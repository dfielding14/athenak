#!/usr/bin/env python3
"""Audit one legacy GOTHAM per-node restart without reading its cell payloads."""

import argparse
import hashlib
import json
from pathlib import Path
import stat
import struct

import numpy as np


_MAX_SPAN_BYTES = 256 * 1024 * 1024
_RANKS_PER_NODE = 8


def _sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _parameters(data):
    values = {}
    block = None
    for raw_line in data.decode().splitlines():
        line = raw_line.split('#', 1)[0].strip()
        if line.startswith('<') and line.endswith('>'):
            block = line[1:-1]
        elif block is not None and '=' in line:
            name, value = line.split('=', 1)
            values[(block, name.strip())] = value.strip()
    return values


def _load_manifest(path):
    data = path.read_bytes()
    marker = b'<par_end>'
    offset = data.index(marker) + len(marker)
    while data[offset:offset + 1] in (b'\r', b'\n'):
        offset += 1
    parameter_dump = data[:offset]

    nmb, root_level = struct.unpack_from('<ii', data, offset)
    offset += struct.calcsize('<2i')
    mesh_size = struct.unpack_from('<9d', data, offset)
    offset += struct.calcsize('<9d')
    mesh_indcs = struct.unpack_from('<19i', data, offset)
    offset += struct.calcsize('<19i')
    mb_indcs = struct.unpack_from('<19i', data, offset)
    offset += struct.calcsize('<19i')
    time, dt = struct.unpack_from('<2d', data, offset)
    offset += struct.calcsize('<2d')
    cycle, nranks, nnodes, version = struct.unpack_from('<4i', data, offset)
    offset += struct.calcsize('<4i')
    if version != 1:
        raise ValueError('unsupported per-node manifest version: {0}'.format(version))

    locations = np.frombuffer(data, dtype='<i4', count=4 * nmb,
                              offset=offset).reshape(nmb, 4).copy()
    offset += 4 * nmb * np.dtype('<i4').itemsize
    costs = np.frombuffer(data, dtype='<f4', count=nmb, offset=offset).copy()
    offset += nmb * np.dtype('<f4').itemsize
    rank_each = np.frombuffer(data, dtype='<i4', count=nmb, offset=offset).copy()
    offset += nmb * np.dtype('<i4').itemsize
    gids_eachrank = np.frombuffer(data, dtype='<i4', count=nranks, offset=offset).copy()
    offset += nranks * np.dtype('<i4').itemsize
    nmb_eachrank = np.frombuffer(data, dtype='<i4', count=nranks, offset=offset).copy()
    offset += nranks * np.dtype('<i4').itemsize
    rank_to_node = np.frombuffer(data, dtype='<i4', count=nranks, offset=offset).copy()
    offset += nranks * np.dtype('<i4').itemsize
    payload_stride = struct.unpack_from('<Q', data, len(data) - 8)[0]

    return {
        'parameter_dump': parameter_dump,
        'parameters': _parameters(parameter_dump),
        'nmb': nmb,
        'root_level': root_level,
        'mesh_size': mesh_size,
        'mesh_indcs': mesh_indcs,
        'mb_indcs': mb_indcs,
        'time': time,
        'dt': dt,
        'cycle': cycle,
        'nranks': nranks,
        'nnodes': nnodes,
        'locations': locations,
        'costs': costs,
        'rank_each': rank_each,
        'gids_eachrank': gids_eachrank,
        'nmb_eachrank': nmb_eachrank,
        'rank_to_node': rank_to_node,
        'payload_stride': payload_stride,
        'state_bytes': len(data) - 8 - offset,
    }


def _particle_header(path):
    with path.open('rb') as stream:
        header = stream.read(struct.calcsize('<d4i'))
    time, cycle, count, nrdata, nidata = struct.unpack('<d4i', header)
    expected_size = struct.calcsize('<d4i') + count * (nrdata * 8 + nidata * 4)
    return {
        'time': time,
        'cycle': cycle,
        'count': count,
        'real_fields': nrdata,
        'integer_fields': nidata,
        'expected_size': expected_size,
        'actual_size': path.stat().st_size,
    }


def _athena_rank_list(nblocks, nranks):
    if nranks > nblocks:
        raise ValueError('rank count exceeds target block count')
    rank_each = np.empty(nblocks, dtype=np.int32)
    total_cost = np.float32(nblocks)
    target_cost = np.float32(total_cost / nranks)
    my_cost = np.float32(0.0)
    rank = nranks - 1
    for gid in range(nblocks - 1, -1, -1):
        my_cost = np.float32(my_cost + 1.0)
        rank_each[gid] = rank
        if my_cost >= target_cost and rank > 0:
            rank -= 1
            total_cost = np.float32(total_cost - my_cost)
            my_cost = np.float32(0.0)
            target_cost = np.float32(total_cost / (rank + 1))
    if rank != 0:
        raise ValueError('AthenaK load-balance emulation left an empty rank')
    return rank_each


def _target_mapping(meta, physical_level):
    target_level = meta['root_level'] + physical_level
    source = meta['locations']
    depth = np.maximum(source[:, 3] - target_level, 0).astype(np.int32)
    mapped = source.copy()
    mapped[:, 0] >>= depth
    mapped[:, 1] >>= depth
    mapped[:, 2] >>= depth
    mapped[:, 3] -= depth

    first = np.ones(meta['nmb'], dtype=bool)
    first[1:] = np.any(mapped[1:] != mapped[:-1], axis=1)
    target_locations = mapped[first]
    target_gid = np.cumsum(first, dtype=np.int32) - 1
    if np.unique(mapped, axis=0).shape[0] != target_locations.shape[0]:
        raise ValueError('mapped target locations are not contiguous in source Z-order')

    coverage = np.bincount(target_gid, weights=np.ldexp(np.ones(meta['nmb']), -3 * depth),
                           minlength=target_locations.shape[0])
    if not np.array_equal(coverage, np.ones_like(coverage)):
        raise ValueError('target coverage contains a hole or overlap')
    fan_in = np.bincount(target_gid, minlength=target_locations.shape[0])
    fan_values, fan_counts = np.unique(fan_in, return_counts=True)
    levels, level_counts = np.unique(target_locations[:, 3] - meta['root_level'],
                                     return_counts=True)
    return {
        'physical_level': physical_level,
        'logical_level': target_level,
        'locations': target_locations,
        'target_gid': target_gid,
        'depth': depth,
        'block_count': int(target_locations.shape[0]),
        'level_histogram': {str(int(level)): int(count)
                            for level, count in zip(levels, level_counts)},
        'fan_in_histogram': {str(int(value)): int(count)
                             for value, count in zip(fan_values, fan_counts)},
        'max_fan_in': int(fan_in.max()),
        'coverage_min': float(coverage.min()),
        'coverage_max': float(coverage.max()),
    }


def _source_chunk_indices(meta):
    prefixes = np.zeros(meta['nranks'], dtype=np.int64)
    node_counts = np.zeros(meta['nnodes'], dtype=np.int64)
    for rank in range(meta['nranks']):
        node = int(meta['rank_to_node'][rank])
        prefixes[rank] = node_counts[node]
        node_counts[node] += int(meta['nmb_eachrank'][rank])
    gids = np.arange(meta['nmb'], dtype=np.int64)
    ranks = meta['rank_each'].astype(np.int64)
    local = gids - meta['gids_eachrank'][ranks]
    return meta['rank_to_node'][ranks].astype(np.int64), prefixes[ranks] + local


def _io_plan(meta, mapping, candidate_ranks, source_nodes, source_chunks):
    rank_each_target = _athena_rank_list(mapping['block_count'], candidate_ranks)
    rank_each_source = rank_each_target[mapping['target_gid']].astype(np.int64)
    key = rank_each_source * meta['nnodes'] + source_nodes
    contiguous = np.ones(meta['nmb'], dtype=bool)
    contiguous[1:] = ((key[1:] == key[:-1]) & (source_chunks[1:] == source_chunks[:-1] + 1))
    starts = np.flatnonzero(~contiguous)
    starts = np.concatenate((np.array([0], dtype=np.int64), starts))
    starts = np.unique(starts)
    lengths = np.diff(np.append(starts, meta['nmb']))
    run_keys = key[starts]
    max_chunks = _MAX_SPAN_BYTES // meta['payload_stride']
    run_spans = (lengths + max_chunks - 1) // max_chunks

    unique_keys, inverse = np.unique(run_keys, return_inverse=True)
    spans_by_rank_shard = np.bincount(inverse, weights=run_spans).astype(np.int64)
    ranks = unique_keys // meta['nnodes']
    source_shards = unique_keys % meta['nnodes']
    requests_per_rank = np.bincount(rank_each_source, minlength=candidate_ranks)

    open_key = (ranks // _RANKS_PER_NODE) * meta['nnodes'] + source_shards
    unique_opens, open_inverse = np.unique(open_key, return_inverse=True)
    collective_reads = np.zeros(unique_opens.size, dtype=np.int64)
    np.maximum.at(collective_reads, open_inverse, spans_by_rank_shard)
    shards_per_rank = np.bincount(ranks, minlength=candidate_ranks)

    target_counts = np.bincount(rank_each_target, minlength=candidate_ranks)
    return {
        'ranks': candidate_ranks,
        'nodes_at_8_ranks_per_node': (candidate_ranks + _RANKS_PER_NODE - 1)
                                     // _RANKS_PER_NODE,
        'target_blocks_per_rank_min': int(target_counts.min()),
        'target_blocks_per_rank_max': int(target_counts.max()),
        'source_blocks_requested_per_rank_min': int(requests_per_rank.min()),
        'source_blocks_requested_per_rank_max': int(requests_per_rank.max()),
        'source_shards_touched_per_rank_min': int(shards_per_rank.min()),
        'source_shards_touched_per_rank_max': int(shards_per_rank.max()),
        'rank_shard_request_groups': int(unique_keys.size),
        'estimated_contiguous_read_spans': int(spans_by_rank_shard.sum()),
        'estimated_collective_read_calls': int(collective_reads.sum()),
        'conversion_node_source_shard_opens': int(unique_opens.size),
        'max_spans_for_one_rank_source_shard': int(spans_by_rank_shard.max()),
        'max_span_bytes': _MAX_SPAN_BYTES,
    }


def _inventory(meta, manifest, inventory_path):
    expected_blocks = np.bincount(meta['rank_to_node'], weights=meta['nmb_eachrank'],
                                  minlength=meta['nnodes']).astype(np.int64)
    lines = ['node\trelative_path\tblocks\texpected_bytes\tactual_bytes\tmtime_ns\t'
             'inode\tdevice\tstatus']
    missing = []
    wrong_size = []
    actual_paths = set()
    total_bytes = 0
    for node in range(meta['nnodes']):
        relative = Path('node_{0:08d}'.format(node)) / manifest.name
        path = manifest.parent / relative
        expected_size = int(expected_blocks[node]) * meta['payload_stride']
        if not path.is_file():
            missing.append(str(relative))
            lines.append('{0}\t{1}\t{2}\t{3}\t-\t-\t-\t-\tmissing'.format(
                node, relative, int(expected_blocks[node]), expected_size))
            continue
        info = path.stat()
        actual_paths.add(path.resolve())
        total_bytes += info.st_size
        status = 'ok' if info.st_size == expected_size else 'wrong_size'
        if status != 'ok':
            wrong_size.append(str(relative))
        lines.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(
            node, relative, int(expected_blocks[node]), expected_size, info.st_size,
            info.st_mtime_ns, info.st_ino, info.st_dev, status))

    globbed = {path.resolve() for path in manifest.parent.glob('node_*/' + manifest.name)}
    extras = sorted(str(path) for path in globbed - actual_paths)
    text = '\n'.join(lines) + '\n'
    inventory_path.write_text(text)
    return {
        'expected_shards': meta['nnodes'],
        'present_expected_shards': meta['nnodes'] - len(missing),
        'missing': missing,
        'wrong_size': wrong_size,
        'extra_matching_payloads': extras,
        'complete': not missing and not wrong_size and not extras,
        'expected_payload_bytes': int(meta['nmb'] * meta['payload_stride']),
        'actual_expected_payload_bytes': total_bytes,
        'inventory_path': str(inventory_path.resolve()),
        'inventory_sha256': hashlib.sha256(text.encode()).hexdigest(),
    }


def audit(args):
    manifest = args.fluid.resolve()
    particles = args.particles.resolve()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    inventory_path = output.with_suffix('.shards.tsv')
    meta = _load_manifest(manifest)
    particle = _particle_header(particles)

    params = meta['parameters']
    expected_stride = (8 + 3) * (meta['mb_indcs'][1] + 2 * meta['mb_indcs'][0]) ** 3 * 8
    checks = {
        'adaptive': params.get(('mesh_refinement', 'refinement')) == 'adaptive',
        'double_precision_gotham_stride': meta['payload_stride'] == expected_stride,
        'inactive_turbulence': (int(params[('turb_driving', 'turb_flag')]) == 1
                                and float(params[('turb_driving', 'dedt')]) == 0.0
                                and float(params[('turb_driving', 'tdriv_start')]) == 100.0),
        'particle_layout': particle['count'] > 0 and particle['real_fields'] == 9
                           and particle['integer_fields'] == 3,
        'particle_size': particle['expected_size'] == particle['actual_size'],
        'fluid_particle_time': particle['time'] == meta['time'],
        'fluid_particle_cycle': particle['cycle'] == meta['cycle'],
        'rank_block_count': int(meta['nmb_eachrank'].sum()) == meta['nmb'],
        'rank_gid_start': int(meta['gids_eachrank'][0]) == 0,
        'rank_gid_contiguous': np.array_equal(meta['gids_eachrank'][1:],
                                              np.cumsum(meta['nmb_eachrank'][:-1])),
        'rank_assignments': np.array_equal(meta['rank_each'],
                                           np.repeat(np.arange(meta['nranks']),
                                                     meta['nmb_eachrank'])),
        'node_ids': (int(meta['rank_to_node'].min()) == 0
                     and int(meta['rank_to_node'].max()) == meta['nnodes'] - 1),
    }

    source_levels, source_level_counts = np.unique(
        meta['locations'][:, 3] - meta['root_level'], return_counts=True)
    source_nodes, source_chunks = _source_chunk_indices(meta)
    targets = []
    for physical_level in args.target_physical_levels:
        mapping = _target_mapping(meta, physical_level)
        mapping['io_plans'] = [
            _io_plan(meta, mapping, ranks, source_nodes, source_chunks)
            for ranks in args.candidate_ranks if ranks <= mapping['block_count']
        ]
        del mapping['locations']
        del mapping['target_gid']
        del mapping['depth']
        targets.append(mapping)

    inventory = _inventory(meta, manifest, inventory_path)
    checks['complete_shards'] = inventory['complete']
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ValueError('audit checks failed: ' + ', '.join(failed))

    manifest_stat = manifest.stat()
    particle_stat = particles.stat()
    report = {
        'label': args.label,
        'audit_scope': 'root manifest, particle sidecar, and payload shard metadata only',
        'payload_content_read': False,
        'fluid': {
            'path': str(manifest),
            'sha256': _sha256(manifest),
            'size': manifest_stat.st_size,
            'mtime_ns': manifest_stat.st_mtime_ns,
            'inode': manifest_stat.st_ino,
            'mode': stat.filemode(manifest_stat.st_mode),
            'time_code': meta['time'],
            'time_gyr': meta['time'] * 0.1,
            'dt_code': meta['dt'],
            'cycle': meta['cycle'],
            'source_blocks': meta['nmb'],
            'source_physical_level_histogram': {
                str(int(level)): int(count)
                for level, count in zip(source_levels, source_level_counts)
            },
            'root_logical_level': meta['root_level'],
            'source_max_logical_level': int(meta['locations'][:, 3].max()),
            'source_max_physical_level': int(meta['locations'][:, 3].max())
                                         - meta['root_level'],
            'original_ranks': meta['nranks'],
            'original_nodes': meta['nnodes'],
            'payload_stride': meta['payload_stride'],
            'payload_bytes': int(meta['nmb'] * meta['payload_stride']),
            'payload_tib': meta['nmb'] * meta['payload_stride'] / 2 ** 40,
        },
        'particles': {
            'path': str(particles),
            'sha256': _sha256(particles),
            'size': particle_stat.st_size,
            'mtime_ns': particle_stat.st_mtime_ns,
            **particle,
        },
        'checks': checks,
        'shards': inventory,
        'targets': targets,
    }
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--label', required=True)
    parser.add_argument('--fluid', required=True, type=Path)
    parser.add_argument('--particles', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--target-physical-levels', type=int, nargs='+', default=[6, 5])
    parser.add_argument('--candidate-ranks', type=int, nargs='+',
                        default=[1024, 2048, 4096, 8192])
    args = parser.parse_args()
    report = audit(args)
    print(json.dumps({
        'label': report['label'],
        'source_blocks': report['fluid']['source_blocks'],
        'time_gyr': report['fluid']['time_gyr'],
        'cycle': report['fluid']['cycle'],
        'payload_tib': report['fluid']['payload_tib'],
        'shards_complete': report['shards']['complete'],
        'targets': [{'physical_level': target['physical_level'],
                     'block_count': target['block_count'],
                     'max_fan_in': target['max_fan_in']}
                    for target in report['targets']],
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
