#!/usr/bin/env python3
"""Apply the Phase-4 acceptance gates to one real GOTHAM conversion."""

import argparse
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess

import numpy as np

from audit_gotham_restart_index import _load_manifest


_FIELDS = ('density', 'momentum_x1', 'momentum_x2', 'momentum_x3',
           'total_energy', 'metallicity', 'dust_small', 'dust_large')


def _sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _key_values(path):
    return dict(line.split('=', 1) for line in path.read_text().splitlines()
                if '=' in line)


def _particles(path):
    data = path.read_bytes()
    time, cycle, count, nrdata, nidata = struct.unpack_from('<d4i', data, 0)
    offset = struct.calcsize('<d4i')
    real = np.frombuffer(data, dtype='<f8', count=count * nrdata,
                         offset=offset).reshape(count, nrdata).copy()
    offset += count * nrdata * 8
    integer = np.frombuffer(data, dtype='<i4', count=count * nidata,
                            offset=offset).reshape(count, nidata).copy()
    if offset + count * nidata * 4 != len(data):
        raise ValueError('particle restart size does not match its header')
    order = np.argsort(integer[:, 1], kind='stable')
    return time, cycle, real[order], integer[order]


def _compare_particles(source, target):
    source_time, source_cycle, source_real, source_int = _particles(source)
    target_time, target_cycle, target_real, target_int = _particles(target)
    if source_time != target_time or source_cycle != target_cycle:
        raise ValueError('particle restart time/cycle changed')
    if not np.array_equal(source_real, target_real):
        raise ValueError('particle real fields changed')
    # Ownership GID is recomputed for the target tree.  Tag, SN counter, and the
    # remaining integer state must be bitwise unchanged.
    if not np.array_equal(source_int[:, 1:], target_int[:, 1:]):
        raise ValueError('particle integer fields other than ownership GID changed')
    return source_real.shape[0]


def _shard_sizes(meta, manifest):
    blocks = np.bincount(meta['rank_to_node'], weights=meta['nmb_eachrank'],
                         minlength=meta['nnodes']).astype(np.int64)
    missing = []
    wrong_size = []
    for node, count in enumerate(blocks):
        path = manifest.parent / ('node_{0:08d}'.format(node)) / manifest.name
        if not path.is_file():
            missing.append(str(path))
        elif path.stat().st_size != int(count) * meta['payload_stride']:
            wrong_size.append(str(path))
    if missing or wrong_size:
        raise ValueError('restart payload shards are missing or incorrectly sized')
    return {'shards': int(meta['nnodes']),
            'payload_bytes': int(blocks.sum()) * meta['payload_stride']}


def _tree_checks(meta, expected_level, expected_blocks):
    locations = meta['locations']
    max_level = meta['root_level'] + expected_level
    if meta['nmb'] != expected_blocks or int(locations[:, 3].max()) != max_level:
        raise ValueError('converted tree block count or maximum level is incorrect')

    root_shape = np.array((meta['mesh_indcs'][3] // meta['mb_indcs'][3],
                           meta['mesh_indcs'][2] // meta['mb_indcs'][2],
                           meta['mesh_indcs'][1] // meta['mb_indcs'][1]),
                          dtype=np.int64)
    scale = 1 << expected_level
    levels = np.full(tuple(root_shape * scale), -1, dtype=np.int8)
    for lx1, lx2, lx3, level in locations:
        width = 1 << (max_level - int(level))
        x0, y0, z0 = int(lx1) * width, int(lx2) * width, int(lx3) * width
        target = levels[z0:z0 + width, y0:y0 + width, x0:x0 + width]
        if np.any(target != -1):
            raise ValueError('converted tree contains overlapping leaves')
        target[...] = int(level)
    if np.any(levels == -1):
        raise ValueError('converted tree contains uncovered volume')

    # Check face, edge, and corner neighbors on the finest target lattice.
    directions = ((1, 0, 0), (0, 1, 0), (0, 0, 1),
                  (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1),
                  (0, 1, 1), (0, 1, -1),
                  (1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1))
    max_jump = 0
    for dz, dy, dx in directions:
        first_slices = tuple(slice(0, -delta) if delta > 0 else
                             slice(-delta, None) if delta < 0 else slice(None)
                             for delta in (dz, dy, dx))
        second_slices = tuple(slice(delta, None) if delta > 0 else
                              slice(0, delta) if delta < 0 else slice(None)
                              for delta in (dz, dy, dx))
        first = levels[first_slices]
        second = levels[second_slices]
        max_jump = max(max_jump, int(np.max(np.abs(first - second))))
    if max_jump > 1:
        raise ValueError('converted tree violates 2:1 balance')
    return {'blocks': int(meta['nmb']), 'maximum_physical_level': expected_level,
            'coverage': 'exact', 'maximum_neighbor_level_jump': max_jump}


def _conversion_metrics(metadata):
    tolerance = Decimal(metadata['conservation_tolerance'])
    residuals = {}
    for field in _FIELDS:
        source = Decimal(metadata['source_integral.' + field])
        target = Decimal(metadata['target_integral.' + field])
        scale = max(Decimal(1), Decimal(metadata['source_abs_integral.' + field]))
        residuals[field] = abs(target - source) / scale
    maximum = max(residuals.values())
    if maximum > tolerance or Decimal(metadata['max_scaled_integral_error']) > tolerance:
        raise ValueError('converted conserved integrals exceed the declared tolerance')
    for name in ('physical_source_invocations', 'particle_updates',
                 'primitive_recovery_floors', 'invalid_target_cells'):
        if int(metadata[name]) != 0:
            raise ValueError(name + ' is nonzero')
    return {
        'tolerance': str(tolerance),
        'maximum_scaled_integral_error': str(maximum),
        'scaled_integral_errors': {key: str(value)
                                   for key, value in residuals.items()},
        'resolved_kinetic_change': metadata['resolved_kinetic_change'],
        'internal_energy_change': metadata['internal_energy_change'],
        'energy_partition_closure': metadata['energy_partition_closure'],
        'positive_finite_target': True,
        'physical_source_invocations': 0,
        'particle_updates': 0,
    }


def _one_file(directory, pattern):
    paths = list(directory.glob(pattern))
    if len(paths) != 1:
        raise ValueError('expected exactly one file matching ' + pattern)
    return paths[0]


def _write_atomic(path, text):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(text)
    os.replace(temporary, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--source-fluid', required=True, type=Path)
    parser.add_argument('--source-particles', required=True, type=Path)
    parser.add_argument('--source-fluid-sha256', required=True)
    parser.add_argument('--source-particle-sha256', required=True)
    parser.add_argument('--executable', required=True, type=Path)
    parser.add_argument('--executable-sha256', required=True)
    parser.add_argument('--source-worktree', required=True, type=Path)
    parser.add_argument('--source-commit', required=True)
    parser.add_argument('--target-physical-level', required=True, type=int)
    parser.add_argument('--target-blocks', required=True, type=int)
    parser.add_argument('--roundtrip', type=Path)
    args = parser.parse_args()

    directory = args.directory.resolve()
    manifest_path = directory / 'conversion_manifest.txt'
    if not (directory / 'CONVERSION_COMPLETE').is_file():
        raise SystemExit('conversion has not been hashed and finalized')
    metadata = _key_values(manifest_path)

    source_fluid = args.source_fluid.resolve()
    source_particles = args.source_particles.resolve()
    executable = args.executable.resolve()
    if _sha256(source_fluid) != args.source_fluid_sha256:
        raise SystemExit('source fluid manifest hash changed')
    if _sha256(source_particles) != args.source_particle_sha256:
        raise SystemExit('source particle hash changed')
    if _sha256(executable) != args.executable_sha256:
        raise SystemExit('conversion executable hash changed')
    head = subprocess.check_output(
        ['git', '-C', args.source_worktree, 'rev-parse', 'HEAD'], text=True).strip()
    status = subprocess.check_output(
        ['git', '-C', args.source_worktree, 'status', '--porcelain'], text=True)
    source_commit_exists = subprocess.run(
        ['git', '-C', args.source_worktree, 'cat-file', '-e',
         args.source_commit + '^{commit}']).returncode == 0
    if not source_commit_exists or status:
        raise SystemExit('conversion source commit is unavailable or worktree is dirty')

    if Path(metadata['source_fluid']) != source_fluid:
        raise SystemExit('conversion manifest names a different source fluid restart')
    if Path(metadata['source_particles']) != source_particles:
        raise SystemExit('conversion manifest names a different particle sidecar')
    if int(metadata['target_max_physical_level']) != args.target_physical_level:
        raise SystemExit('conversion manifest target level is incorrect')
    if int(metadata['target_blocks']) != args.target_blocks:
        raise SystemExit('conversion manifest target block count is incorrect')

    fluid = Path(metadata['converted_fluid'])
    particles = Path(metadata['converted_particles'])
    source_meta = _load_manifest(source_fluid)
    target_meta = _load_manifest(fluid)
    report = {
        'status': 'conversion_validated',
        'source': {
            'fluid': str(source_fluid),
            'fluid_sha256': args.source_fluid_sha256,
            'particles': str(source_particles),
            'particle_sha256': args.source_particle_sha256,
            'time': source_meta['time'],
            'cycle': source_meta['cycle'],
        },
        'implementation': {
            'worktree': str(args.source_worktree.resolve()),
            'commit': args.source_commit,
            'validator_commit': head,
            'diff': 'clean',
            'executable': str(executable),
            'executable_sha256': args.executable_sha256,
        },
        'target': _tree_checks(target_meta, args.target_physical_level,
                               args.target_blocks),
        'shards': _shard_sizes(target_meta, fluid),
        'conversion': _conversion_metrics(metadata),
        'particles': {
            'count': _compare_particles(source_particles, particles),
            'fields_unchanged_except_ownership_gid': True,
            'sn_state_unchanged': True,
        },
        'time_cycle_dt_preserved': (
            (source_meta['time'], source_meta['cycle'], source_meta['dt']) ==
            (target_meta['time'], target_meta['cycle'], target_meta['dt'])
            and float(metadata['time']) == source_meta['time']
            and int(metadata['cycle']) == source_meta['cycle']
            and float(metadata['dt']) == source_meta['dt']),
        'payload_content_read_by_validator': False,
    }
    if not report['time_cycle_dt_preserved']:
        raise SystemExit('converted time, cycle, or timestep changed')

    if args.roundtrip is not None:
        roundtrip_dir = args.roundtrip.resolve()
        roundtrip_fluid = _one_file(roundtrip_dir / 'rst', '*.rst')
        roundtrip_particles = _one_file(roundtrip_dir / 'rst_prtcl', '*.rst_prtcl')
        roundtrip_meta = _load_manifest(roundtrip_fluid)
        if not np.array_equal(target_meta['locations'], roundtrip_meta['locations']):
            raise SystemExit('ordinary round trip changed the target topology')
        if (target_meta['time'], target_meta['cycle']) != (
                roundtrip_meta['time'], roundtrip_meta['cycle']):
            raise SystemExit('ordinary round trip changed time or cycle')
        if not np.isfinite(roundtrip_meta['dt']) or roundtrip_meta['dt'] <= 0.0:
            raise SystemExit('ordinary restart produced an invalid timestep')
        _compare_particles(particles, roundtrip_particles)
        report['roundtrip'] = {
            'directory': str(roundtrip_dir),
            'stored_ranks': target_meta['nranks'],
            'current_ranks': roundtrip_meta['nranks'],
            'topology_identical': True,
            'time_cycle_identical': True,
            'timestep_policy': 'recomputed_by_ordinary_restart_initialization',
            'input_dt': target_meta['dt'],
            'output_dt': roundtrip_meta['dt'],
            'output_to_input_dt_ratio': roundtrip_meta['dt'] / target_meta['dt'],
            'particles_identical_except_ownership_gid': True,
            **_shard_sizes(roundtrip_meta, roundtrip_fluid),
        }
        if target_meta['nranks'] == roundtrip_meta['nranks']:
            raise SystemExit('ordinary round trip did not change MPI decomposition')
        report['status'] = 'accepted'

    report_path = directory / 'phase4_validation.json'
    _write_atomic(report_path, json.dumps(report, indent=2, sort_keys=True) + '\n')
    lines = [line for line in manifest_path.read_text().splitlines()
             if not line.startswith('phase4.')]
    lines.extend((
        'phase4.source_fluid_sha256=' + args.source_fluid_sha256,
        'phase4.source_particles_sha256=' + args.source_particle_sha256,
        'phase4.executable_sha256=' + args.executable_sha256,
        'phase4.source_commit=' + args.source_commit,
        'phase4.source_diff=clean',
        'phase4.target_physical_level=' + str(args.target_physical_level),
        'phase4.target_blocks=' + str(args.target_blocks),
        'phase4.validation=' + report['status'],
        'phase4.validation_report=' + str(report_path),
    ))
    _write_atomic(manifest_path, '\n'.join(lines) + '\n')
    marker = 'PHASE4_ACCEPTED' if args.roundtrip is not None \
        else 'PHASE4_CONVERSION_VALIDATED'
    _write_atomic(directory / marker, report['status'] + '\n')


if __name__ == '__main__':
    main()
