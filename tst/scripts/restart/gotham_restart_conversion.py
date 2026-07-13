import json
import logging
from pathlib import Path
import shutil
import struct
import subprocess

import numpy as np
import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_REPO_DIR = Path(__file__).resolve().parents[3]
_EXE = _REPO_DIR / 'tst' / 'build' / 'src' / 'athena'
_SOURCE_INPUT = _REPO_DIR / 'inputs' / 'tests' / 'gotham_restart_conversion_source.athinput'
_CONVERSION_INPUT = _REPO_DIR / 'inputs' / 'tests' / 'gotham_restart_conversion.athinput'
_WORK = _REPO_DIR / 'tst' / 'build' / 'phase3_gotham_conversion'
_RTOL = 5.0e-13
_ATOL = 5.0e-14


class Restart:
    """Read the exact double-precision, per-node layout used by this GOTHAM test."""

    def __init__(self, manifest):
        self.manifest = Path(manifest)
        data = self.manifest.read_bytes()
        marker = b'<par_end>'
        offset = data.index(marker) + len(marker)
        while data[offset:offset + 1] in (b'\r', b'\n'):
            offset += 1
        self.parameter_dump = data[:offset]
        self.binary_offset = offset

        self.nmb, self.root_level = struct.unpack_from('<ii', data, offset)
        offset += 2 * struct.calcsize('<i')
        self.mesh_size = struct.unpack_from('<9d', data, offset)
        offset += struct.calcsize('<9d')
        self.mesh_indcs = struct.unpack_from('<19i', data, offset)
        offset += struct.calcsize('<19i')
        self.mb_indcs = struct.unpack_from('<19i', data, offset)
        offset += struct.calcsize('<19i')
        self.time, self.dt = struct.unpack_from('<2d', data, offset)
        offset += struct.calcsize('<2d')
        self.cycle = struct.unpack_from('<i', data, offset)[0]
        offset += struct.calcsize('<i')
        self.head_before_nranks = data[self.binary_offset:offset]
        self.nranks, self.nnodes, version = struct.unpack_from('<3i', data, offset)
        offset += struct.calcsize('<3i')
        assert version == 1

        self.locations = np.frombuffer(data, dtype='<i4', count=4 * self.nmb,
                                       offset=offset).reshape(self.nmb, 4).copy()
        offset += 4 * self.nmb * np.dtype('<i4').itemsize
        self.costs = np.frombuffer(data, dtype='<f4', count=self.nmb,
                                   offset=offset).copy()
        offset += self.nmb * np.dtype('<f4').itemsize
        self.rank_each = np.frombuffer(data, dtype='<i4', count=self.nmb,
                                       offset=offset).copy()
        offset += self.nmb * np.dtype('<i4').itemsize
        self.gids_eachrank = np.frombuffer(data, dtype='<i4', count=self.nranks,
                                           offset=offset).copy()
        offset += self.nranks * np.dtype('<i4').itemsize
        self.nmb_eachrank = np.frombuffer(data, dtype='<i4', count=self.nranks,
                                          offset=offset).copy()
        offset += self.nranks * np.dtype('<i4').itemsize
        self.rank_to_node = np.frombuffer(data, dtype='<i4', count=self.nranks,
                                          offset=offset).copy()
        offset += self.nranks * np.dtype('<i4').itemsize
        self.state_tail = data[offset:-8]
        self.stride = struct.unpack_from('<Q', data, len(data) - 8)[0]

        self.ng = self.mb_indcs[0]
        self.nx = tuple(self.mb_indcs[1:4])
        assert self.nx[0] == self.nx[1] == self.nx[2]
        self.ncell = self.nx[0] + 2 * self.ng
        self.hydro_values = 8 * self.ncell ** 3
        self.force_values = 3 * self.ncell ** 3
        assert self.stride == (self.hydro_values + self.force_values) * 8

    def _chunk_index(self, gid):
        rank = int(self.rank_each[gid])
        node = int(self.rank_to_node[rank])
        prefix = sum(int(self.nmb_eachrank[r]) for r in range(rank)
                     if int(self.rank_to_node[r]) == node)
        return node, prefix + gid - int(self.gids_eachrank[rank])

    def active_hydro(self):
        active = np.empty((self.nmb, 8, self.nx[2], self.nx[1], self.nx[0]))
        shards = {}
        for node in range(self.nnodes):
            path = self.manifest.parent / ('node_{0:08d}'.format(node)) / self.manifest.name
            shards[node] = np.memmap(path, mode='r', dtype='u1')
        for gid in range(self.nmb):
            node, chunk = self._chunk_index(gid)
            block = np.frombuffer(shards[node], dtype='<f8', count=self.hydro_values,
                                  offset=chunk * self.stride)
            block = block.reshape(8, self.ncell, self.ncell, self.ncell)
            s = slice(self.ng, self.ng + self.nx[0])
            active[gid] = block[:, s, s, s]
        return active

    def force_nonzero(self):
        nonzero = 0
        shards = {}
        for node in range(self.nnodes):
            path = self.manifest.parent / ('node_{0:08d}'.format(node)) / self.manifest.name
            shards[node] = np.memmap(path, mode='r', dtype='u1')
        for gid in range(self.nmb):
            node, chunk = self._chunk_index(gid)
            force = np.frombuffer(shards[node], dtype='<f8', count=self.force_values,
                                  offset=chunk * self.stride + self.hydro_values * 8)
            nonzero += int(np.count_nonzero(force))
        return nonzero


def _run(command):
    logger.debug('Executing: %s', ' '.join(str(value) for value in command))
    result = subprocess.run([str(value) for value in command], cwd=str(_REPO_DIR),
                            text=True, capture_output=True)
    if result.returncode != 0:
        logger.error('Athena failed:\n%s%s', result.stdout, result.stderr)
        result.check_returncode()


def _balanced_counts(total, count):
    return np.array([total // count + (rank < total % count) for rank in range(count)],
                    dtype='<i4')


def _repack_per_node(source, destination, nranks, nnodes):
    """Re-express a seed with a different stored rank/node decomposition."""
    source = Restart(source)
    destination = Path(destination)
    destination.parent.mkdir(parents=True)
    counts = _balanced_counts(source.nmb, nranks)
    gids = np.zeros(nranks, dtype='<i4')
    gids[1:] = np.cumsum(counts[:-1])
    rank_each = np.repeat(np.arange(nranks, dtype='<i4'), counts)
    rank_to_node = np.array([rank * nnodes // nranks for rank in range(nranks)],
                            dtype='<i4')

    manifest = bytearray(source.parameter_dump)
    manifest.extend(source.head_before_nranks)
    manifest.extend(struct.pack('<3i', nranks, nnodes, 1))
    manifest.extend(source.locations.astype('<i4').tobytes())
    manifest.extend(source.costs.astype('<f4').tobytes())
    manifest.extend(rank_each.tobytes())
    manifest.extend(gids.tobytes())
    manifest.extend(counts.tobytes())
    manifest.extend(rank_to_node.tobytes())
    manifest.extend(source.state_tail)
    manifest.extend(struct.pack('<Q', source.stride))
    destination.write_bytes(manifest)

    source_shard = source.manifest.parent / 'node_00000000' / source.manifest.name
    payload = np.memmap(source_shard, mode='r', dtype='u1')
    for node in range(nnodes):
        node_dir = destination.parent / ('node_{0:08d}'.format(node))
        node_dir.mkdir()
        with (node_dir / destination.name).open('wb') as stream:
            for rank in range(nranks):
                if int(rank_to_node[rank]) != node:
                    continue
                begin = int(gids[rank]) * source.stride
                end = begin + int(counts[rank]) * source.stride
                stream.write(payload[begin:end])


def _move_particles_to_boundaries(path):
    data = bytearray(Path(path).read_bytes())
    time, cycle, count, nrdata, nidata = struct.unpack_from('<d4i', data, 0)
    assert count > 16 and nrdata == 9 and nidata == 3
    rdata = np.frombuffer(data, dtype='<f8', count=count * nrdata,
                          offset=struct.calcsize('<d4i')).reshape(count, nrdata)
    eps = 1.0e-9
    points = [
        (-eps, -eps, -eps), (eps, eps, eps),
        (-eps, eps, -eps), (eps, -eps, eps),
        (50.0 - eps, eps, eps), (50.0 + eps, eps, eps),
        (-50.0 - eps, -eps, eps), (-50.0 + eps, -eps, eps),
    ]
    for particle, (x1, x2, x3) in enumerate(points):
        rdata[particle, 0] = x1
        rdata[particle, 2] = x2
        rdata[particle, 4] = x3
    Path(path).write_bytes(data)
    return {'time': time, 'cycle': cycle, 'count': count,
            'boundary_particles': len(points)}


def _particles(path):
    data = Path(path).read_bytes()
    time, cycle, count, nrdata, nidata = struct.unpack_from('<d4i', data, 0)
    offset = struct.calcsize('<d4i')
    real = np.frombuffer(data, dtype='<f8', count=count * nrdata,
                         offset=offset).reshape(count, nrdata).copy()
    offset += count * nrdata * 8
    integer = np.frombuffer(data, dtype='<i4', count=count * nidata,
                            offset=offset).reshape(count, nidata).copy()
    order = np.argsort(integer[:, 1], kind='stable')
    return time, cycle, real[order], integer[order]


def _restrict_once(values):
    # This is the eight-term cell-centered expression in MeshRefinement::RestrictCC.
    return 0.125 * (
        values[:, 0::2, 0::2, 0::2] + values[:, 0::2, 0::2, 1::2]
        + values[:, 0::2, 1::2, 0::2] + values[:, 0::2, 1::2, 1::2]
        + values[:, 1::2, 0::2, 0::2] + values[:, 1::2, 0::2, 1::2]
        + values[:, 1::2, 1::2, 0::2] + values[:, 1::2, 1::2, 1::2])


def _native_oracle(source, target, target_level, source_values):
    expected_locations = {}
    target_logical_level = source.root_level + target_level
    for location in source.locations:
        depth = max(0, int(location[3]) - target_logical_level)
        mapped = (int(location[0]) >> depth, int(location[1]) >> depth,
                  int(location[2]) >> depth, int(location[3]) - depth)
        expected_locations[mapped] = None
    target_locations = {tuple(int(value) for value in location): gid
                        for gid, location in enumerate(target.locations)}
    assert set(expected_locations) == set(target_locations)

    oracle = np.empty((target.nmb, 8, target.nx[2], target.nx[1], target.nx[0]))
    coverage = np.zeros((target.nmb, target.nx[2], target.nx[1], target.nx[0]),
                        dtype=np.uint8)
    for source_gid, location in enumerate(source.locations):
        depth = max(0, int(location[3]) - target_logical_level)
        mapped = (int(location[0]) >> depth, int(location[1]) >> depth,
                  int(location[2]) >> depth, int(location[3]) - depth)
        target_gid = target_locations[mapped]
        restricted = source_values[source_gid]
        for _ in range(depth):
            restricted = _restrict_once(restricted)
        subnx = source.nx[0] >> depth
        x0 = (int(location[0]) - (mapped[0] << depth)) * subnx
        y0 = (int(location[1]) - (mapped[1] << depth)) * subnx
        z0 = (int(location[2]) - (mapped[2] << depth)) * subnx
        oracle[target_gid, :, z0:z0 + subnx, y0:y0 + subnx, x0:x0 + subnx] = restricted
        coverage[target_gid, z0:z0 + subnx, y0:y0 + subnx, x0:x0 + subnx] += 1
    assert np.all(coverage == 1)
    return oracle


def _integrals(restart, values):
    totals = np.zeros(8, dtype=np.longdouble)
    absolute_inventory = np.zeros(8, dtype=np.longdouble)
    kinetic = np.longdouble(0.0)
    internal = np.longdouble(0.0)
    root_volume = np.longdouble(restart.mesh_size[6] * restart.mesh_size[7]
                                * restart.mesh_size[8])
    for gid, location in enumerate(restart.locations):
        depth = int(location[3]) - restart.root_level
        volume = root_volume / np.longdouble(8 ** depth)
        block = values[gid].astype(np.longdouble)
        totals += block.sum(axis=(1, 2, 3)) * volume
        absolute_inventory += np.abs(block).sum(axis=(1, 2, 3)) * volume
        block_kinetic = 0.5 * (block[1] ** 2 + block[2] ** 2 + block[3] ** 2) / block[0]
        kinetic += block_kinetic.sum() * volume
        internal += (block[4] - block_kinetic).sum() * volume
    return totals, absolute_inventory, kinetic, internal


def _compare_particles(source, converted):
    source_time, source_cycle, source_real, source_int = _particles(source)
    target_time, target_cycle, target_real, target_int = _particles(converted)
    assert source_time == target_time and source_cycle == target_cycle
    np.testing.assert_array_equal(source_real, target_real)
    np.testing.assert_array_equal(source_int[:, 1:], target_int[:, 1:])
    return source_real.shape[0]


def _validate_case(source_manifest, source_particles, target_dir, basename, target_level):
    target_manifest = target_dir / 'rst' / (basename + '.00000.rst')
    target_particles = target_dir / 'rst_prtcl' / (basename + '.00000.rst_prtcl')
    source = Restart(source_manifest)
    target = Restart(target_manifest)
    source_values = source.active_hydro()
    target_values = target.active_hydro()
    oracle = _native_oracle(source, target, target_level, source_values)

    difference = np.abs(target_values - oracle)
    scale = np.maximum(1.0, np.abs(oracle))
    max_scaled = float(np.max(difference / scale))
    assert np.allclose(target_values, oracle, rtol=_RTOL, atol=_ATOL)

    source_totals, source_inventory, source_kinetic, source_internal = \
        _integrals(source, source_values)
    target_totals, _, target_kinetic, target_internal = _integrals(target, target_values)
    # Net momentum nearly cancels in this symmetric state, so normalize conservation
    # residuals by each field's absolute volume inventory rather than its signed total.
    integral_scale = np.maximum(np.longdouble(1.0), source_inventory)
    integral_error = np.max(np.abs(target_totals - source_totals) / integral_scale)
    assert integral_error <= np.longdouble(_RTOL)
    assert target.force_nonzero() == 0
    assert source.time == target.time and source.cycle == target.cycle
    assert source.dt == target.dt
    particle_count = _compare_particles(source_particles, target_particles)

    return {
        'source_blocks': source.nmb,
        'target_blocks': target.nmb,
        'target_physical_level': target_level,
        'max_cell_scaled_error': max_scaled,
        'max_integral_l1_scaled_error': float(integral_error),
        'resolved_kinetic_change': float(target_kinetic - source_kinetic),
        'internal_energy_change': float(target_internal - source_internal),
        'energy_partition_closure': float((target_kinetic - source_kinetic)
                                          + (target_internal - source_internal)),
        'particles': particle_count,
    }


def _generate_source(name, uniform):
    output = _WORK / name
    args = [_EXE, '-i', _SOURCE_INPUT, '-d', output, 'job/basename=' + name]
    if uniform:
        for axis in ('x1', 'x2', 'x3'):
            args.extend(['refinement1/{0}min=-200.0'.format(axis),
                         'refinement1/{0}max=200.0'.format(axis)])
    _run(args)
    manifest = output / 'rst' / (name + '.00000.rst')
    particles = output / 'rst_prtcl' / (name + '.00000.rst_prtcl')
    particle_info = _move_particles_to_boundaries(particles)

    repacked = _WORK / (name + '_source4') / 'rst' / manifest.name
    _repack_per_node(manifest, repacked, nranks=4, nnodes=2)
    return repacked, particles, particle_info


def _convert(source, particles, case, target_level):
    output = _WORK / case
    basename = case + '_converted'
    _run([_EXE, '-r', source, '-i', _CONVERSION_INPUT, '-d', output,
          'job/basename=' + basename,
          'restart_conversion/target_physical_level=' + str(target_level),
          'restart_conversion/particle_restart=' + str(particles)])
    assert (output / 'CONVERSION_WRITTEN').is_file()
    return output, basename


def _roundtrip_different_stored_rank_count(target_dir, basename):
    manifest = target_dir / 'rst' / (basename + '.00000.rst')
    repacked = _WORK / 'roundtrip_source3' / 'rst' / manifest.name
    _repack_per_node(manifest, repacked, nranks=3, nnodes=2)
    output = _WORK / 'roundtrip_current1'
    _run([_EXE, '-r', repacked, '-d', output, 'time/nlim=0', 'time/tlim=0.0'])

    source = Restart(manifest)
    roundtrip_manifests = sorted((output / 'rst').glob('*.rst'))
    assert len(roundtrip_manifests) == 1
    target = Restart(roundtrip_manifests[0])
    assert set(map(tuple, source.locations)) == set(map(tuple, target.locations))
    np.testing.assert_array_equal(source.active_hydro(), target.active_hydro())
    source_particle = target_dir / 'rst_prtcl' / (basename + '.00000.rst_prtcl')
    roundtrip_particles = sorted((output / 'rst_prtcl').glob('*.rst_prtcl'))
    assert len(roundtrip_particles) == 1
    _compare_particles(source_particle, roundtrip_particles[0])
    return {'stored_source_ranks': 3, 'stored_source_nodes': 2, 'current_ranks': 1}


def run(**kwargs):
    logger.debug('Running test %s', __name__)
    if _WORK.is_dir():
        shutil.rmtree(_WORK)
    _WORK.mkdir(parents=True)

    mixed, mixed_particles, particle_info = _generate_source('phase3_mixed', uniform=False)
    uniform, uniform_particles, _ = _generate_source('phase3_uniform', uniform=True)
    results = []
    roundtrip_target = None
    for target_level in (2, 1, 0):
        case = 'mixed_l{0}'.format(target_level)
        output, basename = _convert(mixed, mixed_particles, case, target_level)
        results.append(_validate_case(mixed, mixed_particles, output, basename, target_level))
        if target_level == 1:
            roundtrip_target = (output, basename)

    output, basename = _convert(uniform, uniform_particles, 'uniform_l0', 0)
    results.append(_validate_case(uniform, uniform_particles, output, basename, 0))
    roundtrip = _roundtrip_different_stored_rank_count(*roundtrip_target)
    report = {
        'cell_tolerance': {'rtol': _RTOL, 'atol': _ATOL},
        'source_stored_decomposition': {'ranks': 4, 'nodes': 2},
        'current_conversion_ranks': 1,
        'particle_boundary_cases': particle_info['boundary_particles'],
        'cases': results,
        'roundtrip': roundtrip,
    }
    (_WORK / 'synthetic_results.json').write_text(
        json.dumps(report, indent=2, sort_keys=True) + '\n')
    logger.info('GOTHAM restart-conversion Phase 3 results:\n%s',
                json.dumps(report, indent=2, sort_keys=True))


def analyze():
    return True
