"""WO2 scratch restart comparison for the pinned double-precision Linux ABI.

By default only nine unused mesh_indcs coarse-index integers are normalized. AthenaK
writes RegionIndcs as raw bytes, and these fields are not initialized for the
root mesh. The MeshBlock RegionIndcs, time, diagnostics, and every state byte
remain intact. An explicit opt-in additionally normalizes dormant startup forcing
RNG members after validating their serialized state. Original files are never changed.
"""
import hashlib
from pathlib import Path
import re
import struct


def _dormant_forcing_ranges(raw, header, nmb):
    """Recognize the pinned native v3 forcing record; retain RNG padding and seed."""
    blocks = {}
    block = None
    for line in raw[:header].decode().splitlines():
        line = line.split('#', 1)[0]
        match = re.match(r'\s*<([^>]+)>', line)
        if match:
            block = match[1]
            blocks[block] = {}
        elif block and '=' in line:
            key, value = line.split('=', 1)
            blocks[block][key.strip()] = value.strip()
    info = {'enabled': True, 'applied': False}
    if 'turb_driving' not in blocks:
        return [], dict(info, reason='no forcing record')
    time, _, cycle = struct.unpack_from('=ddi', raw, header + 232)
    if time != 0.0 or cycle != 0:
        return [], dict(info, reason='not initial time and cycle', time=time, cycle=cycle)
    assert 'z4c' not in blocks, 'unsupported metadata before the forcing record'
    start = header + 252 + nmb*(16 + 4)
    integers = struct.unpack_from('=24i', raw, start)
    assert integers[0] == 3 and integers[1] > 0, 'unexpected forcing metadata signature'
    config = blocks['turb_driving']
    keys = ('nlow', 'nhigh', 'driving_type', 'min_kx', 'max_kx', 'min_ky', 'max_ky',
            'min_kz', 'max_kz')
    for index, key in enumerate(keys, 3):
        assert integers[index] == int(config[key]), ('forcing signature', key)
    assert integers[12] == int('npeak' in config), 'forcing npeak signature'
    for index, key in enumerate(('turb_flag', 'tile_nx', 'tile_ny', 'tile_nz'), 13):
        assert integers[index] == int(config[key]), ('forcing signature', key)
    enums = {'normalization': ('edot', 'accel_rms'),
             'localization': ('none', 'include', 'exclude'),
             'spectrum': ('parabolic', 'power_law'),
             'projection_policy': ('solenoidal_compressive',
                                   'mks24_random_unprojected',
                                   'mks24_alfvenic_perpendicular')}
    for index, (key, choices) in enumerate(enums.items(), 17):
        assert integers[index] == choices.index(config[key]), ('forcing signature', key)
    for index, key in enumerate(('physical_k_shell', 'isotropic_power_spectrum',
                                 'record_injected_work'), 21):
        value = config[key].lower()
        expected = int(value in ('1', 'true'))
        assert value in ('0', '1', 'false', 'true') and integers[index] == expected
    rng = start + 24*4 + 19*8
    idum, = struct.unpack_from('=q', raw, rng)
    iset, = struct.unpack_from('=i', raw, rng + 280)
    info.update(metadata_version=integers[0], mode_count=integers[1],
                n_updates=integers[2], time=time, cycle=cycle,
                rng_offset=rng, rng_bytes=296, idum=idum, iset=iset,
                retained_padding_range=[rng + 284, rng + 288])
    if integers[2] != 0 or idum != -1 or iset != 0:
        return [], dict(info, reason='forcing RNG is not the verified dormant seed state')
    ranges = [[rng + 8, rng + 280], [rng + 288, rng + 296]]
    info.update(applied=True, ranges=ranges,
                members=['idum2', 'iy', 'iv[0:32]', 'gset'],
                reason='negative seed initializes shuffle state and zero cache flag initializes gset before use',
                source='src/srcterms/turb_driver.cpp:507; src/utils/random.hpp:55,152; src/outputs/restart.cpp:329')
    return ranges, info


def normalized_restart(path, *, normalize_dormant_startup_forcing=False):
    path = Path(path)
    raw = path.read_bytes()
    marker = b'<par_end>\n'
    header = raw.index(marker) + len(marker)
    # Header: nmb_total, root_level, RegionSize(9 doubles), RegionIndcs(19 ints).
    nmb, root_level = struct.unpack_from('=ii', raw, header)
    assert nmb > 0 and root_level >= 0, 'unexpected restart ABI'
    indcs = header + 2*4 + 9*8
    ng, nx1, nx2, nx3 = struct.unpack_from('=iiii', raw, indcs)
    assert ng > 0 and min(nx1,nx2,nx3) > 0, 'unexpected RegionIndcs layout'
    start = indcs + 10*4
    stop = start + 9*4
    ranges = [[start, stop]]
    forcing = None
    if normalize_dormant_startup_forcing:
        extra, forcing = _dormant_forcing_ranges(raw, header, nmb)
        ranges.extend(extra)
    normalized = bytearray(raw)
    for lo, hi in ranges:
        normalized[lo:hi] = bytes(hi-lo)
    normalized = bytes(normalized)
    metadata = {
        'rule': 'zero only unused mesh_indcs cnx1/cnx2/cnx3/cis/cie/cjs/cje/cks/cke',
        'source': 'src/outputs/restart.cpp serializes RegionIndcs; src/mesh/mesh.cpp initializes only its root fine indices',
        'real_bytes': 8,
        'int_bytes': 4,
        'binary_header_offset': header,
        'normalized_ranges': ranges,
        'normalized_byte_count': sum(hi-lo for lo,hi in ranges),
        'raw_sha256': hashlib.sha256(raw).hexdigest(),
        'normalized_sha256': hashlib.sha256(normalized).hexdigest(),
    }
    if forcing is not None:
        metadata['dormant_startup_forcing'] = forcing
        if forcing['applied']:
            metadata['rule'] += '; explicitly opted-in verified dormant initial forcing RNG members only'
    return normalized, metadata


def restart_hashes(path, *, normalize_dormant_startup_forcing=False):
    return normalized_restart(path,
        normalize_dormant_startup_forcing=normalize_dormant_startup_forcing)[1]
