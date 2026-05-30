import glob
import json
import logging
import os
import re
import shutil
import struct
import subprocess
import tempfile

import numpy as np
import scripts.utils.athena as athena

import sys
sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_restart_safety_guards.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}
_PER_RANK_WATCH = {}
_RESTART_GUARDS = {}
_TRACK_GUARDS = {}
_PUSHER_GUARDS = {}
_EDGE_GUARDS = {}
_SKIPPED_DRILLS = {}
_PIC_RESTART_MAGIC = 0x5049435253543031
_REAL_BYTES = 8
_PSP_INDEX = 2
_MOMENT_COUNT_META_INDEX = 9
_MODEL_INT_COUNT = 31
_MODEL_REAL_COUNT = 37
_EXPECTED_RESTART_SCHEMA = 7

_CASES = {
    'no_mhd': [
        'particles/pic_background_mode=no_mhd',
        'particles/pic_feedback_mode=test_particle',
        'particles/couple_moments_to_mhd=false',
        'particles/couple_moments_momentum_to_mhd=false',
        'particles/couple_moments_energy_to_mhd=false',
    ],
    'passive_mhd': [
        'particles/pic_background_mode=passive_mhd',
        'particles/pic_feedback_mode=test_particle',
        'particles/couple_moments_to_mhd=false',
        'particles/couple_moments_momentum_to_mhd=false',
        'particles/couple_moments_energy_to_mhd=false',
    ],
    'coupled_edge_direct': [
        'particles/pic_background_mode=coupled',
        'particles/pic_feedback_mode=coupled',
        'particles/couple_moments_to_mhd=true',
        'particles/couple_moments_momentum_to_mhd=false',
        'particles/couple_moments_energy_to_mhd=false',
        'particles/couple_j_to_efield_representation=edge_staggered',
        'particles/couple_j_deposition_mode=direct_staggered',
    ],
}

_GUARD_CASES = [
    {
        'tag': 'guard_no_mhd_requires_test_particle',
        'args': [
            'particles/pic_background_mode=no_mhd',
            'particles/pic_feedback_mode=coupled',
            'particles/couple_moments_to_mhd=false',
            'particles/couple_moments_momentum_to_mhd=false',
            'particles/couple_moments_energy_to_mhd=false',
        ],
        'reason': '<particles>/pic_background_mode=no_mhd requires '
                  '<particles>/pic_feedback_mode=test_particle',
    },
    {
        'tag': 'guard_passive_mhd_requires_test_particle',
        'args': [
            'particles/pic_background_mode=passive_mhd',
            'particles/pic_feedback_mode=coupled',
            'particles/couple_moments_to_mhd=false',
            'particles/couple_moments_momentum_to_mhd=false',
            'particles/couple_moments_energy_to_mhd=false',
        ],
        'reason': '<particles>/pic_background_mode=passive_mhd requires '
                  '<particles>/pic_feedback_mode=test_particle unless coupled '
                  'feedback is explicitly implemented for this mode',
    },
    {
        'tag': 'guard_test_particle_rejects_coupling_toggles',
        'args': [
            'particles/pic_background_mode=coupled',
            'particles/pic_feedback_mode=test_particle',
            'particles/couple_moments_to_mhd=true',
            'particles/couple_moments_momentum_to_mhd=false',
            'particles/couple_moments_energy_to_mhd=false',
        ],
        'reason': '<particles>/pic_feedback_mode=test_particle does not support '
                  'particle-to-MHD coupling toggles',
    },
]

_PER_RANK_REASON_HINTS = [
    'coupled restart',
    'particle restart',
    'coupled particle restart state is missing',
    'particle restart state is missing',
    'failed to read coupled restart',
    'failed to read particle restart',
    'restart metadata is inconsistent with coupled particle section layout',
    'restart metadata is inconsistent with particle section layout',
]


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path(input_deck=_INPUT_DECK):
    return '../../' + athena.athena_rel_path + 'inputs/' + input_deck


def _athena_mpi_enabled():
    proc = subprocess.run(['./athena', '-c'], cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Unable to query Athena configuration with -c')
    output = (proc.stdout or '') + (proc.stderr or '')
    return 'MPI parallelism:            ON' in output


def _remove_outputs(basename):
    exe_dir = _athena_exe_dir()
    for pattern in [
            os.path.join(exe_dir, 'bin', basename + '.*.bin'),
            os.path.join(exe_dir, 'rst', basename + '.*.rst*'),
            os.path.join(exe_dir, 'rst', 'rank_*', basename + '.*.rst*'),
            os.path.join(exe_dir, 'trk', basename + '.trk')]:
        for fname in glob.glob(pattern):
            os.remove(fname)


def _build_command(nproc, arguments, restart_file=None, input_deck=_INPUT_DECK):
    command = ['./athena']
    if restart_file is None:
        command += ['-i', _athena_input_path(input_deck)]
    else:
        command += ['-r', restart_file]
    command += list(arguments)
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command
    return command


def _execute(label, nproc, arguments, restart_file=None, input_deck=_INPUT_DECK,
             env=None, timeout=None):
    command = _build_command(nproc, arguments, restart_file=restart_file,
                             input_deck=input_deck)
    logger.info('Executing %s: %s', label, ' '.join(command))
    child_env = os.environ.copy()
    if env is not None:
        child_env.update(env)
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True, env=child_env,
                          timeout=timeout)
    output = (proc.stdout or '') + (proc.stderr or '')
    return proc.returncode, output


def _run_success(label, nproc, arguments, restart_file=None):
    code, output = _execute(label, nproc, arguments, restart_file=restart_file)
    if code != 0:
        raise RuntimeError('Command failed for ' + label + '\n' + output)
    return output


def _run_expect_fail(label, nproc, arguments, reason):
    code, output = _execute(label, nproc, arguments)
    if code == 0:
        raise RuntimeError('Expected failure for ' + label + ', but command passed')
    if reason not in output:
        raise RuntimeError('Unexpected failure reason for ' + label + '\n'
                           'Expected substring: ' + reason + '\n'
                           'Output:\n' + output)


def _latest_restart_path(basename, per_rank=False):
    if per_rank:
        pattern = os.path.join(_athena_exe_dir(), 'rst', 'rank_00000000',
                               basename + '.*.rst')
    else:
        pattern = os.path.join(_athena_exe_dir(), 'rst', basename + '.*.rst')
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise RuntimeError('No restart files found for pattern: ' + pattern)
    return matches[-1], os.path.relpath(matches[-1], _athena_exe_dir())


def _assert_incomplete_restart_only(basename, per_rank=False):
    if per_rank:
        directory = os.path.join(_athena_exe_dir(), 'rst', 'rank_00000000')
    else:
        directory = os.path.join(_athena_exe_dir(), 'rst')
    partials = glob.glob(os.path.join(directory, basename + '.*.rst.partial'))
    published = glob.glob(os.path.join(directory, basename + '.*.rst'))
    if not partials:
        raise RuntimeError('Expected incomplete restart residue for ' + basename)
    if published:
        raise RuntimeError('Incomplete restart unexpectedly published for ' + basename)


def _fault_injector_env(injector, fault):
    preload = injector
    if os.environ.get('LD_PRELOAD'):
        preload += ':' + os.environ['LD_PRELOAD']
    return {
        'ATHENAK_RESTART_FAULT': fault,
        'LD_PRELOAD': preload,
    }


def _build_fault_injector(output_dir):
    source = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'restart_fault_injector.c')
    library = os.path.join(output_dir, 'restart_fault_injector.so')
    compiler = os.environ.get('CC', 'cc')
    command = [compiler, '-shared', '-fPIC', '-O2', '-o', library, source, '-ldl']
    proc = subprocess.run(command, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Unable to build restart fault injector\n' +
                           (proc.stdout or '') + (proc.stderr or ''))
    return library


def _skip_drill(name, reason):
    logger.warning('Skipping %s drill: %s', name, reason)
    _SKIPPED_DRILLS[name] = reason


def _find_particle_restart_section(data, restart_path):
    magic = struct.pack('<Q', _PIC_RESTART_MAGIC)
    section = data.find(magic)
    if section < 0:
        magic = struct.pack('>Q', _PIC_RESTART_MAGIC)
        section = data.find(magic)
    if section < 0:
        raise RuntimeError('Particle restart marker not found in ' + restart_path)
    return section


def _corrupt_particle_restart_metadata_int(restart_path, meta_index, bad_value):
    with open(restart_path, 'rb') as fp:
        data = bytearray(fp.read())

    section = _find_particle_restart_section(data, restart_path)
    offset = section + struct.calcsize('<Q') + meta_index*struct.calcsize('<i')
    struct.pack_into('<i', data, offset, bad_value)

    with open(restart_path, 'wb') as fp:
        fp.write(data)


def _fnv1a64(path):
    value = 14695981039346656037
    with open(path, 'rb') as fp:
        while True:
            payload = fp.read(1024 * 1024)
            if not payload:
                break
            for byte in payload:
                value ^= byte
                value = (value * 1099511628211) & 0xffffffffffffffff
    return value


def _write_completion_marker(path):
    with open(path + '.complete', 'w', encoding='ascii') as fp:
        fp.write('ATHENAK_RESTART_COMPLETE_V1\n')
        fp.write('size=' + str(os.path.getsize(path)) + '\n')
        fp.write('fnv1a64=' + format(_fnv1a64(path), '016x') + '\n')


def _write_shared_restart_publication(path):
    _write_completion_marker(path)
    manifest_path = path + '.manifest'
    manifest = {
        'schema': 'ATHENAK_RESTART_MANIFEST_V1',
        'members': [{
            'path': os.path.relpath(path, _athena_exe_dir()),
            'size': os.path.getsize(path),
            'fnv1a64': format(_fnv1a64(path), '016x'),
        }],
    }
    with open(manifest_path, 'w', encoding='ascii') as fp:
        json.dump(manifest, fp, indent=2, sort_keys=True)
        fp.write('\n')
    _write_completion_marker(manifest_path)


def _corrupt_first_particle_species(restart_path, bad_species):
    with open(restart_path, 'rb') as fp:
        data = bytearray(fp.read())

    section = _find_particle_restart_section(data, restart_path)

    offset = section + struct.calcsize('<Q')
    meta_fmt = '<15i'
    (version, nmb_section, nrdata, nidata, _nout1, _nout2, _nout3,
     _has_moments, _has_edge, _moment_cnt, _edge1_cnt, _edge2_cnt,
     _edge3_cnt, _state_kind, _physical_mode) = struct.unpack_from(
         meta_fmt, data, offset)
    offset += struct.calcsize(meta_fmt)
    offset += _REAL_BYTES  # cr_light_speed
    offset += _MODEL_INT_COUNT * struct.calcsize('<i')
    offset += _MODEL_REAL_COUNT * _REAL_BYTES
    npart_section = struct.unpack_from('<Q', data, offset)[0]
    offset += struct.calcsize('<Q')

    if (version != _EXPECTED_RESTART_SCHEMA or nmb_section <= 0 or nrdata <= 0
            or nidata <= _PSP_INDEX):
        raise RuntimeError('Unexpected particle restart metadata in ' + restart_path)
    if npart_section <= 0:
        raise RuntimeError('Cannot corrupt species in an empty particle restart')

    pr_real_offset = offset + nmb_section*struct.calcsize('<i')
    pr_int_offset = pr_real_offset + npart_section*nrdata*_REAL_BYTES
    first_species_offset = pr_int_offset + _PSP_INDEX*struct.calcsize('<i')
    struct.pack_into('<i', data, first_species_offset, bad_species)

    with open(restart_path, 'wb') as fp:
        fp.write(data)


def _run_direct_inflow_edge_current_bc_guard():
    base = 'pic_rst_safe_direct_inflow_edge_bc'
    _remove_outputs(base)

    args = [
        'job/basename=' + base,
        'time/nlim=2',
        'mesh/ix2_bc=inflow',
        'mesh/ox2_bc=inflow',
        'meshblock/nx1=16',
        'meshblock/nx2=8',
        'meshblock/nx3=8',
        'particles/pusher=drift',
        'particles/deposit_order=2',
        'particles/cr_vx0=0.0',
        'particles/cr_vy0=0.25',
        'particles/cr_vz0=0.125',
    ] + _CASES['coupled_edge_direct']
    expected = ('<particles>/couple_j_deposition_mode=direct_staggered does not '
                'support mesh/ix2_bc=inflow')
    _run_expect_fail('guard_direct_staggered_inflow_boundary', 1, args, expected)
    _EDGE_GUARDS['direct_inflow_rejected'] = True


def _run_corrupt_moment_count_restart_guard():
    base = 'pic_rst_safe_guard_bad_moment_count'
    base_seg = base + '_seg'
    base_rst = base + '_rst'
    base_bad = base + '_corrupt'

    for name in [base_seg, base_rst, base_bad]:
        _remove_outputs(name)

    case_args = _CASES['no_mhd']
    _run_success(base + '_seg_run', 1,
                 ['job/basename=' + base_seg, 'time/nlim=1'] + case_args)

    src_rst_path = os.path.join('rst', base_seg + '.00000.rst')
    dst_rst_path = os.path.join('rst', base_bad + '.00000.rst')
    full_src = os.path.join(_athena_exe_dir(), src_rst_path)
    full_dst = os.path.join(_athena_exe_dir(), dst_rst_path)
    if not os.path.exists(full_src):
        raise RuntimeError('Expected restart file not found: ' + full_src)
    shutil.copyfile(full_src, full_dst)
    _corrupt_particle_restart_metadata_int(full_dst, _MOMENT_COUNT_META_INDEX, 123456)
    _write_shared_restart_publication(full_dst)

    code, output = _execute(base + '_restart_run', 1,
                            ['job/basename=' + base_rst, 'time/nlim=2'] + case_args,
                            restart_file=dst_rst_path)
    expected = 'Particle restart moment size mismatch'
    if code == 0:
        raise RuntimeError('Expected corrupted restart failure, but command passed')
    if expected not in output:
        raise RuntimeError('Unexpected corrupted restart failure reason\n'
                           'Expected substring: ' + expected + '\n'
                           'Output:\n' + output)
    _RESTART_GUARDS['bad_moment_count'] = True


def _run_schema_version_restart_guard():
    base = 'pic_rst_safe_guard_bad_schema'
    base_seg = base + '_seg'
    base_rst = base + '_rst'
    base_bad = base + '_corrupt'

    for name in [base_seg, base_rst, base_bad]:
        _remove_outputs(name)

    case_args = _CASES['no_mhd']
    _run_success(base + '_seg_run', 1,
                 ['job/basename=' + base_seg, 'time/nlim=1'] + case_args)

    src_rst_path = os.path.join('rst', base_seg + '.00000.rst')
    dst_rst_path = os.path.join('rst', base_bad + '.00000.rst')
    full_src = os.path.join(_athena_exe_dir(), src_rst_path)
    full_dst = os.path.join(_athena_exe_dir(), dst_rst_path)
    if not os.path.exists(full_src):
        raise RuntimeError('Expected restart file not found: ' + full_src)
    shutil.copyfile(full_src, full_dst)
    _corrupt_particle_restart_metadata_int(
        full_dst, meta_index=0, bad_value=_EXPECTED_RESTART_SCHEMA + 1)
    _write_shared_restart_publication(full_dst)

    code, output = _execute(base + '_restart_run', 1,
                            ['job/basename=' + base_rst, 'time/nlim=2'] + case_args,
                            restart_file=dst_rst_path)
    expected = 'Unsupported particle restart version'
    if code == 0:
        raise RuntimeError('Expected schema-version restart failure, but command passed')
    if expected not in output:
        raise RuntimeError('Unexpected schema-version restart failure reason\n'
                           'Expected substring: ' + expected + '\n'
                           'Output:\n' + output)
    _RESTART_GUARDS['bad_schema_version'] = True


def _run_checksum_restart_guard():
    base = 'pic_rst_safe_guard_bad_checksum'
    base_seg = base + '_seg'
    base_rst = base + '_rst'
    base_bad = base + '_corrupt'

    for name in [base_seg, base_rst, base_bad]:
        _remove_outputs(name)

    case_args = _CASES['no_mhd']
    _run_success(base + '_seg_run', 1,
                 ['job/basename=' + base_seg, 'time/nlim=1'] + case_args)

    src_rst_path = os.path.join('rst', base_seg + '.00000.rst')
    dst_rst_path = os.path.join('rst', base_bad + '.00000.rst')
    full_src = os.path.join(_athena_exe_dir(), src_rst_path)
    full_dst = os.path.join(_athena_exe_dir(), dst_rst_path)
    if not os.path.exists(full_src):
        raise RuntimeError('Expected restart file not found: ' + full_src)
    shutil.copyfile(full_src, full_dst)
    _write_shared_restart_publication(full_dst)
    with open(full_dst, 'ab') as fp:
        fp.write(b'corrupt-checksum-fixture')

    code, output = _execute(base + '_restart_run', 1,
                            ['job/basename=' + base_rst, 'time/nlim=2'] + case_args,
                            restart_file=dst_rst_path)
    expected = 'restart checksum mismatch for completed artifact'
    if code == 0:
        raise RuntimeError('Expected corrupted restart failure, but command passed')
    if expected not in output:
        raise RuntimeError('Unexpected corrupted restart failure reason\\n'
                           'Expected substring: ' + expected + '\\n'
                           'Output:\\n' + output)
    _RESTART_GUARDS['bad_checksum'] = True


def _run_bound_restart_override_guard():
    base = 'pic_rst_safe_guard_bound_jcoef'
    base_seg = base + '_seg'
    base_rst = base + '_rst'

    for name in [base_seg, base_rst]:
        _remove_outputs(name)

    case_args = _CASES['coupled_edge_direct']
    _run_success(base + '_seg_run', 1,
                 ['job/basename=' + base_seg, 'time/nlim=1'] + case_args)

    rst_path = os.path.join('rst', base_seg + '.00000.rst')
    full_rst_path = os.path.join(_athena_exe_dir(), rst_path)
    if not os.path.exists(full_rst_path):
        raise RuntimeError('Expected restart file not found: ' + full_rst_path)

    code, output = _execute(
        base + '_restart_run', 1,
        ['job/basename=' + base_rst,
         'time/nlim=2'] + case_args + ['particles/couple_j_to_efield_coeff=2.0'],
        restart_file=rst_path)
    expected = 'Particle restart physical-model metadata mismatch'
    if code == 0:
        raise RuntimeError('Expected bound restart override failure, but command passed')
    if expected not in output:
        raise RuntimeError('Unexpected bound restart override failure reason\n'
                           'Expected substring: ' + expected + '\n'
                           'Output:\n' + output)
    _RESTART_GUARDS['bound_j_to_efield_coeff_override'] = True


def _copy_shared_restart_publication(src, dst):
    shutil.copyfile(src, dst)
    _write_shared_restart_publication(dst)


def _rewrite_manifest(path, update):
    manifest_path = path + '.manifest'
    with open(manifest_path, encoding='ascii') as fp:
        manifest = json.load(fp)
    update(manifest)
    with open(manifest_path, 'w', encoding='ascii') as fp:
        json.dump(manifest, fp, indent=2, sort_keys=True)
        fp.write('\n')
    _write_completion_marker(manifest_path)


def _run_restart_expect_fail(label, restart_file, arguments, expected):
    code, output = _execute(label, 1, arguments, restart_file=restart_file)
    if code == 0:
        raise RuntimeError('Expected restart publication failure, but command passed')
    if expected not in output:
        raise RuntimeError('Unexpected restart publication failure reason\n'
                           'Expected substring: ' + expected + '\n'
                           'Output:\n' + output)


def _run_publication_failure_path_guards():
    base = 'pic_rst_safe_guard_publication'
    base_seed = base + '_seed'
    case_args = _CASES['no_mhd']
    _remove_outputs(base_seed)
    _run_success(base + '_seed_run', 1,
                 ['job/basename=' + base_seed, 'time/nlim=1'] + case_args)

    src_rst_path = os.path.join('rst', base_seed + '.00000.rst')
    full_src = os.path.join(_athena_exe_dir(), src_rst_path)
    if not os.path.exists(full_src):
        raise RuntimeError('Expected restart file not found: ' + full_src)

    def prepare_fixture(tag):
        fixture_base = base + '_' + tag
        _remove_outputs(fixture_base)
        dst_rst_path = os.path.join('rst', fixture_base + '.00000.rst')
        full_dst = os.path.join(_athena_exe_dir(), dst_rst_path)
        _copy_shared_restart_publication(full_src, full_dst)
        return fixture_base, dst_rst_path, full_dst

    fixture_base, rst_path, full_dst = prepare_fixture('missing_marker')
    os.remove(full_dst + '.complete')
    _run_restart_expect_fail(
        fixture_base, rst_path,
        ['job/basename=' + fixture_base + '_rst', 'time/nlim=2'] + case_args,
        'restart completion marker is missing')
    _RESTART_GUARDS['missing_completion_marker'] = True

    fixture_base, rst_path, full_dst = prepare_fixture('truncated_marker')
    with open(full_dst + '.complete', 'w', encoding='ascii') as fp:
        fp.write('ATHENAK_RESTART_COMPLETE_V1\nsize=')
    _run_restart_expect_fail(
        fixture_base, rst_path,
        ['job/basename=' + fixture_base + '_rst', 'time/nlim=2'] + case_args,
        'restart completion marker is malformed')
    _RESTART_GUARDS['truncated_completion_marker'] = True

    fixture_base, rst_path, full_dst = prepare_fixture('truncated_payload')
    os.truncate(full_dst, max(1, os.path.getsize(full_dst) // 2))
    _run_restart_expect_fail(
        fixture_base, rst_path,
        ['job/basename=' + fixture_base + '_rst', 'time/nlim=2'] + case_args,
        'restart checksum mismatch for completed artifact')
    _RESTART_GUARDS['truncated_payload'] = True

    fixture_base, rst_path, full_dst = prepare_fixture('missing_manifest_marker')
    os.remove(full_dst + '.manifest.complete')
    _run_restart_expect_fail(
        fixture_base, rst_path,
        ['job/basename=' + fixture_base + '_rst', 'time/nlim=2'] + case_args,
        'restart completion marker is missing')
    _RESTART_GUARDS['missing_manifest_completion_marker'] = True

    fixture_base, rst_path, full_dst = prepare_fixture('manifest_wrong_path')
    payload_name = os.path.basename(full_dst)
    _rewrite_manifest(
        full_dst,
        lambda manifest: manifest['members'][0].update({
            'path': os.path.join('rst', 'rank_00000000', payload_name),
        }))
    _run_restart_expect_fail(
        fixture_base, rst_path,
        ['job/basename=' + fixture_base + '_rst', 'time/nlim=2'] + case_args,
        'restart manifest does not bind requested artifact')
    _RESTART_GUARDS['manifest_wrong_path'] = True

    fixture_base, rst_path, full_dst = prepare_fixture('manifest_wrong_digest')
    _rewrite_manifest(
        full_dst,
        lambda manifest: manifest['members'][0].update({
            'fnv1a64': '0000000000000000',
        }))
    _run_restart_expect_fail(
        fixture_base, rst_path,
        ['job/basename=' + fixture_base + '_rst', 'time/nlim=2'] + case_args,
        'restart manifest digest mismatch for member')
    _RESTART_GUARDS['manifest_wrong_digest'] = True

    published = glob.glob(os.path.join(_athena_exe_dir(), 'rst',
                                       base_seed + '.*.rst'))
    sequence = max(int(path.rsplit('.', 2)[-2]) for path in published) + 1
    newer = os.path.join(_athena_exe_dir(), 'rst',
                         base_seed + '.' + format(sequence, '05d') + '.rst')
    shutil.copyfile(full_src, newer + '.partial')
    with open(newer + '.partial', 'ab') as fp:
        fp.write(b'interrupted-restart-publication')
    with open(newer + '.complete.partial', 'w', encoding='ascii') as fp:
        fp.write('ATHENAK_RESTART_COMPLETE_V1\n')
    with open(newer + '.manifest.partial', 'w', encoding='ascii') as fp:
        fp.write('{"schema": "ATHENAK_RESTART_MANIFEST_V1", "members": [')
    prior_digest = _fnv1a64(full_src)
    _run_success(
        base + '_interrupted_prior_restart', 1,
        ['job/basename=' + base + '_prior_rst', 'time/nlim=2'] + case_args,
        restart_file=src_rst_path)
    if _fnv1a64(full_src) != prior_digest:
        raise RuntimeError(
            'Prior completed restart changed after interrupted publication')
    if os.path.exists(newer):
        raise RuntimeError('Interrupted publication unexpectedly exposed a final restart')
    _RESTART_GUARDS['interrupted_publication_preserves_prior'] = True


def _run_preload_failure_guard(injector, fault, guard, expected):
    base = 'pic_rst_safe_guard_' + guard
    _remove_outputs(base)
    code, output = _execute(
        base, 1,
        ['job/basename=' + base,
         'time/nlim=1',
         'output7/single_file_per_rank=true'] + _CASES['no_mhd'],
        env=_fault_injector_env(injector, fault),
        timeout=30)
    if code == 0:
        raise RuntimeError('Expected injected restart failure for ' + guard)
    if expected not in output:
        raise RuntimeError('Unexpected injected restart failure reason for ' + guard +
                           '\nExpected substring: ' + expected + '\nOutput:\n' + output)
    _assert_incomplete_restart_only(base, per_rank=True)
    _RESTART_GUARDS[guard] = True


def _run_killed_writer_restart_guard(injector):
    base = 'pic_rst_safe_guard_killed_writer'
    base_seed = base + '_seed'
    base_recovered = base + '_recovered'
    for name in [base_seed, base_recovered]:
        _remove_outputs(name)

    common = ['output7/single_file_per_rank=true'] + _CASES['no_mhd']
    _run_success(base + '_seed_run', 1,
                 ['job/basename=' + base_seed, 'time/nlim=1'] + common)
    prior_full, prior_rst = _latest_restart_path(base_seed, per_rank=True)
    prior_digest = _fnv1a64(prior_full)
    prior_published = set(glob.glob(os.path.join(
        _athena_exe_dir(), 'rst', 'rank_00000000', base_seed + '.*.rst')))

    code, _ = _execute(
        base + '_injected_run', 1,
        ['job/basename=' + base_seed, 'time/nlim=2'] + common,
        restart_file=prior_rst,
        env=_fault_injector_env(injector, 'kill_writer'),
        timeout=30)
    if code == 0:
        raise RuntimeError('Expected killed restart writer, but command passed')

    published = set(glob.glob(os.path.join(
        _athena_exe_dir(), 'rst', 'rank_00000000', base_seed + '.*.rst')))
    partials = glob.glob(os.path.join(
        _athena_exe_dir(), 'rst', 'rank_00000000', base_seed + '.*.rst.partial'))
    if published != prior_published:
        raise RuntimeError('Killed writer unexpectedly changed published restart set')
    if not partials or not any(os.path.getsize(path) > 0 for path in partials):
        raise RuntimeError('Killed writer did not leave non-empty incomplete residue')
    if _fnv1a64(prior_full) != prior_digest:
        raise RuntimeError('Killed writer changed the prior completed restart')

    _run_success(
        base + '_prior_restart', 1,
        ['job/basename=' + base_recovered,
         'time/nlim=2',
         'output7/dcycle=0'] + common,
        restart_file=prior_rst)
    _RESTART_GUARDS['killed_writer_preserves_prior'] = True


def _run_preload_restart_guards():
    build_dir = tempfile.mkdtemp(prefix='pic_restart_fault_', dir=_athena_exe_dir())
    try:
        try:
            injector = _build_fault_injector(build_dir)
        except RuntimeError as err:
            reason = str(err)
            for guard in [
                    'short_header_write',
                    'stdio_fseek_failure',
                    'killed_writer_preserves_prior']:
                _skip_drill(guard, reason)
            return
        _run_preload_failure_guard(
            injector, 'short_header_write', 'short_header_write',
            'Failed to write restart header to partial artifact')
        _run_preload_failure_guard(
            injector, 'fseek_failure', 'stdio_fseek_failure',
            'Error seeking before writing data')
        _run_killed_writer_restart_guard(injector)
    finally:
        shutil.rmtree(build_dir)


def _run_full_device_restart_target_guard():
    guard = 'dev_full_publication_failure'
    if not os.path.exists('/dev/full'):
        _skip_drill(guard, '/dev/full is not available on this host')
        return

    base = 'pic_rst_safe_guard_dev_full'
    run_dir = tempfile.mkdtemp(prefix=base + '_', dir=_athena_exe_dir())
    rank_dir = os.path.join(run_dir, 'rst', 'rank_00000000')
    os.makedirs(rank_dir)
    partial = os.path.join(rank_dir, base + '.00000.rst.partial')
    os.symlink('/dev/full', partial)
    try:
        code, output = _execute(
            base, 1,
            ['-d', run_dir,
             'job/basename=' + base,
             'time/nlim=1',
             'output7/single_file_per_rank=true'] + _CASES['no_mhd'],
            timeout=30)
        if code == 0:
            raise RuntimeError('Expected /dev/full restart publication failure')
        expected = [
            'Failed to write restart header to partial artifact',
            'Error seeking before writing data',
            'Failed to sync or close restart partial artifact',
        ]
        if not any(reason in output for reason in expected):
            raise RuntimeError('Unexpected /dev/full restart failure reason\n'
                               'Output:\n' + output)
        if os.path.exists(partial[:-len('.partial')]):
            raise RuntimeError('/dev/full restart unexpectedly published final artifact')
        _RESTART_GUARDS[guard] = True
    finally:
        shutil.rmtree(run_dir)


def _run_unwritable_restart_target_guard():
    base = 'pic_rst_safe_guard_unwritable_target'
    run_dir = tempfile.mkdtemp(prefix=base + '_', dir=_athena_exe_dir())
    rst_dir = os.path.join(run_dir, 'rst')
    os.mkdir(rst_dir)
    os.chmod(rst_dir, 0o555)
    try:
        code, output = _execute(
            base, 1,
            ['-d', run_dir,
             'job/basename=' + base,
             'time/nlim=1'] + _CASES['no_mhd'])
        if code == 0:
            raise RuntimeError('Expected unwritable restart-target failure, but command '
                               'passed')
        if '.rst.partial' not in output or 'could not be opened' not in output:
            raise RuntimeError('Unexpected unwritable restart-target failure reason\n'
                               'Output:\n' + output)
        _RESTART_GUARDS['unwritable_target'] = True
    finally:
        os.chmod(rst_dir, 0o755)
        shutil.rmtree(run_dir)


def _run_corrupt_species_restart_guard():
    base = 'pic_rst_safe_guard_bad_species'
    base_seg = base + '_seg'
    base_rst = base + '_rst'
    base_bad = base + '_corrupt'

    for name in [base_seg, base_rst, base_bad]:
        _remove_outputs(name)

    case_args = _CASES['no_mhd']
    _run_success(base + '_seg_run', 1,
                 ['job/basename=' + base_seg, 'time/nlim=1'] + case_args)

    src_rst_path = os.path.join('rst', base_seg + '.00000.rst')
    dst_rst_path = os.path.join('rst', base_bad + '.00000.rst')
    full_src = os.path.join(_athena_exe_dir(), src_rst_path)
    full_dst = os.path.join(_athena_exe_dir(), dst_rst_path)
    if not os.path.exists(full_src):
        raise RuntimeError('Expected restart file not found: ' + full_src)
    shutil.copyfile(full_src, full_dst)
    _corrupt_first_particle_species(full_dst, bad_species=99)
    _write_shared_restart_publication(full_dst)

    code, output = _execute(base + '_restart_run', 1,
                            ['job/basename=' + base_rst, 'time/nlim=2'] + case_args,
                            restart_file=dst_rst_path)
    expected = 'Restarted cosmic-ray particle species is out of range'
    if code == 0:
        raise RuntimeError('Expected corrupted restart failure, but command passed')
    if expected not in output:
        raise RuntimeError('Unexpected corrupted restart failure reason\n'
                           'Expected substring: ' + expected + '\n'
                           'Output:\n' + output)
    _RESTART_GUARDS['bad_species'] = True


def _run_missing_tracked_particle_guard():
    base = 'pic_rst_safe_guard_missing_tracked_tag'
    _remove_outputs(base)

    case_args = _CASES['no_mhd']
    code, output = _execute(base, 1,
                            ['job/basename=' + base,
                             'time/nlim=1',
                             'output8/nparticles=1025'] + case_args)
    expected = 'Tracked-particle output expected exactly one particle with tag'
    if code == 0:
        raise RuntimeError(
            'Expected missing tracked-particle failure, but command passed')
    if expected not in output:
        raise RuntimeError('Unexpected missing tracked-particle failure reason\n'
                           'Expected substring: ' + expected + '\n'
                           'Output:\n' + output)
    _TRACK_GUARDS['missing_tag'] = True


def _run_tracked_output_requires_particles_guard():
    base = 'pic_rst_safe_guard_trk_without_particles'
    _remove_outputs(base)

    code, output = _execute(
        base, 1,
        ['job/basename=' + base,
         'output1/file_type=trk'],
        input_deck='tests/linear_wave_hydro.athinput')
    expected = 'Tracked-particle output requires an active <particles> block'
    if code == 0:
        raise RuntimeError(
            'Expected tracked-output/no-particles failure, but command passed')
    if expected not in output:
        raise RuntimeError('Unexpected tracked-output/no-particles failure reason\n'
                           'Expected substring: ' + expected + '\n'
                           'Output:\n' + output)
    _TRACK_GUARDS['requires_particles'] = True


def _write_star_particle_file():
    path = os.path.join(_athena_exe_dir(), 'pic_guard_star_particles.txt')
    with open(path, 'w') as fp:
        fp.write('# x y z vx vy vz t_create mass\n')
        fp.write('1.0 1.0 1.0 0.0 0.0 0.0 0.0 1.0\n')
    return path


def _run_pusher_type_guards():
    base = 'pic_rst_safe_guard_cosmic_rk4'
    _remove_outputs(base)
    expected = '<particles>/pusher=rk4_gravity requires <particles>/particle_type=star'
    _run_expect_fail('guard_cosmic_rk4_pusher', 1,
                     ['job/basename=' + base,
                      'time/nlim=0',
                      'particles/pusher=rk4_gravity'],
                     expected)
    _PUSHER_GUARDS['cosmic_rk4'] = True

    _write_star_particle_file()
    base = 'pic_rst_safe_guard_star_boris'
    _remove_outputs(base)
    expected = 'Boris pushers are incompatible with star particles'
    _run_expect_fail('guard_star_boris_pusher', 1,
                     ['job/basename=' + base,
                      'time/nlim=0',
                      'particles/particle_type=star',
                      'particles/star_particle_file=pic_guard_star_particles.txt',
                      'particles/pusher=boris_tsc'],
                     expected)
    _PUSHER_GUARDS['star_boris'] = True


def _run_star_gravity_restart_override_guard():
    _write_star_particle_file()
    base = 'pic_rst_safe_guard_star_gravity'
    base_seg = base + '_seg'
    base_rst = base + '_rst'
    for name in [base_seg, base_rst]:
        _remove_outputs(name)

    common = [
        'particles/particle_type=star',
        'particles/star_particle_file=pic_guard_star_particles.txt',
        'particles/pusher=rk4_gravity',
        'particles/deposit_moments=false',
        'particles/couple_moments_to_mhd=false',
        'particles/couple_moments_momentum_to_mhd=false',
        'particles/couple_moments_energy_to_mhd=false',
        'output1/dcycle=0',
        'output2/dcycle=0',
        'output3/dcycle=0',
        'output4/dcycle=0',
        'output5/dcycle=0',
        'output6/dcycle=0',
        'output8/dcycle=0',
    ]
    _run_success(base + '_segment', 1,
                 ['job/basename=' + base_seg, 'time/nlim=1'] + common)
    rst_path = os.path.join('rst', base_seg + '.00000.rst')
    full_rst_path = os.path.join(_athena_exe_dir(), rst_path)
    if not os.path.exists(full_rst_path):
        raise RuntimeError('Expected restart file not found: ' + full_rst_path)

    _run_success(base + '_unchanged_restart', 1,
                 ['job/basename=' + base_rst, 'time/nlim=2'] + common,
                 restart_file=rst_path)
    code, output = _execute(
        base + '_changed_restart', 1,
        ['job/basename=' + base_rst + '_changed', 'time/nlim=2',
         'potential/mass_gal=2.0'] + common,
        restart_file=rst_path)
    expected = 'Particle restart physical-model metadata mismatch'
    if code == 0:
        raise RuntimeError('Expected star-gravity restart override failure, but command '
                           'passed')
    if expected not in output:
        raise RuntimeError('Unexpected star-gravity restart override failure reason\n'
                           'Expected substring: ' + expected + '\n'
                           'Output:\n' + output)
    _RESTART_GUARDS['star_gravity_bound_override'] = True


def _latest_output_file(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           basename + '.' + file_id + '.*.bin')
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise RuntimeError('No output files found for pattern: ' + pattern)
    return matches[-1]


def _integrate_quantity(dataset, quantity):
    dx1 = np.diff(dataset['x1f'])
    dx2 = np.diff(dataset['x2f'])
    dx3 = np.diff(dataset['x3f'])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    return float(np.sum(dataset[quantity] * dvol))


def _l2_quantity(dataset, quantity):
    dx1 = np.diff(dataset['x1f'])
    dx2 = np.diff(dataset['x2f'])
    dx3 = np.diff(dataset['x3f'])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    return float(np.sqrt(np.sum(dataset[quantity] * dataset[quantity] * dvol)))


def _read_tracked_snapshot(basename, ntrack=16):
    path = os.path.join(_athena_exe_dir(), 'trk', basename + '.trk')
    if not os.path.exists(path):
        raise RuntimeError('Tracked-particle output not found: ' + path)

    with open(path, 'rb') as fp:
        contents = fp.read()

    marker = b'# AthenaK tracked particle data at time='
    header_start = contents.rfind(marker)
    if header_start < 0:
        raise RuntimeError('Tracked-particle header not found in ' + path)
    data_start = contents.find(b'\n', header_start)
    if data_start < 0:
        raise RuntimeError('Tracked-particle header is unterminated in ' + path)
    data_start += 1
    if contents[data_start:data_start + 2] == b' \n':
        data_start += 2

    nvals = 6 * ntrack
    payload = contents[data_start:data_start + 4 * nvals]
    if len(payload) != 4 * nvals:
        raise RuntimeError('Tracked-particle payload is truncated in ' + path)
    return np.frombuffer(payload, dtype='>f4').astype(np.float64).reshape(ntrack, 6)


def _measure_case(basename):
    rho_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_rho'))
    jx_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jx'))
    jy_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jy'))
    jz_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jz'))
    pdens_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_d'))
    bcc_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'mhd_bcc'))

    return {
        'Q': _integrate_quantity(rho_data, 'prtcl_rho'),
        'Jx': _integrate_quantity(jx_data, 'prtcl_jx'),
        'Jy': _integrate_quantity(jy_data, 'prtcl_jy'),
        'Jz': _integrate_quantity(jz_data, 'prtcl_jz'),
        'npart': float(np.sum(pdens_data['pdens'])),
        'rho_l2': _l2_quantity(rho_data, 'prtcl_rho'),
        'jx_l2': _l2_quantity(jx_data, 'prtcl_jx'),
        'jy_l2': _l2_quantity(jy_data, 'prtcl_jy'),
        'jz_l2': _l2_quantity(jz_data, 'prtcl_jz'),
        'pdens_l2': _l2_quantity(pdens_data, 'pdens'),
        'bcc1_l2': _l2_quantity(bcc_data, 'bcc1'),
        'bcc2_l2': _l2_quantity(bcc_data, 'bcc2'),
        'bcc3_l2': _l2_quantity(bcc_data, 'bcc3'),
        'tracked': _read_tracked_snapshot(basename),
    }


def _run_soft_wallclock_continuation_parity():
    guard = 'soft_wallclock_continuation_parity'
    base = 'pic_rst_safe_soft_wallclock'
    base_full = base + '_full'
    base_seg = base + '_seg'
    base_rst = base + '_rst'
    for name in [base_full, base_seg, base_rst]:
        _remove_outputs(name)

    coarse_outputs = ['output' + str(index) + '/dcycle=1000000'
                      for index in range(1, 9)]
    common = ['time/tlim=1000000',
              'time/ndiag=1000000'] + coarse_outputs + _CASES['no_mhd']
    code, output = _execute(
        base + '_segment', 1,
        ['-t', '00:00:01',
         'job/basename=' + base_seg,
         'time/nlim=100000'] + common,
        timeout=30)
    if code != 0:
        raise RuntimeError('Soft wallclock segment failed\n' + output)
    if 'Terminating on wall clock limit' not in output:
        _skip_drill(guard, 'soft -t run did not terminate on the wallclock limit')
        return

    cycles = re.findall(r'time=.* cycle=([0-9]+)', output)
    if not cycles:
        raise RuntimeError('Unable to parse soft wallclock termination cycle')
    target_cycle = int(cycles[-1]) + 2
    _, restart_path = _latest_restart_path(base_seg)

    for label, arguments, restart_file in [
            (base + '_full_run',
             ['job/basename=' + base_full,
              'time/nlim=' + str(target_cycle)] + common,
             None),
            (base + '_restart_run',
             ['job/basename=' + base_rst,
              'time/nlim=' + str(target_cycle)] + common,
             restart_path)]:
        code, output = _execute(label, 1, arguments, restart_file=restart_file,
                                timeout=30)
        if code != 0:
            raise RuntimeError('Command failed for ' + label + '\n' + output)

    _RESULTS[base] = {
        'full': _measure_case(base_full),
        'restart': _measure_case(base_rst),
    }
    _RESTART_GUARDS[guard] = True


def _run_restart_triplet(case_tag, nproc, case_args):
    base = 'pic_rst_safe_' + case_tag + '_np' + str(nproc)
    base_full = base + '_full'
    base_seg = base + '_seg'
    base_rst = base + '_rst'

    for name in [base_full, base_seg, base_rst]:
        _remove_outputs(name)

    _run_success(base + '_full_run', nproc,
                 ['job/basename=' + base_full, 'time/nlim=3'] + case_args)
    full_measured = _measure_case(base_full)

    _run_success(base + '_seg_run', nproc,
                 ['job/basename=' + base_seg, 'time/nlim=1'] + case_args)

    rst_path = os.path.join('rst', base_seg + '.00000.rst')
    full_rst_path = os.path.join(_athena_exe_dir(), rst_path)
    if not os.path.exists(full_rst_path):
        raise RuntimeError('Expected restart file not found: ' + full_rst_path)

    _run_success(base + '_restart_run', nproc,
                 ['job/basename=' + base_rst, 'time/nlim=3'] + case_args,
                 restart_file=rst_path)
    rst_measured = _measure_case(base_rst)

    _RESULTS[base] = {
        'full': full_measured,
        'restart': rst_measured,
    }


def _run_per_rank_watch(nproc, case_args):
    base = 'pic_rst_safe_per_rank_np' + str(nproc)
    base_full = base + '_full'
    base_seg = base + '_seg'
    base_rst = base + '_rst'

    for name in [base_full, base_seg, base_rst]:
        _remove_outputs(name)

    full_args = ['job/basename=' + base_full,
                 'time/nlim=3',
                 'output7/single_file_per_rank=true'] + case_args
    seg_args = ['job/basename=' + base_seg,
                'time/nlim=1',
                'output7/single_file_per_rank=true'] + case_args
    rst_args = ['job/basename=' + base_rst,
                'time/nlim=3'] + case_args

    _run_success(base + '_full_run', nproc, full_args)
    full_measured = _measure_case(base_full)

    _run_success(base + '_seg_run', nproc, seg_args)

    rst_path = os.path.join('rst', 'rank_00000000', base_seg + '.00000.rst')
    full_rst_path = os.path.join(_athena_exe_dir(), rst_path)
    if not os.path.exists(full_rst_path):
        raise RuntimeError('Expected restart file not found: ' + full_rst_path)

    code, output = _execute(base + '_restart_run', nproc,
                            rst_args, restart_file=rst_path)
    if code == 0:
        rst_measured = _measure_case(base_rst)
        _PER_RANK_WATCH['status'] = 'supported'
        _PER_RANK_WATCH['full'] = full_measured
        _PER_RANK_WATCH['restart'] = rst_measured
        return

    lower = output.lower()
    matched = any(hint in lower for hint in _PER_RANK_REASON_HINTS)
    _PER_RANK_WATCH['status'] = 'guarded' if matched else 'unexpected_failure'
    _PER_RANK_WATCH['reason'] = output


def _check_with_tolerance(label, measured, expected, abs_tol, rel_tol):
    abs_err = abs(measured - expected)
    rel_err = abs_err / max(abs(expected), 1.0)
    logger.info('%s measured=% .8e expected=% .8e abs_err=% .8e rel_err=% .8e',
                label, measured, expected, abs_err, rel_err)
    return abs_err <= abs_tol or rel_err <= rel_tol


def _check_lower(label, measured, lower):
    logger.info('%s measured=% .8e lower=% .8e margin=% .8e',
                label, measured, lower, measured - lower)
    return measured >= lower


def _check_tracked_snapshot(label, tracked):
    ok = True
    if not np.all(np.isfinite(tracked)):
        logger.error('%s tracked output contains non-finite values', label)
        return False

    pos = tracked[:, :3]
    vel = tracked[:, 3:]
    lower = np.array([0.0, 0.0, 0.0])
    upper = np.array([16.0, 8.0, 8.0])
    tol = 1.0e-4

    min_pos = np.min(pos, axis=0)
    max_pos = np.max(pos, axis=0)
    span_pos = np.ptp(pos, axis=0)
    min_pos_norm = float(np.min(np.linalg.norm(pos, axis=1)))
    min_vel_norm = float(np.min(np.linalg.norm(vel, axis=1)))
    max_vel_norm = float(np.max(np.linalg.norm(vel, axis=1)))
    logger.info('%s tracked pos_min=%s pos_max=%s pos_span=%s '
                'vel_norm_min=% .8e vel_norm_max=% .8e',
                label, min_pos, max_pos, span_pos, min_vel_norm, max_vel_norm)

    if np.any(min_pos < lower - tol) or np.any(max_pos > upper + tol):
        logger.error('%s tracked positions are outside the test mesh domain', label)
        ok = False
    if min_pos_norm <= 1.0e-3:
        logger.error(
            '%s tracked output contains empty or byte-swapped position rows', label)
        ok = False
    if min_vel_norm <= 1.0e-4 or max_vel_norm <= 1.0e-3:
        logger.error('%s tracked output lacks resolved particle velocities', label)
        ok = False
    if span_pos[0] <= 0.5 or span_pos[1] <= 0.5:
        logger.error('%s tracked output lacks the expected x/y particle spread', label)
        ok = False

    return ok


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    proc_list = [1]
    mpi_enabled = _athena_mpi_enabled()
    if mpi_enabled:
        proc_list.append(2)

    for nproc in proc_list:
        for case_tag, case_args in _CASES.items():
            _run_restart_triplet(case_tag, nproc, case_args)

    for guard in _GUARD_CASES:
        base = 'pic_rst_safe_' + guard['tag']
        _remove_outputs(base)
        _run_expect_fail(guard['tag'], 1,
                         ['job/basename=' + base] + guard['args'],
                         guard['reason'])
    _run_corrupt_species_restart_guard()
    _run_corrupt_moment_count_restart_guard()
    _run_schema_version_restart_guard()
    _run_checksum_restart_guard()
    _run_bound_restart_override_guard()
    _run_publication_failure_path_guards()
    _run_preload_restart_guards()
    _run_full_device_restart_target_guard()
    _run_unwritable_restart_target_guard()
    _run_soft_wallclock_continuation_parity()
    _run_missing_tracked_particle_guard()
    _run_tracked_output_requires_particles_guard()
    _run_pusher_type_guards()
    _run_star_gravity_restart_override_guard()
    _run_direct_inflow_edge_current_bc_guard()

    if mpi_enabled:
        _run_per_rank_watch(2, _CASES['coupled_edge_direct'])


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    abs_tol = 1.0e-8
    rel_tol = 1.0e-8
    quantities = [
        'Q', 'Jx', 'Jy', 'Jz', 'npart',
        'rho_l2', 'jx_l2', 'jy_l2', 'jz_l2', 'pdens_l2',
        'bcc1_l2', 'bcc2_l2', 'bcc3_l2'
    ]

    for case_name, case_data in _RESULTS.items():
        full = case_data['full']
        rst = case_data['restart']
        for quantity in quantities:
            ok = _check_with_tolerance(case_name + ':' + quantity,
                                       rst[quantity], full[quantity],
                                       abs_tol, rel_tol) and ok

        ok = _check_tracked_snapshot(case_name + ':full', full['tracked']) and ok
        ok = _check_tracked_snapshot(case_name + ':restart', rst['tracked']) and ok
        track_err = float(np.max(np.abs(rst['tracked'] - full['tracked'])))
        logger.info('%s:tracked_particle_max_abs_err=% .8e', case_name, track_err)
        ok = (track_err <= 1.0e-5) and ok

        signal = abs(full['Jx']) + abs(full['Jy']) + abs(full['Jz'])
        ok = _check_lower(case_name + ':signal_nonzero', signal, 1.0e-10) and ok

    ok = bool(_RESTART_GUARDS.get('bad_species', False)) and ok
    if not _RESTART_GUARDS.get('bad_species', False):
        logger.error('Missing corrupted restart species guard result')
    ok = bool(_RESTART_GUARDS.get('bad_moment_count', False)) and ok
    if not _RESTART_GUARDS.get('bad_moment_count', False):
        logger.error('Missing corrupted restart moment-count guard result')
    ok = bool(_RESTART_GUARDS.get('bad_schema_version', False)) and ok
    if not _RESTART_GUARDS.get('bad_schema_version', False):
        logger.error('Missing corrupted restart schema-version guard result')
    ok = bool(_RESTART_GUARDS.get('bad_checksum', False)) and ok
    if not _RESTART_GUARDS.get('bad_checksum', False):
        logger.error('Missing corrupted restart checksum guard result')
    ok = bool(_RESTART_GUARDS.get('bound_j_to_efield_coeff_override', False)) and ok
    if not _RESTART_GUARDS.get('bound_j_to_efield_coeff_override', False):
        logger.error('Missing bound restart J-to-E coefficient override guard result')
    publication_guards = [
        'missing_completion_marker',
        'truncated_completion_marker',
        'truncated_payload',
        'missing_manifest_completion_marker',
        'manifest_wrong_path',
        'manifest_wrong_digest',
        'interrupted_publication_preserves_prior',
        'short_header_write',
        'stdio_fseek_failure',
        'dev_full_publication_failure',
        'killed_writer_preserves_prior',
        'unwritable_target',
        'soft_wallclock_continuation_parity',
    ]
    for guard in publication_guards:
        if guard in _SKIPPED_DRILLS:
            logger.warning('Restart publication drill skipped: %s: %s',
                           guard, _SKIPPED_DRILLS[guard])
            continue
        ok = bool(_RESTART_GUARDS.get(guard, False)) and ok
        if not _RESTART_GUARDS.get(guard, False):
            logger.error('Missing restart publication guard result: %s', guard)
    ok = bool(_TRACK_GUARDS.get('missing_tag', False)) and ok
    if not _TRACK_GUARDS.get('missing_tag', False):
        logger.error('Missing tracked-particle missing-tag guard result')
    ok = bool(_TRACK_GUARDS.get('requires_particles', False)) and ok
    if not _TRACK_GUARDS.get('requires_particles', False):
        logger.error('Missing tracked-output requires-particles guard result')
    ok = bool(_PUSHER_GUARDS.get('cosmic_rk4', False)) and ok
    if not _PUSHER_GUARDS.get('cosmic_rk4', False):
        logger.error('Missing cosmic-ray/rk4 pusher guard result')
    ok = bool(_PUSHER_GUARDS.get('star_boris', False)) and ok
    if not _PUSHER_GUARDS.get('star_boris', False):
        logger.error('Missing star/Boris pusher guard result')
    ok = bool(_RESTART_GUARDS.get('star_gravity_bound_override', False)) and ok
    if not _RESTART_GUARDS.get('star_gravity_bound_override', False):
        logger.error('Missing bound star-gravity restart override guard result')
    ok = bool(_EDGE_GUARDS.get('direct_inflow_rejected', False)) and ok
    if not _EDGE_GUARDS.get('direct_inflow_rejected', False):
        logger.error('Missing direct edge-current inflow guard result')

    if _PER_RANK_WATCH:
        status = _PER_RANK_WATCH.get('status', 'missing')
        logger.info('per-rank restart watch status: %s', status)
        if status == 'supported':
            full = _PER_RANK_WATCH['full']
            rst = _PER_RANK_WATCH['restart']
            for quantity in quantities:
                ok = _check_with_tolerance('per_rank:' + quantity,
                                           rst[quantity], full[quantity],
                                           abs_tol, rel_tol) and ok
            ok = _check_tracked_snapshot('per_rank:full', full['tracked']) and ok
            ok = _check_tracked_snapshot('per_rank:restart', rst['tracked']) and ok
            track_err = float(np.max(np.abs(rst['tracked'] - full['tracked'])))
            logger.info('per_rank:tracked_particle_max_abs_err=% .8e', track_err)
            ok = (track_err <= 1.0e-5) and ok
        elif status == 'guarded':
            logger.info('per-rank restart remains guarded with known reason family')
            ok = True and ok
        else:
            logger.error('Unexpected per-rank restart failure:\n%s',
                         _PER_RANK_WATCH.get('reason', '<no output>'))
            ok = False

    return ok
