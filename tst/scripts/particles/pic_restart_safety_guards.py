import glob
import logging
import os
import shutil
import struct
import subprocess

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
_PIC_RESTART_MAGIC = 0x5049435253543031
_REAL_BYTES = 8
_PSP_INDEX = 2
_MOMENT_COUNT_META_INDEX = 9

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
            os.path.join(exe_dir, 'rst', basename + '.*.rst'),
            os.path.join(exe_dir, 'rst', 'rank_*', basename + '.*.rst'),
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


def _execute(label, nproc, arguments, restart_file=None, input_deck=_INPUT_DECK):
    command = _build_command(nproc, arguments, restart_file=restart_file,
                             input_deck=input_deck)
    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
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


def _corrupt_first_particle_species(restart_path, bad_species):
    with open(restart_path, 'rb') as fp:
        data = bytearray(fp.read())

    section = _find_particle_restart_section(data, restart_path)

    offset = section + struct.calcsize('<Q')
    meta_fmt = '<13i'
    (version, nmb_section, nrdata, nidata, _nout1, _nout2, _nout3,
     _has_moments, _has_edge, _moment_cnt, _edge1_cnt, _edge2_cnt,
     _edge3_cnt) = struct.unpack_from(meta_fmt, data, offset)
    offset += struct.calcsize(meta_fmt)
    npart_section = struct.unpack_from('<Q', data, offset)[0]
    offset += struct.calcsize('<Q')

    if version != 1 or nmb_section <= 0 or nrdata <= 0 or nidata <= _PSP_INDEX:
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
    _run_missing_tracked_particle_guard()
    _run_tracked_output_requires_particles_guard()
    _run_pusher_type_guards()
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
