import glob
import logging
import os
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
    'coupled particle restart state is missing',
    'failed to read coupled restart',
    'restart metadata is inconsistent with coupled particle section layout',
]


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


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
            os.path.join(exe_dir, 'rst', 'rank_*', basename + '.*.rst')]:
        for fname in glob.glob(pattern):
            os.remove(fname)


def _build_command(nproc, arguments, restart_file=None):
    command = ['./athena']
    if restart_file is None:
        command += ['-i', _athena_input_path()]
    else:
        command += ['-r', restart_file]
    command += list(arguments)
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command
    return command


def _execute(label, nproc, arguments, restart_file=None):
    command = _build_command(nproc, arguments, restart_file=restart_file)
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
        'bcc1_l2': _l2_quantity(bcc_data, 'bcc1'),
        'bcc2_l2': _l2_quantity(bcc_data, 'bcc2'),
        'bcc3_l2': _l2_quantity(bcc_data, 'bcc3'),
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

    if mpi_enabled:
        _run_per_rank_watch(2, _CASES['coupled_edge_direct'])


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    abs_tol = 1.0e-8
    rel_tol = 1.0e-8
    quantities = ['Q', 'Jx', 'Jy', 'Jz', 'npart', 'bcc1_l2', 'bcc2_l2', 'bcc3_l2']

    for case_name, case_data in _RESULTS.items():
        full = case_data['full']
        rst = case_data['restart']
        for quantity in quantities:
            ok = _check_with_tolerance(case_name + ':' + quantity,
                                       rst[quantity], full[quantity],
                                       abs_tol, rel_tol) and ok

        signal = abs(full['Jx']) + abs(full['Jy']) + abs(full['Jz'])
        ok = _check_lower(case_name + ':signal_nonzero', signal, 1.0e-10) and ok

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
        elif status == 'guarded':
            logger.info('per-rank restart remains guarded with known reason family')
            ok = True and ok
        else:
            logger.error('Unexpected per-rank restart failure:\n%s',
                         _PER_RANK_WATCH.get('reason', '<no output>'))
            ok = False

    return ok
