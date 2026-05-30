import glob
import logging
import os
import subprocess

import numpy as np
import scripts.utils.athena as athena
from scripts.particles.pic_analysis_utils import fit_exponential_growth_windowed
from scripts.particles.pic_analysis_utils import polarization_split

import sys
sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_crsi_deltaf_proxy.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}


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
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.*.bin')
    for fname in glob.glob(pattern):
        os.remove(fname)


def _mode1_complex(dataset, field):
    values = np.asarray(dataset[field], dtype=float)
    x_mode = np.mean(values, axis=(0, 1))
    fluc = x_mode - np.mean(x_mode)

    x1f = np.asarray(dataset['x1f'], dtype=float)
    x1v = np.asarray(dataset['x1v'], dtype=float)
    length = float(x1f[-1] - x1f[0])
    phase = np.exp(-2.0j * np.pi * (x1v - x1f[0]) / length)
    return np.sum(fluc * phase) / x_mode.size


def _fit_branch(time, amp):
    floor = max(1.0e-30, 1.0e-4 * float(np.max(amp)))
    fit = fit_exponential_growth_windowed(
        time, amp, min_points=6, floor=floor, min_growth_factor=1.2
    )
    ratio = float((amp[-1] + floor) / max(amp[0] + floor, 1.0e-30))
    return float(fit['gamma']), float(fit['r2']), float(ratio)


def _load_growth_metrics(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.mhd_bcc.*.bin')
    files = sorted(glob.glob(pattern))
    if len(files) < 8:
        raise RuntimeError('Not enough mhd_bcc outputs for CRSI proxy fit')

    times = []
    by_mode = []
    bz_mode = []
    for fname in files:
        data = bin_convert.read_binary_as_athdf(fname)
        times.append(float(data['Time']))
        by_mode.append(_mode1_complex(data, 'bcc2'))
        bz_mode.append(_mode1_complex(data, 'bcc3'))

    t = np.asarray(times, dtype=float)
    by = np.asarray(by_mode, dtype=complex)
    bz = np.asarray(bz_mode, dtype=complex)
    right, left = polarization_split(by, bz)

    amp_r = np.abs(right)
    amp_l = np.abs(left)

    gamma_r, r2_r, ratio_r = _fit_branch(t, amp_r)
    gamma_l, r2_l, ratio_l = _fit_branch(t, amp_l)

    dominant = 'right' if amp_r[-1] >= amp_l[-1] else 'left'
    if dominant == 'right':
        dom_gamma = gamma_r
        dom_r2 = r2_r
        dom_ratio = ratio_r
        sub_gamma = gamma_l
    else:
        dom_gamma = gamma_l
        dom_r2 = r2_l
        dom_ratio = ratio_l
        sub_gamma = gamma_r

    return {
        'dom_branch': dominant,
        'dom_gamma': float(dom_gamma),
        'dom_r2': float(dom_r2),
        'dom_ratio': float(dom_ratio),
        'sub_gamma': float(sub_gamma),
        'gamma_r': float(gamma_r),
        'gamma_l': float(gamma_l),
    }


def _run_case(label, basename, nproc, deltaf_on):
    args = [
        'job/basename=' + basename,
        'particles/pic_deltaf_mode=' + ('on' if deltaf_on else 'off'),
    ]
    if deltaf_on:
        args.append('particles/pic_deltaf_f0=kappa_iso')

    command = ['./athena', '-i', _athena_input_path()] + args
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command

    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed for ' + label + '\n' + output)


def _check_close(label, measured, expected, abs_tol, rel_tol):
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

    _remove_outputs('pic_crsi_deltaf_serial_off')
    _run_case('serial_off', 'pic_crsi_deltaf_serial_off', 1, False)
    _RESULTS['serial_off'] = _load_growth_metrics('pic_crsi_deltaf_serial_off')

    _remove_outputs('pic_crsi_deltaf_serial_on')
    _run_case('serial_on', 'pic_crsi_deltaf_serial_on', 1, True)
    _RESULTS['serial_on'] = _load_growth_metrics('pic_crsi_deltaf_serial_on')

    if _athena_mpi_enabled():
        _remove_outputs('pic_crsi_deltaf_mpi2_on')
        _run_case('mpi2_on', 'pic_crsi_deltaf_mpi2_on', 2, True)
        _RESULTS['mpi2_on'] = _load_growth_metrics('pic_crsi_deltaf_mpi2_on')

        _remove_outputs('pic_crsi_deltaf_mpi4_on')
        _run_case('mpi4_on', 'pic_crsi_deltaf_mpi4_on', 4, True)
        _RESULTS['mpi4_on'] = _load_growth_metrics('pic_crsi_deltaf_mpi4_on')


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    off = _RESULTS['serial_off']
    on = _RESULTS['serial_on']

    ok = _check_lower('serial_off:dom_ratio_lower', off['dom_ratio'], 10.0) and ok
    ok = _check_lower('serial_off:dom_gamma_lower', off['dom_gamma'], 2.0e-1) and ok
    ok = _check_lower('serial_off:dom_r2_lower', off['dom_r2'], 0.20) and ok
    ok = _check_lower('serial_on:dom_ratio_lower', on['dom_ratio'], 10.0) and ok
    ok = _check_lower('serial_on:dom_gamma_lower', on['dom_gamma'], 2.0e-1) and ok
    ok = _check_lower('serial_on:dom_r2_lower', on['dom_r2'], 0.20) and ok
    ok = _check_lower('serial_on:branch_split',
                      abs(on['dom_gamma'] - on['sub_gamma']), 1.0e-3) and ok

    if 'mpi2_on' in _RESULTS:
        mpi = _RESULTS['mpi2_on']
        ok = _check_lower('mpi2_on:dom_ratio_lower', mpi['dom_ratio'], 10.0) and ok
        ok = _check_lower('mpi2_on:dom_gamma_lower', mpi['dom_gamma'], 2.0e-1) and ok
        ok = _check_lower('mpi2_on:dom_r2_lower', mpi['dom_r2'], 0.20) and ok
        ok = _check_close('serial_vs_mpi2:on_dom_gamma',
                          mpi['dom_gamma'], on['dom_gamma'], 1.0e-8, 1.0e-8) and ok

    if 'mpi4_on' in _RESULTS:
        mpi = _RESULTS['mpi4_on']
        ok = _check_lower('mpi4_on:dom_ratio_lower', mpi['dom_ratio'], 10.0) and ok
        ok = _check_lower('mpi4_on:dom_gamma_lower', mpi['dom_gamma'], 2.0e-1) and ok
        ok = _check_lower('mpi4_on:dom_r2_lower', mpi['dom_r2'], 0.10) and ok
        ok = _check_close('serial_vs_mpi4:on_dom_gamma',
                          mpi['dom_gamma'], on['dom_gamma'], 1.0e-8, 1.0e-8) and ok

    if 'mpi2_on' in _RESULTS and 'mpi4_on' in _RESULTS:
        ok = _check_close('mpi2_vs_mpi4:on_dom_gamma',
                          _RESULTS['mpi4_on']['dom_gamma'],
                          _RESULTS['mpi2_on']['dom_gamma'], 1.0e-8, 1.0e-8) and ok

    return ok
