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

_INPUT_DECK = 'tests/pic_crsi_deltaf_publication.athinput'
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
    return {
        'gamma': float(fit['gamma']),
        'r2': float(fit['r2']),
        'tmin': float(fit['tmin']),
        'tmax': float(fit['tmax']),
        'growth_factor': float(fit['growth_factor']),
    }


def _load_growth_metrics(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.mhd_bcc.*.bin')
    files = sorted(glob.glob(pattern))
    if len(files) < 10:
        raise RuntimeError('Not enough mhd_bcc outputs for CRSI publication fit')

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

    fit_r = _fit_branch(t, amp_r)
    fit_l = _fit_branch(t, amp_l)

    dom = 'right' if amp_r[-1] >= amp_l[-1] else 'left'
    dom_gamma = fit_r['gamma'] if dom == 'right' else fit_l['gamma']

    return {
        'dom_branch': dom,
        'dom_gamma': float(dom_gamma),
        'gamma_r': float(fit_r['gamma']),
        'gamma_l': float(fit_l['gamma']),
        'r2_r': float(fit_r['r2']),
        'r2_l': float(fit_l['r2']),
        'growth_r': float(fit_r['growth_factor']),
        'growth_l': float(fit_l['growth_factor']),
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


def _check_lower(label, measured, lower):
    logger.info('%s measured=% .8e lower=% .8e margin=% .8e',
                label, measured, lower, measured - lower)
    return measured >= lower


def _check_upper(label, measured, upper):
    logger.info('%s measured=% .8e upper=% .8e margin=% .8e',
                label, measured, upper, measured - upper)
    return measured <= upper


def _check_log_ratio_close(label, measured, expected, tol_log):
    if measured <= 0.0 or expected <= 0.0:
        logger.error('%s requires positive values (measured=%e expected=%e)',
                     label, measured, expected)
        return False
    log_ratio = abs(np.log(measured / expected))
    logger.info('%s measured=% .8e expected=% .8e |log-ratio|=% .8e',
                label, measured, expected, log_ratio)
    return log_ratio <= tol_log


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    _remove_outputs('pic_crsi_pub_serial_off')
    _run_case('serial_off', 'pic_crsi_pub_serial_off', 1, False)
    _RESULTS['serial_off'] = _load_growth_metrics('pic_crsi_pub_serial_off')

    _remove_outputs('pic_crsi_pub_serial_on')
    _run_case('serial_on', 'pic_crsi_pub_serial_on', 1, True)
    _RESULTS['serial_on'] = _load_growth_metrics('pic_crsi_pub_serial_on')

    if _athena_mpi_enabled():
        _remove_outputs('pic_crsi_pub_mpi2_on')
        _run_case('mpi2_on', 'pic_crsi_pub_mpi2_on', 2, True)
        _RESULTS['mpi2_on'] = _load_growth_metrics('pic_crsi_pub_mpi2_on')

        _remove_outputs('pic_crsi_pub_mpi4_on')
        _run_case('mpi4_on', 'pic_crsi_pub_mpi4_on', 4, True)
        _RESULTS['mpi4_on'] = _load_growth_metrics('pic_crsi_pub_mpi4_on')


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    off = _RESULTS['serial_off']
    on = _RESULTS['serial_on']

    max_off = max(off['gamma_r'], off['gamma_l'])
    max_on = max(on['gamma_r'], on['gamma_l'])

    ok = _check_lower('serial_on:max_branch_gamma', max_on, 5.0e-2) and ok
    ok = _check_lower('serial_on:max_branch_r2', max(on['r2_r'], on['r2_l']), 0.6) and ok
    ok = _check_upper('serial:deltaf_noise_guard', max_off - max_on, 5.0) and ok

    if 'mpi2_on' in _RESULTS:
        mpi2 = _RESULTS['mpi2_on']
        max2 = max(mpi2['gamma_r'], mpi2['gamma_l'])
        ok = _check_lower('mpi2_on:max_branch_gamma', max2, 5.0e-2) and ok
        ok = _check_lower(
            'mpi2_on:max_branch_r2', max(mpi2['r2_r'], mpi2['r2_l']), 0.5
        ) and ok
        ok = _check_log_ratio_close('serial_vs_mpi2:max_branch_gamma',
                                    max2, max_on, 1.5) and ok

    if 'mpi4_on' in _RESULTS:
        mpi4 = _RESULTS['mpi4_on']
        max4 = max(mpi4['gamma_r'], mpi4['gamma_l'])
        ok = _check_lower('mpi4_on:max_branch_gamma', max4, 5.0e-2) and ok
        ok = _check_lower(
            'mpi4_on:max_branch_r2', max(mpi4['r2_r'], mpi4['r2_l']), 0.4
        ) and ok
        ok = _check_log_ratio_close('serial_vs_mpi4:max_branch_gamma',
                                    max4, max_on, 1.6) and ok

    return ok
