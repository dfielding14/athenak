"""Legacy-named Weibel-like engineering proxy; not reproduction evidence."""

import glob
import logging
import os
import subprocess

import numpy as np
import scripts.utils.athena as athena
from scripts.particles.pic_analysis_utils import fit_exponential_growth_windowed

import sys
sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_weibel_growth_publication.athinput'
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


def _mode1_amplitude(dataset):
    by = np.asarray(dataset['bcc2'], dtype=float)
    by_x = np.mean(by, axis=(0, 1))
    by_fluc = by_x - np.mean(by_x)

    x1f = np.asarray(dataset['x1f'], dtype=float)
    x1v = np.asarray(dataset['x1v'], dtype=float)
    length = float(x1f[-1] - x1f[0])
    phase = np.exp(-2.0j * np.pi * (x1v - x1f[0]) / length)
    mode = np.sum(by_fluc * phase) / by_x.size
    return float(np.abs(mode))


def _load_mode_series(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.mhd_bcc.*.bin')
    files = sorted(glob.glob(pattern))
    if len(files) < 12:
        raise RuntimeError('Not enough output samples for growth fit')

    times = []
    amp = []
    for fname in files:
        data = bin_convert.read_binary_as_athdf(fname)
        times.append(float(data['Time']))
        amp.append(_mode1_amplitude(data))

    t = np.asarray(times, dtype=float)
    a = np.asarray(amp, dtype=float)

    floor = max(1.0e-30, 1.0e-4 * np.max(a))
    fit = fit_exponential_growth_windowed(
        t, a, min_points=6, floor=floor, min_growth_factor=2.0
    )

    ratio = float((np.max(a) + floor) / max(np.min(a) + floor, 1.0e-30))
    return {
        'time': t,
        'amp': a,
        'gamma': float(fit['gamma']),
        'r2': float(fit['r2']),
        'fit_tmin': float(fit['tmin']),
        'fit_tmax': float(fit['tmax']),
        'fit_growth_factor': float(fit['growth_factor']),
        'ratio': ratio,
    }


def _run_case(label, nproc, basename):
    args = ['job/basename=' + basename]
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


def _check_close(label, measured, expected, abs_tol, rel_tol):
    abs_err = abs(measured - expected)
    rel_err = abs_err / max(abs(expected), 1.0)
    logger.info('%s measured=% .8e expected=% .8e abs_err=% .8e rel_err=% .8e',
                label, measured, expected, abs_err, rel_err)
    return abs_err <= abs_tol or rel_err <= rel_tol


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    _remove_outputs('pic_weibel_pub_np1')
    _run_case('serial', 1, 'pic_weibel_pub_np1')
    _RESULTS['np1'] = _load_mode_series('pic_weibel_pub_np1')

    if _athena_mpi_enabled():
        _remove_outputs('pic_weibel_pub_np2')
        _run_case('mpi2', 2, 'pic_weibel_pub_np2')
        _RESULTS['np2'] = _load_mode_series('pic_weibel_pub_np2')


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    g1 = _RESULTS['np1']['gamma']
    r1 = _RESULTS['np1']['r2']
    ok = _check_lower('np1:gamma_lower_bound', g1, 1.0e-1) and ok
    ok = _check_lower('np1:r2_lower_bound', r1, 0.80) and ok
    ok = _check_lower(
        'np1:growth_factor', _RESULTS['np1']['fit_growth_factor'], 2.0
    ) and ok

    if 'np2' in _RESULTS:
        g2 = _RESULTS['np2']['gamma']
        r2 = _RESULTS['np2']['r2']
        ok = _check_lower('np2:gamma_lower_bound', g2, 1.0e-1) and ok
        ok = _check_lower('np2:r2_lower_bound', r2, 0.75) and ok
        amp_linf = float(np.max(np.abs(_RESULTS['np2']['amp'] - _RESULTS['np1']['amp'])))
        time_linf = float(np.max(
            np.abs(_RESULTS['np2']['time'] - _RESULTS['np1']['time'])))
        ok = _check_close('serial_vs_mpi2:time_linf', time_linf, 0.0,
                          1.0e-12, 0.0) and ok
        ok = _check_close('serial_vs_mpi2:amp_linf', amp_linf, 0.0,
                          1.0e-10, 0.0) and ok
        ok = _check_close('serial_vs_mpi2:gamma', g2, g1, 1.0e-10, 1.0e-10) and ok

    return ok
