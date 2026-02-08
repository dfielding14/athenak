import glob
import logging
import os
import subprocess

import numpy as np
import scripts.utils.athena as athena
from scripts.particles.pic_analysis_utils import fit_exponential_growth

import sys
sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_two_stream_growth_proxy.athinput'
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
    rho = np.asarray(dataset['prtcl_rho'], dtype=float)
    rho_x = np.mean(rho, axis=(0, 1))
    rho_fluc = rho_x - np.mean(rho_x)

    x1f = np.asarray(dataset['x1f'], dtype=float)
    x1v = np.asarray(dataset['x1v'], dtype=float)
    length = float(x1f[-1] - x1f[0])
    phase = np.exp(-2.0j * np.pi * (x1v - x1f[0]) / length)
    mode = np.sum(rho_fluc * phase) / rho_x.size
    return float(np.abs(mode))


def _load_mode_series(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.prtcl_rho.*.bin')
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
    floor = max(1.0e-16, 1.0e-3 * np.max(a))
    gamma, _, r2 = fit_exponential_growth(t, a + floor, t[2], t[-3], floor=1.0e-30)
    ratio = float((np.max(a) + floor) / max(np.min(a) + floor, 1.0e-30))
    return {
        'time': t,
        'amp': a,
        'gamma': float(gamma),
        'r2': float(r2),
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


def _check_close(label, measured, expected, abs_tol, rel_tol):
    abs_err = abs(measured - expected)
    rel_err = abs_err / max(abs(expected), 1.0)
    logger.info('%s measured=% .8e expected=% .8e abs_err=% .8e rel_err=% .8e',
                label, measured, expected, abs_err, rel_err)
    return abs_err <= abs_tol or rel_err <= rel_tol


def _check_upper(label, measured, upper, tol):
    logger.info('%s measured=% .8e upper=% .8e margin=% .8e',
                label, measured, upper, measured - upper)
    return measured <= (upper + tol)


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    _remove_outputs('pic_two_stream_growth_proxy_np1')
    _run_case('serial', 1, 'pic_two_stream_growth_proxy_np1')
    _RESULTS['np1'] = _load_mode_series('pic_two_stream_growth_proxy_np1')

    if _athena_mpi_enabled():
        _remove_outputs('pic_two_stream_growth_proxy_np2')
        _run_case('mpi2', 2, 'pic_two_stream_growth_proxy_np2')
        _RESULTS['np2'] = _load_mode_series('pic_two_stream_growth_proxy_np2')


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    g1 = _RESULTS['np1']['gamma']
    ok = _check_upper('np1:gamma_upper_bound', g1, 0.0, 3.0e-2) and ok

    if 'np2' in _RESULTS:
        g2 = _RESULTS['np2']['gamma']
        ok = _check_upper('np2:gamma_upper_bound', g2, 0.0, 3.0e-2) and ok
        ok = _check_close('serial_vs_mpi2:gamma', g2, g1, 6.0e-2, 0.0) and ok

    return ok
