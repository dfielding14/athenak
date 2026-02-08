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

_INPUT_DECK = 'tests/pic_bell_growth_proxy.athinput'
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


def _load_growth_metrics(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.mhd_bcc.*.bin')
    files = sorted(glob.glob(pattern))
    if len(files) < 8:
        raise RuntimeError('Not enough mhd_bcc outputs for Bell fit')

    times = []
    bt = []
    for fname in files:
        data = bin_convert.read_binary_as_athdf(fname)
        bperp = np.sqrt(np.mean(data['bcc2'] * data['bcc2'] +
                                data['bcc3'] * data['bcc3']))
        times.append(float(data['Time']))
        bt.append(float(bperp))

    t = np.asarray(times, dtype=float)
    b = np.asarray(bt, dtype=float)
    gamma, _, _ = fit_exponential_growth(t, b, t[1], t[-1], floor=1.0e-30)
    return {
        'gamma': float(gamma),
        'ratio': float(b[-1] / max(b[0], 1.0e-30)),
    }


def _run_case(label, basename, nproc, coupled):
    args = [
        'job/basename=' + basename,
        'particles/couple_moments_to_mhd=' + ('true' if coupled else 'false'),
        'particles/couple_moments_momentum_to_mhd=false',
        'particles/couple_moments_energy_to_mhd=false',
    ]
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


def _check_upper(label, measured, upper):
    logger.info('%s measured=% .8e upper=% .8e margin=% .8e',
                label, measured, upper, measured - upper)
    return measured <= upper


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    _remove_outputs('pic_bell_proxy_serial_uncoupled')
    _run_case('serial_uncoupled', 'pic_bell_proxy_serial_uncoupled', 1, False)
    _RESULTS['serial_uncoupled'] = _load_growth_metrics('pic_bell_proxy_serial_uncoupled')

    _remove_outputs('pic_bell_proxy_serial_coupled')
    _run_case('serial_coupled', 'pic_bell_proxy_serial_coupled', 1, True)
    _RESULTS['serial_coupled'] = _load_growth_metrics('pic_bell_proxy_serial_coupled')

    if _athena_mpi_enabled():
        _remove_outputs('pic_bell_proxy_mpi2_uncoupled')
        _run_case('mpi2_uncoupled', 'pic_bell_proxy_mpi2_uncoupled', 2, False)
        _RESULTS['mpi2_uncoupled'] = _load_growth_metrics('pic_bell_proxy_mpi2_uncoupled')

        _remove_outputs('pic_bell_proxy_mpi2_coupled')
        _run_case('mpi2_coupled', 'pic_bell_proxy_mpi2_coupled', 2, True)
        _RESULTS['mpi2_coupled'] = _load_growth_metrics('pic_bell_proxy_mpi2_coupled')


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    su = _RESULTS['serial_uncoupled']
    sc = _RESULTS['serial_coupled']

    ok = _check_upper('serial_uncoupled:gamma_upper', su['gamma'], 1.0e-3) and ok
    ok = _check_lower('serial_coupled:gamma_lower', sc['gamma'], 4.0e-3) and ok
    ok = _check_lower('serial_coupled:ratio_lower', sc['ratio'], 1.05) and ok
    ok = _check_lower('serial:gamma_separation',
                      sc['gamma'] - su['gamma'], 4.0e-3) and ok

    if 'mpi2_uncoupled' in _RESULTS and 'mpi2_coupled' in _RESULTS:
        mu = _RESULTS['mpi2_uncoupled']
        mc = _RESULTS['mpi2_coupled']
        ok = _check_upper('mpi2_uncoupled:gamma_upper', mu['gamma'], 1.0e-3) and ok
        ok = _check_lower('mpi2_coupled:gamma_lower', mc['gamma'], 4.0e-3) and ok
        ok = _check_lower('mpi2_coupled:ratio_lower', mc['ratio'], 1.05) and ok
        ok = _check_lower('mpi2:gamma_separation',
                          mc['gamma'] - mu['gamma'], 4.0e-3) and ok

        ok = _check_close('serial_vs_mpi2:uncoupled_gamma',
                          mu['gamma'], su['gamma'], 1.0e-5, 0.0) and ok
        ok = _check_close('serial_vs_mpi2:coupled_gamma',
                          mc['gamma'], sc['gamma'], 1.0e-5, 0.0) and ok
        ok = _check_close('serial_vs_mpi2:coupled_ratio',
                          mc['ratio'], sc['ratio'], 1.0e-6, 0.0) and ok

    return ok
