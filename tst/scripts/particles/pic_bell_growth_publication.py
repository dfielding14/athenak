import glob
import logging
import os
import subprocess

import numpy as np
import scripts.utils.athena as athena
from scripts.particles.pic_analysis_utils import fit_exponential_growth
from scripts.particles.pic_analysis_utils import fit_exponential_growth_windowed

import sys
sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_bell_growth_publication.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_COUPLED_GAMMA_MIN = 4.0
_COUPLED_FIT_GROWTH_MIN = 4.0
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


def _load_growth_metrics(basename, fit_coupled):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.mhd_bcc.*.bin')
    files = sorted(glob.glob(pattern))
    if len(files) < 10:
        raise RuntimeError('Not enough mhd_bcc outputs for Bell fit')

    times = []
    bmode = []
    for fname in files:
        data = bin_convert.read_binary_as_athdf(fname)
        b2_k = _mode1_complex(data, 'bcc2')
        b3_k = _mode1_complex(data, 'bcc3')
        bperp = np.sqrt(np.abs(b2_k) * np.abs(b2_k) + np.abs(b3_k) * np.abs(b3_k))
        times.append(float(data['Time']))
        bmode.append(float(bperp))

    t = np.asarray(times, dtype=float)
    b = np.asarray(bmode, dtype=float)
    floor = max(1.0e-30, 1.0e-6 * np.max(b))

    if fit_coupled:
        fit = fit_exponential_growth_windowed(
            t, b, min_points=8, floor=floor, min_growth_factor=2.0
        )
        gamma = float(fit['gamma'])
        r2 = float(fit['r2'])
    else:
        gamma, _, r2 = fit_exponential_growth(
            t, np.maximum(b, floor), float(t[0]), float(t[-1]), floor=floor
        )
        fit = {
            'gamma': float(gamma),
            'r2': float(r2),
            'tmin': float(t[1]),
            'tmax': float(t[-1]),
            'growth_factor': float(b[-1] / max(b[0], floor)),
        }

    return {
        'gamma': float(fit['gamma']),
        'r2': float(fit['r2']),
        'ratio': float(b[-1] / max(b[0], 1.0e-30)),
        'fit_tmin': float(fit['tmin']),
        'fit_tmax': float(fit['tmax']),
        'fit_growth_factor': float(fit['growth_factor']),
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


def _check_lower(label, measured, lower):
    logger.info('%s measured=% .8e lower=% .8e margin=% .8e',
                label, measured, lower, measured - lower)
    return measured >= lower


def _check_upper(label, measured, upper):
    logger.info('%s measured=% .8e upper=% .8e margin=% .8e',
                label, measured, upper, measured - upper)
    return measured <= upper


def _check_close(label, measured, expected, abs_tol, rel_tol):
    abs_err = abs(measured - expected)
    rel_err = abs_err / max(abs(expected), 1.0)
    logger.info('%s measured=% .8e expected=% .8e abs_err=% .8e rel_err=% .8e',
                label, measured, expected, abs_err, rel_err)
    return abs_err <= abs_tol or rel_err <= rel_tol


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    _remove_outputs('pic_bell_pub_serial_uncoupled')
    _run_case('serial_uncoupled', 'pic_bell_pub_serial_uncoupled', 1, False)
    _RESULTS['serial_uncoupled'] = _load_growth_metrics(
        'pic_bell_pub_serial_uncoupled', fit_coupled=False
    )

    _remove_outputs('pic_bell_pub_serial_coupled')
    _run_case('serial_coupled', 'pic_bell_pub_serial_coupled', 1, True)
    _RESULTS['serial_coupled'] = _load_growth_metrics(
        'pic_bell_pub_serial_coupled', fit_coupled=True
    )

    if _athena_mpi_enabled():
        _remove_outputs('pic_bell_pub_mpi2_uncoupled')
        _run_case('mpi2_uncoupled', 'pic_bell_pub_mpi2_uncoupled', 2, False)
        _RESULTS['mpi2_uncoupled'] = _load_growth_metrics(
            'pic_bell_pub_mpi2_uncoupled', fit_coupled=False
        )

        _remove_outputs('pic_bell_pub_mpi2_coupled')
        _run_case('mpi2_coupled', 'pic_bell_pub_mpi2_coupled', 2, True)
        _RESULTS['mpi2_coupled'] = _load_growth_metrics(
            'pic_bell_pub_mpi2_coupled', fit_coupled=True
        )


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    su = _RESULTS['serial_uncoupled']
    sc = _RESULTS['serial_coupled']

    ok = _check_upper('serial_uncoupled:gamma_upper', abs(su['gamma']), 2.0e-2) and ok
    ok = _check_lower(
        'serial_coupled:gamma_lower', sc['gamma'], _COUPLED_GAMMA_MIN
    ) and ok
    ok = _check_lower('serial_coupled:r2_lower', sc['r2'], 0.85) and ok
    ok = _check_lower(
        'serial_coupled:fit_growth_factor_lower',
        sc['fit_growth_factor'],
        _COUPLED_FIT_GROWTH_MIN,
    ) and ok

    if 'mpi2_uncoupled' in _RESULTS and 'mpi2_coupled' in _RESULTS:
        mu = _RESULTS['mpi2_uncoupled']
        mc = _RESULTS['mpi2_coupled']
        ok = _check_upper('mpi2_uncoupled:gamma_upper', abs(mu['gamma']), 2.0e-2) and ok
        ok = _check_lower(
            'mpi2_coupled:gamma_lower', mc['gamma'], _COUPLED_GAMMA_MIN
        ) and ok
        ok = _check_lower('mpi2_coupled:r2_lower', mc['r2'], 0.85) and ok
        ok = _check_lower(
            'mpi2_coupled:fit_growth_factor_lower',
            mc['fit_growth_factor'],
            _COUPLED_FIT_GROWTH_MIN,
        ) and ok
        ok = _check_close('serial_vs_mpi2:coupled_gamma',
                          mc['gamma'], sc['gamma'], 2.0, 0.5) and ok
        ok = _check_close(
            'serial_vs_mpi2:coupled_fit_growth_factor',
            mc['fit_growth_factor'], sc['fit_growth_factor'], 8.0, 0.6
        ) and ok

    return ok
