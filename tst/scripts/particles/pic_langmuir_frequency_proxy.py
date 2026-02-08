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

_INPUT_DECK = 'tests/pic_langmuir_frequency_proxy.athinput'
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


def _integrate_quantity(dataset, quantity):
    dx1 = np.diff(dataset['x1f'])
    dx2 = np.diff(dataset['x2f'])
    dx3 = np.diff(dataset['x3f'])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    return float(np.sum(dataset[quantity] * dvol))


def _load_series(basename):
    pattern_jx = os.path.join(_athena_exe_dir(), 'bin',
                              basename + '.prtcl_jx.*.bin')
    files = sorted(glob.glob(pattern_jx))
    if len(files) < 8:
        raise RuntimeError('Not enough output samples for frequency fit')

    times = []
    jx_series = []
    jy_series = []

    for jx_file in files:
        stem = os.path.basename(jx_file)
        cycle = stem.split('.')[-2]
        jy_file = os.path.join(_athena_exe_dir(), 'bin',
                               basename + '.prtcl_jy.' + cycle + '.bin')

        jx_data = bin_convert.read_binary_as_athdf(jx_file)
        jy_data = bin_convert.read_binary_as_athdf(jy_file)

        jxt = _integrate_quantity(jx_data, 'prtcl_jx')
        jyt = _integrate_quantity(jy_data, 'prtcl_jy')

        times.append(float(jx_data['Time']))
        jx_series.append(jxt)
        jy_series.append(jyt)

    return np.array(times), np.array(jx_series), np.array(jy_series)


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


def _dominant_frequency_nonuniform(time, signal, fmin, fmax, nsample=4000):
    t = np.asarray(time, dtype=float)
    y = np.asarray(signal, dtype=float)
    if t.size != y.size or t.size < 8:
        raise RuntimeError('insufficient samples for frequency fit')

    y = y - np.mean(y)
    freqs = np.linspace(fmin, fmax, nsample)
    best_f = freqs[0]
    best_amp = -1.0

    for f in freqs:
        omega_t = 2.0 * np.pi * f * t
        design = np.column_stack((np.sin(omega_t), np.cos(omega_t)))
        coeff, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        amp = float(np.sqrt(coeff[0] * coeff[0] + coeff[1] * coeff[1]))
        if amp > best_amp:
            best_amp = amp
            best_f = float(f)

    return best_f, best_amp


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    _remove_outputs('pic_langmuir_freq_proxy_np1')
    _run_case('serial', 1, 'pic_langmuir_freq_proxy_np1')
    _RESULTS['np1'] = _load_series('pic_langmuir_freq_proxy_np1')

    if _athena_mpi_enabled():
        _remove_outputs('pic_langmuir_freq_proxy_np2')
        _run_case('mpi2', 2, 'pic_langmuir_freq_proxy_np2')
        _RESULTS['np2'] = _load_series('pic_langmuir_freq_proxy_np2')


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    f_expected = 1.0 / (2.0 * np.pi)

    t1, vx1, vy1 = _RESULTS['np1']
    fx1, ax1 = _dominant_frequency_nonuniform(t1, vx1, 0.05, 0.40)
    fy1, ay1 = _dominant_frequency_nonuniform(t1, vy1, 0.05, 0.40)

    ok = _check_close('np1:vx_freq', fx1, f_expected, 5.0e-3, 5.0e-2) and ok
    ok = _check_close('np1:vy_freq', fy1, f_expected, 5.0e-3, 5.0e-2) and ok
    ok = _check_close('np1:amp_ratio', ay1 / max(ax1, 1.0e-30), 1.0,
                      5.0e-2, 1.0e-1) and ok

    if 'np2' in _RESULTS:
        t2, vx2, vy2 = _RESULTS['np2']
        fx2, ax2 = _dominant_frequency_nonuniform(t2, vx2, 0.05, 0.40)
        fy2, ay2 = _dominant_frequency_nonuniform(t2, vy2, 0.05, 0.40)

        ok = _check_close('np2:vx_freq', fx2, f_expected,
                          5.0e-3, 5.0e-2) and ok
        ok = _check_close('np2:vy_freq', fy2, f_expected,
                          5.0e-3, 5.0e-2) and ok
        ok = _check_close('np2:amp_ratio', ay2 / max(ax2, 1.0e-30), 1.0,
                          5.0e-2, 1.0e-1) and ok

        ok = _check_close('serial_vs_mpi2:vx_freq', fx2, fx1,
                          1.0e-4, 1.0e-3) and ok
        ok = _check_close('serial_vs_mpi2:vy_freq', fy2, fy1,
                          1.0e-4, 1.0e-3) and ok

    return ok
