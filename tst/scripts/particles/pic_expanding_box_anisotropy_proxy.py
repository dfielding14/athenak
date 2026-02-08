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

_CASES = {
    'expanding': 'tests/pic_expanding_box_proxy.athinput',
    'compressing': 'tests/pic_compressing_box_proxy.athinput',
}
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path(input_rel):
    return '../../' + athena.athena_rel_path + 'inputs/' + input_rel


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
    pattern_x = os.path.join(_athena_exe_dir(), 'bin', basename + '.prtcl_jx.*.bin')
    files_x = sorted(glob.glob(pattern_x))
    if len(files_x) < 12:
        raise RuntimeError('Not enough outputs for anisotropy trend fit')

    times = []
    ratios = []
    jx_series = []
    jy_series = []

    for xfile in files_x:
        stem = os.path.basename(xfile)
        cycle = stem.split('.')[-2]
        yfile = os.path.join(_athena_exe_dir(), 'bin',
                             basename + '.prtcl_jy.' + cycle + '.bin')

        xdat = bin_convert.read_binary_as_athdf(xfile)
        ydat = bin_convert.read_binary_as_athdf(yfile)

        jx = _integrate_quantity(xdat, 'prtcl_jx')
        jy = _integrate_quantity(ydat, 'prtcl_jy')

        times.append(float(xdat['Time']))
        jx_series.append(jx)
        jy_series.append(jy)
        ratios.append(abs(jx) / max(abs(jy), 1.0e-12))

    t = np.asarray(times, dtype=float)
    r = np.asarray(ratios, dtype=float)

    coeff = np.polyfit(t, r, 1)
    slope = float(coeff[0])

    return {
        'slope': slope,
        'ratio_start': float(r[0]),
        'ratio_end': float(r[-1]),
        'jx_start': float(jx_series[0]),
        'jx_end': float(jx_series[-1]),
    }


def _run_case(label, nproc, input_rel, basename):
    args = ['job/basename=' + basename]
    command = ['./athena', '-i', _athena_input_path(input_rel)] + args
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


def _check_sign(label, value, expected_positive):
    sign_ok = (value > 0.0) if expected_positive else (value < 0.0)
    logger.info('%s value=% .8e expected_positive=%s',
                label, value, str(expected_positive))
    return sign_ok


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    for label, input_rel in _CASES.items():
        base = 'pic_box_' + label + '_np1'
        _remove_outputs(base)
        _run_case(label + '_serial', 1, input_rel, base)
        _RESULTS[label + '_np1'] = _load_series(base)

    if _athena_mpi_enabled():
        for label, input_rel in _CASES.items():
            base = 'pic_box_' + label + '_np2'
            _remove_outputs(base)
            _run_case(label + '_mpi2', 2, input_rel, base)
            _RESULTS[label + '_np2'] = _load_series(base)

        for label, input_rel in _CASES.items():
            base = 'pic_box_' + label + '_np4'
            _remove_outputs(base)
            _run_case(label + '_mpi4', 4, input_rel, base)
            _RESULTS[label + '_np4'] = _load_series(base)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    exp1 = _RESULTS['expanding_np1']
    cmp1 = _RESULTS['compressing_np1']

    ok = _check_upper('expanding_np1:slope_upper', exp1['slope'], -1.0e-4) and ok
    ok = _check_lower('compressing_np1:slope_lower', cmp1['slope'], 1.0e-4) and ok
    ok = _check_lower('np1:slope_separation', cmp1['slope'] - exp1['slope'],
                      2.0e-4) and ok

    if 'expanding_np2' in _RESULTS:
        exp2 = _RESULTS['expanding_np2']
        cmp2 = _RESULTS['compressing_np2']

        ok = _check_sign('expanding_np2:slope_sign', exp2['slope'], False) and ok
        ok = _check_sign('compressing_np2:slope_sign', cmp2['slope'], True) and ok

        ok = _check_upper('expanding_np1_vs_np2:slope_abs_diff',
                          abs(exp2['slope'] - exp1['slope']), 2.0e-2) and ok
        ok = _check_upper('compressing_np1_vs_np2:slope_abs_diff',
                          abs(cmp2['slope'] - cmp1['slope']), 2.0e-2) and ok

    if 'expanding_np4' in _RESULTS:
        exp4 = _RESULTS['expanding_np4']
        cmp4 = _RESULTS['compressing_np4']

        ok = _check_sign('expanding_np4:slope_sign', exp4['slope'], False) and ok
        ok = _check_sign('compressing_np4:slope_sign', cmp4['slope'], True) and ok

        ok = _check_upper('expanding_np1_vs_np4:slope_abs_diff',
                          abs(exp4['slope'] - exp1['slope']), 2.0e-2) and ok
        ok = _check_upper('compressing_np1_vs_np4:slope_abs_diff',
                          abs(cmp4['slope'] - cmp1['slope']), 2.0e-2) and ok

    if 'expanding_np2' in _RESULTS and 'expanding_np4' in _RESULTS:
        ok = _check_upper('expanding_np2_vs_np4:slope_abs_diff',
                          abs(_RESULTS['expanding_np4']['slope'] -
                              _RESULTS['expanding_np2']['slope']), 2.0e-2) and ok
        ok = _check_upper('compressing_np2_vs_np4:slope_abs_diff',
                          abs(_RESULTS['compressing_np4']['slope'] -
                              _RESULTS['compressing_np2']['slope']), 2.0e-2) and ok

    return ok
