import logging
import os
import subprocess

import numpy as np
import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_em_vacuum_wave.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESOLUTIONS = [16, 24, 32, 48, 64]
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


def _run_case(label, nproc, arguments):
    command = ['./athena', '-i', _athena_input_path()] + list(arguments)
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command

    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed for ' + label + '\n' + output)


def _read_l1_error(basename):
    err_file = os.path.join(_athena_exe_dir(), basename + '-errs.dat')
    if not os.path.exists(err_file):
        raise RuntimeError('Missing error file: ' + err_file)

    data = np.loadtxt(err_file)
    row = data if data.ndim == 1 else data[-1, :]
    if row.size < 5:
        raise RuntimeError('Unexpected error-file format in ' + err_file)
    return float(row[4])


def _run_resolution(nproc, res):
    base = f'pic_em_vacuum_pub_np{nproc}_n{res}'
    args = [
        'job/basename=' + base,
        'mesh/nx1=' + str(res),
        'mesh/nx2=' + str(res // 2),
        'mesh/nx3=' + str(res // 2),
        'meshblock/nx1=' + str(max(4, res // 4)),
        'meshblock/nx2=' + str(max(4, res // 4)),
        'meshblock/nx3=' + str(max(4, res // 4)),
        'time/tlim=1.0',
        'time/nlim=1200',
        'problem/amp=1.0e-6',
        'problem/wave_flag=1',
        'problem/vflow=0.0',
        'output1/dt=-1.0',
        'output2/dt=-1.0',
        'output3/dt=-1.0',
        'output4/dt=-1.0',
        'output5/dt=-1.0',
    ]
    _run_case(f'np{nproc}_n{res}', nproc, args)
    return _read_l1_error(base)


def _fit_order(l1_by_res):
    x = np.log(np.asarray(sorted(l1_by_res), dtype=float))
    y = np.log(np.asarray([l1_by_res[r] for r in sorted(l1_by_res)], dtype=float))
    coeff = np.polyfit(x, y, 1)
    return float(-coeff[0])


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

    np1 = {}
    for res in _RESOLUTIONS:
        np1[res] = _run_resolution(1, res)
    _RESULTS['np1'] = np1
    _RESULTS['np1_order'] = _fit_order(np1)

    if _athena_mpi_enabled():
        np2 = {}
        for res in _RESOLUTIONS:
            np2[res] = _run_resolution(2, res)
        _RESULTS['np2'] = np2
        _RESULTS['np2_order'] = _fit_order(np2)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    np1 = _RESULTS['np1']
    ok = _check_lower('np1:order', _RESULTS['np1_order'], 1.4) and ok

    res_sorted = sorted(np1)
    for i in range(len(res_sorted) - 1):
        c = np1[res_sorted[i]]
        f = np1[res_sorted[i + 1]]
        ok = _check_upper(f'np1:ratio_n{res_sorted[i]}_to_n{res_sorted[i + 1]}',
                          f / max(c, 1.0e-300), 0.75) and ok

    if 'np2' in _RESULTS:
        np2 = _RESULTS['np2']
        ok = _check_lower('np2:order', _RESULTS['np2_order'], 1.4) and ok
        for res in res_sorted:
            ok = _check_close(f'np1_vs_np2:n{res}', np2[res], np1[res],
                              1.0e-12, 2.0e-6) and ok

    return ok
