import logging
import os
import subprocess

import numpy as np
import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_em_vacuum_wave.athinput'
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
    if data.ndim == 1:
        row = data
    else:
        row = data[-1, :]

    if row.size < 5:
        raise RuntimeError('Unexpected error-file format in ' + err_file)
    return float(row[4])


def _check_ratio(label, coarse, fine, max_ratio):
    ratio = fine / max(coarse, 1.0e-300)
    logger.info('%s coarse=% .8e fine=% .8e ratio=% .8e',
                label, coarse, fine, ratio)
    return ratio <= max_ratio


def _check_close(label, measured, expected, abs_tol, rel_tol):
    abs_err = abs(measured - expected)
    rel_err = abs_err / max(abs(expected), 1.0)
    logger.info('%s measured=% .8e expected=% .8e abs_err=% .8e rel_err=% .8e',
                label, measured, expected, abs_err, rel_err)
    return abs_err <= abs_tol or rel_err <= rel_tol


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    resolutions = [16, 32]
    for res in resolutions:
        base = 'pic_em_vacuum_np1_n' + str(res)
        args = [
            'job/basename=' + base,
            'mesh/nx1=' + str(res),
            'mesh/nx2=' + str(res // 2),
            'mesh/nx3=' + str(res // 2),
            'meshblock/nx1=' + str(max(4, res // 4)),
            'meshblock/nx2=' + str(max(4, res // 4)),
            'meshblock/nx3=' + str(max(4, res // 4)),
            'time/tlim=1.0',
            'time/nlim=1000',
            'problem/amp=1.0e-6',
            'problem/wave_flag=1',
            'problem/vflow=0.0',
            'output1/dt=-1.0',
            'output2/dt=-1.0',
            'output3/dt=-1.0',
            'output4/dt=-1.0',
            'output5/dt=-1.0',
        ]
        _run_case('serial_n' + str(res), 1, args)
        _RESULTS['np1_n' + str(res)] = _read_l1_error(base)

    if _athena_mpi_enabled():
        for res in resolutions:
            base = 'pic_em_vacuum_np2_n' + str(res)
            args = [
                'job/basename=' + base,
                'mesh/nx1=' + str(res),
                'mesh/nx2=' + str(res // 2),
                'mesh/nx3=' + str(res // 2),
                'meshblock/nx1=' + str(max(4, res // 4)),
                'meshblock/nx2=' + str(max(4, res // 4)),
                'meshblock/nx3=' + str(max(4, res // 4)),
                'time/tlim=1.0',
                'time/nlim=1000',
                'problem/amp=1.0e-6',
                'problem/wave_flag=1',
                'problem/vflow=0.0',
                'output1/dt=-1.0',
                'output2/dt=-1.0',
                'output3/dt=-1.0',
                'output4/dt=-1.0',
                'output5/dt=-1.0',
            ]
            _run_case('mpi2_n' + str(res), 2, args)
            _RESULTS['np2_n' + str(res)] = _read_l1_error(base)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    ok = _check_ratio('serial_convergence',
                      _RESULTS['np1_n16'], _RESULTS['np1_n32'], 0.6) and ok

    if 'np2_n16' in _RESULTS and 'np2_n32' in _RESULTS:
        ok = _check_ratio('mpi2_convergence',
                          _RESULTS['np2_n16'], _RESULTS['np2_n32'], 0.6) and ok
        ok = _check_close('serial_vs_mpi2_n32',
                          _RESULTS['np2_n32'], _RESULTS['np1_n32'],
                          1.0e-12, 1.0e-6) and ok

    return ok
