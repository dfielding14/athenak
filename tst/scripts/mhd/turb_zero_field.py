"""Magnetic-field and thermal-energy initialization in the turb pgen."""

import glob
import logging
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena  # noqa: F401 (run_tests.py build entry)

sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT = 'tests/pic_turb_zero_field.athinput'
_MPIEXEC = shlex.split(os.environ.get('MPIEXEC', 'mpiexec'))
_RESULTS = {}
_SKIPPED = False


def _exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _configuration():
    proc = subprocess.run(
        ['./athena', '-c'], cwd=_exe_dir(), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Unable to query Athena configuration')
    return (proc.stdout or '') + (proc.stderr or '')


def _remove_outputs(basename):
    for pattern in (
            os.path.join(_exe_dir(), 'bin', basename + '.*.bin'),
            os.path.join(_exe_dir(), basename + '.*.hst')):
        for path in glob.glob(pattern):
            os.remove(path)


def _mpi_launcher_available():
    return bool(_MPIEXEC and shutil.which(_MPIEXEC[0]))


def _run(basename, nproc, ifield=0, twofluid=False):
    source = Path(__file__).resolve().parents[3] / 'inputs' / _INPUT
    text = source.read_text().replace('ifield = 0', 'ifield = ' + str(ifield))
    text = text.replace('variable = mhd_bcc', 'variable = mhd_w_bcc')
    if twofluid:
        text = text.replace('integrator = rk2', 'integrator = imex2')
        text += ('\n<hydro>\neos = ideal\ngamma = 1.6666666666666667\n'
                 'rsolver = llf\nreconstruct = plm\n'
                 '<ion-neutral>\ndrag_coeff = 0\n')
    input_path = str(Path(_exe_dir()) / (basename + '.athinput'))
    Path(input_path).write_text(text)
    command = ['./athena', '-i', input_path, 'job/basename=' + basename]
    if nproc > 1:
        command = _MPIEXEC + ['-n', str(nproc)] + command
    proc = subprocess.run(
        command, cwd=_exe_dir(), capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('ifield=0 run failed:\n' + output)


def _latest(basename, file_id):
    paths = sorted(glob.glob(os.path.join(
        _exe_dir(), 'bin', basename + '.' + file_id + '.*.bin')))
    if not paths:
        raise RuntimeError('Missing binary output for ' + basename)
    return paths[-1]


def _measure(basename, ifield=0, twofluid=False):
    magnetic = bin_convert.read_binary_as_athdf(
        _latest(basename, 'mhd_bcc'))
    divb = bin_convert.read_binary_as_athdf(
        _latest(basename, 'mhd_divb'))
    user_history = os.path.join(_exe_dir(), basename + '.user.hst')
    if not os.path.exists(user_history):
        raise RuntimeError('Missing user history output for ' + basename)
    with open(user_history, 'r', encoding='utf-8') as stream:
        history_text = stream.read()
    amplitude = 2.0 if twofluid else np.sqrt(2.0/(5.0/3.0))
    expected_bz = np.zeros_like(magnetic['bcc3'])
    if ifield == 1:
        expected_bz[:] = amplitude*np.sin(2.0*np.pi*magnetic['x1v'])
    elif ifield == 2:
        expected_bz[:] = amplitude
    return {
        'field_error': float(np.max(np.abs(np.stack([
            magnetic['bcc1'], magnetic['bcc2'], magnetic['bcc3'] - expected_bz])))),
        'zero_required': ifield == 0,
        'thermal_error': float(np.max(np.abs(magnetic['eint'] - 0.9))),
        'max_b': max(float(np.max(np.abs(magnetic[name])))
                     for name in ('bcc1', 'bcc2', 'bcc3')),
        'max_divb': float(np.max(np.abs(divb['divb']))),
        'history_ok': all(label in history_text for label in (
            'vx_vol', 'divB_max', 'cr_Ekin', 'cr_Pz')),
    }


def run(**kwargs):
    global _SKIPPED
    logger.debug('Running test ' + __name__)
    config = _configuration()
    if 'Problem generator:          turb' not in config:
        logger.info('Skipping: configure this test with -DPROBLEM=turb')
        _SKIPPED = True
        return

    cases = [('np1', 1)]
    if ('MPI parallelism:            ON' in config
            and _mpi_launcher_available()):
        cases.append(('np2', 2))
    elif 'MPI parallelism:            ON' in config:
        logger.info('Skipping MPI case: launcher %s not found',
                    _MPIEXEC[0] if _MPIEXEC else '<empty>')

    for name, nproc in cases:
        basename = 'TurbZeroField_' + name
        _remove_outputs(basename)
        _run(basename, nproc)
        _RESULTS[name] = _measure(basename)
        for ifield in (1, 2):
            for twofluid in (False, True):
                label = name + '_field' + str(ifield) + '_two' + str(twofluid)
                basename = 'TurbField_' + label
                _remove_outputs(basename)
                _run(basename, nproc, ifield, twofluid)
                _RESULTS[label] = _measure(basename, ifield, twofluid)


def analyze():
    if _SKIPPED:
        return True
    ok = True
    for name, result in _RESULTS.items():
        logger.info('%s max_b=% .8e max_divb=% .8e',
                    name, result['max_b'], result['max_divb'])
        logger.info('%s field error=% .8e thermal error=% .8e',
                    name, result['field_error'], result['thermal_error'])
        if result['zero_required']:
            ok = (result['max_b'] == 0.0) and ok
        ok = (result['field_error'] < 2.0e-7) and ok
        ok = (result['thermal_error'] < 2.0e-7) and ok
        ok = (result['max_divb'] == 0.0) and ok
        ok = result['history_ok'] and ok
    return ok
