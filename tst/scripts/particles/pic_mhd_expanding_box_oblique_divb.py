import glob
import logging
import os
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena

sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_mhd_expanding_box_oblique_divb.athinput'
_RATES = np.asarray([0.05, 0.10, 0.15])
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _remove_outputs():
    pattern = os.path.join(_athena_exe_dir(),
                           'bin/pic_mhd_expanding_box_oblique_divb.*.bin')
    for fname in glob.glob(pattern):
        os.remove(fname)


def _endpoints(file_id):
    pattern = os.path.join(
        _athena_exe_dir(),
        'bin/pic_mhd_expanding_box_oblique_divb.' + file_id + '.*.bin')
    paths = sorted(glob.glob(pattern))
    if len(paths) < 2:
        raise RuntimeError('Need initial and final snapshots for: ' + pattern)
    return (bin_convert.read_binary_as_athdf(paths[0]),
            bin_convert.read_binary_as_athdf(paths[-1]))


def _relative_error(actual, expected):
    scale = max(float(np.max(np.abs(expected))), 1.0e-14)
    return float(np.max(np.abs(actual - expected))) / scale


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _RESULTS.clear()
    _remove_outputs()
    command = ['./athena', '-i', _athena_input_path()]
    logger.info('Executing: %s', ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed\n' + output)

    bfirst, blast = _endpoints('mhd_bcc')
    dfirst, dlast = _endpoints('mhd_divb')
    time0 = float(bfirst['Time'])
    time1 = float(blast['Time'])
    a0 = 1.0 + _RATES*time0
    a1 = 1.0 + _RATES*time1
    physical_ratio = np.asarray([
        (a0[1]*a0[2])/(a1[1]*a1[2]),
        (a0[0]*a0[2])/(a1[0]*a1[2]),
        (a0[0]*a0[1])/(a1[0]*a1[1]),
    ])
    _RESULTS.update({
        'time0': time0,
        'time1': time1,
        'initial_divb_max': float(np.max(np.abs(dfirst['divb']))),
        'final_divb_max': float(np.max(np.abs(dlast['divb']))),
        'bcc1_error': _relative_error(
            blast['bcc1'], physical_ratio[0]*bfirst['bcc1']),
        'bcc2_error': _relative_error(
            blast['bcc2'], physical_ratio[1]*bfirst['bcc2']),
        'bcc3_error': _relative_error(
            blast['bcc3'], physical_ratio[2]*bfirst['bcc3']),
    })


def analyze():
    logger.debug('Analyzing test ' + __name__)
    logger.info('oblique expanding-box metrics=%s', _RESULTS)
    return (
        _RESULTS['time1'] > _RESULTS['time0']
        and _RESULTS['initial_divb_max'] <= 1.0e-12
        and _RESULTS['final_divb_max'] <= 1.0e-12
        and _RESULTS['bcc1_error'] <= 5.0e-7
        and _RESULTS['bcc2_error'] <= 5.0e-7
        and _RESULTS['bcc3_error'] <= 5.0e-7
    )
