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

_INPUT_DECK = 'tests/pic_extended_hall_ct_smoke.athinput'
_HALL_MODE = 'current_to_ct_experimental'
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _remove_outputs(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.*.bin')
    for fname in glob.glob(pattern):
        os.remove(fname)


def _latest_output_file(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           basename + '.mhd_bcc.*.bin')
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise RuntimeError('No output files found for pattern: ' + pattern)
    return matches[-1]


def _read_bcc(basename):
    data = bin_convert.read_binary_as_athdf(_latest_output_file(basename))
    return np.concatenate([
        np.asarray(data['bcc1'], dtype=float).ravel(),
        np.asarray(data['bcc2'], dtype=float).ravel(),
        np.asarray(data['bcc3'], dtype=float).ravel(),
    ])


def _run_case(case_name, hall_mode, coefficient):
    basename = 'pic_extended_hall_ct_' + case_name
    args = [
        'job/basename=' + basename,
        'particles/pic_cr_hall_mode=' + hall_mode,
        'particles/couple_j_to_efield_coeff=' + str(coefficient),
    ]
    command = ['./athena', '-i', _athena_input_path()] + args
    logger.info('Executing %s: %s', case_name, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + case_name + '\n' + output)

    induction = ('ideal_mhd_only' if hall_mode == 'off'
                 else 'cr_current_to_ct')
    expected = [
        'physical_mode=extended_mhd_pic',
        'state=momentum_p_over_m',
        'feedback=coupled',
        'induction=' + induction,
    ]
    for token in expected:
        if token not in output:
            raise RuntimeError(case_name + ' missing runtime token: ' + token)
    return basename


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    cases = (
        ('off', 'off', 1.0),
        ('plus', _HALL_MODE, 1.0),
        ('minus', _HALL_MODE, -1.0),
    )
    for case_name, hall_mode, coefficient in cases:
        basename = 'pic_extended_hall_ct_' + case_name
        _remove_outputs(basename)
        _run_case(case_name, hall_mode, coefficient)
        _RESULTS[case_name] = _read_bcc(basename)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    plus_delta = _RESULTS['plus'] - _RESULTS['off']
    minus_delta = _RESULTS['minus'] - _RESULTS['off']
    plus_norm = float(np.linalg.norm(plus_delta))
    minus_norm = float(np.linalg.norm(minus_delta))
    scale = max(plus_norm, minus_norm, 1.0e-30)
    odd_error = float(np.linalg.norm(plus_delta + minus_delta)) / scale
    cosine = float(np.dot(plus_delta, minus_delta)
                   / max(plus_norm * minus_norm, 1.0e-30))
    measured = {
        'plus_norm': plus_norm,
        'minus_norm': minus_norm,
        'odd_error': odd_error,
        'cosine': cosine,
    }
    logger.info('extended Hall CT manufactured-source metrics: %s', measured)
    _RESULTS['metrics'] = measured
    return (plus_norm > 1.0e-8
            and minus_norm > 1.0e-8
            and odd_error <= 5.0e-5
            and cosine <= -0.9999)
