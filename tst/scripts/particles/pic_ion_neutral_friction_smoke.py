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

_INPUT_DECK = 'tests/pic_ion_neutral_friction_smoke.athinput'
_NU_IN = 0.7
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _remove_outputs(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.*.bin')
    for fname in glob.glob(pattern):
        os.remove(fname)


def _output_files(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.mhd_u.*.bin')
    matches = sorted(glob.glob(pattern))
    if len(matches) < 2:
        raise RuntimeError('Need at least two output snapshots for: ' + pattern)
    return matches


def _run_case(case_name, damping_mode, collision_rate):
    basename = 'pic_ion_neutral_friction_' + case_name
    args = [
        'job/basename=' + basename,
        'particles/pic_wave_damping_mode=' + damping_mode,
        'particles/pic_ion_neutral_collision_rate=' + str(collision_rate),
    ]
    command = ['./athena', '-i', _athena_input_path()] + args
    logger.info('Executing %s: %s', case_name, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + case_name + '\n' + output)
    expected = [
        'physical_mode=extended_mhd_pic',
        'wave_damping=' + damping_mode,
        'nu_in=',
    ]
    for token in expected:
        if token not in output:
            raise RuntimeError(case_name + ' missing runtime token: ' + token)
    files = _output_files(basename)
    return (bin_convert.read_binary_as_athdf(files[0]),
            bin_convert.read_binary_as_athdf(files[-1]))


def _run_guard(case_name, args, expected_message):
    command = ['./athena', '-i', _athena_input_path(), 'time/nlim=0'] + args
    logger.info('Executing guard %s: %s', case_name, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode == 0:
        raise RuntimeError('Guard unexpectedly passed: ' + case_name)
    if expected_message not in output:
        raise RuntimeError(case_name + ' missing guard reason\n' + output)


def _relative_error(actual, expected):
    scale = max(float(np.max(np.abs(expected))), 1.0e-14)
    return float(np.max(np.abs(actual - expected))) / scale


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    for case_name in ('off', 'damped'):
        _remove_outputs('pic_ion_neutral_friction_' + case_name)
    _RESULTS['off'] = _run_case('off', 'off', 0.0)
    _RESULTS['damped'] = _run_case('damped', 'ion_neutral_friction', _NU_IN)
    _run_guard(
        'paper_mode_rejects_damping',
        ['particles/pic_physical_mode=paper_mhd_pic',
         'particles/pic_wave_damping_mode=ion_neutral_friction',
         'particles/pic_ion_neutral_collision_rate=0.7'],
        'requires <particles>/pic_physical_mode=extended_mhd_pic')
    _run_guard(
        'friction_requires_positive_rate',
        ['particles/pic_wave_damping_mode=ion_neutral_friction',
         'particles/pic_ion_neutral_collision_rate=0.0'],
        'requires <particles>/pic_ion_neutral_collision_rate > 0')


def analyze():
    logger.debug('Analyzing test ' + __name__)
    _off_first, off = _RESULTS['off']
    damp_first, damp = _RESULTS['damped']
    dt = float(damp['Time'] - damp_first['Time'])
    factor = float(np.exp(-_NU_IN * dt))
    expected_energy = (
        off['ener']
        - 0.5 * (1.0 - factor * factor)
        * (off['mom2'] * off['mom2'] + off['mom3'] * off['mom3'])
        / off['dens'])
    errors = {
        'dens': _relative_error(damp['dens'], off['dens']),
        'mom1': _relative_error(damp['mom1'], off['mom1']),
        'mom2': _relative_error(damp['mom2'], factor * off['mom2']),
        'mom3': _relative_error(damp['mom3'], factor * off['mom3']),
        'ener': _relative_error(damp['ener'], expected_energy),
    }
    transverse_norm = float(np.linalg.norm(off['mom2'])
                            + np.linalg.norm(off['mom3']))
    measured = {
        'dt': dt,
        'factor': factor,
        'transverse_norm': transverse_norm,
        'errors': errors,
    }
    logger.info('ion-neutral friction manufactured-source metrics: %s', measured)
    _RESULTS['metrics'] = measured
    return (dt > 0.0 and transverse_norm > 1.0e-8
            and max(errors.values()) <= 5.0e-7)
