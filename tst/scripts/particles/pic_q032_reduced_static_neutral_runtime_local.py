import glob
import logging
import os
import shlex
import subprocess
import sys
import tempfile

import numpy as np

_SOURCE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, os.path.join(_SOURCE_ROOT, 'vis', 'python'))
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_q032_reduced_static_neutral_runtime_local.athinput'
_NU_IN = 0.7
_RESULTS = {}
_DAMPING_ADMISSION = (
    '<particles>/pic_wave_damping_mode=ion_neutral_friction requires '
    'the Newtonian single-fluid MHD source task path; does not support '
)
_MHD_FIELDS = ('dens', 'mom1', 'mom2', 'mom3', 'ener')


def _athena_exe_dir():
    return os.environ.get(
        'ATHENA_Q032_REDUCED_STATIC_NEUTRAL_EXE_DIR',
        os.path.join(os.getcwd(), 'build', 'src'))


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, 'inputs', _INPUT_DECK)


def _rank_count():
    try:
        nproc = int(os.environ.get('ATHENA_Q032_NPROC', '1'))
    except ValueError as error:
        raise RuntimeError('ATHENA_Q032_NPROC must be 1 or 2') from error
    if nproc not in (1, 2):
        raise RuntimeError('ATHENA_Q032_NPROC must be 1 or 2')
    return nproc


def _launcher_prefix():
    if _rank_count() == 1:
        return []
    launcher = shlex.split(os.environ.get('ATHENA_Q032_LAUNCHER', 'mpiexec'))
    if not launcher:
        raise RuntimeError('ATHENA_Q032_LAUNCHER must not be empty')
    return launcher + ['-n', str(_rank_count())]


def _decomposition_overrides():
    if _rank_count() == 1:
        return []
    return ['meshblock/nx1=16']


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


def _run_case(case_name, arguments, expected_damping_mode):
    basename = 'pic_q032_reduced_static_neutral_runtime_local_' + case_name
    _remove_outputs(basename)
    command = _launcher_prefix() + [
        './athena', '-i', _athena_input_path(), 'job/basename=' + basename,
    ] + arguments + _decomposition_overrides()
    logger.info('Executing %s: %s', case_name, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + case_name + '\n' + output)
    for token in [
            'wave_damping=' + expected_damping_mode,
            'nu_in=']:
        if token not in output:
            raise RuntimeError(case_name + ' missing runtime token: ' + token)
    files = _output_files(basename)
    return (bin_convert.read_binary_as_athdf(files[0]),
            bin_convert.read_binary_as_athdf(files[-1]))


def _guard_input_path(case_name, appended_input):
    if not appended_input:
        return _athena_input_path(), None
    with open(_athena_input_path(), encoding='utf-8') as source:
        deck = source.read()
    with tempfile.NamedTemporaryFile(
            mode='w', encoding='utf-8', prefix='athena_q032_' + case_name + '_',
            suffix='.athinput', delete=False) as temporary:
        temporary.write(deck)
        temporary.write('\n')
        temporary.write(appended_input.strip())
        temporary.write('\n')
        return temporary.name, temporary.name


def _run_guard(case_name, arguments, expected_message, appended_input=''):
    input_path, temporary_input = _guard_input_path(case_name, appended_input)
    command = [
        './athena', '-i', input_path, 'time/nlim=0',
    ] + arguments
    logger.info('Executing guard %s: %s', case_name, ' '.join(command))
    try:
        proc = subprocess.run(command, cwd=_athena_exe_dir(),
                              capture_output=True, text=True)
    finally:
        if temporary_input is not None:
            os.remove(temporary_input)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode == 0:
        raise RuntimeError('Guard unexpectedly passed: ' + case_name)
    if expected_message not in output:
        raise RuntimeError(case_name + ' missing guard reason\n' + output)


def _relative_error(actual, expected):
    scale = max(float(np.max(np.abs(expected))), 1.0e-14)
    return float(np.max(np.abs(actual - expected))) / scale


def _require_equal_snapshots(label, actual, expected):
    if float(actual['Time']) != float(expected['Time']):
        raise RuntimeError(label + ' snapshot times differ')
    for field in _MHD_FIELDS:
        if not np.array_equal(actual[field], expected[field]):
            raise RuntimeError(label + ' snapshot arrays differ for ' + field)


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _RESULTS.clear()
    _RESULTS['control'] = _run_case(
        'control',
        ['problem/pgen_name=linear_wave',
         'particles/pic_wave_damping_mode=off',
         'particles/pic_ion_neutral_collision_rate=0.0'],
        'off')
    _RESULTS['damped'] = _run_case('damped', [], 'ion_neutral_friction')
    _run_guard(
        'rejects_amr',
        ['mesh_refinement/refinement=static'],
        'q032_reduced_static_neutral_local rejects AMR/SMR')
    _run_guard(
        'rejects_nonperiodic_carrier',
        ['mesh/ix1_bc=outflow', 'mesh/ox1_bc=outflow'],
        'requires periodic active carrier boundaries')
    _run_guard(
        'rejects_role_drift',
        ['q032_reduced_static_neutral_local/deck_role=plotnikov_evidence'],
        '<q032_reduced_static_neutral_local>/deck_role does not match')
    _run_guard(
        'rejects_wrapper_special_relativistic',
        ['particles/pic_wave_damping_mode=off',
         'particles/pic_ion_neutral_collision_rate=0.0'],
        'q032_reduced_static_neutral_local requires the exact Newtonian '
        'single-fluid MHD source task path',
        '''
<coord>
special_rel = true
general_rel = false
''')
    _run_guard(
        'rejects_special_relativistic',
        [],
        _DAMPING_ADMISSION + '<coord>/special_rel=true',
        '''
<coord>
special_rel = true
general_rel = false
''')
    _run_guard(
        'rejects_general_relativistic',
        [],
        _DAMPING_ADMISSION + '<coord>/general_rel=true',
        '''
<coord>
special_rel = false
general_rel = true
minkowski   = true
''')
    _run_guard(
        'rejects_radiation_task_path',
        [],
        _DAMPING_ADMISSION + '<radiation> task paths',
        '''
<coord>
special_rel = false
general_rel = true
minkowski   = true

<radiation>
rad_source  = false
fixed_fluid = true
nlevel      = 0
rotate_geo  = false
angular_fluxes = false
''')
    _run_guard(
        'rejects_ion_neutral_task_path',
        [],
        _DAMPING_ADMISSION + '<ion-neutral> alternate task lists',
        '''
<hydro>
eos         = ideal
reconstruct = plm
rsolver     = llf
gamma       = 1.66666666667

<ion-neutral>
drag_coeff  = 1.0
''')
    _run_guard(
        'rejects_dynamical_gr_task_path',
        [],
        _DAMPING_ADMISSION + 'dynamical-GR <adm>/<z4c> task paths',
        '''
<coord>
minkowski = true

<mhd>
dyn_eos   = ideal
dyn_error = reset_floor

<adm>
''')


def analyze():
    logger.debug('Analyzing test ' + __name__)
    _control_first, control = _RESULTS['control']
    damped_first, damped = _RESULTS['damped']
    _require_equal_snapshots('initial control/damped', _control_first, damped_first)
    if float(control['Time']) != float(damped['Time']):
        raise RuntimeError('final control/damped snapshot times differ')
    dt = float(damped['Time'] - damped_first['Time'])
    factor = float(np.exp(-_NU_IN * dt))
    expected_energy = (
        control['ener']
        - 0.5 * (1.0 - factor * factor)
        * (control['mom2'] * control['mom2'] + control['mom3'] * control['mom3'])
        / control['dens'])
    errors = {
        'dens': _relative_error(damped['dens'], control['dens']),
        'mom1': _relative_error(damped['mom1'], control['mom1']),
        'mom2': _relative_error(damped['mom2'], factor * control['mom2']),
        'mom3': _relative_error(damped['mom3'], factor * control['mom3']),
        'ener': _relative_error(damped['ener'], expected_energy),
    }
    transverse_norm = float(np.linalg.norm(control['mom2'])
                            + np.linalg.norm(control['mom3']))
    measured = {
        'dt': dt,
        'factor': factor,
        'transverse_norm': transverse_norm,
        'configured_ranks': _rank_count(),
        'errors': errors,
        'qualification_effect': 'none',
        'plotnikov_qualification': 'not_claimed',
    }
    logger.info('Q032 bounded local mechanics metrics: %s', measured)
    _RESULTS['metrics'] = measured
    return (dt > 0.0 and transverse_norm > 1.0e-8
            and max(errors.values()) <= 5.0e-7)
