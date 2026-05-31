"""Bounded Q008 Appendix-A CPAW exponential-history preparation matrix."""

import glob
import json
import logging
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'vis/python'))
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_mhd_expanding_box_cpaw_history_preparation.athinput'
_BASENAME = 'pic_mhd_expanding_box_cpaw_history'
_RESOLUTIONS = (32, 64, 128)
_PROFILE_RATES = {
    'static_on': 0.0,
    'expanding': 0.01,
    'compressing': -0.01,
}
_OFF_CONTROL = 'static_off_control'
_WAVENUMBER = 2.0*np.pi
_ALFVEN_SPEED = 1.0
_MINIMUM_CROSSING_TIMES = 4.0
_PREPARATION_FIELD_BOUND = 0.25
_PREPARATION_PHASE_ORACLE_BOUND = 0.25
_PREPARATION_CONSTANT_HISTORY_BOUND = 1.0e-6
_PREPARATION_SCALED_ENERGY_BOUND = 0.25
_PREPARATION_ELSASSER_BOUND = 0.01
_PREPARATION_PHASE_LOCK_BOUND = 0.01
_PREPARATION_PARITY_BOUND = 1.0e-6
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        'ATHENA_Q008_CPAW_HISTORY_EXE_DIR',
        os.path.join(os.getcwd(), 'build', 'src'),
    )


def _athena_run_dir():
    return os.environ.get(
        'ATHENA_Q008_CPAW_HISTORY_RUN_DIR',
        _athena_exe_dir(),
    )


def _athena_input_path():
    return os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), '..', '..', '..', 'inputs', _INPUT_DECK
        )
    )


def _basename(profile, nx1):
    return f'{_BASENAME}_{profile}_nx{nx1}'


def _remove_outputs(basename):
    run_dir = _athena_run_dir()
    patterns = [
        os.path.join(run_dir, 'bin', basename + '.*.bin'),
        os.path.join(run_dir, basename + '.mhd.hst'),
        os.path.join(run_dir, basename + '-errs.dat'),
    ]
    for pattern in patterns:
        for path in glob.glob(pattern):
            os.remove(path)


def _run_case(profile, nx1, mode, rate):
    basename = _basename(profile, nx1)
    os.makedirs(_athena_run_dir(), exist_ok=True)
    _remove_outputs(basename)
    rates = (2.0*rate, rate, rate)
    arguments = [
        'job/basename=' + basename,
        'mesh/nx1=' + str(nx1),
        'meshblock/nx1=' + str(nx1),
        'particles/pic_expanding_box_mode=' + mode,
        'particles/pic_expansion_rate_x1=' + str(rates[0]),
        'particles/pic_expansion_rate_x2=' + str(rates[1]),
        'particles/pic_expansion_rate_x3=' + str(rates[2]),
    ]
    command = [
        os.path.join(_athena_exe_dir(), 'athena'),
        '-i',
        _athena_input_path(),
    ] + arguments
    logger.info('Executing %s nx1=%d: %s', profile, nx1, ' '.join(command))
    proc = subprocess.run(
        command, cwd=_athena_run_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed for ' + basename + '\n' + output)


def _snapshots(basename):
    pattern = os.path.join(
        _athena_run_dir(), 'bin', basename + '.mhd_w_bcc.*.bin'
    )
    paths = sorted(glob.glob(pattern))
    if len(paths) < 2:
        raise RuntimeError('Need at least two mhd_w_bcc snapshots for ' + pattern)
    return [bin_convert.read_binary_as_athdf(path) for path in paths]


def _history_rows(basename):
    path = Path(_athena_run_dir()) / (basename + '.mhd.hst')
    rows = np.loadtxt(path)
    if rows.ndim == 1:
        rows = rows[np.newaxis, :]
    if rows.shape[0] < 2:
        raise RuntimeError('Need at least two MHD history rows for ' + str(path))
    return rows


def _scale_factor(rate, time):
    return np.exp(rate*time)


def _phase_oracle(rate, time):
    if rate == 0.0:
        return _WAVENUMBER*_ALFVEN_SPEED*time
    return (
        _WAVENUMBER*_ALFVEN_SPEED*(1.0 - np.exp(-2.0*rate*time))
        / (2.0*rate)
    )


def _mode_coefficient(snapshot, y_name, z_name):
    phase = np.exp(1j*_WAVENUMBER*snapshot['x1v'])[None, None, :]
    return complex(np.mean((snapshot[y_name] + 1j*snapshot[z_name])*phase))


def _relative_drift(values, expected):
    values = np.asarray(values)
    expected = np.asarray(expected)
    return float(np.max(np.abs(values/values[0] - expected)))


def _wrapped_abs(values):
    values = np.asarray(values)
    return np.abs(np.arctan2(np.sin(values), np.cos(values)))


def _field_oracle_metrics(snapshots, rate):
    times = np.asarray([float(snapshot['Time']) for snapshot in snapshots])
    scales = _scale_factor(rate, times)
    rho0 = float(np.mean(snapshots[0]['dens']))
    bx0 = float(np.mean(snapshots[0]['bcc1']))
    u_modes = np.asarray(
        [_mode_coefficient(snapshot, 'vely', 'velz') for snapshot in snapshots]
    )
    b_modes = np.asarray(
        [_mode_coefficient(snapshot, 'bcc2', 'bcc3') for snapshot in snapshots]
    )
    rho_error = max(
        float(np.max(np.abs(snapshot['dens']/rho0 - scale**-4)))
        for snapshot, scale in zip(snapshots, scales)
    )
    bx_error = max(
        float(np.max(np.abs(snapshot['bcc1']/bx0 - scale**-2)))
        for snapshot, scale in zip(snapshots, scales)
    )
    velocity_amplitude_error = np.abs(np.abs(u_modes/u_modes[0]) - scales**-1)
    magnetic_amplitude_error = np.abs(np.abs(b_modes/b_modes[0]) - scales**-3)
    phase_expected = np.asarray([_phase_oracle(rate, time) for time in times])
    velocity_phase = np.unwrap(np.angle(u_modes/u_modes[0]))
    magnetic_phase = np.unwrap(np.angle(b_modes/b_modes[0]))
    velocity_phase_error = np.abs(np.abs(velocity_phase) - phase_expected)
    magnetic_phase_error = np.abs(np.abs(magnetic_phase) - phase_expected)

    rho_means = np.asarray(
        [float(np.mean(snapshot['dens'])) for snapshot in snapshots]
    )
    magnetic_alfven_modes = b_modes/np.sqrt(rho_means)
    zplus = u_modes + magnetic_alfven_modes
    zminus = u_modes - magnetic_alfven_modes
    if abs(zplus[0]) >= abs(zminus[0]):
        desired = zplus
        undesired = zminus
        desired_branch = 'zplus'
    else:
        desired = zminus
        undesired = zplus
        desired_branch = 'zminus'
    phase_lock = _wrapped_abs(np.angle((u_modes/b_modes)/(u_modes[0]/b_modes[0])))
    return {
        'rho_relative_error_max': rho_error,
        'parallel_magnetic_field_relative_error_max': bx_error,
        'velocity_mode_amplitude_relative_error_max': float(
            np.max(velocity_amplitude_error)
        ),
        'velocity_mode_amplitude_endpoint_relative_error': float(
            velocity_amplitude_error[-1]
        ),
        'magnetic_mode_amplitude_relative_error_max': float(
            np.max(magnetic_amplitude_error)
        ),
        'magnetic_mode_amplitude_endpoint_relative_error': float(
            magnetic_amplitude_error[-1]
        ),
        'velocity_phase_absolute_error_max': float(np.max(velocity_phase_error)),
        'velocity_phase_endpoint_absolute_error': float(velocity_phase_error[-1]),
        'magnetic_phase_absolute_error_max': float(np.max(magnetic_phase_error)),
        'magnetic_phase_endpoint_absolute_error': float(magnetic_phase_error[-1]),
        'desired_elsasser_branch': desired_branch,
        'undesired_to_desired_elsasser_ratio_max': float(
            np.max(np.abs(undesired)/np.abs(desired))
        ),
        'velocity_magnetic_phase_lock_absolute_error_max': float(
            np.max(phase_lock)
        ),
    }


def _history_oracle_metrics(history, rate):
    times = history[:, 0]
    scales = _scale_factor(rate, times)
    transverse_kinetic_energy = history[:, 8] + history[:, 9]
    transverse_magnetic_energy = history[:, 11] + history[:, 12]
    return {
        'mass_relative_drift_max': _relative_drift(
            history[:, 2], np.ones_like(times)
        ),
        'parallel_magnetic_energy_relative_drift_max': _relative_drift(
            history[:, 10], np.ones_like(times)
        ),
        'transverse_kinetic_energy_scaled_relative_error_max': _relative_drift(
            transverse_kinetic_energy, scales**-2
        ),
        'transverse_magnetic_energy_scaled_relative_error_max': _relative_drift(
            transverse_magnetic_energy, scales**-2
        ),
    }


def _measure_case(profile, nx1, mode, rate):
    basename = _basename(profile, nx1)
    snapshots = _snapshots(basename)
    history = _history_rows(basename)
    return {
        'profile': profile,
        'nx1': nx1,
        'mode': mode,
        'rates': [2.0*rate, rate, rate],
        'time': float(snapshots[-1]['Time']),
        'snapshot_count': len(snapshots),
        'history_row_count': int(history.shape[0]),
        'field_oracles': _field_oracle_metrics(snapshots, rate),
        'history_oracles': _history_oracle_metrics(history, rate),
    }


def _max_relative_difference(actual, expected):
    scale = max(float(np.max(np.abs(expected))), 1.0e-30)
    return float(np.max(np.abs(actual - expected)))/scale


def _control_parity():
    mode_on_basename = _basename('static_on', 128)
    mode_off_basename = _basename(_OFF_CONTROL, 128)
    mode_on_snapshots = _snapshots(mode_on_basename)
    mode_off_snapshots = _snapshots(mode_off_basename)
    if len(mode_on_snapshots) != len(mode_off_snapshots):
        raise RuntimeError('Static control snapshot count mismatch')
    field_error = 0.0
    for mode_on, mode_off in zip(mode_on_snapshots, mode_off_snapshots):
        for name in ['dens', 'velx', 'vely', 'velz', 'eint',
                     'bcc1', 'bcc2', 'bcc3']:
            field_error = max(
                field_error,
                _max_relative_difference(mode_on[name], mode_off[name]),
            )
    history_error = _max_relative_difference(
        _history_rows(mode_on_basename), _history_rows(mode_off_basename)
    )
    return {
        'mhd_w_bcc_relative_difference_max': field_error,
        'mhd_history_relative_difference_max': history_error,
    }


def _apparent_order(coarse, fine):
    if coarse <= 0.0 or fine <= 0.0:
        return None
    return float(np.log(coarse/fine)/np.log(2.0))


def _convergence_observations(cases):
    observations = {}
    metrics = [
        'velocity_mode_amplitude_endpoint_relative_error',
        'magnetic_mode_amplitude_endpoint_relative_error',
        'velocity_phase_endpoint_absolute_error',
        'magnetic_phase_endpoint_absolute_error',
    ]
    for profile in _PROFILE_RATES:
        profile_cases = [cases[f'{profile}_nx{nx1}'] for nx1 in _RESOLUTIONS]
        observations[profile] = {}
        for metric in metrics:
            residuals = [
                case['field_oracles'][metric] for case in profile_cases
            ]
            observations[profile][metric] = {
                'residuals_nx1_32_64_128': residuals,
                'apparent_orders_32_to_64_and_64_to_128': [
                    _apparent_order(residuals[0], residuals[1]),
                    _apparent_order(residuals[1], residuals[2]),
                ],
            }
    return observations


def _finite(value):
    if isinstance(value, dict):
        return all(_finite(item) for item in value.values())
    if isinstance(value, list):
        return all(_finite(item) for item in value)
    if value is None or isinstance(value, str):
        return True
    return bool(np.isfinite(value))


def _case_passes(case):
    field = case['field_oracles']
    history = case['history_oracles']
    field_metrics = [
        field['rho_relative_error_max'],
        field['parallel_magnetic_field_relative_error_max'],
        field['velocity_mode_amplitude_relative_error_max'],
        field['magnetic_mode_amplitude_relative_error_max'],
    ]
    phase_metrics = [
        field['velocity_phase_absolute_error_max'],
        field['magnetic_phase_absolute_error_max'],
    ]
    constant_history_metrics = [
        history['mass_relative_drift_max'],
        history['parallel_magnetic_energy_relative_drift_max'],
    ]
    scaled_energy_metrics = [
        history['transverse_kinetic_energy_scaled_relative_error_max'],
        history['transverse_magnetic_energy_scaled_relative_error_max'],
    ]
    return (
        _finite(case)
        and case['time'] >= _MINIMUM_CROSSING_TIMES - 1.0e-12
        and case['snapshot_count'] >= 2
        and case['history_row_count'] >= 2
        and max(field_metrics) < _PREPARATION_FIELD_BOUND
        and max(phase_metrics) < _PREPARATION_PHASE_ORACLE_BOUND
        and max(constant_history_metrics) < _PREPARATION_CONSTANT_HISTORY_BOUND
        and max(scaled_energy_metrics) < _PREPARATION_SCALED_ENERGY_BOUND
        and field[
            'undesired_to_desired_elsasser_ratio_max'
        ] < _PREPARATION_ELSASSER_BOUND
        and field[
            'velocity_magnetic_phase_lock_absolute_error_max'
        ] < _PREPARATION_PHASE_LOCK_BOUND
    )


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _RESULTS.clear()
    cases = {}
    for profile, rate in _PROFILE_RATES.items():
        for nx1 in _RESOLUTIONS:
            _run_case(profile, nx1, 'on', rate)
            key = f'{profile}_nx{nx1}'
            cases[key] = _measure_case(profile, nx1, 'on', rate)
    _run_case(_OFF_CONTROL, 128, 'off', 0.0)
    cases[_OFF_CONTROL + '_nx128'] = _measure_case(
        _OFF_CONTROL, 128, 'off', 0.0
    )
    _RESULTS['cases'] = cases
    _RESULTS['static_mode_on_off_parity'] = _control_parity()
    _RESULTS['convergence_observations'] = _convergence_observations(cases)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    for name, case in _RESULTS['cases'].items():
        logger.info('%s metrics=%s', name, case)
    logger.info(
        'static mode-on/off parity=%s', _RESULTS['static_mode_on_off_parity']
    )
    logger.info(
        'convergence observations=%s', _RESULTS['convergence_observations']
    )
    parity = _RESULTS['static_mode_on_off_parity']
    return (
        all(_case_passes(case) for case in _RESULTS['cases'].values())
        and parity['mhd_w_bcc_relative_difference_max'] <= _PREPARATION_PARITY_BOUND
        and parity['mhd_history_relative_difference_max'] <= _PREPARATION_PARITY_BOUND
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    run()
    if not analyze():
        raise SystemExit('pic_mhd_expanding_box_cpaw_history_preparation: FAIL')
    print(json.dumps(_RESULTS, indent=2, sort_keys=True))
