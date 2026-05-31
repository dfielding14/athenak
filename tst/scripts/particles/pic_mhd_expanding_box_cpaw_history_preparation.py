"""Bounded Q008 Appendix-A CPAW expanding-box history preparation matrix."""

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
_LITERAL_PROFILES = {
    'linear_expanding': ('linear', 0.01),
    'linear_compressing': ('linear', -0.01),
    'reciprocal_linear_expanding': ('reciprocal_linear', 0.01),
    'reciprocal_linear_compressing': ('reciprocal_linear', -0.01),
}
_DIMENSIONAL_PROFILES = {
    'exponential_expanding': ('exponential', 0.01),
    'exponential_compressing': ('exponential', -0.01),
}
_AXES = {
    'x1': {
        'index': 0,
        'coordinate': 'x1v',
        'phase_shape': (1, 1, -1),
        'parallel_b': 'bcc1',
        'velocity_transverse': ('vely', 'velz'),
        'magnetic_transverse': ('bcc2', 'bcc3'),
        'history_parallel_me': 10,
        'history_transverse_ke': (8, 9),
        'history_transverse_me': (11, 12),
    },
    'x2': {
        'index': 1,
        'coordinate': 'x2v',
        'phase_shape': (1, -1, 1),
        'parallel_b': 'bcc2',
        'velocity_transverse': ('velz', 'velx'),
        'magnetic_transverse': ('bcc3', 'bcc1'),
        'history_parallel_me': 11,
        'history_transverse_ke': (9, 7),
        'history_transverse_me': (12, 10),
    },
    'x3': {
        'index': 2,
        'coordinate': 'x3v',
        'phase_shape': (-1, 1, 1),
        'parallel_b': 'bcc3',
        'velocity_transverse': ('velx', 'vely'),
        'magnetic_transverse': ('bcc1', 'bcc2'),
        'history_parallel_me': 12,
        'history_transverse_ke': (7, 8),
        'history_transverse_me': (10, 11),
    },
}
_SOLVERS = ('llf', 'hlle', 'hlld')
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
_ROE_REJECTION_TEXT = "<mhd>/rsolver = 'roe' not implemented for dynamic problems"
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


def _spec(profile, law, rate, axis, resolution, solver='llf', mode='on'):
    return {
        'profile': profile,
        'law': law,
        'rate': rate,
        'axis': axis,
        'resolution': resolution,
        'solver': solver,
        'mode': mode,
    }


def _successful_specs():
    specs = []
    for profile, (law, rate) in _LITERAL_PROFILES.items():
        for resolution in _RESOLUTIONS:
            specs.append(_spec(profile, law, rate, 'x1', resolution))
    specs.extend([
        _spec('static_on', 'exponential', 0.0, 'x1', 128),
        _spec('static_off_control', 'exponential', 0.0, 'x1', 128, mode='off'),
    ])
    for profile, (law, rate) in _DIMENSIONAL_PROFILES.items():
        for axis in _AXES:
            specs.append(_spec(profile, law, rate, axis, 64))
    for solver in _SOLVERS[1:]:
        specs.append(
            _spec('exponential_expanding', 'exponential', 0.01, 'x1', 64, solver)
        )
    return specs


def _case_id(spec):
    return '{profile}_{axis}_{solver}_nx{resolution}'.format(**spec)


def _basename(spec):
    return _BASENAME + '_' + _case_id(spec)


def _rates(spec):
    rates = [spec['rate'], spec['rate'], spec['rate']]
    rates[_AXES[spec['axis']]['index']] = 2.0*spec['rate']
    return rates


def _dimensions(spec):
    dimensions = [4, 4, 4]
    dimensions[_AXES[spec['axis']]['index']] = spec['resolution']
    return dimensions


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


def _outputs(basename):
    run_dir = _athena_run_dir()
    patterns = [
        os.path.join(run_dir, 'bin', basename + '.*.bin'),
        os.path.join(run_dir, basename + '.mhd.hst'),
        os.path.join(run_dir, basename + '-errs.dat'),
    ]
    return sorted(
        os.path.relpath(path, run_dir)
        for pattern in patterns
        for path in glob.glob(pattern)
    )


def _command(spec):
    dimensions = _dimensions(spec)
    rates = _rates(spec)
    arguments = [
        'job/basename=' + _basename(spec),
        'mesh/nx1=' + str(dimensions[0]),
        'mesh/nx2=' + str(dimensions[1]),
        'mesh/nx3=' + str(dimensions[2]),
        'meshblock/nx1=' + str(dimensions[0]),
        'meshblock/nx2=' + str(dimensions[1]),
        'meshblock/nx3=' + str(dimensions[2]),
        'mhd/rsolver=' + spec['solver'],
        'particles/pic_expanding_box_mode=' + spec['mode'],
        'particles/pic_expansion_law=' + spec['law'],
        'particles/pic_expansion_rate_x1=' + str(rates[0]),
        'particles/pic_expansion_rate_x2=' + str(rates[1]),
        'particles/pic_expansion_rate_x3=' + str(rates[2]),
        'problem/along_x1=' + str(spec['axis'] == 'x1').lower(),
        'problem/along_x2=' + str(spec['axis'] == 'x2').lower(),
        'problem/along_x3=' + str(spec['axis'] == 'x3').lower(),
    ]
    return [
        os.path.join(_athena_exe_dir(), 'athena'),
        '-i',
        _athena_input_path(),
    ] + arguments


def _run_case(spec):
    basename = _basename(spec)
    os.makedirs(_athena_run_dir(), exist_ok=True)
    _remove_outputs(basename)
    command = _command(spec)
    logger.info('Executing %s: %s', _case_id(spec), ' '.join(command))
    proc = subprocess.run(
        command, cwd=_athena_run_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed for ' + basename + '\n' + output)


def _run_roe_rejection_probe():
    spec = _spec('roe_unsupported_probe', 'exponential', 0.0, 'x1', 32, 'roe')
    basename = _basename(spec)
    _remove_outputs(basename)
    command = _command(spec) + ['time/nlim=0']
    logger.info('Executing Roe rejection probe: %s', ' '.join(command))
    proc = subprocess.run(
        command, cwd=_athena_run_dir(), capture_output=True, text=True
    )
    probe_dir = Path(_athena_run_dir()) / 'roe_unsupported_probe'
    probe_dir.mkdir(parents=True, exist_ok=True)
    (probe_dir / 'command.json').write_text(
        json.dumps(command, indent=2) + '\n', encoding='utf-8'
    )
    (probe_dir / 'stdout.txt').write_text(proc.stdout or '', encoding='utf-8')
    (probe_dir / 'stderr.txt').write_text(proc.stderr or '', encoding='utf-8')
    (probe_dir / 'returncode.txt').write_text(
        str(proc.returncode) + '\n', encoding='utf-8'
    )
    output = (proc.stdout or '') + (proc.stderr or '')
    outputs = _outputs(basename)
    return {
        'profile': spec['profile'],
        'solver': spec['solver'],
        'returncode': proc.returncode,
        'expected_rejection_text': _ROE_REJECTION_TEXT,
        'expected_rejection_text_observed': _ROE_REJECTION_TEXT in output,
        'outputs_after_rejection': outputs,
        'raw_sidecars': {
            'command': 'roe_unsupported_probe/command.json',
            'stdout': 'roe_unsupported_probe/stdout.txt',
            'stderr': 'roe_unsupported_probe/stderr.txt',
            'returncode': 'roe_unsupported_probe/returncode.txt',
        },
        'status': (
            'pass_fail_closed_unsupported'
            if proc.returncode != 0 and _ROE_REJECTION_TEXT in output and not outputs
            else 'fail_unexpected_roe_disposition'
        ),
    }


def _snapshots(spec):
    pattern = os.path.join(
        _athena_run_dir(), 'bin', _basename(spec) + '.mhd_w_bcc.*.bin'
    )
    paths = sorted(glob.glob(pattern))
    if len(paths) < 2:
        raise RuntimeError('Need at least two mhd_w_bcc snapshots for ' + pattern)
    return [bin_convert.read_binary_as_athdf(path) for path in paths]


def _history_rows(spec):
    path = Path(_athena_run_dir()) / (_basename(spec) + '.mhd.hst')
    rows = np.loadtxt(path)
    if rows.ndim == 1:
        rows = rows[np.newaxis, :]
    if rows.shape[0] < 2:
        raise RuntimeError('Need at least two MHD history rows for ' + str(path))
    return rows


def _scale_factor(law, rate, time):
    if law == 'linear':
        return 1.0 + rate*time
    if law == 'reciprocal_linear':
        return 1.0/(1.0 + rate*time)
    if law == 'exponential':
        return np.exp(rate*time)
    raise ValueError('Unsupported expansion law ' + law)


def _geometry(spec, times):
    transverse = _scale_factor(spec['law'], spec['rate'], times)
    parallel = _scale_factor(spec['law'], 2.0*spec['rate'], times)
    return parallel, transverse


def _phase_oracle(law, rate, time):
    if rate == 0.0:
        integral = time
    elif law == 'linear':
        root = np.sqrt(1.0 + 2.0*rate*time)
        integral = 2.0*(np.arctan(root) - np.pi/4.0)/rate
    elif law == 'reciprocal_linear':
        root = np.sqrt(1.0 + 2.0*rate*time)
        integral = ((root**5 - 1.0)/5.0 + (root**3 - 1.0)/3.0)/(2.0*rate)
    elif law == 'exponential':
        integral = (1.0 - np.exp(-2.0*rate*time))/(2.0*rate)
    else:
        raise ValueError('Unsupported expansion law ' + law)
    return _WAVENUMBER*_ALFVEN_SPEED*integral


def _mode_coefficient(snapshot, axis, names):
    metadata = _AXES[axis]
    coordinate = np.asarray(snapshot[metadata['coordinate']])
    phase = np.exp(1j*_WAVENUMBER*coordinate).reshape(metadata['phase_shape'])
    return complex(np.mean((snapshot[names[0]] + 1j*snapshot[names[1]])*phase))


def _relative_drift(values, expected):
    values = np.asarray(values)
    expected = np.asarray(expected)
    return float(np.max(np.abs(values/values[0] - expected)))


def _wrapped_abs(values):
    values = np.asarray(values)
    return np.abs(np.arctan2(np.sin(values), np.cos(values)))


def _field_oracle_metrics(snapshots, spec):
    metadata = _AXES[spec['axis']]
    times = np.asarray([float(snapshot['Time']) for snapshot in snapshots])
    parallel, transverse = _geometry(spec, times)
    rho0 = float(np.mean(snapshots[0]['dens']))
    bparallel0 = float(np.mean(snapshots[0][metadata['parallel_b']]))
    u_modes = np.asarray([
        _mode_coefficient(snapshot, spec['axis'], metadata['velocity_transverse'])
        for snapshot in snapshots
    ])
    b_modes = np.asarray([
        _mode_coefficient(snapshot, spec['axis'], metadata['magnetic_transverse'])
        for snapshot in snapshots
    ])
    rho_error = max(
        float(np.max(np.abs(snapshot['dens']/rho0 - 1.0/(apar*atrans**2))))
        for snapshot, apar, atrans in zip(snapshots, parallel, transverse)
    )
    bparallel_error = max(
        float(np.max(np.abs(
            snapshot[metadata['parallel_b']]/bparallel0 - atrans**-2
        )))
        for snapshot, atrans in zip(snapshots, transverse)
    )
    velocity_amplitude_error = np.abs(np.abs(u_modes/u_modes[0]) - transverse**-1)
    magnetic_amplitude_error = np.abs(
        np.abs(b_modes/b_modes[0]) - 1.0/(parallel*transverse)
    )
    phase_expected = np.asarray([
        _phase_oracle(spec['law'], spec['rate'], time) for time in times
    ])
    velocity_phase = np.unwrap(np.angle(u_modes/u_modes[0]))
    magnetic_phase = np.unwrap(np.angle(b_modes/b_modes[0]))
    velocity_phase_error = np.abs(np.abs(velocity_phase) - phase_expected)
    magnetic_phase_error = np.abs(np.abs(magnetic_phase) - phase_expected)

    rho_means = np.asarray([
        float(np.mean(snapshot['dens'])) for snapshot in snapshots
    ])
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
        'parallel_magnetic_field_relative_error_max': bparallel_error,
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


def _history_oracle_metrics(history, spec):
    metadata = _AXES[spec['axis']]
    times = history[:, 0]
    parallel, transverse = _geometry(spec, times)
    transverse_kinetic_energy = sum(
        history[:, index] for index in metadata['history_transverse_ke']
    )
    transverse_magnetic_energy = sum(
        history[:, index] for index in metadata['history_transverse_me']
    )
    return {
        'mass_relative_drift_max': _relative_drift(history[:, 2], np.ones_like(times)),
        'parallel_magnetic_energy_scaled_relative_error_max': _relative_drift(
            history[:, metadata['history_parallel_me']], parallel/transverse**2
        ),
        'transverse_kinetic_energy_scaled_relative_error_max': _relative_drift(
            transverse_kinetic_energy, transverse**-2
        ),
        'transverse_magnetic_energy_scaled_relative_error_max': _relative_drift(
            transverse_magnetic_energy, parallel**-1
        ),
    }


def _measure_case(spec):
    snapshots = _snapshots(spec)
    history = _history_rows(spec)
    return {
        'profile': spec['profile'],
        'expansion_law': spec['law'],
        'axis': spec['axis'],
        'solver': spec['solver'],
        'mode': spec['mode'],
        'resolution': spec['resolution'],
        'dimensions': _dimensions(spec),
        'rates_x1_x2_x3': _rates(spec),
        'time': float(snapshots[-1]['Time']),
        'snapshot_count': len(snapshots),
        'history_row_count': int(history.shape[0]),
        'field_oracles': _field_oracle_metrics(snapshots, spec),
        'history_oracles': _history_oracle_metrics(history, spec),
    }


def _max_relative_difference(actual, expected):
    scale = max(float(np.max(np.abs(expected))), 1.0e-30)
    return float(np.max(np.abs(actual - expected)))/scale


def _control_parity():
    mode_on = _spec('static_on', 'exponential', 0.0, 'x1', 128)
    mode_off = _spec('static_off_control', 'exponential', 0.0, 'x1', 128, mode='off')
    mode_on_snapshots = _snapshots(mode_on)
    mode_off_snapshots = _snapshots(mode_off)
    if len(mode_on_snapshots) != len(mode_off_snapshots):
        raise RuntimeError('Static control snapshot count mismatch')
    field_error = 0.0
    for actual, expected in zip(mode_on_snapshots, mode_off_snapshots):
        for name in ['dens', 'velx', 'vely', 'velz', 'eint',
                     'bcc1', 'bcc2', 'bcc3']:
            field_error = max(
                field_error,
                _max_relative_difference(actual[name], expected[name]),
            )
    history_error = _max_relative_difference(
        _history_rows(mode_on), _history_rows(mode_off)
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
    for profile, (law, rate) in _LITERAL_PROFILES.items():
        specs = [_spec(profile, law, rate, 'x1', resolution)
                 for resolution in _RESOLUTIONS]
        observations[profile] = {}
        for metric in metrics:
            residuals = [
                cases[_case_id(spec)]['field_oracles'][metric] for spec in specs
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
        and history['mass_relative_drift_max'] < _PREPARATION_CONSTANT_HISTORY_BOUND
        and history[
            'parallel_magnetic_energy_scaled_relative_error_max'
        ] < _PREPARATION_CONSTANT_HISTORY_BOUND
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
    specs = _successful_specs()
    for spec in specs:
        _run_case(spec)
        cases[_case_id(spec)] = _measure_case(spec)
    _RESULTS['execution_contract'] = {
        'direct_serial_host_only': True,
        'mpi_launcher': 'not_used',
        'slurm_launcher': 'not_used',
        'frontier_execution': 'not_used',
        'errs_dat_oracle': 'not_used',
    }
    _RESULTS['successful_case_count'] = len(specs)
    _RESULTS['cases'] = cases
    _RESULTS['roe_disposition'] = _run_roe_rejection_probe()
    _RESULTS['static_mode_on_off_parity'] = _control_parity()
    _RESULTS['literal_profile_convergence_observations'] = (
        _convergence_observations(cases)
    )


def analyze():
    logger.debug('Analyzing test ' + __name__)
    for name, case in _RESULTS['cases'].items():
        logger.info('%s metrics=%s', name, case)
    logger.info('Roe disposition=%s', _RESULTS['roe_disposition'])
    logger.info(
        'static mode-on/off parity=%s', _RESULTS['static_mode_on_off_parity']
    )
    parity = _RESULTS['static_mode_on_off_parity']
    return (
        all(_case_passes(case) for case in _RESULTS['cases'].values())
        and _RESULTS['roe_disposition']['status'] == 'pass_fail_closed_unsupported'
        and parity['mhd_w_bcc_relative_difference_max'] <= _PREPARATION_PARITY_BOUND
        and parity['mhd_history_relative_difference_max'] <= _PREPARATION_PARITY_BOUND
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    run()
    if not analyze():
        raise SystemExit('pic_mhd_expanding_box_cpaw_history_preparation: FAIL')
    print(json.dumps(_RESULTS, indent=2, sort_keys=True))
