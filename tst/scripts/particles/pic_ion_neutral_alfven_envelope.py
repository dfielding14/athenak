import logging
import math
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from scripts.particles import pic_ion_neutral_friction_smoke as smoke

logger = logging.getLogger('athena' + __name__[7:])

_NLIM = 8
_SCAN = [
    ('off', 'off', 0.0),
    ('nu_0p2', 'ion_neutral_friction', 0.2),
    ('nu_0p7', 'ion_neutral_friction', 0.7),
    ('nu_1p4', 'ion_neutral_friction', 1.4),
]
_RESULTS = {}


def _transverse_fluctuation_amplitude(snapshot):
    dens = snapshot['dens']
    vel2 = snapshot['mom2']/dens
    vel3 = snapshot['mom3']/dens
    vel2 = vel2 - np.mean(vel2)
    vel3 = vel3 - np.mean(vel3)
    return float(np.sqrt(np.mean(vel2*vel2 + vel3*vel3)))


def _snapshot_series(basename):
    snapshots = []
    for path in smoke._output_files(basename):
        snapshot = smoke.bin_convert.read_binary_as_athdf(path)
        row = {
            'time': float(snapshot['Time']),
            'amplitude': _transverse_fluctuation_amplitude(snapshot),
        }
        if snapshots and math.isclose(
                row['time'], snapshots[-1]['time'], rel_tol=0.0, abs_tol=1.0e-12):
            if not math.isclose(
                    row['amplitude'], snapshots[-1]['amplitude'],
                    rel_tol=0.0, abs_tol=1.0e-12):
                raise RuntimeError('Duplicate-time snapshots disagree for ' + basename)
            continue
        snapshots.append(row)
    return snapshots


def _run_case(label, damping_mode, collision_rate):
    basename = 'pic_ion_neutral_alfven_envelope_' + label
    smoke._remove_outputs(basename)
    command = [
        './athena', '-i', smoke._athena_input_path(),
        'job/basename=' + basename,
        'time/nlim=' + str(_NLIM),
        'particles/pic_wave_damping_mode=' + damping_mode,
        'particles/pic_ion_neutral_collision_rate=' + str(collision_rate),
    ]
    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=smoke._athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + label + '\n' + output)
    for token in [
            'physical_mode=extended_mhd_pic',
            'wave_damping=' + damping_mode,
            'nu_in=']:
        if token not in output:
            raise RuntimeError(label + ' missing runtime token: ' + token)
    return _snapshot_series(basename)


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _RESULTS.clear()
    for label, damping_mode, collision_rate in _SCAN:
        _RESULTS[label] = _run_case(label, damping_mode, collision_rate)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    labels = [label for label, _mode, _rate in _SCAN]
    collision_rates = np.asarray([rate for _label, _mode, rate in _SCAN])
    times = np.asarray([row['time'] for row in _RESULTS['off']])
    amplitudes = {
        label: np.asarray([row['amplitude'] for row in _RESULTS[label]])
        for label in labels
    }
    final_amplitudes = np.asarray([amplitudes[label][-1] for label in labels])
    excess_rates = -np.log(final_amplitudes/final_amplitudes[0])/times[-1]
    rate_errors = np.abs(excess_rates - collision_rates)
    measured = {
        'times': times.tolist(),
        'final_amplitudes': dict(zip(labels, final_amplitudes.tolist())),
        'excess_damping_rates': dict(zip(labels, excess_rates.tolist())),
        'rate_errors': dict(zip(labels, rate_errors.tolist())),
    }
    logger.info('ion-neutral Alfven envelope metrics: %s', measured)
    _RESULTS['metrics'] = measured
    return (
        len(times) == _NLIM + 1
        and np.all(np.isfinite(times))
        and np.all(np.diff(times) > 0.0)
        and all(len(amplitudes[label]) == len(times) for label in labels)
        and all(np.all(np.isfinite(amplitudes[label])) for label in labels)
        and max(abs(amplitudes[label][0] - amplitudes['off'][0])
                for label in labels) <= 1.0e-12
        and amplitudes['off'][0] > 1.0e-6
        and all(np.all(np.diff(amplitudes[label]) < 0.0) for label in labels)
        and np.all(np.diff(final_amplitudes) < 0.0)
        and np.all(np.diff(excess_rates) > 0.0)
        and np.max(rate_errors) <= 8.0e-2
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    run()
    if not analyze():
        raise SystemExit('pic_ion_neutral_alfven_envelope: FAIL')
    print('pic_ion_neutral_alfven_envelope: PASS')
