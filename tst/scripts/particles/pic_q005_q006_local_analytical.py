"""Bounded serial-host analytical characterization for readiness Q-005/Q-006.

The reused Athena fixtures are engineering proxies, not paper-reproduction
decks.  This test validates its analytical extractors independently and records
the proxy mismatch so the local evidence cannot be mistaken for claim closure.
"""

import glob
import json
import logging
import math
import os
import subprocess

import numpy as np
import scripts.utils.athena as athena
from scripts.particles.pic_analysis_utils import fit_exponential_growth
from scripts.particles.pic_analysis_utils import fit_exponential_growth_windowed

import sys
sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_BELL_DECK = 'tests/pic_bell_growth_proxy.athinput'
_OSCILLATION_DECKS = {
    'uniform': 'tests/pic_multispecies_osc_uniform.athinput',
    'smr': 'tests/pic_multispecies_osc_smr.athinput',
    'amr_proxy': 'tests/pic_multispecies_osc_amr_proxy.athinput',
}
_RESULTS = {}

# Deck-bound normalization for the reused engineering proxies.  The Q-005
# paper deck must still replace this proxy and bind the publication parameters.
_BELL_LENGTH = 16.0
_BELL_PPC = 8.0
_BELL_QSCALE = 0.1
_BELL_CHARGE = 1.0
_BELL_DRIFT = 0.8
_BELL_RHO = 1.0
_BELL_BPAR = 1.0
_BELL_VA = 1.0

# The oscillation deck has ppc=1 split round-robin across two species, rho=1,
# |q|/m=1, and B_parallel=1 from the linear-wave initializer.
_OSC_OMEGA_C = 1.0
_OSC_SPECIES_MASS = 1.0
_OSC_SPECIES_DENSITY = 0.5
_OSC_RHO = 1.0


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path(input_rel):
    return '../../' + athena.athena_rel_path + 'inputs/' + input_rel


def _remove_outputs(basename):
    output_dir = os.path.join(_athena_exe_dir(), 'bin')
    os.makedirs(output_dir, exist_ok=True)
    pattern = os.path.join(output_dir, basename + '.*.bin')
    for fname in glob.glob(pattern):
        os.remove(fname)


def _run_case(label, input_rel, basename, arguments=None):
    command = [
        './athena', '-i', _athena_input_path(input_rel),
        'job/basename=' + basename,
    ]
    if arguments:
        command.extend(arguments)
    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed for ' + label + '\n' + output)


def _check_close(label, measured, expected, abs_tol):
    error = abs(measured - expected)
    logger.info('%s measured=% .8e expected=% .8e error=% .8e tolerance=% .8e',
                label, measured, expected, error, abs_tol)
    return error <= abs_tol


def _check_lower(label, measured, lower):
    logger.info('%s measured=% .8e lower=% .8e margin=% .8e',
                label, measured, lower, measured - lower)
    return measured >= lower


def _check_upper(label, measured, upper):
    logger.info('%s measured=% .8e upper=% .8e margin=% .8e',
                label, measured, upper, measured - upper)
    return measured <= upper


def _bell_no_hall_oracle(k, k0, va):
    """Return the standard unstable no-Hall Bell branch for 0 < k < k0."""
    gamma2 = va * va * k * (k0 - k)
    return {
        'omega_real': 0.0,
        'gamma': math.sqrt(max(0.0, gamma2)),
        'unstable': bool(gamma2 > 0.0),
    }


def _oscillation_angular_frequency(omega_c, species_mass, species_density, rho):
    """Section 5.3 expression: Omega * sqrt(1 + 2 m_e n0 / rho)."""
    return omega_c * math.sqrt(1.0 + 2.0 * species_mass * species_density / rho)


def _single_tone_fit(time, signal, fmin=0.02, fmax=0.65, nsample=6000):
    """Fit a single sinusoid on nonuniform samples using a bounded grid."""
    t = np.asarray(time, dtype=float)
    y = np.asarray(signal, dtype=float)
    if t.size != y.size or t.size < 12:
        raise RuntimeError('insufficient samples for bounded frequency fit')

    y_centered = y - np.mean(y)
    variance = float(np.dot(y_centered, y_centered))
    if variance <= 0.0:
        raise RuntimeError('frequency fit requires a non-constant signal')

    best = None
    for freq in np.linspace(fmin, fmax, nsample):
        phase = 2.0 * np.pi * freq * t
        design = np.column_stack((np.ones(t.size), np.sin(phase), np.cos(phase)))
        coeff, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        residual = y - design.dot(coeff)
        rss = float(np.dot(residual, residual))
        if best is None or rss < best['rss']:
            best = {
                'frequency_cycles': float(freq),
                'omega': float(2.0 * np.pi * freq),
                'amplitude': float(math.hypot(coeff[1], coeff[2])),
                'r2': float(1.0 - rss / variance),
                'rss': rss,
            }
    return best


def _synthetic_oracle_checks():
    k0 = 0.64
    va = 1.0
    midpoint = _bell_no_hall_oracle(0.5 * k0, k0, va)
    left = _bell_no_hall_oracle(0.25 * k0, k0, va)
    right = _bell_no_hall_oracle(0.75 * k0, k0, va)

    omega = _oscillation_angular_frequency(
        _OSC_OMEGA_C, _OSC_SPECIES_MASS, _OSC_SPECIES_DENSITY, _OSC_RHO
    )
    time = np.linspace(0.0, 20.0, 401)
    signal = 3.0 + 0.25 * np.sin(omega * time + 0.37)
    tone = _single_tone_fit(time, signal)

    return {
        'bell': {
            'k0': k0,
            'midpoint_gamma': midpoint['gamma'],
            'expected_midpoint_gamma': 0.5 * k0 * va,
            'left_gamma': left['gamma'],
            'right_gamma': right['gamma'],
        },
        'oscillation': {
            'expected_omega': omega,
            'expected_frequency_cycles': omega / (2.0 * np.pi),
            'measured_frequency_cycles': tone['frequency_cycles'],
            'single_tone_r2': tone['r2'],
        },
    }


def _mode1_complex(dataset, field):
    values = np.asarray(dataset[field], dtype=float)
    x_mode = np.mean(values, axis=(0, 1))
    fluctuation = x_mode - np.mean(x_mode)
    x1f = np.asarray(dataset['x1f'], dtype=float)
    x1v = np.asarray(dataset['x1v'], dtype=float)
    length = float(x1f[-1] - x1f[0])
    phase = np.exp(-2.0j * np.pi * (x1v - x1f[0]) / length)
    return np.sum(fluctuation * phase) / x_mode.size


def _load_bell_metrics(basename, fit_coupled):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.mhd_bcc.*.bin')
    files = sorted(glob.glob(pattern))
    if len(files) < 8:
        raise RuntimeError('Not enough mhd_bcc outputs for Bell fit: ' + basename)

    times = []
    bperp = []
    right = []
    left = []
    for fname in files:
        data = bin_convert.read_binary_as_athdf(fname)
        b2_k = _mode1_complex(data, 'bcc2')
        b3_k = _mode1_complex(data, 'bcc3')
        times.append(float(data['Time']))
        bperp.append(math.hypot(abs(b2_k), abs(b3_k)))
        right.append(0.5 * (b2_k + 1.0j * b3_k))
        left.append(0.5 * (b2_k - 1.0j * b3_k))

    t = np.asarray(times, dtype=float)
    b = np.asarray(bperp, dtype=float)
    floor = max(1.0e-30, 1.0e-6 * np.max(b))
    if fit_coupled:
        fit = fit_exponential_growth_windowed(
            t, b, min_points=8, floor=floor, min_growth_factor=4.0
        )
        mask = (t >= fit['tmin']) & (t <= fit['tmax'])
    else:
        gamma, _, r2 = fit_exponential_growth(
            t, np.maximum(b, floor), float(t[1]), float(t[-1]), floor=floor
        )
        fit = {'gamma': gamma, 'r2': r2, 'tmin': float(t[1]), 'tmax': float(t[-1])}
        mask = (t >= t[1]) & (t <= t[-1])

    right = np.asarray(right)
    left = np.asarray(left)
    dominant = right if abs(right[-1]) >= abs(left[-1]) else left
    phase = np.unwrap(np.angle(dominant))
    phase_slope, _ = np.polyfit(t[mask], phase[mask], 1)

    return {
        'outputs': len(files),
        'time_span': float(t[-1] - t[0]),
        'gamma': float(fit['gamma']),
        'growth_r2': float(fit['r2']),
        'growth_ratio': float(b[-1] / max(b[0], 1.0e-30)),
        'omega_real_measured': float(abs(phase_slope)),
        'fit_tmin': float(fit['tmin']),
        'fit_tmax': float(fit['tmax']),
    }


def _integrate_quantity(dataset, quantity):
    dx1 = np.diff(dataset['x1f'])
    dx2 = np.diff(dataset['x2f'])
    dx3 = np.diff(dataset['x3f'])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    return float(np.sum(dataset[quantity] * dvol))


def _zero_crossings(signal):
    centered = np.asarray(signal, dtype=float) - np.mean(signal)
    signs = np.sign(centered)
    signs = signs[signs != 0.0]
    if signs.size < 2:
        return 0
    return int(np.count_nonzero(signs[1:] != signs[:-1]))


def _load_oscillation_metrics(basename, fit_frequency):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.mhd_u_m2.*.bin')
    files = sorted(glob.glob(pattern))
    if len(files) < 12:
        raise RuntimeError('Not enough mhd_u_m2 outputs: ' + basename)

    times = []
    momentum = []
    energy = []
    for m2_file in files:
        cycle = os.path.basename(m2_file).split('.')[-2]
        e_file = os.path.join(
            _athena_exe_dir(), 'bin', basename + '.mhd_u_e.' + cycle + '.bin'
        )
        m2_data = bin_convert.read_binary_as_athdf(m2_file)
        e_data = bin_convert.read_binary_as_athdf(e_file)
        times.append(float(m2_data['Time']))
        momentum.append(_integrate_quantity(m2_data, 'mom2'))
        energy.append(_integrate_quantity(e_data, 'ener'))

    t = np.asarray(times, dtype=float)
    m2 = np.asarray(momentum, dtype=float)
    ener = np.asarray(energy, dtype=float)
    finite = bool(np.all(np.isfinite(m2)) and np.all(np.isfinite(ener)))
    metrics = {
        'outputs': len(files),
        'time_span': float(t[-1] - t[0]),
        'finite': finite,
        'amplitude': float(np.max(m2) - np.min(m2)),
        'zero_crossings': _zero_crossings(m2),
        'energy_drift': float(np.max(np.abs(ener - ener[0])) / max(abs(ener[0]), 1.0)),
    }
    if fit_frequency:
        metrics['single_tone'] = _single_tone_fit(t, m2)
    return metrics


def _bell_proxy_oracle():
    current_density = _BELL_PPC * _BELL_QSCALE * _BELL_CHARGE * _BELL_DRIFT
    k0 = abs(current_density * _BELL_BPAR) / (_BELL_RHO * _BELL_VA * _BELL_VA)
    k = 2.0 * np.pi / _BELL_LENGTH
    branch = _bell_no_hall_oracle(k, k0, _BELL_VA)
    branch.update({'current_density': current_density, 'k': k, 'k0': k0})
    return branch


def _oscillation_proxy_oracle():
    omega = _oscillation_angular_frequency(
        _OSC_OMEGA_C, _OSC_SPECIES_MASS, _OSC_SPECIES_DENSITY, _OSC_RHO
    )
    return {
        'omega': omega,
        'frequency_cycles': omega / (2.0 * np.pi),
    }


def _summary():
    bell_oracle = _RESULTS['bell']['deck_bound_standard_no_hall_oracle']
    bell_measured = _RESULTS['bell']['coupled']
    osc_oracle = _RESULTS['oscillation']['deck_bound_section_5_3_oracle']
    uniform = _RESULTS['oscillation']['uniform']['single_tone']
    smr = _RESULTS['oscillation']['smr']['single_tone']
    return {
        'classification': 'bounded_serial_host_characterization_not_claim_closure',
        'q005': {
            'fixture': 'engineering_proxy',
            'oracle_gamma': bell_oracle['gamma'],
            'measured_gamma': bell_measured['gamma'],
            'measured_omega_real': bell_measured['omega_real_measured'],
            'paper_qualified': False,
        },
        'q006': {
            'fixture': 'engineering_proxy',
            'oracle_frequency_cycles': osc_oracle['frequency_cycles'],
            'uniform_frequency_cycles': uniform['frequency_cycles'],
            'uniform_single_tone_r2': uniform['r2'],
            'smr_frequency_cycles': smr['frequency_cycles'],
            'smr_single_tone_r2': smr['r2'],
            'uniform_vs_smr_abs_frequency_difference': abs(
                uniform['frequency_cycles'] - smr['frequency_cycles']
            ),
            'amr_proxy_frequency_qualified': False,
            'paper_qualified': False,
        },
    }


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    _RESULTS.clear()
    _RESULTS['synthetic_oracles'] = _synthetic_oracle_checks()

    bell_args = [
        'particles/couple_moments_momentum_to_mhd=false',
        'particles/couple_moments_energy_to_mhd=false',
    ]
    bell = {}
    for tag, coupled in [('uncoupled', False), ('coupled', True)]:
        base = 'pic_q005_local_' + tag
        _remove_outputs(base)
        _run_case(
            'q005_' + tag, _BELL_DECK, base,
            ['particles/couple_moments_to_mhd=' + str(coupled).lower()] + bell_args,
        )
        bell[tag] = _load_bell_metrics(base, fit_coupled=coupled)
    bell['deck_bound_standard_no_hall_oracle'] = _bell_proxy_oracle()
    _RESULTS['bell'] = bell

    oscillation = {}
    for tag, input_rel in _OSCILLATION_DECKS.items():
        base = 'pic_q006_local_' + tag
        _remove_outputs(base)
        _run_case('q006_' + tag, input_rel, base)
        oscillation[tag] = _load_oscillation_metrics(
            base, fit_frequency=(tag != 'amr_proxy')
        )
    oscillation['deck_bound_section_5_3_oracle'] = _oscillation_proxy_oracle()
    _RESULTS['oscillation'] = oscillation

    print(json.dumps(_summary(), indent=2, sort_keys=True))


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    synthetic_bell = _RESULTS['synthetic_oracles']['bell']
    ok = _check_close(
        'synthetic_bell:midpoint_gamma',
        synthetic_bell['midpoint_gamma'],
        synthetic_bell['expected_midpoint_gamma'],
        1.0e-14,
    ) and ok
    ok = _check_close(
        'synthetic_bell:symmetric_gamma',
        synthetic_bell['left_gamma'],
        synthetic_bell['right_gamma'],
        1.0e-14,
    ) and ok

    synthetic_osc = _RESULTS['synthetic_oracles']['oscillation']
    ok = _check_close(
        'synthetic_oscillation:frequency_cycles',
        synthetic_osc['measured_frequency_cycles'],
        synthetic_osc['expected_frequency_cycles'],
        2.0e-4,
    ) and ok
    ok = _check_lower(
        'synthetic_oscillation:single_tone_r2',
        synthetic_osc['single_tone_r2'],
        0.999999,
    ) and ok

    bell = _RESULTS['bell']
    oracle_gamma = bell['deck_bound_standard_no_hall_oracle']['gamma']
    ok = _check_upper('q005_proxy:uncoupled_abs_gamma', abs(bell['uncoupled']['gamma']),
                      5.0e-2) and ok
    ok = _check_lower('q005_proxy:coupled_gamma', bell['coupled']['gamma'], 1.0e-1) and ok
    ok = _check_lower('q005_proxy:coupled_growth_r2', bell['coupled']['growth_r2'],
                      0.5) and ok
    ok = _check_lower('q005_proxy:oracle_gamma_mismatch',
                      abs(bell['coupled']['gamma'] - oracle_gamma), 1.0) and ok

    oscillation = _RESULTS['oscillation']
    oracle_frequency = oscillation['deck_bound_section_5_3_oracle']['frequency_cycles']
    for tag in ['uniform', 'smr']:
        trace = oscillation[tag]
        ok = bool(trace['finite']) and ok
        ok = _check_lower('q006_proxy:' + tag + ':amplitude',
                          trace['amplitude'], 1.0) and ok
        ok = _check_lower('q006_proxy:' + tag + ':zero_crossings',
                          trace['zero_crossings'], 4.0) and ok
        ok = _check_lower('q006_proxy:' + tag + ':oracle_frequency_mismatch',
                          abs(trace['single_tone']['frequency_cycles'] -
                              oracle_frequency), 5.0e-2) and ok

    ok = _check_close(
        'q006_proxy:uniform_vs_smr_frequency',
        oscillation['uniform']['single_tone']['frequency_cycles'],
        oscillation['smr']['single_tone']['frequency_cycles'],
        2.0e-2,
    ) and ok

    amr = oscillation['amr_proxy']
    ok = bool(amr['finite']) and ok
    ok = _check_lower('q006_proxy:amr_proxy:amplitude', amr['amplitude'], 5.0e-2) and ok
    ok = _check_upper('q006_proxy:amr_proxy:time_span', amr['time_span'], 2.0) and ok

    print(json.dumps(_summary(), indent=2, sort_keys=True))
    return ok


if __name__ == '__main__':
    run()
    if not analyze():
        raise SystemExit(1)
