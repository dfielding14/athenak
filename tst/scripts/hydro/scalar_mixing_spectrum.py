# Regression test for scalar_mixing velocity initialization.
#
# This test is only active when AthenaK is configured with
#   --cmake=-DPROBLEM=scalar_mixing
#
# It verifies that projection, 2D stream, and 3D Clebsch initialization produce
# the requested RMS velocity, near-zero mean flow, a discretely divergence-free
# staggered face transport field, and the expected in-band power-law slope.
# Because u0 is reconstructed from that face field by cell-centering averages,
# its FFT-based divergence is small but no longer roundoff-level. The spectrum
# check is performed on an ensemble average over a fixed seed set.

import glob
import json
import logging
import os
import sys

import numpy as np

import scripts.utils.athena as athena

sys.path.insert(0, '../vis/python')
import bin_convert_new  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_enabled = False
_SEEDS = [123, 321, 777, 2026, 12345]

_CASES = [
    {
        'name': 'projection_2d',
        'input': 'tests/scalar_mixing_projection_2d_regression.athinput',
        'basename': 'scalar_mix_proj_2d',
        'nlow': 2,
        'nhigh': 6,
        'expo': 1.6666667,
        'vrms': 1.0,
        'spectral_ndim': 2,
        'expect_vz_zero': True,
        'check_leakage': False,
        'slope_tol': 0.2,
        'check_exact_k_slope': True,
        'div_tol': 2.0e-2,
        'face_div_tol': 1.0e-10,
    },
    {
        'name': 'projection_3d',
        'input': 'tests/scalar_mixing_projection_3d_regression.athinput',
        'basename': 'scalar_mix_proj_3d',
        'nlow': 2,
        'nhigh': 6,
        'expo': 1.6666667,
        'vrms': 1.0,
        'spectral_ndim': 3,
        'expect_vz_zero': False,
        'check_leakage': False,
        'slope_tol': 0.2,
        'check_exact_k_slope': True,
        'div_tol': 2.0e-2,
        'face_div_tol': 1.0e-10,
    },
    {
        'name': 'stream_2d',
        'input': 'tests/scalar_mixing_stream_2d_regression.athinput',
        'basename': 'scalar_mix_stream_2d',
        'nlow': 2,
        'nhigh': 6,
        'expo': 1.6666667,
        'vrms': 1.0,
        'spectral_ndim': 2,
        'expect_vz_zero': True,
        'check_leakage': False,
        'slope_tol': 0.2,
        'check_exact_k_slope': True,
        'div_tol': 2.0e-2,
        'face_div_tol': 1.0e-10,
    },
    {
        'name': 'clebsch_3d',
        'input': 'tests/scalar_mixing_clebsch_3d_regression.athinput',
        'basename': 'scalar_mix_clebsch_3d',
        'nlow': 2,
        'nhigh': 10,
        'expo': 1.6666667,
        'alpha': 1.0/3.0,
        'phi_slope': 11.0/3.0,
        'velocity_slope': 5.0/3.0,
        'vrms': 1.0,
        'spectral_ndim': 3,
        'expect_vz_zero': False,
        'check_leakage': True,
        'slope_tol': 0.7,
        'check_exact_k_slope': False,
        'check_shell_slope': True,
        'div_tol': 5.0e-2,
        'face_div_tol': 1.0e-10,
    },
]


def _compiled_problem():
    cfg_path = os.path.join('build', 'config.hpp')
    if not os.path.exists(cfg_path):
        return None
    with open(cfg_path, 'r', encoding='utf-8') as fobj:
        for line in fobj:
            if line.startswith('#define PROBLEM_GENERATOR'):
                return line.split('"')[1]
    return None


def _run_basename(case, seed):
    return f"{case['basename']}_{seed}"


def _output_path(case, seed):
    pattern = os.path.join('build', 'src', 'bin',
                           _run_basename(case, seed) + '.hydro_w.*.bin')
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise RuntimeError('No hydro_w output found for '
                           + case['name'] + f' seed={seed}')
    return matches[-1]


def _diag_path(case, seed):
    return os.path.join('build', 'src',
                        _run_basename(case, seed) + '.turb_init_diag.json')


def _clear_old_outputs(case):
    for seed in _SEEDS:
        pattern = os.path.join('build', 'src', 'bin',
                               _run_basename(case, seed) + '.hydro_w.*.bin')
        for path in glob.glob(pattern):
            os.remove(path)
        diag_path = _diag_path(case, seed)
        if os.path.exists(diag_path):
            os.remove(diag_path)


def _load_diag(case, seed):
    diag_path = _diag_path(case, seed)
    if not os.path.exists(diag_path):
        raise RuntimeError('No diagnostics JSON found for '
                           + case['name'] + f' seed={seed}: {diag_path}')
    with open(diag_path, 'r', encoding='utf-8') as fobj:
        return json.load(fobj)


def _shell_bins(shape):
    nz, ny, nx = shape
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny if ny > 1 else np.array([0.0])
    kz = np.fft.fftfreq(nz) * nz if nz > 1 else np.array([0.0])
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing='ij')
    k2 = kx3**2 + ky3**2 + kz3**2
    shell = np.ceil(np.sqrt(k2) - 1.0e-12).astype(np.int64)
    return kx3, ky3, kz3, k2, shell


def _mode_group_means(energy, kmag, nlow, nhigh):
    mask = (kmag >= nlow) & (kmag <= nhigh) & (kmag > 0.0)
    groups = {}
    unique_k = np.unique(np.round(kmag[mask], 12))
    for kval in unique_k:
        kval_mask = mask & (np.abs(kmag - kval) < 1.0e-12)
        groups[float(kval)] = float(np.mean(energy[kval_mask]))
    return groups


def _fit_grouped_mode_slope(group_means):
    keys = sorted(group_means)
    if len(keys) < 2:
        return np.nan
    xs = np.log(np.asarray(keys, dtype=np.float64))
    ys = np.log(np.asarray([group_means[key] for key in keys], dtype=np.float64))
    A = np.vstack([xs, np.ones_like(xs)]).T
    slope, _ = np.linalg.lstsq(A, ys, rcond=None)[0]
    return slope


def _fit_shell_slope(shell_energy, nlow, nhigh):
    kvals = np.arange(nlow, nhigh + 1, dtype=np.float64)
    band = np.asarray(shell_energy[nlow:nhigh + 1], dtype=np.float64)
    mask = band > 0.0
    if np.count_nonzero(mask) < 2:
        return np.nan
    xs = np.log(kvals[mask])
    ys = np.log(band[mask])
    A = np.vstack([xs, np.ones_like(xs)]).T
    slope, _ = np.linalg.lstsq(A, ys, rcond=None)[0]
    return slope


def _analyze_output(case, seed):
    data = bin_convert_new.read_binary_as_athdf(
        _output_path(case, seed), quantities=['velx', 'vely', 'velz'])
    vx = np.asarray(data['velx'], dtype=np.float64)
    vy = np.asarray(data['vely'], dtype=np.float64)
    vz = np.asarray(data['velz'], dtype=np.float64)

    means = np.array([vx.mean(), vy.mean(), vz.mean()])
    vrms = np.sqrt(np.mean(vx*vx + vy*vy + vz*vz))

    uhat_x = np.fft.fftn(vx)
    uhat_y = np.fft.fftn(vy)
    uhat_z = np.fft.fftn(vz)
    kx3, ky3, kz3, k2, shell = _shell_bins(vx.shape)
    energy = np.abs(uhat_x)**2 + np.abs(uhat_y)**2 + np.abs(uhat_z)**2
    nonzero = k2 > 0.0

    kdot_u = kx3*uhat_x + ky3*uhat_y + kz3*uhat_z
    div_num = np.sum(np.abs(kdot_u[nonzero])**2)
    div_den = np.sum(k2[nonzero] * energy[nonzero])
    div_ratio = np.sqrt(div_num / div_den) if div_den > 0.0 else 0.0

    shell_energy = np.bincount(shell.ravel(),
                               weights=energy.ravel(),
                               minlength=max(shell.max() + 1, 2*case['nhigh'] + 2))
    out_band = (np.sum(shell_energy[1:case['nlow']]) +
                np.sum(shell_energy[case['nhigh'] + 1:]))
    in_band = np.sum(shell_energy[case['nlow']:case['nhigh'] + 1])
    leakage = out_band / (in_band + out_band) if (in_band + out_band) > 0.0 else 0.0
    grouped_mode_means = _mode_group_means(energy, np.sqrt(k2), case['nlow'],
                                           case['nhigh'])

    return {
        'means': means,
        'vrms': vrms,
        'div_ratio': div_ratio,
        'leakage': leakage,
        'shell_energy': shell_energy,
        'shell_slope': _fit_shell_slope(shell_energy, case['nlow'], case['nhigh']),
        'grouped_mode_means': grouped_mode_means,
        'vz_max': np.max(np.abs(vz)),
    }


def run(**kwargs):
    del kwargs
    global _enabled
    problem = _compiled_problem()
    _enabled = (problem == 'scalar_mixing')
    if not _enabled:
        logger.info('Skipping %s because build/config.hpp reports PROBLEM=%s',
                    __name__, problem)
        return

    for case in _CASES:
        _clear_old_outputs(case)
        for seed in _SEEDS:
            athena.run(case['input'],
                       [f'job/basename={_run_basename(case, seed)}',
                        f'problem/turb_rseed={seed}'])


def analyze():
    if not _enabled:
        return True

    analyze_status = True
    for case in _CASES:
        per_seed = [_analyze_output(case, seed) for seed in _SEEDS]
        diag_per_seed = [_load_diag(case, seed) for seed in _SEEDS]
        grouped_mode_means = {}
        for metrics in per_seed:
            for kval, mean_energy in metrics['grouped_mode_means'].items():
                grouped_mode_means.setdefault(kval, []).append(mean_energy)
        grouped_mode_means = {
            kval: float(np.mean(samples))
            for kval, samples in grouped_mode_means.items()
        }

        exact_k_slope = _fit_grouped_mode_slope(grouped_mode_means)
        target_mode_slope = -(case['expo'] + case['spectral_ndim'] - 1.0)
        logger.info('%s: exact-k slope=%.6f target=%.6f max_div=%.3e max_leak=%.3e',
                    case['name'], exact_k_slope, target_mode_slope,
                    max(metrics['div_ratio'] for metrics in per_seed),
                    max(metrics['leakage'] for metrics in per_seed))
        shell_slope = float(np.nanmean([metrics['shell_slope'] for metrics in per_seed]))

        max_mean_err = max(np.max(np.abs(metrics['means'])) for metrics in per_seed)
        if max_mean_err > 1.0e-9:
            logger.warning('%s mean velocity too large across seeds: %.3e',
                           case['name'], max_mean_err)
            analyze_status = False

        max_vrms_rel_err = max(abs(metrics['vrms'] - case['vrms']) / case['vrms']
                               for metrics in per_seed)
        if max_vrms_rel_err > 1.0e-3:
            logger.warning('%s vrms mismatch across seeds: %.6e',
                           case['name'], max_vrms_rel_err)
            analyze_status = False

        max_div_ratio = max(metrics['div_ratio'] for metrics in per_seed)
        if max_div_ratio > case.get('div_tol', 2.0e-7):
            logger.warning('%s divergence ratio too large across seeds: %.3e',
                           case['name'], max_div_ratio)
            analyze_status = False

        if case['check_exact_k_slope']:
            slope_err = abs(exact_k_slope - target_mode_slope)
            if not np.isfinite(exact_k_slope) or slope_err > case['slope_tol']:
                logger.warning('%s exact-k slope mismatch: slope %.6f target %.6f',
                               case['name'], exact_k_slope, target_mode_slope)
                analyze_status = False

        if case.get('check_shell_slope', False):
            target_shell_slope = -float(case['velocity_slope'])
            slope_err = abs(shell_slope - target_shell_slope)
            if not np.isfinite(shell_slope) or slope_err > case['slope_tol']:
                logger.warning('%s shell slope mismatch: slope %.6f target %.6f',
                               case['name'], shell_slope, target_shell_slope)
                analyze_status = False

        if case['expect_vz_zero']:
            max_vz = max(metrics['vz_max'] for metrics in per_seed)
            if max_vz > 1.0e-12:
                logger.warning('%s expected velz=0 but max(|vz|)=%.3e',
                               case['name'], max_vz)
                analyze_status = False

        if case['check_leakage']:
            max_leakage = max(metrics['leakage'] for metrics in per_seed)
            if max_leakage > 5.0e-2:
                logger.warning('%s leakage too large across seeds: %.3e',
                               case['name'], max_leakage)
                analyze_status = False

        max_face_div = 0.0
        for seed, diag in zip(_SEEDS, diag_per_seed):
            meta = diag.get('metadata', {})
            spectra = diag.get('spectra', {})
            target = meta.get('target', {})
            if meta.get('construction') != case['name']:
                logger.warning('%s wrong construction metadata for seed=%d',
                               case['name'], seed)
                analyze_status = False
            if meta.get('scalar_face_velocity') is not True:
                logger.warning('%s missing scalar_face_velocity=true for seed=%d',
                               case['name'], seed)
                analyze_status = False
            face_div = float(meta.get('face_divergence_ratio', np.inf))
            max_face_div = max(max_face_div, face_div)
            if not np.isfinite(face_div) or face_div > case['face_div_tol']:
                logger.warning('%s face divergence ratio too large for seed=%d: %.3e',
                               case['name'], seed, face_div)
                analyze_status = False

            if case['name'] != 'clebsch_3d':
                continue

            for forbidden_key in ('fit', 'quality_gate', 'sector_response'):
                if forbidden_key in diag:
                    logger.warning('%s stale diagnostics key %s present for seed=%d',
                                   case['name'], forbidden_key, seed)
                    analyze_status = False
            if 'phi1_shell_energy' not in spectra or 'phi2_shell_energy' not in spectra:
                logger.warning('%s missing phi shell spectra for seed=%d',
                               case['name'], seed)
                analyze_status = False
                continue
            if 'velocity_shell_energy' in spectra:
                logger.warning('%s stale retained velocity shell spectrum present for seed=%d',
                               case['name'], seed)
                analyze_status = False
            if abs(float(target.get('alpha', np.nan)) - float(case['alpha'])) > 1.0e-12:
                logger.warning('%s wrong alpha metadata for seed=%d',
                               case['name'], seed)
                analyze_status = False
            if abs(float(target.get('phi_slope', np.nan)) - float(case['phi_slope'])) > 1.0e-12:
                logger.warning('%s wrong phi slope metadata for seed=%d',
                               case['name'], seed)
                analyze_status = False
            if abs(float(target.get('velocity_slope', np.nan)) -
                   float(case['velocity_slope'])) > 1.0e-12:
                logger.warning('%s wrong velocity slope metadata for seed=%d',
                               case['name'], seed)
                analyze_status = False

            for key in ('phi1_shell_energy', 'phi2_shell_energy'):
                shell = np.asarray(spectra[key], dtype=np.float64)
                ref = shell[case['nlow']]
                target_shell = np.array([
                    0.0 if kval == 0 else kval ** (-float(case['phi_slope']))
                    for kval in range(len(shell))
                ], dtype=np.float64)
                target_shell *= ref / target_shell[case['nlow']]
                rel_err = np.max(np.abs(shell[case['nlow']:case['nhigh'] + 1]
                                        - target_shell[case['nlow']:case['nhigh'] + 1])
                                 / target_shell[case['nlow']:case['nhigh'] + 1])
                if rel_err > 1.0e-10:
                    logger.warning('%s %s not exact power law for seed=%d: %.3e',
                                   case['name'], key, seed, rel_err)
                    analyze_status = False

        logger.info('%s: max_face_div=%.3e', case['name'], max_face_div)

    return analyze_status
