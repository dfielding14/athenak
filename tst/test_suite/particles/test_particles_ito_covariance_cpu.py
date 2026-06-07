"""Finite-step multidimensional covariance tests for Ito-2 flux tracers."""

from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from read_prtcl_thermo_history import read_history  # noqa: E402


INPUT = str(ROOT / "inputs/particles/ito_tracers.athinput")
RUN_ROOT = Path("run_particles_ito_covariance")
NPART = 65536
CELL_WIDTHS = np.array([1.0 / 32.0, 1.0 / 8.0, 1.0 / 8.0])


CASES = [
    pytest.param(
        "positive_2d", (1.0, 0.1, 0.0), 2, 0.95, 0.01, 0.99,
        id="positive-2d",
    ),
    pytest.param(
        "mixed_sign_2d", (1.0, -0.1, 0.0), 2, 0.95, 0.01, 0.99,
        id="mixed-sign-2d",
    ),
    pytest.param(
        "stationary_axis_2d", (1.0, 0.0, 0.0), 2, 0.95, 0.01, 0.99,
        id="stationary-axis-2d",
    ),
    pytest.param(
        "near_limit_3d", (1.0, 4.0, -4.0), 3, 0.95, 0.0, 0.99,
        id="near-outgoing-one-3d",
    ),
    pytest.param(
        "rank_two_3d", (1.0, 4.0, -4.0), 3, 0.95, 0.0, 1.0,
        id="correlated-rank-deficient-3d",
    ),
    pytest.param(
        "rank_one_3d", (1.0, 0.0, 0.0), 3, 0.99, 0.01, 0.99,
        id="rank-deficient-3d",
    ),
    pytest.param(
        "zero_flow_2d", (0.0, 0.0, 0.0), 2, 0.95, 1.0, 0.99,
        id="zero-flow-2d",
    ),
]


def _run_case(
    name,
    velocities,
    ndim,
    cfl,
    sound_speed,
    probability_target,
    npart=NPART,
    suffix="",
):
    run_dir = RUN_ROOT / f"{name}{suffix}"
    shutil.rmtree(run_dir, ignore_errors=True)
    run_dir.mkdir(parents=True)
    basename = f"ito_cov_{name}{suffix}"
    flags = [
        "-d",
        str(run_dir),
        f"job/basename={basename}",
        f"time/cfl_number={cfl}",
        f"hydro/iso_sound_speed={sound_speed}",
        f"particles/ito_probability_target={probability_target}",
        f"tracer_seed1/count_per_event={npart}",
        f"problem/vx0={velocities[0]}",
        f"problem/vy0={velocities[1]}",
        f"problem/vz0={velocities[2]}",
        "output2/dt=-1.0",
        "output3/dt=-1.0",
    ]
    if ndim == 3:
        flags.extend(["mesh/nx3=8", "meshblock/nx3=8"])
    assert testutils.run(INPUT, flags)
    history = read_history(
        run_dir / "prtcl_thermo_history" / f"{basename}.prtcl_thermo_history.thp"
    )
    return _displacements(history), np.max(history["time"])


def _displacements(history):
    first_cycle = np.min(history["cycle"])
    last_cycle = np.max(history["cycle"])
    assert first_cycle == 0
    assert last_cycle == 1

    initial = history["cycle"] == first_cycle
    final = history["cycle"] == last_cycle
    initial_order = np.argsort(history["tag"][initial])
    final_order = np.argsort(history["tag"][final])
    initial_tags = history["tag"][initial][initial_order]
    final_tags = history["tag"][final][final_order]
    assert np.unique(initial_tags).size == initial_tags.size
    assert np.unique(final_tags).size == final_tags.size
    np.testing.assert_array_equal(initial_tags, final_tags)

    displacement = np.column_stack(
        [
            history[component][final][final_order]
            - history[component][initial][initial_order]
            for component in ("x1", "x2", "x3")
        ]
    )
    return (displacement + 0.5) % 1.0 - 0.5


def _expected_moments(velocities, dt, ndim):
    velocity = np.asarray(velocities)
    active = np.arange(ndim)
    mean = velocity * dt
    courant = np.zeros(3)
    courant[active] = np.abs(velocity[active]) * dt / CELL_WIDTHS[active]
    covariance = np.diag(CELL_WIDTHS**2 * courant) - np.outer(mean, mean)
    covariance[ndim:, :] = 0.0
    covariance[:, ndim:] = 0.0
    return mean, covariance, np.sum(courant)


def test_spatially_varying_deterministic_drift_has_no_interpolation_noise():
    """Interpolating central covariance does not turn drift gradients into noise."""
    means = np.array(
        [
            [CELL_WIDTHS[0], 0.0, 0.0],
            [-CELL_WIDTHS[0], 0.0, 0.0],
        ]
    )
    covariance = np.zeros((2, 3, 3))
    weights = np.array([0.5, 0.5])

    interpolated_mean = weights @ means
    interpolated_covariance = np.tensordot(weights, covariance, axes=1)
    np.testing.assert_array_equal(interpolated_mean, np.zeros(3))
    np.testing.assert_array_equal(interpolated_covariance, np.zeros((3, 3)))

    raw_second = covariance + np.einsum("ni,nj->nij", means, means)
    old_covariance = (
        np.tensordot(weights, raw_second, axes=1)
        - np.outer(interpolated_mean, interpolated_mean)
    )
    assert old_covariance[0, 0] == CELL_WIDTHS[0] ** 2
    assert np.trace(old_covariance) > 0.0

    source = (ROOT / "src/particles/particles_lagrangian_ito.cpp").read_text(
        encoding="utf-8"
    )
    push_start = source.index("TaskStatus Particles::PushIto2")
    push = source[push_start:]
    assert "q[n] += weight*coeff(m,ITO_Q11+n" in push
    assert "raw_second" not in push


def test_anisotropic_factor_tolerance_uses_local_scales():
    """A strong direction must not mask a negative weak-direction pivot."""
    relative_tolerance = 2.0e-12
    strong_variance = 1.0
    weak_variance = 1.0e-18

    old_global_tolerance = relative_tolerance * strong_variance
    assert weak_variance < old_global_tolerance

    weak_pivot_tolerance = relative_tolerance * weak_variance
    assert weak_variance > weak_pivot_tolerance
    assert -weak_variance < -weak_pivot_tolerance

    source = (ROOT / "src/particles/particles_lagrangian_ito.cpp").read_text(
        encoding="utf-8"
    )
    factor_start = source.index("bool ItoFactorCovariance")
    factor_end = source.index("void FatalIto", factor_start)
    factor = source[factor_start:factor_end]
    assert "pivot_scale = fmax(fabs(a[pivot][pivot]), pivot_correction)" in factor
    assert "residual_scale = fmax(diagonal_scale, residual_correction)" in factor
    assert "rel_tol*scale" not in factor


def _assert_sample_mean(samples, expected, variance, label):
    measured = np.mean(samples)
    standard_error = np.sqrt(max(variance, 0.0) / samples.size)
    tolerance = 5.0 * standard_error + 128.0 * np.finfo(samples.dtype).eps
    residual = (measured - expected) / standard_error if standard_error > 0.0 else 0.0
    print(
        f"{label}: expected={expected:.16e}, measured={measured:.16e}, "
        f"se={standard_error:.3e}, residual={residual:.3f}"
    )
    assert abs(measured - expected) <= tolerance


def _assert_sample_covariance(dx, mean, expected, i, j, label):
    products = (dx[:, i] - mean[i]) * (dx[:, j] - mean[j])
    measured = np.mean(products)
    standard_error = np.std(products, ddof=1) / np.sqrt(products.size)
    scale = max(CELL_WIDTHS[i] * CELL_WIDTHS[j], abs(expected), 1.0e-30)
    tolerance = 5.0 * standard_error + 256.0 * np.finfo(dx.dtype).eps * scale
    residual = (measured - expected) / standard_error if standard_error > 0.0 else 0.0
    print(
        f"{label}: expected={expected:.16e}, measured={measured:.16e}, "
        f"se={standard_error:.3e}, residual={residual:.3f}"
    )
    assert abs(measured - expected) <= tolerance


@pytest.mark.parametrize(
    "name,velocities,ndim,cfl,sound_speed,probability_target", CASES
)
def test_ito2_matches_full_finite_step_covariance_cpu(
    name, velocities, ndim, cfl, sound_speed, probability_target
):
    """All finite-step means and covariance entries match the MC jump kernel."""
    try:
        dx, dt = _run_case(
            name, velocities, ndim, cfl, sound_speed, probability_target
        )
        assert dx.shape == (NPART, 3)
        mean, covariance, outgoing = _expected_moments(velocities, dt, ndim)
        precision = np.finfo(dx.dtype).eps
        probability_tolerance = max(2.0e-12, 8.0 * precision)
        print(f"{name}: dt={dt:.16e}, outgoing_probability={outgoing:.16e}")
        assert outgoing <= 1.0 + probability_tolerance
        if name == "near_limit_3d":
            assert outgoing > 0.98
        if name == "rank_two_3d":
            assert abs(outgoing - 1.0) < probability_tolerance

        for i in range(ndim):
            _assert_sample_mean(
                dx[:, i], mean[i], covariance[i, i], f"{name} mean[{i}]"
            )
        for i in range(ndim):
            for j in range(i, ndim):
                _assert_sample_covariance(
                    dx, mean, covariance[i, j], i, j, f"{name} cov[{i},{j}]"
                )

        covariance_scale = np.max(np.abs(covariance))
        rank_tolerance = max(1.0e-14, 64.0 * precision * covariance_scale)
        eigenvalues = np.linalg.eigvalsh(covariance[:ndim, :ndim])
        assert np.min(eigenvalues) >= -rank_tolerance
        if name == "rank_two_3d":
            values, vectors = np.linalg.eigh(covariance)
            assert np.count_nonzero(values > rank_tolerance) == 2
            null_projection = np.sum((dx - mean) * vectors[:, 0], axis=1)
            projection_tolerance = max(
                1.0e-12, 32.0 * precision * np.max(CELL_WIDTHS)
            )
            assert np.max(np.abs(null_projection)) < projection_tolerance
        if name == "rank_one_3d":
            assert np.count_nonzero(eigenvalues > rank_tolerance) == 1
            assert np.max(np.abs(dx[:, 1:])) == 0.0
        if name == "zero_flow_2d":
            assert np.max(np.abs(dx)) == 0.0
        if ndim == 2:
            assert np.max(np.abs(dx[:, 2])) == 0.0
        if velocities[1] == 0.0:
            assert np.max(np.abs(dx[:, 1])) == 0.0
    finally:
        shutil.rmtree(RUN_ROOT, ignore_errors=True)


def test_ito2_covariance_rng_is_deterministic_cpu():
    """Fixed tags, cycle, and seed give identical correlated displacements."""
    try:
        args = ("deterministic", (1.0, -0.1, 0.0), 2, 0.95, 0.01, 0.99)
        first, first_dt = _run_case(*args, npart=8192, suffix="_a")
        second, second_dt = _run_case(*args, npart=8192, suffix="_b")
        assert first_dt == second_dt
        np.testing.assert_array_equal(first, second)
    finally:
        shutil.rmtree(RUN_ROOT, ignore_errors=True)
