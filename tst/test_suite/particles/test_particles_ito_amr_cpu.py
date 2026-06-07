"""AMR tests for Ito mean and central-covariance coefficients."""

from pathlib import Path
import shutil
import sys

import numpy as np
import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
ITO_SOURCE = ROOT / "src/particles/particles_lagrangian_ito.cpp"
INPUT = str(ROOT / "tst/inputs/particles_ito_amr_moments.athinput")
RUN_DIR = Path("run_particles_ito_amr_moments")

sys.path.insert(0, str(ROOT / "scripts"))
from read_prtcl_thermo_history import read_history  # noqa: E402


def _central_coefficients(probabilities, cell_widths=np.ones(3)):
    probabilities = np.asarray(probabilities)
    plus = probabilities[0::2]
    minus = probabilities[1::2]
    mean = cell_widths * (plus - minus)
    raw_diagonal = cell_widths**2 * (plus + minus)
    covariance = np.diag(raw_diagonal) - np.outer(mean, mean)
    return mean, covariance


def _assert_realizable(_mean, covariance, tolerance=1.0e-14):
    np.testing.assert_allclose(covariance, covariance.T, atol=tolerance)
    assert np.min(np.linalg.eigvalsh(covariance)) >= -tolerance


def test_restriction_does_not_convert_drift_variation_into_covariance():
    """Average central covariance, not raw moments, across fine children."""
    fine_means = np.array(
        [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
        ]
    )
    fine_covariances = np.zeros((2, 3, 3))

    restricted_mean = np.mean(fine_means, axis=0)
    restricted_covariance = np.mean(fine_covariances, axis=0)
    np.testing.assert_array_equal(restricted_mean, np.zeros(3))
    np.testing.assert_array_equal(restricted_covariance, np.zeros((3, 3)))

    raw_second = fine_covariances + np.einsum(
        "ni,nj->nij", fine_means, fine_means
    )
    old_covariance = np.mean(raw_second, axis=0) - np.outer(
        restricted_mean, restricted_mean
    )
    assert old_covariance[0, 0] == 1.0


def test_constant_prolongation_and_central_covariance_averaging_are_realizable():
    """Convex averages of central covariance remain positive semidefinite."""
    rng = np.random.default_rng(918273)
    for _ in range(1000):
        probabilities = rng.dirichlet(np.ones(7), size=8)[:, :6]
        coefficients = [
            _central_coefficients(cell_probabilities)
            for cell_probabilities in probabilities
        ]
        means = np.array([coefficient[0] for coefficient in coefficients])
        covariances = np.array(
            [coefficient[1] for coefficient in coefficients]
        )

        # Restriction and CIC are positive weighted averages of central coefficients.
        weights = rng.random(means.shape[0])
        weights /= np.sum(weights)
        restricted_mean = weights @ means
        restricted_covariance = np.tensordot(weights, covariances, axes=1)
        _assert_realizable(restricted_mean, restricted_covariance)

        cic_weights = rng.random(means.shape[0])
        cic_weights /= np.sum(cic_weights)
        interpolated_mean = cic_weights @ means
        interpolated_covariance = np.tensordot(
            cic_weights, covariances, axes=1
        )
        _assert_realizable(interpolated_mean, interpolated_covariance)

        # Piecewise-constant coarse-to-fine prolongation copies a realizable tuple.
        for _child in range(2):
            np.testing.assert_array_equal(
                restricted_covariance, restricted_covariance.copy()
            )
            _assert_realizable(restricted_mean, restricted_covariance)


def test_ito_uses_particle_specific_amr_prolongation():
    """Ito central coefficients must use piecewise-constant prolongation."""
    source = ITO_SOURCE.read_text(encoding="utf-8")
    start = source.index("TaskStatus Particles::ProlongateItoCoefficients")
    end = source.index("//! \\fn TaskStatus Particles::PushIto2", start)
    body = source[start:end]

    assert "pbval_ito->ProlongateCC" not in body
    assert "ito2_prolong_central_coefficients" in body


def _assert_conditional_moments(displacement, expected_mean, expected_covariance):
    sample_mean = np.mean(displacement, axis=0)
    sample_covariance = np.cov(displacement, rowvar=False, bias=True)
    count = displacement.shape[0]
    for axis in range(2):
        standard_error = np.sqrt(expected_covariance[axis, axis] / count)
        assert abs(sample_mean[axis] - expected_mean[axis]) <= 5.0 * standard_error
    for row in range(2):
        for column in range(row, 2):
            products = (
                (displacement[:, row] - expected_mean[row])
                * (displacement[:, column] - expected_mean[column])
            )
            standard_error = np.std(products, ddof=1) / np.sqrt(count)
            assert (
                abs(sample_covariance[row, column] - expected_covariance[row, column])
                <= 5.0 * standard_error
            )


def test_static_amr_coarse_and_fine_conditional_moments():
    """Actual coarse/fine regions preserve the chosen coefficient-field moments."""
    shutil.rmtree(RUN_DIR, ignore_errors=True)
    try:
        assert testutils.run(INPUT, ["-d", str(RUN_DIR)])
        data = read_history(
            RUN_DIR / "prtcl_thermo_history/ito_amr_moments.prtcl_thermo_history.thp"
        )
        cycles = np.unique(data["cycle"])
        assert cycles.size == 2
        initial = data["cycle"] == cycles[0]
        final = data["cycle"] == cycles[1]
        initial_order = np.argsort(data["tag"][initial])
        final_tags = data["tag"][final]
        _, final_unique = np.unique(final_tags, return_index=True)
        final_indices = np.flatnonzero(final)[final_unique]
        final_order = final_indices[np.argsort(data["tag"][final_indices])]
        np.testing.assert_array_equal(
            data["tag"][initial][initial_order], data["tag"][final_order]
        )
        seed_id = data["seed_id"][initial][initial_order]
        start = np.column_stack(
            (data["x1"][initial][initial_order], data["x2"][initial][initial_order])
        )
        end = np.column_stack(
            (data["x1"][final_order], data["x2"][final_order])
        )
        displacement = end - start
        displacement -= np.rint(displacement)

        dt = float(np.max(data["time"][final]))
        velocity = np.array([0.4, 0.2])
        expected_mean = velocity * dt
        for current_seed, cell_width in ((1, 1.0 / 32.0), (2, 1.0 / 64.0)):
            selected = seed_id == current_seed
            assert np.count_nonzero(selected) == 65536
            raw_diagonal = cell_width * np.abs(velocity) * dt
            expected_covariance = np.diag(raw_diagonal) - np.outer(
                expected_mean, expected_mean
            )
            _assert_conditional_moments(
                displacement[selected], expected_mean, expected_covariance
            )
    finally:
        shutil.rmtree(RUN_DIR, ignore_errors=True)
