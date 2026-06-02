#!/usr/bin/env python3
"""Synthetic tests for the pure Q-011 Section 5.4 particle reducers."""

from __future__ import annotations

import json
import unittest

import numpy as np

try:
    from tst.publication import q011_section54_particles as particles
except ModuleNotFoundError:
    import q011_section54_particles as particles


def _pvtk_velocity_from_chi(chi: np.ndarray) -> np.ndarray:
    momentum_per_mass = particles.UPSTREAM_SPEED_U0 * np.sqrt(chi)
    velocity = momentum_per_mass / np.sqrt(
        1.0 + (momentum_per_mass / particles.PARTICLE_LIGHT_SPEED) ** 2
    )
    return np.column_stack((velocity, np.zeros_like(velocity), np.zeros_like(velocity)))


def _tail_samples(slope: float = -1.5) -> tuple[np.ndarray, np.ndarray]:
    edges = np.asarray(particles.CHI_BIN_EDGES)
    centers = np.sqrt(edges[:-1] * edges[1:])
    selected = (centers >= 20.0) & (centers <= 160.0)
    chi = centers[selected]
    macro_weight = 7.0 * chi**slope * np.diff(edges)[selected]
    return chi, macro_weight


class Q011Section54ParticleReducerTests(unittest.TestCase):
    def test_downstream_filter_archives_disjoint_rejection_census(self) -> None:
        filtered = particles.downstream_filter(
            snapshot_time=2.0,
            x1=[1.0, 1.0, 21.0, 20.0, 19.0],
            cr_source=[0, 1, 1, 1, 1],
            birth_time=[50.0, 44.0, 45.0, 45.0, 45.0],
            macro_weight=[1.0, 2.0, 3.0, 4.0, 5.0],
        )
        np.testing.assert_array_equal(filtered.mask, [False, False, False, False, True])
        self.assertEqual(filtered.record["ideal_surface_x1_c_over_omega_pi"], 20.0)
        census = filtered.record["disjoint_census"]
        self.assertEqual(census["all_particles"], {"particle_count": 5, "macro_weight": 15.0})
        self.assertEqual(
            census["admitted_downstream"], {"particle_count": 1, "macro_weight": 5.0}
        )
        self.assertEqual(
            census["rejected_total"], {"particle_count": 4, "macro_weight": 10.0}
        )
        self.assertEqual(census["rejected_wrong_source"]["particle_count"], 1)
        self.assertEqual(census["rejected_early_birth_time"]["particle_count"], 1)
        self.assertEqual(census["rejected_upstream"]["particle_count"], 1)
        self.assertEqual(census["rejected_on_surface"]["particle_count"], 1)

    def test_chi_reconstruction_inverts_pvtk_physical_velocity_output(self) -> None:
        expected = np.asarray([1.0, 10.0, 1024.0])
        reconstructed = particles.reconstruct_chi_from_pvtk_velocity(
            _pvtk_velocity_from_chi(expected)
        )
        np.testing.assert_allclose(reconstructed, expected, rtol=2.0e-13)

    def test_overflow_contributes_to_normalization_and_overflow_fraction(self) -> None:
        spectrum = particles.weighted_spectrum_record(
            [2.0, 2048.0],
            [1.0, 3.0],
        )
        self.assertEqual(spectrum["macro_weight_in_bins"], 1.0)
        self.assertEqual(spectrum["overflow_macro_weight"], 3.0)
        self.assertEqual(spectrum["total_post_filter_macro_weight"], 4.0)
        self.assertEqual(spectrum["overflow_macro_weight_fraction"], 0.75)
        self.assertFalse(spectrum["overflow_gate"]["passed"])

        nonzero = np.flatnonzero(spectrum["weighted_counts"])
        self.assertEqual(nonzero.size, 1)
        index = int(nonzero[0])
        edges = np.asarray(spectrum["bin_edges"])
        center = np.sqrt(edges[index] * edges[index + 1])
        width = edges[index + 1] - edges[index]
        self.assertAlmostEqual(
            spectrum["normalized_chi_f_chi"][index],
            center * (1.0 / width) / 4.0,
        )
        self.assertNotAlmostEqual(
            spectrum["normalized_chi_f_chi"][index],
            center * (1.0 / width) / 1.0,
        )

    def test_underflow_and_overflow_both_contribute_to_total_post_filter_weight(self) -> None:
        spectrum = particles.weighted_spectrum_record(
            [0.5, 2.0, 2048.0],
            [2.0, 1.0, 3.0],
        )
        self.assertEqual(spectrum["underflow_macro_weight"], 2.0)
        self.assertEqual(spectrum["macro_weight_in_bins"], 1.0)
        self.assertEqual(spectrum["overflow_macro_weight"], 3.0)
        self.assertEqual(spectrum["total_post_filter_macro_weight"], 6.0)
        self.assertEqual(spectrum["overflow_macro_weight_fraction"], 0.5)

    def test_late_slope_requires_eight_positive_bins_and_applies_gate(self) -> None:
        chi, macro_weight = _tail_samples()
        spectrum = particles.weighted_spectrum_record(chi, macro_weight)
        slope = particles.late_slope_record(spectrum["f_chi"])
        self.assertGreaterEqual(slope["positive_fit_bin_count"], 8)
        self.assertAlmostEqual(slope["slope"], -1.5, places=12)
        self.assertTrue(slope["slope_gate_passed"])

        off_target_chi, off_target_weight = _tail_samples(slope=-1.0)
        off_target = particles.late_slope_record(
            particles.weighted_spectrum_record(off_target_chi, off_target_weight)["f_chi"]
        )
        self.assertAlmostEqual(off_target["slope"], -1.0, places=12)
        self.assertFalse(off_target["slope_gate_passed"])

        with self.assertRaisesRegex(
            particles.ParticleReducerError,
            "insufficient positive fit bins",
        ):
            particles.late_slope_record(
                particles.weighted_spectrum_record(chi[:7], macro_weight[:7])["f_chi"]
            )

    def test_snapshot_reduction_returns_deterministic_archive_record(self) -> None:
        chi, macro_weight = _tail_samples()
        count = chi.size
        record = particles.reduce_particle_snapshot(
            snapshot_time=1200.0,
            points=np.column_stack(
                (np.full(count, 100.0), np.zeros(count), np.zeros(count))
            ),
            cr_source=np.ones(count, dtype=np.int64),
            birth_time=np.full(count, 45.0),
            velocity=_pvtk_velocity_from_chi(chi),
            macro_weight=macro_weight,
            evaluate_late_slope=True,
        )
        census = record["particle_filter"]["disjoint_census"]
        self.assertEqual(census["admitted_downstream"]["particle_count"], count)
        self.assertEqual(census["rejected_total"]["particle_count"], 0)
        self.assertTrue(record["late_slope"]["slope_gate_passed"])
        first = particles.canonical_record_bytes(record)
        second = particles.canonical_record_bytes(record)
        self.assertEqual(first, second)
        self.assertEqual(json.loads(first), record)

    def test_reducers_reject_nonfinite_negative_and_shape_drift(self) -> None:
        filter_arguments = {
            "snapshot_time": 2.0,
            "x1": [1.0],
            "cr_source": [1],
            "birth_time": [45.0],
            "macro_weight": [1.0],
        }
        with self.assertRaisesRegex(particles.ParticleReducerError, "x1 must be finite"):
            particles.downstream_filter(**{**filter_arguments, "x1": [np.nan]})
        with self.assertRaisesRegex(
            particles.ParticleReducerError, "macro_weight must be non-negative"
        ):
            particles.downstream_filter(**{**filter_arguments, "macro_weight": [-1.0]})
        with self.assertRaisesRegex(particles.ParticleReducerError, "shape drifted"):
            particles.downstream_filter(**{**filter_arguments, "birth_time": [[45.0]]})
        with self.assertRaisesRegex(
            particles.ParticleReducerError, "decoded integers"
        ):
            particles.downstream_filter(**{**filter_arguments, "cr_source": [1.0]})
        with self.assertRaisesRegex(particles.ParticleReducerError, "velocity must be finite"):
            particles.reconstruct_chi_from_pvtk_velocity([[np.inf, 0.0, 0.0]])
        with self.assertRaisesRegex(
            particles.ParticleReducerError, "below the frozen particle light speed"
        ):
            particles.reconstruct_chi_from_pvtk_velocity(
                [[particles.PARTICLE_LIGHT_SPEED, 0.0, 0.0]]
            )
        with self.assertRaisesRegex(particles.ParticleReducerError, "points shape drifted"):
            particles.reduce_particle_snapshot(
                snapshot_time=1200.0,
                points=[[1.0, 0.0]],
                cr_source=[1],
                birth_time=[45.0],
                velocity=[[1.0, 0.0, 0.0]],
                macro_weight=[1.0],
            )


if __name__ == "__main__":
    unittest.main()
