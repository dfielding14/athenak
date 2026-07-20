#!/usr/bin/env python3
"""Focused tests for the source-local nonlinear Bell campaign foundation."""

from __future__ import annotations

import json
import math
from pathlib import Path
import unittest

import numpy as np

from tst.publication import analyze_q023_prod_bell_nonlinear_foundation as bell


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q023_prod_bell_nonlinear_source_local_foundation_2026-06-06.json"
)


def _synthetic_datasets() -> list[dict[str, object]]:
    nx1 = 64
    nx2 = 64
    l1, l2, _ = bell.EXPECTED_EXTENT
    x1 = (np.arange(nx1, dtype=float) + 0.5) * l1 / nx1
    x2 = (np.arange(nx2, dtype=float) + 0.5) * l2 / nx2
    x2_grid, x1_grid = np.meshgrid(x2, x1, indexing="ij")
    parallel, transverse_a, transverse_b = bell._mode_basis()
    phase = bell.K0 * (parallel[0] * x1_grid + parallel[1] * x2_grid)
    low_k_phase = 0.5 * phase
    datasets = []
    for tau in np.linspace(0.0, 25.0, 51):
        if tau <= 8.0:
            fast_amplitude = 1.0e-4 * math.exp(bell.EXPECTED_LINEAR_GROWTH * tau)
            low_k_amplitude = 0.0
        else:
            fast_amplitude = 0.145 * (1.0 + 0.01 * math.sin(tau))
            low_k_amplitude = 0.06 * (1.0 - math.exp(-(tau - 8.0) / 3.0))
        gas_velocity_parallel = (
            0.0
            if tau <= 8.0
            else 1.4 * (1.0 - math.exp(-(tau - 8.0) / 4.0))
        )
        mean_density = (
            1.0
            if tau <= 8.0
            else 1.0 + 0.02 * (1.0 - math.exp(-(tau - 8.0) / 5.0))
        )
        fast = (
            fast_amplitude * np.cos(phase)[None, ...] * transverse_a[:, None, None]
            - fast_amplitude * np.sin(phase)[None, ...] * transverse_b[:, None, None]
        )
        low_k = (
            low_k_amplitude
            * np.cos(low_k_phase)[None, ...]
            * transverse_a[:, None, None]
            - low_k_amplitude
            * np.sin(low_k_phase)[None, ...]
            * transverse_b[:, None, None]
        )
        magnetic = parallel[:, None, None] + fast + low_k
        velocity = (
            fast_amplitude
            * (
                (
                    -bell.EPSILON * np.cos(phase)
                    + bell.EXPECTED_LINEAR_GROWTH * np.sin(phase)
                )[None, ...]
                * transverse_a[:, None, None]
                + (
                    bell.EXPECTED_LINEAR_GROWTH * np.cos(phase)
                    + bell.EPSILON * np.sin(phase)
                )[None, ...]
                * transverse_b[:, None, None]
            )
            + 0.5 * low_k
            + gas_velocity_parallel * parallel[:, None, None]
        )
        datasets.append(
            {
                "Time": tau / (bell.K0 * bell.U_A),
                "x1v": x1,
                "x2v": x2,
                "x3v": np.asarray([0.5]),
                "dens": mean_density * np.ones((1, nx2, nx1)),
                "bcc1": magnetic[0][None, ...],
                "bcc2": magnetic[1][None, ...],
                "bcc3": magnetic[2][None, ...],
                "velx": velocity[0][None, ...],
                "vely": velocity[1][None, ...],
                "velz": velocity[2][None, ...],
            }
        )
    return datasets


class Q023ProdBellNonlinearFoundationTests(unittest.TestCase):
    def test_production_deck_freezes_vl2_contract_and_no_authority(self) -> None:
        contract = bell.validate_production_deck()
        self.assertEqual(contract["pressure_p0"], 1.0)
        self.assertEqual(contract["epsilon_ua_over_vcr"], 0.4)
        self.assertEqual(contract["carrier_periods_by_periodic_axis"], [8.0, 8.0])
        self.assertEqual(contract["cells_per_carrier_period_by_axis"], [64, 32])
        self.assertTrue(contract["mhd_spatial_spectra_supported"])
        self.assertFalse(contract["particle_spectra_supported"])
        self.assertFalse(contract["launch_authorized"])
        self.assertFalse(contract["scientific_claim_authorized"])

    def test_linear_to_nonlinear_fixture_reports_bounded_diagnostics(self) -> None:
        report = bell.analyze_datasets(
            _synthetic_datasets(), source_kind="synthetic_contract_fixture"
        )
        self.assertTrue(report["foundation_diagnostics_complete"])
        self.assertFalse(report["passed"])
        self.assertFalse(report["launch_authorized"])
        self.assertFalse(report["scientific_claim_authorized"])
        self.assertFalse(report["scientific_scope"]["particle_spectra_available"])
        self.assertFalse(
            report["linear_growth_diagnostic"]["science_acceptance_threshold_frozen"]
        )
        self.assertAlmostEqual(
            report["linear_growth_diagnostic"]["measured_growth_rate_over_k0_ua"],
            bell.EXPECTED_LINEAR_GROWTH,
            places=11,
        )
        self.assertTrue(report["nonlinear_onset_diagnostic"]["detected"])
        self.assertEqual(
            report["scientific_scope"]["anisotropic_pressure_saturation_relation"],
            "not_mapped_not_claimed_normalization_and_CR_pressure_tensor_open",
        )
        self.assertFalse(
            report["scientific_scope"]["evolving_cr_current_or_speed_measured"]
        )
        self.assertGreater(
            report["saturation_window_diagnostic"]["median_delta_b_rms_over_bg"],
            bell.NONLINEAR_ONSET_DELTA_B_OVER_BG,
        )
        final = report["trace"][-1]
        self.assertAlmostEqual(
            final["fixed_cr_stream_parallel_speed_over_ua"],
            bell.FIXED_CR_STREAM_PARALLEL_SPEED / bell.U_A,
        )
        self.assertAlmostEqual(
            final[
                "configured_fixed_cr_minus_gas_relative_drift_parallel_over_ua"
            ],
            final["fixed_cr_stream_parallel_speed_over_ua"]
            - final["density_weighted_mean_gas_velocity_parallel_over_ua"],
        )
        self.assertAlmostEqual(
            final["amplified_total_field_alfven_speed_over_ua"],
            final["b_rms_over_bg"]
            / math.sqrt(final["mean_density_over_rho0"]),
        )
        self.assertAlmostEqual(
            final[
                "v_a_amplified_over_abs_configured_fixed_cr_minus_gas_"
                "relative_drift"
            ],
            final["amplified_total_field_alfven_speed_over_ua"]
            / final[
                "abs_configured_fixed_cr_minus_gas_relative_drift_parallel_"
                "over_ua"
            ],
        )
        mechanism = report["saturation_window_diagnostic"][
            "constant_current_mechanism_diagnostic"
        ]
        self.assertFalse(mechanism["science_acceptance_threshold_frozen"])
        summaries = mechanism["fixed_window_summaries"]
        self.assertGreater(
            summaries[
                "density_weighted_mean_gas_velocity_parallel_over_ua"
            ]["end"],
            summaries[
                "density_weighted_mean_gas_velocity_parallel_over_ua"
            ]["start"],
        )
        self.assertLess(
            summaries[
                "abs_configured_fixed_cr_minus_gas_relative_drift_parallel_"
                "over_ua"
            ]["end"],
            summaries[
                "abs_configured_fixed_cr_minus_gas_relative_drift_parallel_"
                "over_ua"
            ]["start"],
        )
        self.assertGreater(
            summaries[
                "v_a_amplified_over_abs_configured_fixed_cr_minus_gas_"
                "relative_drift"
            ]["end"],
            summaries[
                "v_a_amplified_over_abs_configured_fixed_cr_minus_gas_"
                "relative_drift"
            ]["start"],
        )
        spectra = {
            row["label"]: row for row in report["selected_magnetic_spatial_spectra"]
        }
        self.assertAlmostEqual(spectra["initial"]["dominant_k_over_k0"], 1.0)
        self.assertAlmostEqual(
            spectra["initial"]["dominant_wavelength_over_seed_wavelength"], 1.0
        )
        self.assertGreater(
            spectra["final"]["low_k_power_fraction_k_over_k0_lt_0p75"],
            spectra["initial"]["low_k_power_fraction_k_over_k0_lt_0p75"] + 0.05,
        )

    def test_parallel_gas_velocity_is_density_weighted(self) -> None:
        dataset = _synthetic_datasets()[-1]
        density = np.asarray(dataset["dens"])
        density[:, :, : density.shape[-1] // 2] *= 3.0
        dataset["dens"] = density
        velocity = np.stack(
            [np.asarray(dataset[name]) for name in ("velx", "vely", "velz")]
        )
        parallel, _, _ = bell._mode_basis()
        velocity[:, :, :, : velocity.shape[-1] // 2] += (
            0.5 * parallel[:, None, None, None]
        )
        for name, component in zip(("velx", "vely", "velz"), velocity):
            dataset[name] = component
        parallel_velocity = np.tensordot(parallel, velocity, axes=1)
        expected_density_weighted = float(
            np.sum(density * parallel_velocity) / np.sum(density)
        )
        unweighted = float(np.mean(parallel_velocity))
        metrics = bell._snapshot_metrics(dataset, bell._geometry(dataset))
        self.assertAlmostEqual(
            metrics["density_weighted_mean_gas_velocity_parallel_over_ua"],
            expected_density_weighted / bell.U_A,
        )
        self.assertNotAlmostEqual(expected_density_weighted, unweighted)
        self.assertAlmostEqual(
            metrics["mean_density_over_rho0"],
            float(np.mean(density)) / bell.RHO0,
        )

    def test_fixed_windows_fail_closed_when_run_is_incomplete(self) -> None:
        with self.assertRaisesRegex(
            bell.ContractError, "saturation-diagnostic window"
        ):
            bell.analyze_datasets(
                _synthetic_datasets()[:30], source_kind="synthetic_contract_fixture"
            )

    def test_duplicate_time_nonfinite_field_and_shape_drift_fail_closed(self) -> None:
        datasets = _synthetic_datasets()
        datasets[1]["Time"] = datasets[0]["Time"]
        with self.assertRaisesRegex(bell.ContractError, "unique and increasing"):
            bell.analyze_datasets(
                datasets, source_kind="synthetic_contract_fixture"
            )

        datasets = _synthetic_datasets()
        datasets[0]["bcc3"][0, 0, 0] = np.nan
        with self.assertRaisesRegex(bell.ContractError, "bcc3 shape or values"):
            bell.analyze_datasets(
                datasets, source_kind="synthetic_contract_fixture"
            )

        datasets = _synthetic_datasets()
        datasets[-1]["velz"] = datasets[-1]["velz"][:, :, :-1]
        with self.assertRaisesRegex(bell.ContractError, "velz shape or values"):
            bell.analyze_datasets(
                datasets, source_kind="synthetic_contract_fixture"
            )

    def test_raw_provenance_and_snapshot_count_fail_closed(self) -> None:
        with self.assertRaisesRegex(bell.ContractError, "artifact binding"):
            bell.analyze_datasets(
                _synthetic_datasets(), source_kind="source_local_raw_mhd_w_bcc"
            )
        with self.assertRaisesRegex(bell.ContractError, "bounded snapshot count"):
            bell.analyze_datasets(
                [{}] * (bell.MAX_SNAPSHOTS + 1),
                source_kind="synthetic_contract_fixture",
            )

    def test_source_local_readiness_binds_only_new_non_authorizing_artifacts(
        self,
    ) -> None:
        record = json.loads(READINESS.read_text(encoding="utf-8"))
        self.assertEqual(record["campaign_id"], bell.CAMPAIGN_ID)
        self.assertEqual(
            record["qualification_effect"],
            "none_source_local_foundation_does_not_authorize_launch_or_science_claims",
        )
        self.assertFalse(record["authority"]["launch_authorized"])
        self.assertFalse(record["authority"]["scientific_claim_authorized"])
        self.assertFalse(record["authority"]["publication_authorized"])
        self.assertIn(
            "v_a_amplified_over_abs_configured_fixed_cr_minus_gas_relative_drift",
            record["foundation_scope"]["supported_diagnostics"],
        )
        self.assertEqual(
            record["literature_relationship"][
                "anisotropic_pressure_saturation_relation"
            ],
            "not_mapped_not_claimed_because_deck_normalization_and_CR_"
            "pressure_tensor_are_not_bound",
        )
        self.assertTrue(
            all(
                "Zacharegkas" not in item
                for item in record["literature_relationship"]["nonlinear_context"]
            )
        )
        bindings = {
            row["path"]: row["sha256"] for row in record["artifact_bindings"]
        }
        expected_paths = {
            str(bell.DECK.relative_to(REPO_ROOT)),
            "tst/publication/analyze_q023_prod_bell_nonlinear_foundation.py",
            "tst/publication/test_analyze_q023_prod_bell_nonlinear_foundation.py",
        }
        self.assertEqual(set(bindings), expected_paths)
        self.assertIn(
            "Q022 reference-specific nonlinear Bell mapping and numeric tolerances",
            record["open_dependencies"],
        )


if __name__ == "__main__":
    unittest.main()
