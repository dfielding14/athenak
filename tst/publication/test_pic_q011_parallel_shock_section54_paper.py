#!/usr/bin/env python3
"""Focused tests for the Q-011 Section 5.4 preparation tranche."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tst"))

from scripts.particles import pic_parallel_shock_section54_paper as q011  # noqa: E402

SIDECAR = (
    REPO_ROOT
    / "tst"
    / "publication"
    / "readiness"
    / "q011_parallel_shock_section54_paper_preparation_2026-05-30.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PicQ011ParallelShockSection54PaperTests(unittest.TestCase):
    def test_deck_contract_freezes_paper_values_and_open_items(self) -> None:
        contract = q011.build_preparation_contract()
        self.assertEqual(contract["gate"], "Q-011")
        self.assertEqual(
            contract["qualification_effect"],
            "source_controlled_preparation_only",
        )
        values = contract["paper_text_values"]
        self.assertEqual(values["domain_c_over_omega_pi"], [48000.0, 3120.0])
        self.assertEqual(values["amr_cell_sizes_c_over_omega_pi"], [12.0, 6.0, 3.0])
        self.assertEqual(values["snapshot_times_omega0_inverse"], [500.0, 1200.0])
        self.assertEqual(values["exclude_birth_time_before_omega0_inverse"], 45.0)
        self.assertEqual(values["shock_surface_model"], "ideal_surface")
        self.assertEqual(
            values["injection_distribution"],
            "monoenergetic_full_sphere_isotropic_relative_to_ideal_surface",
        )
        derivation = contract["deck"]["shock_surface_derivation"]
        self.assertEqual(derivation["selected_model"], "ideal_surface")
        self.assertAlmostEqual(derivation["ideal_surface_speed_over_ua0"], 10.0)
        self.assertAlmostEqual(
            derivation["finite_mach_engineering_option_speed_over_ua0"],
            10.007408779453616,
        )
        self.assertAlmostEqual(
            derivation["upstream_relative_sweep_speed_over_ua0"],
            40.00000000005,
        )
        self.assertEqual(
            derivation["shock_surface_carrier_selection"],
            "single_half_open_cell_with_surface_x1",
        )
        self.assertEqual(
            derivation["early_injected_particle_removal"],
            "runtime_state_removal_for_birth_time_below_45",
        )
        self.assertAlmostEqual(
            derivation["ideal_surface_positions_c_over_omega_pi"]["t500"],
            5000.000000025,
        )
        self.assertAlmostEqual(
            derivation["ideal_surface_positions_c_over_omega_pi"]["t1200"],
            12000.00000006,
        )
        self.assertNotIn(
            "finite_mach_shock_speed_estimate_to_ideal_surface_mapping_audit",
            contract["open_items"],
        )
        self.assertIn(
            "executed_shock_surface_injection_distribution_audit",
            contract["open_items"],
        )
        self.assertTrue(contract["open_items"])
        self.assertEqual(contract["result_metrics"], [])

    def test_static_regression_passes_without_campaign_results(self) -> None:
        q011.run()
        self.assertTrue(q011.analyze())

    def test_campaign_artifact_analysis_fails_closed(self) -> None:
        with self.assertRaisesRegex(q011.ContractError, "campaign analysis is blocked"):
            q011.analyze_campaign_artifacts({})

    def test_sidecar_freezes_analyzer_manifest_contract(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        frozen = sidecar["analyzer_contract"]
        contract = q011.build_preparation_contract()["artifact_manifest_contract"]
        for field in (
            "required_candidate_bindings",
            "particle_filters",
            "primary_observables",
            "required_snapshot_times_omega0_inverse",
            "required_grid_variants",
            "required_raw_artifacts_per_snapshot",
            "required_run_artifacts",
            "required_seeds",
        ):
            with self.subTest(field=field):
                self.assertEqual(frozen[field], contract[field])

    def test_sidecar_hashes_and_orion_only_boundary(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-011")
        self.assertEqual(sidecar["qualification_effect"], "preparation_only")
        policy = sidecar["execution_policy"]
        self.assertEqual(policy["frontier_campaign_execution"], "open_not_run")
        self.assertEqual(policy["allowed_bulk_storage_systems"], ["Orion"])
        self.assertEqual(policy["forbidden_storage_systems"], ["Kronos"])
        for artifact in sidecar["artifacts"]:
            path = REPO_ROOT / artifact["path"]
            with self.subTest(path=path):
                self.assertTrue(path.is_file())
                self.assertEqual(artifact["sha256"], _sha256(path))


if __name__ == "__main__":
    unittest.main()
