#!/usr/bin/env python3
"""Focused tests for bounded Q-006 Section 5.3 source-local preparation."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest

from tst.publication import analyze_q006_paper_multispecies_oscillation as osc


REPO_ROOT = Path(__file__).resolve().parents[2]
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q006_paper_multispecies_oscillation_source_local_preparation_2026-05-30.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q006PaperMultispeciesOscillationTests(unittest.TestCase):
    def test_paper_formula_and_code_normalization(self) -> None:
        contract = osc.analytical_contract()
        self.assertEqual(contract["aggregate_ppc"], 128.0)
        self.assertEqual(contract["species_ppc_each"], 64.0)
        self.assertEqual(contract["deposit_qscale"], 0.0234375)
        self.assertEqual(contract["gyro_omega"], 1.0)
        self.assertEqual(contract["oscillation_omega"], 2.0)
        self.assertAlmostEqual(contract["oscillation_frequency_cycles"], 1.0 / math.pi)
        self.assertAlmostEqual(contract["oscillation_period"], math.pi)
        self.assertAlmostEqual(contract["gas_uy"], -0.3)
        self.assertAlmostEqual(contract["initial_kinetic_energy_density"], 0.06)
        self.assertAlmostEqual(contract["ideal_compatibility_pressure"], 0.6)
        self.assertAlmostEqual(
            64.0 * contract["deposit_qscale"] * 1.0,
            contract["species_mass_density_each"],
        )
        self.assertAlmostEqual(
            1.0 * contract["gas_uy"] + 2.0 * 1.5 * 0.1,
            0.0,
        )

    def test_three_decks_freeze_uniform_smr_and_audited_amr_preparation(self) -> None:
        decks = osc.validate_source_local_candidate_decks()
        self.assertEqual(
            [item["grid_setup"] for item in decks],
            ["uniform", "smr", "audited_amr_preparation"],
        )
        self.assertTrue(all(item["cycle_zero_only"] for item in decks))
        self.assertTrue(all(not item["qualifying_evidence"] for item in decks))
        self.assertEqual(decks[1]["smr_refined_volume_fraction"], 0.125)
        self.assertIn("10_refine_60_derefine", decks[2]["amr_policy"])

    def test_audited_amr_policy_is_deterministic_and_bounded(self) -> None:
        first = [
            osc.audited_amr_decision(osc.AMR_SEED, cycle, gid)
            for cycle in range(8)
            for gid in range(32)
        ]
        second = [
            osc.audited_amr_decision(osc.AMR_SEED, cycle, gid)
            for cycle in range(8)
            for gid in range(32)
        ]
        self.assertEqual(first, second)
        self.assertEqual(set(first), {-1, 0, 1})
        summary = osc.audited_amr_policy_summary()
        self.assertAlmostEqual(summary["refine_fraction"], 0.10, delta=0.01)
        self.assertAlmostEqual(summary["derefine_fraction"], 0.60, delta=0.02)
        self.assertAlmostEqual(summary["retain_fraction"], 0.30, delta=0.02)
        self.assertFalse(summary["runtime_exercised"])
        self.assertFalse(summary["true_amr_policy_qualification"])

    def test_candidate_deck_drift_fails_closed(self) -> None:
        mutations = (
            ("uniform", "ppc                                = 128.0",
             "ppc                                = 64.0", "particles/ppc"),
            ("smr", "x1max = 12.0", "x1max = 8.0", "one-eighth region"),
            ("audited_amr_preparation", "amr_seed                     = 60053",
             "amr_seed                     = 1", "AMR seed"),
            ("uniform", "eos         = ideal", "eos         = isothermal", "mhd/eos"),
        )
        for grid_setup, old, new, error in mutations:
            with self.subTest(grid_setup=grid_setup, error=error):
                with tempfile.TemporaryDirectory() as temporary_directory:
                    path = Path(temporary_directory) / "candidate.athinput"
                    original = osc.DECKS[grid_setup].read_text(encoding="utf-8")
                    self.assertIn(old, original)
                    path.write_text(original.replace(old, new, 1), encoding="utf-8")
                    with self.assertRaisesRegex(osc.ContractError, error):
                        osc.validate_candidate_deck(path, grid_setup)

    def test_dedicated_source_contract_and_parent_integrated_registration(self) -> None:
        source = osc.validate_source_contract()
        self.assertEqual(source["pgen_name"], osc.PGEN_NAME)
        self.assertEqual(source["method"], osc.PGEN_METHOD)
        self.assertTrue(source["dedicated_source_contract"])
        self.assertTrue(source["audited_amr_callback_prepared"])
        self.assertTrue(source["shared_registration"]["complete"])
        self.assertTrue(source["shared_registration"]["compilation_unit_registered"])
        self.assertTrue(source["shared_registration"]["method_declared"])
        self.assertTrue(source["shared_registration"]["fresh_dispatch_registered"])
        self.assertTrue(source["shared_registration"]["restart_dispatch_registered"])

    def test_report_separates_freeze_from_open_qualification_work(self) -> None:
        report = osc.build_preparation_report()
        self.assertTrue(report["deck_source_freeze"])
        self.assertFalse(report["qualifying_evidence"])
        self.assertFalse(report["exact_section53_isothermal_runtime_supported"])
        self.assertFalse(report["long_horizon_runtime_evidence"])
        self.assertFalse(report["true_amr_policy_qualification"])
        self.assertFalse(report["mpi_qualification"])
        self.assertFalse(report["gpu_qualification"])
        self.assertFalse(report["frontier_authorization"])
        self.assertFalse(report["external_review"])
        self.assertIn("preparation_only", report["qualification_effect"])
        self.assertEqual(report["parent_registration_steps"], [])

    def test_sidecar_binds_only_new_q006_artifacts_and_nonqualification(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-006")
        self.assertEqual(sidecar["campaign_id"], osc.CAMPAIGN_ID)
        self.assertEqual(sidecar["qualification_effect"], osc.QUALIFICATION_EFFECT)
        self.assertFalse(sidecar["claim_closure"])
        self.assertEqual(sidecar["frontier_authorization"], "not_bound")
        expected_paths = {
            "src/pgen/tests/q006_paper_multispecies_oscillation.cpp",
            "inputs/tests/pic_q006_paper_multispecies_oscillation_uniform_candidate.athinput",
            "inputs/tests/pic_q006_paper_multispecies_oscillation_smr_candidate.athinput",
            "inputs/tests/pic_q006_paper_multispecies_oscillation_audited_amr_candidate.athinput",
            "tst/publication/analyze_q006_paper_multispecies_oscillation.py",
            "tst/publication/test_analyze_q006_paper_multispecies_oscillation.py",
        }
        bindings = sidecar["artifact_bindings"]
        self.assertEqual(set(bindings), expected_paths)
        for relative, expected_sha256 in bindings.items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected_sha256)
        boundary = sidecar["nonqualification_boundary"]
        self.assertTrue(boundary["deck_source_freeze"])
        self.assertFalse(boundary["qualifying_evidence"])
        self.assertFalse(boundary["long_horizon_runtime_evidence"])
        self.assertFalse(boundary["true_amr_policy_qualification"])
        self.assertFalse(boundary["mpi_qualification"])
        self.assertFalse(boundary["gpu_qualification"])
        self.assertFalse(boundary["external_review"])
        self.assertEqual(sidecar["parent_registration_steps"], [])


if __name__ == "__main__":
    unittest.main()
