#!/usr/bin/env python3
"""Focused tests for bounded Q-007 true-delta-f source-local preparation."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest

from tst.publication import analyze_q007_paper_deltaf_linear_preparation as q007


REPO_ROOT = Path(__file__).resolve().parents[2]
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q007_paper_deltaf_linear_source_local_preparation_2026-05-30.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q007PaperDeltaFLinearPreparationTests(unittest.TestCase):
    def test_crsi_static_paper_mapping(self) -> None:
        contract = q007.analytical_contract()
        crsi = contract["crsi"]
        self.assertAlmostEqual(crsi["k0"], 1.0 / 300.0)
        self.assertAlmostEqual(crsi["lambda0"], 600.0 * math.pi)
        self.assertGreater(crsi["q2_at_k0"], 0.0)
        self.assertGreater(crsi["forward_growth_at_k0"], 0.0)
        self.assertLess(crsi["backward_growth_at_minus_k0"], 0.0)
        self.assertFalse(contract["full_q1_dispersion_runtime_oracle"])
        self.assertFalse(contract["qualifying_evidence"])

    def test_crpai_distribution_and_signed_branch_mapping(self) -> None:
        isotropic = q007.isotropic_kappa_distribution(1.0, 0.0, 300.0, 1.75)
        for role, xi, branch in (("prolate", 0.99, "-1"), ("oblate", 1.01, "1")):
            with self.subTest(role=role):
                anisotropic = q007.anisotropic_kappa_distribution(
                    1.0, 0.0, 0.0, 0.0, 300.0, 1.75, xi
                )
                self.assertAlmostEqual(anisotropic / isotropic, xi * xi)
                case = q007.analytical_contract()["crpai"][role]
                self.assertAlmostEqual(case["athenak_transverse_anisotropy_scale"],
                                       1.0 / xi)
                self.assertEqual(case["unstable_signed_branch"], branch)
                self.assertGreater(case["signed_branch_growth"][branch], 0.0)
                self.assertEqual(
                    case["handedness_mapping"],
                    "blocked_pending_manuscript_text_caption_review",
                )
        self.assertFalse(q007.analytical_contract()["crpai_handedness_claimed"])

    def test_three_cycle_zero_decks_freeze_true_deltaf_parser_carriers(self) -> None:
        decks = q007.validate_decks()
        self.assertEqual(
            [item["case"] for item in decks],
            ["crsi", "crpai_prolate", "crpai_oblate"],
        )
        self.assertTrue(all(item["cycle_zero_only"] for item in decks))
        self.assertTrue(all(item["true_deltaf_parser_path"] for item in decks))
        self.assertTrue(all(item["exact_isothermal_mhd"] for item in decks))
        self.assertTrue(
            all(item["source_local_placeholder_particles_total"] == 8 for item in decks)
        )
        self.assertTrue(
            all(item["paper_particles_per_cell_total"] == 2048 for item in decks)
        )
        self.assertTrue(all(not item["physical_loading_implemented"] for item in decks))
        self.assertTrue(all(not item["runtime_evolution_admitted"] for item in decks))
        self.assertTrue(all(not item["qualifying_evidence"] for item in decks))

    def test_deck_drift_fails_closed(self) -> None:
        mutations = (
            ("crsi", "nlim       = 0", "nlim       = 1", "time/nlim"),
            (
                "crsi",
                "couple_moments_energy_to_mhd      = false",
                "couple_moments_energy_to_mhd      = true",
                "couple_moments_energy_to_mhd",
            ),
            (
                "crpai_prolate",
                "pic_deltaf_f0                     = kappa_aniso",
                "pic_deltaf_f0                     = kappa_iso",
                "pic_deltaf_f0",
            ),
            (
                "crpai_oblate",
                "blocked_pending_manuscript_text_caption_review",
                "reviewed_right_handed",
                "handedness_mapping",
            ),
        )
        for case, old, new, expected in mutations:
            with self.subTest(case=case, expected=expected):
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "candidate.athinput"
                    text = q007.DECKS[case].read_text(encoding="utf-8")
                    self.assertIn(old, text)
                    path.write_text(text.replace(old, new, 1), encoding="utf-8")
                    with self.assertRaisesRegex(q007.ContractError, expected):
                        q007.validate_deck(path, case)

    def test_source_contract_is_additive_and_narrow(self) -> None:
        source = q007.validate_source_contract()
        self.assertTrue(source["compilation_unit_registered"])
        self.assertTrue(source["fresh_and_restart_dispatch_registered"])
        self.assertTrue(source["narrow_exact_isothermal_true_deltaf_parser_allowance"])
        self.assertTrue(source["cycle_zero_guarded"])

    def test_report_keeps_qualification_boundaries_open(self) -> None:
        report = q007.build_preparation_report()
        boundary = report["nonqualification_boundary"]
        self.assertTrue(boundary["deck_source_freeze"])
        self.assertTrue(boundary["true_deltaf_parser_path"])
        self.assertTrue(boundary["exact_isothermal_mhd_cycle_zero_startup"])
        for name in (
            "paper_log_bin_weighted_loading",
            "paper_random_phase_four_branch_wave_spectrum",
            "runtime_evolution",
            "full_q1_dispersion_runtime_oracle",
            "crpai_handedness_mapping_reviewed",
            "mpi_qualification",
            "gpu_qualification",
            "frontier_authorization",
            "external_review",
            "qualifying_evidence",
        ):
            self.assertFalse(boundary[name])
        self.assertFalse(report["claim_closure"])
        self.assertFalse(report["frontier_authorization"])

    def test_strict_parser_rejects_duplicate_parameter(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.athinput"
            path.write_text("<mesh>\nnx1 = 32\nnx1 = 32\n", encoding="utf-8")
            with self.assertRaisesRegex(q007.ContractError, "duplicate parameter"):
                q007.parse_athinput(path)

    def test_readiness_sidecar_binds_only_new_q007_artifacts(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-007")
        self.assertEqual(sidecar["qualification_effect"], q007.QUALIFICATION_EFFECT)
        self.assertFalse(sidecar["claim_closure"])
        self.assertEqual(sidecar["frontier_authorization"], "not_bound")
        expected_paths = {
            "src/pgen/tests/q007_paper_deltaf_linear.hpp",
            "src/pgen/tests/q007_paper_deltaf_linear.cpp",
            "inputs/tests/pic_q007_paper_crsi_linear_preparation.athinput",
            "inputs/tests/pic_q007_paper_crpai_linear_prolate_preparation.athinput",
            "inputs/tests/pic_q007_paper_crpai_linear_oblate_preparation.athinput",
            "tst/publication/analyze_q007_paper_deltaf_linear_preparation.py",
            "tst/publication/test_analyze_q007_paper_deltaf_linear_preparation.py",
        }
        bindings = sidecar["artifact_bindings"]
        self.assertEqual(set(bindings), expected_paths)
        for relative, expected in bindings.items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected)
        boundary = sidecar["nonqualification_boundary"]
        self.assertFalse(boundary["qualifying_evidence"])
        self.assertFalse(boundary["runtime_evolution"])
        self.assertFalse(boundary["paper_log_bin_weighted_loading"])
        self.assertFalse(boundary["paper_random_phase_four_branch_wave_spectrum"])


if __name__ == "__main__":
    unittest.main()
