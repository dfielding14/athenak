#!/usr/bin/env python3
"""Focused tests for the Q-032 launch-blocked damped-CRSI candidate."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from tst.publication import analyze_q032_plotnikov_damped_crsi_candidate as damping


REPO_ROOT = Path(__file__).resolve().parents[2]
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q032_plotnikov_damped_crsi_parser_hardening_successor_v2_2026-06-01.json"
)
BOUNDARY_NOTE = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q032_reduced_static_neutral_plotnikov_boundary_2026-05-30.md"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q032PlotnikovDampedCRSICandidateTests(unittest.TestCase):
    def test_exact_source_prerequisite_and_launch_blocked_deck_validate(self) -> None:
        self.assertEqual(damping.validate_source_bindings(), damping._SOURCE_BINDINGS)
        self.assertEqual(damping.validate_prerequisite_boundaries(),
                         damping._PREREQUISITE_BINDINGS)
        deck = damping.validate_candidate_deck()
        self.assertEqual(deck["launch_status"], damping.LAUNCH_STATUS)
        self.assertEqual(deck["qualification_effect"], "none")
        self.assertFalse(deck["qualifying_evidence"])
        self.assertFalse(deck["matched_plotnikov_qualification"])
        self.assertEqual(deck["nu_in_fiducial"], 0.7)
        self.assertEqual(deck["collision_rate_candidate_grid"],
                         [0.0, 0.2, 0.7, 1.4])

    def test_complete_synthetic_contract_grid_is_nonqualifying(self) -> None:
        report = damping.analyze_synthetic_contract_bundle(
            damping.build_synthetic_contract_bundle()
        )
        self.assertTrue(report["synthetic_static_neutral_contract_consistent"])
        self.assertEqual(report["sample_count"], 9)
        self.assertFalse(report["qualifying_evidence"])
        self.assertFalse(report["matched_plotnikov_qualification"])
        self.assertEqual(report["qualification_effect"],
                         damping.QUALIFICATION_EFFECT)
        self.assertIn("not_plotnikov_qualification", report["status"])

    def test_numeric_mismatch_is_reported_without_qualification(self) -> None:
        bundle = damping.build_synthetic_contract_bundle()
        bundle["samples"][0]["energy_after"] += 0.125
        report = damping.analyze_synthetic_contract_bundle(bundle)
        self.assertFalse(report["synthetic_static_neutral_contract_consistent"])
        self.assertFalse(report["qualifying_evidence"])
        self.assertFalse(report["matched_plotnikov_qualification"])
        self.assertIn("mismatch_not_plotnikov_qualification", report["status"])

    def test_incomplete_and_duplicate_synthetic_grids_fail_closed(self) -> None:
        bundle = damping.build_synthetic_contract_bundle()
        bundle["samples"].pop()
        with self.assertRaisesRegex(damping.ContractError, "grid is incomplete"):
            damping.analyze_synthetic_contract_bundle(bundle)

        bundle = damping.build_synthetic_contract_bundle()
        bundle["samples"].append(copy.deepcopy(bundle["samples"][0]))
        with self.assertRaisesRegex(damping.ContractError, "duplicate"):
            damping.analyze_synthetic_contract_bundle(bundle)

    def test_runtime_like_extra_field_and_wrong_role_fail_closed(self) -> None:
        bundle = damping.build_synthetic_contract_bundle()
        bundle["plotnikov_runtime_artifact_path"] = "/tmp/not-admitted"
        with self.assertRaisesRegex(damping.ContractError, "bundle keys"):
            damping.analyze_synthetic_contract_bundle(bundle)

        bundle = damping.build_synthetic_contract_bundle()
        bundle["artifact_role"] = "matched_plotnikov_runtime_evidence"
        with self.assertRaisesRegex(damping.ContractError, "not synthetic"):
            damping.analyze_synthetic_contract_bundle(bundle)

    def test_source_prerequisite_and_deck_drift_fail_closed(self) -> None:
        with mock.patch.dict(
            damping._SOURCE_BINDINGS,
            {"src/mhd/mhd_tasks.cpp": "0" * 64},
        ):
            with self.assertRaisesRegex(damping.ContractError,
                                        "source checksum mismatch"):
                damping.validate_source_bindings()

        prerequisite = next(iter(damping._PREREQUISITE_BINDINGS))
        with mock.patch.dict(damping._PREREQUISITE_BINDINGS,
                             {prerequisite: "0" * 64}):
            with self.assertRaisesRegex(damping.ContractError,
                                        "prerequisite checksum mismatch"):
                damping.validate_prerequisite_boundaries()

        with tempfile.TemporaryDirectory() as directory:
            mutated = Path(directory) / "candidate.athinput"
            mutated.write_text(
                damping.DECK.read_text(encoding="utf-8") + "\n# mutation\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(damping.ContractError,
                                        "deck checksum mismatch"):
                damping.validate_candidate_deck(mutated)

    def test_accidentally_available_generator_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            dispatch = Path(directory) / "pgen.cpp"
            dispatch.write_text(
                'const char *name = "q032_plotnikov_damped_crsi_open";\n',
                encoding="utf-8",
            )
            with mock.patch.object(damping, "_PGEN_DISPATCH", dispatch):
                with self.assertRaisesRegex(damping.ContractError,
                                            "unexpectedly became available"):
                    damping.validate_launch_block()

    def test_readiness_sidecar_binds_q032_only_artifacts(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-032")
        self.assertEqual(sidecar["qualification_effect"],
                         "source_local_launch_blocked_preparation_only")
        self.assertFalse(sidecar["claim_closure"])
        bindings = sidecar["artifact_bindings"]
        expected_paths = {
            "inputs/tests/pic_q032_plotnikov_damped_crsi_candidate.athinput",
            "tst/publication/analyze_q032_plotnikov_damped_crsi_candidate.py",
            "tst/publication/test_analyze_q032_plotnikov_damped_crsi_candidate.py",
            "tst/publication/readiness/"
            "q032_reduced_static_neutral_plotnikov_boundary_2026-05-30.md",
        }
        self.assertEqual(set(bindings), expected_paths)
        for relative, expected_sha256 in bindings.items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected_sha256)
        self.assertEqual(sidecar["implementation_source_bindings"],
                         damping._SOURCE_BINDINGS)
        self.assertEqual(sidecar["prerequisite_record_bindings"],
                         damping._PREREQUISITE_BINDINGS)
        self.assertEqual(
            _sha256(BOUNDARY_NOTE),
            bindings[str(BOUNDARY_NOTE.relative_to(REPO_ROOT))],
        )

    def test_sidecar_retains_explicit_nonqualification_boundary(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        launch_block = sidecar["launch_block"]
        self.assertFalse(launch_block["analyzer_always_emits_qualifying_evidence"])
        self.assertFalse(launch_block["analyzer_can_emit_matched_plotnikov_qualification"])
        self.assertIn("extracted reference dataset", sidecar["explicitly_not_claimed"])
        self.assertIn("unstable-bandwidth agreement",
                      sidecar["explicitly_not_claimed"])


if __name__ == "__main__":
    unittest.main()
