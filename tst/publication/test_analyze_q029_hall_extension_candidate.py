#!/usr/bin/env python3
"""Focused tests for the Q-029 launch-blocked Hall normalization candidate."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from tst.publication import analyze_q029_hall_extension_candidate as hall


REPO_ROOT = Path(__file__).resolve().parents[2]
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q029_hall_extension_parser_hardening_successor_v2_2026-06-01.json"
)
DERIVATION = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q029_experimental_hall_normalization_derivation_2026-05-30.md"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q029HallExtensionCandidateTests(unittest.TestCase):
    def test_exact_source_bindings_and_launch_blocked_deck_validate(self) -> None:
        self.assertEqual(hall.validate_source_bindings(), hall._SOURCE_BINDINGS)
        deck = hall.validate_candidate_deck()
        self.assertEqual(deck["launch_status"], hall.LAUNCH_STATUS)
        self.assertEqual(deck["qualification_effect"], "none")
        self.assertFalse(deck["qualifying_evidence"])
        self.assertEqual(deck["normalization"]["u_a0"], 1.0)
        self.assertEqual(deck["normalization"]["c_e0"], 1.0)
        self.assertEqual(deck["normalization"]["chi_h_fiducial"], 0.5)
        self.assertEqual(deck["normalization"]["chi_h_candidate_grid"],
                         [0.0, 0.25, 0.5, 1.0])

    def test_complete_synthetic_contract_grid_is_nonqualifying(self) -> None:
        report = hall.analyze_synthetic_contract_bundle(
            hall.build_synthetic_contract_bundle()
        )
        self.assertTrue(report["synthetic_normalization_contract_consistent"])
        self.assertEqual(report["sample_count"], 9)
        self.assertFalse(report["qualifying_evidence"])
        self.assertEqual(report["qualification_effect"], hall.QUALIFICATION_EFFECT)
        self.assertIn("not_qualifying_evidence", report["status"])

    def test_numeric_mismatch_is_reported_without_qualification(self) -> None:
        bundle = hall.build_synthetic_contract_bundle()
        bundle["samples"][0]["observed_c_e_over_c_e0"] += 0.125
        report = hall.analyze_synthetic_contract_bundle(bundle)
        self.assertFalse(report["synthetic_normalization_contract_consistent"])
        self.assertFalse(report["qualifying_evidence"])
        self.assertIn("mismatch_not_qualifying_evidence", report["status"])

    def test_incomplete_and_duplicate_synthetic_grids_fail_closed(self) -> None:
        bundle = hall.build_synthetic_contract_bundle()
        bundle["samples"].pop()
        with self.assertRaisesRegex(hall.ContractError, "grid is incomplete"):
            hall.analyze_synthetic_contract_bundle(bundle)

        bundle = hall.build_synthetic_contract_bundle()
        bundle["samples"].append(copy.deepcopy(bundle["samples"][0]))
        with self.assertRaisesRegex(hall.ContractError, "duplicate"):
            hall.analyze_synthetic_contract_bundle(bundle)

    def test_runtime_like_extra_field_and_wrong_role_fail_closed(self) -> None:
        bundle = hall.build_synthetic_contract_bundle()
        bundle["runtime_artifact_path"] = "/tmp/not-admitted"
        with self.assertRaisesRegex(hall.ContractError, "bundle keys"):
            hall.analyze_synthetic_contract_bundle(bundle)

        bundle = hall.build_synthetic_contract_bundle()
        bundle["artifact_role"] = "qualifying_runtime_evidence"
        with self.assertRaisesRegex(hall.ContractError, "not synthetic"):
            hall.analyze_synthetic_contract_bundle(bundle)

    def test_source_and_deck_drift_fail_closed(self) -> None:
        with mock.patch.dict(
            hall._SOURCE_BINDINGS,
            {"src/mhd/mhd_tasks.cpp": "0" * 64},
        ):
            with self.assertRaisesRegex(hall.ContractError, "source checksum mismatch"):
                hall.validate_source_bindings()

        with tempfile.TemporaryDirectory() as directory:
            mutated = Path(directory) / "candidate.athinput"
            mutated.write_text(
                hall.DECK.read_text(encoding="utf-8") + "\n# mutation\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(hall.ContractError, "deck checksum mismatch"):
                hall.validate_candidate_deck(mutated)

    def test_accidentally_available_generator_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            dispatch = Path(directory) / "pgen.cpp"
            dispatch.write_text(
                'const char *name = "q029_extended_hall_normalization_open";\n',
                encoding="utf-8",
            )
            with mock.patch.object(hall, "_PGEN_DISPATCH", dispatch):
                with self.assertRaisesRegex(hall.ContractError,
                                            "unexpectedly became available"):
                    hall.validate_launch_block()

    def test_readiness_sidecar_binds_q029_only_artifacts(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-029")
        self.assertEqual(sidecar["qualification_effect"],
                         "source_local_launch_blocked_preparation_only")
        self.assertFalse(sidecar["claim_closure"])
        bindings = sidecar["artifact_bindings"]
        expected_paths = {
            "inputs/tests/pic_q029_extended_hall_normalization_candidate.athinput",
            "tst/publication/analyze_q029_hall_extension_candidate.py",
            "tst/publication/test_analyze_q029_hall_extension_candidate.py",
            "tst/publication/readiness/"
            "q029_experimental_hall_normalization_derivation_2026-05-30.md",
        }
        self.assertEqual(set(bindings), expected_paths)
        for relative, expected_sha256 in bindings.items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected_sha256)
        self.assertEqual(sidecar["implementation_source_bindings"],
                         hall._SOURCE_BINDINGS)
        self.assertEqual(_sha256(DERIVATION), bindings[str(DERIVATION.relative_to(
            REPO_ROOT
        ))])


if __name__ == "__main__":
    unittest.main()
