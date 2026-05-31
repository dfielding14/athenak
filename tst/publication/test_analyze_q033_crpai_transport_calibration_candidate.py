#!/usr/bin/env python3
"""Focused tests for the Q-033 launch-blocked CRPAI transport candidate."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from tst.publication import analyze_q033_crpai_transport_calibration_candidate as q033


REPO_ROOT = Path(__file__).resolve().parents[2]
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q033_crpai_transport_calibration_source_local_candidate_2026-05-30.json"
)
APPLICABILITY = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q033_crpai_transport_calibration_applicability_2026-05-30.md"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q033CRPAITransportCalibrationCandidateTests(unittest.TestCase):
    def test_exact_boundaries_and_launch_blocked_deck_validate(self) -> None:
        self.assertEqual(q033.validate_source_bindings(), q033._SOURCE_BINDINGS)
        self.assertEqual(
            q033.validate_boundary_artifact_bindings(),
            q033._BOUNDARY_ARTIFACT_BINDINGS,
        )
        deck = q033.validate_candidate_deck()
        self.assertEqual(deck["launch_status"], q033.LAUNCH_STATUS)
        self.assertEqual(deck["qualification_effect"], "none")
        self.assertFalse(deck["qualifying_evidence"])

    def test_complete_synthetic_contract_is_nonqualifying(self) -> None:
        report = q033.analyze_synthetic_contract_bundle(
            q033.build_synthetic_contract_bundle()
        )
        self.assertTrue(report["synthetic_transport_contract_consistent"])
        self.assertEqual(report["series_count"], 4)
        self.assertFalse(report["qualifying_evidence"])
        self.assertFalse(report["physical_calibration_claimed"])
        self.assertEqual(report["qualification_effect"], q033.QUALIFICATION_EFFECT)
        self.assertIn("not_physical_calibration", report["status"])
        self.assertIn("not_qualifying_evidence", report["status"])
        self.assertTrue(
            report["scaling"]["synthetic_oracle_only_not_physical_scaling_law"]
        )

    def test_numeric_mismatch_never_claims_physical_calibration(self) -> None:
        bundle = q033.build_synthetic_contract_bundle()
        bundle["series"][0]["tail_samples"][-1][
            "parallel_displacement_variance"
        ] *= 1.125
        report = q033.analyze_synthetic_contract_bundle(bundle)
        self.assertFalse(report["synthetic_transport_contract_consistent"])
        self.assertFalse(report["qualifying_evidence"])
        self.assertFalse(report["physical_calibration_claimed"])
        self.assertIn("mismatch_not_physical_calibration", report["status"])

    def test_incomplete_and_duplicate_series_grids_fail_closed(self) -> None:
        bundle = q033.build_synthetic_contract_bundle()
        bundle["series"].pop()
        with self.assertRaisesRegex(q033.ContractError, "grid is incomplete"):
            q033.analyze_synthetic_contract_bundle(bundle)

        bundle = q033.build_synthetic_contract_bundle()
        bundle["series"].append(copy.deepcopy(bundle["series"][0]))
        with self.assertRaisesRegex(q033.ContractError, "duplicate"):
            q033.analyze_synthetic_contract_bundle(bundle)

    def test_malformed_spectrum_and_noncanonical_k_grid_fail_closed(self) -> None:
        bundle = q033.build_synthetic_contract_bundle()
        del bundle["series"][0]["tail_samples"][0]["spectrum"][0]["forward_left"]
        with self.assertRaisesRegex(q033.ContractError, "spectrum keys"):
            q033.analyze_synthetic_contract_bundle(bundle)

        bundle = q033.build_synthetic_contract_bundle()
        bundle["series"][0]["tail_samples"][0]["spectrum"][0]["k"] = 0.5
        with self.assertRaisesRegex(q033.ContractError, "spectrum k grid"):
            q033.analyze_synthetic_contract_bundle(bundle)

    def test_runtime_like_extra_field_wrong_role_and_binding_drift_fail_closed(self) -> None:
        bundle = q033.build_synthetic_contract_bundle()
        bundle["runtime_artifact_root"] = "/tmp/not-admitted"
        with self.assertRaisesRegex(q033.ContractError, "bundle keys"):
            q033.analyze_synthetic_contract_bundle(bundle)

        bundle = q033.build_synthetic_contract_bundle()
        bundle["artifact_role"] = "qualifying_runtime_evidence"
        with self.assertRaisesRegex(q033.ContractError, "not synthetic"):
            q033.analyze_synthetic_contract_bundle(bundle)

        bundle = q033.build_synthetic_contract_bundle()
        bundle["boundary_artifact_sha256"] = dict(
            bundle["boundary_artifact_sha256"]
        )
        first = next(iter(bundle["boundary_artifact_sha256"]))
        bundle["boundary_artifact_sha256"][first] = "0" * 64
        with self.assertRaisesRegex(q033.ContractError, "checksums mismatch"):
            q033.analyze_synthetic_contract_bundle(bundle)

    def test_source_boundary_artifact_and_deck_drift_fail_closed(self) -> None:
        with mock.patch.dict(
            q033._SOURCE_BINDINGS,
            {"src/particles/particles.cpp": "0" * 64},
        ):
            with self.assertRaisesRegex(q033.ContractError, "source checksum mismatch"):
                q033.validate_source_bindings()

        with mock.patch.dict(
            q033._BOUNDARY_ARTIFACT_BINDINGS,
            {
                "inputs/tests/pic_mhd_expanding_box_adaptive_damping_smoke.athinput":
                    "0" * 64
            },
        ):
            with self.assertRaisesRegex(
                q033.ContractError, "boundary artifact checksum mismatch"
            ):
                q033.validate_boundary_artifact_bindings()

        with tempfile.TemporaryDirectory() as directory:
            mutated = Path(directory) / "candidate.athinput"
            mutated.write_text(
                q033.DECK.read_text(encoding="utf-8") + "\n# mutation\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(q033.ContractError, "deck checksum mismatch"):
                q033.validate_candidate_deck(mutated)

    def test_accidentally_available_generator_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            dispatch = Path(directory) / "pgen.cpp"
            dispatch.write_text(
                'const char *name = "q033_crpai_transport_calibration_open";\n',
                encoding="utf-8",
            )
            with mock.patch.object(q033, "_PGEN_DISPATCH", dispatch):
                with self.assertRaisesRegex(
                    q033.ContractError, "unexpectedly became available"
                ):
                    q033.validate_launch_block()

    def test_readiness_sidecar_binds_q033_only_artifacts(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-033")
        self.assertEqual(
            sidecar["qualification_effect"],
            "source_local_launch_blocked_preparation_only",
        )
        self.assertFalse(sidecar["claim_closure"])
        bindings = sidecar["artifact_bindings"]
        expected_paths = {
            "inputs/tests/pic_q033_crpai_transport_calibration_candidate.athinput",
            "tst/publication/analyze_q033_crpai_transport_calibration_candidate.py",
            "tst/publication/test_analyze_q033_crpai_transport_calibration_candidate.py",
            "tst/publication/readiness/"
            "q033_crpai_transport_calibration_applicability_2026-05-30.md",
        }
        self.assertEqual(set(bindings), expected_paths)
        for relative, expected_sha256 in bindings.items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected_sha256)
        self.assertEqual(sidecar["implementation_source_bindings"],
                         q033._SOURCE_BINDINGS)
        self.assertEqual(sidecar["bounded_boundary_artifact_bindings"],
                         q033._BOUNDARY_ARTIFACT_BINDINGS)
        relative = str(APPLICABILITY.relative_to(REPO_ROOT))
        self.assertEqual(_sha256(APPLICABILITY), bindings[relative])


if __name__ == "__main__":
    unittest.main()
