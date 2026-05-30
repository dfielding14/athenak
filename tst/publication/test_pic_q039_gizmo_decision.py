#!/usr/bin/env python3
"""Regression tests for the fail-closed PIC Q-039 GIZMO/RSOL decision sidecar."""

from __future__ import annotations

import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
Q039_PATH = READINESS_DIR / "q039_gizmo_rsol_decision_2026-05-30.json"
Q022_PATH = READINESS_DIR / "q022_independent_comparison_preregistration_2026-05-30.json"
CLAIMS_PATH = READINESS_DIR / "claims_registry.json"
ARTIFACTS_PATH = READINESS_DIR / "external_artifacts.json"

CLAIM_ID = "CLAIM-XCODE-GIZMO-RSOL-DECISION-001"
COMPARISON_ID = "XCMP-GIZMO-RSOL-DECISION"
REFERENCE_ID = "ji_hopkins_2022_arxiv_2111.14704"
EXTERNAL_ARTIFACT_ID = "CROSS_CODE_REFERENCE_DATA"


def _load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


class PicQ039GizmoDecisionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.q039 = _load(Q039_PATH)
        self.q022 = _load(Q022_PATH)
        self.claims = _load(CLAIMS_PATH)
        self.artifacts = _load(ARTIFACTS_PATH)

    def test_sidecar_remains_open_and_nonqualifying(self) -> None:
        self.assertEqual(self.q039["gate"], "Q-039")
        self.assertEqual(
            self.q039["qualification_effect"],
            "local_fail_closed_decision_scaffold_only",
        )
        self.assertEqual(self.q039["disposition"], "open")
        self.assertIn("open", self.q039["q039_status"])

    def test_linked_ids_match_frozen_q022_inventory(self) -> None:
        linked = self.q039["linked_records"]
        self.assertEqual(linked["claim_id"], CLAIM_ID)
        self.assertEqual(linked["comparison_id"], COMPARISON_ID)
        self.assertEqual(linked["reference_id"], REFERENCE_ID)
        self.assertEqual(linked["external_artifact_id"], EXTERNAL_ARTIFACT_ID)

        claims = {item["claim_id"]: item for item in self.claims["claims"]}
        self.assertEqual(claims[CLAIM_ID]["disposition"], "open")

        references = {
            item["reference_id"]: item for item in self.q022["reference_families"]
        }
        self.assertEqual(
            references[REFERENCE_ID]["retrieved_artifact_status"],
            "local_staged_checksummed_pending_orion_post_copy_verification_"
            "and_data_extraction",
        )
        comparisons = {
            item["comparison_id"]: item for item in self.q022["comparison_matrix"]
        }
        self.assertEqual(comparisons[COMPARISON_ID]["claim_ids"], [CLAIM_ID])
        self.assertEqual(comparisons[COMPARISON_ID]["reference_ids"], [REFERENCE_ID])

        artifacts = {
            item["artifact_id"]: item for item in self.artifacts["artifacts"]
        }
        self.assertEqual(
            artifacts[EXTERNAL_ARTIFACT_ID]["status"],
            "local_staged_checksummed_reference_sources_pending_orion_copy_"
            "post_copy_verification_and_dataset_extraction",
        )

    def test_sources_remain_staged_pending_orion_copy_and_verification(self) -> None:
        preparation = self.q039["reference_source_preparation"]
        self.assertEqual(
            preparation["status"],
            "local_staged_checksummed_pending_orion_copy_post_copy_verification_"
            "mapping_and_data_extraction",
        )
        self.assertTrue(
            all(
                artifact["status"]
                == "local_staged_checksummed_pending_orion_copy_and_post_copy_"
                "checksum_verification"
                for artifact in preparation["required_artifacts"]
            )
        )

    def test_matrix_accepts_no_exact_overlap_before_verification_and_mapping(self) -> None:
        matrix = self.q039["decision_matrix"]
        dimensions = {item["dimension"] for item in matrix}
        self.assertEqual(len(dimensions), len(matrix))
        self.assertIn("sun_bai_artificial_C_reproduction_contract", dimensions)
        self.assertIn("steady_state_rsol_invariance", dimensions)
        self.assertTrue(all(not item["accepted_exact_overlap"] for item in matrix))
        self.assertEqual(
            self.q039["frozen_open_placeholders"]["accepted_exact_overlap_count"],
            0,
        )

    def test_fail_closed_rule_blocks_manuscript_claims(self) -> None:
        rule = self.q039["fail_closed_rule"]
        self.assertEqual(rule["status"], "active")
        self.assertEqual(rule["manuscript_claim_use"], "prohibited")
        self.assertIn("orion_reference_source_copy_complete",
                      rule["lift_only_after"])
        self.assertIn("orion_post_copy_sha256sums_verified",
                      rule["lift_only_after"])

        exclusion = self.q039["decision_paths"]["documented_exclusion"]
        self.assertEqual(exclusion["status"],
                         "provisional_fail_closed_exclusion_only")
        self.assertIn("not the reviewer-approved documented exclusion",
                      exclusion["rule"])
        self.assertIn("external_reviewer_approval",
                      exclusion["closure_requires"])

    def test_placeholders_and_mapping_prerequisites_stay_open(self) -> None:
        placeholders = self.q039["frozen_open_placeholders"]
        for key in (
            "reference_artifacts",
            "equation_map",
            "normalization_map",
            "bounded_metrics",
            "numerical_tolerances",
        ):
            with self.subTest(key=key):
                self.assertIn("open", placeholders[key])
        prerequisites = self.q039["mapping_prerequisites_before_bounded_comparison"]
        self.assertGreaterEqual(len(prerequisites), 6)
        self.assertTrue(any("checksum" in item for item in prerequisites))
        self.assertTrue(any("tolerances" in item for item in prerequisites))
        self.assertTrue(any("reviewer" in item for item in prerequisites))


if __name__ == "__main__":
    unittest.main()
