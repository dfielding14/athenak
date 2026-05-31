#!/usr/bin/env python3
"""Regression tests for the fail-closed Q-029 Hall Q-022 prerequisite route."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.publication.pic_qualification_manifest import validate_schema  # noqa: E402


READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
SUCCESSOR_PATH = (
    READINESS_DIR / "q029_hall_bell_q022_prerequisite_successor_2026-05-30.json"
)
Q022_PATH = READINESS_DIR / "q022_independent_comparison_preregistration_2026-05-30.json"
PROVENANCE_PATH = READINESS_DIR / "q022_dataset_provenance_manifest_2026-05-30.json"
EQUATION_PATH = READINESS_DIR / "q022_xcmp_ext_hall_bell_equation_map_2026-05-30.json"
TOLERANCE_PATH = (
    READINESS_DIR / "q022_xcmp_ext_hall_bell_tolerance_table_2026-05-30.json"
)
EQUATION_SCHEMA_PATH = (
    READINESS_DIR / "schemas" / "q022_equation_normalization_map.schema.json"
)
TOLERANCE_SCHEMA_PATH = READINESS_DIR / "schemas" / "q022_tolerance_table.schema.json"
Q023_DRAFTS_PATH = READINESS_DIR / "q023_campaign_drafts_2026-05-30.json"


def _load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q029HallQ022PrerequisiteTests(unittest.TestCase):
    def test_dedicated_q022_route_stays_fail_closed(self) -> None:
        successor = _load(SUCCESSOR_PATH)
        q022 = _load(Q022_PATH)
        provenance = _load(PROVENANCE_PATH)
        equation = _load(EQUATION_PATH)
        tolerance = _load(TOLERANCE_PATH)
        comparison = next(
            item for item in q022["comparison_matrix"]
            if item["comparison_id"] == "XCMP-EXT-HALL-BELL"
        )
        dataset = next(
            item for item in provenance["dataset_candidates"]
            if item["comparison_id"] == "XCMP-EXT-HALL-BELL"
        )
        contract = successor["q022_prerequisite_contract"]

        validate_schema(equation, _load(EQUATION_SCHEMA_PATH))
        validate_schema(tolerance, _load(TOLERANCE_SCHEMA_PATH))
        self.assertEqual(comparison["dataset_provenance_id"], dataset["dataset_id"])
        self.assertEqual(comparison["dataset_provenance_id"],
                         equation["dataset_provenance_id"])
        self.assertEqual(comparison["dataset_provenance_id"],
                         tolerance["dataset_provenance_id"])
        self.assertEqual(contract["comparison_id"], comparison["comparison_id"])
        self.assertEqual(contract["equation_map_status"], equation["map_status"])
        self.assertEqual(contract["tolerance_table_status"],
                         tolerance["freeze_status"])
        self.assertEqual(dataset["extraction_status"],
                         "blocked_extraction_input_unavailable")
        self.assertEqual(equation["matched_equations"], [])
        self.assertEqual(equation["unit_map"], {})
        self.assertEqual(equation["normalization_map"], {})
        self.assertEqual(equation["parameter_overlap"], {})
        self.assertEqual(tolerance["rows"], [])
        self.assertFalse(successor["claim_closure"])

    def test_successor_binds_exact_local_and_q022_artifacts(self) -> None:
        successor = _load(SUCCESSOR_PATH)
        for bindings_name in ("artifact_bindings", "q022_prerequisite_bindings"):
            for path, expected_sha256 in successor[bindings_name].items():
                with self.subTest(bindings=bindings_name, path=path):
                    self.assertEqual(_sha256(REPO_ROOT / path), expected_sha256)

    def test_q023_hall_draft_inventories_existing_q029_preparation(self) -> None:
        campaigns = _load(Q023_DRAFTS_PATH)["campaigns"]
        campaign = next(
            item for item in campaigns if item["campaign_id"] == "Q023-EXT-HALL-BELL"
        )
        self.assertTrue({
            "Q023-INPUT-Q029-HALL-BELL-LINEAR-1D-CANDIDATE",
            "Q023-INPUT-Q029-HALL-BELL-LINEAR-2D-CANDIDATE",
            "Q023-INPUT-Q029-HALL-BELL-LINEAR-3D-CANDIDATE",
        } <= set(campaign["local_input_ids"]))
        self.assertIn("Q023-ANALYZER-Q029-HALL-BELL-GRID",
                      campaign["local_analyzer_ids"])


if __name__ == "__main__":
    unittest.main()
