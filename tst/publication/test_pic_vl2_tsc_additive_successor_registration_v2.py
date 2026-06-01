#!/usr/bin/env python3
"""Fail-closed checks for the corrected additive paper VL2/TSC tranche."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
RECORD = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "pic_vl2_tsc_additive_successor_registration_v2_2026-06-01.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PicVL2TSCAdditiveSuccessorRegistrationV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = json.loads(RECORD.read_text(encoding="utf-8"))

    def test_registration_is_nonqualifying_and_supersedes_v1_mechanically(self) -> None:
        self.assertEqual(self.record["schema_version"], 2)
        self.assertEqual(
            self.record["qualification_effect"],
            "source_local_registration_only_no_execution_authorization_no_claim_closure",
        )
        prior = self.record["supersedes_mechanically"]
        self.assertEqual(
            _sha256(REPO_ROOT / prior["path"]),
            prior["sha256"],
        )
        frontier = self.record["frontier_binding"]
        self.assertEqual(frontier["node_hour_cap"], 10000.0)
        self.assertEqual(
            frontier["bulk_artifact_root"],
            "/lustre/orion/ast207/proj-shared/dfielding/PIC",
        )
        self.assertFalse(frontier["new_qualification_launch_authorized"])

    def test_active_successor_files_are_exactly_bound(self) -> None:
        for relative, expected in self.record["active_successor_files"].items():
            with self.subTest(relative=relative):
                self.assertEqual(_sha256(REPO_ROOT / relative), expected)

    def test_identity_split_is_explicit_and_schema_compatible(self) -> None:
        identity = self.record["identity_split"]
        self.assertEqual(identity["historical_mode"], "paper_mhd_pic")
        self.assertEqual(identity["successor_mode"], "paper_mhd_pic_vl2_tsc")
        self.assertEqual(identity["restart_schema_version"], 7)
        self.assertFalse(identity["restart_schema_changed"])
        self.assertTrue(identity["successor_selects_vl2_coefficients"])
        self.assertTrue(identity["successor_requires_tsc_order2"])
        self.assertTrue(identity["successor_deltaf_staging_implemented"])
        self.assertFalse(identity["successor_expanding_box_staging_implemented"])

    def test_runtime_oracle_includes_signed_deltaf_and_exact_remote_routes(self) -> None:
        oracle = self.record["runtime_oracle"]
        self.assertEqual(oracle["serial_cases"], list("abcdefghijklmnop"))
        self.assertEqual(oracle["mpi2_cases"], list("abcdefghijklmnop"))
        self.assertEqual(oracle["mpi3_cases"], ["j", "p"])
        self.assertEqual(
            oracle["exact_remote_receiver_routes"],
            {
                "mpi2_g": "g_periodic_x1",
                "mpi2_n": "g_periodic_x1",
                "mpi3_j": "j_periodic_x1_x3_edge",
                "mpi3_p": "j_periodic_x1_x3_edge",
            },
        )
        self.assertEqual(
            oracle["deltaf_scale_pairs"],
            {
                "m": ["a", 0.5],
                "n": ["g", -0.25],
                "o": ["k", 0.5],
                "p": ["j", -0.25],
            },
        )


if __name__ == "__main__":
    unittest.main()
