#!/usr/bin/env python3
"""Fail-closed checks for the additive paper VL2/TSC source-local tranche."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
RECORD = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "pic_vl2_tsc_additive_successor_registration_2026-06-01.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PicVL2TSCAdditiveSuccessorRegistrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = json.loads(RECORD.read_text(encoding="utf-8"))

    def test_registration_is_nonqualifying_and_uses_required_roots(self) -> None:
        self.assertEqual(
            self.record["qualification_effect"],
            "source_local_registration_only_no_execution_authorization_no_claim_closure",
        )
        self.assertEqual(
            self.record["artifact_root"],
            "/lustre/orion/ast207/proj-shared/dfielding/PIC",
        )
        frontier = self.record["frontier_binding"]
        self.assertEqual(frontier["node_hour_cap"], 10000.0)
        self.assertFalse(frontier["new_qualification_launch_authorized"])
        self.assertEqual(frontier["project_home_root"], "/ccs/home/dfielding/entity")

    def test_historical_files_retain_frozen_bytes(self) -> None:
        for relative, expected in self.record["historical_files_preserved"].items():
            with self.subTest(relative=relative):
                self.assertEqual(_sha256(REPO_ROOT / relative), expected)

    def test_additive_file_bindings_remain_as_v1_chronology(self) -> None:
        bindings = self.record["additive_vl2_tsc_files"]
        self.assertGreater(len(bindings), 0)
        for relative, expected in bindings.items():
            with self.subTest(relative=relative):
                self.assertTrue((REPO_ROOT / relative).is_file())
                self.assertRegex(expected, r"^[0-9a-f]{64}$")

    def test_runtime_oracle_includes_periodic_and_cross_rank_cases(self) -> None:
        oracle = self.record["runtime_oracle"]
        self.assertEqual(oracle["serial_cases"], list("abcdefghijkl"))
        self.assertEqual(oracle["mpi2_cases"], list("abcdefghijkl"))
        self.assertEqual(oracle["periodic_image_case"], "g")
        self.assertEqual(oracle["same_gid_periodic_image_case"], "h")
        self.assertEqual(oracle["physical_boundary_normalization_case"], "i")
        self.assertEqual(
            oracle["physical_boundary_modes"],
            {"i": "reflect", "k": "outflow", "l": "inflow"},
        )
        self.assertEqual(oracle["three_dimensional_periodic_edge_case"], "j")
        self.assertEqual(oracle["cross_rank_split_cases"], ["b", "c", "d", "g"])
        self.assertEqual(
            oracle["expected_raw_totals"],
            {
                "a": 1.0,
                "b": 1.10125,
                "c": 1.15625,
                "d": 0.88375,
                "e": 0.98,
                "f": 1.0,
                "g": 1.14,
                "h": 1.0,
                "i": 1.0,
                "j": 0.86,
            },
        )


if __name__ == "__main__":
    unittest.main()
