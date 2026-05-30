#!/usr/bin/env python3
"""Regression tests for the bounded-local Q-009 inflow lifetime blocker record."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
RECORD_PATH = READINESS_DIR / "q009_coupled_inflow_lifetime_bounded_local_2026-05-30.json"


def _sha256(relative_path: str) -> str:
    return hashlib.sha256((REPO_ROOT / relative_path).read_bytes()).hexdigest()


class Q009CoupledInflowLifetimeReadinessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = json.loads(RECORD_PATH.read_text(encoding="utf-8"))

    def test_archived_reflect_outflow_sibling_is_byte_bound(self) -> None:
        archived = self.record["archived_reflect_outflow_sibling"]
        self.assertTrue(archived["behavior_preserved"])
        for binding in [
            archived["input"],
            archived["harness"],
            archived["readiness_record"],
        ]:
            with self.subTest(path=binding["path"]):
                self.assertEqual(_sha256(binding["path"]), binding["sha256"])

    def test_inflow_sibling_is_byte_bound(self) -> None:
        sibling = self.record["inflow_sibling"]
        for binding in [sibling["input"], sibling["harness"]]:
            with self.subTest(path=binding["path"]):
                self.assertEqual(_sha256(binding["path"]), binding["sha256"])
        self.assertEqual(
            sibling["supported_physical_boundaries"],
            {"x1": "reflect", "x2": "inflow", "x3": "periodic"},
        )
        self.assertEqual(
            sibling["particle_moment_inflow_contract"],
            "zero_valued_moment_ghosts",
        )

    def test_inflow_deck_changes_only_the_reviewed_sibling_fields(self) -> None:
        archived = (
            REPO_ROOT / "inputs/tests/pic_q009_coupled_boundary_lifetime.athinput"
        ).read_text(encoding="utf-8")
        sibling = (
            REPO_ROOT / "inputs/tests/pic_q009_coupled_inflow_lifetime.athinput"
        ).read_text(encoding="utf-8")
        expected = archived.replace("Reflect/outflow", "Reflect/inflow")
        expected = expected.replace(
            "Q-009 bounded coupled PIC/MHD boundary lifetime stress",
            "Q-009 bounded coupled PIC/MHD inflow lifetime stress",
        )
        expected = expected.replace(
            "pic_q009_coupled_boundary_lifetime",
            "pic_q009_coupled_inflow_lifetime",
        )
        expected = expected.replace("ix2_bc    = outflow", "ix2_bc    = inflow")
        expected = expected.replace("ox2_bc    = outflow", "ox2_bc    = inflow")
        self.assertEqual(sibling, expected)

    def test_record_is_fail_closed_without_asan_success_claim(self) -> None:
        self.assertEqual(
            self.record["qualification_effect"],
            "none_fail_closed_characterization_only",
        )
        self.assertEqual(
            self.record["disposition"],
            "blocked_linear_wave_mhd_inflow_reservoir_not_initialized",
        )
        debug = self.record["debug_serial_characterization"]
        self.assertEqual(debug["result"], "fail_closed")
        self.assertFalse(debug["finite_mesh_state"])
        self.assertFalse(debug["subsequent_stage_time_advance"])
        asan = self.record["fresh_host_gnu_asan_attempt"]
        self.assertEqual(asan["compiler_probe"], "pass")
        self.assertEqual(asan["configure_result"], "pass")
        self.assertEqual(asan["build_result"], "pass")
        self.assertEqual(asan["startup_result"],
                         "pass: built_in_pgens double precision MPI OFF OpenMP OFF")
        for harness in [
            asan["archived_reflect_outflow_harness"],
            asan["inflow_sibling_harness"],
        ]:
            self.assertTrue(harness["result"].startswith("fail_"))
            self.assertFalse(harness["asan_success_claimed"])

    def test_source_constraints_are_q009_inflow_specific(self) -> None:
        constraints = self.record["source_constraints"]
        self.assertFalse(constraints["production_source_edited_by_this_slice"])
        self.assertFalse(constraints["archived_reflect_outflow_files_edited"])
        self.assertFalse(constraints["plan_edited_by_this_slice"])
        self.assertEqual(
            constraints["new_q009_inflow_specific_files_only"],
            [
                "inputs/tests/pic_q009_coupled_inflow_lifetime.athinput",
                "tst/scripts/particles/pic_q009_coupled_inflow_lifetime.py",
                "tst/publication/readiness/"
                "q009_coupled_inflow_lifetime_bounded_local_2026-05-30.json",
                "tst/publication/test_q009_coupled_inflow_lifetime_readiness.py",
            ],
        )


if __name__ == "__main__":
    unittest.main()
