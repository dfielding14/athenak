#!/usr/bin/env python3
"""Regression tests for the bounded-local repaired Q-009 inflow successor."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
RECORD_PATH = (
    READINESS_DIR
    / "q009_coupled_inflow_lifetime_repaired_bounded_local_2026-05-30.json"
)


def _sha256(relative_path: str) -> str:
    return hashlib.sha256((REPO_ROOT / relative_path).read_bytes()).hexdigest()


class Q009CoupledInflowLifetimeRepairedReadinessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = json.loads(RECORD_PATH.read_text(encoding="utf-8"))

    def test_predecessor_blocker_is_preserved(self) -> None:
        predecessor = self.record["predecessor_blocker"]
        self.assertTrue(predecessor["preserved_as_immutable_chronology"])
        self.assertEqual(_sha256(predecessor["path"]), predecessor["sha256"])

    def test_repaired_bytes_are_bound(self) -> None:
        for binding in self.record["production_repairs"]:
            with self.subTest(path=binding["path"]):
                self.assertEqual(_sha256(binding["path"]), binding["sha256"])
        sibling = self.record["inflow_sibling"]
        for key in [
            "input",
            "historical_characterization_harness",
            "repaired_successor_harness",
        ]:
            with self.subTest(path=sibling[key]["path"]):
                self.assertEqual(_sha256(sibling[key]["path"]), sibling[key]["sha256"])

    def test_debug_replay_records_advancing_repeated_amr(self) -> None:
        replay = self.record["debug_serial_replay"]
        self.assertEqual(replay["result"], "pass")
        self.assertEqual(replay["restart_continuations"], 5)
        self.assertTrue(replay["finite_mesh_state"])
        self.assertTrue(
            all(after > before for before, after in zip(replay["stage_times"][:-1],
                                                        replay["stage_times"][1:]))
        )
        self.assertEqual(replay["meshblocks"], [8, 2, 8, 2, 8, 2])
        self.assertEqual(len(set(replay["particle_counts"])), 1)
        self.assertTrue(all(count > 0 for count in replay["ownership_changes_per_transition"]))
        self.assertTrue(replay["reflection_observed"])

    def test_sanitizer_replay_is_explicitly_bounded(self) -> None:
        replay = self.record["fresh_host_gnu_asan_ubsan_replay"]
        self.assertEqual(replay["archived_reflect_outflow_result"], "pass")
        self.assertEqual(replay["repaired_inflow_result"], "pass")
        self.assertEqual(replay["leaksanitizer_reported_errors"], 0)
        self.assertEqual(replay["ubsan_reported_errors"], 0)
        self.assertEqual(
            self.record["qualification_effect"],
            "bounded_local_repair_evidence_only",
        )
        self.assertTrue(any("MPI" in gap for gap in self.record["remaining_q009_gaps"]))
        self.assertTrue(any("HIP" in gap for gap in self.record["remaining_q009_gaps"]))


if __name__ == "__main__":
    unittest.main()
