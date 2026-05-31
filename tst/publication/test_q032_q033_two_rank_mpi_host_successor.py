#!/usr/bin/env python3
"""Regression for rejected Q-032/Q-033 login-node parallel chronology."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q032_q033_two_rank_mpi_host_successor_2026-05-30.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _inventory_sha256(root: Path) -> str:
    rows = [
        f"{_sha256(path)}  {path.relative_to(root).as_posix()}\n"
        for path in sorted(path for path in root.rglob("*") if path.is_file())
    ]
    return hashlib.sha256("".join(rows).encode("utf-8")).hexdigest()


class Q032Q033TwoRankMPIHostSuccessorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        cls.artifacts = cls.sidecar["retained_artifacts"]
        cls.root = Path(cls.artifacts["artifact_root"])

    def test_record_remains_narrow_and_nonqualifying(self) -> None:
        self.assertEqual(
            self.sidecar["qualification_effect"],
            "none_rejected_login_node_parallel_launch_chronology_only",
        )
        self.assertFalse(self.sidecar["claim_closure"])
        self.assertEqual(self.sidecar["frontier_authorization"], "not_bound")
        self.assertIn("one Frontier login host", self.sidecar["scope"])
        self.assertIn("not readiness evidence", self.sidecar["scope"].lower())
        self.assertTrue(self.sidecar["explicitly_not_claimed"])
        rejection = self.sidecar["policy_rejection"]
        self.assertEqual(rejection["status"], "rejected_not_readiness_evidence")
        self.assertIn("frontier_user_guide.html", rejection["source"])
        self.assertIn("never launch parallel jobs", rejection["rule"])

    def test_source_bindings_match_current_bytes(self) -> None:
        for relative, expected in self.sidecar["source_bindings"].items():
            with self.subTest(relative=relative):
                self.assertEqual(_sha256(REPO_ROOT / relative), expected)

    def test_retained_tree_is_read_only_and_inventory_bound(self) -> None:
        files = [path for path in self.root.rglob("*") if path.is_file()]
        writable = [
            path
            for path in [self.root, *self.root.rglob("*")]
            if path.stat().st_mode & 0o222
        ]
        self.assertEqual(len(files), self.artifacts["file_count"])
        self.assertEqual(len(writable), self.artifacts["writable_entries"])
        self.assertEqual(_inventory_sha256(self.root), self.artifacts["inventory_sha256"])
        for path_key, sha_key in (
            ("retained_manifest_path", "retained_manifest_sha256"),
            ("executable_path", "executable_sha256"),
            ("cmake_cache_path", "cmake_cache_sha256"),
            ("q032_log_path", "q032_log_sha256"),
            ("q033_log_path", "q033_log_sha256"),
        ):
            self.assertEqual(
                _sha256(self.root / self.artifacts[path_key]),
                self.artifacts[sha_key],
            )

    def test_launcher_and_build_profile_are_rejected_login_host_chronology(self) -> None:
        launcher = self.sidecar["launcher"]
        self.assertEqual(launcher["hostname_probe"], ["login10", "login10"])
        self.assertEqual(launcher["slurm_job_id"], "")
        self.assertEqual(launcher["slurm_step_id"], "")
        self.assertFalse(launcher["scheduler_launch_used"])
        build = self.sidecar["build_profile"]
        self.assertEqual(build["Athena_ENABLE_MPI"], "ON")
        self.assertEqual(build["Kokkos_ENABLE_MPI"], "ON")
        self.assertEqual(build["Kokkos_ENABLE_HIP"], "OFF")

    def test_q032_replay_preserves_bounded_mechanics(self) -> None:
        replay = self.sidecar["q032_two_rank_replay"]
        self.assertEqual(replay["configured_ranks"], 2)
        self.assertEqual(replay["decomposition_override"], "meshblock/nx1=16")
        self.assertEqual(replay["qualification_effect"], "none")
        self.assertEqual(replay["plotnikov_qualification"], "not_claimed")
        self.assertLessEqual(max(replay["relative_errors"].values()), 1.0e-7)
        self.assertEqual(replay["result"], "pass")

    def test_q033_replay_preserves_exact_restart_endpoint(self) -> None:
        replay = self.sidecar["q033_two_rank_replay"]
        self.assertEqual(replay["configured_ranks"], 2)
        self.assertEqual(replay["decomposition_override"], "mesh/nx1=8")
        self.assertEqual(replay["checkpoint_continuations"], 1)
        self.assertEqual(replay["particle_count"], 128)
        self.assertEqual(replay["mhd_time_absolute_error"], 0.0)
        self.assertEqual(replay["mhd_field_absolute_error_max"], 0.0)
        self.assertEqual(replay["history_endpoint_absolute_error_max"], 0.0)
        self.assertTrue(replay["particle_integer_payload_equal"])
        self.assertEqual(replay["particle_float_payload_absolute_error_max"], 0.0)
        self.assertFalse(replay["restart_refit_observed"])


if __name__ == "__main__":
    unittest.main()
