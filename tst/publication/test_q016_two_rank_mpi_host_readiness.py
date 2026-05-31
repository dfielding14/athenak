"""Q-016 rejected login-node parallel chronology regression."""

from __future__ import annotations

import hashlib
import json
import re
import unittest
from pathlib import Path


_REPO_ROOT = Path(__file__).resolve().parents[2]
_RECORD_PATH = (
    _REPO_ROOT
    / "tst"
    / "publication"
    / "readiness"
    / "q016_two_rank_mpi_host_provenance_spectrum_replay_2026-05-30.json"
)
_SHA256 = re.compile(r"[0-9a-f]{64}")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class TestQ016TwoRankMPIHostReadiness(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = json.loads(_RECORD_PATH.read_text(encoding="utf-8"))

    def test_record_stays_narrow_and_fail_closed(self) -> None:
        record = self.record
        self.assertEqual(
            record["qualification_effect"],
            "none_rejected_login_node_parallel_launch_chronology_only",
        )
        self.assertIn("not readiness evidence", record["scope"].lower())
        rejection = record["policy_rejection"]
        self.assertEqual(rejection["status"], "rejected_not_readiness_evidence")
        self.assertIn("frontier_user_guide.html", rejection["source"])
        self.assertIn("never launch parallel jobs", rejection["rule"])
        self.assertTrue(
            (_REPO_ROOT / rejection["valid_compute_node_successor"]).is_file()
        )
        gaps = record["remaining_q016_gaps"]
        self.assertTrue(any("Cross-node MPI" in gap for gap in gaps))
        self.assertTrue(any("Frontier MPI and HIP" in gap for gap in gaps))
        self.assertTrue(any("Repeated-AMR" in gap for gap in gaps))
        self.assertTrue(any("shock-campaign" in gap for gap in gaps))

    def test_source_bindings_match_current_bytes(self) -> None:
        for binding in self.record["source_bindings"]:
            with self.subTest(path=binding["path"]):
                self.assertRegex(binding["sha256"], _SHA256)
                self.assertEqual(
                    binding["sha256"],
                    _sha256(_REPO_ROOT / binding["path"]),
                )

    def test_reused_build_and_local_launcher_are_explicit(self) -> None:
        build = self.record["reused_mpi_host_build"]
        self.assertEqual(build["configuration_query_result"]["mpi_parallelism"], "ON")
        self.assertEqual(build["cmake_cache_contract"]["Athena_ENABLE_MPI"], "ON")
        self.assertEqual(build["cmake_cache_contract"]["Kokkos_ENABLE_MPI"], "ON")
        probes = self.record["launcher_probes"]
        self.assertEqual(probes[0]["status"], "blocked_before_launch")
        self.assertEqual(probes[1]["status"], "pass")
        self.assertEqual(probes[1]["stdout_lines"], ["login10", "login10"])

    def test_replay_metrics_bind_two_rank_restart_and_spectra(self) -> None:
        replay = self.record["replay"]
        self.assertEqual(replay["environment"]["ATHENA_Q016_NPROC"], "2")
        self.assertEqual(replay["environment"]["SLURM_JOB_ID"], "")
        self.assertEqual(len(replay["observed_athena_commands"]), 3)
        self.assertEqual(replay["status"], "pass")
        metrics = replay["metrics"]
        self.assertEqual(metrics["configured_ranks"], 2)
        self.assertEqual(metrics["restart_schema"], 7)
        self.assertTrue(all(metrics["int_metadata_equal_after_restart"].values()))
        self.assertEqual(max(metrics["max_float_errors_after_restart"].values()), 0.0)
        self.assertTrue(metrics["migration_observed"])
        self.assertTrue(metrics["source_metadata_ok"])
        self.assertEqual(metrics["full_f_spectrum_agreement"], 0.0)
        self.assertEqual(metrics["delta_f_spectrum_agreement"], 0.0)

    def test_ephemeral_artifact_bindings_are_checksums_not_closure(self) -> None:
        evidence = self.record["ephemeral_evidence_bindings"]
        self.assertIn("not archived qualification", evidence["retention"])
        self.assertRegex(evidence["log"]["sha256"], _SHA256)
        artifacts = evidence["selected_terminal_artifacts"]
        self.assertEqual(len(artifacts), 4)
        for artifact in artifacts:
            self.assertRegex(artifact["sha256"], _SHA256)
        self.assertEqual(artifacts[0]["sha256"], artifacts[1]["sha256"])


if __name__ == "__main__":
    unittest.main()
