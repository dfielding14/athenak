#!/usr/bin/env python3
"""Validate the bounded Orion-local two-rank Q016 replay record."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
RECORD = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q016_particle_provenance_spectra_mpi2_orion_local_2026-05-30.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q016OrionLocalMPIReadinessTests(unittest.TestCase):
    def setUp(self) -> None:
        self.record = json.loads(RECORD.read_text(encoding="utf-8"))

    def test_source_local_bindings_match_current_bytes(self) -> None:
        expected = {
            "inputs/tests/pic_q016_particle_provenance.athinput",
            "tst/scripts/particles/pic_q016_particle_provenance.py",
            "tst/publication/pvtk_particles.py",
            "tst/publication/q016_particle_spectra.py",
            (
                "tst/publication/readiness/"
                "q016_particle_provenance_spectra_local_2026-05-30.json"
            ),
        }
        self.assertEqual(set(self.record["source_local_bindings"]), expected)
        for relative, digest in self.record["source_local_bindings"].items():
            with self.subTest(path=relative):
                self.assertEqual(_sha256(REPO_ROOT / relative), digest)

    def test_orion_shared_executable_binding_matches_retained_bytes(self) -> None:
        build = self.record["orion_shared_build"]
        self.assertTrue(build["path"].startswith(
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/"
        ))
        self.assertEqual(_sha256(Path(build["executable"])), build["executable_sha256"])

    def test_direct_scheduler_chronology_is_explicit_and_nonqualifying(self) -> None:
        chronology = self.record["manual_scheduler_chronology"]
        jobs = [
            chronology["failed_tmp_mount_probe"],
            *chronology["shared_root_dry_sequence"],
            *self.record["validated_replay"]["jobs"],
        ]
        node_hours = sum(job["elapsed_seconds"] * job["allocated_nodes"] for job in jobs)
        self.assertAlmostEqual(
            chronology["all_manual_q016_allocated_node_hours"], node_hours / 3600.0
        )
        self.assertIn("not inserted", chronology["accounting_boundary"])
        self.assertEqual(
            self.record["qualification_effect"],
            "bounded_controlled_orion_host_mpi_decomposition_evidence_only",
        )


if __name__ == "__main__":
    unittest.main()
