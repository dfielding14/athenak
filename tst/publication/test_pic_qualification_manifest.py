#!/usr/bin/env python3
"""Regression tests for the fail-closed PIC qualification-manifest gate."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.publication.pic_qualification_manifest import (
    freeze_qualification_manifest,
    validate_qualification_manifest,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PicQualificationManifestTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        for name, contents in (
            ("athena", "executable"),
            ("CMakeCache.txt", "cache"),
            ("modules.txt", "modules"),
            ("environment.txt", "environment"),
            ("metrics.json", '{"relative_error": 0.0}\n'),
        ):
            (self.root / name).write_text(contents, encoding="utf-8")
        self.manifest = {
            "schema_version": 1,
            "manifest_id": "qualification-fixture-001",
            "created_utc": "2026-05-30T12:00:00Z",
            "claim_ids": ["CLAIM-PAPER-GYRO-001"],
            "test_id": "pic_relativistic_gyro_paper",
            "evidence_class": "sun_bai_2023_reproduction",
            "physical_mode": "paper_test_particle",
            "git": {
                "commit": "0" * 40,
                "status": [],
                "submodules": [],
            },
            "executable": {
                "path": "athena",
                "sha256": _sha256(self.root / "athena"),
                "cmake_cache": "CMakeCache.txt",
                "modules": "modules.txt",
                "environment_allowlist": "environment.txt",
            },
            "parameters": {},
            "oracle": {
                "kind": "analytic",
                "reference": "bounded fixture",
                "tolerances": {"relative_error": 1.0e-6},
            },
            "metrics": [{"name": "relative_error", "value": 0.0}],
            "resources": {
                "platform": "host",
                "artifact_root": str(self.root),
            },
            "artifacts": [
                {
                    "path": "metrics.json",
                    "sha256": _sha256(self.root / "metrics.json"),
                }
            ],
            "review": {
                "reviewer": "pending external review",
                "disposition": "pending external review",
            },
        }

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def _assert_rejected(self, manifest: dict[str, object]) -> None:
        with self.assertRaises(ValueError):
            validate_qualification_manifest(manifest)

    def test_reviewable_qualification_manifest_is_accepted(self) -> None:
        validate_qualification_manifest(self.manifest)

    def test_proxy_and_unit_evidence_are_rejected(self) -> None:
        for evidence_class in ("engineering_proxy", "unit/regression"):
            with self.subTest(evidence_class=evidence_class):
                manifest = copy.deepcopy(self.manifest)
                manifest["evidence_class"] = evidence_class
                self._assert_rejected(manifest)

    def test_dirty_source_and_unknown_claim_are_rejected(self) -> None:
        dirty = copy.deepcopy(self.manifest)
        dirty["git"]["status"] = [" M src/particles/particles.cpp"]
        self._assert_rejected(dirty)

        unknown = copy.deepcopy(self.manifest)
        unknown["claim_ids"] = ["CLAIM-DOES-NOT-EXIST"]
        self._assert_rejected(unknown)

    def test_escaped_missing_and_checksum_mismatched_files_are_rejected(
        self,
    ) -> None:
        escaped = copy.deepcopy(self.manifest)
        escaped["artifacts"][0]["path"] = "../outside.json"
        self._assert_rejected(escaped)

        missing = copy.deepcopy(self.manifest)
        missing["executable"]["modules"] = "missing-modules.txt"
        self._assert_rejected(missing)

        mismatched = copy.deepcopy(self.manifest)
        mismatched["artifacts"][0]["sha256"] = "0" * 64
        self._assert_rejected(mismatched)

    def test_freeze_writes_canonical_file_once(self) -> None:
        source = self.root / "prepared.json"
        output = self.root / "frozen.json"
        source.write_text(json.dumps(self.manifest), encoding="utf-8")
        freeze_qualification_manifest(source, output)
        self.assertEqual(json.loads(output.read_text(encoding="utf-8")),
                         self.manifest)
        with self.assertRaises(FileExistsError):
            freeze_qualification_manifest(source, output)


if __name__ == "__main__":
    unittest.main()
