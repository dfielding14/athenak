#!/usr/bin/env python3
"""Regression tests for the fail-closed PIC qualification-manifest gate."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.publication.pic_qualification_manifest import (
    _source_bundle_sha256,
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
        nested = self.root / "nested-repo"
        source = self.root / "source-repo"
        for repository in (nested, source):
            subprocess.run(["git", "init", str(repository)], check=True, capture_output=True)
            (repository / "tracked.txt").write_text(f"{repository.name}\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(repository), "add", "tracked.txt"], check=True)
            subprocess.run(
                [
                    "git", "-C", str(repository), "-c", "user.name=PIC Test",
                    "-c", "user.email=pic-test@example.invalid", "commit", "-m", "fixture",
                ],
                check=True, capture_output=True,
            )
        subprocess.run(
            [
                "git", "-c", "protocol.file.allow=always", "-C", str(source),
                "submodule", "add", str(nested), "kokkos",
            ],
            check=True, capture_output=True,
        )
        subprocess.run(["git", "-C", str(source), "add", "."], check=True)
        subprocess.run(
            [
                "git", "-C", str(source), "-c", "user.name=PIC Test",
                "-c", "user.email=pic-test@example.invalid", "commit", "-m", "add submodule",
            ],
            check=True, capture_output=True,
        )
        subprocess.run(
            ["git", "-C", str(source), "archive", "--format=tar",
             f"--output={self.root / 'source.tar'}", "HEAD"],
            check=True,
        )
        subprocess.run(
            ["git", "-C", str(source / "kokkos"), "archive", "--format=tar",
             f"--output={self.root / 'kokkos.tar'}", "HEAD"],
            check=True,
        )
        commit = subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
        ).strip()
        tree = subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "HEAD^{tree}"], text=True
        ).strip()
        kokkos_commit = subprocess.check_output(
            ["git", "-C", str(source / "kokkos"), "rev-parse", "HEAD"], text=True
        ).strip()
        kokkos_tree = subprocess.check_output(
            ["git", "-C", str(source / "kokkos"), "rev-parse", "HEAD^{tree}"], text=True
        ).strip()
        profile_submodules = [
            {
                "path": "kokkos",
                "archive_sha256": _sha256(self.root / "kokkos.tar"),
                "git_commit": kokkos_commit,
                "git_tree": kokkos_tree,
            }
        ]
        source_bundle = _source_bundle_sha256(
            _sha256(self.root / "source.tar"), profile_submodules
        )
        candidate = {
            "schema_version": 2,
            "freeze_id": "03a7bd9a-7d4c-4e37-a12b-46de3817eff2",
            "created_utc": "2026-05-30T12:00:00Z",
            "source": {
                "archive_path": "/original/source.tar",
                "archive_sha256": _sha256(self.root / "source.tar"),
                "source_bundle_sha256": source_bundle,
                "git_commit": commit,
                "git_tree": tree,
                "worktree_status": "clean",
                "submodule_status": "clean_pinned_archived",
                "submodules": [
                    {
                        **profile_submodules[0],
                        "archive_path": "/original/submodules/0000.tar",
                        "worktree_status": "clean",
                    }
                ],
            },
            "build": {
                "profile_id": "fixture",
                "profile_path": "/original/build_profile.json",
                "profile_sha256": "4" * 64,
                "source_archive_sha256": _sha256(self.root / "source.tar"),
                "source_bundle_sha256": source_bundle,
                "toolchain": "fixture",
                "build_command": "fixture",
                "executable_path": "/original/athena",
                "executable_sha256": _sha256(self.root / "athena"),
            },
        }
        (self.root / "clean_candidate_manifest.json").write_text(
            json.dumps(candidate), encoding="utf-8"
        )
        self.manifest = {
            "schema_version": 1,
            "manifest_id": "qualification-fixture-001",
            "created_utc": "2026-05-30T12:00:00Z",
            "claim_ids": ["CLAIM-PAPER-GYRO-001"],
            "test_id": "pic_relativistic_gyro_paper",
            "evidence_class": "sun_bai_2023_reproduction",
            "physical_mode": "paper_test_particle",
            "git": {
                "commit": commit,
                "tree": tree,
                "status": [],
                "source_archive": {
                    "path": "source.tar",
                    "sha256": _sha256(self.root / "source.tar"),
                },
                "source_bundle_sha256": "",
                "submodule_status": "clean_pinned_archived",
                "submodules": [
                    {
                        "path": "kokkos",
                        "archive_path": "kokkos.tar",
                        "archive_sha256": _sha256(self.root / "kokkos.tar"),
                        "git_commit": kokkos_commit,
                        "git_tree": kokkos_tree,
                        "worktree_status": "clean",
                    }
                ],
                "clean_candidate_manifest": {
                    "path": "clean_candidate_manifest.json",
                    "sha256": _sha256(self.root / "clean_candidate_manifest.json"),
                },
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
        self.manifest["git"]["source_bundle_sha256"] = _source_bundle_sha256(
            self.manifest["git"]["source_archive"]["sha256"],
            self.manifest["git"]["submodules"],
        )

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

    def test_submodule_bundle_and_clean_candidate_drift_are_rejected(self) -> None:
        dirty = copy.deepcopy(self.manifest)
        dirty["git"]["submodules"][0]["worktree_status"] = "dirty"
        self._assert_rejected(dirty)

        mismatched = copy.deepcopy(self.manifest)
        mismatched["git"]["source_bundle_sha256"] = "0" * 64
        self._assert_rejected(mismatched)

        candidate = copy.deepcopy(self.manifest)
        candidate["git"]["clean_candidate_manifest"]["sha256"] = "0" * 64
        self._assert_rejected(candidate)

        projected = copy.deepcopy(self.manifest)
        projected["git"]["tree"] = "0" * 40
        self._assert_rejected(projected)

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

        aliased = copy.deepcopy(self.manifest)
        aliased["git"]["source_archive"]["path"] = "./source.tar"
        self._assert_rejected(aliased)

        (self.root / "metrics-link.json").symlink_to("metrics.json")
        symlinked = copy.deepcopy(self.manifest)
        symlinked["artifacts"][0]["path"] = "metrics-link.json"
        self._assert_rejected(symlinked)

    def test_freeze_writes_canonical_file_once(self) -> None:
        source = self.root / "prepared.json"
        output = self.root / "frozen.json"
        source.write_text(json.dumps(self.manifest), encoding="utf-8")
        freeze_qualification_manifest(source, output)
        self.assertEqual(json.loads(output.read_text(encoding="utf-8")),
                         self.manifest)
        self.assertFalse(output.stat().st_mode & 0o222)
        with self.assertRaises(FileExistsError):
            freeze_qualification_manifest(source, output)


if __name__ == "__main__":
    unittest.main()
