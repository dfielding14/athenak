#!/usr/bin/env python3
"""Focused tests for immutable Q-011 Section 5.4 derived artifacts."""

from __future__ import annotations

import json
import os
from pathlib import Path
import stat
import tempfile
import unittest
from unittest import mock

from . import q011_section54_artifacts as artifacts


_RAW_DIGEST = "a" * 64
_ANALYZER_DIGEST = "b" * 64


def _pending() -> dict[str, object]:
    return {
        "status": "pending_external_review",
        "reviewer": None,
        "reviewed_at_utc": None,
        "notes": "Named external review remains explicitly deferred.",
    }


def _manifest() -> dict[str, object]:
    return artifacts.build_derived_manifest(
        campaign_id="q011-section54-fixture",
        raw_artifact_inventory_sha256=_RAW_DIGEST,
        analyzer_bindings={"tst/publication/analyze_fixture.py": _ANALYZER_DIGEST},
        reviewer_disposition=_pending(),
    )


def _make_writable(root: Path) -> None:
    for directory, names, filenames in os.walk(root):
        os.chmod(directory, 0o755)
        for name in names:
            os.chmod(Path(directory) / name, 0o755)
        for name in filenames:
            os.chmod(Path(directory) / name, 0o644)


class Q011Section54DerivedArtifactTests(unittest.TestCase):
    def test_publish_verify_and_refuse_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            receipt = artifacts.publish_derived_bundle(
                root,
                manifest=_manifest(),
                artifacts={
                    "metrics/summary.json": artifacts.canonical_json_bytes({"pass": True}),
                    "figures/morphology.png": b"png-fixture",
                },
            )
            self.assertEqual(
                artifacts.verify_published_derived_bundle(
                    root, expected_inventory_sha256=receipt["inventory_sha256"]
                ),
                _manifest(),
            )
            self.assertFalse(root.stat().st_mode & 0o222)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "already exists"):
                artifacts.publish_derived_bundle(root, manifest=_manifest(), artifacts={})
            _make_writable(root)

    def test_verify_rejects_tamper_and_tree_addition(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            artifacts.publish_derived_bundle(
                root, manifest=_manifest(), artifacts={"metrics/result.json": b"{}\n"}
            )
            _make_writable(root)
            (root / "metrics/result.json").write_bytes(b'{"drift": true}\n')
            for path in root.rglob("*"):
                os.chmod(path, 0o444 if path.is_file() else 0o555)
            os.chmod(root, 0o555)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "checksum drifted"):
                artifacts.verify_published_derived_bundle(root)
            _make_writable(root)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            artifacts.publish_derived_bundle(root, manifest=_manifest(), artifacts={})
            _make_writable(root)
            (root / "undeclared.txt").write_text("unexpected", encoding="utf-8")
            for path in root.rglob("*"):
                os.chmod(path, 0o444 if path.is_file() else 0o555)
            os.chmod(root, 0o555)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "tree closure"):
                artifacts.verify_published_derived_bundle(root)
            _make_writable(root)

    def test_publication_rejects_unsafe_reserved_and_nonfinite_members(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            for name in ("../escape", "/absolute", artifacts.INVENTORY_NAME):
                with self.subTest(name=name), self.assertRaises(artifacts.DerivedArtifactError):
                    artifacts.publish_derived_bundle(
                        Path(directory) / name.replace("/", "_"),
                        manifest=_manifest(),
                        artifacts={name: b"payload"},
                    )
        with self.assertRaises(artifacts.DerivedArtifactError):
            artifacts.canonical_json_bytes({"not_finite": float("nan")})

    def test_reviewer_disposition_requires_named_terminal_reviewer(self) -> None:
        accepted = {
            "status": "accepted",
            "reviewer": "Named Reviewer",
            "reviewed_at_utc": "2026-06-02T00:00:00Z",
            "notes": "Reviewed.",
        }
        self.assertEqual(artifacts.validate_reviewer_disposition(accepted), accepted)
        for disposition in (
            {**accepted, "reviewer": ""},
            {**accepted, "reviewed_at_utc": None},
            {**_pending(), "reviewer": "Premature Name"},
            {**_pending(), "status": "unknown"},
        ):
            with self.subTest(disposition=disposition), self.assertRaises(
                artifacts.DerivedArtifactError
            ):
                artifacts.validate_reviewer_disposition(disposition)

    def test_verify_rejects_duplicate_inventory_keys(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            artifacts.publish_derived_bundle(root, manifest=_manifest(), artifacts={})
            _make_writable(root)
            inventory = root / artifacts.INVENTORY_NAME
            inventory.write_text(
                '{"record_type":"q011_section54_derived_artifact_inventory",'
                '"schema_version":1,"schema_version":1,"members":[]}',
                encoding="utf-8",
            )
            for path in root.rglob("*"):
                os.chmod(path, 0o444 if path.is_file() else 0o555)
            os.chmod(root, 0o555)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "duplicate JSON key"):
                artifacts.verify_published_derived_bundle(root)
            _make_writable(root)

    def test_numeric_schema_aliases_fail_closed(self) -> None:
        with self.assertRaisesRegex(artifacts.DerivedArtifactError, "identity drifted"):
            artifacts.validate_derived_manifest({**_manifest(), "schema_version": True})

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            artifacts.publish_derived_bundle(root, manifest=_manifest(), artifacts={})
            _make_writable(root)
            inventory = root / artifacts.INVENTORY_NAME
            value = json.loads(inventory.read_text(encoding="utf-8"))
            value["schema_version"] = True
            inventory.write_bytes(artifacts.canonical_json_bytes(value))
            for path in root.rglob("*"):
                os.chmod(path, 0o444 if path.is_file() else 0o555)
            os.chmod(root, 0o555)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "schema drifted"):
                artifacts.verify_published_derived_bundle(root)
            _make_writable(root)

    def test_atomic_publish_refuses_competing_destination(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            root = parent / "bundle"
            original_rename = artifacts._rename_no_replace_at

            def collide(
                source_parent_fd: int,
                source_name: str,
                destination_parent_fd: int,
                destination_name: str,
            ) -> None:
                root.mkdir()
                (root / "sentinel.txt").write_text("competing\n", encoding="utf-8")
                original_rename(
                    source_parent_fd, source_name, destination_parent_fd, destination_name
                )

            with mock.patch.object(
                artifacts, "_rename_no_replace_at", side_effect=collide
            ):
                with self.assertRaisesRegex(artifacts.DerivedArtifactError, "already exists"):
                    artifacts.publish_derived_bundle(root, manifest=_manifest(), artifacts={})
            self.assertEqual((root / "sentinel.txt").read_text(encoding="utf-8"), "competing\n")

    def test_private_container_substitution_never_exposes_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            root = parent / "bundle"
            original_rename = artifacts._rename_no_replace_at
            substituted = False
            parked: list[Path] = []
            observed_public_replacement: list[bool] = []
            safe_public_inodes: list[int] = []
            observed_public_inodes: list[int] = []
            observed_public_modes: list[int] = []

            def substitute_then_rename(
                source_parent_fd: int,
                source_name: str,
                destination_parent_fd: int,
                destination_name: str,
            ) -> None:
                nonlocal substituted
                if not substituted:
                    private = next(
                        path
                        for path in parent.iterdir()
                        if path.name.startswith(".bundle.staging-")
                    )
                    parked_root = private.with_name(f"{private.name}.parked")
                    private.rename(parked_root)
                    safe_public_inodes.append(
                        (parked_root / artifacts._STAGING_ROOT_NAME).stat().st_ino
                    )
                    private.mkdir()
                    (private / artifacts._STAGING_ROOT_NAME).mkdir()
                    (private / artifacts._STAGING_ROOT_NAME / "sentinel.txt").write_text(
                        "attacker replacement\n", encoding="utf-8"
                    )
                    parked.append(parked_root)
                    substituted = True
                original_rename(
                    source_parent_fd, source_name, destination_parent_fd, destination_name
                )
                published = root.stat()
                observed_public_inodes.append(published.st_ino)
                observed_public_modes.append(stat.S_IMODE(published.st_mode))
                try:
                    replacement_is_public = (root / "sentinel.txt").exists()
                except PermissionError:
                    replacement_is_public = False
                observed_public_replacement.append(replacement_is_public)

            with mock.patch.object(
                artifacts, "_rename_no_replace_at", side_effect=substitute_then_rename
            ):
                receipt = artifacts.publish_derived_bundle(
                    root, manifest=_manifest(), artifacts={}
                )
            self.assertEqual(receipt["path"], str(root))
            self.assertTrue(root.is_dir())
            self.assertEqual(observed_public_replacement, [False])
            self.assertEqual(observed_public_inodes, safe_public_inodes)
            self.assertEqual(observed_public_modes, [stat.S_IWUSR])
            self.assertTrue(parked)
            self.assertFalse((root / "sentinel.txt").exists())

    def test_rollback_quarantines_moved_original_and_preserves_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            destination = parent / "bundle"
            destination.mkdir()
            (destination / "owned.txt").write_text("owned\n", encoding="utf-8")
            moved = parent / "moved-public-bundle"
            parent_fd = os.open(parent, artifacts._DIRECTORY_FLAGS)
            destination_fd = os.open(
                destination.name, artifacts._DIRECTORY_FLAGS, dir_fd=parent_fd
            )
            original_rename = artifacts._rename_no_replace_at
            raced = False

            def race_then_rename(*args: object, **kwargs: object) -> None:
                nonlocal raced
                if not raced:
                    destination.rename(moved)
                    destination.mkdir()
                    (destination / "sentinel.txt").write_text(
                        "replacement\n", encoding="utf-8"
                    )
                    raced = True
                original_rename(*args, **kwargs)

            try:
                with mock.patch.object(
                    artifacts, "_rename_no_replace_at", side_effect=race_then_rename
                ), self.assertRaisesRegex(
                    artifacts.DerivedArtifactError, "substituted public destination"
                ):
                    artifacts._rollback_published_destination(
                        parent_fd, destination.name, destination_fd
                    )
                self.assertEqual(
                    (destination / "sentinel.txt").read_text(encoding="utf-8"),
                    "replacement\n",
                )
                self.assertFalse(moved.exists())
                self.assertFalse(
                    any(
                        path.name.startswith(f".{destination.name}.rollback-")
                        for path in parent.iterdir()
                    )
                )
            finally:
                os.close(destination_fd)
                os.close(parent_fd)

    def test_verify_rejects_read_only_empty_directory_injection(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            artifacts.publish_derived_bundle(root, manifest=_manifest(), artifacts={})
            _make_writable(root)
            injected = root / "injected-empty"
            injected.mkdir()
            for path in root.rglob("*"):
                os.chmod(path, 0o444 if path.is_file() else 0o555)
            root.chmod(0o555)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "directory closure"):
                artifacts.verify_published_derived_bundle(root)
            _make_writable(root)

    def test_verify_rejects_hardlinked_member(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "bundle"
            artifacts.publish_derived_bundle(
                root, manifest=_manifest(), artifacts={"metrics/result.json": b"{}\n"}
            )
            _make_writable(root)
            os.link(root / "metrics/result.json", root / "metrics/result-alias.json")
            for path in root.rglob("*"):
                os.chmod(path, 0o444 if path.is_file() else 0o555)
            root.chmod(0o555)
            with self.assertRaisesRegex(artifacts.DerivedArtifactError, "hard-linked"):
                artifacts.verify_published_derived_bundle(root)
            _make_writable(root)

    def test_nested_directory_entries_are_fsynced_during_write(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            descriptor = os.open(root, artifacts._DIRECTORY_FLAGS)
            synced_inodes: list[int] = []
            original_fsync = os.fsync

            def record_fsync(fd: int) -> None:
                status = os.fstat(fd)
                if stat.S_ISDIR(status.st_mode):
                    synced_inodes.append(status.st_ino)
                original_fsync(fd)

            try:
                with mock.patch.object(artifacts.os, "fsync", side_effect=record_fsync):
                    artifacts._write_exclusive_at(
                        descriptor, "nested/deep/result.json", b"{}\n"
                    )
                self.assertIn((root / "nested").stat().st_ino, synced_inodes)
                self.assertIn((root / "nested/deep").stat().st_ino, synced_inodes)
            finally:
                os.close(descriptor)


if __name__ == "__main__":
    unittest.main()
