#!/usr/bin/env python3
"""Focused regression tests for immutable Orion tree publication."""

from __future__ import annotations

import os
from pathlib import Path
import stat
import tempfile
import unittest
from unittest import mock

from tst.publication import immutable_orion_tree


_RECEIPT = {
    "schema_version": 1,
    "artifact_role": "focused_test_evidence",
    "qualification_effect": "none",
    "inventory_excludes": immutable_orion_tree.INVENTORY_NAME,
    "freeze_policy": "remove all owner, group and other write bits recursively",
}
_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH


def _make_writable_tree(root: Path) -> None:
    if not root.exists():
        return
    for path in [root, *root.rglob("*")]:
        if not path.is_symlink():
            path.chmod(path.stat().st_mode | stat.S_IWUSR)


def _swap_sibling_directories(first: Path, second: Path) -> None:
    temporary = first.with_name("temporary-swap")
    os.rename(first, temporary)
    os.rename(second, first)
    os.rename(temporary, second)


class ImmutableOrionTreeTests(unittest.TestCase):
    def test_freeze_and_verify_use_recursive_read_only_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            tree = base / "retained"
            (tree / "nested").mkdir(parents=True)
            (tree / "root.txt").write_text("root\n", encoding="utf-8")
            (tree / "nested/payload.txt").write_text("nested\n", encoding="utf-8")
            try:
                report = immutable_orion_tree.freeze_tree(
                    tree,
                    _RECEIPT,
                    authorized_root=base,
                )
                verified = immutable_orion_tree.verify_frozen_tree(
                    tree,
                    report["inventory_sha256"],
                    authorized_root=base,
                )
                self.assertEqual(report["inventoried_file_count"], 3)
                self.assertEqual(verified["inventory_sha256"], report["inventory_sha256"])
                self.assertEqual(verified["writable_entries"], [])
                self.assertTrue(verified["recursively_read_only"])
                for path in [tree, *tree.rglob("*")]:
                    self.assertFalse(path.stat().st_mode & _WRITE_BITS, path)
            finally:
                _make_writable_tree(tree)

    def test_freeze_rejects_same_content_directory_swap_during_inventory_publication(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            tree = base / "retained"
            nested = tree / "nested"
            nested.mkdir(parents=True)
            (nested / "payload.txt").write_text("same bytes\n", encoding="utf-8")
            replacement = tree / "replacement"
            replacement.mkdir()
            (replacement / "payload.txt").write_text("same bytes\n", encoding="utf-8")
            original_write = immutable_orion_tree._write_new_metadata

            def swap_while_publishing(
                root_fd: int, name: str, payload: str, **kwargs: object
            ) -> None:
                if name == immutable_orion_tree.INVENTORY_NAME:
                    _swap_sibling_directories(nested, replacement)
                original_write(root_fd, name, payload, **kwargs)

            try:
                with mock.patch.object(
                    immutable_orion_tree,
                    "_write_new_metadata",
                    side_effect=swap_while_publishing,
                ):
                    with self.assertRaisesRegex(ValueError, "before inventory publication"):
                        immutable_orion_tree.freeze_tree(
                            tree,
                            _RECEIPT,
                            authorized_root=base,
                        )
            finally:
                _make_writable_tree(tree)

    def test_verify_rejects_same_content_directory_swap_between_scans(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            tree = base / "retained"
            nested = tree / "nested"
            nested.mkdir(parents=True)
            (nested / "payload.txt").write_text("same bytes\n", encoding="utf-8")
            replacement = tree / "replacement"
            replacement.mkdir()
            (replacement / "payload.txt").write_text("same bytes\n", encoding="utf-8")
            try:
                report = immutable_orion_tree.freeze_tree(
                    tree,
                    _RECEIPT,
                    authorized_root=base,
                )
                original_scan = immutable_orion_tree._scan_anchored_tree
                scan_count = 0

                def swap_after_baseline_scan(root_fd: int, **kwargs: object):
                    nonlocal scan_count
                    snapshot = original_scan(root_fd, **kwargs)
                    scan_count += 1
                    if scan_count == 1:
                        tree.chmod(tree.stat().st_mode | stat.S_IWUSR)
                        _swap_sibling_directories(nested, replacement)
                        tree.chmod(tree.stat().st_mode & ~_WRITE_BITS)
                    return snapshot

                with mock.patch.object(
                    immutable_orion_tree,
                    "_scan_anchored_tree",
                    side_effect=swap_after_baseline_scan,
                ):
                    with self.assertRaisesRegex(ValueError, "verification hash pass"):
                        immutable_orion_tree.verify_frozen_tree(
                            tree,
                            report["inventory_sha256"],
                            authorized_root=base,
                        )
            finally:
                _make_writable_tree(tree)

    def test_staged_snapshot_rejects_payload_mutation_after_initial_verification(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            tree = base / "retained"
            tree.mkdir()
            payload = tree / "payload.txt"
            payload.write_text("original\n", encoding="utf-8")
            try:
                report = immutable_orion_tree.freeze_tree(
                    tree,
                    _RECEIPT,
                    authorized_root=base,
                )
                original_verify = immutable_orion_tree._verify_frozen_tree_anchored
                verify_count = 0

                def mutate_after_initial_verify(*args: object, **kwargs: object):
                    nonlocal verify_count
                    verified = original_verify(*args, **kwargs)
                    verify_count += 1
                    if verify_count == 1:
                        payload.chmod(payload.stat().st_mode | stat.S_IWUSR)
                        payload.write_text("replacement\n", encoding="utf-8")
                        payload.chmod(payload.stat().st_mode & ~_WRITE_BITS)
                    return verified

                with mock.patch.object(
                    immutable_orion_tree,
                    "_verify_frozen_tree_anchored",
                    side_effect=mutate_after_initial_verify,
                ):
                    with self.assertRaisesRegex(ValueError, "SHA-256 drifted"):
                        with immutable_orion_tree.staged_verified_frozen_tree(
                            tree,
                            report["inventory_sha256"],
                            authorized_root=base,
                        ):
                            self.fail("unsafe post-verification mutation was accepted")
            finally:
                _make_writable_tree(tree)


if __name__ == "__main__":
    unittest.main()
