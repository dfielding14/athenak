#!/usr/bin/env python3
"""Focused regression tests for immutable Orion tree publication."""

from __future__ import annotations

import fcntl
import hashlib
import os
from pathlib import Path
import stat
import tarfile
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
                    with self.assertRaisesRegex(
                        ValueError,
                        "(SHA-256|artifact hash) drifted|sealed snapshot topology capture",
                    ):
                        with immutable_orion_tree.staged_verified_frozen_tree(
                            tree,
                            report["inventory_sha256"],
                            authorized_root=base,
                        ):
                            self.fail("unsafe post-verification mutation was accepted")
            finally:
                _make_writable_tree(tree)

    def test_staged_snapshot_yields_sealed_payload_member(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            tree = base / "retained"
            tree.mkdir()
            (tree / "payload.txt").write_text("verified\n", encoding="utf-8")
            try:
                report = immutable_orion_tree.freeze_tree(
                    tree,
                    _RECEIPT,
                    authorized_root=base,
                )
                with immutable_orion_tree.staged_verified_frozen_tree(
                    tree,
                    report["inventory_sha256"],
                    authorized_root=base,
                ) as (_, snapshot):
                    sealed = snapshot.member_path("payload.txt")
                    self.assertTrue(immutable_orion_tree.is_sealed_snapshot_member(sealed))
                    self.assertEqual(sealed.read_text(encoding="utf-8"), "verified\n")
                    with self.assertRaises(OSError):
                        os.open(sealed, os.O_WRONLY)
                    self.assertEqual(sealed.read_text(encoding="utf-8"), "verified\n")
                    self.assertFalse(snapshot.has_file("../payload.txt"))
                    self.assertFalse(snapshot.has_directory("../retained"))
                    with self.assertRaisesRegex(ValueError, "unsafe sealed snapshot"):
                        snapshot.member_path("../payload.txt")
            finally:
                _make_writable_tree(tree)

    def test_staged_snapshot_rejects_directory_hidden_during_topology_capture(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            tree = base / "retained"
            nested = tree / "nested"
            nested.mkdir(parents=True)
            (tree / "payload.txt").write_text("verified\n", encoding="utf-8")
            try:
                report = immutable_orion_tree.freeze_tree(
                    tree,
                    _RECEIPT,
                    authorized_root=base,
                )
                original_scan = immutable_orion_tree._scan_anchored_tree
                scan_count = 0

                def hide_on_topology_scan(root_fd: int, **kwargs: object):
                    nonlocal scan_count
                    scan_count += 1
                    if scan_count != 4:
                        return original_scan(root_fd, **kwargs)
                    tree.chmod(tree.stat().st_mode | stat.S_IWUSR)
                    nested.rmdir()
                    try:
                        return original_scan(root_fd, **kwargs)
                    finally:
                        nested.mkdir()
                        nested.chmod(nested.stat().st_mode & ~_WRITE_BITS)
                        tree.chmod(tree.stat().st_mode & ~_WRITE_BITS)

                with mock.patch.object(
                    immutable_orion_tree,
                    "_scan_anchored_tree",
                    side_effect=hide_on_topology_scan,
                ):
                    with self.assertRaisesRegex(ValueError, "sealed snapshot topology capture"):
                        with immutable_orion_tree.staged_verified_frozen_tree(
                            tree,
                            report["inventory_sha256"],
                            authorized_root=base,
                        ):
                            self.fail("hidden staged directory was accepted")
            finally:
                _make_writable_tree(tree)

    def test_staged_snapshot_rejects_directory_injected_during_topology_capture(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            tree = base / "retained"
            tree.mkdir()
            (tree / "payload.txt").write_text("verified\n", encoding="utf-8")
            injected = tree / "injected"
            try:
                report = immutable_orion_tree.freeze_tree(
                    tree,
                    _RECEIPT,
                    authorized_root=base,
                )
                original_scan = immutable_orion_tree._scan_anchored_tree
                scan_count = 0

                def inject_on_topology_scan(root_fd: int, **kwargs: object):
                    nonlocal scan_count
                    scan_count += 1
                    if scan_count != 4:
                        return original_scan(root_fd, **kwargs)
                    tree.chmod(tree.stat().st_mode | stat.S_IWUSR)
                    injected.mkdir()
                    try:
                        return original_scan(root_fd, **kwargs)
                    finally:
                        injected.rmdir()
                        tree.chmod(tree.stat().st_mode & ~_WRITE_BITS)

                with mock.patch.object(
                    immutable_orion_tree,
                    "_scan_anchored_tree",
                    side_effect=inject_on_topology_scan,
                ):
                    with self.assertRaisesRegex(ValueError, "sealed snapshot topology capture"):
                        with immutable_orion_tree.staged_verified_frozen_tree(
                            tree,
                            report["inventory_sha256"],
                            authorized_root=base,
                        ):
                            self.fail("injected staged directory was accepted")
            finally:
                _make_writable_tree(tree)

    def test_staged_snapshot_rejects_file_hidden_during_topology_capture(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            tree = base / "retained"
            tree.mkdir()
            payload = tree / "payload.txt"
            payload.write_text("verified\n", encoding="utf-8")
            parked = base / "parked.txt"
            try:
                report = immutable_orion_tree.freeze_tree(
                    tree,
                    _RECEIPT,
                    authorized_root=base,
                )
                original_scan = immutable_orion_tree._scan_anchored_tree
                scan_count = 0

                def hide_on_topology_scan(root_fd: int, **kwargs: object):
                    nonlocal scan_count
                    scan_count += 1
                    if scan_count != 4:
                        return original_scan(root_fd, **kwargs)
                    tree.chmod(tree.stat().st_mode | stat.S_IWUSR)
                    os.rename(payload, parked)
                    try:
                        return original_scan(root_fd, **kwargs)
                    finally:
                        os.rename(parked, payload)
                        tree.chmod(tree.stat().st_mode & ~_WRITE_BITS)

                with mock.patch.object(
                    immutable_orion_tree,
                    "_scan_anchored_tree",
                    side_effect=hide_on_topology_scan,
                ):
                    with self.assertRaisesRegex(ValueError, "sealed snapshot topology capture"):
                        with immutable_orion_tree.staged_verified_frozen_tree(
                            tree,
                            report["inventory_sha256"],
                            authorized_root=base,
                        ):
                            self.fail("hidden staged file was accepted")
            finally:
                _make_writable_tree(tree)

    def test_q006_and_q007_reject_unrelated_sealed_memfd(self) -> None:
        from tst.publication import analyze_q006_paper_multispecies_oscillation_runtime_local
        from tst.publication import analyze_q007_paper_deltaf_linear_preparation

        fd = os.memfd_create("unrelated", flags=os.MFD_ALLOW_SEALING)
        try:
            os.write(fd, b"unrelated\n")
            fcntl.fcntl(fd, fcntl.F_ADD_SEALS, immutable_orion_tree._MEMFD_SEALS)
            sealed = Path("/proc/self/fd") / str(fd)
            with tempfile.TemporaryDirectory() as directory:
                retained_root = Path(directory)
                with self.assertRaisesRegex(ValueError, "remain below retained root"):
                    analyze_q006_paper_multispecies_oscillation_runtime_local._contained_regular_file(
                        retained_root,
                        sealed,
                    )
                with self.assertRaisesRegex(ValueError, "remain below retained root"):
                    analyze_q007_paper_deltaf_linear_preparation._contained_regular_file(
                        retained_root,
                        sealed,
                    )
        finally:
            os.close(fd)

    def test_q006_and_q007_reject_unknown_file_injected_after_snapshot_handoff(self) -> None:
        from tst.publication import analyze_q006_paper_multispecies_oscillation_runtime_local
        from tst.publication import analyze_q007_paper_deltaf_linear_preparation

        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            tree = base / "retained"
            tree.mkdir()
            (tree / "payload.txt").write_text("verified\n", encoding="utf-8")
            try:
                report = immutable_orion_tree.freeze_tree(
                    tree,
                    _RECEIPT,
                    authorized_root=base,
                )
                with immutable_orion_tree.staged_verified_frozen_tree(
                    tree,
                    report["inventory_sha256"],
                    authorized_root=base,
                ) as (_, snapshot):
                    injected = snapshot.staged_root / "injected.txt"
                    injected.write_text("attacker bytes\n", encoding="utf-8")
                    for analyzer in (
                        analyze_q006_paper_multispecies_oscillation_runtime_local,
                        analyze_q007_paper_deltaf_linear_preparation,
                    ):
                        with analyzer._use_staged_tree(tree, snapshot):
                            with self.assertRaisesRegex(ValueError, "snapshot member is absent"):
                                analyzer._contained_regular_file(tree, tree / "injected.txt")
            finally:
                _make_writable_tree(tree)

    def test_archive_and_elf_validators_inspect_sealed_copy_after_path_swap(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            source = base / "source.txt"
            source.write_text("source\n", encoding="utf-8")
            archive = base / "source.tar"
            with tarfile.open(archive, "w") as handle:
                handle.add(source, arcname="source.txt")
            archive_sha256 = hashlib.sha256(archive.read_bytes()).hexdigest()
            archive_replacement = base / "source-replacement.tar"
            archive_replacement.write_bytes(b"not a tar archive\n")
            archived_original = base / "source-original.tar"
            original_copy = immutable_orion_tree._sealed_regular_copy
            swapped = False

            def swap_archive_after_copy(*args: object, **kwargs: object):
                nonlocal swapped
                copied = original_copy(*args, **kwargs)
                if not swapped:
                    os.rename(archive, archived_original)
                    os.rename(archive_replacement, archive)
                    swapped = True
                return copied

            with mock.patch.object(
                immutable_orion_tree,
                "_sealed_regular_copy",
                side_effect=swap_archive_after_copy,
            ):
                report = immutable_orion_tree.validate_source_archive(
                    archive,
                    archive_sha256,
                )
            self.assertTrue(report["passed"])

            executable = base / "athena"
            executable.write_bytes(b"\x7fELFpayload")
            executable.chmod(0o555)
            executable_sha256 = hashlib.sha256(executable.read_bytes()).hexdigest()
            executable_replacement = base / "replacement-athena"
            executable_replacement.write_bytes(b"not-elf")
            executable_replacement.chmod(0o555)
            archived_executable = base / "original-athena"
            swapped = False

            def swap_executable_after_copy(*args: object, **kwargs: object):
                nonlocal swapped
                copied = original_copy(*args, **kwargs)
                if not swapped:
                    os.rename(executable, archived_executable)
                    os.rename(executable_replacement, executable)
                    swapped = True
                return copied

            with mock.patch.object(
                immutable_orion_tree,
                "_sealed_regular_copy",
                side_effect=swap_executable_after_copy,
            ):
                report = immutable_orion_tree.validate_executable_elf(
                    executable,
                    executable_sha256,
                )
            self.assertTrue(report["elf_identity"])

            dependency_archive = base / "dependency-source.tar"
            with tarfile.open(dependency_archive, "w") as handle:
                handle.add(source, arcname="source.txt")
            dependency_replacement = base / "dependency-replacement.tar"
            dependency_replacement.write_bytes(b"not a tar archive\n")
            archived_dependency = base / "dependency-original.tar"
            swapped = False

            def swap_dependency_after_copy(*args: object, **kwargs: object):
                nonlocal swapped
                copied = original_copy(*args, **kwargs)
                if not swapped:
                    os.rename(dependency_archive, archived_dependency)
                    os.rename(dependency_replacement, dependency_archive)
                    swapped = True
                return copied

            with mock.patch.object(
                immutable_orion_tree,
                "_sealed_regular_copy",
                side_effect=swap_dependency_after_copy,
            ):
                report = immutable_orion_tree.validate_source_archive_dependencies(
                    dependency_archive,
                    {"source.txt": hashlib.sha256(source.read_bytes()).hexdigest()},
                    {"source.txt"},
                )
            self.assertTrue(report["passed"])

    def test_multi_pass_validators_inspect_sealed_copy_after_inplace_rewrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            original_source = base / "original.txt"
            original_source.write_text("original\n", encoding="utf-8")
            replacement_source = base / "replacement.txt"
            replacement_source.write_text("replacement\n", encoding="utf-8")
            archive = base / "source.tar"
            with tarfile.open(archive, "w") as handle:
                handle.add(original_source, arcname="original.txt")
            archive_sha256 = hashlib.sha256(archive.read_bytes()).hexdigest()
            replacement_archive = base / "replacement.tar"
            with tarfile.open(replacement_archive, "w") as handle:
                handle.add(replacement_source, arcname="replacement.txt")
                handle.add(original_source, arcname="original.txt")
            executable = base / "athena"
            executable.write_bytes(b"\x7fELForiginal")
            executable.chmod(0o555)
            executable_sha256 = hashlib.sha256(executable.read_bytes()).hexdigest()
            original_copy = immutable_orion_tree._sealed_regular_copy
            rewritten = False

            def rewrite_archive_after_copy(*args: object, **kwargs: object):
                nonlocal rewritten
                copied = original_copy(*args, **kwargs)
                if not rewritten:
                    archive.chmod(0o644)
                    archive.write_bytes(replacement_archive.read_bytes())
                    rewritten = True
                return copied

            with mock.patch.object(
                immutable_orion_tree,
                "_sealed_regular_copy",
                side_effect=rewrite_archive_after_copy,
            ):
                report = immutable_orion_tree.validate_source_archive(
                    archive,
                    archive_sha256,
                )
            self.assertEqual(report["member_count"], 1)
            self.assertEqual(report["regular_file_count"], 1)

            rewritten = False

            def rewrite_executable_after_copy(*args: object, **kwargs: object):
                nonlocal rewritten
                copied = original_copy(*args, **kwargs)
                if not rewritten:
                    executable.chmod(0o755)
                    executable.write_bytes(b"not-elf-replacement")
                    rewritten = True
                return copied

            with mock.patch.object(
                immutable_orion_tree,
                "_sealed_regular_copy",
                side_effect=rewrite_executable_after_copy,
            ):
                report = immutable_orion_tree.validate_executable_elf(
                    executable,
                    executable_sha256,
                )
            self.assertTrue(report["elf_identity"])

            dependency_archive = base / "dependency.tar"
            with tarfile.open(dependency_archive, "w") as handle:
                handle.add(original_source, arcname="source.txt")
            replacement_dependency = base / "replacement-dependency.tar"
            with tarfile.open(replacement_dependency, "w") as handle:
                handle.add(replacement_source, arcname="source.txt")
            rewritten = False

            def rewrite_dependency_after_copy(*args: object, **kwargs: object):
                nonlocal rewritten
                copied = original_copy(*args, **kwargs)
                if not rewritten:
                    dependency_archive.write_bytes(replacement_dependency.read_bytes())
                    rewritten = True
                return copied

            with mock.patch.object(
                immutable_orion_tree,
                "_sealed_regular_copy",
                side_effect=rewrite_dependency_after_copy,
            ):
                report = immutable_orion_tree.validate_source_archive_dependencies(
                    dependency_archive,
                    {"source.txt": hashlib.sha256(original_source.read_bytes()).hexdigest()},
                    {"source.txt"},
                )
            self.assertTrue(report["passed"])

    def test_metadata_text_read_rejects_inplace_rewrite_during_read(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "metadata.json"
            path.write_text('{"state": "original"}\n', encoding="utf-8")
            original_read = immutable_orion_tree.os.read
            rewritten = False

            def rewrite_after_first_chunk(fd: int, count: int) -> bytes:
                nonlocal rewritten
                payload = original_read(fd, count)
                if payload and not rewritten:
                    path.write_text('{"state": "replacement"}\n', encoding="utf-8")
                    rewritten = True
                return payload

            with mock.patch.object(
                immutable_orion_tree.os,
                "read",
                side_effect=rewrite_after_first_chunk,
            ):
                with self.assertRaisesRegex(ValueError, "changed while reading"):
                    immutable_orion_tree._read_regular_text(
                        path,
                        error_type=ValueError,
                        label="metadata test",
                    )


if __name__ == "__main__":
    unittest.main()
