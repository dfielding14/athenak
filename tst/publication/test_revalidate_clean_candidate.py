#!/usr/bin/env python3
"""Focused tests for standalone read-only clean-candidate revalidation."""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


SCRIPT_DIR = Path(__file__).absolute().parent
CONTROL_PLANE_DIR = SCRIPT_DIR / "frontier_control_plane"
sys.path.insert(0, str(CONTROL_PLANE_DIR))
import control_plane_common
import revalidate_clean_candidate as revalidator


FREEZE_ID = "12345678-1234-4234-8234-123456789abc"
CURRENT_VERSION = "a" * 64
HISTORICAL_VERSION = "b" * 64


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


class RevalidateCleanCandidateTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.orion_root = self.base / "orion"
        self.project_home_root = self.base / "project-home"
        self.source_root = self.base / "source"
        self.orion_root.mkdir()
        self.project_home_root.mkdir()
        self.source_root.mkdir()
        self.candidate = self._candidate_tree()

    def tearDown(self) -> None:
        root = Path(self.temporary.name)
        if root.exists():
            for path in sorted(root.rglob("*"), reverse=True):
                if path.is_dir() and not path.is_symlink():
                    path.chmod(0o700)
            root.chmod(0o700)
        self.temporary.cleanup()

    def _candidate_tree(self, *, with_submodule: bool = False) -> Path:
        root = self.orion_root / "clean_candidates" / FREEZE_ID
        provenance = root / "build_provenance"
        provenance.mkdir(parents=True)
        for filename in control_plane_common.BUILD_PROVENANCE_FILENAMES.values():
            (provenance / filename).write_bytes(f"{filename}\n".encode("utf-8"))
        (root / "source.tar").write_bytes(b"source archive\n")
        (root / "source.commit").write_bytes(b"source commit\n")
        (root / "build_profile.json").write_bytes(b'{"profile": "fixture"}\n')
        (root / "profile_receipt.json").write_bytes(
            _json_bytes(
                {
                    "schema_version": 1,
                    "control_plane_version": HISTORICAL_VERSION,
                }
            )
        )
        (root / "athena").write_bytes(b"fixture executable\n")
        submodules: list[dict[str, object]] = []
        if with_submodule:
            submodule_root = root / "submodules"
            submodule_root.mkdir()
            (submodule_root / "0000.tar").write_bytes(b"submodule archive\n")
            (submodule_root / "0000.commit").write_bytes(b"submodule commit\n")
            submodules.append(
                {
                    "archive_path": str(submodule_root / "0000.tar"),
                    "commit_path": str(submodule_root / "0000.commit"),
                }
            )
        manifest = {
            "schema_version": 4,
            "freeze_id": FREEZE_ID,
            "created_utc": "2026-06-03T12:34:56Z",
            "prepared_artifacts": {},
            "source": {
                "archive_path": str(root / "source.tar"),
                "commit_path": str(root / "source.commit"),
                "git_commit": "c" * 40,
                "git_tree": "d" * 40,
                "source_bundle_sha256": "e" * 64,
                "submodules": submodules,
            },
            "build": {
                "executable_path": str(root / "athena"),
                "profile_id": "fixture-profile",
                "profile_path": str(root / "build_profile.json"),
                "profile_receipt_path": str(root / "profile_receipt.json"),
            },
        }
        manifest_path = root / "clean_candidate_manifest.json"
        manifest_path.write_bytes(_json_bytes(manifest))
        for path in root.rglob("*"):
            path.chmod(0o555 if path.is_dir() else 0o444)
        root.chmod(0o555)
        return manifest_path

    def _unseal_candidate(self) -> None:
        self.candidate.parent.chmod(0o755)

    def _snapshot(self) -> dict[str, tuple[int, bytes | None]]:
        snapshot: dict[str, tuple[int, bytes | None]] = {}
        for path in sorted(self.base.rglob("*")):
            relative = str(path.relative_to(self.base))
            mode = path.stat(follow_symlinks=False).st_mode
            payload = path.read_bytes() if path.is_file() else None
            snapshot[relative] = (mode, payload)
        return snapshot

    def _invoke(
        self,
        *,
        expected_manifest_sha256: str | None = None,
        verify_installed: revalidator.InstalledVerifier | None = None,
        verify_historical: revalidator.InstalledVerifier | None = None,
        validate_bundle: revalidator.BundleValidator | None = None,
    ) -> tuple[dict[str, object], list[tuple[Path, Path]], list[tuple[Path, Path]]]:
        current_calls: list[tuple[Path, Path]] = []
        historical_calls: list[tuple[Path, Path]] = []

        def current(path: Path, *, authorized_pic_root: Path) -> dict[str, object]:
            current_calls.append((path, authorized_pic_root))
            return {"version": CURRENT_VERSION}

        def historical(path: Path, *, authorized_pic_root: Path) -> dict[str, object]:
            historical_calls.append((path, authorized_pic_root))
            return {"version": HISTORICAL_VERSION}

        def validate(candidate: dict[str, object], **kwargs: object) -> list[dict[str, str]]:
            self.assertEqual(candidate["freeze_id"], FREEZE_ID)
            self.assertEqual(kwargs["authorized_pic_root"], self.orion_root)
            self.assertEqual(kwargs["authorized_source_root"], self.source_root)
            self.assertEqual(kwargs["executable_sha256"], _sha256(b"fixture executable\n"))
            return []

        result = revalidator.revalidate_clean_candidate(
            self.candidate,
            expected_manifest_sha256=(
                _sha256(self.candidate.read_bytes())
                if expected_manifest_sha256 is None
                else expected_manifest_sha256
            ),
            control_plane_dir=(
                self.orion_root / "control_plane" / CURRENT_VERSION
            ),
            authorized_pic_root=self.orion_root,
            authorized_project_home_root=self.project_home_root,
            authorized_source_root=self.source_root,
            verify_installed=verify_installed or current,
            verify_historical=verify_historical or historical,
            validate_bundle=validate_bundle or validate,
        )
        return result, current_calls, historical_calls

    def test_revalidation_is_read_only_and_checks_both_controller_pairs(self) -> None:
        (self.orion_root / "policy.json").write_bytes(b"policy sentinel\n")
        (self.orion_root / "ledger.jsonl").write_bytes(b"ledger sentinel\n")
        (self.orion_root / "scheduler.txt").write_bytes(b"scheduler sentinel\n")
        before = self._snapshot()

        result, current_calls, historical_calls = self._invoke()

        self.assertEqual(self._snapshot(), before)
        self.assertEqual(
            current_calls,
            [
                (
                    self.orion_root / "control_plane" / CURRENT_VERSION,
                    self.orion_root,
                ),
                (
                    self.project_home_root / "control_plane" / CURRENT_VERSION,
                    self.project_home_root,
                ),
            ],
        )
        self.assertEqual(
            historical_calls,
            [
                (
                    self.orion_root / "control_plane" / HISTORICAL_VERSION,
                    self.orion_root,
                ),
                (
                    self.project_home_root / "control_plane" / HISTORICAL_VERSION,
                    self.project_home_root,
                ),
            ],
        )
        self.assertEqual(result["status"], "passed")
        self.assertEqual(result["current_control_plane_version"], CURRENT_VERSION)
        self.assertEqual(
            result["build"]["receipt_control_plane_version"],
            HISTORICAL_VERSION,
        )
        self.assertEqual(
            result["clean_candidate_manifest"]["path"],
            str(self.candidate),
        )
        self.assertEqual(
            result["clean_candidate_manifest"]["expected_sha256"],
            _sha256(self.candidate.read_bytes()),
        )
        self.assertEqual(
            result["clean_candidate_manifest"]["sha256"],
            _sha256(self.candidate.read_bytes()),
        )

    def test_reader_captures_optional_fixed_submodule_closure(self) -> None:
        self.temporary.cleanup()
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.orion_root = self.base / "orion"
        self.project_home_root = self.base / "project-home"
        self.source_root = self.base / "source"
        self.orion_root.mkdir()
        self.project_home_root.mkdir()
        self.source_root.mkdir()
        self.candidate = self._candidate_tree(with_submodule=True)

        tree = control_plane_common.read_clean_candidate_tree(
            self.candidate,
            authorized_pic_root=self.orion_root,
        )

        self.assertEqual(tree["submodule_archives"], [b"submodule archive\n"])
        self.assertEqual(tree["submodule_commits"], [b"submodule commit\n"])

    def test_reader_rejects_extra_candidate_member(self) -> None:
        self._unseal_candidate()
        (self.candidate.parent / "unexpected.txt").write_bytes(b"extra\n")
        (self.candidate.parent / "unexpected.txt").chmod(0o444)
        self.candidate.parent.chmod(0o555)

        with self.assertRaisesRegex(ValueError, "fixed layout"):
            control_plane_common.read_clean_candidate_tree(
                self.candidate,
                authorized_pic_root=self.orion_root,
            )

    def test_reader_rejects_writable_member(self) -> None:
        source_archive = self.candidate.parent / "source.tar"
        source_archive.chmod(0o644)

        with self.assertRaisesRegex(ValueError, "not read-only"):
            control_plane_common.read_clean_candidate_tree(
                self.candidate,
                authorized_pic_root=self.orion_root,
            )

    def test_reader_rejects_hard_linked_member(self) -> None:
        self._unseal_candidate()
        source_archive = self.candidate.parent / "source.tar"
        source_archive.unlink()
        outside = self.base / "outside-source.tar"
        outside.write_bytes(b"source archive\n")
        outside.chmod(0o444)
        os.link(outside, source_archive)
        self.candidate.parent.chmod(0o555)

        with self.assertRaisesRegex(ValueError, "hard link"):
            control_plane_common.read_clean_candidate_tree(
                self.candidate,
                authorized_pic_root=self.orion_root,
            )

    def test_reader_rejects_symlinked_member(self) -> None:
        self._unseal_candidate()
        source_archive = self.candidate.parent / "source.tar"
        source_archive.unlink()
        outside = self.base / "outside-source.tar"
        outside.write_bytes(b"source archive\n")
        outside.chmod(0o444)
        source_archive.symlink_to(outside)
        self.candidate.parent.chmod(0o555)

        with self.assertRaises(OSError):
            control_plane_common.read_clean_candidate_tree(
                self.candidate,
                authorized_pic_root=self.orion_root,
            )

    def test_reader_closes_watcher_if_candidate_root_is_missing(self) -> None:
        candidate_root = self.orion_root / "clean_candidates"
        held = self.orion_root / "held-clean-candidates"
        candidate_root.rename(held)
        before = len(os.listdir(f"/proc/{os.getpid()}/fd"))
        try:
            with self.assertRaises(FileNotFoundError):
                control_plane_common.read_clean_candidate_tree(
                    self.candidate,
                    authorized_pic_root=self.orion_root,
                )
        finally:
            held.rename(candidate_root)
        after = len(os.listdir(f"/proc/{os.getpid()}/fd"))
        self.assertEqual(after, before)

    def test_reader_rejects_substitution_of_already_read_source_member(self) -> None:
        read_regular = control_plane_common._read_clean_candidate_regular_file_at
        substituted = False

        def substitute_source(
            directory_descriptor: int, name: str, *, label: str
        ) -> bytes:
            nonlocal substituted
            payload = read_regular(directory_descriptor, name, label=label)
            if label == "Clean-candidate source commit object" and not substituted:
                substituted = True
                self._unseal_candidate()
                source_archive = self.candidate.parent / "source.tar"
                source_archive.unlink()
                source_archive.write_bytes(b"coherent substitute source archive\n")
                source_archive.chmod(0o444)
                self.candidate.parent.chmod(0o555)
            return payload

        with self.assertRaisesRegex(ValueError, "retained closure"):
            control_plane_common.read_clean_candidate_tree(
                self.candidate,
                authorized_pic_root=self.orion_root,
                read_regular_file_at=substitute_source,
            )
        self.assertTrue(substituted)

    def test_reader_rejects_substitution_of_already_read_manifest(self) -> None:
        read_regular = control_plane_common._read_clean_candidate_regular_file_at
        substituted = False

        def substitute_manifest(
            directory_descriptor: int, name: str, *, label: str
        ) -> bytes:
            nonlocal substituted
            payload = read_regular(directory_descriptor, name, label=label)
            if label == "Clean-candidate source archive" and not substituted:
                substituted = True
                manifest_payload = self.candidate.read_bytes()
                self._unseal_candidate()
                self.candidate.unlink()
                self.candidate.write_bytes(manifest_payload)
                self.candidate.chmod(0o444)
                self.candidate.parent.chmod(0o555)
            return payload

        with self.assertRaisesRegex(ValueError, "retained closure"):
            control_plane_common.read_clean_candidate_tree(
                self.candidate,
                authorized_pic_root=self.orion_root,
                read_regular_file_at=substitute_manifest,
            )
        self.assertTrue(substituted)

    def test_reader_rejects_manifest_rename_away_and_restore(self) -> None:
        read_regular = control_plane_common._read_clean_candidate_regular_file_at
        renamed = False

        def rename_manifest(
            directory_descriptor: int, name: str, *, label: str
        ) -> bytes:
            nonlocal renamed
            payload = read_regular(directory_descriptor, name, label=label)
            if label == "Clean-candidate source archive" and not renamed:
                renamed = True
                held = self.candidate.with_name("held-manifest.json")
                self._unseal_candidate()
                self.candidate.rename(held)
                held.rename(self.candidate)
                self.candidate.parent.chmod(0o555)
            return payload

        with self.assertRaisesRegex(ValueError, "retained closure"):
            control_plane_common.read_clean_candidate_tree(
                self.candidate,
                authorized_pic_root=self.orion_root,
                read_regular_file_at=rename_manifest,
            )
        self.assertTrue(renamed)

    def test_reader_rejects_persistent_clean_candidates_parent_substitution(self) -> None:
        read_regular = control_plane_common._read_clean_candidate_regular_file_at
        substituted = False

        def substitute_parent(
            directory_descriptor: int, name: str, *, label: str
        ) -> bytes:
            nonlocal substituted
            payload = read_regular(directory_descriptor, name, label=label)
            if label == "Clean-candidate build profile" and not substituted:
                substituted = True
                candidate_root = self.orion_root / "clean_candidates"
                candidate_root.rename(self.orion_root / "held-clean-candidates")
                candidate_root.mkdir()
            return payload

        with self.assertRaisesRegex(ValueError, "retained closure|ancestry"):
            control_plane_common.read_clean_candidate_tree(
                self.candidate,
                authorized_pic_root=self.orion_root,
                read_regular_file_at=substitute_parent,
            )
        self.assertTrue(substituted)

    def test_reader_rejects_clean_candidates_parent_rename_away_and_restore(self) -> None:
        read_regular = control_plane_common._read_clean_candidate_regular_file_at
        renamed = False

        def rename_parent(
            directory_descriptor: int, name: str, *, label: str
        ) -> bytes:
            nonlocal renamed
            payload = read_regular(directory_descriptor, name, label=label)
            if label == "Clean-candidate build profile" and not renamed:
                renamed = True
                candidate_root = self.orion_root / "clean_candidates"
                held = self.orion_root / "held-clean-candidates"
                candidate_root.rename(held)
                held.rename(candidate_root)
            return payload

        with self.assertRaisesRegex(ValueError, "retained closure"):
            control_plane_common.read_clean_candidate_tree(
                self.candidate,
                authorized_pic_root=self.orion_root,
                read_regular_file_at=rename_parent,
            )
        self.assertTrue(renamed)

    def test_reader_rejects_late_manifest_rename_away_and_restore(self) -> None:
        require_same = control_plane_common.PinnedDirectoryAncestry.require_same
        ancestry_checks = 0
        renamed = False

        def rename_manifest_after_ancestry(
            ancestry: control_plane_common.PinnedDirectoryAncestry,
        ) -> None:
            nonlocal ancestry_checks, renamed
            require_same(ancestry)
            ancestry_checks += 1
            if ancestry_checks == 2:
                renamed = True
                held = self.candidate.with_name("held-manifest.json")
                self._unseal_candidate()
                self.candidate.rename(held)
                held.rename(self.candidate)
                self.candidate.parent.chmod(0o555)

        with mock.patch.object(
            control_plane_common.PinnedDirectoryAncestry,
            "require_same",
            rename_manifest_after_ancestry,
        ):
            with self.assertRaisesRegex(ValueError, "retained closure"):
                control_plane_common.read_clean_candidate_tree(
                    self.candidate,
                    authorized_pic_root=self.orion_root,
                )
        self.assertTrue(renamed)

    def test_reader_rejects_late_clean_candidates_parent_rename_away_and_restore(
        self,
    ) -> None:
        require_same = control_plane_common.PinnedDirectoryAncestry.require_same
        ancestry_checks = 0
        renamed = False

        def rename_parent_after_ancestry(
            ancestry: control_plane_common.PinnedDirectoryAncestry,
        ) -> None:
            nonlocal ancestry_checks, renamed
            require_same(ancestry)
            ancestry_checks += 1
            if ancestry_checks == 2:
                renamed = True
                candidate_root = self.orion_root / "clean_candidates"
                held = self.orion_root / "held-clean-candidates"
                candidate_root.rename(held)
                held.rename(candidate_root)

        with mock.patch.object(
            control_plane_common.PinnedDirectoryAncestry,
            "require_same",
            rename_parent_after_ancestry,
        ):
            with self.assertRaisesRegex(ValueError, "retained closure"):
                control_plane_common.read_clean_candidate_tree(
                    self.candidate,
                    authorized_pic_root=self.orion_root,
                )
        self.assertTrue(renamed)

    def test_reader_rejects_manifest_rename_during_retained_descriptor_close(
        self,
    ) -> None:
        close = control_plane_common.os.close
        renamed = False

        def close_and_rename_manifest(descriptor: int) -> None:
            nonlocal renamed
            caller = sys._getframe(1)
            close(descriptor)
            if (
                not renamed
                and caller.f_code
                is control_plane_common._RetainedCleanCandidateClosure.close.__code__
            ):
                renamed = True
                held = self.candidate.with_name("held-manifest.json")
                self._unseal_candidate()
                self.candidate.rename(held)
                held.rename(self.candidate)
                self.candidate.parent.chmod(0o555)

        with mock.patch.object(
            control_plane_common.os,
            "close",
            close_and_rename_manifest,
        ):
            with self.assertRaisesRegex(ValueError, "retained closure"):
                control_plane_common.read_clean_candidate_tree(
                    self.candidate,
                    authorized_pic_root=self.orion_root,
                )
        self.assertTrue(renamed)

    def test_reader_rejects_parent_rename_during_retained_descriptor_close(self) -> None:
        close = control_plane_common.os.close
        renamed = False

        def close_and_rename_parent(descriptor: int) -> None:
            nonlocal renamed
            caller = sys._getframe(1)
            close(descriptor)
            if (
                not renamed
                and caller.f_code
                is control_plane_common._RetainedCleanCandidateClosure.close.__code__
            ):
                renamed = True
                candidate_root = self.orion_root / "clean_candidates"
                held = self.orion_root / "held-clean-candidates"
                candidate_root.rename(held)
                held.rename(candidate_root)

        with mock.patch.object(
            control_plane_common.os,
            "close",
            close_and_rename_parent,
        ):
            with self.assertRaisesRegex(ValueError, "retained closure"):
                control_plane_common.read_clean_candidate_tree(
                    self.candidate,
                    authorized_pic_root=self.orion_root,
                )
        self.assertTrue(renamed)

    def test_reader_rejects_manifest_rename_during_ancestry_descriptor_close(
        self,
    ) -> None:
        close = control_plane_common.os.close
        renamed = False

        def close_and_rename_manifest(descriptor: int) -> None:
            nonlocal renamed
            caller = sys._getframe(1)
            close(descriptor)
            if (
                not renamed
                and caller.f_code
                is control_plane_common.PinnedDirectoryAncestry.close.__code__
            ):
                renamed = True
                held = self.candidate.with_name("held-manifest.json")
                self._unseal_candidate()
                self.candidate.rename(held)
                held.rename(self.candidate)
                self.candidate.parent.chmod(0o555)

        with mock.patch.object(
            control_plane_common.os,
            "close",
            close_and_rename_manifest,
        ):
            with self.assertRaisesRegex(ValueError, "retained closure"):
                control_plane_common.read_clean_candidate_tree(
                    self.candidate,
                    authorized_pic_root=self.orion_root,
                )
        self.assertTrue(renamed)

    def test_reader_rejects_parent_rename_during_ancestry_descriptor_close(self) -> None:
        close = control_plane_common.os.close
        renamed = False

        def close_and_rename_parent(descriptor: int) -> None:
            nonlocal renamed
            caller = sys._getframe(1)
            close(descriptor)
            if (
                not renamed
                and caller.f_code
                is control_plane_common.PinnedDirectoryAncestry.close.__code__
            ):
                renamed = True
                candidate_root = self.orion_root / "clean_candidates"
                held = self.orion_root / "held-clean-candidates"
                candidate_root.rename(held)
                held.rename(candidate_root)

        with mock.patch.object(
            control_plane_common.os,
            "close",
            close_and_rename_parent,
        ):
            with self.assertRaisesRegex(ValueError, "retained closure"):
                control_plane_common.read_clean_candidate_tree(
                    self.candidate,
                    authorized_pic_root=self.orion_root,
                )
        self.assertTrue(renamed)

    def test_revalidation_propagates_missing_historical_controller(self) -> None:
        def missing(path: Path, *, authorized_pic_root: Path) -> dict[str, object]:
            del path, authorized_pic_root
            raise ValueError("Missing historical installed control-plane directory")

        with self.assertRaisesRegex(ValueError, "Missing historical"):
            self._invoke(verify_historical=missing)

    def test_revalidation_requires_lowercase_expected_manifest_sha256(self) -> None:
        verify_installed = mock.Mock()
        for malformed in ["", "A" * 64, "0" * 63, "not-a-digest"]:
            with self.subTest(malformed=malformed):
                with self.assertRaisesRegex(ValueError, "lowercase SHA-256"):
                    self._invoke(
                        expected_manifest_sha256=malformed,
                        verify_installed=verify_installed,
                    )
        verify_installed.assert_not_called()

    def test_revalidation_rejects_expected_manifest_sha256_mismatch_before_semantic_validation(
        self,
    ) -> None:
        validate_bundle = mock.Mock()

        with self.assertRaisesRegex(ValueError, "differs from expected binding"):
            self._invoke(
                expected_manifest_sha256="0" * 64,
                validate_bundle=validate_bundle,
            )

        validate_bundle.assert_not_called()

    def test_cli_has_no_root_override_and_prints_only_canonical_json(self) -> None:
        parser = revalidator._parser()
        expected_sha256 = _sha256(self.candidate.read_bytes())
        parsed = parser.parse_args(
            [
                "--manifest",
                str(self.candidate),
                "--expected-manifest-sha256",
                expected_sha256,
            ]
        )
        self.assertEqual(parsed.candidate_manifest_path, self.candidate)
        self.assertEqual(parsed.expected_manifest_sha256, expected_sha256)
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parser.parse_args(
                    [
                        "--manifest",
                        str(self.candidate),
                        "--expected-manifest-sha256",
                        expected_sha256,
                        "--orion-root",
                        str(self.orion_root),
                    ]
                )

        expected = {"schema_version": 1, "status": "passed"}
        output = io.StringIO()
        with mock.patch.object(
            revalidator,
            "revalidate_clean_candidate",
            return_value=expected,
        ) as revalidate, mock.patch.object(
            sys,
            "_pic_control_plane_bootstrapped",
            True,
            create=True,
        ), mock.patch.object(
            sys,
            "argv",
            [
                revalidator.ENTRYPOINT_NAME,
                "--manifest",
                str(self.candidate),
                "--expected-manifest-sha256",
                expected_sha256,
            ],
        ), redirect_stdout(output):
            revalidator.main()
        revalidate.assert_called_once_with(
            self.candidate,
            expected_manifest_sha256=expected_sha256,
        )
        self.assertEqual(
            output.getvalue(),
            revalidator._canonical_json_bytes(expected).decode("utf-8"),
        )

    def test_cli_rejects_direct_execution_without_authenticated_runner(self) -> None:
        with mock.patch.object(
            sys,
            "_pic_control_plane_bootstrapped",
            False,
            create=True,
        ):
            with self.assertRaisesRegex(SystemExit, "through run_control_plane.py"):
                revalidator.main()


if __name__ == "__main__":
    unittest.main()
