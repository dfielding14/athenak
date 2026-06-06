#!/usr/bin/env python3
"""Focused tests for authenticated fixed-root storage-preflight evidence."""

from __future__ import annotations

from contextlib import contextmanager, redirect_stderr, redirect_stdout
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import stat
import sys
import tempfile
import unittest
from unittest import mock


SCRIPT_DIR = Path(__file__).absolute().parent
CONTROL_PLANE_DIR = SCRIPT_DIR / "frontier_control_plane"
sys.path.insert(0, str(CONTROL_PLANE_DIR))
import capture_storage_preflight_evidence as storage_preflight
import run_control_plane


FIXED_TIME = datetime(2026, 6, 3, 12, 34, 56, tzinfo=timezone.utc)
PROBE_ID = "12345678-1234-4234-8234-123456789abc"
SOURCE_AUTHENTICATION = {
    "common_sha256": "e" * 64,
    "entrypoint_sha256": "a" * 64,
    "git_commit": "b" * 40,
    "runner_sha256": "c" * 64,
    "schema_sha256": "d" * 64,
    "tracked_clean_head_blobs": True,
}


class SourceAuthenticator:
    """Return deterministic reviewed-source records while counting calls."""

    def __init__(self, records: list[dict[str, object]] | None = None) -> None:
        self.records = records or [SOURCE_AUTHENTICATION]
        self.calls = 0

    def __call__(self) -> dict[str, object]:
        record = self.records[min(self.calls, len(self.records) - 1)]
        self.calls += 1
        return dict(record)


class CaptureStoragePreflightEvidenceTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.orion_root = self.base / "orion"
        self.project_home_root = self.base / "project-home"
        self.orion_root.mkdir()
        self.project_home_root.mkdir()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _capture(
        self,
        *,
        authenticate_source: SourceAuthenticator | None = None,
    ) -> dict[str, object]:
        return storage_preflight.capture_storage_preflight_evidence(
            expected_git_commit=str(SOURCE_AUTHENTICATION["git_commit"]),
            orion_root=self.orion_root,
            project_home_root=self.project_home_root,
            probe_id=PROBE_ID,
            now=lambda: FIXED_TIME,
            token_bytes=lambda size: bytes(range(size)),
            authenticate_source=authenticate_source or SourceAuthenticator(),
        )

    def _evidence_paths(self) -> tuple[Path, Path]:
        relative = Path("policy") / "storage_preflight_evidence" / f"{PROBE_ID}.json"
        return self.orion_root / relative, self.project_home_root / relative

    def _recover(
        self,
        *,
        expected_sha256: str,
        expected_existing_role: str = storage_preflight.ORION_ROLE,
        authenticate_source: SourceAuthenticator | None = None,
    ) -> dict[str, object]:
        return storage_preflight.recover_storage_preflight_evidence_pair(
            probe_id=PROBE_ID,
            expected_sha256=expected_sha256,
            expected_existing_role=expected_existing_role,
            expected_git_commit=str(SOURCE_AUTHENTICATION["git_commit"]),
            orion_root=self.orion_root,
            project_home_root=self.project_home_root,
            authenticate_source=authenticate_source or SourceAuthenticator(),
        )

    def _assert_no_probe_files(self) -> None:
        for root in [self.orion_root, self.project_home_root]:
            self.assertEqual(
                [path.name for path in root.iterdir() if path.name.startswith(".pic-")],
                [],
            )

    def test_capture_publishes_byte_identical_read_only_evidence_and_policy_fragment(
        self,
    ) -> None:
        authenticate_source = SourceAuthenticator()
        fragment = self._capture(authenticate_source=authenticate_source)
        orion_path, project_home_path = self._evidence_paths()

        self.assertEqual(authenticate_source.calls, 4)
        self.assertEqual(orion_path.read_bytes(), project_home_path.read_bytes())
        self.assertFalse(orion_path.stat().st_mode & 0o222)
        self.assertFalse(project_home_path.stat().st_mode & 0o222)
        self._assert_no_probe_files()
        payload = orion_path.read_bytes()
        artifact = json.loads(payload)
        self.assertEqual(artifact["schema_version"], 2)
        self.assertEqual(artifact["record_type"], storage_preflight.RECORD_TYPE)
        self.assertEqual(artifact["probe_id"], PROBE_ID)
        self.assertEqual(artifact["method"], storage_preflight.METHOD)
        self.assertEqual(artifact["source_authentication"], SOURCE_AUTHENTICATION)
        self.assertEqual(
            [probe["role"] for probe in artifact["probes"]],
            ["orion_simulation_root", "project_home_mirror_root"],
        )
        for probe in artifact["probes"]:
            self.assertEqual(probe["operations"], storage_preflight.OPERATIONS)
            self.assertEqual(probe["status"], "passed")
            self.assertEqual(probe["payload_bytes"], storage_preflight.PAYLOAD_BYTES)
        self.assertEqual(
            fragment["storage_preflight_evidence"],
            {
                "orion_path": str(orion_path),
                "probe_id": PROBE_ID,
                "project_home_path": str(project_home_path),
                "sha256": hashlib.sha256(payload).hexdigest(),
            },
        )

    def test_capture_loops_until_short_writes_are_complete(self) -> None:
        real_write = os.write

        def short_write(descriptor: int, payload: bytes) -> int:
            return real_write(descriptor, payload[:3])

        with mock.patch.object(storage_preflight.os, "write", side_effect=short_write):
            self._capture()

        orion_path, project_home_path = self._evidence_paths()
        self.assertEqual(orion_path.read_bytes(), project_home_path.read_bytes())
        self._assert_no_probe_files()

    def test_source_authentication_is_required_before_root_mutation(self) -> None:
        malformed = {**SOURCE_AUTHENTICATION, "runner_sha256": "not-a-digest"}
        with self.assertRaisesRegex(ValueError, "digest is malformed"):
            self._capture(authenticate_source=SourceAuthenticator([malformed]))

        self.assertEqual(list(self.orion_root.iterdir()), [])
        self.assertEqual(list(self.project_home_root.iterdir()), [])

    def test_expected_git_commit_is_required_before_root_mutation(self) -> None:
        with self.assertRaisesRegex(ValueError, "differs from expected Git commit"):
            storage_preflight.capture_storage_preflight_evidence(
                expected_git_commit="e" * 40,
                orion_root=self.orion_root,
                project_home_root=self.project_home_root,
                probe_id=PROBE_ID,
                now=lambda: FIXED_TIME,
                token_bytes=lambda size: bytes(range(size)),
                authenticate_source=SourceAuthenticator(),
            )

        self.assertEqual(list(self.orion_root.iterdir()), [])
        self.assertEqual(list(self.project_home_root.iterdir()), [])

    def test_source_change_after_root_probes_fails_before_publication(self) -> None:
        changed = {**SOURCE_AUTHENTICATION, "git_commit": "e" * 40}
        authenticate_source = SourceAuthenticator(
            [SOURCE_AUTHENTICATION, SOURCE_AUTHENTICATION, changed]
        )
        with self.assertRaisesRegex(ValueError, "changed after probing roots"):
            self._capture(authenticate_source=authenticate_source)

        self._assert_no_probe_files()
        self.assertFalse((self.orion_root / "policy").exists())
        self.assertFalse((self.project_home_root / "policy").exists())

    def test_readback_mismatch_fails_closed_and_removes_probe_file(self) -> None:
        with mock.patch.object(storage_preflight, "_read_all", return_value=b"x" * 32):
            with self.assertRaisesRegex(ValueError, "readback differs"):
                self._capture()

        self._assert_no_probe_files()
        self.assertEqual(list(self.orion_root.iterdir()), [])
        self.assertEqual(list(self.project_home_root.iterdir()), [])

    def test_failed_probe_write_removes_exact_created_file(self) -> None:
        with mock.patch.object(
            storage_preflight,
            "_write_all",
            side_effect=OSError("injected write failure"),
        ), self.assertRaisesRegex(OSError, "injected write failure"):
            self._capture()

        self._assert_no_probe_files()
        self.assertEqual(list(self.orion_root.iterdir()), [])
        self.assertEqual(list(self.project_home_root.iterdir()), [])

    def test_probe_cleanup_preserves_substituted_namespace_entry(self) -> None:
        real_require = storage_preflight._require_same_regular_entry
        replacement = b"replacement must be preserved\n"
        substituted = False

        def substitute_before_cleanup(
            directory_descriptor: int,
            name: str,
            expected: os.stat_result,
            *,
            label: str,
        ) -> os.stat_result:
            nonlocal substituted
            if label.endswith("storage probe") and not substituted:
                substituted = True
                os.unlink(name, dir_fd=directory_descriptor)
                descriptor = os.open(
                    name,
                    storage_preflight.FILE_CREATE_FLAGS,
                    0o600,
                    dir_fd=directory_descriptor,
                )
                try:
                    os.write(descriptor, replacement)
                finally:
                    os.close(descriptor)
            return real_require(
                directory_descriptor,
                name,
                expected,
                label=label,
            )

        with mock.patch.object(
            storage_preflight,
            "_require_same_regular_entry",
            side_effect=substitute_before_cleanup,
        ), self.assertRaisesRegex(ValueError, "namespace entry changed"):
            self._capture()

        residue = [
            path
            for path in self.orion_root.iterdir()
            if path.name.startswith(".pic-storage-preflight-")
        ]
        self.assertEqual(len(residue), 1)
        self.assertEqual(residue[0].read_bytes(), replacement)

    def test_symlink_root_component_is_rejected_without_redirected_write(self) -> None:
        outside = self.base / "outside"
        outside.mkdir()
        self.orion_root.rmdir()
        self.orion_root.symlink_to(outside, target_is_directory=True)

        with self.assertRaises(OSError):
            self._capture()

        self.assertEqual(list(outside.iterdir()), [])
        self.assertEqual(list(self.project_home_root.iterdir()), [])

    def test_root_identity_replacement_before_publication_is_rejected(self) -> None:
        moved = self.base / "orion-moved"
        real_publish = storage_preflight._publish_staged_evidence
        replaced = False

        def replace_before_publish(
            root: Path,
            *,
            expected_root_identity: tuple[int, int],
            filename: str,
            payload: bytes,
        ) -> Path:
            nonlocal replaced
            if root == self.orion_root and not replaced:
                replaced = True
                self.orion_root.rename(moved)
                self.orion_root.mkdir()
            return real_publish(
                root,
                expected_root_identity=expected_root_identity,
                filename=filename,
                payload=payload,
            )

        with mock.patch.object(
            storage_preflight,
            "_publish_staged_evidence",
            side_effect=replace_before_publish,
        ):
            with self.assertRaisesRegex(ValueError, "changed after its storage probe"):
                self._capture()

        self.assertFalse((self.orion_root / "policy").exists())
        self.assertFalse((moved / "policy").exists())
        self.assertFalse((self.project_home_root / "policy").exists())

    def test_root_identity_replacement_before_probe_is_rejected_without_write(
        self,
    ) -> None:
        moved = self.base / "orion-moved-before-probe"
        real_probe = storage_preflight._probe_root
        replaced = False

        def replace_before_probe(*args: object, **kwargs: object) -> object:
            nonlocal replaced
            if kwargs["role"] == storage_preflight.ORION_ROLE and not replaced:
                replaced = True
                self.orion_root.rename(moved)
                self.orion_root.mkdir()
            return real_probe(*args, **kwargs)

        with mock.patch.object(
            storage_preflight,
            "_probe_root",
            side_effect=replace_before_probe,
        ), self.assertRaisesRegex(ValueError, "changed after its storage probe"):
            self._capture()

        self.assertEqual(list(self.orion_root.iterdir()), [])
        self.assertEqual(list(moved.iterdir()), [])
        self.assertEqual(list(self.project_home_root.iterdir()), [])

    def test_second_publication_failure_preserves_first_artifact_for_recovery(self) -> None:
        real_publish = storage_preflight._publish_staged_evidence

        def fail_project_home(
            root: Path,
            *,
            expected_root_identity: tuple[int, int],
            filename: str,
            payload: bytes,
        ) -> Path:
            if root == self.project_home_root:
                raise OSError("injected Project Home publication failure")
            return real_publish(
                root,
                expected_root_identity=expected_root_identity,
                filename=filename,
                payload=payload,
            )

        with mock.patch.object(
            storage_preflight,
            "_publish_staged_evidence",
            side_effect=fail_project_home,
        ):
            with self.assertRaisesRegex(OSError, "Project Home publication failure"):
                self._capture()

        orion_path, project_home_path = self._evidence_paths()
        self.assertTrue(orion_path.exists())
        self.assertEqual(stat.S_IMODE(orion_path.stat().st_mode), 0o400)
        self.assertFalse(project_home_path.exists())

    def test_capture_holds_serialization_lock_across_probes_and_publication(self) -> None:
        held = False
        real_lock = storage_preflight._serialization_anchor_lock
        real_probe = storage_preflight._probe_root
        real_publish = storage_preflight._publish_staged_evidence

        @contextmanager
        def record_lock(root: Path):
            nonlocal held
            with real_lock(root):
                held = True
                try:
                    yield
                finally:
                    held = False

        def require_lock_for_probe(*args: object, **kwargs: object) -> object:
            self.assertTrue(held)
            return real_probe(*args, **kwargs)

        def require_lock_for_publish(*args: object, **kwargs: object) -> object:
            self.assertTrue(held)
            return real_publish(*args, **kwargs)

        with mock.patch.object(
            storage_preflight,
            "_serialization_anchor_lock",
            side_effect=record_lock,
        ), mock.patch.object(
            storage_preflight,
            "_probe_root",
            side_effect=require_lock_for_probe,
        ), mock.patch.object(
            storage_preflight,
            "_publish_staged_evidence",
            side_effect=require_lock_for_publish,
        ):
            self._capture()
        self.assertFalse(held)

    def test_capture_rejects_another_one_sided_pair_before_probing(self) -> None:
        other_probe_id = "aaaaaaaa-aaaa-4aaa-8aaa-aaaaaaaaaaaa"
        storage_preflight.capture_storage_preflight_evidence(
            expected_git_commit=str(SOURCE_AUTHENTICATION["git_commit"]),
            orion_root=self.orion_root,
            project_home_root=self.project_home_root,
            probe_id=other_probe_id,
            now=lambda: FIXED_TIME,
            token_bytes=lambda size: bytes(range(size)),
            authenticate_source=SourceAuthenticator(),
        )
        other_project_home = (
            self.project_home_root
            / "policy"
            / "storage_preflight_evidence"
            / f"{other_probe_id}.json"
        )
        other_project_home.unlink()

        with mock.patch.object(storage_preflight, "_probe_root") as probe:
            with self.assertRaisesRegex(ValueError, "another one-sided pair"):
                self._capture()

        probe.assert_not_called()

    def test_existing_evidence_collision_is_not_replaced(self) -> None:
        orion_path, _ = self._evidence_paths()
        orion_path.parent.mkdir(parents=True)
        orion_path.write_bytes(b"existing\n")
        orion_path.chmod(0o400)

        with self.assertRaises(ValueError):
            self._capture()

        self.assertEqual(orion_path.read_bytes(), b"existing\n")

    def test_absent_exact_pair_audit_is_filesystem_read_only(self) -> None:
        before = sorted(
            str(path.relative_to(self.base)) for path in self.base.rglob("*")
        )

        result = storage_preflight.audit_storage_preflight_evidence_pair(
            probe_id=PROBE_ID,
            expected_sha256="a" * 64,
            expected_git_commit=str(SOURCE_AUTHENTICATION["git_commit"]),
            orion_root=self.orion_root,
            project_home_root=self.project_home_root,
            authenticate_source=SourceAuthenticator(),
        )

        self.assertEqual(result["state"], "absent_both")
        self.assertEqual(
            sorted(str(path.relative_to(self.base)) for path in self.base.rglob("*")),
            before,
        )

    def test_exact_recovery_completes_one_sided_evidence_pair(self) -> None:
        fragment = self._capture()
        orion_path, project_home_path = self._evidence_paths()
        payload = orion_path.read_bytes()
        project_home_path.unlink()
        authenticate_source = SourceAuthenticator()

        recovered = self._recover(
            expected_sha256=hashlib.sha256(payload).hexdigest(),
            authenticate_source=authenticate_source,
        )

        self.assertEqual(authenticate_source.calls, 3)
        self.assertEqual(recovered, fragment)
        self.assertEqual(orion_path.read_bytes(), payload)
        self.assertEqual(project_home_path.read_bytes(), payload)
        self.assertFalse(project_home_path.stat().st_mode & 0o222)

    def test_exact_recovery_completes_project_home_sided_pair_and_is_idempotent(
        self,
    ) -> None:
        fragment = self._capture()
        orion_path, project_home_path = self._evidence_paths()
        payload = project_home_path.read_bytes()
        orion_path.unlink()
        recovered = self._recover(
            expected_sha256=hashlib.sha256(payload).hexdigest(),
            expected_existing_role=storage_preflight.PROJECT_HOME_ROLE,
        )
        metadata = tuple(path.stat().st_mtime_ns for path in self._evidence_paths())

        with mock.patch.object(
            storage_preflight, "_publish_staged_evidence"
        ) as publish:
            again = self._recover(
                expected_sha256=hashlib.sha256(payload).hexdigest(),
                expected_existing_role=storage_preflight.PROJECT_HOME_ROLE,
            )

        publish.assert_not_called()
        self.assertEqual(recovered, fragment)
        self.assertEqual(again, fragment)
        self.assertEqual(
            tuple(path.stat().st_mtime_ns for path in self._evidence_paths()),
            metadata,
        )

    def test_exact_pair_audit_reports_absent_one_sided_and_complete_states(self) -> None:
        expected_sha256 = "a" * 64
        audit = lambda: storage_preflight.audit_storage_preflight_evidence_pair(
            probe_id=PROBE_ID,
            expected_sha256=expected_sha256,
            expected_git_commit=str(SOURCE_AUTHENTICATION["git_commit"]),
            orion_root=self.orion_root,
            project_home_root=self.project_home_root,
            authenticate_source=SourceAuthenticator(),
        )
        self.assertEqual(audit()["state"], "absent_both")

        fragment = self._capture()
        expected_sha256 = fragment["storage_preflight_evidence"]["sha256"]
        _, project_home_path = self._evidence_paths()
        project_home_path.unlink()
        self.assertEqual(audit()["state"], "valid_orion_only")
        self._recover(expected_sha256=expected_sha256)
        self.assertEqual(audit()["state"], "valid_identical_pair")

    def test_exact_recovery_rejects_divergent_or_missing_evidence_without_mutation(
        self,
    ) -> None:
        with self.assertRaisesRegex(ValueError, "found no published evidence"):
            self._recover(
                expected_sha256="a" * 64,
            )
        fragment = self._capture()
        orion_path, project_home_path = self._evidence_paths()
        original = orion_path.read_bytes()
        project_home_path.chmod(0o600)
        project_home_path.write_bytes(b"{}\n")
        project_home_path.chmod(0o400)
        with self.assertRaisesRegex(ValueError, "recovery evidence differs"):
            self._recover(
                expected_sha256=fragment["storage_preflight_evidence"]["sha256"],
            )
        self.assertEqual(orion_path.read_bytes(), original)
        self.assertEqual(project_home_path.read_bytes(), b"{}\n")

    def test_exact_recovery_rejects_wrong_role_writable_and_staging_residue(self) -> None:
        fragment = self._capture()
        expected_sha256 = fragment["storage_preflight_evidence"]["sha256"]
        orion_path, project_home_path = self._evidence_paths()
        project_home_path.unlink()
        with self.assertRaisesRegex(ValueError, "existing role differs"):
            self._recover(
                expected_sha256=expected_sha256,
                expected_existing_role=storage_preflight.PROJECT_HOME_ROLE,
            )
        orion_path.chmod(0o600)
        with self.assertRaisesRegex(ValueError, "mode is not exactly 0400"):
            self._recover(expected_sha256=expected_sha256)
        orion_path.chmod(0o400)
        staging = project_home_path.parent / (
            f".{project_home_path.name}.manual{storage_preflight.RECOVERY_STAGING_SUFFIX}"
        )
        project_home_path.parent.mkdir(parents=True, exist_ok=True)
        staging.write_bytes(b"preserved")
        with self.assertRaisesRegex(ValueError, "unresolved recovery residue"):
            self._recover(expected_sha256=expected_sha256)
        self.assertEqual(staging.read_bytes(), b"preserved")
        self.assertFalse(project_home_path.exists())

    def test_exact_recovery_does_not_probe_and_rejects_other_one_sided_pair(self) -> None:
        fragment = self._capture()
        expected_sha256 = fragment["storage_preflight_evidence"]["sha256"]
        _, project_home_path = self._evidence_paths()
        project_home_path.unlink()
        other = (
            self.orion_root
            / "policy"
            / "storage_preflight_evidence"
            / "aaaaaaaa-aaaa-4aaa-8aaa-aaaaaaaaaaaa.json"
        )
        other.write_bytes(b"{}\n")
        other.chmod(0o400)
        with mock.patch.object(storage_preflight, "_probe_root") as probe:
            with self.assertRaisesRegex(ValueError, "another one-sided pair"):
                self._recover(expected_sha256=expected_sha256)
        probe.assert_not_called()

    def test_exact_recovery_accepts_historical_recorded_root_identities(self) -> None:
        self._capture()
        orion_path, project_home_path = self._evidence_paths()
        artifact = json.loads(orion_path.read_bytes())
        for index, probe in enumerate(artifact["probes"], start=1):
            probe["st_dev"] = 1000 + index
            probe["st_ino"] = 2000 + index
        payload = storage_preflight._canonical_json_bytes(artifact)
        for path in [orion_path, project_home_path]:
            path.chmod(0o600)
            path.write_bytes(payload)
            path.chmod(0o400)
        project_home_path.unlink()

        self._recover(expected_sha256=hashlib.sha256(payload).hexdigest())

        self.assertEqual(orion_path.read_bytes(), payload)
        self.assertEqual(project_home_path.read_bytes(), payload)

    def test_exact_recovery_source_drift_before_publication_preserves_one_sided_pair(
        self,
    ) -> None:
        fragment = self._capture()
        expected_sha256 = fragment["storage_preflight_evidence"]["sha256"]
        orion_path, project_home_path = self._evidence_paths()
        payload = orion_path.read_bytes()
        project_home_path.unlink()
        changed = {**SOURCE_AUTHENTICATION, "git_commit": "e" * 40}

        with self.assertRaisesRegex(ValueError, "changed before evidence recovery"):
            self._recover(
                expected_sha256=expected_sha256,
                authenticate_source=SourceAuthenticator(
                    [SOURCE_AUTHENTICATION, changed]
                ),
            )

        self.assertEqual(orion_path.read_bytes(), payload)
        self.assertFalse(project_home_path.exists())

    def test_exact_recovery_source_drift_after_publication_preserves_complete_pair(
        self,
    ) -> None:
        fragment = self._capture()
        expected_sha256 = fragment["storage_preflight_evidence"]["sha256"]
        orion_path, project_home_path = self._evidence_paths()
        payload = orion_path.read_bytes()
        project_home_path.unlink()
        changed = {**SOURCE_AUTHENTICATION, "git_commit": "e" * 40}

        with self.assertRaisesRegex(ValueError, "changed during evidence recovery"):
            self._recover(
                expected_sha256=expected_sha256,
                authenticate_source=SourceAuthenticator(
                    [SOURCE_AUTHENTICATION, SOURCE_AUTHENTICATION, changed]
                ),
            )

        self.assertEqual(orion_path.read_bytes(), payload)
        self.assertEqual(project_home_path.read_bytes(), payload)

    def test_exact_recovery_destination_race_preserves_all_evidence(self) -> None:
        fragment = self._capture()
        expected_sha256 = fragment["storage_preflight_evidence"]["sha256"]
        orion_path, project_home_path = self._evidence_paths()
        project_home_path.unlink()
        real_publish = storage_preflight._publish_staged_evidence

        def race_destination(
            root: Path,
            *,
            expected_root_identity: tuple[int, int],
            filename: str,
            payload: bytes,
        ) -> Path:
            project_home_path.write_bytes(b"racing destination\n")
            project_home_path.chmod(0o400)
            return real_publish(
                root,
                expected_root_identity=expected_root_identity,
                filename=filename,
                payload=payload,
            )

        with mock.patch.object(
            storage_preflight,
            "_publish_staged_evidence",
            side_effect=race_destination,
        ), self.assertRaises(FileExistsError):
            self._recover(expected_sha256=expected_sha256)

        self.assertTrue(orion_path.exists())
        self.assertEqual(project_home_path.read_bytes(), b"racing destination\n")
        staging = list(
            project_home_path.parent.glob(
                f".{project_home_path.name}.*{storage_preflight.RECOVERY_STAGING_SUFFIX}"
            )
        )
        self.assertEqual(len(staging), 1)

    def test_staged_publication_detects_substituted_source_before_cleanup(
        self,
    ) -> None:
        expected_identity = (
            self.orion_root.stat().st_dev,
            self.orion_root.stat().st_ino,
        )
        real_link = storage_preflight.os.link
        replacement = b"replacement must be preserved\n"

        def substitute_source(
            source: str,
            destination: str,
            *,
            src_dir_fd: int,
            dst_dir_fd: int,
            follow_symlinks: bool,
        ) -> None:
            os.unlink(source, dir_fd=src_dir_fd)
            descriptor = os.open(
                source,
                storage_preflight.FILE_CREATE_FLAGS,
                0o400,
                dir_fd=src_dir_fd,
            )
            try:
                os.write(descriptor, replacement)
            finally:
                os.close(descriptor)
            real_link(
                source,
                destination,
                src_dir_fd=src_dir_fd,
                dst_dir_fd=dst_dir_fd,
                follow_symlinks=follow_symlinks,
            )

        with mock.patch.object(
            storage_preflight.os,
            "link",
            side_effect=substitute_source,
        ), self.assertRaisesRegex(ValueError, "namespace entry changed"):
            storage_preflight._publish_staged_evidence(
                self.orion_root,
                expected_root_identity=expected_identity,
                filename=f"{PROBE_ID}.json",
                payload=b"reviewed payload\n",
            )

        parent = self.orion_root.joinpath(*storage_preflight.EVIDENCE_PARENT_PARTS)
        destination = parent / f"{PROBE_ID}.json"
        staging = list(parent.glob(f".{PROBE_ID}.json.*{storage_preflight.RECOVERY_STAGING_SUFFIX}"))
        self.assertEqual(destination.read_bytes(), replacement)
        self.assertEqual(len(staging), 1)
        self.assertEqual(staging[0].read_bytes(), replacement)

    def test_fifo_evidence_entry_fails_without_blocking(self) -> None:
        evidence = self.orion_root.joinpath(*storage_preflight.EVIDENCE_PARENT_PARTS)
        evidence.mkdir(parents=True)
        os.mkfifo(evidence / f"{PROBE_ID}.json", 0o400)

        with self.assertRaisesRegex(ValueError, "not a regular file"):
            storage_preflight.audit_storage_preflight_evidence_pair(
                probe_id=PROBE_ID,
                expected_sha256="a" * 64,
                expected_git_commit=str(SOURCE_AUTHENTICATION["git_commit"]),
                orion_root=self.orion_root,
                project_home_root=self.project_home_root,
                authenticate_source=SourceAuthenticator(),
            )

    def test_exact_recovery_detects_evidence_parent_move_during_publication(self) -> None:
        fragment = self._capture()
        expected_sha256 = fragment["storage_preflight_evidence"]["sha256"]
        _, project_home_path = self._evidence_paths()
        project_home_path.unlink()
        parent = project_home_path.parent
        moved = self.project_home_root / "moved-evidence-parent"
        real_link = storage_preflight.os.link
        moved_once = False

        def move_parent_before_link(*args: object, **kwargs: object) -> None:
            nonlocal moved_once
            if not moved_once:
                moved_once = True
                parent.rename(moved)
                parent.mkdir()
            real_link(*args, **kwargs)

        with mock.patch.object(
            storage_preflight.os,
            "link",
            side_effect=move_parent_before_link,
        ), self.assertRaisesRegex(ValueError, "Evidence parent changed"):
            self._recover(expected_sha256=expected_sha256)

        self.assertFalse(project_home_path.exists())
        self.assertTrue((moved / project_home_path.name).exists())

    def test_exact_recovery_rejects_root_probe_residue_and_symlink_evidence(self) -> None:
        fragment = self._capture()
        expected_sha256 = fragment["storage_preflight_evidence"]["sha256"]
        orion_path, project_home_path = self._evidence_paths()
        project_home_path.unlink()
        residue = self.project_home_root / ".pic-storage-preflight-interrupted"
        residue.write_bytes(b"preserved")
        with self.assertRaisesRegex(ValueError, "probe residue"):
            self._recover(expected_sha256=expected_sha256)
        residue.unlink()
        orion_path.chmod(0o600)
        payload = orion_path.read_bytes()
        orion_path.unlink()
        outside = self.base / "outside-evidence"
        outside.write_bytes(payload)
        outside.chmod(0o400)
        orion_path.symlink_to(outside)
        with self.assertRaises(OSError):
            self._recover(expected_sha256=expected_sha256)
        self.assertEqual(outside.read_bytes(), payload)
        self.assertFalse(project_home_path.exists())

    def test_fsync_is_called_for_regular_files_and_directories(self) -> None:
        real_fsync = os.fsync
        synced_modes: list[int] = []

        def record_fsync(descriptor: int) -> None:
            synced_modes.append(os.fstat(descriptor).st_mode)
            real_fsync(descriptor)

        with mock.patch.object(storage_preflight.os, "fsync", side_effect=record_fsync):
            self._capture()

        self.assertTrue(any(stat.S_ISREG(mode) for mode in synced_modes))
        self.assertTrue(any(stat.S_ISDIR(mode) for mode in synced_modes))

    def test_cli_has_no_root_override_and_prints_canonical_fragment(self) -> None:
        parser = storage_preflight._parser()
        expected_commit = str(SOURCE_AUTHENTICATION["git_commit"])
        self.assertEqual(
            vars(parser.parse_args(["--expected-git-commit", expected_commit])),
            {
                "audit_exact_pair": False,
                "expected_git_commit": expected_commit,
                "expected_evidence_sha256": None,
                "expected_existing_role": None,
                "probe_id": None,
                "recover_exact_pair": False,
            },
        )
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parser.parse_args(
                    [
                        "--expected-git-commit",
                        expected_commit,
                        "--orion-root",
                        str(self.orion_root),
                    ]
                )

        expected = {"policy": "fragment"}
        output = io.StringIO()
        with mock.patch.object(
            storage_preflight,
            "capture_storage_preflight_evidence",
            return_value=expected,
        ) as capture, mock.patch.object(
            storage_preflight,
            "recover_storage_preflight_evidence_pair",
        ), mock.patch.object(sys, "_pic_control_plane_bootstrapped", True, create=True), (
            mock.patch.object(
                sys,
                "argv",
                [
                    storage_preflight.ENTRYPOINT_NAME,
                    "--expected-git-commit",
                    expected_commit,
                ],
            )
        ), redirect_stdout(output):
            storage_preflight.main()
        capture.assert_called_once_with(expected_git_commit=expected_commit)
        self.assertEqual(
            output.getvalue(),
            storage_preflight._canonical_json_bytes(expected).decode(),
        )

    def test_cli_dispatches_exact_recovery_bindings(self) -> None:
        expected_commit = str(SOURCE_AUTHENTICATION["git_commit"])
        expected = {"policy": "fragment"}
        output = io.StringIO()
        with mock.patch.object(
            storage_preflight,
            "recover_storage_preflight_evidence_pair",
            return_value=expected,
        ) as recover, mock.patch.object(
            storage_preflight,
            "capture_storage_preflight_evidence",
        ), mock.patch.object(
            sys, "_pic_control_plane_bootstrapped", True, create=True
        ), mock.patch.object(
            sys,
            "argv",
            [
                storage_preflight.ENTRYPOINT_NAME,
                "--expected-git-commit",
                expected_commit,
                "--recover-exact-pair",
                "--probe-id",
                PROBE_ID,
                "--expected-evidence-sha256",
                "a" * 64,
                "--expected-existing-role",
                storage_preflight.ORION_ROLE,
            ],
        ), redirect_stdout(output):
            storage_preflight.main()
        recover.assert_called_once_with(
            probe_id=PROBE_ID,
            expected_sha256="a" * 64,
            expected_existing_role=storage_preflight.ORION_ROLE,
            expected_git_commit=expected_commit,
        )
        self.assertEqual(
            output.getvalue(),
            storage_preflight._canonical_json_bytes(expected).decode(),
        )

    def test_cli_dispatches_read_only_exact_pair_audit(self) -> None:
        expected_commit = str(SOURCE_AUTHENTICATION["git_commit"])
        expected = {
            "expected_sha256": "a" * 64,
            "probe_id": PROBE_ID,
            "state": "valid_identical_pair",
        }
        output = io.StringIO()
        with mock.patch.object(
            storage_preflight,
            "audit_storage_preflight_evidence_pair",
            return_value=expected,
        ) as audit, mock.patch.object(
            storage_preflight,
            "capture_storage_preflight_evidence",
        ), mock.patch.object(
            storage_preflight,
            "recover_storage_preflight_evidence_pair",
        ), mock.patch.object(
            sys, "_pic_control_plane_bootstrapped", True, create=True
        ), mock.patch.object(
            sys,
            "argv",
            [
                storage_preflight.ENTRYPOINT_NAME,
                "--expected-git-commit",
                expected_commit,
                "--audit-exact-pair",
                "--probe-id",
                PROBE_ID,
                "--expected-evidence-sha256",
                "a" * 64,
            ],
        ), redirect_stdout(output):
            storage_preflight.main()
        audit.assert_called_once_with(
            probe_id=PROBE_ID,
            expected_sha256="a" * 64,
            expected_git_commit=expected_commit,
        )
        self.assertEqual(
            output.getvalue(),
            storage_preflight._canonical_json_bytes(expected).decode(),
        )

    def test_cli_rejects_direct_execution_without_authenticated_runner(self) -> None:
        with mock.patch.object(
            sys,
            "_pic_control_plane_bootstrapped",
            False,
            create=True,
        ):
            with self.assertRaisesRegex(SystemExit, "through run_control_plane.py"):
                storage_preflight.main()

    def test_runner_allows_probe_only_as_authenticated_source_entrypoint(self) -> None:
        self.assertEqual(
            run_control_plane.SOURCE_ONLY_ENTRYPOINTS,
            {
                "capture_storage_preflight_evidence.py",
                "install_control_plane.py",
            },
        )
        self.assertIn(
            "capture_storage_preflight_evidence.py",
            run_control_plane.SOURCE_CONTROL_PLANE_FILES,
        )
        self.assertIn(
            "storage_preflight.schema.json",
            run_control_plane.SOURCE_CONTROL_PLANE_FILES,
        )
        self.assertNotIn(
            "capture_storage_preflight_evidence.py",
            run_control_plane.CONTROL_PLANE_FILES,
        )
        with self.assertRaisesRegex(ValueError, "source-only entrypoints"):
            run_control_plane._verify_source(
                Path("."),
                -1,
                "unreviewed.py",
                expected_git_commit="a" * 40,
            )

    def test_source_runner_requires_and_dispatches_exact_expected_git_commit(
        self,
    ) -> None:
        runner = CONTROL_PLANE_DIR / "run_control_plane.py"
        with mock.patch.object(
            run_control_plane.sys,
            "argv",
            [str(runner), "install_control_plane.py", "--help"],
        ), self.assertRaisesRegex(ValueError, "requires an expected Git commit"):
            run_control_plane.main()

        expected_commit = "a" * 40
        with mock.patch.object(
            run_control_plane.sys,
            "argv",
            [
                str(runner),
                "--expected-git-commit",
                expected_commit,
                "install_control_plane.py",
                "--help",
            ],
        ), mock.patch.object(
            run_control_plane,
            "_verify_source",
            return_value={"install_control_plane.py": b""},
        ) as verify_source, mock.patch.object(
            run_control_plane,
            "_execute_captured",
        ) as execute:
            run_control_plane.main()
        self.assertEqual(
            verify_source.call_args.kwargs["expected_git_commit"],
            expected_commit,
        )
        execute.assert_called_once()

    def test_schema_is_closed_and_binds_production_paths(self) -> None:
        schema = json.loads(
            (CONTROL_PLANE_DIR / storage_preflight.SCHEMA_NAME).read_text(
                encoding="utf-8"
            )
        )
        self.assertFalse(schema["additionalProperties"])
        self.assertEqual(schema["properties"]["schema_version"]["const"], 2)
        self.assertEqual(
            schema["properties"]["method"]["const"],
            storage_preflight.METHOD,
        )
        probes = schema["properties"]["probes"]["prefixItems"]
        self.assertEqual(
            schema["$defs"]["orionProbe"]["allOf"][1]["properties"]["path"]["const"],
            str(storage_preflight.AUTHORIZED_PIC_ROOT),
        )
        self.assertEqual(
            schema["$defs"]["projectHomeProbe"]["allOf"][1]["properties"]["path"][
                "const"
            ],
            str(storage_preflight.AUTHORIZED_PROJECT_HOME_ROOT),
        )
        self.assertEqual(len(probes), 2)


if __name__ == "__main__":
    unittest.main()
