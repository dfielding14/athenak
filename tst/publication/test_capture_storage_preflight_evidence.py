#!/usr/bin/env python3
"""Focused tests for authenticated fixed-root storage-preflight evidence."""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
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

        self.assertEqual(authenticate_source.calls, 3)
        self.assertEqual(orion_path.read_bytes(), project_home_path.read_bytes())
        self.assertFalse(orion_path.stat().st_mode & 0o222)
        self.assertFalse(project_home_path.stat().st_mode & 0o222)
        self._assert_no_probe_files()
        payload = orion_path.read_bytes()
        artifact = json.loads(payload)
        self.assertEqual(artifact["schema_version"], 1)
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

    def test_source_change_after_root_probes_fails_before_publication(self) -> None:
        changed = {**SOURCE_AUTHENTICATION, "git_commit": "e" * 40}
        authenticate_source = SourceAuthenticator([SOURCE_AUTHENTICATION, changed])
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
        real_publish = storage_preflight._publish_evidence
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
            "_publish_evidence",
            side_effect=replace_before_publish,
        ):
            with self.assertRaisesRegex(ValueError, "changed after its storage probe"):
                self._capture()

        self.assertFalse((self.orion_root / "policy").exists())
        self.assertFalse((moved / "policy").exists())
        self.assertFalse((self.project_home_root / "policy").exists())

    def test_second_publication_failure_rolls_back_first_artifact(self) -> None:
        real_publish = storage_preflight._publish_evidence

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
            "_publish_evidence",
            side_effect=fail_project_home,
        ):
            with self.assertRaisesRegex(OSError, "Project Home publication failure"):
                self._capture()

        orion_path, project_home_path = self._evidence_paths()
        self.assertFalse(orion_path.exists())
        self.assertFalse(project_home_path.exists())

    def test_existing_evidence_collision_is_not_replaced(self) -> None:
        orion_path, _ = self._evidence_paths()
        orion_path.parent.mkdir(parents=True)
        orion_path.write_bytes(b"existing\n")

        with self.assertRaises(FileExistsError):
            self._capture()

        self.assertEqual(orion_path.read_bytes(), b"existing\n")

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
        self.assertEqual(vars(parser.parse_args([])), {})
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parser.parse_args(["--orion-root", str(self.orion_root)])

        expected = {"policy": "fragment"}
        output = io.StringIO()
        with mock.patch.object(
            storage_preflight,
            "capture_storage_preflight_evidence",
            return_value=expected,
        ), mock.patch.object(sys, "_pic_control_plane_bootstrapped", True, create=True), (
            mock.patch.object(sys, "argv", [storage_preflight.ENTRYPOINT_NAME])
        ), redirect_stdout(output):
            storage_preflight.main()
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
            run_control_plane._verify_source(Path("."), -1, "unreviewed.py")

    def test_schema_is_closed_and_binds_production_paths(self) -> None:
        schema = json.loads(
            (CONTROL_PLANE_DIR / storage_preflight.SCHEMA_NAME).read_text(
                encoding="utf-8"
            )
        )
        self.assertFalse(schema["additionalProperties"])
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
