#!/usr/bin/env python3
"""Focused tests for the operator-only pre-policy-promotion attestation helper."""

from __future__ import annotations

from datetime import datetime, timezone
import errno
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


SCRIPT_DIR = Path(__file__).absolute().parent
sys.path.insert(0, str(SCRIPT_DIR))
import capture_frontier_pre_policy_promotion_attestation as capture_attestation


FIXED_TIME = datetime(2026, 6, 1, 15, 16, 17, tzinfo=timezone.utc)
CONTROL_PLANE_VERSION = "a" * 64


class ReadOnlyRunner:
    """Provide deterministic ps and squeue output while recording exact argv."""

    def __init__(self, queues: list[str] | None = None) -> None:
        self.queues = list(queues or [""])
        self.calls: list[list[str]] = []

    def __call__(self, argv: list[str]) -> subprocess.CompletedProcess[str]:
        self.calls.append(list(argv))
        if argv[0] == capture_attestation.TRUSTED_PS:
            stdout = "100 1 S /usr/bin/python3 operator-review.py\n"
        elif argv[0] == capture_attestation.TRUSTED_SQUEUE:
            stdout = self.queues.pop(0) if len(self.queues) > 1 else self.queues[0]
        else:
            raise AssertionError(f"Unexpected subprocess: {argv}")
        return subprocess.CompletedProcess(argv, 0, stdout, "")


class CaptureFrontierPrePolicyPromotionAttestationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.pic_root = self.base / "orion"
        self.project_home_root = self.base / "project-home"
        self.archive_root = self.pic_root / "operator_attestations"
        (self.pic_root / "ledger").mkdir(parents=True)
        (self.project_home_root / "ledger").mkdir(parents=True)
        self.archive_root.mkdir()
        self.ledger = self.pic_root / "ledger" / "node_hours.jsonl"
        self.receipts = self.pic_root / "ledger" / "mirror_receipts.jsonl"
        self.mirror = self.project_home_root / "ledger" / "node_hours.jsonl"
        self._write_ledger_lines(1)
        self.records = [
            {
                "event_type": "genesis",
                "event_sha256": "b" * 64,
                "state": "initialized",
            }
        ]

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write_ledger_lines(self, count: int) -> None:
        payload = b"{}\n" * count
        for path in [self.ledger, self.receipts, self.mirror]:
            path.write_bytes(payload)

    def _validator(
        self, paths: capture_attestation.LedgerPaths
    ) -> list[dict[str, object]]:
        self.assertEqual(paths.ledger, self.ledger)
        self.assertEqual(paths.receipts, self.receipts)
        self.assertEqual(paths.mirror, self.mirror)
        return [dict(record) for record in self.records]

    def _capture(
        self,
        *,
        runner: ReadOnlyRunner,
        authorization_id: str = "reviewed-authorization",
        phase: str = capture_attestation.PHASE,
        ledger_validator=None,
    ) -> Path:
        return capture_attestation.capture(
            archive_root=self.archive_root,
            authorization_id=authorization_id,
            control_plane_version=CONTROL_PLANE_VERSION,
            phase=phase,
            pic_root=self.pic_root,
            project_home_root=self.project_home_root,
            runner=runner,
            now=lambda: FIXED_TIME,
            user="dfielding",
            ledger_validator=ledger_validator or self._validator,
        )

    def _seal(
        self,
        staging: Path,
        *,
        runner: ReadOnlyRunner,
        attest_reviewed: bool = True,
        ledger_validator=None,
    ) -> Path:
        return capture_attestation.seal(
            staging,
            attest_reviewed=attest_reviewed,
            runner=runner,
            now=lambda: FIXED_TIME,
            ledger_validator=ledger_validator or self._validator,
        )

    def test_success_uses_only_absolute_read_only_scheduler_queries(self) -> None:
        runner = ReadOnlyRunner(["", ""])
        staging = self._capture(runner=runner)
        self.assertTrue(staging.name.startswith("."))
        final = self._seal(staging, runner=runner)

        self.assertTrue(final.is_dir())
        self.assertFalse(staging.exists())
        self.assertEqual(
            runner.calls,
            [
                [
                    "/usr/bin/ps",
                    "-u",
                    "dfielding",
                    "-o",
                    "pid=,ppid=,state=,args=",
                ],
                [
                    "/usr/bin/squeue",
                    "-u",
                    "dfielding",
                    "-h",
                    "-o",
                    "%i|%a|%P|%q|%T|%j|%k",
                ],
                [
                    "/usr/bin/squeue",
                    "-u",
                    "dfielding",
                    "-h",
                    "-o",
                    "%i|%a|%P|%q|%T|%j|%k",
                ],
            ],
        )
        attestation = json.loads((final / "attestation.json").read_bytes())
        self.assertEqual(attestation["schema_version"], 1)
        self.assertEqual(attestation["phase"], "pre_policy_promotion")
        self.assertEqual(
            attestation["operator_statement"], capture_attestation.OPERATOR_STATEMENT
        )
        self.assertEqual(
            final.name,
            "20260601T151617Z-reviewed-authorization-pre_policy_promotion",
        )
        self.assertEqual(attestation["pending_submission_marker"]["value"], "absent")
        self.assertEqual(
            attestation["pending_manual_accounting_marker"]["value"], "absent"
        )
        self.assertEqual(
            attestation["validated_mirrored_ledger_state"]["state"][
                "active_reservation_count"
            ],
            0,
        )
        self.assertTrue((final / "capture_queue_snapshot.txt").is_file())

    def test_capture_and_seal_encode_selected_execution_phase(self) -> None:
        for phase in ["pre_manifest", "pre_submit_wrapper"]:
            with self.subTest(phase=phase):
                runner = ReadOnlyRunner(["", ""])
                staging = self._capture(
                    runner=runner,
                    authorization_id=f"reviewed-{phase}",
                    phase=phase,
                )
                expected_name = f"20260601T151617Z-reviewed-{phase}-{phase}"
                self.assertTrue(staging.name.startswith(f".{expected_name}.staging-"))
                metadata = json.loads(
                    (staging / capture_attestation.METADATA_FILENAME).read_bytes()
                )
                self.assertEqual(metadata["phase"], phase)

                final = self._seal(staging, runner=runner)
                attestation = json.loads((final / "attestation.json").read_bytes())
                self.assertEqual(final.name, expected_name)
                self.assertEqual(attestation["phase"], phase)
                self.assertEqual(
                    attestation["operator_statement"],
                    capture_attestation.OPERATOR_STATEMENTS[phase],
                )
                self.assertEqual(
                    attestation["pending_submission_marker"]["value"], "absent"
                )
                self.assertEqual(
                    attestation["validated_mirrored_ledger_state"]["state"][
                        "active_reservation_count"
                    ],
                    0,
                )

    def test_capture_cli_phase_choices_default_to_pre_policy_promotion(self) -> None:
        required = [
            "capture",
            "--authorization-id",
            "reviewed-authorization",
            "--control-plane-version",
            CONTROL_PLANE_VERSION,
        ]
        parser = capture_attestation._parser()
        self.assertEqual(parser.parse_args(required).phase, "pre_policy_promotion")
        for phase in capture_attestation.PHASES:
            with self.subTest(phase=phase):
                self.assertEqual(
                    parser.parse_args([*required, "--phase", phase]).phase,
                    phase,
                )

    def test_capture_cli_rejects_invalid_phase(self) -> None:
        with mock.patch.object(sys, "stderr"):
            with self.assertRaises(SystemExit) as raised:
                capture_attestation._parser().parse_args(
                    [
                        "capture",
                        "--authorization-id",
                        "reviewed-authorization",
                        "--control-plane-version",
                        CONTROL_PLANE_VERSION,
                        "--phase",
                        "post_submit_wrapper",
                    ]
                )
        self.assertEqual(raised.exception.code, 2)

    def test_capture_rejects_invalid_phase_before_queries(self) -> None:
        runner = ReadOnlyRunner([""])
        with self.assertRaisesRegex(ValueError, "Phase must be one of"):
            self._capture(runner=runner, phase="post_submit_wrapper")
        self.assertEqual(runner.calls, [])
        self.assertEqual(list(self.archive_root.iterdir()), [])

    def test_seal_accepts_legacy_default_phase_staging_metadata(self) -> None:
        runner = ReadOnlyRunner(["", ""])
        staging = self._capture(runner=runner)
        metadata_path = staging / capture_attestation.METADATA_FILENAME
        metadata = json.loads(metadata_path.read_bytes())
        del metadata["phase"]
        metadata_path.write_text(json.dumps(metadata) + "\n", encoding="utf-8")

        final = self._seal(staging, runner=runner)
        attestation = json.loads((final / "attestation.json").read_bytes())
        self.assertEqual(attestation["phase"], "pre_policy_promotion")
        self.assertTrue(final.name.endswith("-pre_policy_promotion"))

    def test_seal_rejects_nonempty_queue(self) -> None:
        runner = ReadOnlyRunner(["", "123|ast207|batch|debug|RUNNING|job|comment\n"])
        staging = self._capture(runner=runner)
        with self.assertRaisesRegex(ValueError, "non-empty queue"):
            self._seal(staging, runner=runner)
        self.assertTrue(staging.is_dir())

    def test_seal_rejects_each_pending_marker(self) -> None:
        markers = [
            self.pic_root / "ledger" / "pending_submission.json",
            self.pic_root / "ledger" / "pending_manual_accounting.json",
            self.project_home_root / "ledger" / "pending_manual_accounting.json",
        ]
        for index, marker in enumerate(markers):
            with self.subTest(marker=marker):
                runner = ReadOnlyRunner(["", ""])
                staging = self._capture(
                    runner=runner, authorization_id=f"reviewed-marker-{index}"
                )
                marker.write_text("pending\n", encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "blocked by pending"):
                    self._seal(staging, runner=runner)
                marker.unlink()

    def test_seal_rejects_outstanding_reservation(self) -> None:
        self._write_ledger_lines(2)
        self.records.append(
            {
                "event_type": "reservation",
                "event_sha256": "c" * 64,
                "reservation_id": "reservation-one",
                "state": "reserved",
            }
        )
        runner = ReadOnlyRunner(["", ""])
        staging = self._capture(runner=runner)
        with self.assertRaisesRegex(ValueError, "outstanding reservation"):
            self._seal(staging, runner=runner)

    def test_capture_rejects_mirrored_ledger_divergence(self) -> None:
        def divergent_validator(
            paths: capture_attestation.LedgerPaths,
        ) -> list[dict[str, object]]:
            raise ValueError("Local and mirrored PIC ledger records differ")

        runner = ReadOnlyRunner([""])
        with self.assertRaisesRegex(ValueError, "ledger records differ"):
            self._capture(runner=runner, ledger_validator=divergent_validator)
        self.assertEqual(list(self.archive_root.iterdir()), [])

    def test_capture_rejects_unsafe_authorization_ids_before_queries(self) -> None:
        for authorization_id in ["../escape", "has space", "-leading", "trailing-", ""]:
            with self.subTest(authorization_id=authorization_id):
                runner = ReadOnlyRunner([""])
                with self.assertRaisesRegex(ValueError, "Authorization ID"):
                    self._capture(runner=runner, authorization_id=authorization_id)
                self.assertEqual(runner.calls, [])

    def test_seal_requires_explicit_review_attestation(self) -> None:
        runner = ReadOnlyRunner([""])
        staging = self._capture(runner=runner)
        with self.assertRaisesRegex(ValueError, "--attest-reviewed"):
            self._seal(staging, runner=runner, attest_reviewed=False)

    def test_seal_publishes_atomically_with_read_only_tree(self) -> None:
        runner = ReadOnlyRunner(["", ""])
        staging = self._capture(runner=runner)
        original_rename = capture_attestation._rename_no_replace
        observed = {"called": False}

        def inspect_then_rename(source: Path, destination: Path) -> None:
            observed["called"] = True
            self.assertFalse(destination.exists())
            self.assertEqual(stat.S_IMODE(os.lstat(source).st_mode), 0o500)
            for member in source.iterdir():
                self.assertTrue(member.is_file())
                self.assertEqual(stat.S_IMODE(os.lstat(member).st_mode), 0o400)
            self.assertTrue((source / "attestation.json").is_file())
            self.assertFalse((source / capture_attestation.METADATA_FILENAME).exists())
            original_rename(source, destination)

        with mock.patch.object(
            capture_attestation, "_rename_no_replace", side_effect=inspect_then_rename
        ):
            final = self._seal(staging, runner=runner)
        self.assertTrue(observed["called"])
        self.assertTrue(final.is_dir())
        self.assertEqual(stat.S_IMODE(os.lstat(final).st_mode), 0o500)

    def test_lustre_rename_fallback_preserves_identity_and_rejects_collision(self) -> None:
        source = self.archive_root / "source"
        destination = self.archive_root / "destination"
        source.mkdir()
        with mock.patch.object(
            capture_attestation,
            "_renameat2_no_replace",
            side_effect=OSError(errno.EINVAL, os.strerror(errno.EINVAL)),
        ):
            capture_attestation._rename_no_replace(source, destination)
        self.assertFalse(source.exists())
        self.assertTrue(destination.is_dir())

        collision_source = self.archive_root / "collision-source"
        collision_destination = self.archive_root / "collision-destination"
        collision_source.mkdir()
        collision_destination.mkdir()
        with mock.patch.object(
            capture_attestation,
            "_renameat2_no_replace",
            side_effect=OSError(errno.EINVAL, os.strerror(errno.EINVAL)),
        ):
            with self.assertRaisesRegex(ValueError, "collides"):
                capture_attestation._rename_no_replace(
                    collision_source, collision_destination
                )
        self.assertTrue(collision_source.is_dir())
        self.assertTrue(collision_destination.is_dir())

    def test_failed_final_rename_can_publish_sealed_staging(self) -> None:
        runner = ReadOnlyRunner(["", ""])
        staging = self._capture(runner=runner)
        with mock.patch.object(
            capture_attestation,
            "_rename_no_replace",
            side_effect=OSError(errno.EIO, os.strerror(errno.EIO)),
        ):
            with self.assertRaises(OSError):
                self._seal(staging, runner=runner)
        self.assertEqual(stat.S_IMODE(os.lstat(staging).st_mode), 0o500)
        self.assertFalse((staging / capture_attestation.METADATA_FILENAME).exists())

        final = capture_attestation.publish_sealed_staging(staging)
        self.assertFalse(staging.exists())
        self.assertTrue((final / "attestation.json").is_file())

    def test_capture_rejects_archive_root_symlink_alias(self) -> None:
        alias = self.base / "archive-alias"
        alias.symlink_to(self.archive_root, target_is_directory=True)
        runner = ReadOnlyRunner([""])
        with self.assertRaisesRegex(ValueError, "symlink alias"):
            capture_attestation.capture(
                archive_root=alias,
                authorization_id="reviewed-alias",
                control_plane_version=CONTROL_PLANE_VERSION,
                pic_root=self.pic_root,
                project_home_root=self.project_home_root,
                runner=runner,
                now=lambda: FIXED_TIME,
                user="dfielding",
                ledger_validator=self._validator,
            )
        self.assertEqual(runner.calls, [])

    def test_exact_trusted_project_home_mount_alias_is_accepted(self) -> None:
        target = self.base / "trusted-project-home-target"
        target.mkdir()
        alias = self.base / "trusted-project-home-alias"
        alias.symlink_to(target, target_is_directory=True)
        self.assertEqual(
            capture_attestation._canonical_existing_directory(
                alias,
                label="Project Home root",
                trusted_lexical_alias=alias,
            ),
            alias,
        )
        second_alias = self.base / "untrusted-project-home-alias"
        second_alias.symlink_to(target, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "symlink alias"):
            capture_attestation._canonical_existing_directory(
                second_alias,
                label="Project Home root",
                trusted_lexical_alias=alias,
            )

    def test_seal_rejects_final_directory_collision(self) -> None:
        runner = ReadOnlyRunner([""])
        staging = self._capture(runner=runner)
        metadata = json.loads((staging / capture_attestation.METADATA_FILENAME).read_bytes())
        (self.archive_root / metadata["final_directory_name"]).mkdir()
        with self.assertRaisesRegex(ValueError, "collides"):
            self._seal(staging, runner=runner)

    def test_seal_rejects_staged_snapshot_symlink(self) -> None:
        runner = ReadOnlyRunner([""])
        staging = self._capture(runner=runner)
        queue = staging / "queue_snapshot.txt"
        queue.unlink()
        queue.symlink_to(self.ledger)
        with self.assertRaisesRegex(ValueError, "not a regular file"):
            self._seal(staging, runner=runner)

    def test_seal_rejects_unknown_staged_file(self) -> None:
        runner = ReadOnlyRunner([""])
        staging = self._capture(runner=runner)
        (staging / "unexpected.txt").write_text("unexpected\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "members differ"):
            self._seal(staging, runner=runner)


if __name__ == "__main__":
    unittest.main()
