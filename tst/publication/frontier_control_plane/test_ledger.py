#!/usr/bin/env python3
"""Unit tests for the mirrored Frontier PIC ledger."""

from __future__ import annotations

import csv
import fcntl
import json
import os
from pathlib import Path
import shutil
import stat
import tempfile
import unittest
from unittest.mock import patch

from control_plane_common import AUTHORIZED_PIC_ROOT, stable_serialization_anchor
from ledger import accounting, append_primary_event, initialize_ledger
from ledger import genesis_anchor_paths, migrate_existing_genesis_anchors
from ledger import ledger_lock
from ledger import repair_mirrored_state, validate_mirrored_state
from ledger import validate_primary_chain, validate_receipts, write_csv
from ledger import validated_read_only_mirrored_state_snapshot


class LedgerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        self.ledger = root / "orion" / "ledger" / "node_hours.jsonl"
        self.csv = root / "orion" / "ledger" / "node_hours.csv"
        self.receipts = root / "orion" / "ledger" / "mirror_receipts.jsonl"
        self.mirror = root / "project_home" / "ledger" / "node_hours.jsonl"
        initialize_ledger(
            self.ledger,
            self.csv,
            self.receipts,
            self.mirror,
            mirror_transport="filesystem_copy",
            notes="temporary test genesis",
            control_plane_version="temporary-test-control-plane",
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def append(self, event: dict[str, object]) -> dict[str, object]:
        return append_primary_event(
            self.ledger,
            self.csv,
            self.receipts,
            self.mirror,
            event,
            mirror_transport="filesystem_copy",
        )

    def _manual_accounting_event(self) -> dict[str, object]:
        return {
            "event_type": "manual_allocation_reconciliation",
            "job_id": "4746332",
            "control_plane_version": "a" * 64,
            "reconciled_by_control_plane_version": "a" * 64,
            "manual_accounting_authorization_id": "reviewed-allocation",
            "manual_accounting_authorization_path": "/tmp/reviewed-allocation.json",
            "manual_accounting_project_home_authorization_path": (
                "/tmp/project-home/reviewed-allocation.json"
            ),
            "manual_accounting_authorization_sha256": "b" * 64,
            "accounting_scope": "manual_direct_srun_accounting_only",
            "scientific_evidence_eligible": False,
            "active_policy_sha256": "c" * 64,
            "active_promotion_sha256": "d" * 64,
            "partition": "batch",
            "qos": "normal",
            "scheduler_reported_allocated_nodes": 1,
            "billed_nodes": 1,
            "elapsed_seconds": 5,
            "reconciled": True,
            "consumed_node_hours": 1.0 / 720.0,
            "cumulative_consumed_node_hours": 1.0 / 720.0,
            "state": "FAILED",
            "notes": (
                "Reviewed direct-srun accounting only; "
                "ineligible for scientific evidence."
            ),
        }

    def _drop_last_line(self, path: Path) -> None:
        lines = path.read_text(encoding="utf-8").splitlines()
        path.write_text("".join(line + "\n" for line in lines[:-1]), encoding="utf-8")

    def test_reserve_attach_reconcile_and_csv_projection(self) -> None:
        common = {
            "reservation_id": "reservation-1",
            "submission_id": "submission-1",
            "submission_scope": "registered_science",
            "registered_science_authorization_id": "f1-clean-gyro-v1",
            "clean_candidate_manifest_sha256": "a" * 64,
            "requested_nodes": 2,
            "requested_walltime": "00:30:00",
            "reserved_node_hours": 1.0,
        }
        self.append({**common, "event_type": "reservation", "state": "reserved"})
        self.append(
            {
                **common,
                "event_type": "job_id_attached",
                "job_id": "1234",
                "state": "submitted",
            }
        )
        self.assertEqual(accounting(validate_primary_chain(self.ledger)), {
            "cumulative_consumed_node_hours": 0.0,
            "currently_reserved_node_hours": 1.0,
        })
        self.append(
            {
                **common,
                "event_type": "reconciliation",
                "job_id": "1234",
                "state": "COMPLETED",
                "reconciled": True,
                "scheduler_reported_allocated_nodes": 2,
                "billed_nodes": 2,
                "elapsed_seconds": 600,
                "consumed_node_hours": 1.0 / 3.0,
                "cumulative_consumed_node_hours": 1.0 / 3.0,
            }
        )
        totals = accounting(validate_primary_chain(self.ledger))
        self.assertAlmostEqual(totals["cumulative_consumed_node_hours"], 1.0 / 3.0)
        self.assertEqual(totals["currently_reserved_node_hours"], 0.0)
        with self.csv.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 4)
        self.assertTrue(rows[-1]["mirror_ack_sha256"])
        self.assertEqual(rows[-1]["submission_scope"], "registered_science")
        self.assertEqual(
            rows[-1]["registered_science_authorization_id"], "f1-clean-gyro-v1"
        )
        self.assertEqual(rows[-1]["clean_candidate_manifest_sha256"], "a" * 64)

    def test_manual_direct_srun_reconciliation_counts_toward_consumed_budget(
        self,
    ) -> None:
        self.append(self._manual_accounting_event())
        totals = accounting(validate_primary_chain(self.ledger))
        self.assertAlmostEqual(
            totals["cumulative_consumed_node_hours"], 1.0 / 720.0
        )
        self.assertEqual(totals["currently_reserved_node_hours"], 0.0)

    def test_manual_direct_srun_reconciliation_rejects_under_bound_event(
        self,
    ) -> None:
        with self.assertRaisesRegex(ValueError, "semantics are invalid"):
            self.append(
                {
                    "event_type": "manual_allocation_reconciliation",
                    "job_id": "4746332",
                    "accounting_scope": "manual_direct_srun_accounting_only",
                    "scientific_evidence_eligible": False,
                    "reconciled": True,
                    "consumed_node_hours": 1.0 / 720.0,
                }
            )

    def test_manual_direct_srun_reconciliation_rejects_false_cumulative_usage(
        self,
    ) -> None:
        event = self._manual_accounting_event()
        event["cumulative_consumed_node_hours"] = 999.0
        with self.assertRaisesRegex(ValueError, "cumulative usage differs"):
            self.append(event)

    def test_mirror_head_divergence_fails_closed(self) -> None:
        with self.mirror.open("a", encoding="utf-8") as stream:
            stream.write("{}\n")
        with self.assertRaises(ValueError):
            self.append({"event_type": "should_not_append"})

    def test_receipts_are_valid_and_non_recursive(self) -> None:
        primary = validate_primary_chain(self.ledger)
        receipts = validate_receipts(
            self.receipts,
            primary,
            mirror_jsonl=self.mirror,
            mirror_transport="filesystem_copy",
        )
        self.assertEqual(len(primary), 1)
        self.assertEqual(len(receipts), 1)
        self.assertNotIn("sequence_number", receipts[0])

    def test_read_only_snapshot_rejects_receipt_drift_during_use(self) -> None:
        with self.assertRaisesRegex(ValueError, "changed during snapshot use"):
            with validated_read_only_mirrored_state_snapshot(
                self.ledger,
                self.receipts,
                self.mirror,
                ledger_root=self.ledger.parent.parent,
                receipts_root=self.receipts.parent.parent,
                mirror_root=self.mirror.parent.parent,
            ):
                self.receipts.write_text(
                    self.receipts.read_text(encoding="utf-8") + "\n",
                    encoding="utf-8",
                )

    def test_read_only_snapshot_rejects_mirror_parent_replacement(self) -> None:
        parent = self.mirror.parent
        original = parent.with_name("ledger-original")
        try:
            with self.assertRaisesRegex(ValueError, "parent path changed"):
                with validated_read_only_mirrored_state_snapshot(
                    self.ledger,
                    self.receipts,
                    self.mirror,
                    ledger_root=self.ledger.parent.parent,
                    receipts_root=self.receipts.parent.parent,
                    mirror_root=self.mirror.parent.parent,
                ):
                    parent.rename(original)
                    parent.mkdir()
        finally:
            if parent.exists():
                shutil.rmtree(parent)
            if original.exists():
                original.rename(parent)

    def test_read_only_snapshot_rejects_byte_identical_receipt_replacement(self) -> None:
        replacement = self.receipts.with_name("replacement.jsonl")
        replacement.write_bytes(self.receipts.read_bytes())
        with self.assertRaisesRegex(ValueError, "snapshot path changed"):
            with validated_read_only_mirrored_state_snapshot(
                self.ledger,
                self.receipts,
                self.mirror,
                ledger_root=self.ledger.parent.parent,
                receipts_root=self.receipts.parent.parent,
                mirror_root=self.mirror.parent.parent,
            ):
                replacement.replace(self.receipts)

    def test_read_only_snapshot_rejects_genesis_anchor_replacement(self) -> None:
        anchor, _ = genesis_anchor_paths(self.ledger, self.mirror)
        replacement = anchor.with_name("replacement-anchor.json")
        replacement.write_bytes(anchor.read_bytes())
        replacement.chmod(0o444)
        with self.assertRaisesRegex(ValueError, "snapshot path changed"):
            with validated_read_only_mirrored_state_snapshot(
                self.ledger,
                self.receipts,
                self.mirror,
                ledger_root=self.ledger.parent.parent,
                receipts_root=self.receipts.parent.parent,
                mirror_root=self.mirror.parent.parent,
            ):
                replacement.replace(anchor)

    def test_read_only_snapshot_validates_pinned_bytes_not_transient_paths(self) -> None:
        paths = [self.ledger, self.csv, self.receipts, self.mirror]
        original = {path: path.read_bytes() for path in paths}
        self.append(
            {
                "event_type": "reservation",
                "reservation_id": "transient-alternate-reservation",
                "state": "reserved",
                "reserved_node_hours": 0.0,
            }
        )
        alternate = {path: path.read_bytes() for path in paths}
        for path, data in original.items():
            path.write_bytes(data)

        from ledger import _validate_mirrored_state_bytes

        reopened_lengths = []

        def validate_pinned_bytes(
            state: dict[Path, bytes], **kwargs: object
        ) -> list[dict[str, object]]:
            for path in [self.ledger, self.receipts, self.mirror]:
                path.write_bytes(alternate[path])
            try:
                reopened_lengths.append(
                    len(validate_mirrored_state(self.ledger, self.receipts, self.mirror))
                )
                return _validate_mirrored_state_bytes(state, **kwargs)
            finally:
                for path, data in original.items():
                    path.write_bytes(data)

        with patch(
            "ledger._validate_mirrored_state_bytes", side_effect=validate_pinned_bytes
        ):
            with validated_read_only_mirrored_state_snapshot(
                self.ledger,
                self.receipts,
                self.mirror,
                ledger_root=self.ledger.parent.parent,
                receipts_root=self.receipts.parent.parent,
                mirror_root=self.mirror.parent.parent,
            ) as records:
                self.assertEqual(len(records), 1)
        self.assertEqual(reopened_lengths, [2])

    def test_read_only_snapshot_missing_parent_fails_without_creation(self) -> None:
        missing = self.mirror.parent.parent / "missing" / "node_hours.jsonl"
        with self.assertRaises(FileNotFoundError):
            with validated_read_only_mirrored_state_snapshot(
                self.ledger,
                self.receipts,
                missing,
                ledger_root=self.ledger.parent.parent,
                receipts_root=self.receipts.parent.parent,
                mirror_root=missing.parent.parent,
            ):
                pass
        self.assertFalse(missing.parent.exists())

    def test_read_only_snapshot_rechecks_rooted_traversal_for_nested_pin(self) -> None:
        from ledger import _ledger_state_paths, _pinned_parent_directories

        paths = _ledger_state_paths(self.ledger, None, self.receipts, self.mirror)
        unrelated_root = self.ledger.parent.parent / "unrelated"
        unrelated_root.mkdir()
        with _pinned_parent_directories(paths, create_missing=False):
            with self.assertRaisesRegex(ValueError, "outside trusted lexical root"):
                with validated_read_only_mirrored_state_snapshot(
                    self.ledger,
                    self.receipts,
                    self.mirror,
                    ledger_root=unrelated_root,
                    receipts_root=unrelated_root,
                    mirror_root=self.mirror.parent.parent,
                ):
                    pass

    def test_self_consistent_receipt_provenance_rewrite_is_rejected(self) -> None:
        receipt = json.loads(self.receipts.read_text(encoding="utf-8"))
        receipt["mirror_destination"] = "/tmp/substituted-mirror.jsonl"
        receipt["mirror_transport"] = "substituted_transport"
        from ledger import canonical_json, record_sha256
        receipt["mirror_ack_sha256"] = record_sha256(
            receipt, "mirror_ack_sha256"
        )
        self.receipts.write_text(canonical_json(receipt) + "\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            validate_mirrored_state(self.ledger, self.receipts, self.mirror)

    def test_non_genesis_append_rejects_empty_ledgers(self) -> None:
        root = Path(self.temporary.name) / "empty"
        with self.assertRaises(ValueError):
            append_primary_event(
                root / "orion" / "ledger" / "node_hours.jsonl",
                root / "orion" / "ledger" / "node_hours.csv",
                root / "orion" / "ledger" / "mirror_receipts.jsonl",
                root / "project_home" / "ledger" / "node_hours.jsonl",
                {"event_type": "reservation", "state": "reserved"},
                mirror_transport="filesystem_copy",
            )

    def test_corrupt_receipt_rejects_append_before_primary_mutation(self) -> None:
        with self.receipts.open("a", encoding="utf-8") as stream:
            stream.write("{}\n")
        with self.assertRaises(ValueError):
            self.append({"event_type": "must_not_append"})
        self.assertEqual(len(validate_primary_chain(self.ledger)), 1)
        self.assertEqual(len(validate_primary_chain(self.mirror)), 1)

    def test_missing_receipt_rejects_validation(self) -> None:
        self._drop_last_line(self.receipts)
        with self.assertRaises(ValueError):
            validate_mirrored_state(self.ledger, self.receipts, self.mirror)

    def test_duplicate_receipt_rejects_validation(self) -> None:
        line = self.receipts.read_text(encoding="utf-8").splitlines()[0]
        with self.receipts.open("a", encoding="utf-8") as stream:
            stream.write(line + "\n")
        with self.assertRaises(ValueError):
            validate_mirrored_state(self.ledger, self.receipts, self.mirror)

    def test_duplicate_json_key_rejects_validation(self) -> None:
        self.ledger.write_text('{"sequence_number": 0, "sequence_number": 1}\n',
                               encoding="utf-8")
        with self.assertRaises(ValueError):
            validate_primary_chain(self.ledger)

    def test_repair_copies_missing_mirror_suffix_and_receipt(self) -> None:
        self.append({"event_type": "reservation", "state": "reserved"})
        self._drop_last_line(self.mirror)
        self._drop_last_line(self.receipts)
        result = repair_mirrored_state(
            self.ledger,
            self.csv,
            self.receipts,
            self.mirror,
            mirror_transport="filesystem_copy",
        )
        self.assertEqual(
            result,
            {"appended_mirror_records": 1, "appended_receipts": 1},
        )
        self.assertEqual(
            validate_primary_chain(self.ledger),
            validate_mirrored_state(self.ledger, self.receipts, self.mirror),
        )

    def test_repair_appends_missing_receipt_only(self) -> None:
        self.append({"event_type": "reservation", "state": "reserved"})
        self._drop_last_line(self.receipts)
        result = repair_mirrored_state(
            self.ledger,
            self.csv,
            self.receipts,
            self.mirror,
            mirror_transport="filesystem_copy",
        )
        self.assertEqual(
            result,
            {"appended_mirror_records": 0, "appended_receipts": 1},
        )
        validate_mirrored_state(self.ledger, self.receipts, self.mirror)

    def test_repair_rejects_forged_primary_before_mutating_mirror(self) -> None:
        from ledger import canonical_json, record_sha256

        forged = json.loads(self.ledger.read_text(encoding="utf-8"))
        forged["notes"] = "forged replacement genesis"
        forged["event_sha256"] = record_sha256(forged, "event_sha256")
        self.ledger.write_text(canonical_json(forged) + "\n", encoding="utf-8")
        self.mirror.write_bytes(b"")
        self.receipts.write_bytes(b"")
        with self.assertRaises(ValueError):
            repair_mirrored_state(
                self.ledger,
                self.csv,
                self.receipts,
                self.mirror,
                mirror_transport="filesystem_copy",
            )
        self.assertEqual(self.mirror.read_bytes(), b"")
        self.assertEqual(self.receipts.read_bytes(), b"")

    def test_repair_rejects_missing_anchor_before_mirror_preflight(self) -> None:
        local_anchor, _ = genesis_anchor_paths(self.ledger, self.mirror)
        local_anchor.chmod(0o600)
        local_anchor.unlink()
        with patch("ledger.mirror_preflight") as preflight:
            with self.assertRaises(ValueError):
                repair_mirrored_state(
                    self.ledger,
                    self.csv,
                    self.receipts,
                    self.mirror,
                    mirror_transport="filesystem_copy",
                )
        preflight.assert_not_called()

    def test_append_rejects_missing_anchor_before_mirror_preflight(self) -> None:
        local_anchor, _ = genesis_anchor_paths(self.ledger, self.mirror)
        local_anchor.chmod(0o600)
        local_anchor.unlink()
        with patch("ledger.mirror_preflight") as preflight:
            with self.assertRaises(ValueError):
                self.append({"event_type": "must_not_append"})
        preflight.assert_not_called()

    def test_lock_symlink_alias_rejects_before_external_write(self) -> None:
        lock_path = self.ledger.with_suffix(self.ledger.suffix + ".lock")
        lock_path.unlink()
        outside = Path(self.temporary.name) / "outside-lock"
        outside.write_text("preserve external lock target\n", encoding="utf-8")
        lock_path.symlink_to(outside)
        with self.assertRaises(ValueError):
            self.append({"event_type": "must_not_append"})
        self.assertEqual(
            outside.read_text(encoding="utf-8"),
            "preserve external lock target\n",
        )

    def test_replaced_lock_entries_fail_closed_while_held(self) -> None:
        paths = [
            self.ledger.parent.parent / ".ledger.lock",
            self.ledger.with_suffix(self.ledger.suffix + ".lock"),
        ]
        for path in paths:
            with self.subTest(path=path):
                with self.assertRaisesRegex(ValueError, "lock path changed"):
                    with ledger_lock(self.ledger):
                        path.unlink()
                        path.write_text("replacement lock\n", encoding="utf-8")

    def test_production_serialization_anchor_is_outside_replaceable_pic_root(self) -> None:
        self.assertEqual(
            stable_serialization_anchor(AUTHORIZED_PIC_ROOT),
            Path("/lustre/orion/ast207"),
        )
        self.assertEqual(
            stable_serialization_anchor(self.ledger.parent.parent),
            Path(self.temporary.name),
        )

    def test_ledger_lock_holds_stable_serialization_anchor(self) -> None:
        anchor = stable_serialization_anchor(self.ledger.parent.parent)
        with ledger_lock(self.ledger):
            descriptor = os.open(anchor, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
            try:
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            finally:
                os.close(descriptor)

    def test_mirror_parent_swap_fails_closed_while_held(self) -> None:
        mirror_parent = self.mirror.parent
        displaced = mirror_parent.with_name("displaced-mirror")
        with self.assertRaisesRegex(ValueError, "path changed"):
            with ledger_lock(self.ledger, self.mirror):
                mirror_parent.rename(displaced)
                mirror_parent.mkdir()
        self.assertEqual(list(mirror_parent.iterdir()), [])

    def test_ledger_parent_swap_fails_without_writing_replacement(self) -> None:
        import ledger

        ledger_parent = self.ledger.parent
        displaced = ledger_parent.with_name("displaced-ledger")
        real_append = ledger._append_jsonl
        swapped = False

        def swap_parent_then_append(path: Path, record: dict[str, object]) -> None:
            nonlocal swapped
            if path == self.ledger and not swapped:
                swapped = True
                ledger_parent.rename(displaced)
                ledger_parent.mkdir()
            real_append(path, record)

        with patch("ledger._append_jsonl", side_effect=swap_parent_then_append):
            with self.assertRaises(ValueError):
                self.append({"event_type": "must_fail_closed"})
        self.assertTrue(swapped)
        self.assertEqual(list(ledger_parent.iterdir()), [])
        self.assertEqual(
            len(validate_primary_chain(displaced / self.ledger.name)),
            2,
        )

    def test_csv_output_symlink_alias_rejects_before_external_write(self) -> None:
        outside = Path(self.temporary.name) / "outside-csv-output"
        outside.write_text("preserve external CSV target\n", encoding="utf-8")
        self.csv.unlink()
        self.csv.symlink_to(outside)
        with self.assertRaises(ValueError):
            write_csv(
                self.ledger,
                self.receipts,
                self.csv,
                mirror_jsonl=self.mirror,
                mirror_transport="filesystem_copy",
            )
        self.assertEqual(
            outside.read_text(encoding="utf-8"),
            "preserve external CSV target\n",
        )

    def test_initialize_fsyncs_created_ledger_directories(self) -> None:
        root = Path(self.temporary.name) / "durable-genesis"
        fsynced_modes: list[int] = []
        real_fsync = os.fsync

        def record_fsync(descriptor: int) -> None:
            fsynced_modes.append(os.fstat(descriptor).st_mode)
            real_fsync(descriptor)

        with patch("ledger.os.fsync", side_effect=record_fsync):
            initialize_ledger(
                root / "orion" / "ledger" / "node_hours.jsonl",
                root / "orion" / "ledger" / "node_hours.csv",
                root / "orion" / "ledger" / "mirror_receipts.jsonl",
                root / "project_home" / "ledger" / "node_hours.jsonl",
                mirror_transport="filesystem_copy",
                notes="durable test genesis",
                control_plane_version="durable-test-control-plane",
            )
        self.assertGreaterEqual(sum(stat.S_ISDIR(mode) for mode in fsynced_modes), 1)

    def test_csv_replace_fsyncs_parent_directory(self) -> None:
        fsynced_modes: list[int] = []
        real_fsync = os.fsync

        def record_fsync(descriptor: int) -> None:
            fsynced_modes.append(os.fstat(descriptor).st_mode)
            real_fsync(descriptor)

        with patch("ledger.os.fsync", side_effect=record_fsync):
            write_csv(
                self.ledger,
                self.receipts,
                self.csv,
                mirror_jsonl=self.mirror,
                mirror_transport="filesystem_copy",
            )
        self.assertTrue(any(stat.S_ISDIR(mode) for mode in fsynced_modes))

    def test_csv_rejects_missing_genesis_anchors_before_rewrite(self) -> None:
        original = self.csv.read_bytes()
        for anchor in genesis_anchor_paths(self.ledger, self.mirror):
            anchor.chmod(0o600)
            anchor.unlink()
        with self.assertRaises(ValueError):
            write_csv(
                self.ledger,
                self.receipts,
                self.csv,
                mirror_jsonl=self.mirror,
                mirror_transport="filesystem_copy",
            )
        self.assertEqual(self.csv.read_bytes(), original)

    def test_csv_parent_fsync_failure_rolls_back_previous_projection(self) -> None:
        original = self.csv.read_bytes()
        real_fsync = os.fsync
        failed = False

        def fail_first_directory_fsync(descriptor: int) -> None:
            nonlocal failed
            if stat.S_ISDIR(os.fstat(descriptor).st_mode) and not failed:
                failed = True
                raise OSError("directory fsync failed")
            real_fsync(descriptor)

        with patch(
            "control_plane_common.os.fsync",
            side_effect=fail_first_directory_fsync,
        ):
            with self.assertRaises(OSError):
                write_csv(
                    self.ledger,
                    self.receipts,
                    self.csv,
                    mirror_jsonl=self.mirror,
                    mirror_transport="filesystem_copy",
                )
        self.assertTrue(failed)
        self.assertEqual(self.csv.read_bytes(), original)

    def test_initialize_accepts_trusted_mount_alias_above_existing_root(self) -> None:
        root = Path(self.temporary.name) / "trusted-alias-genesis"
        real_parent = root / "real-parent"
        real_mirror_root = real_parent / "project-home"
        real_mirror_root.mkdir(parents=True)
        alias_parent = root / "project-home-alias"
        alias_parent.symlink_to(real_parent, target_is_directory=True)
        mirror = alias_parent / "project-home" / "ledger" / "node_hours.jsonl"
        ledger = root / "orion" / "ledger" / "node_hours.jsonl"
        receipts = root / "orion" / "ledger" / "mirror_receipts.jsonl"
        initialize_ledger(
            ledger,
            root / "orion" / "ledger" / "node_hours.csv",
            receipts,
            mirror,
            mirror_transport="filesystem_copy",
            notes="trusted alias genesis",
            control_plane_version="trusted-alias-control-plane",
        )
        records = validate_mirrored_state(ledger, receipts, mirror)
        self.assertEqual(len(records), 1)
        receipt = json.loads(receipts.read_text(encoding="utf-8"))
        self.assertEqual(receipt["mirror_destination"], str(mirror))

    def test_truncated_ledger_files_cannot_reinitialize_after_genesis(self) -> None:
        for path in [self.ledger, self.csv, self.receipts, self.mirror]:
            path.write_bytes(b"")
        with self.assertRaises(ValueError):
            initialize_ledger(
                self.ledger,
                self.csv,
                self.receipts,
                self.mirror,
                mirror_transport="filesystem_copy",
                notes="must not recreate genesis",
                control_plane_version="replacement-control-plane",
            )
        with self.assertRaises(ValueError):
            validate_mirrored_state(self.ledger, self.receipts, self.mirror)

    def test_interrupted_fresh_genesis_primary_only_is_repairable(self) -> None:
        root = Path(self.temporary.name) / "interrupted-genesis-primary"
        ledger = root / "orion" / "ledger" / "node_hours.jsonl"
        csv_path = root / "orion" / "ledger" / "node_hours.csv"
        receipts = root / "orion" / "ledger" / "mirror_receipts.jsonl"
        mirror = root / "project_home" / "ledger" / "node_hours.jsonl"
        real_append = __import__("ledger")._append_jsonl

        def interrupt_mirror(path: Path, record: dict[str, object]) -> None:
            if path == mirror:
                raise RuntimeError("simulated mirror append interruption")
            real_append(path, record)

        with patch("ledger._append_jsonl", side_effect=interrupt_mirror):
            with self.assertRaises(RuntimeError):
                initialize_ledger(
                    ledger,
                    csv_path,
                    receipts,
                    mirror,
                    mirror_transport="filesystem_copy",
                    notes="recover exact interrupted genesis",
                    control_plane_version="recoverable-control-plane",
                )
        initialize_ledger(
            ledger,
            csv_path,
            receipts,
            mirror,
            mirror_transport="filesystem_copy",
            notes="recover exact interrupted genesis",
            control_plane_version="recoverable-control-plane",
        )
        self.assertEqual(len(validate_mirrored_state(ledger, receipts, mirror)), 1)

    def test_interrupted_fresh_genesis_receipt_only_is_repairable(self) -> None:
        root = Path(self.temporary.name) / "interrupted-genesis-receipt"
        ledger = root / "orion" / "ledger" / "node_hours.jsonl"
        csv_path = root / "orion" / "ledger" / "node_hours.csv"
        receipts = root / "orion" / "ledger" / "mirror_receipts.jsonl"
        mirror = root / "project_home" / "ledger" / "node_hours.jsonl"
        with patch(
            "ledger._write_genesis_anchors",
            side_effect=RuntimeError("simulated anchor publication interruption"),
        ):
            with self.assertRaises(RuntimeError):
                initialize_ledger(
                    ledger,
                    csv_path,
                    receipts,
                    mirror,
                    mirror_transport="filesystem_copy",
                    notes="recover exact interrupted genesis",
                    control_plane_version="recoverable-control-plane",
                )
        initialize_ledger(
            ledger,
            csv_path,
            receipts,
            mirror,
            mirror_transport="filesystem_copy",
            notes="recover exact interrupted genesis",
            control_plane_version="recoverable-control-plane",
        )
        self.assertEqual(len(validate_mirrored_state(ledger, receipts, mirror)), 1)

    def test_interrupted_fresh_genesis_mirror_only_is_repairable(self) -> None:
        root = Path(self.temporary.name) / "interrupted-genesis-mirror"
        ledger = root / "orion" / "ledger" / "node_hours.jsonl"
        csv_path = root / "orion" / "ledger" / "node_hours.csv"
        receipts = root / "orion" / "ledger" / "mirror_receipts.jsonl"
        mirror = root / "project_home" / "ledger" / "node_hours.jsonl"
        with patch(
            "ledger._append_mirror_receipt",
            side_effect=RuntimeError("simulated receipt append interruption"),
        ):
            with self.assertRaises(RuntimeError):
                initialize_ledger(
                    ledger,
                    csv_path,
                    receipts,
                    mirror,
                    mirror_transport="filesystem_copy",
                    notes="recover exact interrupted genesis",
                    control_plane_version="recoverable-control-plane",
                )
        initialize_ledger(
            ledger,
            csv_path,
            receipts,
            mirror,
            mirror_transport="filesystem_copy",
            notes="recover exact interrupted genesis",
            control_plane_version="recoverable-control-plane",
        )
        self.assertEqual(len(validate_mirrored_state(ledger, receipts, mirror)), 1)

    def test_interrupted_fresh_genesis_one_anchor_is_repairable(self) -> None:
        root = Path(self.temporary.name) / "interrupted-genesis-one-anchor"
        ledger = root / "orion" / "ledger" / "node_hours.jsonl"
        csv_path = root / "orion" / "ledger" / "node_hours.csv"
        receipts = root / "orion" / "ledger" / "mirror_receipts.jsonl"
        mirror = root / "project_home" / "ledger" / "node_hours.jsonl"
        real_write = __import__("ledger")._write_one_genesis_anchor

        def interrupt_second_anchor(path: Path, root: Path, anchor: dict[str, object]) -> None:
            if path == mirror.parent / "genesis_anchor.json":
                raise RuntimeError("simulated second anchor interruption")
            real_write(path, root, anchor)

        with patch("ledger._write_one_genesis_anchor", side_effect=interrupt_second_anchor):
            with self.assertRaises(RuntimeError):
                initialize_ledger(
                    ledger,
                    csv_path,
                    receipts,
                    mirror,
                    mirror_transport="filesystem_copy",
                    notes="recover exact interrupted genesis",
                    control_plane_version="recoverable-control-plane",
                )
        initialize_ledger(
            ledger,
            csv_path,
            receipts,
            mirror,
            mirror_transport="filesystem_copy",
            notes="recover exact interrupted genesis",
            control_plane_version="recoverable-control-plane",
        )
        self.assertEqual(len(validate_mirrored_state(ledger, receipts, mirror)), 1)

    def test_interrupted_fresh_genesis_missing_csv_is_repairable(self) -> None:
        root = Path(self.temporary.name) / "interrupted-genesis-csv"
        ledger = root / "orion" / "ledger" / "node_hours.jsonl"
        csv_path = root / "orion" / "ledger" / "node_hours.csv"
        receipts = root / "orion" / "ledger" / "mirror_receipts.jsonl"
        mirror = root / "project_home" / "ledger" / "node_hours.jsonl"
        with patch("ledger.write_csv", side_effect=RuntimeError("simulated CSV interruption")):
            with self.assertRaises(RuntimeError):
                initialize_ledger(
                    ledger,
                    csv_path,
                    receipts,
                    mirror,
                    mirror_transport="filesystem_copy",
                    notes="recover exact interrupted genesis",
                    control_plane_version="recoverable-control-plane",
                )
        initialize_ledger(
            ledger,
            csv_path,
            receipts,
            mirror,
            mirror_transport="filesystem_copy",
            notes="recover exact interrupted genesis",
            control_plane_version="recoverable-control-plane",
        )
        self.assertEqual(len(validate_mirrored_state(ledger, receipts, mirror)), 1)

    def test_interrupted_fresh_genesis_rejects_identity_drift(self) -> None:
        root = Path(self.temporary.name) / "interrupted-genesis-identity"
        ledger = root / "orion" / "ledger" / "node_hours.jsonl"
        csv_path = root / "orion" / "ledger" / "node_hours.csv"
        receipts = root / "orion" / "ledger" / "mirror_receipts.jsonl"
        mirror = root / "project_home" / "ledger" / "node_hours.jsonl"
        real_append = __import__("ledger")._append_jsonl

        def interrupt_mirror(path: Path, record: dict[str, object]) -> None:
            if path == mirror:
                raise RuntimeError("simulated mirror append interruption")
            real_append(path, record)

        with patch("ledger._append_jsonl", side_effect=interrupt_mirror):
            with self.assertRaises(RuntimeError):
                initialize_ledger(
                    ledger,
                    csv_path,
                    receipts,
                    mirror,
                    mirror_transport="filesystem_copy",
                    notes="recover exact interrupted genesis",
                    control_plane_version="recoverable-control-plane",
                )
        for notes, version in [
            ("changed notes", "recoverable-control-plane"),
            ("recover exact interrupted genesis", "changed-control-plane"),
        ]:
            with self.assertRaises(ValueError):
                initialize_ledger(
                    ledger,
                    csv_path,
                    receipts,
                    mirror,
                    mirror_transport="filesystem_copy",
                    notes=notes,
                    control_plane_version=version,
                )

    def test_interrupted_fresh_genesis_rejects_corrupt_receipt(self) -> None:
        root = Path(self.temporary.name) / "interrupted-genesis-corrupt-receipt"
        ledger = root / "orion" / "ledger" / "node_hours.jsonl"
        csv_path = root / "orion" / "ledger" / "node_hours.csv"
        receipts = root / "orion" / "ledger" / "mirror_receipts.jsonl"
        mirror = root / "project_home" / "ledger" / "node_hours.jsonl"
        with patch(
            "ledger._write_genesis_anchors",
            side_effect=RuntimeError("simulated anchor publication interruption"),
        ):
            with self.assertRaises(RuntimeError):
                initialize_ledger(
                    ledger,
                    csv_path,
                    receipts,
                    mirror,
                    mirror_transport="filesystem_copy",
                    notes="recover exact interrupted genesis",
                    control_plane_version="recoverable-control-plane",
                )
        receipt = json.loads(receipts.read_text(encoding="utf-8"))
        receipt["mirrored_event_sha256"] = "0" * 64
        receipts.write_text(json.dumps(receipt) + "\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            initialize_ledger(
                ledger,
                csv_path,
                receipts,
                mirror,
                mirror_transport="filesystem_copy",
                notes="recover exact interrupted genesis",
                control_plane_version="recoverable-control-plane",
            )

    def test_audited_pre_anchor_ledger_migration_is_one_shot(self) -> None:
        local_anchor, mirror_anchor = genesis_anchor_paths(self.ledger, self.mirror)
        event = json.loads(self.ledger.read_text(encoding="utf-8"))
        receipt = json.loads(self.receipts.read_text(encoding="utf-8"))
        local_anchor.chmod(0o600)
        mirror_anchor.chmod(0o600)
        local_anchor.unlink()
        mirror_anchor.unlink()
        anchor = migrate_existing_genesis_anchors(
            self.ledger,
            self.receipts,
            self.mirror,
            expected_event_sha256=event["event_sha256"],
            expected_mirror_ack_sha256=receipt["mirror_ack_sha256"],
        )
        self.assertEqual(anchor["event_sha256"], event["event_sha256"])
        self.assertEqual(len(validate_mirrored_state(
            self.ledger, self.receipts, self.mirror
        )), 1)
        with self.assertRaises(ValueError):
            migrate_existing_genesis_anchors(
                self.ledger,
                self.receipts,
                self.mirror,
                expected_event_sha256="0" * 64,
                expected_mirror_ack_sha256=receipt["mirror_ack_sha256"],
            )

    def test_audited_pre_anchor_migration_repairs_one_missing_anchor(self) -> None:
        local_anchor, mirror_anchor = genesis_anchor_paths(self.ledger, self.mirror)
        event = json.loads(self.ledger.read_text(encoding="utf-8"))
        receipt = json.loads(self.receipts.read_text(encoding="utf-8"))
        mirror_anchor.chmod(0o600)
        mirror_anchor.unlink()
        with self.assertRaises(ValueError):
            validate_mirrored_state(self.ledger, self.receipts, self.mirror)
        anchor = migrate_existing_genesis_anchors(
            self.ledger,
            self.receipts,
            self.mirror,
            expected_event_sha256=event["event_sha256"],
            expected_mirror_ack_sha256=receipt["mirror_ack_sha256"],
        )
        self.assertTrue(local_anchor.is_file())
        self.assertTrue(mirror_anchor.is_file())
        self.assertEqual(anchor["event_sha256"], event["event_sha256"])


if __name__ == "__main__":
    unittest.main()
