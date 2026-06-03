#!/usr/bin/env python3
"""Unit tests for the mirrored Frontier PIC ledger."""

from __future__ import annotations

import csv
import fcntl
import hashlib
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
from ledger import slurm_walltime_seconds
from ledger import transition_payload, validate_primary_chain, validate_receipts, write_csv
from ledger import validated_read_only_mirrored_state_snapshot
from operator_attestation import CAPTURED_SNAPSHOT_FILENAMES, OPERATOR_STATEMENTS


class LedgerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        self.root = root
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

    def _reservation_event(
        self,
        *,
        reservation_id: str = "reservation-1",
        include_operator_attestations: bool = True,
        authorization_id: str = "test-authorization",
    ) -> dict[str, object]:
        event = {
            "event_type": "reservation",
            "reservation_id": reservation_id,
            "submission_id": "submission-1",
            "control_plane_version": "a" * 64,
            "active_policy_sha256": "b" * 64,
            "active_promotion_sha256": "c" * 64,
            "git_commit": "d" * 40,
            "campaign": "test-campaign",
            "test_id": "test-id",
            "manifest_path": str(
                self.root
                / "orion"
                / "manifests"
                / "test-campaign"
                / "submission-1"
                / "manifest.json"
            ),
            "submission_scope": "registered_science",
            "registered_science_authorization_id": authorization_id,
            "clean_candidate_manifest_sha256": "e" * 64,
            "partition": "batch",
            "qos": "normal",
            "qos_selection_reason": "test",
            "queue_snapshot_sha256": "f" * 64,
            "manifest_sha256": "1" * 64,
            "job_script_sha256": "2" * 64,
            "executable_sha256": "3" * 64,
            "site_policy_checked_utc": "2026-05-31T00:00:00Z",
            "state": "reserved",
            "reconciled": False,
            "requested_nodes": 1,
            "requested_walltime": "00:01:00",
            "reserved_node_hours": 1.0 / 60.0,
            "artifact_dir": str(self.root / "artifacts" / "submission-1"),
        }
        if include_operator_attestations:
            event.update(self._operator_attestation_quartet(authorization_id))
        return event

    def _operator_attestation_quartet(
        self, authorization_id: str = "test-authorization"
    ) -> dict[str, str]:
        archive = self.root / "orion" / "operator_attestations"
        bindings = {}
        for index, phase in enumerate(["pre_manifest", "pre_submit_wrapper"]):
            timestamp = f"20260601T00000{index}Z"
            root = archive / f"{timestamp}-{authorization_id}-{phase}"
            if not root.exists():
                root.mkdir(parents=True)
                manual_paths = [
                    self.root / "orion" / "ledger" / "pending_manual_accounting.json",
                    self.root
                    / "project_home"
                    / "ledger"
                    / "pending_manual_accounting.json",
                ]
                ledger_paths = [self.ledger, self.receipts, self.mirror]
                manual = "".join(f"{path} absent\n" for path in manual_paths).encode()
                counts = "".join(f"1 {path}\n" for path in ledger_paths).encode()
                state = {
                    "validation": "coherent",
                    "validator": "existing_read_only_mirrored_state_snapshot",
                    "ledger_record_count": 1,
                    "ledger_tail_event_sha256": None,
                    "active_reservation_count": 0,
                    "active_reservation_ids": [],
                }
                payloads = {
                    "same_account_process_snapshot.txt": b"test process snapshot\n",
                    "queue_snapshot.txt": b"",
                    "pending_submission_marker.txt": b"absent\n",
                    "pending_manual_accounting_marker.txt": manual,
                    "mirrored_ledger_line_counts.txt": counts,
                    "validated_mirrored_ledger_state.json": (
                        json.dumps(state, indent=2, sort_keys=True) + "\n"
                    ).encode(),
                }
                payloads.update(
                    {
                        f"capture_{name}": payload
                        for name, payload in list(payloads.items())
                        if name != "same_account_process_snapshot.txt"
                    }
                )
                payloads["capture_same_account_process_snapshot.txt"] = (
                    b"test captured process snapshot\n"
                )
                for filename, payload in payloads.items():
                    path = root / filename
                    path.write_bytes(payload)
                    path.chmod(0o400)

                def file_record(filename: str) -> dict[str, object]:
                    payload = (root / filename).read_bytes()
                    return {
                        "path": filename,
                        "sha256": hashlib.sha256(payload).hexdigest(),
                    }

                pending = file_record("pending_submission_marker.txt")
                pending["value"] = "absent"
                pending_manual = file_record("pending_manual_accounting_marker.txt")
                pending_manual["value"] = "absent"
                pending_manual["values"] = {
                    str(path): "absent" for path in manual_paths
                }
                line_counts = file_record("mirrored_ledger_line_counts.txt")
                line_counts["counts"] = {str(path): 1 for path in ledger_paths}
                validated_state = file_record("validated_mirrored_ledger_state.json")
                validated_state["state"] = state
                attestation = {
                    "schema_version": 1,
                    "record_type": (
                        "q027_frontier_registered_science_same_account_isolation_attestation"
                    ),
                    "recorded_utc": f"2026-06-01T00:00:0{index}Z",
                    "sealed_utc": f"2026-06-01T00:00:0{index}Z",
                    "registered_science_authorization_id": authorization_id,
                    "control_plane_version": "a" * 64,
                    "phase": phase,
                    "same_account_process_snapshot": file_record(
                        "same_account_process_snapshot.txt"
                    ),
                    "queue_snapshot": file_record("queue_snapshot.txt"),
                    "pending_submission_marker": pending,
                    "pending_manual_accounting_marker": pending_manual,
                    "mirrored_ledger_line_counts": line_counts,
                    "validated_mirrored_ledger_state": validated_state,
                    "captured_snapshots": {
                        filename: file_record(filename)
                        for filename in sorted(CAPTURED_SNAPSHOT_FILENAMES)
                    },
                    "operator_statement": OPERATOR_STATEMENTS[phase],
                }
                attestation_path = root / "attestation.json"
                attestation_path.write_text(
                    json.dumps(attestation, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
                attestation_path.chmod(0o400)
                root.chmod(0o500)
            attestation_path = root / "attestation.json"
            bindings[f"{phase}_attestation_path"] = str(attestation_path)
            bindings[f"{phase}_attestation_sha256"] = hashlib.sha256(
                attestation_path.read_bytes()
            ).hexdigest()
        return bindings

    def _recovery_handoff(
        self,
        reservation: dict[str, object],
        *,
        mode: str = "fresh_scheduler_binding",
    ) -> tuple[Path, str]:
        name = "35181896-6465-4d01-b80f-69841dcf947c.json"
        handoff: dict[str, object] = {
            "schema_version": 1,
            "status": "authorized_terminal_scheduler_job_id_received_recovery",
            "handoff_id": name.removesuffix(".json"),
            "reservation_id": reservation["reservation_id"],
            "submission_id": reservation["submission_id"],
            "job_id": "1234",
            "manifest_path": reservation["manifest_path"],
            "manifest_sha256": reservation["manifest_sha256"],
            "prior_control_plane_version": reservation["control_plane_version"],
            "recovery_control_plane_version": "b" * 64,
            "prior_active_policy_sha256": reservation["active_policy_sha256"],
            "prior_active_promotion_sha256": reservation["active_promotion_sha256"],
            "pending_marker_sha256": "e" * 64,
        }
        if mode == "purged_scontrol_cancelled_zero_execution":
            handoff["scheduler_binding_recovery"] = {
                "mode": mode,
                "job_id": "1234",
                "job_name": "run_installed_control_plane_job.sh",
                "state": "CANCELLED by 123",
                "elapsed_raw": 0,
                "allocated_nodes": 0,
                "comment": "",
                "account": "ast207",
                "submit": "2026-05-30T00:00:00",
                "start": "None",
                "end": "2026-05-30T00:00:00",
                "exit_code": "0:0",
            }
            handoff["reservation_job_binding_attestation"] = {
                "mode": (
                    "reviewed_operator_attestation_for_unprovable_purged_"
                    "reservation_job_binding"
                ),
                "reservation_id": reservation["reservation_id"],
                "job_id": "1234",
            }
        payload = (json.dumps(handoff, sort_keys=True) + "\n").encode("utf-8")
        paths = [
            self.root / "orion" / "policy" / "recovery_handoffs" / name,
            self.root / "project_home" / "policy" / "recovery_handoffs" / name,
        ]
        for path in paths:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
            path.chmod(0o444)
        return paths[0], hashlib.sha256(payload).hexdigest()

    def _recovery_reservation_event(
        self,
        *,
        reservation_id: str = "reservation-1",
        submission_id: str = "submission",
    ) -> dict[str, object]:
        event = self._reservation_event(reservation_id=reservation_id)
        event.update(
            {
                "manifest_path": str(
                    self.root
                    / "orion"
                    / "manifests"
                    / "campaign"
                    / submission_id
                    / "manifest.json"
                ),
                "submission_id": submission_id,
                "manifest_sha256": "f" * 64,
                "control_plane_version": "a" * 64,
                "active_policy_sha256": "c" * 64,
                "active_promotion_sha256": "d" * 64,
            }
        )
        return event

    def _rewrite_recovery_handoff(
        self, handoff_path: Path, mutate: object
    ) -> str:
        handoff = json.loads(handoff_path.read_text(encoding="utf-8"))
        mutate(handoff)
        payload = (json.dumps(handoff, sort_keys=True) + "\n").encode("utf-8")
        for root in [self.root / "orion", self.root / "project_home"]:
            path = root / "policy" / "recovery_handoffs" / handoff_path.name
            path.chmod(0o600)
            path.write_bytes(payload)
            path.chmod(0o444)
        return hashlib.sha256(payload).hexdigest()

    def _drop_last_line(self, path: Path) -> None:
        lines = path.read_text(encoding="utf-8").splitlines()
        path.write_text("".join(line + "\n" for line in lines[:-1]), encoding="utf-8")

    def test_reserve_attach_reconcile_and_csv_projection(self) -> None:
        common = {
            **transition_payload(
                self._reservation_event(authorization_id="f1-clean-gyro-v1")
            ),
            "reservation_id": "reservation-1",
            "submission_id": "submission-1",
            "submission_scope": "registered_science",
            "registered_science_authorization_id": "f1-clean-gyro-v1",
            "clean_candidate_manifest_sha256": "a" * 64,
            "requested_nodes": 2,
            "requested_walltime": "00:30:00",
            "reserved_node_hours": 1.0,
            "reconciled": False,
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

    def test_manual_direct_scheduler_reconciliation_counts_toward_consumed_budget(
        self,
    ) -> None:
        event = self._manual_accounting_event()
        event["accounting_scope"] = "manual_direct_scheduler_accounting_only"
        event["notes"] = (
            "Reviewed direct-scheduler accounting only; "
            "ineligible for scientific evidence."
        )
        self.append(event)
        totals = accounting(validate_primary_chain(self.ledger))
        self.assertAlmostEqual(
            totals["cumulative_consumed_node_hours"], 1.0 / 720.0
        )
        self.assertEqual(totals["currently_reserved_node_hours"], 0.0)

    def test_manual_direct_srun_reconciliation_rejects_under_bound_event(
        self,
    ) -> None:
        with self.assertRaisesRegex(ValueError, "schema is unsupported"):
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

    def test_manual_reconciliation_rejects_nonstring_accounting_scope(self) -> None:
        event = self._manual_accounting_event()
        event["accounting_scope"] = []
        with self.assertRaisesRegex(
            ValueError, "Manual-accounting ledger event semantics are invalid"
        ):
            self.append(event)

    def test_manual_reconciliation_rejects_unknown_accounting_scope(self) -> None:
        event = self._manual_accounting_event()
        event["accounting_scope"] = "unknown_manual_scope"
        with self.assertRaisesRegex(
            ValueError, "Manual-accounting ledger event semantics are invalid"
        ):
            self.append(event)

    def test_registered_reconciliation_rejects_inconsistent_usage(self) -> None:
        common = {
            **transition_payload(
                self._reservation_event(authorization_id="f1-clean-gyro-v1")
            ),
            "reservation_id": "reservation-1",
            "submission_id": "submission-1",
            "submission_scope": "registered_science",
            "registered_science_authorization_id": "f1-clean-gyro-v1",
            "clean_candidate_manifest_sha256": "a" * 64,
            "requested_nodes": 2,
            "requested_walltime": "00:30:00",
            "reserved_node_hours": 1.0,
            "reconciled": False,
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
        event = {
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
        for field, value, message in [
            ("scheduler_reported_allocated_nodes", "2", "usage is invalid"),
            ("billed_nodes", 999, "usage is invalid"),
            ("billed_nodes", 2.0, "usage is invalid"),
            ("elapsed_seconds", -1, "usage is invalid"),
            ("consumed_node_hours", 0.001, "usage differs"),
            ("cumulative_consumed_node_hours", 999.0, "cumulative usage differs"),
        ]:
            with self.subTest(field=field):
                mutated = dict(event)
                mutated[field] = value
                with self.assertRaisesRegex(ValueError, message):
                    self.append(mutated)

    def test_registered_reconciliation_rejects_duplicate_terminal_event(self) -> None:
        common = {
            **transition_payload(self._reservation_event()),
            "reservation_id": "reservation-1",
            "submission_id": "submission-1",
            "requested_nodes": 2,
            "requested_walltime": "00:30:00",
            "reserved_node_hours": 1.0,
            "reconciled": False,
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
        event = {
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
        self.append(event)
        event["cumulative_consumed_node_hours"] = 2.0 / 3.0
        with self.assertRaisesRegex(ValueError, "transition is invalid"):
            self.append(event)

    def test_registered_reconciliation_rejects_nonterminal_state(self) -> None:
        reservation = self.append(self._reservation_event())
        attachment = transition_payload(reservation)
        attachment.update(
            {"event_type": "job_id_attached", "job_id": "1234", "state": "submitted"}
        )
        attachment = self.append(attachment)
        reconciliation = transition_payload(attachment)
        reconciliation.update(
            {
                "event_type": "reconciliation",
                "state": "RUNNING",
                "reconciled": True,
                "scheduler_reported_allocated_nodes": 1,
                "billed_nodes": 1,
                "elapsed_seconds": 60,
                "consumed_node_hours": 1.0 / 60.0,
                "cumulative_consumed_node_hours": 1.0 / 60.0,
            }
        )
        with self.assertRaisesRegex(ValueError, "state is not terminal"):
            self.append(reconciliation)

    def test_registered_attachment_rejects_scheduler_job_id_reuse(self) -> None:
        for reservation_id in ["reservation-1", "reservation-2"]:
            reservation = self.append(
                self._reservation_event(reservation_id=reservation_id)
            )
            attachment = transition_payload(reservation)
            attachment.update(
                {
                    "event_type": "job_id_attached",
                    "job_id": "1234",
                    "state": "submitted",
                }
            )
            if reservation_id == "reservation-1":
                self.append(attachment)
            else:
                with self.assertRaisesRegex(ValueError, "already has a ledger owner"):
                    self.append(attachment)

    def test_registered_lifecycle_rejects_payload_rewrites_and_unknown_events(
        self,
    ) -> None:
        reservation = self.append(self._reservation_event())
        attachment = transition_payload(reservation)
        attachment.update(
            {"event_type": "job_id_attached", "job_id": "1234", "state": "submitted"}
        )
        rewritten_attachment = dict(attachment)
        rewritten_attachment["reserved_node_hours"] = 0.0
        with self.assertRaisesRegex(ValueError, "rewrites immutable fields"):
            self.append(rewritten_attachment)
        attachment = self.append(attachment)
        reconciliation = transition_payload(attachment)
        reconciliation.update(
            {
                "event_type": "reconciliation",
                "state": "COMPLETED",
                "reconciled": True,
                "scheduler_reported_allocated_nodes": 1,
                "billed_nodes": 1,
                "elapsed_seconds": 3600,
                "consumed_node_hours": 1.0,
                "cumulative_consumed_node_hours": 1.0,
            }
        )
        rewritten_reconciliation = dict(reconciliation)
        rewritten_reconciliation["requested_nodes"] = 10
        rewritten_reconciliation["billed_nodes"] = 10
        rewritten_reconciliation["consumed_node_hours"] = 10.0
        rewritten_reconciliation["cumulative_consumed_node_hours"] = 10.0
        with self.assertRaisesRegex(ValueError, "rewrites immutable fields"):
            self.append(rewritten_reconciliation)
        self.append(reconciliation)
        with self.assertRaisesRegex(ValueError, "cancellation transition is invalid"):
            self.append(
                {
                    **transition_payload(reservation),
                    "event_type": "reservation_cancelled",
                    "state": "cancelled",
                }
            )
        with self.assertRaisesRegex(ValueError, "missing or unsupported"):
            self.append(
                {
                    "event_type": "opaque",
                    "reservation_id": "reservation-1",
                    "state": "cancelled",
                }
            )

    def test_registered_attachment_rejects_immutable_numeric_aliases(self) -> None:
        reservation = self.append(self._reservation_event())
        attachment = transition_payload(reservation)
        attachment.update(
            {"event_type": "job_id_attached", "job_id": "1234", "state": "submitted"}
        )
        for value in [1.0, True]:
            with self.subTest(value=value):
                aliased = dict(attachment)
                aliased["requested_nodes"] = value
                with self.assertRaisesRegex(ValueError, "rewrites immutable fields"):
                    self.append(aliased)

    def test_registered_reservation_rejects_false_reserved_usage(self) -> None:
        for field, value in [
            ("requested_walltime", "00:00:00"),
            ("requested_walltime", "invalid"),
            ("reserved_node_hours", 0.0),
        ]:
            with self.subTest(field=field, value=value):
                event = self._reservation_event()
                event[field] = value
                with self.assertRaisesRegex(ValueError, "reservation"):
                    self.append(event)

    def test_registered_reservation_rejects_partial_operator_attestation_quartet(
        self,
    ) -> None:
        quartet = self._operator_attestation_quartet()
        for field, value in quartet.items():
            with self.subTest(field=field):
                event = self._reservation_event(include_operator_attestations=False)
                event[field] = value
                with self.assertRaisesRegex(ValueError, "quartet is incomplete"):
                    self.append(event)

    def test_registered_reservation_accepts_complete_operator_attestation_quartet(
        self,
    ) -> None:
        event = self._reservation_event()
        event.update(self._operator_attestation_quartet())
        reservation = self.append(event)
        self.assertEqual(
            reservation["pre_submit_wrapper_attestation_sha256"],
            self._operator_attestation_quartet()[
                "pre_submit_wrapper_attestation_sha256"
            ],
        )

    def test_registered_reservation_rejects_missing_operator_attestation_quartet(
        self,
    ) -> None:
        with self.assertRaisesRegex(
            ValueError, "requires operator-attestation provenance"
        ):
            self.append(
                self._reservation_event(include_operator_attestations=False)
            )

    def test_historical_registered_reservation_accepts_missing_operator_attestation_quartet(
        self,
    ) -> None:
        event = self._reservation_event(include_operator_attestations=False)
        event["control_plane_version"] = (
            "6002c80e305d6cfd322675b3e27e17e169b1718d8edd722330464cccb0c4fd86"
        )
        self.append(event)

    def test_registered_reservation_rejects_off_root_operator_attestation_quartet(
        self,
    ) -> None:
        event = self._reservation_event()
        for phase in ["pre_manifest", "pre_submit_wrapper"]:
            event[f"{phase}_attestation_path"] = str(
                self.root
                / "forged"
                / "operator_attestations"
                / f"forged-{phase}"
                / "attestation.json"
            )
        with self.assertRaisesRegex(ValueError, "fixed archive layout"):
            self.append(event)

    def test_slurm_walltime_parser_rounds_and_supports_day_forms(self) -> None:
        self.assertEqual(slurm_walltime_seconds("10"), 600)
        self.assertEqual(slurm_walltime_seconds("60"), 3600)
        self.assertEqual(slurm_walltime_seconds("120"), 7200)
        self.assertEqual(slurm_walltime_seconds("10:01"), 660)
        self.assertEqual(slurm_walltime_seconds("60:00"), 3600)
        self.assertEqual(slurm_walltime_seconds("00:10:01"), 660)
        self.assertEqual(slurm_walltime_seconds("1-02"), 93600)
        self.assertEqual(slurm_walltime_seconds("1-02:03"), 93780)
        self.assertEqual(slurm_walltime_seconds("0-01:02"), 3720)

    def test_registered_reservation_rejects_unrounded_reserved_usage(self) -> None:
        event = self._reservation_event()
        event["requested_walltime"] = "00:10:01"
        event["reserved_node_hours"] = 601.0 / 3600.0
        with self.assertRaisesRegex(ValueError, "reservation"):
            self.append(event)

    def test_scheduler_job_ids_are_canonical_and_globally_unique(self) -> None:
        reservation = self.append(self._reservation_event())
        for job_id in ["", " ", "0123", "0"]:
            with self.subTest(job_id=job_id):
                attachment = transition_payload(reservation)
                attachment.update(
                    {
                        "event_type": "job_id_attached",
                        "job_id": job_id,
                        "state": "submitted",
                    }
                )
                with self.assertRaisesRegex(ValueError, "not canonical"):
                    self.append(attachment)

        manual = self._manual_accounting_event()
        manual["job_id"] = "0123"
        with self.assertRaisesRegex(ValueError, "not canonical"):
            self.append(manual)

        manual = self.append(self._manual_accounting_event())
        duplicate = dict(manual)
        duplicate.pop("sequence_number")
        duplicate.pop("previous_event_sha256")
        duplicate.pop("event_sha256")
        duplicate.pop("timestamp")
        with self.assertRaisesRegex(ValueError, "already has a ledger owner"):
            self.append(duplicate)

        reservation = self.append(
            self._reservation_event(reservation_id="reservation-2")
        )
        attachment = transition_payload(reservation)
        attachment.update(
            {
                "event_type": "job_id_attached",
                "job_id": str(manual["job_id"]),
                "state": "submitted",
            }
        )
        with self.assertRaisesRegex(ValueError, "already has a ledger owner"):
            self.append(attachment)

        reservation = self.append(
            self._reservation_event(reservation_id="reservation-3")
        )
        attachment = transition_payload(reservation)
        attachment.update(
            {
                "event_type": "job_id_attached",
                "job_id": "9999",
                "state": "submitted",
            }
        )
        self.append(attachment)
        manual = self._manual_accounting_event()
        manual["job_id"] = "9999"
        with self.assertRaisesRegex(ValueError, "already has a ledger owner"):
            self.append(manual)

    def test_registered_recovery_provenance_is_valid_and_preserved(self) -> None:
        reservation_event = self._reservation_event()
        reservation_event["manifest_path"] = str(
            self.root / "orion" / "manifests" / "campaign" / "submission" / "manifest.json"
        )
        reservation_event["submission_id"] = "submission"
        reservation_event["manifest_sha256"] = "f" * 64
        reservation_event["control_plane_version"] = "a" * 64
        reservation_event["active_policy_sha256"] = "c" * 64
        reservation_event["active_promotion_sha256"] = "d" * 64
        reservation = self.append(reservation_event)
        valid_handoff_path, valid_handoff_sha256 = self._recovery_handoff(reservation)
        for fields in [
            {"terminal_recovery_handoff_path": "/tmp/handoff"},
            {
                "terminal_recovery_handoff_path": "relative",
                "terminal_recovery_handoff_sha256": "a" * 64,
                "terminal_recovery_mode": "fresh_scheduler_binding",
            },
            {
                "terminal_recovery_handoff_path": str(valid_handoff_path),
                "terminal_recovery_handoff_sha256": "not-a-hash",
                "terminal_recovery_mode": "fresh_scheduler_binding",
            },
            {
                "terminal_recovery_handoff_path": str(valid_handoff_path),
                "terminal_recovery_handoff_sha256": "a" * 64,
                "terminal_recovery_mode": "invented",
            },
        ]:
            with self.subTest(fields=fields):
                attachment = transition_payload(reservation)
                attachment.update(
                    {
                        "event_type": "job_id_attached",
                        "job_id": "1234",
                        "state": "submitted",
                        "attached_by_control_plane_version": "b" * 64,
                        **fields,
                    }
                )
                with self.assertRaisesRegex(
                    ValueError, "schema is unsupported|provenance"
                ):
                    self.append(attachment)

        attachment = transition_payload(reservation)
        attachment.update(
            {
                "event_type": "job_id_attached",
                "job_id": "1234",
                "state": "submitted",
                "attached_by_control_plane_version": "b" * 64,
                "terminal_recovery_handoff_path": str(valid_handoff_path),
                "terminal_recovery_handoff_sha256": valid_handoff_sha256,
                "terminal_recovery_mode": "fresh_scheduler_binding",
            }
        )
        attachment = self.append(attachment)
        reconciliation = transition_payload(attachment)
        for field in [
            "terminal_recovery_handoff_path",
            "terminal_recovery_handoff_sha256",
            "terminal_recovery_mode",
        ]:
            reconciliation.pop(field)
        reconciliation.update(
            {
                "event_type": "reconciliation",
                "state": "COMPLETED",
                "reconciled": True,
                "reconciled_by_control_plane_version": "b" * 64,
                "scheduler_reported_allocated_nodes": 1,
                "billed_nodes": 1,
                "elapsed_seconds": 60,
                "consumed_node_hours": 1.0 / 60.0,
                "cumulative_consumed_node_hours": 1.0 / 60.0,
            }
        )
        with self.assertRaisesRegex(ValueError, "rewrites immutable fields"):
            self.append(reconciliation)

    def test_registered_recovery_provenance_requires_mirrored_handoff(self) -> None:
        reservation_event = self._reservation_event()
        reservation_event["manifest_path"] = str(
            self.root / "orion" / "manifests" / "campaign" / "submission" / "manifest.json"
        )
        reservation = self.append(reservation_event)
        attachment = transition_payload(reservation)
        attachment.update(
            {
                "event_type": "job_id_attached",
                "job_id": "1234",
                "state": "submitted",
                "attached_by_control_plane_version": "b" * 64,
                "terminal_recovery_handoff_path": str(
                    self.root
                    / "orion"
                    / "policy"
                    / "recovery_handoffs"
                    / "35181896-6465-4d01-b80f-69841dcf947c.json"
                ),
                "terminal_recovery_handoff_sha256": "a" * 64,
                "terminal_recovery_mode": "fresh_scheduler_binding",
            }
        )
        with self.assertRaisesRegex(ValueError, "mirror is missing"):
            self.append(attachment)

    def test_registered_reconciliation_rejects_late_recovery_provenance(self) -> None:
        reservation = self.append(self._reservation_event())
        attachment = transition_payload(reservation)
        attachment.update(
            {"event_type": "job_id_attached", "job_id": "1234", "state": "submitted"}
        )
        attachment = self.append(attachment)
        reconciliation = transition_payload(attachment)
        reconciliation.update(
            {
                "event_type": "reconciliation",
                "state": "COMPLETED",
                "reconciled": True,
                "scheduler_reported_allocated_nodes": 1,
                "billed_nodes": 1,
                "elapsed_seconds": 60,
                "consumed_node_hours": 1.0 / 60.0,
                "cumulative_consumed_node_hours": 1.0 / 60.0,
                "terminal_recovery_handoff_path": (
                    "/tmp/policy/recovery_handoffs/"
                    "35181896-6465-4d01-b80f-69841dcf947c.json"
                ),
                "terminal_recovery_handoff_sha256": "a" * 64,
                "terminal_recovery_mode": "fresh_scheduler_binding",
            }
        )
        with self.assertRaisesRegex(ValueError, "adds unexpected fields"):
            self.append(reconciliation)

    def test_registered_reconciliation_rejects_false_purged_zero_execution(self) -> None:
        reservation_event = self._reservation_event()
        reservation_event["manifest_path"] = str(
            self.root / "orion" / "manifests" / "campaign" / "submission" / "manifest.json"
        )
        reservation_event["submission_id"] = "submission"
        reservation_event["manifest_sha256"] = "f" * 64
        reservation_event["control_plane_version"] = "a" * 64
        reservation_event["active_policy_sha256"] = "c" * 64
        reservation_event["active_promotion_sha256"] = "d" * 64
        reservation = self.append(reservation_event)
        handoff_path, handoff_sha256 = self._recovery_handoff(
            reservation, mode="purged_scontrol_cancelled_zero_execution"
        )
        attachment = transition_payload(reservation)
        attachment.update(
            {
                "event_type": "job_id_attached",
                "job_id": "1234",
                "state": "submitted",
                "attached_by_control_plane_version": "b" * 64,
                "terminal_recovery_handoff_path": str(handoff_path),
                "terminal_recovery_handoff_sha256": handoff_sha256,
                "terminal_recovery_mode": "purged_scontrol_cancelled_zero_execution",
            }
        )
        attachment = self.append(attachment)
        reconciliation = transition_payload(attachment)
        reconciliation.update(
            {
                "event_type": "reconciliation",
                "state": "COMPLETED",
                "reconciled": True,
                "reconciled_by_control_plane_version": "b" * 64,
                "scheduler_reported_allocated_nodes": 1,
                "billed_nodes": 1,
                "elapsed_seconds": 3600,
                "consumed_node_hours": 1.0,
                "cumulative_consumed_node_hours": 1.0,
            }
        )
        with self.assertRaisesRegex(ValueError, "zero-execution"):
            self.append(reconciliation)

    def test_registered_reconciliation_rejects_purged_zero_numeric_aliases(
        self,
    ) -> None:
        reservation = self.append(self._recovery_reservation_event())
        handoff_path, handoff_sha256 = self._recovery_handoff(
            reservation, mode="purged_scontrol_cancelled_zero_execution"
        )
        attachment = transition_payload(reservation)
        attachment.update(
            {
                "event_type": "job_id_attached",
                "job_id": "1234",
                "state": "submitted",
                "attached_by_control_plane_version": "b" * 64,
                "terminal_recovery_handoff_path": str(handoff_path),
                "terminal_recovery_handoff_sha256": handoff_sha256,
                "terminal_recovery_mode": "purged_scontrol_cancelled_zero_execution",
            }
        )
        attachment = self.append(attachment)
        reconciliation = transition_payload(attachment)
        reconciliation.update(
            {
                "event_type": "reconciliation",
                "state": "CANCELLED",
                "reconciled": True,
                "reconciled_by_control_plane_version": "b" * 64,
                "scheduler_reported_allocated_nodes": 0,
                "billed_nodes": 1,
                "elapsed_seconds": 0,
                "consumed_node_hours": 0.0,
                "cumulative_consumed_node_hours": 0.0,
            }
        )
        for field, value in [
            ("scheduler_reported_allocated_nodes", False),
            ("scheduler_reported_allocated_nodes", 0.0),
            ("elapsed_seconds", False),
            ("elapsed_seconds", 0.0),
            ("consumed_node_hours", False),
        ]:
            with self.subTest(field=field, value=value):
                aliased = dict(reconciliation)
                aliased[field] = value
                with self.assertRaisesRegex(ValueError, "zero-execution"):
                    self.append(aliased)

    def test_registered_recovery_handoff_is_bound_for_every_record(self) -> None:
        reservation = self.append(self._recovery_reservation_event())
        handoff_path, handoff_sha256 = self._recovery_handoff(reservation)
        attachment = transition_payload(reservation)
        attachment.update(
            {
                "event_type": "job_id_attached",
                "job_id": "1234",
                "state": "submitted",
                "attached_by_control_plane_version": "b" * 64,
                "terminal_recovery_handoff_path": str(handoff_path),
                "terminal_recovery_handoff_sha256": handoff_sha256,
                "terminal_recovery_mode": "fresh_scheduler_binding",
            }
        )
        self.append(attachment)

        second = self.append(
            self._recovery_reservation_event(
                reservation_id="reservation-2", submission_id="submission-2"
            )
        )
        attachment = transition_payload(second)
        attachment.update(
            {
                "event_type": "job_id_attached",
                "job_id": "5678",
                "state": "submitted",
                "attached_by_control_plane_version": "b" * 64,
                "terminal_recovery_handoff_path": str(handoff_path),
                "terminal_recovery_handoff_sha256": handoff_sha256,
                "terminal_recovery_mode": "fresh_scheduler_binding",
            }
        )
        with self.assertRaisesRegex(ValueError, "not bound"):
            self.append(attachment)

    def test_registered_recovery_handoff_rejects_matching_null_security_binding(
        self,
    ) -> None:
        event = self._recovery_reservation_event()
        event["manifest_sha256"] = None
        reservation = self.append(event)
        handoff_path, handoff_sha256 = self._recovery_handoff(reservation)
        attachment = transition_payload(reservation)
        attachment.update(
            {
                "event_type": "job_id_attached",
                "job_id": "1234",
                "state": "submitted",
                "attached_by_control_plane_version": "b" * 64,
                "terminal_recovery_handoff_path": str(handoff_path),
                "terminal_recovery_handoff_sha256": handoff_sha256,
                "terminal_recovery_mode": "fresh_scheduler_binding",
            }
        )
        with self.assertRaisesRegex(ValueError, "not bound"):
            self.append(attachment)

    def test_registered_recovery_handoff_rejects_malformed_purged_snapshot(
        self,
    ) -> None:
        reservation = self.append(self._recovery_reservation_event())
        handoff_path, _ = self._recovery_handoff(
            reservation, mode="purged_scontrol_cancelled_zero_execution"
        )
        for field, value in [
            ("state", "CANCELLEDEVIL"),
            ("elapsed_raw", False),
            ("allocated_nodes", False),
            ("account", "OTHER"),
        ]:
            with self.subTest(field=field):
                handoff_sha256 = self._rewrite_recovery_handoff(
                    handoff_path,
                    lambda handoff, field=field, value=value: handoff[
                        "scheduler_binding_recovery"
                    ].__setitem__(field, value),
                )
                attachment = transition_payload(reservation)
                attachment.update(
                    {
                        "event_type": "job_id_attached",
                        "job_id": "1234",
                        "state": "submitted",
                        "attached_by_control_plane_version": "b" * 64,
                        "terminal_recovery_handoff_path": str(handoff_path),
                        "terminal_recovery_handoff_sha256": handoff_sha256,
                        "terminal_recovery_mode": (
                            "purged_scontrol_cancelled_zero_execution"
                        ),
                    }
                )
                with self.assertRaisesRegex(ValueError, "zero-execution handoff"):
                    self.append(attachment)

    def test_manual_direct_srun_reconciliation_rejects_boolean_billed_nodes(
        self,
    ) -> None:
        event = self._manual_accounting_event()
        event["billed_nodes"] = True
        with self.assertRaisesRegex(ValueError, "usage is invalid"):
            self.append(event)

    def test_append_rejects_unknown_event_type(self) -> None:
        original = self.ledger.read_bytes()
        with self.assertRaisesRegex(ValueError, "missing or unsupported"):
            self.append({"event_type": "opaque"})
        self.assertEqual(self.ledger.read_bytes(), original)

    def test_append_rejects_missing_event_type(self) -> None:
        original = self.ledger.read_bytes()
        with self.assertRaisesRegex(ValueError, "missing or unsupported"):
            self.append({})
        self.assertEqual(self.ledger.read_bytes(), original)

    def test_registered_reservation_rejects_under_bound_schema(self) -> None:
        event = self._reservation_event()
        event.pop("artifact_dir")
        with self.assertRaisesRegex(ValueError, "schema is unsupported"):
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

    def test_read_only_snapshot_rejects_transient_receipt_hide_and_restore(self) -> None:
        hidden = self.receipts.with_name("hidden-receipts.jsonl")
        with self.assertRaisesRegex(ValueError, "parent namespace changed"):
            with validated_read_only_mirrored_state_snapshot(
                self.ledger,
                self.receipts,
                self.mirror,
                ledger_root=self.ledger.parent.parent,
                receipts_root=self.receipts.parent.parent,
                mirror_root=self.mirror.parent.parent,
            ):
                self.receipts.rename(hidden)
                hidden.rename(self.receipts)

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

    def test_read_only_snapshot_rejects_retained_operator_attestation_tamper(
        self,
    ) -> None:
        reservation = self.append(self._reservation_event())
        attestation = Path(str(reservation["pre_manifest_attestation_path"]))
        attestation.chmod(0o600)
        attestation.write_bytes(attestation.read_bytes() + b"\n")
        attestation.chmod(0o400)
        with self.assertRaises(ValueError):
            with validated_read_only_mirrored_state_snapshot(
                self.ledger,
                self.receipts,
                self.mirror,
                ledger_root=self.ledger.parent.parent,
                receipts_root=self.receipts.parent.parent,
                mirror_root=self.mirror.parent.parent,
            ):
                pass

    def test_genesis_anchor_rejects_schema_numeric_aliases_before_equality(
        self,
    ) -> None:
        anchors = genesis_anchor_paths(self.ledger, self.mirror)
        original = {path: path.read_bytes() for path in anchors}
        try:
            for value in [1.0, True]:
                with self.subTest(value=value):
                    anchor = json.loads(original[anchors[0]])
                    anchor["schema_version"] = value
                    payload = (
                        json.dumps(anchor, indent=2, sort_keys=True) + "\n"
                    ).encode("utf-8")
                    for path in anchors:
                        path.chmod(0o600)
                        path.write_bytes(payload)
                        path.chmod(0o444)
                    with self.assertRaisesRegex(ValueError, "genesis-anchor schema"):
                        validate_mirrored_state(
                            self.ledger, self.receipts, self.mirror
                        )
        finally:
            for path, payload in original.items():
                path.chmod(0o600)
                path.write_bytes(payload)
                path.chmod(0o444)

    def test_read_only_snapshot_validates_pinned_bytes_not_transient_paths(self) -> None:
        paths = [self.ledger, self.csv, self.receipts, self.mirror]
        original = {path: path.read_bytes() for path in paths}
        self.append(
            self._reservation_event(
                reservation_id="transient-alternate-reservation"
            )
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
        self.assertEqual(reopened_lengths, [2, 2])

    def test_read_only_snapshot_rejects_recovery_handoff_replacement(self) -> None:
        reservation = self.append(self._recovery_reservation_event())
        handoff_path, handoff_sha256 = self._recovery_handoff(reservation)
        attachment = transition_payload(reservation)
        attachment.update(
            {
                "event_type": "job_id_attached",
                "job_id": "1234",
                "state": "submitted",
                "attached_by_control_plane_version": "b" * 64,
                "terminal_recovery_handoff_path": str(handoff_path),
                "terminal_recovery_handoff_sha256": handoff_sha256,
                "terminal_recovery_mode": "fresh_scheduler_binding",
            }
        )
        self.append(attachment)
        mirror_path = (
            self.root
            / "project_home"
            / "policy"
            / "recovery_handoffs"
            / handoff_path.name
        )
        with self.assertRaisesRegex(ValueError, "snapshot path changed"):
            with validated_read_only_mirrored_state_snapshot(
                self.ledger,
                self.receipts,
                self.mirror,
                ledger_root=self.ledger.parent.parent,
                receipts_root=self.receipts.parent.parent,
                mirror_root=self.mirror.parent.parent,
            ):
                for path in [handoff_path, mirror_path]:
                    replacement = path.with_name(f"replacement-{path.name}")
                    replacement.write_bytes(path.read_bytes())
                    replacement.chmod(0o444)
                    replacement.replace(path)

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

    def test_primary_chain_rejects_rehashed_sequence_numeric_aliases(self) -> None:
        from ledger import canonical_json, record_sha256

        self.append(self._reservation_event())
        canonical_records = [
            json.loads(line)
            for line in self.ledger.read_text(encoding="utf-8").splitlines()
        ]
        for index, value in [(0, 0.0), (1, True)]:
            with self.subTest(index=index, value=value):
                records = [dict(record) for record in canonical_records]
                records[index]["sequence_number"] = value
                previous = ""
                for record in records:
                    record["previous_event_sha256"] = previous
                    record["event_sha256"] = record_sha256(record, "event_sha256")
                    previous = str(record["event_sha256"])
                forged = self.root / f"forged-sequence-{index}.jsonl"
                forged.write_text(
                    "".join(canonical_json(record) + "\n" for record in records),
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(ValueError, "Invalid sequence number"):
                    validate_primary_chain(forged)

    def test_primary_chain_rejects_rehashed_genesis_aliases_and_extra_fields(
        self,
    ) -> None:
        from ledger import canonical_json, record_sha256

        canonical = json.loads(self.ledger.read_text(encoding="utf-8"))
        mutations = [
            ("boolean-version", "control_plane_version", True, "nonempty string"),
            ("numeric-version", "control_plane_version", 1, "nonempty string"),
            ("extra-field", "unsupported_extra", "forged", "schema is unsupported"),
        ]
        for name, field, value, message in mutations:
            with self.subTest(name=name):
                forged_record = dict(canonical)
                forged_record[field] = value
                forged_record["event_sha256"] = record_sha256(
                    forged_record, "event_sha256"
                )
                forged = self.root / f"forged-genesis-{name}.jsonl"
                forged.write_text(
                    canonical_json(forged_record) + "\n", encoding="utf-8"
                )
                with self.assertRaisesRegex(ValueError, message):
                    validate_primary_chain(forged)

    def test_receipts_reject_rehashed_extra_field_and_boolean_ack_timestamp(
        self,
    ) -> None:
        from ledger import canonical_json, record_sha256

        primary = validate_primary_chain(self.ledger)
        canonical = json.loads(self.receipts.read_text(encoding="utf-8"))
        mutations = [
            ("extra-field", "unsupported_extra", "forged"),
            ("boolean-acknowledged-utc", "mirror_acknowledged_utc", True),
        ]
        for name, field, value in mutations:
            with self.subTest(name=name):
                forged_receipt = dict(canonical)
                forged_receipt[field] = value
                forged_receipt["mirror_ack_sha256"] = record_sha256(
                    forged_receipt, "mirror_ack_sha256"
                )
                forged = self.root / f"forged-receipt-{name}.jsonl"
                forged.write_text(
                    canonical_json(forged_receipt) + "\n", encoding="utf-8"
                )
                with self.assertRaisesRegex(ValueError, "receipt schema is invalid"):
                    validate_receipts(
                        forged,
                        primary,
                        mirror_jsonl=self.mirror,
                        mirror_transport="filesystem_copy",
                    )

    def test_repair_copies_missing_mirror_suffix_and_receipt(self) -> None:
        self.append(self._reservation_event())
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
        self.append(self._reservation_event())
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

    def test_repair_rejects_retained_operator_attestation_tamper_before_write(
        self,
    ) -> None:
        reservation = self.append(self._reservation_event())
        self._drop_last_line(self.mirror)
        self._drop_last_line(self.receipts)
        mirror_before = self.mirror.read_bytes()
        receipts_before = self.receipts.read_bytes()
        attestation = Path(str(reservation["pre_manifest_attestation_path"]))
        attestation.chmod(0o600)
        attestation.write_bytes(attestation.read_bytes() + b"\n")
        attestation.chmod(0o400)
        with self.assertRaises(ValueError):
            repair_mirrored_state(
                self.ledger,
                self.csv,
                self.receipts,
                self.mirror,
                mirror_transport="filesystem_copy",
            )
        self.assertEqual(self.mirror.read_bytes(), mirror_before)
        self.assertEqual(self.receipts.read_bytes(), receipts_before)

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
                self.append(self._reservation_event())
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

    def test_interrupted_fresh_genesis_rejects_rehashed_boolean_sequence(self) -> None:
        from ledger import canonical_json, record_sha256

        root = Path(self.temporary.name) / "interrupted-genesis-boolean-sequence"
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
        genesis = json.loads(ledger.read_text(encoding="utf-8"))
        genesis["sequence_number"] = False
        genesis["event_sha256"] = record_sha256(genesis, "event_sha256")
        ledger.write_text(canonical_json(genesis) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Invalid sequence number"):
            initialize_ledger(
                ledger,
                csv_path,
                receipts,
                mirror,
                mirror_transport="filesystem_copy",
                notes="recover exact interrupted genesis",
                control_plane_version="recoverable-control-plane",
            )

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

    def test_audited_pre_anchor_migration_rejects_existing_anchor_schema_alias(
        self,
    ) -> None:
        local_anchor, mirror_anchor = genesis_anchor_paths(self.ledger, self.mirror)
        event = json.loads(self.ledger.read_text(encoding="utf-8"))
        receipt = json.loads(self.receipts.read_text(encoding="utf-8"))
        anchor = json.loads(local_anchor.read_text(encoding="utf-8"))
        anchor["schema_version"] = True
        local_anchor.chmod(0o600)
        local_anchor.write_text(
            json.dumps(anchor, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        local_anchor.chmod(0o444)
        mirror_anchor.chmod(0o600)
        mirror_anchor.unlink()
        with self.assertRaisesRegex(ValueError, "genesis-anchor schema"):
            migrate_existing_genesis_anchors(
                self.ledger,
                self.receipts,
                self.mirror,
                expected_event_sha256=event["event_sha256"],
                expected_mirror_ack_sha256=receipt["mirror_ack_sha256"],
            )


if __name__ == "__main__":
    unittest.main()
