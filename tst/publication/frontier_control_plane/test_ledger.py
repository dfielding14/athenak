#!/usr/bin/env python3
"""Unit tests for the mirrored Frontier PIC ledger."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from ledger import accounting, append_primary_event, initialize_ledger
from ledger import repair_mirrored_state, validate_mirrored_state
from ledger import validate_primary_chain, validate_receipts, write_csv


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

    def _drop_last_line(self, path: Path) -> None:
        lines = path.read_text(encoding="utf-8").splitlines()
        path.write_text("".join(line + "\n" for line in lines[:-1]), encoding="utf-8")

    def test_reserve_attach_reconcile_and_csv_projection(self) -> None:
        common = {
            "reservation_id": "reservation-1",
            "submission_id": "submission-1",
            "submission_scope": "registered_science",
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
        self.assertEqual(rows[-1]["clean_candidate_manifest_sha256"], "a" * 64)

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

    def test_csv_temporary_symlink_alias_rejects_before_external_write(self) -> None:
        temporary = self.csv.with_suffix(self.csv.suffix + ".tmp")
        outside = Path(self.temporary.name) / "outside-csv-temporary"
        outside.write_text("preserve external CSV target\n", encoding="utf-8")
        temporary.symlink_to(outside)
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


if __name__ == "__main__":
    unittest.main()
