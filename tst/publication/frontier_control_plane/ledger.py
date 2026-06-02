#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Hash-chained PIC node-hour ledger with a durable non-recursive mirror."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

from contextlib import contextmanager
from contextvars import ContextVar
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import io
import json
import math
import os
from pathlib import Path
import stat
from typing import Iterator, TextIO
import uuid

from control_plane_common import atomic_write_bytes, atomic_write_bytes_at
from control_plane_common import atomic_write_json, atomic_write_json_at
from control_plane_common import durable_mkdir_parents, fsync_directory
from control_plane_common import open_directory_below, read_json_bytes
from control_plane_common import read_stable_regular_file_below
from control_plane_common import scheduler_account_matches_authorized
from control_plane_common import stable_serialization_anchor
from operator_attestation import validate_sealed_operator_attestation


CSV_FIELDS = [
    "sequence_number",
    "previous_event_sha256",
    "event_sha256",
    "event_type",
    "timestamp",
    "reservation_id",
    "submission_id",
    "job_id",
    "control_plane_version",
    "attached_by_control_plane_version",
    "reconciled_by_control_plane_version",
    "terminal_recovery_handoff_path",
    "terminal_recovery_handoff_sha256",
    "terminal_recovery_mode",
    "manual_accounting_authorization_id",
    "manual_accounting_authorization_path",
    "manual_accounting_project_home_authorization_path",
    "manual_accounting_authorization_sha256",
    "accounting_scope",
    "scientific_evidence_eligible",
    "active_policy_sha256",
    "active_promotion_sha256",
    "git_commit",
    "campaign",
    "test_id",
    "submission_scope",
    "registered_science_authorization_id",
    "clean_candidate_manifest_sha256",
    "pre_manifest_attestation_path",
    "pre_manifest_attestation_sha256",
    "pre_submit_wrapper_attestation_path",
    "pre_submit_wrapper_attestation_sha256",
    "partition",
    "qos",
    "qos_selection_reason",
    "queue_snapshot_sha256",
    "manifest_path",
    "manifest_sha256",
    "job_script_sha256",
    "executable_sha256",
    "site_policy_checked_utc",
    "requested_nodes",
    "scheduler_reported_allocated_nodes",
    "billed_nodes",
    "requested_walltime",
    "reserved_node_hours",
    "elapsed_seconds",
    "consumed_node_hours",
    "cumulative_consumed_node_hours",
    "state",
    "reconciled",
    "artifact_dir",
    "mirror_destination",
    "mirror_transport",
    "mirror_acknowledged_utc",
    "mirror_ack_sha256",
    "notes",
]
GENESIS_ANCHOR_FILENAME = "genesis_anchor.json"
INCOMPLETE_MANUAL_ACCOUNTING_MARKER_FILENAME = "pending_manual_accounting.json"
MANUAL_DIRECT_SRUN_ACCOUNTING_SCOPE = "manual_direct_srun_accounting_only"
MANUAL_DIRECT_SCHEDULER_ACCOUNTING_SCOPE = (
    "manual_direct_scheduler_accounting_only"
)
MANUAL_ACCOUNTING_SCOPE_NOTES = {
    MANUAL_DIRECT_SRUN_ACCOUNTING_SCOPE: (
        "Reviewed direct-srun accounting only; ineligible for scientific evidence."
    ),
    MANUAL_DIRECT_SCHEDULER_ACCOUNTING_SCOPE: (
        "Reviewed direct-scheduler accounting only; ineligible for scientific evidence."
    ),
}
MANUAL_ACCOUNTING_TERMINAL_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "COMPLETED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "TIMEOUT",
}
MANUAL_ACCOUNTING_EVENT_FIELDS = {
    "sequence_number",
    "previous_event_sha256",
    "event_sha256",
    "timestamp",
    "event_type",
    "job_id",
    "control_plane_version",
    "reconciled_by_control_plane_version",
    "manual_accounting_authorization_id",
    "manual_accounting_authorization_path",
    "manual_accounting_project_home_authorization_path",
    "manual_accounting_authorization_sha256",
    "accounting_scope",
    "scientific_evidence_eligible",
    "active_policy_sha256",
    "active_promotion_sha256",
    "partition",
    "qos",
    "scheduler_reported_allocated_nodes",
    "billed_nodes",
    "elapsed_seconds",
    "consumed_node_hours",
    "cumulative_consumed_node_hours",
    "state",
    "reconciled",
    "notes",
}
CHAIN_EVENT_FIELDS = {
    "sequence_number",
    "previous_event_sha256",
    "event_sha256",
    "timestamp",
    "event_type",
}
GENESIS_EVENT_FIELDS = CHAIN_EVENT_FIELDS | {
    "state",
    "notes",
    "control_plane_version",
}
RESERVATION_PAYLOAD_BASE_FIELDS = {
    "reservation_id",
    "submission_id",
    "control_plane_version",
    "git_commit",
    "campaign",
    "test_id",
    "partition",
    "qos",
    "qos_selection_reason",
    "queue_snapshot_sha256",
    "site_policy_checked_utc",
    "requested_nodes",
    "requested_walltime",
    "reserved_node_hours",
    "artifact_dir",
    "state",
    "reconciled",
}
RESERVATION_PAYLOAD_BOUND_FIELDS = {
    "active_policy_sha256",
    "active_promotion_sha256",
    "clean_candidate_manifest_sha256",
    "executable_sha256",
    "job_script_sha256",
    "manifest_path",
    "manifest_sha256",
    "submission_scope",
}
RESERVATION_PAYLOAD_OPTIONAL_FIELDS = RESERVATION_PAYLOAD_BOUND_FIELDS | {
    "registered_science_authorization_id",
    "pre_manifest_attestation_path",
    "pre_manifest_attestation_sha256",
    "pre_submit_wrapper_attestation_path",
    "pre_submit_wrapper_attestation_sha256",
}
OPERATOR_ATTESTATION_PROVENANCE_FIELDS = {
    "pre_manifest_attestation_path",
    "pre_manifest_attestation_sha256",
    "pre_submit_wrapper_attestation_path",
    "pre_submit_wrapper_attestation_sha256",
}
LEGACY_OPERATOR_ATTESTATION_FREE_CONTROL_PLANE_VERSIONS = {
    "6002c80e305d6cfd322675b3e27e17e169b1718d8edd722330464cccb0c4fd86",
    "4ccde8dfe557fbb6450c15b2ee68d41d2ac0711d5d56d79e775fa75178177bd3",
    "3d3d0d20ab3cedb1b650d8a19e99a0123b550ead25ebd312882d03d3be3b41d9",
}
RECONCILIATION_PAYLOAD_FIELDS = {
    "scheduler_reported_allocated_nodes",
    "billed_nodes",
    "elapsed_seconds",
    "consumed_node_hours",
    "cumulative_consumed_node_hours",
}
MIRROR_RECEIPT_FIELDS = {
    "mirrored_event_sha256",
    "mirror_destination",
    "mirror_transport",
    "mirror_acknowledged_utc",
    "mirror_ack_sha256",
}
REGISTERED_RESERVATION_EVENT_TYPES = {
    "reservation",
    "job_id_attached",
    "reservation_cancelled",
    "reconciliation",
}
SCHEDULER_JOB_OWNERSHIP_EVENT_TYPES = {
    "job_id_attached",
    "manual_allocation_reconciliation",
}
TERMINAL_RECOVERY_FIELDS = {
    "terminal_recovery_handoff_path",
    "terminal_recovery_handoff_sha256",
    "terminal_recovery_mode",
}
TERMINAL_RECOVERY_MODES = {
    "fresh_scheduler_binding",
    "purged_scontrol_cancelled_zero_execution",
}
REGISTERED_ATTACHMENT_CHANGES = {
    "state",
    "job_id",
    "attached_by_control_plane_version",
    "terminal_recovery_handoff_path",
    "terminal_recovery_handoff_sha256",
    "terminal_recovery_mode",
    "notes",
}
REGISTERED_CANCELLATION_CHANGES = {
    "state",
    "notes",
}
REGISTERED_RECONCILIATION_CHANGES = {
    "state",
    "reconciled",
    "reconciled_by_control_plane_version",
    "scheduler_reported_allocated_nodes",
    "billed_nodes",
    "elapsed_seconds",
    "consumed_node_hours",
    "cumulative_consumed_node_hours",
}
_PINNED_PARENT_DESCRIPTORS: ContextVar[dict[Path, int]] = ContextVar(
    "_PINNED_PARENT_DESCRIPTORS", default={}
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def canonical_json(record: dict[str, object], omit: str | None = None) -> str:
    payload = {key: value for key, value in record.items() if key != omit}
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def record_sha256(record: dict[str, object], omit: str) -> str:
    return hashlib.sha256(canonical_json(record, omit).encode("utf-8")).hexdigest()


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"Non-finite JSON number is not allowed: {value}")


def _reject_duplicate_json_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON object key is not allowed: {key}")
        value[key] = item
    return value


def _parent_descriptor(path: Path) -> int | None:
    return _PINNED_PARENT_DESCRIPTORS.get().get(Path(os.path.abspath(path.parent)))


def _read_regular_bytes(path: Path, *, require_read_only_mode: bool = False) -> bytes:
    parent_descriptor = _parent_descriptor(path)
    descriptor = os.open(
        path.name if parent_descriptor is not None else path,
        os.O_RDONLY | os.O_NOFOLLOW,
        dir_fd=parent_descriptor,
    )
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"Expected a regular file: {path}")
        if require_read_only_mode and metadata.st_mode & 0o222:
            raise ValueError(f"Expected a read-only regular file: {path}")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            return stream.read()
    finally:
        os.close(descriptor)


def _path_exists(path: Path) -> bool:
    parent_descriptor = _parent_descriptor(path)
    if parent_descriptor is None:
        return path.exists()
    try:
        os.stat(path.name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        return False
    return True


def _read_jsonl(path: Path, *, root: Path | None = None) -> list[dict[str, object]]:
    try:
        if root is not None:
            data = read_stable_regular_file_below(path, root)
        else:
            data = _read_regular_bytes(path)
    except FileNotFoundError:
        return []
    return _read_jsonl_bytes(data, path=path)


def _read_jsonl_bytes(data: bytes, *, path: Path) -> list[dict[str, object]]:
    records = []
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"JSONL file is not UTF-8: {path}") from error
    for number, raw in enumerate(text.splitlines(), 1):
        if not raw:
            raise ValueError(f"Blank JSONL line in {path}:{number}")
        record = json.loads(
            raw,
            parse_constant=_reject_json_constant,
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
        if not isinstance(record, dict):
            raise ValueError(f"JSONL record is not an object in {path}:{number}")
        records.append(record)
    return records


def _append_jsonl(path: Path, record: dict[str, object]) -> None:
    parent_descriptor = _parent_descriptor(path)
    if parent_descriptor is None:
        durable_mkdir_parents(path.parent)
    created = False
    flags = os.O_APPEND | os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_NOFOLLOW
    try:
        descriptor = os.open(
            path.name if parent_descriptor is not None else path,
            flags,
            0o600,
            dir_fd=parent_descriptor,
        )
        created = True
    except FileExistsError:
        descriptor = os.open(
            path.name if parent_descriptor is not None else path,
            os.O_APPEND | os.O_WRONLY | os.O_NOFOLLOW,
            dir_fd=parent_descriptor,
        )
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise ValueError(f"Expected a regular JSONL file: {path}")
        stream = os.fdopen(descriptor, "a", encoding="utf-8")
    except BaseException:
        os.close(descriptor)
        raise
    with stream:
        stream.write(canonical_json(record) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    if created:
        if parent_descriptor is not None:
            os.fsync(parent_descriptor)
        else:
            fsync_directory(path.parent)


def validate_primary_chain(
    path: Path, *, root: Path | None = None
) -> list[dict[str, object]]:
    return _validate_primary_records(_read_jsonl(path, root=root), path=path)


def _validate_primary_records(
    records: list[dict[str, object]], *, path: Path
) -> list[dict[str, object]]:
    previous = ""
    for expected_sequence, record in enumerate(records):
        if (
            type(record.get("sequence_number")) is not int
            or record.get("sequence_number") != expected_sequence
        ):
            raise ValueError(f"Invalid sequence number in {path}")
        if record.get("previous_event_sha256") != previous:
            raise ValueError(f"Broken previous-event link in {path}")
        expected_hash = record_sha256(record, "event_sha256")
        if record.get("event_sha256") != expected_hash:
            raise ValueError(f"Invalid event hash in {path}")
        previous = expected_hash
    _validate_accounting_records(records)
    return records


def _validate_receipt_provenance(
    receipt: dict[str, object],
    *,
    mirror_jsonl: Path,
    mirror_transport: str,
) -> None:
    if (
        set(receipt) != MIRROR_RECEIPT_FIELDS
        or not _is_lowercase_sha256(receipt.get("mirrored_event_sha256"))
        or type(receipt.get("mirror_destination")) is not str
        or not receipt["mirror_destination"]
        or type(receipt.get("mirror_transport")) is not str
        or not receipt["mirror_transport"]
        or type(receipt.get("mirror_acknowledged_utc")) is not str
        or not receipt["mirror_acknowledged_utc"]
        or not _is_lowercase_sha256(receipt.get("mirror_ack_sha256"))
    ):
        raise ValueError("Mirror receipt schema is invalid")
    if receipt.get("mirror_destination") != str(mirror_jsonl):
        raise ValueError("Mirror receipt destination does not match configured mirror")
    if receipt.get("mirror_transport") != mirror_transport:
        raise ValueError("Mirror receipt transport does not match configured transport")


def validate_receipts(
    path: Path,
    primary_records: list[dict[str, object]],
    *,
    mirror_jsonl: Path,
    mirror_transport: str,
    root: Path | None = None,
) -> list[dict[str, object]]:
    return _validate_receipt_records(
        _read_jsonl(path, root=root),
        primary_records,
        path=path,
        mirror_jsonl=mirror_jsonl,
        mirror_transport=mirror_transport,
    )


def _validate_receipt_records(
    records: list[dict[str, object]],
    primary_records: list[dict[str, object]],
    *,
    path: Path,
    mirror_jsonl: Path,
    mirror_transport: str,
) -> list[dict[str, object]]:
    primary_hashes = [record["event_sha256"] for record in primary_records]
    for record in records:
        _validate_receipt_provenance(
            record,
            mirror_jsonl=mirror_jsonl,
            mirror_transport=mirror_transport,
        )
        if record.get("mirrored_event_sha256") not in primary_hashes:
            raise ValueError(f"Mirror receipt references an unknown event in {path}")
        expected_hash = record_sha256(record, "mirror_ack_sha256")
        if record.get("mirror_ack_sha256") != expected_hash:
            raise ValueError(f"Invalid mirror receipt hash in {path}")
    receipt_hashes = [record.get("mirrored_event_sha256") for record in records]
    if receipt_hashes != primary_hashes:
        raise ValueError(f"Mirror receipts must acknowledge each primary event exactly once in {path}")
    return records


def chain_head(records: list[dict[str, object]]) -> str:
    return str(records[-1]["event_sha256"]) if records else ""


def require_explicit_genesis(records: list[dict[str, object]]) -> None:
    if not records:
        raise ValueError("Frontier PIC ledger has no explicit genesis event")
    genesis = records[0]
    if (
        genesis.get("event_type") != "genesis"
        or genesis.get("state") != "initialized"
        or type(genesis.get("control_plane_version")) is not str
        or not genesis["control_plane_version"]
    ):
        raise ValueError("Frontier PIC ledger does not start with a valid genesis event")
    if any(record.get("event_type") == "genesis" for record in records[1:]):
        raise ValueError("Frontier PIC ledger contains more than one genesis event")


def _nonnegative_finite_number(
    record: dict[str, object], field: str, *, label: str
) -> float:
    value = record.get(field)
    if type(value) not in {int, float}:
        raise ValueError(f"{label} {field} must be a real number")
    try:
        result = float(value)
    except OverflowError as error:
        raise ValueError(f"{label} {field} must be finite") from error
    if not math.isfinite(result) or result < 0.0:
        raise ValueError(f"{label} {field} must be finite and non-negative")
    return result


def _is_lowercase_sha256(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _registered_event_field_bounds(event_type: str) -> tuple[set[str], set[str]]:
    base = CHAIN_EVENT_FIELDS | RESERVATION_PAYLOAD_BASE_FIELDS
    optional = RESERVATION_PAYLOAD_OPTIONAL_FIELDS
    if event_type == "reservation":
        return base, base | optional
    if event_type == "reservation_cancelled":
        return base, base | optional | {"notes"}
    if event_type == "job_id_attached":
        required = base | {"job_id"}
        return required, (
            required
            | optional
            | {"attached_by_control_plane_version"}
            | TERMINAL_RECOVERY_FIELDS
        )
    if event_type == "reconciliation":
        required = base | {"job_id"} | RECONCILIATION_PAYLOAD_FIELDS
        return required, (
            required
            | optional
            | {
                "attached_by_control_plane_version",
                "reconciled_by_control_plane_version",
            }
            | TERMINAL_RECOVERY_FIELDS
        )
    raise ValueError(f"Unsupported registered ledger event type: {event_type}")


def _require_closed_primary_event_schema(record: dict[str, object]) -> None:
    event_type = record.get("event_type")
    if type(event_type) is not str or event_type not in {
        "genesis",
        "reservation",
        "job_id_attached",
        "reservation_cancelled",
        "reconciliation",
        "manual_allocation_reconciliation",
    }:
        raise ValueError("Ledger event type is missing or unsupported")
    if type(record.get("timestamp")) is not str or not record["timestamp"]:
        raise ValueError("Ledger event timestamp must be a nonempty string")
    if event_type == "genesis":
        required_fields = allowed_fields = GENESIS_EVENT_FIELDS
    elif event_type == "manual_allocation_reconciliation":
        required_fields = allowed_fields = MANUAL_ACCOUNTING_EVENT_FIELDS
    else:
        required_fields, allowed_fields = _registered_event_field_bounds(event_type)
    if not required_fields <= set(record) or not set(record) <= allowed_fields:
        raise ValueError(f"{event_type} ledger event schema is unsupported")
    if event_type == "genesis" and (
        type(record.get("control_plane_version")) is not str
        or not record["control_plane_version"]
    ):
        raise ValueError("Genesis control-plane version must be a nonempty string")


def slurm_walltime_seconds(value: object) -> int:
    if not isinstance(value, str) or not value:
        raise ValueError("Registered reservation walltime is invalid")
    days = 0
    time_value = value
    has_days = "-" in value
    try:
        if has_days:
            day_text, time_value = value.split("-", 1)
            days = int(day_text)
        fields = [int(field) for field in time_value.split(":")]
    except ValueError as error:
        raise ValueError("Registered reservation walltime is invalid") from error
    if has_days and len(fields) == 1:
        hours = fields[0]
        minutes = 0
        seconds = 0
    elif has_days and len(fields) == 2:
        hours, minutes = fields
        seconds = 0
    elif len(fields) == 1:
        hours = 0
        minutes = fields[0]
        seconds = 0
    elif len(fields) == 2:
        hours = 0
        minutes, seconds = fields
    elif len(fields) == 3:
        hours, minutes, seconds = fields
    else:
        raise ValueError("Registered reservation walltime is invalid")
    if (
        days < 0
        or hours < 0
        or minutes < 0
        or seconds not in range(60)
        or ((has_days or len(fields) == 3) and minutes not in range(60))
    ):
        raise ValueError("Registered reservation walltime is invalid")
    seconds = (((days * 24) + hours) * 60 + minutes) * 60 + seconds
    return ((seconds + 59) // 60) * 60


def _accounting_without_validation(
    records: list[dict[str, object]],
) -> dict[str, float]:
    consumed = sum(
        _nonnegative_finite_number(
            record, "consumed_node_hours", label="Manual-accounting ledger event"
        )
        for record in records
        if record.get("event_type") == "manual_allocation_reconciliation"
    )
    reserved = 0.0
    for record in latest_reservations(records).values():
        if record.get("reconciled", False) is True:
            consumed += _nonnegative_finite_number(
                record, "consumed_node_hours", label="Registered reconciliation"
            )
        elif record.get("state") not in {"cancelled", "submission_attach_failed"}:
            reserved += _nonnegative_finite_number(
                record, "reserved_node_hours", label="Ledger reservation"
            )
    return {
        "cumulative_consumed_node_hours": consumed,
        "currently_reserved_node_hours": reserved,
    }


def _require_registered_transition_payload(
    prior: dict[str, object],
    record: dict[str, object],
    *,
    allowed_changes: set[str],
) -> None:
    prior_payload = transition_payload(prior)
    payload = transition_payload(record)
    if set(payload) - set(prior_payload) - allowed_changes:
        raise ValueError("Registered reservation transition adds unexpected fields")
    for key, value in prior_payload.items():
        if (
            key not in allowed_changes
            and (
                key not in payload
                or canonical_json({"value": payload[key]})
                != canonical_json({"value": value})
            )
        ):
            raise ValueError("Registered reservation transition rewrites immutable fields")


def _require_canonical_scheduler_job_id(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value.isascii()
        or not value.isdigit()
        or value.startswith("0")
    ):
        raise ValueError("Scheduler job ID is not canonical")
    return value


def _require_new_scheduler_job_id(
    records: list[dict[str, object]], job_id: str
) -> None:
    if any(
        earlier.get("event_type") in SCHEDULER_JOB_OWNERSHIP_EVENT_TYPES
        and earlier.get("job_id") == job_id
        for earlier in records
    ):
        raise ValueError("Scheduler job ID already has a ledger owner")


def _require_terminal_recovery_provenance(
    record: dict[str, object], *, prior: dict[str, object] | None = None
) -> None:
    present = TERMINAL_RECOVERY_FIELDS & set(record)
    if present and present != TERMINAL_RECOVERY_FIELDS:
        raise ValueError("Terminal-recovery provenance is incomplete")
    if present:
        path = record["terminal_recovery_handoff_path"]
        handoff_path = Path(str(path))
        manifest_path = Path(str(record.get("manifest_path", "")))
        try:
            canonical_handoff_id = str(uuid.UUID(handoff_path.stem))
        except ValueError:
            canonical_handoff_id = ""
        try:
            pic_root = manifest_path.parents[3]
        except IndexError:
            pic_root = Path("")
        if (
            not isinstance(path, str)
            or not os.path.isabs(path)
            or path != os.path.abspath(path)
            or not manifest_path.is_absolute()
            or manifest_path.parents[2] != pic_root / "manifests"
            or handoff_path.parent != pic_root / "policy" / "recovery_handoffs"
            or handoff_path.suffix != ".json"
            or handoff_path.stem != canonical_handoff_id
            or not _is_lowercase_sha256(record["terminal_recovery_handoff_sha256"])
            or record["terminal_recovery_mode"] not in TERMINAL_RECOVERY_MODES
        ):
            raise ValueError("Terminal-recovery provenance is invalid")
    if prior is not None and TERMINAL_RECOVERY_FIELDS & set(prior):
        if any(record.get(field) != prior[field] for field in TERMINAL_RECOVERY_FIELDS):
            raise ValueError("Terminal-recovery provenance is not preserved")


def _require_purged_zero_execution_semantics(record: dict[str, object]) -> None:
    if record.get("terminal_recovery_mode") == "purged_scontrol_cancelled_zero_execution":
        if (
            record.get("state") != "CANCELLED"
            or type(record.get("scheduler_reported_allocated_nodes")) is not int
            or record.get("scheduler_reported_allocated_nodes") != 0
            or type(record.get("elapsed_seconds")) is not int
            or record.get("elapsed_seconds") != 0
            or type(record.get("consumed_node_hours")) not in {int, float}
            or record.get("consumed_node_hours") != 0.0
        ):
            raise ValueError("Purged zero-execution reconciliation semantics differ")


def _require_operator_attestation_provenance(record: dict[str, object]) -> None:
    present = OPERATOR_ATTESTATION_PROVENANCE_FIELDS & set(record)
    if present and present != OPERATOR_ATTESTATION_PROVENANCE_FIELDS:
        raise ValueError("Registered operator-attestation provenance quartet is incomplete")
    if not present:
        if (
            record.get("submission_scope") == "registered_science"
            and record.get("control_plane_version")
            not in LEGACY_OPERATOR_ATTESTATION_FREE_CONTROL_PLANE_VERSIONS
        ):
            raise ValueError(
                "Registered science requires operator-attestation provenance"
            )
        return
    if record.get("submission_scope") != "registered_science":
        raise ValueError("Only registered science may carry operator-attestation provenance")
    for phase in ["pre_manifest", "pre_submit_wrapper"]:
        path = record[f"{phase}_attestation_path"]
        digest = record[f"{phase}_attestation_sha256"]
        if (
            not isinstance(path, str)
            or not os.path.isabs(path)
            or Path(path).name != "attestation.json"
            or Path(path).parent.parent.name != "operator_attestations"
            or not Path(path).parent.name.endswith(f"-{phase}")
            or not _is_lowercase_sha256(digest)
        ):
            raise ValueError("Registered operator-attestation provenance is malformed")


def _validate_operator_attestation_trees(
    records: list[dict[str, object]],
    *,
    ledger_jsonl: Path,
    mirror_jsonl: Path,
) -> None:
    pic_root = Path(os.path.abspath(ledger_jsonl.parent.parent))
    project_home_root = Path(os.path.abspath(mirror_jsonl.parent.parent))
    validated: set[tuple[str, str, str, str, str]] = set()
    for record in records:
        if not OPERATOR_ATTESTATION_PROVENANCE_FIELDS & set(record):
            continue
        authorization_id = record.get("registered_science_authorization_id")
        control_plane_version = record.get("control_plane_version")
        if not isinstance(authorization_id, str) or not authorization_id:
            raise ValueError("Registered operator-attestation authorization ID is invalid")
        if not isinstance(control_plane_version, str):
            raise ValueError("Registered operator-attestation control-plane version is invalid")
        for phase in ["pre_manifest", "pre_submit_wrapper"]:
            path = str(record[f"{phase}_attestation_path"])
            digest = str(record[f"{phase}_attestation_sha256"])
            key = (path, digest, authorization_id, control_plane_version, phase)
            if key in validated:
                continue
            binding = validate_sealed_operator_attestation(
                Path(path),
                authorization_id=authorization_id,
                phase=phase,
                control_plane_version=control_plane_version,
                authorized_pic_root=pic_root,
                authorized_project_home_root=project_home_root,
                enforce_freshness=False,
            )
            if binding != {"path": path, "sha256": digest}:
                raise ValueError("Registered operator-attestation binding differs")
            validated.add(key)


def _validate_terminal_recovery_handoff_mirrors(
    records: list[dict[str, object]],
    *,
    ledger_jsonl: Path,
    mirror_jsonl: Path,
    state: dict[Path, bytes] | None = None,
) -> None:
    pic_root = Path(os.path.abspath(ledger_jsonl.parent.parent))
    project_home_root = Path(os.path.abspath(mirror_jsonl.parent.parent))
    cached_bytes: dict[Path, tuple[bytes, bytes]] = {}
    for record in records:
        if not TERMINAL_RECOVERY_FIELDS & set(record):
            continue
        handoff_path = Path(str(record["terminal_recovery_handoff_path"]))
        if handoff_path.parent != pic_root / "policy" / "recovery_handoffs":
            raise ValueError("Terminal-recovery handoff is outside the PIC root")
        mirror_path = project_home_root / "policy" / "recovery_handoffs" / handoff_path.name
        if handoff_path not in cached_bytes:
            try:
                if state is None:
                    handoff_bytes = read_stable_regular_file_below(
                        handoff_path, pic_root, require_read_only_mode=True
                    )
                    mirror_bytes = read_stable_regular_file_below(
                        mirror_path, project_home_root, require_read_only_mode=True
                    )
                else:
                    handoff_bytes = state[handoff_path]
                    mirror_bytes = state[mirror_path]
            except (FileNotFoundError, KeyError) as error:
                raise ValueError("Terminal-recovery handoff mirror is missing") from error
            cached_bytes[handoff_path] = (handoff_bytes, mirror_bytes)
        handoff_bytes, mirror_bytes = cached_bytes[handoff_path]
        if handoff_bytes != mirror_bytes:
            raise ValueError("Terminal-recovery handoff mirror bytes differ")
        if hashlib.sha256(handoff_bytes).hexdigest() != record[
            "terminal_recovery_handoff_sha256"
        ]:
            raise ValueError("Terminal-recovery handoff SHA-256 differs")
        _validate_terminal_recovery_handoff_bytes(record, handoff_path, handoff_bytes)


def _validate_terminal_recovery_handoff_bytes(
    record: dict[str, object], handoff_path: Path, handoff_bytes: bytes
) -> None:
    handoff = read_json_bytes(handoff_bytes, label=str(handoff_path))
    base_fields = {
        "schema_version",
        "status",
        "handoff_id",
        "reservation_id",
        "submission_id",
        "job_id",
        "manifest_path",
        "manifest_sha256",
        "prior_control_plane_version",
        "recovery_control_plane_version",
        "prior_active_policy_sha256",
        "prior_active_promotion_sha256",
        "pending_marker_sha256",
    }
    mode = record["terminal_recovery_mode"]
    extra_fields = (
        {"scheduler_binding_recovery", "reservation_job_binding_attestation"}
        if mode == "purged_scontrol_cancelled_zero_execution"
        else set()
    )
    if (
        set(handoff) != base_fields | extra_fields
        or type(handoff.get("schema_version")) is not int
        or handoff.get("schema_version") != 1
        or handoff.get("status")
        != "authorized_terminal_scheduler_job_id_received_recovery"
        or handoff.get("handoff_id") != handoff_path.stem
        or not isinstance(handoff.get("reservation_id"), str)
        or not handoff["reservation_id"]
        or handoff.get("reservation_id") != record.get("reservation_id")
        or not isinstance(handoff.get("submission_id"), str)
        or not handoff["submission_id"]
        or handoff.get("submission_id") != record.get("submission_id")
        or not isinstance(handoff.get("job_id"), str)
        or not handoff["job_id"]
        or handoff.get("job_id") != record.get("job_id")
        or not isinstance(handoff.get("manifest_path"), str)
        or not os.path.isabs(handoff["manifest_path"])
        or handoff.get("manifest_path") != record.get("manifest_path")
        or not _is_lowercase_sha256(handoff.get("manifest_sha256"))
        or handoff.get("manifest_sha256") != record.get("manifest_sha256")
        or not _is_lowercase_sha256(handoff.get("prior_control_plane_version"))
        or handoff.get("prior_control_plane_version")
        != record.get("control_plane_version")
        or not _is_lowercase_sha256(handoff.get("recovery_control_plane_version"))
        or handoff.get("recovery_control_plane_version")
        != record.get("attached_by_control_plane_version")
        or not _is_lowercase_sha256(handoff.get("prior_active_policy_sha256"))
        or handoff.get("prior_active_policy_sha256")
        != record.get("active_policy_sha256")
        or not _is_lowercase_sha256(handoff.get("prior_active_promotion_sha256"))
        or handoff.get("prior_active_promotion_sha256")
        != record.get("active_promotion_sha256")
        or not _is_lowercase_sha256(handoff.get("pending_marker_sha256"))
    ):
        raise ValueError("Terminal-recovery handoff bytes are not bound to the ledger")
    if mode == "fresh_scheduler_binding":
        return
    snapshot = handoff["scheduler_binding_recovery"]
    attestation = handoff["reservation_job_binding_attestation"]
    if (
        not isinstance(snapshot, dict)
        or set(snapshot)
        != {
            "mode",
            "job_id",
            "job_name",
            "state",
            "elapsed_raw",
            "allocated_nodes",
            "comment",
            "account",
            "submit",
            "start",
            "end",
            "exit_code",
        }
        or snapshot.get("mode") != mode
        or snapshot.get("job_id") != record.get("job_id")
        or snapshot.get("job_name") != "run_installed_control_plane_job.sh"
        or not isinstance(snapshot.get("state"), str)
        or not snapshot["state"].split()
        or snapshot["state"].split()[0] != "CANCELLED"
        or type(snapshot.get("elapsed_raw")) is not int
        or snapshot.get("elapsed_raw") != 0
        or type(snapshot.get("allocated_nodes")) is not int
        or snapshot.get("allocated_nodes") != 0
        or snapshot.get("comment") != ""
        or not scheduler_account_matches_authorized(snapshot.get("account"))
        or not isinstance(snapshot.get("submit"), str)
        or not snapshot["submit"]
        or snapshot.get("start") != "None"
        or snapshot.get("end") != snapshot["submit"]
        or snapshot.get("exit_code") != "0:0"
        or attestation
        != {
            "mode": (
                "reviewed_operator_attestation_for_unprovable_purged_"
                "reservation_job_binding"
            ),
            "reservation_id": record.get("reservation_id"),
            "job_id": record.get("job_id"),
        }
    ):
        raise ValueError("Purged zero-execution handoff bytes are invalid")


def _validate_accounting_records(records: list[dict[str, object]]) -> None:
    for index, record in enumerate(records):
        _require_closed_primary_event_schema(record)
        _require_operator_attestation_provenance(record)
        event_type = record.get("event_type")
        if "reconciled" in record and type(record["reconciled"]) is not bool:
            raise ValueError("Ledger reconciled status must be boolean")
        if "reserved_node_hours" in record:
            _nonnegative_finite_number(
                record, "reserved_node_hours", label="Ledger reservation"
            )
        if (
            "reservation_id" in record
            or event_type in REGISTERED_RESERVATION_EVENT_TYPES
        ):
            reservation_id = record.get("reservation_id")
            if not isinstance(reservation_id, str) or not reservation_id:
                raise ValueError("Registered reservation identifier is invalid")
            prior = latest_reservations(records[:index]).get(reservation_id)
            if event_type == "reservation":
                requested_nodes = record.get("requested_nodes")
                walltime_seconds = slurm_walltime_seconds(
                    record.get("requested_walltime")
                )
                reserved = _nonnegative_finite_number(
                    record, "reserved_node_hours", label="Registered reservation"
                )
                if (
                    prior is not None
                    or record.get("state") != "reserved"
                    or record.get("reconciled") is not False
                    or type(requested_nodes) is not int
                    or requested_nodes <= 0
                    or walltime_seconds <= 0
                    or reserved != requested_nodes * walltime_seconds / 3600.0
                    or TERMINAL_RECOVERY_FIELDS & set(record)
                ):
                    raise ValueError("Registered reservation event is invalid")
            elif event_type == "job_id_attached":
                job_id = _require_canonical_scheduler_job_id(record.get("job_id"))
                _require_new_scheduler_job_id(records[:index], job_id)
                if (
                    prior is None
                    or prior.get("event_type") != "reservation"
                    or prior.get("state") != "reserved"
                    or prior.get("reconciled") is not False
                    or record.get("state") != "submitted"
                    or record.get("reconciled") is not False
                ):
                    raise ValueError("Registered attachment transition is invalid")
                _require_registered_transition_payload(
                    prior, record, allowed_changes=REGISTERED_ATTACHMENT_CHANGES
                )
                _require_terminal_recovery_provenance(record, prior=prior)
            elif event_type == "reservation_cancelled":
                if (
                    prior is None
                    or prior.get("event_type") != "reservation"
                    or prior.get("state") != "reserved"
                    or prior.get("reconciled") is not False
                    or record.get("state") != "cancelled"
                    or record.get("reconciled") is not False
                ):
                    raise ValueError("Registered cancellation transition is invalid")
                _require_registered_transition_payload(
                    prior, record, allowed_changes=REGISTERED_CANCELLATION_CHANGES
                )
            elif event_type != "reconciliation":
                raise ValueError("Unknown reservation-bearing ledger event")
        if event_type == "manual_allocation_reconciliation":
            job_id = _require_canonical_scheduler_job_id(record.get("job_id"))
            _require_new_scheduler_job_id(records[:index], job_id)
            if (
                set(record) != MANUAL_ACCOUNTING_EVENT_FIELDS
                or record.get("reconciled") is not True
                or not isinstance(record.get("accounting_scope"), str)
                or record["accounting_scope"] not in MANUAL_ACCOUNTING_SCOPE_NOTES
                or record.get("scientific_evidence_eligible") is not False
                or not isinstance(record.get("job_id"), str)
                or not record["job_id"]
                or not _is_lowercase_sha256(record.get("control_plane_version"))
                or record.get("reconciled_by_control_plane_version")
                != record["control_plane_version"]
                or not isinstance(record.get("manual_accounting_authorization_id"), str)
                or not record["manual_accounting_authorization_id"]
                or not isinstance(record.get("manual_accounting_authorization_path"), str)
                or not os.path.isabs(str(record["manual_accounting_authorization_path"]))
                or not isinstance(
                    record.get("manual_accounting_project_home_authorization_path"), str
                )
                or not os.path.isabs(
                    str(record["manual_accounting_project_home_authorization_path"])
                )
                or not _is_lowercase_sha256(
                    record.get("manual_accounting_authorization_sha256")
                )
                or not _is_lowercase_sha256(record.get("active_policy_sha256"))
                or not _is_lowercase_sha256(record.get("active_promotion_sha256"))
                or record.get("partition") != "batch"
                or record.get("qos") not in {"debug", "normal"}
                or record.get("state") not in MANUAL_ACCOUNTING_TERMINAL_STATES
                or not isinstance(record.get("timestamp"), str)
                or not record["timestamp"]
                or record.get("notes")
                != MANUAL_ACCOUNTING_SCOPE_NOTES[record["accounting_scope"]]
            ):
                raise ValueError("Manual-accounting ledger event semantics are invalid")
            allocated_nodes = record.get("scheduler_reported_allocated_nodes")
            billed_nodes = record.get("billed_nodes")
            elapsed_seconds = record.get("elapsed_seconds")
            if (
                type(allocated_nodes) is not int
                or allocated_nodes < 0
                or type(billed_nodes) is not int
                or billed_nodes != allocated_nodes
                or type(elapsed_seconds) is not int
                or elapsed_seconds < 0
            ):
                raise ValueError("Manual-accounting ledger event usage is invalid")
            consumed = _nonnegative_finite_number(
                record, "consumed_node_hours", label="Manual-accounting ledger event"
            )
            if consumed != billed_nodes * elapsed_seconds / 3600.0:
                raise ValueError("Manual-accounting ledger event usage differs")
            cumulative = _nonnegative_finite_number(
                record,
                "cumulative_consumed_node_hours",
                label="Manual-accounting ledger event",
            )
            expected_cumulative = (
                _accounting_without_validation(records[:index])[
                    "cumulative_consumed_node_hours"
                ]
                + consumed
            )
            if not math.isclose(
                cumulative, expected_cumulative, rel_tol=0.0, abs_tol=1.0e-12
            ):
                raise ValueError(
                    "Manual-accounting ledger event cumulative usage differs"
                )
        elif event_type == "reconciliation":
            if record.get("reconciled") is not True:
                raise ValueError("Registered reconciliation must be reconciled")
            if record.get("state") not in MANUAL_ACCOUNTING_TERMINAL_STATES:
                raise ValueError("Registered reconciliation state is not terminal")
            reservation_id = record.get("reservation_id")
            prior = latest_reservations(records[:index]).get(str(reservation_id))
            _require_canonical_scheduler_job_id(record.get("job_id"))
            if (
                prior is None
                or prior.get("event_type") != "job_id_attached"
                or prior.get("reconciled", False) is not False
                or prior.get("job_id") != record.get("job_id")
            ):
                raise ValueError("Registered reconciliation transition is invalid")
            _require_registered_transition_payload(
                prior, record, allowed_changes=REGISTERED_RECONCILIATION_CHANGES
            )
            _require_terminal_recovery_provenance(record, prior=prior)
            _require_purged_zero_execution_semantics(record)
            requested_nodes = record.get("requested_nodes")
            allocated_nodes = record.get("scheduler_reported_allocated_nodes")
            billed_nodes = record.get("billed_nodes")
            elapsed_seconds = record.get("elapsed_seconds")
            if (
                type(requested_nodes) is not int
                or requested_nodes <= 0
                or type(allocated_nodes) is not int
                or allocated_nodes < 0
                or type(billed_nodes) is not int
                or billed_nodes != max(requested_nodes, allocated_nodes)
                or type(elapsed_seconds) is not int
                or elapsed_seconds < 0
            ):
                raise ValueError("Registered reconciliation usage is invalid")
            consumed = _nonnegative_finite_number(
                record, "consumed_node_hours", label="Registered reconciliation"
            )
            if consumed != billed_nodes * elapsed_seconds / 3600.0:
                raise ValueError("Registered reconciliation usage differs")
            cumulative = _nonnegative_finite_number(
                record,
                "cumulative_consumed_node_hours",
                label="Registered reconciliation",
            )
            expected_cumulative = (
                _accounting_without_validation(records[:index])[
                    "cumulative_consumed_node_hours"
                ]
                + consumed
            )
            if not math.isclose(
                cumulative, expected_cumulative, rel_tol=0.0, abs_tol=1.0e-12
            ):
                raise ValueError("Registered reconciliation cumulative usage differs")


def accounting(records: list[dict[str, object]]) -> dict[str, float]:
    _validate_accounting_records(records)
    return _accounting_without_validation(records)


def latest_reservations(
    records: list[dict[str, object]],
) -> dict[str, dict[str, object]]:
    latest: dict[str, dict[str, object]] = {}
    for record in records:
        reservation_id = str(record.get("reservation_id", ""))
        if reservation_id:
            latest[reservation_id] = record
    return latest


def transition_payload(record: dict[str, object]) -> dict[str, object]:
    ignored = {
        "event_sha256",
        "event_type",
        "previous_event_sha256",
        "sequence_number",
        "timestamp",
    }
    return {key: value for key, value in record.items() if key not in ignored}


def _require_canonical_path(path: Path) -> Path:
    lexical = Path(os.path.abspath(path))
    resolved = path.resolve()
    if lexical != resolved:
        raise ValueError(f"Path must not use a symlink alias: {lexical}")
    return lexical


def _open_regular_nofollow(
    path: Path,
    mode: str,
    *,
    flags: int,
    newline: str | None = None,
    dir_fd: int | None = None,
) -> TextIO:
    descriptor = os.open(path, flags | os.O_NOFOLLOW, 0o600, dir_fd=dir_fd)
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise ValueError(f"Expected a regular file: {path}")
        return os.fdopen(descriptor, mode, encoding="utf-8", newline=newline)
    except BaseException:
        os.close(descriptor)
        raise


def _require_same_directory(path: Path, descriptor: int) -> None:
    lexical_descriptor = os.open(
        path,
        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
    )
    try:
        expected = os.fstat(descriptor)
        actual = os.fstat(lexical_descriptor)
        if (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino):
            raise ValueError(f"Writable parent path changed while locked: {path}")
    finally:
        os.close(lexical_descriptor)


def _require_same_regular_file_at(
    parent_descriptor: int, name: str, descriptor: int, *, label: str
) -> None:
    try:
        actual = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError as error:
        raise ValueError(f"{label} path changed while locked: {name}") from error
    expected = os.fstat(descriptor)
    if (
        not stat.S_ISREG(actual.st_mode)
        or (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino)
    ):
        raise ValueError(f"{label} path changed while locked: {name}")


@contextmanager
def _pinned_parent_directories(
    paths: list[Path],
    *,
    create_missing: bool = True,
    roots: dict[Path, Path] | None = None,
) -> Iterator[None]:
    existing = _PINNED_PARENT_DESCRIPTORS.get()
    opened: dict[Path, int] = {}
    try:
        for path in paths:
            parent = Path(os.path.abspath(path.parent))
            if roots is not None and parent not in roots:
                raise ValueError(f"Missing authorized root for pinned parent: {parent}")
            if parent in existing:
                if roots is not None:
                    rooted_descriptor = open_directory_below(parent, root=roots[parent])
                    try:
                        expected = os.fstat(existing[parent])
                        actual = os.fstat(rooted_descriptor)
                        if (actual.st_dev, actual.st_ino) != (
                            expected.st_dev,
                            expected.st_ino,
                        ):
                            raise ValueError(
                                f"Pinned parent differs from authorized-root traversal: {parent}"
                            )
                    finally:
                        os.close(rooted_descriptor)
                continue
            if parent in opened:
                continue
            if create_missing:
                durable_mkdir_parents(parent)
            descriptor = (
                open_directory_below(parent, root=roots[parent])
                if roots is not None and parent in roots
                else os.open(
                    parent,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                )
            )
            try:
                _require_same_directory(parent, descriptor)
            except BaseException:
                os.close(descriptor)
                raise
            opened[parent] = descriptor
    except BaseException:
        for descriptor in opened.values():
            os.close(descriptor)
        raise
    token = _PINNED_PARENT_DESCRIPTORS.set({**existing, **opened})
    try:
        yield
        for parent, descriptor in opened.items():
            _require_same_directory(parent, descriptor)
    finally:
        _PINNED_PARENT_DESCRIPTORS.reset(token)
        for descriptor in opened.values():
            os.close(descriptor)


@contextmanager
def _local_ledger_lock(ledger_jsonl: Path) -> Iterator[None]:
    ledger_parent = Path(os.path.abspath(ledger_jsonl.parent))
    durable_mkdir_parents(ledger_parent)
    lock_path = _require_canonical_path(
        ledger_jsonl.with_suffix(ledger_jsonl.suffix + ".lock")
    )
    stable_parent = _require_canonical_path(ledger_parent.parent)
    anchored_lock_path = _require_canonical_path(stable_parent / ".ledger.lock")
    stable_parent_descriptor = os.open(
        stable_parent,
        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
    )
    try:
        fcntl.flock(stable_parent_descriptor, fcntl.LOCK_EX)
        try:
            _require_same_directory(stable_parent, stable_parent_descriptor)
            ledger_parent_descriptor = os.open(
                ledger_parent.name,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=stable_parent_descriptor,
            )
            try:
                fcntl.flock(ledger_parent_descriptor, fcntl.LOCK_EX)
                try:
                    _require_same_directory(stable_parent, stable_parent_descriptor)
                    _require_same_directory(ledger_parent, ledger_parent_descriptor)
                    with _open_regular_nofollow(
                        Path(anchored_lock_path.name),
                        "a",
                        flags=os.O_APPEND | os.O_CREAT | os.O_WRONLY,
                        dir_fd=stable_parent_descriptor,
                    ) as anchored_lock_stream:
                        fcntl.flock(anchored_lock_stream.fileno(), fcntl.LOCK_EX)
                        try:
                            _require_same_regular_file_at(
                                stable_parent_descriptor,
                                anchored_lock_path.name,
                                anchored_lock_stream.fileno(),
                                label="Anchored ledger lock",
                            )
                            # Retain the historical lock for deployed callers while
                            # pinned directory locks prevent a replacement split.
                            with _open_regular_nofollow(
                                Path(lock_path.name),
                                "a",
                                flags=os.O_APPEND | os.O_CREAT | os.O_WRONLY,
                                dir_fd=ledger_parent_descriptor,
                            ) as lock_stream:
                                fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX)
                                try:
                                    _require_same_directory(
                                        stable_parent, stable_parent_descriptor
                                    )
                                    _require_same_directory(
                                        ledger_parent, ledger_parent_descriptor
                                    )
                                    _require_same_regular_file_at(
                                        stable_parent_descriptor,
                                        anchored_lock_path.name,
                                        anchored_lock_stream.fileno(),
                                        label="Anchored ledger lock",
                                    )
                                    _require_same_regular_file_at(
                                        ledger_parent_descriptor,
                                        lock_path.name,
                                        lock_stream.fileno(),
                                        label="Historical ledger lock",
                                    )
                                    yield
                                finally:
                                    try:
                                        _require_same_directory(
                                            stable_parent, stable_parent_descriptor
                                        )
                                        _require_same_directory(
                                            ledger_parent, ledger_parent_descriptor
                                        )
                                        _require_same_regular_file_at(
                                            stable_parent_descriptor,
                                            anchored_lock_path.name,
                                            anchored_lock_stream.fileno(),
                                            label="Anchored ledger lock",
                                        )
                                        _require_same_regular_file_at(
                                            ledger_parent_descriptor,
                                            lock_path.name,
                                            lock_stream.fileno(),
                                            label="Historical ledger lock",
                                        )
                                    finally:
                                        fcntl.flock(
                                            lock_stream.fileno(), fcntl.LOCK_UN
                                        )
                        finally:
                            _require_same_directory(
                                stable_parent, stable_parent_descriptor
                            )
                            try:
                                _require_same_regular_file_at(
                                    stable_parent_descriptor,
                                    anchored_lock_path.name,
                                    anchored_lock_stream.fileno(),
                                    label="Anchored ledger lock",
                                )
                            finally:
                                fcntl.flock(
                                    anchored_lock_stream.fileno(), fcntl.LOCK_UN
                                )
                finally:
                    try:
                        _require_same_directory(
                            stable_parent, stable_parent_descriptor
                        )
                        _require_same_directory(
                            ledger_parent, ledger_parent_descriptor
                        )
                    finally:
                        fcntl.flock(ledger_parent_descriptor, fcntl.LOCK_UN)
            finally:
                os.close(ledger_parent_descriptor)
        finally:
            try:
                _require_same_directory(stable_parent, stable_parent_descriptor)
            finally:
                fcntl.flock(stable_parent_descriptor, fcntl.LOCK_UN)
    finally:
        os.close(stable_parent_descriptor)


@contextmanager
def _ledger_lock_within_serialization_anchor(
    ledger_jsonl: Path, mirror_jsonl: Path | None = None
) -> Iterator[None]:
    if mirror_jsonl is None:
        with _local_ledger_lock(ledger_jsonl):
            yield
        return

    mirror_parent = Path(os.path.abspath(mirror_jsonl.parent))
    durable_mkdir_parents(mirror_parent)
    mirror_parent_descriptor = os.open(
        mirror_parent,
        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
    )
    try:
        fcntl.flock(mirror_parent_descriptor, fcntl.LOCK_EX)
        try:
            _require_same_directory(mirror_parent, mirror_parent_descriptor)
            with _local_ledger_lock(ledger_jsonl):
                _require_same_directory(mirror_parent, mirror_parent_descriptor)
                yield
                _require_same_directory(mirror_parent, mirror_parent_descriptor)
        finally:
            try:
                _require_same_directory(mirror_parent, mirror_parent_descriptor)
            finally:
                fcntl.flock(mirror_parent_descriptor, fcntl.LOCK_UN)
    finally:
        os.close(mirror_parent_descriptor)


def incomplete_manual_accounting_marker_paths(
    ledger_jsonl: Path, mirror_jsonl: Path | None
) -> tuple[Path, Path | None]:
    return (
        ledger_jsonl.parent / INCOMPLETE_MANUAL_ACCOUNTING_MARKER_FILENAME,
        (
            None
            if mirror_jsonl is None
            else mirror_jsonl.parent / INCOMPLETE_MANUAL_ACCOUNTING_MARKER_FILENAME
        ),
    )


def _incomplete_manual_accounting_marker_bytes(
    marker: dict[str, object],
) -> bytes:
    return (json.dumps(marker, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
        "utf-8"
    )


def _validate_incomplete_manual_accounting_marker(
    marker: dict[str, object],
) -> None:
    if (
        set(marker)
        != {
            "schema_version",
            "state",
            "manual_accounting_authorization_id",
            "manual_accounting_authorization_sha256",
            "pre_tranche_sequence_number",
            "pre_tranche_chain_head",
            "pre_tranche_authorized_job_count",
            "control_plane_version",
            "active_policy_sha256",
            "active_promotion_sha256",
        }
        or type(marker.get("schema_version")) is not int
        or marker.get("schema_version") != 3
        or marker.get("state") != "manual_accounting_incomplete"
        or not isinstance(marker.get("manual_accounting_authorization_id"), str)
        or not isinstance(marker.get("manual_accounting_authorization_sha256"), str)
        or type(marker.get("pre_tranche_sequence_number")) is not int
        or not isinstance(marker.get("pre_tranche_chain_head"), str)
        or type(marker.get("pre_tranche_authorized_job_count")) is not int
    ):
        raise ValueError("Incomplete manual-accounting marker is malformed")
    authorization_id = str(marker["manual_accounting_authorization_id"])
    authorization_sha256 = str(marker["manual_accounting_authorization_sha256"])
    pre_tranche_sequence_number = int(marker["pre_tranche_sequence_number"])
    pre_tranche_chain_head = str(marker["pre_tranche_chain_head"])
    pre_tranche_authorized_job_count = int(marker["pre_tranche_authorized_job_count"])
    if (
        not authorization_id
        or len(authorization_id) > 128
        or not authorization_sha256
        or len(authorization_sha256) != 64
        or any(character not in "0123456789abcdef" for character in authorization_sha256)
        or pre_tranche_sequence_number < 0
        or pre_tranche_authorized_job_count < 0
        or not _is_lowercase_sha256(marker.get("control_plane_version"))
        or not _is_lowercase_sha256(marker.get("active_policy_sha256"))
        or not _is_lowercase_sha256(marker.get("active_promotion_sha256"))
        or (
            pre_tranche_chain_head != ""
            and (
                len(pre_tranche_chain_head) != 64
                or any(
                    character not in "0123456789abcdef"
                    for character in pre_tranche_chain_head
                )
            )
        )
        or (pre_tranche_sequence_number == 0) != (pre_tranche_chain_head == "")
    ):
        raise ValueError("Incomplete manual-accounting marker binding is malformed")


def _read_incomplete_manual_accounting_marker(path: Path) -> dict[str, object]:
    marker = read_json_bytes(
        _read_regular_bytes(path, require_read_only_mode=True),
        label=str(path),
    )
    _validate_incomplete_manual_accounting_marker(marker)
    return marker


def require_no_incomplete_manual_accounting_marker(
    ledger_jsonl: Path, mirror_jsonl: Path | None
) -> None:
    local_marker, mirror_marker = incomplete_manual_accounting_marker_paths(
        ledger_jsonl, mirror_jsonl
    )
    paths = [local_marker] + ([] if mirror_marker is None else [mirror_marker])
    with _pinned_parent_directories(paths):
        if any(_path_exists(path) for path in paths):
            raise ValueError("Ledger mutation is blocked by incomplete manual accounting")


def publish_incomplete_manual_accounting_marker_locked(
    ledger_jsonl: Path,
    mirror_jsonl: Path,
    marker: dict[str, object],
) -> None:
    _validate_incomplete_manual_accounting_marker(marker)
    expected = _incomplete_manual_accounting_marker_bytes(marker)
    local_marker, mirror_marker = incomplete_manual_accounting_marker_paths(
        ledger_jsonl, mirror_jsonl
    )
    assert mirror_marker is not None
    with _pinned_parent_directories([local_marker, mirror_marker]):
        for path in [local_marker, mirror_marker]:
            if _path_exists(path) and _read_regular_bytes(
                path, require_read_only_mode=True
            ) != expected:
                raise ValueError("Incomplete manual-accounting marker differs")
        for path in [local_marker, mirror_marker]:
            if not _path_exists(path):
                parent_descriptor = _parent_descriptor(path)
                assert parent_descriptor is not None
                atomic_write_bytes_at(
                    parent_descriptor,
                    path.name,
                    expected,
                    mode=0o444,
                    replace=False,
                )
        if any(
            _read_regular_bytes(path, require_read_only_mode=True) != expected
            for path in [local_marker, mirror_marker]
        ):
            raise ValueError("Incomplete manual-accounting marker publication differs")


def clear_matching_incomplete_manual_accounting_marker_locked(
    ledger_jsonl: Path,
    mirror_jsonl: Path,
    marker: dict[str, object],
) -> None:
    _validate_incomplete_manual_accounting_marker(marker)
    expected = _incomplete_manual_accounting_marker_bytes(marker)
    local_marker, mirror_marker = incomplete_manual_accounting_marker_paths(
        ledger_jsonl, mirror_jsonl
    )
    assert mirror_marker is not None
    with _pinned_parent_directories([local_marker, mirror_marker]):
        paths = [local_marker, mirror_marker]
        existing = [path for path in paths if _path_exists(path)]
        if not existing:
            return
        if any(
            _read_regular_bytes(path, require_read_only_mode=True) != expected
            for path in existing
        ):
            raise ValueError("Incomplete manual-accounting marker pair differs")
        for path in existing:
            parent_descriptor = _parent_descriptor(path)
            assert parent_descriptor is not None
            os.unlink(path.name, dir_fd=parent_descriptor)
            os.fsync(parent_descriptor)


def _validate_recoverable_manual_accounting_suffix(
    records: list[dict[str, object]],
    *,
    pre_tranche_sequence_number: int,
    authorization_id: str,
    authorization_sha256: str,
    authorization_path: Path,
    project_home_authorization_path: Path,
    reviewed_job_ids: list[str],
    reviewed_scheduler_results: dict[str, dict[str, object]],
    pre_tranche_authorized_job_count: int,
    control_plane_version: str,
    active_policy_sha256: str,
    active_promotion_sha256: str,
    accounting_scope: str,
) -> None:
    suffix = records[pre_tranche_sequence_number:]
    authorized_prefix = [
        record
        for record in records[:pre_tranche_sequence_number]
        if (
            record.get("event_type") == "manual_allocation_reconciliation"
            and record.get("manual_accounting_authorization_id") == authorization_id
        )
    ]
    if (
        pre_tranche_authorized_job_count > len(reviewed_job_ids)
        or [record.get("job_id") for record in authorized_prefix]
        != reviewed_job_ids[:pre_tranche_authorized_job_count]
    ):
        raise ValueError("Incomplete manual-accounting marker authorized prefix differs")
    suffix_job_ids = reviewed_job_ids[pre_tranche_authorized_job_count:]
    if len(suffix) > len(suffix_job_ids):
        raise ValueError("Incomplete manual-accounting Orion suffix exceeds authorization")
    cumulative = accounting(records[:pre_tranche_sequence_number])[
        "cumulative_consumed_node_hours"
    ]
    for record, job_id in zip(suffix, suffix_job_ids):
        scheduler = reviewed_scheduler_results[job_id]
        _require_canonical_scheduler_job_id(job_id)
        _require_new_scheduler_job_id(
            records[:pre_tranche_sequence_number], job_id
        )
        if set(record) != MANUAL_ACCOUNTING_EVENT_FIELDS:
            raise ValueError("Incomplete manual-accounting Orion suffix event is malformed")
        if (
            record.get("event_type") != "manual_allocation_reconciliation"
            or record.get("job_id") != job_id
            or record.get("manual_accounting_authorization_id") != authorization_id
            or record.get("manual_accounting_authorization_sha256")
            != authorization_sha256
            or record.get("manual_accounting_authorization_path")
            != str(authorization_path)
            or record.get("manual_accounting_project_home_authorization_path")
            != str(project_home_authorization_path)
            or record.get("accounting_scope") != accounting_scope
            or record.get("scientific_evidence_eligible") is not False
            or record.get("reconciled") is not True
            or record.get("control_plane_version") != control_plane_version
            or record.get("reconciled_by_control_plane_version")
            != control_plane_version
            or record.get("active_policy_sha256") != active_policy_sha256
            or record.get("active_promotion_sha256") != active_promotion_sha256
            or record.get("partition") != scheduler["partition"]
            or record.get("qos") != scheduler["qos"]
            or record.get("state") != scheduler["state"]
            or record.get("scheduler_reported_allocated_nodes")
            != scheduler["allocated_nodes"]
            or record.get("elapsed_seconds") != scheduler["elapsed_seconds"]
            or not isinstance(record.get("timestamp"), str)
            or not record["timestamp"]
            or record.get("notes")
            != MANUAL_ACCOUNTING_SCOPE_NOTES[accounting_scope]
        ):
            raise ValueError("Incomplete manual-accounting Orion suffix event is invalid")
        allocated_nodes = record.get("scheduler_reported_allocated_nodes")
        billed_nodes = record.get("billed_nodes")
        elapsed_seconds = record.get("elapsed_seconds")
        if (
            type(allocated_nodes) is not int
            or allocated_nodes < 0
            or type(billed_nodes) is not int
            or billed_nodes != allocated_nodes
            or type(elapsed_seconds) is not int
            or elapsed_seconds < 0
        ):
            raise ValueError("Incomplete manual-accounting Orion suffix usage is invalid")
        consumed = _nonnegative_finite_number(
            record, "consumed_node_hours", label="Manual-accounting Orion suffix"
        )
        if consumed != billed_nodes * elapsed_seconds / 3600.0:
            raise ValueError("Incomplete manual-accounting Orion suffix usage differs")
        cumulative += consumed
        if (
            _nonnegative_finite_number(
                record,
                "cumulative_consumed_node_hours",
                label="Manual-accounting Orion suffix",
            )
            != cumulative
        ):
            raise ValueError("Incomplete manual-accounting Orion suffix cumulative usage differs")


def recover_incomplete_manual_accounting_locked(
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    authorization_id: str,
    authorization_sha256: str,
    authorization_path: Path,
    project_home_authorization_path: Path,
    reviewed_job_ids: list[str],
    reviewed_scheduler_results: dict[str, dict[str, object]],
    control_plane_version: str,
    active_policy_sha256: str,
    active_promotion_sha256: str,
    accounting_scope: str,
) -> dict[str, object] | None:
    """Repair only a marker-bound manual-accounting publication from Orion."""
    local_marker, mirror_marker = incomplete_manual_accounting_marker_paths(
        ledger_jsonl, mirror_jsonl
    )
    assert mirror_marker is not None
    paths = [
        *_ledger_state_paths(ledger_jsonl, None, receipts_jsonl, mirror_jsonl),
        local_marker,
        mirror_marker,
    ]
    with _pinned_parent_directories(paths):
        markers = [
            _read_incomplete_manual_accounting_marker(path)
            for path in [local_marker, mirror_marker]
            if _path_exists(path)
        ]
        if not markers:
            return None
        marker = markers[0]
        if any(candidate != marker for candidate in markers[1:]):
            raise ValueError("Incomplete manual-accounting marker pair differs")
        if (
            marker["manual_accounting_authorization_id"] != authorization_id
            or marker["manual_accounting_authorization_sha256"]
            != authorization_sha256
            or marker["control_plane_version"] != control_plane_version
            or marker["active_policy_sha256"] != active_policy_sha256
            or marker["active_promotion_sha256"] != active_promotion_sha256
        ):
            raise ValueError("Incomplete manual-accounting marker binding differs")
        local_records = validate_primary_chain(ledger_jsonl)
        pre_tranche_sequence_number = int(marker["pre_tranche_sequence_number"])
        if (
            pre_tranche_sequence_number > len(local_records)
            or chain_head(local_records[:pre_tranche_sequence_number])
            != marker["pre_tranche_chain_head"]
        ):
            raise ValueError("Incomplete manual-accounting marker Orion boundary differs")
        require_explicit_genesis(local_records)
        _validate_recoverable_manual_accounting_suffix(
            local_records,
            pre_tranche_sequence_number=pre_tranche_sequence_number,
            authorization_id=authorization_id,
            authorization_sha256=authorization_sha256,
            authorization_path=authorization_path,
            project_home_authorization_path=project_home_authorization_path,
            reviewed_job_ids=reviewed_job_ids,
            reviewed_scheduler_results=reviewed_scheduler_results,
            pre_tranche_authorized_job_count=int(
                marker["pre_tranche_authorized_job_count"]
            ),
            control_plane_version=control_plane_version,
            active_policy_sha256=active_policy_sha256,
            active_promotion_sha256=active_promotion_sha256,
            accounting_scope=accounting_scope,
        )
        mirror_records = validate_primary_chain(mirror_jsonl)
        if mirror_records != local_records[:len(mirror_records)]:
            raise ValueError("Project Home ledger is not an exact prefix of canonical Orion")
        receipts = _read_jsonl(receipts_jsonl)
        for receipt in receipts:
            _validate_receipt_provenance(
                receipt,
                mirror_jsonl=mirror_jsonl,
                mirror_transport="filesystem_copy",
            )
            if receipt.get("mirror_ack_sha256") != record_sha256(
                receipt, "mirror_ack_sha256"
            ):
                raise ValueError("Cannot recover a corrupt mirror receipt")
        if (
            len(mirror_records) < pre_tranche_sequence_number
            or len(receipts) < pre_tranche_sequence_number
        ):
            raise ValueError(
                "Incomplete manual-accounting recovery cannot repair pre-tranche publication loss"
            )
        if [receipt.get("mirrored_event_sha256") for receipt in receipts] != [
            record["event_sha256"] for record in mirror_records[:len(receipts)]
        ]:
            raise ValueError("Mirror receipts are not an exact prefix of mirrored events")
        if not receipts:
            raise ValueError("Cannot recover manual accounting without its genesis receipt")
        _validate_genesis_anchors(ledger_jsonl, mirror_jsonl, local_records, receipts)
        _validate_operator_attestation_trees(
            local_records, ledger_jsonl=ledger_jsonl, mirror_jsonl=mirror_jsonl
        )
        publish_incomplete_manual_accounting_marker_locked(
            ledger_jsonl, mirror_jsonl, marker
        )
        mirror_preflight(mirror_jsonl)
        for record in local_records[len(mirror_records):]:
            _append_jsonl(mirror_jsonl, record)
            mirror_records.append(record)
        for record in mirror_records[len(receipts):]:
            _append_mirror_receipt(
                receipts_jsonl, record, mirror_jsonl, "filesystem_copy"
            )
        validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
        return marker


@contextmanager
def ledger_lock(
    ledger_jsonl: Path,
    mirror_jsonl: Path | None = None,
    *,
    allow_incomplete_manual_accounting: bool = False,
) -> Iterator[None]:
    pic_root = Path(os.path.abspath(ledger_jsonl.parent.parent))
    anchor = stable_serialization_anchor(pic_root)
    durable_mkdir_parents(anchor)
    anchor_descriptor = os.open(
        anchor,
        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
    )
    try:
        fcntl.flock(anchor_descriptor, fcntl.LOCK_EX)
        try:
            _require_same_directory(anchor, anchor_descriptor)
            with _ledger_lock_within_serialization_anchor(ledger_jsonl, mirror_jsonl):
                _require_same_directory(anchor, anchor_descriptor)
                if not allow_incomplete_manual_accounting:
                    require_no_incomplete_manual_accounting_marker(
                        ledger_jsonl, mirror_jsonl
                    )
                yield
                _require_same_directory(anchor, anchor_descriptor)
        finally:
            try:
                _require_same_directory(anchor, anchor_descriptor)
            finally:
                fcntl.flock(anchor_descriptor, fcntl.LOCK_UN)
    finally:
        os.close(anchor_descriptor)


def mirror_preflight(mirror_jsonl: Path) -> None:
    parent_descriptor = _parent_descriptor(mirror_jsonl)
    if parent_descriptor is None:
        durable_mkdir_parents(mirror_jsonl.parent)
    probe_name = f".pic-ledger-probe-{os.getpid()}-{uuid.uuid4()}"
    descriptor = os.open(
        probe_name if parent_descriptor is not None else mirror_jsonl.parent / probe_name,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
        0o600,
        dir_fd=parent_descriptor,
    )
    try:
        os.write(descriptor, b"pic-ledger-mirror-preflight\n")
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    if parent_descriptor is not None:
        os.unlink(probe_name, dir_fd=parent_descriptor)
        os.fsync(parent_descriptor)
    else:
        (mirror_jsonl.parent / probe_name).unlink()
        fsync_directory(mirror_jsonl.parent)


def _append_mirror_receipt(
    receipts_jsonl: Path,
    event: dict[str, object],
    mirror_jsonl: Path,
    mirror_transport: str,
) -> dict[str, object]:
    receipt: dict[str, object] = {
        "mirrored_event_sha256": event["event_sha256"],
        "mirror_destination": str(mirror_jsonl),
        "mirror_transport": mirror_transport,
        "mirror_acknowledged_utc": utc_now(),
    }
    receipt["mirror_ack_sha256"] = record_sha256(receipt, "mirror_ack_sha256")
    _append_jsonl(receipts_jsonl, receipt)
    return receipt


def genesis_anchor_paths(ledger_jsonl: Path, mirror_jsonl: Path) -> tuple[Path, Path]:
    return (
        ledger_jsonl.parent / GENESIS_ANCHOR_FILENAME,
        mirror_jsonl.parent / GENESIS_ANCHOR_FILENAME,
    )


def _genesis_anchor(
    first_event: dict[str, object],
    first_receipt: dict[str, object],
    *,
    mirror_jsonl: Path,
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "event_sha256": first_event["event_sha256"],
        "mirror_ack_sha256": first_receipt["mirror_ack_sha256"],
        "mirror_destination": str(mirror_jsonl),
        "mirror_transport": "filesystem_copy",
    }


def _write_genesis_anchors(
    ledger_jsonl: Path,
    mirror_jsonl: Path,
    anchor: dict[str, object],
) -> None:
    local_anchor, mirror_anchor = genesis_anchor_paths(ledger_jsonl, mirror_jsonl)
    _write_one_genesis_anchor(local_anchor, ledger_jsonl.parent.parent, anchor)
    _write_one_genesis_anchor(mirror_anchor, mirror_jsonl.parent.parent, anchor)


def _write_one_genesis_anchor(
    path: Path,
    root: Path,
    anchor: dict[str, object],
) -> None:
    parent_descriptor = _parent_descriptor(path)
    if parent_descriptor is None:
        atomic_write_json(path, anchor, replace=False, root=root)
    else:
        atomic_write_json_at(parent_descriptor, path.name, anchor, replace=False)


def _read_one_genesis_anchor_bytes(path: Path, root: Path) -> bytes:
    return (
        _read_regular_bytes(path, require_read_only_mode=True)
        if _parent_descriptor(path) is not None
        else read_stable_regular_file_below(
            path,
            root,
            require_read_only_mode=True,
        )
    )


def _validate_genesis_anchor_schema(anchor: dict[str, object]) -> None:
    if (
        type(anchor.get("schema_version")) is not int
        or anchor.get("schema_version") != 1
    ):
        raise ValueError("Unsupported Frontier PIC genesis-anchor schema")


def _read_one_genesis_anchor(path: Path, root: Path) -> dict[str, object]:
    anchor = read_json_bytes(_read_one_genesis_anchor_bytes(path, root), label=str(path))
    _validate_genesis_anchor_schema(anchor)
    return anchor


def _validate_genesis_anchors(
    ledger_jsonl: Path,
    mirror_jsonl: Path,
    records: list[dict[str, object]],
    receipts: list[dict[str, object]],
) -> None:
    local_anchor, mirror_anchor = genesis_anchor_paths(ledger_jsonl, mirror_jsonl)
    if not records:
        if _path_exists(local_anchor) or _path_exists(mirror_anchor):
            raise ValueError("Genesis anchor exists without Frontier PIC ledger records")
        return
    require_explicit_genesis(records)
    try:
        local_bytes = _read_one_genesis_anchor_bytes(
            local_anchor, ledger_jsonl.parent.parent
        )
        mirror_bytes = _read_one_genesis_anchor_bytes(
            mirror_anchor, mirror_jsonl.parent.parent
        )
    except FileNotFoundError as error:
        raise ValueError("Frontier PIC ledger is missing a genesis anchor") from error
    _validate_genesis_anchor_bytes(
        local_bytes,
        mirror_bytes,
        ledger_jsonl=ledger_jsonl,
        mirror_jsonl=mirror_jsonl,
        records=records,
        receipts=receipts,
    )


def _validate_genesis_anchor_bytes(
    local_bytes: bytes,
    mirror_bytes: bytes,
    *,
    ledger_jsonl: Path,
    mirror_jsonl: Path,
    records: list[dict[str, object]],
    receipts: list[dict[str, object]],
) -> None:
    require_explicit_genesis(records)
    if local_bytes != mirror_bytes:
        raise ValueError("Orion and Project Home genesis-anchor bytes differ")
    local_anchor, _ = genesis_anchor_paths(ledger_jsonl, mirror_jsonl)
    anchor = read_json_bytes(local_bytes, label=str(local_anchor))
    _validate_genesis_anchor_schema(anchor)
    if anchor != _genesis_anchor(records[0], receipts[0], mirror_jsonl=mirror_jsonl):
        raise ValueError("Frontier PIC genesis anchor differs from ledger genesis")


def validate_mirrored_state(
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    ledger_root: Path | None = None,
    receipts_root: Path | None = None,
    mirror_root: Path | None = None,
) -> list[dict[str, object]]:
    local_records = validate_primary_chain(ledger_jsonl, root=ledger_root)
    mirror_records = validate_primary_chain(mirror_jsonl, root=mirror_root)
    if local_records != mirror_records:
        raise ValueError("Local and mirrored PIC ledger records differ")
    _validate_operator_attestation_trees(
        local_records, ledger_jsonl=ledger_jsonl, mirror_jsonl=mirror_jsonl
    )
    receipts = validate_receipts(
        receipts_jsonl,
        local_records,
        mirror_jsonl=mirror_jsonl,
        mirror_transport="filesystem_copy",
        root=receipts_root,
    )
    _validate_genesis_anchors(ledger_jsonl, mirror_jsonl, local_records, receipts)
    _validate_terminal_recovery_handoff_mirrors(
        local_records, ledger_jsonl=ledger_jsonl, mirror_jsonl=mirror_jsonl
    )
    return local_records


def _validate_mirrored_state_bytes(
    state: dict[Path, bytes],
    *,
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    validate_handoffs: bool = True,
) -> list[dict[str, object]]:
    local_records = _validate_primary_records(
        _read_jsonl_bytes(state[ledger_jsonl], path=ledger_jsonl),
        path=ledger_jsonl,
    )
    mirror_records = _validate_primary_records(
        _read_jsonl_bytes(state[mirror_jsonl], path=mirror_jsonl),
        path=mirror_jsonl,
    )
    if local_records != mirror_records:
        raise ValueError("Local and mirrored PIC ledger records differ")
    _validate_operator_attestation_trees(
        local_records, ledger_jsonl=ledger_jsonl, mirror_jsonl=mirror_jsonl
    )
    receipts = _validate_receipt_records(
        _read_jsonl_bytes(state[receipts_jsonl], path=receipts_jsonl),
        local_records,
        path=receipts_jsonl,
        mirror_jsonl=mirror_jsonl,
        mirror_transport="filesystem_copy",
    )
    local_anchor, mirror_anchor = genesis_anchor_paths(ledger_jsonl, mirror_jsonl)
    _validate_genesis_anchor_bytes(
        state[local_anchor],
        state[mirror_anchor],
        ledger_jsonl=ledger_jsonl,
        mirror_jsonl=mirror_jsonl,
        records=local_records,
        receipts=receipts,
    )
    if validate_handoffs:
        _validate_terminal_recovery_handoff_mirrors(
            local_records,
            ledger_jsonl=ledger_jsonl,
            mirror_jsonl=mirror_jsonl,
            state=state,
        )
    return local_records


def _terminal_recovery_handoff_paths(
    records: list[dict[str, object]],
    *,
    ledger_jsonl: Path,
    mirror_jsonl: Path,
) -> list[Path]:
    pic_root = Path(os.path.abspath(ledger_jsonl.parent.parent))
    project_home_root = Path(os.path.abspath(mirror_jsonl.parent.parent))
    paths: list[Path] = []
    for record in records:
        if not TERMINAL_RECOVERY_FIELDS & set(record):
            continue
        handoff_path = Path(str(record["terminal_recovery_handoff_path"]))
        if handoff_path.parent != pic_root / "policy" / "recovery_handoffs":
            raise ValueError("Terminal-recovery handoff is outside the PIC root")
        mirror_path = project_home_root / "policy" / "recovery_handoffs" / handoff_path.name
        for path in [handoff_path, mirror_path]:
            if path not in paths:
                paths.append(path)
    return paths


def _ledger_state_paths(
    ledger_jsonl: Path,
    ledger_csv: Path | None,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
) -> list[Path]:
    paths = [ledger_jsonl, receipts_jsonl, mirror_jsonl]
    if ledger_csv is not None:
        paths.append(ledger_csv)
    local_anchor, mirror_anchor = genesis_anchor_paths(ledger_jsonl, mirror_jsonl)
    return [*paths, local_anchor, mirror_anchor]


def _read_descriptor_bytes(descriptor: int) -> bytes:
    chunks = []
    offset = 0
    while True:
        chunk = os.pread(descriptor, 1024 * 1024, offset)
        if not chunk:
            return b"".join(chunks)
        chunks.append(chunk)
        offset += len(chunk)


def _require_same_pinned_regular_files(descriptors: dict[Path, int]) -> None:
    for path, descriptor in descriptors.items():
        parent_descriptor = _parent_descriptor(path)
        if parent_descriptor is None:
            raise ValueError(f"Missing pinned parent descriptor for snapshot file: {path}")
        _require_same_regular_file_at(
            parent_descriptor,
            path.name,
            descriptor,
            label="Read-only ledger snapshot",
        )


@contextmanager
def validated_read_only_mirrored_state_snapshot(
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    ledger_root: Path,
    receipts_root: Path,
    mirror_root: Path,
) -> Iterator[list[dict[str, object]]]:
    """Pin and recheck exact authoritative bytes without taking writer locks."""
    paths = _ledger_state_paths(ledger_jsonl, None, receipts_jsonl, mirror_jsonl)
    read_only_paths = set(genesis_anchor_paths(ledger_jsonl, mirror_jsonl))
    roots: dict[Path, Path] = {}
    for path, root in [
        (ledger_jsonl, ledger_root),
        (receipts_jsonl, receipts_root),
        (mirror_jsonl, mirror_root),
    ]:
        parent = Path(os.path.abspath(path.parent))
        root = Path(os.path.abspath(root))
        if parent in roots and roots[parent] != root:
            raise ValueError(f"Snapshot parent has conflicting authorized roots: {parent}")
        roots[parent] = root
    with _pinned_parent_directories(paths, create_missing=False, roots=roots):
        descriptors: dict[Path, int] = {}
        try:
            for path in paths:
                parent_descriptor = _parent_descriptor(path)
                if parent_descriptor is None:
                    raise ValueError(
                        f"Missing pinned parent descriptor for snapshot file: {path}"
                    )
                descriptor = os.open(
                    path.name,
                    os.O_RDONLY | os.O_NOFOLLOW,
                    dir_fd=parent_descriptor,
                )
                metadata = os.fstat(descriptor)
                if not stat.S_ISREG(metadata.st_mode):
                    os.close(descriptor)
                    raise ValueError(f"Expected a regular snapshot file: {path}")
                if path in read_only_paths and metadata.st_mode & 0o222:
                    os.close(descriptor)
                    raise ValueError(f"Expected a read-only snapshot file: {path}")
                descriptors[path] = descriptor
            _require_same_pinned_regular_files(descriptors)
            initial = {
                path: _read_descriptor_bytes(descriptor)
                for path, descriptor in descriptors.items()
            }
            records = _validate_mirrored_state_bytes(
                initial,
                ledger_jsonl=ledger_jsonl,
                receipts_jsonl=receipts_jsonl,
                mirror_jsonl=mirror_jsonl,
                validate_handoffs=False,
            )
            handoff_paths = _terminal_recovery_handoff_paths(
                records, ledger_jsonl=ledger_jsonl, mirror_jsonl=mirror_jsonl
            )
            handoff_roots = {
                Path(os.path.abspath(path.parent)): (
                    Path(os.path.abspath(ledger_root))
                    if path.parent == ledger_jsonl.parent.parent / "policy" / "recovery_handoffs"
                    else Path(os.path.abspath(mirror_root))
                )
                for path in handoff_paths
            }
            with _pinned_parent_directories(
                handoff_paths, create_missing=False, roots=handoff_roots
            ):
                for path in handoff_paths:
                    parent_descriptor = _parent_descriptor(path)
                    if parent_descriptor is None:
                        raise ValueError(
                            f"Missing pinned parent descriptor for snapshot file: {path}"
                        )
                    descriptor = os.open(
                        path.name,
                        os.O_RDONLY | os.O_NOFOLLOW,
                        dir_fd=parent_descriptor,
                    )
                    metadata = os.fstat(descriptor)
                    if not stat.S_ISREG(metadata.st_mode):
                        os.close(descriptor)
                        raise ValueError(f"Expected a regular snapshot file: {path}")
                    if metadata.st_mode & 0o222:
                        os.close(descriptor)
                        raise ValueError(f"Expected a read-only snapshot file: {path}")
                    descriptors[path] = descriptor
                _require_same_pinned_regular_files(descriptors)
                validated = {
                    path: _read_descriptor_bytes(descriptor)
                    for path, descriptor in descriptors.items()
                }
                records = _validate_mirrored_state_bytes(
                    validated,
                    ledger_jsonl=ledger_jsonl,
                    receipts_jsonl=receipts_jsonl,
                    mirror_jsonl=mirror_jsonl,
                )
                if initial != {path: validated[path] for path in initial}:
                    raise ValueError("Authoritative PIC ledger bytes changed during validation")
                yield records
                _require_same_pinned_regular_files(descriptors)
                final = {
                    path: _read_descriptor_bytes(descriptor)
                    for path, descriptor in descriptors.items()
                }
                if validated != final:
                    raise ValueError("Authoritative PIC ledger bytes changed during snapshot use")
        finally:
            for descriptor in descriptors.values():
                os.close(descriptor)


def migrate_existing_genesis_anchors(
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    expected_event_sha256: str,
    expected_mirror_ack_sha256: str,
) -> dict[str, object]:
    """Create paired anchors for one audited pre-anchor ledger exactly once."""
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        with _pinned_parent_directories(
            _ledger_state_paths(ledger_jsonl, None, receipts_jsonl, mirror_jsonl)
        ):
            return _migrate_existing_genesis_anchors_locked(
                ledger_jsonl,
                receipts_jsonl,
                mirror_jsonl,
                expected_event_sha256=expected_event_sha256,
                expected_mirror_ack_sha256=expected_mirror_ack_sha256,
            )


def _migrate_existing_genesis_anchors_locked(
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    expected_event_sha256: str,
    expected_mirror_ack_sha256: str,
) -> dict[str, object]:
    local_anchor, mirror_anchor = genesis_anchor_paths(ledger_jsonl, mirror_jsonl)
    if _path_exists(local_anchor) and _path_exists(mirror_anchor):
        records = validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
        receipts = validate_receipts(
            receipts_jsonl,
            records,
            mirror_jsonl=mirror_jsonl,
            mirror_transport="filesystem_copy",
        )
        anchor = _genesis_anchor(records[0], receipts[0], mirror_jsonl=mirror_jsonl)
    else:
        local_records = validate_primary_chain(ledger_jsonl)
        mirror_records = validate_primary_chain(mirror_jsonl)
        if local_records != mirror_records:
            raise ValueError("Local and mirrored PIC ledger records differ")
        require_explicit_genesis(local_records)
        receipts = validate_receipts(
            receipts_jsonl,
            local_records,
            mirror_jsonl=mirror_jsonl,
            mirror_transport="filesystem_copy",
        )
        anchor = _genesis_anchor(
            local_records[0], receipts[0], mirror_jsonl=mirror_jsonl
        )
        if (
            anchor["event_sha256"] != expected_event_sha256
            or anchor["mirror_ack_sha256"] != expected_mirror_ack_sha256
        ):
            raise ValueError("Audited policy genesis anchor differs from ledger genesis")
        if _path_exists(local_anchor):
            if _read_one_genesis_anchor(
                local_anchor, ledger_jsonl.parent.parent
            ) != anchor:
                raise ValueError("Existing Orion genesis anchor differs from ledger")
        else:
            _write_one_genesis_anchor(
                local_anchor, ledger_jsonl.parent.parent, anchor
            )
        if _path_exists(mirror_anchor):
            if _read_one_genesis_anchor(
                mirror_anchor, mirror_jsonl.parent.parent
            ) != anchor:
                raise ValueError(
                    "Existing Project Home genesis anchor differs from ledger"
                )
        else:
            _write_one_genesis_anchor(
                mirror_anchor, mirror_jsonl.parent.parent, anchor
            )
        validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
    if (
        anchor["event_sha256"] != expected_event_sha256
        or anchor["mirror_ack_sha256"] != expected_mirror_ack_sha256
    ):
        raise ValueError("Existing genesis anchor differs from audited policy")
    return anchor


def _repair_mirrored_state_pinned(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    mirror_transport: str,
) -> dict[str, int]:
    """Repair only missing trailing mirror records or receipts from canonical Orion."""
    if mirror_transport != "filesystem_copy":
        raise ValueError("Only preflighted filesystem_copy transport is implemented")
    local_records = validate_primary_chain(ledger_jsonl)
    mirror_records = validate_primary_chain(mirror_jsonl)
    if mirror_records != local_records[:len(mirror_records)]:
        raise ValueError("Project Home ledger is not an exact prefix of canonical Orion")
    receipts = _read_jsonl(receipts_jsonl)
    for receipt in receipts:
        _validate_receipt_provenance(
            receipt,
            mirror_jsonl=mirror_jsonl,
            mirror_transport=mirror_transport,
        )
        expected_hash = record_sha256(receipt, "mirror_ack_sha256")
        if receipt.get("mirror_ack_sha256") != expected_hash:
            raise ValueError("Cannot repair a corrupt mirror receipt")
    receipt_hashes = [receipt.get("mirrored_event_sha256") for receipt in receipts]
    expected_receipt_prefix = [
        record["event_sha256"] for record in mirror_records[:len(receipts)]
    ]
    if receipt_hashes != expected_receipt_prefix:
        raise ValueError("Mirror receipts are not an exact prefix of mirrored events")
    if not receipts:
        raise ValueError("Cannot repair a mirrored ledger without its genesis receipt")
    _validate_genesis_anchors(ledger_jsonl, mirror_jsonl, local_records, receipts)
    _validate_operator_attestation_trees(
        local_records, ledger_jsonl=ledger_jsonl, mirror_jsonl=mirror_jsonl
    )
    mirror_preflight(mirror_jsonl)

    appended_mirror = 0
    for record in local_records[len(mirror_records):]:
        _append_jsonl(mirror_jsonl, record)
        mirror_records.append(record)
        appended_mirror += 1
    appended_receipts = 0
    for record in mirror_records[len(receipts):]:
        _append_mirror_receipt(receipts_jsonl, record, mirror_jsonl, mirror_transport)
        appended_receipts += 1
    validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
    write_csv(
        ledger_jsonl,
        receipts_jsonl,
        ledger_csv,
        mirror_jsonl=mirror_jsonl,
        mirror_transport=mirror_transport,
    )
    return {
        "appended_mirror_records": appended_mirror,
        "appended_receipts": appended_receipts,
    }


def repair_mirrored_state_locked(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    mirror_transport: str,
) -> dict[str, int]:
    require_no_incomplete_manual_accounting_marker(ledger_jsonl, mirror_jsonl)
    with _pinned_parent_directories(
        _ledger_state_paths(ledger_jsonl, ledger_csv, receipts_jsonl, mirror_jsonl)
    ):
        return _repair_mirrored_state_pinned(
            ledger_jsonl,
            ledger_csv,
            receipts_jsonl,
            mirror_jsonl,
            mirror_transport=mirror_transport,
        )


def repair_mirrored_state(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    mirror_transport: str,
) -> dict[str, int]:
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        return repair_mirrored_state_locked(
            ledger_jsonl,
            ledger_csv,
            receipts_jsonl,
            mirror_jsonl,
            mirror_transport=mirror_transport,
        )


def _append_primary_event_pinned(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    event: dict[str, object],
    *,
    mirror_transport: str,
    allow_incomplete_manual_accounting: bool,
) -> dict[str, object]:
    if mirror_transport != "filesystem_copy":
        raise ValueError("Only preflighted filesystem_copy transport is implemented")

    local_records = validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
    mirror_preflight(mirror_jsonl)
    mirror_records = local_records
    if event.get("event_type") == "genesis":
        if local_records or mirror_records:
            raise ValueError("Refusing to append a second Frontier PIC genesis event")
        if not event.get("control_plane_version"):
            raise ValueError("Genesis event must record the control-plane version")
    else:
        require_explicit_genesis(local_records)
        require_explicit_genesis(mirror_records)

    record = dict(event)
    record["sequence_number"] = len(local_records)
    record["previous_event_sha256"] = chain_head(local_records)
    record.setdefault("timestamp", utc_now())
    record["event_sha256"] = record_sha256(record, "event_sha256")
    _validate_accounting_records([*local_records, record])
    _validate_operator_attestation_trees(
        [*local_records, record], ledger_jsonl=ledger_jsonl, mirror_jsonl=mirror_jsonl
    )
    _validate_terminal_recovery_handoff_mirrors(
        [*local_records, record], ledger_jsonl=ledger_jsonl, mirror_jsonl=mirror_jsonl
    )
    if not allow_incomplete_manual_accounting:
        require_no_incomplete_manual_accounting_marker(ledger_jsonl, mirror_jsonl)
    _append_jsonl(ledger_jsonl, record)
    _append_jsonl(mirror_jsonl, record)
    receipt = _append_mirror_receipt(
        receipts_jsonl, record, mirror_jsonl, mirror_transport
    )
    if event.get("event_type") == "genesis":
        _write_genesis_anchors(
            ledger_jsonl,
            mirror_jsonl,
            _genesis_anchor(record, receipt, mirror_jsonl=mirror_jsonl),
        )
    write_csv(
        ledger_jsonl,
        receipts_jsonl,
        ledger_csv,
        mirror_jsonl=mirror_jsonl,
        mirror_transport=mirror_transport,
        allow_incomplete_manual_accounting=allow_incomplete_manual_accounting,
    )
    return record


def append_primary_event_locked(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    event: dict[str, object],
    *,
    mirror_transport: str,
    allow_incomplete_manual_accounting: bool = False,
) -> dict[str, object]:
    with _pinned_parent_directories(
        _ledger_state_paths(ledger_jsonl, ledger_csv, receipts_jsonl, mirror_jsonl)
    ):
        return _append_primary_event_pinned(
            ledger_jsonl,
            ledger_csv,
            receipts_jsonl,
            mirror_jsonl,
            event,
            mirror_transport=mirror_transport,
            allow_incomplete_manual_accounting=allow_incomplete_manual_accounting,
        )


def append_primary_event(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    event: dict[str, object],
    *,
    mirror_transport: str,
) -> dict[str, object]:
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        return append_primary_event_locked(
            ledger_jsonl,
            ledger_csv,
            receipts_jsonl,
            mirror_jsonl,
            event,
            mirror_transport=mirror_transport,
        )


def write_csv(
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    ledger_csv: Path,
    *,
    mirror_jsonl: Path,
    mirror_transport: str,
    allow_incomplete_manual_accounting: bool = False,
) -> None:
    if not allow_incomplete_manual_accounting:
        require_no_incomplete_manual_accounting_marker(ledger_jsonl, mirror_jsonl)
    primary = validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
    receipts = validate_receipts(
        receipts_jsonl,
        primary,
        mirror_jsonl=mirror_jsonl,
        mirror_transport=mirror_transport,
    )
    receipt_by_event = {
        str(receipt["mirrored_event_sha256"]): receipt for receipt in receipts
    }
    durable_mkdir_parents(ledger_csv.parent)
    ledger_csv = _require_canonical_path(ledger_csv)
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, extrasaction="ignore")
    writer.writeheader()
    for record in primary:
        row = dict(record)
        row.update(receipt_by_event.get(str(record["event_sha256"]), {}))
        writer.writerow({field: row.get(field, "") for field in CSV_FIELDS})
    parent_descriptor = _parent_descriptor(ledger_csv)
    if parent_descriptor is None:
        atomic_write_bytes(ledger_csv, stream.getvalue().encode("utf-8"), mode=0o600)
    else:
        atomic_write_bytes_at(
            parent_descriptor,
            ledger_csv.name,
            stream.getvalue().encode("utf-8"),
            mode=0o600,
        )


def initialize_ledger(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    mirror_transport: str,
    notes: str,
    control_plane_version: str,
) -> dict[str, object]:
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        with _pinned_parent_directories(
            _ledger_state_paths(ledger_jsonl, ledger_csv, receipts_jsonl, mirror_jsonl)
        ):
            anchors = genesis_anchor_paths(ledger_jsonl, mirror_jsonl)
            for path in [ledger_jsonl, ledger_csv, receipts_jsonl, mirror_jsonl, *anchors]:
                if _path_exists(path):
                    return _repair_interrupted_genesis_locked(
                        ledger_jsonl,
                        ledger_csv,
                        receipts_jsonl,
                        mirror_jsonl,
                        mirror_transport=mirror_transport,
                        notes=notes,
                        control_plane_version=control_plane_version,
                    )
            return append_primary_event_locked(
                ledger_jsonl,
                ledger_csv,
                receipts_jsonl,
                mirror_jsonl,
                {
                    "event_type": "genesis",
                    "state": "initialized",
                    "notes": notes,
                    "control_plane_version": control_plane_version,
                },
                mirror_transport=mirror_transport,
            )


def _repair_interrupted_genesis_locked(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    mirror_transport: str,
    notes: str,
    control_plane_version: str,
) -> dict[str, object]:
    """Complete only one exact interrupted fresh-genesis publication."""
    if mirror_transport != "filesystem_copy":
        raise ValueError("Only preflighted filesystem_copy transport is implemented")
    local_records = validate_primary_chain(ledger_jsonl)
    mirror_records = validate_primary_chain(mirror_jsonl)
    if len(local_records) != 1:
        raise ValueError("Refusing to overwrite existing ledger state")
    genesis = local_records[0]
    if (
        set(genesis)
        != {
            "event_type",
            "state",
            "notes",
            "control_plane_version",
            "sequence_number",
            "previous_event_sha256",
            "timestamp",
            "event_sha256",
        }
        or genesis.get("event_type") != "genesis"
        or genesis.get("state") != "initialized"
        or genesis.get("notes") != notes
        or genesis.get("control_plane_version") != control_plane_version
        or type(genesis.get("sequence_number")) is not int
        or genesis.get("sequence_number") != 0
        or genesis.get("previous_event_sha256") != ""
        or record_sha256(genesis, "event_sha256") != genesis.get("event_sha256")
    ):
        raise ValueError("Existing state is not the expected interrupted ledger genesis")
    if mirror_records not in [[], local_records]:
        raise ValueError("Existing Project Home state is not an interrupted genesis prefix")
    receipts = _read_jsonl(receipts_jsonl)
    if len(receipts) > 1 or (receipts and not mirror_records):
        raise ValueError("Existing receipt state is not an interrupted genesis prefix")
    for receipt in receipts:
        _validate_receipt_provenance(
            receipt,
            mirror_jsonl=mirror_jsonl,
            mirror_transport=mirror_transport,
        )
        if (
            receipt.get("mirrored_event_sha256") != genesis["event_sha256"]
            or receipt.get("mirror_ack_sha256")
            != record_sha256(receipt, "mirror_ack_sha256")
        ):
            raise ValueError("Existing genesis receipt is invalid")
    local_anchor, mirror_anchor = genesis_anchor_paths(ledger_jsonl, mirror_jsonl)
    if receipts:
        anchor = _genesis_anchor(genesis, receipts[0], mirror_jsonl=mirror_jsonl)
        for path, root in [
            (local_anchor, ledger_jsonl.parent.parent),
            (mirror_anchor, mirror_jsonl.parent.parent),
        ]:
            if _path_exists(path) and _read_one_genesis_anchor(path, root) != anchor:
                raise ValueError("Existing genesis anchor differs from interrupted genesis")
    elif _path_exists(local_anchor) or _path_exists(mirror_anchor):
        raise ValueError("Genesis anchor exists without a durable mirror receipt")

    mirror_preflight(mirror_jsonl)
    if not mirror_records:
        _append_jsonl(mirror_jsonl, genesis)
    if not receipts:
        receipts = [
            _append_mirror_receipt(
                receipts_jsonl, genesis, mirror_jsonl, mirror_transport
            )
        ]
    anchor = _genesis_anchor(genesis, receipts[0], mirror_jsonl=mirror_jsonl)
    if not _path_exists(local_anchor):
        _write_one_genesis_anchor(local_anchor, ledger_jsonl.parent.parent, anchor)
    if not _path_exists(mirror_anchor):
        _write_one_genesis_anchor(mirror_anchor, mirror_jsonl.parent.parent, anchor)
    validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
    write_csv(
        ledger_jsonl,
        receipts_jsonl,
        ledger_csv,
        mirror_jsonl=mirror_jsonl,
        mirror_transport=mirror_transport,
    )
    return genesis
