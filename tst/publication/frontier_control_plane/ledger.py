#!/opt/cray/pe/python/3.11.7/bin/python3
"""Hash-chained PIC node-hour ledger with a durable non-recursive mirror."""

from __future__ import annotations

from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Iterator, TextIO
import uuid


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
    "git_commit",
    "campaign",
    "test_id",
    "submission_scope",
    "clean_candidate_manifest_sha256",
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


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def canonical_json(record: dict[str, object], omit: str | None = None) -> str:
    payload = {key: value for key, value in record.items() if key != omit}
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def record_sha256(record: dict[str, object], omit: str) -> str:
    return hashlib.sha256(canonical_json(record, omit).encode("utf-8")).hexdigest()


def _read_jsonl(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    records = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw:
            raise ValueError(f"Blank JSONL line in {path}:{number}")
        record = json.loads(raw)
        if not isinstance(record, dict):
            raise ValueError(f"JSONL record is not an object in {path}:{number}")
        records.append(record)
    return records


def _append_jsonl(path: Path, record: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as stream:
        stream.write(canonical_json(record) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def validate_primary_chain(path: Path) -> list[dict[str, object]]:
    records = _read_jsonl(path)
    previous = ""
    for expected_sequence, record in enumerate(records):
        if record.get("sequence_number") != expected_sequence:
            raise ValueError(f"Invalid sequence number in {path}")
        if record.get("previous_event_sha256") != previous:
            raise ValueError(f"Broken previous-event link in {path}")
        expected_hash = record_sha256(record, "event_sha256")
        if record.get("event_sha256") != expected_hash:
            raise ValueError(f"Invalid event hash in {path}")
        previous = expected_hash
    return records


def _validate_receipt_provenance(
    receipt: dict[str, object],
    *,
    mirror_jsonl: Path,
    mirror_transport: str,
) -> None:
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
) -> list[dict[str, object]]:
    records = _read_jsonl(path)
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
        or not genesis.get("control_plane_version")
    ):
        raise ValueError("Frontier PIC ledger does not start with a valid genesis event")
    if any(record.get("event_type") == "genesis" for record in records[1:]):
        raise ValueError("Frontier PIC ledger contains more than one genesis event")


def accounting(records: list[dict[str, object]]) -> dict[str, float]:
    consumed = 0.0
    reserved = 0.0
    for record in latest_reservations(records).values():
        if bool(record.get("reconciled", False)):
            consumed += float(record.get("consumed_node_hours", 0.0))
        elif record.get("state") not in {"cancelled", "submission_attach_failed"}:
            reserved += float(record.get("reserved_node_hours", 0.0))
    return {
        "cumulative_consumed_node_hours": consumed,
        "currently_reserved_node_hours": reserved,
    }


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
) -> TextIO:
    descriptor = os.open(path, flags | os.O_NOFOLLOW, 0o600)
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise ValueError(f"Expected a regular file: {path}")
        return os.fdopen(descriptor, mode, encoding="utf-8", newline=newline)
    except BaseException:
        os.close(descriptor)
        raise


@contextmanager
def ledger_lock(ledger_jsonl: Path) -> Iterator[None]:
    lock_path = _require_canonical_path(
        ledger_jsonl.with_suffix(ledger_jsonl.suffix + ".lock")
    )
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = _require_canonical_path(lock_path)
    with _open_regular_nofollow(
        lock_path,
        "a",
        flags=os.O_APPEND | os.O_CREAT | os.O_WRONLY,
    ) as lock_stream:
        fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)


def mirror_preflight(mirror_jsonl: Path) -> None:
    mirror_jsonl.parent.mkdir(parents=True, exist_ok=True)
    probe = mirror_jsonl.parent / f".pic-ledger-probe-{os.getpid()}-{uuid.uuid4()}"
    probe.write_text("pic-ledger-mirror-preflight\n", encoding="utf-8")
    with probe.open("rb") as stream:
        os.fsync(stream.fileno())
    probe.unlink()


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


def validate_mirrored_state(
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
) -> list[dict[str, object]]:
    local_records = validate_primary_chain(ledger_jsonl)
    mirror_records = validate_primary_chain(mirror_jsonl)
    if local_records != mirror_records:
        raise ValueError("Local and mirrored PIC ledger records differ")
    validate_receipts(
        receipts_jsonl,
        local_records,
        mirror_jsonl=mirror_jsonl,
        mirror_transport="filesystem_copy",
    )
    return local_records


def repair_mirrored_state_locked(
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
    mirror_preflight(mirror_jsonl)
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


def repair_mirrored_state(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    mirror_transport: str,
) -> dict[str, int]:
    with ledger_lock(ledger_jsonl):
        return repair_mirrored_state_locked(
            ledger_jsonl,
            ledger_csv,
            receipts_jsonl,
            mirror_jsonl,
            mirror_transport=mirror_transport,
        )


def append_primary_event_locked(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    event: dict[str, object],
    *,
    mirror_transport: str,
) -> dict[str, object]:
    if mirror_transport != "filesystem_copy":
        raise ValueError("Only preflighted filesystem_copy transport is implemented")

    mirror_preflight(mirror_jsonl)
    local_records = validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
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
    _append_jsonl(ledger_jsonl, record)
    _append_jsonl(mirror_jsonl, record)
    _append_mirror_receipt(receipts_jsonl, record, mirror_jsonl, mirror_transport)
    write_csv(
        ledger_jsonl,
        receipts_jsonl,
        ledger_csv,
        mirror_jsonl=mirror_jsonl,
        mirror_transport=mirror_transport,
    )
    return record


def append_primary_event(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    event: dict[str, object],
    *,
    mirror_transport: str,
) -> dict[str, object]:
    with ledger_lock(ledger_jsonl):
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
) -> None:
    primary = validate_primary_chain(ledger_jsonl)
    receipts = validate_receipts(
        receipts_jsonl,
        primary,
        mirror_jsonl=mirror_jsonl,
        mirror_transport=mirror_transport,
    )
    receipt_by_event = {
        str(receipt["mirrored_event_sha256"]): receipt for receipt in receipts
    }
    ledger_csv.parent.mkdir(parents=True, exist_ok=True)
    ledger_csv = _require_canonical_path(ledger_csv)
    temporary = _require_canonical_path(
        ledger_csv.with_suffix(ledger_csv.suffix + ".tmp")
    )
    with _open_regular_nofollow(
        temporary,
        "w",
        flags=os.O_CREAT | os.O_TRUNC | os.O_WRONLY,
        newline="",
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for record in primary:
            row = dict(record)
            row.update(receipt_by_event.get(str(record["event_sha256"]), {}))
            writer.writerow({field: row.get(field, "") for field in CSV_FIELDS})
        stream.flush()
        os.fsync(stream.fileno())
    temporary = _require_canonical_path(temporary)
    ledger_csv = _require_canonical_path(ledger_csv)
    os.replace(temporary, ledger_csv)


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
    with ledger_lock(ledger_jsonl):
        for path in [ledger_jsonl, ledger_csv, receipts_jsonl, mirror_jsonl]:
            if path.exists() and path.stat().st_size:
                raise ValueError(f"Refusing to overwrite existing ledger state: {path}")
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
