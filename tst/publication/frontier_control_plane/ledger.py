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
import os
from pathlib import Path
import stat
from typing import Iterator, TextIO
import uuid

from control_plane_common import atomic_write_bytes, atomic_write_bytes_at
from control_plane_common import atomic_write_json, atomic_write_json_at
from control_plane_common import durable_mkdir_parents, fsync_directory
from control_plane_common import read_json_bytes
from control_plane_common import read_stable_regular_file_below
from control_plane_common import stable_serialization_anchor


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
    "active_policy_sha256",
    "active_promotion_sha256",
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
GENESIS_ANCHOR_FILENAME = "genesis_anchor.json"
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
    records = _read_jsonl(path, root=root)
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
    root: Path | None = None,
) -> list[dict[str, object]]:
    records = _read_jsonl(path, root=root)
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
def _pinned_parent_directories(paths: list[Path]) -> Iterator[None]:
    existing = _PINNED_PARENT_DESCRIPTORS.get()
    opened: dict[Path, int] = {}
    try:
        for path in paths:
            parent = Path(os.path.abspath(path.parent))
            if parent in existing or parent in opened:
                continue
            durable_mkdir_parents(parent)
            descriptor = os.open(
                parent,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
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


@contextmanager
def ledger_lock(ledger_jsonl: Path, mirror_jsonl: Path | None = None) -> Iterator[None]:
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


def _read_one_genesis_anchor(path: Path, root: Path) -> dict[str, object]:
    return read_json_bytes(_read_one_genesis_anchor_bytes(path, root), label=str(path))


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
    if local_bytes != mirror_bytes:
        raise ValueError("Orion and Project Home genesis-anchor bytes differ")
    anchor = read_json_bytes(local_bytes, label=str(local_anchor))
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
    receipts = validate_receipts(
        receipts_jsonl,
        local_records,
        mirror_jsonl=mirror_jsonl,
        mirror_transport="filesystem_copy",
        root=receipts_root,
    )
    _validate_genesis_anchors(ledger_jsonl, mirror_jsonl, local_records, receipts)
    return local_records


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
) -> None:
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
