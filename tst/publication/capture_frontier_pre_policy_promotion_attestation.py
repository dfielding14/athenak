#!/usr/bin/env python3
"""Capture and seal one reviewed Frontier same-account isolation attestation."""

from __future__ import annotations

import argparse
import ctypes
from dataclasses import dataclass
from datetime import datetime, timezone
import errno
import hashlib
import importlib
import json
import os
from pathlib import Path
import pwd
import re
import secrets
import stat
import subprocess
import sys
from typing import Callable, Sequence


AUTHORIZED_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
AUTHORIZED_PROJECT_HOME_ROOT = Path("/ccs/proj/ast207/proj-shared/PIC")
DEFAULT_ARCHIVE_ROOT = AUTHORIZED_PIC_ROOT / "operator_attestations"
TRUSTED_PS = "/usr/bin/ps"
TRUSTED_SQUEUE = "/usr/bin/squeue"
QUEUE_FORMAT = "%i|%a|%P|%q|%T|%j|%k"
PHASES = ("pre_policy_promotion", "pre_manifest", "pre_submit_wrapper")
PHASE = "pre_policy_promotion"
RECORD_TYPE = "q027_frontier_registered_science_same_account_isolation_attestation"
OPERATOR_STATEMENTS = {
    "pre_policy_promotion": (
        "Reviewed same-account process snapshot immediately before registered policy "
        "promotion; no competing process is authorized to mutate the PIC roots "
        "throughout the paired mirrored-policy replacement."
    ),
    "pre_manifest": (
        "Reviewed same-account process snapshot immediately before registered manifest "
        "creation; no competing process is authorized to mutate the PIC roots "
        "throughout manifest snapshot publication."
    ),
    "pre_submit_wrapper": (
        "Reviewed same-account process snapshot immediately before the registered "
        "submission wrapper; no competing process is authorized to mutate the PIC "
        "roots throughout reservation, scheduler dispatch, attachment and release."
    ),
}
OPERATOR_STATEMENT = OPERATOR_STATEMENTS[PHASE]
PENDING_SUBMISSION_FILENAME = "pending_submission.json"
PENDING_MANUAL_ACCOUNTING_FILENAME = "pending_manual_accounting.json"
METADATA_FILENAME = ".capture_metadata.json"
ACTIVE_RESERVATION_STATES = {"reserved", "submitted"}
AUTHORIZATION_ID_PATTERN = re.compile(r"[A-Za-z0-9](?:[A-Za-z0-9_-]{0,126}[A-Za-z0-9])?")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
CAPTURE_TO_SEAL_MAX_AGE_SECONDS = 15 * 60
AT_FDCWD = -100
RENAME_NOREPLACE = 1
GATE_SNAPSHOT_FILENAMES = {
    "queue_snapshot.txt",
    "pending_submission_marker.txt",
    "pending_manual_accounting_marker.txt",
    "mirrored_ledger_line_counts.txt",
    "validated_mirrored_ledger_state.json",
}
CAPTURED_GATE_SNAPSHOT_FILENAMES = {
    f"capture_{filename}" for filename in GATE_SNAPSHOT_FILENAMES
}
CAPTURED_PROCESS_SNAPSHOT_FILENAME = "capture_same_account_process_snapshot.txt"
CAPTURED_SNAPSHOT_FILENAMES = {
    CAPTURED_PROCESS_SNAPSHOT_FILENAME,
    *CAPTURED_GATE_SNAPSHOT_FILENAMES,
}
REVIEW_STAGING_FILENAMES = {
    METADATA_FILENAME,
    "same_account_process_snapshot.txt",
    *GATE_SNAPSHOT_FILENAMES,
    *CAPTURED_SNAPSHOT_FILENAMES,
}
SEALED_FILENAMES = {
    "attestation.json",
    "same_account_process_snapshot.txt",
    *GATE_SNAPSHOT_FILENAMES,
    *CAPTURED_SNAPSHOT_FILENAMES,
}


@dataclass(frozen=True)
class LedgerPaths:
    """Paths inspected without mutating scheduler or ledger state."""

    ledger: Path
    receipts: Path
    mirror: Path
    pending_submission: Path
    pending_manual_accounting: Path
    mirror_pending_manual_accounting: Path


@dataclass(frozen=True)
class GateSnapshot:
    """One stable read-only scheduler and ledger observation."""

    queue: bytes
    pending_submission_status: str
    pending_manual_accounting_statuses: dict[str, str]
    line_counts: dict[str, int]
    line_counts_payload: bytes
    validated_state: dict[str, object]
    validated_state_payload: bytes


Runner = Callable[[Sequence[str]], subprocess.CompletedProcess[str]]
LedgerValidator = Callable[[LedgerPaths], list[dict[str, object]]]
Clock = Callable[[], datetime]


def _default_runner(argv: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(argv, check=False, capture_output=True, text=True)


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _current_user() -> str:
    return pwd.getpwuid(os.getuid()).pw_name


def _validate_user(user: str) -> str:
    if not user or any(character in user for character in "\x00\r\n"):
        raise ValueError("Same-account user is malformed")
    return user


def _validate_authorization_id(authorization_id: str) -> str:
    if AUTHORIZATION_ID_PATTERN.fullmatch(authorization_id) is None:
        raise ValueError(
            "Authorization ID must contain only ASCII letters, digits, '_' or '-', "
            "must start and end with an alphanumeric character, and must be at most "
            "128 characters"
        )
    return authorization_id


def _validate_phase(phase: str) -> str:
    if phase not in PHASES:
        raise ValueError(f"Phase must be one of: {', '.join(PHASES)}")
    return phase


def _validate_control_plane_version(control_plane_version: str) -> str:
    if SHA256_PATTERN.fullmatch(control_plane_version) is None:
        raise ValueError("Control-plane version must be one lowercase SHA-256 digest")
    return control_plane_version


def _timestamp_strings(now: datetime) -> tuple[str, str]:
    if now.tzinfo is None:
        raise ValueError("Injected clock must return a timezone-aware datetime")
    normalized = now.astimezone(timezone.utc).replace(microsecond=0)
    return (
        normalized.strftime("%Y%m%dT%H%M%SZ"),
        normalized.isoformat().replace("+00:00", "Z"),
    )


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _canonical_existing_directory(
    path: Path, *, label: str, trusted_lexical_alias: Path | None = None
) -> Path:
    lexical = Path(os.path.abspath(path))
    try:
        resolved = lexical.resolve(strict=True)
    except FileNotFoundError as error:
        raise ValueError(f"{label} does not exist: {lexical}") from error
    allowed_alias = (
        trusted_lexical_alias is not None
        and lexical == Path(os.path.abspath(trusted_lexical_alias))
    )
    if lexical != resolved and not allowed_alias:
        raise ValueError(f"{label} must not use a symlink alias: {lexical}")
    metadata = os.stat(lexical) if allowed_alias else os.lstat(lexical)
    if not stat.S_ISDIR(metadata.st_mode):
        raise ValueError(f"{label} is not a directory: {lexical}")
    return lexical


def _ensure_archive_root(path: Path) -> Path:
    lexical = Path(os.path.abspath(path))
    try:
        return _canonical_existing_directory(lexical, label="Archive root")
    except ValueError:
        if os.path.lexists(lexical):
            raise
    parent = _canonical_existing_directory(lexical.parent, label="Archive-root parent")
    try:
        os.mkdir(lexical, 0o700)
    except FileExistsError as error:
        raise ValueError(f"Archive-root path collided during creation: {lexical}") from error
    _fsync_directory(parent)
    return _canonical_existing_directory(lexical, label="Archive root")


def _ledger_paths(pic_root: Path, project_home_root: Path) -> LedgerPaths:
    pic_root = _canonical_existing_directory(pic_root, label="PIC root")
    project_home_root = _canonical_existing_directory(
        project_home_root,
        label="Project Home root",
        trusted_lexical_alias=AUTHORIZED_PROJECT_HOME_ROOT,
    )
    ledger_root = _canonical_existing_directory(pic_root / "ledger", label="PIC ledger root")
    mirror_ledger_root = _canonical_existing_directory(
        project_home_root / "ledger",
        label="Project Home ledger root",
        trusted_lexical_alias=AUTHORIZED_PROJECT_HOME_ROOT / "ledger",
    )
    return LedgerPaths(
        ledger=ledger_root / "node_hours.jsonl",
        receipts=ledger_root / "mirror_receipts.jsonl",
        mirror=mirror_ledger_root / "node_hours.jsonl",
        pending_submission=ledger_root / PENDING_SUBMISSION_FILENAME,
        pending_manual_accounting=ledger_root / PENDING_MANUAL_ACCOUNTING_FILENAME,
        mirror_pending_manual_accounting=(
            mirror_ledger_root / PENDING_MANUAL_ACCOUNTING_FILENAME
        ),
    )


def _run_capture_command(runner: Runner, argv: Sequence[str]) -> bytes:
    result = runner(list(argv))
    if result.returncode != 0:
        stderr = result.stderr.strip() if isinstance(result.stderr, str) else ""
        detail = f": {stderr}" if stderr else ""
        raise ValueError(f"Read-only capture command failed ({argv[0]}){detail}")
    if not isinstance(result.stdout, str):
        raise ValueError(f"Read-only capture command returned non-text output: {argv[0]}")
    return result.stdout.encode("utf-8")


def _capture_same_account_processes(runner: Runner, user: str) -> bytes:
    return _run_capture_command(
        runner,
        [TRUSTED_PS, "-u", user, "-o", "pid=,ppid=,state=,args="],
    )


def _capture_queue(runner: Runner, user: str) -> bytes:
    return _run_capture_command(
        runner,
        [TRUSTED_SQUEUE, "-u", user, "-h", "-o", QUEUE_FORMAT],
    )


def _marker_status(path: Path) -> str:
    try:
        metadata = os.lstat(path)
    except FileNotFoundError:
        return "absent"
    if stat.S_ISLNK(metadata.st_mode):
        raise ValueError(f"Pending-marker path must not be a symlink: {path}")
    return "present"


def _read_regular_bytes(path: Path) -> bytes:
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"Expected a regular file: {path}")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                return b"".join(chunks)
            chunks.append(chunk)
    finally:
        os.close(descriptor)


def _line_count(path: Path) -> int:
    return _read_regular_bytes(path).count(b"\n")


def _line_counts(paths: LedgerPaths) -> dict[str, int]:
    return {
        str(path): _line_count(path)
        for path in [paths.ledger, paths.receipts, paths.mirror]
    }


def _line_counts_payload(counts: dict[str, int], paths: LedgerPaths) -> bytes:
    return "".join(
        f"{counts[str(path)]} {path}\n"
        for path in [paths.ledger, paths.receipts, paths.mirror]
    ).encode("utf-8")


def _load_source_ledger_module() -> object:
    control_plane_dir = Path(__file__).absolute().parent / "frontier_control_plane"
    expected = (control_plane_dir / "ledger.py").resolve()
    sys.path.insert(0, str(control_plane_dir))
    try:
        module = importlib.import_module("ledger")
    finally:
        sys.path.pop(0)
    loaded = Path(str(getattr(module, "__file__", ""))).resolve()
    if loaded != expected:
        raise ValueError(
            f"Refusing unexpected ledger validator module {loaded}; expected {expected}"
        )
    return module


def _validate_with_existing_read_only_ledger_snapshot(
    paths: LedgerPaths,
) -> list[dict[str, object]]:
    ledger_module = _load_source_ledger_module()
    validator = getattr(
        ledger_module, "validated_read_only_mirrored_state_snapshot", None
    )
    if validator is None:
        raise ValueError("Existing read-only mirrored-ledger snapshot validator is unavailable")
    with validator(
        paths.ledger,
        paths.receipts,
        paths.mirror,
        ledger_root=paths.ledger.parent.parent,
        receipts_root=paths.receipts.parent.parent,
        mirror_root=paths.mirror.parent.parent,
    ) as records:
        return [dict(record) for record in records]


def _active_reservation_ids(records: list[dict[str, object]]) -> list[str]:
    latest: dict[str, dict[str, object]] = {}
    for record in records:
        reservation_id = record.get("reservation_id")
        if isinstance(reservation_id, str) and reservation_id:
            latest[reservation_id] = record
    return sorted(
        reservation_id
        for reservation_id, record in latest.items()
        if record.get("state") in ACTIVE_RESERVATION_STATES
    )


def _collect_ledger_snapshot(
    paths: LedgerPaths, ledger_validator: LedgerValidator
) -> tuple[dict[str, int], bytes, dict[str, object], bytes]:
    before = _line_counts(paths)
    records = ledger_validator(paths)
    after = _line_counts(paths)
    if before != after:
        raise ValueError("Mirrored ledger line counts changed during validation")
    if len(set(after.values())) != 1:
        raise ValueError("Mirrored ledger line counts differ")
    if after[str(paths.ledger)] != len(records):
        raise ValueError("Validated mirrored-ledger record count differs from captured lines")
    active_reservation_ids = _active_reservation_ids(records)
    state: dict[str, object] = {
        "validation": "coherent",
        "validator": "existing_read_only_mirrored_state_snapshot",
        "ledger_record_count": len(records),
        "ledger_tail_event_sha256": (
            records[-1].get("event_sha256") if records else None
        ),
        "active_reservation_count": len(active_reservation_ids),
        "active_reservation_ids": active_reservation_ids,
    }
    state_payload = (
        json.dumps(state, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    return after, _line_counts_payload(after, paths), state, state_payload


def _capture_gate_snapshot(
    paths: LedgerPaths,
    *,
    runner: Runner,
    user: str,
    ledger_validator: LedgerValidator,
) -> GateSnapshot:
    queue = _capture_queue(runner, user)
    pending_submission_status = _marker_status(paths.pending_submission)
    pending_manual_accounting_statuses = {
        str(paths.pending_manual_accounting): _marker_status(
            paths.pending_manual_accounting
        ),
        str(paths.mirror_pending_manual_accounting): _marker_status(
            paths.mirror_pending_manual_accounting
        ),
    }
    counts, counts_payload, state, state_payload = _collect_ledger_snapshot(
        paths, ledger_validator
    )
    if pending_submission_status != _marker_status(paths.pending_submission):
        raise ValueError("Pending-submission marker changed during capture")
    if pending_manual_accounting_statuses != {
        str(paths.pending_manual_accounting): _marker_status(
            paths.pending_manual_accounting
        ),
        str(paths.mirror_pending_manual_accounting): _marker_status(
            paths.mirror_pending_manual_accounting
        ),
    }:
        raise ValueError("Pending-manual-accounting marker changed during capture")
    return GateSnapshot(
        queue=queue,
        pending_submission_status=pending_submission_status,
        pending_manual_accounting_statuses=pending_manual_accounting_statuses,
        line_counts=counts,
        line_counts_payload=counts_payload,
        validated_state=state,
        validated_state_payload=state_payload,
    )


def _manual_marker_payload(statuses: dict[str, str]) -> bytes:
    return "".join(f"{path} {status}\n" for path, status in statuses.items()).encode(
        "utf-8"
    )


def _gate_snapshot_files(snapshot: GateSnapshot) -> dict[str, bytes]:
    return {
        "queue_snapshot.txt": snapshot.queue,
        "pending_submission_marker.txt": (
            snapshot.pending_submission_status + "\n"
        ).encode("ascii"),
        "pending_manual_accounting_marker.txt": _manual_marker_payload(
            snapshot.pending_manual_accounting_statuses
        ),
        "mirrored_ledger_line_counts.txt": snapshot.line_counts_payload,
        "validated_mirrored_ledger_state.json": snapshot.validated_state_payload,
    }


def _write_new_file(path: Path, payload: bytes) -> None:
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
        0o600,
    )
    try:
        with os.fdopen(descriptor, "wb", closefd=False) as stream:
            stream.write(payload)
            stream.flush()
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _write_json_new(path: Path, value: dict[str, object]) -> None:
    _write_new_file(
        path,
        (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
            "utf-8"
        ),
    )


def _atomic_replace_file(path: Path, payload: bytes) -> None:
    metadata = os.lstat(path)
    if not stat.S_ISREG(metadata.st_mode):
        raise ValueError(f"Staged snapshot path is not a regular file: {path}")
    temporary = path.parent / f".{path.name}.tmp-{secrets.token_hex(8)}"
    _write_new_file(temporary, payload)
    os.replace(temporary, path)
    _fsync_directory(path.parent)


def _require_absent_entry(path: Path, *, label: str) -> None:
    if os.path.lexists(path):
        raise ValueError(f"{label} collides with an existing path: {path}")


def _new_staging_directory(archive_root: Path, final_name: str) -> Path:
    for _ in range(8):
        staging = archive_root / f".{final_name}.staging-{secrets.token_hex(8)}"
        try:
            os.mkdir(staging, 0o700)
        except FileExistsError:
            continue
        _fsync_directory(archive_root)
        return staging
    raise ValueError("Could not allocate a unique hidden staging directory")


def _metadata_bytes(metadata: dict[str, object]) -> bytes:
    return (
        json.dumps(metadata, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _read_json_regular(path: Path) -> dict[str, object]:
    try:
        value = json.loads(_read_regular_bytes(path))
    except json.JSONDecodeError as error:
        raise ValueError(f"Malformed JSON file: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def capture(
    *,
    archive_root: Path,
    authorization_id: str,
    control_plane_version: str,
    phase: str = PHASE,
    pic_root: Path = AUTHORIZED_PIC_ROOT,
    project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    runner: Runner = _default_runner,
    now: Clock = _utc_now,
    user: str | None = None,
    ledger_validator: LedgerValidator = _validate_with_existing_read_only_ledger_snapshot,
) -> Path:
    """Capture review inputs into one hidden staging directory."""
    authorization_id = _validate_authorization_id(authorization_id)
    control_plane_version = _validate_control_plane_version(control_plane_version)
    phase = _validate_phase(phase)
    user = _validate_user(_current_user() if user is None else user)
    pic_root = _canonical_existing_directory(pic_root, label="PIC root")
    project_home_root = _canonical_existing_directory(
        project_home_root,
        label="Project Home root",
        trusted_lexical_alias=AUTHORIZED_PROJECT_HOME_ROOT,
    )
    archive_root = _ensure_archive_root(archive_root)
    paths = _ledger_paths(pic_root, project_home_root)
    timestamp, recorded_utc = _timestamp_strings(now())
    final_name = f"{timestamp}-{authorization_id}-{phase}"
    final_path = archive_root / final_name
    _require_absent_entry(final_path, label="Final attestation directory")

    process_snapshot = _capture_same_account_processes(runner, user)
    gate_snapshot = _capture_gate_snapshot(
        paths, runner=runner, user=user, ledger_validator=ledger_validator
    )
    staging = _new_staging_directory(archive_root, final_name)
    metadata: dict[str, object] = {
        "schema_version": 1,
        "registered_science_authorization_id": authorization_id,
        "control_plane_version": control_plane_version,
        "phase": phase,
        "recorded_utc": recorded_utc,
        "final_directory_name": final_name,
        "staging_directory_name": staging.name,
        "pic_root": str(pic_root),
        "project_home_root": str(project_home_root),
        "same_account_user": user,
    }
    _write_new_file(staging / METADATA_FILENAME, _metadata_bytes(metadata))
    _write_new_file(staging / "same_account_process_snapshot.txt", process_snapshot)
    _write_new_file(staging / CAPTURED_PROCESS_SNAPSHOT_FILENAME, process_snapshot)
    for filename, payload in _gate_snapshot_files(gate_snapshot).items():
        _write_new_file(staging / filename, payload)
        _write_new_file(staging / f"capture_{filename}", payload)
    _fsync_directory(staging)
    return staging


def _validate_metadata(metadata: dict[str, object], staging: Path) -> str:
    required = {
        "schema_version",
        "registered_science_authorization_id",
        "control_plane_version",
        "recorded_utc",
        "final_directory_name",
        "staging_directory_name",
        "pic_root",
        "project_home_root",
        "same_account_user",
    }
    # Accept legacy hidden staging captures so an upgrade does not strand an
    # already-reviewed pre-policy-promotion capture.
    allowed_fields = [required, required | {"phase"}]
    string_fields = required - {"schema_version"}
    if "phase" in metadata:
        string_fields.add("phase")
    if (
        set(metadata) not in allowed_fields
        or type(metadata.get("schema_version")) is not int
        or metadata.get("schema_version") != 1
        or any(type(metadata.get(field)) is not str for field in string_fields)
    ):
        raise ValueError("Staged capture metadata schema is invalid")
    authorization_id = _validate_authorization_id(
        str(metadata["registered_science_authorization_id"])
    )
    phase = _validate_phase(str(metadata.get("phase", PHASE)))
    _validate_control_plane_version(str(metadata["control_plane_version"]))
    _validate_user(str(metadata["same_account_user"]))
    timestamp = str(metadata["final_directory_name"]).split("-", 1)[0]
    if re.fullmatch(r"[0-9]{8}T[0-9]{6}Z", timestamp) is None:
        raise ValueError("Staged final directory timestamp is malformed")
    try:
        recorded = datetime.strptime(
            str(metadata["recorded_utc"]), "%Y-%m-%dT%H:%M:%SZ"
        )
    except ValueError as error:
        raise ValueError("Staged recorded UTC timestamp is malformed") from error
    if recorded.strftime("%Y%m%dT%H%M%SZ") != timestamp:
        raise ValueError("Staged recorded UTC timestamp differs from directory name")
    expected_final_name = f"{timestamp}-{authorization_id}-{phase}"
    if metadata["final_directory_name"] != expected_final_name:
        raise ValueError("Staged final directory name differs from capture metadata")
    if metadata["staging_directory_name"] != staging.name:
        raise ValueError("Staged directory name differs from capture metadata")
    if not staging.name.startswith(f".{expected_final_name}.staging-"):
        raise ValueError("Staged directory name is malformed")
    return phase


def _require_matching_optional_root(
    captured: Path, override: Path | None, *, label: str
) -> Path:
    trusted_alias = (
        AUTHORIZED_PROJECT_HOME_ROOT if label == "Project Home root" else None
    )
    captured = _canonical_existing_directory(
        captured, label=label, trusted_lexical_alias=trusted_alias
    )
    if override is None:
        return captured
    override = _canonical_existing_directory(
        override, label=label, trusted_lexical_alias=trusted_alias
    )
    if override != captured:
        raise ValueError(f"{label} override differs from staged capture")
    return captured


def _require_sealable(snapshot: GateSnapshot) -> None:
    if snapshot.queue.strip():
        raise ValueError("Isolation attestation is blocked by a non-empty queue")
    if snapshot.pending_submission_status != "absent":
        raise ValueError("Isolation attestation is blocked by pending submission")
    if any(
        status != "absent"
        for status in snapshot.pending_manual_accounting_statuses.values()
    ):
        raise ValueError("Isolation attestation is blocked by pending manual accounting")
    if snapshot.validated_state["active_reservation_count"] != 0:
        raise ValueError("Isolation attestation is blocked by an outstanding reservation")


def _require_fresh_capture(metadata: dict[str, object], sealed: datetime) -> None:
    if sealed.tzinfo is None:
        raise ValueError("Injected clock must return a timezone-aware datetime")
    try:
        recorded = datetime.strptime(
            str(metadata["recorded_utc"]), "%Y-%m-%dT%H:%M:%SZ"
        ).replace(tzinfo=timezone.utc)
    except ValueError as error:
        raise ValueError("Staged recorded UTC timestamp is malformed") from error
    sealed = sealed.astimezone(timezone.utc)
    if (
        sealed < recorded
        or (sealed - recorded).total_seconds() > CAPTURE_TO_SEAL_MAX_AGE_SECONDS
    ):
        raise ValueError("Staged capture-to-seal interval exceeds fifteen minutes")


def _sha256(path: Path) -> str:
    return hashlib.sha256(_read_regular_bytes(path)).hexdigest()


def _file_record(staging: Path, filename: str) -> dict[str, object]:
    return {"path": filename, "sha256": _sha256(staging / filename)}


def _attestation(
    staging: Path,
    metadata: dict[str, object],
    snapshot: GateSnapshot,
    *,
    phase: str,
    sealed_utc: str,
) -> dict[str, object]:
    pending_submission = _file_record(staging, "pending_submission_marker.txt")
    pending_submission["value"] = snapshot.pending_submission_status
    pending_manual_accounting = _file_record(
        staging, "pending_manual_accounting_marker.txt"
    )
    pending_manual_accounting["value"] = "absent"
    pending_manual_accounting["values"] = snapshot.pending_manual_accounting_statuses
    line_counts = _file_record(staging, "mirrored_ledger_line_counts.txt")
    line_counts["counts"] = snapshot.line_counts
    validated_state = _file_record(staging, "validated_mirrored_ledger_state.json")
    validated_state["state"] = snapshot.validated_state
    return {
        "schema_version": 1,
        "record_type": RECORD_TYPE,
        "recorded_utc": metadata["recorded_utc"],
        "sealed_utc": sealed_utc,
        "registered_science_authorization_id": metadata[
            "registered_science_authorization_id"
        ],
        "control_plane_version": metadata["control_plane_version"],
        "phase": phase,
        "same_account_process_snapshot": _file_record(
            staging, "same_account_process_snapshot.txt"
        ),
        "queue_snapshot": _file_record(staging, "queue_snapshot.txt"),
        "pending_submission_marker": pending_submission,
        "pending_manual_accounting_marker": pending_manual_accounting,
        "mirrored_ledger_line_counts": line_counts,
        "validated_mirrored_ledger_state": validated_state,
        "captured_snapshots": {
            filename: _file_record(staging, filename)
            for filename in sorted(CAPTURED_SNAPSHOT_FILENAMES)
        },
        "operator_statement": OPERATOR_STATEMENTS[phase],
    }


def _walk_tree(path: Path) -> list[Path]:
    paths = [path]
    with os.scandir(path) as entries:
        for entry in entries:
            child = Path(entry.path)
            metadata = os.lstat(child)
            if stat.S_ISLNK(metadata.st_mode):
                raise ValueError(f"Attestation tree must not contain symlinks: {child}")
            if stat.S_ISDIR(metadata.st_mode):
                paths.extend(_walk_tree(child))
            elif not stat.S_ISREG(metadata.st_mode):
                raise ValueError(f"Attestation tree contains a non-regular file: {child}")
            else:
                paths.append(child)
    return paths


def _require_exact_regular_files(path: Path, expected: set[str]) -> None:
    actual: set[str] = set()
    with os.scandir(path) as entries:
        for entry in entries:
            child = Path(entry.path)
            metadata = os.lstat(child)
            if not stat.S_ISREG(metadata.st_mode):
                raise ValueError(f"Staged attestation member is not a regular file: {child}")
            actual.add(entry.name)
    if actual != expected:
        raise ValueError(
            "Staged attestation members differ: "
            f"expected={sorted(expected)}, actual={sorted(actual)}"
        )


def _fsync_tree(path: Path) -> None:
    for member in _walk_tree(path):
        metadata = os.lstat(member)
        if stat.S_ISREG(metadata.st_mode):
            descriptor = os.open(member, os.O_RDONLY | os.O_NOFOLLOW)
            try:
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
    for member in reversed(_walk_tree(path)):
        if stat.S_ISDIR(os.lstat(member).st_mode):
            _fsync_directory(member)


def _make_tree_read_only(path: Path) -> None:
    members = _walk_tree(path)
    for member in members:
        if stat.S_ISREG(os.lstat(member).st_mode):
            os.chmod(member, 0o400, follow_symlinks=False)
    for member in reversed(members):
        if stat.S_ISDIR(os.lstat(member).st_mode):
            os.chmod(member, 0o500, follow_symlinks=False)


def _require_tree_read_only(path: Path) -> None:
    for member in _walk_tree(path):
        metadata = os.lstat(member)
        expected_mode = 0o500 if stat.S_ISDIR(metadata.st_mode) else 0o400
        if stat.S_IMODE(metadata.st_mode) != expected_mode:
            raise ValueError(f"Sealed attestation tree member is not read-only: {member}")


def _renameat2_no_replace(source: Path, destination: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise OSError(errno.ENOSYS, os.strerror(errno.ENOSYS), str(destination))
    renameat2.argtypes = [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    ]
    renameat2.restype = ctypes.c_int
    result = renameat2(
        AT_FDCWD,
        os.fsencode(source),
        AT_FDCWD,
        os.fsencode(destination),
        RENAME_NOREPLACE,
    )
    if result == 0:
        return
    error_number = ctypes.get_errno()
    raise OSError(error_number, os.strerror(error_number), str(destination))


def _rename_no_replace(source: Path, destination: Path) -> None:
    try:
        _renameat2_no_replace(source, destination)
        return
    except OSError as error:
        error_number = error.errno
    if error_number == errno.EEXIST:
        raise ValueError(f"Final attestation directory collided during rename: {destination}")
    unsupported = {
        errno.EINVAL,
        errno.ENOSYS,
        getattr(errno, "ENOTSUP", errno.EINVAL),
        getattr(errno, "EOPNOTSUPP", errno.EINVAL),
    }
    if error_number not in unsupported:
        raise

    # Lustre rejects renameat2(RENAME_NOREPLACE). The operator attestation
    # establishes the same-account isolation boundary for this short fallback.
    before = os.lstat(source)
    _require_absent_entry(destination, label="Final attestation directory")
    os.rename(source, destination)
    after = os.lstat(destination)
    if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
        raise RuntimeError("Final attestation directory identity changed during rename")
    if os.path.lexists(source):
        raise RuntimeError("Staged attestation directory still exists after rename")


def _sealed_final_path(staging: Path) -> Path:
    attestation = _read_json_regular(staging / "attestation.json")
    required = {
        "schema_version",
        "record_type",
        "recorded_utc",
        "sealed_utc",
        "registered_science_authorization_id",
        "control_plane_version",
        "phase",
    }
    if (
        not required.issubset(attestation)
        or attestation.get("schema_version") != 1
        or attestation.get("record_type") != RECORD_TYPE
        or any(type(attestation.get(field)) is not str for field in required - {"schema_version"})
    ):
        raise ValueError("Sealed attestation schema is invalid")
    authorization_id = _validate_authorization_id(
        str(attestation["registered_science_authorization_id"])
    )
    phase = _validate_phase(str(attestation["phase"]))
    _validate_control_plane_version(str(attestation["control_plane_version"]))
    try:
        recorded = datetime.strptime(str(attestation["recorded_utc"]), "%Y-%m-%dT%H:%M:%SZ")
        datetime.strptime(str(attestation["sealed_utc"]), "%Y-%m-%dT%H:%M:%SZ")
    except ValueError as error:
        raise ValueError("Sealed attestation UTC timestamp is malformed") from error
    final_name = f"{recorded.strftime('%Y%m%dT%H%M%SZ')}-{authorization_id}-{phase}"
    if re.fullmatch(rf"\.{re.escape(final_name)}\.staging-[0-9a-f]{{16}}", staging.name) is None:
        raise ValueError("Sealed staging directory name is malformed")
    return staging.parent / final_name


def publish_sealed_staging(staging_dir: Path) -> Path:
    """Publish a previously revalidated immutable staging tree after rename failure."""
    staging = _canonical_existing_directory(staging_dir, label="Sealed staging directory")
    archive_root = _canonical_existing_directory(staging.parent, label="Archive root")
    _require_exact_regular_files(staging, SEALED_FILENAMES)
    _require_tree_read_only(staging)
    final = _sealed_final_path(staging)
    _require_absent_entry(final, label="Final attestation directory")
    _fsync_tree(staging)
    _rename_no_replace(staging, final)
    _fsync_directory(archive_root)
    return final


def seal(
    staging_dir: Path,
    *,
    attest_reviewed: bool,
    pic_root: Path | None = None,
    project_home_root: Path | None = None,
    runner: Runner = _default_runner,
    now: Clock = _utc_now,
    user: str | None = None,
    ledger_validator: LedgerValidator = _validate_with_existing_read_only_ledger_snapshot,
) -> Path:
    """Revalidate the selected execution boundary and publish one attestation."""
    if not attest_reviewed:
        raise ValueError("Seal requires explicit --attest-reviewed")
    staging = _canonical_existing_directory(staging_dir, label="Staging directory")
    archive_root = _canonical_existing_directory(staging.parent, label="Archive root")
    metadata = _read_json_regular(staging / METADATA_FILENAME)
    phase = _validate_metadata(metadata, staging)
    _require_exact_regular_files(staging, REVIEW_STAGING_FILENAMES)
    captured_user = _validate_user(str(metadata["same_account_user"]))
    if user is not None and _validate_user(user) != captured_user:
        raise ValueError("Same-account user override differs from staged capture")
    pic_root = _require_matching_optional_root(
        Path(str(metadata["pic_root"])), pic_root, label="PIC root"
    )
    project_home_root = _require_matching_optional_root(
        Path(str(metadata["project_home_root"])),
        project_home_root,
        label="Project Home root",
    )
    paths = _ledger_paths(pic_root, project_home_root)
    final = archive_root / str(metadata["final_directory_name"])
    _require_absent_entry(final, label="Final attestation directory")

    process_snapshot = _capture_same_account_processes(runner, captured_user)
    snapshot = _capture_gate_snapshot(
        paths,
        runner=runner,
        user=captured_user,
        ledger_validator=ledger_validator,
    )
    _require_sealable(snapshot)
    _atomic_replace_file(staging / "same_account_process_snapshot.txt", process_snapshot)
    for filename, payload in _gate_snapshot_files(snapshot).items():
        _atomic_replace_file(staging / filename, payload)
    sealed_now = now()
    _require_fresh_capture(metadata, sealed_now)
    _, sealed_utc = _timestamp_strings(sealed_now)
    _write_json_new(
        staging / "attestation.json",
        _attestation(staging, metadata, snapshot, phase=phase, sealed_utc=sealed_utc),
    )
    os.unlink(staging / METADATA_FILENAME)
    _require_exact_regular_files(staging, SEALED_FILENAMES)
    _fsync_directory(staging)
    _fsync_tree(staging)
    _make_tree_read_only(staging)
    _fsync_tree(staging)
    return publish_sealed_staging(staging)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    capture_parser = commands.add_parser(
        "capture", help="capture review inputs into hidden staging"
    )
    capture_parser.add_argument("--archive-root", type=Path, default=DEFAULT_ARCHIVE_ROOT)
    capture_parser.add_argument("--authorization-id", required=True)
    capture_parser.add_argument("--control-plane-version", required=True)
    capture_parser.add_argument("--phase", choices=PHASES, default=PHASE)
    capture_parser.add_argument("--pic-root", type=Path, default=AUTHORIZED_PIC_ROOT)
    capture_parser.add_argument(
        "--project-home-root", type=Path, default=AUTHORIZED_PROJECT_HOME_ROOT
    )
    seal_parser = commands.add_parser(
        "seal", help="revalidate gates and publish reviewed staging atomically"
    )
    seal_parser.add_argument("--staging-dir", type=Path, required=True)
    seal_parser.add_argument("--attest-reviewed", action="store_true")
    seal_parser.add_argument("--pic-root", type=Path)
    seal_parser.add_argument("--project-home-root", type=Path)
    recover_parser = commands.add_parser(
        "recover-sealed", help="publish an immutable staging tree after rename failure"
    )
    recover_parser.add_argument("--staging-dir", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _parser()
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "capture":
            path = capture(
                archive_root=arguments.archive_root,
                authorization_id=arguments.authorization_id,
                control_plane_version=arguments.control_plane_version,
                phase=arguments.phase,
                pic_root=arguments.pic_root,
                project_home_root=arguments.project_home_root,
            )
        elif arguments.command == "seal":
            path = seal(
                arguments.staging_dir,
                attest_reviewed=arguments.attest_reviewed,
                pic_root=arguments.pic_root,
                project_home_root=arguments.project_home_root,
                user=_current_user(),
            )
        else:
            path = publish_sealed_staging(arguments.staging_dir)
    except (OSError, ValueError) as error:
        parser.exit(1, f"error: {error}\n")
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
