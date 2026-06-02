#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Validate one sealed Frontier same-account isolation attestation."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import stat


RECORD_TYPE = "q027_frontier_registered_science_same_account_isolation_attestation"
CAPTURE_TO_SEAL_MAX_AGE_SECONDS = 15 * 60
PHASE_MAX_AGE_SECONDS = {
    "pre_policy_promotion": 60 * 60,
    "pre_manifest": 60 * 60,
    "pre_submit_wrapper": 15 * 60,
}
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
GATE_FILENAMES = {
    "queue_snapshot.txt",
    "pending_submission_marker.txt",
    "pending_manual_accounting_marker.txt",
    "mirrored_ledger_line_counts.txt",
    "validated_mirrored_ledger_state.json",
}
CAPTURED_SNAPSHOT_FILENAMES = {
    "capture_same_account_process_snapshot.txt",
    *(f"capture_{name}" for name in GATE_FILENAMES),
}
SEALED_FILENAMES = {
    "attestation.json",
    "same_account_process_snapshot.txt",
    *GATE_FILENAMES,
    *CAPTURED_SNAPSHOT_FILENAMES,
}
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
AUTHORIZATION_ID_PATTERN = re.compile(
    r"[A-Za-z0-9](?:[A-Za-z0-9_-]{0,126}[A-Za-z0-9])?"
)


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"Non-finite JSON number is not allowed: {value}")


def _reject_duplicate_json_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON object key is not allowed: {key}")
        value[key] = item
    return value


def _opened_identity(metadata: os.stat_result) -> tuple[int, int, int, int]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
    )


def _read_regular_bytes_at(
    directory_fd: int,
    filename: str,
    *,
    expected_mode: int = 0o400,
    retained_files: dict[str, tuple[int, tuple[int, int, int, int]]] | None = None,
) -> bytes:
    if not filename or filename in {".", ".."} or "/" in filename:
        raise ValueError(f"Sealed attestation member name is invalid: {filename!r}")
    descriptor = (
        retained_files[filename][0]
        if retained_files is not None
        else os.open(
            filename,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_fd,
        )
    )
    close_descriptor = retained_files is None
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or stat.S_IMODE(before.st_mode) != expected_mode
        ):
            raise ValueError(f"Sealed attestation member mode is invalid: {filename}")
        chunks: list[bytes] = []
        offset = 0
        while True:
            chunk = os.pread(descriptor, 1024 * 1024, offset)
            if not chunk:
                break
            chunks.append(chunk)
            offset += len(chunk)
        after = os.fstat(descriptor)
        if _opened_identity(before) != _opened_identity(after):
            raise ValueError(f"Sealed attestation member changed while read: {filename}")
        return b"".join(chunks)
    finally:
        if close_descriptor:
            os.close(descriptor)


def _read_json_object_at(
    directory_fd: int,
    filename: str,
    *,
    retained_files: dict[str, tuple[int, tuple[int, int, int, int]]] | None = None,
) -> tuple[bytes, dict[str, object]]:
    payload = _read_regular_bytes_at(
        directory_fd, filename, retained_files=retained_files
    )
    try:
        value = json.loads(
            payload.decode("utf-8"),
            parse_constant=_reject_json_constant,
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"Sealed attestation JSON is invalid: {filename}") from error
    if not isinstance(value, dict):
        raise ValueError(f"Sealed attestation JSON is not an object: {filename}")
    return payload, value


def _timestamp(value: object, *, label: str) -> datetime:
    if not isinstance(value, str):
        raise ValueError(f"{label} must be a canonical UTC timestamp")
    try:
        return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(
            tzinfo=timezone.utc
        )
    except ValueError as error:
        raise ValueError(f"{label} must be a canonical UTC timestamp") from error


def _file_record(
    root_fd: int,
    record: object,
    *,
    filename: str,
    extra_keys: set[str] | None = None,
    retained_files: dict[str, tuple[int, tuple[int, int, int, int]]] | None = None,
) -> tuple[bytes, dict[str, object]]:
    expected_extra_keys = extra_keys or set()
    if (
        not isinstance(record, dict)
        or set(record) != {"path", "sha256", *expected_extra_keys}
        or record.get("path") != filename
        or not isinstance(record.get("sha256"), str)
        or SHA256_PATTERN.fullmatch(str(record["sha256"])) is None
    ):
        raise ValueError(f"Sealed attestation file record is invalid: {filename}")
    payload = _read_regular_bytes_at(
        root_fd, filename, retained_files=retained_files
    )
    if hashlib.sha256(payload).hexdigest() != record["sha256"]:
        raise ValueError(f"Sealed attestation member checksum differs: {filename}")
    return payload, record


def _require_same_directory(path: Path, descriptor: int, *, label: str) -> None:
    metadata = os.stat(path, follow_symlinks=False)
    retained = os.fstat(descriptor)
    if (
        not stat.S_ISDIR(metadata.st_mode)
        or (metadata.st_dev, metadata.st_ino) != (retained.st_dev, retained.st_ino)
    ):
        raise ValueError(f"{label} changed while the attestation was validated")


def _expected_manual_marker_paths(
    pic_root: Path, project_home_root: Path
) -> set[str]:
    return {
        str(pic_root / "ledger" / "pending_manual_accounting.json"),
        str(project_home_root / "ledger" / "pending_manual_accounting.json"),
    }


def _expected_ledger_paths(pic_root: Path, project_home_root: Path) -> set[str]:
    return {
        str(pic_root / "ledger" / "node_hours.jsonl"),
        str(pic_root / "ledger" / "mirror_receipts.jsonl"),
        str(project_home_root / "ledger" / "node_hours.jsonl"),
    }


def _validate_captured_snapshots(
    root_fd: int,
    records: object,
    *,
    retained_files: dict[str, tuple[int, tuple[int, int, int, int]]],
) -> None:
    if not isinstance(records, dict) or set(records) != CAPTURED_SNAPSHOT_FILENAMES:
        raise ValueError("Captured snapshot binding set differs")
    for filename, record in records.items():
        _file_record(
            root_fd, record, filename=filename, retained_files=retained_files
        )


def _require_same_regular_files(
    directory_fd: int,
    retained_files: dict[str, tuple[int, tuple[int, int, int, int]]],
) -> None:
    for filename, (descriptor, initial) in retained_files.items():
        retained = os.fstat(descriptor)
        current = os.stat(filename, dir_fd=directory_fd, follow_symlinks=False)
        if (
            not stat.S_ISREG(current.st_mode)
            or stat.S_IMODE(current.st_mode) != 0o400
            or _opened_identity(retained) != initial
            or _opened_identity(current) != initial
        ):
            raise ValueError(
                f"Sealed attestation member changed while validated: {filename}"
            )


def validate_sealed_operator_attestation(
    attestation_path: Path,
    *,
    authorization_id: str,
    phase: str,
    control_plane_version: str,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    now: datetime | None = None,
    enforce_freshness: bool = True,
) -> dict[str, str]:
    """Return an immutable path/digest binding after descriptor-rooted validation."""
    if AUTHORIZATION_ID_PATTERN.fullmatch(authorization_id) is None:
        raise ValueError("Operator-attestation authorization ID is malformed")
    if phase not in PHASE_MAX_AGE_SECONDS:
        raise ValueError("Operator-attestation phase is unsupported")
    if SHA256_PATTERN.fullmatch(control_plane_version) is None:
        raise ValueError("Operator-attestation control-plane version is malformed")
    pic_root = Path(os.path.abspath(authorized_pic_root))
    project_home_root = Path(os.path.abspath(authorized_project_home_root))
    archive_root = pic_root / "operator_attestations"
    path = Path(os.path.abspath(attestation_path))
    if path.name != "attestation.json" or path.parent.parent != archive_root:
        raise ValueError("Operator-attestation path is outside the fixed archive layout")
    if pic_root.resolve(strict=True) != pic_root:
        raise ValueError("Operator-attestation PIC root must not use a symlink alias")

    flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    pic_fd = os.open(pic_root, flags)
    archive_fd = -1
    root_fd = -1
    retained_files: dict[str, tuple[int, tuple[int, int, int, int]]] = {}
    try:
        archive_fd = os.open("operator_attestations", flags, dir_fd=pic_fd)
        root_fd = os.open(path.parent.name, flags, dir_fd=archive_fd)
        root_metadata = os.fstat(root_fd)
        if stat.S_IMODE(root_metadata.st_mode) != 0o500:
            raise ValueError("Sealed operator-attestation directory mode is invalid")
        actual_names = set(os.listdir(root_fd))
        if actual_names != SEALED_FILENAMES:
            raise ValueError("Sealed operator-attestation tree closure differs")
        for name in actual_names:
            descriptor = os.open(
                name,
                os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=root_fd,
            )
            metadata = os.fstat(descriptor)
            if (
                not stat.S_ISREG(metadata.st_mode)
                or stat.S_IMODE(metadata.st_mode) != 0o400
            ):
                os.close(descriptor)
                raise ValueError(f"Sealed attestation member mode is invalid: {name}")
            retained_files[name] = (descriptor, _opened_identity(metadata))

        payload, attestation = _read_json_object_at(
            root_fd, "attestation.json", retained_files=retained_files
        )
        required = {
            "schema_version",
            "record_type",
            "recorded_utc",
            "sealed_utc",
            "registered_science_authorization_id",
            "control_plane_version",
            "phase",
            "same_account_process_snapshot",
            "queue_snapshot",
            "pending_submission_marker",
            "pending_manual_accounting_marker",
            "mirrored_ledger_line_counts",
            "validated_mirrored_ledger_state",
            "captured_snapshots",
            "operator_statement",
        }
        if set(attestation) != required:
            raise ValueError("Sealed operator-attestation root schema is invalid")
        if (
            type(attestation.get("schema_version")) is not int
            or attestation["schema_version"] != 1
            or attestation.get("record_type") != RECORD_TYPE
            or attestation.get("registered_science_authorization_id") != authorization_id
            or attestation.get("control_plane_version") != control_plane_version
            or attestation.get("phase") != phase
            or attestation.get("operator_statement") != OPERATOR_STATEMENTS[phase]
        ):
            raise ValueError("Sealed operator-attestation identity differs")
        recorded = _timestamp(attestation["recorded_utc"], label="recorded_utc")
        sealed = _timestamp(attestation["sealed_utc"], label="sealed_utc")
        current = (now or datetime.now(timezone.utc)).astimezone(timezone.utc)
        if not recorded <= sealed <= current:
            raise ValueError("Sealed operator-attestation timestamps are out of order")
        if (sealed - recorded).total_seconds() > CAPTURE_TO_SEAL_MAX_AGE_SECONDS:
            raise ValueError("Sealed operator-attestation capture-to-seal interval is stale")
        if (
            enforce_freshness
            and (current - sealed).total_seconds() > PHASE_MAX_AGE_SECONDS[phase]
        ):
            raise ValueError("Sealed operator-attestation is stale")
        expected_name = (
            f"{recorded.strftime('%Y%m%dT%H%M%SZ')}-{authorization_id}-{phase}"
        )
        if path.parent.name != expected_name:
            raise ValueError("Sealed operator-attestation directory name differs")

        process_snapshot, _ = _file_record(
            root_fd,
            attestation["same_account_process_snapshot"],
            filename="same_account_process_snapshot.txt",
            retained_files=retained_files,
        )
        if not process_snapshot:
            raise ValueError("Sealed operator-attestation process snapshot is empty")
        queue, _ = _file_record(
            root_fd,
            attestation["queue_snapshot"],
            filename="queue_snapshot.txt",
            retained_files=retained_files,
        )
        if queue.strip():
            raise ValueError("Sealed operator-attestation queue snapshot is not empty")
        pending, pending_record = _file_record(
            root_fd,
            attestation["pending_submission_marker"],
            filename="pending_submission_marker.txt",
            extra_keys={"value"},
            retained_files=retained_files,
        )
        if pending != b"absent\n" or pending_record["value"] != "absent":
            raise ValueError("Sealed operator-attestation pending submission marker differs")
        manual, manual_record = _file_record(
            root_fd,
            attestation["pending_manual_accounting_marker"],
            filename="pending_manual_accounting_marker.txt",
            extra_keys={"value", "values"},
            retained_files=retained_files,
        )
        expected_markers = _expected_manual_marker_paths(pic_root, project_home_root)
        if (
            manual_record["value"] != "absent"
            or not isinstance(manual_record["values"], dict)
            or set(manual_record["values"]) != expected_markers
            or any(value != "absent" for value in manual_record["values"].values())
            or sorted(manual.decode("utf-8").splitlines())
            != sorted(f"{name} absent" for name in expected_markers)
        ):
            raise ValueError("Sealed operator-attestation manual-accounting marker differs")
        counts_payload, counts_record = _file_record(
            root_fd,
            attestation["mirrored_ledger_line_counts"],
            filename="mirrored_ledger_line_counts.txt",
            extra_keys={"counts"},
            retained_files=retained_files,
        )
        expected_ledgers = _expected_ledger_paths(pic_root, project_home_root)
        counts = counts_record["counts"]
        if (
            not isinstance(counts, dict)
            or set(counts) != expected_ledgers
            or any(type(value) is not int or value < 0 for value in counts.values())
            or sorted(counts_payload.decode("utf-8").splitlines())
            != sorted(f"{value} {name}" for name, value in counts.items())
            or len(set(counts.values())) != 1
        ):
            raise ValueError("Sealed operator-attestation ledger line counts differ")
        state_payload, state_record = _file_record(
            root_fd,
            attestation["validated_mirrored_ledger_state"],
            filename="validated_mirrored_ledger_state.json",
            extra_keys={"state"},
            retained_files=retained_files,
        )
        try:
            state_from_file = json.loads(
                state_payload.decode("utf-8"),
                parse_constant=_reject_json_constant,
                object_pairs_hook=_reject_duplicate_json_pairs,
            )
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise ValueError("Sealed operator-attestation ledger state is invalid") from error
        state = state_record["state"]
        expected_state_keys = {
            "validation",
            "validator",
            "ledger_record_count",
            "ledger_tail_event_sha256",
            "active_reservation_count",
            "active_reservation_ids",
        }
        tail = state.get("ledger_tail_event_sha256") if isinstance(state, dict) else None
        if (
            not isinstance(state, dict)
            or set(state) != expected_state_keys
            or state_from_file != state
            or state.get("validation") != "coherent"
            or state.get("validator") != "existing_read_only_mirrored_state_snapshot"
            or state.get("ledger_record_count") != next(iter(counts.values()))
            or (tail is not None and (not isinstance(tail, str) or SHA256_PATTERN.fullmatch(tail) is None))
            or state.get("active_reservation_count") != 0
            or state.get("active_reservation_ids") != []
        ):
            raise ValueError("Sealed operator-attestation ledger state differs")
        _validate_captured_snapshots(
            root_fd, attestation["captured_snapshots"], retained_files=retained_files
        )

        _require_same_regular_files(root_fd, retained_files)
        _require_same_directory(path.parent, root_fd, label="Operator-attestation directory")
        _require_same_directory(archive_root, archive_fd, label="Operator-attestation archive")
        _require_same_directory(pic_root, pic_fd, label="Operator-attestation PIC root")
        return {"path": str(path), "sha256": hashlib.sha256(payload).hexdigest()}
    finally:
        for descriptor, _ in retained_files.values():
            os.close(descriptor)
        if root_fd >= 0:
            os.close(root_fd)
        if archive_fd >= 0:
            os.close(archive_fd)
        os.close(pic_fd)
