#!/usr/bin/env python3
"""Publish the reviewed Q-011 Section 5.4 pressure selection without launch authority.

The schema-v3 pressure-selection receipt remains unchanged.  External
acceptance additionally requires a selection-specific success seal that binds
the receipt inode and one immutable, launch-prohibited controller-state
attestation captured while policy mutation is serialized.
"""

from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.abc
import importlib.util
import io
import json
import os
from pathlib import Path
import re
import stat
import sys
import tarfile
import tempfile
from types import ModuleType
from typing import Any, Mapping, Sequence
import uuid

_PUBLICATION_DIR = Path(__file__).resolve().parent
_CONTROL_PLANE_DIR = _PUBLICATION_DIR / "frontier_control_plane"
if str(_CONTROL_PLANE_DIR) not in sys.path:
    sys.path.insert(0, str(_CONTROL_PLANE_DIR))

if __package__:
    from . import publish_q011_section54_pressure_pilot_bundle as pilot_publisher
    from . import q011_section54_pressure_selection as pressure_selection
    from .frontier_control_plane import control_plane_common
    from .frontier_control_plane import ledger
    from .frontier_control_plane import (
        revalidate_clean_candidate as clean_candidate_revalidator,
    )
else:
    sys.path.insert(0, str(_PUBLICATION_DIR))
    import publish_q011_section54_pressure_pilot_bundle as pilot_publisher
    import q011_section54_pressure_selection as pressure_selection
    from frontier_control_plane import control_plane_common
    from frontier_control_plane import ledger
    from frontier_control_plane import (
        revalidate_clean_candidate as clean_candidate_revalidator,
    )


AUTHORIZED_PIC_ROOT = control_plane_common.AUTHORIZED_PIC_ROOT
AUTHORIZED_PROJECT_HOME_ROOT = control_plane_common.AUTHORIZED_PROJECT_HOME_ROOT
CANONICAL_RECEIPT_NAME = "q011_section54_pressure_selection_receipt.json"
PRESSURE_GATE_ATTESTATION_ROOT_NAME = "pressure_gate_attestations"
PRESSURE_GATE_CANDIDATE_ROOT_NAME = "pressure_gate_candidates"
PRESSURE_GATE_HUMAN_DECISION_ROOT_NAME = "pressure_gate_human_decisions"
HUMAN_DECISION_RECORD_TYPE = "q011_section54_pressure_selection_human_decision"
HUMAN_DECISION_STATEMENT = (
    "I reviewed the immutable pressure-pilot packet and the exact sealed "
    "authoritative reanalysis, and I select the bound pressure case."
)
CONTROLLER_STATE_RECORD_TYPE = (
    "q011_section54_pressure_selection_controller_state_attestation"
)
SUCCESS_SEAL_RECORD_TYPE = (
    "q011_section54_pressure_selection_publication_success_seal"
)
RECOVERY_GUARD_RECORD_TYPE = (
    "q011_section54_pressure_selection_publication_recovery_guard"
)
STAGE4_PREPARATION_RECORD_TYPE = (
    "q011_section54_pressure_selection_stage4_preparation_attestation"
)
CANDIDATE_AUTHORIZATION_RECORD_TYPE = (
    "q011_section54_pressure_selection_candidate_publication_authorization"
)
QUALIFICATION_EFFECT = (
    "human_pressure_selection_publication_only_no_science_launch_authority"
)
REVIEWED_SELECTED_CASE = {
    "case_id": "ps_p0_1p00",
    "problem_ps_p0": 1.0,
}
REVIEWED_REVIEWER_ID = "dfielding"
REVIEWED_RATIONALE = (
    "Selected p0=1.0 as the recommended baseline because it explicitly matches "
    "Bai et al. (2015), which uses P0=T0=1 and treats the choice as unimportant "
    "while thermal pressure is much smaller than ram pressure."
)
MAX_JSON_BYTES = 8 * 1024 * 1024
MAX_CANDIDATE_ARCHIVE_RAW_BYTES = 256 * 1024 * 1024
MAX_CANDIDATE_ARCHIVE_MEMBERS = 4096
MAX_CANDIDATE_ARCHIVE_REGULAR_FILES = 4096
MAX_CANDIDATE_ARCHIVE_REGULAR_BYTES = 128 * 1024 * 1024
MAX_CANDIDATE_ARCHIVE_MEMBER_BYTES = 16 * 1024 * 1024
MAX_CANDIDATE_ARCHIVE_MEMBER_NAME_CHARACTERS = 4096
MAX_CANDIDATE_ARCHIVE_PATH_DEPTH = 16
MAX_CANDIDATE_TREE_ARCHIVES = 64
MAX_CANDIDATE_TREE_ARCHIVE_MEMBERS = 16384
MAX_CANDIDATE_TREE_ARCHIVE_REGULAR_BYTES = 512 * 1024 * 1024
MAX_CANDIDATE_TREE_REGULAR_BYTES = 512 * 1024 * 1024
MAX_CANDIDATE_TREE_MEMBER_BYTES = 256 * 1024 * 1024
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_FILE_FLAGS = (
    os.O_RDONLY
    | getattr(os, "O_NOFOLLOW", 0)
    | getattr(os, "O_NONBLOCK", 0)
)
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT_PATTERN = re.compile(r"[0-9a-f]{40}")
_OPERATOR_ID_PATTERN = re.compile(r"[a-z][a-z0-9._-]{2,63}")
PUBLISHER_SOURCE_PATHS = tuple(
    sorted(
        {
            "tst/publication/publish_q011_section54_pressure_selection.py",
            *control_plane_common.Q011_SECTION54_HELPER_SOURCES,
            *(
                f"tst/publication/frontier_control_plane/{name}"
                for name in control_plane_common.CONTROL_PLANE_FILES
            ),
        }
    )
)


class PressureSelectionPublicationError(ValueError):
    """Raised when pressure-selection publication fails closed."""

    def __init__(
        self,
        message: str,
        *,
        recovery: Mapping[str, object] | None = None,
    ) -> None:
        super().__init__(message)
        self.recovery = None if recovery is None else copy.deepcopy(dict(recovery))


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PressureSelectionPublicationError(message)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical_json_bytes(value: object) -> bytes:
    try:
        payload = (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (RecursionError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "pressure-selection publication value is not canonical JSON"
        ) from error
    _require(
        len(payload) <= MAX_JSON_BYTES,
        "pressure-selection publication JSON exceeds the size limit",
    )
    return payload


def _decode_canonical_object(payload: bytes, label: str) -> dict[str, Any]:
    def reject_constant(value: str) -> None:
        raise PressureSelectionPublicationError(
            f"{label} contains forbidden JSON constant: {value}"
        )

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            _require(key not in result, f"{label} contains duplicate JSON key: {key}")
            result[key] = value
        return result

    try:
        value = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as error:
        raise PressureSelectionPublicationError(
            f"{label} is not valid UTF-8 JSON"
        ) from error
    _require(type(value) is dict, f"{label} is not a JSON object")
    _require(payload == _canonical_json_bytes(value), f"{label} is not canonical JSON")
    return value


def _canonical_utc(now: datetime | None = None) -> str:
    current = datetime.now(timezone.utc) if now is None else now
    _require(
        current.tzinfo is not None and current.utcoffset() is not None,
        "pressure-selection publication timestamp must be timezone-aware",
    )
    return current.astimezone(timezone.utc).replace(microsecond=0).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )


def _require_canonical_utc(value: object, label: str) -> str:
    _require(type(value) is str, f"{label} must be canonical UTC text")
    try:
        parsed = datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(
            tzinfo=timezone.utc
        )
    except ValueError as error:
        raise PressureSelectionPublicationError(
            f"{label} must be canonical UTC text"
        ) from error
    _require(_canonical_utc(parsed) == value, f"{label} must be canonical UTC text")
    return value


def _operator_id(value: object, label: str) -> str:
    _require(
        type(value) is str and _OPERATOR_ID_PATTERN.fullmatch(value) is not None,
        f"{label} must be a canonical lowercase operator ID",
    )
    return value


def _lowercase_sha256(value: object, label: str) -> str:
    _require(
        type(value) is str and _SHA256_PATTERN.fullmatch(value) is not None,
        f"{label} must be a lowercase SHA-256",
    )
    return value


def _binding(value: object, label: str) -> dict[str, str]:
    _require(
        type(value) is dict and set(value) == {"path", "sha256"},
        f"{label} binding schema drifted",
    )
    path = value["path"]
    _require(
        type(path) is str and Path(path).is_absolute(),
        f"{label} path must be absolute",
    )
    return {
        "path": path,
        "sha256": _lowercase_sha256(value["sha256"], f"{label} SHA-256"),
    }


def _trusted_publisher_source_archive(expected_git_commit: str) -> tuple[str, bytes]:
    try:
        repository = pilot_publisher.TRUSTED_SOURCE_REPOSITORY.resolve(strict=True)
    except OSError as error:
        raise PressureSelectionPublicationError(
            "trusted Stage-4 source repository is unavailable"
        ) from error
    _require(repository.is_dir(), "trusted Stage-4 source repository is invalid")
    git_environment = control_plane_common.trusted_git_environment()
    git_environment["GIT_NO_REPLACE_OBJECTS"] = "1"
    commit_result = pilot_publisher.subprocess.run(
        control_plane_common.trusted_git_command(
            "--no-replace-objects",
            "-C",
            str(repository),
            "rev-parse",
            "--verify",
            f"{expected_git_commit}^{{commit}}",
        ),
        check=False,
        stdout=pilot_publisher.subprocess.PIPE,
        stderr=pilot_publisher.subprocess.PIPE,
        env=git_environment,
    )
    _require(
        commit_result.returncode == 0,
        "trusted Stage-4 source commit is unavailable",
    )
    commit = commit_result.stdout.decode("ascii").strip()
    _require(
        commit == expected_git_commit
        and _GIT_COMMIT_PATTERN.fullmatch(commit) is not None,
        "trusted Stage-4 source commit differs from the expected commit",
    )
    archive_result = pilot_publisher.subprocess.run(
        control_plane_common.trusted_git_command(
            "--no-replace-objects",
            "-C",
            str(repository),
            "archive",
            "--format=tar",
            commit,
        ),
        check=False,
        stdout=pilot_publisher.subprocess.PIPE,
        stderr=pilot_publisher.subprocess.PIPE,
        env=git_environment,
    )
    _require(
        archive_result.returncode == 0 and bool(archive_result.stdout),
        "trusted Stage-4 source archive is unavailable",
    )
    return commit, archive_result.stdout


def _runtime_publisher_source_authentication(expected_git_commit: str) -> dict[str, object]:
    expected = str(expected_git_commit)
    _require(
        _GIT_COMMIT_PATTERN.fullmatch(expected) is not None,
        "expected publisher Git commit is malformed",
    )
    archive_path_text = os.environ.get("PIC_PRESSURE_PUBLICATION_SOURCE_ARCHIVE_PATH")
    snapshot_root_text = os.environ.get(
        pilot_publisher.WORKER_SOURCE_SNAPSHOT_ROOT_ENV
    )
    _require(
        archive_path_text is not None and snapshot_root_text is not None,
        "pressure-selection mutation requires an authenticated source snapshot",
    )
    archive_path = Path(archive_path_text)
    snapshot_root = Path(snapshot_root_text)
    _require(
        archive_path.is_absolute() and snapshot_root.is_absolute(),
        "pressure-selection authenticated source paths must be absolute",
    )
    try:
        snapshot_root = snapshot_root.resolve(strict=True)
        executing_root = Path(__file__).resolve().parents[2]
        trusted_repository = pilot_publisher.TRUSTED_SOURCE_REPOSITORY.resolve(
            strict=True
        )
    except OSError as error:
        raise PressureSelectionPublicationError(
            "pressure-selection authenticated source snapshot is unavailable"
        ) from error
    _require(
        snapshot_root == executing_root and snapshot_root != trusted_repository,
        "pressure-selection mutation is not executing from its authenticated snapshot",
    )
    try:
        commit, trusted_archive = _trusted_publisher_source_archive(expected)
        archive_payload = pilot_publisher._read_stable_readonly_regular(
            archive_path,
            "pressure-selection authenticated source archive",
            max_bytes=len(trusted_archive),
        )
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "pressure-selection authenticated source archive failed verification"
        ) from error
    _require(
        commit == expected and archive_payload == trusted_archive,
        "pressure-selection authenticated source archive differs from the expected commit",
    )
    try:
        archive_files = control_plane_common._source_archive_regular_files(
            trusted_archive
        )
    except (OSError, RecursionError, tarfile.TarError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "pressure-selection authenticated source archive failed verification"
        ) from error
    closure = []
    for relative in PUBLISHER_SOURCE_PATHS:
        archived = archive_files.get(relative)
        _require(
            archived is not None,
            f"pressure-selection authenticated source archive omits {relative}",
        )
        try:
            snapshot = pilot_publisher._read_stable_readonly_regular(
                snapshot_root / relative,
                f"pressure-selection authenticated source member {relative}",
            )
        except (OSError, TypeError, ValueError) as error:
            raise PressureSelectionPublicationError(
                f"pressure-selection authenticated source member is unavailable: {relative}"
            ) from error
        _require(
            snapshot == archived,
            f"pressure-selection authenticated source member drifted: {relative}",
        )
        closure.append({"path": relative, "sha256": _sha256(archived)})
    return {
        "execution_mode": "worker_extracted_git_archive_expected_commit_verified",
        "git_commit": expected,
        "archive_sha256": _sha256(archive_payload),
        "source_closure_sha256": _sha256(
            json.dumps(
                closure,
                separators=(",", ":"),
                sort_keys=True,
                allow_nan=False,
            ).encode("utf-8")
        ),
        "source_closure": closure,
    }


def _validate_publisher_source_authentication(value: object) -> dict[str, object]:
    _require(
        type(value) is dict
        and set(value)
        == {
            "execution_mode",
            "git_commit",
            "archive_sha256",
            "source_closure_sha256",
            "source_closure",
        }
        and value["execution_mode"]
        == "worker_extracted_git_archive_expected_commit_verified"
        and type(value["git_commit"]) is str
        and _GIT_COMMIT_PATTERN.fullmatch(value["git_commit"]) is not None,
        "publisher source authentication schema drifted",
    )
    _lowercase_sha256(value["archive_sha256"], "publisher source archive SHA-256")
    closure_sha256 = _lowercase_sha256(
        value["source_closure_sha256"], "publisher source closure SHA-256"
    )
    closure = value["source_closure"]
    _require(
        type(closure) is list
        and [record.get("path") for record in closure if type(record) is dict]
        == list(PUBLISHER_SOURCE_PATHS),
        "publisher source closure path set drifted",
    )
    for record in closure:
        _require(
            type(record) is dict
            and set(record) == {"path", "sha256"}
            and type(record["path"]) is str,
            "publisher source closure record drifted",
        )
        _lowercase_sha256(record["sha256"], "publisher source member SHA-256")
    _require(
        closure_sha256
        == _sha256(
            json.dumps(
                closure,
                separators=(",", ":"),
                sort_keys=True,
                allow_nan=False,
            ).encode("utf-8")
        ),
        "publisher source closure SHA-256 drifted",
    )
    return copy.deepcopy(value)


def _read_bounded_descriptor(descriptor: int, label: str) -> bytes:
    metadata = os.fstat(descriptor)
    _require(metadata.st_size <= MAX_JSON_BYTES, f"{label} exceeds the size limit")
    payload = bytearray()
    while True:
        chunk = os.read(descriptor, min(1024 * 1024, MAX_JSON_BYTES + 1 - len(payload)))
        if not chunk:
            return bytes(payload)
        payload.extend(chunk)
        _require(len(payload) <= MAX_JSON_BYTES, f"{label} exceeds the size limit")


def _read_exact_mode_regular_at(
    parent_descriptor: int,
    name: str,
    *,
    mode: int,
    label: str,
) -> tuple[bytes, tuple[int, int]]:
    descriptor = os.open(name, _FILE_FLAGS, dir_fd=parent_descriptor)
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and stat.S_IMODE(before.st_mode) == mode
            and before.st_nlink == 1,
            f"{label} must be one singly linked regular file with mode {mode:04o}",
        )
        payload = _read_bounded_descriptor(descriptor, label)
        after = os.fstat(descriptor)
        current = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
        stable = (
            before.st_dev,
            before.st_ino,
            before.st_mode,
            before.st_nlink,
            before.st_size,
            before.st_mtime_ns,
            before.st_ctime_ns,
        )
        _require(
            stable
            == (
                after.st_dev,
                after.st_ino,
                after.st_mode,
                after.st_nlink,
                after.st_size,
                after.st_mtime_ns,
                after.st_ctime_ns,
            )
            == (
                current.st_dev,
                current.st_ino,
                current.st_mode,
                current.st_nlink,
                current.st_size,
                current.st_mtime_ns,
                current.st_ctime_ns,
            )
            and len(payload) == after.st_size,
            f"{label} changed while reading",
        )
        return payload, (after.st_dev, after.st_ino)
    finally:
        os.close(descriptor)


def _write_controller_member_at(parent_descriptor: int, name: str, payload: bytes) -> None:
    _require("/" not in name and name not in {"", ".", ".."}, "invalid attestation member")
    descriptor = os.open(
        name,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o600,
        dir_fd=parent_descriptor,
    )
    try:
        with os.fdopen(descriptor, "wb", closefd=False) as stream:
            stream.write(payload)
            stream.flush()
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o400)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _unlink_exact_file_at(
    parent_descriptor: int,
    name: str,
    identity: tuple[int, int] | None,
) -> None:
    if identity is None:
        return
    try:
        pilot_publisher._require_same_file_at(
            parent_descriptor, name, identity, "pressure-selection owned staging file"
        )
    except FileNotFoundError:
        return
    os.unlink(name, dir_fd=parent_descriptor)
    pilot_publisher._fsync_descriptor(parent_descriptor)


def _selection_success_seal_name(receipt_name: str) -> str:
    return pilot_publisher._publication_seal_name(receipt_name)


def _selection_recovery_guard_payload(
    *,
    receipt_path: Path,
    receipt_sha256: str,
    controller_state_attestation: object,
    expected_git_commit: str,
) -> bytes:
    return _canonical_json_bytes(
        {
            "schema_version": 1,
            "record_type": RECOVERY_GUARD_RECORD_TYPE,
            "status": "publication_invalid_pending_authenticated_reconciliation",
            "qualification_effect": QUALIFICATION_EFFECT,
            "expected_git_commit": expected_git_commit,
            "receipt": {
                "path": str(receipt_path),
                "sha256": _lowercase_sha256(
                    receipt_sha256, "recovery guard receipt SHA-256"
                ),
            },
            "controller_state_attestation": _binding(
                controller_state_attestation, "recovery guard controller-state attestation"
            ),
        }
    )


def _validate_selection_recovery_guard(
    value: object,
    *,
    receipt_path: Path,
    receipt_sha256: str | None = None,
    controller_state_attestation: object | None = None,
    expected_git_commit: str | None = None,
) -> dict[str, object]:
    _require(
        type(value) is dict
        and set(value)
        == {
            "schema_version",
            "record_type",
            "status",
            "qualification_effect",
            "expected_git_commit",
            "receipt",
            "controller_state_attestation",
        },
        "pressure-selection recovery guard schema drifted",
    )
    receipt = value["receipt"]
    _require(
        type(value["schema_version"]) is int
        and value["schema_version"] == 1
        and value["record_type"] == RECOVERY_GUARD_RECORD_TYPE
        and value["status"]
        == "publication_invalid_pending_authenticated_reconciliation"
        and value["qualification_effect"] == QUALIFICATION_EFFECT
        and type(value["expected_git_commit"]) is str
        and _GIT_COMMIT_PATTERN.fullmatch(value["expected_git_commit"]) is not None
        and type(receipt) is dict
        and set(receipt) == {"path", "sha256"}
        and receipt["path"] == str(receipt_path),
        "pressure-selection recovery guard binding drifted",
    )
    _lowercase_sha256(receipt["sha256"], "recovery guard receipt SHA-256")
    normalized = copy.deepcopy(value)
    normalized["controller_state_attestation"] = _binding(
        value["controller_state_attestation"],
        "recovery guard controller-state attestation",
    )
    if receipt_sha256 is not None:
        _require(
            receipt["sha256"] == receipt_sha256,
            "pressure-selection recovery guard binds another receipt",
        )
    if controller_state_attestation is not None:
        _require(
            normalized["controller_state_attestation"]
            == _binding(
                controller_state_attestation,
                "expected recovery guard controller-state attestation",
            ),
            "pressure-selection recovery guard binds another controller state",
        )
    if expected_git_commit is not None:
        _require(
            value["expected_git_commit"] == expected_git_commit,
            "pressure-selection recovery guard binds another Git commit",
        )
    return normalized


def _read_selection_recovery_guard_at(
    publication_descriptor: int,
    receipt_path: Path,
    *,
    receipt_sha256: str | None = None,
    controller_state_attestation: object | None = None,
    expected_git_commit: str | None = None,
) -> tuple[dict[str, object], tuple[int, int]]:
    payload, identity = _read_exact_mode_regular_at(
        publication_descriptor,
        pilot_publisher._publication_guard_name(receipt_path.name),
        mode=0o444,
        label="pressure-selection recovery guard",
    )
    value = _decode_canonical_object(payload, "pressure-selection recovery guard")
    return (
        _validate_selection_recovery_guard(
            value,
            receipt_path=receipt_path,
            receipt_sha256=receipt_sha256,
            controller_state_attestation=controller_state_attestation,
            expected_git_commit=expected_git_commit,
        ),
        identity,
    )


def _arm_selection_recovery_guard_at(
    publication_descriptor: int,
    receipt_path: Path,
    *,
    receipt_sha256: str,
    controller_state_attestation: object,
    expected_git_commit: str,
) -> tuple[dict[str, object], tuple[int, int]]:
    name = pilot_publisher._publication_guard_name(receipt_path.name)
    pilot_publisher._require_absent_at(
        publication_descriptor, name, "pressure-selection recovery guard"
    )
    payload = _selection_recovery_guard_payload(
        receipt_path=receipt_path,
        receipt_sha256=receipt_sha256,
        controller_state_attestation=controller_state_attestation,
        expected_git_commit=expected_git_commit,
    )
    pilot_publisher._write_exclusive_at(publication_descriptor, name, payload)
    pilot_publisher._fsync_descriptor(publication_descriptor)
    return _read_selection_recovery_guard_at(
        publication_descriptor,
        receipt_path,
        receipt_sha256=receipt_sha256,
        controller_state_attestation=controller_state_attestation,
        expected_git_commit=expected_git_commit,
    )


def _ensure_selection_recovery_guard_at(
    publication_descriptor: int,
    receipt_path: Path,
    *,
    receipt_sha256: str,
    controller_state_attestation: object,
    expected_git_commit: str,
) -> tuple[dict[str, object], tuple[int, int]]:
    try:
        return _read_selection_recovery_guard_at(
            publication_descriptor,
            receipt_path,
            receipt_sha256=receipt_sha256,
            controller_state_attestation=controller_state_attestation,
            expected_git_commit=expected_git_commit,
        )
    except FileNotFoundError:
        return _arm_selection_recovery_guard_at(
            publication_descriptor,
            receipt_path,
            receipt_sha256=receipt_sha256,
            controller_state_attestation=controller_state_attestation,
            expected_git_commit=expected_git_commit,
        )


def _selection_success_seal_payload(
    publication_descriptor: int,
    receipt_name: str,
    receipt_payload: bytes,
    receipt_identity: tuple[int, int],
    controller_state_attestation: object,
    *,
    sealed_utc: str,
) -> bytes:
    controller_binding = _binding(
        controller_state_attestation, "controller-state attestation"
    )
    return _canonical_json_bytes(
        {
            "schema_version": 1,
            "record_type": SUCCESS_SEAL_RECORD_TYPE,
            "qualification_effect": QUALIFICATION_EFFECT,
            "sealed_utc": _require_canonical_utc(
                sealed_utc, "selection success seal sealed_utc"
            ),
            "publication_root_identity": pilot_publisher._directory_identity(
                publication_descriptor
            ),
            "receipt_name": receipt_name,
            "receipt_sha256": _sha256(receipt_payload),
            "receipt_identity": {
                "device": receipt_identity[0],
                "inode": receipt_identity[1],
            },
            "controller_state_attestation": controller_binding,
        }
    )


def _validate_selection_success_seal(
    value: object,
    publication_descriptor: int,
    receipt_name: str,
    receipt_payload: bytes,
    receipt_identity: tuple[int, int],
) -> dict[str, object]:
    expected_keys = {
        "schema_version",
        "record_type",
        "qualification_effect",
        "sealed_utc",
        "publication_root_identity",
        "receipt_name",
        "receipt_sha256",
        "receipt_identity",
        "controller_state_attestation",
    }
    _require(type(value) is dict and set(value) == expected_keys, "selection success seal schema drifted")
    _require(
        type(value["schema_version"]) is int
        and value["schema_version"] == 1
        and value["record_type"] == SUCCESS_SEAL_RECORD_TYPE
        and value["qualification_effect"] == QUALIFICATION_EFFECT
        and value["receipt_name"] == receipt_name
        and value["receipt_sha256"] == _sha256(receipt_payload)
        and value["receipt_identity"]
        == {"device": receipt_identity[0], "inode": receipt_identity[1]},
        "selection success seal binding drifted",
    )
    _require_canonical_utc(value["sealed_utc"], "selection success seal sealed_utc")
    pilot_publisher._require_directory_identity(
        value["publication_root_identity"],
        publication_descriptor,
        "selection success seal publication root",
    )
    normalized = copy.deepcopy(value)
    normalized["controller_state_attestation"] = _binding(
        value["controller_state_attestation"], "controller-state attestation"
    )
    return normalized


def _publish_selection_success_seal_at(
    acceptance_descriptor: int,
    publication_descriptor: int,
    receipt_name: str,
    receipt_payload: bytes,
    receipt_identity: tuple[int, int],
    controller_state_attestation: object,
    *,
    sealed_utc: str,
) -> tuple[int, int]:
    pilot_publisher._require_same_file_at(
        publication_descriptor,
        receipt_name,
        receipt_identity,
        "canonical pressure-selection receipt",
    )
    seal_name = _selection_success_seal_name(receipt_name)
    staging_name = f".{seal_name}.staging-{uuid.uuid4()}"
    pilot_publisher._require_absent_at(
        acceptance_descriptor, seal_name, "pressure-selection durable success seal"
    )
    payload = _selection_success_seal_payload(
        publication_descriptor,
        receipt_name,
        receipt_payload,
        receipt_identity,
        controller_state_attestation,
        sealed_utc=sealed_utc,
    )
    staging_identity: tuple[int, int] | None = None
    committed = False
    try:
        pilot_publisher._write_exclusive_at(acceptance_descriptor, staging_name, payload)
        staging_identity = pilot_publisher._file_identity_at(
            acceptance_descriptor, staging_name, "staged pressure-selection success seal"
        )
        pilot_publisher._fsync_descriptor(acceptance_descriptor)
        try:
            pilot_publisher._rename_no_replace_at(
                acceptance_descriptor, staging_name, seal_name
            )
        except BaseException:
            pilot_publisher._require_absent_at(
                acceptance_descriptor, staging_name, "staged pressure-selection success seal"
            )
            pilot_publisher._require_same_file_at(
                acceptance_descriptor,
                seal_name,
                staging_identity,
                "pressure-selection durable success seal",
            )
            observed, _ = _read_exact_mode_regular_at(
                acceptance_descriptor,
                seal_name,
                mode=0o444,
                label="pressure-selection durable success seal",
            )
            _require(observed == payload, "pressure-selection success seal drifted during commit")
            pilot_publisher._fsync_descriptor(acceptance_descriptor)
            committed = True
            return staging_identity
        pilot_publisher._require_absent_at(
            acceptance_descriptor, staging_name, "staged pressure-selection success seal"
        )
        observed, _ = _read_exact_mode_regular_at(
            acceptance_descriptor,
            seal_name,
            mode=0o444,
            label="pressure-selection durable success seal",
        )
        _, identity = _read_exact_mode_regular_at(
            acceptance_descriptor,
            seal_name,
            mode=0o444,
            label="pressure-selection durable success seal",
        )
        _require(
            identity == staging_identity and observed == payload,
            "pressure-selection success seal drifted after commit",
        )
        pilot_publisher._fsync_descriptor(acceptance_descriptor)
        committed = True
        return staging_identity
    finally:
        if not committed:
            _unlink_exact_file_at(
                acceptance_descriptor, staging_name, staging_identity
            )


def _read_selection_success_seal_at(
    acceptance_descriptor: int,
    publication_descriptor: int,
    receipt_name: str,
    receipt_payload: bytes,
    receipt_identity: tuple[int, int],
) -> tuple[dict[str, object], tuple[int, int]]:
    payload, identity = _read_exact_mode_regular_at(
        acceptance_descriptor,
        _selection_success_seal_name(receipt_name),
        mode=0o444,
        label="pressure-selection durable success seal",
    )
    value = _decode_canonical_object(payload, "pressure-selection durable success seal")
    return (
        _validate_selection_success_seal(
            value,
            publication_descriptor,
            receipt_name,
            receipt_payload,
            receipt_identity,
        ),
        identity,
    )


def build_pressure_selection_receipt(
    *,
    published_pressure_pilot_receipt: Mapping[str, object],
    published_pressure_pilot_review_packet_receipt: Mapping[str, object],
    pilot_bundle_manifest_sha256: str,
    aggregate_pilot_analysis_sha256: str,
    case_descriptors: Sequence[Mapping[str, object]],
    authoritative_reanalysis_attestation: Mapping[str, object],
    reviewer_attestation: Mapping[str, object],
    selected_case: Mapping[str, object] = REVIEWED_SELECTED_CASE,
) -> dict[str, object]:
    """Build the unchanged schema-v3 receipt for the exact reviewed choice."""
    _require(
        dict(selected_case) == REVIEWED_SELECTED_CASE,
        "pressure-selection choice differs from the reviewed p0=1.0 baseline",
    )
    return {
        "schema_version": 3,
        "record_type": pressure_selection.RECORD_TYPE,
        "selection_method": pressure_selection.SELECTION_METHOD,
        "published_pressure_pilot_receipt": copy.deepcopy(
            dict(published_pressure_pilot_receipt)
        ),
        "published_pressure_pilot_review_packet_receipt": copy.deepcopy(
            dict(published_pressure_pilot_review_packet_receipt)
        ),
        "pilot_bundle_manifest_sha256": pilot_bundle_manifest_sha256,
        "aggregate_pilot_analysis_sha256": aggregate_pilot_analysis_sha256,
        "case_descriptors": [copy.deepcopy(dict(item)) for item in case_descriptors],
        "selected_case": dict(REVIEWED_SELECTED_CASE),
        "authoritative_reanalysis_attestation": copy.deepcopy(
            dict(authoritative_reanalysis_attestation)
        ),
        "reviewer_attestation": copy.deepcopy(dict(reviewer_attestation)),
    }


def _validate_reviewed_selection(
    receipt: object,
    *,
    authorized_pic_root: Path,
) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    try:
        normalized = pressure_selection.validate_pressure_selection_receipt(
            receipt, authorized_pic_root=authorized_pic_root
        )
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "schema-v3 pressure-selection receipt failed immutable verification"
        ) from error
    _require(
        normalized["selected_case"] == REVIEWED_SELECTED_CASE,
        "pressure-selection choice differs from the reviewed p0=1.0 baseline",
    )
    verifier = pressure_selection.pressure_review_packet_verifier
    try:
        reanalysis = verifier.consume_sealed_pressure_reanalysis_attestation(
            normalized["authoritative_reanalysis_attestation"],
            aggregate_receipt_binding=normalized["published_pressure_pilot_receipt"],
            packet_receipt_binding=normalized[
                "published_pressure_pilot_review_packet_receipt"
            ],
            pilot_bundle_manifest_sha256=normalized["pilot_bundle_manifest_sha256"],
            aggregate_pilot_analysis_sha256=normalized[
                "aggregate_pilot_analysis_sha256"
            ],
            authorized_pic_root=authorized_pic_root,
        )
        reviewer = verifier.consume_sealed_pressure_reviewer_attestation(
            normalized["reviewer_attestation"],
            aggregate_receipt_binding=normalized["published_pressure_pilot_receipt"],
            packet_receipt_binding=normalized[
                "published_pressure_pilot_review_packet_receipt"
            ],
            reanalysis_verification=reanalysis,
            selected_case=normalized["selected_case"],
            authorized_pic_root=authorized_pic_root,
        )
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "pressure-selection attestations failed immutable verification"
        ) from error
    _require(
        type(reanalysis) is dict
        and reanalysis.get("binding")
        == normalized["authoritative_reanalysis_attestation"],
        "authoritative reanalysis attestation binding drifted",
    )
    _require(
        type(reviewer) is dict
        and reviewer.get("binding") == normalized["reviewer_attestation"]
        and reviewer.get("selected_case") == REVIEWED_SELECTED_CASE,
        "pressure-selection reviewer attestation binding drifted",
    )
    _require(
        reviewer.get("reviewer_id") == REVIEWED_REVIEWER_ID,
        "pressure-selection reviewer differs from the reviewed human",
    )
    _require(
        reviewer.get("rationale") == REVIEWED_RATIONALE,
        "pressure-selection rationale differs from the reviewed rationale",
    )
    return normalized, reanalysis, reviewer


def _read_stable_bounded_regular_file_below(
    path: Path,
    root: Path,
    label: str,
    *,
    max_bytes: int,
) -> bytes:
    lexical_root = Path(os.path.abspath(root))
    lexical_root.resolve(strict=True)
    lexical_path = Path(os.path.abspath(path))
    try:
        relative = lexical_path.relative_to(lexical_root)
    except ValueError as error:
        raise PressureSelectionPublicationError(
            f"{label} is outside the authorized root"
        ) from error
    _require(bool(relative.parts), f"{label} is not below the authorized root")
    root_descriptor = os.open(lexical_root, os.O_RDONLY | os.O_DIRECTORY)
    directory_descriptor = root_descriptor
    descriptor: int | None = None
    try:
        for part in relative.parts[:-1]:
            next_descriptor = os.open(part, _DIRECTORY_FLAGS, dir_fd=directory_descriptor)
            if directory_descriptor != root_descriptor:
                os.close(directory_descriptor)
            directory_descriptor = next_descriptor
        descriptor = os.open(
            relative.parts[-1], _FILE_FLAGS, dir_fd=directory_descriptor
        )
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode) and not before.st_mode & 0o222,
            f"{label} is not a read-only regular file",
        )
        payload = pilot_publisher._read_bounded_descriptor(
            descriptor, label, max_bytes
        )
        after = os.fstat(descriptor)
        stable_fields = (
            "st_dev",
            "st_ino",
            "st_mode",
            "st_size",
            "st_mtime_ns",
            "st_ctime_ns",
        )
        _require(
            all(
                getattr(before, field) == getattr(after, field)
                for field in stable_fields
            )
            and len(payload) == after.st_size,
            f"{label} changed while reading",
        )
        return payload
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if directory_descriptor != root_descriptor:
            os.close(directory_descriptor)
        os.close(root_descriptor)


def _live_file(
    path: Path,
    root: Path,
    label: str,
    *,
    max_bytes: int | None = None,
) -> bytes:
    try:
        canonical = control_plane_common.require_canonical_path_below(path, root)
        if max_bytes is not None:
            return _read_stable_bounded_regular_file_below(
                canonical, root, label, max_bytes=max_bytes
            )
        return control_plane_common.read_stable_regular_file_below(
            canonical, root, require_read_only_mode=True
        )
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(f"{label} is unavailable") from error


def _capture_mirrored_ledger_state(
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> dict[str, object]:
    ledger_root = Path(os.path.abspath(authorized_pic_root)) / "ledger"
    project_ledger_root = (
        control_plane_common.project_home_ledger_root(authorized_project_home_root)
        / "ledger"
    )
    ledger_jsonl = ledger_root / "node_hours.jsonl"
    receipts_jsonl = ledger_root / "mirror_receipts.jsonl"
    mirror_jsonl = project_ledger_root / "node_hours.jsonl"
    pending_submission = ledger_root / "pending_submission.json"
    try:
        os.stat(pending_submission, follow_symlinks=False)
    except FileNotFoundError:
        pass
    except OSError as error:
        raise PressureSelectionPublicationError(
            "pending scheduler-submission marker cannot be inspected"
        ) from error
    else:
        raise PressureSelectionPublicationError(
            "pending scheduler submission blocks pressure selection"
        )
    try:
        ledger.require_no_incomplete_manual_accounting_marker(
            ledger_jsonl, mirror_jsonl
        )
        with ledger.validated_read_only_mirrored_state_snapshot(
            ledger_jsonl,
            receipts_jsonl,
            mirror_jsonl,
            ledger_root=authorized_pic_root,
            receipts_root=authorized_pic_root,
            mirror_root=control_plane_common.project_home_ledger_root(
                authorized_project_home_root
            ),
        ) as records:
            _require(bool(records), "pressure selection requires initialized ledgers")
            latest = ledger.latest_reservations(records)
            active_ids = sorted(
                reservation_id
                for reservation_id, record in latest.items()
                if record.get("state") in {"reserved", "submitted"}
            )
            totals = ledger.accounting(records)
            _require(
                active_ids == []
                and totals.get("currently_reserved_node_hours") == 0.0,
                "outstanding reservation blocks pressure selection",
            )
            cumulative = totals.get("cumulative_consumed_node_hours")
            _require(
                type(cumulative) in {int, float}
                and type(cumulative) is not bool
                and cumulative >= 0.0,
                "mirrored ledger cumulative accounting is malformed",
            )
            tail = records[-1].get("event_sha256")
            _lowercase_sha256(tail, "mirrored ledger tail event SHA-256")
            return {
                "validator": "validated_read_only_mirrored_state_snapshot",
                "orion_ledger_path": str(ledger_jsonl),
                "orion_mirror_receipts_path": str(receipts_jsonl),
                "project_home_ledger_path": str(mirror_jsonl),
                "record_count": len(records),
                "tail_event_sha256": tail,
                "active_reservation_ids": [],
                "currently_reserved_node_hours": 0.0,
                "cumulative_consumed_node_hours": float(cumulative),
                "pending_submission_marker": "absent",
                "pending_manual_accounting_markers": {
                    "orion": "absent",
                    "project_home": "absent",
                },
            }
    except PressureSelectionPublicationError:
        raise
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "mirrored ledger state failed pressure-selection verification"
        ) from error


def _capture_live_controller_state(
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> dict[str, object]:
    pic_root = Path(os.path.abspath(authorized_pic_root))
    project_root = Path(os.path.abspath(authorized_project_home_root))
    try:
        control_plane_common.require_no_active_promotion_transaction(
            authorized_pic_root=pic_root,
            authorized_project_home_root=project_root,
        )
    except (OSError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "active-policy promotion transaction blocks pressure selection"
        ) from error
    policy_path = control_plane_common.canonical_policy_path(pic_root)
    mirror_policy_path = control_plane_common.canonical_policy_path(project_root)
    promotion_path = control_plane_common.active_promotion_path(pic_root)
    mirror_promotion_path = control_plane_common.active_promotion_path(project_root)
    policy_payload = _live_file(policy_path, pic_root, "active Orion policy")
    mirror_policy_payload = _live_file(
        mirror_policy_path, project_root, "active Project Home policy"
    )
    promotion_payload = _live_file(promotion_path, pic_root, "active Orion promotion")
    mirror_promotion_payload = _live_file(
        mirror_promotion_path, project_root, "active Project Home promotion"
    )
    _require(
        policy_payload == mirror_policy_payload,
        "Orion and Project Home active-policy bytes differ",
    )
    _require(
        promotion_payload == mirror_promotion_payload,
        "Orion and Project Home active-promotion bytes differ",
    )
    promotion = _decode_canonical_object(promotion_payload, "active promotion")
    version = _lowercase_sha256(
        promotion.get("control_plane_version"), "active control-plane version"
    )
    try:
        policy, anchors = control_plane_common.require_storage_policy_unlock_snapshot(
            control_plane_version=version,
            authorized_pic_root=pic_root,
            authorized_project_home_root=project_root,
            allow_pending_genesis=True,
        )
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "active storage policy is not unlocked for pressure selection"
        ) from error
    _require(
        type(policy) is dict
        and policy == _decode_canonical_object(policy_payload, "active policy"),
        "validated active policy differs from captured bytes",
    )
    _require(
        type(anchors) is dict
        and anchors
        == {
            "active_policy_sha256": _sha256(policy_payload),
            "active_promotion_sha256": _sha256(promotion_payload),
        },
        "active policy or promotion anchor drifted during capture",
    )
    slices = policy.get("registered_science_slices")
    _require(
        slices == [],
        "pressure-selection publication requires an empty registered-science allowlist",
    )
    admission_smoke = policy.get("frontier_admission_smoke")
    _require(
        admission_smoke
        == {"status": control_plane_common.CLOSED_ADMISSION_SMOKE_STATUS},
        "pressure-selection publication requires closed admission-smoke authority",
    )
    freeze = policy.get("science_submission_freeze")
    _require(
        type(freeze) is dict
        and set(freeze)
        == {
            "status",
            "manifest_path",
            "manifest_sha256",
            "build_profile_control_plane_version",
        }
        and freeze.get("status") == control_plane_common.AUTHORIZED_CLEAN_CANDIDATE_FREEZE
        and freeze.get("build_profile_control_plane_version") == version,
        "pressure-selection publication requires one authorized clean-candidate freeze",
    )
    _lowercase_sha256(freeze.get("manifest_sha256"), "clean-candidate manifest SHA-256")
    storage = policy.get("olcf_side_storage")
    manual_accounting_authorizations = (
        storage.get("manual_accounting_authorizations")
        if type(storage) is dict
        else None
    )
    _require(
        type(storage) is dict
        and storage.get("installed_control_plane_version") == version
        and type(manual_accounting_authorizations) is list,
        "pressure-selection publication requires a validated accounting registry",
    )
    mirrored_ledger_state = _capture_mirrored_ledger_state(
        authorized_pic_root=pic_root,
        authorized_project_home_root=project_root,
    )
    installed_records = []
    inventory_payloads = []
    for role, root in (("orion", pic_root), ("project_home", project_root)):
        directory = root / "control_plane" / version
        try:
            inventory = control_plane_common.verify_installed_control_plane(
                directory, authorized_pic_root=root
            )
        except (OSError, TypeError, ValueError) as error:
            raise PressureSelectionPublicationError(
                f"{role} installed control plane failed verification"
            ) from error
        inventory_path = directory / "inventory.json"
        inventory_payload = _live_file(
            inventory_path, root, f"{role} installed control-plane inventory"
        )
        _require(
            inventory == _decode_canonical_object(
                inventory_payload, f"{role} installed control-plane inventory"
            )
            and inventory.get("version") == version,
            f"{role} installed control-plane inventory drifted",
        )
        inventory_payloads.append(inventory_payload)
        installed_records.append(
            {
                "root_role": role,
                "inventory_path": str(inventory_path),
                "inventory_sha256": _sha256(inventory_payload),
            }
        )
    _require(
        inventory_payloads[0] == inventory_payloads[1],
        "installed control-plane inventory bytes differ",
    )
    manifest_path = Path(str(freeze["manifest_path"]))
    try:
        revalidation = clean_candidate_revalidator.revalidate_clean_candidate(
            manifest_path,
            expected_manifest_sha256=str(freeze["manifest_sha256"]),
            expected_receipt_control_plane_version=version,
            control_plane_dir=pic_root / "control_plane" / version,
            authorized_pic_root=pic_root,
            authorized_project_home_root=project_root,
            read_candidate_tree=_read_bounded_clean_candidate_tree,
        )
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "authorized clean candidate failed read-only revalidation"
        ) from error
    manifest_payload = _live_file(
        manifest_path, pic_root, "authorized clean-candidate manifest"
    )
    manifest = _decode_canonical_object(
        manifest_payload, "authorized clean-candidate manifest"
    )
    source = manifest.get("source")
    _require(
        type(revalidation) is dict,
        "clean-candidate revalidation report drifted",
    )
    revalidated_binding = revalidation.get("clean_candidate_manifest")
    revalidated_source = revalidation.get("source")
    _require(
        type(revalidation) is dict
        and revalidation.get("status") == "passed"
        and revalidation.get("current_control_plane_version") == version
        and type(revalidated_binding) is dict
        and revalidated_binding.get("path") == str(manifest_path)
        and revalidated_binding.get("sha256") == freeze["manifest_sha256"]
        and revalidated_binding.get("expected_sha256") == freeze["manifest_sha256"]
        and type(revalidated_source) is dict
        and type(source) is dict
        and revalidated_source.get("git_commit") == source.get("git_commit"),
        "clean-candidate revalidation report drifted",
    )
    _require(
        _sha256(manifest_payload) == freeze["manifest_sha256"],
        "clean-candidate manifest bytes differ from the active freeze",
    )
    clean_candidate = {
        "manifest_path": str(manifest_path),
        "manifest_sha256": str(freeze["manifest_sha256"]),
        "git_commit": str(source.get("git_commit", "")),
        "source_archive_sha256": str(source.get("archive_sha256", "")),
    }
    _require(
        _GIT_COMMIT_PATTERN.fullmatch(clean_candidate["git_commit"]) is not None,
        "clean-candidate Git commit is malformed",
    )
    _lowercase_sha256(
        clean_candidate["source_archive_sha256"],
        "clean-candidate source archive SHA-256",
    )
    try:
        control_plane_common.require_no_active_promotion_transaction(
            authorized_pic_root=pic_root,
            authorized_project_home_root=project_root,
        )
    except (OSError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "active-policy promotion transaction appeared during pressure selection"
        ) from error
    return {
        "controller_state": {
            "control_plane_version": version,
            "installed_control_planes": installed_records,
            "active_policy": {
                "orion_path": str(policy_path),
                "project_home_path": str(mirror_policy_path),
                "sha256": _sha256(policy_payload),
            },
            "active_promotion": {
                "orion_path": str(promotion_path),
                "project_home_path": str(mirror_promotion_path),
                "sha256": _sha256(promotion_payload),
            },
            "registered_science_slices": [],
            "frontier_admission_smoke": copy.deepcopy(admission_smoke),
            "science_submission_freeze": copy.deepcopy(freeze),
            "clean_candidate": clean_candidate,
            "mirrored_ledger_state": mirrored_ledger_state,
            "pending_submission_marker": "absent",
            "manual_accounting_authorizations": copy.deepcopy(
                manual_accounting_authorizations
            ),
            "pending_manual_accounting_marker": "absent",
            "active_promotion_transaction": "absent",
        },
        "active_policy_payload": policy_payload,
        "active_promotion_payload": promotion_payload,
    }


def _source_binding(
    reanalysis: Mapping[str, object],
    controller_state: Mapping[str, object],
    *,
    candidate_publication_authorization: Mapping[str, object],
    authorized_pic_root: Path,
) -> dict[str, object]:
    authorization = reanalysis.get("source_authorization")
    _require(
        type(authorization) is dict,
        "authoritative reanalysis source authorization is absent",
    )
    git_commit = authorization.get("git_commit")
    _require(
        type(git_commit) is str and _GIT_COMMIT_PATTERN.fullmatch(git_commit) is not None,
        "authoritative reanalysis Git commit is malformed",
    )
    candidate = controller_state["clean_candidate"]
    _require(type(candidate) is dict, "captured clean-candidate binding is malformed")
    _require(
        git_commit == candidate.get("git_commit")
        and authorization.get("source_archive_sha256")
        == candidate.get("source_archive_sha256"),
        "authoritative reanalysis source differs from the authorized clean candidate",
    )
    expected_authorization, _archive_files = _candidate_reanalysis_source_authorization(
        controller_state,
        authorized_pic_root=authorized_pic_root,
    )
    _require(
        authorization == expected_authorization,
        "authoritative reanalysis source closure differs from the active clean-candidate archive",
    )
    try:
        pressure_selection.pressure_review_packet_verifier.validate_pressure_reanalysis_source_snapshot(
            reanalysis,
            git_commit=candidate.get("git_commit"),
            source_archive_sha256=candidate.get("source_archive_sha256"),
            helper_source_closure=expected_authorization["source_closure"],
            authorized_pic_root=authorized_pic_root,
        )
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "authoritative reanalysis source closure differs from the authorized clean candidate"
        ) from error
    return {
        "authoritative_reanalysis_attestation": _binding(
            reanalysis.get("binding"), "authoritative reanalysis attestation"
        ),
        "git_commit": git_commit,
        "source_archive_sha256": _lowercase_sha256(
            authorization.get("source_archive_sha256"),
            "authoritative reanalysis source archive SHA-256",
        ),
        "reanalysis_source_closure_sha256": _lowercase_sha256(
            authorization.get("source_closure_sha256"),
            "authoritative reanalysis source closure SHA-256",
        ),
        "clean_candidate_manifest": {
            "path": candidate["manifest_path"],
            "sha256": candidate["manifest_sha256"],
        },
        "candidate_publication_authorization": copy.deepcopy(
            candidate_publication_authorization
        ),
    }


def build_controller_state_attestation(
    capture: Mapping[str, object],
    *,
    reanalysis_verification: Mapping[str, object],
    publisher_source_authentication: Mapping[str, object],
    candidate_publication_authorization: Mapping[str, object],
    operator_id: str,
    captured_utc: str,
    sealed_utc: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    """Build one launch-prohibited immutable controller-state attestation."""
    controller_state = copy.deepcopy(capture["controller_state"])
    _require(type(controller_state) is dict, "captured controller state is malformed")
    return {
        "schema_version": 1,
        "record_type": CONTROLLER_STATE_RECORD_TYPE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "operator_id": _operator_id(operator_id, "controller-state operator ID"),
        "captured_utc": _require_canonical_utc(
            captured_utc, "controller-state captured_utc"
        ),
        "sealed_utc": _require_canonical_utc(sealed_utc, "controller-state sealed_utc"),
        "publisher_source_authentication": _validate_publisher_source_authentication(
            dict(publisher_source_authentication)
        ),
        "controller_state": controller_state,
        "source_binding": _source_binding(
            reanalysis_verification,
            controller_state,
            candidate_publication_authorization=candidate_publication_authorization,
            authorized_pic_root=authorized_pic_root,
        ),
    }


def _controller_attestation_directory_name(sealed_utc: str, operator_id: str) -> str:
    timestamp = datetime.strptime(sealed_utc, "%Y-%m-%dT%H:%M:%SZ").strftime(
        "%Y%m%dT%H%M%SZ"
    )
    return (
        f"{timestamp}-q011-section54-pressure-selection-controller-state-{operator_id}"
    )


def _remove_owned_partial_attestation_at(
    archive_descriptor: int,
    directory_name: str,
    directory_identity: tuple[int, int],
) -> None:
    descriptor = os.open(directory_name, _DIRECTORY_FLAGS, dir_fd=archive_descriptor)
    try:
        metadata = os.fstat(descriptor)
        _require(
            (metadata.st_dev, metadata.st_ino) == directory_identity,
            "partial controller-state attestation directory changed during cleanup",
        )
        os.fchmod(descriptor, 0o700)
        names = set(os.listdir(descriptor))
        _require(
            names
            <= {"attestation.json", "active_policy.json", "active_promotion.json"},
            "partial controller-state attestation gained an unexpected member",
        )
        for name in names:
            observed = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
            _require(
                stat.S_ISREG(observed.st_mode) and observed.st_nlink == 1,
                "partial controller-state attestation member is unsafe to remove",
            )
            os.unlink(name, dir_fd=descriptor)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    current = os.stat(directory_name, dir_fd=archive_descriptor, follow_symlinks=False)
    _require(
        stat.S_ISDIR(current.st_mode)
        and (current.st_dev, current.st_ino) == directory_identity,
        "partial controller-state attestation directory changed before cleanup",
    )
    os.rmdir(directory_name, dir_fd=archive_descriptor)
    pilot_publisher._fsync_descriptor(archive_descriptor)


def _write_controller_state_attestation_at(
    archive_root: Path,
    archive_descriptor: int,
    value: Mapping[str, object],
    active_policy_payload: bytes,
    active_promotion_payload: bytes,
) -> dict[str, str]:
    operator_id = str(value["operator_id"])
    sealed_utc = str(value["sealed_utc"])
    directory_name = _controller_attestation_directory_name(sealed_utc, operator_id)
    payload = _canonical_json_bytes(value)
    pilot_publisher._require_absent_at(
        archive_descriptor, directory_name, "controller-state attestation"
    )
    os.mkdir(directory_name, mode=0o700, dir_fd=archive_descriptor)
    descriptor: int | None = None
    identity: tuple[int, int] | None = None
    sealed = False
    try:
        pilot_publisher._fsync_descriptor(archive_descriptor)
        descriptor = os.open(directory_name, _DIRECTORY_FLAGS, dir_fd=archive_descriptor)
        metadata = os.fstat(descriptor)
        identity = (metadata.st_dev, metadata.st_ino)
        _write_controller_member_at(descriptor, "active_policy.json", active_policy_payload)
        _write_controller_member_at(
            descriptor, "active_promotion.json", active_promotion_payload
        )
        _write_controller_member_at(descriptor, "attestation.json", payload)
        _require(
            set(os.listdir(descriptor))
            == {"attestation.json", "active_policy.json", "active_promotion.json"},
            "controller-state attestation member closure drifted",
        )
        os.fchmod(descriptor, 0o500)
        os.fsync(descriptor)
        pilot_publisher._fsync_descriptor(archive_descriptor)
        pilot_publisher._require_same_directory_at(
            archive_descriptor,
            directory_name,
            descriptor,
            "controller-state attestation",
        )
        sealed = True
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if not sealed:
            if identity is None:
                os.rmdir(directory_name, dir_fd=archive_descriptor)
                pilot_publisher._fsync_descriptor(archive_descriptor)
            else:
                _remove_owned_partial_attestation_at(
                    archive_descriptor, directory_name, identity
                )
    return {
        "path": str(archive_root / directory_name / "attestation.json"),
        "sha256": _sha256(payload),
    }


def _open_or_create_private_root(
    pic_root: Path,
    name: str,
    label: str,
) -> tuple[Path, int]:
    _require(
        "/" not in name and name not in {"", ".", ".."},
        f"{label} name is invalid",
    )
    root = pic_root / name
    parent_descriptor = pilot_publisher._open_absolute_directory(pic_root)
    descriptor: int | None = None
    try:
        pilot_publisher._require_same_directory(
            pic_root, parent_descriptor, "authorized PIC root"
        )
        try:
            os.mkdir(name, mode=0o700, dir_fd=parent_descriptor)
        except FileExistsError:
            pass
        else:
            pilot_publisher._fsync_descriptor(parent_descriptor)
        descriptor = os.open(name, _DIRECTORY_FLAGS, dir_fd=parent_descriptor)
        metadata = os.fstat(descriptor)
        _require(
            stat.S_ISDIR(metadata.st_mode)
            and stat.S_IMODE(metadata.st_mode) == 0o700
            and metadata.st_uid == os.geteuid(),
            f"{label} must be one private directory with mode 0700",
        )
        pilot_publisher._require_same_directory_at(
            parent_descriptor, name, descriptor, label
        )
        return root, descriptor
    except BaseException:
        if descriptor is not None:
            os.close(descriptor)
        raise
    finally:
        os.close(parent_descriptor)


def _write_single_file_attestation_at(
    archive_root: Path,
    archive_descriptor: int,
    directory_name: str,
    value: Mapping[str, object],
) -> dict[str, str]:
    payload = _canonical_json_bytes(value)
    pilot_publisher._require_absent_at(
        archive_descriptor, directory_name, "pressure-gate attestation"
    )
    os.mkdir(directory_name, mode=0o700, dir_fd=archive_descriptor)
    descriptor: int | None = None
    identity: tuple[int, int] | None = None
    sealed = False
    try:
        pilot_publisher._fsync_descriptor(archive_descriptor)
        descriptor = os.open(directory_name, _DIRECTORY_FLAGS, dir_fd=archive_descriptor)
        metadata = os.fstat(descriptor)
        identity = (metadata.st_dev, metadata.st_ino)
        _write_controller_member_at(descriptor, "attestation.json", payload)
        _require(
            os.listdir(descriptor) == ["attestation.json"],
            "pressure-gate attestation member closure drifted",
        )
        os.fchmod(descriptor, 0o500)
        os.fsync(descriptor)
        pilot_publisher._fsync_descriptor(archive_descriptor)
        pilot_publisher._require_same_directory_at(
            archive_descriptor,
            directory_name,
            descriptor,
            "pressure-gate attestation",
        )
        sealed = True
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if not sealed:
            if identity is None:
                os.rmdir(directory_name, dir_fd=archive_descriptor)
                pilot_publisher._fsync_descriptor(archive_descriptor)
            else:
                _remove_owned_partial_attestation_at(
                    archive_descriptor, directory_name, identity
                )
    return {
        "path": str(archive_root / directory_name / "attestation.json"),
        "sha256": _sha256(payload),
    }


def _validate_stage4_preparation_attestation(value: object) -> dict[str, object]:
    _require(
        type(value) is dict
        and set(value)
        == {
            "schema_version",
            "record_type",
            "qualification_effect",
            "operator_id",
            "sealed_utc",
            "publisher_source_authentication",
            "authoritative_reanalysis_attestation",
        }
        and type(value["schema_version"]) is int
        and value["schema_version"] == 1
        and value["record_type"] == STAGE4_PREPARATION_RECORD_TYPE
        and value["qualification_effect"] == QUALIFICATION_EFFECT,
        "Stage-4 preparation attestation schema drifted",
    )
    normalized = copy.deepcopy(value)
    normalized["operator_id"] = _operator_id(
        value["operator_id"], "Stage-4 preparation operator ID"
    )
    normalized["sealed_utc"] = _require_canonical_utc(
        value["sealed_utc"], "Stage-4 preparation sealed_utc"
    )
    normalized["publisher_source_authentication"] = (
        _validate_publisher_source_authentication(
            value["publisher_source_authentication"]
        )
    )
    normalized["authoritative_reanalysis_attestation"] = _binding(
        value["authoritative_reanalysis_attestation"],
        "Stage-4 preparation authoritative reanalysis attestation",
    )
    return normalized


def _consume_stage4_preparation_attestation(
    attestation_binding: object,
    *,
    authorized_pic_root: Path,
) -> dict[str, object]:
    binding = _binding(attestation_binding, "Stage-4 preparation attestation")
    pic_root = Path(os.path.abspath(authorized_pic_root))
    archive_root = pic_root / PRESSURE_GATE_ATTESTATION_ROOT_NAME
    path = Path(binding["path"])
    _require(
        path.name == "attestation.json"
        and path.parent.parent == archive_root
        and path == Path(os.path.abspath(path)),
        "Stage-4 preparation attestation path is outside the archive",
    )
    try:
        canonical = control_plane_common.require_canonical_path_below(path, archive_root)
    except (OSError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "Stage-4 preparation attestation path is not canonical"
        ) from error
    descriptor = os.open(canonical.parent, _DIRECTORY_FLAGS)
    try:
        metadata = os.fstat(descriptor)
        _require(
            stat.S_IMODE(metadata.st_mode) == 0o500
            and metadata.st_uid == os.geteuid()
            and os.listdir(descriptor) == ["attestation.json"],
            "Stage-4 preparation attestation directory closure drifted",
        )
        payload, _ = _read_exact_mode_regular_at(
            descriptor,
            "attestation.json",
            mode=0o400,
            label="Stage-4 preparation attestation",
        )
        _require(
            _sha256(payload) == binding["sha256"],
            "Stage-4 preparation attestation binding hash mismatch",
        )
        attestation = _validate_stage4_preparation_attestation(
            _decode_canonical_object(payload, "Stage-4 preparation attestation")
        )
        _require(
            canonical.parent.name
            == _pressure_gate_directory_name(
                str(attestation["sealed_utc"]),
                "stage4-preparation",
                str(attestation["operator_id"]),
            )
            and os.listdir(descriptor) == ["attestation.json"],
            "Stage-4 preparation attestation changed during verification",
        )
    finally:
        os.close(descriptor)
    return {"binding": binding, "attestation": attestation}


def _write_candidate_receipt_at(
    candidate_root: Path,
    candidate_descriptor: int,
    receipt: Mapping[str, object],
) -> dict[str, str]:
    payload = pressure_selection.canonical_json_bytes(dict(receipt))
    name = f"{Path(CANONICAL_RECEIPT_NAME).stem}.candidate-{uuid.uuid4()}.json"
    _write_controller_member_at(candidate_descriptor, name, payload)
    pilot_publisher._fsync_descriptor(candidate_descriptor)
    observed, _ = _read_exact_mode_regular_at(
        candidate_descriptor,
        name,
        mode=0o400,
        label="candidate pressure-selection receipt",
    )
    _require(observed == payload, "candidate pressure-selection receipt drifted")
    return {"path": str(candidate_root / name), "sha256": _sha256(payload)}


def _validate_candidate_publication_authorization(value: object) -> dict[str, object]:
    _require(
        type(value) is dict
        and set(value)
        == {
            "schema_version",
            "record_type",
            "qualification_effect",
            "sealed_utc",
            "publisher_source_authentication",
            "stage4_preparation_attestation",
            "human_decision",
            "candidate_pressure_selection_receipt",
            "authoritative_reanalysis_attestation",
            "reviewer_attestation",
        }
        and type(value["schema_version"]) is int
        and value["schema_version"] == 1
        and value["record_type"] == CANDIDATE_AUTHORIZATION_RECORD_TYPE
        and value["qualification_effect"] == QUALIFICATION_EFFECT,
        "candidate publication authorization schema drifted",
    )
    normalized = copy.deepcopy(value)
    normalized["sealed_utc"] = _require_canonical_utc(
        value["sealed_utc"], "candidate publication authorization sealed_utc"
    )
    normalized["publisher_source_authentication"] = (
        _validate_publisher_source_authentication(
            value["publisher_source_authentication"]
        )
    )
    for key, label in (
        ("stage4_preparation_attestation", "Stage-4 preparation attestation"),
        ("human_decision", "human pressure-selection decision"),
        ("candidate_pressure_selection_receipt", "candidate pressure-selection receipt"),
        ("authoritative_reanalysis_attestation", "authoritative reanalysis attestation"),
        ("reviewer_attestation", "reviewer attestation"),
    ):
        normalized[key] = _binding(value[key], label)
    return normalized


def _authorize_candidate_publication(
    value: object,
    *,
    receipt_payload: bytes,
    receipt: Mapping[str, object],
    reanalysis_verification: Mapping[str, object],
    reviewer_verification: Mapping[str, object],
    publisher_source_authentication: Mapping[str, object],
    authorized_pic_root: Path,
) -> dict[str, object]:
    authorization = _validate_candidate_publication_authorization(value)
    _require(
        authorization["publisher_source_authentication"]
        == publisher_source_authentication
        and authorization["candidate_pressure_selection_receipt"]["sha256"]
        == _sha256(receipt_payload)
        and authorization["authoritative_reanalysis_attestation"]
        == receipt["authoritative_reanalysis_attestation"]
        and authorization["reviewer_attestation"] == receipt["reviewer_attestation"],
        "candidate publication authorization differs from the reviewed receipt or source",
    )
    preparation = _consume_stage4_preparation_attestation(
        authorization["stage4_preparation_attestation"],
        authorized_pic_root=authorized_pic_root,
    )
    human_decision, human_decision_binding = _load_human_decision(
        Path(str(authorization["human_decision"]["path"])),
        authorized_pic_root=authorized_pic_root,
    )
    _require(
        human_decision_binding == authorization["human_decision"],
        "candidate publication authorization human-decision binding drifted",
    )
    normalized_decision = _validate_human_decision(
        human_decision, reanalysis_verification=reanalysis_verification
    )
    _require(
        type(reviewer_verification) is dict
        and reviewer_verification.get("reviewer_id")
        == normalized_decision["reviewer_id"]
        and reviewer_verification.get("selected_case")
        == normalized_decision["selected_case"]
        and reviewer_verification.get("rationale") == normalized_decision["rationale"]
        and reviewer_verification.get("reviewed_utc")
        == normalized_decision["reviewed_utc"]
        and reviewer_verification.get("authoritative_reanalysis_attestation")
        == normalized_decision["authoritative_reanalysis_attestation"],
        "reviewer attestation differs from reopened human pressure-selection decision",
    )
    preparation_attestation = preparation["attestation"]
    _require(
        type(preparation_attestation) is dict
        and preparation_attestation["publisher_source_authentication"]
        == publisher_source_authentication
        and preparation_attestation["authoritative_reanalysis_attestation"]
        == receipt["authoritative_reanalysis_attestation"],
        "candidate publication source differs from Stage-4 preparation source",
    )
    _require(
        normalized_decision["stage4_preparation_attestation"]
        == authorization["stage4_preparation_attestation"],
        "candidate publication human decision binds another Stage-4 preparation",
    )
    return authorization


def _write_candidate_authorization_at(
    candidate_root: Path,
    candidate_descriptor: int,
    authorization: Mapping[str, object],
) -> dict[str, str]:
    payload = _canonical_json_bytes(authorization)
    name = f"pressure-selection-candidate-authorization-{uuid.uuid4()}.json"
    _write_controller_member_at(candidate_descriptor, name, payload)
    pilot_publisher._fsync_descriptor(candidate_descriptor)
    observed, _ = _read_exact_mode_regular_at(
        candidate_descriptor,
        name,
        mode=0o400,
        label="candidate publication authorization",
    )
    _require(observed == payload, "candidate publication authorization drifted")
    return {"path": str(candidate_root / name), "sha256": _sha256(payload)}


def _candidate_reanalysis_source_authorization(
    state: Mapping[str, object],
    *,
    authorized_pic_root: Path,
) -> tuple[dict[str, object], dict[str, bytes]]:
    candidate = state.get("clean_candidate")
    _require(type(candidate) is dict, "captured clean candidate is malformed")
    manifest_path = Path(str(candidate.get("manifest_path", "")))
    source_archive_path = manifest_path.parent / "source.tar"
    source_archive = _live_file(
        source_archive_path,
        authorized_pic_root,
        "authorized clean-candidate source archive",
        max_bytes=MAX_CANDIDATE_ARCHIVE_RAW_BYTES,
    )
    _require(
        _sha256(source_archive) == candidate.get("source_archive_sha256"),
        "authorized clean-candidate source archive hash drifted",
    )
    _require_unambiguous_plain_git_archive(source_archive)
    archive_files = _validated_candidate_archive_regular_files(source_archive)
    source_closure = []
    for relative in pressure_selection.pressure_review_packet_verifier.PRESSURE_REANALYSIS_SOURCE_PATHS:
        archived = archive_files.get(relative)
        _require(
            archived is not None,
            f"authorized clean-candidate source archive omits {relative}",
        )
        source_closure.append({"path": relative, "sha256": _sha256(archived)})
    closure_payload = json.dumps(
        source_closure,
        separators=(",", ":"),
        sort_keys=True,
        allow_nan=False,
    ).encode("utf-8")
    return (
        {
            "execution_mode": (
                pressure_selection.pressure_review_packet_verifier.PRESSURE_REANALYSIS_EXECUTION_MODE
            ),
            "git_commit": candidate["git_commit"],
            "source_archive_sha256": candidate["source_archive_sha256"],
            "source_closure_sha256": _sha256(closure_payload),
            "source_closure": source_closure,
            "historical_production_source_authorization": dict(
                pressure_selection.pressure_review_packet_verifier.AUTHORIZED_HISTORICAL_REANALYSIS_SOURCE_AUTHORIZATION
            ),
        },
        archive_files,
    )


def _validate_candidate_archive_resource_limits(
    source_archive: bytes,
) -> tuple[int, int]:
    try:
        _require(
            len(source_archive) <= MAX_CANDIDATE_ARCHIVE_RAW_BYTES,
            "authorized clean-candidate source archive exceeds resource limits",
        )
        member_count = 0
        regular_count = 0
        regular_bytes = 0
        with tarfile.open(fileobj=io.BytesIO(source_archive), mode="r:") as stream:
            for member in stream:
                member_count += 1
                _require(
                    member_count <= MAX_CANDIDATE_ARCHIVE_MEMBERS
                    and len(member.name)
                    <= MAX_CANDIDATE_ARCHIVE_MEMBER_NAME_CHARACTERS
                    and len(Path(member.name).parts)
                    <= MAX_CANDIDATE_ARCHIVE_PATH_DEPTH,
                    "authorized clean-candidate source archive exceeds resource limits",
                )
                if not member.isfile():
                    continue
                regular_count += 1
                regular_bytes += member.size
                _require(
                    regular_count <= MAX_CANDIDATE_ARCHIVE_REGULAR_FILES
                    and regular_bytes <= MAX_CANDIDATE_ARCHIVE_REGULAR_BYTES
                    and member.size <= MAX_CANDIDATE_ARCHIVE_MEMBER_BYTES
                    and member.size >= 0,
                    "authorized clean-candidate source archive exceeds resource limits",
                )
        return member_count, regular_bytes
    except PressureSelectionPublicationError:
        raise
    except (OSError, RecursionError, tarfile.TarError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "authorized clean-candidate source archive is unreadable"
        ) from error


def _require_unambiguous_plain_git_archive(source_archive: bytes) -> None:
    _require(
        len(source_archive) >= tarfile.BLOCKSIZE
        and source_archive.startswith(b"pax_global_header\0")
        and source_archive[156:157] == tarfile.XGLTYPE
        and source_archive[257:263] == b"ustar\0",
        "authorized clean-candidate archive must use canonical uncompressed Git tar format",
    )


def _validated_candidate_archive_regular_files(
    source_archive: bytes,
) -> dict[str, bytes]:
    _validate_candidate_archive_resource_limits(source_archive)
    try:
        records: dict[str, bytes] = {}
        names: set[str] = set()
        with tarfile.open(fileobj=io.BytesIO(source_archive), mode="r:") as stream:
            for member in stream.getmembers():
                path = control_plane_common.canonical_relative_posix_path(
                    member.name, field="candidate source-archive member"
                )
                canonical_name = path.as_posix()
                _require(
                    member.name == canonical_name and canonical_name not in names,
                    f"authorized clean-candidate source archive path drifted: {member.name!r}",
                )
                names.add(canonical_name)
                if not member.isfile():
                    continue
                extracted = stream.extractfile(member)
                _require(
                    extracted is not None,
                    f"authorized clean-candidate source archive member is unreadable: {member.name!r}",
                )
                payload = extracted.read(member.size + 1)
                _require(
                    len(payload) == member.size,
                    f"authorized clean-candidate source archive member size drifted: {member.name!r}",
                )
                records[canonical_name] = payload
        return records
    except PressureSelectionPublicationError:
        raise
    except (OSError, RecursionError, tarfile.TarError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "authorized clean-candidate source archive is unreadable"
        ) from error


# Keep bounded capture Stage-4-local because the shared control-plane bytes are
# already bound by historical registered-science records.
class _BoundedCleanCandidateClosure(
    control_plane_common._RetainedCleanCandidateClosure
):
    def __init__(self) -> None:
        super().__init__()
        self._stage4_archive_count = 0
        self._stage4_archive_members = 0
        self._stage4_archive_regular_bytes = 0
        self._stage4_total_bytes = 0

    def read_regular_file_at(
        self,
        directory_descriptor: int,
        name: str,
        *,
        label: str,
    ) -> bytes:
        self._require_no_watch_events()
        control_plane_common._clean_candidate_fixed_layout_name(name, label=label)
        descriptor = os.open(
            name,
            _FILE_FLAGS,
            dir_fd=directory_descriptor,
        )
        try:
            before = os.fstat(descriptor)
            control_plane_common._require_clean_candidate_regular_metadata(
                before, label=label
            )
            self._watch(descriptor)
            is_archive = name.endswith(".tar")
            _require(
                self._stage4_total_bytes + before.st_size
                <= MAX_CANDIDATE_TREE_REGULAR_BYTES,
                "authorized clean-candidate tree exceeds resource limits",
            )
            if is_archive:
                _require(
                    self._stage4_archive_count < MAX_CANDIDATE_TREE_ARCHIVES,
                    "authorized clean-candidate tree exceeds resource limits",
                )
            payload = pilot_publisher._read_bounded_descriptor(
                descriptor,
                label,
                (
                    MAX_CANDIDATE_ARCHIVE_RAW_BYTES
                    if is_archive
                    else MAX_CANDIDATE_TREE_MEMBER_BYTES
                ),
            )
            after = os.fstat(descriptor)
            metadata = control_plane_common._clean_candidate_stable_metadata(after)
            _require(
                metadata == control_plane_common._clean_candidate_stable_metadata(before)
                and len(payload) == after.st_size,
                f"{label} changed while reading",
            )
            self._stage4_total_bytes += len(payload)
            _require(
                self._stage4_total_bytes <= MAX_CANDIDATE_TREE_REGULAR_BYTES,
                "authorized clean-candidate tree exceeds resource limits",
            )
            if is_archive:
                self._stage4_archive_count += 1
                _require(
                    self._stage4_archive_count <= MAX_CANDIDATE_TREE_ARCHIVES,
                    "authorized clean-candidate tree exceeds resource limits",
                )
                _require_unambiguous_plain_git_archive(payload)
                member_count, regular_bytes = (
                    _validate_candidate_archive_resource_limits(payload)
                )
                self._stage4_archive_members += member_count
                self._stage4_archive_regular_bytes += regular_bytes
                _require(
                    self._stage4_archive_members
                    <= MAX_CANDIDATE_TREE_ARCHIVE_MEMBERS
                    and self._stage4_archive_regular_bytes
                    <= MAX_CANDIDATE_TREE_ARCHIVE_REGULAR_BYTES,
                    "authorized clean-candidate tree exceeds resource limits",
                )
            self._require_no_watch_events()
            self._files.append(
                (directory_descriptor, name, descriptor, metadata, label)
            )
            descriptor = -1
            return payload
        finally:
            if descriptor >= 0:
                os.close(descriptor)


def _read_bounded_clean_candidate_tree(
    candidate_manifest_path: Path,
    *,
    authorized_pic_root: Path,
) -> dict[str, object]:
    try:
        lexical_pic_root = Path(os.path.abspath(authorized_pic_root))
        candidate_root = lexical_pic_root / "clean_candidates"
        candidate_path = control_plane_common.require_canonical_path_below(
            Path(os.path.abspath(candidate_manifest_path)), candidate_root
        )
        _require(
            candidate_path.name == "clean_candidate_manifest.json"
            and candidate_path.parent.parent == candidate_root,
            "clean-candidate manifest must use the fixed candidate layout",
        )
        closure = _BoundedCleanCandidateClosure()
        try:
            candidate_root_ancestry = control_plane_common.PinnedDirectoryAncestry(
                candidate_root, root=lexical_pic_root
            )
        except BaseException:
            closure.close()
            raise
        with candidate_root_ancestry:
            try:
                closure.watch_ancestry(candidate_root_ancestry.descriptors)
                candidate_root_descriptor = candidate_root_ancestry.descriptor
                candidate_descriptor = closure.open_directory_at(
                    candidate_root_descriptor,
                    candidate_path.parent.name,
                    label="Clean-candidate directory",
                )
            except BaseException:
                closure.close(close_ancestry=candidate_root_ancestry.close)
                raise
            try:
                candidate_bytes = closure.read_regular_file_at(
                    candidate_descriptor,
                    "clean_candidate_manifest.json",
                    label="Clean-candidate manifest",
                )
                candidate = control_plane_common.read_json_bytes(
                    candidate_bytes, label="Clean-candidate manifest"
                )
                freeze_id = control_plane_common._clean_candidate_text(
                    candidate, "freeze_id"
                )
                try:
                    parsed_freeze_id = uuid.UUID(freeze_id)
                except ValueError as error:
                    raise ValueError("Clean-candidate freeze ID is malformed") from error
                _require(
                    str(parsed_freeze_id) == freeze_id
                    and candidate_path.parent.name == freeze_id,
                    "clean-candidate freeze ID or path is not canonical",
                )
                control_plane_common.utc_datetime(
                    candidate.get("created_utc"), field="clean_candidate.created_utc"
                )
                source = control_plane_common._clean_candidate_mapping(
                    candidate, "source"
                )
                build = control_plane_common._clean_candidate_mapping(candidate, "build")
                control_plane_common._require_exact_clean_candidate_layout_path(
                    source, "archive_path", candidate_path.parent / "source.tar"
                )
                control_plane_common._require_exact_clean_candidate_layout_path(
                    source, "commit_path", candidate_path.parent / "source.commit"
                )
                control_plane_common._require_exact_clean_candidate_layout_path(
                    build, "profile_path", candidate_path.parent / "build_profile.json"
                )
                control_plane_common._require_exact_clean_candidate_layout_path(
                    build,
                    "profile_receipt_path",
                    candidate_path.parent / "profile_receipt.json",
                )
                control_plane_common._require_exact_clean_candidate_layout_path(
                    build, "executable_path", candidate_path.parent / "athena"
                )
                expected = {
                    "athena",
                    "build_provenance",
                    "build_profile.json",
                    "clean_candidate_manifest.json",
                    "profile_receipt.json",
                    "source.commit",
                    "source.tar",
                }
                if source.get("submodules"):
                    expected.add("submodules")
                closure.require_exact_entries(
                    candidate_descriptor,
                    expected,
                    label="Clean-candidate directory",
                )
                source_archive = closure.read_regular_file_at(
                    candidate_descriptor,
                    "source.tar",
                    label="Clean-candidate source archive",
                )
                source_commit = closure.read_regular_file_at(
                    candidate_descriptor,
                    "source.commit",
                    label="Clean-candidate source commit object",
                )
                submodule_archives, submodule_commits = (
                    control_plane_common._read_clean_candidate_submodules(
                        source,
                        candidate_dir=candidate_path.parent,
                        candidate_descriptor=candidate_descriptor,
                        closure=closure,
                    )
                )
                build_profile = closure.read_regular_file_at(
                    candidate_descriptor,
                    "build_profile.json",
                    label="Clean-candidate build profile",
                )
                build_profile_receipt = closure.read_regular_file_at(
                    candidate_descriptor,
                    "profile_receipt.json",
                    label="Clean-candidate build-profile receipt",
                )
                executable = closure.read_regular_file_at(
                    candidate_descriptor,
                    "athena",
                    label="Clean-candidate executable",
                )
                provenance_descriptor = closure.open_directory_at(
                    candidate_descriptor,
                    "build_provenance",
                    label="Frozen build provenance directory",
                )
                provenance_expected = set(
                    control_plane_common.BUILD_PROVENANCE_FILENAMES.values()
                )
                closure.require_exact_entries(
                    provenance_descriptor,
                    provenance_expected,
                    label="Frozen build provenance directory",
                )
                build_provenance = {
                    label: closure.read_regular_file_at(
                        provenance_descriptor,
                        filename,
                        label=f"Frozen build provenance {label}",
                    )
                    for label, filename in (
                        control_plane_common.BUILD_PROVENANCE_FILENAMES.items()
                    )
                }
                closure.require_exact_entries(
                    provenance_descriptor,
                    provenance_expected,
                    label="Frozen build provenance directory",
                )
                closure.require_exact_entries(
                    candidate_descriptor,
                    expected,
                    label="Clean-candidate directory",
                )
                closure.require_same()
                candidate_root_ancestry.require_same()
                closure.require_same()
                return {
                    "candidate_manifest_path": candidate_path,
                    "candidate_manifest_bytes": candidate_bytes,
                    "candidate": candidate,
                    "source_archive": source_archive,
                    "source_commit": source_commit,
                    "submodule_archives": submodule_archives,
                    "submodule_commits": submodule_commits,
                    "build_profile": build_profile,
                    "build_profile_receipt": build_profile_receipt,
                    "build_provenance": build_provenance,
                    "executable": executable,
                }
            finally:
                closure.close(close_ancestry=candidate_root_ancestry.close)
    except PressureSelectionPublicationError:
        raise
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "authorized clean-candidate tree failed bounded capture"
        ) from error


class _CandidateSourceLoader(importlib.abc.Loader):
    def __init__(self, source: bytes, filename: str) -> None:
        self.source = source
        self.filename = filename

    def create_module(self, spec: object) -> ModuleType | None:
        return None

    def exec_module(self, module: ModuleType) -> None:
        module.__file__ = self.filename
        exec(
            compile(self.source, self.filename, "exec", dont_inherit=True),
            module.__dict__,
        )


class _CandidateSourceFinder(importlib.abc.MetaPathFinder):
    def __init__(
        self,
        package_name: str,
        archive_files: Mapping[str, bytes],
        source_root: Path,
    ) -> None:
        self.package_name = package_name
        self.sources: dict[str, tuple[bytes, str]] = {}
        self.blocked_top_level: set[str] = set()
        for relative in archive_files:
            parts = Path(relative).parts
            if (
                parts[:2] == ("tst", "publication")
                and relative.endswith(".py")
            ):
                self.blocked_top_level.add(parts[-1][:-3])
        for relative in (
            pressure_selection.pressure_review_packet_verifier.PRESSURE_REANALYSIS_SOURCE_PATHS
        ):
            _require(
                relative in archive_files,
                f"active-candidate reanalysis archive omits {relative}",
            )
            parts = Path(relative).parts
            _require(
                parts[:2] == ("tst", "publication") and relative.endswith(".py"),
                f"reanalysis source path is unsupported: {relative}",
            )
            module_parts = list(parts[2:])
            module_parts[-1] = module_parts[-1][:-3]
            aliases = {
                ".".join((package_name, *module_parts)),
                ".".join(("tst", "publication", *module_parts)),
                module_parts[-1],
            }
            if module_parts[0] == "frontier_control_plane":
                aliases.add(".".join(module_parts))
            for fullname in aliases:
                self.sources[fullname] = (
                    archive_files[relative],
                    str(source_root / relative),
                )

    def restricts(self, fullname: str) -> bool:
        return (
            fullname == "tst"
            or fullname.startswith("tst.")
            or fullname == "frontier_control_plane"
            or fullname.startswith("frontier_control_plane.")
            or fullname in self.blocked_top_level
        )

    def find_spec(
        self,
        fullname: str,
        path: object = None,
        target: ModuleType | None = None,
    ) -> object:
        del path, target
        source = self.sources.get(fullname)
        if source is None:
            if self.restricts(fullname):
                raise ImportError(
                    f"active-candidate reanalysis import is outside the authenticated closure: {fullname}"
                )
            return None
        payload, filename = source
        return importlib.util.spec_from_loader(
            fullname,
            _CandidateSourceLoader(payload, filename),
            origin=filename,
        )


def _candidate_package(name: str, path: str) -> ModuleType:
    module = ModuleType(name)
    module.__file__ = path
    module.__package__ = name
    module.__path__ = [path]  # type: ignore[attr-defined]
    return module


def _execute_candidate_reanalysis(archive_files: Mapping[str, bytes]) -> dict[str, str]:
    """Execute the reanalysis from active-candidate archive bytes only."""
    with tempfile.TemporaryDirectory(prefix="q011-stage4-candidate-source-") as directory:
        source_root = Path(directory)
        for relative, payload in archive_files.items():
            path = source_root / relative
            _require(
                path != source_root
                and source_root in path.parents
                and path == Path(os.path.abspath(path)),
                f"active-candidate archive member path is unsafe: {relative}",
            )
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
            path.chmod(0o400)
        package_name = f"_q011_stage4_candidate_{uuid.uuid4().hex}"
        package_path = str(source_root / "tst" / "publication")
        control_plane_name = f"{package_name}.frontier_control_plane"
        finder = _CandidateSourceFinder(package_name, archive_files, source_root)
        saved_modules = {
            name: module
            for name, module in tuple(sys.modules.items())
            if finder.restricts(name)
        }
        for name in saved_modules:
            del sys.modules[name]
        temporary_packages = {
            "tst": _candidate_package("tst", str(source_root / "tst")),
            "tst.publication": _candidate_package("tst.publication", package_path),
            "tst.publication.frontier_control_plane": _candidate_package(
                "tst.publication.frontier_control_plane",
                f"{package_path}/frontier_control_plane",
            ),
            "frontier_control_plane": _candidate_package(
                "frontier_control_plane", f"{package_path}/frontier_control_plane"
            ),
            package_name: _candidate_package(package_name, package_path),
            control_plane_name: _candidate_package(
                control_plane_name, f"{package_path}/frontier_control_plane"
            ),
        }
        sys.modules.update(temporary_packages)
        sys.modules["tst"].publication = sys.modules["tst.publication"]  # type: ignore[attr-defined]
        sys.modules["tst.publication"].frontier_control_plane = sys.modules[  # type: ignore[attr-defined]
            "tst.publication.frontier_control_plane"
        ]
        sys.meta_path.insert(0, finder)
        try:
            historical = importlib.import_module(
                f"{package_name}.q011_section54_historical_pressure_pilot_consumer"
            )
            result = historical.consume_exact_historical_production_pressure_pilot()
            _require(
                type(result) is dict,
                "active-candidate authoritative reanalysis returned a malformed result",
            )
            return copy.deepcopy(result)
        except PressureSelectionPublicationError:
            raise
        except Exception as error:
            raise PressureSelectionPublicationError(
                "active-candidate authoritative historical pressure-pilot recomputation failed"
            ) from error
        finally:
            sys.meta_path.remove(finder)
            for name in tuple(sys.modules):
                if (
                    name == package_name
                    or name.startswith(f"{package_name}.")
                    or finder.restricts(name)
                ):
                    del sys.modules[name]
            sys.modules.update(saved_modules)


def _pressure_gate_directory_name(
    sealed_utc: str,
    role: str,
    operator_id: str,
) -> str:
    timestamp = datetime.strptime(sealed_utc, "%Y-%m-%dT%H:%M:%SZ").strftime(
        "%Y%m%dT%H%M%SZ"
    )
    return f"{timestamp}-q011-section54-pressure-{role}-{operator_id}"


def _pressure_publication_evidence(
    *,
    authorized_pic_root: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    verifier = pressure_selection.pressure_review_packet_verifier
    aggregate_binding = dict(
        pressure_selection.historical_pressure_pilot_consumer.AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING
    )
    packet_binding = dict(
        pressure_selection.historical_pressure_pilot_consumer.AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING
    )
    try:
        packet = verifier.consume_published_pressure_pilot_review_packet(
            packet_binding["path"],
            aggregate_receipt_binding=aggregate_binding,
            authorized_pic_root=authorized_pic_root,
        )
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "published pressure-review packet failed immutable verification"
        ) from error
    aggregate_bundle = packet.get("aggregate_bundle")
    aggregate_analysis = packet.get("aggregate_analysis")
    _require(
        type(aggregate_bundle) is dict
        and set(aggregate_bundle) == {"path", "manifest_sha256"}
        and type(aggregate_analysis) is dict
        and set(aggregate_analysis) == {"path", "sha256"},
        "published pressure-review packet evidence tuple drifted",
    )
    return packet, {
        "published_pressure_pilot_receipt": aggregate_binding,
        "published_pressure_pilot_review_packet_receipt": packet_binding,
        "pilot_bundle_manifest_sha256": aggregate_bundle["manifest_sha256"],
        "aggregate_pilot_analysis_sha256": aggregate_analysis["sha256"],
    }


def _consume_reanalysis_binding(
    binding: object,
    *,
    authorized_pic_root: Path,
    now: datetime | None = None,
) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    packet, evidence = _pressure_publication_evidence(
        authorized_pic_root=authorized_pic_root
    )
    verifier = pressure_selection.pressure_review_packet_verifier
    try:
        reanalysis = verifier.consume_sealed_pressure_reanalysis_attestation(
            binding,
            aggregate_receipt_binding=evidence["published_pressure_pilot_receipt"],
            packet_receipt_binding=evidence[
                "published_pressure_pilot_review_packet_receipt"
            ],
            pilot_bundle_manifest_sha256=evidence["pilot_bundle_manifest_sha256"],
            aggregate_pilot_analysis_sha256=evidence[
                "aggregate_pilot_analysis_sha256"
            ],
            authorized_pic_root=authorized_pic_root,
            now=now,
        )
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "sealed authoritative reanalysis failed immutable verification"
        ) from error
    return packet, evidence, reanalysis


def _case_descriptors(packet: Mapping[str, object]) -> list[dict[str, object]]:
    aggregate_receipt = packet.get("aggregate_receipt")
    _require(type(aggregate_receipt) is dict, "verified aggregate receipt is absent")
    raw_cases = aggregate_receipt.get("raw_cases")
    _require(type(raw_cases) is list, "verified aggregate raw cases are absent")
    descriptors_by_case = {
        item.get("case_id"): item.get("descriptor_sha256")
        for item in raw_cases
        if type(item) is dict
    }
    descriptors = [
        {
            "case_id": case_id,
            "problem_ps_p0": problem_ps_p0,
            "descriptor_sha256": descriptors_by_case.get(case_id),
        }
        for case_id, problem_ps_p0 in pressure_selection.REGISTERED_CASES
    ]
    for descriptor in descriptors:
        _lowercase_sha256(
            descriptor["descriptor_sha256"],
            f"{descriptor['case_id']} descriptor SHA-256",
        )
    return descriptors


def _validate_human_decision(
    value: object,
    *,
    reanalysis_verification: Mapping[str, object],
) -> dict[str, object]:
    expected_keys = {
        "schema_version",
        "record_type",
        "reviewer_id",
        "reviewed_utc",
        "rationale",
        "reviewer_statement",
        "authoritative_reanalysis_attestation",
        "stage4_preparation_attestation",
        "selected_case",
    }
    _require(
        type(value) is dict and set(value) == expected_keys,
        "human pressure-selection decision schema drifted",
    )
    reviewed_utc = _require_canonical_utc(
        value["reviewed_utc"], "human pressure-selection decision reviewed_utc"
    )
    reanalysis_binding = _binding(
        reanalysis_verification.get("binding"),
        "verified authoritative reanalysis attestation",
    )
    reanalysis_sealed_utc = _require_canonical_utc(
        reanalysis_verification.get("sealed_utc"),
        "verified authoritative reanalysis sealed_utc",
    )
    _require(
        type(value["schema_version"]) is int
        and value["schema_version"] == 1
        and value["record_type"] == HUMAN_DECISION_RECORD_TYPE
        and value["reviewer_id"] == REVIEWED_REVIEWER_ID
        and value["rationale"] == REVIEWED_RATIONALE
        and value["reviewer_statement"] == HUMAN_DECISION_STATEMENT
        and value["selected_case"] == REVIEWED_SELECTED_CASE
        and value["authoritative_reanalysis_attestation"]
        == _binding(reanalysis_binding, "authoritative reanalysis attestation"),
        "human pressure-selection decision differs from the required reviewed choice",
    )
    _require(
        reanalysis_sealed_utc < reviewed_utc,
        "human pressure-selection decision must be strictly later than sealed reanalysis",
    )
    normalized = copy.deepcopy(value)
    normalized["reviewed_utc"] = reviewed_utc
    normalized["authoritative_reanalysis_attestation"] = _binding(
        value["authoritative_reanalysis_attestation"],
        "human decision authoritative reanalysis attestation",
    )
    normalized["stage4_preparation_attestation"] = _binding(
        value["stage4_preparation_attestation"],
        "human decision Stage-4 preparation attestation",
    )
    return normalized


def _load_human_decision(
    path: Path,
    *,
    authorized_pic_root: Path,
) -> tuple[dict[str, object], dict[str, str]]:
    decision_root = pilot_publisher._canonical_existing_directory(
        authorized_pic_root / PRESSURE_GATE_HUMAN_DECISION_ROOT_NAME,
        "human pressure-selection decision root",
    )
    lexical = Path(os.path.abspath(path))
    _require(
        lexical == path
        and lexical.parent == decision_root
        and lexical.name not in {"", ".", ".."},
        "human pressure-selection decision must be one direct canonical decision-root child",
    )
    descriptor = pilot_publisher._open_absolute_directory(decision_root)
    try:
        metadata = os.fstat(descriptor)
        _require(
            stat.S_IMODE(metadata.st_mode) == 0o700
            and metadata.st_uid == os.geteuid(),
            "human pressure-selection decision root must be private",
        )
        payload, _ = _read_exact_mode_regular_at(
            descriptor,
            lexical.name,
            mode=0o400,
            label="human pressure-selection decision",
        )
    finally:
        os.close(descriptor)
    return (
        _decode_canonical_object(payload, "human pressure-selection decision"),
        {"path": str(lexical), "sha256": _sha256(payload)},
    )


def prepare_pressure_reanalysis(
    *,
    reanalysis_operator_id: str,
    expected_git_commit: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    now: datetime | None = None,
) -> dict[str, object]:
    """Seal machine reanalysis only; never create a reviewer attestation or receipt."""
    operator_id = _operator_id(reanalysis_operator_id, "reanalysis operator ID")
    pic_root = Path(os.path.abspath(authorized_pic_root))
    project_root = Path(os.path.abspath(authorized_project_home_root))
    _require(
        pic_root == AUTHORIZED_PIC_ROOT and project_root == AUTHORIZED_PROJECT_HOME_ROOT,
        "production pressure-selection preparation requires the authorized roots",
    )
    publisher_source_authentication = _runtime_publisher_source_authentication(
        expected_git_commit
    )
    archive_root, archive_descriptor = _open_or_create_private_root(
        pic_root,
        PRESSURE_GATE_ATTESTATION_ROOT_NAME,
        "pressure-gate attestation archive root",
    )
    decision_root, decision_descriptor = _open_or_create_private_root(
        pic_root,
        PRESSURE_GATE_HUMAN_DECISION_ROOT_NAME,
        "human pressure-selection decision root",
    )
    acceptance_root = pilot_publisher._publication_acceptance_root(pic_root)
    acceptance_descriptor = pilot_publisher._open_absolute_directory(acceptance_root)
    transaction_anchor = pilot_publisher._canonical_existing_directory(
        control_plane_common.stable_serialization_anchor(pic_root),
        "stable publication transaction anchor",
    )
    transaction_descriptor = pilot_publisher._open_absolute_directory(transaction_anchor)
    try:
        pilot_publisher._lock_publication_transaction(
            transaction_descriptor, acceptance_descriptor
        )
        _require(
            os.listdir(decision_descriptor) == [],
            "human pressure-selection decision root must be empty before reanalysis",
        )
        capture = _capture_live_controller_state(
            authorized_pic_root=pic_root,
            authorized_project_home_root=project_root,
        )
        state = capture["controller_state"]
        assert isinstance(state, dict)
        source_authorization, archive_files = _candidate_reanalysis_source_authorization(
            state, authorized_pic_root=pic_root
        )
        recomputed_utc = _canonical_utc(now)
        result = _execute_candidate_reanalysis(archive_files)
        verifier = pressure_selection.pressure_review_packet_verifier
        _packet, evidence = _pressure_publication_evidence(
            authorized_pic_root=pic_root
        )
        sealed_utc = _canonical_utc(now)
        reanalysis_attestation = {
            "schema_version": 1,
            "record_type": verifier.PRESSURE_REANALYSIS_RECORD_TYPE,
            "qualification_effect": verifier.PRESSURE_REANALYSIS_QUALIFICATION_EFFECT,
            "operator_id": operator_id,
            "recomputed_utc": recomputed_utc,
            "sealed_utc": sealed_utc,
            "operator_statement": verifier.PRESSURE_REANALYSIS_OPERATOR_STATEMENT,
            "evidence": evidence,
            "source_authorization": source_authorization,
            "result": result,
        }
        reanalysis_binding = _write_single_file_attestation_at(
            archive_root,
            archive_descriptor,
            _pressure_gate_directory_name(sealed_utc, "reanalysis", operator_id),
            reanalysis_attestation,
        )
        reanalysis = verifier.consume_sealed_pressure_reanalysis_attestation(
            reanalysis_binding,
            aggregate_receipt_binding=evidence["published_pressure_pilot_receipt"],
            packet_receipt_binding=evidence[
                "published_pressure_pilot_review_packet_receipt"
            ],
            pilot_bundle_manifest_sha256=evidence["pilot_bundle_manifest_sha256"],
            aggregate_pilot_analysis_sha256=evidence[
                "aggregate_pilot_analysis_sha256"
            ],
            authorized_pic_root=pic_root,
            expected_result=result,
            now=now,
        )
        stage4_preparation = {
            "schema_version": 1,
            "record_type": STAGE4_PREPARATION_RECORD_TYPE,
            "qualification_effect": QUALIFICATION_EFFECT,
            "operator_id": operator_id,
            "sealed_utc": sealed_utc,
            "publisher_source_authentication": publisher_source_authentication,
            "authoritative_reanalysis_attestation": reanalysis_binding,
        }
        stage4_preparation_binding = _write_single_file_attestation_at(
            archive_root,
            archive_descriptor,
            _pressure_gate_directory_name(
                sealed_utc, "stage4-preparation", operator_id
            ),
            stage4_preparation,
        )
        consumed_preparation = _consume_stage4_preparation_attestation(
            stage4_preparation_binding, authorized_pic_root=pic_root
        )
        _require(
            consumed_preparation["attestation"] == stage4_preparation,
            "Stage-4 preparation attestation drifted after sealing",
        )
        return {
            "authoritative_reanalysis_attestation": reanalysis_binding,
            "stage4_preparation_attestation": stage4_preparation_binding,
            "human_decision_root": str(decision_root),
            "required_human_confirmation": {
                "schema_version": 1,
                "record_type": HUMAN_DECISION_RECORD_TYPE,
                "reviewer_id": REVIEWED_REVIEWER_ID,
                "selected_case": dict(REVIEWED_SELECTED_CASE),
                "rationale": REVIEWED_RATIONALE,
                "reviewer_statement": HUMAN_DECISION_STATEMENT,
            },
            "publisher_source_authentication": publisher_source_authentication,
            "qualification_effect": verifier.PRESSURE_REANALYSIS_QUALIFICATION_EFFECT,
        }
    except PressureSelectionPublicationError:
        raise
    except BaseException as error:
        raise PressureSelectionPublicationError(
            f"pressure-selection preparation failed closed: {error}"
        ) from error
    finally:
        pilot_publisher._close_descriptors(
            (
                transaction_descriptor,
                acceptance_descriptor,
                decision_descriptor,
                archive_descriptor,
            )
        )


def seal_human_pressure_selection(
    *,
    human_decision_path: Path,
    expected_git_commit: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    now: datetime | None = None,
) -> dict[str, object]:
    """Seal one explicit post-reanalysis human decision and create a candidate receipt."""
    pic_root = Path(os.path.abspath(authorized_pic_root))
    project_root = Path(os.path.abspath(authorized_project_home_root))
    _require(
        pic_root == AUTHORIZED_PIC_ROOT and project_root == AUTHORIZED_PROJECT_HOME_ROOT,
        "production human pressure selection requires the authorized roots",
    )
    publisher_source_authentication = _runtime_publisher_source_authentication(
        expected_git_commit
    )
    decision, decision_binding = _load_human_decision(
        human_decision_path, authorized_pic_root=pic_root
    )
    reanalysis_binding = _binding(
        decision.get("authoritative_reanalysis_attestation"),
        "human decision authoritative reanalysis attestation",
    )
    stage4_preparation_binding = _binding(
        decision.get("stage4_preparation_attestation"),
        "human decision Stage-4 preparation attestation",
    )
    archive_root, archive_descriptor = _open_or_create_private_root(
        pic_root,
        PRESSURE_GATE_ATTESTATION_ROOT_NAME,
        "pressure-gate attestation archive root",
    )
    candidate_root, candidate_descriptor = _open_or_create_private_root(
        pic_root,
        PRESSURE_GATE_CANDIDATE_ROOT_NAME,
        "pressure-gate candidate root",
    )
    acceptance_root = pilot_publisher._publication_acceptance_root(pic_root)
    acceptance_descriptor = pilot_publisher._open_absolute_directory(acceptance_root)
    transaction_anchor = pilot_publisher._canonical_existing_directory(
        control_plane_common.stable_serialization_anchor(pic_root),
        "stable publication transaction anchor",
    )
    transaction_descriptor = pilot_publisher._open_absolute_directory(transaction_anchor)
    try:
        pilot_publisher._lock_publication_transaction(
            transaction_descriptor, acceptance_descriptor
        )
        _capture_live_controller_state(
            authorized_pic_root=pic_root,
            authorized_project_home_root=project_root,
        )
        packet, evidence, reanalysis = _consume_reanalysis_binding(
            reanalysis_binding, authorized_pic_root=pic_root, now=now
        )
        stage4_preparation = _consume_stage4_preparation_attestation(
            stage4_preparation_binding, authorized_pic_root=pic_root
        )
        preparation_attestation = stage4_preparation["attestation"]
        _require(
            type(preparation_attestation) is dict
            and preparation_attestation["publisher_source_authentication"]
            == publisher_source_authentication
            and preparation_attestation["authoritative_reanalysis_attestation"]
            == reanalysis_binding,
            "human sealing source differs from Stage-4 preparation source",
        )
        normalized_decision = _validate_human_decision(
            decision, reanalysis_verification=reanalysis
        )
        _require(
            normalized_decision["stage4_preparation_attestation"]
            == stage4_preparation_binding,
            "human decision binds another Stage-4 preparation attestation",
        )
        verifier = pressure_selection.pressure_review_packet_verifier
        sealed_utc = _canonical_utc(now)
        reviewer_attestation = {
            "schema_version": 1,
            "record_type": verifier.PRESSURE_REVIEWER_RECORD_TYPE,
            "qualification_effect": verifier.PRESSURE_REVIEWER_QUALIFICATION_EFFECT,
            "selection_method": pressure_selection.SELECTION_METHOD,
            "reviewer_id": normalized_decision["reviewer_id"],
            "reviewed_utc": normalized_decision["reviewed_utc"],
            "sealed_utc": sealed_utc,
            "rationale": normalized_decision["rationale"],
            "reviewer_statement": verifier.PRESSURE_REVIEWER_STATEMENT,
            "published_pressure_pilot_receipt": evidence[
                "published_pressure_pilot_receipt"
            ],
            "published_pressure_pilot_review_packet_receipt": evidence[
                "published_pressure_pilot_review_packet_receipt"
            ],
            "authoritative_reanalysis_attestation": reanalysis_binding,
            "selected_case": copy.deepcopy(normalized_decision["selected_case"]),
        }
        reviewer_binding = _write_single_file_attestation_at(
            archive_root,
            archive_descriptor,
            _pressure_gate_directory_name(
                sealed_utc, "selection", str(normalized_decision["reviewer_id"])
            ),
            reviewer_attestation,
        )
        verifier.consume_sealed_pressure_reviewer_attestation(
            reviewer_binding,
            aggregate_receipt_binding=evidence["published_pressure_pilot_receipt"],
            packet_receipt_binding=evidence[
                "published_pressure_pilot_review_packet_receipt"
            ],
            reanalysis_verification=reanalysis,
            selected_case=REVIEWED_SELECTED_CASE,
            authorized_pic_root=pic_root,
            now=now,
        )
        receipt = build_pressure_selection_receipt(
            published_pressure_pilot_receipt=evidence[
                "published_pressure_pilot_receipt"
            ],
            published_pressure_pilot_review_packet_receipt=evidence[
                "published_pressure_pilot_review_packet_receipt"
            ],
            pilot_bundle_manifest_sha256=str(
                evidence["pilot_bundle_manifest_sha256"]
            ),
            aggregate_pilot_analysis_sha256=str(
                evidence["aggregate_pilot_analysis_sha256"]
            ),
            case_descriptors=_case_descriptors(packet),
            authoritative_reanalysis_attestation=reanalysis_binding,
            reviewer_attestation=reviewer_binding,
        )
        normalized, _reanalysis, _reviewer = _validate_reviewed_selection(
            receipt, authorized_pic_root=pic_root
        )
        candidate_binding = _write_candidate_receipt_at(
            candidate_root, candidate_descriptor, normalized
        )
        candidate_authorization = _validate_candidate_publication_authorization(
            {
                "schema_version": 1,
                "record_type": CANDIDATE_AUTHORIZATION_RECORD_TYPE,
                "qualification_effect": QUALIFICATION_EFFECT,
                "sealed_utc": sealed_utc,
                "publisher_source_authentication": publisher_source_authentication,
                "stage4_preparation_attestation": stage4_preparation_binding,
                "human_decision": decision_binding,
                "candidate_pressure_selection_receipt": candidate_binding,
                "authoritative_reanalysis_attestation": reanalysis_binding,
                "reviewer_attestation": reviewer_binding,
            }
        )
        candidate_authorization_binding = _write_candidate_authorization_at(
            candidate_root, candidate_descriptor, candidate_authorization
        )
        return {
            "candidate_pressure_selection_receipt": candidate_binding,
            "candidate_publication_authorization": candidate_authorization_binding,
            "human_decision": decision_binding,
            "authoritative_reanalysis_attestation": reanalysis_binding,
            "reviewer_attestation": reviewer_binding,
            "selected_case": dict(REVIEWED_SELECTED_CASE),
            "reviewer_id": REVIEWED_REVIEWER_ID,
            "rationale": REVIEWED_RATIONALE,
            "publisher_source_authentication": publisher_source_authentication,
            "qualification_effect": QUALIFICATION_EFFECT,
        }
    except PressureSelectionPublicationError:
        raise
    except BaseException as error:
        raise PressureSelectionPublicationError(
            f"human pressure-selection sealing failed closed: {error}"
        ) from error
    finally:
        pilot_publisher._close_descriptors(
            (
                transaction_descriptor,
                acceptance_descriptor,
                candidate_descriptor,
                archive_descriptor,
            )
        )


def _validate_controller_state_attestation(value: object) -> dict[str, object]:
    top_keys = {
        "schema_version",
        "record_type",
        "qualification_effect",
        "operator_id",
        "captured_utc",
        "sealed_utc",
        "publisher_source_authentication",
        "controller_state",
        "source_binding",
    }
    _require(type(value) is dict and set(value) == top_keys, "controller-state attestation schema drifted")
    _require(
        type(value["schema_version"]) is int
        and value["schema_version"] == 1
        and value["record_type"] == CONTROLLER_STATE_RECORD_TYPE
        and value["qualification_effect"] == QUALIFICATION_EFFECT,
        "controller-state attestation identity drifted",
    )
    _operator_id(value["operator_id"], "controller-state operator ID")
    normalized_source_authentication = _validate_publisher_source_authentication(
        value["publisher_source_authentication"]
    )
    captured = _require_canonical_utc(value["captured_utc"], "controller-state captured_utc")
    sealed = _require_canonical_utc(value["sealed_utc"], "controller-state sealed_utc")
    _require(captured <= sealed, "controller-state attestation timestamps are out of order")
    state = value["controller_state"]
    state_keys = {
        "control_plane_version",
        "installed_control_planes",
        "active_policy",
        "active_promotion",
        "registered_science_slices",
        "frontier_admission_smoke",
        "science_submission_freeze",
        "clean_candidate",
        "mirrored_ledger_state",
        "pending_submission_marker",
        "manual_accounting_authorizations",
        "pending_manual_accounting_marker",
        "active_promotion_transaction",
    }
    _require(type(state) is dict and set(state) == state_keys, "controller-state snapshot schema drifted")
    version = _lowercase_sha256(state["control_plane_version"], "controller-state version")
    installed = state["installed_control_planes"]
    _require(
        type(installed) is list
        and len(installed) == 2
        and [record.get("root_role") for record in installed if type(record) is dict]
        == ["orion", "project_home"],
        "controller-state installed control-plane pair drifted",
    )
    for record in installed:
        _require(
            type(record) is dict
            and set(record) == {"root_role", "inventory_path", "inventory_sha256"}
            and type(record["inventory_path"]) is str
            and Path(record["inventory_path"]).is_absolute(),
            "controller-state installed control-plane binding drifted",
        )
        _lowercase_sha256(record["inventory_sha256"], "installed inventory SHA-256")
    for key in ("active_policy", "active_promotion"):
        record = state[key]
        _require(
            type(record) is dict
            and set(record) == {"orion_path", "project_home_path", "sha256"}
            and type(record["orion_path"]) is str
            and Path(record["orion_path"]).is_absolute()
            and type(record["project_home_path"]) is str
            and Path(record["project_home_path"]).is_absolute(),
            f"controller-state {key} binding drifted",
        )
        _lowercase_sha256(record["sha256"], f"controller-state {key} SHA-256")
    freeze = state["science_submission_freeze"]
    manual_accounting_authorizations = state["manual_accounting_authorizations"]
    mirrored_ledger_state = state["mirrored_ledger_state"]
    mirrored_ledger_keys = {
        "validator",
        "orion_ledger_path",
        "orion_mirror_receipts_path",
        "project_home_ledger_path",
        "record_count",
        "tail_event_sha256",
        "active_reservation_ids",
        "currently_reserved_node_hours",
        "cumulative_consumed_node_hours",
        "pending_submission_marker",
        "pending_manual_accounting_markers",
    }
    _require(
        state["registered_science_slices"] == []
        and state["frontier_admission_smoke"]
        == {"status": control_plane_common.CLOSED_ADMISSION_SMOKE_STATUS}
        and state["pending_submission_marker"] == "absent"
        and type(manual_accounting_authorizations) is list
        and state["pending_manual_accounting_marker"] == "absent"
        and state["active_promotion_transaction"] == "absent"
        and type(freeze) is dict
        and set(freeze)
        == {
            "status",
            "manifest_path",
            "manifest_sha256",
            "build_profile_control_plane_version",
        }
        and freeze.get("status") == control_plane_common.AUTHORIZED_CLEAN_CANDIDATE_FREEZE
        and freeze.get("build_profile_control_plane_version") == version,
        "controller-state snapshot is not launch-prohibited",
    )
    _require(
        type(mirrored_ledger_state) is dict
        and set(mirrored_ledger_state) == mirrored_ledger_keys
        and mirrored_ledger_state["validator"]
        == "validated_read_only_mirrored_state_snapshot"
        and all(
            type(mirrored_ledger_state[key]) is str
            and Path(mirrored_ledger_state[key]).is_absolute()
            for key in (
                "orion_ledger_path",
                "orion_mirror_receipts_path",
                "project_home_ledger_path",
            )
        )
        and type(mirrored_ledger_state["record_count"]) is int
        and mirrored_ledger_state["record_count"] > 0
        and mirrored_ledger_state["active_reservation_ids"] == []
        and type(mirrored_ledger_state["currently_reserved_node_hours"]) is float
        and mirrored_ledger_state["currently_reserved_node_hours"] == 0.0
        and type(mirrored_ledger_state["cumulative_consumed_node_hours"]) is float
        and mirrored_ledger_state["cumulative_consumed_node_hours"] >= 0.0
        and mirrored_ledger_state["pending_submission_marker"] == "absent"
        and mirrored_ledger_state["pending_manual_accounting_markers"]
        == {"orion": "absent", "project_home": "absent"},
        "controller-state mirrored ledger is not quiescent",
    )
    _lowercase_sha256(
        mirrored_ledger_state["tail_event_sha256"],
        "controller-state mirrored ledger tail SHA-256",
    )
    for authorization in manual_accounting_authorizations:
        _require(
            type(authorization) is dict
            and set(authorization)
            == {"authorization_id", "path", "project_home_path", "sha256"}
            and type(authorization["authorization_id"]) is str
            and bool(authorization["authorization_id"])
            and type(authorization["path"]) is str
            and Path(authorization["path"]).is_absolute()
            and type(authorization["project_home_path"]) is str
            and Path(authorization["project_home_path"]).is_absolute(),
            "controller-state historical accounting binding drifted",
        )
        _lowercase_sha256(
            authorization["sha256"],
            "controller-state historical accounting binding SHA-256",
        )
    _require(
        type(freeze["manifest_path"]) is str
        and Path(freeze["manifest_path"]).is_absolute(),
        "controller-state clean-candidate manifest path drifted",
    )
    _lowercase_sha256(freeze.get("manifest_sha256"), "controller-state clean-candidate SHA-256")
    candidate = state["clean_candidate"]
    _require(
        type(candidate) is dict
        and set(candidate)
        == {
            "manifest_path",
            "manifest_sha256",
            "git_commit",
            "source_archive_sha256",
        }
        and candidate["manifest_path"] == freeze["manifest_path"]
        and candidate["manifest_sha256"] == freeze["manifest_sha256"]
        and type(candidate["git_commit"]) is str
        and _GIT_COMMIT_PATTERN.fullmatch(candidate["git_commit"]) is not None,
        "controller-state clean-candidate binding drifted",
    )
    _lowercase_sha256(
        candidate["source_archive_sha256"],
        "controller-state clean-candidate source archive SHA-256",
    )
    source = value["source_binding"]
    _require(
        type(source) is dict
        and set(source)
        == {
            "authoritative_reanalysis_attestation",
            "git_commit",
            "source_archive_sha256",
            "reanalysis_source_closure_sha256",
            "clean_candidate_manifest",
            "candidate_publication_authorization",
        }
        and type(source["git_commit"]) is str
        and _GIT_COMMIT_PATTERN.fullmatch(source["git_commit"]) is not None,
        "controller-state source binding drifted",
    )
    _binding(source["authoritative_reanalysis_attestation"], "authoritative reanalysis attestation")
    _binding(source["clean_candidate_manifest"], "clean-candidate manifest")
    candidate_authorization = _validate_candidate_publication_authorization(
        source["candidate_publication_authorization"]
    )
    _lowercase_sha256(source["source_archive_sha256"], "source archive SHA-256")
    _lowercase_sha256(source["reanalysis_source_closure_sha256"], "source closure SHA-256")
    _require(
        source["clean_candidate_manifest"]
        == {
            "path": candidate["manifest_path"],
            "sha256": candidate["manifest_sha256"],
        }
        and candidate_authorization["publisher_source_authentication"]
        == normalized_source_authentication
        and candidate_authorization["authoritative_reanalysis_attestation"]
        == source["authoritative_reanalysis_attestation"]
        and source["git_commit"] == candidate["git_commit"]
        and source["source_archive_sha256"] == candidate["source_archive_sha256"],
        "controller-state clean-candidate source binding drifted",
    )
    normalized = copy.deepcopy(value)
    normalized["publisher_source_authentication"] = normalized_source_authentication
    return normalized


def consume_sealed_controller_state_attestation(
    attestation_binding: object,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    """Verify one immutable controller-state attestation without consulting live policy."""
    binding = _binding(attestation_binding, "controller-state attestation")
    pic_root = Path(os.path.abspath(authorized_pic_root))
    archive_root = pic_root / PRESSURE_GATE_ATTESTATION_ROOT_NAME
    path = Path(binding["path"])
    _require(
        path.name == "attestation.json"
        and path.parent.parent == archive_root
        and path == Path(os.path.abspath(path)),
        "controller-state attestation path is outside the archive",
    )
    try:
        canonical = control_plane_common.require_canonical_path_below(path, archive_root)
    except (OSError, ValueError) as error:
        raise PressureSelectionPublicationError(
            "controller-state attestation path is not canonical"
        ) from error
    directory = canonical.parent
    descriptor = os.open(directory, _DIRECTORY_FLAGS)
    try:
        metadata = os.fstat(descriptor)
        _require(
            stat.S_IMODE(metadata.st_mode) == 0o500,
            "controller-state attestation directory mode must be 0500",
        )
        _require(
            set(os.listdir(descriptor))
            == {"attestation.json", "active_policy.json", "active_promotion.json"},
            "controller-state attestation member closure drifted",
        )
        attestation_payload, _ = _read_exact_mode_regular_at(
            descriptor,
            "attestation.json",
            mode=0o400,
            label="controller-state attestation payload",
        )
        policy_payload, _ = _read_exact_mode_regular_at(
            descriptor,
            "active_policy.json",
            mode=0o400,
            label="controller-state active-policy snapshot",
        )
        promotion_payload, _ = _read_exact_mode_regular_at(
            descriptor,
            "active_promotion.json",
            mode=0o400,
            label="controller-state active-promotion snapshot",
        )
        _require(
            _sha256(attestation_payload) == binding["sha256"],
            "controller-state attestation binding hash mismatch",
        )
        attestation = _validate_controller_state_attestation(
            _decode_canonical_object(attestation_payload, "controller-state attestation")
        )
        _require(
            directory.name
            == _controller_attestation_directory_name(
                str(attestation["sealed_utc"]), str(attestation["operator_id"])
            ),
            "controller-state attestation directory name drifted",
        )
        state = attestation["controller_state"]
        assert isinstance(state, dict)
        _require(
            _sha256(policy_payload) == state["active_policy"]["sha256"]
            and _sha256(promotion_payload) == state["active_promotion"]["sha256"],
            "controller-state policy or promotion snapshot hash drifted",
        )
        policy = _decode_canonical_object(
            policy_payload, "controller-state active-policy snapshot"
        )
        promotion = _decode_canonical_object(
            promotion_payload, "controller-state active-promotion snapshot"
        )
        _require(
            promotion.get("control_plane_version") == state["control_plane_version"]
            and promotion.get("policy_sha256") == state["active_policy"]["sha256"]
            and policy.get("registered_science_slices")
            == state["registered_science_slices"]
            and policy.get("frontier_admission_smoke")
            == state["frontier_admission_smoke"]
            and policy.get("science_submission_freeze")
            == state["science_submission_freeze"]
            and type(policy.get("olcf_side_storage")) is dict
            and policy["olcf_side_storage"].get("manual_accounting_authorizations")
            == state["manual_accounting_authorizations"],
            "controller-state policy or promotion snapshot differs from the attestation",
        )
        _require(
            set(os.listdir(descriptor))
            == {"attestation.json", "active_policy.json", "active_promotion.json"},
            "controller-state attestation changed during verification",
        )
    finally:
        os.close(descriptor)
    return {"binding": binding, "attestation": attestation}


def _require_capture_matches_attestation(
    capture: Mapping[str, object],
    consumed: Mapping[str, object],
) -> None:
    attestation = consumed["attestation"]
    _require(
        type(attestation) is dict
        and attestation["controller_state"] == capture["controller_state"],
        "live controller state differs from the sealed controller-state attestation",
    )


def _verify_published_pressure_selection(
    receipt_path: str | Path,
    *,
    authorized_pic_root: Path,
    allow_publication_guard: bool,
    require_publication_seal: bool,
    expected_controller_state_attestation: object | None = None,
) -> dict[str, object]:
    pic_root, publication_root = pilot_publisher._publication_root(authorized_pic_root)
    acceptance_root = pilot_publisher._publication_acceptance_root(pic_root)
    target = pilot_publisher._direct_publication_target(
        receipt_path, publication_root, "pressure-selection receipt"
    )
    _require(target.name == CANONICAL_RECEIPT_NAME, "pressure-selection receipt name drifted")
    publication_descriptor = pilot_publisher._open_absolute_directory(publication_root)
    acceptance_descriptor = pilot_publisher._open_absolute_directory(acceptance_root)
    try:
        pilot_publisher._require_same_directory(
            publication_root, publication_descriptor, "authorized PIC publication root"
        )
        pilot_publisher._require_same_directory(
            acceptance_root,
            acceptance_descriptor,
            "authorized PIC publication acceptance root",
        )
        if not allow_publication_guard:
            pilot_publisher._require_publication_guard_absent_at(
                publication_descriptor, target.name, "pressure-selection receipt"
            )
        receipt_payload, receipt_identity = _read_exact_mode_regular_at(
            publication_descriptor,
            target.name,
            mode=0o444,
            label="pressure-selection receipt",
        )
        receipt = _decode_canonical_object(receipt_payload, "pressure-selection receipt")
        normalized, reanalysis, reviewer = _validate_reviewed_selection(
            receipt, authorized_pic_root=pic_root
        )
        _require(
            receipt_payload == pressure_selection.canonical_json_bytes(normalized),
            "pressure-selection receipt bytes drifted from schema-v3 normalization",
        )
        seal = None
        if require_publication_seal:
            seal, _ = _read_selection_success_seal_at(
                acceptance_descriptor,
                publication_descriptor,
                target.name,
                receipt_payload,
                receipt_identity,
            )
            controller_binding = seal["controller_state_attestation"]
        else:
            _require(
                expected_controller_state_attestation is not None,
                "internal guarded verification requires a controller-state attestation",
            )
            controller_binding = _binding(
                expected_controller_state_attestation, "controller-state attestation"
            )
        controller = consume_sealed_controller_state_attestation(
            controller_binding, authorized_pic_root=pic_root
        )
        controller_attestation = controller["attestation"]
        assert isinstance(controller_attestation, dict)
        source_binding = controller_attestation["source_binding"]
        assert isinstance(source_binding, dict)
        candidate_authorization = _validate_candidate_publication_authorization(
            source_binding["candidate_publication_authorization"]
        )
        _require(
            _authorize_candidate_publication(
                candidate_authorization,
                receipt_payload=receipt_payload,
                receipt=normalized,
                reanalysis_verification=reanalysis,
                reviewer_verification=reviewer,
                publisher_source_authentication=controller_attestation[
                    "publisher_source_authentication"
                ],
                authorized_pic_root=pic_root,
            )
            == candidate_authorization,
            "published pressure-selection candidate authorization drifted",
        )
        preparation = _consume_stage4_preparation_attestation(
            candidate_authorization["stage4_preparation_attestation"],
            authorized_pic_root=pic_root,
        )
        preparation_attestation = preparation["attestation"]
        _require(
            type(preparation_attestation) is dict
            and preparation_attestation["publisher_source_authentication"]
            == controller_attestation["publisher_source_authentication"]
            and preparation_attestation["authoritative_reanalysis_attestation"]
            == normalized["authoritative_reanalysis_attestation"],
            "published pressure-selection source differs from Stage-4 preparation source",
        )
        if seal is not None:
            _require(
                seal["sealed_utc"]
                >= controller["attestation"]["sealed_utc"],
                "pressure-selection success seal predates its controller-state attestation",
            )
        _require(
            controller["attestation"]["source_binding"][
                "authoritative_reanalysis_attestation"
            ]
            == normalized["authoritative_reanalysis_attestation"],
            "controller-state attestation binds a different reanalysis",
        )
        repeated_payload, repeated_identity = _read_exact_mode_regular_at(
            publication_descriptor,
            target.name,
            mode=0o444,
            label="pressure-selection receipt",
        )
        _require(
            repeated_payload == receipt_payload and repeated_identity == receipt_identity,
            "pressure-selection receipt changed during verification",
        )
        if not allow_publication_guard:
            pilot_publisher._require_publication_guard_absent_at(
                publication_descriptor, target.name, "pressure-selection receipt"
            )
        return {
            "receipt_binding": {"path": str(target), "sha256": _sha256(receipt_payload)},
            "receipt": normalized,
            "selected_case": copy.deepcopy(normalized["selected_case"]),
            "reviewer_id": reviewer["reviewer_id"],
            "rationale": reviewer["rationale"],
            "controller_state_attestation": controller,
            "success_seal": seal,
            "qualification_effect": QUALIFICATION_EFFECT,
        }
    finally:
        os.close(acceptance_descriptor)
        os.close(publication_descriptor)


def verify_published_pressure_selection_receipt(
    receipt_path: str | Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    """Re-audit one externally consumable pressure-selection receipt."""
    return _verify_published_pressure_selection(
        receipt_path,
        authorized_pic_root=authorized_pic_root,
        allow_publication_guard=False,
        require_publication_seal=True,
    )


def verify_published_pressure_selection_live_state(
    receipt_path: str | Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    """Re-audit the receipt and require the paired live state to match its seal."""
    verified = verify_published_pressure_selection_receipt(
        receipt_path, authorized_pic_root=authorized_pic_root
    )
    capture = _capture_live_controller_state(
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    _require_capture_matches_attestation(
        capture, verified["controller_state_attestation"]
    )
    result = copy.deepcopy(verified)
    result["live_controller_state_verification"] = {
        "status": "matches_sealed_controller_state_attestation",
        "authorized_pic_root": str(Path(os.path.abspath(authorized_pic_root))),
        "authorized_project_home_root": str(
            Path(os.path.abspath(authorized_project_home_root))
        ),
    }
    return result


def consume_published_pressure_selection_receipt(
    receipt_path: str | Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    """Accept pressure selection only through the canonical receipt and exact seal."""
    return verify_published_pressure_selection_receipt(
        receipt_path, authorized_pic_root=authorized_pic_root
    )


def _same_capture(first: Mapping[str, object], second: Mapping[str, object]) -> bool:
    return (
        first["controller_state"] == second["controller_state"]
        and first["active_policy_payload"] == second["active_policy_payload"]
        and first["active_promotion_payload"] == second["active_promotion_payload"]
    )


def publish_pressure_selection(
    receipt: object,
    *,
    candidate_publication_authorization: object,
    receipt_path: str | Path | None = None,
    controller_operator_id: str,
    expected_git_commit: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    now: datetime | None = None,
) -> dict[str, object]:
    """Publish the exact reviewed p0=1.0 choice without granting launch authority."""
    publisher_source_authentication = _runtime_publisher_source_authentication(
        expected_git_commit
    )
    pic_root, publication_root = pilot_publisher._publication_root(authorized_pic_root)
    acceptance_root = pilot_publisher._publication_acceptance_root(pic_root)
    archive_root = pilot_publisher._canonical_existing_directory(
        pic_root / PRESSURE_GATE_ATTESTATION_ROOT_NAME,
        "pressure-gate attestation archive root",
    )
    transaction_anchor = pilot_publisher._canonical_existing_directory(
        control_plane_common.stable_serialization_anchor(pic_root),
        "stable publication transaction anchor",
    )
    target = pilot_publisher._direct_publication_target(
        publication_root / CANONICAL_RECEIPT_NAME
        if receipt_path is None
        else receipt_path,
        publication_root,
        "pressure-selection receipt",
    )
    _require(target.name == CANONICAL_RECEIPT_NAME, "pressure-selection receipt name drifted")
    operator_id = _operator_id(controller_operator_id, "controller-state operator ID")
    normalized, reanalysis, reviewer = _validate_reviewed_selection(
        receipt, authorized_pic_root=pic_root
    )
    receipt_payload = pressure_selection.canonical_json_bytes(normalized)
    candidate_authorization = _authorize_candidate_publication(
        candidate_publication_authorization,
        receipt_payload=receipt_payload,
        receipt=normalized,
        reanalysis_verification=reanalysis,
        reviewer_verification=reviewer,
        publisher_source_authentication=publisher_source_authentication,
        authorized_pic_root=pic_root,
    )
    publication_descriptor = pilot_publisher._open_absolute_directory(publication_root)
    acceptance_descriptor = pilot_publisher._open_absolute_directory(acceptance_root)
    archive_descriptor = pilot_publisher._open_absolute_directory(archive_root)
    transaction_descriptor = pilot_publisher._open_absolute_directory(transaction_anchor)
    staging_name = f".{target.name}.staging-{uuid.uuid4()}"
    staging_identity: tuple[int, int] | None = None
    receipt_identity: tuple[int, int] | None = None
    guard_identity: tuple[int, int] | None = None
    controller_binding: dict[str, str] | None = None
    guard_armed = False
    receipt_published = False
    seal_committed = False
    result: dict[str, object] | None = None
    try:
        pilot_publisher._lock_publication_transaction(
            transaction_descriptor, acceptance_descriptor
        )
        pilot_publisher._require_same_directory(
            transaction_anchor,
            transaction_descriptor,
            "stable publication transaction anchor",
        )
        pilot_publisher._require_same_directory(
            publication_root, publication_descriptor, "authorized PIC publication root"
        )
        pilot_publisher._require_same_directory(
            acceptance_root,
            acceptance_descriptor,
            "authorized PIC publication acceptance root",
        )
        pilot_publisher._require_same_directory(
            archive_root, archive_descriptor, "pressure-gate attestation archive root"
        )
        archive_metadata = os.fstat(archive_descriptor)
        _require(
            stat.S_IMODE(archive_metadata.st_mode) == 0o700
            and archive_metadata.st_uid == os.geteuid(),
            "pressure-gate attestation archive root is not private",
        )
        pilot_publisher._require_same_account_isolated_parent(publication_descriptor)
        pilot_publisher._require_absent_at(
            publication_descriptor, target.name, "pressure-selection receipt"
        )
        pilot_publisher._require_publication_guard_absent_at(
            publication_descriptor, target.name, "pressure-selection receipt"
        )
        pilot_publisher._require_absent_at(
            acceptance_descriptor,
            _selection_success_seal_name(target.name),
            "pressure-selection durable success seal",
        )
        captured_utc = _canonical_utc(now)
        capture = _capture_live_controller_state(
            authorized_pic_root=pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        sealed_utc = _canonical_utc(now)
        attestation = build_controller_state_attestation(
            capture,
            reanalysis_verification=reanalysis,
            publisher_source_authentication=publisher_source_authentication,
            candidate_publication_authorization=candidate_authorization,
            operator_id=operator_id,
            captured_utc=captured_utc,
            sealed_utc=sealed_utc,
            authorized_pic_root=pic_root,
        )
        controller_binding = _write_controller_state_attestation_at(
            archive_root,
            archive_descriptor,
            attestation,
            capture["active_policy_payload"],
            capture["active_promotion_payload"],
        )
        consumed_controller = consume_sealed_controller_state_attestation(
            controller_binding, authorized_pic_root=pic_root
        )
        _require_capture_matches_attestation(capture, consumed_controller)
        _recovery_guard, guard_identity = _arm_selection_recovery_guard_at(
            publication_descriptor,
            target,
            receipt_sha256=_sha256(receipt_payload),
            controller_state_attestation=controller_binding,
            expected_git_commit=expected_git_commit,
        )
        guard_armed = True
        pilot_publisher._write_exclusive_at(
            publication_descriptor, staging_name, receipt_payload
        )
        staging_identity = pilot_publisher._file_identity_at(
            publication_descriptor, staging_name, "staged pressure-selection receipt"
        )
        pilot_publisher._fsync_descriptor(publication_descriptor)
        pilot_publisher._rename_no_replace_at(
            publication_descriptor, staging_name, target.name
        )
        receipt_identity = staging_identity
        staging_identity = None
        receipt_published = True
        pilot_publisher._fsync_descriptor(publication_descriptor)
        _verify_published_pressure_selection(
            target,
            authorized_pic_root=pic_root,
            allow_publication_guard=True,
            require_publication_seal=False,
            expected_controller_state_attestation=controller_binding,
        )
        recaptured = _capture_live_controller_state(
            authorized_pic_root=pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        _require(
            _same_capture(capture, recaptured),
            "live controller state changed before pressure-selection acceptance",
        )
        _publish_selection_success_seal_at(
            acceptance_descriptor,
            publication_descriptor,
            target.name,
            receipt_payload,
            receipt_identity,
            controller_binding,
            sealed_utc=sealed_utc,
        )
        seal_committed = True
        final_capture = _capture_live_controller_state(
            authorized_pic_root=pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        _require(
            _same_capture(capture, final_capture),
            "live controller state changed before pressure-selection guard removal",
        )
        pilot_publisher._disarm_publication_guard_at(
            publication_descriptor, target.name, guard_identity
        )
        guard_armed = False
        guard_identity = None
        result = verify_published_pressure_selection_receipt(
            target, authorized_pic_root=pic_root
        )
        return result
    except BaseException as publication_error:
        if not guard_armed:
            try:
                _recovery_guard, guard_identity = _read_selection_recovery_guard_at(
                    publication_descriptor,
                    target,
                    receipt_sha256=_sha256(receipt_payload),
                    controller_state_attestation=controller_binding,
                    expected_git_commit=expected_git_commit,
                )
            except BaseException:
                pass
            else:
                guard_armed = True
        if seal_committed:
            try:
                result = verify_published_pressure_selection_receipt(
                    target, authorized_pic_root=pic_root
                )
            except BaseException:
                pass
            else:
                guard_armed = False
                return result
        if guard_armed or receipt_published:
            try:
                recovery_guard, guard_identity = _ensure_selection_recovery_guard_at(
                    publication_descriptor,
                    target,
                    receipt_sha256=_sha256(receipt_payload),
                    controller_state_attestation=controller_binding,
                    expected_git_commit=expected_git_commit,
                )
                guard_armed = True
            except BaseException as error:
                raise PressureSelectionPublicationError(
                    "pressure-selection receipt requires reviewed reconciliation "
                    "and its fail-closed guard could not be assured"
                ) from error
            raise PressureSelectionPublicationError(
                "pressure-selection receipt retained under fail-closed guard; "
                f"reviewed reconciliation required: {publication_error}",
                recovery=recovery_guard,
            ) from publication_error
        if isinstance(publication_error, PressureSelectionPublicationError):
            raise
        raise PressureSelectionPublicationError(
            f"pressure-selection publication failed closed: {publication_error}"
        ) from publication_error
    finally:
        try:
            _unlink_exact_file_at(
                publication_descriptor, staging_name, staging_identity
            )
        finally:
            close_error = pilot_publisher._close_descriptors(
                (
                    archive_descriptor,
                    acceptance_descriptor,
                    publication_descriptor,
                    transaction_descriptor,
                )
            )
            if close_error is not None and not seal_committed:
                raise close_error


def reconcile_pressure_selection_publication(
    receipt_path: str | Path,
    *,
    controller_state_attestation: object,
    expected_git_commit: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    now: datetime | None = None,
) -> dict[str, object]:
    """Complete or confirm one exact guarded pressure-selection publication."""
    runtime_source_authentication = _runtime_publisher_source_authentication(
        expected_git_commit
    )
    pic_root, publication_root = pilot_publisher._publication_root(authorized_pic_root)
    acceptance_root = pilot_publisher._publication_acceptance_root(pic_root)
    transaction_anchor = pilot_publisher._canonical_existing_directory(
        control_plane_common.stable_serialization_anchor(pic_root),
        "stable publication transaction anchor",
    )
    target = pilot_publisher._direct_publication_target(
        receipt_path, publication_root, "pressure-selection receipt"
    )
    _require(target.name == CANONICAL_RECEIPT_NAME, "pressure-selection receipt name drifted")
    controller_binding = _binding(
        controller_state_attestation, "controller-state attestation"
    )
    publication_descriptor = pilot_publisher._open_absolute_directory(publication_root)
    acceptance_descriptor = pilot_publisher._open_absolute_directory(acceptance_root)
    transaction_descriptor = pilot_publisher._open_absolute_directory(transaction_anchor)
    try:
        pilot_publisher._lock_publication_transaction(
            transaction_descriptor, acceptance_descriptor
        )
        try:
            published = verify_published_pressure_selection_receipt(
                target, authorized_pic_root=pic_root
            )
        except (OSError, ValueError):
            pass
        else:
            _require(
                published["controller_state_attestation"]["binding"]
                == controller_binding,
                "published pressure-selection receipt binds another controller state",
            )
            _require(
                published["controller_state_attestation"]["attestation"][
                    "publisher_source_authentication"
                ]
                == runtime_source_authentication,
                "published pressure-selection controller differs from the runtime publisher source",
            )
            return published
        recovery_guard, guard_identity = _read_selection_recovery_guard_at(
            publication_descriptor,
            target,
            controller_state_attestation=controller_binding,
            expected_git_commit=expected_git_commit,
        )
        receipt_sha256 = str(recovery_guard["receipt"]["sha256"])
        try:
            os.stat(target.name, dir_fd=publication_descriptor, follow_symlinks=False)
        except FileNotFoundError:
            pilot_publisher._require_absent_at(
                acceptance_descriptor,
                _selection_success_seal_name(target.name),
                "pressure-selection durable success seal",
            )
            staging_prefix = f".{target.name}.staging-"
            _require(
                not any(
                    name.startswith(staging_prefix)
                    for name in os.listdir(publication_descriptor)
                ),
                "guard-only pressure-selection recovery has retained staging residue",
            )
            controller = consume_sealed_controller_state_attestation(
                controller_binding, authorized_pic_root=pic_root
            )
            _require(
                controller["attestation"]["publisher_source_authentication"]
                == runtime_source_authentication,
                "guarded pressure-selection controller differs from the runtime publisher source",
            )
            capture = _capture_live_controller_state(
                authorized_pic_root=pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )
            _require_capture_matches_attestation(capture, controller)
            recaptured = _capture_live_controller_state(
                authorized_pic_root=pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )
            _require_capture_matches_attestation(recaptured, controller)
            pilot_publisher._disarm_publication_guard_at(
                publication_descriptor, target.name, guard_identity
            )
            return {
                "status": "no_publication_exposed_guard_cleared",
                "recovery_guard": recovery_guard,
                "controller_state_attestation": controller,
                "qualification_effect": QUALIFICATION_EFFECT,
            }
        except OSError as error:
            raise PressureSelectionPublicationError(
                "guarded pressure-selection receipt cannot be inspected"
            ) from error
        guarded = _verify_published_pressure_selection(
            target,
            authorized_pic_root=pic_root,
            allow_publication_guard=True,
            require_publication_seal=False,
            expected_controller_state_attestation=controller_binding,
        )
        controller = guarded["controller_state_attestation"]
        _require(
            controller["attestation"]["publisher_source_authentication"]
            == runtime_source_authentication,
            "guarded pressure-selection controller differs from the runtime publisher source",
        )
        capture = _capture_live_controller_state(
            authorized_pic_root=pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        _require_capture_matches_attestation(capture, controller)
        receipt_payload, receipt_identity = _read_exact_mode_regular_at(
            publication_descriptor,
            target.name,
            mode=0o444,
            label="pressure-selection receipt",
        )
        _require(
            _sha256(receipt_payload) == receipt_sha256,
            "guarded pressure-selection receipt differs from the recovery guard",
        )
        before_seal = _capture_live_controller_state(
            authorized_pic_root=pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        _require_capture_matches_attestation(before_seal, controller)
        try:
            seal, _ = _read_selection_success_seal_at(
                acceptance_descriptor,
                publication_descriptor,
                target.name,
                receipt_payload,
                receipt_identity,
            )
            _require(
                seal["controller_state_attestation"] == controller_binding,
                "existing pressure-selection success seal binds another controller state",
            )
        except FileNotFoundError:
            _publish_selection_success_seal_at(
                acceptance_descriptor,
                publication_descriptor,
                target.name,
                receipt_payload,
                receipt_identity,
                controller_binding,
                sealed_utc=_canonical_utc(now),
            )
        before_guard_removal = _capture_live_controller_state(
            authorized_pic_root=pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        _require_capture_matches_attestation(before_guard_removal, controller)
        pilot_publisher._disarm_publication_guard_at(
            publication_descriptor, target.name, guard_identity
        )
        return verify_published_pressure_selection_receipt(
            target, authorized_pic_root=pic_root
        )
    except PressureSelectionPublicationError:
        raise
    except BaseException as error:
        raise PressureSelectionPublicationError(
            f"pressure-selection reconciliation failed closed: {error}"
        ) from error
    finally:
        pilot_publisher._close_descriptors(
            (acceptance_descriptor, publication_descriptor, transaction_descriptor)
        )


def _load_json(path: Path, label: str) -> dict[str, object]:
    try:
        payload = pilot_publisher._read_stable_readonly_regular(
            path, label, max_bytes=MAX_JSON_BYTES
        )
    except (OSError, TypeError, ValueError) as error:
        raise PressureSelectionPublicationError(f"{label} is unavailable") from error
    return _decode_canonical_object(payload, label)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare-reanalysis")
    prepare_parser.add_argument("--reanalysis-operator-id", required=True)
    human_parser = subparsers.add_parser("seal-human-selection")
    human_parser.add_argument("--human-decision", required=True, type=Path)
    publish_parser = subparsers.add_parser("publish-selection")
    publish_parser.add_argument("--receipt", required=True, type=Path)
    publish_parser.add_argument("--candidate-authorization", required=True, type=Path)
    publish_parser.add_argument("--controller-operator-id", required=True)
    verify_parser = subparsers.add_parser("verify-published")
    verify_parser.add_argument(
        "--receipt",
        type=Path,
        default=AUTHORIZED_PIC_ROOT / "publication" / CANONICAL_RECEIPT_NAME,
    )
    reconcile_parser = subparsers.add_parser("reconcile-published")
    reconcile_parser.add_argument(
        "--receipt",
        type=Path,
        default=AUTHORIZED_PIC_ROOT / "publication" / CANONICAL_RECEIPT_NAME,
    )
    reconcile_parser.add_argument("--controller-state-attestation", required=True, type=Path)
    reconcile_parser.add_argument("--expected-controller-state-sha256", required=True)
    for command_parser in (
        prepare_parser,
        human_parser,
        publish_parser,
        reconcile_parser,
    ):
        command_parser.add_argument("--expected-git-commit", required=True)
    for command_parser in (
        prepare_parser,
        human_parser,
        publish_parser,
        verify_parser,
        reconcile_parser,
    ):
        command_parser.add_argument(
            "--authorized-pic-root", type=Path, default=AUTHORIZED_PIC_ROOT
        )
        command_parser.add_argument(
            "--authorized-project-home-root",
            type=Path,
            default=AUTHORIZED_PROJECT_HOME_ROOT,
        )
    args = parser.parse_args(argv)
    try:
        if args.command == "prepare-reanalysis":
            result = prepare_pressure_reanalysis(
                reanalysis_operator_id=args.reanalysis_operator_id,
                expected_git_commit=args.expected_git_commit,
                authorized_pic_root=args.authorized_pic_root,
                authorized_project_home_root=args.authorized_project_home_root,
            )
        elif args.command == "seal-human-selection":
            result = seal_human_pressure_selection(
                human_decision_path=args.human_decision,
                expected_git_commit=args.expected_git_commit,
                authorized_pic_root=args.authorized_pic_root,
                authorized_project_home_root=args.authorized_project_home_root,
            )
        elif args.command == "publish-selection":
            result = publish_pressure_selection(
                _load_json(args.receipt, "candidate pressure-selection receipt"),
                candidate_publication_authorization=_load_json(
                    args.candidate_authorization,
                    "candidate publication authorization",
                ),
                controller_operator_id=args.controller_operator_id,
                expected_git_commit=args.expected_git_commit,
                authorized_pic_root=args.authorized_pic_root,
                authorized_project_home_root=args.authorized_project_home_root,
            )
        elif args.command == "verify-published":
            result = verify_published_pressure_selection_live_state(
                args.receipt,
                authorized_pic_root=args.authorized_pic_root,
                authorized_project_home_root=args.authorized_project_home_root,
            )
        else:
            result = reconcile_pressure_selection_publication(
                args.receipt,
                controller_state_attestation={
                    "path": str(args.controller_state_attestation),
                    "sha256": args.expected_controller_state_sha256,
                },
                expected_git_commit=args.expected_git_commit,
                authorized_pic_root=args.authorized_pic_root,
                authorized_project_home_root=args.authorized_project_home_root,
            )
    except PressureSelectionPublicationError as error:
        failure: dict[str, object] = {
            "schema_version": 1,
            "record_type": "q011_section54_pressure_selection_command_failure",
            "status": "failed_closed",
            "command": args.command,
            "message": str(error),
            "recovery": error.recovery,
        }
        sys.stderr.write(_canonical_json_bytes(failure).decode("utf-8"))
        return 2
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "CANONICAL_RECEIPT_NAME",
    "CANDIDATE_AUTHORIZATION_RECORD_TYPE",
    "CONTROLLER_STATE_RECORD_TYPE",
    "HUMAN_DECISION_RECORD_TYPE",
    "HUMAN_DECISION_STATEMENT",
    "PressureSelectionPublicationError",
    "QUALIFICATION_EFFECT",
    "REVIEWED_RATIONALE",
    "REVIEWED_REVIEWER_ID",
    "REVIEWED_SELECTED_CASE",
    "SUCCESS_SEAL_RECORD_TYPE",
    "STAGE4_PREPARATION_RECORD_TYPE",
    "build_controller_state_attestation",
    "build_pressure_selection_receipt",
    "consume_published_pressure_selection_receipt",
    "consume_sealed_controller_state_attestation",
    "prepare_pressure_reanalysis",
    "publish_pressure_selection",
    "reconcile_pressure_selection_publication",
    "seal_human_pressure_selection",
    "verify_published_pressure_selection_live_state",
    "verify_published_pressure_selection_receipt",
]
