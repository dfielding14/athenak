#!/usr/bin/env python3
"""Immutable derived-artifact publication helpers for Q-011 Section 5.4."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
from typing import Any, Mapping
import uuid


INVENTORY_NAME = "artifact_inventory.json"
MANIFEST_NAME = "derived_manifest.json"
_RESERVED_NAMES = {INVENTORY_NAME, MANIFEST_NAME}
_SHA256_LENGTH = 64


class DerivedArtifactError(ValueError):
    """Raised when derived-artifact publication or verification fails closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise DerivedArtifactError(message)


def _reject_constant(value: str) -> None:
    raise DerivedArtifactError(f"JSON constant is forbidden: {value}")


def _decode_json(payload: bytes, label: str) -> Any:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result = {}
        for key, value in pairs:
            _require(key not in result, f"{label}: duplicate JSON key {key!r}")
            result[key] = value
        return result

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise DerivedArtifactError(f"{label}: JSON is not UTF-8") from error
    try:
        return json.loads(
            text,
            object_pairs_hook=reject_duplicates,
            parse_constant=_reject_constant,
        )
    except json.JSONDecodeError as error:
        raise DerivedArtifactError(f"{label}: malformed JSON") from error


def canonical_json_bytes(value: object) -> bytes:
    """Serialize one JSON value deterministically and reject non-finite numbers."""
    try:
        return (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise DerivedArtifactError("value is not canonical JSON") from error


def sha256_bytes(payload: bytes) -> str:
    """Return a lowercase SHA-256 digest."""
    return hashlib.sha256(payload).hexdigest()


def _sha256(value: object, label: str) -> str:
    _require(
        type(value) is str
        and len(value) == _SHA256_LENGTH
        and all(character in "0123456789abcdef" for character in value),
        f"{label}: malformed SHA-256",
    )
    return value


def _relative_path(value: object, label: str) -> str:
    _require(type(value) is str and bool(value), f"{label}: expected path text")
    path = PurePosixPath(value)
    _require(
        not path.is_absolute()
        and path.as_posix() == value
        and value != "."
        and all(part not in {"", ".", ".."} for part in path.parts),
        f"{label}: unsafe relative path",
    )
    return value


def validate_reviewer_disposition(value: object) -> dict[str, object]:
    """Validate the explicit human-review boundary for one derived bundle."""
    _require(type(value) is dict, "reviewer disposition must be an object")
    _require(
        set(value)
        == {
            "status",
            "reviewer",
            "reviewed_at_utc",
            "notes",
        },
        "reviewer disposition keys drifted",
    )
    status = value["status"]
    reviewer = value["reviewer"]
    reviewed_at = value["reviewed_at_utc"]
    notes = value["notes"]
    _require(
        status in {"pending_external_review", "accepted", "rejected"},
        "reviewer disposition status is unsupported",
    )
    _require(type(notes) is str, "reviewer disposition notes must be text")
    if status == "pending_external_review":
        _require(reviewer is None, "pending review must not name an unassigned reviewer")
        _require(reviewed_at is None, "pending review must not carry a review timestamp")
    else:
        _require(
            type(reviewer) is str and bool(reviewer.strip()),
            "terminal reviewer disposition requires a named reviewer",
        )
        _require(
            type(reviewed_at) is str
            and reviewed_at.endswith("Z")
            and "T" in reviewed_at,
            "terminal reviewer disposition requires a UTC timestamp",
        )
    return {
        "status": status,
        "reviewer": reviewer,
        "reviewed_at_utc": reviewed_at,
        "notes": notes,
    }


def build_derived_manifest(
    *,
    campaign_id: str,
    raw_artifact_inventory_sha256: str,
    analyzer_bindings: Mapping[str, str],
    reviewer_disposition: object,
) -> dict[str, object]:
    """Build the minimal provenance manifest for one derived-artifact bundle."""
    _require(type(campaign_id) is str and bool(campaign_id), "campaign id is required")
    _sha256(raw_artifact_inventory_sha256, "raw artifact inventory")
    _require(
        type(analyzer_bindings) is dict and bool(analyzer_bindings),
        "analyzer bindings must be a nonempty object",
    )
    parsed_bindings = {}
    for path, digest in sorted(analyzer_bindings.items()):
        parsed_bindings[_relative_path(path, "analyzer binding path")] = _sha256(
            digest, f"analyzer binding {path}"
        )
    return {
        "record_type": "q011_section54_derived_artifact_bundle",
        "schema_version": 1,
        "campaign_id": campaign_id,
        "evidence_boundary": "derived_records_do_not_replace_retained_raw_artifacts",
        "raw_artifact_inventory_sha256": raw_artifact_inventory_sha256,
        "analyzer_bindings": parsed_bindings,
        "reviewer_disposition": validate_reviewer_disposition(reviewer_disposition),
    }


def validate_derived_manifest(value: object) -> dict[str, object]:
    """Validate the exact provenance schema for one derived-artifact bundle."""
    _require(type(value) is dict, "derived manifest must be an object")
    _require(
        set(value)
        == {
            "record_type",
            "schema_version",
            "campaign_id",
            "evidence_boundary",
            "raw_artifact_inventory_sha256",
            "analyzer_bindings",
            "reviewer_disposition",
        },
        "derived manifest keys drifted",
    )
    _require(
        value["record_type"] == "q011_section54_derived_artifact_bundle"
        and type(value["schema_version"]) is int
        and value["schema_version"] == 1,
        "derived manifest identity drifted",
    )
    _require(
        type(value["campaign_id"]) is str and bool(value["campaign_id"]),
        "derived manifest campaign id drifted",
    )
    _require(
        value["evidence_boundary"]
        == "derived_records_do_not_replace_retained_raw_artifacts",
        "derived manifest evidence boundary drifted",
    )
    _sha256(value["raw_artifact_inventory_sha256"], "raw artifact inventory")
    bindings = value["analyzer_bindings"]
    _require(type(bindings) is dict and bool(bindings), "analyzer bindings drifted")
    for path, digest in bindings.items():
        _relative_path(path, "analyzer binding path")
        _sha256(digest, f"analyzer binding {path}")
    validate_reviewer_disposition(value["reviewer_disposition"])
    return value


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _write_exclusive(path: Path, payload: bytes) -> None:
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o600,
    )
    try:
        with os.fdopen(descriptor, "wb", closefd=False) as stream:
            stream.write(payload)
            stream.flush()
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _inventory_record(relative: str, payload: bytes) -> dict[str, object]:
    return {
        "path": relative,
        "size": len(payload),
        "sha256": sha256_bytes(payload),
    }


def _freeze_directories(root: Path) -> None:
    directories = [Path(directory) for directory, _, _ in os.walk(root)]
    for directory in reversed(directories):
        os.chmod(directory, 0o555)
        _fsync_directory(directory)


def _cleanup_staging(path: Path) -> None:
    if not path.exists():
        return
    for directory, names, _ in os.walk(path):
        os.chmod(directory, 0o755)
        for name in names:
            os.chmod(Path(directory) / name, 0o755)
    shutil.rmtree(path)


def publish_derived_bundle(
    output_path: str | Path,
    *,
    manifest: Mapping[str, object],
    artifacts: Mapping[str, bytes],
) -> dict[str, str]:
    """Publish one new immutable bundle by exclusive same-parent rename."""
    target = Path(os.path.abspath(output_path))
    parent = target.parent
    _require(parent.is_dir(), "derived-artifact parent directory is unavailable")
    _require(not target.exists(), "derived-artifact output already exists")
    validate_derived_manifest(manifest)
    canonical_manifest = canonical_json_bytes(dict(manifest))
    payloads: dict[str, bytes] = {}
    for relative, payload in artifacts.items():
        path = _relative_path(relative, "derived-artifact member")
        _require(path not in _RESERVED_NAMES, f"derived-artifact member is reserved: {path}")
        _require(type(payload) is bytes, f"derived-artifact member must be bytes: {path}")
        _require(path not in payloads, f"duplicate derived-artifact member: {path}")
        payloads[path] = payload
    payloads[MANIFEST_NAME] = canonical_manifest
    inventory = {
        "record_type": "q011_section54_derived_artifact_inventory",
        "schema_version": 1,
        "members": [
            _inventory_record(relative, payload)
            for relative, payload in sorted(payloads.items())
        ],
    }
    inventory_payload = canonical_json_bytes(inventory)
    staging = parent / f".{target.name}.staging-{uuid.uuid4()}"
    staging.mkdir(mode=0o700)
    try:
        for relative, payload in sorted(payloads.items()):
            path = staging / relative
            path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
            _write_exclusive(path, payload)
        _write_exclusive(staging / INVENTORY_NAME, inventory_payload)
        _freeze_directories(staging)
        _require(not target.exists(), "derived-artifact output appeared during publication")
        os.rename(staging, target)
        _fsync_directory(parent)
    except BaseException:
        _cleanup_staging(staging)
        raise
    return {
        "path": str(target),
        "manifest_sha256": sha256_bytes(canonical_manifest),
        "inventory_sha256": sha256_bytes(inventory_payload),
    }


def _read_immutable_file(path: Path, label: str) -> bytes:
    _require(not path.is_symlink(), f"{label}: symlink is forbidden")
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode) and not before.st_mode & 0o222,
            f"{label}: expected read-only regular file",
        )
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            payload = stream.read()
        after = os.fstat(descriptor)
        stable = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns")
        _require(
            all(getattr(before, name) == getattr(after, name) for name in stable)
            and len(payload) == after.st_size,
            f"{label}: file changed while reading",
        )
        return payload
    finally:
        os.close(descriptor)


def verify_published_derived_bundle(
    output_path: str | Path, *, expected_inventory_sha256: str | None = None
) -> dict[str, object]:
    """Recompute one immutable bundle inventory and return its manifest."""
    root = Path(os.path.abspath(output_path))
    _require(root.is_dir() and not root.is_symlink(), "derived-artifact bundle is unavailable")
    _require(not root.stat().st_mode & 0o222, "derived-artifact root is writable")
    inventory_payload = _read_immutable_file(root / INVENTORY_NAME, INVENTORY_NAME)
    if expected_inventory_sha256 is not None:
        _require(
            sha256_bytes(inventory_payload)
            == _sha256(expected_inventory_sha256, "expected inventory"),
            "derived-artifact inventory checksum drifted",
        )
    inventory = _decode_json(inventory_payload, INVENTORY_NAME)
    _require(type(inventory) is dict, "derived-artifact inventory must be an object")
    _require(
        set(inventory) == {"record_type", "schema_version", "members"}
        and inventory["record_type"] == "q011_section54_derived_artifact_inventory"
        and type(inventory["schema_version"]) is int
        and inventory["schema_version"] == 1
        and type(inventory["members"]) is list,
        "derived-artifact inventory schema drifted",
    )
    declared = {INVENTORY_NAME}
    for item in inventory["members"]:
        _require(type(item) is dict and set(item) == {"path", "size", "sha256"}, "inventory member schema drifted")
        relative = _relative_path(item["path"], "inventory member path")
        _require(relative not in declared, f"inventory repeats member: {relative}")
        declared.add(relative)
        _require(type(item["size"]) is int and item["size"] >= 0, "inventory member size drifted")
        payload = _read_immutable_file(root / relative, relative)
        _require(
            len(payload) == item["size"]
            and sha256_bytes(payload) == _sha256(item["sha256"], relative),
            f"inventory member checksum drifted: {relative}",
        )
    actual = set()
    for directory, names, filenames in os.walk(root, followlinks=False):
        base = Path(directory)
        _require(not base.is_symlink(), "derived-artifact directory symlink is forbidden")
        _require(not base.stat().st_mode & 0o222, "derived-artifact directory is writable")
        for name in names:
            _require(not (base / name).is_symlink(), "derived-artifact directory symlink is forbidden")
        for name in filenames:
            actual.add((base / name).relative_to(root).as_posix())
    _require(actual == declared, "derived-artifact tree closure drifted")
    manifest = _decode_json(
        _read_immutable_file(root / MANIFEST_NAME, MANIFEST_NAME), MANIFEST_NAME
    )
    return validate_derived_manifest(manifest)


__all__ = [
    "DerivedArtifactError",
    "INVENTORY_NAME",
    "MANIFEST_NAME",
    "build_derived_manifest",
    "canonical_json_bytes",
    "publish_derived_bundle",
    "sha256_bytes",
    "validate_derived_manifest",
    "validate_reviewer_disposition",
    "verify_published_derived_bundle",
]
