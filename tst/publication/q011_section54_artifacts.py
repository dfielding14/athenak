#!/usr/bin/env python3
"""Immutable derived-artifact publication helpers for Q-011 Section 5.4."""

from __future__ import annotations

import ctypes
import errno
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import stat
from typing import Any, Mapping
import uuid

if __package__:
    from . import immutable_orion_tree
else:
    import immutable_orion_tree


INVENTORY_NAME = "artifact_inventory.json"
MANIFEST_NAME = "derived_manifest.json"
_RESERVED_NAMES = {INVENTORY_NAME, MANIFEST_NAME}
_SHA256_LENGTH = 64
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_RENAME_NOREPLACE = 1
_STAGING_ROOT_NAME = "publishable"


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


def _write_exclusive_at(root_fd: int, relative: str, payload: bytes) -> None:
    path = PurePosixPath(relative)
    descriptor = os.dup(root_fd)
    try:
        for part in path.parts[:-1]:
            try:
                os.mkdir(part, mode=0o700, dir_fd=descriptor)
            except FileExistsError:
                pass
            else:
                os.fsync(descriptor)
            child = os.open(part, _DIRECTORY_FLAGS, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        output = os.open(
            path.parts[-1],
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
            0o600,
            dir_fd=descriptor,
        )
    except OSError as error:
        os.close(descriptor)
        raise DerivedArtifactError(f"cannot create derived-artifact member: {relative}") from error
    try:
        offset = 0
        while offset < len(payload):
            written = os.write(output, payload[offset:])
            _require(written > 0, f"short write for derived-artifact member: {relative}")
            offset += written
        os.fsync(output)
        os.fchmod(output, 0o444)
        os.fsync(output)
        os.fsync(descriptor)
    finally:
        os.close(output)
        os.close(descriptor)


def _require_same_directory_at(parent_fd: int, name: str, descriptor: int, label: str) -> None:
    _require("/" not in name, f"{label}: nested descriptor-relative name")
    try:
        actual = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except OSError as error:
        raise DerivedArtifactError(f"{label}: directory is unavailable") from error
    expected = os.fstat(descriptor)
    _require(
        stat.S_ISDIR(actual.st_mode)
        and (actual.st_dev, actual.st_ino) == (expected.st_dev, expected.st_ino),
        f"{label}: directory binding changed",
    )


def _require_same_directory(path: Path, descriptor: int, label: str) -> None:
    try:
        actual = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise DerivedArtifactError(f"{label}: directory is unavailable") from error
    expected = os.fstat(descriptor)
    _require(
        stat.S_ISDIR(actual.st_mode)
        and (actual.st_dev, actual.st_ino) == (expected.st_dev, expected.st_ino),
        f"{label}: directory binding changed",
    )


def _require_absent_at(parent_fd: int, name: str, label: str) -> None:
    try:
        os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        return
    except OSError as error:
        raise DerivedArtifactError(f"{label}: cannot inspect path") from error
    raise DerivedArtifactError(f"{label}: already exists")


def _rename_no_replace_at(
    source_parent_fd: int,
    source_name: str,
    destination_parent_fd: int,
    destination_name: str,
) -> None:
    _require(
        "/" not in source_name and "/" not in destination_name,
        "derived-artifact rename received a nested name",
    )
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    error_number = errno.ENOSYS
    if renameat2 is not None:
        renameat2.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renameat2.restype = ctypes.c_int
        if renameat2(
            source_parent_fd,
            os.fsencode(source_name),
            destination_parent_fd,
            os.fsencode(destination_name),
            _RENAME_NOREPLACE,
        ) == 0:
            return
        error_number = ctypes.get_errno()
    if error_number == errno.EEXIST:
        raise DerivedArtifactError(f"derived-artifact output already exists: {destination_name}")
    unsupported = {
        errno.EINVAL,
        errno.ENOSYS,
        getattr(errno, "ENOTSUP", errno.EINVAL),
        getattr(errno, "EOPNOTSUPP", errno.EINVAL),
    }
    if error_number not in unsupported:
        raise OSError(error_number, os.strerror(error_number), destination_name)
    raise DerivedArtifactError("derived-artifact publication requires atomic no-replace rename")


def _remove_anchored_tree_at(parent_fd: int, name: str, descriptor: int, label: str) -> None:
    """Remove one tree recursively without reopening its ancestor pathname."""
    _require_same_directory_at(parent_fd, name, descriptor, label)

    def remove_members(directory_fd: int) -> None:
        status = os.fstat(directory_fd)
        os.fchmod(directory_fd, stat.S_IMODE(status.st_mode) | 0o700)
        for member_name in os.listdir(directory_fd):
            observed = os.stat(member_name, dir_fd=directory_fd, follow_symlinks=False)
            if stat.S_ISDIR(observed.st_mode):
                child_fd = os.open(member_name, _DIRECTORY_FLAGS, dir_fd=directory_fd)
                try:
                    opened = os.fstat(child_fd)
                    _require(
                        (observed.st_dev, observed.st_ino) == (opened.st_dev, opened.st_ino),
                        f"{label}: tree changed during anchored removal",
                    )
                    remove_members(child_fd)
                    current = os.stat(member_name, dir_fd=directory_fd, follow_symlinks=False)
                    _require(
                        (current.st_dev, current.st_ino) == (opened.st_dev, opened.st_ino),
                        f"{label}: tree changed during anchored removal",
                    )
                finally:
                    os.close(child_fd)
                os.rmdir(member_name, dir_fd=directory_fd)
            else:
                os.unlink(member_name, dir_fd=directory_fd)

    try:
        remove_members(descriptor)
        _require_same_directory_at(parent_fd, name, descriptor, label)
        os.rmdir(name, dir_fd=parent_fd)
    except DerivedArtifactError:
        raise
    except OSError as error:
        raise DerivedArtifactError(f"{label}: cannot remove tree") from error


def _cleanup_staging(parent_fd: int, name: str, descriptor: int) -> None:
    try:
        _remove_anchored_tree_at(parent_fd, name, descriptor, "derived-artifact staging tree")
    except (DerivedArtifactError, OSError):
        return


def _cleanup_private_container(parent_fd: int, name: str, descriptor: int) -> None:
    """Best-effort removal of one empty pinned private staging container."""
    try:
        _require_same_directory_at(
            parent_fd, name, descriptor, "derived-artifact private staging container"
        )
        os.rmdir(name, dir_fd=parent_fd)
        os.fsync(parent_fd)
    except (DerivedArtifactError, OSError):
        return


def _required_directories(paths: set[str]) -> set[str]:
    directories: set[str] = set()
    for relative in paths:
        parent = PurePosixPath(relative).parent
        while parent.as_posix() != ".":
            directories.add(parent.as_posix())
            parent = parent.parent
    return directories


def _verify_anchored_bundle(
    root_fd: int, *, expected_inventory_sha256: str | None = None
) -> dict[str, object]:
    snapshot = immutable_orion_tree._scan_anchored_tree(
        root_fd,
        hash_regular=True,
        capture_paths=frozenset({INVENTORY_NAME, MANIFEST_NAME}),
        error_type=DerivedArtifactError,
        label="Q-011 derived-artifact bundle",
    )
    immutable_orion_tree._require_read_only(
        snapshot,
        include_root=True,
        error_type=DerivedArtifactError,
        label="Q-011 derived-artifact bundle",
    )
    entries = snapshot.by_relative_path()
    inventory_entry = entries.get(INVENTORY_NAME)
    _require(
        inventory_entry is not None
        and inventory_entry.entry_type == "file"
        and inventory_entry.captured_payload is not None,
        "derived-artifact inventory is unavailable",
    )
    inventory_payload = inventory_entry.captured_payload
    if expected_inventory_sha256 is not None:
        _require(
            inventory_entry.sha256 == _sha256(expected_inventory_sha256, "expected inventory"),
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
        _require(
            type(item) is dict and set(item) == {"path", "size", "sha256"},
            "inventory member schema drifted",
        )
        relative = _relative_path(item["path"], "inventory member path")
        _require(relative not in declared, f"inventory repeats member: {relative}")
        declared.add(relative)
        _require(type(item["size"]) is int and item["size"] >= 0, "inventory member size drifted")
        entry = entries.get(relative)
        _require(
            entry is not None
            and entry.entry_type == "file"
            and entry.size == item["size"]
            and entry.sha256 == _sha256(item["sha256"], relative),
            f"inventory member checksum drifted: {relative}",
        )
    measured_files = {
        relative for relative, entry in entries.items() if entry.entry_type == "file"
    }
    measured_directories = {
        relative for relative, entry in entries.items() if entry.entry_type == "directory"
    }
    _require(measured_files == declared, "derived-artifact tree closure drifted")
    _require(
        measured_directories == _required_directories(declared),
        "derived-artifact directory closure drifted",
    )
    manifest_entry = entries.get(MANIFEST_NAME)
    _require(
        manifest_entry is not None and manifest_entry.captured_payload is not None,
        "derived manifest is unavailable",
    )
    final = immutable_orion_tree._scan_anchored_tree(
        root_fd,
        hash_regular=False,
        error_type=DerivedArtifactError,
        label="Q-011 derived-artifact bundle",
    )
    immutable_orion_tree._require_same_snapshot(
        snapshot,
        final,
        error_type=DerivedArtifactError,
        label="Q-011 derived-artifact bundle",
        phase="derived-artifact stability pass",
    )
    return validate_derived_manifest(_decode_json(manifest_entry.captured_payload, MANIFEST_NAME))


def _verify_anchored_bundle_at(
    parent_fd: int,
    name: str,
    descriptor: int,
    *,
    expected_inventory_sha256: str | None = None,
) -> dict[str, object]:
    _require_same_directory_at(parent_fd, name, descriptor, "derived-artifact bundle")
    verified = _verify_anchored_bundle(
        descriptor, expected_inventory_sha256=expected_inventory_sha256
    )
    _require_same_directory_at(parent_fd, name, descriptor, "derived-artifact bundle")
    return verified


def _rollback_published_destination(
    parent_fd: int, destination_name: str, descriptor: int
) -> None:
    """Withdraw only the pinned invalid bundle, never a path replacement."""
    rollback_name = f".{destination_name}.rollback-{uuid.uuid4()}"

    def quarantine_moved_original() -> None:
        expected = os.fstat(descriptor)
        candidates = []
        for name in os.listdir(parent_fd):
            try:
                observed = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
            except FileNotFoundError:
                continue
            if (
                stat.S_ISDIR(observed.st_mode)
                and (observed.st_dev, observed.st_ino)
                == (expected.st_dev, expected.st_ino)
            ):
                candidates.append(name)
        if len(candidates) != 1:
            raise DerivedArtifactError(
                "cannot locate raced derived-artifact publication for quarantine"
            )
        _rename_no_replace_at(parent_fd, candidates[0], parent_fd, rollback_name)
        os.fsync(parent_fd)
        _require_same_directory_at(
            parent_fd, rollback_name, descriptor, "derived-artifact rollback tree"
        )

    raced_replacement = False
    try:
        _rename_no_replace_at(parent_fd, destination_name, parent_fd, rollback_name)
    except FileNotFoundError:
        raced_replacement = True
        os.fsync(parent_fd)
        quarantine_moved_original()
    else:
        os.fsync(parent_fd)
        try:
            _require_same_directory_at(
                parent_fd, rollback_name, descriptor, "derived-artifact rollback tree"
            )
        except DerivedArtifactError:
            raced_replacement = True
            _rename_no_replace_at(parent_fd, rollback_name, parent_fd, destination_name)
            os.fsync(parent_fd)
            quarantine_moved_original()
    if not raced_replacement:
        _require_absent_at(parent_fd, destination_name, "invalid derived-artifact output")
    rollback_fd = os.open(rollback_name, _DIRECTORY_FLAGS, dir_fd=parent_fd)
    try:
        _remove_anchored_tree_at(
            parent_fd, rollback_name, rollback_fd, "derived-artifact rollback tree"
        )
    finally:
        os.close(rollback_fd)
    os.fsync(parent_fd)
    if raced_replacement:
        raise DerivedArtifactError(
            "derived-artifact rollback rejected a substituted public destination"
        )


def _write_exclusive(path: Path, payload: bytes) -> None:
    """Retained compatibility helper for direct path-local unit use."""
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


def publish_derived_bundle(
    output_path: str | Path,
    *,
    manifest: Mapping[str, object],
    artifacts: Mapping[str, bytes],
) -> dict[str, str]:
    """Publish one new immutable bundle by exclusive same-parent rename."""
    target = Path(os.path.abspath(output_path))
    _require(bool(target.name), "derived-artifact output name is unavailable")
    try:
        parent = target.parent.resolve(strict=True)
    except OSError as error:
        raise DerivedArtifactError("derived-artifact parent directory is unavailable") from error
    _require(parent == target.parent, "derived-artifact parent must be canonical")
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
    private_container = parent / f".{target.name}.staging-{uuid.uuid4()}"
    staging = private_container / _STAGING_ROOT_NAME
    parent_fd = os.open(parent, _DIRECTORY_FLAGS)
    private_fd: int | None = None
    staging_fd: int | None = None
    renamed = False
    try:
        _require_same_directory(parent, parent_fd, "derived-artifact parent")
        _require_absent_at(parent_fd, target.name, "derived-artifact output")
        os.mkdir(private_container.name, mode=0o700, dir_fd=parent_fd)
        os.fsync(parent_fd)
        private_fd = os.open(private_container.name, _DIRECTORY_FLAGS, dir_fd=parent_fd)
        os.mkdir(staging.name, mode=0o700, dir_fd=private_fd)
        os.fsync(private_fd)
        staging_fd = os.open(staging.name, _DIRECTORY_FLAGS, dir_fd=private_fd)
        for relative, payload in sorted(payloads.items()):
            _write_exclusive_at(staging_fd, relative, payload)
        _write_exclusive_at(staging_fd, INVENTORY_NAME, inventory_payload)
        immutable_orion_tree._remove_write_bits_below_root(
            staging_fd,
            error_type=DerivedArtifactError,
            label="Q-011 derived-artifact bundle",
        )
        os.fchmod(staging_fd, os.fstat(staging_fd).st_mode & ~0o222)
        _verify_anchored_bundle_at(
            private_fd,
            staging.name,
            staging_fd,
            expected_inventory_sha256=sha256_bytes(inventory_payload),
        )
        _require_same_directory(parent, parent_fd, "derived-artifact parent")
        _require_absent_at(parent_fd, target.name, "derived-artifact output")
        # Orion rejects cross-parent rename of a read-only directory. Descendants
        # remain frozen; make the root owner-write-only and non-traversable for
        # the rename, then restore its exact frozen mode through the pinned fd.
        root_mode = stat.S_IMODE(os.fstat(staging_fd).st_mode)
        os.fchmod(staging_fd, stat.S_IWUSR)
        os.fsync(staging_fd)
        _rename_no_replace_at(private_fd, staging.name, parent_fd, target.name)
        renamed = True
        os.fchmod(staging_fd, root_mode)
        os.fsync(staging_fd)
        os.fsync(private_fd)
        os.fsync(parent_fd)
        _verify_anchored_bundle_at(
            parent_fd,
            target.name,
            staging_fd,
            expected_inventory_sha256=sha256_bytes(inventory_payload),
        )
    except BaseException:
        if staging_fd is not None:
            if renamed:
                try:
                    _rollback_published_destination(parent_fd, target.name, staging_fd)
                except BaseException as rollback_error:
                    raise DerivedArtifactError(
                        "cannot remove invalid derived-artifact output"
                    ) from rollback_error
            elif private_fd is not None:
                _cleanup_staging(private_fd, staging.name, staging_fd)
        raise
    finally:
        if staging_fd is not None:
            os.close(staging_fd)
        if private_fd is not None:
            _cleanup_private_container(parent_fd, private_container.name, private_fd)
            os.close(private_fd)
        os.close(parent_fd)
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
    _require(bool(root.name), "derived-artifact bundle is unavailable")
    parent_fd: int | None = None
    try:
        parent = root.parent.resolve(strict=True)
        _require(parent == root.parent, "derived-artifact parent must be canonical")
        parent_fd = os.open(parent, _DIRECTORY_FLAGS)
        _require_same_directory(parent, parent_fd, "derived-artifact parent")
        root_fd = os.open(root.name, _DIRECTORY_FLAGS, dir_fd=parent_fd)
    except DerivedArtifactError:
        if parent_fd is not None:
            os.close(parent_fd)
        raise
    except OSError as error:
        if parent_fd is not None:
            os.close(parent_fd)
        raise DerivedArtifactError("derived-artifact bundle is unavailable") from error
    try:
        return _verify_anchored_bundle_at(
            parent_fd,
            root.name,
            root_fd,
            expected_inventory_sha256=expected_inventory_sha256,
        )
    finally:
        os.close(root_fd)
        os.close(parent_fd)


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
