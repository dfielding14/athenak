#!/usr/bin/env python3
"""Freeze one source-local Q-011 Section 5.4 raw campaign attempt.

This publisher only retains an already-produced raw attempt.  It does not
submit work, mutate policy, or claim numerical qualification.
"""

from __future__ import annotations

import argparse
import copy
import ctypes
import errno
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Any, Mapping, Sequence
import uuid

if __package__:
    from . import analyze_q011_section54_campaign as campaign
    from . import immutable_orion_tree
else:
    import analyze_q011_section54_campaign as campaign
    import immutable_orion_tree


ATTEMPT_STATUSES = ("completed", "failed")
RESULT_RECORD_TYPE = "q011_section54_campaign_attempt_freeze_result"
FAILED_MANIFEST_RECORD_TYPE = "q011_section54_failed_campaign_attempt_manifest"
_COMPLETED_RECEIPT = dict(campaign._EXPECTED_FREEZE_RECEIPT)
_FAILED_RECEIPT = {
    "schema_version": 1,
    "artifact_role": "q011_section54_failed_campaign_attempt",
    "qualification_effect": "retained_failed_campaign_attempt_no_admission",
    "inventory_excludes": immutable_orion_tree.INVENTORY_NAME,
    "freeze_policy": "remove all owner, group and other write bits recursively",
}
_RESERVED_SOURCE_MEMBERS = {
    immutable_orion_tree.INVENTORY_NAME,
    immutable_orion_tree.FREEZE_RECEIPT_NAME,
}
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_RENAME_NOREPLACE = 1
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_ATTEMPT_PATTERN = re.compile(r"[a-z0-9][a-z0-9._-]{0,127}")
_FAILED_LOG_KINDS = ("stderr", "stdout")
_STAGING_ROOT_NAME = "publishable"


class PublicationError(ValueError):
    """Reject an unsafe or incomplete raw-attempt publication."""


def _canonical_existing_directory(value: str | Path, *, label: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        raise PublicationError(f"{label} must be absolute")
    try:
        resolved = path.resolve(strict=True)
        mode = os.lstat(path).st_mode
    except (OSError, RuntimeError) as error:
        raise PublicationError(f"{label} is unavailable: {error}") from error
    if resolved != path:
        raise PublicationError(f"{label} must be canonical without symlink aliases")
    if not stat.S_ISDIR(mode):
        raise PublicationError(f"{label} must be a real directory")
    return path


def _relative_path(value: object, *, label: str) -> str:
    if not isinstance(value, str):
        raise PublicationError(f"{label} must be a string")
    path = PurePosixPath(value)
    if (
        not value
        or path.is_absolute()
        or value != path.as_posix()
        or any(part in ("", ".", "..") for part in path.parts)
    ):
        raise PublicationError(f"{label} must be a canonical relative path")
    return value


def _regular_identity(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _sha256(value: object, *, label: str) -> str:
    if type(value) is not str or _SHA256_PATTERN.fullmatch(value) is None:
        raise PublicationError(f"{label} must be one lowercase SHA-256")
    return value


def _failed_binding(value: object, *, label: str, include_kind: bool) -> dict[str, str]:
    expected = {"path", "sha256", "kind"} if include_kind else {"path", "sha256"}
    if type(value) is not dict or set(value) != expected:
        raise PublicationError(f"{label} schema drifted")
    binding = {
        "path": _relative_path(value["path"], label=f"{label}/path"),
        "sha256": _sha256(value["sha256"], label=f"{label}/sha256"),
    }
    if include_kind:
        kind = value["kind"]
        if type(kind) is not str or kind not in _FAILED_LOG_KINDS:
            raise PublicationError(f"{label}/kind must be one of {_FAILED_LOG_KINDS}")
        binding["kind"] = kind
    return binding


def _validate_failed_manifest_schema(value: object) -> dict[str, Any]:
    label = "failed campaign manifest"
    expected = {
        "schema_version",
        "record_type",
        "attempt_status",
        "attempt_id",
        "failure_logs",
        "retained_artifacts",
    }
    if type(value) is not dict or set(value) != expected:
        raise PublicationError(f"{label} schema drifted")
    if type(value["schema_version"]) is not int or value["schema_version"] != 1:
        raise PublicationError(f"{label}/schema_version must be integer 1")
    if value["record_type"] != FAILED_MANIFEST_RECORD_TYPE:
        raise PublicationError(f"{label}/record_type drifted")
    if value["attempt_status"] != "failed":
        raise PublicationError(f"{label}/attempt_status must be 'failed'")
    attempt_id = value["attempt_id"]
    if type(attempt_id) is not str or _ATTEMPT_PATTERN.fullmatch(attempt_id) is None:
        raise PublicationError(f"{label}/attempt_id is noncanonical")
    if type(value["failure_logs"]) is not list or not value["failure_logs"]:
        raise PublicationError(f"{label}/failure_logs must retain at least one available log")
    if type(value["retained_artifacts"]) is not list:
        raise PublicationError(f"{label}/retained_artifacts must be an array")
    failure_logs = [
        _failed_binding(item, label=f"{label}/failure_logs[{index}]", include_kind=True)
        for index, item in enumerate(value["failure_logs"])
    ]
    retained_artifacts = [
        _failed_binding(
            item,
            label=f"{label}/retained_artifacts[{index}]",
            include_kind=False,
        )
        for index, item in enumerate(value["retained_artifacts"])
    ]
    kinds = [binding["kind"] for binding in failure_logs]
    if kinds != sorted(set(kinds)):
        raise PublicationError(f"{label}/failure_logs must use unique canonical kind order")
    retained_paths = [binding["path"] for binding in retained_artifacts]
    if retained_paths != sorted(set(retained_paths)):
        raise PublicationError(
            f"{label}/retained_artifacts must use unique canonical path order"
        )
    _binding_paths([*failure_logs, *retained_artifacts])
    return {
        "schema_version": 1,
        "record_type": FAILED_MANIFEST_RECORD_TYPE,
        "attempt_status": "failed",
        "attempt_id": attempt_id,
        "failure_logs": failure_logs,
        "retained_artifacts": retained_artifacts,
    }


def _read_anchored_source_member(
    root_fd: int, relative: str, *, expected_sha256: str
) -> tuple[bytes, int]:
    path = PurePosixPath(_relative_path(relative, label="raw source member path"))
    parent_fd = os.dup(root_fd)
    file_fd: int | None = None
    try:
        for part in path.parts[:-1]:
            child_fd = os.open(
                part,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=parent_fd,
            )
            os.close(parent_fd)
            parent_fd = child_fd
        file_fd = os.open(
            path.parts[-1], os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent_fd
        )
        before = os.fstat(file_fd)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
            raise PublicationError(f"raw source member is not one regular file: {relative}")
        payload = bytearray()
        while chunk := os.read(file_fd, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(file_fd)
        if _regular_identity(before) != _regular_identity(after):
            raise PublicationError(f"raw source member changed while reading: {relative}")
        measured = hashlib.sha256(payload).hexdigest()
        if measured != expected_sha256:
            raise PublicationError(f"raw source member SHA-256 drifted: {relative}")
        return bytes(payload), before.st_mode
    except OSError as error:
        raise PublicationError(f"cannot read raw source member {relative}: {error}") from error
    finally:
        if file_fd is not None:
            os.close(file_fd)
        os.close(parent_fd)


def _write_anchored_member(root_fd: int, relative: str, payload: bytes, mode: int) -> None:
    path = PurePosixPath(_relative_path(relative, label="retained member path"))
    parent_fd = os.dup(root_fd)
    file_fd: int | None = None
    try:
        for part in path.parts[:-1]:
            try:
                child_fd = os.open(
                    part,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                    dir_fd=parent_fd,
                )
            except FileNotFoundError:
                os.mkdir(part, mode=0o755, dir_fd=parent_fd)
                os.fsync(parent_fd)
                child_fd = os.open(
                    part,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                    dir_fd=parent_fd,
                )
            os.close(parent_fd)
            parent_fd = child_fd
        file_fd = os.open(
            path.parts[-1],
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
            0o600,
            dir_fd=parent_fd,
        )
        view = memoryview(payload)
        while view:
            written = os.write(file_fd, view)
            if written <= 0:
                raise OSError("short write")
            view = view[written:]
        os.fchmod(file_fd, stat.S_IMODE(mode) | stat.S_IWUSR)
        os.fsync(file_fd)
        os.fsync(parent_fd)
    except OSError as error:
        raise PublicationError(f"cannot publish retained member {relative}: {error}") from error
    finally:
        if file_fd is not None:
            os.close(file_fd)
        os.close(parent_fd)


def _declared_bindings(parsed: Mapping[str, Any]) -> list[dict[str, str]]:
    return [
        parsed["candidate_binding"]["clean_candidate_manifest"],
        *parsed["artifact_bindings"].values(),
        *parsed["products"],
    ]


def _binding_paths(bindings: Sequence[Mapping[str, str]]) -> list[str]:
    paths = [binding["path"] for binding in bindings]
    if len(paths) != len(set(paths)):
        raise PublicationError("campaign manifest contains duplicate retained paths")
    forbidden = {campaign.MANIFEST_NAME, *_RESERVED_SOURCE_MEMBERS}
    if forbidden.intersection(paths):
        raise PublicationError("campaign manifest declares a reserved metadata path")
    return paths


def _declared_paths(parsed: Mapping[str, Any]) -> list[str]:
    return _binding_paths(_declared_bindings(parsed))


def _failed_declared_bindings(parsed: Mapping[str, Any]) -> list[dict[str, str]]:
    return [*parsed["failure_logs"], *parsed["retained_artifacts"]]


def _validate_completed_contract(
    parsed: Mapping[str, Any],
    payloads: Mapping[str, tuple[bytes, int]],
    *,
    destination: Path,
    authorized_pic_root: Path,
) -> None:
    preregistration = parsed["artifact_bindings"]["preregistration"]
    retained_policy = payloads[preregistration["path"]][0]
    if preregistration["sha256"] != campaign.EXPECTED_PREREGISTRATION_SHA256:
        raise PublicationError("completed attempt preregistration binding drifted")
    if retained_policy != campaign.PREREGISTRATION_PATH.read_bytes():
        raise PublicationError("completed attempt preregistration bytes drifted")
    policy = campaign._load_json_bytes(retained_policy, "preregistration")
    if (
        policy.get("record_type") != "q011_section54_qualifying_campaign_preregistration"
        or type(policy.get("schema_version")) is not int
        or policy["schema_version"] != 1
    ):
        raise PublicationError("completed attempt preregistration identity drifted")
    campaign._exact_match(
        campaign._policy_projection(policy),
        campaign._EXPECTED_POLICY_PROJECTION,
        "preregistration projection",
    )
    campaign._validate_identity(
        parsed["run_identity"], parsed["attempt_identity"], policy
    )
    retained = campaign._select_retained_snapshot_products(parsed["products"], policy)
    campaign._select_products(retained, policy)
    campaign._select_stdout_product(parsed["products"])

    class PayloadSnapshot:
        class Member:
            def __init__(self, payload: bytes) -> None:
                self._payload = payload

            def read_bytes(self) -> bytes:
                return self._payload

        def member_path(self, relative: str) -> "PayloadSnapshot.Member":
            try:
                return self.Member(payloads[relative][0])
            except KeyError as error:
                raise ValueError(f"retained raw member is missing: {relative}") from error

    campaign._validate_restart_publications(parsed["products"], PayloadSnapshot(), policy)
    contract_binding = parsed["artifact_bindings"]["attempt_contract"]
    contract = immutable_orion_tree.loads_json_reject_duplicate_keys(
        payloads[contract_binding["path"]][0].decode("utf-8"),
        error_type=PublicationError,
        label="completed attempt contract",
    )
    if type(contract) is not dict:
        raise PublicationError("completed attempt contract must be a JSON object")
    planner_root_value = contract.get("authorized_orion_attempt_root")
    if type(planner_root_value) is not str:
        raise PublicationError("completed attempt contract lacks its planner-authorized root")
    planner_root = Path(planner_root_value)
    if not planner_root.is_absolute() or Path(os.path.abspath(planner_root)) != planner_root:
        raise PublicationError("completed attempt contract uses a noncanonical authorized root")
    try:
        planner_root.relative_to(authorized_pic_root)
    except ValueError as error:
        raise PublicationError("completed attempt contract escaped the authorized PIC root") from error
    # The current planner schema exposes one deterministic per-attempt root. It
    # therefore binds both launch handoff and retained publication. If those
    # roles diverge later, add a distinct retained root to the planner and
    # admission schemas before relaxing this equality.
    if destination != planner_root:
        raise PublicationError(
            "completed attempt destination differs from planner-authorized retained root"
        )


def _snapshot_files(snapshot: Any) -> dict[str, Any]:
    return {
        entry.relative: entry
        for entry in snapshot.entries
        if entry.entry_type == "file"
    }


def _read_source_attempt(
    source_root: Path, *, attempt_status: str
) -> tuple[dict[str, Any], dict[str, tuple[bytes, int]], dict[str, Any]]:
    label = "Q-011 Section 5.4 raw source tree"
    try:
        source_fd = os.open(source_root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    except OSError as error:
        raise PublicationError(f"cannot open raw source tree: {error}") from error
    try:
        _require_same_directory(source_root, source_fd, label=label)
        before = immutable_orion_tree._scan_anchored_tree(
            source_fd, hash_regular=True, error_type=PublicationError, label=label
        )
        files = _snapshot_files(before)
        if _RESERVED_SOURCE_MEMBERS.intersection(files):
            raise PublicationError("raw source tree already contains reserved freeze metadata")
        try:
            manifest_entry = files[campaign.MANIFEST_NAME]
        except KeyError as error:
            raise PublicationError("raw source tree lacks campaign_manifest.json") from error
        manifest_payload, _ = _read_anchored_source_member(
            source_fd,
            campaign.MANIFEST_NAME,
            expected_sha256=manifest_entry.sha256,
        )
        manifest = campaign._load_json_bytes(manifest_payload, "campaign manifest")
        if attempt_status == "completed":
            parsed = campaign._validate_manifest_schema(manifest, source_root)
            declared_bindings = _declared_bindings(parsed)
        else:
            parsed = _validate_failed_manifest_schema(manifest)
            declared_bindings = _failed_declared_bindings(parsed)
        declared_paths = _binding_paths(declared_bindings)
        expected_files = {campaign.MANIFEST_NAME, *declared_paths}
        if set(files) != expected_files:
            raise PublicationError(
                "raw source tree membership differs from the manifest-declared retained files"
            )
        payloads = {}
        expected_by_path = {
            binding["path"]: binding["sha256"] for binding in declared_bindings
        }
        for relative in declared_paths:
            payloads[relative] = _read_anchored_source_member(
                source_fd, relative, expected_sha256=expected_by_path[relative]
            )
        after = immutable_orion_tree._scan_anchored_tree(
            source_fd, hash_regular=True, error_type=PublicationError, label=label
        )
        immutable_orion_tree._require_same_snapshot(
            before,
            after,
            error_type=PublicationError,
            label=label,
            phase="source-local publication copy",
        )
        _require_same_directory(source_root, source_fd, label=label)
        if attempt_status == "failed" and not any(
            payloads[binding["path"]][0] for binding in parsed["failure_logs"]
        ):
            raise PublicationError("failed campaign manifest retains only empty failure logs")
        return manifest, payloads, parsed
    finally:
        os.close(source_fd)


def _require_same_directory(path: Path, descriptor: int, *, label: str) -> None:
    try:
        actual = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise PublicationError(f"{label} is unavailable: {error}") from error
    expected = os.fstat(descriptor)
    if (
        not stat.S_ISDIR(actual.st_mode)
        or (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino)
    ):
        raise PublicationError(f"{label} changed during publication")


def _require_same_directory_at(
    parent_fd: int, name: str, descriptor: int, *, label: str
) -> None:
    if "/" in name:
        raise PublicationError("descriptor-relative directory check received a nested path")
    try:
        actual = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except OSError as error:
        raise PublicationError(f"{label} is unavailable: {error}") from error
    expected = os.fstat(descriptor)
    if (
        not stat.S_ISDIR(actual.st_mode)
        or (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino)
    ):
        raise PublicationError(f"{label} changed during publication")


def _require_absent_at(parent_fd: int, name: str, *, label: str) -> None:
    try:
        os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        return
    except OSError as error:
        raise PublicationError(f"cannot inspect {label}: {error}") from error
    raise PublicationError(f"{label} already exists")


def _rename_no_replace_at(
    source_parent_fd: int,
    source_name: str,
    destination_parent_fd: int,
    destination_name: str,
) -> None:
    if "/" in source_name or "/" in destination_name:
        raise PublicationError("descriptor-relative rename received a nested path")
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
        raise PublicationError(f"retained attempt destination already exists: {destination_name}")
    unsupported = {
        errno.EINVAL,
        errno.ENOSYS,
        getattr(errno, "ENOTSUP", errno.EINVAL),
        getattr(errno, "EOPNOTSUPP", errno.EINVAL),
    }
    if error_number not in unsupported:
        raise OSError(error_number, os.strerror(error_number), destination_name)
    raise PublicationError("retained attempt publication requires atomic no-replace rename support")


def _verify_frozen_tree_at(
    parent_fd: int,
    name: str,
    descriptor: int,
    expected_inventory_sha256: str,
    *,
    runtime_root: Path,
    authorized_destination_root: Path,
    label: str,
) -> dict[str, Any]:
    """Verify exact inventory while retaining the parent-relative root binding."""
    _require_same_directory_at(parent_fd, name, descriptor, label=label)
    verified = immutable_orion_tree._verify_frozen_tree_anchored(
        runtime_root,
        descriptor,
        expected_inventory_sha256,
        authorized_root=authorized_destination_root,
        error_type=PublicationError,
        label="Q-011 Section 5.4 retained campaign attempt",
    )
    _require_same_directory_at(parent_fd, name, descriptor, label=label)
    return verified


def _remove_anchored_tree_at(
    parent_fd: int, name: str, descriptor: int, *, label: str
) -> None:
    """Remove one quarantined tree without reopening any ancestor path."""
    _require_same_directory_at(parent_fd, name, descriptor, label=label)

    def remove_members(directory_fd: int) -> None:
        try:
            status = os.fstat(directory_fd)
            os.fchmod(
                directory_fd,
                stat.S_IMODE(status.st_mode) | stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR,
            )
            names = os.listdir(directory_fd)
        except OSError as error:
            raise PublicationError(f"cannot prepare {label} for removal: {error}") from error
        for member_name in names:
            try:
                observed = os.stat(member_name, dir_fd=directory_fd, follow_symlinks=False)
                if stat.S_ISDIR(observed.st_mode):
                    child_fd = os.open(member_name, _DIRECTORY_FLAGS, dir_fd=directory_fd)
                    try:
                        opened = os.fstat(child_fd)
                        if (observed.st_dev, observed.st_ino) != (
                            opened.st_dev,
                            opened.st_ino,
                        ):
                            raise PublicationError(
                                f"{label} changed during descriptor-relative removal"
                            )
                        remove_members(child_fd)
                        current = os.stat(
                            member_name, dir_fd=directory_fd, follow_symlinks=False
                        )
                        if (current.st_dev, current.st_ino) != (
                            opened.st_dev,
                            opened.st_ino,
                        ):
                            raise PublicationError(
                                f"{label} changed during descriptor-relative removal"
                            )
                    finally:
                        os.close(child_fd)
                    os.rmdir(member_name, dir_fd=directory_fd)
                else:
                    os.unlink(member_name, dir_fd=directory_fd)
            except PublicationError:
                raise
            except OSError as error:
                raise PublicationError(f"cannot remove member from {label}: {error}") from error

    remove_members(descriptor)
    _require_same_directory_at(parent_fd, name, descriptor, label=label)
    try:
        os.rmdir(name, dir_fd=parent_fd)
    except OSError as error:
        raise PublicationError(f"cannot remove {label}: {error}") from error


def _rollback_published_destination(
    parent_fd: int, destination_name: str, descriptor: int
) -> None:
    """Withdraw only the pinned invalid publication, never a path replacement."""
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
            raise PublicationError(
                "cannot locate raced retained attempt publication for quarantine"
            )
        _rename_no_replace_at(parent_fd, candidates[0], parent_fd, rollback_name)
        os.fsync(parent_fd)
        _require_same_directory_at(
            parent_fd,
            rollback_name,
            descriptor,
            label="retained attempt rollback tree",
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
                parent_fd,
                rollback_name,
                descriptor,
                label="retained attempt rollback tree",
            )
        except PublicationError:
            raced_replacement = True
            _rename_no_replace_at(parent_fd, rollback_name, parent_fd, destination_name)
            os.fsync(parent_fd)
            quarantine_moved_original()
    if not raced_replacement:
        _require_absent_at(
            parent_fd,
            destination_name,
            label="invalid retained attempt destination after rollback",
        )
    rollback_fd = os.open(rollback_name, _DIRECTORY_FLAGS, dir_fd=parent_fd)
    try:
        _remove_anchored_tree_at(
            parent_fd,
            rollback_name,
            rollback_fd,
            label="retained attempt rollback tree",
        )
    finally:
        os.close(rollback_fd)
    os.fsync(parent_fd)
    if raced_replacement:
        raise PublicationError(
            "retained attempt rollback rejected a substituted public destination"
        )


def _cleanup_staging(parent_fd: int, name: str, descriptor: int) -> None:
    """Best-effort removal of the pinned staging inode without reopening its path."""
    try:
        _remove_anchored_tree_at(
            parent_fd,
            name,
            descriptor,
            label="retained attempt staging tree",
        )
    except (OSError, PublicationError):
        return


def _cleanup_private_container(parent_fd: int, name: str, descriptor: int) -> None:
    """Best-effort removal of one empty pinned private staging container."""
    try:
        _require_same_directory_at(
            parent_fd,
            name,
            descriptor,
            label="retained attempt private staging container",
        )
        os.rmdir(name, dir_fd=parent_fd)
        os.fsync(parent_fd)
    except (OSError, PublicationError):
        return


def _validate_staged_completed_admission(
    staging: Path,
    final_destination: Path,
    expected_inventory_sha256: str,
    *,
    authorized_destination_root: Path,
    authorized_pic_root: Path,
) -> dict[str, Any]:
    """Run every campaign admission gate before publishing the intended final path."""
    root = immutable_orion_tree.authorized_tree_root(
        staging,
        authorized_root=authorized_destination_root,
        error_type=campaign.QualificationError,
        label="Q-011 Section 5.4 staged campaign root",
    )
    with immutable_orion_tree.staged_verified_frozen_tree(
        root,
        expected_inventory_sha256,
        authorized_root=authorized_destination_root,
        error_type=campaign.QualificationError,
        label="Q-011 Section 5.4 staged campaign tree",
    ) as (tree_report, snapshot):
        campaign._validate_freeze_receipt_semantics(tree_report["freeze_receipt"])
        manifest = campaign._load_json_bytes(
            snapshot.member_path(campaign.MANIFEST_NAME).read_bytes(),
            "campaign manifest",
        )
        parsed = campaign._validate_manifest_schema(manifest, final_destination)
        campaign._validate_declared_tree_members(parsed, snapshot)
        frozen_candidate, external_candidate_closure = campaign._validate_bound_files(
            parsed, snapshot
        )
        policy = campaign._load_bound_policy(
            snapshot, parsed["artifact_bindings"]["preregistration"]
        )
        campaign._validate_identity(parsed["run_identity"], parsed["attempt_identity"], policy)
        # Keep this sequence aligned with campaign._admit_campaign(). The manifest
        # is checked against its intended final path while the tree is still hidden.
        retained_attempt_snapshot = (
            campaign._validated_retained_attempt_semantics_snapshot(
                parsed,
                snapshot,
                policy,
                authorized_pic_root=authorized_pic_root,
            )
        )
        with retained_attempt_snapshot as retained_attempt_semantics:
            return campaign._complete_campaign_admission(
                parsed,
                snapshot,
                policy,
                retained_attempt_semantics,
                frozen_candidate=frozen_candidate,
                external_candidate_closure=external_candidate_closure,
                tree_report=tree_report,
            )


def run_final_admission_self_check(
    campaign_root: str | Path,
    expected_inventory_sha256: str,
    *,
    authorized_destination_root: str | Path,
    authorized_pic_root: str | Path | None = None,
) -> dict[str, Any]:
    """Run the existing analyzer only for one completed immutable attempt."""
    authorized = _canonical_existing_directory(
        authorized_destination_root, label="authorized destination root"
    )
    pic_root = _canonical_existing_directory(
        campaign.ORION_BULK_ROOT if authorized_pic_root is None else authorized_pic_root,
        label="authorized PIC root",
    )
    report = immutable_orion_tree.verify_frozen_tree(
        campaign_root,
        expected_inventory_sha256,
        authorized_root=authorized,
        error_type=PublicationError,
        label="Q-011 Section 5.4 retained campaign attempt",
    )
    if report["freeze_receipt"] != _COMPLETED_RECEIPT:
        raise PublicationError("final admission self-check requires a completed attempt receipt")
    result = campaign.qualify_campaign(
        campaign_root,
        expected_inventory_sha256,
        authorized_orion_root=pic_root,
    )
    if not result["admitted_for_follow_on_numerical_qualification"]:
        raise PublicationError("completed attempt failed final admission self-check")
    return result


def freeze_campaign_attempt(
    source_root: str | Path,
    destination_parent: str | Path,
    *,
    attempt_status: str,
    run_admission_self_check: bool = False,
    authorized_pic_root: str | Path | None = None,
) -> dict[str, Any]:
    """Exclusively publish and recursively freeze one raw campaign attempt."""
    if attempt_status not in ATTEMPT_STATUSES:
        raise PublicationError(f"attempt status must be one of {ATTEMPT_STATUSES}")
    if run_admission_self_check and attempt_status != "completed":
        raise PublicationError("final admission self-check is only available for completed attempts")
    source = _canonical_existing_directory(source_root, label="raw source root")
    destination_parent_path = _canonical_existing_directory(
        destination_parent, label="destination parent"
    )
    pic_root = _canonical_existing_directory(
        campaign.ORION_BULK_ROOT if authorized_pic_root is None else authorized_pic_root,
        label="authorized PIC root",
    )
    manifest, payloads, source_parsed = _read_source_attempt(
        source, attempt_status=attempt_status
    )
    attempt_id = (
        manifest["run_identity"]["attempt_id"]
        if attempt_status == "completed"
        else source_parsed["attempt_id"]
    )
    destination = destination_parent_path / attempt_id
    retained_manifest = copy.deepcopy(manifest)
    if attempt_status == "completed":
        retained_manifest["authorized_orion_campaign_root"] = str(destination)
        parsed = campaign._validate_manifest_schema(retained_manifest, destination)
        _declared_paths(parsed)
    else:
        parsed = _validate_failed_manifest_schema(retained_manifest)
    if attempt_status == "completed":
        _validate_completed_contract(
            parsed,
            payloads,
            destination=destination,
            authorized_pic_root=pic_root,
        )

    parent_fd = os.open(destination_parent_path, _DIRECTORY_FLAGS)
    private_container = destination_parent_path / f".{destination.name}.staging-{uuid.uuid4()}"
    staging = private_container / _STAGING_ROOT_NAME
    private_fd: int | None = None
    staging_fd: int | None = None
    renamed = False
    admission_self_check = None
    try:
        _require_same_directory(
            destination_parent_path, parent_fd, label="retained attempt destination parent"
        )
        _require_absent_at(parent_fd, destination.name, label="retained attempt destination")
        os.mkdir(private_container.name, mode=0o700, dir_fd=parent_fd)
        os.fsync(parent_fd)
        private_fd = os.open(private_container.name, _DIRECTORY_FLAGS, dir_fd=parent_fd)
        os.mkdir(staging.name, mode=0o700, dir_fd=private_fd)
        os.fsync(private_fd)
        staging_fd = os.open(staging.name, _DIRECTORY_FLAGS, dir_fd=private_fd)
        _require_same_directory(staging, staging_fd, label="retained attempt staging tree")
        for relative in sorted(payloads):
            payload, mode = payloads[relative]
            _write_anchored_member(staging_fd, relative, payload, mode)
        _write_anchored_member(
            staging_fd,
            campaign.MANIFEST_NAME,
            (json.dumps(retained_manifest, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
                "utf-8"
            ),
            0o644,
        )
        _require_same_directory(staging, staging_fd, label="retained attempt staging tree")
        receipt = _COMPLETED_RECEIPT if attempt_status == "completed" else _FAILED_RECEIPT
        frozen = immutable_orion_tree.freeze_tree_anchored(
            staging,
            staging_fd,
            receipt,
            authorized_root=destination_parent_path,
            error_type=PublicationError,
            label="Q-011 Section 5.4 retained campaign attempt",
        )
        _require_same_directory(staging, staging_fd, label="retained attempt staging tree")
        if attempt_status == "completed":
            try:
                _validate_staged_completed_admission(
                    staging,
                    destination,
                    frozen["inventory_sha256"],
                    authorized_destination_root=destination_parent_path,
                    authorized_pic_root=pic_root,
                )
            except (campaign.QualificationError, OSError, ValueError) as error:
                raise PublicationError(
                    "completed attempt failed prepublication campaign admission"
                ) from error
        _require_same_directory(
            destination_parent_path, parent_fd, label="retained attempt destination parent"
        )
        _require_same_directory(staging, staging_fd, label="retained attempt staging tree")
        _require_absent_at(parent_fd, destination.name, label="retained attempt destination")
        verified = _verify_frozen_tree_at(
            private_fd,
            staging.name,
            staging_fd,
            frozen["inventory_sha256"],
            runtime_root=staging,
            authorized_destination_root=destination_parent_path,
            label="retained attempt staging tree",
        )
        # Orion rejects cross-parent rename of a read-only directory. Descendants
        # remain frozen; make the root owner-write-only and non-traversable for
        # the rename, then restore its exact frozen mode through the pinned fd.
        root_mode = stat.S_IMODE(os.fstat(staging_fd).st_mode)
        os.fchmod(staging_fd, stat.S_IWUSR)
        os.fsync(staging_fd)
        _rename_no_replace_at(private_fd, staging.name, parent_fd, destination.name)
        renamed = True
        os.fchmod(staging_fd, root_mode)
        os.fsync(staging_fd)
        os.fsync(private_fd)
        os.fsync(parent_fd)
        verified = _verify_frozen_tree_at(
            parent_fd,
            destination.name,
            staging_fd,
            frozen["inventory_sha256"],
            runtime_root=destination,
            authorized_destination_root=destination_parent_path,
            label="retained attempt published tree",
        )
        if run_admission_self_check:
            admission_self_check = run_final_admission_self_check(
                destination,
                frozen["inventory_sha256"],
                authorized_destination_root=destination_parent_path,
                authorized_pic_root=pic_root,
            )
    except BaseException:
        if renamed:
            try:
                _rollback_published_destination(parent_fd, destination.name, staging_fd)
            except BaseException as rollback_error:
                raise PublicationError(
                    "cannot remove invalid retained attempt destination"
                ) from rollback_error
        elif private_fd is not None and staging_fd is not None:
            _cleanup_staging(private_fd, staging.name, staging_fd)
        raise
    finally:
        if staging_fd is not None:
            os.close(staging_fd)
        if private_fd is not None:
            _cleanup_private_container(parent_fd, private_container.name, private_fd)
            os.close(private_fd)
        os.close(parent_fd)
    result = {
        "schema_version": 1,
        "record_type": RESULT_RECORD_TYPE,
        "attempt_status": attempt_status,
        "campaign_root": str(destination),
        "inventory_sha256": frozen["inventory_sha256"],
        "recursively_read_only": verified["recursively_read_only"],
        "prepublication_admission_passed": attempt_status == "completed",
        "admission_self_check_available": attempt_status == "completed",
        "admission_self_check": admission_self_check,
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_root")
    parser.add_argument("destination_parent")
    parser.add_argument("attempt_status", choices=ATTEMPT_STATUSES)
    parser.add_argument("--run-final-admission-self-check", action="store_true")
    parser.add_argument("--authorized-pic-root", default=str(campaign.ORION_BULK_ROOT))
    args = parser.parse_args(argv)
    result = freeze_campaign_attempt(
        args.source_root,
        args.destination_parent,
        attempt_status=args.attempt_status,
        run_admission_self_check=args.run_final_admission_self_check,
        authorized_pic_root=args.authorized_pic_root,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
