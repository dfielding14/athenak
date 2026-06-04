#!/usr/bin/env python3
"""Provision or recover the fixed Q011 pressure-publication acceptance root."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import stat
from typing import Any


AUTHORIZED_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
ACCEPTANCE_DIRECTORY = "publication_acceptance"
AUTHORIZED_UID = 18664
AUTHORIZED_GID = 31114
REVIEWED_PIC_ROOT_IDENTITY = (135357496, 720587399016531307)
REVIEWED_POLICY_ROOT_IDENTITY = (135357496, 720587399016535607)
REVIEWED_ACCEPTANCE_ROOT_IDENTITY = (135357496, 720587400627193972)
REVIEWED_PIC_ROOT_MODE = 0o2755
REVIEWED_POLICY_ROOT_MODE = 0o2755
REVIEWED_PIC_ROOT_XATTRS = {"lustre.lov"}
REVIEWED_POLICY_ROOT_XATTRS = {"lustre.lov"}
RECOVERY_RECEIPT_NAME = (
    "q011_pressure_publication_acceptance_root_recovery_"
    "1554766c-21e2-48b1-8cfe-b1e7e4e75aa2.json"
)
ACL_XATTRS = {"system.posix_acl_access", "system.posix_acl_default"}
ALLOWED_XATTRS = {"lustre.lov"}
REVIEWED_RECOVERY_XATTRS = {"lustre.lov"}
NOFOLLOW_FLAG = getattr(os, "O_NOFOLLOW", None)
DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | (NOFOLLOW_FLAG or 0)
FILE_FLAGS = os.O_RDONLY | (NOFOLLOW_FLAG or 0)
STAGING_PREFIX = ".q011-acceptance-recovery."


class AcceptanceRootError(ValueError):
    """Raised when the acceptance-root checkpoint is not exactly recoverable."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AcceptanceRootError(message)


def _identity_tuple(status: os.stat_result) -> tuple[int, int]:
    return status.st_dev, status.st_ino


def _identity(status: os.stat_result) -> dict[str, int]:
    return {"device": status.st_dev, "inode": status.st_ino}


def _same_identity(first: os.stat_result, second: os.stat_result) -> bool:
    return _identity_tuple(first) == _identity_tuple(second)


def _require_identity(
    status: os.stat_result, expected: tuple[int, int] | None, label: str
) -> None:
    if expected is not None:
        _require(_identity_tuple(status) == expected, f"{label} identity drifted")


def _open_absolute_directory(path: Path) -> int:
    _require(NOFOLLOW_FLAG is not None, "Platform does not provide O_NOFOLLOW")
    absolute = Path(os.path.abspath(path))
    _require(absolute.is_absolute(), f"Directory path is not absolute: {path}")
    descriptor = os.open("/", DIRECTORY_FLAGS)
    try:
        for component in absolute.parts[1:]:
            child = os.open(component, DIRECTORY_FLAGS, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _require_same_absolute_directory(path: Path, descriptor: int, label: str) -> None:
    reopened = _open_absolute_directory(path)
    try:
        actual = os.fstat(reopened)
        retained = os.fstat(descriptor)
        _require(
            stat.S_ISDIR(actual.st_mode)
            and stat.S_ISDIR(retained.st_mode)
            and _same_identity(actual, retained),
            f"{label} changed during acceptance-root provisioning",
        )
    finally:
        os.close(reopened)


def _child_status(
    parent_descriptor: int, name: str, label: str = "Acceptance root"
) -> os.stat_result:
    try:
        status = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError as error:
        raise AcceptanceRootError(f"{label} is unavailable") from error
    _require(stat.S_ISDIR(status.st_mode), f"{label} is not a directory")
    return status


def _require_same_child(
    parent_descriptor: int, name: str, descriptor: int, label: str
) -> None:
    lexical = _child_status(parent_descriptor, name)
    retained = os.fstat(descriptor)
    _require(
        stat.S_ISDIR(retained.st_mode) and _same_identity(lexical, retained),
        f"{label} changed during acceptance-root provisioning",
    )


def _open_retained_directory_child(
    parent_descriptor: int,
    name: str,
    expected_status: os.stat_result,
    label: str,
) -> int:
    descriptor = os.open(name, DIRECTORY_FLAGS, dir_fd=parent_descriptor)
    try:
        retained = os.fstat(descriptor)
        _require(
            _same_identity(expected_status, retained),
            f"{label} changed before retained open",
        )
        _require_same_child(parent_descriptor, name, descriptor, label)
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _open_retained_child(
    parent_descriptor: int, name: str, expected_status: os.stat_result
) -> int:
    return _open_retained_directory_child(
        parent_descriptor, name, expected_status, "Acceptance root"
    )


def _xattrs(descriptor: int) -> set[str]:
    return set(os.listxattr(descriptor))


def _require_xattr_closure(
    descriptor: int,
    *,
    label: str = "Acceptance root",
    allowed_xattrs: set[str] = ALLOWED_XATTRS,
) -> list[str]:
    names = _xattrs(descriptor)
    _require(not names.intersection(ACL_XATTRS), f"{label} retains ACL xattrs")
    unexpected = names.difference(allowed_xattrs)
    _require(
        not unexpected,
        f"{label} has unexpected xattrs: {sorted(unexpected)}",
    )
    return sorted(names)


def _remove_inherited_acl_xattrs(descriptor: int) -> None:
    for name in sorted(_xattrs(descriptor).intersection(ACL_XATTRS)):
        os.removexattr(descriptor, name)


def _state(descriptor: int) -> dict[str, Any]:
    status = os.fstat(descriptor)
    return {
        "identity": _identity(status),
        "uid": status.st_uid,
        "gid": status.st_gid,
        "mode": f"{stat.S_IMODE(status.st_mode):04o}",
        "entries": sorted(os.listdir(descriptor)),
        "xattrs": _require_xattr_closure(descriptor),
    }


def _require_owner_group(
    status: os.stat_result,
    *,
    expected_uid: int,
    expected_gid: int,
    label: str = "Acceptance root",
) -> None:
    _require(status.st_uid == expected_uid, f"{label} owner drifted")
    _require(status.st_gid == expected_gid, f"{label} group drifted")


def _require_directory_metadata(
    descriptor: int,
    *,
    label: str,
    expected_identity: tuple[int, int],
    expected_mode: int,
    expected_xattrs: set[str],
) -> None:
    status = os.fstat(descriptor)
    _require(stat.S_ISDIR(status.st_mode), f"{label} is not a directory")
    _require_identity(status, expected_identity, label)
    _require_owner_group(
        status,
        expected_uid=AUTHORIZED_UID,
        expected_gid=AUTHORIZED_GID,
        label=label,
    )
    _require(
        stat.S_IMODE(status.st_mode) == expected_mode,
        f"{label} mode is not exact {expected_mode:04o}",
    )
    names = _require_xattr_closure(descriptor, label=label)
    _require(
        set(names) == expected_xattrs,
        f"{label} xattr closure differs from reviewed checkpoint",
    )


def _require_state_shape(
    state: dict[str, Any],
    *,
    require_empty: bool,
    expected_xattrs: set[str] | None = None,
) -> None:
    if require_empty:
        _require(not state["entries"], "Acceptance root is not empty")
    if expected_xattrs is not None:
        _require(
            set(state["xattrs"]) == expected_xattrs,
            "Acceptance root xattr closure differs from reviewed checkpoint",
        )


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_json_bytes(value: dict[str, Any]) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def provision_or_verify(
    *,
    pic_root: Path,
    expected_uid: int,
    expected_gid: int,
    validated_source_commit: str,
    helper_sha256: str,
    expected_pic_root_identity: tuple[int, int] | None = None,
    expected_recovery_identity: tuple[int, int] | None = None,
    expected_recovery_xattrs: set[str] | None = None,
    recover_exact_empty_inherited_setgid_root: bool = False,
    reconcile_exact_empty_normalized_root: bool = False,
    require_empty_recovery_receipt_namespace: bool = False,
) -> dict[str, Any]:
    """Create, verify, recover, or reconcile the fixed sibling acceptance root."""
    _require(
        not (
            recover_exact_empty_inherited_setgid_root
            and reconcile_exact_empty_normalized_root
        ),
        "Select at most one acceptance-root recovery action",
    )
    if require_empty_recovery_receipt_namespace:
        _require_empty_recovery_receipt_namespace()
    pic_root = Path(os.path.abspath(pic_root))
    acceptance_root = pic_root / ACCEPTANCE_DIRECTORY
    pic_descriptor = _open_absolute_directory(pic_root)
    child_descriptor: int | None = None
    try:
        pic_status = os.fstat(pic_descriptor)
        _require_identity(pic_status, expected_pic_root_identity, "Authorized PIC root")
        _require_same_absolute_directory(pic_root, pic_descriptor, "Authorized PIC root")
        created = False
        try:
            child_status = os.stat(
                ACCEPTANCE_DIRECTORY,
                dir_fd=pic_descriptor,
                follow_symlinks=False,
            )
        except FileNotFoundError:
            _require(
                expected_recovery_identity is None,
                "Bound acceptance root is unavailable",
            )
            _require(
                not (
                    recover_exact_empty_inherited_setgid_root
                    or reconcile_exact_empty_normalized_root
                ),
                "Recovery requires the retained acceptance-root checkpoint",
            )
            os.mkdir(ACCEPTANCE_DIRECTORY, mode=0o700, dir_fd=pic_descriptor)
            created = True
            child_status = _child_status(pic_descriptor, ACCEPTANCE_DIRECTORY)
        _require(stat.S_ISDIR(child_status.st_mode), "Acceptance root is not a directory")
        child_descriptor = _open_retained_child(
            pic_descriptor, ACCEPTANCE_DIRECTORY, child_status
        )
        status = os.fstat(child_descriptor)
        _require_owner_group(status, expected_uid=expected_uid, expected_gid=expected_gid)
        if created:
            _remove_inherited_acl_xattrs(child_descriptor)
        before = _state(child_descriptor)
        mode = stat.S_IMODE(status.st_mode)
        if created:
            _require_state_shape(before, require_empty=True)
            os.fchmod(child_descriptor, 0o700)
            action = "fresh_created_and_normalized"
        elif recover_exact_empty_inherited_setgid_root:
            _require_identity(
                status, expected_recovery_identity, "Reviewed acceptance root"
            )
            _require(
                mode == 0o2700,
                "Recovery requires exact inherited-setgid mode 02700",
            )
            _require_state_shape(
                before,
                require_empty=True,
                expected_xattrs=expected_recovery_xattrs,
            )
            os.fchmod(child_descriptor, 0o700)
            action = "recovered_exact_empty_inherited_setgid_root"
        elif reconcile_exact_empty_normalized_root:
            _require_identity(
                status, expected_recovery_identity, "Reviewed acceptance root"
            )
            _require(mode == 0o700, "Reconciliation requires exact normalized mode 0700")
            _require_state_shape(
                before,
                require_empty=True,
                expected_xattrs=expected_recovery_xattrs,
            )
            action = "reconciled_exact_empty_normalized_root"
        else:
            _require_identity(
                status, expected_recovery_identity, "Reviewed acceptance root"
            )
            _require(mode == 0o700, "Acceptance root mode is not exact 0700")
            _require_state_shape(
                before,
                require_empty=True,
                expected_xattrs=expected_recovery_xattrs,
            )
            action = "verified_existing_empty_root"
        _require_same_absolute_directory(pic_root, pic_descriptor, "Authorized PIC root")
        _require_same_child(
            pic_descriptor,
            ACCEPTANCE_DIRECTORY,
            child_descriptor,
            "Acceptance root",
        )
        status = os.fstat(child_descriptor)
        _require_identity(status, _identity_tuple(child_status), "Acceptance root")
        _require_owner_group(status, expected_uid=expected_uid, expected_gid=expected_gid)
        after = _state(child_descriptor)
        _require(after["mode"] == "0700", "Acceptance root mode is not exact 0700")
        _require_state_shape(
            after,
            require_empty=True,
            expected_xattrs=(
                expected_recovery_xattrs
                if action
                in {
                    "recovered_exact_empty_inherited_setgid_root",
                    "reconciled_exact_empty_normalized_root",
                    "verified_existing_empty_root",
                }
                else None
            ),
        )
        os.fsync(child_descriptor)
        os.fsync(pic_descriptor)
        _require_same_absolute_directory(pic_root, pic_descriptor, "Authorized PIC root")
        _require_same_child(
            pic_descriptor,
            ACCEPTANCE_DIRECTORY,
            child_descriptor,
            "Acceptance root",
        )
        status = os.fstat(child_descriptor)
        _require_identity(status, _identity_tuple(child_status), "Acceptance root")
        _require_owner_group(status, expected_uid=expected_uid, expected_gid=expected_gid)
        final = _state(child_descriptor)
        _require(final == after, "Acceptance root changed after durable sync")
        return {
            "schema_version": 1,
            "record_type": "q011_pressure_publication_acceptance_root_checkpoint",
            "action": action,
            "validated_source_commit": validated_source_commit,
            "helper_sha256": helper_sha256,
            "pic_root": str(pic_root),
            "pic_root_identity": _identity(pic_status),
            "acceptance_root": str(acceptance_root),
            "descriptor_relative_operation": True,
            "before": before,
            "after": after,
            "child_and_parent_synced": True,
        }
    finally:
        if child_descriptor is not None:
            os.close(child_descriptor)
        os.close(pic_descriptor)


def _require_same_regular_file_at(
    parent_descriptor: int, name: str, descriptor: int
) -> None:
    lexical = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    retained = os.fstat(descriptor)
    _require(
        stat.S_ISREG(lexical.st_mode)
        and stat.S_ISREG(retained.st_mode)
        and _same_identity(lexical, retained),
        "Recovery receipt changed during publication",
    )


def _read_stable_regular_file_descriptor(descriptor: int) -> bytes:
    before = os.fstat(descriptor)
    _require(stat.S_ISREG(before.st_mode), "Recovery receipt is not a regular file")
    os.lseek(descriptor, 0, os.SEEK_SET)
    with os.fdopen(descriptor, "rb", closefd=False) as stream:
        payload = stream.read()
    after = os.fstat(descriptor)
    _require(_same_identity(before, after), "Recovery receipt identity drifted")
    _require(
        before.st_size == after.st_size == len(payload),
        "Recovery receipt size drifted",
    )
    return payload


def _read_regular_file_at(parent_descriptor: int, name: str) -> bytes:
    descriptor = os.open(name, FILE_FLAGS, dir_fd=parent_descriptor)
    try:
        payload = _read_stable_regular_file_descriptor(descriptor)
        _require_same_regular_file_at(parent_descriptor, name, descriptor)
        return payload
    finally:
        os.close(descriptor)


def _expected_identity(value: tuple[int, int]) -> dict[str, int]:
    return {"device": value[0], "inode": value[1]}


def validate_recovery_receipt(
    payload: bytes, *, validated_source_commit: str, helper_sha256: str
) -> dict[str, Any]:
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise AcceptanceRootError("Recovery receipt is not canonical JSON") from error
    _require(isinstance(value, dict), "Recovery receipt is not an object")
    _require(canonical_json_bytes(value) == payload, "Recovery receipt is not canonical")
    action = value.get("action")
    _require(
        action
        in {
            "recovered_exact_empty_inherited_setgid_root",
            "reconciled_exact_empty_normalized_root",
        },
        "Recovery receipt action is not authorized",
    )
    before_mode = (
        "2700"
        if action == "recovered_exact_empty_inherited_setgid_root"
        else "0700"
    )
    expected_state = {
        "identity": _expected_identity(REVIEWED_ACCEPTANCE_ROOT_IDENTITY),
        "uid": AUTHORIZED_UID,
        "gid": AUTHORIZED_GID,
        "mode": "0700",
        "entries": [],
        "xattrs": sorted(REVIEWED_RECOVERY_XATTRS),
    }
    expected = {
        "schema_version": 1,
        "record_type": "q011_pressure_publication_acceptance_root_checkpoint",
        "action": action,
        "validated_source_commit": validated_source_commit,
        "helper_sha256": helper_sha256,
        "pic_root": str(AUTHORIZED_PIC_ROOT),
        "pic_root_identity": _expected_identity(REVIEWED_PIC_ROOT_IDENTITY),
        "acceptance_root": str(AUTHORIZED_PIC_ROOT / ACCEPTANCE_DIRECTORY),
        "descriptor_relative_operation": True,
        "before": {**expected_state, "mode": before_mode},
        "after": expected_state,
        "child_and_parent_synced": True,
    }
    _require(value == expected, "Recovery receipt differs from exact reviewed schema")
    return value


def _policy_directory() -> Path:
    return AUTHORIZED_PIC_ROOT / "policy"


def _require_policy_descriptors(pic_descriptor: int, policy_descriptor: int) -> None:
    _require_directory_metadata(
        pic_descriptor,
        label="Authorized PIC root",
        expected_identity=REVIEWED_PIC_ROOT_IDENTITY,
        expected_mode=REVIEWED_PIC_ROOT_MODE,
        expected_xattrs=REVIEWED_PIC_ROOT_XATTRS,
    )
    _require_directory_metadata(
        policy_descriptor,
        label="PIC policy root",
        expected_identity=REVIEWED_POLICY_ROOT_IDENTITY,
        expected_mode=REVIEWED_POLICY_ROOT_MODE,
        expected_xattrs=REVIEWED_POLICY_ROOT_XATTRS,
    )
    _require_same_absolute_directory(
        AUTHORIZED_PIC_ROOT, pic_descriptor, "Authorized PIC root"
    )
    _require_same_child(pic_descriptor, "policy", policy_descriptor, "PIC policy root")


def _require_live_acceptance_root(pic_descriptor: int) -> None:
    status = _child_status(pic_descriptor, ACCEPTANCE_DIRECTORY)
    descriptor = _open_retained_child(pic_descriptor, ACCEPTANCE_DIRECTORY, status)
    try:
        _require_directory_metadata(
            descriptor,
            label="Acceptance root",
            expected_identity=REVIEWED_ACCEPTANCE_ROOT_IDENTITY,
            expected_mode=0o700,
            expected_xattrs=REVIEWED_RECOVERY_XATTRS,
        )
        state = _state(descriptor)
        _require_state_shape(
            state,
            require_empty=True,
            expected_xattrs=REVIEWED_RECOVERY_XATTRS,
        )
        _require_same_child(
            pic_descriptor,
            ACCEPTANCE_DIRECTORY,
            descriptor,
            "Acceptance root",
        )
    finally:
        os.close(descriptor)


def _recovery_staging_aliases(policy_descriptor: int) -> list[str]:
    return [
        name
        for name in sorted(os.listdir(policy_descriptor))
        if name.startswith(STAGING_PREFIX)
    ]


def _require_no_recovery_staging_aliases(policy_descriptor: int) -> None:
    _require(
        not _recovery_staging_aliases(policy_descriptor),
        "Recovery receipt staging-prefix namespace is not empty",
    )


def _require_exact_recovery_staging_namespace(
    policy_descriptor: int, staging_name: str
) -> None:
    _require(
        _recovery_staging_aliases(policy_descriptor) == [staging_name],
        "Recovery receipt staging-prefix namespace is not exact",
    )


def _require_empty_recovery_receipt_namespace() -> None:
    pic_descriptor, policy_descriptor = _open_policy_descriptors()
    try:
        _require_policy_descriptors(pic_descriptor, policy_descriptor)
        names = set(os.listdir(policy_descriptor))
        _require(
            RECOVERY_RECEIPT_NAME not in names,
            "Published recovery receipt already exists",
        )
        _require_no_recovery_staging_aliases(policy_descriptor)
        _require_policy_descriptors(pic_descriptor, policy_descriptor)
    finally:
        os.close(policy_descriptor)
        os.close(pic_descriptor)


def _open_policy_descriptors() -> tuple[int, int]:
    pic_descriptor = _open_absolute_directory(AUTHORIZED_PIC_ROOT)
    try:
        pic_status = os.fstat(pic_descriptor)
        _require_identity(pic_status, REVIEWED_PIC_ROOT_IDENTITY, "Authorized PIC root")
        _require_same_absolute_directory(
            AUTHORIZED_PIC_ROOT, pic_descriptor, "Authorized PIC root"
        )
        policy_status = _child_status(pic_descriptor, "policy", "PIC policy root")
        policy_descriptor = _open_retained_directory_child(
            pic_descriptor, "policy", policy_status, "PIC policy root"
        )
        try:
            _require_policy_descriptors(pic_descriptor, policy_descriptor)
            return pic_descriptor, policy_descriptor
        except BaseException:
            os.close(policy_descriptor)
            raise
    except BaseException:
        os.close(pic_descriptor)
        raise


def _require_receipt_file_shape(
    descriptor: int, *, expected_nlink: int, expected_mode: int | None
) -> None:
    status = os.fstat(descriptor)
    _require(stat.S_ISREG(status.st_mode), "Recovery receipt is not a regular file")
    _require(status.st_uid == AUTHORIZED_UID, "Recovery receipt owner drifted")
    _require(status.st_gid == AUTHORIZED_GID, "Recovery receipt group drifted")
    _require(status.st_nlink == expected_nlink, "Recovery receipt link count drifted")
    _require_xattr_closure(descriptor, label="Recovery receipt")
    if expected_mode is None:
        return
    _require(
        stat.S_IMODE(status.st_mode) == expected_mode,
        f"Recovery receipt mode is not exact {expected_mode:04o}",
    )


def _staging_name(path: Path) -> str:
    absolute = Path(os.path.abspath(path))
    _require(absolute.parent == _policy_directory(), "Recovery staging path drifted")
    _require(
        absolute.name.startswith(".q011-acceptance-recovery."),
        "Recovery staging name drifted",
    )
    return absolute.name


def publish_recovery_receipt(
    staging_path: Path, *, validated_source_commit: str, helper_sha256: str
) -> dict[str, Any]:
    pic_descriptor, policy_descriptor = _open_policy_descriptors()
    try:
        _require_policy_descriptors(pic_descriptor, policy_descriptor)
        _require_live_acceptance_root(pic_descriptor)
        staging_name = _staging_name(staging_path)
        _require_exact_recovery_staging_namespace(policy_descriptor, staging_name)
        staging_descriptor = os.open(staging_name, FILE_FLAGS, dir_fd=policy_descriptor)
        try:
            payload = _read_stable_regular_file_descriptor(staging_descriptor)
            _require_same_regular_file_at(
                policy_descriptor, staging_name, staging_descriptor
            )
            _require_receipt_file_shape(
                staging_descriptor, expected_nlink=1, expected_mode=None
            )
            value = validate_recovery_receipt(
                payload,
                validated_source_commit=validated_source_commit,
                helper_sha256=helper_sha256,
            )
            os.fchmod(staging_descriptor, 0o400)
            os.fsync(staging_descriptor)
            _require_receipt_file_shape(
                staging_descriptor, expected_nlink=1, expected_mode=0o400
            )
            _require_same_regular_file_at(
                policy_descriptor, staging_name, staging_descriptor
            )
            _require_policy_descriptors(pic_descriptor, policy_descriptor)
            _require_live_acceptance_root(pic_descriptor)
            _require_exact_recovery_staging_namespace(policy_descriptor, staging_name)
            os.link(
                staging_name,
                RECOVERY_RECEIPT_NAME,
                src_dir_fd=policy_descriptor,
                dst_dir_fd=policy_descriptor,
                follow_symlinks=False,
            )
            final_descriptor = os.open(
                RECOVERY_RECEIPT_NAME, FILE_FLAGS, dir_fd=policy_descriptor
            )
            try:
                _require(
                    _same_identity(
                        os.fstat(staging_descriptor), os.fstat(final_descriptor)
                    ),
                    "Published recovery receipt identity drifted",
                )
                _require_receipt_file_shape(
                    final_descriptor, expected_nlink=2, expected_mode=0o400
                )
            finally:
                os.close(final_descriptor)
            os.unlink(staging_name, dir_fd=policy_descriptor)
            os.fsync(policy_descriptor)
            _require_policy_descriptors(pic_descriptor, policy_descriptor)
            _require_live_acceptance_root(pic_descriptor)
            _require_no_recovery_staging_aliases(policy_descriptor)
            _require_receipt_file_shape(
                staging_descriptor, expected_nlink=1, expected_mode=0o400
            )
            _require_same_regular_file_at(
                policy_descriptor, RECOVERY_RECEIPT_NAME, staging_descriptor
            )
            retained = _read_regular_file_at(policy_descriptor, RECOVERY_RECEIPT_NAME)
            _require(retained == payload, "Published recovery receipt bytes drifted")
            _require_same_regular_file_at(
                policy_descriptor, RECOVERY_RECEIPT_NAME, staging_descriptor
            )
            validate_recovery_receipt(
                retained,
                validated_source_commit=validated_source_commit,
                helper_sha256=helper_sha256,
            )
            _require_live_acceptance_root(pic_descriptor)
            return value
        finally:
            os.close(staging_descriptor)
    finally:
        os.close(policy_descriptor)
        os.close(pic_descriptor)


def verify_recovery_receipt(
    path: Path, *, validated_source_commit: str, helper_sha256: str
) -> dict[str, Any]:
    absolute = Path(os.path.abspath(path))
    _require(
        absolute == _policy_directory() / RECOVERY_RECEIPT_NAME,
        "Published recovery receipt path drifted",
    )
    pic_descriptor, policy_descriptor = _open_policy_descriptors()
    try:
        _require_policy_descriptors(pic_descriptor, policy_descriptor)
        _require_live_acceptance_root(pic_descriptor)
        _require_no_recovery_staging_aliases(policy_descriptor)
        descriptor = os.open(RECOVERY_RECEIPT_NAME, FILE_FLAGS, dir_fd=policy_descriptor)
        try:
            payload = _read_stable_regular_file_descriptor(descriptor)
            _require_same_regular_file_at(
                policy_descriptor, RECOVERY_RECEIPT_NAME, descriptor
            )
            _require_receipt_file_shape(
                descriptor, expected_nlink=1, expected_mode=0o400
            )
            value = validate_recovery_receipt(
                payload,
                validated_source_commit=validated_source_commit,
                helper_sha256=helper_sha256,
            )
            os.fsync(descriptor)
            os.fsync(policy_descriptor)
            _require_policy_descriptors(pic_descriptor, policy_descriptor)
            _require_live_acceptance_root(pic_descriptor)
            _require_no_recovery_staging_aliases(policy_descriptor)
            retained = _read_stable_regular_file_descriptor(descriptor)
            _require(retained == payload, "Published recovery receipt bytes drifted")
            _require_same_regular_file_at(
                policy_descriptor, RECOVERY_RECEIPT_NAME, descriptor
            )
            _require_receipt_file_shape(
                descriptor, expected_nlink=1, expected_mode=0o400
            )
            return value
        finally:
            os.close(descriptor)
    finally:
        os.close(policy_descriptor)
        os.close(pic_descriptor)


def reconcile_linked_recovery_receipt(
    *, validated_source_commit: str, helper_sha256: str
) -> dict[str, Any]:
    """Remove one exact retained staging alias after an interrupted hard-link publish."""
    pic_descriptor, policy_descriptor = _open_policy_descriptors()
    try:
        _require_policy_descriptors(pic_descriptor, policy_descriptor)
        _require_live_acceptance_root(pic_descriptor)
        descriptor = os.open(RECOVERY_RECEIPT_NAME, FILE_FLAGS, dir_fd=policy_descriptor)
        try:
            payload = _read_stable_regular_file_descriptor(descriptor)
            _require_same_regular_file_at(
                policy_descriptor, RECOVERY_RECEIPT_NAME, descriptor
            )
            _require_receipt_file_shape(
                descriptor, expected_nlink=2, expected_mode=0o400
            )
            value = validate_recovery_receipt(
                payload,
                validated_source_commit=validated_source_commit,
                helper_sha256=helper_sha256,
            )
            receipt_status = os.fstat(descriptor)
            aliases = _recovery_staging_aliases(policy_descriptor)
            _require(
                len(aliases) == 1,
                "Recovery receipt staging-prefix namespace is not exact",
            )
            alias_status = os.stat(
                aliases[0], dir_fd=policy_descriptor, follow_symlinks=False
            )
            _require(
                stat.S_ISREG(alias_status.st_mode)
                and _same_identity(receipt_status, alias_status),
                "Recovery receipt does not have one exact retained staging alias",
            )
            os.unlink(aliases[0], dir_fd=policy_descriptor)
            os.fsync(descriptor)
            os.fsync(policy_descriptor)
            _require_policy_descriptors(pic_descriptor, policy_descriptor)
            _require_live_acceptance_root(pic_descriptor)
            _require_no_recovery_staging_aliases(policy_descriptor)
            retained = _read_stable_regular_file_descriptor(descriptor)
            _require(retained == payload, "Published recovery receipt bytes drifted")
            _require_same_regular_file_at(
                policy_descriptor, RECOVERY_RECEIPT_NAME, descriptor
            )
            _require_receipt_file_shape(
                descriptor, expected_nlink=1, expected_mode=0o400
            )
            return value
        finally:
            os.close(descriptor)
    finally:
        os.close(policy_descriptor)
        os.close(pic_descriptor)


def _self_authenticate(
    *, validated_source_commit: str, expected_helper_sha256: str
) -> str:
    _require(
        re.fullmatch(r"[0-9a-f]{40}", validated_source_commit) is not None,
        "Validated source commit is not one lowercase Git SHA-1",
    )
    _require(
        re.fullmatch(r"[0-9a-f]{64}", expected_helper_sha256) is not None,
        "Expected helper SHA-256 is malformed",
    )
    actual = _sha256(Path(__file__).read_bytes())
    _require(actual == expected_helper_sha256, "Executing helper SHA-256 drifted")
    return actual


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "--recover-exact-empty-inherited-setgid-root",
        action="store_true",
        help="Normalize only the retained exact empty 02700 checkpoint.",
    )
    action.add_argument(
        "--reconcile-exact-empty-normalized-root",
        action="store_true",
        help="Reconcile an interrupted recovery only after exact 0700 normalization.",
    )
    action.add_argument("--publish-recovery-receipt", type=Path)
    action.add_argument("--verify-recovery-receipt", type=Path)
    action.add_argument("--reconcile-linked-recovery-receipt", action="store_true")
    parser.add_argument("--validated-source-commit", required=True)
    parser.add_argument("--expected-helper-sha256", required=True)
    args = parser.parse_args(argv)
    helper_sha256 = _self_authenticate(
        validated_source_commit=args.validated_source_commit,
        expected_helper_sha256=args.expected_helper_sha256,
    )
    common = {
        "validated_source_commit": args.validated_source_commit,
        "helper_sha256": helper_sha256,
    }
    if args.publish_recovery_receipt is not None:
        result = publish_recovery_receipt(args.publish_recovery_receipt, **common)
    elif args.verify_recovery_receipt is not None:
        result = verify_recovery_receipt(args.verify_recovery_receipt, **common)
    elif args.reconcile_linked_recovery_receipt:
        result = reconcile_linked_recovery_receipt(**common)
    else:
        result = provision_or_verify(
            pic_root=AUTHORIZED_PIC_ROOT,
            expected_uid=AUTHORIZED_UID,
            expected_gid=AUTHORIZED_GID,
            expected_pic_root_identity=REVIEWED_PIC_ROOT_IDENTITY,
            expected_recovery_identity=REVIEWED_ACCEPTANCE_ROOT_IDENTITY,
            expected_recovery_xattrs=REVIEWED_RECOVERY_XATTRS,
            recover_exact_empty_inherited_setgid_root=(
                args.recover_exact_empty_inherited_setgid_root
            ),
            reconcile_exact_empty_normalized_root=(
                args.reconcile_exact_empty_normalized_root
            ),
            require_empty_recovery_receipt_namespace=(
                args.recover_exact_empty_inherited_setgid_root
                or args.reconcile_exact_empty_normalized_root
            ),
            **common,
        )
    print(canonical_json_bytes(result).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
