#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Atomically promote one reviewed policy into the mirrored active-policy anchor."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import re
import stat
from typing import Callable, Iterator
import uuid

from control_plane_common import AUTHORIZED_ACCOUNT, AUTHORIZED_PIC_ROOT
from control_plane_common import AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import ACTIVE_PROMOTION_ROLLBACK_ANCHOR_PREFIXES
from control_plane_common import active_promotion_path, active_promotion_transaction_path
from control_plane_common import atomic_write_bytes_at
from control_plane_common import atomic_write_json_at
from control_plane_common import canonical_policy_path, durable_mkdir_parents
from control_plane_common import open_directory_below, require_same_directory
from control_plane_common import project_home_ledger_root
from control_plane_common import read_json_bytes
from control_plane_common import read_stable_regular_file, sha256_bytes
from control_plane_common import require_no_symlink_components_below
from control_plane_common import require_policy_predecessor_snapshot_for_promotion
from control_plane_common import require_storage_policy_unlock_snapshot
from control_plane_common import stable_serialization_anchor
from control_plane_common import validate_storage_policy, verify_installed_control_plane
from control_plane_common import verify_historical_installed_control_plane
from ledger import _path_exists, _pinned_parent_directories
from ledger import latest_reservations, validate_mirrored_state, validate_receipts
from ledger import require_no_incomplete_manual_accounting_marker
from operator_attestation import validate_sealed_operator_attestation
from revalidate_clean_candidate import revalidate_clean_candidate


SCRIPT_DIR = Path(__file__).absolute().parent
Q011_JOB_SCRIPT_SHA256 = (
    "3048493d376dfa7954595e586c7cebe6740460a7455110d0fc4020bedf697999"
)
Q011_INPUT_DECK_SHA256 = (
    "0b1cbd62d54027ec81a5f4f5c88d5ee56b86b8cc0cb018c3fbebfb37a11be7b1"
)
Q011_LAUNCH_CONTRACT_SHA256 = {
    "2413a91247d32fb6d93d4903bddab65ed518badc1be6c514ad29c03ca8eafb7b",
    "f9f615bfaa4cc18479dcd02ed4f688f3723fc3dd30a50e754e9f69e10616a8ea",
    "8ab4528b55bdf71e4a13047971aa5ccf8da35c28896ba0bf875a7dca368dfd2f",
    "d84b9217c8810f33ffda995f9611b44f7eeb519bb2e7667c4691ab5833895dd4",
}
PROMOTION_TRANSACTION_RECORD_TYPE = "frontier_pic_active_policy_promotion_transaction"
PROMOTION_TRANSACTION_SCHEMA_VERSION = 2


def _lowercase_sha256(value: object, *, label: str) -> str:
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError(f"{label} is not a lowercase SHA-256 digest")
    return value


def _canonical_uuid(value: object, *, label: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{label} is not a canonical UUID")
    try:
        parsed = uuid.UUID(value)
    except ValueError as error:
        raise ValueError(f"{label} is not a canonical UUID") from error
    if str(parsed) != value:
        raise ValueError(f"{label} is not a canonical UUID")
    return value


def _read_regular_file_at(
    parent_descriptor: int, name: str, *, label: str
) -> bytes:
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=parent_descriptor,
    )
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_mode & 0o222:
            raise ValueError(f"{label} is not one read-only regular file")
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        if (
            (before.st_dev, before.st_ino, before.st_mode, before.st_size)
            != (after.st_dev, after.st_ino, after.st_mode, after.st_size)
            or len(payload) != after.st_size
        ):
            raise ValueError(f"{label} changed while reading")
        return bytes(payload)
    finally:
        os.close(descriptor)


def _optional_regular_file_at(
    parent_descriptor: int, name: str, *, label: str
) -> bytes | None:
    try:
        return _read_regular_file_at(parent_descriptor, name, label=label)
    except FileNotFoundError:
        return None


def _unlink_if_exists_at(parent_descriptor: int, name: str) -> None:
    try:
        os.unlink(name, dir_fd=parent_descriptor)
    except FileNotFoundError:
        return


@contextmanager
def _locked_mirror_policy_parent(
    mirror_policy_parent: Path, *, authorized_project_home_root: Path
) -> Iterator[int]:
    descriptor = open_directory_below(
        mirror_policy_parent, root=authorized_project_home_root
    )
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        try:
            require_same_directory(
                mirror_policy_parent,
                descriptor,
                root=authorized_project_home_root,
            )
            yield descriptor
            require_same_directory(
                mirror_policy_parent,
                descriptor,
                root=authorized_project_home_root,
            )
        finally:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
    finally:
        os.close(descriptor)


def _validate_promotion_transaction_marker(
    marker: dict[str, object],
) -> dict[str, object]:
    if set(marker) != {
        "schema_version",
        "record_type",
        "transaction_id",
        "state",
        "predecessor_state",
        "anchors",
    } or type(marker.get("schema_version")) is not int or marker.get(
        "schema_version"
    ) != PROMOTION_TRANSACTION_SCHEMA_VERSION or marker.get(
        "record_type"
    ) != PROMOTION_TRANSACTION_RECORD_TYPE:
        raise ValueError("Active-policy promotion transaction marker is malformed")
    transaction_id = _canonical_uuid(
        marker.get("transaction_id"),
        label="Active-policy promotion transaction ID",
    )
    anchors = marker.get("anchors")
    expected = [
        ("orion", "storage_policy.json"),
        ("orion", "active_promotion.json"),
        ("project_home", "storage_policy.json"),
        ("project_home", "active_promotion.json"),
    ]
    if not isinstance(anchors, list) or len(anchors) != len(expected):
        raise ValueError("Active-policy promotion transaction anchors are malformed")
    state = marker.get("state")
    predecessor_state = marker.get("predecessor_state")
    if state not in {"prepared", "committed"} or predecessor_state not in {
        "absent",
        "complete",
    }:
        raise ValueError("Active-policy promotion transaction state is malformed")
    predecessor_digests: list[object] = []
    successor_digests: list[str] = []
    for anchor, (root_role, name) in zip(anchors, expected):
        rollback_name = f".{name}.transaction-rollback-{transaction_id}"
        if (
            not isinstance(anchor, dict)
            or set(anchor)
            != {
                "root_role",
                "name",
                "rollback_name",
                "predecessor_sha256",
                "successor_sha256",
            }
            or anchor.get("root_role") != root_role
            or anchor.get("name") != name
            or anchor.get("rollback_name") != rollback_name
            or (
                anchor.get("predecessor_sha256") is not None
                and (
                    not isinstance(anchor.get("predecessor_sha256"), str)
                    or len(str(anchor["predecessor_sha256"])) != 64
                    or any(
                        character not in "0123456789abcdef"
                        for character in str(anchor["predecessor_sha256"])
                    )
                )
            )
            or not isinstance(anchor.get("successor_sha256"), str)
            or len(str(anchor["successor_sha256"])) != 64
            or any(
                character not in "0123456789abcdef"
                for character in str(anchor["successor_sha256"])
            )
        ):
            raise ValueError("Active-policy promotion transaction anchor is malformed")
        predecessor_digests.append(anchor["predecessor_sha256"])
        successor_digests.append(str(anchor["successor_sha256"]))
    if (
        predecessor_state == "absent"
        and any(digest is not None for digest in predecessor_digests)
    ) or (
        predecessor_state == "complete"
        and any(digest is None for digest in predecessor_digests)
    ):
        raise ValueError("Active-policy promotion predecessor state is malformed")
    if (
        successor_digests[0] != successor_digests[2]
        or successor_digests[1] != successor_digests[3]
        or (
            predecessor_state == "complete"
            and (
                predecessor_digests[0] != predecessor_digests[2]
                or predecessor_digests[1] != predecessor_digests[3]
            )
        )
    ):
        raise ValueError("Active-policy promotion transaction mirror digests differ")
    return marker


def _read_promotion_transaction_marker(
    policy_descriptor: int,
    mirror_policy_descriptor: int,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> dict[str, object] | None:
    orion_name = active_promotion_transaction_path(authorized_pic_root).name
    mirror_name = active_promotion_transaction_path(
        authorized_project_home_root
    ).name
    orion_bytes = _optional_regular_file_at(
        policy_descriptor,
        orion_name,
        label="Orion active-policy promotion transaction marker",
    )
    mirror_bytes = _optional_regular_file_at(
        mirror_policy_descriptor,
        mirror_name,
        label="Project Home active-policy promotion transaction marker",
    )
    if orion_bytes is None and mirror_bytes is None:
        return None
    markers = [
        _validate_promotion_transaction_marker(
            read_json_bytes(payload, label="active-policy promotion transaction marker")
        )
        for payload in [orion_bytes, mirror_bytes]
        if payload is not None
    ]
    assert markers
    expected = {**markers[0], "state": None}
    if any({**marker, "state": None} != expected for marker in markers[1:]):
        raise ValueError("Mirrored active-policy promotion transaction markers differ")
    state = (
        "committed"
        if any(marker["state"] == "committed" for marker in markers)
        else "prepared"
    )
    return {**markers[0], "state": state}


def _require_visible_prepared_promotion_transaction_marker(
    marker: dict[str, object],
    policy_descriptor: int,
    mirror_policy_descriptor: int,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    allow_absent: bool = False,
) -> None:
    if marker.get("state") != "prepared":
        raise ValueError("Committed active-policy promotion transaction cannot roll back")
    current = _read_promotion_transaction_marker(
        policy_descriptor,
        mirror_policy_descriptor,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    if current is None and allow_absent:
        return
    if current != marker:
        raise ValueError(
            "Prepared active-policy promotion transaction marker changed before rollback"
        )


def _promotion_transaction_marker_cleanup_order(
    marker: dict[str, object],
    policy_descriptor: int,
    mirror_policy_descriptor: int,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    required_state: str | None = None,
) -> list[tuple[int, str]]:
    expected = {**marker, "state": None}
    markers: list[tuple[int, int, str]] = []
    for descriptor, marker_name in [
        (
            mirror_policy_descriptor,
            active_promotion_transaction_path(authorized_project_home_root).name,
        ),
        (
            policy_descriptor,
            active_promotion_transaction_path(authorized_pic_root).name,
        ),
    ]:
        payload = _optional_regular_file_at(
            descriptor,
            marker_name,
            label="Active-policy promotion transaction marker during cleanup",
        )
        if payload is None:
            continue
        current = _validate_promotion_transaction_marker(
            read_json_bytes(
                payload,
                label="active-policy promotion transaction marker during cleanup",
            )
        )
        if {**current, "state": None} != expected:
            raise ValueError(
                "Active-policy promotion transaction marker changed during cleanup"
            )
        if required_state is not None and current["state"] != required_state:
            raise ValueError(
                "Active-policy promotion transaction marker state changed during cleanup"
            )
        markers.append(
            (
                0 if current["state"] == "committed" else 1,
                descriptor,
                marker_name,
            )
        )
    return [(descriptor, marker_name) for _, descriptor, marker_name in sorted(markers)]


def _promotion_transaction_publication_prefix_length(
    marker: dict[str, object],
    descriptors: dict[str, int],
) -> int:
    anchors = marker["anchors"]
    assert isinstance(anchors, list)
    current_digests: list[str | None] = []
    predecessor_digests: list[str | None] = []
    successor_digests: list[str] = []
    for anchor in anchors:
        assert isinstance(anchor, dict)
        active_bytes = _optional_regular_file_at(
            descriptors[str(anchor["root_role"])],
            str(anchor["name"]),
            label="Prepared active-policy promotion current anchor",
        )
        current_digests.append(
            sha256_bytes(active_bytes) if active_bytes is not None else None
        )
        predecessor_sha256 = anchor["predecessor_sha256"]
        assert predecessor_sha256 is None or isinstance(predecessor_sha256, str)
        predecessor_digests.append(predecessor_sha256)
        successor_digests.append(str(anchor["successor_sha256"]))
    publication_order = [0, 2, 1, 3]
    matching_prefix_lengths = [
        publication_prefix_length
        for publication_prefix_length in range(len(publication_order) + 1)
        if all(
            current_digests[index]
            == (
                successor_digests[index]
                if index in publication_order[:publication_prefix_length]
                else predecessor_digests[index]
            )
            for index in range(len(anchors))
        )
    ]
    if not matching_prefix_lengths:
        raise ValueError(
            "Prepared active-policy promotion current anchors are not one reachable "
            "publication-prefix state"
        )
    return max(matching_prefix_lengths)


def _rollback_promotion_transaction(
    marker: dict[str, object],
    policy_descriptor: int,
    mirror_policy_descriptor: int,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    post_publish_check: Callable[[], None],
    validate_predecessor: Callable[
        [bytes | None, bytes | None],
        tuple[dict[str, object], dict[str, str]] | None,
    ],
    allow_absent_predecessor: bool = False,
    allow_absent_prepared_marker: bool = False,
) -> None:
    if marker.get("state") != "prepared":
        raise ValueError("Committed active-policy promotion transaction cannot roll back")
    descriptors = {
        "orion": policy_descriptor,
        "project_home": mirror_policy_descriptor,
    }
    anchors = marker["anchors"]
    assert isinstance(anchors, list)
    predecessor_state = marker["predecessor_state"]
    if predecessor_state == "absent" and not allow_absent_predecessor:
        raise ValueError(
            "Absent-predecessor promotion transaction requires explicit manual recovery"
        )
    post_publish_check()
    predecessor_digests: list[str | None] = []
    for anchor in anchors:
        assert isinstance(anchor, dict)
        predecessor_sha256 = anchor["predecessor_sha256"]
        assert predecessor_sha256 is None or isinstance(predecessor_sha256, str)
        predecessor_digests.append(predecessor_sha256)
    publication_order = [0, 2, 1, 3]
    publication_prefix_length = _promotion_transaction_publication_prefix_length(
        marker, descriptors
    )
    if publication_prefix_length == len(publication_order):
        raise ValueError(
            "Complete active-policy promotion successor cannot roll back"
        )
    rollback_payloads: list[bytes | None] = []
    for anchor in anchors:
        assert isinstance(anchor, dict)
        descriptor = descriptors[str(anchor["root_role"])]
        rollback_name = str(anchor["rollback_name"])
        predecessor_sha256 = anchor["predecessor_sha256"]
        if predecessor_sha256 is None:
            if (
                _optional_regular_file_at(
                    descriptor,
                    rollback_name,
                    label=f"Absent-predecessor rollback anchor {rollback_name}",
                )
                is not None
            ):
                raise ValueError(
                    "Absent-predecessor promotion transaction has a rollback anchor"
                )
            rollback_payloads.append(None)
            continue
        rollback_bytes = _read_regular_file_at(
            descriptor,
            rollback_name,
            label=f"Active-policy promotion rollback anchor {rollback_name}",
        )
        if sha256_bytes(rollback_bytes) != predecessor_sha256:
            raise ValueError("Active-policy promotion rollback anchor digest differs")
        rollback_payloads.append(rollback_bytes)
    expected_predecessor_snapshot = (
        None
        if predecessor_state == "absent"
        else {
            "active_policy_sha256": str(predecessor_digests[0]),
            "active_promotion_sha256": str(predecessor_digests[1]),
        }
    )
    if predecessor_state == "complete":
        policy_bytes = rollback_payloads[0]
        promotion_bytes = rollback_payloads[1]
        mirror_policy_bytes = rollback_payloads[2]
        mirror_promotion_bytes = rollback_payloads[3]
        assert isinstance(policy_bytes, bytes)
        assert isinstance(promotion_bytes, bytes)
        assert isinstance(mirror_policy_bytes, bytes)
        assert isinstance(mirror_promotion_bytes, bytes)
        if (
            policy_bytes != mirror_policy_bytes
            or promotion_bytes != mirror_promotion_bytes
        ):
            raise ValueError(
                "Prepared active-policy promotion rollback mirrors differ"
            )
        validated_predecessor = validate_predecessor(
            policy_bytes,
            promotion_bytes,
        )
        if (
            validated_predecessor is None
            or validated_predecessor[1] != expected_predecessor_snapshot
        ):
            raise ValueError(
                "Prepared active-policy promotion semantic predecessor snapshot differs"
            )
    post_publish_check()
    _require_visible_prepared_promotion_transaction_marker(
        marker,
        policy_descriptor,
        mirror_policy_descriptor,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        allow_absent=allow_absent_prepared_marker,
    )
    for index in reversed(publication_order[:publication_prefix_length]):
        _require_visible_prepared_promotion_transaction_marker(
            marker,
            policy_descriptor,
            mirror_policy_descriptor,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
            allow_absent=allow_absent_prepared_marker,
        )
        anchor = anchors[index]
        rollback_bytes = rollback_payloads[index]
        assert isinstance(anchor, dict)
        descriptor = descriptors[str(anchor["root_role"])]
        name = str(anchor["name"])
        if rollback_bytes is None:
            _unlink_if_exists_at(descriptor, name)
            os.fsync(descriptor)
            continue
        atomic_write_bytes_at(
            descriptor,
            name,
            rollback_bytes,
            post_publish_check=post_publish_check,
        )
        _require_visible_prepared_promotion_transaction_marker(
            marker,
            policy_descriptor,
            mirror_policy_descriptor,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
            allow_absent=allow_absent_prepared_marker,
        )
    post_publish_check()
    restored_predecessor = validate_predecessor(None, None)
    restored_predecessor_snapshot = (
        restored_predecessor[1] if restored_predecessor is not None else None
    )
    if restored_predecessor_snapshot != expected_predecessor_snapshot:
        raise ValueError(
            "Prepared active-policy promotion restored predecessor snapshot differs"
        )
    post_publish_check()
    marker_cleanup_order = _promotion_transaction_marker_cleanup_order(
        marker,
        policy_descriptor,
        mirror_policy_descriptor,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        required_state="prepared",
    )
    for descriptor, marker_name in marker_cleanup_order:
        _unlink_if_exists_at(descriptor, marker_name)
        os.fsync(descriptor)
    for anchor in anchors:
        assert isinstance(anchor, dict)
        descriptor = descriptors[str(anchor["root_role"])]
        _unlink_if_exists_at(descriptor, str(anchor["rollback_name"]))
        os.fsync(descriptor)


def _require_promotion_transaction_successor(
    marker: dict[str, object],
    descriptors: dict[str, int],
) -> None:
    anchors = marker["anchors"]
    assert isinstance(anchors, list)
    for anchor in anchors:
        assert isinstance(anchor, dict)
        active_bytes = _read_regular_file_at(
            descriptors[str(anchor["root_role"])],
            str(anchor["name"]),
            label="Committed active-policy promotion successor anchor",
        )
        if sha256_bytes(active_bytes) != anchor["successor_sha256"]:
            raise ValueError(
                "Committed active-policy promotion successor anchor digest differs"
        )


def _require_valid_complete_promotion_transaction_successor(
    marker: dict[str, object],
    descriptors: dict[str, int],
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    authorized_account: str,
    post_publish_check: Callable[[], None],
) -> tuple[dict[str, object], str]:
    _require_promotion_transaction_successor(marker, descriptors)
    promotion_bytes = _read_regular_file_at(
        descriptors["orion"],
        active_promotion_path(authorized_pic_root).name,
        label="Committed active-policy promotion successor record",
    )
    promotion = read_json_bytes(
        promotion_bytes,
        label="committed active-policy promotion successor record",
    )
    version = promotion.get("control_plane_version")
    if (
        not isinstance(version, str)
        or len(version) != 64
        or any(character not in "0123456789abcdef" for character in version)
    ):
        raise ValueError(
            "Committed active-policy promotion successor controller is malformed"
        )
    inventory: dict[str, object] | None = None
    for root in [authorized_pic_root, authorized_project_home_root]:
        current = verify_historical_installed_control_plane(
            Path(os.path.abspath(root)) / "control_plane" / version,
            authorized_pic_root=root,
        )
        if inventory is not None and current != inventory:
            raise ValueError(
                "Committed active-policy promotion successor controller mirrors differ"
            )
        inventory = current
    post_publish_check()
    policy, snapshot = require_storage_policy_unlock_snapshot(
        control_plane_version=version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        allow_pending_genesis=True,
        allow_active_promotion_transaction=True,
    )
    anchors = marker["anchors"]
    assert isinstance(anchors, list)
    policy_anchor = anchors[0]
    promotion_anchor = anchors[1]
    assert isinstance(policy_anchor, dict)
    assert isinstance(promotion_anchor, dict)
    expected_snapshot = {
        "active_policy_sha256": str(policy_anchor["successor_sha256"]),
        "active_promotion_sha256": str(promotion_anchor["successor_sha256"]),
    }
    if snapshot != expected_snapshot:
        raise ValueError(
            "Committed active-policy promotion semantic successor snapshot differs"
        )
    post_publish_check()
    return policy, version


def _finalize_complete_promotion_transaction_successor(
    marker: dict[str, object],
    policy_descriptor: int,
    mirror_policy_descriptor: int,
    descriptors: dict[str, int],
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    authorized_account: str,
    post_publish_check: Callable[[], None],
    validate_active_successor: Callable[[dict[str, object], str], None],
) -> None:
    anchors = marker["anchors"]
    assert isinstance(anchors, list)
    policy, version = _require_valid_complete_promotion_transaction_successor(
        marker,
        descriptors,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        post_publish_check=post_publish_check,
    )
    validate_active_successor(policy, version)
    _require_promotion_transaction_successor(marker, descriptors)
    if (
        _read_promotion_transaction_marker(
            policy_descriptor,
            mirror_policy_descriptor,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        != marker
    ):
        raise ValueError(
            "Complete active-policy promotion transaction marker changed during "
            "recovery"
        )
    post_publish_check()
    _remove_promotion_transaction_rollback_anchors(
        anchors,
        descriptors,
    )
    for descriptor, marker_name in [
        (
            mirror_policy_descriptor,
            active_promotion_transaction_path(authorized_project_home_root).name,
        ),
        (
            policy_descriptor,
            active_promotion_transaction_path(authorized_pic_root).name,
        ),
    ]:
        try:
            _unlink_if_exists_at(descriptor, marker_name)
            os.fsync(descriptor)
        except OSError as error:
            if (
                _read_promotion_transaction_marker(
                    policy_descriptor,
                    mirror_policy_descriptor,
                    authorized_pic_root=authorized_pic_root,
                    authorized_project_home_root=authorized_project_home_root,
                )
                is not None
            ):
                raise
            print(
                "warning: recovered complete active-policy promotion with "
                f"marker cleanup error: {error}",
                file=_sys.stderr,
            )
            break


def _require_authorized_clean_candidate_freeze_revalidation(
    policy: dict[str, object],
    *,
    control_plane_version: str,
    control_plane_dir: Path,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> None:
    science_freeze = policy.get("science_submission_freeze")
    if (
        not isinstance(science_freeze, dict)
        or science_freeze.get("status") != "authorized"
    ):
        return
    receipt_control_plane_version = science_freeze.get(
        "build_profile_control_plane_version"
    )
    if not isinstance(receipt_control_plane_version, str):
        raise ValueError(
            "Authorized clean-candidate freeze build controller is malformed"
        )
    revalidation = revalidate_clean_candidate(
        Path(str(science_freeze["manifest_path"])),
        expected_manifest_sha256=str(science_freeze["manifest_sha256"]),
        expected_receipt_control_plane_version=receipt_control_plane_version,
        control_plane_dir=control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    build = revalidation.get("build")
    if (
        revalidation.get("status") != "passed"
        or revalidation.get("current_control_plane_version")
        != control_plane_version
        or not isinstance(build, dict)
        or build.get("receipt_control_plane_version")
        != receipt_control_plane_version
    ):
        raise ValueError(
            "Authorized clean-candidate freeze revalidation is not bound to the "
            "active and build-receipt control planes"
        )


def _recover_promotion_transaction(
    policy_descriptor: int,
    mirror_policy_descriptor: int,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    authorized_account: str,
    post_publish_check: Callable[[], None],
    validate_predecessor: Callable[
        [bytes | None, bytes | None],
        tuple[dict[str, object], dict[str, str]] | None,
    ],
    validate_active_successor: Callable[[dict[str, object], str], None],
) -> None:
    marker = _read_promotion_transaction_marker(
        policy_descriptor,
        mirror_policy_descriptor,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    descriptors = {
        "orion": policy_descriptor,
        "project_home": mirror_policy_descriptor,
    }
    reserved_anchor_names = {
        root_role: {
            name
            for name in os.listdir(descriptor)
            if any(
                name.startswith(prefix)
                for prefix in ACTIVE_PROMOTION_ROLLBACK_ANCHOR_PREFIXES
            )
        }
        for root_role, descriptor in descriptors.items()
    }
    if marker is None:
        if any(reserved_anchor_names.values()):
            raise ValueError(
                "Markerless active-policy promotion rollback anchors require "
                "reviewed manual recovery"
            )
        return
    expected_anchor_names: dict[str, set[str]] = {
        "orion": set(),
        "project_home": set(),
    }
    anchors = marker["anchors"]
    assert isinstance(anchors, list)
    for anchor in anchors:
        assert isinstance(anchor, dict)
        if anchor["predecessor_sha256"] is not None:
            expected_anchor_names[str(anchor["root_role"])].add(
                str(anchor["rollback_name"])
            )
    if any(
        reserved_anchor_names[root_role] - expected_anchor_names[root_role]
        for root_role in descriptors
    ):
        raise ValueError(
            "Unexpected active-policy promotion rollback anchors require reviewed "
            "manual recovery"
        )
    if marker is not None:
        if (
            marker["state"] == "committed"
            or _promotion_transaction_publication_prefix_length(marker, descriptors)
            == 4
        ):
            _finalize_complete_promotion_transaction_successor(
                marker,
                policy_descriptor,
                mirror_policy_descriptor,
                descriptors,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
                authorized_account=authorized_account,
                post_publish_check=post_publish_check,
                validate_active_successor=validate_active_successor,
            )
        else:
            _rollback_promotion_transaction(
                marker,
                policy_descriptor,
                mirror_policy_descriptor,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
                post_publish_check=post_publish_check,
                validate_predecessor=validate_predecessor,
            )


def _remove_promotion_transaction_rollback_anchors(
    anchors: list[object],
    descriptors: dict[str, int],
) -> None:
    for anchor in reversed(anchors):
        assert isinstance(anchor, dict)
        descriptor = descriptors[str(anchor["root_role"])]
        _unlink_if_exists_at(descriptor, str(anchor["rollback_name"]))
        os.fsync(descriptor)


@contextmanager
def _active_policy_transaction(
    policy_descriptor: int,
    mirror_policy_descriptor: int,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    post_publish_check: Callable[[], None],
    expected_predecessor_snapshot: dict[str, str] | None,
    expected_successor_snapshot: dict[str, str],
    validate_predecessor: Callable[
        [bytes | None, bytes | None],
        tuple[dict[str, object], dict[str, str]] | None,
    ],
    validate_successor: Callable[[], None],
) -> Iterator[None]:
    transaction_id = str(uuid.uuid4())
    descriptors = {
        "orion": policy_descriptor,
        "project_home": mirror_policy_descriptor,
    }
    marker: dict[str, object] = {
        "schema_version": PROMOTION_TRANSACTION_SCHEMA_VERSION,
        "record_type": PROMOTION_TRANSACTION_RECORD_TYPE,
        "transaction_id": transaction_id,
        "state": "prepared",
        "predecessor_state": (
            "absent" if expected_predecessor_snapshot is None else "complete"
        ),
        "anchors": [],
    }
    anchors = marker["anchors"]
    assert isinstance(anchors, list)
    try:
        for root_role, name in [
            ("orion", "storage_policy.json"),
            ("orion", "active_promotion.json"),
            ("project_home", "storage_policy.json"),
            ("project_home", "active_promotion.json"),
        ]:
            descriptor = descriptors[root_role]
            rollback_name = f".{name}.transaction-rollback-{transaction_id}"
            predecessor = _optional_regular_file_at(
                descriptor,
                name,
                label=f"Active-policy predecessor anchor {root_role}/{name}",
            )
            expected_predecessor_sha256 = (
                None
                if expected_predecessor_snapshot is None
                else expected_predecessor_snapshot[
                    (
                        "active_policy_sha256"
                        if name == "storage_policy.json"
                        else "active_promotion_sha256"
                    )
                ]
            )
            predecessor_sha256 = (
                sha256_bytes(predecessor) if predecessor is not None else None
            )
            if predecessor_sha256 != expected_predecessor_sha256:
                raise ValueError(
                    "Active-policy predecessor changed before transaction setup"
                )
            successor_sha256 = expected_successor_snapshot[
                (
                    "active_policy_sha256"
                    if name == "storage_policy.json"
                    else "active_promotion_sha256"
                )
            ]
            anchors.append(
                {
                    "root_role": root_role,
                    "name": name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": predecessor_sha256,
                    "successor_sha256": successor_sha256,
                }
            )
            if predecessor is not None:
                os.link(
                    name,
                    rollback_name,
                    src_dir_fd=descriptor,
                    dst_dir_fd=descriptor,
                    follow_symlinks=False,
                )
            os.fsync(descriptor)
        for anchor in anchors:
            assert isinstance(anchor, dict)
            if anchor["predecessor_sha256"] is None:
                continue
            descriptor = descriptors[str(anchor["root_role"])]
            rollback_bytes = _read_regular_file_at(
                descriptor,
                str(anchor["rollback_name"]),
                label="Active-policy promotion setup rollback anchor",
            )
            if sha256_bytes(rollback_bytes) != anchor["predecessor_sha256"]:
                raise ValueError(
                    "Active-policy predecessor changed during transaction setup"
                )
    except BaseException:
        try:
            _remove_promotion_transaction_rollback_anchors(anchors, descriptors)
        except BaseException as cleanup_error:
            raise RuntimeError(
                "Failed to clean up active-policy promotion transaction setup"
            ) from cleanup_error
        raise
    marker_names = [
        (
            policy_descriptor,
            active_promotion_transaction_path(authorized_pic_root).name,
        ),
        (
            mirror_policy_descriptor,
            active_promotion_transaction_path(authorized_project_home_root).name,
        ),
    ]
    commit_started = False
    committed = False
    committed_marker = {**marker, "state": "committed"}
    try:
        for descriptor, marker_name in marker_names:
            atomic_write_json_at(
                descriptor,
                marker_name,
                marker,
                mode=0o400,
                replace=False,
                post_publish_check=post_publish_check,
            )
        yield
        if _read_promotion_transaction_marker(
            policy_descriptor,
            mirror_policy_descriptor,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        ) != marker:
            raise ValueError("Active-policy promotion transaction marker changed")
        _require_promotion_transaction_successor(marker, descriptors)
        validate_successor()
        _require_promotion_transaction_successor(marker, descriptors)
        post_publish_check()
        commit_started = True
        for descriptor, marker_name in marker_names:
            atomic_write_json_at(
                descriptor,
                marker_name,
                committed_marker,
                mode=0o400,
                replace=True,
                post_publish_check=post_publish_check,
            )
        if (
            _read_promotion_transaction_marker(
                policy_descriptor,
                mirror_policy_descriptor,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )
            != committed_marker
        ):
            raise ValueError(
                "Committed active-policy promotion transaction marker changed"
            )
        _require_promotion_transaction_successor(marker, descriptors)
        validate_successor()
        _require_promotion_transaction_successor(marker, descriptors)
        if (
            _read_promotion_transaction_marker(
                policy_descriptor,
                mirror_policy_descriptor,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )
            != committed_marker
        ):
            raise ValueError(
                "Committed active-policy promotion transaction marker changed"
            )
        post_publish_check()
        committed = True
    except BaseException as error:
        current_marker: dict[str, object] | None = None
        try:
            current_marker = _read_promotion_transaction_marker(
                policy_descriptor,
                mirror_policy_descriptor,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )
            publication_prefix_length = (
                _promotion_transaction_publication_prefix_length(marker, descriptors)
            )
        except BaseException as recovery_state_error:
            raise RuntimeError(
                "Active-policy promotion commit state requires locked recovery"
                if commit_started
                else "Active-policy promotion state requires locked recovery"
            ) from recovery_state_error
        same_marker_generation = current_marker is not None and {
            **current_marker,
            "state": None,
        } == {**marker, "state": None}
        if commit_started and not same_marker_generation:
            raise RuntimeError(
                "Active-policy promotion commit state requires locked recovery"
            ) from error
        if publication_prefix_length == 4:
            if not same_marker_generation:
                raise RuntimeError(
                    "Complete active-policy promotion requires locked recovery"
                ) from error
            try:
                _require_promotion_transaction_successor(marker, descriptors)
                validate_successor()
                _require_promotion_transaction_successor(marker, descriptors)
                if (
                    _read_promotion_transaction_marker(
                        policy_descriptor,
                        mirror_policy_descriptor,
                        authorized_pic_root=authorized_pic_root,
                        authorized_project_home_root=authorized_project_home_root,
                    )
                    != current_marker
                ):
                    raise ValueError(
                        "Complete active-policy promotion transaction marker changed"
                    )
                post_publish_check()
                committed = True
            except BaseException as validation_error:
                raise RuntimeError(
                    "Committed active-policy promotion requires locked recovery"
                ) from validation_error
        elif commit_started or (
            current_marker is not None
            and (
                not same_marker_generation
                or current_marker["state"] == "committed"
            )
        ) or (current_marker is None and publication_prefix_length != 0):
            raise RuntimeError(
                "Active-policy promotion state requires locked recovery"
            ) from error
        if committed:
            print(
                f"warning: active-policy promotion committed before cleanup error: {error}",
                file=_sys.stderr,
            )
        else:
            try:
                _rollback_promotion_transaction(
                    marker,
                    policy_descriptor,
                    mirror_policy_descriptor,
                    authorized_pic_root=authorized_pic_root,
                    authorized_project_home_root=authorized_project_home_root,
                    post_publish_check=post_publish_check,
                    validate_predecessor=validate_predecessor,
                    allow_absent_predecessor=True,
                    allow_absent_prepared_marker=(
                        current_marker is None and publication_prefix_length == 0
                    ),
                )
            except BaseException as rollback_error:
                raise RuntimeError(
                    "Failed to roll back active-policy promotion transaction"
                ) from rollback_error
            raise
    if committed:
        try:
            _remove_promotion_transaction_rollback_anchors(anchors, descriptors)
        except OSError as error:
            print(
                f"warning: active-policy promotion committed with retained rollback "
                f"anchor and marker: {error}",
                file=_sys.stderr,
            )
        else:
            try:
                for descriptor, marker_name in reversed(marker_names):
                    _unlink_if_exists_at(descriptor, marker_name)
                    os.fsync(descriptor)
            except OSError as error:
                print(
                    f"warning: active-policy promotion committed with retained marker: "
                    f"{error}",
                    file=_sys.stderr,
                )


def _requires_q011_launch_prohibited_baseline(policy: dict[str, object]) -> bool:
    slices = policy.get("registered_science_slices")
    return isinstance(slices, list) and any(
        isinstance(record, dict)
        and (
            str(record.get("authorization_id", "")).startswith(
                "q011-section54-pressure-"
            )
            or str(record.get("campaign", "")).startswith(
                "q011_section54_pressure_"
            )
            or str(record.get("test_id", "")).startswith(
                "pic_parallel_shock_section54_pressure_"
            )
            or record.get("job_script_sha256") == Q011_JOB_SCRIPT_SHA256
            or record.get("input_deck_sha256") == Q011_INPUT_DECK_SHA256
            or record.get("launch_contract_sha256") in Q011_LAUNCH_CONTRACT_SHA256
        )
        for record in slices
    )


def _require_q011_launch_prohibited_baseline(
    policy: dict[str, object],
    *,
    control_plane_version: str,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    authorized_account: str,
) -> None:
    if not _requires_q011_launch_prohibited_baseline(policy):
        return
    active_policy, _ = require_storage_policy_unlock_snapshot(
        control_plane_version=control_plane_version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        allow_pending_genesis=True,
    )
    if active_policy.get("registered_science_slices") != []:
        raise ValueError(
            "Q011 pressure-slice promotion requires one active launch-prohibited "
            "same-controller baseline"
        )


def _require_no_outstanding_submissions(
    policy: dict[str, object],
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> None:
    ledger_root = Path(os.path.abspath(authorized_pic_root)) / "ledger"
    project_ledger_root = (
        project_home_ledger_root(authorized_project_home_root) / "ledger"
    )
    marker = ledger_root / "pending_submission.json"
    ledger = ledger_root / "node_hours.jsonl"
    receipts = ledger_root / "mirror_receipts.jsonl"
    mirror = project_ledger_root / "node_hours.jsonl"
    paths = [ledger, receipts, mirror]
    resolved_mirror = mirror.resolve()
    with _pinned_parent_directories([marker, *paths, resolved_mirror]):
        if _path_exists(marker):
            raise ValueError("Policy promotion is blocked by a pending scheduler submission")
        existing = [_path_exists(path) for path in paths]
        if not any(existing):
            storage = policy["olcf_side_storage"]
            if (
                not isinstance(storage, dict)
                or storage.get("ledger_genesis_allowed") is not True
                or storage.get("ledger_genesis") is not None
            ):
                raise ValueError(
                    "Fresh-root policy promotion requires pending ledger genesis"
                )
            return
        if not all(existing):
            raise ValueError("Policy promotion is blocked by incomplete mirrored ledger state")
        validated_mirror = mirror
        try:
            records = validate_mirrored_state(ledger, receipts, mirror)
        except ValueError as lexical_error:
            if resolved_mirror == mirror:
                raise
            try:
                records = validate_mirrored_state(ledger, receipts, resolved_mirror)
                validated_mirror = resolved_mirror
            except ValueError:
                raise lexical_error
        outstanding = [
            record for record in latest_reservations(records).values()
            if record.get("state") in {"reserved", "submitted"}
        ]
        if outstanding:
            raise ValueError("Policy promotion is blocked by an outstanding reservation")
        receipts_records = validate_receipts(
            receipts,
            records,
            mirror_jsonl=validated_mirror,
            mirror_transport="filesystem_copy",
        )
        storage = policy["olcf_side_storage"]
        expected_genesis = {
            "status": "initialized",
            "timestamp": records[0]["timestamp"],
            "control_plane_version": records[0]["control_plane_version"],
            "event_sha256": records[0]["event_sha256"],
            "mirror_transport": receipts_records[0]["mirror_transport"],
            "mirror_ack_sha256": receipts_records[0]["mirror_ack_sha256"],
        }
        if (
            not isinstance(storage, dict)
            or storage.get("ledger_genesis_allowed") is not False
            or storage.get("ledger_genesis") != expected_genesis
        ):
            raise ValueError(
                "Policy promotion ledger genesis does not match retained ledger bindings"
            )


def _require_same_directory(path: Path, descriptor: int) -> None:
    lexical_descriptor = os.open(
        path,
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
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
def _promotion_lock_within_serialization_anchor(
    authorized_pic_root: Path,
) -> Iterator[int]:
    lexical_root = Path(os.path.abspath(authorized_pic_root))
    policy_parent = lexical_root / "policy"
    path = lexical_root / ".promotion.lock"
    durable_mkdir_parents(policy_parent, root=lexical_root)
    require_no_symlink_components_below(path, lexical_root)
    require_no_symlink_components_below(policy_parent, lexical_root)
    root_descriptor = os.open(
        lexical_root,
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        fcntl.flock(root_descriptor, fcntl.LOCK_EX)
        try:
            _require_same_directory(lexical_root, root_descriptor)
            policy_descriptor = os.open(
                policy_parent.name,
                os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=root_descriptor,
            )
            try:
                fcntl.flock(policy_descriptor, fcntl.LOCK_EX)
                try:
                    _require_same_directory(lexical_root, root_descriptor)
                    _require_same_directory(policy_parent, policy_descriptor)
                    descriptor = os.open(
                        path.name,
                        os.O_APPEND
                        | os.O_CREAT
                        | os.O_WRONLY
                        | getattr(os, "O_NOFOLLOW", 0),
                        0o600,
                        dir_fd=root_descriptor,
                    )
                    try:
                        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                            raise ValueError("Policy promotion lock is not a regular file")
                        fcntl.flock(descriptor, fcntl.LOCK_EX)
                        try:
                            _require_same_regular_file_at(
                                root_descriptor,
                                path.name,
                                descriptor,
                                label="Policy promotion lock",
                            )
                            yield policy_descriptor
                        finally:
                            try:
                                _require_same_directory(lexical_root, root_descriptor)
                                _require_same_directory(
                                    policy_parent, policy_descriptor
                                )
                                _require_same_regular_file_at(
                                    root_descriptor,
                                    path.name,
                                    descriptor,
                                    label="Policy promotion lock",
                                )
                            finally:
                                fcntl.flock(descriptor, fcntl.LOCK_UN)
                    finally:
                        os.close(descriptor)
                finally:
                    try:
                        _require_same_directory(lexical_root, root_descriptor)
                        _require_same_directory(policy_parent, policy_descriptor)
                    finally:
                        fcntl.flock(policy_descriptor, fcntl.LOCK_UN)
            finally:
                os.close(policy_descriptor)
        finally:
            try:
                _require_same_directory(lexical_root, root_descriptor)
            finally:
                fcntl.flock(root_descriptor, fcntl.LOCK_UN)
    finally:
        os.close(root_descriptor)


@contextmanager
def _promotion_lock(authorized_pic_root: Path) -> Iterator[int]:
    anchor = stable_serialization_anchor(authorized_pic_root)
    durable_mkdir_parents(anchor)
    anchor_descriptor = os.open(
        anchor,
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        fcntl.flock(anchor_descriptor, fcntl.LOCK_EX)
        try:
            _require_same_directory(anchor, anchor_descriptor)
            with _promotion_lock_within_serialization_anchor(
                authorized_pic_root
            ) as policy_descriptor:
                _require_same_directory(anchor, anchor_descriptor)
                yield policy_descriptor
                _require_same_directory(anchor, anchor_descriptor)
        finally:
            try:
                _require_same_directory(anchor, anchor_descriptor)
            finally:
                fcntl.flock(anchor_descriptor, fcntl.LOCK_UN)
    finally:
        os.close(anchor_descriptor)


@contextmanager
def _serialization_anchor_lock(authorized_pic_root: Path) -> Iterator[None]:
    anchor = stable_serialization_anchor(authorized_pic_root)
    anchor_descriptor = os.open(
        anchor,
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        fcntl.flock(anchor_descriptor, fcntl.LOCK_EX)
        try:
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


def verify_active_launch_prohibited_generation(
    *,
    expected_control_plane_version: str,
    expected_active_policy_sha256: str,
    expected_active_promotion_sha256: str,
    expected_authorized_freeze_manifest: Path,
    expected_authorized_freeze_manifest_sha256: str,
    expected_authorized_freeze_build_controller: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
) -> dict[str, object]:
    """Revalidate one exact launch-prohibited active generation without mutation."""
    control_plane_dir = SCRIPT_DIR
    version = _lowercase_sha256(
        expected_control_plane_version,
        label="Expected active control-plane version",
    )
    expected_snapshot = {
        "active_policy_sha256": _lowercase_sha256(
            expected_active_policy_sha256,
            label="Expected active-policy SHA-256",
        ),
        "active_promotion_sha256": _lowercase_sha256(
            expected_active_promotion_sha256,
            label="Expected active-promotion SHA-256",
        ),
    }
    expected_manifest = Path(expected_authorized_freeze_manifest)
    if expected_manifest != Path(os.path.abspath(expected_manifest)):
        raise ValueError("Expected authorized freeze manifest is not absolute")
    expected_freeze = {
        "status": "authorized",
        "manifest_path": str(expected_manifest),
        "manifest_sha256": _lowercase_sha256(
            expected_authorized_freeze_manifest_sha256,
            label="Expected authorized freeze manifest SHA-256",
        ),
        "build_profile_control_plane_version": _lowercase_sha256(
            expected_authorized_freeze_build_controller,
            label="Expected authorized freeze build controller",
        ),
    }

    executing_inventory = verify_installed_control_plane(
        control_plane_dir,
        authorized_pic_root=authorized_pic_root,
    )
    if executing_inventory.get("version") != version:
        raise ValueError("Executing control-plane inventory differs from expected")
    installed_inventory: dict[str, object] | None = None
    for root in [authorized_pic_root, authorized_project_home_root]:
        lexical_root = Path(os.path.abspath(root))
        inventory = verify_historical_installed_control_plane(
            lexical_root / "control_plane" / version,
            authorized_pic_root=root,
        )
        if inventory.get("version") != version:
            raise ValueError("Installed control-plane inventory differs from expected")
        if installed_inventory is not None and inventory != installed_inventory:
            raise ValueError(
                "Installed Orion and Project Home control-plane inventories differ"
            )
        installed_inventory = inventory
    if installed_inventory != executing_inventory:
        raise ValueError("Executing and installed control-plane inventories differ")

    def require_quiescent(policy: dict[str, object]) -> None:
        require_no_incomplete_manual_accounting_marker(
            Path(os.path.abspath(authorized_pic_root)) / "ledger" / "node_hours.jsonl",
            project_home_ledger_root(authorized_project_home_root)
            / "ledger"
            / "node_hours.jsonl",
        )
        _require_no_outstanding_submissions(
            policy, authorized_pic_root, authorized_project_home_root
        )

    def require_expected_active() -> dict[str, object]:
        policy, snapshot = require_storage_policy_unlock_snapshot(
            control_plane_version=version,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
            authorized_account=authorized_account,
            allow_pending_genesis=True,
        )
        if snapshot != expected_snapshot:
            raise ValueError("Active launch-prohibited generation differs from expected")
        if policy.get("registered_science_slices") != []:
            raise ValueError(
                "Active launch-prohibited generation has registered science slices"
            )
        if policy.get("science_submission_freeze") != expected_freeze:
            raise ValueError(
                "Active launch-prohibited authorized freeze differs from expected"
            )
        return policy

    with _serialization_anchor_lock(authorized_pic_root):
        policy = require_expected_active()
        require_quiescent(policy)
        _require_authorized_clean_candidate_freeze_revalidation(
            policy,
            control_plane_version=version,
            control_plane_dir=control_plane_dir,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        if require_expected_active() != policy:
            raise ValueError("Active launch-prohibited policy changed during verification")
        require_quiescent(policy)
        return {
            **expected_snapshot,
            "control_plane_version": version,
            "record_type": (
                "frontier_pic_active_launch_prohibited_generation_verification"
            ),
            "schema_version": 1,
            "science_submission_freeze": expected_freeze,
            "status": "passed",
        }


def promote(
    reviewed_policy: Path,
    *,
    pre_policy_promotion_attestation: Path | None = None,
    pre_policy_promotion_authorization_id: str | None = None,
    retire_historical_storage_preflight_predecessor: bool = False,
    migrate_exact_reviewed_storage_preflight_predecessor: bool = False,
    replace_exact_authorized_clean_candidate_freeze: bool = False,
    retire_completed_q043_registered_slices: bool = False,
    retire_completed_q023_registered_slices: bool = False,
    q043_registered_matrix: Path | None = None,
    q043_registered_matrix_sha256: str | None = None,
    q023_registered_matrix: Path | None = None,
    q023_registered_matrix_sha256: str | None = None,
    expected_active_policy_sha256: str | None = None,
    expected_active_promotion_sha256: str | None = None,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
) -> dict[str, object]:
    if (
        sum(
            [
                retire_historical_storage_preflight_predecessor,
                migrate_exact_reviewed_storage_preflight_predecessor,
                replace_exact_authorized_clean_candidate_freeze,
                retire_completed_q043_registered_slices,
                retire_completed_q023_registered_slices,
            ]
        )
        > 1
    ):
        raise ValueError("Policy predecessor transition modes are exclusive")
    inventory = verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    version = str(inventory["version"])
    verify_installed_control_plane(
        Path(os.path.abspath(authorized_project_home_root)) / "control_plane" / version,
        authorized_pic_root=authorized_project_home_root,
    )
    reviewed_bytes = read_stable_regular_file(reviewed_policy)
    policy = read_json_bytes(reviewed_bytes, label=str(reviewed_policy))
    validate_storage_policy(
        policy,
        control_plane_version=version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        allow_pending_genesis=True,
    )

    policy_path = canonical_policy_path(authorized_pic_root)
    mirror_policy_path = canonical_policy_path(authorized_project_home_root)
    require_no_symlink_components_below(policy_path, authorized_pic_root)
    require_no_symlink_components_below(
        mirror_policy_path, authorized_project_home_root
    )
    record = {
        "schema_version": 2,
        "promotion_id": str(uuid.uuid4()),
        "control_plane_version": version,
        "policy_path": str(policy_path),
        "project_home_policy_path": str(mirror_policy_path),
        "policy_sha256": sha256_bytes(reviewed_bytes),
    }
    registered_science_slices = policy.get("registered_science_slices")
    if isinstance(registered_science_slices, list) and registered_science_slices:
        if (
            pre_policy_promotion_attestation is None
            or pre_policy_promotion_authorization_id is None
        ):
            raise ValueError(
                "Registered-science policy promotion requires a sealed "
                "pre-policy-promotion attestation"
            )
        attestation = validate_sealed_operator_attestation(
            pre_policy_promotion_attestation,
            authorization_id=pre_policy_promotion_authorization_id,
            phase="pre_policy_promotion",
            control_plane_version=version,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=project_home_ledger_root(
                authorized_project_home_root
            ),
        )
        record.update(
            {
                "pre_policy_promotion_attestation_authorization_id": (
                    pre_policy_promotion_authorization_id
                ),
                "pre_policy_promotion_attestation_path": attestation["path"],
                "pre_policy_promotion_attestation_sha256": attestation["sha256"],
            }
        )
    elif (
        pre_policy_promotion_attestation is not None
        or pre_policy_promotion_authorization_id is not None
    ):
        raise ValueError(
            "Launch-prohibited policy promotion must not claim a "
            "pre-policy-promotion attestation"
        )
    promotion_bytes = (
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    successor_snapshot = {
        "active_policy_sha256": sha256_bytes(reviewed_bytes),
        "active_promotion_sha256": sha256_bytes(promotion_bytes),
    }
    if retire_historical_storage_preflight_predecessor and (
        registered_science_slices != []
        or policy.get("science_submission_freeze")
        != {"status": "pending_clean_candidate_freeze"}
    ):
        raise ValueError(
            "Storage-preflight predecessor migration requires one launch-prohibited "
            "pending-freeze replacement"
        )
    if (
        migrate_exact_reviewed_storage_preflight_predecessor
        and registered_science_slices != []
    ):
        raise ValueError(
            "Exact reviewed storage-preflight predecessor migration requires one "
            "empty-allowlist replacement"
        )
    if (
        replace_exact_authorized_clean_candidate_freeze
        and registered_science_slices != []
    ):
        raise ValueError(
            "Exact authorized clean-candidate freeze replacement requires one "
            "empty-allowlist replacement"
        )
    if retire_completed_q043_registered_slices and registered_science_slices != []:
        raise ValueError(
            "Completed Q043 retirement requires one empty-allowlist replacement"
        )
    if retire_completed_q023_registered_slices and registered_science_slices != []:
        raise ValueError(
            "Completed Q023 retirement requires one empty-allowlist replacement"
        )
    mirror_policy_parent = Path(os.path.abspath(authorized_project_home_root)) / "policy"
    durable_mkdir_parents(mirror_policy_parent, root=authorized_project_home_root)
    with _promotion_lock(authorized_pic_root) as policy_descriptor, (
        _locked_mirror_policy_parent(
            mirror_policy_parent,
            authorized_project_home_root=authorized_project_home_root,
        )
    ) as mirror_policy_descriptor:
        def require_pinned_policy_parents() -> None:
            require_same_directory(
                policy_path.parent, policy_descriptor, root=authorized_pic_root
            )
            require_same_directory(
                mirror_policy_parent,
                mirror_policy_descriptor,
                root=authorized_project_home_root,
            )

        def require_predecessor_snapshot(
            *,
            predecessor_policy_bytes: bytes | None = None,
            predecessor_promotion_bytes: bytes | None = None,
            allow_active_promotion_transaction: bool = False,
        ) -> tuple[dict[str, object], dict[str, str]] | None:
            return require_policy_predecessor_snapshot_for_promotion(
                successor_policy=policy,
                successor_control_plane_version=version,
                permit_historical_retirement_predecessor=(
                    retire_historical_storage_preflight_predecessor
                ),
                permit_exact_reviewed_storage_preflight_predecessor=(
                    migrate_exact_reviewed_storage_preflight_predecessor
                ),
                permit_exact_authorized_clean_candidate_freeze_replacement=(
                    replace_exact_authorized_clean_candidate_freeze
                ),
                permit_completed_q043_registered_slice_retirement=(
                    retire_completed_q043_registered_slices
                ),
                permit_completed_q023_registered_slice_retirement=(
                    retire_completed_q023_registered_slices
                ),
                q043_registered_matrix_path=q043_registered_matrix,
                q043_registered_matrix_sha256=q043_registered_matrix_sha256,
                q023_registered_matrix_path=q023_registered_matrix,
                q023_registered_matrix_sha256=q023_registered_matrix_sha256,
                expected_active_policy_sha256=expected_active_policy_sha256,
                expected_active_promotion_sha256=expected_active_promotion_sha256,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
                authorized_account=authorized_account,
                predecessor_policy_bytes=predecessor_policy_bytes,
                predecessor_promotion_bytes=predecessor_promotion_bytes,
                allow_active_promotion_transaction=allow_active_promotion_transaction,
            )

        def validate_transaction_predecessor(
            predecessor_policy_bytes: bytes | None,
            predecessor_promotion_bytes: bytes | None,
        ) -> tuple[dict[str, object], dict[str, str]] | None:
            return require_predecessor_snapshot(
                predecessor_policy_bytes=predecessor_policy_bytes,
                predecessor_promotion_bytes=predecessor_promotion_bytes,
                allow_active_promotion_transaction=True,
            )

        def validate_active_successor(
            active_policy: dict[str, object],
            active_control_plane_version: str,
        ) -> None:
            _require_authorized_clean_candidate_freeze_revalidation(
                active_policy,
                control_plane_version=active_control_plane_version,
                control_plane_dir=(
                    Path(os.path.abspath(authorized_pic_root))
                    / "control_plane"
                    / active_control_plane_version
                ),
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )

        def validate_reviewed_successor_candidate() -> None:
            _require_authorized_clean_candidate_freeze_revalidation(
                policy,
                control_plane_version=version,
                control_plane_dir=control_plane_dir,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )

        def validate_transaction_successor() -> None:
            active_policy, active_snapshot = require_storage_policy_unlock_snapshot(
                control_plane_version=version,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
                authorized_account=authorized_account,
                allow_pending_genesis=True,
                allow_active_promotion_transaction=True,
            )
            if active_snapshot != successor_snapshot or active_policy != policy:
                raise ValueError(
                    "Active-policy promotion semantic successor snapshot differs"
                )
            _require_authorized_clean_candidate_freeze_revalidation(
                active_policy,
                control_plane_version=version,
                control_plane_dir=control_plane_dir,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )

        require_pinned_policy_parents()
        _recover_promotion_transaction(
            policy_descriptor,
            mirror_policy_descriptor,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
            authorized_account=authorized_account,
            post_publish_check=require_pinned_policy_parents,
            validate_predecessor=validate_transaction_predecessor,
            validate_active_successor=validate_active_successor,
        )

        predecessor = require_predecessor_snapshot()
        predecessor_snapshot = predecessor[1] if predecessor is not None else None
        if (
            migrate_exact_reviewed_storage_preflight_predecessor
            or replace_exact_authorized_clean_candidate_freeze
            or retire_completed_q043_registered_slices
            or retire_completed_q023_registered_slices
        ):
            science_freeze = policy["science_submission_freeze"]
            assert isinstance(science_freeze, dict)
            validate_reviewed_successor_candidate()
        _require_q011_launch_prohibited_baseline(
            policy,
            control_plane_version=version,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
            authorized_account=authorized_account,
        )
        require_no_incomplete_manual_accounting_marker(
            Path(os.path.abspath(authorized_pic_root)) / "ledger" / "node_hours.jsonl",
            project_home_ledger_root(authorized_project_home_root)
            / "ledger"
            / "node_hours.jsonl",
        )
        _require_no_outstanding_submissions(
            policy, authorized_pic_root, authorized_project_home_root
        )
        if isinstance(registered_science_slices, list) and registered_science_slices:
            assert pre_policy_promotion_attestation is not None
            assert pre_policy_promotion_authorization_id is not None
            attestation = validate_sealed_operator_attestation(
                pre_policy_promotion_attestation,
                authorization_id=pre_policy_promotion_authorization_id,
                phase="pre_policy_promotion",
                control_plane_version=version,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=project_home_ledger_root(
                    authorized_project_home_root
                ),
            )
            if attestation != {
                "path": record["pre_policy_promotion_attestation_path"],
                "sha256": record["pre_policy_promotion_attestation_sha256"],
            }:
                raise ValueError(
                    "Pre-policy-promotion attestation changed while acquiring lock"
                )
        require_pinned_policy_parents()
        if isinstance(registered_science_slices, list) and registered_science_slices:
            assert pre_policy_promotion_attestation is not None
            assert pre_policy_promotion_authorization_id is not None
            attestation = validate_sealed_operator_attestation(
                pre_policy_promotion_attestation,
                authorization_id=pre_policy_promotion_authorization_id,
                phase="pre_policy_promotion",
                control_plane_version=version,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=project_home_ledger_root(
                    authorized_project_home_root
                ),
            )
            if attestation != {
                "path": record["pre_policy_promotion_attestation_path"],
                "sha256": record["pre_policy_promotion_attestation_sha256"],
            }:
                raise ValueError(
                    "Pre-policy-promotion attestation changed while acquiring "
                    "mirror lock"
                )
        current_predecessor = require_predecessor_snapshot()
        current_predecessor_snapshot = (
            current_predecessor[1] if current_predecessor is not None else None
        )
        if current_predecessor_snapshot != predecessor_snapshot:
            raise ValueError(
                "Active-policy predecessor changed before transaction publication"
            )
        with _active_policy_transaction(
            policy_descriptor,
            mirror_policy_descriptor,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
            post_publish_check=require_pinned_policy_parents,
            expected_predecessor_snapshot=predecessor_snapshot,
            expected_successor_snapshot=successor_snapshot,
            validate_predecessor=validate_transaction_predecessor,
            validate_successor=validate_transaction_successor,
        ):
            atomic_write_bytes_at(
                policy_descriptor,
                policy_path.name,
                reviewed_bytes,
                post_publish_check=require_pinned_policy_parents,
            )
            atomic_write_bytes_at(
                mirror_policy_descriptor,
                mirror_policy_path.name,
                reviewed_bytes,
                post_publish_check=require_pinned_policy_parents,
            )
            atomic_write_json_at(
                policy_descriptor,
                active_promotion_path(authorized_pic_root).name,
                record,
                post_publish_check=require_pinned_policy_parents,
            )
            atomic_write_json_at(
                mirror_policy_descriptor,
                active_promotion_path(authorized_project_home_root).name,
                record,
                post_publish_check=require_pinned_policy_parents,
            )
            require_storage_policy_unlock_snapshot(
                control_plane_version=version,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
                authorized_account=authorized_account,
                allow_pending_genesis=True,
                allow_active_promotion_transaction=True,
            )
            require_pinned_policy_parents()
    print(record["policy_sha256"])
    return record


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--reviewed-policy", type=Path)
    mode.add_argument(
        "--verify-active-launch-prohibited-generation",
        action="store_true",
    )
    parser.add_argument("--pre-policy-promotion-attestation", type=Path)
    parser.add_argument("--pre-policy-promotion-authorization-id")
    parser.add_argument(
        "--retire-historical-storage-preflight-predecessor",
        action="store_true",
    )
    parser.add_argument(
        "--migrate-exact-reviewed-storage-preflight-predecessor",
        action="store_true",
    )
    parser.add_argument(
        "--replace-exact-authorized-clean-candidate-freeze",
        action="store_true",
    )
    parser.add_argument(
        "--retire-completed-q043-registered-slices",
        action="store_true",
    )
    parser.add_argument(
        "--retire-completed-q023-registered-slices",
        action="store_true",
    )
    parser.add_argument("--q043-registered-matrix", type=Path)
    parser.add_argument("--q043-registered-matrix-sha256")
    parser.add_argument("--q023-registered-matrix", type=Path)
    parser.add_argument("--q023-registered-matrix-sha256")
    parser.add_argument("--expected-active-policy-sha256")
    parser.add_argument("--expected-active-promotion-sha256")
    parser.add_argument("--expected-control-plane-version")
    parser.add_argument("--expected-authorized-freeze-manifest", type=Path)
    parser.add_argument("--expected-authorized-freeze-manifest-sha256")
    parser.add_argument("--expected-authorized-freeze-build-controller")
    return parser


def main() -> None:
    parser = _parser()
    args = parser.parse_args()
    if args.verify_active_launch_prohibited_generation:
        if any(
            value is not None
            for value in (
                args.pre_policy_promotion_attestation,
                args.pre_policy_promotion_authorization_id,
                args.q043_registered_matrix,
                args.q043_registered_matrix_sha256,
                args.q023_registered_matrix,
                args.q023_registered_matrix_sha256,
            )
        ) or any(
            (
                args.retire_historical_storage_preflight_predecessor,
                args.migrate_exact_reviewed_storage_preflight_predecessor,
                args.replace_exact_authorized_clean_candidate_freeze,
                args.retire_completed_q043_registered_slices,
                args.retire_completed_q023_registered_slices,
            )
        ):
            parser.error("Active-generation verification does not accept promotion flags")
        required = {
            "expected_active_policy_sha256": args.expected_active_policy_sha256,
            "expected_active_promotion_sha256": args.expected_active_promotion_sha256,
            "expected_control_plane_version": args.expected_control_plane_version,
            "expected_authorized_freeze_manifest": (
                args.expected_authorized_freeze_manifest
            ),
            "expected_authorized_freeze_manifest_sha256": (
                args.expected_authorized_freeze_manifest_sha256
            ),
            "expected_authorized_freeze_build_controller": (
                args.expected_authorized_freeze_build_controller
            ),
        }
        missing = [name for name, value in required.items() if value is None]
        if missing:
            parser.error(
                "Active-generation verification requires: " + ", ".join(missing)
            )
        result = verify_active_launch_prohibited_generation(**required)
        print(
            json.dumps(
                result,
                allow_nan=False,
                ensure_ascii=True,
                separators=(",", ":"),
                sort_keys=True,
            )
        )
        return
    if any(
        value is not None
        for value in (
            args.expected_control_plane_version,
            args.expected_authorized_freeze_manifest,
            args.expected_authorized_freeze_manifest_sha256,
            args.expected_authorized_freeze_build_controller,
        )
    ):
        parser.error("Policy promotion does not accept verifier-only bindings")
    assert args.reviewed_policy is not None
    promote(
        args.reviewed_policy,
        pre_policy_promotion_attestation=args.pre_policy_promotion_attestation,
        pre_policy_promotion_authorization_id=(
            args.pre_policy_promotion_authorization_id
        ),
        retire_historical_storage_preflight_predecessor=(
            args.retire_historical_storage_preflight_predecessor
        ),
        migrate_exact_reviewed_storage_preflight_predecessor=(
            args.migrate_exact_reviewed_storage_preflight_predecessor
        ),
        replace_exact_authorized_clean_candidate_freeze=(
            args.replace_exact_authorized_clean_candidate_freeze
        ),
        retire_completed_q043_registered_slices=(
            args.retire_completed_q043_registered_slices
        ),
        retire_completed_q023_registered_slices=(
            args.retire_completed_q023_registered_slices
        ),
        q043_registered_matrix=args.q043_registered_matrix,
        q043_registered_matrix_sha256=args.q043_registered_matrix_sha256,
        q023_registered_matrix=args.q023_registered_matrix,
        q023_registered_matrix_sha256=args.q023_registered_matrix_sha256,
        expected_active_policy_sha256=args.expected_active_policy_sha256,
        expected_active_promotion_sha256=args.expected_active_promotion_sha256,
    )


if __name__ == "__main__":
    main()
