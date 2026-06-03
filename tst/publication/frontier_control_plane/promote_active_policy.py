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
import os
from pathlib import Path
import stat
from typing import Iterator

from control_plane_common import AUTHORIZED_ACCOUNT, AUTHORIZED_PIC_ROOT
from control_plane_common import AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import active_promotion_path, atomic_write_bytes_at
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
from ledger import _path_exists, _pinned_parent_directories
from ledger import latest_reservations, validate_mirrored_state, validate_receipts
from ledger import require_no_incomplete_manual_accounting_marker
from operator_attestation import validate_sealed_operator_attestation


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


def promote(
    reviewed_policy: Path,
    *,
    pre_policy_promotion_attestation: Path | None = None,
    pre_policy_promotion_authorization_id: str | None = None,
    retire_historical_storage_preflight_predecessor: bool = False,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
) -> dict[str, object]:
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
        "schema_version": 1,
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
    if retire_historical_storage_preflight_predecessor and (
        registered_science_slices != []
        or policy.get("science_submission_freeze")
        != {"status": "pending_clean_candidate_freeze"}
    ):
        raise ValueError(
            "Historical storage-preflight retirement requires one launch-prohibited "
            "pending-freeze replacement"
        )
    mirror_policy_parent = Path(os.path.abspath(authorized_project_home_root)) / "policy"
    durable_mkdir_parents(mirror_policy_parent, root=authorized_project_home_root)
    with _promotion_lock(authorized_pic_root) as policy_descriptor:
        require_policy_predecessor_snapshot_for_promotion(
            successor_policy=policy,
            successor_control_plane_version=version,
            permit_historical_retirement_predecessor=(
                retire_historical_storage_preflight_predecessor
            ),
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
            authorized_account=authorized_account,
        )
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
        mirror_policy_descriptor = open_directory_below(
            mirror_policy_parent, root=authorized_project_home_root
        )
        try:
            fcntl.flock(mirror_policy_descriptor, fcntl.LOCK_EX)
            try:
                require_same_directory(
                    mirror_policy_parent,
                    mirror_policy_descriptor,
                    root=authorized_project_home_root,
                )
                def require_pinned_policy_parents() -> None:
                    require_same_directory(
                        policy_path.parent, policy_descriptor, root=authorized_pic_root
                    )
                    require_same_directory(
                        mirror_policy_parent,
                        mirror_policy_descriptor,
                        root=authorized_project_home_root,
                    )

                require_pinned_policy_parents()
                if (
                    isinstance(registered_science_slices, list)
                    and registered_science_slices
                ):
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
                            "Pre-policy-promotion attestation changed while "
                            "acquiring mirror lock"
                        )
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
                )
                require_pinned_policy_parents()
            finally:
                fcntl.flock(mirror_policy_descriptor, fcntl.LOCK_UN)
        finally:
            os.close(mirror_policy_descriptor)
    print(record["policy_sha256"])
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reviewed-policy", required=True, type=Path)
    parser.add_argument("--pre-policy-promotion-attestation", type=Path)
    parser.add_argument("--pre-policy-promotion-authorization-id")
    parser.add_argument(
        "--retire-historical-storage-preflight-predecessor",
        action="store_true",
    )
    args = parser.parse_args()
    promote(
        args.reviewed_policy,
        pre_policy_promotion_attestation=args.pre_policy_promotion_attestation,
        pre_policy_promotion_authorization_id=(
            args.pre_policy_promotion_authorization_id
        ),
        retire_historical_storage_preflight_predecessor=(
            args.retire_historical_storage_preflight_predecessor
        ),
    )


if __name__ == "__main__":
    main()
