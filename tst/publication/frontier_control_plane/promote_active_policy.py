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
from control_plane_common import read_json_bytes
from control_plane_common import read_stable_regular_file, sha256_bytes
from control_plane_common import require_no_symlink_components_below
from control_plane_common import require_storage_policy_unlock_snapshot
from control_plane_common import stable_serialization_anchor
from control_plane_common import validate_storage_policy, verify_installed_control_plane
from ledger import _path_exists, _pinned_parent_directories
from ledger import latest_reservations, validate_mirrored_state
from ledger import require_no_incomplete_manual_accounting_marker


SCRIPT_DIR = Path(__file__).absolute().parent


def _require_no_outstanding_submissions(
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> None:
    ledger_root = Path(os.path.abspath(authorized_pic_root)) / "ledger"
    project_ledger_root = Path(os.path.abspath(authorized_project_home_root)) / "ledger"
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
            return
        if not all(existing):
            raise ValueError("Policy promotion is blocked by incomplete mirrored ledger state")
        try:
            records = validate_mirrored_state(ledger, receipts, mirror)
        except ValueError as lexical_error:
            if resolved_mirror == mirror:
                raise
            try:
                records = validate_mirrored_state(ledger, receipts, resolved_mirror)
            except ValueError:
                raise lexical_error
        outstanding = [
            record for record in latest_reservations(records).values()
            if record.get("state") in {"reserved", "submitted"}
        ]
        if outstanding:
            raise ValueError("Policy promotion is blocked by an outstanding reservation")


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
    mirror_policy_parent = Path(os.path.abspath(authorized_project_home_root)) / "policy"
    durable_mkdir_parents(mirror_policy_parent, root=authorized_project_home_root)
    with _promotion_lock(authorized_pic_root) as policy_descriptor:
        require_no_incomplete_manual_accounting_marker(
            Path(os.path.abspath(authorized_pic_root)) / "ledger" / "node_hours.jsonl",
            Path(os.path.abspath(authorized_project_home_root))
            / "ledger"
            / "node_hours.jsonl",
        )
        _require_no_outstanding_submissions(
            authorized_pic_root, authorized_project_home_root
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
    args = parser.parse_args()
    promote(args.reviewed_policy)


if __name__ == "__main__":
    main()
