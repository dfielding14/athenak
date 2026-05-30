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
from control_plane_common import active_promotion_path, atomic_write_bytes, atomic_write_json
from control_plane_common import canonical_policy_path, durable_mkdir_parents
from control_plane_common import read_json_bytes
from control_plane_common import read_stable_regular_file, sha256_bytes
from control_plane_common import require_no_symlink_components_below
from control_plane_common import require_storage_policy_unlock_snapshot
from control_plane_common import validate_storage_policy, verify_installed_control_plane


SCRIPT_DIR = Path(__file__).absolute().parent


@contextmanager
def _promotion_lock(authorized_pic_root: Path) -> Iterator[None]:
    lexical_root = Path(os.path.abspath(authorized_pic_root))
    path = lexical_root / "policy" / ".promotion.lock"
    durable_mkdir_parents(path.parent, root=lexical_root)
    require_no_symlink_components_below(path, lexical_root)
    descriptor = os.open(
        path,
        os.O_APPEND | os.O_CREAT | os.O_WRONLY | getattr(os, "O_NOFOLLOW", 0),
        0o600,
    )
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise ValueError("Policy promotion lock is not a regular file")
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
    finally:
        os.close(descriptor)


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
    with _promotion_lock(authorized_pic_root):
        atomic_write_bytes(
            policy_path, reviewed_bytes, root=authorized_pic_root
        )
        atomic_write_bytes(
            mirror_policy_path,
            reviewed_bytes,
            root=authorized_project_home_root,
        )
        atomic_write_json(
            active_promotion_path(authorized_pic_root),
            record,
            root=authorized_pic_root,
        )
        atomic_write_json(
            active_promotion_path(authorized_project_home_root),
            record,
            root=authorized_project_home_root,
        )
        require_storage_policy_unlock_snapshot(
            control_plane_version=version,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
            authorized_account=authorized_account,
            allow_pending_genesis=True,
        )
    print(record["policy_sha256"])
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reviewed-policy", required=True, type=Path)
    args = parser.parse_args()
    promote(args.reviewed_policy)


if __name__ == "__main__":
    main()
