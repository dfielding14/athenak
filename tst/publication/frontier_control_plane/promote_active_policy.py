#!/opt/cray/pe/python/3.11.7/bin/python3
"""Atomically promote one reviewed policy into the mirrored active-policy anchor."""

from __future__ import annotations

import argparse
from pathlib import Path

from control_plane_common import AUTHORIZED_ACCOUNT, AUTHORIZED_PIC_ROOT
from control_plane_common import AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import active_promotion_path, atomic_write_bytes, atomic_write_json
from control_plane_common import canonical_policy_path, read_json, sha256
from control_plane_common import require_no_symlink_components_below
from control_plane_common import validate_storage_policy, verify_installed_control_plane


SCRIPT_DIR = Path(__file__).absolute().parent


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
        authorized_project_home_root.resolve() / "control_plane" / version,
        authorized_pic_root=authorized_project_home_root,
    )
    reviewed_bytes = reviewed_policy.read_bytes()
    policy = read_json(reviewed_policy)
    validate_storage_policy(
        policy,
        control_plane_version=version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
    )

    policy_path = canonical_policy_path(authorized_pic_root)
    mirror_policy_path = canonical_policy_path(authorized_project_home_root)
    require_no_symlink_components_below(policy_path, authorized_pic_root.resolve())
    require_no_symlink_components_below(
        mirror_policy_path, authorized_project_home_root.resolve()
    )
    atomic_write_bytes(policy_path, reviewed_bytes)
    atomic_write_bytes(mirror_policy_path, reviewed_bytes)
    record = {
        "schema_version": 1,
        "control_plane_version": version,
        "policy_path": str(policy_path),
        "project_home_policy_path": str(mirror_policy_path),
        "policy_sha256": sha256(policy_path),
    }
    atomic_write_json(active_promotion_path(authorized_pic_root), record)
    atomic_write_json(active_promotion_path(authorized_project_home_root), record)
    print(record["policy_sha256"])
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reviewed-policy", required=True, type=Path)
    args = parser.parse_args()
    promote(args.reviewed_policy)


if __name__ == "__main__":
    main()
