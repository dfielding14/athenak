#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Create the explicit mirrored genesis event for the Frontier PIC ledger."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
from pathlib import Path
import re

from control_plane_common import AUTHORIZED_ACCOUNT, AUTHORIZED_PIC_ROOT
from control_plane_common import AUTHORIZED_PROJECT_HOME_ROOT, require_ledger_paths
from control_plane_common import active_promotion_path, canonical_policy_path
from control_plane_common import read_json_bytes, read_stable_regular_file_below
from control_plane_common import sha256_bytes
from control_plane_common import require_storage_policy_unlock
from control_plane_common import verify_installed_control_plane
from ledger import initialize_ledger, migrate_existing_genesis_anchors


SCRIPT_DIR = Path(__file__).absolute().parent


def initialize_from_policy(
    *,
    ledger_jsonl: Path,
    ledger_csv: Path,
    mirror_receipts: Path,
    mirror_jsonl: Path,
    mirror_transport: str,
    notes: str,
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
        authorized_project_home_root / "control_plane" / version,
        authorized_pic_root=authorized_project_home_root,
    )
    require_ledger_paths(
        ledger_jsonl,
        ledger_csv,
        mirror_receipts,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    policy = require_storage_policy_unlock(
        control_plane_version=version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        ledger_mirror_transport=mirror_transport,
        allow_pending_genesis=True,
    )
    storage = policy["olcf_side_storage"]
    if (
        not isinstance(storage, dict)
        or storage.get("ledger_genesis_allowed") is not True
        or storage.get("ledger_genesis") is not None
    ):
        raise ValueError("Active storage policy does not authorize ledger genesis")
    return initialize_ledger(
        ledger_jsonl,
        ledger_csv,
        mirror_receipts,
        mirror_jsonl,
        mirror_transport=mirror_transport,
        notes=notes,
        control_plane_version=version,
    )


def migrate_anchor_from_active_policy(
    *,
    ledger_jsonl: Path,
    ledger_csv: Path,
    mirror_receipts: Path,
    mirror_jsonl: Path,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    """Migrate one audited predecessor ledger into paired immutable anchors."""
    inventory = verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    verify_installed_control_plane(
        authorized_project_home_root / "control_plane" / str(inventory["version"]),
        authorized_pic_root=authorized_project_home_root,
    )
    require_ledger_paths(
        ledger_jsonl,
        ledger_csv,
        mirror_receipts,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    policy_path = canonical_policy_path(authorized_pic_root)
    mirror_policy_path = canonical_policy_path(authorized_project_home_root)
    promotion_path = active_promotion_path(authorized_pic_root)
    mirror_promotion_path = active_promotion_path(authorized_project_home_root)
    policy_bytes = read_stable_regular_file_below(
        policy_path, authorized_pic_root, require_read_only_mode=True
    )
    mirror_policy_bytes = read_stable_regular_file_below(
        mirror_policy_path,
        authorized_project_home_root,
        require_read_only_mode=True,
    )
    promotion_bytes = read_stable_regular_file_below(
        promotion_path, authorized_pic_root, require_read_only_mode=True
    )
    mirror_promotion_bytes = read_stable_regular_file_below(
        mirror_promotion_path,
        authorized_project_home_root,
        require_read_only_mode=True,
    )
    if policy_bytes != mirror_policy_bytes or promotion_bytes != mirror_promotion_bytes:
        raise ValueError("Active predecessor policy mirrors differ")
    promotion = read_json_bytes(promotion_bytes, label=str(promotion_path))
    if (
        set(promotion) != {
            "schema_version",
            "control_plane_version",
            "policy_path",
            "project_home_policy_path",
            "policy_sha256",
        }
        or promotion.get("schema_version") != 1
        or not re.fullmatch(
            r"[0-9a-f]{64}", str(promotion.get("control_plane_version", ""))
        )
        or Path(str(promotion.get("policy_path", ""))).resolve()
        != policy_path.resolve()
        or Path(str(promotion.get("project_home_policy_path", ""))).resolve()
        != mirror_policy_path.resolve()
        or promotion.get("policy_sha256") != sha256_bytes(policy_bytes)
    ):
        raise ValueError("Active predecessor promotion anchor is malformed")
    policy = read_json_bytes(policy_bytes, label=str(policy_path))
    storage = policy.get("olcf_side_storage")
    genesis = storage.get("ledger_genesis") if isinstance(storage, dict) else None
    if (
        not isinstance(genesis, dict)
        or genesis.get("status") != "initialized"
        or genesis.get("mirror_transport") != "filesystem_copy"
        or not re.fullmatch(r"[0-9a-f]{64}", str(genesis.get("event_sha256", "")))
        or not re.fullmatch(r"[0-9a-f]{64}", str(genesis.get("mirror_ack_sha256", "")))
    ):
        raise ValueError("Active predecessor policy has no audited genesis anchor")
    return migrate_existing_genesis_anchors(
        ledger_jsonl,
        mirror_receipts,
        mirror_jsonl,
        expected_event_sha256=str(genesis["event_sha256"]),
        expected_mirror_ack_sha256=str(genesis["mirror_ack_sha256"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger-jsonl", required=True, type=Path)
    parser.add_argument("--ledger-csv", required=True, type=Path)
    parser.add_argument("--mirror-receipts", required=True, type=Path)
    parser.add_argument("--mirror-jsonl", required=True, type=Path)
    parser.add_argument("--mirror-transport", default="filesystem_copy")
    parser.add_argument("--notes")
    parser.add_argument("--migrate-existing-anchor", action="store_true")
    args = parser.parse_args()
    if args.migrate_existing_anchor:
        anchor = migrate_anchor_from_active_policy(
            ledger_jsonl=args.ledger_jsonl,
            ledger_csv=args.ledger_csv,
            mirror_receipts=args.mirror_receipts,
            mirror_jsonl=args.mirror_jsonl,
        )
        print(anchor["event_sha256"])
    else:
        if not args.notes:
            parser.error("--notes is required unless --migrate-existing-anchor is selected")
        event = initialize_from_policy(
            ledger_jsonl=args.ledger_jsonl,
            ledger_csv=args.ledger_csv,
            mirror_receipts=args.mirror_receipts,
            mirror_jsonl=args.mirror_jsonl,
            mirror_transport=args.mirror_transport,
            notes=args.notes,
        )
        print(event["event_sha256"])


if __name__ == "__main__":
    main()
