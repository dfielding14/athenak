#!/opt/cray/pe/python/3.11.7/bin/python3
"""Create the explicit mirrored genesis event for the Frontier PIC ledger."""

from __future__ import annotations

import argparse
from pathlib import Path

from control_plane_common import AUTHORIZED_ACCOUNT, AUTHORIZED_PIC_ROOT
from control_plane_common import AUTHORIZED_PROJECT_HOME_ROOT, require_ledger_paths
from control_plane_common import require_storage_policy_unlock
from control_plane_common import verify_installed_control_plane
from ledger import initialize_ledger


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
    require_storage_policy_unlock(
        control_plane_version=version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        ledger_mirror_transport=mirror_transport,
    )
    return initialize_ledger(
        ledger_jsonl,
        ledger_csv,
        mirror_receipts,
        mirror_jsonl,
        mirror_transport=mirror_transport,
        notes=notes,
        control_plane_version=version,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger-jsonl", required=True, type=Path)
    parser.add_argument("--ledger-csv", required=True, type=Path)
    parser.add_argument("--mirror-receipts", required=True, type=Path)
    parser.add_argument("--mirror-jsonl", required=True, type=Path)
    parser.add_argument("--mirror-transport", default="filesystem_copy")
    parser.add_argument("--notes", required=True)
    args = parser.parse_args()
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
