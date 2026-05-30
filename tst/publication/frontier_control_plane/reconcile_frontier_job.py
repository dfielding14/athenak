#!/opt/cray/pe/python/3.11.7/bin/python3
"""Reconcile a completed Frontier PIC job into the node-hour ledger."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess

from control_plane_common import AUTHORIZED_ACCOUNT, AUTHORIZED_PIC_ROOT
from control_plane_common import AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import TRUSTED_SACCT
from control_plane_common import require_ledger_paths, verify_installed_control_plane
from ledger import accounting, append_primary_event_locked, ledger_lock
from ledger import latest_reservations, require_explicit_genesis, transition_payload
from ledger import validate_mirrored_state


ACTIVE_STATES = {"CONFIGURING", "COMPLETING", "PENDING", "RUNNING"}
TERMINAL_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "COMPLETED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "TIMEOUT",
}
SCRIPT_DIR = Path(__file__).absolute().parent


def _scheduler_result(job_id: str, reservation_id: str) -> tuple[str, int, int]:
    output = subprocess.check_output(
        [
            TRUSTED_SACCT,
            "-j",
            job_id,
            "--allocations",
            "--noheader",
            "--parsable2",
            "--format=JobIDRaw,State,ElapsedRaw,AllocNodes,Comment,Account",
        ],
        text=True,
    )
    for line in output.splitlines():
        fields = line.split("|")
        if len(fields) >= 6 and fields[0] == job_id:
            if fields[4] != f"pic-reservation={reservation_id}":
                raise ValueError("Slurm accounting comment does not bind the PIC reservation")
            if fields[5] != AUTHORIZED_ACCOUNT:
                raise ValueError("Slurm accounting account does not match the PIC account")
            return fields[1].split()[0], int(fields[2]), int(fields[3])
    raise ValueError(f"No Slurm allocation record found for job {job_id}")


def reconcile(
    *,
    job_id: str,
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    require_ledger_paths(
        ledger_jsonl,
        ledger_csv,
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    with ledger_lock(ledger_jsonl):
        inventory = verify_installed_control_plane(
            control_plane_dir, authorized_pic_root=authorized_pic_root
        )
        records = validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
        require_explicit_genesis(records)
        latest = [
            record for record in latest_reservations(records).values()
            if record.get("job_id") == job_id
        ]
        if len(latest) != 1 or latest[0].get("state") != "submitted":
            raise ValueError(f"Expected one submitted reservation for job {job_id}")
        state, elapsed_seconds, allocated_nodes = _scheduler_result(
            job_id, str(latest[0]["reservation_id"])
        )
        state = state.upper()
        if state in ACTIVE_STATES:
            raise ValueError(f"Job is not in a terminal Slurm state: {state}")
        if state not in TERMINAL_STATES:
            raise ValueError(f"Unrecognized terminal Slurm state: {state}")
        if elapsed_seconds < 0 or allocated_nodes < 0:
            raise ValueError("Scheduler accounting values must not be negative")

        submitted = latest[0]
        if submitted.get("control_plane_version") != inventory["version"]:
            raise ValueError("Reservation belongs to another control-plane version")
        requested_nodes = int(submitted["requested_nodes"])
        billed_nodes = max(requested_nodes, allocated_nodes)
        consumed = billed_nodes * elapsed_seconds / 3600.0
        cumulative = accounting(records)["cumulative_consumed_node_hours"] + consumed
        event = transition_payload(submitted)
        event.update(
            {
                "event_type": "reconciliation",
                "state": state,
                "reconciled": True,
                "scheduler_reported_allocated_nodes": allocated_nodes,
                "billed_nodes": billed_nodes,
                "elapsed_seconds": elapsed_seconds,
                "consumed_node_hours": consumed,
                "cumulative_consumed_node_hours": cumulative,
            }
        )
        return append_primary_event_locked(
            ledger_jsonl,
            ledger_csv,
            receipts_jsonl,
            mirror_jsonl,
            event,
            mirror_transport="filesystem_copy",
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--ledger-jsonl", required=True, type=Path)
    parser.add_argument("--ledger-csv", required=True, type=Path)
    parser.add_argument("--receipts-jsonl", required=True, type=Path)
    parser.add_argument("--mirror-jsonl", required=True, type=Path)
    args = parser.parse_args()
    event = reconcile(
        job_id=args.job_id,
        ledger_jsonl=args.ledger_jsonl,
        ledger_csv=args.ledger_csv,
        receipts_jsonl=args.receipts_jsonl,
        mirror_jsonl=args.mirror_jsonl,
    )
    print(event["consumed_node_hours"])


if __name__ == "__main__":
    main()
