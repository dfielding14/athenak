#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Reconcile a completed Frontier PIC job into the node-hour ledger."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess

from control_plane_common import AUTHORIZED_PIC_ROOT
from control_plane_common import AUTHORIZED_SLURM_CLUSTER
from control_plane_common import AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import TRUSTED_SACCT
from control_plane_common import PinnedDirectoryAncestry, atomic_write_bytes_at
from control_plane_common import read_json, read_json_bytes, record_for_role
from control_plane_common import read_stable_regular_file_below, require_ledger_paths
from control_plane_common import require_canonical_path_below
from control_plane_common import scheduler_account_matches_authorized
from control_plane_common import validate_planner_retention_binding
from control_plane_common import trusted_slurm_environment
from control_plane_common import verify_historical_installed_control_plane
from control_plane_common import verify_installed_control_plane
from ledger import accounting, append_primary_event_locked, ledger_lock
from ledger import latest_reservations, require_explicit_genesis, transition_payload
from ledger import validate_mirrored_state
from terminal_recovery_handoff import PURGED_CANCELLED_ZERO_EXECUTION_MODE
from terminal_recovery_handoff import require_closed_received_marker
from terminal_recovery_handoff import require_purged_cancelled_zero_execution_snapshot
from terminal_recovery_handoff import verify_terminal_recovery_handoff
from validate_and_reserve_frontier_job import _clear_matching_pending_marker
from validate_and_reserve_frontier_job import _matching_pending_marker
from validate_and_reserve_frontier_job import _pending_marker_path
from validate_and_reserve_frontier_job import _require_current_reservation_marker
from validate_and_reserve_frontier_job import _require_reservation_policy_snapshot
from validate_and_reserve_frontier_job import _require_run_artifact_dir
from validate_and_reserve_frontier_job import _verify_scheduler_job_binding


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
REGISTERED_EXECUTION_RECEIPT_NAME = "q011_section54_registered_execution_receipt.json"
REGISTERED_EXECUTION_RECEIPT_RECORD_TYPE = (
    "q011_section54_reconciled_registered_execution_receipt"
)
REGISTERED_EXECUTION_RECEIPT_ROLE = "immutable_reconciled_registered_execution"


def _scheduler_result(
    job_id: str,
    reservation_id: str,
    *,
    allow_purged_cancelled_zero_execution: bool = False,
) -> tuple[str, int, int, str]:
    output = subprocess.check_output(
        [
            TRUSTED_SACCT,
            "-j",
            job_id,
            "--allocations",
            f"--clusters={AUTHORIZED_SLURM_CLUSTER}",
            "--noheader",
            "--parsable2",
            "--format=JobIDRaw,State,ElapsedRaw,AllocNodes,Comment,Account,ExitCode",
        ],
        text=True,
        env=trusted_slurm_environment(),
    )
    matching = [
        fields
        for line in output.splitlines()
        if (fields := line.split("|")) and fields[0] == job_id
    ]
    if len(matching) != 1 or len(matching[0]) != 7:
        raise ValueError(f"Expected one exact Slurm allocation record for job {job_id}")
    fields = matching[0]
    if fields[4]:
        if fields[4] != f"pic-reservation={reservation_id}":
            raise ValueError("Slurm accounting comment does not bind the PIC reservation")
    else:
        try:
            _verify_scheduler_job_binding(job_id, reservation_id)
        except subprocess.CalledProcessError:
            if not allow_purged_cancelled_zero_execution:
                raise
            snapshot = require_purged_cancelled_zero_execution_snapshot(job_id)
            if (
                fields[1] != snapshot["state"]
                or int(fields[2]) != snapshot["elapsed_raw"]
                or int(fields[3]) != snapshot["allocated_nodes"]
                or fields[5] != snapshot["account"]
                or fields[6] != snapshot["exit_code"]
            ):
                raise ValueError("Slurm accounting changed during terminal recovery")
    if not scheduler_account_matches_authorized(fields[5]):
        raise ValueError("Slurm accounting account does not match the PIC account")
    if re.fullmatch(r"[0-9]+:[0-9]+", fields[6]) is None:
        raise ValueError("Slurm accounting exit code is malformed")
    return fields[1].split()[0], int(fields[2]), int(fields[3]), fields[6]


def _verify_reservation_control_plane_pair(
    reservation: dict[str, object],
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> None:
    version = reservation.get("control_plane_version")
    if not isinstance(version, str) or re.fullmatch(r"[0-9a-f]{64}", version) is None:
        raise ValueError("Reservation has an invalid control-plane version")
    for root in [authorized_pic_root, authorized_project_home_root]:
        inventory = verify_historical_installed_control_plane(
            root / "control_plane" / version,
            authorized_pic_root=root,
        )
        if inventory["version"] != version:
            raise ValueError("Reservation control-plane inventory digest mismatch")


def _json_bytes(value: dict[str, object]) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _registered_execution_receipt_payload(
    event: dict[str, object], *, authorized_pic_root: Path
) -> tuple[Path, bytes] | None:
    """Derive one Q011 receipt solely from immutable manifest and reconciliation."""
    retention = event.get("planner_retention")
    if retention is None or event.get("state") != "COMPLETED":
        return None
    planner_retention = validate_planner_retention_binding(
        retention,
        authorized_pic_root=authorized_pic_root,
        expected_clean_candidate_manifest_sha256=str(
            event["clean_candidate_manifest_sha256"]
        ),
    )
    manifest_path = require_canonical_path_below(
        Path(str(event.get("manifest_path", ""))),
        Path(os.path.abspath(authorized_pic_root)) / "manifests",
    )
    manifest_bytes = read_stable_regular_file_below(
        manifest_path,
        Path(os.path.abspath(authorized_pic_root)),
        require_read_only_mode=True,
    )
    if hashlib.sha256(manifest_bytes).hexdigest() != event.get("manifest_sha256"):
        raise ValueError("Registered-execution manifest digest differs from ledger")
    manifest = read_json_bytes(manifest_bytes, label="registered-execution manifest")
    manifest_retention = validate_planner_retention_binding(
        manifest.get("planner_retention"),
        authorized_pic_root=authorized_pic_root,
        expected_clean_candidate_manifest_sha256=str(
            manifest["clean_candidate_manifest_sha256"]
        ),
    )
    artifact_dir = _require_run_artifact_dir(manifest)
    if (
        planner_retention != manifest_retention
        or manifest.get("submission_scope") != "registered_science"
        or event.get("submission_scope") != "registered_science"
        or event.get("reconciled") is not True
        or event.get("submission_id") != manifest.get("submission_id")
        or event.get("control_plane_version") != manifest.get("control_plane_version")
        or event.get("git_commit") != manifest.get("git_commit")
        or event.get("artifact_dir") != str(artifact_dir)
    ):
        raise ValueError("Registered-execution reconciliation differs from manifest")
    receipt = {
        "record_type": REGISTERED_EXECUTION_RECEIPT_RECORD_TYPE,
        "schema_version": 1,
        "receipt_role": REGISTERED_EXECUTION_RECEIPT_ROLE,
        "registration_scope": "registered_science",
        "reconciled": True,
        "reservation_id": event["reservation_id"],
        "submission_id": event["submission_id"],
        "reconciliation_event_sha256": event["event_sha256"],
        "attempt_id": planner_retention["attempt_id"],
        "source_commit": manifest["git_commit"],
        "executable_sha256": record_for_role(manifest, "executable")["sha256"],
        "deck_sha256": record_for_role(manifest, "input-deck")["sha256"],
        "environment_sha256": record_for_role(manifest, "environment-profile")[
            "sha256"
        ],
        "control_plane_version": manifest["control_plane_version"],
        "argv": planner_retention["argv"],
        "slurm_job_id": event["job_id"],
        "slurm_terminal_state": event["state"],
        "raw_output_root": planner_retention["authorized_orion_raw_root"],
        "artifact_dir": str(artifact_dir),
        "planner_retention": planner_retention,
        "pre_submit_manifest_sha256": event["manifest_sha256"],
    }
    return artifact_dir / "analysis" / REGISTERED_EXECUTION_RECEIPT_NAME, _json_bytes(
        receipt
    )


def _read_pinned_read_only_receipt(parent_fd: int, name: str) -> bytes:
    descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent_fd)
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_mode & 0o222
        ):
            raise ValueError("Registered-execution receipt is not one read-only file")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            payload = stream.read()
        after = os.fstat(descriptor)
        current = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        identity = lambda value: (
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_nlink,
            value.st_size,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
        if (
            identity(before) != identity(after)
            or (after.st_dev, after.st_ino) != (current.st_dev, current.st_ino)
        ):
            raise ValueError("Registered-execution receipt changed while reading")
        return payload
    finally:
        os.close(descriptor)


def _publish_registered_execution_receipt(
    event: dict[str, object], *, authorized_pic_root: Path
) -> Path | None:
    """Publish or verify the post-mirror Q011 receipt inside pinned /runs analysis."""
    derived = _registered_execution_receipt_payload(
        event, authorized_pic_root=authorized_pic_root
    )
    if derived is None:
        return None
    path, payload = derived
    analysis_dir = path.parent
    with PinnedDirectoryAncestry(
        analysis_dir, root=Path(os.path.abspath(authorized_pic_root))
    ) as ancestry:
        metadata = os.fstat(ancestry.descriptor)
        if not stat.S_ISDIR(metadata.st_mode) or stat.S_IMODE(metadata.st_mode) != 0o700:
            raise ValueError("Registered-execution analysis directory is not private")
        try:
            atomic_write_bytes_at(
                ancestry.descriptor,
                path.name,
                payload,
                mode=0o444,
                replace=False,
                post_publish_check=ancestry.require_same,
            )
        except FileExistsError:
            pass
        ancestry.require_same()
        if _read_pinned_read_only_receipt(ancestry.descriptor, path.name) != payload:
            raise ValueError("Registered-execution receipt retry bytes differ")
        ancestry.require_same()
    return path


def _require_cross_generation_recovery(
    reservation: dict[str, object],
    *,
    job_id: str,
    executing_version: str,
    terminal_recovery_handoff: Path | None,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    require_marker: bool = True,
    refresh_scheduler_binding_recovery: bool = True,
) -> tuple[Path, str, str, dict[str, object]] | None:
    origin_version = str(reservation.get("control_plane_version", ""))
    if origin_version == executing_version:
        if terminal_recovery_handoff is not None:
            raise ValueError("Terminal-recovery handoff was supplied for a current reservation")
        return None
    if terminal_recovery_handoff is None:
        raise ValueError("Cross-generation recovery requires an immutable mirrored handoff")
    handoff, handoff_sha256 = verify_terminal_recovery_handoff(
        terminal_recovery_handoff,
        reservation=reservation,
        job_id=job_id,
        recovery_control_plane_version=executing_version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        require_marker=require_marker,
        refresh_scheduler_binding_recovery=refresh_scheduler_binding_recovery,
    )
    scheduler_binding_recovery = handoff.get("scheduler_binding_recovery")
    mode = (
        str(scheduler_binding_recovery["mode"])
        if isinstance(scheduler_binding_recovery, dict)
        else "fresh_scheduler_binding"
    )
    return terminal_recovery_handoff, handoff_sha256, mode, handoff


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
    terminal_recovery_handoff: Path | None = None,
) -> dict[str, object]:
    require_ledger_paths(
        ledger_jsonl,
        ledger_csv,
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    inventory = verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    verify_installed_control_plane(
        authorized_project_home_root / "control_plane" / str(inventory["version"]),
        authorized_pic_root=authorized_project_home_root,
    )
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        inventory = verify_installed_control_plane(
            control_plane_dir, authorized_pic_root=authorized_pic_root
        )
        verify_installed_control_plane(
            authorized_project_home_root / "control_plane" / str(inventory["version"]),
            authorized_pic_root=authorized_project_home_root,
        )
        records = validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
        require_explicit_genesis(records)
        latest_by_reservation = latest_reservations(records)
        latest = [
            record for record in latest_by_reservation.values()
            if record.get("job_id") == job_id
        ]
        if (
            len(latest) == 1
            and latest[0].get("event_type") == "reconciliation"
            and latest[0].get("reconciled") is True
        ):
            _verify_reservation_control_plane_pair(
                latest[0],
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )
            if latest[0].get("control_plane_version") == inventory["version"]:
                _require_reservation_policy_snapshot(
                    latest[0],
                    authorized_pic_root=authorized_pic_root,
                    authorized_project_home_root=authorized_project_home_root,
                )
            marker_path = _pending_marker_path(authorized_pic_root)
            recovery = _require_cross_generation_recovery(
                latest[0],
                job_id=job_id,
                executing_version=str(inventory["version"]),
                terminal_recovery_handoff=terminal_recovery_handoff,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
                require_marker=marker_path.is_file(),
                refresh_scheduler_binding_recovery=False,
            )
            if recovery is not None and (
                latest[0].get("reconciled_by_control_plane_version")
                != inventory["version"]
                or latest[0].get("terminal_recovery_handoff_path") != str(recovery[0])
                or latest[0].get("terminal_recovery_handoff_sha256") != recovery[1]
                or latest[0].get("terminal_recovery_mode") != recovery[2]
            ):
                raise ValueError("Cross-generation retry belongs to another recovery handoff")
            _publish_registered_execution_receipt(
                latest[0], authorized_pic_root=authorized_pic_root
            )
            if marker_path.is_file():
                if recovery is not None:
                    marker = _matching_pending_marker(
                        marker_path, str(latest[0]["reservation_id"])
                    )
                    if marker is None:
                        raise ValueError("Cross-generation retry is missing its pending marker")
                    require_closed_received_marker(marker, latest[0], job_id=job_id)
                _clear_matching_pending_marker(
                    marker_path, str(latest[0]["reservation_id"])
                )
            return latest[0]
        pending_attachment = False
        pending_requires_scheduler_binding = False
        marker: dict[str, object] | None = None
        if not latest:
            marker_path = _pending_marker_path(authorized_pic_root)
            read_marker = read_json(marker_path) if marker_path.is_file() else None
            marker = (
                _matching_pending_marker(marker_path, str(read_marker["reservation_id"]))
                if read_marker is not None
                else None
            )
            if (
                marker is None
                or marker.get("state")
                not in {"scheduler_job_id_received", "submitted_not_attached"}
                or marker.get("job_id") != job_id
            ):
                raise ValueError(f"Expected one submitted reservation for job {job_id}")
            reserved = latest_by_reservation.get(str(marker["reservation_id"]))
            if reserved is None or reserved.get("state") != "reserved":
                raise ValueError("Pending scheduler attachment is not a reserved job")
            if marker.get("state") == "scheduler_job_id_received":
                if reserved.get("control_plane_version") == inventory["version"]:
                    _require_current_reservation_marker(
                        marker,
                        reserved,
                        control_plane_version=str(inventory["version"]),
                    )
                require_closed_received_marker(marker, reserved, job_id=job_id)
            elif reserved.get("control_plane_version") == inventory["version"]:
                _require_current_reservation_marker(
                    marker,
                    reserved,
                    control_plane_version=str(inventory["version"]),
                )
            latest = [reserved]
            pending_attachment = True
            pending_requires_scheduler_binding = (
                marker.get("state") == "scheduler_job_id_received"
            )
        if len(latest) != 1 or latest[0].get("state") not in {"reserved", "submitted"}:
            raise ValueError(f"Expected one submitted reservation for job {job_id}")
        _verify_reservation_control_plane_pair(
            latest[0],
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        if latest[0].get("control_plane_version") == inventory["version"]:
            _require_reservation_policy_snapshot(
                latest[0],
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )
        elif marker is None:
            marker = _matching_pending_marker(
                _pending_marker_path(authorized_pic_root),
                str(latest[0]["reservation_id"]),
            )
        recovery_resume_candidate = (
            terminal_recovery_handoff is not None
            and latest[0].get("event_type") == "job_id_attached"
            and latest[0].get("state") == "submitted"
            and latest[0].get("attached_by_control_plane_version")
            == inventory["version"]
            and latest[0].get("terminal_recovery_handoff_path")
            == str(terminal_recovery_handoff)
            and latest[0].get("terminal_recovery_mode")
            == PURGED_CANCELLED_ZERO_EXECUTION_MODE
        )
        recovery = _require_cross_generation_recovery(
            latest[0],
            job_id=job_id,
            executing_version=str(inventory["version"]),
            terminal_recovery_handoff=terminal_recovery_handoff,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
            refresh_scheduler_binding_recovery=not recovery_resume_candidate,
        )
        recovery_resume = (
            recovery is not None
            and latest[0].get("event_type") == "job_id_attached"
            and latest[0].get("state") == "submitted"
            and latest[0].get("attached_by_control_plane_version")
            == inventory["version"]
            and latest[0].get("terminal_recovery_handoff_path") == str(recovery[0])
            and latest[0].get("terminal_recovery_handoff_sha256") == recovery[1]
            and latest[0].get("terminal_recovery_mode") == recovery[2]
        )
        if recovery is not None and (
            not (pending_attachment or recovery_resume)
            or marker is None
            or marker.get("state") != "scheduler_job_id_received"
        ):
            raise ValueError(
                "Cross-generation recovery requires a scheduler-job-id-received marker"
            )
        if pending_requires_scheduler_binding and (
            recovery is None or recovery[2] != PURGED_CANCELLED_ZERO_EXECUTION_MODE
        ):
            _verify_scheduler_job_binding(job_id, str(latest[0]["reservation_id"]))
        if recovery_resume and recovery[2] == PURGED_CANCELLED_ZERO_EXECUTION_MODE:
            snapshot = recovery[3]["scheduler_binding_recovery"]
            if not isinstance(snapshot, dict):
                raise ValueError("Exceptional terminal recovery snapshot is missing")
            state = str(snapshot["state"]).split()[0]
            elapsed_seconds = int(snapshot["elapsed_raw"])
            allocated_nodes = int(snapshot["allocated_nodes"])
            scheduler_exit_code = str(snapshot["exit_code"])
        else:
            (
                state,
                elapsed_seconds,
                allocated_nodes,
                scheduler_exit_code,
            ) = _scheduler_result(
                job_id,
                str(latest[0]["reservation_id"]),
                allow_purged_cancelled_zero_execution=(
                    recovery is not None
                    and recovery[2] == PURGED_CANCELLED_ZERO_EXECUTION_MODE
                ),
            )
        state = state.upper()
        if state in ACTIVE_STATES:
            raise ValueError(f"Job is not in a terminal Slurm state: {state}")
        if state not in TERMINAL_STATES:
            raise ValueError(f"Unrecognized terminal Slurm state: {state}")
        if elapsed_seconds < 0 or allocated_nodes < 0:
            raise ValueError("Scheduler accounting values must not be negative")

        submitted = latest[0]
        requested_nodes = int(submitted["requested_nodes"])
        billed_nodes = max(requested_nodes, allocated_nodes)
        consumed = billed_nodes * elapsed_seconds / 3600.0
        if pending_attachment:
            attachment = transition_payload(submitted)
            attachment.update(
                {
                    "event_type": "job_id_attached",
                    "state": "submitted",
                    "job_id": job_id,
                    "attached_by_control_plane_version": str(inventory["version"]),
                }
            )
            if recovery is not None:
                attachment.update(
                    {
                        "terminal_recovery_handoff_path": str(recovery[0]),
                        "terminal_recovery_handoff_sha256": recovery[1],
                        "terminal_recovery_mode": recovery[2],
                    }
                )
            submitted = append_primary_event_locked(
                ledger_jsonl,
                ledger_csv,
                receipts_jsonl,
                mirror_jsonl,
                attachment,
                mirror_transport="filesystem_copy",
            )
        cumulative = accounting(records)["cumulative_consumed_node_hours"] + consumed
        event = transition_payload(submitted)
        event.update(
            {
                "event_type": "reconciliation",
                "state": state,
                "reconciled": True,
                "reconciled_by_control_plane_version": str(inventory["version"]),
                "scheduler_reported_allocated_nodes": allocated_nodes,
                "billed_nodes": billed_nodes,
                "elapsed_seconds": elapsed_seconds,
                "scheduler_exit_code": scheduler_exit_code,
                "consumed_node_hours": consumed,
                "cumulative_consumed_node_hours": cumulative,
            }
        )
        if recovery is not None:
            event.update(
                {
                    "terminal_recovery_handoff_path": str(recovery[0]),
                    "terminal_recovery_handoff_sha256": recovery[1],
                    "terminal_recovery_mode": recovery[2],
                }
            )
        result = append_primary_event_locked(
            ledger_jsonl,
            ledger_csv,
            receipts_jsonl,
            mirror_jsonl,
            event,
            mirror_transport="filesystem_copy",
        )
        _publish_registered_execution_receipt(
            result, authorized_pic_root=authorized_pic_root
        )
        marker_path = _pending_marker_path(authorized_pic_root)
        if marker_path.is_file():
            _clear_matching_pending_marker(
                marker_path, str(result["reservation_id"])
            )
        return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--ledger-jsonl", required=True, type=Path)
    parser.add_argument("--ledger-csv", required=True, type=Path)
    parser.add_argument("--receipts-jsonl", required=True, type=Path)
    parser.add_argument("--mirror-jsonl", required=True, type=Path)
    parser.add_argument("--terminal-recovery-handoff", type=Path)
    args = parser.parse_args()
    event = reconcile(
        job_id=args.job_id,
        ledger_jsonl=args.ledger_jsonl,
        ledger_csv=args.ledger_csv,
        receipts_jsonl=args.receipts_jsonl,
        mirror_jsonl=args.mirror_jsonl,
        terminal_recovery_handoff=args.terminal_recovery_handoff,
    )
    print(event["consumed_node_hours"])


if __name__ == "__main__":
    main()
