#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Publish and verify a mirrored immutable terminal-recovery authorization."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import uuid

from control_plane_common import AUTHORIZED_PIC_ROOT, AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import AUTHORIZED_SLURM_CLUSTER
from control_plane_common import TRUSTED_SACCT, TRUSTED_SCONTROL, TRUSTED_SQUEUE
from control_plane_common import atomic_write_bytes, durable_mkdir_parents
from control_plane_common import read_json_bytes, read_stable_regular_file_below
from control_plane_common import require_canonical_path_below, require_ledger_paths
from control_plane_common import require_storage_policy_unlock_snapshot
from control_plane_common import scheduler_account_matches_authorized
from control_plane_common import sha256_bytes, verify_installed_control_plane
from control_plane_common import verify_historical_installed_control_plane
from control_plane_common import trusted_slurm_environment
from ledger import latest_reservations, ledger_lock, require_explicit_genesis
from ledger import validate_mirrored_state


SCRIPT_DIR = Path(__file__).absolute().parent
HANDOFF_STATUS = "authorized_terminal_scheduler_job_id_received_recovery"
PURGED_CANCELLED_ZERO_EXECUTION_MODE = "purged_scontrol_cancelled_zero_execution"
REVIEWED_PURGED_BINDING_ATTESTATION_MODE = (
    "reviewed_operator_attestation_for_unprovable_purged_reservation_job_binding"
)
TRUSTED_SUBMISSION_WRAPPER_JOB_NAME = "run_installed_control_plane_job.sh"


def _validate_purged_cancelled_zero_execution_snapshot(
    snapshot: object, *, job_id: str
) -> dict[str, object]:
    if not isinstance(snapshot, dict) or set(snapshot) != {
        "mode",
        "job_id",
        "job_name",
        "state",
        "elapsed_raw",
        "allocated_nodes",
        "comment",
        "account",
        "submit",
        "start",
        "end",
        "exit_code",
    }:
        raise ValueError("Purged scheduler-recovery snapshot is malformed")
    state = snapshot.get("state")
    submit = snapshot.get("submit")
    if (
        snapshot.get("mode") != PURGED_CANCELLED_ZERO_EXECUTION_MODE
        or snapshot.get("job_id") != job_id
        or snapshot.get("job_name") != TRUSTED_SUBMISSION_WRAPPER_JOB_NAME
        or not isinstance(state, str)
        or not state.split()
        or state.split()[0] != "CANCELLED"
        or type(snapshot.get("elapsed_raw")) is not int
        or snapshot.get("elapsed_raw") != 0
        or type(snapshot.get("allocated_nodes")) is not int
        or snapshot.get("allocated_nodes") != 0
        or snapshot.get("comment") != ""
        or not scheduler_account_matches_authorized(snapshot.get("account"))
        or not isinstance(submit, str)
        or not submit
        or snapshot.get("start") != "None"
        or snapshot.get("end") != submit
        or snapshot.get("exit_code") != "0:0"
    ):
        raise ValueError("Slurm accounting is not a purged zero-execution cancellation")
    return snapshot


def require_purged_cancelled_zero_execution_snapshot(job_id: str) -> dict[str, object]:
    if re.fullmatch(r"[0-9]+", job_id) is None:
        raise ValueError("Terminal-recovery scheduler job ID is malformed")
    try:
        subprocess.check_output(
            [TRUSTED_SCONTROL, "show", "job", "--oneliner", job_id],
            text=True,
            stderr=subprocess.PIPE,
            env=trusted_slurm_environment(),
        )
    except subprocess.CalledProcessError as error:
        if (
            error.returncode != 1
            or str(error.stderr).strip()
            != "slurm_load_jobs error: Invalid job id specified"
        ):
            raise ValueError("scontrol failed for a reason other than scheduler-record purge") from error
    else:
        raise ValueError("Purged scheduler recovery requires an absent scontrol record")
    try:
        queued = subprocess.check_output(
            [TRUSTED_SQUEUE, "--jobs", job_id, "--noheader", "--format=%i"],
            text=True,
            stderr=subprocess.PIPE,
            env=trusted_slurm_environment(),
        )
    except subprocess.CalledProcessError as error:
        if (
            error.returncode != 1
            or str(error.stderr).strip()
            != "slurm_load_jobs error: Invalid job id specified"
        ):
            raise ValueError("squeue failed for a reason other than scheduler-record purge") from error
        queued = ""
    if queued.strip():
        raise ValueError("Purged scheduler recovery requires an empty live queue")
    output = subprocess.check_output(
        [
            TRUSTED_SACCT,
            "-j",
            job_id,
            "--allocations",
            f"--clusters={AUTHORIZED_SLURM_CLUSTER}",
            "--noheader",
            "--parsable2",
            "--format=JobIDRaw,JobName,State,ElapsedRaw,AllocNodes,Comment,Account,"
            "Submit,Start,End,ExitCode",
        ],
        text=True,
        env=trusted_slurm_environment(),
    )
    matching = [
        fields
        for line in output.splitlines()
        if (fields := line.split("|")) and fields[0] == job_id
    ]
    if len(matching) != 1 or len(matching[0]) != 11:
        raise ValueError(f"Expected one Slurm allocation record for job {job_id}")
    fields = matching[0]
    return _validate_purged_cancelled_zero_execution_snapshot(
        {
            "mode": PURGED_CANCELLED_ZERO_EXECUTION_MODE,
            "job_id": fields[0],
            "job_name": fields[1],
            "state": fields[2],
            "elapsed_raw": int(fields[3]),
            "allocated_nodes": int(fields[4]),
            "comment": fields[5],
            "account": fields[6],
            "submit": fields[7],
            "start": fields[8],
            "end": fields[9],
            "exit_code": fields[10],
        },
        job_id=job_id,
    )


def _verify_control_plane_pair(
    version: str,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    historical: bool = False,
) -> None:
    if re.fullmatch(r"[0-9a-f]{64}", version) is None:
        raise ValueError("Terminal-recovery control-plane version is malformed")
    for root in [authorized_pic_root, authorized_project_home_root]:
        verifier = (
            verify_historical_installed_control_plane
            if historical
            else verify_installed_control_plane
        )
        inventory = verifier(
            root / "control_plane" / version,
            authorized_pic_root=root,
        )
        if inventory["version"] != version:
            raise ValueError("Terminal-recovery control-plane inventory digest mismatch")


def _pending_marker_bytes(authorized_pic_root: Path) -> tuple[Path, bytes]:
    root = Path(os.path.abspath(authorized_pic_root))
    path = root / "ledger" / "pending_submission.json"
    require_canonical_path_below(path, root)
    return path, read_stable_regular_file_below(path, root, require_read_only_mode=True)


def require_closed_received_marker(
    marker: dict[str, object],
    reservation: dict[str, object],
    *,
    job_id: str,
) -> None:
    schema_version = marker.get("schema_version")
    if type(schema_version) is not int or schema_version not in {1, 2}:
        raise ValueError("Terminal-recovery marker schema is not supported")
    expected = {
        "schema_version": schema_version,
        "state": "scheduler_job_id_received",
        "reservation_id": reservation["reservation_id"],
        "submission_id": reservation["submission_id"],
        "manifest_path": reservation["manifest_path"],
        "manifest_sha256": reservation["manifest_sha256"],
        "job_id": job_id,
    }
    if schema_version == 2:
        expected["control_plane_version"] = reservation["control_plane_version"]
    if marker != expected:
        raise ValueError("Terminal-recovery marker does not match the reservation")


def _require_prior_policy_snapshot(
    reservation: dict[str, object],
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> dict[str, str]:
    _, snapshot = require_storage_policy_unlock_snapshot(
        control_plane_version=str(reservation["control_plane_version"]),
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    for field in ["active_policy_sha256", "active_promotion_sha256"]:
        if reservation.get(field) != snapshot[field]:
            raise ValueError("Terminal recovery reservation policy snapshot is not active")
    return snapshot


def _handoff_path(root: Path, handoff_id: str) -> Path:
    if str(uuid.UUID(handoff_id)) != handoff_id:
        raise ValueError("Terminal-recovery handoff ID must be a canonical UUID")
    lexical_root = Path(os.path.abspath(root))
    return lexical_root / "policy" / "recovery_handoffs" / f"{handoff_id}.json"


def _publish_copy(path: Path, payload: bytes, *, root: Path) -> None:
    durable_mkdir_parents(path.parent, root=root)
    if path.exists():
        existing = read_stable_regular_file_below(
            path, Path(os.path.abspath(root)), require_read_only_mode=True
        )
        if existing != payload:
            raise ValueError("Existing terminal-recovery handoff bytes differ")
        return
    atomic_write_bytes(path, payload, mode=0o400, replace=False, root=root)


def verify_terminal_recovery_handoff(
    handoff_path: Path,
    *,
    reservation: dict[str, object],
    job_id: str,
    recovery_control_plane_version: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    require_marker: bool = True,
    refresh_scheduler_binding_recovery: bool = True,
) -> tuple[dict[str, object], str]:
    lexical_root = Path(os.path.abspath(authorized_pic_root))
    handoff_path = require_canonical_path_below(handoff_path, lexical_root)
    if handoff_path.parent != lexical_root / "policy" / "recovery_handoffs":
        raise ValueError("Terminal-recovery handoff is outside the authorized directory")
    handoff_id = handoff_path.stem
    if _handoff_path(authorized_pic_root, handoff_id) != handoff_path:
        raise ValueError("Terminal-recovery handoff path is not canonical")
    mirror_path = _handoff_path(authorized_project_home_root, handoff_id)
    primary_bytes = read_stable_regular_file_below(
        handoff_path, lexical_root, require_read_only_mode=True
    )
    mirror_bytes = read_stable_regular_file_below(
        mirror_path,
        Path(os.path.abspath(authorized_project_home_root)),
        require_read_only_mode=True,
    )
    if primary_bytes != mirror_bytes:
        raise ValueError("Terminal-recovery handoff mirror bytes differ")
    marker_path, marker_bytes = _pending_marker_bytes(authorized_pic_root) if require_marker else (
        Path(""),
        b"",
    )
    marker_sha256 = ""
    if require_marker:
        marker = read_json_bytes(marker_bytes, label=str(marker_path))
        require_closed_received_marker(marker, reservation, job_id=job_id)
        marker_sha256 = sha256_bytes(marker_bytes)
    _verify_control_plane_pair(
        str(reservation["control_plane_version"]),
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        historical=True,
    )
    _verify_control_plane_pair(
        recovery_control_plane_version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    snapshot = _require_prior_policy_snapshot(
        reservation,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    handoff = read_json_bytes(primary_bytes, label=str(handoff_path))
    scheduler_binding_recovery = handoff.get("scheduler_binding_recovery")
    if scheduler_binding_recovery is not None:
        expected_scheduler_binding_recovery = _validate_purged_cancelled_zero_execution_snapshot(
            (
                require_purged_cancelled_zero_execution_snapshot(job_id)
                if refresh_scheduler_binding_recovery
                else scheduler_binding_recovery
            ),
            job_id=job_id,
        )
    else:
        expected_scheduler_binding_recovery = None
    expected = {
        "schema_version": 1,
        "status": HANDOFF_STATUS,
        "handoff_id": handoff_id,
        "reservation_id": reservation["reservation_id"],
        "submission_id": reservation["submission_id"],
        "job_id": job_id,
        "manifest_path": reservation["manifest_path"],
        "manifest_sha256": reservation["manifest_sha256"],
        "prior_control_plane_version": reservation["control_plane_version"],
        "recovery_control_plane_version": recovery_control_plane_version,
        "prior_active_policy_sha256": snapshot["active_policy_sha256"],
        "prior_active_promotion_sha256": snapshot["active_promotion_sha256"],
        "pending_marker_sha256": marker_sha256 or handoff.get("pending_marker_sha256"),
    }
    if expected_scheduler_binding_recovery is not None:
        expected["scheduler_binding_recovery"] = expected_scheduler_binding_recovery
        expected["reservation_job_binding_attestation"] = {
            "mode": REVIEWED_PURGED_BINDING_ATTESTATION_MODE,
            "reservation_id": reservation["reservation_id"],
            "job_id": job_id,
        }
    if type(handoff.get("schema_version")) is not int or handoff != expected:
        raise ValueError("Terminal-recovery handoff does not match the authorized recovery")
    return handoff, sha256_bytes(primary_bytes)


def create_handoff(
    *,
    job_id: str,
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    handoff_id: str | None = None,
    authorize_purged_cancelled_zero_execution: bool = False,
    attest_reviewed_purged_reservation_job_binding: bool = False,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> Path:
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
    recovery_version = str(inventory["version"])
    _verify_control_plane_pair(
        recovery_version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        records = validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
        require_explicit_genesis(records)
        marker_path, marker_bytes = _pending_marker_bytes(authorized_pic_root)
        marker = read_json_bytes(marker_bytes, label=str(marker_path))
        reservation = latest_reservations(records).get(str(marker.get("reservation_id", "")))
        if reservation is None or reservation.get("state") != "reserved":
            raise ValueError("Terminal-recovery marker is not bound to a reserved job")
        require_closed_received_marker(marker, reservation, job_id=job_id)
        if reservation["control_plane_version"] == recovery_version:
            raise ValueError("Terminal recovery handoff requires a successor control plane")
        _verify_control_plane_pair(
            str(reservation["control_plane_version"]),
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
            historical=True,
        )
        snapshot = _require_prior_policy_snapshot(
            reservation,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        selected_id = handoff_id or str(uuid.uuid4())
        primary_path = _handoff_path(authorized_pic_root, selected_id)
        payload = {
            "schema_version": 1,
            "status": HANDOFF_STATUS,
            "handoff_id": selected_id,
            "reservation_id": reservation["reservation_id"],
            "submission_id": reservation["submission_id"],
            "job_id": job_id,
            "manifest_path": reservation["manifest_path"],
            "manifest_sha256": reservation["manifest_sha256"],
            "prior_control_plane_version": reservation["control_plane_version"],
            "recovery_control_plane_version": recovery_version,
            "prior_active_policy_sha256": snapshot["active_policy_sha256"],
            "prior_active_promotion_sha256": snapshot["active_promotion_sha256"],
            "pending_marker_sha256": sha256_bytes(marker_bytes),
        }
        if (
            attest_reviewed_purged_reservation_job_binding
            and not authorize_purged_cancelled_zero_execution
        ):
            raise ValueError(
                "Purged reservation-job binding attestation requires exceptional recovery"
            )
        if authorize_purged_cancelled_zero_execution:
            if not attest_reviewed_purged_reservation_job_binding:
                raise ValueError(
                    "Purged recovery requires reviewed reservation-job binding attestation"
                )
            payload["scheduler_binding_recovery"] = (
                require_purged_cancelled_zero_execution_snapshot(job_id)
            )
            payload["reservation_job_binding_attestation"] = {
                "mode": REVIEWED_PURGED_BINDING_ATTESTATION_MODE,
                "reservation_id": reservation["reservation_id"],
                "job_id": job_id,
            }
        serialized = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
        _publish_copy(primary_path, serialized, root=authorized_pic_root)
        _publish_copy(
            _handoff_path(authorized_project_home_root, selected_id),
            serialized,
            root=authorized_project_home_root,
        )
        verify_terminal_recovery_handoff(
            primary_path,
            reservation=reservation,
            job_id=job_id,
            recovery_control_plane_version=recovery_version,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        return primary_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--ledger-jsonl", required=True, type=Path)
    parser.add_argument("--ledger-csv", required=True, type=Path)
    parser.add_argument("--receipts-jsonl", required=True, type=Path)
    parser.add_argument("--mirror-jsonl", required=True, type=Path)
    parser.add_argument("--handoff-id")
    parser.add_argument("--authorize-purged-cancelled-zero-execution", action="store_true")
    parser.add_argument(
        "--attest-reviewed-purged-reservation-job-binding",
        action="store_true",
    )
    args = parser.parse_args()
    print(create_handoff(**vars(args)))


if __name__ == "__main__":
    main()
