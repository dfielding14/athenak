#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Reconcile reviewed direct-srun Frontier allocations into the node-hour ledger."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
import os
from pathlib import Path
import pwd
import re
import subprocess

from control_plane_common import AUTHORIZED_PARTITION, AUTHORIZED_PIC_ROOT
from control_plane_common import AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import AUTHORIZED_SLURM_CLUSTER
from control_plane_common import TRUSTED_SACCT, TRUSTED_SQUEUE
from control_plane_common import read_json_bytes, read_stable_regular_file_below
from control_plane_common import require_canonical_path_below, require_ledger_paths
from control_plane_common import require_storage_policy_unlock_snapshot
from control_plane_common import scheduler_account_matches_authorized, sha256_bytes
from control_plane_common import trusted_slurm_environment
from control_plane_common import verify_installed_control_plane
from ledger import accounting, append_primary_event_locked, chain_head
from ledger import clear_matching_incomplete_manual_accounting_marker_locked
from ledger import ledger_lock, publish_incomplete_manual_accounting_marker_locked
from ledger import recover_incomplete_manual_accounting_locked
from ledger import latest_reservations, require_explicit_genesis
from ledger import validate_mirrored_state, validate_primary_chain, write_csv


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
MANUAL_ACCOUNTING_SCOPE = "manual_direct_srun_accounting_only"
SCRIPT_DIR = Path(__file__).absolute().parent


def _verify_installed_control_plane_pair(
    control_plane_dir: Path,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> dict[str, object]:
    inventory = verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    verify_installed_control_plane(
        Path(os.path.abspath(authorized_project_home_root))
        / "control_plane"
        / str(inventory["version"]),
        authorized_pic_root=authorized_project_home_root,
    )
    return inventory


def _authorization(
    path: Path,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> tuple[dict[str, object], str, Path, Path]:
    root = (
        Path(os.path.abspath(authorized_pic_root))
        / "policy"
        / "manual_accounting_authorizations"
    )
    path = require_canonical_path_below(path, root)
    data = read_stable_regular_file_below(
        path, root, require_read_only_mode=True
    )
    mirror_root = (
        Path(os.path.abspath(authorized_project_home_root))
        / "policy"
        / "manual_accounting_authorizations"
    )
    mirror_path = mirror_root / path.name
    mirror_data = read_stable_regular_file_below(
        mirror_path, mirror_root, require_read_only_mode=True
    )
    if data != mirror_data:
        raise ValueError("Manual-accounting authorization mirror bytes differ")
    authorization = read_json_bytes(data, label=str(path))
    if set(authorization) != {
        "schema_version",
        "authorization_id",
        "accounting_scope",
        "scientific_evidence_eligible",
        "jobs",
    }:
        raise ValueError("Manual-accounting authorization has an unsupported shape")
    if (
        type(authorization.get("schema_version")) is not int
        or authorization.get("schema_version") != 1
    ):
        raise ValueError("Manual-accounting authorization schema is unsupported")
    authorization_id = authorization.get("authorization_id")
    if (
        not isinstance(authorization_id, str)
        or re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,127}", authorization_id) is None
        or path.name != f"{authorization_id}.json"
    ):
        raise ValueError("Manual-accounting authorization ID is invalid")
    if authorization.get("accounting_scope") != MANUAL_ACCOUNTING_SCOPE:
        raise ValueError("Manual-accounting authorization scope is invalid")
    if authorization.get("scientific_evidence_eligible") is not False:
        raise ValueError("Manual-accounting authorization must reject scientific evidence use")
    jobs = authorization.get("jobs")
    if not isinstance(jobs, list) or not 1 <= len(jobs) <= 256:
        raise ValueError("Manual-accounting authorization jobs are invalid")
    job_ids: set[str] = set()
    for job in jobs:
        if not isinstance(job, dict) or set(job) != {"job_id", "expected_qos"}:
            raise ValueError("Manual-accounting authorization job is malformed")
        job_id = job.get("job_id")
        if (
            not isinstance(job_id, str)
            or re.fullmatch(r"[1-9][0-9]*", job_id) is None
            or job_id in job_ids
        ):
            raise ValueError("Manual-accounting authorization job ID is invalid")
        job_ids.add(job_id)
        if job.get("expected_qos") not in {"debug", "normal"}:
            raise ValueError("Manual-accounting authorization QoS is invalid")
    return authorization, sha256_bytes(data), path, mirror_path


def _require_empty_queue() -> None:
    output = subprocess.check_output(
        [
            TRUSTED_SQUEUE,
            "-u",
            pwd.getpwuid(os.getuid()).pw_name,
            "-h",
            "-o",
            "%i",
        ],
        text=True,
        env=trusted_slurm_environment(),
    )
    if output.strip():
        raise ValueError("Manual accounting requires an empty trusted Frontier queue")


def _scheduler_results(
    authorized_jobs: list[dict[str, object]],
) -> dict[str, dict[str, object]]:
    job_ids = [str(job["job_id"]) for job in authorized_jobs]
    output = subprocess.check_output(
        [
            TRUSTED_SACCT,
            "-j",
            ",".join(job_ids),
            "--allocations",
            f"--clusters={AUTHORIZED_SLURM_CLUSTER}",
            "--noheader",
            "--parsable2",
            "--format=JobIDRaw,State,ElapsedRaw,AllocNodes,Comment,Account,Partition,QOS",
        ],
        text=True,
        env=trusted_slurm_environment(),
    )
    expected_qos = {
        str(job["job_id"]): str(job["expected_qos"]) for job in authorized_jobs
    }
    rows: dict[str, dict[str, object]] = {}
    for line in output.splitlines():
        fields = line.split("|")
        if len(fields) != 8 or fields[0] not in expected_qos:
            raise ValueError("Slurm accounting returned an unexpected allocation record")
        job_id, raw_state, raw_elapsed, raw_nodes, comment, account, partition, qos = fields
        if job_id in rows:
            raise ValueError(f"Expected one exact Slurm allocation record for job {job_id}")
        state = raw_state.split()[0].upper()
        if state not in TERMINAL_STATES:
            raise ValueError(f"Job is not in a recognized terminal Slurm state: {state}")
        try:
            elapsed_seconds = int(raw_elapsed)
            allocated_nodes = int(raw_nodes)
        except ValueError as error:
            raise ValueError("Scheduler accounting values must be integers") from error
        if elapsed_seconds < 0 or allocated_nodes < 0:
            raise ValueError("Scheduler accounting values must not be negative")
        if comment:
            raise ValueError("Manual direct-srun allocation must have an empty Slurm comment")
        if not scheduler_account_matches_authorized(account):
            raise ValueError("Slurm accounting account does not match the PIC account")
        if partition != AUTHORIZED_PARTITION:
            raise ValueError("Slurm accounting partition does not match the PIC partition")
        if qos != expected_qos[job_id]:
            raise ValueError("Slurm accounting QoS does not match the reviewed authorization")
        rows[job_id] = {
            "job_id": job_id,
            "state": state,
            "elapsed_seconds": elapsed_seconds,
            "allocated_nodes": allocated_nodes,
            "partition": partition,
            "qos": qos,
        }
    if set(rows) != set(job_ids):
        raise ValueError("Expected one exact Slurm allocation record for each authorized job")
    return rows


def _event_payload(
    scheduler: dict[str, object],
    *,
    authorization_id: str,
    authorization_path: Path,
    project_home_authorization_path: Path,
    authorization_sha256: str,
    control_plane_version: str,
    active_policy_sha256: str,
    active_promotion_sha256: str,
    cumulative: float,
) -> dict[str, object]:
    allocated_nodes = int(scheduler["allocated_nodes"])
    elapsed_seconds = int(scheduler["elapsed_seconds"])
    consumed = allocated_nodes * elapsed_seconds / 3600.0
    return {
        "event_type": "manual_allocation_reconciliation",
        "job_id": scheduler["job_id"],
        "control_plane_version": control_plane_version,
        "reconciled_by_control_plane_version": control_plane_version,
        "manual_accounting_authorization_id": authorization_id,
        "manual_accounting_authorization_path": str(authorization_path),
        "manual_accounting_project_home_authorization_path": str(
            project_home_authorization_path
        ),
        "manual_accounting_authorization_sha256": authorization_sha256,
        "accounting_scope": MANUAL_ACCOUNTING_SCOPE,
        "scientific_evidence_eligible": False,
        "active_policy_sha256": active_policy_sha256,
        "active_promotion_sha256": active_promotion_sha256,
        "partition": scheduler["partition"],
        "qos": scheduler["qos"],
        "scheduler_reported_allocated_nodes": allocated_nodes,
        "billed_nodes": allocated_nodes,
        "elapsed_seconds": elapsed_seconds,
        "consumed_node_hours": consumed,
        "cumulative_consumed_node_hours": cumulative + consumed,
        "state": scheduler["state"],
        "reconciled": True,
        "notes": "Reviewed direct-srun accounting only; ineligible for scientific evidence.",
    }


def _stable_event_payload(record: dict[str, object]) -> dict[str, object]:
    omitted = {
        "event_sha256",
        "previous_event_sha256",
        "sequence_number",
        "timestamp",
    }
    return {key: value for key, value in record.items() if key not in omitted}


def reconcile_manual_allocations(
    *,
    authorization: Path,
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> list[dict[str, object]]:
    require_ledger_paths(
        ledger_jsonl,
        ledger_csv,
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    reviewed, authorization_sha256, authorization, project_home_authorization = (
        _authorization(
            authorization,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
    )
    _verify_installed_control_plane_pair(
        control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    with ledger_lock(
        ledger_jsonl,
        mirror_jsonl,
        allow_incomplete_manual_accounting=True,
    ):
        inventory = _verify_installed_control_plane_pair(
            control_plane_dir,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        version = str(inventory["version"])
        policy, policy_snapshot = require_storage_policy_unlock_snapshot(
            control_plane_version=version,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        policy_binding = {
            "authorization_id": reviewed["authorization_id"],
            "path": str(authorization),
            "project_home_path": str(project_home_authorization),
            "sha256": authorization_sha256,
        }
        storage = policy["olcf_side_storage"]
        assert isinstance(storage, dict)
        if policy_binding not in storage["manual_accounting_authorizations"]:
            raise ValueError(
                "Manual-accounting authorization is not bound by the promoted policy"
            )
        jobs = reviewed["jobs"]
        assert isinstance(jobs, list)
        job_ids = [str(job["job_id"]) for job in jobs]
        authorization_id = str(reviewed["authorization_id"])
        local_records = validate_primary_chain(ledger_jsonl)
        require_explicit_genesis(local_records)
        pending_marker = (
            Path(os.path.abspath(authorized_pic_root)) / "ledger" / "pending_submission.json"
        )
        require_canonical_path_below(pending_marker, authorized_pic_root)
        if pending_marker.exists():
            raise ValueError("Manual accounting is blocked by a pending scheduler submission")
        totals = accounting(local_records)
        outstanding = [
            record for record in latest_reservations(local_records).values()
            if record.get("state") in {"reserved", "submitted"}
        ]
        if outstanding or totals["currently_reserved_node_hours"]:
            raise ValueError("Manual accounting is blocked by an active reservation")
        _require_empty_queue()
        scheduler = _scheduler_results(jobs)
        incomplete_marker = recover_incomplete_manual_accounting_locked(
            ledger_jsonl,
            receipts_jsonl,
            mirror_jsonl,
            authorization_id=authorization_id,
            authorization_sha256=authorization_sha256,
            authorization_path=authorization,
            project_home_authorization_path=project_home_authorization,
            reviewed_job_ids=job_ids,
            reviewed_scheduler_results=scheduler,
            control_plane_version=version,
            active_policy_sha256=policy_snapshot["active_policy_sha256"],
            active_promotion_sha256=policy_snapshot["active_promotion_sha256"],
        )
        records = validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
        require_explicit_genesis(records)
        existing = [
            record for record in records
            if record.get("event_type") == "manual_allocation_reconciliation"
            and record.get("manual_accounting_authorization_id") == authorization_id
        ]
        if [record.get("job_id") for record in existing] != job_ids[:len(existing)]:
            raise ValueError("Existing manual-accounting events are not an authorized prefix")
        if len(existing) > len(job_ids):
            raise ValueError("Existing manual-accounting events exceed the authorization")
        existing_positions = [
            index for index, record in enumerate(records) if record in existing
        ]
        if existing_positions and existing_positions != list(
            range(existing_positions[0], existing_positions[0] + len(existing))
        ):
            raise ValueError("Existing manual-accounting events are not contiguous")
        if (
            existing
            and len(existing) < len(job_ids)
            and existing_positions[-1] != len(records) - 1
        ):
            raise ValueError(
                "Partial manual-accounting tranche is not the terminal ledger suffix"
            )
        existing_by_job = {str(record["job_id"]): record for record in existing}
        if len(existing_by_job) != len(existing):
            raise ValueError("Existing manual-accounting events contain duplicate jobs")
        prior_by_job: dict[str, list[dict[str, object]]] = {
            job_id: [
                record for record in records if record.get("job_id") == job_id
            ]
            for job_id in job_ids
        }
        for job_id, prior in prior_by_job.items():
            if prior and prior != ([existing_by_job[job_id]] if job_id in existing_by_job else []):
                raise ValueError(f"Authorized job already has a prior ledger event: {job_id}")
        preceding_records = (
            records[:existing_positions[0]] if existing_positions else records
        )
        cumulative = accounting(preceding_records)["cumulative_consumed_node_hours"]
        if incomplete_marker is None:
            incomplete_marker = {
                "schema_version": 3,
                "state": "manual_accounting_incomplete",
                "manual_accounting_authorization_id": authorization_id,
                "manual_accounting_authorization_sha256": authorization_sha256,
                "pre_tranche_sequence_number": len(records),
                "pre_tranche_chain_head": chain_head(records),
                "pre_tranche_authorized_job_count": len(existing),
                "control_plane_version": version,
                "active_policy_sha256": policy_snapshot["active_policy_sha256"],
                "active_promotion_sha256": policy_snapshot["active_promotion_sha256"],
            }
            marker_published = False
        else:
            marker_published = True
        results: list[dict[str, object]] = []
        for index, job_id in enumerate(job_ids):
            historical = existing[index] if index < len(existing) else None
            payload = _event_payload(
                scheduler[job_id],
                authorization_id=authorization_id,
                authorization_path=authorization,
                project_home_authorization_path=project_home_authorization,
                authorization_sha256=authorization_sha256,
                control_plane_version=(
                    str(historical["control_plane_version"])
                    if historical is not None
                    else version
                ),
                active_policy_sha256=(
                    str(historical["active_policy_sha256"])
                    if historical is not None
                    else policy_snapshot["active_policy_sha256"]
                ),
                active_promotion_sha256=(
                    str(historical["active_promotion_sha256"])
                    if historical is not None
                    else policy_snapshot["active_promotion_sha256"]
                ),
                cumulative=cumulative,
            )
            if historical is not None:
                result = historical
                if _stable_event_payload(result) != payload:
                    raise ValueError(
                        f"Existing manual-accounting event differs from reviewed scheduler data: {job_id}"
                    )
            else:
                if not marker_published:
                    publish_incomplete_manual_accounting_marker_locked(
                        ledger_jsonl, mirror_jsonl, incomplete_marker
                    )
                    marker_published = True
                result = append_primary_event_locked(
                    ledger_jsonl,
                    ledger_csv,
                    receipts_jsonl,
                    mirror_jsonl,
                    payload,
                    mirror_transport="filesystem_copy",
                    allow_incomplete_manual_accounting=True,
                )
            cumulative = float(payload["cumulative_consumed_node_hours"])
            results.append(result)
        if not marker_published:
            publish_incomplete_manual_accounting_marker_locked(
                ledger_jsonl, mirror_jsonl, incomplete_marker
            )
        write_csv(
            ledger_jsonl,
            receipts_jsonl,
            ledger_csv,
            mirror_jsonl=mirror_jsonl,
            mirror_transport="filesystem_copy",
            allow_incomplete_manual_accounting=True,
        )
        clear_matching_incomplete_manual_accounting_marker_locked(
            ledger_jsonl, mirror_jsonl, incomplete_marker
        )
        return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--authorization", required=True, type=Path)
    parser.add_argument("--ledger-jsonl", required=True, type=Path)
    parser.add_argument("--ledger-csv", required=True, type=Path)
    parser.add_argument("--receipts-jsonl", required=True, type=Path)
    parser.add_argument("--mirror-jsonl", required=True, type=Path)
    args = parser.parse_args()
    events = reconcile_manual_allocations(
        authorization=args.authorization,
        ledger_jsonl=args.ledger_jsonl,
        ledger_csv=args.ledger_csv,
        receipts_jsonl=args.receipts_jsonl,
        mirror_jsonl=args.mirror_jsonl,
    )
    print(events[-1]["cumulative_consumed_node_hours"])


if __name__ == "__main__":
    main()
