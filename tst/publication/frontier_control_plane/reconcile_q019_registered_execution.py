#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Reconcile Q019 and publish deterministic registered-execution evidence."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re

from control_plane_common import AUTHORIZED_PIC_ROOT, AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import launch_contract_sha256, project_home_ledger_root
from control_plane_common import record_for_role, require_canonical_path_below
from control_plane_common import validate_launch_contract
from control_plane_common import verify_installed_control_plane
from ledger import require_explicit_genesis, validate_mirrored_state
from reconcile_frontier_job import reconcile
import reconcile_q043_registered_execution as secure
from validate_and_reserve_frontier_job import _require_run_artifact_dir


SCRIPT_DIR = Path(__file__).absolute().parent
ENTRYPOINT_NAME = "reconcile_q019_registered_execution.py"
TRAMPOLINE_ENTRYPOINT = "launch_trampoline.py"
REGISTERED_CAMPAIGN = "q019_nonlinear_bell_registered_successor_v1"
CAMPAIGN_ID = "Q019-PHYSICS-FIRST-NONLINEAR-BELL"
AUTHORIZATION_PREFIX = "q019-"
EXECUTION_RECEIPT_NAME = "q019_registered_execution_receipt.json"
EXECUTION_RECEIPT_RECORD_TYPE = "q019_reconciled_registered_execution_receipt"
TERMINAL_RECEIPT_NAME = "q019_terminal_receipt.json"
TERMINAL_RECEIPT_RECORD_TYPE = "q019_registered_execution_terminal_receipt"
PROJECT_HOME_MIRROR_NAMESPACE = Path("ledger/q019_registered_execution_receipts")
REQUIRED_FIELDS = (
    "mhd_w_bcc",
    "prtcl_rho",
    "prtcl_jx",
    "prtcl_jy",
    "prtcl_jz",
    "prtcl_dedt",
    "prtcl_dpxdt",
    "prtcl_dpydt",
    "prtcl_dpzdt",
    "prtcl_ebdot",
)
_SHA256 = re.compile(r"[0-9a-f]{64}")
_UUID = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}"
)
_JOB_ID = re.compile(r"[1-9][0-9]*")
_MEMBER_ID = re.compile(r"q019-[a-z0-9][a-z0-9_-]{0,123}")
_OUTPUT_INDEX = r"([0-9]{5,})"


def _json_bytes(value: dict[str, object]) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _canonical_sha256(value: object) -> str:
    payload = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _inventory_file_digest(inventory: dict[str, object], name: str) -> str:
    records = inventory.get("files")
    if not isinstance(records, list):
        raise ValueError("Installed control-plane inventory files are malformed")
    matches = [
        record
        for record in records
        if isinstance(record, dict) and record.get("path") == name
    ]
    if len(matches) != 1:
        raise ValueError(f"Installed control plane lacks required entrypoint: {name}")
    digest = matches[0].get("sha256")
    if not isinstance(digest, str) or _SHA256.fullmatch(digest) is None:
        raise ValueError(f"Installed control-plane digest is malformed: {name}")
    return digest


def _producer_binding(inventory: dict[str, object]) -> dict[str, object]:
    version = inventory.get("version")
    if not isinstance(version, str) or _SHA256.fullmatch(version) is None:
        raise ValueError("Q019 producer control-plane version is malformed")
    return {
        "entrypoint": ENTRYPOINT_NAME,
        "entrypoint_sha256": _inventory_file_digest(inventory, ENTRYPOINT_NAME),
        "launch_trampoline_sha256": _inventory_file_digest(
            inventory, TRAMPOLINE_ENTRYPOINT
        ),
        "control_plane_version": version,
    }


def _require_bootstrapped_installed_producer(
    control_plane_dir: Path,
    inventory: dict[str, object],
    *,
    authorized_pic_root: Path,
) -> None:
    if not getattr(_sys, "_pic_control_plane_bootstrapped", False):
        raise ValueError("Q019 reconciliation requires bootstrapped installed execution")
    captured = globals().get("_PIC_CAPTURED_CONTROL_PLANE_BINDING")
    producer = _producer_binding(inventory)
    if (
        not isinstance(captured, dict)
        or set(captured) != {"directory", "version", "filename", "sha256"}
        or captured.get("directory") != str(Path(os.path.abspath(control_plane_dir)))
        or captured.get("directory") != str(Path(os.path.abspath(SCRIPT_DIR)))
        or captured.get("version") != inventory.get("version")
        or captured.get("filename") != ENTRYPOINT_NAME
        or captured.get("sha256") != producer["entrypoint_sha256"]
    ):
        raise ValueError(
            "Q019 captured executing entrypoint differs from verified generation"
        )
    require_canonical_path_below(
        Path(os.path.abspath(control_plane_dir)),
        Path(os.path.abspath(authorized_pic_root)),
    )


def _project_home_mirror_paths(
    submission_id: str, *, authorized_project_home_root: Path
) -> tuple[Path, Path]:
    root = secure._canonical_project_home_root(authorized_project_home_root)
    mirror_root = root / PROJECT_HOME_MIRROR_NAMESPACE / submission_id
    return (
        mirror_root / TERMINAL_RECEIPT_NAME,
        mirror_root / EXECUTION_RECEIPT_NAME,
    )


def _q019_event(
    records: list[dict[str, object]],
    *,
    job_id: str,
    producer_control_plane_version: str,
) -> dict[str, object]:
    matches = [
        record
        for record in records
        if record.get("event_type") == "reconciliation"
        and record.get("job_id") == job_id
        and record.get("campaign") == REGISTERED_CAMPAIGN
    ]
    if len(matches) != 1:
        raise ValueError("Expected one Q019 reconciliation event for the scheduler job")
    event = matches[0]
    authorization = event.get("registered_science_authorization_id")
    if (
        event.get("submission_scope") != "registered_science"
        or not isinstance(authorization, str)
        or not authorization.startswith(AUTHORIZATION_PREFIX)
        or event.get("state") != "COMPLETED"
        or event.get("scheduler_exit_code") != "0:0"
        or event.get("reconciled") is not True
        or event.get("reconciled_by_control_plane_version")
        != producer_control_plane_version
        or not isinstance(event.get("event_sha256"), str)
        or _SHA256.fullmatch(str(event["event_sha256"])) is None
    ):
        raise ValueError("Q019 reconciliation event is not a successful registered execution")
    member_id = event.get("test_id")
    if not isinstance(member_id, str) or _MEMBER_ID.fullmatch(member_id) is None:
        raise ValueError("Q019 reconciliation member identity is malformed")
    return event


def _trusted_wrapper_evidence(
    payload: bytes, *, member_id: str, tasks: int, stdout_sha256: str
) -> dict[str, object]:
    try:
        lines = payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("Q019 retained stdout is not UTF-8") from error
    rank_line = (
        f"Q019_REGISTERED_EXECUTION case_id={member_id} "
        f"mpi_world_size={tasks} rank_ids={','.join(str(rank) for rank in range(tasks))}"
    )
    exit_line = "Q019_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0"
    if lines.count(rank_line) != 1 or lines.count(exit_line) != 1:
        raise ValueError("Q019 retained stdout lacks exact trusted rank/exit evidence")
    finite = r"[+-]?(?:[0-9]+(?:[.][0-9]*)?|[.][0-9]+)(?:[eE][+-]?[0-9]+)?"
    terminal = [
        (float(match.group(1)), int(match.group(2)))
        for line in lines
        if (match := re.fullmatch(rf"time=({finite}) cycle=([0-9]+)", line))
        is not None
    ]
    limits = [
        (float(match.group(1)), int(match.group(2)))
        for line in lines
        if (match := re.fullmatch(rf"tlim=({finite}) nlim=([0-9]+)", line))
        is not None
    ]
    termination_lines = {
        "Terminating on time limit",
        "Terminating on cycle limit",
        "Terminating on user request",
    }
    observed_termination = [line for line in lines if line in termination_lines]
    if (
        len(observed_termination) != 1
        or len(terminal) != 1
        or len(limits) != 1
        or not math.isfinite(terminal[0][0])
        or terminal[0][0] < 0.0
        or terminal[0][1] <= 0
        or not math.isfinite(limits[0][0])
        or limits[0][0] <= 0.0
        or limits[0][1] <= 0
        or terminal[0][0] > limits[0][0] + 1.0e-12
    ):
        raise ValueError("Q019 retained stdout lacks exact terminal-time evidence")
    status_lines = [
        line
        for line in lines
        if line.startswith("Q019_FINAL_EVIDENCE_STATUS=")
    ]
    eligibility_lines = [
        line
        for line in lines
        if line.startswith("Q019_SATURATION_EVIDENCE_ELIGIBLE=")
    ]
    if len(status_lines) != 1 or len(eligibility_lines) != 1:
        raise ValueError("Q019 retained stdout lacks exact final-evidence status")
    return {
        "source": "installed_trampoline_retained_stdout_exact_bytes",
        "stdout_sha256": stdout_sha256,
        "required_exact_rank_line": rank_line,
        "required_exact_exit_line": exit_line,
        "observed_world_size": tasks,
        "observed_rank_ids": list(range(tasks)),
        "exit_code": 0,
        "signal": 0,
        "terminal_time": terminal[0][0],
        "terminal_cycle": terminal[0][1],
        "configured_tlim": limits[0][0],
        "configured_nlim": limits[0][1],
        "termination_reason": observed_termination[0],
        "problem_final_evidence_status": status_lines[0].split("=", 1)[1],
        "problem_saturation_evidence_eligible": eligibility_lines[0].split(
            "=", 1
        )[1],
    }


def _contiguous_indices(values: set[int], *, label: str) -> tuple[int, ...]:
    if not values or min(values) != 0 or values != set(range(max(values) + 1)):
        raise ValueError(f"Q019 {label} indices are not contiguous from zero")
    return tuple(sorted(values))


def _structured_launch_evidence(
    manifest: dict[str, object],
    sealed_inventory: dict[str, object],
    *,
    member_id: str,
    artifact_dir: Path,
    producer: dict[str, object],
) -> dict[str, object]:
    contract = validate_launch_contract(manifest.get("launch_contract"))
    actions = contract["actions"]
    if len(actions) != 1:
        raise ValueError("Q019 requires one trusted trampoline Athena action")
    action = actions[0]
    arguments = action["arguments"]
    raw_directories = [
        str(arguments[index + 1]["artifact_directory"])
        for index, argument in enumerate(arguments[:-1])
        if argument == {"literal": "-d"}
        and isinstance(arguments[index + 1], dict)
        and set(arguments[index + 1]) == {"artifact_directory"}
    ]
    input_decks = sum(
        1
        for index, argument in enumerate(arguments[:-1])
        if argument == {"literal": "-i"}
        and arguments[index + 1] == {"snapshot_role": "input-deck"}
    )
    if len(raw_directories) != 1 or input_decks != 1:
        raise ValueError("Q019 launch contract lacks one raw directory and input deck")
    raw_relative = secure._safe_relative(
        raw_directories[0], label="Q019 raw artifact directory"
    )
    raw_root = artifact_dir.joinpath(*PurePosixPath(raw_relative).parts)
    tasks = action["resources"]["tasks"]
    if type(tasks) is not int or not 1 <= tasks <= 4096:
        raise ValueError("Q019 launch contract has unsupported MPI task count")
    records = sealed_inventory["records"]
    if not isinstance(records, dict):
        raise ValueError("Q019 sealed artifact records are malformed")
    stdout_record = records.get("athena_stdout.txt")
    if not isinstance(stdout_record, dict):
        raise ValueError("Q019 sealed artifact records lack trusted stdout")
    wrapper = _trusted_wrapper_evidence(
        sealed_inventory["athena_stdout_payload"],
        member_id=member_id,
        tasks=tasks,
        stdout_sha256=str(stdout_record["sha256"]),
    )
    basename = member_id
    prefix = raw_relative + "/"
    raw_paths = sorted(
        relative[len(prefix) :]
        for relative in records
        if relative.startswith(prefix)
    )
    binary_pattern = re.compile(
        rf"bin/{re.escape(basename)}[.]({'|'.join(REQUIRED_FIELDS)})[.]"
        rf"{_OUTPUT_INDEX}[.]bin"
    )
    vtk_pattern = re.compile(
        rf"pvtk/{re.escape(basename)}[.]prtcl_all[.]"
        rf"{_OUTPUT_INDEX}[.]part[.]vtk"
    )
    restart_pattern = re.compile(
        rf"rst/{re.escape(basename)}[.]"
        rf"{_OUTPUT_INDEX}[.]rst(?:[.]complete|[.]manifest(?:[.]complete)?)?"
    )
    histories = {
        f"{basename}.mhd.hst",
        f"{basename}.user.hst",
    }
    binary_by_field = {field: set() for field in REQUIRED_FIELDS}
    vtk_indices: set[int] = set()
    restart_indices: set[int] = set()
    observed_histories: set[str] = set()
    unsupported: list[str] = []
    for relative in raw_paths:
        if match := binary_pattern.fullmatch(relative):
            binary_by_field[match.group(1)].add(int(match.group(2)))
        elif match := vtk_pattern.fullmatch(relative):
            vtk_indices.add(int(match.group(1)))
        elif match := restart_pattern.fullmatch(relative):
            restart_indices.add(int(match.group(1)))
        elif relative in histories:
            observed_histories.add(relative)
        else:
            unsupported.append(relative)
    if unsupported:
        raise ValueError("Q019 sealed raw tree contains unsupported output")
    binary_indices = _contiguous_indices(
        binary_by_field[REQUIRED_FIELDS[0]], label="binary snapshot"
    )
    if len(binary_indices) < 2 or any(
        binary_by_field[field] != set(binary_indices)
        for field in REQUIRED_FIELDS
    ):
        raise ValueError("Q019 ten-product binary chronology is incomplete")
    checkpoint_indices = _contiguous_indices(
        restart_indices, label="restart checkpoint"
    )
    if (
        len(checkpoint_indices) < 2
        or vtk_indices != set(checkpoint_indices)
        or observed_histories != histories
    ):
        raise ValueError("Q019 checkpoint, particle-VTK, or history inventory drifted")
    expected_paths = set(histories)
    for index in binary_indices:
        expected_paths.update(
            f"bin/{basename}.{field}.{index:05d}.bin"
            for field in REQUIRED_FIELDS
        )
    for index in checkpoint_indices:
        expected_paths.add(f"pvtk/{basename}.prtcl_all.{index:05d}.part.vtk")
        stem = f"rst/{basename}.{index:05d}.rst"
        expected_paths.update(
            {stem, stem + ".complete", stem + ".manifest", stem + ".manifest.complete"}
        )
    if set(raw_paths) != expected_paths:
        raise ValueError("Q019 sealed raw output differs from exact expected inventory")
    raw_inventory = [
        {
            "path": relative,
            "sha256": records[prefix + relative]["sha256"],
            "byte_count": records[prefix + relative]["size"],
            "member_id": member_id,
            "artifact_kind": (
                "binary_snapshot"
                if relative.startswith("bin/")
                else "particle_vtk"
                if relative.startswith("pvtk/")
                else "restart_publication"
                if relative.startswith("rst/")
                else "history"
            ),
        }
        for relative in raw_paths
    ]
    return {
        "raw_root": raw_root,
        "raw_inventory": raw_inventory,
        "raw_inventory_sha256": _canonical_sha256(raw_inventory),
        "command_evidence": {
            "source": "trusted_pre_submit_manifest_and_installed_trampoline",
            "executor": contract["executor"],
            "action": action,
            "launch_contract_sha256": launch_contract_sha256(contract),
            "launch_trampoline_entrypoint": TRAMPOLINE_ENTRYPOINT,
            "launch_trampoline_sha256": producer["launch_trampoline_sha256"],
            "trusted_wrapper_evidence": wrapper,
        },
        "mpi_evidence": {
            "source": "trusted_pre_submit_manifest_and_installed_trampoline_stdout",
            **action["resources"],
            "observed_world_size": wrapper["observed_world_size"],
            "observed_rank_ids": wrapper["observed_rank_ids"],
        },
        "terminal_time": wrapper["terminal_time"],
        "terminal_cycle": wrapper["terminal_cycle"],
        "binary_output_indices": list(binary_indices),
        "checkpoint_output_indices": list(checkpoint_indices),
    }


def derive_q019_registered_execution_evidence(
    event: dict[str, object],
    mirror_ack: dict[str, object],
    inventory: dict[str, object],
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> tuple[Path, bytes, Path, bytes, Path, Path]:
    """Derive exact Q019 receipts from controller-owned immutable evidence."""
    root = Path(os.path.abspath(authorized_pic_root))
    producer = _producer_binding(inventory)
    if (
        event.get("reconciled_by_control_plane_version")
        != producer["control_plane_version"]
        or mirror_ack.get("mirrored_event_sha256") != event.get("event_sha256")
    ):
        raise ValueError("Q019 producer, reconciliation, and mirror acknowledgment differ")
    manifest_path = require_canonical_path_below(
        Path(str(event.get("manifest_path", ""))), root / "manifests"
    )
    manifest, manifest_payload = secure._read_json(
        manifest_path, root=root, label="Q019 pre-submit manifest"
    )
    if hashlib.sha256(manifest_payload).hexdigest() != event.get("manifest_sha256"):
        raise ValueError("Q019 pre-submit manifest digest differs from reconciliation")
    artifact_dir = _require_run_artifact_dir(manifest)
    if (
        manifest.get("submission_scope") != "registered_science"
        or manifest.get("campaign") != REGISTERED_CAMPAIGN
        or manifest.get("submission_id") != event.get("submission_id")
        or manifest.get("test_id") != event.get("test_id")
        or manifest.get("control_plane_version") != event.get("control_plane_version")
        or manifest.get("git_commit") != event.get("git_commit")
        or manifest.get("registered_science_authorization_id")
        != event.get("registered_science_authorization_id")
        or event.get("artifact_dir") != str(artifact_dir)
    ):
        raise ValueError("Q019 pre-submit manifest differs from reconciliation")
    member_id = str(event["test_id"])
    submission_id = str(event["submission_id"])
    job_id = str(event["job_id"])
    if _UUID.fullmatch(submission_id) is None or _JOB_ID.fullmatch(job_id) is None:
        raise ValueError("Q019 reconciliation has malformed scheduler identity")
    trampoline_completion = secure._trusted_trampoline_completion(
        manifest,
        event,
        authorized_pic_root=root,
        authorized_project_home_root=authorized_project_home_root,
    )
    sealed_inventory = secure._sealed_artifact_inventory(
        artifact_dir,
        authorized_pic_root=root,
        trampoline_completion=trampoline_completion,
    )
    launch = _structured_launch_evidence(
        manifest,
        sealed_inventory,
        member_id=member_id,
        artifact_dir=artifact_dir,
        producer=producer,
    )
    terminal_mirror_path, receipt_mirror_path = _project_home_mirror_paths(
        submission_id,
        authorized_project_home_root=authorized_project_home_root,
    )
    terminal = {
        "schema_version": 1,
        "record_type": TERMINAL_RECEIPT_RECORD_TYPE,
        "campaign_id": CAMPAIGN_ID,
        "member_id": member_id,
        "submission_id": submission_id,
        "slurm_job_id": job_id,
        "slurm_terminal_state": event["state"],
        "slurm_exit_code": event["scheduler_exit_code"],
        "terminal_cycle": launch["terminal_cycle"],
        "terminal_time": launch["terminal_time"],
        "binary_output_indices": launch["binary_output_indices"],
        "checkpoint_output_indices": launch["checkpoint_output_indices"],
        "registered_mpi_tasks": launch["mpi_evidence"]["tasks"],
        "artifact_inventory_sha256": sealed_inventory["sha256"],
        "trampoline_completion_receipt_sha256": trampoline_completion["sha256"],
        "raw_inventory_sha256": launch["raw_inventory_sha256"],
        "reconciliation_event_sha256": event["event_sha256"],
        "reconciliation_mirror_ack_sha256": mirror_ack["mirror_ack_sha256"],
        "project_home_mirror_path": str(terminal_mirror_path),
        "producer": producer,
    }
    terminal_path = artifact_dir / "analysis" / TERMINAL_RECEIPT_NAME
    terminal_payload = _json_bytes(terminal)
    candidate_manifest_path = Path(str(manifest["clean_candidate_manifest_path"]))
    candidate, candidate_payload = secure._read_json(
        candidate_manifest_path, root=root, label="Q019 clean-candidate manifest"
    )
    source = candidate.get("source")
    if (
        hashlib.sha256(candidate_payload).hexdigest()
        != manifest.get("clean_candidate_manifest_sha256")
        or not isinstance(source, dict)
        or source.get("git_commit") != event.get("git_commit")
    ):
        raise ValueError("Q019 clean candidate differs from pre-submit manifest")
    receipt = {
        "schema_version": 1,
        "record_type": EXECUTION_RECEIPT_RECORD_TYPE,
        "receipt_role": "immutable_reconciled_registered_execution",
        "registration_scope": "registered_science",
        "reconciled": True,
        "campaign_id": CAMPAIGN_ID,
        "member_id": member_id,
        "reservation_id": event["reservation_id"],
        "submission_id": submission_id,
        "reconciliation_event_sha256": event["event_sha256"],
        "reconciliation_mirror_ack_sha256": mirror_ack["mirror_ack_sha256"],
        "control_plane_version": event["control_plane_version"],
        "project_home_mirrors": {
            "registered_execution_receipt_path": str(receipt_mirror_path),
            "terminal_receipt_path": str(terminal_mirror_path),
        },
        "producer": producer,
        "registered_science_authorization_id": event[
            "registered_science_authorization_id"
        ],
        "source_commit": event["git_commit"],
        "source_bundle_sha256": source["source_bundle_sha256"],
        "source_archive_sha256": source["archive_sha256"],
        "clean_candidate_manifest_sha256": manifest[
            "clean_candidate_manifest_sha256"
        ],
        "executable_sha256": record_for_role(manifest, "executable")["sha256"],
        "environment_sha256": record_for_role(manifest, "environment-profile")[
            "sha256"
        ],
        "deck_sha256": record_for_role(manifest, "input-deck")["sha256"],
        "command_evidence": launch["command_evidence"],
        "mpi_evidence": launch["mpi_evidence"],
        "slurm_job_id": job_id,
        "slurm_terminal_state": event["state"],
        "slurm_exit_code": event["scheduler_exit_code"],
        "terminal_cycle": launch["terminal_cycle"],
        "terminal_time": launch["terminal_time"],
        "binary_output_indices": launch["binary_output_indices"],
        "checkpoint_output_indices": launch["checkpoint_output_indices"],
        "raw_output_root": str(launch["raw_root"]),
        "artifact_dir": str(artifact_dir),
        "artifact_inventory": {
            key: sealed_inventory[key] for key in ("path", "sha256", "byte_count")
        },
        "trampoline_completion_receipt": {
            key: trampoline_completion[key]
            for key in ("orion_path", "project_home_path", "sha256", "byte_count")
        },
        "terminal_receipt_sha256": hashlib.sha256(terminal_payload).hexdigest(),
        "pre_submit_manifest_path": str(manifest_path),
        "pre_submit_manifest_sha256": event["manifest_sha256"],
        "raw_inventory": launch["raw_inventory"],
        "raw_inventory_sha256": launch["raw_inventory_sha256"],
    }
    receipt_path = artifact_dir / "analysis" / EXECUTION_RECEIPT_NAME
    return (
        terminal_path,
        terminal_payload,
        receipt_path,
        _json_bytes(receipt),
        terminal_mirror_path,
        receipt_mirror_path,
    )


def publish_q019_registered_execution_evidence(
    event: dict[str, object],
    mirror_ack: dict[str, object],
    inventory: dict[str, object],
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    derived = derive_q019_registered_execution_evidence(
        event,
        mirror_ack,
        inventory,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    paths = (derived[0], derived[2], derived[4], derived[5])
    payloads = (derived[1], derived[3], derived[1], derived[3])
    roots = (
        authorized_pic_root,
        authorized_pic_root,
        authorized_project_home_root,
        authorized_project_home_root,
    )
    identities = [
        secure._publish_exact(path, payload, authorized_root=root)
        for path, payload, root in zip(paths, payloads, roots)
    ]
    if identities[0] == identities[2] or identities[1] == identities[3]:
        raise ValueError("Q019 Orion and Project Home evidence reuse one filesystem object")
    os.chmod(paths[3].parent, 0o500, follow_symlinks=False)
    receipt = json.loads(derived[3])
    return {
        "terminal_receipt": {
            "orion_path": str(paths[0]),
            "project_home_mirror_path": str(paths[2]),
            "sha256": hashlib.sha256(derived[1]).hexdigest(),
        },
        "registered_execution_receipt": {
            "orion_path": str(paths[1]),
            "project_home_mirror_path": str(paths[3]),
            "sha256": hashlib.sha256(derived[3]).hexdigest(),
        },
        "artifact_inventory_sha256": receipt["artifact_inventory"]["sha256"],
        "raw_inventory_sha256": receipt["raw_inventory_sha256"],
        "reconciliation_mirror_ack_sha256": mirror_ack["mirror_ack_sha256"],
        "producer": _producer_binding(inventory),
    }


def reconcile_q019(
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
    inventory = verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    _require_bootstrapped_installed_producer(
        control_plane_dir, inventory, authorized_pic_root=authorized_pic_root
    )
    paired = verify_installed_control_plane(
        Path(os.path.abspath(authorized_project_home_root))
        / "control_plane"
        / str(inventory["version"]),
        authorized_pic_root=authorized_project_home_root,
    )
    if paired != inventory:
        raise ValueError("Q019 producer installed control-plane pair differs")
    before = validate_mirrored_state(
        ledger_jsonl,
        receipts_jsonl,
        mirror_jsonl,
        ledger_root=authorized_pic_root,
        receipts_root=authorized_pic_root,
        mirror_root=project_home_ledger_root(authorized_project_home_root),
    )
    require_explicit_genesis(before)
    completion_events = [
        record
        for record in before
        if record.get("event_type") == "trampoline_completion"
        and record.get("job_id") == job_id
        and record.get("campaign") == REGISTERED_CAMPAIGN
    ]
    if len(completion_events) != 1:
        raise ValueError(
            "Q019 reconciliation requires one canonical trampoline-completion anchor"
        )
    completion_event = completion_events[0]
    secure._mirror_ack(
        completion_event,
        before,
        receipts_jsonl=receipts_jsonl,
        mirror_jsonl=mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
    )
    completion_manifest_path = require_canonical_path_below(
        Path(str(completion_event.get("manifest_path", ""))),
        Path(os.path.abspath(authorized_pic_root)) / "manifests",
    )
    completion_manifest, completion_payload = secure._read_json(
        completion_manifest_path,
        root=Path(os.path.abspath(authorized_pic_root)),
        label="Q019 trampoline-completion pre-submit manifest",
    )
    if (
        hashlib.sha256(completion_payload).hexdigest()
        != completion_event.get("manifest_sha256")
    ):
        raise ValueError("Q019 trampoline-completion manifest digest differs")
    secure._trusted_trampoline_completion(
        completion_manifest,
        completion_event,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    event = reconcile(
        job_id=job_id,
        ledger_jsonl=ledger_jsonl,
        ledger_csv=ledger_csv,
        receipts_jsonl=receipts_jsonl,
        mirror_jsonl=mirror_jsonl,
        control_plane_dir=control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    records = validate_mirrored_state(
        ledger_jsonl,
        receipts_jsonl,
        mirror_jsonl,
        ledger_root=authorized_pic_root,
        receipts_root=authorized_pic_root,
        mirror_root=project_home_ledger_root(authorized_project_home_root),
    )
    require_explicit_genesis(records)
    trusted_event = _q019_event(
        records,
        job_id=job_id,
        producer_control_plane_version=str(inventory["version"]),
    )
    if trusted_event != event:
        raise ValueError("Q019 reconciler result differs from canonical mirrored ledger")
    mirror_ack = secure._mirror_ack(
        event,
        records,
        receipts_jsonl=receipts_jsonl,
        mirror_jsonl=mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
    )
    evidence = publish_q019_registered_execution_evidence(
        event,
        mirror_ack,
        inventory,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    return {"reconciliation": event, "evidence": evidence}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--ledger-jsonl", required=True, type=Path)
    parser.add_argument("--ledger-csv", required=True, type=Path)
    parser.add_argument("--receipts-jsonl", required=True, type=Path)
    parser.add_argument("--mirror-jsonl", required=True, type=Path)
    args = parser.parse_args()
    result = reconcile_q019(
        job_id=args.job_id,
        ledger_jsonl=args.ledger_jsonl,
        ledger_csv=args.ledger_csv,
        receipts_jsonl=args.receipts_jsonl,
        mirror_jsonl=args.mirror_jsonl,
    )
    print(json.dumps(result["evidence"], sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
