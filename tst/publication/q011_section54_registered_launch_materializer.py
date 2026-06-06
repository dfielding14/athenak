#!/usr/bin/env python3
"""Materialize non-authorizing Q-011 registered-launch review candidates.

The tool consumes one closed selected-bindings record and exact immutable
planner-retention records.  It emits policy-slice and pre-submit-config review
candidates only.  It never edits live policy, captures operator attestations,
invokes a scheduler, promotes a controller, or submits work.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import shutil
import stat
import sys
from typing import Any, Mapping

CONTROL_PLANE_SOURCE_DIR = Path(__file__).resolve().parent / "frontier_control_plane"
if str(CONTROL_PLANE_SOURCE_DIR) not in sys.path:
    sys.path.insert(0, str(CONTROL_PLANE_SOURCE_DIR))
from control_plane_common import (  # type: ignore[import-not-found]
    AUTHORIZED_ACCOUNT,
    AUTHORIZED_NODE_HOUR_CAP,
    AUTHORIZED_PIC_ROOT,
    launch_contract_sha256,
    validate_launch_contract,
    validate_planner_retention_binding,
    verify_installed_control_plane,
)
from ledger import slurm_walltime_seconds  # type: ignore[import-not-found]


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_CONTRACT = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_registered_launch_orchestration_preregistration_2026-06-06.json"
)
SELECTED_BINDING_RECORD_TYPE = "q011_section54_registered_launch_selected_bindings"
POLICY_FRAGMENT_RECORD_TYPE = "q011_section54_registered_launch_policy_slice_fragment"
CONFIG_CANDIDATE_RECORD_TYPE = (
    "q011_section54_registered_launch_pre_submit_config_candidate"
)
RECEIPT_RECORD_TYPE = "q011_section54_registered_launch_materialization_receipt"
PHYSICAL_MODE = "paper_mhd_pic_vl2_tsc"
EVIDENCE_CLASS = "q011_section54_qualifying_production"
RUNTIME_PROFILE = "frontier_minimum_supported"
CAMPAIGN = "q011_section54_qualifying"
TRUSTED_EXECUTOR = "trusted_trampoline_athena_argv_v1"
REGISTERED_SLICE_KEYS = {
    "authorization_id",
    "status",
    "campaign",
    "test_id",
    "evidence_class",
    "physical_mode",
    "runtime_profile",
    "selected_qos",
    "registered_short_nonproduction",
    "maximum_nodes",
    "maximum_walltime_seconds",
    "maximum_attempts",
    "job_script_sha256",
    "input_deck_sha256",
    "environment_profile_sha256",
    "analysis_script_sha256",
    "executable_sha256",
    "launch_contract_sha256",
    "clean_candidate_manifest_sha256",
}
_SHA256 = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT = re.compile(r"[0-9a-f]{40}")
_ATTEMPT_ID = re.compile(r"baseline-(?P<index>[0-9]{3})-[a-z0-9._-]+")
_SAFE_AUTHORIZATION = re.compile(r"[a-z0-9][a-z0-9_-]{0,63}")
_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH


class RegisteredLaunchMaterializationError(ValueError):
    """Raised when registered-launch review materialization fails closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RegisteredLaunchMaterializationError(message)


def _json_bytes(value: object) -> bytes:
    try:
        return (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise RegisteredLaunchMaterializationError(
            "registered-launch artifact contains noncanonical JSON values"
        ) from error


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _decode_json(payload: bytes, *, label: str) -> Any:
    def reject_constant(value: str) -> None:
        raise RegisteredLaunchMaterializationError(
            f"{label} contains non-finite JSON value {value}"
        )

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            _require(key not in result, f"{label} contains duplicate key {key!r}")
            result[key] = value
        return result

    try:
        return json.loads(
            payload.decode("utf-8"),
            parse_constant=reject_constant,
            object_pairs_hook=reject_duplicates,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RegisteredLaunchMaterializationError(
            f"{label} is not valid UTF-8 JSON"
        ) from error


def _stable_regular_bytes(
    path: Path, *, label: str, require_read_only: bool = True
) -> tuple[Path, bytes]:
    lexical = Path(os.path.abspath(path))
    _require(lexical == path, f"{label} path must be absolute and normalized")
    _require(path.resolve() == path, f"{label} path must not contain symlinks")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise RegisteredLaunchMaterializationError(f"{label} cannot be opened") from error
    try:
        before = os.fstat(descriptor)
        _require(stat.S_ISREG(before.st_mode), f"{label} must be a regular file")
        if require_read_only:
            _require(not before.st_mode & _WRITE_BITS, f"{label} must be read-only")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            payload = stream.read()
        after = os.fstat(descriptor)
        stable_fields = (
            "st_dev",
            "st_ino",
            "st_mode",
            "st_size",
            "st_mtime_ns",
            "st_ctime_ns",
        )
        _require(
            all(getattr(before, field) == getattr(after, field) for field in stable_fields)
            and len(payload) == after.st_size,
            f"{label} changed while reading",
        )
        return lexical, payload
    finally:
        os.close(descriptor)


def _read_canonical_json(
    path: Path, *, label: str, require_read_only: bool = True
) -> tuple[Path, bytes, dict[str, Any]]:
    lexical, payload = _stable_regular_bytes(
        path, label=label, require_read_only=require_read_only
    )
    value = _decode_json(payload, label=label)
    _require(type(value) is dict, f"{label} must be one JSON object")
    _require(payload == _json_bytes(value), f"{label} must use canonical JSON bytes")
    return lexical, payload, value


def _object(value: object, keys: set[str], *, label: str) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == keys, f"{label} schema is malformed")
    return value


def _text(value: object, *, label: str) -> str:
    _require(type(value) is str and bool(value.strip()), f"{label} must be non-empty text")
    return value.strip()


def _digest(value: object, *, label: str) -> str:
    text = _text(value, label=label)
    _require(_SHA256.fullmatch(text) is not None, f"{label} must be lowercase SHA-256")
    return text


def _positive_int(value: object, *, label: str) -> int:
    _require(type(value) is int and value > 0, f"{label} must be a positive integer")
    return value


def _binding(
    value: object, *, label: str, require_read_only: bool = True
) -> tuple[dict[str, str], bytes]:
    record = _object(value, {"path", "sha256"}, label=label)
    path = Path(_text(record["path"], label=f"{label}/path"))
    digest = _digest(record["sha256"], label=f"{label}/sha256")
    lexical, payload = _stable_regular_bytes(
        path, label=label, require_read_only=require_read_only
    )
    _require(_sha256_bytes(payload) == digest, f"{label} SHA-256 drifted")
    return {"path": str(lexical), "sha256": digest}, payload


def _readiness_contract(
    path: Path, *, require_read_only: bool = True
) -> tuple[dict[str, str], dict[str, Any]]:
    lexical, payload, value = _read_canonical_json(
        path,
        label="registered-launch readiness contract",
        require_read_only=require_read_only,
    )
    expected = {
        "record_type": "q011_section54_registered_launch_orchestration_preregistration",
        "schema_version": 1,
        "date": "2026-06-06",
        "status": "source_local_review_materializer_only_no_launch_authority",
        "qualification_effect": (
            "none_review_candidates_only_no_launch_authority_no_claim_closure"
        ),
        "selected_binding_record_type": SELECTED_BINDING_RECORD_TYPE,
        "required_attempt_scope": (
            "one_or_more_unique_ordered_baseline_attempts_each_with_one_exact_"
            "revalidated_planner_retention_record"
        ),
        "resource_ceiling_policy": {
            "maximum_attempts_per_registered_slice": 1,
            "maximum_total_node_hours": 10000.0,
            "selected_qos": "normal",
            "qos_selection_reason": "normal_required_by_registered_campaign",
            "registered_short_nonproduction": False,
        },
        "registered_slice_required_fields": [
            "authorization_id",
            "status",
            "campaign",
            "test_id",
            "evidence_class",
            "physical_mode",
            "runtime_profile",
            "selected_qos",
            "registered_short_nonproduction",
            "maximum_nodes",
            "maximum_walltime_seconds",
            "maximum_attempts",
            "job_script_sha256",
            "input_deck_sha256",
            "environment_profile_sha256",
            "analysis_script_sha256",
            "executable_sha256",
            "launch_contract_sha256",
            "clean_candidate_manifest_sha256",
        ],
        "pre_submit_candidate_late_bindings": [
            "fresh_submission_uuid",
            "fresh_pre_manifest_operator_attestation_after_policy_promotion",
            "fresh_timeout_margin_artifact",
            "fresh_six_field_queue_snapshot",
            "fresh_site_policy_checked_utc",
            "artifact_dir_derived_from_campaign_and_submission_uuid",
        ],
        "outputs": [
            "registered_policy_slice_fragment.json",
            "pre_submit_config_candidates/<attempt_id>.json",
            "materialization_receipt.json",
        ],
        "non_authorization_guards": {
            "claim_closure_authorized": False,
            "frontier_execution_authorized": False,
            "live_policy_mutation_authorized": False,
            "scheduler_calls_authorized": False,
            "scheduler_submission_authorized": False,
        },
        "required_next_gates": [
            "review_the_registered_policy_slice_fragment",
            "freeze_source_and_install_the_exact_selected_controller_generation",
            "authorize_the_exact_selected_clean_candidate_in_a_launch_prohibited_policy",
            "verify_remaining_ledger_budget_against_the_declared_batch_ceiling",
            "capture_a_fresh_pre_policy_promotion_operator_attestation_with_an_empty_user_queue",
            "promote_a_separately_reviewed_policy_successor_containing_only_the_selected_slices",
            "for_each_attempt_capture_fresh_late_bindings_and_materialize_one_exact_pre_submit_config",
            "capture_a_fresh_pre_submit_wrapper_attestation_and_use_only_the_installed_submission_wrapper",
            "wait_for_terminal_state_and_reconcile_before_materializing_another_attempt",
        ],
    }
    _require(value == expected, "registered-launch readiness contract drifted")
    return {"path": str(lexical), "sha256": _sha256_bytes(payload)}, value


def _walltime_seconds(value: str) -> int:
    try:
        result = slurm_walltime_seconds(value)
    except ValueError as error:
        raise RegisteredLaunchMaterializationError(
            "selected job-script walltime is malformed"
        ) from error
    _require(result > 0, "selected job-script walltime must be positive")
    return result


def _job_directives(payload: bytes) -> dict[str, str]:
    try:
        lines = payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise RegisteredLaunchMaterializationError(
            "selected job script is not UTF-8"
        ) from error
    aliases = {
        "-A": "account",
        "--account": "account",
        "-p": "partition",
        "--partition": "partition",
        "-q": "qos",
        "--qos": "qos",
        "-N": "nodes",
        "--nodes": "nodes",
        "-t": "time",
        "--time": "time",
        "-o": "output",
        "--output": "output",
    }
    options: dict[str, str] = {}
    for line in lines:
        if not line.startswith("#SBATCH"):
            continue
        remainder = line[len("#SBATCH") :].strip()
        option, separator, attached = remainder.partition("=")
        if separator:
            key = aliases.get(option)
            value = attached
        else:
            try:
                tokens = shlex.split(remainder)
            except ValueError as error:
                raise RegisteredLaunchMaterializationError(
                    "selected job script has malformed Slurm quoting"
                ) from error
            key = aliases.get(tokens[0]) if tokens else None
            value = tokens[1] if len(tokens) == 2 else ""
        if key is None:
            continue
        _require(bool(value), f"selected job script has missing {key} directive value")
        _require(key not in options, f"selected job script has duplicate {key} directive")
        options[key] = value
    required = {"account", "partition", "qos", "nodes", "time", "output"}
    _require(set(options) == required, "selected job script required directive schema drifted")
    return options


def _plan_summary(
    planner_root: Path, *, inventory_sha256: str, plan_id: str
) -> dict[str, Any]:
    _require(planner_root.is_absolute(), "selected planner root must be absolute")
    _require(planner_root.resolve() == planner_root, "selected planner root must not contain symlinks")
    _require(
        planner_root.name == f"q011-section54-qualifying-campaign-plan-{plan_id}",
        "selected planner root name differs from selected plan ID",
    )
    _, inventory_payload = _stable_regular_bytes(
        planner_root / "artifact_inventory.sha256",
        label="immutable planner inventory",
    )
    _require(
        _sha256_bytes(inventory_payload) == inventory_sha256,
        "selected planner inventory SHA-256 drifted",
    )
    try:
        text = inventory_payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise RegisteredLaunchMaterializationError(
            "immutable planner inventory is not UTF-8"
        ) from error
    inventory: dict[str, str] = {}
    for line in text.splitlines(keepends=True):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)\n", line)
        _require(match is not None, "immutable planner inventory is malformed")
        relative = PurePosixPath(match.group(2))
        _require(
            not relative.is_absolute()
            and relative.parts
            and all(part not in {"", ".", ".."} for part in relative.parts),
            "immutable planner inventory member path is unsafe",
        )
        key = relative.as_posix()
        _require(key not in inventory, "immutable planner inventory contains duplicates")
        inventory[key] = match.group(1)
    _require(
        text == "".join(f"{inventory[path]}  {path}\n" for path in sorted(inventory)),
        "immutable planner inventory is noncanonical",
    )
    _require("campaign_plan.json" in inventory, "immutable planner omits campaign plan")
    _, plan_payload = _stable_regular_bytes(
        planner_root / "campaign_plan.json", label="immutable planner campaign plan"
    )
    _require(
        _sha256_bytes(plan_payload) == inventory["campaign_plan.json"],
        "immutable planner campaign-plan SHA-256 drifted",
    )
    plan = _decode_json(plan_payload, label="immutable planner campaign plan")
    _require(type(plan) is dict, "immutable planner campaign plan must be an object")
    _require(
        plan.get("record_type") == "q011_section54_qualifying_campaign_execution_plan"
        and plan.get("schema_version") == 1
        and plan.get("plan_id") == plan_id
        and plan.get("status") == "source_local_immutable_review_plan_only"
        and plan.get("execution_boundary")
        == {
            "mutates_live_policy": False,
            "scheduler_calls": False,
            "submits_jobs": False,
            "infers_pressure_selection": False,
            "launch_authorized": False,
            "frontier_execution_authorized": False,
            "claim_closure_authorized": False,
        },
        "immutable planner campaign-plan identity drifted",
    )
    return plan


def _clean_candidate(
    value: object, *, authorized_pic_root: Path
) -> tuple[dict[str, Any], dict[str, str], dict[str, str], str]:
    selected = _object(
        value, {"manifest", "executable", "git_commit"}, label="selected clean candidate"
    )
    manifest_binding, manifest_payload = _binding(
        selected["manifest"], label="selected clean-candidate manifest"
    )
    executable_binding, _ = _binding(
        selected["executable"], label="selected clean-candidate executable"
    )
    git_commit = _text(selected["git_commit"], label="selected clean-candidate Git commit")
    _require(_GIT_COMMIT.fullmatch(git_commit) is not None, "selected clean-candidate Git commit is malformed")
    manifest_path = Path(manifest_binding["path"])
    candidate_root = authorized_pic_root / "clean_candidates"
    _require(
        manifest_path.parent.parent == candidate_root
        and manifest_path.name == "clean_candidate_manifest.json",
        "selected clean-candidate manifest path is not one level below clean_candidates",
    )
    _require(
        Path(executable_binding["path"]) == manifest_path.parent / "athena",
        "selected clean-candidate executable is not adjacent to its manifest",
    )
    manifest = _decode_json(manifest_payload, label="selected clean-candidate manifest")
    _require(type(manifest) is dict, "selected clean-candidate manifest must be an object")
    source = manifest.get("source")
    build = manifest.get("build")
    _require(type(source) is dict and type(build) is dict, "selected clean-candidate manifest objects drifted")
    _require(source.get("git_commit") == git_commit, "selected clean-candidate Git commit drifted")
    _require(
        build.get("executable_path") == executable_binding["path"]
        and build.get("executable_sha256") == executable_binding["sha256"],
        "selected clean-candidate executable binding drifted",
    )
    return manifest, manifest_binding, executable_binding, git_commit


def _resource_profile(value: object, *, attempt_count: int) -> dict[str, Any]:
    profile = _object(
        value,
        {
            "selected_qos",
            "qos_selection_reason",
            "registered_short_nonproduction",
            "maximum_nodes",
            "maximum_walltime_seconds",
            "maximum_attempts_per_slice",
            "maximum_batch_node_hours",
            "launch_resources",
        },
        label="selected resource profile",
    )
    _require(profile["selected_qos"] == "normal", "Q011 qualifying production requires selected_qos=normal")
    _require(
        profile["qos_selection_reason"] == "normal_required_by_registered_campaign",
        "Q011 qualifying production QoS reason drifted",
    )
    _require(
        profile["registered_short_nonproduction"] is False,
        "Q011 qualifying production cannot be registered short nonproduction",
    )
    maximum_nodes = _positive_int(profile["maximum_nodes"], label="resource maximum_nodes")
    maximum_walltime = _positive_int(
        profile["maximum_walltime_seconds"], label="resource maximum_walltime_seconds"
    )
    _require(
        profile["maximum_attempts_per_slice"] == 1
        and type(profile["maximum_attempts_per_slice"]) is int,
        "each Q011 registered slice must have maximum_attempts_per_slice=1",
    )
    resources = _object(
        profile["launch_resources"],
        {"nodes", "tasks", "cpus_per_task", "gpus_per_task", "gpu_bind"},
        label="selected launch resources",
    )
    for field in ["nodes", "tasks", "cpus_per_task", "gpus_per_task"]:
        _positive_int(resources[field], label=f"selected launch resources/{field}")
    _require(resources["gpu_bind"] == "closest", "selected launch resources gpu_bind must be closest")
    _require(resources["nodes"] <= maximum_nodes, "launch resources exceed selected node ceiling")
    declared = profile["maximum_batch_node_hours"]
    _require(
        isinstance(declared, (int, float))
        and not isinstance(declared, bool)
        and math.isfinite(float(declared))
        and float(declared) > 0,
        "selected maximum_batch_node_hours is invalid",
    )
    computed = attempt_count * maximum_nodes * maximum_walltime / 3600.0
    _require(
        float(declared) == computed,
        "selected maximum_batch_node_hours differs from attempt/resource ceilings",
    )
    _require(
        computed <= float(AUTHORIZED_NODE_HOUR_CAP),
        "selected batch ceiling exceeds the authorized total node-hour cap",
    )
    return profile


def _attempt_identity(attempt_id: str, *, plan_id: str) -> dict[str, str]:
    match = _ATTEMPT_ID.fullmatch(attempt_id)
    _require(match is not None, "selected attempt ID is not one baseline planner attempt")
    index = match.group("index")
    authorization_id = f"q011-s54-{index}-{plan_id[:12]}-v1"
    _require(
        _SAFE_AUTHORIZATION.fullmatch(authorization_id) is not None,
        "derived registered-science authorization ID is malformed",
    )
    return {
        "authorization_id": authorization_id,
        "campaign": CAMPAIGN,
        "test_id": f"q011_s54_baseline_{index}",
    }


def _runtime_launch_contract(
    retention: Mapping[str, Any], *, resources: Mapping[str, Any]
) -> dict[str, Any]:
    argv = retention.get("argv")
    _require(
        type(argv) is list
        and len(argv) >= 5
        and all(type(value) is str and value for value in argv),
        "planner-retention argv is malformed",
    )
    _require(
        argv[:4]
        == [
            "-i",
            "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
            "-d",
            retention["authorized_orion_raw_root"],
        ],
        "planner-retention argv input/output prefix drifted",
    )
    attempt_id = str(retention["attempt_id"])
    index = _ATTEMPT_ID.fullmatch(attempt_id)
    _require(index is not None, "planner-retention attempt ID is malformed")
    contract = {
        "schema_version": 1,
        "executor": TRUSTED_EXECUTOR,
        "pre_actions": [],
        "actions": [
            {
                "action_id": f"q011-s54-{index.group('index')}",
                "kind": "athena",
                "resources": dict(resources),
                "arguments": [
                    {"literal": "-i"},
                    {"snapshot_role": "input-deck"},
                    {"literal": "-d"},
                    {"artifact_directory": "raw"},
                    *({"literal": token} for token in argv[4:]),
                ],
                "stdout_artifact": "athena_stdout.txt",
                "stderr_artifact": "athena_stderr.txt",
            }
        ],
        "post_actions": [
            {
                "action_id": "require-athena-stdout",
                "kind": "artifact_nonempty",
                "artifact": "athena_stdout.txt",
            },
            {
                "action_id": "sha-athena-stdout",
                "kind": "artifact_sha256",
                "artifact": "athena_stdout.txt",
                "output_artifact": "athena_stdout.sha256",
            },
        ],
    }
    try:
        return validate_launch_contract(contract)
    except ValueError as error:
        raise RegisteredLaunchMaterializationError(
            f"derived runtime launch contract is invalid: {error}"
        ) from error


def _registered_slice(
    *,
    identity: Mapping[str, str],
    resource_profile: Mapping[str, Any],
    runtime_artifacts: Mapping[str, Any],
    environment_profile: Mapping[str, str],
    executable: Mapping[str, str],
    clean_manifest: Mapping[str, str],
    launch_contract: Mapping[str, Any],
) -> dict[str, Any]:
    result = {
        "authorization_id": identity["authorization_id"],
        "status": "authorized",
        "campaign": identity["campaign"],
        "test_id": identity["test_id"],
        "evidence_class": EVIDENCE_CLASS,
        "physical_mode": PHYSICAL_MODE,
        "runtime_profile": RUNTIME_PROFILE,
        "selected_qos": resource_profile["selected_qos"],
        "registered_short_nonproduction": resource_profile[
            "registered_short_nonproduction"
        ],
        "maximum_nodes": resource_profile["maximum_nodes"],
        "maximum_walltime_seconds": resource_profile["maximum_walltime_seconds"],
        "maximum_attempts": resource_profile["maximum_attempts_per_slice"],
        "job_script_sha256": runtime_artifacts["job_script"]["sha256"],
        "input_deck_sha256": runtime_artifacts["input_deck"]["sha256"],
        "environment_profile_sha256": environment_profile["sha256"],
        "analysis_script_sha256": [
            record["sha256"] for record in runtime_artifacts["analysis_scripts"]
        ],
        "executable_sha256": executable["sha256"],
        "launch_contract_sha256": launch_contract_sha256(launch_contract),
        "clean_candidate_manifest_sha256": clean_manifest["sha256"],
    }
    _require(set(result) == REGISTERED_SLICE_KEYS, "derived registered-science slice schema drifted")
    return result


def _pre_submit_candidate(
    *,
    attempt_id: str,
    identity: Mapping[str, str],
    policy_slice: Mapping[str, Any],
    retention: Mapping[str, Any],
    resource_profile: Mapping[str, Any],
    runtime_artifacts: Mapping[str, Any],
    clean_manifest: Mapping[str, str],
    executable: Mapping[str, str],
    controller: Mapping[str, Any],
    environment_profile: Mapping[str, str],
    git_commit: str,
    launch_contract: Mapping[str, Any],
    authorized_pic_root: Path,
) -> dict[str, Any]:
    static = {
        "pic_root": str(authorized_pic_root),
        "campaign": identity["campaign"],
        "test_id": identity["test_id"],
        "submission_scope": "registered_science",
        "registered_science_authorization_id": identity["authorization_id"],
        "git_commit": git_commit,
        "evidence_class": EVIDENCE_CLASS,
        "physical_mode": PHYSICAL_MODE,
        "selected_qos": resource_profile["selected_qos"],
        "qos_selection_reason": resource_profile["qos_selection_reason"],
        "registered_short_nonproduction": resource_profile[
            "registered_short_nonproduction"
        ],
        "job_script_executable_env": "PIC_EXECUTABLE",
        "job_script": runtime_artifacts["job_script"]["path"],
        "executable": executable["path"],
        "input_deck": runtime_artifacts["input_deck"]["path"],
        "environment_profile": environment_profile["path"],
        "analysis_scripts": [
            record["path"] for record in runtime_artifacts["analysis_scripts"]
        ],
        "prior_case_closures": [],
        "clean_candidate_manifest": clean_manifest["path"],
        "launch_contract": dict(launch_contract),
        "planner_retention": dict(retention),
    }
    return {
        "record_type": CONFIG_CANDIDATE_RECORD_TYPE,
        "schema_version": 1,
        "status": "review_candidate_incomplete_not_pre_submit_config",
        "qualification_effect": "none_no_launch_authority_no_claim_closure",
        "attempt_id": attempt_id,
        "authorization_id": identity["authorization_id"],
        "registered_policy_slice_sha256": _sha256_bytes(_json_bytes(policy_slice)),
        "selected_controller_version": controller["version"],
        "declared_resource_ceiling": {
            "maximum_nodes": resource_profile["maximum_nodes"],
            "maximum_walltime_seconds": resource_profile["maximum_walltime_seconds"],
            "maximum_attempts": resource_profile["maximum_attempts_per_slice"],
            "maximum_node_hours": (
                resource_profile["maximum_nodes"]
                * resource_profile["maximum_walltime_seconds"]
                / 3600.0
            ),
        },
        "static_pre_submit_config": static,
        "required_fresh_late_bindings": {
            "submission_id": "fresh_uuid_required_at_final_config_materialization",
            "pre_manifest_attestation": {
                "phase": "pre_manifest",
                "authorization_id": identity["authorization_id"],
                "control_plane_version": controller["version"],
                "capture_boundary": "after_exact_registered_policy_promotion",
            },
            "timeout_margin_artifact": {
                "scheduler_walltime_seconds": resource_profile[
                    "maximum_walltime_seconds"
                ],
                "environment_profile_sha256": environment_profile["sha256"],
                "athena_walltime_rule": "positive_integer_strictly_below_scheduler_walltime",
                "freshness_rule": "must_be_valid_at_reservation_time",
            },
            "queue_snapshot": {
                "format": "%i|%P|%q|%T|%j|%k",
                "freshness_rule": "must_equal_fresh_reservation_boundary_queue_output",
            },
            "site_policy_checked_utc": "fresh_canonical_utc_required",
            "artifact_dir": {
                "parent": str(authorized_pic_root / "runs" / identity["campaign"]),
                "rule": "parent/fresh_submission_uuid",
            },
        },
        "execution_boundary": {
            "complete_pre_submit_config": False,
            "live_policy_mutation_authorized": False,
            "scheduler_calls_authorized": False,
            "scheduler_submission_authorized": False,
            "frontier_execution_authorized": False,
            "claim_closure_authorized": False,
        },
    }


def _write_new(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags, 0o444)
    try:
        offset = 0
        while offset < len(payload):
            offset += os.write(descriptor, payload[offset:])
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _freeze_tree(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o555 if path.is_dir() else 0o444)
    root.chmod(0o555)


def _protected_live_output(output_root: Path, authorized_pic_root: Path) -> bool:
    try:
        relative = output_root.relative_to(authorized_pic_root)
    except ValueError:
        return False
    return not relative.parts or relative.parts[0] != "review_candidates"


def materialize_registered_launch_review_bundle(
    *,
    selected_bindings: Path,
    output_root: Path,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    readiness_contract: Path = READINESS_CONTRACT,
) -> dict[str, Any]:
    """Write one recursively read-only, non-authorizing review bundle."""
    authorized_pic_root = Path(os.path.abspath(authorized_pic_root))
    _require(
        authorized_pic_root.is_dir() and authorized_pic_root.resolve() == authorized_pic_root,
        "authorized PIC root must be one existing canonical directory",
    )
    output_root = Path(os.path.abspath(output_root))
    _require(
        output_root.parent.is_dir() and output_root.parent.resolve() == output_root.parent,
        "review-bundle output parent must be one existing canonical directory",
    )
    _require(
        not _protected_live_output(output_root, authorized_pic_root),
        "review bundle cannot be materialized in a live PIC control namespace",
    )
    readiness_binding, readiness = _readiness_contract(
        readiness_contract,
        # Git does not preserve non-executable permission bits.  The exact
        # built-in source contract is instead protected by the full schema
        # equality check below; caller-supplied contracts remain read-only.
        require_read_only=readiness_contract != READINESS_CONTRACT,
    )
    selected_path, selected_payload, selected = _read_canonical_json(
        selected_bindings, label="selected registered-launch bindings"
    )
    _object(
        selected,
        {
            "record_type",
            "schema_version",
            "authorized_pic_root",
            "planner",
            "clean_candidate",
            "controller",
            "runtime_artifacts",
            "resource_profile",
            "attempts",
        },
        label="selected registered-launch bindings",
    )
    _require(
        selected["record_type"] == SELECTED_BINDING_RECORD_TYPE
        and selected["schema_version"] == 1
        and type(selected["schema_version"]) is int,
        "selected registered-launch binding identity drifted",
    )
    _require(
        selected["authorized_pic_root"] == str(authorized_pic_root),
        "selected registered-launch authorized PIC root drifted",
    )
    planner = _object(
        selected["planner"],
        {"root", "inventory_sha256", "plan_id"},
        label="selected planner binding",
    )
    planner_root = Path(_text(planner["root"], label="selected planner root"))
    planner_inventory = _digest(
        planner["inventory_sha256"], label="selected planner inventory SHA-256"
    )
    plan_id = _digest(planner["plan_id"], label="selected planner plan ID")
    _require(
        planner_root.parent == authorized_pic_root / "plans",
        "selected planner root is outside the authorized plans root",
    )
    attempts = selected["attempts"]
    _require(type(attempts) is list and attempts, "selected attempts must be a non-empty list")
    normalized_attempts: list[dict[str, Any]] = []
    prior_attempt_id = ""
    seen_indices: set[str] = set()
    for value in attempts:
        attempt = _object(
            value, {"attempt_id", "planner_retention"}, label="selected attempt binding"
        )
        attempt_id = _text(attempt["attempt_id"], label="selected attempt ID")
        match = _ATTEMPT_ID.fullmatch(attempt_id)
        _require(match is not None, "selected attempt ID is not a baseline attempt")
        _require(attempt_id > prior_attempt_id, "selected attempt bindings must be uniquely ordered")
        _require(match.group("index") not in seen_indices, "selected attempt baseline index is duplicated")
        prior_attempt_id = attempt_id
        seen_indices.add(match.group("index"))
        retention_binding, retention_payload = _binding(
            attempt["planner_retention"], label=f"{attempt_id} planner-retention record"
        )
        retention = _decode_json(
            retention_payload, label=f"{attempt_id} planner-retention record"
        )
        _require(type(retention) is dict, f"{attempt_id} planner-retention record must be an object")
        _require(retention.get("attempt_id") == attempt_id, f"{attempt_id} planner-retention attempt binding drifted")
        _require(retention.get("planner_root") == str(planner_root), f"{attempt_id} planner-retention planner root drifted")
        _require(retention.get("planner_inventory_sha256") == planner_inventory, f"{attempt_id} planner-retention inventory drifted")
        _require(retention.get("planner_plan_id") == plan_id, f"{attempt_id} planner-retention plan ID drifted")
        normalized_attempts.append(
            {
                "attempt_id": attempt_id,
                "planner_retention_binding": retention_binding,
                "planner_retention": retention,
            }
        )
    resource_profile = _resource_profile(
        selected["resource_profile"], attempt_count=len(normalized_attempts)
    )
    _, clean_manifest, executable, git_commit = _clean_candidate(
        selected["clean_candidate"], authorized_pic_root=authorized_pic_root
    )
    controller = _object(
        selected["controller"],
        {"path", "version", "environment_profile"},
        label="selected controller binding",
    )
    controller_path = Path(_text(controller["path"], label="selected controller path"))
    controller_version = _digest(controller["version"], label="selected controller version")
    _require(
        controller_path == authorized_pic_root / "control_plane" / controller_version,
        "selected controller path/version binding drifted",
    )
    try:
        inventory = verify_installed_control_plane(
            controller_path, authorized_pic_root=authorized_pic_root
        )
    except ValueError as error:
        raise RegisteredLaunchMaterializationError(
            f"selected installed controller is invalid: {error}"
        ) from error
    _require(inventory.get("version") == controller_version, "selected installed controller version drifted")
    environment_profile, _ = _binding(
        controller["environment_profile"], label="selected controller environment profile"
    )
    _require(
        Path(environment_profile["path"]) == controller_path / "frontier_pic_environment.sh",
        "selected environment profile is not the installed controller profile",
    )
    runtime = _object(
        selected["runtime_artifacts"],
        {"job_script", "input_deck", "analysis_scripts"},
        label="selected runtime artifacts",
    )
    job_script, job_payload = _binding(runtime["job_script"], label="selected job script")
    input_deck, _ = _binding(runtime["input_deck"], label="selected input deck")
    analysis = runtime["analysis_scripts"]
    _require(type(analysis) is list and 1 <= len(analysis) <= 16, "selected analysis scripts must contain one to sixteen records")
    analysis_bindings = [
        _binding(value, label=f"selected analysis script {index:03d}")[0]
        for index, value in enumerate(analysis)
    ]
    _require(
        len({record["path"] for record in analysis_bindings}) == len(analysis_bindings)
        and len({record["sha256"] for record in analysis_bindings}) == len(analysis_bindings),
        "selected analysis scripts contain duplicate paths or bytes",
    )
    runtime_artifacts = {
        "job_script": job_script,
        "input_deck": input_deck,
        "analysis_scripts": analysis_bindings,
    }
    directives = _job_directives(job_payload)
    _require(directives["account"] == AUTHORIZED_ACCOUNT, "selected job script account is not authorized")
    _require(directives["partition"] == "batch", "selected job script partition is not batch")
    _require(directives["qos"] == resource_profile["selected_qos"], "selected job script QoS differs from resource profile")
    try:
        directive_nodes = int(directives["nodes"])
    except ValueError as error:
        raise RegisteredLaunchMaterializationError(
            "selected job script nodes directive is malformed"
        ) from error
    _require(
        directive_nodes > 0
        and directive_nodes == resource_profile["maximum_nodes"],
        "selected job script nodes differ from resource ceiling",
    )
    _require(
        _walltime_seconds(directives["time"])
        == resource_profile["maximum_walltime_seconds"],
        "selected job script walltime differs from resource ceiling",
    )
    _require(
        directives["output"] == str(authorized_pic_root / "logs/slurm/%x.%j.log"),
        "selected job script output path is not the dedicated PIC log path",
    )
    plan = _plan_summary(
        planner_root, inventory_sha256=planner_inventory, plan_id=plan_id
    )
    candidate = plan.get("candidate_binding")
    source_bindings = plan.get("source_bindings")
    _require(type(candidate) is dict and type(source_bindings) is dict, "immutable planner candidate/source bindings drifted")
    candidate_environment = candidate.get("environment_profile")
    _require(
        type(candidate_environment) is dict,
        "immutable planner environment-profile binding drifted",
    )
    _require(
        candidate.get("clean_candidate_manifest") == clean_manifest
        and candidate.get("executable") == executable
        and candidate_environment
        == {
            "path": environment_profile["path"],
            "sha256": environment_profile["sha256"],
            "control_plane_version": controller_version,
            "reviewed_source": candidate_environment.get("reviewed_source"),
        }
        and candidate.get("git_commit") == git_commit,
        "selected candidate/controller bindings differ from immutable planner",
    )
    paper_deck = source_bindings.get("paper_deck")
    _require(
        type(paper_deck) is dict and paper_deck.get("sha256") == input_deck["sha256"],
        "selected input deck differs from immutable planner paper deck",
    )
    validated_attempts: list[dict[str, Any]] = []
    for attempt in normalized_attempts:
        try:
            validated = validate_planner_retention_binding(
                attempt["planner_retention"],
                authorized_pic_root=authorized_pic_root,
                expected_clean_candidate_manifest_sha256=clean_manifest["sha256"],
            )
        except ValueError as error:
            raise RegisteredLaunchMaterializationError(
                f"{attempt['attempt_id']} planner retention is invalid: {error}"
            ) from error
        _require(
            validated == attempt["planner_retention"],
            f"{attempt['attempt_id']} planner retention validator returned drifted data",
        )
        validated_attempts.append(attempt)
    policy_slices: list[dict[str, Any]] = []
    config_candidates: list[tuple[str, dict[str, Any]]] = []
    attempt_records: list[dict[str, Any]] = []
    for attempt in validated_attempts:
        attempt_id = attempt["attempt_id"]
        identity = _attempt_identity(attempt_id, plan_id=plan_id)
        launch_contract = _runtime_launch_contract(
            attempt["planner_retention"],
            resources=resource_profile["launch_resources"],
        )
        policy_slice = _registered_slice(
            identity=identity,
            resource_profile=resource_profile,
            runtime_artifacts=runtime_artifacts,
            environment_profile=environment_profile,
            executable=executable,
            clean_manifest=clean_manifest,
            launch_contract=launch_contract,
        )
        candidate_record = _pre_submit_candidate(
            attempt_id=attempt_id,
            identity=identity,
            policy_slice=policy_slice,
            retention=attempt["planner_retention"],
            resource_profile=resource_profile,
            runtime_artifacts=runtime_artifacts,
            clean_manifest=clean_manifest,
            executable=executable,
            controller={"path": str(controller_path), "version": controller_version},
            environment_profile=environment_profile,
            git_commit=git_commit,
            launch_contract=launch_contract,
            authorized_pic_root=authorized_pic_root,
        )
        policy_slices.append(policy_slice)
        config_candidates.append((attempt_id, candidate_record))
        attempt_records.append(
            {
                "attempt_id": attempt_id,
                "authorization_id": identity["authorization_id"],
                "planner_retention": attempt["planner_retention_binding"],
                "launch_contract_sha256": policy_slice["launch_contract_sha256"],
                "maximum_attempts": 1,
                "maximum_nodes": resource_profile["maximum_nodes"],
                "maximum_walltime_seconds": resource_profile[
                    "maximum_walltime_seconds"
                ],
                "maximum_node_hours": (
                    resource_profile["maximum_nodes"]
                    * resource_profile["maximum_walltime_seconds"]
                    / 3600.0
                ),
            }
        )
    selected_binding = {"path": str(selected_path), "sha256": _sha256_bytes(selected_payload)}
    fragment = {
        "record_type": POLICY_FRAGMENT_RECORD_TYPE,
        "schema_version": 1,
        "status": "review_fragment_only_not_live_policy",
        "qualification_effect": readiness["qualification_effect"],
        "selected_bindings": selected_binding,
        "readiness_contract": readiness_binding,
        "planner": {
            "root": str(planner_root),
            "inventory_sha256": planner_inventory,
            "plan_id": plan_id,
        },
        "selected_clean_candidate": clean_manifest,
        "selected_controller": {
            "path": str(controller_path),
            "version": controller_version,
            "environment_profile": environment_profile,
        },
        "declared_batch_ceiling": {
            "included_attempts": len(validated_attempts),
            "maximum_registered_attempts": len(validated_attempts),
            "maximum_nodes_per_attempt": resource_profile["maximum_nodes"],
            "maximum_walltime_seconds_per_attempt": resource_profile[
                "maximum_walltime_seconds"
            ],
            "maximum_batch_node_hours": resource_profile["maximum_batch_node_hours"],
            "authorized_total_node_hour_cap": float(AUTHORIZED_NODE_HOUR_CAP),
        },
        "attempts": attempt_records,
        "registered_science_slices": policy_slices,
        "slice_status_semantics": (
            "status_authorized_is_required_policy_schema_only_and_has_no_effect_until_"
            "a_separately_reviewed_full_policy_successor_is_promoted"
        ),
        "required_next_gates": readiness["required_next_gates"],
        "execution_boundary": {
            "complete_storage_policy": False,
            "live_policy_mutation_authorized": False,
            "scheduler_calls_authorized": False,
            "scheduler_submission_authorized": False,
            "frontier_execution_authorized": False,
            "claim_closure_authorized": False,
        },
    }
    _require(not output_root.exists(), "review-bundle output root already exists")
    try:
        output_root.mkdir(mode=0o700)
        fragment_payload = _json_bytes(fragment)
        _write_new(output_root / "registered_policy_slice_fragment.json", fragment_payload)
        candidate_bindings: list[dict[str, str]] = []
        for attempt_id, candidate_record in config_candidates:
            relative = f"pre_submit_config_candidates/{attempt_id}.json"
            payload = _json_bytes(candidate_record)
            _write_new(output_root / relative, payload)
            candidate_bindings.append({"path": relative, "sha256": _sha256_bytes(payload)})
        receipt = {
            "record_type": RECEIPT_RECORD_TYPE,
            "schema_version": 1,
            "status": "materialized_review_candidates_only_no_launch_authority",
            "qualification_effect": readiness["qualification_effect"],
            "selected_bindings": selected_binding,
            "readiness_contract": readiness_binding,
            "registered_policy_slice_fragment": {
                "path": "registered_policy_slice_fragment.json",
                "sha256": _sha256_bytes(fragment_payload),
            },
            "pre_submit_config_candidates": candidate_bindings,
            "included_attempt_count": len(validated_attempts),
            "declared_batch_ceiling": fragment["declared_batch_ceiling"],
            "execution_boundary": fragment["execution_boundary"],
        }
        receipt_payload = _json_bytes(receipt)
        _write_new(output_root / "materialization_receipt.json", receipt_payload)
        _freeze_tree(output_root)
        return {
            "output_root": str(output_root),
            "materialization_receipt_sha256": _sha256_bytes(receipt_payload),
            "registered_policy_slice_fragment_sha256": _sha256_bytes(fragment_payload),
            "included_attempt_count": len(validated_attempts),
            "maximum_batch_node_hours": resource_profile["maximum_batch_node_hours"],
            "recursively_read_only": all(
                not path.stat().st_mode & _WRITE_BITS
                for path in [output_root, *output_root.rglob("*")]
            ),
            "launch_authorized": False,
        }
    except BaseException:
        if output_root.exists():
            for path in [output_root, *output_root.rglob("*")]:
                if not path.is_symlink():
                    path.chmod(path.stat().st_mode | stat.S_IWUSR)
            shutil.rmtree(output_root)
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected-bindings", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    arguments = parser.parse_args()
    result = materialize_registered_launch_review_bundle(
        selected_bindings=arguments.selected_bindings,
        output_root=arguments.output_root,
    )
    print(_json_bytes(result).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
