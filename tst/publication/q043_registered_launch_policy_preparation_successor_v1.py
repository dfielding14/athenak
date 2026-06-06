#!/usr/bin/env python3
"""Prepare non-authorizing Q043 registered-launch and policy review candidates.

This source-local successor binds the exact checked-in 132-case deposited-current
matrix to deterministic launch-contract, policy-slice, budget, output-inventory,
and trusted-reconciliation review candidates.  It never mutates live policy,
calls a scheduler, consumes a final clean candidate, or grants authority.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import sys
from typing import Any, Mapping, Sequence

from tst.publication import (
    q043_bell_current_volume_aware_deposited_current_oracle as oracle,
)


CONTROL_PLANE_SOURCE_DIR = Path(__file__).resolve().parent / "frontier_control_plane"
CONTROL_PLANE_SCHEMA_SOURCE = CONTROL_PLANE_SOURCE_DIR / "control_plane_common.py"
if str(CONTROL_PLANE_SOURCE_DIR) not in sys.path:
    sys.path.insert(0, str(CONTROL_PLANE_SOURCE_DIR))
from control_plane_common import (  # type: ignore[import-not-found]
    TRUSTED_LAUNCH_EXECUTOR,
    launch_contract_sha256,
    validate_launch_contract,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_ROOT = REPO_ROOT / "tst/publication/readiness"
AUTHORIZED_ORION_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
CANONICAL_PROJECT_HOME_ROOT = Path("/autofs/nccs-svm1_proj/ast207/proj-shared/PIC")
RUN_NAMESPACE = "runs/q043_registered_execution_raw_oracle_successor_v1"
PROJECT_HOME_RECEIPT_NAMESPACE = (
    "receipts/q043_registered_execution_raw_oracle_successor_v1"
)
DECK_ROOT = (
    REPO_ROOT / "inputs/tests/q043_bell_current_volume_aware_deposited_current_oracle"
)
DECK_MANIFEST = DECK_ROOT / "deck_manifest.json"
DECK_MANIFEST_SHA256 = (
    "77d6ec862314d308c41490dde1137f6c384d1da24b6b72ad43e357545cff98b0"
)

SCHEMA_VERSION = 1
SUCCESSOR_ID = "q043_registered_launch_policy_preparation_successor_v1"
MANIFEST_RECORD_TYPE = "q043_registered_launch_policy_preparation_manifest"
LAUNCH_RECORD_TYPE = "q043_registered_launch_review_candidate"
POLICY_RECORD_TYPE = "q043_registered_policy_slice_review_candidate"
POLICY_FRAGMENT_RECORD_TYPE = "q043_registered_policy_slice_review_fragment"
FINAL_BINDING_RECORD_TYPE = "q043_registered_launch_final_digest_bindings"
QUALIFICATION_EFFECT = (
    "source_local_launch_policy_preparation_only_no_launch_no_policy_no_q043_"
    "no_q023_no_q019_no_science_no_publication_authority"
)
CAMPAIGN = "q043_registered_execution_raw_oracle_successor_v1"
EVIDENCE_CLASS = "q043_registered_execution_raw_oracle_prerequisite"
PHYSICAL_MODE = "q043_bell_deposited_current_volume_aware_vl2_tsc"
RUNTIME_PROFILE = "frontier_minimum_supported"

EXPECTED_CASE_COUNT = 132
MAXIMUM_NODES_PER_CASE = 1
MAXIMUM_WALLTIME_SECONDS_PER_CASE = 600
MAXIMUM_ATTEMPTS_PER_CASE = 1
MAXIMUM_RETRIES_PER_CASE = 0
MAXIMUM_LIVE_Q043_SUBMISSIONS = 1
AUTHORIZED_TOTAL_NODE_HOUR_CAP = 10000.0
MAXIMUM_RAW_ARTIFACT_BYTES = 64 * 1024 * 1024
MAXIMUM_NON_RAW_CASE_BYTES = 64 * 1024 * 1024
MAXIMUM_BATCH_STORAGE_BYTES = 256 * 1024 * 1024 * 1024

SUBMISSION_ID_TEMPLATE = "{submission_id}"
RESERVATION_ID_TEMPLATE = "{reservation_id}"
FINAL_PLACEHOLDER_PREFIX = "PENDING_FINAL_"
FINAL_BINDING_KEYS = {
    "record_type",
    "schema_version",
    "source_commit",
    "source_bundle_sha256",
    "source_archive_path",
    "source_archive_sha256",
    "clean_candidate_manifest_path",
    "clean_candidate_manifest_sha256",
    "executable_path",
    "executable_sha256",
    "installed_control_plane_version",
    "orion_installed_control_plane_root",
    "project_home_installed_control_plane_root",
    "environment_profile_path",
    "environment_profile_sha256",
    "job_script_path",
    "job_script_sha256",
    "analysis_script_paths",
    "analysis_script_sha256",
    "reconcile_q043_registered_execution_path",
    "reconcile_q043_registered_execution_sha256",
}
PENDING_FINAL_BINDINGS = {
    "record_type": FINAL_BINDING_RECORD_TYPE,
    "schema_version": SCHEMA_VERSION,
    "source_commit": "PENDING_FINAL_SOURCE_COMMIT",
    "source_bundle_sha256": "PENDING_FINAL_SOURCE_BUNDLE_SHA256",
    "source_archive_path": "PENDING_FINAL_SOURCE_ARCHIVE_PATH",
    "source_archive_sha256": "PENDING_FINAL_SOURCE_ARCHIVE_SHA256",
    "clean_candidate_manifest_path": "PENDING_FINAL_CLEAN_CANDIDATE_MANIFEST_PATH",
    "clean_candidate_manifest_sha256": "PENDING_FINAL_CLEAN_CANDIDATE_MANIFEST_SHA256",
    "executable_path": "PENDING_FINAL_EXECUTABLE_PATH",
    "executable_sha256": "PENDING_FINAL_EXECUTABLE_SHA256",
    "installed_control_plane_version": "PENDING_FINAL_INSTALLED_CONTROL_PLANE_VERSION",
    "orion_installed_control_plane_root": (
        "PENDING_FINAL_ORION_INSTALLED_CONTROL_PLANE_ROOT"
    ),
    "project_home_installed_control_plane_root": (
        "PENDING_FINAL_PROJECT_HOME_INSTALLED_CONTROL_PLANE_ROOT"
    ),
    "environment_profile_path": "PENDING_FINAL_ENVIRONMENT_PROFILE_PATH",
    "environment_profile_sha256": "PENDING_FINAL_ENVIRONMENT_PROFILE_SHA256",
    "job_script_path": "PENDING_FINAL_JOB_SCRIPT_PATH",
    "job_script_sha256": "PENDING_FINAL_JOB_SCRIPT_SHA256",
    "analysis_script_paths": ["PENDING_FINAL_ANALYSIS_SCRIPT_PATHS"],
    "analysis_script_sha256": ["PENDING_FINAL_ANALYSIS_SCRIPT_SHA256"],
    "reconcile_q043_registered_execution_path": (
        "PENDING_FINAL_RECONCILE_Q043_REGISTERED_EXECUTION_PATH"
    ),
    "reconcile_q043_registered_execution_sha256": (
        "PENDING_FINAL_RECONCILE_Q043_REGISTERED_EXECUTION_SHA256"
    ),
}

REQUIRED_BLOCKERS = (
    "empty_user_queue_proved_at_fresh_reservation_boundary",
    "final_clean_candidate_manifest_source_archive_and_executable_bound",
    "final_paired_installed_control_plane_generation_bound_and_independently_verified",
    "installed_reconcile_q043_registered_execution_py_producer_bound",
    "exact_registered_policy_slice_promoted_by_installed_control_plane",
    "fresh_live_ledger_budget_and_storage_preflight_passed",
    "fresh_pre_manifest_pre_submit_wrapper_and_timeout_margin_attestations_bound",
    "fresh_submission_id_and_reservation_id_bound_once",
    "trusted_stdout_wrapper_evidence_support_installed",
    "paired_byte_identical_orion_and_canonical_project_home_receipts_required",
)

AUTHORIZATION_BOUNDARY = {
    "launch_authorized": False,
    "scheduler_submission_authorized": False,
    "policy_mutation_authorized": False,
    "frontier_execution_authorized": False,
    "q043_qualification_authorized": False,
    "q023_qualification_authorized": False,
    "q019_qualification_authorized": False,
    "scientific_claim_authorized": False,
    "publication_authorized": False,
}

_SHA256 = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT = re.compile(r"[0-9a-f]{40}")
_SAFE_ID = re.compile(r"[a-z0-9][a-z0-9_-]{0,127}")
_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH


class PreparationError(ValueError):
    """Reject incomplete, drifted, over-budget, or authority-bearing candidates."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PreparationError(message)


def _strict_equal(left: object, right: object) -> bool:
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return set(left) == set(right) and all(
            _strict_equal(left[key], right[key]) for key in left
        )
    if isinstance(left, list):
        return len(left) == len(right) and all(
            _strict_equal(a, b) for a, b in zip(left, right)
        )
    return left == right


def _json_bytes(value: object) -> bytes:
    try:
        return (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise PreparationError("candidate contains noncanonical JSON values") from error


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _binding(path: str, payload: bytes) -> dict[str, object]:
    return {"path": path, "sha256": _sha256_bytes(payload), "byte_count": len(payload)}


def _is_absolute_normalized_path(value: object) -> bool:
    if type(value) is not str or not value:
        return False
    path = PurePosixPath(value)
    return path.is_absolute() and path.as_posix() == value and ".." not in path.parts


def _contains_final_placeholder(value: object) -> bool:
    if isinstance(value, str):
        return FINAL_PLACEHOLDER_PREFIX in value
    if isinstance(value, dict):
        return any(_contains_final_placeholder(item) for item in value.values())
    if isinstance(value, list):
        return any(_contains_final_placeholder(item) for item in value)
    return False


def validate_final_bindings(value: object) -> dict[str, object]:
    """Validate one exact future final-digest binding without reading live state."""
    _require(type(value) is dict and set(value) == FINAL_BINDING_KEYS, "final binding schema drifted")
    result = dict(value)
    _require(
        result["record_type"] == FINAL_BINDING_RECORD_TYPE
        and type(result["schema_version"]) is int
        and result["schema_version"] == SCHEMA_VERSION,
        "final binding identity drifted",
    )
    _require(not _contains_final_placeholder(result), "final bindings retain a final-digest placeholder")
    for key in {
        "source_bundle_sha256",
        "source_archive_sha256",
        "clean_candidate_manifest_sha256",
        "executable_sha256",
        "installed_control_plane_version",
        "environment_profile_sha256",
        "job_script_sha256",
        "reconcile_q043_registered_execution_sha256",
    }:
        _require(
            type(result[key]) is str and _SHA256.fullmatch(str(result[key])) is not None,
            f"final binding {key} is not a lowercase SHA-256 digest",
        )
    _require(
        type(result["source_commit"]) is str
        and _GIT_COMMIT.fullmatch(str(result["source_commit"])) is not None,
        "final source commit is malformed",
    )
    for key in {
        "source_archive_path",
        "clean_candidate_manifest_path",
        "executable_path",
        "orion_installed_control_plane_root",
        "project_home_installed_control_plane_root",
        "environment_profile_path",
        "job_script_path",
        "reconcile_q043_registered_execution_path",
    }:
        _require(_is_absolute_normalized_path(result[key]), f"final binding {key} path is malformed")
    analysis_paths = result["analysis_script_paths"]
    analysis_sha256 = result["analysis_script_sha256"]
    _require(
        type(analysis_paths) is list
        and type(analysis_sha256) is list
        and 1 <= len(analysis_paths) == len(analysis_sha256) <= 16
        and all(_is_absolute_normalized_path(path) for path in analysis_paths)
        and all(type(digest) is str and _SHA256.fullmatch(digest) is not None for digest in analysis_sha256),
        "final analysis-script bindings are malformed",
    )
    version = str(result["installed_control_plane_version"])
    orion_controller = AUTHORIZED_ORION_ROOT / "control_plane" / version
    project_controller = CANONICAL_PROJECT_HOME_ROOT / "control_plane" / version
    _require(
        result["orion_installed_control_plane_root"] == str(orion_controller)
        and result["project_home_installed_control_plane_root"] == str(project_controller),
        "final paired installed-control-plane roots drifted",
    )
    _require(
        result["environment_profile_path"] == str(orion_controller / "frontier_pic_environment.sh")
        and result["job_script_path"] == str(orion_controller / "frontier_job.sh")
        and result["reconcile_q043_registered_execution_path"]
        == str(orion_controller / "reconcile_q043_registered_execution.py"),
        "final installed-control-plane member path drifted",
    )
    clean_manifest = PurePosixPath(str(result["clean_candidate_manifest_path"]))
    executable = PurePosixPath(str(result["executable_path"]))
    source_archive = PurePosixPath(str(result["source_archive_path"]))
    clean_root = PurePosixPath(str(AUTHORIZED_ORION_ROOT / "clean_candidates"))
    _require(
        clean_manifest.parent.parent == clean_root
        and clean_manifest.name == "clean_candidate_manifest.json"
        and executable.parent == clean_manifest.parent
        and source_archive.parent == clean_manifest.parent,
        "final clean-candidate path family drifted",
    )
    return result


def _selected_final_bindings(
    final_bindings: Mapping[str, object] | None,
) -> tuple[str, dict[str, object], list[str]]:
    if final_bindings is None:
        return (
            "source_local_preparation_pending_final_digests",
            dict(PENDING_FINAL_BINDINGS),
            sorted(FINAL_BINDING_KEYS - {"record_type", "schema_version"}),
        )
    return (
        "final_digests_bound_review_candidate_still_blocked",
        validate_final_bindings(final_bindings),
        [],
    )


def validate_case_matrix(cases: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    """Require the exact checked-in 132-case matrix, order, labels, and decompositions."""
    expected = [dict(case) for case in oracle.expected_cases()]
    contract_keys = set(expected[0])
    projected = [
        {key: case[key] for key in contract_keys}
        if type(case) is dict and contract_keys <= set(case)
        else None
        for case in cases
    ]
    _require(
        type(cases) in {list, tuple}
        and len(cases) == EXPECTED_CASE_COUNT
        and _strict_equal(projected, expected),
        "Q043 matrix must be the exact ordered 132-case checked-in contract",
    )
    return expected


def _checked_in_deck_manifest() -> dict[str, object]:
    _require(
        _sha256_file(DECK_MANIFEST) == DECK_MANIFEST_SHA256,
        "checked-in Q043 deck-manifest digest drifted",
    )
    try:
        manifest = oracle.validate_checked_in_decks()
    except oracle.ContractError as error:
        raise PreparationError("checked-in Q043 deck matrix failed validation") from error
    validate_case_matrix(manifest["cases"])
    return manifest


def _case_root_template(case_id: str) -> str:
    return str(AUTHORIZED_ORION_ROOT / RUN_NAMESPACE / case_id / SUBMISSION_ID_TEMPLATE)


def _raw_root_template(case_id: str) -> str:
    return f"{_case_root_template(case_id)}/raw"


def _artifact_root_template(case_id: str) -> str:
    return f"{_case_root_template(case_id)}/artifacts"


def _project_home_receipt_template(case_id: str) -> str:
    return str(
        CANONICAL_PROJECT_HOME_ROOT
        / PROJECT_HOME_RECEIPT_NAMESPACE
        / case_id
        / SUBMISSION_ID_TEMPLATE
        / "registered_execution_receipt.json"
    )


def _orion_receipt_template(case_id: str) -> str:
    return f"{_artifact_root_template(case_id)}/registration/registered_execution_receipt.json"


def _raw_relative_path(case: Mapping[str, object], field: str, cycle: int, rank: int) -> str:
    basename = str(case["case_id"]).replace("-", "_")
    filename = f"{basename}.{field}.{cycle:05d}.bin"
    if int(case["mpi_ranks"]) > 1:
        return f"bin/rank_{rank:08d}/{filename}"
    return f"bin/{filename}"


def _output_topology(case: Mapping[str, object], deck_text: str) -> dict[str, object]:
    blocks = oracle.parse_athinput_text(deck_text)
    ranks = int(case["mpi_ranks"])
    expected_per_rank = "true" if ranks > 1 else "false"
    for index, field in enumerate(oracle.FIELDS, 1):
        output = blocks.get(f"output{index}", {})
        _require(
            output.get("variable") == field
            and output.get("file_type") == "bin"
            and output.get("dcycle") == str(oracle.RAW_ORACLE_DCYCLE)
            and output.get("single_file_per_rank") == expected_per_rank,
            f"{case['case_id']}: checked-in output topology drifted",
        )
    inventory = [
        _raw_relative_path(case, field, cycle, rank)
        for cycle in (0, 1)
        for field in oracle.FIELDS
        for rank in range(ranks)
    ]
    return {
        "one_meshblock_per_rank": True,
        "single_file_per_rank": ranks > 1,
        "single_rank_file_topology": ranks == 1,
        "required_cycles": [0, 1],
        "required_fields": list(oracle.FIELDS),
        "mpi_rank_count": ranks,
        "expected_raw_artifact_count": len(inventory),
        "expected_relative_raw_paths": inventory,
        "maximum_bytes_per_raw_artifact": MAXIMUM_RAW_ARTIFACT_BYTES,
        "maximum_non_raw_case_bytes": MAXIMUM_NON_RAW_CASE_BYTES,
        "maximum_case_storage_bytes": (
            len(inventory) * MAXIMUM_RAW_ARTIFACT_BYTES + MAXIMUM_NON_RAW_CASE_BYTES
        ),
    }


def _launch_contract(case: Mapping[str, object]) -> dict[str, object]:
    case_id = str(case["case_id"])
    ranks = int(case["mpi_ranks"])
    contract = {
        "schema_version": 1,
        "executor": TRUSTED_LAUNCH_EXECUTOR,
        "pre_actions": [],
        "actions": [
            {
                "action_id": case_id,
                "kind": "athena",
                "resources": {
                    "nodes": MAXIMUM_NODES_PER_CASE,
                    "tasks": ranks,
                    "cpus_per_task": 1,
                    "gpus_per_task": 1,
                    "gpu_bind": "closest",
                },
                "arguments": [
                    {"literal": "-i"},
                    {"snapshot_role": "input-deck"},
                    {"literal": "-d"},
                    {"artifact_directory": "raw"},
                ],
                "stdout_artifact": "athena_stdout.txt",
                "stderr_artifact": "athena_stderr.txt",
            }
        ],
        "post_actions": [
            {
                "action_id": "require-stdout",
                "kind": "artifact_nonempty",
                "artifact": "athena_stdout.txt",
            },
            {
                "action_id": "sha-stdout",
                "kind": "artifact_sha256",
                "artifact": "athena_stdout.txt",
                "output_artifact": "athena_stdout.sha256",
            },
        ],
    }
    try:
        return validate_launch_contract(contract)
    except ValueError as error:
        raise PreparationError(f"{case_id}: installed launch-contract schema rejected candidate") from error


def _expected_command(
    case: Mapping[str, object], final_bindings: Mapping[str, object]
) -> list[str]:
    ranks = int(case["mpi_ranks"])
    artifact_root = _artifact_root_template(str(case["case_id"]))
    return [
        "srun",
        "--nodes=1",
        f"--ntasks={ranks}",
        f"--ntasks-per-node={ranks}",
        "--cpus-per-task=1",
        "--gpus-per-task=1",
        "--gpu-bind=closest",
        str(final_bindings["executable_path"]),
        "-i",
        f"{artifact_root}/registration/input_deck.athinput",
        "-d",
        _raw_root_template(str(case["case_id"])),
    ]


def _wrapper_evidence(case: Mapping[str, object]) -> dict[str, object]:
    ranks = int(case["mpi_ranks"])
    return {
        "required_exact_rank_line": (
            f"Q043_REGISTERED_EXECUTION case_id={case['case_id']} "
            f"mpi_world_size={ranks} rank_ids={','.join(str(rank) for rank in range(ranks))}"
        ),
        "required_exact_exit_line": "Q043_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0",
        "required_exact_termination_line": "Terminating on cycle limit",
        "required_terminal_regex_lines": [
            "time=<positive-finite> cycle=1",
            "tlim=<finite> nlim=1",
        ],
        "trusted_wrapper_must_emit_each_exact_line_once": True,
        "athena_stdout_alone_without_trusted_wrapper_evidence_sufficient": False,
    }


def _reconciliation_interface(
    case: Mapping[str, object], final_bindings: Mapping[str, object]
) -> dict[str, object]:
    case_id = str(case["case_id"])
    return {
        "producer_name": "reconcile_q043_registered_execution.py",
        "producer_path": final_bindings["reconcile_q043_registered_execution_path"],
        "producer_sha256": final_bindings["reconcile_q043_registered_execution_sha256"],
        "installed_generation_only": True,
        "run_only_after_canonical_mirrored_ledger_reconciliation": True,
        "required_receipt_record_type": "q043_reconciled_registered_execution_receipt",
        "required_receipt_role": "immutable_reconciled_registered_execution",
        "orion_receipt_path_template": _orion_receipt_template(case_id),
        "canonical_project_home_receipt_path_template": _project_home_receipt_template(case_id),
        "paired_receipts_required": True,
        "paired_receipts_must_be_byte_identical": True,
        "paired_receipts_must_be_recursively_read_only": True,
        "receipt_sha256": "fresh_at_reconciliation_not_assumed_by_preparation",
        "raw_inventory_captured_at_reconciliation": True,
        "required_bound_evidence": [
            "source_commit_and_source_bundle",
            "clean_candidate_manifest_source_archive_and_executable",
            "checked_in_and_immutable_launch_deck",
            "installed_environment_profile",
            "exact_launch_contract_and_command",
            "actual_mpi_world_rank_and_meshblock_topology",
            "slurm_job_identity_terminal_state_and_exit_code",
            "stdout_wrapper_terminal_receipt_and_pre_submit_manifest",
            "exact_raw_filesystem_inventory_with_sha256_and_byte_counts",
            "canonical_mirrored_ledger_reconciliation_event",
        ],
    }


def _deck_binding(case: Mapping[str, object], manifest_record: Mapping[str, object]) -> dict[str, object]:
    relative = (
        "inputs/tests/q043_bell_current_volume_aware_deposited_current_oracle/"
        f"{case['case_id']}.athinput"
    )
    path = REPO_ROOT / relative
    payload = path.read_bytes()
    expected = oracle.render_oracle_deck(case).encode("utf-8")
    _require(payload == expected, f"{case['case_id']}: checked-in deck bytes drifted")
    _require(
        manifest_record["deck_path"] == path.name
        and manifest_record["deck_sha256"] == _sha256_bytes(payload),
        f"{case['case_id']}: deck-manifest cross-link drifted",
    )
    return {
        "path": relative,
        "sha256": _sha256_bytes(payload),
        "byte_count": len(payload),
    }


def _case_identity(case: Mapping[str, object], index: int) -> dict[str, object]:
    case_id = str(case["case_id"])
    attempt_id = f"q043-re-{index:03d}-{case_id}"
    authorization_id = f"q043-re-{index:03d}-v1"
    _require(
        _SAFE_ID.fullmatch(attempt_id) is not None
        and _SAFE_ID.fullmatch(authorization_id) is not None,
        f"{case_id}: derived identity is unsafe",
    )
    return {
        "case_index": index,
        "case_id": case_id,
        "attempt_id": attempt_id,
        "authorization_id": authorization_id,
        "maximum_registered_attempts": MAXIMUM_ATTEMPTS_PER_CASE,
        "maximum_retries": MAXIMUM_RETRIES_PER_CASE,
        "fresh_submission_id_template": SUBMISSION_ID_TEMPLATE,
        "fresh_reservation_id_template": RESERVATION_ID_TEMPLATE,
    }


def _launch_candidate(
    *,
    case: Mapping[str, object],
    index: int,
    deck: Mapping[str, object],
    deck_text: str,
    binding_stage: str,
    final_bindings: Mapping[str, object],
    unresolved_final_bindings: Sequence[str],
) -> dict[str, object]:
    identity = _case_identity(case, index)
    case_id = str(case["case_id"])
    contract = _launch_contract(case)
    topology = _output_topology(case, deck_text)
    candidate = {
        "record_type": LAUNCH_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "status": "source_local_review_candidate_incomplete_launch_prohibited",
        "qualification_effect": QUALIFICATION_EFFECT,
        "campaign": CAMPAIGN,
        "identity": identity,
        "case_contract": dict(case),
        "checked_in_deck": dict(deck),
        "binding_stage": binding_stage,
        "selected_final_bindings": dict(final_bindings),
        "unresolved_final_bindings": list(unresolved_final_bindings),
        "environment": {
            "profile_path": final_bindings["environment_profile_path"],
            "profile_sha256": final_bindings["environment_profile_sha256"],
            "installed_control_plane_version": final_bindings["installed_control_plane_version"],
            "caller_environment_authoritative": False,
            "exact_installed_profile_required": True,
        },
        "roots": {
            "authorized_orion_case_root_template": _case_root_template(case_id),
            "raw_output_root_template": _raw_root_template(case_id),
            "artifact_root_template": _artifact_root_template(case_id),
            "canonical_project_home_receipt_root_template": str(
                CANONICAL_PROJECT_HOME_ROOT
                / PROJECT_HOME_RECEIPT_NAMESPACE
                / case_id
                / SUBMISSION_ID_TEMPLATE
            ),
        },
        "command_template": _expected_command(case, final_bindings),
        "launch_contract_candidate": contract,
        "launch_contract_sha256": launch_contract_sha256(contract),
        "mpi_and_decomposition": {
            "mpi_ranks": case["mpi_ranks"],
            "requested_nodes": MAXIMUM_NODES_PER_CASE,
            "requested_tasks": case["mpi_ranks"],
            "meshblock_grid": case["meshblock_grid"],
            "meshblock_nx": case["meshblock_nx"],
            "global_nx": case["global_nx"],
            "partitioned_axes": case["partitioned_axes"],
            "one_meshblock_per_rank": True,
            "rank_meshblock_ids": [
                {"rank": rank, "meshblock_ids": [rank]}
                for rank in range(int(case["mpi_ranks"]))
            ],
        },
        "output_topology": topology,
        "stdout_wrapper_evidence": _wrapper_evidence(case),
        "trusted_reconciliation_interface": _reconciliation_interface(case, final_bindings),
        "resource_and_retry_ceiling": {
            "selected_qos": "normal",
            "registered_short_nonproduction": False,
            "maximum_nodes": MAXIMUM_NODES_PER_CASE,
            "maximum_walltime_seconds": MAXIMUM_WALLTIME_SECONDS_PER_CASE,
            "maximum_attempts": MAXIMUM_ATTEMPTS_PER_CASE,
            "maximum_retries": MAXIMUM_RETRIES_PER_CASE,
            "maximum_node_hours": (
                MAXIMUM_NODES_PER_CASE
                * MAXIMUM_WALLTIME_SECONDS_PER_CASE
                * MAXIMUM_ATTEMPTS_PER_CASE
                / 3600.0
            ),
            "maximum_storage_bytes": topology["maximum_case_storage_bytes"],
        },
        "required_blockers": list(REQUIRED_BLOCKERS),
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }
    _require(
        binding_stage != "final_digests_bound_review_candidate_still_blocked"
        or not _contains_final_placeholder(candidate),
        f"{case_id}: final-digest candidate retains a placeholder",
    )
    return candidate


def _policy_candidate(launch: Mapping[str, object]) -> dict[str, object]:
    identity = launch["identity"]
    deck = launch["checked_in_deck"]
    final_bindings = launch["selected_final_bindings"]
    ceiling = launch["resource_and_retry_ceiling"]
    result = {
        "record_type": POLICY_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "status": "review_required_not_authorized_not_live_policy",
        "qualification_effect": QUALIFICATION_EFFECT,
        "campaign": CAMPAIGN,
        "case_id": identity["case_id"],
        "attempt_id": identity["attempt_id"],
        "authorization_id": identity["authorization_id"],
        "policy_slice_candidate": {
            "authorization_id": identity["authorization_id"],
            "status": "review_required_not_authorized",
            "campaign": CAMPAIGN,
            "test_id": str(identity["case_id"]).replace("-", "_"),
            "evidence_class": EVIDENCE_CLASS,
            "physical_mode": PHYSICAL_MODE,
            "runtime_profile": RUNTIME_PROFILE,
            "selected_qos": ceiling["selected_qos"],
            "registered_short_nonproduction": ceiling["registered_short_nonproduction"],
            "maximum_nodes": ceiling["maximum_nodes"],
            "maximum_walltime_seconds": ceiling["maximum_walltime_seconds"],
            "maximum_attempts": ceiling["maximum_attempts"],
            "job_script_sha256": final_bindings["job_script_sha256"],
            "input_deck_sha256": deck["sha256"],
            "environment_profile_sha256": final_bindings["environment_profile_sha256"],
            "analysis_script_sha256": list(final_bindings["analysis_script_sha256"]),
            "executable_sha256": final_bindings["executable_sha256"],
            "launch_contract_sha256": launch["launch_contract_sha256"],
            "clean_candidate_manifest_sha256": final_bindings[
                "clean_candidate_manifest_sha256"
            ],
        },
        "launch_candidate_sha256": _sha256_bytes(_json_bytes(launch)),
        "unresolved_final_bindings": launch["unresolved_final_bindings"],
        "fresh_runtime_bindings_required": [
            "submission_id",
            "reservation_id",
            "empty_user_queue_snapshot",
            "live_ledger_budget_snapshot",
            "pre_manifest_attestation",
            "pre_submit_wrapper_attestation",
            "timeout_margin_artifact",
            "paired_reconciliation_receipt_sha256",
        ],
        "required_blockers": list(REQUIRED_BLOCKERS),
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }
    _require(
        launch["binding_stage"] != "final_digests_bound_review_candidate_still_blocked"
        or not _contains_final_placeholder(result),
        f"{identity['case_id']}: final-digest policy candidate retains a placeholder",
    )
    return result


def _budget_accounting(launches: Sequence[Mapping[str, object]]) -> dict[str, object]:
    total_node_hours = (
        len(launches)
        * MAXIMUM_NODES_PER_CASE
        * MAXIMUM_WALLTIME_SECONDS_PER_CASE
        * MAXIMUM_ATTEMPTS_PER_CASE
        / 3600.0
    )
    total_storage = sum(
        int(item["resource_and_retry_ceiling"]["maximum_storage_bytes"])
        for item in launches
    )
    result = {
        "case_count": len(launches),
        "maximum_live_q043_submissions": MAXIMUM_LIVE_Q043_SUBMISSIONS,
        "maximum_attempts_per_case": MAXIMUM_ATTEMPTS_PER_CASE,
        "maximum_retries_per_case": MAXIMUM_RETRIES_PER_CASE,
        "maximum_nodes_per_case": MAXIMUM_NODES_PER_CASE,
        "maximum_walltime_seconds_per_case": MAXIMUM_WALLTIME_SECONDS_PER_CASE,
        "maximum_registered_attempts": len(launches) * MAXIMUM_ATTEMPTS_PER_CASE,
        "maximum_batch_node_hours": total_node_hours,
        "authorized_total_node_hour_cap": AUTHORIZED_TOTAL_NODE_HOUR_CAP,
        "minimum_live_ledger_remaining_node_hours_required": total_node_hours,
        "live_ledger_remaining_node_hours": "fresh_at_policy_promotion_not_assumed",
        "maximum_batch_storage_bytes": total_storage,
        "maximum_storage_cap_bytes": MAXIMUM_BATCH_STORAGE_BYTES,
        "fresh_storage_preflight_required": True,
        "ledger_mutation_authorized": False,
    }
    return validate_budget_accounting(result, launches)


def _aggregate_policy_fragment(
    *,
    policies: Sequence[Mapping[str, object]],
    budget: Mapping[str, object],
    binding_stage: str,
    unresolved_final_bindings: Sequence[str],
) -> dict[str, object]:
    return {
        "record_type": POLICY_FRAGMENT_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "status": "aggregate_review_fragment_only_not_live_policy",
        "qualification_effect": QUALIFICATION_EFFECT,
        "campaign": CAMPAIGN,
        "binding_stage": binding_stage,
        "case_count": len(policies),
        "registered_science_slice_candidates": [
            dict(policy["policy_slice_candidate"]) for policy in policies
        ],
        "declared_batch_ceiling": dict(budget),
        "unresolved_final_bindings": list(unresolved_final_bindings),
        "required_blockers": list(REQUIRED_BLOCKERS),
        "execution_boundary": {
            "complete_storage_policy": False,
            "live_policy_mutation_authorized": False,
            "scheduler_submission_authorized": False,
            "frontier_execution_authorized": False,
            "scientific_claim_authorized": False,
            "publication_authorized": False,
        },
    }


def validate_budget_accounting(
    value: object, launches: Sequence[Mapping[str, object]]
) -> dict[str, object]:
    """Require exact batch ceilings and reject node-hour or storage overflow."""
    _require(type(value) is dict, "budget accounting must be an object")
    for item in launches:
        ceiling = item.get("resource_and_retry_ceiling")
        topology = item.get("output_topology")
        _require(
            type(ceiling) is dict
            and type(topology) is dict
            and ceiling.get("maximum_nodes") == MAXIMUM_NODES_PER_CASE
            and ceiling.get("maximum_walltime_seconds")
            == MAXIMUM_WALLTIME_SECONDS_PER_CASE
            and ceiling.get("maximum_attempts") == MAXIMUM_ATTEMPTS_PER_CASE
            and ceiling.get("maximum_retries") == MAXIMUM_RETRIES_PER_CASE
            and ceiling.get("maximum_node_hours")
            == (
                MAXIMUM_NODES_PER_CASE
                * MAXIMUM_WALLTIME_SECONDS_PER_CASE
                * MAXIMUM_ATTEMPTS_PER_CASE
                / 3600.0
            )
            and ceiling.get("maximum_storage_bytes")
            == topology.get("maximum_case_storage_bytes"),
            "Q043 per-case resource, retry, or storage ceiling drifted",
        )
    expected_node_hours = (
        len(launches)
        * MAXIMUM_NODES_PER_CASE
        * MAXIMUM_WALLTIME_SECONDS_PER_CASE
        * MAXIMUM_ATTEMPTS_PER_CASE
        / 3600.0
    )
    expected_storage = sum(
        int(item["resource_and_retry_ceiling"]["maximum_storage_bytes"])
        for item in launches
    )
    _require(
        value.get("case_count") == EXPECTED_CASE_COUNT
        and value.get("maximum_registered_attempts") == EXPECTED_CASE_COUNT
        and value.get("maximum_attempts_per_case") == 1
        and value.get("maximum_retries_per_case") == 0
        and value.get("maximum_live_q043_submissions") == 1,
        "budget accounting attempt or concurrency ceiling drifted",
    )
    _require(
        value.get("maximum_batch_node_hours") == expected_node_hours
        and expected_node_hours <= AUTHORIZED_TOTAL_NODE_HOUR_CAP,
        "Q043 batch node-hour ceiling overflow or drift",
    )
    _require(
        value.get("maximum_batch_storage_bytes") == expected_storage
        and expected_storage <= MAXIMUM_BATCH_STORAGE_BYTES,
        "Q043 batch storage ceiling overflow or drift",
    )
    _require(
        value.get("live_ledger_remaining_node_hours")
        == "fresh_at_policy_promotion_not_assumed"
        and value.get("fresh_storage_preflight_required") is True
        and value.get("ledger_mutation_authorized") is False,
        "budget accounting acquired authority or assumed fresh live state",
    )
    return dict(value)


def build_materialization(
    final_bindings: Mapping[str, object] | None = None,
) -> tuple[dict[str, object], dict[str, bytes]]:
    """Build deterministic source-local launch and policy review candidates."""
    deck_manifest = _checked_in_deck_manifest()
    cases = validate_case_matrix(deck_manifest["cases"])
    records_by_id = {
        str(record["case_id"]): record for record in deck_manifest["cases"]
    }
    binding_stage, selected_bindings, unresolved = _selected_final_bindings(final_bindings)
    files: dict[str, bytes] = {}
    launches: list[dict[str, object]] = []
    policies: list[dict[str, object]] = []
    case_records: list[dict[str, object]] = []
    for index, case in enumerate(cases, 1):
        case_id = str(case["case_id"])
        deck = _deck_binding(case, records_by_id[case_id])
        deck_text = (REPO_ROOT / str(deck["path"])).read_text(encoding="utf-8")
        launch = _launch_candidate(
            case=case,
            index=index,
            deck=deck,
            deck_text=deck_text,
            binding_stage=binding_stage,
            final_bindings=selected_bindings,
            unresolved_final_bindings=unresolved,
        )
        policy = _policy_candidate(launch)
        launch_path = f"launch_candidates/{case_id}.json"
        policy_path = f"policy_candidates/{case_id}.json"
        launch_payload = _json_bytes(launch)
        policy_payload = _json_bytes(policy)
        files[launch_path] = launch_payload
        files[policy_path] = policy_payload
        launches.append(launch)
        policies.append(policy)
        case_records.append(
            {
                "case_index": index,
                "case_id": case_id,
                "dimension": case["dimension"],
                "resolution": case["resolution"],
                "ppc": case["ppc"],
                "decomposition": case["decomposition"],
                "mpi_ranks": case["mpi_ranks"],
                "meshblock_grid": case["meshblock_grid"],
                "checked_in_deck": deck,
                "launch_candidate": _binding(launch_path, launch_payload),
                "policy_candidate": _binding(policy_path, policy_payload),
            }
        )
    budget = _budget_accounting(launches)
    policy_fragment = _aggregate_policy_fragment(
        policies=policies,
        budget=budget,
        binding_stage=binding_stage,
        unresolved_final_bindings=unresolved,
    )
    reconciliation = {
        "producer_interface": "reconcile_q043_registered_execution.py",
        "producer_path": selected_bindings["reconcile_q043_registered_execution_path"],
        "producer_sha256": selected_bindings[
            "reconcile_q043_registered_execution_sha256"
        ],
        "installed_paired_generation_required": True,
        "orion_root": str(AUTHORIZED_ORION_ROOT),
        "canonical_project_home_root": str(CANONICAL_PROJECT_HOME_ROOT),
        "paired_byte_identical_receipts_required_for_every_case": True,
        "receipt_digest_assumed_by_this_materializer": False,
    }
    budget_payload = _json_bytes(budget)
    policy_fragment_payload = _json_bytes(policy_fragment)
    reconciliation_payload = _json_bytes(reconciliation)
    files["batch_budget_accounting_input.json"] = budget_payload
    files["registered_policy_slice_review_fragment.json"] = policy_fragment_payload
    files["trusted_reconciliation_producer_interface.json"] = reconciliation_payload
    manifest = {
        "record_type": MANIFEST_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "status": "source_local_preparation_complete_registration_and_launch_blocked",
        "qualification_effect": QUALIFICATION_EFFECT,
        "campaign": CAMPAIGN,
        "binding_stage": binding_stage,
        "selected_final_bindings": selected_bindings,
        "unresolved_final_bindings": unresolved,
        "source_bindings": {
            "checked_in_deck_manifest": {
                "path": str(DECK_MANIFEST.relative_to(REPO_ROOT)),
                "sha256": DECK_MANIFEST_SHA256,
            },
            "source_local_oracle": {
                "path": str(Path(oracle.__file__).resolve().relative_to(REPO_ROOT)),
                "sha256": _sha256_file(Path(oracle.__file__).resolve()),
            },
            "checked_in_launch_contract_schema_source": {
                "path": str(CONTROL_PLANE_SCHEMA_SOURCE.relative_to(REPO_ROOT)),
                "sha256": _sha256_file(CONTROL_PLANE_SCHEMA_SOURCE),
            },
        },
        "launch_contract_schema_compatibility": {
            "checked_in_control_plane_schema_source_validation": (
                "passed_for_all_132_candidates"
            ),
            "final_installed_generation_schema_validation_required": True,
            "final_installed_generation_assumed": False,
        },
        "exact_matrix_contract": {
            "case_count": EXPECTED_CASE_COUNT,
            "one_registered_attempt_per_case": True,
            "one_meshblock_per_rank": True,
            "dimensions": list(oracle.DIMENSIONS),
            "resolutions": list(oracle.RESOLUTIONS),
            "ppc": list(oracle.PPC_VALUES),
            "decompositions_by_dimension": {
                str(key): list(value)
                for key, value in oracle.DECOMPOSITIONS_BY_DIMENSION.items()
            },
            "artificial_c_over_v_cr": list(oracle.ARTIFICIAL_C_OVER_V_CR_VALUES),
        },
        "case_records": case_records,
        "case_count": len(case_records),
        "launch_candidate_count": len(launches),
        "policy_candidate_count": len(policies),
        "batch_budget_accounting_input": _binding(
            "batch_budget_accounting_input.json", budget_payload
        ),
        "registered_policy_slice_review_fragment": _binding(
            "registered_policy_slice_review_fragment.json", policy_fragment_payload
        ),
        "trusted_reconciliation_producer_interface": _binding(
            "trusted_reconciliation_producer_interface.json", reconciliation_payload
        ),
        "required_blockers": list(REQUIRED_BLOCKERS),
        "execution_boundary": {
            "empty_user_queue_assumed": False,
            "final_clean_candidate_assumed": False,
            "final_installed_control_plane_assumed": False,
            "live_policy_complete": False,
            "scheduler_calls_authorized": False,
            **AUTHORIZATION_BOUNDARY,
        },
    }
    _require(
        binding_stage != "final_digests_bound_review_candidate_still_blocked"
        or (
            not _contains_final_placeholder(manifest)
            and not any(_contains_final_placeholder(json.loads(payload)) for payload in files.values())
        ),
        "post-final-digest materialization retains a final placeholder",
    )
    return manifest, files


def validate_materialization(
    manifest: object,
    files: Mapping[str, bytes],
    final_bindings: Mapping[str, object] | None = None,
) -> tuple[dict[str, object], dict[str, bytes]]:
    """Rebuild and exactly validate all source-local preparation bytes."""
    expected_manifest, expected_files = build_materialization(final_bindings)
    _require(_strict_equal(manifest, expected_manifest), "Q043 preparation manifest drifted")
    _require(set(files) == set(expected_files), "Q043 preparation file inventory drifted")
    for path, payload in expected_files.items():
        _require(type(files[path]) is bytes and files[path] == payload, f"Q043 preparation file drifted: {path}")
    for record in expected_manifest["case_records"]:
        launch = json.loads(files[record["launch_candidate"]["path"]])
        policy = json.loads(files[record["policy_candidate"]["path"]])
        _require(
            validate_launch_contract(launch["launch_contract_candidate"])
            == launch["launch_contract_candidate"],
            f"{record['case_id']}: installed launch-contract schema validation drifted",
        )
        _require(
            launch["authorization"] == AUTHORIZATION_BOUNDARY
            and policy["authorization"] == AUTHORIZATION_BOUNDARY,
            f"{record['case_id']}: candidate acquired authority",
        )
    budget = json.loads(files["batch_budget_accounting_input.json"])
    policy_fragment = json.loads(files["registered_policy_slice_review_fragment.json"])
    launches = [
        json.loads(files[record["launch_candidate"]["path"]])
        for record in expected_manifest["case_records"]
    ]
    validate_budget_accounting(budget, launches)
    _require(
        policy_fragment
        == _aggregate_policy_fragment(
            policies=[
                json.loads(files[record["policy_candidate"]["path"]])
                for record in expected_manifest["case_records"]
            ],
            budget=budget,
            binding_stage=expected_manifest["binding_stage"],
            unresolved_final_bindings=expected_manifest["unresolved_final_bindings"],
        ),
        "Q043 aggregate policy review fragment drifted or acquired authority",
    )
    return expected_manifest, expected_files


def _write_new(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags, 0o444)
    try:
        remaining = memoryview(payload)
        while remaining:
            written = os.write(descriptor, remaining)
            _require(written > 0, f"short write while materializing {path}")
            remaining = remaining[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _freeze_tree(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o555 if path.is_dir() else 0o444)
    root.chmod(0o555)


def materialize_review_bundle(
    output_root: Path,
    final_bindings: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Write one recursively read-only review bundle outside live control namespaces."""
    output_root = Path(os.path.abspath(output_root))
    _require(
        output_root.parent.is_dir() and output_root.parent.resolve() == output_root.parent,
        "output parent must be one existing canonical directory",
    )
    _require(not output_root.exists(), "output root already exists")
    for protected in (AUTHORIZED_ORION_ROOT, CANONICAL_PROJECT_HOME_ROOT):
        try:
            output_root.relative_to(protected)
        except ValueError:
            continue
        raise PreparationError("review bundle cannot be materialized in a live PIC namespace")
    manifest, files = build_materialization(final_bindings)
    validate_materialization(manifest, files, final_bindings)
    try:
        output_root.mkdir(mode=0o700)
        bindings = []
        for relative, payload in sorted(files.items()):
            _write_new(output_root / relative, payload)
            bindings.append(_binding(relative, payload))
        manifest = {**manifest, "materialized_files": bindings}
        manifest_payload = _json_bytes(manifest)
        _write_new(output_root / "materialization_manifest.json", manifest_payload)
        _freeze_tree(output_root)
        return {
            "output_root": str(output_root),
            "materialization_manifest_sha256": _sha256_bytes(manifest_payload),
            "case_count": EXPECTED_CASE_COUNT,
            "binding_stage": manifest["binding_stage"],
            "maximum_batch_node_hours": json.loads(
                files["batch_budget_accounting_input.json"]
            )["maximum_batch_node_hours"],
            "recursively_read_only": all(
                not path.stat().st_mode & _WRITE_BITS
                for path in [output_root, *output_root.rglob("*")]
            ),
            "launch_authorized": False,
            "policy_mutation_authorized": False,
            "publication_authorized": False,
        }
    except BaseException:
        if output_root.exists():
            for path in [output_root, *output_root.rglob("*")]:
                if not path.is_symlink():
                    path.chmod(path.stat().st_mode | stat.S_IWUSR)
            shutil.rmtree(output_root)
        raise


def _read_final_bindings(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PreparationError("final bindings file is not readable UTF-8 JSON") from error
    return validate_final_bindings(value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--final-bindings", type=Path)
    arguments = parser.parse_args()
    final_bindings = (
        _read_final_bindings(arguments.final_bindings)
        if arguments.final_bindings is not None
        else None
    )
    result = materialize_review_bundle(arguments.output_root, final_bindings)
    print(_json_bytes(result).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
