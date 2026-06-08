#!/usr/bin/env python3
"""Prepare non-authorizing Q023 registered launch and policy candidates."""

from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
from typing import Any, Mapping, Sequence

from tst.publication import analyze_q023_paper_bell_linear_joverc as bell
from tst.publication import q043_registered_launch_policy_preparation_successor_v1 as q043


CONTROL_PLANE_SOURCE_DIR = Path(__file__).resolve().parent / "frontier_control_plane"
if str(CONTROL_PLANE_SOURCE_DIR) not in os.sys.path:
    os.sys.path.insert(0, str(CONTROL_PLANE_SOURCE_DIR))
from control_plane_common import (  # type: ignore[import-not-found]
    TRUSTED_LAUNCH_EXECUTOR,
    launch_contract_sha256,
    validate_launch_contract,
    validate_storage_policy,
)
from revalidate_clean_candidate import (  # type: ignore[import-not-found]  # noqa: E402
    revalidate_clean_candidate,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
AUTHORIZED_ORION_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
CANONICAL_PROJECT_HOME_ROOT = Path("/autofs/nccs-svm1_proj/ast207/proj-shared/PIC")
RUN_NAMESPACE = "runs/q023_paper_bell_linear_joverc_registered_successor_v1"
PROJECT_HOME_RECEIPT_NAMESPACE = "ledger/q023_registered_execution_receipts"
CAMPAIGN = "q023_paper_bell_linear_joverc_registered_successor_v1"
EVIDENCE_CLASS = "q023_corrected_linear_bell_registered_qualification"
PHYSICAL_MODE = "q023_corrected_joverc_linear_bell_vl2_tsc"
RUNTIME_PROFILE = "frontier_minimum_supported"
SUCCESSOR_ID = "q023_registered_launch_policy_preparation_successor_v1"
SCHEMA_VERSION = 1
MANIFEST_RECORD_TYPE = "q023_registered_launch_policy_preparation_manifest"
LAUNCH_RECORD_TYPE = "q023_registered_launch_review_candidate"
POLICY_RECORD_TYPE = "q023_registered_policy_slice_review_candidate"
POLICY_FRAGMENT_RECORD_TYPE = "q023_registered_policy_slice_review_fragment"
FINAL_BINDING_RECORD_TYPE = "q023_registered_launch_final_digest_bindings"
EXPECTED_CASE_COUNT = 55
MAXIMUM_NODES_PER_CASE = 1
MAXIMUM_WALLTIME_SECONDS_PER_CASE = 3600
ATHENA_WALLTIME_SECONDS = 3300
MAXIMUM_ATTEMPTS_PER_CASE = 1
MAXIMUM_RETRIES_PER_CASE = 0
MAXIMUM_LIVE_Q023_SUBMISSIONS = 1
MAXIMUM_NON_RAW_CASE_BYTES = 512 * 1024 * 1024
MAXIMUM_CASE_STORAGE_BYTES = 4 * 1024 * 1024 * 1024
MAXIMUM_BATCH_STORAGE_BYTES = 256 * 1024 * 1024 * 1024
MAXIMUM_BATCH_NODE_HOURS = 55.0
REQUIRED_OUTPUT_INDICES = tuple(range(89))
REQUIRED_OUTPUT_FIELDS = tuple(bell._RAW_OUTPUT_VARIABLES)
MAXIMUM_RAW_ARTIFACT_BYTES_BY_FIELD = {
    "mhd_w_bcc": 16 * 1024 * 1024,
    "prtcl_rho": 4 * 1024 * 1024,
    "prtcl_jx": 4 * 1024 * 1024,
    "prtcl_jy": 4 * 1024 * 1024,
    "prtcl_jz": 4 * 1024 * 1024,
}
SUBMISSION_ID_TEMPLATE = "{submission_id}"
RESERVATION_ID_TEMPLATE = "{reservation_id}"
FINAL_PLACEHOLDER_PREFIX = "PENDING_FINAL_"
QUALIFICATION_EFFECT = (
    "source_local_q023_launch_policy_preparation_only_no_launch_no_policy_"
    "no_linear_claim_no_q019_no_science_no_publication_authority"
)
AUTHORIZATION_BOUNDARY = {
    "launch_authorized": False,
    "scheduler_submission_authorized": False,
    "policy_mutation_authorized": False,
    "frontier_execution_authorized": False,
    "q023_qualification_authorized": False,
    "q019_qualification_authorized": False,
    "scientific_claim_authorized": False,
    "publication_authorized": False,
}
REQUIRED_BLOCKERS = (
    "exact_registered_q043_matrix_prerequisite_bound_and_revalidated",
    "q023_source_decks_rebound_to_exact_q043_matrix_before_execution",
    "empty_user_queue_proved_at_fresh_reservation_boundary",
    "final_clean_candidate_manifest_source_archive_and_executable_bound",
    "final_paired_installed_control_plane_generation_bound_and_independently_verified",
    "installed_reconcile_q023_registered_execution_py_producer_bound",
    "exact_registered_policy_slice_promoted_by_installed_control_plane",
    "fresh_live_ledger_budget_and_storage_preflight_passed",
    "fresh_pre_manifest_pre_submit_wrapper_and_timeout_margin_attestations_bound",
    "fresh_submission_id_and_reservation_id_bound_once",
    "trusted_stdout_wrapper_evidence_support_installed",
    "paired_byte_identical_orion_and_canonical_project_home_receipts_required",
)
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
    "reconcile_q023_registered_execution_path",
    "reconcile_q023_registered_execution_sha256",
    "q043_registered_matrix_path",
    "q043_registered_matrix_sha256",
    "q043_registered_matrix_record_type",
    "q043_registered_matrix_case_bindings_sha256",
    "q043_registered_dependency_sha256",
}
PENDING_FINAL_BINDINGS: dict[str, object] = {
    "record_type": FINAL_BINDING_RECORD_TYPE,
    "schema_version": SCHEMA_VERSION,
    **{
        key: f"{FINAL_PLACEHOLDER_PREFIX}{key.upper()}"
        for key in FINAL_BINDING_KEYS
        - {"record_type", "schema_version", "analysis_script_paths", "analysis_script_sha256"}
    },
    "analysis_script_paths": ["PENDING_FINAL_ANALYSIS_SCRIPT_PATHS"],
    "analysis_script_sha256": ["PENDING_FINAL_ANALYSIS_SCRIPT_SHA256"],
}
_SHA256 = re.compile(r"[0-9a-f]{64}")
_COMMIT = re.compile(r"[0-9a-f]{40}")
_SUBMISSION = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}"
)
_STORAGE_POLICY_REQUIRED_KEYS = {
    "schema_version",
    "frontier",
    "science_submission_freeze",
    "registered_science_slices",
    "frontier_admission_smoke",
    "olcf_side_storage",
    "long_term_storage",
}
_STORAGE_POLICY_ALLOWED_KEYS = _STORAGE_POLICY_REQUIRED_KEYS | {
    "authorization_date",
    "authorized_by",
    "reviewer",
}


class PreparationError(ValueError):
    """Reject incomplete, drifted, over-budget, or authority-bearing candidates."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PreparationError(message)


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _write_json_exclusive(path: Path, value: object) -> None:
    payload = _json_bytes(value)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o444,
    )
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError(f"short write: {path}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


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


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _binding(path: str, payload: bytes) -> dict[str, object]:
    return {"path": path, "sha256": _sha256_bytes(payload), "byte_count": len(payload)}


def _contains_placeholder(value: object) -> bool:
    if isinstance(value, str):
        return value.startswith(FINAL_PLACEHOLDER_PREFIX)
    if isinstance(value, Mapping):
        return any(_contains_placeholder(item) for item in value.values())
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return any(_contains_placeholder(item) for item in value)
    return False


def _stable_binding(
    path: object,
    digest: object,
    *,
    label: str,
    executable: bool = False,
) -> str:
    _require(
        isinstance(path, str) and Path(path).is_absolute(),
        f"{label} path is not absolute",
    )
    _require(
        isinstance(digest, str) and _SHA256.fullmatch(digest) is not None,
        f"{label} digest is malformed",
    )
    _, payload = q043._stable_regular_bytes(
        Path(path), label=label, require_executable=executable
    )
    _require(_sha256_bytes(payload) == digest, f"{label} digest drifted")
    return path


def validate_final_bindings(value: object) -> dict[str, object]:
    """Validate exact final identities without assuming mutable live files."""
    _require(
        type(value) is dict and set(value) == FINAL_BINDING_KEYS,
        "Q023 final-binding keys drifted",
    )
    final = dict(value)
    _require(
        final["record_type"] == FINAL_BINDING_RECORD_TYPE
        and type(final["schema_version"]) is int
        and final["schema_version"] == SCHEMA_VERSION
        and not _contains_placeholder(final),
        "Q023 final-binding identity drifted or retains placeholders",
    )
    _require(
        type(final["source_commit"]) is str
        and _COMMIT.fullmatch(final["source_commit"]) is not None,
        "Q023 source commit is malformed",
    )
    for name in (
        "source_bundle_sha256",
        "source_archive_sha256",
        "clean_candidate_manifest_sha256",
        "executable_sha256",
        "installed_control_plane_version",
        "environment_profile_sha256",
        "job_script_sha256",
        "reconcile_q023_registered_execution_sha256",
        "q043_registered_matrix_sha256",
        "q043_registered_matrix_case_bindings_sha256",
        "q043_registered_dependency_sha256",
    ):
        _require(
            type(final[name]) is str and _SHA256.fullmatch(final[name]) is not None,
            f"Q023 final binding {name} is malformed",
        )
    for name in (
        "source_archive_path",
        "clean_candidate_manifest_path",
        "executable_path",
        "orion_installed_control_plane_root",
        "project_home_installed_control_plane_root",
        "environment_profile_path",
        "job_script_path",
        "reconcile_q023_registered_execution_path",
        "q043_registered_matrix_path",
    ):
        _require(
            type(final[name]) is str
            and PurePosixPath(final[name]).is_absolute()
            and PurePosixPath(final[name]).as_posix() == final[name]
            and ".." not in PurePosixPath(final[name]).parts,
            f"Q023 final binding {name} path is malformed",
        )
    version = str(final["installed_control_plane_version"])
    orion_controller = AUTHORIZED_ORION_ROOT / "control_plane" / version
    project_controller = CANONICAL_PROJECT_HOME_ROOT / "control_plane" / version
    _require(
        final["orion_installed_control_plane_root"] == str(orion_controller)
        and final["project_home_installed_control_plane_root"] == str(project_controller)
        and final["environment_profile_path"]
        == str(orion_controller / "frontier_pic_environment.sh")
        and final["job_script_path"] == str(orion_controller / "frontier_job.sh")
        and final["reconcile_q023_registered_execution_path"]
        == str(orion_controller / "reconcile_q023_registered_execution.py"),
        "Q023 paired installed-control-plane path family drifted",
    )
    _require(
        final["analysis_script_paths"]
        == [final["reconcile_q023_registered_execution_path"]]
        and final["analysis_script_sha256"]
        == [final["reconcile_q023_registered_execution_sha256"]],
        "Q023 analysis support must be the immutable installed reconciler",
    )
    candidate_manifest = PurePosixPath(str(final["clean_candidate_manifest_path"]))
    executable = PurePosixPath(str(final["executable_path"]))
    source_archive = PurePosixPath(str(final["source_archive_path"]))
    _require(
        candidate_manifest.parent.parent
        == PurePosixPath(str(AUTHORIZED_ORION_ROOT / "clean_candidates"))
        and candidate_manifest.name == "clean_candidate_manifest.json"
        and executable.parent == candidate_manifest.parent
        and source_archive.parent == candidate_manifest.parent,
        "Q023 final clean-candidate path family drifted",
    )
    _require(
        final["q043_registered_matrix_record_type"]
        == "q043_registered_execution_raw_oracle_matrix_qualification",
        "Q023 final Q043 matrix record type drifted",
    )
    return final


def validate_final_binding_files(value: object) -> dict[str, object]:
    """Verify immutable files, paired controller members, and exact Q043 matrix."""
    final = validate_final_bindings(value)
    for path_name, digest_name, executable in (
        ("source_archive_path", "source_archive_sha256", False),
        ("clean_candidate_manifest_path", "clean_candidate_manifest_sha256", False),
        ("executable_path", "executable_sha256", True),
        ("environment_profile_path", "environment_profile_sha256", False),
        ("job_script_path", "job_script_sha256", False),
        (
            "reconcile_q023_registered_execution_path",
            "reconcile_q023_registered_execution_sha256",
            False,
        ),
        ("q043_registered_matrix_path", "q043_registered_matrix_sha256", False),
    ):
        _stable_binding(
            final[path_name],
            final[digest_name],
            label=path_name,
            executable=executable,
        )
    orion = Path(str(final["orion_installed_control_plane_root"]))
    project = Path(str(final["project_home_installed_control_plane_root"]))
    _require(orion.is_dir() and project.is_dir(), "Q023 installed controller pair is absent")
    for name, digest_name in (
        ("frontier_pic_environment.sh", "environment_profile_sha256"),
        ("frontier_job.sh", "job_script_sha256"),
        (
            "reconcile_q023_registered_execution.py",
            "reconcile_q023_registered_execution_sha256",
        ),
    ):
        _, orion_payload = q043._stable_regular_bytes(
            orion / name, label=f"Q023 Orion installed {name}"
        )
        _, project_payload = q043._stable_regular_bytes(
            project / name, label=f"Q023 Project Home installed {name}"
        )
        _require(
            orion_payload == project_payload
            and _sha256_bytes(orion_payload) == final[digest_name],
            f"Q023 paired installed {name} differs or drifted",
        )
    try:
        dependency = bell.registered_q043_raw_oracle_dependency(
            Path(str(final["q043_registered_matrix_path"])),
            artifact_root=AUTHORIZED_ORION_ROOT,
        )
    except bell.ContractError as error:
        raise PreparationError("Q023 final binding lacks a valid Q043 matrix") from error
    _require(
        dependency["registered_matrix_sha256"]
        == final["q043_registered_matrix_sha256"],
        "Q023 final Q043 matrix digest drifted",
    )
    _require(
        dependency["registered_matrix_record_type"]
        == final["q043_registered_matrix_record_type"]
        and dependency["registered_matrix_case_bindings_sha256"]
        == final["q043_registered_matrix_case_bindings_sha256"],
        "Q023 final Q043 matrix schema or case bindings drifted",
    )
    _require(
        bell._dependency_digest(dependency)
        == final["q043_registered_dependency_sha256"],
        "Q023 canonical Q043 dependency digest drifted",
    )
    return final


def materialize_q023_final_bindings(
    *,
    clean_candidate_manifest: Path,
    clean_candidate_manifest_sha256: str,
    expected_source_commit: str,
    installed_control_plane_version: str,
    q043_registered_matrix: Path,
) -> dict[str, object]:
    """Derive one exact Q023 launch-binding record from immutable evidence."""
    _require(
        _SHA256.fullmatch(clean_candidate_manifest_sha256) is not None,
        "Q023 clean-candidate manifest digest is malformed",
    )
    _require(
        _COMMIT.fullmatch(expected_source_commit) is not None,
        "Q023 expected source commit is malformed",
    )
    _require(
        _SHA256.fullmatch(installed_control_plane_version) is not None,
        "Q023 installed control-plane version is malformed",
    )
    report = revalidate_clean_candidate(
        clean_candidate_manifest,
        expected_manifest_sha256=clean_candidate_manifest_sha256,
        expected_git_commit=expected_source_commit,
        expected_receipt_control_plane_version=installed_control_plane_version,
        control_plane_dir=CONTROL_PLANE_SOURCE_DIR,
        authorized_pic_root=AUTHORIZED_ORION_ROOT,
        authorized_project_home_root=CANONICAL_PROJECT_HOME_ROOT,
    )
    _require(
        report.get("status") == "passed"
        and report.get("current_control_plane_version")
        == installed_control_plane_version
        and isinstance(report.get("build"), dict)
        and report["build"].get("receipt_control_plane_version")
        == installed_control_plane_version
        and isinstance(report.get("source"), dict)
        and report["source"].get("git_commit") == expected_source_commit,
        "Q023 clean-candidate revalidation report drifted",
    )

    _, manifest_payload = q043._stable_regular_bytes(
        clean_candidate_manifest,
        label="Q023 clean-candidate manifest",
    )
    _require(
        _sha256_bytes(manifest_payload) == clean_candidate_manifest_sha256,
        "Q023 clean-candidate manifest digest drifted",
    )
    try:
        candidate = json.loads(manifest_payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PreparationError(
            "Q023 clean-candidate manifest is not valid JSON"
        ) from error
    _require(
        type(candidate) is dict
        and type(candidate.get("source")) is dict
        and type(candidate.get("build")) is dict,
        "Q023 clean-candidate manifest structure drifted",
    )
    source = candidate["source"]
    build = candidate["build"]
    _require(
        source.get("git_commit") == expected_source_commit,
        "Q023 clean-candidate source commit drifted",
    )
    candidate_dir = clean_candidate_manifest.parent
    _require(
        clean_candidate_manifest
        == AUTHORIZED_ORION_ROOT
        / "clean_candidates"
        / candidate_dir.name
        / "clean_candidate_manifest.json",
        "Q023 clean-candidate manifest path is not canonical",
    )
    source_archive = Path(str(source.get("archive_path", "")))
    executable = Path(str(build.get("executable_path", "")))
    _require(
        source_archive == candidate_dir / "source.tar"
        and executable == candidate_dir / "athena",
        "Q023 clean-candidate artifact path family drifted",
    )
    source_archive_sha256 = str(source.get("archive_sha256", ""))
    source_bundle_sha256 = str(source.get("source_bundle_sha256", ""))
    executable_sha256 = str(build.get("executable_sha256", ""))
    _stable_binding(
        str(source_archive),
        source_archive_sha256,
        label="Q023 source archive",
    )
    _stable_binding(
        str(executable),
        executable_sha256,
        label="Q023 executable",
        executable=True,
    )

    orion_controller = (
        AUTHORIZED_ORION_ROOT / "control_plane" / installed_control_plane_version
    )
    project_controller = (
        CANONICAL_PROJECT_HOME_ROOT
        / "control_plane"
        / installed_control_plane_version
    )
    controller_bindings: dict[str, tuple[str, str]] = {}
    for name in (
        "frontier_pic_environment.sh",
        "frontier_job.sh",
        "reconcile_q023_registered_execution.py",
    ):
        require_executable = name.endswith(".sh")
        _, orion_payload = q043._stable_regular_bytes(
            orion_controller / name,
            label=f"Q023 Orion installed {name}",
            require_executable=require_executable,
        )
        _, project_payload = q043._stable_regular_bytes(
            project_controller / name,
            label=f"Q023 Project Home installed {name}",
            require_executable=require_executable,
        )
        _require(
            orion_payload == project_payload,
            f"Q023 paired installed {name} differs",
        )
        controller_bindings[name] = (
            str(orion_controller / name),
            _sha256_bytes(orion_payload),
        )

    try:
        dependency = bell.registered_q043_raw_oracle_dependency(
            q043_registered_matrix,
            artifact_root=AUTHORIZED_ORION_ROOT,
        )
    except bell.ContractError as error:
        raise PreparationError("Q023 final binding lacks a valid Q043 matrix") from error
    _require(
        dependency["binding_kind"] == "registered_matrix_qualification",
        "Q023 final binding requires a registered Q043 matrix",
    )
    environment_path, environment_sha256 = controller_bindings[
        "frontier_pic_environment.sh"
    ]
    job_path, job_sha256 = controller_bindings["frontier_job.sh"]
    reconciler_path, reconciler_sha256 = controller_bindings[
        "reconcile_q023_registered_execution.py"
    ]
    final = {
        "record_type": FINAL_BINDING_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "source_commit": expected_source_commit,
        "source_bundle_sha256": source_bundle_sha256,
        "source_archive_path": str(source_archive),
        "source_archive_sha256": source_archive_sha256,
        "clean_candidate_manifest_path": str(clean_candidate_manifest),
        "clean_candidate_manifest_sha256": clean_candidate_manifest_sha256,
        "executable_path": str(executable),
        "executable_sha256": executable_sha256,
        "installed_control_plane_version": installed_control_plane_version,
        "orion_installed_control_plane_root": str(orion_controller),
        "project_home_installed_control_plane_root": str(project_controller),
        "environment_profile_path": environment_path,
        "environment_profile_sha256": environment_sha256,
        "job_script_path": job_path,
        "job_script_sha256": job_sha256,
        "analysis_script_paths": [reconciler_path],
        "analysis_script_sha256": [reconciler_sha256],
        "reconcile_q023_registered_execution_path": reconciler_path,
        "reconcile_q023_registered_execution_sha256": reconciler_sha256,
        "q043_registered_matrix_path": str(q043_registered_matrix),
        "q043_registered_matrix_sha256": dependency[
            "registered_matrix_sha256"
        ],
        "q043_registered_matrix_record_type": dependency[
            "registered_matrix_record_type"
        ],
        "q043_registered_matrix_case_bindings_sha256": dependency[
            "registered_matrix_case_bindings_sha256"
        ],
        "q043_registered_dependency_sha256": bell._dependency_digest(dependency),
    }
    return validate_final_binding_files(final)


def _selected_bindings(
    final_bindings: Mapping[str, object] | None,
) -> tuple[str, dict[str, object], list[str]]:
    if final_bindings is None:
        return (
            "pending_final_digests_launch_prohibited",
            dict(PENDING_FINAL_BINDINGS),
            sorted(FINAL_BINDING_KEYS - {"record_type", "schema_version"}),
        )
    final = validate_final_bindings(final_bindings)
    return "final_digests_bound_review_candidate_still_blocked", final, []


def _members(
    final_bindings: Mapping[str, object] | None = None,
) -> list[dict[str, object]]:
    dependency = None
    artifact_root = None
    measured_manifest = json.loads(bell.DECK_MANIFEST.read_text(encoding="utf-8"))
    registered_manifest = (
        measured_manifest.get("foundational_registered_admission_binding_status")
        == bell.Q043_REGISTERED_MATRIX_BINDING_STATUS
    )
    if final_bindings is not None and registered_manifest:
        artifact_root = AUTHORIZED_ORION_ROOT
        dependency = bell.registered_q043_raw_oracle_dependency(
            Path(str(final_bindings["q043_registered_matrix_path"])),
            artifact_root=artifact_root,
        )
    manifest = bell.validate_checked_in_decks(
        q043_registered_raw_oracle_dependency=dependency,
        q043_artifact_root=artifact_root,
    )
    cases = [dict(case) for case in manifest["cases"]]
    _require(
        len(cases) == EXPECTED_CASE_COUNT
        and len({case["member_id"] for case in cases}) == EXPECTED_CASE_COUNT,
        "Q023 checked-in matrix is incomplete or duplicated",
    )
    expected = [dict(case) for case in bell.expected_deck_members()]
    projected = [
        {key: case[key] for key in expected[0]}
        if set(expected[0]) <= set(case)
        else None
        for case in cases
    ]
    _require(
        _strict_equal(projected, expected),
        "Q023 checked-in matrix differs from the exact ordered 55-case contract",
    )
    if final_bindings is not None:
        required = {
            "foundational_registered_admission_binding_status": (
                bell.Q043_REGISTERED_MATRIX_BINDING_STATUS
            ),
            "foundational_registered_admission_digest_bound": True,
            "foundational_registered_admission_schema_bound": True,
            "foundational_registered_matrix_sha256": final_bindings[
                "q043_registered_matrix_sha256"
            ],
            "foundational_registered_matrix_record_type": final_bindings[
                "q043_registered_matrix_record_type"
            ],
            "foundational_registered_matrix_case_bindings_sha256": final_bindings[
                "q043_registered_matrix_case_bindings_sha256"
            ],
            "foundational_registered_dependency_sha256": final_bindings[
                "q043_registered_dependency_sha256"
            ],
        }
        _require(
            all(manifest.get(key) == expected_value for key, expected_value in required.items()),
            "Q023 final materialization requires decks rebound to the exact Q043 matrix",
        )
    return cases


def _rank_count(member: Mapping[str, object]) -> int:
    return math.prod(int(value) for value in member["decomposition_splits"])


def _authorization_id(index: int) -> str:
    return f"q023-linear-{index:03d}-v1"


def _case_root_template(member_id: str) -> str:
    _require(
        re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,127}", member_id) is not None,
        "Q023 member ID is unsafe",
    )
    return str(AUTHORIZED_ORION_ROOT / RUN_NAMESPACE / SUBMISSION_ID_TEMPLATE)


def _raw_root_template(member_id: str) -> str:
    return f"{_case_root_template(member_id)}/raw"


def _input_deck_snapshot_template(member_id: str) -> str:
    return str(
        AUTHORIZED_ORION_ROOT
        / "manifests"
        / CAMPAIGN
        / SUBMISSION_ID_TEMPLATE
        / "snapshot"
        / f"q023-joverc-{member_id}.athinput"
    )


def _orion_receipt_template(member_id: str) -> str:
    return (
        f"{_case_root_template(member_id)}/analysis/"
        "q023_registered_execution_receipt.json"
    )


def _project_home_receipt_template(member_id: str) -> str:
    _require(member_id, "Q023 member ID is empty")
    return str(
        CANONICAL_PROJECT_HOME_ROOT
        / PROJECT_HOME_RECEIPT_NAMESPACE
        / SUBMISSION_ID_TEMPLATE
        / "q023_registered_execution_receipt.json"
    )


def _raw_relative_path(member_id: str, field: str, output_index: int) -> str:
    basename = "q023_joverc_" + member_id.replace("-", "_")
    return f"bin/{basename}.{field}.{output_index:05d}.bin"


def _output_topology(
    member: Mapping[str, object], deck_text: str
) -> dict[str, object]:
    blocks = bell.parse_athinput_text(deck_text)
    for index, field in enumerate(REQUIRED_OUTPUT_FIELDS, 1):
        output = blocks.get(f"output{index}", {})
        _require(
            output.get("file_type") == "bin"
            and output.get("variable") == field
            and float(output.get("dt", "nan")) == bell.LINEAR_OUTPUT_DT
            and output.get("single_file_per_rank") == "false",
            f"{member['member_id']}: Q023 shared-MPI output topology drifted",
        )
    restart = blocks.get("output6", {})
    _require(
        restart.get("file_type") == "rst"
        and float(restart.get("dt", "nan")) == bell.LINEAR_RUNTIME_TLIM
        and restart.get("single_file_per_rank") == "false",
        f"{member['member_id']}: Q023 restart topology drifted",
    )
    member_id = str(member["member_id"])
    inventory = [
        _raw_relative_path(member_id, field, output_index)
        for output_index in REQUIRED_OUTPUT_INDICES
        for field in REQUIRED_OUTPUT_FIELDS
    ]
    bounded_raw_bytes = len(REQUIRED_OUTPUT_INDICES) * sum(
        MAXIMUM_RAW_ARTIFACT_BYTES_BY_FIELD.values()
    )
    _require(
        bounded_raw_bytes + MAXIMUM_NON_RAW_CASE_BYTES
        <= MAXIMUM_CASE_STORAGE_BYTES,
        f"{member_id}: Q023 declared raw topology exceeds case storage ceiling",
    )
    return {
        "shared_mpi_io": True,
        "single_file_per_rank": False,
        "required_output_indices": list(REQUIRED_OUTPUT_INDICES),
        "required_fields": list(REQUIRED_OUTPUT_FIELDS),
        "expected_raw_artifact_count": len(inventory),
        "expected_relative_raw_paths": inventory,
        "maximum_bytes_per_raw_artifact_by_field": dict(
            MAXIMUM_RAW_ARTIFACT_BYTES_BY_FIELD
        ),
        "maximum_bounded_raw_bytes": bounded_raw_bytes,
        "maximum_non_raw_case_bytes": MAXIMUM_NON_RAW_CASE_BYTES,
        "maximum_case_storage_bytes": MAXIMUM_CASE_STORAGE_BYTES,
        "restart_output_sealed_by_artifact_inventory": True,
    }


def _expected_command(
    member: Mapping[str, object], final_bindings: Mapping[str, object]
) -> list[str]:
    ranks = _rank_count(member)
    member_id = str(member["member_id"])
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
        _input_deck_snapshot_template(member_id),
        "-d",
        _raw_root_template(member_id),
    ]


def _wrapper_evidence(member: Mapping[str, object]) -> dict[str, object]:
    ranks = _rank_count(member)
    return {
        "required_exact_rank_line": (
            f"Q023_REGISTERED_EXECUTION case_id={member['member_id']} "
            f"mpi_world_size={ranks} "
            f"rank_ids={','.join(str(rank) for rank in range(ranks))}"
        ),
        "required_exact_exit_line": (
            "Q023_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0"
        ),
        "required_exact_termination_line": "Terminating on time limit",
        "required_terminal_time": bell.LINEAR_RUNTIME_TLIM,
        "required_nlim": bell.LINEAR_NLIM,
        "trusted_wrapper_must_emit_each_exact_line_once": True,
        "athena_stdout_alone_without_trusted_wrapper_evidence_sufficient": False,
    }


def _reconciliation_interface(
    member: Mapping[str, object], final_bindings: Mapping[str, object]
) -> dict[str, object]:
    member_id = str(member["member_id"])
    return {
        "producer_name": "reconcile_q023_registered_execution.py",
        "producer_path": final_bindings["reconcile_q023_registered_execution_path"],
        "producer_sha256": final_bindings[
            "reconcile_q023_registered_execution_sha256"
        ],
        "installed_generation_only": True,
        "required_receipt_record_type": (
            "q023_reconciled_registered_execution_receipt"
        ),
        "required_receipt_role": "immutable_reconciled_registered_execution",
        "orion_receipt_path_template": _orion_receipt_template(member_id),
        "canonical_project_home_receipt_path_template": (
            _project_home_receipt_template(member_id)
        ),
        "paired_receipts_required": True,
        "paired_receipts_must_be_byte_identical": True,
        "raw_inventory_captured_at_reconciliation": True,
        "required_raw_artifact_count": len(REQUIRED_OUTPUT_INDICES)
        * len(REQUIRED_OUTPUT_FIELDS),
    }


def _deck_binding(member: Mapping[str, object]) -> dict[str, object]:
    relative = (
        bell.DECK_ROOT.relative_to(REPO_ROOT) / str(member["deck_path"])
    ).as_posix()
    path = REPO_ROOT / relative
    _, payload = q043._stable_regular_bytes(
        path, label=f"Q023 checked-in deck {member['member_id']}", require_read_only=False
    )
    _require(
        _sha256_bytes(payload) == member["deck_sha256"],
        f"{member['member_id']}: deck digest drifted",
    )
    return {
        "path": relative,
        "sha256": member["deck_sha256"],
        "byte_count": len(payload),
    }


def _launch_contract(member: Mapping[str, object]) -> dict[str, object]:
    contract = {
        "schema_version": 1,
        "executor": TRUSTED_LAUNCH_EXECUTOR,
        "pre_actions": [],
        "actions": [
            {
                "action_id": member["member_id"],
                "kind": "athena",
                "resources": {
                    "nodes": 1,
                    "tasks": _rank_count(member),
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
        raise PreparationError("Q023 launch contract failed schema validation") from error


def build_materialization(
    final_bindings: Mapping[str, object] | None = None,
) -> tuple[dict[str, object], dict[str, bytes]]:
    stage, selected, unresolved = _selected_bindings(final_bindings)
    files: dict[str, bytes] = {}
    case_records: list[dict[str, object]] = []
    launches: list[dict[str, object]] = []
    policy_candidates: list[dict[str, object]] = []
    policy_slices: list[dict[str, object]] = []
    members = _members(selected if final_bindings is not None else None)
    for index, member in enumerate(members, 1):
        member_id = str(member["member_id"])
        deck = _deck_binding(member)
        deck_text = (REPO_ROOT / str(deck["path"])).read_text(encoding="utf-8")
        contract = _launch_contract(member)
        topology = _output_topology(member, deck_text)
        authorization_id = _authorization_id(index)
        launch = {
            "record_type": LAUNCH_RECORD_TYPE,
            "schema_version": SCHEMA_VERSION,
            "successor_id": SUCCESSOR_ID,
            "status": "source_local_review_candidate_incomplete_launch_prohibited",
            "qualification_effect": QUALIFICATION_EFFECT,
            "campaign": CAMPAIGN,
            "binding_stage": stage,
            "member_id": member_id,
            "authorization_id": authorization_id,
            "attempt_id": f"q023-linear-{index:03d}-{member_id}",
            "case_contract": member,
            "checked_in_deck": deck,
            "selected_final_bindings": selected,
            "unresolved_final_bindings": unresolved,
            "roots": {
                "authorized_orion_case_root_template": _case_root_template(
                    member_id
                ),
                "raw_output_root_template": _raw_root_template(member_id),
                "canonical_project_home_receipt_root_template": str(
                    CANONICAL_PROJECT_HOME_ROOT
                    / PROJECT_HOME_RECEIPT_NAMESPACE
                    / SUBMISSION_ID_TEMPLATE
                ),
            },
            "command_template": _expected_command(member, selected),
            "launch_contract": contract,
            "launch_contract_sha256": launch_contract_sha256(contract),
            "mpi_and_decomposition": {
                "mpi_ranks": _rank_count(member),
                "requested_nodes": MAXIMUM_NODES_PER_CASE,
                "requested_tasks": _rank_count(member),
                "decomposition_splits": member["decomposition_splits"],
                "global_nx": member["global_nx"],
                "meshblock_nx": member["meshblock_nx"],
                "one_meshblock_per_rank": True,
            },
            "output_topology": topology,
            "stdout_wrapper_evidence": _wrapper_evidence(member),
            "trusted_reconciliation_interface": _reconciliation_interface(
                member, selected
            ),
            "resource_ceiling": {
                "selected_qos": "normal",
                "registered_short_nonproduction": False,
                "maximum_nodes": MAXIMUM_NODES_PER_CASE,
                "maximum_walltime_seconds": MAXIMUM_WALLTIME_SECONDS_PER_CASE,
                "maximum_attempts": MAXIMUM_ATTEMPTS_PER_CASE,
                "maximum_retries": MAXIMUM_RETRIES_PER_CASE,
                "maximum_node_hours": 1.0,
                "maximum_storage_bytes": MAXIMUM_CASE_STORAGE_BYTES,
            },
            "required_blockers": list(REQUIRED_BLOCKERS),
            "authorization": dict(AUTHORIZATION_BOUNDARY),
        }
        policy_slice = {
            "authorization_id": authorization_id,
            "status": "review_required_not_authorized",
            "campaign": CAMPAIGN,
            "test_id": member_id,
            "evidence_class": EVIDENCE_CLASS,
            "physical_mode": PHYSICAL_MODE,
            "runtime_profile": RUNTIME_PROFILE,
            "selected_qos": "normal",
            "registered_short_nonproduction": False,
            "maximum_nodes": MAXIMUM_NODES_PER_CASE,
            "maximum_walltime_seconds": MAXIMUM_WALLTIME_SECONDS_PER_CASE,
            "maximum_attempts": MAXIMUM_ATTEMPTS_PER_CASE,
            "job_script_sha256": selected["job_script_sha256"],
            "input_deck_sha256": deck["sha256"],
            "environment_profile_sha256": selected["environment_profile_sha256"],
            "analysis_script_sha256": [
                selected["reconcile_q023_registered_execution_sha256"]
            ],
            "executable_sha256": selected["executable_sha256"],
            "launch_contract_sha256": launch["launch_contract_sha256"],
            "clean_candidate_manifest_sha256": selected[
                "clean_candidate_manifest_sha256"
            ],
        }
        policy_candidate = {
            "record_type": POLICY_RECORD_TYPE,
            "schema_version": SCHEMA_VERSION,
            "successor_id": SUCCESSOR_ID,
            "status": "review_required_not_authorized_not_live_policy",
            "qualification_effect": QUALIFICATION_EFFECT,
            "campaign": CAMPAIGN,
            "binding_stage": stage,
            "member_id": member_id,
            "policy_slice_candidate": policy_slice,
            "launch_candidate_sha256": _sha256_bytes(_json_bytes(launch)),
            "unresolved_final_bindings": unresolved,
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
        launch_path = f"launch_candidates/{member_id}.json"
        policy_path = f"policy_candidates/{member_id}.json"
        launch_payload = _json_bytes(launch)
        _require(
            policy_candidate["launch_candidate_sha256"]
            == _sha256_bytes(launch_payload),
            f"{member_id}: Q023 launch-candidate digest drifted",
        )
        policy_payload = _json_bytes(policy_candidate)
        files[launch_path] = launch_payload
        files[policy_path] = policy_payload
        launches.append(launch)
        policy_candidates.append(policy_candidate)
        policy_slices.append(policy_slice)
        case_records.append(
            {
                "case_index": index,
                "member_id": member_id,
                "dimension": member["dimension"],
                "epsilon": member["epsilon"],
                "resolution": member["resolution"],
                "decomposition": member["decomposition"],
                "mpi_ranks": _rank_count(member),
                "checked_in_deck": deck,
                "launch_candidate": _binding(launch_path, launch_payload),
                "policy_candidate": _binding(policy_path, policy_payload),
            }
        )
    budget = validate_budget_accounting(
        {
        "case_count": EXPECTED_CASE_COUNT,
        "maximum_live_q023_submissions": MAXIMUM_LIVE_Q023_SUBMISSIONS,
        "maximum_attempts_per_case": MAXIMUM_ATTEMPTS_PER_CASE,
        "maximum_retries_per_case": MAXIMUM_RETRIES_PER_CASE,
        "maximum_nodes_per_case": MAXIMUM_NODES_PER_CASE,
        "maximum_walltime_seconds_per_case": MAXIMUM_WALLTIME_SECONDS_PER_CASE,
        "maximum_registered_attempts": EXPECTED_CASE_COUNT,
        "maximum_batch_node_hours": MAXIMUM_BATCH_NODE_HOURS,
        "authorized_total_node_hour_cap": MAXIMUM_BATCH_NODE_HOURS,
        "minimum_live_ledger_remaining_node_hours_required": (
            MAXIMUM_BATCH_NODE_HOURS
        ),
        "live_ledger_remaining_node_hours": (
            "fresh_at_policy_promotion_not_assumed"
        ),
        "maximum_batch_storage_bytes": (
            EXPECTED_CASE_COUNT * MAXIMUM_CASE_STORAGE_BYTES
        ),
        "maximum_storage_cap_bytes": MAXIMUM_BATCH_STORAGE_BYTES,
        "fresh_storage_preflight_required": True,
        "ledger_mutation_authorized": False,
        },
        launches,
    )
    policy_fragment = {
        "record_type": POLICY_FRAGMENT_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "status": "aggregate_review_fragment_only_not_live_policy",
        "qualification_effect": QUALIFICATION_EFFECT,
        "campaign": CAMPAIGN,
        "binding_stage": stage,
        "case_count": EXPECTED_CASE_COUNT,
        "registered_science_slice_candidates": policy_slices,
        "declared_batch_ceiling": budget,
        "unresolved_final_bindings": unresolved,
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
    reconciliation = {
        "producer_interface": "reconcile_q023_registered_execution.py",
        "producer_path": selected["reconcile_q023_registered_execution_path"],
        "producer_sha256": selected[
            "reconcile_q023_registered_execution_sha256"
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
        "status": "source_local_preparation_complete_launch_and_policy_blocked",
        "qualification_effect": QUALIFICATION_EFFECT,
        "campaign": CAMPAIGN,
        "binding_stage": stage,
        "selected_final_bindings": selected,
        "unresolved_final_bindings": unresolved,
        "case_count": len(case_records),
        "case_records": case_records,
        "policy_slices": policy_slices,
        "launch_candidate_count": len(launches),
        "policy_candidate_count": len(policy_candidates),
        "batch_budget_accounting_input": _binding(
            "batch_budget_accounting_input.json", budget_payload
        ),
        "registered_policy_slice_review_fragment": _binding(
            "registered_policy_slice_review_fragment.json",
            policy_fragment_payload,
        ),
        "trusted_reconciliation_producer_interface": _binding(
            "trusted_reconciliation_producer_interface.json",
            reconciliation_payload,
        ),
        "exact_matrix_contract": {
            "case_count": EXPECTED_CASE_COUNT,
            "total_mpi_ranks": sum(_rank_count(member) for member in members),
            "rank_distribution": {"2": 45, "4": 5, "8": 5},
            "raw_artifacts_per_case": len(REQUIRED_OUTPUT_INDICES)
            * len(REQUIRED_OUTPUT_FIELDS),
            "raw_artifacts_total": (
                EXPECTED_CASE_COUNT
                * len(REQUIRED_OUTPUT_INDICES)
                * len(REQUIRED_OUTPUT_FIELDS)
            ),
            "one_registered_attempt_per_case": True,
            "one_meshblock_per_rank": True,
            "shared_mpi_output": True,
        },
        "required_blockers": list(REQUIRED_BLOCKERS),
        "execution_boundary": {
            "exact_registered_q043_matrix_required": True,
            "exact_final_clean_candidate_required": True,
            "exact_installed_control_plane_required": True,
            "empty_user_queue_required_for_every_submission": True,
            **AUTHORIZATION_BOUNDARY,
        },
    }
    _require(
        stage != "final_digests_bound_review_candidate_still_blocked"
        or not _contains_placeholder(manifest),
        "Q023 final-digest materialization retained a placeholder",
    )
    return manifest, files


def validate_budget_accounting(
    value: object, launches: Sequence[Mapping[str, object]]
) -> dict[str, object]:
    _require(type(value) is dict, "Q023 budget accounting must be an object")
    expected_node_hours = (
        len(launches)
        * MAXIMUM_NODES_PER_CASE
        * MAXIMUM_WALLTIME_SECONDS_PER_CASE
        * MAXIMUM_ATTEMPTS_PER_CASE
        / 3600.0
    )
    expected_storage = len(launches) * MAXIMUM_CASE_STORAGE_BYTES
    _require(
        len(launches) == EXPECTED_CASE_COUNT
        and value.get("case_count") == EXPECTED_CASE_COUNT
        and value.get("maximum_registered_attempts") == EXPECTED_CASE_COUNT
        and value.get("maximum_live_q023_submissions")
        == MAXIMUM_LIVE_Q023_SUBMISSIONS
        and value.get("maximum_attempts_per_case") == MAXIMUM_ATTEMPTS_PER_CASE
        and value.get("maximum_retries_per_case") == MAXIMUM_RETRIES_PER_CASE,
        "Q023 budget attempt or concurrency ceiling drifted",
    )
    for launch in launches:
        ceiling = launch.get("resource_ceiling")
        _require(
            type(ceiling) is dict
            and ceiling.get("maximum_nodes") == MAXIMUM_NODES_PER_CASE
            and ceiling.get("maximum_walltime_seconds")
            == MAXIMUM_WALLTIME_SECONDS_PER_CASE
            and ceiling.get("maximum_attempts") == MAXIMUM_ATTEMPTS_PER_CASE
            and ceiling.get("maximum_retries") == MAXIMUM_RETRIES_PER_CASE
            and ceiling.get("maximum_node_hours") == 1.0
            and ceiling.get("maximum_storage_bytes")
            == MAXIMUM_CASE_STORAGE_BYTES,
            "Q023 per-case resource ceiling drifted",
        )
    _require(
        value.get("maximum_batch_node_hours") == expected_node_hours
        and expected_node_hours == MAXIMUM_BATCH_NODE_HOURS,
        "Q023 batch node-hour ceiling drifted",
    )
    _require(
        value.get("maximum_batch_storage_bytes") == expected_storage
        and value.get("maximum_storage_cap_bytes")
        == MAXIMUM_BATCH_STORAGE_BYTES
        and expected_storage <= MAXIMUM_BATCH_STORAGE_BYTES,
        "Q023 batch storage ceiling overflow or drift",
    )
    _require(
        value.get("live_ledger_remaining_node_hours")
        == "fresh_at_policy_promotion_not_assumed"
        and value.get("fresh_storage_preflight_required") is True
        and value.get("ledger_mutation_authorized") is False,
        "Q023 budget accounting assumed live authority",
    )
    return dict(value)


def validate_materialization(
    manifest: object,
    files: Mapping[str, bytes],
    final_bindings: Mapping[str, object] | None = None,
) -> tuple[dict[str, object], dict[str, bytes]]:
    expected_manifest, expected_files = build_materialization(final_bindings)
    _require(
        _strict_equal(manifest, expected_manifest),
        "Q023 preparation manifest drifted",
    )
    _require(set(files) == set(expected_files), "Q023 preparation files drifted")
    for path, payload in expected_files.items():
        _require(
            type(files[path]) is bytes and files[path] == payload,
            f"Q023 preparation file drifted: {path}",
        )
    launches = []
    for record in expected_manifest["case_records"]:
        launch = json.loads(files[record["launch_candidate"]["path"]])
        policy = json.loads(files[record["policy_candidate"]["path"]])
        _require(
            validate_launch_contract(launch["launch_contract"])
            == launch["launch_contract"],
            f"{record['member_id']}: Q023 launch-contract schema drifted",
        )
        _require(
            launch["authorization"] == AUTHORIZATION_BOUNDARY
            and policy["authorization"] == AUTHORIZATION_BOUNDARY,
            f"{record['member_id']}: Q023 candidate acquired authority",
        )
        launches.append(launch)
    validate_budget_accounting(
        json.loads(files["batch_budget_accounting_input.json"]), launches
    )
    return expected_manifest, expected_files


def materialize_q023_timeout_margin(
    *,
    final_bindings: Mapping[str, object],
    measured_utc: str,
    expires_utc: str,
) -> dict[str, object]:
    final = validate_final_bindings(final_bindings)
    measured = q043._utc_datetime(measured_utc, label="Q023 timeout measured UTC")
    expires = q043._utc_datetime(expires_utc, label="Q023 timeout expires UTC")
    _require(
        measured < expires and (expires - measured).total_seconds() <= 24 * 3600,
        "Q023 timeout validity interval is invalid",
    )
    return {
        "athena_walltime_seconds": ATHENA_WALLTIME_SECONDS,
        "scheduler_walltime_seconds": MAXIMUM_WALLTIME_SECONDS_PER_CASE,
        "environment_profile_sha256": final["environment_profile_sha256"],
        "measured_utc": measured_utc,
        "expires_utc": expires_utc,
    }


def validate_q023_timeout_margin_artifact(
    path: Path,
    *,
    final_bindings: Mapping[str, object],
    now: datetime | None = None,
) -> str:
    final = validate_final_bindings(final_bindings)
    path, payload = q043._stable_regular_bytes(
        path, label="Q023 timeout-margin artifact"
    )
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PreparationError(
            "Q023 timeout-margin artifact is not UTF-8 JSON"
        ) from error
    required = {
        "athena_walltime_seconds",
        "scheduler_walltime_seconds",
        "environment_profile_sha256",
        "measured_utc",
        "expires_utc",
    }
    _require(
        type(value) is dict and set(value) == required,
        "Q023 timeout-margin artifact schema drifted",
    )
    _require(
        _strict_equal(
            value,
            materialize_q023_timeout_margin(
                final_bindings=final,
                measured_utc=value["measured_utc"],
                expires_utc=value["expires_utc"],
            ),
        ),
        "Q023 timeout-margin artifact walltime or environment drifted",
    )
    measured = q043._utc_datetime(
        value["measured_utc"], label="Q023 timeout measured UTC"
    )
    expires = q043._utc_datetime(
        value["expires_utc"], label="Q023 timeout expires UTC"
    )
    current = datetime.now(timezone.utc) if now is None else now
    _require(
        measured <= current < expires,
        "Q023 timeout-margin artifact is stale or premature",
    )
    return str(path)


def materialize_q023_pre_submit_config(
    *,
    member_id: str,
    submission_id: str,
    final_bindings: Mapping[str, object],
    pre_manifest_attestation: Path,
    timeout_margin_artifact: Path,
    queue_snapshot: Path,
    site_policy_checked_utc: str,
    now: datetime | None = None,
) -> dict[str, object]:
    final = validate_final_binding_files(final_bindings)
    members = _members(final)
    matches = [
        (index, member)
        for index, member in enumerate(members, 1)
        if member["member_id"] == member_id
    ]
    _require(len(matches) == 1, "Q023 selected member is unknown")
    index, member = matches[0]
    _require(
        _SUBMISSION.fullmatch(submission_id) is not None,
        "Q023 submission ID is malformed",
    )
    current = datetime.now(timezone.utc) if now is None else now
    checked = q043._utc_datetime(
        site_policy_checked_utc, label="Q023 site-policy checked UTC"
    )
    _require(
        0 <= (current - checked).total_seconds() <= 24 * 3600,
        "Q023 site-policy check is stale or in the future",
    )
    timeout = validate_q023_timeout_margin_artifact(
        timeout_margin_artifact, final_bindings=final, now=current
    )
    queue = q043._validate_empty_queue_snapshot(queue_snapshot)
    attestation = q043._validate_pre_manifest_attestation_path(
        pre_manifest_attestation
    )
    deck = _deck_binding(member)
    return {
        "pic_root": str(AUTHORIZED_ORION_ROOT),
        "campaign": CAMPAIGN,
        "test_id": member_id,
        "submission_scope": "registered_science",
        "registered_science_authorization_id": _authorization_id(index),
        "pre_manifest_attestation": attestation,
        "submission_id": submission_id,
        "git_commit": final["source_commit"],
        "evidence_class": EVIDENCE_CLASS,
        "physical_mode": PHYSICAL_MODE,
        "selected_qos": "normal",
        "qos_selection_reason": "normal_required_by_registered_campaign",
        "site_policy_checked_utc": site_policy_checked_utc,
        "registered_short_nonproduction": False,
        "artifact_dir": str(AUTHORIZED_ORION_ROOT / RUN_NAMESPACE / submission_id),
        "job_script_executable_env": "PIC_EXECUTABLE",
        "job_script": final["job_script_path"],
        "executable": final["executable_path"],
        "input_deck": str(REPO_ROOT / str(deck["path"])),
        "environment_profile": final["environment_profile_path"],
        "timeout_margin_artifact": timeout,
        "analysis_scripts": [final["reconcile_q023_registered_execution_path"]],
        "queue_snapshot": queue,
        "prior_case_closures": [],
        "clean_candidate_manifest": final["clean_candidate_manifest_path"],
        "launch_contract": _launch_contract(member),
    }


def materialize_q023_promotable_policy(
    *,
    baseline_policy: Mapping[str, object],
    final_bindings: Mapping[str, object],
) -> dict[str, object]:
    final = validate_final_binding_files(final_bindings)
    _require(
        isinstance(baseline_policy, Mapping)
        and _STORAGE_POLICY_REQUIRED_KEYS <= set(baseline_policy)
        and set(baseline_policy) <= _STORAGE_POLICY_ALLOWED_KEYS
        and baseline_policy.get("schema_version") == SCHEMA_VERSION,
        "Q023 baseline policy schema drifted",
    )
    _require(
        baseline_policy.get("registered_science_slices") == [],
        "Q023 promotable policy requires an empty registered allowlist",
    )
    _require(
        baseline_policy.get("frontier_admission_smoke")
        == {"status": "closed_after_pass"},
        "Q023 promotable policy requires closed Frontier admission smoke",
    )
    storage = baseline_policy.get("olcf_side_storage")
    _require(
        type(storage) is dict
        and storage.get("installed_control_plane_version")
        == final["installed_control_plane_version"]
        and storage.get("staged_control_plane_candidate_version")
        == final["installed_control_plane_version"]
        and storage.get("installed_control_plane_lifecycle")
        == "paired_installed_reviewed_generation",
        "Q023 baseline policy does not bind the selected controller generation",
    )
    try:
        validate_storage_policy(
            copy.deepcopy(dict(baseline_policy)),
            control_plane_version=str(final["installed_control_plane_version"]),
            authorized_pic_root=AUTHORIZED_ORION_ROOT,
            authorized_project_home_root=CANONICAL_PROJECT_HOME_ROOT,
        )
    except ValueError as error:
        raise PreparationError("Q023 baseline storage policy is invalid") from error
    manifest, _ = build_materialization(final)
    successor: dict[str, Any] = copy.deepcopy(dict(baseline_policy))
    successor["science_submission_freeze"] = {
        "status": "authorized",
        "manifest_path": final["clean_candidate_manifest_path"],
        "manifest_sha256": final["clean_candidate_manifest_sha256"],
        "build_profile_control_plane_version": final["installed_control_plane_version"],
    }
    successor["registered_science_slices"] = [
        {**slice_, "status": "authorized"} for slice_ in manifest["policy_slices"]
    ]
    try:
        validate_storage_policy(
            copy.deepcopy(successor),
            control_plane_version=str(final["installed_control_plane_version"]),
            authorized_pic_root=AUTHORIZED_ORION_ROOT,
            authorized_project_home_root=CANONICAL_PROJECT_HOME_ROOT,
        )
    except ValueError as error:
        raise PreparationError("Q023 successor policy failed validation") from error
    return successor


# Compatibility aliases retained for source-local callers created during preparation.
materialize_timeout_margin = materialize_q023_timeout_margin
materialize_pre_submit_config = materialize_q023_pre_submit_config
materialize_promotable_policy = materialize_q023_promotable_policy


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--materialize-final-bindings", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--clean-candidate-manifest", type=Path)
    parser.add_argument("--clean-candidate-manifest-sha256")
    parser.add_argument("--expected-source-commit")
    parser.add_argument("--installed-control-plane-version")
    parser.add_argument("--q043-registered-matrix", type=Path)
    return parser


def main() -> None:
    args = _parser().parse_args()
    if args.materialize_final_bindings:
        required = {
            "--output": args.output,
            "--clean-candidate-manifest": args.clean_candidate_manifest,
            "--clean-candidate-manifest-sha256": (
                args.clean_candidate_manifest_sha256
            ),
            "--expected-source-commit": args.expected_source_commit,
            "--installed-control-plane-version": (
                args.installed_control_plane_version
            ),
            "--q043-registered-matrix": args.q043_registered_matrix,
        }
        missing = [name for name, value in required.items() if value is None]
        _require(not missing, f"Q023 final binding arguments missing: {missing}")
        result = materialize_q023_final_bindings(
            clean_candidate_manifest=args.clean_candidate_manifest,
            clean_candidate_manifest_sha256=(
                args.clean_candidate_manifest_sha256
            ),
            expected_source_commit=args.expected_source_commit,
            installed_control_plane_version=args.installed_control_plane_version,
            q043_registered_matrix=args.q043_registered_matrix,
        )
        _write_json_exclusive(args.output, result)
    else:
        optional_values = (
            args.output,
            args.clean_candidate_manifest,
            args.clean_candidate_manifest_sha256,
            args.expected_source_commit,
            args.installed_control_plane_version,
            args.q043_registered_matrix,
        )
        _require(
            all(value is None for value in optional_values),
            "Q023 final binding arguments require --materialize-final-bindings",
        )
        manifest, files = build_materialization()
        result = {
            "manifest": manifest,
            "file_count": len(files),
            "manifest_sha256": _sha256_bytes(_json_bytes(manifest)),
        }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
