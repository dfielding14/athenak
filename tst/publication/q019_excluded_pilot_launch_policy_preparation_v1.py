#!/usr/bin/env python3
"""Prepare non-authorizing Q019 runtime-controller excluded pilots."""

from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
from typing import Any, Mapping, Sequence

from tst.publication import q019_nonlinear_bell_runtime_controller_v1 as controller
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as design
from tst.publication import (
    q023_registered_execution_linear_qualification_successor_v1 as q023,
)
from tst.publication import (
    q043_registered_execution_raw_oracle_qualification_successor_v1 as q043,
)
from tst.publication import (
    q043_registered_launch_policy_preparation_successor_v1 as q043_prep,
)


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
RUN_NAMESPACE = "runs/q019_nonlinear_bell_registered_successor_v1"
PROJECT_HOME_RECEIPT_NAMESPACE = "ledger/q019_registered_execution_receipts"
CAMPAIGN = "q019_nonlinear_bell_registered_successor_v1"
EVIDENCE_CLASS = "q019_excluded_resource_science_window_pilot"
PHYSICAL_MODE = "q019_finite_rigidity_nonlinear_bell_excluded_pilot"
RUNTIME_PROFILE = "frontier_minimum_supported"
SUCCESSOR_ID = "q019_excluded_pilot_launch_policy_preparation_v1"
SCHEMA_VERSION = 1
MANIFEST_RECORD_TYPE = "q019_excluded_pilot_launch_policy_preparation_manifest"
LAUNCH_RECORD_TYPE = "q019_excluded_pilot_launch_review_candidate"
POLICY_RECORD_TYPE = "q019_excluded_pilot_policy_slice_review_candidate"
POLICY_FRAGMENT_RECORD_TYPE = "q019_excluded_pilot_policy_slice_review_fragment"
FINAL_BINDING_RECORD_TYPE = "q019_excluded_pilot_launch_final_digest_bindings"
EXPECTED_ATTEMPT_COUNT = 4
MAXIMUM_ATTEMPTS_PER_PILOT = 1
MAXIMUM_RETRIES_PER_PILOT = 0
MAXIMUM_LIVE_Q019_SUBMISSIONS = 1
TASKS_PER_NODE = 8
MAXIMUM_BATCH_NODE_HOURS = 36.0
MAXIMUM_BATCH_STORAGE_BYTES = 640 * 1024**3
SUBMISSION_ID_TEMPLATE = "{submission_id}"
FINAL_PLACEHOLDER_PREFIX = "PENDING_FINAL_"
PILOT_RESOURCES = {
    2: {
        "nodes": 4,
        "tasks": 32,
        "scheduler_walltime_seconds": 1800,
        "athena_walltime_seconds": 1500,
        "maximum_storage_bytes": 64 * 1024**3,
    },
    3: {
        "nodes": 16,
        "tasks": 128,
        "scheduler_walltime_seconds": 3600,
        "athena_walltime_seconds": 3300,
        "maximum_storage_bytes": 256 * 1024**3,
    },
}
REQUIRED_OUTPUT_FIELDS = (
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
QUALIFICATION_EFFECT = (
    "resource_and_runtime_controller_pilot_only_excluded_from_q019_saturation_"
    "qualification_no_science_no_publication_authority"
)
AUTHORIZATION_BOUNDARY = {
    "launch_authorized": False,
    "scheduler_submission_authorized": False,
    "policy_mutation_authorized": False,
    "frontier_execution_authorized": False,
    "q019_qualification_authorized": False,
    "nonlinear_saturation_claim_authorized": False,
    "scientific_claim_authorized": False,
    "publication_authorized": False,
}
REQUIRED_BLOCKERS = (
    "exact_registered_q043_matrix_prerequisite_bound_and_revalidated",
    "exact_registered_q023_linear_matrix_prerequisite_bound_and_revalidated",
    "empty_user_queue_proved_at_fresh_reservation_boundary",
    "final_clean_candidate_manifest_source_archive_and_executable_bound",
    "final_paired_installed_control_plane_generation_bound_and_independently_verified",
    "installed_reconcile_q019_registered_execution_py_producer_bound",
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
    "reconcile_q019_registered_execution_path",
    "reconcile_q019_registered_execution_sha256",
    "q043_registered_matrix_path",
    "q043_registered_matrix_sha256",
    "q043_registered_matrix_record_type",
    "q043_registered_matrix_case_bindings_sha256",
    "q043_registered_dependency_sha256",
    "q023_registered_matrix_path",
    "q023_registered_matrix_sha256",
    "q023_registered_matrix_record_type",
    "q023_registered_matrix_case_bindings_sha256",
    "q023_registered_dependency_sha256",
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
_SAFE_ID = re.compile(r"[a-z0-9][a-z0-9_-]{0,127}")
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


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _binding(path: str, payload: bytes) -> dict[str, object]:
    return {"path": path, "sha256": _sha256_bytes(payload), "byte_count": len(payload)}


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


def _contains_placeholder(value: object) -> bool:
    if isinstance(value, str):
        return value.startswith(FINAL_PLACEHOLDER_PREFIX)
    if isinstance(value, Mapping):
        return any(_contains_placeholder(item) for item in value.values())
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return any(_contains_placeholder(item) for item in value)
    return False


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


def _stable_binding(
    path: object,
    digest: object,
    *,
    label: str,
    executable: bool = False,
) -> bytes:
    _require(
        isinstance(path, str)
        and Path(path).is_absolute()
        and isinstance(digest, str)
        and _SHA256.fullmatch(digest) is not None,
        f"{label} binding is malformed",
    )
    _, payload = q043_prep._stable_regular_bytes(
        Path(path), label=label, require_executable=executable
    )
    _require(_sha256_bytes(payload) == digest, f"{label} digest drifted")
    return payload


def validate_final_bindings(value: object) -> dict[str, object]:
    _require(
        type(value) is dict and set(value) == FINAL_BINDING_KEYS,
        "Q019 pilot final-binding keys drifted",
    )
    final = dict(value)
    _require(
        final["record_type"] == FINAL_BINDING_RECORD_TYPE
        and type(final["schema_version"]) is int
        and final["schema_version"] == SCHEMA_VERSION
        and not _contains_placeholder(final),
        "Q019 pilot final-binding identity drifted or retains placeholders",
    )
    _require(
        type(final["source_commit"]) is str
        and _COMMIT.fullmatch(final["source_commit"]) is not None,
        "Q019 pilot source commit is malformed",
    )
    digest_names = {
        key
        for key in FINAL_BINDING_KEYS
        if key.endswith("_sha256") and key not in {"analysis_script_sha256"}
    } | {"installed_control_plane_version"}
    for name in digest_names:
        _require(
            type(final[name]) is str and _SHA256.fullmatch(final[name]) is not None,
            f"Q019 pilot final binding {name} is malformed",
        )
    _require(
        type(final["analysis_script_paths"]) is list
        and type(final["analysis_script_sha256"]) is list
        and len(final["analysis_script_paths"]) == 1
        and len(final["analysis_script_sha256"]) == 1
        and _SHA256.fullmatch(final["analysis_script_sha256"][0]) is not None,
        "Q019 pilot analysis binding is malformed",
    )
    path_names = {
        key
        for key in FINAL_BINDING_KEYS
        if key.endswith("_path") or key.endswith("_root")
    } | {"analysis_script_paths"}
    for name in path_names:
        values = final[name] if isinstance(final[name], list) else [final[name]]
        _require(
            all(
                isinstance(item, str)
                and PurePosixPath(item).is_absolute()
                and PurePosixPath(item).as_posix() == item
                and ".." not in PurePosixPath(item).parts
                for item in values
            ),
            f"Q019 pilot final binding {name} path is malformed",
        )
    version = str(final["installed_control_plane_version"])
    orion = AUTHORIZED_ORION_ROOT / "control_plane" / version
    project = CANONICAL_PROJECT_HOME_ROOT / "control_plane" / version
    reconciler = orion / "reconcile_q019_registered_execution.py"
    _require(
        final["orion_installed_control_plane_root"] == str(orion)
        and final["project_home_installed_control_plane_root"] == str(project)
        and final["environment_profile_path"]
        == str(orion / "frontier_pic_environment.sh")
        and final["job_script_path"] == str(orion / "frontier_job.sh")
        and final["reconcile_q019_registered_execution_path"] == str(reconciler)
        and final["analysis_script_paths"] == [str(reconciler)]
        and final["analysis_script_sha256"]
        == [final["reconcile_q019_registered_execution_sha256"]],
        "Q019 pilot paired installed-control-plane path family drifted",
    )
    candidate_manifest = PurePosixPath(str(final["clean_candidate_manifest_path"]))
    _require(
        candidate_manifest.parent.parent
        == PurePosixPath(str(AUTHORIZED_ORION_ROOT / "clean_candidates"))
        and candidate_manifest.name == "clean_candidate_manifest.json"
        and PurePosixPath(str(final["executable_path"])).parent
        == candidate_manifest.parent
        and PurePosixPath(str(final["source_archive_path"])).parent
        == candidate_manifest.parent,
        "Q019 pilot clean-candidate path family drifted",
    )
    _require(
        final["q043_registered_matrix_record_type"]
        == q043.MATRIX_RECORD_TYPE
        and final["q023_registered_matrix_record_type"] == q023.MATRIX_RECORD_TYPE,
        "Q019 pilot predecessor matrix record type drifted",
    )
    return final


def _read_json_binding(path: str, digest: str, *, label: str) -> dict[str, object]:
    payload = _stable_binding(path, digest, label=label)
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PreparationError(f"{label} is not UTF-8 JSON") from error
    _require(type(value) is dict, f"{label} is not a JSON object")
    return value


def validate_final_binding_files(value: object) -> dict[str, object]:
    final = validate_final_bindings(value)
    for path_name, digest_name, executable in (
        ("source_archive_path", "source_archive_sha256", False),
        ("clean_candidate_manifest_path", "clean_candidate_manifest_sha256", False),
        ("executable_path", "executable_sha256", True),
        ("environment_profile_path", "environment_profile_sha256", False),
        ("job_script_path", "job_script_sha256", False),
        (
            "reconcile_q019_registered_execution_path",
            "reconcile_q019_registered_execution_sha256",
            False,
        ),
    ):
        _stable_binding(
            final[path_name],
            final[digest_name],
            label=path_name,
            executable=executable,
        )
    orion = Path(str(final["orion_installed_control_plane_root"]))
    project = Path(str(final["project_home_installed_control_plane_root"]))
    _require(orion.is_dir() and project.is_dir(), "Q019 pilot controller pair is absent")
    for name, digest_name in (
        ("frontier_pic_environment.sh", "environment_profile_sha256"),
        ("frontier_job.sh", "job_script_sha256"),
        (
            "reconcile_q019_registered_execution.py",
            "reconcile_q019_registered_execution_sha256",
        ),
    ):
        _, orion_payload = q043_prep._stable_regular_bytes(
            orion / name, label=f"Q019 pilot Orion installed {name}"
        )
        _, project_payload = q043_prep._stable_regular_bytes(
            project / name, label=f"Q019 pilot Project Home installed {name}"
        )
        _require(
            orion_payload == project_payload
            and _sha256_bytes(orion_payload) == final[digest_name],
            f"Q019 pilot paired installed {name} differs or drifted",
        )
    q043_record = _read_json_binding(
        str(final["q043_registered_matrix_path"]),
        str(final["q043_registered_matrix_sha256"]),
        label="Q019 pilot Q043 matrix",
    )
    q023_record = _read_json_binding(
        str(final["q023_registered_matrix_path"]),
        str(final["q023_registered_matrix_sha256"]),
        label="Q019 pilot Q023 matrix",
    )
    try:
        q043_validated = q043.validate_downstream_q023_q019_prerequisite(q043_record)
        q023_validated = q023.validate_downstream_q019_prerequisite(
            q023_record,
            q043_artifact_root=AUTHORIZED_ORION_ROOT,
            authorized_orion_root=AUTHORIZED_ORION_ROOT,
            authorized_project_home_root=CANONICAL_PROJECT_HOME_ROOT,
        )
    except (q043.AdmissionError, q023.QualificationError) as error:
        raise PreparationError("Q019 pilot predecessor matrix is not passing") from error
    _require(
        q043_validated["case_bindings_sha256"]
        == final["q043_registered_matrix_case_bindings_sha256"]
        and q043.canonical_sha256(q043_validated)
        == final["q043_registered_dependency_sha256"]
        and q023_validated["case_bindings_sha256"]
        == final["q023_registered_matrix_case_bindings_sha256"]
        and q023.canonical_sha256(q023_validated)
        == final["q023_registered_dependency_sha256"],
        "Q019 pilot predecessor matrix derived binding drifted",
    )
    return final


def materialize_q019_final_bindings(
    *,
    clean_candidate_manifest: Path,
    clean_candidate_manifest_sha256: str,
    expected_source_commit: str,
    installed_control_plane_version: str,
    q043_registered_matrix: Path,
    q023_registered_matrix: Path,
) -> dict[str, object]:
    _require(
        _SHA256.fullmatch(clean_candidate_manifest_sha256) is not None
        and _COMMIT.fullmatch(expected_source_commit) is not None
        and _SHA256.fullmatch(installed_control_plane_version) is not None,
        "Q019 pilot final-binding identity is malformed",
    )
    report = revalidate_clean_candidate(
        clean_candidate_manifest,
        expected_manifest_sha256=clean_candidate_manifest_sha256,
        expected_git_commit=expected_source_commit,
        expected_receipt_control_plane_version=installed_control_plane_version,
        control_plane_dir=(
            AUTHORIZED_ORION_ROOT / "control_plane" / installed_control_plane_version
        ),
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
        "Q019 pilot clean-candidate revalidation report drifted",
    )
    manifest_payload = _stable_binding(
        str(clean_candidate_manifest),
        clean_candidate_manifest_sha256,
        label="Q019 pilot clean-candidate manifest",
    )
    try:
        candidate = json.loads(manifest_payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PreparationError("Q019 pilot clean-candidate manifest is invalid") from error
    _require(
        type(candidate) is dict
        and type(candidate.get("source")) is dict
        and type(candidate.get("build")) is dict,
        "Q019 pilot clean-candidate manifest structure drifted",
    )
    source = candidate["source"]
    build = candidate["build"]
    candidate_dir = clean_candidate_manifest.parent
    source_archive = Path(str(source.get("archive_path", "")))
    executable = Path(str(build.get("executable_path", "")))
    _require(
        source.get("git_commit") == expected_source_commit
        and clean_candidate_manifest
        == AUTHORIZED_ORION_ROOT
        / "clean_candidates"
        / candidate_dir.name
        / "clean_candidate_manifest.json"
        and source_archive == candidate_dir / "source.tar"
        and executable == candidate_dir / "athena",
        "Q019 pilot clean-candidate lineage or path family drifted",
    )
    source_archive_sha256 = str(source.get("archive_sha256", ""))
    executable_sha256 = str(build.get("executable_sha256", ""))
    _stable_binding(
        str(source_archive), source_archive_sha256, label="Q019 pilot source archive"
    )
    _stable_binding(
        str(executable),
        executable_sha256,
        label="Q019 pilot executable",
        executable=True,
    )
    orion = AUTHORIZED_ORION_ROOT / "control_plane" / installed_control_plane_version
    project = (
        CANONICAL_PROJECT_HOME_ROOT / "control_plane" / installed_control_plane_version
    )
    controller_bindings: dict[str, tuple[str, str]] = {}
    for name in (
        "frontier_pic_environment.sh",
        "frontier_job.sh",
        "reconcile_q019_registered_execution.py",
    ):
        _, orion_payload = q043_prep._stable_regular_bytes(
            orion / name,
            label=f"Q019 pilot Orion installed {name}",
            require_executable=name.endswith(".sh"),
        )
        _, project_payload = q043_prep._stable_regular_bytes(
            project / name,
            label=f"Q019 pilot Project Home installed {name}",
            require_executable=name.endswith(".sh"),
        )
        _require(orion_payload == project_payload, f"Q019 pilot paired {name} differs")
        controller_bindings[name] = (str(orion / name), _sha256_bytes(orion_payload))
    _, q043_payload = q043_prep._stable_regular_bytes(
        q043_registered_matrix, label="Q019 pilot Q043 matrix"
    )
    _, q023_payload = q043_prep._stable_regular_bytes(
        q023_registered_matrix, label="Q019 pilot Q023 matrix"
    )
    try:
        q043_record = q043.validate_downstream_q023_q019_prerequisite(
            json.loads(q043_payload)
        )
        q023_record = q023.validate_downstream_q019_prerequisite(
            json.loads(q023_payload),
            q043_artifact_root=AUTHORIZED_ORION_ROOT,
            authorized_orion_root=AUTHORIZED_ORION_ROOT,
            authorized_project_home_root=CANONICAL_PROJECT_HOME_ROOT,
        )
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        q043.AdmissionError,
        q023.QualificationError,
    ) as error:
        raise PreparationError("Q019 pilot prerequisite matrix is invalid") from error
    environment_path, environment_sha256 = controller_bindings[
        "frontier_pic_environment.sh"
    ]
    job_path, job_sha256 = controller_bindings["frontier_job.sh"]
    reconciler_path, reconciler_sha256 = controller_bindings[
        "reconcile_q019_registered_execution.py"
    ]
    final = {
        "record_type": FINAL_BINDING_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "source_commit": expected_source_commit,
        "source_bundle_sha256": source["source_bundle_sha256"],
        "source_archive_path": str(source_archive),
        "source_archive_sha256": source_archive_sha256,
        "clean_candidate_manifest_path": str(clean_candidate_manifest),
        "clean_candidate_manifest_sha256": clean_candidate_manifest_sha256,
        "executable_path": str(executable),
        "executable_sha256": executable_sha256,
        "installed_control_plane_version": installed_control_plane_version,
        "orion_installed_control_plane_root": str(orion),
        "project_home_installed_control_plane_root": str(project),
        "environment_profile_path": environment_path,
        "environment_profile_sha256": environment_sha256,
        "job_script_path": job_path,
        "job_script_sha256": job_sha256,
        "analysis_script_paths": [reconciler_path],
        "analysis_script_sha256": [reconciler_sha256],
        "reconcile_q019_registered_execution_path": reconciler_path,
        "reconcile_q019_registered_execution_sha256": reconciler_sha256,
        "q043_registered_matrix_path": str(q043_registered_matrix),
        "q043_registered_matrix_sha256": _sha256_bytes(q043_payload),
        "q043_registered_matrix_record_type": q043_record["record_type"],
        "q043_registered_matrix_case_bindings_sha256": q043_record[
            "case_bindings_sha256"
        ],
        "q043_registered_dependency_sha256": q043.canonical_sha256(q043_record),
        "q023_registered_matrix_path": str(q023_registered_matrix),
        "q023_registered_matrix_sha256": _sha256_bytes(q023_payload),
        "q023_registered_matrix_record_type": q023_record["record_type"],
        "q023_registered_matrix_case_bindings_sha256": q023_record[
            "case_bindings_sha256"
        ],
        "q023_registered_dependency_sha256": q023.canonical_sha256(q023_record),
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
    return (
        "final_digests_bound_review_candidate_still_blocked",
        validate_final_bindings(final_bindings),
        [],
    )


def _pilot_members() -> list[dict[str, object]]:
    base_manifest = design.validate_checked_in_decks()
    runtime_manifest = controller.build_manifest()
    _require(
        json.loads(controller.CHECKED_IN_MANIFEST.read_text(encoding="utf-8"))
        == runtime_manifest,
        "Q019 runtime-controller manifest drifted",
    )
    base_by_id = {item["case_id"]: item for item in base_manifest["decks"]}
    pilots = []
    for overlay in runtime_manifest["artifacts"]:
        if overlay["authority"] != "excluded_pilot_only":
            continue
        artifact_id = str(overlay["artifact_id"])
        source_case_id = str(overlay["source_case_id"])
        _require(
            artifact_id in {
                "q019-controller-pilot-2d-baseline",
                "q019-controller-pilot-2d-instrumented",
                "q019-controller-pilot-3d-baseline",
                "q019-controller-pilot-3d-instrumented",
            }
            and source_case_id in base_by_id,
            "Q019 excluded-pilot inventory drifted",
        )
        path = controller.CHECKED_IN_ROOT / str(overlay["filename"])
        _, payload = q043_prep._stable_regular_bytes(
            path,
            label=f"Q019 excluded-pilot deck {artifact_id}",
            require_read_only=False,
        )
        _require(
            payload.decode("utf-8") == controller.render_overlay(overlay)
            and _sha256_bytes(payload) == overlay["rendered_sha256"],
            f"{artifact_id}: Q019 excluded-pilot deck drifted",
        )
        base = dict(base_by_id[source_case_id])
        pilots.append(
            {
                "artifact_id": artifact_id,
                "source_case_id": source_case_id,
                "pair_dimension": int(base["dimension"]),
                "pair_role": (
                    "instrumented" if artifact_id.endswith("-instrumented") else "baseline"
                ),
                "overlay": dict(overlay),
                "base_case": base,
                "deck_path": str(path.relative_to(REPO_ROOT)),
                "deck_sha256": overlay["rendered_sha256"],
                "deck_byte_count": len(payload),
            }
        )
    _require(
        len(pilots) == EXPECTED_ATTEMPT_COUNT
        and [item["artifact_id"] for item in pilots]
        == [
            "q019-controller-pilot-2d-baseline",
            "q019-controller-pilot-2d-instrumented",
            "q019-controller-pilot-3d-baseline",
            "q019-controller-pilot-3d-instrumented",
        ],
        "Q019 excluded-pilot ordered matrix drifted",
    )
    return pilots


def _authorization_id(index: int, member: Mapping[str, object]) -> str:
    role = str(member["pair_role"])
    dimension = int(member["pair_dimension"])
    return f"q019-pilot-{index:02d}-{dimension}d-{role}-v1"


def _resources(member: Mapping[str, object]) -> dict[str, int]:
    return dict(PILOT_RESOURCES[int(member["pair_dimension"])])


def _case_root_template() -> str:
    return str(AUTHORIZED_ORION_ROOT / RUN_NAMESPACE / SUBMISSION_ID_TEMPLATE)


def _input_deck_snapshot_template(member: Mapping[str, object]) -> str:
    return str(
        AUTHORIZED_ORION_ROOT
        / "manifests"
        / CAMPAIGN
        / SUBMISSION_ID_TEMPLATE
        / "snapshot"
        / f"{member['artifact_id']}.athinput"
    )


def _deck_binding(member: Mapping[str, object]) -> dict[str, object]:
    path = REPO_ROOT / str(member["deck_path"])
    _, payload = q043_prep._stable_regular_bytes(
        path,
        label=f"Q019 excluded-pilot deck {member['artifact_id']}",
        require_read_only=False,
    )
    _require(
        _sha256_bytes(payload) == member["deck_sha256"],
        f"{member['artifact_id']}: Q019 excluded-pilot deck digest drifted",
    )
    return {
        "path": str(member["deck_path"]),
        "sha256": str(member["deck_sha256"]),
        "byte_count": len(payload),
    }


def _output_topology(member: Mapping[str, object]) -> dict[str, object]:
    text = (REPO_ROOT / str(member["deck_path"])).read_text(encoding="utf-8")
    blocks = design.parse_athinput_text(text)
    for index, field in enumerate(REQUIRED_OUTPUT_FIELDS, 1):
        output = blocks.get(f"output{index}", {})
        _require(
            output.get("file_type") == "bin"
            and output.get("variable") == field
            and float(output.get("dt", "nan")) == 0.1
            and output.get("single_file_per_rank") == "false",
            f"{member['artifact_id']}: Q019 binary output topology drifted",
        )
    _require(
        blocks.get("output11", {}).get("file_type") == "pvtk"
        and blocks["output11"].get("variable") == "prtcl_all"
        and float(blocks["output11"].get("dt", "nan")) == 0.5
        and blocks.get("output12", {}).get("file_type") == "rst"
        and float(blocks["output12"].get("dt", "nan")) == 0.5
        and blocks["output12"].get("single_file_per_rank") == "false"
        and blocks.get("output13", {}).get("file_type") == "hst",
        f"{member['artifact_id']}: Q019 checkpoint or history topology drifted",
    )
    return {
        "shared_mpi_binary_and_restart_io": True,
        "binary_fields": list(REQUIRED_OUTPUT_FIELDS),
        "binary_cadence": 0.1,
        "particle_vtk_cadence": 0.5,
        "restart_cadence": 0.5,
        "minimum_binary_output_indices": 2,
        "minimum_checkpoint_output_indices": 2,
        "contiguous_indices_from_zero_required": True,
        "particle_vtk_indices_must_equal_restart_indices": True,
        "both_mhd_and_user_histories_required": True,
        "maximum_case_storage_bytes": _resources(member)["maximum_storage_bytes"],
    }


def _launch_contract(member: Mapping[str, object]) -> dict[str, object]:
    resources = _resources(member)
    contract = {
        "schema_version": 1,
        "executor": TRUSTED_LAUNCH_EXECUTOR,
        "pre_actions": [],
        "actions": [
            {
                "action_id": member["artifact_id"],
                "kind": "athena",
                "resources": {
                    "nodes": resources["nodes"],
                    "tasks": resources["tasks"],
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
        raise PreparationError("Q019 pilot launch contract failed validation") from error


def _expected_command(
    member: Mapping[str, object], final_bindings: Mapping[str, object]
) -> list[str]:
    resources = _resources(member)
    return [
        "srun",
        f"--nodes={resources['nodes']}",
        f"--ntasks={resources['tasks']}",
        f"--ntasks-per-node={TASKS_PER_NODE}",
        "--cpus-per-task=1",
        "--gpus-per-task=1",
        "--gpu-bind=closest",
        str(final_bindings["executable_path"]),
        "-i",
        _input_deck_snapshot_template(member),
        "-d",
        f"{_case_root_template()}/raw",
    ]


def _wrapper_evidence(member: Mapping[str, object]) -> dict[str, object]:
    resources = _resources(member)
    case_id = str(member["source_case_id"])
    return {
        "required_exact_rank_line": (
            f"Q019_REGISTERED_EXECUTION case_id={case_id} "
            f"mpi_world_size={resources['tasks']} "
            f"rank_ids={','.join(str(rank) for rank in range(resources['tasks']))}"
        ),
        "required_exact_exit_line": "Q019_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0",
        "required_exact_termination_line": "Terminating on user request",
        "required_final_evidence_status": "completed_not_acceptance_eligible",
        "required_saturation_evidence_eligible": "false",
        "required_controller_trigger_reason": 1903,
        "required_controller_trigger_cycle": 20,
        "trusted_wrapper_must_emit_each_exact_line_once": True,
    }


def validate_budget_accounting(
    value: object, launches: Sequence[Mapping[str, object]]
) -> dict[str, object]:
    _require(type(value) is dict, "Q019 pilot budget accounting must be an object")
    node_hours = sum(
        launch["resource_ceiling"]["maximum_node_hours"] for launch in launches
    )
    storage = sum(
        launch["resource_ceiling"]["maximum_storage_bytes"] for launch in launches
    )
    _require(
        len(launches) == EXPECTED_ATTEMPT_COUNT
        and value.get("pilot_attempt_count") == EXPECTED_ATTEMPT_COUNT
        and value.get("maximum_live_q019_submissions") == MAXIMUM_LIVE_Q019_SUBMISSIONS
        and value.get("maximum_attempts_per_pilot") == MAXIMUM_ATTEMPTS_PER_PILOT
        and value.get("maximum_retries_per_pilot") == MAXIMUM_RETRIES_PER_PILOT
        and value.get("maximum_registered_attempts") == EXPECTED_ATTEMPT_COUNT
        and value.get("maximum_batch_node_hours") == node_hours
        and node_hours == MAXIMUM_BATCH_NODE_HOURS
        and value.get("maximum_batch_storage_bytes") == storage
        and storage == MAXIMUM_BATCH_STORAGE_BYTES
        and value.get("live_ledger_remaining_node_hours")
        == "fresh_at_policy_promotion_not_assumed"
        and value.get("fresh_storage_preflight_required") is True
        and value.get("ledger_mutation_authorized") is False,
        "Q019 pilot budget or attempt ceiling drifted",
    )
    return dict(value)


def build_materialization(
    final_bindings: Mapping[str, object] | None = None,
) -> tuple[dict[str, object], dict[str, bytes]]:
    stage, selected, unresolved = _selected_bindings(final_bindings)
    members = _pilot_members()
    files: dict[str, bytes] = {}
    records = []
    launches = []
    slices = []
    for index, member in enumerate(members, 1):
        artifact_id = str(member["artifact_id"])
        source_case_id = str(member["source_case_id"])
        authorization_id = _authorization_id(index, member)
        deck = _deck_binding(member)
        resources = _resources(member)
        contract = _launch_contract(member)
        launch = {
            "record_type": LAUNCH_RECORD_TYPE,
            "schema_version": SCHEMA_VERSION,
            "successor_id": SUCCESSOR_ID,
            "status": "source_local_review_candidate_incomplete_launch_prohibited",
            "qualification_effect": QUALIFICATION_EFFECT,
            "campaign": CAMPAIGN,
            "binding_stage": stage,
            "artifact_id": artifact_id,
            "test_id": source_case_id,
            "authorization_id": authorization_id,
            "pair_dimension": member["pair_dimension"],
            "pair_role": member["pair_role"],
            "checked_in_overlay_deck": deck,
            "base_case_binding": {
                key: member["base_case"][key]
                for key in (
                    "case_id",
                    "path",
                    "sha256",
                    "dimension",
                    "nx",
                    "meshblock_nx",
                    "ppc",
                    "matrix_identity_fingerprint",
                )
            },
            "controller_contract": {
                "controller_identity_fingerprint": member["overlay"][
                    "controller_parameters"
                ]["controller_identity_fingerprint"],
                "monitor_dt": member["overlay"]["monitor_dt"],
                "pilot_cycle_limit": member["overlay"]["pilot_cycle_limit"],
                "expected_stop_reason": member["overlay"]["expected_stop_reason"],
                "box_edge_monitor_enabled": member["overlay"][
                    "box_edge_monitor_enabled"
                ],
                "resolution_monitor_enabled": member["overlay"][
                    "resolution_monitor_enabled"
                ],
                "diagnostic_failure_stop_armed": member["overlay"][
                    "diagnostic_failure_stop_armed"
                ],
                "saturation_evidence_eligible": False,
            },
            "selected_final_bindings": selected,
            "unresolved_final_bindings": unresolved,
            "roots": {
                "authorized_orion_case_root_template": _case_root_template(),
                "raw_output_root_template": f"{_case_root_template()}/raw",
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
                "mpi_ranks": resources["tasks"],
                "requested_nodes": resources["nodes"],
                "tasks_per_node": TASKS_PER_NODE,
                "global_nx": member["base_case"]["nx"],
                "meshblock_nx": member["base_case"]["meshblock_nx"],
                "meshblock_count": (
                    int(member["base_case"]["nx"][0])
                    // int(member["base_case"]["meshblock_nx"][0])
                )
                * (
                    int(member["base_case"]["nx"][1])
                    // int(member["base_case"]["meshblock_nx"][1])
                )
                * (
                    int(member["base_case"]["nx"][2])
                    // int(member["base_case"]["meshblock_nx"][2])
                ),
            },
            "output_topology": _output_topology(member),
            "stdout_wrapper_evidence": _wrapper_evidence(member),
            "trusted_reconciliation_interface": {
                "producer_name": "reconcile_q019_registered_execution.py",
                "producer_path": selected[
                    "reconcile_q019_registered_execution_path"
                ],
                "producer_sha256": selected[
                    "reconcile_q019_registered_execution_sha256"
                ],
                "installed_generation_only": True,
                "required_receipt_record_type": (
                    "q019_reconciled_registered_execution_receipt"
                ),
                "paired_receipts_required": True,
                "paired_receipts_must_be_byte_identical": True,
                "raw_inventory_captured_at_reconciliation": True,
            },
            "resource_ceiling": {
                "selected_qos": "normal",
                "registered_short_nonproduction": False,
                "maximum_nodes": resources["nodes"],
                "maximum_walltime_seconds": resources[
                    "scheduler_walltime_seconds"
                ],
                "athena_walltime_seconds": resources["athena_walltime_seconds"],
                "maximum_attempts": MAXIMUM_ATTEMPTS_PER_PILOT,
                "maximum_retries": MAXIMUM_RETRIES_PER_PILOT,
                "maximum_node_hours": (
                    resources["nodes"]
                    * resources["scheduler_walltime_seconds"]
                    / 3600.0
                ),
                "maximum_storage_bytes": resources["maximum_storage_bytes"],
            },
            "exclusion_boundary": {
                "registered_execution": True,
                "resource_and_controller_measurement_only": True,
                "saturation_evidence_eligible": False,
                "production_promotion_authorized": False,
                "paired_baseline_instrumented_comparison_required": True,
            },
            "required_blockers": list(REQUIRED_BLOCKERS),
            "authorization": dict(AUTHORIZATION_BOUNDARY),
        }
        policy_slice = {
            "authorization_id": authorization_id,
            "status": "review_required_not_authorized",
            "campaign": CAMPAIGN,
            "test_id": source_case_id,
            "evidence_class": EVIDENCE_CLASS,
            "physical_mode": PHYSICAL_MODE,
            "runtime_profile": RUNTIME_PROFILE,
            "selected_qos": "normal",
            "registered_short_nonproduction": False,
            "maximum_nodes": resources["nodes"],
            "maximum_walltime_seconds": resources[
                "scheduler_walltime_seconds"
            ],
            "maximum_attempts": MAXIMUM_ATTEMPTS_PER_PILOT,
            "job_script_sha256": selected["job_script_sha256"],
            "input_deck_sha256": deck["sha256"],
            "environment_profile_sha256": selected["environment_profile_sha256"],
            "analysis_script_sha256": [
                selected["reconcile_q019_registered_execution_sha256"]
            ],
            "executable_sha256": selected["executable_sha256"],
            "launch_contract_sha256": launch["launch_contract_sha256"],
            "clean_candidate_manifest_sha256": selected[
                "clean_candidate_manifest_sha256"
            ],
        }
        policy = {
            "record_type": POLICY_RECORD_TYPE,
            "schema_version": SCHEMA_VERSION,
            "successor_id": SUCCESSOR_ID,
            "status": "review_required_not_authorized_not_live_policy",
            "qualification_effect": QUALIFICATION_EFFECT,
            "campaign": CAMPAIGN,
            "binding_stage": stage,
            "artifact_id": artifact_id,
            "test_id": source_case_id,
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
        launch_path = f"launch_candidates/{artifact_id}.json"
        policy_path = f"policy_candidates/{artifact_id}.json"
        launch_payload = _json_bytes(launch)
        policy_payload = _json_bytes(policy)
        files[launch_path] = launch_payload
        files[policy_path] = policy_payload
        launches.append(launch)
        slices.append(policy_slice)
        records.append(
            {
                "attempt_index": index,
                "artifact_id": artifact_id,
                "test_id": source_case_id,
                "authorization_id": authorization_id,
                "pair_dimension": member["pair_dimension"],
                "pair_role": member["pair_role"],
                "checked_in_overlay_deck": deck,
                "launch_candidate": _binding(launch_path, launch_payload),
                "policy_candidate": _binding(policy_path, policy_payload),
            }
        )
    budget = validate_budget_accounting(
        {
            "pilot_attempt_count": EXPECTED_ATTEMPT_COUNT,
            "maximum_live_q019_submissions": MAXIMUM_LIVE_Q019_SUBMISSIONS,
            "maximum_attempts_per_pilot": MAXIMUM_ATTEMPTS_PER_PILOT,
            "maximum_retries_per_pilot": MAXIMUM_RETRIES_PER_PILOT,
            "maximum_registered_attempts": EXPECTED_ATTEMPT_COUNT,
            "maximum_batch_node_hours": MAXIMUM_BATCH_NODE_HOURS,
            "authorized_total_node_hour_cap": MAXIMUM_BATCH_NODE_HOURS,
            "minimum_live_ledger_remaining_node_hours_required": (
                MAXIMUM_BATCH_NODE_HOURS
            ),
            "live_ledger_remaining_node_hours": (
                "fresh_at_policy_promotion_not_assumed"
            ),
            "maximum_batch_storage_bytes": MAXIMUM_BATCH_STORAGE_BYTES,
            "maximum_storage_cap_bytes": MAXIMUM_BATCH_STORAGE_BYTES,
            "fresh_storage_preflight_required": True,
            "ledger_mutation_authorized": False,
        },
        launches,
    )
    fragment = {
        "record_type": POLICY_FRAGMENT_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "status": "aggregate_review_fragment_only_not_live_policy",
        "qualification_effect": QUALIFICATION_EFFECT,
        "campaign": CAMPAIGN,
        "binding_stage": stage,
        "pilot_attempt_count": EXPECTED_ATTEMPT_COUNT,
        "registered_science_slice_candidates": slices,
        "declared_batch_ceiling": budget,
        "unresolved_final_bindings": unresolved,
        "required_blockers": list(REQUIRED_BLOCKERS),
        "execution_boundary": {
            "complete_storage_policy": False,
            "live_policy_mutation_authorized": False,
            "scheduler_submission_authorized": False,
            "frontier_execution_authorized": False,
            "saturation_evidence_eligible": False,
            "scientific_claim_authorized": False,
            "publication_authorized": False,
        },
    }
    budget_payload = _json_bytes(budget)
    fragment_payload = _json_bytes(fragment)
    files["batch_budget_accounting_input.json"] = budget_payload
    files["registered_policy_slice_review_fragment.json"] = fragment_payload
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
        "pilot_attempt_count": len(records),
        "attempt_records": records,
        "policy_slices": slices,
        "batch_budget_accounting_input": _binding(
            "batch_budget_accounting_input.json", budget_payload
        ),
        "registered_policy_slice_review_fragment": _binding(
            "registered_policy_slice_review_fragment.json", fragment_payload
        ),
        "paired_pilot_contract": {
            "dimensions": [2, 3],
            "pair_roles": ["baseline", "instrumented"],
            "one_attempt_per_pair_member": True,
            "same_source_case_within_each_pair": True,
            "only_controller_instrumentation_differs_within_each_pair": True,
            "pilot_cycle_limit": 20,
            "expected_stop_reason": 1903,
            "saturation_evidence_eligible": False,
        },
        "required_blockers": list(REQUIRED_BLOCKERS),
        "execution_boundary": {
            "exact_registered_q043_matrix_required": True,
            "exact_registered_q023_matrix_required": True,
            "exact_final_clean_candidate_required": True,
            "exact_installed_control_plane_required": True,
            "empty_user_queue_required_for_every_submission": True,
            **AUTHORIZATION_BOUNDARY,
        },
    }
    _require(
        stage != "final_digests_bound_review_candidate_still_blocked"
        or not _contains_placeholder(manifest),
        "Q019 pilot final-digest materialization retained a placeholder",
    )
    return manifest, files


def validate_materialization(
    manifest: object,
    files: Mapping[str, bytes],
    final_bindings: Mapping[str, object] | None = None,
) -> tuple[dict[str, object], dict[str, bytes]]:
    expected_manifest, expected_files = build_materialization(final_bindings)
    _require(
        _strict_equal(manifest, expected_manifest),
        "Q019 pilot preparation manifest drifted",
    )
    _require(set(files) == set(expected_files), "Q019 pilot preparation files drifted")
    for path, payload in expected_files.items():
        _require(
            type(files[path]) is bytes and files[path] == payload,
            f"Q019 pilot preparation file drifted: {path}",
        )
    launches = []
    for record in expected_manifest["attempt_records"]:
        launch = json.loads(files[record["launch_candidate"]["path"]])
        policy = json.loads(files[record["policy_candidate"]["path"]])
        _require(
            validate_launch_contract(launch["launch_contract"])
            == launch["launch_contract"]
            and launch["authorization"] == AUTHORIZATION_BOUNDARY
            and policy["authorization"] == AUTHORIZATION_BOUNDARY
            and launch["exclusion_boundary"]["saturation_evidence_eligible"] is False,
            f"{record['artifact_id']}: Q019 pilot candidate acquired authority",
        )
        launches.append(launch)
    validate_budget_accounting(
        json.loads(files["batch_budget_accounting_input.json"]), launches
    )
    return expected_manifest, expected_files


def materialize_q019_timeout_margin(
    *,
    artifact_id: str,
    final_bindings: Mapping[str, object],
    measured_utc: str,
    expires_utc: str,
) -> dict[str, object]:
    final = validate_final_bindings(final_bindings)
    matches = [item for item in _pilot_members() if item["artifact_id"] == artifact_id]
    _require(len(matches) == 1, "Q019 pilot timeout artifact ID is unknown")
    resources = _resources(matches[0])
    measured = q043_prep._utc_datetime(measured_utc, label="Q019 pilot timeout measured UTC")
    expires = q043_prep._utc_datetime(expires_utc, label="Q019 pilot timeout expires UTC")
    _require(
        measured < expires and (expires - measured).total_seconds() <= 24 * 3600,
        "Q019 pilot timeout validity interval is invalid",
    )
    return {
        "artifact_id": artifact_id,
        "athena_walltime_seconds": resources["athena_walltime_seconds"],
        "scheduler_walltime_seconds": resources["scheduler_walltime_seconds"],
        "environment_profile_sha256": final["environment_profile_sha256"],
        "measured_utc": measured_utc,
        "expires_utc": expires_utc,
    }


def validate_q019_timeout_margin_artifact(
    path: Path,
    *,
    artifact_id: str,
    final_bindings: Mapping[str, object],
    now: datetime | None = None,
) -> str:
    final = validate_final_bindings(final_bindings)
    _, payload = q043_prep._stable_regular_bytes(
        path, label="Q019 pilot timeout-margin artifact"
    )
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PreparationError("Q019 pilot timeout artifact is invalid") from error
    _require(
        type(value) is dict
        and _strict_equal(
            value,
            materialize_q019_timeout_margin(
                artifact_id=artifact_id,
                final_bindings=final,
                measured_utc=value.get("measured_utc"),
                expires_utc=value.get("expires_utc"),
            ),
        ),
        "Q019 pilot timeout artifact schema or walltime drifted",
    )
    measured = q043_prep._utc_datetime(
        value["measured_utc"], label="Q019 pilot timeout measured UTC"
    )
    expires = q043_prep._utc_datetime(
        value["expires_utc"], label="Q019 pilot timeout expires UTC"
    )
    current = datetime.now(timezone.utc) if now is None else now
    _require(measured <= current < expires, "Q019 pilot timeout artifact is stale")
    return str(path)


def materialize_q019_pre_submit_config(
    *,
    artifact_id: str,
    submission_id: str,
    final_bindings: Mapping[str, object],
    pre_manifest_attestation: Path,
    timeout_margin_artifact: Path,
    queue_snapshot: Path,
    site_policy_checked_utc: str,
    now: datetime | None = None,
) -> dict[str, object]:
    final = validate_final_binding_files(final_bindings)
    matches = [
        (index, item)
        for index, item in enumerate(_pilot_members(), 1)
        if item["artifact_id"] == artifact_id
    ]
    _require(len(matches) == 1, "Q019 selected pilot artifact is unknown")
    index, member = matches[0]
    _require(
        _SUBMISSION.fullmatch(submission_id) is not None,
        "Q019 pilot submission ID is malformed",
    )
    current = datetime.now(timezone.utc) if now is None else now
    checked = q043_prep._utc_datetime(
        site_policy_checked_utc, label="Q019 pilot site-policy checked UTC"
    )
    _require(
        0 <= (current - checked).total_seconds() <= 24 * 3600,
        "Q019 pilot site-policy check is stale or in the future",
    )
    timeout = validate_q019_timeout_margin_artifact(
        timeout_margin_artifact,
        artifact_id=artifact_id,
        final_bindings=final,
        now=current,
    )
    queue = q043_prep._validate_empty_queue_snapshot(queue_snapshot)
    attestation = q043_prep._validate_pre_manifest_attestation_path(
        pre_manifest_attestation
    )
    deck = _deck_binding(member)
    return {
        "pic_root": str(AUTHORIZED_ORION_ROOT),
        "campaign": CAMPAIGN,
        "test_id": member["source_case_id"],
        "submission_scope": "registered_science",
        "registered_science_authorization_id": _authorization_id(index, member),
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
        "analysis_scripts": [final["reconcile_q019_registered_execution_path"]],
        "queue_snapshot": queue,
        "prior_case_closures": [],
        "clean_candidate_manifest": final["clean_candidate_manifest_path"],
        "launch_contract": _launch_contract(member),
    }


def materialize_q019_promotable_policy(
    *,
    baseline_policy: Mapping[str, object],
    final_bindings: Mapping[str, object],
) -> dict[str, object]:
    final = validate_final_binding_files(final_bindings)
    _require(
        isinstance(baseline_policy, Mapping)
        and _STORAGE_POLICY_REQUIRED_KEYS <= set(baseline_policy)
        and set(baseline_policy) <= _STORAGE_POLICY_ALLOWED_KEYS
        and baseline_policy.get("schema_version") == SCHEMA_VERSION
        and baseline_policy.get("registered_science_slices") == []
        and baseline_policy.get("frontier_admission_smoke")
        == {"status": "closed_after_pass"},
        "Q019 pilot baseline policy is not an empty-allowlist closed-smoke policy",
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
        "Q019 pilot baseline policy does not bind the selected controller",
    )
    try:
        validate_storage_policy(
            copy.deepcopy(dict(baseline_policy)),
            control_plane_version=str(final["installed_control_plane_version"]),
            authorized_pic_root=AUTHORIZED_ORION_ROOT,
            authorized_project_home_root=CANONICAL_PROJECT_HOME_ROOT,
        )
    except ValueError as error:
        raise PreparationError("Q019 pilot baseline storage policy is invalid") from error
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
        raise PreparationError("Q019 pilot successor policy failed validation") from error
    return successor


materialize_timeout_margin = materialize_q019_timeout_margin
materialize_pre_submit_config = materialize_q019_pre_submit_config
materialize_promotable_policy = materialize_q019_promotable_policy


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--materialize-final-bindings", action="store_true")
    parser.add_argument("--materialize-promotable-policy", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--baseline-policy", type=Path)
    parser.add_argument("--final-bindings", type=Path)
    parser.add_argument("--clean-candidate-manifest", type=Path)
    parser.add_argument("--clean-candidate-manifest-sha256")
    parser.add_argument("--expected-source-commit")
    parser.add_argument("--installed-control-plane-version")
    parser.add_argument("--q043-registered-matrix", type=Path)
    parser.add_argument("--q023-registered-matrix", type=Path)
    return parser


def main() -> None:
    args = _parser().parse_args()
    _require(
        int(args.materialize_final_bindings)
        + int(args.materialize_promotable_policy)
        <= 1,
        "select at most one Q019 pilot materialization operation",
    )
    if args.materialize_final_bindings:
        required = {
            "--output": args.output,
            "--clean-candidate-manifest": args.clean_candidate_manifest,
            "--clean-candidate-manifest-sha256": args.clean_candidate_manifest_sha256,
            "--expected-source-commit": args.expected_source_commit,
            "--installed-control-plane-version": args.installed_control_plane_version,
            "--q043-registered-matrix": args.q043_registered_matrix,
            "--q023-registered-matrix": args.q023_registered_matrix,
        }
        missing = [name for name, value in required.items() if value is None]
        _require(not missing, f"Q019 pilot final binding arguments missing: {missing}")
        result = materialize_q019_final_bindings(
            clean_candidate_manifest=args.clean_candidate_manifest,
            clean_candidate_manifest_sha256=args.clean_candidate_manifest_sha256,
            expected_source_commit=args.expected_source_commit,
            installed_control_plane_version=args.installed_control_plane_version,
            q043_registered_matrix=args.q043_registered_matrix,
            q023_registered_matrix=args.q023_registered_matrix,
        )
        _write_json_exclusive(args.output, result)
    elif args.materialize_promotable_policy:
        required = {
            "--output": args.output,
            "--baseline-policy": args.baseline_policy,
            "--final-bindings": args.final_bindings,
        }
        missing = [name for name, value in required.items() if value is None]
        _require(not missing, f"Q019 pilot policy arguments missing: {missing}")
        _, baseline_payload = q043_prep._stable_regular_bytes(
            args.baseline_policy, label="Q019 pilot baseline policy"
        )
        _, final_payload = q043_prep._stable_regular_bytes(
            args.final_bindings, label="Q019 pilot final bindings"
        )
        try:
            baseline = json.loads(baseline_payload)
            final = json.loads(final_payload)
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise PreparationError("Q019 pilot policy inputs are invalid") from error
        result = materialize_q019_promotable_policy(
            baseline_policy=baseline, final_bindings=final
        )
        _write_json_exclusive(args.output, result)
    else:
        optional = (
            args.output,
            args.baseline_policy,
            args.final_bindings,
            args.clean_candidate_manifest,
            args.clean_candidate_manifest_sha256,
            args.expected_source_commit,
            args.installed_control_plane_version,
            args.q043_registered_matrix,
            args.q023_registered_matrix,
        )
        _require(
            all(value is None for value in optional),
            "Q019 pilot binding arguments require an explicit materialization operation",
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
