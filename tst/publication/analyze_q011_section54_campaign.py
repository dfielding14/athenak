#!/usr/bin/env python3
"""Fail-closed admission gate for one Q-011 Section 5.4 campaign run.

This campaign-bound admission slice validates immutable provenance and retained
raw-artifact structure.  It intentionally does not claim that the later
numerical or qualitative Section 5.4 gates have passed.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
from typing import Any, Iterator, Mapping, Sequence
import uuid

import numpy as np

if __package__:
    from . import analyze_q011_section54_outputs as output_primitives
    from . import q011_section54_model as frozen_model
    from . import q011_section54_pressure_selection as pressure_selection
    from . import q011_section54_qualifying_campaign_execution as campaign_planner
    from .frontier_control_plane.control_plane_common import (
        AUTHORIZED_PROJECT_HOME_ROOT,
        BUILD_PROVENANCE_FILENAMES,
        _planner_expected_policy_fragment,
        project_home_ledger_root,
    )
    from .frontier_control_plane.control_plane_common import (
        read_json_bytes as read_control_plane_json_bytes,
    )
    from .frontier_control_plane.control_plane_common import (
        read_stable_regular_file_below,
    )
    from .frontier_control_plane.control_plane_common import (
        validate_clean_candidate_bundle,
    )
    from .frontier_control_plane.control_plane_common import (
        validate_planner_retention_binding,
    )
    from .frontier_control_plane import ledger as frontier_ledger
    from .immutable_orion_tree import authorized_tree_root
    from .immutable_orion_tree import loads_json_reject_duplicate_keys
    from .immutable_orion_tree import require_exact_primitive_types
    from .immutable_orion_tree import staged_verified_frozen_tree
    from .immutable_orion_tree import validate_executable_elf
    from .pvtk_particles import ParticleVTKData, read_particle_vtk
else:
    import analyze_q011_section54_outputs as output_primitives
    import q011_section54_model as frozen_model
    import q011_section54_pressure_selection as pressure_selection
    import q011_section54_qualifying_campaign_execution as campaign_planner
    from frontier_control_plane import ledger as frontier_ledger
    from frontier_control_plane.control_plane_common import (
        AUTHORIZED_PROJECT_HOME_ROOT,
        BUILD_PROVENANCE_FILENAMES,
        _planner_expected_policy_fragment,
        project_home_ledger_root,
    )
    from frontier_control_plane.control_plane_common import (
        read_json_bytes as read_control_plane_json_bytes,
    )
    from frontier_control_plane.control_plane_common import read_stable_regular_file_below
    from frontier_control_plane.control_plane_common import validate_clean_candidate_bundle
    from frontier_control_plane.control_plane_common import validate_planner_retention_binding
    from immutable_orion_tree import authorized_tree_root
    from immutable_orion_tree import loads_json_reject_duplicate_keys
    from immutable_orion_tree import require_exact_primitive_types
    from immutable_orion_tree import staged_verified_frozen_tree
    from immutable_orion_tree import validate_executable_elf
    from pvtk_particles import ParticleVTKData, read_particle_vtk


REPO_ROOT = Path(__file__).resolve().parents[2]
ANALYZER_PATH = Path(__file__).resolve()
PREREGISTRATION_PATH = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_qualifying_campaign_preregistration_successor_v2_2026-06-01.json"
)
RESTART_PREREGISTRATION_PATH = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_restart_continuation_preregistration_2026-06-01.json"
)
EXPECTED_PREREGISTRATION_SHA256 = (
    "6fd9ebbc247b6cace69f0ff61553cf198241577b457410d57d26afcbf27cdc35"
)
EXPECTED_RESTART_PREREGISTRATION_SHA256 = (
    "c3360694dc90d391c5ccf7a0620ae576733e87beea3fa974c69602c82dd866ab"
)
ORION_BULK_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
ACTIVE_DECK_SOURCE_PATH = (
    "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
)
CAMPAIGN_ANALYZER_SOURCE_PATH = (
    "tst/publication/analyze_q011_section54_campaign.py"
)
PREPARED_ARTIFACT_INVENTORY_SOURCE_PATH = (
    "tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json"
)
CAMPAIGN_ID = "Q011-SECTION54-QUALIFYING-CAMPAIGN"
MANIFEST_NAME = "campaign_manifest.json"
INVENTORY_NAME = "artifact_inventory.sha256"
FREEZE_RECEIPT_NAME = "freeze_receipt.json"
ARTIFACT_ROLE = "q011_section54_qualifying_campaign_attempt"
QUALIFICATION_SCOPE = "artifact_admission_only_no_final_claim_closure"
RESULT_RECORD_TYPE = "q011_section54_campaign_admission_result"
REGISTERED_EXECUTION_RECEIPT_NAME = "q011_section54_registered_execution_receipt.json"
_REGISTERED_EXECUTION_RECEIPT_RECORD_TYPE = (
    "q011_section54_reconciled_registered_execution_receipt"
)
_REGISTERED_EXECUTION_RECEIPT_ROLE = "immutable_reconciled_registered_execution"
_REGISTERED_EXECUTION_SCOPE = "registered_science"
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_COMMIT_PATTERN = re.compile(r"[0-9a-f]{40}")
_ATTEMPT_PATTERN = re.compile(r"[a-z0-9][a-z0-9._-]{0,127}")
_SLURM_JOB_ID_PATTERN = re.compile(r"[1-9][0-9]*")
_MODEL_OVERRIDE_PATTERN = re.compile(r"([^/=\s]+)/([^/=\s]+)=([^\s=]+)")
_RFC3339_UTC_PATTERN = re.compile(
    r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}"
    r"(?:\.[0-9]{1,6})?Z"
)
_ENDPOINT_PRODUCT_KINDS = ("rho", "bmag", "prtcl_jx", "j2", "prtcl_all")
_RESTART_PRODUCT_KINDS = (
    "restart",
    "restart_manifest",
    "restart_complete",
    "restart_manifest_complete",
)
_RUN_PRODUCT_KINDS = ("stdout",)
_SNAPSHOT_PRODUCT_KINDS = frozenset((*_ENDPOINT_PRODUCT_KINDS, *_RESTART_PRODUCT_KINDS))
_ALL_PRODUCT_KINDS = frozenset((*_ENDPOINT_PRODUCT_KINDS, *_RUN_PRODUCT_KINDS))
_ALL_PRODUCT_KINDS = frozenset((*_ALL_PRODUCT_KINDS, *_RESTART_PRODUCT_KINDS))
_BINDING_NAMES = (
    "clean_candidate_manifest",
    "executable",
    "deck",
    "analyzer",
    "preregistration",
    "campaign_plan",
    "planner_materialization_receipt",
    "attempt_contract",
    "selected_pressure_receipt",
    "analyzer_helper_source_closure_manifest",
    "registered_execution_receipt",
)
_ATTEMPT_BINDING_SHA256_FIELDS = {
    "campaign_plan": "campaign_plan_sha256",
    "planner_materialization_receipt": "planner_materialization_receipt_sha256",
    "attempt_contract": "attempt_contract_sha256",
    "selected_pressure_receipt": "selected_pressure_receipt_sha256",
    "analyzer_helper_source_closure_manifest": (
        "analyzer_helper_source_closure_manifest_sha256"
    ),
    "registered_execution_receipt": "registered_execution_receipt_sha256",
}
_SEED_OVERRIDE_NAMES = (
    "particles/pic_random_seed",
    "problem/ps_inject_seed",
    "problem/ps_seed_noise_seed",
)
_MODEL_VARIANT_ORDER = (
    "coarse_uniform_dx12",
    "three_level_amr_root_dx12_finest_dx3",
    "fine_uniform_dx3",
)
_EXPECTED_MODEL_LAUNCH_OVERRIDES = {
    "coarse_uniform_dx12": (
        "mesh_refinement/refinement=none",
        "mesh_refinement/num_levels=1",
        "problem/ps_enable_curvature_amr=false",
    ),
    "three_level_amr_root_dx12_finest_dx3": (),
    "fine_uniform_dx3": (
        "mesh/nx1=16000",
        "mesh/nx2=1040",
        "mesh_refinement/refinement=none",
        "mesh_refinement/num_levels=1",
        "problem/ps_enable_curvature_amr=false",
    ),
}
_EXPECTED_SECTION54_RUNTIME_PROJECTION = (
    ("physical_mode", "paper_mhd_pic_vl2_tsc"),
    ("state", "momentum_p_over_m"),
    ("C", "10000"),
    ("background", "coupled"),
    ("feedback", "coupled"),
    ("induction", "ideal_mhd_only"),
    ("deposition", "tsc"),
    ("deltaf", "off"),
    ("deltaf_adapt", "off"),
    ("deltaf_adapt_interval", "0"),
    ("expanding_box", "off"),
    ("expansion_law", "linear"),
    ("wave_damping", "off"),
    ("nu_in", "0"),
    ("lb_cost_per_particle", "0"),
    ("max_cell_cross", "2"),
    ("theta_max", "0.3"),
    ("restart_schema", "7"),
)
_EXPECTED_HELPER_SOURCE_PATHS = (
    "tst/publication/q011_section54_model.py",
    "tst/publication/q011_section54_pressure_pilot_execution.py",
    "tst/publication/q011_section54_pressure_selection.py",
    "tst/publication/q011_section54_restart.py",
    "tst/publication/analyze_q011_section54_outputs.py",
    "tst/publication/analyze_q011_section54_campaign.py",
    "tst/publication/analyze_q011_section54_numerical_qualification.py",
    "tst/publication/q011_section54_particles.py",
    "tst/publication/q011_section54_spatial.py",
    "tst/publication/q011_section54_artifacts.py",
    "tst/publication/publish_q011_section54_pressure_pilot_bundle.py",
    "tst/publication/analyze_q011_section54_pressure_pilot.py",
    "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
    "tst/publication/frontier_f1_structured_artifacts.py",
    "tst/publication/q011_section54_attempt_manifest_materializer.py",
    "tst/publication/publish_q011_section54_campaign_attempt.py",
    "tst/publication/immutable_orion_tree.py",
    "tst/publication/pvtk_particles.py",
    "tst/publication/q011_parallel_shock_storage_estimator.py",
    "tst/publication/frontier_control_plane/control_plane_common.py",
    "tst/publication/frontier_control_plane/ledger.py",
    "tst/publication/frontier_control_plane/operator_attestation.py",
    "tst/publication/q011_section54_qualifying_campaign_execution.py",
)
_QUALIFYING_CAMPAIGN_PLAN_RECORD_TYPE = (
    "q011_section54_qualifying_campaign_execution_plan"
)
_BASELINE_ATTEMPT_CONTRACT_RECORD_TYPE = "q011_section54_launch_prohibited_handoff_contract"
_BASELINE_ATTEMPT_COUNT = 24
_CAMPAIGN_PLAN_ARTIFACT_ROLE = (
    "q011_section54_source_local_immutable_qualifying_campaign_plan"
)
_CAMPAIGN_PLAN_QUALIFICATION_EFFECT = (
    "plan_only_no_execution_authorization_no_claim_closure"
)
_PLANNER_MATERIALIZATION_RECEIPT_RECORD_TYPE = (
    "q011_section54_qualifying_campaign_plan_materialization_receipt"
)
_PLANNER_MATERIALIZATION_RECEIPT_NAME = "materialization_receipt.json"
_PLANNER_MATERIALIZED_MEMBER_INVENTORY_ALGORITHM = (
    "sha256 of '<file_sha256>  <root-relative-path>\\n' entries ordered "
    "lexically by root-relative path"
)
_PLANNER_MATERIALIZED_MEMBER_INVENTORY_SCOPE = (
    "all materialized campaign-plan members before this receipt "
    "and recursive-freeze metadata"
)
_EXPECTED_PLANNER_FREEZE_RECEIPT = {
    "schema_version": 1,
    "artifact_role": _CAMPAIGN_PLAN_ARTIFACT_ROLE,
    "qualification_effect": _CAMPAIGN_PLAN_QUALIFICATION_EFFECT,
    "inventory_excludes": INVENTORY_NAME,
    "freeze_policy": "remove all owner, group and other write bits recursively",
}
_BASELINE_REQUIRED_SEPARATE_BOUNDARY = (
    "review_and_promote_a_registered_frontier_submission_policy_then_use_"
    "the_installed_control_plane_wrapper"
)
_PVTK_SCALARS = frozenset(
    {
        "gid",
        "ptag",
        "species",
        "cr_source",
        "macro_weight",
        "birth_time",
        "deltaf_f0",
        "deltaf_weight",
    }
)
_PVTK_INTEGER_SCALARS = frozenset({"gid", "ptag", "species", "cr_source"})
_PVTK_REAL_SCALARS = _PVTK_SCALARS - _PVTK_INTEGER_SCALARS
_EXPECTED_BIN_FIELDS = {
    "rho": ("dens",),
    "bmag": ("bmag",),
    "prtcl_jx": ("prtcl_jx",),
    "j2": ("j2",),
}
_PVTK_EXECUTION_PATTERN = re.compile(
    rb"^# vtk DataFile Version 2\.0\n"
    rb"# AthenaK particle data at time= ([^ \n]+)  "
    rb"nranks= (0|[1-9][0-9]*)  cycle=(0|[1-9][0-9]*)  variables=([^\n]+)\n"
)
_RESTART_MARKER_PATTERN = re.compile(
    rb"ATHENAK_RESTART_COMPLETE_V1\n"
    rb"size=(0|[1-9][0-9]*)\n"
    rb"fnv1a64=([0-9a-f]{16})\n"
)
_Q017_TELEMETRY_PATTERN = re.compile(
    r"^q017\.telemetry\.([A-Za-z0-9_.]+)=([^\s]+)$", re.MULTILINE
)
_Q017_REQUIRED_NAMES = frozenset(
    {
        "schema_version",
        "mpi.ranks",
        "timer.driver.seconds_rank_max",
        "timer.task_lists.seconds_rank_max",
        "timer.task_lists.seconds_rank_mean",
        "timer.task_lists.calls_rank_max",
        "timer.task_list.before_timeintegrator.seconds_rank_max",
        "timer.task_list.before_stagen.seconds_rank_max",
        "timer.task_list.stagen.seconds_rank_max",
        "timer.task_list.after_stagen.seconds_rank_max",
        "timer.task_list.after_timeintegrator.seconds_rank_max",
        "timer.output_publication.seconds_rank_max",
        "timer.output_publication.calls_rank_max",
        "timer.amr_load_balance.seconds_rank_max",
        "timer.amr_load_balance.calls_rank_max",
        "timer.particle.adaptive_deltaf.seconds_rank_max",
        "timer.particle.adaptive_deltaf.calls_rank_max",
        "timer.particle.push.seconds_rank_max",
        "timer.particle.push.calls_rank_max",
        "timer.particle.deposition.seconds_rank_max",
        "timer.particle.deposition.calls_rank_max",
        "timer.particle.migration.seconds_rank_max",
        "timer.particle.migration.calls_rank_max",
        "cycles",
        "meshblocks.total",
        "meshblocks.rank_min",
        "meshblocks.rank_max",
        "mesh.cells_per_meshblock",
        "mesh.active_cells",
        "particles.total",
        "particles.rank_min",
        "particles.rank_max",
        "updates.meshblock_cycles",
        "updates.particle_updates",
        "throughput.zone_cycles_per_second",
        "throughput.particle_updates_per_second",
        "load.meshblock_efficiency",
        "load.particle_efficiency",
        "load.cost.total",
        "load.cost.rank_min",
        "load.cost.rank_max",
        "load.cost.efficiency",
        "load.cost.invalid_meshblocks",
        "amr.enabled",
        "amr.meshblocks_created",
        "amr.meshblocks_deleted",
        "amr.meshblocks_communicated",
        "particle_memory.sync_kernel_timers",
        "particle_memory.record_bytes",
        "particle_memory.root_level",
        "particle_memory.max_level",
        "particle_memory.resident_records.bytes_total",
        "particle_memory.direct_views.allocated_snapshot_bytes_total",
        "particle_memory.direct_views.allocated_snapshot_bytes_rank_max",
        "particle_memory.athenak_owned_tracked_kokkos_views.allocated_snapshot_bytes_total",
        "particle_memory.athenak_owned_tracked_kokkos_views.allocated_snapshot_bytes_rank_max",
        "particle_memory.athenak_owned_tracked_kokkos_views.allocated_high_water_bytes_rank_sum",
        "particle_memory.athenak_owned_tracked_kokkos_views.allocated_high_water_bytes_rank_max",
        "particle_memory.paper_smooth_host_transport.allocated_snapshot_bytes_total",
        "particle_memory.paper_smooth_host_transport.allocated_snapshot_bytes_rank_max",
        "particle_memory.paper_smooth_host_transport.allocated_high_water_bytes_rank_sum",
        "particle_memory.paper_smooth_host_transport.allocated_high_water_bytes_rank_max",
        "particle_memory.invalid_records",
        "particle_memory.species.0.count",
        "particle_memory.species.0.resident_bytes",
        "particle_memory.species.1.count",
        "particle_memory.species.1.resident_bytes",
        "particle_memory.level.1.count",
        "particle_memory.level.1.resident_bytes",
        "particle_memory.level.2.count",
        "particle_memory.level.2.resident_bytes",
        "particle_memory.species.0.level.1.count",
        "particle_memory.species.0.level.2.count",
        "particle_memory.species.1.level.1.count",
        "particle_memory.species.1.level.2.count",
    }
)
_EXPECTED_FREEZE_RECEIPT = {
    "schema_version": 1,
    "artifact_role": ARTIFACT_ROLE,
    "qualification_effect": "retained_qualifying_campaign_attempt",
    "inventory_excludes": INVENTORY_NAME,
    "freeze_policy": "remove all owner, group and other write bits recursively",
}
_EXPECTED_POLICY_PROJECTION = {
    "campaign_matrix": {
        "physical_mode": "paper_mhd_pic_vl2_tsc",
        "grid_variants": [
            "coarse_uniform_dx12",
            "three_level_amr_root_dx12_finest_dx3",
            "fine_uniform_dx3",
        ],
        "qualifying_seeds": [
            23050101,
            23050102,
            23050103,
            23050104,
            23050105,
            23050106,
            23050107,
            23050108,
        ],
        "expected_baseline_attempts": 24,
        "paired_seed_rule": (
            "Use the same qualifying seed for coarse-uniform, AMR and "
            "fine-uniform variants."
        ),
    },
    "snapshot_selection": {
        "time_unit": "omega0_inverse",
        "required_times": [500.0, 1200.0],
        "absolute_match_tolerance": 1.0e-06,
        "missing_or_ambiguous_snapshot": "fail_endpoint",
    },
    "raw_bin_ids": ["rho", "bmag", "prtcl_jx", "j2"],
    "raw_pvtk_ids": ["prtcl_all"],
    "required_output_times_omega0_inverse": [
        0.0,
        100.0,
        200.0,
        300.0,
        400.0,
        500.0,
        600.0,
        700.0,
        800.0,
        900.0,
        1000.0,
        1100.0,
        1200.0,
    ],
    "raw_rst_policy": "retain_each_emitted_restart_checkpoint_and_completion_metadata",
}


class QualificationError(ValueError):
    """Raised when one campaign attempt cannot enter numerical qualification."""

    def __init__(self, message: str, *, code: str = "invalid_campaign_tree"):
        super().__init__(message)
        self.code = code


def _fail(code: str, message: str) -> None:
    raise QualificationError(message, code=code)


def _require(condition: bool, code: str, message: str) -> None:
    if not condition:
        _fail(code, message)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_path(path: Path) -> str:
    try:
        return _sha256_bytes(path.read_bytes())
    except OSError as error:
        _fail("unreadable_bound_file", f"unable to read bound file {path}: {error}")


def _object(value: object, expected: set[str], label: str) -> dict[str, Any]:
    _require(type(value) is dict, "schema_type_error", f"{label}: expected object")
    mapping = value
    _require(set(mapping) == expected, "schema_key_error", f"{label}: keys drifted")
    return mapping


def _object_with_required_keys(
    value: object, required: set[str], label: str
) -> dict[str, Any]:
    """Accept additive receipt fields while requiring the planner's stable core."""
    _require(type(value) is dict, "schema_type_error", f"{label}: expected object")
    mapping = value
    _require(
        required <= set(mapping),
        "schema_key_error",
        f"{label}: required keys are missing",
    )
    return mapping


def _list(value: object, label: str) -> list[Any]:
    _require(type(value) is list, "schema_type_error", f"{label}: expected list")
    return value


def _text(value: object, label: str) -> str:
    _require(
        type(value) is str and bool(value),
        "schema_type_error",
        f"{label}: expected text",
    )
    return value


def _exact_int(value: object, label: str) -> int:
    _require(type(value) is int, "schema_type_error", f"{label}: expected integer")
    return value


def _finite_float(value: object, label: str) -> float:
    _require(type(value) is float, "schema_type_error", f"{label}: expected float")
    _require(math.isfinite(value), "schema_type_error", f"{label}: expected finite float")
    return value


def _sha256(value: object, label: str) -> str:
    text = _text(value, label)
    _require(
        _SHA256_PATTERN.fullmatch(text) is not None,
        "schema_type_error",
        f"{label}: expected 64 lowercase hexadecimal digits",
    )
    return text


def _relative_path(value: object, label: str) -> str:
    text = _text(value, label)
    candidate = PurePosixPath(text)
    _require(
        not candidate.is_absolute()
        and text == candidate.as_posix()
        and text != "."
        and not any(character.isspace() for character in text)
        and all(part not in ("", ".", "..") for part in candidate.parts),
        "unsafe_relative_path",
        f"{label}: unsafe root-relative path {text!r}",
    )
    return text


def _git_sha1(value: object, label: str) -> str:
    text = _text(value, label)
    _require(
        _COMMIT_PATTERN.fullmatch(text) is not None,
        "schema_type_error",
        f"{label}: expected 40 lowercase hexadecimal digits",
    )
    return text


def _uuid(value: object, label: str) -> str:
    text = _text(value, label)
    try:
        parsed = uuid.UUID(text)
    except ValueError:
        _fail("schema_type_error", f"{label}: expected canonical UUID")
    _require(
        str(parsed) == text, "schema_type_error", f"{label}: expected canonical UUID"
    )
    return text


def _utc_timestamp(value: object, label: str) -> str:
    text = _text(value, label)
    _require(
        _RFC3339_UTC_PATTERN.fullmatch(text) is not None,
        "schema_type_error",
        f"{label}: expected canonical RFC-3339 UTC timestamp",
    )
    try:
        datetime.fromisoformat(text.removesuffix("Z") + "+00:00")
    except ValueError:
        _fail("schema_type_error", f"{label}: expected valid UTC timestamp")
    return text


def _frozen_candidate_path(value: object, freeze_id: str, suffix: str, label: str) -> str:
    text = _text(value, label)
    expected = (ORION_BULK_ROOT / "clean_candidates" / freeze_id / suffix).as_posix()
    _require(
        text == expected,
        "clean_candidate_schema_error",
        f"{label}: expected frozen clean-candidate path {expected!r}",
    )
    return text


def _load_json_bytes(payload: bytes, label: str) -> dict[str, Any]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        _fail("invalid_json", f"{label}: JSON is not UTF-8: {error}")
    decoded = loads_json_reject_duplicate_keys(
        text,
        error_type=QualificationError,
        label=label,
    )
    _require(type(decoded) is dict, "schema_type_error", f"{label}: expected object")
    return decoded


def _exact_match(actual: object, expected: object, label: str) -> None:
    require_exact_primitive_types(
        actual,
        expected,
        error_type=QualificationError,
        label=label,
    )
    _require(actual == expected, "policy_drift", f"{label}: value drifted")


def _policy_projection(policy: Mapping[str, Any]) -> dict[str, Any]:
    criteria = _object(
        policy.get("athenak_selected_release_criteria"),
        {
            "criteria_provenance",
            "campaign_matrix",
            "snapshot_selection",
            "ideal_injection_surface_classifier",
            "particle_filter",
            "spectrum",
            "shock_front",
            "upstream_magnetic_amplification",
            "amr_vs_fine_uniform_residuals",
            "qualitative_figure_requirements",
        },
        "preregistration/athenak_selected_release_criteria",
    )
    retention = _object(
        policy.get("artifact_retention_policy"),
        {
            "retention_scope",
            "required_output_times_omega0_inverse",
            "raw_bin_ids",
            "raw_pvtk_ids",
            "raw_rst_policy",
            "inventory_policy",
            "derived_artifact_policy",
        },
        "preregistration/artifact_retention_policy",
    )
    return {
        "campaign_matrix": criteria["campaign_matrix"],
        "snapshot_selection": criteria["snapshot_selection"],
        "raw_bin_ids": retention["raw_bin_ids"],
        "raw_pvtk_ids": retention["raw_pvtk_ids"],
        "required_output_times_omega0_inverse": retention[
            "required_output_times_omega0_inverse"
        ],
        "raw_rst_policy": retention["raw_rst_policy"],
    }


def _load_bound_policy(snapshot: Any, binding: Mapping[str, str]) -> dict[str, Any]:
    retained_path = snapshot.member_path(binding["path"])
    retained_payload = retained_path.read_bytes()
    _require(
        _sha256_bytes(retained_payload) == binding["sha256"],
        "hash_drift",
        "preregistration SHA-256 drifted",
    )
    _require(
        binding["sha256"] == EXPECTED_PREREGISTRATION_SHA256,
        "hash_drift",
        "retained preregistration does not match the frozen complete-policy SHA-256",
    )
    _require(
        _sha256_path(PREREGISTRATION_PATH) == binding["sha256"],
        "hash_drift",
        "invoked preregistration source SHA-256 drifted",
    )
    _require(
        retained_payload == PREREGISTRATION_PATH.read_bytes(),
        "hash_drift",
        "retained preregistration bytes differ from invoked frozen policy",
    )
    policy = _load_json_bytes(retained_payload, "preregistration")
    _require(
        policy.get("record_type") == "q011_section54_qualifying_campaign_preregistration"
        and type(policy.get("schema_version")) is int
        and policy["schema_version"] == 1,
        "policy_drift",
        "preregistration identity drifted",
    )
    projection = _policy_projection(policy)
    _exact_match(projection, _EXPECTED_POLICY_PROJECTION, "preregistration projection")
    return policy


def _validate_binding(value: object, label: str) -> dict[str, str]:
    binding = _object(value, {"path", "sha256"}, label)
    return {
        "path": _relative_path(binding["path"], f"{label}/path"),
        "sha256": _sha256(binding["sha256"], f"{label}/sha256"),
    }


def _validate_model_launch_overrides(value: object, label: str) -> list[str]:
    records = []
    parameters = set()
    for index, raw_record in enumerate(_list(value, label)):
        record_label = f"{label}[{index}]"
        record = _text(raw_record, record_label)
        match = _MODEL_OVERRIDE_PATTERN.fullmatch(record)
        _require(
            match is not None,
            "schema_type_error",
            f"{record_label}: malformed model override",
        )
        parameter = match.group(1, 2)
        _require(
            parameter not in parameters,
            "duplicate_attempt_binding",
            f"{record_label}: duplicate model override {parameter[0]}/{parameter[1]}",
        )
        parameters.add(parameter)
        records.append(record)
    return records


def _validate_seed_overrides(value: object, label: str) -> dict[str, int]:
    overrides = _object(value, set(_SEED_OVERRIDE_NAMES), label)
    return {
        name: _exact_int(overrides[name], f"{label}/{name}")
        for name in _SEED_OVERRIDE_NAMES
    }


def _validate_attempt_identity(
    value: object, artifact_bindings: Mapping[str, Mapping[str, str]]
) -> dict[str, Any]:
    label = "campaign manifest/attempt_identity"
    item = _object(
        value,
        {
            "campaign_plan_sha256",
            "planner_materialization_receipt_sha256",
            "attempt_contract_sha256",
            "selected_pressure_receipt_sha256",
            "model_launch_overrides",
            "seed_overrides",
            "analyzer_helper_source_closure_manifest_sha256",
            "registered_execution_receipt_sha256",
        },
        label,
    )
    parsed = {
        "campaign_plan_sha256": _sha256(
            item["campaign_plan_sha256"], f"{label}/campaign_plan_sha256"
        ),
        "planner_materialization_receipt_sha256": _sha256(
            item["planner_materialization_receipt_sha256"],
            f"{label}/planner_materialization_receipt_sha256",
        ),
        "attempt_contract_sha256": _sha256(
            item["attempt_contract_sha256"], f"{label}/attempt_contract_sha256"
        ),
        "selected_pressure_receipt_sha256": _sha256(
            item["selected_pressure_receipt_sha256"],
            f"{label}/selected_pressure_receipt_sha256",
        ),
        "model_launch_overrides": _validate_model_launch_overrides(
            item["model_launch_overrides"], f"{label}/model_launch_overrides"
        ),
        "seed_overrides": _validate_seed_overrides(
            item["seed_overrides"], f"{label}/seed_overrides"
        ),
        "analyzer_helper_source_closure_manifest_sha256": _sha256(
            item["analyzer_helper_source_closure_manifest_sha256"],
            f"{label}/analyzer_helper_source_closure_manifest_sha256",
        ),
        "registered_execution_receipt_sha256": _sha256(
            item["registered_execution_receipt_sha256"],
            f"{label}/registered_execution_receipt_sha256",
        ),
    }
    for binding_name, digest_name in _ATTEMPT_BINDING_SHA256_FIELDS.items():
        _require(
            parsed[digest_name] == artifact_bindings[binding_name]["sha256"],
            "attempt_binding_drift",
            f"{label}/{digest_name}: differs from retained {binding_name} binding",
        )
    return parsed


def _validate_manifest_schema(manifest: object, root: Path) -> dict[str, Any]:
    item = _object(
        manifest,
        {
            "schema_version",
            "record_type",
            "campaign_id",
            "qualification_scope",
            "authorized_orion_campaign_root",
            "run_identity",
            "candidate_binding",
            "artifact_bindings",
            "attempt_identity",
            "products",
        },
        "campaign manifest",
    )
    _require(
        _exact_int(item["schema_version"], "campaign manifest/schema_version") == 1,
        "schema_value_error",
        "campaign manifest/schema_version: expected 1",
    )
    _require(
        _text(item["record_type"], "campaign manifest/record_type")
        == "q011_section54_campaign_run_manifest",
        "schema_value_error",
        "campaign manifest/record_type drifted",
    )
    _require(
        _text(item["campaign_id"], "campaign manifest/campaign_id") == CAMPAIGN_ID,
        "schema_value_error",
        "campaign manifest/campaign_id drifted",
    )
    _require(
        _text(item["qualification_scope"], "campaign manifest/qualification_scope")
        == QUALIFICATION_SCOPE,
        "schema_value_error",
        "campaign manifest/qualification_scope drifted",
    )
    _require(
        _text(
            item["authorized_orion_campaign_root"],
            "campaign manifest/authorized_orion_campaign_root",
        )
        == str(root),
        "orion_root_binding_drift",
        "campaign manifest authorized Orion root binding drifted",
    )

    identity = _object(
        item["run_identity"],
        {"variant", "seed", "physical_mode", "attempt_id"},
        "campaign manifest/run_identity",
    )
    attempt_id = _text(
        identity["attempt_id"], "campaign manifest/run_identity/attempt_id"
    )
    _require(
        _ATTEMPT_PATTERN.fullmatch(attempt_id) is not None,
        "schema_value_error",
        "campaign manifest/run_identity/attempt_id is noncanonical",
    )
    parsed_identity = {
        "variant": _text(identity["variant"], "campaign manifest/run_identity/variant"),
        "seed": _exact_int(identity["seed"], "campaign manifest/run_identity/seed"),
        "physical_mode": _text(
            identity["physical_mode"], "campaign manifest/run_identity/physical_mode"
        ),
        "attempt_id": attempt_id,
    }

    candidate = _object(
        item["candidate_binding"],
        {"git_commit", "source_bundle_sha256", "clean_candidate_manifest"},
        "campaign manifest/candidate_binding",
    )
    commit = _text(
        candidate["git_commit"], "campaign manifest/candidate_binding/git_commit"
    )
    _require(
        _COMMIT_PATTERN.fullmatch(commit) is not None,
        "schema_type_error",
        "campaign manifest/candidate_binding/git_commit must be 40 lowercase hex digits",
    )
    parsed_candidate = {
        "git_commit": commit,
        "source_bundle_sha256": _sha256(
            candidate["source_bundle_sha256"],
            "campaign manifest/candidate_binding/source_bundle_sha256",
        ),
        "clean_candidate_manifest": _validate_binding(
            candidate["clean_candidate_manifest"],
            "campaign manifest/candidate_binding/clean_candidate_manifest",
        ),
    }

    raw_bindings = _object(
        item["artifact_bindings"],
        set(_BINDING_NAMES[1:]),
        "campaign manifest/artifact_bindings",
    )
    bindings = {
        name: _validate_binding(
            raw_bindings[name], f"campaign manifest/artifact_bindings/{name}"
        )
        for name in _BINDING_NAMES[1:]
    }
    attempt_identity = _validate_attempt_identity(item["attempt_identity"], bindings)
    products = [
        _validate_product(product, index)
        for index, product in enumerate(
            _list(item["products"], "campaign manifest/products")
        )
    ]
    _require(bool(products), "missing_product", "campaign manifest/products is empty")
    return {
        "authorized_orion_campaign_root": str(root),
        "run_identity": parsed_identity,
        "candidate_binding": parsed_candidate,
        "artifact_bindings": bindings,
        "attempt_identity": attempt_identity,
        "products": products,
    }


def _validate_product(value: object, index: int) -> dict[str, Any]:
    label = f"campaign manifest/products[{index}]"
    product = _object(value, {"kind", "path", "sha256", "snapshot_time"}, label)
    kind = _text(product["kind"], f"{label}/kind")
    _require(
        kind in _ALL_PRODUCT_KINDS,
        "unknown_product",
        f"{label}: unknown kind {kind!r}",
    )
    snapshot_time = product["snapshot_time"]
    if kind in _SNAPSHOT_PRODUCT_KINDS:
        snapshot_time = _finite_float(snapshot_time, f"{label}/snapshot_time")
    else:
        _require(
            snapshot_time is None,
            "schema_type_error",
            f"{label}/snapshot_time: run-level products require null",
        )
    return {
        "kind": kind,
        "path": _relative_path(product["path"], f"{label}/path"),
        "sha256": _sha256(product["sha256"], f"{label}/sha256"),
        "snapshot_time": snapshot_time,
    }


def _validate_prepared_records(value: object, label: str) -> list[dict[str, str]]:
    raw_records = _list(value, label)
    _require(
        1 <= len(raw_records) <= 4096,
        "clean_candidate_schema_error",
        f"{label}: expected between 1 and 4096 records",
    )
    records = [
        _validate_binding(record, f"{label}[{index}]")
        for index, record in enumerate(raw_records)
    ]
    paths = [record["path"] for record in records]
    _require(
        paths == sorted(set(paths)),
        "clean_candidate_schema_error",
        f"{label}: paths must be unique and canonically ordered",
    )
    return records


def _validate_clean_candidate_submodules(
    value: object, freeze_id: str
) -> list[dict[str, str]]:
    raw_records = _list(value, "clean candidate/source/submodules")
    records: list[dict[str, str]] = []
    for index, raw_record in enumerate(raw_records):
        label = f"clean candidate/source/submodules[{index}]"
        record = _object(
            raw_record,
            {
                "path",
                "archive_path",
                "archive_sha256",
                "commit_path",
                "commit_sha256",
                "git_commit",
                "git_tree",
                "worktree_status",
            },
            label,
        )
        _require(
            _text(record["worktree_status"], f"{label}/worktree_status") == "clean",
            "clean_candidate_schema_error",
            f"{label}/worktree_status: expected 'clean'",
        )
        records.append(
            {
                "path": _relative_path(record["path"], f"{label}/path"),
                "archive_path": _frozen_candidate_path(
                    record["archive_path"],
                    freeze_id,
                    f"submodules/{index:04d}.tar",
                    f"{label}/archive_path",
                ),
                "archive_sha256": _sha256(
                    record["archive_sha256"], f"{label}/archive_sha256"
                ),
                "commit_path": _frozen_candidate_path(
                    record["commit_path"],
                    freeze_id,
                    f"submodules/{index:04d}.commit",
                    f"{label}/commit_path",
                ),
                "commit_sha256": _sha256(
                    record["commit_sha256"], f"{label}/commit_sha256"
                ),
                "git_commit": _git_sha1(record["git_commit"], f"{label}/git_commit"),
                "git_tree": _git_sha1(record["git_tree"], f"{label}/git_tree"),
            }
        )
    paths = [record["path"] for record in records]
    _require(
        paths == sorted(set(paths)),
        "clean_candidate_schema_error",
        "clean candidate/source/submodules: paths must be unique and canonically ordered",
    )
    return records


def _validate_clean_candidate_manifest(payload: bytes) -> dict[str, Any]:
    label = "clean candidate"
    item = _object(
        _load_json_bytes(payload, label),
        {
            "schema_version",
            "freeze_id",
            "created_utc",
            "prepared_artifacts",
            "source",
            "build",
        },
        label,
    )
    _require(
        _exact_int(item["schema_version"], f"{label}/schema_version") == 4,
        "clean_candidate_schema_error",
        f"{label}/schema_version: expected 4",
    )
    freeze_id = _uuid(item["freeze_id"], f"{label}/freeze_id")
    created_utc = _utc_timestamp(item["created_utc"], f"{label}/created_utc")

    prepared = _object(
        item["prepared_artifacts"],
        {"inventory_path", "inventory_sha256", "paper_decks", "analyzers"},
        f"{label}/prepared_artifacts",
    )
    inventory_path = _relative_path(
        prepared["inventory_path"], f"{label}/prepared_artifacts/inventory_path"
    )
    _require(
        inventory_path == PREPARED_ARTIFACT_INVENTORY_SOURCE_PATH,
        "clean_candidate_schema_error",
        "clean candidate prepared-artifact inventory path drifted",
    )
    paper_decks = _validate_prepared_records(
        prepared["paper_decks"], f"{label}/prepared_artifacts/paper_decks"
    )
    analyzers = _validate_prepared_records(
        prepared["analyzers"], f"{label}/prepared_artifacts/analyzers"
    )
    prepared_paths = [record["path"] for record in (*paper_decks, *analyzers)]
    _require(
        inventory_path not in prepared_paths
        and len(prepared_paths) == len(set(prepared_paths)),
        "clean_candidate_schema_error",
        "clean candidate prepared-artifact paths and inventory must be distinct",
    )

    source = _object(
        item["source"],
        {
            "archive_path",
            "archive_sha256",
            "commit_path",
            "commit_sha256",
            "source_bundle_sha256",
            "git_commit",
            "git_tree",
            "worktree_status",
            "submodule_status",
            "submodules",
        },
        f"{label}/source",
    )
    _require(
        _text(source["worktree_status"], f"{label}/source/worktree_status") == "clean",
        "clean_candidate_schema_error",
        f"{label}/source/worktree_status: expected 'clean'",
    )
    submodule_status = _text(
        source["submodule_status"], f"{label}/source/submodule_status"
    )
    submodules = _validate_clean_candidate_submodules(source["submodules"], freeze_id)
    expected_submodule_status = "clean_pinned_archived" if submodules else "absent"
    _require(
        submodule_status == expected_submodule_status,
        "clean_candidate_schema_error",
        f"{label}/source/submodule_status: does not match retained submodules",
    )
    parsed_source = {
        "archive_path": _frozen_candidate_path(
            source["archive_path"],
            freeze_id,
            "source.tar",
            f"{label}/source/archive_path",
        ),
        "archive_sha256": _sha256(
            source["archive_sha256"], f"{label}/source/archive_sha256"
        ),
        "commit_path": _frozen_candidate_path(
            source["commit_path"],
            freeze_id,
            "source.commit",
            f"{label}/source/commit_path",
        ),
        "commit_sha256": _sha256(
            source["commit_sha256"], f"{label}/source/commit_sha256"
        ),
        "source_bundle_sha256": _sha256(
            source["source_bundle_sha256"], f"{label}/source/source_bundle_sha256"
        ),
        "git_commit": _git_sha1(source["git_commit"], f"{label}/source/git_commit"),
        "git_tree": _git_sha1(source["git_tree"], f"{label}/source/git_tree"),
        "submodules": submodules,
    }

    build = _object(
        item["build"],
        {
            "profile_id",
            "profile_path",
            "profile_sha256",
            "profile_receipt_path",
            "profile_receipt_sha256",
            "source_archive_sha256",
            "source_commit_sha256",
            "source_bundle_sha256",
            "toolchain",
            "build_invocations_sha256",
            "executable_path",
            "executable_sha256",
        },
        f"{label}/build",
    )
    profile_id = _text(build["profile_id"], f"{label}/build/profile_id")
    _require(
        profile_id == "hip-mpi-release-paper-pic",
        "clean_candidate_schema_error",
        f"{label}/build/profile_id: expected production PIC build profile",
    )
    toolchain = _text(build["toolchain"], f"{label}/build/toolchain")
    _require(
        bool(toolchain.strip()),
        "clean_candidate_schema_error",
        f"{label}/build/toolchain: expected nonblank text",
    )
    parsed_build = {
        "profile_id": profile_id,
        "profile_path": _frozen_candidate_path(
            build["profile_path"],
            freeze_id,
            "build_profile.json",
            f"{label}/build/profile_path",
        ),
        "profile_sha256": _sha256(
            build["profile_sha256"], f"{label}/build/profile_sha256"
        ),
        "profile_receipt_path": _frozen_candidate_path(
            build["profile_receipt_path"],
            freeze_id,
            "profile_receipt.json",
            f"{label}/build/profile_receipt_path",
        ),
        "profile_receipt_sha256": _sha256(
            build["profile_receipt_sha256"], f"{label}/build/profile_receipt_sha256"
        ),
        "source_archive_sha256": _sha256(
            build["source_archive_sha256"], f"{label}/build/source_archive_sha256"
        ),
        "source_commit_sha256": _sha256(
            build["source_commit_sha256"], f"{label}/build/source_commit_sha256"
        ),
        "source_bundle_sha256": _sha256(
            build["source_bundle_sha256"], f"{label}/build/source_bundle_sha256"
        ),
        "toolchain": toolchain,
        "build_invocations_sha256": _sha256(
            build["build_invocations_sha256"],
            f"{label}/build/build_invocations_sha256",
        ),
        "executable_path": _frozen_candidate_path(
            build["executable_path"],
            freeze_id,
            "athena",
            f"{label}/build/executable_path",
        ),
        "executable_sha256": _sha256(
            build["executable_sha256"], f"{label}/build/executable_sha256"
        ),
    }
    _require(
        parsed_build["source_archive_sha256"] == parsed_source["archive_sha256"]
        and parsed_build["source_commit_sha256"] == parsed_source["commit_sha256"]
        and parsed_build["source_bundle_sha256"]
        == parsed_source["source_bundle_sha256"],
        "clean_candidate_binding_drift",
        "clean candidate build is not bound to its frozen source",
    )
    return {
        "schema_version": 4,
        "freeze_id": freeze_id,
        "created_utc": created_utc,
        "prepared_artifacts": {
            "inventory_path": inventory_path,
            "inventory_sha256": _sha256(
                prepared["inventory_sha256"],
                f"{label}/prepared_artifacts/inventory_sha256",
            ),
            "paper_decks": paper_decks,
            "analyzers": analyzers,
        },
        "source": parsed_source,
        "build": parsed_build,
    }


def _validate_identity(
    identity: Mapping[str, Any],
    attempt_identity_or_policy: Mapping[str, Any],
    policy: Mapping[str, Any] | None = None,
) -> None:
    attempt_identity = None if policy is None else attempt_identity_or_policy
    if policy is None:
        policy = attempt_identity_or_policy
    matrix = policy["athenak_selected_release_criteria"]["campaign_matrix"]
    _require(
        identity["variant"] in matrix["grid_variants"],
        "invalid_run_identity",
        "run identity variant is not preregistered",
    )
    _require(
        identity["seed"] in matrix["qualifying_seeds"],
        "invalid_run_identity",
        "run identity seed is not preregistered",
    )
    _require(
        identity["physical_mode"] == matrix["physical_mode"],
        "invalid_run_identity",
        "run identity physical mode is not preregistered",
    )
    try:
        binding = frozen_model.variant_binding(identity["variant"])
    except frozen_model.ModelContractError as error:
        _fail("invalid_run_identity", f"run identity variant is not canonical: {error}")
    if attempt_identity is None:
        return
    _require(
        tuple(attempt_identity["model_launch_overrides"])
        == _EXPECTED_MODEL_LAUNCH_OVERRIDES[binding.variant],
        "invalid_attempt_identity",
        "attempt model launch overrides differ from the canonical campaign binding",
    )
    _require(
        all(
            seed == identity["seed"]
            for seed in attempt_identity["seed_overrides"].values()
        ),
        "invalid_attempt_identity",
        "attempt seed overrides must exactly match the preregistered run seed",
    )


def _member_payload(snapshot: Any, binding: Mapping[str, str], label: str) -> bytes:
    try:
        path = snapshot.member_path(binding["path"])
        payload = path.read_bytes()
    except (OSError, ValueError) as error:
        _fail("missing_bound_file", f"{label}: retained member is unavailable: {error}")
    measured = _sha256_bytes(payload)
    _require(measured == binding["sha256"], "hash_drift", f"{label}: SHA-256 drifted")
    _require(bool(payload), "empty_bound_file", f"{label}: retained member is empty")
    return payload


def _member_sha256(snapshot: Any, binding: Mapping[str, str], label: str) -> str:
    return _sha256_bytes(_member_payload(snapshot, binding, label))


def _absolute_binding(value: object, label: str) -> dict[str, str]:
    binding = _object(value, {"path", "sha256"}, label)
    path = _text(binding["path"], f"{label}/path")
    _require(Path(path).is_absolute(), "schema_type_error", f"{label}/path: expected absolute path")
    return {"path": path, "sha256": _sha256(binding["sha256"], f"{label}/sha256")}


def _validate_selected_pressure_receipt(
    payload: bytes, *, authorized_pic_root: Path = ORION_BULK_ROOT
) -> dict[str, Any]:
    try:
        return pressure_selection.validate_pressure_selection_receipt_bytes(
            payload,
            authorized_pic_root=authorized_pic_root,
        )
    except pressure_selection.PressureSelectionReceiptError as error:
        _fail(
            "selected_pressure_receipt_drift",
            f"selected pressure receipt failed immutable validation: {error}",
        )


def _validate_planner_freeze_receipt_semantics(receipt: object) -> None:
    require_exact_primitive_types(
        receipt,
        _EXPECTED_PLANNER_FREEZE_RECEIPT,
        error_type=QualificationError,
        label="planner materialization freeze receipt",
    )
    _require(
        receipt == _EXPECTED_PLANNER_FREEZE_RECEIPT,
        "planner_materialization_receipt_drift",
        "planner materialization freeze receipt role or effect semantics drifted",
    )


def _validate_planner_materialization_receipt(
    payload: bytes,
    *,
    retained_campaign_plan_payload: bytes,
    retained_campaign_plan_sha256: str,
    authorized_pic_root: Path,
) -> dict[str, Any]:
    """Reopen the immutable planner tree and bind its emitted graph root bytes."""
    label = "planner materialization receipt"
    item = _object_with_required_keys(
        _load_json_bytes(payload, label),
        {
            "plan_root",
            "plan_id",
            "campaign_plan_sha256",
            "materialization_receipt",
            "materialized_member_inventory_sha256",
            "inventory_sha256",
            "inventoried_file_count",
            "baseline_attempt_count",
            "restart_continuation_carrier_count",
            "recursively_read_only",
        },
        label,
    )
    plan_id = _sha256(item["plan_id"], f"{label}/plan_id")
    plan_root = Path(_text(item["plan_root"], f"{label}/plan_root"))
    _require(
        plan_root.is_absolute()
        and plan_root.name == f"q011-section54-qualifying-campaign-plan-{plan_id}",
        "planner_materialization_receipt_drift",
        f"{label}/plan_root: expected deterministic absolute planner output root",
    )
    campaign_plan_sha256 = _sha256(
        item["campaign_plan_sha256"], f"{label}/campaign_plan_sha256"
    )
    receipt_binding = _validate_binding(
        item["materialization_receipt"], f"{label}/materialization_receipt"
    )
    _require(
        receipt_binding["path"] == _PLANNER_MATERIALIZATION_RECEIPT_NAME,
        "planner_materialization_receipt_drift",
        f"{label}: internal planner receipt path drifted",
    )
    materialized_member_inventory_sha256 = _sha256(
        item["materialized_member_inventory_sha256"],
        f"{label}/materialized_member_inventory_sha256",
    )
    inventory_sha256 = _sha256(item["inventory_sha256"], f"{label}/inventory_sha256")
    inventoried_file_count = _exact_int(
        item["inventoried_file_count"], f"{label}/inventoried_file_count"
    )
    _require(
        item["recursively_read_only"] is True
        and _exact_int(item["baseline_attempt_count"], f"{label}/baseline_attempt_count")
        == _BASELINE_ATTEMPT_COUNT
        and _exact_int(
            item["restart_continuation_carrier_count"],
            f"{label}/restart_continuation_carrier_count",
        )
        == 1,
        "planner_materialization_receipt_drift",
        f"{label}: immutable planner result summary drifted",
    )
    _require(
        campaign_plan_sha256 == retained_campaign_plan_sha256
        and _sha256_bytes(retained_campaign_plan_payload) == campaign_plan_sha256,
        "planner_materialization_receipt_drift",
        f"{label}: retained campaign plan differs from planner result digest",
    )
    try:
        with staged_verified_frozen_tree(
            plan_root,
            inventory_sha256,
            authorized_root=authorized_pic_root,
            error_type=QualificationError,
            label="Q-011 Section 5.4 immutable planner materialization tree",
        ) as (tree_report, snapshot):
            _validate_planner_freeze_receipt_semantics(tree_report["freeze_receipt"])
            _require(
                tree_report["inventoried_file_count"] == inventoried_file_count
                and tree_report["recursively_read_only"] is True,
                "planner_materialization_receipt_drift",
                f"{label}: immutable planner-tree report drifted",
            )
            emitted_payload = snapshot.member_path("campaign_plan.json").read_bytes()
            internal_receipt_payload = snapshot.member_path(
                receipt_binding["path"]
            ).read_bytes()
            frozen_inventory_payload = snapshot.member_path(INVENTORY_NAME).read_bytes()
            tree_member_payloads = {
                relative: snapshot.member_path(relative).read_bytes()
                for relative in snapshot.relative_files()
            }
            tree_directories = snapshot.relative_directories()
    except QualificationError as error:
        if error.code == "planner_materialization_receipt_drift":
            raise
        _fail(
            "planner_materialization_receipt_drift",
            f"{label}: immutable planner tree failed validation: {error}",
        )
    except (OSError, ValueError) as error:
        _fail(
            "planner_materialization_receipt_drift",
            f"{label}: emitted campaign plan is unavailable: {error}",
        )
    _require(
        emitted_payload == retained_campaign_plan_payload
        and _sha256_bytes(emitted_payload) == campaign_plan_sha256,
        "planner_materialization_receipt_drift",
        f"{label}: retained campaign plan was not emitted by the immutable planner tree",
    )
    _require(
        _sha256_bytes(internal_receipt_payload) == receipt_binding["sha256"],
        "planner_materialization_receipt_drift",
        f"{label}: internal planner receipt SHA-256 drifted",
    )
    internal = _object_with_required_keys(
        _load_json_bytes(internal_receipt_payload, f"{label}/internal"),
        {
            "record_type",
            "schema_version",
            "plan_id",
            "campaign_plan",
            "helper_source_closure",
            "tree_inventory",
        },
        f"{label}/internal",
    )
    internal_campaign_plan = _validate_binding(
        internal["campaign_plan"], f"{label}/internal/campaign_plan"
    )
    helper_source_closure = _validate_binding(
        internal["helper_source_closure"],
        f"{label}/internal/helper_source_closure",
    )
    tree_inventory = _object_with_required_keys(
        internal["tree_inventory"],
        {"algorithm", "scope", "excludes", "sha256", "inventoried_file_count"},
        f"{label}/internal/tree_inventory",
    )
    excludes = _list(tree_inventory["excludes"], f"{label}/internal/tree_inventory/excludes")
    _require(
        internal["record_type"] == _PLANNER_MATERIALIZATION_RECEIPT_RECORD_TYPE
        and _exact_int(internal["schema_version"], f"{label}/internal/schema_version")
        == 1
        and _sha256(internal["plan_id"], f"{label}/internal/plan_id") == plan_id
        and internal_campaign_plan
        == {"path": "campaign_plan.json", "sha256": campaign_plan_sha256}
        and tree_inventory["algorithm"]
        == _PLANNER_MATERIALIZED_MEMBER_INVENTORY_ALGORITHM
        and tree_inventory["scope"] == _PLANNER_MATERIALIZED_MEMBER_INVENTORY_SCOPE
        and excludes
        == [
            _PLANNER_MATERIALIZATION_RECEIPT_NAME,
            FREEZE_RECEIPT_NAME,
            INVENTORY_NAME,
        ],
        "planner_materialization_receipt_drift",
        f"{label}: internal planner receipt identity drifted",
    )
    try:
        lines = frozen_inventory_payload.decode("utf-8").splitlines(keepends=True)
    except UnicodeDecodeError as error:
        _fail(
            "planner_materialization_receipt_drift",
            f"{label}: frozen planner inventory is not UTF-8: {error}",
        )
    materialized_lines = [
        line
        for line in lines
        if line.removesuffix("\n").partition("  ")[2] not in set(excludes)
    ]
    _require(
        all(line.endswith("\n") and line.count("  ") == 1 for line in lines)
        and _sha256_bytes("".join(materialized_lines).encode("utf-8"))
        == materialized_member_inventory_sha256
        == _sha256(tree_inventory["sha256"], f"{label}/internal/tree_inventory/sha256")
        and len(materialized_lines)
        == _exact_int(
            tree_inventory["inventoried_file_count"],
            f"{label}/internal/tree_inventory/inventoried_file_count",
        ),
        "planner_materialization_receipt_drift",
        f"{label}: internal materialized-member inventory drifted",
    )
    return {
        "plan_root": str(plan_root),
        "plan_id": plan_id,
        "campaign_plan_sha256": campaign_plan_sha256,
        "materialization_receipt": receipt_binding,
        "materialized_member_inventory_sha256": materialized_member_inventory_sha256,
        "helper_source_closure": helper_source_closure,
        "inventory_sha256": inventory_sha256,
        "inventoried_file_count": inventoried_file_count,
        "baseline_attempt_count": _BASELINE_ATTEMPT_COUNT,
        "restart_continuation_carrier_count": 1,
        "recursively_read_only": True,
        "_tree_member_payloads": tree_member_payloads,
        "_tree_directories": tree_directories,
    }


def _validate_helper_source_closure(payload: bytes) -> dict[str, Any]:
    label = "analyzer helper source closure manifest"
    item = _object(
        _load_json_bytes(payload, label),
        {"record_type", "schema_version", "plan_id", "sources"},
        label,
    )
    _require(
        item["record_type"] == "q011_section54_helper_source_closure"
        and _exact_int(item["schema_version"], f"{label}/schema_version") == 1,
        "helper_source_closure_drift",
        f"{label}: identity drifted",
    )
    plan_id = _sha256(item["plan_id"], f"{label}/plan_id")
    records = [
        _validate_binding(raw_record, f"{label}/sources[{index}]")
        for index, raw_record in enumerate(_list(item["sources"], f"{label}/sources"))
    ]
    _require(
        tuple(record["path"] for record in records) == _EXPECTED_HELPER_SOURCE_PATHS,
        "helper_source_closure_drift",
        f"{label}: expected helper snapshot membership or order drifted",
    )
    for record in records:
        _require(
            _sha256_path(REPO_ROOT / record["path"]) == record["sha256"],
            "helper_source_closure_drift",
            f"{label}: invoked helper bytes drifted for {record['path']}",
        )
    return {"plan_id": plan_id, "sources": records}


def _validate_plan_candidate(
    value: object,
    parsed: Mapping[str, Any],
) -> dict[str, Any]:
    label = "campaign plan/candidate_binding"
    item = _object(
        value,
        {
            "clean_candidate_manifest",
            "freeze_id",
            "git_commit",
            "git_tree",
            "source_archive_sha256",
            "source_commit_sha256",
            "source_bundle_sha256",
            "prepared_artifact_inventory_sha256",
            "validated_submodules",
            "build_profile",
            "build_profile_receipt",
            "build_invocations_sha256",
            "executable",
            "environment_profile",
        },
        label,
    )
    clean_manifest = _absolute_binding(
        item["clean_candidate_manifest"], f"{label}/clean_candidate_manifest"
    )
    executable = _absolute_binding(item["executable"], f"{label}/executable")
    environment_item = _object(
        item["environment_profile"],
        {"path", "sha256", "control_plane_version", "reviewed_source"},
        f"{label}/environment_profile",
    )
    environment_path = _text(
        environment_item["path"], f"{label}/environment_profile/path"
    )
    _require(
        Path(environment_path).is_absolute(),
        "unsafe_absolute_path",
        f"{label}/environment_profile/path: expected absolute path",
    )
    environment = {
        "path": environment_path,
        "sha256": _sha256(
            environment_item["sha256"], f"{label}/environment_profile/sha256"
        ),
        "control_plane_version": _sha256(
            environment_item["control_plane_version"],
            f"{label}/environment_profile/control_plane_version",
        ),
        "reviewed_source": _validate_binding(
            environment_item["reviewed_source"],
            f"{label}/environment_profile/reviewed_source",
        ),
    }
    _require(
        environment["reviewed_source"]["path"]
        == "tst/publication/frontier_control_plane/frontier_pic_environment.sh"
        and _sha256_path(REPO_ROOT / environment["reviewed_source"]["path"])
        == environment["reviewed_source"]["sha256"],
        "campaign_plan_crosslink_drift",
        f"{label}: reviewed environment-profile source drifted",
    )
    _absolute_binding(item["build_profile"], f"{label}/build_profile")
    _absolute_binding(item["build_profile_receipt"], f"{label}/build_profile_receipt")
    _uuid(item["freeze_id"], f"{label}/freeze_id")
    _git_sha1(item["git_commit"], f"{label}/git_commit")
    _git_sha1(item["git_tree"], f"{label}/git_tree")
    for name in (
        "source_archive_sha256",
        "source_commit_sha256",
        "source_bundle_sha256",
        "prepared_artifact_inventory_sha256",
        "build_invocations_sha256",
    ):
        _sha256(item[name], f"{label}/{name}")
    validated_submodules = _list(
        item["validated_submodules"], f"{label}/validated_submodules"
    )
    for index, raw_record in enumerate(validated_submodules):
        submodule_label = f"{label}/validated_submodules[{index}]"
        record = _object(
            raw_record,
            {"path", "archive_sha256", "commit_sha256", "git_commit", "git_tree"},
            submodule_label,
        )
        _relative_path(record["path"], f"{submodule_label}/path")
        _sha256(record["archive_sha256"], f"{submodule_label}/archive_sha256")
        _sha256(record["commit_sha256"], f"{submodule_label}/commit_sha256")
        _git_sha1(record["git_commit"], f"{submodule_label}/git_commit")
        _git_sha1(record["git_tree"], f"{submodule_label}/git_tree")
    _require(
        [record["path"] for record in validated_submodules]
        == sorted({record["path"] for record in validated_submodules}),
        "campaign_plan_crosslink_drift",
        f"{label}: validated submodule order or uniqueness drifted",
    )
    _require(
        clean_manifest["sha256"]
        == parsed["candidate_binding"]["clean_candidate_manifest"]["sha256"]
        and item["git_commit"] == parsed["candidate_binding"]["git_commit"]
        and item["source_bundle_sha256"]
        == parsed["candidate_binding"]["source_bundle_sha256"]
        and executable["sha256"] == parsed["artifact_bindings"]["executable"]["sha256"],
        "campaign_plan_crosslink_drift",
        f"{label}: candidate cross-link drifted",
    )
    return {"executable": executable, "environment_profile": environment}


def _planner_json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _planner_binding(path: str, payload: bytes) -> dict[str, str]:
    return {"path": path, "sha256": _sha256_bytes(payload)}


def _planner_exact_match(actual: object, expected: object, label: str) -> None:
    try:
        require_exact_primitive_types(
            actual,
            expected,
            error_type=QualificationError,
            label=label,
        )
    except QualificationError as error:
        _fail("planner_materialization_receipt_drift", str(error))
    _require(
        actual == expected,
        "planner_materialization_receipt_drift",
        f"{label}: retained planner member drifted",
    )


def _planner_member_payload(
    planner_receipt: Mapping[str, Any],
    binding: Mapping[str, str],
    label: str,
) -> bytes:
    members = planner_receipt["_tree_member_payloads"]
    payload = members.get(binding["path"])
    _require(
        type(payload) is bytes and _sha256_bytes(payload) == binding["sha256"],
        "planner_materialization_receipt_drift",
        f"{label}: retained planner member is missing or hash-drifted",
    )
    return payload


def _planner_json_member(
    planner_receipt: Mapping[str, Any],
    binding: Mapping[str, str],
    label: str,
) -> dict[str, Any]:
    payload = _planner_member_payload(planner_receipt, binding, label)
    try:
        return _load_json_bytes(payload, label)
    except QualificationError as error:
        _fail("planner_materialization_receipt_drift", str(error))


def _expected_planner_directories(paths: Sequence[str]) -> set[str]:
    directories: set[str] = set()
    for path in paths:
        parent = PurePosixPath(path).parent
        while parent.as_posix() != ".":
            directories.add(parent.as_posix())
            parent = parent.parent
    return directories


def _validate_materialized_planner_graph(
    item: Mapping[str, Any],
    parsed: Mapping[str, Any],
    planner_receipt: Mapping[str, Any],
    helper_closure: Mapping[str, Any],
) -> None:
    """Reconstruct and verify the exact graph emitted by the production planner."""
    plan_id = item["plan_id"]
    campaign_root = Path(item["authorized_orion_campaign_root"])
    source_bindings = item["source_bindings"]
    candidate = item["candidate_binding"]
    selected_ps_p0 = item["selected_pressure"]["selected_case"]["problem_ps_p0"]
    expected_source_paths = {
        "pressure_selection_receipt": "bindings/human_pressure_selection_receipt.json",
        "clean_candidate_manifest": "bindings/clean_candidate_manifest.json",
        "environment_profile": "bindings/environment_profile.sh",
        "qualifying_preregistration": (
            "bindings/q011_section54_qualifying_campaign_preregistration.json"
        ),
        "restart_preregistration": (
            "bindings/q011_section54_restart_continuation_preregistration.json"
        ),
        "paper_deck": "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
    }
    expected_members: dict[str, bytes] = {}
    source_payloads: dict[str, bytes] = {}
    for name, expected_path in expected_source_paths.items():
        binding = source_bindings[name]
        _require(
            binding["path"] == expected_path,
            "planner_materialization_receipt_drift",
            f"campaign plan/source_bindings/{name}: production path drifted",
        )
        payload = _planner_member_payload(
            planner_receipt, binding, f"planner source binding {name}"
        )
        expected_members[expected_path] = payload
        source_payloads[name] = payload
    _require(
        source_bindings["pressure_selection_receipt"]["sha256"]
        == parsed["artifact_bindings"]["selected_pressure_receipt"]["sha256"]
        and source_bindings["clean_candidate_manifest"]["sha256"]
        == parsed["candidate_binding"]["clean_candidate_manifest"]["sha256"]
        and source_bindings["environment_profile"]["sha256"]
        == candidate["environment_profile"]["sha256"]
        and source_bindings["qualifying_preregistration"]["sha256"]
        == parsed["artifact_bindings"]["preregistration"]["sha256"]
        and source_bindings["restart_preregistration"]["sha256"]
        == EXPECTED_RESTART_PREREGISTRATION_SHA256
        == _sha256_path(RESTART_PREREGISTRATION_PATH)
        and source_bindings["paper_deck"]["sha256"]
        == parsed["artifact_bindings"]["deck"]["sha256"],
        "campaign_plan_crosslink_drift",
        "campaign plan retained authoritative source bytes drifted",
    )
    expected_plan_id = campaign_planner._digest_value(
        {
            "record_type": campaign_planner.PLAN_RECORD_TYPE,
            "schema_version": 1,
            "pressure_selection_receipt_sha256": source_bindings[
                "pressure_selection_receipt"
            ]["sha256"],
            "selected_case": item["selected_pressure"]["selected_case"],
            "candidate_binding": candidate,
            "source_binding_sha256": {
                name: binding["sha256"] for name, binding in source_bindings.items()
            },
            "helper_source_closure": helper_closure["sources"],
            "campaign_matrix": item["campaign_matrix"],
            "authorized_orion_root": str(ORION_BULK_ROOT),
        }
    )
    _require(
        plan_id == expected_plan_id,
        "planner_materialization_receipt_drift",
        "campaign plan ID differs from the production planner input digest",
    )
    try:
        restart_policy = campaign_planner.restart.decode_preregistration(
            source_payloads["restart_preregistration"].decode("utf-8")
        )
        campaign_planner.restart.validate_preregistration(restart_policy)
    except (UnicodeDecodeError, campaign_planner.restart.RestartPolicyError) as error:
        _fail(
            "planner_materialization_receipt_drift",
            f"planner restart preregistration is invalid: {error}",
        )

    contract_bindings = []
    descriptor_bindings = []
    descriptors = []
    index = 0
    try:
        for variant_id in _EXPECTED_POLICY_PROJECTION["campaign_matrix"]["grid_variants"]:
            variant = campaign_planner._variant_binding(variant_id)
            for seed in _EXPECTED_POLICY_PROJECTION["campaign_matrix"]["qualifying_seeds"]:
                index += 1
                attempt_id = campaign_planner._attempt_id(index, variant.variant, seed)
                artifact_root = campaign_root / "baseline" / attempt_id
                contract_path = f"launch_contracts/baseline/{attempt_id}.json"
                contract = campaign_planner._baseline_launch_contract(
                    attempt_id=attempt_id,
                    variant=variant,
                    seed=seed,
                    selected_ps_p0=selected_ps_p0,
                    candidate=candidate,
                    paper_deck_binding=source_bindings["paper_deck"],
                    artifact_root=artifact_root,
                )
                contract_payload = _planner_json_bytes(contract)
                contract_binding = _planner_binding(contract_path, contract_payload)
                _planner_exact_match(
                    _planner_json_member(
                        planner_receipt, contract_binding, f"baseline contract {attempt_id}"
                    ),
                    contract,
                    f"baseline contract {attempt_id}",
                )
                expected_members[contract_path] = contract_payload
                contract_bindings.append(contract_binding)
                descriptor = campaign_planner._attempt_descriptor(
                    index=index,
                    attempt_id=attempt_id,
                    variant=variant,
                    seed=seed,
                    selected_ps_p0=selected_ps_p0,
                    candidate=candidate,
                    artifact_root=artifact_root,
                    contract_path=contract_path,
                    contract_payload=contract_payload,
                )
                descriptor_path = f"attempts/baseline/{attempt_id}.json"
                descriptor_payload = _planner_json_bytes(descriptor)
                descriptor_binding = _planner_binding(descriptor_path, descriptor_payload)
                _planner_exact_match(
                    _planner_json_member(
                        planner_receipt,
                        descriptor_binding,
                        f"baseline descriptor {attempt_id}",
                    ),
                    descriptor,
                    f"baseline descriptor {attempt_id}",
                )
                expected_members[descriptor_path] = descriptor_payload
                descriptor_bindings.append(descriptor_binding)
                descriptors.append(descriptor)
    except campaign_planner.CampaignPlanError as error:
        _fail(
            "planner_materialization_receipt_drift",
            f"production planner graph reconstruction failed: {error}",
        )
    _planner_exact_match(
        item["baseline_attempt_descriptors"],
        descriptor_bindings,
        "campaign plan/baseline_attempt_descriptors",
    )

    source_attempt = next(
        descriptor
        for descriptor in descriptors
        if descriptor["variant"] == "three_level_amr_root_dx12_finest_dx3"
        and descriptor["qualifying_seed"]
        == _EXPECTED_POLICY_PROJECTION["campaign_matrix"]["qualifying_seeds"][0]
    )
    carrier_id = campaign_planner._restart_carrier_id(
        source_attempt["qualifying_seed"]
    )
    restart_artifact_root = campaign_root / "restart_continuation" / carrier_id
    restart_contract_path = f"launch_contracts/restart_continuation/{carrier_id}.json"
    restart_contract = campaign_planner._restart_launch_contract(
        carrier_id=carrier_id,
        source_attempt=source_attempt,
        restart_preregistration=restart_policy,
        candidate=candidate,
        paper_deck_binding=source_bindings["paper_deck"],
        artifact_root=restart_artifact_root,
    )
    restart_contract_payload = _planner_json_bytes(restart_contract)
    restart_contract_binding = _planner_binding(
        restart_contract_path, restart_contract_payload
    )
    _planner_exact_match(
        _planner_json_member(
            planner_receipt, restart_contract_binding, "restart launch contract"
        ),
        restart_contract,
        "restart launch contract",
    )
    expected_members[restart_contract_path] = restart_contract_payload
    restart_carrier = {
        "record_type": campaign_planner.RESTART_CARRIER_RECORD_TYPE,
        "schema_version": 1,
        "carrier_id": carrier_id,
        "status": "planned_not_authorized",
        "source_baseline_attempt_id": source_attempt["attempt_id"],
        "variant": source_attempt["variant"],
        "qualifying_seed": source_attempt["qualifying_seed"],
        "selected_problem_ps_p0": selected_ps_p0,
        "authorized_orion_attempt_root": str(restart_artifact_root),
        "restart_preregistration": source_bindings["restart_preregistration"],
        "checkpoint_time_omega0_inverse": restart_policy["continuation_contract"][
            "checkpoint_time_omega0_inverse"
        ],
        "retained_output_schedule_after_checkpoint_omega0_inverse": restart_policy[
            "continuation_contract"
        ]["retained_output_schedule_after_checkpoint_omega0_inverse"],
        "comparison_tolerances_max_absolute_difference": restart_policy[
            "continuation_contract"
        ]["comparison_tolerances_max_absolute_difference"],
        "launch_contract": restart_contract_binding,
    }
    restart_carrier_payload = _planner_json_bytes(restart_carrier)
    restart_carrier_binding = _planner_binding(
        "restart_continuation/amr_restart_continuation_carrier.json",
        restart_carrier_payload,
    )
    _planner_exact_match(
        _planner_json_member(planner_receipt, restart_carrier_binding, "restart carrier"),
        restart_carrier,
        "restart carrier",
    )
    expected_members[restart_carrier_binding["path"]] = restart_carrier_payload
    _planner_exact_match(
        item["restart_continuation_carrier"],
        restart_carrier_binding,
        "campaign plan/restart_continuation_carrier",
    )

    helper_binding = _planner_binding(
        "helper_source_closure.json",
        _planner_json_bytes(
            {
                "record_type": "q011_section54_helper_source_closure",
                "schema_version": 1,
                "plan_id": plan_id,
                "sources": helper_closure["sources"],
            }
        ),
    )
    helper_payload = _planner_member_payload(
        planner_receipt, helper_binding, "planner helper source closure"
    )
    expected_members[helper_binding["path"]] = helper_payload
    _planner_exact_match(item["helper_source_closure"], helper_binding, "campaign plan/helper closure")

    recompute = campaign_planner._independent_recompute_plan(
        plan_id=plan_id,
        campaign_root=campaign_root,
        qualifying_preregistration_binding=source_bindings["qualifying_preregistration"],
    )
    recompute_payload = _planner_json_bytes(recompute)
    recompute_binding = _planner_binding(
        "independent_raw_artifact_recompute_plan.json", recompute_payload
    )
    _planner_exact_match(
        _planner_json_member(planner_receipt, recompute_binding, "independent recompute plan"),
        recompute,
        "independent recompute plan",
    )
    expected_members[recompute_binding["path"]] = recompute_payload
    _planner_exact_match(
        item["independent_raw_artifact_recompute_plan"],
        recompute_binding,
        "campaign plan/independent recompute plan",
    )

    fragment = _planner_expected_policy_fragment(
        plan_id=plan_id,
        pic_root=ORION_BULK_ROOT,
        campaign_root=campaign_root,
        candidate=candidate,
        pressure_receipt_binding=source_bindings["pressure_selection_receipt"],
        contract_bindings=contract_bindings,
        restart_contract_binding=restart_contract_binding,
    )
    fragment_payload = _planner_json_bytes(fragment)
    fragment_binding = _planner_binding("nonauthorizing_policy_fragment.json", fragment_payload)
    _planner_exact_match(
        _planner_json_member(planner_receipt, fragment_binding, "nonauthorizing policy fragment"),
        fragment,
        "nonauthorizing policy fragment",
    )
    expected_members[fragment_binding["path"]] = fragment_payload
    _planner_exact_match(
        item["nonauthorizing_policy_fragment"],
        fragment_binding,
        "campaign plan/nonauthorizing policy fragment",
    )

    canonical_campaign_plan = _planner_json_bytes(item)
    _require(
        _planner_member_payload(
            planner_receipt,
            {"path": "campaign_plan.json", "sha256": _sha256_bytes(canonical_campaign_plan)},
            "canonical campaign plan",
        )
        == canonical_campaign_plan,
        "planner_materialization_receipt_drift",
        "campaign plan bytes are not the canonical production-planner encoding",
    )
    expected_members["campaign_plan.json"] = canonical_campaign_plan
    expected_files = set(expected_members) | {
        _PLANNER_MATERIALIZATION_RECEIPT_NAME,
        FREEZE_RECEIPT_NAME,
        INVENTORY_NAME,
    }
    _require(
        set(planner_receipt["_tree_member_payloads"]) == expected_files
        and planner_receipt["_tree_directories"]
        == _expected_planner_directories(tuple(expected_files)),
        "planner_materialization_receipt_drift",
        "planner materialization tree contains missing or unexpected members",
    )


def _expected_attempt_id(identity: Mapping[str, Any]) -> str:
    try:
        variant_index = _MODEL_VARIANT_ORDER.index(identity["variant"])
        seed_index = _EXPECTED_POLICY_PROJECTION["campaign_matrix"][
            "qualifying_seeds"
        ].index(identity["seed"])
    except ValueError:
        _fail("campaign_plan_crosslink_drift", "run identity is absent from campaign matrix")
    index = variant_index * 8 + seed_index + 1
    return f"baseline-{index:03d}-{identity['variant']}-seed-{identity['seed']}"


def _validate_campaign_plan(
    payload: bytes,
    parsed: Mapping[str, Any],
    planner_receipt: Mapping[str, Any],
    pressure_receipt: Mapping[str, Any],
    helper_closure: Mapping[str, Any],
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    label = "campaign plan"
    item = _object(
        _load_json_bytes(payload, label),
        {
            "record_type",
            "schema_version",
            "plan_id",
            "artifact_role",
            "qualification_effect",
            "status",
            "authorized_orion_root",
            "authorized_orion_campaign_root",
            "selected_pressure",
            "candidate_binding",
            "source_bindings",
            "helper_source_closure",
            "campaign_matrix",
            "baseline_attempt_count",
            "baseline_attempt_descriptors",
            "restart_continuation_carrier",
            "independent_raw_artifact_recompute_plan",
            "nonauthorizing_policy_fragment",
            "execution_boundary",
            "preregistration_execution_boundary",
        },
        label,
    )
    _require(
        item["record_type"] == _QUALIFYING_CAMPAIGN_PLAN_RECORD_TYPE
        and _exact_int(item["schema_version"], f"{label}/schema_version") == 1,
        "campaign_plan_crosslink_drift",
        f"{label}: identity drifted",
    )
    plan_id = _sha256(item["plan_id"], f"{label}/plan_id")
    _require(
        plan_id == planner_receipt["plan_id"] == helper_closure["plan_id"],
        "campaign_plan_crosslink_drift",
        f"{label}: planner receipt or helper closure plan ID drifted",
    )
    _require(
        item["artifact_role"] == _CAMPAIGN_PLAN_ARTIFACT_ROLE
        and item["qualification_effect"] == _CAMPAIGN_PLAN_QUALIFICATION_EFFECT
        and item["status"] == "source_local_immutable_review_plan_only"
        and item["authorized_orion_root"] == str(ORION_BULK_ROOT),
        "campaign_plan_crosslink_drift",
        f"{label}: review-only role or Orion root drifted",
    )
    campaign_root = _text(
        item["authorized_orion_campaign_root"],
        f"{label}/authorized_orion_campaign_root",
    )
    expected_campaign_root = str(
        ORION_BULK_ROOT / "campaigns" / f"q011-section54-{plan_id}"
    )
    _require(
        campaign_root == expected_campaign_root,
        "campaign_plan_crosslink_drift",
        f"{label}: authorized campaign root drifted",
    )
    selected_pressure = _object(
        item["selected_pressure"],
        {"selection_method", "selected_case", "receipt"},
        f"{label}/selected_pressure",
    )
    selected_receipt = _validate_binding(
        selected_pressure["receipt"], f"{label}/selected_pressure/receipt"
    )
    _exact_match(
        selected_pressure["selected_case"],
        pressure_receipt["selected_case"],
        f"{label}/selected_pressure/selected_case",
    )
    _require(
        selected_pressure["selection_method"] == pressure_receipt["selection_method"]
        and selected_receipt["sha256"]
        == parsed["artifact_bindings"]["selected_pressure_receipt"]["sha256"],
        "campaign_plan_crosslink_drift",
        f"{label}: selected pressure cross-link drifted",
    )
    source_bindings = _object(
        item["source_bindings"],
        {
            "pressure_selection_receipt",
            "clean_candidate_manifest",
            "environment_profile",
            "qualifying_preregistration",
            "restart_preregistration",
            "paper_deck",
        },
        f"{label}/source_bindings",
    )
    source_bindings = {
        name: _validate_binding(binding, f"{label}/source_bindings/{name}")
        for name, binding in source_bindings.items()
    }
    closure_binding = _validate_binding(
        item["helper_source_closure"], f"{label}/helper_source_closure"
    )
    _require(
        source_bindings["pressure_selection_receipt"]["sha256"]
        == parsed["artifact_bindings"]["selected_pressure_receipt"]["sha256"]
        and source_bindings["clean_candidate_manifest"]["sha256"]
        == parsed["candidate_binding"]["clean_candidate_manifest"]["sha256"]
        and source_bindings["qualifying_preregistration"]["sha256"]
        == parsed["artifact_bindings"]["preregistration"]["sha256"]
        and source_bindings["environment_profile"]["sha256"]
        == item["candidate_binding"]["environment_profile"]["sha256"]
        and source_bindings["restart_preregistration"]["sha256"]
        == EXPECTED_RESTART_PREREGISTRATION_SHA256
        == _sha256_path(RESTART_PREREGISTRATION_PATH)
        and source_bindings["paper_deck"]["sha256"]
        == parsed["artifact_bindings"]["deck"]["sha256"]
        and closure_binding["sha256"]
        == parsed["artifact_bindings"]["analyzer_helper_source_closure_manifest"][
            "sha256"
        ]
        and closure_binding == planner_receipt["helper_source_closure"],
        "campaign_plan_crosslink_drift",
        f"{label}: retained source cross-link drifted",
    )
    _validate_plan_candidate(item["candidate_binding"], parsed)
    _exact_match(
        item["campaign_matrix"],
        _EXPECTED_POLICY_PROJECTION["campaign_matrix"],
        f"{label}/campaign_matrix",
    )
    _require(
        _exact_int(item["baseline_attempt_count"], f"{label}/baseline_attempt_count")
        == _BASELINE_ATTEMPT_COUNT,
        "campaign_plan_crosslink_drift",
        f"{label}: baseline attempt count drifted",
    )
    descriptors = [
        _validate_binding(record, f"{label}/baseline_attempt_descriptors[{index}]")
        for index, record in enumerate(
            _list(item["baseline_attempt_descriptors"], f"{label}/baseline_attempt_descriptors")
        )
    ]
    _require(
        len(descriptors) == _BASELINE_ATTEMPT_COUNT
        and len({record["path"] for record in descriptors}) == _BASELINE_ATTEMPT_COUNT,
        "campaign_plan_crosslink_drift",
        f"{label}: baseline attempt descriptor set drifted",
    )
    for name in (
        "restart_continuation_carrier",
        "independent_raw_artifact_recompute_plan",
        "nonauthorizing_policy_fragment",
    ):
        _validate_binding(item[name], f"{label}/{name}")
    execution = _object(
        item["execution_boundary"],
        {
            "mutates_live_policy",
            "scheduler_calls",
            "submits_jobs",
            "infers_pressure_selection",
            "launch_authorized",
            "frontier_execution_authorized",
            "claim_closure_authorized",
        },
        f"{label}/execution_boundary",
    )
    _require(
        all(value is False for value in execution.values()),
        "campaign_plan_crosslink_drift",
        f"{label}: execution boundary must remain launch-prohibited",
    )
    _exact_match(
        item["preregistration_execution_boundary"],
        policy["qualifying_execution_bindings"],
        f"{label}/preregistration_execution_boundary",
    )
    expected_attempt_id = _expected_attempt_id(parsed["run_identity"])
    _require(
        parsed["run_identity"]["attempt_id"] == expected_attempt_id,
        "campaign_plan_crosslink_drift",
        f"{label}: run attempt ID differs from deterministic baseline identity",
    )
    _validate_materialized_planner_graph(
        item,
        parsed,
        planner_receipt,
        helper_closure,
    )
    return {
        "plan_id": plan_id,
        "authorized_orion_campaign_root": campaign_root,
        "selected_problem_ps_p0": pressure_receipt["selected_case"]["problem_ps_p0"],
        "source_bindings": source_bindings,
        "candidate_binding": item["candidate_binding"],
    }


def _validate_attempt_contract(
    payload: bytes,
    parsed: Mapping[str, Any],
    plan: Mapping[str, Any],
) -> dict[str, Any]:
    label = "attempt contract"
    item = _object(
        _load_json_bytes(payload, label),
        {
            "record_type",
            "schema_version",
            "contract_role",
            "launch_authorized",
            "scheduler_submission_authorized",
            "live_policy_mutation_authorized",
            "attempt_id",
            "variant",
            "qualifying_seed",
            "selected_problem_ps_p0",
            "executable",
            "environment_profile",
            "paper_deck",
            "authorized_orion_attempt_root",
            "argv",
            "required_separate_boundary",
        },
        label,
    )
    identity = parsed["run_identity"]
    _require(
        item["record_type"] == _BASELINE_ATTEMPT_CONTRACT_RECORD_TYPE
        and _exact_int(item["schema_version"], f"{label}/schema_version") == 1
        and item["contract_role"] == "source_local_review_handoff_only"
        and item["launch_authorized"] is False
        and item["scheduler_submission_authorized"] is False
        and item["live_policy_mutation_authorized"] is False,
        "attempt_contract_crosslink_drift",
        f"{label}: launch-prohibited identity drifted",
    )
    pressure = _finite_float(
        item["selected_problem_ps_p0"], f"{label}/selected_problem_ps_p0"
    )
    executable = _absolute_binding(item["executable"], f"{label}/executable")
    environment_item = _object(
        item["environment_profile"],
        {"path", "sha256", "control_plane_version", "reviewed_source"},
        f"{label}/environment_profile",
    )
    environment_path = _text(
        environment_item["path"], f"{label}/environment_profile/path"
    )
    _require(
        Path(environment_path).is_absolute(),
        "unsafe_absolute_path",
        f"{label}/environment_profile/path: expected absolute path",
    )
    environment = {
        "path": environment_path,
        "sha256": _sha256(
            environment_item["sha256"], f"{label}/environment_profile/sha256"
        ),
        "control_plane_version": _sha256(
            environment_item["control_plane_version"],
            f"{label}/environment_profile/control_plane_version",
        ),
        "reviewed_source": _validate_binding(
            environment_item["reviewed_source"],
            f"{label}/environment_profile/reviewed_source",
        ),
    }
    paper_deck = _validate_binding(item["paper_deck"], f"{label}/paper_deck")
    attempt_root = _text(
        item["authorized_orion_attempt_root"],
        f"{label}/authorized_orion_attempt_root",
    )
    expected_attempt_root = str(
        Path(plan["authorized_orion_campaign_root"]) / "baseline" / identity["attempt_id"]
    )
    _require(
        item["attempt_id"] == identity["attempt_id"]
        and item["variant"] == identity["variant"]
        and _exact_int(item["qualifying_seed"], f"{label}/qualifying_seed")
        == identity["seed"]
        and pressure == plan["selected_problem_ps_p0"]
        and executable["sha256"] == parsed["artifact_bindings"]["executable"]["sha256"]
        and environment == plan["candidate_binding"]["environment_profile"]
        and paper_deck == plan["source_bindings"]["paper_deck"]
        and attempt_root == expected_attempt_root,
        "attempt_contract_crosslink_drift",
        f"{label}: run or plan cross-link drifted",
    )
    seed = identity["seed"]
    expected_argv = [
        "-i",
        "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
        "-d",
        f"{attempt_root}/raw",
        f"job/basename={identity['attempt_id']}",
        f"problem/ps_p0={pressure!r}",
        f"particles/pic_random_seed={seed}",
        f"problem/ps_inject_seed={seed}",
        f"problem/ps_seed_noise_seed={seed}",
        *parsed["attempt_identity"]["model_launch_overrides"],
    ]
    _exact_match(item["argv"], expected_argv, f"{label}/argv")
    _require(
        item["required_separate_boundary"] == _BASELINE_REQUIRED_SEPARATE_BOUNDARY,
        "attempt_contract_crosslink_drift",
        f"{label}: required separate launch boundary drifted",
    )
    return {
        "record_type": item["record_type"],
        "attempt_id": identity["attempt_id"],
        "launch_authorized": False,
        "argv": expected_argv,
    }


def _validate_registered_execution_receipt(
    payload: bytes,
    *,
    expected_attempt_id: str,
    expected_source_commit: str,
    expected_executable_sha256: str,
    expected_deck_sha256: str,
    expected_environment_sha256: str,
    expected_control_plane_version: str,
    expected_argv: Sequence[str],
    expected_raw_output_root: str,
    expected_artifact_dir: str | None = None,
) -> dict[str, Any]:
    label = "registered execution receipt"
    item = _object(
        _load_json_bytes(payload, label),
        {
            "record_type",
            "schema_version",
            "receipt_role",
            "registration_scope",
            "reconciled",
            "reservation_id",
            "submission_id",
            "reconciliation_event_sha256",
            "attempt_id",
            "source_commit",
            "executable_sha256",
            "deck_sha256",
            "environment_sha256",
            "control_plane_version",
            "argv",
            "slurm_job_id",
            "slurm_terminal_state",
            "raw_output_root",
            "artifact_dir",
            "planner_retention",
            "pre_submit_manifest_sha256",
        },
        label,
    )
    argv = [
        _text(value, f"{label}/argv[{index}]")
        for index, value in enumerate(_list(item["argv"], f"{label}/argv"))
    ]
    slurm_job_id = _text(item["slurm_job_id"], f"{label}/slurm_job_id")
    raw_output_root = _text(item["raw_output_root"], f"{label}/raw_output_root")
    artifact_dir = _text(item["artifact_dir"], f"{label}/artifact_dir")
    try:
        planner_retention = validate_planner_retention_binding(
            item["planner_retention"], authorized_pic_root=ORION_BULK_ROOT
        )
    except ValueError as error:
        _fail(
            "registered_execution_receipt_drift",
            f"{label}: planner-retention binding drifted: {error}",
        )
    _require(
        item["record_type"] == _REGISTERED_EXECUTION_RECEIPT_RECORD_TYPE
        and _exact_int(item["schema_version"], f"{label}/schema_version") == 1
        and item["receipt_role"] == _REGISTERED_EXECUTION_RECEIPT_ROLE
        and item["registration_scope"] == _REGISTERED_EXECUTION_SCOPE
        and item["reconciled"] is True
        and _uuid(item["reservation_id"], f"{label}/reservation_id")
        == item["reservation_id"]
        and _uuid(item["submission_id"], f"{label}/submission_id")
        == item["submission_id"]
        and _sha256(
            item["reconciliation_event_sha256"],
            f"{label}/reconciliation_event_sha256",
        )
        == item["reconciliation_event_sha256"]
        and item["attempt_id"] == expected_attempt_id
        and _git_sha1(item["source_commit"], f"{label}/source_commit")
        == expected_source_commit
        and _sha256(item["executable_sha256"], f"{label}/executable_sha256")
        == expected_executable_sha256
        and _sha256(item["deck_sha256"], f"{label}/deck_sha256")
        == expected_deck_sha256
        and _sha256(item["environment_sha256"], f"{label}/environment_sha256")
        == expected_environment_sha256
        and _sha256(item["control_plane_version"], f"{label}/control_plane_version")
        == expected_control_plane_version
        and argv == list(expected_argv)
        and _SLURM_JOB_ID_PATTERN.fullmatch(slurm_job_id) is not None
        and item["slurm_terminal_state"] == "COMPLETED"
        and raw_output_root == expected_raw_output_root
        and artifact_dir != raw_output_root
        and (
            expected_artifact_dir is None
            or artifact_dir == expected_artifact_dir
        )
        and planner_retention["attempt_id"] == expected_attempt_id
        and planner_retention["authorized_orion_raw_root"] == expected_raw_output_root
        and planner_retention["argv"] == list(expected_argv)
        and _sha256(
            item["pre_submit_manifest_sha256"],
            f"{label}/pre_submit_manifest_sha256",
        )
        == item["pre_submit_manifest_sha256"],
        "registered_execution_receipt_drift",
        f"{label}: immutable registered-execution binding drifted",
    )
    return {
        "record_type": item["record_type"],
        "receipt_role": item["receipt_role"],
        "registration_scope": item["registration_scope"],
        "reservation_id": item["reservation_id"],
        "submission_id": item["submission_id"],
        "reconciliation_event_sha256": item["reconciliation_event_sha256"],
        "attempt_id": expected_attempt_id,
        "slurm_job_id": slurm_job_id,
        "slurm_terminal_state": item["slurm_terminal_state"],
        "raw_output_root": raw_output_root,
        "artifact_dir": artifact_dir,
        "planner_retention": planner_retention,
        "pre_submit_manifest_sha256": item["pre_submit_manifest_sha256"],
        "source_commit": expected_source_commit,
        "executable_sha256": expected_executable_sha256,
        "control_plane_version": expected_control_plane_version,
    }


@contextmanager
def _validated_registered_execution_receipt_ledger_snapshot(
    receipt: Mapping[str, Any],
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> Iterator[dict[str, str]]:
    """Require one externally mirrored completed registered-science reconciliation."""
    pic_root = Path(os.path.abspath(authorized_pic_root))
    project_home_root = project_home_ledger_root(authorized_project_home_root)
    ledger_jsonl = pic_root / "ledger/node_hours.jsonl"
    receipts_jsonl = pic_root / "ledger/mirror_receipts.jsonl"
    mirror_jsonl = project_home_root / "ledger/node_hours.jsonl"
    snapshot = frontier_ledger.validated_read_only_mirrored_state_snapshot(
        ledger_jsonl,
        receipts_jsonl,
        mirror_jsonl,
        ledger_root=pic_root,
        receipts_root=pic_root,
        mirror_root=project_home_root,
    )
    try:
        records = snapshot.__enter__()
    except (OSError, ValueError) as error:
        _fail(
            "registered_execution_receipt_ledger_drift",
            f"registered execution receipt ledger validation failed: {error}",
        )
    try:
        frontier_ledger.require_explicit_genesis(records)
        matches = [
            record
            for record in records
            if record.get("event_sha256") == receipt["reconciliation_event_sha256"]
        ]
        _require(
            len(matches) == 1,
            "registered_execution_receipt_ledger_drift",
            "registered execution receipt lacks one mirrored reconciliation event",
        )
        event = matches[0]
        _require(
            event.get("event_type") == "reconciliation"
            and event.get("submission_scope") == _REGISTERED_EXECUTION_SCOPE
            and event.get("reconciled") is True
            and event.get("state") == "COMPLETED"
            and event.get("reservation_id") == receipt["reservation_id"]
            and event.get("submission_id") == receipt["submission_id"]
            and event.get("job_id") == receipt["slurm_job_id"]
            and event.get("git_commit") == receipt["source_commit"]
            and event.get("executable_sha256") == receipt["executable_sha256"]
            and event.get("control_plane_version") == receipt["control_plane_version"]
            and event.get("manifest_sha256") == receipt["pre_submit_manifest_sha256"]
            and event.get("artifact_dir") == receipt["artifact_dir"]
            and event.get("planner_retention") == receipt["planner_retention"],
            "registered_execution_receipt_ledger_drift",
            "registered execution receipt differs from mirrored reconciliation event",
        )
    except BaseException as error:
        snapshot.__exit__(type(error), error, error.__traceback__)
        raise
    binding = {
        "reservation_id": receipt["reservation_id"],
        "submission_id": receipt["submission_id"],
        "reconciliation_event_sha256": receipt["reconciliation_event_sha256"],
    }
    try:
        yield binding
    except BaseException as error:
        snapshot.__exit__(type(error), error, error.__traceback__)
        raise
    else:
        try:
            snapshot.__exit__(None, None, None)
        except (OSError, ValueError) as error:
            _fail(
                "registered_execution_receipt_ledger_drift",
                f"registered execution receipt ledger validation failed: {error}",
            )


def _validate_registered_execution_receipt_ledger_binding(
    receipt: Mapping[str, Any],
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> dict[str, str]:
    with _validated_registered_execution_receipt_ledger_snapshot(
        receipt,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    ) as binding:
        return binding


@contextmanager
def _validated_retained_attempt_semantics_snapshot(
    parsed: Mapping[str, Any],
    snapshot: Any,
    policy: Mapping[str, Any],
    *,
    authorized_pic_root: Path = ORION_BULK_ROOT,
) -> Iterator[dict[str, Any]]:
    bindings = parsed["artifact_bindings"]
    pressure_receipt = _validate_selected_pressure_receipt(
        _member_payload(
            snapshot, bindings["selected_pressure_receipt"], "selected pressure receipt"
        ),
        authorized_pic_root=authorized_pic_root,
    )
    helper_closure = _validate_helper_source_closure(
        _member_payload(
            snapshot,
            bindings["analyzer_helper_source_closure_manifest"],
            "analyzer helper source closure manifest",
        )
    )
    campaign_plan_payload = _member_payload(
        snapshot, bindings["campaign_plan"], "campaign plan"
    )
    planner_receipt = _validate_planner_materialization_receipt(
        _member_payload(
            snapshot,
            bindings["planner_materialization_receipt"],
            "planner materialization receipt",
        ),
        retained_campaign_plan_payload=campaign_plan_payload,
        retained_campaign_plan_sha256=bindings["campaign_plan"]["sha256"],
        authorized_pic_root=authorized_pic_root,
    )
    plan = _validate_campaign_plan(
        campaign_plan_payload,
        parsed,
        planner_receipt,
        pressure_receipt,
        helper_closure,
        policy,
    )
    contract = _validate_attempt_contract(
        _member_payload(snapshot, bindings["attempt_contract"], "attempt contract"),
        parsed,
        plan,
    )
    expected_retained_root = str(
        Path(plan["authorized_orion_campaign_root"])
        / "baseline"
        / parsed["run_identity"]["attempt_id"]
    )
    _require(
        parsed["authorized_orion_campaign_root"] == expected_retained_root,
        "retained_root_binding_drift",
        "retained campaign root differs from deterministic planner-authorized "
        "destination",
    )
    registered_execution_receipt = _validate_registered_execution_receipt(
        _member_payload(
            snapshot,
            bindings["registered_execution_receipt"],
            "registered execution receipt",
        ),
        expected_attempt_id=parsed["run_identity"]["attempt_id"],
        expected_source_commit=parsed["candidate_binding"]["git_commit"],
        expected_executable_sha256=bindings["executable"]["sha256"],
        expected_deck_sha256=bindings["deck"]["sha256"],
        expected_environment_sha256=plan["candidate_binding"]["environment_profile"][
            "sha256"
        ],
        expected_control_plane_version=plan["candidate_binding"][
            "environment_profile"
        ]["control_plane_version"],
        expected_argv=contract["argv"],
        expected_raw_output_root=f"{expected_retained_root}/raw",
    )
    with _validated_registered_execution_receipt_ledger_snapshot(
        registered_execution_receipt,
        authorized_pic_root=authorized_pic_root,
    ) as registered_execution_ledger_binding:
        yield {
            "campaign_plan_id": plan["plan_id"],
            "planner_materialization_receipt": {
                name: value
                for name, value in planner_receipt.items()
                if not name.startswith("_")
            },
            "selected_pressure": pressure_receipt["selected_case"],
            "attempt_contract": contract,
            "registered_execution_receipt": registered_execution_receipt,
            "registered_execution_ledger_binding": registered_execution_ledger_binding,
            "helper_source_count": len(helper_closure["sources"]),
        }


def _validate_retained_attempt_semantics(
    parsed: Mapping[str, Any],
    snapshot: Any,
    policy: Mapping[str, Any],
    *,
    authorized_pic_root: Path = ORION_BULK_ROOT,
) -> dict[str, Any]:
    """Validate one retained attempt where no later tree reads need ledger pinning."""
    with _validated_retained_attempt_semantics_snapshot(
        parsed,
        snapshot,
        policy,
        authorized_pic_root=authorized_pic_root,
    ) as semantics:
        return semantics


def _validate_prepared_binding(
    records: Sequence[Mapping[str, str]],
    required_path: str,
    binding: Mapping[str, str],
    label: str,
) -> None:
    matches = [record for record in records if record["path"] == required_path]
    _require(
        len(matches) == 1,
        "clean_candidate_binding_drift",
        f"{label}: frozen prepared-artifact path is missing or ambiguous",
    )
    _require(
        matches[0]["sha256"] == binding["sha256"],
        "clean_candidate_binding_drift",
        f"{label}: SHA-256 differs from frozen prepared artifact",
    )


def _cross_bind_clean_candidate(
    parsed: Mapping[str, Any], frozen: Mapping[str, Any]
) -> None:
    candidate = parsed["candidate_binding"]
    bindings = parsed["artifact_bindings"]
    source = frozen["source"]
    build = frozen["build"]
    prepared = frozen["prepared_artifacts"]
    _require(
        candidate["git_commit"] == source["git_commit"],
        "clean_candidate_binding_drift",
        "campaign Git commit differs from frozen clean candidate",
    )
    _require(
        candidate["source_bundle_sha256"] == source["source_bundle_sha256"]
        and candidate["source_bundle_sha256"] == build["source_bundle_sha256"],
        "clean_candidate_binding_drift",
        "campaign source bundle differs from frozen clean candidate",
    )
    _require(
        bindings["executable"]["sha256"] == build["executable_sha256"],
        "clean_candidate_binding_drift",
        "campaign executable differs from frozen clean candidate",
    )
    _validate_prepared_binding(
        prepared["paper_decks"],
        ACTIVE_DECK_SOURCE_PATH,
        bindings["deck"],
        "campaign active VL2 deck",
    )
    _validate_prepared_binding(
        prepared["analyzers"],
        CAMPAIGN_ANALYZER_SOURCE_PATH,
        bindings["analyzer"],
        "campaign admission analyzer",
    )


def _validate_elf_identity(
    path: Path, payload: bytes, expected_sha256: str, label: str
) -> dict[str, str | bool]:
    try:
        report = validate_executable_elf(
            path,
            expected_sha256,
            error_type=QualificationError,
            label=label,
        )
    except ValueError as error:
        _fail("invalid_executable_elf", f"{label}: {error}")
    _require(
        len(payload) >= 64
        and payload[:4] == b"\x7fELF"
        and payload[4] in (1, 2)
        and payload[5] in (1, 2)
        and payload[6] == 1,
        "invalid_executable_elf",
        f"{label}: expected retained ELF executable identity",
    )
    return report


def _read_external_candidate_member(path: str, label: str) -> bytes:
    try:
        return read_stable_regular_file_below(
            Path(path),
            ORION_BULK_ROOT,
            require_read_only_mode=True,
        )
    except (OSError, ValueError) as error:
        _fail(
            "external_clean_candidate_closure_error",
            f"{label}: unable to read frozen clean-candidate member: {error}",
        )


def _validate_external_clean_candidate_closure(
    frozen: Mapping[str, Any],
    retained_manifest_payload: bytes,
    executable_sha256: str,
) -> dict[str, Any]:
    """Revalidate the referenced external source/build closure from frozen bytes."""
    candidate_root = ORION_BULK_ROOT / "clean_candidates" / frozen["freeze_id"]
    external_manifest_payload = _read_external_candidate_member(
        str(candidate_root / "clean_candidate_manifest.json"),
        "external clean candidate manifest",
    )
    _require(
        external_manifest_payload == retained_manifest_payload,
        "external_clean_candidate_closure_error",
        "retained clean-candidate manifest differs from referenced external freeze",
    )
    try:
        candidate = read_control_plane_json_bytes(
            external_manifest_payload,
            label="external clean candidate manifest",
        )
        source = frozen["source"]
        build = frozen["build"]
        build_profile = _read_external_candidate_member(
            build["profile_path"], "external clean-candidate build profile"
        )
        profile = read_control_plane_json_bytes(
            build_profile,
            label="external clean-candidate build profile",
        )
        provenance_inputs = profile.get("provenance_inputs")
        _require(
            type(provenance_inputs) is dict
            and set(provenance_inputs) == set(BUILD_PROVENANCE_FILENAMES),
            "external_clean_candidate_closure_error",
            "external clean-candidate build provenance input set is incomplete",
        )
        build_provenance: dict[str, bytes] = {}
        for name in BUILD_PROVENANCE_FILENAMES:
            record = provenance_inputs[name]
            _require(
                type(record) is dict
                and set(record) == {"path", "sha256"}
                and type(record["path"]) is str,
                "external_clean_candidate_closure_error",
                f"external clean-candidate build provenance record is malformed: {name}",
            )
            build_provenance[name] = _read_external_candidate_member(
                record["path"], f"external clean-candidate build provenance {name}"
            )
        external_executable = _read_external_candidate_member(
            build["executable_path"], "external clean-candidate executable"
        )
        executable_report = _validate_elf_identity(
            Path(build["executable_path"]),
            external_executable,
            executable_sha256,
            "external clean-candidate executable",
        )
        _require(
            _sha256_bytes(external_executable) == executable_sha256,
            "external_clean_candidate_closure_error",
            "external clean-candidate executable differs from retained executable",
        )
        submodules = source["submodules"]
        validated_submodules = validate_clean_candidate_bundle(
            candidate,
            source_archive=_read_external_candidate_member(
                source["archive_path"], "external clean-candidate source archive"
            ),
            source_commit=_read_external_candidate_member(
                source["commit_path"], "external clean-candidate source commit"
            ),
            submodule_archives=[
                _read_external_candidate_member(
                    record["archive_path"],
                    f"external clean-candidate submodule archive {index}",
                )
                for index, record in enumerate(submodules)
            ],
            submodule_commits=[
                _read_external_candidate_member(
                    record["commit_path"],
                    f"external clean-candidate submodule commit {index}",
                )
                for index, record in enumerate(submodules)
            ],
            build_profile=build_profile,
            build_profile_receipt=_read_external_candidate_member(
                build["profile_receipt_path"],
                "external clean-candidate build-profile receipt",
            ),
            build_provenance=build_provenance,
            executable_sha256=executable_sha256,
            authorized_pic_root=ORION_BULK_ROOT,
        )
    except QualificationError:
        raise
    except (OSError, ValueError) as error:
        _fail(
            "external_clean_candidate_closure_error",
            f"external clean-candidate closure validation failed: {error}",
        )
    return {
        "freeze_id": frozen["freeze_id"],
        "referenced_manifest_sha256": _sha256_bytes(external_manifest_payload),
        "retained_executable_sha256": executable_sha256,
        "executable": executable_report,
        "validated_submodule_count": len(validated_submodules),
        "validation": "control_plane_validate_clean_candidate_bundle",
    }


def _validate_bound_files(
    parsed: Mapping[str, Any], snapshot: Any
) -> tuple[dict[str, Any], dict[str, Any]]:
    candidate = parsed["candidate_binding"]["clean_candidate_manifest"]
    retained_candidate_payload = _member_payload(
        snapshot, candidate, "clean candidate manifest"
    )
    frozen = _validate_clean_candidate_manifest(retained_candidate_payload)
    retained_artifacts = {}
    for name, binding in parsed["artifact_bindings"].items():
        retained_artifacts[name] = _member_payload(snapshot, binding, name)
    _cross_bind_clean_candidate(parsed, frozen)
    _validate_elf_identity(
        snapshot.member_path(parsed["artifact_bindings"]["executable"]["path"]),
        retained_artifacts["executable"],
        parsed["artifact_bindings"]["executable"]["sha256"],
        "retained campaign executable",
    )
    external_closure = _validate_external_clean_candidate_closure(
        frozen,
        retained_candidate_payload,
        parsed["artifact_bindings"]["executable"]["sha256"],
    )
    _require(
        parsed["artifact_bindings"]["analyzer"]["sha256"] == _sha256_path(ANALYZER_PATH),
        "hash_drift",
        "invoked analyzer source SHA-256 drifted",
    )
    return frozen, external_closure


def _required_retained_times(policy: Mapping[str, Any]) -> list[float]:
    return policy["artifact_retention_policy"]["required_output_times_omega0_inverse"]


def _snapshot_tolerance(policy: Mapping[str, Any]) -> float:
    snapshot_policy = policy["athenak_selected_release_criteria"]["snapshot_selection"]
    return snapshot_policy["absolute_match_tolerance"]


def _required_time_key(snapshot_time: float, policy: Mapping[str, Any]) -> str:
    matches = [
        required_time
        for required_time in _required_retained_times(policy)
        if abs(snapshot_time - required_time) <= _snapshot_tolerance(policy)
    ]
    _require(
        len(matches) == 1,
        "snapshot_cadence_drift",
        f"snapshot time {snapshot_time!r} does not identify exactly one retained cadence time",
    )
    return f"{matches[0]:.1f}"


def _select_retained_snapshot_products(
    products: Sequence[Mapping[str, Any]], policy: Mapping[str, Any]
) -> dict[str, dict[str, dict[str, Any]]]:
    retained: dict[str, dict[str, dict[str, Any]]] = {}
    raw_products = [product for product in products if product["kind"] in _ENDPOINT_PRODUCT_KINDS]
    for product in raw_products:
        _required_time_key(product["snapshot_time"], policy)
    for required_time in _required_retained_times(policy):
        snapshot: dict[str, dict[str, Any]] = {}
        for kind in _ENDPOINT_PRODUCT_KINDS:
            matches = [
                product
                for product in raw_products
                if product["kind"] == kind
                and abs(product["snapshot_time"] - required_time)
                <= _snapshot_tolerance(policy)
            ]
            _require(
                len(matches) == 1,
                "ambiguous_snapshot_cadence" if matches else "missing_snapshot_cadence",
                f"{kind} snapshot at t={required_time:.1f} has {len(matches)} matches",
            )
            snapshot[kind] = dict(matches[0])
        retained[f"{required_time:.1f}"] = snapshot
    return retained


def _select_products(
    retained: Mapping[str, Mapping[str, Mapping[str, Any]]],
    policy: Mapping[str, Any],
) -> dict[str, dict[str, dict[str, Any]]]:
    snapshot_policy = policy["athenak_selected_release_criteria"]["snapshot_selection"]
    selected: dict[str, dict[str, dict[str, Any]]] = {}
    for required_time in snapshot_policy["required_times"]:
        key = f"{required_time:.1f}"
        _require(key in retained, "missing_endpoint", f"missing endpoint snapshot t={key}")
        selected[key] = {kind: dict(product) for kind, product in retained[key].items()}
    return selected


def _select_stdout_product(products: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    matches = [product for product in products if product["kind"] == "stdout"]
    _require(
        len(matches) == 1,
        "duplicate_product" if matches else "missing_product",
        f"stdout product has {len(matches)} matches",
    )
    return dict(matches[0])


def _fnv1a64(payload: bytes) -> int:
    digest = 14695981039346656037
    for value in payload:
        digest ^= value
        digest = (digest * 1099511628211) & ((1 << 64) - 1)
    return digest


def _validate_restart_marker(marker_payload: bytes, artifact_payload: bytes, label: str) -> None:
    match = _RESTART_MARKER_PATTERN.fullmatch(marker_payload)
    _require(
        match is not None,
        "restart_marker_corruption",
        f"{label}: malformed restart completion marker",
    )
    size = int(match.group(1))
    digest = int(match.group(2), 16)
    _require(
        size == len(artifact_payload) and digest == _fnv1a64(artifact_payload),
        "restart_marker_corruption",
        f"{label}: completion marker does not bind retained artifact bytes",
    )


def _validate_restart_payload_identity(payload: bytes, label: str) -> None:
    marker = b"<par_end>\n"
    location = payload.find(marker)
    _require(
        0 < location <= 160 * 1024 and location + len(marker) < len(payload),
        "invalid_restart_payload",
        f"{label}: retained restart lacks serialized parameter header and binary body",
    )
    try:
        header = payload[:location].decode("utf-8")
    except UnicodeDecodeError as error:
        _fail("invalid_restart_payload", f"{label}: parameter header is not UTF-8: {error}")
    _require(
        header.startswith("<") and "=" in header,
        "invalid_restart_payload",
        f"{label}: serialized parameter header is malformed",
    )


def _parse_restart_manifest(payload: bytes, manifest_path: str) -> list[dict[str, Any]]:
    item = _object(
        _load_json_bytes(payload, f"restart manifest {manifest_path}"),
        {"schema", "members"},
        f"restart manifest {manifest_path}",
    )
    _require(
        item["schema"] == "ATHENAK_RESTART_MANIFEST_V1",
        "invalid_restart_manifest",
        f"{manifest_path}: unsupported restart manifest schema",
    )
    members = []
    for index, raw_member in enumerate(_list(item["members"], f"{manifest_path}/members")):
        label = f"{manifest_path}/members[{index}]"
        member = _object(raw_member, {"path", "size", "fnv1a64"}, label)
        path = _relative_path(member["path"], f"{label}/path")
        size = _exact_int(member["size"], f"{label}/size")
        digest = _text(member["fnv1a64"], f"{label}/fnv1a64")
        _require(
            size >= 0 and re.fullmatch(r"[0-9a-f]{16}", digest) is not None,
            "invalid_restart_manifest",
            f"{label}: invalid restart member digest",
        )
        members.append({"path": path, "size": size, "fnv1a64": digest})
    _require(
        bool(members) and len({member["path"] for member in members}) == len(members),
        "invalid_restart_manifest",
        f"{manifest_path}: restart members are empty or duplicated",
    )
    return members


def _validate_restart_manifest_layout(
    manifest_path: str, members: Sequence[Mapping[str, Any]]
) -> str:
    manifest = PurePosixPath(manifest_path)
    _require(
        manifest.parent == PurePosixPath("rst") and manifest.name.endswith(".rst.manifest"),
        "restart_publication_layout_drift",
        f"{manifest_path}: restart manifest must be directly below rst/",
    )
    payload_name = manifest.name.removesuffix(".manifest")
    ranked = []
    for member in members:
        path = PurePosixPath(member["path"])
        _require(
            path.name == payload_name and path.parts[0] == "rst",
            "restart_publication_layout_drift",
            f"{manifest_path}: restart member filename or rst/ root drifted",
        )
        parent = path.parent.name
        match = re.fullmatch(r"rank_([0-9]{8})", parent)
        _require(
            path.parent == PurePosixPath("rst")
            or (match is not None and path.parent.parent == PurePosixPath("rst")),
            "restart_publication_layout_drift",
            f"{manifest_path}: restart member directory layout drifted",
        )
        ranked.append(None if match is None else int(match.group(1)))
    _require(
        all(rank is None for rank in ranked) or all(rank is not None for rank in ranked),
        "restart_publication_layout_drift",
        f"{manifest_path}: restart manifest mixes shared and per-rank members",
    )
    if ranked[0] is None:
        _require(
            len(ranked) == 1,
            "restart_publication_layout_drift",
            f"{manifest_path}: shared restart manifest must contain one member",
        )
        return "shared_mpi_io"
    _require(
        sorted(ranked) == list(range(len(ranked))),
        "restart_publication_layout_drift",
        f"{manifest_path}: per-rank restart members are not contiguous",
    )
    return "per_rank_shards"


def _validate_restart_publications(
    products: Sequence[Mapping[str, Any]], snapshot: Any, policy: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    restart_products = [
        product for product in products if product["kind"] in _RESTART_PRODUCT_KINDS
    ]
    for product in restart_products:
        _required_time_key(product["snapshot_time"], policy)
    actual = {(product["kind"], product["path"]) for product in restart_products}
    expected: set[tuple[str, str]] = set()
    report: dict[str, dict[str, Any]] = {}
    for required_time in _required_retained_times(policy):
        key = f"{required_time:.1f}"
        manifest_matches = [
            product
            for product in restart_products
            if product["kind"] == "restart_manifest"
            and abs(product["snapshot_time"] - required_time) <= _snapshot_tolerance(policy)
        ]
        _require(
            len(manifest_matches) == 1,
            "ambiguous_restart_publication" if manifest_matches else "missing_restart_publication",
            f"restart manifest at t={key} has {len(manifest_matches)} matches",
        )
        manifest_product = manifest_matches[0]
        manifest_path = manifest_product["path"]
        manifest_payload = _member_payload(
            snapshot, manifest_product, f"restart manifest at t={key}"
        )
        members = _parse_restart_manifest(manifest_payload, manifest_path)
        layout = _validate_restart_manifest_layout(manifest_path, members)
        manifest_marker_path = manifest_path + ".complete"
        expected.add(("restart_manifest", manifest_path))
        expected.add(("restart_manifest_complete", manifest_marker_path))
        marker_matches = [
            product
            for product in restart_products
            if product["kind"] == "restart_manifest_complete"
            and product["path"] == manifest_marker_path
            and abs(product["snapshot_time"] - required_time) <= _snapshot_tolerance(policy)
        ]
        _require(
            len(marker_matches) == 1,
            "restart_publication_layout_drift",
            f"{manifest_path}: manifest completion marker is missing or ambiguous",
        )
        _validate_restart_marker(
            _member_payload(snapshot, marker_matches[0], f"{manifest_path} marker"),
            manifest_payload,
            f"{manifest_path} marker",
        )
        for member in members:
            member_path = member["path"]
            marker_path = member_path + ".complete"
            expected.add(("restart", member_path))
            expected.add(("restart_complete", marker_path))
            payload_matches = [
                product
                for product in restart_products
                if product["kind"] == "restart"
                and product["path"] == member_path
                and abs(product["snapshot_time"] - required_time)
                <= _snapshot_tolerance(policy)
            ]
            marker_matches = [
                product
                for product in restart_products
                if product["kind"] == "restart_complete"
                and product["path"] == marker_path
                and abs(product["snapshot_time"] - required_time)
                <= _snapshot_tolerance(policy)
            ]
            _require(
                len(payload_matches) == 1 and len(marker_matches) == 1,
                "restart_publication_layout_drift",
                f"{member_path}: payload or completion marker is missing or ambiguous",
            )
            member_payload = _member_payload(
                snapshot, payload_matches[0], f"restart payload {member_path}"
            )
            _validate_restart_payload_identity(member_payload, member_path)
            _validate_restart_marker(
                _member_payload(snapshot, marker_matches[0], f"{member_path} marker"),
                member_payload,
                f"{member_path} marker",
            )
            _require(
                member["size"] == len(member_payload)
                and int(member["fnv1a64"], 16) == _fnv1a64(member_payload),
                "invalid_restart_manifest",
                f"{manifest_path}: member digest drifted for {member_path}",
            )
        report[key] = {
            "manifest_path": manifest_path,
            "layout": layout,
            "member_count": len(members),
        }
    _require(
        actual == expected,
        "restart_publication_layout_drift",
        "restart product inventory contains missing or unexpected publication members",
    )
    return report


def _validate_declared_tree_members(parsed: Mapping[str, Any], snapshot: Any) -> None:
    declared = {MANIFEST_NAME, INVENTORY_NAME, FREEZE_RECEIPT_NAME}
    declared.add(parsed["candidate_binding"]["clean_candidate_manifest"]["path"])
    declared.update(binding["path"] for binding in parsed["artifact_bindings"].values())
    declared.update(product["path"] for product in parsed["products"])
    _require(
        len(declared)
        == 3 + 1 + len(parsed["artifact_bindings"]) + len(parsed["products"]),
        "duplicate_declared_path",
        "campaign manifest contains duplicate retained paths",
    )
    measured = snapshot.relative_files()
    _require(
        measured == declared,
        "artifact_inventory_drift",
        "campaign tree membership does not match declared retained files",
    )


def _validate_product_hashes(
    products: Sequence[Mapping[str, Any]], snapshot: Any
) -> None:
    for product in products:
        _member_sha256(snapshot, product, f"{product['kind']} product {product['path']}")


def _parse_pvtk_execution_header(payload: bytes, label: str) -> dict[str, Any]:
    match = _PVTK_EXECUTION_PATTERN.match(payload[:4096])
    _require(match is not None, "invalid_prtcl_all", f"{label}: malformed execution header")
    try:
        time = float(match.group(1))
        nranks = int(match.group(2))
        cycle = int(match.group(3))
        variables = match.group(4).decode("ascii")
    except (UnicodeDecodeError, ValueError) as error:
        _fail("invalid_prtcl_all", f"{label}: invalid execution header: {error}")
    _require(
        math.isfinite(time) and nranks > 0 and cycle >= 0 and variables == "prtcl_all",
        "invalid_prtcl_all",
        f"{label}: execution header binding drifted",
    )
    return {"time": time, "nranks": nranks, "cycle": cycle, "variables": variables}


def _validate_particle_snapshot(path: Path, payload: bytes, label: str) -> dict[str, Any]:
    execution = _parse_pvtk_execution_header(payload, label)
    try:
        data: ParticleVTKData = read_particle_vtk(path)
    except (OSError, ValueError) as error:
        _fail("invalid_prtcl_all", f"{label}: particle VTK decode failed: {error}")
    _require(
        set(data.scalars) == _PVTK_SCALARS,
        "invalid_prtcl_all",
        f"{label}: particle scalar inventory drifted",
    )
    _require(
        set(data.vectors) == {"vel"},
        "invalid_prtcl_all",
        f"{label}: particle vector inventory drifted",
    )
    _require(
        data.points.shape == data.vectors["vel"].shape
        and data.points.ndim == 2
        and data.points.shape[1] == 3,
        "invalid_prtcl_all",
        f"{label}: particle point/vector shape drifted",
    )
    _require(
        data.points.shape[0] > 0,
        "invalid_prtcl_all",
        f"{label}: particle payload is empty",
    )
    _require(
        all(
            np.issubdtype(data.scalars[name].dtype, np.integer)
            for name in _PVTK_INTEGER_SCALARS
        )
        and all(
            np.issubdtype(data.scalars[name].dtype, np.floating)
            for name in _PVTK_REAL_SCALARS
        ),
        "invalid_prtcl_all",
        f"{label}: particle scalar integer/real typing drifted",
    )
    _require(
        np.all(np.isfinite(data.points))
        and np.all(np.isfinite(data.vectors["vel"]))
        and all(np.all(np.isfinite(data.scalars[name])) for name in _PVTK_REAL_SCALARS),
        "invalid_prtcl_all",
        f"{label}: particle payload contains non-finite values",
    )
    _require(
        np.all(data.scalars["macro_weight"] >= 0.0)
        and np.all(data.scalars["gid"] >= 0)
        and np.all(data.scalars["ptag"] >= 0)
        and np.all(data.scalars["species"] >= 0)
        and np.all(np.isin(data.scalars["cr_source"], (0, 1))),
        "invalid_prtcl_all",
        f"{label}: particle provenance or macro weight is invalid",
    )
    _require(
        np.unique(data.scalars["ptag"]).size == data.points.shape[0],
        "invalid_prtcl_all",
        f"{label}: particle ptag values are not unique",
    )
    return {
        "particle_count": int(data.points.shape[0]),
        "execution_header": execution,
        "storage_layout": "single_mpi_io_pvtk_file_with_nranks_header",
    }


def _validate_snapshot_payloads(
    retained: Mapping[str, Mapping[str, Mapping[str, Any]]],
    snapshot: Any,
    policy: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    report = {}
    tolerance = _snapshot_tolerance(policy)
    for time, products in retained.items():
        expected_time = float(time)
        cycles: dict[str, int] = {}
        mesh_report = {}
        for kind in ("rho", "bmag", "prtcl_jx", "j2"):
            product = products[kind]
            payload = _member_payload(snapshot, product, f"{kind} snapshot at t={time}")
            try:
                dataset = output_primitives.parse_athenak_binary_bytes(
                    payload, source=product["path"]
                )
            except ValueError as error:
                _fail("invalid_mesh_bin", f"{kind} snapshot at t={time}: {error}")
            _require(
                abs(dataset.time - expected_time) <= tolerance
                and abs(dataset.time - product["snapshot_time"]) <= tolerance,
                "snapshot_metadata_drift",
                f"{kind} snapshot at t={time}: embedded time is {dataset.time!r}",
            )
            _require(
                dataset.variable_names == _EXPECTED_BIN_FIELDS[kind],
                "invalid_mesh_bin",
                f"{kind} snapshot at t={time}: binary field inventory drifted",
            )
            cycles[kind] = dataset.cycle
            mesh_report[kind] = {
                "cycle": dataset.cycle,
                "time": dataset.time,
                "meshblock_count": len(dataset.blocks),
            }
        particle = products["prtcl_all"]
        particle_payload = _member_payload(
            snapshot, particle, f"prtcl_all snapshot at t={time}"
        )
        particle_report = _validate_particle_snapshot(
            snapshot.member_path(particle["path"]),
            particle_payload,
            f"prtcl_all snapshot at t={time}",
        )
        header = particle_report["execution_header"]
        _require(
            abs(header["time"] - expected_time) <= tolerance
            and abs(header["time"] - particle["snapshot_time"]) <= tolerance,
            "snapshot_metadata_drift",
            f"prtcl_all snapshot at t={time}: embedded time is {header['time']!r}",
        )
        cycles["prtcl_all"] = header["cycle"]
        _require(
            len(set(cycles.values())) == 1,
            "snapshot_metadata_drift",
            f"snapshot at t={time}: embedded cycles differ across retained raw products",
        )
        report[time] = {
            "cycle": next(iter(cycles.values())),
            "mesh_bins": mesh_report,
            "particles": particle_report,
        }
    return report


def _validate_stdout_telemetry(
    payload: bytes, run_identity: Mapping[str, Any]
) -> dict[str, Any]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        _fail("invalid_q017_telemetry", f"stdout telemetry is not UTF-8: {error}")
    telemetry = {}
    for match in _Q017_TELEMETRY_PATTERN.finditer(text):
        name = match.group(1)
        _require(
            name not in telemetry,
            "invalid_q017_telemetry",
            f"duplicate Q-017 stdout telemetry name: {name}",
        )
        try:
            telemetry[name] = float(match.group(2))
        except ValueError:
            _fail("invalid_q017_telemetry", f"non-numeric Q-017 telemetry value: {name}")
    missing = sorted(_Q017_REQUIRED_NAMES - set(telemetry))
    _require(
        not missing,
        "invalid_q017_telemetry",
        f"stdout is missing required Q-017 telemetry: {', '.join(missing)}",
    )
    _require(
        all(math.isfinite(value) and value >= 0.0 for value in telemetry.values())
        and telemetry["schema_version"] == 2.0,
        "invalid_q017_telemetry",
        "stdout Q-017 telemetry contains invalid values or schema version",
    )
    expected_runtime_line = "PIC runtime model: " + " ".join(
        f"{name}={value}" for name, value in _EXPECTED_SECTION54_RUNTIME_PROJECTION
    )
    runtime_lines = [
        line for line in text.splitlines() if line.startswith("PIC runtime model:")
    ]
    _require(
        runtime_lines == [expected_runtime_line],
        "invalid_runtime_identity",
        "stdout PIC runtime scientific projection is incomplete, reordered, or drifted",
    )
    try:
        runtime_identity = frozen_model.parse_runtime_identity(text)
    except frozen_model.ModelContractError as error:
        _fail("invalid_runtime_identity", f"stdout PIC runtime identity is invalid: {error}")
    _require(
        runtime_identity.physical_mode == run_identity["physical_mode"],
        "runtime_identity_mismatch",
        "stdout PIC runtime identity differs from the retained run identity",
    )
    return {
        "schema_version": telemetry["schema_version"],
        "mpi_ranks": telemetry["mpi.ranks"],
        "retained_name_count": len(telemetry),
        "pic_runtime_identity": {
            "physical_mode": runtime_identity.physical_mode,
            "state": runtime_identity.state,
            "light_speed": runtime_identity.light_speed,
            "restart_schema": runtime_identity.restart_schema,
        },
    }


def _pairing_schema(
    identity: Mapping[str, Any], policy: Mapping[str, Any]
) -> dict[str, Any]:
    matrix = policy["athenak_selected_release_criteria"]["campaign_matrix"]
    variant = identity["variant"]
    seed = identity["seed"]
    if variant == "three_level_amr_root_dx12_finest_dx3":
        counterpart = {"variant": "fine_uniform_dx3", "seed": seed}
    elif variant == "fine_uniform_dx3":
        counterpart = {"variant": "three_level_amr_root_dx12_finest_dx3", "seed": seed}
    else:
        counterpart = None
    return {
        "paired_seed_rule": matrix["paired_seed_rule"],
        "pair_key": {"seed": seed},
        "amr_fine_uniform_counterpart": counterpart,
        "comparison_status": "schema_wired_not_evaluated_by_artifact_admission_slice",
    }


def spectrum_wiring(policy: Mapping[str, Any]) -> dict[str, Any]:
    """Expose the frozen weighted-spectrum arithmetic contract without inventing data."""
    spectrum = policy["athenak_selected_release_criteria"]["spectrum"]
    return {
        "histogram_primitive": "analyze_q011_section54_outputs.fixed_histogram",
        "slope_fit_primitive": (
            "analyze_q011_section54_outputs.fit_fixed_bin_loglog_slope"
        ),
        "bin_edges": spectrum["bin_edges"],
        "histogram_weight": spectrum["histogram_weight"],
        "slope_fit_window_chi": spectrum["slope_fit_window_chi"],
        "minimum_positive_fit_bins": spectrum["minimum_positive_fit_bins"],
        "runtime_projection_status": (
            "wired_requires_reviewed_chi_values_not_inferred_from_pvtk_velocity"
        ),
    }


def weighted_spectrum_from_chi(
    chi: object, macro_weight: object, policy: Mapping[str, Any]
) -> dict[str, Any]:
    """Apply the preregistered fixed weighted histogram to caller-reviewed chi values."""
    spectrum = policy["athenak_selected_release_criteria"]["spectrum"]
    histogram = output_primitives.fixed_histogram(
        chi,
        spectrum["bin_edges"],
        weights=macro_weight,
    )
    widths = np.diff(histogram.bin_edges)
    f_chi = histogram.weighted_counts / widths
    macro_weight_in_bins = float(np.sum(histogram.weighted_counts))
    admitted_weight = float(
        macro_weight_in_bins + histogram.underflow_weight + histogram.overflow_weight
    )
    plotted = np.zeros_like(f_chi)
    if admitted_weight > 0.0:
        plotted = histogram.geometric_bin_centers * f_chi / admitted_weight
    overflow_fraction = (
        histogram.overflow_weight / admitted_weight if admitted_weight > 0.0 else 0.0
    )
    return {
        "bin_edges": histogram.bin_edges.tolist(),
        "counts": histogram.counts.tolist(),
        "weighted_counts": histogram.weighted_counts.tolist(),
        "f_chi": f_chi.tolist(),
        "normalized_chi_f_chi": plotted.tolist(),
        "underflow_count": histogram.underflow_count,
        "overflow_count": histogram.overflow_count,
        "underflow_weight": histogram.underflow_weight,
        "overflow_weight": histogram.overflow_weight,
        "macro_weight_in_bins": macro_weight_in_bins,
        "admitted_weight": admitted_weight,
        "overflow_macro_weight_fraction": overflow_fraction,
        "overflow_gate": {
            "maximum_macro_weight_fraction": spectrum[
                "max_overflow_macro_weight_fraction"
            ],
            "passed": bool(
                overflow_fraction <= spectrum["max_overflow_macro_weight_fraction"]
            ),
        },
    }


def _validate_freeze_receipt_semantics(receipt: object) -> None:
    require_exact_primitive_types(
        receipt,
        _EXPECTED_FREEZE_RECEIPT,
        error_type=QualificationError,
        label="campaign freeze receipt",
    )
    _require(
        receipt == _EXPECTED_FREEZE_RECEIPT,
        "freeze_receipt_semantics_drift",
        "campaign freeze receipt role or effect semantics drifted",
    )


def _complete_campaign_admission(
    parsed: Mapping[str, Any],
    snapshot: Any,
    policy: Mapping[str, Any],
    retained_attempt_semantics: Mapping[str, Any],
    *,
    frozen_candidate: Mapping[str, Any],
    external_candidate_closure: Mapping[str, Any],
    tree_report: Mapping[str, Any],
) -> dict[str, Any]:
    """Read retained products while the registered-execution ledger stays pinned."""
    retained_snapshot_products = _select_retained_snapshot_products(
        parsed["products"], policy
    )
    endpoint_products = _select_products(retained_snapshot_products, policy)
    stdout_product = _select_stdout_product(parsed["products"])
    _validate_product_hashes(parsed["products"], snapshot)
    snapshot_payloads = _validate_snapshot_payloads(
        retained_snapshot_products, snapshot, policy
    )
    restart_publications = _validate_restart_publications(
        parsed["products"], snapshot, policy
    )
    stdout_telemetry = _validate_stdout_telemetry(
        _member_payload(snapshot, stdout_product, "stdout product"),
        parsed["run_identity"],
    )
    particle_endpoints = {
        time: snapshot_payloads[time]["particles"] for time in endpoint_products
    }
    return {
        "run_identity": parsed["run_identity"],
        "candidate_binding": parsed["candidate_binding"],
        "frozen_clean_candidate": frozen_candidate,
        "external_clean_candidate_closure": external_candidate_closure,
        "artifact_bindings": parsed["artifact_bindings"],
        "attempt_identity": parsed["attempt_identity"],
        "retained_attempt_semantics": retained_attempt_semantics,
        "preregistration_binding": {
            "sha256": parsed["artifact_bindings"]["preregistration"]["sha256"],
            "expected_sha256": EXPECTED_PREREGISTRATION_SHA256,
            "binding_scope": "complete_retained_bytes_equal_invoked_frozen_policy",
        },
        "retained_snapshot_products": retained_snapshot_products,
        "endpoint_products": endpoint_products,
        "run_products": {"stdout": stdout_product},
        "snapshot_payloads": snapshot_payloads,
        "restart_publications": restart_publications,
        "stdout_telemetry": stdout_telemetry,
        "particle_endpoints": particle_endpoints,
        "immutable_tree": {
            "inventory_sha256": tree_report["inventory_sha256"],
            "inventoried_file_count": tree_report["inventoried_file_count"],
            "recursively_read_only": tree_report["recursively_read_only"],
        },
        "weighted_spectrum": spectrum_wiring(policy),
        "amr_pairing": _pairing_schema(parsed["run_identity"], policy),
        "numerical_qualification_status": (
            "not_evaluated_by_artifact_admission_slice"
        ),
    }


def _admit_campaign(
    campaign_root: str | Path,
    expected_inventory_sha256: str,
    *,
    authorized_orion_root: Path,
) -> dict[str, Any]:
    root = authorized_tree_root(
        campaign_root,
        authorized_root=authorized_orion_root,
        error_type=QualificationError,
        label="Q-011 Section 5.4 campaign root",
    )
    with staged_verified_frozen_tree(
        root,
        expected_inventory_sha256,
        authorized_root=authorized_orion_root,
        error_type=QualificationError,
        label="Q-011 Section 5.4 campaign tree",
    ) as (tree_report, snapshot):
        _validate_freeze_receipt_semantics(tree_report["freeze_receipt"])
        manifest = _load_json_bytes(
            snapshot.member_path(MANIFEST_NAME).read_bytes(),
            "campaign manifest",
        )
        parsed = _validate_manifest_schema(manifest, root)
        _validate_declared_tree_members(parsed, snapshot)
        frozen_candidate, external_candidate_closure = _validate_bound_files(parsed, snapshot)
        policy = _load_bound_policy(
            snapshot, parsed["artifact_bindings"]["preregistration"]
        )
        _validate_identity(parsed["run_identity"], parsed["attempt_identity"], policy)
        retained_attempt_snapshot = _validated_retained_attempt_semantics_snapshot(
            parsed,
            snapshot,
            policy,
            authorized_pic_root=authorized_orion_root,
        )
        with retained_attempt_snapshot as retained_attempt_semantics:
            return _complete_campaign_admission(
                parsed,
                snapshot,
                policy,
                retained_attempt_semantics,
                frozen_candidate=frozen_candidate,
                external_candidate_closure=external_candidate_closure,
                tree_report=tree_report,
            )


def qualify_campaign(
    campaign_root: str | Path,
    expected_inventory_sha256: str,
    *,
    authorized_orion_root: Path = ORION_BULK_ROOT,
) -> dict[str, Any]:
    """Return one deterministic admission result without partial success."""
    result = {
        "schema_version": 1,
        "record_type": RESULT_RECORD_TYPE,
        "campaign_id": CAMPAIGN_ID,
        "qualification_scope": QUALIFICATION_SCOPE,
        "admitted_for_follow_on_numerical_qualification": False,
        "final_claim_closure": False,
        "status": "rejected",
        "failure_reasons": [],
        "admission": None,
    }
    try:
        admission = _admit_campaign(
            campaign_root,
            expected_inventory_sha256,
            authorized_orion_root=authorized_orion_root,
        )
    except (QualificationError, OSError, ValueError) as error:
        code = getattr(error, "code", "invalid_campaign")
        result["failure_reasons"] = [{"code": code, "message": str(error)}]
        return result
    result["admitted_for_follow_on_numerical_qualification"] = True
    result["status"] = "admitted_for_follow_on_numerical_qualification"
    result["admission"] = admission
    return result


def _write_json(value: Mapping[str, Any]) -> None:
    print(json.dumps(value, indent=2, sort_keys=True))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign_root")
    parser.add_argument("expected_inventory_sha256")
    args = parser.parse_args(argv)
    result = qualify_campaign(args.campaign_root, args.expected_inventory_sha256)
    _write_json(result)
    return 0 if result["admitted_for_follow_on_numerical_qualification"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
