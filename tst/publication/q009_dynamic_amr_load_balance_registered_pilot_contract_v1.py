#!/usr/bin/env python3
"""Validate the non-authorizing Q009 registered-pilot preregistration.

This module validates checked-in planning bytes only.  It cannot materialize a
deck, mutate policy, inspect runtime output, call a scheduler, or admit evidence.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = (
    REPO_ROOT
    / "tst/publication/q009_dynamic_amr_load_balance_registered_pilot_contract_v1.json"
)
READINESS_PATH = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q009_dynamic_amr_load_balance_registered_pilot_contract_v1_2026-06-06.json"
)

RECORD_TYPE = "q009_dynamic_amr_load_balance_registered_pilot_contract"
CONTRACT_ID = "q009_dynamic_amr_load_balance_registered_pilot_contract_v1"
QUALIFICATION_EFFECT = (
    "registered_engineering_pilot_preregistration_only_no_launch_no_policy_"
    "no_science_no_publication_authority"
)
EXPECTED_CASE_IDS = (
    "q009-amr-migration-lb-off-v1",
    "q009-amr-migration-lb-on-v1",
    "q009-shock-restart-post45-lb-off-v1",
    "q009-shock-restart-post45-lb-on-v1",
    "q009-shock-heldout-full-geometry-lb-on-v1",
)
EXPECTED_CASE_SUMMARIES = {
    "q009-amr-migration-lb-off-v1": (
        1,
        "migration_pair_v1",
        "q009_coupled_boundary_lifetime",
        "reduced_lifetime",
        "repeated_refine_derefine",
        "none",
        "telemetry_only",
        "off",
        0.0,
        2,
        16,
        1,
        3600,
        2.0,
        4.0,
    ),
    "q009-amr-migration-lb-on-v1": (
        2,
        "migration_pair_v1",
        "q009_coupled_boundary_lifetime",
        "reduced_lifetime",
        "repeated_refine_derefine",
        "none",
        "telemetry_only",
        "particle_weight_0p001",
        0.001,
        2,
        16,
        1,
        3600,
        2.0,
        4.0,
    ),
    "q009-shock-restart-post45-lb-off-v1": (
        3,
        "restart_post45_pair_v1",
        "q011_section54_production_preparation",
        "reduced_transverse_shock",
        "shock_curvature_dynamic_amr",
        "checkpoint_before45_resume_after45",
        "telemetry_and_bound_restart",
        "off",
        0.0,
        8,
        64,
        2,
        3600,
        16.0,
        32.0,
    ),
    "q009-shock-restart-post45-lb-on-v1": (
        4,
        "restart_post45_pair_v1",
        "q011_section54_production_preparation",
        "reduced_transverse_shock",
        "shock_curvature_dynamic_amr",
        "checkpoint_before45_resume_after45",
        "telemetry_and_bound_restart",
        "particle_weight_0p001",
        0.001,
        8,
        64,
        2,
        3600,
        16.0,
        32.0,
    ),
    "q009-shock-heldout-full-geometry-lb-on-v1": (
        5,
        None,
        "q011_section54_production_preparation",
        "full_section54_geometry",
        "shock_curvature_dynamic_amr",
        "none",
        "telemetry_only",
        "particle_weight_0p001",
        0.001,
        32,
        256,
        1,
        7200,
        64.0,
        128.0,
    ),
}
EXPECTED_GLOBAL_NODE_HOURS = 100.0
EXPECTED_GLOBAL_STORAGE_GIB = 200.0
EXPECTED_GATES = (
    "final_clean_candidate_exact_source_archive_executable_and_environment_bound",
    "installed_registered_control_plane_generation_exactly_bound",
    "separate_reviewed_registered_science_policy_slice_for_each_case",
    "project_wide_node_hour_ledger_and_other_campaign_reserves_bound",
    "orion_only_artifact_root_and_storage_reservation_bound",
    "user_queue_empty_immediately_before_every_submission",
    "directive_only_normal_qos_job_template_and_fresh_operator_attestations_bound",
    "serial_case_order_and_one_attempt_per_case_enforced",
    "required_q009_telemetry_and_restart_receipt_producer_installed_and_bound",
    "no_runtime_output_bytes_inspected_before_this_preregistration_is_frozen",
)
EXPECTED_ALLOWED_INSPECTION = (
    "scheduler_state_elapsed_and_charged_node_hours",
    "artifact_byte_counts_and_publication_times",
    "aggregate_meshblock_cell_particle_and_migration_counters",
    "aggregate_load_cost_imbalance_and_timing_telemetry",
    "aggregate_conservation_accounting_terms_and_residuals",
    "restart_identity_continuity_and_terminal_time_receipts",
    "hip_memory_high_water_and_failure_telemetry",
)
EXPECTED_FORBIDDEN_INSPECTION = (
    "cell_or_particle_field_values",
    "shock_morphology_or_profiles",
    "particle_spectra_or_acceleration_efficiency",
    "magnetic_amplification_or_instability_growth",
    "scientific_residual_threshold_selection",
    "publication_figures_or_claim_metrics",
    "qualifying_campaign_outputs",
)
_SHA256 = re.compile(r"[0-9a-f]{64}")


class ContractError(ValueError):
    """Reject incomplete, drifted, or authority-bearing Q009 contracts."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def _object(value: object, keys: set[str], *, label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label}: expected object")
    _require(set(value) == keys, f"{label}: keys drifted")
    return value


def _array(value: object, *, label: str) -> list[Any]:
    _require(type(value) is list, f"{label}: expected array")
    return value


def _relative_path(value: object, *, label: str) -> str:
    _require(type(value) is str and bool(value), f"{label}: expected path")
    path = PurePosixPath(value)
    _require(
        not path.is_absolute()
        and path.as_posix() == value
        and value != "."
        and all(part not in ("", ".", "..") for part in path.parts),
        f"{label}: unsafe repository-relative path",
    )
    return value


def _sha256(value: object, *, label: str) -> str:
    _require(
        type(value) is str and _SHA256.fullmatch(value) is not None,
        f"{label}: malformed SHA-256",
    )
    return value


def _stable_regular_bytes(path: Path, *, label: str) -> bytes:
    lexical = Path(os.path.abspath(path))
    _require(path.is_absolute(), f"{label}: expected absolute path")
    try:
        resolved = lexical.resolve(strict=True)
    except OSError as error:
        raise ContractError(f"{label}: unavailable") from error
    _require(resolved == lexical, f"{label}: symlink or path alias forbidden")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(lexical, flags)
    except OSError as error:
        raise ContractError(f"{label}: not an openable regular file") from error
    try:
        before = os.fstat(descriptor)
        _require(stat.S_ISREG(before.st_mode), f"{label}: expected regular file")
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        current = os.stat(lexical, follow_symlinks=False)
        identity = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns")
        _require(
            all(getattr(before, key) == getattr(after, key) for key in identity)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and len(payload) == after.st_size,
            f"{label}: changed while reading",
        )
        return bytes(payload)
    finally:
        os.close(descriptor)


def _decode_json(payload: bytes, *, label: str) -> Any:
    def reject_constant(value: str) -> None:
        raise ContractError(f"{label}: forbidden JSON constant {value}")

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            _require(key not in result, f"{label}: duplicate key {key!r}")
            result[key] = value
        return result

    try:
        return json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractError(f"{label}: invalid UTF-8 JSON") from error


def _validate_binding(
    value: object, *, label: str, verify_repo_bindings: bool
) -> dict[str, Any]:
    binding = _object(value, {"path", "sha256", "role"}, label=label)
    relative = _relative_path(binding["path"], label=f"{label}.path")
    digest = _sha256(binding["sha256"], label=f"{label}.sha256")
    _require(type(binding["role"]) is str and bool(binding["role"]), f"{label}: role missing")
    if verify_repo_bindings:
        observed = hashlib.sha256(
            _stable_regular_bytes(REPO_ROOT / relative, label=f"{label}.path")
        ).hexdigest()
        _require(observed == digest, f"{label}: SHA-256 drifted")
    return binding


def _validate_profiles(value: object) -> None:
    profiles = _object(
        value,
        {
            "geometry_profiles",
            "dynamic_amr_profiles",
            "restart_profiles",
            "output_profiles",
        },
        label="profile_definitions",
    )
    geometry = _object(
        profiles["geometry_profiles"],
        {"reduced_lifetime", "reduced_transverse_shock", "full_section54_geometry"},
        label="geometry_profiles",
    )
    reduced_lifetime = _object(
        geometry["reduced_lifetime"],
        {"classification", "template_geometry", "purpose"},
        label="reduced_lifetime",
    )
    _require(
        reduced_lifetime["classification"] == "reduced_engineering_geometry"
        and reduced_lifetime["template_geometry"] == "unchanged",
        "reduced lifetime geometry drifted",
    )
    reduced_shock = _object(
        geometry["reduced_transverse_shock"],
        {"classification", "exact_overrides", "purpose"},
        label="reduced_transverse_shock",
    )
    _require(
        reduced_shock["classification"] == "reduced_engineering_geometry"
        and reduced_shock["exact_overrides"]
        == {"mesh/nx2": 20, "mesh/x2max": 240.0, "time/tlim": 50.0},
        "reduced shock geometry drifted",
    )
    full = _object(
        geometry["full_section54_geometry"],
        {"classification", "exact_geometry", "exact_overrides", "held_out"},
        label="full_section54_geometry",
    )
    _require(
        full["classification"] == "full_production_geometry_engineering_pilot"
        and full["exact_geometry"]
        == {
            "mesh/nx1": 4000,
            "mesh/nx2": 260,
            "mesh/x1max": 48000.0,
            "mesh/x2max": 3120.0,
            "meshblock/nx1": 20,
            "meshblock/nx2": 20,
        }
        and full["exact_overrides"] == {"time/nlim": 512}
        and full["held_out"] is True,
        "full held-out geometry drifted",
    )

    amr = _object(
        profiles["dynamic_amr_profiles"],
        {"repeated_refine_derefine", "shock_curvature_dynamic_amr"},
        label="dynamic_amr_profiles",
    )
    repeated = _object(
        amr["repeated_refine_derefine"],
        {
            "dynamic_amr_required",
            "minimum_refine_events",
            "minimum_derefine_events",
            "inter_rank_particle_migration_required",
            "rank_ownership_change_required",
            "load_balance_event_required",
        },
        label="repeated_refine_derefine",
    )
    _require(
        repeated
        == {
            "dynamic_amr_required": True,
            "minimum_refine_events": 3,
            "minimum_derefine_events": 3,
            "inter_rank_particle_migration_required": True,
            "rank_ownership_change_required": True,
            "load_balance_event_required": True,
        },
        "repeated AMR profile drifted",
    )
    shock = _object(
        amr["shock_curvature_dynamic_amr"],
        {
            "dynamic_amr_required",
            "minimum_refine_events",
            "minimum_derefine_events",
            "inter_rank_particle_migration_required",
            "rank_ownership_change_required",
            "load_balance_event_required",
        },
        label="shock_curvature_dynamic_amr",
    )
    _require(
        shock["dynamic_amr_required"] is True
        and shock["minimum_refine_events"] >= 1
        and shock["minimum_derefine_events"] >= 0
        and shock["inter_rank_particle_migration_required"] is True
        and shock["rank_ownership_change_required"] is True
        and shock["load_balance_event_required"] is True,
        "shock dynamic-AMR profile drifted",
    )

    restart = _object(
        profiles["restart_profiles"],
        {"none", "checkpoint_before45_resume_after45"},
        label="restart_profiles",
    )
    _require(
        restart["none"] == {"restart_required": False},
        "no-restart profile drifted",
    )
    post45 = _object(
        restart["checkpoint_before45_resume_after45"],
        {
            "restart_required",
            "checkpoint_time_min_inclusive",
            "checkpoint_time_max_exclusive",
            "resumed_terminal_time_min_exclusive",
            "resumed_terminal_time_max_inclusive",
            "distinct_scheduler_job_ids_required",
            "checkpoint_and_restart_receipts_byte_bound",
            "restart_boundary_state_fingerprint_exact",
            "restart_boundary_integer_counters_exact",
        },
        label="checkpoint_before45_resume_after45",
    )
    _require(
        post45["restart_required"] is True
        and 0.0 <= post45["checkpoint_time_min_inclusive"]
        < post45["checkpoint_time_max_exclusive"]
        <= 45.0
        <= post45["resumed_terminal_time_min_exclusive"]
        <= post45["resumed_terminal_time_max_inclusive"]
        and all(
            post45[key] is True
            for key in (
                "distinct_scheduler_job_ids_required",
                "checkpoint_and_restart_receipts_byte_bound",
                "restart_boundary_state_fingerprint_exact",
                "restart_boundary_integer_counters_exact",
            )
        ),
        "post-45 restart profile drifted",
    )

    output = _object(
        profiles["output_profiles"],
        {"telemetry_only", "telemetry_and_bound_restart"},
        label="output_profiles",
    )
    for name, profile in output.items():
        checked = _object(
            profile,
            {
                "physical_field_outputs",
                "particle_dump_outputs",
                "aggregate_telemetry",
                "restart_output",
                "artifact_inventory_required",
            },
            label=f"output_profiles.{name}",
        )
        _require(
            checked["physical_field_outputs"] == "disabled"
            and checked["particle_dump_outputs"] == "disabled"
            and checked["aggregate_telemetry"] == "required"
            and checked["artifact_inventory_required"] is True,
            f"{name}: output inspection boundary drifted",
        )
    _require(
        output["telemetry_only"]["restart_output"] == "disabled"
        and output["telemetry_and_bound_restart"]["restart_output"]
        == "required_immutable_and_cross_linked",
        "restart output profiles drifted",
    )


def _case_signature(case: dict[str, Any]) -> tuple[object, ...]:
    return tuple(
        case[key]
        for key in (
            "role",
            "template_id",
            "geometry_profile_id",
            "dynamic_amr_profile_id",
            "restart_profile_id",
            "output_profile_id",
            "random_seed",
            "resources",
            "ceilings",
            "execution_boundary",
        )
    )


def validate_contract(value: object, *, verify_repo_bindings: bool = True) -> dict[str, Any]:
    """Validate one in-memory Q009 contract and return it unchanged."""

    contract = _object(
        value,
        {
            "record_type",
            "schema_version",
            "contract_id",
            "date",
            "gate",
            "status",
            "qualification_effect",
            "scope",
            "predecessor_evidence",
            "source_templates",
            "scientific_boundary",
            "registered_execution_gates",
            "inspection_policy",
            "required_telemetry",
            "profile_definitions",
            "budget_and_storage_envelope",
            "pilot_matrix",
            "pairwise_acceptance",
            "held_out_acceptance",
            "authority_boundary",
            "open_dependencies",
        },
        label="contract",
    )
    _require(contract["record_type"] == RECORD_TYPE, "record type drifted")
    _require(contract["schema_version"] == 1, "schema version drifted")
    _require(contract["contract_id"] == CONTRACT_ID, "contract id drifted")
    _require(contract["date"] == "2026-06-06" and contract["gate"] == "Q-009", "identity drifted")
    _require(
        contract["status"] == "preregistered_non_authorizing_execution_blocked"
        and contract["qualification_effect"] == QUALIFICATION_EFFECT,
        "status or qualification effect drifted",
    )
    _require(type(contract["scope"]) is str and "HIP" in contract["scope"], "scope drifted")

    predecessors = _object(
        contract["predecessor_evidence"],
        {"amr_policy", "coupled_boundary_local", "coupled_inflow_repaired_local"},
        label="predecessor_evidence",
    )
    for name, binding in predecessors.items():
        _validate_binding(binding, label=f"predecessor_evidence.{name}", verify_repo_bindings=verify_repo_bindings)
    templates = _object(
        contract["source_templates"],
        {"q009_coupled_boundary_lifetime", "q011_section54_production_preparation"},
        label="source_templates",
    )
    for name, binding in templates.items():
        _validate_binding(binding, label=f"source_templates.{name}", verify_repo_bindings=verify_repo_bindings)

    boundary = _object(
        contract["scientific_boundary"],
        {
            "uniform_mesh_exact_conservation_qualifies_dynamic_amr",
            "paper_smooth_amr_individually_conservative_claimed",
            "engineering_pilot_can_close_scientific_amr",
            "full_geometry_pilot_is_production_science",
            "required_interpretation",
        },
        label="scientific_boundary",
    )
    _require(
        all(boundary[key] is False for key in boundary if key != "required_interpretation")
        and "does not qualify AMR" in boundary["required_interpretation"],
        "scientific boundary drifted",
    )

    gates = tuple(_array(contract["registered_execution_gates"], label="registered_execution_gates"))
    _require(gates == EXPECTED_GATES, "registered execution gates drifted")
    inspection = _object(
        contract["inspection_policy"],
        {
            "pre_preregistration_runtime_output_inspection",
            "post_preregistration_allowed_inspection",
            "forbidden_inspection",
            "same_contract_output_driven_tuning",
            "held_out_opening_rule",
        },
        label="inspection_policy",
    )
    allowed = tuple(
        _array(
            inspection["post_preregistration_allowed_inspection"],
            label="post_preregistration_allowed_inspection",
        )
    )
    forbidden = tuple(_array(inspection["forbidden_inspection"], label="forbidden_inspection"))
    _require(
        inspection["pre_preregistration_runtime_output_inspection"] == "forbidden"
        and allowed == EXPECTED_ALLOWED_INSPECTION
        and forbidden == EXPECTED_FORBIDDEN_INSPECTION
        and set(allowed).isdisjoint(forbidden)
        and inspection["same_contract_output_driven_tuning"] == "forbidden_requires_new_version"
        and inspection["held_out_opening_rule"]
        == "open_aggregate_telemetry_only_after_terminal_receipt_and_artifact_seal",
        "inspection policy drifted",
    )

    telemetry = _object(
        contract["required_telemetry"],
        {
            "execution_identity",
            "amr_and_migration",
            "load_balance",
            "conservation_accounting",
            "restart_continuity",
            "hip_resource",
            "fail_closed_rules",
        },
        label="required_telemetry",
    )
    for key, entries in telemetry.items():
        values = _array(entries, label=f"required_telemetry.{key}")
        _require(values and len(values) == len(set(values)), f"required_telemetry.{key}: drifted")
    _require(
        "unexplained_particle_loss_or_duplicate_count_must_equal_zero"
        in telemetry["fail_closed_rules"]
        and "missing_or_nonfinite_conservation_accounting_fails"
        in telemetry["fail_closed_rules"]
        and "no_scientific_conservation_tolerance_is_selected_by_this_pilot"
        in telemetry["fail_closed_rules"],
        "conservation fail-closed rules drifted",
    )
    _validate_profiles(contract["profile_definitions"])

    cases = _array(contract["pilot_matrix"], label="pilot_matrix")
    _require(
        tuple(case.get("case_id") for case in cases if type(case) is dict) == EXPECTED_CASE_IDS,
        "pilot matrix case ids or order drifted",
    )
    observed: dict[str, dict[str, Any]] = {}
    total_node_hours = 0.0
    total_storage_gib = 0.0
    for case in cases:
        checked = _object(
            case,
            {
                "case_id",
                "sequence",
                "role",
                "pair_id",
                "template_id",
                "geometry_profile_id",
                "dynamic_amr_profile_id",
                "restart_profile_id",
                "output_profile_id",
                "random_seed",
                "load_balance",
                "resources",
                "ceilings",
                "execution_boundary",
            },
            label="pilot case",
        )
        case_id = checked["case_id"]
        _require(case_id in EXPECTED_CASE_SUMMARIES and case_id not in observed, "unexpected case")
        load_balance = _object(
            checked["load_balance"], {"mode", "pic_load_balance_cost_per_particle"}, label=f"{case_id}.load_balance"
        )
        resources = _object(
            checked["resources"],
            {
                "platform",
                "accelerator",
                "qos",
                "nodes",
                "mpi_ranks",
                "gpus_per_node",
                "gpus_per_rank",
                "cpus_per_rank",
                "scheduler_jobs_per_attempt",
                "walltime_seconds_per_scheduler_job",
            },
            label=f"{case_id}.resources",
        )
        ceilings = _object(
            checked["ceilings"],
            {"maximum_attempts", "maximum_node_hours", "maximum_artifact_storage_gib"},
            label=f"{case_id}.ceilings",
        )
        execution = _object(
            checked["execution_boundary"],
            {
                "submission_scope",
                "launch_state",
                "policy_slice_state",
                "serial_order_required",
                "science_output_inspection_authorized",
                "qualification_authority",
            },
            label=f"{case_id}.execution_boundary",
        )
        summary = (
            checked["sequence"],
            checked["pair_id"],
            checked["template_id"],
            checked["geometry_profile_id"],
            checked["dynamic_amr_profile_id"],
            checked["restart_profile_id"],
            checked["output_profile_id"],
            load_balance["mode"],
            load_balance["pic_load_balance_cost_per_particle"],
            resources["nodes"],
            resources["mpi_ranks"],
            resources["scheduler_jobs_per_attempt"],
            resources["walltime_seconds_per_scheduler_job"],
            ceilings["maximum_node_hours"],
            ceilings["maximum_artifact_storage_gib"],
        )
        _require(summary == EXPECTED_CASE_SUMMARIES[case_id], f"{case_id}: summary drifted")
        _require(
            resources["platform"] == "Frontier"
            and resources["accelerator"] == "HIP"
            and resources["qos"] == "normal"
            and resources["gpus_per_node"] == 8
            and resources["gpus_per_rank"] == 1
            and resources["cpus_per_rank"] == 7
            and resources["mpi_ranks"] == resources["nodes"] * resources["gpus_per_node"]
            and resources["mpi_ranks"] > 1,
            f"{case_id}: not an exact multi-rank Frontier HIP allocation",
        )
        expected_jobs = 2 if checked["restart_profile_id"] == "checkpoint_before45_resume_after45" else 1
        _require(
            resources["scheduler_jobs_per_attempt"] == expected_jobs,
            f"{case_id}: scheduler jobs per logical attempt drifted",
        )
        derived_node_hours = (
            resources["nodes"]
            * resources["scheduler_jobs_per_attempt"]
            * resources["walltime_seconds_per_scheduler_job"]
            / 3600.0
        )
        _require(
            ceilings["maximum_attempts"] == 1
            and math.isclose(ceilings["maximum_node_hours"], derived_node_hours)
            and ceilings["maximum_artifact_storage_gib"] > 0.0,
            f"{case_id}: ceilings drifted",
        )
        _require(
            execution
            == {
                "submission_scope": "registered_science",
                "launch_state": "blocked_pending_separate_admission",
                "policy_slice_state": "absent_not_authorized_by_this_contract",
                "serial_order_required": True,
                "science_output_inspection_authorized": False,
                "qualification_authority": False,
            },
            f"{case_id}: execution boundary drifted",
        )
        observed[case_id] = checked
        total_node_hours += ceilings["maximum_node_hours"]
        total_storage_gib += ceilings["maximum_artifact_storage_gib"]

    migration_off, migration_on, restart_off, restart_on, held_out = (
        observed[case_id] for case_id in EXPECTED_CASE_IDS
    )
    _require(
        _case_signature(migration_off) == _case_signature(migration_on)
        and _case_signature(restart_off) == _case_signature(restart_on),
        "paired cases differ outside case identity and load-balance setting",
    )
    _require(
        held_out["pair_id"] is None
        and held_out["geometry_profile_id"] == "full_section54_geometry"
        and held_out["sequence"] == len(cases),
        "held-out case drifted",
    )

    envelope = _object(
        contract["budget_and_storage_envelope"],
        {
            "project_frontier_cap_node_hours",
            "contract_hard_cap_node_hours",
            "contract_hard_cap_artifact_storage_gib",
            "derived_sum_case_node_hour_ceilings",
            "derived_sum_case_storage_ceilings_gib",
            "maximum_attempts_per_case",
            "overrun_action",
            "other_registered_campaign_reserve_rule",
            "project_cap_admission_equation",
        },
        label="budget_and_storage_envelope",
    )
    _require(
        envelope["project_frontier_cap_node_hours"] == 10000.0
        and envelope["contract_hard_cap_node_hours"] == EXPECTED_GLOBAL_NODE_HOURS
        and envelope["contract_hard_cap_artifact_storage_gib"] == EXPECTED_GLOBAL_STORAGE_GIB
        and math.isclose(total_node_hours, EXPECTED_GLOBAL_NODE_HOURS)
        and math.isclose(total_storage_gib, EXPECTED_GLOBAL_STORAGE_GIB)
        and envelope["derived_sum_case_node_hour_ceilings"] == total_node_hours
        and envelope["derived_sum_case_storage_ceilings_gib"] == total_storage_gib
        and envelope["maximum_attempts_per_case"] == 1
        and envelope["overrun_action"] == "stop_fail_closed_no_retry"
        and envelope["other_registered_campaign_reserve_rule"]
        == "must_be_bound_before_first_q009_submission"
        and "10000" in envelope["project_cap_admission_equation"],
        "budget or storage envelope drifted",
    )

    pairwise = _object(
        contract["pairwise_acceptance"],
        {"required_pair_ids", "common_requirements", "load_balance_pair_requirements"},
        label="pairwise_acceptance",
    )
    _require(
        pairwise["required_pair_ids"] == ["migration_pair_v1", "restart_post45_pair_v1"]
        and "same_candidate_executable_deck_seed_resources_and_stage_plan"
        in pairwise["common_requirements"]
        and "only_load_balance_mode_and_case_identity_may_differ"
        in pairwise["common_requirements"]
        and "lb_on_must_change_at_least_one_post_amr_rank_ownership_assignment"
        in pairwise["load_balance_pair_requirements"],
        "pairwise acceptance drifted",
    )
    held = _object(
        contract["held_out_acceptance"],
        {
            "case_id",
            "launch_preconditions",
            "same_contract_retuning_after_pair_telemetry",
            "allowed_postrun_read",
            "qualification_effect",
        },
        label="held_out_acceptance",
    )
    _require(
        held["case_id"] == EXPECTED_CASE_IDS[-1]
        and held["same_contract_retuning_after_pair_telemetry"] == "forbidden"
        and held["allowed_postrun_read"] == "aggregate_preregistered_telemetry_only"
        and held["qualification_effect"] == "engineering_feasibility_only",
        "held-out acceptance drifted",
    )

    authority = _object(
        contract["authority_boundary"],
        {
            "launch_authorized",
            "scheduler_submission_authorized",
            "policy_mutation_authorized",
            "live_policy_slice_created",
            "q009_qualified",
            "dynamic_amr_science_qualified",
            "exact_amr_conservation_claim_authorized",
            "q011_production_authorized",
            "scientific_claim_authorized",
            "publication_authorized",
        },
        label="authority_boundary",
    )
    _require(all(value is False for value in authority.values()), "authority must remain false")
    dependencies = _array(contract["open_dependencies"], label="open_dependencies")
    _require(
        "implement_and_review_exact_registered_deck_materializer_and_admission"
        in dependencies
        and "implement_and_review_required_amr_conservation_and_migration_telemetry"
        in dependencies
        and "create_separate_live_registered_science_policy_slices_after_review"
        in dependencies,
        "open dependencies drifted",
    )
    return contract


def load_contract(path: Path = CONTRACT_PATH) -> dict[str, Any]:
    """Load and validate the checked-in contract without granting authority."""

    payload = _stable_regular_bytes(path, label="Q009 contract")
    return validate_contract(_decode_json(payload, label="Q009 contract"))


def main() -> int:
    contract = load_contract()
    print(
        json.dumps(
            {
                "contract_id": contract["contract_id"],
                "case_count": len(contract["pilot_matrix"]),
                "contract_hard_cap_node_hours": contract[
                    "budget_and_storage_envelope"
                ]["contract_hard_cap_node_hours"],
                "contract_hard_cap_artifact_storage_gib": contract[
                    "budget_and_storage_envelope"
                ]["contract_hard_cap_artifact_storage_gib"],
                "launch_authorized": contract["authority_boundary"]["launch_authorized"],
                "publication_authorized": contract["authority_boundary"][
                    "publication_authorized"
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
