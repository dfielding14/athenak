#!/usr/bin/env python3
"""Materialize non-authorizing Q-011 resource-scaling pilot review artifacts.

The emitted decks and registered-control-plane config fragments are engineering
review inputs only.  This tool never mutates live policy, authorizes a launch,
calls a scheduler, inspects physical outputs, or closes a scientific claim.
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
from typing import Any, Mapping

if __package__:
    from . import q011_section54_model as model
else:
    import q011_section54_model as model


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_ROOT = REPO_ROOT / "tst/publication/readiness"
RESOURCE_PLAN = READINESS_ROOT / "q011_resource_scaling_preproduction_plan_2026-06-06.json"
PRODUCTION_DECK = REPO_ROOT / "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
SELECTED_PRESSURE_RECEIPT = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/PIC/publication/"
    "q011_section54_pressure_selection_receipt.json"
)
AUTHORIZED_OUTPUT_PARENT = REPO_ROOT / "tst/.codex"
DEFAULT_OUTPUT_ROOT = AUTHORIZED_OUTPUT_PARENT / "q011-resource-scaling-preproduction-pilots"

RESOURCE_PLAN_SHA256 = "f6b059183cd3dd1cfbca55ef2d02300d597bac4e9ac5851aecc3ae79abdebfb5"
PRODUCTION_DECK_SHA256 = "0b1cbd62d54027ec81a5f4f5c88d5ee56b86b8cc0cb018c3fbebfb37a11be7b1"
SELECTED_PRESSURE_RECEIPT_SHA256 = (
    "e5c1492cfc67d5cfad0110e7d772bf75f9d6d2e1fd90b0d5d8dd3338965905cc"
)

RECORD_TYPE = "q011_resource_scaling_preproduction_pilot_materialization"
CONFIG_RECORD_TYPE = "q011_resource_scaling_engineering_pilot_registered_config_fragment"
QUALIFICATION_EFFECT = "engineering_resource_measurement_only_no_execution_authorization_no_claim_closure"
AUTHORIZATION_STATE = "non_authorizing_review_fragment_only"
EVIDENCE_CLASS = "q011_section54_resource_scaling_engineering_only"
PHYSICAL_MODE = model.PAPER_PHYSICAL_MODE
RUNTIME_PROFILE = "frontier_minimum_supported"
SELECTED_QOS = "normal"
HARD_PILOT_CAP_NODE_HOURS = 500
FRONTIER_GPUS_PER_NODE = 8
CPUS_PER_GPU_TASK = 7

VARIANTS = (
    "coarse_uniform_dx12",
    "three_level_amr_root_dx12_finest_dx3",
    "fine_uniform_dx3",
)
VARIANT_SHORT = {
    "coarse_uniform_dx12": "coarse",
    "three_level_amr_root_dx12_finest_dx3": "amr",
    "fine_uniform_dx3": "fine",
}
QUALIFYING_SEEDS = tuple(range(23050101, 23050109))
ENGINEERING_SEEDS = (24060601, 24060602, 24060603, 24060604)
PHASE_SEEDS = {
    "phase_1_ladder": 24060601,
    "phase_1_selected_repeat": 24060602,
    "phase_2_reduced_transverse": 24060604,
    "phase_3_held_out": 24060603,
}
NODE_LADDERS = {
    "coarse_uniform_dx12": (4, 8, 16),
    "three_level_amr_root_dx12_finest_dx3": (8, 16, 32),
    "fine_uniform_dx3": (32, 64, 128),
}
ALLOWED_INSPECTION = (
    "scheduler_elapsed_and_node_hours",
    "runtime_timing_telemetry",
    "particle_counts_and_particle_update_counts",
    "active_cell_and_meshblock_counts",
    "memory_high_water_marks",
    "artifact_byte_counts_and_output_publication_times",
)
FORBIDDEN_INSPECTION = (
    "shock_morphology",
    "spectra",
    "magnetic_amplification",
    "scientific_profiles",
    "scientific_residual_metrics",
    "qualifying_campaign_outputs",
)
COMMON_RUNTIME_GATES = (
    "separate_reviewed_registered_science_policy_slice",
    "final_clean_candidate_and_exact_executable_binding",
    "installed_control_plane_and_environment_profile_binding",
    "directive_only_normal_qos_job_template_with_reviewed_nodes_and_walltime",
    "fresh_timeout_margin_and_queue_snapshot",
    "fresh_pre_manifest_and_pre_submit_wrapper_operator_attestations",
    "valid_q011_planner_retention_binding",
    "serial_submission_and_500_consumed_node_hour_cap_enforcement",
)

_HEADER = (
    "# Q-011 RESOURCE-SCALING ENGINEERING PILOT DECK.\n"
    "# NON-AUTHORIZING: separate registered policy and installed-control-plane review required.\n"
    "# PHYSICAL-OUTPUT INSPECTION IS FORBIDDEN; retain resource telemetry and byte counts only.\n"
)
_SAFE_ID = re.compile(r"[a-z0-9][a-z0-9_-]{0,127}")
_ASSIGNMENT = re.compile(r"^(\s*)([A-Za-z0-9_]+)(\s*=\s*)([^#\n]*?)(\s*(?:#.*)?)$")
_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH


class MaterializationError(ValueError):
    """Reject drifted inputs or an unsafe materialization request."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise MaterializationError(message)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json_bytes(value: object) -> bytes:
    try:
        return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
            "utf-8"
        )
    except (TypeError, ValueError) as error:
        raise MaterializationError("materialization contains noncanonical JSON") from error


def _decode_json(payload: bytes, *, label: str) -> Any:
    def reject_constant(value: str) -> None:
        raise MaterializationError(f"{label} contains forbidden JSON constant {value}")

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            _require(key not in result, f"{label} contains duplicate key {key!r}")
            result[key] = value
        return result

    try:
        return json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise MaterializationError(f"{label} is not valid UTF-8 JSON") from error


def _stable_regular_bytes(path: Path, *, label: str) -> bytes:
    lexical = Path(os.path.abspath(path))
    _require(path.is_absolute(), f"{label} must be absolute")
    try:
        resolved = lexical.resolve(strict=True)
    except OSError as error:
        raise MaterializationError(f"{label} is unavailable") from error
    _require(resolved == lexical, f"{label} must not use a symlink or path alias")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(lexical, flags)
    except OSError as error:
        raise MaterializationError(f"{label} is not an openable regular file") from error
    try:
        before = os.fstat(descriptor)
        _require(stat.S_ISREG(before.st_mode), f"{label} must be a regular file")
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        current = os.stat(lexical, follow_symlinks=False)
        stable = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns")
        _require(
            all(getattr(before, field) == getattr(after, field) for field in stable)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and len(payload) == after.st_size,
            f"{label} changed while reading",
        )
        return bytes(payload)
    finally:
        os.close(descriptor)


def _binding(path: str, payload: bytes) -> dict[str, str]:
    return {"path": path, "sha256": _sha256_bytes(payload)}


def _source_inputs(
    *,
    resource_plan: Path = RESOURCE_PLAN,
    selected_pressure_receipt: Path = SELECTED_PRESSURE_RECEIPT,
    production_deck: Path = PRODUCTION_DECK,
) -> tuple[dict[str, Any], bytes, dict[str, Any], bytes, str, bytes]:
    plan_payload = _stable_regular_bytes(resource_plan, label="Q011 resource-scaling plan")
    receipt_payload = _stable_regular_bytes(
        selected_pressure_receipt, label="selected Q011 pressure receipt"
    )
    deck_payload = _stable_regular_bytes(production_deck, label="Q011 production deck")
    _require(_sha256_bytes(plan_payload) == RESOURCE_PLAN_SHA256, "resource-scaling plan SHA-256 drifted")
    _require(
        _sha256_bytes(receipt_payload) == SELECTED_PRESSURE_RECEIPT_SHA256,
        "selected pressure receipt SHA-256 drifted",
    )
    _require(_sha256_bytes(deck_payload) == PRODUCTION_DECK_SHA256, "production deck SHA-256 drifted")

    plan = _decode_json(plan_payload, label="Q011 resource-scaling plan")
    receipt = _decode_json(receipt_payload, label="selected Q011 pressure receipt")
    _require(type(plan) is dict, "resource-scaling plan must be an object")
    _require(type(receipt) is dict, "selected pressure receipt must be an object")
    selected = receipt.get("selected_case")
    _require(
        receipt.get("record_type") == "q011_section54_pressure_selection_receipt"
        and receipt.get("selection_method") == "human_review_only"
        and selected == {"case_id": "ps_p0_1p00", "problem_ps_p0": 1.0},
        "selected pressure receipt does not bind ps_p0_1p00 / problem_ps_p0=1.0",
    )

    bindings = plan.get("source_and_evidence_bindings")
    design = plan.get("preproduction_scaling_pilot_design")
    matrix = plan.get("current_campaign_matrix")
    boundaries = plan.get("historical_and_evidence_boundaries")
    _require(type(bindings) is dict and type(design) is dict, "resource-scaling plan structure drifted")
    _require(type(matrix) is dict and type(boundaries) is dict, "resource-scaling plan boundary drifted")
    _require(
        bindings.get("paper_deck")
        == {
            "path": "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
            "sha256": PRODUCTION_DECK_SHA256,
        },
        "resource-scaling plan production-deck binding drifted",
    )
    _require(
        bindings.get("selected_pressure_receipt")
        == {
            "path": str(SELECTED_PRESSURE_RECEIPT),
            "sha256": SELECTED_PRESSURE_RECEIPT_SHA256,
            "selected_case_id": "ps_p0_1p00",
            "selected_problem_ps_p0": 1.0,
        },
        "resource-scaling plan selected-pressure binding drifted",
    )
    _require(tuple(matrix.get("grid_variants", ())) == VARIANTS, "resource-scaling variants drifted")
    _require(
        tuple(matrix.get("qualifying_seeds", ())) == QUALIFYING_SEEDS,
        "resource-scaling qualifying seeds drifted",
    )
    seed_policy = design.get("engineering_seed_policy")
    _require(type(seed_policy) is dict, "engineering seed policy is missing")
    _require(tuple(seed_policy.get("seeds", ())) == ENGINEERING_SEEDS, "engineering seed policy drifted")
    _require(set(ENGINEERING_SEEDS).isdisjoint(QUALIFYING_SEEDS), "engineering seeds overlap qualifying seeds")
    _require(design.get("hard_consumed_node_hour_cap") == HARD_PILOT_CAP_NODE_HOURS, "pilot cap drifted")
    _require(
        boundaries.get("qualifying_output_inspection") == "forbidden"
        and boundaries.get("engineering_pilot_physical_field_inspection")
        == "forbidden_resource_telemetry_scheduler_accounting_and_artifact_byte_counts_only",
        "physical-output inspection boundary drifted",
    )
    phase1 = design.get("phase_1_full_geometry_short_step_strong_scaling")
    phase2 = design.get("phase_2_reduced_transverse_full_time_calibration")
    phase3 = design.get("phase_3_held_out_cost_model_validation")
    _require(type(phase1) is dict and type(phase2) is dict and type(phase3) is dict, "pilot phase design drifted")
    _require(phase1.get("cycle_limit") == 512, "Phase 1 cycle limit drifted")
    _require(
        {key: tuple(value) for key, value in phase1.get("node_ladders", {}).items()} == NODE_LADDERS,
        "Phase 1 node ladders drifted",
    )
    _require("engineering seed 24060604" in str(phase2.get("time_and_outputs")), "Phase 2 seed drifted")
    _require("engineering seed 24060603" in str(phase3.get("run")), "Phase 3 seed drifted")
    try:
        deck_text = deck_payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise MaterializationError("production deck is not UTF-8") from error
    base = model.parse_deck_contract(deck_text)
    _require(base.variant == "three_level_amr_root_dx12_finest_dx3", "production deck is not the AMR base")
    return plan, plan_payload, receipt, receipt_payload, deck_text, deck_payload


def _parse_deck(text: str) -> dict[str, dict[str, str]]:
    try:
        return model._parse_athinput(text)
    except model.ModelContractError as error:
        raise MaterializationError(f"materialized deck is malformed: {error}") from error


def _render_deck(base_text: str, overrides: Mapping[tuple[str, str], str]) -> str:
    remaining = dict(overrides)
    current: str | None = None
    rendered: list[str] = []
    for raw in base_text.splitlines():
        stripped = raw.split("#", 1)[0].strip()
        if stripped.startswith("<") and stripped.endswith(">"):
            current = stripped[1:-1].strip()
            rendered.append(raw)
            continue
        match = _ASSIGNMENT.fullmatch(raw)
        if match is None or current is None:
            rendered.append(raw)
            continue
        prefix, name, separator, _, suffix = match.groups()
        key = (current, name)
        if key not in remaining:
            rendered.append(raw)
            continue
        rendered.append(f"{prefix}{name}{separator}{remaining.pop(key)}{suffix}")
    _require(not remaining, f"deck overrides reference missing parameters: {sorted(remaining)}")
    return _HEADER + "\n".join(rendered) + "\n"


def _phase_overrides(phase: str, variant: str, seed: int, basename: str) -> dict[tuple[str, str], str]:
    binding = model.variant_binding(variant)
    overrides: dict[tuple[str, str], str] = {
        ("job", "basename"): basename,
        ("problem", "ps_p0"): "1.0",
        ("particles", "pic_random_seed"): str(seed),
        ("problem", "ps_inject_seed"): str(seed),
        ("problem", "ps_seed_noise_seed"): str(seed),
    }
    for value in binding.model_launch_overrides:
        key, replacement = value.split("=", 1)
        block, name = key.split("/", 1)
        overrides[(block, name)] = replacement
    if phase in {"phase_1_ladder", "phase_1_selected_repeat", "phase_3_held_out"}:
        overrides[("time", "nlim")] = "512"
        for output in range(1, 7):
            overrides[(f"output{output}", "dt")] = "-1.0"
    elif phase == "phase_2_reduced_transverse":
        overrides[("mesh", "x2max")] = "240.0"
        overrides[("mesh", "nx2")] = "80" if variant == "fine_uniform_dx3" else "20"
    else:
        raise MaterializationError(f"unsupported pilot phase: {phase}")
    return overrides


def _validate_rendered_deck(
    base_text: str,
    rendered: str,
    *,
    phase: str,
    variant: str,
    seed: int,
    basename: str,
) -> None:
    base = _parse_deck(base_text)
    actual = _parse_deck(rendered)
    expected = {block: dict(values) for block, values in base.items()}
    overrides = _phase_overrides(phase, variant, seed, basename)
    for (block, name), value in overrides.items():
        expected[block][name] = value
    _require(actual == expected, f"{phase}/{variant} deck differs outside exact engineering overrides")
    _require(actual["problem"]["ps_p0"] == "1.0", "materialized deck pressure drifted")
    actual_seeds = {
        int(actual["particles"]["pic_random_seed"]),
        int(actual["problem"]["ps_inject_seed"]),
        int(actual["problem"]["ps_seed_noise_seed"]),
    }
    _require(actual_seeds == {seed}, "materialized deck engineering seed drifted")
    _require(seed not in QUALIFYING_SEEDS, "materialized deck uses a qualifying seed")
    if phase in {"phase_1_ladder", "phase_1_selected_repeat", "phase_3_held_out"}:
        _require(actual["time"]["nlim"] == "512", "short-step deck nlim drifted")
        _require(
            all(float(actual[f"output{index}"]["dt"]) <= 0.0 for index in range(1, 7)),
            "short-step bulk outputs are not disabled",
        )
        _require(actual["mesh"]["x2max"] == "3120.0", "short-step deck is not full geometry")
    else:
        _require(actual["time"]["tlim"] == "1200.0", "reduced-transverse deck tlim drifted")
        _require(actual["time"]["nlim"] == "200000", "reduced-transverse deck nlim drifted")
        _require(
            all(actual[f"output{index}"]["dt"] == "100.0" for index in range(1, 7)),
            "reduced-transverse deck does not retain production output cadence",
        )
        _require(actual["mesh"]["x2max"] == "240.0", "reduced-transverse x2 extent drifted")


def _deck_id(phase: str, variant: str) -> str:
    short = VARIANT_SHORT[variant]
    value = f"{phase.replace('phase_', 'p').replace('_', '-')}-{short}-seed-{PHASE_SEEDS[phase]}"
    _require(_SAFE_ID.fullmatch(value) is not None, "generated deck ID is unsafe")
    return value


def _case_id(phase_token: str, variant: str, nodes: int | None = None) -> str:
    value = f"{phase_token}-{VARIANT_SHORT[variant]}"
    if nodes is not None:
        value += f"-n{nodes:03d}"
    _require(_SAFE_ID.fullmatch(value) is not None, "generated case ID is unsafe")
    return value


def _launch_contract(case_id: str, nodes: int) -> dict[str, object]:
    _require(type(nodes) is int and nodes > 0, "launch-contract nodes must be positive")
    return {
        "schema_version": 1,
        "executor": "trusted_trampoline_athena_argv_v1",
        "pre_actions": [],
        "actions": [
            {
                "action_id": case_id,
                "kind": "athena",
                "resources": {
                    "nodes": nodes,
                    "tasks": nodes * FRONTIER_GPUS_PER_NODE,
                    "cpus_per_task": CPUS_PER_GPU_TASK,
                    "gpus_per_task": 1,
                    "gpu_bind": "closest",
                },
                "arguments": [
                    {"literal": "-i"},
                    {"snapshot_role": "input-deck"},
                    {"literal": "-d"},
                    {"artifact_directory": "output"},
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


def _config_fragment(
    *,
    case_id: str,
    phase: str,
    variant: str,
    seed: int,
    deck_binding: Mapping[str, str],
    nodes: int | None,
    maximum_attempts: int,
    execution_condition: Mapping[str, object],
) -> dict[str, object]:
    authorization_id = f"q011-rs-{case_id}-v1"
    campaign = f"q011_rs_{case_id.replace('-', '_')}"
    test_id = f"pic_parallel_shock_section54_rs_{case_id.replace('-', '_')}"
    concrete_contract = _launch_contract(case_id, nodes) if nodes is not None else None
    node_binding: dict[str, object]
    if nodes is None:
        node_binding = {
            "status": "blocked_pending_phase_1_selection_and_memory_meshblock_check",
            "nodes": None,
            "tasks": None,
            "selection_rule": (
                "bind selected Phase 1 node count; Phase 2 may scale down only after "
                "a separately reviewed memory and MeshBlock-per-rank check"
            ),
        }
    else:
        node_binding = {
            "status": "concrete_review_candidate",
            "nodes": nodes,
            "tasks": nodes * FRONTIER_GPUS_PER_NODE,
            "gpus_per_task": 1,
            "cpus_per_task": CPUS_PER_GPU_TASK,
        }
    policy_slice = {
        "authorization_id": authorization_id,
        "status": "review_required_not_authorized",
        "campaign": campaign,
        "test_id": test_id,
        "evidence_class": EVIDENCE_CLASS,
        "physical_mode": PHYSICAL_MODE,
        "runtime_profile": RUNTIME_PROFILE,
        "selected_qos": SELECTED_QOS,
        "registered_short_nonproduction": False,
        "maximum_nodes": nodes,
        "maximum_walltime_seconds": None,
        "maximum_attempts": maximum_attempts,
        "input_deck_sha256": deck_binding["sha256"],
        "launch_contract_sha256": (
            _sha256_bytes(
                json.dumps(
                    concrete_contract, sort_keys=True, separators=(",", ":"), allow_nan=False
                ).encode("utf-8")
            )
            if concrete_contract is not None
            else None
        ),
        "unresolved_policy_bindings": [
            "job_script_sha256",
            "environment_profile_sha256",
            "analysis_script_sha256",
            "executable_sha256",
            "clean_candidate_manifest_sha256",
            "maximum_walltime_seconds",
            *(["maximum_nodes", "launch_contract_sha256"] if nodes is None else []),
        ],
    }
    return {
        "record_type": CONFIG_RECORD_TYPE,
        "schema_version": 1,
        "qualification_effect": QUALIFICATION_EFFECT,
        "authorization_state": AUTHORIZATION_STATE,
        "case_id": case_id,
        "phase": phase,
        "variant": variant,
        "engineering_seed": seed,
        "qualifying_seed": False,
        "hard_consumed_node_hour_cap": HARD_PILOT_CAP_NODE_HOURS,
        "submission_scope": "registered_science",
        "selected_qos": SELECTED_QOS,
        "physical_output_inspection": "forbidden",
        "allowed_inspection": list(ALLOWED_INSPECTION),
        "forbidden_inspection": list(FORBIDDEN_INSPECTION),
        "input_deck": dict(deck_binding),
        "node_binding": node_binding,
        "execution_condition": dict(execution_condition),
        "registered_policy_slice_fragment": policy_slice,
        "launch_contract_candidate": concrete_contract,
        "launch_contract_template_gate": (
            None
            if concrete_contract is not None
            else "materialize_and_review_only_after_phase_1_node_selection"
        ),
        "required_runtime_gates": list(COMMON_RUNTIME_GATES),
        "planner_retention": {
            "required_by_installed_q011_registered_control_plane": True,
            "materialized_by_this_artifact": False,
            "gate": "separate_reviewed_pilot_planner_retention_binding_required",
        },
        "execution_boundary": {
            "mutates_live_policy": False,
            "scheduler_calls": False,
            "submits_jobs": False,
            "launch_authorized": False,
            "frontier_execution_authorized": False,
            "physical_output_inspection_authorized": False,
            "qualifying_evidence": False,
            "claim_closure_authorized": False,
        },
    }


def build_materialization(
    *,
    resource_plan: Path = RESOURCE_PLAN,
    selected_pressure_receipt: Path = SELECTED_PRESSURE_RECEIPT,
    production_deck: Path = PRODUCTION_DECK,
) -> tuple[dict[str, object], dict[str, bytes]]:
    """Build deterministic manifest and member bytes without writing them."""
    plan, plan_payload, _, receipt_payload, base_text, deck_payload = _source_inputs(
        resource_plan=resource_plan,
        selected_pressure_receipt=selected_pressure_receipt,
        production_deck=production_deck,
    )
    files: dict[str, bytes] = {
        "bindings/q011_resource_scaling_preproduction_plan.json": plan_payload,
        "bindings/q011_section54_pressure_selection_receipt.json": receipt_payload,
        "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput": deck_payload,
    }
    deck_bindings: dict[tuple[str, str], dict[str, str]] = {}
    for phase in PHASE_SEEDS:
        for variant in VARIANTS:
            identifier = _deck_id(phase, variant)
            rendered = _render_deck(
                base_text,
                _phase_overrides(phase, variant, PHASE_SEEDS[phase], identifier),
            )
            _validate_rendered_deck(
                base_text,
                rendered,
                phase=phase,
                variant=variant,
                seed=PHASE_SEEDS[phase],
                basename=identifier,
            )
            relative = f"decks/{identifier}.athinput"
            payload = rendered.encode("utf-8")
            files[relative] = payload
            deck_bindings[(phase, variant)] = _binding(relative, payload)

    configs: list[dict[str, str]] = []

    def add_config(value: dict[str, object]) -> None:
        relative = f"config_fragments/{value['case_id']}.json"
        payload = _json_bytes(value)
        _require(relative not in files, f"duplicate materialized config path: {relative}")
        files[relative] = payload
        configs.append(_binding(relative, payload))

    for variant in VARIANTS:
        for nodes in NODE_LADDERS[variant]:
            add_config(
                _config_fragment(
                    case_id=_case_id("p1-ladder", variant, nodes),
                    phase="phase_1_full_geometry_short_step_strong_scaling",
                    variant=variant,
                    seed=PHASE_SEEDS["phase_1_ladder"],
                    deck_binding=deck_bindings[("phase_1_ladder", variant)],
                    nodes=nodes,
                    maximum_attempts=1,
                    execution_condition={"status": "initial_phase_1_ladder_point"},
                )
            )
            add_config(
                _config_fragment(
                    case_id=_case_id("p1-repeat", variant, nodes),
                    phase="phase_1_selected_node_repeat",
                    variant=variant,
                    seed=PHASE_SEEDS["phase_1_selected_repeat"],
                    deck_binding=deck_bindings[("phase_1_selected_repeat", variant)],
                    nodes=nodes,
                    maximum_attempts=3,
                    execution_condition={
                        "status": "conditional",
                        "run_only_if_phase_1_selected_node_count_equals": nodes,
                        "selection_rule": (
                            "smallest node count whose median node-hours per completed cycle "
                            "is within 10 percent of the observed minimum and passes memory "
                            "and projected-walltime gates"
                        ),
                    },
                )
            )
        add_config(
            _config_fragment(
                case_id=_case_id("p2-reduced", variant),
                phase="phase_2_reduced_transverse_full_time_calibration",
                variant=variant,
                seed=PHASE_SEEDS["phase_2_reduced_transverse"],
                deck_binding=deck_bindings[("phase_2_reduced_transverse", variant)],
                nodes=None,
                maximum_attempts=1,
                execution_condition={
                    "status": "blocked",
                    "requires_passed_gates": ["Q011-RS-0", "Q011-RS-1"],
                },
            )
        )
        add_config(
            _config_fragment(
                case_id=_case_id("p3-heldout", variant),
                phase="phase_3_held_out_cost_model_validation",
                variant=variant,
                seed=PHASE_SEEDS["phase_3_held_out"],
                deck_binding=deck_bindings[("phase_3_held_out", variant)],
                nodes=None,
                maximum_attempts=1,
                execution_condition={
                    "status": "blocked",
                    "requires_passed_gates": ["Q011-RS-0", "Q011-RS-1", "Q011-RS-2"],
                    "held_out_from_cost_model_fit": True,
                },
            )
        )

    deck_records = [
        {"phase": phase, "variant": variant, **deck_bindings[(phase, variant)]}
        for phase in PHASE_SEEDS
        for variant in VARIANTS
    ]
    manifest = {
        "record_type": RECORD_TYPE,
        "schema_version": 1,
        "qualification_effect": QUALIFICATION_EFFECT,
        "authorization_state": AUTHORIZATION_STATE,
        "status": "source_local_review_artifact_only",
        "source_bindings": {
            "resource_plan": _binding(
                "bindings/q011_resource_scaling_preproduction_plan.json", plan_payload
            ),
            "selected_pressure_receipt": _binding(
                "bindings/q011_section54_pressure_selection_receipt.json", receipt_payload
            ),
            "production_deck": _binding(
                "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput", deck_payload
            ),
        },
        "selected_pressure": {"case_id": "ps_p0_1p00", "problem_ps_p0": 1.0},
        "engineering_seed_allocation": dict(PHASE_SEEDS),
        "qualifying_seeds": list(QUALIFYING_SEEDS),
        "engineering_seeds_disjoint_from_qualifying_seeds": True,
        "hard_consumed_node_hour_cap": HARD_PILOT_CAP_NODE_HOURS,
        "physical_output_inspection": "forbidden",
        "allowed_inspection": list(ALLOWED_INSPECTION),
        "forbidden_inspection": list(FORBIDDEN_INSPECTION),
        "phase_1_node_ladders": {key: list(value) for key, value in NODE_LADDERS.items()},
        "deck_count": len(deck_records),
        "decks": deck_records,
        "config_fragment_count": len(configs),
        "config_fragments": configs,
        "required_runtime_gates": list(COMMON_RUNTIME_GATES),
        "decision_gates": plan["decision_gates"],
        "execution_boundary": {
            "historical_preregistration_mutation": False,
            "mutates_live_policy": False,
            "scheduler_calls": False,
            "submits_jobs": False,
            "launch_authorized": False,
            "frontier_execution_authorized": False,
            "physical_output_inspection_authorized": False,
            "qualifying_evidence": False,
            "claim_closure_authorized": False,
        },
    }
    return manifest, files


def _inventory_payload(files: Mapping[str, bytes]) -> bytes:
    lines = [f"{_sha256_bytes(files[path])}  {path}\n" for path in sorted(files)]
    return "".join(lines).encode("utf-8")


def _safe_output_root(output_root: Path) -> Path:
    root = Path(output_root)
    _require(root.is_absolute(), "materialization output root must be absolute")
    try:
        authorized_parent = AUTHORIZED_OUTPUT_PARENT.resolve(strict=True)
        parent = root.parent.resolve(strict=True)
    except OSError as error:
        raise MaterializationError("materialization output parent is unavailable") from error
    _require(parent == authorized_parent and root.parent == authorized_parent, "output root must be a direct tst/.codex child")
    _require(not os.path.lexists(root), "materialization output root already exists")
    return root


def _write_new_file(path: Path, payload: bytes) -> None:
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags, 0o600)
    try:
        written = 0
        while written < len(payload):
            written += os.write(descriptor, payload[written:])
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _freeze_tree(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        mode = path.stat(follow_symlinks=False).st_mode
        if stat.S_ISREG(mode):
            path.chmod(stat.S_IMODE(mode) & ~_WRITE_BITS)
        elif stat.S_ISDIR(mode):
            path.chmod(0o555)
        else:
            raise MaterializationError(f"unexpected materialized member type: {path}")
    root.chmod(0o555)


def materialize_pilot_artifacts(
    output_root: Path,
    *,
    resource_plan: Path = RESOURCE_PLAN,
    selected_pressure_receipt: Path = SELECTED_PRESSURE_RECEIPT,
    production_deck: Path = PRODUCTION_DECK,
) -> dict[str, object]:
    """Write one deterministic, recursively read-only review tree."""
    root = _safe_output_root(output_root)
    manifest, files = build_materialization(
        resource_plan=resource_plan,
        selected_pressure_receipt=selected_pressure_receipt,
        production_deck=production_deck,
    )
    files = dict(files)
    manifest_payload = _json_bytes(manifest)
    files["materialization_manifest.json"] = manifest_payload
    files["inventory.sha256"] = _inventory_payload(files)
    try:
        root.mkdir(mode=0o700)
        for relative, payload in sorted(files.items()):
            canonical = PurePosixPath(relative)
            _require(
                not canonical.is_absolute()
                and canonical.as_posix() == relative
                and all(part not in {"", ".", ".."} for part in canonical.parts),
                f"unsafe materialized member path: {relative}",
            )
            _write_new_file(root / relative, payload)
        _freeze_tree(root)
    except BaseException:
        if root.exists():
            for path in [root, *root.rglob("*")]:
                try:
                    path.chmod(path.stat().st_mode | stat.S_IWUSR)
                except OSError:
                    pass
            shutil.rmtree(root, ignore_errors=True)
        raise
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--print-manifest", action="store_true")
    args = parser.parse_args()
    if args.print_manifest:
        manifest, _ = build_materialization()
        print(_json_bytes(manifest).decode("utf-8"), end="")
        return
    manifest = materialize_pilot_artifacts(args.output_root)
    print(_json_bytes(manifest).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
