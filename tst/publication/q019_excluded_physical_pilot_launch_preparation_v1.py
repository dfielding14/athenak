#!/usr/bin/env python3
"""Prepare one measured, cumulative-gated Q019 physical-pilot stage.

This module deliberately contains no default node, walltime, or storage
allocation.  A launch candidate can be built only from immutable engineering
qualification, measured resource-model, live budget, and prior-stage evidence.
The resulting packet remains non-authorizing until the installed control plane
performs the usual policy promotion and fresh submission checks.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import stat
from typing import Mapping, Sequence

from tst.publication import q019_excluded_physical_window_preregistration_v1 as prereg
from tst.publication import q019_excluded_physical_window_qualification_v1 as physical_qual
from tst.publication import q019_excluded_pilot_campaign_driver_v1 as engineering_campaign
from tst.publication import q019_excluded_pilot_engineering_qualification_v1 as engineering_qual
from tst.publication import q019_excluded_pilot_launch_policy_preparation_v1 as engineering_prep
from tst.publication import q019_nonlinear_bell_runtime_controller_v1 as controller
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as design
from tst.publication import (
    q043_registered_execution_raw_oracle_qualification_successor_v1 as q043,
)


SCHEMA_VERSION = 1
SUCCESSOR_ID = "q019_excluded_physical_pilot_launch_preparation_v1"
RESOURCE_MODEL_RECORD_TYPE = "q019_measured_physical_pilot_resource_model_v1"
BUDGET_SNAPSHOT_RECORD_TYPE = "q019_live_unreserved_budget_snapshot_v1"
MANIFEST_RECORD_TYPE = "q019_excluded_physical_pilot_stage_preparation_v1"
LAUNCH_RECORD_TYPE = "q019_excluded_physical_pilot_launch_candidate_v1"
MAXIMUM_LIVE_SUBMISSIONS = 1
TASKS_PER_NODE = engineering_prep.TASKS_PER_NODE
MINIMUM_TIMEOUT_MARGIN_SECONDS = 300
MAXIMUM_BUDGET_SNAPSHOT_AGE = timedelta(minutes=15)
PHYSICAL_AUTHORIZATION_PREFIX = "q019-physical-pilot-"
LEDGER_PATH = engineering_prep.AUTHORIZED_ORION_ROOT / "ledger/node_hours.jsonl"
RECEIPTS_PATH = (
    engineering_prep.AUTHORIZED_ORION_ROOT / "ledger/mirror_receipts.jsonl"
)
MIRROR_PATH = (
    Path("/autofs/nccs-svm1_proj/ast207/proj-shared/PIC")
    / "ledger/node_hours.jsonl"
)
POLICY_PATH = engineering_prep.AUTHORIZED_ORION_ROOT / "policy/storage_policy.json"
PROMOTION_PATH = (
    engineering_prep.AUTHORIZED_ORION_ROOT / "policy/active_promotion.json"
)
AUTHORIZATION = {
    "launch_authorized": False,
    "scheduler_submission_authorized": False,
    "policy_mutation_authorized": False,
    "frontier_execution_authorized": False,
    "production_resource_freeze_authorized": False,
    "production_deck_freeze_authorized": False,
    "q019_qualification_authorized": False,
    "nonlinear_saturation_claim_authorized": False,
    "scientific_claim_authorized": False,
    "publication_authorized": False,
}
RESOURCE_MODEL_KEYS = {
    "schema_version",
    "record_type",
    "status",
    "engineering_qualification_canonical_sha256",
    "engineering_execution_index_canonical_sha256",
    "rows",
    "rows_sha256",
    "default_resources_permitted",
    "resource_freeze_authorized",
    "authorization",
}
RESOURCE_ROW_KEYS = {
    "case_id",
    "artifact_id",
    "stage",
    "dimension",
    "measurement_artifact_id",
    "measurement_elapsed_seconds",
    "measurement_allocated_nodes",
    "measurement_cycle_count",
    "method",
    "source_case_id",
    "source_work_units",
    "target_work_units",
    "measurement_work_units_per_node",
    "target_work_units_per_node",
    "target_meshblock_count",
    "target_initial_cfl_limited_cycle_count",
    "nonlinear_timestep_design_B_over_B0",
    "target_estimated_cycle_count",
    "measured_seconds_per_cycle",
    "projected_seconds_at_measurement_nodes",
    "projected_seconds_at_selected_nodes",
    "safety_factor",
    "required_athena_walltime_seconds",
    "nodes",
    "tasks",
    "scheduler_walltime_seconds",
    "athena_walltime_seconds",
    "measurement_artifact_payload_bytes",
    "measurement_output_slot_equivalent",
    "target_output_slot_equivalent",
    "storage_scaling_equivalent",
    "projected_storage_bytes",
    "storage_safety_factor",
    "required_storage_bytes",
    "maximum_storage_bytes",
}
BUDGET_SNAPSHOT_KEYS = {
    "schema_version",
    "record_type",
    "source",
    "measured_utc",
    "expires_utc",
    "installed_control_plane_version",
    "ledger",
    "receipts",
    "mirror",
    "policy",
    "promotion",
    "policy_maximum_node_hours",
    "cumulative_consumed_node_hours",
    "currently_reserved_node_hours",
    "then_unreserved_project_node_hours",
    "prior_q019_physical_node_hours",
    "ledger_mutation_authorized",
    "authorization",
}


class PreparationError(ValueError):
    """Reject unmeasured resources, stale budgets, or bypassed stage gates."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PreparationError(message)


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def _canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


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


def _stable_read_only_json(binding: object, *, label: str) -> dict[str, object]:
    _require(
        type(binding) is dict and set(binding) == {"path", "sha256"},
        f"{label} binding is malformed",
    )
    path = Path(str(binding["path"]))
    digest = str(binding["sha256"])
    _require(path.is_absolute(), f"{label} path is not absolute")
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and not before.st_mode & 0o222,
            f"{label} is not one read-only regular file",
        )
        chunks = []
        while chunk := os.read(descriptor, 1024 * 1024):
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
        current = path.stat(follow_symlinks=False)
        _require(
            (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
            == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and len(payload) == after.st_size
            and _sha256(payload) == digest,
            f"{label} changed while reading or its digest drifted",
        )
    finally:
        os.close(descriptor)
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PreparationError(f"{label} is not UTF-8 JSON") from error
    _require(type(value) is dict, f"{label} is not a JSON object")
    return value


def _stage_case_ids(stage: int) -> tuple[str, ...]:
    _require(
        type(stage) is int and stage in {1, 2, 3},
        "Q019 physical pilot stage is invalid",
    )
    return (prereg.PREDECESSOR_CORE, prereg.WINDOW_CORE, prereg.WINDOW_REPLICATION)[
        stage - 1
    ]


def _all_case_ids() -> tuple[str, ...]:
    return prereg.PREDECESSOR_CORE + prereg.WINDOW_CORE + prereg.WINDOW_REPLICATION


def _root_cells(case: Mapping[str, object]) -> int:
    return math.prod(int(value) for value in case["nx"])


def _work_units(case: Mapping[str, object]) -> int:
    return _root_cells(case) * (1 + int(case["ppc"]))


def _meshblock_count(case: Mapping[str, object]) -> int:
    return math.prod(
        int(cells) // int(block)
        for cells, block in zip(case["nx"], case["meshblock_nx"])
    )


def _cycle_count(case: Mapping[str, object], *, nonlinear: bool) -> int:
    cfl = float(case["cfl"])
    cell_bound = float(case["initial_particle_cell_crossing_dt_bound"])
    gyro_bound = float(case["initial_particle_gyro_dt_bound"])
    magnetic_bound = (
        max(1.0, float(case["absolute_total_energy_B_over_B0_bound"]))
        if nonlinear
        else 1.0
    )
    timestep = cfl * min(cell_bound, gyro_bound / magnetic_bound)
    _require(
        math.isfinite(cfl)
        and 0.0 < cfl <= 1.0
        and math.isfinite(magnetic_bound)
        and magnetic_bound >= 1.0
        and math.isfinite(timestep)
        and timestep > 0.0,
        f"{case['case_id']}: estimated timestep is invalid",
    )
    return math.ceil(float(case["terminal_time"]) / timestep)


def _estimated_cycle_count(case: Mapping[str, object]) -> int:
    return _cycle_count(case, nonlinear=True)


def _initial_cfl_limited_cycle_count(case: Mapping[str, object]) -> int:
    return _cycle_count(case, nonlinear=False)


def _physical_authorization_id(case_id: str) -> str:
    digest = hashlib.sha256(case_id.encode("utf-8")).hexdigest()[:20]
    return f"{PHYSICAL_AUTHORIZATION_PREFIX}{digest}-v1"


def _physical_authorization_ids() -> frozenset[str]:
    return frozenset(_physical_authorization_id(case_id) for case_id in _all_case_ids())


def _physical_members() -> dict[str, dict[str, object]]:
    base = {str(item["case_id"]): item for item in design.validate_checked_in_decks()["decks"]}
    overlays = {
        str(item["source_case_id"]): item
        for item in controller.build_manifest()["artifacts"]
        if str(item["artifact_id"]).startswith("q019-physical-pilot-")
    }
    _require(
        tuple(overlays) == _all_case_ids() and set(overlays) == set(_all_case_ids()),
        "Q019 physical controller inventory drifted",
    )
    result = {}
    for case_id in _all_case_ids():
        overlay = overlays[case_id]
        path = controller.CHECKED_IN_ROOT / str(overlay["filename"])
        payload = path.read_bytes()
        _require(
            payload.decode() == controller.render_overlay(overlay)
            and _sha256(payload) == overlay["rendered_sha256"],
            f"{case_id}: Q019 physical controller deck drifted",
        )
        result[case_id] = {
            "case": dict(base[case_id]),
            "overlay": dict(overlay),
            "deck_path": str(path.relative_to(engineering_prep.REPO_ROOT)),
            "deck_sha256": _sha256(payload),
            "deck_byte_count": len(payload),
        }
    return result


def _engineering_evidence(binding: object) -> tuple[dict[str, object], dict[str, object]]:
    qualification = _stable_read_only_json(binding, label="Q019 engineering qualification")
    try:
        qualification = engineering_qual.validate_qualification(qualification)
    except engineering_qual.QualificationError as error:
        raise PreparationError("Q019 engineering qualification did not validate") from error
    _require(
        qualification["status"] == engineering_qual.STATUS_PASS
        and qualification["decision"]["engineering_gate_pass"] is True
        and qualification["decision"]["production_resource_freeze_recommended"] is True
        and qualification["decision"]["production_resource_freeze_authorized"] is False
        and all(value is False for value in qualification["authorization"].values()),
        "Q019 engineering gate is not passing and non-authorizing",
    )
    index_binding = qualification["execution_index"]
    execution = _stable_read_only_json(
        {"path": index_binding["path"], "sha256": index_binding["sha256"]},
        label="Q019 engineering execution index",
    )
    try:
        execution = engineering_campaign.validate_execution_index(execution)
    except engineering_campaign.DriverError as error:
        raise PreparationError("Q019 engineering execution index did not validate") from error
    return qualification, execution


def validate_resource_model(
    value: object,
    *,
    engineering_qualification: Mapping[str, object],
    engineering_execution_index: Mapping[str, object],
) -> dict[str, object]:
    _require(
        type(value) is dict and set(value) == RESOURCE_MODEL_KEYS,
        "Q019 measured resource model keys drifted",
    )
    rows = value.get("rows")
    _require(
        value.get("schema_version") == SCHEMA_VERSION
        and value.get("record_type") == RESOURCE_MODEL_RECORD_TYPE
        and value.get("status") == "measured_engineering_calibrated_resource_model"
        and type(rows) is list
        and len(rows) == 19,
        "Q019 measured resource model identity or row count drifted",
    )
    _require(
        value.get("engineering_qualification_canonical_sha256")
        == _canonical_sha256(engineering_qualification)
        and value.get("engineering_execution_index_canonical_sha256")
        == _canonical_sha256(engineering_execution_index),
        "Q019 measured resource model engineering binding drifted",
    )
    attempts = {
        str(item["artifact_id"]): item for item in engineering_execution_index["attempts"]
    }
    members = _physical_members()
    expected_sources = {
        2: "q019-controller-pilot-2d-instrumented",
        3: "q019-controller-pilot-3d-instrumented",
    }
    design_cases = {str(item["case_id"]): item for item in design.expected_cases()}
    for case_id, row in zip(_all_case_ids(), rows):
        _require(
            type(row) is dict and set(row) == RESOURCE_ROW_KEYS,
            f"{case_id}: resource row keys drifted",
        )
        case = members[case_id]["case"]
        source_id = expected_sources[int(case["dimension"])]
        source = attempts[source_id]
        source_case = design_cases[str(source["source_case_id"])]
        nodes = row.get("nodes")
        tasks = row.get("tasks")
        scheduler = row.get("scheduler_walltime_seconds")
        athena = row.get("athena_walltime_seconds")
        storage = row.get("maximum_storage_bytes")
        measured_seconds_per_cycle = (
            float(source["reconciliation_event"]["elapsed_seconds"])
            / engineering_qual.PILOT_CYCLE_COUNT
        )
        target_cycles = _estimated_cycle_count(case)
        initial_cycles = _initial_cfl_limited_cycle_count(case)
        source_work_units = _work_units(source_case)
        target_work_units = _work_units(case)
        measured_nodes = int(
            source["reconciliation_event"]["scheduler_reported_allocated_nodes"]
        )
        selected_nodes = math.ceil(
            measured_nodes * target_work_units / source_work_units
        )
        measurement_work_units_per_node = source_work_units / measured_nodes
        target_work_units_per_node = target_work_units / selected_nodes
        target_meshblocks = _meshblock_count(case)
        safety_factor = float(row.get("safety_factor", math.nan))
        projected_at_measured_nodes = (
            measured_seconds_per_cycle
            * target_cycles
            * target_work_units
            / source_work_units
        )
        projected_at_selected_nodes = (
            projected_at_measured_nodes
            * measured_nodes
            / int(nodes)
            if type(nodes) is int and nodes > 0
            else math.inf
        )
        required_athena_seconds = math.ceil(
            safety_factor * projected_at_selected_nodes
        ) if math.isfinite(safety_factor * projected_at_selected_nodes) else -1
        measured_bytes = int(source["artifact_payload_bytes"])
        target_output_slots = (
            10 * (math.floor(float(case["terminal_time"]) / 0.1) + 1)
            + 2 * (math.floor(float(case["terminal_time"]) / 0.5) + 1)
            + 2
        )
        measurement_output_slots = 14
        storage_scaling_equivalent = max(
            target_output_slots / measurement_output_slots,
            target_cycles / engineering_qual.PILOT_CYCLE_COUNT,
        )
        projected_storage = math.ceil(
            measured_bytes
            * storage_scaling_equivalent
            * target_work_units
            / source_work_units
        )
        storage_safety_factor = float(row.get("storage_safety_factor", math.nan))
        required_storage = math.ceil(storage_safety_factor * projected_storage)
        if not math.isfinite(storage_safety_factor * projected_storage):
            required_storage = -1
        _require(
            row.get("case_id") == case_id
            and row.get("artifact_id") == members[case_id]["overlay"]["artifact_id"]
            and row.get("stage")
            == (1 if case_id in prereg.PREDECESSOR_CORE else 2 if case_id in prereg.WINDOW_CORE else 3)
            and row.get("dimension") == case["dimension"]
            and row.get("measurement_artifact_id") == source_id
            and row.get("measurement_elapsed_seconds")
            == source["reconciliation_event"]["elapsed_seconds"]
            and row.get("measurement_allocated_nodes")
            == source["reconciliation_event"]["scheduler_reported_allocated_nodes"]
            and row.get("measurement_cycle_count") == engineering_qual.PILOT_CYCLE_COUNT
            and row.get("method")
            == "linear_cycle_work_density_preserving_frontier_extrapolation_v2"
            and row.get("source_case_id") == source_case["case_id"]
            and row.get("source_work_units") == source_work_units
            and row.get("target_work_units") == target_work_units
            and math.isclose(
                float(row.get("measurement_work_units_per_node", math.nan)),
                measurement_work_units_per_node,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            and math.isclose(
                float(row.get("target_work_units_per_node", math.nan)),
                target_work_units_per_node,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            and target_work_units_per_node <= measurement_work_units_per_node
            and row.get("target_meshblock_count") == target_meshblocks
            and row.get("target_initial_cfl_limited_cycle_count") == initial_cycles
            and math.isclose(
                float(row.get("nonlinear_timestep_design_B_over_B0", math.nan)),
                max(1.0, float(case["absolute_total_energy_B_over_B0_bound"])),
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            and row.get("target_estimated_cycle_count") == target_cycles
            and math.isclose(
                float(row.get("measured_seconds_per_cycle", math.nan)),
                measured_seconds_per_cycle,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            and math.isclose(
                float(row.get("projected_seconds_at_measurement_nodes", math.nan)),
                projected_at_measured_nodes,
                rel_tol=1.0e-12,
                abs_tol=1.0e-9,
            )
            and math.isclose(
                float(row.get("projected_seconds_at_selected_nodes", math.nan)),
                projected_at_selected_nodes,
                rel_tol=1.0e-12,
                abs_tol=1.0e-9,
            )
            and type(row.get("safety_factor")) in {int, float}
            and math.isfinite(safety_factor)
            and safety_factor >= 1.5
            and row.get("required_athena_walltime_seconds")
            == required_athena_seconds
            and type(nodes) is int
            and nodes == selected_nodes
            and type(tasks) is int
            and tasks == nodes * TASKS_PER_NODE
            and tasks <= target_meshblocks
            and type(scheduler) is int
            and type(athena) is int
            and scheduler >= athena + MINIMUM_TIMEOUT_MARGIN_SECONDS
            and athena == required_athena_seconds
            and row.get("measurement_artifact_payload_bytes") == measured_bytes
            and row.get("measurement_output_slot_equivalent")
            == measurement_output_slots
            and row.get("target_output_slot_equivalent") == target_output_slots
            and math.isclose(
                float(row.get("storage_scaling_equivalent", math.nan)),
                storage_scaling_equivalent,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            and row.get("projected_storage_bytes") == projected_storage
            and type(row.get("storage_safety_factor")) in {int, float}
            and math.isfinite(storage_safety_factor)
            and storage_safety_factor >= 1.5
            and row.get("required_storage_bytes") == required_storage
            and type(storage) is int
            and storage == required_storage,
            f"{case_id}: resource row is not exact, measured, or internally consistent",
        )
    _require(
        value.get("rows_sha256") == _canonical_sha256(rows)
        and value.get("default_resources_permitted") is False
        and value.get("resource_freeze_authorized") is False
        and _strict_equal(value.get("authorization"), AUTHORIZATION),
        "Q019 measured resource model acquired authority or drifted",
    )
    return dict(value)


def _prior_stage_evidence(stage: int, binding: object | None) -> dict[str, object] | None:
    if stage == 1:
        _require(binding is None, "Q019 stage 1 must not carry a prior-stage qualification")
        return None
    _require(binding is not None, "Q019 prior physical qualification is required")
    _require(
        type(binding) is dict
        and set(binding) == {"qualification", "analysis_reports"}
        and type(binding["analysis_reports"]) is list,
        "Q019 prior physical qualification bundle is malformed",
    )
    value = _stable_read_only_json(
        binding["qualification"], label="Q019 prior physical qualification"
    )
    expected_ids = (
        list(prereg.PREDECESSOR_CORE)
        if stage == 2
        else list(prereg.PREDECESSOR_CORE + prereg.WINDOW_CORE)
    )
    _require(
        len(binding["analysis_reports"]) == len(expected_ids),
        "Q019 prior physical analysis-report count drifted",
    )
    reports = [
        _stable_read_only_json(
            report_binding,
            label=f"Q019 prior physical analysis report {case_id}",
        )
        for case_id, report_binding in zip(expected_ids, binding["analysis_reports"])
    ]
    _require(
        [report.get("case_id") for report in reports] == expected_ids,
        "Q019 prior physical analysis-report order or identity drifted",
    )
    try:
        value = physical_qual.validate_qualification(value, reports)
    except physical_qual.QualificationError as error:
        raise PreparationError(
            "Q019 prior physical qualification did not rederive from registered reports"
        ) from error
    _require(
        value.get("record_type") == physical_qual.RECORD_TYPE
        and value.get("status") == physical_qual.STATUS_PASS
        and value.get("cumulative_stage") == stage - 1
        and value.get("expected_case_ids") == expected_ids
        and value.get("decision", {}).get("current_stage_gate_passed") is True
        and value.get("decision", {}).get("next_excluded_stage_submission_recommended")
        is True
        and value.get("decision", {}).get("next_excluded_stage_submission_authorized")
        is False
        and _strict_equal(value.get("authorization"), physical_qual.AUTHORIZATION),
        "Q019 prior physical stage is not a passing non-authorizing gate",
    )
    return value


def _physical_node_hours(records: Sequence[Mapping[str, object]]) -> float:
    latest = {}
    for record in records:
        reservation_id = str(record.get("reservation_id", ""))
        if reservation_id:
            latest[reservation_id] = record
    total = 0.0
    for record in latest.values():
        authorization_id = str(
            record.get("registered_science_authorization_id", "")
        )
        if authorization_id not in _physical_authorization_ids():
            continue
        if record.get("reconciled") is True:
            total += float(record["consumed_node_hours"])
        elif record.get("state") not in {"cancelled", "submission_attach_failed"}:
            total += float(record["reserved_node_hours"])
    return total


def _rederive_live_budget_snapshot(
    value: Mapping[str, object],
    *,
    final_bindings: Mapping[str, object],
) -> dict[str, float]:
    expected_paths = {
        "ledger": LEDGER_PATH,
        "receipts": RECEIPTS_PATH,
        "mirror": MIRROR_PATH,
        "policy": POLICY_PATH,
        "promotion": PROMOTION_PATH,
    }
    bound_payloads = {}
    bound_digests = {}
    for name, expected_path in expected_paths.items():
        binding = value[name]
        _require(
            type(binding) is dict
            and set(binding) == {"path", "sha256"}
            and binding["path"] == str(expected_path),
            f"Q019 live budget {name} binding drifted",
        )
        _, payload = engineering_prep.q043_prep._stable_regular_bytes(
            expected_path,
            label=f"Q019 live budget {name}",
            require_read_only=name in {"policy", "promotion"},
        )
        digest = _sha256(payload)
        _require(
            digest == binding["sha256"],
            f"Q019 live budget {name} digest drifted",
        )
        bound_payloads[name] = payload
        bound_digests[name] = digest
    version = str(final_bindings["installed_control_plane_version"])
    try:
        with q043._installed_control_plane_modules(
            version, ("ledger.py", "control_plane_common.py")
        ) as (modules, _pair):
            ledger = modules["ledger.py"]
            common = modules["control_plane_common.py"]
            policy, policy_snapshot = common.require_storage_policy_unlock_snapshot(
                control_plane_version=version,
                authorized_pic_root=engineering_prep.AUTHORIZED_ORION_ROOT,
                authorized_project_home_root=engineering_prep.CANONICAL_PROJECT_HOME_ROOT,
            )
            _require(
                policy_snapshot
                == {
                    "active_policy_sha256": bound_digests["policy"],
                    "active_promotion_sha256": bound_digests["promotion"],
                },
                "Q019 live budget active policy bindings drifted",
            )
            with ledger.validated_read_only_mirrored_state_snapshot(
                LEDGER_PATH,
                RECEIPTS_PATH,
                MIRROR_PATH,
                ledger_root=engineering_prep.AUTHORIZED_ORION_ROOT,
                receipts_root=engineering_prep.AUTHORIZED_ORION_ROOT,
                mirror_root=MIRROR_PATH.parents[1],
            ) as records:
                bound_receipts = ledger._read_jsonl_bytes(
                    bound_payloads["receipts"], path=RECEIPTS_PATH
                )
                _require(
                    bound_payloads["mirror"] == bound_payloads["ledger"]
                    and
                    records
                    == ledger._read_jsonl_bytes(
                        bound_payloads["ledger"], path=LEDGER_PATH
                    ),
                    "Q019 bound ledger or mirror bytes differ from validated snapshot",
                )
                ledger.require_explicit_genesis(records)
                validated_receipts = ledger.validate_receipts(
                    RECEIPTS_PATH,
                    records,
                    mirror_jsonl=MIRROR_PATH,
                    mirror_transport="filesystem_copy",
                    root=engineering_prep.AUTHORIZED_ORION_ROOT,
                )
                _require(
                    validated_receipts
                    == ledger._validate_receipt_records(
                        bound_receipts,
                        records,
                        path=RECEIPTS_PATH,
                        mirror_jsonl=MIRROR_PATH,
                        mirror_transport="filesystem_copy",
                    ),
                    "Q019 bound receipt bytes differ from validated snapshot",
                )
                totals = ledger.accounting(records)
                prior_physical = _physical_node_hours(records)
            _, final_policy_snapshot = common.require_storage_policy_unlock_snapshot(
                control_plane_version=version,
                authorized_pic_root=engineering_prep.AUTHORIZED_ORION_ROOT,
                authorized_project_home_root=engineering_prep.CANONICAL_PROJECT_HOME_ROOT,
            )
            _require(
                final_policy_snapshot == policy_snapshot,
                "Q019 active policy changed during budget validation",
            )
    except (OSError, ValueError, q043.AdmissionError) as error:
        raise PreparationError(
            "Q019 live mirrored budget state did not validate"
        ) from error
    for name in expected_paths:
        _, payload = engineering_prep.q043_prep._stable_regular_bytes(
            expected_paths[name],
            label=f"Q019 live budget post-validation {name}",
            require_read_only=name in {"policy", "promotion"},
        )
        _require(
            _sha256(payload) == bound_digests[name],
            f"Q019 live budget {name} changed during validation",
        )
    cap = float(policy["frontier"]["maximum_node_hours"])
    consumed = float(totals["cumulative_consumed_node_hours"])
    reserved = float(totals["currently_reserved_node_hours"])
    return {
        "policy_maximum_node_hours": cap,
        "cumulative_consumed_node_hours": consumed,
        "currently_reserved_node_hours": reserved,
        "then_unreserved_project_node_hours": cap - consumed - reserved,
        "prior_q019_physical_node_hours": prior_physical,
    }


def _budget_snapshot(
    binding: object,
    *,
    now: datetime,
    final_bindings: Mapping[str, object],
) -> dict[str, object]:
    value = _stable_read_only_json(binding, label="Q019 live unreserved budget snapshot")
    try:
        measured = datetime.fromisoformat(str(value["measured_utc"]).replace("Z", "+00:00"))
        expires = datetime.fromisoformat(str(value["expires_utc"]).replace("Z", "+00:00"))
    except (KeyError, ValueError) as error:
        raise PreparationError("Q019 live budget timestamps are malformed") from error
    unreserved = value.get("then_unreserved_project_node_hours")
    _require(
        set(value) == BUDGET_SNAPSHOT_KEYS
        and type(value.get("schema_version")) is int
        and value.get("schema_version") == SCHEMA_VERSION
        and value.get("record_type") == BUDGET_SNAPSHOT_RECORD_TYPE
        and value.get("source") == "installed_control_plane_fresh_live_ledger_preflight"
        and value.get("installed_control_plane_version")
        == final_bindings["installed_control_plane_version"]
        and measured.tzinfo is not None
        and expires.tzinfo is not None
        and measured <= now <= expires
        and now - measured <= MAXIMUM_BUDGET_SNAPSHOT_AGE
        and expires - measured <= MAXIMUM_BUDGET_SNAPSHOT_AGE
        and type(unreserved) in {int, float}
        and type(unreserved) is not bool
        and math.isfinite(float(unreserved))
        and float(unreserved) > 0.0
        and value.get("ledger_mutation_authorized") is False
        and _strict_equal(value.get("authorization"), AUTHORIZATION),
        "Q019 live unreserved budget snapshot is stale or malformed",
    )
    derived = _rederive_live_budget_snapshot(value, final_bindings=final_bindings)
    _require(
        all(
            type(value[name]) in {int, float}
            and type(value[name]) is not bool
            and math.isclose(
                float(value[name]),
                expected,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            for name, expected in derived.items()
        ),
        "Q019 live budget arithmetic differs from mirrored ledger and policy",
    )
    return value


def _launch_contract(artifact_id: str, resources: Mapping[str, object]) -> dict[str, object]:
    action_id = f"q019-physical-{hashlib.sha256(artifact_id.encode()).hexdigest()[:16]}"
    contract = {
        "schema_version": 1,
        "executor": engineering_prep.TRUSTED_LAUNCH_EXECUTOR,
        "pre_actions": [],
        "actions": [
            {
                "action_id": action_id,
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
        return engineering_prep.validate_launch_contract(contract)
    except ValueError as error:
        raise PreparationError("Q019 physical launch contract failed validation") from error


def build_materialization(
    *,
    stage: int,
    final_bindings: Mapping[str, object],
    engineering_qualification_binding: object,
    resource_model_binding: object,
    budget_snapshot_binding: object,
    prior_stage_qualification_binding: object | None = None,
    now: datetime | None = None,
) -> tuple[dict[str, object], dict[str, bytes]]:
    final = engineering_prep.validate_final_binding_files(final_bindings)
    engineering, execution = _engineering_evidence(engineering_qualification_binding)
    resource_model = validate_resource_model(
        _stable_read_only_json(resource_model_binding, label="Q019 measured resource model"),
        engineering_qualification=engineering,
        engineering_execution_index=execution,
    )
    prior = _prior_stage_evidence(stage, prior_stage_qualification_binding)
    current_time = now or datetime.now(timezone.utc)
    budget_snapshot = _budget_snapshot(
        budget_snapshot_binding,
        now=current_time,
        final_bindings=final,
    )
    members = _physical_members()
    rows = {str(row["case_id"]): row for row in resource_model["rows"]}
    files: dict[str, bytes] = {}
    launches = []
    records = []
    for index, case_id in enumerate(_stage_case_ids(stage), 1):
        member = members[case_id]
        row = rows[case_id]
        artifact_id = str(member["overlay"]["artifact_id"])
        contract = _launch_contract(artifact_id, row)
        launch = {
            "schema_version": SCHEMA_VERSION,
            "record_type": LAUNCH_RECORD_TYPE,
            "successor_id": SUCCESSOR_ID,
            "status": "measured_stage_gated_review_candidate_not_authorized",
            "stage": stage,
            "attempt_index_within_stage": index,
            "artifact_id": artifact_id,
            "case_id": case_id,
            "registered_science_authorization_id": _physical_authorization_id(case_id),
            "checked_in_overlay_deck": {
                "path": member["deck_path"],
                "sha256": member["deck_sha256"],
                "byte_count": member["deck_byte_count"],
            },
            "controller_contract": {
                "physical_pilot_stage": member["overlay"]["physical_pilot_stage"],
                "controller_identity_fingerprint": member["overlay"][
                    "controller_parameters"
                ]["controller_identity_fingerprint"],
                "expected_stop_reason": None,
                "monitor_only_to_configured_time_horizon": True,
                "diagnostic_failure_stop_armed": True,
                "saturation_evidence_eligible": False,
            },
            "selected_final_bindings": final,
            "measured_resource_row": row,
            "launch_contract": contract,
            "launch_contract_sha256": engineering_prep.launch_contract_sha256(contract),
            "resource_ceiling": {
                "selected_qos": "normal",
                "maximum_nodes": row["nodes"],
                "maximum_walltime_seconds": row["scheduler_walltime_seconds"],
                "athena_walltime_seconds": row["athena_walltime_seconds"],
                "maximum_storage_bytes": row["maximum_storage_bytes"],
                "maximum_attempts": prereg.MAXIMUM_ATTEMPTS_PER_CASE,
                "maximum_retries": prereg.MAXIMUM_RETRIES_PER_CASE,
                "maximum_node_hours": row["nodes"]
                * row["scheduler_walltime_seconds"]
                / 3600.0,
            },
            "required_terminal_evidence": {
                "termination_reason": "Terminating on time limit",
                "runtime_controller_never_triggered": True,
                "problem_final_evidence_status": "completed_not_acceptance_eligible",
                "problem_saturation_evidence_eligible": "false",
            },
            "authorization": dict(AUTHORIZATION),
        }
        path = f"launch_candidates/stage_{stage}/{artifact_id}.json"
        payload = _json_bytes(launch)
        files[path] = payload
        launches.append(launch)
        records.append(
            {
                "case_id": case_id,
                "artifact_id": artifact_id,
                "launch_candidate": {
                    "path": path,
                    "sha256": _sha256(payload),
                    "byte_count": len(payload),
                },
            }
        )
    requested = sum(item["resource_ceiling"]["maximum_node_hours"] for item in launches)
    unreserved = float(budget_snapshot["then_unreserved_project_node_hours"])
    prior_physical = float(budget_snapshot["prior_q019_physical_node_hours"])
    available_fraction = prereg.MAXIMUM_UNRESERVED_BUDGET_FRACTION * unreserved
    _require(
        prior_physical + requested <= prereg.MAXIMUM_PILOT_NODE_HOURS
        and prior_physical + requested <= available_fraction,
        "Q019 physical stage exceeds cumulative 500 or live 10-percent node-hour cap",
    )
    budget = {
        "requested_stage_node_hours": requested,
        "then_unreserved_project_node_hours": unreserved,
        "prior_q019_physical_node_hours": prior_physical,
        "absolute_cap_node_hours": prereg.MAXIMUM_PILOT_NODE_HOURS,
        "unreserved_fraction_cap": prereg.MAXIMUM_UNRESERVED_BUDGET_FRACTION,
        "live_fractional_campaign_ceiling": available_fraction,
        "cumulative_campaign_node_hours_after_stage": prior_physical + requested,
        "within_cumulative_and_live_fractional_ceilings": True,
        "maximum_live_q019_submissions": MAXIMUM_LIVE_SUBMISSIONS,
        "maximum_attempts_per_case": prereg.MAXIMUM_ATTEMPTS_PER_CASE,
        "maximum_retries_per_case": prereg.MAXIMUM_RETRIES_PER_CASE,
        "ledger_mutation_authorized": False,
    }
    budget_payload = _json_bytes(budget)
    files[f"stage_{stage}_budget_accounting.json"] = budget_payload
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "record_type": MANIFEST_RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "status": "stage_preparation_complete_policy_and_launch_not_authorized",
        "stage": stage,
        "stage_case_ids": list(_stage_case_ids(stage)),
        "engineering_qualification": dict(engineering_qualification_binding),
        "engineering_qualification_canonical_sha256": _canonical_sha256(engineering),
        "engineering_execution_index_canonical_sha256": _canonical_sha256(execution),
        "resource_model": dict(resource_model_binding),
        "resource_model_canonical_sha256": _canonical_sha256(resource_model),
        "prior_stage_qualification": (
            None
            if prior is None
            else {
                "qualification": dict(
                    prior_stage_qualification_binding["qualification"]
                ),
                "analysis_reports": [
                    dict(item)
                    for item in prior_stage_qualification_binding[
                        "analysis_reports"
                    ]
                ],
                "canonical_sha256": _canonical_sha256(prior),
            }
        ),
        "budget_snapshot": dict(budget_snapshot_binding),
        "budget_accounting": {
            "path": f"stage_{stage}_budget_accounting.json",
            "sha256": _sha256(budget_payload),
            "byte_count": len(budget_payload),
        },
        "attempt_count": len(records),
        "attempt_records": records,
        "serial_submission_required": True,
        "stage_n_plus_one_requires_passing_stage_n": True,
        "all_attempts_and_failures_retained": True,
        "saturation_evidence_eligible": False,
        "authorization": dict(AUTHORIZATION),
    }
    return manifest, files


def validate_materialization(
    manifest: object,
    files: Mapping[str, bytes],
    **build_arguments: object,
) -> tuple[dict[str, object], dict[str, bytes]]:
    expected_manifest, expected_files = build_materialization(**build_arguments)
    _require(
        _strict_equal(manifest, expected_manifest),
        "Q019 physical stage manifest drifted",
    )
    _require(
        _strict_equal(files, expected_files),
        "Q019 physical stage file inventory drifted",
    )
    _require(
        all(value is False for value in expected_manifest["authorization"].values()),
        "Q019 physical stage preparation acquired authority",
    )
    launch_records = [
        json.loads(payload)
        for path, payload in expected_files.items()
        if path.startswith("launch_candidates/")
    ]
    _require(
        len(launch_records) == expected_manifest["attempt_count"]
        and all(
            _strict_equal(record.get("authorization"), AUTHORIZATION)
            and all(value is False for value in record["authorization"].values())
            for record in launch_records
        ),
        "Q019 physical launch candidate acquired authority",
    )
    return expected_manifest, expected_files


__all__ = [
    "AUTHORIZATION",
    "BUDGET_SNAPSHOT_RECORD_TYPE",
    "MANIFEST_RECORD_TYPE",
    "PreparationError",
    "RESOURCE_MODEL_RECORD_TYPE",
    "build_materialization",
    "validate_materialization",
    "validate_resource_model",
]
