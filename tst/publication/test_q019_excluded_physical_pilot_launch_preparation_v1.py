#!/usr/bin/env python3
"""Tests for measured and cumulative-gated Q019 physical-pilot preparation."""

from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
import hashlib
import json
from unittest.mock import patch

import pytest

from tst.publication import q019_excluded_physical_pilot_launch_preparation_v1 as prep
from tst.publication import q019_excluded_physical_window_preregistration_v1 as prereg
from tst.publication import q019_excluded_pilot_engineering_qualification_v1 as engineering_qual
from tst.publication.test_q019_excluded_physical_window_qualification_v1 import (
    _reports as physical_reports,
)


def _engineering() -> tuple[dict[str, object], dict[str, object]]:
    execution = {
        "attempts": [
            {
                "artifact_id": "q019-controller-pilot-2d-instrumented",
                "source_case_id": "q019-fr-grid-k8-rho1em05-s0",
                "artifact_inventory": {"byte_count": 1024},
                "artifact_payload_bytes": 1024,
                "reconciliation_event": {
                    "elapsed_seconds": 40,
                    "scheduler_reported_allocated_nodes": 4,
                },
            },
            {
                "artifact_id": "q019-controller-pilot-3d-instrumented",
                "source_case_id": "q019-fr-3d-onset-small-s0",
                "artifact_inventory": {"byte_count": 2048},
                "artifact_payload_bytes": 2048,
                "reconciliation_event": {
                    "elapsed_seconds": 80,
                    "scheduler_reported_allocated_nodes": 16,
                },
            },
        ]
    }
    qualification = {
        "record_type": engineering_qual.RECORD_TYPE,
        "status": engineering_qual.STATUS_PASS,
    }
    return qualification, execution


def _resource_model(
    qualification: dict[str, object], execution: dict[str, object]
) -> dict[str, object]:
    members = prep._physical_members()
    rows = []
    for case_id in prep._all_case_ids():
        case = members[case_id]["case"]
        dimension = int(case["dimension"])
        source = (
            "q019-controller-pilot-2d-instrumented"
            if dimension == 2
            else "q019-controller-pilot-3d-instrumented"
        )
        elapsed, measured_nodes = (40, 4) if dimension == 2 else (80, 16)
        source_case_id = (
            "q019-fr-grid-k8-rho1em05-s0"
            if dimension == 2
            else "q019-fr-3d-onset-small-s0"
        )
        cases = {item["case_id"]: item for item in prep.design.expected_cases()}
        source_case = cases[source_case_id]
        source_work_units = prep._work_units(source_case)
        target_work_units = prep._work_units(case)
        target_cycles = prep._estimated_cycle_count(case)
        initial_cycles = prep._initial_cfl_limited_cycle_count(case)
        selected_nodes = prep.math.ceil(
            measured_nodes * target_work_units / source_work_units
        )
        measurement_work_units_per_node = source_work_units / measured_nodes
        target_work_units_per_node = target_work_units / selected_nodes
        target_meshblocks = prep._meshblock_count(case)
        measured_seconds_per_cycle = elapsed / engineering_qual.PILOT_CYCLE_COUNT
        projected_measured = (
            measured_seconds_per_cycle * target_cycles * target_work_units / source_work_units
        )
        nodes = selected_nodes
        projected_selected = projected_measured * measured_nodes / nodes
        safety_factor = 1.5
        required_athena = prep.math.ceil(safety_factor * projected_selected)
        target_output_slots = (
            10 * (prep.math.floor(float(case["terminal_time"]) / 0.1) + 1)
            + 2 * (prep.math.floor(float(case["terminal_time"]) / 0.5) + 1)
            + 2
        )
        measured_bytes = 1024 if dimension == 2 else 2048
        storage_scaling_equivalent = max(
            target_output_slots / 14,
            target_cycles / engineering_qual.PILOT_CYCLE_COUNT,
        )
        projected_storage = prep.math.ceil(
            measured_bytes
            * storage_scaling_equivalent
            * target_work_units
            / source_work_units
        )
        storage_safety_factor = 1.5
        required_storage = prep.math.ceil(storage_safety_factor * projected_storage)
        rows.append(
            {
                "case_id": case_id,
                "artifact_id": members[case_id]["overlay"]["artifact_id"],
                "stage": (
                    1
                    if case_id in prereg.PREDECESSOR_CORE
                    else 2
                    if case_id in prereg.WINDOW_CORE
                    else 3
                ),
                "dimension": dimension,
                "measurement_artifact_id": source,
                "measurement_elapsed_seconds": elapsed,
                "measurement_allocated_nodes": measured_nodes,
                "measurement_cycle_count": engineering_qual.PILOT_CYCLE_COUNT,
                "method": (
                    "linear_cycle_work_density_preserving_frontier_extrapolation_v2"
                ),
                "source_case_id": source_case_id,
                "source_work_units": source_work_units,
                "target_work_units": target_work_units,
                "measurement_work_units_per_node": measurement_work_units_per_node,
                "target_work_units_per_node": target_work_units_per_node,
                "target_meshblock_count": target_meshblocks,
                "target_initial_cfl_limited_cycle_count": initial_cycles,
                "nonlinear_timestep_design_B_over_B0": max(
                    1.0, float(case["absolute_total_energy_B_over_B0_bound"])
                ),
                "target_estimated_cycle_count": target_cycles,
                "measured_seconds_per_cycle": measured_seconds_per_cycle,
                "projected_seconds_at_measurement_nodes": projected_measured,
                "projected_seconds_at_selected_nodes": projected_selected,
                "safety_factor": safety_factor,
                "required_athena_walltime_seconds": required_athena,
                "nodes": nodes,
                "tasks": nodes * prep.TASKS_PER_NODE,
                "scheduler_walltime_seconds": required_athena + 300,
                "athena_walltime_seconds": required_athena,
                "measurement_artifact_payload_bytes": measured_bytes,
                "measurement_output_slot_equivalent": 14,
                "target_output_slot_equivalent": target_output_slots,
                "storage_scaling_equivalent": storage_scaling_equivalent,
                "projected_storage_bytes": projected_storage,
                "storage_safety_factor": storage_safety_factor,
                "required_storage_bytes": required_storage,
                "maximum_storage_bytes": required_storage,
            }
        )
    return {
        "schema_version": 1,
        "record_type": prep.RESOURCE_MODEL_RECORD_TYPE,
        "status": "measured_engineering_calibrated_resource_model",
        "engineering_qualification_canonical_sha256": prep._canonical_sha256(
            qualification
        ),
        "engineering_execution_index_canonical_sha256": prep._canonical_sha256(
            execution
        ),
        "rows": rows,
        "rows_sha256": prep._canonical_sha256(rows),
        "default_resources_permitted": False,
        "resource_freeze_authorized": False,
        "authorization": dict(prep.AUTHORIZATION),
    }


def _final_bindings() -> dict[str, object]:
    return {"installed_control_plane_version": "a" * 64}


def test_engineering_evidence_reopens_file_backed_execution_artifacts() -> None:
    execution = {"attempts": []}
    qualification = {
        "status": engineering_qual.STATUS_PASS,
        "decision": {
            "engineering_gate_pass": True,
            "production_resource_freeze_recommended": True,
            "production_resource_freeze_authorized": False,
        },
        "authorization": dict(engineering_qual.AUTHORIZATION_BOUNDARY),
        "execution_index": {"path": "/execution-index", "sha256": "a" * 64},
    }

    def read_binding(_binding: object, *, label: str) -> dict[str, object]:
        if label == "Q019 engineering qualification":
            return qualification
        if label == "Q019 engineering execution index":
            return execution
        raise AssertionError(label)

    with patch.object(
        prep, "_stable_read_only_json", side_effect=read_binding
    ), patch.object(
        prep.engineering_qual,
        "validate_qualification",
        return_value=qualification,
    ), patch.object(
        prep.engineering_campaign,
        "validate_execution_index_files",
        return_value=execution,
    ) as validate_files:
        assert prep._engineering_evidence({}) == (qualification, execution)

    validate_files.assert_called_once_with(execution)


def test_prior_case_evidence_is_rebuilt_from_registered_raw_bundle() -> None:
    case_id = prereg.PREDECESSOR_CORE[0]
    report = {"case_id": case_id, "status": "exact"}
    admission = {"case_id": case_id}
    bundle = {
        "snapshots": ["snapshot"],
        "particle_states": ["particle-state"],
        "completion_record": {"status": "complete"},
    }
    item = {
        "case_id": case_id,
        "analysis_report": {"path": "/report", "sha256": "a" * 64},
        "admission": {"path": "/admission", "sha256": "b" * 64},
    }

    def read_binding(_binding: object, *, label: str) -> dict[str, object]:
        if label.startswith("Q019 prior physical analysis report"):
            return report
        if label.startswith("Q019 prior physical admission"):
            return admission
        raise AssertionError(label)

    with patch.object(
        prep, "_stable_read_only_json", side_effect=read_binding
    ), patch.object(
        prep.registered_admission,
        "validate_analysis_bundle",
        return_value=(admission, bundle),
    ) as validate_bundle, patch.object(
        prep.physics_analysis,
        "analyze_snapshots",
        return_value=report,
    ) as analyze:
        assert prep._reopen_prior_case_evidence(
            item, expected_case_id=case_id
        ) == report

    validate_bundle.assert_called_once_with(admission)
    analyze.assert_called_once_with(
        case_id,
        bundle["snapshots"],
        bundle["particle_states"],
        source_kind="raw_registered_bundle",
        provenance=admission,
        completion_record=bundle["completion_record"],
    )


def test_prior_case_evidence_rejects_report_raw_reanalysis_drift() -> None:
    case_id = prereg.PREDECESSOR_CORE[0]
    report = {"case_id": case_id, "status": "stored"}
    admission = {"case_id": case_id}
    bundle = {
        "snapshots": [],
        "particle_states": [],
        "completion_record": {},
    }
    item = {
        "case_id": case_id,
        "analysis_report": {"path": "/report", "sha256": "a" * 64},
        "admission": {"path": "/admission", "sha256": "b" * 64},
    }

    def read_binding(_binding: object, *, label: str) -> dict[str, object]:
        if label.startswith("Q019 prior physical analysis report"):
            return report
        if label.startswith("Q019 prior physical admission"):
            return admission
        raise AssertionError(label)

    with patch.object(
        prep, "_stable_read_only_json", side_effect=read_binding
    ), patch.object(
        prep.registered_admission,
        "validate_analysis_bundle",
        return_value=(admission, bundle),
    ), patch.object(
        prep.physics_analysis,
        "analyze_snapshots",
        return_value={"case_id": case_id, "status": "rebuilt"},
    ):
        with pytest.raises(
            prep.PreparationError,
            match="differs from exact raw reanalysis",
        ):
            prep._reopen_prior_case_evidence(item, expected_case_id=case_id)


def _budget(
    now: datetime,
    unreserved: float = 1000.0,
    *,
    prior_physical: float = 0.0,
    measured_minutes_ago: int = 1,
    expires_minutes_from_now: int = 10,
) -> dict[str, object]:
    policy_maximum = 10000.0
    consumed = 0.0
    reserved = policy_maximum - unreserved
    return {
        "schema_version": 1,
        "record_type": prep.BUDGET_SNAPSHOT_RECORD_TYPE,
        "source": "installed_control_plane_fresh_live_ledger_preflight",
        "measured_utc": (
            now - timedelta(minutes=measured_minutes_ago)
        ).isoformat().replace("+00:00", "Z"),
        "expires_utc": (
            now + timedelta(minutes=expires_minutes_from_now)
        ).isoformat().replace("+00:00", "Z"),
        "installed_control_plane_version": "a" * 64,
        "ledger": {"path": str(prep.LEDGER_PATH), "sha256": "1" * 64},
        "receipts": {"path": str(prep.RECEIPTS_PATH), "sha256": "2" * 64},
        "mirror": {"path": str(prep.MIRROR_PATH), "sha256": "3" * 64},
        "policy": {"path": str(prep.POLICY_PATH), "sha256": "4" * 64},
        "promotion": {"path": str(prep.PROMOTION_PATH), "sha256": "5" * 64},
        "policy_maximum_node_hours": policy_maximum,
        "cumulative_consumed_node_hours": consumed,
        "currently_reserved_node_hours": reserved,
        "then_unreserved_project_node_hours": unreserved,
        "prior_q019_physical_node_hours": prior_physical,
        "ledger_mutation_authorized": False,
        "authorization": dict(prep.AUTHORIZATION),
    }


def _derived_budget(snapshot: dict[str, object]) -> dict[str, float]:
    return {
        name: float(snapshot[name])
        for name in (
            "policy_maximum_node_hours",
            "cumulative_consumed_node_hours",
            "currently_reserved_node_hours",
            "then_unreserved_project_node_hours",
            "prior_q019_physical_node_hours",
        )
    }


def test_resource_model_is_exact_measured_and_non_authorizing() -> None:
    qualification, execution = _engineering()
    model = _resource_model(qualification, execution)
    assert (
        prep.validate_resource_model(
            model,
            engineering_qualification=qualification,
            engineering_execution_index=execution,
        )
        == model
    )
    model["rows"][0]["measurement_elapsed_seconds"] = 41
    with pytest.raises(prep.PreparationError, match="resource row"):
        prep.validate_resource_model(
            model,
            engineering_qualification=qualification,
            engineering_execution_index=execution,
        )


def test_stage_one_materialization_is_exact_serial_and_bounded() -> None:
    qualification, execution = _engineering()
    now = datetime(2026, 6, 8, 12, 0, tzinfo=timezone.utc)
    budget = _budget(now)
    with patch.object(prep, "_estimated_cycle_count", return_value=20), patch.object(
        prep, "_work_units", return_value=1
    ):
        model = _resource_model(qualification, execution)

        def read_binding(_binding: object, *, label: str) -> dict[str, object]:
            if label == "Q019 measured resource model":
                return model
            if label == "Q019 live unreserved budget snapshot":
                return budget
            raise AssertionError(label)

        arguments = {
            "stage": 1,
            "final_bindings": _final_bindings(),
            "engineering_qualification_binding": {"path": "/qual", "sha256": "a"},
            "resource_model_binding": {"path": "/model", "sha256": "b"},
            "budget_snapshot_binding": {"path": "/budget", "sha256": "c"},
            "now": now,
        }
        with patch.object(
            prep.engineering_prep,
            "validate_final_binding_files",
            return_value=arguments["final_bindings"],
        ), patch.object(
            prep, "_engineering_evidence", return_value=(qualification, execution)
        ), patch.object(
            prep, "_stable_read_only_json", side_effect=read_binding
        ), patch.object(
            prep, "_rederive_live_budget_snapshot", return_value=_derived_budget(budget)
        ):
            manifest, files = prep.build_materialization(**arguments)
            prep.validate_materialization(manifest, files, **arguments)

    assert manifest["stage_case_ids"] == list(prereg.PREDECESSOR_CORE)
    assert manifest["attempt_count"] == 7
    assert manifest["serial_submission_required"] is True
    assert all(value is False for value in manifest["authorization"].values())
    assert len(files) == 8
    launch_records = [
        json.loads(payload)
        for path, payload in files.items()
        if path.startswith("launch_candidates/")
    ]
    assert {
        record["registered_science_authorization_id"] for record in launch_records
    } == {
        prep._physical_authorization_id(case_id)
        for case_id in prereg.PREDECESSOR_CORE
    }


def test_later_stage_cannot_bypass_prior_gate() -> None:
    qualification, execution = _engineering()
    model = _resource_model(qualification, execution)
    now = datetime(2026, 6, 8, 12, 0, tzinfo=timezone.utc)

    def read_binding(_binding: object, *, label: str) -> dict[str, object]:
        if label == "Q019 measured resource model":
            return model
        if label == "Q019 live unreserved budget snapshot":
            return _budget(now)
        raise AssertionError(label)

    with patch.object(
        prep.engineering_prep,
        "validate_final_binding_files",
        return_value=_final_bindings(),
    ), patch.object(
        prep, "_engineering_evidence", return_value=(qualification, execution)
    ), patch.object(prep, "_stable_read_only_json", side_effect=read_binding):
        with pytest.raises(prep.PreparationError, match="prior physical qualification"):
            prep.build_materialization(
                stage=2,
                final_bindings=_final_bindings(),
                engineering_qualification_binding={},
                resource_model_binding={},
                budget_snapshot_binding={},
                prior_stage_qualification_binding=None,
                now=now,
            )


def test_forged_prior_stage_qualification_is_rederived_and_rejected() -> None:
    expected_ids = list(prereg.PREDECESSOR_CORE)
    reports = [{"case_id": case_id} for case_id in expected_ids]
    forged = {
        "record_type": prep.physical_qual.RECORD_TYPE,
        "status": prep.physical_qual.STATUS_PASS,
        "cumulative_stage": 1,
        "expected_case_ids": expected_ids,
        "decision": {
            "current_stage_gate_passed": True,
            "next_excluded_stage_submission_recommended": True,
            "next_excluded_stage_submission_authorized": False,
        },
        "authorization": dict(prep.physical_qual.AUTHORIZATION),
    }
    bundle = {
        "qualification": {"path": "/qualification", "sha256": "a"},
        "case_evidence": [
            {
                "case_id": case_id,
                "analysis_report": {
                    "path": f"/reports/{case_id}",
                    "sha256": "b",
                },
                "admission": {
                    "path": f"/admissions/{case_id}",
                    "sha256": "c",
                },
            }
            for case_id in expected_ids
        ],
    }

    def read_binding(binding: object, *, label: str) -> dict[str, object]:
        del binding
        if label == "Q019 prior physical qualification":
            return forged
        raise AssertionError(label)

    with patch.object(
        prep, "_stable_read_only_json", side_effect=read_binding
    ), patch.object(
        prep,
        "_reopen_prior_case_evidence",
        side_effect=lambda _item, *, expected_case_id: reports[
            expected_ids.index(expected_case_id)
        ],
    ):
        with pytest.raises(prep.PreparationError, match="did not rederive"):
            prep._prior_stage_evidence(2, bundle)


def test_live_fractional_budget_cap_is_enforced() -> None:
    qualification, execution = _engineering()
    model = _resource_model(qualification, execution)
    now = datetime(2026, 6, 8, 12, 0, tzinfo=timezone.utc)
    budget = _budget(now, unreserved=1.0)

    def read_binding(_binding: object, *, label: str) -> dict[str, object]:
        if label == "Q019 measured resource model":
            return model
        if label == "Q019 live unreserved budget snapshot":
            return budget
        raise AssertionError(label)

    with patch.object(
        prep.engineering_prep,
        "validate_final_binding_files",
        return_value=_final_bindings(),
    ), patch.object(
        prep, "_engineering_evidence", return_value=(qualification, execution)
    ), patch.object(
        prep, "_stable_read_only_json", side_effect=read_binding
    ), patch.object(
        prep, "_rederive_live_budget_snapshot", return_value=_derived_budget(budget)
    ):
        with pytest.raises(
            prep.PreparationError, match="cumulative 500 or live 10-percent"
        ):
            prep.build_materialization(
                stage=1,
                final_bindings=_final_bindings(),
                engineering_qualification_binding={},
                resource_model_binding={},
                budget_snapshot_binding={},
                now=now,
            )


def test_stale_long_lived_budget_snapshot_is_rejected() -> None:
    now = datetime(2026, 6, 8, 12, 0, tzinfo=timezone.utc)
    budget = _budget(
        now,
        measured_minutes_ago=16,
        expires_minutes_from_now=60,
    )
    with patch.object(prep, "_stable_read_only_json", return_value=budget), patch.object(
        prep, "_rederive_live_budget_snapshot"
    ) as rederive:
        with pytest.raises(prep.PreparationError, match="stale or malformed"):
            prep._budget_snapshot(
                {},
                now=now,
                final_bindings=_final_bindings(),
            )
    rederive.assert_not_called()


def test_cumulative_physical_budget_cap_is_enforced() -> None:
    qualification, execution = _engineering()
    now = datetime(2026, 6, 8, 12, 0, tzinfo=timezone.utc)
    budget = _budget(now, unreserved=10000.0, prior_physical=499.0)
    with patch.object(prep, "_estimated_cycle_count", return_value=20), patch.object(
        prep, "_work_units", return_value=1
    ):
        model = _resource_model(qualification, execution)

        def read_binding(_binding: object, *, label: str) -> dict[str, object]:
            if label == "Q019 measured resource model":
                return model
            if label == "Q019 live unreserved budget snapshot":
                return budget
            raise AssertionError(label)

        with patch.object(
            prep.engineering_prep,
            "validate_final_binding_files",
            return_value=_final_bindings(),
        ), patch.object(
            prep, "_engineering_evidence", return_value=(qualification, execution)
        ), patch.object(
            prep, "_stable_read_only_json", side_effect=read_binding
        ), patch.object(
            prep,
            "_rederive_live_budget_snapshot",
            return_value=_derived_budget(budget),
        ):
            with pytest.raises(
                prep.PreparationError, match="cumulative 500 or live 10-percent"
            ):
                prep.build_materialization(
                    stage=1,
                    final_bindings=_final_bindings(),
                    engineering_qualification_binding={},
                    resource_model_binding={},
                    budget_snapshot_binding={},
                    now=now,
                )


def test_cumulative_live_fractional_campaign_cap_is_enforced() -> None:
    qualification, execution = _engineering()
    now = datetime(2026, 6, 8, 12, 0, tzinfo=timezone.utc)
    budget = _budget(now, unreserved=100.0, prior_physical=9.0)
    with patch.object(prep, "_estimated_cycle_count", return_value=20), patch.object(
        prep, "_work_units", return_value=1
    ):
        model = _resource_model(qualification, execution)

        def read_binding(_binding: object, *, label: str) -> dict[str, object]:
            if label == "Q019 measured resource model":
                return model
            if label == "Q019 live unreserved budget snapshot":
                return budget
            raise AssertionError(label)

        with patch.object(
            prep.engineering_prep,
            "validate_final_binding_files",
            return_value=_final_bindings(),
        ), patch.object(
            prep, "_engineering_evidence", return_value=(qualification, execution)
        ), patch.object(
            prep, "_stable_read_only_json", side_effect=read_binding
        ), patch.object(
            prep,
            "_rederive_live_budget_snapshot",
            return_value=_derived_budget(budget),
        ):
            with pytest.raises(
                prep.PreparationError, match="cumulative 500 or live 10-percent"
            ):
                prep.build_materialization(
                    stage=1,
                    final_bindings=_final_bindings(),
                    engineering_qualification_binding={},
                    resource_model_binding={},
                    budget_snapshot_binding={},
                    now=now,
                )


@pytest.mark.parametrize("location", ["model", "row"])
def test_resource_model_rejects_unknown_authority_fields(location: str) -> None:
    qualification, execution = _engineering()
    model = _resource_model(qualification, execution)
    if location == "model":
        model["unregistered_authority"] = True
    else:
        model["rows"][0]["unregistered_authority"] = True
    with pytest.raises(prep.PreparationError, match="keys drifted"):
        prep.validate_resource_model(
            model,
            engineering_qualification=qualification,
            engineering_execution_index=execution,
        )


def test_boolean_stage_is_rejected() -> None:
    with pytest.raises(prep.PreparationError, match="stage is invalid"):
        prep._stage_case_ids(True)


def test_cycle_estimate_applies_cfl_and_nonlinear_safety_factor() -> None:
    case = {
        "case_id": "fixture",
        "cfl": 0.2,
        "initial_particle_cell_crossing_dt_bound": 0.5,
        "initial_particle_gyro_dt_bound": 1.0,
        "absolute_total_energy_B_over_B0_bound": 10.0,
        "terminal_time": 1.0,
    }
    assert prep._initial_cfl_limited_cycle_count(case) == 10
    assert prep._estimated_cycle_count(case) == 50


def test_live_budget_rederivation_binds_exact_policy_ledger_and_receipts() -> None:
    authorization_id = prep._physical_authorization_id(prep._all_case_ids()[0])
    records = [
        {
            "reservation_id": "physical",
            "registered_science_authorization_id": authorization_id,
            "state": "complete",
            "reconciled": True,
            "consumed_node_hours": 3.5,
        }
    ]
    receipts = [{"mirrored_event_sha256": "a" * 64}]
    payloads = {
        "ledger": (json.dumps(records[0]) + "\n").encode(),
        "receipts": (json.dumps(receipts[0]) + "\n").encode(),
        "mirror": (json.dumps(records[0]) + "\n").encode(),
        "policy": b"{}\n",
        "promotion": b"{}\n",
    }
    paths = {
        "ledger": prep.LEDGER_PATH,
        "receipts": prep.RECEIPTS_PATH,
        "mirror": prep.MIRROR_PATH,
        "policy": prep.POLICY_PATH,
        "promotion": prep.PROMOTION_PATH,
    }
    bindings = {
        name: {
            "path": str(paths[name]),
            "sha256": hashlib.sha256(payload).hexdigest(),
        }
        for name, payload in payloads.items()
    }
    snapshot = {
        **bindings,
        "installed_control_plane_version": "a" * 64,
    }

    class Ledger:
        @staticmethod
        @contextmanager
        def validated_read_only_mirrored_state_snapshot(*_args: object, **_kwargs: object):
            yield records

        @staticmethod
        def _read_jsonl_bytes(payload: bytes, *, path: object) -> list[dict[str, object]]:
            del path
            return [json.loads(line) for line in payload.decode().splitlines()]

        @staticmethod
        def _validate_receipt_records(
            values: list[dict[str, object]],
            _records: list[dict[str, object]],
            **_kwargs: object,
        ) -> list[dict[str, object]]:
            return values

        @staticmethod
        def require_explicit_genesis(_records: object) -> None:
            return None

        @staticmethod
        def validate_receipts(*_args: object, **_kwargs: object) -> list[dict[str, object]]:
            return receipts

        @staticmethod
        def accounting(_records: object) -> dict[str, float]:
            return {
                "cumulative_consumed_node_hours": 100.0,
                "currently_reserved_node_hours": 25.0,
            }

    policy_snapshot = {
        "active_policy_sha256": bindings["policy"]["sha256"],
        "active_promotion_sha256": bindings["promotion"]["sha256"],
    }

    class Common:
        @staticmethod
        def require_storage_policy_unlock_snapshot(
            **_kwargs: object,
        ) -> tuple[dict[str, object], dict[str, str]]:
            return {"frontier": {"maximum_node_hours": 10000.0}}, policy_snapshot

    @contextmanager
    def installed_modules(*_args: object, **_kwargs: object):
        yield {
            "ledger.py": Ledger,
            "control_plane_common.py": Common,
        }, {}

    def stable_bytes(path: object, **_kwargs: object) -> tuple[object, bytes]:
        name = next(name for name, expected in paths.items() if expected == path)
        return path, payloads[name]

    with patch.object(
        prep.engineering_prep.q043_prep,
        "_stable_regular_bytes",
        side_effect=stable_bytes,
    ), patch.object(
        prep.q043,
        "_installed_control_plane_modules",
        side_effect=installed_modules,
    ):
        derived = prep._rederive_live_budget_snapshot(
            snapshot,
            final_bindings=_final_bindings(),
        )
    assert derived == {
        "policy_maximum_node_hours": 10000.0,
        "cumulative_consumed_node_hours": 100.0,
        "currently_reserved_node_hours": 25.0,
        "then_unreserved_project_node_hours": 9875.0,
        "prior_q019_physical_node_hours": 3.5,
    }


def test_physical_node_hours_uses_latest_reservation_state() -> None:
    first, second, third = prep._all_case_ids()[:3]
    records = [
        {
            "reservation_id": "reserved",
            "registered_science_authorization_id": prep._physical_authorization_id(first),
            "state": "reserved",
            "reserved_node_hours": 10.0,
        },
        {
            "reservation_id": "reconciled",
            "registered_science_authorization_id": prep._physical_authorization_id(second),
            "state": "reserved",
            "reserved_node_hours": 20.0,
        },
        {
            "reservation_id": "reconciled",
            "registered_science_authorization_id": prep._physical_authorization_id(second),
            "state": "complete",
            "reconciled": True,
            "consumed_node_hours": 7.5,
        },
        {
            "reservation_id": "cancelled",
            "registered_science_authorization_id": prep._physical_authorization_id(third),
            "state": "cancelled",
            "reserved_node_hours": 30.0,
        },
        {
            "reservation_id": "other",
            "registered_science_authorization_id": "q023-linear",
            "state": "reserved",
            "reserved_node_hours": 40.0,
        },
    ]
    assert prep._physical_node_hours(records) == 17.5


@pytest.mark.parametrize("stage, expected_count", [(2, 8), (3, 4)])
def test_later_stage_materialization_accepts_rederived_gate(
    stage: int, expected_count: int
) -> None:
    qualification, execution = _engineering()
    now = datetime(2026, 6, 8, 12, 0, tzinfo=timezone.utc)
    budget = _budget(now, unreserved=10000.0)
    reports = physical_reports(stage - 1)
    prior = prep.physical_qual.build_qualification(reports, upto_stage=stage - 1)
    bundle = {
        "qualification": {"path": "/prior", "sha256": "a"},
        "case_evidence": [
            {
                "case_id": report["case_id"],
                "analysis_report": {
                    "path": f"/reports/{report['case_id']}",
                    "sha256": "b",
                },
                "admission": {
                    "path": f"/admissions/{report['case_id']}",
                    "sha256": "c",
                },
            }
            for report in reports
        ],
    }
    with patch.object(prep, "_estimated_cycle_count", return_value=20), patch.object(
        prep, "_work_units", return_value=1
    ):
        model = _resource_model(qualification, execution)

        def read_binding(_binding: object, *, label: str) -> dict[str, object]:
            if label == "Q019 measured resource model":
                return model
            if label == "Q019 live unreserved budget snapshot":
                return budget
            if label == "Q019 prior physical qualification":
                return prior
            raise AssertionError(label)

        with patch.object(
            prep.engineering_prep,
            "validate_final_binding_files",
            return_value=_final_bindings(),
        ), patch.object(
            prep, "_engineering_evidence", return_value=(qualification, execution)
        ), patch.object(
            prep, "_stable_read_only_json", side_effect=read_binding
        ), patch.object(
            prep,
            "_reopen_prior_case_evidence",
            side_effect=lambda _item, *, expected_case_id: next(
                report for report in reports if report["case_id"] == expected_case_id
            ),
        ), patch.object(
            prep,
            "_rederive_live_budget_snapshot",
            return_value=_derived_budget(budget),
        ):
            manifest, files = prep.build_materialization(
                stage=stage,
                final_bindings=_final_bindings(),
                engineering_qualification_binding={},
                resource_model_binding={},
                budget_snapshot_binding={},
                prior_stage_qualification_binding=bundle,
                now=now,
            )
    assert manifest["attempt_count"] == expected_count
    assert manifest["stage_case_ids"] == list(prep._stage_case_ids(stage))
    assert len([path for path in files if path.startswith("launch_candidates/")]) == expected_count
