#!/usr/bin/env python3
"""Focused tests for the six-run Q019 carrier calibration lifecycle."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from tst.publication import (
    q019_q023_carrier_calibration_campaign_driver_v1 as driver,
)
from tst.publication import (
    q019_q023_carrier_calibration_engineering_qualification_v1 as qualification,
)
from tst.publication import (
    q019_q023_carrier_calibration_launch_policy_preparation_v1 as preparation,
)
from tst.publication import (
    q019_q023_carrier_calibration_retirement_v1 as retirement,
)


def _provenance(index: int, event_sha256: str) -> dict[str, object]:
    digest = lambda offset: f"{(index + offset) % 16:x}" * 64
    return {
        "source_commit": f"{index % 16:x}" * 40,
        "source_bundle_sha256": digest(1),
        "source_archive_sha256": digest(2),
        "executable_sha256": digest(3),
        "environment_sha256": digest(4),
        "deck_sha256": digest(5),
        "pre_submit_manifest": {"sha256": digest(6)},
        "clean_candidate_manifest": {"sha256": digest(7)},
        "source_archive": {"sha256": digest(2)},
        "executable_snapshot": {"sha256": digest(3)},
        "execution_receipt": {"sha256": digest(8)},
        "terminal_receipt": {"sha256": digest(9)},
        "installed_producer": {"entrypoint_sha256": digest(10)},
        "reconciliation_event_sha256": event_sha256,
        "reconciliation_mirror_ack_sha256": digest(11),
        "raw_inventory_sha256": digest(12),
        "retained_raw_bindings": [],
        "reduction_binding": {"record_type": "fixture"},
    }


def _memory_rows(job_id: str, index: int) -> list[dict[str, object]]:
    return [
        {
            "job_id_raw": job_id,
            "job_name": "frontier_job.sh",
            "state": "COMPLETED",
            "exit_code": "0:0",
            "elapsed_seconds": 200,
            "allocated_nodes": 1,
            "task_count": 0,
            "max_rss": "",
            "max_rss_node": "",
            "max_rss_task": "",
            "tres_usage_in_max": "",
            "tres_usage_in_max_node": "",
            "tres_usage_in_max_task": "",
        },
        {
            "job_id_raw": f"{job_id}.0",
            "job_name": "athena",
            "state": "COMPLETED",
            "exit_code": "0:0",
            "elapsed_seconds": 190,
            "allocated_nodes": 1,
            "task_count": 8 if index <= 3 else 16,
            "max_rss": f"{100000 + index}K",
            "max_rss_node": f"frontier{index:05d}",
            "max_rss_task": "0",
            "tres_usage_in_max": f"mem={100000 + index}K",
            "tres_usage_in_max_node": f"mem=frontier{index:05d}",
            "tres_usage_in_max_task": "mem=0",
        },
    ]


def _attempt(index: int, artifact_id: str) -> dict[str, object]:
    case_id = driver.EXPECTED_CASES[artifact_id]
    cycle = driver.EXPECTED_CYCLES[artifact_id]
    nodes = driver.EXPECTED_NODES[artifact_id]
    tasks = driver.EXPECTED_TASKS[artifact_id]
    elapsed_by_cycle = {32: 40, 256: 200, 16: 80, 64: 240}
    elapsed = elapsed_by_cycle[cycle]
    if artifact_id.endswith("-instrumented"):
        elapsed += 20
    submission = f"0000000{index}-0000-4000-8000-00000000000{index}"
    job_id = str(6000000 + index)
    authorization = driver.EXPECTED_AUTHORIZATIONS[artifact_id]
    artifact_root = str(
        driver.PIC_ROOT / driver.preparation.RUN_NAMESPACE / submission
    )
    event = {
        "sequence_number": index,
        "previous_event_sha256": f"{index + 4:02x}" * 32,
        "timestamp": f"2026-06-08T12:0{index}:00Z",
        "event_type": "reconciliation",
        "campaign": driver.preparation.CAMPAIGN,
        "test_id": case_id,
        "registered_science_authorization_id": authorization,
        "submission_id": submission,
        "job_id": job_id,
        "artifact_dir": artifact_root,
        "state": "COMPLETED",
        "reconciled": True,
        "requested_nodes": nodes,
        "scheduler_reported_allocated_nodes": nodes,
        "scheduler_exit_code": "0:0",
        "elapsed_seconds": elapsed,
        "billed_nodes": nodes,
        "consumed_node_hours": nodes * elapsed / 3600.0,
    }
    event["event_sha256"] = driver.legacy._ledger_event_sha256(event)
    measurement = driver.build_resource_measurement(
        job_id=job_id, tasks=tasks, rows=_memory_rows(job_id, index)
    )
    provenance = _provenance(index, str(event["event_sha256"]))
    role = "instrumented" if artifact_id.endswith("-instrumented") else "baseline"
    dimension = 2 if "-2d-" in artifact_id else 3
    raw_bytes = dimension * 512 * 1024**2
    artifact_bytes = raw_bytes + 64 * 1024**2
    output_slots = 14 if cycle in {16, 32} else 42
    return {
        "attempt_index": index,
        "artifact_id": artifact_id,
        "source_case_id": case_id,
        "authorization_id": authorization,
        "submission_id": submission,
        "job_id": job_id,
        "artifact_root": artifact_root,
        "admission": {
            "path": f"/analysis/{artifact_id}.json",
            "sha256": f"{index}" * 64,
            "byte_count": 1000 + index,
        },
        "admission_facts": {
            "record_type": driver.registered_admission.RECORD_TYPE,
            "case_id": case_id,
            "artifact_root": artifact_root,
            "registered_execution_identity": {
                "submission_id": submission,
                "registered_science_authorization_id": authorization,
                "slurm_job_id": job_id,
                "reconciliation_event_sha256": event["event_sha256"],
                "reconciliation_mirror_ack_sha256": provenance[
                    "reconciliation_mirror_ack_sha256"
                ],
            },
            "execution_profile": {
                "kind": "runtime_controller_overlay",
                "artifact_id": artifact_id,
                "source_case_id": case_id,
                "authority": "excluded_pilot_only",
                "expected_stop_reason": 1903,
                "saturation_evidence_eligible": False,
            },
            "runtime_completion": {
                "run_completion_status": "completed_not_acceptance_eligible",
                "problem_stop_requested": True,
                "stop_reason_code": "1903",
                "runtime_controller_trigger_cycle": cycle,
                "process_exit_code": 0,
                "scheduler_terminal_state": "COMPLETED",
                "trusted_execution_binding_present": True,
            },
            "raw_science_admission_eligible": True,
            "problem_reported_saturation_evidence_eligible": False,
            "saturation_evidence_eligible": False,
            "authorization": {
                "scientific_claim_authorized": False,
                "publication_authorized": False,
            },
        },
        "execution_provenance": provenance,
        "artifact_inventory": {
            "path": f"{artifact_root}/artifact_inventory.json",
            "sha256": f"{index + 4}" * 64,
            "byte_count": 1024,
        },
        "artifact_payload_bytes": artifact_bytes,
        "raw_payload_bytes": raw_bytes,
        "output_slot_count": output_slots,
        "completed_cycles": cycle,
        "resource_measurement": measurement,
        "peak_resident_memory_bytes": measurement[
            "peak_resident_memory_bytes"
        ],
        "reconciliation_event": event,
        "fixture_role": role,
    }


def _attempts() -> list[dict[str, object]]:
    attempts = [
        _attempt(index, artifact_id)
        for index, artifact_id in enumerate(driver.EXPECTED_ARTIFACTS, 1)
    ]
    for item in attempts:
        item.pop("fixture_role")
    return attempts


def test_preparation_is_exact_six_run_non_authorizing_matrix() -> None:
    manifest, files = preparation.build_materialization()
    preparation.validate_materialization(manifest, files)
    assert manifest["pilot_attempt_count"] == 6
    assert len(files) == 14
    assert [item["cycle_limit"] for item in manifest["attempt_records"]] == [
        32,
        256,
        256,
        16,
        64,
        64,
    ]
    assert [item["pair_role"] for item in manifest["attempt_records"]] == [
        "baseline",
        "baseline",
        "instrumented",
        "baseline",
        "baseline",
        "instrumented",
    ]
    assert [
        json.loads(files[item["launch_candidate"]["path"]])["launch_contract"][
            "actions"
        ][0]["resources"]["nodes"]
        for item in manifest["attempt_records"]
    ] == [1, 1, 1, 2, 2, 2]
    assert all(value is False for value in preparation.AUTHORIZATION_BOUNDARY.values())


def test_memory_measurement_is_attributed_and_conservative() -> None:
    value = driver.build_resource_measurement(
        job_id="6000001", tasks=8, rows=_memory_rows("6000001", 1)
    )
    assert value["maximum_observed_task_rss_bytes"] == 100001 * 1024
    assert value["peak_resident_memory_bytes"] == 8 * 100001 * 1024
    assert value["peak_resident_memory_exact_aggregate"] is False
    assert value["trustworthy_for_capacity_planning"] is True
    missing = copy.deepcopy(_memory_rows("6000001", 1))
    missing[1]["max_rss"] = ""
    with pytest.raises(driver.DriverError, match="no attributable MaxRSS"):
        driver.build_resource_measurement(job_id="6000001", tasks=8, rows=missing)
    ambiguous = copy.deepcopy(_memory_rows("6000001", 1))
    ambiguous[1]["max_rss_task"] = ""
    with pytest.raises(driver.DriverError, match="node/task-count attribution"):
        driver.build_resource_measurement(
            job_id="6000001", tasks=8, rows=ambiguous
        )


def test_execution_index_captures_all_required_measurements() -> None:
    index = driver.build_execution_index(_attempts())
    driver.validate_execution_index(index)
    assert index["attempt_count"] == 6
    assert [item["dimension"] for item in index["startup_measurements"]] == [2, 3]
    assert len(index["steady_cycle_pairs"]) == 2
    assert index["steady_cycle_pairs"][0]["incremental_cycles"] == 224
    assert index["steady_cycle_pairs"][1]["incremental_cycles"] == 48
    assert index["total_raw_bytes"] > 0
    assert index["total_output_slots"] == 196
    assert index["production_resource_freeze_authorized"] is False
    assert all(value is False for value in index["authorization"].values())


def test_wrong_cycle_duplicate_identity_and_memory_claim_fail_closed() -> None:
    wrong_cycle = _attempts()
    wrong_cycle[0]["completed_cycles"] = 31
    with pytest.raises(driver.DriverError, match="required measurements drifted"):
        driver.build_execution_index(wrong_cycle)
    duplicate = _attempts()
    duplicate[1]["job_id"] = duplicate[0]["job_id"]
    duplicate[1]["admission_facts"]["registered_execution_identity"][
        "slurm_job_id"
    ] = duplicate[0]["job_id"]
    duplicate[1]["reconciliation_event"]["job_id"] = duplicate[0]["job_id"]
    duplicate[1]["reconciliation_event"]["event_sha256"] = (
        driver.legacy._ledger_event_sha256(duplicate[1]["reconciliation_event"])
    )
    duplicate[1]["execution_provenance"]["reconciliation_event_sha256"] = duplicate[
        1
    ]["reconciliation_event"]["event_sha256"]
    duplicate[1]["resource_measurement"] = driver.build_resource_measurement(
        job_id=duplicate[0]["job_id"],
        tasks=driver.EXPECTED_TASKS[duplicate[1]["artifact_id"]],
        rows=_memory_rows(duplicate[0]["job_id"], 2),
    )
    duplicate[1]["peak_resident_memory_bytes"] = duplicate[1][
        "resource_measurement"
    ]["peak_resident_memory_bytes"]
    with pytest.raises(driver.DriverError, match="identity reused: job_id"):
        driver.build_execution_index(duplicate)
    forged = _attempts()
    forged[2]["peak_resident_memory_bytes"] += 1
    with pytest.raises(driver.DriverError, match="required measurements drifted"):
        driver.build_execution_index(forged)


def _build_qualification(
    attempts: list[dict[str, object]],
) -> dict[str, object]:
    index = driver.build_execution_index(attempts)
    path = (
        driver.PIC_ROOT
        / "analysis"
        / driver.preparation.CAMPAIGN
        / "q019_q023_carrier_calibration_execution_index_v1.json"
    )
    with patch.object(driver, "validate_execution_index_files", return_value=index):
        return qualification.build_qualification(
            index,
            execution_index_path=path,
            execution_index_sha256="a" * 64,
        )


def test_qualification_builds_non_authorizing_resource_model() -> None:
    record = _build_qualification(_attempts())
    assert record["status"] == qualification.STATUS_PASS
    assert record["decision"]["engineering_gate_pass"] is True
    assert record["decision"]["production_resource_freeze_recommended"] is True
    assert record["decision"]["production_resource_freeze_authorized"] is False
    assert set(record["decision"]["measured_resource_model"]) == {"2", "3"}
    assert all(
        item["peak_resident_memory_is_conservative_upper_bound"]
        for item in record["dimension_qualifications"]
    )


def test_qualification_materialization_is_atomic_exact_and_idempotent(
    tmp_path: Path,
) -> None:
    index = driver.build_execution_index(_attempts())
    index_path = (
        tmp_path
        / "analysis"
        / driver.preparation.CAMPAIGN
        / "q019_q023_carrier_calibration_execution_index_v1.json"
    )
    index_path.parent.mkdir(parents=True)
    index_path.write_text(
        json.dumps(index, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    index_path.chmod(0o444)
    output_path = (
        index_path.parent
        / "q019_q023_carrier_calibration_engineering_qualification_v1.json"
    )
    with (
        patch.object(driver, "PIC_ROOT", tmp_path),
        patch.object(
            driver,
            "validate_execution_index_files",
            return_value=index,
        ),
    ):
        first = qualification.materialize_qualification(
            execution_index_path=index_path,
            output_path=output_path,
        )
        second = qualification.materialize_qualification(
            execution_index_path=index_path,
            output_path=output_path,
        )
    assert first == second
    assert output_path.stat().st_mode & 0o222 == 0
    assert json.loads(output_path.read_text(encoding="utf-8")) == first


def test_qualification_rejects_overhead_and_memory_without_gaining_authority() -> None:
    slow = _attempts()
    event = slow[2]["reconciliation_event"]
    event["elapsed_seconds"] = 280
    event["consumed_node_hours"] = (
        event["billed_nodes"] * event["elapsed_seconds"] / 3600.0
    )
    event["event_sha256"] = driver.legacy._ledger_event_sha256(event)
    slow[2]["execution_provenance"]["reconciliation_event_sha256"] = event[
        "event_sha256"
    ]
    failed = _build_qualification(slow)
    assert failed["status"] == qualification.STATUS_FAIL
    assert failed["decision"]["engineering_gate_pass"] is False
    assert failed["decision"]["production_resource_freeze_authorized"] is False

    memory = _attempts()
    short = memory[0]
    capacity = qualification.ALLOCATED_MEMORY_BYTES_PER_NODE
    rows = _memory_rows(short["job_id"], 1)
    rows[1]["max_rss"] = f"{capacity // driver.EXPECTED_TASKS[short['artifact_id']] // 1024}K"
    short["resource_measurement"] = driver.build_resource_measurement(
        job_id=short["job_id"],
        tasks=driver.EXPECTED_TASKS[short["artifact_id"]],
        rows=rows,
    )
    short["peak_resident_memory_bytes"] = short["resource_measurement"][
        "peak_resident_memory_bytes"
    ]
    record = _build_qualification(memory)
    assert record["status"] == qualification.STATUS_FAIL

    long_memory = _attempts()
    long_run = long_memory[1]
    rows = _memory_rows(long_run["job_id"], 2)
    rows[1]["max_rss"] = (
        f"{capacity // driver.EXPECTED_TASKS[long_run['artifact_id']] // 1024}K"
    )
    long_run["resource_measurement"] = driver.build_resource_measurement(
        job_id=long_run["job_id"],
        tasks=driver.EXPECTED_TASKS[long_run["artifact_id"]],
        rows=rows,
    )
    long_run["peak_resident_memory_bytes"] = long_run["resource_measurement"][
        "peak_resident_memory_bytes"
    ]
    record = _build_qualification(long_memory)
    assert record["status"] == qualification.STATUS_FAIL


def test_retirement_requires_exact_six_slice_policy(tmp_path: Path) -> None:
    manifest, _ = preparation.build_materialization()
    active = {
        "registered_science_slices": [
            {**item, "status": "authorized"} for item in manifest["policy_slices"]
        ],
        "olcf_side_storage": {"installed_control_plane_version": "a" * 64},
    }
    index = driver.build_execution_index(_attempts())
    index_path = tmp_path / "index.json"
    index_path.write_text("{}\n", encoding="utf-8")
    qualified = {
        "status": qualification.STATUS_PASS,
        "decision": {
            "engineering_gate_pass": True,
            "production_resource_freeze_authorized": False,
        },
        "execution_index": {"path": str(index_path)},
    }
    preflight = tmp_path / "preflight.json"
    preflight.write_text("{}\n", encoding="utf-8")
    successor = {**active, "registered_science_slices": []}
    with (
        patch.object(
            retirement.qualification,
            "validate_qualification",
            return_value=qualified,
        ),
        patch.object(
            retirement.preparation,
            "validate_final_binding_files",
            return_value={"record_type": "fixture"},
        ),
        patch.object(
            retirement.preparation,
            "build_materialization",
            return_value=(manifest, {}),
        ),
        patch.object(
            retirement.campaign,
            "validate_execution_index_files",
            return_value=index,
        ),
        patch.object(retirement, "validate_storage_policy"),
        patch.object(
            retirement.execution,
            "_advance_control_plane_fields",
            return_value=successor,
        ),
    ):
        assert retirement.materialize_retired_policy(
            active_policy=active,
            engineering_qualification=qualified,
            final_bindings={"record_type": "fixture"},
            successor_control_plane_version="b" * 64,
            storage_preflight_binding=preflight,
        )["registered_science_slices"] == []
        drifted = copy.deepcopy(active)
        drifted["registered_science_slices"][0]["test_id"] = "wrong"
        with pytest.raises(retirement.RetirementError, match="exact authorized"):
            retirement.materialize_retired_policy(
                active_policy=drifted,
                engineering_qualification=qualified,
                final_bindings={"record_type": "fixture"},
                successor_control_plane_version="b" * 64,
                storage_preflight_binding=preflight,
            )


def test_historical_four_run_packet_files_remain_at_head_bytes() -> None:
    root = Path(__file__).resolve().parents[2]
    paths = (
        "tst/publication/q019_excluded_pilot_launch_policy_preparation_v1.py",
        "tst/publication/q019_excluded_pilot_campaign_driver_v1.py",
        "tst/publication/q019_excluded_pilot_engineering_qualification_v1.py",
        "inputs/publication/q019_nonlinear_bell_runtime_controller_v1/deck_manifest.json",
    )
    for relative in paths:
        current = (root / relative).read_bytes()
        expected = __import__("subprocess").check_output(
            ["git", "show", f"HEAD:{relative}"], cwd=root
        )
        assert hashlib.sha256(current).digest() == hashlib.sha256(expected).digest()
