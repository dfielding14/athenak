#!/usr/bin/env python3
"""Focused tests for the Q019 excluded-pilot campaign driver."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from tst.publication import q019_excluded_pilot_campaign_driver_v1 as driver


def _attempt(index: int, artifact_id: str) -> dict[str, object]:
    case_id = driver.EXPECTED_CASES[artifact_id]
    dimension = 2 if "-2d-" in artifact_id else 3
    role = "baseline" if artifact_id.endswith("-baseline") else "instrumented"
    elapsed = (100 if dimension == 2 else 300) + (10 if role == "instrumented" else 0)
    nodes = 4 if dimension == 2 else 16
    submission = f"0000000{index}-0000-4000-8000-00000000000{index}"
    authorization = driver.EXPECTED_AUTHORIZATIONS[artifact_id]
    artifact_root = (
        "/lustre/orion/ast207/proj-shared/dfielding/PIC/"
        f"runs/q019_nonlinear_bell_registered_successor_v1/{submission}"
    )
    event = {
        "sequence_number": index,
        "previous_event_sha256": f"{index + 4:02x}" * 32,
        "timestamp": f"2026-06-08T10:0{index}:00Z",
        "event_type": "reconciliation",
        "campaign": "q019_nonlinear_bell_registered_successor_v1",
        "test_id": case_id,
        "registered_science_authorization_id": authorization,
        "submission_id": submission,
        "job_id": str(5000000 + index),
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
    event["event_sha256"] = driver._ledger_event_sha256(event)
    return {
        "attempt_index": index,
        "artifact_id": artifact_id,
        "source_case_id": case_id,
        "authorization_id": authorization,
        "submission_id": submission,
        "job_id": str(5000000 + index),
        "artifact_root": artifact_root,
        "admission": {
            "path": f"/analysis/{artifact_id}.json",
            "sha256": f"{index}" * 64,
            "byte_count": 1000 + index,
        },
        "admission_facts": {
            "record_type": (
                "q019_hardened_installed_control_plane_registered_admission_v1"
            ),
            "case_id": case_id,
            "artifact_root": artifact_root,
            "registered_execution_identity": {
                "submission_id": submission,
                "registered_science_authorization_id": authorization,
                "slurm_job_id": str(5000000 + index),
                "reconciliation_event_sha256": event["event_sha256"],
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
                "runtime_controller_trigger_cycle": 20,
                "process_exit_code": 0,
                "scheduler_terminal_state": "COMPLETED",
                "trusted_execution_binding_present": True,
            },
            "raw_science_admission_eligible": True,
            "problem_reported_saturation_evidence_eligible": False,
            "saturation_evidence_eligible": False,
            "authorization": {"scientific_claim_authorized": False},
        },
        "artifact_inventory": {
            "path": f"{artifact_root}/artifact_inventory.json",
            "sha256": f"{index + 4}" * 64,
            "byte_count": 1024,
        },
        "artifact_payload_bytes": 1024**3 * dimension,
        "reconciliation_event": event,
    }


def _attempts() -> list[dict[str, object]]:
    return [
        _attempt(index, artifact_id)
        for index, artifact_id in enumerate(driver.EXPECTED_ARTIFACTS, 1)
    ]


def test_execution_index_is_exact_paired_and_non_authorizing() -> None:
    record = driver.build_execution_index(_attempts())
    driver.validate_execution_index(record)
    assert record["attempt_count"] == 4
    assert [item["dimension"] for item in record["pairs"]] == [2, 3]
    assert record["pairs"][0]["instrumentation_elapsed_ratio"] == 1.1
    assert record["pairs"][1]["instrumentation_elapsed_ratio"] == 310 / 300
    assert record["saturation_evidence_eligible"] is False
    assert record["production_resource_freeze_authorized"] is False
    assert all(value is False for value in record["authorization"].values())


def test_duplicate_source_cases_are_allowed_but_execution_ids_are_unique() -> None:
    attempts = _attempts()
    assert attempts[0]["source_case_id"] == attempts[1]["source_case_id"]
    assert attempts[2]["source_case_id"] == attempts[3]["source_case_id"]
    driver.build_execution_index(attempts)
    duplicate = copy.deepcopy(attempts)
    duplicate[1]["job_id"] = duplicate[0]["job_id"]
    duplicate[1]["admission_facts"]["registered_execution_identity"][
        "slurm_job_id"
    ] = duplicate[0]["job_id"]
    duplicate[1]["reconciliation_event"]["job_id"] = duplicate[0]["job_id"]
    duplicate[1]["reconciliation_event"]["event_sha256"] = (
        driver._ledger_event_sha256(duplicate[1]["reconciliation_event"])
    )
    duplicate[1]["admission_facts"]["registered_execution_identity"][
        "reconciliation_event_sha256"
    ] = duplicate[1]["reconciliation_event"]["event_sha256"]
    with pytest.raises(driver.DriverError, match="identity reused: job_id"):
        driver.build_execution_index(duplicate)


def test_artifact_order_case_binding_and_storage_ceiling_fail_closed() -> None:
    attempts = _attempts()
    swapped = [attempts[1], attempts[0], attempts[2], attempts[3]]
    with pytest.raises(driver.DriverError, match="execution identity drifted"):
        driver.build_execution_index(swapped)
    wrong_case = copy.deepcopy(attempts)
    wrong_case[0]["source_case_id"] = "q019-fr-3d-onset-small-s0"
    with pytest.raises(driver.DriverError, match="execution identity drifted"):
        driver.build_execution_index(wrong_case)
    oversized = copy.deepcopy(attempts)
    oversized[2]["artifact_payload_bytes"] = 257 * 1024**3
    with pytest.raises(driver.DriverError, match="storage binding drifted"):
        driver.build_execution_index(oversized)


def test_reconciliation_event_must_be_exact_and_successful() -> None:
    attempts = _attempts()
    failed = copy.deepcopy(attempts)
    failed[3]["reconciliation_event"]["state"] = "TIMEOUT"
    with pytest.raises(driver.DriverError, match="reconciliation event drifted"):
        driver.build_execution_index(failed)
    wrong_billing = copy.deepcopy(attempts)
    wrong_billing[0]["reconciliation_event"]["consumed_node_hours"] = 1.0
    with pytest.raises(driver.DriverError, match="reconciliation event drifted"):
        driver.build_execution_index(wrong_billing)


def test_exact_authorization_allocation_and_event_hash_are_required() -> None:
    attempts = _attempts()
    wrong_authorization = copy.deepcopy(attempts)
    wrong_authorization[0]["authorization_id"] = "q019-pilot-99-2d-baseline-v1"
    with pytest.raises(driver.DriverError, match="execution identity drifted"):
        driver.build_execution_index(wrong_authorization)
    wrong_nodes = copy.deepcopy(attempts)
    wrong_nodes[2]["reconciliation_event"]["billed_nodes"] = 15
    with pytest.raises(driver.DriverError, match="reconciliation event drifted"):
        driver.build_execution_index(wrong_nodes)
    tampered_event = copy.deepcopy(attempts)
    tampered_event[1]["reconciliation_event"]["timestamp"] = "2026-06-08T11:00:00Z"
    with pytest.raises(driver.DriverError, match="hash binding drifted"):
        driver.build_execution_index(tampered_event)


def test_admission_must_prove_excluded_controller_stop_contract() -> None:
    attempts = _attempts()
    wrong_cycle = copy.deepcopy(attempts)
    wrong_cycle[3]["admission_facts"]["runtime_completion"][
        "runtime_controller_trigger_cycle"
    ] = 19
    with pytest.raises(driver.DriverError, match="controller facts drifted"):
        driver.build_execution_index(wrong_cycle)
    saturation_authority = copy.deepcopy(attempts)
    saturation_authority[0]["admission_facts"]["saturation_evidence_eligible"] = True
    with pytest.raises(driver.DriverError, match="controller facts drifted"):
        driver.build_execution_index(saturation_authority)


def test_artifact_payload_bytes_sum_inventory_members(tmp_path: Path) -> None:
    inventory = {
        "schema_version": 1,
        "files": [
            {"path": "raw/a.bin", "sha256": "a" * 64, "size": 100},
            {"path": "athena_stdout.txt", "sha256": "b" * 64, "size": 23},
        ],
    }
    payload = (json.dumps(inventory, sort_keys=True) + "\n").encode()
    path = tmp_path / "artifact_inventory.json"
    path.write_bytes(payload)
    path.chmod(0o444)
    assert (
        driver._artifact_payload_bytes(
            {
                "path": str(path),
                "sha256": hashlib.sha256(payload).hexdigest(),
                "byte_count": len(payload),
            },
            artifact_root=tmp_path,
        )
        == 123
    )


def test_execution_index_file_validation_rederives_payload_bytes() -> None:
    attempts = _attempts()
    index = driver.build_execution_index(attempts)
    by_root = {
        attempt["artifact_root"]: attempt["artifact_payload_bytes"]
        for attempt in attempts
    }
    with patch.object(
        driver,
        "_artifact_payload_bytes",
        side_effect=lambda _binding, *, artifact_root: by_root[str(artifact_root)],
    ):
        assert driver.validate_execution_index_files(index) == index
    with patch.object(driver, "_artifact_payload_bytes", return_value=1):
        with pytest.raises(driver.DriverError, match="payload-byte total drifted"):
            driver.validate_execution_index_files(index)
