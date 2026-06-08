#!/opt/cray/pe/python/3.11.7/bin/python3
"""Run, measure, admit, and index six Q019 carrier calibration cases."""

from __future__ import annotations

import argparse
from datetime import timedelta
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import time
import uuid
from typing import Mapping, Sequence

from tst.publication import q019_excluded_pilot_campaign_driver_v1 as legacy
from tst.publication import (
    q019_hardened_installed_control_plane_registered_admission_v1
    as registered_admission,
)
from tst.publication import (
    q019_q023_carrier_calibration_launch_policy_preparation_v1 as preparation,
)
from tst.publication import (
    q019_q023_carrier_resource_calibration_runtime_controller_v1 as calibration,
)


PIC_ROOT = legacy.PIC_ROOT
PYTHON = legacy.PYTHON
LEDGER_JSONL = legacy.LEDGER_JSONL
LEDGER_CSV = legacy.LEDGER_CSV
RECEIPTS_JSONL = legacy.RECEIPTS_JSONL
MIRROR_JSONL = legacy.MIRROR_JSONL
JOB_ID_PATTERN = legacy.JOB_ID_PATTERN
SHA256_PATTERN = legacy.SHA256_PATTERN
COMMIT_PATTERN = legacy.COMMIT_PATTERN
UUID_PATTERN = legacy.UUID_PATTERN
SCHEMA_VERSION = 1
EXECUTION_INDEX_RECORD_TYPE = (
    "q019_q023_carrier_calibration_execution_index_v1"
)
RESOURCE_MEASUREMENT_RECORD_TYPE = (
    "q019_q023_carrier_calibration_resource_measurement_v1"
)
MEMORY_SEMANTICS = (
    "conservative_sum_of_slurm_step_max_task_rss_times_reported_step_tasks"
)
ACCOUNTING_CAPTURE_ATTEMPTS = 30
ACCOUNTING_CAPTURE_POLL_SECONDS = 10
AUTHORIZATION_BOUNDARY = {
    "launch_authorized": False,
    "scheduler_submission_authorized": False,
    "policy_mutation_authorized": False,
    "production_resource_freeze_authorized": False,
    "q019_qualification_authorized": False,
    "nonlinear_saturation_claim_authorized": False,
    "scientific_claim_authorized": False,
    "publication_authorized": False,
}

_OVERLAYS = tuple(calibration.expected_overlays())
EXPECTED_ARTIFACTS = tuple(str(item["artifact_id"]) for item in _OVERLAYS)
EXPECTED_CASES = {
    str(item["artifact_id"]): str(item["source_case_id"]) for item in _OVERLAYS
}
EXPECTED_CYCLES = {
    str(item["artifact_id"]): int(item["cycle_limit"]) for item in _OVERLAYS
}
EXPECTED_NODES = {
    str(item["artifact_id"]): int(item["nodes"]) for item in _OVERLAYS
}
EXPECTED_TASKS = {
    str(item["artifact_id"]): int(item["tasks"]) for item in _OVERLAYS
}
EXPECTED_AUTHORIZATIONS = {
    str(item["artifact_id"]): (
        f"q019-carrier-cal-{index:02d}-{int(item['dimension'])}d-"
        f"c{int(item['cycle_limit']):04d}-{item['instrumentation']}-v1"
    )
    for index, item in enumerate(_OVERLAYS, 1)
}
MAXIMUM_ELAPSED_SECONDS = {
    artifact_id: preparation.PILOT_RESOURCES[
        int(next(item["dimension"] for item in _OVERLAYS if item["artifact_id"] == artifact_id))
    ]["scheduler_walltime_seconds"]
    for artifact_id in EXPECTED_ARTIFACTS
}
MAXIMUM_STORAGE_BYTES = {
    artifact_id: preparation.PILOT_RESOURCES[
        int(next(item["dimension"] for item in _OVERLAYS if item["artifact_id"] == artifact_id))
    ]["maximum_storage_bytes"]
    for artifact_id in EXPECTED_ARTIFACTS
}


class DriverError(RuntimeError):
    """Reject drifted carrier calibration evidence or lifecycle state."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise DriverError(message)


def _canonical_sha256(value: object) -> str:
    payload = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _strict_equal(left: object, right: object) -> bool:
    return legacy._strict_equal(left, right)


def _memory_kib(value: object, *, label: str) -> int:
    _require(type(value) is str, f"{label}: Slurm memory value is not text")
    match = re.fullmatch(r"([1-9][0-9]*)K", value)
    _require(match is not None, f"{label}: Slurm memory value is absent or not KiB")
    return int(match.group(1))


def build_resource_measurement(
    *,
    job_id: str,
    tasks: int,
    rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Build a conservative job-memory bound from attributed Slurm step rows."""
    _require(
        type(job_id) is str
        and job_id.isdigit()
        and type(tasks) is int
        and tasks > 0
        and type(rows) in {list, tuple}
        and bool(rows),
        "carrier calibration accounting identity is malformed",
    )
    normalized = []
    observed = []
    for raw in rows:
        _require(type(raw) is dict, "carrier calibration sacct row is malformed")
        row = dict(raw)
        _require(
            set(row)
            == {
                "job_id_raw",
                "job_name",
                "state",
                "exit_code",
                "elapsed_seconds",
                "allocated_nodes",
                "task_count",
                "max_rss",
                "max_rss_node",
                "max_rss_task",
                "tres_usage_in_max",
                "tres_usage_in_max_node",
                "tres_usage_in_max_task",
            }
            and type(row["job_id_raw"]) is str
            and (
                row["job_id_raw"] == job_id
                or row["job_id_raw"].startswith(f"{job_id}.")
            )
            and type(row["job_name"]) is str
            and type(row["state"]) is str
            and bool(row["state"])
            and type(row["exit_code"]) is str
            and type(row["elapsed_seconds"]) is int
            and row["elapsed_seconds"] >= 0
            and type(row["allocated_nodes"]) is int
            and row["allocated_nodes"] >= 0
            and type(row["task_count"]) is int
            and row["task_count"] >= 0,
            "carrier calibration sacct row schema drifted",
        )
        if row["max_rss"]:
            kib = _memory_kib(
                row["max_rss"], label=f"{row['job_id_raw']}/MaxRSS"
            )
            _require(
                type(row["max_rss_node"]) is str
                and bool(row["max_rss_node"])
                and type(row["max_rss_task"]) is str
                and row["max_rss_task"].isdigit()
                and 0 < row["task_count"] <= tasks,
                "carrier calibration MaxRSS lacks exact node/task-count attribution",
            )
            max_task_rss_bytes = kib * 1024
            observed.append(
                {
                    "job_id_raw": row["job_id_raw"],
                    "reported_step_tasks": row["task_count"],
                    "max_task_rss_bytes": max_task_rss_bytes,
                    "step_memory_upper_bound_bytes": (
                        max_task_rss_bytes * row["task_count"]
                    ),
                    "node": row["max_rss_node"],
                    "task": int(row["max_rss_task"]),
                }
            )
        normalized.append(row)
    _require(observed, "carrier calibration sacct rows contain no attributable MaxRSS")
    maximum = max(item["max_task_rss_bytes"] for item in observed)
    peak_upper_bound = sum(
        item["step_memory_upper_bound_bytes"] for item in observed
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RESOURCE_MEASUREMENT_RECORD_TYPE,
        "job_id": job_id,
        "tasks": tasks,
        "sacct_units": "K",
        "sacct_rows": normalized,
        "sacct_rows_sha256": _canonical_sha256(normalized),
        "maximum_observed_task_rss_bytes": maximum,
        "maximum_observed_task_rss_locations": [
            item for item in observed if item["max_task_rss_bytes"] == maximum
        ],
        "slurm_step_memory_upper_bounds": observed,
        "peak_resident_memory_bytes": peak_upper_bound,
        "peak_resident_memory_semantics": MEMORY_SEMANTICS,
        "peak_resident_memory_exact_aggregate": False,
        "peak_resident_memory_is_conservative_upper_bound": True,
        "trustworthy_for_capacity_planning": True,
    }


def _capture_resource_measurement(job_id: str, tasks: int) -> dict[str, object]:
    last_error: Exception | None = None
    for attempt in range(ACCOUNTING_CAPTURE_ATTEMPTS):
        output = legacy._run(
            [
                "/usr/bin/sacct",
                "--clusters=frontier",
                "-j",
                job_id,
                "--units=K",
                (
                    "--format=JobIDRaw,JobName,State,ExitCode,ElapsedRaw,AllocNodes,"
                    "NTasks,MaxRSS,MaxRSSNode,MaxRSSTask,TRESUsageInMax,"
                    "TRESUsageInMaxNode,TRESUsageInMaxTask"
                ),
                "-n",
                "-P",
            ]
        )
        try:
            rows = []
            for line in output.splitlines():
                if not line.strip():
                    continue
                fields = line.split("|")
                if fields and fields[-1] == "":
                    fields.pop()
                _require(
                    len(fields) == 13,
                    "carrier calibration sacct column count drifted",
                )
                rows.append(
                    {
                        "job_id_raw": fields[0],
                        "job_name": fields[1],
                        "state": fields[2],
                        "exit_code": fields[3],
                        "elapsed_seconds": int(fields[4]),
                        "allocated_nodes": int(fields[5]),
                        "task_count": int(fields[6]) if fields[6] else 0,
                        "max_rss": fields[7],
                        "max_rss_node": fields[8],
                        "max_rss_task": fields[9],
                        "tres_usage_in_max": fields[10],
                        "tres_usage_in_max_node": fields[11],
                        "tres_usage_in_max_task": fields[12],
                    }
                )
            return build_resource_measurement(
                job_id=job_id, tasks=tasks, rows=rows
            )
        except (DriverError, ValueError) as error:
            last_error = error
            if attempt + 1 < ACCOUNTING_CAPTURE_ATTEMPTS:
                time.sleep(ACCOUNTING_CAPTURE_POLL_SECONDS)
    raise DriverError(
        "carrier calibration accounting did not expose attributable MaxRSS "
        f"within {ACCOUNTING_CAPTURE_ATTEMPTS * ACCOUNTING_CAPTURE_POLL_SECONDS} "
        "seconds"
    ) from last_error


def _validate_measurement(
    measurement: object, *, artifact_id: str, job_id: str
) -> dict[str, object]:
    _require(
        type(measurement) is dict
        and measurement.get("record_type") == RESOURCE_MEASUREMENT_RECORD_TYPE,
        f"{artifact_id}: carrier calibration resource measurement is absent",
    )
    rebuilt = build_resource_measurement(
        job_id=job_id,
        tasks=EXPECTED_TASKS[artifact_id],
        rows=measurement.get("sacct_rows"),
    )
    _require(
        _strict_equal(measurement, rebuilt),
        f"{artifact_id}: carrier calibration resource measurement drifted",
    )
    return rebuilt


def _raw_payload_bytes(receipt: Mapping[str, object]) -> int:
    inventory = receipt.get("raw_inventory")
    _require(type(inventory) is list and inventory, "carrier raw inventory is absent")
    _require(
        all(
            type(item) is dict
            and type(item.get("byte_count")) is int
            and item["byte_count"] > 0
            for item in inventory
        ),
        "carrier raw inventory byte counts are malformed",
    )
    return sum(int(item["byte_count"]) for item in inventory)


def _output_slot_count(receipt: Mapping[str, object]) -> int:
    binary = receipt.get("binary_output_indices")
    checkpoint = receipt.get("checkpoint_output_indices")
    _require(
        type(binary) is list
        and type(checkpoint) is list
        and binary == list(range(len(binary)))
        and checkpoint == list(range(len(checkpoint))),
        "carrier output chronology is not contiguous from zero",
    )
    return 10 * len(binary) + 2 * len(checkpoint) + 2


def build_execution_index(
    attempts: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    _require(
        type(attempts) in {list, tuple}
        and len(attempts) == len(EXPECTED_ARTIFACTS),
        "carrier calibration execution index requires six attempts",
    )
    normalized = []
    identities = {key: set() for key in (
        "submission_id",
        "job_id",
        "authorization_id",
        "artifact_root",
        "admission_sha256",
        "reconciliation_event_sha256",
    )}
    for offset, raw in enumerate(attempts):
        _require(type(raw) is dict, "carrier calibration attempt must be an object")
        attempt = dict(raw)
        artifact_id = EXPECTED_ARTIFACTS[offset]
        case_id = EXPECTED_CASES[artifact_id]
        authorization_id = EXPECTED_AUTHORIZATIONS[artifact_id]
        cycle_limit = EXPECTED_CYCLES[artifact_id]
        submission_id = attempt.get("submission_id")
        job_id = attempt.get("job_id")
        _require(
            attempt.get("attempt_index") == offset + 1
            and attempt.get("artifact_id") == artifact_id
            and attempt.get("source_case_id") == case_id
            and attempt.get("authorization_id") == authorization_id
            and type(submission_id) is str
            and UUID_PATTERN.fullmatch(submission_id) is not None
            and type(job_id) is str
            and job_id.isdigit()
            and attempt.get("artifact_root")
            == str(PIC_ROOT / preparation.RUN_NAMESPACE / submission_id),
            f"{artifact_id}: carrier calibration execution identity drifted",
        )
        admission = attempt.get("admission")
        facts = attempt.get("admission_facts")
        provenance = attempt.get("execution_provenance")
        inventory = attempt.get("artifact_inventory")
        event = attempt.get("reconciliation_event")
        measurement = _validate_measurement(
            attempt.get("resource_measurement"),
            artifact_id=artifact_id,
            job_id=job_id,
        )
        _require(
            type(admission) is dict
            and set(admission) == {"path", "sha256", "byte_count"}
            and SHA256_PATTERN.fullmatch(str(admission.get("sha256"))) is not None
            and type(admission.get("byte_count")) is int
            and admission["byte_count"] > 0
            and type(facts) is dict
            and type(provenance) is dict
            and set(provenance) == legacy.EXECUTION_PROVENANCE_KEYS
            and type(inventory) is dict
            and set(inventory) == {"path", "sha256", "byte_count"}
            and type(event) is dict,
            f"{artifact_id}: carrier admission or provenance binding drifted",
        )
        identity = facts.get("registered_execution_identity")
        profile = facts.get("execution_profile")
        completion = facts.get("runtime_completion")
        _require(
            facts.get("record_type")
            == registered_admission.RECORD_TYPE
            and facts.get("case_id") == case_id
            and facts.get("artifact_root") == attempt["artifact_root"]
            and type(identity) is dict
            and identity.get("submission_id") == submission_id
            and identity.get("registered_science_authorization_id")
            == authorization_id
            and identity.get("slurm_job_id") == job_id
            and type(profile) is dict
            and profile.get("kind") == "runtime_controller_overlay"
            and profile.get("artifact_id") == artifact_id
            and profile.get("source_case_id") == case_id
            and profile.get("authority") == "excluded_pilot_only"
            and profile.get("expected_stop_reason") == 1903
            and profile.get("saturation_evidence_eligible") is False
            and type(completion) is dict
            and completion.get("run_completion_status")
            == "completed_not_acceptance_eligible"
            and completion.get("problem_stop_requested") is True
            and completion.get("stop_reason_code") == "1903"
            and completion.get("runtime_controller_trigger_cycle") == cycle_limit
            and completion.get("process_exit_code") == 0
            and completion.get("scheduler_terminal_state") == "COMPLETED"
            and completion.get("trusted_execution_binding_present") is True
            and facts.get("raw_science_admission_eligible") is True
            and facts.get("saturation_evidence_eligible") is False
            and all(value is False for value in facts["authorization"].values()),
            f"{artifact_id}: carrier admitted controller facts drifted",
        )
        elapsed = event.get("elapsed_seconds")
        _require(
            event.get("event_type") == "reconciliation"
            and event.get("campaign") == preparation.CAMPAIGN
            and event.get("test_id") == case_id
            and event.get("registered_science_authorization_id")
            == authorization_id
            and event.get("submission_id") == submission_id
            and event.get("job_id") == job_id
            and event.get("artifact_dir") == attempt["artifact_root"]
            and event.get("state") == "COMPLETED"
            and event.get("reconciled") is True
            and event.get("requested_nodes") == EXPECTED_NODES[artifact_id]
            and event.get("scheduler_reported_allocated_nodes")
            == EXPECTED_NODES[artifact_id]
            and event.get("scheduler_exit_code") == "0:0"
            and type(elapsed) is int
            and 0 < elapsed <= MAXIMUM_ELAPSED_SECONDS[artifact_id]
            and event.get("billed_nodes") == EXPECTED_NODES[artifact_id]
            and math.isclose(
                float(event.get("consumed_node_hours")),
                EXPECTED_NODES[artifact_id] * elapsed / 3600.0,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            and event.get("event_sha256") == legacy._ledger_event_sha256(event),
            f"{artifact_id}: carrier reconciliation event drifted",
        )
        raw_bytes = attempt.get("raw_payload_bytes")
        artifact_bytes = attempt.get("artifact_payload_bytes")
        output_slots = attempt.get("output_slot_count")
        completed_cycles = attempt.get("completed_cycles")
        _require(
            completed_cycles == cycle_limit
            and type(raw_bytes) is int
            and raw_bytes > 0
            and type(artifact_bytes) is int
            and raw_bytes <= artifact_bytes <= MAXIMUM_STORAGE_BYTES[artifact_id]
            and type(output_slots) is int
            and output_slots >= 14
            and attempt.get("peak_resident_memory_bytes")
            == measurement["peak_resident_memory_bytes"],
            f"{artifact_id}: carrier required measurements drifted",
        )
        for key, value in (
            ("submission_id", submission_id),
            ("job_id", job_id),
            ("authorization_id", authorization_id),
            ("artifact_root", attempt["artifact_root"]),
            ("admission_sha256", admission["sha256"]),
            ("reconciliation_event_sha256", event["event_sha256"]),
        ):
            _require(
                value not in identities[key],
                f"{artifact_id}: carrier execution identity reused: {key}",
            )
            identities[key].add(value)
        normalized.append(attempt)

    by_id = {str(item["artifact_id"]): item for item in normalized}
    startup = []
    pairs = []
    for dimension, short_id, base_id, instrumented_id in (
        (2, EXPECTED_ARTIFACTS[0], EXPECTED_ARTIFACTS[1], EXPECTED_ARTIFACTS[2]),
        (3, EXPECTED_ARTIFACTS[3], EXPECTED_ARTIFACTS[4], EXPECTED_ARTIFACTS[5]),
    ):
        short = by_id[short_id]
        baseline = by_id[base_id]
        instrumented = by_id[instrumented_id]
        short_cycles = int(short["completed_cycles"])
        long_cycles = int(baseline["completed_cycles"])
        short_elapsed = int(short["reconciliation_event"]["elapsed_seconds"])
        baseline_elapsed = int(baseline["reconciliation_event"]["elapsed_seconds"])
        instrumented_elapsed = int(
            instrumented["reconciliation_event"]["elapsed_seconds"]
        )
        _require(
            short["source_case_id"]
            == baseline["source_case_id"]
            == instrumented["source_case_id"]
            and short_cycles < long_cycles
            and long_cycles == instrumented["completed_cycles"],
            f"carrier {dimension}D calibration grouping drifted",
        )
        baseline_increment = baseline_elapsed - short_elapsed
        instrumented_increment = instrumented_elapsed - short_elapsed
        _require(
            baseline_increment > 0 and instrumented_increment > 0,
            f"carrier {dimension}D steady-cycle timing is not positive",
        )
        startup.append(
            {
                "dimension": dimension,
                "artifact_id": short_id,
                "completed_cycles": short_cycles,
                "elapsed_seconds": short_elapsed,
                "peak_resident_memory_bytes": short[
                    "peak_resident_memory_bytes"
                ],
                "maximum_observed_task_rss_bytes": short[
                    "resource_measurement"
                ]["maximum_observed_task_rss_bytes"],
                "memory_is_conservative_upper_bound": True,
            }
        )
        pairs.append(
            {
                "dimension": dimension,
                "source_case_id": short["source_case_id"],
                "startup_artifact_id": short_id,
                "baseline_artifact_id": base_id,
                "instrumented_artifact_id": instrumented_id,
                "incremental_cycles": long_cycles - short_cycles,
                "baseline_incremental_elapsed_seconds": baseline_increment,
                "instrumented_incremental_elapsed_seconds": instrumented_increment,
                "baseline_seconds_per_incremental_cycle": (
                    baseline_increment / (long_cycles - short_cycles)
                ),
                "instrumented_seconds_per_incremental_cycle": (
                    instrumented_increment / (long_cycles - short_cycles)
                ),
                "instrumentation_elapsed_ratio": (
                    instrumented_increment / baseline_increment
                ),
                "resource_measurement_only": True,
                "scientific_selection_authorized": False,
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": EXECUTION_INDEX_RECORD_TYPE,
        "status": (
            "complete_registered_carrier_calibration_execution_index_"
            "non_authorizing"
        ),
        "campaign": preparation.CAMPAIGN,
        "attempt_count": len(normalized),
        "attempts": normalized,
        "attempt_bindings_sha256": _canonical_sha256(normalized),
        "startup_measurements": startup,
        "startup_measurements_sha256": _canonical_sha256(startup),
        "steady_cycle_pairs": pairs,
        "steady_cycle_pairs_sha256": _canonical_sha256(pairs),
        "total_consumed_node_hours": sum(
            float(item["reconciliation_event"]["consumed_node_hours"])
            for item in normalized
        ),
        "total_artifact_bytes": sum(
            int(item["artifact_payload_bytes"]) for item in normalized
        ),
        "total_raw_bytes": sum(int(item["raw_payload_bytes"]) for item in normalized),
        "total_output_slots": sum(
            int(item["output_slot_count"]) for item in normalized
        ),
        "saturation_evidence_eligible": False,
        "production_resource_freeze_authorized": False,
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def validate_execution_index(value: object) -> dict[str, object]:
    _require(
        type(value) is dict
        and value.get("record_type") == EXECUTION_INDEX_RECORD_TYPE,
        "carrier calibration execution index identity drifted",
    )
    rebuilt = build_execution_index(value.get("attempts"))
    _require(
        _strict_equal(value, rebuilt),
        "carrier calibration execution index derived fields drifted",
    )
    return rebuilt


def validate_execution_index_files(value: object) -> dict[str, object]:
    rebuilt = validate_execution_index(value)
    for attempt in rebuilt["attempts"]:
        admission = legacy._reopen_attempt_admission(attempt)
        legacy._validate_attempt_admission_projection(attempt, admission)
        artifact_root = Path(str(attempt["artifact_root"]))
        _require(
            legacy._artifact_payload_bytes(
                attempt["artifact_inventory"], artifact_root=artifact_root
            )
            == attempt["artifact_payload_bytes"],
            f"{attempt['artifact_id']}: carrier artifact bytes drifted",
        )
        receipt_path = (
            artifact_root
            / "analysis"
            / registered_admission.EXECUTION_RECEIPT_NAME
        )
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        _require(
            _raw_payload_bytes(receipt) == attempt["raw_payload_bytes"]
            and _output_slot_count(receipt) == attempt["output_slot_count"]
            and receipt.get("terminal_cycle") == attempt["completed_cycles"],
            f"{attempt['artifact_id']}: carrier receipt measurement drifted",
        )
    return rebuilt


class Driver(legacy.Driver):
    """Live serialized carrier calibration driver."""

    def __init__(self, arguments: argparse.Namespace) -> None:
        self.source_root = arguments.source_root.resolve(strict=True)
        self.source_commit = arguments.source_commit
        self.control_plane_version = arguments.control_plane_version
        self.final_bindings_path = arguments.final_bindings.resolve(strict=True)
        self.final_bindings_sha256 = arguments.final_bindings_sha256
        self.reviewed_policy_path = arguments.reviewed_policy.resolve(strict=True)
        self.reviewed_policy_sha256 = arguments.reviewed_policy_sha256
        self.baseline_policy_sha256 = arguments.baseline_policy_sha256
        self.baseline_promotion_sha256 = arguments.baseline_promotion_sha256
        self.control_plane_root = (
            PIC_ROOT / "control_plane" / self.control_plane_version
        )
        self.run_control_plane = self.control_plane_root / "run_control_plane.py"
        self.attestation_helper = (
            self.source_root
            / "tst/publication/capture_frontier_pre_policy_promotion_attestation.py"
        )
        self._validate_cli_identities()
        sys.path.insert(0, str(self.source_root))
        self.preparation = preparation
        self.admission = registered_admission
        self.final = self._read_bound_json(
            self.final_bindings_path,
            expected_sha256=self.final_bindings_sha256,
            label="carrier calibration final bindings",
        )
        self.reviewed_policy = self._read_bound_json(
            self.reviewed_policy_path,
            expected_sha256=self.reviewed_policy_sha256,
            label="carrier calibration reviewed policy",
        )
        self.preparation.validate_final_binding_files(self.final)
        _require(
            self.final["source_commit"] == self.source_commit
            and self.final["installed_control_plane_version"]
            == self.control_plane_version,
            "carrier calibration final bindings differ from driver identities",
        )
        self.members = self.preparation._pilot_members()
        self.campaign = self.preparation.CAMPAIGN
        self.work_root = PIC_ROOT / "jobs" / self.campaign
        self.admission_root = (
            PIC_ROOT / "analysis" / self.campaign / "case_admissions"
        )
        self.index_path = (
            PIC_ROOT
            / "analysis"
            / self.campaign
            / "q019_q023_carrier_calibration_execution_index_v1.json"
        )

    def _state_path(self, artifact_id: str) -> Path:
        _require(
            artifact_id in EXPECTED_ARTIFACTS,
            "carrier calibration state key is not an exact artifact ID",
        )
        return self.work_root / "state" / f"{artifact_id}.json"

    def _promote_policy_if_needed(self) -> None:
        policy, policy_sha256, promotion_sha256 = self._active_policy()
        if policy_sha256 == self.reviewed_policy_sha256:
            _require(
                policy == self.reviewed_policy
                and len(policy["registered_science_slices"]) == len(self.members),
                "active carrier calibration policy differs from reviewed policy",
            )
            return
        _require(
            policy_sha256 == self.baseline_policy_sha256
            and promotion_sha256 == self.baseline_promotion_sha256
            and policy["registered_science_slices"] == [],
            "active policy is neither carrier predecessor nor successor",
        )
        expected = self.preparation.materialize_q019_promotable_policy(
            baseline_policy=policy, final_bindings=self.final
        )
        _require(
            expected == self.reviewed_policy,
            "reviewed carrier policy is not the exact materialized successor",
        )
        self._wait_for_empty_queue()
        authorization_id = f"q019-carrier-cal-policy-{self.source_commit[:8]}-v1"
        attestation = self._capture_attestation(
            authorization_id, "pre_policy_promotion"
        )
        policy, policy_sha256, promotion_sha256 = self._active_policy()
        _require(
            policy_sha256 == self.baseline_policy_sha256
            and promotion_sha256 == self.baseline_promotion_sha256
            and policy["registered_science_slices"] == [],
            "carrier predecessor changed during promotion preparation",
        )
        self._controller(
            "promote_active_policy.py",
            "--reviewed-policy",
            self.reviewed_policy_path,
            "--pre-policy-promotion-attestation",
            attestation,
            "--pre-policy-promotion-authorization-id",
            authorization_id,
        )
        policy, policy_sha256, _ = self._active_policy()
        _require(
            policy_sha256 == self.reviewed_policy_sha256
            and len(policy["registered_science_slices"]) == len(self.members),
            "carrier calibration policy promotion produced the wrong successor",
        )

    def _attempt_record(
        self,
        *,
        index: int,
        member: Mapping[str, object],
        state: Mapping[str, object],
        admission_path: Path,
    ) -> dict[str, object]:
        record = super()._attempt_record(
            index=index,
            member=member,
            state=state,
            admission_path=admission_path,
        )
        artifact_root = Path(str(record["artifact_root"]))
        receipt_path = (
            artifact_root
            / "analysis"
            / self.admission.EXECUTION_RECEIPT_NAME
        )
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        measurement = _capture_resource_measurement(
            str(state["job_id"]), int(member["overlay"]["tasks"])
        )
        record.update(
            {
                "completed_cycles": int(receipt["terminal_cycle"]),
                "raw_payload_bytes": _raw_payload_bytes(receipt),
                "output_slot_count": _output_slot_count(receipt),
                "resource_measurement": measurement,
                "peak_resident_memory_bytes": measurement[
                    "peak_resident_memory_bytes"
                ],
            }
        )
        return record

    def _run_attempt(
        self, index: int, member: Mapping[str, object]
    ) -> dict[str, object]:
        artifact_id = str(member["artifact_id"])
        admission_path = self.admission_root / f"{artifact_id}.json"
        state = self._load_state(artifact_id)
        if admission_path.exists():
            _require(state is not None, f"{artifact_id}: admission lacks resume state")
            record = self.admission.validate_admission(
                json.loads(admission_path.read_text(encoding="utf-8"))
            )
            self._validate_admitted_resume_state(
                index=index,
                member=member,
                state=state,
                admission_path=admission_path,
                admission_record=record,
            )
            return self._attempt_record(
                index=index,
                member=member,
                state=state,
                admission_path=admission_path,
            )
        if state is None:
            self._wait_for_empty_queue()
            submission_id = str(uuid.uuid4())
            authorization_id = self.preparation._authorization_id(index, member)
            state = {
                "artifact_id": artifact_id,
                "source_case_id": member["source_case_id"],
                "attempt_index": index,
                "submission_id": submission_id,
                "authorization_id": authorization_id,
                "status": "materializing",
            }
            self._save_state(artifact_id, state)
            manifest = self._materialize_attempt_inputs(
                artifact_id=artifact_id,
                submission_id=submission_id,
                authorization_id=authorization_id,
            )
            state["manifest_path"] = str(manifest)
            state["status"] = "manifest_created"
            self._save_state(artifact_id, state)
        job_id = str(state.get("job_id", ""))
        if not job_id:
            recovered = self._ledger_job_for_submission(str(state["submission_id"]))
            if recovered:
                job_id = recovered
                state["job_id"] = job_id
                state["status"] = "submitted"
                self._save_state(artifact_id, state)
            else:
                _require(
                    state.get("status") == "manifest_created",
                    f"{artifact_id}: incomplete pre-submission state requires review",
                )
                self._wait_for_empty_queue()
                pre_submit = self._capture_attestation(
                    str(state["authorization_id"]), "pre_submit_wrapper"
                )
                output = legacy._run(
                    [
                        self.control_plane_root / "submit_frontier_job.sh",
                        Path(str(state["manifest_path"])),
                        pre_submit,
                    ]
                )
                match = JOB_ID_PATTERN.search(output)
                _require(match is not None, f"could not parse submission: {output}")
                job_id = match.group(1)
                state["job_id"] = job_id
                state["status"] = "submitted"
                self._save_state(artifact_id, state)
                print(f"{artifact_id}: {output}", flush=True)
        self._wait_for_terminal(job_id)
        evidence = json.loads(
            self._controller(
                "reconcile_q019_registered_execution.py",
                "--job-id",
                job_id,
                "--ledger-jsonl",
                LEDGER_JSONL,
                "--ledger-csv",
                LEDGER_CSV,
                "--receipts-jsonl",
                RECEIPTS_JSONL,
                "--mirror-jsonl",
                MIRROR_JSONL,
            )
        )
        _require(type(evidence) is dict, f"{artifact_id}: malformed reconciler evidence")
        artifact_root = (
            PIC_ROOT / self.preparation.RUN_NAMESPACE / str(state["submission_id"])
        )
        admission, _ = self.admission.derive_case_bundle(
            case_id=str(member["source_case_id"]),
            artifact_root=artifact_root,
            q043_qualification_path=Path(
                str(self.final["q043_registered_matrix_path"])
            ),
            q043_artifact_root=PIC_ROOT,
            q023_qualification_path=Path(
                str(self.final["q023_registered_matrix_path"])
            ),
        )
        _require(
            admission["execution_profile"]["artifact_id"] == artifact_id
            and admission["runtime_completion"]["stop_reason_code"] == "1903"
            and admission["runtime_completion"]["runtime_controller_trigger_cycle"]
            == int(member["cycle_limit"])
            and admission["saturation_evidence_eligible"] is False,
            f"{artifact_id}: admitted carrier stop contract drifted",
        )
        legacy._atomic_json(admission_path, admission)
        state["status"] = "admitted"
        state["admission_path"] = str(admission_path)
        state["admission_sha256"] = legacy._sha256(admission_path)
        self._save_state(artifact_id, state)
        self._wait_for_empty_queue()
        return self._attempt_record(
            index=index,
            member=member,
            state=state,
            admission_path=admission_path,
        )

    def run(self) -> int:
        self._promote_policy_if_needed()
        self.admission_root.mkdir(parents=True, exist_ok=True)
        attempts = [
            self._run_attempt(index, member)
            for index, member in enumerate(self.members, 1)
        ]
        execution_index = build_execution_index(attempts)
        if self.index_path.exists():
            existing = json.loads(self.index_path.read_text(encoding="utf-8"))
            validate_execution_index_files(existing)
            _require(
                existing == execution_index,
                "existing carrier calibration index differs",
            )
        else:
            legacy._atomic_json(self.index_path, execution_index)
        from tst.publication import (
            q019_q023_carrier_calibration_engineering_qualification_v1
            as engineering_qualification,
        )

        qualified = engineering_qualification.materialize_qualification(
            execution_index_path=self.index_path
        )
        qualification_path = (
            PIC_ROOT
            / "analysis"
            / self.campaign
            / "q019_q023_carrier_calibration_engineering_qualification_v1.json"
        )
        print(
            json.dumps(
                {
                    "attempt_count": len(attempts),
                    "execution_index_path": str(self.index_path),
                    "execution_index_sha256": legacy._sha256(self.index_path),
                    "total_consumed_node_hours": execution_index[
                        "total_consumed_node_hours"
                    ],
                    "engineering_qualification_path": str(qualification_path),
                    "engineering_qualification_sha256": legacy._sha256(
                        qualification_path
                    ),
                    "engineering_qualification_status": qualified["status"],
                    "production_resource_freeze_authorized": False,
                },
                sort_keys=True,
            ),
            flush=True,
        )
        return 0


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--control-plane-version", required=True)
    parser.add_argument("--final-bindings", required=True, type=Path)
    parser.add_argument("--final-bindings-sha256", required=True)
    parser.add_argument("--reviewed-policy", required=True, type=Path)
    parser.add_argument("--reviewed-policy-sha256", required=True)
    parser.add_argument("--baseline-policy-sha256", required=True)
    parser.add_argument("--baseline-promotion-sha256", required=True)
    return parser


def main() -> int:
    return Driver(_parser().parse_args()).run()


if __name__ == "__main__":
    raise SystemExit(main())
