#!/usr/bin/env python3
"""Qualify measured Q019 carrier calibration resources without production authority."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import stat
from typing import Mapping

from tst.publication import (
    q019_q023_carrier_calibration_campaign_driver_v1 as campaign,
)


SCHEMA_VERSION = 1
RECORD_TYPE = "q019_q023_carrier_calibration_engineering_qualification_v1"
STATUS_PASS = "passed_carrier_calibration_engineering_gate_non_authorizing"
STATUS_FAIL = "failed_carrier_calibration_engineering_gate_repeat_required"
MAXIMUM_INSTRUMENTATION_RATIO = 1.25
MAXIMUM_INSTRUMENTATION_ADDED_SECONDS = 30
MAXIMUM_RESOURCE_FRACTION = 0.80
ALLOCATED_MEMORY_BYTES_PER_NODE = 500 * 1024**3
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


class QualificationError(ValueError):
    """Reject drifted calibration evidence or authority-bearing decisions."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise QualificationError(message)


def _canonical_sha256(value: object) -> str:
    payload = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _strict_equal(left: object, right: object) -> bool:
    return campaign._strict_equal(left, right)


def _stable_read_only(path: Path) -> bytes:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and not before.st_mode & 0o222,
            "carrier calibration index is not one read-only regular file",
        )
        chunks = []
        while chunk := os.read(descriptor, 1024 * 1024):
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
        current = path.stat(follow_symlinks=False)
        _require(
            (
                before.st_dev,
                before.st_ino,
                before.st_mode,
                before.st_nlink,
                before.st_size,
                before.st_mtime_ns,
                before.st_ctime_ns,
            )
            == (
                after.st_dev,
                after.st_ino,
                after.st_mode,
                after.st_nlink,
                after.st_size,
                after.st_mtime_ns,
                after.st_ctime_ns,
            )
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and len(payload) == after.st_size,
            "carrier calibration index changed while reading",
        )
        return payload
    finally:
        os.close(descriptor)


def _canonical_index_path() -> Path:
    return (
        campaign.PIC_ROOT
        / "analysis"
        / campaign.preparation.CAMPAIGN
        / "q019_q023_carrier_calibration_execution_index_v1.json"
    )


def _canonical_qualification_path() -> Path:
    return (
        campaign.PIC_ROOT
        / "analysis"
        / campaign.preparation.CAMPAIGN
        / "q019_q023_carrier_calibration_engineering_qualification_v1.json"
    )


def build_qualification(
    execution_index: Mapping[str, object],
    *,
    execution_index_path: Path,
    execution_index_sha256: str,
) -> dict[str, object]:
    index = campaign.validate_execution_index_files(execution_index)
    expected_path = _canonical_index_path()
    _require(
        execution_index_path.is_absolute()
        and execution_index_path == expected_path
        and campaign.SHA256_PATTERN.fullmatch(execution_index_sha256) is not None,
        "carrier calibration qualification index binding is malformed",
    )
    attempts = {str(item["artifact_id"]): item for item in index["attempts"]}
    dimension_records = []
    overall_pass = True
    for startup, pair in zip(
        index["startup_measurements"], index["steady_cycle_pairs"]
    ):
        dimension = int(pair["dimension"])
        baseline = attempts[str(pair["baseline_artifact_id"])]
        instrumented = attempts[str(pair["instrumented_artifact_id"])]
        baseline_increment = int(pair["baseline_incremental_elapsed_seconds"])
        instrumented_increment = int(
            pair["instrumented_incremental_elapsed_seconds"]
        )
        ratio = float(pair["instrumentation_elapsed_ratio"])
        added_seconds = instrumented_increment - baseline_increment
        overhead_pass = (
            ratio <= MAXIMUM_INSTRUMENTATION_RATIO
            or added_seconds <= MAXIMUM_INSTRUMENTATION_ADDED_SECONDS
        )
        elapsed_fraction = max(
            int(baseline["reconciliation_event"]["elapsed_seconds"])
            / campaign.MAXIMUM_ELAPSED_SECONDS[str(baseline["artifact_id"])],
            int(instrumented["reconciliation_event"]["elapsed_seconds"])
            / campaign.MAXIMUM_ELAPSED_SECONDS[str(instrumented["artifact_id"])],
        )
        storage_fraction = max(
            int(baseline["artifact_payload_bytes"])
            / campaign.MAXIMUM_STORAGE_BYTES[str(baseline["artifact_id"])],
            int(instrumented["artifact_payload_bytes"])
            / campaign.MAXIMUM_STORAGE_BYTES[str(instrumented["artifact_id"])],
        )
        memory_capacity = (
            campaign.EXPECTED_NODES[str(startup["artifact_id"])]
            * ALLOCATED_MEMORY_BYTES_PER_NODE
        )
        dimension_attempts = (
            attempts[str(startup["artifact_id"])],
            baseline,
            instrumented,
        )
        maximum_memory_attempt = max(
            dimension_attempts,
            key=lambda item: int(item["peak_resident_memory_bytes"]),
        )
        maximum_peak_memory = int(
            maximum_memory_attempt["peak_resident_memory_bytes"]
        )
        maximum_task_rss = max(
            int(item["resource_measurement"]["maximum_observed_task_rss_bytes"])
            for item in dimension_attempts
        )
        memory_fraction = maximum_peak_memory / memory_capacity
        headroom_pass = (
            elapsed_fraction <= MAXIMUM_RESOURCE_FRACTION
            and storage_fraction <= MAXIMUM_RESOURCE_FRACTION
            and memory_fraction <= MAXIMUM_RESOURCE_FRACTION
        )
        pair_pass = overhead_pass and headroom_pass
        overall_pass = overall_pass and pair_pass
        bytes_per_output_slot = max(
            int(baseline["artifact_payload_bytes"])
            / int(baseline["output_slot_count"]),
            int(instrumented["artifact_payload_bytes"])
            / int(instrumented["output_slot_count"]),
        )
        dimension_records.append(
            {
                "dimension": dimension,
                "source_case_id": pair["source_case_id"],
                "startup_artifact_id": startup["artifact_id"],
                "baseline_artifact_id": baseline["artifact_id"],
                "instrumented_artifact_id": instrumented["artifact_id"],
                "startup_completed_cycles": startup["completed_cycles"],
                "startup_elapsed_seconds": startup["elapsed_seconds"],
                "peak_resident_memory_upper_bound_bytes": startup[
                    "peak_resident_memory_bytes"
                ],
                "maximum_observed_task_rss_bytes": startup[
                    "maximum_observed_task_rss_bytes"
                ],
                "peak_resident_memory_is_exact_aggregate": False,
                "peak_resident_memory_is_conservative_upper_bound": True,
                "maximum_peak_resident_memory_artifact_id": (
                    maximum_memory_attempt["artifact_id"]
                ),
                "maximum_peak_resident_memory_upper_bound_bytes": (
                    maximum_peak_memory
                ),
                "maximum_observed_task_rss_across_dimension_bytes": (
                    maximum_task_rss
                ),
                "incremental_cycles": pair["incremental_cycles"],
                "baseline_seconds_per_incremental_cycle": pair[
                    "baseline_seconds_per_incremental_cycle"
                ],
                "instrumented_seconds_per_incremental_cycle": pair[
                    "instrumented_seconds_per_incremental_cycle"
                ],
                "instrumentation_added_seconds": added_seconds,
                "instrumentation_elapsed_ratio": ratio,
                "maximum_allowed_ratio": MAXIMUM_INSTRUMENTATION_RATIO,
                "maximum_allowed_added_seconds": (
                    MAXIMUM_INSTRUMENTATION_ADDED_SECONDS
                ),
                "overhead_gate_pass": overhead_pass,
                "maximum_elapsed_fraction_of_ceiling": elapsed_fraction,
                "maximum_storage_fraction_of_ceiling": storage_fraction,
                "maximum_memory_fraction_of_allocated_capacity": memory_fraction,
                "maximum_allowed_resource_fraction": MAXIMUM_RESOURCE_FRACTION,
                "resource_headroom_gate_pass": headroom_pass,
                "artifact_bytes_per_output_slot_upper_estimate": math.ceil(
                    bytes_per_output_slot
                ),
                "dimension_engineering_gate_pass": pair_pass,
            }
        )
    measured_resource_model = {
        str(item["dimension"]): {
            "startup_elapsed_seconds": item["startup_elapsed_seconds"],
            "instrumented_seconds_per_cycle": item[
                "instrumented_seconds_per_incremental_cycle"
            ],
            "peak_resident_memory_upper_bound_bytes": item[
                "maximum_peak_resident_memory_upper_bound_bytes"
            ],
            "peak_resident_memory_source_artifact_id": item[
                "maximum_peak_resident_memory_artifact_id"
            ],
            "artifact_bytes_per_output_slot_upper_estimate": item[
                "artifact_bytes_per_output_slot_upper_estimate"
            ],
        }
        for item in dimension_records
    }
    decision = {
        "engineering_gate_pass": overall_pass,
        "measured_resource_model_complete": len(dimension_records) == 2,
        "measured_resource_model": measured_resource_model,
        "production_resource_freeze_recommended": overall_pass,
        "production_resource_freeze_authorized": False,
        "retire_calibration_policy_before_resource_freeze": True,
        "repeat_calibration_required": not overall_pass,
        "threshold_classification": "preregistered_engineering_not_literature",
        "scientific_selection_authorized": False,
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "status": STATUS_PASS if overall_pass else STATUS_FAIL,
        "campaign": campaign.preparation.CAMPAIGN,
        "execution_index": {
            "path": str(execution_index_path),
            "sha256": execution_index_sha256,
            "canonical_record_sha256": _canonical_sha256(index),
        },
        "thresholds": {
            "maximum_instrumentation_elapsed_ratio": (
                MAXIMUM_INSTRUMENTATION_RATIO
            ),
            "maximum_instrumentation_added_seconds": (
                MAXIMUM_INSTRUMENTATION_ADDED_SECONDS
            ),
            "maximum_elapsed_storage_or_memory_fraction": (
                MAXIMUM_RESOURCE_FRACTION
            ),
            "allocated_memory_bytes_per_node": ALLOCATED_MEMORY_BYTES_PER_NODE,
            "logical_overhead_rule": "ratio_pass_or_added_seconds_pass",
            "both_dimensions_must_pass": True,
            "threshold_source": (
                "preregistered_conservative_engineering_headroom_policy_"
                "not_a_scientific_or_literature_threshold"
            ),
        },
        "dimension_qualifications": dimension_records,
        "dimension_qualifications_sha256": _canonical_sha256(dimension_records),
        "decision": decision,
        "saturation_evidence_eligible": False,
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def validate_qualification(value: object) -> dict[str, object]:
    _require(
        type(value) is dict and value.get("record_type") == RECORD_TYPE,
        "carrier calibration qualification identity drifted",
    )
    binding = value.get("execution_index")
    _require(type(binding) is dict, "carrier execution-index binding is absent")
    path = Path(str(binding.get("path", "")))
    payload = _stable_read_only(path)
    _require(
        hashlib.sha256(payload).hexdigest() == binding.get("sha256"),
        "carrier execution-index digest drifted",
    )
    try:
        index = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise QualificationError("carrier execution index is invalid") from error
    rebuilt = build_qualification(
        index,
        execution_index_path=path,
        execution_index_sha256=str(binding["sha256"]),
    )
    _require(
        _strict_equal(value, rebuilt),
        "carrier qualification derived fields or authority drifted",
    )
    return rebuilt


def materialize_qualification(
    *,
    execution_index_path: Path | None = None,
    output_path: Path | None = None,
) -> dict[str, object]:
    index_path = (
        _canonical_index_path()
        if execution_index_path is None
        else execution_index_path.resolve(strict=True)
    )
    destination = (
        _canonical_qualification_path()
        if output_path is None
        else output_path.resolve(strict=False)
    )
    _require(
        index_path == _canonical_index_path()
        and destination == _canonical_qualification_path(),
        "carrier qualification materialization requires canonical paths",
    )
    payload = _stable_read_only(index_path)
    digest = hashlib.sha256(payload).hexdigest()
    try:
        index = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise QualificationError("carrier execution index is invalid") from error
    record = build_qualification(
        index,
        execution_index_path=index_path,
        execution_index_sha256=digest,
    )
    encoded = (
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    if destination.exists():
        existing = _stable_read_only(destination)
        _require(
            existing == encoded,
            "existing carrier qualification differs from exact rederivation",
        )
        return record
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(
        destination,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o444,
    )
    try:
        view = memoryview(encoded)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError(f"short write: {destination}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    directory = os.open(destination.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)
    _require(
        _stable_read_only(destination) == encoded,
        "materialized carrier qualification changed after write",
    )
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    record = materialize_qualification()
    path = _canonical_qualification_path()
    print(
        json.dumps(
            {
                "engineering_gate_pass": record["decision"][
                    "engineering_gate_pass"
                ],
                "qualification_path": str(path),
                "qualification_sha256": hashlib.sha256(
                    path.read_bytes()
                ).hexdigest(),
                "status": record["status"],
                "production_resource_freeze_authorized": False,
            },
            sort_keys=True,
        )
    )
    return 0


__all__ = [
    "AUTHORIZATION_BOUNDARY",
    "QualificationError",
    "RECORD_TYPE",
    "STATUS_FAIL",
    "STATUS_PASS",
    "build_qualification",
    "materialize_qualification",
    "validate_qualification",
]


if __name__ == "__main__":
    raise SystemExit(main())
