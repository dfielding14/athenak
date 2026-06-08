#!/usr/bin/env python3
"""Qualify Q019 excluded pilots for controller and resource engineering only."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Mapping

from tst.publication import q019_excluded_pilot_campaign_driver_v1 as campaign
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as design


SCHEMA_VERSION = 1
RECORD_TYPE = "q019_excluded_pilot_engineering_qualification_v1"
STATUS_PASS = "passed_engineering_pilot_gate_non_authorizing"
STATUS_FAIL = "failed_engineering_pilot_gate_redesign_required_non_authorizing"
MAXIMUM_INSTRUMENTATION_RATIO = 1.25
MAXIMUM_INSTRUMENTATION_ADDED_SECONDS = 30
MAXIMUM_RESOURCE_FRACTION = 0.80
PILOT_CYCLE_COUNT = 20
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
    """Reject drifted pilot evidence or authority-bearing engineering decisions."""


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


def _stable_read_only(path: Path) -> bytes:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and not before.st_mode & 0o222,
            "Q019 engineering execution index is not one read-only regular file",
        )
        chunks = []
        while chunk := os.read(descriptor, 1024 * 1024):
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
        current = path.stat(follow_symlinks=False)
        identity = lambda item: (
            item.st_dev,
            item.st_ino,
            item.st_mode,
            item.st_nlink,
            item.st_size,
            item.st_mtime_ns,
            item.st_ctime_ns,
        )
        _require(
            identity(before) == identity(after)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and len(payload) == after.st_size,
            "Q019 engineering execution index changed while reading",
        )
        return payload
    finally:
        os.close(descriptor)


def build_qualification(
    execution_index: Mapping[str, object],
    *,
    execution_index_path: Path,
    execution_index_sha256: str,
) -> dict[str, object]:
    index = campaign.validate_execution_index_files(execution_index)
    _require(
        execution_index_path.is_absolute()
        and execution_index_path
        == Path(
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/analysis/"
            "q019_nonlinear_bell_registered_successor_v1/"
            "q019_excluded_pilot_execution_index_v1.json"
        )
        and campaign.SHA256_PATTERN.fullmatch(execution_index_sha256) is not None,
        "Q019 engineering qualification index binding is malformed",
    )
    attempts = {item["artifact_id"]: item for item in index["attempts"]}
    pair_records = []
    overall_pass = True
    for pair in index["pairs"]:
        dimension = int(pair["dimension"])
        baseline = attempts[str(pair["baseline_artifact_id"])]
        instrumented = attempts[str(pair["instrumented_artifact_id"])]
        baseline_event = baseline["reconciliation_event"]
        instrumented_event = instrumented["reconciliation_event"]
        baseline_elapsed = int(baseline_event["elapsed_seconds"])
        instrumented_elapsed = int(instrumented_event["elapsed_seconds"])
        added_seconds = instrumented_elapsed - baseline_elapsed
        ratio = float(pair["instrumentation_elapsed_ratio"])
        overhead_pass = (
            ratio <= MAXIMUM_INSTRUMENTATION_RATIO
            or added_seconds <= MAXIMUM_INSTRUMENTATION_ADDED_SECONDS
        )
        elapsed_ceiling = max(
            campaign.MAXIMUM_ELAPSED_SECONDS[str(baseline["artifact_id"])],
            campaign.MAXIMUM_ELAPSED_SECONDS[str(instrumented["artifact_id"])],
        )
        elapsed_fraction = max(baseline_elapsed, instrumented_elapsed) / elapsed_ceiling
        storage_ceiling = max(
            campaign.MAXIMUM_STORAGE_BYTES[str(baseline["artifact_id"])],
            campaign.MAXIMUM_STORAGE_BYTES[str(instrumented["artifact_id"])],
        )
        storage_fraction = max(
            int(baseline["artifact_payload_bytes"]),
            int(instrumented["artifact_payload_bytes"]),
        ) / storage_ceiling
        headroom_pass = (
            elapsed_fraction <= MAXIMUM_RESOURCE_FRACTION
            and storage_fraction <= MAXIMUM_RESOURCE_FRACTION
        )
        pair_pass = overhead_pass and headroom_pass
        overall_pass = overall_pass and pair_pass
        pair_records.append(
            {
                "dimension": dimension,
                "source_case_id": pair["source_case_id"],
                "baseline_artifact_id": baseline["artifact_id"],
                "instrumented_artifact_id": instrumented["artifact_id"],
                "pilot_cycle_count": PILOT_CYCLE_COUNT,
                "baseline_elapsed_seconds": baseline_elapsed,
                "instrumented_elapsed_seconds": instrumented_elapsed,
                "instrumentation_added_seconds": added_seconds,
                "instrumentation_elapsed_ratio": ratio,
                "maximum_allowed_ratio": MAXIMUM_INSTRUMENTATION_RATIO,
                "maximum_allowed_added_seconds": (
                    MAXIMUM_INSTRUMENTATION_ADDED_SECONDS
                ),
                "overhead_gate_pass": overhead_pass,
                "maximum_elapsed_fraction_of_ceiling": elapsed_fraction,
                "maximum_storage_fraction_of_ceiling": storage_fraction,
                "maximum_allowed_resource_fraction": MAXIMUM_RESOURCE_FRACTION,
                "resource_headroom_gate_pass": headroom_pass,
                "pair_engineering_gate_pass": pair_pass,
            }
        )
    decision = {
        "engineering_gate_pass": overall_pass,
        "selected_runtime_box_edge_monitor_dt": (
            design.BOX_EDGE_MONITOR_DT if overall_pass else None
        ),
        "selected_monitoring_disposition": (
            "retain_instrumented_controller_at_frozen_production_cadence"
            if overall_pass
            else "redesign_or_reduce_controller_cadence_and_repeat_excluded_pilots"
        ),
        "production_resource_freeze_recommended": overall_pass,
        "production_resource_freeze_authorized": False,
        "threshold_classification": "preregistered_engineering_not_literature",
        "literature_supplies_numeric_overhead_threshold": False,
        "scientific_selection_authorized": False,
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "status": STATUS_PASS if overall_pass else STATUS_FAIL,
        "campaign": "q019_nonlinear_bell_registered_successor_v1",
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
            "maximum_elapsed_or_storage_fraction_of_ceiling": (
                MAXIMUM_RESOURCE_FRACTION
            ),
            "logical_overhead_rule": "ratio_pass_or_added_seconds_pass",
            "both_dimension_pairs_must_pass": True,
            "threshold_source": (
                "preregistered_conservative_engineering_headroom_policy_"
                "not_a_scientific_or_literature_threshold"
            ),
        },
        "pair_qualifications": pair_records,
        "pair_qualifications_sha256": _canonical_sha256(pair_records),
        "decision": decision,
        "saturation_evidence_eligible": False,
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def validate_qualification(value: object) -> dict[str, object]:
    _require(
        type(value) is dict and value.get("record_type") == RECORD_TYPE,
        "Q019 engineering qualification identity drifted",
    )
    binding = value.get("execution_index")
    _require(type(binding) is dict, "Q019 engineering index binding is absent")
    path = Path(str(binding.get("path", "")))
    payload = _stable_read_only(path)
    _require(
        hashlib.sha256(payload).hexdigest() == binding.get("sha256"),
        "Q019 engineering execution index digest drifted",
    )
    try:
        index = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise QualificationError("Q019 engineering execution index is invalid") from error
    rebuilt = build_qualification(
        index,
        execution_index_path=path,
        execution_index_sha256=str(binding["sha256"]),
    )
    _require(
        _strict_equal(value, rebuilt),
        "Q019 engineering qualification derived fields or authority drifted",
    )
    return rebuilt


__all__ = [
    "AUTHORIZATION_BOUNDARY",
    "QualificationError",
    "RECORD_TYPE",
    "build_qualification",
    "validate_qualification",
]
