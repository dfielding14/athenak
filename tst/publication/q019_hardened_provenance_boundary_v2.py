#!/usr/bin/env python3
"""Fail-closed Q019 boundary to future hardened registered provenance.

No source-local receipt or caller-supplied scheduler record is accepted. A
future integration must provide an installed-control-plane adapter that binds
the immutable registered execution, artifact inventory, Q043 admission, and
the complete runtime-source closure without reopening lexical paths.
"""

from __future__ import annotations

from typing import Mapping


RECORD_TYPE = "q019_hardened_provenance_boundary_v2"
REQUIRED_ADAPTER_RECORD_TYPE = (
    "q019_hardened_installed_control_plane_registered_admission_v1"
)
RUNTIME_COMPLETION_RECORD_TYPE = "q019_runtime_completion_status_v1"


class ProvenanceBoundaryError(ValueError):
    """Raised whenever raw evidence lacks the future hardened adapter."""


def contract() -> dict[str, object]:
    return {
        "schema_version": 2,
        "record_type": RECORD_TYPE,
        "status": "hardened_installed_control_plane_adapter_not_installed",
        "self_attested_receipts_accepted": False,
        "caller_supplied_scheduler_facts_accepted": False,
        "source_local_raw_science_analysis_enabled": False,
        "required_adapter_record_type": REQUIRED_ADAPTER_RECORD_TYPE,
        "required_bindings": [
            "immutable_registered_execution_identity",
            "scheduler_completion_and_zero_exit_from_hardened_control_plane",
            "descriptor_bound_immutable_artifact_inventory",
            "passed_Q043_registered_admission",
            "complete_runtime_source_closure",
            "exact_deck_and_analysis_closure",
            "installed_runtime_resolution_stop_controller_and_stop_event_binding",
            "structured_runtime_completion_status_and_problem_stop_reason",
            "incomplete_rejected_problem_stop_forces_nonzero_process_exit",
            "passed_loading_noise_convergence_and_applicability_gates",
        ],
        "authority": {
            "launch_authorized": False,
            "policy_authorized": False,
            "qualification_authorized": False,
            "claim_authorized": False,
            "raw_production_authorized": False,
            "nonlinear_saturation_claim_authorized": False,
        },
    }


def classify_runtime_completion(record: Mapping[str, object]) -> dict[str, object]:
    """Classify structured completion without granting raw-science admission."""
    if type(record) is not dict or set(record) != {
        "record_type",
        "run_completion_status",
        "problem_stop_requested",
        "stop_reason_code",
        "process_exit_code",
        "scheduler_terminal_state",
        "trusted_execution_binding_present",
    }:
        raise ProvenanceBoundaryError("runtime completion record shape drifted")
    if record["record_type"] != RUNTIME_COMPLETION_RECORD_TYPE:
        raise ProvenanceBoundaryError("runtime completion record identity drifted")
    status = record["run_completion_status"]
    if type(status) is not str or not status:
        raise ProvenanceBoundaryError("runtime completion status is invalid")
    incomplete = status.startswith("incomplete_rejected_")
    if incomplete and (
        record["problem_stop_requested"] is not True
        or type(record["stop_reason_code"]) is not str
        or not record["stop_reason_code"]
        or type(record["process_exit_code"]) is not int
        or record["process_exit_code"] == 0
    ):
        raise ProvenanceBoundaryError(
            "incomplete-rejected stop lacks problem-stop reason or nonzero process exit"
        )
    return {
        "run_completion_status": status,
        "incomplete_rejected": incomplete,
        "evidence_disposition": (
            "rejected_incomplete" if incomplete else "not_admitted_pending_hardened_adapter"
        ),
        "saturation_evidence_eligible": False,
        "raw_science_admission_eligible": False,
        "trusted_execution_binding_present": (
            record["trusted_execution_binding_present"] is True
        ),
    }


def validate_raw_science_admission(_: Mapping[str, object]) -> dict[str, object]:
    """Reject every raw packet until the hardened adapter is actually installed."""
    raise ProvenanceBoundaryError(
        "raw Q019 science analysis is disabled until the hardened installed "
        "control-plane/Q043 admission adapter is integrated and independently reviewed"
    )
