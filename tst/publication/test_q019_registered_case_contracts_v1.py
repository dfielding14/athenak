#!/usr/bin/env python3
"""Cross-version tests for registered Q019 case and controller resolution."""

from __future__ import annotations

import copy
import hashlib
import json

import pytest

from tst.publication import (
    analyze_q019_physics_first_nonlinear_bell_successor_v2 as analysis,
)
from tst.publication import (
    q019_hardened_installed_control_plane_registered_admission_v1 as admission,
)
from tst.publication import q019_nonlinear_bell_runtime_controller_v1 as controller
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as decks
from tst.publication import q019_q023_carrier_nonlinear_bell_redesign_v1 as carrier
from tst.publication import (
    q019_q023_carrier_resource_calibration_runtime_controller_v1 as calibration,
)
from tst.publication import q019_registered_case_contracts_v1 as contracts
from tst.publication import q019_registered_raw_reduction_v1 as reduction


HISTORICAL_CASE_ID = "q019-fr-runtime-initializer-ppc24-s0"
CARRIER_CASE_ID = "q019-q023-carrier-s1-onset-s0"
CARRIER_OVERLAY_ID = "q019-carrier-calibration-2d-onset-c0032-baseline"


def _carrier_payloads(case: dict[str, object]) -> dict[str, bytes]:
    return {
        relative: (admission.REPO_ROOT / relative).read_bytes()
        for relative in admission._case_required_source_paths(case)
    }


def _runtime_parameters(artifact_id: str) -> dict[str, dict[str, str]]:
    overlay = next(
        item
        for item in calibration.expected_overlays()
        if item["artifact_id"] == artifact_id
    )
    parameters = decks.parse_athinput_text(calibration.render_overlay(overlay))
    parameters[controller.CONTROLLER_BLOCK].update(
        {
            "runtime_resolution_samples": "0",
            "runtime_resolution_last_cycle": "0",
            "runtime_resolution_last_time": "0",
            "runtime_resolution_last_B_over_B0": "-1",
            "runtime_resolution_max_B_over_B0": "-1",
            "runtime_controller_triggered": "false",
            "runtime_controller_trigger_failure": "false",
            "runtime_controller_trigger_reason": "0",
            "runtime_controller_trigger_cycle": "0",
            "runtime_controller_trigger_time": "0",
            "runtime_controller_trigger_metric": "-1",
        }
    )
    return parameters


def test_registry_selects_one_exact_family_per_case_id() -> None:
    historical = contracts.resolve_case(HISTORICAL_CASE_ID)
    carrier_contract = contracts.resolve_case(CARRIER_CASE_ID)
    assert historical["family_id"] == contracts.HISTORICAL_FAMILY_ID
    assert carrier_contract["family_id"] == contracts.CARRIER_FAMILY_ID
    assert historical["base_manifest_path"] == contracts.HISTORICAL_MANIFEST_PATH
    assert carrier_contract["base_manifest_path"] == contracts.CARRIER_MANIFEST_PATH
    assert len(contracts.expected_cases()) == (
        len(decks.expected_cases()) + len(carrier.expected_cases())
    )
    assert analysis._case_map()[CARRIER_CASE_ID]["case_id"] == CARRIER_CASE_ID


def test_carrier_admission_closure_selects_exact_manifest_and_overlay() -> None:
    case = admission._case(CARRIER_CASE_ID)
    payloads = _carrier_payloads(case)
    overlay = next(
        item
        for item in calibration.build_manifest()["artifacts"]
        if item["artifact_id"] == CARRIER_OVERLAY_ID
    )
    closure = admission._candidate_analysis_closure(
        payloads,
        case=case,
        execution_deck_sha256=str(overlay["rendered_sha256"]),
    )
    assert closure["deck_manifest"] == {
        "path": contracts.CARRIER_MANIFEST_PATH,
        "sha256": hashlib.sha256(
            payloads[contracts.CARRIER_MANIFEST_PATH]
        ).hexdigest(),
        "byte_count": len(payloads[contracts.CARRIER_MANIFEST_PATH]),
        "family_id": contracts.CARRIER_FAMILY_ID,
    }
    assert closure["registered_deck"]["matrix_identity_fingerprint"] == case[
        "matrix_identity_fingerprint"
    ]
    assert closure["runtime_controller_manifest"]["packet_id"] == (
        contracts.CARRIER_CALIBRATION_PACKET_ID
    )
    assert closure["execution_deck"]["artifact_id"] == CARRIER_OVERLAY_ID
    assert closure["execution_deck"]["source_deck_sha256"] == case["sha256"]


def test_carrier_archived_matrix_fingerprint_and_deck_sha_tampering_fail() -> None:
    case = admission._case(CARRIER_CASE_ID)
    payloads = _carrier_payloads(case)
    manifest = json.loads(payloads[contracts.CARRIER_MANIFEST_PATH])
    row = next(
        item for item in manifest["decks"] if item["case_id"] == CARRIER_CASE_ID
    )
    row["matrix_identity_fingerprint"] = "0" * 64
    tampered = dict(payloads)
    tampered[contracts.CARRIER_MANIFEST_PATH] = (
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    with pytest.raises(
        admission.RegisteredAdmissionError,
        match="selected base manifest drifted",
    ):
        admission._candidate_analysis_closure(tampered, case=case)

    deck_tampered = dict(payloads)
    deck_tampered[str(case["path"])] += b"\n# substituted\n"
    with pytest.raises(
        admission.RegisteredAdmissionError,
        match="fingerprint or SHA-256 drifted",
    ):
        admission._candidate_analysis_closure(deck_tampered, case=case)


def test_carrier_raw_runtime_profile_binds_exact_calibration_packet() -> None:
    case = next(
        item for item in carrier.expected_cases() if item["case_id"] == CARRIER_CASE_ID
    )
    base, profile = reduction._execution_profile(
        _runtime_parameters(CARRIER_OVERLAY_ID),
        case=case,
    )
    assert decks.deck_semantics_payload(base) == decks.matrix_identity_payload(case)
    assert profile["base_family_id"] == contracts.CARRIER_FAMILY_ID
    assert profile["controller_packet_id"] == (
        contracts.CARRIER_CALIBRATION_PACKET_ID
    )
    assert profile["artifact_id"] == CARRIER_OVERLAY_ID
    assert profile["source_matrix_identity_fingerprint"] == case[
        "matrix_identity_fingerprint"
    ]

    hostile = _runtime_parameters(CARRIER_OVERLAY_ID)
    hostile[controller.CONTROLLER_BLOCK]["pilot_cycle_limit"] = "33"
    with pytest.raises(
        reduction.RawReductionError,
        match="exact checked-in contract",
    ):
        reduction._execution_profile(hostile, case=case)


def test_carrier_overlay_from_another_exact_case_is_rejected() -> None:
    case = next(
        item for item in carrier.expected_cases() if item["case_id"] == CARRIER_CASE_ID
    )
    hostile = copy.deepcopy(
        _runtime_parameters("q019-carrier-calibration-3d-fiducial-c0016-baseline")
    )
    with pytest.raises(
        reduction.RawReductionError,
        match="exact checked-in contract",
    ):
        reduction._execution_profile(hostile, case=case)
