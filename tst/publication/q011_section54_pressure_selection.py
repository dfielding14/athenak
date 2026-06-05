#!/usr/bin/env python3
"""Validate one human-authored Q-011 Section 5.4 pressure-selection receipt.

This source-local validator does not create a receipt, rank pressure cases, or
automate selection.  It only fails closed unless a human review record binds
the immutable four-case pressure-pilot bundle and names exactly one registered
``problem/ps_p0`` choice.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
from typing import Any

if __package__:
    from . import (
        q011_section54_historical_pressure_pilot_consumer
        as historical_pressure_pilot_consumer,
    )
    from .frontier_control_plane import (
        q011_pressure_review_packet_verifier as pressure_review_packet_verifier,
    )
else:
    import q011_section54_historical_pressure_pilot_consumer as historical_pressure_pilot_consumer
    from frontier_control_plane import (
        q011_pressure_review_packet_verifier as pressure_review_packet_verifier,
    )


RECORD_TYPE = "q011_section54_pressure_selection_receipt"
SELECTION_METHOD = "human_review_only"
REGISTERED_CASES = (
    ("ps_p0_1p00", 1.0),
    ("ps_p0_0p05", 0.05),
    ("ps_p0_0p10", 0.1),
    ("ps_p0_0p20", 0.2),
)
_CASE_BY_ID = dict(REGISTERED_CASES)
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
MAX_PRESSURE_SELECTION_RECEIPT_BYTES = 1024 * 1024


class PressureSelectionReceiptError(ValueError):
    """Raised when a human pressure-selection receipt fails closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PressureSelectionReceiptError(message)


def _object(value: object, expected: set[str], label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label}: expected object")
    mapping = value
    _require(set(mapping) == expected, f"{label}: keys drifted")
    return mapping


def _list(value: object, label: str) -> list[Any]:
    _require(type(value) is list, f"{label}: expected list")
    return value


def _text(value: object, label: str) -> str:
    _require(
        type(value) is str and bool(value.strip()),
        f"{label}: expected nonempty text",
    )
    return value


def _sha256(value: object, label: str) -> str:
    _require(
        type(value) is str and _SHA256_PATTERN.fullmatch(value) is not None,
        f"{label}: malformed SHA-256",
    )
    return value


def _problem_ps_p0(value: object, label: str) -> float:
    _require(type(value) is float and math.isfinite(value), f"{label}: expected JSON float")
    return value


def _decode_json_bytes(payload: bytes, label: str) -> Any:
    def reject_constant(value: str) -> None:
        raise PressureSelectionReceiptError(f"{label}: forbidden JSON constant: {value}")

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result = {}
        for key, value in pairs:
            _require(key not in result, f"{label}: duplicate JSON key: {key!r}")
            result[key] = value
        return result

    try:
        return json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as error:
        raise PressureSelectionReceiptError(f"{label}: not valid UTF-8 JSON") from error


def _published_pressure_pilot_receipt_binding(value: object) -> dict[str, object]:
    binding = _object(
        value,
        {"path", "sha256"},
        "published_pressure_pilot_receipt",
    )
    path = _text(binding["path"], "published_pressure_pilot_receipt/path")
    expected_sha256 = _sha256(
        binding["sha256"], "published_pressure_pilot_receipt/sha256"
    )
    return {
        "path": path,
        "sha256": expected_sha256,
    }


def _pressure_gate_attestation_binding(value: object, label: str) -> dict[str, object]:
    binding = _object(value, {"path", "sha256"}, label)
    return {
        "path": _text(binding["path"], f"{label}/path"),
        "sha256": _sha256(binding["sha256"], f"{label}/sha256"),
    }


def _published_pressure_pilot_review_packet_receipt(
    value: object,
    *,
    aggregate_receipt_binding: dict[str, object],
    authorized_pic_root: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    binding = _object(
        value,
        {"path", "sha256"},
        "published_pressure_pilot_review_packet_receipt",
    )
    normalized = {
        "path": _text(
            binding["path"],
            "published_pressure_pilot_review_packet_receipt/path",
        ),
        "sha256": _sha256(
            binding["sha256"],
            "published_pressure_pilot_review_packet_receipt/sha256",
        ),
    }
    try:
        verified = (
            pressure_review_packet_verifier.consume_published_pressure_pilot_review_packet(
                normalized["path"],
                aggregate_receipt_binding=aggregate_receipt_binding,
                authorized_pic_root=authorized_pic_root,
            )
        )
    except (
        pressure_review_packet_verifier.PressureReviewPacketVerificationError,
        OSError,
    ) as error:
        raise PressureSelectionReceiptError(
            "published pressure-pilot review-packet receipt failed immutable verification"
        ) from error
    verified = _object(
        verified,
        {
            "receipt_binding",
            "aggregate_receipt_binding",
            "packet_receipt",
            "aggregate_receipt",
            "aggregate_bundle",
            "aggregate_analysis",
            "source_bindings",
            "inventory",
        },
        "verified pressure-pilot review packet",
    )
    _require(
        verified["receipt_binding"] == normalized,
        "published pressure-pilot review-packet receipt binding drifted",
    )
    _require(
        verified["aggregate_receipt_binding"] == aggregate_receipt_binding,
        "published pressure-pilot review packet is bound to a different aggregate receipt",
    )
    return normalized, verified


def _authoritative_pressure_pilot_publication(
    aggregate_receipt_binding: dict[str, object],
    review_packet_receipt_binding: dict[str, object],
    packet_publication: dict[str, object],
    *,
    authorized_pic_root: Path,
) -> dict[str, object]:
    try:
        verified = (
            historical_pressure_pilot_consumer.consume_exact_historical_production_pressure_pilot()
        )
    except (
        historical_pressure_pilot_consumer.HistoricalPressurePilotConsumerError,
        OSError,
        TypeError,
        ValueError,
    ) as error:
        raise PressureSelectionReceiptError(
            "published pressure-pilot receipt failed authoritative verification"
        ) from error
    verified = _object(
        verified,
        {
            "packet_receipt_sha256",
            "aggregate_receipt_sha256",
            "manifest_sha256",
            "analysis_result_sha256",
            "status",
        },
        "authoritative pressure-pilot publication",
    )
    packet_bundle = _object(
        packet_publication["aggregate_bundle"],
        {"path", "manifest_sha256"},
        "verified aggregate bundle",
    )
    packet_analysis = _object(
        packet_publication["aggregate_analysis"],
        {"path", "sha256"},
        "verified aggregate analysis",
    )
    _require(
        verified["packet_receipt_sha256"] == review_packet_receipt_binding["sha256"],
        "authoritative pressure-pilot packet receipt SHA-256 differs from immutable publication binding",
    )
    _require(
        verified["aggregate_receipt_sha256"] == aggregate_receipt_binding["sha256"],
        "authoritative pressure-pilot receipt SHA-256 differs from immutable publication binding",
    )
    _require(
        verified["manifest_sha256"] == packet_bundle["manifest_sha256"],
        "authoritative pressure-pilot manifest SHA-256 differs from immutable publication binding",
    )
    _require(
        verified["analysis_result_sha256"] == packet_analysis["sha256"],
        "authoritative pressure-pilot analysis SHA-256 differs from immutable publication binding",
    )
    _require(
        verified["status"] == "pass_engineering_calibration_only",
        "authoritative pressure-pilot aggregate status is not pass_engineering_calibration_only",
    )
    return verified


def _case_descriptor(value: object, index: int) -> dict[str, object]:
    label = f"case_descriptors[{index}]"
    item = _object(
        value,
        {"case_id", "problem_ps_p0", "descriptor_sha256"},
        label,
    )
    case_id = _text(item["case_id"], f"{label}/case_id")
    return {
        "case_id": case_id,
        "problem_ps_p0": _problem_ps_p0(
            item["problem_ps_p0"], f"{label}/problem_ps_p0"
        ),
        "descriptor_sha256": _sha256(
            item["descriptor_sha256"], f"{label}/descriptor_sha256"
        ),
    }


def validate_pressure_selection_receipt(
    value: object,
    *,
    authorized_pic_root: Path = (
        historical_pressure_pilot_consumer.AUTHORIZED_PRODUCTION_PIC_ROOT
    ),
) -> dict[str, object]:
    """Validate and normalize one human-authored pressure-selection receipt."""
    receipt = _object(
        value,
        {
            "schema_version",
            "record_type",
            "selection_method",
            "published_pressure_pilot_receipt",
            "published_pressure_pilot_review_packet_receipt",
            "pilot_bundle_manifest_sha256",
            "aggregate_pilot_analysis_sha256",
            "case_descriptors",
            "selected_case",
            "authoritative_reanalysis_attestation",
            "reviewer_attestation",
        },
        "pressure-selection receipt",
    )
    _require(receipt["schema_version"] == 3, "receipt schema version drifted")
    _require(type(receipt["schema_version"]) is int, "receipt schema version must be integer")
    _require(receipt["record_type"] == RECORD_TYPE, "receipt record type drifted")
    _require(receipt["selection_method"] == SELECTION_METHOD, "selection must remain human-only")

    published_binding = _published_pressure_pilot_receipt_binding(
        receipt["published_pressure_pilot_receipt"]
    )
    review_packet_binding, publication = _published_pressure_pilot_review_packet_receipt(
        receipt["published_pressure_pilot_review_packet_receipt"],
        aggregate_receipt_binding=published_binding,
        authorized_pic_root=authorized_pic_root,
    )
    authoritative_publication = _authoritative_pressure_pilot_publication(
        published_binding,
        review_packet_binding,
        publication,
        authorized_pic_root=authorized_pic_root,
    )
    descriptors = [
        _case_descriptor(item, index)
        for index, item in enumerate(
            _list(receipt["case_descriptors"], "case_descriptors")
        )
    ]
    _require(len(descriptors) == len(REGISTERED_CASES), "descriptor count drifted")
    for index, (descriptor, registered) in enumerate(zip(descriptors, REGISTERED_CASES)):
        case_id, problem_ps_p0 = registered
        _require(descriptor["case_id"] == case_id, f"case_descriptors[{index}]: case drifted")
        _require(
            descriptor["problem_ps_p0"] == problem_ps_p0,
            f"case_descriptors[{index}]: problem/ps_p0 drifted",
        )
    descriptor_sha256 = [str(item["descriptor_sha256"]) for item in descriptors]
    _require(
        len(set(descriptor_sha256)) == len(REGISTERED_CASES),
        "case descriptor SHA-256 values must be unique",
    )
    published = publication["aggregate_receipt"]
    _require(type(published) is dict, "verified aggregate receipt: expected object")
    published_bundle = _object(
        published["aggregate_bundle"],
        {"path", "manifest_sha256"},
        "published pressure-pilot receipt/aggregate_bundle",
    )
    published_analysis = _object(
        published["aggregate_analysis"],
        {"path", "sha256"},
        "published pressure-pilot receipt/aggregate_analysis",
    )
    published_cases = _list(
        published["raw_cases"],
        "published pressure-pilot receipt/raw_cases",
    )
    verified_bundle = _object(
        publication["aggregate_bundle"],
        {"path", "manifest_sha256"},
        "verified aggregate bundle",
    )
    verified_analysis = _object(
        publication["aggregate_analysis"],
        {"path", "sha256"},
        "verified aggregate analysis",
    )
    _require(
        published_bundle == verified_bundle,
        "verified aggregate bundle differs from immutable publication receipt",
    )
    _require(
        published_analysis == verified_analysis,
        "verified aggregate analysis differs from immutable publication receipt",
    )
    _require(
        len(published_cases) == len(REGISTERED_CASES),
        "published pressure-pilot receipt raw-case count drifted",
    )
    for index, (descriptor, published_case) in enumerate(
        zip(descriptors, published_cases)
    ):
        case = _object(
            published_case,
            {
                "case_id",
                "artifact_dir",
                "descriptor_path",
                "descriptor_sha256",
                "artifact_inventory_sha256",
                "runtime_artifacts",
            },
            f"published pressure-pilot receipt/raw_cases[{index}]",
        )
        _require(
            case["case_id"] == descriptor["case_id"]
            and case["descriptor_sha256"] == descriptor["descriptor_sha256"],
            f"case_descriptors[{index}]: published descriptor SHA-256 drifted",
        )
    pilot_bundle_manifest_sha256 = _sha256(
        receipt["pilot_bundle_manifest_sha256"],
        "pilot_bundle_manifest_sha256",
    )
    aggregate_pilot_analysis_sha256 = _sha256(
        receipt["aggregate_pilot_analysis_sha256"],
        "aggregate_pilot_analysis_sha256",
    )
    _require(
        pilot_bundle_manifest_sha256
        == published_bundle["manifest_sha256"]
        == verified_bundle["manifest_sha256"]
        == authoritative_publication["manifest_sha256"],
        "pilot bundle manifest SHA-256 differs from immutable publication receipt",
    )
    _require(
        aggregate_pilot_analysis_sha256
        == published_analysis["sha256"]
        == verified_analysis["sha256"]
        == authoritative_publication["analysis_result_sha256"],
        "aggregate pilot-analysis SHA-256 differs from immutable publication receipt",
    )
    reanalysis_binding = _pressure_gate_attestation_binding(
        receipt["authoritative_reanalysis_attestation"],
        "authoritative_reanalysis_attestation",
    )
    try:
        reanalysis = (
            pressure_review_packet_verifier.consume_sealed_pressure_reanalysis_attestation(
                reanalysis_binding,
                aggregate_receipt_binding=published_binding,
                packet_receipt_binding=review_packet_binding,
                pilot_bundle_manifest_sha256=pilot_bundle_manifest_sha256,
                aggregate_pilot_analysis_sha256=aggregate_pilot_analysis_sha256,
                authorized_pic_root=authorized_pic_root,
                expected_result=authoritative_publication,
            )
        )
    except (
        pressure_review_packet_verifier.PressureReviewPacketVerificationError,
        OSError,
    ) as error:
        raise PressureSelectionReceiptError(
            "authoritative pressure-pilot reanalysis attestation failed immutable verification"
        ) from error
    _require(
        reanalysis["binding"] == reanalysis_binding,
        "authoritative reanalysis attestation binding drifted",
    )

    selected = _object(
        receipt["selected_case"],
        {"case_id", "problem_ps_p0"},
        "selected_case",
    )
    selected_case_id = _text(selected["case_id"], "selected_case/case_id")
    _require(selected_case_id in _CASE_BY_ID, "selected_case is not registered")
    selected_ps_p0 = _problem_ps_p0(
        selected["problem_ps_p0"], "selected_case/problem_ps_p0"
    )
    _require(
        selected_ps_p0 == _CASE_BY_ID[selected_case_id],
        "selected_case problem/ps_p0 does not match its registered case",
    )
    selected_case = {
        "case_id": selected_case_id,
        "problem_ps_p0": selected_ps_p0,
    }
    reviewer_binding = _pressure_gate_attestation_binding(
        receipt["reviewer_attestation"],
        "reviewer_attestation",
    )
    try:
        reviewer = (
            pressure_review_packet_verifier.consume_sealed_pressure_reviewer_attestation(
                reviewer_binding,
                aggregate_receipt_binding=published_binding,
                packet_receipt_binding=review_packet_binding,
                reanalysis_verification=reanalysis,
                selected_case=selected_case,
                authorized_pic_root=authorized_pic_root,
            )
        )
    except (
        pressure_review_packet_verifier.PressureReviewPacketVerificationError,
        OSError,
    ) as error:
        raise PressureSelectionReceiptError(
            "human pressure-selection reviewer attestation failed immutable verification"
        ) from error
    _require(
        reviewer["binding"] == reviewer_binding,
        "human pressure-selection reviewer attestation binding drifted",
    )

    return {
        "schema_version": 3,
        "record_type": RECORD_TYPE,
        "selection_method": SELECTION_METHOD,
        "published_pressure_pilot_receipt": published_binding,
        "published_pressure_pilot_review_packet_receipt": review_packet_binding,
        "pilot_bundle_manifest_sha256": pilot_bundle_manifest_sha256,
        "aggregate_pilot_analysis_sha256": aggregate_pilot_analysis_sha256,
        "case_descriptors": descriptors,
        "selected_case": selected_case,
        "authoritative_reanalysis_attestation": reanalysis_binding,
        "reviewer_attestation": reviewer_binding,
    }


def canonical_json_bytes(value: object) -> bytes:
    """Serialize one receipt deterministically and reject non-finite values."""
    try:
        payload = (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (RecursionError, TypeError, ValueError) as error:
        raise PressureSelectionReceiptError("receipt is not canonical JSON") from error
    _require(
        len(payload) <= MAX_PRESSURE_SELECTION_RECEIPT_BYTES,
        "receipt exceeds the pressure-selection receipt size limit",
    )
    return payload


def validate_pressure_selection_receipt_bytes(
    payload: bytes,
    *,
    authorized_pic_root: Path = (
        historical_pressure_pilot_consumer.AUTHORIZED_PRODUCTION_PIC_ROOT
    ),
) -> dict[str, object]:
    """Decode canonical JSON bytes and validate one immutable human receipt."""
    _require(type(payload) is bytes, "receipt: expected immutable bytes")
    _require(
        len(payload) <= MAX_PRESSURE_SELECTION_RECEIPT_BYTES,
        "receipt exceeds the pressure-selection receipt size limit",
    )
    decoded = _decode_json_bytes(payload, "receipt")
    receipt = validate_pressure_selection_receipt(
        decoded,
        authorized_pic_root=authorized_pic_root,
    )
    _require(payload == canonical_json_bytes(receipt), "receipt is not canonical JSON")
    return receipt


def validate_pressure_selection_source_snapshot(
    receipt: object,
    *,
    git_commit: object,
    source_archive_sha256: object,
    helper_source_closure: object,
    authorized_pic_root: Path = (
        historical_pressure_pilot_consumer.AUTHORIZED_PRODUCTION_PIC_ROOT
    ),
) -> None:
    """Bind an accepted selection receipt to one frozen qualifying source snapshot."""
    normalized = _object(
        receipt,
        {
            "schema_version",
            "record_type",
            "selection_method",
            "published_pressure_pilot_receipt",
            "published_pressure_pilot_review_packet_receipt",
            "pilot_bundle_manifest_sha256",
            "aggregate_pilot_analysis_sha256",
            "case_descriptors",
            "selected_case",
            "authoritative_reanalysis_attestation",
            "reviewer_attestation",
        },
        "pressure-selection receipt",
    )
    _require(
        type(normalized["schema_version"]) is int
        and normalized["schema_version"] == 3
        and normalized["record_type"] == RECORD_TYPE
        and normalized["selection_method"] == SELECTION_METHOD,
        "pressure-selection receipt identity drifted",
    )
    try:
        reanalysis = (
            pressure_review_packet_verifier.consume_sealed_pressure_reanalysis_attestation(
                normalized["authoritative_reanalysis_attestation"],
                aggregate_receipt_binding=normalized["published_pressure_pilot_receipt"],
                packet_receipt_binding=normalized[
                    "published_pressure_pilot_review_packet_receipt"
                ],
                pilot_bundle_manifest_sha256=normalized["pilot_bundle_manifest_sha256"],
                aggregate_pilot_analysis_sha256=normalized[
                    "aggregate_pilot_analysis_sha256"
                ],
                authorized_pic_root=authorized_pic_root,
            )
        )
        pressure_review_packet_verifier.validate_pressure_reanalysis_source_snapshot(
            reanalysis,
            git_commit=git_commit,
            source_archive_sha256=source_archive_sha256,
            helper_source_closure=helper_source_closure,
            authorized_pic_root=authorized_pic_root,
        )
    except (
        pressure_review_packet_verifier.PressureReviewPacketVerificationError,
        OSError,
        TypeError,
        ValueError,
    ) as error:
        raise PressureSelectionReceiptError(
            "pressure-selection reanalysis differs from the frozen qualifying source snapshot"
        ) from error


__all__ = [
    "PressureSelectionReceiptError",
    "MAX_PRESSURE_SELECTION_RECEIPT_BYTES",
    "RECORD_TYPE",
    "REGISTERED_CASES",
    "SELECTION_METHOD",
    "canonical_json_bytes",
    "validate_pressure_selection_receipt",
    "validate_pressure_selection_receipt_bytes",
    "validate_pressure_selection_source_snapshot",
]
