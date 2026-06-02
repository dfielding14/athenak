#!/usr/bin/env python3
"""Validate one human-authored Q-011 Section 5.4 pressure-selection receipt.

This source-local validator does not create a receipt, rank pressure cases, or
automate selection.  It only fails closed unless a human review record binds
the immutable four-case pressure-pilot bundle and names exactly one registered
``problem/ps_p0`` choice.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
from typing import Any

if __package__:
    from . import publish_q011_section54_pressure_pilot_bundle as pressure_pilot_publisher
else:
    import publish_q011_section54_pressure_pilot_bundle as pressure_pilot_publisher


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
_UTC_PATTERN = re.compile(
    r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}"
    r"(?:\.[0-9]{1,6})?Z"
)
_FILE_FLAGS = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)


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


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _problem_ps_p0(value: object, label: str) -> float:
    _require(type(value) is float and math.isfinite(value), f"{label}: expected JSON float")
    return value


def _reviewed_utc(value: object) -> str:
    text = _text(value, "reviewed_utc")
    _require(_UTC_PATTERN.fullmatch(text) is not None, "reviewed_utc: expected canonical UTC")
    try:
        parsed = datetime.fromisoformat(text[:-1] + "+00:00")
    except ValueError as error:
        raise PressureSelectionReceiptError("reviewed_utc: invalid UTC timestamp") from error
    _require(parsed.tzinfo == timezone.utc, "reviewed_utc: expected UTC")
    return text


def _stable_readonly_regular_bytes(path: str | Path, label: str) -> bytes:
    target = Path(path)
    _require(target.is_absolute(), f"{label}: expected absolute path")
    try:
        descriptor = os.open(target, _FILE_FLAGS)
    except OSError as error:
        raise PressureSelectionReceiptError(f"{label}: cannot open immutable file") from error
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode) and not before.st_mode & 0o222,
            f"{label}: expected read-only regular file",
        )
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        stable = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns")
        _require(
            all(getattr(before, field) == getattr(after, field) for field in stable)
            and len(payload) == after.st_size,
            f"{label}: file changed while reading",
        )
        return bytes(payload)
    finally:
        os.close(descriptor)


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
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PressureSelectionReceiptError(f"{label}: not valid UTF-8 JSON") from error


def _published_pressure_pilot_receipt(
    value: object,
    *,
    authorized_pic_root: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    binding = _object(
        value,
        {"path", "sha256"},
        "published_pressure_pilot_receipt",
    )
    path = _text(binding["path"], "published_pressure_pilot_receipt/path")
    expected_sha256 = _sha256(
        binding["sha256"], "published_pressure_pilot_receipt/sha256"
    )
    payload = _stable_readonly_regular_bytes(
        path, "published_pressure_pilot_receipt/path"
    )
    _require(
        _sha256_bytes(payload) == expected_sha256,
        "published pressure-pilot receipt SHA-256 drifted",
    )
    try:
        verified = pressure_pilot_publisher.verify_published_pressure_pilot_receipt(
            path,
            authorized_pic_root=authorized_pic_root,
        )
    except (pressure_pilot_publisher.PressurePilotPublicationError, OSError) as error:
        raise PressureSelectionReceiptError(
            "published pressure-pilot receipt failed immutable verification"
        ) from error
    _require(
        verified["receipt_sha256"] == expected_sha256,
        "published pressure-pilot receipt verifier SHA-256 drifted",
    )
    published = _decode_json_bytes(payload, "published pressure-pilot receipt")
    _require(type(published) is dict, "published pressure-pilot receipt: expected object")
    _require(
        payload == canonical_json_bytes(published),
        "published pressure-pilot receipt is not canonical JSON",
    )
    return {
        "path": path,
        "sha256": expected_sha256,
    }, {
        "verified": verified,
        "published": published,
    }


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
    authorized_pic_root: Path = pressure_pilot_publisher.AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    """Validate and normalize one human-authored pressure-selection receipt."""
    receipt = _object(
        value,
        {
            "schema_version",
            "record_type",
            "selection_method",
            "published_pressure_pilot_receipt",
            "pilot_bundle_manifest_sha256",
            "aggregate_pilot_analysis_sha256",
            "case_descriptors",
            "selected_case",
            "reviewer_identity",
            "reviewed_utc",
            "rationale",
        },
        "pressure-selection receipt",
    )
    _require(receipt["schema_version"] == 1, "receipt schema version drifted")
    _require(type(receipt["schema_version"]) is int, "receipt schema version must be integer")
    _require(receipt["record_type"] == RECORD_TYPE, "receipt record type drifted")
    _require(receipt["selection_method"] == SELECTION_METHOD, "selection must remain human-only")

    published_binding, publication = _published_pressure_pilot_receipt(
        receipt["published_pressure_pilot_receipt"],
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
    published = publication["published"]
    verified = publication["verified"]
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
        == verified["manifest_sha256"],
        "pilot bundle manifest SHA-256 differs from immutable publication receipt",
    )
    _require(
        aggregate_pilot_analysis_sha256
        == published_analysis["sha256"]
        == verified["analysis_result_sha256"],
        "aggregate pilot-analysis SHA-256 differs from immutable publication receipt",
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

    return {
        "schema_version": 1,
        "record_type": RECORD_TYPE,
        "selection_method": SELECTION_METHOD,
        "published_pressure_pilot_receipt": published_binding,
        "pilot_bundle_manifest_sha256": pilot_bundle_manifest_sha256,
        "aggregate_pilot_analysis_sha256": aggregate_pilot_analysis_sha256,
        "case_descriptors": descriptors,
        "selected_case": {
            "case_id": selected_case_id,
            "problem_ps_p0": selected_ps_p0,
        },
        "reviewer_identity": _text(receipt["reviewer_identity"], "reviewer_identity"),
        "reviewed_utc": _reviewed_utc(receipt["reviewed_utc"]),
        "rationale": _text(receipt["rationale"], "rationale"),
    }


def canonical_json_bytes(value: object) -> bytes:
    """Serialize one receipt deterministically and reject non-finite values."""
    try:
        return (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise PressureSelectionReceiptError("receipt is not canonical JSON") from error


def validate_pressure_selection_receipt_bytes(
    payload: bytes,
    *,
    authorized_pic_root: Path = pressure_pilot_publisher.AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    """Decode canonical JSON bytes and validate one immutable human receipt."""
    decoded = _decode_json_bytes(payload, "receipt")
    receipt = validate_pressure_selection_receipt(
        decoded,
        authorized_pic_root=authorized_pic_root,
    )
    _require(payload == canonical_json_bytes(receipt), "receipt is not canonical JSON")
    return receipt


__all__ = [
    "PressureSelectionReceiptError",
    "RECORD_TYPE",
    "REGISTERED_CASES",
    "SELECTION_METHOD",
    "canonical_json_bytes",
    "validate_pressure_selection_receipt",
    "validate_pressure_selection_receipt_bytes",
]
