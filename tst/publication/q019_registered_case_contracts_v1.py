#!/usr/bin/env python3
"""Resolve versioned Q019 base matrices and runtime-controller packets."""

from __future__ import annotations

import copy
from functools import lru_cache
import hashlib
from pathlib import Path
from typing import Mapping

from tst.publication import q019_nonlinear_bell_runtime_controller_v1 as historical_controller
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as historical
from tst.publication import q019_q023_carrier_nonlinear_bell_redesign_v1 as carrier
from tst.publication import (
    q019_q023_carrier_resource_calibration_runtime_controller_v1
    as carrier_calibration,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_SOURCE_PATH = (
    "tst/publication/q019_registered_case_contracts_v1.py"
)
HISTORICAL_FAMILY_ID = "q019_physics_first_nonlinear_bell_successor_v2"
CARRIER_FAMILY_ID = "q019_q023_carrier_nonlinear_bell_redesign_v1"
HISTORICAL_MANIFEST_PATH = (
    "inputs/publication/q019_physics_first_nonlinear_bell_successor_v2/"
    "deck_manifest.json"
)
CARRIER_MANIFEST_PATH = (
    "inputs/publication/q019_q023_carrier_nonlinear_bell_redesign_v1/"
    "deck_manifest.json"
)
HISTORICAL_CONTROLLER_PACKET_ID = "q019_nonlinear_bell_runtime_controller_v1"
HISTORICAL_CONTROLLER_MANIFEST_PATH = (
    "inputs/publication/q019_nonlinear_bell_runtime_controller_v1/"
    "deck_manifest.json"
)
HISTORICAL_CONTROLLER_ROOT = (
    "inputs/publication/q019_nonlinear_bell_runtime_controller_v1"
)
CARRIER_CALIBRATION_PACKET_ID = (
    "q019_q023_carrier_resource_calibration_runtime_controller_v1"
)
CARRIER_CALIBRATION_MANIFEST_PATH = (
    "inputs/publication/"
    "q019_q023_carrier_resource_calibration_runtime_controller_v1/"
    "deck_manifest.json"
)
CARRIER_CALIBRATION_ROOT = (
    "inputs/publication/"
    "q019_q023_carrier_resource_calibration_runtime_controller_v1"
)


class CaseContractError(ValueError):
    """Reject ambiguous or drifted Q019 case and controller contracts."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise CaseContractError(message)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _strict_equal(left: object, right: object) -> bool:
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return set(left) == set(right) and all(
            _strict_equal(left[key], right[key]) for key in left
        )
    if type(left) is list:
        return len(left) == len(right) and all(
            _strict_equal(a, b) for a, b in zip(left, right)
        )
    return left == right


@lru_cache(maxsize=1)
def _expected_cases() -> tuple[dict[str, object], ...]:
    cases = tuple(
        dict(case)
        for case in (
            *historical.expected_cases(),
            *carrier.expected_cases(),
        )
    )
    case_ids = [str(case["case_id"]) for case in cases]
    _require(
        len(case_ids) == len(set(case_ids)),
        "Q019 registered case IDs collide across versioned matrices",
    )
    return cases


def expected_cases() -> tuple[dict[str, object], ...]:
    return tuple(copy.deepcopy(case) for case in _expected_cases())


def _family(case_id: str) -> tuple[str, tuple[dict[str, object], ...]]:
    historical_ids = {
        str(case["case_id"]) for case in historical.expected_cases()
    }
    families = (
        (
            HISTORICAL_FAMILY_ID,
            tuple(
                case
                for case in _expected_cases()
                if str(case["case_id"]) in historical_ids
            ),
        ),
        (
            CARRIER_FAMILY_ID,
            tuple(
                case
                for case in _expected_cases()
                if str(case["case_id"]) not in historical_ids
            ),
        ),
    )
    matches = [
        (family_id, cases)
        for family_id, cases in families
        if any(str(case["case_id"]) == case_id for case in cases)
    ]
    _require(
        len(matches) == 1,
        "Q019 registered case identity is unknown or ambiguous",
    )
    return matches[0]


@lru_cache(maxsize=2)
def _base_manifest(family_id: str) -> dict[str, object]:
    if family_id == HISTORICAL_FAMILY_ID:
        return historical.build_deck_manifest()[0]
    _require(
        family_id == CARRIER_FAMILY_ID,
        "Q019 registered base family is unknown",
    )
    return carrier.build_deck_manifest()[0]


@lru_cache(maxsize=None)
def _resolve_case(case_id: str) -> dict[str, object]:
    family_id, cases = _family(case_id)
    case_matches = [
        dict(case) for case in cases if str(case["case_id"]) == case_id
    ]
    _require(len(case_matches) == 1, "Q019 registered case row is ambiguous")
    case = case_matches[0]
    rendered = historical.render_deck(case).encode("utf-8")
    deck_sha256 = _sha256(rendered)
    matrix_fingerprint = historical.matrix_identity_fingerprint(case)
    _require(
        case.get("matrix_identity_fingerprint") == matrix_fingerprint,
        f"{case_id}: Q019 matrix identity fingerprint drifted",
    )
    manifest = _base_manifest(family_id)
    manifest_matches = [
        item
        for item in manifest.get("decks", [])
        if type(item) is dict and item.get("case_id") == case_id
    ]
    _require(
        len(manifest_matches) == 1,
        f"{case_id}: Q019 base manifest row is absent or ambiguous",
    )
    manifest_record = manifest_matches[0]
    manifest_path = (
        HISTORICAL_MANIFEST_PATH
        if family_id == HISTORICAL_FAMILY_ID
        else CARRIER_MANIFEST_PATH
    )
    deck_root = str(Path(manifest_path).parent)
    deck_path = f"{deck_root}/{case_id}.athinput"
    _require(
        manifest_record.get("path") == deck_path
        and manifest_record.get("sha256") == deck_sha256
        and manifest_record.get("matrix_identity_fingerprint")
        == matrix_fingerprint
        and all(
            key in manifest_record
            and _strict_equal(manifest_record[key], value)
            for key, value in case.items()
        ),
        f"{case_id}: Q019 base manifest fingerprint or deck binding drifted",
    )
    return {
        "family_id": family_id,
        "case": case,
        "base_manifest_path": manifest_path,
        "base_manifest": manifest,
        "base_manifest_record": dict(manifest_record),
        "deck_path": deck_path,
        "deck_sha256": deck_sha256,
        "deck_byte_count": len(rendered),
        "matrix_identity_fingerprint": matrix_fingerprint,
    }


def resolve_case(case_id: str) -> dict[str, object]:
    return copy.deepcopy(_resolve_case(case_id))


@lru_cache(maxsize=None)
def _controller_packets(case_id: str) -> tuple[dict[str, object], ...]:
    contract = _resolve_case(case_id)
    if contract["family_id"] == HISTORICAL_FAMILY_ID:
        packet_definitions = (
            (
                HISTORICAL_CONTROLLER_PACKET_ID,
                HISTORICAL_CONTROLLER_MANIFEST_PATH,
                HISTORICAL_CONTROLLER_ROOT,
                historical_controller.build_manifest(),
                historical_controller.render_overlay,
            ),
        )
    else:
        packet_definitions = (
            (
                CARRIER_CALIBRATION_PACKET_ID,
                CARRIER_CALIBRATION_MANIFEST_PATH,
                CARRIER_CALIBRATION_ROOT,
                carrier_calibration.build_manifest(),
                carrier_calibration.render_overlay,
            ),
        )
    packets = []
    for packet_id, manifest_path, root, manifest, renderer in packet_definitions:
        artifacts = []
        for item in manifest["artifacts"]:
            if item["source_case_id"] != case_id:
                continue
            rendered = renderer(item).encode("utf-8")
            _require(
                item.get("source_deck") == contract["deck_path"]
                and item.get("source_deck_sha256")
                == contract["deck_sha256"]
                and (
                    "source_matrix_identity_fingerprint" not in item
                    or item["source_matrix_identity_fingerprint"]
                    == contract["matrix_identity_fingerprint"]
                )
                and item.get("rendered_sha256") == _sha256(rendered),
                f"{case_id}: Q019 controller source or rendered binding drifted",
            )
            artifacts.append(dict(item))
        packets.append(
            {
                "packet_id": packet_id,
                "manifest_path": manifest_path,
                "root": root,
                "manifest": manifest,
                "artifacts": tuple(artifacts),
            }
        )
    return tuple(packets)


def controller_packets(case_id: str) -> tuple[dict[str, object], ...]:
    return copy.deepcopy(_controller_packets(case_id))


def match_controller_overlay(
    case_id: str,
    *,
    immutable_parameters: Mapping[str, str] | None = None,
    rendered_sha256: str | None = None,
) -> dict[str, object]:
    _require(
        (immutable_parameters is None) != (rendered_sha256 is None),
        "Q019 controller lookup requires exactly one identity selector",
    )
    matches = []
    for packet in _controller_packets(case_id):
        for artifact in packet["artifacts"]:
            selected = (
                artifact["controller_parameters"] == immutable_parameters
                if immutable_parameters is not None
                else artifact["rendered_sha256"] == rendered_sha256
            )
            if selected:
                matches.append(
                    copy.deepcopy(
                        {
                        **artifact,
                        "controller_packet_id": packet["packet_id"],
                        "controller_manifest_path": packet["manifest_path"],
                        "controller_root": packet["root"],
                        }
                    )
                )
    _require(
        len(matches) == 1,
        "Q019 execution deck is not one exact checked-in contract or "
        "runtime-controller overlay",
    )
    return matches[0]


def required_controller_paths(case_id: str) -> frozenset[str]:
    paths: set[str] = set()
    for packet in _controller_packets(case_id):
        paths.add(str(packet["manifest_path"]))
        paths.update(
            f"{packet['root']}/{artifact['filename']}"
            for artifact in packet["artifacts"]
        )
    return frozenset(paths)


__all__ = [
    "CARRIER_FAMILY_ID",
    "CARRIER_MANIFEST_PATH",
    "CaseContractError",
    "HISTORICAL_FAMILY_ID",
    "HISTORICAL_MANIFEST_PATH",
    "REGISTRY_SOURCE_PATH",
    "controller_packets",
    "expected_cases",
    "match_controller_overlay",
    "required_controller_paths",
    "resolve_case",
]
