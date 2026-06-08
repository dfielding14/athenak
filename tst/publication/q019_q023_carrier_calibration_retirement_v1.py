#!/usr/bin/env python3
"""Retire an exact completed six-slice Q019 carrier calibration allowlist."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Mapping

from tst.publication import q011_section54_pressure_pilot_execution as execution
from tst.publication import (
    q019_q023_carrier_calibration_campaign_driver_v1 as campaign,
)
from tst.publication import (
    q019_q023_carrier_calibration_engineering_qualification_v1 as qualification,
)
from tst.publication import (
    q019_q023_carrier_calibration_launch_policy_preparation_v1 as preparation,
)


CONTROL_PLANE_SOURCE_DIR = Path(__file__).resolve().parent / "frontier_control_plane"
if str(CONTROL_PLANE_SOURCE_DIR) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(CONTROL_PLANE_SOURCE_DIR))
from control_plane_common import validate_storage_policy  # type: ignore[import-not-found]  # noqa: E402


SCHEMA_VERSION = 1
RECORD_TYPE = "q019_q023_carrier_calibration_retirement_candidate_v1"
AUTHORIZATION_BOUNDARY = {
    "policy_mutation_authorized": False,
    "production_resource_freeze_authorized": False,
    "q019_qualification_authorized": False,
    "nonlinear_saturation_claim_authorized": False,
    "scientific_claim_authorized": False,
    "publication_authorized": False,
}


class RetirementError(ValueError):
    """Reject incomplete calibration evidence or a drifted active allowlist."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RetirementError(message)


def _canonical_sha256(value: object) -> str:
    payload = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _strict_equal(left: object, right: object) -> bool:
    return campaign._strict_equal(left, right)


def materialize_retired_policy(
    *,
    active_policy: Mapping[str, object],
    engineering_qualification: Mapping[str, object],
    final_bindings: Mapping[str, object],
    successor_control_plane_version: str,
    storage_preflight_binding: Path,
) -> dict[str, object]:
    """Remove exactly six completed calibration slices into a fresh generation."""
    try:
        qualified = qualification.validate_qualification(
            engineering_qualification
        )
    except qualification.QualificationError as error:
        raise RetirementError(
            "carrier retirement requires re-openable engineering qualification"
        ) from error
    _require(
        qualified["status"] == qualification.STATUS_PASS
        and qualified["decision"]["engineering_gate_pass"] is True
        and qualified["decision"]["production_resource_freeze_authorized"] is False,
        "carrier retirement requires a passing non-authorizing qualification",
    )
    _require(isinstance(active_policy, Mapping), "active policy must be an object")
    policy = copy.deepcopy(dict(active_policy))
    storage = policy.get("olcf_side_storage")
    _require(
        type(storage) is dict
        and type(storage.get("installed_control_plane_version")) is str,
        "active carrier policy lacks one installed controller binding",
    )
    predecessor_version = str(storage["installed_control_plane_version"])
    try:
        validate_storage_policy(
            copy.deepcopy(policy),
            control_plane_version=predecessor_version,
            authorized_pic_root=preparation.AUTHORIZED_ORION_ROOT,
            authorized_project_home_root=preparation.CANONICAL_PROJECT_HOME_ROOT,
        )
    except ValueError as error:
        raise RetirementError("active carrier policy failed validation") from error
    final = preparation.validate_final_binding_files(final_bindings)
    manifest, _ = preparation.build_materialization(final)
    expected_slices = [
        {**item, "status": "authorized"} for item in manifest["policy_slices"]
    ]
    observed_slices = policy.get("registered_science_slices")
    _require(
        _strict_equal(observed_slices, expected_slices),
        "active policy is not the exact authorized six-slice carrier allowlist",
    )
    index_path = Path(str(qualified["execution_index"]["path"]))
    index_payload = json.loads(index_path.read_text(encoding="utf-8"))
    index = campaign.validate_execution_index_files(index_payload)
    observed_execution_pairs = [
        (item["authorization_id"], item["source_case_id"])
        for item in index["attempts"]
    ]
    expected_pairs = [
        (item["authorization_id"], item["test_id"])
        for item in manifest["attempt_records"]
    ]
    _require(
        observed_execution_pairs == expected_pairs,
        "carrier execution index differs from active policy slices",
    )
    successor = copy.deepcopy(policy)
    successor["registered_science_slices"] = []
    try:
        successor = execution._advance_control_plane_fields(
            successor,
            control_plane_version=successor_control_plane_version,
            storage_preflight_binding=storage_preflight_binding,
            require_fresh_preflight=True,
            require_new_control_plane=True,
        )
    except (OSError, ValueError) as error:
        raise RetirementError(
            "carrier retirement requires a new controller and fresh preflight"
        ) from error
    try:
        validate_storage_policy(
            copy.deepcopy(successor),
            control_plane_version=successor_control_plane_version,
            authorized_pic_root=preparation.AUTHORIZED_ORION_ROOT,
            authorized_project_home_root=preparation.CANONICAL_PROJECT_HOME_ROOT,
        )
    except ValueError as error:
        raise RetirementError("retired carrier policy failed validation") from error
    return successor


def build_retirement_candidate(
    *,
    active_policy: Mapping[str, object],
    engineering_qualification: Mapping[str, object],
    final_bindings: Mapping[str, object],
    successor_control_plane_version: str,
    storage_preflight_binding: Path,
) -> dict[str, object]:
    successor = materialize_retired_policy(
        active_policy=active_policy,
        engineering_qualification=engineering_qualification,
        final_bindings=final_bindings,
        successor_control_plane_version=successor_control_plane_version,
        storage_preflight_binding=storage_preflight_binding,
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "status": "review_candidate_policy_mutation_not_authorized",
        "campaign": preparation.CAMPAIGN,
        "predecessor_policy_sha256": _canonical_sha256(dict(active_policy)),
        "engineering_qualification_sha256": _canonical_sha256(
            dict(engineering_qualification)
        ),
        "final_bindings_sha256": _canonical_sha256(dict(final_bindings)),
        "successor_control_plane_version": successor_control_plane_version,
        "storage_preflight_binding": {
            "path": str(storage_preflight_binding),
            "sha256": hashlib.sha256(
                storage_preflight_binding.read_bytes()
            ).hexdigest(),
        },
        "retired_slice_count": len(campaign.EXPECTED_ARTIFACTS),
        "successor_registered_science_slice_count": 0,
        "successor_policy": successor,
        "successor_policy_sha256": _canonical_sha256(successor),
        "production_resource_freeze_may_be_prepared_after_live_retirement": (
            True
        ),
        "production_resource_freeze_authorized": False,
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


__all__ = [
    "AUTHORIZATION_BOUNDARY",
    "RECORD_TYPE",
    "RetirementError",
    "build_retirement_candidate",
    "materialize_retired_policy",
]
