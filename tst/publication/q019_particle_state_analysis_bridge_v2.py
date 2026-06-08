#!/usr/bin/env python3
"""Bridge a Q019 schema-7 particle reduction into the v2 analyzer schema.

This is a structural bridge only. It does not attest to extraction, execution,
artifact provenance, or scientific validity.
"""

from __future__ import annotations

import math
from typing import Mapping


RECORD_TYPE = "q019_particle_state_analysis_bridge_v2"


class BridgeError(ValueError):
    """Raised when a particle reduction cannot enter the analyzer."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise BridgeError(message)


def build_record(
    reduction: Mapping[str, object],
    *,
    case_id: str,
    campaign_id: str,
    cycle: int,
    time: float,
    raw_restart_binding: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Return one non-authorizing analyzer record from an existing reduction."""
    _require(type(reduction) is dict, "reduction must be an object")
    _require(type(case_id) is str and case_id, "case_id is invalid")
    _require(type(campaign_id) is str and campaign_id, "campaign_id is invalid")
    _require(type(cycle) is int and cycle >= 0, "cycle is invalid")
    _require(type(time) is float and math.isfinite(time), "time is invalid")
    _require(
        reduction.get("record_type")
        == "q019_nonlinear_bell_schema7_particle_state_reduction"
        and reduction.get("qualification_effect") == "none"
        and reduction.get("launch_authorized") is False
        and reduction.get("scientific_claim_authorized") is False,
        "particle reduction identity or authority drifted",
    )
    momentum = reduction.get("particle_momentum")
    kinetic = reduction.get("particle_kinetic_energy")
    mass = reduction.get("particle_mass")
    current = reduction.get("current")
    gyroradius = reduction.get("gyroradius")
    conservation = reduction.get("conservation")
    _require(type(momentum) is dict, "particle momentum is missing")
    _require(type(kinetic) is dict, "particle kinetic energy is missing")
    _require(type(mass) is dict, "particle mass diagnostics are missing")
    _require(type(current) is dict, "particle current diagnostics are missing")
    _require(type(gyroradius) is dict and gyroradius, "particle gyroradius is missing")
    _require(type(conservation) is dict, "particle conservation is missing")
    if raw_restart_binding is not None:
        _require(
            type(raw_restart_binding) is dict
            and set(raw_restart_binding) == {"path", "sha256"},
            "raw restart binding is invalid",
        )
    return {
        "schema_version": 2,
        "record_type": RECORD_TYPE,
        "case_id": case_id,
        "campaign_id": campaign_id,
        "cycle": cycle,
        "time": time,
        "raw_restart_binding": None
        if raw_restart_binding is None
        else dict(raw_restart_binding),
        "particle_bulk_velocity": list(mass["bulk_velocity"]),
        "particle_momentum": {
            "volume_integrated_vector": list(momentum["volume_integrated_vector"]),
            "momentum_flux_tensor": momentum["momentum_flux_tensor"],
            "velocity_pressure_tensor": momentum["velocity_pressure_tensor"],
        },
        "particle_kinetic_energy": kinetic["volume_integrated"],
        "current": {
            "lab_j_over_c": list(current["deposited_current_j_over_c"]),
            "gas_frame_j_over_c": list(
                current["gas_frame_deposited_current_j_over_c"]
            ),
        },
        "gyroradius": dict(gyroradius),
        "conservation": dict(conservation),
        "authority": {
            "launch_authorized": False,
            "policy_authorized": False,
            "qualification_authorized": False,
            "claim_authorized": False,
            "raw_production_authorized": False,
            "nonlinear_saturation_claim_authorized": False,
        },
        "limitations": [
            "structural bridge only",
            "does not attest to extraction execution or artifact provenance",
        ],
    }
