#!/usr/bin/env python3
"""Pure fail-closed helpers for Q-011 Section 5.4 restart continuation.

This module intentionally performs no execution, scheduling, artifact
discovery, or output mutation. Callers supply retained restart bytes and
deterministically ordered extracted output values.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import struct
from typing import Mapping

import numpy as np


PREREGISTRATION = (
    Path(__file__).resolve().parent
    / "readiness"
    / "q011_section54_restart_continuation_preregistration_successor_v2_2026-06-06.json"
)

PIC_RESTART_MAGIC = 0x5049435253543031
EXPECTED_RESTART_SCHEMA = 7
EXPECTED_SHOCK_LEDGER_SCHEMA = 3
EXPECTED_ESCAPE_LEDGER_SCHEMA = 2
_PIC_RESTART_MARKER = struct.pack("<Q", PIC_RESTART_MAGIC)
_PIC_METADATA_FORMAT = "<15i"
_MODEL_INTEGER_COUNT = 31
_MODEL_REAL_COUNT = 37
_REAL_BYTES = 8
_MAX_LAYOUT_COUNT = 1 << 31
_MAX_SIGNED_INT = (1 << 31) - 1
_BINARY64_EPSILON = 2.220446049250313e-16
_Q011_STATE_KIND = 1
_Q011_PHYSICAL_MODE = 4
_Q011_REAL_FIELDS_PER_PARTICLE = 26
_Q011_INTEGER_FIELDS_PER_PARTICLE = 4
_Q011_RESTART_FINGERPRINT_SCHEMA = "athenak_pic_parallel_shock_restart_controls_v2"
_RESTART_FINGERPRINT_PATTERN = re.compile(r"v1:[0-9a-f]{16}")
_Q011_IPM = 6
_Q011_IPWT = 22
_Q011_IPT_BIRTH = 25
_Q011_PGID = 0
_Q011_PTAG = 1
_Q011_PSP = 2
_Q011_PCRSOURCE = 3
_Q011_SOURCE_INITIAL = 0
_Q011_SOURCE_SHOCK_INJECTED = 1
_HISTORICAL_PREREGISTRATION_CANONICAL_SHA256 = (
    "bb514439ac95d541559cccdf8533f56d8cc6bc15f476638d82462803dbe3e3bf"
)

_INTEGER_LEDGER_FIELDS = (
    "ps_cr_ledger_schema",
    "ps_injection_tag_floor",
    "ps_next_tag",
)
_BOOLEAN_LEDGER_FIELDS = (
    "ps_cr_ledger_complete",
    "ps_removed_excluded_early_cohort",
    "ps_tag_seeded",
)
_REAL_LEDGER_FIELDS = (
    "ps_mass_reservoir_global",
    "ps_injected_cr_count_global",
    "ps_injected_cr_mass_global",
    "ps_injected_cr_momentum_x1_global",
    "ps_injected_cr_momentum_x2_global",
    "ps_injected_cr_momentum_x3_global",
    "ps_injected_cr_energy_global",
    "ps_removed_cr_count_global",
    "ps_removed_cr_mass_global",
    "ps_removed_cr_momentum_x1_global",
    "ps_removed_cr_momentum_x2_global",
    "ps_removed_cr_momentum_x3_global",
    "ps_removed_cr_energy_global",
)
STARTUP_SHOCK_LEDGER_FIELDS = (
    "ps_cr_ledger_schema",
    "ps_cr_ledger_complete",
    "ps_mass_reservoir_global",
    "ps_injected_cr_count_global",
    "ps_injected_cr_mass_global",
    "ps_injected_cr_momentum_x1_global",
    "ps_injected_cr_momentum_x2_global",
    "ps_injected_cr_momentum_x3_global",
    "ps_injected_cr_energy_global",
    "ps_removed_excluded_early_cohort",
    "ps_removed_cr_count_global",
    "ps_removed_cr_mass_global",
    "ps_removed_cr_momentum_x1_global",
    "ps_removed_cr_momentum_x2_global",
    "ps_removed_cr_momentum_x3_global",
    "ps_removed_cr_energy_global",
    "ps_tag_seeded",
    "ps_injection_tag_floor",
    "ps_next_tag",
)
_NONNEGATIVE_REAL_LEDGER_FIELDS = (
    "ps_mass_reservoir_global",
    "ps_injected_cr_count_global",
    "ps_injected_cr_mass_global",
    "ps_injected_cr_energy_global",
    "ps_removed_cr_count_global",
    "ps_removed_cr_mass_global",
    "ps_removed_cr_energy_global",
)
_COUNT_LEDGER_FIELDS = (
    "ps_injected_cr_count_global",
    "ps_removed_cr_count_global",
)
_ESCAPE_INTEGER_LEDGER_FIELDS = (
    "ps_escape_ledger_schema",
    "ps_escape_audit_calls",
)
_ESCAPE_BOOLEAN_LEDGER_FIELDS = ("ps_escape_ledger_complete",)
_ESCAPE_REAL_LEDGER_FIELDS = (
    "ps_escape_last_audit_time",
    "ps_escaped_injected_cr_count_global",
    "ps_escaped_injected_cr_mass_global",
    "ps_escaped_injected_cr_momentum_x1_global",
    "ps_escaped_injected_cr_momentum_x2_global",
    "ps_escaped_injected_cr_momentum_x3_global",
    "ps_escaped_injected_cr_energy_global",
    "ps_escaped_initial_cr_count_global",
    "ps_escaped_injected_cr_term_count_global",
    "ps_escaped_injected_cr_abs_mass_global",
    "ps_escaped_injected_cr_abs_momentum_x1_global",
    "ps_escaped_injected_cr_abs_momentum_x2_global",
    "ps_escaped_injected_cr_abs_momentum_x3_global",
    "ps_escaped_injected_cr_abs_energy_global",
)
ESCAPE_LEDGER_FIELDS = (
    *_ESCAPE_INTEGER_LEDGER_FIELDS,
    *_ESCAPE_BOOLEAN_LEDGER_FIELDS,
    *_ESCAPE_REAL_LEDGER_FIELDS,
)
_INTEGER_PATTERN = re.compile(r"-?[0-9]+")
_REAL_PATTERN = re.compile(r"-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
_RESTART_CONTROL_KINDS = (
    ("ps_rho0", "real"),
    ("ps_p0", "real"),
    ("ps_u0", "real"),
    ("ps_b0", "real"),
    ("ps_eta", "real"),
    ("ps_vinj_over_u0", "real"),
    ("ps_inject_half_width_cells", "real"),
    ("ps_inject_t_start", "real"),
    ("ps_inject_t_stop", "real"),
    ("ps_remove_birth_time_before", "real"),
    ("ps_shock_speed_model", "integer"),
    ("ps_refine_curv", "real"),
    ("ps_derefine_curv", "real"),
    ("ps_rho_floor_frac", "real"),
    ("ps_p_floor_frac", "real"),
    ("ps_enable_injection", "integer"),
    ("ps_enable_subtraction", "integer"),
    ("ps_enable_curvature_amr", "integer"),
    ("ps_test_source_transaction_terms_override", "integer"),
    ("ps_test_source_transaction_terms", "real"),
    ("ps_inject_species", "integer"),
    ("ps_inject_seed", "integer"),
    ("ps_particle_mass", "real"),
    ("ps_particle_charge", "real"),
    ("ps_particle_q_over_m", "real"),
    ("ps_particle_macro_mass", "real"),
    ("ps_particle_momentum_state", "integer"),
    ("ps_particle_light_speed", "real"),
    ("ps_enable_frame_tracking", "integer"),
    ("ps_frame_mode", "integer"),
    ("ps_frame_t_start", "real"),
    ("ps_frame_t_ramp", "real"),
    ("ps_frame_vfrac", "real"),
    ("ps_frame_dv_max", "real"),
    ("ps_frame_apply_to_particles", "integer"),
    ("ps_frame_apply_to_inflow", "integer"),
    ("ps_frame_require_uniform", "integer"),
    ("ps_recenter_x_target", "real"),
    ("ps_recenter_x_trigger", "real"),
    ("ps_recenter_dx1", "real"),
    ("ps_recenter_shift_cells", "integer"),
    ("ps_recenter_vshock_model", "real"),
    ("ps_shock_speed", "real"),
    ("ps_xshock0", "real"),
    ("ps_use_2d3v", "integer"),
)

_COMPARISON_TOLERANCES = {
    "rho_bin": 1.0e-12,
    "bmag_bin": 1.0e-12,
    "prtcl_jx_bin": 1.0e-12,
    "j2_bin": 1.0e-12,
    "prtcl_all_pvtk_integer_payload": 0.0,
    "prtcl_all_pvtk_float_payload": 1.0e-06,
}
_COMPARISON_FIELD_KINDS = {
    "rho_bin": "float",
    "bmag_bin": "float",
    "prtcl_jx_bin": "float",
    "j2_bin": "float",
    "prtcl_all_pvtk_integer_payload": "integer",
    "prtcl_all_pvtk_float_payload": "float",
}
_RETAINED_OUTPUT_NOMINAL_SLOTS = [
    600.0,
    700.0,
    800.0,
    900.0,
    1000.0,
    1100.0,
    1200.0,
]

_EXPECTED_PREREGISTRATION = {
    "record_type": "q011_section54_restart_continuation_preregistration",
    "schema_version": 3,
    "date": "2026-06-06",
    "gate": "Q-011",
    "claim_id": "CLAIM-PAPER-SHOCK-001",
    "qualification_effect": (
        "policy_freeze_only_no_execution_authorization_no_claim_closure"
    ),
    "predecessor_record": (
        "tst/publication/readiness/"
        "q011_section54_restart_continuation_preregistration_successor_2026-06-06.json"
    ),
    "predecessor_sha256": (
        "3ad3e05a34e7fedf578bfea3424fa98a62ece52d4cda827d010d9840a715376f"
    ),
    "scope": (
        "Versioned bounded successor for one future Q-011 Section 5.4 "
        "restart-continuation carrier. It selects the checkpoint and paired "
        "continued outputs by preregistered nominal slots while binding exact "
        "canonical observed committed cycle and time metadata, raw schema-7 "
        "particle model bytes, the recomputed C++ continuation-control "
        "fingerprint, active particle source cohorts, and cancellation-aware "
        "schema-2 physical-escape metadata. Payload tolerances remain "
        "payload-only. This record contains no result and authorizes no "
        "scheduler call."
    ),
    "schema_contract": {
        "schema_style": "self_contained_exact_key_policy",
        "unknown_keys": "reject",
        "missing_keys": "reject",
        "numeric_aliases": "reject_boolean_nonfinite_and_numeric_string_aliases",
        "policy_drift": "requires_versioned_successor_before_execution",
    },
    "helper_binding": {
        "module": "tst/publication/q011_section54_restart.py",
        "role": "pure_fail_closed_policy_helper_only",
        "scheduler_calls": "forbidden",
        "artifact_discovery": "forbidden_callers_supply_retained_bytes_and_values",
    },
    "restart_payload_probe": {
        "pic_restart_magic_hex": "0x5049435253543031",
        "byte_order": "little",
        "restart_schema": 7,
        "metadata_integer_count": 15,
        "model_integer_count": 31,
        "model_real_count": 37,
        "real_bytes": 8,
        "minimum_meshblock_count": 1,
        "minimum_real_fields_per_particle": 1,
        "minimum_integer_fields_per_particle": 4,
        "minimum_particle_count": 1,
        "q011_exact_real_fields_per_particle": 26,
        "q011_exact_integer_fields_per_particle": 4,
        "q011_state_kind": 1,
        "q011_physical_mode": 4,
        "raw_model_policy": (
            "bind_exact_state_kind_physical_mode_light_speed_model_integer_"
            "array_and_model_real_array"
        ),
        "active_particle_policy": (
            "derive_and_validate_every_raw_particle_source_tag_species_q_over_m_"
            "weight_and_birth_cohort_require_zero_active_initial_particles_after_"
            "startup_removal_and_reconcile_injected_removed_and_escape_counts"
        ),
    },
    "startup_shock_ledger": {
        "parameter_block": "problem",
        "ledger_schema": 3,
        "accepted_boolean_literals": ["0", "1", "false", "true"],
        "required_fields": list(STARTUP_SHOCK_LEDGER_FIELDS),
        "checkpoint_state_requirements": {
            "ps_cr_ledger_complete": True,
            "ps_removed_excluded_early_cohort": True,
            "ps_tag_seeded": True,
        },
        "identity_policy": (
            "require_exact_typed_identity_between_uninterrupted_and_continued_"
            "comparison_bindings"
        ),
    },
    "particle_escape_ledger": {
        "parameter_block": "problem",
        "ledger_schema": 2,
        "required_fields": list(ESCAPE_LEDGER_FIELDS),
        "comparison_metadata_policy": (
            "retain_per_component_absolute_contribution_sums_and_term_count_for_"
            "cancellation_aware_future_generic_vs_reason_coded_escape_crosscheck"
        ),
    },
    "restart_control_binding": {
        "fingerprint_schema": _Q011_RESTART_FINGERPRINT_SCHEMA,
        "policy": (
            "recompute_the_cpp_fingerprint_from_the_complete_typed_control_set_"
            "and_bind_the_controls_in_addition_to_the_digest"
        ),
    },
    "continuation_contract": {
        "checkpoint_nominal_slot_omega0_inverse": 500.0,
        "checkpoint_selection_policy": (
            "select_the_committed_checkpoint_assigned_to_nominal_slot_500_and_"
            "bind_its_exact_canonical_observed_committed_cycle_and_time"
        ),
        "retained_output_nominal_slots_after_checkpoint_omega0_inverse": list(
            _RETAINED_OUTPUT_NOMINAL_SLOTS
        ),
        "retained_output_pairing_policy": (
            "pair_uninterrupted_and_continued_outputs_by_nominal_slot_and_"
            "require_exact_observed_committed_cycle_and_time_parity_before_"
            "payload_comparison"
        ),
        "observed_commit_metadata_contract": {
            "cycle": "canonical_nonnegative_integer",
            "time_omega0_inverse": "canonical_nonnegative_finite_float",
            "sequence": (
                "strictly_increasing_after_checkpoint_and_across_retained_"
                "nominal_slots"
            ),
            "parity": (
                "exact_typed_identity_between_uninterrupted_and_continued_"
                "paired_outputs"
            ),
        },
        "comparison_tolerances_max_absolute_difference": dict(
            _COMPARISON_TOLERANCES
        ),
        "comparison_field_kinds": dict(_COMPARISON_FIELD_KINDS),
        "comparison_value_order": (
            "caller_supplied_deterministically_ordered_flat_values_per_field"
        ),
        "comparison_policy": (
            "require_identical_checkpoint_observed_commit_restart_schema_"
            "startup_cohort_ledger_escape_ledger_raw_particle_model_restart_"
            "controls_active_particle_cohort_retained_nominal_slots_and_"
            "tolerances_then_exact_paired_output_observed_commit_parity_before_"
            "payload_field_comparison"
        ),
        "tolerance_boundary": (
            "These are AthenaK deterministic restart-continuation payload-only "
            "release screens selected before execution, not manuscript "
            "tolerances. They never apply to nominal slots, observed committed "
            "cycles or observed committed times."
        ),
    },
    "execution_policy": {
        "status": (
            "blocked_until_immutable_execution_record_binds_continuation_"
            "execution_fields"
        ),
        "required_before_continuation_execution": [
            "checkpoint_nominal_slot_omega0_inverse",
            "checkpoint_observed_committed_cycle",
            "checkpoint_observed_committed_time_omega0_inverse",
            "restart_schema",
            "startup_shock_ledger",
            "particle_escape_ledger",
            "raw_particle_model",
            "restart_control_binding",
            "active_particle_cohort",
            "retained_output_nominal_slots_after_checkpoint_omega0_inverse",
            "comparison_tolerances_max_absolute_difference",
            "clean_candidate_git_commit",
            "clean_frontier_executable_sha256",
            "qualifying_input_deck_sha256",
            "authorized_orion_campaign_root",
            "registered_frontier_submission_policy",
        ],
        "required_after_execution_before_parity_result": [
            "paired_output_observed_committed_cycle_and_time",
        ],
        "campaign_results_inspected": False,
        "scheduler_calls_authorized_by_this_record": False,
        "frontier_execution_authorized_by_this_record": False,
    },
    "limitations": [
        "This tranche is a non-executing policy and helper freeze only.",
        "It does not create a scheduler wrapper or authorize a Frontier submission.",
        "It does not close Q-011 qualification, independent recompute or external review.",
        (
            "The runtime still requires one integer MPI event-probe reduction per "
            "paper-VL2 stage; the larger 14-real payload reduction is event-only "
            "and production performance remains pilot-gated."
        ),
    ],
}

_BINDING_KEYS = {
    "checkpoint_nominal_slot_omega0_inverse",
    "checkpoint_observed_committed_cycle",
    "checkpoint_observed_committed_time_omega0_inverse",
    "restart_schema",
    "startup_shock_ledger",
    "particle_escape_ledger",
    "raw_particle_model",
    "restart_control_binding",
    "active_particle_cohort",
    "retained_output_nominal_slots_after_checkpoint_omega0_inverse",
    "comparison_tolerances_max_absolute_difference",
}
_RAW_PARTICLE_MODEL_KEYS = {
    "state_kind",
    "physical_mode",
    "cr_light_speed",
    "model_ints",
    "model_reals",
}
_RESTART_CONTROL_BINDING_KEYS = {
    "fingerprint_schema",
    "stored_fingerprint",
    "controls",
}
_ACTIVE_PARTICLE_COHORT_KEYS = {
    "particle_count",
    "initial_count",
    "shock_injected_count",
    "minimum_tag",
    "maximum_tag",
}
_OBSERVATION_KEYS = {"binding", "outputs_after_checkpoint"}
_OUTPUT_KEYS = {
    "nominal_slot_omega0_inverse",
    "observed_committed_cycle",
    "observed_committed_time_omega0_inverse",
    "fields",
}


class RestartPolicyError(ValueError):
    """Raised when Q-011 restart policy or retained payloads fail closed."""


@dataclass(frozen=True)
class RestartPayloadProbe:
    """Validated schema-7 particle restart layout."""

    restart_schema: int
    meshblock_count: int
    real_fields_per_particle: int
    integer_fields_per_particle: int
    state_kind: int
    physical_mode: int
    cr_light_speed: float
    model_ints: tuple[int, ...]
    model_reals: tuple[float, ...]
    particle_count: int
    particle_real_offset: int
    particle_integer_offset: int
    payload_end_offset: int


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RestartPolicyError(message)


def _keys(value: object, expected: set[str], label: str) -> dict[str, object]:
    _require(type(value) is dict, f"{label}: expected object")
    mapping = value
    _require(set(mapping) == expected, f"{label}: schema drift")
    return mapping


def _strict_equal(actual: object, expected: object, label: str) -> None:
    _require(type(actual) is type(expected), f"{label}: scalar type drift")
    if type(expected) is dict:
        actual_mapping = _keys(actual, set(expected), label)
        for key, expected_value in expected.items():
            _strict_equal(actual_mapping[key], expected_value, f"{label}/{key}")
    elif type(expected) is list:
        _require(len(actual) == len(expected), f"{label}: list length drift")
        for index, (actual_value, expected_value) in enumerate(zip(actual, expected)):
            _strict_equal(actual_value, expected_value, f"{label}[{index}]")
    elif type(expected) is float:
        _require(math.isfinite(actual), f"{label}: non-finite number")
        _require(actual == expected, f"{label}: numeric drift")
    else:
        _require(actual == expected, f"{label}: value drift")


def _canonical_observed_cycle(value: object, label: str) -> int:
    _require(type(value) is int, f"{label}: expected canonical integer")
    _require(value >= 0, f"{label}: expected nonnegative cycle")
    return value


def _canonical_observed_time(value: object, label: str) -> float:
    _require(type(value) is float, f"{label}: expected canonical float")
    _require(math.isfinite(value), f"{label}: expected finite float")
    _require(value >= 0.0, f"{label}: expected nonnegative time")
    return value


def _reject_json_constant(value: str) -> None:
    raise RestartPolicyError(f"JSON constant is forbidden: {value}")


def _reject_duplicate_json_keys(
    pairs: list[tuple[str, object]],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        _require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def decode_preregistration(payload: str) -> dict[str, object]:
    """Decode strict preregistration JSON without aliases or duplicate keys."""
    _require(type(payload) is str, "preregistration JSON: expected text")
    try:
        decoded = json.loads(
            payload,
            parse_constant=_reject_json_constant,
            object_pairs_hook=_reject_duplicate_json_keys,
        )
    except json.JSONDecodeError as exc:
        raise RestartPolicyError("invalid preregistration JSON") from exc
    _require(type(decoded) is dict, "preregistration: expected object")
    return decoded


def validate_preregistration(preregistration: object) -> None:
    """Validate the active successor or an exact protected historical packet."""
    if type(preregistration) is dict and preregistration.get("schema_version") == 2:
        try:
            canonical = (json.dumps(preregistration, indent=2) + "\n").encode("utf-8")
        except (TypeError, ValueError) as exc:
            raise RestartPolicyError("historical preregistration is not canonical") from exc
        _require(
            hashlib.sha256(canonical).hexdigest()
            == _HISTORICAL_PREREGISTRATION_CANONICAL_SHA256,
            "historical preregistration drifted",
        )
        return
    _strict_equal(preregistration, _EXPECTED_PREREGISTRATION, "preregistration")


def load_preregistration(path: Path = PREREGISTRATION) -> dict[str, object]:
    """Load and validate the checked-in frozen preregistration."""
    preregistration = decode_preregistration(path.read_text(encoding="utf-8"))
    validate_preregistration(preregistration)
    return preregistration


def frozen_preregistration() -> dict[str, object]:
    """Return an isolated copy of the frozen preregistration."""
    return copy.deepcopy(_EXPECTED_PREREGISTRATION)


def _validated_preregistration(
    preregistration: Mapping[str, object] | None,
) -> dict[str, object]:
    if preregistration is None:
        return frozen_preregistration()
    _strict_equal(
        preregistration,
        _EXPECTED_PREREGISTRATION,
        "active restart-continuation preregistration",
    )
    return copy.deepcopy(preregistration)


def _unpack_from(
    format_string: str, payload: bytes, offset: int, label: str, source: str
) -> tuple[object, ...]:
    size = struct.calcsize(format_string)
    _require(offset >= 0 and offset + size <= len(payload), f"{source}: truncated {label}")
    return struct.unpack_from(format_string, payload, offset)


def probe_schema7_restart_payload(
    payload: bytes, *, source: str = "<restart-bytes>"
) -> RestartPayloadProbe:
    """Probe and validate the schema-7 particle section of retained restart bytes."""
    _require(type(payload) is bytes, f"{source}: restart payload must be bytes")
    marker_offset = payload.find(_PIC_RESTART_MARKER)
    _require(marker_offset >= 0, f"{source}: particle restart marker not found")
    _require(
        payload.find(_PIC_RESTART_MARKER, marker_offset + 1) < 0,
        f"{source}: multiple particle restart markers found",
    )
    offset = marker_offset + len(_PIC_RESTART_MARKER)
    metadata = _unpack_from(_PIC_METADATA_FORMAT, payload, offset, "PIC metadata", source)
    offset += struct.calcsize(_PIC_METADATA_FORMAT)
    (
        restart_schema,
        meshblock_count,
        real_fields,
        integer_fields,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        state_kind,
        physical_mode,
    ) = metadata
    _require(restart_schema == EXPECTED_RESTART_SCHEMA, f"{source}: restart schema is not 7")
    _require(
        type(meshblock_count) is int and 1 <= meshblock_count < _MAX_LAYOUT_COUNT,
        f"{source}: invalid particle restart meshblock count",
    )
    _require(
        type(real_fields) is int and 1 <= real_fields < _MAX_LAYOUT_COUNT,
        f"{source}: invalid particle real-field count",
    )
    _require(
        type(integer_fields) is int and 4 <= integer_fields < _MAX_LAYOUT_COUNT,
        f"{source}: invalid particle integer-field count",
    )

    cr_light_speed = _unpack_from("<d", payload, offset, "particle light speed", source)[0]
    offset += _REAL_BYTES
    model_ints = _unpack_from(
        f"<{_MODEL_INTEGER_COUNT}i", payload, offset, "PIC model integers", source
    )
    offset += _MODEL_INTEGER_COUNT * struct.calcsize("<i")
    model_reals = _unpack_from(
        f"<{_MODEL_REAL_COUNT}d", payload, offset, "PIC model reals", source
    )
    offset += _MODEL_REAL_COUNT * _REAL_BYTES
    _require(state_kind in {0, 1}, f"{source}: invalid particle state kind")
    _require(0 <= physical_mode <= 4, f"{source}: invalid PIC physical mode")
    _require(
        type(cr_light_speed) is float
        and math.isfinite(cr_light_speed)
        and cr_light_speed > 0.0,
        f"{source}: invalid particle light speed",
    )
    _require(
        all(type(value) is int for value in model_ints),
        f"{source}: invalid PIC model integer metadata",
    )
    _require(
        all(type(value) is float and math.isfinite(value) for value in model_reals),
        f"{source}: invalid PIC model real metadata",
    )
    particle_count = _unpack_from("<Q", payload, offset, "particle count", source)[0]
    offset += struct.calcsize("<Q")
    _require(
        type(particle_count) is int and 1 <= particle_count < _MAX_LAYOUT_COUNT,
        f"{source}: invalid particle count",
    )
    meshblock_counts = _unpack_from(
        f"<{meshblock_count}i", payload, offset, "particle MeshBlock counts", source
    )
    _require(
        all(count >= 0 for count in meshblock_counts),
        f"{source}: particle MeshBlock count table contains a negative count",
    )
    _require(
        sum(meshblock_counts) == particle_count,
        f"{source}: particle MeshBlock count table is inconsistent",
    )
    particle_real_offset = offset + meshblock_count * struct.calcsize("<i")
    particle_integer_offset = (
        particle_real_offset + particle_count * real_fields * _REAL_BYTES
    )
    payload_end_offset = (
        particle_integer_offset
        + particle_count * integer_fields * struct.calcsize("<i")
    )
    _require(
        payload_end_offset <= len(payload),
        f"{source}: truncated particle restart payload",
    )
    return RestartPayloadProbe(
        restart_schema=restart_schema,
        meshblock_count=meshblock_count,
        real_fields_per_particle=real_fields,
        integer_fields_per_particle=integer_fields,
        state_kind=state_kind,
        physical_mode=physical_mode,
        cr_light_speed=cr_light_speed,
        model_ints=tuple(model_ints),
        model_reals=tuple(model_reals),
        particle_count=particle_count,
        particle_real_offset=particle_real_offset,
        particle_integer_offset=particle_integer_offset,
        payload_end_offset=payload_end_offset,
    )


def _parameter_blocks(payload: bytes, source: str) -> dict[str, dict[str, str]]:
    _require(type(payload) is bytes, f"{source}: restart payload must be bytes")
    parameter_end = payload.find(b"<par_end>")
    _require(parameter_end >= 0, f"{source}: restart header is missing <par_end>")
    try:
        text = payload[:parameter_end].decode("ascii")
    except UnicodeDecodeError as exc:
        raise RestartPolicyError(f"{source}: restart parameter header is not ASCII") from exc
    active_block: str | None = None
    blocks: dict[str, dict[str, str]] = {}
    for lineno, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<"):
            _require(
                line.endswith(">") and line.count("<") == 1 and line.count(">") == 1,
                f"{source}: malformed parameter block at line {lineno}",
            )
            active_block = line[1:-1]
            _require(active_block, f"{source}: empty parameter block at line {lineno}")
            blocks.setdefault(active_block, {})
            continue
        _require(active_block is not None, f"{source}: parameter outside a block at line {lineno}")
        _require(
            line.count("=") == 1,
            f"{source}: malformed parameter at line {lineno}",
        )
        key, value = (item.strip() for item in line.split("=", 1))
        _require(key, f"{source}: empty parameter name at line {lineno}")
        parameters = blocks[active_block]
        _require(
            key not in parameters,
            f"{source}: duplicate {active_block} parameter: {key}",
        )
        parameters[key] = value
    return blocks


def _problem_parameters(payload: bytes, source: str) -> dict[str, str]:
    blocks = _parameter_blocks(payload, source)
    _require("problem" in blocks, f"{source}: restart header is missing problem block")
    return blocks["problem"]


def _parse_integer(value: str, label: str) -> int:
    _require(_INTEGER_PATTERN.fullmatch(value) is not None, f"{label}: expected integer")
    return int(value)


def _parse_boolean(value: str, label: str) -> bool:
    _require(
        value in {"0", "1", "false", "true"},
        f"{label}: expected 0, 1, false or true",
    )
    return value in {"1", "true"}


def _parse_real(value: str, label: str) -> float:
    _require(_REAL_PATTERN.fullmatch(value) is not None, f"{label}: expected finite real")
    parsed = float(value)
    _require(math.isfinite(parsed), f"{label}: expected finite real")
    return parsed


def _header_parameter(
    blocks: Mapping[str, Mapping[str, str]], block: str, key: str, source: str
) -> str:
    _require(block in blocks, f"{source}: restart header is missing {block} block")
    _require(key in blocks[block], f"{source}: restart header is missing {block}/{key}")
    return blocks[block][key]


def _header_integer(
    blocks: Mapping[str, Mapping[str, str]], block: str, key: str, source: str
) -> int:
    return _parse_integer(
        _header_parameter(blocks, block, key, source), f"{source}/{block}/{key}"
    )


def _header_boolean(
    blocks: Mapping[str, Mapping[str, str]], block: str, key: str, source: str
) -> bool:
    return _parse_boolean(
        _header_parameter(blocks, block, key, source), f"{source}/{block}/{key}"
    )


def _header_real(
    blocks: Mapping[str, Mapping[str, str]], block: str, key: str, source: str
) -> float:
    return _parse_real(
        _header_parameter(blocks, block, key, source), f"{source}/{block}/{key}"
    )


def _fnv1a64_update(value: int, payload: bytes) -> int:
    for byte in payload:
        value ^= byte
        value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return value


def parallel_shock_restart_control_fingerprint(controls: object) -> str:
    """Recompute the C++ Q011 continuation-control fingerprint."""
    mapping = _keys(
        controls,
        {name for name, _ in _RESTART_CONTROL_KINDS},
        "parallel-shock restart controls",
    )
    value = _fnv1a64_update(
        14695981039346656037,
        _Q011_RESTART_FINGERPRINT_SCHEMA.encode("ascii") + b"\0",
    )
    for name, kind in _RESTART_CONTROL_KINDS:
        control = mapping[name]
        if kind == "integer":
            _require(
                type(control) is int,
                f"parallel-shock restart controls/{name}: expected integer",
            )
            encoded = struct.pack("<i", control)
        else:
            _require(
                type(control) is float,
                f"parallel-shock restart controls/{name}: expected float",
            )
            _require(
                math.isfinite(control),
                f"parallel-shock restart controls/{name}: expected finite float",
            )
            encoded = struct.pack("<d", control)
        value = _fnv1a64_update(value, name.encode("ascii") + b"\0" + encoded)
    return f"v1:{value:016x}"


def _species_config_hash(
    blocks: Mapping[str, Mapping[str, str]], nspecies: int, source: str
) -> int:
    value = _fnv1a64_update(14695981039346656037, struct.pack("<i", nspecies))
    for species in range(nspecies):
        block = f"species{species}"
        value = _fnv1a64_update(
            value, struct.pack("<d", _header_real(blocks, block, "mass", source))
        )
        value = _fnv1a64_update(
            value, struct.pack("<d", _header_real(blocks, block, "charge", source))
        )
    return value


def _validate_raw_particle_model(model: object, label: str) -> dict[str, object]:
    mapping = _keys(model, _RAW_PARTICLE_MODEL_KEYS, label)
    _require(mapping["state_kind"] == _Q011_STATE_KIND, f"{label}: Q011 state kind drift")
    _require(mapping["physical_mode"] == _Q011_PHYSICAL_MODE, f"{label}: Q011 physical mode drift")
    light_speed = mapping["cr_light_speed"]
    _require(
        type(light_speed) is float and math.isfinite(light_speed) and light_speed > 0.0,
        f"{label}: invalid particle light speed",
    )
    model_ints = mapping["model_ints"]
    model_reals = mapping["model_reals"]
    _require(
        type(model_ints) is list
        and len(model_ints) == _MODEL_INTEGER_COUNT
        and all(type(value) is int for value in model_ints),
        f"{label}: invalid PIC model integer binding",
    )
    _require(
        type(model_reals) is list
        and len(model_reals) == _MODEL_REAL_COUNT
        and all(type(value) is float and math.isfinite(value) for value in model_reals),
        f"{label}: invalid PIC model real binding",
    )
    return mapping


def _raw_q011_particle_model(
    probe: RestartPayloadProbe,
    blocks: Mapping[str, Mapping[str, str]],
    source: str,
) -> dict[str, object]:
    _require(
        probe.real_fields_per_particle == _Q011_REAL_FIELDS_PER_PARTICLE
        and probe.integer_fields_per_particle == _Q011_INTEGER_FIELDS_PER_PARTICLE,
        f"{source}: Q011 particle payload layout drift",
    )
    _require(probe.state_kind == _Q011_STATE_KIND, f"{source}: Q011 requires momentum state")
    _require(
        probe.physical_mode == _Q011_PHYSICAL_MODE,
        f"{source}: Q011 requires paper_mhd_pic_vl2_tsc physical mode",
    )
    _require(
        _header_parameter(blocks, "particles", "pic_physical_mode", source)
        == "paper_mhd_pic_vl2_tsc",
        f"{source}: particles/pic_physical_mode drift",
    )
    _require(
        _header_parameter(blocks, "particles", "pic_cr_initial_state", source)
        == "momentum",
        f"{source}: particles/pic_cr_initial_state drift",
    )
    _strict_equal(
        _header_real(blocks, "particles", "pic_cr_light_speed", source),
        probe.cr_light_speed,
        f"{source}/raw particle light speed",
    )

    ints = list(probe.model_ints)
    reals = list(probe.model_reals)
    nspecies = _header_integer(blocks, "particles", "nspecies", source)
    _require(nspecies > 0, f"{source}: invalid particles/nspecies")
    species_hash = _species_config_hash(blocks, nspecies, source)
    def exact_particle_choice(key: str, expected: str, encoded: int) -> int:
        _strict_equal(
            _header_parameter(blocks, "particles", key, source),
            expected,
            f"{source}/particles/{key}",
        )
        return encoded

    expected_ints = [
        exact_particle_choice("pic_deltaf_mode", "off", 0),
        0,
        exact_particle_choice("pic_expanding_box_mode", "off", 0),
        exact_particle_choice("pic_expansion_law", "linear", 0),
        exact_particle_choice("pic_wave_damping_mode", "off", 0),
        exact_particle_choice("pic_deltaf_adapt_mode", "off", 0),
        exact_particle_choice("particle_type", "cosmic_ray", 0),
        exact_particle_choice("pusher", "boris_tsc", 6),
        nspecies,
        int(_header_boolean(blocks, "particles", "track_displacement", source)),
        int(_header_boolean(blocks, "particles", "deposit_moments", source)),
        _header_integer(blocks, "particles", "deposit_order", source),
        int(_header_boolean(blocks, "particles", "couple_moments_to_mhd", source)),
        exact_particle_choice("couple_j_to_efield_representation", "cell_centered", 0),
        exact_particle_choice("couple_j_deposition_mode", "cc_convert", 0),
        exact_particle_choice("couple_fluid_feedback_order", "mhd_src_terms", 0),
        int(
            _header_boolean(
                blocks, "particles", "couple_moments_momentum_to_mhd", source
            )
        ),
        int(
            _header_boolean(blocks, "particles", "couple_moments_energy_to_mhd", source)
        ),
        exact_particle_choice("pic_background_mode", "coupled", 0),
        exact_particle_choice("pic_feedback_mode", "coupled", 0),
        exact_particle_choice("pic_cr_hall_mode", "off", 0),
        exact_particle_choice("pic_cr_initial_state", "momentum", 1),
        exact_particle_choice("pic_interp_scheme", "tsc", 0),
        int(_header_boolean(blocks, "particles", "pic_enable_2d3v", source)),
        exact_particle_choice("pic_intermediate_arrays", "auto", 0),
        _header_integer(blocks, "particles", "pic_max_cell_cross", source),
        _header_integer(blocks, "particles", "pic_sort_interval", source),
        _header_integer(blocks, "particles", "pic_random_seed", source),
        species_hash & 0x3FFFFF,
        (species_hash >> 22) & 0x3FFFFF,
        (species_hash >> 44) & 0xFFFFF,
    ]
    for index, expected in enumerate(expected_ints):
        _require(ints[index] == expected, f"{source}: PIC model integer {index} drift")
    expected_reals = [
        _header_real(blocks, "particles", "pic_expansion_rate_x1", source),
        _header_real(blocks, "particles", "pic_expansion_rate_x2", source),
        _header_real(blocks, "particles", "pic_expansion_rate_x3", source),
        _header_real(blocks, "particles", "pic_deltaf_p0", source),
        _header_real(blocks, "particles", "pic_deltaf_kappa", source),
        _header_real(blocks, "particles", "pic_deltaf_drift_x1", source),
        _header_real(blocks, "particles", "pic_deltaf_drift_x2", source),
        _header_real(blocks, "particles", "pic_deltaf_drift_x3", source),
        _header_real(blocks, "particles", "pic_deltaf_aniso_x1", source),
        _header_real(blocks, "particles", "pic_deltaf_aniso_x2", source),
        _header_real(blocks, "particles", "pic_deltaf_aniso_x3", source),
        _header_real(blocks, "particles", "pic_deltaf_background_rho", source),
        _header_real(blocks, "particles", "pic_deltaf_background_jx", source),
        _header_real(blocks, "particles", "pic_deltaf_background_jy", source),
        _header_real(blocks, "particles", "pic_deltaf_background_jz", source),
        _header_real(blocks, "particles", "pic_no_mhd_bx", source),
        _header_real(blocks, "particles", "pic_no_mhd_by", source),
        _header_real(blocks, "particles", "pic_no_mhd_bz", source),
        _header_real(blocks, "particles", "pic_ion_neutral_collision_rate", source),
        _header_real(blocks, "particles", "pic_deltaf_adapt_interval", source),
        _header_real(blocks, "particles", "deposit_qscale", source),
        _header_real(blocks, "particles", "couple_j_to_efield_coeff", source),
        _header_real(
            blocks, "particles", "couple_moments_momentum_coeff", source
        ),
        _header_real(blocks, "particles", "couple_moments_energy_coeff", source),
        _header_real(blocks, "particles", "pic_theta_max", source),
        _header_real(
            blocks, "particles", "pic_load_balance_cost_per_particle", source
        ),
    ]
    for index, expected in enumerate(expected_reals):
        _strict_equal(reals[index], expected, f"{source}/PIC model real {index}")
    return dict(
        _validate_raw_particle_model(
            {
                "state_kind": probe.state_kind,
                "physical_mode": probe.physical_mode,
                "cr_light_speed": probe.cr_light_speed,
                "model_ints": ints,
                "model_reals": reals,
            },
            f"{source}/raw_particle_model",
        )
    )


def _parallel_shock_restart_controls(
    blocks: Mapping[str, Mapping[str, str]],
    probe: RestartPayloadProbe,
    source: str,
) -> dict[str, int | float]:
    _require(
        _header_parameter(blocks, "problem", "pgen_name", source) == "pic_parallel_shock",
        f"{source}: problem/pgen_name drift",
    )
    problem_real = lambda key: _header_real(blocks, "problem", key, source)
    problem_int = lambda key: _header_integer(blocks, "problem", key, source)
    problem_bool = lambda key: int(_header_boolean(blocks, "problem", key, source))
    shock_model_name = _header_parameter(blocks, "problem", "ps_shock_speed_model", source)
    _require(
        shock_model_name in {"finite_mach", "ideal_surface"},
        f"{source}: shock-speed model drift",
    )
    shock_model = 0 if shock_model_name == "finite_mach" else 1
    frame_mode_name = _header_parameter(blocks, "problem", "ps_frame_mode", source)
    _require(frame_mode_name in {"velocity", "recenter"}, f"{source}: frame mode drift")
    frame_mode = 0 if frame_mode_name == "velocity" else 1
    inject_species = problem_int("ps_inject_species")
    species_block = f"species{inject_species}"
    particle_mass = _header_real(blocks, species_block, "mass", source)
    particle_charge = _header_real(blocks, species_block, "charge", source)
    deposit_qscale = _header_real(blocks, "particles", "deposit_qscale", source)
    rho0 = problem_real("ps_rho0")
    p0 = problem_real("ps_p0")
    u0 = problem_real("ps_u0")
    gamma = _header_real(blocks, "mhd", "gamma", source)
    if shock_model == 1:
        shock_speed = 0.5 * (gamma - 1.0) * u0
    else:
        cs2 = gamma * p0 / rho0
        ms2 = u0 * u0 / cs2
        compression = ((gamma + 1.0) * ms2) / ((gamma - 1.0) * ms2 + 2.0)
        shock_speed = u0 / (compression - 1.0)
    x1min = _header_real(blocks, "mesh", "x1min", source)
    x1max = _header_real(blocks, "mesh", "x1max", source)
    nx1 = _header_integer(blocks, "mesh", "nx1", source)
    nx2 = _header_integer(blocks, "mesh", "nx2", source)
    nx3 = _header_integer(blocks, "mesh", "nx3", source)
    _require(nx1 > 0 and nx2 > 0 and nx3 > 0 and x1max > x1min, f"{source}: mesh controls drift")
    controls: dict[str, int | float] = {
        "ps_rho0": rho0,
        "ps_p0": p0,
        "ps_u0": u0,
        "ps_b0": problem_real("ps_b0"),
        "ps_eta": problem_real("ps_eta"),
        "ps_vinj_over_u0": problem_real("ps_vinj_over_u0"),
        "ps_inject_half_width_cells": problem_real("ps_inject_half_width_cells"),
        "ps_inject_t_start": problem_real("ps_inject_t_start"),
        "ps_inject_t_stop": problem_real("ps_inject_t_stop"),
        "ps_remove_birth_time_before": problem_real("ps_remove_birth_time_before"),
        "ps_shock_speed_model": shock_model,
        "ps_refine_curv": problem_real("ps_refine_curv"),
        "ps_derefine_curv": problem_real("ps_derefine_curv"),
        "ps_rho_floor_frac": problem_real("ps_rho_floor_frac"),
        "ps_p_floor_frac": problem_real("ps_p_floor_frac"),
        "ps_enable_injection": problem_bool("ps_enable_injection"),
        "ps_enable_subtraction": problem_bool("ps_enable_gas_subtraction"),
        "ps_enable_curvature_amr": problem_bool("ps_enable_curvature_amr"),
        "ps_test_source_transaction_terms_override": problem_bool(
            "ps_test_source_transaction_terms_override"
        ),
        "ps_test_source_transaction_terms": problem_real(
            "ps_test_source_transaction_terms"
        ),
        "ps_inject_species": inject_species,
        "ps_inject_seed": problem_int("ps_inject_seed"),
        "ps_particle_mass": particle_mass,
        "ps_particle_charge": particle_charge,
        "ps_particle_q_over_m": particle_charge / particle_mass,
        "ps_particle_macro_mass": deposit_qscale * particle_mass,
        "ps_particle_momentum_state": probe.state_kind,
        "ps_particle_light_speed": probe.cr_light_speed,
        "ps_enable_frame_tracking": problem_bool("ps_enable_frame_tracking"),
        "ps_frame_mode": frame_mode,
        "ps_frame_t_start": problem_real("ps_frame_t_start"),
        "ps_frame_t_ramp": problem_real("ps_frame_t_ramp"),
        "ps_frame_vfrac": problem_real("ps_frame_vfrac"),
        "ps_frame_dv_max": problem_real("ps_frame_dv_max"),
        "ps_frame_apply_to_particles": problem_bool("ps_frame_apply_to_particles"),
        "ps_frame_apply_to_inflow": problem_bool("ps_frame_apply_to_inflow"),
        "ps_frame_require_uniform": problem_bool("ps_frame_require_uniform"),
        "ps_recenter_x_target": problem_real("ps_recenter_x_target"),
        "ps_recenter_x_trigger": problem_real("ps_recenter_x_trigger"),
        "ps_recenter_dx1": (x1max - x1min) / nx1,
        "ps_recenter_shift_cells": problem_int("ps_recenter_shift_cells"),
        "ps_recenter_vshock_model": problem_real("ps_recenter_vshock_model"),
        "ps_shock_speed": shock_speed,
        "ps_xshock0": x1min,
        "ps_use_2d3v": int(
            nx2 > 1
            and nx3 == 1
            and _header_boolean(blocks, "particles", "pic_enable_2d3v", source)
        ),
    }
    parallel_shock_restart_control_fingerprint(controls)
    return controls


def _validate_restart_control_binding(binding: object, label: str) -> dict[str, object]:
    mapping = _keys(binding, _RESTART_CONTROL_BINDING_KEYS, label)
    _strict_equal(
        mapping["fingerprint_schema"],
        _Q011_RESTART_FINGERPRINT_SCHEMA,
        f"{label}/fingerprint_schema",
    )
    stored = mapping["stored_fingerprint"]
    _require(
        type(stored) is str and _RESTART_FINGERPRINT_PATTERN.fullmatch(stored) is not None,
        f"{label}: invalid stored restart-control fingerprint",
    )
    _strict_equal(
        parallel_shock_restart_control_fingerprint(mapping["controls"]),
        stored,
        f"{label}/recomputed restart-control fingerprint",
    )
    return mapping


def _restart_control_binding(
    blocks: Mapping[str, Mapping[str, str]],
    probe: RestartPayloadProbe,
    source: str,
) -> dict[str, object]:
    controls = _parallel_shock_restart_controls(blocks, probe, source)
    stored = _header_parameter(
        blocks, "problem", "ps_restart_control_fingerprint", source
    )
    binding = {
        "fingerprint_schema": _Q011_RESTART_FINGERPRINT_SCHEMA,
        "stored_fingerprint": stored,
        "controls": controls,
    }
    return dict(
        _validate_restart_control_binding(binding, f"{source}/restart_control_binding")
    )


def _summation_tolerance(
    lhs: float, rhs: float, absolute_contributions: float, accumulated_terms: float
) -> float:
    _require(
        all(
            math.isfinite(value)
            for value in (lhs, rhs, absolute_contributions, accumulated_terms)
        )
        and absolute_contributions >= 0.0
        and accumulated_terms >= 0.0,
        "invalid cancellation-aware summation metadata",
    )
    relative_bound = max(accumulated_terms, 1.0) * _BINARY64_EPSILON
    _require(relative_bound < 0.5, "unreliable cancellation-aware summation metadata")
    summation_bound = relative_bound / (1.0 - relative_bound)
    return (8.0 * summation_bound + 64.0 * _BINARY64_EPSILON) * max(
        abs(lhs), abs(rhs), absolute_contributions, 1.0
    )


def _validate_q011_population_energy(
    ledger: Mapping[str, object], prefix: str, label: str, light_speed: float
) -> None:
    count = ledger[prefix + "_count_global"]
    mass = ledger[prefix + "_mass_global"]
    momentum = [
        ledger[prefix + "_momentum_x1_global"],
        ledger[prefix + "_momentum_x2_global"],
        ledger[prefix + "_momentum_x3_global"],
    ]
    energy = ledger[prefix + "_energy_global"]
    if count == 0.0:
        _require(
            mass == 0.0
            and all(component == 0.0 for component in momentum)
            and energy == 0.0,
            f"{label}: empty {prefix} population contains accumulated state",
        )
        return
    _require(mass > 0.0, f"{label}: nonempty {prefix} population has nonpositive mass")
    _require(
        type(light_speed) is float and math.isfinite(light_speed) and light_speed > 0.0,
        f"{label}: invalid particle light speed",
    )
    state_components = [component / mass for component in momentum]
    state_squared = sum(component * component for component in state_components)
    _require(
        math.isfinite(state_squared),
        f"{label}: {prefix} aggregate momentum is not finite",
    )
    specific_lower_bound = state_squared / (
        math.sqrt(1.0 + state_squared / (light_speed * light_speed)) + 1.0
    )
    lower_bound = mass * specific_lower_bound
    relative_bound = max(count, 1.0) * _BINARY64_EPSILON
    _require(
        math.isfinite(lower_bound) and relative_bound < 0.5,
        f"{label}: {prefix} energy-momentum admissibility bound is invalid",
    )
    summation_bound = relative_bound / (1.0 - relative_bound)
    tolerance = (8.0 * summation_bound + 64.0 * _BINARY64_EPSILON) * max(
        abs(energy), abs(lower_bound), abs(mass), 1.0
    )
    _require(
        energy + tolerance >= lower_bound,
        f"{label}: {prefix} energy is below its aggregate momentum lower bound",
    )


def _validate_startup_shock_ledger(
    ledger: object, label: str, *, light_speed: float
) -> dict[str, object]:
    ledger_mapping = _keys(ledger, set(STARTUP_SHOCK_LEDGER_FIELDS), label)
    for field in _INTEGER_LEDGER_FIELDS:
        _require(type(ledger_mapping[field]) is int, f"{label}/{field}: expected integer")
    for field in _BOOLEAN_LEDGER_FIELDS:
        _require(type(ledger_mapping[field]) is bool, f"{label}/{field}: expected boolean")
    for field in _REAL_LEDGER_FIELDS:
        value = ledger_mapping[field]
        _require(type(value) is float, f"{label}/{field}: expected float")
        _require(math.isfinite(value), f"{label}/{field}: expected finite float")

    _require(
        ledger_mapping["ps_cr_ledger_schema"] == EXPECTED_SHOCK_LEDGER_SCHEMA,
        f"{label}: shock ledger schema is not 3",
    )
    _require(ledger_mapping["ps_cr_ledger_complete"], f"{label}: shock ledger is incomplete")
    _require(
        ledger_mapping["ps_removed_excluded_early_cohort"],
        f"{label}: startup cohort removal is incomplete",
    )
    _require(ledger_mapping["ps_tag_seeded"], f"{label}: shock particle tags are unseeded")
    for field in _NONNEGATIVE_REAL_LEDGER_FIELDS:
        _require(ledger_mapping[field] >= 0.0, f"{label}/{field}: negative ledger value")
    for field in _COUNT_LEDGER_FIELDS:
        _require(
            ledger_mapping[field].is_integer(),
            f"{label}/{field}: non-integral particle count",
        )
    _require(
        ledger_mapping["ps_injected_cr_count_global"]
        >= ledger_mapping["ps_removed_cr_count_global"],
        f"{label}: removed count exceeds injected count",
    )
    _require(
        ledger_mapping["ps_injected_cr_mass_global"]
        >= ledger_mapping["ps_removed_cr_mass_global"],
        f"{label}: removed mass exceeds injected mass",
    )
    _require(
        ledger_mapping["ps_injected_cr_energy_global"]
        >= ledger_mapping["ps_removed_cr_energy_global"],
        f"{label}: removed energy exceeds injected energy",
    )
    for prefix in ("ps_injected_cr", "ps_removed_cr"):
        if ledger_mapping[prefix + "_count_global"] == 0.0:
            for suffix in (
                "_momentum_x1_global",
                "_momentum_x2_global",
                "_momentum_x3_global",
                "_energy_global",
            ):
                _require(
                    ledger_mapping[prefix + suffix] == 0.0,
                    f"{label}: empty {prefix} ledger contains accumulated state",
                )
        _validate_q011_population_energy(ledger_mapping, prefix, label, light_speed)
    tag_floor = ledger_mapping["ps_injection_tag_floor"]
    injected_count = int(ledger_mapping["ps_injected_cr_count_global"])
    next_tag = ledger_mapping["ps_next_tag"]
    _require(0 <= tag_floor <= _MAX_SIGNED_INT, f"{label}: invalid injection tag floor")
    _require(
        injected_count <= _MAX_SIGNED_INT - tag_floor
        and next_tag == tag_floor + injected_count,
        f"{label}: invalid next-tag progression",
    )
    return ledger_mapping


def extract_startup_shock_ledger(
    payload: bytes, *, source: str = "<restart-bytes>"
) -> dict[str, object]:
    """Extract the complete schema-3 startup shock ledger from restart bytes."""
    probe = probe_schema7_restart_payload(payload, source=source)
    parameters = _problem_parameters(payload, source)
    missing = [field for field in STARTUP_SHOCK_LEDGER_FIELDS if field not in parameters]
    _require(not missing, f"{source}: startup shock ledger is missing entries {missing!r}")
    ledger: dict[str, object] = {}
    for field in STARTUP_SHOCK_LEDGER_FIELDS:
        label = f"{source}/{field}"
        if field in _INTEGER_LEDGER_FIELDS:
            ledger[field] = _parse_integer(parameters[field], label)
        elif field in _BOOLEAN_LEDGER_FIELDS:
            ledger[field] = _parse_boolean(parameters[field], label)
        else:
            ledger[field] = _parse_real(parameters[field], label)
    return dict(
        _validate_startup_shock_ledger(
            ledger,
            f"{source}/startup_shock_ledger",
            light_speed=probe.cr_light_speed,
        )
    )


def _validate_particle_escape_ledger(
    ledger: object,
    label: str,
    *,
    startup_shock_ledger: object | None = None,
    committed_cycle: object | None = None,
    committed_time: object | None = None,
    light_speed: float,
) -> dict[str, object]:
    ledger_mapping = _keys(ledger, set(ESCAPE_LEDGER_FIELDS), label)
    for field in _ESCAPE_INTEGER_LEDGER_FIELDS:
        _require(type(ledger_mapping[field]) is int, f"{label}/{field}: expected integer")
    for field in _ESCAPE_BOOLEAN_LEDGER_FIELDS:
        _require(type(ledger_mapping[field]) is bool, f"{label}/{field}: expected boolean")
    for field in _ESCAPE_REAL_LEDGER_FIELDS:
        value = ledger_mapping[field]
        _require(type(value) is float, f"{label}/{field}: expected float")
        _require(math.isfinite(value), f"{label}/{field}: expected finite float")

    _require(
        ledger_mapping["ps_escape_ledger_schema"] == EXPECTED_ESCAPE_LEDGER_SCHEMA,
        f"{label}: escape ledger schema is not 2",
    )
    _require(
        ledger_mapping["ps_escape_ledger_complete"],
        f"{label}: escape ledger is incomplete",
    )
    _require(
        ledger_mapping["ps_escape_audit_calls"] >= 0,
        f"{label}: negative escape-audit call count",
    )
    for field in (
        "ps_escape_last_audit_time",
        "ps_escaped_injected_cr_count_global",
        "ps_escaped_injected_cr_mass_global",
        "ps_escaped_injected_cr_energy_global",
        "ps_escaped_initial_cr_count_global",
        "ps_escaped_injected_cr_term_count_global",
        "ps_escaped_injected_cr_abs_mass_global",
        "ps_escaped_injected_cr_abs_momentum_x1_global",
        "ps_escaped_injected_cr_abs_momentum_x2_global",
        "ps_escaped_injected_cr_abs_momentum_x3_global",
        "ps_escaped_injected_cr_abs_energy_global",
    ):
        _require(ledger_mapping[field] >= 0.0, f"{label}/{field}: negative ledger value")
    for field in (
        "ps_escaped_injected_cr_count_global",
        "ps_escaped_initial_cr_count_global",
        "ps_escaped_injected_cr_term_count_global",
    ):
        _require(
            ledger_mapping[field].is_integer(),
            f"{label}/{field}: non-integral particle count",
        )
    _require(
        ledger_mapping["ps_escaped_initial_cr_count_global"] == 0.0,
        f"{label}: initial-particle physical escape is unaccounted",
    )
    if ledger_mapping["ps_escaped_injected_cr_count_global"] == 0.0:
        for field in (
            "ps_escaped_injected_cr_momentum_x1_global",
            "ps_escaped_injected_cr_momentum_x2_global",
            "ps_escaped_injected_cr_momentum_x3_global",
            "ps_escaped_injected_cr_energy_global",
            "ps_escaped_injected_cr_term_count_global",
            "ps_escaped_injected_cr_abs_mass_global",
            "ps_escaped_injected_cr_abs_momentum_x1_global",
            "ps_escaped_injected_cr_abs_momentum_x2_global",
            "ps_escaped_injected_cr_abs_momentum_x3_global",
            "ps_escaped_injected_cr_abs_energy_global",
        ):
            _require(
                ledger_mapping[field] == 0.0,
                f"{label}: empty escape ledger contains accumulated state",
            )
    _validate_q011_population_energy(
        ledger_mapping, "ps_escaped_injected_cr", label, light_speed
    )
    term_count = ledger_mapping["ps_escaped_injected_cr_term_count_global"]
    escaped_count = ledger_mapping["ps_escaped_injected_cr_count_global"]
    _require(
        term_count == escaped_count,
        f"{label}: escape comparison term count differs from escaped count",
    )
    for net_field, absolute_field in (
        (
            "ps_escaped_injected_cr_mass_global",
            "ps_escaped_injected_cr_abs_mass_global",
        ),
        (
            "ps_escaped_injected_cr_energy_global",
            "ps_escaped_injected_cr_abs_energy_global",
        ),
    ):
        net = ledger_mapping[net_field]
        absolute = ledger_mapping[absolute_field]
        _require(
            abs(net - absolute)
            <= _summation_tolerance(net, absolute, absolute, term_count),
            f"{label}: {absolute_field} is inconsistent with positive contributions",
        )
    for net_field, absolute_field in (
        (
            "ps_escaped_injected_cr_momentum_x1_global",
            "ps_escaped_injected_cr_abs_momentum_x1_global",
        ),
        (
            "ps_escaped_injected_cr_momentum_x2_global",
            "ps_escaped_injected_cr_abs_momentum_x2_global",
        ),
        (
            "ps_escaped_injected_cr_momentum_x3_global",
            "ps_escaped_injected_cr_abs_momentum_x3_global",
        ),
    ):
        net = ledger_mapping[net_field]
        absolute = ledger_mapping[absolute_field]
        _require(
            abs(net)
            <= absolute + _summation_tolerance(net, absolute, absolute, term_count),
            f"{label}: {absolute_field} is below the net escaped momentum",
        )
    if ledger_mapping["ps_escape_audit_calls"] == 0:
        for field in _ESCAPE_REAL_LEDGER_FIELDS:
            _require(
                ledger_mapping[field] == 0.0,
                f"{label}: zero-call escape ledger contains accumulated state",
            )

    if startup_shock_ledger is not None:
        startup = _validate_startup_shock_ledger(
            startup_shock_ledger,
            f"{label}/startup_shock_ledger",
            light_speed=light_speed,
        )
        escaped_mass = ledger_mapping["ps_escaped_injected_cr_mass_global"]
        injected_count = startup["ps_injected_cr_count_global"]
        injected_mass = startup["ps_injected_cr_mass_global"]
        removed_count = startup["ps_removed_cr_count_global"]
        removed_mass = startup["ps_removed_cr_mass_global"]
        _require(
            escaped_count + removed_count <= injected_count,
            f"{label}: removed plus escaped count exceeds injected count",
        )
        _require(
            escaped_mass + removed_mass <= injected_mass,
            f"{label}: removed plus escaped mass exceeds injected mass",
        )
        if injected_count == 0.0:
            _require(
                escaped_count == 0.0 and escaped_mass == 0.0,
                f"{label}: escaped population exists without injected population",
            )
        else:
            expected_escaped_mass = escaped_count * injected_mass / injected_count
            _require(
                math.isclose(
                    escaped_mass,
                    expected_escaped_mass,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-12,
                ),
                f"{label}: escaped mass is inconsistent with injected macro-mass",
            )

    if committed_cycle is not None or committed_time is not None:
        cycle = _canonical_observed_cycle(committed_cycle, f"{label}/committed_cycle")
        time = _canonical_observed_time(committed_time, f"{label}/committed_time")
        expected_calls = 2 * cycle
        _require(
            expected_calls <= _MAX_SIGNED_INT
            and ledger_mapping["ps_escape_audit_calls"] == expected_calls,
            f"{label}: paper-VL2 escape-audit call chronology is inconsistent",
        )
        expected_last_time = 0.0 if expected_calls == 0 else time
        _strict_equal(
            ledger_mapping["ps_escape_last_audit_time"],
            expected_last_time,
            f"{label}/paper-VL2 last escape-audit time",
        )
    return ledger_mapping


def extract_particle_escape_ledger(
    payload: bytes, *, source: str = "<restart-bytes>"
) -> dict[str, object]:
    """Extract the complete schema-2 Q-011 physical-boundary escape ledger."""
    probe = probe_schema7_restart_payload(payload, source=source)
    parameters = _problem_parameters(payload, source)
    missing = [field for field in ESCAPE_LEDGER_FIELDS if field not in parameters]
    _require(not missing, f"{source}: particle escape ledger is missing entries {missing!r}")
    ledger: dict[str, object] = {}
    for field in ESCAPE_LEDGER_FIELDS:
        label = f"{source}/{field}"
        if field in _ESCAPE_INTEGER_LEDGER_FIELDS:
            ledger[field] = _parse_integer(parameters[field], label)
        elif field in _ESCAPE_BOOLEAN_LEDGER_FIELDS:
            ledger[field] = _parse_boolean(parameters[field], label)
        else:
            ledger[field] = _parse_real(parameters[field], label)
    return dict(
        _validate_particle_escape_ledger(
            ledger,
            f"{source}/particle_escape_ledger",
            light_speed=probe.cr_light_speed,
        )
    )


def _validate_active_particle_cohort(
    cohort: object,
    label: str,
    *,
    startup_shock_ledger: object,
    particle_escape_ledger: object,
    light_speed: float,
) -> dict[str, object]:
    mapping = _keys(cohort, _ACTIVE_PARTICLE_COHORT_KEYS, label)
    for field in _ACTIVE_PARTICLE_COHORT_KEYS:
        _require(type(mapping[field]) is int, f"{label}/{field}: expected integer")
    _require(mapping["particle_count"] > 0, f"{label}: particle cohort is empty")
    _require(
        mapping["initial_count"] >= 0 and mapping["shock_injected_count"] >= 0,
        f"{label}: negative particle-source count",
    )
    _require(
        mapping["initial_count"] + mapping["shock_injected_count"]
        == mapping["particle_count"],
        f"{label}: particle-source counts do not cover the active cohort",
    )
    _require(
        mapping["initial_count"] == 0,
        f"{label}: active initial particle remains after startup-cohort removal",
    )
    _require(
        0 <= mapping["minimum_tag"] <= mapping["maximum_tag"] <= _MAX_SIGNED_INT,
        f"{label}: invalid active particle tag range",
    )
    startup = _validate_startup_shock_ledger(
        startup_shock_ledger, f"{label}/startup_shock_ledger", light_speed=light_speed
    )
    escape = _validate_particle_escape_ledger(
        particle_escape_ledger,
        f"{label}/particle_escape_ledger",
        startup_shock_ledger=startup,
        light_speed=light_speed,
    )
    _require(
        mapping["shock_injected_count"]
        + int(startup["ps_removed_cr_count_global"])
        + int(escape["ps_escaped_injected_cr_count_global"])
        == int(startup["ps_injected_cr_count_global"]),
        f"{label}: active plus removed plus escaped injected count is inconsistent",
    )
    return mapping


def _active_particle_cohort_from_raw(
    payload: bytes,
    probe: RestartPayloadProbe,
    restart_control_binding: Mapping[str, object],
    startup_shock_ledger: Mapping[str, object],
    particle_escape_ledger: Mapping[str, object],
    committed_time: float,
    source: str,
) -> dict[str, object]:
    controls = _validate_restart_control_binding(
        restart_control_binding, f"{source}/restart_control_binding"
    )["controls"]
    real_values = np.frombuffer(
        payload,
        dtype="<f8",
        count=probe.particle_count * probe.real_fields_per_particle,
        offset=probe.particle_real_offset,
    ).reshape(probe.particle_count, probe.real_fields_per_particle)
    integer_values = np.frombuffer(
        payload,
        dtype="<i4",
        count=probe.particle_count * probe.integer_fields_per_particle,
        offset=probe.particle_integer_offset,
    ).reshape(probe.particle_count, probe.integer_fields_per_particle)
    _require(np.all(np.isfinite(real_values)), f"{source}: active particle payload is nonfinite")
    gids = integer_values[:, _Q011_PGID]
    tags = integer_values[:, _Q011_PTAG]
    species = integer_values[:, _Q011_PSP]
    sources = integer_values[:, _Q011_PCRSOURCE]
    nspecies = probe.model_ints[8]
    _require(np.all(gids >= 0), f"{source}: active particle gid is negative")
    _require(np.all(tags >= 0), f"{source}: active particle tag is negative")
    _require(
        np.unique(tags).size == probe.particle_count,
        f"{source}: active particle tags are not unique",
    )
    _require(
        np.all((species >= 0) & (species < nspecies)),
        f"{source}: active particle species is out of range",
    )
    _require(
        np.all(
            (sources == _Q011_SOURCE_INITIAL)
            | (sources == _Q011_SOURCE_SHOCK_INJECTED)
        ),
        f"{source}: active particle source is invalid",
    )
    initial = sources == _Q011_SOURCE_INITIAL
    injected = sources == _Q011_SOURCE_SHOCK_INJECTED
    tag_floor = int(startup_shock_ledger["ps_injection_tag_floor"])
    next_tag = int(startup_shock_ledger["ps_next_tag"])
    _require(
        np.all(tags[initial] < tag_floor),
        f"{source}: active initial particle overlaps the injected tag range",
    )
    _require(
        np.all((tags[injected] >= tag_floor) & (tags[injected] < next_tag)),
        f"{source}: active injected particle tag is outside the injected tag range",
    )
    _require(
        np.all(species[injected] == controls["ps_inject_species"]),
        f"{source}: active injected particle species drift",
    )
    _require(
        np.all(real_values[injected, _Q011_IPM] == controls["ps_particle_q_over_m"]),
        f"{source}: active injected particle q/m drift",
    )
    _require(
        np.all(real_values[injected, _Q011_IPWT] == 1.0),
        f"{source}: active injected particle macro-weight drift",
    )
    births = real_values[injected, _Q011_IPT_BIRTH]
    _require(
        np.all(
            (births >= controls["ps_inject_t_start"])
            & (births <= controls["ps_inject_t_stop"])
            & (births <= committed_time)
            & (births >= controls["ps_remove_birth_time_before"])
        ),
        f"{source}: active injected particle birth-time cohort is inconsistent",
    )
    cohort = {
        "particle_count": int(probe.particle_count),
        "initial_count": int(np.count_nonzero(initial)),
        "shock_injected_count": int(np.count_nonzero(injected)),
        "minimum_tag": int(np.min(tags)),
        "maximum_tag": int(np.max(tags)),
    }
    return dict(
        _validate_active_particle_cohort(
            cohort,
            f"{source}/active_particle_cohort",
            startup_shock_ledger=startup_shock_ledger,
            particle_escape_ledger=particle_escape_ledger,
            light_speed=probe.cr_light_speed,
        )
    )


def validate_checkpoint_binding(
    binding: object,
    *,
    preregistration: Mapping[str, object] | None = None,
) -> None:
    """Validate one frozen checkpoint binding before any execution."""
    policy = _validated_preregistration(preregistration)
    binding_mapping = _keys(binding, _BINDING_KEYS, "checkpoint binding")
    contract = policy["continuation_contract"]
    _strict_equal(
        binding_mapping["checkpoint_nominal_slot_omega0_inverse"],
        contract["checkpoint_nominal_slot_omega0_inverse"],
        "checkpoint binding/checkpoint_nominal_slot_omega0_inverse",
    )
    _canonical_observed_cycle(
        binding_mapping["checkpoint_observed_committed_cycle"],
        "checkpoint binding/checkpoint_observed_committed_cycle",
    )
    _canonical_observed_time(
        binding_mapping["checkpoint_observed_committed_time_omega0_inverse"],
        "checkpoint binding/checkpoint_observed_committed_time_omega0_inverse",
    )
    _strict_equal(
        binding_mapping["restart_schema"],
        policy["restart_payload_probe"]["restart_schema"],
        "checkpoint binding/restart_schema",
    )
    raw_model = _validate_raw_particle_model(
        binding_mapping["raw_particle_model"], "checkpoint binding/raw_particle_model"
    )
    restart_controls = _validate_restart_control_binding(
        binding_mapping["restart_control_binding"],
        "checkpoint binding/restart_control_binding",
    )
    _strict_equal(
        restart_controls["controls"]["ps_particle_light_speed"],
        raw_model["cr_light_speed"],
        "checkpoint binding/particle light-speed parity",
    )
    _validate_startup_shock_ledger(
        binding_mapping["startup_shock_ledger"],
        "checkpoint binding/startup_shock_ledger",
        light_speed=raw_model["cr_light_speed"],
    )
    _validate_particle_escape_ledger(
        binding_mapping["particle_escape_ledger"],
        "checkpoint binding/particle_escape_ledger",
        startup_shock_ledger=binding_mapping["startup_shock_ledger"],
        committed_cycle=binding_mapping["checkpoint_observed_committed_cycle"],
        committed_time=binding_mapping[
            "checkpoint_observed_committed_time_omega0_inverse"
        ],
        light_speed=raw_model["cr_light_speed"],
    )
    _validate_active_particle_cohort(
        binding_mapping["active_particle_cohort"],
        "checkpoint binding/active_particle_cohort",
        startup_shock_ledger=binding_mapping["startup_shock_ledger"],
        particle_escape_ledger=binding_mapping["particle_escape_ledger"],
        light_speed=raw_model["cr_light_speed"],
    )
    _strict_equal(
        binding_mapping["retained_output_nominal_slots_after_checkpoint_omega0_inverse"],
        contract["retained_output_nominal_slots_after_checkpoint_omega0_inverse"],
        "checkpoint binding/retained_output_nominal_slots_after_checkpoint_omega0_inverse",
    )
    _strict_equal(
        binding_mapping["comparison_tolerances_max_absolute_difference"],
        contract["comparison_tolerances_max_absolute_difference"],
        "checkpoint binding/comparison_tolerances_max_absolute_difference",
    )


def bind_checkpoint_for_continuation(
    restart_payload: bytes,
    *,
    checkpoint_nominal_slot_omega0_inverse: object,
    checkpoint_observed_committed_cycle: object,
    checkpoint_observed_committed_time_omega0_inverse: object,
    retained_output_nominal_slots_after_checkpoint_omega0_inverse: object,
    comparison_tolerances_max_absolute_difference: object,
    preregistration: Mapping[str, object] | None = None,
    source: str = "<restart-bytes>",
) -> dict[str, object]:
    """Probe retained bytes and build a validated pre-execution checkpoint binding."""
    policy = _validated_preregistration(preregistration)
    probe = probe_schema7_restart_payload(restart_payload, source=source)
    blocks = _parameter_blocks(restart_payload, source)
    raw_model = _raw_q011_particle_model(probe, blocks, source)
    restart_controls = _restart_control_binding(blocks, probe, source)
    ledger = extract_startup_shock_ledger(restart_payload, source=source)
    escape_ledger = extract_particle_escape_ledger(restart_payload, source=source)
    committed_time = _canonical_observed_time(
        checkpoint_observed_committed_time_omega0_inverse,
        "checkpoint observed committed time",
    )
    active_cohort = _active_particle_cohort_from_raw(
        restart_payload,
        probe,
        restart_controls,
        ledger,
        escape_ledger,
        committed_time,
        source,
    )
    binding = {
        "checkpoint_nominal_slot_omega0_inverse": checkpoint_nominal_slot_omega0_inverse,
        "checkpoint_observed_committed_cycle": checkpoint_observed_committed_cycle,
        "checkpoint_observed_committed_time_omega0_inverse": (
            checkpoint_observed_committed_time_omega0_inverse
        ),
        "restart_schema": probe.restart_schema,
        "startup_shock_ledger": ledger,
        "particle_escape_ledger": escape_ledger,
        "raw_particle_model": raw_model,
        "restart_control_binding": restart_controls,
        "active_particle_cohort": active_cohort,
        "retained_output_nominal_slots_after_checkpoint_omega0_inverse": (
            retained_output_nominal_slots_after_checkpoint_omega0_inverse
        ),
        "comparison_tolerances_max_absolute_difference": (
            comparison_tolerances_max_absolute_difference
        ),
    }
    validate_checkpoint_binding(binding, preregistration=policy)
    return copy.deepcopy(binding)


def _validated_values(values: object, kind: object, label: str) -> list[int | float]:
    _require(type(values) is list and values, f"{label}: expected non-empty flat list")
    _require(kind in {"integer", "float"}, f"{label}: unsupported comparison kind")
    for index, value in enumerate(values):
        if kind == "integer":
            _require(type(value) is int, f"{label}[{index}]: expected integer")
        else:
            _require(type(value) is float, f"{label}[{index}]: expected float")
            _require(math.isfinite(value), f"{label}[{index}]: expected finite float")
    return values


def _validate_observation(
    observation: object,
    policy: Mapping[str, object],
    label: str,
) -> dict[str, object]:
    observation_mapping = _keys(observation, _OBSERVATION_KEYS, label)
    validate_checkpoint_binding(
        observation_mapping["binding"], preregistration=policy
    )
    contract = policy["continuation_contract"]
    nominal_slots = contract[
        "retained_output_nominal_slots_after_checkpoint_omega0_inverse"
    ]
    outputs = observation_mapping["outputs_after_checkpoint"]
    _require(type(outputs) is list, f"{label}/outputs_after_checkpoint: expected list")
    _require(
        len(outputs) == len(nominal_slots),
        f"{label}/outputs_after_checkpoint: schedule length drift",
    )
    field_kinds = contract["comparison_field_kinds"]
    expected_fields = set(field_kinds)
    previous_cycle = observation_mapping["binding"][
        "checkpoint_observed_committed_cycle"
    ]
    previous_time = observation_mapping["binding"][
        "checkpoint_observed_committed_time_omega0_inverse"
    ]
    for output_index, (output, nominal_slot) in enumerate(zip(outputs, nominal_slots)):
        output_label = f"{label}/outputs_after_checkpoint[{output_index}]"
        output_mapping = _keys(output, _OUTPUT_KEYS, output_label)
        _strict_equal(
            output_mapping["nominal_slot_omega0_inverse"],
            nominal_slot,
            f"{output_label}/nominal_slot_omega0_inverse",
        )
        observed_cycle = _canonical_observed_cycle(
            output_mapping["observed_committed_cycle"],
            f"{output_label}/observed_committed_cycle",
        )
        observed_time = _canonical_observed_time(
            output_mapping["observed_committed_time_omega0_inverse"],
            f"{output_label}/observed_committed_time_omega0_inverse",
        )
        _require(
            observed_cycle > previous_cycle,
            f"{output_label}/observed_committed_cycle: sequence is not strictly increasing",
        )
        _require(
            observed_time > previous_time,
            f"{output_label}/observed_committed_time_omega0_inverse: "
            "sequence is not strictly increasing",
        )
        previous_cycle = observed_cycle
        previous_time = observed_time
        fields = _keys(output_mapping["fields"], expected_fields, f"{output_label}/fields")
        for field, kind in field_kinds.items():
            _validated_values(fields[field], kind, f"{output_label}/fields/{field}")
    return observation_mapping


def compare_deterministic_continuation_parity(
    uninterrupted: object,
    continued: object,
    *,
    preregistration: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Compare retained post-checkpoint values under one identical frozen binding."""
    policy = _validated_preregistration(preregistration)
    reference = _validate_observation(uninterrupted, policy, "uninterrupted")
    candidate = _validate_observation(continued, policy, "continued")
    _strict_equal(candidate["binding"], reference["binding"], "comparison binding identity")

    contract = policy["continuation_contract"]
    tolerances = contract["comparison_tolerances_max_absolute_difference"]
    maxima = {field: 0.0 for field in tolerances}
    for output_index, (expected_output, actual_output) in enumerate(
        zip(reference["outputs_after_checkpoint"], candidate["outputs_after_checkpoint"])
    ):
        _strict_equal(
            actual_output["nominal_slot_omega0_inverse"],
            expected_output["nominal_slot_omega0_inverse"],
            f"paired outputs[{output_index}]/nominal slot parity",
        )
        _strict_equal(
            actual_output["observed_committed_cycle"],
            expected_output["observed_committed_cycle"],
            f"paired outputs[{output_index}]/observed committed cycle parity",
        )
        _strict_equal(
            actual_output["observed_committed_time_omega0_inverse"],
            expected_output["observed_committed_time_omega0_inverse"],
            f"paired outputs[{output_index}]/observed committed time parity",
        )
        for field, tolerance in tolerances.items():
            expected_values = expected_output["fields"][field]
            actual_values = actual_output["fields"][field]
            _require(
                len(actual_values) == len(expected_values),
                f"continued/outputs_after_checkpoint[{output_index}]/fields/{field}: "
                "value-count drift",
            )
            maximum = max(
                (
                    abs(actual - expected)
                    for actual, expected in zip(actual_values, expected_values)
                ),
                default=0.0,
            )
            maxima[field] = max(maxima[field], float(maximum))
            _require(
                maximum <= tolerance,
                f"continued/outputs_after_checkpoint[{output_index}]/fields/{field}: "
                f"tolerance exceeded ({maximum!r} > {tolerance!r})",
            )
    return {
        "result": "pass_deterministic_continuation_parity",
        "checkpoint_nominal_slot_omega0_inverse": reference["binding"][
            "checkpoint_nominal_slot_omega0_inverse"
        ],
        "checkpoint_observed_committed_cycle": reference["binding"][
            "checkpoint_observed_committed_cycle"
        ],
        "checkpoint_observed_committed_time_omega0_inverse": reference["binding"][
            "checkpoint_observed_committed_time_omega0_inverse"
        ],
        "retained_output_nominal_slots_after_checkpoint_omega0_inverse": copy.deepcopy(
            reference["binding"][
                "retained_output_nominal_slots_after_checkpoint_omega0_inverse"
            ]
        ),
        "paired_output_observed_commits": [
            {
                "nominal_slot_omega0_inverse": output["nominal_slot_omega0_inverse"],
                "observed_committed_cycle": output["observed_committed_cycle"],
                "observed_committed_time_omega0_inverse": output[
                    "observed_committed_time_omega0_inverse"
                ],
            }
            for output in reference["outputs_after_checkpoint"]
        ],
        "maximum_absolute_difference_by_field": maxima,
    }
