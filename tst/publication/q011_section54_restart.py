#!/usr/bin/env python3
"""Pure fail-closed helpers for Q-011 Section 5.4 restart continuation.

This module intentionally performs no execution, scheduling, artifact
discovery, or output mutation. Callers supply retained restart bytes and
deterministically ordered extracted output values.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
import struct
from typing import Mapping


PREREGISTRATION = (
    Path(__file__).resolve().parent
    / "readiness"
    / "q011_section54_restart_continuation_preregistration_2026-06-01.json"
)

PIC_RESTART_MAGIC = 0x5049435253543031
EXPECTED_RESTART_SCHEMA = 7
EXPECTED_SHOCK_LEDGER_SCHEMA = 3
_PIC_RESTART_MARKER = struct.pack("<Q", PIC_RESTART_MAGIC)
_PIC_METADATA_FORMAT = "<15i"
_MODEL_INTEGER_COUNT = 31
_MODEL_REAL_COUNT = 37
_REAL_BYTES = 8
_MAX_LAYOUT_COUNT = 1 << 31
_MAX_SIGNED_INT = (1 << 31) - 1

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
_INTEGER_PATTERN = re.compile(r"-?[0-9]+")
_REAL_PATTERN = re.compile(r"-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")

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
_RETAINED_OUTPUT_SCHEDULE = [
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
    "schema_version": 1,
    "date": "2026-06-01",
    "gate": "Q-011",
    "claim_id": "CLAIM-PAPER-SHOCK-001",
    "qualification_effect": (
        "policy_freeze_only_no_execution_authorization_no_claim_closure"
    ),
    "predecessor_record": (
        "tst/publication/readiness/"
        "q011_section54_qualifying_campaign_preregistration_successor_v2_2026-06-01.json"
    ),
    "scope": (
        "Bounded additive preregistration for one future Q-011 Section 5.4 "
        "restart-continuation carrier. This record freezes the restart payload "
        "probe, startup shock-ledger extraction, checkpoint, retained "
        "post-checkpoint schedule and deterministic parity tolerances before "
        "any execution. It contains no result and authorizes no scheduler call."
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
    "continuation_contract": {
        "checkpoint_time_omega0_inverse": 500.0,
        "retained_output_schedule_after_checkpoint_omega0_inverse": list(
            _RETAINED_OUTPUT_SCHEDULE
        ),
        "comparison_tolerances_max_absolute_difference": dict(
            _COMPARISON_TOLERANCES
        ),
        "comparison_field_kinds": dict(_COMPARISON_FIELD_KINDS),
        "comparison_value_order": (
            "caller_supplied_deterministically_ordered_flat_values_per_field"
        ),
        "comparison_policy": (
            "require_identical_restart_schema_startup_cohort_ledger_retained_"
            "output_schedule_and_tolerances_before_field_comparison"
        ),
        "tolerance_boundary": (
            "These are AthenaK deterministic restart-continuation release "
            "screens selected before execution, not manuscript tolerances."
        ),
    },
    "execution_policy": {
        "status": "blocked_until_immutable_execution_record_binds_all_fields",
        "required_before_execution": [
            "checkpoint_time_omega0_inverse",
            "restart_schema",
            "startup_shock_ledger",
            "retained_output_schedule_after_checkpoint_omega0_inverse",
            "comparison_tolerances_max_absolute_difference",
            "clean_candidate_git_commit",
            "clean_frontier_executable_sha256",
            "qualifying_input_deck_sha256",
            "authorized_orion_campaign_root",
            "registered_frontier_submission_policy",
        ],
        "campaign_results_inspected": False,
        "scheduler_calls_authorized_by_this_record": False,
        "frontier_execution_authorized_by_this_record": False,
    },
    "limitations": [
        "This tranche is a non-executing policy and helper freeze only.",
        "It does not create a scheduler wrapper or authorize a Frontier submission.",
        "It does not close Q-011 qualification, independent recompute or external review.",
    ],
}

_BINDING_KEYS = {
    "checkpoint_time_omega0_inverse",
    "restart_schema",
    "startup_shock_ledger",
    "retained_output_schedule_after_checkpoint_omega0_inverse",
    "comparison_tolerances_max_absolute_difference",
}
_OBSERVATION_KEYS = {"binding", "outputs_after_checkpoint"}
_OUTPUT_KEYS = {"time_omega0_inverse", "fields"}


class RestartPolicyError(ValueError):
    """Raised when Q-011 restart policy or retained payloads fail closed."""


@dataclass(frozen=True)
class RestartPayloadProbe:
    """Validated schema-7 particle restart layout."""

    restart_schema: int
    meshblock_count: int
    real_fields_per_particle: int
    integer_fields_per_particle: int
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
    """Require exact identity with the frozen non-executing preregistration."""
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
    validate_preregistration(preregistration)
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
    restart_schema, meshblock_count, real_fields, integer_fields, *_ = metadata
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

    offset += _REAL_BYTES
    offset += _MODEL_INTEGER_COUNT * struct.calcsize("<i")
    offset += _MODEL_REAL_COUNT * _REAL_BYTES
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
        particle_count=particle_count,
        particle_real_offset=particle_real_offset,
        particle_integer_offset=particle_integer_offset,
        payload_end_offset=payload_end_offset,
    )


def _problem_parameters(payload: bytes, source: str) -> dict[str, str]:
    _require(type(payload) is bytes, f"{source}: restart payload must be bytes")
    parameter_end = payload.find(b"<par_end>")
    _require(parameter_end >= 0, f"{source}: restart header is missing <par_end>")
    try:
        text = payload[:parameter_end].decode("ascii")
    except UnicodeDecodeError as exc:
        raise RestartPolicyError(f"{source}: restart parameter header is not ASCII") from exc
    active_block: str | None = None
    parameters: dict[str, str] = {}
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
            continue
        if active_block != "problem":
            continue
        _require(
            line.count("=") == 1,
            f"{source}: malformed problem parameter at line {lineno}",
        )
        key, value = (item.strip() for item in line.split("=", 1))
        _require(key and value, f"{source}: empty problem parameter at line {lineno}")
        _require(
            key not in parameters,
            f"{source}: duplicate problem parameter: {key}",
        )
        parameters[key] = value
    return parameters


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


def _validate_startup_shock_ledger(ledger: object, label: str) -> dict[str, object]:
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
    return dict(_validate_startup_shock_ledger(ledger, f"{source}/startup_shock_ledger"))


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
        binding_mapping["checkpoint_time_omega0_inverse"],
        contract["checkpoint_time_omega0_inverse"],
        "checkpoint binding/checkpoint_time_omega0_inverse",
    )
    _strict_equal(
        binding_mapping["restart_schema"],
        policy["restart_payload_probe"]["restart_schema"],
        "checkpoint binding/restart_schema",
    )
    _validate_startup_shock_ledger(
        binding_mapping["startup_shock_ledger"],
        "checkpoint binding/startup_shock_ledger",
    )
    _strict_equal(
        binding_mapping["retained_output_schedule_after_checkpoint_omega0_inverse"],
        contract["retained_output_schedule_after_checkpoint_omega0_inverse"],
        "checkpoint binding/retained_output_schedule_after_checkpoint_omega0_inverse",
    )
    _strict_equal(
        binding_mapping["comparison_tolerances_max_absolute_difference"],
        contract["comparison_tolerances_max_absolute_difference"],
        "checkpoint binding/comparison_tolerances_max_absolute_difference",
    )


def bind_checkpoint_for_continuation(
    restart_payload: bytes,
    *,
    checkpoint_time_omega0_inverse: object,
    retained_output_schedule_after_checkpoint_omega0_inverse: object,
    comparison_tolerances_max_absolute_difference: object,
    preregistration: Mapping[str, object] | None = None,
    source: str = "<restart-bytes>",
) -> dict[str, object]:
    """Probe retained bytes and build a validated pre-execution checkpoint binding."""
    policy = _validated_preregistration(preregistration)
    probe = probe_schema7_restart_payload(restart_payload, source=source)
    ledger = extract_startup_shock_ledger(restart_payload, source=source)
    binding = {
        "checkpoint_time_omega0_inverse": checkpoint_time_omega0_inverse,
        "restart_schema": probe.restart_schema,
        "startup_shock_ledger": ledger,
        "retained_output_schedule_after_checkpoint_omega0_inverse": (
            retained_output_schedule_after_checkpoint_omega0_inverse
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
    schedule = contract["retained_output_schedule_after_checkpoint_omega0_inverse"]
    outputs = observation_mapping["outputs_after_checkpoint"]
    _require(type(outputs) is list, f"{label}/outputs_after_checkpoint: expected list")
    _require(
        len(outputs) == len(schedule),
        f"{label}/outputs_after_checkpoint: schedule length drift",
    )
    field_kinds = contract["comparison_field_kinds"]
    expected_fields = set(field_kinds)
    for output_index, (output, expected_time) in enumerate(zip(outputs, schedule)):
        output_label = f"{label}/outputs_after_checkpoint[{output_index}]"
        output_mapping = _keys(output, _OUTPUT_KEYS, output_label)
        _strict_equal(
            output_mapping["time_omega0_inverse"],
            expected_time,
            f"{output_label}/time_omega0_inverse",
        )
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
        "checkpoint_time_omega0_inverse": reference["binding"][
            "checkpoint_time_omega0_inverse"
        ],
        "retained_output_schedule_after_checkpoint_omega0_inverse": copy.deepcopy(
            reference["binding"][
                "retained_output_schedule_after_checkpoint_omega0_inverse"
            ]
        ),
        "maximum_absolute_difference_by_field": maxima,
    }
