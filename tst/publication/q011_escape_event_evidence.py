"""Strict reducer for Q011 raw physical-boundary escape evidence."""

from __future__ import annotations

import json
import math
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


class EscapeEvidenceError(ValueError):
    """Raised when raw escape evidence is incomplete or inconsistent."""


_BINARY_HEADER_MAGIC = b"Q011ESC1"
_BINARY_TRAILER_MAGIC = b"Q011END1"
_BINARY_HEADER_FIXED = struct.Struct("<8sIIIIIIIIQII")
_BINARY_HEADER_DOUBLES = struct.Struct("<26d")
_BINARY_EVENT = struct.Struct("<Q8I17d")
_BINARY_TRAILER = struct.Struct("<8sQQ")


_HEADER_FIELDS = {
    "record_type",
    "schema_version",
    "control_fingerprint",
    "basename",
    "rank",
    "nranks",
    "segment_start_cycle",
    "segment_start_time",
    "state_kind",
    "light_speed",
    "injected_species",
    "species_mass",
    "species_charge",
    "q_over_m",
    "macro_mass",
    "field_interpolation",
    "mesh_x1min",
    "mesh_x1max",
    "mesh_x2min",
    "mesh_x2max",
    "mesh_x3min",
    "mesh_x3max",
    "escape_audit_calls_at_start",
    "escape_last_audit_time_at_start",
    "global_escape_count_at_start",
    "global_escape_mass_at_start",
    "global_escape_momentum_x1_at_start",
    "global_escape_momentum_x2_at_start",
    "global_escape_momentum_x3_at_start",
    "global_escape_energy_at_start",
    "global_initial_escape_count_at_start",
    "global_escape_term_count_at_start",
    "global_escape_abs_mass_at_start",
    "global_escape_abs_momentum_x1_at_start",
    "global_escape_abs_momentum_x2_at_start",
    "global_escape_abs_momentum_x3_at_start",
    "global_escape_abs_energy_at_start",
}

_EVENT_FIELDS = {
    "record_type",
    "schema_version",
    "rank",
    "rank_event_index",
    "cycle",
    "stage",
    "audit_time",
    "tag",
    "source",
    "species",
    "destruction_reason",
    "physical_boundary_mask",
    "parent_gid",
    "x1",
    "x2",
    "x3",
    "state_x1",
    "state_x2",
    "state_x3",
    "state_kind",
    "light_speed",
    "species_mass",
    "species_charge",
    "q_over_m",
    "macro_weight",
    "macro_mass",
    "kinetic_energy_per_mass",
    "b1",
    "b2",
    "b3",
    "fluid_v1",
    "fluid_v2",
    "fluid_v3",
    "field_interpolation",
    "frame_velocity_x1",
    "shock_speed",
}

_TRAILER_FIELDS = {
    "record_type",
    "schema_version",
    "rank",
    "rank_event_count",
    "prefix_fnv1a64",
}

_HEADER_FLOATS = {
    field
    for field in _HEADER_FIELDS
    if field
    not in {
        "record_type",
        "schema_version",
        "control_fingerprint",
        "basename",
        "rank",
        "nranks",
        "segment_start_cycle",
        "state_kind",
        "injected_species",
        "field_interpolation",
        "escape_audit_calls_at_start",
    }
}

_EVENT_FLOATS = {
    field
    for field in _EVENT_FIELDS
    if field
    not in {
        "record_type",
        "schema_version",
        "rank",
        "rank_event_index",
        "cycle",
        "stage",
        "tag",
        "source",
        "species",
        "destruction_reason",
        "physical_boundary_mask",
        "parent_gid",
        "state_kind",
        "field_interpolation",
    }
}

_START_TO_LEDGER = {
    "global_escape_count_at_start": "count",
    "global_escape_mass_at_start": "mass",
    "global_escape_momentum_x1_at_start": "momentum_x1",
    "global_escape_momentum_x2_at_start": "momentum_x2",
    "global_escape_momentum_x3_at_start": "momentum_x3",
    "global_escape_energy_at_start": "energy",
    "global_initial_escape_count_at_start": "initial_count",
    "global_escape_term_count_at_start": "term_count",
    "global_escape_abs_mass_at_start": "abs_mass",
    "global_escape_abs_momentum_x1_at_start": "abs_momentum_x1",
    "global_escape_abs_momentum_x2_at_start": "abs_momentum_x2",
    "global_escape_abs_momentum_x3_at_start": "abs_momentum_x3",
    "global_escape_abs_energy_at_start": "abs_energy",
}


@dataclass(frozen=True)
class EscapeStream:
    path: str
    header: Mapping[str, Any]
    events: tuple[Mapping[str, Any], ...]
    trailer: Mapping[str, Any]


def _reject_duplicate_keys(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise EscapeEvidenceError(f"duplicate JSON field {key!r}")
        result[key] = value
    return result


def _parse_json_object(raw: bytes, *, path: str, line_number: int) -> dict[str, Any]:
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError as exc:
        raise EscapeEvidenceError(f"{path}:{line_number}: non-ASCII JSON") from exc
    try:
        value = json.loads(text, object_pairs_hook=_reject_duplicate_keys)
    except (json.JSONDecodeError, EscapeEvidenceError) as exc:
        raise EscapeEvidenceError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
    if not isinstance(value, dict):
        raise EscapeEvidenceError(f"{path}:{line_number}: record is not an object")
    return value


def _require_exact_fields(
    record: Mapping[str, Any], expected: set[str], *, context: str
) -> None:
    actual = set(record)
    missing = sorted(expected - actual)
    extra = sorted(actual - expected)
    if missing or extra:
        raise EscapeEvidenceError(
            f"{context}: schema mismatch; missing={missing!r}, extra={extra!r}"
        )


def _require_int(value: Any, *, context: str, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise EscapeEvidenceError(f"{context}: expected integer")
    if minimum is not None and value < minimum:
        raise EscapeEvidenceError(f"{context}: value is below {minimum}")
    return value


def _require_finite(value: Any, *, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise EscapeEvidenceError(f"{context}: expected number")
    converted = float(value)
    if not math.isfinite(converted):
        raise EscapeEvidenceError(f"{context}: non-finite number")
    return converted


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=2.0e-14, abs_tol=2.0e-14)


def _fnv1a64(payload: bytes) -> int:
    value = 14695981039346656037
    for byte in payload:
        value ^= byte
        value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return value


def _specific_energy(event: Mapping[str, Any]) -> float:
    state_squared = math.fsum(
        float(event[field]) ** 2 for field in ("state_x1", "state_x2", "state_x3")
    )
    if event["state_kind"] == "velocity":
        return 0.5 * state_squared
    light_speed = float(event["light_speed"])
    gamma = math.sqrt(1.0 + state_squared / (light_speed * light_speed))
    return state_squared / (gamma + 1.0)


def _validate_header(header: Mapping[str, Any], *, context: str) -> None:
    _require_exact_fields(header, _HEADER_FIELDS, context=context)
    if header["record_type"] != "q011_escape_event_stream_header":
        raise EscapeEvidenceError(f"{context}: wrong record_type")
    if header["schema_version"] != 1:
        raise EscapeEvidenceError(f"{context}: unsupported schema_version")
    rank = _require_int(header["rank"], context=f"{context}.rank", minimum=0)
    nranks = _require_int(header["nranks"], context=f"{context}.nranks", minimum=1)
    if rank >= nranks:
        raise EscapeEvidenceError(f"{context}: rank is outside nranks")
    _require_int(
        header["segment_start_cycle"],
        context=f"{context}.segment_start_cycle",
        minimum=0,
    )
    _require_int(
        header["injected_species"],
        context=f"{context}.injected_species",
        minimum=0,
    )
    _require_int(
        header["escape_audit_calls_at_start"],
        context=f"{context}.escape_audit_calls_at_start",
        minimum=0,
    )
    if not isinstance(header["control_fingerprint"], str) or not header[
        "control_fingerprint"
    ].startswith("v1:"):
        raise EscapeEvidenceError(f"{context}: invalid control_fingerprint")
    if not isinstance(header["basename"], str) or not header["basename"]:
        raise EscapeEvidenceError(f"{context}: invalid basename")
    if header["state_kind"] not in {"velocity", "momentum_per_mass"}:
        raise EscapeEvidenceError(f"{context}: invalid state_kind")
    if header["field_interpolation"] != "tsc_bcc0_w0":
        raise EscapeEvidenceError(f"{context}: invalid field_interpolation")
    for field in _HEADER_FLOATS:
        _require_finite(header[field], context=f"{context}.{field}")
    for field in ("light_speed", "species_mass", "macro_mass"):
        if float(header[field]) <= 0.0:
            raise EscapeEvidenceError(f"{context}.{field}: must be positive")
    for lower, upper in (
        ("mesh_x1min", "mesh_x1max"),
        ("mesh_x2min", "mesh_x2max"),
        ("mesh_x3min", "mesh_x3max"),
    ):
        if float(header[upper]) < float(header[lower]):
            raise EscapeEvidenceError(f"{context}: invalid mesh bounds")
    if not _close(
        float(header["q_over_m"]),
        float(header["species_charge"]) / float(header["species_mass"]),
    ):
        raise EscapeEvidenceError(f"{context}: q_over_m is inconsistent")
    count = float(header["global_escape_count_at_start"])
    mass = float(header["global_escape_mass_at_start"])
    term_count = float(header["global_escape_term_count_at_start"])
    if count < 0.0 or count != math.floor(count):
        raise EscapeEvidenceError(f"{context}: invalid starting escape count")
    if term_count != count:
        raise EscapeEvidenceError(f"{context}: starting term count mismatch")
    if not _close(mass, count * float(header["macro_mass"])):
        raise EscapeEvidenceError(f"{context}: starting mass mismatch")
    if float(header["global_initial_escape_count_at_start"]) != 0.0:
        raise EscapeEvidenceError(f"{context}: initial-particle escape is nonzero")
    for field in (
        "global_escape_energy_at_start",
        "global_escape_abs_mass_at_start",
        "global_escape_abs_momentum_x1_at_start",
        "global_escape_abs_momentum_x2_at_start",
        "global_escape_abs_momentum_x3_at_start",
        "global_escape_abs_energy_at_start",
    ):
        if float(header[field]) < 0.0:
            raise EscapeEvidenceError(f"{context}.{field}: must be nonnegative")


def _validate_event(
    event: Mapping[str, Any],
    header: Mapping[str, Any],
    *,
    context: str,
    expected_index: int,
) -> None:
    _require_exact_fields(event, _EVENT_FIELDS, context=context)
    if event["record_type"] != "q011_escape_event" or event["schema_version"] != 1:
        raise EscapeEvidenceError(f"{context}: wrong event schema")
    integer_minima = {
        "rank": 0,
        "rank_event_index": 0,
        "cycle": 0,
        "stage": 1,
        "tag": 0,
        "source": 0,
        "species": 0,
        "destruction_reason": 0,
        "physical_boundary_mask": 0,
        "parent_gid": 0,
    }
    for field, minimum in integer_minima.items():
        _require_int(event[field], context=f"{context}.{field}", minimum=minimum)
    for field in _EVENT_FLOATS:
        _require_finite(event[field], context=f"{context}.{field}")
    if event["rank"] != header["rank"] or event["rank_event_index"] != expected_index:
        raise EscapeEvidenceError(f"{context}: rank or rank_event_index mismatch")
    if event["cycle"] < header["segment_start_cycle"] or event["stage"] not in {1, 2}:
        raise EscapeEvidenceError(f"{context}: invalid cycle/stage chronology")
    if float(event["audit_time"]) < float(header["segment_start_time"]):
        raise EscapeEvidenceError(f"{context}: audit_time predates segment")
    if event["source"] != 1 or event["species"] != header["injected_species"]:
        raise EscapeEvidenceError(f"{context}: source/species mismatch")
    if event["destruction_reason"] != 1 or event["physical_boundary_mask"] != 2:
        raise EscapeEvidenceError(f"{context}: event is not outer-x1 physical escape")
    if float(event["x1"]) < float(header["mesh_x1max"]):
        raise EscapeEvidenceError(f"{context}: x1 does not cross the outer boundary")
    for field in (
        "state_kind",
        "light_speed",
        "species_mass",
        "species_charge",
        "q_over_m",
        "macro_mass",
        "field_interpolation",
    ):
        if event[field] != header[field]:
            raise EscapeEvidenceError(f"{context}: {field} differs from header")
    if float(event["macro_weight"]) != 1.0:
        raise EscapeEvidenceError(f"{context}: macro_weight is not unity")
    if float(event["shock_speed"]) <= 0.0:
        raise EscapeEvidenceError(f"{context}: shock_speed is not positive")
    expected_energy = _specific_energy(event)
    if not _close(float(event["kinetic_energy_per_mass"]), expected_energy):
        raise EscapeEvidenceError(f"{context}: kinetic energy is not recomputable")


def _parse_binary_stream_bytes(payload: bytes, *, path: str) -> EscapeStream:
    minimum = (
        _BINARY_HEADER_FIXED.size
        + _BINARY_HEADER_DOUBLES.size
        + _BINARY_TRAILER.size
    )
    if len(payload) < minimum:
        raise EscapeEvidenceError(f"{path}: binary stream is truncated")
    fixed = _BINARY_HEADER_FIXED.unpack_from(payload, 0)
    (
        magic,
        schema,
        header_bytes,
        event_bytes,
        rank,
        nranks,
        start_cycle,
        injected_species,
        state_kind_code,
        audit_calls,
        basename_bytes,
        fingerprint_bytes,
    ) = fixed
    if magic != _BINARY_HEADER_MAGIC or schema != 1:
        raise EscapeEvidenceError(f"{path}: unsupported binary header")
    if event_bytes != _BINARY_EVENT.size:
        raise EscapeEvidenceError(f"{path}: binary event size drifted")
    expected_header = (
        _BINARY_HEADER_FIXED.size
        + _BINARY_HEADER_DOUBLES.size
        + basename_bytes
        + fingerprint_bytes
    )
    if header_bytes != expected_header or header_bytes > len(payload) - _BINARY_TRAILER.size:
        raise EscapeEvidenceError(f"{path}: binary header size drifted")
    doubles = _BINARY_HEADER_DOUBLES.unpack_from(
        payload, _BINARY_HEADER_FIXED.size
    )
    strings_offset = _BINARY_HEADER_FIXED.size + _BINARY_HEADER_DOUBLES.size
    try:
        basename = payload[
            strings_offset : strings_offset + basename_bytes
        ].decode("ascii")
        fingerprint = payload[
            strings_offset + basename_bytes : header_bytes
        ].decode("ascii")
    except UnicodeDecodeError as exc:
        raise EscapeEvidenceError(f"{path}: binary metadata is not ASCII") from exc
    (
        start_time,
        light_speed,
        species_mass,
        species_charge,
        q_over_m,
        macro_mass,
        mesh_x1min,
        mesh_x1max,
        mesh_x2min,
        mesh_x2max,
        mesh_x3min,
        mesh_x3max,
        last_audit_time,
        start_count,
        start_mass,
        start_momentum_x1,
        start_momentum_x2,
        start_momentum_x3,
        start_energy,
        start_initial_count,
        start_term_count,
        start_abs_mass,
        start_abs_momentum_x1,
        start_abs_momentum_x2,
        start_abs_momentum_x3,
        start_abs_energy,
    ) = doubles
    if state_kind_code not in {0, 1}:
        raise EscapeEvidenceError(f"{path}: invalid binary state kind")
    state_kind = "momentum_per_mass" if state_kind_code == 1 else "velocity"
    header = {
        "record_type": "q011_escape_event_stream_header",
        "schema_version": schema,
        "control_fingerprint": fingerprint,
        "basename": basename,
        "rank": rank,
        "nranks": nranks,
        "segment_start_cycle": start_cycle,
        "segment_start_time": start_time,
        "state_kind": state_kind,
        "light_speed": light_speed,
        "injected_species": injected_species,
        "species_mass": species_mass,
        "species_charge": species_charge,
        "q_over_m": q_over_m,
        "macro_mass": macro_mass,
        "field_interpolation": "tsc_bcc0_w0",
        "mesh_x1min": mesh_x1min,
        "mesh_x1max": mesh_x1max,
        "mesh_x2min": mesh_x2min,
        "mesh_x2max": mesh_x2max,
        "mesh_x3min": mesh_x3min,
        "mesh_x3max": mesh_x3max,
        "escape_audit_calls_at_start": audit_calls,
        "escape_last_audit_time_at_start": last_audit_time,
        "global_escape_count_at_start": start_count,
        "global_escape_mass_at_start": start_mass,
        "global_escape_momentum_x1_at_start": start_momentum_x1,
        "global_escape_momentum_x2_at_start": start_momentum_x2,
        "global_escape_momentum_x3_at_start": start_momentum_x3,
        "global_escape_energy_at_start": start_energy,
        "global_initial_escape_count_at_start": start_initial_count,
        "global_escape_term_count_at_start": start_term_count,
        "global_escape_abs_mass_at_start": start_abs_mass,
        "global_escape_abs_momentum_x1_at_start": start_abs_momentum_x1,
        "global_escape_abs_momentum_x2_at_start": start_abs_momentum_x2,
        "global_escape_abs_momentum_x3_at_start": start_abs_momentum_x3,
        "global_escape_abs_energy_at_start": start_abs_energy,
    }
    trailer_magic, event_count, prefix_hash = _BINARY_TRAILER.unpack_from(
        payload, len(payload) - _BINARY_TRAILER.size
    )
    if trailer_magic != _BINARY_TRAILER_MAGIC:
        raise EscapeEvidenceError(f"{path}: binary trailer magic drifted")
    expected_size = (
        header_bytes + event_count * event_bytes + _BINARY_TRAILER.size
    )
    if len(payload) != expected_size:
        raise EscapeEvidenceError(f"{path}: binary event count/size mismatch")
    actual_hash = _fnv1a64(payload[: -_BINARY_TRAILER.size])
    if actual_hash != prefix_hash:
        raise EscapeEvidenceError(f"{path}: prefix FNV-1a mismatch")
    events = []
    offset = header_bytes
    for event_number in range(event_count):
        values = _BINARY_EVENT.unpack_from(payload, offset)
        offset += event_bytes
        (
            event_index,
            cycle,
            stage,
            tag,
            source,
            species,
            destruction_reason,
            boundary_mask,
            parent_gid,
            audit_time,
            x1,
            x2,
            x3,
            state_x1,
            state_x2,
            state_x3,
            event_q_over_m,
            macro_weight,
            b1,
            b2,
            b3,
            fluid_v1,
            fluid_v2,
            fluid_v3,
            frame_velocity_x1,
            shock_speed,
        ) = values
        event = {
            "record_type": "q011_escape_event",
            "schema_version": 1,
            "rank": rank,
            "rank_event_index": event_index,
            "cycle": cycle,
            "stage": stage,
            "audit_time": audit_time,
            "tag": tag,
            "source": source,
            "species": species,
            "destruction_reason": destruction_reason,
            "physical_boundary_mask": boundary_mask,
            "parent_gid": parent_gid,
            "x1": x1,
            "x2": x2,
            "x3": x3,
            "state_x1": state_x1,
            "state_x2": state_x2,
            "state_x3": state_x3,
            "state_kind": state_kind,
            "light_speed": light_speed,
            "species_mass": species_mass,
            "species_charge": species_charge,
            "q_over_m": event_q_over_m,
            "macro_weight": macro_weight,
            "macro_mass": macro_mass,
            "kinetic_energy_per_mass": 0.0,
            "b1": b1,
            "b2": b2,
            "b3": b3,
            "fluid_v1": fluid_v1,
            "fluid_v2": fluid_v2,
            "fluid_v3": fluid_v3,
            "field_interpolation": "tsc_bcc0_w0",
            "frame_velocity_x1": frame_velocity_x1,
            "shock_speed": shock_speed,
        }
        event["kinetic_energy_per_mass"] = _specific_energy(event)
        _validate_event(
            event,
            header,
            context=f"{path}:event[{event_number}]",
            expected_index=event_number,
        )
        events.append(event)
    _validate_header(header, context=f"{path}:header")
    trailer = {
        "record_type": "q011_escape_event_stream_trailer",
        "schema_version": 1,
        "rank": rank,
        "rank_event_count": event_count,
        "prefix_fnv1a64": f"{prefix_hash:016x}",
    }
    return EscapeStream(
        path=path, header=header, events=tuple(events), trailer=trailer
    )


def parse_stream_bytes(payload: bytes, *, path: str = "<bytes>") -> EscapeStream:
    if payload.startswith(_BINARY_HEADER_MAGIC):
        return _parse_binary_stream_bytes(payload, path=path)
    if not payload or not payload.endswith(b"\n"):
        raise EscapeEvidenceError(f"{path}: stream is empty or lacks final newline")
    lines = payload.splitlines(keepends=True)
    if len(lines) < 2:
        raise EscapeEvidenceError(f"{path}: header/trailer are incomplete")
    if any(not line.endswith(b"\n") or line == b"\n" for line in lines):
        raise EscapeEvidenceError(f"{path}: blank or unterminated record")
    records = [
        _parse_json_object(line[:-1], path=path, line_number=index + 1)
        for index, line in enumerate(lines)
    ]
    header = records[0]
    trailer = records[-1]
    events = records[1:-1]
    _validate_header(header, context=f"{path}:header")
    _require_exact_fields(trailer, _TRAILER_FIELDS, context=f"{path}:trailer")
    if (
        trailer["record_type"] != "q011_escape_event_stream_trailer"
        or trailer["schema_version"] != 1
    ):
        raise EscapeEvidenceError(f"{path}: invalid trailer schema")
    _require_int(trailer["rank"], context=f"{path}:trailer.rank", minimum=0)
    _require_int(
        trailer["rank_event_count"],
        context=f"{path}:trailer.rank_event_count",
        minimum=0,
    )
    if trailer["rank"] != header["rank"] or trailer["rank_event_count"] != len(events):
        raise EscapeEvidenceError(f"{path}: trailer count/rank mismatch")
    expected_hash = f"{_fnv1a64(b''.join(lines[:-1])):016x}"
    if trailer["prefix_fnv1a64"] != expected_hash:
        raise EscapeEvidenceError(f"{path}: prefix FNV-1a mismatch")
    previous: tuple[int, int, float] | None = None
    for index, event in enumerate(events):
        _validate_event(
            event,
            header,
            context=f"{path}:event[{index}]",
            expected_index=index,
        )
        chronology = (event["cycle"], event["stage"], float(event["audit_time"]))
        if previous is not None and chronology < previous:
            raise EscapeEvidenceError(f"{path}: event chronology regressed")
        previous = chronology
    return EscapeStream(path=path, header=header, events=tuple(events), trailer=trailer)


def parse_stream(path: str | Path) -> EscapeStream:
    resolved = Path(path)
    return parse_stream_bytes(resolved.read_bytes(), path=str(resolved))


def _starting_ledger(header: Mapping[str, Any]) -> dict[str, float]:
    return {
        ledger_name: float(header[header_name])
        for header_name, ledger_name in _START_TO_LEDGER.items()
    }


def _event_increment(events: Iterable[Mapping[str, Any]]) -> dict[str, float]:
    event_list = list(events)
    masses = [float(event["macro_mass"]) for event in event_list]
    momenta = {
        axis: [
            float(event["macro_mass"]) * float(event[f"state_x{axis}"])
            for event in event_list
        ]
        for axis in (1, 2, 3)
    }
    energies = [
        float(event["macro_mass"]) * _specific_energy(event) for event in event_list
    ]
    return {
        "count": float(len(event_list)),
        "mass": math.fsum(masses),
        "momentum_x1": math.fsum(momenta[1]),
        "momentum_x2": math.fsum(momenta[2]),
        "momentum_x3": math.fsum(momenta[3]),
        "energy": math.fsum(energies),
        "initial_count": 0.0,
        "term_count": float(len(event_list)),
        "abs_mass": math.fsum(abs(value) for value in masses),
        "abs_momentum_x1": math.fsum(abs(value) for value in momenta[1]),
        "abs_momentum_x2": math.fsum(abs(value) for value in momenta[2]),
        "abs_momentum_x3": math.fsum(abs(value) for value in momenta[3]),
        "abs_energy": math.fsum(abs(value) for value in energies),
    }


def _add_ledgers(left: Mapping[str, float], right: Mapping[str, float]) -> dict[str, float]:
    return {key: float(left[key]) + float(right[key]) for key in left}


def _ledgers_agree(left: Mapping[str, float], right: Mapping[str, float]) -> bool:
    return set(left) == set(right) and all(_close(left[key], right[key]) for key in left)


def reduce_streams(streams: Iterable[EscapeStream]) -> dict[str, Any]:
    materialized = list(streams)
    if not materialized:
        raise EscapeEvidenceError("no escape streams supplied")
    groups: dict[tuple[Any, ...], list[EscapeStream]] = {}
    for stream in materialized:
        header = stream.header
        key = (
            header["control_fingerprint"],
            header["basename"],
            header["segment_start_cycle"],
            header["segment_start_time"],
        )
        groups.setdefault(key, []).append(stream)
    ordered_groups = sorted(
        groups.values(),
        key=lambda group: (
            group[0].header["segment_start_cycle"],
            group[0].header["segment_start_time"],
            group[0].header["basename"],
        ),
    )
    previous_final: dict[str, float] | None = None
    all_tags: set[int] = set()
    segment_results: list[dict[str, Any]] = []
    seen_start_slots: set[tuple[int, float]] = set()
    for group_index, group in enumerate(ordered_groups):
        representative = group[0].header
        slot = (
            representative["segment_start_cycle"],
            float(representative["segment_start_time"]),
        )
        if slot in seen_start_slots:
            raise EscapeEvidenceError("ambiguous duplicate segment start slot")
        seen_start_slots.add(slot)
        nranks = representative["nranks"]
        ranks = {stream.header["rank"] for stream in group}
        if len(group) != nranks or ranks != set(range(nranks)):
            raise EscapeEvidenceError(
                f"segment {representative['basename']!r}: incomplete rank set"
            )
        reference = dict(representative)
        reference.pop("rank")
        for stream in group[1:]:
            candidate = dict(stream.header)
            candidate.pop("rank")
            if candidate != reference:
                raise EscapeEvidenceError(
                    f"segment {representative['basename']!r}: rank headers differ"
                )
        expected_audits = 2 * int(representative["segment_start_cycle"])
        if representative["escape_audit_calls_at_start"] != expected_audits:
            raise EscapeEvidenceError("segment start VL2 audit count is inconsistent")
        expected_last_time = 0.0 if expected_audits == 0 else float(
            representative["segment_start_time"]
        )
        if not _close(
            float(representative["escape_last_audit_time_at_start"]),
            expected_last_time,
        ):
            raise EscapeEvidenceError("segment start audit time is inconsistent")
        start = _starting_ledger(representative)
        if group_index == 0:
            if slot != (0, 0.0):
                raise EscapeEvidenceError(
                    "escape evidence chain does not begin at cycle/time zero"
                )
            if not _ledgers_agree(start, {key: 0.0 for key in start}):
                raise EscapeEvidenceError(
                    "escape evidence chain begins with a nonzero ledger"
                )
        if previous_final is not None and not _ledgers_agree(start, previous_final):
            raise EscapeEvidenceError("continuation ledger does not join prior segment")
        events = [event for stream in group for event in stream.events]
        for event in events:
            tag = int(event["tag"])
            if tag in all_tags:
                raise EscapeEvidenceError(f"duplicate escaped particle tag {tag}")
            all_tags.add(tag)
        increment = _event_increment(events)
        final = _add_ledgers(start, increment)
        segment_results.append(
            {
                "basename": representative["basename"],
                "start_cycle": representative["segment_start_cycle"],
                "start_time": representative["segment_start_time"],
                "nranks": nranks,
                "event_count": len(events),
                "starting_ledger": start,
                "increment": increment,
                "derived_final_ledger": final,
            }
        )
        previous_final = final
    assert previous_final is not None
    return {
        "record_type": "q011_escape_event_independent_reduction",
        "schema_version": 1,
        "control_fingerprint": ordered_groups[0][0].header["control_fingerprint"],
        "stream_count": len(materialized),
        "segment_count": len(segment_results),
        "unique_particle_tags": len(all_tags),
        "segments": segment_results,
        "derived_final_ledger": previous_final,
    }


def reduce_paths(paths: Iterable[str | Path]) -> dict[str, Any]:
    return reduce_streams(parse_stream(path) for path in paths)


def assert_matches_restart_ledger(
    reduction: Mapping[str, Any], restart_parameters: Mapping[str, Any]
) -> None:
    expected = {
        "count": float(restart_parameters["ps_escaped_injected_cr_count_global"]),
        "mass": float(restart_parameters["ps_escaped_injected_cr_mass_global"]),
        "momentum_x1": float(
            restart_parameters["ps_escaped_injected_cr_momentum_x1_global"]
        ),
        "momentum_x2": float(
            restart_parameters["ps_escaped_injected_cr_momentum_x2_global"]
        ),
        "momentum_x3": float(
            restart_parameters["ps_escaped_injected_cr_momentum_x3_global"]
        ),
        "energy": float(restart_parameters["ps_escaped_injected_cr_energy_global"]),
        "initial_count": float(
            restart_parameters["ps_escaped_initial_cr_count_global"]
        ),
        "term_count": float(
            restart_parameters["ps_escaped_injected_cr_term_count_global"]
        ),
        "abs_mass": float(
            restart_parameters["ps_escaped_injected_cr_abs_mass_global"]
        ),
        "abs_momentum_x1": float(
            restart_parameters["ps_escaped_injected_cr_abs_momentum_x1_global"]
        ),
        "abs_momentum_x2": float(
            restart_parameters["ps_escaped_injected_cr_abs_momentum_x2_global"]
        ),
        "abs_momentum_x3": float(
            restart_parameters["ps_escaped_injected_cr_abs_momentum_x3_global"]
        ),
        "abs_energy": float(
            restart_parameters["ps_escaped_injected_cr_abs_energy_global"]
        ),
    }
    actual = reduction["derived_final_ledger"]
    if not _ledgers_agree(actual, expected):
        raise EscapeEvidenceError(
            f"raw-event reduction does not match restart ledger: "
            f"actual={actual!r}, expected={expected!r}"
        )
