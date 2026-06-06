#!/usr/bin/env python3
"""Pure supplied-artifact MHD+CR conservation closure for the Q011 successor."""

from __future__ import annotations

import hashlib
import json
import math
import re
import struct
from typing import Any, Mapping, Sequence

import numpy as np

if __package__:
    from . import q011_section54_restart as restart_layout
    from . import q019_nonlinear_bell_particle_state as particle_state
else:
    import q011_section54_restart as restart_layout
    import q019_nonlinear_bell_particle_state as particle_state


SCHEMA_VERSION = 1
RECORD_TYPE = "q011_section54_exact_conservation_closure_successor_v1"
_LINEAGE_KEYS = {
    "attempt_id",
    "attempt_sha256",
    "candidate_commit",
    "candidate_sha256",
    "deck_sha256",
    "executable_sha256",
}
_ARTIFACT_KEYS = {"name", "payload", "sha256", "byte_count", "lineage"}
_MHD_LABELS = (
    "time", "dt", "mass", "1-mom", "2-mom", "3-mom", "tot-E",
    "1-KE", "2-KE", "3-KE", "1-ME", "2-ME", "3-ME",
)
_USER_LABELS = (
    "time", "dt", "ext_mass", "ext_mom1", "ext_mom2", "ext_mom3",
    "ext_etot", "mhd_bmass", "mhd_bmom1", "mhd_bmom2", "mhd_bmom3",
    "mhd_betot", "cr_bmass", "cr_betot",
)
_COMPONENTS = ("mass", "momentum_x1", "momentum_x2", "momentum_x3", "energy")
_CONS_PREFIXES = (
    "ps_cons_mhd_boundary",
    "ps_cons_particle_reflect",
    "ps_cons_particle_escape",
    "ps_cons_gas_subtracted",
)
_ESCAPE_INTEGER_FIELDS = (
    "ps_escape_ledger_schema",
    "ps_escape_audit_calls",
)
_ESCAPE_BOOLEAN_FIELDS = ("ps_escape_ledger_complete",)
_ESCAPE_REAL_FIELDS = (
    "ps_escape_last_audit_time",
    "ps_escaped_injected_cr_count_global",
    "ps_escaped_injected_cr_mass_global",
    "ps_escaped_injected_cr_momentum_x1_global",
    "ps_escaped_injected_cr_momentum_x2_global",
    "ps_escaped_injected_cr_momentum_x3_global",
    "ps_escaped_injected_cr_energy_global",
    "ps_escaped_initial_cr_count_global",
)
_MESH_METADATA_MAGIC = 0x4154484B4D455348
_MESH_METADATA_VERSION = 2
_REGION_INDCS_COUNT = 19
_LOGICAL_LOCATION_COUNT = 4
_IDEAL_MHD_CONSERVED_COUNT = 5
_HEX40 = re.compile(r"[0-9a-f]{40}")
_HEX64 = re.compile(r"[0-9a-f]{64}")
_INT = re.compile(r"-?[0-9]+")
_REAL = re.compile(r"-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")


class ConservationClosureError(ValueError):
    """Raised when exact closure inputs are incomplete, inconsistent, or ambiguous."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ConservationClosureError(message)


def _validate_fixed_uniform_topology(
    *,
    logical_locations: Sequence[tuple[int, int, int, int]],
    deck_mesh_cells: Sequence[int],
    deck_block_cells: Sequence[int],
    root_level: int,
    nmb_total: int,
    label: str,
) -> None:
    """Require one complete, unique root-level tiling of the configured mesh."""
    _require(
        len(deck_mesh_cells) == 3
        and len(deck_block_cells) == 3
        and all(value > 0 for value in (*deck_mesh_cells, *deck_block_cells)),
        f"{label}: invalid fixed-uniform mesh geometry",
    )
    _require(
        all(
            mesh_cells % block_cells == 0
            for mesh_cells, block_cells in zip(deck_mesh_cells, deck_block_cells)
        ),
        f"{label}: MeshBlocks do not tile the fixed-uniform mesh",
    )
    root_counts = tuple(
        mesh_cells // block_cells
        for mesh_cells, block_cells in zip(deck_mesh_cells, deck_block_cells)
    )
    expected_root_level = 0
    while (1 << expected_root_level) < max(root_counts):
        expected_root_level += 1
    _require(
        root_level == expected_root_level,
        f"{label}: root level is inconsistent with the fixed-uniform root grid",
    )
    expected_count = math.prod(root_counts)
    _require(
        nmb_total == expected_count and len(logical_locations) == nmb_total,
        f"{label}: fixed-uniform MeshBlock count is incomplete",
    )
    occupied: set[tuple[int, int, int]] = set()
    for location in logical_locations:
        coordinate = location[:3]
        _require(
            location[3] == root_level
            and all(
                0 <= value < extent
                for value, extent in zip(coordinate, root_counts)
            ),
            f"{label}: fixed-uniform logical coordinate is out of range",
        )
        _require(
            coordinate not in occupied,
            f"{label}: fixed-uniform logical coordinate is duplicated",
        )
        occupied.add(coordinate)
    _require(
        len(occupied) == expected_count,
        f"{label}: fixed-uniform logical coordinates do not completely tile the mesh",
    )


def canonical_json_bytes(value: object) -> bytes:
    """Serialize a result deterministically without permitting non-finite aliases."""
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")


def _mapping(value: object, keys: set[str], label: str) -> Mapping[str, object]:
    _require(type(value) is dict, f"{label} must be an object")
    _require(set(value) == keys, f"{label} keys drifted")
    return value


def _lineage(value: object) -> dict[str, str]:
    item = _mapping(value, _LINEAGE_KEYS, "lineage")
    for key, member in item.items():
        _require(type(member) is str and member, f"lineage/{key} must be a string")
    _require(_HEX40.fullmatch(item["candidate_commit"]) is not None,
             "lineage/candidate_commit must be lowercase git hex")
    for key in ("attempt_sha256", "candidate_sha256", "deck_sha256", "executable_sha256"):
        _require(_HEX64.fullmatch(item[key]) is not None,
                 f"lineage/{key} must be lowercase sha256")
    _require(
        re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", item["attempt_id"]) is not None,
        "lineage/attempt_id is not canonical",
    )
    return dict(item)


def _artifact(value: object, lineage: Mapping[str, str], label: str) -> dict[str, object]:
    item = _mapping(value, _ARTIFACT_KEYS, label)
    _require(type(item["name"]) is str and item["name"], f"{label}/name must be nonempty")
    _require(type(item["payload"]) is bytes, f"{label}/payload must be bytes")
    _require(type(item["byte_count"]) is int and item["byte_count"] >= 0,
             f"{label}/byte_count must be a nonnegative integer")
    _require(type(item["sha256"]) is str and _HEX64.fullmatch(item["sha256"]) is not None,
             f"{label}/sha256 must be lowercase sha256")
    _require(item["lineage"] == lineage, f"{label}/lineage differs from bound lineage")
    payload = item["payload"]
    _require(len(payload) == item["byte_count"], f"{label}/byte_count mismatch")
    _require(hashlib.sha256(payload).hexdigest() == item["sha256"],
             f"{label}/sha256 mismatch")
    return {
        "name": item["name"],
        "payload": payload,
        "sha256": item["sha256"],
        "byte_count": item["byte_count"],
    }


def _artifact_record(item: Mapping[str, object]) -> dict[str, object]:
    return {key: item[key] for key in ("name", "sha256", "byte_count")}


def _parse_input(
    payload: bytes, label: str, *, allow_empty_values: bool = False
) -> dict[str, dict[str, str]]:
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ConservationClosureError(f"{label} is not ASCII") from exc
    blocks: dict[str, dict[str, str]] = {}
    active: str | None = None
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<"):
            _require(line.endswith(">") and line.count("<") == 1 and line.count(">") == 1,
                     f"{label}: malformed block at line {lineno}")
            active = line[1:-1]
            _require(active not in blocks, f"{label}: duplicate block {active}")
            blocks[active] = {}
            continue
        _require(active is not None and line.count("=") == 1,
                 f"{label}: malformed parameter at line {lineno}")
        key, value = (part.strip() for part in line.split("=", 1))
        _require(key and (value or allow_empty_values) and key not in blocks[active],
                 f"{label}: duplicate or empty parameter at line {lineno}")
        blocks[active][key] = value
    return blocks


def _deck_physics(payload: bytes) -> dict[str, object]:
    blocks = _parse_input(payload, "deck")
    for block in ("mesh", "time", "mhd", "coord", "particles", "problem"):
        _require(block in blocks, f"deck is missing <{block}>")
    particles = blocks["particles"]
    problem = blocks["problem"]
    required = {
        ("mesh", "ix1_bc", "reflect"),
        ("mesh", "ix2_bc", "periodic"),
        ("mesh", "ox2_bc", "periodic"),
        ("mesh", "nx3", "1"),
        ("mesh", "ix3_bc", "periodic"),
        ("mesh", "ox3_bc", "periodic"),
        ("time", "integrator", "rk2"),
        ("mhd", "eos", "ideal"),
        ("coord", "special_rel", "false"),
        ("coord", "general_rel", "false"),
        ("particles", "particle_type", "cosmic_ray"),
        ("particles", "pusher", "boris_tsc"),
        ("particles", "pic_enable_2d3v", "true"),
        ("particles", "deposit_moments", "true"),
        ("particles", "deposit_order", "2"),
        ("particles", "couple_moments_to_mhd", "true"),
        ("particles", "couple_fluid_feedback_order", "mhd_src_terms"),
        ("particles", "couple_moments_momentum_to_mhd", "true"),
        ("particles", "couple_moments_energy_to_mhd", "true"),
        ("particles", "couple_moments_momentum_coeff", "1.0"),
        ("particles", "couple_moments_energy_coeff", "1.0"),
        ("particles", "pic_boundary_conservation_ledger", "true"),
        ("particles", "pic_physical_mode", "paper_mhd_pic_vl2_tsc"),
        ("particles", "pic_background_mode", "coupled"),
        ("particles", "pic_feedback_mode", "coupled"),
        ("particles", "pic_interp_scheme", "tsc"),
        ("particles", "pic_cr_initial_state", "momentum"),
        ("particles", "pic_cr_hall_mode", "off"),
        ("particles", "pic_wave_damping_mode", "off"),
        ("particles", "pic_deltaf_mode", "off"),
        ("particles", "pic_expanding_box_mode", "off"),
        ("particles", "pic_load_balance_cost_per_particle", "0.0"),
        ("problem", "pgen_name", "pic_parallel_shock"),
        ("problem", "ps_enable_conservation_ledger", "true"),
        ("problem", "ps_enable_injection", "true"),
        ("problem", "ps_enable_gas_subtraction", "true"),
        ("problem", "ps_enable_curvature_amr", "false"),
        ("problem", "ps_enable_frame_tracking", "false"),
        ("problem", "ps_frame_mode", "velocity"),
        ("problem", "ps_inject_t_start", "0.0"),
        ("problem", "ps_remove_birth_time_before", "45.0"),
        ("problem", "user_hist", "true"),
        ("problem", "ps_p0", "1.0"),
    }
    for block, key, expected in required:
        _require(blocks[block].get(key) == expected,
                 f"deck <{block}>/{key} does not equal {expected}")
    try:
        nx1 = int(blocks["mesh"]["nx1"])
        nx2 = int(blocks["mesh"]["nx2"])
        nx3 = int(blocks["mesh"]["nx3"])
    except (KeyError, ValueError) as error:
        raise ConservationClosureError("deck mesh dimensionality is invalid") from error
    _require(nx1 > 1 and nx2 > 1 and nx3 == 1,
             "deck exact closure requires true 2D x1-x2 geometry")
    _require(blocks["mesh"].get("ox1_bc") in {"inflow", "outflow"},
             "deck <mesh>/ox1_bc is not an exact-ledger physical boundary")
    refinement = blocks.get("mesh_refinement", {})
    _require(refinement.get("refinement", "none") == "none",
             "deck exact closure requires a fixed uniform mesh")
    _require(refinement.get("num_levels", "1") == "1",
             "deck exact closure requires one mesh level")
    _require(particles.get("ppc") in {"0", "0.0"},
             "deck must bind an initially empty CR population")
    forbidden_blocks = {
        "hydro", "radiation", "ion_neutral", "turb_driving", "initial_turb",
        "shearing_box", "adm", "z4c", "dyngr",
    }
    _require(not (forbidden_blocks & set(blocks)),
             "deck enables an untracked physics or source block")
    forbidden_mhd_parameters = {
        "viscosity", "ohmic_resistivity", "conductivity", "tdep_conductivity",
        "const_accel", "ism_cooling", "cgm_cooling", "rel_cooling", "beam_source",
        "nscalars", "dfloor", "pfloor", "tfloor", "sfloor",
    }
    _require(not (forbidden_mhd_parameters & set(blocks["mhd"])),
             "deck enables an untracked MHD source, scalar, diffusion, or floor control")
    _require(
        refinement.get("prolong_primitives", "false") == "false",
        "deck enables non-conservative primitive-variable AMR prolongation",
    )
    try:
        qscale = float(particles["deposit_qscale"])
        light_speed = float(particles["pic_cr_light_speed"])
    except (KeyError, ValueError) as exc:
        raise ConservationClosureError("deck particle normalization is missing") from exc
    _require(math.isfinite(qscale) and qscale > 0.0, "deck deposit_qscale is invalid")
    _require(math.isfinite(light_speed) and light_speed > 0.0,
             "deck artificial light speed is invalid")
    species: list[dict[str, float]] = []
    index = 0
    while f"species{index}" in blocks:
        block = blocks[f"species{index}"]
        try:
            mass = float(block["mass"])
            charge = float(block["charge"])
        except (KeyError, ValueError) as exc:
            raise ConservationClosureError(f"deck species{index} is incomplete") from exc
        _require(math.isfinite(mass) and mass > 0.0 and math.isfinite(charge),
                 f"deck species{index} is invalid")
        species.append({"mass": mass, "charge": charge, "q_over_mc": charge / mass})
        index += 1
    _require(species, "deck has no contiguous species table")
    history_outputs = [
        block for name, block in blocks.items()
        if name.startswith("output") and block.get("file_type") == "hst"
    ]
    _require(len(history_outputs) == 1, "deck must define exactly one history output")
    _require(history_outputs[0].get("data_format") == "%24.16e",
             "deck history output must bind %24.16e")
    _require(history_outputs[0].get("dt") == "100.0",
             "deck history output must bind the 100.0 checkpoint cadence")
    _require(history_outputs[0].get("user_hist_only", "false") == "false",
             "deck history output must include both MHD and user histories")
    restart_outputs = [
        block for name, block in blocks.items()
        if name.startswith("output") and block.get("file_type") == "rst"
    ]
    _require(len(restart_outputs) == 1 and restart_outputs[0].get("dt") == "100.0",
             "deck must define one restart output at the 100.0 history cadence")
    try:
        tlim = float(blocks["time"]["tlim"])
        checkpoint_cadence = float(restart_outputs[0]["dt"])
        inject_species = int(problem["ps_inject_species"])
    except (KeyError, ValueError) as exc:
        raise ConservationClosureError(
            "deck checkpoint or injected-species contract is incomplete"
        ) from exc
    _require(
        math.isfinite(tlim)
        and math.isfinite(checkpoint_cadence)
        and tlim > 0.0
        and checkpoint_cadence > 0.0
        and tlim / checkpoint_cadence < 1.0e6,
        "deck checkpoint schedule is invalid",
    )
    _require(0 <= inject_species < len(species), "deck injected species is invalid")
    return {
        "deposit_qscale": qscale,
        "artificial_light_speed": light_speed,
        "species": species,
        "injected_macro_mass": qscale * species[inject_species]["mass"],
        "tlim": tlim,
        "checkpoint_cadence": checkpoint_cadence,
        "expected_checkpoint_times": _expected_checkpoint_times(
            tlim, checkpoint_cadence
        ),
        "input_blocks": blocks,
    }


def _expected_checkpoint_times(tlim: float, cadence: float) -> list[float]:
    times = _expected_cadence_times(tlim, cadence)
    tolerance = 64.0 * math.ulp(max(abs(tlim), abs(cadence), 1.0))
    if not times or abs(times[-1] - tlim) > tolerance:
        times.append(tlim)
    else:
        times[-1] = tlim
    return times


def _expected_cadence_times(tlim: float, cadence: float) -> list[float]:
    tolerance = 64.0 * math.ulp(max(abs(tlim), abs(cadence), 1.0))
    count = int(math.floor((tlim + tolerance) / cadence))
    times = [cadence * index for index in range(1, count + 1)]
    if times and abs(times[-1] - tlim) <= tolerance:
        times[-1] = tlim
    return times


def _validate_checkpoint_sequence(
    observed: Sequence[float],
    expected_slots: Sequence[float],
    *,
    include_initial: bool,
    label: str,
) -> None:
    values = list(observed)
    if include_initial:
        _require(values and values[0] == 0.0, f"{label}: initial time is not zero")
        values = values[1:]
    _require(
        len(values) == len(expected_slots),
        f"{label}: does not contain the complete configured checkpoint sequence",
    )
    for index, (time, slot) in enumerate(zip(values, expected_slots)):
        _require(math.isfinite(time), f"{label}: checkpoint time is nonfinite")
        if index == len(expected_slots) - 1:
            _require(time == slot, f"{label}: terminal checkpoint time drifted")
        else:
            _require(
                time >= slot and time < expected_slots[index + 1],
                f"{label}: checkpoint does not bind its configured nominal slot",
            )


def _parameter_values_equal(expected: str, observed: str) -> bool:
    bool_aliases = {"0": False, "1": True, "false": False, "true": True}
    if expected in bool_aliases and observed in bool_aliases:
        return bool_aliases[expected] is bool_aliases[observed]
    if _REAL.fullmatch(expected) is not None and _REAL.fullmatch(observed) is not None:
        expected_real = float(expected)
        observed_real = float(observed)
        return (
            math.isfinite(expected_real)
            and math.isfinite(observed_real)
            and expected_real == observed_real
        )
    return expected == observed


def _bind_restart_to_deck(
    payload: bytes, deck_blocks: Mapping[str, Mapping[str, str]], label: str
) -> None:
    marker = b"<par_end>\n"
    end = payload.find(marker)
    _require(end >= 0 and payload.find(marker, end + 1) < 0,
             f"{label}: missing or duplicate <par_end> marker")
    observed = _parse_input(
        payload[:end], f"{label}/effective_input", allow_empty_values=True
    )
    for block, parameters in deck_blocks.items():
        if block in {"comment", "job"}:
            continue
        _require(block in observed, f"{label}: effective input is missing <{block}>")
        for key, expected in parameters.items():
            _require(
                key in observed[block]
                and _parameter_values_equal(expected, observed[block][key]),
                f"{label}: restart/deck control drift at <{block}>/{key}",
            )


def _parse_history(
    payload: bytes, expected: Sequence[str], label: str
) -> tuple[list[dict[str, float]], dict[str, int]]:
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise ConservationClosureError(f"{label} is not ASCII") from exc
    header: tuple[str, ...] | None = None
    awaiting_segment_header = False
    segment_headers = 0
    identical_duplicate_times = 0
    rows: list[dict[str, float]] = []
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            if line == "# Athena++ history data":
                _require(not awaiting_segment_header,
                         f"{label}: duplicate history segment banner")
                awaiting_segment_header = True
                continue
            labels = tuple(re.findall(r"\[[0-9]+\]=([^\s]+)", line))
            if labels:
                _require(awaiting_segment_header,
                         f"{label}: history segment header lacks banner")
                _require(labels == tuple(expected), f"{label}: history header drifted")
                header = labels
                segment_headers += 1
                awaiting_segment_header = False
            else:
                _require(False, f"{label}: unrecognized history comment")
            continue
        _require(header is not None and not awaiting_segment_header,
                 f"{label}: data preceded or interrupted a history header")
        fields = line.split()
        _require(len(fields) == len(header), f"{label}: row width drift at line {lineno}")
        values: list[float] = []
        for field in fields:
            _require(_REAL.fullmatch(field) is not None,
                     f"{label}: noncanonical real at line {lineno}")
            parsed = float(field)
            _require(math.isfinite(parsed), f"{label}: nonfinite value at line {lineno}")
            values.append(parsed)
        _require(values[0] >= 0.0 and values[1] > 0.0,
                 f"{label}: invalid time or timestep at line {lineno}")
        if rows:
            if values[0] == rows[-1]["time"]:
                candidate = dict(zip(header, values))
                _require(
                    all(candidate[key] == rows[-1][key] for key in header[2:]),
                    f"{label}: duplicate time has state drift",
                )
                identical_duplicate_times += 1
                continue
            _require(values[0] > rows[-1]["time"],
                     f"{label}: nonmonotonic time")
        rows.append(dict(zip(header, values)))
    _require(header is not None and rows and not awaiting_segment_header,
             f"{label}: missing or incomplete header or rows")
    return rows, {
        "segment_headers": segment_headers,
        "identical_duplicate_times_collapsed": identical_duplicate_times,
    }


def _restart_mhd_state(
    payload: bytes, deck: Mapping[str, object], label: str
) -> tuple[dict[str, object], list[float]]:
    marker = b"<par_end>\n"
    start = payload.find(marker)
    _require(start >= 0 and payload.find(marker, start + 1) < 0,
             f"{label}: missing or duplicate <par_end> marker")
    offset = start + len(marker)
    # Restart schema 7 in this campaign is double precision: RegionSize=9d,
    # RegionIndcs=19i, followed by time, dt, cycle and original rank count.
    minimum = offset + 2 * 4 + 9 * 8 + 2 * 19 * 4 + 2 * 8 + 2 * 4
    _require(minimum <= len(payload), f"{label}: truncated mesh restart header")
    nmb_total, root_level = struct.unpack_from("<ii", payload, offset)
    offset += 8
    mesh_size = struct.unpack_from("<9d", payload, offset)
    offset += 9 * 8
    mesh_indcs = struct.unpack_from(f"<{_REGION_INDCS_COUNT}i", payload, offset)
    offset += _REGION_INDCS_COUNT * 4
    mb_indcs = struct.unpack_from(f"<{_REGION_INDCS_COUNT}i", payload, offset)
    offset += _REGION_INDCS_COUNT * 4
    time, dt, cycle, nranks = struct.unpack_from("<ddii", payload, offset)
    offset += 2 * 8 + 2 * 4
    _require(nmb_total > 0 and root_level >= 0 and nranks > 0 and nranks <= nmb_total,
             f"{label}: invalid mesh/rank header")
    _require(math.isfinite(time) and time >= 0.0 and math.isfinite(dt) and dt > 0.0,
             f"{label}: invalid time header")
    _require(cycle >= 0, f"{label}: invalid cycle header")
    lengths = (
        mesh_size[3] - mesh_size[0],
        mesh_size[4] - mesh_size[1],
        mesh_size[5] - mesh_size[2],
    )
    _require(all(math.isfinite(length) and length > 0.0 for length in lengths),
             f"{label}: invalid mesh extents")
    _require(all(math.isfinite(value) and value > 0.0 for value in mesh_size[6:9]),
             f"{label}: invalid mesh cell widths")

    blocks = deck["input_blocks"]
    mesh = blocks["mesh"]
    meshblock = blocks["meshblock"]
    try:
        deck_lengths = tuple(
            float(mesh[f"x{axis}max"]) - float(mesh[f"x{axis}min"])
            for axis in (1, 2, 3)
        )
        deck_mesh_cells = tuple(int(mesh[f"nx{axis}"]) for axis in (1, 2, 3))
        deck_block_cells = tuple(int(meshblock[f"nx{axis}"]) for axis in (1, 2, 3))
        deck_nghost = int(mesh["nghost"])
    except (KeyError, ValueError) as exc:
        raise ConservationClosureError(f"{label}: deck mesh geometry is incomplete") from exc
    _require(lengths == deck_lengths, f"{label}: restart/deck mesh extents drifted")
    _require(tuple(mesh_indcs[1:4]) == deck_mesh_cells,
             f"{label}: restart/deck mesh cell counts drifted")
    _require(tuple(mb_indcs[1:4]) == deck_block_cells and mb_indcs[0] == deck_nghost,
             f"{label}: restart/deck MeshBlock geometry drifted")

    logical_locations: list[tuple[int, int, int, int]] = []
    for _ in range(nmb_total):
        _require(offset + _LOGICAL_LOCATION_COUNT * 4 <= len(payload),
                 f"{label}: truncated logical-location list")
        location = struct.unpack_from(f"<{_LOGICAL_LOCATION_COUNT}i", payload, offset)
        offset += _LOGICAL_LOCATION_COUNT * 4
        _require(location[3] >= root_level, f"{label}: invalid MeshBlock level")
        logical_locations.append(location)
    _require(offset + nmb_total * 4 <= len(payload),
             f"{label}: truncated MeshBlock cost list")
    costs = struct.unpack_from(f"<{nmb_total}f", payload, offset)
    offset += nmb_total * 4
    _require(all(math.isfinite(value) and value > 0.0 for value in costs),
             f"{label}: invalid MeshBlock cost list")
    _require(all(value == 1.0 for value in costs),
             f"{label}: exact closure rejects non-unit load-balance costs")
    _require(offset + nmb_total * 4 + 2 * nranks * 4 <= len(payload),
             f"{label}: truncated rank layout")
    rank_eachmb = struct.unpack_from(f"<{nmb_total}i", payload, offset)
    offset += nmb_total * 4
    gids_eachrank = struct.unpack_from(f"<{nranks}i", payload, offset)
    offset += nranks * 4
    nmb_eachrank = struct.unpack_from(f"<{nranks}i", payload, offset)
    offset += nranks * 4
    partition_end = 0
    for rank, (gid, count) in enumerate(zip(gids_eachrank, nmb_eachrank)):
        _require(gid == partition_end and count > 0 and count <= nmb_total - partition_end,
                 f"{label}: invalid contiguous rank partition")
        _require(all(value == rank for value in rank_eachmb[gid:gid + count]),
                 f"{label}: MeshBlock rank assignment drifted")
        partition_end += count
    _require(partition_end == nmb_total, f"{label}: rank partition is incomplete")

    _require(offset + 16 <= len(payload), f"{label}: missing mesh metadata")
    mesh_magic, mesh_version, has_refinement_cooldown = struct.unpack_from(
        "<Qii", payload, offset
    )
    offset += 16
    _require(mesh_magic == _MESH_METADATA_MAGIC and
             mesh_version == _MESH_METADATA_VERSION and
             has_refinement_cooldown in {0, 1},
             f"{label}: mesh metadata schema drifted")
    _require(offset + 8 <= len(payload), f"{label}: missing checkpoint nonce")
    checkpoint_nonce = struct.unpack_from("<Q", payload, offset)[0]
    offset += 8
    _require(checkpoint_nonce != 0, f"{label}: invalid checkpoint nonce")
    if has_refinement_cooldown:
        _require(offset + nmb_total * 4 <= len(payload),
                 f"{label}: truncated refinement-cooldown list")
        cooldowns = struct.unpack_from(f"<{nmb_total}i", payload, offset)
        offset += nmb_total * 4
        _require(all(value >= 0 for value in cooldowns),
                 f"{label}: invalid refinement-cooldown list")
    _require(not has_refinement_cooldown,
             f"{label}: exact closure rejects adaptive restart metadata")
    _require(all(location[3] == root_level for location in logical_locations),
             f"{label}: exact closure rejects refined restart topology")
    _validate_fixed_uniform_topology(
        logical_locations=logical_locations,
        deck_mesh_cells=deck_mesh_cells,
        deck_block_cells=deck_block_cells,
        root_level=root_level,
        nmb_total=nmb_total,
        label=label,
    )
    _require(offset + 8 <= len(payload), f"{label}: missing restart data size")
    data_size = struct.unpack_from("<Q", payload, offset)[0]
    offset += 8
    data_offset = offset

    ng, nx1, nx2, nx3, is_, ie, js, je, ks, ke = mb_indcs[:10]
    _require(
        nx1 > 0 and nx2 > 0 and nx3 > 0 and
        ie - is_ + 1 == nx1 and je - js + 1 == nx2 and ke - ks + 1 == nx3,
        f"{label}: invalid active MeshBlock indices",
    )
    nout1 = nx1 + 2 * ng if nx1 > 1 else 1
    nout2 = nx2 + 2 * ng if nx2 > 1 else 1
    nout3 = nx3 + 2 * ng if nx3 > 1 else 1
    mhd_cell_count = _IDEAL_MHD_CONSERVED_COUNT * nout3 * nout2 * nout1
    expected_data_size = 8 * (
        mhd_cell_count
        + nout3 * nout2 * (nout1 + 1)
        + nout3 * (nout2 + 1) * nout1
        + (nout3 + 1) * nout2 * nout1
    )
    _require(data_size == expected_data_size,
             f"{label}: restart contains an untracked or malformed physics payload")
    pic_magic = struct.pack("<Q", restart_layout.PIC_RESTART_MAGIC)
    pic_offset = payload.find(pic_magic)
    _require(pic_offset == data_offset + nmb_total * data_size,
             f"{label}: restart MHD/PIC payload boundary drifted")

    state = np.zeros(_IDEAL_MHD_CONSERVED_COUNT, dtype=np.float64)
    active_dimensions = sum(value > 1 for value in mesh_indcs[1:4])
    _require(active_dimensions > 0, f"{label}: restart mesh has no active dimension")
    for gid, location in enumerate(logical_locations):
        level_scale = 2 ** (location[3] - root_level)
        volume = float(np.prod([
            mesh_size[5 + axis] / level_scale
            if mesh_indcs[axis] > 1 else mesh_size[5 + axis]
            for axis in (1, 2, 3)
        ]))
        _require(math.isfinite(volume) and volume > 0.0,
                 f"{label}: invalid MeshBlock cell volume")
        values = np.frombuffer(
            payload,
            dtype="<f8",
            count=mhd_cell_count,
            offset=data_offset + gid * data_size,
        ).reshape(_IDEAL_MHD_CONSERVED_COUNT, nout3, nout2, nout1)
        active = values[:, ks:ke + 1, js:je + 1, is_:ie + 1]
        _require(np.all(np.isfinite(active)), f"{label}: nonfinite restart MHD state")
        state += np.sum(active, axis=(1, 2, 3)) * volume

    spatial_payload = payload[data_offset:pic_offset]
    return {
        "time": time,
        "dt": dt,
        "cycle": cycle,
        "domain_volume": float(lengths[0] * lengths[1] * lengths[2]),
        "nmb_total": nmb_total,
        "nranks": nranks,
        "has_refinement_cooldown": bool(has_refinement_cooldown),
        "checkpoint_nonce": checkpoint_nonce,
        "spatial_mhd_and_face_field_byte_count": len(spatial_payload),
        "spatial_mhd_and_face_field_sha256": hashlib.sha256(spatial_payload).hexdigest(),
    }, state.astype(float).tolist()


def _parse_int(value: str, label: str) -> int:
    _require(_INT.fullmatch(value) is not None, f"{label}: expected integer")
    return int(value)


def _parse_real(value: str, label: str) -> float:
    _require(_REAL.fullmatch(value) is not None, f"{label}: expected finite real")
    parsed = float(value)
    _require(math.isfinite(parsed), f"{label}: expected finite real")
    return parsed


def _parse_bool(value: str, label: str) -> bool:
    _require(value in {"0", "1", "false", "true"}, f"{label}: expected boolean")
    return value in {"1", "true"}


def _ledger_values_agree(lhs: float, rhs: float, accumulated_terms: int) -> bool:
    scale = max(abs(lhs), abs(rhs), 1.0)
    tolerance = (
        256.0
        * max(accumulated_terms, 1)
        * np.finfo(np.float64).eps
        * scale
    )
    return math.isfinite(tolerance) and abs(lhs - rhs) <= tolerance


def _escape_ledger(
    *,
    parameters: Mapping[str, str],
    startup: Mapping[str, object],
    generic_escape: Sequence[float],
    committed_cycle: int,
    committed_time: float,
    deck: Mapping[str, object],
    label: str,
) -> dict[str, object]:
    fields = set(_ESCAPE_INTEGER_FIELDS + _ESCAPE_BOOLEAN_FIELDS + _ESCAPE_REAL_FIELDS)
    missing = sorted(fields - set(parameters))
    _require(not missing, f"{label}: physical escape ledger is missing {missing!r}")
    values: dict[str, object] = {}
    for field in _ESCAPE_INTEGER_FIELDS:
        values[field] = _parse_int(parameters[field], f"{label}/{field}")
    for field in _ESCAPE_BOOLEAN_FIELDS:
        values[field] = _parse_bool(parameters[field], f"{label}/{field}")
    for field in _ESCAPE_REAL_FIELDS:
        values[field] = _parse_real(parameters[field], f"{label}/{field}")

    _require(values["ps_escape_ledger_schema"] == 1,
             f"{label}: physical escape schema is not 1")
    _require(values["ps_escape_ledger_complete"],
             f"{label}: physical escape ledger is incomplete")
    expected_audit_calls = 2 * committed_cycle
    _require(values["ps_escape_audit_calls"] == expected_audit_calls,
             f"{label}: physical escape audit-call chronology is incomplete")
    _require(
        _ledger_values_agree(
            float(values["ps_escape_last_audit_time"]),
            committed_time,
            expected_audit_calls,
        ),
        f"{label}: physical escape audit time differs from committed time",
    )
    for field in (
        "ps_escaped_injected_cr_count_global",
        "ps_escaped_injected_cr_mass_global",
        "ps_escaped_injected_cr_energy_global",
        "ps_escaped_initial_cr_count_global",
    ):
        _require(float(values[field]) >= 0.0, f"{label}/{field}: negative value")
    for field in (
        "ps_escaped_injected_cr_count_global",
        "ps_escaped_initial_cr_count_global",
    ):
        _require(float(values[field]).is_integer(), f"{label}/{field}: non-integral count")
    _require(values["ps_escaped_initial_cr_count_global"] == 0.0,
             f"{label}: initially seeded CR escape is not fully accounted")
    escaped_count = float(values["ps_escaped_injected_cr_count_global"])
    escaped_mass = float(values["ps_escaped_injected_cr_mass_global"])
    _require(
        _ledger_values_agree(
            escaped_mass,
            escaped_count * float(deck["injected_macro_mass"]),
            expected_audit_calls,
        ),
        f"{label}: escaped injected count/mass relation is inconsistent",
    )
    _require(
        escaped_count + float(startup["ps_removed_cr_count_global"])
        <= float(startup["ps_injected_cr_count_global"]),
        f"{label}: removed plus escaped injected count exceeds injected count",
    )
    _require(
        escaped_mass + float(startup["ps_removed_cr_mass_global"])
        <= float(startup["ps_injected_cr_mass_global"])
        + 256.0 * np.finfo(np.float64).eps
        * max(float(startup["ps_injected_cr_mass_global"]), 1.0),
        f"{label}: removed plus escaped injected mass exceeds injected mass",
    )
    reason_coded = [
        escaped_mass,
        float(values["ps_escaped_injected_cr_momentum_x1_global"]),
        float(values["ps_escaped_injected_cr_momentum_x2_global"]),
        float(values["ps_escaped_injected_cr_momentum_x3_global"]),
        float(values["ps_escaped_injected_cr_energy_global"]),
    ]
    for component, generic, observed in zip(_COMPONENTS, generic_escape, reason_coded):
        _require(
            _ledger_values_agree(generic, -observed, expected_audit_calls),
            f"{label}: generic and reason-coded physical escape differ at {component}",
        )
    values["reason_coded_escape_vector"] = reason_coded
    return values


def _conservation_ledger(
    payload: bytes,
    deck: Mapping[str, object],
    label: str,
    *,
    expected_cycle: int | None = None,
    expected_time: float | None = None,
) -> dict[str, object]:
    try:
        parameters = restart_layout._problem_parameters(payload, label)
        startup = restart_layout.extract_startup_shock_ledger(payload, source=label)
    except restart_layout.RestartPolicyError as exc:
        raise ConservationClosureError(
            f"{label}: restart parameter ledger failed: {exc}"
        ) from exc
    fields = {
        "ps_conservation_ledger_schema",
        "ps_conservation_ledger_complete",
        "ps_conservation_committed_cycles",
        "ps_conservation_committed_time",
    }
    for prefix in _CONS_PREFIXES:
        for component in _COMPONENTS:
            fields.add(f"{prefix}_{component}_global")
    missing = sorted(fields - set(parameters))
    _require(not missing, f"{label}: conservation ledger is missing {missing!r}")
    _require(_parse_int(parameters["ps_conservation_ledger_schema"], label) == 1,
             f"{label}: conservation schema is not 1")
    _require(_parse_bool(parameters["ps_conservation_ledger_complete"], label),
             f"{label}: conservation ledger is incomplete")
    vectors: dict[str, list[float]] = {}
    for prefix in _CONS_PREFIXES:
        vectors[prefix] = [
            _parse_real(parameters[f"{prefix}_{component}_global"],
                        f"{label}/{prefix}/{component}")
            for component in _COMPONENTS
        ]
    _require(vectors["ps_cons_particle_reflect"][0] == 0.0 and
             vectors["ps_cons_particle_reflect"][4] == 0.0,
             f"{label}: reflecting-wall mass or energy delta is nonzero")
    _require(vectors["ps_cons_particle_escape"][0] <= 0.0 and
             vectors["ps_cons_particle_escape"][4] <= 0.0,
             f"{label}: particle-escape mass or energy delta is positive")
    _require(vectors["ps_cons_gas_subtracted"][0] >= 0.0 and
             vectors["ps_cons_gas_subtracted"][4] >= 0.0,
             f"{label}: gas-subtraction mass or energy is negative")
    committed_cycle = _parse_int(parameters["ps_conservation_committed_cycles"], label)
    committed_time = _parse_real(parameters["ps_conservation_committed_time"], label)
    if expected_cycle is not None or expected_time is not None:
        _require(
            committed_cycle == expected_cycle and committed_time == expected_time,
            f"{label}: restart header/ledger commit discontinuity",
        )
    escape = _escape_ledger(
        parameters=parameters,
        startup=startup,
        generic_escape=vectors["ps_cons_particle_escape"],
        committed_cycle=committed_cycle,
        committed_time=committed_time,
        deck=deck,
        label=label,
    )
    return {
        "committed_cycle": committed_cycle,
        "committed_time": committed_time,
        "vectors": vectors,
        "startup": startup,
        "escape": escape,
    }


def _particle_budget(payload: bytes, deck: Mapping[str, object], label: str) -> list[float]:
    try:
        extracted = particle_state.extract_schema7_particle_state(payload, source=label)
    except particle_state.ParticleStateError as exc:
        raise ConservationClosureError(f"{label}: schema-7 particle extraction failed") from exc
    _require(extracted["deposit_qscale"] == deck["deposit_qscale"],
             f"{label}: restart/deck deposit_qscale drift")
    _require(extracted["artificial_light_speed"] == deck["artificial_light_speed"],
             f"{label}: restart/deck artificial light-speed drift")
    species = np.asarray(extracted["species"], dtype=np.int64)
    table = deck["species"]
    _require(np.all((species >= 0) & (species < len(table))),
             f"{label}: particle species is outside deck table")
    masses = np.asarray([entry["mass"] for entry in table], dtype=np.float64)
    qom = np.asarray([entry["q_over_mc"] for entry in table], dtype=np.float64)
    particle_qom = np.asarray(extracted["particle_q_over_mc"], dtype=np.float64)
    _require(np.array_equal(particle_qom, qom[species]),
             f"{label}: particle q/(mc) differs from deck species table")
    momentum = np.asarray(extracted["momentum_per_mass"], dtype=np.float64)
    weights = np.asarray(extracted["macro_weight"], dtype=np.float64)
    _require(np.all(np.isfinite(momentum)) and np.all(np.isfinite(weights)) and
             np.all(weights > 0.0), f"{label}: invalid particle state")
    macro_mass = extracted["deposit_qscale"] * weights * masses[species]
    light_speed = extracted["artificial_light_speed"]
    gamma = np.sqrt(1.0 + np.sum(momentum * momentum, axis=1) / light_speed**2)
    return [
        float(np.sum(macro_mass)),
        *np.sum(macro_mass[:, None] * momentum, axis=0).astype(float).tolist(),
        float(np.sum(macro_mass * (gamma - 1.0) * light_speed**2)),
    ]


def _mhd_budget(row: Mapping[str, float]) -> tuple[list[float], dict[str, float]]:
    state = [row["mass"], row["1-mom"], row["2-mom"], row["3-mom"], row["tot-E"]]
    kinetic = row["1-KE"] + row["2-KE"] + row["3-KE"]
    magnetic = row["1-ME"] + row["2-ME"] + row["3-ME"]
    thermal = row["tot-E"] - kinetic - magnetic
    return state, {
        "kinetic_energy": kinetic,
        "magnetic_energy": magnetic,
        "thermal_energy_by_difference": thermal,
        "total_energy": row["tot-E"],
    }


def _external_delta(ledger: Mapping[str, object]) -> tuple[list[float], dict[str, list[float]]]:
    startup = ledger["startup"]
    terms = dict(ledger["vectors"])
    terms["injected_cr"] = [
        startup["ps_injected_cr_mass_global"],
        startup["ps_injected_cr_momentum_x1_global"],
        startup["ps_injected_cr_momentum_x2_global"],
        startup["ps_injected_cr_momentum_x3_global"],
        startup["ps_injected_cr_energy_global"],
    ]
    terms["removed_cr"] = [
        startup["ps_removed_cr_mass_global"],
        startup["ps_removed_cr_momentum_x1_global"],
        startup["ps_removed_cr_momentum_x2_global"],
        startup["ps_removed_cr_momentum_x3_global"],
        startup["ps_removed_cr_energy_global"],
    ]
    external = []
    for n in range(5):
        external.append(
            terms["ps_cons_mhd_boundary"][n]
            + terms["ps_cons_particle_reflect"][n]
            + terms["ps_cons_particle_escape"][n]
            + terms["injected_cr"][n]
            - terms["ps_cons_gas_subtracted"][n]
            - terms["removed_cr"][n]
        )
    return external, terms


def _normalized(residual: Sequence[float], denominator: Sequence[float]) -> list[float | None]:
    return [
        None if scale == 0.0 else value / scale
        for value, scale in zip(residual, denominator)
    ]


def _checkpoint_closure(
    *,
    time: float,
    cycle: int,
    initial: Sequence[float],
    mhd_history_state: Sequence[float],
    mhd_restart_state: Sequence[float],
    mhd_partition: Mapping[str, float],
    cr_state: Sequence[float],
    external: Sequence[float],
    terms: Mapping[str, Sequence[float]],
) -> dict[str, object]:
    history_actual = [mhd_history_state[n] + cr_state[n] for n in range(5)]
    actual = [mhd_restart_state[n] + cr_state[n] for n in range(5)]
    mhd_post_history_delta = [
        mhd_restart_state[n] - mhd_history_state[n] for n in range(5)
    ]
    expected = [initial[n] + external[n] for n in range(5)]
    history_signed_residual = [
        history_actual[n] - expected[n] for n in range(5)
    ]
    signed_residual = [actual[n] - expected[n] for n in range(5)]
    absolute_residual = [abs(value) for value in signed_residual]
    transport_l1 = [
        sum(abs(term[n]) for term in terms.values())
        for n in range(5)
    ]
    signed_contributions = {
        name: [
            (-value if name in {"ps_cons_gas_subtracted", "removed_cr"} else value)
            for value in term
        ]
        for name, term in terms.items()
    }
    denominators = {
        "initial_state_abs": [abs(value) for value in initial],
        "expected_final_state_abs": [abs(value) for value in expected],
        "max_state_abs": [
            max(abs(initial[n]), abs(actual[n]), abs(expected[n])) for n in range(5)
        ],
        "cumulative_transport_l1": transport_l1,
    }
    signed_normalized = {
        key: _normalized(signed_residual, value) for key, value in denominators.items()
    }
    absolute_normalized = {
        key: _normalized(absolute_residual, value) for key, value in denominators.items()
    }
    momentum_residual = signed_residual[1:4]
    momentum_denominators = {
        "initial_state_l2": float(np.linalg.norm(initial[1:4])),
        "expected_final_state_l2": float(np.linalg.norm(expected[1:4])),
        "max_state_l2": max(
            float(np.linalg.norm(initial[1:4])),
            float(np.linalg.norm(actual[1:4])),
            float(np.linalg.norm(expected[1:4])),
        ),
        "cumulative_transport_l2": float(np.linalg.norm(transport_l1[1:4])),
    }
    momentum_l2 = float(np.linalg.norm(momentum_residual))
    return {
        "cycle": cycle,
        "time": time,
        "instantaneous_partitions": {
            "is_conservation_closure": False,
            "mhd_history_pre_restart": dict(mhd_partition),
            "mhd_restart_state": dict(zip(_COMPONENTS, mhd_restart_state)),
            "cr": {
                "mass": cr_state[0],
                "momentum": list(cr_state[1:4]),
                "kinetic_energy": cr_state[4],
            },
            "combined_state": dict(zip(_COMPONENTS, actual)),
        },
        "conservation_closure": {
            "is_instantaneous_partition": False,
            "initial_combined_state": dict(zip(_COMPONENTS, initial)),
            "cumulative_ledger_terms": {
                key: dict(zip(_COMPONENTS, value)) for key, value in terms.items()
            },
            "signed_external_contributions": {
                key: dict(zip(_COMPONENTS, value))
                for key, value in signed_contributions.items()
            },
            "expected_external_delta": dict(zip(_COMPONENTS, external)),
            "expected_combined_state": dict(zip(_COMPONENTS, expected)),
            "history_combined_state": dict(zip(_COMPONENTS, history_actual)),
            "history_signed_residual": dict(
                zip(_COMPONENTS, history_signed_residual)
            ),
            "mhd_state_change_after_history_before_restart": dict(
                zip(_COMPONENTS, mhd_post_history_delta)
            ),
            "actual_combined_state": dict(zip(_COMPONENTS, actual)),
            "signed_residual": dict(zip(_COMPONENTS, signed_residual)),
            "absolute_residual": dict(zip(_COMPONENTS, absolute_residual)),
            "denominators": {
                key: dict(zip(_COMPONENTS, value)) for key, value in denominators.items()
            },
            "signed_normalized_residual": {
                key: dict(zip(_COMPONENTS, value))
                for key, value in signed_normalized.items()
            },
            "absolute_normalized_residual": {
                key: dict(zip(_COMPONENTS, value))
                for key, value in absolute_normalized.items()
            },
            "momentum_l2": {
                "absolute_residual": momentum_l2,
                "denominators": momentum_denominators,
                "normalized_residual": {
                    key: None if value == 0.0 else momentum_l2 / value
                    for key, value in momentum_denominators.items()
                },
            },
        },
    }


def reduce_exact_conservation_closure(
    *,
    lineage: object,
    candidate_artifact: object,
    deck_artifact: object,
    attempt_artifact: object,
    mhd_history_artifact: object,
    user_history_artifact: object,
    restart_artifacts: object,
) -> dict[str, Any]:
    """Bind supplied artifacts and compute exact combined MHD+CR closure residuals."""
    bound_lineage = _lineage(lineage)
    candidate = _artifact(candidate_artifact, bound_lineage, "candidate_artifact")
    deck_art = _artifact(deck_artifact, bound_lineage, "deck_artifact")
    attempt = _artifact(attempt_artifact, bound_lineage, "attempt_artifact")
    mhd_art = _artifact(mhd_history_artifact, bound_lineage, "mhd_history_artifact")
    user_art = _artifact(user_history_artifact, bound_lineage, "user_history_artifact")
    _require(candidate["sha256"] == bound_lineage["candidate_sha256"],
             "candidate artifact differs from lineage")
    _require(deck_art["sha256"] == bound_lineage["deck_sha256"],
             "deck artifact differs from lineage")
    _require(attempt["sha256"] == bound_lineage["attempt_sha256"],
             "attempt artifact differs from lineage")
    _require(mhd_art["name"].endswith(".mhd.hst"), "MHD history name drifted")
    _require(user_art["name"].endswith(".user.hst"), "user history name drifted")
    _require(type(restart_artifacts) is list and len(restart_artifacts) >= 2,
             "at least two restart checkpoints are required")
    restarts = [
        _artifact(value, bound_lineage, f"restart_artifacts[{index}]")
        for index, value in enumerate(restart_artifacts)
    ]
    _require(all(item["name"].endswith(".rst") for item in restarts),
             "restart artifact name drifted")

    deck = _deck_physics(deck_art["payload"])
    mhd_rows, mhd_history_audit = _parse_history(
        mhd_art["payload"], _MHD_LABELS, "mhd_history"
    )
    user_rows, user_history_audit = _parse_history(
        user_art["payload"], _USER_LABELS, "user_history"
    )
    _require([row["time"] for row in mhd_rows] == [row["time"] for row in user_rows],
             "MHD and user history time grids differ")
    _validate_checkpoint_sequence(
        [row["time"] for row in mhd_rows],
        deck["expected_checkpoint_times"],
        include_initial=True,
        label="history",
    )
    _require(mhd_rows[0]["time"] == 0.0, "history lacks the t=0 initial state")
    _require(all(user_rows[0][key] == 0.0 for key in _USER_LABELS[2:]),
             "initial user conservation history is nonzero")
    initial_mhd, _ = _mhd_budget(mhd_rows[0])
    initial_combined = initial_mhd
    mhd_by_time = {row["time"]: row for row in mhd_rows}
    user_by_time = {row["time"]: row for row in user_rows}

    checkpoints: list[dict[str, object]] = []
    previous_cycle = -1
    previous_time = -1.0
    previous_ledger: Mapping[str, object] | None = None
    domain_volume: float | None = None
    observed_restart_times: list[float] = []
    for index, artifact in enumerate(restarts):
        label = f"restart_artifacts[{index}]"
        payload = artifact["payload"]
        _bind_restart_to_deck(payload, deck["input_blocks"], label)
        try:
            restart_layout.probe_schema7_restart_payload(payload, source=label)
        except restart_layout.RestartPolicyError as exc:
            raise ConservationClosureError(f"{label}: schema-7 probe failed") from exc
        header, restart_mhd = _restart_mhd_state(payload, deck, label)
        ledger = _conservation_ledger(
            payload,
            deck,
            label,
            expected_cycle=header["cycle"],
            expected_time=header["time"],
        )
        _require(header["cycle"] > previous_cycle and header["time"] > previous_time,
                 f"{label}: restart sequence is duplicate or nonmonotonic")
        _require(header["time"] in mhd_by_time, f"{label}: MHD history time is missing")
        _require(header["time"] in user_by_time, f"{label}: user history time is missing")
        observed_restart_times.append(header["time"])
        if domain_volume is None:
            domain_volume = header["domain_volume"]
        _require(header["domain_volume"] == domain_volume,
                 f"{label}: restart domain volume drifted")

        if previous_ledger is not None:
            current_startup = ledger["startup"]
            prior_startup = previous_ledger["startup"]
            for field in (
                "ps_injected_cr_count_global", "ps_injected_cr_mass_global",
                "ps_injected_cr_energy_global", "ps_removed_cr_count_global",
                "ps_removed_cr_mass_global", "ps_removed_cr_energy_global",
            ):
                _require(current_startup[field] >= prior_startup[field],
                         f"{label}: startup ledger regressed at {field}")
            for n in (0, 4):
                _require(
                    ledger["vectors"]["ps_cons_gas_subtracted"][n]
                    >= previous_ledger["vectors"]["ps_cons_gas_subtracted"][n],
                    f"{label}: gas-subtraction ledger regressed",
                )
                _require(
                    ledger["vectors"]["ps_cons_particle_escape"][n]
                    <= previous_ledger["vectors"]["ps_cons_particle_escape"][n],
                    f"{label}: particle-escape ledger regressed",
                )
            for field in (
                "ps_escape_audit_calls",
                "ps_escape_last_audit_time",
                "ps_escaped_injected_cr_count_global",
                "ps_escaped_injected_cr_mass_global",
                "ps_escaped_injected_cr_energy_global",
            ):
                _require(
                    ledger["escape"][field] >= previous_ledger["escape"][field],
                    f"{label}: reason-coded escape ledger regressed at {field}",
                )

        external, terms = _external_delta(ledger)
        user = user_by_time[header["time"]]
        observed_external = [
            user["ext_mass"], user["ext_mom1"], user["ext_mom2"],
            user["ext_mom3"], user["ext_etot"],
        ]
        _require(observed_external == external,
                 f"{label}: user-history external delta differs from restart ledger")
        _require(
            [user["mhd_bmass"], user["mhd_bmom1"], user["mhd_bmom2"],
             user["mhd_bmom3"], user["mhd_betot"]]
            == ledger["vectors"]["ps_cons_mhd_boundary"],
            f"{label}: user-history MHD boundary delta differs from restart ledger",
        )
        _require(
            user["cr_bmass"] ==
            ledger["vectors"]["ps_cons_particle_reflect"][0]
            + ledger["vectors"]["ps_cons_particle_escape"][0]
            and user["cr_betot"] ==
            ledger["vectors"]["ps_cons_particle_reflect"][4]
            + ledger["vectors"]["ps_cons_particle_escape"][4],
            f"{label}: user-history CR boundary delta differs from restart ledger",
        )
        cr = _particle_budget(payload, deck, label)
        mhd, partition = _mhd_budget(mhd_by_time[header["time"]])
        checkpoints.append(
            _checkpoint_closure(
                time=header["time"],
                cycle=header["cycle"],
                initial=initial_combined,
                mhd_history_state=mhd,
                mhd_restart_state=restart_mhd,
                mhd_partition=partition,
                cr_state=cr,
                external=external,
                terms=terms,
            )
        )
        previous_cycle = header["cycle"]
        previous_time = header["time"]
        previous_ledger = ledger

    _validate_checkpoint_sequence(
        observed_restart_times,
        deck["expected_checkpoint_times"],
        include_initial=False,
        label="restart artifacts",
    )

    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "status":
            "supplied_artifact_self_consistency_reduction_only_not_execution_provenance",
        "qualification_effect": "none",
        "authority": {
            "execution_authorized": False,
            "launch_authorized": False,
            "policy_mutation_authorized": False,
            "scientific_claim_authorized": False,
            "publication_authorized": False,
        },
        "lineage": bound_lineage,
        "artifact_bindings": {
            "candidate": _artifact_record(candidate),
            "deck": _artifact_record(deck_art),
            "attempt": _artifact_record(attempt),
            "mhd_history": _artifact_record(mhd_art),
            "user_history": _artifact_record(user_art),
            "restarts": [_artifact_record(item) for item in restarts],
        },
        "history_publication_audit": {
            "mhd": mhd_history_audit,
            "user": user_history_audit,
            "duplicate_time_policy":
                "collapse_only_when_all_state_columns_are_exactly_identical",
        },
        "restart_deck_control_binding":
            "all_nominal_deck_parameters_except_comment_and_job_match_each_restart",
        "checkpoint_timing_contract": {
            "mhd_history_state": "end_of_cycle_before_restart_output",
            "mhd_and_particle_restart_state": "same_committed_time_fixed_uniform_topology",
            "interpretation":
                "closure_uses_restart_state_and_separately_reports_mhd_state_change_after_history",
        },
        "domain_volume": domain_volume,
        "initial_state_contract": {
            "time": 0.0,
            "initial_cr_population": "empty_bound_by_deck_ppc_zero",
            "combined_state": dict(zip(_COMPONENTS, initial_combined)),
        },
        "checkpoints": checkpoints,
        "limitations": [
            "This pure reducer accepts supplied bytes only and grants no execution authority.",
            "It does not establish build, scheduler, command, or independent-reducer execution provenance.",
            "Instantaneous energy partitions are diagnostics and are not conservation closure.",
            "History and restart state at the same committed time are reported separately.",
            "The numerical candidate is limited to the fixed-uniform Section 5.4 boundary geometry.",
            "AMR, SMR, and runtime load balancing are rejected by this bounded successor.",
            "Other post-history MHD state change is measured directly from restart state.",
        ],
    }
