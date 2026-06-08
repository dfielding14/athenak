#!/usr/bin/env python3
"""Pure Q019 raw-artifact reduction for registered nonlinear Bell runs.

This module accepts bytes that a separate trusted boundary has already retained.
It performs no discovery, scheduling, receipt validation, publication, or
authorization.  Its output is suitable for the Q019 physics analyzer only after
the installed control plane independently binds the exact byte inventory.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import PurePosixPath
from typing import Mapping, Sequence

import numpy as np

from tst.publication import analyze_q011_section54_outputs as binary
from tst.publication import q019_nonlinear_bell_particle_state as particle_reducer
from tst.publication import q019_particle_state_analysis_bridge_v2 as particle_bridge
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as decks


SCHEMA_VERSION = 1
RECORD_TYPE = "q019_registered_raw_reduction_v1"
MHD_FIELDS = (
    "dens",
    "eint",
    "velx",
    "vely",
    "velz",
    "bcc1",
    "bcc2",
    "bcc3",
)
MOMENT_FIELDS = (
    "prtcl_rho",
    "prtcl_jx",
    "prtcl_jy",
    "prtcl_jz",
    "prtcl_dedt",
    "prtcl_dpxdt",
    "prtcl_dpydt",
    "prtcl_dpzdt",
    "prtcl_ebdot",
)
REQUIRED_BINARY_PRODUCTS = ("mhd_w_bcc", *MOMENT_FIELDS)
_MUTABLE_OUTPUT_PARAMETERS = frozenset({"file_number", "last_time"})


class RawReductionError(ValueError):
    """Raised when retained Q019 raw bytes cannot be reduced exactly."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RawReductionError(message)


def _case(case_id: str) -> dict[str, object]:
    matches = [row for row in decks.expected_cases() if row["case_id"] == case_id]
    _require(len(matches) == 1, "Q019 raw reduction case identity is unknown")
    return dict(matches[0])


def _canonical_nonnegative_integer(value: str, *, label: str) -> None:
    try:
        parsed = int(value)
    except ValueError as error:
        raise RawReductionError(f"{label} is not an integer") from error
    _require(parsed >= 0 and value == str(parsed), f"{label} is not canonical")


def _finite_output_last_time(value: str, *, label: str) -> None:
    try:
        parsed = float(value)
    except ValueError as error:
        raise RawReductionError(f"{label} is not numeric") from error
    _require(
        math.isfinite(parsed) and (parsed == -1.0 or parsed >= 0.0),
        f"{label} is neither the unwritten sentinel nor a nonnegative time",
    )


def _normalized_runtime_parameters(
    parameters: Mapping[str, Mapping[str, str]],
) -> dict[str, dict[str, str]]:
    normalized: dict[str, dict[str, str]] = {}
    for block_name, block in parameters.items():
        _require(
            isinstance(block_name, str) and isinstance(block, Mapping),
            "Q019 raw runtime parameter block is malformed",
        )
        normalized_block: dict[str, str] = {}
        for name, value in block.items():
            _require(
                isinstance(name, str) and isinstance(value, str),
                "Q019 raw runtime parameter value is malformed",
            )
            if block_name.startswith("output") and name in _MUTABLE_OUTPUT_PARAMETERS:
                if name == "file_number":
                    _canonical_nonnegative_integer(
                        value, label=f"{block_name}/{name}"
                    )
                else:
                    _finite_output_last_time(value, label=f"{block_name}/{name}")
                continue
            normalized_block[name] = value
        normalized[block_name] = normalized_block
    return normalized


def _parse_products(
    products: Mapping[str, bytes],
) -> dict[str, binary.AthenaBinaryDataset]:
    _require(
        type(products) is dict and set(products) == set(REQUIRED_BINARY_PRODUCTS),
        "Q019 raw binary product inventory drifted",
    )
    parsed: dict[str, binary.AthenaBinaryDataset] = {}
    for product in REQUIRED_BINARY_PRODUCTS:
        payload = products[product]
        _require(type(payload) is bytes and payload, f"{product} payload is invalid")
        try:
            parsed[product] = binary.parse_athenak_binary_bytes(
                payload, source=f"<retained:{product}>"
            )
        except binary.AnalysisError as error:
            raise RawReductionError(
                f"{product} failed strict AthenaK binary parsing"
            ) from error
    return parsed


def _validate_cross_product_metadata(
    datasets: Mapping[str, binary.AthenaBinaryDataset],
    *,
    case: Mapping[str, object],
) -> binary.AthenaBinaryDataset:
    reference = datasets["mhd_w_bcc"]
    _require(
        len(reference.variable_names) == len(MHD_FIELDS)
        and set(reference.variable_names) == set(MHD_FIELDS),
        "Q019 mhd_w_bcc field inventory drifted",
    )
    reference_parameters = _normalized_runtime_parameters(
        reference.input_parameters
    )
    _require(
        decks.deck_semantics_payload(reference_parameters)
        == decks.matrix_identity_payload(case),
        "Q019 raw runtime parameters differ from the immutable matrix row",
    )
    for product in MOMENT_FIELDS:
        dataset = datasets[product]
        _require(
            dataset.variable_names == (product,),
            f"Q019 {product} field inventory drifted",
        )
        _require(
            dataset.cycle == reference.cycle
            and dataset.time == reference.time
            and dataset.location_size == reference.location_size
            and dataset.variable_size == reference.variable_size
            and dataset.root_grid_shape == reference.root_grid_shape
            and dataset.meshblock_shape == reference.meshblock_shape
            and dataset.nghost == reference.nghost
            and dataset.domain_bounds == reference.domain_bounds
            and _normalized_runtime_parameters(dataset.input_parameters)
            == reference_parameters,
            f"Q019 {product} metadata differs from mhd_w_bcc",
        )
    return reference


def _composite_fields(
    datasets: Mapping[str, binary.AthenaBinaryDataset],
) -> tuple[dict[str, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    reference_grid = binary.compose_leaf_field(datasets["mhd_w_bcc"], MHD_FIELDS[0])
    fields: dict[str, np.ndarray] = {}
    for field in MHD_FIELDS:
        product = "mhd_w_bcc"
        composite = binary.compose_leaf_field(
            datasets[product], field, target_level=reference_grid.target_level
        )
        _require(
            np.array_equal(composite.source_levels, reference_grid.source_levels)
            and np.array_equal(composite.x1_faces, reference_grid.x1_faces)
            and np.array_equal(composite.x2_faces, reference_grid.x2_faces)
            and np.array_equal(composite.x3_faces, reference_grid.x3_faces),
            f"Q019 {field} composite geometry drifted",
        )
        fields[field] = np.asarray(composite.values, dtype=np.float64)
    for field in MOMENT_FIELDS:
        composite = binary.compose_leaf_field(
            datasets[field], field, target_level=reference_grid.target_level
        )
        _require(
            np.array_equal(composite.source_levels, reference_grid.source_levels)
            and np.array_equal(composite.x1_faces, reference_grid.x1_faces)
            and np.array_equal(composite.x2_faces, reference_grid.x2_faces)
            and np.array_equal(composite.x3_faces, reference_grid.x3_faces),
            f"Q019 {field} composite geometry drifted",
        )
        fields[field] = np.asarray(composite.values, dtype=np.float64)
    return fields, (
        reference_grid.x1_faces,
        reference_grid.x2_faces,
        reference_grid.x3_faces,
    )


def compose_snapshot(
    case_id: str,
    products: Mapping[str, bytes],
) -> dict[str, object]:
    """Compose one exact matched ten-product snapshot."""
    case = _case(case_id)
    datasets = _parse_products(products)
    reference = _validate_cross_product_metadata(datasets, case=case)
    fields, faces = _composite_fields(datasets)
    _require(np.all(fields["dens"] > 0.0), "Q019 density is not positive")
    _require(np.all(fields["eint"] >= 0.0), "Q019 internal energy is negative")
    return {
        "cycle": reference.cycle,
        "time": float(reference.time),
        "x1_faces": faces[0],
        "x2_faces": faces[1],
        "x3_faces": faces[2],
        "fields": fields,
    }


def _cell_volumes(snapshot: Mapping[str, object]) -> np.ndarray:
    return (
        np.diff(np.asarray(snapshot["x3_faces"], dtype=np.float64))[:, None, None]
        * np.diff(np.asarray(snapshot["x2_faces"], dtype=np.float64))[None, :, None]
        * np.diff(np.asarray(snapshot["x1_faces"], dtype=np.float64))[None, None, :]
    )


def mhd_budget(snapshot: Mapping[str, object]) -> dict[str, object]:
    """Integrate the exact MHD momentum and energy carried by one snapshot."""
    fields = snapshot["fields"]
    _require(type(fields) is dict, "Q019 snapshot fields are malformed")
    volume = _cell_volumes(snapshot)
    density = np.asarray(fields["dens"], dtype=np.float64)
    velocity = np.stack(
        [np.asarray(fields[name], dtype=np.float64) for name in ("velx", "vely", "velz")]
    )
    magnetic = np.stack(
        [np.asarray(fields[name], dtype=np.float64) for name in ("bcc1", "bcc2", "bcc3")]
    )
    momentum = np.sum(density[None, ...] * velocity * volume[None, ...], axis=(1, 2, 3))
    kinetic = float(
        np.sum(0.5 * density * np.sum(velocity * velocity, axis=0) * volume)
    )
    thermal = float(np.sum(np.asarray(fields["eint"]) * volume))
    magnetic_energy = float(
        np.sum(0.5 * np.sum(magnetic * magnetic, axis=0) * volume)
    )
    return {
        "gas_momentum": momentum.tolist(),
        "gas_kinetic_energy": kinetic,
        "gas_thermal_energy": thermal,
        "magnetic_energy": magnetic_energy,
        "mhd_total_energy": kinetic + thermal + magnetic_energy,
    }


def _gas_bulk_velocity(snapshot: Mapping[str, object]) -> np.ndarray:
    fields = snapshot["fields"]
    volume = _cell_volumes(snapshot)
    mass = np.asarray(fields["dens"], dtype=np.float64) * volume
    total_mass = float(np.sum(mass))
    _require(total_mass > 0.0, "Q019 gas mass is not positive")
    return np.asarray(
        [
            np.sum(mass * np.asarray(fields[name], dtype=np.float64)) / total_mass
            for name in ("velx", "vely", "velz")
        ]
    )


def _restart_binding(path: str, payload: bytes) -> dict[str, str]:
    relative = PurePosixPath(path)
    _require(
        not relative.is_absolute()
        and ".." not in relative.parts
        and relative.as_posix() == path
        and path.endswith(".rst"),
        "Q019 restart binding path is malformed",
    )
    return {"path": path, "sha256": hashlib.sha256(payload).hexdigest()}


def _reference_budget(reduction: Mapping[str, object]) -> dict[str, object]:
    conservation = reduction["conservation"]
    _require(type(conservation) is dict, "Q019 particle conservation record is malformed")
    return {
        "total_momentum": list(conservation["gas_plus_cr_total_momentum"]),
        "total_energy": float(conservation["gas_plus_cr_total_energy"]),
    }


def _reduce_checkpoint(
    *,
    case: Mapping[str, object],
    checkpoint: Mapping[str, object],
    reference_budget: Mapping[str, object] | None,
) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    _require(
        type(checkpoint) is dict
        and set(checkpoint)
        == {"cycle", "time", "binary_products", "restart_path", "restart_payload"},
        "Q019 checkpoint schema drifted",
    )
    cycle = checkpoint["cycle"]
    time = checkpoint["time"]
    restart_payload = checkpoint["restart_payload"]
    _require(type(cycle) is int and cycle >= 0, "Q019 checkpoint cycle is invalid")
    _require(
        type(time) is float and math.isfinite(time) and time >= 0.0,
        "Q019 checkpoint time is invalid",
    )
    _require(
        type(restart_payload) is bytes and restart_payload,
        "Q019 checkpoint restart payload is invalid",
    )
    snapshot = compose_snapshot(str(case["case_id"]), checkpoint["binary_products"])
    _require(
        snapshot["cycle"] == cycle and snapshot["time"] == time,
        "Q019 checkpoint projection differs from embedded binary cycle/time",
    )
    budget = mhd_budget(snapshot)
    reduction = particle_reducer.reduce_schema7_restart_payload(
        restart_payload,
        source=str(checkpoint["restart_path"]),
        species_mass=[float(case["species_mass"])],
        species_q_over_mc=[
            float(case["species_charge"]) / float(case["species_mass"])
        ],
        domain_volume=math.prod(float(value) for value in case["extents"]),
        gas_bulk_velocity=_gas_bulk_velocity(snapshot),
        guide_field_direction=[1.0, 0.0, 0.0],
        mhd_budget=budget,
        reference_budget=reference_budget,
    )
    particle_state = particle_bridge.build_record(
        reduction,
        case_id=str(case["case_id"]),
        campaign_id=str(case["campaign_id"]),
        cycle=cycle,
        time=time,
        raw_restart_binding=_restart_binding(
            str(checkpoint["restart_path"]), restart_payload
        ),
    )
    return snapshot, particle_state, reduction


def reduce_matched_checkpoints(
    case_id: str,
    checkpoints: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Reduce a strict cycle-zero-starting sequence of matched checkpoints."""
    case = _case(case_id)
    _require(
        isinstance(checkpoints, Sequence)
        and not isinstance(checkpoints, (str, bytes))
        and len(checkpoints) >= 2,
        "Q019 matched checkpoint inventory is incomplete",
    )
    cycles = [checkpoint.get("cycle") for checkpoint in checkpoints]
    times = [checkpoint.get("time") for checkpoint in checkpoints]
    _require(
        cycles[0] == 0
        and times[0] == 0.0
        and all(
            type(cycles[index]) is int
            and cycles[index] > cycles[index - 1]
            and type(times[index]) is float
            and times[index] > times[index - 1]
            for index in range(1, len(checkpoints))
        ),
        "Q019 matched checkpoint chronology drifted",
    )

    _, _, initial = _reduce_checkpoint(
        case=case, checkpoint=checkpoints[0], reference_budget=None
    )
    reference = _reference_budget(initial)
    snapshots: list[dict[str, object]] = []
    particle_states: list[dict[str, object]] = []
    reductions: list[dict[str, object]] = []
    for checkpoint in checkpoints:
        snapshot, particle_state, reduction = _reduce_checkpoint(
            case=case,
            checkpoint=checkpoint,
            reference_budget=reference,
        )
        snapshots.append(snapshot)
        particle_states.append(particle_state)
        reductions.append(reduction)
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "case_id": case_id,
        "campaign_id": case["campaign_id"],
        "matched_checkpoint_count": len(checkpoints),
        "chronology": [
            {"cycle": int(cycle), "time": float(time)}
            for cycle, time in zip(cycles, times)
        ],
        "reference_budget": reference,
        "snapshots": snapshots,
        "particle_states": particle_states,
        "particle_reductions": reductions,
        "authority": {
            "launch_authorized": False,
            "policy_authorized": False,
            "qualification_authorized": False,
            "claim_authorized": False,
            "raw_production_authorized": False,
            "nonlinear_saturation_claim_authorized": False,
        },
        "limitations": [
            "pure byte reduction only",
            "retained-artifact provenance and scheduler evidence are external",
            "checkpoint cycle/time projections require installed-control-plane binding",
        ],
    }


__all__ = [
    "MHD_FIELDS",
    "MOMENT_FIELDS",
    "REQUIRED_BINARY_PRODUCTS",
    "RawReductionError",
    "compose_snapshot",
    "mhd_budget",
    "reduce_matched_checkpoints",
]
