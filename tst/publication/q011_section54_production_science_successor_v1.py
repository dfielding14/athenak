#!/usr/bin/env python3
"""Pure source-local Q-011 production-science diagnostics successor.

This module consumes already decoded AthenaK ``mhd_w_bcc``, deposited-current,
``mhd_j2``, and particle snapshots.  It does not launch work, mutate policy,
bind campaign artifacts, or close a scientific claim.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import wraps
import json
import math
from numbers import Real
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

if __package__:
    from . import analyze_q011_section54_outputs as output_primitives
    from . import q011_section54_particles as particle_primitives
    from . import q011_section54_model as frozen_model
else:
    import analyze_q011_section54_outputs as output_primitives
    import q011_section54_particles as particle_primitives
    import q011_section54_model as frozen_model


SCHEMA_VERSION = 1
SUCCESSOR_ID = "q011_section54_production_science_diagnostic_successor_v1"
REPO_ROOT = Path(__file__).resolve().parents[2]
BASELINE_DECK = (
    REPO_ROOT / "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
)
SUCCESSOR_DECK = (
    REPO_ROOT
    / "inputs/publication/"
    "pic_parallel_shock_section54_production_science_successor_v1_vl2_tsc.athinput"
)

MHD_PRIMITIVE_FIELDS = (
    "dens",
    "velx",
    "vely",
    "velz",
    "eint",
    "bcc1",
    "bcc2",
    "bcc3",
)
CURRENT_PRODUCT_FIELDS: Mapping[str, str] = MappingProxyType(
    {
        "prtcl_rho": "prtcl_rho",
        "prtcl_jx": "prtcl_jx",
        "prtcl_jy": "prtcl_jy",
        "prtcl_jz": "prtcl_jz",
        "mhd_j2": "j2",
    }
)
DERIVED_PROFILE_FIELDS = ("pressure", "bmag")
REFERENCE_B0 = 1.0
GAS_GAMMA = 5.0 / 3.0
PARTICLE_LIGHT_SPEED = frozen_model.LIGHT_SPEED
UPSTREAM_SPEED_U0 = frozen_model.UPSTREAM_SPEED_U0
PARTICLE_MACRO_MASS = 9.0e-4
SHOCK_INJECTED_SOURCE = 1
BIRTH_TIME_MINIMUM = 45.0
SHOCK_SEARCH_OFFSETS_C_OVER_OMEGA_PI = (-1200.0, 1200.0)
MAX_FRONT_OFFSET_C_OVER_OMEGA_PI = 600.0
DOWNSTREAM_STATE_OFFSETS_C_OVER_OMEGA_PI = (-1200.0, -120.0)
UPSTREAM_STATE_OFFSETS_C_OVER_OMEGA_PI = (120.0, 1200.0)
SHOCK_ENERGY_OFFSETS_C_OVER_OMEGA_PI = (-1200.0, 1200.0)
UPSTREAM_B_OFFSETS_C_OVER_OMEGA_PI = (120.0, 1200.0)
UPSTREAM_B0_QUALIFYING_GATE = (1.2, 3.5)
MAGNETIC_QUANTILES = (0.5, 0.9, 0.99)
MAGNETIC_AREA_FRACTION_THRESHOLDS_OVER_B0 = (2.0, 3.0, 4.0)
PARTICLE_ENERGY_QUANTILES = (0.99, 0.999)
PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES = 1000
DOWNSTREAM_SPECTRUM_MINIMUM_PARTICLE_SAMPLES = 1000
PREREGISTERED_SPECTRUM_NOMINAL_TIMES = (500.0, 1200.0)
HISTORY_LINEAR_FIT_MINIMUM_SNAPSHOTS = 3
MATCHED_SNAPSHOT_FIRST_ELIGIBLE_NOMINAL_TIME = 200.0
_MUTABLE_OUTPUT_PARAMETER_NAMES = frozenset({"file_number", "last_time"})

_SUCCESSOR_BLOCK = MappingProxyType(
    {
        "contract_id": SUCCESSOR_ID,
        "deck_role": "source_local_diagnostic_successor_not_authorized",
        "qualification_effect": "none_no_launch_no_policy_authorization_no_claim_closure",
        "full_state_output": "mhd_w_bcc",
        "particle_energy_scope": "shock_injected_birth_time_ge_45",
        "shock_front_gradient": "negative",
        "shock_front_preregistration_status": "requires_supersession",
    }
)
_SUCCESSOR_OUTPUTS: Mapping[str, Mapping[str, str]] = MappingProxyType(
    {
        "output1": MappingProxyType(
            {
                "file_type": "bin",
                "variable": "mhd_w_bcc",
                "id": "mhd_w_bcc",
                "dt": "100.0",
                "ghost_zones": "false",
            }
        ),
        "output2": MappingProxyType(
            {
                "file_type": "bin",
                "variable": "prtcl_rho",
                "id": "prtcl_rho",
                "dt": "100.0",
                "ghost_zones": "false",
            }
        ),
        "output3": MappingProxyType(
            {
                "file_type": "bin",
                "variable": "prtcl_jx",
                "id": "prtcl_jx",
                "dt": "100.0",
                "ghost_zones": "false",
            }
        ),
        "output4": MappingProxyType(
            {
                "file_type": "bin",
                "variable": "prtcl_jy",
                "id": "prtcl_jy",
                "dt": "100.0",
                "ghost_zones": "false",
            }
        ),
        "output5": MappingProxyType(
            {
                "file_type": "bin",
                "variable": "prtcl_jz",
                "id": "prtcl_jz",
                "dt": "100.0",
                "ghost_zones": "false",
            }
        ),
        "output6": MappingProxyType(
            {
                "file_type": "bin",
                "variable": "mhd_j2",
                "id": "mhd_j2",
                "dt": "100.0",
                "ghost_zones": "false",
            }
        ),
        "output7": MappingProxyType(
            {
                "file_type": "pvtk",
                "variable": "prtcl_all",
                "id": "prtcl_all",
                "dt": "100.0",
            }
        ),
        "output8": MappingProxyType({"file_type": "hst", "dt": "10.0"}),
        "output9": MappingProxyType({"file_type": "rst", "dt": "100.0"}),
    }
)


class ProductionScienceError(ValueError):
    """Raised when a production-science diagnostic cannot be trusted."""


_UNDERLYING_CONTRACT_EXCEPTIONS = (
    output_primitives.AnalysisError,
    particle_primitives.ParticleReducerError,
    frozen_model.ModelContractError,
    OSError,
    UnicodeError,
    KeyError,
    IndexError,
    AttributeError,
    TypeError,
    ValueError,
    OverflowError,
    FloatingPointError,
    RecursionError,
)


def _public_contract(label: str):
    """Normalize expected parser, analysis, and malformed-input failures."""

    def decorate(function):
        @wraps(function)
        def wrapped(*args: object, **kwargs: object):
            try:
                return function(*args, **kwargs)
            except ProductionScienceError:
                raise
            except _UNDERLYING_CONTRACT_EXCEPTIONS as error:
                raise ProductionScienceError(f"{label} failed: {error}") from error

        return wrapped

    return decorate


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ProductionScienceError(message)


def _finite_scalar(value: object, label: str, *, minimum: float | None = None) -> float:
    _require(
        isinstance(value, Real) and not isinstance(value, (bool, np.bool_)),
        f"{label} must be a real scalar",
    )
    result = float(value)
    _require(math.isfinite(result), f"{label} must be finite")
    if minimum is not None:
        _require(result >= minimum, f"{label} must be at least {minimum}")
    return result


def _finite_array(values: object, label: str, *, ndim: int | None = None) -> np.ndarray:
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ProductionScienceError(f"{label} must be numeric") from error
    if ndim is not None:
        _require(array.ndim == ndim, f"{label} must be {ndim}-dimensional")
    _require(array.size > 0, f"{label} must not be empty")
    _require(np.all(np.isfinite(array)), f"{label} must be finite")
    return array


def _integer_vector(values: object, label: str) -> np.ndarray:
    array = np.asarray(values)
    _require(array.ndim == 1, f"{label} must be one-dimensional")
    _require(array.dtype.kind in "iu", f"{label} must contain decoded integers")
    return array


def _require_same_length(reference: np.ndarray, **arrays: np.ndarray) -> None:
    for label, array in arrays.items():
        _require(
            array.shape[0] == reference.shape[0],
            f"{label} particle count disagrees",
        )


def _normalized_runtime_parameters(
    parameters: object, *, label: str
) -> dict[str, dict[str, str]]:
    """Remove only AthenaK's sequential output-state counters before comparison."""
    _require(isinstance(parameters, Mapping), f"{label} input parameters must be a mapping")
    normalized: dict[str, dict[str, str]] = {}
    for block_name, block_values in parameters.items():
        _require(
            isinstance(block_name, str) and isinstance(block_values, Mapping),
            f"{label} input parameter block is malformed",
        )
        normalized_block: dict[str, str] = {}
        for name, value in block_values.items():
            _require(
                isinstance(name, str) and isinstance(value, str),
                f"{label} input parameter value is malformed",
            )
            if block_name in _SUCCESSOR_OUTPUTS and name in _MUTABLE_OUTPUT_PARAMETER_NAMES:
                if name == "file_number":
                    try:
                        counter = int(value)
                    except ValueError as error:
                        raise ProductionScienceError(
                            f"{label} output file_number must be an integer"
                        ) from error
                    _require(
                        counter >= 0 and value == str(counter),
                        f"{label} output file_number must be a canonical non-negative integer",
                    )
                else:
                    try:
                        last_time = float(value)
                    except ValueError as error:
                        raise ProductionScienceError(
                            f"{label} output last_time must be numeric"
                        ) from error
                    _require(
                        math.isfinite(last_time) and last_time >= 0.0,
                        f"{label} output last_time must be finite and non-negative",
                    )
                continue
            normalized_block[name] = value
        normalized[block_name] = normalized_block
    return normalized


def _mesh_time_projection(observed_committed_time: float) -> float:
    return float(format(observed_committed_time, ".6g"))


def _translated_window(center: float, offsets: tuple[float, float]) -> tuple[float, float]:
    window = (center + offsets[0], center + offsets[1])
    _require(all(math.isfinite(value) for value in window), "translated window is invalid")
    return window


def _parse_athinput(path: Path) -> dict[str, dict[str, str]]:
    blocks: dict[str, dict[str, str]] = {}
    current: str | None = None
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), 1
    ):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            current = line[1:-1].strip()
            _require(
                bool(current) and current not in blocks,
                f"{path}:{line_number}: invalid or duplicate block",
            )
            blocks[current] = {}
            continue
        _require(
            current is not None and "=" in line,
            f"{path}:{line_number}: malformed Athena input",
        )
        name, value = (part.strip() for part in line.split("=", 1))
        _require(
            bool(name) and bool(value) and name not in blocks[current],
            f"{path}:{line_number}: invalid or duplicate parameter",
        )
        blocks[current][name] = value
    return blocks


@_public_contract("successor deck validation")
def validate_successor_deck(path: Path = SUCCESSOR_DECK) -> dict[str, Any]:
    """Validate the additive deck while explicitly refusing launch authority."""
    baseline = _parse_athinput(BASELINE_DECK)
    successor = _parse_athinput(path)
    ignored = {"comment", "job", *(f"output{index}" for index in range(1, 10))}
    baseline_physics = {
        block: values for block, values in baseline.items() if block not in ignored
    }
    successor_physics = {
        block: values
        for block, values in successor.items()
        if block not in ignored
        and block != "q011_section54_production_science_successor_v1"
    }
    _require(
        successor_physics == baseline_physics,
        "successor deck changed the historical Q011 physical baseline",
    )
    _require(
        successor.get("job", {}).get("basename")
        == "pic_parallel_shock_section54_production_science_successor_v1_vl2_tsc",
        "successor deck basename drifted",
    )
    _require(
        successor.get("q011_section54_production_science_successor_v1")
        == dict(_SUCCESSOR_BLOCK),
        "successor deck diagnostic contract block drifted",
    )
    for block, expected in _SUCCESSOR_OUTPUTS.items():
        _require(successor.get(block) == dict(expected), f"successor deck {block} drifted")
    frozen_model.parse_deck_contract(path.read_text(encoding="utf-8"))
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_production_science_successor_deck_validation",
        "successor_id": SUCCESSOR_ID,
        "path": str(path.relative_to(REPO_ROOT)),
        "historical_physical_baseline_preserved": True,
        "output_ids": [
            values["id"] for values in _SUCCESSOR_OUTPUTS.values() if "id" in values
        ],
        "history_cadence_omega0_inverse": 10.0,
        "launch_authorized": False,
        "policy_mutation_authorized": False,
        "claim_closure_authorized": False,
    }


@dataclass(frozen=True)
class ComposedMHDState:
    """One complete 2D finest-level composite of an ``mhd_w_bcc`` snapshot."""

    nominal_slot_time: float
    observed_committed_time: float
    x1_faces: np.ndarray
    x2_faces: np.ndarray
    x3_faces: np.ndarray
    source_levels_y_x: np.ndarray
    target_level: int
    fields_y_x: Mapping[str, np.ndarray]


@_public_contract("full MHD state composition")
def compose_full_mhd_state(
    dataset: output_primitives.AthenaBinaryDataset,
    *,
    nominal_slot_time: object | None = None,
    observed_committed_time: object | None = None,
    target_level: int | None = None,
) -> ComposedMHDState:
    """Compose and validate the exact unordered ``mhd_w_bcc`` field inventory."""
    _require(
        isinstance(dataset, output_primitives.AthenaBinaryDataset),
        "full MHD state requires an AthenaBinaryDataset",
    )
    _require(
        len(dataset.variable_names) == len(MHD_PRIMITIVE_FIELDS)
        and set(dataset.variable_names) == set(MHD_PRIMITIVE_FIELDS),
        "mhd_w_bcc variable inventory drifted",
    )
    _require(
        dataset.root_grid_shape[2] == 1 and dataset.meshblock_shape[2] == 1,
        "production-science Q011 state requires a singleton x3 axis",
    )
    if nominal_slot_time is None and observed_committed_time is None:
        nominal = _finite_scalar(dataset.time, "nominal slot time", minimum=0.0)
        observed = nominal
    else:
        _require(
            nominal_slot_time is not None and observed_committed_time is not None,
            "nominal and observed times must be supplied together",
        )
        nominal = _finite_scalar(nominal_slot_time, "nominal slot time", minimum=0.0)
        observed = _finite_scalar(
            observed_committed_time, "observed committed time", minimum=0.0
        )
    _require(
        dataset.time == _mesh_time_projection(observed),
        "mhd_w_bcc embedded time differs from observed committed time projection",
    )

    composites = {
        field: output_primitives.compose_leaf_field(
            dataset, field, target_level=target_level
        )
        for field in MHD_PRIMITIVE_FIELDS
    }
    reference = composites[MHD_PRIMITIVE_FIELDS[0]]
    _require(
        reference.values.shape[0] == 1 and reference.x3_faces.size == 2,
        "composed Q011 state did not preserve singleton x3",
    )
    fields: dict[str, np.ndarray] = {}
    for field in MHD_PRIMITIVE_FIELDS:
        composite = composites[field]
        _require(
            composite.target_level == reference.target_level
            and np.array_equal(composite.x1_faces, reference.x1_faces)
            and np.array_equal(composite.x2_faces, reference.x2_faces)
            and np.array_equal(composite.x3_faces, reference.x3_faces)
            and np.array_equal(composite.source_levels, reference.source_levels),
            "mhd_w_bcc composite grids disagree",
        )
        values = np.array(composite.values[0], dtype=np.float64, copy=True)
        _require(np.all(np.isfinite(values)), f"mhd_w_bcc {field} must be finite")
        values.setflags(write=False)
        fields[field] = values
    _require(np.all(fields["dens"] > 0.0), "mhd_w_bcc density must remain positive")
    _require(
        np.all(fields["eint"] >= 0.0),
        "mhd_w_bcc internal energy must remain non-negative",
    )
    levels = np.array(reference.source_levels[0], dtype=np.int64, copy=True)
    levels.setflags(write=False)
    return ComposedMHDState(
        nominal_slot_time=nominal,
        observed_committed_time=observed,
        x1_faces=np.array(reference.x1_faces, copy=True),
        x2_faces=np.array(reference.x2_faces, copy=True),
        x3_faces=np.array(reference.x3_faces, copy=True),
        source_levels_y_x=levels,
        target_level=reference.target_level,
        fields_y_x=MappingProxyType(fields),
    )


def _compose_matched_current_fields(
    mhd_dataset: output_primitives.AthenaBinaryDataset,
    state: ComposedMHDState,
    current_datasets: Mapping[str, output_primitives.AthenaBinaryDataset],
) -> Mapping[str, np.ndarray]:
    _require(
        isinstance(current_datasets, Mapping),
        "current products must be supplied as a mapping",
    )
    _require(
        set(current_datasets) == set(CURRENT_PRODUCT_FIELDS),
        "current product inventory drifted",
    )
    reference_parameters = _normalized_runtime_parameters(
        mhd_dataset.input_parameters, label="mhd_w_bcc"
    )
    fields: dict[str, np.ndarray] = {}
    for product, field in CURRENT_PRODUCT_FIELDS.items():
        dataset = current_datasets[product]
        product_parameters = _normalized_runtime_parameters(
            dataset.input_parameters, label=product
        )
        _require(
            isinstance(dataset, output_primitives.AthenaBinaryDataset),
            f"{product} requires an AthenaBinaryDataset",
        )
        _require(
            dataset.variable_names == (field,),
            f"{product} variable inventory drifted",
        )
        _require(
            dataset.time == _mesh_time_projection(state.observed_committed_time)
            and dataset.cycle == mhd_dataset.cycle
            and dataset.location_size == mhd_dataset.location_size
            and dataset.variable_size == mhd_dataset.variable_size
            and dataset.root_grid_shape == mhd_dataset.root_grid_shape
            and dataset.meshblock_shape == mhd_dataset.meshblock_shape
            and dataset.domain_bounds == mhd_dataset.domain_bounds
            and product_parameters == reference_parameters,
            f"{product} metadata disagrees with matched mhd_w_bcc",
        )
        composite = output_primitives.compose_leaf_field(
            dataset, field, target_level=state.target_level
        )
        _require(
            composite.target_level == state.target_level
            and np.array_equal(composite.x1_faces, state.x1_faces)
            and np.array_equal(composite.x2_faces, state.x2_faces)
            and np.array_equal(composite.x3_faces, state.x3_faces)
            and np.array_equal(composite.source_levels[0], state.source_levels_y_x),
            f"{product} composite grid disagrees with matched mhd_w_bcc",
        )
        values = np.array(composite.values[0], dtype=np.float64, copy=True)
        _require(np.all(np.isfinite(values)), f"{product} must be finite")
        values.setflags(write=False)
        fields[product] = values
    _require(np.all(fields["prtcl_rho"] >= 0.0), "prtcl_rho must be non-negative")
    _require(np.all(fields["mhd_j2"] >= 0.0), "mhd_j2 must be non-negative")
    return MappingProxyType(fields)


def _x1_centers(state: ComposedMHDState) -> np.ndarray:
    return 0.5 * (state.x1_faces[:-1] + state.x1_faces[1:])


def _cell_areas(state: ComposedMHDState) -> np.ndarray:
    return np.diff(state.x2_faces)[:, None] * np.diff(state.x1_faces)[None, :]


def _cell_volumes(state: ComposedMHDState) -> np.ndarray:
    _require(state.x3_faces.size == 2, "state x3 thickness is ambiguous")
    return _cell_areas(state) * float(state.x3_faces[1] - state.x3_faces[0])


def _window_mask(
    state: ComposedMHDState,
    window: tuple[float, float],
    *,
    label: str,
    require_contained: bool = True,
    closed: bool = False,
) -> np.ndarray:
    lower, upper = window
    _require(
        math.isfinite(lower) and math.isfinite(upper) and upper > lower,
        f"{label} must be finite and increasing",
    )
    if require_contained:
        _require(
            lower >= state.x1_faces[0] and upper <= state.x1_faces[-1],
            f"{label} escaped the retained domain",
        )
    if closed:
        selected_x = (_x1_centers(state) >= lower) & (_x1_centers(state) <= upper)
    else:
        selected_x = (_x1_centers(state) > lower) & (_x1_centers(state) < upper)
    mask = np.broadcast_to(selected_x[None, :], state.fields_y_x["dens"].shape)
    _require(np.any(mask), f"{label} selects no cells")
    return mask


def _weighted_mean(values: np.ndarray, weights: np.ndarray, label: str) -> float:
    selected = _finite_array(values, f"{label} values").reshape(-1)
    selected_weights = _finite_array(weights, f"{label} weights").reshape(-1)
    _require(selected.size == selected_weights.size, f"{label} sizes disagree")
    _require(np.all(selected_weights >= 0.0), f"{label} weights must be non-negative")
    total = float(np.sum(selected_weights))
    _require(total > 0.0 and math.isfinite(total), f"{label} total weight must be positive")
    return float(np.sum(selected * selected_weights) / total)


def _weighted_quantile(
    values: np.ndarray, weights: np.ndarray, quantile: float, *, label: str
) -> float:
    q = _finite_scalar(quantile, f"{label} quantile")
    _require(0.0 <= q <= 1.0, f"{label} quantile must lie in [0, 1]")
    samples = _finite_array(values, f"{label} values").reshape(-1)
    sample_weights = _finite_array(weights, f"{label} weights").reshape(-1)
    _require(samples.size == sample_weights.size, f"{label} sizes disagree")
    positive = sample_weights > 0.0
    _require(np.any(positive), f"{label} requires positive total weight")
    order = np.argsort(samples[positive], kind="stable")
    ordered = samples[positive][order]
    ordered_weights = sample_weights[positive][order]
    cumulative = np.cumsum(ordered_weights)
    threshold = q * float(cumulative[-1])
    index = int(np.searchsorted(cumulative, threshold, side="left"))
    return float(ordered[min(index, ordered.size - 1)])


def _area_weighted_profile(values_y_x: np.ndarray, areas_y_x: np.ndarray) -> list[float]:
    column_areas = np.sum(areas_y_x, axis=0)
    _require(np.all(column_areas > 0.0), "profile column areas must be positive")
    return (np.sum(values_y_x * areas_y_x, axis=0) / column_areas).tolist()


def _scalar_statistics(
    values_y_x: np.ndarray,
    areas_y_x: np.ndarray,
    mask_y_x: np.ndarray,
    *,
    label: str,
) -> dict[str, Any]:
    selected = values_y_x[mask_y_x]
    weights = areas_y_x[mask_y_x]
    return {
        "selected_cell_count": int(np.count_nonzero(mask_y_x)),
        "selected_area": float(np.sum(weights)),
        "area_weighted_mean": _weighted_mean(selected, weights, label),
        "area_weighted_rms": math.sqrt(
            _weighted_mean(selected * selected, weights, f"{label} squared")
        ),
        "local_minimum": float(np.min(selected)),
        "local_maximum": float(np.max(selected)),
        "area_weighted_quantiles": {
            f"q{int(round(100.0 * quantile)):02d}": _weighted_quantile(
                selected, weights, quantile, label=label
            )
            for quantile in MAGNETIC_QUANTILES
        },
    }


def _vector_current_statistics(
    vector_y_x: np.ndarray,
    areas_y_x: np.ndarray,
    mask_y_x: np.ndarray,
    *,
    label: str,
) -> dict[str, Any]:
    _require(vector_y_x.shape[0] == 3, f"{label} must have three components")
    magnitude = np.sqrt(np.sum(vector_y_x * vector_y_x, axis=0))
    weights = areas_y_x[mask_y_x]
    return {
        "selected_cell_count": int(np.count_nonzero(mask_y_x)),
        "selected_area": float(np.sum(weights)),
        "area_weighted_mean_vector": [
            _weighted_mean(
                vector_y_x[index][mask_y_x],
                weights,
                f"{label} component {component}",
            )
            for index, component in enumerate(("x", "y", "z"))
        ],
        "area_weighted_mean_magnitude": _weighted_mean(
            magnitude[mask_y_x], weights, f"{label} magnitude"
        ),
        "area_weighted_rms_magnitude": math.sqrt(
            _weighted_mean(
                magnitude[mask_y_x] ** 2,
                weights,
                f"{label} magnitude squared",
            )
        ),
        "local_maximum_magnitude": float(np.max(magnitude[mask_y_x])),
        "area_weighted_quantiles_magnitude": {
            f"q{int(round(100.0 * quantile)):02d}": _weighted_quantile(
                magnitude[mask_y_x],
                weights,
                quantile,
                label=f"{label} magnitude",
            )
            for quantile in MAGNETIC_QUANTILES
        },
    }


def _magnetic_amplification_statistics(
    state: ComposedMHDState,
    bmag_y_x: np.ndarray,
    window: tuple[float, float],
    *,
    label: str,
    qualification_role: str,
    closed_window: bool = False,
) -> dict[str, Any]:
    mask = _window_mask(state, window, label=label, closed=closed_window)
    areas = _cell_areas(state)[mask]
    amplification = bmag_y_x[mask] / REFERENCE_B0
    selected_area = float(np.sum(areas))
    return {
        "qualification_role": qualification_role,
        "window_c_over_omega_pi": list(window),
        "reference_b0": REFERENCE_B0,
        "selected_cell_count": int(np.count_nonzero(mask)),
        "selected_area": selected_area,
        "area_weighted_mean_abs_b_over_b0": _weighted_mean(
            amplification, areas, label
        ),
        "area_weighted_rms_abs_b_over_b0": math.sqrt(
            _weighted_mean(amplification * amplification, areas, f"{label} squared")
        ),
        "local_max_abs_b_over_b0": float(np.max(amplification)),
        "area_weighted_quantiles_abs_b_over_b0": {
            f"q{int(round(100.0 * quantile)):02d}": _weighted_quantile(
                amplification, areas, quantile, label=label
            )
            for quantile in MAGNETIC_QUANTILES
        },
        "area_fractions": {
            f"area_fraction_abs_b_over_b0_ge_{threshold:g}": float(
                np.sum(areas[amplification >= threshold]) / selected_area
            )
            for threshold in MAGNETIC_AREA_FRACTION_THRESHOLDS_OVER_B0
        },
    }


def _mhd_energy_components(
    state: ComposedMHDState, mask: np.ndarray
) -> dict[str, float]:
    volume = _cell_volumes(state)
    rho = state.fields_y_x["dens"]
    velocity_squared = sum(
        state.fields_y_x[field] ** 2 for field in ("velx", "vely", "velz")
    )
    magnetic_squared = sum(
        state.fields_y_x[field] ** 2 for field in ("bcc1", "bcc2", "bcc3")
    )
    gas_kinetic = float(np.sum(0.5 * rho[mask] * velocity_squared[mask] * volume[mask]))
    gas_internal = float(np.sum(state.fields_y_x["eint"][mask] * volume[mask]))
    magnetic = float(np.sum(0.5 * magnetic_squared[mask] * volume[mask]))
    _require(
        all(math.isfinite(value) and value >= 0.0 for value in (gas_kinetic, gas_internal, magnetic)),
        "MHD energy components must be finite and non-negative",
    )
    return {
        "gas_kinetic_energy": gas_kinetic,
        "gas_internal_energy": gas_internal,
        "gas_energy": gas_kinetic + gas_internal,
        "magnetic_energy": magnetic,
    }


@_public_contract("MHD snapshot reduction")
def reduce_mhd_snapshot(
    dataset: output_primitives.AthenaBinaryDataset,
    *,
    nominal_slot_time: object | None = None,
    observed_committed_time: object | None = None,
    target_level: int | None = None,
) -> dict[str, Any]:
    """Reduce one full-state snapshot into shock, field, and MHD-energy records."""
    state = compose_full_mhd_state(
        dataset,
        nominal_slot_time=nominal_slot_time,
        observed_committed_time=observed_committed_time,
        target_level=target_level,
    )
    x1 = _x1_centers(state)
    areas = _cell_areas(state)
    density_profile = np.asarray(
        _area_weighted_profile(state.fields_y_x["dens"], areas), dtype=np.float64
    )
    x_ideal = frozen_model.x_ideal(state.observed_committed_time)
    front = output_primitives.detect_shock_front(
        x1,
        density_profile,
        search_window=_translated_window(
            x_ideal, SHOCK_SEARCH_OFFSETS_C_OVER_OMEGA_PI
        ),
        gradient_sign="negative",
    )
    offset = front.x1 - x_ideal
    _require(
        abs(offset) <= MAX_FRONT_OFFSET_C_OVER_OMEGA_PI,
        "detected shock front exceeds its fixed ideal-surface offset bound",
    )

    downstream_window = _translated_window(
        front.x1, DOWNSTREAM_STATE_OFFSETS_C_OVER_OMEGA_PI
    )
    upstream_window = _translated_window(front.x1, UPSTREAM_STATE_OFFSETS_C_OVER_OMEGA_PI)
    shock_energy_window = _translated_window(front.x1, SHOCK_ENERGY_OFFSETS_C_OVER_OMEGA_PI)
    downstream_mask = _window_mask(state, downstream_window, label="downstream state window")
    upstream_mask = _window_mask(state, upstream_window, label="upstream state window")
    shock_energy_mask = _window_mask(
        state, shock_energy_window, label="shock energy window"
    )
    all_cells = np.ones_like(state.fields_y_x["dens"], dtype=bool)
    downstream_density = _weighted_mean(
        state.fields_y_x["dens"][downstream_mask],
        areas[downstream_mask],
        "downstream density",
    )
    upstream_density = _weighted_mean(
        state.fields_y_x["dens"][upstream_mask],
        areas[upstream_mask],
        "upstream density",
    )
    _require(upstream_density > 0.0, "upstream density must be positive")

    bmag = np.sqrt(
        sum(state.fields_y_x[field] ** 2 for field in ("bcc1", "bcc2", "bcc3"))
    )
    ideal_upstream_window = _translated_window(
        x_ideal, UPSTREAM_B_OFFSETS_C_OVER_OMEGA_PI
    )
    qualifying_amplification = _magnetic_amplification_statistics(
        state,
        bmag,
        ideal_upstream_window,
        label="ideal-surface-relative upstream magnetic amplification",
        qualification_role=(
            "preregistered_qualifying_metric_at_nominal_t500_"
            "detected_front_substitution_forbidden"
        ),
        closed_window=True,
    )
    qualifying_mean = qualifying_amplification["area_weighted_mean_abs_b_over_b0"]
    qualifying_amplification["acceptance_gate"] = {
        "nominal_snapshot_omega0_inverse": 500.0,
        "acceptance_range": list(UPSTREAM_B0_QUALIFYING_GATE),
        "evaluated": state.nominal_slot_time == 500.0,
        "passed": (
            UPSTREAM_B0_QUALIFYING_GATE[0]
            <= qualifying_mean
            <= UPSTREAM_B0_QUALIFYING_GATE[1]
            if state.nominal_slot_time == 500.0
            else None
        ),
    }
    supplemental_amplification = _magnetic_amplification_statistics(
        state,
        bmag,
        upstream_window,
        label="detected-front-relative upstream magnetic amplification",
        qualification_role=(
            "supplemental_only_not_a_substitute_for_the_preregistered_"
            "ideal_surface_relative_metric"
        ),
    )

    levels, level_counts = np.unique(state.source_levels_y_x, return_counts=True)
    profiles = {
        field: _area_weighted_profile(state.fields_y_x[field], areas)
        for field in MHD_PRIMITIVE_FIELDS
    }
    profiles["pressure"] = _area_weighted_profile(
        (GAS_GAMMA - 1.0) * state.fields_y_x["eint"], areas
    )
    profiles["bmag"] = _area_weighted_profile(bmag, areas)
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_production_science_mhd_snapshot_v1",
        "successor_id": SUCCESSOR_ID,
        "nominal_slot_time": state.nominal_slot_time,
        "observed_committed_time": state.observed_committed_time,
        "ideal_surface_x1_c_over_omega_pi": x_ideal,
        "detected_front": {
            "detector": "unique_strongest_negative_density_gradient",
            "physical_orientation": (
                "downstream_density_is_higher_at_lower_x1_so_the_outward_front_"
                "is_a_negative_density_gradient"
            ),
            "preregistration_status": (
                "corrected_detector_requires_preregistration_supersession_"
                "before_qualifying_use"
            ),
            "qualifying_use_authorized": False,
            "x_front_c_over_omega_pi": front.x1,
            "density_gradient": front.density_gradient,
            "offset_from_ideal_surface_c_over_omega_pi": offset,
            "maximum_absolute_offset_c_over_omega_pi": MAX_FRONT_OFFSET_C_OVER_OMEGA_PI,
        },
        "shock_state": {
            "downstream_definition": "x1 < detected_front",
            "upstream_definition": "x1 > detected_front",
            "downstream_window_c_over_omega_pi": list(downstream_window),
            "upstream_window_c_over_omega_pi": list(upstream_window),
            "downstream_area_weighted_density": downstream_density,
            "upstream_area_weighted_density": upstream_density,
            "compression_ratio_downstream_over_upstream": (
                downstream_density / upstream_density
            ),
        },
        "upstream_magnetic_amplification": {
            "qualifying_ideal_surface_relative": qualifying_amplification,
            "supplemental_detected_front_relative": supplemental_amplification,
        },
        "mhd_energy_components": {
            "interpretation": (
                "instantaneous_wall_frame_cell_centered_diagnostic_not_"
                "conserved_energy_closure"
            ),
            "excluded_closure_terms": [
                "boundary_fluxes",
                "particle_injection",
                "particle_removal",
                "gas_subtraction",
                "face_centered_magnetic_energy",
                "conserved_total_energy_history",
            ],
            "closure_successor_required": True,
            "full_domain": _mhd_energy_components(state, all_cells),
            "shock_centered_window": {
                "window_c_over_omega_pi": list(shock_energy_window),
                **_mhd_energy_components(state, shock_energy_mask),
            },
        },
        "full_state": {
            "source_product": "mhd_w_bcc",
            "raw_primitive_field_inventory": list(MHD_PRIMITIVE_FIELDS),
            "derived_profile_field_inventory": list(DERIVED_PROFILE_FIELDS),
            "x1_centers_c_over_omega_pi": x1.tolist(),
            "y_area_weighted_profiles": profiles,
            "target_composite_level": state.target_level,
            "source_level_cell_counts": [
                {"source_level": int(level), "cell_count": int(count)}
                for level, count in zip(levels, level_counts)
            ],
        },
    }


@_public_contract("CR current snapshot reduction")
def reduce_cr_current_snapshot(
    mhd_dataset: output_primitives.AthenaBinaryDataset,
    current_datasets: Mapping[str, output_primitives.AthenaBinaryDataset],
    *,
    nominal_slot_time: object,
    observed_committed_time: object,
    detected_front_x1_c_over_omega_pi: object,
    target_level: int | None = None,
) -> dict[str, Any]:
    """Reduce matched deposited moments into lab- and gas-frame CR currents."""
    state = compose_full_mhd_state(
        mhd_dataset,
        nominal_slot_time=nominal_slot_time,
        observed_committed_time=observed_committed_time,
        target_level=target_level,
    )
    fields = _compose_matched_current_fields(mhd_dataset, state, current_datasets)
    front_x = _finite_scalar(
        detected_front_x1_c_over_omega_pi, "detected front x1"
    )
    x_ideal = frozen_model.x_ideal(state.observed_committed_time)
    _require(
        abs(front_x - x_ideal) <= MAX_FRONT_OFFSET_C_OVER_OMEGA_PI,
        "current diagnostic detected front exceeds its ideal-surface offset bound",
    )
    areas = _cell_areas(state)
    rho_cr = fields["prtcl_rho"]
    lab_current = np.stack(
        [fields[f"prtcl_j{component}"] for component in ("x", "y", "z")],
        axis=0,
    )
    gas_velocity = np.stack(
        [state.fields_y_x[f"vel{component}"] for component in ("x", "y", "z")],
        axis=0,
    )
    gas_current = lab_current - rho_cr[None, :, :] * gas_velocity
    _require(np.all(np.isfinite(gas_current)), "gas-frame CR current must be finite")
    transform_residual = gas_current + rho_cr[None, :, :] * gas_velocity - lab_current
    _require(
        np.all(np.isfinite(transform_residual))
        and float(np.max(np.abs(transform_residual))) <= 1.0e-12,
        "gas-frame CR current transform failed closure",
    )

    ideal_upstream_window = _translated_window(
        x_ideal, UPSTREAM_B_OFFSETS_C_OVER_OMEGA_PI
    )
    detected_upstream_window = _translated_window(
        front_x, UPSTREAM_STATE_OFFSETS_C_OVER_OMEGA_PI
    )
    masks = {
        "full_domain": np.ones_like(rho_cr, dtype=bool),
        "qualifying_ideal_surface_relative_upstream": _window_mask(
            state,
            ideal_upstream_window,
            label="ideal-surface-relative current window",
            closed=True,
        ),
        "supplemental_detected_front_relative_upstream": _window_mask(
            state,
            detected_upstream_window,
            label="detected-front-relative current window",
        ),
    }
    statistics = {
        region: {
            "window_c_over_omega_pi": (
                None
                if region == "full_domain"
                else list(
                    ideal_upstream_window
                    if region == "qualifying_ideal_surface_relative_upstream"
                    else detected_upstream_window
                )
            ),
            "prtcl_rho": _scalar_statistics(
                rho_cr, areas, mask, label=f"{region} prtcl_rho"
            ),
            "j_cr_lab": _vector_current_statistics(
                lab_current, areas, mask, label=f"{region} lab-frame CR current"
            ),
            "j_cr_gas": _vector_current_statistics(
                gas_current, areas, mask, label=f"{region} gas-frame CR current"
            ),
            "mhd_current_squared": _scalar_statistics(
                fields["mhd_j2"],
                areas,
                mask,
                label=f"{region} MHD current squared",
            ),
        }
        for region, mask in masks.items()
    }
    profiles = {
        "prtcl_rho": _area_weighted_profile(rho_cr, areas),
        "j_cr_lab_x": _area_weighted_profile(lab_current[0], areas),
        "j_cr_lab_y": _area_weighted_profile(lab_current[1], areas),
        "j_cr_lab_z": _area_weighted_profile(lab_current[2], areas),
        "j_cr_lab_magnitude": _area_weighted_profile(
            np.sqrt(np.sum(lab_current * lab_current, axis=0)), areas
        ),
        "j_cr_gas_x": _area_weighted_profile(gas_current[0], areas),
        "j_cr_gas_y": _area_weighted_profile(gas_current[1], areas),
        "j_cr_gas_z": _area_weighted_profile(gas_current[2], areas),
        "j_cr_gas_magnitude": _area_weighted_profile(
            np.sqrt(np.sum(gas_current * gas_current, axis=0)), areas
        ),
        "mhd_current_squared": _area_weighted_profile(fields["mhd_j2"], areas),
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_production_science_cr_current_snapshot_v1",
        "successor_id": SUCCESSOR_ID,
        "nominal_slot_time": state.nominal_slot_time,
        "observed_committed_time": state.observed_committed_time,
        "source_products": {
            **dict(CURRENT_PRODUCT_FIELDS),
            "mhd_w_bcc_gas_velocity_fields": ["velx", "vely", "velz"],
        },
        "matched_metadata_contract": {
            "immutable_parameter_comparison": (
                "exact_after_validating_and_removing_only_sequential_output_state"
            ),
            "normalized_output_blocks": list(_SUCCESSOR_OUTPUTS),
            "normalized_output_state_fields": sorted(_MUTABLE_OUTPUT_PARAMETER_NAMES),
        },
        "source_representation": {
            "prtcl_rho": "deposited_rho_CR_over_c",
            "prtcl_j_components": "deposited_J_CR_over_c",
            "reported_j_cr_shorthand": (
                "code_normalized_deposited_current_representation_only"
            ),
            "shared_c_inverse_factor": (
                "prtcl_rho_and_prtcl_j_each_carry_one_factor_of_c_inverse"
            ),
            "physical_j_cr_conversion": (
                "external_normalization_required_not_applied_by_this_successor"
            ),
        },
        "frame_transform": {
            "formula": "(J_CR/c)_gas = (J_CR/c)_lab - (rho_CR/c) * u_gas",
            "gas_velocity_source": "matched_mhd_w_bcc_cell_centered_velocity",
            "maximum_absolute_reconstruction_residual": float(
                np.max(np.abs(transform_residual))
            ),
        },
        "mhd_current_squared_interpretation": {
            "source_product": "mhd_j2",
            "source_field": "j2",
            "meaning": "cell_centered_squared_magnitude_of_MHD_curl_B_current",
            "not_particle_current_squared": True,
        },
        "region_roles": {
            "qualifying_ideal_surface_relative_upstream": (
                "preserves_preregistered_ideal_surface_relative_window"
            ),
            "supplemental_detected_front_relative_upstream": "supplemental_only",
        },
        "statistics_by_region": statistics,
        "y_area_weighted_profiles": {
            "x1_centers_c_over_omega_pi": _x1_centers(state).tolist(),
            **profiles,
        },
    }


def _decoded_particle_arrays(
    *,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
) -> dict[str, np.ndarray]:
    decoded_points = _finite_array(points, "particle points", ndim=2)
    decoded_velocity = _finite_array(velocity, "particle velocity", ndim=2)
    _require(
        decoded_points.shape[1:] == (3,) and decoded_velocity.shape[1:] == (3,),
        "particle points and velocity must have shape (nparticle, 3)",
    )
    sources = _integer_vector(cr_source, "particle cr_source")
    births = _finite_array(birth_time, "particle birth_time", ndim=1)
    weights = _finite_array(macro_weight, "particle macro_weight", ndim=1)
    _require(np.all(weights >= 0.0), "particle macro_weight must be non-negative")
    _require_same_length(
        decoded_points,
        velocity=decoded_velocity,
        cr_source=sources,
        birth_time=births,
        macro_weight=weights,
    )
    projected_velocity = decoded_velocity.astype(np.float32)
    restored_velocity = projected_velocity.astype(np.float64)
    _require(
        np.array_equal(decoded_velocity, restored_velocity),
        "particle velocity must be exactly representable decoded float32 PVTK data",
    )
    previous_velocity = np.nextafter(
        projected_velocity, np.float32(-np.inf)
    ).astype(np.float64)
    next_velocity = np.nextafter(
        projected_velocity, np.float32(np.inf)
    ).astype(np.float64)
    velocity_lower = restored_velocity - 0.5 * (
        restored_velocity - previous_velocity
    )
    velocity_upper = restored_velocity + 0.5 * (
        next_velocity - restored_velocity
    )
    _require(
        np.all(np.isfinite(velocity_lower)) and np.all(np.isfinite(velocity_upper)),
        "float32 PVTK velocity rounding interval must be finite",
    )
    speed_squared = np.sum(decoded_velocity * decoded_velocity, axis=1)
    minimum_component_magnitude = np.where(
        (velocity_lower <= 0.0) & (velocity_upper >= 0.0),
        0.0,
        np.minimum(np.abs(velocity_lower), np.abs(velocity_upper)),
    )
    maximum_component_magnitude = np.maximum(
        np.abs(velocity_lower), np.abs(velocity_upper)
    )
    speed_squared_lower = np.sum(minimum_component_magnitude**2, axis=1)
    speed_squared_upper = np.sum(maximum_component_magnitude**2, axis=1)
    _require(
        np.all(speed_squared_upper < PARTICLE_LIGHT_SPEED**2),
        "particle float32 uncertainty interval must remain below the Q011 artificial light speed",
    )
    gamma_squared = 1.0 / (1.0 - speed_squared / PARTICLE_LIGHT_SPEED**2)
    momentum_squared = gamma_squared * speed_squared
    specific_kinetic_energy = momentum_squared / (
        np.sqrt(1.0 + momentum_squared / PARTICLE_LIGHT_SPEED**2) + 1.0
    )
    chi = momentum_squared / UPSTREAM_SPEED_U0**2
    momentum_squared_lower = speed_squared_lower / (
        1.0 - speed_squared_lower / PARTICLE_LIGHT_SPEED**2
    )
    momentum_squared_upper = speed_squared_upper / (
        1.0 - speed_squared_upper / PARTICLE_LIGHT_SPEED**2
    )
    specific_kinetic_energy_lower = momentum_squared_lower / (
        np.sqrt(1.0 + momentum_squared_lower / PARTICLE_LIGHT_SPEED**2) + 1.0
    )
    specific_kinetic_energy_upper = momentum_squared_upper / (
        np.sqrt(1.0 + momentum_squared_upper / PARTICLE_LIGHT_SPEED**2) + 1.0
    )
    chi_lower = momentum_squared_lower / UPSTREAM_SPEED_U0**2
    chi_upper = momentum_squared_upper / UPSTREAM_SPEED_U0**2
    _require(
        np.all(np.isfinite(specific_kinetic_energy))
        and np.all(specific_kinetic_energy >= 0.0)
        and np.all(np.isfinite(chi))
        and np.all(chi >= 0.0),
        "particle energy reconstruction failed",
    )
    _require(
        np.all(np.isfinite(specific_kinetic_energy_lower))
        and np.all(np.isfinite(specific_kinetic_energy_upper))
        and np.all(specific_kinetic_energy_lower <= specific_kinetic_energy)
        and np.all(specific_kinetic_energy <= specific_kinetic_energy_upper)
        and np.all(np.isfinite(chi_lower))
        and np.all(np.isfinite(chi_upper))
        and np.all(chi_lower <= chi)
        and np.all(chi <= chi_upper),
        "float32 PVTK particle-energy uncertainty propagation failed",
    )
    active = (sources == SHOCK_INJECTED_SOURCE) & (births >= BIRTH_TIME_MINIMUM)
    energetic = active & (weights > 0.0)
    _require(np.any(energetic), "particle snapshot has no positive-weight active Q011 CRs")
    return {
        "points": decoded_points,
        "velocity": decoded_velocity,
        "sources": sources,
        "births": births,
        "weights": weights,
        "specific_kinetic_energy": specific_kinetic_energy,
        "specific_kinetic_energy_lower": specific_kinetic_energy_lower,
        "specific_kinetic_energy_upper": specific_kinetic_energy_upper,
        "chi": chi,
        "chi_lower": chi_lower,
        "chi_upper": chi_upper,
        "maximum_float32_velocity_component_interval_width": float(
            np.max(velocity_upper - velocity_lower)
        ),
        "active": active,
        "energetic": energetic,
    }


def _uncertainty_interval(
    lower: object, nominal: object, upper: object, *, label: str
) -> dict[str, float]:
    low = _finite_scalar(lower, f"{label} lower", minimum=0.0)
    center = _finite_scalar(nominal, f"{label} nominal", minimum=0.0)
    high = _finite_scalar(upper, f"{label} upper", minimum=0.0)
    _require(low <= center <= high, f"{label} nominal value escaped its uncertainty interval")
    return {
        "lower": low,
        "nominal": center,
        "upper": high,
        "absolute_width": high - low,
    }


@_public_contract("particle-energy snapshot reduction")
def reduce_particle_energy_snapshot(
    *,
    snapshot_time: object,
    nominal_slot_time: object | None = None,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
) -> dict[str, Any]:
    """Reduce one particle snapshot into maximum-energy and CR-energy records."""
    time = _finite_scalar(snapshot_time, "particle snapshot time", minimum=0.0)
    nominal = (
        time
        if nominal_slot_time is None
        else _finite_scalar(nominal_slot_time, "particle nominal slot time", minimum=0.0)
    )
    arrays = _decoded_particle_arrays(
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )
    active = arrays["active"]
    energetic = arrays["energetic"]
    weights = arrays["weights"]
    chi = arrays["chi"]
    specific_energy = arrays["specific_kinetic_energy"]
    positive_weight_count = int(np.count_nonzero(energetic))
    _require(
        positive_weight_count >= PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES,
        "particle tail quantiles require at least "
        f"{PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES} positive-weight active samples",
    )
    chi_quantiles = {
        f"q{int(round(1000.0 * quantile)):03d}": _weighted_quantile(
            chi[energetic],
            weights[energetic],
            quantile,
            label="active particle chi",
        )
        for quantile in PARTICLE_ENERGY_QUANTILES
    }
    specific_energy_quantiles = {
        f"q{int(round(1000.0 * quantile)):03d}": _weighted_quantile(
            specific_energy[energetic],
            weights[energetic],
            quantile,
            label="active particle specific kinetic energy",
        )
        for quantile in PARTICLE_ENERGY_QUANTILES
    }
    chi_quantile_intervals = {
        label: _uncertainty_interval(
            _weighted_quantile(
                arrays["chi_lower"][energetic],
                weights[energetic],
                quantile,
                label=f"active particle chi {label} lower",
            ),
            chi_quantiles[label],
            _weighted_quantile(
                arrays["chi_upper"][energetic],
                weights[energetic],
                quantile,
                label=f"active particle chi {label} upper",
            ),
            label=f"active particle chi {label}",
        )
        for label, quantile in zip(("q990", "q999"), PARTICLE_ENERGY_QUANTILES)
    }
    specific_energy_quantile_intervals = {
        label: _uncertainty_interval(
            _weighted_quantile(
                arrays["specific_kinetic_energy_lower"][energetic],
                weights[energetic],
                quantile,
                label=f"active particle specific energy {label} lower",
            ),
            specific_energy_quantiles[label],
            _weighted_quantile(
                arrays["specific_kinetic_energy_upper"][energetic],
                weights[energetic],
                quantile,
                label=f"active particle specific energy {label} upper",
            ),
            label=f"active particle specific energy {label}",
        )
        for label, quantile in zip(("q990", "q999"), PARTICLE_ENERGY_QUANTILES)
    }
    maximum_chi = float(np.max(chi[energetic]))
    maximum_specific_energy = float(np.max(specific_energy[energetic]))
    record: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_production_science_particle_energy_snapshot_v1",
        "successor_id": SUCCESSOR_ID,
        "nominal_slot_time": nominal,
        "snapshot_time_omega0_inverse": time,
        "selection": {
            "source_rule": "cr_source == 1",
            "birth_time_rule": "birth_time >= 45",
            "spatial_rule": "none_for_acceleration_and_full_domain_cr_energy",
            "maximum_energy_population": "selected_particles_with_positive_macro_weight",
        },
        "particle_census": {
            "all_particle_count": int(active.size),
            "active_particle_count": int(np.count_nonzero(active)),
            "active_positive_weight_particle_count": int(np.count_nonzero(energetic)),
            "active_macro_weight": float(np.sum(weights[active])),
            "active_cr_macro_mass": float(PARTICLE_MACRO_MASS * np.sum(weights[active])),
        },
        "energy_reconstruction": {
            "input_vector": "pvtk physical velocity",
            "momentum_formula": "p_over_m = gamma(v) * v",
            "specific_kinetic_energy_formula": (
                "p_over_m_squared / (sqrt(1 + p_over_m_squared / C_squared) + 1)"
            ),
            "chi_formula": "p_over_m_squared / u0_squared",
            "particle_macro_mass": PARTICLE_MACRO_MASS,
            "particle_light_speed": PARTICLE_LIGHT_SPEED,
            "u0_over_u_a0": UPSTREAM_SPEED_U0,
        },
        "maximum_particle_energy": {
            "maximum_chi": maximum_chi,
            "maximum_specific_kinetic_energy": maximum_specific_energy,
            "macro_weighted_chi_quantiles": chi_quantiles,
            "macro_weighted_specific_kinetic_energy_quantiles": (
                specific_energy_quantiles
            ),
            "tail_quantile_sample_requirement": {
                "minimum_positive_weight_particle_count": (
                    PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES
                ),
                "observed_positive_weight_particle_count": positive_weight_count,
                "passed": True,
            },
        },
        "float32_pvtk_velocity_uncertainty": {
            "input_contract": (
                "velocity_components_are_exactly_decoded_float32_values"
            ),
            "component_interval": (
                "closed_round_to_nearest_interval_bounded_by_midpoints_to_"
                "adjacent_float32_values"
            ),
            "propagation": (
                "monotonic_component_speed_squared_bounds_to_momentum_squared_"
                "chi_and_relativistic_specific_kinetic_energy"
            ),
            "maximum_velocity_component_interval_width": arrays[
                "maximum_float32_velocity_component_interval_width"
            ],
            "maximum_particle_energy_intervals": {
                "maximum_chi": _uncertainty_interval(
                    np.max(arrays["chi_lower"][energetic]),
                    maximum_chi,
                    np.max(arrays["chi_upper"][energetic]),
                    label="maximum chi",
                ),
                "maximum_specific_kinetic_energy": _uncertainty_interval(
                    np.max(arrays["specific_kinetic_energy_lower"][energetic]),
                    maximum_specific_energy,
                    np.max(arrays["specific_kinetic_energy_upper"][energetic]),
                    label="maximum specific kinetic energy",
                ),
            },
            "macro_weighted_tail_quantile_intervals": {
                "chi": chi_quantile_intervals,
                "specific_kinetic_energy": specific_energy_quantile_intervals,
            },
            "all_nominal_values_contained": True,
        },
        "active_cr_kinetic_energy": float(
            PARTICLE_MACRO_MASS * np.sum(weights[active] * specific_energy[active])
        ),
    }
    if nominal in PREREGISTERED_SPECTRUM_NOMINAL_TIMES:
        spectrum = particle_primitives.reduce_particle_snapshot(
            snapshot_time=time,
            points=points,
            cr_source=cr_source,
            birth_time=birth_time,
            velocity=velocity,
            macro_weight=macro_weight,
            evaluate_late_slope=nominal == particle_primitives.LATE_SLOPE_SNAPSHOT_TIME,
        )
        downstream_spectrum = spectrum["weighted_spectrum"]
        _require(
            int(downstream_spectrum["particle_count_post_filter"])
            >= DOWNSTREAM_SPECTRUM_MINIMUM_PARTICLE_SAMPLES
            and float(downstream_spectrum["total_post_filter_macro_weight"]) > 0.0,
            "preregistered downstream spectrum has insufficient admitted samples",
        )
        record["preregistered_downstream_spectrum"] = {
            "reuse_contract": "q011_section54_particles.reduce_particle_snapshot",
            "nominal_slot_time": nominal,
            "observed_committed_time": time,
            "downstream_classifier": "ideal_injection_surface_detected_front_substitution_forbidden",
            "minimum_positive_fit_bins": (
                particle_primitives.LATE_SLOPE_MINIMUM_POSITIVE_BINS
            ),
            "minimum_admitted_particle_count": (
                DOWNSTREAM_SPECTRUM_MINIMUM_PARTICLE_SAMPLES
            ),
            "reduction": spectrum,
        }
    return record


def _energy_partition(
    mhd_components: Mapping[str, object], cr_kinetic_energy: object
) -> dict[str, Any]:
    gas_kinetic = _finite_scalar(
        mhd_components["gas_kinetic_energy"], "gas kinetic energy", minimum=0.0
    )
    gas_internal = _finite_scalar(
        mhd_components["gas_internal_energy"], "gas internal energy", minimum=0.0
    )
    magnetic = _finite_scalar(
        mhd_components["magnetic_energy"], "magnetic energy", minimum=0.0
    )
    cr = _finite_scalar(cr_kinetic_energy, "CR kinetic energy", minimum=0.0)
    gas = gas_kinetic + gas_internal
    total = gas + magnetic + cr
    _require(total > 0.0 and math.isfinite(total), "partition total energy must be positive")
    return {
        "interpretation": (
            "instantaneous_wall_frame_partition_not_a_conserved_energy_budget"
        ),
        "closure_status": {
            "conserved_energy_closure_evaluated": False,
            "boundary_fluxes_included": False,
            "particle_injection_included": False,
            "particle_removal_included": False,
            "gas_subtraction_included": False,
            "separate_successor_required": True,
        },
        "gas_kinetic_energy": gas_kinetic,
        "gas_internal_energy": gas_internal,
        "gas_energy": gas,
        "magnetic_energy": magnetic,
        "cr_kinetic_energy": cr,
        "partition_total_energy": total,
        "fractions": {
            "gas": gas / total,
            "magnetic": magnetic / total,
            "cr": cr / total,
        },
    }


@_public_contract("matched production-science snapshot reduction")
def reduce_production_science_snapshot(
    mhd_dataset: output_primitives.AthenaBinaryDataset,
    current_datasets: Mapping[str, output_primitives.AthenaBinaryDataset],
    *,
    nominal_slot_time: object,
    observed_committed_time: object,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
    target_level: int | None = None,
) -> dict[str, Any]:
    """Reduce one matched full-state and particle snapshot under successor v1."""
    nominal = _finite_scalar(
        nominal_slot_time, "matched snapshot nominal slot time", minimum=0.0
    )
    _require(
        nominal >= MATCHED_SNAPSHOT_FIRST_ELIGIBLE_NOMINAL_TIME,
        "matched production-science snapshots begin at nominal t=200; "
        "raw t=0 and t=100 products are retained for initial-state and conservation work",
    )
    mhd = reduce_mhd_snapshot(
        mhd_dataset,
        nominal_slot_time=nominal,
        observed_committed_time=observed_committed_time,
        target_level=target_level,
    )
    particle = reduce_particle_energy_snapshot(
        snapshot_time=mhd["observed_committed_time"],
        nominal_slot_time=mhd["nominal_slot_time"],
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )
    current = reduce_cr_current_snapshot(
        mhd_dataset,
        current_datasets,
        nominal_slot_time=mhd["nominal_slot_time"],
        observed_committed_time=mhd["observed_committed_time"],
        detected_front_x1_c_over_omega_pi=mhd["detected_front"][
            "x_front_c_over_omega_pi"
        ],
        target_level=target_level,
    )
    arrays = _decoded_particle_arrays(
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )
    domain = mhd_dataset.domain_bounds[0:2]
    _require(
        np.all(arrays["points"][arrays["active"], 0] >= domain[0])
        and np.all(arrays["points"][arrays["active"], 0] <= domain[1]),
        "active Q011 particle escaped the retained x1 domain",
    )
    shock_window = tuple(
        float(value)
        for value in mhd["mhd_energy_components"]["shock_centered_window"][
            "window_c_over_omega_pi"
        ]
    )
    in_shock_window = (
        arrays["active"]
        & (arrays["points"][:, 0] > shock_window[0])
        & (arrays["points"][:, 0] < shock_window[1])
    )
    cr_shock_energy = float(
        PARTICLE_MACRO_MASS
        * np.sum(
            arrays["weights"][in_shock_window]
            * arrays["specific_kinetic_energy"][in_shock_window]
        )
    )
    full_partition = _energy_partition(
        mhd["mhd_energy_components"]["full_domain"],
        particle["active_cr_kinetic_energy"],
    )
    shock_partition = _energy_partition(
        mhd["mhd_energy_components"]["shock_centered_window"], cr_shock_energy
    )
    shock_partition["window_c_over_omega_pi"] = list(shock_window)
    shock_partition["active_cr_particle_count"] = int(np.count_nonzero(in_shock_window))
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_production_science_snapshot_v1",
        "successor_id": SUCCESSOR_ID,
        "nominal_slot_time": mhd["nominal_slot_time"],
        "observed_committed_time": mhd["observed_committed_time"],
        "mhd": mhd,
        "cr_current": current,
        "particles": particle,
        "energy_partition": {
            "interpretation": (
                "instantaneous_wall_frame_partitions_only_not_conserved_energy_closure"
            ),
            "closure_successor_required": True,
            "full_domain": full_partition,
            "shock_centered_window": shock_partition,
        },
    }


def _linear_fit(times: np.ndarray, values: np.ndarray, *, label: str) -> dict[str, Any]:
    _require(
        times.size >= HISTORY_LINEAR_FIT_MINIMUM_SNAPSHOTS
        and times.size == values.size,
        f"{label} fit requires at least "
        f"{HISTORY_LINEAR_FIT_MINIMUM_SNAPSHOTS} snapshots",
    )
    centered = times - float(np.mean(times))
    denominator = float(np.sum(centered * centered))
    _require(denominator > 0.0, f"{label} fit times must span a finite interval")
    slope = float(np.sum(centered * (values - np.mean(values))) / denominator)
    intercept = float(np.mean(values) - slope * np.mean(times))
    fitted = slope * times + intercept
    residual_sum = float(np.sum((values - fitted) ** 2))
    total_sum = float(np.sum((values - np.mean(values)) ** 2))
    r_squared = 1.0 if total_sum == 0.0 and residual_sum == 0.0 else (
        1.0 - residual_sum / total_sum if total_sum > 0.0 else 0.0
    )
    _require(
        all(math.isfinite(value) for value in (slope, intercept, r_squared)),
        f"{label} fit is non-finite",
    )
    return {
        "sample_count": int(times.size),
        "minimum_sample_count": HISTORY_LINEAR_FIT_MINIMUM_SNAPSHOTS,
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_squared,
    }


def _ordered_snapshot_records(
    records: Sequence[Mapping[str, object]], *, expected_record_type: str
) -> list[Mapping[str, object]]:
    _require(
        isinstance(records, Sequence) and not isinstance(records, (str, bytes)),
        "snapshot history must be a sequence",
    )
    _require(len(records) >= 2, "snapshot history requires at least two snapshots")
    ordered = sorted(
        records,
        key=lambda record: _finite_scalar(
            record["observed_committed_time"], "history observed committed time"
        ),
    )
    times = [
        _finite_scalar(record["observed_committed_time"], "history observed committed time")
        for record in ordered
    ]
    _require(
        all(record.get("record_type") == expected_record_type for record in ordered),
        "snapshot history record type drifted",
    )
    _require(np.all(np.diff(times) > 0.0), "snapshot history times must be unique and increase")
    return ordered


def _validate_production_snapshot_record(record: Mapping[str, object]) -> None:
    _require(
        record.get("schema_version") == SCHEMA_VERSION
        and record.get("record_type") == "q011_section54_production_science_snapshot_v1"
        and record.get("successor_id") == SUCCESSOR_ID,
        "production-science snapshot identity drifted",
    )
    nominal = _finite_scalar(record["nominal_slot_time"], "snapshot nominal slot time")
    observed = _finite_scalar(
        record["observed_committed_time"], "snapshot observed committed time"
    )
    _require(
        nominal >= MATCHED_SNAPSHOT_FIRST_ELIGIBLE_NOMINAL_TIME,
        "matched production-science history contains a pre-t200 snapshot",
    )
    _require(
        float(record["mhd"]["nominal_slot_time"]) == nominal
        and float(record["mhd"]["observed_committed_time"]) == observed
        and float(record["cr_current"]["nominal_slot_time"]) == nominal
        and float(record["cr_current"]["observed_committed_time"]) == observed
        and float(record["particles"]["nominal_slot_time"]) == nominal
        and float(record["particles"]["snapshot_time_omega0_inverse"]) == observed,
        "production-science snapshot component times disagree",
    )
    _require(
        record["mhd"]["detected_front"]["qualifying_use_authorized"] is False
        and record["mhd"]["detected_front"]["preregistration_status"]
        == (
            "corrected_detector_requires_preregistration_supersession_"
            "before_qualifying_use"
        ),
        "negative-gradient detector supersession status drifted",
    )
    _require(
        set(record["mhd"]["upstream_magnetic_amplification"])
        == {
            "qualifying_ideal_surface_relative",
            "supplemental_detected_front_relative",
        },
        "magnetic-amplification metric roles drifted",
    )
    _require(
        record["cr_current"]["frame_transform"]["formula"]
        == "(J_CR/c)_gas = (J_CR/c)_lab - (rho_CR/c) * u_gas"
        and record["cr_current"]["source_representation"]["prtcl_rho"]
        == "deposited_rho_CR_over_c"
        and record["cr_current"]["source_representation"]["prtcl_j_components"]
        == "deposited_J_CR_over_c"
        and set(record["cr_current"]["statistics_by_region"])
        == {
            "full_domain",
            "qualifying_ideal_surface_relative_upstream",
            "supplemental_detected_front_relative_upstream",
        },
        "CR current diagnostic contract drifted",
    )
    _require(
        record["energy_partition"]["interpretation"]
        == "instantaneous_wall_frame_partitions_only_not_conserved_energy_closure"
        and record["energy_partition"]["closure_successor_required"] is True,
        "instantaneous energy-partition scope drifted",
    )
    uncertainty = record["particles"]["float32_pvtk_velocity_uncertainty"]
    _require(
        uncertainty["all_nominal_values_contained"] is True
        and uncertainty["input_contract"]
        == "velocity_components_are_exactly_decoded_float32_values",
        "float32 PVTK velocity uncertainty contract drifted",
    )
    _finite_scalar(
        uncertainty["maximum_velocity_component_interval_width"],
        "maximum float32 velocity component interval width",
        minimum=0.0,
    )
    maximum_energy = record["particles"]["maximum_particle_energy"]
    for observable in ("maximum_chi", "maximum_specific_kinetic_energy"):
        interval = uncertainty["maximum_particle_energy_intervals"][observable]
        validated = _uncertainty_interval(
            interval["lower"],
            interval["nominal"],
            interval["upper"],
            label=f"history {observable} uncertainty",
        )
        _require(
            validated["nominal"] == float(maximum_energy[observable]),
            f"history {observable} uncertainty nominal drifted",
        )
    for family, nominal_key in (
        ("chi", "macro_weighted_chi_quantiles"),
        (
            "specific_kinetic_energy",
            "macro_weighted_specific_kinetic_energy_quantiles",
        ),
    ):
        for label in ("q990", "q999"):
            interval = uncertainty["macro_weighted_tail_quantile_intervals"][family][
                label
            ]
            validated = _uncertainty_interval(
                interval["lower"],
                interval["nominal"],
                interval["upper"],
                label=f"history {family} {label} uncertainty",
            )
            _require(
                validated["nominal"] == float(maximum_energy[nominal_key][label]),
                f"history {family} {label} uncertainty nominal drifted",
            )


@_public_contract("production-science history reduction")
def reduce_production_science_history(
    records: Sequence[Mapping[str, object]],
) -> dict[str, Any]:
    """Reduce matched successor snapshots into acceleration and shock histories."""
    ordered = _ordered_snapshot_records(
        records, expected_record_type="q011_section54_production_science_snapshot_v1"
    )
    for record in ordered:
        _validate_production_snapshot_record(record)
    times = np.asarray(
        [float(record["observed_committed_time"]) for record in ordered], dtype=np.float64
    )
    max_chi = np.asarray(
        [
            float(record["particles"]["maximum_particle_energy"]["maximum_chi"])
            for record in ordered
        ],
        dtype=np.float64,
    )
    max_specific_energy = np.asarray(
        [
            float(
                record["particles"]["maximum_particle_energy"][
                    "maximum_specific_kinetic_energy"
                ]
            )
            for record in ordered
        ],
        dtype=np.float64,
    )
    tail_quantiles = {
        label: np.asarray(
            [
                float(
                    record["particles"]["maximum_particle_energy"][
                        "macro_weighted_chi_quantiles"
                    ][label]
                )
                for record in ordered
            ],
            dtype=np.float64,
        )
        for label in ("q990", "q999")
    }
    tail_specific_energy_quantiles = {
        label: np.asarray(
            [
                float(
                    record["particles"]["maximum_particle_energy"][
                        "macro_weighted_specific_kinetic_energy_quantiles"
                    ][label]
                )
                for record in ordered
            ],
            dtype=np.float64,
        )
        for label in ("q990", "q999")
    }
    tail_sample_counts = np.asarray(
        [
            int(
                record["particles"]["maximum_particle_energy"][
                    "tail_quantile_sample_requirement"
                ]["observed_positive_weight_particle_count"]
            )
            for record in ordered
        ],
        dtype=np.int64,
    )
    fronts = np.asarray(
        [
            float(record["mhd"]["detected_front"]["x_front_c_over_omega_pi"])
            for record in ordered
        ],
        dtype=np.float64,
    )
    compression = np.asarray(
        [
            float(
                record["mhd"]["shock_state"][
                    "compression_ratio_downstream_over_upstream"
                ]
            )
            for record in ordered
        ],
        dtype=np.float64,
    )
    _require(
        np.all(np.isfinite(max_chi))
        and np.all(max_chi >= 0.0)
        and np.all(np.isfinite(max_specific_energy))
        and np.all(max_specific_energy >= 0.0)
        and all(
            np.all(np.isfinite(values)) and np.all(values >= 0.0)
            for values in tail_quantiles.values()
        )
        and all(
            np.all(np.isfinite(values)) and np.all(values >= 0.0)
            for values in tail_specific_energy_quantiles.values()
        )
        and np.all(tail_sample_counts >= PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES)
        and np.all(np.isfinite(fronts))
        and np.all(np.isfinite(compression))
        and np.all(compression > 0.0),
        "production-science history contains invalid diagnostics",
    )
    dt = np.diff(times)
    interval_rows = [
        {
            "time_start_omega0_inverse": float(times[index]),
            "time_stop_omega0_inverse": float(times[index + 1]),
            "delta_time_omega0_inverse": float(dt[index]),
            "delta_maximum_chi_over_delta_time": float(
                (max_chi[index + 1] - max_chi[index]) / dt[index]
            ),
            "delta_maximum_specific_kinetic_energy_over_delta_time": float(
                (max_specific_energy[index + 1] - max_specific_energy[index]) / dt[index]
            ),
            "measured_shock_speed_over_u_a0": float(
                (fronts[index + 1] - fronts[index]) / dt[index]
            ),
        }
        for index in range(dt.size)
    ]
    spectrum_time_series = []
    seen_spectrum_slots: set[float] = set()
    for record in ordered:
        nominal = _finite_scalar(record["nominal_slot_time"], "history nominal slot time")
        expected = nominal in PREREGISTERED_SPECTRUM_NOMINAL_TIMES
        present = "preregistered_downstream_spectrum" in record["particles"]
        _require(
            present == expected,
            "preregistered downstream spectrum presence drifted from its nominal slots",
        )
        if not present:
            continue
        _require(nominal not in seen_spectrum_slots, "duplicate preregistered spectrum slot")
        seen_spectrum_slots.add(nominal)
        spectrum = record["particles"]["preregistered_downstream_spectrum"]
        _require(
            spectrum["reuse_contract"]
            == "q011_section54_particles.reduce_particle_snapshot"
            and float(spectrum["nominal_slot_time"]) == nominal,
            "preregistered downstream spectrum reducer binding drifted",
        )
        reduction = spectrum["reduction"]
        _require(
            reduction["record_type"] == "q011_section54_particle_snapshot_reduction",
            "preregistered downstream spectrum record type drifted",
        )
        _require(
            int(reduction["weighted_spectrum"]["particle_count_post_filter"])
            >= DOWNSTREAM_SPECTRUM_MINIMUM_PARTICLE_SAMPLES
            and float(reduction["weighted_spectrum"]["total_post_filter_macro_weight"])
            > 0.0,
            "preregistered downstream spectrum history has insufficient admitted samples",
        )
        if nominal == particle_primitives.LATE_SLOPE_SNAPSHOT_TIME:
            late_slope = reduction["late_slope"]
            _require(
                int(late_slope["positive_fit_bin_count"])
                >= particle_primitives.LATE_SLOPE_MINIMUM_POSITIVE_BINS,
                "late downstream spectrum slope has insufficient positive fit bins",
            )
        spectrum_time_series.append(
            {
                "observed_committed_time": float(record["observed_committed_time"]),
                **spectrum,
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_production_science_history_v1",
        "successor_id": SUCCESSOR_ID,
        "observed_committed_times_omega0_inverse": times.tolist(),
        "fit_requirements": {
            "linear_fit_minimum_snapshot_count": (
                HISTORY_LINEAR_FIT_MINIMUM_SNAPSHOTS
            ),
            "adjacent_interval_rates_are_two_snapshot_secants_not_regression_fits": True,
        },
        "particle_acceleration": {
            "maximum_chi": max_chi.tolist(),
            "maximum_specific_kinetic_energy": max_specific_energy.tolist(),
            "maximum_chi_linear_fit": _linear_fit(
                times, max_chi, label="maximum chi acceleration"
            ),
            "maximum_specific_kinetic_energy_linear_fit": _linear_fit(
                times, max_specific_energy, label="maximum specific energy acceleration"
            ),
            "macro_weighted_chi_tail_quantile_history": {
                **{label: values.tolist() for label, values in tail_quantiles.items()},
                "minimum_positive_weight_particle_count_per_snapshot": (
                    PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES
                ),
                "observed_positive_weight_particle_counts": tail_sample_counts.tolist(),
            },
            "macro_weighted_specific_kinetic_energy_tail_quantile_history": {
                **{
                    label: values.tolist()
                    for label, values in tail_specific_energy_quantiles.items()
                },
                "minimum_positive_weight_particle_count_per_snapshot": (
                    PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES
                ),
                "observed_positive_weight_particle_counts": tail_sample_counts.tolist(),
            },
            "float32_pvtk_velocity_uncertainty_time_series": [
                {
                    "observed_committed_time": float(record["observed_committed_time"]),
                    **record["particles"]["float32_pvtk_velocity_uncertainty"],
                }
                for record in ordered
            ],
            "nondecreasing_interval_fraction": float(
                np.count_nonzero(np.diff(max_chi) >= 0.0) / dt.size
            ),
        },
        "shock_kinematics": {
            "negative_gradient_detector_preregistration_status": (
                "requires_preregistration_supersession_before_qualifying_use"
            ),
            "detected_front_x1_c_over_omega_pi": fronts.tolist(),
            "compression_ratio_downstream_over_upstream": compression.tolist(),
            "measured_wall_frame_shock_speed_linear_fit": _linear_fit(
                times, fronts, label="measured shock speed"
            ),
        },
        "interval_rates": interval_rows,
        "energy_partition_time_series": [
            {
                "observed_committed_time": float(record["observed_committed_time"]),
                "full_domain": record["energy_partition"]["full_domain"],
                "shock_centered_window": record["energy_partition"][
                    "shock_centered_window"
                ],
            }
            for record in ordered
        ],
        "upstream_magnetic_amplification_time_series": [
            {
                "observed_committed_time": float(record["observed_committed_time"]),
                **record["mhd"]["upstream_magnetic_amplification"],
            }
            for record in ordered
        ],
        "cr_current_statistics_time_series": [
            {
                "observed_committed_time": float(record["observed_committed_time"]),
                "statistics_by_region": record["cr_current"]["statistics_by_region"],
            }
            for record in ordered
        ],
        "preregistered_downstream_spectrum_time_series": spectrum_time_series,
        "spectrum_requirements": {
            "required_nominal_slots_omega0_inverse": list(
                PREREGISTERED_SPECTRUM_NOMINAL_TIMES
            ),
            "reducer": "q011_section54_particles.reduce_particle_snapshot",
            "downstream_classifier": "ideal_injection_surface",
            "minimum_positive_fit_bins_at_t1200": (
                particle_primitives.LATE_SLOPE_MINIMUM_POSITIVE_BINS
            ),
            "minimum_admitted_particle_count_per_spectrum": (
                DOWNSTREAM_SPECTRUM_MINIMUM_PARTICLE_SAMPLES
            ),
        },
    }


@_public_contract("canonical diagnostic serialization")
def canonical_record_bytes(record: Mapping[str, object]) -> bytes:
    """Return deterministic JSON bytes for one successor diagnostic record."""
    _require(isinstance(record, Mapping), "diagnostic record must be a mapping")
    return (
        json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


__all__ = [
    "BASELINE_DECK",
    "BIRTH_TIME_MINIMUM",
    "ComposedMHDState",
    "CURRENT_PRODUCT_FIELDS",
    "DOWNSTREAM_SPECTRUM_MINIMUM_PARTICLE_SAMPLES",
    "HISTORY_LINEAR_FIT_MINIMUM_SNAPSHOTS",
    "MATCHED_SNAPSHOT_FIRST_ELIGIBLE_NOMINAL_TIME",
    "MHD_PRIMITIVE_FIELDS",
    "PARTICLE_MACRO_MASS",
    "PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES",
    "ProductionScienceError",
    "REFERENCE_B0",
    "SUCCESSOR_DECK",
    "SUCCESSOR_ID",
    "canonical_record_bytes",
    "compose_full_mhd_state",
    "reduce_cr_current_snapshot",
    "reduce_mhd_snapshot",
    "reduce_particle_energy_snapshot",
    "reduce_production_science_history",
    "reduce_production_science_snapshot",
    "validate_successor_deck",
]
