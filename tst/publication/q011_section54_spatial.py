#!/usr/bin/env python3
"""Pure spatial reducers for the bounded Q-011 Section 5.4 t=500 tranche.

This module accepts already parsed AthenaK mesh datasets.  It does not bind
campaign artifacts, load policy files, or decide whether a campaign closes a
qualification claim.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math
from numbers import Real
from types import MappingProxyType
from typing import Any

import numpy as np

if __package__:
    from . import analyze_q011_section54_outputs as output_primitives
else:
    import analyze_q011_section54_outputs as output_primitives


AnalysisError = output_primitives.AnalysisError

MESH_QUANTITY_FIELDS: Mapping[str, str] = MappingProxyType(
    {
        "rho": "dens",
        "bmag": "bmag",
        "prtcl_jx": "prtcl_jx",
        "j2": "j2",
    }
)
REQUIRED_MESH_QUANTITIES = tuple(MESH_QUANTITY_FIELDS)
T500_OMEGA0_INVERSE = 500.0
SHOCK_SEARCH_OFFSETS_C_OVER_OMEGA_PI = (-1200.0, 1200.0)
MAX_FRONT_OFFSET_C_OVER_OMEGA_PI = 600.0
UPSTREAM_B_OFFSETS_C_OVER_OMEGA_PI = (120.0, 1200.0)
REFERENCE_B0 = 1.0
UPSTREAM_B0_GATE = (1.2, 3.5)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AnalysisError(message)


def _finite_real(value: object, label: str) -> float:
    _require(
        isinstance(value, Real) and not isinstance(value, (bool, np.bool_)),
        f"{label} must be a real number",
    )
    result = float(value)
    _require(math.isfinite(result), f"{label} must be finite")
    return result


def _snapshot_times(
    nominal_slot_time: object,
    observed_committed_time: object,
    *,
    label: str,
) -> tuple[float, float]:
    nominal = _finite_real(nominal_slot_time, f"{label} nominal slot time")
    observed = _finite_real(
        observed_committed_time, f"{label} observed committed time"
    )
    _require(
        nominal >= 0.0 and observed >= 0.0,
        f"{label} times must be non-negative",
    )
    return nominal, observed


def _mesh_time_projection(observed_committed_time: float) -> float:
    return float(format(observed_committed_time, ".6g"))


def _readonly_finite_array(values: object, label: str) -> np.ndarray:
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise AnalysisError(f"{label} must be numeric") from error
    _require(array.size > 0, f"{label} must not be empty")
    _require(np.all(np.isfinite(array)), f"{label} must be finite")
    result = np.array(array, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result


def _readonly_source_levels(values: object, label: str) -> np.ndarray:
    try:
        array = np.asarray(values)
    except (TypeError, ValueError) as error:
        raise AnalysisError(f"{label} must be an integer array") from error
    _require(array.size > 0, f"{label} must not be empty")
    _require(
        np.issubdtype(array.dtype, np.integer)
        and not np.issubdtype(array.dtype, np.bool_),
        f"{label} must be an integer array",
    )
    result = np.array(array, dtype=np.int64, copy=True)
    result.setflags(write=False)
    return result


def _readonly_faces(values: object, label: str) -> np.ndarray:
    faces = _readonly_finite_array(values, label)
    _require(faces.ndim == 1 and faces.size >= 2, f"{label} shape is invalid")
    _require(np.all(np.diff(faces) > 0.0), f"{label} must increase")
    return faces


def mesh_field_name(quantity: object) -> str:
    """Return the one permitted AthenaK field for a retained mesh quantity."""
    _require(type(quantity) is str, "mesh quantity must be text")
    _require(quantity in MESH_QUANTITY_FIELDS, f"unknown mesh quantity {quantity!r}")
    return MESH_QUANTITY_FIELDS[quantity]


def _finite_window(value: object, label: str) -> tuple[float, float]:
    _require(
        isinstance(value, tuple) and len(value) == 2,
        f"{label} must be a two-item tuple",
    )
    lower = _finite_real(value[0], f"{label} lower bound")
    upper = _finite_real(value[1], f"{label} upper bound")
    _require(upper > lower, f"{label} must increase")
    return (lower, upper)


@dataclass(frozen=True)
class CartesianXYRaster:
    """One finite Cartesian x-y raster after collapsing the singleton x3 axis."""

    quantity: str
    source_field: str
    nominal_slot_time: float
    observed_committed_time: float
    x1_faces_c_over_omega_pi: np.ndarray
    x2_faces_c_over_omega_pi: np.ndarray
    collapsed_x3_faces: np.ndarray
    values_y_x: np.ndarray
    source_levels_y_x: np.ndarray
    target_level: int

    def __post_init__(self) -> None:
        expected_field = mesh_field_name(self.quantity)
        _require(
            self.source_field == expected_field,
            f"{self.quantity}: source field must be {expected_field!r}",
        )
        nominal, observed = _snapshot_times(
            self.nominal_slot_time,
            self.observed_committed_time,
            label="raster",
        )
        x1_faces = _readonly_faces(
            self.x1_faces_c_over_omega_pi, "raster x1 faces"
        )
        x2_faces = _readonly_faces(
            self.x2_faces_c_over_omega_pi, "raster x2 faces"
        )
        x3_faces = _readonly_faces(self.collapsed_x3_faces, "raster x3 faces")
        _require(x3_faces.size == 2, "raster x3 axis must contain exactly one cell")
        values = _readonly_finite_array(self.values_y_x, "raster values")
        levels = _readonly_source_levels(
            self.source_levels_y_x, "raster source levels"
        )
        expected_shape = (x2_faces.size - 1, x1_faces.size - 1)
        _require(values.ndim == 2, "raster values must be a two-dimensional y-x array")
        _require(values.shape == expected_shape, "raster values shape disagrees with faces")
        _require(levels.shape == values.shape, "raster source-level shape disagrees")
        _require(
            isinstance(self.target_level, int)
            and not isinstance(self.target_level, bool)
            and self.target_level >= 0,
            "raster target level must be a non-negative integer",
        )
        _require(np.all(levels >= 0), "raster source levels must be non-negative")
        _require(
            np.all(levels <= self.target_level),
            "raster source level exceeds target level",
        )
        object.__setattr__(self, "nominal_slot_time", nominal)
        object.__setattr__(self, "observed_committed_time", observed)
        object.__setattr__(self, "x1_faces_c_over_omega_pi", x1_faces)
        object.__setattr__(self, "x2_faces_c_over_omega_pi", x2_faces)
        object.__setattr__(self, "collapsed_x3_faces", x3_faces)
        object.__setattr__(self, "values_y_x", values)
        object.__setattr__(self, "source_levels_y_x", levels)


@dataclass(frozen=True)
class YAreaWeightedProfile:
    """One finite x profile formed by conservative y-area weighting."""

    quantity: str
    source_field: str
    nominal_slot_time: float
    observed_committed_time: float
    x1_centers_c_over_omega_pi: np.ndarray
    values_x: np.ndarray
    column_areas: np.ndarray

    def __post_init__(self) -> None:
        expected_field = mesh_field_name(self.quantity)
        _require(
            self.source_field == expected_field,
            f"{self.quantity}: profile source field must be {expected_field!r}",
        )
        nominal, observed = _snapshot_times(
            self.nominal_slot_time,
            self.observed_committed_time,
            label="profile",
        )
        x1 = _readonly_finite_array(
            self.x1_centers_c_over_omega_pi, "profile x1 centers"
        )
        values = _readonly_finite_array(self.values_x, "profile values")
        areas = _readonly_finite_array(self.column_areas, "profile column areas")
        _require(
            x1.ndim == values.ndim == areas.ndim == 1,
            "profile arrays must be one-dimensional",
        )
        _require(
            x1.size == values.size == areas.size,
            "profile array sizes disagree",
        )
        _require(np.all(np.diff(x1) > 0.0), "profile x1 centers must increase")
        _require(np.all(areas > 0.0), "profile column areas must be positive")
        object.__setattr__(self, "nominal_slot_time", nominal)
        object.__setattr__(self, "observed_committed_time", observed)
        object.__setattr__(self, "x1_centers_c_over_omega_pi", x1)
        object.__setattr__(self, "values_x", values)
        object.__setattr__(self, "column_areas", areas)


@dataclass(frozen=True)
class DetectedFrontRecord:
    """A unique strongest positive density gradient near the ideal surface."""

    nominal_slot_time: float
    observed_committed_time: float
    x_ideal_c_over_omega_pi: float
    search_window_c_over_omega_pi: tuple[float, float]
    front_index: int
    x_front_c_over_omega_pi: float
    density_gradient: float
    offset_from_x_ideal_c_over_omega_pi: float
    max_absolute_offset_c_over_omega_pi: float

    def __post_init__(self) -> None:
        nominal, observed = _snapshot_times(
            self.nominal_slot_time,
            self.observed_committed_time,
            label="detected front",
        )
        x_ideal = _finite_real(self.x_ideal_c_over_omega_pi, "ideal shock position")
        window = _finite_window(
            self.search_window_c_over_omega_pi, "shock search window"
        )
        _require(
            window
            == tuple(
                x_ideal + offset
                for offset in SHOCK_SEARCH_OFFSETS_C_OVER_OMEGA_PI
            ),
            "shock search window differs from the fixed x_ideal-relative window",
        )
        _require(
            isinstance(self.front_index, int)
            and not isinstance(self.front_index, bool)
            and self.front_index >= 0,
            "detected front index must be a non-negative integer",
        )
        x_front = _finite_real(self.x_front_c_over_omega_pi, "detected front position")
        gradient = _finite_real(self.density_gradient, "detected density gradient")
        _require(gradient > 0.0, "detected density gradient must be positive")
        offset = _finite_real(
            self.offset_from_x_ideal_c_over_omega_pi,
            "detected front offset",
        )
        _require(offset == x_front - x_ideal, "detected front offset is inconsistent")
        max_offset = _finite_real(
            self.max_absolute_offset_c_over_omega_pi,
            "maximum detected front offset",
        )
        _require(
            max_offset == MAX_FRONT_OFFSET_C_OVER_OMEGA_PI,
            "maximum detected front offset differs from the fixed bound",
        )
        _require(abs(offset) <= max_offset, "detected shock front exceeds its offset bound")
        object.__setattr__(self, "nominal_slot_time", nominal)
        object.__setattr__(self, "observed_committed_time", observed)
        object.__setattr__(self, "x_ideal_c_over_omega_pi", x_ideal)
        object.__setattr__(self, "search_window_c_over_omega_pi", window)
        object.__setattr__(self, "x_front_c_over_omega_pi", x_front)
        object.__setattr__(self, "density_gradient", gradient)
        object.__setattr__(self, "offset_from_x_ideal_c_over_omega_pi", offset)
        object.__setattr__(self, "max_absolute_offset_c_over_omega_pi", max_offset)


@dataclass(frozen=True)
class UpstreamBAmplificationRecord:
    """Area-weighted upstream magnetic amplification and its fixed t=500 gate."""

    nominal_slot_time: float
    observed_committed_time: float
    x_ideal_c_over_omega_pi: float
    upstream_window_c_over_omega_pi: tuple[float, float]
    selected_cell_count: int
    selected_area: float
    mean_magnetic_magnitude: float
    reference_b0: float
    amplification_over_b0: float
    acceptance_range: tuple[float, float]
    passes_gate: bool

    def __post_init__(self) -> None:
        nominal, observed = _snapshot_times(
            self.nominal_slot_time,
            self.observed_committed_time,
            label="amplification",
        )
        _require(
            nominal == T500_OMEGA0_INVERSE,
            "upstream magnetic-amplification gate is defined only for nominal t=500",
        )
        x_ideal = _finite_real(self.x_ideal_c_over_omega_pi, "ideal shock position")
        window = _finite_window(
            self.upstream_window_c_over_omega_pi, "upstream magnetic window"
        )
        _require(
            window
            == tuple(
                x_ideal + offset
                for offset in UPSTREAM_B_OFFSETS_C_OVER_OMEGA_PI
            ),
            "upstream magnetic window differs from the fixed x_ideal-relative window",
        )
        _require(
            isinstance(self.selected_cell_count, int)
            and not isinstance(self.selected_cell_count, bool)
            and self.selected_cell_count > 0,
            "upstream selected-cell count must be a positive integer",
        )
        selected_area = _finite_real(self.selected_area, "upstream selected area")
        _require(selected_area > 0.0, "upstream selected area must be positive")
        mean = _finite_real(
            self.mean_magnetic_magnitude, "upstream mean magnetic magnitude"
        )
        _require(mean >= 0.0, "upstream mean magnetic magnitude must be non-negative")
        reference = _finite_real(self.reference_b0, "reference B0")
        _require(reference == REFERENCE_B0, "reference B0 differs from the fixed value")
        amplification = _finite_real(self.amplification_over_b0, "upstream B amplification")
        _require(
            amplification == mean / reference,
            "upstream B amplification is inconsistent with its mean and reference",
        )
        gate = _finite_window(self.acceptance_range, "upstream B0 gate")
        _require(gate == UPSTREAM_B0_GATE, "upstream B0 gate differs from the fixed range")
        expected_pass = gate[0] <= amplification <= gate[1]
        _require(type(self.passes_gate) is bool, "upstream B0 gate result must be boolean")
        _require(self.passes_gate == expected_pass, "upstream B0 gate result is inconsistent")
        object.__setattr__(self, "nominal_slot_time", nominal)
        object.__setattr__(self, "observed_committed_time", observed)
        object.__setattr__(self, "x_ideal_c_over_omega_pi", x_ideal)
        object.__setattr__(self, "upstream_window_c_over_omega_pi", window)
        object.__setattr__(self, "selected_area", selected_area)
        object.__setattr__(self, "mean_magnetic_magnitude", mean)
        object.__setattr__(self, "reference_b0", reference)
        object.__setattr__(self, "amplification_over_b0", amplification)
        object.__setattr__(self, "acceptance_range", gate)


def compose_xy_quantity(
    dataset: output_primitives.AthenaBinaryDataset,
    quantity: object,
    *,
    nominal_slot_time: object | None = None,
    observed_committed_time: object | None = None,
    target_level: int | None = None,
) -> CartesianXYRaster:
    """Compose one exact scalar mesh product and collapse its singleton x3 axis."""
    field = mesh_field_name(quantity)
    _require(
        isinstance(dataset, output_primitives.AthenaBinaryDataset),
        f"{quantity}: expected an AthenaBinaryDataset",
    )
    _require(
        len(dataset.root_grid_shape) == 3 and len(dataset.meshblock_shape) == 3,
        f"{quantity}: mesh shape metadata must contain exactly three axes",
    )
    _require(
        dataset.variable_names == (field,),
        f"{quantity}: retained scalar product must contain exactly field {field!r}",
    )
    _require(
        dataset.root_grid_shape[2] == 1 and dataset.meshblock_shape[2] == 1,
        f"{quantity}: x3 must contain exactly one cell",
    )
    if nominal_slot_time is None and observed_committed_time is None:
        nominal_slot_time = dataset.time
        observed_committed_time = dataset.time
    else:
        _require(
            nominal_slot_time is not None and observed_committed_time is not None,
            f"{quantity}: nominal and observed snapshot times must be supplied together",
        )
    nominal, observed = _snapshot_times(
        nominal_slot_time,
        observed_committed_time,
        label=f"{quantity} snapshot",
    )
    _require(
        dataset.time == _mesh_time_projection(observed),
        f"{quantity}: embedded mesh time differs from observed committed time projection",
    )
    composite = output_primitives.compose_leaf_field(
        dataset, field, target_level=target_level
    )
    _require(
        composite.values.ndim == 3
        and composite.values.shape[0] == 1
        and composite.source_levels.shape == composite.values.shape,
        f"{quantity}: composed x3 axis must contain exactly one cell",
    )
    return CartesianXYRaster(
        quantity=quantity,
        source_field=field,
        nominal_slot_time=nominal,
        observed_committed_time=observed,
        x1_faces_c_over_omega_pi=composite.x1_faces,
        x2_faces_c_over_omega_pi=composite.x2_faces,
        collapsed_x3_faces=composite.x3_faces,
        values_y_x=composite.values[0],
        source_levels_y_x=composite.source_levels[0],
        target_level=composite.target_level,
    )


def compose_xy_snapshot(
    datasets: Mapping[str, output_primitives.AthenaBinaryDataset],
    *,
    nominal_slot_time: object | None = None,
    observed_committed_time: object | None = None,
    target_level: int | None = None,
) -> dict[str, CartesianXYRaster]:
    """Compose the exact four retained scalar products onto one Cartesian grid."""
    _require(isinstance(datasets, Mapping), "snapshot datasets must be a mapping")
    _require(
        set(datasets) == set(REQUIRED_MESH_QUANTITIES),
        "snapshot mesh quantities must be exactly rho, bmag, prtcl_jx, and j2",
    )
    rasters = {
        quantity: compose_xy_quantity(
            datasets[quantity],
            quantity,
            nominal_slot_time=nominal_slot_time,
            observed_committed_time=observed_committed_time,
            target_level=target_level,
        )
        for quantity in REQUIRED_MESH_QUANTITIES
    }
    reference = rasters[REQUIRED_MESH_QUANTITIES[0]]
    for quantity in REQUIRED_MESH_QUANTITIES[1:]:
        raster = rasters[quantity]
        _require(
            raster.nominal_slot_time == reference.nominal_slot_time
            and raster.observed_committed_time == reference.observed_committed_time,
            "snapshot mesh-product times disagree",
        )
        _require(
            raster.target_level == reference.target_level
            and np.array_equal(
                raster.x1_faces_c_over_omega_pi,
                reference.x1_faces_c_over_omega_pi,
            )
            and np.array_equal(
                raster.x2_faces_c_over_omega_pi,
                reference.x2_faces_c_over_omega_pi,
            )
            and np.array_equal(raster.collapsed_x3_faces, reference.collapsed_x3_faces)
            and np.array_equal(raster.source_levels_y_x, reference.source_levels_y_x),
            "snapshot mesh-product Cartesian grids disagree",
        )
    return rasters


def x1_cell_centers(raster: CartesianXYRaster) -> np.ndarray:
    """Return finite x1 centers for one validated Cartesian raster."""
    _require(isinstance(raster, CartesianXYRaster), "expected a CartesianXYRaster")
    return 0.5 * (
        raster.x1_faces_c_over_omega_pi[:-1]
        + raster.x1_faces_c_over_omega_pi[1:]
    )


def xy_cell_areas(raster: CartesianXYRaster) -> np.ndarray:
    """Return Cartesian cell areas in the raster's y-x order."""
    _require(isinstance(raster, CartesianXYRaster), "expected a CartesianXYRaster")
    return np.diff(raster.x2_faces_c_over_omega_pi)[:, None] * np.diff(
        raster.x1_faces_c_over_omega_pi
    )[None, :]


def y_area_weighted_profile(raster: CartesianXYRaster) -> YAreaWeightedProfile:
    """Collapse one x-y raster to a conservative area-weighted x profile."""
    _require(isinstance(raster, CartesianXYRaster), "expected a CartesianXYRaster")
    areas = xy_cell_areas(raster)
    column_areas = np.sum(areas, axis=0)
    profile = np.sum(raster.values_y_x * areas, axis=0) / column_areas
    return YAreaWeightedProfile(
        quantity=raster.quantity,
        source_field=raster.source_field,
        nominal_slot_time=raster.nominal_slot_time,
        observed_committed_time=raster.observed_committed_time,
        x1_centers_c_over_omega_pi=x1_cell_centers(raster),
        values_x=profile,
        column_areas=column_areas,
    )


def _translated_window(
    x_ideal_c_over_omega_pi: object,
    offsets: tuple[float, float],
) -> tuple[float, float]:
    x_ideal = _finite_real(x_ideal_c_over_omega_pi, "ideal shock position")
    window = (x_ideal + offsets[0], x_ideal + offsets[1])
    _require(all(math.isfinite(value) for value in window), "translated window is not finite")
    return window


def ideal_shock_search_window(
    x_ideal_c_over_omega_pi: object,
) -> tuple[float, float]:
    """Return the fixed detected-front search window around ``x_ideal``."""
    return _translated_window(
        x_ideal_c_over_omega_pi, SHOCK_SEARCH_OFFSETS_C_OVER_OMEGA_PI
    )


def upstream_b_window(x_ideal_c_over_omega_pi: object) -> tuple[float, float]:
    """Return the fixed open-upstream magnetic-amplification window."""
    return _translated_window(
        x_ideal_c_over_omega_pi, UPSTREAM_B_OFFSETS_C_OVER_OMEGA_PI
    )


def detect_positive_gradient_front(
    rho_profile: YAreaWeightedProfile,
    *,
    x_ideal_c_over_omega_pi: object,
) -> DetectedFrontRecord:
    """Detect and bound the unique strongest positive density gradient."""
    _require(
        isinstance(rho_profile, YAreaWeightedProfile),
        "expected a YAreaWeightedProfile",
    )
    _require(
        rho_profile.quantity == "rho" and rho_profile.source_field == "dens",
        "detected front requires the rho -> dens profile",
    )
    x_ideal = _finite_real(x_ideal_c_over_omega_pi, "ideal shock position")
    window = ideal_shock_search_window(x_ideal)
    front = output_primitives.detect_shock_front(
        rho_profile.x1_centers_c_over_omega_pi,
        rho_profile.values_x,
        search_window=window,
        gradient_sign="positive",
    )
    offset = front.x1 - x_ideal
    _require(
        abs(offset) <= MAX_FRONT_OFFSET_C_OVER_OMEGA_PI,
        "detected shock front exceeds the maximum absolute offset from x_ideal",
    )
    return DetectedFrontRecord(
        nominal_slot_time=rho_profile.nominal_slot_time,
        observed_committed_time=rho_profile.observed_committed_time,
        x_ideal_c_over_omega_pi=x_ideal,
        search_window_c_over_omega_pi=window,
        front_index=front.index,
        x_front_c_over_omega_pi=front.x1,
        density_gradient=front.density_gradient,
        offset_from_x_ideal_c_over_omega_pi=offset,
        max_absolute_offset_c_over_omega_pi=MAX_FRONT_OFFSET_C_OVER_OMEGA_PI,
    )


def reduce_upstream_b_amplification_at_t500(
    bmag_raster: CartesianXYRaster,
    *,
    x_ideal_c_over_omega_pi: object,
) -> UpstreamBAmplificationRecord:
    """Reduce ``|B| / B0`` in the fixed t=500 open-upstream window."""
    _require(isinstance(bmag_raster, CartesianXYRaster), "expected a CartesianXYRaster")
    _require(
        bmag_raster.quantity == "bmag" and bmag_raster.source_field == "bmag",
        "upstream amplification requires the bmag -> bmag raster",
    )
    _require(
        bmag_raster.nominal_slot_time == T500_OMEGA0_INVERSE,
        "upstream magnetic-amplification gate is defined only for nominal t=500",
    )
    x_ideal = _finite_real(x_ideal_c_over_omega_pi, "ideal shock position")
    window = upstream_b_window(x_ideal)
    estimate = output_primitives.estimate_upstream_amplification(
        x1_cell_centers(bmag_raster),
        bmag_raster.values_y_x,
        upstream_window=window,
        reference_magnetic_field=REFERENCE_B0,
        cell_areas=xy_cell_areas(bmag_raster),
    )
    lower, upper = UPSTREAM_B0_GATE
    return UpstreamBAmplificationRecord(
        nominal_slot_time=bmag_raster.nominal_slot_time,
        observed_committed_time=bmag_raster.observed_committed_time,
        x_ideal_c_over_omega_pi=x_ideal,
        upstream_window_c_over_omega_pi=window,
        selected_cell_count=estimate.selected_cell_count,
        selected_area=estimate.selected_area,
        mean_magnetic_magnitude=estimate.mean_magnetic_magnitude,
        reference_b0=REFERENCE_B0,
        amplification_over_b0=estimate.amplification,
        acceptance_range=UPSTREAM_B0_GATE,
        passes_gate=lower <= estimate.amplification <= upper,
    )


def morphology_raster_record(raster: CartesianXYRaster) -> dict[str, Any]:
    """Return a deterministic JSON-ready raster with source-level overlay data."""
    _require(isinstance(raster, CartesianXYRaster), "expected a CartesianXYRaster")
    levels, counts = np.unique(raster.source_levels_y_x, return_counts=True)
    return {
        "schema_version": 1,
        "record_type": "q011_section54_morphology_raster",
        "quantity": raster.quantity,
        "source_field": raster.source_field,
        "nominal_slot_time": raster.nominal_slot_time,
        "observed_committed_time": raster.observed_committed_time,
        "cartesian_axes": ["x2", "x1"],
        "shape_y_x": list(raster.values_y_x.shape),
        "x1_faces_c_over_omega_pi": raster.x1_faces_c_over_omega_pi.tolist(),
        "x2_faces_c_over_omega_pi": raster.x2_faces_c_over_omega_pi.tolist(),
        "collapsed_x3": {
            "cell_count": 1,
            "faces": raster.collapsed_x3_faces.tolist(),
        },
        "values_y_x": raster.values_y_x.tolist(),
        "source_level_overlay": {
            "encoding": "physical_refinement_level_per_cell",
            "target_composite_level": raster.target_level,
            "levels_y_x": raster.source_levels_y_x.tolist(),
            "level_cell_counts": [
                {"source_level": int(level), "cell_count": int(count)}
                for level, count in zip(levels, counts)
            ],
        },
    }


def y_area_weighted_profile_record(profile: YAreaWeightedProfile) -> dict[str, Any]:
    """Return a deterministic JSON-ready y-area-weighted profile."""
    _require(isinstance(profile, YAreaWeightedProfile), "expected a YAreaWeightedProfile")
    return {
        "quantity": profile.quantity,
        "source_field": profile.source_field,
        "nominal_slot_time": profile.nominal_slot_time,
        "observed_committed_time": profile.observed_committed_time,
        "x1_centers_c_over_omega_pi": profile.x1_centers_c_over_omega_pi.tolist(),
        "values_x": profile.values_x.tolist(),
        "column_areas": profile.column_areas.tolist(),
    }


def detected_front_record(front: DetectedFrontRecord) -> dict[str, Any]:
    """Return a deterministic JSON-ready detected-front record."""
    _require(isinstance(front, DetectedFrontRecord), "expected a DetectedFrontRecord")
    return {
        "nominal_slot_time": front.nominal_slot_time,
        "observed_committed_time": front.observed_committed_time,
        "x_ideal_c_over_omega_pi": front.x_ideal_c_over_omega_pi,
        "search_window_c_over_omega_pi": list(front.search_window_c_over_omega_pi),
        "detector": "unique_strongest_positive_density_gradient",
        "front_index": front.front_index,
        "x_front_c_over_omega_pi": front.x_front_c_over_omega_pi,
        "density_gradient": front.density_gradient,
        "offset_from_x_ideal_c_over_omega_pi": (
            front.offset_from_x_ideal_c_over_omega_pi
        ),
        "max_absolute_offset_c_over_omega_pi": (
            front.max_absolute_offset_c_over_omega_pi
        ),
    }


def upstream_b_amplification_record(
    amplification: UpstreamBAmplificationRecord,
) -> dict[str, Any]:
    """Return a deterministic JSON-ready upstream-amplification record."""
    _require(
        isinstance(amplification, UpstreamBAmplificationRecord),
        "expected an UpstreamBAmplificationRecord",
    )
    return {
        "nominal_slot_time": amplification.nominal_slot_time,
        "observed_committed_time": amplification.observed_committed_time,
        "x_ideal_c_over_omega_pi": amplification.x_ideal_c_over_omega_pi,
        "upstream_window_c_over_omega_pi": list(
            amplification.upstream_window_c_over_omega_pi
        ),
        "selected_cell_count": amplification.selected_cell_count,
        "selected_area": amplification.selected_area,
        "mean_magnetic_magnitude": amplification.mean_magnetic_magnitude,
        "reference_b0": amplification.reference_b0,
        "amplification_over_b0": amplification.amplification_over_b0,
        "acceptance_range": list(amplification.acceptance_range),
        "passes_gate": amplification.passes_gate,
    }


def reduce_t500_spatial_snapshot(
    datasets: Mapping[str, output_primitives.AthenaBinaryDataset],
    *,
    nominal_slot_time: object,
    observed_committed_time: object,
    x_ideal_c_over_omega_pi: object,
    target_level: int | None = None,
) -> dict[str, Any]:
    """Reduce one exact t=500 four-product spatial snapshot to archive records."""
    nominal, observed = _snapshot_times(
        nominal_slot_time,
        observed_committed_time,
        label="bounded spatial snapshot",
    )
    rasters = compose_xy_snapshot(
        datasets,
        nominal_slot_time=nominal,
        observed_committed_time=observed,
        target_level=target_level,
    )
    _require(
        rasters["rho"].nominal_slot_time == T500_OMEGA0_INVERSE,
        "bounded spatial snapshot reducer is defined only for nominal t=500",
    )
    profiles = {
        quantity: y_area_weighted_profile(raster)
        for quantity, raster in rasters.items()
    }
    front = detect_positive_gradient_front(
        profiles["rho"], x_ideal_c_over_omega_pi=x_ideal_c_over_omega_pi
    )
    amplification = reduce_upstream_b_amplification_at_t500(
        rasters["bmag"], x_ideal_c_over_omega_pi=x_ideal_c_over_omega_pi
    )
    return {
        "schema_version": 1,
        "record_type": "q011_section54_t500_spatial_reduction",
        "nominal_slot_time": nominal,
        "observed_committed_time": observed,
        "x_ideal_c_over_omega_pi": _finite_real(
            x_ideal_c_over_omega_pi, "ideal shock position"
        ),
        "mesh_quantity_fields": dict(MESH_QUANTITY_FIELDS),
        "detected_front": detected_front_record(front),
        "upstream_b_amplification": upstream_b_amplification_record(amplification),
        "y_area_weighted_profiles": {
            quantity: y_area_weighted_profile_record(profiles[quantity])
            for quantity in REQUIRED_MESH_QUANTITIES
        },
        "morphology_rasters": {
            quantity: morphology_raster_record(rasters[quantity])
            for quantity in REQUIRED_MESH_QUANTITIES
        },
    }


__all__ = [
    "AnalysisError",
    "CartesianXYRaster",
    "DetectedFrontRecord",
    "MAX_FRONT_OFFSET_C_OVER_OMEGA_PI",
    "MESH_QUANTITY_FIELDS",
    "REFERENCE_B0",
    "REQUIRED_MESH_QUANTITIES",
    "SHOCK_SEARCH_OFFSETS_C_OVER_OMEGA_PI",
    "T500_OMEGA0_INVERSE",
    "UPSTREAM_B0_GATE",
    "UPSTREAM_B_OFFSETS_C_OVER_OMEGA_PI",
    "UpstreamBAmplificationRecord",
    "YAreaWeightedProfile",
    "compose_xy_quantity",
    "compose_xy_snapshot",
    "detect_positive_gradient_front",
    "detected_front_record",
    "ideal_shock_search_window",
    "mesh_field_name",
    "morphology_raster_record",
    "reduce_t500_spatial_snapshot",
    "reduce_upstream_b_amplification_at_t500",
    "upstream_b_amplification_record",
    "upstream_b_window",
    "x1_cell_centers",
    "xy_cell_areas",
    "y_area_weighted_profile",
    "y_area_weighted_profile_record",
]
