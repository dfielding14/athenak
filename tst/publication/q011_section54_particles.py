"""Pure fail-closed particle reducers for the Q-011 Section 5.4 campaign."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from numbers import Real
from typing import Mapping

import numpy as np

if __package__:
    from . import analyze_q011_section54_outputs as arithmetic
    from . import q011_section54_model as frozen_model
else:
    import analyze_q011_section54_outputs as arithmetic
    import q011_section54_model as frozen_model


SCHEMA_VERSION = 1
GAS_GAMMA = 5.0 / 3.0
FROZEN_MODEL_SOURCE = "q011_section54_model"
UPSTREAM_SPEED_U0 = frozen_model.UPSTREAM_SPEED_U0
IDEAL_SURFACE_SPEED = frozen_model.IDEAL_SURFACE_SPEED
PARTICLE_LIGHT_SPEED = frozen_model.LIGHT_SPEED
BIRTH_TIME_MINIMUM = 45.0
SHOCK_INJECTED_SOURCE = 1
CHI_BIN_EDGES = tuple(2.0 ** (quarter_octave / 4.0) for quarter_octave in range(41))
LATE_SLOPE_SNAPSHOT_TIME = 1200.0
LATE_SLOPE_FIT_WINDOW = (20.0, 160.0)
LATE_SLOPE_MINIMUM_POSITIVE_BINS = 8
LATE_SLOPE_TARGET = -1.5
LATE_SLOPE_ABSOLUTE_TOLERANCE = 0.2
MAX_OVERFLOW_MACRO_WEIGHT_FRACTION = 0.001


class ParticleReducerError(ValueError):
    """Raised when decoded particle arrays cannot produce a bounded record."""


@dataclass(frozen=True)
class DownstreamFilterResult:
    """One downstream selection mask and its JSON-ready disjoint census."""

    mask: np.ndarray
    record: dict[str, object]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ParticleReducerError(message)


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


def _finite_vector(values: object, label: str, *, allow_empty: bool = True) -> np.ndarray:
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ParticleReducerError(f"{label} must be a numeric array") from error
    _require(array.ndim == 1, f"{label} shape drifted: expected (nparticle,)")
    _require(allow_empty or array.size > 0, f"{label} must not be empty")
    _require(np.all(np.isfinite(array)), f"{label} must be finite")
    return array


def _finite_triplets(values: object, label: str) -> np.ndarray:
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ParticleReducerError(f"{label} must be a numeric array") from error
    _require(
        array.ndim == 2 and array.shape[1:] == (3,),
        f"{label} shape drifted: expected (nparticle, 3)",
    )
    _require(np.all(np.isfinite(array)), f"{label} must be finite")
    return array


def _integer_vector(values: object, label: str) -> np.ndarray:
    array = np.asarray(values)
    _require(array.ndim == 1, f"{label} shape drifted: expected (nparticle,)")
    _require(array.dtype.kind in "iu", f"{label} must contain decoded integers")
    return array


def _macro_weights(values: object) -> np.ndarray:
    weights = _finite_vector(values, "macro_weight")
    _require(np.all(weights >= 0.0), "macro_weight must be non-negative")
    return weights


def _require_same_length(reference: np.ndarray, **arrays: np.ndarray) -> None:
    for label, array in arrays.items():
        _require(
            array.shape[0] == reference.shape[0],
            f"{label} shape drifted: particle counts disagree",
        )


def _weighted_census(mask: np.ndarray, weights: np.ndarray) -> dict[str, object]:
    return {
        "particle_count": int(np.count_nonzero(mask)),
        "macro_weight": float(np.sum(weights[mask])),
    }


def ideal_surface_x1(snapshot_time: object) -> float:
    """Return the frozen ideal injection-surface position at one snapshot time."""
    try:
        return frozen_model.x_ideal(snapshot_time)
    except frozen_model.ModelContractError as error:
        raise ParticleReducerError(f"snapshot_time is invalid: {error}") from error


def downstream_filter(
    *,
    snapshot_time: object,
    x1: object,
    cr_source: object,
    birth_time: object,
    macro_weight: object,
) -> DownstreamFilterResult:
    """Select downstream shock-injected particles and archive a disjoint census."""
    positions = _finite_vector(x1, "x1")
    sources = _integer_vector(cr_source, "cr_source")
    births = _finite_vector(birth_time, "birth_time")
    weights = _macro_weights(macro_weight)
    _require_same_length(
        positions,
        cr_source=sources,
        birth_time=births,
        macro_weight=weights,
    )
    surface = ideal_surface_x1(snapshot_time)

    rejected_wrong_source = sources != SHOCK_INJECTED_SOURCE
    source_admitted = ~rejected_wrong_source
    rejected_early_birth_time = source_admitted & (births < BIRTH_TIME_MINIMUM)
    provenance_admitted = source_admitted & ~rejected_early_birth_time
    rejected_upstream = provenance_admitted & (positions > surface)
    rejected_on_surface = provenance_admitted & (positions == surface)
    admitted_downstream = provenance_admitted & (positions < surface)

    partitions = np.stack(
        (
            rejected_wrong_source,
            rejected_early_birth_time,
            rejected_upstream,
            rejected_on_surface,
            admitted_downstream,
        ),
        axis=0,
    )
    _require(
        np.all(np.sum(partitions, axis=0) == 1),
        "internal downstream-filter census is not a disjoint complete partition",
    )
    rejected = ~admitted_downstream
    record = {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_downstream_particle_filter",
        "snapshot_time_omega0_inverse": float(snapshot_time),
        "ideal_surface_x1_c_over_omega_pi": surface,
        "selection": {
            "source_rule": "cr_source == 1",
            "birth_time_rule": "birth_time >= 45",
            "spatial_rule": "x1 < x_ideal(t)",
        },
        "disjoint_census": {
            "all_particles": _weighted_census(np.ones(positions.size, dtype=bool), weights),
            "admitted_downstream": _weighted_census(admitted_downstream, weights),
            "rejected_total": _weighted_census(rejected, weights),
            "rejected_wrong_source": _weighted_census(rejected_wrong_source, weights),
            "rejected_early_birth_time": _weighted_census(
                rejected_early_birth_time, weights
            ),
            "rejected_upstream": _weighted_census(rejected_upstream, weights),
            "rejected_on_surface": _weighted_census(rejected_on_surface, weights),
        },
    }
    return DownstreamFilterResult(mask=admitted_downstream, record=record)


def reconstruct_chi_from_pvtk_velocity(velocity: object) -> np.ndarray:
    """Reconstruct frozen chi from PVTK physical velocity output.

    AthenaK stores momentum per mass in paper mode but emits physical velocity.
    The inverse output transform is p/m = gamma(v) * v.
    """
    physical_velocity = _finite_triplets(velocity, "velocity")
    speed_squared = np.sum(physical_velocity * physical_velocity, axis=1)
    _require(np.all(np.isfinite(speed_squared)), "velocity magnitude must be finite")
    _require(
        np.all(speed_squared < PARTICLE_LIGHT_SPEED**2),
        "velocity magnitude must remain below the frozen particle light speed",
    )
    gamma_squared = 1.0 / (1.0 - speed_squared / PARTICLE_LIGHT_SPEED**2)
    chi = gamma_squared * speed_squared / (UPSTREAM_SPEED_U0**2)
    _require(np.all(np.isfinite(chi)), "reconstructed chi must be finite")
    _require(np.all(chi >= 0.0), "reconstructed chi must be non-negative")
    return chi


def weighted_spectrum_record(chi: object, macro_weight: object) -> dict[str, object]:
    """Return a fixed-bin weighted spectrum with bounded overflow accounting."""
    values = _finite_vector(chi, "chi")
    weights = _macro_weights(macro_weight)
    _require_same_length(values, macro_weight=weights)
    try:
        histogram = arithmetic.fixed_histogram(values, CHI_BIN_EDGES, weights=weights)
    except arithmetic.AnalysisError as error:
        raise ParticleReducerError(f"fixed chi histogram failed: {error}") from error

    widths = np.diff(histogram.bin_edges)
    f_chi = histogram.weighted_counts / widths
    total_post_filter_macro_weight = float(np.sum(weights))
    _require(
        math.isfinite(total_post_filter_macro_weight),
        "total post-filter macro weight must be finite",
    )
    accounted_weight = float(
        np.sum(histogram.weighted_counts)
        + histogram.underflow_weight
        + histogram.overflow_weight
    )
    _require(
        math.isclose(
            accounted_weight,
            total_post_filter_macro_weight,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        ),
        "fixed chi histogram macro-weight accounting drifted",
    )
    normalized = np.zeros_like(f_chi)
    overflow_fraction = 0.0
    if total_post_filter_macro_weight > 0.0:
        normalized = (
            histogram.geometric_bin_centers * f_chi / total_post_filter_macro_weight
        )
        overflow_fraction = histogram.overflow_weight / total_post_filter_macro_weight
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_fixed_weighted_chi_spectrum",
        "bin_edges": histogram.bin_edges.tolist(),
        "counts": histogram.counts.tolist(),
        "weighted_counts": histogram.weighted_counts.tolist(),
        "f_chi": f_chi.tolist(),
        "normalized_chi_f_chi": normalized.tolist(),
        "particle_count_post_filter": int(values.size),
        "particle_count_in_bins": int(np.sum(histogram.counts)),
        "macro_weight_in_bins": float(np.sum(histogram.weighted_counts)),
        "total_post_filter_macro_weight": total_post_filter_macro_weight,
        "underflow_count": histogram.underflow_count,
        "underflow_macro_weight": histogram.underflow_weight,
        "overflow_count": histogram.overflow_count,
        "overflow_macro_weight": histogram.overflow_weight,
        "overflow_macro_weight_fraction": overflow_fraction,
        "overflow_gate": {
            "maximum_macro_weight_fraction": MAX_OVERFLOW_MACRO_WEIGHT_FRACTION,
            "passed": bool(overflow_fraction <= MAX_OVERFLOW_MACRO_WEIGHT_FRACTION),
        },
        "normalization": (
            "chi * f_chi / total_post_filter_macro_weight, including fixed-bin "
            "underflow and overflow macro weight"
        ),
    }


def late_slope_record(f_chi: object) -> dict[str, object]:
    """Fit and gate the frozen late-time downstream energy-tail slope."""
    values = _finite_vector(f_chi, "f_chi", allow_empty=False)
    _require(
        values.shape == (len(CHI_BIN_EDGES) - 1,),
        "f_chi shape drifted: expected one value per frozen chi bin",
    )
    _require(np.all(values >= 0.0), "f_chi must be non-negative")
    try:
        fit = arithmetic.fit_fixed_bin_loglog_slope(
            CHI_BIN_EDGES,
            values,
            fit_window=LATE_SLOPE_FIT_WINDOW,
            minimum_bins=LATE_SLOPE_MINIMUM_POSITIVE_BINS,
        )
    except arithmetic.AnalysisError as error:
        raise ParticleReducerError(f"late chi slope fit failed: {error}") from error
    absolute_error = abs(fit.slope - LATE_SLOPE_TARGET)
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_late_chi_slope",
        "fit_window_chi": list(LATE_SLOPE_FIT_WINDOW),
        "minimum_positive_fit_bins": LATE_SLOPE_MINIMUM_POSITIVE_BINS,
        "selected_bin_indices": list(fit.selected_bin_indices),
        "positive_fit_bin_count": len(fit.selected_bin_indices),
        "slope": fit.slope,
        "intercept": fit.intercept,
        "target_slope": LATE_SLOPE_TARGET,
        "absolute_tolerance": LATE_SLOPE_ABSOLUTE_TOLERANCE,
        "absolute_error": absolute_error,
        "slope_gate_passed": bool(absolute_error <= LATE_SLOPE_ABSOLUTE_TOLERANCE),
    }


def reduce_particle_snapshot(
    *,
    snapshot_time: object,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
    evaluate_late_slope: bool = False,
) -> dict[str, object]:
    """Reduce one decoded particle snapshot into deterministic archive records."""
    _require(
        isinstance(evaluate_late_slope, bool),
        "evaluate_late_slope must be a boolean",
    )
    time = _finite_scalar(snapshot_time, "snapshot_time", minimum=0.0)
    decoded_points = _finite_triplets(points, "points")
    decoded_velocity = _finite_triplets(velocity, "velocity")
    sources = _integer_vector(cr_source, "cr_source")
    births = _finite_vector(birth_time, "birth_time")
    weights = _macro_weights(macro_weight)
    _require_same_length(
        decoded_points,
        velocity=decoded_velocity,
        cr_source=sources,
        birth_time=births,
        macro_weight=weights,
    )
    selected = downstream_filter(
        snapshot_time=time,
        x1=decoded_points[:, 0],
        cr_source=sources,
        birth_time=births,
        macro_weight=weights,
    )
    chi = reconstruct_chi_from_pvtk_velocity(decoded_velocity[selected.mask])
    spectrum = weighted_spectrum_record(chi, weights[selected.mask])
    record: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_particle_snapshot_reduction",
        "snapshot_time_omega0_inverse": time,
        "chi_reconstruction": {
            "input_vector": "pvtk physical velocity",
            "formula": "chi = gamma(v)^2 * |v|^2 / u0^2",
            "frozen_model_source": FROZEN_MODEL_SOURCE,
            "particle_light_speed": PARTICLE_LIGHT_SPEED,
            "u0_over_u_a0": UPSTREAM_SPEED_U0,
        },
        "particle_filter": selected.record,
        "weighted_spectrum": spectrum,
    }
    if evaluate_late_slope:
        _require(
            time == LATE_SLOPE_SNAPSHOT_TIME,
            "late slope may only be evaluated at the frozen t=1200 snapshot",
        )
        record["late_slope"] = late_slope_record(spectrum["f_chi"])
    return record


def canonical_record_bytes(record: Mapping[str, object]) -> bytes:
    """Encode one reducer record as deterministic archive JSON bytes."""
    _require(isinstance(record, Mapping), "archive record must be a mapping")
    try:
        payload = json.dumps(
            record,
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError) as error:
        raise ParticleReducerError("archive record is not finite JSON data") from error
    return (payload + "\n").encode("utf-8")


__all__ = [
    "BIRTH_TIME_MINIMUM",
    "CHI_BIN_EDGES",
    "DownstreamFilterResult",
    "FROZEN_MODEL_SOURCE",
    "IDEAL_SURFACE_SPEED",
    "LATE_SLOPE_ABSOLUTE_TOLERANCE",
    "LATE_SLOPE_FIT_WINDOW",
    "LATE_SLOPE_MINIMUM_POSITIVE_BINS",
    "LATE_SLOPE_SNAPSHOT_TIME",
    "LATE_SLOPE_TARGET",
    "MAX_OVERFLOW_MACRO_WEIGHT_FRACTION",
    "PARTICLE_LIGHT_SPEED",
    "ParticleReducerError",
    "SHOCK_INJECTED_SOURCE",
    "UPSTREAM_SPEED_U0",
    "canonical_record_bytes",
    "downstream_filter",
    "ideal_surface_x1",
    "late_slope_record",
    "reconstruct_chi_from_pvtk_velocity",
    "reduce_particle_snapshot",
    "weighted_spectrum_record",
]
