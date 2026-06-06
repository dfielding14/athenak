#!/usr/bin/env python3
"""Fail-closed source-local physical-applicability diagnostics for Q011.

This additive successor consumes already decoded Q011 production-science mesh,
deposited-moment, and particle products.  It computes snapshot applicability
diagnostics and validates a future runtime time/escape record.  It does not
launch work, mutate policy, authorize qualifying-output inspection, or close
any scientific claim.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import wraps
import math
from numbers import Real
from types import MappingProxyType
from typing import Any

import numpy as np

if __package__:
    from . import q011_section54_production_science_successor_v1 as science
else:
    import q011_section54_production_science_successor_v1 as science


SCHEMA_VERSION = 1
SUCCESSOR_ID = "q011_section54_physical_applicability_successor_v1"
SNAPSHOT_RECORD_TYPE = "q011_section54_physical_applicability_snapshot_v1"
HISTORY_RECORD_TYPE = "q011_section54_physical_applicability_history_v1"
RUNTIME_RECORD_TYPE = "q011_section54_runtime_time_escape_applicability_v1"
QUALIFICATION_EFFECT = (
    "source_local_non_authorizing_diagnostic_and_gate_only_no_launch_no_policy_"
    "mutation_no_qualifying_output_inspection_no_claim_closure"
)

AUTHORIZATION: Mapping[str, bool] = MappingProxyType(
    {
        "launch_authorized": False,
        "policy_mutation_authorized": False,
        "qualifying_output_inspection_authorized": False,
        "claim_closure_authorized": False,
    }
)

# These values and representation labels must be supplied exactly.  The
# successor intentionally does not infer normalization from deck arithmetic.
EXACT_NORMALIZATION: Mapping[str, object] = MappingProxyType(
    {
        "selected_species_count": 1,
        "selected_species_index": 0,
        "selected_species_charge_sign": 1,
        "selected_species_abs_q_over_m": 1.0,
        "reference_rho0": 1.0,
        "reference_b0": 1.0,
        "particle_light_speed": 10000.0,
        "gas_density_representation": "ion_mass_density_rho_g",
        "prtcl_rho_representation": "deposited_rho_CR_over_c_named_rho_q",
        "prtcl_j_representation": "deposited_J_CR_over_c_named_J_q",
        "alfven_speed_formula": "abs_B_over_sqrt_rho_g",
        "gas_frame_current_formula": "J_q_minus_rho_q_times_v_g",
        "particle_momentum_formula": "gamma_of_v_times_v",
        "particle_gyroradius_formula": "abs_p_over_abs_q_over_m_times_abs_B_local",
    }
)

# Bai et al. require R << 1, Lambda << 1 for the no-Hall approximation, and
# MHD-PIC scales much larger than d_i, but do not prescribe these numeric
# thresholds.  They are conservative AthenaK-selected applicability bounds.
R_MAXIMUM = 0.01
LAMBDA_MAXIMUM = 0.1
S_DELTA_MINIMUM = 1.0
LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM = 10.0
SUB_10DI_POWER_FRACTION_MAXIMUM = 0.05

# Sun & Bai state that the transverse domain should contain several
# high-energy gyroradii but do not prescribe these numeric thresholds.  They
# are conservative AthenaK-selected containment bounds.
RG_Q999_OVER_LY_MAXIMUM = 1.0 / 8.0
RG_MAXIMUM_OVER_LY_MAXIMUM = 1.0 / 2.0
RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM = 1.0e-3

STARTUP_REMOVAL_TIME = 45.0
EXPECTED_TERMINAL_TIME = 1200.0
PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES = 1000
SHOCK_TRANSITION_HALF_WIDTH = 120.0
DETECTED_FRONT_REGIONS: Mapping[str, tuple[float, float] | None] = MappingProxyType(
    {
        "full_domain": None,
        "detected_front_downstream": (-1200.0, -120.0),
        "detected_front_precursor": (120.0, 1200.0),
        "detected_front_far_upstream": (1200.0, 2400.0),
    }
)
DI_MAGNETIC_SPECTRUM_REGION = "detected_front_precursor"

CLAIM_REJECTION_RULES: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        "Q011-APP-NORM": (
            "all_physical_MHD_PIC_Bell_shock_and_DSA_claims",
        ),
        "Q011-APP-R": (
            "all_physical_MHD_PIC_Bell_shock_and_DSA_claims",
        ),
        "Q011-APP-LAMBDA": (
            "Hall_negligible_claim",
            "Bell_mechanism_claim",
            "physical_magnetic_amplification_claim",
            "physical_DSA_scattering_claim",
        ),
        "Q011-APP-DI": (
            "physical_precursor_turbulence_claim",
            "physical_scattering_interpretation",
        ),
        "Q011-APP-RG": (
            "Emax_claim",
            "high_energy_slope_or_cutoff_claim",
            "acceleration_rate_claim",
            "acceleration_efficiency_claim",
        ),
        "Q011-APP-TIME": (
            "all_history_and_global_applicability_claims",
        ),
    }
)
PERMANENT_CLAIM_EXCLUSIONS = (
    "microscopic_shock_structure_claim",
    "self_consistent_injection_claim",
)


class PhysicalApplicabilityError(ValueError):
    """Raised when Q011 physical applicability cannot be established."""


_UNDERLYING_EXCEPTIONS = (
    science.ProductionScienceError,
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
    def decorate(function):
        @wraps(function)
        def wrapped(*args: object, **kwargs: object):
            try:
                return function(*args, **kwargs)
            except PhysicalApplicabilityError:
                raise
            except _UNDERLYING_EXCEPTIONS as error:
                raise PhysicalApplicabilityError(f"{label} failed: {error}") from error

        return wrapped

    return decorate


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PhysicalApplicabilityError(message)


def _exact_keys(value: object, expected: set[str], label: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), f"{label} must be a mapping")
    _require(set(value) == expected, f"{label} keys drifted")
    return value


def _finite_scalar(
    value: object, label: str, *, minimum: float | None = None
) -> float:
    _require(
        isinstance(value, Real) and not isinstance(value, (bool, np.bool_)),
        f"{label} must be a real scalar",
    )
    result = float(value)
    _require(math.isfinite(result), f"{label} must be finite")
    if minimum is not None:
        _require(result >= minimum, f"{label} must be at least {minimum}")
    return result


def _nonnegative_int(value: object, label: str) -> int:
    _require(type(value) is int and value >= 0, f"{label} must be a non-negative integer")
    return value


def _strict_bool(value: object, label: str) -> bool:
    _require(type(value) is bool, f"{label} must be a boolean")
    return value


def _exact_numeric(value: object, expected: float, label: str) -> float:
    decoded = _finite_scalar(value, label)
    _require(decoded == expected, f"{label} drifted")
    return decoded


def _validate_authorization(value: object, label: str) -> dict[str, bool]:
    authorization = _exact_keys(value, set(AUTHORIZATION), label)
    for key, expected in AUTHORIZATION.items():
        _require(
            _strict_bool(authorization[key], f"{label} {key}") is expected,
            f"{label} {key} drifted",
        )
    return dict(AUTHORIZATION)


def _finite_array(values: object, label: str, *, ndim: int | None = None) -> np.ndarray:
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise PhysicalApplicabilityError(f"{label} must be numeric") from error
    if ndim is not None:
        _require(array.ndim == ndim, f"{label} must be {ndim}-dimensional")
    _require(array.size > 0, f"{label} must not be empty")
    _require(np.all(np.isfinite(array)), f"{label} must be finite")
    return array


def _readonly(values: np.ndarray) -> np.ndarray:
    result = np.array(values, copy=True)
    result.setflags(write=False)
    return result


def _weighted_quantile(
    values: np.ndarray, weights: np.ndarray, quantile: float, *, label: str
) -> float:
    samples = _finite_array(values, f"{label} values").reshape(-1)
    sample_weights = _finite_array(weights, f"{label} weights").reshape(-1)
    _require(samples.size == sample_weights.size, f"{label} sizes disagree")
    _require(np.all(sample_weights >= 0.0), f"{label} weights must be non-negative")
    positive = sample_weights > 0.0
    _require(np.any(positive), f"{label} requires positive total weight")
    order = np.argsort(samples[positive], kind="stable")
    ordered = samples[positive][order]
    ordered_weights = sample_weights[positive][order]
    cumulative = np.cumsum(ordered_weights)
    threshold = quantile * float(cumulative[-1])
    index = int(np.searchsorted(cumulative, threshold, side="left"))
    return float(ordered[min(index, ordered.size - 1)])


def _weighted_statistics(
    values: np.ndarray,
    area_weights: np.ndarray,
    current_weights: np.ndarray,
    *,
    label: str,
) -> dict[str, Any]:
    samples = _finite_array(values, f"{label} values").reshape(-1)
    area = _finite_array(area_weights, f"{label} area weights").reshape(-1)
    current = _finite_array(current_weights, f"{label} current weights").reshape(-1)
    _require(samples.size == area.size == current.size, f"{label} sizes disagree")
    _require(np.all(area > 0.0), f"{label} area weights must be positive")
    _require(np.all(current >= 0.0), f"{label} current weights must be non-negative")

    def summarize(weights: np.ndarray, weight_label: str) -> dict[str, Any]:
        total = float(np.sum(weights))
        if total == 0.0:
            return {
                "available": False,
                "reason": f"zero_total_{weight_label}_weight",
                "total_weight": 0.0,
                "weighted_mean": None,
                "weighted_quantiles": None,
            }
        return {
            "available": True,
            "reason": None,
            "total_weight": total,
            "weighted_mean": float(np.sum(samples * weights) / total),
            "weighted_quantiles": {
                "q500": _weighted_quantile(samples, weights, 0.5, label=label),
                "q900": _weighted_quantile(samples, weights, 0.9, label=label),
                "q990": _weighted_quantile(samples, weights, 0.99, label=label),
                "q999": _weighted_quantile(samples, weights, 0.999, label=label),
            },
        }

    return {
        "cell_count": int(samples.size),
        "local_minimum": float(np.min(samples)),
        "local_maximum": float(np.max(samples)),
        "area_weighted": summarize(area, "area"),
        "gas_frame_current_weighted": summarize(current, "gas_frame_current"),
    }


def _validate_exact_normalization(value: object) -> dict[str, object]:
    normalization = _exact_keys(
        value, set(EXACT_NORMALIZATION), "Q011 exact normalization"
    )
    for key, expected in EXACT_NORMALIZATION.items():
        _require(
            type(normalization[key]) is type(expected) and normalization[key] == expected,
            f"Q011 exact normalization field {key!r} drifted",
        )
    return dict(EXACT_NORMALIZATION)


def _uniform_spacing(faces: np.ndarray, label: str) -> float:
    widths = np.diff(_finite_array(faces, f"{label} faces", ndim=1))
    _require(np.all(widths > 0.0), f"{label} faces must increase")
    reference = float(widths[0])
    _require(
        np.allclose(widths, reference, rtol=0.0, atol=1.0e-12 * max(1.0, reference)),
        f"{label} composite spacing must be uniform",
    )
    return reference


def _cell_areas(state: science.ComposedMHDState) -> np.ndarray:
    return np.diff(state.x2_faces)[:, None] * np.diff(state.x1_faces)[None, :]


def _x1_centers(state: science.ComposedMHDState) -> np.ndarray:
    return 0.5 * (state.x1_faces[:-1] + state.x1_faces[1:])


def _detected_front(
    mhd_dataset: object,
    *,
    nominal_slot_time: object,
    observed_committed_time: object,
    target_level: int | None,
) -> tuple[science.ComposedMHDState, float]:
    state = science.compose_full_mhd_state(
        mhd_dataset,
        nominal_slot_time=nominal_slot_time,
        observed_committed_time=observed_committed_time,
        target_level=target_level,
    )
    reduced = science.reduce_mhd_snapshot(
        mhd_dataset,
        nominal_slot_time=nominal_slot_time,
        observed_committed_time=observed_committed_time,
        target_level=target_level,
    )
    return state, float(reduced["detected_front"]["x_front_c_over_omega_pi"])


def _region_masks(
    state: science.ComposedMHDState, front_x: float
) -> dict[str, np.ndarray]:
    x = _x1_centers(state)
    shape = state.fields_y_x["dens"].shape
    masks: dict[str, np.ndarray] = {}
    for name, offsets in DETECTED_FRONT_REGIONS.items():
        if offsets is None:
            masks[name] = np.ones(shape, dtype=bool)
            continue
        lower = front_x + offsets[0]
        upper = front_x + offsets[1]
        _require(
            lower >= state.x1_faces[0] and upper <= state.x1_faces[-1],
            f"{name} escaped the retained x1 domain",
        )
        selected_x = (x > lower) & (x < upper)
        _require(np.count_nonzero(selected_x) >= 4, f"{name} requires at least four x1 cells")
        masks[name] = np.broadcast_to(selected_x[None, :], shape)
    return masks


def _actual_leaf_spacing_maps(
    state: science.ComposedMHDState,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dx1_target = _uniform_spacing(state.x1_faces, "x1")
    dx2_target = _uniform_spacing(state.x2_faces, "x2")
    level_factor = np.power(
        2.0, state.target_level - state.source_levels_y_x.astype(np.float64)
    )
    _require(np.all(level_factor >= 1.0), "source level exceeds target composite level")
    dx1_leaf = dx1_target * level_factor
    dx2_leaf = dx2_target * level_factor
    delta_leaf = np.minimum(dx1_leaf, dx2_leaf)
    return dx1_leaf, dx2_leaf, delta_leaf


def _magnetic_spectrum(
    state: science.ComposedMHDState,
    mask: np.ndarray,
    *,
    di_max: float,
) -> dict[str, Any]:
    selected_x = np.any(mask, axis=0)
    selected_y = np.any(mask, axis=1)
    _require(
        np.all(mask == np.outer(selected_y, selected_x)),
        "magnetic spectrum region must be rectangular",
    )
    nx = int(np.count_nonzero(selected_x))
    ny = int(np.count_nonzero(selected_y))
    _require(nx >= 4 and ny >= 4, "magnetic spectrum requires at least four cells per axis")
    dx1 = _uniform_spacing(state.x1_faces, "x1")
    dx2 = _uniform_spacing(state.x2_faces, "x2")
    magnetic = np.stack(
        [
            state.fields_y_x[field][np.ix_(selected_y, selected_x)]
            for field in ("bcc1", "bcc2", "bcc3")
        ],
        axis=0,
    )
    fluctuation = magnetic - np.mean(magnetic, axis=(1, 2), keepdims=True)
    window = np.outer(np.hanning(ny), np.hanning(nx))
    transformed = np.fft.fftn(fluctuation * window[None, :, :], axes=(-2, -1))
    power = np.sum(np.abs(transformed) ** 2, axis=0)
    power[0, 0] = 0.0
    total_power = float(np.sum(power))
    _require(
        math.isfinite(total_power) and total_power > np.finfo(np.float64).tiny,
        "precursor magnetic spectrum has no resolvable fluctuation power",
    )
    k1 = 2.0 * math.pi * np.fft.fftfreq(nx, d=dx1)
    k2 = 2.0 * math.pi * np.fft.fftfreq(ny, d=dx2)
    k1_grid, k2_grid = np.meshgrid(k1, k2, indexing="xy")
    kmag = np.sqrt(k1_grid * k1_grid + k2_grid * k2_grid)
    nonzero = kmag > 0.0
    mean_k = float(np.sum(kmag[nonzero] * power[nonzero]) / total_power)
    _require(mean_k > 0.0 and math.isfinite(mean_k), "characteristic magnetic k is invalid")
    lambda_char = 2.0 * math.pi / mean_k
    sub_10di = kmag > 2.0 * math.pi / (10.0 * di_max)
    sub_fraction = float(np.sum(power[sub_10di]) / total_power)
    return {
        "region": DI_MAGNETIC_SPECTRUM_REGION,
        "method": (
            "two_dimensional_Hann_windowed_vector_delta_B_power_weighted_mean_k"
        ),
        "target_composite_dx1": dx1,
        "target_composite_dx2": dx2,
        "selected_nx1": nx,
        "selected_nx2": ny,
        "total_windowed_delta_B_power": total_power,
        "local_di_maximum": di_max,
        "lambda_B_characteristic": lambda_char,
        "lambda_B_characteristic_over_local_di_maximum": lambda_char / di_max,
        "sub_10di_magnetic_power_fraction": sub_fraction,
        "shock_transition_excluded": True,
    }


def _tsc_axis_indices_weights(
    coordinates: np.ndarray,
    faces: np.ndarray,
    *,
    periodic: bool,
    label: str,
) -> tuple[np.ndarray, np.ndarray]:
    dx = _uniform_spacing(faces, label)
    count = faces.size - 1
    normalized = (coordinates - faces[0]) / dx - 0.5
    center = np.floor(normalized + 0.5).astype(np.int64)
    indices = center[:, None] + np.asarray([-1, 0, 1], dtype=np.int64)[None, :]
    distance = np.abs(normalized[:, None] - indices.astype(np.float64))
    weights = np.where(
        distance < 0.5,
        0.75 - distance * distance,
        np.where(distance < 1.5, 0.5 * (1.5 - distance) ** 2, 0.0),
    )
    _require(
        np.allclose(np.sum(weights, axis=1), 1.0, rtol=0.0, atol=2.0e-15),
        f"{label} TSC weights failed unity closure",
    )
    if periodic:
        indices %= count
    else:
        _require(
            np.all((indices >= 0) & (indices < count)),
            f"{label} particle TSC stencil escaped the retained nonperiodic domain",
        )
    return indices, weights


def _tsc_sample_magnetic_field(
    state: science.ComposedMHDState, points: np.ndarray
) -> np.ndarray:
    _require(points.ndim == 2 and points.shape[1] == 3, "particle points shape drifted")
    _require(
        np.all(points[:, 0] >= state.x1_faces[0])
        and np.all(points[:, 0] <= state.x1_faces[-1])
        and np.all(points[:, 1] >= state.x2_faces[0])
        and np.all(points[:, 1] <= state.x2_faces[-1]),
        "particle escaped retained x1-x2 domain",
    )
    x_indices, x_weights = _tsc_axis_indices_weights(
        points[:, 0], state.x1_faces, periodic=False, label="x1"
    )
    y_indices, y_weights = _tsc_axis_indices_weights(
        points[:, 1], state.x2_faces, periodic=True, label="x2"
    )
    sampled = np.zeros((points.shape[0], 3), dtype=np.float64)
    for component, field in enumerate(("bcc1", "bcc2", "bcc3")):
        values = state.fields_y_x[field]
        for iy in range(3):
            for ix in range(3):
                sampled[:, component] += (
                    y_weights[:, iy]
                    * x_weights[:, ix]
                    * values[y_indices[:, iy], x_indices[:, ix]]
                )
    _require(np.all(np.isfinite(sampled)), "TSC sampled magnetic field must be finite")
    return sampled


def _tsc_sample_scalar_field(
    state: science.ComposedMHDState,
    values_y_x: np.ndarray,
    points: np.ndarray,
    *,
    label: str,
) -> np.ndarray:
    _require(
        values_y_x.shape == state.fields_y_x["dens"].shape,
        f"{label} scalar field shape drifted",
    )
    x_indices, x_weights = _tsc_axis_indices_weights(
        points[:, 0], state.x1_faces, periodic=False, label="x1"
    )
    y_indices, y_weights = _tsc_axis_indices_weights(
        points[:, 1], state.x2_faces, periodic=True, label="x2"
    )
    sampled = np.zeros(points.shape[0], dtype=np.float64)
    for iy in range(3):
        for ix in range(3):
            sampled += (
                y_weights[:, iy]
                * x_weights[:, ix]
                * values_y_x[y_indices[:, iy], x_indices[:, ix]]
            )
    _require(np.all(np.isfinite(sampled)), f"TSC sampled {label} must be finite")
    return sampled


def _particle_scalar_exposure_statistics(
    values: np.ndarray,
    macro_weights: np.ndarray,
    energy_weights: np.ndarray,
    *,
    threshold: float,
    label: str,
) -> dict[str, Any]:
    samples = _finite_array(values, f"{label} samples", ndim=1)
    macro = _finite_array(macro_weights, f"{label} macro weights", ndim=1)
    energy = _finite_array(energy_weights, f"{label} energy weights", ndim=1)
    _require(samples.size == macro.size == energy.size, f"{label} sizes disagree")
    _require(
        np.all(macro >= 0.0) and np.all(energy >= 0.0),
        f"{label} particle weights must be non-negative",
    )

    def weighted(weights: np.ndarray, weight_label: str) -> dict[str, Any]:
        total = float(np.sum(weights))
        _require(total > 0.0, f"{label} {weight_label} total weight must be positive")
        return {
            "weighted_mean": float(np.sum(samples * weights) / total),
            "weighted_quantiles": {
                "q500": _weighted_quantile(samples, weights, 0.5, label=label),
                "q900": _weighted_quantile(samples, weights, 0.9, label=label),
                "q990": _weighted_quantile(samples, weights, 0.99, label=label),
                "q999": _weighted_quantile(samples, weights, 0.999, label=label),
            },
            "exceedance_fraction": float(np.sum(weights[samples > threshold]) / total),
        }

    return {
        "particle_count": int(samples.size),
        "local_minimum": float(np.min(samples)),
        "local_maximum": float(np.max(samples)),
        "threshold": threshold,
        "macro_weighted": weighted(macro, "macro"),
        "CR_kinetic_energy_weighted": weighted(energy, "CR kinetic energy"),
    }


def _particle_R_Lambda_exposure(
    state: science.ComposedMHDState,
    front_x: float,
    r_map: np.ndarray,
    lambda_map: np.ndarray,
    *,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
) -> dict[str, Any]:
    arrays = science._decoded_particle_arrays(
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )
    energetic = arrays["energetic"]
    _require(np.any(energetic), "particle exposure requires positive-weight active CRs")
    selected_points = arrays["points"][energetic]
    macro = arrays["weights"][energetic]
    energy = macro * arrays["specific_kinetic_energy"][energetic]
    _require(float(np.sum(energy)) > 0.0, "particle exposure energy weight must be positive")
    sampled_r = _tsc_sample_scalar_field(
        state, r_map, selected_points, label="particle R exposure"
    )
    sampled_lambda = _tsc_sample_scalar_field(
        state, lambda_map, selected_points, label="particle Lambda exposure"
    )
    high_energy_threshold = _weighted_quantile(
        arrays["specific_kinetic_energy"][energetic],
        macro,
        0.99,
        label="particle high-energy-tail threshold",
    )
    populations = {
        "all_active": np.ones(sampled_r.size, dtype=bool),
        "detected_front_upstream": (
            selected_points[:, 0] > front_x + SHOCK_TRANSITION_HALF_WIDTH
        ),
        "high_energy_tail": (
            arrays["specific_kinetic_energy"][energetic] >= high_energy_threshold
        ),
    }
    result: dict[str, Any] = {}
    for name, selected in populations.items():
        if not np.any(selected):
            result[name] = {
                "available": False,
                "reason": "population_contains_no_positive_weight_active_particles",
                "particle_count": 0,
                "R": None,
                "Lambda": None,
            }
            continue
        result[name] = {
            "available": True,
            "reason": None,
            "particle_count": int(np.count_nonzero(selected)),
            "R": _particle_scalar_exposure_statistics(
                sampled_r[selected],
                macro[selected],
                energy[selected],
                threshold=R_MAXIMUM,
                label=f"{name} particle R exposure",
            ),
            "Lambda": _particle_scalar_exposure_statistics(
                sampled_lambda[selected],
                macro[selected],
                energy[selected],
                threshold=LAMBDA_MAXIMUM,
                label=f"{name} particle Lambda exposure",
            ),
        }
    return {
        "sampling": (
            "TSC_on_matched_finest_composite_periodic_x2_nonperiodic_x1_"
            "stencil_must_remain_retained"
        ),
        "population_definitions": {
            "all_active": "cr_source_eq_1_birth_time_ge_45_positive_macro_weight",
            "detected_front_upstream": "all_active_and_x1_gt_detected_front_plus_120",
            "high_energy_tail": (
                "all_active_with_specific_kinetic_energy_ge_macro_weighted_q990"
            ),
        },
        "high_energy_tail_specific_kinetic_energy_threshold": high_energy_threshold,
        "populations": result,
    }


def _particle_gyroradius(
    state: science.ComposedMHDState,
    *,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
) -> dict[str, Any]:
    arrays = science._decoded_particle_arrays(
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )
    selected = arrays["energetic"]
    count = int(np.count_nonzero(selected))
    _require(
        count >= PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES,
        "particle gyroradius q999 requires at least 1000 positive-weight active CRs",
    )
    selected_points = arrays["points"][selected]
    sampled_b = _tsc_sample_magnetic_field(state, selected_points)
    bmag = np.linalg.norm(sampled_b, axis=1)
    _require(np.all(bmag > 0.0), "particle local magnetic field must be nonzero")
    speed_squared = np.sum(arrays["velocity"][selected] ** 2, axis=1)
    gamma = 1.0 / np.sqrt(1.0 - speed_squared / science.PARTICLE_LIGHT_SPEED**2)
    momentum_magnitude = gamma * np.sqrt(speed_squared)
    rg = momentum_magnitude / bmag
    _require(np.all(np.isfinite(rg)) and np.all(rg >= 0.0), "particle gyroradius is invalid")
    macro = arrays["weights"][selected]
    energy = macro * arrays["specific_kinetic_energy"][selected]
    _require(float(np.sum(energy)) > 0.0, "particle energy weights must be positive")
    ly = float(state.x2_faces[-1] - state.x2_faces[0])
    macro_q999 = _weighted_quantile(rg, macro, 0.999, label="particle gyroradius macro")
    energy_q999 = _weighted_quantile(rg, energy, 0.999, label="particle gyroradius energy")
    maximum = float(np.max(rg))
    energy_fraction = float(np.sum(energy[rg > ly / 4.0]) / np.sum(energy))
    return {
        "population": "cr_source_eq_1_birth_time_ge_45_positive_macro_weight",
        "local_B_sampling": (
            "TSC_on_matched_finest_composite_periodic_x2_nonperiodic_x1_"
            "stencil_must_remain_retained"
        ),
        "gyroradius_formula": "abs_gamma_v_over_abs_q_over_m_times_abs_B_local",
        "selected_particle_count": count,
        "transverse_domain_size_Ly": ly,
        "macro_weighted": {
            "q500": _weighted_quantile(rg, macro, 0.5, label="particle gyroradius macro"),
            "q900": _weighted_quantile(rg, macro, 0.9, label="particle gyroradius macro"),
            "q990": _weighted_quantile(rg, macro, 0.99, label="particle gyroradius macro"),
            "q999": macro_q999,
        },
        "energy_weighted": {
            "q500": _weighted_quantile(rg, energy, 0.5, label="particle gyroradius energy"),
            "q900": _weighted_quantile(rg, energy, 0.9, label="particle gyroradius energy"),
            "q990": _weighted_quantile(rg, energy, 0.99, label="particle gyroradius energy"),
            "q999": energy_q999,
        },
        "maximum": maximum,
        "macro_q999_over_Ly": macro_q999 / ly,
        "energy_q999_over_Ly": energy_q999 / ly,
        "maximum_over_Ly": maximum / ly,
        "energy_fraction_with_rg_above_Ly_over_4": energy_fraction,
    }


def _claim_rejections(gates: Mapping[str, Mapping[str, Any]], maximum_lambda: float) -> list[str]:
    rejected = set(PERMANENT_CLAIM_EXCLUSIONS)
    for gate, claims in CLAIM_REJECTION_RULES.items():
        if not bool(gates[gate]["pass"]):
            rejected.update(claims)
    if maximum_lambda >= 1.0:
        rejected.add("no_Hall_run_approximates_target_plasma_claim")
    return sorted(rejected)


@dataclass(frozen=True)
class ApplicabilitySnapshot:
    """One JSON-ready summary plus immutable per-cell applicability maps."""

    record: Mapping[str, Any]
    cell_maps: Mapping[str, np.ndarray]


@_public_contract("Q011 physical-applicability snapshot reduction")
def reduce_physical_applicability_snapshot(
    mhd_dataset: object,
    current_datasets: Mapping[str, object],
    *,
    normalization: object,
    nominal_slot_time: object,
    observed_committed_time: object,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
    target_level: int | None = None,
) -> ApplicabilitySnapshot:
    """Reduce one matched snapshot into fail-closed physical-applicability gates."""
    normalized = _validate_exact_normalization(normalization)
    state, front_x = _detected_front(
        mhd_dataset,
        nominal_slot_time=nominal_slot_time,
        observed_committed_time=observed_committed_time,
        target_level=target_level,
    )
    currents = science._compose_matched_current_fields(mhd_dataset, state, current_datasets)
    rho_g = state.fields_y_x["dens"]
    rho_q = currents["prtcl_rho"]
    gas_velocity = np.stack(
        [state.fields_y_x[f"vel{component}"] for component in ("x", "y", "z")],
        axis=0,
    )
    j_lab = np.stack(
        [currents[f"prtcl_j{component}"] for component in ("x", "y", "z")],
        axis=0,
    )
    j_gas = j_lab - rho_q[None, :, :] * gas_velocity
    j_gas_magnitude = np.linalg.norm(j_gas, axis=0)
    bmag = np.sqrt(
        sum(state.fields_y_x[field] ** 2 for field in ("bcc1", "bcc2", "bcc3"))
    )
    _require(np.all(bmag > 0.0), "cell magnetic-field magnitude must be positive")
    va = bmag / np.sqrt(rho_g)
    denominator = rho_g + rho_q
    _require(np.all(denominator > 0.0), "rho_g plus rho_q must be positive")
    r_map = rho_q / denominator
    lambda_map = j_gas_magnitude / (denominator * va)
    di_map = 1.0 / np.sqrt(rho_g)
    dx1_leaf, dx2_leaf, delta_leaf = _actual_leaf_spacing_maps(state)
    s_delta_map = delta_leaf / di_map
    for name, values in {
        "R": r_map,
        "Lambda": lambda_map,
        "d_i": di_map,
        "S_delta": s_delta_map,
    }.items():
        _require(
            np.all(np.isfinite(values)) and np.all(values >= 0.0),
            f"{name} cell map must be finite and non-negative",
        )

    areas = _cell_areas(state)
    masks = _region_masks(state, front_x)
    region_statistics: dict[str, Any] = {}
    for region, mask in masks.items():
        area_weights = areas[mask]
        current_weights = area_weights * j_gas_magnitude[mask]
        offsets = DETECTED_FRONT_REGIONS[region]
        region_statistics[region] = {
            "detected_front_offsets_c_over_omega_pi": (
                None if offsets is None else list(offsets)
            ),
            "R": _weighted_statistics(
                r_map[mask], area_weights, current_weights, label=f"{region} R"
            ),
            "Lambda": _weighted_statistics(
                lambda_map[mask],
                area_weights,
                current_weights,
                label=f"{region} Lambda",
            ),
        }

    shock_transition = np.abs(_x1_centers(state) - front_x) <= SHOCK_TRANSITION_HALF_WIDTH
    di_mask = np.broadcast_to((~shock_transition)[None, :], rho_g.shape)
    _require(np.any(di_mask), "shock-transition exclusion removed every cell")
    s_delta_min = float(np.min(s_delta_map[di_mask]))
    precursor_mask = masks[DI_MAGNETIC_SPECTRUM_REGION]
    precursor_di_max = float(np.max(di_map[precursor_mask]))
    spectrum = _magnetic_spectrum(state, precursor_mask, di_max=precursor_di_max)
    particle = _particle_gyroradius(
        state,
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )
    particle_exposure = _particle_R_Lambda_exposure(
        state,
        front_x,
        r_map,
        lambda_map,
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )

    maximum_r = float(np.max(r_map))
    maximum_lambda = float(np.max(lambda_map))
    gates: dict[str, dict[str, Any]] = {
        "Q011-APP-NORM": {
            "pass": True,
            "threshold_provenance": "exact_equation_and_normalization_identity",
            "observed": "exact_match",
            "required": "exact_match",
        },
        "Q011-APP-R": {
            "pass": maximum_r <= R_MAXIMUM,
            "threshold_provenance": "AthenaK_selected_Bai_states_R_much_less_than_one",
            "observed_maximum": maximum_r,
            "required_maximum": R_MAXIMUM,
        },
        "Q011-APP-LAMBDA": {
            "pass": maximum_lambda <= LAMBDA_MAXIMUM,
            "threshold_provenance": (
                "AthenaK_selected_Bai_states_Lambda_much_less_than_one_for_no_Hall"
            ),
            "observed_maximum": maximum_lambda,
            "required_maximum": LAMBDA_MAXIMUM,
        },
        "Q011-APP-DI": {
            "pass": (
                s_delta_min >= S_DELTA_MINIMUM
                and spectrum["lambda_B_characteristic_over_local_di_maximum"]
                >= LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM
                and spectrum["sub_10di_magnetic_power_fraction"]
                <= SUB_10DI_POWER_FRACTION_MAXIMUM
            ),
            "threshold_provenance": (
                "AthenaK_selected_decade_scale_separation_Bai_is_qualitative"
            ),
            "observed_S_delta_minimum_excluding_shock_transition": s_delta_min,
            "required_S_delta_minimum": S_DELTA_MINIMUM,
            "observed_lambda_B_characteristic_over_local_di_maximum": spectrum[
                "lambda_B_characteristic_over_local_di_maximum"
            ],
            "required_lambda_B_characteristic_over_local_di_maximum": (
                LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM
            ),
            "observed_sub_10di_magnetic_power_fraction": spectrum[
                "sub_10di_magnetic_power_fraction"
            ],
            "required_sub_10di_magnetic_power_fraction_maximum": (
                SUB_10DI_POWER_FRACTION_MAXIMUM
            ),
        },
        "Q011-APP-RG": {
            "pass": (
                particle["macro_q999_over_Ly"] <= RG_Q999_OVER_LY_MAXIMUM
                and particle["energy_q999_over_Ly"] <= RG_Q999_OVER_LY_MAXIMUM
                and particle["maximum_over_Ly"] <= RG_MAXIMUM_OVER_LY_MAXIMUM
                and particle["energy_fraction_with_rg_above_Ly_over_4"]
                <= RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM
            ),
            "threshold_provenance": (
                "AthenaK_selected_interpretation_of_Sun_and_Bai_several_gyroradii"
            ),
            "observed_macro_q999_over_Ly": particle["macro_q999_over_Ly"],
            "observed_energy_q999_over_Ly": particle["energy_q999_over_Ly"],
            "required_q999_over_Ly_maximum": RG_Q999_OVER_LY_MAXIMUM,
            "observed_maximum_over_Ly": particle["maximum_over_Ly"],
            "required_maximum_over_Ly": RG_MAXIMUM_OVER_LY_MAXIMUM,
            "observed_energy_fraction_with_rg_above_Ly_over_4": particle[
                "energy_fraction_with_rg_above_Ly_over_4"
            ],
            "required_energy_fraction_with_rg_above_Ly_over_4_maximum": (
                RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM
            ),
        },
        "Q011-APP-TIME": {
            "pass": False,
            "threshold_provenance": "engineering_and_analysis_completeness",
            "observed": "not_available_from_snapshot",
            "required": (
                "complete_per_cycle_post_startup_runtime_extrema_particle_exposure_"
                "and_boundary_escape_ledger_through_t1200"
            ),
        },
    }
    record: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "record_type": SNAPSHOT_RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "authorization": dict(AUTHORIZATION),
        "nominal_slot_time": state.nominal_slot_time,
        "observed_committed_time": state.observed_committed_time,
        "detected_front_x1_c_over_omega_pi": front_x,
        "exact_normalization": normalized,
        "formulae": {
            "R": "rho_q_over_rho_g_plus_rho_q",
            "v_A": "abs_B_over_sqrt_rho_g",
            "Lambda": "abs_J_q_minus_rho_q_v_g_over_rho_g_plus_rho_q_times_v_A",
            "d_i": "rho_g_to_the_minus_one_half",
            "S_delta": "minimum_actual_leaf_dx1_dx2_over_local_d_i",
        },
        "cell_map_contract": {
            "returned_separately_as_immutable_numpy_arrays": True,
            "retention_required_for_physical_applicability_evidence": True,
            "shape_y_x": list(rho_g.shape),
            "target_composite_level": state.target_level,
            "x1_faces_c_over_omega_pi": state.x1_faces.tolist(),
            "x2_faces_c_over_omega_pi": state.x2_faces.tolist(),
            "maps": [
                "R",
                "Lambda",
                "d_i",
                "S_delta",
                "actual_leaf_dx1",
                "actual_leaf_dx2",
                "gas_frame_current_magnitude",
            ],
            "per_cell_maxima_required_no_rare_cell_waiver": True,
        },
        "regional_statistics": region_statistics,
        "ion_scale_separation": {
            "actual_leaf_spacing_used": True,
            "shock_transition_exclusion": {
                "center": "detected_front",
                "half_width_c_over_omega_pi": SHOCK_TRANSITION_HALF_WIDTH,
                "microscopic_shock_structure_claim_permanently_excluded": True,
                "self_consistent_injection_claim_permanently_excluded": True,
            },
            "S_delta_minimum_excluding_shock_transition": s_delta_min,
            "precursor_magnetic_spectrum": spectrum,
        },
        "particle_R_Lambda_exposure": particle_exposure,
        "particle_gyroradius_containment": particle,
        "gates": gates,
        "snapshot_gate_pass_excluding_time_completeness": all(
            gates[name]["pass"]
            for name in ("Q011-APP-NORM", "Q011-APP-R", "Q011-APP-LAMBDA", "Q011-APP-DI", "Q011-APP-RG")
        ),
        "claim_rejection_rules": {
            gate: list(claims) for gate, claims in CLAIM_REJECTION_RULES.items()
        },
        "permanent_claim_exclusions": list(PERMANENT_CLAIM_EXCLUSIONS),
        "claim_rejections": _claim_rejections(gates, maximum_lambda),
        "runtime_time_escape_evidence_required": True,
    }
    maps = MappingProxyType(
        {
            "R": _readonly(r_map),
            "Lambda": _readonly(lambda_map),
            "d_i": _readonly(di_map),
            "S_delta": _readonly(s_delta_map),
            "actual_leaf_dx1": _readonly(dx1_leaf),
            "actual_leaf_dx2": _readonly(dx2_leaf),
            "gas_frame_current_magnitude": _readonly(j_gas_magnitude),
        }
    )
    return ApplicabilitySnapshot(record=MappingProxyType(record), cell_maps=maps)


_COVERAGE_KEYS = {
    "complete",
    "sampling_mode",
    "post_startup_removal_start_time",
    "terminal_time",
    "first_cycle",
    "last_cycle",
    "covered_cycle_count",
    "expected_contiguous_cycle_count",
    "gap_count",
    "restart_segment_count",
    "restart_segments_complete",
}
_EXTREMA_KEYS = {
    "R_maximum",
    "Lambda_maximum",
    "S_delta_minimum_excluding_shock_transition",
    "lambda_B_characteristic_over_local_di_maximum_minimum",
    "sub_10di_magnetic_power_fraction_maximum",
    "macro_q999_over_Ly_maximum",
    "energy_q999_over_Ly_maximum",
    "particle_rg_maximum_over_Ly",
    "energy_fraction_with_rg_above_Ly_over_4_maximum",
}
_PARTICLE_EXPOSURE_KEYS = {
    "complete",
    "method",
    "active_particle_updates_included",
    "pre_destruction_boundary_events_included",
    "escaped_particles_included",
    "observation_count",
    "maximum_sampled_R",
    "maximum_sampled_Lambda",
    "cumulative_macro_weighted_R_exceedance_fraction",
    "cumulative_CR_energy_weighted_R_exceedance_fraction",
    "cumulative_macro_weighted_Lambda_exceedance_fraction",
    "cumulative_CR_energy_weighted_Lambda_exceedance_fraction",
}
_FACE_KEYS = {
    "escaped_particle_count",
    "escaped_macro_weight",
    "escaped_kinetic_energy",
}
_ESCAPE_KEYS = {
    "complete",
    "scope",
    "nonperiodic_faces",
    "periodic_faces",
    "total_escaped_particle_count",
    "total_escaped_macro_weight",
    "total_escaped_kinetic_energy",
    "unaccounted_particle_count",
    "unaccounted_macro_weight",
    "unaccounted_kinetic_energy",
    "source_census_particle_count_residual",
    "source_census_macro_weight_residual",
    "escaped_particles_included_in_exposure",
}
_RUNTIME_KEYS = {
    "schema_version",
    "record_type",
    "successor_id",
    "qualification_effect",
    "authorization",
    "exact_normalization",
    "cycle_coverage",
    "all_cycle_extrema",
    "particle_exposure",
    "boundary_escape_ledger",
}


def _validate_runtime_time_escape(value: object) -> dict[str, Any]:
    runtime = _exact_keys(value, _RUNTIME_KEYS, "runtime time/escape evidence")
    _require(
        type(runtime["schema_version"]) is int
        and runtime["schema_version"] == SCHEMA_VERSION,
        "runtime schema version drifted",
    )
    _require(
        type(runtime["record_type"]) is str
        and runtime["record_type"] == RUNTIME_RECORD_TYPE,
        "runtime record type drifted",
    )
    _require(
        type(runtime["successor_id"]) is str
        and runtime["successor_id"] == SUCCESSOR_ID,
        "runtime successor id drifted",
    )
    _require(
        type(runtime["qualification_effect"]) is str
        and runtime["qualification_effect"] == QUALIFICATION_EFFECT,
        "runtime effect drifted",
    )
    authorization = _validate_authorization(runtime["authorization"], "runtime authorization")
    normalization = _validate_exact_normalization(runtime["exact_normalization"])

    coverage = _exact_keys(runtime["cycle_coverage"], _COVERAGE_KEYS, "cycle coverage")
    complete = _strict_bool(coverage["complete"], "cycle coverage complete")
    _require(
        coverage["sampling_mode"]
        == "every_integrator_cycle_post_startup_removal_through_terminal_time",
        "cycle coverage sampling mode drifted",
    )
    start = _finite_scalar(
        coverage["post_startup_removal_start_time"], "cycle coverage start", minimum=0.0
    )
    terminal = _finite_scalar(coverage["terminal_time"], "cycle coverage terminal", minimum=0.0)
    first_cycle = _nonnegative_int(coverage["first_cycle"], "first cycle")
    last_cycle = _nonnegative_int(coverage["last_cycle"], "last cycle")
    covered = _nonnegative_int(coverage["covered_cycle_count"], "covered cycle count")
    expected = _nonnegative_int(
        coverage["expected_contiguous_cycle_count"], "expected contiguous cycle count"
    )
    gaps = _nonnegative_int(coverage["gap_count"], "cycle gap count")
    segments = _nonnegative_int(coverage["restart_segment_count"], "restart segment count")
    segments_complete = _strict_bool(
        coverage["restart_segments_complete"], "restart segments complete"
    )
    _require(
        start == STARTUP_REMOVAL_TIME and terminal == EXPECTED_TERMINAL_TIME,
        "runtime coverage must span exact post-startup interval through t=1200",
    )
    _require(last_cycle >= first_cycle, "runtime cycle interval is reversed")
    _require(
        covered == expected == last_cycle - first_cycle + 1,
        "runtime cycle counts do not prove contiguous coverage",
    )
    _require(segments >= 1, "runtime coverage requires at least one restart segment")

    extrema = _exact_keys(runtime["all_cycle_extrema"], _EXTREMA_KEYS, "all-cycle extrema")
    decoded_extrema = {
        key: _finite_scalar(extrema[key], f"all-cycle {key}", minimum=0.0)
        for key in _EXTREMA_KEYS
    }

    exposure = _exact_keys(
        runtime["particle_exposure"], _PARTICLE_EXPOSURE_KEYS, "particle exposure"
    )
    exposure_complete = _strict_bool(exposure["complete"], "particle exposure complete")
    _require(
        exposure["method"]
        == "every_particle_update_and_pre_destruction_boundary_event",
        "particle exposure method drifted",
    )
    update_included = _strict_bool(
        exposure["active_particle_updates_included"], "active particle updates included"
    )
    boundary_included = _strict_bool(
        exposure["pre_destruction_boundary_events_included"],
        "pre-destruction boundary events included",
    )
    escaped_included = _strict_bool(
        exposure["escaped_particles_included"], "escaped particles included"
    )
    observations = _nonnegative_int(exposure["observation_count"], "particle observations")
    _require(observations > 0, "particle exposure requires observations")
    exposure_statistics = {
        key: _finite_scalar(exposure[key], f"particle exposure {key}", minimum=0.0)
        for key in (
            "maximum_sampled_R",
            "maximum_sampled_Lambda",
            "cumulative_macro_weighted_R_exceedance_fraction",
            "cumulative_CR_energy_weighted_R_exceedance_fraction",
            "cumulative_macro_weighted_Lambda_exceedance_fraction",
            "cumulative_CR_energy_weighted_Lambda_exceedance_fraction",
        )
    }
    for key in (
        "cumulative_macro_weighted_R_exceedance_fraction",
        "cumulative_CR_energy_weighted_R_exceedance_fraction",
        "cumulative_macro_weighted_Lambda_exceedance_fraction",
        "cumulative_CR_energy_weighted_Lambda_exceedance_fraction",
    ):
        _require(exposure_statistics[key] <= 1.0, f"particle exposure {key} exceeds unity")
    r_roundoff = 1.0e-12 * max(1.0, decoded_extrema["R_maximum"])
    lambda_roundoff = 1.0e-12 * max(1.0, decoded_extrema["Lambda_maximum"])
    _require(
        exposure_statistics["maximum_sampled_R"]
        <= decoded_extrema["R_maximum"] + r_roundoff
        and exposure_statistics["maximum_sampled_Lambda"]
        <= decoded_extrema["Lambda_maximum"] + lambda_roundoff,
        "particle exposure maxima exceed all-cycle cell maxima",
    )
    if decoded_extrema["R_maximum"] <= R_MAXIMUM:
        _require(
            exposure_statistics["cumulative_macro_weighted_R_exceedance_fraction"] == 0.0
            and exposure_statistics[
                "cumulative_CR_energy_weighted_R_exceedance_fraction"
            ]
            == 0.0,
            "particle R exposure reports exceedance below the all-cycle cell maximum bound",
        )
    if decoded_extrema["Lambda_maximum"] <= LAMBDA_MAXIMUM:
        _require(
            exposure_statistics[
                "cumulative_macro_weighted_Lambda_exceedance_fraction"
            ]
            == 0.0
            and exposure_statistics[
                "cumulative_CR_energy_weighted_Lambda_exceedance_fraction"
            ]
            == 0.0,
            "particle Lambda exposure reports exceedance below the all-cycle cell maximum bound",
        )

    escape = _exact_keys(
        runtime["boundary_escape_ledger"], _ESCAPE_KEYS, "boundary escape ledger"
    )
    escape_complete = _strict_bool(escape["complete"], "escape ledger complete")
    _require(
        escape["scope"] == "all_Q011_shock_injected_particles_including_startup_removal",
        "escape ledger scope drifted",
    )
    faces = _exact_keys(
        escape["nonperiodic_faces"], {"ix1", "ox1"}, "nonperiodic escape faces"
    )
    decoded_faces: dict[str, dict[str, float | int]] = {}
    for face in ("ix1", "ox1"):
        values = _exact_keys(faces[face], _FACE_KEYS, f"{face} escape ledger")
        decoded_faces[face] = {
            "escaped_particle_count": _nonnegative_int(
                values["escaped_particle_count"], f"{face} escaped count"
            ),
            "escaped_macro_weight": _finite_scalar(
                values["escaped_macro_weight"], f"{face} escaped macro weight", minimum=0.0
            ),
            "escaped_kinetic_energy": _finite_scalar(
                values["escaped_kinetic_energy"], f"{face} escaped energy", minimum=0.0
            ),
        }
    _require(
        escape["periodic_faces"] == ["ix2", "ox2"],
        "escape ledger periodic face inventory drifted",
    )
    total_count = _nonnegative_int(
        escape["total_escaped_particle_count"], "total escaped particle count"
    )
    total_weight = _finite_scalar(
        escape["total_escaped_macro_weight"], "total escaped macro weight", minimum=0.0
    )
    total_energy = _finite_scalar(
        escape["total_escaped_kinetic_energy"], "total escaped energy", minimum=0.0
    )
    _require(
        total_count == sum(int(decoded_faces[face]["escaped_particle_count"]) for face in decoded_faces),
        "escape face counts do not close to total",
    )
    _require(
        total_weight == sum(float(decoded_faces[face]["escaped_macro_weight"]) for face in decoded_faces)
        and total_energy
        == sum(float(decoded_faces[face]["escaped_kinetic_energy"]) for face in decoded_faces),
        "escape face weights or energy do not close to total",
    )
    unaccounted_count = _nonnegative_int(
        escape["unaccounted_particle_count"], "unaccounted particle count"
    )
    unaccounted_weight = _finite_scalar(
        escape["unaccounted_macro_weight"], "unaccounted macro weight", minimum=0.0
    )
    unaccounted_energy = _finite_scalar(
        escape["unaccounted_kinetic_energy"], "unaccounted energy", minimum=0.0
    )
    census_count_residual = _finite_scalar(
        escape["source_census_particle_count_residual"], "source census count residual"
    )
    census_weight_residual = _finite_scalar(
        escape["source_census_macro_weight_residual"], "source census weight residual"
    )
    escape_exposure = _strict_bool(
        escape["escaped_particles_included_in_exposure"],
        "escape ledger particles included in exposure",
    )
    time_pass = (
        complete
        and gaps == 0
        and segments_complete
        and exposure_complete
        and update_included
        and boundary_included
        and escaped_included
        and escape_complete
        and escape_exposure
        and unaccounted_count == 0
        and unaccounted_weight == 0.0
        and unaccounted_energy == 0.0
        and census_count_residual == 0.0
        and census_weight_residual == 0.0
    )
    return {
        "exact_normalization": normalization,
        "authorization": authorization,
        "cycle_coverage": dict(coverage),
        "all_cycle_extrema": decoded_extrema,
        "particle_exposure": {
            **dict(exposure),
            **exposure_statistics,
        },
        "boundary_escape_ledger": {
            **dict(escape),
            "nonperiodic_faces": decoded_faces,
        },
        "time_and_escape_complete": time_pass,
    }


def _validate_snapshot_record(value: object) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), "snapshot applicability record must be a mapping")
    _require(
        type(value.get("schema_version")) is int
        and value.get("schema_version") == SCHEMA_VERSION,
        "snapshot schema version drifted",
    )
    _require(
        type(value.get("record_type")) is str
        and value.get("record_type") == SNAPSHOT_RECORD_TYPE,
        "snapshot record type drifted",
    )
    _require(
        type(value.get("successor_id")) is str
        and value.get("successor_id") == SUCCESSOR_ID,
        "snapshot successor id drifted",
    )
    _require(
        type(value.get("qualification_effect")) is str
        and value.get("qualification_effect") == QUALIFICATION_EFFECT,
        "snapshot qualification effect drifted",
    )
    _validate_authorization(value.get("authorization"), "snapshot authorization")
    _validate_exact_normalization(value.get("exact_normalization"))
    gates = _exact_keys(
        value.get("gates"),
        {
            "Q011-APP-NORM",
            "Q011-APP-R",
            "Q011-APP-LAMBDA",
            "Q011-APP-DI",
            "Q011-APP-RG",
            "Q011-APP-TIME",
        },
        "snapshot gates",
    )
    _require(gates["Q011-APP-NORM"]["pass"] is True, "snapshot normalization gate drifted")
    _require(gates["Q011-APP-TIME"]["pass"] is False, "snapshot time gate must remain false")
    r_gate = gates["Q011-APP-R"]
    lambda_gate = gates["Q011-APP-LAMBDA"]
    di_gate = gates["Q011-APP-DI"]
    rg_gate = gates["Q011-APP-RG"]
    expected = {
        "Q011-APP-R": (
            _finite_scalar(r_gate["observed_maximum"], "snapshot R maximum", minimum=0.0)
            <= R_MAXIMUM
            and _exact_numeric(
                r_gate["required_maximum"], R_MAXIMUM, "snapshot required R maximum"
            )
            == R_MAXIMUM
        ),
        "Q011-APP-LAMBDA": (
            _finite_scalar(
                lambda_gate["observed_maximum"], "snapshot Lambda maximum", minimum=0.0
            )
            <= LAMBDA_MAXIMUM
            and _exact_numeric(
                lambda_gate["required_maximum"],
                LAMBDA_MAXIMUM,
                "snapshot required Lambda maximum",
            )
            == LAMBDA_MAXIMUM
        ),
        "Q011-APP-DI": (
            _finite_scalar(
                di_gate["observed_S_delta_minimum_excluding_shock_transition"],
                "snapshot S_delta minimum",
                minimum=0.0,
            )
            >= S_DELTA_MINIMUM
            and _finite_scalar(
                di_gate["observed_lambda_B_characteristic_over_local_di_maximum"],
                "snapshot lambda_B over d_i maximum",
                minimum=0.0,
            )
            >= LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM
            and _finite_scalar(
                di_gate["observed_sub_10di_magnetic_power_fraction"],
                "snapshot sub-10di power fraction",
                minimum=0.0,
            )
            <= SUB_10DI_POWER_FRACTION_MAXIMUM
            and _exact_numeric(
                di_gate["required_S_delta_minimum"],
                S_DELTA_MINIMUM,
                "snapshot required S_delta minimum",
            )
            == S_DELTA_MINIMUM
            and _exact_numeric(
                di_gate["required_lambda_B_characteristic_over_local_di_maximum"],
                LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM,
                "snapshot required lambda_B over d_i maximum",
            )
            == LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM
            and _exact_numeric(
                di_gate["required_sub_10di_magnetic_power_fraction_maximum"],
                SUB_10DI_POWER_FRACTION_MAXIMUM,
                "snapshot required sub-10di power fraction maximum",
            )
            == SUB_10DI_POWER_FRACTION_MAXIMUM
        ),
        "Q011-APP-RG": (
            _finite_scalar(
                rg_gate["observed_macro_q999_over_Ly"],
                "snapshot macro q999 over Ly",
                minimum=0.0,
            )
            <= RG_Q999_OVER_LY_MAXIMUM
            and _finite_scalar(
                rg_gate["observed_energy_q999_over_Ly"],
                "snapshot energy q999 over Ly",
                minimum=0.0,
            )
            <= RG_Q999_OVER_LY_MAXIMUM
            and _finite_scalar(
                rg_gate["observed_maximum_over_Ly"],
                "snapshot maximum gyroradius over Ly",
                minimum=0.0,
            )
            <= RG_MAXIMUM_OVER_LY_MAXIMUM
            and _finite_scalar(
                rg_gate["observed_energy_fraction_with_rg_above_Ly_over_4"],
                "snapshot large-gyroradius energy fraction",
                minimum=0.0,
            )
            <= RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM
            and _exact_numeric(
                rg_gate["required_q999_over_Ly_maximum"],
                RG_Q999_OVER_LY_MAXIMUM,
                "snapshot required q999 over Ly maximum",
            )
            == RG_Q999_OVER_LY_MAXIMUM
            and _exact_numeric(
                rg_gate["required_maximum_over_Ly"],
                RG_MAXIMUM_OVER_LY_MAXIMUM,
                "snapshot required maximum gyroradius over Ly",
            )
            == RG_MAXIMUM_OVER_LY_MAXIMUM
            and _exact_numeric(
                rg_gate["required_energy_fraction_with_rg_above_Ly_over_4_maximum"],
                RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM,
                "snapshot required large-gyroradius energy fraction maximum",
            )
            == RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM
        ),
    }
    for name, expected_pass in expected.items():
        _require(
            type(gates[name]["pass"]) is bool and gates[name]["pass"] is expected_pass,
            f"snapshot {name} pass flag disagrees with bound observables",
        )
    return value


@_public_contract("Q011 physical-applicability history reduction")
def reduce_physical_applicability_history(
    snapshots: Sequence[Mapping[str, Any]], runtime_time_escape_evidence: object
) -> dict[str, Any]:
    """Combine snapshot gates with future all-cycle time/escape evidence."""
    _require(
        isinstance(snapshots, Sequence) and not isinstance(snapshots, (str, bytes)),
        "snapshot applicability history must be a sequence",
    )
    _require(len(snapshots) > 0, "snapshot applicability history must not be empty")
    records = [_validate_snapshot_record(record) for record in snapshots]
    times = [
        _finite_scalar(record["observed_committed_time"], "snapshot observed time", minimum=0.0)
        for record in records
    ]
    _require(
        all(right > left for left, right in zip(times, times[1:])),
        "snapshot applicability history times must strictly increase",
    )
    runtime = _validate_runtime_time_escape(runtime_time_escape_evidence)
    extrema = runtime["all_cycle_extrema"]
    _require(
        times[0] >= runtime["cycle_coverage"]["post_startup_removal_start_time"]
        and times[-1] <= runtime["cycle_coverage"]["terminal_time"],
        "snapshot applicability history escaped runtime coverage",
    )
    snapshot_extrema = {
        "R_maximum": max(
            float(record["gates"]["Q011-APP-R"]["observed_maximum"]) for record in records
        ),
        "Lambda_maximum": max(
            float(record["gates"]["Q011-APP-LAMBDA"]["observed_maximum"])
            for record in records
        ),
        "S_delta_minimum_excluding_shock_transition": min(
            float(
                record["gates"]["Q011-APP-DI"][
                    "observed_S_delta_minimum_excluding_shock_transition"
                ]
            )
            for record in records
        ),
        "lambda_B_characteristic_over_local_di_maximum_minimum": min(
            float(
                record["gates"]["Q011-APP-DI"][
                    "observed_lambda_B_characteristic_over_local_di_maximum"
                ]
            )
            for record in records
        ),
        "sub_10di_magnetic_power_fraction_maximum": max(
            float(
                record["gates"]["Q011-APP-DI"][
                    "observed_sub_10di_magnetic_power_fraction"
                ]
            )
            for record in records
        ),
        "macro_q999_over_Ly_maximum": max(
            float(record["gates"]["Q011-APP-RG"]["observed_macro_q999_over_Ly"])
            for record in records
        ),
        "energy_q999_over_Ly_maximum": max(
            float(record["gates"]["Q011-APP-RG"]["observed_energy_q999_over_Ly"])
            for record in records
        ),
        "particle_rg_maximum_over_Ly": max(
            float(record["gates"]["Q011-APP-RG"]["observed_maximum_over_Ly"])
            for record in records
        ),
        "energy_fraction_with_rg_above_Ly_over_4_maximum": max(
            float(
                record["gates"]["Q011-APP-RG"][
                    "observed_energy_fraction_with_rg_above_Ly_over_4"
                ]
            )
            for record in records
        ),
    }
    _require(
        extrema["R_maximum"] >= snapshot_extrema["R_maximum"]
        and extrema["Lambda_maximum"] >= snapshot_extrema["Lambda_maximum"]
        and extrema["S_delta_minimum_excluding_shock_transition"]
        <= snapshot_extrema["S_delta_minimum_excluding_shock_transition"]
        and extrema["lambda_B_characteristic_over_local_di_maximum_minimum"]
        <= snapshot_extrema["lambda_B_characteristic_over_local_di_maximum_minimum"]
        and extrema["sub_10di_magnetic_power_fraction_maximum"]
        >= snapshot_extrema["sub_10di_magnetic_power_fraction_maximum"]
        and extrema["macro_q999_over_Ly_maximum"]
        >= snapshot_extrema["macro_q999_over_Ly_maximum"]
        and extrema["energy_q999_over_Ly_maximum"]
        >= snapshot_extrema["energy_q999_over_Ly_maximum"]
        and extrema["particle_rg_maximum_over_Ly"]
        >= snapshot_extrema["particle_rg_maximum_over_Ly"]
        and extrema["energy_fraction_with_rg_above_Ly_over_4_maximum"]
        >= snapshot_extrema["energy_fraction_with_rg_above_Ly_over_4_maximum"],
        "all-cycle runtime extrema do not conservatively dominate supplied snapshots",
    )
    snapshot_gate_names = ("Q011-APP-NORM", "Q011-APP-R", "Q011-APP-LAMBDA", "Q011-APP-DI", "Q011-APP-RG")
    snapshot_pass = {
        name: all(bool(record["gates"][name]["pass"]) for record in records)
        for name in snapshot_gate_names
    }
    gates: dict[str, dict[str, Any]] = {
        "Q011-APP-NORM": {
            "pass": snapshot_pass["Q011-APP-NORM"],
            "basis": "exact_normalization_in_every_snapshot_and_runtime_record",
        },
        "Q011-APP-R": {
            "pass": snapshot_pass["Q011-APP-R"] and extrema["R_maximum"] <= R_MAXIMUM,
            "all_cycle_observed_maximum": extrema["R_maximum"],
            "required_maximum": R_MAXIMUM,
        },
        "Q011-APP-LAMBDA": {
            "pass": (
                snapshot_pass["Q011-APP-LAMBDA"]
                and extrema["Lambda_maximum"] <= LAMBDA_MAXIMUM
            ),
            "all_cycle_observed_maximum": extrema["Lambda_maximum"],
            "required_maximum": LAMBDA_MAXIMUM,
        },
        "Q011-APP-DI": {
            "pass": (
                snapshot_pass["Q011-APP-DI"]
                and extrema["S_delta_minimum_excluding_shock_transition"]
                >= S_DELTA_MINIMUM
                and extrema["lambda_B_characteristic_over_local_di_maximum_minimum"]
                >= LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM
                and extrema["sub_10di_magnetic_power_fraction_maximum"]
                <= SUB_10DI_POWER_FRACTION_MAXIMUM
            ),
            "all_cycle_extrema": {
                key: extrema[key]
                for key in (
                    "S_delta_minimum_excluding_shock_transition",
                    "lambda_B_characteristic_over_local_di_maximum_minimum",
                    "sub_10di_magnetic_power_fraction_maximum",
                )
            },
        },
        "Q011-APP-RG": {
            "pass": (
                snapshot_pass["Q011-APP-RG"]
                and extrema["macro_q999_over_Ly_maximum"] <= RG_Q999_OVER_LY_MAXIMUM
                and extrema["energy_q999_over_Ly_maximum"] <= RG_Q999_OVER_LY_MAXIMUM
                and extrema["particle_rg_maximum_over_Ly"]
                <= RG_MAXIMUM_OVER_LY_MAXIMUM
                and extrema["energy_fraction_with_rg_above_Ly_over_4_maximum"]
                <= RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM
            ),
            "all_cycle_extrema": {
                key: extrema[key]
                for key in (
                    "macro_q999_over_Ly_maximum",
                    "energy_q999_over_Ly_maximum",
                    "particle_rg_maximum_over_Ly",
                    "energy_fraction_with_rg_above_Ly_over_4_maximum",
                )
            },
        },
        "Q011-APP-TIME": {
            "pass": runtime["time_and_escape_complete"],
            "basis": (
                "complete_contiguous_per_cycle_post_startup_coverage_particle_"
                "exposure_and_closed_boundary_escape_ledger"
            ),
        },
    }
    maximum_lambda = extrema["Lambda_maximum"]
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": HISTORY_RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "authorization": dict(AUTHORIZATION),
        "snapshot_count": len(records),
        "first_observed_committed_time": times[0],
        "last_observed_committed_time": times[-1],
        "runtime_time_escape_evidence": runtime,
        "snapshot_extrema": snapshot_extrema,
        "gates": gates,
        "all_physical_applicability_gates_pass": all(
            bool(gate["pass"]) for gate in gates.values()
        ),
        "claim_rejection_rules": {
            gate: list(claims) for gate, claims in CLAIM_REJECTION_RULES.items()
        },
        "claim_rejections": _claim_rejections(gates, maximum_lambda),
        "permanent_claim_exclusions": list(PERMANENT_CLAIM_EXCLUSIONS),
        "authorization_effect": "none_even_if_all_gates_pass",
    }


__all__ = [
    "AUTHORIZATION",
    "ApplicabilitySnapshot",
    "CLAIM_REJECTION_RULES",
    "DETECTED_FRONT_REGIONS",
    "EXACT_NORMALIZATION",
    "HISTORY_RECORD_TYPE",
    "LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM",
    "LAMBDA_MAXIMUM",
    "PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES",
    "PERMANENT_CLAIM_EXCLUSIONS",
    "PhysicalApplicabilityError",
    "QUALIFICATION_EFFECT",
    "RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM",
    "RG_MAXIMUM_OVER_LY_MAXIMUM",
    "RG_Q999_OVER_LY_MAXIMUM",
    "RUNTIME_RECORD_TYPE",
    "R_MAXIMUM",
    "SCHEMA_VERSION",
    "SNAPSHOT_RECORD_TYPE",
    "STARTUP_REMOVAL_TIME",
    "SUB_10DI_POWER_FRACTION_MAXIMUM",
    "SUCCESSOR_ID",
    "S_DELTA_MINIMUM",
    "reduce_physical_applicability_history",
    "reduce_physical_applicability_snapshot",
]
