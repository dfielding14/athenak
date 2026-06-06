#!/usr/bin/env python3
"""Fail-closed Q019 particle-state diagnostics from trusted schema-7 extraction.

This reducer intentionally does not parse AthenaK restart bytes.  A separate
trusted extractor must bind the exact retained schema-7 artifact and provide
the arrays consumed here.  The reducer has no launch or scientific authority.
"""

from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np


SCHEMA_VERSION = 1
RECORD_TYPE = "q019_nonlinear_bell_schema7_particle_state_reduction"
CAMPAIGN_ID = "Q019-HR-VOLUME-AWARE-NOHALL"
CLAIM_ID = "CLAIM-PROD-BELL-NONLINEAR-NOHALL-001"
RESTART_SCHEMA = 7
STATE_KIND = "mass_normalized_momentum"


class ParticleStateError(ValueError):
    """Raised when a Q019 particle-state reduction cannot be trusted."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ParticleStateError(message)


def _positive_float(value: object, *, label: str) -> float:
    _require(
        type(value) is float and math.isfinite(value) and value > 0.0,
        f"{label} must be a finite positive float",
    )
    return value


def _finite_vector(value: object, *, label: str, length: int | None = None) -> np.ndarray:
    array = np.asarray(value)
    _require(
        np.issubdtype(array.dtype, np.floating),
        f"{label} must use a floating dtype",
    )
    array = array.astype(np.float64, copy=False)
    _require(array.ndim == 1, f"{label} must be one-dimensional")
    if length is not None:
        _require(array.shape == (length,), f"{label} has the wrong length")
    _require(np.all(np.isfinite(array)), f"{label} must be finite")
    return array


def _finite_vectors(value: object, *, label: str, count: int | None = None) -> np.ndarray:
    array = np.asarray(value)
    _require(
        np.issubdtype(array.dtype, np.floating),
        f"{label} must use a floating dtype",
    )
    array = array.astype(np.float64, copy=False)
    _require(array.ndim == 2 and array.shape[1] == 3, f"{label} must have shape (n, 3)")
    if count is not None:
        _require(array.shape[0] == count, f"{label} has the wrong particle count")
    _require(np.all(np.isfinite(array)), f"{label} must be finite")
    return array


def _species_indices(value: object, *, count: int, nspecies: int) -> np.ndarray:
    species = np.asarray(value)
    _require(np.issubdtype(species.dtype, np.integer), "species must use an integer dtype")
    species = species.astype(np.int64, copy=False)
    _require(species.shape == (count,), "species has the wrong particle count")
    _require(np.all((species >= 0) & (species < nspecies)), "species index is out of range")
    return species


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values, kind="stable")
    sorted_values = values[order]
    cumulative = np.cumsum(weights[order])
    target = 0.5 * float(cumulative[-1])
    return float(sorted_values[int(np.searchsorted(cumulative, target, side="left"))])


def _mhd_budget(value: object) -> dict[str, object]:
    _require(type(value) is dict, "mhd_budget must be an object")
    expected = {
        "gas_momentum",
        "gas_kinetic_energy",
        "gas_thermal_energy",
        "magnetic_energy",
        "mhd_total_energy",
    }
    _require(set(value) == expected, "mhd_budget keys drifted")
    momentum = _finite_vector(value["gas_momentum"], label="mhd_budget/gas_momentum", length=3)
    energies = {}
    for key in expected - {"gas_momentum"}:
        parsed = value[key]
        _require(
            type(parsed) is float and math.isfinite(parsed) and parsed >= 0.0,
            f"mhd_budget/{key} must be a finite non-negative float",
        )
        energies[key] = parsed
    component_total = (
        energies["gas_kinetic_energy"]
        + energies["gas_thermal_energy"]
        + energies["magnetic_energy"]
    )
    _require(
        math.isclose(
            energies["mhd_total_energy"],
            component_total,
            rel_tol=1.0e-12,
            abs_tol=1.0e-14,
        ),
        "mhd_budget total energy is inconsistent with its components",
    )
    return {"gas_momentum": momentum.tolist(), **energies}


def _reference_budget(value: object) -> dict[str, object]:
    _require(type(value) is dict, "reference_budget must be an object")
    _require(
        set(value) == {"total_momentum", "total_energy"},
        "reference_budget keys drifted",
    )
    momentum = _finite_vector(
        value["total_momentum"], label="reference_budget/total_momentum", length=3
    )
    energy = value["total_energy"]
    _require(
        type(energy) is float and math.isfinite(energy) and energy >= 0.0,
        "reference_budget/total_energy must be a finite non-negative float",
    )
    return {"total_momentum": momentum.tolist(), "total_energy": energy}


def reduce_schema7_particle_state(
    *,
    restart_schema: int,
    state_kind: str,
    deltaf_mode: str,
    momentum_per_mass: object,
    sampled_magnetic_field: object,
    macro_weight: object,
    species: object,
    particle_q_over_mc: object,
    species_mass: object,
    species_q_over_mc: object,
    deposit_qscale: float,
    artificial_light_speed: float,
    domain_volume: float,
    gas_bulk_velocity: object,
    guide_field_direction: object,
    mhd_budget: object,
    reference_budget: object | None = None,
) -> dict[str, Any]:
    """Reduce a trusted full-f schema-7 particle state into physical diagnostics."""
    _require(type(restart_schema) is int and restart_schema == RESTART_SCHEMA,
             "only restart schema 7 is accepted")
    _require(state_kind == STATE_KIND, "particle state must be mass-normalized momentum")
    _require(deltaf_mode == "off", "Q019 full-f reducer rejects delta-f particle state")
    qscale = _positive_float(deposit_qscale, label="deposit_qscale")
    light_speed = _positive_float(artificial_light_speed, label="artificial_light_speed")
    volume = _positive_float(domain_volume, label="domain_volume")

    momentum = _finite_vectors(momentum_per_mass, label="momentum_per_mass")
    count = momentum.shape[0]
    _require(count > 0, "particle state must not be empty")
    sampled_b = _finite_vectors(
        sampled_magnetic_field, label="sampled_magnetic_field", count=count
    )
    weights = _finite_vector(macro_weight, label="macro_weight", length=count)
    _require(np.all(weights > 0.0), "macro_weight must be strictly positive")
    masses = _finite_vector(species_mass, label="species_mass")
    q_over_mc_by_species = _finite_vector(
        species_q_over_mc, label="species_q_over_mc", length=masses.size
    )
    _require(np.all(masses > 0.0), "species_mass must be strictly positive")
    _require(np.all(q_over_mc_by_species != 0.0), "species_q_over_mc must be nonzero")
    species_index = _species_indices(species, count=count, nspecies=masses.size)
    particle_qom = _finite_vector(
        particle_q_over_mc, label="particle_q_over_mc", length=count
    )
    _require(
        np.array_equal(particle_qom, q_over_mc_by_species[species_index]),
        "particle q/(mc) does not match the bound species table",
    )
    gas_velocity = _finite_vector(gas_bulk_velocity, label="gas_bulk_velocity", length=3)
    guide = _finite_vector(guide_field_direction, label="guide_field_direction", length=3)
    guide_norm = float(np.linalg.norm(guide))
    _require(guide_norm > 0.0, "guide_field_direction must be nonzero")
    guide /= guide_norm
    parsed_mhd = _mhd_budget(mhd_budget)

    momentum_squared = np.sum(momentum * momentum, axis=1)
    gamma = np.sqrt(1.0 + momentum_squared / (light_speed * light_speed))
    velocity = momentum / gamma[:, None]
    macro_mass = qscale * weights * masses[species_index]
    total_macro_mass = float(np.sum(macro_mass))
    mass_density = total_macro_mass / volume
    bulk_velocity = np.sum(macro_mass[:, None] * velocity, axis=0) / total_macro_mass

    particle_momentum = np.sum(macro_mass[:, None] * momentum, axis=0)
    kinetic_per_mass = (gamma - 1.0) * light_speed * light_speed
    particle_kinetic_energy = float(np.sum(macro_mass * kinetic_per_mass))
    momentum_flux = np.einsum("n,ni,nj->ij", macro_mass, momentum, velocity) / volume
    centered_velocity = velocity - bulk_velocity
    velocity_pressure = (
        np.einsum("n,ni,nj->ij", macro_mass, centered_velocity, centered_velocity)
        / volume
    )

    charge_over_c = macro_mass * particle_qom
    charge_density_over_c = float(np.sum(charge_over_c) / volume)
    current_over_c = np.sum(charge_over_c[:, None] * velocity, axis=0) / volume
    gas_frame_current_over_c = current_over_c - charge_density_over_c * gas_velocity
    charge_density = light_speed * charge_density_over_c
    lab_current = light_speed * current_over_c
    gas_frame_current = lab_current - charge_density * gas_velocity

    magnetic_norm = np.linalg.norm(sampled_b, axis=1)
    _require(np.all(magnetic_norm > 0.0), "sampled magnetic field must be nonzero")
    magnetic_unit = sampled_b / magnetic_norm[:, None]
    momentum_parallel_b = np.sum(momentum * magnetic_unit, axis=1)
    momentum_perpendicular = momentum - momentum_parallel_b[:, None] * magnetic_unit
    gyroradius = (
        np.linalg.norm(momentum_perpendicular, axis=1)
        / (np.abs(particle_qom) * magnetic_norm)
    )
    _require(np.all(np.isfinite(gyroradius)), "computed gyroradius must be finite")

    parallel_flux = float(guide @ momentum_flux @ guide)
    perpendicular_flux = float((np.trace(momentum_flux) - parallel_flux) / 2.0)
    total_momentum = np.asarray(parsed_mhd["gas_momentum"]) + particle_momentum
    total_energy = parsed_mhd["mhd_total_energy"] + particle_kinetic_energy

    conservation: dict[str, object] = {
        "gas_plus_cr_total_momentum": total_momentum.tolist(),
        "gas_plus_cr_total_energy": total_energy,
        "reference_bound": reference_budget is not None,
    }
    if reference_budget is not None:
        reference = _reference_budget(reference_budget)
        reference_momentum = np.asarray(reference["total_momentum"])
        momentum_residual = total_momentum - reference_momentum
        energy_residual = total_energy - reference["total_energy"]
        conservation.update(
            {
                "reference_total_momentum": reference["total_momentum"],
                "reference_total_energy": reference["total_energy"],
                "momentum_residual": momentum_residual.tolist(),
                "momentum_residual_l2": float(np.linalg.norm(momentum_residual)),
                "energy_residual": energy_residual,
                "relative_energy_residual": (
                    None
                    if reference["total_energy"] == 0.0
                    else energy_residual / reference["total_energy"]
                ),
            }
        )

    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "campaign_id": CAMPAIGN_ID,
        "claim_id": CLAIM_ID,
        "status": "source_local_reduction_only_not_qualified",
        "qualification_effect": "none",
        "launch_authorized": False,
        "scientific_claim_authorized": False,
        "restart_schema": RESTART_SCHEMA,
        "state_kind": STATE_KIND,
        "particle_count": count,
        "domain_volume": volume,
        "particle_mass": {
            "total_macro_mass": total_macro_mass,
            "volume_averaged_mass_density": mass_density,
            "bulk_velocity": bulk_velocity.tolist(),
        },
        "particle_momentum": {
            "volume_integrated_vector": particle_momentum.tolist(),
            "momentum_flux_tensor": momentum_flux.tolist(),
            "velocity_pressure_tensor": velocity_pressure.tolist(),
            "guide_parallel_momentum_flux": parallel_flux,
            "mean_guide_perpendicular_momentum_flux": perpendicular_flux,
        },
        "particle_kinetic_energy": {
            "volume_integrated": particle_kinetic_energy,
            "volume_averaged_density": particle_kinetic_energy / volume,
        },
        "current": {
            "charge_density_over_c": charge_density_over_c,
            "deposited_current_j_over_c": current_over_c.tolist(),
            "gas_frame_deposited_current_j_over_c": gas_frame_current_over_c.tolist(),
            "reconstructed_charge_density": charge_density,
            "reconstructed_lab_current": lab_current.tolist(),
            "reconstructed_gas_frame_current": gas_frame_current.tolist(),
            "lab_current_parallel_to_initial_guide": float(lab_current @ guide),
            "gas_frame_current_parallel_to_initial_guide": float(gas_frame_current @ guide),
        },
        "gyroradius": {
            "definition": "|p_perp/m|/(|q/(mc)| |B_sampled|)",
            "macro_mass_weighted_mean": float(np.average(gyroradius, weights=macro_mass)),
            "macro_mass_weighted_median": _weighted_median(gyroradius, macro_mass),
            "minimum": float(np.min(gyroradius)),
            "maximum": float(np.max(gyroradius)),
        },
        "mhd_budget": parsed_mhd,
        "conservation": conservation,
        "limitations": [
            "trusted schema-7 extraction and exact retained-artifact binding are external prerequisites",
            "sampled particle magnetic field is used for gyroradius reconstruction",
            "full-f deltaf-off states only",
            "no launch, qualification, mechanism classification, or scientific claim authority",
        ],
    }
