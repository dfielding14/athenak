#!/usr/bin/env python3
"""Analyze source-local Q019 physics-first contract fixtures.

Raw science analysis is deliberately disabled until a hardened installed
control-plane/Q043 admission adapter exists. Synthetic fixtures exercise the
diagnostic and fail-closed gate logic only.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from tst.publication import q019_finite_rigidity_early_time_physics_predecessor_v2
from tst.publication import q019_hardened_provenance_boundary_v2 as provenance_boundary
from tst.publication import q019_particle_state_analysis_bridge_v2 as particle_bridge
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as decks
from tst.publication import q019_registered_case_contracts_v1 as case_contracts


RECORD_TYPE = "q019_physics_first_nonlinear_bell_successor_v2_analysis"
CURRENT_AGREEMENT_RELATIVE_TOLERANCE = 1.0e-8
MINIMUM_CHARACTERISTIC_SHELL_RL_OVER_DX = decks.MIN_CHARACTERISTIC_SHELL_RL_OVER_DX
REQUIRED_GRID_FIELDS = (
    "dens",
    "eint",
    "velx",
    "vely",
    "velz",
    "bcc1",
    "bcc2",
    "bcc3",
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


class ContractError(ValueError):
    """Raised when Q019 analysis inputs or gates drift."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def _case_map() -> dict[str, dict[str, object]]:
    try:
        return {
            str(case["case_id"]): case for case in case_contracts.expected_cases()
        }
    except case_contracts.CaseContractError as error:
        raise ContractError(str(error)) from error


def validate_raw_provenance(provenance: Mapping[str, object]) -> dict[str, object]:
    """Validate one hardened registered raw-science admission."""
    try:
        return provenance_boundary.validate_raw_science_admission(provenance)
    except provenance_boundary.ProvenanceBoundaryError as error:
        raise ContractError(str(error)) from error


def _finite_array(
    value: object, *, label: str, shape: tuple[int, ...] | None = None
) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    _require(array.size > 0 and np.all(np.isfinite(array)), f"{label} is invalid")
    if shape is not None:
        _require(array.shape == shape, f"{label} shape drifted")
    return array


def _faces(value: object, *, label: str, count: int) -> np.ndarray:
    faces = _finite_array(value, label=label)
    _require(faces.ndim == 1 and faces.size == count + 1, f"{label} count drifted")
    spacing = np.diff(faces)
    _require(
        np.all(spacing > 0.0)
        and np.allclose(spacing, spacing[0], rtol=1.0e-12, atol=1.0e-14),
        f"{label} must be uniform and increasing",
    )
    return faces


def _snapshot(snapshot: Mapping[str, object]) -> dict[str, object]:
    _require(type(snapshot) is dict, "snapshot must be an object")
    _require(
        set(snapshot) == {
            "cycle",
            "time",
            "x1_faces",
            "x2_faces",
            "x3_faces",
            "fields",
        },
        "snapshot keys drifted",
    )
    cycle = snapshot["cycle"]
    time = snapshot["time"]
    _require(type(cycle) is int and cycle >= 0, "snapshot cycle is invalid")
    _require(type(time) is float and math.isfinite(time), "snapshot time is invalid")
    fields = snapshot["fields"]
    _require(type(fields) is dict and set(fields) == set(REQUIRED_GRID_FIELDS), "field inventory drifted")
    first = _finite_array(fields["dens"], label="dens")
    _require(first.ndim == 3, "fields must use (x3,x2,x1)")
    parsed = {"dens": first}
    for name in REQUIRED_GRID_FIELDS[1:]:
        parsed[name] = _finite_array(fields[name], label=name, shape=first.shape)
    _require(np.all(parsed["dens"] > 0.0), "density must remain positive")
    faces = (
        _faces(snapshot["x1_faces"], label="x1_faces", count=first.shape[2]),
        _faces(snapshot["x2_faces"], label="x2_faces", count=first.shape[1]),
        _faces(snapshot["x3_faces"], label="x3_faces", count=first.shape[0]),
    )
    return {
        "cycle": cycle,
        "time": time,
        "fields": parsed,
        "faces": faces,
        "shape": first.shape,
    }


def _particle_state(record: Mapping[str, object], *, case: Mapping[str, object]) -> dict[str, object]:
    _require(type(record) is dict, "particle state must be an object")
    _require(
        record.get("record_type") == particle_bridge.RECORD_TYPE
        and record.get("case_id") == case["case_id"]
        and record.get("campaign_id") == case["campaign_id"],
        "particle state identity drifted",
    )
    _require(
        record.get("authority")
        == {
            "launch_authorized": False,
            "policy_authorized": False,
            "qualification_authorized": False,
            "claim_authorized": False,
            "raw_production_authorized": False,
            "nonlinear_saturation_claim_authorized": False,
        },
        "particle state authority drifted",
    )
    time = record.get("time")
    cycle = record.get("cycle")
    _require(type(cycle) is int and cycle >= 0, "particle cycle is invalid")
    _require(type(time) is float and math.isfinite(time), "particle time is invalid")
    bulk = _finite_array(record.get("particle_bulk_velocity"), label="particle bulk velocity")
    _require(bulk.shape == (3,), "particle bulk velocity shape drifted")
    momentum = record.get("particle_momentum")
    current = record.get("current")
    gyroradius = record.get("gyroradius")
    conservation = record.get("conservation")
    _require(type(momentum) is dict and type(current) is dict, "particle diagnostics missing")
    _require(type(gyroradius) is dict and gyroradius, "particle gyroradius missing")
    _require(type(conservation) is dict, "particle conservation diagnostics missing")
    vector = _finite_array(momentum["volume_integrated_vector"], label="particle momentum")
    _require(vector.shape == (3,), "particle momentum shape drifted")
    momentum_flux = _finite_array(
        momentum["momentum_flux_tensor"], label="particle momentum flux", shape=(3, 3)
    )
    velocity_pressure = _finite_array(
        momentum["velocity_pressure_tensor"],
        label="particle velocity pressure",
        shape=(3, 3),
    )
    lab_current = _finite_array(current["lab_j_over_c"], label="particle lab current")
    gas_current = _finite_array(current["gas_frame_j_over_c"], label="particle gas current")
    _require(lab_current.shape == gas_current.shape == (3,), "particle current shape drifted")
    return {
        "cycle": cycle,
        "time": time,
        "particle_bulk_velocity": bulk,
        "particle_momentum": vector,
        "particle_momentum_flux_tensor": momentum_flux,
        "particle_velocity_pressure_tensor": velocity_pressure,
        "particle_kinetic_energy": float(record["particle_kinetic_energy"]),
        "lab_j_over_c": lab_current,
        "gas_frame_j_over_c": gas_current,
        "gyroradius": dict(gyroradius),
        "conservation": dict(conservation),
    }


def _available_mode_extents(
    case: Mapping[str, object], case_map: Mapping[str, Mapping[str, object]]
) -> tuple[float, float, float]:
    pair_id = str(case["box_pair_id"])
    if pair_id == "not_applicable":
        return tuple(float(value) for value in case["extents"])
    pair = [
        row
        for row in case_map.values()
        if row["box_pair_id"] == pair_id and row["dimension"] == case["dimension"]
    ]
    _require(pair, "box-pair identity has no matrix rows")
    return tuple(
        min(float(row["extents"][axis]) for row in pair) for axis in range(3)
    )


def _spectrum(
    fields: Mapping[str, np.ndarray],
    dx: tuple[float, float, float],
    *,
    k0: float,
    available_mode_extents: tuple[float, float, float],
) -> dict[str, object]:
    total = np.abs(np.fft.fftn(fields["bcc2"])) ** 2 + np.abs(np.fft.fftn(fields["bcc3"])) ** 2
    total[(0,) * total.ndim] = 0.0
    power = float(np.sum(total))
    _require(power > np.finfo(float).tiny, "transverse magnetic spectrum has no power")
    axes = [
        2.0 * math.pi * np.fft.fftfreq(count, d=spacing)
        for count, spacing in zip(total.shape, reversed(dx))
    ]
    grids = np.meshgrid(*axes, indexing="ij")
    kmag = np.sqrt(sum(grid * grid for grid in grids))
    dominant = np.unravel_index(int(np.argmax(total)), total.shape)
    shared_components = [
        np.isclose(
            grid * small_extent / (2.0 * math.pi),
            np.rint(grid * small_extent / (2.0 * math.pi)),
            rtol=0.0,
            atol=1.0e-10,
        )
        for grid, small_extent in zip(grids, reversed(available_mode_extents))
    ]
    shared_mask = np.logical_and.reduce(shared_components) & (kmag > 0.0)
    unavailable_mask = (kmag > 0.0) & ~shared_mask
    return {
        "dominant_k_over_k0": float(kmag[dominant] / k0),
        "dominant_wavelength": float(2.0 * math.pi / kmag[dominant]),
        "available_mode_extents": list(available_mode_extents),
        "large_box_only_unavailable_mode_power_fraction": float(
            np.sum(total[unavailable_mask]) / power
        ),
        "shared_small_box_mode_power_fraction": float(np.sum(total[shared_mask]) / power),
        "shared_unavailable_power_partition_residual": float(
            (np.sum(total[unavailable_mask]) + np.sum(total[shared_mask])) / power - 1.0
        ),
    }


def _morphology(fields: Mapping[str, np.ndarray]) -> dict[str, object]:
    density = fields["dens"]
    bperp2 = fields["bcc2"] ** 2 + fields["bcc3"] ** 2
    density_mean = float(np.mean(density))
    bperp2_mean = float(np.mean(bperp2))
    density_std = float(np.std(density))
    bperp2_std = float(np.std(bperp2))
    pearson = (
        None
        if density_std == 0.0 or bperp2_std == 0.0
        else float(np.corrcoef(density.ravel(), bperp2.ravel())[0, 1])
    )
    cavity = density < 0.5 * density_mean
    filament = bperp2 > 4.0 * max(bperp2_mean, np.finfo(float).tiny)
    return {
        "density_Bperp2_pearson": pearson,
        "density_quantiles": np.quantile(density, [0.01, 0.1, 0.5, 0.9, 0.99]).tolist(),
        "Bperp2_quantiles": np.quantile(bperp2, [0.01, 0.1, 0.5, 0.9, 0.99]).tolist(),
        "cavity_volume_fraction_rho_below_half_mean": float(np.mean(cavity)),
        "filament_volume_fraction_Bperp2_above_four_mean": float(np.mean(filament)),
        "joint_cavity_filament_volume_fraction": float(np.mean(cavity & filament)),
        "numeric_morphology_acceptance_thresholds": None,
    }


def _energy_transfer(fields: Mapping[str, np.ndarray]) -> dict[str, float]:
    speed2 = fields["velx"] ** 2 + fields["vely"] ** 2 + fields["velz"] ** 2
    b2 = fields["bcc1"] ** 2 + fields["bcc2"] ** 2 + fields["bcc3"] ** 2
    bperp2 = fields["bcc2"] ** 2 + fields["bcc3"] ** 2
    force2 = (
        fields["prtcl_dpxdt"] ** 2
        + fields["prtcl_dpydt"] ** 2
        + fields["prtcl_dpzdt"] ** 2
    )
    return {
        "volume_mean_gas_kinetic_energy": float(np.mean(0.5 * fields["dens"] * speed2)),
        "volume_mean_gas_internal_energy_proxy": float(np.mean(fields["eint"])),
        "volume_mean_total_magnetic_energy": float(np.mean(0.5 * b2)),
        "volume_mean_transverse_magnetic_energy": float(np.mean(0.5 * bperp2)),
        "volume_mean_particle_energy_transfer_rate": float(np.mean(fields["prtcl_dedt"])),
        "volume_mean_parallel_particle_force": float(np.mean(fields["prtcl_dpxdt"])),
        "volume_rms_particle_force": float(np.sqrt(np.mean(force2))),
        "volume_mean_particle_E_dot_B_diagnostic": float(np.mean(fields["prtcl_ebdot"])),
    }


def _initial_rho_jx_spatial_noise(fields: Mapping[str, np.ndarray]) -> dict[str, float]:
    rho = fields["prtcl_rho"]
    jx = fields["prtcl_jx"]
    rho_mean = float(np.mean(rho))
    jx_mean = float(np.mean(jx))
    rho_scale = max(abs(rho_mean), np.finfo(float).tiny)
    jx_scale = max(abs(jx_mean), np.finfo(float).tiny)
    transverse_rms = float(
        np.sqrt(np.mean(fields["prtcl_jy"] ** 2 + fields["prtcl_jz"] ** 2))
    )
    return {
        "deposited_charge_density_mean": rho_mean,
        "deposited_charge_density_coefficient_of_variation": float(
            np.std(rho) / rho_scale
        ),
        "deposited_charge_density_max_absolute_fractional_deviation": float(
            np.max(np.abs(rho - rho_mean)) / rho_scale
        ),
        "jx_over_c_mean": jx_mean,
        "jx_over_c_coefficient_of_variation": float(np.std(jx) / jx_scale),
        "jx_over_c_max_absolute_fractional_deviation": float(
            np.max(np.abs(jx - jx_mean)) / jx_scale
        ),
        "transverse_current_rms_over_abs_parallel_mean": transverse_rms / jx_scale,
    }


def _complex_k0_mode(
    fields: Mapping[str, np.ndarray], x1_faces: np.ndarray, *, k0: float
) -> tuple[complex, complex]:
    x1 = 0.5 * (x1_faces[:-1] + x1_faces[1:])
    phase = np.exp(-1j * k0 * x1)[None, None, :]
    bmode = np.mean((fields["bcc2"] + 1j * fields["bcc3"]) * phase)
    jmode = np.mean((fields["prtcl_jy"] + 1j * fields["prtcl_jz"]) * phase)
    return complex(bmode), complex(jmode)


def _deposited_particle_moment_state(
    fields: Mapping[str, np.ndarray], *, cycle: int, time: float
) -> dict[str, object]:
    rho_q = fields["prtcl_rho"]
    currents = (fields["prtcl_jx"], fields["prtcl_jy"], fields["prtcl_jz"])
    exactly_empty = bool(
        np.all(rho_q == 0.0) and all(np.all(value == 0.0) for value in currents)
    )
    if exactly_empty:
        _require(
            cycle == 0 and time == 0.0,
            "empty deposited particle moments are permitted only at exact cycle-zero startup",
        )
        return {
            "available": False,
            "status": "natural_cycle_zero_pre_deposition_output",
            "positive_charge_density_cell_fraction": 0.0,
        }
    _require(np.all(rho_q >= 0.0), "deposited charge density must be nonnegative")
    positive = rho_q > 0.0
    _require(np.any(positive), "nonempty deposited particle moments lack charge density")
    _require(
        all(np.all(value[~positive] == 0.0) for value in currents),
        "deposited current exists where deposited charge density is zero",
    )
    return {
        "available": True,
        "status": "valid_single_species_full_f_deposited_moments",
        "positive_charge_density_cell_fraction": float(np.mean(positive)),
    }


def _grid_metrics(
    parsed: Mapping[str, object],
    *,
    case: Mapping[str, object],
    available_mode_extents: tuple[float, float, float],
) -> dict[str, object]:
    fields = parsed["fields"]
    dx = tuple(float(np.diff(axis)[0]) for axis in parsed["faces"])
    measured_extents = tuple(
        float(axis[-1] - axis[0]) for axis in parsed["faces"]
    )
    _require(
        all(
            math.isclose(measured, float(expected), rel_tol=1.0e-12, abs_tol=1.0e-14)
            for measured, expected in zip(measured_extents, case["extents"])
        ),
        "snapshot physical extents drifted from matrix row",
    )
    gas_mass = float(np.sum(fields["dens"]))
    gas_bulk_velocity = np.asarray(
        [
            np.sum(fields["dens"] * fields[name]) / gas_mass
            for name in ("velx", "vely", "velz")
        ],
        dtype=np.float64,
    )
    deposited = np.asarray(
        [
            np.mean(fields["prtcl_jx"]),
            np.mean(fields["prtcl_jy"]),
            np.mean(fields["prtcl_jz"]),
        ],
        dtype=np.float64,
    )
    moment_state = _deposited_particle_moment_state(
        fields, cycle=int(parsed["cycle"]), time=float(parsed["time"])
    )
    deposited_charge_density = float(np.mean(fields["prtcl_rho"]))
    deposited_mass_density = (
        deposited_charge_density
        * float(case["species_mass"])
        / float(case["species_charge"])
    )
    bmode, jmode = _complex_k0_mode(fields, parsed["faces"][0], k0=float(case["k0"]))
    return {
        "cycle": parsed["cycle"],
        "time": parsed["time"],
        "dx": list(dx),
        "physical_extents": list(measured_extents),
        "gas_bulk_velocity": gas_bulk_velocity,
        "deposited_grid_j_over_c": deposited,
        "deposited_particle_moment_state": moment_state,
        "deposited_grid_charge_density": deposited_charge_density,
        "deposited_grid_cr_mass_density": deposited_mass_density,
        "mean_Bperp_rms_over_B0": float(
            np.sqrt(np.mean(fields["bcc2"] ** 2 + fields["bcc3"] ** 2))
            / float(case["b_g"])
        ),
        "complex_Bperp_k0": bmode,
        "complex_deposited_Jperp_k0": jmode,
        "initial_rho_jx_spatial_noise": (
            _initial_rho_jx_spatial_noise(fields) if moment_state["available"] else None
        ),
        "energy_transfer": _energy_transfer(fields),
        "spectrum": _spectrum(
            fields,
            dx,
            k0=float(case["k0"]),
            available_mode_extents=available_mode_extents,
        ),
        "morphology": _morphology(fields),
    }


def _complex_record(value: complex) -> dict[str, float]:
    return {"real": float(value.real), "imag": float(value.imag), "amplitude": abs(value), "phase": float(np.angle(value))}


def _current_agreement(
    grid: Sequence[Mapping[str, object]], particle: Sequence[Mapping[str, object]]
) -> dict[str, object]:
    rows = []
    for g, p in zip(grid, particle):
        direct = np.asarray(g["deposited_grid_j_over_c"])
        reconstructed = np.asarray(p["lab_j_over_c"])
        absolute = float(np.linalg.norm(direct - reconstructed))
        scale = max(float(np.linalg.norm(direct)), float(np.linalg.norm(reconstructed)), 1.0)
        rows.append(
            {
                "cycle": g["cycle"],
                "time": g["time"],
                "absolute_l2_error": absolute,
                "relative_l2_error": absolute / scale,
            }
        )
    maximum = max(float(row["relative_l2_error"]) for row in rows)
    return {
        "trace": rows,
        "relative_l2_tolerance": CURRENT_AGREEMENT_RELATIVE_TOLERANCE,
        "maximum_relative_l2_error": maximum,
        "gate_passed": maximum <= CURRENT_AGREEMENT_RELATIVE_TOLERANCE,
        "qualification_effect": "none",
    }


def _evolving_resolution_gate(
    case: Mapping[str, object],
    grid: Sequence[Mapping[str, object]],
    particle: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    applicable = case["branch"] != "high_rigidity_current_retention_candidate"
    if not applicable:
        return {
            "applicable": False,
            "high_rigidity_validity_uses_rg_over_dominant_wavelength_instead": True,
        }
    rows = []
    first_violation = None
    for g, p in zip(grid, particle):
        dx = max(float(value) for value in g["dx"][: int(case["dimension"])])
        characteristic = float(p["gyroradius"]["macro_mass_weighted_median"]) / dx
        stop_required = characteristic < MINIMUM_CHARACTERISTIC_SHELL_RL_OVER_DX
        if stop_required and first_violation is None:
            first_violation = float(g["time"])
        rows.append(
            {
                "cycle": g["cycle"],
                "time": g["time"],
                "characteristic_median_particle_rl_over_dx": characteristic,
                "mean_particle_rl_over_dx": (
                    float(p["gyroradius"]["macro_mass_weighted_mean"]) / dx
                ),
                "minimum_pitch_angle_sensitive_particle_rl_over_dx_report_only": (
                    float(p["gyroradius"]["minimum"]) / dx
                ),
                "maximum_particle_rl_over_dx": float(p["gyroradius"]["maximum"]) / dx,
                "stop_required": stop_required,
            }
        )
    resolution_floor_satisfied = first_violation is None
    return {
        "applicable": True,
        "minimum_characteristic_shell_rl_over_dx": (
            MINIMUM_CHARACTERISTIC_SHELL_RL_OVER_DX
        ),
        "minimum_sampled_pitch_angle_guard_used": False,
        "distribution_percentile_and_PPC_resolution_convergence_required": True,
        "trace": rows,
        "first_postprocessing_violation_time": first_violation,
        "postprocessing_resolution_floor_satisfied": resolution_floor_satisfied,
        "postprocessing_resolution_gate_required": True,
        "postprocessing_resolution_gate_evaluated": True,
        "runtime_stop_controller_installed": False,
        "runtime_resolution_guard_pilot_qualified": False,
        "gate_passed": False,
        "acceptance_eligible": False,
        "fail_closed_reason": (
            "runtime controller removed and postprocessing distribution-aware "
            "convergence and dominant-scale thresholds remain unset"
            if resolution_floor_satisfied
            else "characteristic finite-rigidity resolution floor crossed"
        ),
    }


def _high_rigidity_validity_gate(
    case: Mapping[str, object],
    grid: Sequence[Mapping[str, object]],
    particle: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    if case["branch"] != "high_rigidity_current_retention_candidate":
        return {"applicable": False, "fixed_current_like_label_authorized": False}
    initial_j = float(grid[0]["deposited_grid_j_over_c"][0])
    initial_p = float(particle[0]["particle_momentum"][0])
    initial_cr_energy = float(particle[0]["particle_kinetic_energy"])
    rows = []
    for g, p in zip(grid, particle):
        dominant_wavelength = float(g["spectrum"]["dominant_wavelength"])
        rows.append(
            {
                "time": g["time"],
                "lab_parallel_current_retention": float(g["deposited_grid_j_over_c"][0]) / initial_j,
                "CR_parallel_momentum_fractional_change": (
                    float(p["particle_momentum"][0]) - initial_p
                ) / abs(initial_p),
                "gas_parallel_acceleration_from_initial": float(
                    g["gas_bulk_velocity"][0] - grid[0]["gas_bulk_velocity"][0]
                ),
                "CR_kinetic_energy_fractional_change": (
                    float(p["particle_kinetic_energy"]) - initial_cr_energy
                )
                / abs(initial_cr_energy),
                "relative_parallel_drift": float(
                    p["particle_bulk_velocity"][0] - g["gas_bulk_velocity"][0]
                ),
                "mean_particle_rg_over_dominant_wavelength": float(
                    p["gyroradius"]["macro_mass_weighted_mean"]
                )
                / dominant_wavelength,
            }
        )
    return {
        "applicable": True,
        "trace": rows,
        "required_measurements": [
            "lab_parallel_current_retention",
            "CR_parallel_momentum_fractional_change",
            "CR_kinetic_energy_fractional_change",
            "gas_parallel_acceleration_from_initial",
            "relative_parallel_drift",
            "mean_particle_rg_over_dominant_wavelength",
        ],
        "required_measurements_complete": True,
        "numeric_validity_thresholds": None,
        "threshold_source": "excluded_pilots_and_external_review",
        "gate_passed": None,
        "all_numeric_thresholds_frozen": False,
        "all_required_gates_passed": False,
        "fixed_current_like_label_authorized": False,
    }


def _high_rigidity_saturation_mechanism_discriminants(
    case: Mapping[str, object],
    grid: Sequence[Mapping[str, object]],
    particle: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    if case["branch"] != "high_rigidity_current_retention_candidate":
        return {
            "applicable": False,
            "mechanism_classification_authorized": False,
            "nonlinear_saturation_claim_authorized": False,
        }
    initial_j = float(grid[0]["deposited_grid_j_over_c"][0])
    initial_p = float(particle[0]["particle_momentum"][0])
    initial_cr_energy = float(particle[0]["particle_kinetic_energy"])
    trace = []
    for g, p in zip(grid, particle):
        trace.append(
            {
                "time": g["time"],
                "Bperp_rms_over_B0": g["mean_Bperp_rms_over_B0"],
                "dominant_k_over_k0": g["spectrum"]["dominant_k_over_k0"],
                "large_box_only_unavailable_mode_power_fraction": g["spectrum"][
                    "large_box_only_unavailable_mode_power_fraction"
                ],
                "density_Bperp2_pearson": g["morphology"]["density_Bperp2_pearson"],
                "cavity_volume_fraction": g["morphology"][
                    "cavity_volume_fraction_rho_below_half_mean"
                ],
                "filament_volume_fraction": g["morphology"][
                    "filament_volume_fraction_Bperp2_above_four_mean"
                ],
                "lab_parallel_current_retention": float(
                    g["deposited_grid_j_over_c"][0]
                )
                / initial_j,
                "CR_parallel_momentum_fractional_change": (
                    float(p["particle_momentum"][0]) - initial_p
                )
                / abs(initial_p),
                "CR_kinetic_energy_fractional_change": (
                    float(p["particle_kinetic_energy"]) - initial_cr_energy
                )
                / abs(initial_cr_energy),
                "gas_parallel_acceleration_from_initial": float(
                    g["gas_bulk_velocity"][0] - grid[0]["gas_bulk_velocity"][0]
                ),
                "relative_parallel_drift": float(
                    p["particle_bulk_velocity"][0] - g["gas_bulk_velocity"][0]
                ),
                "mean_particle_rg_over_dominant_wavelength": float(
                    p["gyroradius"]["macro_mass_weighted_mean"]
                )
                / float(g["spectrum"]["dominant_wavelength"]),
                "energy_transfer": dict(g["energy_transfer"]),
            }
        )
    return {
        "applicable": True,
        "case_role": case["role"],
        "high_sampling_mode": case["high_sampling_mode"],
        "position_noise_reference_case_id": case["high_position_noise_reference_case_id"],
        "first_valid_deposited_moment_cycle": grid[0]["cycle"],
        "first_valid_deposited_moment_time": grid[0]["time"],
        "first_snapshot_rho_jx_spatial_noise": dict(
            grid[0]["initial_rho_jx_spatial_noise"]
        ),
        "first_snapshot_noise_is_first_valid_deposited_moment": True,
        "trace": trace,
        "required_discriminants_complete": True,
        "required_convergence_axes": [
            "PPC",
            "resolution",
            "timestep",
            "cell_centered_vs_stochastic_position_sampling",
            "multiple_stochastic_particle_seeds",
        ],
        "convergence_controls_runtime_complete": False,
        "numeric_mechanism_thresholds": None,
        "mechanism_classification": None,
        "mechanism_classification_authorized": False,
        "nonlinear_saturation_claim_authorized": False,
    }


def _initial_rho_jx_noise_pair_gate(
    case: Mapping[str, object], grid: Sequence[Mapping[str, object]]
) -> dict[str, object]:
    if case["branch"] == "high_rigidity_current_retention_candidate":
        return {
            "applicable": False,
            "gate_passed": None,
            "physical_fiducial_eligible": False,
        }
    return {
        "applicable": True,
        "finite_sampling_mode": case["finite_sampling_mode"],
        "packet_grouping_contract": case["finite_packet_grouping_contract"],
        "packet_grouping_against_actual_initializer_proved": case[
            "finite_packet_grouping_against_actual_initializer_proved"
        ],
        "density_parallel_current_quiet_by_construction": case[
            "finite_density_parallel_current_quiet_by_construction"
        ],
        "first_valid_deposited_moment_cycle": grid[0]["cycle"],
        "first_valid_deposited_moment_time": grid[0]["time"],
        "first_snapshot_measurement": dict(grid[0]["initial_rho_jx_spatial_noise"]),
        "first_snapshot_measurement_is_first_valid_deposited_moment": True,
        "paired_mode_case_id": case["finite_initial_rho_jx_noise_pair_case_id"],
        "paired_mode_comparison_required": True,
        "numeric_acceptance_thresholds": None,
        "threshold_source": "excluded_pilots_and_external_review",
        "gate_passed": None,
        "physical_fiducial_eligible": False,
    }


def _finite_predecessor_measurements(
    case: Mapping[str, object],
    grid: Sequence[Mapping[str, object]],
    current_agreement: Mapping[str, object],
    initial_noise_gate: Mapping[str, object],
    resolution_gate: Mapping[str, object],
) -> dict[str, object] | None:
    if case["branch"] != "finite_rigidity_early_time_predecessor":
        return None
    times = np.asarray([row["time"] for row in grid], dtype=np.float64)
    bmodes = np.asarray([row["complex_Bperp_k0"] for row in grid], dtype=np.complex128)
    jmodes = np.asarray([row["complex_deposited_Jperp_k0"] for row in grid], dtype=np.complex128)
    growth = None
    frequency = None
    if len(times) >= 3 and np.all(np.abs(bmodes) > 0.0):
        growth = float(np.polyfit(times, np.log(np.abs(bmodes)), 1)[0])
        frequency = float(np.polyfit(times, np.unwrap(np.angle(bmodes)), 1)[0])
    return {
        "signed_complex_Bperp_k0_trace": [_complex_record(value) for value in bmodes],
        "signed_complex_deposited_Jperp_k0_trace": [_complex_record(value) for value in jmodes],
        "complex_exponential_growth_and_frequency_fit": {
            "growth_rate": growth,
            "angular_frequency": frequency,
        },
        "polarization_handedness": "measured_from_signed_complex_Bperp_trace",
        "Bperp_Jperp_phase_relation": [
            float(np.angle(jmode / bmode)) if abs(bmode) > 0.0 else None
            for bmode, jmode in zip(bmodes, jmodes)
        ],
        "parallel_current_retention": [
            float(row["deposited_grid_j_over_c"][0] / decks.EXPECTED_J_OVER_C)
            for row in grid
        ],
        "transverse_force_noise": [
            float(abs(row["complex_deposited_Jperp_k0"])) for row in grid
        ],
        "deposited_grid_vs_reconstructed_particle_current_agreement": dict(current_agreement),
        "initial_rho_jx_noise_pair_gate": dict(initial_noise_gate),
        "evolving_rl_over_dx": dict(resolution_gate),
        "fit_window": None,
        "numeric_acceptance_tolerances": None,
        "gate_passed": None,
        "runtime_predecessor_complete": False,
    }


def _distribution_summary(values: np.ndarray) -> dict[str, object]:
    finite = _finite_array(values, label="local applicability diagnostic")
    return {
        "minimum": float(np.min(finite)),
        "maximum": float(np.max(finite)),
        "mean": float(np.mean(finite)),
        "quantiles_01_10_50_90_99": np.quantile(
            finite, [0.01, 0.1, 0.5, 0.9, 0.99]
        ).tolist(),
    }


def _evolving_local_hall_applicability_diagnostics(
    case: Mapping[str, object],
    parsed_grid: Sequence[Mapping[str, object]],
    grid: Sequence[Mapping[str, object]],
    particle: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    trace = []
    species_mass = float(case["species_mass"])
    species_charge = float(case["species_charge"])
    background_qom = float(case["background_q_over_mc_reference"])
    _require(
        species_mass > 0.0 and species_charge > 0.0 and background_qom > 0.0,
        "single-species charge/mass applicability inputs are invalid",
    )
    all_cells_complete = True
    for parsed, g, p in zip(parsed_grid, grid, particle):
        fields = parsed["fields"]
        deposited_charge_density = fields["prtcl_rho"]
        local_b = np.sqrt(
            fields["bcc1"] ** 2 + fields["bcc2"] ** 2 + fields["bcc3"] ** 2
        )
        charge_valid = deposited_charge_density > 0.0
        lambda_valid = charge_valid & (local_b > np.finfo(float).tiny)
        _require(
            np.any(lambda_valid),
            "local Lambda diagnostics have no positive-charge, nonzero-B cells",
        )
        drift_coverage = float(np.mean(charge_valid))
        lambda_coverage = float(np.mean(lambda_valid))
        all_cells_complete = all_cells_complete and lambda_coverage == 1.0
        rho_cr_mass = deposited_charge_density * species_mass / species_charge
        rho_cr_over_rho_gas = rho_cr_mass / fields["dens"]
        raw_charge_ratio = deposited_charge_density / (
            fields["dens"] * background_qom
        )
        local_bai_r = raw_charge_ratio / (1.0 + raw_charge_ratio)
        local_va = local_b / np.sqrt(fields["dens"])
        local_di = 1.0 / (background_qom * np.sqrt(fields["dens"]))
        local_cr_velocity = np.stack(
            [
                np.divide(
                    fields[name],
                    deposited_charge_density,
                    out=np.zeros_like(deposited_charge_density),
                    where=charge_valid,
                )
                for name in ("prtcl_jx", "prtcl_jy", "prtcl_jz")
            ]
        )
        local_gas_velocity = np.stack(
            [fields[name] for name in ("velx", "vely", "velz")]
        )
        local_relative_drift = local_cr_velocity - local_gas_velocity
        local_lambda = np.divide(
            local_bai_r * local_relative_drift[0],
            local_va,
            out=np.zeros_like(local_va),
            where=lambda_valid,
        )
        global_relative_drift = (
            np.asarray(p["particle_bulk_velocity"])
            - np.asarray(g["gas_bulk_velocity"])
        )
        active_dx = np.asarray(g["dx"][: int(case["dimension"])], dtype=np.float64)
        trace.append(
            {
                "cycle": g["cycle"],
                "time": g["time"],
                "deposited_prtcl_rho_semantics": "charge_density",
                "species_mass": species_mass,
                "species_charge": species_charge,
                "background_q_over_mc": background_qom,
                "local_concentration_and_Bai_R_cell_coverage_fraction": 1.0,
                "local_cr_relative_drift_cell_coverage_fraction": drift_coverage,
                "local_background_ion_inertial_length_and_dx_over_di_cell_coverage_fraction": 1.0,
                "local_Lambda_cell_coverage_fraction": lambda_coverage,
                "local_cr_relative_drift_unavailable_only_where": (
                    "deposited_charge_density_is_nonpositive"
                ),
                "local_Lambda_unavailable_only_where": (
                    "deposited_charge_density_is_nonpositive_or_local_B_is_zero"
                ),
                "local_deposited_charge_density": _distribution_summary(
                    deposited_charge_density
                ),
                "local_cr_mass_density": _distribution_summary(rho_cr_mass),
                "local_rho_cr_over_rho_gas": _distribution_summary(
                    rho_cr_over_rho_gas
                ),
                "local_Bai_R": _distribution_summary(local_bai_r),
                "local_cr_minus_gas_relative_drift": {
                    axis: _distribution_summary(
                        local_relative_drift[index][charge_valid]
                    )
                    for index, axis in enumerate(("x1", "x2", "x3"))
                },
                "global_particle_minus_gas_relative_drift": (
                    global_relative_drift.tolist()
                ),
                "global_parallel_relative_drift": float(global_relative_drift[0]),
                "local_alfven_speed": _distribution_summary(local_va),
                "local_background_ion_inertial_length": _distribution_summary(
                    local_di
                ),
                "minimum_active_dx_over_local_di": _distribution_summary(
                    np.min(active_dx) / local_di
                ),
                "maximum_active_dx_over_local_di": _distribution_summary(
                    np.max(active_dx) / local_di
                ),
                "signed_local_Lambda": _distribution_summary(
                    local_lambda[lambda_valid]
                ),
                "absolute_local_Lambda": _distribution_summary(
                    np.abs(local_lambda[lambda_valid])
                ),
                "local_cr_relative_drift_field_available": True,
                "local_background_ion_inertial_length_available": True,
                "local_dx_over_di_available": True,
                "exact_local_Lambda_available": True,
            }
        )
    return {
        "status": decks.NONLINEAR_NO_HALL_APPLICABILITY_STATUS,
        "derivation_scope": "single_species_full_f_positive_charge",
        "derivations": {
            "cr_mass_density": "prtcl_rho_times_species_mass_over_species_charge",
            "cr_bulk_velocity": "prtcl_j_vector_divided_by_prtcl_rho",
            "cr_minus_gas_relative_drift": "cr_bulk_velocity_minus_local_gas_velocity",
            "background_ion_inertial_length": (
                "one_over_background_q_over_mc_times_sqrt_local_gas_mass_density"
            ),
            "Bai_R": (
                "raw_charge_ratio_over_one_plus_raw_charge_ratio_with_raw_charge_ratio_"
                "equal_prtcl_rho_over_gas_density_times_background_q_over_mc"
            ),
            "local_Lambda": "Bai_R_times_local_parallel_relative_drift_over_local_alfven_speed",
        },
        "trace": trace,
        "evolving_local_cr_concentration_and_R_available": True,
        "evolving_global_relative_drift_available": True,
        "local_cr_relative_drift_field_available": True,
        "local_background_ion_inertial_length_available": True,
        "local_dx_over_di_available": True,
        "exact_local_Lambda_available": True,
        "unavailable_local_cr_drift_limited_to_zero_charge_cells": True,
        "unavailable_local_Lambda_limited_to_zero_charge_or_zero_B_cells": True,
        "all_valid_snapshot_cells_have_local_diagnostics": all_cells_complete,
        "evolving_local_hall_applicability_diagnostics_complete": all_cells_complete,
        "registered_evolving_local_hall_diagnostics_bound": False,
        "numeric_no_hall_acceptance_thresholds_frozen": False,
        "gate_passed": False,
        "no_hall_applicability_accepted": False,
        "claim_authorized": False,
    }


def _resolved_scale_and_applicability_gate(
    case: Mapping[str, object],
    evolving_local: Mapping[str, object],
) -> dict[str, object]:
    return {
        "background_q_over_mc": case["background_q_over_mc_reference"],
        "background_ion_inertial_length": case["background_ion_inertial_length"],
        "k0_di": case["k0_background_ion_inertial_length"],
        "minimum_materialized_dx_over_di": case["minimum_active_dx_over_background_di"],
        "charge_density_ratio_R": case["charge_density_ratio_equal_background_qom"],
        "Lambda_Hall": case["hall_parameter_from_current"],
        "Bai_linear_factor": case["bai_hall_linear_factor"],
        "Bai_growth_rate_fractional_shift": case[
            "bai_hall_growth_rate_fractional_shift"
        ],
        "Bai_wavenumber_fractional_shift": case[
            "bai_hall_wavenumber_fractional_shift"
        ],
        "Bai_growth_rate_reduction_factor": case[
            "bai_hall_growth_rate_reduction_factor"
        ],
        "Bai_wavenumber_reduction_factor": case[
            "bai_hall_wavenumber_reduction_factor"
        ],
        "no_subion_cell_scale_envelope_satisfied": case[
            "no_subion_cell_scale_envelope_satisfied"
        ],
        "mhd_resolved_scale_applicability_accepted": False,
        "R_much_less_than_one_applicability_accepted": False,
        "bounded_hall_omission_candidate": True,
        "bounded_hall_omission_review_complete": False,
        "no_hall_applicability_accepted": False,
        "nonlinear_no_hall_applicability_status": case[
            "nonlinear_no_hall_applicability_status"
        ],
        "evolving_local_hall_applicability_diagnostics": dict(evolving_local),
        "numeric_much_larger_and_much_less_acceptance_thresholds": None,
        "gate_passed": False,
        "claim_authorized": False,
    }


def _finite_onset_and_saturation_prerequisite_gate(
    case: Mapping[str, object],
    grid: Sequence[Mapping[str, object]],
    particle: Sequence[Mapping[str, object]],
    resolution_gate: Mapping[str, object],
) -> dict[str, object]:
    if case["branch"] == "high_rigidity_current_retention_candidate":
        return {"applicable": False, "saturation_claim_authorized": False}
    trace = []
    initial_current = float(grid[0]["deposited_grid_j_over_c"][0])
    initial_drift = float(
        particle[0]["particle_bulk_velocity"][0] - grid[0]["gas_bulk_velocity"][0]
    )
    for g, p in zip(grid, particle):
        flux = np.asarray(p["particle_momentum_flux_tensor"])
        trace.append(
            {
                "time": g["time"],
                "Bperp_rms_over_B0": g["mean_Bperp_rms_over_B0"],
                "dominant_wavelength": g["spectrum"]["dominant_wavelength"],
                "dominant_wavelength_over_minimum_active_extent": (
                    float(g["spectrum"]["dominant_wavelength"])
                    / min(float(value) for value in case["extents"][: int(case["dimension"])])
                ),
                "parallel_current_retention": (
                    float(g["deposited_grid_j_over_c"][0]) / initial_current
                ),
                "relative_parallel_drift_change": (
                    float(p["particle_bulk_velocity"][0] - g["gas_bulk_velocity"][0])
                    - initial_drift
                ),
                "measured_anisotropic_momentum_flux": float(
                    flux[0, 0] - 0.5 * (flux[1, 1] + flux[2, 2])
                ),
                "conservation": dict(p["conservation"]),
                "energy_transfer": dict(g["energy_transfer"]),
            }
        )
    return {
        "applicable": True,
        "current_design_role": case["role"],
        "current_row_saturation_candidate": False,
        "target_B_over_B0": None,
        "preregistered_nonlinear_onset_Bperp_rms_over_B0": case[
            "preregistered_nonlinear_onset_Bperp_rms_over_B0"
        ],
        "resolution_design_maximum_sampled_B_over_B0": case[
            "resolution_design_maximum_sampled_B_over_B0"
        ],
        "common_nonlinear_onset_resolution_reachable": case[
            "common_nonlinear_onset_resolution_reachable"
        ],
        "common_nonlinear_onset_reachability_status": case[
            "common_nonlinear_onset_reachability_status"
        ],
        "reachability_inferred_from_resolution_design_envelope": False,
        "coarse_resolution_stop_or_intermittency_limitation": case[
            "coarse_resolution_stop_or_intermittency_limitation"
        ],
        "initial_anisotropic_momentum_flux_predictor_input": case[
            "zacharegkas_saturation_predictor_input"
        ],
        "absolute_total_energy_B_over_B0_bound": case[
            "absolute_total_energy_B_over_B0_bound"
        ],
        "trace": trace,
        "characteristic_resolution_gate": dict(resolution_gate),
        "required_gates": [
            "nonlinear_amplitude",
            "plateau_log_slope",
            "sustained_post_plateau_window",
            "gas_plus_CR_energy_conservation",
            "gas_plus_CR_momentum_conservation",
            "current_retention",
            "relative_drift_change",
            "backreaction_energy_transfer",
            "characteristic_shell_resolution_and_distribution_convergence",
            "dominant_scale_vs_box",
            "matched_shared_spectrum_box_convergence",
        ],
        "numeric_thresholds": None,
        "threshold_source": "excluded_pilots_and_preregistered_review",
        "dominant_scale_runtime_stop_controller_installed": False,
        "all_required_gates_passed": False,
        "saturation_claim_authorized": False,
    }


def _completion_gate(
    source_kind: str, completion_record: Mapping[str, object] | None
) -> dict[str, object]:
    if completion_record is not None:
        try:
            return provenance_boundary.classify_runtime_completion(completion_record)
        except provenance_boundary.ProvenanceBoundaryError as error:
            raise ContractError(str(error)) from error
    _require(
        source_kind == "synthetic_contract_fixture",
        "raw registered bundle requires structured runtime completion status",
    )
    return {
        "run_completion_status": "synthetic_fixture_not_executed",
        "incomplete_rejected": False,
        "evidence_disposition": "not_evidence_synthetic_fixture",
        "saturation_evidence_eligible": False,
        "raw_science_admission_eligible": False,
        "trusted_execution_binding_present": False,
    }


def analyze_snapshots(
    case_id: str,
    snapshots: Sequence[Mapping[str, object]],
    particle_states: Sequence[Mapping[str, object]],
    *,
    source_kind: str,
    provenance: Mapping[str, object] | None = None,
    completion_record: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Analyze matched fixture snapshots without granting authority."""
    _require(source_kind in {"synthetic_contract_fixture", "raw_registered_bundle"}, "source kind is invalid")
    raw_admission: dict[str, object] | None = None
    if source_kind == "raw_registered_bundle":
        _require(provenance is not None, "raw registered bundle requires provenance")
        _require(
            completion_record is not None,
            "raw registered bundle requires structured runtime completion status",
        )
        try:
            raw_admission = provenance_boundary.validate_raw_science_bundle(
                provenance,
                snapshots=snapshots,
                particle_states=particle_states,
                completion_record=completion_record,
            )
        except provenance_boundary.ProvenanceBoundaryError as error:
            raise ContractError(str(error)) from error
    else:
        _require(provenance is None, "synthetic fixtures cannot carry provenance")
    case_map = _case_map()
    case = case_map.get(case_id)
    _require(case is not None, "unknown Q019 case")
    _require(2 <= len(snapshots) == len(particle_states), "snapshot inventory drifted")
    available_mode_extents = _available_mode_extents(case, case_map)
    parsed_grid = [_snapshot(row) for row in snapshots]
    grid = [
        _grid_metrics(
            row, case=case, available_mode_extents=available_mode_extents
        )
        for row in parsed_grid
    ]
    particle = [_particle_state(row, case=case) for row in particle_states]
    cycles = np.asarray([row["cycle"] for row in grid], dtype=np.int64)
    times = np.asarray([row["time"] for row in grid], dtype=np.float64)
    _require(
        cycles[0] == 0 and times[0] == 0.0,
        "chronology must begin at exact cycle zero/time zero",
    )
    _require(np.all(np.diff(cycles) > 0), "grid cycles must increase")
    _require(np.all(np.diff(times) > 0.0), "grid times must increase")
    _require(
        np.array_equal(
            cycles, np.asarray([row["cycle"] for row in particle], dtype=np.int64)
        )
        and np.array_equal(
            times, np.asarray([row["time"] for row in particle], dtype=np.float64)
        ),
        "particle and grid cycle/time chronology drifted",
    )
    valid_moment_indexes = [
        index
        for index, row in enumerate(grid)
        if row["deposited_particle_moment_state"]["available"]
    ]
    _require(
        len(valid_moment_indexes) >= 1,
        "at least one valid deposited-particle-moment snapshot is required",
    )
    first_valid_index = valid_moment_indexes[0]
    _require(
        valid_moment_indexes == list(range(first_valid_index, len(grid))),
        "valid deposited particle moments must form a continuous suffix",
    )
    valid_parsed_grid = [parsed_grid[index] for index in valid_moment_indexes]
    valid_grid = [grid[index] for index in valid_moment_indexes]
    valid_particle = [particle[index] for index in valid_moment_indexes]
    moment_chronology = {
        "strict_grid_and_particle_cycle_time_chronology": True,
        "chronology": [
            {
                "cycle": row["cycle"],
                "time": row["time"],
                "deposited_particle_moment_state": dict(
                    row["deposited_particle_moment_state"]
                ),
            }
            for row in grid
        ],
        "natural_cycle_zero_pre_deposition_output_present": (
            not grid[0]["deposited_particle_moment_state"]["available"]
        ),
        "first_valid_deposited_moment_cycle": valid_grid[0]["cycle"],
        "first_valid_deposited_moment_time": valid_grid[0]["time"],
        "deposited_moment_denominator_and_noise_baseline_cycle": valid_grid[0][
            "cycle"
        ],
        "deposited_moment_denominator_and_noise_baseline_time": valid_grid[0]["time"],
        "valid_deposited_moment_snapshot_count": len(valid_grid),
        "cycle_zero_excluded_only_from_deposited_moment_dependent_diagnostics": True,
    }
    current_agreement = _current_agreement(valid_grid, valid_particle)
    resolution_gate = _evolving_resolution_gate(case, grid, particle)
    initial_noise_gate = _initial_rho_jx_noise_pair_gate(case, valid_grid)
    evolving_local_hall = _evolving_local_hall_applicability_diagnostics(
        case, valid_parsed_grid, valid_grid, valid_particle
    )
    completion_gate = _completion_gate(source_kind, completion_record)
    if raw_admission is not None:
        completion_gate = {
            **completion_gate,
            "evidence_disposition": (
                "rejected_incomplete"
                if completion_gate["incomplete_rejected"]
                else "admitted_registered_raw_analysis"
            ),
            "saturation_evidence_eligible": (
                not completion_gate["incomplete_rejected"]
                and raw_admission["saturation_evidence_eligible"] is True
            ),
            "raw_science_admission_eligible": (
                not completion_gate["incomplete_rejected"]
                and raw_admission["raw_science_admission_eligible"] is True
            ),
        }
    report = {
        "schema_version": 3,
        "record_type": RECORD_TYPE,
        "case_id": case_id,
        "campaign_id": case["campaign_id"],
        "branch": case["branch"],
        "analysis_input_kind": source_kind,
        "authority": {
            "launch_authorized": False,
            "policy_authorized": False,
            "qualification_authorized": False,
            "claim_authorized": False,
            "raw_production_authorized": False,
            "nonlinear_saturation_claim_authorized": False,
        },
        "loading_and_no_hall_accounting": {
            key: case[key]
            for key in (
                "rho_cr_over_rho0",
                "n_cr",
                "species_mass",
                "species_charge",
                "charge_density_ratio_equal_background_qom",
                "hall_parameter_equal_background_qom",
                "hall_parameter_from_current",
                "background_q_over_mc_reference",
                "background_ion_inertial_length",
                "k0_background_ion_inertial_length",
                "minimum_active_dx_over_background_di",
                "bai_hall_linear_factor",
                "bai_hall_growth_rate_fractional_shift",
                "bai_hall_wavenumber_fractional_shift",
                "bai_hall_growth_rate_reduction_factor",
                "bai_hall_wavenumber_reduction_factor",
                "background_q_over_mc_at_lambda_equal_one_reference",
                "hall_order_unity_reference",
                "hall_order_unity_reference_margin",
                "cr_inertia_parameter",
                "cr_momentum_loading_parameter",
                "cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy",
                "energy_loading_regime",
                "energy_loading_gate_status",
                "energy_loading_grid_review_complete",
                "universal_saturation_inference_authorized",
                "feedback_force_parameter",
                "positive_finite_loading_accounting_satisfied",
                "applicability_scope",
                "strong_shock_applicability_authorized",
                "no_hall_applicability_accepted",
                "nonlinear_no_hall_applicability_status",
                "evolving_local_hall_applicability_diagnostics_complete",
                "bounded_hall_omission_candidate",
                "bounded_hall_omission_review_complete",
                "mhd_resolved_scale_applicability_accepted",
                "R_much_less_than_one_applicability_accepted",
            )
        },
        "structured_runtime_completion_gate": completion_gate,
        "actual_timestep_history_gate": {
            "configured_limiter_is_actual_timestep_evidence": False,
            "per_cycle_history_retention_configured": True,
            "registered_history_binding_present": False,
            "actual_timestep_convergence_gate_passed": None,
            "claim_authorized": False,
        },
        "mhd_scale_R_and_hall_applicability_gate": (
            _resolved_scale_and_applicability_gate(case, evolving_local_hall)
        ),
        "evolving_local_hall_applicability_diagnostics": evolving_local_hall,
        "independent_prerequisite_gate": {
            "q043_independent_raw_cycle_one_oracle_id": case[
                "q043_independent_raw_cycle_one_oracle_id"
            ],
            "q043_independent_raw_cycle_one_oracle_bound": raw_admission is not None,
            "q023_independent_linear_predecessor_id": case[
                "q023_independent_linear_predecessor_id"
            ],
            "q023_independent_linear_predecessor_bound": raw_admission is not None,
            "independent_prerequisites_complete": raw_admission is not None,
            "nonlinear_execution_prerequisites_passed": raw_admission is not None,
            "matrix_authorizes_execution": False,
            "gate_passed": raw_admission is not None,
        },
        "energy_loading_gate": {
            "cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy": case[
                "cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy"
            ],
            "regime": case["energy_loading_regime"],
            "gate_status": case["energy_loading_gate_status"],
            "grid_and_conservation_response_review_complete": False,
            "numeric_acceptance_thresholds": None,
            "gate_passed": None,
            "universal_saturation_inference_authorized": False,
        },
        "deposited_grid_vs_reconstructed_particle_current_gate": current_agreement,
        "deposited_particle_moment_chronology": moment_chronology,
        "initial_rho_jx_noise_pair_gate": initial_noise_gate,
        "evolving_finite_rigidity_resolution_stop_gate": resolution_gate,
        "finite_nonlinear_onset_and_saturation_prerequisite_gate": (
            _finite_onset_and_saturation_prerequisite_gate(
                case, valid_grid, valid_particle, resolution_gate
            )
        ),
        "high_rigidity_fixed_current_like_validity_gate": _high_rigidity_validity_gate(
            case, valid_grid, valid_particle
        ),
        "high_rigidity_saturation_mechanism_discriminants": (
            _high_rigidity_saturation_mechanism_discriminants(
                case, valid_grid, valid_particle
            )
        ),
        "finite_rigidity_early_time_physics_predecessor": _finite_predecessor_measurements(
            case, valid_grid, current_agreement, initial_noise_gate, resolution_gate
        ),
        "relative_drift_trace": [
            {
                "time": g["time"],
                "particle_minus_gas_bulk_velocity": (
                    p["particle_bulk_velocity"] - g["gas_bulk_velocity"]
                ).tolist(),
            }
            for g, p in zip(grid, particle)
        ],
        "nested_spectrum_trace": [
            {"time": row["time"], **row["spectrum"]} for row in grid
        ],
        "morphology_trace": [
            {"time": row["time"], **row["morphology"]} for row in grid
        ],
        "grid_trace": [
            {
                **row,
                "gas_bulk_velocity": row["gas_bulk_velocity"].tolist(),
                "deposited_grid_j_over_c": row["deposited_grid_j_over_c"].tolist(),
                "complex_Bperp_k0": _complex_record(row["complex_Bperp_k0"]),
                "complex_deposited_Jperp_k0": _complex_record(
                    row["complex_deposited_Jperp_k0"]
                ),
            }
            for row in grid
        ],
        "particle_state_trace": [
            {
                **row,
                "particle_bulk_velocity": row["particle_bulk_velocity"].tolist(),
                "particle_momentum": row["particle_momentum"].tolist(),
                "particle_momentum_flux_tensor": row[
                    "particle_momentum_flux_tensor"
                ].tolist(),
                "particle_velocity_pressure_tensor": row[
                    "particle_velocity_pressure_tensor"
                ].tolist(),
                "lab_j_over_c": row["lab_j_over_c"].tolist(),
                "gas_frame_j_over_c": row["gas_frame_j_over_c"].tolist(),
            }
            for row in particle
        ],
        "limitations": [
            "raw science analysis disabled pending hardened installed-control-plane adapter",
            "bounded Hall omission, R<<1, and MHD resolved-scale applicability remain review-pending and unaccepted",
            "single-species full-f local drift, d_i, dx/d_i, and Lambda are derived, but registered binding, thresholds, and external review remain absent so nonlinear no-Hall applicability fails closed",
            "ideal-MHD carrier rows are model demonstrations, not strong-shock mappings",
            "energy-loading grid and conservation response are unreviewed; universal saturation inference is prohibited",
            "finite initial rho/Jx noise thresholds and paired-mode acceptance remain unset",
            "the unbounded per-cycle runtime particle-scan resolution controller is removed; the replacement postprocessing resolution gate remains unqualified and production remains prohibited",
            "current 3D rows are nonlinear-onset/resource pilots, not saturation candidates",
            "common nonlinear-onset reachability is pilot-pending and is not inferred from the max-sampled-B resolution envelope; the coarse row may stop early or sample onset intermittently",
            "independent Q043 raw-cycle-one and corrected Q023 linear-predecessor prerequisites are unbound, so the matrix is non-authorizing",
            "fixed-current-like label remains prohibited until measured gate thresholds freeze and pass",
            "finite predecessor fit windows and tolerances remain pilot-derived and unset",
            "no launch qualification claim or publication authority",
        ],
    }
    if case["branch"] == "finite_rigidity_early_time_predecessor":
        q019_finite_rigidity_early_time_physics_predecessor_v2.validate_analysis_report(report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--synthetic-fixture", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    _require(args.synthetic_fixture is not None and args.output is not None, "fixture and output required")
    fixture = json.loads(args.synthetic_fixture.read_text(encoding="utf-8"))
    report = analyze_snapshots(
        fixture["case_id"],
        fixture["snapshots"],
        fixture["particle_states"],
        source_kind="synthetic_contract_fixture",
    )
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
