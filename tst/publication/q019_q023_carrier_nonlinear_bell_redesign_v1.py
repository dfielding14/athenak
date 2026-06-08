#!/usr/bin/env python3
"""Versioned, non-authorizing Q019 high-rigidity pilot redesign.

The superseded finite-rigidity pilot matrix remains immutable chronology.  This
successor uses the corrected Q023 low-gyrofrequency carrier as an explicitly
large-inertia external-current surrogate.  It prepares only excluded pilots;
production horizons and scientific thresholds remain open until those pilots
are measured.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
from typing import Mapping, Sequence

from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as q019


REPO_ROOT = Path(__file__).resolve().parents[2]
CHECKED_IN_DECK_ROOT = (
    REPO_ROOT / "inputs/publication/q019_q023_carrier_nonlinear_bell_redesign_v1"
)
CHECKED_IN_MANIFEST = CHECKED_IN_DECK_ROOT / "deck_manifest.json"
SCHEMA_VERSION = 1
SUCCESSOR_ID = "q019_q023_carrier_nonlinear_bell_redesign_v1"
CAMPAIGN_ID = "Q019-HR-JOVERC-Q023-CARRIER-FIXED-CURRENT-LIKE-NOHALL-V1"
BRANCH = "high_rigidity_q023_carrier_candidate"
QUALIFICATION_EFFECT = (
    "none_excluded_pilot_design_only_no_saturation_or_publication_authority"
)

RHO0 = 1.0
B0 = 1.0
PRESSURE0 = 1.0
U_A = 1.0
K0 = 2.0 * math.pi
WAVELENGTH = 1.0
EPSILON = 0.4
STREAM_SPEED = U_A / EPSILON
EXPECTED_J_OVER_C = 2.0 * B0 * K0
ARTIFICIAL_LIGHT_SPEED = 2500.0
BACKGROUND_Q_OVER_MC = 10000.0
FIDUCIAL_CR_Q_OVER_MC = K0 * 1.0e-6
RIGIDITY_CONTROL_Q_OVER_MC = (K0 * 1.0e-4, K0 * 1.0e-5)
PPC = 1
GAMMA_ADIABATIC = 5.0 / 3.0
EIGENMODE_AMPLITUDE = 1.0e-4
BROADBAND_AMPLITUDE = 2.5e-5
PILOT_NONLINEAR_B_OVER_B0_SAFETY_ENVELOPE = 10.0
PILOT_NODE_HOUR_CAP = 500.0
TASKS_PER_NODE = 8
MAXIMUM_STRUCTURAL_CYCLE_COUNT = 20000
MAXIMUM_PILOT_ROOT_CELLS = 4_194_304

SOURCE_LINEAGE = (
    "hardened_q043_then_registered_q023_joverc_then_q019_q023_carrier_v1"
)
MATRIX_SCOPE = "q023_carrier_excluded_physical_pilot_design_v1"
DOMAIN_TIME_STATUS = (
    "versioned_redesign_after_finite_rigidity_resource_overrun_execution_prohibited"
)
APPLICABILITY_SCOPE = (
    "q023_low_gyrofrequency_large_inertia_external_current_surrogate_"
    "periodic_ideal_mhd_no_finite_rigidity_or_physical_cr_density_claim"
)
ENERGY_LOADING_REGIME = (
    "large_inertial_carrier_reservoir_explicitly_reported_not_physical_cr_density"
)
ENERGY_LOADING_GATE_STATUS = (
    "external_current_surrogate_only_requires_measured_current_momentum_energy_"
    "invariance_and_conservation"
)
NO_HALL_MAPPING_BASIS = (
    "q023_cr_qom_separate_from_similarity_mapped_background_ion_qom_with_exact_"
    "Bai_R_k0di_Lambda_accounting_and_no_hall_claim_withheld"
)
PHYSICAL_PILOT_GATE = (
    "explicit_independent_q043_then_registered_q023_then_versioned_resource_"
    "redesign_then_excluded_q023_carrier_pilots"
)

AUTHORIZATION = {
    "launch_authorized": False,
    "policy_authorized": False,
    "physical_pilot_authorized": False,
    "qualification_authorized": False,
    "nonlinear_saturation_claim_authorized": False,
    "publication_authorized": False,
}


class RedesignError(ValueError):
    """Reject a physically inconsistent or computationally unbounded redesign."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RedesignError(message)


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _root_cell_volume(extents: Sequence[float], nx: Sequence[int]) -> float:
    return math.prod(float(extent) / int(count) for extent, count in zip(extents, nx))


def _base_template() -> dict[str, object]:
    return copy.deepcopy(
        next(
            case
            for case in q019.expected_cases()
            if case["case_id"] == "q019-hr-current-retention-s0"
        )
    )


def _charge_density_ratio(
    rho_cr_over_rho0: float, species_q_over_mc: float
) -> float:
    ratio = (
        rho_cr_over_rho0 * species_q_over_mc / BACKGROUND_Q_OVER_MC
    )
    return ratio / (1.0 + ratio)


def _case(
    *,
    suffix: str,
    stage: int,
    role: str,
    extents: tuple[float, float, float],
    nx: tuple[int, int, int],
    meshblock_nx: tuple[int, int, int],
    terminal_tau: float,
    field_seed: int,
    species_q_over_mc: float = FIDUCIAL_CR_Q_OVER_MC,
    reconstruct: str = "plm",
    rsolver: str = "llf",
    cfl: float = 0.2,
    pic_max_cell_cross: int = 2,
    matched_control_family: str = "q023_carrier_fiducial",
    control_interpretation: str = "excluded_pilot_fiducial",
    box_pair_id: str = "not_applicable",
) -> dict[str, object]:
    _require(stage in {1, 2, 3}, "Q019 Q023-carrier stage is invalid")
    _require(species_q_over_mc > 0.0, "CR q/(mc) must be positive")
    _require(terminal_tau > 0.0, "terminal normalized time must be positive")
    dimension = 3 if nx[2] > 1 else 2
    _require(
        math.isclose(extents[0] / extents[1], 2.0)
        and (dimension == 2 or math.isclose(extents[1], extents[2])),
        "Q019 Q023-carrier box geometry must preserve the 2:1 aspect ratio",
    )
    root_cell_volume = _root_cell_volume(extents, nx)
    rho_cr_over_rho0 = EXPECTED_J_OVER_C / (
        species_q_over_mc * STREAM_SPEED * RHO0
    )
    deposit_qscale = rho_cr_over_rho0 * RHO0 * root_cell_volume / PPC
    omega = species_q_over_mc * B0
    nominal_rg0 = STREAM_SPEED / omega
    k0_rg0 = K0 * nominal_rg0
    active_dx = [extents[index] / nx[index] for index in range(dimension)]
    minimum_dx = min(active_dx)
    maximum_dx = max(active_dx)
    background_gyrofrequency = BACKGROUND_Q_OVER_MC * B0
    background_di = U_A / background_gyrofrequency
    charge_density_ratio = _charge_density_ratio(
        rho_cr_over_rho0, species_q_over_mc
    )
    hall_parameter = charge_density_ratio * STREAM_SPEED / U_A
    hall_from_current = (
        EXPECTED_J_OVER_C
        / (RHO0 * BACKGROUND_Q_OVER_MC * U_A)
        * (1.0 - charge_density_ratio)
    )
    _require(
        math.isclose(hall_parameter, hall_from_current, rel_tol=1.0e-13),
        "Q019 Q023-carrier Hall accounting identity drifted",
    )
    lorentz = 1.0 / math.sqrt(
        1.0 - (STREAM_SPEED / ARTIFICIAL_LIGHT_SPEED) ** 2
    )
    cr_kinetic = (
        rho_cr_over_rho0
        * RHO0
        * ARTIFICIAL_LIGHT_SPEED**2
        * (lorentz - 1.0)
    )
    seed_bounds = (
        1.6 * BROADBAND_AMPLITUDE,
        EIGENMODE_AMPLITUDE + 1.5 * BROADBAND_AMPLITUDE,
        EIGENMODE_AMPLITUDE + 0.7 * BROADBAND_AMPLITUDE,
    )
    magnetic_upper = 0.5 * (
        (B0 + seed_bounds[0]) ** 2
        + seed_bounds[1] ** 2
        + seed_bounds[2] ** 2
    )
    gas_internal = PRESSURE0 / (GAMMA_ADIABATIC - 1.0)
    total_energy = (
        gas_internal
        + 0.5 * EIGENMODE_AMPLITUDE**2
        + magnetic_upper
        + cr_kinetic
    )
    particle_cell_dt = pic_max_cell_cross * minimum_dx / STREAM_SPEED
    particle_gyro_dt = q019.PIC_THETA_MAX / abs(omega)
    case = _base_template()
    case.update(
        {
            "case_id": f"q019-q023-carrier-{suffix}",
            "pilot_stage": stage,
            "branch": BRANCH,
            "campaign_id": CAMPAIGN_ID,
            "role": role,
            "dimension": dimension,
            "extents": list(extents),
            "nx": list(nx),
            "meshblock_nx": list(meshblock_nx),
            "reconstruct": reconstruct,
            "rsolver": rsolver,
            "nghost": 4 if reconstruct == "wenoz" else 2,
            "cfl": cfl,
            "pic_max_cell_cross": pic_max_cell_cross,
            "ppc": PPC,
            "distribution": "center",
            "finite_sampling_mode": "not_applicable",
            "high_sampling_mode": q019.HIGH_CENTERED_SAMPLING_MODE,
            "high_position_noise_reference_case_id": "not_applicable",
            "field_seed": field_seed,
            "particle_seed": 0,
            "seed_topology": "shared_spectrum_only",
            "separated_long_mode_amplitude": 0.0,
            "rho": RHO0,
            "pressure": PRESSURE0,
            "b_g": B0,
            "u_a": U_A,
            "wavelength": WAVELENGTH,
            "k0": K0,
            "expected_j_over_c": EXPECTED_J_OVER_C,
            "guide_parallel_stream_speed": STREAM_SPEED,
            "epsilon": EPSILON,
            "species_charge": species_q_over_mc,
            "species_q_over_mc_matches_background": False,
            "omega": omega,
            "finite_shell_speed": 0.0,
            "characteristic_shell_p_iso_over_m": 0.0,
            "maximum_initial_speed_over_artificial_c": (
                STREAM_SPEED / ARTIFICIAL_LIGHT_SPEED
            ),
            "rms_initial_speed_over_artificial_c": (
                STREAM_SPEED / ARTIFICIAL_LIGHT_SPEED
            ),
            "artificial_light_speed": ARTIFICIAL_LIGHT_SPEED,
            "k0_rg0": k0_rg0,
            "rho_cr_over_rho0": rho_cr_over_rho0,
            "n_cr": rho_cr_over_rho0,
            "charge_density_ratio_equal_background_qom": charge_density_ratio,
            "hall_parameter_equal_background_qom": hall_parameter,
            "hall_parameter_from_current": hall_from_current,
            "background_q_over_mc_reference": BACKGROUND_Q_OVER_MC,
            "background_ion_gyrofrequency": background_gyrofrequency,
            "background_ion_inertial_length": background_di,
            "k0_background_ion_inertial_length": K0 * background_di,
            "hall_parameter_over_twice_k0_di": (
                hall_from_current / (2.0 * K0 * background_di)
            ),
            "hall_order_unity_reference_margin": 1.0 / hall_from_current,
            "background_q_over_mc_at_lambda_equal_one_reference": (
                EXPECTED_J_OVER_C
                / (RHO0 * U_A)
                * (1.0 - charge_density_ratio)
            ),
            "cr_inertia_parameter": rho_cr_over_rho0,
            "cr_momentum_loading_parameter": (
                rho_cr_over_rho0 * STREAM_SPEED / U_A
            ),
            "cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy": (
                cr_kinetic / (0.5 * B0 * B0)
            ),
            "energy_loading_regime": ENERGY_LOADING_REGIME,
            "energy_loading_gate_status": ENERGY_LOADING_GATE_STATUS,
            "applicability_scope": APPLICABILITY_SCOPE,
            "no_hall_mapping_basis": NO_HALL_MAPPING_BASIS,
            "nominal_rg0": nominal_rg0,
            "initial_nominal_rg0_over_max_active_dx": nominal_rg0 / maximum_dx,
            "required_initial_rg0_over_max_active_dx": 0.0,
            "preregistered_nonlinear_onset_Bperp_rms_over_B0": 0.0,
            "resolution_design_maximum_sampled_B_over_B0": 0.0,
            "maximum_sampled_B_over_B0_before_resolution_stop": 0.0,
            "common_nonlinear_onset_resolution_reachable": False,
            "common_nonlinear_onset_reachability_status": "not_applicable",
            "coarse_resolution_stop_or_intermittency_limitation": False,
            "postprocessing_resolution_gate_required": False,
            "runtime_box_edge_monitor_unique_mode_count": (
                q019.box_edge_unique_mode_count(dimension)
            ),
            "runtime_box_edge_monitor_cell_passes_per_sample": (
                q019.box_edge_cell_passes_per_sample(dimension)
            ),
            "runtime_box_edge_monitor_global_reductions_per_sample": (
                q019.box_edge_cell_passes_per_sample(dimension)
            ),
            "bai_hall_linear_factor": 1.0 + 0.25 * hall_from_current**2,
            "bai_hall_growth_rate_fractional_shift": (
                1.0 / math.sqrt(1.0 + 0.25 * hall_from_current**2) - 1.0
            ),
            "bai_hall_wavenumber_fractional_shift": (
                1.0 / (1.0 + 0.25 * hall_from_current**2) - 1.0
            ),
            "bai_hall_growth_rate_reduction_factor": (
                1.0 / math.sqrt(1.0 + 0.25 * hall_from_current**2)
            ),
            "bai_hall_wavenumber_reduction_factor": (
                1.0 / (1.0 + 0.25 * hall_from_current**2)
            ),
            "minimum_active_dx_over_background_di": minimum_dx / background_di,
            "initial_anisotropic_momentum_flux_proxy": (
                rho_cr_over_rho0 * STREAM_SPEED**2
            ),
            "zacharegkas_saturation_predictor_input": (
                rho_cr_over_rho0 * STREAM_SPEED**2
            ),
            "absolute_total_energy_B_over_B0_bound": (
                math.sqrt(2.0 * total_energy) / B0
            ),
            "initial_cr_momentum_flux_tensor": [
                [rho_cr_over_rho0 * STREAM_SPEED**2, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ],
            "initial_cr_kinetic_energy_density_rms_proxy": cr_kinetic,
            "conservative_cr_kinetic_energy_density_upper_bound": cr_kinetic,
            "initial_seed_delta_b_component_bounds": list(seed_bounds),
            "initial_magnetic_energy_density_upper_bound": magnetic_upper,
            "initial_gas_kinetic_energy_density_upper_bound": (
                0.5 * EIGENMODE_AMPLITUDE**2
            ),
            "matched_control_family": matched_control_family,
            "control_interpretation": control_interpretation,
            "box_pair_id": box_pair_id,
            "saturation_candidate": False,
            "spectral_sensitivity_control": False,
            "minimum_active_extent_over_nominal_rg0": (
                min(extents[:dimension]) / nominal_rg0
            ),
            "minimum_active_extent_over_seed_wavelength": (
                min(extents[:dimension]) / WAVELENGTH
            ),
            "initial_particle_cell_crossing_dt_bound": particle_cell_dt,
            "initial_particle_gyro_dt_bound": particle_gyro_dt,
            "initial_particle_timestep_limiter": (
                "cell_crossing"
                if particle_cell_dt < particle_gyro_dt
                else "gyro_angle"
            ),
            "root_cell_volume": root_cell_volume,
            "deposit_qscale": deposit_qscale,
            "terminal_time": terminal_tau / (K0 * U_A),
            "terminal_tau": terminal_tau,
            "cycle_limit": -1,
            "eigenmode_amplitude": EIGENMODE_AMPLITUDE,
            "broadband_amplitude": BROADBAND_AMPLITUDE,
            "resource_classification": (
                "frontier_medium_excluded_3d_pilot"
                if dimension == 3
                else "frontier_small_excluded_2d_pilot"
            ),
            "source_lineage": SOURCE_LINEAGE,
            "matrix_scope": MATRIX_SCOPE,
            "domain_time_status": DOMAIN_TIME_STATUS,
            "physical_pilot_gate": PHYSICAL_PILOT_GATE,
        }
    )
    case["matrix_identity_fingerprint"] = q019.matrix_identity_fingerprint(case)
    return case


def expected_cases() -> tuple[dict[str, object], ...]:
    cases = [
        _case(
            suffix=f"s1-onset-s{seed_index}",
            stage=1,
            role="q023_carrier_linear_onset_pilot",
            extents=(4.0, 2.0, 1.0),
            nx=(128, 64, 1),
            meshblock_nx=(32, 32, 1),
            terminal_tau=12.0,
            field_seed=23050091 + seed_index,
            control_interpretation="excluded_linear_onset_seed_pair",
        )
        for seed_index in range(2)
    ]
    fiducial = {
        "stage": 2,
        "role": "q023_carrier_nonlinear_window_pilot",
        "extents": (8.0, 4.0, 1.0),
        "nx": (256, 128, 1),
        "meshblock_nx": (32, 32, 1),
        "terminal_tau": 30.0,
    }
    for seed_index in range(2):
        cases.append(
            _case(
                suffix=f"s2-window-fiducial-s{seed_index}",
                field_seed=23050091 + seed_index,
                control_interpretation="excluded_nonlinear_window_seed_pair",
                **fiducial,
            )
        )
    controls = (
        ("resolution-coarse", {"nx": (128, 64, 1)}),
        ("resolution-fine", {"nx": (512, 256, 1)}),
        ("cfl-small", {"cfl": 0.1}),
        ("cell-cross-one", {"pic_max_cell_cross": 1}),
        ("riemann-hlld", {"rsolver": "hlld"}),
        ("reconstruct-wenoz", {"reconstruct": "wenoz"}),
    )
    for name, changes in controls:
        arguments = dict(fiducial)
        arguments.update(changes)
        arguments["role"] = "q023_carrier_numerical_control"
        cases.append(
            _case(
                suffix=f"s2-{name}-s0",
                field_seed=23050091,
                control_interpretation=f"excluded_single_axis_{name}_control",
                **arguments,
            )
        )
    for exponent, species_q_over_mc in zip((4, 5), RIGIDITY_CONTROL_Q_OVER_MC):
        cases.append(
            _case(
                suffix=f"s2-qom1em{exponent}-s0",
                field_seed=23050091,
                species_q_over_mc=species_q_over_mc,
                control_interpretation=(
                    "fixed_current_vary_cr_qom_and_inertial_reservoir_"
                    "toward_q023_fixed_current_limit"
                ),
                matched_control_family="q023_carrier_rigidity_ladder",
                **{**fiducial, "role": "q023_carrier_rigidity_control"},
            )
        )
    cases.extend(
        [
            _case(
                suffix="s3-3d-fiducial-s0",
                stage=3,
                role="q023_carrier_3d_window_pilot",
                extents=(8.0, 4.0, 4.0),
                nx=(128, 64, 64),
                meshblock_nx=(32, 32, 32),
                terminal_tau=12.0,
                field_seed=23050091,
                matched_control_family="q023_carrier_3d_short_window",
                control_interpretation="excluded_3d_short_window_pilot",
                box_pair_id="q023-carrier-3d-short-window-s0",
            ),
            _case(
                suffix="s3-3d-large-s0",
                stage=3,
                role="q023_carrier_3d_box_control",
                extents=(16.0, 8.0, 8.0),
                nx=(256, 128, 128),
                meshblock_nx=(32, 32, 32),
                terminal_tau=6.0,
                field_seed=23050091,
                matched_control_family="q023_carrier_3d_short_window",
                control_interpretation="excluded_3d_large_box_cost_and_mode_control",
                box_pair_id="q023-carrier-3d-short-window-s0",
            ),
        ]
    )
    _require(
        len({str(case["case_id"]) for case in cases}) == len(cases),
        "Q019 Q023-carrier case IDs collided",
    )
    return tuple(cases)


def estimated_cycle_count(
    case: Mapping[str, object],
    *,
    nonlinear_b_over_b0: float = PILOT_NONLINEAR_B_OVER_B0_SAFETY_ENVELOPE,
) -> int:
    _require(nonlinear_b_over_b0 >= 1.0, "nonlinear B envelope is invalid")
    dimension = int(case["dimension"])
    extents = [float(value) for value in case["extents"]]
    nx = [int(value) for value in case["nx"]]
    minimum_dx = min(extents[index] / nx[index] for index in range(dimension))
    cfl = float(case["cfl"])
    fast_speed = math.sqrt(
        GAMMA_ADIABATIC * PRESSURE0 / RHO0
        + nonlinear_b_over_b0**2 * B0**2 / RHO0
    )
    mhd_bound = minimum_dx / fast_speed
    particle_cell_bound = (
        int(case["pic_max_cell_cross"]) * minimum_dx / STREAM_SPEED
    )
    particle_gyro_bound = q019.PIC_THETA_MAX / (
        abs(float(case["species_charge"])) * nonlinear_b_over_b0
    )
    timestep = cfl * min(mhd_bound, particle_cell_bound, particle_gyro_bound)
    _require(timestep > 0.0 and math.isfinite(timestep), "timestep bound is invalid")
    return math.ceil(float(case["terminal_time"]) / timestep)


def selected_nodes(case: Mapping[str, object]) -> int:
    blocks = math.prod(
        int(global_nx) // int(block_nx)
        for global_nx, block_nx in zip(case["nx"], case["meshblock_nx"])
    )
    return max(1, math.ceil(blocks / TASKS_PER_NODE))


def feasibility_summary(
    cases: Sequence[Mapping[str, object]] | None = None,
) -> dict[str, object]:
    selected = tuple(expected_cases() if cases is None else cases)
    rows = []
    for case in selected:
        cycles = estimated_cycle_count(case)
        nodes = selected_nodes(case)
        root_cells = math.prod(int(value) for value in case["nx"])
        rows.append(
            {
                "case_id": case["case_id"],
                "stage": case["pilot_stage"],
                "estimated_cycles_at_B_over_B0_10": cycles,
                "root_cells": root_cells,
                "macro_particles": root_cells * int(case["ppc"]),
                "selected_nodes_pre_measurement": nodes,
                "node_cycles": nodes * cycles,
            }
        )
    total_node_cycles = sum(int(row["node_cycles"]) for row in rows)
    maximum_seconds_per_cycle_to_fit_cap = (
        PILOT_NODE_HOUR_CAP * 3600.0 / total_node_cycles
    )
    return {
        "record_type": "q019_q023_carrier_structural_feasibility_v1",
        "nonlinear_B_over_B0_safety_envelope": (
            PILOT_NONLINEAR_B_OVER_B0_SAFETY_ENVELOPE
        ),
        "case_count": len(rows),
        "rows": rows,
        "maximum_case_cycles": max(int(row["estimated_cycles_at_B_over_B0_10"]) for row in rows),
        "maximum_root_cells": max(int(row["root_cells"]) for row in rows),
        "total_node_cycles": total_node_cycles,
        "pilot_node_hour_cap": PILOT_NODE_HOUR_CAP,
        "maximum_seconds_per_cycle_to_fit_cap": (
            maximum_seconds_per_cycle_to_fit_cap
        ),
        "measured_resource_authority": False,
        "execution_authorized": False,
    }


def validate_cases(cases: Sequence[Mapping[str, object]]) -> None:
    _require(len(cases) == 14, "Q019 Q023-carrier pilot inventory drifted")
    _require(
        {int(case["pilot_stage"]) for case in cases} == {1, 2, 3},
        "Q019 Q023-carrier stage inventory drifted",
    )
    for case in cases:
        _require(case["campaign_id"] == CAMPAIGN_ID, "campaign identity drifted")
        _require(case["branch"] == BRANCH, "branch identity drifted")
        _require(case["ppc"] == 1, "cold centered carrier must use one PPC")
        _require(
            math.isclose(
                q019.configured_j_over_c(case),
                EXPECTED_J_OVER_C,
                rel_tol=1.0e-13,
            ),
            f"{case['case_id']}: volume-aware current closure drifted",
        )
        _require(
            case["species_q_over_mc_matches_background"] is False,
            f"{case['case_id']}: CR/background q/(mc) separation drifted",
        )
        text = q019.render_deck(case)
        q019.validate_rendered_deck(case, text)
    feasibility = feasibility_summary(cases)
    _require(
        feasibility["maximum_case_cycles"] <= MAXIMUM_STRUCTURAL_CYCLE_COUNT,
        "Q019 Q023-carrier pilot exceeds the structural cycle ceiling",
    )
    _require(
        feasibility["maximum_root_cells"] <= MAXIMUM_PILOT_ROOT_CELLS,
        "Q019 Q023-carrier pilot exceeds the root-cell ceiling",
    )
    _require(
        feasibility["maximum_seconds_per_cycle_to_fit_cap"] >= 1.0,
        "Q019 Q023-carrier design would require an implausible resource rate",
    )


def build_deck_manifest() -> tuple[dict[str, object], dict[str, str]]:
    cases = expected_cases()
    validate_cases(cases)
    rendered = {
        f"{case['case_id']}.athinput": q019.render_deck(case) for case in cases
    }
    records = []
    for case in cases:
        filename = f"{case['case_id']}.athinput"
        payload = rendered[filename].encode("utf-8")
        records.append(
            {
                **case,
                "path": (
                    "inputs/publication/q019_q023_carrier_nonlinear_bell_"
                    f"redesign_v1/{filename}"
                ),
                "sha256": _sha256(payload),
                "byte_count": len(payload),
            }
        )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q019_q023_carrier_nonlinear_bell_redesign_v1_manifest",
        "successor_id": SUCCESSOR_ID,
        "status": "versioned_redesign_complete_execution_prohibited",
        "supersedes_for_execution": (
            "q019_excluded_physical_window_preregistration_v1_finite_"
            "rigidity_matrix_due_to_measured_resource_overrun"
        ),
        "qualification_effect": QUALIFICATION_EFFECT,
        "authorization": dict(AUTHORIZATION),
        "physics": {
            "rho0": RHO0,
            "B0": B0,
            "P0": PRESSURE0,
            "U_A": U_A,
            "k0": K0,
            "epsilon": EPSILON,
            "v_CR": STREAM_SPEED,
            "gamma0": K0 * U_A * math.sqrt(1.0 - EPSILON**2),
            "Omega0": FIDUCIAL_CR_Q_OVER_MC * B0,
            "gamma0_over_Omega0": (
                K0 * U_A * math.sqrt(1.0 - EPSILON**2)
                / (FIDUCIAL_CR_Q_OVER_MC * B0)
            ),
            "fiducial_k0_r_g0": K0 * STREAM_SPEED / FIDUCIAL_CR_Q_OVER_MC,
            "expected_J_CR_over_c": EXPECTED_J_OVER_C,
            "fiducial_rho_CR_over_rho0": (
                EXPECTED_J_OVER_C / (FIDUCIAL_CR_Q_OVER_MC * STREAM_SPEED)
            ),
            "artificial_light_speed": ARTIFICIAL_LIGHT_SPEED,
            "background_q_over_mc_mapping": BACKGROUND_Q_OVER_MC,
            "external_current_surrogate_not_physical_cr_density": True,
            "exact_locked_current": False,
            "finite_rigidity_claim_authorized": False,
            "no_hall_plasma_claim_authorized": False,
        },
        "structural_feasibility": feasibility_summary(cases),
        "required_next_evidence": [
            "passed_registered_Q043_matrix",
            "passed_registered_Q023_linear_matrix",
            "measured_Frontier_seconds_per_cycle_and_memory",
            "registered_excluded_pilot_policy_and_admission",
            "current_momentum_energy_invariance_threshold_freeze",
            "production_window_and_resource_freeze_before_qualifying_output",
        ],
        "decks": records,
    }
    return manifest, rendered


def materialize_checked_in_decks(*, replace: bool = False) -> dict[str, object]:
    manifest, rendered = build_deck_manifest()
    if CHECKED_IN_DECK_ROOT.exists():
        _require(replace, "checked-in redesign root exists; pass --replace")
        shutil.rmtree(CHECKED_IN_DECK_ROOT)
    CHECKED_IN_DECK_ROOT.mkdir(parents=True)
    for filename, text in rendered.items():
        (CHECKED_IN_DECK_ROOT / filename).write_text(text, encoding="utf-8")
    CHECKED_IN_MANIFEST.write_bytes(_json_bytes(manifest))
    return manifest


def validate_checked_in_decks() -> dict[str, object]:
    manifest, rendered = build_deck_manifest()
    _require(CHECKED_IN_MANIFEST.is_file(), "checked-in redesign manifest is absent")
    _require(
        json.loads(CHECKED_IN_MANIFEST.read_text(encoding="utf-8")) == manifest,
        "checked-in redesign manifest drifted",
    )
    expected = set(rendered) | {"deck_manifest.json"}
    actual = {path.name for path in CHECKED_IN_DECK_ROOT.iterdir() if path.is_file()}
    _require(actual == expected, "checked-in redesign deck inventory drifted")
    for filename, text in rendered.items():
        _require(
            (CHECKED_IN_DECK_ROOT / filename).read_text(encoding="utf-8") == text,
            f"checked-in redesign deck drifted: {filename}",
        )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--materialize-checked-in-decks", action="store_true")
    parser.add_argument("--validate-checked-in-decks", action="store_true")
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()
    _require(
        args.materialize_checked_in_decks != args.validate_checked_in_decks,
        "choose exactly one checked-in redesign operation",
    )
    result = (
        materialize_checked_in_decks(replace=args.replace)
        if args.materialize_checked_in_decks
        else validate_checked_in_decks()
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
