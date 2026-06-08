#!/usr/bin/env python3
"""Materialize the fail-closed physics-first Q019 nonlinear Bell design.

This module defines source-local candidate decks and their exact physics
contracts. It cannot launch, modify policy, qualify evidence, or authorize a
scientific claim. Numerical science thresholds not supplied by the literature
remain explicitly pilot-derived and unset.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
from typing import Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[2]
CHECKED_IN_DECK_ROOT = (
    REPO_ROOT / "inputs/publication/q019_physics_first_nonlinear_bell_successor_v2"
)
CHECKED_IN_MANIFEST = CHECKED_IN_DECK_ROOT / "deck_manifest.json"

SCHEMA_VERSION = 4
PGEN_NAME = "q019_physics_first_nonlinear_bell_successor_v2"
PGEN_BLOCK = PGEN_NAME
HIGH_RIGIDITY_CAMPAIGN = (
    "Q019-HR-SIMILARITY-MAPPED-HALL-OMISSION-CANDIDATE-NONSHOCK-V4"
)
FINITE_RIGIDITY_CAMPAIGN = (
    "Q019-FR-ISOTROPIC-SHELL-ONSET-CANDIDATE-NONSHOCK-V4"
)
FINITE_PREDECESSOR_CAMPAIGN = "Q019-FR-ISOTROPIC-SHELL-EARLY-TIME-PREDECESSOR-V4"
QUALIFICATION_EFFECT = "none_source_local_design_only"
MATRIX_SCOPE = "physics_first_preproduction_design_successor_v4"
DOMAIN_TIME_STATUS = "engineering_candidate_pending_excluded_resource_and_window_pilots"
NUMERIC_TOLERANCE_STATUS = (
    "machine_scale_moment_tolerance_and_pilot_science_thresholds_unset"
)
COMMON_ONSET_REACHABILITY_STATUS = (
    "pilot_pending_not_inferred_from_resolution_design_envelope"
)
NONLINEAR_NO_HALL_APPLICABILITY_STATUS = (
    "fail_closed_pending_registered_local_diagnostics_thresholds_and_external_review"
)
RUNTIME_IDENTITY_CHECKSUM_STATUS = (
    "sha256_cryptographic_immutable_runtime_semantics_bound_by_compiled_registry"
)
RUNTIME_RESOLUTION_STOP_CONTROLLER_STATUS = (
    "not_installed_unbounded_per_cycle_scan_removed_pending_excluded_pilot_"
    "benchmark_and_freeze"
)
DEPOSITED_PRTCL_RHO_SEMANTICS = "single_species_charge_density"
DEPOSITED_CR_MASS_DENSITY_DERIVATION = (
    "prtcl_rho_times_species_mass_over_species_charge"
)
Q043_INDEPENDENT_ORACLE_ID = (
    "Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE_raw_cycle_one"
)
Q023_INDEPENDENT_PREDECESSOR_ID = (
    "Q023-PAPER-BELL-LINEAR-JOVERC_after_passed_Q043_oracle"
)

BELL_2004_DOI = "10.1111/j.1365-2966.2004.08097.x"
RIQUELME_SPITKOVSKY_2009_DOI = "10.1088/0004-637X/694/1/626"
GARGATE_2010_DOI = "10.1088/2041-8205/711/2/L127"
ZACHAREGKAS_2024_DOI = "10.3847/1538-4357/ad3960"
SUN_BAI_2023_ARXIV = "2304.10568v1"
BAI_2015_ARXIV = "1412.1087"
RETAINED_REFERENCE_MAP = (
    "tst/publication/readiness/"
    "q022_xcmp_corrected_nonlinear_bell_reference_map_successor_v2_2026-06-07.json"
)

B_G = 1.0
RHO0 = 1.0
PRESSURE0 = 1.0
U_A = B_G / math.sqrt(RHO0)
K0 = 2.0 * math.pi
WAVELENGTH = 1.0
EXPECTED_J_OVER_C = 2.0 * B_G * K0
ARTIFICIAL_LIGHT_SPEED = 100000000.0
SIMILARITY_SCALED_ION_Q_OVER_MC = 10000.0
BACKGROUND_Q_OVER_MC_MAPPING = SIMILARITY_SCALED_ION_Q_OVER_MC
BACKGROUND_ION_GYROFREQUENCY = BACKGROUND_Q_OVER_MC_MAPPING * B_G
BACKGROUND_ION_INERTIAL_LENGTH = U_A / BACKGROUND_ION_GYROFREQUENCY
K0_BACKGROUND_ION_INERTIAL_LENGTH = K0 * BACKGROUND_ION_INERTIAL_LENGTH
HALL_ORDER_UNITY_REFERENCE = 1.0
FIDUCIAL_EXACT_BAI_R = 0.000003 / (1.0 + 0.000003)
FIDUCIAL_HALL_PARAMETER = (
    EXPECTED_J_OVER_C / (RHO0 * BACKGROUND_Q_OVER_MC_MAPPING * U_A)
    * (1.0 - FIDUCIAL_EXACT_BAI_R)
)

FINITE_K0_RG0_GRID = (4.0, 8.0, 16.0)
FINITE_RHO_CR_GRID = (0.000001, 0.000003, 0.00001)
FINITE_PPC_LADDER = (24, 48, 96)
FINITE_FIDUCIAL_K0_RG0 = 8.0
FINITE_FIDUCIAL_RHO_CR = 0.000003
FINITE_FIDUCIAL_PPC = 48

HIGH_K0_RG0 = 256.0
HIGH_PPC_LADDER = (1, 8, 32)
HIGH_FIDUCIAL_PPC = 8
HIGH_STREAM_SPEED = HIGH_K0_RG0 * SIMILARITY_SCALED_ION_Q_OVER_MC * B_G / K0
HIGH_FIDUCIAL_RHO_CR = EXPECTED_J_OVER_C / (
    SIMILARITY_SCALED_ION_Q_OVER_MC * HIGH_STREAM_SPEED
)
HIGH_CENTERED_SAMPLING_MODE = "cell_centered_cold_beam"
HIGH_STOCHASTIC_SAMPLING_MODE = "seeded_random_cold_beam_noise_control"

FIELD_SEEDS = (19001, 19003, 19007)
PARTICLE_SEEDS = (29001, 29003, 29009)
SEED_PAIRS = tuple(zip(FIELD_SEEDS, PARTICLE_SEEDS))

EIGENMODE_AMPLITUDE = 1.0e-6
BROADBAND_AMPLITUDE = 2.5e-7
SEPARATED_LONG_MODE_AMPLITUDE = 2.5e-7
PREDECESSOR_EIGENMODE_AMPLITUDE = 1.0e-7

MIN_CHARACTERISTIC_SHELL_RL_OVER_DX = 8.0
REQUIRED_INITIAL_RL_OVER_DX = MIN_CHARACTERISTIC_SHELL_RL_OVER_DX
PREREGISTERED_NONLINEAR_ONSET_BPERP_RMS_OVER_B0 = 1.0
PREREGISTERED_ONSET_RESOLUTION_MAXIMUM_SAMPLED_B_OVER_B0 = 2.0
PIC_THETA_MAX = 0.3
ENERGY_LOADING_REGIME = "isotropic_shell_plus_drift_energy_accounted_periodic_system"
ENERGY_LOADING_GATE_STATUS = (
    "unqualified_requires_coupled_response_grid_conservation_and_saturation_review"
)
APPLICABILITY_SCOPE = (
    "similarity_scaled_equal_cr_background_qom_mhd_scale_and_R_explicit_"
    "applicability_review_pending_not_strong_shock_mapping"
)

ENGINEERING_TERMINAL_TIME = 12.0
PREDECESSOR_TERMINAL_TIME = 1.5
FIELD_OUTPUT_DT = 0.1
PARTICLE_OUTPUT_DT = 0.5
RESTART_OUTPUT_DT = 0.5
BOX_EDGE_MONITOR_DT = 0.1
BOX_EDGE_STOP_PPM = -1
BOX_EDGE_MONITOR_SCHEMA = 2
BOX_EDGE_MONITOR_SELECTION_RULE = (
    "physical_wave_number_ball_abs_k_le_sqrt_dimension_times_2pi_over_"
    "shortest_active_extent"
)
BOX_EDGE_MONITOR_GEOMETRY = "frozen_x1_to_transverse_aspect_ratio_2"
BOX_EDGE_MONITOR_DIAGNOSTIC_AUTHORITY = "none_quarantined_excluded_pilot_only"

FINITE_2D_EXTENTS = (8.0, 4.0, 1.0)
FINITE_2D_NX = (1024, 512, 1)
FINITE_3D_SMALL_EXTENTS = (16.0, 8.0, 8.0)
FINITE_3D_SMALL_NX = (256, 128, 128)
FINITE_3D_LARGE_EXTENTS = (32.0, 16.0, 16.0)
FINITE_3D_LARGE_NX = (512, 256, 256)
HIGH_2D_NX = (512, 256, 1)
FINITE_3D_FIDUCIAL_K0_RG0 = 8.0
FINITE_3D_FIDUCIAL_RHO_CR = FINITE_FIDUCIAL_RHO_CR
PIC_MAX_CELL_CROSS_BASELINE = 2
PIC_MAX_CELL_CROSS_CONTROL = 1

FINITE_STRATIFIED_SAMPLING_MODE = "cell_centered_nested_haar_octahedral_packets"
FINITE_QUIET_POSITION_SAMPLING = (
    "one_cell_centered_collocated_nested_haar_octahedral_packets_per_cell"
)
FINITE_NOISE_POSITION_SAMPLING = "seeded_random_unique_positions_noise_sensitivity"
FINITE_QUIET_GYROPHASE_POSITION_COUPLING = (
    "nested_haar_rotated_antipodal_packets_exact_zero_rest_mean_isotropic_second_"
    "moment_and_local_tsc_current_quiet_start"
)
FINITE_NOISE_GYROPHASE_POSITION_COUPLING = (
    "nested_haar_rotated_octahedral_packets_independent_positions"
)


def box_edge_unique_mode_count(dimension: int) -> int:
    _require(dimension in {2, 3}, "box-edge monitor dimension is invalid")
    return 7 if dimension == 2 else 23


def box_edge_reduction_group_count(dimension: int) -> int:
    return math.ceil(box_edge_unique_mode_count(dimension) / 5)


def box_edge_cell_passes_per_sample(dimension: int) -> int:
    return 1 + box_edge_reduction_group_count(dimension)


def box_edge_complex_multiplications_per_cell(dimension: int) -> int:
    _require(dimension in {2, 3}, "box-edge monitor dimension is invalid")
    return 28 if dimension == 2 else 106
FINITE_PACKET_GROUPING_CONTRACT = (
    "ppc_divisible_by_six_nested_packets_bound_to_global_root_cell_identity_and_"
    "within_cell_packet_index_and_field_particle_seeds"
)
FINITE_QUIET_PACKET_GROUPING_CONTRACT = (
    FINITE_PACKET_GROUPING_CONTRACT
    + "_with_collocated_packet_spatial_current_quiet_path_proved"
)
FINITE_NOISE_PACKET_GROUPING_CONTRACT = (
    FINITE_PACKET_GROUPING_CONTRACT
    + "_velocity_sequence_only_independent_positions_not_spatially_current_quiet"
)
FINITE_INITIAL_NOISE_GATE = (
    "required_quiet_isotropic_shell_vs_independent_positions_first_snapshot_rho_jx"
)

PHYSICAL_PILOT_GATE_ORDER = (
    f"independent_prerequisite:{Q043_INDEPENDENT_ORACLE_ID}",
    f"independent_prerequisite:{Q023_INDEPENDENT_PREDECESSOR_ID}",
    "q019_finite_rigidity_early_time_physics_predecessor",
    "q019_source_local_compatibility",
    "excluded_resource_and_science_window_pilots",
    "reviewed_registered_execution_preregistration",
)

REQUIRED_OUTPUTS = (
    ("mhd_w_bcc", "bin", FIELD_OUTPUT_DT),
    ("prtcl_rho", "bin", FIELD_OUTPUT_DT),
    ("prtcl_jx", "bin", FIELD_OUTPUT_DT),
    ("prtcl_jy", "bin", FIELD_OUTPUT_DT),
    ("prtcl_jz", "bin", FIELD_OUTPUT_DT),
    ("prtcl_dedt", "bin", FIELD_OUTPUT_DT),
    ("prtcl_dpxdt", "bin", FIELD_OUTPUT_DT),
    ("prtcl_dpydt", "bin", FIELD_OUTPUT_DT),
    ("prtcl_dpzdt", "bin", FIELD_OUTPUT_DT),
    ("prtcl_ebdot", "bin", FIELD_OUTPUT_DT),
    ("prtcl_all", "pvtk", PARTICLE_OUTPUT_DT),
    ("rst", "rst", RESTART_OUTPUT_DT),
)

# These parameters are deterministically added before UserProblem() when absent
# from the rendered deck. They are part of the exact runtime semantics and hence
# part of the restart-stable matrix SHA-256 identity.
RUNTIME_DEFAULT_SEMANTICS = {
    "mesh_refinement": {
        "refinement": "none",
    },
    "time": {"start_time": "0"},
    "coord": {
        "special_rel": "0",
        "general_rel": "0",
    },
    "mhd": {
        "dfloor": "1.17549e-38",
        "pfloor": "1.17549e-38",
        "tfloor": "1.17549e-38",
        "sfloor": "1.17549e-38",
        "nscalars": "0",
        "const_accel": "0",
        "cooling_dt_factor": "1",
        "t_start_ism_cooling": "0",
        "ism_cooling": "0",
        "cgm_cooling": "0",
        "beam_source": "0",
        "rel_cooling": "0",
        "fofc": "0",
    },
    "particles": {
        "track_displacement": "0",
        "pic_ion_neutral_collision_rate": "0",
        "pic_deltaf_f0": "",
        "pic_deltaf_p0": "1",
        "pic_deltaf_kappa": "1.25",
        "pic_deltaf_drift_x1": "0",
        "pic_deltaf_drift_x2": "0",
        "pic_deltaf_drift_x3": "0",
        "pic_deltaf_aniso_x1": "1",
        "pic_deltaf_aniso_x2": "1",
        "pic_deltaf_aniso_x3": "1",
        "pic_deltaf_background_rho": "0",
        "pic_deltaf_background_jx": "0",
        "pic_deltaf_background_jy": "0",
        "pic_deltaf_background_jz": "0",
        "pic_deltaf_adapt_mode": "off",
        "pic_deltaf_adapt_interval": "0",
        "pic_q017_sync_kernel_timers": "0",
        "pic_expansion_law": "linear",
        "pic_expansion_rate_x1": "0",
        "pic_expansion_rate_x2": "0",
        "pic_expansion_rate_x3": "0",
        "pic_no_mhd_bx": "0",
        "pic_no_mhd_by": "0",
        "pic_no_mhd_bz": "0",
        "assign_tag": "index_order",
    },
    "problem": {
        "user_srcs": "0",
        "user_hist": "0",
    },
}

ANALYSIS_BINDING_PATHS = (
    "tst/publication/q019_physics_first_nonlinear_bell_successor_v2.py",
    "tst/publication/analyze_q019_physics_first_nonlinear_bell_successor_v2.py",
    "tst/publication/q019_finite_rigidity_early_time_physics_predecessor_v2.py",
    "tst/publication/q019_nonlinear_bell_particle_state.py",
    "tst/publication/q019_particle_state_analysis_bridge_v2.py",
    "tst/publication/q019_registered_raw_reduction_v1.py",
    "tst/publication/q019_hardened_provenance_boundary_v2.py",
)


class ContractError(ValueError):
    """Raised when the source-local Q019 design drifts."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


MUTABLE_BOX_EDGE_MONITOR_STATE_PARAMETERS = (
    "runtime_box_edge_monitor_next_nominal_time",
    "runtime_box_edge_monitor_last_nominal_time",
    "runtime_box_edge_monitor_last_cycle",
    "runtime_box_edge_monitor_completed_slots",
    "runtime_box_edge_monitor_valid_samples",
    "runtime_box_edge_monitor_skipped_slots",
    "runtime_box_edge_monitor_last_status",
    "runtime_box_edge_monitor_status_mask",
    "runtime_box_edge_monitor_last_prior_time",
    "runtime_box_edge_monitor_last_time",
    "runtime_box_edge_monitor_last_power_fraction",
    "runtime_box_edge_monitor_max_power_fraction",
    "runtime_box_edge_monitor_last_fluctuation_mean",
)


def _mutable_runtime_bookkeeping_parameter(block_name: str, name: str) -> bool:
    return (
        block_name.startswith("output") and name in {"file_number", "last_time"}
    ) or (
        block_name == PGEN_BLOCK
        and name in MUTABLE_BOX_EDGE_MONITOR_STATE_PARAMETERS
    )


def deck_semantics_payload(blocks: Mapping[str, Mapping[str, str]]) -> str:
    """Canonicalize immutable inputs for the restart-stable SHA-256 binding."""
    entries = []
    for block_name, parameters in blocks.items():
        if block_name == "comment":
            continue
        for name, value in parameters.items():
            if block_name == PGEN_BLOCK and name == "matrix_identity_fingerprint":
                continue
            if _mutable_runtime_bookkeeping_parameter(block_name, name):
                continue
            entries.append((block_name, name, value))
    for block_name, parameters in RUNTIME_DEFAULT_SEMANTICS.items():
        for name, value in parameters.items():
            if block_name not in blocks or name not in blocks[block_name]:
                entries.append((block_name, name, value))
    return "".join(
        f"{block_name}/{name}={value}\n"
        for block_name, name, value in sorted(entries)
    )


def matrix_identity_payload(case: Mapping[str, object]) -> str:
    text = _render_deck(case, matrix_identity_fingerprint_value="excluded")
    return deck_semantics_payload(parse_athinput_text(text))


def matrix_identity_fingerprint(case: Mapping[str, object]) -> str:
    return _sha256_bytes(matrix_identity_payload(case).encode("utf-8"))


def _float_token(value: float) -> str:
    return format(value, ".17g")


def _id_float(value: float) -> str:
    return format(value, "g").replace(".", "p").replace("-", "m").replace("+", "")


def _root_cell_volume(extents: Sequence[float], nx: Sequence[int]) -> float:
    _require(len(extents) == len(nx) == 3, "root-cell geometry must be 3D")
    _require(all(value > 0.0 for value in extents), "domain extents must be positive")
    _require(all(type(value) is int and value > 0 for value in nx), "root nx is invalid")
    return math.prod(extent / count for extent, count in zip(extents, nx))


def _runtime_source_paths() -> tuple[str, ...]:
    """Return the deliberately broad future registered runtime-source closure."""
    paths = ["CMakeLists.txt"]
    paths.extend(
        str(path.relative_to(REPO_ROOT))
        for path in (REPO_ROOT / "src").rglob("*")
        if path.is_file()
    )
    return tuple(sorted(paths))


def runtime_source_closure() -> dict[str, object]:
    bindings = [
        {"path": path, "sha256": _sha256_file(REPO_ROOT / path)}
        for path in _runtime_source_paths()
    ]
    return {
        "closure_kind": "all_regular_files_under_src_plus_top_level_CMakeLists",
        "file_count": len(bindings),
        "tree_sha256": _sha256_bytes(_canonical_json_bytes(bindings)),
        "bindings": bindings,
    }


def finite_stream_speed(
    rho_cr_over_rho0: float,
    *,
    expected_j_over_c: float = EXPECTED_J_OVER_C,
    rho0: float = RHO0,
) -> float:
    return expected_j_over_c / (
        rho_cr_over_rho0 * rho0 * SIMILARITY_SCALED_ION_Q_OVER_MC
    )


def finite_shell_speed(
    k0_rg0: float, *, b_g: float = B_G, k0: float = K0
) -> float:
    return k0_rg0 * SIMILARITY_SCALED_ION_Q_OVER_MC * b_g / k0


def background_ion_inertial_length(
    *, rho0: float, b_g: float, u_a: float
) -> float:
    _require(rho0 > 0.0 and b_g > 0.0 and u_a > 0.0, "background mapping is invalid")
    return u_a / (BACKGROUND_Q_OVER_MC_MAPPING * b_g)


def exact_bai_charge_density_ratio(
    *, rho_cr_over_rho0: float, species_q_over_mc: float, background_q_over_mc: float
) -> float:
    cr_to_background_charge_density = (
        rho_cr_over_rho0 * species_q_over_mc / background_q_over_mc
    )
    return cr_to_background_charge_density / (1.0 + cr_to_background_charge_density)


def required_deposit_qscale(
    *, root_cell_volume: float, ppc: int, rho_cr_over_rho0: float, rho0: float = RHO0
) -> float:
    _require(root_cell_volume > 0.0, "root-cell volume must be positive")
    _require(type(ppc) is int and ppc > 0, "PPC must be a positive integer")
    _require(rho_cr_over_rho0 > 0.0, "rho_CR/rho0 must be positive")
    return rho_cr_over_rho0 * rho0 * root_cell_volume / ppc


def configured_j_over_c(case: Mapping[str, object]) -> float:
    return (
        int(case["ppc"])
        * float(case["deposit_qscale"])
        * float(case["species_charge"])
        * float(case["guide_parallel_stream_speed"])
        / float(case["root_cell_volume"])
    )


def configured_rho_cr_over_rho0(case: Mapping[str, object]) -> float:
    return (
        int(case["ppc"])
        * float(case["deposit_qscale"])
        / float(case["root_cell_volume"])
        / float(case["rho"])
    )


def finite_shell_velocity_moments(
    ppc: int, *, guide_parallel_stream_speed: float, shell_speed: float
) -> dict[str, list[list[float]] | list[float]]:
    _require(ppc in FINITE_PPC_LADDER, "finite PPC is outside the design ladder")
    _require(ppc % 6 == 0, "finite PPC must contain complete octahedral packets")
    centered_covariance = [
        [shell_speed * shell_speed / 3.0 if i == j else 0.0 for j in range(3)]
        for i in range(3)
    ]
    second = [
        [
            centered_covariance[i][j]
            + (
                guide_parallel_stream_speed * guide_parallel_stream_speed
                if i == 0 and j == 0
                else 0.0
            )
            for j in range(3)
        ]
        for i in range(3)
    ]
    return {
        "mean": [guide_parallel_stream_speed, 0.0, 0.0],
        "second_moment": second,
        "centered_covariance": centered_covariance,
    }


def _finite_noise_pair_case_id(
    *,
    branch: str,
    k0_rg0: float | None,
    rho_cr_over_rho0: float | None,
    ppc: int | None,
    nx: Sequence[int],
    finite_sampling_mode: str,
) -> str:
    if branch == "hr":
        return "not_applicable"
    if not (
        k0_rg0 == FINITE_FIDUCIAL_K0_RG0
        and rho_cr_over_rho0 == FINITE_FIDUCIAL_RHO_CR
    ):
        return "not_materialized_exact_pair"
    if branch == "fr" and ppc == FINITE_FIDUCIAL_PPC and tuple(nx) == FINITE_2D_NX:
        return (
            "q019-fr-grid-k8-rho3em06-s0"
            if finite_sampling_mode == "independent_position_noise_seeded"
            else "q019-fr-fiducial-noise-seeded-s0"
        )
    if branch == "fr_predecessor" and ppc == max(FINITE_PPC_LADDER) and tuple(nx) == (
        1536,
        768,
        1,
    ):
        return (
            "q019-fr-predecessor-k8-rho3em06-s0"
            if finite_sampling_mode == "independent_position_noise_seeded"
            else "q019-fr-predecessor-fiducial-noise-seeded-s0"
        )
    return "not_materialized_exact_pair"


def _base_case(
    *,
    branch: str,
    suffix: str,
    role: str,
    field_seed: int,
    particle_seed: int,
    dimension: int = 2,
    extents: tuple[float, float, float] = FINITE_2D_EXTENTS,
    nx: tuple[int, int, int] = FINITE_2D_NX,
    meshblock_nx: tuple[int, int, int] = (32, 32, 1),
    k0_rg0: float | None = None,
    rho_cr_over_rho0: float | None = None,
    ppc: int | None = None,
    seed_topology: str = "shared_spectrum_only",
    finite_sampling_mode: str = FINITE_STRATIFIED_SAMPLING_MODE,
    high_sampling_mode: str = HIGH_CENTERED_SAMPLING_MODE,
    reconstruct: str = "plm",
    rsolver: str = "llf",
    cfl: float = 0.2,
    pic_max_cell_cross: int = PIC_MAX_CELL_CROSS_BASELINE,
    terminal_time: float = ENGINEERING_TERMINAL_TIME,
    cycle_limit: int = -1,
    eigenmode_amplitude: float = EIGENMODE_AMPLITUDE,
    broadband_amplitude: float = BROADBAND_AMPLITUDE,
    b_g: float = B_G,
    k0: float = K0,
    finite_shell_speed_override: float | None = None,
    matched_control_family: str = "none",
    control_interpretation: str = "single_case_or_numerical_control",
    box_pair_id: str = "not_applicable",
    saturation_candidate: bool = False,
    spectral_sensitivity_control: bool = False,
) -> dict[str, object]:
    _require(branch in {"hr", "fr", "fr_predecessor"}, "unknown branch")
    _require(
        seed_topology
        in {"shared_spectrum_only", "shared_spectrum_plus_separated_long_modes"},
        "unknown seed topology",
    )
    high = branch == "hr"
    finite = not high
    rho0 = RHO0
    pressure0 = PRESSURE0
    u_a = b_g / math.sqrt(rho0)
    wavelength = 2.0 * math.pi / k0
    expected_j_over_c = 2.0 * b_g * k0
    campaign_id = (
        HIGH_RIGIDITY_CAMPAIGN
        if high
        else FINITE_PREDECESSOR_CAMPAIGN
        if branch == "fr_predecessor"
        else FINITE_RIGIDITY_CAMPAIGN
    )
    branch_name = (
        "high_rigidity_current_retention_candidate"
        if high
        else "finite_rigidity_early_time_predecessor"
        if branch == "fr_predecessor"
        else "finite_rigidity_self_consistent"
    )
    if high:
        k0_rg0 = HIGH_K0_RG0
        _require(
            math.isclose(
                float(rho_cr_over_rho0), HIGH_FIDUCIAL_RHO_CR, rel_tol=1.0e-13
            ),
            "high-rigidity rho_CR/rho0 is not the single closed physical point",
        )
        ppc = 1 if ppc is None else ppc
        _require(ppc in HIGH_PPC_LADDER, "high-rigidity PPC is outside the ladder")
        _require(
            high_sampling_mode
            in {HIGH_CENTERED_SAMPLING_MODE, HIGH_STOCHASTIC_SAMPLING_MODE},
            "high-rigidity sampling mode is invalid",
        )
        if high_sampling_mode == HIGH_CENTERED_SAMPLING_MODE:
            distribution = "center"
            particle_seed = 0
        else:
            _require(
                ppc >= HIGH_FIDUCIAL_PPC and particle_seed > 0,
                "stochastic high-rigidity controls require fiducial PPC and seed",
            )
            distribution = "random"
        finite_sampling_mode = "not_applicable"
        _require(
            math.isclose(b_g, B_G) and math.isclose(k0, K0),
            "high-rigidity controls use the fiducial background mapping",
        )
        stream_speed = HIGH_STREAM_SPEED
        shell_speed = 0.0
        species_charge = SIMILARITY_SCALED_ION_Q_OVER_MC
        rigidity_momentum_per_mass = stream_speed
    else:
        high_sampling_mode = "not_applicable"
        distribution = "random"
        _require(k0_rg0 in FINITE_K0_RG0_GRID, "finite k0*r_g0 is outside the grid")
        _require(
            rho_cr_over_rho0 in FINITE_RHO_CR_GRID,
            "finite rho_CR/rho0 is outside the grid",
        )
        ppc = FINITE_FIDUCIAL_PPC if ppc is None else ppc
        _require(ppc in FINITE_PPC_LADDER, "finite PPC is outside the ladder")
        _require(
            finite_sampling_mode
            in {FINITE_STRATIFIED_SAMPLING_MODE, "independent_position_noise_seeded"},
            "finite sampling mode is invalid",
        )
        species_charge = SIMILARITY_SCALED_ION_Q_OVER_MC
        stream_speed = finite_stream_speed(
            rho_cr_over_rho0, expected_j_over_c=expected_j_over_c, rho0=rho0
        )
        shell_speed = (
            k0_rg0 * species_charge * b_g / k0
            if finite_shell_speed_override is None
            else finite_shell_speed_override
        )
        rigidity_momentum_per_mass = shell_speed
        _require(
            math.isclose(
                k0 * shell_speed / (species_charge * b_g),
                k0_rg0,
                rel_tol=1.0e-13,
            ),
            "finite rigidity mapping drifted",
        )
    omega = species_charge * b_g
    nominal_rg0 = rigidity_momentum_per_mass / omega
    active_dx = [extents[index] / nx[index] for index in range(dimension)]
    initial_rl_over_dx = nominal_rg0 / max(active_dx)
    three_d_onset_family = finite and dimension == 3
    preregistered_nonlinear_onset_bperp_rms_over_b0 = (
        PREREGISTERED_NONLINEAR_ONSET_BPERP_RMS_OVER_B0
        if three_d_onset_family
        else 0.0
    )
    resolution_design_maximum_sampled_b_over_b0 = (
        PREREGISTERED_ONSET_RESOLUTION_MAXIMUM_SAMPLED_B_OVER_B0
        if three_d_onset_family
        else 0.0
    )
    required_initial = 0.0 if high else REQUIRED_INITIAL_RL_OVER_DX
    maximum_sampled_b_over_b0_before_resolution_stop = (
        0.0
        if high
        else initial_rl_over_dx / MIN_CHARACTERISTIC_SHELL_RL_OVER_DX
    )
    common_nonlinear_onset_resolution_reachable = False
    root_cell_volume = _root_cell_volume(extents, nx)
    long_mode_amplitude = (
        SEPARATED_LONG_MODE_AMPLITUDE
        if seed_topology == "shared_spectrum_plus_separated_long_modes"
        else 0.0
    )
    species_mass = 1.0
    n_cr = rho_cr_over_rho0 * rho0 / species_mass
    background_ion_gyrofrequency = BACKGROUND_Q_OVER_MC_MAPPING * b_g
    background_di = background_ion_inertial_length(rho0=rho0, b_g=b_g, u_a=u_a)
    k0_di = k0 * background_di
    charge_density_ratio_equal_background_qom = exact_bai_charge_density_ratio(
        rho_cr_over_rho0=rho_cr_over_rho0,
        species_q_over_mc=species_charge,
        background_q_over_mc=BACKGROUND_Q_OVER_MC_MAPPING,
    )
    hall_parameter_equal_background_qom = (
        charge_density_ratio_equal_background_qom * stream_speed / u_a
    )
    hall_parameter_from_current = expected_j_over_c / (
        rho0 * BACKGROUND_Q_OVER_MC_MAPPING * u_a
    ) * (1.0 - charge_density_ratio_equal_background_qom)
    maximum_speed = stream_speed + shell_speed
    rms_speed = math.sqrt(stream_speed * stream_speed + shell_speed**2)
    total_speed = maximum_speed
    _require(
        maximum_speed < ARTIFICIAL_LIGHT_SPEED,
        "artificial light speed does not bound the shell-plus-drift support",
    )
    rms_gamma = 1.0 / math.sqrt(1.0 - (rms_speed / ARTIFICIAL_LIGHT_SPEED) ** 2)
    maximum_gamma = 1.0 / math.sqrt(
        1.0 - (maximum_speed / ARTIFICIAL_LIGHT_SPEED) ** 2
    )
    cr_inertia_parameter = rho_cr_over_rho0
    cr_momentum_loading_parameter = rho_cr_over_rho0 * stream_speed / u_a
    cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy = (
        rho_cr_over_rho0
        * rho0
        * ARTIFICIAL_LIGHT_SPEED**2
        * (rms_gamma - 1.0)
        / (0.5 * b_g * b_g)
    )
    feedback_force_parameter = (
        expected_j_over_c * b_g / (rho0 * u_a * u_a * k0)
    )
    positive_finite_loading_accounting_satisfied = all(
        math.isfinite(value) and value > 0.0
        for value in (cr_inertia_parameter, cr_momentum_loading_parameter)
    )
    rho_cr = rho_cr_over_rho0 * rho0
    shell_moments = finite_shell_velocity_moments(
        ppc,
        guide_parallel_stream_speed=stream_speed,
        shell_speed=shell_speed,
    ) if finite else None
    isotropic_shell_pressure_proxy = (
        rho_cr * shell_speed * shell_speed / 3.0 if finite else 0.0
    )
    anisotropic_momentum_flux_proxy = rho_cr * stream_speed * stream_speed
    initial_cr_momentum_flux_tensor = [
        [isotropic_shell_pressure_proxy + anisotropic_momentum_flux_proxy, 0.0, 0.0],
        [0.0, isotropic_shell_pressure_proxy, 0.0],
        [0.0, 0.0, isotropic_shell_pressure_proxy],
    ]
    initial_cr_kinetic_energy_density_rms_proxy = (
        rho_cr * ARTIFICIAL_LIGHT_SPEED**2 * (rms_gamma - 1.0)
    )
    conservative_cr_kinetic_energy_density_upper_bound = (
        rho_cr * ARTIFICIAL_LIGHT_SPEED**2 * (maximum_gamma - 1.0)
    )
    gas_internal_energy_density = pressure0 / (5.0 / 3.0 - 1.0)
    seed_delta_b_component_bounds = (
        1.6 * broadband_amplitude + 1.1 * long_mode_amplitude,
        eigenmode_amplitude + 1.5 * broadband_amplitude + 1.7 * long_mode_amplitude,
        eigenmode_amplitude + 0.7 * broadband_amplitude + 0.8 * long_mode_amplitude,
    )
    initial_magnetic_energy_density_upper_bound = 0.5 * (
        (b_g + seed_delta_b_component_bounds[0]) ** 2
        + seed_delta_b_component_bounds[1] ** 2
        + seed_delta_b_component_bounds[2] ** 2
    )
    initial_gas_kinetic_energy_density_upper_bound = 0.5 * eigenmode_amplitude**2
    initial_total_energy_density = (
        gas_internal_energy_density
        + initial_gas_kinetic_energy_density_upper_bound
        + initial_magnetic_energy_density_upper_bound
        + conservative_cr_kinetic_energy_density_upper_bound
    )
    absolute_total_energy_B_over_B0_bound = (
        math.sqrt(2.0 * initial_total_energy_density) / b_g
    )
    minimum_active_extent = min(extents[:dimension])
    minimum_active_extent_over_nominal_rg0 = minimum_active_extent / nominal_rg0
    minimum_active_extent_over_seed_wavelength = minimum_active_extent / wavelength
    minimum_active_dx = min(active_dx)
    initial_particle_cell_crossing_dt_bound = (
        pic_max_cell_cross * minimum_active_dx / total_speed
    )
    initial_particle_gyro_dt_bound = 0.3 / abs(omega)
    finite_noise_pair_case_id = _finite_noise_pair_case_id(
        branch=branch,
        k0_rg0=k0_rg0,
        rho_cr_over_rho0=rho_cr_over_rho0,
        ppc=ppc,
        nx=nx,
        finite_sampling_mode=finite_sampling_mode,
    )
    quiet_isotropic_shell_packet = finite_sampling_mode == FINITE_STRATIFIED_SAMPLING_MODE
    coarse_stop_or_intermittency_limitation = (
        three_d_onset_family
        and initial_rl_over_dx
        < (FINITE_3D_FIDUCIAL_K0_RG0 / K0)
        / (FINITE_3D_SMALL_EXTENTS[0] / FINITE_3D_SMALL_NX[0])
    )
    high_position_noise_reference_case_id = (
        "q019-hr-fiducial-ppc8-centered-s0"
        if high and high_sampling_mode == HIGH_STOCHASTIC_SAMPLING_MODE
        else "not_applicable"
    )
    macro_particles = math.prod(nx) * ppc
    if role == "finite_rigidity_source_local_initializer_regression":
        resource_classification = "source_local_runtime_regression"
    elif dimension == 3 and macro_particles >= 500_000_000:
        resource_classification = "frontier_large_excluded_3d_pilot"
    elif dimension == 3:
        resource_classification = "frontier_medium_excluded_3d_pilot"
    else:
        resource_classification = "source_local_or_frontier_staged_nonproduction"
    case = {
        "case_id": f"q019-{suffix}",
        "branch": branch_name,
        "campaign_id": campaign_id,
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
        "ppc": ppc,
        "nspecies": 1,
        "distribution": distribution,
        "finite_sampling_mode": finite_sampling_mode,
        "high_sampling_mode": high_sampling_mode,
        "high_position_noise_reference_case_id": high_position_noise_reference_case_id,
        "field_seed": field_seed,
        "particle_seed": particle_seed,
        "seed_topology": seed_topology,
        "separated_long_mode_amplitude": long_mode_amplitude,
        "rho": rho0,
        "pressure": pressure0,
        "b_g": b_g,
        "u_a": u_a,
        "wavelength": wavelength,
        "k0": k0,
        "expected_j_over_c": expected_j_over_c,
        "guide_parallel_stream_speed": stream_speed,
        "epsilon": u_a / stream_speed,
        "species_charge": species_charge,
        "species_q_over_mc_matches_background": math.isclose(
            species_charge, BACKGROUND_Q_OVER_MC_MAPPING, rel_tol=0.0, abs_tol=0.0
        ),
        "omega": omega,
        "finite_shell_speed": shell_speed,
        "characteristic_shell_p_iso_over_m": shell_speed,
        "maximum_initial_speed_over_artificial_c": (
            maximum_speed / ARTIFICIAL_LIGHT_SPEED
        ),
        "rms_initial_speed_over_artificial_c": rms_speed / ARTIFICIAL_LIGHT_SPEED,
        "k0_rg0": k0_rg0,
        "rho_cr_over_rho0": rho_cr_over_rho0,
        "species_mass": species_mass,
        "deposited_prtcl_rho_semantics": DEPOSITED_PRTCL_RHO_SEMANTICS,
        "deposited_cr_mass_density_derivation": (
            DEPOSITED_CR_MASS_DENSITY_DERIVATION
        ),
        "species_mass_and_charge_bound_in_immutable_payload": True,
        "n_cr": n_cr,
        "charge_density_ratio_equal_background_qom": (
            charge_density_ratio_equal_background_qom
        ),
        "hall_parameter_equal_background_qom": hall_parameter_equal_background_qom,
        "hall_parameter_from_current": hall_parameter_from_current,
        "background_q_over_mc_reference": BACKGROUND_Q_OVER_MC_MAPPING,
        "background_ion_gyrofrequency": background_ion_gyrofrequency,
        "background_ion_inertial_length": background_di,
        "k0_background_ion_inertial_length": k0_di,
        "hall_parameter_over_twice_k0_di": hall_parameter_from_current / (2.0 * k0_di),
        "hall_order_unity_reference_margin": (
            HALL_ORDER_UNITY_REFERENCE / hall_parameter_from_current
        ),
        "no_hall_mapping_basis": (
            "equal_similarity_scaled_cr_background_qom_with_exact_Bai_R_k0di_Lambda_"
            "resolved_scale_and_signed_Bai_linear_reduction"
        ),
        "hall_order_unity_reference": HALL_ORDER_UNITY_REFERENCE,
        "background_q_over_mc_at_lambda_equal_one_reference": (
            expected_j_over_c
            / (rho0 * u_a)
            * (1.0 - charge_density_ratio_equal_background_qom)
        ),
        "cr_inertia_parameter": cr_inertia_parameter,
        "cr_momentum_loading_parameter": cr_momentum_loading_parameter,
        "cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy": (
            cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy
        ),
        "feedback_force_parameter": feedback_force_parameter,
        "positive_finite_loading_accounting_satisfied": (
            positive_finite_loading_accounting_satisfied
        ),
        "energy_loading_regime": ENERGY_LOADING_REGIME,
        "energy_loading_gate_status": ENERGY_LOADING_GATE_STATUS,
        "energy_loading_grid_review_complete": False,
        "universal_saturation_inference_authorized": False,
        "applicability_scope": APPLICABILITY_SCOPE,
        "strong_shock_applicability_authorized": False,
        "nominal_rg0": nominal_rg0,
        "initial_nominal_rg0_over_max_active_dx": initial_rl_over_dx,
        "required_initial_rg0_over_max_active_dx": required_initial,
        "target_nonlinear_B_over_B0": None,
        "preregistered_nonlinear_onset_Bperp_rms_over_B0": (
            preregistered_nonlinear_onset_bperp_rms_over_b0
        ),
        "resolution_design_maximum_sampled_B_over_B0": (
            resolution_design_maximum_sampled_b_over_b0
        ),
        "maximum_sampled_B_over_B0_before_resolution_stop": (
            maximum_sampled_b_over_b0_before_resolution_stop
        ),
        "common_nonlinear_onset_resolution_reachable": (
            common_nonlinear_onset_resolution_reachable
        ),
        "common_nonlinear_onset_reachability_status": (
            COMMON_ONSET_REACHABILITY_STATUS
            if three_d_onset_family
            else "not_applicable"
        ),
        "coarse_resolution_stop_or_intermittency_limitation": (
            coarse_stop_or_intermittency_limitation
        ),
        "minimum_characteristic_shell_rl_over_dx": (
            MIN_CHARACTERISTIC_SHELL_RL_OVER_DX
        ),
        "resolution_envelope_satisfied_by_design": high
        or initial_rl_over_dx >= required_initial,
        "runtime_identity_checksum_status": RUNTIME_IDENTITY_CHECKSUM_STATUS,
        "external_sha256_execution_receipt_required": True,
        "external_sha256_execution_receipt_bound": False,
        "runtime_resolution_stop_controller_installed": False,
        "runtime_resolution_stop_controller_status": (
            RUNTIME_RESOLUTION_STOP_CONTROLLER_STATUS
        ),
        "runtime_resolution_guard_pilot_qualified": False,
        "postprocessing_resolution_gate_required": finite,
        "runtime_box_edge_monitor_installed": True,
        "runtime_box_edge_monitor_enabled": False,
        "runtime_box_edge_stop_boundary_frozen": False,
        "runtime_box_edge_stop_armed": False,
        "runtime_box_edge_monitor_passive": True,
        "runtime_box_edge_monitor_selection_rule": BOX_EDGE_MONITOR_SELECTION_RULE,
        "runtime_box_edge_monitor_geometry": BOX_EDGE_MONITOR_GEOMETRY,
        "runtime_box_edge_monitor_diagnostic_authority": (
            BOX_EDGE_MONITOR_DIAGNOSTIC_AUTHORITY
        ),
        "runtime_box_edge_monitor_unique_mode_count": box_edge_unique_mode_count(
            dimension
        ),
        "runtime_box_edge_monitor_cell_passes_per_sample": (
            box_edge_cell_passes_per_sample(dimension)
        ),
        "runtime_box_edge_monitor_global_reductions_per_sample": (
            box_edge_cell_passes_per_sample(dimension)
        ),
        "runtime_box_edge_monitor_dt": BOX_EDGE_MONITOR_DT,
        "runtime_box_edge_stop_ppm": BOX_EDGE_STOP_PPM,
        "runtime_dominant_scale_stop_controller_installed": False,
        "raw_production_authorized": False,
        "nonlinear_saturation_claim_authorized": False,
        "no_hall_applicability_accepted": False,
        "nonlinear_no_hall_applicability_status": (
            NONLINEAR_NO_HALL_APPLICABILITY_STATUS
        ),
        "evolving_local_hall_applicability_diagnostics_complete": False,
        "bounded_hall_omission_candidate": True,
        "bounded_hall_omission_review_complete": False,
        "bai_hall_linear_factor": 1.0 + 0.25 * hall_parameter_from_current**2,
        "bai_hall_growth_rate_fractional_shift": (
            1.0 / math.sqrt(1.0 + 0.25 * hall_parameter_from_current**2) - 1.0
        ),
        "bai_hall_wavenumber_fractional_shift": (
            1.0 / (1.0 + 0.25 * hall_parameter_from_current**2) - 1.0
        ),
        "bai_hall_growth_rate_reduction_factor": (
            1.0 / math.sqrt(1.0 + 0.25 * hall_parameter_from_current**2)
        ),
        "bai_hall_wavenumber_reduction_factor": (
            1.0 / (1.0 + 0.25 * hall_parameter_from_current**2)
        ),
        "mhd_resolved_scale_applicability_accepted": False,
        "R_much_less_than_one_applicability_accepted": False,
        "minimum_active_dx_over_background_di": minimum_active_dx / background_di,
        "no_subion_cell_scale_envelope_satisfied": minimum_active_dx > background_di,
        "initial_shell_velocity_moments": shell_moments,
        "isotropic_shell_pressure_proxy": isotropic_shell_pressure_proxy,
        "initial_anisotropic_momentum_flux_proxy": anisotropic_momentum_flux_proxy,
        "zacharegkas_saturation_predictor_input": anisotropic_momentum_flux_proxy,
        "absolute_total_energy_B_over_B0_bound": absolute_total_energy_B_over_B0_bound,
        "saturation_plateau_thresholds": None,
        "initial_cr_momentum_flux_tensor": initial_cr_momentum_flux_tensor,
        "initial_cr_kinetic_energy_density_rms_proxy": (
            initial_cr_kinetic_energy_density_rms_proxy
        ),
        "conservative_cr_kinetic_energy_density_upper_bound": (
            conservative_cr_kinetic_energy_density_upper_bound
        ),
        "initial_seed_delta_b_component_bounds": list(seed_delta_b_component_bounds),
        "initial_magnetic_energy_density_upper_bound": (
            initial_magnetic_energy_density_upper_bound
        ),
        "initial_gas_kinetic_energy_density_upper_bound": (
            initial_gas_kinetic_energy_density_upper_bound
        ),
        "matched_control_family": matched_control_family,
        "control_interpretation": control_interpretation,
        "isolated_variable_claim_authorized": False,
        "box_pair_id": box_pair_id,
        "saturation_candidate": saturation_candidate,
        "spectral_sensitivity_control": spectral_sensitivity_control,
        "minimum_active_extent_over_nominal_rg0": (
            minimum_active_extent_over_nominal_rg0
        ),
        "minimum_active_extent_over_seed_wavelength": (
            minimum_active_extent_over_seed_wavelength
        ),
        "initial_particle_cell_crossing_dt_bound": (
            initial_particle_cell_crossing_dt_bound
        ),
        "initial_particle_gyro_dt_bound": initial_particle_gyro_dt_bound,
        "pic_theta_max": PIC_THETA_MAX,
        "initial_particle_timestep_limiter": (
            "cell_crossing"
            if initial_particle_cell_crossing_dt_bound
            < initial_particle_gyro_dt_bound
            else "gyro_angle"
        ),
        "finite_packet_grouping_contract": (
            "not_applicable"
            if high
            else FINITE_QUIET_PACKET_GROUPING_CONTRACT
            if quiet_isotropic_shell_packet
            else FINITE_NOISE_PACKET_GROUPING_CONTRACT
        ),
        "finite_packet_grouping_against_actual_initializer_proved": (
            quiet_isotropic_shell_packet
        ),
        "finite_packet_independent_spatial_anchors_per_cell_expected": (
            0 if high else (1 if quiet_isotropic_shell_packet else ppc)
        ),
        "finite_density_parallel_current_quiet_by_construction": (
            quiet_isotropic_shell_packet
        ),
        "finite_density_parallel_current_quiet_measured_and_accepted": False,
        "finite_initial_rho_jx_noise_pair_gate": (
            "not_applicable" if high else FINITE_INITIAL_NOISE_GATE
        ),
        "finite_initial_rho_jx_noise_pair_case_id": finite_noise_pair_case_id,
        "root_cell_volume": root_cell_volume,
        "deposit_qscale": required_deposit_qscale(
            root_cell_volume=root_cell_volume,
            ppc=ppc,
            rho_cr_over_rho0=rho_cr_over_rho0,
            rho0=rho0,
        ),
        "terminal_time": terminal_time,
        "cycle_limit": cycle_limit,
        "eigenmode_amplitude": eigenmode_amplitude,
        "broadband_amplitude": broadband_amplitude,
        "resource_classification": resource_classification,
        "q043_independent_raw_cycle_one_oracle_id": Q043_INDEPENDENT_ORACLE_ID,
        "q043_independent_raw_cycle_one_oracle_bound": False,
        "q023_independent_linear_predecessor_id": Q023_INDEPENDENT_PREDECESSOR_ID,
        "q023_independent_linear_predecessor_bound": False,
        "independent_prerequisites_complete": False,
        "nonlinear_execution_prerequisites_passed": False,
    }
    case["matrix_identity_fingerprint"] = matrix_identity_fingerprint(case)
    return case


def _seed_pair(index: int) -> tuple[int, int]:
    return SEED_PAIRS[index]


def expected_cases() -> tuple[dict[str, object], ...]:
    """Return the exact non-authorizing physics-first candidate matrix."""
    cases: list[dict[str, object]] = []
    for seed_index, (field_seed, particle_seed) in enumerate(SEED_PAIRS):
        cases.append(
            _base_case(
                branch="hr",
                suffix=f"hr-current-retention-s{seed_index}",
                role="high_rigidity_current_retention_seed_ensemble",
                field_seed=field_seed,
                particle_seed=particle_seed,
                nx=HIGH_2D_NX,
                rho_cr_over_rho0=HIGH_FIDUCIAL_RHO_CR,
                control_interpretation=(
                    "single_closed_high_rigidity_current_retention_point_not_loading_scan"
                ),
            )
        )
    field_seed, particle_seed = _seed_pair(0)
    for name, changes in (
        ("ppc8-centered", {"ppc": HIGH_FIDUCIAL_PPC}),
        ("ppc32-centered", {"ppc": max(HIGH_PPC_LADDER)}),
        ("resolution-coarse", {"ppc": HIGH_FIDUCIAL_PPC, "nx": (256, 128, 1)}),
        ("resolution-fine", {"ppc": HIGH_FIDUCIAL_PPC, "nx": (1024, 512, 1)}),
        (
            "particle-step-small",
            {
                "ppc": HIGH_FIDUCIAL_PPC,
                "pic_max_cell_cross": PIC_MAX_CELL_CROSS_CONTROL,
            },
        ),
    ):
        control_changes = {"nx": HIGH_2D_NX, **changes}
        cases.append(
            _base_case(
                branch="hr",
                suffix=f"hr-fiducial-{name}-s0",
                role="high_rigidity_staged_convergence_control",
                field_seed=field_seed,
                particle_seed=particle_seed,
                rho_cr_over_rho0=HIGH_FIDUCIAL_RHO_CR,
                control_interpretation=(
                    "matched_high_rigidity_numerical_or_position_sampling_control"
                ),
                **control_changes,
            )
        )
    for seed_index, particle_seed in enumerate(PARTICLE_SEEDS):
        cases.append(
            _base_case(
                branch="hr",
                suffix=f"hr-fiducial-stochastic-position-s{seed_index}",
                role="high_rigidity_staged_convergence_control",
                field_seed=FIELD_SEEDS[0],
                particle_seed=particle_seed,
                rho_cr_over_rho0=HIGH_FIDUCIAL_RHO_CR,
                ppc=HIGH_FIDUCIAL_PPC,
                nx=HIGH_2D_NX,
                high_sampling_mode=HIGH_STOCHASTIC_SAMPLING_MODE,
                control_interpretation=(
                    "matched_high_rigidity_numerical_or_position_sampling_control"
                ),
            )
        )
    for k0_rg0 in FINITE_K0_RG0_GRID:
        for rho_cr in FINITE_RHO_CR_GRID:
            for seed_index, (field_seed, particle_seed) in enumerate(SEED_PAIRS):
                cases.append(
                    _base_case(
                        branch="fr",
                        suffix=(
                            f"fr-grid-k{_id_float(k0_rg0)}-"
                            f"rho{_id_float(rho_cr)}-s{seed_index}"
                        ),
                        role="finite_rigidity_density_drift_coupled_response_ensemble",
                        field_seed=field_seed,
                        particle_seed=particle_seed,
                        k0_rg0=k0_rg0,
                        rho_cr_over_rho0=rho_cr,
                        control_interpretation=(
                            "coupled_rigidity_density_drift_momentum_flux_and_kinetic_"
                            "loading_response_not_an_isolated_variable_scan"
                        ),
                    )
                )
    fixed_shell_speed = finite_shell_speed(FINITE_3D_FIDUCIAL_K0_RG0)
    field_seed, particle_seed = _seed_pair(0)
    for target_k0_rg0 in FINITE_K0_RG0_GRID:
        b_g = math.sqrt(FINITE_3D_FIDUCIAL_K0_RG0 / target_k0_rg0)
        k0 = K0 / b_g
        wavelength = 2.0 * math.pi / k0
        cases.append(
            _base_case(
                branch="fr",
                suffix=f"fr-rigidity-isolation-k{_id_float(target_k0_rg0)}-s0",
                role="finite_rigidity_isolation_fixed_cr_distribution_control",
                field_seed=field_seed,
                particle_seed=particle_seed,
                extents=(8.0 * wavelength, 4.0 * wavelength, wavelength),
                nx=FINITE_2D_NX,
                k0_rg0=target_k0_rg0,
                rho_cr_over_rho0=FINITE_FIDUCIAL_RHO_CR,
                b_g=b_g,
                k0=k0,
                finite_shell_speed_override=fixed_shell_speed,
                matched_control_family="finite_rigidity_isolation_fixed_cr_distribution",
                control_interpretation=(
                    "fixed_full_cr_distribution_current_momentum_flux_and_kinetic_"
                    "loading_while_background_B0_k0_and_UA_change_not_universal_pure_"
                    "rigidity_scaling"
                ),
            )
        )
    field_seed, particle_seed = _seed_pair(0)
    controls = (
        ("ppc24", {"ppc": 24}),
        ("ppc96", {"ppc": 96}),
        ("resolution-coarse", {"nx": (768, 384, 1)}),
        ("resolution-fine", {"nx": (1536, 768, 1)}),
        (
            "particle-step-small",
            {"pic_max_cell_cross": PIC_MAX_CELL_CROSS_CONTROL},
        ),
        ("noise-seeded", {"finite_sampling_mode": "independent_position_noise_seeded"}),
        ("riemann-hlld", {"rsolver": "hlld"}),
        ("reconstruct-wenoz", {"reconstruct": "wenoz"}),
    )
    for name, changes in controls:
        cases.append(
            _base_case(
                branch="fr",
                suffix=f"fr-fiducial-{name}-s0",
                role="finite_rigidity_numerical_or_noise_control",
                field_seed=field_seed,
                particle_seed=particle_seed,
                k0_rg0=FINITE_FIDUCIAL_K0_RG0,
                rho_cr_over_rho0=FINITE_FIDUCIAL_RHO_CR,
                control_interpretation="matched_finite_2d_numerical_or_noise_control",
                **changes,
            )
        )
    for ppc in (min(FINITE_PPC_LADDER), max(FINITE_PPC_LADDER)):
        cases.append(
            _base_case(
                branch="fr",
                suffix=f"fr-runtime-initializer-ppc{ppc}-s0",
                role="finite_rigidity_source_local_initializer_regression",
                field_seed=field_seed,
                particle_seed=particle_seed,
                dimension=2,
                extents=(2.0, 1.0, 1.0),
                nx=(32, 16, 1),
                meshblock_nx=(16, 16, 1),
                k0_rg0=FINITE_FIDUCIAL_K0_RG0,
                rho_cr_over_rho0=FINITE_FIDUCIAL_RHO_CR,
                ppc=ppc,
                terminal_time=1.0e-4,
                cycle_limit=1,
                matched_control_family="q019_source_local_initializer_regression",
                control_interpretation=(
                    "source_local_nested_haar_initializer_runtime_regression_not_science"
                ),
            )
        )
    for seed_index, (field_seed, particle_seed) in enumerate(SEED_PAIRS):
        for size, extents, nx in (
            ("small", FINITE_3D_SMALL_EXTENTS, FINITE_3D_SMALL_NX),
            ("large", FINITE_3D_LARGE_EXTENTS, FINITE_3D_LARGE_NX),
        ):
            cases.append(
                _base_case(
                    branch="fr",
                    suffix=f"fr-3d-onset-{size}-s{seed_index}",
                    role=f"finite_rigidity_3d_nonlinear_onset_box_{size}",
                    field_seed=field_seed,
                    particle_seed=particle_seed,
                    dimension=3,
                    extents=extents,
                    nx=nx,
                    meshblock_nx=(16, 16, 16),
                    k0_rg0=FINITE_3D_FIDUCIAL_K0_RG0,
                    rho_cr_over_rho0=FINITE_3D_FIDUCIAL_RHO_CR,
                    ppc=min(FINITE_PPC_LADDER),
                    seed_topology="shared_spectrum_only",
                    matched_control_family="finite_3d_nonlinear_onset",
                    box_pair_id=f"finite-3d-onset-s{seed_index}",
                    saturation_candidate=False,
                    control_interpretation=(
                        "matched_finite_3d_onset_resource_and_box_pilot_not_saturation_"
                        "candidate"
                    ),
                )
            )
    field_seed, particle_seed = _seed_pair(0)
    for name, changes in (
        ("ppc48", {"ppc": 48}),
        ("resolution-coarse", {"nx": (224, 112, 112)}),
        ("resolution-fine", {"nx": (512, 256, 256)}),
        (
            "particle-step-small",
            {"pic_max_cell_cross": PIC_MAX_CELL_CROSS_CONTROL},
        ),
    ):
        control_changes = {
            "nx": FINITE_3D_SMALL_NX,
            "ppc": min(FINITE_PPC_LADDER),
            **changes,
        }
        cases.append(
            _base_case(
                branch="fr",
                suffix=f"fr-3d-onset-small-{name}-s0",
                role="finite_rigidity_3d_onset_matched_convergence_control",
                field_seed=field_seed,
                particle_seed=particle_seed,
                dimension=3,
                extents=FINITE_3D_SMALL_EXTENTS,
                meshblock_nx=(16, 16, 16),
                k0_rg0=FINITE_3D_FIDUCIAL_K0_RG0,
                rho_cr_over_rho0=FINITE_3D_FIDUCIAL_RHO_CR,
                matched_control_family="finite_3d_nonlinear_onset",
                box_pair_id="finite-3d-onset-s0",
                control_interpretation="matched_finite_3d_numerical_control",
                **control_changes,
            )
        )
    cases.append(
        _base_case(
            branch="fr",
            suffix="fr-3d-onset-large-long-mode-sensitivity-s0",
            role="finite_rigidity_3d_spectral_sensitivity_control",
            field_seed=field_seed,
            particle_seed=particle_seed,
            dimension=3,
            extents=FINITE_3D_LARGE_EXTENTS,
            nx=FINITE_3D_LARGE_NX,
            meshblock_nx=(16, 16, 16),
            k0_rg0=FINITE_3D_FIDUCIAL_K0_RG0,
            rho_cr_over_rho0=FINITE_3D_FIDUCIAL_RHO_CR,
            ppc=min(FINITE_PPC_LADDER),
            seed_topology="shared_spectrum_plus_separated_long_modes",
            matched_control_family="finite_3d_spectral_sensitivity",
            box_pair_id="finite-3d-onset-s0",
            spectral_sensitivity_control=True,
            control_interpretation=(
                "separate_long_mode_spectral_sensitivity_not_box_convergence"
            ),
        )
    )
    for k0_rg0 in FINITE_K0_RG0_GRID:
        for rho_cr in FINITE_RHO_CR_GRID:
            cases.append(
                _base_case(
                    branch="fr_predecessor",
                    suffix=(
                        f"fr-predecessor-k{_id_float(k0_rg0)}-"
                        f"rho{_id_float(rho_cr)}-s0"
                    ),
                    role="finite_rigidity_early_time_physics_reference_grid",
                    field_seed=field_seed,
                    particle_seed=particle_seed,
                    k0_rg0=k0_rg0,
                    rho_cr_over_rho0=rho_cr,
                    ppc=max(FINITE_PPC_LADDER),
                    nx=(1536, 768, 1),
                    terminal_time=PREDECESSOR_TERMINAL_TIME,
                    eigenmode_amplitude=PREDECESSOR_EIGENMODE_AMPLITUDE,
                    broadband_amplitude=0.0,
                    control_interpretation=(
                        "finite_early_time_reference_grid_not_nonlinear_saturation_evidence"
                    ),
                )
            )
    for name, changes in (
        ("ppc24", {"ppc": 24}),
        ("ppc48", {"ppc": 48, "nx": (1536, 768, 1)}),
        ("resolution-coarse", {"nx": (768, 384, 1)}),
        ("resolution-fiducial", {"nx": (1024, 512, 1)}),
        (
            "particle-step-small",
            {"pic_max_cell_cross": PIC_MAX_CELL_CROSS_CONTROL},
        ),
        ("noise-seeded", {"finite_sampling_mode": "independent_position_noise_seeded"}),
    ):
        predecessor_changes = {
            "ppc": max(FINITE_PPC_LADDER),
            "nx": (1536, 768, 1),
            **changes,
        }
        cases.append(
            _base_case(
                branch="fr_predecessor",
                suffix=f"fr-predecessor-fiducial-{name}-s0",
                role="finite_rigidity_early_time_convergence_control",
                field_seed=field_seed,
                particle_seed=particle_seed,
                k0_rg0=FINITE_FIDUCIAL_K0_RG0,
                rho_cr_over_rho0=FINITE_FIDUCIAL_RHO_CR,
                terminal_time=PREDECESSOR_TERMINAL_TIME,
                eigenmode_amplitude=PREDECESSOR_EIGENMODE_AMPLITUDE,
                broadband_amplitude=0.0,
                control_interpretation="matched_finite_early_time_numerical_or_noise_control",
                **predecessor_changes,
            )
        )
    _validate_case_matrix(cases)
    return tuple(cases)


def _validate_case_matrix(cases: Sequence[Mapping[str, object]]) -> None:
    _require(len(cases) == 77, "candidate matrix size drifted")
    _require(len({case["case_id"] for case in cases}) == len(cases), "case IDs collided")
    for case in cases:
        nx = tuple(int(value) for value in case["nx"])
        mb = tuple(int(value) for value in case["meshblock_nx"])
        _require(all(n % block == 0 for n, block in zip(nx, mb)), "decomposition drifted")
        _require(
            case["matrix_identity_fingerprint"] == matrix_identity_fingerprint(case),
            "runtime matrix identity SHA-256 drifted",
        )
        _require(
            math.isclose(
                configured_j_over_c(case),
                float(case["expected_j_over_c"]),
                rel_tol=1.0e-13,
            ),
            "deposited J_CR/c closure drifted",
        )
        _require(
            math.isclose(
                configured_rho_cr_over_rho0(case),
                float(case["rho_cr_over_rho0"]),
                rel_tol=1.0e-13,
            ),
            "rho_CR/rho0 closure drifted",
        )
        _require(
            not case["no_hall_applicability_accepted"]
            and case["bounded_hall_omission_candidate"]
            and not case["bounded_hall_omission_review_complete"]
            and not case["mhd_resolved_scale_applicability_accepted"]
            and not case["R_much_less_than_one_applicability_accepted"]
            and not case["evolving_local_hall_applicability_diagnostics_complete"]
            and case["nonlinear_no_hall_applicability_status"]
            == NONLINEAR_NO_HALL_APPLICABILITY_STATUS,
            "applicability review authority leaked",
        )
        _require(
            case["species_q_over_mc_matches_background"]
            and math.isclose(
                float(case["species_charge"]),
                SIMILARITY_SCALED_ION_Q_OVER_MC,
                rel_tol=0.0,
                abs_tol=0.0,
            )
            and math.isclose(
                float(case["background_q_over_mc_reference"]),
                SIMILARITY_SCALED_ION_Q_OVER_MC,
                rel_tol=0.0,
                abs_tol=0.0,
            ),
            "CR/background similarity-scaled ion q/(mc) equality drifted",
        )
        _require(
            not case["strong_shock_applicability_authorized"]
            and not case["energy_loading_grid_review_complete"]
            and not case["universal_saturation_inference_authorized"]
            and not case["runtime_resolution_stop_controller_installed"]
            and case["runtime_resolution_stop_controller_status"]
            == RUNTIME_RESOLUTION_STOP_CONTROLLER_STATUS
            and not case["runtime_resolution_guard_pilot_qualified"]
            and bool(case["postprocessing_resolution_gate_required"])
            == (case["branch"] != "high_rigidity_current_retention_candidate")
            and case["runtime_identity_checksum_status"]
            == RUNTIME_IDENTITY_CHECKSUM_STATUS
            and case["deposited_prtcl_rho_semantics"]
            == DEPOSITED_PRTCL_RHO_SEMANTICS
            and case["deposited_cr_mass_density_derivation"]
            == DEPOSITED_CR_MASS_DENSITY_DERIVATION
            and case["species_mass_and_charge_bound_in_immutable_payload"]
            and case["external_sha256_execution_receipt_required"]
            and not case["external_sha256_execution_receipt_bound"]
            and case["runtime_box_edge_monitor_installed"]
            and not case["runtime_box_edge_monitor_enabled"]
            and not case["runtime_box_edge_stop_boundary_frozen"]
            and not case["runtime_box_edge_stop_armed"]
            and case["runtime_box_edge_monitor_passive"]
            and case["runtime_box_edge_monitor_selection_rule"]
            == BOX_EDGE_MONITOR_SELECTION_RULE
            and case["runtime_box_edge_monitor_geometry"] == BOX_EDGE_MONITOR_GEOMETRY
            and case["runtime_box_edge_monitor_diagnostic_authority"]
            == BOX_EDGE_MONITOR_DIAGNOSTIC_AUTHORITY
            and int(case["runtime_box_edge_monitor_unique_mode_count"])
            == box_edge_unique_mode_count(int(case["dimension"]))
            and int(case["runtime_box_edge_monitor_cell_passes_per_sample"])
            == box_edge_cell_passes_per_sample(int(case["dimension"]))
            and int(case["runtime_box_edge_monitor_global_reductions_per_sample"])
            == box_edge_cell_passes_per_sample(int(case["dimension"]))
            and math.isclose(
                float(case["runtime_box_edge_monitor_dt"]),
                BOX_EDGE_MONITOR_DT,
                rel_tol=0.0,
                abs_tol=0.0,
            )
            and int(case["runtime_box_edge_stop_ppm"]) == BOX_EDGE_STOP_PPM
            and not case["runtime_dominant_scale_stop_controller_installed"]
            and not case["raw_production_authorized"]
            and not case["nonlinear_saturation_claim_authorized"]
            and case["q043_independent_raw_cycle_one_oracle_id"]
            == Q043_INDEPENDENT_ORACLE_ID
            and not case["q043_independent_raw_cycle_one_oracle_bound"]
            and case["q023_independent_linear_predecessor_id"]
            == Q023_INDEPENDENT_PREDECESSOR_ID
            and not case["q023_independent_linear_predecessor_bound"]
            and not case["independent_prerequisites_complete"]
            and not case["nonlinear_execution_prerequisites_passed"],
            "physics or production authority leaked",
        )
        _require(
            case["positive_finite_loading_accounting_satisfied"],
            "CR loading accounting is not finite and positive",
        )
        _require(
            math.isclose(
                float(case["hall_parameter_equal_background_qom"]),
                float(case["hall_parameter_from_current"]),
                rel_tol=1.0e-13,
            )
            and math.isclose(
                float(case["hall_parameter_from_current"]),
                2.0
                * float(case["k0_background_ion_inertial_length"])
                * (1.0 - float(case["charge_density_ratio_equal_background_qom"])),
                rel_tol=1.0e-13,
            )
            and math.isclose(
                float(case["background_ion_inertial_length"]),
                float(case["u_a"])
                / (
                    BACKGROUND_Q_OVER_MC_MAPPING
                    * float(case["b_g"])
                ),
                rel_tol=1.0e-13,
            ),
            "explicit d_i/k0/Lambda_Hall mapping drifted",
        )
        raw_charge_ratio = (
            float(case["rho_cr_over_rho0"])
            * float(case["species_charge"])
            / float(case["background_q_over_mc_reference"])
        )
        _require(
            math.isclose(
                float(case["charge_density_ratio_equal_background_qom"]),
                raw_charge_ratio / (1.0 + raw_charge_ratio),
                rel_tol=1.0e-13,
            )
            and float(case["bai_hall_growth_rate_fractional_shift"]) < 0.0
            and float(case["bai_hall_wavenumber_fractional_shift"]) < 0.0
            and 0.0 < float(case["bai_hall_growth_rate_reduction_factor"]) < 1.0
            and 0.0 < float(case["bai_hall_wavenumber_reduction_factor"]) < 1.0,
            "exact Bai R or signed reduction semantics drifted",
        )
        _require(
            float(case["hall_parameter_from_current"])
            < HALL_ORDER_UNITY_REFERENCE,
            "mapped Hall-omission candidate crossed the order-unity reference",
        )
        _require(
            float(case["charge_density_ratio_equal_background_qom"]) <= 1.0e-5
            and float(case["guide_parallel_stream_speed"]) > float(case["u_a"])
            and float(case["minimum_active_dx_over_background_di"]) > 1.0
            and case["no_subion_cell_scale_envelope_satisfied"],
            "explicit R, super-Alfvenic drift, or resolved-scale design envelope failed",
        )
        if case["branch"] != "high_rigidity_current_retention_candidate":
            _require(case["resolution_envelope_satisfied_by_design"], "resolution envelope failed")
            quiet = case["finite_sampling_mode"] == FINITE_STRATIFIED_SAMPLING_MODE
            _require(
                bool(case["finite_packet_grouping_against_actual_initializer_proved"])
                == quiet
                and int(case["finite_packet_independent_spatial_anchors_per_cell_expected"])
                == (
                    1
                    if quiet
                    else int(case["ppc"])
                )
                and case["finite_packet_grouping_contract"]
                == (
                    FINITE_QUIET_PACKET_GROUPING_CONTRACT
                    if quiet
                    else FINITE_NOISE_PACKET_GROUPING_CONTRACT
                ),
                "finite packet grouping contract drifted",
            )
            moments = finite_shell_velocity_moments(
                int(case["ppc"]),
                guide_parallel_stream_speed=float(case["guide_parallel_stream_speed"]),
                shell_speed=float(case["finite_shell_speed"]),
            )
            _require(
                moments == case["initial_shell_velocity_moments"]
                and all(
                    math.isclose(
                        float(moments["centered_covariance"][i][j]),
                        float(case["finite_shell_speed"]) ** 2 / 3.0
                        if i == j
                        else 0.0,
                        rel_tol=1.0e-13,
                        abs_tol=1.0e-10,
                    )
                    for i in range(3)
                    for j in range(3)
                ),
                "finite isotropic-shell moment contract drifted",
            )
    grid = [
        case
        for case in cases
        if case["role"] == "finite_rigidity_density_drift_coupled_response_ensemble"
    ]
    _require(len(grid) == 27, "finite physics grid inventory drifted")
    _require(
        {
            (case["k0_rg0"], case["rho_cr_over_rho0"])
            for case in grid
        }
        == {
            (k0_rg0, rho_cr)
            for k0_rg0 in FINITE_K0_RG0_GRID
            for rho_cr in FINITE_RHO_CR_GRID
        },
        "finite physics grid drifted",
    )
    _require(
        {
            (case["field_seed"], case["particle_seed"])
            for case in grid
        }
        == set(SEED_PAIRS),
        "finite deterministic seed ensemble drifted",
    )
    _require(
        all(
            not case["isolated_variable_claim_authorized"]
            and "not_an_isolated_variable_scan" in str(case["control_interpretation"])
            for case in grid
        ),
        "coupled finite response grid was mislabeled as an isolated scan",
    )
    high_seed_ensemble = [
        case
        for case in cases
        if case["role"] == "high_rigidity_current_retention_seed_ensemble"
    ]
    _require(
        len(high_seed_ensemble) == len(SEED_PAIRS)
        and {case["rho_cr_over_rho0"] for case in high_seed_ensemble}
        == {HIGH_FIDUCIAL_RHO_CR}
        and {case["species_charge"] for case in high_seed_ensemble}
        == {SIMILARITY_SCALED_ION_Q_OVER_MC},
        "high-rigidity branch silently became a loading or q/(mc) scan",
    )
    high_controls = [
        case for case in cases if case["role"] == "high_rigidity_staged_convergence_control"
    ]
    _require(len(high_controls) == 8, "high-rigidity convergence controls drifted")
    _require(
        {int(case["ppc"]) for case in high_controls}
        == {HIGH_FIDUCIAL_PPC, max(HIGH_PPC_LADDER)},
        "high-rigidity PPC convergence ladder drifted",
    )
    stochastic = [
        case
        for case in high_controls
        if case["high_sampling_mode"] == HIGH_STOCHASTIC_SAMPLING_MODE
    ]
    _require(
        len(stochastic) == len(SEED_PAIRS)
        and {int(case["particle_seed"]) for case in stochastic} == set(PARTICLE_SEEDS),
        "high-rigidity stochastic position seed ensemble drifted",
    )
    by_id = {str(case["case_id"]): case for case in cases}
    high_centered = by_id["q019-hr-fiducial-ppc8-centered-s0"]
    for case in stochastic:
        _require(
            all(
                case[key] == high_centered[key]
                for key in (
                    "field_seed",
                    "rho_cr_over_rho0",
                    "k0_rg0",
                    "ppc",
                    "nx",
                    "extents",
                    "deposit_qscale",
                    "root_cell_volume",
                )
            )
            and case["high_position_noise_reference_case_id"]
            == high_centered["case_id"],
            "high-rigidity stochastic control is not matched to the centered reference",
        )
    for quiet_id, independent_id in (
        ("q019-fr-grid-k8-rho3em06-s0", "q019-fr-fiducial-noise-seeded-s0"),
        (
            "q019-fr-predecessor-k8-rho3em06-s0",
            "q019-fr-predecessor-fiducial-noise-seeded-s0",
        ),
    ):
        quiet = by_id[quiet_id]
        independent = by_id[independent_id]
        _require(
            all(
                quiet[key] == independent[key]
                for key in (
                    "field_seed",
                    "particle_seed",
                    "rho_cr_over_rho0",
                    "k0_rg0",
                    "ppc",
                    "nx",
                    "extents",
                    "deposit_qscale",
                )
            ),
            "finite quiet-shell/independent-position comparison is not an exact pair",
        )
    onset = [
        case
        for case in cases
        if str(case["role"]).startswith("finite_rigidity_3d_nonlinear_onset_box_")
    ]
    _require(len(onset) == 6, "3D nonlinear-onset seed/box matrix drifted")
    _require(
        all(not case["saturation_candidate"] for case in cases),
        "a current row was mislabeled as a saturation candidate",
    )
    for pair_id in {str(case["box_pair_id"]) for case in onset}:
        pair = [case for case in onset if case["box_pair_id"] == pair_id]
        _require(len(pair) == 2, "3D nonlinear-onset box pair drifted")
        small = next(case for case in pair if str(case["role"]).endswith("_small"))
        large = next(case for case in pair if str(case["role"]).endswith("_large"))
        _require(
            small["seed_topology"] == "shared_spectrum_only"
            and large["seed_topology"] == "shared_spectrum_only"
            and small["field_seed"] == large["field_seed"]
            and small["particle_seed"] == large["particle_seed"]
            and all(2 * int(a) == int(b) for a, b in zip(small["nx"], large["nx"]))
            and all(
                math.isclose(2.0 * float(a), float(b))
                for a, b in zip(small["extents"], large["extents"])
            ),
            "3D onset box convergence changed the shared initial spectrum",
        )
    spectral = [
        case
        for case in cases
        if case["role"] == "finite_rigidity_3d_spectral_sensitivity_control"
    ]
    _require(
        len(spectral) == 1
        and spectral[0]["seed_topology"]
        == "shared_spectrum_plus_separated_long_modes"
        and spectral[0]["spectral_sensitivity_control"]
        and spectral[0]["box_pair_id"] == "finite-3d-onset-s0",
        "separate long-mode spectral-sensitivity control drifted",
    )
    convergence = [
        case
        for case in cases
        if case["role"] == "finite_rigidity_3d_onset_matched_convergence_control"
    ]
    _require(len(convergence) == 4, "matched 3D convergence controls drifted")
    baseline = next(
        case for case in onset if case["case_id"] == "q019-fr-3d-onset-small-s0"
    )
    for case in convergence:
        for key in (
            "field_seed",
            "particle_seed",
            "rho_cr_over_rho0",
            "k0_rg0",
            "extents",
            "b_g",
            "k0",
            "initial_cr_momentum_flux_tensor",
            "initial_cr_kinetic_energy_density_rms_proxy",
            "conservative_cr_kinetic_energy_density_upper_bound",
        ):
            _require(case[key] == baseline[key], "matched 3D convergence control drifted")
    all_three_d = onset + convergence + spectral
    _require(
        all(
            not case["common_nonlinear_onset_resolution_reachable"]
            and case["common_nonlinear_onset_reachability_status"]
            == COMMON_ONSET_REACHABILITY_STATUS
            and math.isclose(
                float(case["preregistered_nonlinear_onset_Bperp_rms_over_B0"]),
                PREREGISTERED_NONLINEAR_ONSET_BPERP_RMS_OVER_B0,
            )
            and math.isclose(
                float(case["resolution_design_maximum_sampled_B_over_B0"]),
                PREREGISTERED_ONSET_RESOLUTION_MAXIMUM_SAMPLED_B_OVER_B0,
            )
            for case in all_three_d
        ),
        "3D onset reachability was not kept explicitly pilot-pending",
    )
    _require(
        next(
            case
            for case in convergence
            if case["case_id"] == "q019-fr-3d-onset-small-resolution-coarse-s0"
        )["coarse_resolution_stop_or_intermittency_limitation"],
        "coarse 3D resolution stop/intermittency limitation was not explicit",
    )
    runtime_regressions = [
        case
        for case in cases
        if case["role"] == "finite_rigidity_source_local_initializer_regression"
    ]
    _require(
        len(runtime_regressions) == 2
        and {int(case["ppc"]) for case in runtime_regressions}
        == {min(FINITE_PPC_LADDER), max(FINITE_PPC_LADDER)}
        and all(
            case["resource_classification"] == "source_local_runtime_regression"
            for case in runtime_regressions
        ),
        "finite initializer runtime regression inventory drifted",
    )
    isolation = [
        case
        for case in cases
        if case["role"] == "finite_rigidity_isolation_fixed_cr_distribution_control"
    ]
    _require(
        len(isolation) == len(FINITE_K0_RG0_GRID)
        and {case["k0_rg0"] for case in isolation} == set(FINITE_K0_RG0_GRID),
        "rigidity-isolation control inventory drifted",
    )
    for case in isolation[1:]:
        for key in (
            "rho_cr_over_rho0",
            "guide_parallel_stream_speed",
            "finite_shell_speed",
            "species_charge",
            "expected_j_over_c",
            "initial_cr_kinetic_energy_density_rms_proxy",
            "conservative_cr_kinetic_energy_density_upper_bound",
        ):
            _require(
                math.isclose(
                    float(case[key]), float(isolation[0][key]), rel_tol=1.0e-13
                ),
                "rigidity-isolation matching drifted",
            )
        for row, reference_row in zip(
            case["initial_cr_momentum_flux_tensor"],
            isolation[0]["initial_cr_momentum_flux_tensor"],
        ):
            _require(
                all(
                    math.isclose(float(value), float(reference), rel_tol=1.0e-13)
                    for value, reference in zip(row, reference_row)
                ),
                "rigidity-isolation momentum-flux matching drifted",
            )
    _require(
        all(
            not case["isolated_variable_claim_authorized"]
            and "not_universal_pure_rigidity_scaling"
            in str(case["control_interpretation"])
            for case in isolation
        ),
        "fixed-CR-distribution rigidity controls overclaim isolation",
    )


def _species_rows(case: Mapping[str, object]) -> list[dict[str, float]]:
    return [
        {
            "mass": 1.0,
            "charge": float(case["species_charge"]),
            "vx0": float(case["guide_parallel_stream_speed"]),
            "vy0": 0.0,
            "vz0": 0.0,
        }
    ]


def _render_deck(
    case: Mapping[str, object], *, matrix_identity_fingerprint_value: str
) -> str:
    """Render one exact non-authorizing candidate deck."""
    extents = [float(value) for value in case["extents"]]
    nx = [int(value) for value in case["nx"]]
    mb = [int(value) for value in case["meshblock_nx"]]
    high = case["branch"] == "high_rigidity_current_retention_candidate"
    quiet = case["finite_sampling_mode"] == FINITE_STRATIFIED_SAMPLING_MODE
    lines = [
        "# Q019 physics-first nonlinear Bell candidate. No execution authority.",
        "",
        "<comment>",
        f"problem = {case['case_id']} source-local candidate not authorized",
        "",
        "<job>",
        f"basename = {case['case_id']}",
        "",
        "<mesh>",
        f"nghost = {case['nghost']}",
    ]
    for axis in range(3):
        lines.extend(
            [
                f"nx{axis + 1} = {nx[axis]}",
                f"x{axis + 1}min = 0",
                f"x{axis + 1}max = {_float_token(extents[axis])}",
                f"ix{axis + 1}_bc = periodic",
                f"ox{axis + 1}_bc = periodic",
            ]
        )
    lines.extend(["", "<mesh_refinement>", "refinement = none", "", "<meshblock>"])
    lines.extend(f"nx{axis + 1} = {mb[axis]}" for axis in range(3))
    lines.extend(
        [
            "",
            "<time>",
            "evolution = dynamic",
            "integrator = rk2",
            f"cfl_number = {_float_token(float(case['cfl']))}",
            f"nlim = {int(case['cycle_limit'])}",
            f"tlim = {_float_token(float(case['terminal_time']))}",
            "ndiag = 1",
            "",
            "<mhd>",
            "eos = ideal",
            f"reconstruct = {case['reconstruct']}",
            f"rsolver = {case['rsolver']}",
            "gamma = 1.6666666666666667",
            "",
            "<particles>",
            "particle_type = cosmic_ray",
            f"ppc = {case['ppc']}",
            "pusher = boris_tsc",
            "nspecies = 1",
            f"cr_distribution = {case['distribution']}",
            "deposit_moments = true",
            "deposit_order = 2",
            f"deposit_qscale = {_float_token(float(case['deposit_qscale']))}",
            "couple_moments_to_mhd = true",
            "couple_j_to_efield_coeff = 1",
            "couple_j_to_efield_representation = cell_centered",
            "couple_j_deposition_mode = cc_convert",
            "couple_moments_momentum_to_mhd = true",
            "couple_moments_energy_to_mhd = true",
            "couple_fluid_feedback_order = mhd_src_terms",
            "couple_moments_momentum_coeff = 1",
            "couple_moments_energy_coeff = 1",
            f"cr_vx0 = {_float_token(float(case['guide_parallel_stream_speed']))}",
            "cr_vy0 = 0",
            "cr_vz0 = 0",
            "pic_physical_mode = paper_mhd_pic_vl2_tsc",
            "pic_background_mode = coupled",
            "pic_feedback_mode = coupled",
            "pic_interp_scheme = tsc",
            "pic_enable_2d3v = true",
            f"pic_cr_light_speed = {_float_token(ARTIFICIAL_LIGHT_SPEED)}",
            "pic_cr_initial_state = velocity",
            "pic_cr_hall_mode = off",
            "pic_wave_damping_mode = off",
            "pic_deltaf_mode = off",
            "pic_expanding_box_mode = off",
            f"pic_max_cell_cross = {case['pic_max_cell_cross']}",
            f"pic_theta_max = {_float_token(float(case['pic_theta_max']))}",
            "pic_load_balance_cost_per_particle = 1",
            "pic_sort_interval = 10",
            "pic_intermediate_arrays = auto",
            f"pic_random_seed = {case['particle_seed']}",
        ]
    )
    for index, species in enumerate(_species_rows(case)):
        lines.extend(
            [
                "",
                f"<species{index}>",
                f"mass = {_float_token(species['mass'])}",
                f"charge = {_float_token(species['charge'])}",
                f"vx0 = {_float_token(species['vx0'])}",
                f"vy0 = {_float_token(species['vy0'])}",
                f"vz0 = {_float_token(species['vz0'])}",
            ]
        )
    lines.extend(
        [
            "",
            "<problem>",
            f"pgen_name = {PGEN_NAME}",
            "user_work_in_loop = true",
            "user_hist = true",
            "",
            f"<{PGEN_BLOCK}>",
            f"case_id = {case['case_id']}",
            f"campaign_id = {case['campaign_id']}",
            f"branch = {case['branch']}",
            f"role = {case['role']}",
            f"dimension = {case['dimension']}",
            f"seed_topology = {case['seed_topology']}",
            f"finite_sampling_mode = {case['finite_sampling_mode']}",
            f"high_sampling_mode = {case['high_sampling_mode']}",
            "high_position_noise_reference_case_id = "
            f"{case['high_position_noise_reference_case_id']}",
            f"matrix_identity_fingerprint = {matrix_identity_fingerprint_value}",
            f"runtime_identity_checksum_status = {case['runtime_identity_checksum_status']}",
            "external_sha256_execution_receipt_required = true",
            "external_sha256_execution_receipt_bound = false",
            "deck_role = source_local_physics_first_candidate_not_authorized",
            "source_lineage = hardened_control_plane_q043_then_corrected_q023_then_q019_v4",
            f"matrix_scope = {MATRIX_SCOPE}",
            f"domain_time_status = {DOMAIN_TIME_STATUS}",
            "axis_alignment = B0_parallel_JCR_parallel_x1",
            f"rho = {_float_token(float(case['rho']))}",
            f"pressure = {_float_token(float(case['pressure']))}",
            f"b_g = {_float_token(float(case['b_g']))}",
            f"u_a = {_float_token(float(case['u_a']))}",
            f"wavelength = {_float_token(float(case['wavelength']))}",
            f"k0 = {_float_token(float(case['k0']))}",
            f"epsilon = {_float_token(float(case['epsilon']))}",
            f"omega = {_float_token(float(case['omega']))}",
            f"k0_rg0 = {_float_token(float(case['k0_rg0']))}",
            f"rho_cr_over_rho0 = {_float_token(float(case['rho_cr_over_rho0']))}",
            f"species_mass = {_float_token(float(case['species_mass']))}",
            f"n_cr = {_float_token(float(case['n_cr']))}",
            f"applicability_scope = {case['applicability_scope']}",
            "strong_shock_applicability_authorized = false",
            "no_hall_applicability = bounded_linear_correction_candidate_pending_review_not_accepted",
            "no_hall_applicability_accepted = false",
            "nonlinear_no_hall_applicability_status = "
            f"{case['nonlinear_no_hall_applicability_status']}",
            "evolving_local_hall_applicability_diagnostics_complete = false",
            "bounded_hall_omission_candidate = true",
            "bounded_hall_omission_review_complete = false",
            "mhd_resolved_scale_applicability_accepted = false",
            "R_much_less_than_one_applicability_accepted = false",
            f"no_hall_mapping_basis = {case['no_hall_mapping_basis']}",
            "species_q_over_mc_matches_background = true",
            "charge_density_ratio_equal_background_qom = "
            f"{_float_token(float(case['charge_density_ratio_equal_background_qom']))}",
            "hall_parameter_equal_background_qom = "
            f"{_float_token(float(case['hall_parameter_equal_background_qom']))}",
            "hall_parameter_from_current = "
            f"{_float_token(float(case['hall_parameter_from_current']))}",
            "background_q_over_mc_reference = "
            f"{_float_token(float(case['background_q_over_mc_reference']))}",
            "background_ion_gyrofrequency = "
            f"{_float_token(float(case['background_ion_gyrofrequency']))}",
            "background_ion_inertial_length = "
            f"{_float_token(float(case['background_ion_inertial_length']))}",
            "k0_background_ion_inertial_length = "
            f"{_float_token(float(case['k0_background_ion_inertial_length']))}",
            "hall_parameter_over_twice_k0_di = "
            f"{_float_token(float(case['hall_parameter_over_twice_k0_di']))}",
            "hall_order_unity_reference_margin = "
            f"{_float_token(float(case['hall_order_unity_reference_margin']))}",
            "hall_order_unity_reference = "
            f"{_float_token(float(case['hall_order_unity_reference']))}",
            "bai_hall_linear_factor = "
            f"{_float_token(float(case['bai_hall_linear_factor']))}",
            "bai_hall_growth_rate_fractional_shift = "
            f"{_float_token(float(case['bai_hall_growth_rate_fractional_shift']))}",
            "bai_hall_wavenumber_fractional_shift = "
            f"{_float_token(float(case['bai_hall_wavenumber_fractional_shift']))}",
            "bai_hall_growth_rate_reduction_factor = "
            f"{_float_token(float(case['bai_hall_growth_rate_reduction_factor']))}",
            "bai_hall_wavenumber_reduction_factor = "
            f"{_float_token(float(case['bai_hall_wavenumber_reduction_factor']))}",
            "minimum_active_dx_over_background_di = "
            f"{_float_token(float(case['minimum_active_dx_over_background_di']))}",
            "no_subion_cell_scale_envelope_satisfied = true",
            "background_q_over_mc_at_lambda_equal_one_reference = "
            f"{_float_token(float(case['background_q_over_mc_at_lambda_equal_one_reference']))}",
            f"cr_inertia_parameter = {_float_token(float(case['cr_inertia_parameter']))}",
            "cr_momentum_loading_parameter = "
            f"{_float_token(float(case['cr_momentum_loading_parameter']))}",
            "cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy = "
            f"{_float_token(float(case['cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy']))}",
            f"energy_loading_regime = {case['energy_loading_regime']}",
            f"energy_loading_gate_status = {case['energy_loading_gate_status']}",
            "energy_loading_grid_review_complete = false",
            "universal_saturation_inference_authorized = false",
            "feedback_force_parameter = "
            f"{_float_token(float(case['feedback_force_parameter']))}",
            "positive_finite_loading_accounting_satisfied = true",
            f"eigenmode_amplitude = {_float_token(float(case['eigenmode_amplitude']))}",
            f"broadband_amplitude = {_float_token(float(case['broadband_amplitude']))}",
            "separated_long_mode_amplitude = "
            f"{_float_token(float(case['separated_long_mode_amplitude']))}",
            f"finite_shell_speed = {_float_token(float(case['finite_shell_speed']))}",
            "characteristic_shell_p_iso_over_m = "
            f"{_float_token(float(case['characteristic_shell_p_iso_over_m']))}",
            f"nominal_rg0 = {_float_token(float(case['nominal_rg0']))}",
            "initial_nominal_rg0_over_max_active_dx = "
            f"{_float_token(float(case['initial_nominal_rg0_over_max_active_dx']))}",
            "required_initial_rg0_over_max_active_dx = "
            f"{_float_token(float(case['required_initial_rg0_over_max_active_dx']))}",
            "preregistered_nonlinear_onset_Bperp_rms_over_B0 = "
            f"{_float_token(float(case['preregistered_nonlinear_onset_Bperp_rms_over_B0']))}",
            "resolution_design_maximum_sampled_B_over_B0 = "
            f"{_float_token(float(case['resolution_design_maximum_sampled_B_over_B0']))}",
            "maximum_sampled_B_over_B0_before_resolution_stop = "
            f"{_float_token(float(case['maximum_sampled_B_over_B0_before_resolution_stop']))}",
            "common_nonlinear_onset_resolution_reachable = "
            f"{str(bool(case['common_nonlinear_onset_resolution_reachable'])).lower()}",
            "common_nonlinear_onset_reachability_status = "
            f"{case['common_nonlinear_onset_reachability_status']}",
            "coarse_resolution_stop_or_intermittency_limitation = "
            f"{str(bool(case['coarse_resolution_stop_or_intermittency_limitation'])).lower()}",
            "minimum_characteristic_shell_rl_over_dx = "
            f"{_float_token(MIN_CHARACTERISTIC_SHELL_RL_OVER_DX)}",
            "resolution_envelope_satisfied_by_design = true",
            "runtime_resolution_stop_controller_installed = "
            f"{str(bool(case['runtime_resolution_stop_controller_installed'])).lower()}",
            "runtime_resolution_stop_controller_status = "
            f"{case['runtime_resolution_stop_controller_status']}",
            "runtime_resolution_guard_pilot_qualified = false",
            "postprocessing_resolution_gate_required = "
            f"{str(bool(case['postprocessing_resolution_gate_required'])).lower()}",
            "runtime_box_edge_monitor_installed = true",
            "runtime_box_edge_monitor_enabled = false",
            "runtime_box_edge_stop_boundary_frozen = false",
            "runtime_box_edge_stop_armed = false",
            "runtime_box_edge_monitor_passive = true",
            "runtime_box_edge_monitor_selection_rule = "
            f"{case['runtime_box_edge_monitor_selection_rule']}",
            "runtime_box_edge_monitor_geometry = "
            f"{case['runtime_box_edge_monitor_geometry']}",
            "runtime_box_edge_monitor_diagnostic_authority = "
            f"{case['runtime_box_edge_monitor_diagnostic_authority']}",
            "runtime_box_edge_monitor_unique_mode_count = "
            f"{case['runtime_box_edge_monitor_unique_mode_count']}",
            "runtime_box_edge_monitor_cell_passes_per_sample = "
            f"{case['runtime_box_edge_monitor_cell_passes_per_sample']}",
            "runtime_box_edge_monitor_global_reductions_per_sample = "
            f"{case['runtime_box_edge_monitor_global_reductions_per_sample']}",
            "runtime_box_edge_monitor_dt = "
            f"{_float_token(float(case['runtime_box_edge_monitor_dt']))}",
            f"runtime_box_edge_stop_ppm = {case['runtime_box_edge_stop_ppm']}",
            f"runtime_box_edge_monitor_schema = {BOX_EDGE_MONITOR_SCHEMA}",
            "runtime_box_edge_monitor_next_nominal_time = "
            f"{_float_token(float(case['runtime_box_edge_monitor_dt']))}",
            "runtime_box_edge_monitor_last_nominal_time = -1",
            "runtime_box_edge_monitor_last_cycle = 0",
            "runtime_box_edge_monitor_completed_slots = 0",
            "runtime_box_edge_monitor_valid_samples = 0",
            "runtime_box_edge_monitor_skipped_slots = 0",
            "runtime_box_edge_monitor_last_status = 0",
            "runtime_box_edge_monitor_status_mask = 0",
            "runtime_box_edge_monitor_last_prior_time = 0",
            "runtime_box_edge_monitor_last_time = 0",
            "runtime_box_edge_monitor_last_power_fraction = -1",
            "runtime_box_edge_monitor_max_power_fraction = -1",
            "runtime_box_edge_monitor_last_fluctuation_mean = -1",
            "runtime_dominant_scale_stop_controller_installed = false",
            "raw_production_authorized = false",
            "nonlinear_saturation_claim_authorized = false",
            f"pic_max_cell_cross = {case['pic_max_cell_cross']}",
            f"pic_theta_max = {_float_token(float(case['pic_theta_max']))}",
            f"matched_control_family = {case['matched_control_family']}",
            f"control_interpretation = {case['control_interpretation']}",
            "isolated_variable_claim_authorized = false",
            f"box_pair_id = {case['box_pair_id']}",
            "saturation_candidate = "
            f"{str(bool(case['saturation_candidate'])).lower()}",
            "spectral_sensitivity_control = "
            f"{str(bool(case['spectral_sensitivity_control'])).lower()}",
            "minimum_active_extent_over_nominal_rg0 = "
            f"{_float_token(float(case['minimum_active_extent_over_nominal_rg0']))}",
            "minimum_active_extent_over_seed_wavelength = "
            f"{_float_token(float(case['minimum_active_extent_over_seed_wavelength']))}",
            "initial_particle_cell_crossing_dt_bound = "
            f"{_float_token(float(case['initial_particle_cell_crossing_dt_bound']))}",
            "initial_particle_gyro_dt_bound = "
            f"{_float_token(float(case['initial_particle_gyro_dt_bound']))}",
            f"initial_particle_timestep_limiter = {case['initial_particle_timestep_limiter']}",
            f"field_seed = {case['field_seed']}",
            (
                "finite_distribution = not_applicable_high_rigidity_cold_beam"
                if high
                else "finite_distribution = single_species_equal_weight_isotropic_"
                "shell_plus_guide_drift"
            ),
            (
                "finite_position_sampling = not_applicable_high_rigidity_branch"
                if high
                else "finite_position_sampling = "
                f"{FINITE_QUIET_POSITION_SAMPLING if quiet else FINITE_NOISE_POSITION_SAMPLING}"
            ),
            (
                "finite_gyrophase_position_coupling = not_applicable_high_rigidity_cold_beam"
                if high
                else "finite_gyrophase_position_coupling = "
                f"{FINITE_QUIET_GYROPHASE_POSITION_COUPLING if quiet else FINITE_NOISE_GYROPHASE_POSITION_COUPLING}"
            ),
            "high_position_sampling = "
            f"{case['high_sampling_mode'] if high else 'not_applicable'}",
            "finite_packet_grouping_contract = "
            f"{case['finite_packet_grouping_contract']}",
            "finite_packet_grouping_against_actual_initializer_proved = "
            f"{str(bool(case['finite_packet_grouping_against_actual_initializer_proved'])).lower()}",
            "finite_packet_independent_spatial_anchors_per_cell_expected = "
            f"{case['finite_packet_independent_spatial_anchors_per_cell_expected']}",
            "finite_density_parallel_current_quiet_by_construction = "
            f"{str(bool(case['finite_density_parallel_current_quiet_by_construction'])).lower()}",
            "finite_density_parallel_current_quiet_measured_and_accepted = false",
            "finite_initial_rho_jx_noise_pair_gate = "
            f"{case['finite_initial_rho_jx_noise_pair_gate']}",
            "finite_initial_rho_jx_noise_pair_case_id = "
            f"{case['finite_initial_rho_jx_noise_pair_case_id']}",
            "nested_spectrum_contract = box_convergence_uses_identical_shared_spectrum_"
            "and_separate_long_mode_control",
            "evolving_resolution_stop_contract = postprocessing_fail_closed_gate_required_"
            "runtime_controller_removed_pending_excluded_pilot_benchmark_and_freeze",
            "high_rigidity_validity_gate = measured_not_assumed",
            "high_rigidity_convergence_envelope = staged_2d_ppc_resolution_timestep_"
            "and_stochastic_position_controls_required_before_3d_nonlinear_use",
            "finite_rigidity_predecessor_gate = measured_complex_mode_response_and_"
            "convergence_required",
            "current_agreement_gate = deposited_grid_vs_reconstructed_particle_required",
            "morphology_and_relative_drift_diagnostics = required",
            "physical_pilot_gate = explicit_independent_q043_then_q023_then_finite_"
            "predecessor_source_compatibility_then_excluded_pilots",
            "q043_independent_raw_cycle_one_oracle_id = "
            f"{case['q043_independent_raw_cycle_one_oracle_id']}",
            "q043_independent_raw_cycle_one_oracle_bound = false",
            "q023_independent_linear_predecessor_id = "
            f"{case['q023_independent_linear_predecessor_id']}",
            "q023_independent_linear_predecessor_bound = false",
            "independent_prerequisites_complete = false",
            "nonlinear_execution_prerequisites_passed = false",
            "current_normalization = deposited_j_over_c_equals_2_b_g_k0",
            "deposit_qscale_semantics = root_cell_macro_mass_volume_aware",
            f"deposited_prtcl_rho_semantics = {case['deposited_prtcl_rho_semantics']}",
            "deposited_cr_mass_density_derivation = "
            f"{case['deposited_cr_mass_density_derivation']}",
            "species_mass_and_charge_bound_in_immutable_payload = true",
            f"configured_root_cell_volume = {_float_token(float(case['root_cell_volume']))}",
            f"configured_volume_mean_j_over_c = {_float_token(configured_j_over_c(case))}",
            f"configured_rho_cr_over_rho0 = {_float_token(configured_rho_cr_over_rho0(case))}",
            f"qualification_effect = {QUALIFICATION_EFFECT}",
            f"numeric_tolerance_status = {NUMERIC_TOLERANCE_STATUS}",
            f"resource_classification = {case['resource_classification']}",
            "launch_authorized = false",
            "policy_authorized = false",
            "claim_authorized = false",
            "physical_pilot_authorized = false",
        ]
    )
    for index, (variable, file_type, dt) in enumerate(REQUIRED_OUTPUTS, 1):
        lines.extend(["", f"<output{index}>", f"file_type = {file_type}"])
        if file_type != "rst":
            lines.append(f"variable = {variable}")
        lines.extend(
            [
                f"id = {variable}",
                f"dt = {_float_token(dt)}",
                "ghost_zones = false",
                "gid = -1",
                "data_format = %12.5e",
            ]
        )
        if file_type == "bin":
            lines.append("single_file_per_rank = false")
        elif file_type == "rst":
            lines.append("single_file_per_rank = false")
    lines.extend(
        [
            "",
            f"<output{len(REQUIRED_OUTPUTS) + 1}>",
            "file_type = hst",
            "dcycle = 1",
            "ghost_zones = false",
            "gid = -1",
            "user_hist_only = false",
            "data_format = %24.16e",
        ]
    )
    return "\n".join(lines) + "\n"


def render_deck(case: Mapping[str, object]) -> str:
    """Render one exact non-authorizing candidate deck."""
    fingerprint = str(case["matrix_identity_fingerprint"])
    _require(
        fingerprint == matrix_identity_fingerprint(case),
        "runtime matrix identity SHA-256 drifted before rendering",
    )
    return _render_deck(case, matrix_identity_fingerprint_value=fingerprint)


def parse_athinput_text(text: str) -> dict[str, dict[str, str]]:
    blocks: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    for line_number, raw in enumerate(text.splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            name = line[1:-1].strip()
            _require(name and name not in blocks, f"line {line_number}: invalid block")
            current = {}
            blocks[name] = current
            continue
        _require(current is not None and line.count("=") == 1, f"line {line_number}: malformed")
        name, value = (part.strip() for part in line.split("=", 1))
        _require(name and value and name not in current, f"line {line_number}: invalid parameter")
        current[name] = value
    return blocks


def validate_rendered_deck(case: Mapping[str, object], text: str) -> dict[str, object]:
    _require(text.count("ndiag = 1") == 1, "time diagnostic key inventory drifted")
    blocks = parse_athinput_text(text)
    actual_payload = deck_semantics_payload(blocks)
    expected_payload = matrix_identity_payload(case)
    _require(actual_payload == expected_payload, "immutable runtime deck semantics drifted")
    immutable_payload_sha256 = _sha256_bytes(actual_payload.encode("utf-8"))
    _require(
        immutable_payload_sha256 == str(case["matrix_identity_fingerprint"]),
        "immutable runtime deck SHA-256 drifted",
    )
    _require(blocks["problem"]["pgen_name"] == PGEN_NAME, "pgen identity drifted")
    _require(
        blocks["problem"]["user_work_in_loop"] == "true"
        and blocks["problem"]["user_hist"] == "true"
        and blocks["mesh_refinement"]["refinement"] == "none",
        "box-edge runtime callback or uniform-mesh contract drifted",
    )
    contract = blocks[PGEN_BLOCK]
    for key in (
        "case_id",
        "campaign_id",
        "branch",
        "role",
        "seed_topology",
        "finite_sampling_mode",
        "high_sampling_mode",
        "box_pair_id",
        "matched_control_family",
        "control_interpretation",
        "resource_classification",
        "matrix_identity_fingerprint",
        "runtime_identity_checksum_status",
    ):
        _require(contract[key] == str(case[key]), f"{key} drifted")
    for block_name, key in (
        ("mhd", "reconstruct"),
        ("mhd", "rsolver"),
    ):
        _require(
            blocks[block_name][key] == str(case[key]),
            f"{block_name}/{key} drifted",
        )
    _require(
        int(blocks["particles"]["pic_max_cell_cross"])
        == int(case["pic_max_cell_cross"]),
        "particles/pic_max_cell_cross drifted",
    )
    _require(
        blocks["job"]["basename"] == str(case["case_id"]),
        "runtime job/case identity drifted",
    )
    _require(
        int(contract["dimension"]) == int(case["dimension"])
        and int(contract["field_seed"]) == int(case["field_seed"])
        and int(blocks["particles"]["pic_random_seed"]) == int(case["particle_seed"])
        and math.isclose(
            float(blocks["particles"]["pic_theta_max"]),
            float(case["pic_theta_max"]),
            rel_tol=0.0,
            abs_tol=0.0,
        )
        and math.isclose(
            float(contract["pic_theta_max"]),
            float(case["pic_theta_max"]),
            rel_tol=0.0,
            abs_tol=0.0,
        )
        and math.isclose(
            float(contract["minimum_characteristic_shell_rl_over_dx"]),
            float(case["minimum_characteristic_shell_rl_over_dx"]),
            rel_tol=0.0,
            abs_tol=0.0,
        )
        and (contract["saturation_candidate"] == "true")
        == bool(case["saturation_candidate"])
        and (contract["spectral_sensitivity_control"] == "true")
        == bool(case["spectral_sensitivity_control"])
        and contract["runtime_box_edge_monitor_installed"] == "true"
        and contract["runtime_box_edge_monitor_enabled"] == "false"
        and contract["runtime_box_edge_stop_boundary_frozen"] == "false"
        and contract["runtime_box_edge_stop_armed"] == "false"
        and contract["runtime_box_edge_monitor_passive"] == "true"
        and contract["runtime_box_edge_monitor_selection_rule"]
        == BOX_EDGE_MONITOR_SELECTION_RULE
        and contract["runtime_box_edge_monitor_geometry"] == BOX_EDGE_MONITOR_GEOMETRY
        and contract["runtime_box_edge_monitor_diagnostic_authority"]
        == BOX_EDGE_MONITOR_DIAGNOSTIC_AUTHORITY
        and int(contract["runtime_box_edge_monitor_unique_mode_count"])
        == box_edge_unique_mode_count(int(case["dimension"]))
        and int(contract["runtime_box_edge_monitor_cell_passes_per_sample"])
        == box_edge_cell_passes_per_sample(int(case["dimension"]))
        and int(contract["runtime_box_edge_monitor_global_reductions_per_sample"])
        == box_edge_cell_passes_per_sample(int(case["dimension"]))
        and math.isclose(
            float(contract["runtime_box_edge_monitor_dt"]),
            float(case["runtime_box_edge_monitor_dt"]),
            rel_tol=0.0,
            abs_tol=0.0,
        )
        and int(contract["runtime_box_edge_stop_ppm"])
        == int(case["runtime_box_edge_stop_ppm"]),
        "runtime discrete matrix identity drifted",
    )
    _require(
        contract["launch_authorized"]
        == contract["policy_authorized"]
        == contract["claim_authorized"]
        == contract["physical_pilot_authorized"]
        == contract["strong_shock_applicability_authorized"]
        == contract["energy_loading_grid_review_complete"]
        == contract["universal_saturation_inference_authorized"]
        == contract["runtime_box_edge_stop_boundary_frozen"]
        == contract["runtime_box_edge_stop_armed"]
        == contract["runtime_dominant_scale_stop_controller_installed"]
        == contract["raw_production_authorized"]
        == contract["nonlinear_saturation_claim_authorized"]
        == contract["isolated_variable_claim_authorized"]
        == contract["finite_density_parallel_current_quiet_measured_and_accepted"]
        == contract["no_hall_applicability_accepted"]
        == contract["evolving_local_hall_applicability_diagnostics_complete"]
        == contract["bounded_hall_omission_review_complete"]
        == contract["mhd_resolved_scale_applicability_accepted"]
        == contract["R_much_less_than_one_applicability_accepted"]
        == contract["q043_independent_raw_cycle_one_oracle_bound"]
        == contract["q023_independent_linear_predecessor_bound"]
        == contract["independent_prerequisites_complete"]
        == contract["nonlinear_execution_prerequisites_passed"]
        == "false",
        "authority must remain false",
    )
    _require(
        contract["numeric_tolerance_status"] == NUMERIC_TOLERANCE_STATUS,
        "numeric threshold status drifted",
    )
    _require(
        contract["bounded_hall_omission_candidate"] == "true"
        and contract["no_subion_cell_scale_envelope_satisfied"] == "true"
        and contract["positive_finite_loading_accounting_satisfied"] == "true"
        and contract["deposited_prtcl_rho_semantics"]
        == DEPOSITED_PRTCL_RHO_SEMANTICS
        and contract["deposited_cr_mass_density_derivation"]
        == DEPOSITED_CR_MASS_DENSITY_DERIVATION
        and contract["species_mass_and_charge_bound_in_immutable_payload"] == "true"
        and contract["runtime_identity_checksum_status"]
        == RUNTIME_IDENTITY_CHECKSUM_STATUS
        and contract["external_sha256_execution_receipt_required"] == "true"
        and contract["external_sha256_execution_receipt_bound"] == "false"
        and contract["runtime_resolution_stop_controller_installed"] == "false"
        and contract["runtime_resolution_stop_controller_status"]
        == RUNTIME_RESOLUTION_STOP_CONTROLLER_STATUS
        and contract["runtime_resolution_guard_pilot_qualified"] == "false"
        and (contract["postprocessing_resolution_gate_required"] == "true")
        == bool(case["postprocessing_resolution_gate_required"])
        and contract["runtime_box_edge_monitor_installed"] == "true"
        and contract["runtime_box_edge_monitor_enabled"] == "false"
        and int(contract["runtime_box_edge_monitor_schema"]) == BOX_EDGE_MONITOR_SCHEMA
        and math.isclose(
            float(contract["runtime_box_edge_monitor_next_nominal_time"]),
            BOX_EDGE_MONITOR_DT,
            rel_tol=0.0,
            abs_tol=0.0,
        )
        and float(contract["runtime_box_edge_monitor_last_nominal_time"]) == -1.0
        and int(contract["runtime_box_edge_monitor_last_cycle"]) == 0
        and int(contract["runtime_box_edge_monitor_completed_slots"]) == 0
        and int(contract["runtime_box_edge_monitor_valid_samples"]) == 0
        and int(contract["runtime_box_edge_monitor_skipped_slots"]) == 0
        and int(contract["runtime_box_edge_monitor_last_status"]) == 0
        and int(contract["runtime_box_edge_monitor_status_mask"]) == 0
        and float(contract["runtime_box_edge_monitor_last_prior_time"]) == 0.0
        and float(contract["runtime_box_edge_monitor_last_time"]) == 0.0
        and float(contract["runtime_box_edge_monitor_last_power_fraction"]) == -1.0
        and float(contract["runtime_box_edge_monitor_max_power_fraction"]) == -1.0
        and float(contract["runtime_box_edge_monitor_last_fluctuation_mean"]) == -1.0,
        "mapping-candidate or runtime resolution-controller status drifted",
    )
    _require(
        int(blocks["time"]["nlim"]) == int(case["cycle_limit"])
        and int(blocks["time"]["nlim"]) != 0,
        "candidate deck cycle horizon drifted",
    )
    _require(float(blocks["time"]["tlim"]) > 0.0, "candidate deck has no positive horizon")
    output_blocks = [block for name, block in blocks.items() if name.startswith("output")]
    history_blocks = [block for block in output_blocks if block["file_type"] == "hst"]
    _require(
        len(output_blocks) == len(REQUIRED_OUTPUTS) + 1
        and len(history_blocks) == 1
        and history_blocks[0]["dcycle"] == "1"
        and history_blocks[0]["data_format"] == "%24.16e",
        "per-cycle actual-timestep history retention drifted",
    )
    nx = tuple(int(blocks["mesh"][f"nx{axis}"]) for axis in (1, 2, 3))
    meshblock_nx = tuple(
        int(blocks["meshblock"][f"nx{axis}"]) for axis in (1, 2, 3)
    )
    extents = tuple(
        float(blocks["mesh"][f"x{axis}max"]) - float(blocks["mesh"][f"x{axis}min"])
        for axis in (1, 2, 3)
    )
    _require(
        nx == tuple(int(value) for value in case["nx"])
        and meshblock_nx == tuple(int(value) for value in case["meshblock_nx"])
        and all(
            math.isclose(value, float(expected), rel_tol=1.0e-13)
            for value, expected in zip(extents, case["extents"])
        ),
        "runtime mesh matrix identity drifted",
    )
    root_cell_volume = _root_cell_volume(extents, nx)
    ppc = int(blocks["particles"]["ppc"])
    qscale = float(blocks["particles"]["deposit_qscale"])
    species = blocks["species0"]
    measured_j = (
        ppc
        * qscale
        * float(species["charge"])
        * float(species["vx0"])
        / root_cell_volume
    )
    measured_rho = ppc * qscale / root_cell_volume / float(case["rho"])
    _require(
        math.isclose(
            float(species["charge"]),
            float(contract["background_q_over_mc_reference"]),
            rel_tol=0.0,
            abs_tol=0.0,
        )
        and contract["species_q_over_mc_matches_background"] == "true",
        "rendered CR/background q/(mc) equality drifted",
    )
    _require(
        math.isclose(measured_j, float(case["expected_j_over_c"]), rel_tol=1.0e-13),
        "J_CR/c drifted",
    )
    _require(
        math.isclose(measured_rho, float(case["rho_cr_over_rho0"]), rel_tol=1.0e-13),
        "rho_CR/rho0 drifted",
    )
    initial_rg0_over_dx = float(contract["initial_nominal_rg0_over_max_active_dx"])
    required_initial_rg0_over_dx = float(contract["required_initial_rg0_over_max_active_dx"])
    _require(
        math.isclose(
            initial_rg0_over_dx,
            float(case["initial_nominal_rg0_over_max_active_dx"]),
            rel_tol=1.0e-13,
        )
        and math.isclose(
            required_initial_rg0_over_dx,
            float(case["required_initial_rg0_over_max_active_dx"]),
            rel_tol=1.0e-13,
        )
        and initial_rg0_over_dx >= required_initial_rg0_over_dx,
        "finite-resolution envelope drifted",
    )
    _require(
        math.isclose(
            float(contract["preregistered_nonlinear_onset_Bperp_rms_over_B0"]),
            float(case["preregistered_nonlinear_onset_Bperp_rms_over_B0"]),
        )
        and math.isclose(
            float(contract["resolution_design_maximum_sampled_B_over_B0"]),
            float(case["resolution_design_maximum_sampled_B_over_B0"]),
        )
        and math.isclose(
            float(contract["maximum_sampled_B_over_B0_before_resolution_stop"]),
            float(case["maximum_sampled_B_over_B0_before_resolution_stop"]),
        )
        and (
            contract["common_nonlinear_onset_resolution_reachable"] == "true"
        )
        == bool(case["common_nonlinear_onset_resolution_reachable"])
        and contract["common_nonlinear_onset_reachability_status"]
        == str(case["common_nonlinear_onset_reachability_status"])
        and (
            contract["coarse_resolution_stop_or_intermittency_limitation"] == "true"
        )
        == bool(case["coarse_resolution_stop_or_intermittency_limitation"]),
        "common nonlinear-onset resolution contract drifted",
    )
    _require(
        contract["nonlinear_no_hall_applicability_status"]
        == NONLINEAR_NO_HALL_APPLICABILITY_STATUS
        and contract["q043_independent_raw_cycle_one_oracle_id"]
        == Q043_INDEPENDENT_ORACLE_ID
        and contract["q023_independent_linear_predecessor_id"]
        == Q023_INDEPENDENT_PREDECESSOR_ID,
        "nonlinear applicability or independent prerequisite identity drifted",
    )
    return {
        "case_id": case["case_id"],
        "configured_volume_mean_j_over_c": measured_j,
        "configured_rho_cr_over_rho0": measured_rho,
        "resolution_envelope_satisfied_by_design": True,
        "no_hall_applicability_accepted": False,
        "bounded_hall_omission_candidate": True,
        "mhd_resolved_scale_applicability_accepted": False,
        "no_subion_cell_scale_envelope_satisfied": True,
        "species_q_over_mc_matches_background": True,
        "per_cycle_history_retains_actual_completed_step_dt": True,
        "registered_actual_timestep_history_binding_complete": False,
        "complete_runtime_deck_semantics_bound": False,
        "exact_immutable_runtime_deck_semantics_payload_matches": True,
        "immutable_runtime_deck_semantics_sha256_matches": True,
        "runtime_identity_checksum_is_cryptographic_integrity_binding": True,
        "immutable_runtime_deck_semantics_sha256": immutable_payload_sha256,
        "mutable_output_bookkeeping_excluded_from_immutable_payload": [
            "file_number",
            "last_time",
        ],
        "mutable_box_edge_monitor_state_excluded_from_immutable_payload": list(
            MUTABLE_BOX_EDGE_MONITOR_STATE_PARAMETERS
        ),
        "immutable_output_defaults_materialized_and_bound": True,
        "external_sha256_execution_receipt_required": True,
        "external_sha256_execution_receipt_bound": False,
        "external_sha256_receipt_must_bind_exact_immutable_payload": False,
        "compiled_case_registry_binds_immutable_payload_sha256": True,
        "deposited_prtcl_rho_semantics": "single_species_charge_density",
        "deposited_cr_mass_density_derivation": "prtcl_rho_times_species_mass_over_species_charge",
        "species_mass_and_charge_bound_in_immutable_payload": True,
        "independent_prerequisites_complete": False,
        "runtime_box_edge_monitor_installed": True,
        "runtime_box_edge_monitor_enabled": False,
        "runtime_box_edge_stop_boundary_frozen": False,
        "runtime_box_edge_stop_armed": False,
        "run_horizon_advances": True,
        "raw_production_authorized": False,
        "authority": {
            "launch": False,
            "policy": False,
            "qualification": False,
            "claim": False,
        },
    }


def build_deck_manifest() -> tuple[dict[str, object], dict[str, str]]:
    rendered: dict[str, str] = {}
    records = []
    for case in expected_cases():
        text = render_deck(case)
        filename = f"{case['case_id']}.athinput"
        rendered[filename] = text
        records.append(
            {
                **case,
                "path": (
                    "inputs/publication/q019_physics_first_nonlinear_bell_successor_v2/"
                    f"{filename}"
                ),
                "sha256": _sha256_bytes(text.encode("utf-8")),
                "validation": validate_rendered_deck(case, text),
            }
        )
    analysis_bindings = [
        {"path": path, "sha256": _sha256_file(REPO_ROOT / path)}
        for path in ANALYSIS_BINDING_PATHS
    ]
    finite_energy = [
        float(
            case[
                "cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy"
            ]
        )
        for case in records
        if case["branch"] != "high_rigidity_current_retention_candidate"
    ]
    high_energy = [
        float(
            case[
                "cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy"
            ]
        )
        for case in records
        if case["branch"] == "high_rigidity_current_retention_candidate"
    ]
    largest = max(records, key=lambda case: math.prod(case["nx"]) * int(case["ppc"]))
    largest_root_cells = math.prod(largest["nx"])
    largest_terminal_time = float(largest["terminal_time"])
    largest_monitor_tolerance = (
        64.0
        * 2.220446049250313e-16
        * max(1.0, abs(largest_terminal_time), BOX_EDGE_MONITOR_DT)
    )
    largest_monitor_slots = math.floor(
        (largest_terminal_time + largest_monitor_tolerance) / BOX_EDGE_MONITOR_DT
    )
    largest_monitor_groups = box_edge_reduction_group_count(int(largest["dimension"]))
    largest_monitor_passes = box_edge_cell_passes_per_sample(int(largest["dimension"]))
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q019_physics_first_nonlinear_bell_successor_v2_deck_manifest",
        "status": "source_local_physics_first_design_complete_execution_prohibited",
        "qualification_effect": QUALIFICATION_EFFECT,
        "authority": {
            "launch_authorized": False,
            "frontier_submission_authorized": False,
            "physical_pilot_authorized": False,
            "qualification_authorized": False,
            "scientific_claim_authorized": False,
            "publication_authorized": False,
            "raw_production_authorized": False,
            "nonlinear_saturation_claim_authorized": False,
        },
        "physics_grid": {
            "finite_k0_rg0": list(FINITE_K0_RG0_GRID),
            "finite_rho_cr_over_rho0": list(FINITE_RHO_CR_GRID),
            "finite_seed_pairs": [list(pair) for pair in SEED_PAIRS],
            "finite_ppc_ladder": list(FINITE_PPC_LADDER),
            "high_k0_rg0": HIGH_K0_RG0,
            "high_rho_cr_over_rho0": HIGH_FIDUCIAL_RHO_CR,
            "no_hall_applicability_accepted": False,
            "bounded_hall_omission_candidate": True,
            "bounded_hall_omission_review_complete": False,
            "mhd_resolved_scale_applicability_accepted": False,
            "R_much_less_than_one_applicability_accepted": False,
            "all_cr_species_q_over_mc": SIMILARITY_SCALED_ION_Q_OVER_MC,
            "background_ion_q_over_mc_mapping": BACKGROUND_Q_OVER_MC_MAPPING,
            "cr_and_background_q_over_mc_equal_in_every_case": True,
            "fiducial_background_ion_gyrofrequency": BACKGROUND_ION_GYROFREQUENCY,
            "fiducial_background_ion_inertial_length": BACKGROUND_ION_INERTIAL_LENGTH,
            "fiducial_k0_background_ion_inertial_length": (
                K0_BACKGROUND_ION_INERTIAL_LENGTH
            ),
            "Hall_order_unity_reference_not_acceptance_threshold": (
                HALL_ORDER_UNITY_REFERENCE
            ),
            "maximum_mapped_Lambda_Hall": max(
                float(case["hall_parameter_from_current"]) for case in records
            ),
            "minimum_order_unity_reference_over_mapped_Lambda_Hall": min(
                float(case["hall_order_unity_reference_margin"]) for case in records
            ),
            "maximum_charge_density_ratio_R": max(
                float(case["charge_density_ratio_equal_background_qom"])
                for case in records
            ),
            "minimum_materialized_dx_over_di": min(
                float(case["minimum_active_dx_over_background_di"])
                for case in records
            ),
            "maximum_absolute_Bai_growth_rate_fractional_shift": max(
                abs(float(case["bai_hall_growth_rate_fractional_shift"]))
                for case in records
            ),
            "maximum_absolute_Bai_wavenumber_fractional_shift": max(
                abs(float(case["bai_hall_wavenumber_fractional_shift"]))
                for case in records
            ),
            "no_hall_accounting": (
                "use one equal similarity-scaled CR/background q/(mc)=10000 mapping; "
                "derive d_i, resolved dx/d_i, exact Bai R=n_CR/(n_i+n_CR), "
                "Lambda=R*u_d/U_A=2*k0*d_i*(1-R), and signed Bai growth/wavenumber "
                "reductions. The bounded omission and MHD-scale mapping remain "
                "review-pending because the retained literature does not provide "
                "numeric much-less/much-larger acceptance cutoffs"
            ),
            "positive_finite_loading_accounting_satisfied_for_all_rows": True,
            "applicability_scope": APPLICABILITY_SCOPE,
            "strong_shock_applicability_authorized": False,
        },
        "high_rigidity_staged_convergence_envelope": {
            "ppc_ladder": list(HIGH_PPC_LADDER),
            "fiducial_ppc": HIGH_FIDUCIAL_PPC,
            "required_axes": [
                "PPC",
                "resolution",
                "timestep",
                "cell_centered_vs_stochastic_position_sampling",
                "multiple_particle_seeds_for_stochastic_sampling",
            ],
            "control_case_ids": [
                str(case["case_id"])
                for case in expected_cases()
                if case["role"] == "high_rigidity_staged_convergence_control"
            ],
            "scope": "2d_current_retention_only_no_high_rigidity_saturation_claim",
            "exact_volume_aware_j_over_c_and_rho_cr_closure_required": True,
            "deposited_prtcl_rho_semantics": "single_species_charge_density",
            "deposited_cr_mass_density_derivation": (
                "prtcl_rho_times_species_mass_over_species_charge"
            ),
            "species_mass_and_charge_bound_in_immutable_payload": True,
            "runtime_controls_materialized": False,
            "nonlinear_saturation_claim_authorized": False,
        },
        "rigidity_isolation_controls": {
            "case_ids": [
                str(case["case_id"])
                for case in records
                if case["role"]
                == "finite_rigidity_isolation_fixed_cr_distribution_control"
            ],
            "fixed_quantities": [
                "rho_CR/rho0",
                "full_initial_CR_velocity_distribution",
                "species_q_over_mc",
                "deposited_J_CR/c",
                "absolute_CR_momentum_flux_tensor",
                "initial_CR_kinetic_energy_density_RMS_proxy",
                "conservative_CR_kinetic_energy_density_upper_bound",
            ],
            "varied_quantities": ["B0", "k0", "k0_r_g0"],
            "also_changes": ["U_A", "magnetic_energy_density", "k0_d_i"],
            "interpretation": (
                "fixed-CR-distribution and fixed-absolute-loading background-field "
                "rigidity control, not a universal pure-rigidity scan"
            ),
            "isolated_variable_claim_authorized": False,
            "universal_rigidity_scaling_authorized": False,
        },
        "coupled_response_grid": {
            "case_ids": [
                str(case["case_id"])
                for case in records
                if case["role"]
                == "finite_rigidity_density_drift_coupled_response_ensemble"
            ],
            "varied_inputs": ["k0_r_g0", "rho_CR/rho0"],
            "coupled_changes": [
                "guide_parallel_drift",
                "isotropic_shell_speed_and_pressure",
                "initial_anisotropic_momentum_flux",
                "kinetic_loading",
            ],
            "isolated_loading_or_rigidity_inference_authorized": False,
        },
        "energy_loading_boundary": {
            "regime": ENERGY_LOADING_REGIME,
            "reported_dimension": (
                "cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy"
            ),
            "finite_branch_range": [min(finite_energy), max(finite_energy)],
            "high_rigidity_branch_range": [min(high_energy), max(high_energy)],
            "gate_status": ENERGY_LOADING_GATE_STATUS,
            "grid_and_conservation_response_review_complete": False,
            "universal_saturation_inference_authorized": False,
        },
        "nested_spectrum_box_contract": {
            "shared_realization": "identical fixed physical shared modes and phases",
            "box_convergence_large_box_addition": "none",
            "long_mode_addition": "separate spectral-sensitivity control only",
            "shared_and_long_mode_power_must_be_reported_separately": True,
            "box_sufficiency_claim_authorized": False,
        },
        "finite_resolution_envelope": {
            "target_nonlinear_B_over_B0": None,
            "preregistered_nonlinear_onset_Bperp_rms_over_B0": (
                PREREGISTERED_NONLINEAR_ONSET_BPERP_RMS_OVER_B0
            ),
            "preregistered_resolution_design_maximum_sampled_B_over_B0": (
                PREREGISTERED_ONSET_RESOLUTION_MAXIMUM_SAMPLED_B_OVER_B0
            ),
            "common_nonlinear_onset_reachability_status": (
                COMMON_ONSET_REACHABILITY_STATUS
            ),
            "all_3d_rows_reach_common_nonlinear_onset_resolution_envelope": False,
            "reachability_inferred_from_maximum_sampled_B_envelope": False,
            "coarse_resolution_stop_or_intermittency_limitation_explicit": True,
            "rejected_predecessor_target_B_over_B0": 8.0,
            "rejected_predecessor_audit_absolute_energy_bound_B_over_B0": 4.63,
            "rejected_predecessor_target_reason": (
                "energetically_unreachable_and_resolution_guard_incompatible"
            ),
            "minimum_characteristic_shell_rl_over_dx": (
                MIN_CHARACTERISTIC_SHELL_RL_OVER_DX
            ),
            "required_initial_nominal_rg0_over_dx_for_3d_onset": (
                REQUIRED_INITIAL_RL_OVER_DX
            ),
            "predictive_amplification_times_resolution_acceptance_rule": None,
            "predictive_rule_reason": (
                "max sampled B/B0 is a stop-envelope diagnostic, not proof that a "
                "coarse or intermittent row can reach a common nonlinear-onset amplitude"
            ),
            "stop_rule": (
                "postprocessing fails closed on characteristic-shell rL/dx and retained "
                "actual-timestep gyro-angle evidence; the unbounded per-cycle particle "
                "scan and global reductions were removed, and any future runtime "
                "controller requires excluded-pilot benchmarking and a frozen cadence"
            ),
            "runtime_stop_final_status": None,
            "runtime_stop_process_exit": None,
            "stopped_run_saturation_evidence_eligible": False,
            "acceptance_authorized": False,
            "runtime_resolution_stop_controller_installed": False,
            "runtime_resolution_stop_controller_status": (
                RUNTIME_RESOLUTION_STOP_CONTROLLER_STATUS
            ),
            "runtime_resolution_guard_pilot_qualified": False,
            "postprocessing_resolution_gate_required": True,
            "runtime_box_edge_monitor_installed": True,
            "runtime_box_edge_monitor_enabled": False,
            "runtime_box_edge_monitor_enabled_case_ids": [],
            "runtime_box_edge_monitor_definition": (
                "fraction of mean-subtracted transverse magnetic power in every "
                "nonzero discrete mode inside the frozen physical-wave-number ball "
                "|k| <= sqrt(dimension) * 2pi / shortest_active_extent; conjugate "
                "pairs are accumulated once, giving 7 unique modes in 2D and 23 in "
                "3D for the frozen 2:1 x1-to-transverse production geometry"
            ),
            "runtime_box_edge_monitor_selection_rule": BOX_EDGE_MONITOR_SELECTION_RULE,
            "runtime_box_edge_monitor_geometry": BOX_EDGE_MONITOR_GEOMETRY,
            "runtime_box_edge_monitor_dt": BOX_EDGE_MONITOR_DT,
            "runtime_box_edge_monitor_schema": BOX_EDGE_MONITOR_SCHEMA,
            "runtime_box_edge_monitor_passive": True,
            "runtime_box_edge_monitor_diagnostic_authority": (
                BOX_EDGE_MONITOR_DIAGNOSTIC_AUTHORITY
            ),
            "runtime_box_edge_monitor_unavailable_statuses": [
                "cadence_skipped",
                "no_global_cells",
                "zero_fluctuation",
                "numerical_unavailable",
            ],
            "runtime_box_edge_monitor_unavailable_metric_is_fail_closed": True,
            "runtime_box_edge_monitor_can_request_stop": False,
            "runtime_box_edge_monitor_restart_binds_prior_completed_time": True,
            "runtime_box_edge_stop_boundary_frozen": False,
            "runtime_box_edge_stop_armed": False,
            "runtime_dominant_scale_stop_controller_installed": False,
            "raw_production_blocked_until_reviewed_box_edge_stop_controller": True,
            "exact_next_runtime_implementation": (
                "benchmark and physically review the passive monitor in excluded "
                "pilots, then separately implement and review any dominant-scale stop "
                "controller; this quarantined monitor remains unable to stop a run"
            ),
        },
        "high_rigidity_validity_gate": {
            "required_measurements": [
                "lab_frame_current_retention",
                "CR_parallel_momentum_fractional_change",
                "CR_kinetic_energy_fractional_change",
                "gas_parallel_acceleration",
                "relative_CR_gas_drift",
                "sampled_particle_rg_over_dominant_wavelength",
            ],
            "numeric_thresholds": None,
            "threshold_source": "excluded_pilots_and_external_review",
            "fixed_current_like_label_authorized": False,
        },
        "high_rigidity_saturation_mechanism_diagnostic_contract": {
            "required_discriminants": [
                "Bperp_amplification",
                "dominant_wavenumber_and_long_mode_transfer",
                "density_Bperp2_correlation_cavities_and_filaments",
                "current_retention_CR_momentum_and_CR_energy_change",
                "gas_acceleration_and_relative_drift",
                "particle_rg_over_dominant_wavelength",
                "gas_magnetic_and_particle_feedback_energy_transfer",
                "initial_rho_and_current_spatial_noise",
            ],
            "numeric_thresholds": None,
            "mechanism_classification_authorized": False,
            "nonlinear_saturation_claim_authorized": False,
        },
        "finite_3d_nonlinear_onset_prerequisite_contract": {
            "box_pair_ids": sorted(
                {
                    str(case["box_pair_id"])
                    for case in records
                    if str(case["role"]).startswith(
                        "finite_rigidity_3d_nonlinear_onset_box_"
                    )
                }
            ),
            "matched_convergence_axes": [
                "three_deterministic_seeds",
                "PPC",
                "resolution",
                "configured_particle_step_control",
                "actual_completed_step_dt_history",
                "identical_shared_spectrum_box_size",
                "separate_long_mode_spectral_sensitivity",
            ],
            "required_gates": [
                "nonlinear_amplitude",
                "plateau_log_slope",
                "sustained_post_plateau_window",
                "gas_plus_CR_energy_conservation",
                "gas_plus_CR_momentum_conservation",
                "current_retention",
                "relative_drift_change",
                "backreaction_energy_transfer",
                "evolving_rL_over_dx",
                "dominant_scale_vs_box",
                "matched_box_convergence",
            ],
            "numeric_thresholds": None,
            "threshold_source": "excluded_pilots_and_preregistered_review",
            "current_rows_are_saturation_candidates": False,
            "target_B_over_B0": None,
            "preregistered_nonlinear_onset_Bperp_rms_over_B0": (
                PREREGISTERED_NONLINEAR_ONSET_BPERP_RMS_OVER_B0
            ),
            "resolution_design_maximum_sampled_B_over_B0": (
                PREREGISTERED_ONSET_RESOLUTION_MAXIMUM_SAMPLED_B_OVER_B0
            ),
            "common_nonlinear_onset_reachability_status": (
                COMMON_ONSET_REACHABILITY_STATUS
            ),
            "all_rows_reach_common_nonlinear_onset_resolution_envelope": False,
            "coarse_resolution_stop_or_intermittency_limitation_explicit": True,
            "nonlinear_saturation_claim_authorized": False,
        },
        "finite_rigidity_predecessor": {
            "reference": (
                "highest-resolution/highest-PPC nested-Haar octahedral isotropic-shell "
                "complex Fourier-mode response for this exact MHD-PIC system"
            ),
            "measurements": [
                "signed_complex_Bperp_mode",
                "signed_complex_deposited_Jperp_mode",
                "growth_and_frequency_fit",
                "polarization",
                "current_retention",
                "force_noise",
                "resolution_PPC_timestep_and_noise_convergence",
            ],
            "fit_windows_and_numeric_tolerances": None,
            "threshold_source": "excluded_predecessor_pilots_where_literature_is_silent",
            "runtime_predecessor_complete": False,
        },
        "diagnostic_contract": {
            "deposited_grid_vs_reconstructed_particle_current_agreement_required": True,
            "particle_momentum_flux_tensor_retained": True,
            "gas_plus_CR_conservation_reducer_data_retained": True,
            "morphology_metrics_required": True,
            "relative_drift_required": True,
            "multiple_deterministic_seeds_required": True,
            "configured_particle_step_limiter_is_actual_timestep_evidence": False,
            "per_cycle_history_retains_time_and_actual_completed_step_dt": True,
            "registered_actual_timestep_history_binding_complete": False,
            "actual_timestep_convergence_gate_passed": False,
            "finite_quiet_packet_grouping_contract": (
                FINITE_QUIET_PACKET_GROUPING_CONTRACT
            ),
            "finite_independent_position_packet_grouping_contract": (
                FINITE_NOISE_PACKET_GROUPING_CONTRACT
            ),
            "finite_packet_grouping_against_actual_initializer_proved_scope": (
                "quiet_collocated_path_only"
            ),
            "finite_quiet_isotropic_shell_sampling_mode": FINITE_STRATIFIED_SAMPLING_MODE,
            "finite_density_parallel_current_quiet_measured_and_accepted": False,
            "initial_rho_jx_noise_pair_gate": FINITE_INITIAL_NOISE_GATE,
            "exact_matched_independent_position_pairs": [
                [
                    "q019-fr-grid-k8-rho3em06-s0",
                    "q019-fr-fiducial-noise-seeded-s0",
                ],
                [
                    "q019-fr-predecessor-k8-rho3em06-s0",
                    "q019-fr-predecessor-fiducial-noise-seeded-s0",
                ],
            ],
            "numeric_acceptance_thresholds": None,
        },
        "future_provenance_boundary": {
            "self_attested_trusted_receipts_accepted": False,
            "caller_supplied_scheduler_facts_accepted": False,
            "raw_science_analysis_enabled": False,
            "runtime_identity_checksum_algorithm": "SHA-256",
            "runtime_identity_checksum_is_cryptographic_integrity_binding": True,
            "mutable_output_bookkeeping_excluded_from_immutable_payload": [
                "file_number",
                "last_time",
            ],
            "immutable_output_defaults_materialized_and_bound": True,
            "external_sha256_execution_receipt_required": True,
            "external_sha256_execution_receipt_bound": False,
            "external_sha256_receipt_must_bind_exact_immutable_payload": False,
            "external_sha256_receipt_binding_target": (
                "trusted source archive, executable, launch deck, and runtime artifacts"
            ),
            "required_path": (
                "future hardened installed control plane and Q043 registered-admission "
                "adapter with immutable execution and artifact bindings"
            ),
            "runtime_source_closure": runtime_source_closure(),
        },
        "independent_prerequisite_boundary": {
            "q043_independent_raw_cycle_one_oracle_id": Q043_INDEPENDENT_ORACLE_ID,
            "q043_independent_raw_cycle_one_oracle_bound": False,
            "q023_independent_linear_predecessor_id": Q023_INDEPENDENT_PREDECESSOR_ID,
            "q023_independent_linear_predecessor_bound": False,
            "independent_prerequisites_complete": False,
            "nonlinear_execution_prerequisites_passed": False,
            "matrix_authorizes_execution": False,
        },
        "nonlinear_no_hall_applicability_boundary": {
            "status": NONLINEAR_NO_HALL_APPLICABILITY_STATUS,
            "evolving_local_hall_applicability_diagnostics_complete": False,
            "no_hall_applicability_accepted": False,
            "required_evolving_local_diagnostics": [
                "rho_CR_over_rho_gas",
                "Bai_R",
                "CR_minus_gas_relative_drift",
                "background_ion_inertial_length",
                "dx_over_background_ion_inertial_length",
                "Lambda",
            ],
        },
        "staged_resource_boundary": {
            "largest_case_id": largest["case_id"],
            "largest_case_root_cells": largest_root_cells,
            "largest_case_macro_particles": largest_root_cells * int(largest["ppc"]),
            "box_edge_monitor_cost_model": {
                "physical_low_k_unique_modes": {"2d": 7, "3d": 23},
                "mode_groups_of_five": {"2d": 2, "3d": 5},
                "full_grid_passes_per_valid_sample": {"2d": 3, "3d": 6},
                "global_sum_reductions_per_valid_sample": {"2d": 3, "3d": 6},
                "trigonometric_evaluations_per_cell_per_mode_group": 6,
                "mode_component_accumulations_per_cell_per_valid_sample": {
                    "2d": 4 * 7,
                    "3d": 4 * 23,
                },
                "complex_multiplications_per_cell_per_valid_sample": {
                    "2d": box_edge_complex_multiplications_per_cell(2),
                    "3d": box_edge_complex_multiplications_per_cell(3),
                },
                "largest_case_nominal_slots_through_configured_terminal_time": (
                    largest_monitor_slots
                ),
                "matrix_enabled_case_ids": [],
                "matrix_enabled_case_count": 0,
                "matrix_actual_full_grid_cell_visits": 0,
                "matrix_actual_trigonometric_evaluations": 0,
                "matrix_actual_mode_component_accumulations": 0,
                "matrix_actual_complex_multiplications": 0,
                "largest_case_upper_bound_full_grid_cell_visits": (
                    largest_root_cells * largest_monitor_passes * largest_monitor_slots
                ),
                "largest_case_upper_bound_trigonometric_evaluations": (
                    largest_root_cells
                    * largest_monitor_groups
                    * 6
                    * largest_monitor_slots
                ),
                "largest_case_upper_bound_mode_component_accumulations": (
                    largest_root_cells
                    * 4
                    * box_edge_unique_mode_count(int(largest["dimension"]))
                    * largest_monitor_slots
                ),
                "largest_case_upper_bound_complex_multiplications": (
                    largest_root_cells
                    * box_edge_complex_multiplications_per_cell(
                        int(largest["dimension"])
                    )
                    * largest_monitor_slots
                ),
                "counts_exclude_kokkos_and_mpi_implementation_overheads": True,
                "excluded_pilot_benchmark_required": True,
                "performance_acceptance_thresholds": None,
                "production_promotion_blocked_pending_benchmark": True,
            },
            "three_d_resource_classifications": {
                str(case["case_id"]): case["resource_classification"]
                for case in records
                if int(case["dimension"]) == 3
            },
            "source_local_initializer_regression_case_ids": [
                str(case["case_id"])
                for case in records
                if case["resource_classification"] == "source_local_runtime_regression"
            ],
            "high_rigidity_convergence_controls_are_2d_staged": True,
            "resource_model_frozen": False,
            "execution_authorized": False,
        },
        "literature_context": {
            "bell_2004_doi": BELL_2004_DOI,
            "riquelme_spitkovsky_2009_doi": RIQUELME_SPITKOVSKY_2009_DOI,
            "gargate_2010_doi": GARGATE_2010_DOI,
            "zacharegkas_2024_doi": ZACHAREGKAS_2024_DOI,
            "sun_bai_2023_arxiv": SUN_BAI_2023_ARXIV,
            "bai_2015_arxiv": BAI_2015_ARXIV,
            "retained_primary_reference_map": RETAINED_REFERENCE_MAP,
            "hall_order_unity_reference_is_accepted_threshold": False,
            "exact_finite_isotropic_shell_numeric_tolerances_supplied_by_literature": False,
        },
        "physical_pilot_gate_order": list(PHYSICAL_PILOT_GATE_ORDER),
        "analysis_bindings": analysis_bindings,
        "decks": records,
    }
    return manifest, rendered


def materialize_checked_in_decks(*, replace: bool = False) -> dict[str, object]:
    manifest, rendered = build_deck_manifest()
    if CHECKED_IN_DECK_ROOT.exists():
        _require(replace, "checked-in deck root exists; pass --replace")
        shutil.rmtree(CHECKED_IN_DECK_ROOT)
    CHECKED_IN_DECK_ROOT.mkdir(parents=True)
    for filename, text in rendered.items():
        (CHECKED_IN_DECK_ROOT / filename).write_text(text, encoding="utf-8")
    CHECKED_IN_MANIFEST.write_bytes(_canonical_json_bytes(manifest))
    return manifest


def validate_checked_in_decks() -> dict[str, object]:
    manifest, rendered = build_deck_manifest()
    _require(CHECKED_IN_MANIFEST.is_file(), "checked-in manifest is missing")
    actual = json.loads(CHECKED_IN_MANIFEST.read_text(encoding="utf-8"))
    _require(actual == manifest, "checked-in manifest drifted")
    expected_names = set(rendered) | {"deck_manifest.json"}
    actual_names = {path.name for path in CHECKED_IN_DECK_ROOT.iterdir() if path.is_file()}
    _require(actual_names == expected_names, "checked-in deck inventory drifted")
    for filename, text in rendered.items():
        _require(
            (CHECKED_IN_DECK_ROOT / filename).read_text(encoding="utf-8") == text,
            f"checked-in deck drifted: {filename}",
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
        "choose exactly one checked-in deck operation",
    )
    result = (
        materialize_checked_in_decks(replace=args.replace)
        if args.materialize_checked_in_decks
        else validate_checked_in_decks()
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
