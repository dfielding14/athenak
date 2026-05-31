"""Bounded Q-011 Section 5.4 paper-shock preparation contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
from pathlib import Path

logger = logging.getLogger("athena" + __name__[7:])

_REPO_ROOT = Path(__file__).resolve().parents[3]
_INPUT_DECK = (
    _REPO_ROOT
    / "inputs"
    / "publication"
    / "pic_parallel_shock_section54_paper.athinput"
)
_PGEN_SOURCE = _REPO_ROOT / "src" / "pgen" / "tests" / "pic_parallel_shock.cpp"
_PUSHER_SOURCE = _REPO_ROOT / "src" / "particles" / "particles_pushers.cpp"
_RESULTS = {}

_EXPECTED_VALUES = {
    ("mesh", "nx1"): "4000",
    ("mesh", "x1min"): "0.0",
    ("mesh", "x1max"): "48000.0",
    ("mesh", "ix1_bc"): "reflect",
    ("mesh", "ox1_bc"): "inflow",
    ("mesh", "nx2"): "260",
    ("mesh", "x2min"): "0.0",
    ("mesh", "x2max"): "3120.0",
    ("mesh", "ix2_bc"): "periodic",
    ("mesh", "ox2_bc"): "periodic",
    ("mesh", "nx3"): "1",
    ("mesh", "x3min"): "0.0",
    ("mesh", "x3max"): "1.0",
    ("meshblock", "nx1"): "20",
    ("meshblock", "nx2"): "20",
    ("meshblock", "nx3"): "1",
    ("mesh_refinement", "refinement"): "adaptive",
    ("mesh_refinement", "num_levels"): "3",
    ("time", "tlim"): "1200.0",
    ("mhd", "eos"): "ideal",
    ("mhd", "gamma"): "1.66666666667",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "ppc"): "0.0",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "pic_enable_2d3v"): "true",
    ("particles", "nspecies"): "1",
    ("particles", "deposit_moments"): "true",
    ("particles", "deposit_qscale"): "9.0e-4",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "couple_j_to_efield_representation"): "cell_centered",
    ("particles", "couple_j_deposition_mode"): "cc_convert",
    ("particles", "couple_fluid_feedback_order"): "mhd_src_terms",
    ("particles", "couple_moments_momentum_to_mhd"): "true",
    ("particles", "couple_moments_energy_to_mhd"): "true",
    ("particles", "pic_physical_mode"): "paper_mhd_pic",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_cr_light_speed"): "10000.0",
    ("particles", "pic_cr_initial_state"): "momentum",
    ("particles", "pic_cr_hall_mode"): "off",
    ("particles", "pic_wave_damping_mode"): "off",
    ("particles", "pic_theta_max"): "0.3",
    ("particles", "pic_load_balance_cost_per_particle"): "0.0",
    ("problem", "pgen_name"): "pic_parallel_shock",
    ("problem", "ps_rho0"): "1.0",
    ("problem", "ps_p0"): "0.10",
    ("problem", "ps_u0"): "30.0",
    ("problem", "ps_b0"): "1.0",
    ("problem", "ps_eta"): "1.0e-3",
    ("problem", "ps_vinj_over_u0"): "3.16227766017",
    ("problem", "ps_shock_speed_model"): "ideal_surface",
    ("problem", "ps_inject_t_start"): "0.0",
    ("problem", "ps_remove_birth_time_before"): "45.0",
    ("problem", "ps_enable_injection"): "true",
    ("problem", "ps_enable_gas_subtraction"): "true",
    ("problem", "ps_enable_curvature_amr"): "true",
    ("problem", "ps_refine_curv"): "1.0",
    ("problem", "ps_derefine_curv"): "0.1",
    ("problem", "ps_inject_species"): "0",
    ("problem", "ps_inject_seed"): "23050101",
    ("problem", "ps_enable_frame_tracking"): "false",
    ("output1", "file_type"): "bin",
    ("output1", "variable"): "mhd_w_d",
    ("output1", "id"): "rho",
    ("output1", "dt"): "100.0",
    ("output2", "file_type"): "bin",
    ("output2", "variable"): "mhd_bmag",
    ("output2", "id"): "bmag",
    ("output2", "dt"): "100.0",
    ("output3", "file_type"): "bin",
    ("output3", "variable"): "mhd_j2",
    ("output3", "id"): "j2",
    ("output3", "dt"): "100.0",
    ("output4", "file_type"): "pvtk",
    ("output4", "variable"): "prtcl_all",
    ("output4", "id"): "prtcl_all",
    ("output4", "dt"): "100.0",
    ("output5", "file_type"): "rst",
    ("output5", "dt"): "100.0",
    ("species0", "mass"): "1.0",
    ("species0", "charge"): "1.0",
}

_OPEN_ITEMS = [
    "executed_shock_surface_injection_distribution_audit",
    "gas_pressure_thermodynamic_normalization_audit",
    "frontier_load_balance_cost_tuning",
    "snapshot_time_selection_tolerance",
    "downstream_spectrum_fit_energy_interval",
    "amr_vs_fine_uniform_residual_tolerances",
    "clean_candidate_frontier_executable_orion_root_and_campaign_execution",
    "independent_raw_artifact_recompute_and_external_review",
]


class ContractError(ValueError):
    """Raised when the frozen Q-011 preparation contract is not satisfied."""


def _require_close(label: str, measured: float, expected: float) -> None:
    if not math.isfinite(measured) or not math.isfinite(expected):
        raise ContractError(f"{label}: values must be finite")
    if not math.isclose(measured, expected, rel_tol=1.0e-11, abs_tol=1.0e-12):
        raise ContractError(f"{label}: expected {expected!r}, measured {measured!r}")


def _parse_positive_float(label: str, value: str) -> float:
    measured = _parse_finite_float(label, value)
    if measured <= 0.0:
        raise ContractError(f"{label}: expected a positive finite number")
    return measured


def _parse_finite_float(label: str, value: str) -> float:
    try:
        measured = float(value)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{label}: expected a finite number") from exc
    if not math.isfinite(measured):
        raise ContractError(f"{label}: expected a finite number")
    return measured


def _parse_positive_int(label: str, value: str) -> int:
    try:
        measured = int(value)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{label}: expected a positive integer") from exc
    if measured <= 0:
        raise ContractError(f"{label}: expected a positive integer")
    return measured


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_athinput(path: Path = _INPUT_DECK) -> dict[str, dict[str, str]]:
    """Parse the strict subset of Athena input syntax used by the frozen deck."""
    blocks: dict[str, dict[str, str]] = {}
    current = None
    for lineno, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<"):
            if not line.endswith(">"):
                raise ContractError(f"{path}:{lineno}: malformed block header")
            current = line[1:-1].strip()
            if not current:
                raise ContractError(f"{path}:{lineno}: empty block name")
            blocks.setdefault(current, {})
            continue
        if current is None or "=" not in line:
            raise ContractError(f"{path}:{lineno}: malformed parameter line")
        name, value = (item.strip() for item in line.split("=", 1))
        if not name or not value:
            raise ContractError(f"{path}:{lineno}: empty parameter name or value")
        if name in blocks[current]:
            raise ContractError(f"{path}:{lineno}: duplicate {current}/{name}")
        blocks[current][name] = value
    return blocks


def _require_exact_values(blocks: dict[str, dict[str, str]]) -> None:
    mismatches = []
    for (block, name), expected in _EXPECTED_VALUES.items():
        measured = blocks.get(block, {}).get(name)
        if measured != expected:
            mismatches.append(
                f"{block}/{name}: expected {expected!r}, measured {measured!r}"
            )
    if mismatches:
        raise ContractError("Q-011 deck contract mismatch:\n" + "\n".join(mismatches))


def _require_source_contract() -> None:
    required_fragments = {
        _PGEN_SOURCE: [
            "enum class PSShockSpeedModel { finite_mach, ideal_surface };",
            'pin->GetOrAddString(\n      "problem", "ps_shock_speed_model", "finite_mach")',
            "IdealSurfaceShockSpeed(gamma, ps_u0)",
            'pin->GetOrAddReal(\n      "problem", "ps_remove_birth_time_before", -1.0)',
            "const Real sweep_speed = ps_u0 + ps_shock_speed;",
            "const Real swept_mass = ps_eta*ps_rho0*sweep_speed*pm->dt*global_running_area;",
            "ps_particle_q_over_m = ps_particle_charge/ps_particle_mass;",
            "ps_particle_macro_mass = qscale*ps_particle_mass;",
            "h_pr_new(IPM, n) = ps_particle_q_over_m;",
            "h_pr_new(IPWT, n) = 1.0;",
            "xshock >= x1c - half_width && xshock < x1c + half_width",
            "part.x1 = xshock;",
            "RemoveExcludedEarlyInjectedParticles(pm);",
            "const Real pinj = ps_vinj_over_u0*ps_u0;",
            "const Real vinj = VelocityMagnitudeFromMomentumMagnitude(ppart, pinj);",
            "const Real mu = 2.0*TaggedUniform01(tag, 1) - 1.0;",
            "const Real phi = 2.0*M_PI*TaggedUniform01(tag, 2);",
            "BoostRelativeVelocityFromSurface(ppart, surface_vx, vinj*dirx, vinj*diry,",
        ],
        _PUSHER_SOURCE: [
            "Real q_over_m = pr(IPM, p);",
            "Real qdt_2m = q_over_m*dt_half;",
            "Real tx = qdt_2m*Bx*inv_gamma_minus;",
        ],
    }
    missing = []
    for path, fragments in required_fragments.items():
        source = path.read_text(encoding="utf-8")
        missing.extend(
            f"{path.relative_to(_REPO_ROOT)}: {fragment}"
            for fragment in fragments
            if fragment not in source
        )
    if missing:
        raise ContractError(
            "Q-011 source-local shock normalization contract mismatch:\n"
            + "\n".join(missing)
        )


def derive_normalization_and_macro_particle_calibration(
    blocks: dict[str, dict[str, str]],
) -> dict[str, object]:
    """Derive the paper-unit mapping and ideal downstream ppc calibration."""
    rho0 = _parse_positive_float("problem/ps_rho0", blocks["problem"]["ps_rho0"])
    b0 = _parse_positive_float("problem/ps_b0", blocks["problem"]["ps_b0"])
    u0 = _parse_positive_float("problem/ps_u0", blocks["problem"]["ps_u0"])
    eta = _parse_positive_float("problem/ps_eta", blocks["problem"]["ps_eta"])
    gamma = _parse_positive_float("mhd/gamma", blocks["mhd"]["gamma"])
    if gamma <= 1.0:
        raise ContractError("mhd/gamma: expected a value greater than one")
    numerical_light_speed = _parse_positive_float(
        "particles/pic_cr_light_speed",
        blocks["particles"]["pic_cr_light_speed"],
    )
    qscale = _parse_positive_float(
        "particles/deposit_qscale",
        blocks["particles"]["deposit_qscale"],
    )
    species_mass = _parse_positive_float(
        "species0/mass",
        blocks["species0"]["mass"],
    )
    species_charge = _parse_positive_float(
        "species0/charge",
        blocks["species0"]["charge"],
    )

    # Sun & Bai write Omega0 = B0*q/(m*c). AthenaK stores species charge/mass
    # directly in IPM and the Boris rotation uses IPM*B without another C
    # division. Under this deck convention, species charge/mass is therefore
    # the normalized paper q/(m*c), while pic_cr_light_speed is the separate
    # artificial numerical C used for relativistic kinematics. Dividing IPM by
    # pic_cr_light_speed again would double-apply the c normalization.
    normalized_q_over_mc = species_charge / species_mass
    alfven_speed_u_a0 = b0 / math.sqrt(rho0)
    cyclotron_frequency_omega0 = b0 * normalized_q_over_mc
    ion_inertial_length = 1.0 / (normalized_q_over_mc * math.sqrt(rho0))
    ion_inertial_length_via_u_a0_over_omega0 = (
        alfven_speed_u_a0 / cyclotron_frequency_omega0
    )
    rejected_raw_q_over_m_divided_by_artificial_c_omega0 = (
        species_charge * b0 / (species_mass * numerical_light_speed)
    )
    numerical_light_speed_over_u_a0 = numerical_light_speed / alfven_speed_u_a0

    nx1 = _parse_positive_int("mesh/nx1", blocks["mesh"]["nx1"])
    nx2 = _parse_positive_int("mesh/nx2", blocks["mesh"]["nx2"])
    nx3 = _parse_positive_int("mesh/nx3", blocks["mesh"]["nx3"])
    levels = _parse_positive_int(
        "mesh_refinement/num_levels",
        blocks["mesh_refinement"]["num_levels"],
    )
    root_dx = (
        _parse_finite_float("mesh/x1max", blocks["mesh"]["x1max"])
        - _parse_finite_float("mesh/x1min", blocks["mesh"]["x1min"])
    ) / nx1
    root_dy = (
        _parse_finite_float("mesh/x2max", blocks["mesh"]["x2max"])
        - _parse_finite_float("mesh/x2min", blocks["mesh"]["x2min"])
    ) / nx2
    collapsed_dz = (
        _parse_finite_float("mesh/x3max", blocks["mesh"]["x3max"])
        - _parse_finite_float("mesh/x3min", blocks["mesh"]["x3min"])
    ) / nx3
    for label, measured in (
        ("root dx", root_dx),
        ("root dy", root_dy),
        ("collapsed dz", collapsed_dz),
    ):
        if not math.isfinite(measured) or measured <= 0.0:
            raise ContractError(f"{label}: expected a positive finite cell size")

    cell_sizes = []
    for level in range(levels):
        refinement = 2**level
        dx = root_dx / refinement
        dy = root_dy / refinement
        cell_sizes.append(
            {
                "level": level,
                "dx_c_over_omega_pi": dx / ion_inertial_length,
                "dy_c_over_omega_pi": dy / ion_inertial_length,
                "collapsed_dz_c_over_omega_pi": collapsed_dz / ion_inertial_length,
                "effective_2d_cell_volume": dx * dy * collapsed_dz,
            }
        )

    ideal_shock_speed = 0.5 * (gamma - 1.0) * u0
    swept_speed = u0 + ideal_shock_speed
    macro_particle_mass = qscale * species_mass
    downstream_compression = swept_speed / ideal_shock_speed
    downstream_macro_particle_density = (
        eta * rho0 * downstream_compression / macro_particle_mass
    )
    expected_downstream_ppc = [
        downstream_macro_particle_density * item["effective_2d_cell_volume"]
        for item in cell_sizes
    ]
    macro_mass_from_coarse_ppc = (
        eta
        * rho0
        * downstream_compression
        * cell_sizes[0]["effective_2d_cell_volume"]
        / 640.0
    )
    macro_mass_from_fine_ppc = (
        eta
        * rho0
        * downstream_compression
        * cell_sizes[-1]["effective_2d_cell_volume"]
        / 40.0
    )

    _require_close("U_A0", alfven_speed_u_a0, 1.0)
    _require_close("normalized q/(m*c)", normalized_q_over_mc, 1.0)
    _require_close("Omega0", cyclotron_frequency_omega0, 1.0)
    _require_close(
        "c/omega_pi identity",
        ion_inertial_length,
        ion_inertial_length_via_u_a0_over_omega0,
    )
    _require_close("c/omega_pi", ion_inertial_length, 1.0)
    _require_close("C/U_A0", numerical_light_speed_over_u_a0, 10000.0)
    if len(cell_sizes) != 3:
        raise ContractError("AMR ladder: expected exactly three levels")
    for level, (item, expected_size) in enumerate(
        zip(cell_sizes, (12.0, 6.0, 3.0))
    ):
        _require_close(f"level {level} dx", item["dx_c_over_omega_pi"], expected_size)
        _require_close(f"level {level} dy", item["dy_c_over_omega_pi"], expected_size)
    _require_close("collapsed dz", collapsed_dz, 1.0)
    _require_close("ideal shock speed/U_A0", ideal_shock_speed / alfven_speed_u_a0, 10.0)
    _require_close("swept speed/U_A0", swept_speed / alfven_speed_u_a0, 40.0)
    _require_close("qscale to macro mass", macro_particle_mass, 9.0e-4)
    _require_close("ideal downstream compression", downstream_compression, 4.0)
    for level, (measured, expected) in enumerate(
        zip(expected_downstream_ppc, (640.0, 160.0, 40.0))
    ):
        _require_close(f"level {level} downstream ppc", measured, expected)
    _require_close(
        "coarse reverse-calibrated macro mass",
        macro_mass_from_coarse_ppc,
        macro_particle_mass,
    )
    _require_close(
        "fine reverse-calibrated macro mass",
        macro_mass_from_fine_ppc,
        macro_particle_mass,
    )

    return {
        "convention_ambiguity": {
            "status": "resolved_for_the_frozen_athenak_deck",
            "paper_formula": "Omega0 = B0 * q / (m * c)",
            "athenak_deck_mapping": (
                "species charge/mass is the normalized q/(m*c) coefficient "
                "stored in IPM; pic_cr_light_speed is the separate artificial C"
            ),
            "applicable_formula": "Omega0 = B0 * (species_charge / species_mass)",
            "rejected_double_division_formula": (
                "Omega0 = B0 * species_charge / "
                "(species_mass * pic_cr_light_speed)"
            ),
            "rejected_double_division_omega0": (
                rejected_raw_q_over_m_divided_by_artificial_c_omega0
            ),
        },
        "formulas": {
            "alfven_speed_u_a0": "B0 / sqrt(rho0)",
            "cyclotron_frequency_omega0": "B0 * normalized_q_over_mc",
            "ion_inertial_length_c_over_omega_pi": (
                "1 / (normalized_q_over_mc * sqrt(rho0)) = U_A0 / Omega0"
            ),
            "numerical_light_speed_over_u_a0": "pic_cr_light_speed / U_A0",
            "ideal_shock_speed": "(gamma - 1) * u0 / 2",
            "upstream_relative_swept_speed": "u0 + ideal_shock_speed",
            "macro_particle_mass": "deposit_qscale * species_mass * IPWT",
            "ideal_downstream_ppc": (
                "eta * rho0 * swept_speed / ideal_shock_speed "
                "* effective_2d_cell_volume / macro_particle_mass"
            ),
        },
        "alfven_speed_u_a0": alfven_speed_u_a0,
        "normalized_q_over_mc": normalized_q_over_mc,
        "cyclotron_frequency_omega0": cyclotron_frequency_omega0,
        "ion_inertial_length_c_over_omega_pi": ion_inertial_length,
        "ion_inertial_length_via_u_a0_over_omega0": (
            ion_inertial_length_via_u_a0_over_omega0
        ),
        "numerical_light_speed_c": numerical_light_speed,
        "numerical_light_speed_over_u_a0": numerical_light_speed_over_u_a0,
        "cell_sizes": cell_sizes,
        "ideal_shock_speed_over_u_a0": ideal_shock_speed / alfven_speed_u_a0,
        "upstream_relative_swept_speed_over_u_a0": swept_speed / alfven_speed_u_a0,
        "macro_particle": {
            "formula": "deposit_qscale * species_mass * IPWT",
            "injected_ipwt": 1.0,
            "deposit_qscale": qscale,
            "species_mass": species_mass,
            "mass": macro_particle_mass,
        },
        "ideal_downstream_calibration": {
            "scope": (
                "ideal mean immediately-downstream injection calibration; "
                "not a campaign measurement"
            ),
            "effective_dimension": "2d_with_collapsed_x3_thickness",
            "compression_ratio": downstream_compression,
            "macro_particle_density": downstream_macro_particle_density,
            "target_ppc_by_level": [640.0, 160.0, 40.0],
            "expected_ppc_by_level": expected_downstream_ppc,
            "expected_coarse_ppc": expected_downstream_ppc[0],
            "expected_fine_ppc": expected_downstream_ppc[-1],
            "macro_mass_from_coarse_ppc": macro_mass_from_coarse_ppc,
            "macro_mass_from_fine_ppc": macro_mass_from_fine_ppc,
        },
    }


def _shock_surface_derivation(
    blocks: dict[str, dict[str, str]],
    calibration: dict[str, object],
) -> dict[str, object]:
    gamma = float(blocks["mhd"]["gamma"])
    u0 = float(blocks["problem"]["ps_u0"])
    rho0 = float(blocks["problem"]["ps_rho0"])
    p0 = float(blocks["problem"]["ps_p0"])
    ideal_surface_speed = 0.5 * (gamma - 1.0) * u0
    sound_speed_squared = gamma * p0 / rho0
    mach_squared = u0 * u0 / sound_speed_squared
    compression_ratio = (
        (gamma + 1.0) * mach_squared
        / ((gamma - 1.0) * mach_squared + 2.0)
    )
    finite_mach_speed = u0 / (compression_ratio - 1.0)
    upstream_relative_sweep_speed = u0 + ideal_surface_speed
    alfven_speed = calibration["alfven_speed_u_a0"]
    ion_inertial_length = calibration["ion_inertial_length_c_over_omega_pi"]
    return {
        "selected_model": "ideal_surface",
        "paper_equation": "u_sh_prime = (Gamma - 1) * u0 / 2",
        "ideal_surface_speed_over_ua0": ideal_surface_speed / alfven_speed,
        "finite_mach_engineering_option_speed_over_ua0": (
            finite_mach_speed / alfven_speed
        ),
        "upstream_relative_sweep_speed_over_ua0": (
            upstream_relative_sweep_speed / alfven_speed
        ),
        "ideal_surface_positions_c_over_omega_pi": {
            "t500": ideal_surface_speed * 500.0 / ion_inertial_length,
            "t1200": ideal_surface_speed * 1200.0 / ion_inertial_length,
        },
        "injection_distribution": "monoenergetic_full_sphere_isotropic_relative_to_ideal_surface",
        "shock_surface_carrier_selection": "single_half_open_cell_with_surface_x1",
        "early_injected_particle_removal": "runtime_state_removal_for_birth_time_below_45",
    }


def validate_deck(path: Path = _INPUT_DECK) -> dict[str, object]:
    """Validate sourced values and explicit AthenaK preparation choices."""
    blocks = parse_athinput(path)
    _require_exact_values(blocks)
    _require_source_contract()
    calibration = derive_normalization_and_macro_particle_calibration(blocks)
    cell_sizes = [
        item["dx_c_over_omega_pi"] for item in calibration["cell_sizes"]
    ]
    ua0 = calibration["alfven_speed_u_a0"]
    mach_alfven = float(blocks["problem"]["ps_u0"]) / ua0
    light_speed_over_ua0 = calibration["numerical_light_speed_over_u_a0"]
    if mach_alfven != 30.0 or light_speed_over_ua0 != 10000.0:
        raise ContractError("Q-011 dimensionless normalization contract mismatch")
    return {
        "deck": str(path.relative_to(_REPO_ROOT)),
        "deck_sha256": _sha256(path),
        "root_cell_size_c_over_omega_pi": cell_sizes[0],
        "amr_cell_sizes_c_over_omega_pi": cell_sizes,
        "mach_alfven": mach_alfven,
        "light_speed_over_ua0": light_speed_over_ua0,
        "normalization_and_macro_particle_calibration": calibration,
        "shock_surface_derivation": _shock_surface_derivation(blocks, calibration),
        "pgen_source_sha256": _sha256(_PGEN_SOURCE),
        "particles_pusher_source_sha256": _sha256(_PUSHER_SOURCE),
    }


def build_preparation_contract() -> dict[str, object]:
    """Return the frozen contract without inspecting campaign output."""
    return {
        "schema_version": 2,
        "gate": "Q-011",
        "claim_id": "CLAIM-PAPER-SHOCK-001",
        "qualification_effect": "source_controlled_preparation_only",
        "deck": validate_deck(),
        "paper_text_values": {
            "geometry": "reflecting_left_wall_parallel_Bx_periodic_y",
            "domain_c_over_omega_pi": [48000.0, 3120.0],
            "dimensions": "2D3V",
            "mach_alfven": 30.0,
            "gas_gamma": 5.0 / 3.0,
            "injection_efficiency": 1.0e-3,
            "injection_p_over_m_over_u0": 10.0**0.5,
            "light_speed_over_ua0": 10000.0,
            "exclude_birth_time_before_omega0_inverse": 45.0,
            "amr_cell_sizes_c_over_omega_pi": [12.0, 6.0, 3.0],
            "expected_downstream_ppc": {"coarse": 640.0, "fine": 40.0},
            "amr_refine_curvature_threshold": 1.0,
            "amr_derefine_curvature_threshold": 0.1,
            "snapshot_times_omega0_inverse": [500.0, 1200.0],
            "shock_surface_model": "ideal_surface",
            "injection_distribution": (
                "monoenergetic_full_sphere_isotropic_relative_to_ideal_surface"
            ),
        },
        "artifact_manifest_contract": {
            "required_candidate_bindings": [
                "clean_candidate_git_commit",
                "clean_frontier_executable_sha256",
                "input_deck_sha256",
                "analyzer_sha256",
                "authorized_orion_artifact_root",
            ],
            "required_seeds": [
                23050101,
                23050102,
                23050103,
                23050104,
                23050105,
                23050106,
                23050107,
                23050108,
            ],
            "required_grid_variants": [
                "coarse_uniform",
                "three_level_amr_root_12_finest_3",
                "fine_uniform",
            ],
            "required_snapshot_times_omega0_inverse": [500.0, 1200.0],
            "required_raw_artifacts_per_snapshot": [
                "rho_bin",
                "bmag_bin",
                "j2_bin",
                "prtcl_all_pvtk",
            ],
            "required_run_artifacts": [
                "stdout_with_q017_telemetry",
                "restart_checkpoint",
                "attempt_status_and_failure_artifacts",
            ],
            "particle_filters": {
                "source": "shock_injected",
                "birth_time_min_omega0_inverse": 45.0,
                "region": "downstream",
            },
            "primary_observables": [
                "shock_position",
                "upstream_magnetic_amplification_at_t500",
                "downstream_weighted_energy_spectra_at_t500_and_t1200",
                "late_energy_spectrum_power_law_slope_at_t1200",
                "morphology",
            ],
        },
        "campaign_execution": {
            "platform": "Frontier",
            "bulk_artifact_root": "/lustre/orion/ast207/proj-shared/dfielding/PIC",
            "status": "open_not_executed_by_preparation_tranche",
            "forbidden_storage_systems": ["Kronos"],
        },
        "open_items": list(_OPEN_ITEMS),
        "result_metrics": [],
    }


def analyze_campaign_artifacts(_manifest: dict[str, object]) -> None:
    """Reject output inspection until the remaining preregistration fields close."""
    raise ContractError(
        "Q-011 campaign analysis is blocked: Frontier execution and pre-run "
        "contract fields remain open"
    )


def run(**kwargs) -> None:
    """Run the local static deck and fail-closed analyzer-contract regression."""
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    contract = build_preparation_contract()
    _RESULTS["deck_contract"] = bool(contract["deck"]["deck_sha256"])
    _RESULTS["frontier_campaign_execution_open"] = (
        contract["campaign_execution"]["status"]
        == "open_not_executed_by_preparation_tranche"
    )
    try:
        analyze_campaign_artifacts({})
    except ContractError:
        _RESULTS["artifact_analysis_fail_closed"] = True
    else:
        _RESULTS["artifact_analysis_fail_closed"] = False


def analyze() -> bool:
    """Report only preparation-regression status, never fabricated science results."""
    logger.info("Q-011 Section 5.4 preparation contract: %s", _RESULTS)
    return _RESULTS == {
        "deck_contract": True,
        "frontier_campaign_execution_open": True,
        "artifact_analysis_fail_closed": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-contract", action="store_true")
    parser.add_argument("--artifact-manifest", type=Path)
    args = parser.parse_args()
    contract = build_preparation_contract()
    if args.artifact_manifest is not None:
        manifest = json.loads(args.artifact_manifest.read_text(encoding="utf-8"))
        analyze_campaign_artifacts(manifest)
    print(json.dumps(contract, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
