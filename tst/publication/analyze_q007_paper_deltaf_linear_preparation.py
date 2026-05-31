#!/usr/bin/env python3
"""Nonqualifying Q-007 Sun-Bai true-delta-f linear preparation audit."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
QUALIFICATION_EFFECT = "none_source_local_true_deltaf_preparation_only"
ARTIFACT_ROLE = "section55_section56_true_deltaf_preparation_not_qualifying_evidence"
PGEN_METHODS = {
    "crsi": "Q007PaperCRSILinearPreparation",
    "crpai": "Q007PaperCRPAILinearPreparation",
}
DECKS = {
    "crsi": REPO_ROOT / "inputs/tests/pic_q007_paper_crsi_linear_preparation.athinput",
    "crpai_prolate": (
        REPO_ROOT
        / "inputs/tests/pic_q007_paper_crpai_linear_prolate_preparation.athinput"
    ),
    "crpai_oblate": (
        REPO_ROOT
        / "inputs/tests/pic_q007_paper_crpai_linear_oblate_preparation.athinput"
    ),
}
SOURCE = REPO_ROOT / "src/pgen/tests/q007_paper_deltaf_linear.cpp"
HEADER = REPO_ROOT / "src/pgen/tests/q007_paper_deltaf_linear.hpp"
PARTICLES = REPO_ROOT / "src/particles/particles.cpp"
CMAKE = REPO_ROOT / "src/CMakeLists.txt"
PGEN_HEADER = REPO_ROOT / "src/pgen/pgen.hpp"
PGEN_DISPATCH = REPO_ROOT / "src/pgen/pgen.cpp"

_COMMON = {
    ("mesh", "nx1"): "32",
    ("mesh", "nx2"): "4",
    ("mesh", "nx3"): "1",
    ("meshblock", "nx1"): "32",
    ("meshblock", "nx2"): "4",
    ("meshblock", "nx3"): "1",
    ("mesh_refinement", "refinement"): "none",
    ("time", "evolution"): "dynamic",
    ("time", "integrator"): "rk2",
    ("time", "cfl_number"): "0.1",
    ("time", "nlim"): "0",
    ("time", "tlim"): "0.0",
    ("mhd", "eos"): "isothermal",
    ("mhd", "iso_sound_speed"): "1.0",
    ("mhd", "reconstruct"): "plm",
    ("mhd", "rsolver"): "llf",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "ppc"): "0.0625",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "nspecies"): "8",
    ("particles", "cr_distribution"): "center",
    ("particles", "deposit_moments"): "true",
    ("particles", "deposit_order"): "1",
    ("particles", "deposit_qscale"): "1.0",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "couple_j_to_efield_coeff"): "1.0",
    ("particles", "couple_j_to_efield_representation"): "cell_centered",
    ("particles", "couple_j_deposition_mode"): "cc_convert",
    ("particles", "couple_moments_momentum_to_mhd"): "true",
    ("particles", "couple_moments_energy_to_mhd"): "false",
    ("particles", "couple_fluid_feedback_order"): "mhd_src_terms",
    ("particles", "pic_physical_mode"): "paper_mhd_pic",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_interp_scheme"): "tsc",
    ("particles", "pic_enable_2d3v"): "true",
    ("particles", "pic_cr_initial_state"): "momentum",
    ("particles", "pic_cr_hall_mode"): "off",
    ("particles", "pic_wave_damping_mode"): "off",
    ("particles", "pic_deltaf_mode"): "physical",
    ("particles", "pic_deltaf_background_rho"): "1.0e-4",
    ("particles", "pic_deltaf_background_jx"): "0.0",
    ("particles", "pic_deltaf_background_jy"): "0.0",
    ("particles", "pic_deltaf_background_jz"): "0.0",
    ("particles", "pic_deltaf_adapt_mode"): "off",
    ("particles", "pic_expanding_box_mode"): "off",
}
_METADATA_COMMON = {
    "deck_role": "cycle_zero_true_deltaf_parser_and_analytic_mapping_preparation_only",
    "qualification_effect": "none",
    "frontier_authorization": "not_bound",
    "runtime_evolution": "blocked_cycle_zero_only",
    "physical_loading": "blocked_missing_paper_log_bin_weighted_loader",
    "initial_wave_spectrum": "blocked_missing_random_phase_four_branch_loader",
    "theory_runtime_comparison":
        "blocked_missing_runtime_artifacts_and_independent_review",
    "paper_momentum_bin_count": "8",
    "paper_particles_per_cell_per_bin": "256",
    "source_local_placeholder_particles_total": "8",
    "paper_rho0": "1.0",
    "paper_b0": "1.0",
    "paper_ua": "1.0",
    "paper_mncr_over_rho0": "1.0e-4",
    "paper_domain_x1": "96000.0",
}
_CASES = {
    "crsi": {
        "block": "q007_paper_crsi_linear_preparation",
        "pgen_name": "q007_paper_crsi_linear_preparation",
        "campaign_id": "Q007-PAPER-CRSI-LINEAR-PREPARATION",
        "section_anchor": "sun_bai_2023_section_5_5_1_crsi",
        "x1max": "320.0",
        "paper_dx": "10.0",
        "light_speed": "300.0",
        "p0": "300.0",
        "kappa": "1.25",
        "xi": "1.0",
        "f0": "kappa_iso",
        "aniso_transverse": "1.0",
        "case_role": "crsi_isotropic_kappa_gas_drift_minus_vd",
        "paper_vd": "2.0",
        "source_local_gas_vx": "-2.0",
        "handedness_mapping":
            "paper_both_forward_polarizations_static_mapping_only",
    },
    "crpai_prolate": {
        "block": "q007_paper_crpai_linear_preparation",
        "pgen_name": "q007_paper_crpai_linear_preparation",
        "campaign_id": "Q007-PAPER-CRPAI-LINEAR-PREPARATION",
        "section_anchor": "sun_bai_2023_section_5_5_2_crpai",
        "x1max": "640.0",
        "paper_dx": "20.0",
        "light_speed": "30000.0",
        "p0": "300.0",
        "kappa": "1.75",
        "xi": "0.99",
        "f0": "kappa_aniso",
        "aniso_transverse": "1.0101010101010102",
        "case_role": "crpai_prolate_xi_0p99_signed_branch_mapping_only",
        "source_local_gas_vx": "0.0",
        "handedness_mapping": "blocked_pending_manuscript_text_caption_review",
    },
    "crpai_oblate": {
        "block": "q007_paper_crpai_linear_preparation",
        "pgen_name": "q007_paper_crpai_linear_preparation",
        "campaign_id": "Q007-PAPER-CRPAI-LINEAR-PREPARATION",
        "section_anchor": "sun_bai_2023_section_5_5_2_crpai",
        "x1max": "640.0",
        "paper_dx": "20.0",
        "light_speed": "30000.0",
        "p0": "300.0",
        "kappa": "1.75",
        "xi": "1.01",
        "f0": "kappa_aniso",
        "aniso_transverse": "0.9900990099009901",
        "case_role": "crpai_oblate_xi_1p01_signed_branch_mapping_only",
        "source_local_gas_vx": "0.0",
        "handedness_mapping": "blocked_pending_manuscript_text_caption_review",
    },
}


class ContractError(ValueError):
    """Raised when a bounded Q-007 source-local contract fails closed."""


def parse_athinput(path: Path) -> dict[str, dict[str, str]]:
    """Parse the strict Athena input subset used by the Q-007 preparation decks."""
    blocks: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    for lineno, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<"):
            if not line.endswith(">"):
                raise ContractError(f"{path}:{lineno}: malformed block header")
            name = line[1:-1].strip()
            if not name:
                raise ContractError(f"{path}:{lineno}: empty block name")
            current = blocks.setdefault(name, {})
            continue
        if current is None or "=" not in line:
            raise ContractError(f"{path}:{lineno}: malformed parameter line")
        name, value = (item.strip() for item in line.split("=", 1))
        if not name or not value:
            raise ContractError(f"{path}:{lineno}: empty parameter name or value")
        if name in current:
            raise ContractError(f"{path}:{lineno}: duplicate parameter {name}")
        current[name] = value
    return blocks


def _finite(label: str, value: float) -> float:
    measured = float(value)
    if not math.isfinite(measured):
        raise ContractError(f"{label} must be finite")
    return measured


def _positive(label: str, value: float) -> float:
    measured = _finite(label, value)
    if measured <= 0.0:
        raise ContractError(f"{label} must be positive")
    return measured


def kappa_normalization(number_density: float, p0: float, kappa: float) -> float:
    """Return the paper's normalized isotropic kappa prefactor."""
    density = _positive("number density", number_density)
    scale = _positive("p0", p0)
    index = _positive("kappa", kappa)
    if index <= 0.5:
        raise ContractError("kappa must exceed 0.5")
    return (
        density
        * math.gamma(index + 1.0)
        / ((math.pi * index * scale * scale) ** 1.5 * math.gamma(index - 0.5))
    )


def isotropic_kappa_distribution(
    number_density: float, momentum: float, p0: float, kappa: float
) -> float:
    """Return Sun-Bai equation kappa_dist_iso."""
    shape = 1.0 + _finite("momentum", momentum) ** 2 / (kappa * p0 * p0)
    return kappa_normalization(number_density, p0, kappa) * shape ** (-kappa - 1.0)


def anisotropic_kappa_distribution(
    number_density: float,
    px: float,
    py: float,
    pz: float,
    p0: float,
    kappa: float,
    xi: float,
) -> float:
    """Return Sun-Bai equation kappa_dist_aniso in Cartesian momentum form."""
    anisotropy = _positive("xi", xi)
    xi2 = anisotropy * anisotropy
    shape = 1.0 + (
        _finite("px", px) ** 2
        + xi2 * (_finite("py", py) ** 2 + _finite("pz", pz) ** 2)
    ) / (kappa * p0 * p0)
    return (
        xi2
        * kappa_normalization(number_density, p0, kappa)
        * shape ** (-kappa - 1.0)
    )


def resonant_wavenumber(mass: float, omega0: float, p0: float) -> float:
    """Return the paper-literal resonant peak scale k0=m*Omega0/p0."""
    return _positive("mass", mass) * _positive("omega0", omega0) / _positive("p0", p0)


def gyroresonance_q2(
    k: float, mass: float, omega0: float, p0: float, kappa: float
) -> float:
    """Return the paper's closed Q2 gyroresonance factor."""
    wave_number = abs(_finite("k", k))
    if wave_number == 0.0:
        raise ContractError("k must be nonzero")
    index = _positive("kappa", kappa)
    ratio = _positive("mass", mass) * _positive("omega0", omega0)
    ratio /= wave_number * _positive("p0", p0)
    prefactor = (
        math.sqrt(math.pi)
        / index**1.5
        * math.gamma(index + 1.0)
        / math.gamma(index - 0.5)
    )
    return prefactor * ratio * (1.0 + ratio * ratio / index) ** (-index)


def crsi_low_density_growth_rate(
    k: float,
    mass: float,
    omega0: float,
    p0: float,
    kappa: float,
    mncr_over_rho0: float,
    vd_over_ua: float,
) -> float:
    """Return the paper's low-density CRSI approximation for one signed k."""
    direction = 1.0 if _finite("k", k) >= 0.0 else -1.0
    return (
        -0.5
        * _positive("mncr/rho0", mncr_over_rho0)
        * _positive("omega0", omega0)
        * (1.0 - _positive("vd/ua", vd_over_ua) * direction)
        * gyroresonance_q2(k, mass, omega0, p0, kappa)
    )


def crpai_low_density_growth_rate(
    k: float,
    mass: float,
    omega0: float,
    p0: float,
    kappa: float,
    mncr_over_rho0: float,
    ua: float,
    xi: float,
    signed_branch: int,
) -> float:
    """Return the paper's CRPAI approximation without assigning handedness labels."""
    if signed_branch not in (-1, 1):
        raise ContractError("signed branch must be -1 or +1")
    wave_number = abs(_finite("k", k))
    if wave_number == 0.0:
        raise ContractError("k must be nonzero")
    anisotropy = _positive("xi", xi)
    xi2 = anisotropy * anisotropy
    branch_factor = (
        1.0
        + signed_branch
        * (1.0 - xi2)
        / xi2
        * _positive("omega0", omega0)
        / (wave_number * _positive("ua", ua))
    )
    return (
        -0.5
        * _positive("mncr/rho0", mncr_over_rho0)
        * omega0
        * branch_factor
        * gyroresonance_q2(k, mass, omega0, p0, kappa)
        / xi2
    )


def validate_deck(path: Path, case: str) -> dict[str, Any]:
    """Require a cycle-zero parser carrier and paper-literal static metadata."""
    expected_case = _CASES[case]
    blocks = parse_athinput(path)
    for (block, name), expected in _COMMON.items():
        measured = blocks.get(block, {}).get(name)
        if measured != expected:
            raise ContractError(f"{path}: {block}/{name}: expected {expected!r}")
    for name in ("ix1_bc", "ox1_bc", "ix2_bc", "ox2_bc"):
        if blocks["mesh"].get(name) != "periodic":
            raise ContractError(f"{path}: mesh/{name} must remain periodic")
    if blocks["mesh"].get("x1max") != expected_case["x1max"]:
        raise ContractError(f"{path}: mesh/x1max does not match bounded carrier")
    if blocks["problem"].get("pgen_name") != expected_case["pgen_name"]:
        raise ContractError(f"{path}: problem/pgen_name does not match case")
    particles = blocks["particles"]
    for name, expected in (
        ("pic_cr_light_speed", expected_case["light_speed"]),
        ("pic_deltaf_f0", expected_case["f0"]),
        ("pic_deltaf_p0", expected_case["p0"]),
        ("pic_deltaf_kappa", expected_case["kappa"]),
        ("pic_deltaf_aniso_x1", "1.0"),
        ("pic_deltaf_aniso_x2", expected_case["aniso_transverse"]),
        ("pic_deltaf_aniso_x3", expected_case["aniso_transverse"]),
    ):
        if particles.get(name) != expected:
            raise ContractError(f"{path}: particles/{name}: expected {expected!r}")
    for species in range(8):
        block = blocks.get(f"species{species}", {})
        if block != {
            "mass": "1.0",
            "charge": "1.0",
            "vx0": "0.0",
            "vy0": "0.0",
            "vz0": "0.0",
        }:
            raise ContractError(f"{path}: species{species} placeholder contract drifted")
    metadata = blocks[expected_case["block"]]
    for name, expected in _METADATA_COMMON.items():
        if metadata.get(name) != expected:
            raise ContractError(f"{path}: {expected_case['block']}/{name} drifted")
    for name, expected_name in (
        ("campaign_id", "campaign_id"),
        ("section_anchor", "section_anchor"),
        ("paper_dx", "paper_dx"),
        ("paper_light_speed", "light_speed"),
        ("paper_p0", "p0"),
        ("paper_kappa", "kappa"),
        ("paper_xi", "xi"),
        ("case_role", "case_role"),
        ("source_local_gas_vx", "source_local_gas_vx"),
        ("handedness_mapping", "handedness_mapping"),
    ):
        if metadata.get(name) != expected_case[expected_name]:
            raise ContractError(f"{path}: {expected_case['block']}/{name} drifted")
    if case == "crsi" and metadata.get("paper_vd") != expected_case["paper_vd"]:
        raise ContractError(f"{path}: CRSI paper drift mapping drifted")
    placeholder_count = float(particles["ppc"]) * 32 * 4
    if placeholder_count != 8.0:
        raise ContractError(f"{path}: placeholder count must remain exactly eight")
    return {
        "case": case,
        "path": str(path.relative_to(REPO_ROOT)),
        "pgen_name": expected_case["pgen_name"],
        "cycle_zero_only": True,
        "true_deltaf_parser_path": True,
        "exact_isothermal_mhd": True,
        "source_local_placeholder_particles_total": int(placeholder_count),
        "paper_particles_per_cell_total": 8 * 256,
        "physical_loading_implemented": False,
        "random_phase_wave_spectrum_implemented": False,
        "runtime_evolution_admitted": False,
        "qualifying_evidence": False,
    }


def validate_decks() -> list[dict[str, Any]]:
    """Validate all separately named CRSI and CRPAI source-local decks."""
    return [validate_deck(path, case) for case, path in DECKS.items()]


def validate_source_contract() -> dict[str, Any]:
    """Require additive registration and the narrow parser allowance."""
    source = SOURCE.read_text(encoding="utf-8")
    header = HEADER.read_text(encoding="utf-8")
    particles = PARTICLES.read_text(encoding="utf-8")
    cmake = CMAKE.read_text(encoding="utf-8")
    declarations = PGEN_HEADER.read_text(encoding="utf-8")
    dispatch = PGEN_DISPATCH.read_text(encoding="utf-8")
    if "pgen/tests/q007_paper_deltaf_linear.cpp" not in cmake:
        raise ContractError("Q-007 source is absent from CMake")
    for pgen_name, method in (
        ("q007_paper_crsi_linear_preparation", PGEN_METHODS["crsi"]),
        ("q007_paper_crpai_linear_preparation", PGEN_METHODS["crpai"]),
    ):
        if dispatch.count(f'compare("{pgen_name}")') != 2:
            raise ContractError(f"{pgen_name} fresh/restart dispatch is incomplete")
        if f"void {method}(ParameterInput *pin, const bool restart);" not in declarations:
            raise ContractError(f"{method} declaration is absent")
        if f"ProblemGenerator::{method}" not in source:
            raise ContractError(f"{method} implementation is absent")
    for snippet in (
        "exact_isothermal_deltaf_paper_feedback",
        "UsesDeltaF() && !couple_moments_energy_to_mhd",
        "exact isothermal paper delta-f uses momentum-only",
    ):
        if snippet not in particles:
            raise ContractError("narrow exact-isothermal parser allowance drifted")
    for snippet in (
        "blocked_missing_paper_log_bin_weighted_loader",
        "blocked_missing_random_phase_four_branch_loader",
        "blocked_pending_manuscript_text_caption_review",
        "Q007RequireInteger(pin, \"time\", \"nlim\", 0);",
    ):
        if snippet not in source:
            raise ContractError(f"Q-007 source is missing fail-closed marker {snippet!r}")
    for snippet in (
        "IsotropicKappaDistribution",
        "AnisotropicKappaDistribution",
        "GyroresonanceQ2",
        "CRSILowDensityGrowthRate",
        "CRPAILowDensityGrowthRate",
    ):
        if snippet not in header:
            raise ContractError(f"Q-007 header is missing analytical helper {snippet!r}")
    return {
        "compilation_unit_registered": True,
        "fresh_and_restart_dispatch_registered": True,
        "narrow_exact_isothermal_true_deltaf_parser_allowance": True,
        "cycle_zero_guarded": True,
    }


def analytical_contract() -> dict[str, Any]:
    """Build static paper mappings without claiming a runtime theory oracle."""
    k_crsi = resonant_wavenumber(1.0, 1.0, 300.0)
    q2_crsi = gyroresonance_q2(k_crsi, 1.0, 1.0, 300.0, 1.25)
    crsi_forward = crsi_low_density_growth_rate(
        k_crsi, 1.0, 1.0, 300.0, 1.25, 1.0e-4, 2.0
    )
    crsi_backward = crsi_low_density_growth_rate(
        -k_crsi, 1.0, 1.0, 300.0, 1.25, 1.0e-4, 2.0
    )
    crpai = {}
    for role, xi in (("prolate", 0.99), ("oblate", 1.01)):
        growth = {
            str(branch): crpai_low_density_growth_rate(
                k_crsi, 1.0, 1.0, 300.0, 1.75, 1.0e-4, 1.0, xi, branch
            )
            for branch in (-1, 1)
        }
        crpai[role] = {
            "xi": xi,
            "athenak_transverse_anisotropy_scale": 1.0 / xi,
            "signed_branch_growth": growth,
            "unstable_signed_branch": max(growth, key=growth.get),
            "handedness_mapping":
                "blocked_pending_manuscript_text_caption_review",
        }
    return {
        "paper_reference": "sun_bai_2023_arxiv_2304.10568v1",
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "full_q1_dispersion_runtime_oracle": False,
        "crsi": {
            "k0": k_crsi,
            "lambda0": 2.0 * math.pi / k_crsi,
            "q2_at_k0": q2_crsi,
            "forward_growth_at_k0": crsi_forward,
            "backward_growth_at_minus_k0": crsi_backward,
        },
        "crpai": crpai,
        "crpai_handedness_claimed": False,
        "crpai_handedness_boundary":
            "manuscript prose and figure caption require independent review",
    }


def build_preparation_report() -> dict[str, Any]:
    """Return the bounded Q-007 source-local preparation report."""
    return {
        "schema_version": 1,
        "gate": "Q-007",
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "claim_closure": False,
        "frontier_authorization": False,
        "source_contract": validate_source_contract(),
        "decks": validate_decks(),
        "analytical_contract": analytical_contract(),
        "nonqualification_boundary": {
            "deck_source_freeze": True,
            "true_deltaf_parser_path": True,
            "exact_isothermal_mhd_cycle_zero_startup": True,
            "paper_log_bin_weighted_loading": False,
            "paper_random_phase_four_branch_wave_spectrum": False,
            "runtime_evolution": False,
            "full_q1_dispersion_runtime_oracle": False,
            "crpai_handedness_mapping_reviewed": False,
            "mpi_qualification": False,
            "gpu_qualification": False,
            "frontier_authorization": False,
            "external_review": False,
            "qualifying_evidence": False,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_preparation_report()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
