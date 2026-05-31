#!/usr/bin/env python3
"""Nonqualifying Q-006 Sun-Bai Section 5.3 source-local preparation audit."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ID = "Q006-PAPER-MULTISPECIES-OSCILLATION"
PGEN_NAME = "q006_paper_multispecies_oscillation"
PGEN_METHOD = "Q006PaperMultispeciesOscillation"
QUALIFICATION_EFFECT = "none_source_local_deck_and_source_freeze_preparation_only"
ARTIFACT_ROLE = "section53_source_local_preparation_not_qualifying_evidence"
AMR_SEED = 60053
AMR_REFINE_PROBABILITY = 0.10
AMR_DEREFINE_PROBABILITY = 0.60

DECKS = {
    "uniform": (
        REPO_ROOT
        / "inputs/tests/pic_q006_paper_multispecies_oscillation_uniform_candidate.athinput"
    ),
    "smr": (
        REPO_ROOT
        / "inputs/tests/pic_q006_paper_multispecies_oscillation_smr_candidate.athinput"
    ),
    "audited_amr_preparation": (
        REPO_ROOT
        / "inputs/tests/"
        "pic_q006_paper_multispecies_oscillation_audited_amr_candidate.athinput"
    ),
}

SOURCE = REPO_ROOT / "src/pgen/tests/q006_paper_multispecies_oscillation.cpp"

_EXPECTED_COMMON = {
    ("time", "integrator"): "rk2",
    ("time", "cfl_number"): "0.1",
    ("time", "nlim"): "0",
    ("time", "tlim"): "0.0",
    ("mhd", "eos"): "ideal",
    ("mhd", "reconstruct"): "plm",
    ("mhd", "rsolver"): "llf",
    ("mhd", "gamma"): "1.66666666667",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "ppc"): "128.0",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "nspecies"): "2",
    ("particles", "cr_distribution"): "center",
    ("particles", "deposit_moments"): "true",
    ("particles", "deposit_order"): "1",
    ("particles", "deposit_qscale"): "0.0234375",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "couple_j_to_efield_coeff"): "1.0",
    ("particles", "couple_j_to_efield_representation"): "cell_centered",
    ("particles", "couple_j_deposition_mode"): "cc_convert",
    ("particles", "couple_moments_momentum_to_mhd"): "true",
    ("particles", "couple_moments_energy_to_mhd"): "true",
    ("particles", "couple_fluid_feedback_order"): "mhd_src_terms",
    ("particles", "couple_moments_momentum_coeff"): "1.0",
    ("particles", "couple_moments_energy_coeff"): "1.0",
    ("particles", "pic_physical_mode"): "paper_mhd_pic",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_interp_scheme"): "tsc",
    ("particles", "pic_cr_light_speed"): "1000.0",
    ("particles", "pic_cr_initial_state"): "velocity",
    ("particles", "pic_cr_hall_mode"): "off",
    ("particles", "pic_wave_damping_mode"): "off",
    ("particles", "pic_max_cell_cross"): "1",
    ("particles", "pic_theta_max"): "0.1",
    ("particles", "pic_deltaf_mode"): "off",
    ("particles", "pic_expanding_box_mode"): "off",
    ("species0", "mass"): "1.0",
    ("species0", "charge"): "-1.0",
    ("species0", "vx0"): "0.0",
    ("species0", "vy0"): "0.1",
    ("species0", "vz0"): "0.0",
    ("species1", "mass"): "1.0",
    ("species1", "charge"): "1.0",
    ("species1", "vx0"): "0.0",
    ("species1", "vy0"): "0.1",
    ("species1", "vz0"): "0.0",
    ("problem", "pgen_name"): PGEN_NAME,
    ("q006_paper_multispecies_oscillation", "campaign_id"): CAMPAIGN_ID,
    ("q006_paper_multispecies_oscillation", "qualification_effect"): "none",
    ("q006_paper_multispecies_oscillation", "frontier_authorization"): "not_bound",
    ("q006_paper_multispecies_oscillation", "section53_runtime_evidence"):
        "not_claimed",
    ("q006_paper_multispecies_oscillation", "true_amr_policy_qualification"):
        "not_claimed",
    ("q006_paper_multispecies_oscillation", "mpi_qualification"): "not_claimed",
    ("q006_paper_multispecies_oscillation", "gpu_qualification"): "not_claimed",
    ("q006_paper_multispecies_oscillation", "external_review"): "not_claimed",
    ("q006_paper_multispecies_oscillation", "paper_eos_target"):
        "isothermal_cs_1",
    ("q006_paper_multispecies_oscillation", "runtime_eos_compatibility"):
        "ideal_energy_feedback_preparation_only",
    ("q006_paper_multispecies_oscillation", "ppc_semantics"):
        "aggregate_128_round_robin_yields_64_per_species_per_cell",
    ("q006_paper_multispecies_oscillation", "rho"): "1.0",
    ("q006_paper_multispecies_oscillation", "paper_cs"): "1.0",
    ("q006_paper_multispecies_oscillation", "compatibility_pressure"): "0.6",
    ("q006_paper_multispecies_oscillation", "b_g"): "1.0",
    ("q006_paper_multispecies_oscillation", "omega"): "1.0",
    ("q006_paper_multispecies_oscillation", "species_mass_density_ratio"): "1.5",
    ("q006_paper_multispecies_oscillation", "gas_uy"): "-0.3",
    ("q006_paper_multispecies_oscillation", "species_vy"): "0.1",
    ("q006_paper_multispecies_oscillation", "artificial_c"): "1000.0",
    ("q006_paper_multispecies_oscillation", "dt_omega_target"): "0.1",
    ("output1", "file_type"): "bin",
    ("output1", "variable"): "mhd_w_bcc",
    ("output1", "id"): "mhd_w_bcc",
    ("output2", "file_type"): "pvtk",
    ("output2", "variable"): "prtcl_all",
    ("output2", "id"): "prtcl_all",
    ("output3", "file_type"): "rst",
    ("output3", "id"): "rst",
}

_EXPECTED_MESH = {
    "nx": (16, 8, 8),
    "bounds": ((0.0, 16.0), (0.0, 8.0), (0.0, 8.0)),
    "meshblock_nx": (4, 4, 4),
}

_EXPECTED_GRID = {
    "uniform": {
        "deck_role": "source_local_uniform_deck_freeze_preparation_only_not_authorized",
        "amr_policy": "not_applicable_uniform",
        "refinement": "none",
        "num_levels": 1,
    },
    "smr": {
        "deck_role": "source_local_smr_deck_freeze_preparation_only_not_authorized",
        "amr_policy": "paper_static_one_eighth_region",
        "refinement": "static",
        "num_levels": 2,
    },
    "audited_amr_preparation": {
        "deck_role":
            "source_local_audited_amr_deck_freeze_preparation_only_not_authorized",
        "amr_policy":
            "deterministic_audited_randomized_preparation_10_refine_60_derefine_not_qualified",
        "refinement": "adaptive",
        "num_levels": 2,
    },
}


class ContractError(ValueError):
    """Raised when the bounded Q-006 source-local preparation drifts."""


def parse_athinput(path: Path) -> dict[str, dict[str, str]]:
    """Parse the strict Athena input subset used by the preparation decks."""
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


def _parse_float(label: str, value: str) -> float:
    try:
        measured = float(value)
    except ValueError as exc:
        raise ContractError(f"{label}: expected a finite number") from exc
    if not math.isfinite(measured):
        raise ContractError(f"{label}: expected a finite number")
    return measured


def _parse_int(label: str, value: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise ContractError(f"{label}: expected an integer") from exc


def _require_close(label: str, measured: float, expected: float) -> None:
    if not math.isclose(measured, expected, rel_tol=1.0e-14, abs_tol=1.0e-14):
        raise ContractError(f"{label}: expected {expected!r}, measured {measured!r}")


def aggregate_ppc_for_round_robin_species(species_ppc: float, nspecies: int) -> float:
    """Return the aggregate AthenaK ppc needed for round-robin species assignment."""
    if species_ppc <= 0.0 or nspecies <= 0:
        raise ContractError("species ppc and species count must be positive")
    return species_ppc * nspecies


def deposit_qscale_for_species_density(
    species_density: float, species_ppc: float, species_mass: float
) -> float:
    """Return macro-particle scale yielding one species mass density."""
    if species_density <= 0.0 or species_ppc <= 0.0 or species_mass <= 0.0:
        raise ContractError("density, ppc, and species mass must be positive")
    return species_density / (species_ppc * species_mass)


def gyro_frequency(q_over_mc: float, magnetic_field: float) -> float:
    """Return the non-relativistic gyro-frequency magnitude in code units."""
    return abs(q_over_mc * magnetic_field)


def oscillation_angular_frequency(omega: float, species_density_ratio: float) -> float:
    """Return Section 5.3 Omega*sqrt(1 + 2*m_e*n0/rho)."""
    if omega <= 0.0 or species_density_ratio <= 0.0:
        raise ContractError("Omega and species mass-density ratio must be positive")
    return omega * math.sqrt(1.0 + 2.0 * species_density_ratio)


def center_of_mass_gas_velocity(
    rho: float, species_density: float, species_velocity: float
) -> float:
    """Return gas velocity balancing two identical species in the COM frame."""
    if rho <= 0.0 or species_density <= 0.0:
        raise ContractError("gas and species densities must be positive")
    return -2.0 * species_density * species_velocity / rho


def ideal_pressure_for_sound_speed(rho: float, sound_speed: float, gamma: float) -> float:
    """Return ideal-MHD pressure with the requested initial adiabatic sound speed."""
    if rho <= 0.0 or sound_speed <= 0.0 or gamma <= 1.0:
        raise ContractError("ideal compatibility pressure inputs are invalid")
    return rho * sound_speed * sound_speed / gamma


def initial_kinetic_energy_density(
    rho: float, gas_velocity: float, species_density: float, species_velocity: float
) -> float:
    """Return Section 5.3 gas-plus-two-species kinetic energy density."""
    return (
        0.5 * rho * gas_velocity * gas_velocity
        + 0.5 * species_density * (species_velocity * species_velocity) * 2.0
    )


def _splitmix64(value: int) -> int:
    mask = (1 << 64) - 1
    value = (value + 0x9E3779B97F4A7C15) & mask
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & mask
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & mask
    return (value ^ (value >> 31)) & mask


def audited_amr_uniform01(seed: int, cycle: int, gid: int) -> float:
    """Mirror the deterministic source-local AMR audit draw."""
    mask = (1 << 64) - 1
    key = seed & mask
    key ^= ((cycle + 1) * 0xBF58476D1CE4E5B9) & mask
    key ^= ((gid + 1) * 0xD2B74407B1CE6E93) & mask
    return (_splitmix64(key) >> 11) / 9007199254740992.0


def audited_amr_decision(seed: int, cycle: int, gid: int) -> int:
    """Return +1 refine, -1 derefine, or 0 retain for the frozen audit draw."""
    draw = audited_amr_uniform01(seed, cycle, gid)
    if draw < AMR_REFINE_PROBABILITY:
        return 1
    if draw < AMR_REFINE_PROBABILITY + AMR_DEREFINE_PROBABILITY:
        return -1
    return 0


def audited_amr_policy_summary(cycles: int = 128, gids: int = 256) -> dict[str, Any]:
    """Summarize a bounded deterministic policy sample without claiming runtime evidence."""
    counts = {-1: 0, 0: 0, 1: 0}
    for cycle in range(cycles):
        for gid in range(gids):
            counts[audited_amr_decision(AMR_SEED, cycle, gid)] += 1
    total = cycles * gids
    return {
        "seed": AMR_SEED,
        "sample_count": total,
        "refine_fraction": counts[1] / total,
        "derefine_fraction": counts[-1] / total,
        "retain_fraction": counts[0] / total,
        "runtime_exercised": False,
        "true_amr_policy_qualification": False,
    }


def _validate_geometry(path: Path, blocks: dict[str, dict[str, str]]) -> None:
    nx = tuple(
        _parse_int(f"{path}: mesh/nx{axis}", blocks["mesh"][f"nx{axis}"])
        for axis in (1, 2, 3)
    )
    mb_nx = tuple(
        _parse_int(f"{path}: meshblock/nx{axis}", blocks["meshblock"][f"nx{axis}"])
        for axis in (1, 2, 3)
    )
    bounds = tuple(
        (
            _parse_float(f"{path}: mesh/x{axis}min", blocks["mesh"][f"x{axis}min"]),
            _parse_float(f"{path}: mesh/x{axis}max", blocks["mesh"][f"x{axis}max"]),
        )
        for axis in (1, 2, 3)
    )
    if nx != _EXPECTED_MESH["nx"] or mb_nx != _EXPECTED_MESH["meshblock_nx"]:
        raise ContractError(f"{path}: unexpected 3D carrier cell-count geometry")
    if bounds != _EXPECTED_MESH["bounds"]:
        raise ContractError(f"{path}: unexpected 3D carrier bounds")
    for axis in (1, 2, 3):
        mesh = blocks["mesh"]
        if mesh[f"ix{axis}_bc"] != "periodic" or mesh[f"ox{axis}_bc"] != "periodic":
            raise ContractError(f"{path}: all carrier boundaries must be periodic")


def validate_candidate_deck(path: Path, grid_setup: str) -> dict[str, Any]:
    """Validate one cycle-zero non-authorized source-local preparation deck."""
    blocks = parse_athinput(path)
    for (block, name), expected in _EXPECTED_COMMON.items():
        measured = blocks.get(block, {}).get(name)
        if measured != expected:
            raise ContractError(
                f"{path}: {block}/{name}: expected {expected!r}, measured {measured!r}"
            )
    _validate_geometry(path, blocks)

    expected = _EXPECTED_GRID[grid_setup]
    metadata = blocks["q006_paper_multispecies_oscillation"]
    for name in ("deck_role", "amr_policy"):
        if metadata.get(name) != expected[name]:
            raise ContractError(f"{path}: unexpected {name} for {grid_setup}")
    if metadata.get("grid_setup") != grid_setup:
        raise ContractError(f"{path}: unexpected grid_setup")

    refinement = blocks["mesh_refinement"]
    if refinement.get("refinement") != expected["refinement"]:
        raise ContractError(f"{path}: unexpected mesh-refinement mode")
    if _parse_int(f"{path}: mesh_refinement/num_levels", refinement["num_levels"]) != (
        expected["num_levels"]
    ):
        raise ContractError(f"{path}: unexpected refinement-level count")
    for name in ("ncycle_check", "refinement_interval"):
        if refinement.get(name) != "1":
            raise ContractError(f"{path}: mesh_refinement/{name} must remain 1")

    smr_fraction = None
    if grid_setup == "smr":
        region = blocks.get("refinement1")
        if region is None:
            raise ContractError(f"{path}: SMR preparation requires refinement1")
        expected_region = {
            "level": "1",
            "x1min": "4.0",
            "x1max": "12.0",
            "x2min": "2.0",
            "x2max": "6.0",
            "x3min": "2.0",
            "x3max": "6.0",
        }
        if region != expected_region:
            raise ContractError(f"{path}: SMR one-eighth region drifted")
        smr_fraction = (8.0 * 4.0 * 4.0) / (16.0 * 8.0 * 8.0)
        _require_close(f"{path}: SMR refined-volume fraction", smr_fraction, 0.125)
    elif "refinement1" in blocks:
        raise ContractError(f"{path}: only the SMR deck may include refinement1")

    if grid_setup == "audited_amr_preparation":
        if refinement.get("max_nmb_per_rank") != "512":
            raise ContractError(f"{path}: audited AMR max_nmb_per_rank drifted")
        if metadata.get("amr_seed") != str(AMR_SEED):
            raise ContractError(f"{path}: audited AMR seed drifted")

    rho = _parse_float(f"{path}: rho", metadata["rho"])
    species_density_ratio = _parse_float(
        f"{path}: species_mass_density_ratio", metadata["species_mass_density_ratio"]
    )
    species_ppc = 64.0
    qscale = _parse_float(
        f"{path}: particles/deposit_qscale", blocks["particles"]["deposit_qscale"]
    )
    _require_close(
        f"{path}: aggregate ppc",
        _parse_float(f"{path}: particles/ppc", blocks["particles"]["ppc"]),
        aggregate_ppc_for_round_robin_species(species_ppc, 2),
    )
    _require_close(
        f"{path}: species macro density",
        species_ppc * qscale * _parse_float(f"{path}: species0/mass", blocks["species0"]["mass"]),
        species_density_ratio * rho,
    )
    _require_close(
        f"{path}: ideal compatibility pressure",
        _parse_float(f"{path}: compatibility_pressure", metadata["compatibility_pressure"]),
        ideal_pressure_for_sound_speed(rho, 1.0, 5.0 / 3.0),
    )
    return {
        "grid_setup": grid_setup,
        "deck_path": path.relative_to(REPO_ROOT).as_posix(),
        "deck_role": expected["deck_role"],
        "amr_policy": expected["amr_policy"],
        "smr_refined_volume_fraction": smr_fraction,
        "cycle_zero_only": True,
        "qualifying_evidence": False,
    }


def validate_source_local_candidate_decks() -> list[dict[str, Any]]:
    """Validate all three bounded source-local preparation decks."""
    return [validate_candidate_deck(path, grid) for grid, path in DECKS.items()]


def validate_source_contract() -> dict[str, Any]:
    """Validate the dedicated source while reporting parent-owned registration state."""
    source = SOURCE.read_text(encoding="utf-8")
    required_snippets = (
        f"void ProblemGenerator::{PGEN_METHOD}",
        "Q006AuditedAMRRefinement",
        "const Real refine_probability = static_cast<Real>(0.10);",
        "const Real derefine_probability = static_cast<Real>(0.60);",
        "user_ref_func = Q006AuditedAMRRefinement;",
        'Q006RequireString(pin, "particles", "couple_j_deposition_mode", "cc_convert");',
        'Q006RequireReal(pin, "particles", "couple_moments_energy_coeff", 1.0);',
        '"ideal_energy_feedback_preparation_only"',
        '"deterministic_audited_randomized_preparation_10_refine_60_derefine_"',
        '"not_qualified"',
        "w0(m, IVY, k, j, i) = gas_uy;",
        "b0.x3f(m, k, j, i) = b_g;",
        "pmbp->pmhd->peos->PrimToCons",
    )
    missing = [snippet for snippet in required_snippets if snippet not in source]
    if missing:
        raise ContractError(f"Q-006 dedicated source contract is incomplete: {missing!r}")

    cmake = (REPO_ROOT / "src/CMakeLists.txt").read_text(encoding="utf-8")
    declarations = (REPO_ROOT / "src/pgen/pgen.hpp").read_text(encoding="utf-8")
    dispatch = (REPO_ROOT / "src/pgen/pgen.cpp").read_text(encoding="utf-8")
    registration = {
        "compilation_unit_registered":
            "pgen/tests/q006_paper_multispecies_oscillation.cpp" in cmake,
        "method_declared":
            f"void {PGEN_METHOD}(ParameterInput *pin, const bool restart);" in declarations,
        "fresh_dispatch_registered": f"{PGEN_METHOD}(pin, false);" in dispatch,
        "restart_dispatch_registered": f"{PGEN_METHOD}(pin, true);" in dispatch,
    }
    registration["complete"] = all(registration.values())
    return {
        "source": SOURCE.relative_to(REPO_ROOT).as_posix(),
        "pgen_name": PGEN_NAME,
        "method": PGEN_METHOD,
        "dedicated_source_contract": True,
        "audited_amr_callback_prepared": True,
        "shared_registration": registration,
    }


def analytical_contract() -> dict[str, float]:
    """Return the frozen Section 5.3 analytical preparation values."""
    rho = 1.0
    species_density = 1.5
    species_velocity = 0.1
    omega = gyro_frequency(1.0, 1.0)
    oscillation_omega = oscillation_angular_frequency(omega, species_density / rho)
    gas_velocity = center_of_mass_gas_velocity(rho, species_density, species_velocity)
    return {
        "rho": rho,
        "species_mass_density_each": species_density,
        "species_ppc_each": 64.0,
        "aggregate_ppc": aggregate_ppc_for_round_robin_species(64.0, 2),
        "deposit_qscale": deposit_qscale_for_species_density(species_density, 64.0, 1.0),
        "gyro_omega": omega,
        "oscillation_omega": oscillation_omega,
        "oscillation_frequency_cycles": oscillation_omega / (2.0 * math.pi),
        "oscillation_period": 2.0 * math.pi / oscillation_omega,
        "gas_uy": gas_velocity,
        "initial_kinetic_energy_density": initial_kinetic_energy_density(
            rho, gas_velocity, species_density, species_velocity
        ),
        "ideal_compatibility_pressure": ideal_pressure_for_sound_speed(
            rho, 1.0, 5.0 / 3.0
        ),
    }


def parent_registration_steps() -> list[str]:
    """Return any missing shared registration edits."""
    source = validate_source_contract()["shared_registration"]
    if source["complete"]:
        return []
    return [
        "Add pgen/tests/q006_paper_multispecies_oscillation.cpp to src/CMakeLists.txt.",
        "Declare ProblemGenerator::Q006PaperMultispeciesOscillation in src/pgen/pgen.hpp.",
        "Dispatch problem/pgen_name=q006_paper_multispecies_oscillation for fresh and "
        "restart construction in src/pgen/pgen.cpp.",
    ]


def build_preparation_report() -> dict[str, Any]:
    """Build the bounded report; no runtime artifacts are accepted or emitted."""
    decks = validate_source_local_candidate_decks()
    source = validate_source_contract()
    return {
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "status": "source_local_deck_and_source_freeze_preparation_only",
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "deck_source_freeze": True,
        "exact_section53_isothermal_runtime_supported": False,
        "runtime_eos_compatibility":
            "ideal_gamma_5_3_pressure_0p6_energy_feedback_preparation_only",
        "long_horizon_runtime_evidence": False,
        "true_amr_policy_qualification": False,
        "mpi_qualification": False,
        "gpu_qualification": False,
        "frontier_authorization": False,
        "external_review": False,
        "decks": decks,
        "source_contract": source,
        "analytical_contract": analytical_contract(),
        "audited_amr_policy_preparation": audited_amr_policy_summary(),
        "parent_registration_steps": parent_registration_steps(),
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indent", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    print(json.dumps(build_preparation_report(), indent=args.indent, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
