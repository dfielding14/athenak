#!/usr/bin/env python3
"""Nonqualifying Q-007 Sun-Bai true-delta-f linear preparation audit."""

from __future__ import annotations

import argparse
import cmath
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import stat
import sys
from typing import Any

import numpy as np
from scipy.integrate import quad

if __package__:
    from .immutable_orion_tree import authorized_tree_root
    from .immutable_orion_tree import freeze_tree as freeze_immutable_tree
    from .immutable_orion_tree import validate_executable_elf
    from .immutable_orion_tree import validate_serial_host_build_evidence
    from .immutable_orion_tree import validate_source_archive
    from .immutable_orion_tree import validate_source_archive_dependencies
    from .immutable_orion_tree import verify_frozen_tree as verify_immutable_tree
else:
    from immutable_orion_tree import authorized_tree_root
    from immutable_orion_tree import freeze_tree as freeze_immutable_tree
    from immutable_orion_tree import validate_executable_elf
    from immutable_orion_tree import validate_serial_host_build_evidence
    from immutable_orion_tree import validate_source_archive
    from immutable_orion_tree import validate_source_archive_dependencies
    from immutable_orion_tree import verify_frozen_tree as verify_immutable_tree


REPO_ROOT = Path(__file__).resolve().parents[2]
QUALIFICATION_EFFECT = "none_source_local_true_deltaf_preparation_only"
ARTIFACT_ROLE = "section55_section56_true_deltaf_preparation_not_qualifying_evidence"
ORION_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
MOMENTUM_BIN_COUNT = 8
PARTICLES_PER_CELL_PER_BIN = 256
PARTICLES_PER_CELL_TOTAL = MOMENTUM_BIN_COUNT * PARTICLES_PER_CELL_PER_BIN
SOURCE_LOCAL_CELL_COUNT = 32 * 4
SOURCE_LOCAL_PARTICLE_TOTAL = PARTICLES_PER_CELL_TOTAL * SOURCE_LOCAL_CELL_COUNT
WAVE_MODE_COUNT = 8
CRSI_REPLAY_CYCLE_LIMIT = 2
MASK64 = (1 << 64) - 1
INVENTORY_NAME = "artifact_inventory.sha256"
FREEZE_RECEIPT_NAME = "freeze_receipt.json"
PINNED_EXECUTABLE_BINDING_NAME = "pinned_executable_binding.json"
PIN_PROVENANCE_RECEIPT_NAME = "build/provenance_receipt.json"
PIN_RETAINED_EVIDENCE = {
    "build/CMakeCache.txt",
    "build/CMakeConfigureLog.yaml",
    "build/athena.flags.make",
    "build/athena.link.txt",
    "build/build.log",
    "build/build_command.txt",
    "build/build_profile.json",
    "build/clean_directory_preflight.json",
    "build/compile_commands.json",
    "build/compiler_identity.txt",
    "build/config.hpp",
    "build/configure.log",
    "build/configure_command.txt",
    "build/configure_compile_commands_refresh.log",
    "build/configure_compile_commands_refresh_command.txt",
    "build/verbose_clean_rebuild.log",
    "build/verbose_clean_rebuild_command.txt",
    "runtime/analysis_environment.json",
    "runtime/dynamic_dependencies.json",
    "runtime/ldd.txt",
    "runtime/module_environment_receipt.json",
    "runtime/modules.txt",
    "runtime/readelf_dynamic.txt",
    "runtime/readelf_header.txt",
    "runtime/readelf_program_headers.txt",
    "source/archive_validation.json",
    "source/exact_build_worktree_source.tar.gz",
    "source/head_commit.txt",
    "source/shared_dependency_sha256.json",
    "source/tracked_worktree.diff",
    "source/worktree_status.txt",
}
PIN_BUILD_DEPENDENCY_PATHS = {
    "src/CMakeLists.txt",
    "src/particles/particles.cpp",
    "src/pgen/pgen.cpp",
    "src/pgen/pgen.hpp",
    "src/pgen/tests/q006_paper_multispecies_oscillation_runtime_local.cpp",
    "src/pgen/tests/q007_paper_deltaf_linear.cpp",
    "src/pgen/tests/q007_paper_deltaf_linear.hpp",
    "tst/publication/analyze_q006_paper_multispecies_oscillation_runtime_local.py",
    "tst/publication/analyze_q007_paper_deltaf_linear_preparation.py",
    "tst/publication/analyze_q011_injection_distribution_runtime_local.py",
    "tst/publication/immutable_orion_tree.py",
    "tst/scripts/particles/pic_mhd_expanding_box_cpaw_history_preparation.py",
    "tst/scripts/particles/pic_parser_contract_guards.py",
}
PINNED_EXECUTABLE_BINDING_KEYS = {
    "artifact_role",
    "frontier_used",
    "kronos_used",
    "mpi_used",
    "negative_crpai_nlim1_disposition",
    "pinned_executable_path",
    "pinned_executable_root",
    "pinned_executable_root_inventory_sha256",
    "pinned_executable_sha256",
    "qualification_effect",
    "schema_version",
    "slurm_used",
}
RUNTIME_INVOCATION_KEYS = {
    "argv",
    "cwd",
    "deck",
    "executable_realpath",
    "executable_sha256",
    "frontier_used",
    "kronos_used",
    "launch_style",
    "mpi_used",
    "overrides",
    "pinned_executable_root",
    "pinned_executable_root_inventory_sha256",
    "produced_file_sha256",
    "returncode",
    "schema_version",
    "selected_environment",
    "slurm_used",
}
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
    ("mesh", "nghost"): "2",
    ("mesh", "nx1"): "32",
    ("mesh", "x1min"): "0.0",
    ("mesh", "ix1_bc"): "periodic",
    ("mesh", "ox1_bc"): "periodic",
    ("mesh", "nx2"): "4",
    ("mesh", "x2min"): "0.0",
    ("mesh", "x2max"): "4.0",
    ("mesh", "ix2_bc"): "periodic",
    ("mesh", "ox2_bc"): "periodic",
    ("mesh", "nx3"): "1",
    ("mesh", "x3min"): "0.0",
    ("mesh", "x3max"): "1.0",
    ("mesh", "ix3_bc"): "periodic",
    ("mesh", "ox3_bc"): "periodic",
    ("meshblock", "nx1"): "32",
    ("meshblock", "nx2"): "4",
    ("meshblock", "nx3"): "1",
    ("mesh_refinement", "refinement"): "none",
    ("mesh_refinement", "num_levels"): "1",
    ("time", "evolution"): "dynamic",
    ("time", "integrator"): "rk2",
    ("time", "cfl_number"): "0.1",
    ("time", "ndiag"): "1",
    ("mhd", "eos"): "isothermal",
    ("mhd", "iso_sound_speed"): "1.0",
    ("mhd", "reconstruct"): "plm",
    ("mhd", "rsolver"): "llf",
    ("coord", "special_rel"): "false",
    ("coord", "general_rel"): "false",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "ppc"): "2048.0",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "nspecies"): "8",
    ("particles", "cr_distribution"): "center",
    ("particles", "deposit_moments"): "true",
    ("particles", "deposit_order"): "1",
    ("particles", "deposit_qscale"): "1.0e-4",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "couple_j_to_efield_coeff"): "1.0",
    ("particles", "couple_j_to_efield_representation"): "cell_centered",
    ("particles", "couple_j_deposition_mode"): "cc_convert",
    ("particles", "couple_moments_momentum_to_mhd"): "true",
    ("particles", "couple_moments_energy_to_mhd"): "false",
    ("particles", "couple_moments_momentum_coeff"): "1.0",
    ("particles", "couple_moments_energy_coeff"): "1.0",
    ("particles", "couple_fluid_feedback_order"): "mhd_src_terms",
    ("particles", "cr_vx0"): "0.0",
    ("particles", "cr_vy0"): "0.0",
    ("particles", "cr_vz0"): "0.0",
    ("particles", "pic_physical_mode"): "paper_mhd_pic",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_interp_scheme"): "tsc",
    ("particles", "pic_enable_2d3v"): "true",
    ("particles", "pic_cr_initial_state"): "momentum",
    ("particles", "pic_cr_hall_mode"): "off",
    ("particles", "pic_wave_damping_mode"): "off",
    ("particles", "pic_max_cell_cross"): "1",
    ("particles", "pic_theta_max"): "0.3",
    ("particles", "pic_deltaf_mode"): "physical",
    ("particles", "pic_deltaf_drift_x1"): "0.0",
    ("particles", "pic_deltaf_drift_x2"): "0.0",
    ("particles", "pic_deltaf_drift_x3"): "0.0",
    ("particles", "pic_deltaf_background_rho"): "1.0e-4",
    ("particles", "pic_deltaf_background_jx"): "0.0",
    ("particles", "pic_deltaf_background_jy"): "0.0",
    ("particles", "pic_deltaf_background_jz"): "0.0",
    ("particles", "pic_deltaf_adapt_mode"): "off",
    ("particles", "pic_deltaf_adapt_interval"): "0.0",
    ("particles", "pic_sort_interval"): "0",
    ("particles", "pic_intermediate_arrays"): "auto",
    ("particles", "pic_expanding_box_mode"): "off",
    ("output1", "file_type"): "bin",
    ("output1", "variable"): "mhd_w_bcc",
    ("output1", "id"): "mhd_w_bcc",
    ("output1", "dcycle"): "1",
    ("output1", "ghost_zones"): "false",
    ("output2", "file_type"): "pvtk",
    ("output2", "variable"): "prtcl_all",
    ("output2", "id"): "prtcl_all",
    ("output2", "dcycle"): "1",
}
_METADATA_COMMON = {
    "qualification_effect": "none",
    "frontier_authorization": "not_bound",
    "physical_loading": "implemented_source_local_eight_log_bin_ipwt_quadrature",
    "initial_wave_spectrum":
        "implemented_source_local_deterministic_four_branch_carrier",
    "theory_runtime_comparison": "q1_q2_oracle_preparation_only_no_growth_fit",
    "paper_momentum_bin_count": "8",
    "paper_particles_per_cell_per_bin": "256",
    "source_local_particles_per_cell_total": "2048",
    "source_local_particles_total": "262144",
    "weight_encoding": "ipwt_macro_multiplicity_equivalent_equal_q_over_mc",
    "loading_quadrature":
        "geometric_center_shell_midpoint_normalized_p0_over_500_to_500_p0",
    "wave_discrete_normalization": "branch_amplitude_A_over_sqrt_mode",
    "wave_mode_count": "8",
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
        "comment": "Q007 paper CRSI linear true-delta-f source-local preparation",
        "basename": "pic_q007_paper_crsi_linear_preparation",
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
        "nlim": "2",
        "tlim": "1.0",
        "deck_role": "bounded_two_cycle_serial_runtime_replay_preparation_only",
        "runtime_evolution": "admitted_two_cycle_serial_replay_only",
        "momentum_seed": "700701",
        "wave_seed": "700702",
        "runtime_replay_cycle_limit": "2",
        "paper_seed_amplitude": "1.0e-3",
    },
    "crpai_prolate": {
        "block": "q007_paper_crpai_linear_preparation",
        "pgen_name": "q007_paper_crpai_linear_preparation",
        "comment": "Q007 paper CRPAI prolate true-delta-f source-local preparation",
        "basename": "pic_q007_paper_crpai_linear_prolate_preparation",
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
        "nlim": "0",
        "tlim": "0.0",
        "deck_role": "cycle_zero_weighted_loading_wave_oracle_preparation_only",
        "runtime_evolution": "blocked_cycle_zero_only",
        "momentum_seed": "700703",
        "wave_seed": "700704",
        "runtime_replay_cycle_limit": "0",
        "paper_seed_amplitude": "1.0e-3",
    },
    "crpai_oblate": {
        "block": "q007_paper_crpai_linear_preparation",
        "pgen_name": "q007_paper_crpai_linear_preparation",
        "comment": "Q007 paper CRPAI oblate true-delta-f source-local preparation",
        "basename": "pic_q007_paper_crpai_linear_oblate_preparation",
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
        "nlim": "0",
        "tlim": "0.0",
        "deck_role": "cycle_zero_weighted_loading_wave_oracle_preparation_only",
        "runtime_evolution": "blocked_cycle_zero_only",
        "momentum_seed": "700703",
        "wave_seed": "700704",
        "runtime_replay_cycle_limit": "0",
        "paper_seed_amplitude": "1.0e-3",
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


def _allowed_deck_parameters(case: str) -> dict[str, set[str]]:
    """Return the exact allowed block and parameter names for one Q-007 deck."""
    expected_case = _CASES[case]
    allowed: dict[str, set[str]] = {}

    def add(block: str, *names: str) -> None:
        allowed.setdefault(block, set()).update(names)

    for block, name in _COMMON:
        add(block, name)
    add("comment", "problem")
    add("job", "basename")
    add("mesh", "x1max")
    add("time", "nlim", "tlim")
    add("problem", "pgen_name")
    add(
        "particles",
        "pic_cr_light_speed",
        "pic_deltaf_f0",
        "pic_deltaf_p0",
        "pic_deltaf_kappa",
        "pic_deltaf_aniso_x1",
        "pic_deltaf_aniso_x2",
        "pic_deltaf_aniso_x3",
    )
    for species in range(MOMENTUM_BIN_COUNT):
        add(f"species{species}", "mass", "charge", "vx0", "vy0", "vz0")
    add(expected_case["block"], *_METADATA_COMMON)
    add(
        expected_case["block"],
        "campaign_id",
        "section_anchor",
        "paper_dx",
        "source_local_x1",
        "paper_light_speed",
        "paper_p0",
        "paper_kappa",
        "paper_xi",
        "paper_seed_amplitude",
        "case_role",
        "source_local_gas_vx",
        "handedness_mapping",
        "deck_role",
        "runtime_evolution",
        "momentum_seed",
        "wave_seed",
        "runtime_replay_cycle_limit",
    )
    if case == "crsi":
        add(expected_case["block"], "paper_vd")
    return allowed


def _require_exact_deck_map(
    path: Path, blocks: dict[str, dict[str, str]], case: str
) -> None:
    """Reject every block or parameter outside the bounded Q-007 deck surface."""
    allowed = _allowed_deck_parameters(case)
    missing_blocks = sorted(set(allowed) - set(blocks))
    unexpected_blocks = sorted(set(blocks) - set(allowed))
    if missing_blocks or unexpected_blocks:
        raise ContractError(
            f"{path}: exact allowed block map drifted: "
            f"missing={missing_blocks}, unexpected={unexpected_blocks}"
        )
    for block, expected_names in allowed.items():
        measured_names = set(blocks[block])
        missing = sorted(expected_names - measured_names)
        unexpected = sorted(measured_names - expected_names)
        if missing or unexpected:
            raise ContractError(
                f"{path}: <{block}> exact allowed parameter map drifted: "
                f"missing={missing}, unexpected={unexpected}"
            )


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


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def splitmix64(value: int) -> int:
    """Match the Q-007 source-local uint64 SplitMix64 implementation."""
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def deterministic_uniform01(seed: int, item: int, component: int) -> float:
    """Return one deterministic source-local replay draw."""
    key = seed & MASK64
    key ^= ((item + 1) * 0xBF58476D1CE4E5B9) & MASK64
    key ^= ((component + 1) * 0xD2B74407B1CE6E93) & MASK64
    return (splitmix64(key) >> 11) / 9007199254740992.0


def log_bin_edge(p0: float, edge: int) -> float:
    """Return one of the nine paper-range logarithmic momentum-bin edges."""
    if edge < 0 or edge > MOMENTUM_BIN_COUNT:
        raise ContractError("log-bin edge index is outside the bounded contract")
    return _positive("p0", p0) / 500.0 * 250000.0 ** (
        edge / MOMENTUM_BIN_COUNT
    )


def log_bin_center(p0: float, bin_index: int) -> float:
    """Return the geometric center of one bounded momentum bin."""
    if bin_index < 0 or bin_index >= MOMENTUM_BIN_COUNT:
        raise ContractError("log-bin index is outside the bounded contract")
    return math.sqrt(log_bin_edge(p0, bin_index) * log_bin_edge(p0, bin_index + 1))


def log_bin_raw_shell_weight(p0: float, kappa: float, bin_index: int) -> float:
    """Return the source-local shell-midpoint quadrature weight before normalization."""
    lower = log_bin_edge(p0, bin_index)
    upper = log_bin_edge(p0, bin_index + 1)
    center = log_bin_center(p0, bin_index)
    shape = (1.0 + center * center / (kappa * p0 * p0)) ** (-kappa - 1.0)
    return (upper**3 - lower**3) * shape


def log_bin_fraction(p0: float, kappa: float, bin_index: int) -> float:
    """Return one normalized source-local IPWT shell fraction."""
    weights = [
        log_bin_raw_shell_weight(p0, kappa, candidate)
        for candidate in range(MOMENTUM_BIN_COUNT)
    ]
    return weights[bin_index] / sum(weights)


def wave_phase(seed: int, mode: int, direction: int, polarization: int) -> float:
    """Return one deterministic local branch phase."""
    if mode < 1 or mode > WAVE_MODE_COUNT:
        raise ContractError("wave mode is outside the bounded replay contract")
    if direction not in (-1, 1) or polarization not in (-1, 1):
        raise ContractError("wave direction and polarization must be signed branches")
    branch = 2 * int(direction > 0) + int(polarization > 0)
    return 2.0 * math.pi * deterministic_uniform01(seed, mode, branch)


def wave_branch_amplitude(amplitude: float, mode: int) -> float:
    """Freeze the local discrete interpretation of I(k)=A^2/|k|."""
    return _positive("wave amplitude", amplitude) / math.sqrt(mode)


def four_branch_alfven_state(
    x1: float, length: float, amplitude: float, seed: int
) -> tuple[float, float, float, float]:
    """Return the deterministic local By, Bz, Uy, Uz four-branch carrier."""
    by = bz = uy = uz = 0.0
    for mode in range(1, WAVE_MODE_COUNT + 1):
        branch_amplitude = wave_branch_amplitude(amplitude, mode)
        wave_number = 2.0 * math.pi * mode / _positive("length", length)
        for direction in (-1, 1):
            for polarization in (-1, 1):
                phase = wave_number * x1 + wave_phase(
                    seed, mode, direction, polarization
                )
                branch_by = branch_amplitude * math.cos(phase)
                branch_bz = (
                    polarization * branch_amplitude * math.sin(phase)
                )
                by += branch_by
                bz += branch_bz
                uy -= direction * branch_by
                uz -= direction * branch_bz
    return by, bz, uy, uz


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


@lru_cache(maxsize=None)
def gyroresonance_q1(
    k: float, mass: float, omega0: float, p0: float, kappa: float
) -> float:
    """Numerically evaluate the paper's principal-value-free Q1 log integral."""
    wave_number = abs(_finite("k", k))
    if wave_number == 0.0:
        raise ContractError("k must be nonzero")
    index = _positive("kappa", kappa)
    ratio = _positive("mass", mass) * _positive("omega0", omega0)
    ratio /= wave_number * _positive("p0", p0)

    def integrand(s: float) -> float:
        if s == 1.0:
            return 0.0
        logarithm = math.log(abs((1.0 + s) / (1.0 - s)))
        return logarithm * (1.0 + s * s * ratio * ratio / index) ** (
            -index - 1.0
        ) * s

    integral = quad(integrand, 0.0, 1.0, epsabs=1.0e-12, epsrel=1.0e-12)[0]
    integral += quad(integrand, 1.0, math.inf, epsabs=1.0e-12, epsrel=1.0e-12)[0]
    prefactor = (
        2.0
        / (math.sqrt(math.pi) * index**1.5)
        * math.gamma(index + 1.0)
        / math.gamma(index - 0.5)
    )
    return prefactor * ratio**3 * integral


def _ordered_dispersion_roots(
    linear: complex, constant: complex, k: float, ua: float
) -> dict[str, complex]:
    discriminant = linear * linear - 4.0 * constant
    roots = (
        (-linear + cmath.sqrt(discriminant)) / 2.0,
        (-linear - cmath.sqrt(discriminant)) / 2.0,
    )
    forward_target = abs(k) * ua
    forward = min(roots, key=lambda root: abs(root - forward_target))
    backward = roots[0] if roots[1] == forward else roots[1]
    return {"forward": forward, "backward": backward}


def crsi_dispersion_roots(
    k: float,
    mass: float,
    omega0: float,
    p0: float,
    kappa: float,
    mncr_over_rho0: float,
    ua: float,
    vd: float,
    signed_polarization: int,
) -> dict[str, complex]:
    """Return full CRSI Q1+Q2 roots without assigning handedness labels."""
    if signed_polarization not in (-1, 1):
        raise ContractError("signed polarization must be -1 or +1")
    wave_number = _positive("|k|", abs(_finite("k", k)))
    alfven_speed = _positive("ua", ua)
    alpha_omega = _positive("mncr/rho0", mncr_over_rho0) * _positive(
        "omega0", omega0
    )
    q1 = gyroresonance_q1(wave_number, mass, omega0, p0, kappa)
    q2 = gyroresonance_q2(wave_number, mass, omega0, p0, kappa)
    response = 1.0 - q1 + signed_polarization * 1.0j * q2
    linear = signed_polarization * alpha_omega * response
    constant = (
        -wave_number * wave_number * alfven_speed * alfven_speed
        - signed_polarization * alpha_omega * response * wave_number * vd
    )
    return _ordered_dispersion_roots(linear, constant, wave_number, alfven_speed)


def crpai_dispersion_roots(
    k: float,
    mass: float,
    omega0: float,
    p0: float,
    kappa: float,
    mncr_over_rho0: float,
    ua: float,
    xi: float,
    signed_polarization: int,
) -> dict[str, complex]:
    """Return full CRPAI Q1+Q2 roots without assigning handedness labels."""
    if signed_polarization not in (-1, 1):
        raise ContractError("signed polarization must be -1 or +1")
    wave_number = _positive("|k|", abs(_finite("k", k)))
    alfven_speed = _positive("ua", ua)
    anisotropy = _positive("xi", xi)
    xi2 = anisotropy * anisotropy
    alpha_omega = _positive("mncr/rho0", mncr_over_rho0) * _positive(
        "omega0", omega0
    )
    omega = _positive("omega0", omega0)
    anisotropy_factor = (1.0 - xi2) / xi2
    q1 = gyroresonance_q1(wave_number, mass, omega0, p0, kappa)
    q2 = gyroresonance_q2(wave_number, mass, omega0, p0, kappa)
    response = 1.0 - q1 / xi2 + signed_polarization * 1.0j * q2 / xi2
    linear = signed_polarization * alpha_omega * response
    constant = (
        -wave_number * wave_number * alfven_speed * alfven_speed
        + alpha_omega * omega * anisotropy_factor * response
    )
    return _ordered_dispersion_roots(linear, constant, wave_number, alfven_speed)


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
    _require_exact_deck_map(path, blocks, case)
    for (block, name), expected in _COMMON.items():
        measured = blocks.get(block, {}).get(name)
        if measured != expected:
            raise ContractError(f"{path}: {block}/{name}: expected {expected!r}")
    if blocks["comment"].get("problem") != expected_case["comment"]:
        raise ContractError(f"{path}: comment/problem does not match bounded carrier")
    if blocks["job"].get("basename") != expected_case["basename"]:
        raise ContractError(f"{path}: job/basename does not match bounded carrier")
    if blocks["mesh"].get("x1max") != expected_case["x1max"]:
        raise ContractError(f"{path}: mesh/x1max does not match bounded carrier")
    for name in ("nlim", "tlim"):
        if blocks["time"].get(name) != expected_case[name]:
            raise ContractError(f"{path}: time/{name} does not match bounded carrier")
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
            raise ContractError(f"{path}: species{species} log-bin contract drifted")
    metadata = blocks[expected_case["block"]]
    for name, expected in _METADATA_COMMON.items():
        if metadata.get(name) != expected:
            raise ContractError(f"{path}: {expected_case['block']}/{name} drifted")
    for name, expected_name in (
        ("campaign_id", "campaign_id"),
        ("section_anchor", "section_anchor"),
        ("paper_dx", "paper_dx"),
        ("source_local_x1", "x1max"),
        ("paper_light_speed", "light_speed"),
        ("paper_p0", "p0"),
        ("paper_kappa", "kappa"),
        ("paper_xi", "xi"),
        ("paper_seed_amplitude", "paper_seed_amplitude"),
        ("case_role", "case_role"),
        ("source_local_gas_vx", "source_local_gas_vx"),
        ("handedness_mapping", "handedness_mapping"),
        ("deck_role", "deck_role"),
        ("runtime_evolution", "runtime_evolution"),
        ("momentum_seed", "momentum_seed"),
        ("wave_seed", "wave_seed"),
        ("runtime_replay_cycle_limit", "runtime_replay_cycle_limit"),
    ):
        if metadata.get(name) != expected_case[expected_name]:
            raise ContractError(f"{path}: {expected_case['block']}/{name} drifted")
    if case == "crsi" and metadata.get("paper_vd") != expected_case["paper_vd"]:
        raise ContractError(f"{path}: CRSI paper drift mapping drifted")
    particle_count = int(float(particles["ppc"]) * SOURCE_LOCAL_CELL_COUNT)
    if particle_count != SOURCE_LOCAL_PARTICLE_TOTAL:
        raise ContractError(f"{path}: weighted-loading particle count drifted")
    try:
        rendered_path = str(path.relative_to(REPO_ROOT))
    except ValueError:
        rendered_path = str(path)
    return {
        "case": case,
        "path": rendered_path,
        "pgen_name": expected_case["pgen_name"],
        "cycle_zero_only": expected_case["nlim"] == "0",
        "true_deltaf_parser_path": True,
        "exact_isothermal_mhd": True,
        "source_local_particles_total": particle_count,
        "paper_particles_per_cell_total": PARTICLES_PER_CELL_TOTAL,
        "physical_loading_implemented": True,
        "random_phase_wave_spectrum_implemented": True,
        "runtime_evolution_admitted": expected_case["nlim"] != "0",
        "runtime_cycle_limit": int(expected_case["nlim"]),
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
        "implemented_source_local_eight_log_bin_ipwt_quadrature",
        "implemented_source_local_deterministic_four_branch_carrier",
        "blocked_pending_manuscript_text_caption_review",
        "Q007RejectEffectfulOptionalMHDControls",
        "couple_moments_momentum_coeff",
        "couple_moments_energy_coeff",
        '"viscosity"',
        '"ohmic_resistivity"',
        '"shearing_box"',
        "pgen_q007_weighted_log_bin_antipodal_momenta",
        "pgen_q007_four_branch_isothermal_paper_state",
    ):
        if snippet not in source:
            raise ContractError(f"Q-007 source is missing fail-closed marker {snippet!r}")
    for snippet in (
        "IsotropicKappaDistribution",
        "AnisotropicKappaDistribution",
        "GyroresonanceQ2",
        "CRSILowDensityGrowthRate",
        "CRPAILowDensityGrowthRate",
        "LogBinFraction",
        "FourBranchAlfvenState",
    ):
        if snippet not in header:
            raise ContractError(f"Q-007 header is missing analytical helper {snippet!r}")
    return {
        "compilation_unit_registered": True,
        "fresh_and_restart_dispatch_registered": True,
        "narrow_exact_isothermal_true_deltaf_parser_allowance": True,
        "bounded_serial_replay_guarded": True,
        "effectful_optional_mhd_controls_rejected": True,
    }


def _complex_payload(value: complex) -> dict[str, float]:
    return {"real": float(value.real), "imag": float(value.imag)}


def loading_contract(p0: float, kappa: float) -> dict[str, Any]:
    """Describe the bounded source-local shell-midpoint weighted loading."""
    bins = []
    for bin_index in range(MOMENTUM_BIN_COUNT):
        fraction = log_bin_fraction(p0, kappa, bin_index)
        bins.append({
            "bin": bin_index,
            "lower": log_bin_edge(p0, bin_index),
            "center": log_bin_center(p0, bin_index),
            "upper": log_bin_edge(p0, bin_index + 1),
            "shell_fraction": fraction,
            "particle_macro_weight": fraction / PARTICLES_PER_CELL_PER_BIN,
        })
    return {
        "paper_range": "p0_over_500_to_500_p0",
        "bin_count": MOMENTUM_BIN_COUNT,
        "particles_per_cell_per_bin": PARTICLES_PER_CELL_PER_BIN,
        "particles_per_cell_total": PARTICLES_PER_CELL_TOTAL,
        "source_local_particles_total": SOURCE_LOCAL_PARTICLE_TOTAL,
        "weight_encoding": "ipwt_macro_multiplicity_equivalent_equal_q_over_mc",
        "quadrature":
            "geometric_center_shell_midpoint_normalized_over_truncated_range",
        "shell_fraction_sum": sum(item["shell_fraction"] for item in bins),
        "bins": bins,
        "paper_literal_sampling_algorithm_claimed": False,
    }


def wave_contract(length: float, amplitude: float, seed: int) -> dict[str, Any]:
    """Describe the deterministic local discrete four-branch wave carrier."""
    modes = []
    for mode in range(1, WAVE_MODE_COUNT + 1):
        branches = []
        for direction in (-1, 1):
            for polarization in (-1, 1):
                branches.append({
                    "direction": direction,
                    "signed_polarization": polarization,
                    "phase": wave_phase(seed, mode, direction, polarization),
                })
        modes.append({
            "mode": mode,
            "k": 2.0 * math.pi * mode / length,
            "branch_amplitude": wave_branch_amplitude(amplitude, mode),
            "branches": branches,
        })
    return {
        "seed": seed,
        "mode_count": WAVE_MODE_COUNT,
        "discrete_normalization": "branch_amplitude_A_over_sqrt_mode",
        "paper_intensity_mapping": "branch_amplitude_squared_over_delta_k_is_A2_over_abs_k",
        "modes": modes,
        "paper_literal_seed_or_discrete_mode_set_claimed": False,
    }


def _crsi_dispersion_oracle_table(length: float) -> dict[str, Any]:
    """Build CRSI Q1+Q2 roots on one bounded carrier's discrete modes."""
    modes = []
    for mode in range(1, WAVE_MODE_COUNT + 1):
        k = 2.0 * math.pi * mode / length
        crsi = {}
        for polarization in (-1, 1):
            crsi[str(polarization)] = {
                direction: _complex_payload(root)
                for direction, root in crsi_dispersion_roots(
                    k, 1.0, 1.0, 300.0, 1.25, 1.0e-4, 1.0, 2.0,
                    polarization
                ).items()
            }
        modes.append({
            "mode": mode,
            "k": k,
            "q1_crsi": gyroresonance_q1(k, 1.0, 1.0, 300.0, 1.25),
            "q2_crsi": gyroresonance_q2(k, 1.0, 1.0, 300.0, 1.25),
            "crsi_signed_polarizations": crsi,
        })
    return {
        "mode_count": WAVE_MODE_COUNT,
        "carrier_length": length,
        "modes": modes,
    }


def _crpai_dispersion_oracle_table(length: float) -> dict[str, Any]:
    """Build CRPAI Q1+Q2 roots on one bounded carrier's discrete modes."""
    modes = []
    for mode in range(1, WAVE_MODE_COUNT + 1):
        k = 2.0 * math.pi * mode / length
        crpai = {}
        for role, xi in (("prolate", 0.99), ("oblate", 1.01)):
            crpai[role] = {}
            for polarization in (-1, 1):
                crpai[role][str(polarization)] = {
                    direction: _complex_payload(root)
                    for direction, root in crpai_dispersion_roots(
                        k, 1.0, 1.0, 300.0, 1.75, 1.0e-4, 1.0, xi,
                        polarization
                    ).items()
                }
        modes.append({
            "mode": mode,
            "k": k,
            "q1_crpai": gyroresonance_q1(k, 1.0, 1.0, 300.0, 1.75),
            "q2_crpai": gyroresonance_q2(k, 1.0, 1.0, 300.0, 1.75),
            "crpai_signed_polarizations": crpai,
        })
    return {
        "mode_count": WAVE_MODE_COUNT,
        "carrier_length": length,
        "modes": modes,
    }


def full_q1_q2_dispersion_oracle() -> dict[str, Any]:
    """Build bounded carrier-specific Q1+Q2 roots for CRSI and CRPAI."""
    return {
        "mode_count": WAVE_MODE_COUNT,
        "root_semantics": "full_complex_quadratic_roots_near_forward_and_backward_alfven_modes",
        "assigns_crpai_handedness_labels": False,
        "runtime_growth_fit_claimed": False,
        "tables": {
            "crsi": _crsi_dispersion_oracle_table(320.0),
            "crpai": _crpai_dispersion_oracle_table(640.0),
        },
    }


def analytical_contract() -> dict[str, Any]:
    """Build static paper mappings and a bounded nonqualifying Q1+Q2 oracle."""
    k_crsi = resonant_wavenumber(1.0, 1.0, 300.0)
    q1_crsi = gyroresonance_q1(k_crsi, 1.0, 1.0, 300.0, 1.25)
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
        "full_q1_q2_dispersion_oracle_prepared": True,
        "runtime_growth_fit_claimed": False,
        "crsi": {
            "k0": k_crsi,
            "lambda0": 2.0 * math.pi / k_crsi,
            "q1_at_k0": q1_crsi,
            "q2_at_k0": q2_crsi,
            "forward_growth_at_k0": crsi_forward,
            "backward_growth_at_minus_k0": crsi_backward,
        },
        "crpai": crpai,
        "crpai_handedness_claimed": False,
        "crpai_handedness_boundary":
            "manuscript prose and figure caption require independent review",
        "loading_contracts": {
            "crsi": loading_contract(300.0, 1.25),
            "crpai": loading_contract(300.0, 1.75),
        },
        "wave_contracts": {
            "crsi": wave_contract(320.0, 1.0e-3, 700702),
            "crpai": wave_contract(640.0, 1.0e-3, 700704),
        },
        "full_q1_q2_dispersion_oracle": full_q1_q2_dispersion_oracle(),
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
            "exact_isothermal_mhd_startup": True,
            "source_local_log_bin_weighted_loading": True,
            "source_local_deterministic_four_branch_wave_spectrum": True,
            "bounded_serial_crsi_runtime_replay": True,
            "full_q1_q2_dispersion_oracle_prepared": True,
            "paper_literal_sampling_algorithm_reviewed": False,
            "paper_literal_wave_seed_and_discrete_mode_set_reviewed": False,
            "runtime_growth_fit": False,
            "crpai_handedness_mapping_reviewed": False,
            "mpi_qualification": False,
            "gpu_qualification": False,
            "frontier_authorization": False,
            "external_review": False,
            "qualifying_evidence": False,
        },
    }


def _authorized_orion_runtime_root(path: Path) -> Path:
    """Require one canonical real runtime root strictly below Orion PIC."""
    return authorized_tree_root(
        path,
        authorized_root=ORION_PIC_ROOT,
        error_type=ContractError,
        label="Q-007 runtime artifact root",
    )


def _contained_regular_file(root: Path, path: Path) -> Path:
    """Require one self-contained regular artifact below its retained root."""
    status = os.lstat(path)
    if not stat.S_ISREG(status.st_mode):
        raise ContractError(f"runtime artifact must be regular: {path}")
    if status.st_nlink != 1:
        raise ContractError(f"runtime artifact must have one link: {path}")
    resolved = path.resolve(strict=True)
    if resolved != path:
        raise ContractError(f"runtime artifact path must be canonical: {path}")
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise ContractError(f"runtime artifact must remain below retained root: {path}") from error
    return resolved


def freeze_runtime_tree(runtime_root: Path) -> dict[str, Any]:
    """Inventory a retained Orion runtime tree and recursively remove write bits."""
    receipt = {
        "schema_version": 1,
        "gate": "Q-007",
        "artifact_role": "bounded_serial_crsi_runtime_replay_mechanics_only",
        "qualification_effect": QUALIFICATION_EFFECT,
        "inventory_excludes": INVENTORY_NAME,
        "freeze_policy": "remove all owner, group and other write bits recursively",
    }
    return freeze_immutable_tree(
        runtime_root,
        receipt,
        authorized_root=ORION_PIC_ROOT,
        error_type=ContractError,
        label="Q-007 retained runtime tree",
    )


def _validate_runtime_freeze_receipt(root: Path) -> str:
    """Require the retained Q-007 tree to carry its exact nonqualifying role."""
    path = _contained_regular_file(root, root / FREEZE_RECEIPT_NAME)
    receipt = json.loads(path.read_text(encoding="utf-8"))
    expected = {
        "schema_version": 1,
        "gate": "Q-007",
        "artifact_role": "bounded_serial_crsi_runtime_replay_mechanics_only",
        "qualification_effect": QUALIFICATION_EFFECT,
        "inventory_excludes": INVENTORY_NAME,
        "freeze_policy": "remove all owner, group and other write bits recursively",
    }
    if receipt != expected:
        raise ContractError("Q-007 freeze receipt drifted")
    return _sha256(path)


def verify_frozen_runtime_tree(
    runtime_root: Path, expected_inventory_sha256: str
) -> dict[str, Any]:
    """Verify exact payload membership, hashes, and recursive read-only modes."""
    root = _authorized_orion_runtime_root(runtime_root)
    report = verify_immutable_tree(
        root,
        expected_inventory_sha256,
        authorized_root=ORION_PIC_ROOT,
        error_type=ContractError,
        label="Q-007 retained runtime tree",
    )
    report["freeze_receipt_sha256"] = _validate_runtime_freeze_receipt(root)
    return report


def _read_mhd_w_bcc(path: Path) -> dict[str, Any]:
    sys.path.insert(0, str(REPO_ROOT / "vis/python"))
    import bin_convert_new as bin_convert  # noqa: PLC0415

    return bin_convert.read_binary_as_athdf(str(path))


def _read_particle_vtk(path: Path) -> Any:
    try:
        from .pvtk_particles import read_particle_vtk  # noqa: PLC0415
    except ImportError:
        from pvtk_particles import read_particle_vtk  # type: ignore[no-redef]  # noqa: PLC0415

    return read_particle_vtk(path)


def _x1_average(label: str, payload: Any) -> np.ndarray:
    array = np.asarray(payload, dtype=np.float64)
    if array.ndim < 1 or array.shape[-1] != 32 or not np.all(np.isfinite(array)):
        raise ContractError(f"{label} must retain a finite 32-cell x1 payload")
    axes = tuple(range(array.ndim - 1))
    return np.mean(array, axis=axes) if axes else array


def _branch_spectrum(mhd: dict[str, Any]) -> list[dict[str, Any]]:
    fields = {
        name: _x1_average(name, mhd[name])
        for name in ("bcc2", "bcc3", "vely", "velz")
        if name in mhd
    }
    if set(fields) != {"bcc2", "bcc3", "vely", "velz"}:
        raise ContractError("mhd_w_bcc replay artifact is missing transverse fields")
    spectra = []
    for direction in (-1, 1):
        by = 0.5 * (fields["bcc2"] - direction * fields["vely"])
        bz = 0.5 * (fields["bcc3"] - direction * fields["velz"])
        by_fft = np.fft.rfft(by) / by.size
        bz_fft = np.fft.rfft(bz) / bz.size
        for polarization in (-1, 1):
            for mode in range(1, WAVE_MODE_COUNT + 1):
                coefficient = by_fft[mode] + 1.0j * polarization * bz_fft[mode]
                spectra.append({
                    "mode": mode,
                    "direction": direction,
                    "signed_polarization": polarization,
                    "power": float(abs(coefficient) ** 2),
                })
    return spectra


def _validate_initial_mhd_payload(mhd: dict[str, Any], case: str) -> dict[str, float]:
    """Validate every cycle-zero MHD cell before reducing to branch power."""
    expected_case = _CASES[case]
    length = float(expected_case["x1max"])
    expected_shape = (1, 4, 32)
    expected_keys = {
        "MaxLevel", "NumCycles", "Time",
        "bcc1", "bcc2", "bcc3", "dens", "velx", "vely", "velz",
        "x1f", "x1v", "x2f", "x2v", "x3f", "x3v",
    }
    if set(mhd) != expected_keys:
        raise ContractError("initial mhd_w_bcc replay field topology drifted")
    if int(mhd.get("MaxLevel", -1)) != 0:
        raise ContractError("initial mhd_w_bcc replay refinement level drifted")
    x1v = (np.arange(32, dtype=np.float64) + 0.5) * (length / 32.0)
    coordinates = {
        "x1f": np.linspace(0.0, length, 33),
        "x1v": x1v,
        "x2f": np.arange(5, dtype=np.float64),
        "x2v": np.arange(4, dtype=np.float64) + 0.5,
        "x3f": np.arange(2, dtype=np.float64),
        "x3v": np.array([0.5]),
    }
    state = np.array([
        four_branch_alfven_state(
            x1, length, float(expected_case["paper_seed_amplitude"]),
            int(expected_case["wave_seed"]),
        )
        for x1 in x1v
    ])
    expected_fields = {
        "dens": np.ones(32),
        "velx": np.full(32, float(expected_case["source_local_gas_vx"])),
        "vely": state[:, 2],
        "velz": state[:, 3],
        "bcc1": np.ones(32),
        "bcc2": state[:, 0],
        "bcc3": state[:, 1],
    }
    if float(mhd.get("Time", math.nan)) != 0.0 or int(mhd.get("NumCycles", -1)) != 0:
        raise ContractError("initial mhd_w_bcc replay time or cycle drifted")
    for name, expected in coordinates.items():
        if not np.array_equal(np.asarray(mhd.get(name)), expected):
            raise ContractError(f"initial mhd_w_bcc replay {name} coordinates drifted")
    errors = {}
    for name, expected_x1 in expected_fields.items():
        if name not in mhd:
            raise ContractError(f"initial mhd_w_bcc replay is missing {name}")
        measured = np.asarray(mhd[name], dtype=np.float64)
        if measured.shape != expected_shape:
            raise ContractError(f"initial mhd_w_bcc replay {name} shape drifted")
        expected = np.broadcast_to(expected_x1, measured.shape)
        if not np.all(np.isfinite(measured)) or not np.allclose(
            measured, expected, rtol=0.0, atol=5.0e-10
        ):
            raise ContractError(f"initial mhd_w_bcc replay {name} cells drifted")
        errors[f"{name}_max_abs_error"] = float(np.max(np.abs(measured - expected)))
    return errors


def _validate_initial_branch_spectrum(spectrum: list[dict[str, Any]]) -> None:
    for branch in spectrum:
        expected = wave_branch_amplitude(1.0e-3, branch["mode"]) ** 2
        if not math.isclose(branch["power"], expected, rel_tol=1.0e-4,
                            abs_tol=1.0e-12):
            raise ContractError("initial four-branch wave spectrum drifted")


@lru_cache(maxsize=None)
def _expected_startup_loading(case: str) -> dict[str, np.ndarray]:
    """Return vectorized source-local startup values for one retained case."""
    expected_case = _CASES[case]
    tags = np.arange(SOURCE_LOCAL_PARTICLE_TOTAL, dtype=np.int64)
    species = tags % MOMENTUM_BIN_COUNT
    position_index = tags // MOMENTUM_BIN_COUNT
    angular_sample = position_index // (SOURCE_LOCAL_CELL_COUNT)
    pair = angular_sample // 2
    signs = np.where(angular_sample % 2 == 0, 1.0, -1.0)
    draws = np.array([
        [
            deterministic_uniform01(int(expected_case["momentum_seed"]), index, 0),
            deterministic_uniform01(int(expected_case["momentum_seed"]), index, 1),
        ]
        for index in range(PARTICLES_PER_CELL_PER_BIN // 2)
    ])
    mu = 2.0 * draws[pair, 0] - 1.0
    phi = 2.0 * math.pi * draws[pair, 1]
    p0 = float(expected_case["p0"])
    kappa = float(expected_case["kappa"])
    xi = float(expected_case["xi"])
    light_speed = float(expected_case["light_speed"])
    shell_centers = np.array([
        log_bin_center(p0, index) for index in range(MOMENTUM_BIN_COUNT)
    ])[species]
    perpendicular = shell_centers * np.sqrt(1.0 - mu * mu)
    states = np.column_stack((
        signs * shell_centers * mu,
        signs * perpendicular * np.cos(phi) / xi,
        signs * perpendicular * np.sin(phi) / xi,
    ))
    gamma = np.sqrt(1.0 + np.sum(states * states, axis=1) / (light_speed * light_speed))
    velocities = states / gamma[:, np.newaxis]
    weights = np.array([
        log_bin_fraction(p0, kappa, index) / PARTICLES_PER_CELL_PER_BIN
        for index in range(MOMENTUM_BIN_COUNT)
    ])[species]
    shape = 1.0 + (
        states[:, 0] * states[:, 0]
        + xi * xi * (states[:, 1] * states[:, 1] + states[:, 2] * states[:, 2])
    ) / (kappa * p0 * p0)
    f0 = np.maximum(shape ** (-kappa - 1.0), 1.0e-30)
    cell_linear = position_index % SOURCE_LOCAL_CELL_COUNT
    points = np.column_stack((
        (cell_linear % 32 + 0.5) * (float(expected_case["x1max"]) / 32.0),
        ((cell_linear // 32) % 4 + 0.5),
        np.zeros(SOURCE_LOCAL_PARTICLE_TOTAL),
    ))
    return {
        "tags": tags,
        "species": species,
        "states": states,
        "velocities": velocities,
        "weights": weights,
        "f0": f0,
        "shell_centers": shell_centers,
        "points": points,
    }


def _validate_startup_particle_payload(particles: Any, case: str) -> dict[str, float]:
    """Validate retained cycle-zero positions, angular sampling, shells, and f0."""
    expected = _expected_startup_loading(case)
    points = np.asarray(particles.points, dtype=np.float64)
    velocities = np.asarray(particles.vectors["vel"], dtype=np.float64)
    f0 = np.asarray(particles.scalars["deltaf_f0"], dtype=np.float64)
    if points.shape != expected["points"].shape or not np.all(np.isfinite(points)):
        raise ContractError("initial particle replay points must retain the bounded layout")
    if velocities.shape != expected["velocities"].shape or not np.all(np.isfinite(velocities)):
        raise ContractError("initial particle replay velocities must remain finite 3-vectors")
    expected_vtk_velocities = expected["velocities"].astype(np.float32).astype(np.float64)
    if not np.allclose(points, expected["points"], rtol=0.0, atol=1.0e-6):
        raise ContractError("initial particle replay center-distribution points drifted")
    if not np.allclose(velocities, expected_vtk_velocities, rtol=1.0e-6, atol=1.0e-6):
        raise ContractError("initial particle replay deterministic angular sampler drifted")
    if not np.allclose(f0, expected["f0"], rtol=1.0e-6, atol=1.0e-30):
        raise ContractError("initial particle replay delta-f f0 shape drifted")

    light_speed = float(_CASES[case]["light_speed"])
    speed2 = np.sum(velocities * velocities, axis=1)
    if np.any(speed2 >= light_speed * light_speed):
        raise ContractError("initial particle replay velocity exceeds the artificial light speed")
    states = velocities / np.sqrt(
        1.0 - speed2 / (light_speed * light_speed)
    )[:, np.newaxis]
    expected_vtk_speed2 = np.sum(expected_vtk_velocities * expected_vtk_velocities, axis=1)
    expected_vtk_states = expected_vtk_velocities / np.sqrt(
        1.0 - expected_vtk_speed2 / (light_speed * light_speed)
    )[:, np.newaxis]
    xi = float(_CASES[case]["xi"])
    shell_norms = np.sqrt(
        states[:, 0] * states[:, 0]
        + xi * xi * (states[:, 1] * states[:, 1] + states[:, 2] * states[:, 2])
    )
    expected_vtk_shell_norms = np.sqrt(
        expected_vtk_states[:, 0] * expected_vtk_states[:, 0]
        + xi * xi * (
            expected_vtk_states[:, 1] * expected_vtk_states[:, 1]
            + expected_vtk_states[:, 2] * expected_vtk_states[:, 2]
        )
    )
    if not np.allclose(shell_norms, expected_vtk_shell_norms, rtol=1.0e-5, atol=1.0e-5):
        raise ContractError("initial particle replay transformed shell serialization drifted")
    paired = velocities.reshape(
        PARTICLES_PER_CELL_PER_BIN,
        SOURCE_LOCAL_CELL_COUNT,
        MOMENTUM_BIN_COUNT,
        3,
    )
    pair_sums = paired[0::2] + paired[1::2]
    if not np.allclose(pair_sums, 0.0, rtol=0.0, atol=1.0e-8):
        raise ContractError("initial particle replay antipodal pairing drifted")
    return {
        "point_layout_max_abs_error": float(np.max(np.abs(points - expected["points"]))),
        "velocity_max_abs_error": float(
            np.max(np.abs(velocities - expected["velocities"]))
        ),
        "vtk_serialized_velocity_max_abs_error": float(
            np.max(np.abs(velocities - expected_vtk_velocities))
        ),
        "deltaf_f0_max_abs_error": float(np.max(np.abs(f0 - expected["f0"]))),
        "transformed_shell_serialization_max_abs_error": float(
            np.max(np.abs(shell_norms - expected_vtk_shell_norms))
        ),
        "vtk_inferred_shell_center_max_abs_error": float(
            np.max(np.abs(shell_norms - expected["shell_centers"]))
        ),
        "antipodal_velocity_pair_sum_max_abs": float(np.max(np.abs(pair_sums))),
    }


def _particle_payload(
    path: Path,
    *,
    case: str,
    require_initial_df_zero: bool,
    validate_startup_loading: bool,
) -> dict[str, Any]:
    particles = _read_particle_vtk(path)
    required = {"ptag", "species", "macro_weight", "deltaf_f0", "deltaf_weight"}
    if not required.issubset(particles.scalars) or "vel" not in particles.vectors:
        raise ContractError("particle replay artifact is missing Q-007 fields")
    tags = np.asarray(particles.scalars["ptag"])
    species = np.asarray(particles.scalars["species"])
    weights = np.asarray(particles.scalars["macro_weight"], dtype=np.float64)
    f0 = np.asarray(particles.scalars["deltaf_f0"], dtype=np.float64)
    df = np.asarray(particles.scalars["deltaf_weight"], dtype=np.float64)
    expected = _expected_startup_loading(case)
    if tags.shape != expected["tags"].shape or not np.array_equal(tags, expected["tags"]):
        raise ContractError("particle replay tags do not match bounded serial loading")
    if species.shape != expected["species"].shape or not np.array_equal(
        species, expected["species"]
    ):
        raise ContractError("particle replay species do not retain round-robin log bins")
    if weights.shape != expected["weights"].shape or f0.shape != expected["f0"].shape:
        raise ContractError("particle replay delta-f scalar shapes drifted")
    if df.shape != expected["f0"].shape:
        raise ContractError("particle replay delta-f weight shape drifted")
    if not np.all(np.isfinite(weights)) or not np.all(np.isfinite(f0)):
        raise ContractError("particle replay weights must remain finite")
    if np.any(weights <= 0.0) or np.any(f0 <= 0.0) or not np.all(np.isfinite(df)):
        raise ContractError("particle replay delta-f payload must remain positive and finite")
    if require_initial_df_zero and np.max(np.abs(df)) != 0.0:
        raise ContractError("initial particle replay delta-f weights must be zero")
    counts = [int(np.count_nonzero(species == index)) for index in range(MOMENTUM_BIN_COUNT)]
    if counts != [32768] * 8:
        raise ContractError("particle replay species counts do not match eight log bins")
    if not np.allclose(weights, expected["weights"], rtol=1.0e-6, atol=1.0e-12):
        raise ContractError("particle replay IPWT bin weights drifted")
    if not math.isclose(float(np.sum(weights)), 128.0, rel_tol=1.0e-6):
        raise ContractError("particle replay IPWT normalization drifted")
    payload = {
        "count": int(tags.size),
        "species_counts": counts,
        "macro_weight_sum": float(np.sum(weights)),
        "macro_weight_range": [float(np.min(weights)), float(np.max(weights))],
        "deltaf_f0_range": [float(np.min(f0)), float(np.max(f0))],
        "deltaf_weight_range": [float(np.min(df)), float(np.max(df))],
    }
    if validate_startup_loading:
        payload["startup_loading_validation"] = _validate_startup_particle_payload(
            particles, case
        )
    return payload


def _runtime_replay_tree_root(path: Path) -> Path:
    """Accept the complete retained tree or its historical CRSI child path."""
    root = _authorized_orion_runtime_root(path)
    required_children = {"crsi", "crpai-prolate", "crpai-oblate"}
    if required_children.issubset(child.name for child in root.iterdir() if child.is_dir()):
        return root
    if root.name == "crsi":
        parent = _authorized_orion_runtime_root(root.parent)
        if required_children.issubset(
            child.name for child in parent.iterdir() if child.is_dir()
        ):
            return parent
    raise ContractError("Q-007 runtime replay requires CRSI and both CRPAI startup trees")


def _relative_file_hashes(
    root: Path, subtree: Path, *, excluded: tuple[str, ...] = ()
) -> dict[str, str]:
    """Hash every regular payload below one retained subtree."""
    if not subtree.is_dir():
        raise ContractError(f"retained payload subtree is missing: {subtree}")
    payload = {}
    for path in sorted(subtree.rglob("*")):
        if path.is_file():
            relative = path.relative_to(subtree).as_posix()
            if relative not in excluded:
                payload[relative] = _sha256(_contained_regular_file(root, path))
    return payload


def _validate_pinned_executable_binding(root: Path) -> dict[str, Any]:
    """Validate the immutable source-archived executable used for the replay."""
    path = _contained_regular_file(root, root / PINNED_EXECUTABLE_BINDING_NAME)
    binding = json.loads(path.read_text(encoding="utf-8"))
    if (
        not isinstance(binding, dict)
        or set(binding) != PINNED_EXECUTABLE_BINDING_KEYS
        or binding.get("schema_version") != 1
    ):
        raise ContractError("Q-007 pinned executable binding schema drifted")
    expected = {
        "artifact_role": "bounded_serial_crsi_runtime_replay_mechanics_only",
        "qualification_effect": QUALIFICATION_EFFECT,
        "negative_crpai_nlim1_disposition":
            "fail_closed_time_nlim_contract_error_pending_handedness_review",
    }
    for name, value in expected.items():
        if binding.get(name) != value:
            raise ContractError(f"Q-007 pinned executable binding {name} drifted")
    for name in ("pinned_executable_root", "pinned_executable_root_inventory_sha256",
                 "pinned_executable_path", "pinned_executable_sha256"):
        if not isinstance(binding.get(name), str) or not binding[name]:
            raise ContractError(f"Q-007 pinned executable binding requires {name}")
    for name in ("mpi_used", "slurm_used", "frontier_used", "kronos_used"):
        if binding.get(name) is not False:
            raise ContractError(f"Q-007 pinned executable binding requires {name}=false")
    pinned_root = _authorized_orion_runtime_root(Path(binding["pinned_executable_root"]))
    if binding["pinned_executable_path"] != str(pinned_root / "bin/athena"):
        raise ContractError("Q-007 pinned executable binding path drifted")
    verify_immutable_tree(
        pinned_root,
        binding["pinned_executable_root_inventory_sha256"],
        authorized_root=ORION_PIC_ROOT,
        error_type=ContractError,
        label="Q-007 pinned executable tree",
    )
    executable = _contained_regular_file(pinned_root, Path(binding["pinned_executable_path"]))
    validate_executable_elf(
        executable,
        binding["pinned_executable_sha256"],
        error_type=ContractError,
        label="Q-007 pinned executable",
    )
    receipt_path = _contained_regular_file(pinned_root, pinned_root / FREEZE_RECEIPT_NAME)
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    expected_receipt = {
        "schema_version": 1,
        "artifact_role": "integrated_serial_host_binary_clean_build_provenance_pin",
        "qualification_effect": "bounded_serial_host_provenance_only",
        "inventory_excludes": INVENTORY_NAME,
        "freeze_policy": "remove all owner, group and other write bits recursively",
        "mpi_used": False,
        "slurm_used": False,
        "frontier_used": False,
        "kronos_used": False,
    }
    if receipt != expected_receipt:
        raise ContractError("Q-007 pinned executable receipt drifted")
    provenance_path = _contained_regular_file(
        pinned_root, pinned_root / PIN_PROVENANCE_RECEIPT_NAME
    )
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    source = provenance.get("source_archive")
    if not isinstance(source, dict) or source.get("path") != (
        "source/exact_build_worktree_source.tar.gz"
    ):
        raise ContractError("Q-007 pinned executable provenance source archive path drifted")
    archive = _contained_regular_file(pinned_root, pinned_root / source["path"])
    if source.get("sha256") != _sha256(archive):
        raise ContractError("Q-007 pinned executable provenance source archive SHA-256 drifted")
    evidence = provenance.get("retained_evidence")
    if not isinstance(evidence, list) or set(evidence) != PIN_RETAINED_EVIDENCE:
        raise ContractError("Q-007 pinned executable provenance retained evidence drifted")
    expected_provenance = {
        "schema_version": 1,
        "artifact_role": expected_receipt["artifact_role"],
        "qualification_effect": expected_receipt["qualification_effect"],
        "clean_empty_directory_configure": True,
        "configure_returncode": 0,
        "initial_build_returncode": 0,
        "compile_commands_refresh_returncode": 0,
        "verbose_clean_rebuild_returncode": 0,
        "executable": {"path": "bin/athena", "sha256": binding["pinned_executable_sha256"]},
        "source_archive": {
            "path": "source/exact_build_worktree_source.tar.gz",
            "sha256": _sha256(archive),
        },
        "retained_evidence": sorted(PIN_RETAINED_EVIDENCE),
    }
    if provenance != expected_provenance:
        raise ContractError("Q-007 pinned executable provenance receipt drifted")
    for relative in evidence:
        if not isinstance(relative, str) or not relative:
            raise ContractError("Q-007 pinned executable provenance evidence path is invalid")
        _contained_regular_file(pinned_root, pinned_root / relative)
    expected_pin_files = PIN_RETAINED_EVIDENCE | {
        INVENTORY_NAME, FREEZE_RECEIPT_NAME, PIN_PROVENANCE_RECEIPT_NAME, "bin/athena",
    }
    measured_pin_files = {
        item.relative_to(pinned_root).as_posix()
        for item in pinned_root.rglob("*") if item.is_file()
    }
    if measured_pin_files != expected_pin_files:
        raise ContractError("Q-007 pinned executable tree file topology drifted")
    measured_pin_directories = {
        item.relative_to(pinned_root).as_posix()
        for item in pinned_root.rglob("*") if item.is_dir()
    }
    if measured_pin_directories != {"bin", "build", "runtime", "source"}:
        raise ContractError("Q-007 pinned executable tree directory topology drifted")
    archive_report = validate_source_archive(
        archive,
        source["sha256"],
        error_type=ContractError,
        label="Q-007 pinned executable source archive",
    )
    archive_validation = json.loads(_contained_regular_file(
        pinned_root, pinned_root / "source/archive_validation.json"
    ).read_text(encoding="utf-8"))
    expected_archive_validation = {
        **archive_report,
        "schema_version": 1,
        "policy": (
            "canonical_relative_regular_or_directory_members_no_duplicates_no_links_"
            "no_specials_no_bytecode_cache_no_line_separator_names"
        ),
    }
    if archive_validation != expected_archive_validation:
        raise ContractError("Q-007 pinned executable source archive validation receipt drifted")
    dependency_manifest = json.loads(_contained_regular_file(
        pinned_root, pinned_root / "source/shared_dependency_sha256.json"
    ).read_text(encoding="utf-8"))
    validate_source_archive_dependencies(
        archive,
        dependency_manifest,
        PIN_BUILD_DEPENDENCY_PATHS,
        error_type=ContractError,
        label="Q-007 pinned executable archived dependencies",
    )
    validate_serial_host_build_evidence(
        pinned_root,
        error_type=ContractError,
        label="Q-007 pinned executable serial-host build evidence",
    )
    preflight = json.loads(_contained_regular_file(
        pinned_root, pinned_root / "build/clean_directory_preflight.json"
    ).read_text(encoding="utf-8"))
    if not (
        isinstance(preflight, dict)
        and isinstance(preflight.get("build_directory"), str)
        and Path(preflight["build_directory"]).is_absolute()
        and {name: value for name, value in preflight.items() if name != "build_directory"}
        == {
            "schema_version": 1,
            "existed_before_creation": False,
            "entries_immediately_after_creation": [],
            "clean_empty_directory_configure": True,
        }
    ):
        raise ContractError("Q-007 pinned executable clean-directory preflight receipt drifted")
    return {
        **binding,
        "binding_sha256": _sha256(path),
        "freeze_receipt_sha256": _sha256(receipt_path),
        "clean_build_provenance_receipt_sha256": _sha256(provenance_path),
    }


def _validate_runtime_invocations(root: Path, pinned: dict[str, Any]) -> dict[str, Any]:
    """Validate the three positive commands and the retained CRPAI negative control."""
    executable = pinned["pinned_executable_path"]
    deck_root = root / "decks"
    cases = {
        "crsi": {
            "case": "crsi",
            "deck": deck_root / "pic_q007_paper_crsi_linear_preparation.athinput",
            "overrides": [],
            "returncode": 0,
            "output_cycles": range(4),
        },
        "crpai-prolate": {
            "case": "crpai_prolate",
            "deck": deck_root / "pic_q007_paper_crpai_linear_prolate_preparation.athinput",
            "overrides": [],
            "returncode": 0,
            "output_cycles": range(2),
        },
        "crpai-oblate": {
            "case": "crpai_oblate",
            "deck": deck_root / "pic_q007_paper_crpai_linear_oblate_preparation.athinput",
            "overrides": [],
            "returncode": 0,
            "output_cycles": range(2),
        },
        "negative-crpai-nlim1": {
            "case": "crpai_prolate",
            "deck": deck_root / "pic_q007_paper_crpai_linear_prolate_preparation.athinput",
            "overrides": ["time/nlim=1"],
            "returncode": 1,
            "output_cycles": range(0),
        },
    }
    reports = {}
    for label, expected in cases.items():
        run = root / label
        path = _contained_regular_file(root, run / "invocation.json")
        invocation = json.loads(path.read_text(encoding="utf-8"))
        if (
            not isinstance(invocation, dict)
            or set(invocation) != RUNTIME_INVOCATION_KEYS
            or invocation.get("schema_version") != 1
        ):
            raise ContractError(f"{label}: Q-007 invocation schema drifted")
        if invocation.get("launch_style") != "direct_serial_host_execution":
            raise ContractError(f"{label}: Q-007 invocation launch style drifted")
        for name in ("mpi_used", "slurm_used", "frontier_used", "kronos_used"):
            if invocation.get(name) is not False:
                raise ContractError(f"{label}: Q-007 invocation requires {name}=false")
        identity = {
            "pinned_executable_root": pinned["pinned_executable_root"],
            "pinned_executable_root_inventory_sha256":
                pinned["pinned_executable_root_inventory_sha256"],
            "executable_realpath": executable,
            "executable_sha256": pinned["pinned_executable_sha256"],
        }
        for name, value in identity.items():
            if invocation.get(name) != value:
                raise ContractError(f"{label}: Q-007 invocation {name} drifted")
        deck = _contained_regular_file(root, expected["deck"])
        source_deck = DECKS[expected["case"]]
        if deck.read_bytes() != source_deck.read_bytes():
            raise ContractError(f"{label}: Q-007 retained deck drifted from workspace source")
        validate_deck(deck, expected["case"])
        overrides = expected["overrides"]
        expected_argv = [executable, "-i", str(deck), *overrides]
        if invocation.get("argv") != expected_argv:
            raise ContractError(f"{label}: Q-007 invocation argv drifted")
        if invocation.get("cwd") != str(run):
            raise ContractError(f"{label}: Q-007 invocation cwd drifted")
        if invocation.get("deck") != {"path": str(deck), "sha256": _sha256(deck)}:
            raise ContractError(f"{label}: Q-007 invocation deck binding drifted")
        if invocation.get("overrides") != overrides:
            raise ContractError(f"{label}: Q-007 invocation overrides drifted")
        if invocation.get("returncode") != expected["returncode"]:
            raise ContractError(f"{label}: Q-007 invocation returncode drifted")
        returncode_path = _contained_regular_file(root, run / "returncode.txt")
        if returncode_path.read_text(encoding="utf-8") != f"{expected['returncode']}\n":
            raise ContractError(f"{label}: Q-007 retained returncode text drifted")
        if invocation.get("selected_environment") != {"PYTHONDONTWRITEBYTECODE": "1"}:
            raise ContractError(f"{label}: Q-007 selected environment drifted")
        basename = _CASES[expected["case"]]["basename"]
        expected_payload = {"athena.stderr.txt", "athena.stdout.txt", "returncode.txt"}
        for cycle in expected["output_cycles"]:
            expected_payload.add(f"bin/{basename}.mhd_w_bcc.{cycle:05d}.bin")
            expected_payload.add(f"pvtk/{basename}.prtcl_all.{cycle:05d}.part.vtk")
        measured_payload = _relative_file_hashes(root, run, excluded=("invocation.json",))
        if set(measured_payload) != expected_payload:
            raise ContractError(f"{label}: Q-007 produced payload topology drifted")
        if invocation.get("produced_file_sha256") != measured_payload:
            raise ContractError(f"{label}: Q-007 produced payload membership or hash drifted")
        if label == "negative-crpai-nlim1":
            output = _contained_regular_file(root, run / "athena.stdout.txt").read_text(
                encoding="utf-8", errors="replace"
            )
            if "<time>/nlim does not match the Q-007 preparation contract" not in output:
                raise ContractError("Q-007 retained negative-control diagnostic drifted")
        reports[label] = {"invocation_sha256": _sha256(path)}
    return reports


def _validate_runtime_tree_topology(root: Path) -> None:
    """Reject retained files or directories outside the four closed replay invocations."""
    expected = {
        INVENTORY_NAME,
        FREEZE_RECEIPT_NAME,
        PINNED_EXECUTABLE_BINDING_NAME,
        "decks/pic_q007_paper_crsi_linear_preparation.athinput",
        "decks/pic_q007_paper_crpai_linear_prolate_preparation.athinput",
        "decks/pic_q007_paper_crpai_linear_oblate_preparation.athinput",
    }
    layouts = {
        "crsi": ("crsi", range(4)),
        "crpai-prolate": ("crpai_prolate", range(2)),
        "crpai-oblate": ("crpai_oblate", range(2)),
        "negative-crpai-nlim1": ("crpai_prolate", range(0)),
    }
    for label, (case, cycles) in layouts.items():
        expected.update({
            f"{label}/athena.stderr.txt",
            f"{label}/athena.stdout.txt",
            f"{label}/invocation.json",
            f"{label}/returncode.txt",
        })
        basename = _CASES[case]["basename"]
        for cycle in cycles:
            expected.add(f"{label}/bin/{basename}.mhd_w_bcc.{cycle:05d}.bin")
            expected.add(f"{label}/pvtk/{basename}.prtcl_all.{cycle:05d}.part.vtk")
    measured = {
        item.relative_to(root).as_posix()
        for item in root.rglob("*") if item.is_file()
    }
    if measured != expected:
        raise ContractError("Q-007 retained runtime tree file topology drifted")
    expected_directories = {
        "crsi", "crsi/bin", "crsi/pvtk",
        "crpai-prolate", "crpai-prolate/bin", "crpai-prolate/pvtk",
        "crpai-oblate", "crpai-oblate/bin", "crpai-oblate/pvtk",
        "negative-crpai-nlim1", "decks",
    }
    measured_directories = {
        item.relative_to(root).as_posix()
        for item in root.rglob("*") if item.is_dir()
    }
    if measured_directories != expected_directories:
        raise ContractError("Q-007 retained runtime tree directory topology drifted")


def extract_runtime_replay(
    runtime_root: Path, expected_inventory_sha256: str
) -> dict[str, Any]:
    """Extract retained startup semantics plus the bounded serial CRSI replay."""
    root = _runtime_replay_tree_root(runtime_root)
    tree = verify_frozen_runtime_tree(root, expected_inventory_sha256)
    pinned = _validate_pinned_executable_binding(root)
    invocations = _validate_runtime_invocations(root, pinned)
    _validate_runtime_tree_topology(root)
    crsi = "pic_q007_paper_crsi_linear_preparation"
    prolate = "pic_q007_paper_crpai_linear_prolate_preparation"
    oblate = "pic_q007_paper_crpai_linear_oblate_preparation"
    paths = {
        "crsi_initial_mhd": root / "crsi" / "bin" / f"{crsi}.mhd_w_bcc.00000.bin",
        "crsi_final_mhd": root / "crsi" / "bin" / f"{crsi}.mhd_w_bcc.00002.bin",
        "crsi_initial_particles":
            root / "crsi" / "pvtk" / f"{crsi}.prtcl_all.00000.part.vtk",
        "crsi_final_particles":
            root / "crsi" / "pvtk" / f"{crsi}.prtcl_all.00002.part.vtk",
        "crpai_prolate_initial_mhd":
            root / "crpai-prolate" / "bin" / f"{prolate}.mhd_w_bcc.00000.bin",
        "crpai_prolate_initial_particles":
            root / "crpai-prolate" / "pvtk" / f"{prolate}.prtcl_all.00000.part.vtk",
        "crpai_oblate_initial_mhd":
            root / "crpai-oblate" / "bin" / f"{oblate}.mhd_w_bcc.00000.bin",
        "crpai_oblate_initial_particles":
            root / "crpai-oblate" / "pvtk" / f"{oblate}.prtcl_all.00000.part.vtk",
    }
    for label, path in paths.items():
        try:
            paths[label] = _contained_regular_file(root, path)
        except FileNotFoundError as error:
            raise ContractError(f"Q-007 runtime replay artifact is missing: {label}") from error
    startup_cases = {}
    for case in ("crsi", "crpai_prolate", "crpai_oblate"):
        mhd = _read_mhd_w_bcc(paths[f"{case}_initial_mhd"])
        mhd_validation = _validate_initial_mhd_payload(mhd, case)
        spectrum = _branch_spectrum(mhd)
        _validate_initial_branch_spectrum(spectrum)
        startup_cases[case] = {
            "initial_time": float(mhd["Time"]),
            "initial_mhd_validation": mhd_validation,
            "initial_branch_spectrum": spectrum,
            "initial_particles": _particle_payload(
                paths[f"{case}_initial_particles"],
                case=case,
                require_initial_df_zero=True,
                validate_startup_loading=True,
            ),
        }
    final_mhd = _read_mhd_w_bcc(paths["crsi_final_mhd"])
    return {
        "schema_version": 1,
        "gate": "Q-007",
        "artifact_role": "bounded_serial_crsi_runtime_replay_mechanics_only",
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "runtime_growth_fit_claimed": False,
        "artifact_root": str(root),
        "tree_freeze": tree,
        "pinned_executable": pinned,
        "runtime_invocations": invocations,
        "artifacts": {
            label: {"path": str(path.relative_to(root)), "sha256": _sha256(path)}
            for label, path in paths.items()
        },
        "startup_cases": startup_cases,
        "initial_time": startup_cases["crsi"]["initial_time"],
        "final_time": float(final_mhd["Time"]),
        "initial_branch_spectrum": startup_cases["crsi"]["initial_branch_spectrum"],
        "final_branch_spectrum": _branch_spectrum(final_mhd),
        "initial_particles": startup_cases["crsi"]["initial_particles"],
        "final_particles": _particle_payload(
            paths["crsi_final_particles"],
            case="crsi",
            require_initial_df_zero=False,
            validate_startup_loading=False,
        ),
        "dispersion_oracle": full_q1_q2_dispersion_oracle(),
        "interpretation":
            "two_cycle_serial_replay_only_without_growth_fit_or_qualification_credit",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--runtime-root", type=Path)
    parser.add_argument("--freeze-runtime-tree", type=Path)
    parser.add_argument("--verify-frozen-runtime-tree", type=Path)
    parser.add_argument("--expected-inventory-sha256")
    args = parser.parse_args()
    actions = [
        args.runtime_root,
        args.freeze_runtime_tree,
        args.verify_frozen_runtime_tree,
    ]
    if sum(action is not None for action in actions) > 1:
        raise ContractError("Q-007 analyzer accepts only one runtime-tree action")
    if args.freeze_runtime_tree is not None:
        report = freeze_runtime_tree(args.freeze_runtime_tree)
    elif args.verify_frozen_runtime_tree is not None:
        if args.expected_inventory_sha256 is None:
            raise ContractError("Q-007 frozen-tree verification requires an anchored digest")
        report = verify_frozen_runtime_tree(
            args.verify_frozen_runtime_tree, args.expected_inventory_sha256
        )
    elif args.runtime_root is not None:
        if args.expected_inventory_sha256 is None:
            raise ContractError("Q-007 runtime extraction requires an anchored digest")
        report = extract_runtime_replay(args.runtime_root, args.expected_inventory_sha256)
    else:
        report = build_preparation_report()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
