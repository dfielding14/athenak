#!/usr/bin/env python3
"""Bounded serial-host Q-006 exact-isothermal runtime-local mechanics audit."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import contextvars
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
import sys
from typing import Any

import numpy as np

if __package__:
    from .immutable_orion_tree import authorized_tree_root
    from .immutable_orion_tree import freeze_tree as freeze_immutable_tree
    from .immutable_orion_tree import is_sealed_snapshot_member
    from .immutable_orion_tree import require_exact_primitive_types
    from .immutable_orion_tree import staged_verified_frozen_tree
    from .immutable_orion_tree import validate_executable_elf
    from .immutable_orion_tree import validate_serial_host_build_evidence
    from .immutable_orion_tree import validate_source_archive
    from .immutable_orion_tree import validate_source_archive_dependencies
    from .immutable_orion_tree import verify_frozen_tree as verify_immutable_tree
    from .pvtk_particles import ParticleVTKData, read_particle_vtk
else:
    from immutable_orion_tree import authorized_tree_root
    from immutable_orion_tree import freeze_tree as freeze_immutable_tree
    from immutable_orion_tree import is_sealed_snapshot_member
    from immutable_orion_tree import require_exact_primitive_types
    from immutable_orion_tree import staged_verified_frozen_tree
    from immutable_orion_tree import validate_executable_elf
    from immutable_orion_tree import validate_serial_host_build_evidence
    from immutable_orion_tree import validate_source_archive
    from immutable_orion_tree import validate_source_archive_dependencies
    from immutable_orion_tree import verify_frozen_tree as verify_immutable_tree
    from pvtk_particles import ParticleVTKData, read_particle_vtk


REPO_ROOT = Path(__file__).resolve().parents[2]
ORION_BULK_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
CAMPAIGN_ID = "Q006-PAPER-MULTISPECIES-OSCILLATION-RUNTIME-LOCAL"
PGEN_NAME = "q006_paper_multispecies_oscillation_runtime_local"
PGEN_METHOD = "Q006PaperMultispeciesOscillationRuntimeLocal"
ARTIFACT_ROLE = "bounded_serial_host_exact_isothermal_mechanics_only"
QUALIFICATION_EFFECT = "none"
INVENTORY_NAME = "artifact_inventory.sha256"
FREEZE_RECEIPT_NAME = "freeze_receipt.json"
PINNED_EXECUTABLE_BINDING_NAME = "pinned_executable_binding.json"
PIN_PROVENANCE_RECEIPT_NAME = "build/provenance_receipt.json"
PARSER_CONTRACT_SUMMARY_NAME = "parser_contract_suite/summary.json"
PROBE_SUMMARY_NAME = "reports/probe_summary.json"
PARSER_GUARDS_SOURCE = REPO_ROOT / "tst/scripts/particles/pic_parser_contract_guards.py"
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
    "pinned_executable_root",
    "pinned_executable_root_inventory_sha256",
    "produced_file_sha256",
    "restart",
    "returncode",
    "schema_version",
    "selected_environment",
    "slurm_used",
}
PARSER_INVOCATION_KEYS = {
    "argv",
    "cwd",
    "deck",
    "executable_realpath",
    "executable_sha256",
    "frontier_used",
    "generated_file_sha256",
    "kronos_used",
    "launch_style",
    "mpi_used",
    "pinned_executable_root",
    "pinned_executable_root_inventory_sha256",
    "returncode",
    "schema_version",
    "slurm_used",
}
_PVTK_EXECUTION_PATTERN = re.compile(
    rb"^# vtk DataFile Version 2\.0\r?\n"
    rb"# AthenaK particle data at time=\s*([^ \r\n]+)\s+"
    rb"nranks=\s*([0-9]+)\s+cycle=([0-9]+)\s+variables=([^\r\n]+)\r?\n"
)
SOURCE = REPO_ROOT / "src/pgen/tests/q006_paper_multispecies_oscillation_runtime_local.cpp"
PARTICLES_SOURCE = REPO_ROOT / "src/particles/particles.cpp"
DISPATCH = REPO_ROOT / "src/pgen/pgen.cpp"
HEADER = REPO_ROOT / "src/pgen/pgen.hpp"
CMAKE = REPO_ROOT / "src/CMakeLists.txt"
PREPARATION_SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q006_paper_multispecies_oscillation_source_local_preparation_2026-05-30.json"
)
DECKS = {
    "uniform": (
        REPO_ROOT
        / "inputs/tests/pic_q006_paper_multispecies_oscillation_uniform_runtime_local.athinput"
    ),
    "smr": (
        REPO_ROOT
        / "inputs/tests/pic_q006_paper_multispecies_oscillation_smr_runtime_local.athinput"
    ),
    "audited_amr_runtime_local": (
        REPO_ROOT
        / "inputs/tests/"
        "pic_q006_paper_multispecies_oscillation_audited_amr_runtime_local.athinput"
    ),
}
EXPECTED_DECK_SHA256 = {
    "uniform": "5b8d81aa4914247ba82f065c791d736e7f041ebcb56bc409fd312e17b39f8f01",
    "smr": "ec94c9ac5dddcdfc600e0831b2d6a2a39238b34e42da09203948d1eb3c3fb5eb",
    "audited_amr_runtime_local":
        "9baa6cad1995c359f86401be2ce69580d9b5292a5959925b9abadb1d4f4b3847",
}

_EXPECTED_COMMON = {
    ("time", "integrator"): "rk2",
    ("time", "cfl_number"): "0.1",
    ("time", "nlim"): "2",
    ("time", "tlim"): "1.0",
    ("mhd", "eos"): "isothermal",
    ("mhd", "iso_sound_speed"): "1.0",
    ("mhd", "reconstruct"): "plm",
    ("mhd", "rsolver"): "llf",
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
    ("particles", "couple_moments_energy_to_mhd"): "false",
    ("particles", "couple_fluid_feedback_order"): "mhd_src_terms",
    ("particles", "couple_moments_momentum_coeff"): "1.0",
    ("particles", "couple_moments_energy_coeff"): "0.0",
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
    ("species0", "vy0"): "0.1",
    ("species1", "mass"): "1.0",
    ("species1", "charge"): "1.0",
    ("species1", "vy0"): "0.1",
    ("problem", "pgen_name"): PGEN_NAME,
    ("q006_paper_multispecies_oscillation_runtime_local", "campaign_id"):
        CAMPAIGN_ID,
    ("q006_paper_multispecies_oscillation_runtime_local", "qualification_effect"):
        QUALIFICATION_EFFECT,
    ("q006_paper_multispecies_oscillation_runtime_local", "frontier_authorization"):
        "not_bound",
    ("q006_paper_multispecies_oscillation_runtime_local", "runtime_scope"):
        "bounded_serial_host_mechanics_only",
    ("q006_paper_multispecies_oscillation_runtime_local", "long_horizon_qualification"):
        "not_claimed",
    ("q006_paper_multispecies_oscillation_runtime_local", "true_amr_policy_qualification"):
        "not_claimed",
    ("q006_paper_multispecies_oscillation_runtime_local", "mpi_qualification"):
        "not_claimed",
    ("q006_paper_multispecies_oscillation_runtime_local", "gpu_qualification"):
        "not_claimed",
    ("q006_paper_multispecies_oscillation_runtime_local", "external_review"):
        "not_claimed",
    ("q006_paper_multispecies_oscillation_runtime_local", "runtime_eos_contract"):
        "exact_isothermal_cs_1_fullf_momentum_only",
    ("q006_paper_multispecies_oscillation_runtime_local", "rho"): "1.0",
    ("q006_paper_multispecies_oscillation_runtime_local", "paper_cs"): "1.0",
    ("q006_paper_multispecies_oscillation_runtime_local", "b_g"): "1.0",
    ("q006_paper_multispecies_oscillation_runtime_local", "omega"): "1.0",
    ("q006_paper_multispecies_oscillation_runtime_local", "species_mass_density_ratio"):
        "1.5",
    ("q006_paper_multispecies_oscillation_runtime_local", "gas_uy"): "-0.3",
    ("q006_paper_multispecies_oscillation_runtime_local", "species_vy"): "0.1",
    ("q006_paper_multispecies_oscillation_runtime_local", "artificial_c"): "1000.0",
    ("q006_paper_multispecies_oscillation_runtime_local", "dt_omega_target"): "0.1",
    ("output1", "variable"): "mhd_w_bcc",
    ("output2", "variable"): "prtcl_all",
    ("output3", "file_type"): "rst",
}
_EXPECTED_GRID = {
    "uniform": {
        "deck_role": "bounded_serial_host_uniform_mechanics_only",
        "amr_policy": "not_applicable_uniform",
        "refinement": "none",
        "num_levels": "1",
    },
    "smr": {
        "deck_role": "bounded_serial_host_smr_mechanics_only",
        "amr_policy": "paper_static_one_eighth_region",
        "refinement": "static",
        "num_levels": "2",
    },
    "audited_amr_runtime_local": {
        "deck_role": "bounded_serial_host_audited_amr_mechanics_only",
        "amr_policy":
            "deterministic_audited_randomized_runtime_local_10_refine_60_derefine_"
            "not_qualified",
        "refinement": "adaptive",
        "num_levels": "2",
    },
}
_RETAINED_PROBE_CASES = {
    "uniform_startup_cycle0": {
        "run": "uniform_startup",
        "basename": "q006_uniform_startup",
        "file_index": 0,
        "cycle": 0,
        "meshblocks": 16,
        "minimum_level": 0,
        "maximum_level": 0,
        "particle_count": 131072,
        "startup_anchor": True,
    },
    "uniform_evolution_cycle2": {
        "run": "uniform_evolution",
        "basename": "q006_uniform_evolution",
        "file_index": 3,
        "cycle": 2,
        "meshblocks": 16,
        "minimum_level": 0,
        "maximum_level": 0,
        "particle_count": 131072,
    },
    "smr_evolution_cycle2": {
        "run": "smr_evolution",
        "basename": "q006_smr_evolution",
        "file_index": 3,
        "cycle": 2,
        "meshblocks": 72,
        "minimum_level": 0,
        "maximum_level": 1,
        "particle_count": 589824,
    },
    "audited_amr_evolution_cycle2": {
        "run": "audited_amr_evolution",
        "basename": "q006_audited_amr_evolution",
        "file_index": 3,
        "cycle": 2,
        "meshblocks": 44,
        "minimum_level": 0,
        "maximum_level": 1,
        "particle_count": 131072,
    },
    "audited_amr_restart_cycle3": {
        "run": "audited_amr_restart",
        "basename": "q006_audited_amr_restart",
        "file_index": 5,
        "cycle": 3,
        "meshblocks": 58,
        "minimum_level": 0,
        "maximum_level": 1,
        "particle_count": 131072,
    },
}


class AuditError(ValueError):
    """Raised when a bounded Q-006 runtime-local contract fails closed."""


_STAGED_TREES: contextvars.ContextVar[tuple[tuple[Path, Any], ...]] = (
    contextvars.ContextVar("q006_staged_trees", default=())
)


@contextmanager
def _use_staged_tree(logical_root: Path, staged_tree: Any):
    """Route retained-tree reads through one descriptor-anchored sealed snapshot."""
    token = _STAGED_TREES.set((*_STAGED_TREES.get(), (logical_root, staged_tree)))
    try:
        yield
    finally:
        _STAGED_TREES.reset(token)


def _tree_file_io_path(path: Path) -> Path:
    """Route one retained regular file through the active sealed snapshot."""
    for logical_root, staged_tree in reversed(_STAGED_TREES.get()):
        routed = staged_tree.file_io_path(logical_root, path)
        if routed != path:
            return routed
    return path


def _tree_directory_io_path(path: Path) -> Path:
    """Route one retained directory through the active topology snapshot."""
    for logical_root, staged_tree in reversed(_STAGED_TREES.get()):
        routed = staged_tree.directory_io_path(logical_root, path)
        if routed != path:
            return routed
    return path


def _has_staged_tree(root: Path) -> bool:
    return any(logical_root == root for logical_root, _ in _STAGED_TREES.get())


def _staged_tree_for(path: Path) -> tuple[Path, Any, Path] | None:
    for logical_root, staged_tree in reversed(_STAGED_TREES.get()):
        try:
            return logical_root, staged_tree, path.relative_to(logical_root)
        except ValueError:
            pass
        try:
            return logical_root, staged_tree, path.relative_to(staged_tree.staged_root)
        except ValueError:
            pass
    return None


def _tree_relative_files(subtree: Path) -> set[str]:
    staged = _staged_tree_for(subtree)
    if staged is None:
        return {
            path.relative_to(subtree).as_posix()
            for path in subtree.rglob("*")
            if path.is_file()
        }
    _, staged_tree, relative = staged
    return staged_tree.relative_files(relative)


def _tree_relative_directories(subtree: Path) -> set[str]:
    staged = _staged_tree_for(subtree)
    if staged is None:
        return {
            path.relative_to(subtree).as_posix()
            for path in subtree.rglob("*")
            if path.is_dir()
        }
    _, staged_tree, relative = staged
    return staged_tree.relative_directories(relative)


def _tree_immediate_directories(subtree: Path) -> set[str]:
    staged = _staged_tree_for(subtree)
    if staged is None:
        return {path.name for path in subtree.iterdir() if path.is_dir()}
    _, staged_tree, relative = staged
    return staged_tree.immediate_directories(relative)


def _tree_has_directory(subtree: Path) -> bool:
    staged = _staged_tree_for(subtree)
    if staged is None:
        return subtree.is_dir()
    _, staged_tree, relative = staged
    return staged_tree.has_directory(relative)


def _restore_logical_paths(value: Any) -> Any:
    """Restore Orion path labels after parsing private staged bytes."""
    if isinstance(value, dict):
        return {key: _restore_logical_paths(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_restore_logical_paths(item) for item in value]
    if isinstance(value, str):
        for logical_root, staged_tree in reversed(_STAGED_TREES.get()):
            restored = staged_tree.logical_path_for_io(logical_root, value)
            if restored is not None:
                return str(restored)
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def _require_exact_match(actual: Any, expected: Any, label: str) -> None:
    """Reject primitive aliases before requiring an exact retained JSON value."""
    require_exact_primitive_types(actual, expected, error_type=AuditError, label=label)
    _require(actual == expected, f"{label} drifted")


def _require_exact_int(actual: Any, expected: int, label: str) -> None:
    """Require one JSON integer without accepting boolean or float aliases."""
    _require(
        type(actual) is int and type(expected) is int and actual == expected,
        f"{label} drifted",
    )


def _require_canonical_returncode_sidecar(path: Path, expected: int, label: str) -> None:
    """Require the one canonical textual encoding for a retained return code."""
    _require(
        type(expected) is int and path.read_text(encoding="utf-8") == f"{expected}\n",
        f"{label} sidecar drifted",
    )


def _require_schema_version(value: Any, expected: int, label: str) -> None:
    """Reject Python boolean and float aliases for JSON schema integers."""
    _require(
        isinstance(value, dict)
        and type(value.get("schema_version")) is int
        and value["schema_version"] == expected,
        f"{label}: schema drifted",
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_athinput(path: Path) -> dict[str, dict[str, str]]:
    """Parse the strict Athena input subset used by the runtime-local decks."""
    blocks: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    for lineno, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<"):
            _require(line.endswith(">"), f"{path}:{lineno}: malformed block header")
            name = line[1:-1].strip()
            _require(bool(name), f"{path}:{lineno}: empty block name")
            current = blocks.setdefault(name, {})
            continue
        _require(current is not None and "=" in line,
                 f"{path}:{lineno}: malformed parameter line")
        name, value = (item.strip() for item in line.split("=", 1))
        _require(bool(name) and bool(value),
                 f"{path}:{lineno}: empty parameter name or value")
        _require(name not in current, f"{path}:{lineno}: duplicate parameter {name}")
        current[name] = value
    return blocks


def validate_deck(grid_setup: str, path: Path | None = None) -> dict[str, Any]:
    """Require one exact bounded successor deck without assigning qualification credit."""
    selected = DECKS[grid_setup] if path is None else path
    blocks = parse_athinput(selected)
    for (block, name), expected in _EXPECTED_COMMON.items():
        measured = blocks.get(block, {}).get(name)
        _require(
            measured == expected,
            f"{selected}: {block}/{name}: expected {expected!r}, measured {measured!r}",
        )
    mesh = blocks["mesh"]
    for name, expected in (
        ("nx1", "16"), ("nx2", "8"), ("nx3", "8"),
        ("x1min", "0.0"), ("x1max", "16.0"),
        ("x2min", "0.0"), ("x2max", "8.0"),
        ("x3min", "0.0"), ("x3max", "8.0"),
    ):
        _require(mesh.get(name) == expected, f"{selected}: mesh/{name} drifted")
    for name in ("ix1_bc", "ox1_bc", "ix2_bc", "ox2_bc", "ix3_bc", "ox3_bc"):
        _require(mesh.get(name) == "periodic", f"{selected}: mesh/{name} must be periodic")
    _require(
        tuple(blocks["meshblock"].get(f"nx{axis}") for axis in (1, 2, 3))
        == ("4", "4", "4"),
        f"{selected}: meshblock geometry drifted",
    )
    expected_grid = _EXPECTED_GRID[grid_setup]
    metadata = blocks["q006_paper_multispecies_oscillation_runtime_local"]
    refinement = blocks["mesh_refinement"]
    for name in ("deck_role", "amr_policy"):
        _require(metadata.get(name) == expected_grid[name],
                 f"{selected}: {name} drifted")
    _require(metadata.get("grid_setup") == grid_setup,
             f"{selected}: grid_setup drifted")
    for name in ("refinement", "num_levels"):
        _require(refinement.get(name) == expected_grid[name],
                 f"{selected}: mesh_refinement/{name} drifted")
    if grid_setup == "smr":
        _require(
            blocks.get("refinement1") == {
                "level": "1",
                "x1min": "4.0", "x1max": "12.0",
                "x2min": "2.0", "x2max": "6.0",
                "x3min": "2.0", "x3max": "6.0",
            },
            f"{selected}: one-eighth SMR region drifted",
        )
    if grid_setup == "audited_amr_runtime_local":
        _require(refinement.get("max_nmb_per_rank") == "512",
                 f"{selected}: AMR capacity drifted")
        _require(metadata.get("amr_seed") == "60053",
                 f"{selected}: AMR seed drifted")
    _require(_sha256(selected) == EXPECTED_DECK_SHA256[grid_setup],
             f"{selected}: exact bounded Q-006 deck SHA-256 drifted")
    return {
        "grid_setup": grid_setup,
        "path": str(selected.relative_to(REPO_ROOT)),
        "sha256": _sha256(selected),
        "deck_role": metadata["deck_role"],
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
    }


def validate_decks() -> list[dict[str, Any]]:
    """Validate all three bounded serial-host successor carriers."""
    return [validate_deck(grid) for grid in DECKS]


def validate_historical_preparation_bindings() -> dict[str, Any]:
    """Require every committed preparation artifact to retain its frozen hash."""
    sidecar = json.loads(PREPARATION_SIDECAR.read_text(encoding="utf-8"))
    bindings = sidecar["artifact_bindings"]
    for relative, expected in bindings.items():
        _require(_sha256(REPO_ROOT / relative) == expected,
                 f"historical Q-006 preparation hash drifted: {relative}")
    return {"preserved": True, "artifact_count": len(bindings)}


def validate_registration() -> dict[str, Any]:
    """Require additive registration and the two narrow isothermal parser allowances."""
    source = SOURCE.read_text(encoding="utf-8")
    particles = PARTICLES_SOURCE.read_text(encoding="utf-8")
    dispatch = DISPATCH.read_text(encoding="utf-8")
    header = HEADER.read_text(encoding="utf-8")
    cmake = CMAKE.read_text(encoding="utf-8")
    _require(dispatch.count(f'compare("{PGEN_NAME}")') == 2,
             "Q-006 runtime-local dispatch must be additive for fresh and restart")
    _require(re.search(
        rf"void\s+{re.escape(PGEN_METHOD)}\(\s*ParameterInput\s+\*pin,\s*"
        rf"const\s+bool\s+restart\s*\);", header) is not None,
             "Q-006 runtime-local declaration is absent")
    _require("pgen/tests/q006_paper_multispecies_oscillation_runtime_local.cpp" in cmake,
             "Q-006 runtime-local source is absent from CMake")
    for snippet in (
        "exact_isothermal_deltaf_paper_feedback",
        "UsesDeltaF() && !couple_moments_energy_to_mhd",
        "exact isothermal paper delta-f uses momentum-only",
        "exact_isothermal_fullf_paper_feedback",
        "!UsesDeltaF() && !couple_moments_energy_to_mhd",
        '"q006_paper_multispecies_oscillation_runtime_local"',
        "full-f uses the same momentum-only feedback contract only in",
    ):
        _require(snippet in particles, f"Q-006 parser allowance is missing {snippet!r}")
    for snippet in (
        "Q006PaperMultispeciesOscillationRuntimeLocal",
        "Q006RuntimeLocalAuditedAMRRefinement",
        "bounded_serial_host_mechanics_only",
        "exact_isothermal_cs_1_fullf_momentum_only",
        "Q006RuntimeLocalRequireBoolean(pin, \"particles\", "
        "\"couple_moments_energy_to_mhd\", false);",
        "if (restart) return;",
        "pgen_q006_runtime_local_uniform_isothermal_state",
    ):
        _require(snippet in source, f"Q-006 runtime-local source is missing {snippet!r}")
    _require("w0(m, IEN" not in source,
             "Q-006 exact-isothermal initialization must not write an energy slot")
    return {
        "pgen_name": PGEN_NAME,
        "fresh_dispatch": True,
        "restart_dispatch": True,
        "deltaf_allowance_preserved": True,
        "fullf_allowance_added": True,
    }


def analytical_contract() -> dict[str, float]:
    """Return the manuscript-local Section 5.3 normalization."""
    omega = 1.0 * math.sqrt(1.0 + 2.0 * 1.5)
    return {
        "rho": 1.0,
        "species_mass_density_each": 1.5,
        "species_ppc_each": 64.0,
        "aggregate_ppc": 128.0,
        "deposit_qscale": 0.0234375,
        "gas_uy": -0.3,
        "species_vy": 0.1,
        "oscillation_omega": omega,
        "oscillation_frequency_cycles": omega / (2.0 * math.pi),
        "oscillation_period": 2.0 * math.pi / omega,
        "initial_kinetic_energy_density": 0.06,
    }


def _finite_array(label: str, value: Any) -> np.ndarray:
    array = np.asarray(value)
    _require(np.all(np.isfinite(array)), f"{label} must contain only finite values")
    return array


def analyze_snapshot_arrays(
    mhd_density: Any,
    mhd_velocity: Any,
    cell_volume: Any,
    particle_velocity: Any,
    macro_weight: Any,
    particle_species: Any,
    particle_tags: Any,
) -> dict[str, Any]:
    """Recompute bounded weighted mechanics from one serial-host snapshot."""
    density = _finite_array("MHD density", mhd_density).astype(np.float64)
    velocity = _finite_array("MHD velocity", mhd_velocity).astype(np.float64)
    volumes = _finite_array("cell volume", cell_volume).astype(np.float64)
    _require(velocity.shape == density.shape + (3,), "MHD velocity shape mismatch")
    _require(volumes.shape == density.shape, "cell-volume shape mismatch")
    _require(np.all(density > 0.0), "MHD density must remain positive")
    _require(np.all(volumes > 0.0), "cell volumes must remain positive")
    physical_volume = float(np.sum(volumes))
    _require(math.isclose(physical_volume, 1024.0, rel_tol=0.0, abs_tol=1.0e-10),
             "carrier physical volume drifted")

    particle_v = _finite_array("particle velocity", particle_velocity).astype(np.float64)
    weights = _finite_array("macro weight", macro_weight).astype(np.float64)
    species = np.asarray(particle_species)
    tags = np.asarray(particle_tags)
    _require(particle_v.ndim == 2 and particle_v.shape[1] == 3,
             "particle velocity shape mismatch")
    _require(weights.shape == (particle_v.shape[0],), "macro-weight shape mismatch")
    _require(species.shape == (particle_v.shape[0],), "particle species shape mismatch")
    _require(tags.shape == (particle_v.shape[0],), "particle tag shape mismatch")
    _require(np.all(weights > 0.0), "macro weights must remain positive")
    _require(np.issubdtype(species.dtype, np.integer), "particle species must be integral")
    _require(np.issubdtype(tags.dtype, np.integer), "particle tags must be integral")
    _require(np.all((species == 0) | (species == 1)), "unexpected particle species")
    _require(np.unique(tags).size == tags.size, "particle tags must remain unique")

    fluid_mass_weights = density * volumes
    fluid_momentum = np.sum(fluid_mass_weights[..., None] * velocity, axis=tuple(
        range(density.ndim)
    ))
    macro_mass = analytical_contract()["deposit_qscale"] * weights
    particle_momentum = np.sum(macro_mass[:, None] * particle_v, axis=0)
    total_momentum = fluid_momentum + particle_momentum
    fluid_ke = float(np.sum(
        0.5 * fluid_mass_weights * np.sum(velocity * velocity, axis=-1)
    ))
    particle_ke = float(np.sum(0.5 * macro_mass * np.sum(particle_v * particle_v, axis=1)))
    species_mass_density = [
        float(np.sum(macro_mass[species == index]) / physical_volume) for index in (0, 1)
    ]
    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "long_horizon_qualification_claimed": False,
        "true_amr_policy_qualification_claimed": False,
        "mpi_qualification_claimed": False,
        "gpu_qualification_claimed": False,
        "frontier_qualification_claimed": False,
        "external_review_claimed": False,
        "particle_count": int(tags.size),
        "physical_volume": physical_volume,
        "volume_averaged_gas_velocity": (
            np.sum(fluid_mass_weights[..., None] * velocity, axis=tuple(range(density.ndim)))
            / np.sum(fluid_mass_weights)
        ).tolist(),
        "species_mass_density": species_mass_density,
        "fluid_momentum": fluid_momentum.tolist(),
        "particle_momentum": particle_momentum.tolist(),
        "total_momentum": total_momentum.tolist(),
        "maximum_absolute_total_momentum": float(np.max(np.abs(total_momentum))),
        "fluid_kinetic_energy": fluid_ke,
        "particle_kinetic_energy": particle_ke,
        "total_kinetic_energy": fluid_ke + particle_ke,
        "total_kinetic_energy_density": (fluid_ke + particle_ke) / physical_volume,
    }


def _read_mhd_binary(path: Path) -> dict[str, Any]:
    sys.path.insert(0, str(REPO_ROOT / "vis/python"))
    import bin_convert_new as bin_convert  # noqa: PLC0415

    return bin_convert.read_binary(str(path))


def _read_pvtk_execution_metadata(path: Path) -> dict[str, Any]:
    """Require the serial prtcl_all execution header emitted by AthenaK."""
    with path.open("rb") as stream:
        header = stream.read(4096)
    match = _PVTK_EXECUTION_PATTERN.match(header)
    _require(match is not None, "PVTK execution header is missing or malformed")
    time = float(match.group(1))
    nranks = int(match.group(2))
    cycle = int(match.group(3))
    variables = match.group(4).decode("ascii").strip()
    _require(math.isfinite(time), "PVTK execution time must be finite")
    _require(nranks == 1, "bounded Q-006 PVTK artifact must report nranks=1")
    _require(variables == "prtcl_all",
             "bounded Q-006 PVTK artifact must report variables=prtcl_all")
    return {"time": time, "nranks": nranks, "cycle": cycle, "variables": variables}


def _contained_regular_file(root: Path, path: Path) -> Path:
    """Require one self-contained regular artifact below its retained root."""
    io_root = _tree_directory_io_path(root)
    io_path = _tree_file_io_path(path)
    if is_sealed_snapshot_member(io_path):
        try:
            path.relative_to(root)
        except ValueError as error:
            raise AuditError(f"runtime artifact must remain below retained root: {path}") from error
        _require(io_path != path, f"sealed runtime artifact must be routed by active snapshot: {path}")
        return io_path
    status = os.lstat(io_path)
    _require(stat.S_ISREG(status.st_mode), f"runtime artifact must be regular: {path}")
    _require(status.st_nlink == 1, f"runtime artifact must have one link: {path}")
    resolved = io_path.resolve(strict=True)
    _require(resolved == io_path, f"runtime artifact path must be canonical: {path}")
    try:
        resolved.relative_to(io_root)
    except ValueError as error:
        raise AuditError(f"runtime artifact must remain below retained root: {path}") from error
    return resolved


def _authorized_retained_root(root: str | Path) -> Path:
    return authorized_tree_root(
        root,
        authorized_root=ORION_BULK_ROOT,
        error_type=AuditError,
        label="Q-006 retained probe root",
    )


def _validate_runtime_freeze_receipt(root: Path) -> str:
    """Require the retained Q-006 tree to carry its exact nonqualifying role."""
    path = _contained_regular_file(root, root / FREEZE_RECEIPT_NAME)
    receipt = json.loads(path.read_text(encoding="utf-8"))
    _require_schema_version(receipt, 1, "Q-006 freeze receipt")
    expected = {
        "schema_version": 1,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "inventory_excludes": INVENTORY_NAME,
        "freeze_policy": "remove all owner, group and other write bits recursively",
    }
    for name, value in expected.items():
        _require(receipt.get(name) == value, f"Q-006 freeze receipt {name} drifted")
    return _sha256(path)


def extract_runtime_snapshot(
    mhd_path: str | Path,
    pvtk_path: str | Path,
    *,
    retained_root: str | Path | None = None,
) -> dict[str, Any]:
    """Decode and hash one immutable raw MHD plus particle snapshot pair."""
    if retained_root is None:
        mhd_file = Path(mhd_path).resolve()
        particle_file = Path(pvtk_path).resolve()
    else:
        root = _authorized_retained_root(retained_root)
        mhd_file = _contained_regular_file(root, Path(mhd_path))
        particle_file = _contained_regular_file(root, Path(pvtk_path))
    mhd = _read_mhd_binary(mhd_file)
    pvtk_execution = _read_pvtk_execution_metadata(particle_file)
    particles: ParticleVTKData = read_particle_vtk(particle_file)
    for name in ("dens", "velx", "vely", "velz", "bcc1", "bcc2", "bcc3"):
        _require(name in mhd["mb_data"], f"MHD artifact is missing {name}")
    for name in ("macro_weight", "species", "ptag"):
        _require(name in particles.scalars, f"particle artifact is missing {name}")
    _require("vel" in particles.vectors, "particle artifact is missing vel")
    points = _finite_array("particle position", particles.points).astype(np.float64)
    _require(points.ndim == 2 and points.shape[1] == 3,
             "particle position shape mismatch")
    _require(
        np.all((points[:, 0] >= float(mhd["x1min"])) & (points[:, 0] <= float(mhd["x1max"])))
        and np.all((points[:, 1] >= float(mhd["x2min"])) & (points[:, 1] <= float(mhd["x2max"])))
        and np.all((points[:, 2] >= float(mhd["x3min"])) & (points[:, 2] <= float(mhd["x3max"]))),
        "particle positions left the retained carrier domain",
    )
    density = np.asarray(mhd["mb_data"]["dens"], dtype=np.float64)
    mhd_velocity = np.stack(
        [np.asarray(mhd["mb_data"][name], dtype=np.float64)
         for name in ("velx", "vely", "velz")],
        axis=-1,
    )
    levels = np.asarray(mhd["mb_logical"], dtype=np.int64)[:, 3]
    magnetic = np.stack(
        [np.asarray(mhd["mb_data"][name], dtype=np.float64)
         for name in ("bcc1", "bcc2", "bcc3")],
        axis=-1,
    )
    _require(np.all(np.isfinite(magnetic)), "MHD magnetic field must remain finite")
    _require(density.shape[0] == levels.size, "MHD MeshBlock payload shape mismatch")
    block_volumes = np.power(0.5, 3 * levels, dtype=np.float64)
    cell_volume = np.broadcast_to(
        block_volumes.reshape((-1,) + (1,) * (density.ndim - 1)),
        density.shape,
    )
    report = analyze_snapshot_arrays(
        density,
        mhd_velocity,
        cell_volume,
        particles.vectors["vel"],
        particles.scalars["macro_weight"],
        particles.scalars["species"],
        particles.scalars["ptag"],
    )
    report["snapshot"] = {
        "mhd_time": float(mhd["time"]),
        "mhd_cycle": int(mhd["cycle"]),
        "meshblock_count": int(levels.size),
        "minimum_level": int(np.min(levels)),
        "maximum_level": int(np.max(levels)),
    }
    _require(report["snapshot"]["mhd_cycle"] == pvtk_execution["cycle"],
             "MHD and PVTK cycles do not match")
    _require(math.isclose(report["snapshot"]["mhd_time"], pvtk_execution["time"],
                          rel_tol=0.0, abs_tol=1.0e-12),
             "MHD and PVTK times do not match")
    report["pvtk_execution_metadata"] = pvtk_execution
    report["magnetic_field"] = {
        "component_minima": np.min(magnetic, axis=tuple(range(magnetic.ndim - 1))).tolist(),
        "component_maxima": np.max(magnetic, axis=tuple(range(magnetic.ndim - 1))).tolist(),
    }
    report["immutable_artifacts"] = {
        "mhd_w_bcc": {"path": str(mhd_file), "sha256": _sha256(mhd_file)},
        "prtcl_all": {"path": str(particle_file), "sha256": _sha256(particle_file)},
    }
    return report


def _validate_pinned_executable_binding(root: Path) -> dict[str, Any]:
    """Verify the source-archived executable lineage retained beside the probe."""
    binding_path = _contained_regular_file(root, root / PINNED_EXECUTABLE_BINDING_NAME)
    initial_binding = json.loads(binding_path.read_text(encoding="utf-8"))
    initial_pinned_root = _authorized_retained_root(initial_binding["pinned_executable_root"])
    if not _has_staged_tree(initial_pinned_root):
        with staged_verified_frozen_tree(
            initial_pinned_root,
            initial_binding["pinned_executable_root_inventory_sha256"],
            authorized_root=ORION_BULK_ROOT,
            error_type=AuditError,
            label="Q-006 pinned executable tree",
        ) as (_, staged_root), _use_staged_tree(initial_pinned_root, staged_root):
            return _validate_pinned_executable_binding(root)
    path = _contained_regular_file(root, root / PINNED_EXECUTABLE_BINDING_NAME)
    binding = json.loads(path.read_text(encoding="utf-8"))
    _require(isinstance(binding, dict) and set(binding) == PINNED_EXECUTABLE_BINDING_KEYS,
             "pinned executable binding schema drifted")
    _require_schema_version(binding, 1, "pinned executable binding")
    expected_binding = {
        "schema_version": 1,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
    }
    for name, value in expected_binding.items():
        _require(binding.get(name) == value, f"pinned executable binding {name} drifted")
    for name in ("pinned_executable_root", "pinned_executable_root_inventory_sha256",
                 "pinned_executable_path", "pinned_executable_sha256"):
        _require(isinstance(binding.get(name), str) and binding[name],
                 f"pinned executable binding requires {name}")
    for name in ("mpi_used", "slurm_used", "frontier_used", "kronos_used"):
        _require(binding.get(name) is False, f"pinned executable binding requires {name}=false")
    pinned_root = _authorized_retained_root(binding["pinned_executable_root"])
    _require(binding["pinned_executable_path"] == str(pinned_root / "bin/athena"),
             "pinned executable binding path drifted")
    executable = _contained_regular_file(pinned_root, Path(binding["pinned_executable_path"]))
    validate_executable_elf(
        executable,
        binding["pinned_executable_sha256"],
        error_type=AuditError,
        label="Q-006 pinned executable",
    )
    receipt_path = _contained_regular_file(pinned_root, pinned_root / FREEZE_RECEIPT_NAME)
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    _require_schema_version(receipt, 1, "pinned executable receipt")
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
    _require_exact_match(receipt, expected_receipt, "pinned executable receipt")
    provenance_path = _contained_regular_file(
        pinned_root, pinned_root / PIN_PROVENANCE_RECEIPT_NAME
    )
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    _require_schema_version(provenance, 1, "pinned executable provenance receipt")
    source = provenance.get("source_archive")
    _require(isinstance(source, dict) and source.get("path")
             == "source/exact_build_worktree_source.tar.gz",
             "pinned executable provenance source archive path drifted")
    archive = _contained_regular_file(pinned_root, pinned_root / source["path"])
    _require(source.get("sha256") == _sha256(archive),
             "pinned executable provenance source archive SHA-256 drifted")
    evidence = provenance.get("retained_evidence")
    _require(isinstance(evidence, list) and set(evidence) == PIN_RETAINED_EVIDENCE,
             "pinned executable provenance retained evidence drifted")
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
    _require_exact_match(provenance, expected_provenance, "pinned executable provenance receipt")
    for relative in evidence:
        _require(isinstance(relative, str) and relative,
                 "pinned executable provenance evidence path is invalid")
        _contained_regular_file(pinned_root, pinned_root / relative)
    expected_pin_files = PIN_RETAINED_EVIDENCE | {
        INVENTORY_NAME, FREEZE_RECEIPT_NAME, PIN_PROVENANCE_RECEIPT_NAME, "bin/athena",
    }
    io_pinned_root = _tree_directory_io_path(pinned_root)
    measured_pin_files = _tree_relative_files(pinned_root)
    _require(measured_pin_files == expected_pin_files,
             "pinned executable tree file topology drifted")
    measured_pin_directories = _tree_relative_directories(pinned_root)
    _require(measured_pin_directories == {"bin", "build", "runtime", "source"},
             "pinned executable tree directory topology drifted")
    archive_report = validate_source_archive(
        archive,
        source["sha256"],
        error_type=AuditError,
        label="Q-006 pinned executable source archive",
    )
    archive_validation = json.loads(_contained_regular_file(
        pinned_root, pinned_root / "source/archive_validation.json"
    ).read_text(encoding="utf-8"))
    _require_schema_version(archive_validation, 1, "pinned executable archive validation")
    expected_archive_validation = {
        **archive_report,
        "schema_version": 1,
        "policy": (
            "canonical_relative_regular_or_directory_members_no_duplicates_no_links_"
            "no_specials_no_bytecode_cache_no_line_separator_names"
        ),
    }
    _require_exact_match(
        archive_validation,
        expected_archive_validation,
        "pinned executable source archive validation receipt",
    )
    dependency_manifest = json.loads(_contained_regular_file(
        pinned_root, pinned_root / "source/shared_dependency_sha256.json"
    ).read_text(encoding="utf-8"))
    validate_source_archive_dependencies(
        archive,
        dependency_manifest,
        PIN_BUILD_DEPENDENCY_PATHS,
        error_type=AuditError,
        label="Q-006 pinned executable archived dependencies",
    )
    validate_serial_host_build_evidence(
        io_pinned_root,
        file_resolver=lambda relative: _tree_file_io_path(pinned_root / relative),
        error_type=AuditError,
        label="Q-006 pinned executable serial-host build evidence",
    )
    preflight = json.loads(_contained_regular_file(
        pinned_root, pinned_root / "build/clean_directory_preflight.json"
    ).read_text(encoding="utf-8"))
    _require_schema_version(preflight, 1, "pinned executable clean-directory preflight")
    _require(
        isinstance(preflight, dict)
        and isinstance(preflight.get("build_directory"), str)
        and Path(preflight["build_directory"]).is_absolute(),
        "pinned executable clean-directory preflight receipt drifted",
    )
    _require_exact_match(
        {name: value for name, value in preflight.items() if name != "build_directory"},
        {
            "schema_version": 1,
            "existed_before_creation": False,
            "entries_immediately_after_creation": [],
            "clean_empty_directory_configure": True,
        },
        "pinned executable clean-directory preflight receipt",
    )
    return {
        **binding,
        "binding_sha256": _sha256(path),
        "freeze_receipt_sha256": _sha256(receipt_path),
        "clean_build_provenance_receipt_sha256": _sha256(provenance_path),
    }


def _load_rooted_json(root: Path, path: Path) -> tuple[Path, Any]:
    """Load one root-contained JSON sidecar after enforcing regular-file policy."""
    rooted = _contained_regular_file(root, path)
    return rooted, json.loads(rooted.read_text(encoding="utf-8"))


def _relative_file_hashes(
    root: Path, subtree: Path, *, excluded: tuple[str, ...] = ()
) -> dict[str, str]:
    """Hash every regular payload below one retained subtree."""
    relative_files = _tree_relative_files(subtree)
    _require(_tree_has_directory(subtree), f"retained payload subtree is missing: {subtree}")
    payload = {}
    for relative in sorted(relative_files):
        if relative not in excluded:
            payload[relative] = _sha256(_contained_regular_file(root, subtree / relative))
    return payload


def _validate_serial_invocation_identity(
    invocation: dict[str, Any], pinned: dict[str, Any], *, label: str
) -> None:
    """Require direct serial-host execution of the exact pinned executable."""
    _require_schema_version(invocation, 1, f"{label}: invocation")
    _require(invocation.get("launch_style") == "direct_serial_host_execution",
             f"{label}: invocation launch style drifted")
    for name in ("mpi_used", "slurm_used", "frontier_used", "kronos_used"):
        _require(invocation.get(name) is False, f"{label}: invocation requires {name}=false")
    expected = {
        "pinned_executable_root": pinned["pinned_executable_root"],
        "pinned_executable_root_inventory_sha256":
            pinned["pinned_executable_root_inventory_sha256"],
        "executable_realpath": pinned["pinned_executable_path"],
        "executable_sha256": pinned["pinned_executable_sha256"],
    }
    for name, value in expected.items():
        _require(invocation.get(name) == value, f"{label}: invocation {name} drifted")


def _validate_runtime_invocations(root: Path, pinned: dict[str, Any]) -> dict[str, Any]:
    """Validate all five fixed runtime command manifests and their raw payloads."""
    decks = {
        "uniform": root / "decks/pic_q006_paper_multispecies_oscillation_uniform_runtime_local.athinput",
        "smr": root / "decks/pic_q006_paper_multispecies_oscillation_smr_runtime_local.athinput",
        "audited_amr": root / "decks/pic_q006_paper_multispecies_oscillation_audited_amr_runtime_local.athinput",
    }
    restart = (
        root
        / "runs/audited_amr_evolution/rst/q006_audited_amr_evolution.00003.rst"
    )
    executable = pinned["pinned_executable_path"]
    cases = {
        "uniform_startup": {
            "deck": decks["uniform"],
            "argv": [executable, "-i", str(decks["uniform"]), "time/nlim=0",
                     "job/basename=q006_uniform_startup"],
        },
        "uniform_evolution": {
            "deck": decks["uniform"],
            "argv": [executable, "-i", str(decks["uniform"]),
                     "job/basename=q006_uniform_evolution"],
        },
        "smr_evolution": {
            "deck": decks["smr"],
            "argv": [executable, "-i", str(decks["smr"]),
                     "job/basename=q006_smr_evolution"],
        },
        "audited_amr_evolution": {
            "deck": decks["audited_amr"],
            "argv": [executable, "-i", str(decks["audited_amr"]),
                     "job/basename=q006_audited_amr_evolution"],
        },
        "audited_amr_restart": {
            "restart": restart,
            "argv": [executable, "-r", str(restart), "time/nlim=3",
                     "job/basename=q006_audited_amr_restart"],
        },
    }
    reports = {}
    run_root = root / "runs"
    _require(
        _tree_immediate_directories(run_root) == set(cases),
        "Q-006 runtime run membership drifted",
    )
    for label, expected in cases.items():
        run = root / "runs" / label
        path, invocation = _load_rooted_json(root, run / "invocation.json")
        _require(isinstance(invocation, dict), f"{label}: invocation must be an object")
        _require(set(invocation) == RUNTIME_INVOCATION_KEYS,
                 f"{label}: invocation schema drifted")
        _validate_serial_invocation_identity(invocation, pinned, label=label)
        _require(invocation.get("cwd") == str(run), f"{label}: invocation cwd drifted")
        _require(invocation.get("argv") == expected["argv"], f"{label}: invocation argv drifted")
        _require_exact_int(invocation.get("returncode"), 0, f"{label}: invocation returncode")
        _require(invocation.get("selected_environment") == {"PYTHONDONTWRITEBYTECODE": "1"},
                 f"{label}: selected environment drifted")
        produced = _relative_file_hashes(root, run, excluded=("invocation.json",))
        _require(invocation.get("produced_file_sha256") == produced,
            f"{label}: produced payload membership or hash drifted",
        )
        _require(set(produced) == _expected_runtime_payload_paths(label),
                 f"{label}: produced payload topology drifted")
        if "deck" in expected:
            deck = _contained_regular_file(root, expected["deck"])
            _require(invocation.get("restart") is None, f"{label}: unexpected restart record")
            _require(invocation.get("deck") == {
                "path": str(expected["deck"]),
                "sha256": _sha256(deck),
            },
                     f"{label}: deck binding drifted")
        else:
            source = _contained_regular_file(root, expected["restart"])
            _require(invocation.get("deck") is None, f"{label}: restart deck must be null")
            _require(
                invocation.get("restart") == {
                    "path": str(expected["restart"]),
                    "sha256": _sha256(source),
                },
                f"{label}: restart binding drifted",
            )
        reports[label] = {"invocation_sha256": _sha256(path)}
    return reports


def _expected_runtime_payload_paths(label: str) -> set[str]:
    """Return the exact bounded output topology for one fixed runtime replay."""
    setup = {
        "uniform_startup": ("q006_uniform_startup", range(2)),
        "uniform_evolution": ("q006_uniform_evolution", range(4)),
        "smr_evolution": ("q006_smr_evolution", range(4)),
        "audited_amr_evolution": ("q006_audited_amr_evolution", range(4)),
        "audited_amr_restart": ("q006_audited_amr_restart", range(4, 6)),
    }
    basename, indices = setup[label]
    payload = {"command.txt", "stderr.log", "stdout.log"}
    for index in indices:
        suffix = f"{index:05d}"
        payload.add(f"bin/{basename}.mhd_w_bcc.{suffix}.bin")
        payload.add(f"pvtk/{basename}.prtcl_all.{suffix}.part.vtk")
        payload.update({
            f"rst/{basename}.{suffix}.rst",
            f"rst/{basename}.{suffix}.rst.complete",
            f"rst/{basename}.{suffix}.rst.manifest",
            f"rst/{basename}.{suffix}.rst.manifest.complete",
        })
    return payload


def _parser_contract_oracle(root: Path, pinned: dict[str, Any]) -> list[dict[str, Any]]:
    """Build the exact parser matrix from the source-owned regression oracle."""
    sys.path.insert(0, str(REPO_ROOT / "tst"))
    from scripts.particles import pic_parser_contract_guards as guards  # noqa: PLC0415

    generic = root / "decks/pic_parser_contract_guards.athinput"
    q006 = root / "decks/pic_q006_paper_multispecies_oscillation_uniform_runtime_local.athinput"
    oracle = []
    for positive, cases in ((True, guards._POSITIVE_CASES), (False, guards._REJECTION_CASES)):
        for label, arguments, expected in cases:
            deck = q006 if label in guards._Q006_RUNTIME_LOCAL_CASES else generic
            oracle.append({
                "label": label,
                "positive": positive,
                "expected": expected,
                "deck": {"path": str(deck), "sha256": _sha256(deck)},
                "argv": [
                    pinned["pinned_executable_path"], "-i", str(deck), "time/nlim=0",
                    *arguments,
                ],
            })
    return oracle


def _expected_parser_generated_paths(label: str, positive: bool) -> set[str]:
    if not positive:
        return set()
    if label != "paper_mhd_pic_isothermal_fullf_momentum_only":
        return {"pic_parser_contract_guards-errs.dat"}
    basename = "pic_q006_paper_multispecies_oscillation_uniform_runtime_local"
    payload = set()
    for index in range(2):
        suffix = f"{index:05d}"
        payload.add(f"bin/{basename}.mhd_w_bcc.{suffix}.bin")
        payload.add(f"pvtk/{basename}.prtcl_all.{suffix}.part.vtk")
        payload.update({
            f"rst/{basename}.{suffix}.rst",
            f"rst/{basename}.{suffix}.rst.complete",
            f"rst/{basename}.{suffix}.rst.manifest",
            f"rst/{basename}.{suffix}.rst.manifest.complete",
        })
    return payload


def _validate_parser_contract_suite(root: Path, pinned: dict[str, Any]) -> dict[str, Any]:
    """Validate frozen parser command manifests, sidecars, and generated payloads."""
    path, summary = _load_rooted_json(root, root / PARSER_CONTRACT_SUMMARY_NAME)
    _require(isinstance(summary, dict), "parser-contract summary must be an object")
    positive = summary.get("positive_cases")
    negative = summary.get("rejection_cases")
    _require(isinstance(positive, list) and len(positive) == 6,
             "parser-contract positive case count drifted")
    _require(isinstance(negative, list) and len(negative) == 31,
             "parser-contract negative case count drifted")
    _require_exact_int(summary.get("accepted_case_count"), 6,
                       "parser-contract accepted case count")
    _require_exact_int(summary.get("rejected_case_count"), 31,
                       "parser-contract rejected case count")
    _require(summary.get("all_passed") is True, "parser-contract suite did not pass")
    _require(
        summary.get("pinned_executable_root_inventory_sha256")
        == pinned["pinned_executable_root_inventory_sha256"],
        "parser-contract pinned executable inventory drifted",
    )
    _require(summary.get("pinned_executable_sha256") == pinned["pinned_executable_sha256"],
             "parser-contract pinned executable SHA-256 drifted")
    oracle = _parser_contract_oracle(root, pinned)
    expected_by_label = {item["label"]: item for item in oracle}
    labels = []
    for item, expected_positive in [
        *((entry, True) for entry in positive),
        *((entry, False) for entry in negative),
    ]:
        _require(isinstance(item, dict), "parser-contract case summary must be an object")
        label = item.get("label")
        _require(isinstance(label, str) and label and "/" not in label and label != "..",
                 "parser-contract case label is unsafe")
        labels.append(label)
        _require(label in expected_by_label, f"{label}: parser-contract label drifted")
        expected_case = expected_by_label[label]
        case = root / "parser_contract_suite/cases" / label
        invocation_path, invocation = _load_rooted_json(root, case / "invocation.json")
        _require(isinstance(invocation, dict), f"{label}: parser invocation must be an object")
        _require(set(invocation) == PARSER_INVOCATION_KEYS,
                 f"{label}: parser invocation schema drifted")
        _validate_serial_invocation_identity(invocation, pinned, label=label)
        work = case / "work"
        _require(invocation.get("cwd") == str(work), f"{label}: parser invocation cwd drifted")
        _require(item.get("argv") == expected_case["argv"],
                 f"{label}: parser summary argv drifted")
        _require(invocation.get("argv") == expected_case["argv"],
                 f"{label}: parser invocation argv drifted")
        argv = invocation.get("argv")
        _require(isinstance(argv, list) and len(argv) >= 4,
                 f"{label}: parser invocation argv is incomplete")
        deck = _contained_regular_file(root, Path(argv[2]))
        _require(argv[:2] == [pinned["pinned_executable_path"], "-i"],
                 f"{label}: parser invocation prefix drifted")
        _require(invocation.get("deck") == expected_case["deck"],
                 f"{label}: parser deck binding drifted")
        _require(invocation.get("deck") == item.get("deck"),
                 f"{label}: parser summary deck binding drifted")
        generated = _relative_file_hashes(root, work)
        _require(invocation.get("generated_file_sha256") == generated,
                 f"{label}: parser generated payload drifted")
        _require(item.get("generated_file_sha256") == generated,
                 f"{label}: parser summary generated payload drifted")
        _require(set(generated) == _expected_parser_generated_paths(label, expected_positive),
                 f"{label}: parser generated payload topology drifted")
        returncode = invocation.get("returncode")
        _require_exact_int(returncode, item.get("returncode"), f"{label}: parser returncode")
        _require_exact_int(
            returncode,
            0 if expected_positive else 1,
            f"{label}: parser positive/negative disposition",
        )
        _require(item.get("positive") is expected_positive and item.get("passed") is True,
                 f"{label}: parser summary disposition drifted")
        stdout = _contained_regular_file(root, case / "stdout.txt").read_text(
            encoding="utf-8", errors="replace"
        )
        stderr = _contained_regular_file(root, case / "stderr.txt").read_text(
            encoding="utf-8", errors="replace"
        )
        _require(item.get("expected") == expected_case["expected"],
                 f"{label}: parser expected diagnostic or token definition drifted")
        expected_text = expected_case["expected"]
        tokens = expected_text if isinstance(expected_text, list) else [expected_text]
        _require(all(isinstance(token, str) and token in stdout + stderr for token in tokens),
                 f"{label}: parser expected diagnostic or token drifted")
        _require(
            json.loads(_contained_regular_file(root, case / "command.json").read_text(
                encoding="utf-8"
            )) == argv,
            f"{label}: parser command sidecar drifted",
        )
        _require_canonical_returncode_sidecar(
            _contained_regular_file(root, case / "returncode.txt"),
            returncode,
            f"{label}: parser returncode",
        )
        _require(
            set(_relative_file_hashes(root, case)) - {"invocation.json"}
            == {"command.json", "returncode.txt", "stderr.txt", "stdout.txt"}
            | {f"work/{name}" for name in generated},
            f"{label}: parser case payload membership drifted",
        )
        _require(_sha256(invocation_path), f"{label}: parser invocation hash is empty")
    _require(len(labels) == len(set(labels)), "parser-contract labels must be unique")
    _require(set(labels) == set(expected_by_label), "parser-contract labels drifted")
    case_root = root / "parser_contract_suite/cases"
    observed = _tree_immediate_directories(case_root)
    _require(observed == set(labels), "parser-contract case membership drifted")
    return {
        "summary_sha256": _sha256(path),
        "positive_case_count": len(positive),
        "negative_case_count": len(negative),
        "all_passed": True,
    }


def _validate_retained_topology(root: Path, pinned: dict[str, Any]) -> None:
    """Reject unexpected retained directories and descriptive-artifact drift."""
    for name, source in {
        "pic_q006_paper_multispecies_oscillation_uniform_runtime_local.athinput":
            DECKS["uniform"],
        "pic_q006_paper_multispecies_oscillation_smr_runtime_local.athinput": DECKS["smr"],
        "pic_q006_paper_multispecies_oscillation_audited_amr_runtime_local.athinput":
            DECKS["audited_amr_runtime_local"],
        "pic_parser_contract_guards.athinput":
            REPO_ROOT / "inputs/tests/pic_parser_contract_guards.athinput",
    }.items():
        retained = _contained_regular_file(root, root / "decks" / name)
        _require(retained.read_bytes() == source.read_bytes(), f"retained deck drifted: {name}")
    convenience = _contained_regular_file(root, root / "bin/athena")
    _require(_sha256(convenience) == pinned["pinned_executable_sha256"],
             "Q-006 convenience executable drifted")
    static_path, static = _load_rooted_json(root, root / "reports/static_descriptor.json")
    _require_exact_match(static, static_descriptor(), "Q-006 retained static descriptor")
    _require(_sha256(static_path), "Q-006 static descriptor hash is empty")
    allowed = {
        "bin", "decks", "reports", "runs", "parser_contract_suite",
        "parser_contract_suite/cases",
    }
    for label in _RETAINED_PROBE_CASES.values():
        run = f"runs/{label['run']}"
        allowed.update({run, f"{run}/bin", f"{run}/pvtk", f"{run}/rst"})
    for item in _parser_contract_oracle(root, pinned):
        case = f"parser_contract_suite/cases/{item['label']}"
        allowed.update({case, f"{case}/work"})
        for relative in _expected_parser_generated_paths(item["label"], item["positive"]):
            parent = (Path(case) / "work" / relative).parent.as_posix()
            while parent != f"{case}/work":
                allowed.add(parent)
                parent = Path(parent).parent.as_posix()
    measured = _tree_relative_directories(root)
    _require(measured == allowed, "Q-006 retained directory topology drifted")
    allowed_files = {
        INVENTORY_NAME,
        FREEZE_RECEIPT_NAME,
        PINNED_EXECUTABLE_BINDING_NAME,
        PROBE_SUMMARY_NAME,
        "reports/static_descriptor.json",
        "bin/athena",
        "decks/pic_q006_paper_multispecies_oscillation_uniform_runtime_local.athinput",
        "decks/pic_q006_paper_multispecies_oscillation_smr_runtime_local.athinput",
        "decks/pic_q006_paper_multispecies_oscillation_audited_amr_runtime_local.athinput",
        "decks/pic_parser_contract_guards.athinput",
        PARSER_CONTRACT_SUMMARY_NAME,
    }
    for label in {
        "uniform_startup", "uniform_evolution", "smr_evolution",
        "audited_amr_evolution", "audited_amr_restart",
    }:
        allowed_files.add(f"runs/{label}/invocation.json")
        allowed_files.update(f"runs/{label}/{relative}"
                             for relative in _expected_runtime_payload_paths(label))
    for item in _parser_contract_oracle(root, pinned):
        case = f"parser_contract_suite/cases/{item['label']}"
        allowed_files.update({
            f"{case}/command.json",
            f"{case}/invocation.json",
            f"{case}/returncode.txt",
            f"{case}/stderr.txt",
            f"{case}/stdout.txt",
        })
        allowed_files.update(
            f"{case}/work/{relative}"
            for relative in _expected_parser_generated_paths(item["label"], item["positive"])
        )
    measured_files = _tree_relative_files(root)
    _require(measured_files == allowed_files, "Q-006 retained file topology drifted")


def build_retained_probe_report(root: str | Path) -> dict[str, Any]:
    """Regenerate the fixed bounded carrier report from root-contained raw pairs."""
    retained = _authorized_retained_root(root)
    reports = {}
    for label, expected in _RETAINED_PROBE_CASES.items():
        prefix = retained / "runs" / expected["run"]
        index = expected["file_index"]
        basename = expected["basename"]
        report = extract_runtime_snapshot(
            prefix / "bin" / f"{basename}.mhd_w_bcc.{index:05d}.bin",
            prefix / "pvtk" / f"{basename}.prtcl_all.{index:05d}.part.vtk",
            retained_root=retained,
        )
        snapshot = report["snapshot"]
        for name in ("cycle", "meshblocks", "minimum_level", "maximum_level"):
            measured_name = "mhd_cycle" if name == "cycle" else (
                "meshblock_count" if name == "meshblocks" else name
            )
            _require(snapshot[measured_name] == expected[name],
                     f"{label}: retained {name} drifted")
        _require(report["particle_count"] == expected["particle_count"],
                 f"{label}: retained particle count drifted")
        _require(np.allclose(report["species_mass_density"], [1.5, 1.5],
                             rtol=0.0, atol=1.0e-12),
                 f"{label}: retained species density drifted")
        if expected.get("startup_anchor"):
            _require(np.allclose(report["volume_averaged_gas_velocity"], [0.0, -0.3, 0.0],
                                 rtol=0.0, atol=1.0e-7),
                     "uniform startup gas velocity drifted")
            _require(report["maximum_absolute_total_momentum"] <= 1.0e-4,
                     "uniform startup total momentum drifted")
            _require(math.isclose(report["total_kinetic_energy_density"], 0.06,
                                  rel_tol=0.0, abs_tol=1.0e-7),
                     "uniform startup kinetic-energy anchor drifted")
            _require(np.allclose(report["magnetic_field"]["component_minima"],
                                 [0.0, 0.0, 1.0], rtol=0.0, atol=1.0e-12),
                     "uniform startup magnetic minima drifted")
            _require(np.allclose(report["magnetic_field"]["component_maxima"],
                                 [0.0, 0.0, 1.0], rtol=0.0, atol=1.0e-12),
                     "uniform startup magnetic maxima drifted")
        reports[label] = report
    pinned = _validate_pinned_executable_binding(retained)
    return _restore_logical_paths({
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "scope": "bounded local mechanics only",
        "pinned_executable": pinned,
        "runtime_invocations": _validate_runtime_invocations(retained, pinned),
        "parser_contract_suite": _validate_parser_contract_suite(retained, pinned),
        "snapshots": reports,
        "open_items": [
            "long_horizon", "true_amr_policy_qualification", "mpi", "gpu",
            "frontier", "external_review",
        ],
    })


def verify_retained_probe(
    root: str | Path, expected_inventory_sha256: str
) -> dict[str, Any]:
    """Verify one frozen tree and exact regeneration of its fixed probe summary."""
    retained = _authorized_retained_root(root)
    with staged_verified_frozen_tree(
        retained,
        expected_inventory_sha256,
        authorized_root=ORION_BULK_ROOT,
        error_type=AuditError,
        label="Q-006 retained runtime tree",
    ) as (tree, staged_root), _use_staged_tree(retained, staged_root):
        tree["freeze_receipt_sha256"] = _validate_runtime_freeze_receipt(retained)
        _validate_retained_topology(retained, _validate_pinned_executable_binding(retained))
        measured = build_retained_probe_report(retained)
        summary = _contained_regular_file(retained, retained / PROBE_SUMMARY_NAME)
        expected = json.loads(summary.read_text(encoding="utf-8"))
        _require_schema_version(expected, 1, "retained Q-006 probe summary")
        _require_exact_match(expected, measured, "retained Q-006 probe summary")
        return {"tree_freeze": tree, "summary_sha256": _sha256(summary), "summary": measured}


def static_descriptor() -> dict[str, Any]:
    """Describe the bounded successor without claiming long-horizon evidence."""
    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "long_horizon_qualification_claimed": False,
        "true_amr_policy_qualification_claimed": False,
        "mpi_qualification_claimed": False,
        "gpu_qualification_claimed": False,
        "frontier_qualification_claimed": False,
        "external_review_claimed": False,
        "decks": validate_decks(),
        "registration": validate_registration(),
        "historical_preparation": validate_historical_preparation_bindings(),
        "analytical_contract": analytical_contract(),
    }


def verify_frozen_tree(root: str | Path, expected_inventory_sha256: str) -> dict[str, Any]:
    """Verify exact inventory hashes, tree membership, and recursive read-only modes."""
    retained = _authorized_retained_root(root)
    with staged_verified_frozen_tree(
        retained,
        expected_inventory_sha256,
        authorized_root=ORION_BULK_ROOT,
        error_type=AuditError,
        label="Q-006 retained runtime tree",
    ) as (report, staged_root), _use_staged_tree(retained, staged_root):
        report["freeze_receipt_sha256"] = _validate_runtime_freeze_receipt(retained)
        return report


def freeze_tree(root: str | Path) -> dict[str, Any]:
    """Write an exact inventory and recursively remove write bits from the tree."""
    receipt = {
        "schema_version": 1,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "inventory_excludes": INVENTORY_NAME,
        "freeze_policy": "remove all owner, group and other write bits recursively",
    }
    return freeze_immutable_tree(
        root,
        receipt,
        authorized_root=ORION_BULK_ROOT,
        error_type=AuditError,
        label="Q-006 retained runtime tree",
    )


def _write_json(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("static")
    snapshot = subparsers.add_parser("snapshot")
    snapshot.add_argument("mhd")
    snapshot.add_argument("particles")
    retained = subparsers.add_parser("build-retained-probe-report")
    retained.add_argument("root")
    retained.add_argument("output")
    retained_verify = subparsers.add_parser("verify-retained-probe")
    retained_verify.add_argument("root")
    retained_verify.add_argument("--expected-inventory-sha256", required=True)
    freeze = subparsers.add_parser("freeze-tree")
    freeze.add_argument("root")
    verify = subparsers.add_parser("verify-frozen-tree")
    verify.add_argument("root")
    verify.add_argument("--expected-inventory-sha256", required=True)
    args = parser.parse_args()
    if args.command == "static":
        _write_json(static_descriptor())
    elif args.command == "snapshot":
        _write_json(extract_runtime_snapshot(args.mhd, args.particles))
    elif args.command == "build-retained-probe-report":
        Path(args.output).write_text(
            json.dumps(build_retained_probe_report(args.root), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    elif args.command == "verify-retained-probe":
        _write_json(verify_retained_probe(args.root, args.expected_inventory_sha256))
    elif args.command == "freeze-tree":
        _write_json(freeze_tree(args.root))
    else:
        _write_json(verify_frozen_tree(args.root, args.expected_inventory_sha256))


if __name__ == "__main__":
    main()
