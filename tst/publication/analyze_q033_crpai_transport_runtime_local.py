#!/usr/bin/env python3
"""Bounded source-local Q-033 CRPAI runtime diagnostic extractor."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ID = "Q033-CRPAI-TRANSPORT-RUNTIME-LOCAL"
PGEN_NAME = "q033_crpai_transport_runtime_local"
RESERVED_BLOCKED_PGEN_NAME = "q033_crpai_transport_calibration_open"
ARTIFACT_ROLE = "bounded_source_local_runtime_diagnostic_only"
QUALIFICATION_EFFECT = "none"
DECK = REPO_ROOT / "inputs/tests/pic_q033_crpai_transport_runtime_local.athinput"
SOURCE = REPO_ROOT / "src/pgen/tests/q033_crpai_transport_runtime_local.cpp"
DISPATCH = REPO_ROOT / "src/pgen/pgen.cpp"
HEADER = REPO_ROOT / "src/pgen/pgen.hpp"
CMAKE = REPO_ROOT / "src/CMakeLists.txt"
WAVE_MODE_COUNT = 3
MASK64 = (1 << 64) - 1

_EXPECTED_DECK_VALUES = {
    ("mesh", "nx1"): "32",
    ("mesh", "nx2"): "4",
    ("mesh", "nx3"): "1",
    ("meshblock", "nx1"): "32",
    ("meshblock", "nx2"): "4",
    ("meshblock", "nx3"): "1",
    ("mesh_refinement", "refinement"): "none",
    ("time", "integrator"): "rk1",
    ("time", "cfl_number"): "0.1",
    ("time", "nlim"): "1",
    ("mhd", "eos"): "ideal",
    ("mhd", "reconstruct"): "plm",
    ("mhd", "rsolver"): "llf",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "ppc"): "4.0",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "nspecies"): "1",
    ("particles", "cr_distribution"): "center",
    ("particles", "deposit_moments"): "true",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "couple_moments_momentum_to_mhd"): "true",
    ("particles", "couple_moments_energy_to_mhd"): "true",
    ("particles", "pic_physical_mode"): "extended_mhd_pic",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_enable_2d3v"): "true",
    ("particles", "pic_cr_initial_state"): "momentum",
    ("particles", "pic_cr_hall_mode"): "off",
    ("particles", "pic_wave_damping_mode"): "ion_neutral_friction",
    ("particles", "pic_ion_neutral_collision_rate"): "1.0e-4",
    ("particles", "pic_deltaf_mode"): "physical",
    ("particles", "pic_deltaf_f0"): "kappa_aniso",
    ("particles", "pic_deltaf_adapt_mode"):
        "global_bikappa_moments_experimental",
    ("particles", "pic_expanding_box_mode"): "on",
    ("particles", "pic_expansion_law"): "exponential",
    ("particles", "track_displacement"): "true",
    ("problem", "pgen_name"): PGEN_NAME,
    ("q033_crpai_transport_runtime_local", "campaign_id"): CAMPAIGN_ID,
    ("q033_crpai_transport_runtime_local", "deck_role"):
        "bounded_source_local_runtime_diagnostic_only_not_transport_calibration",
    ("q033_crpai_transport_runtime_local", "qualification_effect"): "none",
    ("q033_crpai_transport_runtime_local", "frontier_authorization"): "not_bound",
    ("q033_crpai_transport_runtime_local", "q022_closure"): "not_claimed",
    ("q033_crpai_transport_runtime_local", "physical_calibration"): "not_claimed",
    ("q033_crpai_transport_runtime_local", "external_review"): "not_claimed",
    ("q033_crpai_transport_runtime_local", "momentum_distribution"):
        "deterministic_antipodal_bounded_prolate_lattice",
    ("q033_crpai_transport_runtime_local", "momentum_seed"): "330034",
    ("q033_crpai_transport_runtime_local", "momentum_p0"): "1.0",
    ("q033_crpai_transport_runtime_local", "momentum_xi"): "2.0",
    ("q033_crpai_transport_runtime_local", "wave_spectrum"):
        "deterministic_seeded_three_mode_transverse_ct_carrier",
    ("q033_crpai_transport_runtime_local", "wave_seed"): "330033",
    ("q033_crpai_transport_runtime_local", "wave_mode_count"): "3",
    ("q033_crpai_transport_runtime_local", "wave_amplitude"): "1.0e-4",
    ("output2", "variable"): "mhd_bcc",
    ("output3", "variable"): "prtcl_all",
}


class ContractError(ValueError):
    """Raised when a bounded Q-033 runtime-local contract fails closed."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_athinput(path: Path = DECK) -> dict[str, dict[str, str]]:
    """Parse the strict Athena input subset used by the bounded local deck."""
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


def _finite_float(label: str, value: str) -> float:
    try:
        measured = float(value)
    except ValueError as error:
        raise ContractError(f"{label} must be a finite number") from error
    if not math.isfinite(measured):
        raise ContractError(f"{label} must be a finite number")
    return measured


def validate_deck(path: Path = DECK) -> dict[str, Any]:
    """Require the exact bounded local role without assigning qualification credit."""
    blocks = parse_athinput(path)
    for (block, name), expected in _EXPECTED_DECK_VALUES.items():
        measured = blocks.get(block, {}).get(name)
        if measured != expected:
            raise ContractError(
                f"{path}: {block}/{name}: expected {expected!r}, measured {measured!r}"
            )
    for name in ("ix1_bc", "ox1_bc", "ix2_bc", "ox2_bc"):
        if blocks["mesh"].get(name) != "periodic":
            raise ContractError(f"{path}: mesh/{name} must remain periodic")
    metadata = blocks["q033_crpai_transport_runtime_local"]
    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "sha256": _sha256(path),
        "campaign_id": metadata["campaign_id"],
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "physical_calibration_claimed": False,
        "frontier_authorization": metadata["frontier_authorization"],
        "q022_closure": metadata["q022_closure"],
        "external_review": metadata["external_review"],
        "wave_seed": int(metadata["wave_seed"]),
        "wave_amplitude": _finite_float("wave_amplitude", metadata["wave_amplitude"]),
        "momentum_seed": int(metadata["momentum_seed"]),
        "momentum_p0": _finite_float("momentum_p0", metadata["momentum_p0"]),
        "momentum_xi": _finite_float("momentum_xi", metadata["momentum_xi"]),
    }


def validate_registration() -> dict[str, Any]:
    """Require the local successor while preserving the reserved open-name block."""
    dispatch = DISPATCH.read_text(encoding="utf-8")
    header = HEADER.read_text(encoding="utf-8")
    cmake = CMAKE.read_text(encoding="utf-8")
    source = SOURCE.read_text(encoding="utf-8")
    if dispatch.count(f'compare("{PGEN_NAME}")') != 2:
        raise ContractError("Q-033 runtime-local dispatch must be additive for fresh/restart")
    if RESERVED_BLOCKED_PGEN_NAME in dispatch:
        raise ContractError("reserved blocked Q-033 candidate unexpectedly became runnable")
    if "void Q033CRPAITransportRuntimeLocal(ParameterInput *pin," not in header:
        raise ContractError("Q-033 runtime-local declaration is absent")
    if "pgen/tests/q033_crpai_transport_runtime_local.cpp" not in cmake:
        raise ContractError("Q-033 runtime-local source is absent from CMake")
    for snippet in (
        "deterministic_antipodal_bounded_prolate_lattice",
        "deterministic_seeded_three_mode_transverse_ct_carrier",
        "Q033WaveBy",
        "Q033WaveBz",
        "pgen_q033_local_antipodal_prolate_momenta",
        "if (restart) return;",
    ):
        if snippet not in source:
            raise ContractError(f"Q-033 runtime-local source is missing {snippet!r}")
    return {
        "pgen_name": PGEN_NAME,
        "fresh_dispatch": True,
        "restart_dispatch": True,
        "reserved_blocked_pgen_name": RESERVED_BLOCKED_PGEN_NAME,
        "reserved_blocked_pgen_registered": False,
    }


def splitmix64(value: int) -> int:
    """Match the source-local uint64 SplitMix64 implementation."""
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def deterministic_uniform01(seed: int, item: int, component: int) -> float:
    """Match the source-local deterministic diagnostic carrier hash."""
    key = seed & MASK64
    key ^= ((item + 1) * 0xBF58476D1CE4E5B9) & MASK64
    key ^= ((component + 1) * 0xD2B74407B1CE6E93) & MASK64
    return (splitmix64(key) >> 11) / 9007199254740992.0


def build_static_descriptor() -> dict[str, Any]:
    """Describe the deterministic local initialization without physical calibration."""
    deck = validate_deck()
    registration = validate_registration()
    modes = []
    for mode in range(1, WAVE_MODE_COUNT + 1):
        modes.append({
            "mode": mode,
            "phase": 2.0 * math.pi * deterministic_uniform01(
                deck["wave_seed"], mode, 0
            ),
            "amplitude": deck["wave_amplitude"] / (mode * mode),
            "handedness": "positive" if mode % 2 else "negative",
        })
    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "physical_calibration_claimed": False,
        "q022_closure_claimed": False,
        "external_review_claimed": False,
        "frontier_qualification_claimed": False,
        "launch_authorization_claimed": False,
        "deck": deck,
        "registration": registration,
        "wave_seed_contract": {
            "semantics": "deterministic_seeded_three_mode_transverse_ct_carrier",
            "modes": modes,
        },
        "momentum_seed_contract": {
            "semantics": "deterministic_antipodal_bounded_prolate_lattice",
            "seed": deck["momentum_seed"],
            "p0": deck["momentum_p0"],
            "xi": deck["momentum_xi"],
            "formula": {
                "parallel": "sign*p0*xi*(0.25 + 0.75*u_parallel)",
                "perpendicular_radius": "sign*p0*sqrt(u_perp/xi)",
                "azimuth": "2*pi*u_phi",
                "pairing": "adjacent sorted particle tags carry antipodal states",
            },
        },
        "immutable_extraction_path": {
            "mhd": "bin/<basename>.mhd_bcc.<cycle>.bin",
            "particles": "pvtk/<basename>.prtcl_all.<cycle>.part.vtk",
            "artifact_integrity": "sha256_bound_per_input_artifact",
            "diagnostic_scope":
                "cycle_local_seeded_wave_spectrum_and_particle_pair_payload_only",
        },
    }


def _finite_array(label: str, value: Any) -> np.ndarray:
    array = np.asarray(value)
    if not np.all(np.isfinite(array)):
        raise ContractError(f"{label} must contain only finite values")
    return array


def analyze_cycle_local_arrays(
    bcc2: Any,
    bcc3: Any,
    particle_tags: Any,
    particle_velocity: Any,
    macro_weight: Any,
    deltaf_f0: Any,
    deltaf_weight: Any,
) -> dict[str, Any]:
    """Extract bounded diagnostics without comparison thresholds or calibration claims."""
    by = _finite_array("bcc2", bcc2).astype(np.float64)
    bz = _finite_array("bcc3", bcc3).astype(np.float64)
    if by.shape != bz.shape or by.ndim < 1 or by.shape[-1] != 32:
        raise ContractError("mhd_bcc payload must retain the thin-carrier x1 extent")
    reduce_axes = tuple(range(by.ndim - 1))
    by_x1 = np.mean(by, axis=reduce_axes) if reduce_axes else by
    bz_x1 = np.mean(bz, axis=reduce_axes) if reduce_axes else bz
    by_fft = np.fft.rfft(by_x1 - np.mean(by_x1)) / by_x1.size
    bz_fft = np.fft.rfft(bz_x1 - np.mean(bz_x1)) / bz_x1.size
    spectrum = []
    for mode in range(1, WAVE_MODE_COUNT + 1):
        positive = 0.5 * (by_fft[mode] - 1.0j * bz_fft[mode])
        negative = 0.5 * (by_fft[mode] + 1.0j * bz_fft[mode])
        spectrum.append({
            "mode": mode,
            "positive_helicity_power": float(abs(positive) ** 2),
            "negative_helicity_power": float(abs(negative) ** 2),
        })

    tags = np.asarray(particle_tags)
    if tags.ndim != 1 or tags.size == 0 or tags.size % 2:
        raise ContractError("particle tag payload must contain non-empty adjacent pairs")
    if not np.issubdtype(tags.dtype, np.integer):
        raise ContractError("particle tag payload must be integral")
    order = np.argsort(tags)
    sorted_tags = tags[order]
    if not np.array_equal(sorted_tags, np.arange(sorted_tags.size)):
        raise ContractError("particle tags must be the exact bounded serial tag sequence")
    velocity = _finite_array("particle velocity", particle_velocity).astype(np.float64)
    if velocity.shape != (tags.size, 3):
        raise ContractError("particle velocity payload shape mismatch")
    velocity = velocity[order]
    pair_residual = velocity[0::2] + velocity[1::2]
    weights = _finite_array("macro_weight", macro_weight).astype(np.float64)[order]
    f0 = _finite_array("deltaf_f0", deltaf_f0).astype(np.float64)[order]
    df = _finite_array("deltaf_weight", deltaf_weight).astype(np.float64)[order]
    for label, array in (("macro_weight", weights), ("deltaf_f0", f0),
                         ("deltaf_weight", df)):
        if array.shape != (tags.size,):
            raise ContractError(f"{label} payload shape mismatch")

    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "physical_calibration_claimed": False,
        "q022_closure_claimed": False,
        "external_review_claimed": False,
        "frontier_qualification_claimed": False,
        "launch_authorization_claimed": False,
        "wave_spectrum": spectrum,
        "particle_payload": {
            "count": int(tags.size),
            "adjacent_antipodal_pair_count": int(tags.size // 2),
            "maximum_velocity_pair_residual": float(np.max(np.abs(pair_residual))),
            "mean_parallel_velocity_squared": float(np.mean(velocity[:, 0] ** 2)),
            "mean_perpendicular_velocity_squared": float(
                np.mean(velocity[:, 1] ** 2 + velocity[:, 2] ** 2)
            ),
            "macro_weight_range": [float(np.min(weights)), float(np.max(weights))],
            "deltaf_f0_range": [float(np.min(f0)), float(np.max(f0))],
            "deltaf_weight_range": [float(np.min(df)), float(np.max(df))],
        },
        "interpretation":
            "bounded_cycle_local_runtime_diagnostics_without_reference_thresholds",
    }


def _read_mhd_bcc(path: Path) -> dict[str, Any]:
    sys.path.insert(0, str(REPO_ROOT / "vis/python"))
    import bin_convert_new as bin_convert  # noqa: PLC0415

    return bin_convert.read_binary_as_athdf(str(path))


def _read_particle_vtk(path: Path) -> Any:
    from tst.publication.pvtk_particles import read_particle_vtk  # noqa: PLC0415

    return read_particle_vtk(path)


def extract_runtime_artifacts(mhd_bcc_path: Path, particle_vtk_path: Path) -> dict[str, Any]:
    """Read and hash one immutable cycle-local artifact pair."""
    mhd = _read_mhd_bcc(mhd_bcc_path)
    particles = _read_particle_vtk(particle_vtk_path)
    for field in ("bcc2", "bcc3"):
        if field not in mhd:
            raise ContractError(f"mhd_bcc artifact missing {field}")
    for field in ("ptag", "macro_weight", "deltaf_f0", "deltaf_weight"):
        if field not in particles.scalars:
            raise ContractError(f"particle VTK artifact missing {field}")
    if "vel" not in particles.vectors:
        raise ContractError("particle VTK artifact missing vel")
    report = analyze_cycle_local_arrays(
        mhd["bcc2"],
        mhd["bcc3"],
        particles.scalars["ptag"],
        particles.vectors["vel"],
        particles.scalars["macro_weight"],
        particles.scalars["deltaf_f0"],
        particles.scalars["deltaf_weight"],
    )
    report["immutable_artifacts"] = {
        "mhd_bcc": {"path": str(mhd_bcc_path), "sha256": _sha256(mhd_bcc_path)},
        "prtcl_all": {"path": str(particle_vtk_path), "sha256": _sha256(particle_vtk_path)},
    }
    report["mhd_time"] = float(mhd["Time"])
    return report


def _write_json(payload: dict[str, Any], output: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if output is None:
        print(text, end="")
    else:
        output.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mhd-bcc", type=Path)
    parser.add_argument("--particles", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if (args.mhd_bcc is None) != (args.particles is None):
        parser.error("--mhd-bcc and --particles must be provided together")
    payload = (
        build_static_descriptor()
        if args.mhd_bcc is None
        else extract_runtime_artifacts(args.mhd_bcc, args.particles)
    )
    _write_json(payload, args.output)


if __name__ == "__main__":
    main()
