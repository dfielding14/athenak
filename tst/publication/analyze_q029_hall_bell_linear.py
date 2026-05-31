#!/usr/bin/env python3
"""Q-029 experimental Hall-Bell source-local launch-preparation analyzer."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ID = "Q029-HALL-BELL-LINEAR"
PGEN_NAME = "q029_hall_bell_linear"
ARTIFACT_ROLE = "source_local_hall_bell_launch_preparation_grid_only"
QUALIFICATION_EFFECT = "none_source_local_launch_preparation_only"
LAUNCH_STATUS = "source_local_candidate_only_not_authorized"
EPSILON_VALUES = (0.1, 0.2, 0.4, 0.6, 0.8)
CHI_H_VALUES = (0.25, 0.5, 1.0)
DIMENSIONS = (1, 2, 3)

DECKS = {
    1: REPO_ROOT / "inputs/tests/pic_q029_hall_bell_linear_1d_candidate.athinput",
    2: REPO_ROOT / "inputs/tests/pic_q029_hall_bell_linear_2d_candidate.athinput",
    3: REPO_ROOT / "inputs/tests/pic_q029_hall_bell_linear_3d_candidate.athinput",
}

_EXPECTED_GEOMETRY = {
    1: {
        "nx": (32, 4, 1),
        "meshblock_nx": (32, 4, 1),
        "bounds": ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        "active_dx": (1.0 / 32.0,),
    },
    2: {
        "nx": (64, 32, 1),
        "meshblock_nx": (32, 32, 1),
        "bounds": ((0.0, math.sqrt(5.0)), (0.0, math.sqrt(1.25)), (0.0, 1.0)),
        "active_dx": (math.sqrt(5.0) / 64.0, math.sqrt(1.25) / 32.0),
    },
    3: {
        "nx": (128, 64, 32),
        "meshblock_nx": (32, 32, 32),
        "bounds": (
            (0.0, math.sqrt(21.0)),
            (0.0, math.sqrt(5.25)),
            (0.0, math.sqrt(1.3125)),
        ),
        "active_dx": (
            math.sqrt(21.0) / 128.0,
            math.sqrt(5.25) / 64.0,
            math.sqrt(1.3125) / 32.0,
        ),
    },
}

_EXPECTED_DECK_VALUES = {
    ("time", "nlim"): "0",
    ("time", "tlim"): "0.0",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "ppc"): "2.0",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "nspecies"): "1",
    ("particles", "cr_distribution"): "center",
    ("particles", "deposit_moments"): "true",
    ("particles", "deposit_order"): "1",
    ("particles", "deposit_qscale"): "1.0e9",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "couple_j_to_efield_coeff"): "1.591549430918954e-5",
    ("particles", "couple_j_to_efield_representation"): "cell_centered",
    ("particles", "couple_j_deposition_mode"): "cc_convert",
    ("particles", "couple_moments_momentum_to_mhd"): "true",
    ("particles", "couple_moments_energy_to_mhd"): "true",
    ("particles", "couple_fluid_feedback_order"): "mhd_src_terms",
    ("particles", "pic_physical_mode"): "extended_mhd_pic",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_enable_2d3v"): "true",
    ("particles", "pic_cr_light_speed"): "2500.0",
    ("particles", "pic_cr_initial_state"): "velocity",
    ("particles", "pic_cr_hall_mode"): "current_to_ct_experimental",
    ("particles", "pic_wave_damping_mode"): "off",
    ("species0", "mass"): "1.0",
    ("species0", "charge"): "6.283185307179586e-6",
    ("problem", "pgen_name"): PGEN_NAME,
    ("q029_hall_bell_linear", "campaign_id"): CAMPAIGN_ID,
    ("q029_hall_bell_linear", "epsilon_default"): "0.4",
    ("q029_hall_bell_linear", "epsilon"): "0.4",
    ("q029_hall_bell_linear", "epsilon_grid"): "0.1,0.2,0.4,0.6,0.8",
    ("q029_hall_bell_linear", "chi_h_default"): "0.5",
    ("q029_hall_bell_linear", "chi_h"): "0.5",
    ("q029_hall_bell_linear", "chi_h_grid"): "0.25,0.5,1.0",
    ("q029_hall_bell_linear", "alpha_h_parameter"):
        "particles/couple_j_to_efield_coeff",
    ("q029_hall_bell_linear", "chi_h_definition"): "alpha_h*j_cr/(u_a*b_g)",
    ("q029_hall_bell_linear", "chi_h_grid_semantics"):
        "positive_launch_preparation_only_no_hall_bell_qualification",
    ("q029_hall_bell_linear", "rho"): "1.0",
    ("q029_hall_bell_linear", "pressure"): "1.0",
    ("q029_hall_bell_linear", "amplitude"): "1.0e-6",
    ("q029_hall_bell_linear", "b_g"): "1.0",
    ("q029_hall_bell_linear", "u_a"): "1.0",
    ("q029_hall_bell_linear", "wavelength"): "1.0",
    ("q029_hall_bell_linear", "k0"): "6.283185307179586",
    ("q029_hall_bell_linear", "omega"): "6.283185307179586e-6",
    ("q029_hall_bell_linear", "c_over_v_cr"): "1000.0",
    ("q029_hall_bell_linear", "initial_eigenmode"):
        "q023_section52_seed_carrier_only_not_hall_dispersion_oracle",
    ("q029_hall_bell_linear", "qualification_effect"): "none",
    ("q029_hall_bell_linear", "linear_hall_bell"): "open_not_claimed",
    ("q029_hall_bell_linear", "nonlinear_hall_bell"): "open_not_claimed",
    ("q029_hall_bell_linear", "timestep"):
        "open_clean_candidate_timestep_freeze",
    ("output1", "file_type"): "bin",
    ("output1", "variable"): "mhd_bcc",
    ("output1", "id"): "mhd_bcc",
    ("output1", "dcycle"): "1",
    ("output1", "ghost_zones"): "false",
    ("output2", "file_type"): "rst",
    ("output2", "id"): "rst",
    ("output2", "dcycle"): "1",
    ("output2", "single_file_per_rank"): "false",
}

_BUNDLE_KEYS = {"schema_version", "campaign_id", "artifact_role", "variants"}
_VARIANT_KEYS = {
    "variant_id",
    "deck_path",
    "dimension",
    "epsilon",
    "chi_h",
    "stream_velocity",
    "pic_cr_light_speed",
    "j_cr",
    "alpha_h",
    "athena_overrides",
    "variant_role",
}
_VARIANT_ROLE = "source_local_materialization_pending_not_authorized"


class ContractError(ValueError):
    """Raised when the Q-029 launch-preparation contract fails closed."""


def parse_athinput(path: Path) -> dict[str, dict[str, str]]:
    """Parse the strict Athena input subset used by the source-local decks."""
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


def _require_close(label: str, measured: float, expected: float) -> None:
    if not math.isfinite(measured) or not math.isfinite(expected):
        raise ContractError(f"{label}: values must be finite")
    if not math.isclose(measured, expected, rel_tol=1.0e-14, abs_tol=1.0e-14):
        raise ContractError(f"{label}: expected {expected!r}, measured {measured!r}")


def _parse_finite_float(label: str, value: str) -> float:
    try:
        measured = float(value)
    except ValueError as exc:
        raise ContractError(f"{label}: expected a finite number") from exc
    if not math.isfinite(measured):
        raise ContractError(f"{label}: expected a finite number")
    return measured


def _parse_integer(label: str, value: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise ContractError(f"{label}: expected an integer") from exc


def _csv_floats(value: str) -> tuple[float, ...]:
    return tuple(_parse_finite_float("CSV value", item) for item in value.split(","))


def _mode_basis(dimension: int) -> tuple[float, float, float]:
    if dimension not in DIMENSIONS:
        raise ContractError("dimension must be one of the source-local integer values")
    raw = (1.0, 2.0 if dimension >= 2 else 0.0, 4.0 if dimension >= 3 else 0.0)
    norm = math.sqrt(sum(value * value for value in raw))
    return tuple(value / norm for value in raw)


def _current_density(blocks: dict[str, dict[str, str]]) -> tuple[float, float]:
    particles = blocks["particles"]
    stream = tuple(
        _parse_finite_float(f"particles/cr_v{axis}0", particles[f"cr_v{axis}0"])
        for axis in ("x", "y", "z")
    )
    speed = math.sqrt(sum(value * value for value in stream))
    current = (
        _parse_finite_float("particles/ppc", particles["ppc"])
        * _parse_finite_float("particles/deposit_qscale", particles["deposit_qscale"])
        * _parse_finite_float("species0/charge", blocks["species0"]["charge"])
        * speed
    )
    if not math.isfinite(speed) or not math.isfinite(current):
        raise ContractError("Q-029 stream speed and j_CR must be finite")
    if speed <= 0.0 or current <= 0.0:
        raise ContractError("Q-029 stream speed and j_CR must be positive")
    return speed, current


def validate_candidate_deck(path: Path, expected_dimension: int) -> dict[str, Any]:
    """Validate one non-authorized Hall-Bell launch-preparation deck."""
    blocks = parse_athinput(path)
    for (block, name), expected in _EXPECTED_DECK_VALUES.items():
        measured = blocks.get(block, {}).get(name)
        if measured != expected:
            raise ContractError(
                f"{path}: {block}/{name}: expected {expected!r}, measured {measured!r}"
            )

    metadata = blocks["q029_hall_bell_linear"]
    if metadata.get("dimension") != str(expected_dimension):
        raise ContractError(f"{path}: unexpected launch-preparation dimension")
    expected_role = (
        "source_local_runnable_thin_2d3v_launch_preparation_not_authorized"
        if expected_dimension == 1
        else "source_local_runnable_launch_preparation_not_authorized"
    )
    if metadata.get("deck_role") != expected_role:
        raise ContractError(f"{path}: unexpected source-local launch role")

    chi_h_grid = _csv_floats(metadata["chi_h_grid"])
    if chi_h_grid != CHI_H_VALUES or not all(value > 0.0 for value in chi_h_grid):
        raise ContractError(f"{path}: chi_H grid must contain only prepared positives")

    geometry = _EXPECTED_GEOMETRY[expected_dimension]
    nx = tuple(
        _parse_integer(f"{path}: mesh/nx{axis}", blocks["mesh"][f"nx{axis}"])
        for axis in (1, 2, 3)
    )
    meshblock_nx = tuple(
        _parse_integer(
            f"{path}: meshblock/nx{axis}", blocks["meshblock"][f"nx{axis}"]
        )
        for axis in (1, 2, 3)
    )
    if nx != geometry["nx"]:
        raise ContractError(f"{path}: global mesh cell-count contract mismatch")
    if meshblock_nx != geometry["meshblock_nx"]:
        raise ContractError(f"{path}: meshblock cell-count contract mismatch")
    bounds = []
    for axis, expected_bounds in enumerate(geometry["bounds"], 1):
        measured_bounds = (
            _parse_finite_float(
                f"{path}: mesh/x{axis}min", blocks["mesh"][f"x{axis}min"]
            ),
            _parse_finite_float(
                f"{path}: mesh/x{axis}max", blocks["mesh"][f"x{axis}max"]
            ),
        )
        for label, measured, expected in zip(
            ("min", "max"), measured_bounds, expected_bounds
        ):
            _require_close(f"{path}: x{axis}{label}", measured, expected)
        bounds.append(measured_bounds)
    extent = tuple(upper - lower for lower, upper in bounds)
    active_dx = tuple(extent[axis] / nx[axis] for axis in range(expected_dimension))
    for axis, (measured, expected) in enumerate(
        zip(active_dx, geometry["active_dx"]), 1
    ):
        _require_close(f"{path}: dx{axis}", measured, expected)

    rho = _parse_finite_float(f"{path}: rho", metadata["rho"])
    b_g = _parse_finite_float(f"{path}: b_g", metadata["b_g"])
    u_a = _parse_finite_float(f"{path}: u_a", metadata["u_a"])
    wavelength = _parse_finite_float(f"{path}: wavelength", metadata["wavelength"])
    k0 = _parse_finite_float(f"{path}: k0", metadata["k0"])
    omega = _parse_finite_float(f"{path}: omega", metadata["omega"])
    epsilon = _parse_finite_float(f"{path}: epsilon", metadata["epsilon"])
    chi_h = _parse_finite_float(f"{path}: chi_h", metadata["chi_h"])
    c_over_v_cr = _parse_finite_float(
        f"{path}: c_over_v_cr", metadata["c_over_v_cr"]
    )
    alpha_h = _parse_finite_float(
        f"{path}: particles/couple_j_to_efield_coeff",
        blocks["particles"]["couple_j_to_efield_coeff"],
    )
    speed, j_cr = _current_density(blocks)
    stream = tuple(
        _parse_finite_float(
            f"{path}: particles/cr_v{axis}0",
            blocks["particles"][f"cr_v{axis}0"],
        )
        for axis in ("x", "y", "z")
    )
    light_speed = _parse_finite_float(
        f"{path}: particles/pic_cr_light_speed",
        blocks["particles"]["pic_cr_light_speed"],
    )
    ppc = _parse_finite_float(f"{path}: particles/ppc", blocks["particles"]["ppc"])
    qscale = _parse_finite_float(
        f"{path}: particles/deposit_qscale", blocks["particles"]["deposit_qscale"]
    )
    mass = _parse_finite_float(f"{path}: species0/mass", blocks["species0"]["mass"])
    charge = _parse_finite_float(
        f"{path}: species0/charge", blocks["species0"]["charge"]
    )

    if chi_h not in CHI_H_VALUES or chi_h <= 0.0:
        raise ContractError(f"{path}: fiducial chi_H is not a positive prepared value")
    _require_close(f"{path}: U_A", u_a, b_g / math.sqrt(rho))
    _require_close(f"{path}: k0", k0, 2.0 * math.pi / wavelength)
    _require_close(f"{path}: Omega", omega, 1.0e-6 * k0 * u_a)
    _require_close(f"{path}: epsilon", epsilon, u_a / speed)
    _require_close(f"{path}: C", light_speed, c_over_v_cr * speed)
    _require_close(f"{path}: species q/m", charge / mass, omega / b_g)
    _require_close(f"{path}: j_CR", j_cr, 2.0 * b_g * light_speed * k0)
    _require_close(f"{path}: chi_H", chi_h, alpha_h * j_cr / (u_a * b_g))
    for axis, (measured, expected) in enumerate(
        zip(stream, (speed * value for value in _mode_basis(expected_dimension))), 1
    ):
        _require_close(f"{path}: diagonal seed-carrier stream x{axis}", measured, expected)

    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "dimension": expected_dimension,
        "launch_status": LAUNCH_STATUS,
        "qualification_effect": "none",
        "qualifying_evidence": False,
        "hall_bell_qualification": False,
        "carrier_semantics": (
            "physical_1d_transverse_invariant_thin_2d3v_nx2_4"
            if expected_dimension == 1
            else "native_mesh_dimension"
        ),
        "global_nx": list(nx),
        "meshblock_nx": list(meshblock_nx),
        "bounds": [list(axis_bounds) for axis_bounds in bounds],
        "active_dx": list(active_dx),
        "pic_enable_2d3v": blocks["particles"]["pic_enable_2d3v"] == "true",
        "epsilon_default": _parse_finite_float(
            f"{path}: epsilon_default", metadata["epsilon_default"]
        ),
        "chi_h_default": _parse_finite_float(
            f"{path}: chi_h_default", metadata["chi_h_default"]
        ),
        "chi_h_grid": list(chi_h_grid),
        "u_a": u_a,
        "b_g": b_g,
        "c_over_v_cr": c_over_v_cr,
        "ppc": ppc,
        "deposit_qscale": qscale,
        "charge": charge,
        "stream_speed": speed,
        "pic_cr_light_speed": light_speed,
        "alpha_h_fiducial": alpha_h,
        "j_cr": j_cr,
    }


def validate_source_local_candidate_decks() -> list[dict[str, Any]]:
    """Validate all three Q-029 decks without launching AthenaK."""
    return [
        validate_candidate_deck(DECKS[dimension], dimension)
        for dimension in DIMENSIONS
    ]


def validate_generator_registration() -> dict[str, Any]:
    """Require the separate built-in identity and guarded Q-023 reuse boundary."""
    dispatch = (REPO_ROOT / "src/pgen/pgen.cpp").read_text(encoding="utf-8")
    declarations = (REPO_ROOT / "src/pgen/pgen.hpp").read_text(encoding="utf-8")
    cmake = (REPO_ROOT / "src/CMakeLists.txt").read_text(encoding="utf-8")
    generator = (
        REPO_ROOT / "src/pgen/tests/q029_hall_bell_linear.cpp"
    ).read_text(encoding="utf-8")
    q023_generator = (
        REPO_ROOT / "src/pgen/tests/q023_paper_bell_linear.cpp"
    ).read_text(encoding="utf-8")
    carrier = (
        REPO_ROOT / "src/pgen/tests/q023_paper_bell_linear.hpp"
    ).read_text(encoding="utf-8")
    if dispatch.count(f'compare("{PGEN_NAME}")') != 2:
        raise ContractError("Q-029 dispatch must cover fresh and restart constructors")
    if dispatch.count("Q029HallBellLinear(pin, false);") != 1:
        raise ContractError("Q-029 fresh constructor dispatch is missing")
    if dispatch.count("Q029HallBellLinear(pin, true);") != 1:
        raise ContractError("Q-029 restart constructor dispatch is missing")
    if "void Q029HallBellLinear(ParameterInput *pin, const bool restart);" not in declarations:
        raise ContractError("Q-029 generator declaration is missing")
    if "pgen/tests/q029_hall_bell_linear.cpp" not in cmake:
        raise ContractError("Q-029 generator compilation unit is missing")
    required_source_snippets = (
        '#include "q023_paper_bell_linear.hpp"',
        '"extended_mhd_pic"',
        '"current_to_ct_experimental"',
        '"positive_launch_preparation_only_no_hall_bell_qualification"',
        '"open_not_claimed"',
        'Q029RequireBoolean(pin, "particles", "pic_enable_2d3v", true);',
        "std::isfinite",
    )
    if any(snippet not in generator for snippet in required_source_snippets):
        raise ContractError("Q-029 generator identity or nonqualification boundary is absent")
    if '#include "q023_paper_bell_linear.cpp"' in generator:
        raise ContractError("Q-029 must not textually include the Q-023 compilation unit")
    if '#include "q023_paper_bell_linear.hpp"' not in q023_generator:
        raise ContractError("Q-023 must consume its guarded seed-carrier header")
    required_carrier_snippets = (
        "#ifndef PGEN_TESTS_Q023_PAPER_BELL_LINEAR_HPP_",
        "#define Q023_INLINE KOKKOS_INLINE_FUNCTION",
        "#define Q023_INLINE inline",
        "EigenmodeAtPhase",
        "VectorPotentialAt",
    )
    if any(snippet not in carrier for snippet in required_carrier_snippets):
        raise ContractError("guarded Q-023 accelerator and host carrier contract is absent")
    return {
        "pgen_name": PGEN_NAME,
        "fresh_dispatch": True,
        "restart_dispatch": True,
        "q023_seed_carrier_reused_from_guarded_header": True,
    }


def _variant_token(value: float) -> str:
    return repr(float(value)).replace(".", "p")


def _variant_id(dimension: int, epsilon: float, chi_h: float) -> str:
    return (
        f"Q029-SOURCE-LOCAL-PREPARATION-{dimension}D-"
        f"EPSILON-{_variant_token(epsilon)}-CHI-H-{_variant_token(chi_h)}"
    )


def _basename(dimension: int, epsilon: float, chi_h: float) -> str:
    return (
        f"pic_q029_hall_bell_linear_{dimension}d_"
        f"epsilon_{_variant_token(epsilon)}_chi_h_{_variant_token(chi_h)}_candidate"
    )


def _athena_override(name: str, value: Any) -> str:
    rendered = value if isinstance(value, str) else repr(float(value))
    return f"{name}={rendered}"


def _materialized_variant(
    dimension: int,
    epsilon: float,
    chi_h: float,
    deck: dict[str, Any],
) -> dict[str, Any]:
    stream_speed = deck["u_a"] / epsilon
    stream_velocity = [
        stream_speed * component for component in _mode_basis(dimension)
    ]
    light_speed = deck["c_over_v_cr"] * stream_speed
    j_cr = deck["ppc"] * deck["deposit_qscale"] * deck["charge"] * stream_speed
    alpha_h = chi_h * deck["u_a"] * deck["b_g"] / j_cr
    basename = _basename(dimension, epsilon, chi_h)
    overrides = [
        _athena_override("job/basename", basename),
        _athena_override("q029_hall_bell_linear/epsilon", epsilon),
        _athena_override("q029_hall_bell_linear/chi_h", chi_h),
        _athena_override("particles/cr_vx0", stream_velocity[0]),
        _athena_override("particles/cr_vy0", stream_velocity[1]),
        _athena_override("particles/cr_vz0", stream_velocity[2]),
        _athena_override("particles/pic_cr_light_speed", light_speed),
        _athena_override("particles/couple_j_to_efield_coeff", alpha_h),
    ]
    return {
        "variant_id": _variant_id(dimension, epsilon, chi_h),
        "deck_path": deck["path"],
        "dimension": dimension,
        "epsilon": epsilon,
        "chi_h": chi_h,
        "stream_velocity": stream_velocity,
        "pic_cr_light_speed": light_speed,
        "j_cr": j_cr,
        "alpha_h": alpha_h,
        "athena_overrides": overrides,
        "variant_role": _VARIANT_ROLE,
    }


def build_launch_preparation_bundle() -> dict[str, Any]:
    """Build the exact positive chi_H variant-materialization grid."""
    decks = {
        deck["dimension"]: deck for deck in validate_source_local_candidate_decks()
    }
    variants = []
    for dimension in DIMENSIONS:
        for epsilon in EPSILON_VALUES:
            for chi_h in CHI_H_VALUES:
                variants.append(
                    _materialized_variant(dimension, epsilon, chi_h, decks[dimension])
                )
    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "variants": variants,
    }


def _finite_value(label: str, value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"{label} must be a finite number")
    measured = float(value)
    if not math.isfinite(measured):
        raise ContractError(f"{label} must be a finite number")
    return measured


def _finite_number(variant: dict[str, Any], key: str) -> float:
    return _finite_value(key, variant[key])


def analyze_launch_preparation_bundle(bundle: dict[str, Any]) -> dict[str, Any]:
    """Analyze only variant preparation and never emit Hall-Bell qualification."""
    if not isinstance(bundle, dict) or set(bundle) != _BUNDLE_KEYS:
        raise ContractError("Q-029 launch-preparation bundle keys do not match")
    if bundle["schema_version"] != 1 or bundle["campaign_id"] != CAMPAIGN_ID:
        raise ContractError("Q-029 launch-preparation bundle identity mismatch")
    if bundle["artifact_role"] != ARTIFACT_ROLE:
        raise ContractError("Q-029 bundle is not source-local launch preparation")
    if not isinstance(bundle["variants"], list):
        raise ContractError("Q-029 launch-preparation variants must be a list")

    decks = validate_source_local_candidate_decks()
    registration = validate_generator_registration()
    expected_keys = {
        (dimension, epsilon, chi_h)
        for dimension in DIMENSIONS
        for epsilon in EPSILON_VALUES
        for chi_h in CHI_H_VALUES
    }
    measured_keys = []
    variants = []
    decks_by_dimension = {deck["dimension"]: deck for deck in decks}
    for variant in bundle["variants"]:
        if not isinstance(variant, dict) or set(variant) != _VARIANT_KEYS:
            raise ContractError("Q-029 launch-preparation variant keys do not match")
        dimension = variant["dimension"]
        if type(dimension) is not int or dimension not in DIMENSIONS:
            raise ContractError("Q-029 launch-preparation dimension is invalid")
        epsilon = _finite_number(variant, "epsilon")
        chi_h = _finite_number(variant, "chi_h")
        stream_velocity = variant["stream_velocity"]
        if not isinstance(stream_velocity, list) or len(stream_velocity) != 3:
            raise ContractError("Q-029 launch-preparation stream velocity is invalid")
        stream_velocity = [
            _finite_value(f"stream_velocity[{axis}]", value)
            for axis, value in enumerate(stream_velocity)
        ]
        light_speed = _finite_number(variant, "pic_cr_light_speed")
        j_cr = _finite_number(variant, "j_cr")
        alpha_h = _finite_number(variant, "alpha_h")
        if epsilon not in EPSILON_VALUES:
            raise ContractError("Q-029 launch-preparation epsilon is outside the grid")
        if chi_h not in CHI_H_VALUES or chi_h <= 0.0:
            raise ContractError("Q-029 launch-preparation chi_H must be prepared positive")
        if alpha_h <= 0.0:
            raise ContractError("Q-029 launch-preparation alpha_H must be positive")
        if variant["variant_role"] != _VARIANT_ROLE:
            raise ContractError("Q-029 launch-preparation variant role mismatch")
        expected = _materialized_variant(
            dimension, epsilon, chi_h, decks_by_dimension[dimension]
        )
        if variant["variant_id"] != expected["variant_id"]:
            raise ContractError("Q-029 launch-preparation variant id mismatch")
        if variant["deck_path"] != expected["deck_path"]:
            raise ContractError("Q-029 launch-preparation deck path mismatch")
        for axis, (measured, expected_value) in enumerate(
            zip(stream_velocity, expected["stream_velocity"])
        ):
            _require_close(
                f"Q-029 stream velocity x{axis + 1} materialization",
                measured,
                expected_value,
            )
        _require_close(
            "Q-029 PIC light-speed materialization",
            light_speed,
            expected["pic_cr_light_speed"],
        )
        _require_close("Q-029 j_CR materialization", j_cr, expected["j_cr"])
        _require_close("Q-029 alpha_H materialization", alpha_h, expected["alpha_h"])
        if variant["athena_overrides"] != expected["athena_overrides"]:
            raise ContractError("Q-029 launch-preparation Athena overrides mismatch")
        key = (dimension, epsilon, chi_h)
        if key in measured_keys:
            raise ContractError(f"duplicate Q-029 launch-preparation variant {key!r}")
        measured_keys.append(key)
        variants.append(dict(variant))
    if set(measured_keys) != expected_keys:
        raise ContractError(
            "Q-029 launch-preparation grid is incomplete or contains extra variants"
        )

    variants.sort(key=lambda item: (item["dimension"], item["epsilon"], item["chi_h"]))
    return {
        "schema_version": 1,
        "campaign_id": CAMPAIGN_ID,
        "artifact_role": ARTIFACT_ROLE,
        "launch_status": LAUNCH_STATUS,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualifying_evidence": False,
        "hall_bell_qualification": False,
        "linear_hall_bell_qualification": False,
        "nonlinear_hall_bell_qualification": False,
        "generator_registration": registration,
        "deck_contracts": decks,
        "variant_count": len(variants),
        "variants": variants,
        "status": "source_local_launch_preparation_consistent_not_qualifying_evidence",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-bundle", action="store_true")
    parser.add_argument("bundle", type=Path, nargs="?")
    args = parser.parse_args()
    if args.build_bundle:
        if args.bundle is not None:
            parser.error("--build-bundle does not accept a bundle path")
        result = build_launch_preparation_bundle()
    else:
        if args.bundle is None:
            parser.error("bundle analysis requires one launch-preparation bundle")
        result = analyze_launch_preparation_bundle(
            json.loads(args.bundle.read_text(encoding="utf-8"))
        )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
