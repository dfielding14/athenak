#!/usr/bin/env python3
"""Strict source-local Q-043 volume-aware Bell current preparation analyzer."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ID = "Q043-BELL-CURRENT-VOLUME-AWARE"
SUPERSEDES_CAMPAIGN_ID = "Q023-PAPER-BELL-LINEAR"
PGEN_NAME = "q043_bell_current_volume_aware"
PGEN_BLOCK = PGEN_NAME
READINESS_RECORD_PREFIX = "q043_bell_current_volume_aware_successor"
SUPERSESSION_GATE_ID = "Q043-BELL-CURRENT-NORMALIZATION-SUPERSESSION"
SUPERSESSION_IDENTITY_ROLE = "authoritative_foundational_current_lineage"
SUPERSESSION_CONTRACT_PATH = (
    "tst/publication/readiness/"
    "q043_bell_current_normalization_supersession_2026-06-06.json"
)
INTEGRATED_RAW_CURRENT_ORACLE_CAMPAIGN_ID = (
    "Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE"
)
INTEGRATED_RAW_CURRENT_ORACLE_IDENTITY = PGEN_NAME
INTEGRATED_RAW_CURRENT_ORACLE_BINDING_STATUS = (
    "source_local_ready_runtime_observation_pending"
)
CURRENT_NORMALIZATION = "deposited_j_over_c_equals_2_b_g_k0"
DEPOSITION_MEASURE = (
    "ppc_times_deposit_qscale_times_species_charge_times_v_cr_over_root_cell_volume"
)
ROOT_CELL_VOLUME_DEFINITION = "global_root_mesh_extents_over_global_root_mesh_counts"
DEPOSIT_QSCALE_SEMANTICS = "root_cell_macro_charge_volume_aware"
FIXED_QSCALE_DISPOSITION = "forbidden_across_dimensions_and_resolutions"
HISTORICAL_INVALID_NORMALIZATION = (
    "omits_root_cell_volume_and_multiplies_by_artificial_C_invalid"
)
SUPERSESSION_EFFECT = (
    "historical_c_multiplied_normalization_invalid_for_qualification"
)
QUALIFICATION_EFFECT = "none_source_local_preparation_only"
LAUNCH_STATUS = "source_local_preparation_only_not_authorized"
OBSERVATION_SCHEMA_VERSION = 1
DIMENSIONS = (1, 2, 3)
OBSERVATION_RELATIVE_TOLERANCE = 5.0e-6
OBSERVATION_ABSOLUTE_TOLERANCE = 1.0e-6

DECKS = {
    dimension: REPO_ROOT
    / (
        "inputs/tests/pic_q043_bell_current_volume_aware_"
        f"{dimension}d_candidate_vl2_tsc.athinput"
    )
    for dimension in DIMENSIONS
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
    ("particles", "ppc"): "1.0",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "nspecies"): "1",
    ("particles", "cr_distribution"): "center",
    ("particles", "deposit_moments"): "true",
    ("particles", "deposit_order"): "2",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "couple_j_to_efield_coeff"): "1.0",
    ("particles", "couple_j_to_efield_representation"): "cell_centered",
    ("particles", "couple_j_deposition_mode"): "cc_convert",
    ("particles", "couple_moments_momentum_to_mhd"): "true",
    ("particles", "couple_moments_energy_to_mhd"): "true",
    ("particles", "couple_fluid_feedback_order"): "mhd_src_terms",
    ("particles", "pic_physical_mode"): "paper_mhd_pic_vl2_tsc",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_interp_scheme"): "tsc",
    ("particles", "pic_enable_2d3v"): "true",
    ("particles", "pic_cr_light_speed"): "2500.0",
    ("particles", "pic_cr_initial_state"): "velocity",
    ("particles", "pic_cr_hall_mode"): "off",
    ("particles", "pic_wave_damping_mode"): "off",
    ("species0", "mass"): "1.0",
    ("species0", "charge"): "6.283185307179586e-6",
    ("problem", "pgen_name"): PGEN_NAME,
    (PGEN_BLOCK, "campaign_id"): CAMPAIGN_ID,
    (PGEN_BLOCK, "supersedes_campaign_id"): SUPERSEDES_CAMPAIGN_ID,
    (PGEN_BLOCK, "source_mode"): "corrected_linear_eigenmode",
    (PGEN_BLOCK, "epsilon_default"): "0.4",
    (PGEN_BLOCK, "epsilon"): "0.4",
    (PGEN_BLOCK, "epsilon_grid"): "0.1,0.2,0.4,0.6,0.8",
    (PGEN_BLOCK, "rho"): "1.0",
    (PGEN_BLOCK, "pressure"): "1.0",
    (PGEN_BLOCK, "amplitude"): "1.0e-6",
    (PGEN_BLOCK, "b_g"): "1.0",
    (PGEN_BLOCK, "u_a"): "1.0",
    (PGEN_BLOCK, "wavelength"): "1.0",
    (PGEN_BLOCK, "k0"): "6.283185307179586",
    (PGEN_BLOCK, "omega"): "6.283185307179586e-6",
    (PGEN_BLOCK, "c_over_v_cr"): "1000.0",
    (PGEN_BLOCK, "initial_eigenmode"): "section52_right_polarized_eigenmode",
    (PGEN_BLOCK, "current_normalization"): CURRENT_NORMALIZATION,
    (PGEN_BLOCK, "deposition_measure"): DEPOSITION_MEASURE,
    (PGEN_BLOCK, "root_cell_volume"): ROOT_CELL_VOLUME_DEFINITION,
    (PGEN_BLOCK, "deposit_qscale_semantics"): DEPOSIT_QSCALE_SEMANTICS,
    (PGEN_BLOCK, "fixed_qscale"): FIXED_QSCALE_DISPOSITION,
    (
        PGEN_BLOCK,
        "q_over_mc_representation",
    ): "species_charge_equals_q_over_mc_only_with_explicit_species_mass_equals_1",
    (PGEN_BLOCK, "supersession_effect"): SUPERSESSION_EFFECT,
    (PGEN_BLOCK, "qualification_effect"): QUALIFICATION_EFFECT,
    (PGEN_BLOCK, "launch_authorized"): "false",
    (PGEN_BLOCK, "timestep"): "open_clean_candidate_timestep_freeze",
    ("output1", "file_type"): "bin",
    ("output1", "variable"): "mhd_w_bcc",
    ("output1", "id"): "mhd_w_bcc",
    ("output1", "dcycle"): "1",
    ("output1", "ghost_zones"): "false",
    ("output2", "file_type"): "bin",
    ("output2", "variable"): "prtcl_rho",
    ("output2", "id"): "prtcl_rho",
    ("output2", "dcycle"): "1",
    ("output2", "ghost_zones"): "false",
    ("output3", "file_type"): "bin",
    ("output3", "variable"): "prtcl_jx",
    ("output3", "id"): "prtcl_jx",
    ("output3", "dcycle"): "1",
    ("output3", "ghost_zones"): "false",
    ("output4", "file_type"): "bin",
    ("output4", "variable"): "prtcl_jy",
    ("output4", "id"): "prtcl_jy",
    ("output4", "dcycle"): "1",
    ("output4", "ghost_zones"): "false",
    ("output5", "file_type"): "bin",
    ("output5", "variable"): "prtcl_jz",
    ("output5", "id"): "prtcl_jz",
    ("output5", "dcycle"): "1",
    ("output5", "ghost_zones"): "false",
    ("output6", "file_type"): "rst",
    ("output6", "id"): "rst",
    ("output6", "dcycle"): "1",
    ("output6", "single_file_per_rank"): "false",
}

_OBSERVATION_BUNDLE_KEYS = {
    "schema_version",
    "campaign_id",
    "current_normalization",
    "records",
}
_OBSERVATION_RECORD_KEYS = {
    "dimension",
    "artificial_light_speed",
    "volume_mean_deposited_j_over_c",
}


class ContractError(ValueError):
    """Raised when the corrected Q-043 normalization contract fails closed."""


def parse_athinput(path: Path) -> dict[str, dict[str, str]]:
    """Parse the strict Athena input subset used by the successor decks."""
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


def _finite_float(label: str, value: Any) -> float:
    if isinstance(value, bool):
        raise ContractError(f"{label}: expected a finite number")
    try:
        measured = float(value)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{label}: expected a finite number") from exc
    if not math.isfinite(measured):
        raise ContractError(f"{label}: expected a finite number")
    return measured


def _integer(label: str, value: Any) -> int:
    if isinstance(value, bool):
        raise ContractError(f"{label}: expected an integer")
    try:
        measured = int(value)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{label}: expected an integer") from exc
    if str(measured) != str(value):
        raise ContractError(f"{label}: expected an integer")
    return measured


def _require_close(
    label: str,
    measured: float,
    expected: float,
    *,
    relative_tolerance: float = 1.0e-13,
    absolute_tolerance: float = 1.0e-13,
) -> None:
    if not math.isclose(
        measured,
        expected,
        rel_tol=relative_tolerance,
        abs_tol=absolute_tolerance,
    ):
        raise ContractError(f"{label}: expected {expected!r}, measured {measured!r}")


def _dimension(value: Any) -> int:
    if type(value) is not int or value not in DIMENSIONS:
        raise ContractError("dimension must be one of the source-local integer values")
    return value


def _mode_basis(dimension: int) -> tuple[float, float, float]:
    dimension = _dimension(dimension)
    raw = (1.0, 2.0 if dimension >= 2 else 0.0, 4.0 if dimension >= 3 else 0.0)
    norm = math.sqrt(sum(value * value for value in raw))
    return tuple(value / norm for value in raw)


def _vector_magnitude(vector: Sequence[float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def root_cell_volume(
    bounds: Sequence[Sequence[float]], root_nx: Sequence[int]
) -> float:
    """Return the Cartesian root-cell volume used by the particle loader."""
    if len(bounds) != 3 or len(root_nx) != 3:
        raise ContractError("root-cell geometry must contain exactly three axes")
    volume = 1.0
    for axis, (axis_bounds, count) in enumerate(zip(bounds, root_nx), 1):
        if len(axis_bounds) != 2:
            raise ContractError(f"root-cell axis {axis} bounds are malformed")
        lower = _finite_float(f"root-cell x{axis}min", axis_bounds[0])
        upper = _finite_float(f"root-cell x{axis}max", axis_bounds[1])
        if type(count) is not int or count <= 0 or upper <= lower:
            raise ContractError(f"root-cell axis {axis} geometry is invalid")
        volume *= (upper - lower) / count
    if not math.isfinite(volume) or volume <= 0.0:
        raise ContractError("root-cell volume must be finite and positive")
    return volume


def required_deposit_qscale(
    ppc: float,
    species_charge: float,
    stream_speed: float,
    root_volume: float,
    b_g: float,
    k0: float,
) -> float:
    """Return qscale required by the volume-aware current-density closure."""
    values = {
        "ppc": ppc,
        "species_charge": species_charge,
        "stream_speed": stream_speed,
        "root_cell_volume": root_volume,
        "b_g": b_g,
        "k0": k0,
    }
    parsed = {name: _finite_float(name, value) for name, value in values.items()}
    if any(value <= 0.0 for value in parsed.values()):
        raise ContractError("volume-aware qscale factors must be positive")
    return (
        2.0 * parsed["b_g"] * parsed["k0"] * parsed["root_cell_volume"]
        / (
            parsed["ppc"]
            * parsed["species_charge"]
            * parsed["stream_speed"]
        )
    )


def deposited_j_over_c_vector(
    ppc: float,
    deposit_qscale: float,
    species_charge: float,
    stream_velocity: Sequence[float],
    root_cell_volume: float,
) -> tuple[float, float, float]:
    """Return deposited J_CR/c using the decomposition-invariant root cell volume."""
    ppc = _finite_float("ppc", ppc)
    deposit_qscale = _finite_float("deposit_qscale", deposit_qscale)
    species_charge = _finite_float("species_charge", species_charge)
    root_cell_volume = _finite_float("root_cell_volume", root_cell_volume)
    if (
        ppc <= 0.0
        or deposit_qscale <= 0.0
        or species_charge <= 0.0
        or root_cell_volume <= 0.0
    ):
        raise ContractError("deposited-current factors must be positive")
    if len(stream_velocity) != 3:
        raise ContractError("stream velocity must contain exactly three components")
    velocity = tuple(
        _finite_float(f"stream_velocity[{axis}]", value)
        for axis, value in enumerate(stream_velocity)
    )
    scale = ppc * deposit_qscale * species_charge / root_cell_volume
    return tuple(scale * value for value in velocity)


def _historical_disposition() -> dict[str, Any]:
    return {
        "campaign_id": SUPERSEDES_CAMPAIGN_ID,
        "normalization": HISTORICAL_INVALID_NORMALIZATION,
        "qualification_eligible": False,
        "disposition": "invalid_current_normalization_must_not_be_used_for_qualification",
        "reason": (
            "AthenaK moment deposition multiplies macro weight by configured "
            "species_charge and divides by cell volume. Omitting "
            "V_root_cell and multiplying by artificial C changes the Bell mode."
        ),
        "fixed_qscale_disposition": FIXED_QSCALE_DISPOSITION,
    }


def _integration_binding() -> dict[str, Any]:
    return {
        "selected_successor_campaign_id": CAMPAIGN_ID,
        "selected_generator_name": PGEN_NAME,
        "selected_parameter_block": PGEN_BLOCK,
        "selected_readiness_record_prefix": READINESS_RECORD_PREFIX,
        "supersession_gate_id": SUPERSESSION_GATE_ID,
        "supersession_identity_role": SUPERSESSION_IDENTITY_ROLE,
        "supersession_contract_path": SUPERSESSION_CONTRACT_PATH,
        "supersession_contract_binding_status": (
            "authoritative_q043_foundational_lineage_bound"
        ),
        "integrated_raw_current_oracle_campaign_id": (
            INTEGRATED_RAW_CURRENT_ORACLE_CAMPAIGN_ID
        ),
        "integrated_raw_current_oracle_identity": (
            INTEGRATED_RAW_CURRENT_ORACLE_IDENTITY
        ),
        "integrated_raw_current_oracle_binding_status": (
            INTEGRATED_RAW_CURRENT_ORACLE_BINDING_STATUS
        ),
        "integrated_raw_current_oracle_identity_compatible": True,
        "integrated_raw_current_oracle_required": True,
        "integrated_raw_current_oracle_runtime_observation_complete": False,
        "integrated_raw_current_oracle_accepted_as_qualification_evidence": False,
        "integrated_raw_current_oracle_source_local_ready": True,
        "duplicate_q043_generator_or_oracle_required": False,
        "q023_corrected_linear_role": "downstream_non_authorizing_chronology",
    }


def validate_candidate_deck(path: Path, expected_dimension: int) -> dict[str, Any]:
    """Validate one corrected non-authorizing Section 5.2 preparation deck."""
    expected_dimension = _dimension(expected_dimension)
    blocks = parse_athinput(path)
    for (block, name), expected in _EXPECTED_DECK_VALUES.items():
        measured = blocks.get(block, {}).get(name)
        if measured != expected:
            raise ContractError(
                f"{path}: {block}/{name}: expected {expected!r}, measured {measured!r}"
            )

    metadata = blocks[PGEN_BLOCK]
    if metadata.get("dimension") != str(expected_dimension):
        raise ContractError(f"{path}: unexpected corrected-preparation dimension")
    expected_role = (
        "source_local_thin_2d3v_corrected_preparation_not_authorized"
        if expected_dimension == 1
        else "source_local_corrected_preparation_not_authorized"
    )
    if metadata.get("deck_role") != expected_role:
        raise ContractError(f"{path}: unexpected corrected-preparation role")

    geometry = _EXPECTED_GEOMETRY[expected_dimension]
    nx = tuple(
        _integer(f"{path}: mesh/nx{axis}", blocks["mesh"][f"nx{axis}"])
        for axis in (1, 2, 3)
    )
    meshblock_nx = tuple(
        _integer(
            f"{path}: meshblock/nx{axis}",
            blocks["meshblock"][f"nx{axis}"],
        )
        for axis in (1, 2, 3)
    )
    if nx != geometry["nx"]:
        raise ContractError(f"{path}: mesh cell-count contract mismatch")
    if meshblock_nx != geometry["meshblock_nx"]:
        raise ContractError(f"{path}: meshblock cell-count contract mismatch")
    bounds = []
    for axis, expected_bounds in enumerate(geometry["bounds"], 1):
        measured_bounds = (
            _finite_float(f"{path}: mesh/x{axis}min", blocks["mesh"][f"x{axis}min"]),
            _finite_float(f"{path}: mesh/x{axis}max", blocks["mesh"][f"x{axis}max"]),
        )
        for label, measured, expected in zip(
            ("min", "max"), measured_bounds, expected_bounds
        ):
            _require_close(f"{path}: x{axis}{label}", measured, expected)
        bounds.append(measured_bounds)
    extent = tuple(upper - lower for lower, upper in bounds)
    root_cell_widths = tuple(extent[axis] / nx[axis] for axis in range(3))
    root_volume = root_cell_volume(bounds, nx)
    _require_close(f"{path}: root-cell volume", root_volume, math.prod(root_cell_widths))
    active_dx = root_cell_widths[:expected_dimension]
    for axis, (measured, expected) in enumerate(
        zip(active_dx, geometry["active_dx"]), 1
    ):
        _require_close(f"{path}: dx{axis}", measured, expected)

    particles = blocks["particles"]
    species = blocks["species0"]
    stream_velocity = tuple(
        _finite_float(f"{path}: particles/cr_v{axis}0", particles[f"cr_v{axis}0"])
        for axis in ("x", "y", "z")
    )
    stream_speed = _vector_magnitude(stream_velocity)
    rho = _finite_float(f"{path}: rho", metadata["rho"])
    b_g = _finite_float(f"{path}: b_g", metadata["b_g"])
    u_a = _finite_float(f"{path}: u_a", metadata["u_a"])
    wavelength = _finite_float(f"{path}: wavelength", metadata["wavelength"])
    k0 = _finite_float(f"{path}: k0", metadata["k0"])
    omega = _finite_float(f"{path}: omega", metadata["omega"])
    epsilon = _finite_float(f"{path}: epsilon", metadata["epsilon"])
    c_over_v_cr = _finite_float(f"{path}: c_over_v_cr", metadata["c_over_v_cr"])
    light_speed = _finite_float(
        f"{path}: particles/pic_cr_light_speed",
        particles["pic_cr_light_speed"],
    )
    ppc = _finite_float(f"{path}: particles/ppc", particles["ppc"])
    qscale = _finite_float(
        f"{path}: particles/deposit_qscale", particles["deposit_qscale"]
    )
    mass = _finite_float(f"{path}: species0/mass", species["mass"])
    charge = _finite_float(f"{path}: species0/charge", species["charge"])

    _require_close(f"{path}: U_A", u_a, b_g / math.sqrt(rho))
    _require_close(f"{path}: k0", k0, 2.0 * math.pi / wavelength)
    _require_close(f"{path}: Omega", omega, 1.0e-6 * k0 * u_a)
    _require_close(f"{path}: epsilon", epsilon, u_a / stream_speed)
    _require_close(f"{path}: artificial C", light_speed, c_over_v_cr * stream_speed)
    _require_close(f"{path}: species mass required by q/(mc) shorthand", mass, 1.0)
    _require_close(
        f"{path}: species charge under explicit mass-one q/(mc) shorthand",
        charge,
        omega / b_g,
    )
    for axis, (measured, expected) in enumerate(
        zip(
            stream_velocity,
            (stream_speed * value for value in _mode_basis(expected_dimension)),
        ),
        1,
    ):
        _require_close(f"{path}: diagonal stream x{axis}", measured, expected)

    expected_deposited_j_over_c = 2.0 * b_g * k0
    required_qscale = required_deposit_qscale(
        ppc, charge, stream_speed, root_volume, b_g, k0
    )
    _require_close(
        f"{path}: root-cell-volume-aware deposit_qscale",
        qscale,
        required_qscale,
    )
    current_vector = deposited_j_over_c_vector(
        ppc, qscale, charge, stream_velocity, root_volume
    )
    deposited_j_over_c = _vector_magnitude(current_vector)
    unit_root_volume_qscale = required_deposit_qscale(
        ppc, charge, stream_speed, 1.0, b_g, k0
    )
    unit_root_volume_qscale_value = _vector_magnitude(
        deposited_j_over_c_vector(
            ppc, unit_root_volume_qscale, charge, stream_velocity, root_volume
        )
    )
    historical_invalid_value = (
        expected_deposited_j_over_c * light_speed / root_volume
    )
    _require_close(
        f"{path}: deposited J_CR/c",
        deposited_j_over_c,
        expected_deposited_j_over_c,
    )
    if math.isclose(
        deposited_j_over_c,
        historical_invalid_value,
        rel_tol=1.0e-13,
        abs_tol=1.0e-13,
    ):
        raise ContractError(f"{path}: historical C-multiplied normalization is invalid")

    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "dimension": expected_dimension,
        "campaign_id": CAMPAIGN_ID,
        "supersedes_campaign_id": SUPERSEDES_CAMPAIGN_ID,
        "pgen_name": PGEN_NAME,
        "integration_binding": _integration_binding(),
        "pgen_source_implementation_bound_by_this_preparation": True,
        "pgen_registration_bound_by_this_preparation": True,
        "launch_status": LAUNCH_STATUS,
        "launch_authorized": False,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualification_eligible": False,
        "current_normalization": CURRENT_NORMALIZATION,
        "supersession_effect": SUPERSESSION_EFFECT,
        "carrier_semantics": (
            "physical_1d_transverse_invariant_thin_2d3v_nx2_4"
            if expected_dimension == 1
            else "native_mesh_dimension"
        ),
        "global_nx": list(nx),
        "meshblock_nx": list(meshblock_nx),
        "bounds": [list(item) for item in bounds],
        "active_dx": list(active_dx),
        "root_cell_widths": list(root_cell_widths),
        "root_cell_volume": root_volume,
        "ppc": ppc,
        "deposit_qscale": qscale,
        "ppc_times_deposit_qscale": ppc * qscale,
        "required_deposit_qscale": required_qscale,
        "deposit_qscale_per_root_cell_volume": qscale / root_volume,
        "species_mass": mass,
        "species_charge": charge,
        "species_charge_over_mass": charge / mass,
        "stream_velocity": list(stream_velocity),
        "stream_speed": stream_speed,
        "artificial_light_speed": light_speed,
        "deposited_j_over_c_vector": list(current_vector),
        "deposited_j_over_c": deposited_j_over_c,
        "expected_deposited_j_over_c": expected_deposited_j_over_c,
        "deposited_j_over_c_absolute_residual": abs(
            deposited_j_over_c - expected_deposited_j_over_c
        ),
        "unit_root_volume_qscale": unit_root_volume_qscale,
        "unit_root_volume_qscale_deposited_j_over_c": unit_root_volume_qscale_value,
        "unit_root_volume_qscale_to_corrected_current_ratio": (
            unit_root_volume_qscale_value / expected_deposited_j_over_c
        ),
        "historical_invalid_c_multiplied_value": historical_invalid_value,
        "historical_to_corrected_current_ratio": (
            historical_invalid_value / expected_deposited_j_over_c
        ),
        "deposited_current_output_variables": ["prtcl_jx", "prtcl_jy", "prtcl_jz"],
    }


def validate_source_local_candidate_decks() -> list[dict[str, Any]]:
    """Validate all three corrected decks without launching AthenaK."""
    return [
        validate_candidate_deck(DECKS[dimension], dimension)
        for dimension in DIMENSIONS
    ]


def build_preparation_report() -> dict[str, Any]:
    """Report the corrected closure while retaining the non-authorizing boundary."""
    decks = validate_source_local_candidate_decks()
    return {
        "schema_version": 1,
        "record_type": "q043_bell_current_volume_aware_source_local_preparation",
        "campaign_id": CAMPAIGN_ID,
        "supersedes_campaign_id": SUPERSEDES_CAMPAIGN_ID,
        "current_normalization": CURRENT_NORMALIZATION,
        "integration_binding": _integration_binding(),
        "normalization_contract_pass": True,
        "artificial_c_in_formula": False,
        "root_cell_volume_in_formula": True,
        "meshblock_decomposition_in_formula": False,
        "fixed_qscale_valid": False,
        "dimension_and_resolution_aware_qscale_required": True,
        "artificial_c_invariance_runtime_observation_pending": True,
        "historical_campaign_disposition": _historical_disposition(),
        "launch_status": LAUNCH_STATUS,
        "launch_authorized": False,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualification_eligible": False,
        "scientific_claim_authorized": False,
        "passed": False,
        "deck_contracts": decks,
        "status": (
            "corrected_source_local_preparation_consistent_"
            "not_authorized_not_qualifying_evidence"
        ),
    }


def _finite_vector(label: str, value: Any) -> tuple[float, float, float]:
    if not isinstance(value, list) or len(value) != 3:
        raise ContractError(f"{label} must contain exactly three components")
    return tuple(
        _finite_float(f"{label}[{axis}]", component)
        for axis, component in enumerate(value)
    )


def analyze_observed_current_bundle(bundle: Any) -> dict[str, Any]:
    """Validate direct deposited-current observations at multiple artificial C."""
    if not isinstance(bundle, dict) or set(bundle) != _OBSERVATION_BUNDLE_KEYS:
        raise ContractError("observed-current bundle keys do not match")
    if bundle["schema_version"] != OBSERVATION_SCHEMA_VERSION:
        raise ContractError("observed-current schema version mismatch")
    if bundle["campaign_id"] != CAMPAIGN_ID:
        raise ContractError("observed-current campaign id mismatch")
    if bundle["current_normalization"] != CURRENT_NORMALIZATION:
        raise ContractError("observed-current normalization id mismatch")
    if not isinstance(bundle["records"], list):
        raise ContractError("observed-current records must be a list")

    deck_contracts = {
        deck["dimension"]: deck for deck in validate_source_local_candidate_decks()
    }
    records_by_dimension: dict[int, list[dict[str, Any]]] = {
        dimension: [] for dimension in DIMENSIONS
    }
    measured_keys: set[tuple[int, float]] = set()
    for record in bundle["records"]:
        if not isinstance(record, dict) or set(record) != _OBSERVATION_RECORD_KEYS:
            raise ContractError("observed-current record keys do not match")
        dimension = _dimension(record["dimension"])
        artificial_c = _finite_float(
            "observed-current artificial_light_speed",
            record["artificial_light_speed"],
        )
        if artificial_c <= deck_contracts[dimension]["stream_speed"]:
            raise ContractError("observed-current artificial light speed must exceed v_CR")
        key = (dimension, artificial_c)
        if key in measured_keys:
            raise ContractError(f"duplicate observed-current record {key!r}")
        measured_keys.add(key)
        vector = _finite_vector(
            "observed-current volume_mean_deposited_j_over_c",
            record["volume_mean_deposited_j_over_c"],
        )
        expected_vector = deck_contracts[dimension]["deposited_j_over_c_vector"]
        for axis, (measured, expected) in enumerate(zip(vector, expected_vector), 1):
            _require_close(
                f"observed-current dimension {dimension} x{axis}",
                measured,
                expected,
                relative_tolerance=OBSERVATION_RELATIVE_TOLERANCE,
                absolute_tolerance=OBSERVATION_ABSOLUTE_TOLERANCE,
            )
        records_by_dimension[dimension].append(
            {
                "dimension": dimension,
                "artificial_light_speed": artificial_c,
                "volume_mean_deposited_j_over_c": list(vector),
                "deposited_j_over_c": _vector_magnitude(vector),
            }
        )

    invariance = []
    for dimension in DIMENSIONS:
        records = records_by_dimension[dimension]
        if len(records) < 2:
            raise ContractError(
                f"dimension {dimension} requires at least two artificial-C observations"
            )
        records.sort(key=lambda item: item["artificial_light_speed"])
        values = [record["deposited_j_over_c"] for record in records]
        _require_close(
            f"dimension {dimension} observed artificial-C invariance",
            max(values),
            min(values),
            relative_tolerance=OBSERVATION_RELATIVE_TOLERANCE,
            absolute_tolerance=OBSERVATION_ABSOLUTE_TOLERANCE,
        )
        invariance.append(
            {
                "dimension": dimension,
                "artificial_light_speeds": [
                    record["artificial_light_speed"] for record in records
                ],
                "deposited_j_over_c_values": values,
                "expected_deposited_j_over_c": deck_contracts[dimension][
                    "expected_deposited_j_over_c"
                ],
                "maximum_absolute_spread": max(values) - min(values),
                "artificial_c_invariance_pass": True,
            }
        )

    return {
        "schema_version": OBSERVATION_SCHEMA_VERSION,
        "record_type": "q043_bell_current_volume_aware_observed_current_analysis",
        "campaign_id": CAMPAIGN_ID,
        "supersedes_campaign_id": SUPERSEDES_CAMPAIGN_ID,
        "current_normalization": CURRENT_NORMALIZATION,
        "integration_binding": _integration_binding(),
        "normalization_contract_pass": True,
        "observed_deposited_current_contract_pass": True,
        "artificial_c_invariance_pass": True,
        "root_cell_volume_in_formula": True,
        "fixed_qscale_valid": False,
        "historical_campaign_disposition": _historical_disposition(),
        "observation_relative_tolerance": OBSERVATION_RELATIVE_TOLERANCE,
        "observation_absolute_tolerance": OBSERVATION_ABSOLUTE_TOLERANCE,
        "invariance_by_dimension": invariance,
        "launch_status": LAUNCH_STATUS,
        "launch_authorized": False,
        "qualification_effect": QUALIFICATION_EFFECT,
        "qualification_eligible": False,
        "scientific_claim_authorized": False,
        "passed": False,
        "status": (
            "observed_current_bundle_analysis_passed_"
            "source_local_non_authorizing_nonqualifying"
        ),
    }


def synthetic_observed_current_bundle(
    artificial_light_speeds: Iterable[float],
) -> dict[str, Any]:
    """Build a test-only direct-current fixture; never qualifying evidence."""
    light_speeds = [
        _finite_float("synthetic artificial light speed", value)
        for value in artificial_light_speeds
    ]
    if len(set(light_speeds)) < 2:
        raise ContractError("synthetic fixture requires at least two artificial C values")
    records = []
    for deck in validate_source_local_candidate_decks():
        if any(value <= deck["stream_speed"] for value in light_speeds):
            raise ContractError("synthetic artificial light speed must exceed v_CR")
        for light_speed in light_speeds:
            records.append(
                {
                    "dimension": deck["dimension"],
                    "artificial_light_speed": light_speed,
                    "volume_mean_deposited_j_over_c": deck[
                        "deposited_j_over_c_vector"
                    ],
                }
            )
    return {
        "schema_version": OBSERVATION_SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "current_normalization": CURRENT_NORMALIZATION,
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-report", action="store_true")
    parser.add_argument("observed_current_bundle", type=Path, nargs="?")
    args = parser.parse_args()
    if args.build_report:
        if args.observed_current_bundle is not None:
            parser.error("--build-report does not accept an observed-current bundle")
        report = build_preparation_report()
    else:
        if args.observed_current_bundle is None:
            parser.error("analysis requires an observed-current bundle or --build-report")
        report = analyze_observed_current_bundle(
            json.loads(args.observed_current_bundle.read_text(encoding="utf-8"))
        )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
