#!/usr/bin/env python3
"""Source-local Q-023 Section 5.2 Bell analytical-dispersion candidate."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
from typing import Any, Callable, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN_ID = "Q023-PAPER-BELL-LINEAR"
TRACE_SCHEMA_VERSION = 2
EPSILON_VALUES = (0.1, 0.2, 0.4, 0.6, 0.8)
DIMENSIONS = (1, 2, 3)
QUALIFYING_SEEDS = (
    23050101,
    23050102,
    23050103,
    23050104,
    23050105,
    23050106,
    23050107,
    23050108,
)
PHASE_INTERVAL = math.pi
MIN_GROWTH_SNAPSHOTS = 8
MIN_RIGHT_TO_LEFT_RATIO = 10.0
ABSOLUTE_TOLERANCE = 0.02
RELATIVE_TOLERANCE = 0.05
SOURCE_LOCAL_VELOCITY_MAGNETIC_RATIO_ABSOLUTE_TOLERANCE = 0.05
RAW_GEOMETRY_ABSOLUTE_TOLERANCE = 1.0e-6
K0 = 2.0 * math.pi
U_A = 1.0
RETAINED_OBSERVABLE_CONTRACT = (
    "combined_mhd_w_bcc_magnetic_and_fluid_velocity_mode_preparation_"
    "section52_qualification_still_blocked"
)

DECKS = {
    1: REPO_ROOT / "inputs/tests/pic_q023_paper_bell_linear_1d_candidate_vl2_tsc.athinput",
    2: REPO_ROOT / "inputs/tests/pic_q023_paper_bell_linear_2d_candidate_vl2_tsc.athinput",
    3: REPO_ROOT / "inputs/tests/pic_q023_paper_bell_linear_3d_candidate_vl2_tsc.athinput",
}

_EXPECTED_GEOMETRY = {
    1: {
        "nx": (32, 4, 1),
        "xmin": (0.0, 0.0, 0.0),
        "extent": (1.0, 1.0, 1.0),
        "active_dx": (1.0 / 32.0,),
    },
    2: {
        "nx": (64, 32, 1),
        "xmin": (0.0, 0.0, 0.0),
        "extent": (math.sqrt(5.0), math.sqrt(1.25), 1.0),
        "active_dx": (math.sqrt(5.0) / 64.0, math.sqrt(1.25) / 32.0),
    },
    3: {
        "nx": (128, 64, 32),
        "xmin": (0.0, 0.0, 0.0),
        "extent": (math.sqrt(21.0), math.sqrt(5.25), math.sqrt(1.3125)),
        "active_dx": (
            math.sqrt(21.0) / 128.0,
            math.sqrt(5.25) / 64.0,
            math.sqrt(1.3125) / 32.0,
        ),
    },
}

_APPROVED_SOURCE_LOCAL_RAW_VARIANTS = {
    "Q023-SOURCE-LOCAL-BASELINE-1D-EPSILON-0P4": {
        "dimension": 1,
        "epsilon": 0.4,
        "deck": DECKS[1],
        "deck_sha256": "b77689de8e683d918cbfda2d13adac77ca7aa85b7f5cb7f998fb1dfc02ce1cba",
    },
    "Q023-SOURCE-LOCAL-BASELINE-2D-EPSILON-0P4": {
        "dimension": 2,
        "epsilon": 0.4,
        "deck": DECKS[2],
        "deck_sha256": "1877445adcecd24dbc91aeda80e5f11a4529fd4ab68049ef000e7cdda79b3524",
    },
    "Q023-SOURCE-LOCAL-BASELINE-3D-EPSILON-0P4": {
        "dimension": 3,
        "epsilon": 0.4,
        "deck": DECKS[3],
        "deck_sha256": "e30e44d1566e4acaf1201c46db8cf37537101c79c2d325a48dbe4ae9288556ad",
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
    ("particles", "deposit_order"): "2",
    ("particles", "deposit_qscale"): "1.0e9",
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
    ("particles", "pic_enable_2d3v"): "true",
    ("particles", "pic_cr_light_speed"): "2500.0",
    ("particles", "pic_cr_initial_state"): "velocity",
    ("particles", "pic_cr_hall_mode"): "off",
    ("particles", "pic_wave_damping_mode"): "off",
    ("species0", "mass"): "1.0",
    ("species0", "charge"): "6.283185307179586e-6",
    ("problem", "pgen_name"): "q023_paper_bell_linear",
    ("q023_paper_bell_linear", "campaign_id"): CAMPAIGN_ID,
    ("q023_paper_bell_linear", "epsilon_default"): "0.4",
    ("q023_paper_bell_linear", "epsilon"): "0.4",
    ("q023_paper_bell_linear", "epsilon_grid"): "0.1,0.2,0.4,0.6,0.8",
    ("q023_paper_bell_linear", "rho"): "1.0",
    ("q023_paper_bell_linear", "pressure"): "1.0",
    ("q023_paper_bell_linear", "amplitude"): "1.0e-6",
    ("q023_paper_bell_linear", "b_g"): "1.0",
    ("q023_paper_bell_linear", "u_a"): "1.0",
    ("q023_paper_bell_linear", "wavelength"): "1.0",
    ("q023_paper_bell_linear", "k0"): "6.283185307179586",
    ("q023_paper_bell_linear", "omega"): "6.283185307179586e-6",
    ("q023_paper_bell_linear", "c_over_v_cr"): "1000.0",
    ("q023_paper_bell_linear", "initial_eigenmode"):
        "section52_right_polarized_eigenmode",
    ("q023_paper_bell_linear", "timestep"):
        "open_clean_candidate_timestep_freeze",
    ("output1", "file_type"): "bin",
    ("output1", "variable"): "mhd_w_bcc",
    ("output1", "id"): "mhd_w_bcc",
    ("output1", "dcycle"): "1",
    ("output1", "ghost_zones"): "false",
    ("output2", "file_type"): "rst",
    ("output2", "id"): "rst",
    ("output2", "dcycle"): "1",
    ("output2", "single_file_per_rank"): "false",
}

_TRACE_KEYS = {
    "dimension",
    "epsilon",
    "raw_provenance",
    "normalized_time",
    "right_mode_real",
    "right_mode_imag",
    "left_mode_real",
    "left_mode_imag",
    "velocity_right_mode_real",
    "velocity_right_mode_imag",
    "velocity_left_mode_real",
    "velocity_left_mode_imag",
    "phase_interval",
    "phase_change",
    "velocity_phase_interval",
    "velocity_phase_change",
    "paper_delta_u_y_sine_fit_real",
    "paper_delta_u_y_sine_fit_imag",
    "paper_delta_u_y_phase_interval",
    "paper_delta_u_y_phase_change",
    "paper_volume_averaged_abs_delta_u",
}
_RAW_PROVENANCE_KEYS = {
    "kind",
    "variant_id",
    "deck_path",
    "deck_sha256",
    "raw_geometry",
    "raw_artifacts",
    "retained_observable_contract",
}
_RAW_GEOMETRY_KEYS = {"nx", "xmin", "extent"}
_RAW_ARTIFACT_KEYS = {"path", "sha256"}
_SHA256 = re.compile(r"[0-9a-f]{64}")


class ContractError(ValueError):
    """Raised when a candidate deck or extracted trace violates the contract."""


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
    if not math.isclose(measured, expected, rel_tol=1.0e-14, abs_tol=1.0e-14):
        raise ContractError(f"{label}: expected {expected!r}, measured {measured!r}")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_dimension(dimension: Any) -> int:
    if type(dimension) is not int or dimension not in DIMENSIONS:
        raise ContractError("dimension must be one of the preregistered integer values")
    return dimension


def _raw_geometry(dimension: int) -> dict[str, list[float] | list[int]]:
    dimension = _require_dimension(dimension)
    geometry = _EXPECTED_GEOMETRY[dimension]
    return {
        "nx": list(geometry["nx"]),
        "xmin": list(geometry["xmin"]),
        "extent": list(geometry["extent"]),
    }


def _validate_raw_geometry(
    geometry: Any, dimension: int
) -> dict[str, list[float] | list[int]]:
    dimension = _require_dimension(dimension)
    if (
        not isinstance(geometry, dict)
        or set(geometry) != _RAW_GEOMETRY_KEYS
        or not isinstance(geometry["nx"], list)
        or len(geometry["nx"]) != 3
        or any(type(value) is not int for value in geometry["nx"])
        or not isinstance(geometry["xmin"], list)
        or len(geometry["xmin"]) != 3
        or any(type(value) is not float for value in geometry["xmin"])
        or not isinstance(geometry["extent"], list)
        or len(geometry["extent"]) != 3
        or any(type(value) is not float for value in geometry["extent"])
        or geometry != _raw_geometry(dimension)
    ):
        raise ContractError(
            "raw Bell provenance geometry is not an approved materialized variant"
        )
    return geometry


def synthetic_contract_provenance(dimension: int, epsilon: float) -> dict[str, Any]:
    """Return explicitly non-qualifying provenance for analytical unit fixtures."""
    theoretical_dispersion(epsilon)
    dimension = _require_dimension(dimension)
    return {
        "kind": "synthetic_contract_fixture",
        "variant_id": f"Q023-SYNTHETIC-CONTRACT-{dimension}D-EPSILON-{epsilon:g}",
        "deck_path": "",
        "deck_sha256": "",
        "raw_geometry": _raw_geometry(dimension),
        "raw_artifacts": [],
        "retained_observable_contract": RETAINED_OBSERVABLE_CONTRACT,
    }


def _validated_variant_deck(variant: dict[str, Any]) -> Path:
    deck = Path(variant["deck"])
    expected_digest = variant["deck_sha256"]
    if not isinstance(expected_digest, str) or _SHA256.fullmatch(expected_digest) is None:
        raise ContractError("approved raw Bell deck digest pin is malformed")
    try:
        digest = _sha256(deck)
    except OSError as error:
        raise ContractError("approved raw Bell deck is missing") from error
    if digest != expected_digest:
        raise ContractError("approved raw Bell deck digest does not match the pinned value")
    return deck


def _authorized_artifact_root(artifact_root: Path | None) -> Path:
    if artifact_root is None:
        raise ContractError("raw Bell materialized provenance requires an explicit artifact root")
    root = Path(artifact_root)
    if not root.is_absolute():
        raise ContractError("raw Bell authorized artifact root must be absolute")
    try:
        root = root.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ContractError("raw Bell authorized artifact root is missing") from error
    if not root.is_dir():
        raise ContractError("raw Bell authorized artifact root must be a directory")
    return root


def _normalized_artifact_file(
    path: Path, artifact_root: Path, *, provenance_path: bool
) -> tuple[Path, str]:
    if provenance_path and path.is_absolute():
        raise ContractError("raw Bell artifact path must be normalized root-relative")
    candidate = path if path.is_absolute() else artifact_root / path
    try:
        resolved = candidate.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ContractError("raw Bell artifact path is missing") from error
    try:
        relative = resolved.relative_to(artifact_root)
    except ValueError as error:
        raise ContractError("raw Bell artifact path is outside the authorized root") from error
    if not resolved.is_file():
        raise ContractError("raw Bell artifact path must identify a file")
    normalized = relative.as_posix()
    if provenance_path and path.as_posix() != normalized:
        raise ContractError("raw Bell artifact path must be normalized root-relative")
    return resolved, normalized


def _source_local_materialized_provenance(
    dimension: int,
    epsilon: float,
    variant_id: str,
    paths: Sequence[Path],
    *,
    artifact_root: Path,
) -> dict[str, Any]:
    dimension = _require_dimension(dimension)
    variant = _APPROVED_SOURCE_LOCAL_RAW_VARIANTS.get(variant_id)
    if variant is None:
        raise ContractError("raw Bell variant is not an approved materialized source-local ID")
    if variant["dimension"] != dimension or variant["epsilon"] != epsilon:
        raise ContractError("raw Bell variant does not match the requested dimension and epsilon")
    deck = _validated_variant_deck(variant)
    root = _authorized_artifact_root(artifact_root)
    if not paths:
        raise ContractError("raw Bell source-local extraction requires retained artifact files")
    artifacts = []
    retained_paths = set()
    for path in paths:
        resolved, relative = _normalized_artifact_file(
            Path(path), root, provenance_path=False
        )
        if relative in retained_paths:
            raise ContractError("raw Bell artifact path is duplicated")
        retained_paths.add(relative)
        artifacts.append({"path": relative, "sha256": _sha256(resolved)})
    return {
        "kind": "source_local_materialized_variant",
        "variant_id": variant_id,
        "deck_path": str(deck.relative_to(REPO_ROOT)),
        "deck_sha256": variant["deck_sha256"],
        "raw_geometry": _raw_geometry(dimension),
        "raw_artifacts": artifacts,
        "retained_observable_contract": RETAINED_OBSERVABLE_CONTRACT,
    }


def _validate_raw_provenance(
    provenance: dict[str, Any],
    dimension: int,
    epsilon: float,
    *,
    artifact_root: Path | None = None,
) -> dict[str, Any]:
    dimension = _require_dimension(dimension)
    if not isinstance(provenance, dict) or set(provenance) != _RAW_PROVENANCE_KEYS:
        raise ContractError("raw Bell provenance keys do not match the contract")
    if provenance["retained_observable_contract"] != RETAINED_OBSERVABLE_CONTRACT:
        raise ContractError("raw Bell retained-observable contract mismatch")
    _validate_raw_geometry(provenance["raw_geometry"], dimension)
    artifacts = provenance["raw_artifacts"]
    if not isinstance(artifacts, list):
        raise ContractError("raw Bell provenance artifacts must be a list")

    kind = provenance["kind"]
    if kind == "synthetic_contract_fixture":
        expected = synthetic_contract_provenance(dimension, epsilon)
        if provenance != expected:
            raise ContractError("synthetic Bell contract fixture provenance mismatch")
        return provenance
    if kind != "source_local_materialized_variant":
        raise ContractError("raw Bell provenance kind is not recognized")

    variant = _APPROVED_SOURCE_LOCAL_RAW_VARIANTS.get(provenance["variant_id"])
    if variant is None:
        raise ContractError("raw Bell variant is not an approved materialized source-local ID")
    if variant["dimension"] != dimension or variant["epsilon"] != epsilon:
        raise ContractError("raw Bell variant does not match the retained record")
    deck = _validated_variant_deck(variant)
    if provenance["deck_path"] != str(deck.relative_to(REPO_ROOT)):
        raise ContractError("raw Bell provenance deck path mismatch")
    if provenance["deck_sha256"] != variant["deck_sha256"]:
        raise ContractError("raw Bell provenance deck digest mismatch")
    if not artifacts:
        raise ContractError("raw Bell source-local provenance requires retained artifacts")
    root = _authorized_artifact_root(artifact_root)
    paths = set()
    for artifact in artifacts:
        if not isinstance(artifact, dict) or set(artifact) != _RAW_ARTIFACT_KEYS:
            raise ContractError("raw Bell artifact provenance keys do not match the contract")
        if not isinstance(artifact["path"], str) or not artifact["path"]:
            raise ContractError("raw Bell artifact path is malformed")
        path, relative = _normalized_artifact_file(
            Path(artifact["path"]), root, provenance_path=True
        )
        if relative in paths:
            raise ContractError("raw Bell artifact path is duplicated")
        paths.add(relative)
        digest = artifact["sha256"]
        if not isinstance(digest, str) or _SHA256.fullmatch(digest) is None:
            raise ContractError("raw Bell artifact digest is malformed")
        if digest != _sha256(path):
            raise ContractError("raw Bell artifact digest mismatch")
    return provenance


def _mode_basis(dimension: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dimension = _require_dimension(dimension)
    raw = np.array(
        [1.0, 2.0 if dimension >= 2 else 0.0, 4.0 if dimension >= 3 else 0.0]
    )
    parallel = raw / np.linalg.norm(raw)
    transverse_a = (
        np.array([0.0, 1.0, 0.0])
        if dimension == 1
        else np.array([-parallel[1], parallel[0], 0.0])
    )
    transverse_a /= np.linalg.norm(transverse_a)
    transverse_b = np.cross(parallel, transverse_a)
    return parallel, transverse_a, transverse_b


def validate_candidate_deck(path: Path, expected_dimension: int) -> dict[str, Any]:
    """Validate paper values and retain the non-authorized local-run boundary."""
    expected_dimension = _require_dimension(expected_dimension)
    blocks = parse_athinput(path)
    for (block, name), expected in _EXPECTED_DECK_VALUES.items():
        measured = blocks.get(block, {}).get(name)
        if measured != expected:
            raise ContractError(
                f"{path}: {block}/{name}: expected {expected!r}, measured {measured!r}"
            )

    metadata = blocks["q023_paper_bell_linear"]
    if metadata.get("dimension") != str(expected_dimension):
        raise ContractError(f"{path}: unexpected paper dimension")
    expected_role = (
        "source_local_runnable_thin_2d3v_carrier_preparation_not_authorized"
        if expected_dimension == 1
        else "source_local_runnable_preparation_not_authorized"
    )
    if metadata.get("deck_role") != expected_role:
        raise ContractError(f"{path}: unexpected source-local preparation role")
    geometry = _EXPECTED_GEOMETRY[expected_dimension]
    nx = tuple(int(blocks["mesh"][f"nx{axis}"]) for axis in (1, 2, 3))
    extent = tuple(
        float(blocks["mesh"][f"x{axis}max"]) - float(blocks["mesh"][f"x{axis}min"])
        for axis in (1, 2, 3)
    )
    if nx != geometry["nx"]:
        raise ContractError(f"{path}: mesh cell-count contract mismatch")
    for axis, (measured, expected) in enumerate(
        zip(extent, geometry["extent"]), 1
    ):
        _require_close(f"{path}: x{axis} extent", measured, expected)
    active_dx = tuple(extent[axis] / nx[axis] for axis in range(expected_dimension))
    for axis, (measured, expected) in enumerate(
        zip(active_dx, geometry["active_dx"]), 1
    ):
        _require_close(f"{path}: dx{axis}", measured, expected)

    rho = float(metadata["rho"])
    bg = float(metadata["b_g"])
    ua = float(metadata["u_a"])
    wavelength = float(metadata["wavelength"])
    k0 = float(metadata["k0"])
    omega = float(metadata["omega"])
    epsilon = float(metadata["epsilon_default"])
    configured_epsilon = float(metadata["epsilon"])
    stream = np.array(
        [
            float(blocks["particles"]["cr_vx0"]),
            float(blocks["particles"]["cr_vy0"]),
            float(blocks["particles"]["cr_vz0"]),
        ]
    )
    vcr = float(np.linalg.norm(stream))
    light_speed = float(blocks["particles"]["pic_cr_light_speed"])
    q_over_m = float(blocks["species0"]["charge"]) / float(blocks["species0"]["mass"])
    ppc = float(blocks["particles"]["ppc"])
    qscale = float(blocks["particles"]["deposit_qscale"])
    jcr = ppc * qscale * float(blocks["species0"]["charge"]) * vcr

    _require_close(f"{path}: U_A", ua, bg / math.sqrt(rho))
    _require_close(f"{path}: k0", k0, 2.0 * math.pi / wavelength)
    _require_close(f"{path}: Omega", omega, 1.0e-6 * k0 * ua)
    _require_close(f"{path}: species q/m", q_over_m, omega / bg)
    _require_close(f"{path}: epsilon", epsilon, ua / vcr)
    _require_close(f"{path}: configured epsilon", configured_epsilon, epsilon)
    _require_close(f"{path}: C", light_speed, 1.0e3 * vcr)
    _require_close(f"{path}: j_CR", jcr, 2.0 * bg * light_speed * k0)
    expected_stream = vcr * _mode_basis(expected_dimension)[0]
    for axis, (measured, expected) in enumerate(zip(stream, expected_stream), 1):
        _require_close(f"{path}: diagonal CR stream x{axis}", measured, float(expected))

    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "dimension": expected_dimension,
        "epsilon_default": epsilon,
        "active_dx": list(active_dx),
        "launch_status": (
            "source_local_runnable_thin_2d3v_carrier_preparation_only_not_authorized"
            if expected_dimension == 1
            else "source_local_runnable_preparation_only_not_authorized"
        ),
        "carrier_semantics": (
            "physical_1d_transverse_invariant_thin_2d3v_nx2_4"
            if expected_dimension == 1
            else "native_mesh_dimension"
        ),
    }


def validate_source_local_candidate_decks() -> list[dict[str, Any]]:
    """Validate all three paper-value candidate decks without launching AthenaK."""
    return [
        validate_candidate_deck(DECKS[dimension], dimension)
        for dimension in DIMENSIONS
    ]


def theoretical_dispersion(epsilon: float) -> tuple[float, float]:
    """Return normalized (phase frequency, growth rate) at k=k0."""
    if epsilon not in EPSILON_VALUES:
        raise ContractError(f"epsilon {epsilon!r} is outside the preregistered grid")
    return epsilon, math.sqrt(1.0 - epsilon * epsilon)


def _fixed_interval_phase_trace(
    normalized_time: np.ndarray, right_mode: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    if normalized_time.size != right_mode.size or normalized_time.size < 2:
        raise ContractError("phase extraction requires matching trace arrays")
    if np.any(np.diff(normalized_time) <= 0.0):
        raise ContractError("phase extraction requires strictly increasing time")
    first = math.ceil(normalized_time[0] / PHASE_INTERVAL) * PHASE_INTERVAL
    boundaries = np.arange(
        first, normalized_time[-1] + 1.0e-12, PHASE_INTERVAL, dtype=float
    )
    if boundaries.size < 2:
        raise ContractError("trace does not cover one fixed paper phase interval")
    phase = np.interp(boundaries, normalized_time, np.unwrap(np.angle(right_mode)))
    return np.diff(boundaries), np.diff(phase)


def _spatial_modes_from_components(
    dataset: dict[str, Any],
    dimension: int,
    raw_geometry: dict[str, Any],
    component_names: Sequence[str],
    *,
    label: str,
) -> tuple[complex, complex]:
    parallel, transverse_a, transverse_b = _mode_basis(dimension)
    coordinates = [
        np.asarray(dataset[f"x{axis}v"], dtype=float) for axis in (1, 2, 3)
    ]
    if any(values.ndim != 1 or not np.all(np.isfinite(values))
           for values in coordinates):
        raise ContractError(
            "raw mhd_w_bcc coordinates must be finite one-dimensional arrays"
        )
    nx = tuple(raw_geometry["nx"])
    xmin = tuple(raw_geometry["xmin"])
    extent = tuple(raw_geometry["extent"])
    if tuple(values.size for values in coordinates) != nx:
        raise ContractError("raw mhd_w_bcc geometry does not match the approved variant")
    for axis, values in enumerate(coordinates):
        expected = xmin[axis] + (np.arange(nx[axis], dtype=float) + 0.5) * (
            extent[axis] / nx[axis]
        )
        if not np.allclose(
            values, expected, rtol=0.0, atol=RAW_GEOMETRY_ABSOLUTE_TOLERANCE
        ):
            raise ContractError(
                "raw mhd_w_bcc geometry does not match the approved variant"
            )
    expected_shape = tuple(values.size for values in reversed(coordinates))
    if any(name not in dataset for name in component_names):
        raise ContractError(f"raw mhd_w_bcc {label} components are incomplete")
    vector = np.stack(
        [np.asarray(dataset[name], dtype=float) for name in component_names]
    )
    if vector.shape[1:] != expected_shape or not np.all(np.isfinite(vector)):
        raise ContractError(f"raw mhd_w_bcc {label} shape or values are invalid")
    x3, x2, x1 = np.meshgrid(
        coordinates[2], coordinates[1], coordinates[0], indexing="ij"
    )
    phase = K0 * (parallel[0]*x1 + parallel[1]*x2 + parallel[2]*x3)
    fourier_weight = np.exp(-1.0j * phase)
    mode_a = np.mean(np.tensordot(transverse_a, vector, axes=1) * fourier_weight)
    mode_b = np.mean(np.tensordot(transverse_b, vector, axes=1) * fourier_weight)
    return 0.5*(mode_a - 1.0j*mode_b), 0.5*(mode_a + 1.0j*mode_b)


def _spatial_modes_from_dataset(
    dataset: dict[str, Any], dimension: int, raw_geometry: dict[str, Any]
) -> tuple[complex, complex, complex, complex]:
    magnetic = _spatial_modes_from_components(
        dataset, dimension, raw_geometry, ("bcc1", "bcc2", "bcc3"), label="magnetic"
    )
    velocity = _spatial_modes_from_components(
        dataset, dimension, raw_geometry, ("velx", "vely", "velz"), label="velocity"
    )
    return magnetic + velocity


def _paper_literal_velocity_observables_from_dataset(
    dataset: dict[str, Any], dimension: int, raw_geometry: dict[str, Any]
) -> tuple[complex, float]:
    """Return the Section 5.2 delta-u-y sine fit and volume-averaged |delta u|."""
    parallel, _, _ = _mode_basis(dimension)
    coordinates = [
        np.asarray(dataset[f"x{axis}v"], dtype=float) for axis in (1, 2, 3)
    ]
    nx = tuple(raw_geometry["nx"])
    xmin = tuple(raw_geometry["xmin"])
    extent = tuple(raw_geometry["extent"])
    expected_shape = tuple(values.size for values in reversed(coordinates))
    if tuple(values.size for values in coordinates) != nx:
        raise ContractError("raw mhd_w_bcc geometry does not match the approved variant")
    for axis, values in enumerate(coordinates):
        expected = xmin[axis] + (np.arange(nx[axis], dtype=float) + 0.5) * (
            extent[axis] / nx[axis]
        )
        if not np.allclose(
            values, expected, rtol=0.0, atol=RAW_GEOMETRY_ABSOLUTE_TOLERANCE
        ):
            raise ContractError(
                "raw mhd_w_bcc geometry does not match the approved variant"
            )
    velocity = np.stack(
        [np.asarray(dataset[name], dtype=float) for name in ("velx", "vely", "velz")]
    )
    if velocity.shape[1:] != expected_shape or not np.all(np.isfinite(velocity)):
        raise ContractError("raw mhd_w_bcc velocity shape or values are invalid")
    x3, x2, x1 = np.meshgrid(
        coordinates[2], coordinates[1], coordinates[0], indexing="ij"
    )
    spatial_phase = K0 * (
        parallel[0]*x1 + parallel[1]*x2 + parallel[2]*x3
    )
    design = np.column_stack(
        (np.cos(spatial_phase).ravel(), np.sin(spatial_phase).ravel())
    )
    coefficients, _, rank, _ = np.linalg.lstsq(
        design, velocity[1].ravel(), rcond=None
    )
    if rank != 2 or not np.all(np.isfinite(coefficients)):
        raise ContractError("paper-literal delta_u_y spatial sine fit is singular")
    phase_fit = complex(coefficients[0], -coefficients[1])
    if abs(phase_fit) <= np.finfo(float).tiny:
        raise ContractError("paper-literal delta_u_y spatial sine fit has zero amplitude")
    volume_averaged_abs_delta_u = float(
        np.mean(np.sqrt(np.sum(velocity*velocity, axis=0)))
    )
    if (
        not math.isfinite(volume_averaged_abs_delta_u)
        or volume_averaged_abs_delta_u <= 0.0
    ):
        raise ContractError("paper-literal volume-averaged |delta u| must be positive")
    return phase_fit, volume_averaged_abs_delta_u


def extract_trace_record_from_datasets(
    dimension: int,
    epsilon: float,
    datasets: Sequence[dict[str, Any]],
    *,
    raw_provenance: dict[str, Any],
    artifact_root: Path | None = None,
) -> dict[str, Any]:
    """Extract one paper Bell trace record from raw mhd_w_bcc snapshot datasets."""
    theoretical_dispersion(epsilon)
    provenance = _validate_raw_provenance(
        raw_provenance, dimension, epsilon, artifact_root=artifact_root
    )
    rows = []
    for dataset in datasets:
        if "Time" not in dataset:
            raise ContractError("raw mhd_w_bcc dataset is missing Time")
        time = float(dataset["Time"])
        if not math.isfinite(time):
            raise ContractError("raw mhd_w_bcc dataset Time must be finite")
        right, left, velocity_right, velocity_left = _spatial_modes_from_dataset(
            dataset, dimension, provenance["raw_geometry"]
        )
        paper_phase_fit, paper_abs_delta_u = (
            _paper_literal_velocity_observables_from_dataset(
                dataset, dimension, provenance["raw_geometry"]
            )
        )
        rows.append(
            (
                time*K0*U_A,
                right,
                left,
                velocity_right,
                velocity_left,
                paper_phase_fit,
                paper_abs_delta_u,
            )
        )
    rows.sort(key=lambda item: item[0])
    if not rows:
        raise ContractError("raw mhd_w_bcc extraction requires snapshots")
    time = np.asarray([item[0] for item in rows], dtype=float)
    right = np.asarray([item[1] for item in rows], dtype=complex)
    left = np.asarray([item[2] for item in rows], dtype=complex)
    velocity_right = np.asarray([item[3] for item in rows], dtype=complex)
    velocity_left = np.asarray([item[4] for item in rows], dtype=complex)
    paper_phase_fit = np.asarray([item[5] for item in rows], dtype=complex)
    paper_abs_delta_u = np.asarray([item[6] for item in rows], dtype=float)
    interval, change = _fixed_interval_phase_trace(time, right)
    velocity_interval, velocity_change = _fixed_interval_phase_trace(
        time, velocity_right
    )
    paper_interval, paper_change = _fixed_interval_phase_trace(time, paper_phase_fit)
    return {
        "dimension": dimension,
        "epsilon": epsilon,
        "raw_provenance": provenance,
        "normalized_time": time.tolist(),
        "right_mode_real": right.real.tolist(),
        "right_mode_imag": right.imag.tolist(),
        "left_mode_real": left.real.tolist(),
        "left_mode_imag": left.imag.tolist(),
        "velocity_right_mode_real": velocity_right.real.tolist(),
        "velocity_right_mode_imag": velocity_right.imag.tolist(),
        "velocity_left_mode_real": velocity_left.real.tolist(),
        "velocity_left_mode_imag": velocity_left.imag.tolist(),
        "phase_interval": interval.tolist(),
        "phase_change": change.tolist(),
        "velocity_phase_interval": velocity_interval.tolist(),
        "velocity_phase_change": velocity_change.tolist(),
        "paper_delta_u_y_sine_fit_real": paper_phase_fit.real.tolist(),
        "paper_delta_u_y_sine_fit_imag": paper_phase_fit.imag.tolist(),
        "paper_delta_u_y_phase_interval": paper_interval.tolist(),
        "paper_delta_u_y_phase_change": paper_change.tolist(),
        "paper_volume_averaged_abs_delta_u": paper_abs_delta_u.tolist(),
    }


def extract_trace_record_from_binary_files(
    dimension: int,
    epsilon: float,
    paths: Sequence[Path],
    variant_id: str,
    reader: Callable[[str], dict[str, Any]] | None = None,
    *,
    artifact_root: Path,
) -> dict[str, Any]:
    """Read raw Athena binary snapshots and extract one paper Bell trace record."""
    provenance = _source_local_materialized_provenance(
        dimension, epsilon, variant_id, paths, artifact_root=artifact_root
    )
    root = _authorized_artifact_root(artifact_root)
    normalized_paths = [root / artifact["path"]
                        for artifact in provenance["raw_artifacts"]]
    if reader is None:
        module_path = REPO_ROOT / "vis/python/bin_convert_new.py"
        spec = importlib.util.spec_from_file_location("q023_bin_convert_new", module_path)
        if spec is None or spec.loader is None:
            raise ContractError("unable to load the raw Athena binary reader")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        reader = module.read_binary_as_athdf
    return extract_trace_record_from_datasets(
        dimension,
        epsilon,
        [reader(str(path)) for path in normalized_paths],
        raw_provenance=provenance,
        artifact_root=root,
    )


def _finite_array(record: dict[str, Any], key: str) -> np.ndarray:
    values = np.asarray(record[key], dtype=float)
    if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
        raise ContractError(f"{key} must be a finite nonempty one-dimensional array")
    return values


def _within_both_tolerances(measured: float, expected: float) -> bool:
    absolute_error = abs(measured - expected)
    relative_error = absolute_error / abs(expected)
    return absolute_error <= ABSOLUTE_TOLERANCE and relative_error <= RELATIVE_TOLERANCE


def _fit_growth(
    time: np.ndarray, amplitude: np.ndarray, expected_growth: float, *, label: str
) -> tuple[float, float, np.ndarray]:
    fit_mask = (expected_growth * time >= 1.0) & (expected_growth * time <= 5.0)
    if np.count_nonzero(fit_mask) < MIN_GROWTH_SNAPSHOTS:
        raise ContractError("the fixed theory growth window requires >= 8 snapshots")
    if np.any(amplitude[fit_mask] <= 0.0):
        raise ContractError(f"{label} right-polarized mode amplitude must remain positive")
    fit_time = time[fit_mask]
    fit_log_amplitude = np.log(amplitude[fit_mask])
    coeff = np.polyfit(fit_time, fit_log_amplitude, 1)
    fitted = np.polyval(coeff, fit_time)
    residual = fit_log_amplitude - fitted
    ss_res = float(np.sum(residual * residual))
    ss_tot = float(np.sum((fit_log_amplitude - np.mean(fit_log_amplitude)) ** 2))
    growth_r2 = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    return float(coeff[0]), growth_r2, fit_mask


def _analyze_record(
    record: dict[str, Any], *, artifact_root: Path | None = None
) -> dict[str, Any]:
    if set(record) != _TRACE_KEYS:
        raise ContractError("Bell extracted-trace record keys do not match the contract")
    dimension = record["dimension"]
    if type(dimension) is not int or dimension not in DIMENSIONS:
        raise ContractError("dimension must be one of the preregistered integer values")
    if type(record["epsilon"]) not in (int, float):
        raise ContractError("epsilon must be a preregistered JSON number")
    epsilon = float(record["epsilon"])
    expected_phase, expected_growth = theoretical_dispersion(epsilon)
    provenance = _validate_raw_provenance(
        record["raw_provenance"], dimension, epsilon, artifact_root=artifact_root
    )

    time = _finite_array(record, "normalized_time")
    if time.size < MIN_GROWTH_SNAPSHOTS or np.any(np.diff(time) <= 0.0):
        raise ContractError("normalized_time must be strictly increasing with >= 8 rows")
    magnetic_mode_arrays = [
        _finite_array(record, key)
        for key in (
            "right_mode_real",
            "right_mode_imag",
            "left_mode_real",
            "left_mode_imag",
        )
    ]
    velocity_mode_arrays = [
        _finite_array(record, key)
        for key in (
            "velocity_right_mode_real",
            "velocity_right_mode_imag",
            "velocity_left_mode_real",
            "velocity_left_mode_imag",
        )
    ]
    if any(values.size != time.size
           for values in magnetic_mode_arrays + velocity_mode_arrays):
        raise ContractError("mode arrays must match normalized_time length")
    right = magnetic_mode_arrays[0] + 1.0j * magnetic_mode_arrays[1]
    left = magnetic_mode_arrays[2] + 1.0j * magnetic_mode_arrays[3]
    velocity_right = velocity_mode_arrays[0] + 1.0j * velocity_mode_arrays[1]
    velocity_left = velocity_mode_arrays[2] + 1.0j * velocity_mode_arrays[3]
    right_amplitude = np.abs(right)
    left_amplitude = np.abs(left)
    velocity_right_amplitude = np.abs(velocity_right)
    velocity_left_amplitude = np.abs(velocity_left)
    paper_abs_delta_u = _finite_array(record, "paper_volume_averaged_abs_delta_u")
    if paper_abs_delta_u.size != time.size:
        raise ContractError(
            "paper-literal volume-averaged |delta u| must match normalized_time length"
        )
    paper_delta_u_y_sine_fit = (
        _finite_array(record, "paper_delta_u_y_sine_fit_real")
        + 1.0j * _finite_array(record, "paper_delta_u_y_sine_fit_imag")
    )
    if paper_delta_u_y_sine_fit.size != time.size:
        raise ContractError(
            "paper-literal delta_u_y sine-fit trace must match normalized_time length"
        )
    measured_growth, growth_r2, fit_mask = _fit_growth(
        time, right_amplitude, expected_growth, label="magnetic"
    )
    velocity_measured_growth, velocity_growth_r2, velocity_fit_mask = _fit_growth(
        time, velocity_right_amplitude, expected_growth, label="fluid-velocity"
    )
    paper_measured_growth, paper_growth_r2, paper_fit_mask = _fit_growth(
        time, paper_abs_delta_u, expected_growth, label="paper-literal fluid-velocity"
    )
    if not np.array_equal(fit_mask, velocity_fit_mask):
        raise ContractError("magnetic and fluid-velocity growth windows must match")
    if not np.array_equal(fit_mask, paper_fit_mask):
        raise ContractError("magnetic and paper-literal growth windows must match")

    phase_interval = _finite_array(record, "phase_interval")
    phase_change = _finite_array(record, "phase_change")
    if phase_interval.size != phase_change.size:
        raise ContractError("phase interval and phase change arrays must have equal length")
    if not np.allclose(phase_interval, PHASE_INTERVAL, rtol=0.0, atol=1.0e-14):
        raise ContractError("phase intervals must equal pi/(k0*U_A)")
    extracted_interval, extracted_change = _fixed_interval_phase_trace(time, right)
    if phase_interval.size != extracted_interval.size or not np.allclose(
        phase_interval, extracted_interval, rtol=0.0, atol=1.0e-12
    ):
        raise ContractError("phase intervals do not match the retained mode trace")
    if phase_change.size != extracted_change.size or not np.allclose(
        phase_change, extracted_change, rtol=0.0, atol=1.0e-10
    ):
        raise ContractError("phase changes do not match the retained mode trace")
    measured_phase = float(np.mean(extracted_change / extracted_interval))
    velocity_phase_interval = _finite_array(record, "velocity_phase_interval")
    velocity_phase_change = _finite_array(record, "velocity_phase_change")
    if velocity_phase_interval.size != velocity_phase_change.size:
        raise ContractError(
            "fluid-velocity phase interval and phase change arrays must have equal length"
        )
    if not np.allclose(
        velocity_phase_interval, PHASE_INTERVAL, rtol=0.0, atol=1.0e-14
    ):
        raise ContractError("fluid-velocity phase intervals must equal pi/(k0*U_A)")
    extracted_velocity_interval, extracted_velocity_change = (
        _fixed_interval_phase_trace(time, velocity_right)
    )
    if (
        velocity_phase_interval.size != extracted_velocity_interval.size
        or not np.allclose(
            velocity_phase_interval,
            extracted_velocity_interval,
            rtol=0.0,
            atol=1.0e-12,
        )
    ):
        raise ContractError(
            "fluid-velocity phase intervals do not match the retained mode trace"
        )
    if (
        velocity_phase_change.size != extracted_velocity_change.size
        or not np.allclose(
            velocity_phase_change,
            extracted_velocity_change,
            rtol=0.0,
            atol=1.0e-10,
        )
    ):
        raise ContractError(
            "fluid-velocity phase changes do not match the retained mode trace"
        )
    velocity_measured_phase = float(
        np.mean(extracted_velocity_change / extracted_velocity_interval)
    )
    paper_phase_interval = _finite_array(record, "paper_delta_u_y_phase_interval")
    paper_phase_change = _finite_array(record, "paper_delta_u_y_phase_change")
    if paper_phase_interval.size != paper_phase_change.size:
        raise ContractError(
            "paper-literal delta_u_y phase interval and phase change arrays "
            "must have equal length"
        )
    if not np.allclose(
        paper_phase_interval, PHASE_INTERVAL, rtol=0.0, atol=1.0e-14
    ):
        raise ContractError(
            "paper-literal delta_u_y phase intervals must equal pi/(k0*U_A)"
        )
    extracted_paper_interval, extracted_paper_change = _fixed_interval_phase_trace(
        time, paper_delta_u_y_sine_fit
    )
    if paper_phase_interval.size != extracted_paper_interval.size or not np.allclose(
        paper_phase_interval, extracted_paper_interval, rtol=0.0, atol=1.0e-12
    ):
        raise ContractError(
            "paper-literal delta_u_y phase intervals do not match the retained sine-fit trace"
        )
    if paper_phase_change.size != extracted_paper_change.size or not np.allclose(
        paper_phase_change, extracted_paper_change, rtol=0.0, atol=1.0e-10
    ):
        raise ContractError(
            "paper-literal delta_u_y phase changes do not match the retained sine-fit trace"
        )
    paper_measured_phase = float(
        np.mean(extracted_paper_change / extracted_paper_interval)
    )
    final_fit_index = int(np.flatnonzero(fit_mask)[-1])
    polarization_ratio = float(
        right_amplitude[final_fit_index]
        / max(float(left_amplitude[final_fit_index]), np.finfo(float).tiny)
    )
    velocity_polarization_ratio = float(
        velocity_right_amplitude[final_fit_index]
        / max(float(velocity_left_amplitude[final_fit_index]), np.finfo(float).tiny)
    )
    expected_velocity_to_magnetic_ratio = complex(-epsilon, -expected_growth)
    velocity_to_magnetic_ratio = velocity_right[fit_mask] / right[fit_mask]
    velocity_magnetic_ratio_max_absolute_error = float(
        np.max(np.abs(velocity_to_magnetic_ratio - expected_velocity_to_magnetic_ratio))
    )

    growth_pass = _within_both_tolerances(measured_growth, expected_growth)
    phase_pass = _within_both_tolerances(measured_phase, expected_phase)
    polarization_pass = polarization_ratio >= MIN_RIGHT_TO_LEFT_RATIO
    velocity_growth_pass = _within_both_tolerances(
        velocity_measured_growth, expected_growth
    )
    velocity_phase_pass = _within_both_tolerances(
        velocity_measured_phase, expected_phase
    )
    velocity_polarization_pass = (
        velocity_polarization_ratio >= MIN_RIGHT_TO_LEFT_RATIO
    )
    velocity_magnetic_ratio_pass = (
        velocity_magnetic_ratio_max_absolute_error
        <= SOURCE_LOCAL_VELOCITY_MAGNETIC_RATIO_ABSOLUTE_TOLERANCE
    )
    paper_growth_pass = _within_both_tolerances(
        paper_measured_growth, expected_growth
    )
    paper_phase_pass = _within_both_tolerances(paper_measured_phase, expected_phase)
    return {
        "dimension": dimension,
        "epsilon": epsilon,
        "raw_provenance_kind": provenance["kind"],
        "raw_variant_id": provenance["variant_id"],
        "expected_growth_rate_over_k0_ua": expected_growth,
        "measured_growth_rate_over_k0_ua": measured_growth,
        "growth_fit_r2": growth_r2,
        "expected_phase_frequency_over_k0_ua": expected_phase,
        "measured_phase_frequency_over_k0_ua": measured_phase,
        "right_to_left_amplitude_ratio": polarization_ratio,
        "velocity_measured_growth_rate_over_k0_ua": velocity_measured_growth,
        "velocity_growth_fit_r2": velocity_growth_r2,
        "velocity_measured_phase_frequency_over_k0_ua": velocity_measured_phase,
        "velocity_right_to_left_amplitude_ratio": velocity_polarization_ratio,
        "expected_velocity_to_magnetic_ratio_real":
            expected_velocity_to_magnetic_ratio.real,
        "expected_velocity_to_magnetic_ratio_imag":
            expected_velocity_to_magnetic_ratio.imag,
        "velocity_magnetic_ratio_max_absolute_error":
            velocity_magnetic_ratio_max_absolute_error,
        "velocity_magnetic_ratio_diagnostic_limit":
            SOURCE_LOCAL_VELOCITY_MAGNETIC_RATIO_ABSOLUTE_TOLERANCE,
        "paper_literal_measured_growth_rate_over_k0_ua": paper_measured_growth,
        "paper_literal_growth_fit_r2": paper_growth_r2,
        "paper_literal_measured_phase_frequency_over_k0_ua": paper_measured_phase,
        "growth_pass": growth_pass,
        "phase_pass": phase_pass,
        "polarization_pass": polarization_pass,
        "velocity_growth_pass": velocity_growth_pass,
        "velocity_phase_pass": velocity_phase_pass,
        "velocity_polarization_pass": velocity_polarization_pass,
        "velocity_magnetic_ratio_pass": velocity_magnetic_ratio_pass,
        "paper_literal_growth_pass": paper_growth_pass,
        "paper_literal_phase_pass": paper_phase_pass,
        "velocity_cross_check_qualification_effect":
            "nonqualifying_source_local_diagnostic_only",
        "retained_velocity_observable_status":
            "paper_literal_delta_u_y_phase_and_volume_averaged_abs_delta_u_frozen_"
            "projected_mode_cross_check_retained",
        "section52_qualification_eligible": False,
        "passed": (
            growth_pass
            and phase_pass
            and polarization_pass
            and velocity_growth_pass
            and velocity_phase_pass
            and velocity_polarization_pass
            and velocity_magnetic_ratio_pass
            and paper_growth_pass
            and paper_phase_pass
        ),
    }


def analyze_trace_bundle(
    bundle: dict[str, Any], *, artifact_root: Path | None = None
) -> dict[str, Any]:
    """Analyze an exact extracted-trace grid without making a qualification claim."""
    if set(bundle) != {"schema_version", "campaign_id", "qualifying_seed", "records"}:
        raise ContractError("Bell extracted-trace bundle keys do not match the contract")
    if (
        type(bundle["schema_version"]) is not int
        or bundle["schema_version"] != TRACE_SCHEMA_VERSION
        or bundle["campaign_id"] != CAMPAIGN_ID
    ):
        raise ContractError("Bell extracted-trace bundle identity mismatch")
    qualifying_seed = bundle["qualifying_seed"]
    if type(qualifying_seed) is not int or qualifying_seed not in QUALIFYING_SEEDS:
        raise ContractError("Bell extracted-trace bundle seed is not preregistered")
    if not isinstance(bundle["records"], list):
        raise ContractError("Bell extracted-trace records must be a list")

    expected_keys = {(dimension, epsilon) for dimension in DIMENSIONS
                     for epsilon in EPSILON_VALUES}
    measured_keys = []
    reports = []
    for record in bundle["records"]:
        if not isinstance(record, dict):
            raise ContractError("Bell extracted-trace rows must be objects")
        report = _analyze_record(record, artifact_root=artifact_root)
        key = (report["dimension"], report["epsilon"])
        if key in measured_keys:
            raise ContractError(f"duplicate Bell extracted-trace row {key!r}")
        measured_keys.append(key)
        reports.append(report)
    if set(measured_keys) != expected_keys:
        raise ContractError("Bell extracted-trace grid is incomplete or contains extra rows")

    reports.sort(key=lambda item: (item["dimension"], item["epsilon"]))
    provenance_kinds = {report["raw_provenance_kind"] for report in reports}
    if len(provenance_kinds) != 1:
        raise ContractError("Bell extracted-trace bundle cannot mix provenance kinds")
    analysis_input_kind = next(iter(provenance_kinds))
    synthetic_fixture_analysis = analysis_input_kind == "synthetic_contract_fixture"
    scientific_contract_pass = all(report["passed"] for report in reports)
    materialized_source_local_candidate_pass = (
        not synthetic_fixture_analysis and scientific_contract_pass
    )
    qualification_effect = (
        "none_synthetic_contract_fixture_test_only_not_a_materialized_source_local_"
        "candidate_pass"
        if synthetic_fixture_analysis
        else (
            "source_local_combined_magnetic_and_fluid_velocity_mode_candidate_"
            "cross_check_only_section52_qualification_blocked_clean_candidate_"
            "binding_convergence_registered_frontier_execution_independent_recompute_"
            "and_external_review_required"
        )
    )
    return {
        "schema_version": TRACE_SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "qualifying_seed": qualifying_seed,
        "analysis_input_kind": analysis_input_kind,
        "analysis_scope": (
            "synthetic_contract_fixture_test_only"
            if synthetic_fixture_analysis
            else "source_local_materialized_combined_magnetic_and_fluid_velocity_"
                 "mode_candidate"
        ),
        "synthetic_fixture_analysis": synthetic_fixture_analysis,
        "scientific_contract_pass": scientific_contract_pass,
        "materialized_source_local_candidate_pass": materialized_source_local_candidate_pass,
        "qualification_effect": qualification_effect,
        "section52_qualification_eligible": False,
        "deck_contracts": validate_source_local_candidate_decks(),
        "record_count": len(reports),
        "records": reports,
        "passed": materialized_source_local_candidate_pass,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extract-record", action="store_true")
    parser.add_argument("--dimension", type=int)
    parser.add_argument("--epsilon", type=float)
    parser.add_argument("--variant-id")
    parser.add_argument("--artifact-root", type=Path)
    parser.add_argument("paths", type=Path, nargs="+")
    args = parser.parse_args()
    if args.extract_record:
        if (args.dimension is None or args.epsilon is None or
                args.variant_id is None or args.artifact_root is None):
            parser.error(
                "--extract-record requires --dimension, --epsilon, --variant-id "
                "and --artifact-root"
            )
        result = extract_trace_record_from_binary_files(
            args.dimension,
            args.epsilon,
            args.paths,
            args.variant_id,
            artifact_root=args.artifact_root,
        )
    else:
        if (len(args.paths) != 1 or args.dimension is not None or
                args.epsilon is not None or args.variant_id is not None):
            parser.error("bundle analysis requires exactly one trace-bundle path")
        bundle = json.loads(args.paths[0].read_text(encoding="utf-8"))
        result = analyze_trace_bundle(bundle, artifact_root=args.artifact_root)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
