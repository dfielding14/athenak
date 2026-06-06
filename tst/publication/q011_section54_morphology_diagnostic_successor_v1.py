#!/usr/bin/env python3
"""Publication-bound Q-011 Section 5.4 morphology diagnostics at nominal t=500.

This additive source-local successor consumes already decoded production-science
mesh and particle products.  It does not launch work, mutate policy, inspect an
unbound qualifying campaign, or close a scientific claim.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from functools import wraps
import hashlib
import json
import math
from numbers import Real
from pathlib import Path, PurePosixPath
import re
from types import MappingProxyType
from typing import Any

import numpy as np

if __package__:
    from . import q011_section54_artifacts as artifact_helpers
    from . import q011_section54_production_science_successor_v1 as science
else:
    import q011_section54_artifacts as artifact_helpers
    import q011_section54_production_science_successor_v1 as science


SCHEMA_VERSION = 1
SUCCESSOR_ID = "q011_section54_morphology_diagnostic_successor_v1"
NOMINAL_SLOT_TIME = 500.0
FIXED_PROFILE_CELL_SIZE_C_OVER_OMEGA_PI = 12.0
REPO_ROOT = Path(__file__).resolve().parents[2]

QUALIFICATION_EFFECT = (
    "source_local_publication_bound_diagnostic_only_no_launch_no_policy_"
    "authorization_no_claim_closure"
)
_MANIFEST_ARTIFACT_ROLE = "publication_bound_t500_morphology_diagnostic_successor"
_MANIFEST_EVIDENCE_BOUNDARY = (
    "manifest_and_profiles_bind_but_do_not_replace_immutable_raw_artifacts"
)
_PROFILE_EVIDENCE_BOUNDARY = "derived_profile_does_not_replace_retained_raw_artifacts"
_MANIFEST_PROFILE_CONTRACT = {
    "grid_matching_rule": (
        "conservatively_restrict_to_12_c_over_omega_pi_then_full_y_average"
    ),
    "particle_distribution_weight": "active_cr_kinetic_energy",
    "current_representation": "deposited_J_CR_over_c_in_lab_and_gas_frames",
    "magnetic_profile": "cell_centered_bcc_components_and_magnitude",
    "refinement_overlay": "exact_leaf_block_rectangles_and_source_level_area_fractions",
}
AUTHORIZATION: Mapping[str, bool] = MappingProxyType(
    {
        "launch_authorized": False,
        "policy_mutation_authorized": False,
        "claim_closure_authorized": False,
        "qualifying_output_inspection_authorized": False,
    }
)

EXPECTED_RAW_PRODUCTS = (
    "mhd_w_bcc",
    "prtcl_rho",
    "prtcl_jx",
    "prtcl_jy",
    "prtcl_jz",
    "mhd_j2",
    "prtcl_all",
)
EXPECTED_PVTK_SCALARS = (
    "birth_time",
    "cr_source",
    "deltaf_f0",
    "deltaf_weight",
    "gid",
    "macro_weight",
    "ptag",
    "species",
)
EXPECTED_PVTK_VECTORS = ("vel",)

MANIFEST_MEMBER = "morphology_manifest.json"
INVENTORY_MEMBER = "artifact_inventory.json"
PARTICLE_PROFILE_MEMBER = "profiles/particle_energy_weighted_spatial_profile_t500.json"
CURRENT_PROFILE_MEMBER = "profiles/deposited_current_spatial_profile_t500.json"
MAGNETIC_PROFILE_MEMBER = "profiles/magnetic_field_spatial_profile_t500.json"
REFINEMENT_PROFILE_MEMBER = "profiles/refinement_spatial_profile_t500.json"
REFINEMENT_OVERLAY_MEMBER = "overlays/refinement_leaf_block_overlay_t500.json"
DERIVED_MEMBERS = (
    MANIFEST_MEMBER,
    PARTICLE_PROFILE_MEMBER,
    CURRENT_PROFILE_MEMBER,
    MAGNETIC_PROFILE_MEMBER,
    REFINEMENT_PROFILE_MEMBER,
    REFINEMENT_OVERLAY_MEMBER,
)
ALL_MEMBERS = (*DERIVED_MEMBERS, INVENTORY_MEMBER)

_SOURCE_PATHS = (
    "inputs/publication/"
    "pic_parallel_shock_section54_production_science_successor_v1_vl2_tsc.athinput",
    "tst/publication/analyze_q011_section54_outputs.py",
    "tst/publication/q011_section54_artifacts.py",
    "tst/publication/q011_section54_model.py",
    "tst/publication/q011_section54_morphology_diagnostic_successor_v1.py",
    "tst/publication/q011_section54_particles.py",
    "tst/publication/q011_section54_production_science_successor_v1.py",
)
_PRODUCT_PATH_TOKENS: Mapping[str, str] = MappingProxyType(
    {
        "mhd_w_bcc": ".mhd_w_bcc.",
        "prtcl_rho": ".prtcl_rho.",
        "prtcl_jx": ".prtcl_jx.",
        "prtcl_jy": ".prtcl_jy.",
        "prtcl_jz": ".prtcl_jz.",
        "mhd_j2": ".mhd_j2.",
        "prtcl_all": ".prtcl_all.",
    }
)
_SHA256 = re.compile(r"[0-9a-f]{64}")
_ATTEMPT_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:-]{0,255}")


class MorphologyDiagnosticError(ValueError):
    """Raised when a t=500 morphology artifact cannot be trusted."""


_UNDERLYING_EXCEPTIONS = (
    artifact_helpers.DerivedArtifactError,
    science.ProductionScienceError,
    OSError,
    UnicodeError,
    KeyError,
    IndexError,
    AttributeError,
    TypeError,
    ValueError,
    OverflowError,
    FloatingPointError,
    RecursionError,
)


def _public_contract(label: str):
    def decorate(function):
        @wraps(function)
        def wrapped(*args: object, **kwargs: object):
            try:
                return function(*args, **kwargs)
            except MorphologyDiagnosticError:
                raise
            except _UNDERLYING_EXCEPTIONS as error:
                raise MorphologyDiagnosticError(f"{label} failed: {error}") from error

        return wrapped

    return decorate


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise MorphologyDiagnosticError(message)


def _finite_scalar(value: object, label: str, *, minimum: float | None = None) -> float:
    _require(
        isinstance(value, Real) and not isinstance(value, (bool, np.bool_)),
        f"{label} must be a real scalar",
    )
    result = float(value)
    _require(math.isfinite(result), f"{label} must be finite")
    if minimum is not None:
        _require(result >= minimum, f"{label} must be at least {minimum}")
    return result


def _exact_nonnegative_int(value: object, label: str) -> int:
    _require(
        type(value) is int and value >= 0,
        f"{label} must be a non-negative integer",
    )
    return value


def _safe_relative_path(value: object, label: str) -> str:
    _require(type(value) is str and bool(value), f"{label} must be path text")
    path = PurePosixPath(value)
    _require(
        not path.is_absolute()
        and path.as_posix() == value
        and value != "."
        and all(part not in {"", ".", ".."} for part in path.parts),
        f"{label} must be a safe canonical relative path",
    )
    return value


def _sha256(value: object, label: str) -> str:
    _require(
        type(value) is str and _SHA256.fullmatch(value) is not None,
        f"{label} must be a lowercase SHA-256 digest",
    )
    return value


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical_json_bytes(value: object) -> bytes:
    return artifact_helpers.canonical_json_bytes(value)


def _decode_json(payload: object, label: str) -> dict[str, Any]:
    _require(type(payload) is bytes, f"{label} must be bytes")

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            _require(key not in result, f"{label} contains duplicate key {key!r}")
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise MorphologyDiagnosticError(f"{label} contains forbidden constant {value}")

    try:
        decoded = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise MorphologyDiagnosticError(f"{label} is not canonical JSON") from error
    _require(type(decoded) is dict, f"{label} must decode to an object")
    _require(_canonical_json_bytes(decoded) == payload, f"{label} is not canonical JSON")
    return decoded


def source_bindings() -> dict[str, dict[str, str]]:
    """Return the exact source closure used to build morphology artifacts."""
    result: dict[str, dict[str, str]] = {}
    for relative in _SOURCE_PATHS:
        path = REPO_ROOT / relative
        payload = path.read_bytes()
        result[relative] = {"path": relative, "sha256": _sha256_bytes(payload)}
    return result


def _validate_source_binding_map(value: object) -> dict[str, dict[str, str]]:
    _require(
        type(value) is dict and set(value) == set(_SOURCE_PATHS),
        "source binding inventory drifted",
    )
    parsed: dict[str, dict[str, str]] = {}
    for name, raw_binding in value.items():
        relative = _safe_relative_path(name, "source binding key")
        _require(
            type(raw_binding) is dict and set(raw_binding) == {"path", "sha256"},
            f"source binding {relative} schema drifted",
        )
        _require(
            raw_binding["path"] == relative,
            f"source binding {relative} path drifted",
        )
        parsed[relative] = {
            "path": relative,
            "sha256": _sha256(raw_binding["sha256"], f"source binding {relative}"),
        }
    return parsed


def validate_raw_artifact_bindings(
    value: object,
    *,
    observed_committed_time: object,
    cycle: object,
) -> dict[str, dict[str, object]]:
    """Validate exact product, path, checksum, nominal-time, and cycle bindings."""
    observed = _finite_scalar(
        observed_committed_time, "raw binding observed committed time", minimum=0.0
    )
    expected_cycle = _exact_nonnegative_int(cycle, "raw binding cycle")
    _require(
        type(value) is dict and set(value) == set(EXPECTED_RAW_PRODUCTS),
        "raw artifact binding inventory drifted",
    )
    parsed: dict[str, dict[str, object]] = {}
    paths: set[str] = set()
    for product in EXPECTED_RAW_PRODUCTS:
        binding = value[product]
        _require(
            type(binding) is dict
            and set(binding)
            == {
                "product",
                "path",
                "sha256",
                "nominal_slot_time",
                "observed_committed_time",
                "cycle",
            },
            f"{product} raw artifact binding schema drifted",
        )
        _require(binding["product"] == product, f"{product} raw product identity drifted")
        path = _safe_relative_path(binding["path"], f"{product} raw artifact path")
        _require(
            _PRODUCT_PATH_TOKENS[product] in PurePosixPath(path).name,
            f"{product} raw artifact path does not bind its product identity",
        )
        expected_suffix = ".part.vtk" if product == "prtcl_all" else ".bin"
        _require(path.endswith(expected_suffix), f"{product} raw artifact suffix drifted")
        _require(path not in paths, "raw artifact paths must be unique")
        paths.add(path)
        nominal = _finite_scalar(binding["nominal_slot_time"], f"{product} nominal slot")
        bound_observed = _finite_scalar(
            binding["observed_committed_time"], f"{product} observed committed time"
        )
        bound_cycle = _exact_nonnegative_int(binding["cycle"], f"{product} cycle")
        _require(nominal == NOMINAL_SLOT_TIME, f"{product} must bind nominal t=500")
        _require(
            bound_observed == observed,
            f"{product} observed committed time disagrees",
        )
        _require(bound_cycle == expected_cycle, f"{product} cycle disagrees")
        parsed[product] = {
            "product": product,
            "path": path,
            "sha256": _sha256(binding["sha256"], f"{product} raw artifact"),
            "nominal_slot_time": nominal,
            "observed_committed_time": bound_observed,
            "cycle": bound_cycle,
        }
    return parsed


def validate_particle_snapshot_metadata(
    value: object,
    *,
    raw_binding: Mapping[str, object],
    observed_committed_time: object,
    cycle: object,
) -> dict[str, object]:
    """Validate the decoded PVTK identity carried alongside particle arrays."""
    _require(
        type(value) is dict
        and set(value)
        == {
            "source",
            "nominal_slot_time",
            "observed_committed_time",
            "cycle",
            "scalar_fields",
            "vector_fields",
            "storage_layout",
        },
        "particle snapshot metadata schema drifted",
    )
    observed = _finite_scalar(observed_committed_time, "particle observed committed time")
    expected_cycle = _exact_nonnegative_int(cycle, "particle cycle")
    _require(value["source"] == raw_binding["path"], "particle source binding drifted")
    _require(
        _finite_scalar(value["nominal_slot_time"], "particle nominal slot")
        == NOMINAL_SLOT_TIME,
        "particle metadata must bind nominal t=500",
    )
    _require(
        _finite_scalar(value["observed_committed_time"], "particle observed time")
        == observed,
        "particle metadata observed committed time drifted",
    )
    _require(
        _exact_nonnegative_int(value["cycle"], "particle metadata cycle") == expected_cycle,
        "particle metadata cycle drifted",
    )
    _require(
        type(value["scalar_fields"]) is list
        and tuple(value["scalar_fields"]) == EXPECTED_PVTK_SCALARS,
        "particle scalar inventory drifted",
    )
    _require(
        type(value["vector_fields"]) is list
        and tuple(value["vector_fields"]) == EXPECTED_PVTK_VECTORS,
        "particle vector inventory drifted",
    )
    _require(
        value["storage_layout"] == "single_mpi_io_pvtk_file_with_nranks_header",
        "particle storage layout drifted",
    )
    return dict(value)


def _validate_dataset_sources(
    mhd_dataset: object,
    current_datasets: object,
    bindings: Mapping[str, Mapping[str, object]],
    *,
    cycle: int,
) -> None:
    _require(
        isinstance(mhd_dataset, science.output_primitives.AthenaBinaryDataset),
        "mhd_w_bcc requires an AthenaBinaryDataset",
    )
    _require(
        mhd_dataset.source == bindings["mhd_w_bcc"]["path"],
        "mhd_w_bcc decoded source does not match its raw binding",
    )
    _require(mhd_dataset.cycle == cycle, "mhd_w_bcc decoded cycle drifted")
    _require(
        isinstance(current_datasets, Mapping)
        and set(current_datasets) == set(science.CURRENT_PRODUCT_FIELDS),
        "current dataset inventory drifted",
    )
    for product in science.CURRENT_PRODUCT_FIELDS:
        dataset = current_datasets[product]
        _require(
            isinstance(dataset, science.output_primitives.AthenaBinaryDataset),
            f"{product} requires an AthenaBinaryDataset",
        )
        _require(
            dataset.source == bindings[product]["path"],
            f"{product} decoded source does not match its raw binding",
        )
        _require(dataset.cycle == cycle, f"{product} decoded cycle drifted")


def _uniform_spacing(faces: np.ndarray, label: str) -> float:
    differences = np.diff(np.asarray(faces, dtype=np.float64))
    _require(
        differences.ndim == 1
        and differences.size > 0
        and np.all(np.isfinite(differences))
        and np.all(differences > 0.0),
        f"{label} faces must be finite and increasing",
    )
    spacing = float(differences[0])
    _require(
        np.allclose(differences, spacing, rtol=0.0, atol=1.0e-12),
        f"{label} composite grid must be uniform",
    )
    return spacing


def _fixed_faces(lower: float, upper: float, label: str) -> np.ndarray:
    extent = upper - lower
    count = int(round(extent / FIXED_PROFILE_CELL_SIZE_C_OVER_OMEGA_PI))
    _require(
        count > 0
        and math.isclose(
            lower + count * FIXED_PROFILE_CELL_SIZE_C_OVER_OMEGA_PI,
            upper,
            rel_tol=0.0,
            abs_tol=1.0e-9,
        ),
        f"{label} domain is not exactly divisible by the fixed 12 c/omega_pi grid",
    )
    return lower + FIXED_PROFILE_CELL_SIZE_C_OVER_OMEGA_PI * np.arange(
        count + 1, dtype=np.float64
    )


def _restriction_shape(
    state: science.ComposedMHDState,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    x1_faces = _fixed_faces(float(state.x1_faces[0]), float(state.x1_faces[-1]), "x1")
    x2_faces = _fixed_faces(float(state.x2_faces[0]), float(state.x2_faces[-1]), "x2")
    source_dx1 = _uniform_spacing(state.x1_faces, "x1")
    source_dx2 = _uniform_spacing(state.x2_faces, "x2")
    factor_x1 = int(round(FIXED_PROFILE_CELL_SIZE_C_OVER_OMEGA_PI / source_dx1))
    factor_x2 = int(round(FIXED_PROFILE_CELL_SIZE_C_OVER_OMEGA_PI / source_dx2))
    _require(
        factor_x1 > 0
        and factor_x2 > 0
        and math.isclose(
            factor_x1 * source_dx1,
            FIXED_PROFILE_CELL_SIZE_C_OVER_OMEGA_PI,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            factor_x2 * source_dx2,
            FIXED_PROFILE_CELL_SIZE_C_OVER_OMEGA_PI,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "finest composite grid does not exactly subdivide the fixed 12 c/omega_pi grid",
    )
    _require(
        (state.x1_faces.size - 1) == (x1_faces.size - 1) * factor_x1
        and (state.x2_faces.size - 1) == (x2_faces.size - 1) * factor_x2
        and np.allclose(state.x1_faces[::factor_x1], x1_faces, rtol=0.0, atol=1.0e-9)
        and np.allclose(state.x2_faces[::factor_x2], x2_faces, rtol=0.0, atol=1.0e-9),
        "fixed 12 c/omega_pi grid is not aligned with the finest composite grid",
    )
    return x1_faces, x2_faces, factor_x1, factor_x2


def _restrict_mean(
    values_y_x: object,
    state: science.ComposedMHDState,
    *,
    label: str,
) -> np.ndarray:
    values = np.asarray(values_y_x, dtype=np.float64)
    _require(
        values.shape == state.source_levels_y_x.shape
        and np.all(np.isfinite(values)),
        f"{label} must be one finite finest-composite y-x field",
    )
    x1_faces, x2_faces, factor_x1, factor_x2 = _restriction_shape(state)
    ny = x2_faces.size - 1
    nx = x1_faces.size - 1
    restricted = values.reshape(ny, factor_x2, nx, factor_x1).mean(axis=(1, 3))
    _require(
        restricted.shape == (ny, nx) and np.all(np.isfinite(restricted)),
        f"{label} fixed-grid restriction failed",
    )
    return restricted


def _profile_coordinates(
    state: science.ComposedMHDState, snapshot: Mapping[str, object]
) -> dict[str, list[float]]:
    x1_faces, _, _, _ = _restriction_shape(state)
    centers = 0.5 * (x1_faces[:-1] + x1_faces[1:])
    ideal = float(snapshot["mhd"]["ideal_surface_x1_c_over_omega_pi"])
    front = float(snapshot["mhd"]["detected_front"]["x_front_c_over_omega_pi"])
    return {
        "x1_centers_c_over_omega_pi": centers.tolist(),
        "x1_minus_ideal_surface_c_over_omega_pi": (centers - ideal).tolist(),
        "x1_minus_detected_front_c_over_omega_pi": (centers - front).tolist(),
    }


def _profile(values_y_x: object, state: science.ComposedMHDState, *, label: str) -> list[float]:
    return np.mean(_restrict_mean(values_y_x, state, label=label), axis=0).tolist()


def _profile_grid_record(state: science.ComposedMHDState) -> dict[str, object]:
    x1_faces, x2_faces, factor_x1, factor_x2 = _restriction_shape(state)
    return {
        "restriction": (
            "conservative_equal_area_mean_from_finest_composite_then_full_y_average"
        ),
        "fixed_cell_size_c_over_omega_pi": FIXED_PROFILE_CELL_SIZE_C_OVER_OMEGA_PI,
        "x1_faces_c_over_omega_pi": x1_faces.tolist(),
        "x2_faces_c_over_omega_pi": x2_faces.tolist(),
        "finest_composite_cells_per_fixed_x1_cell": factor_x1,
        "finest_composite_cells_per_fixed_x2_cell": factor_x2,
        "target_composite_level": state.target_level,
    }


def _profile_provenance(
    manifest_sha256: str, source_products: Sequence[str]
) -> dict[str, object]:
    _require(
        tuple(source_products)
        == tuple(product for product in EXPECTED_RAW_PRODUCTS if product in source_products),
        "profile source-product order drifted",
    )
    return {
        "morphology_manifest": {
            "path": MANIFEST_MEMBER,
            "sha256": manifest_sha256,
        },
        "source_products": list(source_products),
        "evidence_boundary": _PROFILE_EVIDENCE_BOUNDARY,
    }


def _particle_profile(
    state: science.ComposedMHDState,
    snapshot: Mapping[str, object],
    arrays: Mapping[str, np.ndarray],
    *,
    manifest_sha256: str,
) -> dict[str, object]:
    points = arrays["points"]
    _require(
        np.array_equal(points, points.astype(np.float32).astype(np.float64)),
        "particle points must be exactly representable decoded float32 PVTK data",
    )
    active = arrays["active"]
    energetic = arrays["energetic"]
    for axis, faces in enumerate((state.x1_faces, state.x2_faces, state.x3_faces)):
        selected = points[active, axis]
        _require(
            np.all(selected >= faces[0]) and np.all(selected <= faces[-1]),
            f"active particle escaped retained x{axis + 1} domain",
        )
    x1_faces, _, _, _ = _restriction_shape(state)
    active_x = points[active, 0]
    energetic_x = points[energetic, 0]
    active_count, _ = np.histogram(active_x, bins=x1_faces)
    positive_count, _ = np.histogram(energetic_x, bins=x1_faces)
    macro_weight, _ = np.histogram(
        points[active, 0], bins=x1_faces, weights=arrays["weights"][active]
    )
    kinetic_weights = (
        science.PARTICLE_MACRO_MASS
        * arrays["weights"][energetic]
        * arrays["specific_kinetic_energy"][energetic]
    )
    kinetic_energy, _ = np.histogram(energetic_x, bins=x1_faces, weights=kinetic_weights)
    total_energy = float(np.sum(kinetic_weights))
    _require(
        math.isfinite(total_energy) and total_energy > 0.0,
        "active positive-weight Q011 CR kinetic energy must be positive",
    )
    energy_fraction = kinetic_energy / total_energy
    probability_density = energy_fraction / np.diff(x1_faces)
    _require(
        int(np.sum(active_count)) == int(np.count_nonzero(active))
        and int(np.sum(positive_count)) == int(np.count_nonzero(energetic)),
        "particle fixed-grid count histogram failed closure",
    )
    macro_weight_total = float(np.sum(arrays["weights"][active]))
    _require(
        math.isclose(float(np.sum(macro_weight)), macro_weight_total, rel_tol=1.0e-13)
        and math.isclose(float(np.sum(kinetic_energy)), total_energy, rel_tol=1.0e-13)
        and math.isclose(float(np.sum(energy_fraction)), 1.0, rel_tol=1.0e-13)
        and math.isclose(
            float(np.sum(probability_density * np.diff(x1_faces))),
            1.0,
            rel_tol=1.0e-13,
        ),
        "particle fixed-grid energy-weighted profile failed normalization closure",
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_t500_particle_energy_weighted_spatial_profile_v1",
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "nominal_slot_time": NOMINAL_SLOT_TIME,
        "observed_committed_time": snapshot["observed_committed_time"],
        "provenance": _profile_provenance(manifest_sha256, ("prtcl_all",)),
        "selection": {
            "source_rule": "cr_source == 1",
            "birth_time_rule": "birth_time >= 45",
            "energy_weight_population": "selected_particles_with_positive_macro_weight",
            "energy_weight": (
                "particle_macro_mass_times_macro_weight_times_reconstructed_"
                "relativistic_specific_kinetic_energy"
            ),
        },
        "profile_grid": _profile_grid_record(state),
        "coordinates": _profile_coordinates(state, snapshot),
        "profiles": {
            "active_particle_count": active_count.astype(np.int64).tolist(),
            "active_positive_weight_particle_count": positive_count.astype(np.int64).tolist(),
            "active_macro_weight": macro_weight.tolist(),
            "active_cr_kinetic_energy": kinetic_energy.tolist(),
            "active_cr_kinetic_energy_fraction": energy_fraction.tolist(),
            "energy_weighted_probability_density_per_c_over_omega_pi": (
                probability_density.tolist()
            ),
        },
        "normalization_closure": {
            "active_particle_count": int(np.sum(active_count)),
            "active_positive_weight_particle_count": int(np.sum(positive_count)),
            "active_macro_weight": macro_weight_total,
            "active_cr_kinetic_energy": total_energy,
            "energy_fraction_sum": float(np.sum(energy_fraction)),
            "probability_density_integral": float(
                np.sum(probability_density * np.diff(x1_faces))
            ),
        },
    }


def _current_profile(
    mhd_dataset: science.output_primitives.AthenaBinaryDataset,
    current_datasets: Mapping[str, science.output_primitives.AthenaBinaryDataset],
    state: science.ComposedMHDState,
    snapshot: Mapping[str, object],
    *,
    manifest_sha256: str,
) -> dict[str, object]:
    fields = science._compose_matched_current_fields(mhd_dataset, state, current_datasets)
    lab = np.stack([fields[f"prtcl_j{component}"] for component in ("x", "y", "z")])
    velocity = np.stack(
        [state.fields_y_x[f"vel{component}"] for component in ("x", "y", "z")]
    )
    gas = lab - fields["prtcl_rho"][None, :, :] * velocity
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_t500_deposited_current_spatial_profile_v1",
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "nominal_slot_time": NOMINAL_SLOT_TIME,
        "observed_committed_time": snapshot["observed_committed_time"],
        "provenance": _profile_provenance(
            manifest_sha256,
            ("mhd_w_bcc", "prtcl_rho", "prtcl_jx", "prtcl_jy", "prtcl_jz", "mhd_j2"),
        ),
        "profile_grid": _profile_grid_record(state),
        "coordinates": _profile_coordinates(state, snapshot),
        "source_representation": snapshot["cr_current"]["source_representation"],
        "frame_transform": snapshot["cr_current"]["frame_transform"],
        "profiles": {
            "deposited_rho_cr_over_c": _profile(
                fields["prtcl_rho"], state, label="deposited rho_CR/c"
            ),
            "deposited_j_cr_over_c_lab_x": _profile(lab[0], state, label="lab current x"),
            "deposited_j_cr_over_c_lab_y": _profile(lab[1], state, label="lab current y"),
            "deposited_j_cr_over_c_lab_z": _profile(lab[2], state, label="lab current z"),
            "deposited_j_cr_over_c_lab_magnitude": _profile(
                np.sqrt(np.sum(lab * lab, axis=0)), state, label="lab current magnitude"
            ),
            "deposited_j_cr_over_c_gas_x": _profile(gas[0], state, label="gas current x"),
            "deposited_j_cr_over_c_gas_y": _profile(gas[1], state, label="gas current y"),
            "deposited_j_cr_over_c_gas_z": _profile(gas[2], state, label="gas current z"),
            "deposited_j_cr_over_c_gas_magnitude": _profile(
                np.sqrt(np.sum(gas * gas, axis=0)), state, label="gas current magnitude"
            ),
            "mhd_current_squared": _profile(
                fields["mhd_j2"], state, label="MHD current squared"
            ),
        },
    }


def _magnetic_profile(
    state: science.ComposedMHDState,
    snapshot: Mapping[str, object],
    *,
    manifest_sha256: str,
) -> dict[str, object]:
    magnetic = np.stack(
        [state.fields_y_x[f"bcc{component}"] for component in (1, 2, 3)]
    )
    magnitude = np.sqrt(np.sum(magnetic * magnetic, axis=0))
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_t500_magnetic_field_spatial_profile_v1",
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "nominal_slot_time": NOMINAL_SLOT_TIME,
        "observed_committed_time": snapshot["observed_committed_time"],
        "provenance": _profile_provenance(manifest_sha256, ("mhd_w_bcc",)),
        "profile_grid": _profile_grid_record(state),
        "coordinates": _profile_coordinates(state, snapshot),
        "profiles": {
            "bcc1": _profile(magnetic[0], state, label="bcc1"),
            "bcc2": _profile(magnetic[1], state, label="bcc2"),
            "bcc3": _profile(magnetic[2], state, label="bcc3"),
            "bmag": _profile(magnitude, state, label="magnetic magnitude"),
            "bmag_over_b0": _profile(
                magnitude / science.REFERENCE_B0, state, label="magnetic amplification"
            ),
            "gas_density_context": _profile(
                state.fields_y_x["dens"], state, label="gas density context"
            ),
        },
        "t500_upstream_magnetic_amplification_context": snapshot["mhd"][
            "upstream_magnetic_amplification"
        ],
    }


def _refinement_profile(
    state: science.ComposedMHDState,
    snapshot: Mapping[str, object],
    *,
    manifest_sha256: str,
) -> dict[str, object]:
    levels = np.unique(state.source_levels_y_x)
    _require(
        levels.size > 0
        and np.issubdtype(levels.dtype, np.integer)
        and np.all(levels >= 0)
        and np.all(levels <= state.target_level),
        "source refinement levels are invalid",
    )
    fractions = {
        f"source_level_{int(level)}_area_fraction": _profile(
            (state.source_levels_y_x == level).astype(np.float64),
            state,
            label=f"source level {int(level)} area fraction",
        )
        for level in levels
    }
    fraction_sum = np.sum(np.asarray(list(fractions.values()), dtype=np.float64), axis=0)
    _require(
        np.allclose(fraction_sum, 1.0, rtol=0.0, atol=1.0e-15),
        "refinement area-fraction profiles failed closure",
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_t500_refinement_spatial_profile_v1",
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "nominal_slot_time": NOMINAL_SLOT_TIME,
        "observed_committed_time": snapshot["observed_committed_time"],
        "provenance": _profile_provenance(manifest_sha256, ("mhd_w_bcc",)),
        "profile_grid": _profile_grid_record(state),
        "coordinates": _profile_coordinates(state, snapshot),
        "encoding": "full_y_area_fraction_by_physical_source_refinement_level",
        "profiles": fractions,
        "area_fraction_sum_by_x1": fraction_sum.tolist(),
    }


def _refinement_overlay(
    mhd_dataset: science.output_primitives.AthenaBinaryDataset,
    state: science.ComposedMHDState,
    snapshot: Mapping[str, object],
    *,
    manifest_sha256: str,
) -> dict[str, object]:
    rows = []
    for block in mhd_dataset.blocks:
        geometry = tuple(float(value) for value in block.geometry)
        _require(
            len(geometry) == 6 and all(math.isfinite(value) for value in geometry),
            "refinement overlay block geometry is invalid",
        )
        rows.append(
            {
                "source_level": int(block.level),
                "logical_location": [int(value) for value in block.logical_location],
                "index_bounds": [int(value) for value in block.index_bounds],
                "x1_bounds_c_over_omega_pi": [geometry[0], geometry[1]],
                "x2_bounds_c_over_omega_pi": [geometry[2], geometry[3]],
            }
        )
    rows.sort(
        key=lambda row: (
            row["source_level"],
            row["logical_location"],
            row["x1_bounds_c_over_omega_pi"],
            row["x2_bounds_c_over_omega_pi"],
        )
    )
    _require(bool(rows), "refinement overlay requires at least one leaf block")
    levels, counts = np.unique(state.source_levels_y_x, return_counts=True)
    block_levels, block_counts = np.unique(
        np.asarray([row["source_level"] for row in rows], dtype=np.int64),
        return_counts=True,
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_t500_refinement_leaf_block_overlay_v1",
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "nominal_slot_time": NOMINAL_SLOT_TIME,
        "observed_committed_time": snapshot["observed_committed_time"],
        "provenance": _profile_provenance(manifest_sha256, ("mhd_w_bcc",)),
        "rendering_contract": {
            "overlay": "draw_each_leaf_block_x1_x2_rectangle",
            "color_key": "source_level",
            "physical_refinement_level_not_target_composite_level": True,
        },
        "target_composite_level": state.target_level,
        "leaf_block_rectangles": rows,
        "leaf_block_count_by_source_level": [
            {"source_level": int(level), "leaf_block_count": int(count)}
            for level, count in zip(block_levels, block_counts)
        ],
        "finest_composite_cell_count_by_physical_source_level": [
            {"source_level": int(level), "cell_count": int(count)}
            for level, count in zip(levels, counts)
        ],
    }


def _manifest_record(
    *,
    attempt_id: str,
    observed_committed_time: float,
    cycle: int,
    raw_artifact_bindings: Mapping[str, Mapping[str, object]],
    source_binding_map: Mapping[str, Mapping[str, str]],
) -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_t500_morphology_publication_manifest_v1",
        "successor_id": SUCCESSOR_ID,
        "artifact_role": _MANIFEST_ARTIFACT_ROLE,
        "qualification_effect": QUALIFICATION_EFFECT,
        "authorization": dict(AUTHORIZATION),
        "attempt_id": attempt_id,
        "nominal_slot_time": NOMINAL_SLOT_TIME,
        "observed_committed_time": observed_committed_time,
        "cycle": cycle,
        "raw_artifact_bindings": dict(raw_artifact_bindings),
        "source_bindings": dict(source_binding_map),
        "derived_member_paths": list(DERIVED_MEMBERS[1:]),
        "profile_contract": dict(_MANIFEST_PROFILE_CONTRACT),
        "evidence_boundary": _MANIFEST_EVIDENCE_BOUNDARY,
    }


@_public_contract("t=500 morphology artifact construction")
def build_morphology_artifacts(
    *,
    attempt_id: object,
    mhd_dataset: science.output_primitives.AthenaBinaryDataset,
    current_datasets: Mapping[str, science.output_primitives.AthenaBinaryDataset],
    raw_artifact_bindings: object,
    particle_snapshot_metadata: object,
    observed_committed_time: object,
    cycle: object,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
    target_level: int | None = None,
) -> dict[str, bytes]:
    """Build and self-verify the deterministic publication-bound artifact set."""
    _require(
        type(attempt_id) is str and _ATTEMPT_ID.fullmatch(attempt_id) is not None,
        "attempt id is invalid",
    )
    observed = _finite_scalar(
        observed_committed_time, "morphology observed committed time", minimum=0.0
    )
    exact_cycle = _exact_nonnegative_int(cycle, "morphology cycle")
    bindings = validate_raw_artifact_bindings(
        raw_artifact_bindings,
        observed_committed_time=observed,
        cycle=exact_cycle,
    )
    validate_particle_snapshot_metadata(
        particle_snapshot_metadata,
        raw_binding=bindings["prtcl_all"],
        observed_committed_time=observed,
        cycle=exact_cycle,
    )
    _validate_dataset_sources(mhd_dataset, current_datasets, bindings, cycle=exact_cycle)
    science.validate_successor_deck()
    snapshot = science.reduce_production_science_snapshot(
        mhd_dataset,
        current_datasets,
        nominal_slot_time=NOMINAL_SLOT_TIME,
        observed_committed_time=observed,
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
        target_level=target_level,
    )
    _require(
        snapshot["nominal_slot_time"] == NOMINAL_SLOT_TIME
        and snapshot["observed_committed_time"] == observed,
        "production-science snapshot identity drifted",
    )
    state = science.compose_full_mhd_state(
        mhd_dataset,
        nominal_slot_time=NOMINAL_SLOT_TIME,
        observed_committed_time=observed,
        target_level=target_level,
    )
    _restriction_shape(state)
    arrays = science._decoded_particle_arrays(
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )
    source_binding_map = _validate_source_binding_map(source_bindings())
    manifest = _manifest_record(
        attempt_id=attempt_id,
        observed_committed_time=observed,
        cycle=exact_cycle,
        raw_artifact_bindings=bindings,
        source_binding_map=source_binding_map,
    )
    manifest_payload = _canonical_json_bytes(manifest)
    manifest_sha256 = _sha256_bytes(manifest_payload)
    members = {
        MANIFEST_MEMBER: manifest_payload,
        PARTICLE_PROFILE_MEMBER: _canonical_json_bytes(
            _particle_profile(state, snapshot, arrays, manifest_sha256=manifest_sha256)
        ),
        CURRENT_PROFILE_MEMBER: _canonical_json_bytes(
            _current_profile(
                mhd_dataset,
                current_datasets,
                state,
                snapshot,
                manifest_sha256=manifest_sha256,
            )
        ),
        MAGNETIC_PROFILE_MEMBER: _canonical_json_bytes(
            _magnetic_profile(state, snapshot, manifest_sha256=manifest_sha256)
        ),
        REFINEMENT_PROFILE_MEMBER: _canonical_json_bytes(
            _refinement_profile(state, snapshot, manifest_sha256=manifest_sha256)
        ),
        REFINEMENT_OVERLAY_MEMBER: _canonical_json_bytes(
            _refinement_overlay(
                mhd_dataset,
                state,
                snapshot,
                manifest_sha256=manifest_sha256,
            )
        ),
    }
    inventory = {
        "schema_version": SCHEMA_VERSION,
        "record_type": "q011_section54_t500_morphology_artifact_inventory_v1",
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "authorization": dict(AUTHORIZATION),
        "morphology_manifest": {
            "path": MANIFEST_MEMBER,
            "sha256": manifest_sha256,
        },
        "members": [
            {
                "path": path,
                "bytes": len(payload),
                "sha256": _sha256_bytes(payload),
            }
            for path, payload in sorted(members.items())
        ],
    }
    members[INVENTORY_MEMBER] = _canonical_json_bytes(inventory)
    result = {path: members[path] for path in ALL_MEMBERS}
    validate_morphology_artifacts(result)
    return result


@_public_contract("t=500 morphology artifact validation")
def validate_morphology_artifacts(value: object) -> dict[str, object]:
    """Validate exact membership, canonical JSON, provenance, and checksums."""
    _require(
        type(value) is dict and set(value) == set(ALL_MEMBERS),
        "morphology artifact member inventory drifted",
    )
    records = {path: _decode_json(value[path], path) for path in ALL_MEMBERS}
    manifest = records[MANIFEST_MEMBER]
    _require(
        set(manifest)
        == {
            "schema_version",
            "record_type",
            "successor_id",
            "artifact_role",
            "qualification_effect",
            "authorization",
            "attempt_id",
            "nominal_slot_time",
            "observed_committed_time",
            "cycle",
            "raw_artifact_bindings",
            "source_bindings",
            "derived_member_paths",
            "profile_contract",
            "evidence_boundary",
        },
        "morphology manifest schema drifted",
    )
    _require(
        manifest["schema_version"] == SCHEMA_VERSION
        and type(manifest["schema_version"]) is int
        and manifest["record_type"]
        == "q011_section54_t500_morphology_publication_manifest_v1"
        and manifest["successor_id"] == SUCCESSOR_ID
        and manifest["artifact_role"] == _MANIFEST_ARTIFACT_ROLE
        and manifest["qualification_effect"] == QUALIFICATION_EFFECT
        and manifest["authorization"] == dict(AUTHORIZATION)
        and manifest["nominal_slot_time"] == NOMINAL_SLOT_TIME
        and manifest["derived_member_paths"] == list(DERIVED_MEMBERS[1:])
        and manifest["profile_contract"] == _MANIFEST_PROFILE_CONTRACT
        and manifest["evidence_boundary"] == _MANIFEST_EVIDENCE_BOUNDARY,
        "morphology manifest identity, contract, evidence boundary, or authority drifted",
    )
    _require(
        type(manifest["attempt_id"]) is str
        and _ATTEMPT_ID.fullmatch(manifest["attempt_id"]) is not None,
        "morphology manifest attempt id drifted",
    )
    observed = _finite_scalar(manifest["observed_committed_time"], "manifest observed time")
    cycle = _exact_nonnegative_int(manifest["cycle"], "manifest cycle")
    validate_raw_artifact_bindings(
        manifest["raw_artifact_bindings"],
        observed_committed_time=observed,
        cycle=cycle,
    )
    _validate_source_binding_map(manifest["source_bindings"])
    manifest_sha256 = _sha256_bytes(value[MANIFEST_MEMBER])
    expected_types = {
        PARTICLE_PROFILE_MEMBER: "q011_section54_t500_particle_energy_weighted_spatial_profile_v1",
        CURRENT_PROFILE_MEMBER: "q011_section54_t500_deposited_current_spatial_profile_v1",
        MAGNETIC_PROFILE_MEMBER: "q011_section54_t500_magnetic_field_spatial_profile_v1",
        REFINEMENT_PROFILE_MEMBER: "q011_section54_t500_refinement_spatial_profile_v1",
        REFINEMENT_OVERLAY_MEMBER: "q011_section54_t500_refinement_leaf_block_overlay_v1",
    }
    expected_products = {
        PARTICLE_PROFILE_MEMBER: ["prtcl_all"],
        CURRENT_PROFILE_MEMBER: [
            "mhd_w_bcc",
            "prtcl_rho",
            "prtcl_jx",
            "prtcl_jy",
            "prtcl_jz",
            "mhd_j2",
        ],
        MAGNETIC_PROFILE_MEMBER: ["mhd_w_bcc"],
        REFINEMENT_PROFILE_MEMBER: ["mhd_w_bcc"],
        REFINEMENT_OVERLAY_MEMBER: ["mhd_w_bcc"],
    }
    common_profile_keys = {
        "schema_version",
        "record_type",
        "successor_id",
        "qualification_effect",
        "nominal_slot_time",
        "observed_committed_time",
        "provenance",
    }
    expected_keys = {
        PARTICLE_PROFILE_MEMBER: common_profile_keys
        | {
            "selection",
            "profile_grid",
            "coordinates",
            "profiles",
            "normalization_closure",
        },
        CURRENT_PROFILE_MEMBER: common_profile_keys
        | {
            "profile_grid",
            "coordinates",
            "source_representation",
            "frame_transform",
            "profiles",
        },
        MAGNETIC_PROFILE_MEMBER: common_profile_keys
        | {
            "profile_grid",
            "coordinates",
            "profiles",
            "t500_upstream_magnetic_amplification_context",
        },
        REFINEMENT_PROFILE_MEMBER: common_profile_keys
        | {
            "profile_grid",
            "coordinates",
            "encoding",
            "profiles",
            "area_fraction_sum_by_x1",
        },
        REFINEMENT_OVERLAY_MEMBER: common_profile_keys
        | {
            "rendering_contract",
            "target_composite_level",
            "leaf_block_rectangles",
            "leaf_block_count_by_source_level",
            "finest_composite_cell_count_by_physical_source_level",
        },
    }
    for path, record_type in expected_types.items():
        record = records[path]
        _require(
            set(record) == expected_keys[path]
            and record.get("schema_version") == SCHEMA_VERSION
            and type(record["schema_version"]) is int
            and record.get("record_type") == record_type
            and record.get("successor_id") == SUCCESSOR_ID
            and record.get("qualification_effect") == QUALIFICATION_EFFECT
            and record.get("nominal_slot_time") == NOMINAL_SLOT_TIME
            and record.get("observed_committed_time") == observed,
            f"{path} identity drifted",
        )
        provenance = record.get("provenance")
        _require(
            type(provenance) is dict
            and set(provenance)
            == {"morphology_manifest", "source_products", "evidence_boundary"}
            and provenance.get("morphology_manifest")
            == {"path": MANIFEST_MEMBER, "sha256": manifest_sha256},
            f"{path} morphology-manifest provenance drifted",
        )
        _require(
            provenance["source_products"] == expected_products[path]
            and provenance["evidence_boundary"] == _PROFILE_EVIDENCE_BOUNDARY,
            f"{path} source-product provenance or evidence boundary drifted",
        )
    inventory = records[INVENTORY_MEMBER]
    _require(
        set(inventory)
        == {
            "schema_version",
            "record_type",
            "successor_id",
            "qualification_effect",
            "authorization",
            "morphology_manifest",
            "members",
        }
        and inventory["schema_version"] == SCHEMA_VERSION
        and type(inventory["schema_version"]) is int
        and inventory["record_type"]
        == "q011_section54_t500_morphology_artifact_inventory_v1"
        and inventory["successor_id"] == SUCCESSOR_ID
        and inventory["qualification_effect"] == QUALIFICATION_EFFECT
        and inventory["authorization"] == dict(AUTHORIZATION)
        and inventory["morphology_manifest"]
        == {"path": MANIFEST_MEMBER, "sha256": manifest_sha256},
        "morphology inventory identity or authority drifted",
    )
    expected_inventory = [
        {
            "path": path,
            "bytes": len(value[path]),
            "sha256": _sha256_bytes(value[path]),
        }
        for path in sorted(DERIVED_MEMBERS)
    ]
    _require(
        inventory["members"] == expected_inventory,
        "morphology artifact inventory checksum or membership drifted",
    )
    return manifest


__all__ = [
    "ALL_MEMBERS",
    "AUTHORIZATION",
    "CURRENT_PROFILE_MEMBER",
    "EXPECTED_PVTK_SCALARS",
    "EXPECTED_PVTK_VECTORS",
    "EXPECTED_RAW_PRODUCTS",
    "FIXED_PROFILE_CELL_SIZE_C_OVER_OMEGA_PI",
    "INVENTORY_MEMBER",
    "MAGNETIC_PROFILE_MEMBER",
    "MANIFEST_MEMBER",
    "MorphologyDiagnosticError",
    "NOMINAL_SLOT_TIME",
    "PARTICLE_PROFILE_MEMBER",
    "QUALIFICATION_EFFECT",
    "REFINEMENT_OVERLAY_MEMBER",
    "REFINEMENT_PROFILE_MEMBER",
    "SUCCESSOR_ID",
    "build_morphology_artifacts",
    "source_bindings",
    "validate_morphology_artifacts",
    "validate_particle_snapshot_metadata",
    "validate_raw_artifact_bindings",
]
