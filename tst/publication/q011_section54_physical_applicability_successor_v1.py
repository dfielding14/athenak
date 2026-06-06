#!/usr/bin/env python3
"""Fail-closed source-local physical-applicability diagnostics for Q011.

This additive successor verifies decoded Q011 production-science mesh,
deposited-moment, and particle products against opened immutable raw bytes. It
computes snapshot applicability diagnostics and validates byte-bound future
runtime time/escape records. It does not launch work, mutate policy, authorize
qualifying-output inspection, or close any scientific claim.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import wraps
import hashlib
import json
import math
from numbers import Real
import os
from pathlib import Path, PurePosixPath
import re
import stat
from types import MappingProxyType
from typing import Any

import numpy as np

if __package__:
    from . import analyze_q011_section54_outputs as output_primitives
    from . import pvtk_particles
    from . import q011_section54_production_science_successor_v1 as science
else:
    import analyze_q011_section54_outputs as output_primitives
    import pvtk_particles
    import q011_section54_production_science_successor_v1 as science


SCHEMA_VERSION = 1
SUCCESSOR_ID = "q011_section54_physical_applicability_successor_v1"
SNAPSHOT_RECORD_TYPE = "q011_section54_physical_applicability_snapshot_v1"
HISTORY_RECORD_TYPE = "q011_section54_physical_applicability_history_v1"
RUNTIME_RECORD_TYPE = "q011_section54_runtime_time_escape_applicability_v1"
NORMALIZATION_RECORD_TYPE = "q011_section54_bound_normalization_evidence_v1"
SOURCE_MANIFEST_RECORD_TYPE = "q011_section54_source_manifest_v1"
SNAPSHOT_PROVENANCE_RECORD_TYPE = "q011_section54_snapshot_provenance_manifest_v1"
CYCLE_TELEMETRY_RECORD_TYPE = "q011_section54_bound_cycle_telemetry_v1"
SLOPE_CUTOFF_ESCAPE_RECORD_TYPE = "q011_section54_bound_slope_cutoff_escape_evidence_v1"
QUALIFICATION_EFFECT = (
    "source_local_non_authorizing_diagnostic_and_gate_only_no_launch_no_policy_"
    "mutation_no_qualifying_output_inspection_no_claim_closure"
)

AUTHORIZATION: Mapping[str, bool] = MappingProxyType(
    {
        "launch_authorized": False,
        "policy_mutation_authorized": False,
        "qualifying_output_inspection_authorized": False,
        "claim_closure_authorized": False,
    }
)

# These values and representation labels must be supplied exactly.  The
# successor intentionally does not infer normalization from deck arithmetic.
EXACT_NORMALIZATION: Mapping[str, object] = MappingProxyType(
    {
        "selected_species_count": 1,
        "selected_species_index": 0,
        "selected_species_charge_sign": 1,
        "selected_species_abs_q_over_m": 1.0,
        "reference_rho0": 1.0,
        "reference_b0": 1.0,
        "particle_light_speed": 10000.0,
        "gas_density_representation": "ion_mass_density_rho_g",
        "prtcl_rho_representation": "deposited_rho_CR_over_c_named_rho_q",
        "prtcl_j_representation": "deposited_J_CR_over_c_named_J_q",
        "alfven_speed_formula": "abs_B_over_sqrt_rho_g",
        "gas_frame_current_formula": "J_q_minus_rho_q_times_v_g",
        "particle_momentum_formula": "gamma_of_v_times_v",
        "particle_gyroradius_formula": "abs_p_over_abs_q_over_m_times_abs_B_local",
    }
)

# Bai et al. require R << 1, Lambda << 1 for the no-Hall approximation, and
# MHD-PIC scales much larger than d_i, but do not prescribe these numeric
# thresholds.  They are conservative AthenaK-selected applicability bounds.
R_MAXIMUM = 0.01
LAMBDA_MAXIMUM = 0.1
S_DELTA_MINIMUM = 1.0
LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM = 10.0
SUB_10DI_POWER_FRACTION_MAXIMUM = 0.05
DELTA_B_RMS_OVER_B0_MINIMUM = 0.1

# Sun & Bai state that the transverse domain should contain several
# high-energy gyroradii but do not prescribe these numeric thresholds.  They
# are conservative AthenaK-selected containment bounds.
RG_MAXIMUM_OVER_LY_MAXIMUM = 1.0 / 8.0
RG_HIGH_ENERGY_MAXIMUM_OVER_LY_MAXIMUM = 1.0 / 8.0
RG_Q999_OVER_LY_MAXIMUM = 1.0 / 8.0
RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM = 1.0e-3
ESCAPED_PARTICLE_COUNT_FRACTION_MAXIMUM = 1.0e-3
ESCAPED_MACRO_WEIGHT_FRACTION_MAXIMUM = 1.0e-3
ESCAPED_KINETIC_ENERGY_FRACTION_MAXIMUM = 1.0e-3
ESCAPE_LEDGER_RELATIVE_RESIDUAL_MAXIMUM = 1.0e-12
PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT = "d614e5c84aad3a541dc96af68ef1178dabc66f71"
TRUSTED_Q011_RUNTIME_SOURCE_COMMIT = PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT
PS_ESCAPE_LEDGER_SCHEMA = 1
PS_CR_LEDGER_SCHEMA = 3
PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE = 2
SLOPE_CUTOFF_BIN_ESCAPE_FRACTION_MAXIMUM = 1.0e-3
SLOPE_CUTOFF_TAIL_ESCAPE_FRACTION_MAXIMUM = 1.0e-3
MINIMUM_SLOPE_CUTOFF_HIGH_ENERGY_BINS = 4
INNER_X1_ESCAPE_REASON = "inner_x1_escape_forbidden"
OUTER_X1_ESCAPE_REASON = "outer_x1_physical_boundary_escape"
REQUIRED_PS_ESCAPE_CHECKPOINT_NOMINAL_TIMES = tuple(
    float(value) for value in range(100, 1201, 100)
)
CLAIM_SPECIFIC_ESCAPE_BOUNDS: Mapping[str, Mapping[str, float | bool]] = (
    MappingProxyType(
        {
            "Emax_claim": MappingProxyType(
                {
                    "particle_count_fraction_maximum": 1.0e-3,
                    "macro_weight_fraction_maximum": 1.0e-3,
                    "kinetic_energy_fraction_maximum": 1.0e-3,
                    "requires_escaped_maximum_inclusion": True,
                }
            ),
            "high_energy_slope_or_cutoff_claim": MappingProxyType(
                {
                    "particle_count_fraction_maximum": 1.0e-3,
                    "macro_weight_fraction_maximum": 1.0e-3,
                    "kinetic_energy_fraction_maximum": 1.0e-3,
                    "requires_escaped_maximum_inclusion": False,
                }
            ),
            "acceleration_rate_claim": MappingProxyType(
                {
                    "particle_count_fraction_maximum": 1.0e-3,
                    "macro_weight_fraction_maximum": 1.0e-3,
                    "kinetic_energy_fraction_maximum": 1.0e-3,
                    "requires_escaped_maximum_inclusion": True,
                }
            ),
            "acceleration_efficiency_claim": MappingProxyType(
                {
                    "particle_count_fraction_maximum": 1.0e-3,
                    "macro_weight_fraction_maximum": 1.0e-3,
                    "kinetic_energy_fraction_maximum": 1.0e-3,
                    "requires_escaped_maximum_inclusion": False,
                }
            ),
        }
    )
)

STARTUP_REMOVAL_TIME = 45.0
EXPECTED_TERMINAL_TIME = 1200.0
MAXIMUM_RUNTIME_CYCLE_SPAN = 1.0
MINIMUM_RUNTIME_CYCLE_INVENTORY_COUNT = 1000
PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES = 1000
SHOCK_TRANSITION_HALF_WIDTH = 120.0
DETECTED_FRONT_REGIONS: Mapping[str, tuple[float, float] | None] = MappingProxyType(
    {
        "full_domain": None,
        "detected_front_downstream": (-1200.0, -120.0),
        "detected_front_precursor": (120.0, 1200.0),
        "detected_front_far_upstream": (1200.0, 2400.0),
    }
)
DI_MAGNETIC_SPECTRUM_REGION = "detected_front_precursor"
REQUIRED_RAW_PRODUCTS = (
    "mhd_w_bcc",
    "prtcl_rho",
    "prtcl_jx",
    "prtcl_jy",
    "prtcl_jz",
    "mhd_j2",
    "prtcl_all",
)
CELL_MAP_NAMES = (
    "R",
    "Lambda",
    "d_i",
    "S_delta",
    "actual_leaf_dx1",
    "actual_leaf_dx2",
    "gas_frame_current_magnitude",
)
_SHA256 = re.compile(r"[0-9a-f]{64}")
_SOURCE_COMMIT = re.compile(r"[0-9a-f]{40}")
_ATTEMPT_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:-]{0,255}")
_PVTK_EXECUTION_PATTERN = re.compile(
    rb"^# vtk DataFile Version 2\.0\n"
    rb"# AthenaK particle data at time= ([^ \n]+)  "
    rb"nranks= (0|[1-9][0-9]*)  cycle=(0|[1-9][0-9]*)  variables=([^\n]+)\n"
)
_PVTK_SCALARS = {
    "gid",
    "ptag",
    "species",
    "cr_source",
    "macro_weight",
    "birth_time",
    "deltaf_f0",
    "deltaf_weight",
}

CLAIM_REJECTION_RULES: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        "Q011-APP-NORM": (
            "all_physical_MHD_PIC_Bell_shock_and_DSA_claims",
        ),
        "Q011-APP-R": (
            "all_physical_MHD_PIC_Bell_shock_and_DSA_claims",
        ),
        "Q011-APP-LAMBDA": (
            "Hall_negligible_claim",
            "Bell_mechanism_claim",
            "physical_magnetic_amplification_claim",
            "physical_DSA_scattering_claim",
            "Emax_claim",
            "high_energy_slope_or_cutoff_claim",
            "acceleration_rate_claim",
            "acceleration_efficiency_claim",
        ),
        "Q011-APP-DI": (
            "physical_precursor_turbulence_claim",
            "physical_scattering_interpretation",
            "Emax_claim",
            "high_energy_slope_or_cutoff_claim",
            "acceleration_rate_claim",
            "acceleration_efficiency_claim",
        ),
        "Q011-APP-RG": (
            "Emax_claim",
            "high_energy_slope_or_cutoff_claim",
            "acceleration_rate_claim",
            "acceleration_efficiency_claim",
        ),
        "Q011-APP-TIME": (
            "all_history_and_global_applicability_claims",
            "Emax_claim",
            "high_energy_slope_or_cutoff_claim",
            "acceleration_rate_claim",
            "acceleration_efficiency_claim",
        ),
    }
)
PERMANENT_CLAIM_EXCLUSIONS = (
    "microscopic_shock_structure_claim",
    "self_consistent_injection_claim",
)


class PhysicalApplicabilityError(ValueError):
    """Raised when Q011 physical applicability cannot be established."""


_UNDERLYING_EXCEPTIONS = (
    science.ProductionScienceError,
    KeyError,
    IndexError,
    AttributeError,
    TypeError,
    ValueError,
    OverflowError,
    FloatingPointError,
    RecursionError,
    OSError,
)


def _public_contract(label: str):
    def decorate(function):
        @wraps(function)
        def wrapped(*args: object, **kwargs: object):
            try:
                return function(*args, **kwargs)
            except PhysicalApplicabilityError:
                raise
            except _UNDERLYING_EXCEPTIONS as error:
                raise PhysicalApplicabilityError(f"{label} failed: {error}") from error

        return wrapped

    return decorate


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PhysicalApplicabilityError(message)


def _exact_keys(value: object, expected: set[str], label: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), f"{label} must be a mapping")
    _require(set(value) == expected, f"{label} keys drifted")
    return value


def _finite_scalar(
    value: object, label: str, *, minimum: float | None = None
) -> float:
    _require(
        isinstance(value, Real) and not isinstance(value, (bool, np.bool_)),
        f"{label} must be a real scalar",
    )
    result = float(value)
    _require(math.isfinite(result), f"{label} must be finite")
    if minimum is not None:
        _require(result >= minimum, f"{label} must be at least {minimum}")
    return result


def _nonnegative_int(value: object, label: str) -> int:
    _require(type(value) is int and value >= 0, f"{label} must be a non-negative integer")
    return value


def _strict_bool(value: object, label: str) -> bool:
    _require(type(value) is bool, f"{label} must be a boolean")
    return value


def _exact_numeric(value: object, expected: float, label: str) -> float:
    decoded = _finite_scalar(value, label)
    _require(decoded == expected, f"{label} drifted")
    return decoded


def _validate_authorization(value: object, label: str) -> dict[str, bool]:
    authorization = _exact_keys(value, set(AUTHORIZATION), label)
    for key, expected in AUTHORIZATION.items():
        _require(
            _strict_bool(authorization[key], f"{label} {key}") is expected,
            f"{label} {key} drifted",
        )
    return dict(AUTHORIZATION)


def _canonical_json_bytes(value: object) -> bytes:
    try:
        return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
            "utf-8"
        )
    except (TypeError, ValueError, OverflowError) as error:
        raise PhysicalApplicabilityError("value is not canonical-JSON encodable") from error


def _decode_canonical_json(payload: bytes, label: str) -> dict[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            _require(key not in result, f"{label} contains duplicate key {key!r}")
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise PhysicalApplicabilityError(f"{label} contains forbidden constant {value}")

    try:
        decoded = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PhysicalApplicabilityError(f"{label} is not canonical JSON") from error
    _require(type(decoded) is dict, f"{label} must decode to an object")
    _require(_canonical_json_bytes(decoded) == payload, f"{label} is not canonical JSON")
    return decoded


def _restart_problem_parameters(payload: bytes, label: str) -> dict[str, str]:
    marker = b"<par_end>"
    _require(marker in payload, f"{label} is missing <par_end>")
    try:
        header = payload.split(marker, 1)[0].decode("latin1")
    except UnicodeDecodeError as error:
        raise PhysicalApplicabilityError(f"{label} header cannot be decoded") from error
    active_block: str | None = None
    parameters: dict[str, str] = {}
    for raw_line in header.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("<") and line.endswith(">"):
            active_block = line[1:-1]
            continue
        if active_block != "problem" or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        _require(key not in parameters, f"{label} contains duplicate problem/{key}")
        parameters[key] = value.split("#", 1)[0].strip()
    return parameters


def _restart_bool(value: str, label: str) -> bool:
    normalized = value.lower()
    _require(normalized in {"true", "false", "1", "0"}, f"{label} is not boolean")
    return normalized in {"true", "1"}


def _sha256_text(value: object, label: str) -> str:
    _require(
        type(value) is str and _SHA256.fullmatch(value) is not None,
        f"{label} must be a lowercase SHA-256 digest",
    )
    return value


def _source_commit(value: object, label: str) -> str:
    _require(
        type(value) is str and _SOURCE_COMMIT.fullmatch(value) is not None,
        f"{label} must be a lowercase full source commit",
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


def _read_bound_artifact(
    evidence_root: object, binding_value: object, *, expected_role: str
) -> tuple[dict[str, object], bytes]:
    _require(isinstance(evidence_root, Path), "evidence root must be a pathlib.Path")
    root = evidence_root.resolve(strict=True)
    _require(root.is_dir(), "evidence root must be a directory")
    binding = _exact_keys(
        binding_value, {"role", "path", "sha256", "byte_count"}, f"{expected_role} binding"
    )
    _require(
        type(binding["role"]) is str and binding["role"] == expected_role,
        f"{expected_role} binding role drifted",
    )
    relative = _safe_relative_path(binding["path"], f"{expected_role} binding path")
    expected_sha = _sha256_text(binding["sha256"], f"{expected_role} binding sha256")
    expected_size = _nonnegative_int(binding["byte_count"], f"{expected_role} byte count")
    path = root / relative
    before = path.lstat()
    _require(stat.S_ISREG(before.st_mode), f"{expected_role} artifact must be regular")
    _require(not path.is_symlink(), f"{expected_role} artifact must not be a symlink")
    with path.open("rb") as stream:
        descriptor_before = os.fstat(stream.fileno())
        payload = stream.read()
        descriptor_after = os.fstat(stream.fileno())
    after = path.lstat()
    identity = lambda status: (
        status.st_dev,
        status.st_ino,
        status.st_mode,
        status.st_size,
        status.st_mtime_ns,
        status.st_ctime_ns,
    )
    _require(
        identity(before)
        == identity(descriptor_before)
        == identity(descriptor_after)
        == identity(after),
        f"{expected_role} artifact changed while read",
    )
    _require(len(payload) == expected_size, f"{expected_role} artifact byte count drifted")
    _require(
        hashlib.sha256(payload).hexdigest() == expected_sha,
        f"{expected_role} artifact digest drifted",
    )
    return {
        "role": expected_role,
        "path": relative,
        "sha256": expected_sha,
        "byte_count": expected_size,
    }, payload


def _array_sha256(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype="<f8"))
    digest = hashlib.sha256()
    digest.update(_canonical_json_bytes({"dtype": "<f8", "shape": list(array.shape)}))
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _typed_array_sha256(values: object, label: str) -> str:
    try:
        array = np.ascontiguousarray(np.asarray(values))
    except (TypeError, ValueError, OverflowError) as error:
        raise PhysicalApplicabilityError(f"{label} must be a numeric array") from error
    _require(
        array.dtype.kind in {"b", "i", "u", "f"} and array.size > 0,
        f"{label} must be a nonempty numeric array",
    )
    if array.dtype.kind == "f":
        _require(np.all(np.isfinite(array)), f"{label} must be finite")
    digest = hashlib.sha256()
    digest.update(
        _canonical_json_bytes({"dtype": array.dtype.str, "shape": list(array.shape)})
    )
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _dataset_sha256(value: object, label: str) -> str:
    required_attributes = (
        "source",
        "time",
        "cycle",
        "location_size",
        "variable_size",
        "variable_names",
        "input_parameters",
        "root_grid_shape",
        "meshblock_shape",
        "nghost",
        "domain_bounds",
        "blocks",
    )
    _require(
        all(hasattr(value, attribute) for attribute in required_attributes),
        f"{label} decoded dataset contract drifted",
    )
    blocks = []
    for index, block in enumerate(value.blocks):
        _require(hasattr(block, "fields"), f"{label} block {index} fields are missing")
        blocks.append(
            {
                "index_bounds": list(block.index_bounds),
                "logical_location": list(block.logical_location),
                "level": block.level,
                "geometry": list(block.geometry),
                "fields": {
                    field: _typed_array_sha256(
                        block.fields[field], f"{label} block {index} field {field}"
                    )
                    for field in sorted(block.fields)
                },
            }
        )
    payload = {
        "source": value.source,
        "time": value.time,
        "cycle": value.cycle,
        "location_size": value.location_size,
        "variable_size": value.variable_size,
        "variable_names": list(value.variable_names),
        "input_parameters": dict(value.input_parameters),
        "root_grid_shape": list(value.root_grid_shape),
        "meshblock_shape": list(value.meshblock_shape),
        "nghost": value.nghost,
        "domain_bounds": list(value.domain_bounds),
        "blocks": blocks,
    }
    return hashlib.sha256(_canonical_json_bytes(payload)).hexdigest()


def _particle_payload_sha256(
    *,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
) -> str:
    payload = {
        "points": _typed_array_sha256(points, "particle points"),
        "cr_source": _typed_array_sha256(cr_source, "particle cr_source"),
        "birth_time": _typed_array_sha256(birth_time, "particle birth time"),
        "velocity": _typed_array_sha256(velocity, "particle velocity"),
        "macro_weight": _typed_array_sha256(macro_weight, "particle macro weight"),
    }
    return hashlib.sha256(_canonical_json_bytes(payload)).hexdigest()


def _decode_particle_vtk_bytes(payload: bytes, label: str) -> dict[str, Any]:
    match = _PVTK_EXECUTION_PATTERN.match(payload[:4096])
    _require(match is not None, f"{label} lacks canonical prtcl_all execution metadata")
    try:
        header_time = float(match.group(1))
        header_nranks = int(match.group(2))
        header_cycle = int(match.group(3))
        header_variables = match.group(4).decode("ascii")
    except (UnicodeDecodeError, ValueError) as error:
        raise PhysicalApplicabilityError(
            f"{label} contains invalid prtcl_all execution metadata"
        ) from error
    _require(
        math.isfinite(header_time)
        and header_nranks > 0
        and header_cycle >= 0
        and header_variables == "prtcl_all",
        f"{label} prtcl_all execution metadata drifted",
    )
    try:
        descriptor = os.memfd_create(
            "q011-bound-particle-payload",
            flags=getattr(os, "MFD_CLOEXEC", 0),
        )
        try:
            remaining = memoryview(payload)
            while remaining:
                written = os.write(descriptor, remaining)
                _require(written > 0, f"{label} could not stage bound particle bytes")
                remaining = remaining[written:]
            particles = pvtk_particles.read_particle_vtk(Path(f"/proc/self/fd/{descriptor}"))
        finally:
            os.close(descriptor)
    except (OSError, ValueError) as error:
        raise PhysicalApplicabilityError(
            f"{label} bound prtcl_all payload cannot be decoded"
        ) from error
    _require(
        set(particles.scalars) == _PVTK_SCALARS and set(particles.vectors) == {"vel"},
        f"{label} prtcl_all field inventory drifted",
    )
    particle_count = particles.points.shape[0]
    _require(
        particle_count > 0
        and particles.points.shape == (particle_count, 3)
        and particles.vectors["vel"].shape == (particle_count, 3),
        f"{label} prtcl_all shape or population is invalid",
    )
    for name in _PVTK_SCALARS:
        _require(
            particles.scalars[name].shape == (particle_count,),
            f"{label} prtcl_all scalar shape drifted",
        )
    _require(
        np.all(np.isfinite(particles.points))
        and np.all(np.isfinite(particles.vectors["vel"]))
        and all(
            np.all(np.isfinite(particles.scalars[name]))
            for name in ("macro_weight", "birth_time", "deltaf_f0", "deltaf_weight")
        )
        and np.all(particles.scalars["macro_weight"] >= 0.0)
        and np.all(particles.scalars["gid"] >= 0)
        and np.all(particles.scalars["ptag"] >= 0)
        and np.all(particles.scalars["species"] >= 0)
        and np.all(np.isin(particles.scalars["cr_source"], (0, 1)))
        and np.unique(particles.scalars["ptag"]).size == particle_count,
        f"{label} prtcl_all values or provenance are invalid",
    )
    return {
        "execution_header": {
            "observed_committed_time": header_time,
            "nranks": header_nranks,
            "cycle": header_cycle,
            "variables": header_variables,
        },
        "points": particles.points,
        "cr_source": particles.scalars["cr_source"],
        "birth_time": particles.scalars["birth_time"],
        "velocity": particles.vectors["vel"],
        "macro_weight": particles.scalars["macro_weight"],
    }


def _validate_bound_normalization_evidence(
    evidence_root: Path,
    value: object,
    *,
    runtime_input_parameters: Mapping[str, Any],
) -> dict[str, Any]:
    evidence = _exact_keys(
        value,
        {
            "runtime_normalization_record",
            "deck",
            "source_manifest",
            "source_archive",
            "executable",
        },
        "bound normalization evidence",
    )
    bindings: dict[str, dict[str, object]] = {}
    payloads: dict[str, bytes] = {}
    for role in (
        "runtime_normalization_record",
        "deck",
        "source_manifest",
        "source_archive",
        "executable",
    ):
        bindings[role], payloads[role] = _read_bound_artifact(
            evidence_root, evidence[role], expected_role=role
        )
        _require(len(payloads[role]) > 0, f"{role} artifact must not be empty")
    runtime = _decode_canonical_json(
        payloads["runtime_normalization_record"], "runtime normalization record"
    )
    runtime = _exact_keys(
        runtime,
        {
            "schema_version",
            "record_type",
            "source_commit",
            "normalization",
            "runtime_input_parameters_sha256",
            "deck_sha256",
            "source_manifest_sha256",
            "source_archive_sha256",
            "executable_sha256",
        },
        "runtime normalization record",
    )
    _require(
        type(runtime["schema_version"]) is int and runtime["schema_version"] == SCHEMA_VERSION,
        "runtime normalization schema version drifted",
    )
    _require(
        runtime["record_type"] == NORMALIZATION_RECORD_TYPE,
        "runtime normalization record type drifted",
    )
    source_commit = _source_commit(runtime["source_commit"], "runtime normalization source commit")
    _require(
        source_commit == TRUSTED_Q011_RUNTIME_SOURCE_COMMIT,
        "runtime normalization is not bound to the exact trusted Q011 source identity",
    )
    normalization = _validate_exact_normalization(runtime["normalization"])
    source_manifest = _decode_canonical_json(
        payloads["source_manifest"], "bound source manifest"
    )
    source_manifest = _exact_keys(
        source_manifest,
        {
            "schema_version",
            "record_type",
            "source_commit",
            "deck_sha256",
            "source_archive_sha256",
            "executable_sha256",
        },
        "bound source manifest",
    )
    _require(
        type(source_manifest["schema_version"]) is int
        and source_manifest["schema_version"] == SCHEMA_VERSION
        and source_manifest["record_type"] == SOURCE_MANIFEST_RECORD_TYPE,
        "bound source manifest identity drifted",
    )
    _require(
        _source_commit(source_manifest["source_commit"], "bound source manifest commit")
        == source_commit,
        "bound source manifest commit disagrees with runtime normalization",
    )
    expected_input_sha = hashlib.sha256(
        _canonical_json_bytes(dict(runtime_input_parameters))
    ).hexdigest()
    _require(
        _sha256_text(
            runtime["runtime_input_parameters_sha256"],
            "runtime input-parameter digest",
        )
        == expected_input_sha,
        "runtime normalization is not bound to decoded runtime input parameters",
    )
    for role in ("deck", "source_manifest", "source_archive", "executable"):
        _require(
            _sha256_text(runtime[f"{role}_sha256"], f"runtime {role} digest")
            == bindings[role]["sha256"],
            f"runtime normalization {role} digest disagrees with bound artifact",
        )
    for role in ("deck", "source_archive", "executable"):
        _require(
            _sha256_text(
                source_manifest[f"{role}_sha256"], f"source manifest {role} digest"
            )
            == bindings[role]["sha256"],
            f"source manifest {role} digest disagrees with bound artifact",
        )
    return {
        "source_commit": source_commit,
        "normalization": normalization,
        "runtime_input_parameters_sha256": expected_input_sha,
        "source_manifest": dict(source_manifest),
        "bindings": bindings,
    }


def _validate_snapshot_provenance(
    evidence_root: Path,
    value: object,
    *,
    normalization_evidence: Mapping[str, Any],
    mhd_dataset: object,
    current_datasets: Mapping[str, object],
    particle_source: object,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
    nominal_slot_time: float,
    observed_committed_time: float,
) -> dict[str, Any]:
    envelope = _exact_keys(value, {"snapshot_manifest"}, "snapshot provenance evidence")
    manifest_binding, manifest_payload = _read_bound_artifact(
        evidence_root, envelope["snapshot_manifest"], expected_role="snapshot_manifest"
    )
    manifest = _decode_canonical_json(manifest_payload, "snapshot provenance manifest")
    manifest = _exact_keys(
        manifest,
        {
            "schema_version",
            "record_type",
            "attempt_id",
            "source_commit",
            "executable_sha256",
            "runtime_normalization_sha256",
            "nominal_slot_time",
            "observed_committed_time",
            "cycle",
            "raw_products",
            "decoded_product_sha256",
        },
        "snapshot provenance manifest",
    )
    _require(
        type(manifest["schema_version"]) is int and manifest["schema_version"] == SCHEMA_VERSION,
        "snapshot provenance schema version drifted",
    )
    _require(
        manifest["record_type"] == SNAPSHOT_PROVENANCE_RECORD_TYPE,
        "snapshot provenance record type drifted",
    )
    _require(
        type(manifest["attempt_id"]) is str
        and _ATTEMPT_ID.fullmatch(manifest["attempt_id"]) is not None,
        "snapshot provenance attempt id drifted",
    )
    _require(
        _source_commit(manifest["source_commit"], "snapshot source commit")
        == normalization_evidence["source_commit"],
        "snapshot source commit disagrees with normalization evidence",
    )
    _require(
        _sha256_text(manifest["executable_sha256"], "snapshot executable digest")
        == normalization_evidence["bindings"]["executable"]["sha256"],
        "snapshot executable disagrees with normalization evidence",
    )
    _require(
        _sha256_text(
            manifest["runtime_normalization_sha256"],
            "snapshot runtime normalization digest",
        )
        == normalization_evidence["bindings"]["runtime_normalization_record"]["sha256"],
        "snapshot runtime normalization digest drifted",
    )
    _require(
        _finite_scalar(manifest["nominal_slot_time"], "snapshot provenance nominal time")
        == nominal_slot_time
        and _finite_scalar(
            manifest["observed_committed_time"], "snapshot provenance observed time"
        )
        == observed_committed_time,
        "snapshot provenance times drifted",
    )
    _require(
        type(manifest["cycle"]) is int
        and manifest["cycle"] == mhd_dataset.cycle,
        "snapshot provenance cycle drifted",
    )
    raw_products = _exact_keys(
        manifest["raw_products"], set(REQUIRED_RAW_PRODUCTS), "snapshot raw products"
    )
    supplied_datasets = {
        "mhd_w_bcc": mhd_dataset,
        **{product: current_datasets[product] for product in science.CURRENT_PRODUCT_FIELDS},
    }
    validated_products: dict[str, dict[str, object]] = {}
    raw_payloads: dict[str, bytes] = {}
    for product in REQUIRED_RAW_PRODUCTS:
        binding, raw_payloads[product] = _read_bound_artifact(
            evidence_root, raw_products[product], expected_role=product
        )
        validated_products[product] = binding
    decoded_expected: dict[str, str] = {}
    for product, supplied in supplied_datasets.items():
        _require(
            supplied.source == validated_products[product]["path"],
            f"{product} decoded source is not bound to raw artifact",
        )
        try:
            parsed = output_primitives.parse_athenak_binary_bytes(
                raw_payloads[product], source=str(validated_products[product]["path"])
            )
        except output_primitives.AnalysisError as error:
            raise PhysicalApplicabilityError(
                f"{product} bound raw product is not a trusted Athena binary: {error}"
            ) from error
        parsed_sha = _dataset_sha256(parsed, f"trusted decoded {product}")
        _require(
            parsed_sha == _dataset_sha256(supplied, f"supplied decoded {product}"),
            f"{product} supplied decoded dataset disagrees with bound raw bytes",
        )
        decoded_expected[product] = parsed_sha
    _require(
        particle_source == validated_products["prtcl_all"]["path"],
        "prtcl_all decoded source is not bound to raw artifact",
    )
    decoded_particle = _decode_particle_vtk_bytes(
        raw_payloads["prtcl_all"], "snapshot bound prtcl_all"
    )
    _require(
        decoded_particle["execution_header"]["observed_committed_time"]
        == observed_committed_time
        and decoded_particle["execution_header"]["cycle"] == mhd_dataset.cycle,
        "snapshot bound prtcl_all cycle/time binding drifted",
    )
    decoded_expected["prtcl_all"] = _particle_payload_sha256(
        points=decoded_particle["points"],
        cr_source=decoded_particle["cr_source"],
        birth_time=decoded_particle["birth_time"],
        velocity=decoded_particle["velocity"],
        macro_weight=decoded_particle["macro_weight"],
    )
    _require(
        decoded_expected["prtcl_all"]
        == _particle_payload_sha256(
            points=points,
            cr_source=cr_source,
            birth_time=birth_time,
            velocity=velocity,
            macro_weight=macro_weight,
        ),
        "prtcl_all supplied decoded payload disagrees with bound raw bytes",
    )
    decoded = _exact_keys(
        manifest["decoded_product_sha256"],
        set(REQUIRED_RAW_PRODUCTS),
        "snapshot decoded-product digests",
    )
    for product, expected_sha in decoded_expected.items():
        _require(
            _sha256_text(decoded[product], f"snapshot decoded {product} digest")
            == expected_sha,
            f"snapshot decoded {product} payload disagrees with trusted raw-byte reduction",
        )
    return {
        "manifest_binding": manifest_binding,
        "attempt_id": manifest["attempt_id"],
        "source_commit": manifest["source_commit"],
        "executable_sha256": manifest["executable_sha256"],
        "runtime_normalization_sha256": manifest["runtime_normalization_sha256"],
        "nominal_slot_time": nominal_slot_time,
        "observed_committed_time": observed_committed_time,
        "cycle": manifest["cycle"],
        "raw_products": validated_products,
        "decoded_product_sha256": decoded_expected,
    }


def _finite_array(values: object, label: str, *, ndim: int | None = None) -> np.ndarray:
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise PhysicalApplicabilityError(f"{label} must be numeric") from error
    if ndim is not None:
        _require(array.ndim == ndim, f"{label} must be {ndim}-dimensional")
    _require(array.size > 0, f"{label} must not be empty")
    _require(np.all(np.isfinite(array)), f"{label} must be finite")
    return array


def _readonly(values: np.ndarray) -> np.ndarray:
    result = np.array(values, copy=True)
    result.setflags(write=False)
    return result


def _weighted_quantile(
    values: np.ndarray, weights: np.ndarray, quantile: float, *, label: str
) -> float:
    samples = _finite_array(values, f"{label} values").reshape(-1)
    sample_weights = _finite_array(weights, f"{label} weights").reshape(-1)
    _require(samples.size == sample_weights.size, f"{label} sizes disagree")
    _require(np.all(sample_weights >= 0.0), f"{label} weights must be non-negative")
    positive = sample_weights > 0.0
    _require(np.any(positive), f"{label} requires positive total weight")
    order = np.argsort(samples[positive], kind="stable")
    ordered = samples[positive][order]
    ordered_weights = sample_weights[positive][order]
    cumulative = np.cumsum(ordered_weights)
    threshold = quantile * float(cumulative[-1])
    index = int(np.searchsorted(cumulative, threshold, side="left"))
    return float(ordered[min(index, ordered.size - 1)])


def _weighted_statistics(
    values: np.ndarray,
    area_weights: np.ndarray,
    current_weights: np.ndarray,
    *,
    label: str,
) -> dict[str, Any]:
    samples = _finite_array(values, f"{label} values").reshape(-1)
    area = _finite_array(area_weights, f"{label} area weights").reshape(-1)
    current = _finite_array(current_weights, f"{label} current weights").reshape(-1)
    _require(samples.size == area.size == current.size, f"{label} sizes disagree")
    _require(np.all(area > 0.0), f"{label} area weights must be positive")
    _require(np.all(current >= 0.0), f"{label} current weights must be non-negative")

    def summarize(weights: np.ndarray, weight_label: str) -> dict[str, Any]:
        total = float(np.sum(weights))
        if total == 0.0:
            return {
                "available": False,
                "reason": f"zero_total_{weight_label}_weight",
                "total_weight": 0.0,
                "weighted_mean": None,
                "weighted_quantiles": None,
            }
        return {
            "available": True,
            "reason": None,
            "total_weight": total,
            "weighted_mean": float(np.sum(samples * weights) / total),
            "weighted_quantiles": {
                "q500": _weighted_quantile(samples, weights, 0.5, label=label),
                "q900": _weighted_quantile(samples, weights, 0.9, label=label),
                "q990": _weighted_quantile(samples, weights, 0.99, label=label),
                "q999": _weighted_quantile(samples, weights, 0.999, label=label),
            },
        }

    return {
        "cell_count": int(samples.size),
        "local_minimum": float(np.min(samples)),
        "local_maximum": float(np.max(samples)),
        "area_weighted": summarize(area, "area"),
        "gas_frame_current_weighted": summarize(current, "gas_frame_current"),
    }


def _validate_exact_normalization(value: object) -> dict[str, object]:
    normalization = _exact_keys(
        value, set(EXACT_NORMALIZATION), "Q011 exact normalization"
    )
    for key, expected in EXACT_NORMALIZATION.items():
        _require(
            type(normalization[key]) is type(expected) and normalization[key] == expected,
            f"Q011 exact normalization field {key!r} drifted",
        )
    return dict(EXACT_NORMALIZATION)


def _uniform_spacing(faces: np.ndarray, label: str) -> float:
    widths = np.diff(_finite_array(faces, f"{label} faces", ndim=1))
    _require(np.all(widths > 0.0), f"{label} faces must increase")
    reference = float(widths[0])
    _require(
        np.allclose(widths, reference, rtol=0.0, atol=1.0e-12 * max(1.0, reference)),
        f"{label} composite spacing must be uniform",
    )
    return reference


def _cell_areas(state: science.ComposedMHDState) -> np.ndarray:
    return np.diff(state.x2_faces)[:, None] * np.diff(state.x1_faces)[None, :]


def _x1_centers(state: science.ComposedMHDState) -> np.ndarray:
    return 0.5 * (state.x1_faces[:-1] + state.x1_faces[1:])


def _detected_front(
    mhd_dataset: object,
    *,
    nominal_slot_time: object,
    observed_committed_time: object,
    target_level: int | None,
) -> tuple[science.ComposedMHDState, float]:
    state = science.compose_full_mhd_state(
        mhd_dataset,
        nominal_slot_time=nominal_slot_time,
        observed_committed_time=observed_committed_time,
        target_level=target_level,
    )
    reduced = science.reduce_mhd_snapshot(
        mhd_dataset,
        nominal_slot_time=nominal_slot_time,
        observed_committed_time=observed_committed_time,
        target_level=target_level,
    )
    return state, float(reduced["detected_front"]["x_front_c_over_omega_pi"])


def _region_masks(
    state: science.ComposedMHDState, front_x: float
) -> dict[str, np.ndarray]:
    x = _x1_centers(state)
    shape = state.fields_y_x["dens"].shape
    masks: dict[str, np.ndarray] = {}
    for name, offsets in DETECTED_FRONT_REGIONS.items():
        if offsets is None:
            masks[name] = np.ones(shape, dtype=bool)
            continue
        lower = front_x + offsets[0]
        upper = front_x + offsets[1]
        _require(
            lower >= state.x1_faces[0] and upper <= state.x1_faces[-1],
            f"{name} escaped the retained x1 domain",
        )
        selected_x = (x > lower) & (x < upper)
        _require(np.count_nonzero(selected_x) >= 4, f"{name} requires at least four x1 cells")
        masks[name] = np.broadcast_to(selected_x[None, :], shape)
    return masks


def _actual_leaf_spacing_maps(
    state: science.ComposedMHDState,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dx1_target = _uniform_spacing(state.x1_faces, "x1")
    dx2_target = _uniform_spacing(state.x2_faces, "x2")
    level_factor = np.power(
        2.0, state.target_level - state.source_levels_y_x.astype(np.float64)
    )
    _require(np.all(level_factor >= 1.0), "source level exceeds target composite level")
    dx1_leaf = dx1_target * level_factor
    dx2_leaf = dx2_target * level_factor
    delta_leaf = np.minimum(dx1_leaf, dx2_leaf)
    return dx1_leaf, dx2_leaf, delta_leaf


def _leaf_aware_magnetic_spectrum(
    state: science.ComposedMHDState,
    mask: np.ndarray,
    *,
    di_max: float,
) -> dict[str, Any]:
    selected_x = np.any(mask, axis=0)
    selected_y = np.any(mask, axis=1)
    _require(
        np.all(mask == np.outer(selected_y, selected_x)),
        "magnetic spectrum region must be rectangular",
    )
    target_dx1 = _uniform_spacing(state.x1_faces, "x1")
    target_dx2 = _uniform_spacing(state.x2_faces, "x2")
    source_levels = state.source_levels_y_x[np.ix_(selected_y, selected_x)]
    levels_present = sorted(int(level) for level in np.unique(source_levels))
    analysis_level = min(levels_present)
    factor = 2 ** (state.target_level - analysis_level)
    _require(factor >= 1, "leaf-aware spectrum restriction factor is invalid")

    def complete_groups(selected: np.ndarray, count: int, label: str) -> np.ndarray:
        groups = []
        for group in np.unique(np.flatnonzero(selected) // factor):
            indices = np.arange(group * factor, (group + 1) * factor)
            if indices[-1] < count and np.all(selected[indices]):
                groups.append(int(group))
        result = np.asarray(groups, dtype=np.int64)
        _require(
            result.size >= 4 and np.all(np.diff(result) == 1),
            f"leaf-aware spectrum requires at least four contiguous complete {label} groups",
        )
        return result

    x_groups = complete_groups(selected_x, selected_x.size, "x1")
    y_groups = complete_groups(selected_y, selected_y.size, "x2")
    x_indices = np.arange(x_groups[0] * factor, (x_groups[-1] + 1) * factor)
    y_indices = np.arange(y_groups[0] * factor, (y_groups[-1] + 1) * factor)
    magnetic_target = np.stack(
        [
            state.fields_y_x[field][np.ix_(y_indices, x_indices)]
            for field in ("bcc1", "bcc2", "bcc3")
        ],
        axis=0,
    )
    ny = y_groups.size
    nx = x_groups.size
    magnetic_restricted = magnetic_target.reshape(3, ny, factor, nx, factor).mean(
        axis=(2, 4)
    )
    mean_b = np.mean(magnetic_target, axis=(1, 2), keepdims=True)
    target_fluctuation = magnetic_target - mean_b
    restricted_fluctuation = magnetic_restricted - mean_b
    delta_b_rms = math.sqrt(
        float(np.mean(np.sum(target_fluctuation * target_fluctuation, axis=0)))
    )
    prolonged = np.repeat(
        np.repeat(magnetic_restricted, factor, axis=1), factor, axis=2
    )
    subgrid_residual = magnetic_target - prolonged
    window = np.outer(np.hanning(ny), np.hanning(nx))
    transformed = np.fft.fftn(
        restricted_fluctuation * window[None, :, :], axes=(-2, -1)
    )
    power = np.sum(np.abs(transformed) ** 2, axis=0) / float(nx * ny)
    power[0, 0] = 0.0
    resolved_power_scale = float(factor * factor)
    resolved_power = float(np.sum(power) * resolved_power_scale)
    fine_window = np.repeat(np.repeat(window, factor, axis=0), factor, axis=1)
    subgrid_power = float(
        np.sum(
            subgrid_residual
            * subgrid_residual
            * fine_window[None, :, :]
            * fine_window[None, :, :]
        )
    )
    total_power = resolved_power + subgrid_power
    _require(
        math.isfinite(total_power) and total_power > np.finfo(np.float64).tiny,
        "precursor leaf-aware magnetic estimator has no resolvable fluctuation power",
    )
    analysis_dx1 = target_dx1 * factor
    analysis_dx2 = target_dx2 * factor
    k1 = 2.0 * math.pi * np.fft.fftfreq(nx, d=analysis_dx1)
    k2 = 2.0 * math.pi * np.fft.fftfreq(ny, d=analysis_dx2)
    k1_grid, k2_grid = np.meshgrid(k1, k2, indexing="xy")
    kmag = np.sqrt(k1_grid * k1_grid + k2_grid * k2_grid)
    nonzero = kmag > 0.0
    minimum_actual_leaf_spacing = min(
        target_dx1 * 2 ** (state.target_level - max(levels_present)),
        target_dx2 * 2 ** (state.target_level - max(levels_present)),
    )
    subgrid_k_upper = math.pi / minimum_actual_leaf_spacing
    mean_k = float(
        (
            np.sum(kmag[nonzero] * power[nonzero]) * resolved_power_scale
            + subgrid_k_upper * subgrid_power
        )
        / total_power
    )
    _require(mean_k > 0.0 and math.isfinite(mean_k), "characteristic magnetic k is invalid")
    lambda_char = 2.0 * math.pi / mean_k
    sub_10di = kmag > 2.0 * math.pi / (10.0 * di_max)
    sub_fraction_upper = float(
        (
            np.sum(power[sub_10di]) * resolved_power_scale
            + subgrid_power
        )
        / total_power
    )
    return {
        "region": DI_MAGNETIC_SPECTRUM_REGION,
        "method": (
            "actual_leaf_aware_coarsest_present_level_conservative_restriction_"
            "FFT_plus_all_subgrid_residual_as_high_k_upper_bound"
        ),
        "actual_leaf_aware": True,
        "finest_composite_FFT_rejected": True,
        "source_levels_present": levels_present,
        "analysis_source_level": analysis_level,
        "target_composite_level": state.target_level,
        "restriction_factor": factor,
        "analysis_dx1": analysis_dx1,
        "analysis_dx2": analysis_dx2,
        "selected_analysis_nx1": nx,
        "selected_analysis_nx2": ny,
        "delta_B_rms_over_B0": delta_b_rms / float(EXACT_NORMALIZATION["reference_b0"]),
        "total_windowed_delta_B_power": total_power,
        "resolved_restricted_power": resolved_power,
        "subgrid_residual_power_upper_bound": subgrid_power,
        "subgrid_residual_power_fraction_upper_bound": subgrid_power / total_power,
        "local_di_maximum": di_max,
        "lambda_B_characteristic": lambda_char,
        "lambda_B_characteristic_over_local_di_maximum": lambda_char / di_max,
        "sub_10di_magnetic_power_fraction_upper_bound": sub_fraction_upper,
        "shock_transition_excluded": True,
    }


def _tsc_axis_indices_weights(
    coordinates: np.ndarray,
    faces: np.ndarray,
    *,
    periodic: bool,
    label: str,
) -> tuple[np.ndarray, np.ndarray]:
    dx = _uniform_spacing(faces, label)
    count = faces.size - 1
    normalized = (coordinates - faces[0]) / dx - 0.5
    center = np.floor(normalized + 0.5).astype(np.int64)
    indices = center[:, None] + np.asarray([-1, 0, 1], dtype=np.int64)[None, :]
    distance = np.abs(normalized[:, None] - indices.astype(np.float64))
    weights = np.where(
        distance < 0.5,
        0.75 - distance * distance,
        np.where(distance < 1.5, 0.5 * (1.5 - distance) ** 2, 0.0),
    )
    _require(
        np.allclose(np.sum(weights, axis=1), 1.0, rtol=0.0, atol=2.0e-15),
        f"{label} TSC weights failed unity closure",
    )
    if periodic:
        indices %= count
    else:
        _require(
            np.all((indices >= 0) & (indices < count)),
            f"{label} particle TSC stencil escaped the retained nonperiodic domain",
        )
    return indices, weights


def _tsc_sample_magnetic_field(
    state: science.ComposedMHDState, points: np.ndarray
) -> np.ndarray:
    _require(points.ndim == 2 and points.shape[1] == 3, "particle points shape drifted")
    _require(
        np.all(points[:, 0] >= state.x1_faces[0])
        and np.all(points[:, 0] <= state.x1_faces[-1])
        and np.all(points[:, 1] >= state.x2_faces[0])
        and np.all(points[:, 1] <= state.x2_faces[-1]),
        "particle escaped retained x1-x2 domain",
    )
    x_indices, x_weights = _tsc_axis_indices_weights(
        points[:, 0], state.x1_faces, periodic=False, label="x1"
    )
    y_indices, y_weights = _tsc_axis_indices_weights(
        points[:, 1], state.x2_faces, periodic=True, label="x2"
    )
    sampled = np.zeros((points.shape[0], 3), dtype=np.float64)
    for component, field in enumerate(("bcc1", "bcc2", "bcc3")):
        values = state.fields_y_x[field]
        for iy in range(3):
            for ix in range(3):
                sampled[:, component] += (
                    y_weights[:, iy]
                    * x_weights[:, ix]
                    * values[y_indices[:, iy], x_indices[:, ix]]
                )
    _require(np.all(np.isfinite(sampled)), "TSC sampled magnetic field must be finite")
    return sampled


def _tsc_sample_scalar_field(
    state: science.ComposedMHDState,
    values_y_x: np.ndarray,
    points: np.ndarray,
    *,
    label: str,
) -> np.ndarray:
    _require(
        values_y_x.shape == state.fields_y_x["dens"].shape,
        f"{label} scalar field shape drifted",
    )
    x_indices, x_weights = _tsc_axis_indices_weights(
        points[:, 0], state.x1_faces, periodic=False, label="x1"
    )
    y_indices, y_weights = _tsc_axis_indices_weights(
        points[:, 1], state.x2_faces, periodic=True, label="x2"
    )
    sampled = np.zeros(points.shape[0], dtype=np.float64)
    for iy in range(3):
        for ix in range(3):
            sampled += (
                y_weights[:, iy]
                * x_weights[:, ix]
                * values_y_x[y_indices[:, iy], x_indices[:, ix]]
            )
    _require(np.all(np.isfinite(sampled)), f"TSC sampled {label} must be finite")
    return sampled


def _particle_scalar_exposure_statistics(
    values: np.ndarray,
    macro_weights: np.ndarray,
    energy_weights: np.ndarray,
    *,
    threshold: float,
    label: str,
) -> dict[str, Any]:
    samples = _finite_array(values, f"{label} samples", ndim=1)
    macro = _finite_array(macro_weights, f"{label} macro weights", ndim=1)
    energy = _finite_array(energy_weights, f"{label} energy weights", ndim=1)
    _require(samples.size == macro.size == energy.size, f"{label} sizes disagree")
    _require(
        np.all(macro >= 0.0) and np.all(energy >= 0.0),
        f"{label} particle weights must be non-negative",
    )

    def weighted(weights: np.ndarray, weight_label: str) -> dict[str, Any]:
        total = float(np.sum(weights))
        _require(total > 0.0, f"{label} {weight_label} total weight must be positive")
        return {
            "weighted_mean": float(np.sum(samples * weights) / total),
            "weighted_quantiles": {
                "q500": _weighted_quantile(samples, weights, 0.5, label=label),
                "q900": _weighted_quantile(samples, weights, 0.9, label=label),
                "q990": _weighted_quantile(samples, weights, 0.99, label=label),
                "q999": _weighted_quantile(samples, weights, 0.999, label=label),
            },
            "exceedance_fraction": float(np.sum(weights[samples > threshold]) / total),
        }

    return {
        "particle_count": int(samples.size),
        "local_minimum": float(np.min(samples)),
        "local_maximum": float(np.max(samples)),
        "threshold": threshold,
        "macro_weighted": weighted(macro, "macro"),
        "CR_kinetic_energy_weighted": weighted(energy, "CR kinetic energy"),
    }


def _particle_R_Lambda_exposure(
    state: science.ComposedMHDState,
    front_x: float,
    r_map: np.ndarray,
    lambda_map: np.ndarray,
    *,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
) -> dict[str, Any]:
    arrays = science._decoded_particle_arrays(
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )
    energetic = arrays["energetic"]
    _require(np.any(energetic), "particle exposure requires positive-weight active CRs")
    selected_points = arrays["points"][energetic]
    macro = arrays["weights"][energetic]
    energy = macro * arrays["specific_kinetic_energy"][energetic]
    _require(float(np.sum(energy)) > 0.0, "particle exposure energy weight must be positive")
    sampled_r = _tsc_sample_scalar_field(
        state, r_map, selected_points, label="particle R exposure"
    )
    sampled_lambda = _tsc_sample_scalar_field(
        state, lambda_map, selected_points, label="particle Lambda exposure"
    )
    high_energy_threshold = _weighted_quantile(
        arrays["specific_kinetic_energy"][energetic],
        macro,
        0.99,
        label="particle high-energy-tail threshold",
    )
    populations = {
        "all_active": np.ones(sampled_r.size, dtype=bool),
        "detected_front_upstream": (
            selected_points[:, 0] > front_x + SHOCK_TRANSITION_HALF_WIDTH
        ),
        "high_energy_tail": (
            arrays["specific_kinetic_energy"][energetic] >= high_energy_threshold
        ),
    }
    result: dict[str, Any] = {}
    for name, selected in populations.items():
        if not np.any(selected):
            result[name] = {
                "available": False,
                "reason": "population_contains_no_positive_weight_active_particles",
                "particle_count": 0,
                "R": None,
                "Lambda": None,
            }
            continue
        result[name] = {
            "available": True,
            "reason": None,
            "particle_count": int(np.count_nonzero(selected)),
            "R": _particle_scalar_exposure_statistics(
                sampled_r[selected],
                macro[selected],
                energy[selected],
                threshold=R_MAXIMUM,
                label=f"{name} particle R exposure",
            ),
            "Lambda": _particle_scalar_exposure_statistics(
                sampled_lambda[selected],
                macro[selected],
                energy[selected],
                threshold=LAMBDA_MAXIMUM,
                label=f"{name} particle Lambda exposure",
            ),
        }
    return {
        "sampling": (
            "TSC_on_matched_finest_composite_periodic_x2_nonperiodic_x1_"
            "stencil_must_remain_retained"
        ),
        "population_definitions": {
            "all_active": "cr_source_eq_1_birth_time_ge_45_positive_macro_weight",
            "detected_front_upstream": "all_active_and_x1_gt_detected_front_plus_120",
            "high_energy_tail": (
                "all_active_with_specific_kinetic_energy_ge_macro_weighted_q990"
            ),
        },
        "high_energy_tail_specific_kinetic_energy_threshold": high_energy_threshold,
        "populations": result,
    }


def _particle_gyroradius(
    state: science.ComposedMHDState,
    *,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
) -> dict[str, Any]:
    arrays = science._decoded_particle_arrays(
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )
    selected = arrays["energetic"]
    count = int(np.count_nonzero(selected))
    _require(
        count >= PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES,
        "particle gyroradius q999 requires at least 1000 positive-weight active CRs",
    )
    selected_points = arrays["points"][selected]
    sampled_b = _tsc_sample_magnetic_field(state, selected_points)
    bmag = np.linalg.norm(sampled_b, axis=1)
    _require(np.all(bmag > 0.0), "particle local magnetic field must be nonzero")
    speed_squared = np.sum(arrays["velocity"][selected] ** 2, axis=1)
    gamma = 1.0 / np.sqrt(1.0 - speed_squared / science.PARTICLE_LIGHT_SPEED**2)
    momentum_magnitude = gamma * np.sqrt(speed_squared)
    rg = momentum_magnitude / bmag
    _require(np.all(np.isfinite(rg)) and np.all(rg >= 0.0), "particle gyroradius is invalid")
    macro = arrays["weights"][selected]
    specific_energy = arrays["specific_kinetic_energy"][selected]
    energy = macro * specific_energy
    _require(float(np.sum(energy)) > 0.0, "particle energy weights must be positive")
    ly = float(state.x2_faces[-1] - state.x2_faces[0])
    macro_q999 = _weighted_quantile(rg, macro, 0.999, label="particle gyroradius macro")
    energy_q999 = _weighted_quantile(rg, energy, 0.999, label="particle gyroradius energy")
    high_energy_threshold = _weighted_quantile(
        specific_energy,
        macro,
        0.99,
        label="particle gyroradius high-energy-tail threshold",
    )
    high_energy_tail = specific_energy >= high_energy_threshold
    _require(np.any(high_energy_tail), "particle gyroradius high-energy tail is empty")
    maximum = float(np.max(rg))
    high_energy_maximum = float(np.max(rg[high_energy_tail]))
    energy_fraction_above_ly_over_4 = float(
        np.sum(energy[rg > ly / 4.0]) / np.sum(energy)
    )
    return {
        "population": "cr_source_eq_1_birth_time_ge_45_positive_macro_weight",
        "local_B_sampling": (
            "TSC_on_matched_finest_composite_periodic_x2_nonperiodic_x1_"
            "stencil_must_remain_retained"
        ),
        "gyroradius_formula": "abs_gamma_v_over_abs_q_over_m_times_abs_B_local",
        "selected_particle_count": count,
        "transverse_domain_size_Ly": ly,
        "macro_weighted": {
            "q500": _weighted_quantile(rg, macro, 0.5, label="particle gyroradius macro"),
            "q900": _weighted_quantile(rg, macro, 0.9, label="particle gyroradius macro"),
            "q990": _weighted_quantile(rg, macro, 0.99, label="particle gyroradius macro"),
            "q999": macro_q999,
        },
        "energy_weighted": {
            "q500": _weighted_quantile(rg, energy, 0.5, label="particle gyroradius energy"),
            "q900": _weighted_quantile(rg, energy, 0.9, label="particle gyroradius energy"),
            "q990": _weighted_quantile(rg, energy, 0.99, label="particle gyroradius energy"),
            "q999": energy_q999,
        },
        "maximum": maximum,
        "maximum_specific_kinetic_energy": float(np.max(specific_energy)),
        "high_energy_tail_definition": (
            "specific_kinetic_energy_ge_macro_weighted_q990"
        ),
        "high_energy_tail_specific_kinetic_energy_threshold": high_energy_threshold,
        "high_energy_tail_particle_count": int(np.count_nonzero(high_energy_tail)),
        "high_energy_tail_maximum": high_energy_maximum,
        "macro_q999_over_Ly": macro_q999 / ly,
        "energy_q999_over_Ly": energy_q999 / ly,
        "maximum_over_Ly": maximum / ly,
        "high_energy_tail_maximum_over_Ly": high_energy_maximum / ly,
        "energy_fraction_with_rg_above_Ly_over_4": energy_fraction_above_ly_over_4,
    }


def _claim_rejections(
    gates: Mapping[str, Mapping[str, Any]],
    maximum_lambda: float,
    claim_specific_escape_applicability: Mapping[str, Mapping[str, Any]] | None = None,
) -> list[str]:
    rejected = set(PERMANENT_CLAIM_EXCLUSIONS)
    for gate, claims in CLAIM_REJECTION_RULES.items():
        if not bool(gates[gate]["pass"]):
            rejected.update(claims)
    if maximum_lambda >= 1.0:
        rejected.add("no_Hall_run_approximates_target_plasma_claim")
    if claim_specific_escape_applicability is not None:
        for claim, applicability in claim_specific_escape_applicability.items():
            if not bool(applicability["pass"]):
                rejected.add(claim)
    return sorted(rejected)


@dataclass(frozen=True)
class ApplicabilitySnapshot:
    """One JSON-ready summary plus immutable per-cell applicability maps."""

    record: Mapping[str, Any]
    cell_maps: Mapping[str, np.ndarray]
    evidence_root: Path


@_public_contract("Q011 physical-applicability snapshot reduction")
def reduce_physical_applicability_snapshot(
    mhd_dataset: object,
    current_datasets: Mapping[str, object],
    *,
    evidence_root: Path,
    normalization_evidence: object,
    snapshot_provenance: object,
    particle_source: object,
    nominal_slot_time: object,
    observed_committed_time: object,
    points: object,
    cr_source: object,
    birth_time: object,
    velocity: object,
    macro_weight: object,
    target_level: int | None = None,
) -> ApplicabilitySnapshot:
    """Reduce one matched snapshot into fail-closed physical-applicability gates."""
    nominal = _finite_scalar(nominal_slot_time, "nominal slot time", minimum=0.0)
    observed = _finite_scalar(
        observed_committed_time, "observed committed time", minimum=0.0
    )
    bound_normalization = _validate_bound_normalization_evidence(
        evidence_root,
        normalization_evidence,
        runtime_input_parameters=mhd_dataset.input_parameters,
    )
    normalized = bound_normalization["normalization"]
    provenance = _validate_snapshot_provenance(
        evidence_root,
        snapshot_provenance,
        normalization_evidence=bound_normalization,
        mhd_dataset=mhd_dataset,
        current_datasets=current_datasets,
        particle_source=particle_source,
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
        nominal_slot_time=nominal,
        observed_committed_time=observed,
    )
    state, front_x = _detected_front(
        mhd_dataset,
        nominal_slot_time=nominal,
        observed_committed_time=observed,
        target_level=target_level,
    )
    currents = science._compose_matched_current_fields(mhd_dataset, state, current_datasets)
    rho_g = state.fields_y_x["dens"]
    rho_q = currents["prtcl_rho"]
    gas_velocity = np.stack(
        [state.fields_y_x[f"vel{component}"] for component in ("x", "y", "z")],
        axis=0,
    )
    j_lab = np.stack(
        [currents[f"prtcl_j{component}"] for component in ("x", "y", "z")],
        axis=0,
    )
    j_gas = j_lab - rho_q[None, :, :] * gas_velocity
    j_gas_magnitude = np.linalg.norm(j_gas, axis=0)
    bmag = np.sqrt(
        sum(state.fields_y_x[field] ** 2 for field in ("bcc1", "bcc2", "bcc3"))
    )
    _require(np.all(bmag > 0.0), "cell magnetic-field magnitude must be positive")
    va = bmag / np.sqrt(rho_g)
    denominator = rho_g + rho_q
    _require(np.all(denominator > 0.0), "rho_g plus rho_q must be positive")
    r_map = rho_q / denominator
    lambda_map = j_gas_magnitude / (denominator * va)
    di_map = 1.0 / np.sqrt(rho_g)
    dx1_leaf, dx2_leaf, delta_leaf = _actual_leaf_spacing_maps(state)
    s_delta_map = delta_leaf / di_map
    for name, values in {
        "R": r_map,
        "Lambda": lambda_map,
        "d_i": di_map,
        "S_delta": s_delta_map,
    }.items():
        _require(
            np.all(np.isfinite(values)) and np.all(values >= 0.0),
            f"{name} cell map must be finite and non-negative",
        )

    areas = _cell_areas(state)
    masks = _region_masks(state, front_x)
    region_statistics: dict[str, Any] = {}
    for region, mask in masks.items():
        area_weights = areas[mask]
        current_weights = area_weights * j_gas_magnitude[mask]
        offsets = DETECTED_FRONT_REGIONS[region]
        region_statistics[region] = {
            "detected_front_offsets_c_over_omega_pi": (
                None if offsets is None else list(offsets)
            ),
            "R": _weighted_statistics(
                r_map[mask], area_weights, current_weights, label=f"{region} R"
            ),
            "Lambda": _weighted_statistics(
                lambda_map[mask],
                area_weights,
                current_weights,
                label=f"{region} Lambda",
            ),
        }

    shock_transition = np.abs(_x1_centers(state) - front_x) <= SHOCK_TRANSITION_HALF_WIDTH
    di_mask = np.broadcast_to((~shock_transition)[None, :], rho_g.shape)
    _require(np.any(di_mask), "shock-transition exclusion removed every cell")
    s_delta_min = float(np.min(s_delta_map[di_mask]))
    precursor_mask = masks[DI_MAGNETIC_SPECTRUM_REGION]
    precursor_di_max = float(np.max(di_map[precursor_mask]))
    spectrum = _leaf_aware_magnetic_spectrum(
        state, precursor_mask, di_max=precursor_di_max
    )
    particle = _particle_gyroradius(
        state,
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )
    particle_exposure = _particle_R_Lambda_exposure(
        state,
        front_x,
        r_map,
        lambda_map,
        points=points,
        cr_source=cr_source,
        birth_time=birth_time,
        velocity=velocity,
        macro_weight=macro_weight,
    )

    maximum_r = float(np.max(r_map))
    maximum_lambda = float(np.max(lambda_map))
    gates: dict[str, dict[str, Any]] = {
        "Q011-APP-NORM": {
            "pass": True,
            "threshold_provenance": "exact_equation_and_normalization_identity",
            "observed": "exact_match",
            "required": "exact_match",
        },
        "Q011-APP-R": {
            "pass": maximum_r <= R_MAXIMUM,
            "threshold_provenance": "AthenaK_selected_Bai_states_R_much_less_than_one",
            "observed_maximum": maximum_r,
            "required_maximum": R_MAXIMUM,
        },
        "Q011-APP-LAMBDA": {
            "pass": maximum_lambda <= LAMBDA_MAXIMUM,
            "threshold_provenance": (
                "AthenaK_selected_Bai_states_Lambda_much_less_than_one_for_no_Hall"
            ),
            "observed_maximum": maximum_lambda,
            "required_maximum": LAMBDA_MAXIMUM,
        },
        "Q011-APP-DI": {
            "pass": (
                s_delta_min >= S_DELTA_MINIMUM
                and spectrum["lambda_B_characteristic_over_local_di_maximum"]
                >= LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM
                and spectrum["sub_10di_magnetic_power_fraction_upper_bound"]
                <= SUB_10DI_POWER_FRACTION_MAXIMUM
                and spectrum["delta_B_rms_over_B0"]
                >= DELTA_B_RMS_OVER_B0_MINIMUM
            ),
            "threshold_provenance": (
                "AthenaK_selected_scale_separation_and_scattering_amplitude_"
                "floors_Bai_is_qualitative"
            ),
            "observed_S_delta_minimum_excluding_shock_transition": s_delta_min,
            "required_S_delta_minimum": S_DELTA_MINIMUM,
            "observed_lambda_B_characteristic_over_local_di_maximum": spectrum[
                "lambda_B_characteristic_over_local_di_maximum"
            ],
            "required_lambda_B_characteristic_over_local_di_maximum": (
                LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM
            ),
            "observed_sub_10di_magnetic_power_fraction_upper_bound": spectrum[
                "sub_10di_magnetic_power_fraction_upper_bound"
            ],
            "required_sub_10di_magnetic_power_fraction_upper_bound_maximum": (
                SUB_10DI_POWER_FRACTION_MAXIMUM
            ),
            "observed_delta_B_rms_over_B0": spectrum["delta_B_rms_over_B0"],
            "required_delta_B_rms_over_B0_minimum": DELTA_B_RMS_OVER_B0_MINIMUM,
        },
        "Q011-APP-RG": {
            "pass": (
                particle["macro_q999_over_Ly"] <= RG_Q999_OVER_LY_MAXIMUM
                and particle["energy_q999_over_Ly"] <= RG_Q999_OVER_LY_MAXIMUM
                and particle["energy_fraction_with_rg_above_Ly_over_4"]
                <= RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM
                and particle["maximum_over_Ly"] <= RG_MAXIMUM_OVER_LY_MAXIMUM
                and particle["high_energy_tail_maximum_over_Ly"]
                <= RG_HIGH_ENERGY_MAXIMUM_OVER_LY_MAXIMUM
            ),
            "threshold_provenance": (
                "AthenaK_selected_interpretation_of_Sun_and_Bai_several_gyroradii"
            ),
            "observed_maximum_over_Ly": particle["maximum_over_Ly"],
            "required_maximum_over_Ly": RG_MAXIMUM_OVER_LY_MAXIMUM,
            "observed_macro_q999_over_Ly": particle["macro_q999_over_Ly"],
            "observed_energy_q999_over_Ly": particle["energy_q999_over_Ly"],
            "required_q999_over_Ly_maximum": RG_Q999_OVER_LY_MAXIMUM,
            "observed_energy_fraction_with_rg_above_Ly_over_4": particle[
                "energy_fraction_with_rg_above_Ly_over_4"
            ],
            "required_energy_fraction_with_rg_above_Ly_over_4_maximum": (
                RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM
            ),
            "observed_high_energy_tail_maximum_over_Ly": particle[
                "high_energy_tail_maximum_over_Ly"
            ],
            "required_high_energy_tail_maximum_over_Ly": (
                RG_HIGH_ENERGY_MAXIMUM_OVER_LY_MAXIMUM
            ),
        },
        "Q011-APP-TIME": {
            "pass": False,
            "threshold_provenance": "engineering_and_analysis_completeness",
            "observed": "not_available_from_snapshot",
            "required": (
                "complete_per_cycle_post_startup_runtime_extrema_particle_exposure_"
                "and_boundary_escape_ledger_through_t1200"
            ),
        },
    }
    maps = MappingProxyType(
        {
            "R": _readonly(r_map),
            "Lambda": _readonly(lambda_map),
            "d_i": _readonly(di_map),
            "S_delta": _readonly(s_delta_map),
            "actual_leaf_dx1": _readonly(dx1_leaf),
            "actual_leaf_dx2": _readonly(dx2_leaf),
            "gas_frame_current_magnitude": _readonly(j_gas_magnitude),
        }
    )
    map_bindings = {
        name: {
            "shape_y_x": list(maps[name].shape),
            "dtype": "float64",
            "sha256": _array_sha256(maps[name]),
        }
        for name in CELL_MAP_NAMES
    }
    record: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "record_type": SNAPSHOT_RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "authorization": dict(AUTHORIZATION),
        "nominal_slot_time": state.nominal_slot_time,
        "observed_committed_time": state.observed_committed_time,
        "detected_front_x1_c_over_omega_pi": front_x,
        "exact_normalization": normalized,
        "bound_normalization_evidence": {
            "source_commit": bound_normalization["source_commit"],
            "runtime_input_parameters_sha256": bound_normalization[
                "runtime_input_parameters_sha256"
            ],
            "source_manifest": bound_normalization["source_manifest"],
            "bindings": bound_normalization["bindings"],
        },
        "snapshot_provenance": provenance,
        "formulae": {
            "R": "rho_q_over_rho_g_plus_rho_q",
            "v_A": "abs_B_over_sqrt_rho_g",
            "Lambda": "abs_J_q_minus_rho_q_v_g_over_rho_g_plus_rho_q_times_v_A",
            "d_i": "rho_g_to_the_minus_one_half",
            "S_delta": "minimum_actual_leaf_dx1_dx2_over_local_d_i",
        },
        "cell_map_contract": {
            "returned_separately_as_immutable_numpy_arrays": True,
            "retention_required_for_physical_applicability_evidence": True,
            "shape_y_x": list(rho_g.shape),
            "target_composite_level": state.target_level,
            "x1_faces_c_over_omega_pi": state.x1_faces.tolist(),
            "x2_faces_c_over_omega_pi": state.x2_faces.tolist(),
            "maps": list(CELL_MAP_NAMES),
            "map_bindings": map_bindings,
            "per_cell_maxima_required_no_rare_cell_waiver": True,
        },
        "regional_statistics": region_statistics,
        "ion_scale_separation": {
            "actual_leaf_spacing_used": True,
            "shock_transition_exclusion": {
                "center": "detected_front",
                "half_width_c_over_omega_pi": SHOCK_TRANSITION_HALF_WIDTH,
                "microscopic_shock_structure_claim_permanently_excluded": True,
                "self_consistent_injection_claim_permanently_excluded": True,
            },
            "S_delta_minimum_excluding_shock_transition": s_delta_min,
            "precursor_magnetic_spectrum": spectrum,
        },
        "particle_R_Lambda_exposure": particle_exposure,
        "particle_gyroradius_containment": particle,
        "gates": gates,
        "snapshot_gate_pass_excluding_time_completeness": all(
            gates[name]["pass"]
            for name in ("Q011-APP-NORM", "Q011-APP-R", "Q011-APP-LAMBDA", "Q011-APP-DI", "Q011-APP-RG")
        ),
        "claim_rejection_rules": {
            gate: list(claims) for gate, claims in CLAIM_REJECTION_RULES.items()
        },
        "permanent_claim_exclusions": list(PERMANENT_CLAIM_EXCLUSIONS),
        "claim_rejections": _claim_rejections(gates, maximum_lambda),
        "runtime_time_escape_evidence_required": True,
    }
    return ApplicabilitySnapshot(
        record=MappingProxyType(record),
        cell_maps=maps,
        evidence_root=evidence_root.resolve(strict=True),
    )


_CYCLE_EXTREMA_KEYS = (
    "R_maximum",
    "Lambda_maximum",
    "S_delta_minimum_excluding_shock_transition",
    "lambda_B_characteristic_over_local_di_maximum_minimum",
    "sub_10di_magnetic_power_fraction_upper_bound_maximum",
    "delta_B_rms_over_B0_minimum",
    "particle_rg_maximum_over_Ly",
    "high_energy_tail_rg_maximum_over_Ly",
    "maximum_particle_specific_kinetic_energy",
    "escaped_particle_rg_maximum_over_Ly",
    "escaped_high_energy_tail_rg_maximum_over_Ly",
    "escaped_particle_specific_kinetic_energy_maximum",
)
_CYCLE_MINIMUM_KEYS = {
    "S_delta_minimum_excluding_shock_transition",
    "lambda_B_characteristic_over_local_di_maximum_minimum",
    "delta_B_rms_over_B0_minimum",
}
_CYCLE_TELEMETRY_KEYS = set(_CYCLE_EXTREMA_KEYS) | {
    "schema_version",
    "record_type",
    "attempt_id",
    "source_commit",
    "executable_sha256",
    "runtime_normalization_sha256",
    "ps_escape_accounting_source_commit",
    "cycle",
    "previous_committed_cycle",
    "previous_committed_time",
    "start_time",
    "end_time",
}
_PARTICLE_EXPOSURE_KEYS = {
    "complete",
    "method",
    "active_particle_updates_included",
    "pre_destruction_boundary_events_included",
    "escaped_particles_included",
    "escaped_particles_included_in_high_energy_and_gyroradius_extrema",
    "observation_count",
    "maximum_sampled_R",
    "maximum_sampled_Lambda",
    "cumulative_macro_weighted_R_exceedance_fraction",
    "cumulative_CR_energy_weighted_R_exceedance_fraction",
    "cumulative_macro_weighted_Lambda_exceedance_fraction",
    "cumulative_CR_energy_weighted_Lambda_exceedance_fraction",
}
_STATE_VECTOR_KEYS = {"particle_count", "macro_weight", "kinetic_energy", "momentum"}
_FACE_ESCAPE_KEYS = _STATE_VECTOR_KEYS | {"reason_code"}
_ESCAPE_KEYS = {
    "complete",
    "scope",
    "nonperiodic_faces",
    "periodic_faces",
    "accumulated_escaped",
    "terminal_active",
    "startup_removed",
    "injected_particle_count",
    "injected_macro_weight",
    "escaped_particles_included_in_exposure",
    "escaped_particles_included_in_high_energy_and_gyroradius_extrema",
}
_PS_ESCAPE_REAL_LEDGER_KEYS = (
    "ps_escape_last_audit_time",
    "ps_escaped_injected_cr_count_global",
    "ps_escaped_injected_cr_mass_global",
    "ps_escaped_injected_cr_momentum_x1_global",
    "ps_escaped_injected_cr_momentum_x2_global",
    "ps_escaped_injected_cr_momentum_x3_global",
    "ps_escaped_injected_cr_energy_global",
    "ps_escaped_initial_cr_count_global",
    "ps_injected_cr_count_global",
    "ps_injected_cr_mass_global",
    "ps_injected_cr_momentum_x1_global",
    "ps_injected_cr_momentum_x2_global",
    "ps_injected_cr_momentum_x3_global",
    "ps_injected_cr_energy_global",
    "ps_removed_cr_count_global",
    "ps_removed_cr_mass_global",
    "ps_removed_cr_momentum_x1_global",
    "ps_removed_cr_momentum_x2_global",
    "ps_removed_cr_momentum_x3_global",
    "ps_removed_cr_energy_global",
)
_PS_ESCAPE_LEDGER_KEYS = set(_PS_ESCAPE_REAL_LEDGER_KEYS) | {
    "ps_cr_ledger_schema",
    "ps_cr_ledger_complete",
    "ps_removed_excluded_early_cohort",
    "ps_escape_ledger_schema",
    "ps_escape_ledger_complete",
    "ps_escape_audit_calls",
}
_PS_ESCAPE_CHECKPOINT_KEYS = {
    "nominal_checkpoint_time",
    "observed_committed_time",
    "cycle",
    "previous_committed_cycle",
    "previous_committed_time",
    "restart_artifact",
    "particle_checkpoint_artifact",
    "ps_escape_ledger",
    "active_injected_cr_count_global",
    "active_injected_cr_mass_global",
    "active_injected_cr_kinetic_energy_global",
    "escaped_injected_max_specific_kinetic_energy_global",
    "escaped_injected_max_rg_over_Ly_global",
}
_CLAIM_ESCAPE_ACCOUNTING = {
    "Emax_claim": "active_plus_escaped_all_cycle_maximum_and_preregistered_fraction_bounds",
    "high_energy_slope_or_cutoff_claim": (
        "bound_complete_binwise_and_high_energy_tail_active_plus_escaped_bias_bounds"
    ),
    "acceleration_rate_claim": (
        "active_plus_escaped_all_cycle_maximum_and_preregistered_fraction_bounds"
    ),
    "acceleration_efficiency_claim": (
        "escaped_energy_explicitly_bounded_by_preregistered_energy_fraction"
    ),
}
_RUNTIME_KEYS = {
    "schema_version",
    "record_type",
    "successor_id",
    "qualification_effect",
    "authorization",
    "attempt_id",
    "source_commit",
    "executable_sha256",
    "runtime_normalization_sha256",
    "ps_escape_accounting_source_commit",
    "per_cycle_inventory",
    "ps_escape_checkpoints",
    "claim_escape_accounting",
    "escaped_slope_cutoff_evidence",
    "particle_exposure",
    "boundary_escape_ledger",
}
_SLOPE_CUTOFF_ESCAPE_KEYS = {
    "schema_version",
    "record_type",
    "attempt_id",
    "source_commit",
    "executable_sha256",
    "runtime_normalization_sha256",
    "ps_escape_accounting_source_commit",
    "complete",
    "unavailable_reason",
    "high_energy_tail_threshold",
    "energy_bin_edges",
    "active_particle_count_by_bin",
    "active_macro_weight_by_bin",
    "active_kinetic_energy_by_bin",
    "escaped_particle_count_by_bin",
    "escaped_macro_weight_by_bin",
    "escaped_kinetic_energy_by_bin",
}



def _closure_residual(
    observed: float, expected: float, contributions: Sequence[float]
) -> float:
    scale = max(1.0, abs(observed), sum(abs(value) for value in contributions))
    return abs(observed - expected) / scale


def _decode_state_vector(value: object, label: str) -> dict[str, Any]:
    state = _exact_keys(value, _STATE_VECTOR_KEYS, label)
    momentum = state["momentum"]
    _require(
        type(momentum) is list and len(momentum) == 3,
        f"{label} momentum must be a three-component list",
    )
    return {
        "particle_count": _nonnegative_int(state["particle_count"], f"{label} count"),
        "macro_weight": _finite_scalar(
            state["macro_weight"], f"{label} macro weight", minimum=0.0
        ),
        "kinetic_energy": _finite_scalar(
            state["kinetic_energy"], f"{label} kinetic energy", minimum=0.0
        ),
        "momentum": [
            _finite_scalar(component, f"{label} momentum component {index}")
            for index, component in enumerate(momentum)
        ],
    }


def _decode_face_escape(
    value: object, label: str, *, expected_reason: str
) -> dict[str, Any]:
    face = _exact_keys(value, _FACE_ESCAPE_KEYS, label)
    _require(face["reason_code"] == expected_reason, f"{label} reason code drifted")
    return {
        "reason_code": expected_reason,
        **_decode_state_vector(
            {key: face[key] for key in _STATE_VECTOR_KEYS}, f"{label} state"
        ),
    }


def _validate_cycle_telemetry_inventory(
    evidence_root: Path,
    value: object,
    *,
    attempt_id: str,
    source_commit: str,
    executable_sha256: str,
    normalization_sha256: str,
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    _require(type(value) is list, "runtime per-cycle inventory must be a list")
    _require(
        len(value) >= MINIMUM_RUNTIME_CYCLE_INVENTORY_COUNT,
        "runtime per-cycle inventory is too short to establish complete coverage",
    )
    decoded: list[dict[str, Any]] = []
    digests: set[str] = set()
    paths: set[str] = set()
    for index, raw_binding in enumerate(value):
        binding, payload = _read_bound_artifact(
            evidence_root, raw_binding, expected_role="cycle_telemetry_record"
        )
        _require(
            binding["sha256"] not in digests and binding["path"] not in paths,
            "runtime per-cycle telemetry artifact was reused",
        )
        digests.add(str(binding["sha256"]))
        paths.add(str(binding["path"]))
        entry = _exact_keys(
            _decode_canonical_json(payload, f"runtime cycle telemetry {index}"),
            _CYCLE_TELEMETRY_KEYS,
            f"runtime cycle telemetry {index}",
        )
        _require(
            type(entry["schema_version"]) is int
            and entry["schema_version"] == SCHEMA_VERSION
            and entry["record_type"] == CYCLE_TELEMETRY_RECORD_TYPE,
            f"runtime cycle telemetry {index} identity drifted",
        )
        _require(
            entry["attempt_id"] == attempt_id
            and _source_commit(
                entry["source_commit"], f"runtime cycle telemetry {index} source commit"
            )
            == source_commit
            == TRUSTED_Q011_RUNTIME_SOURCE_COMMIT
            and _sha256_text(
                entry["executable_sha256"],
                f"runtime cycle telemetry {index} executable digest",
            )
            == executable_sha256
            and _sha256_text(
                entry["runtime_normalization_sha256"],
                f"runtime cycle telemetry {index} normalization digest",
            )
            == normalization_sha256
            and _source_commit(
                entry["ps_escape_accounting_source_commit"],
                f"runtime cycle telemetry {index} escape implementation commit",
            )
            == PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT,
            f"runtime cycle telemetry {index} trusted source/escape identity drifted",
        )
        cycle = _nonnegative_int(entry["cycle"], f"runtime cycle telemetry {index} cycle")
        previous_cycle = _nonnegative_int(
            entry["previous_committed_cycle"],
            f"runtime cycle telemetry {index} previous cycle",
        )
        previous_time = _finite_scalar(
            entry["previous_committed_time"],
            f"runtime cycle telemetry {index} previous time",
            minimum=0.0,
        )
        start = _finite_scalar(
            entry["start_time"], f"runtime cycle telemetry {index} start time", minimum=0.0
        )
        end = _finite_scalar(
            entry["end_time"], f"runtime cycle telemetry {index} end time", minimum=0.0
        )
        _require(
            previous_cycle + 1 == cycle and previous_time <= start < end,
            f"runtime cycle telemetry {index} chronology drifted",
        )
        _require(
            end - start <= MAXIMUM_RUNTIME_CYCLE_SPAN
            and start - previous_time <= MAXIMUM_RUNTIME_CYCLE_SPAN,
            f"runtime cycle telemetry {index} exceeds the maximum chronology span",
        )
        metrics = {
            key: _finite_scalar(
                entry[key], f"runtime cycle telemetry {index} {key}", minimum=0.0
            )
            for key in _CYCLE_EXTREMA_KEYS
        }
        _require(
            metrics["particle_rg_maximum_over_Ly"]
            >= metrics["escaped_particle_rg_maximum_over_Ly"],
            "runtime particle gyroradius maximum omits escaped particles",
        )
        _require(
            metrics["high_energy_tail_rg_maximum_over_Ly"]
            >= metrics["escaped_high_energy_tail_rg_maximum_over_Ly"],
            "runtime high-energy gyroradius maximum omits escaped particles",
        )
        _require(
            metrics["maximum_particle_specific_kinetic_energy"]
            >= metrics["escaped_particle_specific_kinetic_energy_maximum"],
            "runtime maximum particle energy omits escaped particles",
        )
        decoded.append(
            {
                "binding": binding,
                "cycle": cycle,
                "previous_committed_cycle": previous_cycle,
                "previous_committed_time": previous_time,
                "start_time": start,
                "end_time": end,
                **metrics,
            }
        )
    first = decoded[0]
    _require(
        first["previous_committed_time"] < STARTUP_REMOVAL_TIME
        <= first["start_time"]
        and first["previous_committed_cycle"] + 1 == first["cycle"],
        "runtime first retained cycle is not the first committed cycle starting across t=45",
    )
    _require(
        decoded[-1]["end_time"] == EXPECTED_TERMINAL_TIME,
        "runtime per-cycle inventory must bind exact actual terminal endpoint t=1200",
    )
    for left, right in zip(decoded, decoded[1:]):
        _require(
            right["cycle"] == left["cycle"] + 1
            and right["previous_committed_cycle"] == left["cycle"],
            "runtime per-cycle inventory has a cycle gap or previous-cycle drift",
        )
        _require(
            right["previous_committed_time"] == left["end_time"]
            and right["start_time"] == left["end_time"],
            "runtime per-cycle inventory has a time gap, overlap, or previous-time drift",
        )
    extrema = {
        key: (
            min(float(entry[key]) for entry in decoded)
            if key in _CYCLE_MINIMUM_KEYS
            else max(float(entry[key]) for entry in decoded)
        )
        for key in _CYCLE_EXTREMA_KEYS
    }
    coverage = {
        "sampling_mode": "opened_immutable_byte_bound_every_integrator_cycle_telemetry",
        "previous_committed_cycle_before_startup_crossing": first[
            "previous_committed_cycle"
        ],
        "previous_committed_time_before_startup_crossing": first[
            "previous_committed_time"
        ],
        "post_startup_removal_start_time": first["start_time"],
        "terminal_time": decoded[-1]["end_time"],
        "first_cycle": first["cycle"],
        "last_cycle": decoded[-1]["cycle"],
        "covered_cycle_count": len(decoded),
        "maximum_cycle_span": max(entry["end_time"] - entry["start_time"] for entry in decoded),
        "complete_contiguous_actual_endpoints": True,
        "trusted_source_commit": TRUSTED_Q011_RUNTIME_SOURCE_COMMIT,
        "trusted_escape_implementation_commit": PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT,
    }
    return decoded, extrema, coverage


def _nonnegative_numeric_list(
    value: object, label: str, *, integer: bool = False
) -> list[float | int]:
    _require(type(value) is list, f"{label} must be a list")
    if integer:
        return [_nonnegative_int(item, f"{label} item {index}") for index, item in enumerate(value)]
    return [
        _finite_scalar(item, f"{label} item {index}", minimum=0.0)
        for index, item in enumerate(value)
    ]


def _validate_slope_cutoff_escape_evidence(
    evidence_root: Path,
    value: object,
    *,
    attempt_id: str,
    source_commit: str,
    executable_sha256: str,
    normalization_sha256: str,
    active: Mapping[str, Any],
    escaped: Mapping[str, Any],
    maximum_specific_energy: float,
) -> dict[str, Any]:
    binding, payload = _read_bound_artifact(
        evidence_root, value, expected_role="slope_cutoff_escape_evidence"
    )
    record = _exact_keys(
        _decode_canonical_json(payload, "slope/cutoff escape evidence"),
        _SLOPE_CUTOFF_ESCAPE_KEYS,
        "slope/cutoff escape evidence",
    )
    _require(
        type(record["schema_version"]) is int
        and record["schema_version"] == SCHEMA_VERSION
        and record["record_type"] == SLOPE_CUTOFF_ESCAPE_RECORD_TYPE,
        "slope/cutoff escape evidence identity drifted",
    )
    _require(
        record["attempt_id"] == attempt_id
        and _source_commit(record["source_commit"], "slope/cutoff source commit")
        == source_commit
        == TRUSTED_Q011_RUNTIME_SOURCE_COMMIT
        and _sha256_text(record["executable_sha256"], "slope/cutoff executable digest")
        == executable_sha256
        and _sha256_text(
            record["runtime_normalization_sha256"], "slope/cutoff normalization digest"
        )
        == normalization_sha256
        and _source_commit(
            record["ps_escape_accounting_source_commit"],
            "slope/cutoff escape implementation commit",
        )
        == PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT,
        "slope/cutoff trusted source/escape identity drifted",
    )
    complete = _strict_bool(record["complete"], "slope/cutoff evidence complete")
    array_keys = (
        "energy_bin_edges",
        "active_particle_count_by_bin",
        "active_macro_weight_by_bin",
        "active_kinetic_energy_by_bin",
        "escaped_particle_count_by_bin",
        "escaped_macro_weight_by_bin",
        "escaped_kinetic_energy_by_bin",
    )
    if not complete:
        _require(
            type(record["unavailable_reason"]) is str
            and bool(record["unavailable_reason"])
            and record["high_energy_tail_threshold"] is None
            and all(record[key] == [] for key in array_keys),
            "incomplete slope/cutoff evidence must explicitly bind an unavailable reason and no bins",
        )
        return {
            "binding": binding,
            "complete": False,
            "pass": False,
            "unavailable_reason": record["unavailable_reason"],
            "AthenaK_selected_bin_escape_fraction_maximum": (
                SLOPE_CUTOFF_BIN_ESCAPE_FRACTION_MAXIMUM
            ),
            "AthenaK_selected_tail_escape_fraction_maximum": (
                SLOPE_CUTOFF_TAIL_ESCAPE_FRACTION_MAXIMUM
            ),
        }
    _require(record["unavailable_reason"] is None, "complete slope/cutoff evidence has a reason")
    threshold = _finite_scalar(
        record["high_energy_tail_threshold"], "slope/cutoff high-energy threshold", minimum=0.0
    )
    edges = _nonnegative_numeric_list(record["energy_bin_edges"], "slope/cutoff bin edges")
    _require(
        len(edges) >= MINIMUM_SLOPE_CUTOFF_HIGH_ENERGY_BINS + 1
        and all(right > left for left, right in zip(edges, edges[1:]))
        and threshold in edges[:-1]
        and edges[-1] >= maximum_specific_energy,
        "slope/cutoff energy-bin coverage is insufficient or non-monotonic",
    )
    bin_count = len(edges) - 1
    arrays = {
        "active_particle_count_by_bin": _nonnegative_numeric_list(
            record["active_particle_count_by_bin"],
            "slope/cutoff active particle counts",
            integer=True,
        ),
        "active_macro_weight_by_bin": _nonnegative_numeric_list(
            record["active_macro_weight_by_bin"], "slope/cutoff active macro weights"
        ),
        "active_kinetic_energy_by_bin": _nonnegative_numeric_list(
            record["active_kinetic_energy_by_bin"], "slope/cutoff active kinetic energies"
        ),
        "escaped_particle_count_by_bin": _nonnegative_numeric_list(
            record["escaped_particle_count_by_bin"],
            "slope/cutoff escaped particle counts",
            integer=True,
        ),
        "escaped_macro_weight_by_bin": _nonnegative_numeric_list(
            record["escaped_macro_weight_by_bin"], "slope/cutoff escaped macro weights"
        ),
        "escaped_kinetic_energy_by_bin": _nonnegative_numeric_list(
            record["escaped_kinetic_energy_by_bin"], "slope/cutoff escaped kinetic energies"
        ),
    }
    _require(
        all(len(values) == bin_count for values in arrays.values()),
        "slope/cutoff bin array lengths drifted",
    )
    high_bins = [index for index, lower in enumerate(edges[:-1]) if lower >= threshold]
    _require(
        len(high_bins) >= MINIMUM_SLOPE_CUTOFF_HIGH_ENERGY_BINS
        and all(
            arrays["active_particle_count_by_bin"][index]
            + arrays["escaped_particle_count_by_bin"][index]
            > 0
            and arrays["active_macro_weight_by_bin"][index]
            + arrays["escaped_macro_weight_by_bin"][index]
            > 0.0
            and arrays["active_kinetic_energy_by_bin"][index]
            + arrays["escaped_kinetic_energy_by_bin"][index]
            > 0.0
            for index in high_bins
        ),
        "slope/cutoff high-energy binwise evidence is insufficient",
    )
    _require(
        sum(arrays["active_particle_count_by_bin"]) == active["particle_count"]
        and sum(arrays["escaped_particle_count_by_bin"]) == escaped["particle_count"],
        "slope/cutoff particle-count bins disagree with active/escaped ledgers",
    )
    for key, expected, label in (
        ("active_macro_weight_by_bin", active["macro_weight"], "active macro weight"),
        ("active_kinetic_energy_by_bin", active["kinetic_energy"], "active kinetic energy"),
        ("escaped_macro_weight_by_bin", escaped["macro_weight"], "escaped macro weight"),
        ("escaped_kinetic_energy_by_bin", escaped["kinetic_energy"], "escaped kinetic energy"),
    ):
        _require(
            _closure_residual(float(expected), float(sum(arrays[key])), [float(v) for v in arrays[key]])
            <= ESCAPE_LEDGER_RELATIVE_RESIDUAL_MAXIMUM,
            f"slope/cutoff {label} bins disagree with ledger",
        )

    def fractions(active_key: str, escaped_key: str) -> list[float]:
        result = []
        for active_value, escaped_value in zip(arrays[active_key], arrays[escaped_key]):
            denominator = float(active_value) + float(escaped_value)
            result.append(float(escaped_value) / denominator if denominator > 0.0 else 0.0)
        return result

    bin_fractions = {
        "particle_count": fractions(
            "active_particle_count_by_bin", "escaped_particle_count_by_bin"
        ),
        "macro_weight": fractions(
            "active_macro_weight_by_bin", "escaped_macro_weight_by_bin"
        ),
        "kinetic_energy": fractions(
            "active_kinetic_energy_by_bin", "escaped_kinetic_energy_by_bin"
        ),
    }
    tail_fractions: dict[str, float] = {}
    for name, active_key, escaped_key in (
        ("particle_count", "active_particle_count_by_bin", "escaped_particle_count_by_bin"),
        ("macro_weight", "active_macro_weight_by_bin", "escaped_macro_weight_by_bin"),
        ("kinetic_energy", "active_kinetic_energy_by_bin", "escaped_kinetic_energy_by_bin"),
    ):
        active_tail = sum(float(arrays[active_key][index]) for index in high_bins)
        escaped_tail = sum(float(arrays[escaped_key][index]) for index in high_bins)
        tail_fractions[name] = (
            escaped_tail / (active_tail + escaped_tail)
            if active_tail + escaped_tail > 0.0
            else 0.0
        )
    passes = all(
        bin_fractions[name][index] <= SLOPE_CUTOFF_BIN_ESCAPE_FRACTION_MAXIMUM
        for name in bin_fractions
        for index in high_bins
    ) and all(
        value <= SLOPE_CUTOFF_TAIL_ESCAPE_FRACTION_MAXIMUM
        for value in tail_fractions.values()
    )
    return {
        "binding": binding,
        "complete": True,
        "pass": passes,
        "unavailable_reason": None,
        "high_energy_tail_threshold": threshold,
        "energy_bin_edges": edges,
        **arrays,
        "high_energy_bin_indices": high_bins,
        "escaped_fraction_by_bin": bin_fractions,
        "escaped_high_energy_tail_fractions": tail_fractions,
        "AthenaK_selected_bin_escape_fraction_maximum": (
            SLOPE_CUTOFF_BIN_ESCAPE_FRACTION_MAXIMUM
        ),
        "AthenaK_selected_tail_escape_fraction_maximum": (
            SLOPE_CUTOFF_TAIL_ESCAPE_FRACTION_MAXIMUM
        ),
    }


def _integer_valued_real(value: object, label: str) -> float:
    decoded = _finite_scalar(value, label, minimum=0.0)
    _require(decoded == math.floor(decoded), f"{label} must be integer-valued")
    return decoded


def _decode_ps_escape_ledger(value: object, label: str) -> dict[str, Any]:
    ledger = _exact_keys(value, _PS_ESCAPE_LEDGER_KEYS, label)
    _require(
        type(ledger["ps_cr_ledger_schema"]) is int
        and ledger["ps_cr_ledger_schema"] == PS_CR_LEDGER_SCHEMA,
        f"{label} must bind complete schema-3 CR source accounting",
    )
    _require(
        type(ledger["ps_escape_ledger_schema"]) is int
        and ledger["ps_escape_ledger_schema"] == PS_ESCAPE_LEDGER_SCHEMA,
        f"{label} must be complete schema-1 ps_escape",
    )
    _require(
        _strict_bool(ledger["ps_cr_ledger_complete"], f"{label} CR ledger complete"),
        f"{label} CR source accounting is incomplete",
    )
    _require(
        _strict_bool(ledger["ps_escape_ledger_complete"], f"{label} complete"),
        f"{label} is incomplete",
    )
    decoded: dict[str, Any] = {
        "ps_cr_ledger_schema": PS_CR_LEDGER_SCHEMA,
        "ps_cr_ledger_complete": True,
        "ps_removed_excluded_early_cohort": _strict_bool(
            ledger["ps_removed_excluded_early_cohort"],
            f"{label} startup-removal completion",
        ),
        "ps_escape_ledger_schema": PS_ESCAPE_LEDGER_SCHEMA,
        "ps_escape_ledger_complete": True,
        "ps_escape_audit_calls": _nonnegative_int(
            ledger["ps_escape_audit_calls"], f"{label} audit calls"
        ),
    }
    count_fields = {
        "ps_escaped_injected_cr_count_global",
        "ps_escaped_initial_cr_count_global",
        "ps_injected_cr_count_global",
        "ps_removed_cr_count_global",
    }
    nonnegative_fields = count_fields | {
        "ps_escape_last_audit_time",
        "ps_escaped_injected_cr_mass_global",
        "ps_escaped_injected_cr_energy_global",
        "ps_injected_cr_mass_global",
        "ps_injected_cr_energy_global",
        "ps_removed_cr_mass_global",
        "ps_removed_cr_energy_global",
    }
    for key in _PS_ESCAPE_REAL_LEDGER_KEYS:
        if key in count_fields:
            decoded[key] = _integer_valued_real(ledger[key], f"{label} {key}")
        else:
            decoded[key] = _finite_scalar(
                ledger[key],
                f"{label} {key}",
                minimum=0.0 if key in nonnegative_fields else None,
            )
    _require(
        decoded["ps_escaped_initial_cr_count_global"] == 0.0,
        f"{label} production initial-CR escape must be exactly zero",
    )
    return decoded


def _ps_escape_restart_ledger(payload: bytes, label: str) -> dict[str, Any]:
    parameters = _restart_problem_parameters(payload, label)
    _require(
        _PS_ESCAPE_LEDGER_KEYS <= set(parameters),
        f"{label} lacks complete schema-1 ps_escape ledger fields",
    )
    try:
        decoded: dict[str, Any] = {
            "ps_cr_ledger_schema": int(parameters["ps_cr_ledger_schema"]),
            "ps_cr_ledger_complete": _restart_bool(
                parameters["ps_cr_ledger_complete"],
                f"{label} ps_cr_ledger_complete",
            ),
            "ps_removed_excluded_early_cohort": _restart_bool(
                parameters["ps_removed_excluded_early_cohort"],
                f"{label} ps_removed_excluded_early_cohort",
            ),
            "ps_escape_ledger_schema": int(parameters["ps_escape_ledger_schema"]),
            "ps_escape_ledger_complete": _restart_bool(
                parameters["ps_escape_ledger_complete"],
                f"{label} ps_escape_ledger_complete",
            ),
            "ps_escape_audit_calls": int(parameters["ps_escape_audit_calls"]),
            **{
                key: float(parameters[key])
                for key in _PS_ESCAPE_REAL_LEDGER_KEYS
            },
        }
    except ValueError as error:
        raise PhysicalApplicabilityError(
            f"{label} contains invalid ps_escape numeric metadata"
        ) from error
    return _decode_ps_escape_ledger(decoded, label)


def _decode_bound_particle_checkpoint(
    payload: bytes,
    label: str,
    *,
    observed_time: float,
    cycle: int,
    ledger: Mapping[str, Any],
) -> dict[str, Any]:
    match = _PVTK_EXECUTION_PATTERN.match(payload[:4096])
    _require(match is not None, f"{label} lacks canonical prtcl_all execution metadata")
    try:
        header_time = float(match.group(1))
        header_nranks = int(match.group(2))
        header_cycle = int(match.group(3))
        header_variables = match.group(4).decode("ascii")
    except (UnicodeDecodeError, ValueError) as error:
        raise PhysicalApplicabilityError(
            f"{label} contains invalid prtcl_all execution metadata"
        ) from error
    _require(
        math.isfinite(header_time)
        and header_nranks > 0
        and header_cycle == cycle
        and header_time == observed_time
        and header_variables == "prtcl_all",
        f"{label} prtcl_all cycle/time binding drifted",
    )
    try:
        descriptor = os.memfd_create(
            "q011-ps-escape-particle-checkpoint",
            flags=getattr(os, "MFD_CLOEXEC", 0),
        )
        try:
            remaining = memoryview(payload)
            while remaining:
                written = os.write(descriptor, remaining)
                _require(written > 0, f"{label} could not stage bound particle bytes")
                remaining = remaining[written:]
            particles = pvtk_particles.read_particle_vtk(
                Path(f"/proc/self/fd/{descriptor}")
            )
        finally:
            os.close(descriptor)
    except (OSError, ValueError) as error:
        raise PhysicalApplicabilityError(
            f"{label} bound prtcl_all payload cannot be decoded"
        ) from error
    _require(
        set(particles.scalars) == _PVTK_SCALARS and set(particles.vectors) == {"vel"},
        f"{label} prtcl_all field inventory drifted",
    )
    particle_count = particles.points.shape[0]
    _require(
        particles.points.shape == (particle_count, 3)
        and particles.vectors["vel"].shape == (particle_count, 3)
        and particle_count > 0,
        f"{label} prtcl_all shape or population is invalid",
    )
    for name in _PVTK_SCALARS:
        _require(
            particles.scalars[name].shape == (particle_count,),
            f"{label} prtcl_all scalar shape drifted",
        )
    _require(
        np.all(particles.scalars["gid"] >= 0)
        and np.all(particles.scalars["ptag"] >= 0)
        and np.unique(particles.scalars["ptag"]).size == particle_count
        and np.all(np.isin(particles.scalars["cr_source"], (0, 1)))
        and np.all(particles.scalars["macro_weight"] >= 0.0),
        f"{label} prtcl_all provenance or macro weight is invalid",
    )
    arrays = science._decoded_particle_arrays(
        points=particles.points,
        cr_source=particles.scalars["cr_source"],
        birth_time=particles.scalars["birth_time"],
        velocity=particles.vectors["vel"],
        macro_weight=particles.scalars["macro_weight"],
    )
    active = arrays["active"]
    injected = particles.scalars["cr_source"] == science.SHOCK_INJECTED_SOURCE
    _require(
        np.all(particles.scalars["species"][injected] == EXACT_NORMALIZATION["selected_species_index"])
        and np.all(particles.scalars["birth_time"][injected] >= STARTUP_REMOVAL_TIME)
        and np.all(particles.scalars["macro_weight"][injected] == 1.0),
        f"{label} retains invalid post-startup shock-injected particles",
    )
    active_count = int(np.count_nonzero(active))
    _require(active_count > 0, f"{label} contains no active injected particles")
    injected_count = ledger["ps_injected_cr_count_global"]
    _require(injected_count > 0.0, f"{label} sealed injected source is empty")
    particle_macro_mass = ledger["ps_injected_cr_mass_global"] / injected_count
    _require(
        math.isfinite(particle_macro_mass) and particle_macro_mass > 0.0,
        f"{label} sealed particle macro mass is invalid",
    )
    mass_residuals = {
        "injected": _closure_residual(
            ledger["ps_injected_cr_mass_global"],
            injected_count * particle_macro_mass,
            [injected_count * particle_macro_mass],
        ),
        "startup_removed": _closure_residual(
            ledger["ps_removed_cr_mass_global"],
            ledger["ps_removed_cr_count_global"] * particle_macro_mass,
            [ledger["ps_removed_cr_count_global"] * particle_macro_mass],
        ),
        "escaped": _closure_residual(
            ledger["ps_escaped_injected_cr_mass_global"],
            ledger["ps_escaped_injected_cr_count_global"] * particle_macro_mass,
            [ledger["ps_escaped_injected_cr_count_global"] * particle_macro_mass],
        ),
    }
    _require(
        max(mass_residuals.values()) <= ESCAPE_LEDGER_RELATIVE_RESIDUAL_MAXIMUM,
        f"{label} sealed count-to-mass accounting does not close",
    )
    active_mass = active_count * particle_macro_mass
    active_energy = particle_macro_mass * float(
        np.sum(arrays["specific_kinetic_energy"][active])
    )
    active_energy_lower = particle_macro_mass * float(
        np.sum(arrays["specific_kinetic_energy_lower"][active])
    )
    active_energy_upper = particle_macro_mass * float(
        np.sum(arrays["specific_kinetic_energy_upper"][active])
    )
    return {
        "execution_header": {
            "observed_committed_time": header_time,
            "nranks": header_nranks,
            "cycle": header_cycle,
            "variables": header_variables,
        },
        "all_particle_count": particle_count,
        "active_injected_cr_count_global": active_count,
        "active_injected_cr_mass_global": active_mass,
        "active_injected_cr_kinetic_energy_global": active_energy,
        "active_injected_cr_kinetic_energy_lower_bound": active_energy_lower,
        "active_injected_cr_kinetic_energy_upper_bound": active_energy_upper,
        "active_injected_cr_max_specific_kinetic_energy_upper_bound": float(
            np.max(arrays["specific_kinetic_energy_upper"][active])
        ),
        "particle_macro_mass": particle_macro_mass,
        "sealed_count_to_mass_relative_residuals": mass_residuals,
    }


def _validate_ps_escape_checkpoints(
    evidence_root: Path,
    value: object,
    *,
    inventory: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    _require(type(value) is list, "ps_escape checkpoints must be a list")
    _require(
        len(value) == len(REQUIRED_PS_ESCAPE_CHECKPOINT_NOMINAL_TIMES),
        "ps_escape checkpoint cadence is incomplete",
    )
    decoded: list[dict[str, Any]] = []
    for index, (raw_checkpoint, required_nominal) in enumerate(
        zip(value, REQUIRED_PS_ESCAPE_CHECKPOINT_NOMINAL_TIMES)
    ):
        checkpoint = _exact_keys(
            raw_checkpoint, _PS_ESCAPE_CHECKPOINT_KEYS, f"ps_escape checkpoint {index}"
        )
        nominal = _finite_scalar(
            checkpoint["nominal_checkpoint_time"],
            f"ps_escape checkpoint {index} nominal time",
            minimum=0.0,
        )
        _require(nominal == required_nominal, "ps_escape checkpoint cadence drifted")
        observed = _finite_scalar(
            checkpoint["observed_committed_time"],
            f"ps_escape checkpoint {index} observed time",
            minimum=0.0,
        )
        cycle = _nonnegative_int(checkpoint["cycle"], f"ps_escape checkpoint {index} cycle")
        previous_cycle = _nonnegative_int(
            checkpoint["previous_committed_cycle"],
            f"ps_escape checkpoint {index} previous committed cycle",
        )
        previous_time = _finite_scalar(
            checkpoint["previous_committed_time"],
            f"ps_escape checkpoint {index} previous committed time",
            minimum=0.0,
        )
        crossings = [entry for entry in inventory if entry["end_time"] >= nominal]
        _require(bool(crossings), "ps_escape checkpoint nominal slot was never crossed")
        first_crossing = crossings[0]
        _require(
            first_crossing["previous_committed_time"] < nominal
            <= first_crossing["end_time"]
            and first_crossing["cycle"] == cycle
            and first_crossing["end_time"] == observed
            and first_crossing["previous_committed_cycle"] == previous_cycle
            and first_crossing["previous_committed_time"] == previous_time,
            "ps_escape checkpoint is not the first committed cycle crossing its nominal slot",
        )
        if nominal == EXPECTED_TERMINAL_TIME:
            _require(observed == nominal, "ps_escape endpoint checkpoint time drifted")
        restart_binding, restart_payload = _read_bound_artifact(
            evidence_root,
            checkpoint["restart_artifact"],
            expected_role="ps_escape_restart_checkpoint",
        )
        particle_binding, particle_payload = _read_bound_artifact(
            evidence_root,
            checkpoint["particle_checkpoint_artifact"],
            expected_role="ps_escape_particle_checkpoint",
        )
        recorded_ledger = _decode_ps_escape_ledger(
            checkpoint["ps_escape_ledger"], f"ps_escape checkpoint {index} recorded ledger"
        )
        parsed_ledger = _ps_escape_restart_ledger(
            restart_payload, f"ps_escape checkpoint {index} restart"
        )
        _require(
            _canonical_json_bytes(recorded_ledger) == _canonical_json_bytes(parsed_ledger),
            "ps_escape checkpoint ledger disagrees with bound restart header",
        )
        _require(
            parsed_ledger["ps_escape_audit_calls"]
            == PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE * cycle,
            "ps_escape checkpoint audit-call coverage is incomplete",
        )
        _require(
            parsed_ledger["ps_escape_last_audit_time"] == observed,
            "ps_escape checkpoint last audit time disagrees with actual checkpoint time",
        )
        _require(
            parsed_ledger["ps_removed_excluded_early_cohort"],
            "ps_escape checkpoint predates completed startup cohort removal",
        )
        matches = [first_crossing]
        particle_inventory = _decode_bound_particle_checkpoint(
            particle_payload,
            f"ps_escape checkpoint {index} particle artifact",
            observed_time=observed,
            cycle=cycle,
            ledger=parsed_ledger,
        )
        active_count = _integer_valued_real(
            checkpoint["active_injected_cr_count_global"],
            f"ps_escape checkpoint {index} active count",
        )
        active_mass = _finite_scalar(
            checkpoint["active_injected_cr_mass_global"],
            f"ps_escape checkpoint {index} active mass",
            minimum=0.0,
        )
        active_energy = _finite_scalar(
            checkpoint["active_injected_cr_kinetic_energy_global"],
            f"ps_escape checkpoint {index} active energy",
            minimum=0.0,
        )
        escaped_max_energy = _finite_scalar(
            checkpoint["escaped_injected_max_specific_kinetic_energy_global"],
            f"ps_escape checkpoint {index} escaped maximum energy",
            minimum=0.0,
        )
        escaped_max_rg = _finite_scalar(
            checkpoint["escaped_injected_max_rg_over_Ly_global"],
            f"ps_escape checkpoint {index} escaped maximum gyroradius",
            minimum=0.0,
        )
        _require(
            active_count == particle_inventory["active_injected_cr_count_global"],
            "ps_escape active count disagrees with bound particle checkpoint",
        )
        _require(
            _closure_residual(
                active_mass,
                particle_inventory["active_injected_cr_mass_global"],
                [particle_inventory["active_injected_cr_mass_global"]],
            )
            <= ESCAPE_LEDGER_RELATIVE_RESIDUAL_MAXIMUM,
            "ps_escape active mass disagrees with bound particle checkpoint",
        )
        _require(
            particle_inventory["active_injected_cr_kinetic_energy_lower_bound"]
            <= active_energy
            <= particle_inventory["active_injected_cr_kinetic_energy_upper_bound"],
            "ps_escape active energy disagrees with bound particle checkpoint",
        )
        _require(
            matches[0]["maximum_particle_specific_kinetic_energy"]
            >= particle_inventory[
                "active_injected_cr_max_specific_kinetic_energy_upper_bound"
            ],
            "ps_escape per-cycle energy maximum omits bound active particles",
        )
        cumulative_inventory = [entry for entry in inventory if entry["cycle"] <= cycle]
        _require(
            escaped_max_energy
            == max(
                entry["escaped_particle_specific_kinetic_energy_maximum"]
                for entry in cumulative_inventory
            )
            and escaped_max_rg
            == max(
                entry["escaped_particle_rg_maximum_over_Ly"]
                for entry in cumulative_inventory
            ),
            "ps_escape checkpoint escaped maxima disagree with bound per-cycle inventory",
        )
        _require(
            active_count
            + parsed_ledger["ps_removed_cr_count_global"]
            + parsed_ledger["ps_escaped_injected_cr_count_global"]
            == parsed_ledger["ps_injected_cr_count_global"],
            "ps_escape active + startup_removed + escaped_injected count does not equal injected",
        )
        mass_residual = _closure_residual(
            parsed_ledger["ps_injected_cr_mass_global"],
            active_mass
            + parsed_ledger["ps_removed_cr_mass_global"]
            + parsed_ledger["ps_escaped_injected_cr_mass_global"],
            [
                active_mass,
                parsed_ledger["ps_removed_cr_mass_global"],
                parsed_ledger["ps_escaped_injected_cr_mass_global"],
            ],
        )
        _require(
            mass_residual <= ESCAPE_LEDGER_RELATIVE_RESIDUAL_MAXIMUM,
            "ps_escape active + startup_removed + escaped_injected mass does not equal injected",
        )
        escaped_count = parsed_ledger["ps_escaped_injected_cr_count_global"]
        if escaped_count == 0.0:
            _require(
                parsed_ledger["ps_escaped_injected_cr_mass_global"] == 0.0
                and parsed_ledger["ps_escaped_injected_cr_energy_global"] == 0.0
                and escaped_max_energy == 0.0
                and escaped_max_rg == 0.0,
                "zero ps_escape count requires zero escaped mass, energy, and maxima",
            )
        decoded.append(
            {
                "nominal_checkpoint_time": nominal,
                "observed_committed_time": observed,
                "cycle": cycle,
                "previous_committed_cycle": previous_cycle,
                "previous_committed_time": previous_time,
                "restart_artifact": restart_binding,
                "particle_checkpoint_artifact": particle_binding,
                "ps_escape_ledger": parsed_ledger,
                "active_injected_cr_count_global": active_count,
                "active_injected_cr_mass_global": active_mass,
                "active_injected_cr_kinetic_energy_global": active_energy,
                "active_particle_inventory": particle_inventory,
                "escaped_injected_max_specific_kinetic_energy_global": escaped_max_energy,
                "escaped_injected_max_rg_over_Ly_global": escaped_max_rg,
                "particle_checkpoint_sha256": particle_binding["sha256"],
                "source_mass_relative_residual": mass_residual,
            }
        )
    monotonic_keys = (
        "ps_escape_audit_calls",
        "ps_escape_last_audit_time",
        "ps_escaped_injected_cr_count_global",
        "ps_escaped_injected_cr_mass_global",
        "ps_escaped_injected_cr_energy_global",
        "ps_injected_cr_count_global",
        "ps_injected_cr_mass_global",
        "ps_injected_cr_energy_global",
        "ps_removed_cr_count_global",
        "ps_removed_cr_mass_global",
        "ps_removed_cr_energy_global",
    )
    for left, right in zip(decoded, decoded[1:]):
        _require(
            right["cycle"] > left["cycle"]
            and right["observed_committed_time"] > left["observed_committed_time"],
            "ps_escape checkpoint cycle/time coverage is not strictly increasing",
        )
        for key in monotonic_keys:
            _require(
                right["ps_escape_ledger"][key] >= left["ps_escape_ledger"][key],
                f"ps_escape cumulative ledger {key} is not monotonic",
            )
        _require(
            right["escaped_injected_max_specific_kinetic_energy_global"]
            >= left["escaped_injected_max_specific_kinetic_energy_global"]
            and right["escaped_injected_max_rg_over_Ly_global"]
            >= left["escaped_injected_max_rg_over_Ly_global"],
            "ps_escape cumulative escaped-particle maxima are not monotonic",
        )
    return {
        "source_commit": PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT,
        "schema": PS_ESCAPE_LEDGER_SCHEMA,
        "required_nominal_checkpoint_times": list(
            REQUIRED_PS_ESCAPE_CHECKPOINT_NOMINAL_TIMES
        ),
        "paper_vl2_escape_audits_per_cycle": PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE,
        "checkpoints": decoded,
        "terminal": decoded[-1],
        "complete_checkpoint_cadence": True,
    }


def _validate_runtime_time_escape(
    evidence_root: Path, value: object
) -> dict[str, Any]:
    envelope = _exact_keys(
        value, {"runtime_time_escape_record"}, "runtime time/escape evidence envelope"
    )
    binding, payload = _read_bound_artifact(
        evidence_root,
        envelope["runtime_time_escape_record"],
        expected_role="runtime_time_escape_record",
    )
    runtime = _decode_canonical_json(payload, "runtime time/escape record")
    runtime = _exact_keys(runtime, _RUNTIME_KEYS, "runtime time/escape record")
    _require(
        type(runtime["schema_version"]) is int
        and runtime["schema_version"] == SCHEMA_VERSION,
        "runtime schema version drifted",
    )
    _require(runtime["record_type"] == RUNTIME_RECORD_TYPE, "runtime record type drifted")
    _require(runtime["successor_id"] == SUCCESSOR_ID, "runtime successor id drifted")
    _require(
        runtime["qualification_effect"] == QUALIFICATION_EFFECT,
        "runtime effect drifted",
    )
    authorization = _validate_authorization(runtime["authorization"], "runtime authorization")
    _require(
        type(runtime["attempt_id"]) is str
        and _ATTEMPT_ID.fullmatch(runtime["attempt_id"]) is not None,
        "runtime attempt id drifted",
    )
    source_commit = _source_commit(runtime["source_commit"], "runtime source commit")
    _require(
        source_commit == TRUSTED_Q011_RUNTIME_SOURCE_COMMIT,
        "runtime source commit is not the exact trusted Q011 source identity",
    )
    executable_sha256 = _sha256_text(
        runtime["executable_sha256"], "runtime executable digest"
    )
    normalization_sha256 = _sha256_text(
        runtime["runtime_normalization_sha256"], "runtime normalization digest"
    )
    _require(
        _source_commit(
            runtime["ps_escape_accounting_source_commit"],
            "runtime ps_escape accounting source commit",
        )
        == PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT,
        "runtime ps_escape accounting source commit drifted",
    )
    claim_escape_accounting = _exact_keys(
        runtime["claim_escape_accounting"],
        set(_CLAIM_ESCAPE_ACCOUNTING),
        "runtime claim escape accounting",
    )
    _require(
        dict(claim_escape_accounting) == _CLAIM_ESCAPE_ACCOUNTING,
        "runtime claim escape accounting modes drifted",
    )

    decoded_inventory, extrema, coverage = _validate_cycle_telemetry_inventory(
        evidence_root,
        runtime["per_cycle_inventory"],
        attempt_id=runtime["attempt_id"],
        source_commit=source_commit,
        executable_sha256=executable_sha256,
        normalization_sha256=normalization_sha256,
    )
    ps_escape = _validate_ps_escape_checkpoints(
        evidence_root, runtime["ps_escape_checkpoints"], inventory=decoded_inventory
    )

    exposure = _exact_keys(
        runtime["particle_exposure"], _PARTICLE_EXPOSURE_KEYS, "particle exposure"
    )
    exposure_complete = _strict_bool(exposure["complete"], "particle exposure complete")
    _require(
        exposure["method"]
        == "every_particle_update_and_pre_destruction_boundary_event",
        "particle exposure method drifted",
    )
    update_included = _strict_bool(
        exposure["active_particle_updates_included"], "active particle updates included"
    )
    boundary_included = _strict_bool(
        exposure["pre_destruction_boundary_events_included"],
        "pre-destruction boundary events included",
    )
    escaped_included = _strict_bool(
        exposure["escaped_particles_included"], "escaped particles included"
    )
    escaped_extrema_included = _strict_bool(
        exposure["escaped_particles_included_in_high_energy_and_gyroradius_extrema"],
        "escaped particles included in high-energy and gyroradius extrema",
    )
    observations = _nonnegative_int(exposure["observation_count"], "particle observations")
    _require(observations > 0, "particle exposure requires observations")
    exposure_statistics = {
        key: _finite_scalar(exposure[key], f"particle exposure {key}", minimum=0.0)
        for key in (
            "maximum_sampled_R",
            "maximum_sampled_Lambda",
            "cumulative_macro_weighted_R_exceedance_fraction",
            "cumulative_CR_energy_weighted_R_exceedance_fraction",
            "cumulative_macro_weighted_Lambda_exceedance_fraction",
            "cumulative_CR_energy_weighted_Lambda_exceedance_fraction",
        )
    }
    for key in (
        "cumulative_macro_weighted_R_exceedance_fraction",
        "cumulative_CR_energy_weighted_R_exceedance_fraction",
        "cumulative_macro_weighted_Lambda_exceedance_fraction",
        "cumulative_CR_energy_weighted_Lambda_exceedance_fraction",
    ):
        _require(exposure_statistics[key] <= 1.0, f"particle exposure {key} exceeds unity")
    _require(
        exposure_statistics["maximum_sampled_R"]
        <= extrema["R_maximum"] + 1.0e-12 * max(1.0, extrema["R_maximum"])
        and exposure_statistics["maximum_sampled_Lambda"]
        <= extrema["Lambda_maximum"] + 1.0e-12 * max(1.0, extrema["Lambda_maximum"]),
        "particle exposure maxima exceed all-cycle cell maxima",
    )
    if extrema["R_maximum"] <= R_MAXIMUM:
        _require(
            exposure_statistics["cumulative_macro_weighted_R_exceedance_fraction"] == 0.0
            and exposure_statistics[
                "cumulative_CR_energy_weighted_R_exceedance_fraction"
            ]
            == 0.0,
            "particle R exposure reports exceedance below the all-cycle cell maximum bound",
        )
    if extrema["Lambda_maximum"] <= LAMBDA_MAXIMUM:
        _require(
            exposure_statistics[
                "cumulative_macro_weighted_Lambda_exceedance_fraction"
            ]
            == 0.0
            and exposure_statistics[
                "cumulative_CR_energy_weighted_Lambda_exceedance_fraction"
            ]
            == 0.0,
            "particle Lambda exposure reports exceedance below the all-cycle cell maximum bound",
        )

    escape = _exact_keys(
        runtime["boundary_escape_ledger"], _ESCAPE_KEYS, "boundary escape ledger"
    )
    escape_complete = _strict_bool(escape["complete"], "escape ledger complete")
    _require(
        escape["scope"] == "all_Q011_shock_injected_particles_including_startup_removal",
        "escape ledger scope drifted",
    )
    faces = _exact_keys(
        escape["nonperiodic_faces"], {"inner_x1", "outer_x1"}, "nonperiodic escape faces"
    )
    decoded_faces = {
        "inner_x1": _decode_face_escape(
            faces["inner_x1"], "inner_x1 escaped state", expected_reason=INNER_X1_ESCAPE_REASON
        ),
        "outer_x1": _decode_face_escape(
            faces["outer_x1"], "outer_x1 escaped state", expected_reason=OUTER_X1_ESCAPE_REASON
        ),
    }
    _require(
        decoded_faces["inner_x1"]["particle_count"] == 0
        and decoded_faces["inner_x1"]["macro_weight"] == 0.0
        and decoded_faces["inner_x1"]["kinetic_energy"] == 0.0
        and decoded_faces["inner_x1"]["momentum"] == [0.0, 0.0, 0.0],
        "any inner_x1 escape rejects physical applicability",
    )
    _require(
        escape["periodic_faces"] == ["ix2", "ox2"],
        "escape ledger periodic face inventory drifted",
    )
    accumulated = _decode_state_vector(
        escape["accumulated_escaped"], "accumulated escaped state"
    )
    active = _decode_state_vector(escape["terminal_active"], "terminal active state")
    startup = _decode_state_vector(escape["startup_removed"], "startup removed state")
    injected_count = _nonnegative_int(
        escape["injected_particle_count"], "injected particle count"
    )
    injected_weight = _finite_scalar(
        escape["injected_macro_weight"], "injected macro weight", minimum=0.0
    )
    _require(injected_count > 0 and injected_weight > 0.0, "escape ledger source is empty")
    face_counts = [
        decoded_faces[face]["particle_count"] for face in ("inner_x1", "outer_x1")
    ]
    _require(
        accumulated["particle_count"] == sum(face_counts),
        "escape face particle counts do not close to accumulated escape",
    )
    residuals = {
        "macro_weight": _closure_residual(
            accumulated["macro_weight"],
            sum(decoded_faces[face]["macro_weight"] for face in ("inner_x1", "outer_x1")),
            [decoded_faces[face]["macro_weight"] for face in ("inner_x1", "outer_x1")],
        ),
        "kinetic_energy": _closure_residual(
            accumulated["kinetic_energy"],
            sum(decoded_faces[face]["kinetic_energy"] for face in ("inner_x1", "outer_x1")),
            [decoded_faces[face]["kinetic_energy"] for face in ("inner_x1", "outer_x1")],
        ),
        "momentum": [
            _closure_residual(
                accumulated["momentum"][component],
                sum(
                    decoded_faces[face]["momentum"][component]
                    for face in ("inner_x1", "outer_x1")
                ),
                [
                    decoded_faces[face]["momentum"][component]
                    for face in ("inner_x1", "outer_x1")
                ],
            )
            for component in range(3)
        ],
    }
    _require(
        residuals["macro_weight"] <= ESCAPE_LEDGER_RELATIVE_RESIDUAL_MAXIMUM
        and residuals["kinetic_energy"] <= ESCAPE_LEDGER_RELATIVE_RESIDUAL_MAXIMUM
        and max(residuals["momentum"]) <= ESCAPE_LEDGER_RELATIVE_RESIDUAL_MAXIMUM,
        "escape face weight, kinetic-energy, or momentum residual exceeds limit",
    )
    _require(
        injected_count
        == active["particle_count"]
        + startup["particle_count"]
        + accumulated["particle_count"],
        "escape ledger source particle census does not close",
    )
    source_weight_residual = _closure_residual(
        injected_weight,
        active["macro_weight"] + startup["macro_weight"] + accumulated["macro_weight"],
        [active["macro_weight"], startup["macro_weight"], accumulated["macro_weight"]],
    )
    _require(
        source_weight_residual <= ESCAPE_LEDGER_RELATIVE_RESIDUAL_MAXIMUM,
        "escape ledger source macro-weight census residual exceeds limit",
    )
    if accumulated["particle_count"] == 0:
        _require(
            accumulated["macro_weight"] == 0.0
            and accumulated["kinetic_energy"] == 0.0
            and extrema["escaped_particle_rg_maximum_over_Ly"] == 0.0
            and extrema["escaped_high_energy_tail_rg_maximum_over_Ly"] == 0.0
            and extrema["escaped_particle_specific_kinetic_energy_maximum"] == 0.0,
            "zero escaped-particle count requires zero escaped applicability inventory",
        )
    else:
        _require(
            accumulated["macro_weight"] > 0.0
            and accumulated["kinetic_energy"] > 0.0
            and extrema["escaped_particle_rg_maximum_over_Ly"] > 0.0
            and extrema["escaped_particle_specific_kinetic_energy_maximum"]
            >= accumulated["kinetic_energy"] / accumulated["macro_weight"],
            "positive escaped census lacks mass, energy, gyroradius, or maximum-energy evidence",
        )
    escape_exposure = _strict_bool(
        escape["escaped_particles_included_in_exposure"],
        "escape ledger particles included in exposure",
    )
    escape_extrema = _strict_bool(
        escape["escaped_particles_included_in_high_energy_and_gyroradius_extrema"],
        "escape ledger particles included in high-energy and gyroradius extrema",
    )
    terminal_ps = ps_escape["terminal"]
    terminal_particle_inventory = terminal_ps["active_particle_inventory"]
    energy_denominator = (
        accumulated["kinetic_energy"]
        + terminal_particle_inventory[
            "active_injected_cr_kinetic_energy_lower_bound"
        ]
    )
    escape_fractions = {
        "particle_count": accumulated["particle_count"] / injected_count,
        "macro_weight": accumulated["macro_weight"] / injected_weight,
        "kinetic_energy": (
            accumulated["kinetic_energy"] / energy_denominator
            if energy_denominator > 0.0
            else 0.0
        ),
    }
    terminal_ledger = terminal_ps["ps_escape_ledger"]
    _require(
        terminal_ledger["ps_escaped_injected_cr_count_global"]
        == accumulated["particle_count"]
        and terminal_ledger["ps_injected_cr_count_global"] == injected_count
        and terminal_ledger["ps_removed_cr_count_global"] == startup["particle_count"]
        and terminal_ps["active_injected_cr_count_global"] == active["particle_count"],
        "terminal schema-1 ps_escape counts disagree with runtime escape ledger",
    )
    sealed_residuals = {
        "escaped_mass": _closure_residual(
            terminal_ledger["ps_escaped_injected_cr_mass_global"],
            accumulated["macro_weight"],
            [accumulated["macro_weight"]],
        ),
        "injected_mass": _closure_residual(
            terminal_ledger["ps_injected_cr_mass_global"],
            injected_weight,
            [injected_weight],
        ),
        "startup_removed_mass": _closure_residual(
            terminal_ledger["ps_removed_cr_mass_global"],
            startup["macro_weight"],
            [startup["macro_weight"]],
        ),
        "active_mass": _closure_residual(
            terminal_ps["active_injected_cr_mass_global"],
            active["macro_weight"],
            [active["macro_weight"]],
        ),
        "active_energy": _closure_residual(
            terminal_ps["active_injected_cr_kinetic_energy_global"],
            active["kinetic_energy"],
            [active["kinetic_energy"]],
        ),
        "escaped_energy": _closure_residual(
            terminal_ledger["ps_escaped_injected_cr_energy_global"],
            accumulated["kinetic_energy"],
            [accumulated["kinetic_energy"]],
        ),
        "escaped_momentum": [
            _closure_residual(
                terminal_ledger[f"ps_escaped_injected_cr_momentum_x{component + 1}_global"],
                accumulated["momentum"][component],
                [accumulated["momentum"][component]],
            )
            for component in range(3)
        ],
    }
    _require(
        max(
            [
                sealed_residuals["escaped_mass"],
                sealed_residuals["injected_mass"],
                sealed_residuals["startup_removed_mass"],
                sealed_residuals["active_mass"],
                sealed_residuals["active_energy"],
                sealed_residuals["escaped_energy"],
                *sealed_residuals["escaped_momentum"],
            ]
        )
        <= ESCAPE_LEDGER_RELATIVE_RESIDUAL_MAXIMUM,
        "terminal schema-1 ps_escape quantities disagree with runtime escape ledger",
    )
    _require(
        terminal_ps["escaped_injected_max_specific_kinetic_energy_global"]
        == extrema["escaped_particle_specific_kinetic_energy_maximum"]
        and terminal_ps["escaped_injected_max_rg_over_Ly_global"]
        == extrema["escaped_particle_rg_maximum_over_Ly"],
        "terminal schema-1 ps_escape maxima disagree with all-cycle escaped extrema",
    )
    slope_cutoff_escape = _validate_slope_cutoff_escape_evidence(
        evidence_root,
        runtime["escaped_slope_cutoff_evidence"],
        attempt_id=runtime["attempt_id"],
        source_commit=source_commit,
        executable_sha256=executable_sha256,
        normalization_sha256=normalization_sha256,
        active=active,
        escaped=accumulated,
        maximum_specific_energy=extrema["maximum_particle_specific_kinetic_energy"],
    )
    escaped_applicability_complete = (
        exposure_complete
        and update_included
        and boundary_included
        and escape_complete
        and escaped_included
        and escaped_extrema_included
        and escape_exposure
        and escape_extrema
    )
    claim_specific_escape_applicability = {
        claim: {
            "pass": (
                escaped_applicability_complete
                and escape_fractions["particle_count"]
                <= float(bounds["particle_count_fraction_maximum"])
                and escape_fractions["macro_weight"]
                <= float(bounds["macro_weight_fraction_maximum"])
                and escape_fractions["kinetic_energy"]
                <= float(bounds["kinetic_energy_fraction_maximum"])
                and (
                    not bool(bounds["requires_escaped_maximum_inclusion"])
                    or escaped_extrema_included
                )
                and (
                    claim != "high_energy_slope_or_cutoff_claim"
                    or slope_cutoff_escape["pass"]
                )
            ),
            "accounting_mode": claim_escape_accounting[claim],
            "observed_escape_fractions": dict(escape_fractions),
            "AthenaK_selected_preregistered_bounds": dict(bounds),
            "slope_cutoff_binwise_tail_evidence_required": (
                claim == "high_energy_slope_or_cutoff_claim"
            ),
        }
        for claim, bounds in CLAIM_SPECIFIC_ESCAPE_BOUNDS.items()
    }
    time_pass = (
        exposure_complete
        and update_included
        and boundary_included
        and escape_complete
        and escaped_applicability_complete
    )
    return {
        "binding": binding,
        "authorization": authorization,
        "attempt_id": runtime["attempt_id"],
        "source_commit": source_commit,
        "executable_sha256": executable_sha256,
        "runtime_normalization_sha256": normalization_sha256,
        "cycle_coverage": coverage,
        "per_cycle_inventory": decoded_inventory,
        "all_cycle_extrema": extrema,
        "ps_escape_accounting": ps_escape,
        "claim_escape_accounting": dict(claim_escape_accounting),
        "claim_specific_escape_applicability": claim_specific_escape_applicability,
        "slope_cutoff_escape_evidence": slope_cutoff_escape,
        "particle_exposure": {**dict(exposure), **exposure_statistics},
        "boundary_escape_ledger": {
            **dict(escape),
            "nonperiodic_faces": decoded_faces,
            "accumulated_escaped": accumulated,
            "terminal_active": active,
            "startup_removed": startup,
            "face_to_accumulated_relative_residuals": residuals,
            "source_macro_weight_relative_residual": source_weight_residual,
            "sealed_ps_escape_relative_residuals": sealed_residuals,
            "escape_fractions": escape_fractions,
        },
        "escaped_particle_applicability_complete": escaped_applicability_complete,
        "time_and_escape_complete": time_pass,
    }




def _validate_recorded_artifact_binding(
    value: object, role: str, label: str
) -> Mapping[str, Any]:
    binding = _exact_keys(value, {"role", "path", "sha256", "byte_count"}, label)
    _require(binding["role"] == role, f"{label} role drifted")
    _safe_relative_path(binding["path"], f"{label} path")
    _sha256_text(binding["sha256"], f"{label} sha256")
    _nonnegative_int(binding["byte_count"], f"{label} byte count")
    return binding


def _revalidate_retained_snapshot_artifacts(
    snapshot: ApplicabilitySnapshot, record: Mapping[str, Any]
) -> None:
    _require(
        isinstance(snapshot.evidence_root, Path),
        "snapshot retained evidence root must be a pathlib.Path",
    )
    root = snapshot.evidence_root.resolve(strict=True)
    _require(root.is_dir(), "snapshot retained evidence root must remain a directory")
    provenance = record["snapshot_provenance"]
    raw_products = provenance["raw_products"]
    raw_payloads: dict[str, bytes] = {}
    for product in REQUIRED_RAW_PRODUCTS:
        binding, raw_payloads[product] = _read_bound_artifact(
            root, raw_products[product], expected_role=product
        )
        _require(
            binding == dict(raw_products[product]),
            f"retained snapshot {product} binding drifted",
        )
    parsed_datasets: dict[str, Any] = {}
    for product in ("mhd_w_bcc", *science.CURRENT_PRODUCT_FIELDS):
        try:
            parsed_datasets[product] = output_primitives.parse_athenak_binary_bytes(
                raw_payloads[product], source=str(raw_products[product]["path"])
            )
        except output_primitives.AnalysisError as error:
            raise PhysicalApplicabilityError(
                f"retained snapshot {product} is not a trusted Athena binary: {error}"
            ) from error
        _require(
            _dataset_sha256(parsed_datasets[product], f"retained snapshot {product}")
            == provenance["decoded_product_sha256"][product],
            f"retained snapshot {product} decoded digest drifted",
        )
    particle = _decode_particle_vtk_bytes(
        raw_payloads["prtcl_all"], "retained snapshot prtcl_all"
    )
    _require(
        particle["execution_header"]["cycle"] == provenance["cycle"]
        and particle["execution_header"]["observed_committed_time"]
        == provenance["observed_committed_time"],
        "retained snapshot prtcl_all cycle/time binding drifted",
    )
    _require(
        _particle_payload_sha256(
            points=particle["points"],
            cr_source=particle["cr_source"],
            birth_time=particle["birth_time"],
            velocity=particle["velocity"],
            macro_weight=particle["macro_weight"],
        )
        == provenance["decoded_product_sha256"]["prtcl_all"],
        "retained snapshot prtcl_all decoded digest drifted",
    )
    normalization = record["bound_normalization_evidence"]
    revalidated_normalization = _validate_bound_normalization_evidence(
        root,
        normalization["bindings"],
        runtime_input_parameters=parsed_datasets["mhd_w_bcc"].input_parameters,
    )
    _require(
        revalidated_normalization["source_commit"] == normalization["source_commit"]
        and revalidated_normalization["runtime_input_parameters_sha256"]
        == normalization["runtime_input_parameters_sha256"]
        and revalidated_normalization["source_manifest"]
        == normalization["source_manifest"]
        and revalidated_normalization["bindings"] == normalization["bindings"],
        "retained snapshot normalization/source evidence drifted",
    )
    manifest_binding, manifest_payload = _read_bound_artifact(
        root, provenance["manifest_binding"], expected_role="snapshot_manifest"
    )
    manifest = _exact_keys(
        _decode_canonical_json(
            manifest_payload, "retained snapshot provenance manifest"
        ),
        {
            "schema_version",
            "record_type",
            "attempt_id",
            "source_commit",
            "executable_sha256",
            "runtime_normalization_sha256",
            "nominal_slot_time",
            "observed_committed_time",
            "cycle",
            "raw_products",
            "decoded_product_sha256",
        },
        "retained snapshot provenance manifest",
    )
    _require(
        manifest_binding == dict(provenance["manifest_binding"])
        and manifest["attempt_id"] == provenance["attempt_id"]
        and manifest["source_commit"] == provenance["source_commit"]
        and manifest["executable_sha256"] == provenance["executable_sha256"]
        and manifest["runtime_normalization_sha256"]
        == provenance["runtime_normalization_sha256"]
        and manifest["nominal_slot_time"] == provenance["nominal_slot_time"]
        and manifest["observed_committed_time"]
        == provenance["observed_committed_time"]
        and manifest["cycle"] == provenance["cycle"]
        and manifest["raw_products"] == provenance["raw_products"]
        and manifest["decoded_product_sha256"] == provenance["decoded_product_sha256"],
        "retained snapshot provenance manifest drifted",
    )


def _validate_snapshot(snapshot: object) -> Mapping[str, Any]:
    _require(
        isinstance(snapshot, ApplicabilitySnapshot),
        "history requires complete ApplicabilitySnapshot objects with retained maps",
    )
    record = _exact_keys(
        snapshot.record,
        {
            "schema_version",
            "record_type",
            "successor_id",
            "qualification_effect",
            "authorization",
            "nominal_slot_time",
            "observed_committed_time",
            "detected_front_x1_c_over_omega_pi",
            "exact_normalization",
            "bound_normalization_evidence",
            "snapshot_provenance",
            "formulae",
            "cell_map_contract",
            "regional_statistics",
            "ion_scale_separation",
            "particle_R_Lambda_exposure",
            "particle_gyroradius_containment",
            "gates",
            "snapshot_gate_pass_excluding_time_completeness",
            "claim_rejection_rules",
            "permanent_claim_exclusions",
            "claim_rejections",
            "runtime_time_escape_evidence_required",
        },
        "snapshot applicability record",
    )
    _require(
        type(record["schema_version"]) is int and record["schema_version"] == SCHEMA_VERSION,
        "snapshot schema version drifted",
    )
    _require(record["record_type"] == SNAPSHOT_RECORD_TYPE, "snapshot record type drifted")
    _require(record["successor_id"] == SUCCESSOR_ID, "snapshot successor id drifted")
    _require(
        record["qualification_effect"] == QUALIFICATION_EFFECT,
        "snapshot qualification effect drifted",
    )
    _validate_authorization(record["authorization"], "snapshot authorization")
    _validate_exact_normalization(record["exact_normalization"])
    nominal = _finite_scalar(record["nominal_slot_time"], "snapshot nominal time", minimum=0.0)
    observed = _finite_scalar(
        record["observed_committed_time"], "snapshot observed time", minimum=0.0
    )
    front_x = _finite_scalar(
        record["detected_front_x1_c_over_omega_pi"], "snapshot detected front"
    )
    _require(
        _strict_bool(
            record["runtime_time_escape_evidence_required"],
            "snapshot runtime time/escape requirement",
        ),
        "snapshot must require runtime time/escape evidence",
    )

    normalization = _exact_keys(
        record["bound_normalization_evidence"],
        {"source_commit", "runtime_input_parameters_sha256", "source_manifest", "bindings"},
        "snapshot bound normalization evidence",
    )
    _source_commit(normalization["source_commit"], "snapshot normalization source commit")
    _sha256_text(
        normalization["runtime_input_parameters_sha256"],
        "snapshot runtime input-parameter digest",
    )
    normalization_bindings = _exact_keys(
        normalization["bindings"],
        {
            "runtime_normalization_record",
            "deck",
            "source_manifest",
            "source_archive",
            "executable",
        },
        "snapshot normalization bindings",
    )
    for role, binding in normalization_bindings.items():
        _validate_recorded_artifact_binding(
            binding, role, f"snapshot normalization {role} binding"
        )
    source_manifest = _exact_keys(
        normalization["source_manifest"],
        {
            "schema_version",
            "record_type",
            "source_commit",
            "deck_sha256",
            "source_archive_sha256",
            "executable_sha256",
        },
        "snapshot source manifest",
    )
    _require(
        source_manifest["schema_version"] == SCHEMA_VERSION
        and source_manifest["record_type"] == SOURCE_MANIFEST_RECORD_TYPE
        and source_manifest["source_commit"] == normalization["source_commit"],
        "snapshot source manifest identity drifted",
    )
    for role in ("deck", "source_archive", "executable"):
        _require(
            source_manifest[f"{role}_sha256"] == normalization_bindings[role]["sha256"],
            f"snapshot source manifest {role} binding drifted",
        )

    provenance = _exact_keys(
        record["snapshot_provenance"],
        {
            "manifest_binding",
            "attempt_id",
            "source_commit",
            "executable_sha256",
            "runtime_normalization_sha256",
            "nominal_slot_time",
            "observed_committed_time",
            "cycle",
            "raw_products",
            "decoded_product_sha256",
        },
        "snapshot provenance",
    )
    _validate_recorded_artifact_binding(
        provenance["manifest_binding"], "snapshot_manifest", "snapshot manifest binding"
    )
    _require(
        type(provenance["attempt_id"]) is str
        and _ATTEMPT_ID.fullmatch(provenance["attempt_id"]) is not None,
        "snapshot attempt id drifted",
    )
    _require(
        _source_commit(provenance["source_commit"], "snapshot provenance source commit")
        == normalization["source_commit"],
        "snapshot provenance source commit disagrees with normalization evidence",
    )
    _require(
        _sha256_text(provenance["executable_sha256"], "snapshot executable digest")
        == normalization_bindings["executable"]["sha256"],
        "snapshot executable digest disagrees with normalization evidence",
    )
    _require(
        _sha256_text(
            provenance["runtime_normalization_sha256"],
            "snapshot runtime normalization digest",
        )
        == normalization_bindings["runtime_normalization_record"]["sha256"],
        "snapshot runtime normalization digest disagrees with normalization evidence",
    )
    _require(
        _finite_scalar(provenance["nominal_slot_time"], "provenance nominal time")
        == nominal
        and _finite_scalar(
            provenance["observed_committed_time"], "provenance observed time"
        )
        == observed,
        "snapshot provenance times drifted",
    )
    _nonnegative_int(provenance["cycle"], "snapshot provenance cycle")
    raw_products = _exact_keys(
        provenance["raw_products"], set(REQUIRED_RAW_PRODUCTS), "snapshot raw products"
    )
    for role, binding in raw_products.items():
        _validate_recorded_artifact_binding(binding, role, f"snapshot raw {role} binding")
    decoded_products = _exact_keys(
        provenance["decoded_product_sha256"],
        set(REQUIRED_RAW_PRODUCTS),
        "snapshot decoded-product digests",
    )
    for role, digest in decoded_products.items():
        _sha256_text(digest, f"snapshot decoded {role} digest")

    maps = snapshot.cell_maps
    _require(isinstance(maps, Mapping), "snapshot cell maps must be a mapping")
    _require(set(maps) == set(CELL_MAP_NAMES), "snapshot cell map inventory drifted")
    contract = _exact_keys(
        record["cell_map_contract"],
        {
            "returned_separately_as_immutable_numpy_arrays",
            "retention_required_for_physical_applicability_evidence",
            "shape_y_x",
            "target_composite_level",
            "x1_faces_c_over_omega_pi",
            "x2_faces_c_over_omega_pi",
            "maps",
            "map_bindings",
            "per_cell_maxima_required_no_rare_cell_waiver",
        },
        "snapshot cell-map contract",
    )
    for key in (
        "returned_separately_as_immutable_numpy_arrays",
        "retention_required_for_physical_applicability_evidence",
        "per_cell_maxima_required_no_rare_cell_waiver",
    ):
        _require(_strict_bool(contract[key], f"snapshot cell-map {key}"), f"{key} drifted")
    _nonnegative_int(contract["target_composite_level"], "snapshot target composite level")
    _require(contract["maps"] == list(CELL_MAP_NAMES), "snapshot cell-map list drifted")
    shape = contract["shape_y_x"]
    _require(
        type(shape) is list
        and len(shape) == 2
        and all(type(value) is int and value > 0 for value in shape),
        "snapshot cell-map shape drifted",
    )
    x1_faces = _finite_array(
        contract["x1_faces_c_over_omega_pi"], "snapshot x1 faces", ndim=1
    )
    x2_faces = _finite_array(
        contract["x2_faces_c_over_omega_pi"], "snapshot x2 faces", ndim=1
    )
    _require(
        x1_faces.size == shape[1] + 1
        and x2_faces.size == shape[0] + 1
        and np.all(np.diff(x1_faces) > 0.0)
        and np.all(np.diff(x2_faces) > 0.0),
        "snapshot cell-map faces disagree with shape",
    )
    map_bindings = _exact_keys(
        contract["map_bindings"], set(CELL_MAP_NAMES), "snapshot cell-map bindings"
    )
    decoded_maps: dict[str, np.ndarray] = {}
    for name in CELL_MAP_NAMES:
        values = _finite_array(maps[name], f"snapshot {name} map", ndim=2)
        _require(values.shape == tuple(shape), f"snapshot {name} map shape drifted")
        _require(not values.flags.writeable, f"snapshot {name} map must be immutable")
        _require(np.all(values >= 0.0), f"snapshot {name} map must be non-negative")
        binding = _exact_keys(
            map_bindings[name],
            {"shape_y_x", "dtype", "sha256"},
            f"snapshot {name} map binding",
        )
        _require(
            binding["shape_y_x"] == shape and binding["dtype"] == "float64",
            f"snapshot {name} map binding metadata drifted",
        )
        _require(
            _sha256_text(binding["sha256"], f"snapshot {name} map digest")
            == _array_sha256(values),
            f"snapshot {name} map digest drifted",
        )
        decoded_maps[name] = values

    regions = _exact_keys(
        record["regional_statistics"], set(DETECTED_FRONT_REGIONS), "snapshot regions"
    )
    x1_centers = 0.5 * (x1_faces[:-1] + x1_faces[1:])
    areas = np.diff(x2_faces)[:, None] * np.diff(x1_faces)[None, :]
    for name, offsets in DETECTED_FRONT_REGIONS.items():
        region = _exact_keys(
            regions[name],
            {"detected_front_offsets_c_over_omega_pi", "R", "Lambda"},
            f"snapshot region {name}",
        )
        expected_offsets = None if offsets is None else list(offsets)
        _require(
            region["detected_front_offsets_c_over_omega_pi"] == expected_offsets,
            f"snapshot region {name} offsets drifted",
        )
        selected_x = (
            np.ones(x1_centers.shape, dtype=bool)
            if offsets is None
            else (x1_centers > front_x + offsets[0])
            & (x1_centers < front_x + offsets[1])
        )
        mask = np.broadcast_to(selected_x[None, :], tuple(shape))
        _require(np.any(mask), f"snapshot region {name} is empty")
        expected_current = areas[mask] * decoded_maps["gas_frame_current_magnitude"][mask]
        for field in ("R", "Lambda"):
            recomputed = _weighted_statistics(
                decoded_maps[field][mask],
                areas[mask],
                expected_current,
                label=f"{name} {field}",
            )
            _require(
                _canonical_json_bytes(region[field]) == _canonical_json_bytes(recomputed),
                f"snapshot region {name} {field} statistics disagree with retained maps",
            )

    ion = _exact_keys(
        record["ion_scale_separation"],
        {
            "actual_leaf_spacing_used",
            "shock_transition_exclusion",
            "S_delta_minimum_excluding_shock_transition",
            "precursor_magnetic_spectrum",
        },
        "snapshot ion-scale separation",
    )
    _require(
        _strict_bool(ion["actual_leaf_spacing_used"], "snapshot actual leaf spacing used"),
        "snapshot must use actual leaf spacing",
    )
    spectrum = _exact_keys(
        ion["precursor_magnetic_spectrum"],
        {
            "region",
            "method",
            "actual_leaf_aware",
            "finest_composite_FFT_rejected",
            "source_levels_present",
            "analysis_source_level",
            "target_composite_level",
            "restriction_factor",
            "analysis_dx1",
            "analysis_dx2",
            "selected_analysis_nx1",
            "selected_analysis_nx2",
            "delta_B_rms_over_B0",
            "total_windowed_delta_B_power",
            "resolved_restricted_power",
            "subgrid_residual_power_upper_bound",
            "subgrid_residual_power_fraction_upper_bound",
            "local_di_maximum",
            "lambda_B_characteristic",
            "lambda_B_characteristic_over_local_di_maximum",
            "sub_10di_magnetic_power_fraction_upper_bound",
            "shock_transition_excluded",
        },
        "snapshot precursor magnetic spectrum",
    )
    _require(
        spectrum["region"] == DI_MAGNETIC_SPECTRUM_REGION
        and _strict_bool(spectrum["actual_leaf_aware"], "snapshot leaf-aware spectrum")
        and _strict_bool(
            spectrum["finest_composite_FFT_rejected"],
            "snapshot finest-composite FFT rejection",
        )
        and _strict_bool(
            spectrum["shock_transition_excluded"], "snapshot spectrum shock exclusion"
        ),
        "snapshot magnetic spectrum method drifted",
    )
    _require(
        type(spectrum["source_levels_present"]) is list
        and len(spectrum["source_levels_present"]) > 0
        and all(type(level) is int and level >= 0 for level in spectrum["source_levels_present"]),
        "snapshot magnetic spectrum source-level inventory drifted",
    )
    for key in set(spectrum) - {
        "region",
        "method",
        "actual_leaf_aware",
        "finest_composite_FFT_rejected",
        "source_levels_present",
        "shock_transition_excluded",
    }:
        if key in {
            "analysis_source_level",
            "target_composite_level",
            "restriction_factor",
            "selected_analysis_nx1",
            "selected_analysis_nx2",
        }:
            _nonnegative_int(spectrum[key], f"snapshot spectrum {key}")
        else:
            _finite_scalar(spectrum[key], f"snapshot spectrum {key}", minimum=0.0)

    exposure = _exact_keys(
        record["particle_R_Lambda_exposure"],
        {
            "sampling",
            "population_definitions",
            "high_energy_tail_specific_kinetic_energy_threshold",
            "populations",
        },
        "snapshot particle R/Lambda exposure",
    )
    _finite_scalar(
        exposure["high_energy_tail_specific_kinetic_energy_threshold"],
        "snapshot high-energy-tail threshold",
        minimum=0.0,
    )
    _exact_keys(
        exposure["population_definitions"],
        {"all_active", "detected_front_upstream", "high_energy_tail"},
        "snapshot particle population definitions",
    )
    populations = _exact_keys(
        exposure["populations"],
        {"all_active", "detected_front_upstream", "high_energy_tail"},
        "snapshot particle exposure populations",
    )
    for name, population in populations.items():
        population = _exact_keys(
            population,
            {"available", "reason", "particle_count", "R", "Lambda"},
            f"snapshot particle exposure population {name}",
        )
        _require(
            _strict_bool(population["available"], f"snapshot {name} availability"),
            f"snapshot requires complete available particle population {name}",
        )
        _require(population["reason"] is None, f"snapshot {name} reason must be null")
        _require(
            _nonnegative_int(population["particle_count"], f"snapshot {name} count") > 0,
            f"snapshot {name} particle population is empty",
        )
        for field in ("R", "Lambda"):
            stats = _exact_keys(
                population[field],
                {
                    "particle_count",
                    "local_minimum",
                    "local_maximum",
                    "threshold",
                    "macro_weighted",
                    "CR_kinetic_energy_weighted",
                },
                f"snapshot {name} {field} particle statistics",
            )
            _require(
                _nonnegative_int(stats["particle_count"], f"snapshot {name} {field} count")
                == population["particle_count"],
                f"snapshot {name} {field} particle count drifted",
            )
            for scalar in ("local_minimum", "local_maximum", "threshold"):
                _finite_scalar(stats[scalar], f"snapshot {name} {field} {scalar}", minimum=0.0)
            for weighting in ("macro_weighted", "CR_kinetic_energy_weighted"):
                weighted = _exact_keys(
                    stats[weighting],
                    {"weighted_mean", "weighted_quantiles", "exceedance_fraction"},
                    f"snapshot {name} {field} {weighting}",
                )
                _finite_scalar(
                    weighted["weighted_mean"],
                    f"snapshot {name} {field} {weighting} mean",
                    minimum=0.0,
                )
                quantiles = _exact_keys(
                    weighted["weighted_quantiles"],
                    {"q500", "q900", "q990", "q999"},
                    f"snapshot {name} {field} {weighting} quantiles",
                )
                for quantile, value in quantiles.items():
                    _finite_scalar(
                        value,
                        f"snapshot {name} {field} {weighting} {quantile}",
                        minimum=0.0,
                    )
                fraction = _finite_scalar(
                    weighted["exceedance_fraction"],
                    f"snapshot {name} {field} {weighting} exceedance",
                    minimum=0.0,
                )
                _require(fraction <= 1.0, "snapshot particle exceedance exceeds unity")

    particle = _exact_keys(
        record["particle_gyroradius_containment"],
        {
            "population",
            "local_B_sampling",
            "gyroradius_formula",
            "selected_particle_count",
            "transverse_domain_size_Ly",
            "macro_weighted",
            "energy_weighted",
            "maximum",
            "maximum_specific_kinetic_energy",
            "high_energy_tail_definition",
            "high_energy_tail_specific_kinetic_energy_threshold",
            "high_energy_tail_particle_count",
            "high_energy_tail_maximum",
            "macro_q999_over_Ly",
            "energy_q999_over_Ly",
            "maximum_over_Ly",
            "high_energy_tail_maximum_over_Ly",
            "energy_fraction_with_rg_above_Ly_over_4",
        },
        "snapshot particle gyroradius containment",
    )
    ly = _finite_scalar(
        particle["transverse_domain_size_Ly"], "snapshot transverse domain size", minimum=0.0
    )
    _require(ly > 0.0, "snapshot transverse domain size must be positive")
    for key in (
        "maximum",
        "maximum_specific_kinetic_energy",
        "high_energy_tail_specific_kinetic_energy_threshold",
        "high_energy_tail_maximum",
        "macro_q999_over_Ly",
        "energy_q999_over_Ly",
        "maximum_over_Ly",
        "high_energy_tail_maximum_over_Ly",
        "energy_fraction_with_rg_above_Ly_over_4",
    ):
        _finite_scalar(particle[key], f"snapshot particle {key}", minimum=0.0)
    _require(
        particle["maximum_over_Ly"] == particle["maximum"] / ly
        and particle["high_energy_tail_maximum_over_Ly"]
        == particle["high_energy_tail_maximum"] / ly,
        "snapshot particle gyroradius ratios drifted",
    )
    _require(
        _nonnegative_int(
            particle["selected_particle_count"], "snapshot selected particle count"
        )
        >= PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES
        and _nonnegative_int(
            particle["high_energy_tail_particle_count"],
            "snapshot high-energy-tail particle count",
        )
        > 0,
        "snapshot particle containment population is incomplete",
    )
    for weighting in ("macro_weighted", "energy_weighted"):
        quantiles = _exact_keys(
            particle[weighting],
            {"q500", "q900", "q990", "q999"},
            f"snapshot gyroradius {weighting} quantiles",
        )
        for quantile, value in quantiles.items():
            _finite_scalar(
                value, f"snapshot gyroradius {weighting} {quantile}", minimum=0.0
            )

    gates = _exact_keys(
        record["gates"],
        {
            "Q011-APP-NORM",
            "Q011-APP-R",
            "Q011-APP-LAMBDA",
            "Q011-APP-DI",
            "Q011-APP-RG",
            "Q011-APP-TIME",
        },
        "snapshot gates",
    )
    shock = np.abs(x1_centers - front_x) <= SHOCK_TRANSITION_HALF_WIDTH
    di_mask = np.broadcast_to((~shock)[None, :], tuple(shape))
    expected_observables = {
        "R_maximum": float(np.max(decoded_maps["R"])),
        "Lambda_maximum": float(np.max(decoded_maps["Lambda"])),
        "S_delta_minimum_excluding_shock_transition": float(
            np.min(decoded_maps["S_delta"][di_mask])
        ),
        "lambda_B_characteristic_over_local_di_maximum": float(
            spectrum["lambda_B_characteristic_over_local_di_maximum"]
        ),
        "sub_10di_magnetic_power_fraction_upper_bound": float(
            spectrum["sub_10di_magnetic_power_fraction_upper_bound"]
        ),
        "delta_B_rms_over_B0": float(spectrum["delta_B_rms_over_B0"]),
        "particle_rg_maximum_over_Ly": float(particle["maximum_over_Ly"]),
        "high_energy_tail_rg_maximum_over_Ly": float(
            particle["high_energy_tail_maximum_over_Ly"]
        ),
        "macro_q999_over_Ly": float(particle["macro_q999_over_Ly"]),
        "energy_q999_over_Ly": float(particle["energy_q999_over_Ly"]),
        "energy_fraction_with_rg_above_Ly_over_4": float(
            particle["energy_fraction_with_rg_above_Ly_over_4"]
        ),
    }
    _require(
        gates["Q011-APP-NORM"]["pass"] is True
        and gates["Q011-APP-TIME"]["pass"] is False,
        "snapshot normalization/time gate drifted",
    )
    _require(
        gates["Q011-APP-R"]["observed_maximum"] == expected_observables["R_maximum"]
        and gates["Q011-APP-LAMBDA"]["observed_maximum"]
        == expected_observables["Lambda_maximum"]
        and gates["Q011-APP-DI"]["observed_S_delta_minimum_excluding_shock_transition"]
        == expected_observables["S_delta_minimum_excluding_shock_transition"]
        and gates["Q011-APP-DI"][
            "observed_lambda_B_characteristic_over_local_di_maximum"
        ]
        == expected_observables["lambda_B_characteristic_over_local_di_maximum"]
        and gates["Q011-APP-DI"][
            "observed_sub_10di_magnetic_power_fraction_upper_bound"
        ]
        == expected_observables["sub_10di_magnetic_power_fraction_upper_bound"]
        and gates["Q011-APP-DI"]["observed_delta_B_rms_over_B0"]
        == expected_observables["delta_B_rms_over_B0"]
        and gates["Q011-APP-RG"]["observed_maximum_over_Ly"]
        == expected_observables["particle_rg_maximum_over_Ly"]
        and gates["Q011-APP-RG"]["observed_high_energy_tail_maximum_over_Ly"]
        == expected_observables["high_energy_tail_rg_maximum_over_Ly"]
        and gates["Q011-APP-RG"]["observed_macro_q999_over_Ly"]
        == expected_observables["macro_q999_over_Ly"]
        and gates["Q011-APP-RG"]["observed_energy_q999_over_Ly"]
        == expected_observables["energy_q999_over_Ly"]
        and gates["Q011-APP-RG"]["observed_energy_fraction_with_rg_above_Ly_over_4"]
        == expected_observables["energy_fraction_with_rg_above_Ly_over_4"],
        "snapshot gate observables disagree with retained full evidence",
    )
    expected_pass = {
        "Q011-APP-NORM": True,
        "Q011-APP-R": expected_observables["R_maximum"] <= R_MAXIMUM,
        "Q011-APP-LAMBDA": expected_observables["Lambda_maximum"] <= LAMBDA_MAXIMUM,
        "Q011-APP-DI": (
            expected_observables["S_delta_minimum_excluding_shock_transition"]
            >= S_DELTA_MINIMUM
            and expected_observables["lambda_B_characteristic_over_local_di_maximum"]
            >= LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM
            and expected_observables["sub_10di_magnetic_power_fraction_upper_bound"]
            <= SUB_10DI_POWER_FRACTION_MAXIMUM
            and expected_observables["delta_B_rms_over_B0"]
            >= DELTA_B_RMS_OVER_B0_MINIMUM
        ),
        "Q011-APP-RG": (
            expected_observables["macro_q999_over_Ly"] <= RG_Q999_OVER_LY_MAXIMUM
            and expected_observables["energy_q999_over_Ly"] <= RG_Q999_OVER_LY_MAXIMUM
            and expected_observables["energy_fraction_with_rg_above_Ly_over_4"]
            <= RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM
            and expected_observables["particle_rg_maximum_over_Ly"]
            <= RG_MAXIMUM_OVER_LY_MAXIMUM
            and expected_observables["high_energy_tail_rg_maximum_over_Ly"]
            <= RG_HIGH_ENERGY_MAXIMUM_OVER_LY_MAXIMUM
        ),
        "Q011-APP-TIME": False,
    }
    for name, expected in expected_pass.items():
        _require(
            type(gates[name]["pass"]) is bool and gates[name]["pass"] is expected,
            f"snapshot {name} pass flag disagrees with bound observables",
        )
    _require(
        _strict_bool(
            record["snapshot_gate_pass_excluding_time_completeness"],
            "snapshot aggregate gate",
        )
        is all(expected_pass[name] for name in expected_pass if name != "Q011-APP-TIME"),
        "snapshot aggregate gate drifted",
    )
    _require(
        record["claim_rejection_rules"]
        == {gate: list(claims) for gate, claims in CLAIM_REJECTION_RULES.items()}
        and record["permanent_claim_exclusions"] == list(PERMANENT_CLAIM_EXCLUSIONS)
        and record["claim_rejections"]
        == _claim_rejections(gates, expected_observables["Lambda_maximum"]),
        "snapshot claim-rejection contract drifted",
    )
    _revalidate_retained_snapshot_artifacts(snapshot, record)
    return record


@_public_contract("Q011 physical-applicability history reduction")
def reduce_physical_applicability_history(
    snapshots: Sequence[ApplicabilitySnapshot],
    runtime_time_escape_evidence: object,
    *,
    evidence_root: Path,
) -> dict[str, Any]:
    """Combine complete bound snapshots with a bound actual per-cycle runtime record."""
    _require(
        isinstance(snapshots, Sequence) and not isinstance(snapshots, (str, bytes)),
        "snapshot applicability history must be a sequence",
    )
    _require(len(snapshots) > 0, "snapshot applicability history must not be empty")
    records = [_validate_snapshot(snapshot) for snapshot in snapshots]
    times = [
        _finite_scalar(record["observed_committed_time"], "snapshot observed time", minimum=0.0)
        for record in records
    ]
    _require(
        all(right > left for left, right in zip(times, times[1:])),
        "snapshot applicability history times must strictly increase",
    )
    runtime = _validate_runtime_time_escape(evidence_root, runtime_time_escape_evidence)
    extrema = runtime["all_cycle_extrema"]
    _require(
        times[0] >= runtime["cycle_coverage"]["post_startup_removal_start_time"]
        and times[-1] <= runtime["cycle_coverage"]["terminal_time"],
        "snapshot applicability history escaped runtime coverage",
    )
    for record in records:
        provenance = record["snapshot_provenance"]
        _require(
            provenance["attempt_id"] == runtime["attempt_id"]
            and provenance["source_commit"] == runtime["source_commit"]
            and provenance["executable_sha256"] == runtime["executable_sha256"]
            and provenance["runtime_normalization_sha256"]
            == runtime["runtime_normalization_sha256"],
            "snapshot provenance disagrees with bound runtime evidence",
        )
    snapshot_extrema = {
        "R_maximum": max(record["gates"]["Q011-APP-R"]["observed_maximum"] for record in records),
        "Lambda_maximum": max(
            record["gates"]["Q011-APP-LAMBDA"]["observed_maximum"] for record in records
        ),
        "S_delta_minimum_excluding_shock_transition": min(
            record["gates"]["Q011-APP-DI"][
                "observed_S_delta_minimum_excluding_shock_transition"
            ]
            for record in records
        ),
        "lambda_B_characteristic_over_local_di_maximum_minimum": min(
            record["gates"]["Q011-APP-DI"][
                "observed_lambda_B_characteristic_over_local_di_maximum"
            ]
            for record in records
        ),
        "sub_10di_magnetic_power_fraction_upper_bound_maximum": max(
            record["gates"]["Q011-APP-DI"][
                "observed_sub_10di_magnetic_power_fraction_upper_bound"
            ]
            for record in records
        ),
        "delta_B_rms_over_B0_minimum": min(
            record["gates"]["Q011-APP-DI"]["observed_delta_B_rms_over_B0"]
            for record in records
        ),
        "particle_rg_maximum_over_Ly": max(
            record["gates"]["Q011-APP-RG"]["observed_maximum_over_Ly"]
            for record in records
        ),
        "high_energy_tail_rg_maximum_over_Ly": max(
            record["gates"]["Q011-APP-RG"]["observed_high_energy_tail_maximum_over_Ly"]
            for record in records
        ),
        "maximum_particle_specific_kinetic_energy": max(
            record["particle_gyroradius_containment"]["maximum_specific_kinetic_energy"]
            for record in records
        ),
    }
    _require(
        extrema["R_maximum"] >= snapshot_extrema["R_maximum"]
        and extrema["Lambda_maximum"] >= snapshot_extrema["Lambda_maximum"]
        and extrema["S_delta_minimum_excluding_shock_transition"]
        <= snapshot_extrema["S_delta_minimum_excluding_shock_transition"]
        and extrema["lambda_B_characteristic_over_local_di_maximum_minimum"]
        <= snapshot_extrema["lambda_B_characteristic_over_local_di_maximum_minimum"]
        and extrema["sub_10di_magnetic_power_fraction_upper_bound_maximum"]
        >= snapshot_extrema["sub_10di_magnetic_power_fraction_upper_bound_maximum"]
        and extrema["delta_B_rms_over_B0_minimum"]
        <= snapshot_extrema["delta_B_rms_over_B0_minimum"]
        and extrema["particle_rg_maximum_over_Ly"]
        >= snapshot_extrema["particle_rg_maximum_over_Ly"]
        and extrema["high_energy_tail_rg_maximum_over_Ly"]
        >= snapshot_extrema["high_energy_tail_rg_maximum_over_Ly"]
        and extrema["maximum_particle_specific_kinetic_energy"]
        >= snapshot_extrema["maximum_particle_specific_kinetic_energy"],
        "all-cycle runtime extrema do not conservatively dominate supplied snapshots",
    )
    snapshot_pass = {
        name: all(record["gates"][name]["pass"] for record in records)
        for name in (
            "Q011-APP-NORM",
            "Q011-APP-R",
            "Q011-APP-LAMBDA",
            "Q011-APP-DI",
            "Q011-APP-RG",
        )
    }
    escape_fractions = runtime["boundary_escape_ledger"]["escape_fractions"]
    claim_escape = runtime["claim_specific_escape_applicability"]
    gates: dict[str, dict[str, Any]] = {
        "Q011-APP-NORM": {
            "pass": snapshot_pass["Q011-APP-NORM"],
            "basis": "bound_runtime_deck_source_executable_and_raw_snapshot_provenance",
        },
        "Q011-APP-R": {
            "pass": snapshot_pass["Q011-APP-R"] and extrema["R_maximum"] <= R_MAXIMUM,
            "all_cycle_observed_maximum": extrema["R_maximum"],
            "required_maximum": R_MAXIMUM,
        },
        "Q011-APP-LAMBDA": {
            "pass": (
                snapshot_pass["Q011-APP-LAMBDA"]
                and extrema["Lambda_maximum"] <= LAMBDA_MAXIMUM
            ),
            "all_cycle_observed_maximum": extrema["Lambda_maximum"],
            "required_maximum": LAMBDA_MAXIMUM,
        },
        "Q011-APP-DI": {
            "pass": (
                snapshot_pass["Q011-APP-DI"]
                and extrema["S_delta_minimum_excluding_shock_transition"]
                >= S_DELTA_MINIMUM
                and extrema["lambda_B_characteristic_over_local_di_maximum_minimum"]
                >= LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM
                and extrema["sub_10di_magnetic_power_fraction_upper_bound_maximum"]
                <= SUB_10DI_POWER_FRACTION_MAXIMUM
                and extrema["delta_B_rms_over_B0_minimum"]
                >= DELTA_B_RMS_OVER_B0_MINIMUM
            ),
            "all_cycle_extrema": {
                key: extrema[key]
                for key in (
                    "S_delta_minimum_excluding_shock_transition",
                    "lambda_B_characteristic_over_local_di_maximum_minimum",
                    "sub_10di_magnetic_power_fraction_upper_bound_maximum",
                    "delta_B_rms_over_B0_minimum",
                )
            },
        },
        "Q011-APP-RG": {
            "pass": (
                snapshot_pass["Q011-APP-RG"]
                and runtime["escaped_particle_applicability_complete"]
                and extrema["particle_rg_maximum_over_Ly"] <= RG_MAXIMUM_OVER_LY_MAXIMUM
                and extrema["high_energy_tail_rg_maximum_over_Ly"]
                <= RG_HIGH_ENERGY_MAXIMUM_OVER_LY_MAXIMUM
                and escape_fractions["particle_count"]
                <= ESCAPED_PARTICLE_COUNT_FRACTION_MAXIMUM
                and escape_fractions["macro_weight"]
                <= ESCAPED_MACRO_WEIGHT_FRACTION_MAXIMUM
                and escape_fractions["kinetic_energy"]
                <= ESCAPED_KINETIC_ENERGY_FRACTION_MAXIMUM
                and all(applicability["pass"] for applicability in claim_escape.values())
            ),
            "all_cycle_extrema": {
                key: extrema[key]
                for key in (
                    "particle_rg_maximum_over_Ly",
                    "high_energy_tail_rg_maximum_over_Ly",
                    "maximum_particle_specific_kinetic_energy",
                    "escaped_particle_rg_maximum_over_Ly",
                    "escaped_high_energy_tail_rg_maximum_over_Ly",
                    "escaped_particle_specific_kinetic_energy_maximum",
                )
            },
            "escape_fractions": escape_fractions,
            "required_escape_fraction_maxima": {
                "particle_count": ESCAPED_PARTICLE_COUNT_FRACTION_MAXIMUM,
                "macro_weight": ESCAPED_MACRO_WEIGHT_FRACTION_MAXIMUM,
                "kinetic_energy": ESCAPED_KINETIC_ENERGY_FRACTION_MAXIMUM,
            },
            "claim_specific_escape_applicability": claim_escape,
        },
        "Q011-APP-TIME": {
            "pass": runtime["time_and_escape_complete"],
            "basis": (
                "bound_actual_complete_contiguous_per_cycle_inventory_and_closed_"
                "count_weight_energy_momentum_escape_ledger"
            ),
        },
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": HISTORY_RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "authorization": dict(AUTHORIZATION),
        "snapshot_count": len(records),
        "first_observed_committed_time": times[0],
        "last_observed_committed_time": times[-1],
        "runtime_time_escape_evidence": runtime,
        "snapshot_extrema": snapshot_extrema,
        "gates": gates,
        "all_physical_applicability_gates_pass": all(gate["pass"] for gate in gates.values()),
        "claim_rejection_rules": {
            gate: list(claims) for gate, claims in CLAIM_REJECTION_RULES.items()
        },
        "claim_rejections": _claim_rejections(
            gates, extrema["Lambda_maximum"], claim_escape
        ),
        "permanent_claim_exclusions": list(PERMANENT_CLAIM_EXCLUSIONS),
        "authorization_effect": "none_even_if_all_gates_pass",
    }


__all__ = [
    "AUTHORIZATION",
    "ApplicabilitySnapshot",
    "CELL_MAP_NAMES",
    "CLAIM_REJECTION_RULES",
    "CLAIM_SPECIFIC_ESCAPE_BOUNDS",
    "CYCLE_TELEMETRY_RECORD_TYPE",
    "DELTA_B_RMS_OVER_B0_MINIMUM",
    "DETECTED_FRONT_REGIONS",
    "EXACT_NORMALIZATION",
    "ESCAPE_LEDGER_RELATIVE_RESIDUAL_MAXIMUM",
    "ESCAPED_KINETIC_ENERGY_FRACTION_MAXIMUM",
    "ESCAPED_MACRO_WEIGHT_FRACTION_MAXIMUM",
    "ESCAPED_PARTICLE_COUNT_FRACTION_MAXIMUM",
    "EXPECTED_TERMINAL_TIME",
    "HISTORY_RECORD_TYPE",
    "INNER_X1_ESCAPE_REASON",
    "LAMBDA_B_CHAR_OVER_DI_MAX_MINIMUM",
    "LAMBDA_MAXIMUM",
    "MAXIMUM_RUNTIME_CYCLE_SPAN",
    "MINIMUM_RUNTIME_CYCLE_INVENTORY_COUNT",
    "NORMALIZATION_RECORD_TYPE",
    "OUTER_X1_ESCAPE_REASON",
    "PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE",
    "PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES",
    "PERMANENT_CLAIM_EXCLUSIONS",
    "PS_CR_LEDGER_SCHEMA",
    "PS_ESCAPE_ACCOUNTING_SOURCE_COMMIT",
    "PS_ESCAPE_LEDGER_SCHEMA",
    "PhysicalApplicabilityError",
    "QUALIFICATION_EFFECT",
    "REQUIRED_PS_ESCAPE_CHECKPOINT_NOMINAL_TIMES",
    "REQUIRED_RAW_PRODUCTS",
    "RG_HIGH_ENERGY_MAXIMUM_OVER_LY_MAXIMUM",
    "RG_ENERGY_FRACTION_ABOVE_LY_OVER_4_MAXIMUM",
    "RG_MAXIMUM_OVER_LY_MAXIMUM",
    "RG_Q999_OVER_LY_MAXIMUM",
    "RUNTIME_RECORD_TYPE",
    "R_MAXIMUM",
    "SCHEMA_VERSION",
    "SNAPSHOT_RECORD_TYPE",
    "SNAPSHOT_PROVENANCE_RECORD_TYPE",
    "SLOPE_CUTOFF_BIN_ESCAPE_FRACTION_MAXIMUM",
    "SLOPE_CUTOFF_ESCAPE_RECORD_TYPE",
    "SLOPE_CUTOFF_TAIL_ESCAPE_FRACTION_MAXIMUM",
    "SOURCE_MANIFEST_RECORD_TYPE",
    "STARTUP_REMOVAL_TIME",
    "SUB_10DI_POWER_FRACTION_MAXIMUM",
    "SUCCESSOR_ID",
    "S_DELTA_MINIMUM",
    "TRUSTED_Q011_RUNTIME_SOURCE_COMMIT",
    "reduce_physical_applicability_history",
    "reduce_physical_applicability_snapshot",
]
