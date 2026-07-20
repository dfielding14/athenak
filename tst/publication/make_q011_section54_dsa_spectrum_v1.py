#!/usr/bin/env python3
"""Build the Q011 Section 5.4 downstream DSA spectrum comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Mapping, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

if __package__:
    from . import analyze_q011_section54_outputs as output_primitives
    from . import pvtk_particles as pvtk_module
    from . import q011_section54_model as model
    from . import q011_section54_particles as particle_primitives
    from .pvtk_particles import ParticleVTKData, read_particle_vtk
else:
    import analyze_q011_section54_outputs as output_primitives
    import pvtk_particles as pvtk_module
    import q011_section54_model as model
    import q011_section54_particles as particle_primitives
    from pvtk_particles import ParticleVTKData, read_particle_vtk


SCHEMA_VERSION = 1
RECORD_TYPE = "q011_section54_dsa_spectrum_metrics_v1"
MANIFEST_RECORD_TYPE = "q011_section54_dsa_spectrum_manifest_v1"
MHD_FIELDS = ("dens", "velx", "vely", "velz", "eint", "bcc1", "bcc2", "bcc3")
PVTK_SCALARS = frozenset(
    {
        "gid",
        "ptag",
        "species",
        "cr_source",
        "macro_weight",
        "birth_time",
        "deltaf_f0",
        "deltaf_weight",
    }
)
PVTK_VECTORS = frozenset({"vel"})
NOMINAL_TIMES = (500.0, 1200.0)
EARLY_MAXIMUM_LATENESS = 0.1
DEFAULT_ENERGY_BIN_COUNT = 64
DEFAULT_MINIMUM_PARTICLES = 1000
TAIL_FIT_CHI_WINDOW = particle_primitives.LATE_SLOPE_FIT_WINDOW
TAIL_FIT_MINIMUM_POSITIVE_BINS = particle_primitives.LATE_SLOPE_MINIMUM_POSITIVE_BINS
REFERENCE_F_CHI_EXPONENT = -1.5
# Retained for callers that imported the original public constant.
REFERENCE_F_EPSILON_EXPONENT = REFERENCE_F_CHI_EXPONENT
PARTICLE_MACRO_MASS = 9.0e-4
CANONICAL_QUALIFICATION_EFFECT = "derived_figure_only_no_claim_closure"
EXPLORATORY_QUALIFICATION_EFFECT = (
    "exploratory_reduced_ppc_noncanonical_deposit_qscale_no_claim_closure"
)
_MUTABLE_OUTPUT_PARAMETERS = frozenset({"file_number", "last_time"})
_MUTABLE_EXECUTION_PARAMETERS = frozenset({("time", "tlim")})
_SAFE_STEM = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}")
_PVTK_HEADER = re.compile(
    rb"\A# vtk DataFile Version 2[.]0\n"
    rb"# AthenaK particle data at time= ([^ \n]+)  nranks= (0|[1-9][0-9]*)  "
    rb"cycle=(0|[1-9][0-9]*)  variables=([A-Za-z0-9_]+)\n"
)

_EXPECTED_TEXT_PARAMETERS = {
    ("time", "integrator"): model.INTEGRATOR,
    ("mhd", "eos"): "ideal",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_interp_scheme"): "tsc",
    ("particles", "pic_cr_initial_state"): model.DECK_INITIAL_STATE,
    ("particles", "pic_deltaf_mode"): "off",
    ("problem", "pgen_name"): "pic_parallel_shock",
    ("problem", "ps_shock_speed_model"): "ideal_surface",
}
_EXPECTED_BOOLEAN_PARAMETERS = {
    ("particles", "deposit_moments"): True,
    ("particles", "couple_moments_to_mhd"): True,
    ("particles", "couple_moments_momentum_to_mhd"): True,
    ("particles", "couple_moments_energy_to_mhd"): True,
    ("problem", "ps_enable_injection"): True,
    ("problem", "ps_enable_gas_subtraction"): True,
}
_EXPECTED_INTEGER_PARAMETERS = {
    ("particles", "nspecies"): 1,
    ("particles", "deposit_order"): 2,
    ("problem", "ps_inject_species"): 0,
}
_EXPECTED_FLOAT_PARAMETERS = {
    ("mhd", "gamma"): 5.0 / 3.0,
    ("particles", "couple_moments_momentum_coeff"): 1.0,
    ("particles", "couple_moments_energy_coeff"): 1.0,
    ("particles", "pic_cr_light_speed"): model.LIGHT_SPEED,
    ("species0", "mass"): 1.0,
    ("problem", "ps_rho0"): 1.0,
    ("problem", "ps_p0"): 1.0,
    ("problem", "ps_u0"): model.UPSTREAM_SPEED_U0,
    ("problem", "ps_b0"): 1.0,
    ("problem", "ps_eta"): 1.0e-3,
    ("problem", "ps_vinj_over_u0"): math.sqrt(10.0),
    ("problem", "ps_remove_birth_time_before"): 45.0,
}


class DSASpectrumError(ValueError):
    """Reject incomplete, mismatched, or non-finite spectrum inputs."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise DSASpectrumError(message)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": _sha256(path),
        "byte_count": path.stat().st_size,
    }


def _canonical_json(value: Mapping[str, object]) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _validated_expected_deposit_qscale(value: object) -> float:
    _require(
        not isinstance(value, (bool, np.bool_)),
        "expected deposit_qscale must be a finite positive scalar",
    )
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise DSASpectrumError(
            "expected deposit_qscale must be a finite positive scalar"
        ) from error
    _require(
        math.isfinite(result) and result > 0.0,
        "expected deposit_qscale must be a finite positive scalar",
    )
    return result


def _qualification_effect(expected_deposit_qscale: float) -> str:
    if expected_deposit_qscale == PARTICLE_MACRO_MASS:
        return CANONICAL_QUALIFICATION_EFFECT
    return EXPLORATORY_QUALIFICATION_EFFECT


def _parameter_text(
    dataset: output_primitives.AthenaBinaryDataset, block: str, name: str
) -> str:
    try:
        value = dataset.input_parameters[block][name]
    except KeyError as error:
        raise DSASpectrumError(f"mhd_w_bcc header is missing {block}/{name}") from error
    _require(isinstance(value, str) and bool(value), f"{block}/{name} is invalid")
    return value


def _parameter_float(
    dataset: output_primitives.AthenaBinaryDataset, block: str, name: str
) -> float:
    value = _parameter_text(dataset, block, name)
    try:
        result = float(value)
    except ValueError as error:
        raise DSASpectrumError(f"{block}/{name} is not a real scalar") from error
    _require(math.isfinite(result), f"{block}/{name} is non-finite")
    return result


def _parameter_integer(
    dataset: output_primitives.AthenaBinaryDataset, block: str, name: str
) -> int:
    value = _parameter_text(dataset, block, name)
    _require(re.fullmatch(r"-?(?:0|[1-9][0-9]*)", value) is not None,
             f"{block}/{name} is not an integer")
    return int(value)


def _parameter_boolean(
    dataset: output_primitives.AthenaBinaryDataset, block: str, name: str
) -> bool:
    value = _parameter_text(dataset, block, name).lower()
    _require(value in {"true", "false", "1", "0"}, f"{block}/{name} is not boolean")
    return value in {"true", "1"}


def _validate_runtime_model(
    dataset: output_primitives.AthenaBinaryDataset,
    expected_deposit_qscale: float,
) -> dict[str, object]:
    _require(
        set(dataset.variable_names) == set(MHD_FIELDS),
        "mhd_w_bcc field inventory is incomplete or contains unexpected fields",
    )
    for (block, name), expected in _EXPECTED_TEXT_PARAMETERS.items():
        _require(
            _parameter_text(dataset, block, name) == expected,
            f"{block}/{name} is not the Q011 Section 5.4 value {expected!r}",
        )
    for (block, name), expected in _EXPECTED_BOOLEAN_PARAMETERS.items():
        _require(
            _parameter_boolean(dataset, block, name) is expected,
            f"{block}/{name} does not enable the required Q011 coupling",
        )
    for (block, name), expected in _EXPECTED_INTEGER_PARAMETERS.items():
        _require(
            _parameter_integer(dataset, block, name) == expected,
            f"{block}/{name} is not the Q011 Section 5.4 value {expected}",
        )
    measured_floats: dict[str, float] = {}
    for (block, name), expected in _EXPECTED_FLOAT_PARAMETERS.items():
        measured = _parameter_float(dataset, block, name)
        _require(
            math.isclose(measured, expected, rel_tol=1.0e-10, abs_tol=1.0e-12),
            f"{block}/{name} is not the Q011 Section 5.4 value {expected!r}",
        )
        measured_floats[f"{block}/{name}"] = measured
    observed_deposit_qscale = _parameter_float(
        dataset, "particles", "deposit_qscale"
    )
    _require(
        observed_deposit_qscale == expected_deposit_qscale,
        "particles/deposit_qscale does not exactly match the explicitly requested "
        f"value {expected_deposit_qscale!r}; observed {observed_deposit_qscale!r}",
    )
    return {
        "integrator": _parameter_text(dataset, "time", "integrator"),
        "particle_light_speed": measured_floats["particles/pic_cr_light_speed"],
        "upstream_speed_u0": measured_floats["problem/ps_u0"],
        "particle_macro_mass": observed_deposit_qscale,
        "injection_efficiency": measured_floats["problem/ps_eta"],
        "injection_momentum_over_u0": measured_floats["problem/ps_vinj_over_u0"],
        "birth_time_minimum": measured_floats[
            "problem/ps_remove_birth_time_before"
        ],
        "backreaction": {
            "deposit_moments": True,
            "couple_moments_to_mhd": True,
            "momentum_feedback": True,
            "energy_feedback": True,
            "momentum_feedback_coefficient": measured_floats[
                "particles/couple_moments_momentum_coeff"
            ],
            "energy_feedback_coefficient": measured_floats[
                "particles/couple_moments_energy_coeff"
            ],
            "gas_injection_subtraction": True,
        },
    }


def _normalized_input_parameters(
    dataset: output_primitives.AthenaBinaryDataset,
) -> dict[str, dict[str, str]]:
    normalized: dict[str, dict[str, str]] = {}
    for block, parameters in dataset.input_parameters.items():
        normalized[block] = {
            name: value
            for name, value in parameters.items()
            if not (block.startswith("output") and name in _MUTABLE_OUTPUT_PARAMETERS)
            and (block, name) not in _MUTABLE_EXECUTION_PARAMETERS
        }
    return normalized


def _model_identity(dataset: output_primitives.AthenaBinaryDataset) -> dict[str, object]:
    parameters = _normalized_input_parameters(dataset)
    contract: dict[str, object] = {
        "normalized_input_parameters": parameters,
        "root_grid_shape": list(dataset.root_grid_shape),
        "meshblock_shape": list(dataset.meshblock_shape),
        "nghost": int(dataset.nghost),
        "domain_bounds": list(dataset.domain_bounds),
        "location_size": int(dataset.location_size),
        "variable_size": int(dataset.variable_size),
        "variable_names": list(dataset.variable_names),
    }
    digest = hashlib.sha256(_canonical_json(contract).encode("utf-8")).hexdigest()
    return {"sha256": digest, "contract": contract}


def _particle_header(path: Path) -> dict[str, object]:
    with path.open("rb") as stream:
        prefix = stream.read(4096)
    match = _PVTK_HEADER.match(prefix)
    _require(match is not None, "prtcl_all PVTK execution header is malformed")
    assert match is not None
    try:
        time = float(match.group(1))
        nranks = int(match.group(2))
        cycle = int(match.group(3))
        variables = match.group(4).decode("ascii")
    except (UnicodeDecodeError, ValueError) as error:
        raise DSASpectrumError("prtcl_all PVTK execution header is invalid") from error
    _require(math.isfinite(time) and time >= 0.0, "prtcl_all time is invalid")
    _require(nranks > 0 and cycle >= 0, "prtcl_all rank or cycle metadata is invalid")
    _require(variables == "prtcl_all", "particle product is not prtcl_all")
    return {"time": time, "nranks": nranks, "cycle": cycle, "variables": variables}


def _validate_particles(data: ParticleVTKData) -> dict[str, np.ndarray]:
    _require(set(data.scalars) == PVTK_SCALARS, "prtcl_all scalar inventory drifted")
    _require(set(data.vectors) == PVTK_VECTORS, "prtcl_all vector inventory drifted")
    points = np.asarray(data.points, dtype=np.float64)
    velocity = np.asarray(data.vectors["vel"], dtype=np.float64)
    _require(
        points.ndim == 2 and points.shape[1:] == (3,) and velocity.shape == points.shape,
        "prtcl_all points and velocity must have shape (nparticle, 3)",
    )
    count = points.shape[0]
    _require(count > 0, "prtcl_all contains no particles")
    _require(np.all(np.isfinite(points)) and np.all(np.isfinite(velocity)),
             "prtcl_all points or velocity are non-finite")
    _require(
        np.array_equal(points, points.astype(np.float32).astype(np.float64))
        and np.array_equal(velocity, velocity.astype(np.float32).astype(np.float64)),
        "prtcl_all positions and velocities must be exact decoded float32 values",
    )
    arrays = {"points": points, "velocity": velocity}
    for name, raw in data.scalars.items():
        values = np.asarray(raw)
        _require(values.shape == (count,), f"prtcl_all scalar {name} shape drifted")
        _require(np.all(np.isfinite(values)), f"prtcl_all scalar {name} is non-finite")
        arrays[name] = values
    _require(
        arrays["gid"].dtype.kind in "iu"
        and arrays["ptag"].dtype.kind in "iu"
        and arrays["species"].dtype.kind in "iu"
        and arrays["cr_source"].dtype.kind in "iu",
        "prtcl_all identity and provenance fields must be decoded integers",
    )
    _require(
        np.all(arrays["gid"] >= 0)
        and np.all(arrays["ptag"] >= 0)
        and np.unique(arrays["ptag"]).size == count,
        "prtcl_all gid/ptag identity is invalid",
    )
    _require(np.all(arrays["species"] == 0), "prtcl_all species inventory is not Q011")
    _require(np.all(arrays["macro_weight"] >= 0.0), "prtcl_all macro weights are negative")
    return arrays


def _weighted_census(mask: np.ndarray, weights: np.ndarray) -> dict[str, object]:
    return {
        "particle_count": int(np.count_nonzero(mask)),
        "macro_weight": float(np.sum(weights[mask])),
    }


def _filter_particles(
    arrays: Mapping[str, np.ndarray], snapshot_time: float
) -> tuple[np.ndarray, dict[str, object]]:
    source = arrays["cr_source"] == particle_primitives.SHOCK_INJECTED_SOURCE
    early = source & (arrays["birth_time"] < particle_primitives.BIRTH_TIME_MINIMUM)
    provenance = source & ~early
    nonpositive = provenance & (arrays["macro_weight"] <= 0.0)
    positive = provenance & ~nonpositive
    surface = model.x_ideal(snapshot_time)
    upstream = positive & (arrays["points"][:, 0] > surface)
    on_surface = positive & (arrays["points"][:, 0] == surface)
    selected = positive & (arrays["points"][:, 0] < surface)
    wrong_source = ~source
    partitions = np.stack(
        (wrong_source, early, nonpositive, upstream, on_surface, selected), axis=0
    )
    _require(
        np.all(np.sum(partitions, axis=0) == 1),
        "internal particle-filter census is not a disjoint complete partition",
    )
    weights = np.asarray(arrays["macro_weight"], dtype=np.float64)
    census = {
        "all_particles": _weighted_census(np.ones(selected.size, dtype=bool), weights),
        "rejected_wrong_source": _weighted_census(wrong_source, weights),
        "rejected_early_birth_time": _weighted_census(early, weights),
        "rejected_nonpositive_macro_weight": _weighted_census(nonpositive, weights),
        "rejected_upstream": _weighted_census(upstream, weights),
        "rejected_on_ideal_surface": _weighted_census(on_surface, weights),
        "admitted_downstream": _weighted_census(selected, weights),
    }
    return selected, {
        "selection": [
            "cr_source == 1",
            "birth_time >= 45",
            "macro_weight > 0",
            "x1 < x_ideal(t)",
        ],
        "ideal_surface_x1_c_over_omega_pi": surface,
        "disjoint_census": census,
    }


def _validate_nominal_time(observed: float, nominal: float) -> None:
    if nominal == NOMINAL_TIMES[-1]:
        _require(observed == nominal, "terminal t=1200 particle snapshot must be exact")
    else:
        _require(
            observed >= nominal and observed <= nominal + EARLY_MAXIMUM_LATENESS,
            "t=500 particle snapshot exceeds the frozen nominal-slot lateness bound",
        )


def _specific_energy_from_chi(
    chi: object, upstream_speed: float, light_speed: float
) -> np.ndarray:
    values = np.asarray(chi, dtype=np.float64)
    momentum_squared = values * upstream_speed**2
    return momentum_squared / (
        np.sqrt(1.0 + momentum_squared / light_speed**2) + 1.0
    )


def _read_snapshot(
    mhd_path: Path,
    particle_path: Path,
    nominal_time: float,
    expected_deposit_qscale: float,
) -> dict[str, object]:
    try:
        dataset = output_primitives.read_athenak_binary(mhd_path)
    except (OSError, output_primitives.AnalysisError) as error:
        raise DSASpectrumError(f"unable to read mhd_w_bcc: {error}") from error
    runtime_model = _validate_runtime_model(dataset, expected_deposit_qscale)
    header = _particle_header(particle_path)
    _require(dataset.cycle == header["cycle"], "mhd_w_bcc and prtcl_all cycles disagree")
    projected_time = float(format(float(header["time"]), ".6g"))
    _require(
        dataset.time == projected_time,
        "mhd_w_bcc time is not the six-significant-digit projection of PVTK time",
    )
    observed_time = float(header["time"])
    _validate_nominal_time(observed_time, nominal_time)
    try:
        particle_data = read_particle_vtk(particle_path)
    except (OSError, ValueError) as error:
        raise DSASpectrumError(f"unable to read prtcl_all: {error}") from error
    arrays = _validate_particles(particle_data)
    selected, particle_filter = _filter_particles(arrays, observed_time)
    _require(np.any(selected), "prtcl_all has no admitted positive-weight downstream particles")
    try:
        chi = particle_primitives.reconstruct_chi_from_pvtk_velocity(
            arrays["velocity"][selected]
        )
    except particle_primitives.ParticleReducerError as error:
        raise DSASpectrumError(f"particle momentum reconstruction failed: {error}") from error
    u0 = float(runtime_model["upstream_speed_u0"])
    light_speed = float(runtime_model["particle_light_speed"])
    specific_energy = _specific_energy_from_chi(chi, u0, light_speed)
    weights = np.asarray(arrays["macro_weight"][selected], dtype=np.float64)
    _require(
        np.all(np.isfinite(specific_energy))
        and np.all(specific_energy >= 0.0)
        and np.all(np.isfinite(weights))
        and np.all(weights > 0.0),
        "admitted particle energies must be finite and non-negative; macro weights "
        "must be finite and positive",
    )
    return {
        "dataset": dataset,
        "model_identity": _model_identity(dataset),
        "runtime_model": runtime_model,
        "snapshot": {
            "nominal_time_omega0_inverse": nominal_time,
            "particle_time_omega0_inverse": observed_time,
            "mhd_time_omega0_inverse": float(dataset.time),
            "cycle": int(dataset.cycle),
            "particle_nranks": int(header["nranks"]),
            "mesh_time_projection": projected_time,
        },
        "particle_filter": particle_filter,
        "specific_energy": specific_energy,
        "macro_weight": weights,
        "chi": chi,
    }


def _tail_fit(fixed_spectrum: Mapping[str, object]) -> dict[str, object]:
    edges = np.asarray(fixed_spectrum["bin_edges"], dtype=np.float64)
    centers = np.sqrt(edges[:-1] * edges[1:])
    f_chi = np.asarray(fixed_spectrum["f_chi"], dtype=np.float64)
    selected = (
        (centers >= TAIL_FIT_CHI_WINDOW[0])
        & (centers <= TAIL_FIT_CHI_WINDOW[1])
        & np.isfinite(f_chi)
        & (f_chi > 0.0)
    )
    indices = np.flatnonzero(selected)
    record: dict[str, object] = {
        "fit_window_chi": list(TAIL_FIT_CHI_WINDOW),
        "minimum_positive_bins": TAIL_FIT_MINIMUM_POSITIVE_BINS,
        "selected_bin_indices": indices.tolist(),
        "positive_bin_count": int(indices.size),
        "evaluated": bool(indices.size >= TAIL_FIT_MINIMUM_POSITIVE_BINS),
        "slope_f_chi": None,
        "intercept_log_f_chi": None,
        "r_squared": None,
        "reducer_record": None,
    }
    if indices.size < TAIL_FIT_MINIMUM_POSITIVE_BINS:
        return record

    try:
        reduced = particle_primitives.late_slope_record(f_chi)
    except particle_primitives.ParticleReducerError as error:
        raise DSASpectrumError(f"fixed chi tail fit failed: {error}") from error
    x = np.log(centers[selected])
    y = np.log(f_chi[selected])
    slope = float(reduced["slope"])
    intercept = float(reduced["intercept"])
    predicted = slope * x + intercept
    residual = float(np.sum((y - predicted) ** 2))
    total = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 if total == 0.0 else 1.0 - residual / total
    record.update(
        {
            "slope_f_chi": slope,
            "intercept_log_f_chi": intercept,
            "r_squared": float(r_squared),
            "reducer_record": reduced,
        }
    )
    return record


def _spectrum_record(snapshot: Mapping[str, object]) -> dict[str, object]:
    chi = np.asarray(snapshot["chi"], dtype=np.float64)
    energy = np.asarray(snapshot["specific_energy"], dtype=np.float64)
    weights = np.asarray(snapshot["macro_weight"], dtype=np.float64)
    try:
        fixed_spectrum = particle_primitives.weighted_spectrum_record(chi, weights)
    except particle_primitives.ParticleReducerError as error:
        raise DSASpectrumError(f"fixed chi spectrum reduction failed: {error}") from error

    edges = np.asarray(fixed_spectrum["bin_edges"], dtype=np.float64)
    frozen_edges = np.asarray(particle_primitives.CHI_BIN_EDGES, dtype=np.float64)
    _require(
        np.array_equal(edges, frozen_edges),
        "particle spectrum reducer did not return the frozen chi bin edges",
    )
    centers = np.sqrt(edges[:-1] * edges[1:])
    weighted_counts = np.asarray(fixed_spectrum["weighted_counts"], dtype=np.float64)
    total_weight = float(fixed_spectrum["total_post_filter_macro_weight"])
    in_bin_count = int(fixed_spectrum["particle_count_in_bins"])
    in_bin_weight = float(fixed_spectrum["macro_weight_in_bins"])
    underflow_count = int(fixed_spectrum["underflow_count"])
    underflow_weight = float(fixed_spectrum["underflow_macro_weight"])
    overflow_count = int(fixed_spectrum["overflow_count"])
    overflow_weight = float(fixed_spectrum["overflow_macro_weight"])
    accounted_count = in_bin_count + underflow_count + overflow_count
    accounted_weight = in_bin_weight + underflow_weight + overflow_weight
    _require(
        accounted_count == chi.size
        and math.isclose(
            accounted_weight, total_weight, rel_tol=1.0e-12, abs_tol=1.0e-12
        ),
        "fixed chi histogram underflow/in-bin/overflow accounting does not close",
    )

    u0 = float(snapshot["runtime_model"]["upstream_speed_u0"])
    light_speed = float(snapshot["runtime_model"]["particle_light_speed"])
    energy_edges = _specific_energy_from_chi(edges, u0, light_speed)
    energy_centers = _specific_energy_from_chi(centers, u0, light_speed)
    f_epsilon = weighted_counts / total_weight / np.diff(energy_edges)
    epsilon2_f_epsilon = energy_centers**2 * f_epsilon
    fit_energy = _specific_energy_from_chi(TAIL_FIT_CHI_WINDOW, u0, light_speed)
    tail_fit = _tail_fit(fixed_spectrum)
    tail_fit["fit_window_specific_energy"] = fit_energy.tolist()
    macro_mass = float(snapshot["runtime_model"]["particle_macro_mass"])
    return {
        "snapshot": snapshot["snapshot"],
        "particle_filter": snapshot["particle_filter"],
        "weighted_spectrum": fixed_spectrum,
        "chi_bin_edges": fixed_spectrum["bin_edges"],
        "chi_bin_centers": centers.tolist(),
        "f_chi": fixed_spectrum["f_chi"],
        "normalized_chi_f_chi": fixed_spectrum["normalized_chi_f_chi"],
        "specific_energy_bin_edges": energy_edges.tolist(),
        "specific_energy_bin_centers": energy_centers.tolist(),
        "counts": fixed_spectrum["counts"],
        "macro_weighted_counts": fixed_spectrum["weighted_counts"],
        "f_epsilon": f_epsilon.tolist(),
        "epsilon_squared_f_epsilon": epsilon2_f_epsilon.tolist(),
        "selected_particle_count": int(chi.size),
        "selected_macro_weight": total_weight,
        "selected_physical_macro_mass": macro_mass * total_weight,
        "selected_physical_kinetic_energy": float(macro_mass * np.sum(weights * energy)),
        "selected_specific_energy_minimum": float(np.min(energy)),
        "selected_specific_energy_maximum": float(np.max(energy)),
        "histogrammed_particle_count": in_bin_count,
        "histogrammed_macro_weight": in_bin_weight,
        "underflow_count": underflow_count,
        "underflow_macro_weight": underflow_weight,
        "underflow_macro_weight_fraction": underflow_weight / total_weight,
        "overflow_count": overflow_count,
        "overflow_macro_weight": overflow_weight,
        "overflow_macro_weight_fraction": fixed_spectrum[
            "overflow_macro_weight_fraction"
        ],
        "overflow_gate": fixed_spectrum["overflow_gate"],
        "accounted_particle_count": accounted_count,
        "accounted_macro_weight": accounted_weight,
        "in_bin_macro_weight_fraction": in_bin_weight / total_weight,
        "histogram_closure_fraction": accounted_weight / total_weight,
        "tail_fit": tail_fit,
    }


def _reference_record(late_spectrum: Mapping[str, object]) -> dict[str, object]:
    centers = np.asarray(late_spectrum["chi_bin_centers"], dtype=np.float64)
    f_chi = np.asarray(late_spectrum["f_chi"], dtype=np.float64)
    displayed = np.asarray(late_spectrum["normalized_chi_f_chi"], dtype=np.float64)
    fit_window = late_spectrum["tail_fit"]["fit_window_chi"]
    candidates = np.flatnonzero(
        (f_chi > 0.0) & (centers >= fit_window[0]) & (centers <= fit_window[1])
    )
    if candidates.size == 0:
        candidates = np.flatnonzero(f_chi > 0.0)
    _require(candidates.size > 0, "late spectrum contains no positive reference bin")
    target = math.sqrt(float(centers[candidates[0]] * centers[candidates[-1]]))
    index = int(candidates[np.argmin(np.abs(np.log(centers[candidates] / target)))])
    return {
        "f_chi_exponent": REFERENCE_F_CHI_EXPONENT,
        "displayed_chi_f_chi_exponent": 1.0 + REFERENCE_F_CHI_EXPONENT,
        "label": "f(chi) proportional to chi^(-3/2)",
        "anchor_bin_index": index,
        "anchor_chi": float(centers[index]),
        "anchor_f_chi": float(f_chi[index]),
        "anchor_normalized_chi_f_chi": float(displayed[index]),
        "line_chi_range": list(fit_window),
        # Compatibility aliases for consumers of the original metrics record.
        "f_epsilon_exponent": REFERENCE_F_CHI_EXPONENT,
        "displayed_epsilon_squared_f_exponent": (
            2.0 + REFERENCE_F_EPSILON_EXPONENT
        ),
        "anchor_specific_energy": late_spectrum["specific_energy_bin_centers"][index],
        "anchor_f_epsilon": late_spectrum["f_epsilon"][index],
        "line_specific_energy_range": late_spectrum["tail_fit"][
            "fit_window_specific_energy"
        ],
    }


def analyze_dsa_spectrum(
    early_mhd_path: Path,
    early_particle_path: Path,
    late_mhd_path: Path,
    late_particle_path: Path,
    *,
    energy_bin_count: int = DEFAULT_ENERGY_BIN_COUNT,
    minimum_particles: int = DEFAULT_MINIMUM_PARTICLES,
    expected_deposit_qscale: float = PARTICLE_MACRO_MASS,
) -> dict[str, object]:
    """Reduce exact t=500 and t=1200 same-model MHD/PVTK pairs."""
    _require(
        isinstance(energy_bin_count, int)
        and not isinstance(energy_bin_count, (bool, np.bool_))
        and energy_bin_count >= 8,
        "energy bin count must be at least 8",
    )
    _require(minimum_particles > 0, "minimum particle count must be positive")
    requested_qscale = _validated_expected_deposit_qscale(expected_deposit_qscale)
    snapshots = (
        _read_snapshot(
            early_mhd_path,
            early_particle_path,
            NOMINAL_TIMES[0],
            requested_qscale,
        ),
        _read_snapshot(
            late_mhd_path,
            late_particle_path,
            NOMINAL_TIMES[1],
            requested_qscale,
        ),
    )
    _require(
        snapshots[0]["model_identity"] == snapshots[1]["model_identity"],
        "t=500 and t=1200 snapshots do not have an exact same-model identity",
    )
    _require(
        snapshots[0]["runtime_model"] == snapshots[1]["runtime_model"],
        "t=500 and t=1200 runtime Q011 model parameters disagree",
    )
    for snapshot in snapshots:
        count = int(np.asarray(snapshot["chi"]).size)
        _require(
            count >= minimum_particles,
            f"downstream spectrum has {count} particles; minimum is {minimum_particles}",
        )
    u0 = float(snapshots[0]["runtime_model"]["upstream_speed_u0"])
    light_speed = float(snapshots[0]["runtime_model"]["particle_light_speed"])
    fit_energy = _specific_energy_from_chi(TAIL_FIT_CHI_WINDOW, u0, light_speed)
    spectra = [_spectrum_record(snapshot) for snapshot in snapshots]
    reference = _reference_record(spectra[-1])
    observed_qscale = float(snapshots[0]["runtime_model"]["particle_macro_mass"])
    qscale_contract = {
        "canonical": PARTICLE_MACRO_MASS,
        "requested": requested_qscale,
        "observed": observed_qscale,
        "exact_match": observed_qscale == requested_qscale,
        "canonical_run": requested_qscale == PARTICLE_MACRO_MASS,
        "explicit_noncanonical_override": requested_qscale != PARTICLE_MACRO_MASS,
    }
    qualification_effect = _qualification_effect(requested_qscale)
    metrics: dict[str, object] = {
        "record_type": RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "qualification_effect": qualification_effect,
        "deposit_qscale": qscale_contract,
        "same_model_validation": {
            "passed": True,
            "method": (
                "exact normalized Athena input parameters, mesh identity, field inventory, "
                "numeric storage widths, and required Q011 full-backreaction runtime values"
            ),
            "mutable_output_parameters_ignored": sorted(_MUTABLE_OUTPUT_PARAMETERS),
            "mutable_execution_parameters_ignored": [
                f"{block}/{name}"
                for block, name in sorted(_MUTABLE_EXECUTION_PARAMETERS)
            ],
            "model_identity_sha256": snapshots[0]["model_identity"]["sha256"],
        },
        "runtime_model": snapshots[0]["runtime_model"],
        "spectrum_definition": {
            "population": (
                "shock-injected, birth_time>=45, positive-macro-weight particles "
                "strictly downstream of the ideal injection surface"
            ),
            "momentum_reconstruction": "p_over_m = gamma(v) * v from PVTK physical velocity",
            "dimensionless_coordinate": "chi = p_over_m_squared / u0_squared",
            "specific_kinetic_energy": (
                "epsilon = p_over_m_squared / "
                "(sqrt(1 + p_over_m_squared / C_squared) + 1)"
            ),
            "distribution": (
                "f_chi = sum(macro_weight in fixed chi bin) / delta_chi"
            ),
            "displayed_quantity": (
                "chi * f_chi / total_post_filter_macro_weight"
            ),
            "displayed_field": "normalized_chi_f_chi",
            "displayed_units": "dimensionless",
            "physical_particle_weight": "deposit_qscale * macro_weight",
            "fixed_dimensionless_chi_bins": True,
            "fixed_bin_source": "q011_section54_particles.CHI_BIN_EDGES",
            "fixed_bin_edges": list(particle_primitives.CHI_BIN_EDGES),
            "common_data_bound_energy_bins": False,
            "legacy_energy_bin_count_argument": {
                "requested": energy_bin_count,
                "effect": "accepted for CLI compatibility; fixed chi bins are authoritative",
            },
            "underflow_policy": "archive count and macro weight; include weight in normalization",
            "overflow_policy": (
                "archive count and macro weight; include weight in normalization and evaluate "
                "the frozen overflow gate"
            ),
            "float32_velocity_caveat": model.FLOAT32_PROJECTION_UNCERTAINTY,
        },
        "tail_fit_definition": {
            "source_chi_window": list(TAIL_FIT_CHI_WINDOW),
            "physical_specific_energy_window": fit_energy.tolist(),
            "minimum_positive_bins": TAIL_FIT_MINIMUM_POSITIVE_BINS,
            "fit_quantity": "q011_section54_particles.late_slope_record on fixed-bin f_chi",
        },
        "reference_power_law": reference,
        "spectra": spectra,
    }
    return {"metrics": metrics}


def _style() -> None:
    plt.style.use("default")
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def render_dsa_spectrum(
    reduction: Mapping[str, object], output_paths: Sequence[Path], dpi: int
) -> None:
    """Render one validated same-model spectrum reduction."""
    _require(dpi > 0, "dpi must be positive")
    _style()
    metrics = reduction["metrics"]
    spectra = metrics["spectra"]
    fig, axis = plt.subplots(figsize=(6.7, 4.8), constrained_layout=True)
    colors = ("#0072B2", "#D55E00")
    styles = ("-", "--")
    for spectrum, color, style in zip(spectra, colors, styles):
        centers = np.asarray(spectrum["chi_bin_centers"], dtype=np.float64)
        displayed = np.asarray(spectrum["normalized_chi_f_chi"], dtype=np.float64)
        positive = displayed > 0.0
        time = spectrum["snapshot"]["particle_time_omega0_inverse"]
        axis.loglog(
            centers[positive],
            displayed[positive],
            color=color,
            ls=style,
            lw=1.8,
            marker="o",
            ms=2.2,
            markevery=max(1, int(np.count_nonzero(positive) / 14)),
            label=rf"$\Omega_0 t={time:.6g}$",
        )
    reference = metrics["reference_power_law"]
    reference_chi = np.geomspace(*reference["line_chi_range"], 100)
    anchor_chi = float(reference["anchor_chi"])
    anchor_displayed = float(reference["anchor_normalized_chi_f_chi"])
    reference_displayed = anchor_displayed * (
        reference_chi / anchor_chi
    ) ** float(reference["displayed_chi_f_chi_exponent"])
    axis.loglog(
        reference_chi,
        reference_displayed,
        color="black",
        ls=":",
        lw=1.4,
        label=(
            r"reference $f(\chi)\propto\chi^{-3/2}$ "
            r"($\chi f\propto\chi^{-1/2}$)"
        ),
    )
    edges = np.asarray(spectra[0]["chi_bin_edges"], dtype=np.float64)
    axis.set_xlim(edges[0], edges[-1])
    axis.set_xlabel(r"dimensionless momentum-energy coordinate $\chi=(p/m)^2/U_0^2$")
    axis.set_ylabel(r"normalized $\chi f(\chi)$")
    axis.set_title("Q011 Section 5.4 downstream DSA spectrum")
    axis.grid(which="major", color="0.88", lw=0.5)
    axis.legend(frameon=False, fontsize=8.5, loc="best")
    axis.text(
        0.02,
        0.03,
        r"$cr\_source=1$, $t_{birth}\geq45$, $w>0$, $x_1<x_{ideal}(t)$",
        transform=axis.transAxes,
        fontsize=8,
    )
    accounting = []
    for spectrum in spectra:
        time = spectrum["snapshot"]["particle_time_omega0_inverse"]
        accounting.append(
            rf"$t={time:.6g}$: under={spectrum['underflow_count']} "
            rf"({100.0 * spectrum['underflow_macro_weight_fraction']:.3g}\%), "
            rf"over={spectrum['overflow_count']} "
            rf"({100.0 * spectrum['overflow_macro_weight_fraction']:.3g}\%)"
        )
    axis.text(
        0.98,
        0.97,
        "\n".join(accounting),
        transform=axis.transAxes,
        fontsize=7,
        ha="right",
        va="top",
    )
    if not metrics["deposit_qscale"]["canonical_run"]:
        axis.text(
            0.98,
            0.03,
            f"exploratory: deposit_qscale={metrics['deposit_qscale']['observed']:.6g}",
            transform=axis.transAxes,
            fontsize=8,
            ha="right",
        )
    for path in output_paths:
        if path.suffix.lower() == ".png":
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
        else:
            fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def make_dsa_spectrum(
    early_mhd_path: Path,
    early_particle_path: Path,
    late_mhd_path: Path,
    late_particle_path: Path,
    output_dir: Path,
    *,
    stem: str = "q011_section54_dsa_spectrum_v1",
    energy_bin_count: int = DEFAULT_ENERGY_BIN_COUNT,
    minimum_particles: int = DEFAULT_MINIMUM_PARTICLES,
    expected_deposit_qscale: float = PARTICLE_MACRO_MASS,
    dpi: int = 300,
) -> list[Path]:
    """Analyze exact inputs and emit PNG, PDF, JSON, and a hashed manifest."""
    _require(_SAFE_STEM.fullmatch(stem) is not None, "output stem is unsafe")
    input_paths = tuple(
        path.resolve(strict=True)
        for path in (
            early_mhd_path,
            early_particle_path,
            late_mhd_path,
            late_particle_path,
        )
    )
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    reduction = analyze_dsa_spectrum(
        *input_paths,
        energy_bin_count=energy_bin_count,
        minimum_particles=minimum_particles,
        expected_deposit_qscale=expected_deposit_qscale,
    )
    png_path = output_dir / f"{stem}.png"
    pdf_path = output_dir / f"{stem}.pdf"
    metrics_path = output_dir / f"{stem}.json"
    manifest_path = output_dir / f"{stem}_manifest.json"
    render_dsa_spectrum(reduction, (png_path, pdf_path), dpi)
    metrics_path.write_text(_canonical_json(reduction["metrics"]), encoding="utf-8")
    dependencies = [
        Path(__file__),
        Path(output_primitives.__file__),
        Path(model.__file__),
        Path(particle_primitives.__file__),
        Path(pvtk_module.__file__),
    ]
    manifest = {
        "record_type": MANIFEST_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "qualification_effect": reduction["metrics"]["qualification_effect"],
        "deposit_qscale": reduction["metrics"]["deposit_qscale"],
        "generator": _artifact(Path(__file__)),
        "source_dependencies": [_artifact(path) for path in dependencies[1:]],
        "inputs": [_artifact(path) for path in input_paths],
        "outputs": [_artifact(path) for path in (png_path, pdf_path, metrics_path)],
        "same_model_validation": reduction["metrics"]["same_model_validation"],
        "snapshots": [item["snapshot"] for item in reduction["metrics"]["spectra"]],
        "analysis_parameters": {
            "energy_bin_count": energy_bin_count,
            "energy_bin_count_effect": (
                "compatibility argument only; the frozen chi bin edges are authoritative"
            ),
            "fixed_chi_bin_count": len(particle_primitives.CHI_BIN_EDGES) - 1,
            "fixed_chi_bin_edges": list(particle_primitives.CHI_BIN_EDGES),
            "minimum_admitted_particles_per_snapshot": minimum_particles,
            "expected_deposit_qscale": reduction["metrics"]["deposit_qscale"][
                "requested"
            ],
            "nominal_times_omega0_inverse": list(NOMINAL_TIMES),
            "early_maximum_lateness_omega0_inverse": EARLY_MAXIMUM_LATENESS,
            "terminal_time_exact": True,
            "particle_filter": [
                "cr_source == 1",
                "birth_time >= 45",
                "macro_weight > 0",
                "x1 < x_ideal(t)",
            ],
            "displayed_quantity": (
                "chi * f_chi / total_post_filter_macro_weight"
            ),
            "reference_f_chi_exponent": REFERENCE_F_CHI_EXPONENT,
            "reference_f_epsilon_exponent": REFERENCE_F_EPSILON_EXPONENT,
            "underflow_overflow_accounting": True,
            "dpi": dpi,
            "cross_cycle_substitution": False,
        },
    }
    manifest_path.write_text(_canonical_json(manifest), encoding="utf-8")
    return [png_path, pdf_path, metrics_path, manifest_path]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--early-mhd", type=Path, required=True)
    parser.add_argument("--early-particles", type=Path, required=True)
    parser.add_argument("--late-mhd", type=Path, required=True)
    parser.add_argument("--late-particles", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stem", default="q011_section54_dsa_spectrum_v1")
    parser.add_argument(
        "--energy-bin-count",
        type=int,
        default=DEFAULT_ENERGY_BIN_COUNT,
        help=(
            "Deprecated compatibility option; validated but ignored because the "
            "frozen dimensionless chi bins are authoritative"
        ),
    )
    parser.add_argument("--minimum-particles", type=int, default=DEFAULT_MINIMUM_PARTICLES)
    parser.add_argument(
        "--expected-deposit-qscale",
        type=float,
        default=PARTICLE_MACRO_MASS,
        help=(
            "Required particles/deposit_qscale; noncanonical values are accepted only "
            "as explicit exploratory overrides"
        ),
    )
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args(argv)
    outputs = make_dsa_spectrum(
        args.early_mhd,
        args.early_particles,
        args.late_mhd,
        args.late_particles,
        args.output_dir,
        stem=args.stem,
        energy_bin_count=args.energy_bin_count,
        minimum_particles=args.minimum_particles,
        expected_deposit_qscale=args.expected_deposit_qscale,
        dpi=args.dpi,
    )
    print(_canonical_json({"outputs": [str(path) for path in outputs]}), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
