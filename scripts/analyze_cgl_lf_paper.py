#!/usr/bin/env python3
"""Generate MKS24-oriented diagnostics from archived CGL-LF snapshots."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import sys
import tempfile
import types

import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR / "vis" / "python"))
sys.path.insert(0, str(ROOT_DIR / "scripts" / "frontier"))
import bin_convert  # noqa: E402
import cgl_lf_stage_i_validate_segment as stage_i_validator  # noqa: E402


HEADER_LABEL = re.compile(r"\[\d+\]=(\S+)")
REQUIRED_FIELDS = ("dens", "velx", "vely", "velz", "eint", "p_perp",
                   "bcc1", "bcc2", "bcc3")
SQRT_TWO_OVER_PI = 0.7978845608028654
SQRT_EIGHT_OVER_PI = 1.5957691216057308
SQRT_TWO_PI = 2.5066282746310002
SQRT_EIGHT_PI = 5.013256549262000
THREE_PI_MINUS_EIGHT = 1.4247779607693793
BACKUP_COLLISION_RATE = 1.0e10
REFERENCE_CURVE_SCHEMA_VERSION = 1
STAGE_I_PANEL_SCHEMA_VERSION = 2
STAGE_I_PANEL_STATUSES = (
    "external_model",
    "blocked_reference",
    "not_run",
    "failed",
    "passed",
)
STAGE_I_PRODUCTION_WORKFLOW = "paper-mks24-stage-i-production"
EDDY_ANGLE_DEGREES = 15.0
FIGURE_13_FIREHOSE_SCHEMA_VERSION = 2
FIGURE_13_FIREHOSE_PRODUCT_NAME = "figure_13_alternate_firehose_occupancy.json"
ANALYSIS_PUBLICATION_SCHEMA_VERSION = 1
ANALYSIS_PUBLICATION_CURRENT_NAME = "analysis-current"
ANALYSIS_PUBLICATION_GENERATIONS_NAME = ".analysis-generations"
ANALYSIS_PUBLICATION_MANIFEST_NAME = "publication-manifest.json"
MECHANISM_DIAGNOSTIC_SCHEMA_VERSION = 1
MECHANISM_JOINT_SCHEMA_VERSION = 1
DESCRIPTIVE_UNCERTAINTY_SCHEMA_VERSION = 1
DESCRIPTIVE_BLOCK_COUNT_MAX = 4
MECHANISM_FIELD_DEFINITIONS = {
    "delta_p": "pressure anisotropy Delta p = p_perp - p_parallel",
    "bb_grad_velocity": (
        "b b : grad(u), with b = B/sqrt(max(B^2, machine tiny))"
    ),
    "div_velocity": "div(u)",
    "dln_b_dt": (
        "D ln(B)/Dt reconstructed from the ideal-induction identity "
        "b b : grad(u) - div(u), with "
        "b = B/sqrt(max(B^2, machine tiny))"
    ),
    "b_grad_delta_p": (
        "b . grad(Delta p), with b = B/sqrt(max(B^2, machine tiny)) and "
        "Delta p = p_perp - p_parallel"
    ),
    "signed_pressure_stress_power_density": (
        "-Delta p (b b : grad(u)); positive values correspond to local "
        "kinetic-energy gain if the anisotropic pressure-stress force is applied"
    ),
}
MECHANISM_JOINT_FIELDS = {
    "b_grad_delta_p_vs_delta_p": ("delta_p", "b_grad_delta_p"),
    "bb_grad_velocity_vs_delta_p": ("delta_p", "bb_grad_velocity"),
    "bb_grad_velocity_vs_b_grad_delta_p": (
        "b_grad_delta_p",
        "bb_grad_velocity",
    ),
    "dln_b_dt_vs_bb_grad_velocity": ("bb_grad_velocity", "dln_b_dt"),
    "signed_pressure_stress_power_density_vs_bb_grad_velocity": (
        "bb_grad_velocity",
        "signed_pressure_stress_power_density",
    ),
    "signed_pressure_stress_power_density_vs_delta_p": (
        "delta_p",
        "signed_pressure_stress_power_density",
    ),
}
FIGURE_13_FIREHOSE_CASES = ("R03", "R07", "R14", "R15")
FIGURE_13_EXECUTION_EPOCH = "E03-forcing-policy"
FIGURE_13_ACCEPTED_FINAL_TIME = 10.0
STAGE_I_CANONICAL_ROOT = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL"
)
STAGE_I_RUNS_RELATIVE = Path("runs/mks24-stage-i")
FIGURE_13_MATRIX_PATH = "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
FIGURE_13_MATRIX_SHA256 = (
    "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
)
FIGURE_13_WINDOW_START = 8.0
FIGURE_13_WINDOW_END = 10.0
FIGURE_13_SNAPSHOT_CADENCE = 0.25
FIGURE_13_TIME_TOLERANCE_FRACTION = 1.0e-6
FIGURE_13_TIME_TOLERANCE = (
    FIGURE_13_SNAPSHOT_CADENCE * FIGURE_13_TIME_TOLERANCE_FRACTION
)
MAX_RESTART_PARAMETER_DUMP_BYTES = 11 * 4096 + 1
FIGURE_13_QUALIFIED_RESTART_LOAD_EVIDENCE_RELATIVE = Path(
    "runs/qualification-9e075422-e03/"
    "g027-e03-modal-rankio-natural-restart-identity-8gpu/analysis/"
    "g027_e03_modal_rankio_restart_identity_evidence.json"
)
FIGURE_13_QUALIFIED_RESTART_LOAD_ROLE = (
    "natural_post_ou_refresh_restart_identity_against_g026"
)
FIGURE_13_QUALIFIED_RESTART_BINARY_ABIS = (
    stage_i_validator.QUALIFIED_RESTART_BINARY_ABIS
)
FIGURE_13_GRID_SHAPE_Z_Y_X = (384, 192, 192)
FIGURE_13_GRID_LENGTHS_X_Y_Z = (1.0, 1.0, 2.0)
FIGURE_13_SNAPSHOT_TIMES = tuple(
    FIGURE_13_WINDOW_START + FIGURE_13_SNAPSHOT_CADENCE * index
    for index in range(
        int(
            round(
                (FIGURE_13_WINDOW_END - FIGURE_13_WINDOW_START)
                / FIGURE_13_SNAPSHOT_CADENCE
            )
        )
        + 1
    )
)
FIGURE_13_CASE_CONTRACTS = {
    "R03": {
        "name": "paper_standard_active_alfvenic_beta100",
        "input": (
            "inputs/cgl_lf_paper/"
            "cgl_lf_paper_standard_active_alfvenic_beta100.athinput"
        ),
        "input_sha256": (
            "997f449abb3c2e4d1de50509efffa127c6230e3ecb2632612200010b8fa6c0b0"
        ),
        "model_choices_sha256": (
            "31ca3300d12f664886b3a7ccc5141332eb030b204a10d6ced24defba45b09110"
        ),
        "resolution": "192x192x384",
    },
    "R07": {
        "name": "paper_standard_passive_alfvenic_beta100",
        "input": (
            "inputs/cgl_lf_paper/"
            "cgl_lf_paper_standard_passive_alfvenic_beta100.athinput"
        ),
        "input_sha256": (
            "72ed6a3b38342d00855f34a7c7eb67102b44cd22b6c12ff2448ddb6f141e5145"
        ),
        "model_choices_sha256": (
            "a23903850f2b975f309da597e2468a8258b991cad4471396679a22afacebbd9f"
        ),
        "resolution": "192x192x384",
    },
    "R14": {
        "name": "paper_nulim_beta100_20",
        "input": "inputs/cgl_lf_paper/cgl_lf_paper_nulim_beta100_20.athinput",
        "input_sha256": (
            "2b8d5837f8a7f3070ca2ef56f8b44e53d048839918185a8eef155cb084736578"
        ),
        "model_choices_sha256": (
            "a8a32525185804dbf51e2c47e95db0c6f8aba3f77541f628d87b76a387aebe96"
        ),
        "resolution": "192x192x384",
    },
    "R15": {
        "name": "paper_nulim_beta100_200",
        "input": "inputs/cgl_lf_paper/cgl_lf_paper_nulim_beta100_200.athinput",
        "input_sha256": (
            "9a698a60bef4c558ccee4635d69c3acbf3fea478401bd633d09bbc0526d943d8"
        ),
        "model_choices_sha256": (
            "32bc861a143d30fd045bed3e6553d16232324d28466a82cb000f4566cb1e4c39"
        ),
        "resolution": "192x192x384",
    },
}
FIGURE_13_FIREHOSE_CASE_NAMES = {
    str(contract["name"]): case_id
    for case_id, contract in FIGURE_13_CASE_CONTRACTS.items()
}
FIREHOSE_THRESHOLD_DEFINITIONS = {
    "parallel": {
        "beta_delta_threshold": -2.0,
        "paniso_over_b2_threshold": -1.0,
    },
    "oblique": {
        "beta_delta_threshold": -1.4,
        "paniso_over_b2_threshold": -0.7,
    },
}
MIRROR_BETA_DELTA_THRESHOLD = 1.0
INSTABILITY_OCCUPANCY_SEMANTICS = (
    "strict_manuscript",
    "inclusive_solver",
)
INSTABILITY_OCCUPANCY_COMPONENTS = (
    "mirror",
    "parallel_firehose",
    "oblique_firehose",
    "parallel_total_unstable",
    "oblique_total_unstable",
    "oblique_only_reclassification",
)
ANALYSIS_CONFIGURATION_KEYS = (
    "snapshots",
    "bundle",
    "history",
    "lf_history",
    "output_dir",
    "pdf_bins",
    "alignment_shells",
    "eddy_samples",
    "eddy_bins",
    "eddy_seed",
    "time_start",
    "time_end",
    "reference_curves",
    "stage_i_manifest",
    "allow_partial_reference_cases",
    "figure_13_only",
    "synthetic_test",
)
BIN_CONVERT_FUNCTIONS = (
    "read_binary",
    "read_binary_as_athdf",
    "read_all_ranks_binary",
    "read_all_ranks_binary_as_athdf",
)
STAGE_I_VALIDATOR_FUNCTIONS = (
    "parse_athinput",
    "require_input_contract",
    "require_product_parameter_contract",
    "restart_binary_evidence",
    "recheck_profiles",
)


def trapezoidal_integral(values: list[float], times: list[float]) -> float:
    """Integrate with the available NumPy trapezoidal-rule entry point."""

    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(values, times))
    return float(np.trapz(values, times))


def figure_13_time_close(first: float, second: float) -> bool:
    """Compare physical times at negligible scale relative to snapshot cadence."""

    return math.isclose(
        float(first),
        float(second),
        rel_tol=0.0,
        abs_tol=FIGURE_13_TIME_TOLERANCE,
    )


def validate_figure_13_snapshot_times(times: list[float]) -> dict[str, object]:
    """Require every nominal late-window cadence slot without duplicate times."""

    expected = list(FIGURE_13_SNAPSHOT_TIMES)
    if len(times) != len(expected) or not all(math.isfinite(value) for value in times):
        raise ValueError("Figure 13 snapshot cadence is incomplete or nonfinite")
    intervals = np.diff(np.asarray(times, dtype=float))
    if np.any(intervals <= FIGURE_13_TIME_TOLERANCE):
        raise ValueError("Figure 13 snapshot cadence contains duplicates or reversals")
    if any(
        not figure_13_time_close(interval, FIGURE_13_SNAPSHOT_CADENCE)
        for interval in intervals
    ) or any(
        not figure_13_time_close(actual, nominal)
        for actual, nominal in zip(times, expected)
    ):
        raise ValueError("Figure 13 snapshot cadence contains a gap or off-cadence time")
    return {
        "nominal_snapshot_times": expected,
        "time_tolerance": FIGURE_13_TIME_TOLERANCE,
        "time_tolerance_fraction_of_cadence": FIGURE_13_TIME_TOLERANCE_FRACTION,
        "maximum_nominal_time_offset": max(
            abs(actual - nominal) for actual, nominal in zip(times, expected)
        ),
    }


def parse_history(path: Path) -> dict[str, np.ndarray]:
    """Read a history file into named one-dimensional arrays."""

    labels: list[str] | None = None
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            found = HEADER_LABEL.findall(line)
            if found:
                labels = found
        elif line.strip():
            rows.append([float(value) for value in line.split()])
    if labels is None:
        raise ValueError(f"history header is missing: {path}")
    values = np.asarray(rows, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(labels):
        raise ValueError(f"history columns do not match header: {path}")
    return {label: values[:, index] for index, label in enumerate(labels)}


def time_mask(times: np.ndarray, time_start: float | None,
              time_end: float | None) -> np.ndarray:
    """Select data retained inside an optional inclusive analysis window."""

    selected = np.ones(times.shape, dtype=bool)
    if time_start is not None:
        selected &= times >= time_start
    if time_end is not None:
        selected &= times <= time_end
    return selected


def summarize_history(path: Path, time_start: float | None = None,
                      time_end: float | None = None) -> dict[str, object]:
    """Summarize reduced paper history quantities retained by the pgen."""

    data = parse_history(path)
    result: dict[str, object] = {
        "path": str(path),
        "rows": len(data["time"]),
        "time_final": float(data["time"][-1]),
        "finite": bool(all(np.isfinite(values).all() for values in data.values())),
    }
    if "volume" not in data or data["volume"][-1] == 0.0:
        return result
    volume = data["volume"][-1]
    for name in ("mass", "kinetic", "magnetic", "therm_cgl", "delta_p",
                 "abs_dp", "beta", "nu_eff", "force_pwr"):
        if name in data:
            result[f"{name}_mean_final"] = float(data[name][-1] / volume)
    for name in ("mirror_vol", "fire_vol", "hard_vol"):
        if name in data:
            result[f"{name}_fraction_final"] = float(data[name][-1] / volume)
    volumes = data["volume"]
    if np.all(volumes != 0.0):
        # Athena writes an unchanged final history row at termination; retain
        # its last value once so reference-curve coordinates remain ordered.
        series_rows = np.concatenate((np.diff(data["time"]) > 0.0, [True]))
        time_series: dict[str, object] = {
            "time": data["time"][series_rows].tolist()
        }
        fraction_names = {
            "mirror_vol": "mirror_fraction",
            "fire_vol": "firehose_fraction",
            "hard_vol": "hard_bound_fraction",
        }
        for source, target in fraction_names.items():
            if source in data:
                time_series[target] = (
                    data[source][series_rows] / volumes[series_rows]
                ).tolist()
        if "mirror_vol" in data and "fire_vol" in data:
            time_series["unstable_fraction"] = (
                (data["mirror_vol"][series_rows] + data["fire_vol"][series_rows])
                / volumes[series_rows]
            ).tolist()
        result["time_series"] = time_series
    if "force_work" in data:
        result["force_work_final"] = float(data["force_work"][-1])
    if "b2" in data and "b4" in data and data["b2"][-1] != 0.0:
        result["c_b2_final"] = float(
            data["b4"][-1] * volume / (data["b2"][-1] ** 2) - 1.0
        )
    selected = time_mask(data["time"], time_start, time_end)
    window: dict[str, object] = {
        "time_start": time_start,
        "time_end": time_end,
        "rows_selected": int(np.count_nonzero(selected)),
    }
    if np.any(selected):
        window["first_time"] = float(data["time"][selected][0])
        window["last_time"] = float(data["time"][selected][-1])
        volumes = data["volume"][selected]
        for name in ("mass", "kinetic", "magnetic", "therm_cgl", "delta_p",
                     "abs_dp", "beta", "nu_eff", "force_pwr"):
            if name in data:
                window[f"{name}_mean"] = float(np.mean(data[name][selected] / volumes))
        for name in ("mirror_vol", "fire_vol", "hard_vol"):
            if name in data:
                window[f"{name}_fraction_mean"] = float(
                    np.mean(data[name][selected] / volumes)
                )
        if "b2" in data and "b4" in data:
            b2 = data["b2"][selected]
            valid = b2 != 0.0
            if np.any(valid):
                window["c_b2_mean"] = float(np.mean(
                    data["b4"][selected][valid] * volumes[valid]
                    / (b2[valid] ** 2) - 1.0
                ))
    result["analysis_window"] = window
    return result


def summarize_forcing_energy_budget(user_path: Path, mhd_path: Path,
                                    model: dict[str, object] | None,
                                    time_start: float | None = None,
                                    time_end: float | None = None
                                    ) -> dict[str, object]:
    """Compare active CGL energy change to RK-integrated forcing work."""

    user = parse_history(user_path)
    mhd = parse_history(mhd_path)
    result: dict[str, object] = {
        "user_history": str(user_path),
        "mhd_history": str(mhd_path),
    }
    if "force_work" not in user or "tot-E" not in mhd:
        result.update({
            "available": False,
            "reason": "force_work and conserved tot-E histories are required",
        })
        return result
    if len(user["time"]) != len(mhd["time"]) or not np.allclose(
        user["time"], mhd["time"], rtol=0.0, atol=1.0e-14
    ):
        result.update({
            "available": False,
            "reason": "user and MHD history sampling times do not agree",
        })
        return result
    selected = time_mask(user["time"], time_start, time_end)
    rows = np.flatnonzero(selected)
    result["analysis_window"] = {
        "time_start": time_start,
        "time_end": time_end,
        "rows_selected": int(len(rows)),
    }
    if len(rows) < 2:
        result.update({
            "available": False,
            "reason": "at least two synchronized history rows are required",
        })
        return result
    injected = float(user["force_work"][rows[-1]] - user["force_work"][rows[0]])
    result["applied_forcing_work"] = injected
    if model is None:
        result.update({
            "available": False,
            "reason": "case model choices are required for budget interpretation",
        })
        return result
    if model_bool(model, "passive_delta", False):
        result.update({
            "available": False,
            "reason": "passive-Delta flow is not an active CGL total-energy budget",
        })
        return result
    energy_delta = float(mhd["tot-E"][rows[-1]] - mhd["tot-E"][rows[0]])
    residual = energy_delta - injected
    scale = max(abs(energy_delta), abs(injected), np.finfo(float).tiny)
    result.update({
        "available": True,
        "definition": (
            "Delta conserved tot-E - Delta RK-integrated cumulative applied "
            "force source work"
        ),
        "total_energy_delta": energy_delta,
        "residual": residual,
        "relative_residual": abs(residual) / scale,
    })
    return result


def summarize_lf_history(path: Path, time_start: float | None = None,
                         time_end: float | None = None) -> dict[str, object]:
    """Summarize cumulative LF face-cap and applied-face work diagnostics."""

    data = parse_history(path)
    result: dict[str, object] = {"path": str(path)}
    required = ("lf_qface", "lf_qprcap", "lf_qpr10", "lf_qpecap", "lf_qpe10")
    if any(name not in data for name in required):
        result["available"] = False
        return result
    selected = time_mask(data["time"], time_start, time_end)
    rows = np.flatnonzero(selected)
    result["analysis_window"] = {
        "time_start": time_start,
        "time_end": time_end,
        "rows_selected": int(len(rows)),
    }
    if len(rows) < 2:
        result["available"] = False
        result["reason"] = "at least two history rows are needed for interval counters"
        return result
    differences = {
        name: float(data[name][rows[-1]] - data[name][rows[0]])
        for name in required
    }
    qfaces = differences["lf_qface"]
    result.update({
        "available": qfaces > 0.0,
        "first_time": float(data["time"][rows[0]]),
        "last_time": float(data["time"][rows[-1]]),
        "face_evaluations": qfaces,
        "counter_differences": differences,
    })
    if qfaces > 0.0:
        result["heat_flux_cap_fractions"] = {
            "parallel_over_1": differences["lf_qprcap"] / qfaces,
            "parallel_over_10": differences["lf_qpr10"] / qfaces,
            "perpendicular_over_1": differences["lf_qpecap"] / qfaces,
            "perpendicular_over_10": differences["lf_qpe10"] / qfaces,
        }
    work_names = ("lf_qprwrk", "lf_qpewrk")
    if all(name in data for name in work_names):
        qpar_work = float(data["lf_qprwrk"][rows[-1]] - data["lf_qprwrk"][rows[0]])
        qperp_work = float(data["lf_qpewrk"][rows[-1]] - data["lf_qpewrk"][rows[0]])
        result["applied_heat_flux_work"] = {
            "definition": (
                "Signed RKL2-applied owned-face discrete contraction of "
                "-q_parallel dT_parallel and -q_perp dT_perp"
            ),
            "amr_face_ownership": (
                "coarse/fine interfaces are accumulated from the fine-side flux"
            ),
            "interpretation": (
                "operator contraction; not required to be positive or equal "
                "the snapshot reconstruction proxy"
            ),
            "signed": True,
            "parallel": qpar_work,
            "perpendicular": qperp_work,
            "total": qpar_work + qperp_work,
        }
    pressure_names = ("lf_cpwrk", "lf_cawrk")
    if all(name in data for name in pressure_names):
        pressure_work = float(data["lf_cpwrk"][rows[-1]] - data["lf_cpwrk"][rows[0]])
        anisotropic_work = float(data["lf_cawrk"][rows[-1]] - data["lf_cawrk"][rows[0]])
        result["applied_pressure_work"] = {
            "definition": (
                "Signed explicit-RK-applied cell contraction of velocity with "
                "the AMR-corrected retained CGL pressure-traction divergence"
            ),
            "traction_split": (
                "total CGL pressure traction and its Delta-p anisotropic component"
            ),
            "interpretation": "applied hyperbolic momentum-feedback ledger",
            "signed": True,
            "total": pressure_work,
            "anisotropic": anisotropic_work,
        }
    return result


def periodic_gradient(field: np.ndarray,
                      lengths: tuple[float, float, float]) -> np.ndarray:
    """Return Cartesian central differences as (x,y,z,z-index,y-index,x-index)."""

    lx, ly, lz = lengths
    nx, ny, nz = field.shape[2], field.shape[1], field.shape[0]
    derivatives = []
    for axis, spacing in ((2, lx / nx), (1, ly / ny), (0, lz / nz)):
        derivatives.append(
            (np.roll(field, -1, axis=axis) - np.roll(field, 1, axis=axis))
            / (2.0 * spacing)
        )
    return np.asarray(derivatives)


def projected_gradient(field: np.ndarray, bhat: list[np.ndarray],
                       lengths: tuple[float, float, float]
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Return local-field parallel and perpendicular gradient magnitudes."""

    gradient = periodic_gradient(field, lengths)
    parallel = sum(bhat[index] * gradient[index] for index in range(3))
    perpendicular = np.sqrt(np.maximum(
        sum(component * component for component in gradient) - parallel ** 2,
        0.0,
    ))
    return parallel, perpendicular


def model_bool(model: dict[str, object], name: str, default: bool) -> bool:
    """Convert an archived model-choice boolean with Athena input spelling."""

    value = model.get(name, default)
    if isinstance(value, bool):
        return value
    lowered = str(value).lower()
    if lowered in ("true", "1"):
        return True
    if lowered in ("false", "0"):
        return False
    raise ValueError(f"invalid boolean model choice {name}={value}")


def model_float(model: dict[str, object], name: str, default: float | None
                ) -> float | None:
    """Read one archived closure scalar, using a fallback when not archived."""

    value = model.get(name, "unspecified")
    if value in (None, "unspecified"):
        return default
    result = float(str(value))
    if not math.isfinite(result):
        raise ValueError(f"invalid finite model choice {name}={value}")
    return result


def heat_flux_transport_proxy(fields: dict[str, np.ndarray],
                              lengths: tuple[float, float, float],
                              model: dict[str, object] | None
                              ) -> dict[str, object]:
    """Reconstruct a cell-centered LF temperature-smoothing power proxy.

    This follows the implemented closure and algebraic heat-flux cap, but it
    evaluates cell-centered periodic gradients from retained snapshots. It is
    not the discrete face flux applied during time integration.
    """

    if model is None:
        return {
            "available": False,
            "reason": "archived LF model choices are required",
        }
    required_choices = (
        "lf_k_parallel",
        "lf_coefficient_mode",
        "nu_coll",
        "mirror_limiter",
        "firehose_limiter",
        "cgl_firehose_threshold",
        "limiter_nu_coll",
        "backup_limiters",
        "dfloor",
        "pfloor",
        "tfloor",
        "bfloor",
    )
    missing_choices = [name for name in required_choices if name not in model]
    if missing_choices:
        return {
            "available": False,
            "reason": (
                "archived LF model choices are incomplete: "
                + ", ".join(missing_choices)
            ),
        }
    kpar = model_float(model, "lf_k_parallel", None)
    if kpar is None or kpar <= 0.0:
        return {
            "available": False,
            "reason": "positive archived lf_k_parallel is required",
        }
    coefficient_mode = str(model.get("lf_coefficient_mode", "local"))
    if coefficient_mode not in ("local", "background"):
        raise ValueError(f"invalid lf_coefficient_mode={coefficient_mode}")
    cparallel0 = model_float(model, "lf_c_parallel0", None)
    if coefficient_mode == "background" and (
        cparallel0 is None or cparallel0 <= 0.0
    ):
        return {
            "available": False,
            "reason": "background closure requires archived lf_c_parallel0",
        }
    tiny = float(np.finfo(np.float32).tiny)
    dfloor = model_float(model, "dfloor", tiny)
    pfloor = model_float(model, "pfloor", tiny)
    tfloor = model_float(model, "tfloor", tiny)
    bfloor = model_float(model, "bfloor", math.sqrt(1024.0 * tiny))
    assert dfloor is not None and pfloor is not None
    assert tfloor is not None and bfloor is not None
    rho = np.maximum(fields["dens"], dfloor)
    ppar = np.maximum(fields["eint"], pfloor)
    pperp = np.maximum(fields["p_perp"], pfloor)
    magnetic = [fields["bcc1"], fields["bcc2"], fields["bcc3"]]
    bsqr = sum(component * component for component in magnetic)
    bmag = np.sqrt(bsqr)
    valid = bmag > bfloor
    safe_bmag = np.where(valid, bmag, 1.0)
    bhat = [
        np.where(valid, component / safe_bmag, 0.0) for component in magnetic
    ]
    tpar = ppar / rho
    tperp = pperp / rho
    grad_tpar, _ = projected_gradient(tpar, bhat, lengths)
    grad_tperp, _ = projected_gradient(tperp, bhat, lengths)
    grad_bmag, _ = projected_gradient(bmag, bhat, lengths)
    if coefficient_mode == "local":
        cparallel = np.sqrt(np.maximum(ppar / rho, tfloor))
    else:
        assert cparallel0 is not None
        cparallel = np.full(rho.shape, cparallel0)

    nu_limiter = np.zeros(rho.shape, dtype=float)
    paniso = pperp - ppar
    limiter_rate = max(model_float(model, "limiter_nu_coll", 0.0) or 0.0, 0.0)
    backup = model_bool(model, "backup_limiters", False)
    firehose_policy = str(model.get("cgl_firehose_threshold", "oblique"))
    if firehose_policy not in FIREHOSE_THRESHOLD_DEFINITIONS:
        raise ValueError(f"invalid cgl_firehose_threshold={firehose_policy}")
    firehose_threshold = float(
        FIREHOSE_THRESHOLD_DEFINITIONS[firehose_policy][
            "paniso_over_b2_threshold"
        ]
    )
    if model_bool(model, "firehose_limiter", False):
        active = paniso <= firehose_threshold * bsqr
        hard = paniso <= -1.5 * bsqr
        rate = np.where(active, limiter_rate, 0.0)
        if backup:
            rate = np.where(hard, BACKUP_COLLISION_RATE, rate)
        nu_limiter = np.maximum(nu_limiter, rate)
    if model_bool(model, "mirror_limiter", False):
        active = paniso >= 0.5 * bsqr
        hard = paniso >= bsqr
        rate = np.where(active, limiter_rate, 0.0)
        if backup:
            rate = np.where(hard, BACKUP_COLLISION_RATE, rate)
        nu_limiter = np.maximum(nu_limiter, rate)
    nu = max(model_float(model, "nu_coll", 0.0) or 0.0, 0.0) + nu_limiter
    denom_perp = SQRT_TWO_PI * cparallel * kpar + nu
    denom_parallel = SQRT_EIGHT_PI * cparallel * kpar + THREE_PI_MINUS_EIGHT * nu
    chi_perp = np.where(denom_perp > 0.0, 2.0 * cparallel ** 2 / denom_perp, 0.0)
    chi_parallel = np.where(
        denom_parallel > 0.0, 8.0 * cparallel ** 2 / denom_parallel, 0.0
    )
    qpar_unlimited = -chi_parallel * rho * grad_tpar
    qperp_unlimited = -chi_perp * (
        rho * grad_tperp
        - pperp * (1.0 - pperp / ppar) * grad_bmag / safe_bmag
    )
    qpar_max = SQRT_EIGHT_OVER_PI * cparallel * ppar
    qperp_max = SQRT_TWO_OVER_PI * cparallel * pperp
    qpar = np.where(
        qpar_max > 0.0,
        qpar_unlimited * qpar_max / (qpar_max + np.abs(qpar_unlimited)),
        0.0,
    )
    qperp = np.where(
        qperp_max > 0.0,
        qperp_unlimited * qperp_max / (qperp_max + np.abs(qperp_unlimited)),
        0.0,
    )
    qpar = np.where(valid, qpar, 0.0)
    qperp = np.where(valid, qperp, 0.0)
    qpar_unlimited = np.where(valid, qpar_unlimited, 0.0)
    qperp_unlimited = np.where(valid, qperp_unlimited, 0.0)
    qpar_ratio = np.where(
        valid & (qpar_max > 0.0), np.abs(qpar_unlimited) / qpar_max, 0.0
    )
    qperp_ratio = np.where(
        valid & (qperp_max > 0.0), np.abs(qperp_unlimited) / qperp_max, 0.0
    )
    ppar_work = -qpar * grad_tpar
    pperp_work = -qperp * grad_tperp
    ppar_unlimited_work = -qpar_unlimited * grad_tpar
    pperp_unlimited_work = -qperp_unlimited * grad_tperp
    volume = lengths[0] * lengths[1] * lengths[2]
    cell_volume = volume / rho.size

    def integral(values: np.ndarray) -> float:
        return float(cell_volume * np.sum(values))

    parallel_capped = valid & (qpar_ratio > 1.0)
    perpendicular_capped = valid & (qperp_ratio > 1.0)
    regularized_parallel = integral(ppar_work)
    regularized_perpendicular = integral(pperp_work)
    unlimited_parallel = integral(ppar_unlimited_work)
    unlimited_perpendicular = integral(pperp_unlimited_work)
    choices_used = {name: model[name] for name in required_choices}
    choices_used["limiter_hardwall"] = model.get("limiter_hardwall", "false")
    if coefficient_mode == "background":
        choices_used["lf_c_parallel0"] = model.get("lf_c_parallel0", "unspecified")
    return {
        "available": True,
        "definition": (
            "integral[-q_parallel b.grad(T_parallel) "
            "- q_perp b.grad(T_perp)] dV"
        ),
        "discretization": (
            "cell-centered periodic-gradient reconstruction; "
            "not applied finite-volume face flux"
        ),
        "closure_model_choices": choices_used,
        "valid_volume_fraction": float(np.mean(valid)),
        "regularized_parallel_power": regularized_parallel,
        "regularized_perpendicular_power": regularized_perpendicular,
        "regularized_total_power": regularized_parallel + regularized_perpendicular,
        "unlimited_parallel_power": unlimited_parallel,
        "unlimited_perpendicular_power": unlimited_perpendicular,
        "unlimited_total_power": unlimited_parallel + unlimited_perpendicular,
        "parallel_cap_active_volume_fraction": float(np.mean(parallel_capped)),
        "perpendicular_cap_active_volume_fraction": float(
            np.mean(perpendicular_capped)
        ),
        "regularized_parallel_power_on_cap_active_cells": integral(
            np.where(parallel_capped, ppar_work, 0.0)
        ),
        "regularized_perpendicular_power_on_cap_active_cells": integral(
            np.where(perpendicular_capped, pperp_work, 0.0)
        ),
    }


def velocity_gradient_products(velocity: list[np.ndarray], bhat: list[np.ndarray],
                               lengths: tuple[float, float, float]
                               ) -> dict[str, np.ndarray]:
    """Construct local-field flow-gradient products used in MKS24 diagnostics."""

    velocity_parallel = sum(
        velocity[index] * bhat[index] for index in range(3)
    )
    velocity_perp = [
        velocity[index] - bhat[index] * velocity_parallel for index in range(3)
    ]
    parallel_parallel, perpendicular_parallel = projected_gradient(
        velocity_parallel, bhat, lengths
    )
    perp_gradients = [
        projected_gradient(component, bhat, lengths) for component in velocity_perp
    ]
    parallel_perp = np.sqrt(sum(products[0] ** 2 for products in perp_gradients))
    perpendicular_perp = np.sqrt(sum(products[1] ** 2 for products in perp_gradients))
    component_gradients = [periodic_gradient(component, lengths)
                           for component in velocity]
    parallel_component = [
        sum(bhat[direction] * gradient[direction] for direction in range(3))
        for gradient in component_gradients
    ]
    strain_parallel = sum(
        bhat[index] * parallel_component[index] for index in range(3)
    )
    divergence = sum(
        component_gradients[index][index] for index in range(3)
    )
    return {
        "grad_parallel_velocity_parallel": parallel_parallel,
        "grad_perp_velocity_parallel": perpendicular_parallel,
        "grad_parallel_velocity_perp": parallel_perp,
        "grad_perp_velocity_perp": perpendicular_perp,
        "bb_grad_velocity": strain_parallel,
        "div_velocity": divergence,
        "dln_b_dt": strain_parallel - divergence,
    }


def pressure_work_decomposition(fields: dict[str, np.ndarray],
                                lengths: tuple[float, float, float],
                                strain_parallel: np.ndarray,
                                model: dict[str, object] | None,
                                divergence: np.ndarray | None = None,
                                ) -> dict[str, object]:
    """Reconstruct cell-centered CGL pressure-force work from one snapshot.

    For P = p_perp I - Delta p b b, periodic integration of the pressure
    force contribution to kinetic energy is integral[P : grad(u)] dV.
    """

    if divergence is None:
        velocity = [fields["velx"], fields["vely"], fields["velz"]]
        divergence = sum(
            periodic_gradient(velocity[index], lengths)[index] for index in range(3)
        )
    pperp = fields["p_perp"]
    delta_p = pperp - fields["eint"]
    isotropic_density = pperp * divergence
    anisotropic_density = -delta_p * strain_parallel
    volume = lengths[0] * lengths[1] * lengths[2]
    cell_volume = volume / pperp.size

    def integral(values: np.ndarray) -> float:
        return float(cell_volume * np.sum(values))

    applied_to_flow: bool | None = None
    if model is not None and "passive_delta" in model:
        applied_to_flow = not model_bool(model, "passive_delta", False)
    if applied_to_flow is True:
        interpretation = "active pressure-feedback diagnostic"
    elif applied_to_flow is False:
        interpretation = "passive-Delta diagnostic; not applied to flow evolution"
    else:
        interpretation = "model feedback scope was not archived"
    anisotropic_power = integral(anisotropic_density)
    isotropic_power = integral(isotropic_density)
    return {
        "available": True,
        "definition": (
            "integral[p_perp div(u) - Delta p (b b : grad(u))] dV"
        ),
        "anisotropic_stress_definition": (
            "integral[-Delta p (b b : grad(u))] dV"
        ),
        "discretization": (
            "cell-centered periodic-gradient reconstruction from one snapshot; "
            "not a time-integrated energy budget"
        ),
        "sign_convention": (
            "positive power is kinetic-energy gain from the CGL pressure force"
        ),
        "applied_to_flow": applied_to_flow,
        "interpretation": interpretation,
        "isotropic_perpendicular_pressure_power": isotropic_power,
        "anisotropic_stress_power": anisotropic_power,
        "total_cgl_pressure_power": isotropic_power + anisotropic_power,
        "parallel_strain_rms": float(np.sqrt(np.mean(strain_parallel ** 2))),
        "anisotropic_power_density_rms": float(
            np.sqrt(np.mean(anisotropic_density ** 2))
        ),
    }


def wavenumbers(shape: tuple[int, int, int],
                lengths: tuple[float, float, float]) -> tuple[np.ndarray, ...]:
    """Return physical FFT wavenumber grids in array indexing order."""

    nz, ny, nx = shape
    lx, ly, lz = lengths
    kx = 2.0 * math.pi * np.fft.fftfreq(nx, d=lx / nx)
    ky = 2.0 * math.pi * np.fft.fftfreq(ny, d=ly / ny)
    kz = 2.0 * math.pi * np.fft.fftfreq(nz, d=lz / nz)
    return np.meshgrid(kx, ky, kz, indexing="xy")[0].transpose(2, 0, 1), \
        np.meshgrid(kx, ky, kz, indexing="xy")[1].transpose(2, 0, 1), \
        np.meshgrid(kx, ky, kz, indexing="xy")[2].transpose(2, 0, 1)


def shell_indices(shape: tuple[int, int, int], lengths: tuple[float, float, float],
                  dk: float, perpendicular: bool) -> np.ndarray:
    """Return integer shell assignments for Fourier modes."""

    kx, ky, kz = wavenumbers(shape, lengths)
    magnitude = np.sqrt(kx * kx + ky * ky) if perpendicular else np.sqrt(
        kx * kx + ky * ky + kz * kz
    )
    return np.floor(magnitude / dk + 1.0e-12).astype(int)


def shell_spectrum(fields: list[np.ndarray], lengths: tuple[float, float, float],
                   dk: float, perpendicular: bool = True,
                   field_definition: str | None = None) -> dict[str, object]:
    """Compute an MKS24-style bin-summed spectrum for scalar or vector fields."""

    shell = shell_indices(fields[0].shape, lengths, dk, perpendicular)
    power = np.zeros(fields[0].shape, dtype=float)
    normalizer = float(fields[0].size)
    for field in fields:
        transformed = np.fft.fftn(field - np.mean(field)) / normalizer
        power += np.abs(transformed) ** 2
    binned = np.bincount(shell.ravel(), weights=power.ravel())
    result: dict[str, object] = {
        "dk": dk,
        "perpendicular": perpendicular,
        "k": (np.arange(len(binned), dtype=float) * dk).tolist(),
        "power_per_dk": (binned / dk).tolist(),
        "normalization_definition": (
            "E_chi(k) = sum_shell |FFT(chi - <chi>)/N|^2 / dk"
        ),
    }
    if field_definition is not None:
        result["field_definition"] = field_definition
    return result


def compressive_velocity_spectrum(
    velocity: list[np.ndarray], lengths: tuple[float, float, float], dk: float
) -> dict[str, object]:
    """Compute the MKS24 compressive-flow spectrum of khat dot u_k."""

    shell = shell_indices(velocity[0].shape, lengths, dk, True)
    k_components = wavenumbers(velocity[0].shape, lengths)
    magnitude = np.sqrt(sum(component * component for component in k_components))
    valid = magnitude > 0.0
    normalizer = float(velocity[0].size)
    projected = np.zeros(velocity[0].shape, dtype=complex)
    for field, component in zip(velocity, k_components):
        transformed = np.fft.fftn(field - np.mean(field)) / normalizer
        projected[valid] += transformed[valid] * component[valid] / magnitude[valid]
    binned = np.bincount(shell.ravel(), weights=np.abs(projected).ravel() ** 2)
    return {
        "dk": dk,
        "perpendicular": True,
        "k": (np.arange(len(binned), dtype=float) * dk).tolist(),
        "power_per_dk": (binned / dk).tolist(),
        "field_definition": (
            "compressive velocity khat dot u_k projected with the full "
            "three-dimensional Fourier wavevector and binned by k_perp"
        ),
        "normalization_definition": (
            "E_khat_dot_u(k_perp) = sum_shell |khat dot FFT(u)/N|^2 / dk"
        ),
    }


def pdf(values: np.ndarray, bins: int, value_range: tuple[float, float] | None = None
        ) -> dict[str, object]:
    """Return a density-normalized histogram."""

    counts, edges = np.histogram(
        values.ravel(), bins=bins, range=value_range, density=True
    )
    return {"edges": edges.tolist(), "density": counts.tolist()}


def joint_pdf(x_values: np.ndarray, y_values: np.ndarray, bins: int,
              ranges: tuple[tuple[float, float], tuple[float, float]] | None = None,
              include_counts: bool = False,
              ) -> dict[str, object]:
    """Return a density-normalized two-dimensional histogram."""

    if include_counts:
        counts, x_edges, y_edges = np.histogram2d(
            x_values.ravel(), y_values.ravel(), bins=bins, range=ranges, density=False
        )
        binned_count = float(np.sum(counts))
        if binned_count <= 0.0:
            raise ValueError("joint histogram contains no samples in the requested range")
        density = counts / binned_count
        density /= np.diff(x_edges)[:, np.newaxis]
        density /= np.diff(y_edges)[np.newaxis, :]
    else:
        density, x_edges, y_edges = np.histogram2d(
            x_values.ravel(), y_values.ravel(), bins=bins, range=ranges, density=True
        )
        counts = None
    result: dict[str, object] = {
        "x_edges": x_edges.tolist(),
        "y_edges": y_edges.tolist(),
        "density": density.tolist(),
    }
    if counts is not None:
        result.update({
            "bin_counts": counts.astype(np.int64).tolist(),
            "sample_count": int(x_values.size),
            "binned_sample_count": int(np.sum(counts)),
        })
    return result


def conditional_profile_from_joint_pdf(product: dict[str, object]) -> dict[str, object]:
    """Derive descriptive y-given-x summaries from one joint histogram."""

    x_edges = np.asarray(product["x_edges"], dtype=float)
    y_edges = np.asarray(product["y_edges"], dtype=float)
    density = np.asarray(product["density"], dtype=float)
    if (
        x_edges.ndim != 1
        or y_edges.ndim != 1
        or density.shape != (len(x_edges) - 1, len(y_edges) - 1)
        or not np.isfinite(x_edges).all()
        or not np.isfinite(y_edges).all()
        or not np.isfinite(density).all()
        or np.any(np.diff(x_edges) <= 0.0)
        or np.any(np.diff(y_edges) <= 0.0)
        or np.any(density < 0.0)
    ):
        raise ValueError("joint histogram is invalid for conditional summaries")
    x_centers = 0.5 * (x_edges[1:] + x_edges[:-1])
    y_centers = 0.5 * (y_edges[1:] + y_edges[:-1])
    x_widths = np.diff(x_edges)
    y_widths = np.diff(y_edges)
    negative_bin_fractions = np.clip(
        (np.minimum(y_edges[1:], 0.0) - y_edges[:-1]) / y_widths,
        0.0,
        1.0,
    )
    positive_bin_fractions = np.clip(
        (y_edges[1:] - np.maximum(y_edges[:-1], 0.0)) / y_widths,
        0.0,
        1.0,
    )
    y_weights = density * y_widths[np.newaxis, :]
    x_marginal_density = np.sum(y_weights, axis=1)
    x_bin_probability = x_marginal_density * x_widths
    bin_counts_value = product.get("bin_counts", product.get("bin_counts_sum"))
    bin_counts = (
        np.asarray(bin_counts_value, dtype=np.int64)
        if bin_counts_value is not None else None
    )
    if bin_counts is not None and bin_counts.shape != density.shape:
        raise ValueError("joint histogram bin counts differ from density shape")

    summaries: dict[str, list[float | int | None]] = {
        "sample_count": [],
        "response_mean": [],
        "response_rms": [],
        "response_quantile_16": [],
        "response_median": [],
        "response_quantile_84": [],
        "response_negative_fraction": [],
        "response_positive_fraction": [],
    }
    for index, marginal in enumerate(x_marginal_density):
        count = int(np.sum(bin_counts[index])) if bin_counts is not None else None
        summaries["sample_count"].append(count)
        if marginal <= 0.0:
            for name in summaries:
                if name != "sample_count":
                    summaries[name].append(None)
            continue
        weights = y_weights[index] / marginal
        cumulative = np.cumsum(weights)

        def quantile(probability: float) -> float:
            position = min(
                int(np.searchsorted(cumulative, probability, side="left")),
                len(y_centers) - 1,
            )
            return float(y_centers[position])

        summaries["response_mean"].append(float(np.sum(weights * y_centers)))
        summaries["response_rms"].append(
            float(np.sqrt(np.sum(weights * y_centers ** 2)))
        )
        summaries["response_quantile_16"].append(quantile(0.16))
        summaries["response_median"].append(quantile(0.5))
        summaries["response_quantile_84"].append(quantile(0.84))
        summaries["response_negative_fraction"].append(
            float(np.sum(weights * negative_bin_fractions))
        )
        summaries["response_positive_fraction"].append(
            float(np.sum(weights * positive_bin_fractions))
        )
    return {
        "definition": (
            "descriptive y-given-x summaries derived from the joint-PDF bins"
        ),
        "quantile_definition": (
            "first y-bin center whose conditional cumulative probability reaches "
            "the requested quantile"
        ),
        "bin_reconstruction": (
            "response moments and quantiles use y-bin centers; response-sign "
            "fractions allocate bins spanning zero in proportion to bin width"
        ),
        "scope": (
            "histogram-resolution descriptive relationship; not an independent-"
            "sample inference"
        ),
        "x_edges": x_edges.tolist(),
        "x_bin_centers": x_centers.tolist(),
        "x_marginal_density": x_marginal_density.tolist(),
        "x_bin_probability": x_bin_probability.tolist(),
        **summaries,
    }


def mechanism_joint_coordinates(
    mechanism_fields: dict[str, np.ndarray],
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Return named x/y fields for descriptive mechanism relationships."""

    return {
        name: (mechanism_fields[x_name], mechanism_fields[y_name])
        for name, (x_name, y_name) in MECHANISM_JOINT_FIELDS.items()
        if x_name in mechanism_fields and y_name in mechanism_fields
    }


def mechanism_joint_diagnostics(
    mechanism_fields: dict[str, np.ndarray],
    bins: int,
    ranges: dict[
        str, tuple[tuple[float, float], tuple[float, float]]
    ] | None = None,
    applied_to_flow: bool | None = None,
    interpretation: str = "model feedback scope was not archived",
) -> dict[str, object]:
    """Construct descriptive joint and conditional mechanism products."""

    products: dict[str, object] = {}
    for name, (x_values, y_values) in mechanism_joint_coordinates(
        mechanism_fields
    ).items():
        x_name, y_name = MECHANISM_JOINT_FIELDS[name]
        histogram = joint_pdf(
            x_values,
            y_values,
            bins,
            None if ranges is None else ranges.get(name),
            include_counts=True,
        )
        products[name] = {
            "x_field": x_name,
            "y_field": y_name,
            "x_definition": MECHANISM_FIELD_DEFINITIONS[x_name],
            "y_definition": MECHANISM_FIELD_DEFINITIONS[y_name],
            "joint_pdf": histogram,
            "conditional_y_given_x": conditional_profile_from_joint_pdf(histogram),
        }
    return {
        "schema_version": MECHANISM_JOINT_SCHEMA_VERSION,
        "scope": (
            "descriptive retained-snapshot relationships among local mechanism "
            "diagnostics"
        ),
        "discretization": (
            "cell-centered periodic-gradient reconstruction; local signed "
            "pressure-stress exchange is not applied stage accounting"
        ),
        "applied_to_flow": applied_to_flow,
        "interpretation": interpretation,
        "local_signed_exchange_field": "signed_pressure_stress_power_density",
        "scale_resolved_signed_transfer_path": "pressure_transfer.signed_transfer",
        "products": products,
    }


def filter_kperp(field: np.ndarray, shell: np.ndarray, selected: int) -> np.ndarray:
    """Filter one real field to a single perpendicular Fourier shell."""

    return np.fft.ifftn(np.fft.fftn(field) * (shell == selected)).real


def pressure_transfer(rho: np.ndarray, velocity: list[np.ndarray],
                      magnetic: list[np.ndarray], delta_p: np.ndarray,
                      lengths: tuple[float, float, float], dk: float
                      ) -> dict[str, object]:
    """Compute the MKS24 pressure-stress transfer shell partition.

    This implements integral <sqrt(rho) u>_k dot [(B/sqrt(rho)) dot grad
    ((Delta p/B^2) B)], with filtering in k_perp shells.
    MKS24 normalizes this diagnostic by the Kolmogorov estimate
    T_total ~= E_K (2 pi u_rms / L_perp).
    """

    if not np.isfinite(rho).all() or np.any(rho <= 0.0):
        raise ValueError("pressure-transfer normalization requires positive density")
    bsqr = sum(component * component for component in magnetic)
    safe_bsqr = np.maximum(bsqr, np.finfo(float).tiny)
    root_rho = np.sqrt(rho)
    weighted_velocity = [root_rho * component for component in velocity]
    stress_vector = [delta_p * component / safe_bsqr for component in magnetic]
    directional = []
    for component in stress_vector:
        gradient = periodic_gradient(component, lengths)
        directional.append(
            sum(magnetic[index] * gradient[index] for index in range(3)) / root_rho
        )
    volume = lengths[0] * lengths[1] * lengths[2]
    direct = volume * float(np.mean(sum(
        weighted_velocity[index] * directional[index] for index in range(3)
    )))
    shells = shell_indices(rho.shape, lengths, dk, True)
    values = []
    for selected in range(int(shells.max()) + 1):
        filtered = [filter_kperp(field, shells, selected) for field in weighted_velocity]
        values.append(volume * float(np.mean(sum(
            filtered[index] * directional[index] for index in range(3)
        ))))
    kinetic_energy = volume * float(np.mean(
        0.5 * rho * sum(component * component for component in velocity)
    ))
    velocity_rms = float(np.sqrt(np.mean(
        sum(component * component for component in velocity)
    )))
    lperp = math.sqrt(lengths[0] * lengths[1])
    total_transfer_rate = kinetic_energy * (2.0 * math.pi * velocity_rms / lperp)
    normalization_available = bool(
        math.isfinite(total_transfer_rate)
        and total_transfer_rate > np.finfo(float).tiny
    )
    return {
        "dk": dk,
        "k_perp": (np.arange(len(values), dtype=float) * dk).tolist(),
        "transfer": values,
        "signed_transfer": list(values),
        "definition": (
            "perpendicular-shell partition of integral[sqrt(rho) u . "
            "((B/sqrt(rho)) . grad((Delta p/B^2) B))] dV"
        ),
        "filter_definition": (
            "sqrt(rho) u is filtered into nonoverlapping k_perp shells; "
            "the pressure-stress force factor is unfiltered"
        ),
        "sign_convention": (
            "positive signed transfer corresponds to kinetic-energy gain if the "
            "anisotropic CGL pressure-stress force is applied"
        ),
        "interpretation": (
            "signed pressure-stress exchange reconstructed from one retained "
            "snapshot; model feedback scope is not specified here"
        ),
        "normalization_available": normalization_available,
        "normalization_definition": (
            "T_total ~= E_K (2 pi u_rms / L_perp), with "
            "E_K = integral[0.5 rho |u|^2] dV, "
            "u_rms = sqrt(<|u|^2>), and L_perp = sqrt(Lx Ly)"
        ),
        "kinetic_energy": kinetic_energy,
        "velocity_rms": velocity_rms,
        "perpendicular_outer_scale": lperp,
        "total_transfer_rate": total_transfer_rate,
        "transfer_normalized_by_total": (
            [value / total_transfer_rate for value in values]
            if normalization_available else None
        ),
        "signed_transfer_normalized_by_total": (
            [value / total_transfer_rate for value in values]
            if normalization_available else None
        ),
        "direct_real_space": direct,
        "shell_sum": float(sum(values)),
        "closure_error": float(sum(values) - direct),
    }


def alignment_histograms(velocity: list[np.ndarray], magnetic: list[np.ndarray],
                         lengths: tuple[float, float, float], dk: float,
                         selected_shells: list[int], bins: int) -> dict[str, object]:
    """Compute stretching-eigenvector alignment PDFs for selected k_perp shells."""

    bsqr = sum(component * component for component in magnetic)
    bhat = [component / np.sqrt(np.maximum(bsqr, np.finfo(float).tiny))
            for component in magnetic]
    shells = shell_indices(velocity[0].shape, lengths, dk, True)
    output: dict[str, object] = {}
    for selected in selected_shells:
        if not np.any(shells == selected):
            continue
        filtered = [filter_kperp(field, shells, selected) for field in velocity]
        gradients = [periodic_gradient(field, lengths) for field in filtered]
        strain = np.empty(velocity[0].shape + (3, 3), dtype=float)
        for row in range(3):
            for column in range(3):
                strain[..., row, column] = 0.5 * (
                    gradients[row][column] + gradients[column][row]
                )
        _, vectors = np.linalg.eigh(strain)
        stretching = vectors[..., :, 2]
        cosine = np.abs(sum(stretching[..., index] * bhat[index] for index in range(3)))
        output[str(selected)] = pdf(cosine, bins, (0.0, 1.0))
    return output


def eddy_anisotropy_curve(bin_centers: np.ndarray, perpendicular: np.ndarray,
                          parallel: np.ndarray) -> dict[str, object]:
    """Invert perpendicular and parallel structure functions at common power."""

    parallel_valid = np.isfinite(parallel) & (parallel > 0.0)
    perpendicular_valid = np.isfinite(perpendicular) & (perpendicular > 0.0)
    parallel_length = bin_centers[parallel_valid]
    parallel_power = parallel[parallel_valid]
    monotonic: list[int] = []
    maximum = -math.inf
    for index, value in enumerate(parallel_power):
        if value > maximum:
            monotonic.append(index)
            maximum = float(value)
    if len(monotonic) < 2:
        return {
            "available": False,
            "reason": "parallel structure function has fewer than two increasing bins",
        }
    parallel_length = parallel_length[monotonic]
    parallel_power = parallel_power[monotonic]
    selected = perpendicular_valid & (perpendicular >= parallel_power[0]) & (
        perpendicular <= parallel_power[-1]
    )
    if np.count_nonzero(selected) < 2:
        return {
            "available": False,
            "reason": "structure functions do not overlap on at least two bins",
        }
    ell_perp = bin_centers[selected]
    ell_parallel = np.exp(np.interp(
        np.log(perpendicular[selected]),
        np.log(parallel_power),
        np.log(parallel_length),
    ))
    return {
        "available": True,
        "ell_perp_over_lperp": ell_perp.tolist(),
        "ell_parallel_over_lperp": ell_parallel.tolist(),
    }


def local_field_eddy_anisotropy(
    velocity: list[np.ndarray], magnetic: list[np.ndarray],
    lengths: tuple[float, float, float], samples: int, bins: int, seed: int
) -> dict[str, object]:
    """Estimate local-field-conditioned three-point eddy anisotropy curves."""

    definition = (
        "solve S2(phi; ell_perp) = S2(phi; ell_parallel), where "
        "S2 = <|phi(x+ell) - 2 phi(x) + phi(x-ell)|^2>"
    )
    if samples <= 0:
        return {
            "computed": False,
            "available": False,
            "definition": definition,
            "reason": "eddy anisotropy analysis is opt-in; set --eddy-samples",
        }
    nz, ny, nx = velocity[0].shape
    lx, ly, lz = lengths
    lperp = math.sqrt(lx * ly)
    spacing_xyz = np.asarray((lx / nx, ly / ny, lz / nz), dtype=float)
    minimum = max(float(np.min(spacing_xyz)), np.finfo(float).tiny)
    maximum = 0.5 * min(lx, ly)
    if maximum <= minimum:
        return {
            "computed": False,
            "available": False,
            "definition": definition,
            "reason": "snapshot has insufficient perpendicular scale separation",
        }
    edges = np.geomspace(minimum, maximum, bins + 1)
    centers = np.sqrt(edges[:-1] * edges[1:]) / lperp
    per_bin = max(1, int(math.ceil(samples / bins)))
    generated = per_bin * bins
    generator = np.random.default_rng(seed)
    source_bins = np.repeat(np.arange(bins), per_bin)
    radius = np.exp(generator.uniform(
        np.log(edges[source_bins]), np.log(edges[source_bins + 1])
    ))
    directions = generator.normal(size=(generated, 3))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    offsets_xyz = np.rint(
        directions * radius[:, None] / spacing_xyz[None, :]
    ).astype(int)
    separation_xyz = offsets_xyz * spacing_xyz[None, :]
    separation = np.linalg.norm(separation_xyz, axis=1)
    retained = (separation >= edges[0]) & (separation <= edges[-1])
    retained &= np.any(offsets_xyz != 0, axis=1)
    offsets_xyz = offsets_xyz[retained]
    separation_xyz = separation_xyz[retained]
    separation = separation[retained]
    sample_bins = np.searchsorted(edges, separation, side="right") - 1
    valid_bins = (sample_bins >= 0) & (sample_bins < bins)
    offsets_xyz = offsets_xyz[valid_bins]
    separation_xyz = separation_xyz[valid_bins]
    separation = separation[valid_bins]
    sample_bins = sample_bins[valid_bins]
    count = len(sample_bins)
    center_z = generator.integers(0, nz, size=count)
    center_y = generator.integers(0, ny, size=count)
    center_x = generator.integers(0, nx, size=count)
    offset_x = offsets_xyz[:, 0]
    offset_y = offsets_xyz[:, 1]
    offset_z = offsets_xyz[:, 2]
    plus = (
        (center_z + offset_z) % nz,
        (center_y + offset_y) % ny,
        (center_x + offset_x) % nx,
    )
    center = (center_z, center_y, center_x)
    minus = (
        (center_z - offset_z) % nz,
        (center_y - offset_y) % ny,
        (center_x - offset_x) % nx,
    )
    local_field = [
        (component[plus] + component[center] + component[minus]) / 3.0
        for component in magnetic
    ]
    field_norm = np.sqrt(sum(component * component for component in local_field))
    valid_field = field_norm > np.finfo(float).tiny
    bhat = [
        component / np.maximum(field_norm, np.finfo(float).tiny)
        for component in local_field
    ]
    separation_hat = separation_xyz / separation[:, None]
    cosine = np.abs(sum(
        separation_hat[:, index] * bhat[index] for index in range(3)
    ))
    radians = math.radians(EDDY_ANGLE_DEGREES)
    parallel_selected = valid_field & (cosine >= math.cos(radians))
    perpendicular_selected = valid_field & (cosine <= math.sin(radians))

    def second_order_perp(vector: list[np.ndarray]) -> np.ndarray:
        sampled = [
            [component[location] for component in vector]
            for location in (plus, center, minus)
        ]
        perpendicular_values = []
        for values in sampled:
            field_parallel = sum(values[index] * bhat[index] for index in range(3))
            perpendicular_values.append([
                values[index] - field_parallel * bhat[index] for index in range(3)
            ])
        return sum(
            (
                perpendicular_values[0][index]
                - 2.0 * perpendicular_values[1][index]
                + perpendicular_values[2][index]
            ) ** 2
            for index in range(3)
        )

    def conditioned_mean(values: np.ndarray, selected: np.ndarray
                         ) -> tuple[np.ndarray, np.ndarray]:
        counts = np.bincount(sample_bins[selected], minlength=bins)
        sums = np.bincount(
            sample_bins[selected], weights=values[selected], minlength=bins
        )
        mean = np.zeros(bins, dtype=float)
        populated = counts > 0
        mean[populated] = sums[populated] / counts[populated]
        return mean, counts

    output: dict[str, object] = {
        "computed": True,
        "available": True,
        "definition": definition,
        "conditioning": (
            "three-point local mean magnetic field; separation vectors within "
            f"{EDDY_ANGLE_DEGREES:g} degrees of parallel or perpendicular"
        ),
        "sampling": (
            "deterministic random lattice separations, logarithmically balanced "
            "over normalized separation bins"
        ),
        "separation_coordinate": (
            "|ell|/L_perp binned within each angular cone; the selected "
            "parallel or perpendicular projection differs by at most "
            "1 - cos(15 degrees)"
        ),
        "normalization_definition": "lengths divided by L_perp = sqrt(Lx Ly)",
        "lperp": lperp,
        "angle_degrees": EDDY_ANGLE_DEGREES,
        "samples_requested": samples,
        "samples_retained": int(count),
        "seed": seed,
        "bins": bins,
        "bin_centers_over_lperp": centers.tolist(),
    }
    for name, vector in (
        ("velocity_perp", velocity),
        ("magnetic_perp", magnetic),
    ):
        values = second_order_perp(vector)
        perpendicular, perpendicular_counts = conditioned_mean(
            values, perpendicular_selected
        )
        parallel, parallel_counts = conditioned_mean(values, parallel_selected)
        product = eddy_anisotropy_curve(centers, perpendicular, parallel)
        product.update({
            "perpendicular_structure_function": perpendicular.tolist(),
            "parallel_structure_function": parallel.tolist(),
            "perpendicular_sample_counts": perpendicular_counts.tolist(),
            "parallel_sample_counts": parallel_counts.tolist(),
        })
        output[name] = product
    output["available"] = bool(
        output["velocity_perp"]["available"]
        or output["magnetic_perp"]["available"]
    )
    if not output["available"]:
        output["reason"] = (
            "no eddy-anisotropy product has overlapping structure functions"
        )
    return output


def binary_metadata_equal(first: object, second: object) -> bool:
    """Return whether two rank-file metadata values are exactly compatible."""

    if isinstance(first, np.ndarray) or isinstance(second, np.ndarray):
        return bool(np.array_equal(np.asarray(first), np.asarray(second)))
    return first == second


def read_exact_rank_set_binary(rank_files: list[Path]) -> dict[str, object]:
    """Read and combine exactly one authenticated contiguous rank-file set."""

    if not rank_files:
        raise ValueError("rank-local snapshot requires at least one rank file")
    payloads = [bin_convert.read_binary(str(path)) for path in rank_files]
    reference = payloads[0]
    metadata_keys = (
        "header",
        "time",
        "cycle",
        "var_names",
        "Nx1",
        "Nx2",
        "Nx3",
        "nvars",
        "x1min",
        "x1max",
        "x2min",
        "x2max",
        "x3min",
        "x3max",
        "nx1_mb",
        "nx2_mb",
        "nx3_mb",
        "nx1_out_mb",
        "nx2_out_mb",
        "nx3_out_mb",
    )
    for rank, payload in enumerate(payloads[1:], start=1):
        for key in metadata_keys:
            if key not in reference or key not in payload or not binary_metadata_equal(
                reference[key], payload[key]
            ):
                raise ValueError(
                    "rank-local snapshot metadata differs across exact rank set: "
                    f"rank={rank}, key={key}"
                )
    combined = reference.copy()
    for key in ("mb_index", "mb_logical", "mb_geometry"):
        combined[key] = np.concatenate(
            [np.asarray(payload[key]) for payload in payloads], axis=0
        )
    combined["mb_data"] = {
        name: np.concatenate(
            [np.asarray(payload["mb_data"][name]) for payload in payloads],
            axis=0,
        )
        for name in reference["var_names"]
    }
    combined["n_mbs"] = len(combined["mb_index"])
    validate_exact_rank_set_meshblocks(combined)
    return combined


def validate_exact_rank_set_meshblocks(raw: dict[str, object]) -> None:
    """Reject duplicate or incomplete meshblock inventories before reconstruction."""

    count = int(raw["n_mbs"])
    logical = np.asarray(raw["mb_logical"])
    indices = np.asarray(raw["mb_index"])
    geometry = np.asarray(raw["mb_geometry"])
    if (
        count <= 0
        or logical.shape != (count, 4)
        or indices.shape != (count, 6)
        or geometry.shape != (count, 6)
        or not np.issubdtype(logical.dtype, np.integer)
        or not np.issubdtype(indices.dtype, np.integer)
        or np.any(logical < 0)
        or np.any(indices < 0)
        or np.any(indices[:, 0::2] > indices[:, 1::2])
    ):
        raise ValueError("rank-local snapshot meshblock inventory is malformed")
    logical_rows = [tuple(int(value) for value in row) for row in logical]
    if len(set(logical_rows)) != count:
        raise ValueError("rank-local snapshot contains duplicate logical meshblocks")
    if not np.isfinite(geometry).all() or any(
        np.any(geometry[:, lower] >= geometry[:, lower + 1])
        for lower in (0, 2, 4)
    ):
        raise ValueError("rank-local snapshot meshblock geometry is invalid")
    for name in raw["var_names"]:
        values = np.asarray(raw["mb_data"][name])
        if values.shape[0] != count:
            raise ValueError(
                f"rank-local snapshot meshblock data count differs: {name}"
            )

    # Figure 13 uses full, uniform, level-zero output. For such output the exact
    # Cartesian logical inventory is unambiguous, so do not permit holes.
    root_sizes = tuple(int(raw[f"Nx{axis}"]) for axis in (1, 2, 3))
    block_sizes = tuple(int(raw[f"nx{axis}_mb"]) for axis in (1, 2, 3))
    output_block_sizes = tuple(
        int(raw[f"nx{axis}_out_mb"]) for axis in (1, 2, 3)
    )
    if any(value <= 0 for value in (*root_sizes, *block_sizes, *output_block_sizes)):
        raise ValueError("rank-local snapshot meshblock dimensions are invalid")
    full_uniform = bool(
        np.all(logical[:, 3] == 0)
        and output_block_sizes == block_sizes
        and all(root % block == 0 for root, block in zip(root_sizes, block_sizes))
    )
    if full_uniform:
        block_counts = tuple(
            root // block for root, block in zip(root_sizes, block_sizes)
        )
        expected = {
            (i, j, k, 0)
            for k in range(block_counts[2])
            for j in range(block_counts[1])
            for i in range(block_counts[0])
        }
        if set(logical_rows) != expected:
            raise ValueError(
                "rank-local snapshot uniform meshblock coverage is incomplete"
            )


def exact_rank_set_as_athdf(
    path: Path, raw: dict[str, object],
) -> dict[str, np.ndarray]:
    """Convert one exact combined rank set without permitting wildcard discovery."""

    original = bin_convert.read_all_ranks_binary_as_athdf
    if not isinstance(original, types.FunctionType):
        raise ValueError("exact rank-set converter is not a Python function")
    isolated_globals = original.__globals__.copy()
    isolated_globals["read_all_ranks_binary"] = lambda _filename: raw
    converter = types.FunctionType(
        original.__code__,
        isolated_globals,
        name=original.__name__,
        argdefs=original.__defaults__,
        closure=original.__closure__,
    )
    converter.__kwdefaults__ = original.__kwdefaults__
    return converter(
        str(path), quantities=list(REQUIRED_FIELDS), dtype=np.float64
    )


def read_snapshot(
    path: Path, rank_files: list[Path] | None = None,
) -> tuple[dict[str, np.ndarray], tuple[float, float, float], float]:
    """Read one current CGL primitive binary snapshot from an exact file set."""

    rank_local = path.parent.name == "rank_00000000"
    if rank_local:
        if rank_files is None:
            rank_files = snapshot_sibling_paths(path)
        if (
            not rank_files
            or rank_files[0] != path
            or rank_files != snapshot_sibling_paths(path, len(rank_files))
        ):
            raise ValueError(
                f"rank-local snapshot exact rank set is invalid: {path}"
            )
        raw = read_exact_rank_set_binary(rank_files)
    else:
        if rank_files not in (None, [path]):
            raise ValueError(f"shared snapshot received an invalid exact file set: {path}")
        raw = bin_convert.read_binary(str(path))
    missing = sorted(set(REQUIRED_FIELDS) - set(raw["var_names"]))
    if missing:
        raise ValueError(
            f"{path} lacks CGL paper fields {missing}; regenerate with current "
            "mhd_w_bcc output (legacy eint is p_parallel)"
        )
    values = (
        exact_rank_set_as_athdf(path, raw)
        if rank_local else
        bin_convert.read_binary_as_athdf(
            str(path), quantities=list(REQUIRED_FIELDS), dtype=np.float64
        )
    )
    lengths = (
        float(raw["x1max"] - raw["x1min"]),
        float(raw["x2max"] - raw["x2min"]),
        float(raw["x3max"] - raw["x3min"]),
    )
    return values, lengths, float(raw["time"])


def snapshot_time(path: Path) -> float:
    """Read the time in one Athena binary preheader without loading field data."""

    with path.open("rb") as stream:
        code_header = stream.readline().split()
        if not code_header or code_header[0] != b"Athena":
            raise ValueError(f"binary snapshot has invalid header: {path}")
        pheader_count = int(stream.readline().split(b"=")[-1])
        values: dict[str, str] = {}
        for _ in range(pheader_count - 1):
            key, value = [
                token.strip()
                for token in stream.readline().decode("utf-8").split("=", 1)
            ]
            values[key] = value
    if "time" not in values:
        raise ValueError(f"binary snapshot has no time field: {path}")
    return float(values["time"])


def snapshot_sibling_paths(
    path: Path, expected_ranks: int | None = None,
) -> list[Path]:
    """Return the exact contiguous file set represented by one snapshot path."""

    if expected_ranks is not None and (
        type(expected_ranks) is not int or expected_ranks <= 0
    ):
        raise ValueError(f"snapshot expected rank count is invalid: {expected_ranks}")
    if path.parent.name != "rank_00000000":
        if expected_ranks not in (None, 1):
            raise ValueError(
                f"shared snapshot cannot satisfy expected rank count {expected_ranks}: "
                f"{path}"
            )
        siblings = [path]
    else:
        rank_root = path.parent.parent
        rank_entries = sorted(
            candidate for candidate in rank_root.iterdir()
            if candidate.name.startswith("rank_")
        )
        invalid = [
            candidate for candidate in rank_entries
            if not candidate.is_dir()
            or re.fullmatch(r"rank_\d{8}", candidate.name) is None
        ]
        if invalid:
            raise ValueError(
                "rank-local snapshot inventory contains injected rank entries: "
                + ", ".join(str(candidate) for candidate in invalid)
        )
        if expected_ranks is None:
            expected_ranks = len(rank_entries)
        expected_names = [
            f"rank_{rank:08d}" for rank in range(expected_ranks)
        ]
        actual_names = [candidate.name for candidate in rank_entries]
        if actual_names != expected_names:
            raise ValueError(
                f"rank-local snapshot set is not the exact expected contiguous set "
                f"{expected_names}: {path}"
            )
        siblings = [rank_root / name / path.name for name in expected_names]
    if not siblings or path not in siblings:
        raise ValueError(f"snapshot provenance lacks representative file: {path}")
    for sibling in siblings:
        if not sibling.is_file():
            raise ValueError(f"snapshot provenance lacks rank sibling: {sibling}")
    return siblings


def snapshot_digest_provenance(
    path: Path, expected_ranks: int | None = None,
    siblings: list[Path] | None = None,
) -> dict[str, object]:
    """Bind one logical snapshot to all rank-local sibling bytes."""

    exact_siblings = snapshot_sibling_paths(path, expected_ranks)
    if siblings is not None and siblings != exact_siblings:
        raise ValueError(f"snapshot exact rank set changed before hashing: {path}")
    files = []
    retained_identities = []
    for sibling in exact_siblings:
        symlink_target = os.readlink(sibling) if sibling.is_symlink() else None
        before = sibling.stat()
        digest = sha256_file(sibling)
        after = sibling.stat()
        retained_stat = (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
        )
        if retained_stat != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
        ) or symlink_target != (
            os.readlink(sibling) if sibling.is_symlink() else None
        ):
            raise ValueError(f"snapshot changed while hashing provenance: {sibling}")
        retained_identities.append((retained_stat, symlink_target))
        files.append({
            "path": str(sibling),
            "symlink_target": symlink_target,
            "size_bytes": before.st_size,
            "sha256": digest,
        })
    time = snapshot_time(path)
    if not math.isfinite(time):
        raise ValueError(f"snapshot provenance time is nonfinite: {path}")
    if snapshot_sibling_paths(path, expected_ranks) != exact_siblings:
        raise ValueError(f"snapshot exact rank set changed while hashing: {path}")
    for sibling, (retained_stat, symlink_target) in zip(
        exact_siblings, retained_identities
    ):
        final = sibling.stat()
        if retained_stat != (
            final.st_dev,
            final.st_ino,
            final.st_size,
            final.st_mtime_ns,
        ) or symlink_target != (
            os.readlink(sibling) if sibling.is_symlink() else None
        ):
            raise ValueError(f"snapshot changed while hashing provenance: {sibling}")
    return {
        "representative_path": str(path),
        "layout": (
            "rank_local_siblings"
            if path.parent.name == "rank_00000000"
            else "single_file"
        ),
        "expected_rank_count": len(exact_siblings),
        "rank_directory_names": [
            sibling.parent.name for sibling in exact_siblings
            if path.parent.name == "rank_00000000"
        ],
        "files": files,
        "aggregate_sha256": canonical_json_sha256(files),
        "snapshot_time": time,
    }


def firehose_threshold_definitions() -> dict[str, dict[str, object]]:
    """Return inclusive solver firehose thresholds in beta-Delta conventions."""

    return {
        name: {
            "beta_delta_threshold": float(values["beta_delta_threshold"]),
            "paniso_over_b2_threshold": float(
                values["paniso_over_b2_threshold"]
            ),
            "comparison_operator": "<=",
            "definition": (
                "beta_delta = 2 (p_perp - p_parallel) / B^2 "
                f"<= {values['beta_delta_threshold']}"
            ),
            "cgl_firehose_threshold_policy": name,
        }
        for name, values in FIREHOSE_THRESHOLD_DEFINITIONS.items()
    }


def instability_threshold_definitions() -> dict[str, object]:
    """Return strict manuscript and inclusive solver occupancy conventions."""

    return {
        "beta_delta_definition": (
            "beta_delta = 2 (p_perp - p_parallel) / B^2"
        ),
        "semantics": {
            "strict_manuscript": {
                "definition": (
                    "strict instability regions stated in the MKS24 manuscript"
                ),
                "mirror": {
                    "beta_delta_threshold": MIRROR_BETA_DELTA_THRESHOLD,
                    "comparison_operator": ">",
                },
                "parallel_firehose": {
                    "beta_delta_threshold": float(
                        FIREHOSE_THRESHOLD_DEFINITIONS["parallel"][
                            "beta_delta_threshold"
                        ]
                    ),
                    "comparison_operator": "<",
                },
                "oblique_firehose": {
                    "beta_delta_threshold": float(
                        FIREHOSE_THRESHOLD_DEFINITIONS["oblique"][
                            "beta_delta_threshold"
                        ]
                    ),
                    "comparison_operator": "<",
                },
            },
            "inclusive_solver": {
                "definition": (
                    "inclusive threshold occupancy matching limiter activation"
                ),
                "mirror": {
                    "beta_delta_threshold": MIRROR_BETA_DELTA_THRESHOLD,
                    "comparison_operator": ">=",
                },
                "parallel_firehose": {
                    "beta_delta_threshold": float(
                        FIREHOSE_THRESHOLD_DEFINITIONS["parallel"][
                            "beta_delta_threshold"
                        ]
                    ),
                    "comparison_operator": "<=",
                },
                "oblique_firehose": {
                    "beta_delta_threshold": float(
                        FIREHOSE_THRESHOLD_DEFINITIONS["oblique"][
                            "beta_delta_threshold"
                        ]
                    ),
                    "comparison_operator": "<=",
                },
            },
        },
    }


def beta_delta_field(fields: dict[str, np.ndarray]) -> np.ndarray:
    """Return beta-Delta after validating the fields needed by its definition."""

    required = ("eint", "p_perp", "bcc1", "bcc2", "bcc3")
    missing = [name for name in required if name not in fields]
    if missing:
        raise ValueError(
            "firehose occupancy requires retained snapshot fields: "
            + ", ".join(missing)
        )
    values = {name: np.asarray(fields[name], dtype=float) for name in required}
    shape = values[required[0]].shape
    if not shape or any(values[name].shape != shape for name in required[1:]):
        raise ValueError("firehose occupancy fields have inconsistent grids")
    for name, field in values.items():
        if not np.isfinite(field).all():
            raise ValueError(f"firehose occupancy field is nonfinite: {name}")
    if np.any(values["eint"] <= 0.0) or np.any(values["p_perp"] <= 0.0):
        raise ValueError(
            "firehose occupancy requires positive p_parallel and p_perp"
        )
    with np.errstate(over="ignore", invalid="ignore"):
        bsqr = (
            values["bcc1"] ** 2
            + values["bcc2"] ** 2
            + values["bcc3"] ** 2
        )
    if not np.isfinite(bsqr).all() or np.any(bsqr <= 0.0):
        raise ValueError(
            "firehose occupancy beta_delta requires finite positive magnetic "
            "field strength"
        )
    beta_delta = 2.0 * (values["p_perp"] - values["eint"]) / bsqr
    if not np.isfinite(beta_delta).all():
        raise ValueError("firehose occupancy beta_delta is nonfinite")
    return beta_delta


def firehose_threshold_occupancy_from_beta_delta(
    beta_delta: np.ndarray,
) -> dict[str, object]:
    """Measure mirror-or-firehose occupancy under both threshold semantics."""

    beta_delta = np.asarray(beta_delta, dtype=float)
    if not beta_delta.shape or not np.isfinite(beta_delta).all():
        raise ValueError("firehose occupancy requires finite beta_delta cells")
    total = int(beta_delta.size)
    if total <= 0:
        raise ValueError("firehose occupancy requires at least one snapshot cell")
    parallel_threshold = float(
        FIREHOSE_THRESHOLD_DEFINITIONS["parallel"]["beta_delta_threshold"]
    )
    oblique_threshold = float(
        FIREHOSE_THRESHOLD_DEFINITIONS["oblique"]["beta_delta_threshold"]
    )
    masks_by_semantics: dict[str, dict[str, np.ndarray]] = {}
    for semantics in INSTABILITY_OCCUPANCY_SEMANTICS:
        if semantics == "strict_manuscript":
            mirror = beta_delta > MIRROR_BETA_DELTA_THRESHOLD
            parallel_firehose = beta_delta < parallel_threshold
            oblique_firehose = beta_delta < oblique_threshold
        else:
            mirror = beta_delta >= MIRROR_BETA_DELTA_THRESHOLD
            parallel_firehose = beta_delta <= parallel_threshold
            oblique_firehose = beta_delta <= oblique_threshold
        if np.any(parallel_firehose & ~oblique_firehose):
            raise ValueError("parallel firehose mask is not contained in oblique mask")
        if np.any(mirror & oblique_firehose):
            raise ValueError("mirror and firehose masks overlap")
        masks_by_semantics[semantics] = {
            "mirror": mirror,
            "parallel_firehose": parallel_firehose,
            "oblique_firehose": oblique_firehose,
            "parallel_total_unstable": mirror | parallel_firehose,
            "oblique_total_unstable": mirror | oblique_firehose,
            "oblique_only_reclassification": (
                oblique_firehose & ~parallel_firehose
            ),
        }
    return {
        "schema_version": FIGURE_13_FIREHOSE_SCHEMA_VERSION,
        "definition": (
            "instantaneous retained-snapshot mirror-or-firehose occupancy under "
            "parallel and alternative oblique firehose definitions"
        ),
        "normalization": {
            "kind": "equal_cell_volume_fraction",
            "denominator": "all cells in the reconstructed uniform snapshot grid",
            "total_cell_count": total,
        },
        "threshold_definitions": instability_threshold_definitions(),
        "beta_delta_minimum": float(np.min(beta_delta)),
        "beta_delta_maximum": float(np.max(beta_delta)),
        "occupancy": {
            semantics: {
                name: {
                    "cell_count": int(np.count_nonzero(mask)),
                    "volume_fraction": float(np.count_nonzero(mask) / total),
                }
                for name, mask in masks.items()
            }
            for semantics, masks in masks_by_semantics.items()
        },
    }


def firehose_threshold_occupancy(
    fields: dict[str, np.ndarray],
) -> dict[str, object]:
    """Measure microinstability occupancy from retained snapshot fields."""

    return firehose_threshold_occupancy_from_beta_delta(beta_delta_field(fields))


def validate_instantaneous_firehose_occupancy(
    product: object, expected_total_cells: int, context: str,
) -> dict[str, object]:
    """Validate one instantaneous occupancy product against its exact counts."""

    if not isinstance(product, dict):
        raise ValueError(f"{context} occupancy product is malformed")
    if (
        product.get("schema_version") != FIGURE_13_FIREHOSE_SCHEMA_VERSION
        or product.get("threshold_definitions") != instability_threshold_definitions()
    ):
        raise ValueError(f"{context} occupancy threshold definitions do not match")
    normalization = product.get("normalization")
    if (
        not isinstance(normalization, dict)
        or normalization.get("kind") != "equal_cell_volume_fraction"
        or normalization.get("total_cell_count") != expected_total_cells
    ):
        raise ValueError(f"{context} occupancy normalization does not match its grid")
    try:
        beta_delta_minimum = float(product["beta_delta_minimum"])
        beta_delta_maximum = float(product["beta_delta_maximum"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{context} occupancy beta_delta bounds are malformed") from error
    if (
        not math.isfinite(beta_delta_minimum)
        or not math.isfinite(beta_delta_maximum)
        or beta_delta_minimum > beta_delta_maximum
    ):
        raise ValueError(f"{context} occupancy beta_delta bounds are invalid")

    occupancy = product.get("occupancy")
    if not isinstance(occupancy, dict) or set(occupancy) != set(
        INSTABILITY_OCCUPANCY_SEMANTICS
    ):
        raise ValueError(f"{context} occupancy semantics are malformed")
    for semantics in INSTABILITY_OCCUPANCY_SEMANTICS:
        measured = occupancy[semantics]
        if not isinstance(measured, dict) or set(measured) != set(
            INSTABILITY_OCCUPANCY_COMPONENTS
        ):
            raise ValueError(f"{context} occupancy components are malformed: {semantics}")
        counts: dict[str, int] = {}
        for name in INSTABILITY_OCCUPANCY_COMPONENTS:
            record = measured[name]
            if not isinstance(record, dict):
                raise ValueError(
                    f"{context} occupancy component is malformed: {semantics}.{name}"
                )
            count = record.get("cell_count")
            fraction = record.get("volume_fraction")
            if (
                isinstance(count, bool)
                or not isinstance(count, int)
                or count < 0
                or count > expected_total_cells
                or not isinstance(fraction, (int, float))
                or not math.isfinite(float(fraction))
                or not math.isclose(
                    float(fraction),
                    count / expected_total_cells,
                    rel_tol=0.0,
                    abs_tol=1.0e-15,
                )
            ):
                raise ValueError(
                    f"{context} occupancy count/fraction differs: {semantics}.{name}"
                )
            counts[name] = count
        if (
            counts["oblique_firehose"] - counts["parallel_firehose"]
            != counts["oblique_only_reclassification"]
            or counts["parallel_total_unstable"]
            != counts["mirror"] + counts["parallel_firehose"]
            or counts["oblique_total_unstable"]
            != counts["mirror"] + counts["oblique_firehose"]
        ):
            raise ValueError(f"{context} occupancy mask counts are inconsistent")
    return product


def pdf_fields(fields: dict[str, np.ndarray],
               lengths: tuple[float, float, float] | None = None,
               mechanism_fields: dict[str, np.ndarray] | None = None,
               ) -> dict[str, np.ndarray]:
    """Construct the scalar fields used for paper PDF products."""

    rho = fields["dens"]
    ppar = fields["eint"]
    pperp = fields["p_perp"]
    bsqr = fields["bcc1"] ** 2 + fields["bcc2"] ** 2 + fields["bcc3"] ** 2
    values = {
        "density_fluctuation": rho / np.mean(rho) - 1.0,
        "p_parallel_fluctuation": ppar / np.mean(ppar) - 1.0,
        "p_perp_fluctuation": pperp / np.mean(pperp) - 1.0,
        "beta_delta": beta_delta_field(fields),
    }
    if lengths is not None:
        if mechanism_fields is None:
            magnetic = [fields["bcc1"], fields["bcc2"], fields["bcc3"]]
            velocity = [fields["velx"], fields["vely"], fields["velz"]]
            bhat = [
                component / np.sqrt(np.maximum(bsqr, np.finfo(float).tiny))
                for component in magnetic
            ]
            mechanism_fields = velocity_gradient_products(velocity, bhat, lengths)
            mechanism_fields["b_grad_delta_p"] = projected_gradient(
                pperp - ppar, bhat, lengths
            )[0]
            mechanism_fields["delta_p"] = pperp - ppar
            mechanism_fields["signed_pressure_stress_power_density"] = (
                -(pperp - ppar) * mechanism_fields["bb_grad_velocity"]
            )
        values.update({
            name: mechanism_fields[name] for name in MECHANISM_FIELD_DEFINITIONS
        })
    return values


def pressure_density_fields(
    fields: dict[str, np.ndarray],
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Construct Figure 2(a)-style pressure-versus-density coordinates."""

    rho = fields["dens"]
    ppar = fields["eint"]
    pperp = fields["p_perp"]
    mean_pressure = float(np.mean((2.0 * pperp + ppar) / 3.0))
    density_coordinate = mean_pressure * (rho / np.mean(rho) - 1.0)
    return {
        "parallel": (density_coordinate, ppar - np.mean(ppar)),
        "perpendicular": (density_coordinate, pperp - np.mean(pperp)),
    }


def analyze_fields(fields: dict[str, np.ndarray], lengths: tuple[float, float, float],
                   time: float, bins: int, alignment_shells: list[int],
                   pdf_ranges: dict[str, tuple[float, float]] | None = None,
                   model_choices: dict[str, object] | None = None,
                   eddy_samples: int = 0, eddy_bins: int = 20,
                   eddy_seed: int = 0,
                   joint_ranges: dict[
                       str, tuple[tuple[float, float], tuple[float, float]]
                   ] | None = None
                   ) -> dict[str, object]:
    """Analyze a single snapshot already represented as field arrays."""

    rho = fields["dens"]
    velocity = [fields["velx"], fields["vely"], fields["velz"]]
    magnetic = [fields["bcc1"], fields["bcc2"], fields["bcc3"]]
    ppar = fields["eint"]
    pperp = fields["p_perp"]
    bsqr = sum(component * component for component in magnetic)
    delta_p = pperp - ppar
    dk = 2.0 * math.pi / lengths[2]
    bhat = [component / np.sqrt(np.maximum(bsqr, np.finfo(float).tiny))
            for component in magnetic]
    gradient_parallel, gradient_perp = projected_gradient(delta_p, bhat, lengths)
    velocity_products = velocity_gradient_products(velocity, bhat, lengths)
    mechanism_fields = {
        **velocity_products,
        "delta_p": delta_p,
        "b_grad_delta_p": gradient_parallel,
        "signed_pressure_stress_power_density": (
            -delta_p * velocity_products["bb_grad_velocity"]
        ),
    }
    transfer = pressure_transfer(rho, velocity, magnetic, delta_p, lengths, dk)
    pressure_work = pressure_work_decomposition(
        fields,
        lengths,
        velocity_products["bb_grad_velocity"],
        model_choices,
        velocity_products["div_velocity"],
    )
    transfer.update({
        "applied_to_flow": pressure_work["applied_to_flow"],
        "interpretation": (
            f"{pressure_work['interpretation']}; signed pressure-stress exchange "
            "reconstructed from one retained snapshot"
        ),
    })
    anisotropic_power = float(pressure_work["anisotropic_stress_power"])
    transfer_direct = float(transfer["direct_real_space"])
    transfer_difference = transfer_direct - anisotropic_power
    pressure_work.update({
        "mks24_transfer_direct_real_space": transfer_direct,
        "mks24_transfer_minus_anisotropic_stress_power": transfer_difference,
        "mks24_transfer_relative_difference": float(
            transfer_difference / max(
                abs(transfer_direct), abs(anisotropic_power), np.finfo(float).tiny
            )
        ),
    })
    pdf_values = pdf_fields(fields, lengths, mechanism_fields)
    pressure_density_values = pressure_density_fields(fields)
    parallel_gradient_spectrum = shell_spectrum(
        [gradient_parallel],
        lengths,
        dk,
        field_definition=MECHANISM_FIELD_DEFINITIONS["b_grad_delta_p"],
    )
    spectra = {
        "velocity": shell_spectrum(
            velocity, lengths, dk, field_definition="velocity vector u"
        ),
        "compressive_velocity": compressive_velocity_spectrum(
            velocity, lengths, dk
        ),
        "magnetic_fluctuation": shell_spectrum(
            magnetic, lengths, dk, field_definition="magnetic field B"
        ),
        "density": shell_spectrum(
            [rho], lengths, dk, field_definition="mass density rho"
        ),
        "density_fluctuation": shell_spectrum(
            [rho / np.mean(rho) - 1.0], lengths, dk,
            field_definition="normalized density fluctuation rho/<rho> - 1",
        ),
        "p_parallel": shell_spectrum(
            [ppar], lengths, dk, field_definition="parallel thermal pressure p_parallel"
        ),
        "p_perp": shell_spectrum(
            [pperp], lengths, dk, field_definition="perpendicular thermal pressure p_perp"
        ),
        "magnetic_pressure": shell_spectrum(
            [0.5 * bsqr], lengths, dk,
            field_definition="magnetic pressure B^2/2 in AthenaK units",
        ),
        "delta_p": shell_spectrum(
            [delta_p],
            lengths,
            dk,
            field_definition=MECHANISM_FIELD_DEFINITIONS["delta_p"],
        ),
        "signed_pressure_stress_power_density": shell_spectrum(
            [mechanism_fields["signed_pressure_stress_power_density"]],
            lengths,
            dk,
            field_definition=MECHANISM_FIELD_DEFINITIONS[
                "signed_pressure_stress_power_density"
            ],
        ),
        "grad_parallel_delta_p": dict(parallel_gradient_spectrum),
        "b_grad_delta_p": dict(parallel_gradient_spectrum),
        "grad_perp_delta_p": shell_spectrum(
            [gradient_perp], lengths, dk,
            field_definition="magnitude of the gradient of Delta p perpendicular to b",
        ),
    }
    spectra.update({
        name: shell_spectrum(
            [values],
            lengths,
            dk,
            field_definition=MECHANISM_FIELD_DEFINITIONS.get(name),
        )
        for name, values in velocity_products.items()
    })
    return {
        "time": time,
        "shape_z_y_x": list(rho.shape),
        "lengths_x_y_z": list(lengths),
        "pdf": {
            name: pdf(
                values, bins, None if pdf_ranges is None else pdf_ranges.get(name)
            )
            for name, values in pdf_values.items()
        },
        "pressure_density_joint": {
            "definition": (
                "joint PDF in MKS24 Figure 2(a) panel coordinates: "
                "x = <p> delta rho/<rho>, y = delta p_parallel or delta p_perp"
            ),
            "mean_pressure_definition": "<p> = <(2 p_perp + p_parallel)/3>",
            "reference_scope": (
                "AthenaK-coordinate diagnostic; direct paper-raster comparison "
                "requires a qualified donor pressure/color-density transform"
            ),
            **{
                name: joint_pdf(
                    x_values, y_values, bins,
                    None if joint_ranges is None else joint_ranges[name],
                )
                for name, (x_values, y_values) in pressure_density_values.items()
            },
        },
        "spectra": spectra,
        "pressure_transfer": transfer,
        "mechanism_diagnostics": mechanism_diagnostic_index(
            normalized_transfer_available=transfer["normalization_available"]
        ),
        "mechanism_joint_diagnostics": mechanism_joint_diagnostics(
            mechanism_fields,
            bins,
            joint_ranges,
            pressure_work["applied_to_flow"],
            pressure_work["interpretation"],
        ),
        "pressure_work_decomposition": pressure_work,
        "alignment": alignment_histograms(velocity, magnetic, lengths, dk,
                                          alignment_shells, bins),
        "eddy_anisotropy": local_field_eddy_anisotropy(
            velocity, magnetic, lengths, eddy_samples, eddy_bins, eddy_seed
        ),
        "heat_flux_transport_proxy": heat_flux_transport_proxy(
            fields, lengths, model_choices
        ),
        "firehose_threshold_occupancy": firehose_threshold_occupancy_from_beta_delta(
            pdf_values["beta_delta"]
        ),
    }


def mean_distribution(records: list[dict[str, object]]) -> dict[str, object]:
    """Average compatible histogram products with already shared edges."""

    return {
        "edges": records[0]["edges"],
        "density": np.mean(
            [np.asarray(record["density"], dtype=float) for record in records], axis=0
        ).tolist(),
    }


def mean_joint_distribution(
    records: list[dict[str, object]],
    snapshot_times: list[float] | np.ndarray | None = None,
) -> dict[str, object]:
    """Average compatible two-dimensional histogram products."""

    densities = [
        np.asarray(record["density"], dtype=float) for record in records
    ]
    x_edges = np.asarray(records[0]["x_edges"], dtype=float)
    y_edges = np.asarray(records[0]["y_edges"], dtype=float)
    if any(
        not np.array_equal(np.asarray(record["x_edges"], dtype=float), x_edges)
        or not np.array_equal(np.asarray(record["y_edges"], dtype=float), y_edges)
        for record in records[1:]
    ):
        raise ValueError("joint distributions require identical shared bin edges")
    result: dict[str, object] = {
        "x_edges": x_edges.tolist(),
        "y_edges": y_edges.tolist(),
        "density": np.mean(densities, axis=0).tolist(),
    }
    if all("bin_counts" in record for record in records):
        bin_counts = [
            np.asarray(record["bin_counts"], dtype=np.int64) for record in records
        ]
        result.update({
            "bin_counts_sum": np.sum(bin_counts, axis=0).tolist(),
            "sample_count_sum": int(sum(
                int(record.get("sample_count", np.sum(counts)))
                for record, counts in zip(records, bin_counts)
            )),
            "binned_sample_count_sum": int(sum(
                int(record.get("binned_sample_count", np.sum(counts)))
                for record, counts in zip(records, bin_counts)
            )),
        })
    if snapshot_times is not None:
        result["uncertainty"] = descriptive_snapshot_block_uncertainty(
            densities, snapshot_times
        )
    return result


def mechanism_diagnostic_index(
    uncertainty_available: bool = False,
    ensemble: bool = False,
    normalized_transfer_available: bool = True,
) -> dict[str, object]:
    """Describe stable output paths for retained-snapshot mechanism diagnostics."""

    fields = {
        name: {
            "definition": definition,
            "pdf_path": f"pdf.{name}",
            "scale_resolved_path": f"spectra.{name}",
            **({
                "scale_resolved_uncertainty_path": f"spectra.{name}.uncertainty",
            } if uncertainty_available else {}),
        }
        for name, definition in MECHANISM_FIELD_DEFINITIONS.items()
    }
    transfer = {
        "curve_path": "pressure_transfer.signed_transfer",
        "normalized_curve_path": (
            "pressure_transfer.signed_transfer_normalized_by_total"
        ),
        "normalized_curve_available": normalized_transfer_available,
        "applied_to_flow_path": "pressure_transfer.applied_to_flow",
        "direct_real_space_path": (
            "pressure_transfer.direct_real_space_mean"
            if ensemble else "pressure_transfer.direct_real_space"
        ),
        "sign_convention_path": "pressure_transfer.sign_convention",
    }
    if uncertainty_available:
        transfer["curve_uncertainty_path"] = (
            "pressure_transfer.uncertainty.signed_transfer"
        )
        if normalized_transfer_available:
            transfer["normalized_curve_uncertainty_path"] = (
                "pressure_transfer.uncertainty.signed_transfer_normalized_by_total"
            )
        transfer["direct_real_space_uncertainty_path"] = (
            "pressure_transfer.uncertainty.direct_real_space"
        )
    return {
        "schema_version": MECHANISM_DIAGNOSTIC_SCHEMA_VERSION,
        "scope": "descriptive diagnostics reconstructed from retained snapshots",
        "discretization": (
            "cell-centered periodic gradients and perpendicular Fourier-shell "
            "partitions; not applied finite-volume face or stage accounting"
        ),
        "fields": fields,
        "signed_pressure_stress_transfer": transfer,
        "joint_conditional_products": {
            name: {
                "x_field": x_name,
                "y_field": y_name,
                "path": f"mechanism_joint_diagnostics.products.{name}",
            }
            for name, (x_name, y_name) in MECHANISM_JOINT_FIELDS.items()
        },
    }


def descriptive_snapshot_block_uncertainty(
    values: list[object] | np.ndarray,
    snapshot_times: list[float] | np.ndarray,
) -> dict[str, object]:
    """Summarize within-realization snapshot and contiguous-block variability."""

    times = np.asarray(snapshot_times, dtype=float)
    arrays = [np.asarray(value, dtype=float) for value in values]
    if len(arrays) != len(times) or not arrays:
        raise ValueError("uncertainty values and snapshot times must be nonempty and match")
    if not np.isfinite(times).all() or (
        len(times) > 1 and np.any(np.diff(times) < 0.0)
    ):
        raise ValueError("uncertainty snapshot times must be finite and nondecreasing")
    try:
        stacked = np.stack(arrays, axis=0)
    except ValueError as error:
        raise ValueError("uncertainty values have incompatible shapes") from error
    if not np.isfinite(stacked).all():
        raise ValueError("uncertainty values must be finite")
    result: dict[str, object] = {
        "schema_version": DESCRIPTIVE_UNCERTAINTY_SCHEMA_VERSION,
        "scope": (
            "descriptive within-realization temporal variability; not "
            "realization-to-realization or population uncertainty"
        ),
        "snapshot_weighting": "equal weight per retained snapshot",
        "snapshot_count": len(times),
        "snapshot_times": times.tolist(),
        "unique_snapshot_time_count": len(np.unique(times)),
        "duplicate_snapshot_times_present": len(np.unique(times)) != len(times),
        "available": len(times) >= 2,
    }
    if len(times) < 2:
        result.update({
            "reason": "at least two retained snapshots are required",
            "equal_snapshot_standard_deviation": None,
            "equal_snapshot_standard_error": None,
            "contiguous_blocks": {
                "available": False,
                "reason": "at least two retained snapshots are required",
            },
        })
        return result

    snapshot_standard_deviation = np.std(stacked, axis=0, ddof=1)
    block_count = min(
        DESCRIPTIVE_BLOCK_COUNT_MAX,
        max(2, len(times) // 2),
    )
    block_indices = np.array_split(np.arange(len(times)), block_count)
    block_means = np.stack(
        [np.mean(stacked[indices], axis=0) for indices in block_indices],
        axis=0,
    )
    block_standard_deviation = np.std(block_means, axis=0, ddof=1)
    result.update({
        "equal_snapshot_standard_deviation": snapshot_standard_deviation.tolist(),
        "equal_snapshot_standard_error": (
            snapshot_standard_deviation / math.sqrt(len(times))
        ).tolist(),
        "contiguous_blocks": {
            "available": True,
            "method": (
                "deterministic partition into up to four nonoverlapping contiguous "
                "equal-snapshot-count blocks, retaining at least two snapshots per "
                "block when the snapshot count permits; block sizes differ by at "
                "most one"
            ),
            "interpretation": (
                "descriptive sensitivity to contiguous temporal aggregation; "
                "not an autocorrelation-calibrated confidence interval"
            ),
            "maximum_block_count": DESCRIPTIVE_BLOCK_COUNT_MAX,
            "block_count": block_count,
            "block_snapshot_counts": [len(indices) for indices in block_indices],
            "block_time_ranges": [
                [float(times[indices[0]]), float(times[indices[-1]])]
                for indices in block_indices
            ],
            "block_mean_standard_deviation": block_standard_deviation.tolist(),
            "block_mean_standard_error": (
                block_standard_deviation / math.sqrt(block_count)
            ).tolist(),
        },
    })
    return result


def mean_spectrum(
    records: list[dict[str, object]],
    snapshot_times: list[float] | np.ndarray | None = None,
) -> dict[str, object]:
    """Average compatible shell-summed spectra."""

    powers = [
        np.asarray(record["power_per_dk"], dtype=float) for record in records
    ]
    result: dict[str, object] = {
        "dk": records[0]["dk"],
        "perpendicular": records[0]["perpendicular"],
        "k": records[0]["k"],
        "power_per_dk": np.mean(powers, axis=0).tolist(),
    }
    for key in ("field_definition", "normalization_definition"):
        if key in records[0]:
            result[key] = records[0][key]
    if snapshot_times is not None:
        result["uncertainty"] = descriptive_snapshot_block_uncertainty(
            powers, snapshot_times
        )
    return result


def mean_mechanism_joint_diagnostics(
    records: list[dict[str, object]],
    snapshot_times: list[float] | np.ndarray,
) -> dict[str, object]:
    """Average compatible mechanism joint PDFs and derive conditional summaries."""

    if not records or any(
        record.get("schema_version") != MECHANISM_JOINT_SCHEMA_VERSION
        for record in records
    ):
        raise ValueError("mechanism joint diagnostics are missing or incompatible")
    first_products = records[0].get("products")
    if not isinstance(first_products, dict) or set(first_products) != set(
        MECHANISM_JOINT_FIELDS
    ):
        raise ValueError("mechanism joint diagnostic product inventory differs")
    applied_values = [record.get("applied_to_flow") for record in records]
    applied_to_flow = (
        applied_values[0]
        if all(value == applied_values[0] for value in applied_values[1:])
        else None
    )
    interpretations = [str(record.get("interpretation", "")) for record in records]
    interpretation = (
        interpretations[0]
        if all(value == interpretations[0] for value in interpretations[1:])
        else "mixed model feedback scopes"
    )
    products: dict[str, object] = {}
    for name, (x_name, y_name) in MECHANISM_JOINT_FIELDS.items():
        source_products = [
            record["products"][name] for record in records
        ]
        if any(
            source["x_field"] != x_name or source["y_field"] != y_name
            for source in source_products
        ):
            raise ValueError(f"mechanism joint diagnostic fields differ: {name}")
        histogram = mean_joint_distribution(
            [source["joint_pdf"] for source in source_products],
            snapshot_times,
        )
        products[name] = {
            "x_field": x_name,
            "y_field": y_name,
            "x_definition": MECHANISM_FIELD_DEFINITIONS[x_name],
            "y_definition": MECHANISM_FIELD_DEFINITIONS[y_name],
            "joint_pdf": histogram,
            "conditional_y_given_x": conditional_profile_from_joint_pdf(histogram),
        }
    return {
        "schema_version": MECHANISM_JOINT_SCHEMA_VERSION,
        "scope": (
            "equal-snapshot-mean descriptive relationships among local mechanism "
            "diagnostics"
        ),
        "discretization": records[0].get(
            "discretization",
            "cell-centered periodic-gradient reconstruction",
        ),
        "applied_to_flow": applied_to_flow,
        "interpretation": interpretation,
        "local_signed_exchange_field": records[0].get(
            "local_signed_exchange_field",
            "signed_pressure_stress_power_density",
        ),
        "scale_resolved_signed_transfer_path": records[0].get(
            "scale_resolved_signed_transfer_path",
            "pressure_transfer.signed_transfer",
        ),
        "snapshot_count": len(records),
        "snapshot_times": np.asarray(snapshot_times, dtype=float).tolist(),
        "products": products,
    }


def mean_eddy_anisotropy(records: list[dict[str, object]]) -> dict[str, object]:
    """Average sampled structure functions and invert the ensemble result."""

    computed = [record for record in records if record.get("computed", False)]
    if not computed:
        return {
            "available": False,
            "reason": "eddy anisotropy was not computed for selected snapshots",
        }
    centers = np.asarray(computed[0]["bin_centers_over_lperp"], dtype=float)
    output = {
        name: computed[0][name] for name in (
            "definition", "conditioning", "sampling", "separation_coordinate",
            "normalization_definition", "lperp", "angle_degrees",
            "samples_requested", "seed", "bins",
        )
    }
    output["computed"] = True
    output["snapshot_count"] = len(computed)
    output["samples_retained"] = int(sum(
        int(record["samples_retained"]) for record in computed
    ))
    output["bin_centers_over_lperp"] = centers.tolist()
    for name in ("velocity_perp", "magnetic_perp"):
        product = {}
        for direction in ("perpendicular", "parallel"):
            count_name = f"{direction}_sample_counts"
            value_name = f"{direction}_structure_function"
            counts = np.sum([
                np.asarray(record[name][count_name], dtype=float)
                for record in computed
            ], axis=0)
            weighted = np.sum([
                np.nan_to_num(np.asarray(record[name][value_name], dtype=float))
                * np.asarray(record[name][count_name], dtype=float)
                for record in computed
            ], axis=0)
            values = np.zeros(counts.shape, dtype=float)
            populated = counts > 0.0
            values[populated] = weighted[populated] / counts[populated]
            product[value_name] = values
            product[count_name] = counts
        curve = eddy_anisotropy_curve(
            centers,
            product["perpendicular_structure_function"],
            product["parallel_structure_function"],
        )
        curve.update({
            key: value.tolist() for key, value in product.items()
        })
        output[name] = curve
    output["available"] = bool(
        output["velocity_perp"]["available"]
        or output["magnetic_perp"]["available"]
    )
    if not output["available"]:
        output["reason"] = "no ensemble eddy-anisotropy curve could be inverted"
    return output


def average_firehose_threshold_occupancy(
    samples: list[dict[str, object]],
) -> dict[str, object]:
    """Time-average compatible microinstability occupancy measurements."""

    if not samples:
        raise ValueError("microinstability occupancy requires snapshots")
    times = np.asarray([float(sample["time"]) for sample in samples], dtype=float)
    if not np.isfinite(times).all():
        raise ValueError("microinstability occupancy snapshot times are nonfinite")
    if len(times) > 1 and np.any(np.diff(times) <= 0.0):
        raise ValueError(
            "microinstability occupancy snapshot times must be strictly increasing"
        )
    shapes = [tuple(int(value) for value in sample["shape_z_y_x"]) for sample in samples]
    lengths = [
        tuple(float(value) for value in sample["lengths_x_y_z"])
        for sample in samples
    ]
    if any(shape != shapes[0] for shape in shapes[1:]) or any(
        length != lengths[0] for length in lengths[1:]
    ):
        raise ValueError("microinstability occupancy snapshots have inconsistent grids")
    if any(value <= 0.0 or not math.isfinite(value) for value in lengths[0]):
        raise ValueError("microinstability occupancy snapshot grid lengths are invalid")
    total_cells = math.prod(shapes[0])
    products = [
        validate_instantaneous_firehose_occupancy(
            sample["firehose_threshold_occupancy"],
            total_cells,
            f"microinstability snapshot {index}",
        )
        for index, sample in enumerate(samples)
    ]
    thresholds = products[0]["threshold_definitions"]
    physical_volume = math.prod(lengths[0])
    spatial_normalization = {
        **products[0]["normalization"],
        "weighting": (
            "equal physical cell volumes on the reconstructed uniform "
            "Cartesian snapshot grid"
        ),
        "physical_volume": physical_volume,
        "cell_volume": physical_volume / total_cells,
    }

    provenance_values = [sample.get("snapshot_provenance") for sample in samples]
    if any(value is not None for value in provenance_values) and not all(
        isinstance(value, dict) for value in provenance_values
    ):
        raise ValueError("microinstability occupancy snapshot provenance is incomplete")
    snapshot_provenance = [
        {**value, "snapshot_time": float(sample["time"])}
        for sample, value in zip(samples, provenance_values)
        if isinstance(value, dict)
    ]

    occupancy: dict[str, dict[str, object]] = {}
    for semantics in INSTABILITY_OCCUPANCY_SEMANTICS:
        measured: dict[str, object] = {}
        for name in INSTABILITY_OCCUPANCY_COMPONENTS:
            fractions = np.asarray(
                [
                    float(product["occupancy"][semantics][name]["volume_fraction"])
                    for product in products
                ],
                dtype=float,
            )
            if not np.isfinite(fractions).all() or np.any(
                (fractions < 0.0) | (fractions > 1.0)
            ):
                raise ValueError(
                    f"microinstability occupancy fractions are invalid: "
                    f"{semantics}.{name}"
                )
            if len(times) > 1:
                time_average = (
                    trapezoidal_integral(fractions.tolist(), times.tolist())
                    / float(times[-1] - times[0])
                )
                time_average_method = "trapezoidal_time_average"
            else:
                time_average = float(fractions[0])
                time_average_method = "single_snapshot"
            equal_snapshot_mean = float(np.mean(fractions))
            measured[name] = {
                "snapshot_volume_fractions": fractions.tolist(),
                "snapshot_fraction_minimum": float(np.min(fractions)),
                "snapshot_fraction_maximum": float(np.max(fractions)),
                "comparison_volume_fraction": float(time_average),
                "time_average_volume_fraction": float(time_average),
                "equal_snapshot_mean_volume_fraction": equal_snapshot_mean,
            }
        for parallel_name, oblique_name in (
            ("parallel_firehose", "oblique_firehose"),
            ("parallel_total_unstable", "oblique_total_unstable"),
        ):
            parallel = np.asarray(
                measured[parallel_name]["snapshot_volume_fractions"], dtype=float
            )
            oblique = np.asarray(
                measured[oblique_name]["snapshot_volume_fractions"], dtype=float
            )
            reclassified = np.asarray(
                measured["oblique_only_reclassification"][
                    "snapshot_volume_fractions"
                ],
                dtype=float,
            )
            if not np.allclose(
                oblique - parallel, reclassified, rtol=0.0, atol=1.0e-15
            ):
                raise ValueError(
                    "microinstability occupancy reclassification fractions are "
                    f"inconsistent: {semantics}.{oblique_name}"
                )
        mirror = np.asarray(
            measured["mirror"]["snapshot_volume_fractions"], dtype=float
        )
        for total_name, firehose_name in (
            ("parallel_total_unstable", "parallel_firehose"),
            ("oblique_total_unstable", "oblique_firehose"),
        ):
            total_fraction = np.asarray(
                measured[total_name]["snapshot_volume_fractions"], dtype=float
            )
            firehose_fraction = np.asarray(
                measured[firehose_name]["snapshot_volume_fractions"], dtype=float
            )
            if not np.allclose(
                total_fraction,
                mirror + firehose_fraction,
                rtol=0.0,
                atol=1.0e-15,
            ):
                raise ValueError(
                    "microinstability total occupancy is inconsistent: "
                    f"{semantics}.{total_name}"
                )
        occupancy[semantics] = measured

    result: dict[str, object] = {
        "schema_version": FIGURE_13_FIREHOSE_SCHEMA_VERSION,
        "definition": (
            "late-window retained-snapshot mirror-or-firehose occupancy under "
            "parallel and alternative oblique firehose definitions"
        ),
        "normalization": {
            "spatial": spatial_normalization,
            "temporal": {
                "primary_comparison_kind": time_average_method,
                "comparison_definition": (
                    "trapezoidal time average of instantaneous equal-cell-volume "
                    "fractions over retained snapshots"
                    if len(times) > 1 else
                    "instantaneous equal-cell-volume fraction from one snapshot"
                ),
                "secondary_kind": "equal_weight_retained_snapshot_mean",
                "secondary_definition": (
                    "arithmetic mean of instantaneous retained-snapshot fractions"
                ),
            },
        },
        "threshold_definitions": thresholds,
        "grid": {
            "shape_z_y_x": list(shapes[0]),
            "lengths_x_y_z": list(lengths[0]),
        },
        "analysis_window": {
            "selected_time_first": float(times[0]),
            "selected_time_last": float(times[-1]),
            "snapshot_count": len(samples),
            "snapshot_times": times.tolist(),
        },
        "occupancy": occupancy,
    }
    if snapshot_provenance:
        result["snapshot_provenance"] = snapshot_provenance
        result["snapshot_provenance_sha256"] = canonical_json_sha256(
            snapshot_provenance
        )
    return result


def average_snapshot_records(records: dict[str, dict[str, object]]) -> dict[str, object]:
    """Average spatial products selected from one case and one time window."""

    if not records:
        return {"snapshot_count": 0}
    samples = sorted(records.values(), key=lambda sample: float(sample["time"]))
    times = np.asarray([float(sample["time"]) for sample in samples], dtype=float)
    can_integrate = bool(len(times) >= 2 and times[-1] > times[0])
    pdf_names = samples[0]["pdf"].keys()
    pressure_density = [
        sample["pressure_density_joint"] for sample in samples
    ]
    joint_names = ("parallel", "perpendicular")
    spectrum_names = samples[0]["spectra"].keys()
    alignment_names = set(samples[0]["alignment"].keys())
    for sample in samples[1:]:
        alignment_names.intersection_update(sample["alignment"].keys())
    transfer = [sample["pressure_transfer"] for sample in samples]
    normalized_transfer_available = all(
        item.get("normalization_available", False) for item in transfer
    )
    pressure_work = [sample["pressure_work_decomposition"] for sample in samples]
    heat_flux = [
        sample["heat_flux_transport_proxy"] for sample in samples
        if sample["heat_flux_transport_proxy"].get("available", False)
    ]
    shared_work_scope = pressure_work[0]["interpretation"]
    if any(sample["interpretation"] != shared_work_scope for sample in pressure_work[1:]):
        shared_work_scope = "mixed model feedback scopes"
    pressure_work_ensemble: dict[str, object] = {
        "available": True,
        "snapshot_count": len(pressure_work),
        "definition": pressure_work[0]["definition"],
        "anisotropic_stress_definition": pressure_work[0][
            "anisotropic_stress_definition"
        ],
        "discretization": pressure_work[0]["discretization"],
        "sign_convention": pressure_work[0]["sign_convention"],
        "applied_to_flow": pressure_work[0]["applied_to_flow"] if all(
            sample["applied_to_flow"] == pressure_work[0]["applied_to_flow"]
            for sample in pressure_work[1:]
        ) else None,
        "interpretation": shared_work_scope,
    }
    for name in (
        "isotropic_perpendicular_pressure_power",
        "anisotropic_stress_power",
        "total_cgl_pressure_power",
        "parallel_strain_rms",
        "anisotropic_power_density_rms",
        "mks24_transfer_direct_real_space",
        "mks24_transfer_minus_anisotropic_stress_power",
        "mks24_transfer_relative_difference",
    ):
        pressure_work_ensemble[f"{name}_mean"] = float(
            np.mean([sample[name] for sample in pressure_work])
        )
    pressure_work_ensemble["uncertainty"] = {
        name: descriptive_snapshot_block_uncertainty(
            [sample[name] for sample in pressure_work], times
        )
        for name in (
            "isotropic_perpendicular_pressure_power",
            "anisotropic_stress_power",
            "total_cgl_pressure_power",
            "parallel_strain_rms",
            "anisotropic_power_density_rms",
            "mks24_transfer_direct_real_space",
        )
    }
    pressure_integral: dict[str, object] = {
        "available": can_integrate,
        "definition": (
            "trapezoidal time integral of snapshot-reconstructed CGL "
            "pressure-power terms"
        ),
        "discretization": (
            "sparse retained-snapshot quadrature; not applied stage accounting"
        ),
        "snapshot_count": len(times),
        "time_first": float(times[0]),
        "time_last": float(times[-1]),
    }
    if can_integrate:
        for name in (
            "isotropic_perpendicular_pressure_power",
            "anisotropic_stress_power",
            "total_cgl_pressure_power",
            "mks24_transfer_direct_real_space",
            "mks24_transfer_minus_anisotropic_stress_power",
        ):
            pressure_integral[f"{name}_integral"] = trapezoidal_integral(
                [sample[name] for sample in pressure_work], times
            )
    pressure_work_ensemble["time_integral_estimate"] = pressure_integral
    heat_flux_ensemble: dict[str, object] = {
        "available": bool(heat_flux),
        "snapshot_count": len(heat_flux),
    }
    if heat_flux:
        heat_flux_ensemble.update({
            "definition": heat_flux[0]["definition"],
            "discretization": heat_flux[0]["discretization"],
            "closure_model_choices": heat_flux[0]["closure_model_choices"],
        })
        for name in (
            "valid_volume_fraction",
            "regularized_parallel_power",
            "regularized_perpendicular_power",
            "regularized_total_power",
            "unlimited_parallel_power",
            "unlimited_perpendicular_power",
            "unlimited_total_power",
            "parallel_cap_active_volume_fraction",
            "perpendicular_cap_active_volume_fraction",
            "regularized_parallel_power_on_cap_active_cells",
            "regularized_perpendicular_power_on_cap_active_cells",
        ):
            heat_flux_ensemble[f"{name}_mean"] = float(
                np.mean([sample[name] for sample in heat_flux])
            )
        heat_flux_times = np.asarray([
            float(sample["time"]) for sample in samples
            if sample["heat_flux_transport_proxy"].get("available", False)
        ])
        heat_flux_can_integrate = bool(
            len(heat_flux_times) >= 2 and heat_flux_times[-1] > heat_flux_times[0]
        )
        heat_flux_integral: dict[str, object] = {
            "available": heat_flux_can_integrate,
            "definition": (
                "trapezoidal time integral of snapshot-reconstructed LF "
                "heat-flux smoothing-power terms"
            ),
            "discretization": (
                "sparse retained-snapshot quadrature; not applied face accounting"
            ),
            "snapshot_count": len(heat_flux_times),
            "time_first": float(heat_flux_times[0]),
            "time_last": float(heat_flux_times[-1]),
        }
        if heat_flux_can_integrate:
            for name in (
                "regularized_parallel_power",
                "regularized_perpendicular_power",
                "regularized_total_power",
                "unlimited_total_power",
            ):
                heat_flux_integral[f"{name}_integral"] = trapezoidal_integral(
                    [sample[name] for sample in heat_flux], heat_flux_times
                )
        heat_flux_ensemble["time_integral_estimate"] = heat_flux_integral
    spectra_ensemble = {
        name: mean_spectrum(
            [sample["spectra"][name] for sample in samples], times
        )
        for name in spectrum_names
    }
    transfer_values = [
        np.asarray(item["transfer"], dtype=float) for item in transfer
    ]
    signed_transfer_values = [
        np.asarray(item.get("signed_transfer", item["transfer"]), dtype=float)
        for item in transfer
    ]
    transfer_mean = np.mean(transfer_values, axis=0).tolist()
    signed_transfer_mean = np.mean(signed_transfer_values, axis=0).tolist()
    pressure_transfer_ensemble: dict[str, object] = {
        "dk": transfer[0]["dk"],
        "k_perp": transfer[0]["k_perp"],
        "transfer": transfer_mean,
        "signed_transfer": list(signed_transfer_mean),
        "definition": transfer[0].get(
            "definition",
            "signed perpendicular-shell pressure-stress transfer partition",
        ),
        "filter_definition": transfer[0].get(
            "filter_definition",
            "sqrt(rho) u filtered into nonoverlapping k_perp shells",
        ),
        "sign_convention": transfer[0].get(
            "sign_convention",
            "positive transfer is kinetic-energy gain from the pressure force",
        ),
        "interpretation": (
            f"{shared_work_scope}; signed pressure-stress exchange reconstructed "
            "from retained snapshots"
        ),
        "applied_to_flow": transfer[0].get("applied_to_flow") if all(
            item.get("applied_to_flow") == transfer[0].get("applied_to_flow")
            for item in transfer[1:]
        ) else None,
        "normalization_available": normalized_transfer_available,
        "normalization_definition": transfer[0]["normalization_definition"],
        "kinetic_energy_mean": float(np.mean(
            [item["kinetic_energy"] for item in transfer]
        )),
        "velocity_rms_mean": float(np.mean(
            [item["velocity_rms"] for item in transfer]
        )),
        "perpendicular_outer_scale": transfer[0]["perpendicular_outer_scale"],
        "total_transfer_rate_mean": float(np.mean(
            [item["total_transfer_rate"] for item in transfer]
        )),
        "transfer_normalized_by_total": (
            np.mean(
                [
                    np.asarray(item["transfer_normalized_by_total"], dtype=float)
                    for item in transfer
                ],
                axis=0,
            ).tolist()
            if normalized_transfer_available else None
        ),
        "signed_transfer_normalized_by_total": (
            np.mean(
                [
                    np.asarray(
                        item.get(
                            "signed_transfer_normalized_by_total",
                            item["transfer_normalized_by_total"],
                        ),
                        dtype=float,
                    )
                    for item in transfer
                ],
                axis=0,
            ).tolist()
            if normalized_transfer_available else None
        ),
        "direct_real_space_mean": float(np.mean(
            [item["direct_real_space"] for item in transfer]
        )),
        "shell_sum_mean": float(np.mean([item["shell_sum"] for item in transfer])),
        "closure_error_mean": float(np.mean(
            [item["closure_error"] for item in transfer]
        )),
        "uncertainty": {
            "signed_transfer": descriptive_snapshot_block_uncertainty(
                signed_transfer_values, times
            ),
            "direct_real_space": descriptive_snapshot_block_uncertainty(
                [item["direct_real_space"] for item in transfer], times
            ),
            "shell_sum": descriptive_snapshot_block_uncertainty(
                [item["shell_sum"] for item in transfer], times
            ),
        },
    }
    if normalized_transfer_available:
        pressure_transfer_ensemble["uncertainty"][
            "signed_transfer_normalized_by_total"
        ] = descriptive_snapshot_block_uncertainty(
            [
                item.get(
                    "signed_transfer_normalized_by_total",
                    item["transfer_normalized_by_total"],
                )
                for item in transfer
            ],
            times,
        )
    mechanism_ensemble: dict[str, object] = {}
    if all(
        isinstance(sample.get("mechanism_joint_diagnostics"), dict)
        for sample in samples
    ):
        mechanism_ensemble = {
            "mechanism_diagnostics": mechanism_diagnostic_index(
                uncertainty_available=True,
                ensemble=True,
                normalized_transfer_available=normalized_transfer_available,
            ),
            "mechanism_joint_diagnostics": mean_mechanism_joint_diagnostics(
                [sample["mechanism_joint_diagnostics"] for sample in samples],
                times,
            ),
        }
    return {
        "snapshot_count": len(samples),
        "time_first": min(float(sample["time"]) for sample in samples),
        "time_last": max(float(sample["time"]) for sample in samples),
        "snapshot_times": times.tolist(),
        "pdf": {
            name: mean_distribution([sample["pdf"][name] for sample in samples])
            for name in pdf_names
        },
        "pressure_density_joint": {
            "definition": pressure_density[0]["definition"],
            "mean_pressure_definition": pressure_density[0][
                "mean_pressure_definition"
            ],
            "reference_scope": pressure_density[0]["reference_scope"],
            **{
                name: mean_joint_distribution([
                    sample[name] for sample in pressure_density
                ], times)
                for name in joint_names
            },
        },
        "spectra": spectra_ensemble,
        "pressure_transfer": pressure_transfer_ensemble,
        **mechanism_ensemble,
        "pressure_work_decomposition": pressure_work_ensemble,
        "alignment": {
            name: mean_distribution([sample["alignment"][name] for sample in samples])
            for name in sorted(alignment_names)
        },
        "eddy_anisotropy": mean_eddy_anisotropy([
            sample["eddy_anisotropy"] for sample in samples
        ]),
        "heat_flux_transport_proxy": heat_flux_ensemble,
        "firehose_threshold_occupancy": average_firehose_threshold_occupancy(samples),
    }


def analyze_snapshot_paths(paths: list[Path], bins: int, alignment_shells: list[int],
                           time_start: float | None = None,
                           time_end: float | None = None,
                           model_choices: dict[str, object] | None = None,
                           eddy_samples: int = 0, eddy_bins: int = 20,
                           eddy_seed: int = 0,
                           expected_ranks_by_path: dict[str, int] | None = None,
                           ) -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    """Analyze and time-average selected snapshots with common PDF bin edges."""

    selected: list[Path] = []
    extrema: dict[str, list[float]] = {}
    joint_extrema: dict[str, list[list[float]]] = {}
    expected_ranks_by_path = expected_ranks_by_path or {}
    candidate_rank_sets = {
        str(path): snapshot_sibling_paths(
            path, expected_ranks_by_path.get(str(path))
        )
        for path in paths
    }
    for path in paths:
        time = snapshot_time(path)
        if not time_mask(np.asarray([time]), time_start, time_end)[0]:
            continue
        selected.append(path)
    if len({str(path) for path in selected}) != len(selected):
        raise ValueError("selected snapshot paths must be unique")
    exact_rank_sets = {
        str(path): candidate_rank_sets[str(path)]
        for path in selected
    }
    snapshot_provenance = {
        str(path): snapshot_digest_provenance(
            path,
            expected_ranks_by_path.get(str(path)),
            exact_rank_sets[str(path)],
        )
        for path in selected
    }

    def read_selected(path: Path) -> tuple[
        dict[str, np.ndarray], tuple[float, float, float], float
    ]:
        if path.parent.name == "rank_00000000":
            return read_snapshot(path, exact_rank_sets[str(path)])
        return read_snapshot(path)

    for path in selected:
        fields, lengths, _ = read_selected(path)
        scalar_fields = pdf_fields(fields, lengths)
        for name, values in scalar_fields.items():
            finite = values[np.isfinite(values)]
            if finite.size == 0:
                continue
            low, high = float(np.min(finite)), float(np.max(finite))
            if name in extrema:
                extrema[name][0] = min(extrema[name][0], low)
                extrema[name][1] = max(extrema[name][1], high)
            else:
                extrema[name] = [low, high]
        joint_coordinates = {
            **pressure_density_fields(fields),
            **mechanism_joint_coordinates(scalar_fields),
        }
        for name, (x_values, y_values) in joint_coordinates.items():
            coordinates = (x_values, y_values)
            if name not in joint_extrema:
                joint_extrema[name] = [
                    [float(np.min(values)), float(np.max(values))]
                    for values in coordinates
                ]
                continue
            for index, values in enumerate(coordinates):
                joint_extrema[name][index][0] = min(
                    joint_extrema[name][index][0], float(np.min(values))
                )
                joint_extrema[name][index][1] = max(
                    joint_extrema[name][index][1], float(np.max(values))
                )
    ranges: dict[str, tuple[float, float]] = {}
    for name, (low, high) in extrema.items():
        if high <= low:
            delta = max(abs(low), 1.0) * 1.0e-12
            low -= delta
            high += delta
        ranges[name] = (low, high)
    joint_ranges: dict[
        str, tuple[tuple[float, float], tuple[float, float]]
    ] = {}
    for name, coordinates in joint_extrema.items():
        padded: list[tuple[float, float]] = []
        for low, high in coordinates:
            if high <= low:
                delta = max(abs(low), 1.0) * 1.0e-12
                low -= delta
                high += delta
            padded.append((low, high))
        joint_ranges[name] = (padded[0], padded[1])
    records: dict[str, dict[str, object]] = {}
    for path in selected:
        fields, lengths, time = read_selected(path)
        record = analyze_fields(
            fields, lengths, time, bins, alignment_shells, ranges, model_choices,
            eddy_samples, eddy_bins, eddy_seed, joint_ranges
        )
        record["snapshot_provenance"] = snapshot_provenance[str(path)]
        records[str(path)] = record
    for path in selected:
        if snapshot_digest_provenance(
            path,
            expected_ranks_by_path.get(str(path)),
            exact_rank_sets[str(path)],
        ) != snapshot_provenance[str(path)]:
            raise ValueError(f"snapshot changed while being analyzed: {path}")
    ensemble = average_snapshot_records(records)
    ensemble["time_start"] = time_start
    ensemble["time_end"] = time_end
    if "firehose_threshold_occupancy" in ensemble:
        ensemble["firehose_threshold_occupancy"]["analysis_window"].update({
            "requested_time_start": time_start,
            "requested_time_end": time_end,
        })
    return records, ensemble


def combine_firehose_occupancy_products(
    products: list[dict[str, object]], expected_total_cells: int,
) -> dict[str, object]:
    """Combine disjoint meshblock occupancy products by exact integer counts."""

    if not products:
        raise ValueError("streamed Figure 13 occupancy has no meshblock products")
    combined = json.loads(json.dumps(products[0]))
    total = sum(int(product["normalization"]["total_cell_count"]) for product in products)
    if total != expected_total_cells:
        raise ValueError("streamed Figure 13 occupancy does not cover the full grid")
    combined["normalization"]["total_cell_count"] = total
    combined["beta_delta_minimum"] = min(
        float(product["beta_delta_minimum"]) for product in products
    )
    combined["beta_delta_maximum"] = max(
        float(product["beta_delta_maximum"]) for product in products
    )
    for semantics in INSTABILITY_OCCUPANCY_SEMANTICS:
        for name in INSTABILITY_OCCUPANCY_COMPONENTS:
            count = sum(
                int(product["occupancy"][semantics][name]["cell_count"])
                for product in products
            )
            combined["occupancy"][semantics][name] = {
                "cell_count": count,
                "volume_fraction": count / total,
            }
    return validate_instantaneous_firehose_occupancy(
        combined, total, "streamed Figure 13 snapshot"
    )


def streamed_figure_13_snapshot_occupancy(
    path: Path, rank_files: list[Path],
) -> dict[str, object]:
    """Measure Figure 13 occupancy from production blocks without full-grid FFTs."""

    metadata_keys = (
        "header", "time", "cycle", "var_names", "Nx1", "Nx2", "Nx3", "nvars",
        "x1min", "x1max", "x2min", "x2max", "x3min", "x3max",
        "nx1_mb", "nx2_mb", "nx3_mb", "nx1_out_mb", "nx2_out_mb", "nx3_out_mb",
    )
    occupancy_products: list[dict[str, object]] = []
    logical: list[np.ndarray] = []
    indices: list[np.ndarray] = []
    geometry: list[np.ndarray] = []
    reference: dict[str, object] | None = None
    for rank, rank_file in enumerate(rank_files):
        raw = bin_convert.read_binary(str(rank_file))
        if reference is None:
            reference = raw
            missing = sorted(set(REQUIRED_FIELDS) - set(raw["var_names"]))
            if missing:
                raise ValueError(
                    f"{path} lacks CGL paper fields {missing}; regenerate with "
                    "mhd_w_bcc output"
                )
        else:
            for key in metadata_keys:
                if not binary_metadata_equal(reference[key], raw[key]):
                    raise ValueError(
                        "rank-local snapshot metadata differs across exact rank set: "
                        f"rank={rank}, key={key}"
                    )
        fields = {
            name: np.asarray(raw["mb_data"][name])
            for name in ("eint", "p_perp", "bcc1", "bcc2", "bcc3")
        }
        occupancy_products.append(firehose_threshold_occupancy(fields))
        logical.append(np.asarray(raw["mb_logical"]))
        indices.append(np.asarray(raw["mb_index"]))
        geometry.append(np.asarray(raw["mb_geometry"]))
    assert reference is not None
    combined_inventory = reference.copy()
    combined_inventory["mb_logical"] = np.concatenate(logical, axis=0)
    combined_inventory["mb_index"] = np.concatenate(indices, axis=0)
    combined_inventory["mb_geometry"] = np.concatenate(geometry, axis=0)
    combined_inventory["n_mbs"] = len(combined_inventory["mb_logical"])
    combined_inventory["mb_data"] = {
        name: np.empty((combined_inventory["n_mbs"], 0))
        for name in reference["var_names"]
    }
    validate_exact_rank_set_meshblocks(combined_inventory)
    shape = (
        int(reference["Nx3"]),
        int(reference["Nx2"]),
        int(reference["Nx1"]),
    )
    lengths = (
        float(reference["x1max"] - reference["x1min"]),
        float(reference["x2max"] - reference["x2min"]),
        float(reference["x3max"] - reference["x3min"]),
    )
    return {
        "time": float(reference["time"]),
        "shape_z_y_x": list(shape),
        "lengths_x_y_z": list(lengths),
        "firehose_threshold_occupancy": combine_firehose_occupancy_products(
            occupancy_products, math.prod(shape)
        ),
    }


def analyze_figure_13_snapshot_paths(
    paths: list[Path],
    time_start: float | None,
    time_end: float | None,
    expected_ranks_by_path: dict[str, int],
    provenance_cache: dict[str, dict[str, object]],
    occupancy_cache: dict[str, dict[str, object]],
    revalidated_cache: set[str],
) -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    """Run the comparison-ready Figure 13-only streamed CLI analysis."""

    selected = [
        path for path in paths
        if time_mask(np.asarray([snapshot_time(path)]), time_start, time_end)[0]
    ]
    if len({str(path) for path in selected}) != len(selected):
        raise ValueError("selected snapshot paths must be unique")
    records: dict[str, dict[str, object]] = {}
    for path in selected:
        rank_files = snapshot_sibling_paths(
            path, expected_ranks_by_path.get(str(path))
        )
        cache_key = canonical_json_sha256({
            "representative": str(path),
            "rank_files": [str(value) for value in rank_files],
            "expected_ranks": expected_ranks_by_path.get(str(path)),
        })
        provenance = provenance_cache.get(cache_key)
        if provenance is None:
            provenance = snapshot_digest_provenance(
                path, expected_ranks_by_path.get(str(path)), rank_files
            )
            provenance_cache[cache_key] = provenance
        aggregate = str(provenance["aggregate_sha256"])
        record = occupancy_cache.get(aggregate)
        if record is None:
            record = streamed_figure_13_snapshot_occupancy(path, rank_files)
            occupancy_cache[aggregate] = record
        record = json.loads(json.dumps(record))
        record["snapshot_provenance"] = provenance
        records[str(path)] = record
        if aggregate not in revalidated_cache:
            current = snapshot_digest_provenance(
                path, expected_ranks_by_path.get(str(path)), rank_files
            )
            if current != provenance:
                raise ValueError(f"snapshot changed while being analyzed: {path}")
            revalidated_cache.add(aggregate)
    ensemble = {
        "snapshot_count": len(records),
        "time_first": min(float(record["time"]) for record in records.values()),
        "time_last": max(float(record["time"]) for record in records.values()),
        "firehose_threshold_occupancy": average_firehose_threshold_occupancy(
            list(records.values())
        ),
    }
    ensemble["firehose_threshold_occupancy"]["analysis_window"].update({
        "requested_time_start": time_start,
        "requested_time_end": time_end,
    })
    return records, ensemble


def synthetic_test() -> dict[str, object]:
    """Exercise binning, projected gradients, transfer closure, and alignment bounds."""

    lengths = (1.0, 1.0, 2.0)
    shape = (32, 16, 16)
    z = (np.arange(shape[0]) + 0.5) * lengths[2] / shape[0]
    y = (np.arange(shape[1]) + 0.5) * lengths[1] / shape[1]
    x = (np.arange(shape[2]) + 0.5) * lengths[0] / shape[2]
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    fields = {
        "dens": np.ones(shape),
        "velx": np.sin(2.0 * math.pi * xx / lengths[0]),
        "vely": np.zeros(shape),
        "velz": np.zeros(shape),
        "eint": np.ones(shape),
        "p_perp": 1.0 + 0.1 * np.sin(2.0 * math.pi * zz / lengths[2]),
        "bcc1": np.zeros(shape),
        "bcc2": np.zeros(shape),
        "bcc3": np.ones(shape),
    }
    model = {
        "lf_k_parallel": str(math.pi),
        "lf_coefficient_mode": "local",
        "nu_coll": "0.0",
        "mirror_limiter": "false",
        "firehose_limiter": "false",
        "cgl_firehose_threshold": "parallel",
        "limiter_nu_coll": "0.0",
        "backup_limiters": "false",
        "dfloor": "1.0e-12",
        "pfloor": "1.0e-12",
        "tfloor": "1.0e-12",
        "bfloor": "1.0e-10",
        "passive_delta": "false",
    }
    record = analyze_fields(fields, lengths, 0.0, 16, [2], model_choices=model)
    velocity_power = np.asarray(record["spectra"]["velocity"]["power_per_dk"])
    compressive_power = np.asarray(
        record["spectra"]["compressive_velocity"]["power_per_dk"]
    )
    peak = int(np.argmax(velocity_power))
    transverse_fields = {
        name: np.array(values, copy=True) for name, values in fields.items()
    }
    transverse_fields["velx"] = np.zeros(shape)
    transverse_fields["vely"] = np.sin(2.0 * math.pi * xx / lengths[0])
    transverse_record = analyze_fields(
        transverse_fields, lengths, 0.0, 16, [2], model_choices=model
    )
    transverse_compressive_power = np.asarray(
        transverse_record["spectra"]["compressive_velocity"]["power_per_dk"]
    )
    derivative = periodic_gradient(fields["p_perp"] - fields["eint"], lengths)[2]
    exact = 0.1 * math.pi * np.cos(math.pi * zz)
    relative_gradient_error = float(
        np.max(np.abs(derivative - exact)) / np.max(np.abs(exact))
    )
    transfer = record["pressure_transfer"]
    pressure_work = record["pressure_work_decomposition"]
    alignment = record["alignment"].get("2", {})
    strain = np.asarray(record["pdf"]["bb_grad_velocity"]["density"])
    pressure_density = record["pressure_density_joint"]
    heat_flux = record["heat_flux_transport_proxy"]
    finite_alignment = bool(
        alignment and np.isfinite(np.asarray(alignment["density"])).all()
    )
    finite_pressure_density = bool(all(
        np.isfinite(np.asarray(pressure_density[name]["density"], dtype=float)).all()
        for name in ("parallel", "perpendicular")
    ))
    joint_fields = {
        name: np.array(values, copy=True) for name, values in fields.items()
    }
    density_fluctuation = 0.05 * np.sin(2.0 * math.pi * xx / lengths[0])
    joint_fields["dens"] = 1.0 + density_fluctuation
    joint_fields["eint"] = 1.0 + (5.0 / 3.0) * density_fluctuation
    joint_fields["p_perp"] = 1.0 + (5.0 / 3.0) * density_fluctuation
    joint_coordinates = pressure_density_fields(joint_fields)["parallel"]
    joint_guide_error = float(np.max(np.abs(
        joint_coordinates[1] - (5.0 / 3.0) * joint_coordinates[0]
    )))
    correlated_fields = {
        name: np.array(values, copy=True) for name, values in fields.items()
    }
    correlated_fields["velx"] = np.zeros(shape)
    correlated_fields["velz"] = np.sin(math.pi * zz)
    correlated_fields["p_perp"] = 1.0 + 0.1 * np.cos(math.pi * zz)
    correlated_record = analyze_fields(
        correlated_fields, lengths, 0.0, 16, [2], model_choices=model
    )
    correlated = correlated_record["pressure_work_decomposition"]
    correlated_transfer = correlated_record["pressure_transfer"]
    passive_model = dict(model)
    passive_model["passive_delta"] = "true"
    passive_work = analyze_fields(
        correlated_fields, lengths, 0.0, 16, [2], model_choices=passive_model
    )["pressure_work_decomposition"]
    repeated_correlated = analyze_fields(
        correlated_fields, lengths, 2.0, 16, [2], model_choices=model
    )
    correlated_quadrature = average_snapshot_records({
        "early": correlated_record,
        "late": repeated_correlated,
    })["pressure_work_decomposition"]["time_integral_estimate"]
    eddy_fields = {
        name: np.array(values, copy=True) for name, values in fields.items()
    }
    eddy_fields["velx"] = (
        np.sin(2.0 * math.pi * xx) + np.sin(math.pi * zz)
    )
    eddy_fields["vely"] = (
        np.sin(2.0 * math.pi * yy) + np.sin(math.pi * zz)
    )
    eddy_fields["bcc1"] = 0.1 * (
        np.sin(2.0 * math.pi * xx) + np.sin(math.pi * zz)
    )
    eddy_fields["bcc2"] = 0.1 * (
        np.sin(2.0 * math.pi * yy) + np.sin(math.pi * zz)
    )
    eddy_record = analyze_fields(
        eddy_fields, lengths, 0.0, 16, [2], model_choices=model,
        eddy_samples=200000, eddy_bins=10, eddy_seed=731
    )
    eddy = eddy_record["eddy_anisotropy"]
    finite_eddy_anisotropy = bool(
        eddy["available"]
        and eddy["velocity_perp"]["available"]
        and eddy["magnetic_perp"]["available"]
        and np.isfinite(np.asarray(
            eddy["velocity_perp"]["ell_parallel_over_lperp"]
        )).all()
        and np.isfinite(np.asarray(
            eddy["magnetic_perp"]["ell_parallel_over_lperp"]
        )).all()
    )
    passed = bool(
        peak == 2
        and float(np.max(compressive_power)) > 0.0
        and float(np.max(np.abs(transverse_compressive_power))) < 1.0e-28
        and relative_gradient_error < 7.0e-3
        and abs(float(transfer["direct_real_space"])) < 1.0e-14
        and abs(float(transfer["closure_error"])) < 1.0e-14
        and transfer["normalization_available"]
        and np.isfinite(
            np.asarray(transfer["transfer_normalized_by_total"], dtype=float)
        ).all()
        and finite_alignment
        and finite_pressure_density
        and joint_guide_error < 1.0e-14
        and np.isfinite(strain).all()
        and heat_flux["available"]
        and heat_flux["regularized_perpendicular_power"] > 0.0
        and abs(float(heat_flux["regularized_parallel_power"])) < 1.0e-14
        and abs(float(pressure_work["anisotropic_stress_power"])) < 1.0e-14
        and correlated["anisotropic_stress_power"] < 0.0
        and correlated_transfer["normalization_available"]
        and abs(float(
            correlated["mks24_transfer_minus_anisotropic_stress_power"]
        )) < 1.0e-14
        and passive_work["applied_to_flow"] is False
        and correlated_quadrature["available"]
        and abs(float(correlated_quadrature["anisotropic_stress_power_integral"])
                - 2.0 * float(correlated["anisotropic_stress_power"])) < 1.0e-14
        and finite_eddy_anisotropy
    )
    return {
        "passed": passed,
        "velocity_peak_kperp_bin": peak,
        "expected_velocity_peak_kperp_bin": 2,
        "positive_longitudinal_compressive_spectrum": float(
            np.max(compressive_power)
        ),
        "zero_transverse_compressive_spectrum": float(
            np.max(np.abs(transverse_compressive_power))
        ),
        "parallel_gradient_relative_error": relative_gradient_error,
        "zero_transfer": transfer["direct_real_space"],
        "transfer_closure_error": transfer["closure_error"],
        "finite_normalized_transfer": bool(
            transfer["normalization_available"]
            and np.isfinite(
                np.asarray(transfer["transfer_normalized_by_total"], dtype=float)
            ).all()
        ),
        "positive_correlated_transfer_normalization": (
            correlated_transfer["total_transfer_rate"]
        ),
        "finite_alignment": finite_alignment,
        "finite_pressure_density_joint_pdf": finite_pressure_density,
        "pressure_density_five_thirds_coordinate_error": joint_guide_error,
        "finite_strain_pdf": bool(np.isfinite(strain).all()),
        "positive_perpendicular_heat_flux_proxy": (
            heat_flux["regularized_perpendicular_power"]
        ),
        "zero_parallel_heat_flux_proxy": heat_flux["regularized_parallel_power"],
        "zero_anisotropic_pressure_work": pressure_work["anisotropic_stress_power"],
        "negative_correlated_anisotropic_pressure_work": (
            correlated["anisotropic_stress_power"]
        ),
        "correlated_transfer_work_difference": (
            correlated["mks24_transfer_minus_anisotropic_stress_power"]
        ),
        "passive_pressure_work_is_diagnostic_only": (
            passive_work["applied_to_flow"] is False
        ),
        "constant_power_quadrature_error": (
            correlated_quadrature["anisotropic_stress_power_integral"]
            - 2.0 * correlated["anisotropic_stress_power"]
        ),
        "finite_eddy_anisotropy": finite_eddy_anisotropy,
    }


def optional_float(value: object) -> float | None:
    """Convert manifest values to optional finite floating-point numbers."""

    if value in (None, "unspecified"):
        return None
    result = float(str(value))
    if not math.isfinite(result):
        raise ValueError(f"invalid analysis-window value: {value}")
    return result


def case_window(case: dict[str, object], time_start: float | None,
                time_end: float | None) -> tuple[float | None, float | None]:
    """Select explicit CLI analysis limits or per-case archived model choices."""

    model = case.get("model_choices", {})
    if not isinstance(model, dict):
        model = {}
    start = time_start if time_start is not None else optional_float(
        model.get("analysis_t_start")
    )
    end = time_end if time_end is not None else optional_float(
        model.get("analysis_t_end")
    )
    return start, end


def bundle_manifest(bundle: Path) -> dict[str, object]:
    """Read one retained bundle manifest."""

    with (bundle / "manifest.json").open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    if not isinstance(manifest, dict):
        raise ValueError(f"bundle manifest must be an object: {bundle}")
    return manifest


def bundle_cases(bundle: Path,
                 manifest: dict[str, object] | None = None
                 ) -> list[dict[str, object]]:
    """Read case metadata used to preserve per-case analysis windows."""

    if manifest is None:
        manifest = bundle_manifest(bundle)
    cases = manifest.get("cases", [])
    if not isinstance(cases, list):
        raise ValueError(f"manifest cases must be a list: {bundle}")
    if not all(isinstance(case, dict) for case in cases):
        raise ValueError(f"manifest cases must contain objects: {bundle}")
    return cases


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of one file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_regular_file_sha256(path: Path, context: str) -> str:
    """Hash one retained regular file while rejecting concurrent mutation."""

    if not path.is_file() or path.is_symlink():
        raise ValueError(f"{context} is not a regular retained file: {path}")
    before = path.stat()
    digest = sha256_file(path)
    after = path.stat()
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise ValueError(f"{context} changed while being hashed: {path}")
    return digest


def stable_regular_file_bytes(path: Path, context: str) -> bytes:
    """Read one retained regular file while rejecting concurrent mutation."""

    if not path.is_file() or path.is_symlink():
        raise ValueError(f"{context} is not a regular retained file: {path}")
    before = path.stat()
    payload = path.read_bytes()
    after = path.stat()
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise ValueError(f"{context} changed while being read: {path}")
    return payload


def canonical_json_sha256(value: object) -> str:
    """Return a deterministic SHA-256 digest of one JSON-compatible value."""

    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def json_compatible(value: object) -> object:
    """Return one deterministic JSON-compatible representation."""

    if isinstance(value, dict):
        return {
            str(key): json_compatible(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [json_compatible(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted(json_compatible(item) for item in value)
    if isinstance(value, Path):
        return str(value)
    return value


def fsync_directory(path: Path) -> None:
    """Persist one atomic directory-entry update before returning."""

    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_write_json(destination: Path, value: object) -> None:
    """Atomically replace one named JSON output with fully persisted content."""

    atomic_write_json_generation([(destination, value)])


def atomic_write_json_generation(
    outputs: list[tuple[Path, object]],
) -> None:
    """Atomically replace one standalone JSON file."""

    if len(outputs) != 1:
        raise ValueError(
            "standalone atomic JSON replacement accepts exactly one output; "
            "multi-output publications require an immutable generation"
        )
    destination, value = outputs[0]
    parent = destination.parent
    parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink() or (
        destination.exists() and not destination.is_file()
    ):
        raise ValueError(
            f"atomic JSON destination is not a regular file: {destination}"
        )
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=parent,
            prefix=f".{destination.name}.stage.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            staged = Path(stream.name)
            stream.write(
                (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")
            )
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(staged, destination)
        staged = None
        fsync_directory(parent)
    finally:
        if staged is not None:
            staged.unlink(missing_ok=True)


def bin_convert_implementation_provenance() -> dict[str, object]:
    """Authenticate the imported binary reconstruction implementation."""

    expected_path = (ROOT_DIR / "vis" / "python" / "bin_convert.py").resolve()
    imported_path_value = getattr(bin_convert, "__file__", None)
    if not isinstance(imported_path_value, str):
        raise ValueError("Figure 13 bin_convert import has no source path")
    imported_path = Path(imported_path_value).resolve()
    if imported_path != expected_path:
        raise ValueError("Figure 13 imported a different bin_convert implementation")
    functions: dict[str, object] = {}
    for name in BIN_CONVERT_FUNCTIONS:
        function = getattr(bin_convert, name, None)
        code = getattr(function, "__code__", None)
        if (
            not callable(function)
            or getattr(function, "__module__", None) != "bin_convert"
            or code is None
            or Path(str(code.co_filename)).resolve() != expected_path
        ):
            raise ValueError(
                f"Figure 13 imported bin_convert function identity differs: {name}"
            )
        functions[name] = {
            "module": function.__module__,
            "source_path": str(Path(str(code.co_filename)).resolve()),
            "first_line": int(code.co_firstlineno),
        }
    return {
        "source_path": str(expected_path),
        "source_sha256": sha256_file(expected_path),
        "functions": functions,
        "functions_sha256": canonical_json_sha256(functions),
    }


def stage_i_validator_implementation_provenance() -> dict[str, object]:
    """Authenticate the imported closed Stage I restart-contract implementation."""

    expected_path = (
        ROOT_DIR / "scripts" / "frontier" / "cgl_lf_stage_i_validate_segment.py"
    ).resolve()
    imported_path_value = getattr(stage_i_validator, "__file__", None)
    if not isinstance(imported_path_value, str):
        raise ValueError("Figure 13 Stage I validator import has no source path")
    imported_path = Path(imported_path_value).resolve()
    if imported_path != expected_path:
        raise ValueError("Figure 13 imported a different Stage I validator")
    functions: dict[str, object] = {}
    for name in STAGE_I_VALIDATOR_FUNCTIONS:
        function = getattr(stage_i_validator, name, None)
        code = getattr(function, "__code__", None)
        if (
            not callable(function)
            or getattr(function, "__module__", None)
            != "cgl_lf_stage_i_validate_segment"
            or code is None
            or Path(str(code.co_filename)).resolve() != expected_path
        ):
            raise ValueError(
                f"Figure 13 imported Stage I validator function identity differs: "
                f"{name}"
            )
        functions[name] = {
            "module": function.__module__,
            "source_path": str(Path(str(code.co_filename)).resolve()),
            "first_line": int(code.co_firstlineno),
        }
    return {
        "source_path": str(expected_path),
        "source_sha256": sha256_file(expected_path),
        "functions": functions,
        "functions_sha256": canonical_json_sha256(functions),
    }


def analysis_invocation_provenance(
    configuration: dict[str, object],
) -> dict[str, object]:
    """Bind analyzer, binary reconstruction source, and normalized CLI."""

    analyzer_path = Path(__file__).resolve()
    return {
        "analyzer_path": str(analyzer_path),
        "analyzer_sha256": sha256_file(analyzer_path),
        "bin_convert": bin_convert_implementation_provenance(),
        "stage_i_validator": stage_i_validator_implementation_provenance(),
        "configuration": configuration,
        "configuration_sha256": canonical_json_sha256(configuration),
    }


def validate_analysis_invocation_provenance(
    provenance: dict[str, object],
) -> dict[str, object]:
    """Fail closed unless analyzer and CLI provenance match this invocation."""

    configuration = provenance.get("configuration")
    if not isinstance(configuration, dict):
        raise ValueError("Figure 13 analysis provenance requires configuration")
    if set(configuration) != set(ANALYSIS_CONFIGURATION_KEYS):
        raise ValueError(
            "Figure 13 analysis provenance does not bind the complete CLI "
            "configuration"
        )
    if provenance.get("configuration_sha256") != canonical_json_sha256(configuration):
        raise ValueError("Figure 13 analysis configuration checksum does not match")
    analyzer_path = Path(__file__).resolve()
    if provenance.get("analyzer_path") != str(analyzer_path):
        raise ValueError("Figure 13 analysis provenance names a different analyzer")
    if provenance.get("analyzer_sha256") != sha256_file(analyzer_path):
        raise ValueError("Figure 13 analyzer source checksum does not match")
    if provenance.get("bin_convert") != bin_convert_implementation_provenance():
        raise ValueError(
            "Figure 13 binary reconstruction implementation provenance does not match"
        )
    if (
        provenance.get("stage_i_validator")
        != stage_i_validator_implementation_provenance()
    ):
        raise ValueError(
            "Figure 13 Stage I restart-contract implementation provenance does not "
            "match"
        )
    return configuration


def figure_13_matrix_provenance() -> dict[str, object]:
    """Validate the tracked R03/R07/R14/R15 matrix and source inputs exactly."""

    matrix_path = ROOT_DIR / FIGURE_13_MATRIX_PATH
    matrix_digest = sha256_file(matrix_path)
    if matrix_digest != FIGURE_13_MATRIX_SHA256:
        raise ValueError("Figure 13 tracked Stage I matrix checksum does not match")
    with matrix_path.open(encoding="utf-8") as stream:
        matrix = json.load(stream)
    matrix_cases = matrix.get("cases", []) if isinstance(matrix, dict) else []
    if not isinstance(matrix_cases, list):
        raise ValueError("Figure 13 tracked Stage I matrix cases are malformed")
    indexed = {
        str(case.get("id")): case for case in matrix_cases if isinstance(case, dict)
    }
    retained: dict[str, object] = {}
    for case_id, contract in FIGURE_13_CASE_CONTRACTS.items():
        record = indexed.get(case_id)
        if not isinstance(record, dict):
            raise ValueError(f"Figure 13 tracked matrix lacks {case_id}")
        for key in ("name", "input", "resolution"):
            if record.get(key) != contract[key]:
                raise ValueError(
                    f"Figure 13 tracked matrix identity mismatch: {case_id}.{key}"
                )
        if "13" not in [str(value) for value in record.get("figure_roles", [])]:
            raise ValueError(f"Figure 13 tracked matrix lacks role 13 for {case_id}")
        source_input = ROOT_DIR / str(contract["input"])
        source_digest = sha256_file(source_input)
        if source_digest != contract["input_sha256"]:
            raise ValueError(f"Figure 13 source input checksum mismatch: {case_id}")
        retained[case_id] = {
            "matrix_identity": {
                "id": case_id,
                "name": contract["name"],
                "input": contract["input"],
                "resolution": contract["resolution"],
                "figure_role": "13",
            },
            "source_input_path": str(source_input),
            "source_input_sha256": source_digest,
        }
    return {
        "tracked_matrix_path": str(matrix_path),
        "tracked_matrix_sha256": matrix_digest,
        "cases": retained,
    }


def path_is_beneath(path: Path, root: Path) -> bool:
    """Return whether one resolved path is contained by another."""

    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return True


def bundle_relative_member_path(bundle: Path, value: str, context: str) -> Path:
    """Return one lexically relative bundle member with a contained parent."""

    relative = Path(value)
    if (
        relative.is_absolute()
        or relative == Path(".")
        or ".." in relative.parts
        or not path_is_beneath(bundle / relative.parent, bundle)
    ):
        raise ValueError(f"{context} path escapes its bundle")
    return bundle / relative


def stable_json_file(path: Path, context: str) -> tuple[dict[str, object], str]:
    """Read one regular JSON file while rejecting mutation during the read."""

    if not path.is_file() or path.is_symlink():
        raise ValueError(f"{context} is not a regular retained file: {path}")
    before = path.stat()
    payload = path.read_bytes()
    after = path.stat()
    before_identity = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    )
    after_identity = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    )
    if before_identity != after_identity:
        raise ValueError(f"{context} changed while being read: {path}")
    try:
        value = json.loads(payload)
    except json.JSONDecodeError as error:
        raise ValueError(f"{context} is not valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{context} must contain a JSON object: {path}")
    return value, hashlib.sha256(payload).hexdigest()


def authenticate_regular_file(
    path: Path,
    expected_sha256: object,
    expected_size: object,
    context: str,
    root: Path | None = None,
) -> dict[str, object]:
    """Authenticate one exact live regular file against retained byte evidence."""

    digest = str(expected_sha256)
    if (
        not path.is_absolute()
        or path != path.absolute()
        or path.is_symlink()
        or not path.is_file()
        or (root is not None and not path_is_beneath(path, root))
        or re.fullmatch(r"[0-9a-f]{64}", digest) is None
        or (
            expected_size is not None
            and (
                not isinstance(expected_size, int)
                or isinstance(expected_size, bool)
                or expected_size <= 0
            )
        )
    ):
        raise ValueError(f"{context} retained file evidence is invalid")
    actual_digest = stable_regular_file_sha256(path, context)
    actual_size = path.stat().st_size
    if actual_digest != digest or (
        expected_size is not None and actual_size != expected_size
    ):
        raise ValueError(f"{context} size or checksum does not match")
    return {
        "path": str(path.resolve()),
        "sha256": actual_digest,
        "size_bytes": actual_size,
    }


def restart_parameter_dump(path: Path) -> tuple[str, int]:
    """Read one Athena restart parameter dump and binary-payload offset."""

    marker = b"<par_end>\n"
    payload = b""
    with path.open("rb") as stream:
        while len(payload) < MAX_RESTART_PARAMETER_DUMP_BYTES:
            block = stream.read(
                min(4096, MAX_RESTART_PARAMETER_DUMP_BYTES - len(payload))
            )
            if not block:
                break
            payload += block
            end = payload.find(marker)
            if end >= 0:
                try:
                    return payload[:end].decode("utf-8"), end + len(marker)
                except UnicodeDecodeError as error:
                    raise ValueError(
                        f"restart parameter dump is not UTF-8: {path}"
                    ) from error
    raise ValueError(f"restart lacks loadable <par_end> terminator: {path}")


def restart_marker_time(text: str, path: Path) -> float:
    """Return the unique finite time/restart_time marker from a parameter dump."""

    block = ""
    markers: list[str] = []
    for original in text.splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            continue
        if block == "time" and "=" in line:
            key, value = line.split("=", 1)
            if key.strip() == "restart_time":
                markers.append(value.strip())
    if len(markers) != 1:
        raise ValueError(
            f"restart parameter dump must contain one time/restart_time marker: {path}"
        )
    try:
        marker_time = float(markers[0])
    except ValueError as error:
        raise ValueError(f"restart time marker is not numeric: {path}") from error
    if not math.isfinite(marker_time):
        raise ValueError(f"restart time marker is not finite: {path}")
    return marker_time


def turbulence_restart_metadata_contract(
    input_contract: dict[str, object],
) -> tuple[list[str], list[int | None], list[str], list[float]]:
    """Reconstruct Athena's exact TurbulenceRestartMetadata configuration."""

    blocks = dict(input_contract["archived_parameter_blocks"])
    driving = dict(blocks.get("turb_driving", {}))
    driving.update({
        key: str(value)
        for (block, key), value in (
            stage_i_validator.QUALIFIED_PRODUCT_PARAMETER_ADDITIONS.items()
        )
        if block == "turb_driving" and value is not None
    })

    def text(name: str, default: str) -> str:
        return str(driving.get(name, default))

    def integer(name: str, default: str) -> int:
        return int(text(name, default))

    def real(name: str, default: str) -> float:
        value = float(text(name, default))
        if not math.isfinite(value):
            raise ValueError(f"turbulence runtime parameter is nonfinite: {name}")
        return value

    def boolean(name: str, default: str) -> int:
        value = text(name, default)
        if value not in {"true", "false"}:
            raise ValueError(f"turbulence runtime parameter is not boolean: {name}")
        return int(value == "true")

    nlow = integer("nlow", "1")
    nhigh = integer("nhigh", "3")
    driving_type = integer("driving_type", "0")
    limits = {
        name: integer(name, "0" if name.startswith("min_") else str(nhigh))
        for name in ("min_kx", "max_kx", "min_ky", "max_ky", "min_kz", "max_kz")
    }
    tile = {
        name: integer(name, "1") for name in ("tile_nx", "tile_ny", "tile_nz")
    }
    enum_values = {
        "normalization": {"edot": 0, "accel_rms": 1},
        "localization": {"none": 0, "include": 1, "exclude": 2},
        "spectrum": {"parabolic": 0, "power_law": 1},
        "projection_policy": {
            "solenoidal_compressive": 0,
            "mks24_random_unprojected": 1,
            "mks24_alfvenic_perpendicular": 2,
        },
    }
    enums: dict[str, int] = {}
    for name, values in enum_values.items():
        value = text(name, next(iter(values)))
        if value not in values:
            raise ValueError(f"turbulence runtime enum is invalid: {name}")
        enums[name] = values[value]

    domain = tuple(
        float(upper) - float(lower)
        for lower, upper in input_contract["domain"]
    )
    tile_lengths = tuple(
        length / tile[name]
        for length, name in zip(domain, ("tile_nx", "tile_ny", "tile_nz"))
    )
    physical_shell = boolean("physical_k_shell", "false")
    k_shell_unit = real("k_shell_unit", "0")
    mode_count = 0
    for nkx in range(limits["min_kx"], limits["max_kx"] + 1):
        for nky in range(limits["min_ky"], limits["max_ky"] + 1):
            for nkz in range(limits["min_kz"], limits["max_kz"] + 1):
                if nkx == nky == nkz == 0:
                    continue
                if physical_shell:
                    dkx = 2.0 * math.pi / tile_lengths[0]
                    dky = 2.0 * math.pi / tile_lengths[1]
                    dkz = 2.0 * math.pi / tile_lengths[2]
                    normalized = (
                        (dkx * nkx) ** 2
                        + (dky * nky) ** 2
                        + (dkz * nkz) ** 2
                    ) / (k_shell_unit ** 2)
                    selected = nlow * nlow <= normalized <= nhigh * nhigh
                elif driving_type == 0:
                    shell = nkx * nkx + nky * nky + nkz * nkz
                    selected = nlow * nlow <= shell <= nhigh * nhigh
                else:
                    perpendicular = nkx * nkx + nky * nky
                    parallel = nkz * nkz
                    selected = (
                        nlow * nlow <= perpendicular <= nhigh * nhigh
                        and nlow * nlow <= parallel <= nhigh * nhigh
                    )
                mode_count += int(selected)
    use_npeak = int("npeak" in driving)
    npeak = real("npeak", "0")
    kpeak = (
        npeak * 2.0 * math.pi / tile_lengths[0]
        if use_npeak else real("kpeak", str(4.0 * math.pi))
    )
    turb_flag = integer("turb_flag", "2")
    tdriv_duration = (
        real("tdriv_duration", text("tcorr", "0"))
        if turb_flag == 1
        else float(np.finfo(np.float32).max)
    )
    int_names = [
        "version", "mode_count", "n_updates", "nlow", "nhigh", "driving_type",
        "min_kx", "max_kx", "min_ky", "max_ky", "min_kz", "max_kz",
        "use_npeak", "turb_flag", "tile_nx", "tile_ny", "tile_nz",
        "normalization", "localization", "spectrum", "projection_policy",
        "physical_k_shell", "isotropic_power_spectrum", "record_injected_work",
    ]
    int_values: list[int | None] = [
        3, mode_count, None, nlow, nhigh, driving_type,
        limits["min_kx"], limits["max_kx"], limits["min_ky"], limits["max_ky"],
        limits["min_kz"], limits["max_kz"], use_npeak, turb_flag,
        tile["tile_nx"], tile["tile_ny"], tile["tile_nz"],
        enums["normalization"], enums["localization"], enums["spectrum"],
        enums["projection_policy"], physical_shell,
        boolean("isotropic_power_spectrum", "false"),
        boolean("record_injected_work", "false"),
    ]
    real_names = [
        "tcorr", "dt_update", "dedt", "accel_rms", "sol_fraction", "kpeak",
        "npeak", "expo", "exp_prp", "exp_prl", "tdriv_duration", "tdriv_start",
        "sigma_x1", "sigma_x2", "sigma_x3", "center_x1", "center_x2",
        "center_x3", "k_shell_unit",
    ]
    real_values = [
        real("tcorr", "0"), real("dt_update", "0.01"),
        real("dedt", "0") if enums["normalization"] == 0 else 0.0,
        real("accel_rms", "0") if enums["normalization"] == 1 else 0.0,
        real("sol_fraction", "1"), kpeak, npeak, real("expo", str(5.0 / 3.0)),
        real("exp_prp", str(5.0 / 3.0)), real("exp_prl", "0"),
        tdriv_duration, real("tdriv_start", "0"), real("sigma_x1", "-1"),
        real("sigma_x2", "-1"), real("sigma_x3", "-1"),
        real("center_x1", "0"), real("center_x2", "0"), real("center_x3", "0"),
        k_shell_unit,
    ]
    return int_names, int_values, real_names, real_values


def validate_restart_turbulence_metadata(
    path: Path,
    restart_binary_abi: dict[str, object],
    input_contract: dict[str, object],
) -> dict[str, object]:
    """Bind every serialized turbulence configuration field to runtime input."""

    contract = restart_binary_abi["contract"]
    _, parameter_size = restart_parameter_dump(path)
    expected_blocks = int(input_contract["expected_meshblocks"])
    metadata_offset = (
        parameter_size
        + int(contract["mesh_header_size"])
        + expected_blocks * int(contract["logical_location_size"])
        + expected_blocks * int(contract["cost_size"])
    )
    before = path.stat()
    with path.open("rb") as stream:
        stream.seek(metadata_offset)
        payload = stream.read(int(contract["turbulence_metadata_size"]))
    after = path.stat()
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise ValueError(f"restart changed during turbulence metadata validation: {path}")
    if len(payload) != int(contract["turbulence_metadata_size"]):
        raise ValueError(f"restart turbulence metadata is truncated: {path}")
    actual_ints = list(struct.unpack_from("<24i", payload))
    actual_reals = list(struct.unpack_from("<19d", payload, 96))
    int_names, expected_ints, real_names, expected_reals = (
        turbulence_restart_metadata_contract(input_contract)
    )
    mismatches = [
        name for name, actual, expected in zip(int_names, actual_ints, expected_ints)
        if expected is not None and actual != expected
    ]
    mismatches.extend(
        name for name, actual, expected in zip(real_names, actual_reals, expected_reals)
        if actual != expected
    )
    if mismatches:
        raise ValueError(
            f"restart turbulence metadata differs from exact runtime configuration "
            f"for {mismatches}: {path}"
        )
    return {
        "validation": "all-runtime-configuration-fields-exact",
        "configuration_fields": {
            **{
                name: actual for name, actual in zip(int_names, actual_ints)
                if name != "n_updates"
            },
            **dict(zip(real_names, actual_reals)),
        },
        "dynamic_state": {"n_updates": actual_ints[2]},
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def restart_loadability_evidence(
    path: Path,
    expected_time: float,
    restart_binary_abi: dict[str, object],
    input_contract: dict[str, object],
    runtime_contract: dict[str, object],
) -> dict[str, object]:
    """Apply the independent validator's exact production restart contract."""

    contract = restart_binary_abi.get("contract")
    if (
        not isinstance(contract, dict)
        or restart_binary_abi.get("contract_sha256")
        != canonical_json_sha256(contract)
    ):
        raise ValueError(f"restart binary ABI evidence is malformed: {path}")
    profiles: dict[Path, tuple[int, ...]] = {}
    try:
        evidence = stage_i_validator.restart_binary_evidence(
            path, contract, input_contract, runtime_contract, profiles
        )
        stage_i_validator.recheck_profiles(profiles)
    except stage_i_validator.ValidationError as error:
        raise ValueError(f"restart exact production loadability failed: {error}") from error
    binary_time = float(evidence["time"])
    if not figure_13_time_close(binary_time, expected_time):
        raise ValueError(f"restart binary time differs from expected time: {path}")
    metadata = validate_restart_turbulence_metadata(
        path, restart_binary_abi, input_contract
    )
    text, parameter_size = restart_parameter_dump(path)
    return {
        "binary_layout_validation": "independent-validator-exact-production-contract",
        "restart_binary_abi_sha256": restart_binary_abi["contract_sha256"],
        "parameter_dump_size": parameter_size,
        "parameter_dump_sha256": hashlib.sha256(
            (text + "<par_end>\n").encode("utf-8")
        ).hexdigest(),
        "meshblock_count": int(input_contract["expected_meshblocks"]),
        "local_meshblock_count": int(evidence["local_blocks"]),
        "binary_time": binary_time,
        "dt": float(evidence["dt"]),
        "cycle": int(evidence["cycle"]),
        "marker_mode": str(evidence["marker_mode"]),
        "logical_location_count": len(evidence["locations"]),
        "mode_count": int(evidence["mode_count"]),
        "variable_data_size": int(input_contract["restart_data_size"]),
        "logical_location_inventory_sha256": evidence[
            "logical_location_inventory_sha256"
        ],
        "cost_inventory_sha256": evidence["cost_inventory_sha256"],
        "turbulence_metadata_sha256": evidence["turbulence_metadata_sha256"],
        "turbulence_metadata_contract": metadata,
        "rng_lifecycle": evidence["rng_lifecycle"],
        "rng_n_updates": int(evidence["rng_n_updates"]),
        "rng_canonical_continuation_state_sha256": evidence[
            "rng_canonical_continuation_state_sha256"
        ],
        "rng_authenticated_field_count": int(evidence["rng_authenticated_field_count"]),
        "turbulence_amplitudes_sha256": evidence["turbulence_amplitudes_sha256"],
        "injected_work": float(evidence["injected_work"]),
        "lf_diagnostics": [float(value) for value in evidence["lf_diagnostics"]],
        "mesh_region_contract": evidence["mesh_region_contract"],
    }


def authenticate_restart_rank_files(
    rank_files: object,
    expected_ranks: int,
    root: Path,
    expected_time: float,
    restart_binary_abi: dict[str, object],
    input_contract: dict[str, object],
    runtime_contract: dict[str, object],
    context: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Authenticate exact restart siblings, bytes, rank identity, and loadability."""

    if not isinstance(rank_files, list) or len(rank_files) != expected_ranks:
        raise ValueError(f"{context} lacks the exact expected restart rank count")
    records: list[dict[str, object]] = []
    loadability: list[dict[str, object]] = []
    for record in rank_files:
        if not isinstance(record, dict):
            raise ValueError(f"{context} has malformed restart rank evidence")
        authenticated = authenticate_regular_file(
            Path(str(record.get("path", ""))),
            record.get("sha256"),
            record.get("size_bytes"),
            context,
            root,
        )
        records.append(authenticated)
        loadability.append(
            restart_loadability_evidence(
                Path(authenticated["path"]),
                expected_time,
                restart_binary_abi,
                input_contract,
                runtime_contract,
            )
        )
        if authenticate_regular_file(
            Path(authenticated["path"]),
            authenticated["sha256"],
            authenticated["size_bytes"],
            context,
            root,
        ) != authenticated:
            raise ValueError(f"{context} changed during exact loadability validation")
    paths = [Path(str(record["path"])) for record in records]
    expected_names = [f"rank_{rank:08d}" for rank in range(expected_ranks)]
    if (
        len(set(paths)) != expected_ranks
        or [path.parent.name for path in paths] != expected_names
        or len({path.name for path in paths}) != 1
        or any(
            any(
                evidence[key] != loadability[0][key]
                for key in (
                    "binary_time",
                    "meshblock_count",
                    "dt",
                    "cycle",
                    "logical_location_count",
                    "mode_count",
                    "variable_data_size",
                    "restart_binary_abi_sha256",
                    "parameter_dump_sha256",
                    "logical_location_inventory_sha256",
                    "cost_inventory_sha256",
                    "turbulence_metadata_sha256",
                    "turbulence_metadata_contract",
                    "rng_lifecycle",
                    "rng_n_updates",
                    "rng_canonical_continuation_state_sha256",
                    "rng_authenticated_field_count",
                    "turbulence_amplitudes_sha256",
                )
            )
            for evidence in loadability[1:]
        )
        or sum(
            int(evidence["local_meshblock_count"]) for evidence in loadability
        ) != int(loadability[0]["meshblock_count"])
    ):
        raise ValueError(f"{context} restart rank siblings are incomplete or ambiguous")
    return records, loadability


def qualified_restart_binary_abi(
    git_revision: str, executable_sha256: str, context: str,
) -> dict[str, object]:
    """Return the exact native-restart ABI qualified for one executable."""

    contract = FIGURE_13_QUALIFIED_RESTART_BINARY_ABIS.get(
        (git_revision, executable_sha256)
    )
    if not isinstance(contract, dict):
        raise ValueError(f"{context} executable lacks a qualified restart binary ABI")
    retained: dict[str, object] = {
        "executable_revision": git_revision,
        "executable_sha256": executable_sha256,
        "contract": json_compatible(contract),
    }
    retained["contract_sha256"] = canonical_json_sha256(retained["contract"])
    return retained


def authenticate_qualified_restart_load_evidence(
    approval: dict[str, object],
    canonical_root: Path,
    executable: dict[str, object],
    git_revision: str,
    context: str,
) -> dict[str, object]:
    """Authenticate the approved natural-restart continuation qualification."""

    notes = approval.get("review_notes")
    if not isinstance(notes, str):
        raise ValueError(f"{context} approval lacks qualified restart-load evidence")
    match = re.search(r"(?:^|[; ])g027 ([0-9a-f]{64})(?:;|$)", notes)
    if match is None:
        raise ValueError(f"{context} approval lacks exact g027 restart-load evidence")
    path = canonical_root / FIGURE_13_QUALIFIED_RESTART_LOAD_EVIDENCE_RELATIVE
    evidence, digest = stable_json_file(path, f"{context} restart-load evidence")
    acceptance = evidence.get("acceptance")
    execution = evidence.get("execution")
    provenance = evidence.get("provenance")
    review = evidence.get("review")
    checks = evidence.get("checks")
    executable_evidence = (
        provenance.get("executable") if isinstance(provenance, dict) else None
    )
    staged = (
        acceptance.get("staged_restart_boundary")
        if isinstance(acceptance, dict) else None
    )
    terminal = (
        acceptance.get("terminal_identity")
        if isinstance(acceptance, dict) else None
    )
    terminal_metadata = (
        checks.get("terminal_restart_metadata") if isinstance(checks, dict) else None
    )
    try:
        staged_time = float(staged["restart_time"])
        terminal_time = float(terminal_metadata["restart_time"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            f"{context} qualified restart-load evidence is malformed"
        ) from error
    if (
        digest != match.group(1)
        or evidence.get("schema_version") != 2
        or not isinstance(execution, dict)
        or execution.get("git_revision") != git_revision
        or execution.get("result") != "passed"
        or not isinstance(executable_evidence, dict)
        or executable_evidence.get("path") != executable["path"]
        or executable_evidence.get("sha256") != executable["sha256"]
        or not isinstance(acceptance, dict)
        or acceptance.get("qualification_role")
        != FIGURE_13_QUALIFIED_RESTART_LOAD_ROLE
        or not isinstance(staged, dict)
        or staged.get("content_set_matches_g026_selected_boundary") is not True
        or not isinstance(terminal, dict)
        or terminal.get("within_tolerance") is not True
        or not isinstance(review, dict)
        or review.get("derived_by_read_only_inspection") is not True
        or review.get("local_review") != "passed"
        or not math.isfinite(staged_time)
        or not math.isfinite(terminal_time)
        or staged_time <= 0.0
        or terminal_time <= staged_time
    ):
        raise ValueError(f"{context} qualified restart-load evidence does not bind")
    return {
        "path": str(path.resolve()),
        "sha256": digest,
        "qualification_role": FIGURE_13_QUALIFIED_RESTART_LOAD_ROLE,
        "staged_restart_time": staged_time,
        "terminal_restart_time": terminal_time,
        "terminal_identity_within_tolerance": True,
    }


def authenticate_approved_executable(
    command: dict[str, object],
    canonical_root: Path,
    git_revision: str,
    context: str,
) -> dict[str, object]:
    """Bind live executable bytes to the retained qualification approval."""

    executable = Path(str(command.get("executable", "")))
    executable_record = authenticate_regular_file(
        executable,
        command.get("executable_sha256"),
        None,
        f"{context} executable",
        canonical_root / "build",
    )
    if not os.access(executable, os.X_OK):
        raise ValueError(f"{context} executable is not executable")
    qualification = command.get("qualification_approval")
    if not isinstance(qualification, dict):
        raise ValueError(f"{context} lacks qualification approval")
    expected_approval_path = (
        canonical_root / "accounting"
        / "mks24_stage_i_E03_forcing_policy_qualification_approval.json"
    )
    approval_path = Path(str(qualification.get("path", "")))
    if (
        not approval_path.is_absolute()
        or approval_path != expected_approval_path
        or approval_path.is_symlink()
    ):
        raise ValueError(f"{context} qualification approval path differs")
    approval, approval_digest = stable_json_file(
        approval_path, f"{context} qualification approval"
    )
    approved_executable = Path(str(approval.get("approved_executable", "")))
    if (
        qualification.get("sha256") != approval_digest
        or qualification.get("token") != approval
        or approval.get("schema_version") != 1
        or approval.get("execution_epoch") != FIGURE_13_EXECUTION_EPOCH
        or approval.get("approved_executable_revision") != git_revision
        or approval.get("approved_executable_sha256") != executable_record["sha256"]
        or approved_executable != executable
        or approved_executable.resolve() != Path(executable_record["path"])
        or qualification.get("execution_epoch") != FIGURE_13_EXECUTION_EPOCH
        or qualification.get("approved_executable_revision") != git_revision
        or qualification.get("approved_executable_sha256")
        != executable_record["sha256"]
        or command.get("build_manifest") != approval.get("build_manifest")
    ):
        raise ValueError(f"{context} qualification approval does not bind executable")
    restart_binary_abi = qualified_restart_binary_abi(
        git_revision, str(executable_record["sha256"]), context
    )
    restart_load_evidence = authenticate_qualified_restart_load_evidence(
        approval, canonical_root, executable_record, git_revision, context
    )
    return {
        "executable": executable_record,
        "qualification_approval_path": str(approval_path),
        "qualification_approval_sha256": approval_digest,
        "restart_binary_abi": restart_binary_abi,
        "restart_load_evidence": restart_load_evidence,
    }


def unavailable_execution_authentication(reason: str) -> dict[str, object]:
    """Return an explicit non-authorizing execution-authentication result."""

    return {
        "available": False,
        "reason": reason,
        "canonical_root": str(STAGE_I_CANONICAL_ROOT),
        "required_execution_epoch": FIGURE_13_EXECUTION_EPOCH,
        "required_final_time": FIGURE_13_ACCEPTED_FINAL_TIME,
        "snapshot_expected_ranks": {},
        "snapshot_controller_inventory": {},
        "controller_artifacts": {},
    }


def validate_controller_snapshot_inventory(
    artifact: dict[str, object], context: str, case_root: Path,
) -> dict[str, dict[str, object]]:
    """Return exact controller-inspected source rank sets."""

    allocation = artifact.get("allocation")
    inspection = artifact.get("scientific_inspection")
    if not isinstance(allocation, dict) or not isinstance(inspection, dict):
        raise ValueError(f"{context} lacks allocation or scientific inspection")
    try:
        expected_ranks = (
            int(allocation["nodes"]) * int(allocation["ranks_per_node"])
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{context} has invalid rank allocation") from error
    if expected_ranks <= 0:
        raise ValueError(f"{context} has nonpositive expected rank count")
    snapshots = inspection.get("snapshots")
    if not isinstance(snapshots, list):
        raise ValueError(f"{context} lacks inspected snapshot inventory")
    retained: dict[str, dict[str, object]] = {}
    for snapshot in snapshots:
        if not isinstance(snapshot, dict):
            raise ValueError(f"{context} has malformed inspected snapshot")
        rank_files = snapshot.get("rank_files")
        if not isinstance(rank_files, list) or len(rank_files) != expected_ranks:
            raise ValueError(
                f"{context} inspected snapshot lacks the exact expected rank count"
            )
        paths: list[Path] = []
        records: list[dict[str, object]] = []
        for record in rank_files:
            if not isinstance(record, dict):
                raise ValueError(f"{context} has malformed inspected rank file")
            path = Path(str(record.get("path", "")))
            digest = str(record.get("sha256", ""))
            size = record.get("size_bytes")
            if (
                not path.is_absolute()
                or path.is_symlink()
                or not path.is_file()
                or not path_is_beneath(path, case_root)
                or not re.fullmatch(r"[0-9a-f]{64}", digest)
                or not isinstance(size, int)
                or size <= 0
            ):
                raise ValueError(f"{context} has invalid inspected rank-file evidence")
            resolved = path.resolve()
            paths.append(resolved)
            records.append({
                "path": str(resolved),
                "sha256": digest,
                "size_bytes": size,
            })
        expected_names = [
            f"rank_{rank:08d}" for rank in range(expected_ranks)
        ]
        if [path.parent.name for path in paths] != expected_names:
            raise ValueError(
                f"{context} inspected snapshot ranks are not exactly contiguous"
            )
        if str(snapshot.get("path", "")) != str(paths[0]):
            raise ValueError(f"{context} inspected snapshot representative differs")
        if len({path.name for path in paths}) != 1:
            raise ValueError(f"{context} inspected rank files name different snapshots")
        source = str(paths[0])
        value = {
            "expected_rank_count": expected_ranks,
            "rank_files": records,
        }
        previous = retained.get(source)
        if previous is not None and previous != value:
            raise ValueError(f"{context} gives conflicting inspected snapshot evidence")
        retained[source] = value
    return retained


def controller_restart_contract(
    command: dict[str, object],
    run: dict[str, object],
    manifest_path: Path,
    case_id: str,
    case_name: str,
    input_sha256: str,
    context: str,
) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    """Authenticate and reconstruct the exact production restart-load contract."""

    input_path = Path(str(command.get("input_file", "")))
    expected_input_path = manifest_path.parent / "submitted_input.athinput"
    if input_path != expected_input_path:
        raise ValueError(f"{context} input_file does not name its submitted input")
    input_record = authenticate_regular_file(
        input_path, input_sha256, None, f"{context} submitted input", manifest_path.parent
    )
    input_payload = stable_regular_file_bytes(
        input_path, f"{context} submitted input"
    )
    if hashlib.sha256(input_payload).hexdigest() != input_record["sha256"]:
        raise ValueError(f"{context} submitted input changed during authentication")
    try:
        blocks = stage_i_validator.parse_athinput(
            input_payload, f"{context} submitted input",
        )
        case_contract = {
            "name": case_name,
            "input": FIGURE_13_CASE_CONTRACTS[case_id]["input"],
            "resolution": FIGURE_13_CASE_CONTRACTS[case_id]["resolution"],
        }
        if run.get("resolution") != case_contract["resolution"]:
            raise ValueError(f"{context} run resolution differs from its case contract")
        input_contract = stage_i_validator.require_input_contract(
            blocks, case_id, case_contract
        )
    except stage_i_validator.ValidationError as error:
        raise ValueError(f"{context} submitted input contract is invalid: {error}") from error
    try:
        target_time = float(command["time_tlim_target"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{context} lacks a finite runtime target") from error
    run_basename = run.get("run_basename")
    if (
        not isinstance(run_basename, str)
        or not run_basename
        or not math.isfinite(target_time)
        or target_time <= 0.0
    ):
        raise ValueError(f"{context} runtime contract is malformed")
    runtime_contract = {
        "run_basename": run_basename,
        "target_time": target_time,
    }
    retained = {
        "validator": stage_i_validator_implementation_provenance(),
        "submitted_input": input_record,
        "case_id": case_id,
        "case_name": case_name,
        "resolution": case_contract["resolution"],
        "runtime_contract": runtime_contract,
        "mesh": list(input_contract["mesh"]),
        "meshblock": list(input_contract["meshblock"]),
        "nghost": int(input_contract["nghost"]),
        "expected_meshblocks": int(input_contract["expected_meshblocks"]),
        "restart_data_size": int(input_contract["restart_data_size"]),
        "storage": json_compatible(input_contract["storage"]),
        "archived_parameter_blocks_sha256": canonical_json_sha256(
            input_contract["archived_parameter_blocks"]
        ),
        "turbulence_metadata_contract_sha256": canonical_json_sha256(
            {
                "ints": turbulence_restart_metadata_contract(input_contract)[1],
                "reals": turbulence_restart_metadata_contract(input_contract)[3],
            }
        ),
    }
    return input_contract, runtime_contract, retained


def retained_controller_restart_contract(
    record: dict[str, object], context: str,
) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    """Rebuild and compare a retained controller restart contract."""

    retained = record.get("restart_parameter_contract")
    if not isinstance(retained, dict):
        raise ValueError(f"{context} lacks retained restart parameter contract")
    submitted = retained.get("submitted_input")
    runtime = retained.get("runtime_contract")
    if not isinstance(submitted, dict) or not isinstance(runtime, dict):
        raise ValueError(f"{context} retained restart parameter contract is malformed")
    synthetic_command = {
        "input_file": submitted.get("path"),
        "time_tlim_target": runtime.get("target_time"),
    }
    synthetic_run = {
        "resolution": retained.get("resolution"),
        "run_basename": runtime.get("run_basename"),
    }
    manifest_path = Path(str(record.get("path", "")))
    input_contract, runtime_contract, current = controller_restart_contract(
        synthetic_command,
        synthetic_run,
        manifest_path,
        str(record.get("case_id", "")),
        str(record.get("case_name", "")),
        str(record.get("input_sha256", "")),
        context,
    )
    if current != retained:
        raise ValueError(f"{context} restart parameter contract changed")
    return input_contract, runtime_contract, current


def validate_controller_terminal_restart(
    artifact: dict[str, object],
    context: str,
    case_root: Path,
    final_time: float,
    restart_binary_abi: dict[str, object],
    input_contract: dict[str, object],
    runtime_contract: dict[str, object],
) -> dict[str, object]:
    """Return the exact authenticated terminal-restart identity for one segment."""

    allocation = artifact.get("allocation")
    inspection = artifact.get("scientific_inspection")
    if not isinstance(allocation, dict) or not isinstance(inspection, dict):
        raise ValueError(f"{context} lacks allocation or scientific inspection")
    try:
        expected_ranks = (
            int(allocation["nodes"]) * int(allocation["ranks_per_node"])
        )
        terminal_time = float(inspection["terminal_restart_time"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{context} has invalid terminal-restart metadata") from error
    terminal = inspection.get("terminal_restart")
    rank_files = terminal.get("rank_files") if isinstance(terminal, dict) else None
    if (
        expected_ranks <= 0
        or not figure_13_time_close(terminal_time, final_time)
        or not isinstance(rank_files, list)
        or len(rank_files) != expected_ranks
    ):
        raise ValueError(f"{context} lacks one exact terminal restart")
    restarts = inspection.get("restarts")
    restart_times = inspection.get("restart_times")
    if (
        not isinstance(restarts, list)
        or not isinstance(restart_times, list)
        or len(restarts) != len(restart_times)
        or len({
            str(record.get("path", ""))
            for record in restarts if isinstance(record, dict)
        }) != len(restarts)
    ):
        raise ValueError(f"{context} restart inventory is incomplete or ambiguous")
    try:
        parsed_restart_times = [float(value) for value in restart_times]
    except (TypeError, ValueError) as error:
        raise ValueError(f"{context} restart times are malformed") from error
    if not all(math.isfinite(value) for value in parsed_restart_times):
        raise ValueError(f"{context} restart times are nonfinite")
    terminal_matches = [
        index for index, (record, value) in enumerate(
            zip(restarts, parsed_restart_times)
        )
        if isinstance(record, dict)
        and figure_13_time_close(value, final_time)
        and record == terminal
    ]
    if len(terminal_matches) != 1:
        raise ValueError(f"{context} terminal restart is not unique")
    records, loadability = authenticate_restart_rank_files(
        rank_files,
        expected_ranks,
        case_root,
        final_time,
        restart_binary_abi,
        input_contract,
        runtime_contract,
        f"{context} terminal restart",
    )
    first = records[0]
    if (
        terminal.get("path") != first["path"]
        or terminal.get("sha256") != first["sha256"]
        or terminal.get("size_bytes") != first["size_bytes"]
    ):
        raise ValueError(f"{context} terminal restart identity is inconsistent")
    return {
        "time": terminal_time,
        "sha256": first["sha256"],
        "rank_files": records,
        "loadability": loadability,
    }


def validate_controller_artifact(
    path: Path,
    canonical_root: Path,
    case_id: str,
    case_name: str,
    input_sha256: str,
    git_revision: str,
    source_bundle_digests: dict[str, str],
) -> tuple[
    dict[str, object], dict[str, object], dict[str, dict[str, object]]
]:
    """Authenticate one canonical recorded Stage I controller artifact."""

    context = f"Figure 13 controller artifact {case_id}"
    expected_parent = (
        canonical_root / STAGE_I_RUNS_RELATIVE / FIGURE_13_EXECUTION_EPOCH
        / case_id
    )
    resolved = path.resolve()
    if (
        not path.is_absolute()
        or path.is_symlink()
        or not path_is_beneath(resolved, expected_parent)
        or resolved.name != "prepared_run.json"
        or resolved.parent.name != "manifest"
    ):
        raise ValueError(f"{context} path is outside its canonical case lineage")
    artifact, artifact_digest = stable_json_file(resolved, context)
    run = artifact.get("run")
    command = artifact.get("command")
    accounting = artifact.get("accounting")
    inspection = artifact.get("scientific_inspection")
    if not all(
        isinstance(value, dict)
        for value in (run, command, accounting, inspection)
    ):
        raise ValueError(f"{context} lacks recorded controller sections")
    segment = str(run.get("segment", ""))
    result = accounting.get("result")
    expected_manifest = str(resolved)
    if (
        artifact.get("project_root") != str(canonical_root)
        or artifact.get("execution_epoch") != FIGURE_13_EXECUTION_EPOCH
        or artifact.get("state") != "recorded"
        or run.get("case_id") != case_id
        or run.get("case_name") != case_name
        or resolved.parents[1].name != segment
        or accounting.get("case_id") != case_id
        or accounting.get("case_name") != case_name
        or accounting.get("segment") != segment
        or accounting.get("state") != "COMPLETED"
        or accounting.get("exit_code") != "0:0"
        or inspection.get("case_id") != case_id
        or inspection.get("segment") != segment
        or inspection.get("manifest") != expected_manifest
        or inspection.get("execution_epoch") != FIGURE_13_EXECUTION_EPOCH
        or result not in {"accepted", "clean_partial"}
        or (
            result == "accepted"
            and inspection.get("accepted") is not True
        )
        or (
            result == "clean_partial"
            and inspection.get("clean_for_continuation") is not True
        )
    ):
        raise ValueError(f"{context} is not a qualifying recorded execution")
    if (
        command.get("input_revision") != git_revision
        or accounting.get("input_revision") != git_revision
        or command.get("executable_revision") != git_revision
        or accounting.get("executable_revision") != git_revision
        or command.get("input_sha256") != input_sha256
        or command.get("matrix_sha256") != FIGURE_13_MATRIX_SHA256
    ):
        raise ValueError(f"{context} revision or input identity does not match")
    approved_executable = authenticate_approved_executable(
        command, canonical_root, git_revision, context
    )
    input_contract, runtime_contract, restart_parameter_contract = (
        controller_restart_contract(
            command,
            run,
            resolved,
            case_id,
            case_name,
            input_sha256,
            context,
        )
    )
    executable_digest = str(approved_executable["executable"]["sha256"])
    if (
        accounting.get("execution_epoch") != FIGURE_13_EXECUTION_EPOCH
        or accounting.get("executable_sha256") != executable_digest
    ):
        raise ValueError(f"{context} execution code identity does not match")
    source_bundle = command.get("source_bundle")
    if not isinstance(source_bundle, dict):
        raise ValueError(f"{context} lacks authenticated source-bundle evidence")
    source_path = Path(str(source_bundle.get("path", "")))
    source_digest = str(source_bundle.get("sha256", ""))
    revisions = source_bundle.get("verified_revisions")
    if (
        not source_path.is_absolute()
        or not path_is_beneath(source_path, canonical_root / "source-archives")
        or source_path.is_symlink()
        or not source_path.is_file()
        or not re.fullmatch(r"[0-9a-f]{64}", source_digest)
        or not isinstance(revisions, list)
        or git_revision not in revisions
    ):
        raise ValueError(f"{context} source-bundle evidence is invalid")
    actual_source_digest = source_bundle_digests.get(str(source_path.resolve()))
    if actual_source_digest is None:
        actual_source_digest = stable_regular_file_sha256(
            source_path, f"{context} source bundle"
        )
        source_bundle_digests[str(source_path.resolve())] = actual_source_digest
    if actual_source_digest != source_digest:
        raise ValueError(f"{context} source-bundle checksum does not match")
    try:
        final_time = float(inspection["final_time"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{context} has invalid accepted final time") from error
    if not math.isfinite(final_time) or final_time <= 0.0:
        raise ValueError(f"{context} has nonpositive accepted final time")
    snapshots = validate_controller_snapshot_inventory(
        artifact, context, expected_parent
    )
    terminal_restart = validate_controller_terminal_restart(
        artifact,
        context,
        expected_parent,
        final_time,
        approved_executable["restart_binary_abi"],
        input_contract,
        runtime_contract,
    )
    retained = {
        "path": str(resolved),
        "sha256": artifact_digest,
        "canonical_root": str(canonical_root.resolve()),
        "controller_case_root": str(expected_parent.resolve()),
        "case_id": case_id,
        "case_name": case_name,
        "segment": segment,
        "result": result,
        "final_time": final_time,
        "expected_rank_count": (
            int(artifact["allocation"]["nodes"])
            * int(artifact["allocation"]["ranks_per_node"])
        ),
        "input_revision": git_revision,
        "input_sha256": input_sha256,
        "restart_parameter_contract": restart_parameter_contract,
        "executable_revision": git_revision,
        "executable_sha256": executable_digest,
        "executable_path": approved_executable["executable"]["path"],
        "executable_size_bytes": approved_executable["executable"]["size_bytes"],
        "qualification_approval_path": approved_executable[
            "qualification_approval_path"
        ],
        "qualification_approval_sha256": approved_executable[
            "qualification_approval_sha256"
        ],
        "qualified_restart_load_evidence": approved_executable[
            "restart_load_evidence"
        ],
        "restart_binary_abi": approved_executable["restart_binary_abi"],
        "source_bundle_path": str(source_path.resolve()),
        "source_bundle_sha256": source_digest,
        "terminal_restart_sha256": terminal_restart["sha256"],
        "terminal_restart_rank_files": terminal_restart["rank_files"],
        "terminal_restart_loadability": terminal_restart["loadability"],
        "terminal_restart_files": [
            record["path"] for record in terminal_restart["rank_files"]
        ],
    }
    return artifact, retained, snapshots


def authenticate_child_restart_archive(
    command: dict[str, object],
    child_manifest: Path,
    parent_record: dict[str, object],
    context: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Bind a child's archived restart bytes and paths to its exact parent."""

    parent_rank_files = parent_record.get("terminal_restart_rank_files")
    if not isinstance(parent_rank_files, list) or not parent_rank_files:
        raise ValueError(f"{context} parent terminal restart evidence is malformed")
    input_contract, runtime_contract, _ = retained_controller_restart_contract(
        parent_record, f"{context} parent"
    )
    archive_root = child_manifest.parent / "submitted_restart"
    archive_records, loadability = authenticate_restart_rank_files(
        command.get("restart_files"),
        len(parent_rank_files),
        archive_root,
        float(parent_record["final_time"]),
        parent_record["restart_binary_abi"],
        input_contract,
        runtime_contract,
        f"{context} archived restart",
    )
    source_restart = Path(str(command.get("source_restart_file", "")))
    restart_file = Path(str(command.get("restart_file", "")))
    parent_paths = [str(record["path"]) for record in parent_rank_files]
    parent_identities = [
        (record["sha256"], record["size_bytes"]) for record in parent_rank_files
    ]
    archive_identities = [
        (record["sha256"], record["size_bytes"]) for record in archive_records
    ]
    if (
        not source_restart.is_absolute()
        or source_restart != Path(parent_paths[0])
        or not restart_file.is_absolute()
        or restart_file != Path(archive_records[0]["path"])
        or command.get("restart_sha256") != parent_record["terminal_restart_sha256"]
        or parent_identities != archive_identities
        or [Path(path).name for path in parent_paths]
        != [Path(str(record["path"])).name for record in archive_records]
    ):
        raise ValueError(f"{context} archived restart does not bind its parent")
    return archive_records, loadability


def validate_controller_case_lineage(
    case_id: str,
    artifacts: dict[Path, dict[str, object]],
    retained: dict[Path, dict[str, object]],
) -> None:
    """Require one complete authenticated lineage ending in accepted t=10."""

    terminals = [
        path for path, record in retained.items()
        if record["result"] == "accepted"
        and figure_13_time_close(
            float(record["final_time"]), FIGURE_13_ACCEPTED_FINAL_TIME
        )
    ]
    if len(terminals) != 1:
        raise ValueError(
            f"Figure 13 controller evidence for {case_id} lacks one accepted t=10 "
            "terminal"
        )
    visited: set[Path] = set()
    current = terminals[0]
    while True:
        if current in visited:
            raise ValueError(f"Figure 13 controller lineage for {case_id} contains a cycle")
        visited.add(current)
        artifact = artifacts[current]
        command = artifact["command"]
        parent = command.get("parent_segment")
        if parent is None:
            if (
                command.get("source_restart_file") is not None
                or command.get("restart_file") is not None
                or command.get("restart_sha256") is not None
                or command.get("restart_files", []) != []
            ):
                raise ValueError(
                    f"Figure 13 controller lineage for {case_id} has an "
                    "unauthenticated root restart"
                )
            break
        if not isinstance(parent, dict):
            raise ValueError(f"Figure 13 controller lineage for {case_id} is malformed")
        parent_path = Path(str(parent.get("manifest", ""))).resolve()
        parent_record = retained.get(parent_path)
        try:
            parent_final_time = float(parent["final_time"])
            parent_restart_time = float(parent["restart_time"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(
                f"Figure 13 controller lineage for {case_id} has an "
                "unauthenticated parent"
            ) from error
        child_restart_files = command.get("restart_files")
        parent_restart_rank_files = (
            parent_record.get("terminal_restart_rank_files")
            if isinstance(parent_record, dict) else None
        )
        child_restart_identity = (
            [
                (record.get("sha256"), record.get("size_bytes"))
                for record in child_restart_files
                if isinstance(record, dict)
            ]
            if isinstance(child_restart_files, list) else None
        )
        parent_restart_identity = (
            [
                (record.get("sha256"), record.get("size_bytes"))
                for record in parent_restart_rank_files
                if isinstance(record, dict)
            ]
            if isinstance(parent_restart_rank_files, list) else None
        )
        source_restart_file = Path(str(command.get("source_restart_file", "")))
        archive_records: list[dict[str, object]] | None = None
        archive_loadability: list[dict[str, object]] | None = None
        if parent_record is not None:
            try:
                archive_records, archive_loadability = authenticate_child_restart_archive(
                    command,
                    current,
                    parent_record,
                    f"Figure 13 controller lineage for {case_id}",
                )
            except ValueError as error:
                raise ValueError(
                    f"Figure 13 controller lineage for {case_id} has an "
                    "unauthenticated parent"
                ) from error
        if (
            parent_path not in artifacts
            or parent_record is None
            or parent.get("case_id") != case_id
            or parent.get("execution_epoch") != FIGURE_13_EXECUTION_EPOCH
            or parent.get("segment") != parent_record["segment"]
            or parent.get("result") != parent_record["result"]
            or parent.get("input_sha256") != parent_record["input_sha256"]
            or parent.get("executable_sha256")
            != parent_record["executable_sha256"]
            or parent.get("restart_sha256")
            != parent_record["terminal_restart_sha256"]
            or parent.get("restart_files")
            != parent_record["terminal_restart_files"]
            or not source_restart_file.is_absolute()
            or source_restart_file.resolve()
            != Path(parent_record["terminal_restart_files"][0]).resolve()
            or command.get("restart_sha256")
            != parent_record["terminal_restart_sha256"]
            or child_restart_identity != parent_restart_identity
            or not figure_13_time_close(
                parent_final_time, float(parent_record["final_time"])
            )
            or not figure_13_time_close(
                parent_restart_time, float(parent_record["final_time"])
            )
            or float(retained[current]["final_time"])
            <= float(parent_record["final_time"]) + FIGURE_13_TIME_TOLERANCE
        ):
            raise ValueError(
                f"Figure 13 controller lineage for {case_id} has an unauthenticated parent"
            )
        retained[current]["restart_archive_rank_files"] = archive_records
        retained[current]["restart_archive_loadability"] = archive_loadability
        retained[current]["restart_archive_time"] = float(parent_record["final_time"])
        retained[current]["restart_archive_parent_manifest"] = parent_record["path"]
        retained[current]["restart_archive_parameter_contract"] = parent_record[
            "restart_parameter_contract"
        ]
        current = parent_path
    if visited != set(artifacts):
        raise ValueError(
            f"Figure 13 controller evidence for {case_id} contains unlinked segments"
        )


def figure_13_execution_authentication(
    bundle: Path, manifest: dict[str, object], cases: list[dict[str, object]],
) -> dict[str, object]:
    """Authenticate bundle execution identity against canonical controller artifacts."""

    canonical_root = STAGE_I_CANONICAL_ROOT.resolve()
    bundle_root = (
        canonical_root / STAGE_I_RUNS_RELATIVE / FIGURE_13_EXECUTION_EPOCH
        / "bundles"
    )
    if not path_is_beneath(bundle, bundle_root):
        return unavailable_execution_authentication(
            "bundle is outside the canonical retained Stage I bundle namespace"
        )
    segment_values = manifest.get("production_segment_manifests")
    if not isinstance(segment_values, list) or not segment_values:
        return unavailable_execution_authentication(
            "bundle lacks retained production controller artifacts"
        )
    git_revision = manifest.get("git_revision")
    if (
        manifest.get("status") != "accepted_for_analysis"
        or manifest.get("execution_epoch") != FIGURE_13_EXECUTION_EPOCH
        or not isinstance(git_revision, str)
        or re.fullmatch(r"[0-9a-f]{40}", git_revision) is None
    ):
        raise ValueError(
            "Figure 13 bundle execution epoch, revision, or accepted status is invalid"
        )
    selected: dict[str, dict[str, object]] = {}
    for case in cases:
        case_id = figure_13_case_id(case, manifest, len(cases))
        if case_id is not None:
            selected[case_id] = case
    if not selected:
        return unavailable_execution_authentication(
            "bundle contains no Figure 13 cases to authenticate"
        )
    paths = [Path(str(value)) for value in segment_values]
    if len({str(path.resolve()) for path in paths}) != len(paths):
        raise ValueError("Figure 13 bundle repeats production controller artifacts")
    source_bundle_digests: dict[str, str] = {}
    artifacts_by_case: dict[str, dict[Path, dict[str, object]]] = {
        case_id: {} for case_id in selected
    }
    retained_by_case: dict[str, dict[Path, dict[str, object]]] = {
        case_id: {} for case_id in selected
    }
    source_snapshots_by_case: dict[str, dict[str, dict[str, object]]] = {
        case_id: {} for case_id in selected
    }
    for path in paths:
        try:
            relative = path.resolve().relative_to(
                canonical_root / STAGE_I_RUNS_RELATIVE / FIGURE_13_EXECUTION_EPOCH
            )
        except ValueError:
            raise ValueError(
                "Figure 13 bundle names a controller artifact outside canonical Stage I"
            ) from None
        if len(relative.parts) < 4:
            raise ValueError("Figure 13 bundle names a malformed controller artifact")
        case_id = relative.parts[0]
        if case_id not in selected:
            continue
        case = selected[case_id]
        contract = FIGURE_13_CASE_CONTRACTS[case_id]
        artifact, retained, snapshots = validate_controller_artifact(
            path,
            canonical_root,
            case_id,
            str(contract["name"]),
            str(contract["input_sha256"]),
            git_revision,
            source_bundle_digests,
        )
        resolved = path.resolve()
        artifacts_by_case[case_id][resolved] = artifact
        retained_by_case[case_id][resolved] = retained
        for source, snapshot in snapshots.items():
            previous = source_snapshots_by_case[case_id].get(source)
            if previous is not None and previous != snapshot:
                raise ValueError(
                    f"Figure 13 controller artifacts disagree on {case_id} snapshots"
                )
            source_snapshots_by_case[case_id][source] = snapshot
    snapshot_expected_ranks: dict[str, int] = {}
    snapshot_controller_inventory: dict[str, dict[str, object]] = {}
    retained_output: dict[str, list[dict[str, object]]] = {}
    for case_id, case in selected.items():
        if not artifacts_by_case[case_id]:
            return unavailable_execution_authentication(
                f"bundle lacks canonical controller artifacts for {case_id}"
            )
        validate_controller_case_lineage(
            case_id, artifacts_by_case[case_id], retained_by_case[case_id]
        )
        outputs = case.get("outputs")
        declared = outputs.get("snapshot_paths") if isinstance(outputs, dict) else None
        if not isinstance(declared, list) or not declared:
            raise ValueError(f"Figure 13 bundle case {case_id} lacks snapshot paths")
        for value in declared:
            bundle_path = bundle_relative_member_path(
                bundle,
                str(value),
                f"Figure 13 bundle case {case_id} snapshot",
            )
            if not bundle_path.is_file():
                raise ValueError(
                    f"Figure 13 bundle case {case_id} snapshot is missing: {bundle_path}"
                )
            source = str(bundle_path.resolve())
            inspected = source_snapshots_by_case[case_id].get(source)
            if inspected is None:
                raise ValueError(
                    f"Figure 13 bundle case {case_id} snapshot is not retained by its "
                    "authenticated controller lineage"
                )
            expected_ranks = int(inspected["expected_rank_count"])
            bundle_rank_files = snapshot_sibling_paths(bundle_path, expected_ranks)
            inspected_rank_files = inspected["rank_files"]
            if (
                not isinstance(inspected_rank_files, list)
                or [str(path.resolve()) for path in bundle_rank_files]
                != [
                    str(record.get("path", ""))
                    for record in inspected_rank_files
                    if isinstance(record, dict)
                ]
            ):
                raise ValueError(
                    f"Figure 13 bundle case {case_id} rank siblings do not match "
                    "the authenticated controller inventory"
                )
            snapshot_expected_ranks[str(bundle_path)] = expected_ranks
            snapshot_controller_inventory[str(bundle_path)] = inspected
        retained_output[case_id] = [
            retained_by_case[case_id][path]
            for path in sorted(retained_by_case[case_id], key=str)
        ]
    return {
        "available": True,
        "canonical_root": str(canonical_root),
        "required_execution_epoch": FIGURE_13_EXECUTION_EPOCH,
        "required_final_time": FIGURE_13_ACCEPTED_FINAL_TIME,
        "authenticated_execution_epoch": manifest["execution_epoch"],
        "authenticated_git_revision": git_revision,
        "authenticated_bundle_status": manifest["status"],
        "snapshot_expected_ranks": snapshot_expected_ranks,
        "snapshot_expected_ranks_sha256": canonical_json_sha256(
            snapshot_expected_ranks
        ),
        "snapshot_controller_inventory": snapshot_controller_inventory,
        "snapshot_controller_inventory_sha256": canonical_json_sha256(
            snapshot_controller_inventory
        ),
        "controller_artifacts": retained_output,
        "controller_artifacts_sha256": canonical_json_sha256(retained_output),
    }


def require_manifest_text(record: dict[str, object], key: str, context: str) -> str:
    """Return one required nonempty manifest text value."""

    value = record.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context} requires nonempty {key}")
    return value


def figure_13_case_id(case: dict[str, object], manifest: dict[str, object],
                      case_count: int) -> str | None:
    """Resolve one Figure 13 case without silently accepting conflicting IDs."""

    name = str(case.get("name", ""))
    candidates: list[str] = []
    explicit = case.get("case_id")
    if explicit is not None:
        if not isinstance(explicit, str) or not explicit:
            raise ValueError(f"bundle case has invalid case_id: {name}")
        candidates.append(explicit)
    if case_count == 1 and manifest.get("production_case_id") is not None:
        production = manifest["production_case_id"]
        if not isinstance(production, str) or not production:
            raise ValueError("bundle has invalid production_case_id")
        candidates.append(production)
    if name in FIGURE_13_FIREHOSE_CASE_NAMES:
        candidates.append(FIGURE_13_FIREHOSE_CASE_NAMES[name])
    if len(set(candidates)) > 1:
        raise ValueError(f"bundle case has conflicting Figure 13 identities: {name}")
    case_id = candidates[0] if candidates else None
    return case_id if case_id in FIGURE_13_FIREHOSE_CASES else None


def figure_13_case_provenance(
    bundle: Path,
    case_id: str,
    case: dict[str, object],
    matrix_provenance: dict[str, object],
) -> dict[str, object]:
    """Validate one Figure 13 case against its exact tracked physics identity."""

    contract = FIGURE_13_CASE_CONTRACTS[case_id]
    context = f"Figure 13 case {case_id}"
    name = require_manifest_text(case, "name", context)
    input_name = require_manifest_text(case, "input", context)
    if name != contract["name"] or input_name != contract["input"]:
        raise ValueError(f"{context} does not match its exact matrix identity")
    if case.get("status") != "passed" or case.get("lf_active") is not True or (
        case.get("amr") is not False
    ):
        raise ValueError(f"{context} does not have exact accepted LF case metadata")
    model = case.get("model_choices")
    if not isinstance(model, dict):
        raise ValueError(f"{context} requires model_choices")
    model_digest = canonical_json_sha256(model)
    if model_digest != contract["model_choices_sha256"]:
        raise ValueError(f"{context} model choices checksum does not match")
    execution_input_name = require_manifest_text(case, "execution_input", context)
    execution_input = bundle_relative_member_path(
        bundle, execution_input_name, f"{context} execution input"
    )
    if (
        not execution_input.is_file()
        or execution_input.is_symlink()
        or not path_is_beneath(execution_input, bundle)
    ):
        raise ValueError(f"{context} execution input is not a regular file")
    execution_digest = sha256_file(execution_input)
    if execution_digest != contract["input_sha256"]:
        raise ValueError(f"{context} execution input checksum does not match")
    analysis_start = optional_float(model["analysis_t_start"])
    analysis_end = optional_float(model["analysis_t_end"])
    if analysis_start != FIGURE_13_WINDOW_START or analysis_end != FIGURE_13_WINDOW_END:
        raise ValueError(f"{context} does not archive the exact t=8..10 window")
    if str(model["output2_dt"]) != str(FIGURE_13_SNAPSHOT_CADENCE):
        raise ValueError(f"{context} does not archive the exact 0.25 snapshot cadence")
    source = matrix_provenance["cases"][case_id]
    return {
        "case_id": case_id,
        "case_name": name,
        "input": input_name,
        "input_sha256": contract["input_sha256"],
        "execution_input": execution_input_name,
        "execution_input_path": str(execution_input),
        "execution_input_sha256": execution_digest,
        "matrix_identity": source["matrix_identity"],
        "source_input_path": source["source_input_path"],
        "source_input_sha256": source["source_input_sha256"],
        "model_choices_sha256": model_digest,
        "accepted_case_metadata": {
            "status": case["status"],
            "lf_active": case["lf_active"],
            "amr": case["amr"],
        },
        "physics_identity": {
            "beta0": model["beta0"],
            "resolution": (
                f"{model['mesh_nx1']}x{model['mesh_nx2']}x{model['mesh_nx3']}"
            ),
            "forcing_mode": model["forcing_mode"],
            "forcing_seed": model["forcing_seed"],
            "forcing_tcorr": model["forcing_tcorr"],
            "forcing_dedt": model["forcing_dedt"],
            "passive_delta": model["passive_delta"],
            "cgl_firehose_threshold": model["cgl_firehose_threshold"],
            "mirror_limiter": model["mirror_limiter"],
            "firehose_limiter": model["firehose_limiter"],
            "limiter_hardwall": model["limiter_hardwall"],
            "limiter_nu_coll": model["limiter_nu_coll"],
            "output2_dt": model["output2_dt"],
            "output2_single_file_per_rank": model[
                "output2_single_file_per_rank"
            ],
        },
        "archived_analysis_window": {
            "time_start": analysis_start,
            "time_end": analysis_end,
            "snapshot_cadence": FIGURE_13_SNAPSHOT_CADENCE,
        },
    }


def validate_snapshot_provenance_record(
    provenance: object, context: str,
) -> dict[str, object]:
    """Validate one deterministic logical-snapshot provenance record."""

    if not isinstance(provenance, dict):
        raise ValueError(f"{context} snapshot provenance is malformed")
    files = provenance.get("files")
    if not isinstance(files, list) or not files:
        raise ValueError(f"{context} snapshot provenance lacks files")
    layout = provenance.get("layout")
    expected_rank_count = provenance.get("expected_rank_count")
    rank_directory_names = provenance.get("rank_directory_names")
    if (
        layout not in {"rank_local_siblings", "single_file"}
        or type(expected_rank_count) is not int
        or expected_rank_count <= 0
        or not isinstance(rank_directory_names, list)
        or len(files) != expected_rank_count
    ):
        raise ValueError(f"{context} snapshot provenance lacks exact rank-set identity")
    snapshot_time_value = provenance.get("snapshot_time")
    if type(snapshot_time_value) not in (int, float) or not math.isfinite(
        float(snapshot_time_value)
    ):
        raise ValueError(f"{context} snapshot provenance lacks finite time")
    for record in files:
        if not isinstance(record, dict):
            raise ValueError(f"{context} snapshot provenance file is malformed")
        if not isinstance(record.get("path"), str) or not record["path"]:
            raise ValueError(f"{context} snapshot provenance file lacks path")
        if "symlink_target" not in record or (
            record["symlink_target"] is not None
            and not isinstance(record["symlink_target"], str)
        ):
            raise ValueError(
                f"{context} snapshot provenance file has invalid symlink target"
            )
        if not isinstance(record.get("size_bytes"), int) or record["size_bytes"] < 0:
            raise ValueError(f"{context} snapshot provenance file has invalid size")
        if not re.fullmatch(r"[0-9a-f]{64}", str(record.get("sha256", ""))):
            raise ValueError(f"{context} snapshot provenance file has invalid checksum")
    representative = provenance.get("representative_path")
    if not isinstance(representative, str) or representative != files[0]["path"]:
        raise ValueError(f"{context} snapshot provenance representative differs")
    if layout == "rank_local_siblings":
        expected_names = [
            f"rank_{rank:08d}" for rank in range(expected_rank_count)
        ]
        actual_names = [Path(str(record["path"])).parent.name for record in files]
        if rank_directory_names != expected_names or actual_names != expected_names:
            raise ValueError(
                f"{context} snapshot provenance rank set is not exactly contiguous"
            )
    elif expected_rank_count != 1 or rank_directory_names:
        raise ValueError(f"{context} shared snapshot provenance rank identity is invalid")
    if provenance.get("aggregate_sha256") != canonical_json_sha256(files):
        raise ValueError(f"{context} snapshot provenance checksum does not match")
    return provenance


def revalidate_snapshot_provenance_record(
    provenance: dict[str, object], expected_ranks: int | None, context: str,
) -> None:
    """Rehash one bound logical snapshot and require exact byte identity."""

    validated = validate_snapshot_provenance_record(provenance, context)
    representative = Path(str(validated["representative_path"]))
    try:
        current = snapshot_digest_provenance(representative, expected_ranks)
    except (OSError, ValueError) as error:
        raise ValueError(
            f"{context} snapshot bytes changed before publication"
        ) from error
    retained_time = float(validated["snapshot_time"])
    current_time = float(current["snapshot_time"])
    retained_bytes = {
        key: value for key, value in validated.items() if key != "snapshot_time"
    }
    current_bytes = {
        key: value for key, value in current.items() if key != "snapshot_time"
    }
    if (
        current_bytes != retained_bytes
        or not figure_13_time_close(current_time, retained_time)
    ):
        raise ValueError(f"{context} snapshot bytes changed before publication")


def revalidate_figure_13_snapshot_bytes(
    cases: dict[str, dict[str, object]],
    execution_authentication: dict[str, object],
) -> None:
    """Revalidate every bound snapshot immediately before comparison use."""

    expected_by_path = execution_authentication.get("snapshot_expected_ranks")
    if not isinstance(expected_by_path, dict):
        raise ValueError("Figure 13 execution authentication lacks expected rank sets")
    revalidated: set[str] = set()
    for case_id, case in cases.items():
        if case.get("available") is not True:
            continue
        records = case.get("snapshot_provenance")
        if not isinstance(records, list):
            raise ValueError(f"Figure 13 case {case_id} lacks bound snapshots")
        for record in records:
            if not isinstance(record, dict):
                raise ValueError(f"Figure 13 case {case_id} snapshot evidence is malformed")
            representative = str(record.get("representative_path", ""))
            expected_ranks = expected_by_path.get(representative)
            if not isinstance(expected_ranks, int) or expected_ranks <= 0:
                raise ValueError(
                    f"Figure 13 case {case_id} snapshot lacks authenticated expected ranks"
                )
            aggregate = str(record.get("aggregate_sha256", ""))
            if aggregate in revalidated:
                continue
            revalidate_snapshot_provenance_record(
                record, expected_ranks, f"Figure 13 case {case_id}"
            )
            revalidated.add(aggregate)


def figure_13_semantics_measurement(
    case_id: str,
    semantics: str,
    measured: dict[str, object],
    snapshot_times: list[float],
) -> dict[str, object]:
    """Validate and extract one strict/inclusive Figure 13 measurement."""

    components: dict[str, float] = {}
    equal_snapshot_means: dict[str, float] = {}
    series: dict[str, list[float]] = {}
    for name in INSTABILITY_OCCUPANCY_COMPONENTS:
        record = measured.get(name)
        required = (
            "comparison_volume_fraction",
            "time_average_volume_fraction",
            "equal_snapshot_mean_volume_fraction",
            "snapshot_volume_fractions",
            "snapshot_fraction_minimum",
            "snapshot_fraction_maximum",
        )
        if not isinstance(record, dict) or any(key not in record for key in required):
            raise ValueError(
                f"Figure 13 case {case_id} lacks {semantics} measurement {name}"
            )
        comparison = float(record["comparison_volume_fraction"])
        time_average = float(record["time_average_volume_fraction"])
        equal_snapshot = float(record["equal_snapshot_mean_volume_fraction"])
        values = [float(value) for value in record["snapshot_volume_fractions"]]
        minimum = float(record["snapshot_fraction_minimum"])
        maximum = float(record["snapshot_fraction_maximum"])
        if not all(
            math.isfinite(value) and 0.0 <= value <= 1.0
            for value in (
                comparison,
                time_average,
                equal_snapshot,
                minimum,
                maximum,
                *values,
            )
        ):
            raise ValueError(
                f"Figure 13 case {case_id} has invalid {semantics} measurement {name}"
            )
        if len(values) != len(snapshot_times):
            raise ValueError(
                f"Figure 13 case {case_id} {semantics} series lacks exact cadence"
            )
        expected_time_average = (
            trapezoidal_integral(values, snapshot_times)
            / float(snapshot_times[-1] - snapshot_times[0])
            if len(snapshot_times) > 1
            else values[0]
        )
        expected_equal_snapshot = float(np.mean(values))
        if not all((
            math.isclose(
                comparison, expected_time_average, rel_tol=0.0, abs_tol=1.0e-15
            ),
            math.isclose(
                time_average, expected_time_average, rel_tol=0.0, abs_tol=1.0e-15
            ),
            math.isclose(
                equal_snapshot,
                expected_equal_snapshot,
                rel_tol=0.0,
                abs_tol=1.0e-15,
            ),
            math.isclose(minimum, min(values), rel_tol=0.0, abs_tol=1.0e-15),
            math.isclose(maximum, max(values), rel_tol=0.0, abs_tol=1.0e-15),
        )):
            raise ValueError(
                f"Figure 13 case {case_id} {semantics} measurement {name} "
                "does not match its snapshot series"
            )
        if not math.isclose(
            comparison, time_average, rel_tol=0.0, abs_tol=1.0e-15
        ):
            raise ValueError(
                f"Figure 13 case {case_id} does not use time average as primary"
            )
        components[name] = comparison
        equal_snapshot_means[name] = equal_snapshot
        series[name] = values
    if len({len(values) for values in series.values()}) != 1:
        raise ValueError(
            f"Figure 13 case {case_id} {semantics} series lengths are inconsistent"
        )

    def validate_relationships(values: dict[str, object], label: str) -> None:
        for parallel_name, oblique_name in (
            ("parallel_firehose", "oblique_firehose"),
            ("parallel_total_unstable", "oblique_total_unstable"),
        ):
            if not np.allclose(
                np.asarray(values[oblique_name], dtype=float)
                - np.asarray(values[parallel_name], dtype=float),
                np.asarray(values["oblique_only_reclassification"], dtype=float),
                rtol=0.0,
                atol=1.0e-15,
            ):
                raise ValueError(
                    f"Figure 13 case {case_id} {semantics} {label} masks "
                    "are inconsistent"
                )
        for total_name, firehose_name in (
            ("parallel_total_unstable", "parallel_firehose"),
            ("oblique_total_unstable", "oblique_firehose"),
        ):
            if not np.allclose(
                np.asarray(values[total_name], dtype=float),
                np.asarray(values["mirror"], dtype=float)
                + np.asarray(values[firehose_name], dtype=float),
                rtol=0.0,
                atol=1.0e-15,
            ):
                raise ValueError(
                    f"Figure 13 case {case_id} {semantics} {label} total is "
                    f"not mirror OR {firehose_name}"
                )

    validate_relationships(components, "time-average")
    validate_relationships(equal_snapshot_means, "equal-snapshot")
    validate_relationships(series, "snapshot-series")
    return {
        "components": components,
        "equal_snapshot_mean_components": equal_snapshot_means,
        "snapshot_fraction_series": series,
    }


def figure_13_case_occupancy(
    bundle: Path,
    case_id: str,
    case: dict[str, object],
    analyzed_case: dict[str, object],
    matrix_provenance: dict[str, object],
    analysis_configuration: dict[str, object],
    execution_authentication: dict[str, object],
) -> dict[str, object]:
    """Build one fail-closed, comparison-ready Figure 13 case record."""

    provenance = figure_13_case_provenance(
        bundle, case_id, case, matrix_provenance
    )
    ensemble = analyzed_case.get("snapshot_ensemble")
    if not isinstance(ensemble, dict):
        raise ValueError(f"Figure 13 case {case_id} snapshot ensemble is malformed")
    occupancy = ensemble.get("firehose_threshold_occupancy")
    if not isinstance(occupancy, dict):
        return {
            "available": False,
            "reason": "no retained snapshots selected in the requested late-time window",
            "case_provenance": provenance,
        }
    normalization = occupancy.get("normalization")
    temporal_normalization = (
        normalization.get("temporal") if isinstance(normalization, dict) else None
    )
    if (
        occupancy.get("schema_version") != FIGURE_13_FIREHOSE_SCHEMA_VERSION
        or occupancy.get("threshold_definitions") != instability_threshold_definitions()
        or not isinstance(temporal_normalization, dict)
        or temporal_normalization.get("primary_comparison_kind")
        != "trapezoidal_time_average"
        or temporal_normalization.get("secondary_kind")
        != "equal_weight_retained_snapshot_mean"
    ):
        raise ValueError(
            f"Figure 13 case {case_id} occupancy definitions do not match"
        )
    window = occupancy.get("analysis_window")
    grid = occupancy.get("grid")
    measured = occupancy.get("occupancy")
    if not isinstance(window, dict) or not isinstance(grid, dict) or not isinstance(
        measured, dict
    ):
        raise ValueError(f"Figure 13 case {case_id} occupancy product is malformed")
    required_window = (
        "requested_time_start",
        "requested_time_end",
        "selected_time_first",
        "selected_time_last",
        "snapshot_count",
        "snapshot_times",
    )
    if any(name not in window for name in required_window):
        raise ValueError(f"Figure 13 case {case_id} occupancy window is incomplete")
    requested_start = optional_float(window["requested_time_start"])
    requested_end = optional_float(window["requested_time_end"])
    selected_start = float(window["selected_time_first"])
    selected_end = float(window["selected_time_last"])
    snapshot_count = int(window["snapshot_count"])
    snapshot_times = [float(value) for value in window["snapshot_times"]]
    if not all(math.isfinite(value) for value in (
        selected_start, selected_end, *snapshot_times
    )):
        raise ValueError(f"Figure 13 case {case_id} occupancy window is nonfinite")
    try:
        cadence_validation = validate_figure_13_snapshot_times(snapshot_times)
    except ValueError:
        cadence_validation = None
    exact_window = bool(
        requested_start is not None
        and requested_end is not None
        and figure_13_time_close(requested_start, FIGURE_13_WINDOW_START)
        and figure_13_time_close(requested_end, FIGURE_13_WINDOW_END)
        and snapshot_count == len(FIGURE_13_SNAPSHOT_TIMES)
        and snapshot_count == len(snapshot_times)
        and cadence_validation is not None
        and figure_13_time_close(selected_start, snapshot_times[0])
        and figure_13_time_close(selected_end, snapshot_times[-1])
    )
    analysis_window = {
        "requested_time_start": requested_start,
        "requested_time_end": requested_end,
        "selected_time_first": selected_start,
        "selected_time_last": selected_end,
        "snapshot_count": snapshot_count,
        "snapshot_times": snapshot_times,
        "cadence_validation": cadence_validation,
    }
    if analysis_configuration.get("time_start") is not None or (
        analysis_configuration.get("time_end") is not None
    ):
        return {
            "available": False,
            "reason": (
                "comparison-ready Figure 13 status forbids CLI time-window overrides"
            ),
            "case_provenance": provenance,
            "analysis_window": analysis_window,
            "grid": grid,
        }
    if not exact_window:
        return {
            "available": False,
            "reason": (
                "retained snapshots do not exactly cover t=8..10 at complete "
                "0.25 cadence"
            ),
            "case_provenance": provenance,
            "analysis_window": analysis_window,
            "grid": grid,
        }
    try:
        shape = tuple(int(value) for value in grid["shape_z_y_x"])
        grid_lengths = tuple(float(value) for value in grid["lengths_x_y_z"])
        spatial_normalization = occupancy["normalization"]["spatial"]
        total_cell_count = int(spatial_normalization["total_cell_count"])
        physical_volume = float(spatial_normalization["physical_volume"])
        cell_volume = float(spatial_normalization["cell_volume"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            f"Figure 13 case {case_id} grid normalization is malformed"
        ) from error
    expected_cell_count = math.prod(FIGURE_13_GRID_SHAPE_Z_Y_X)
    expected_physical_volume = math.prod(FIGURE_13_GRID_LENGTHS_X_Y_Z)
    exact_grid = bool(
        shape == FIGURE_13_GRID_SHAPE_Z_Y_X
        and len(grid_lengths) == len(FIGURE_13_GRID_LENGTHS_X_Y_Z)
        and all(
            math.isclose(actual, expected, rel_tol=0.0, abs_tol=1.0e-12)
            for actual, expected in zip(
                grid_lengths, FIGURE_13_GRID_LENGTHS_X_Y_Z
            )
        )
        and total_cell_count == expected_cell_count
        and math.isclose(
            physical_volume, expected_physical_volume, rel_tol=0.0, abs_tol=1.0e-12
        )
        and math.isclose(
            cell_volume,
            expected_physical_volume / expected_cell_count,
            rel_tol=0.0,
            abs_tol=1.0e-24,
        )
    )
    if not exact_grid:
        return {
            "available": False,
            "reason": (
                "retained snapshots do not use the exact 192x192x384 "
                "uniform Cartesian Figure 13 grid"
            ),
            "case_provenance": provenance,
            "analysis_window": analysis_window,
            "grid": grid,
        }
    snapshot_provenance = occupancy.get("snapshot_provenance")
    if not isinstance(snapshot_provenance, list) or len(
        snapshot_provenance
    ) != len(FIGURE_13_SNAPSHOT_TIMES):
        return {
            "available": False,
            "reason": "complete per-snapshot and rank-sibling provenance is required",
            "case_provenance": provenance,
            "analysis_window": analysis_window,
            "grid": grid,
        }
    validated_snapshot_provenance = [
        validate_snapshot_provenance_record(value, f"Figure 13 case {case_id}")
        for value in snapshot_provenance
    ]
    if not all(
        figure_13_time_close(float(record["snapshot_time"]), expected)
        for record, expected in zip(
            validated_snapshot_provenance, snapshot_times
        )
    ):
        raise ValueError(
            f"Figure 13 case {case_id} snapshot provenance times do not match"
        )
    outputs = case.get("outputs")
    declared_paths = (
        outputs.get("snapshot_paths") if isinstance(outputs, dict) else None
    )
    if not isinstance(declared_paths, list) or not all(
        isinstance(path, str) and path for path in declared_paths
    ):
        raise ValueError(
            f"Figure 13 case {case_id} lacks declared snapshot-path provenance"
        )
    declared = {str(bundle / path) for path in declared_paths}
    representatives = [
        str(record.get("representative_path", ""))
        for record in validated_snapshot_provenance
    ]
    if len(set(representatives)) != len(representatives) or any(
        representative not in declared for representative in representatives
    ):
        raise ValueError(
            f"Figure 13 case {case_id} analyzed snapshots do not match its manifest"
        )
    expected_provenance_digest = canonical_json_sha256(
        validated_snapshot_provenance
    )
    if occupancy.get("snapshot_provenance_sha256") != expected_provenance_digest:
        raise ValueError(
            f"Figure 13 case {case_id} snapshot provenance checksum does not match"
        )
    execution_authenticated = execution_authentication.get("available") is True
    expected_by_path = execution_authentication.get("snapshot_expected_ranks", {})
    controller_inventory = execution_authentication.get(
        "snapshot_controller_inventory", {}
    )
    if execution_authenticated:
        if not isinstance(expected_by_path, dict) or not isinstance(
            controller_inventory, dict
        ):
            raise ValueError(
                f"Figure 13 case {case_id} lacks authenticated snapshot rank sets"
            )
        for record in validated_snapshot_provenance:
            representative = str(record["representative_path"])
            expected_ranks = expected_by_path.get(representative)
            inspected = controller_inventory.get(representative)
            if (
                not isinstance(expected_ranks, int)
                or record["expected_rank_count"] != expected_ranks
                or not isinstance(inspected, dict)
                or inspected.get("expected_rank_count") != expected_ranks
            ):
                raise ValueError(
                    f"Figure 13 case {case_id} snapshot rank set does not match "
                    "authenticated controller allocation"
                )
            inspected_files = inspected.get("rank_files")
            if not isinstance(inspected_files, list) or len(
                inspected_files
            ) != expected_ranks:
                raise ValueError(
                    f"Figure 13 case {case_id} lacks authenticated rank-file evidence"
                )
            for observed, expected in zip(record["files"], inspected_files):
                if (
                    not isinstance(expected, dict)
                    or Path(str(observed["path"])).resolve()
                    != Path(str(expected.get("path", "")))
                    or observed["size_bytes"] != expected.get("size_bytes")
                    or observed["sha256"] != expected.get("sha256")
                ):
                    raise ValueError(
                        f"Figure 13 case {case_id} analyzed rank bytes do not match "
                        "authenticated controller inventory"
                    )
    semantics_measurements: dict[str, object] = {}
    for semantics in INSTABILITY_OCCUPANCY_SEMANTICS:
        semantic_product = measured.get(semantics)
        if not isinstance(semantic_product, dict):
            raise ValueError(
                f"Figure 13 case {case_id} lacks occupancy semantics {semantics}"
            )
        semantics_measurements[semantics] = figure_13_semantics_measurement(
            case_id, semantics, semantic_product, snapshot_times
        )
        if any(
            len(values) != len(FIGURE_13_SNAPSHOT_TIMES)
            for values in semantics_measurements[semantics][
                "snapshot_fraction_series"
            ].values()
        ):
            raise ValueError(
                f"Figure 13 case {case_id} {semantics} series lacks exact cadence"
            )
    strict = semantics_measurements["strict_manuscript"]["components"]
    return {
        "available": True,
        "comparison_ready": execution_authenticated,
        "execution_authenticated": execution_authenticated,
        "case_provenance": provenance,
        "analysis_window": analysis_window,
        "grid": grid,
        "normalization": occupancy["normalization"],
        "primary_semantics": "strict_manuscript",
        "primary_temporal_aggregation": "trapezoidal_time_average",
        "parallel_total_unstable_volume_fraction": strict[
            "parallel_total_unstable"
        ],
        "oblique_total_unstable_volume_fraction": strict["oblique_total_unstable"],
        "mirror_volume_fraction": strict["mirror"],
        "parallel_firehose_volume_fraction": strict["parallel_firehose"],
        "oblique_firehose_volume_fraction": strict["oblique_firehose"],
        "oblique_only_reclassified_volume_fraction": strict[
            "oblique_only_reclassification"
        ],
        "parallel_total_unstable_percent": (
            100.0 * strict["parallel_total_unstable"]
        ),
        "oblique_total_unstable_percent": (
            100.0 * strict["oblique_total_unstable"]
        ),
        "oblique_only_reclassified_percentage_points": (
            100.0 * strict["oblique_only_reclassification"]
        ),
        "semantics": semantics_measurements,
        "snapshot_provenance": validated_snapshot_provenance,
        "snapshot_provenance_sha256": expected_provenance_digest,
        "measurement_provenance_sha256": canonical_json_sha256({
            "case_provenance": provenance,
            "analysis_window": analysis_window,
            "grid": grid,
            "snapshot_provenance_sha256": expected_provenance_digest,
            "semantics": semantics_measurements,
        }),
    }


def figure_13_analysis_windows_agree(
    first: dict[str, object], second: dict[str, object],
) -> bool:
    """Return whether two accepted windows occupy the same cadence slots."""

    try:
        if int(first["snapshot_count"]) != int(second["snapshot_count"]):
            return False
        for key in (
            "requested_time_start",
            "requested_time_end",
            "selected_time_first",
            "selected_time_last",
        ):
            if not figure_13_time_close(float(first[key]), float(second[key])):
                return False
        first_times = [float(value) for value in first["snapshot_times"]]
        second_times = [float(value) for value in second["snapshot_times"]]
    except (KeyError, TypeError, ValueError):
        return False
    return len(first_times) == len(second_times) and all(
        figure_13_time_close(a, b) for a, b in zip(first_times, second_times)
    )


def figure_13_alternate_firehose_occupancy(
    bundle: Path, manifest: dict[str, object], cases: list[dict[str, object]],
    analyzed_cases: dict[str, object],
    analysis_provenance: dict[str, object],
    execution_authentication: dict[str, object] | None = None,
) -> dict[str, object] | None:
    """Build the named R03/R07/R14/R15 alternate-firehose comparison product."""

    selected_cases: dict[str, dict[str, object]] = {}
    for case in cases:
        case_id = figure_13_case_id(case, manifest, len(cases))
        if case_id is None:
            continue
        if case_id in selected_cases:
            raise ValueError(f"bundle contains duplicate Figure 13 case {case_id}")
        selected_cases[case_id] = case
    if not selected_cases:
        return None

    analysis_configuration = validate_analysis_invocation_provenance(
        analysis_provenance
    )
    current_manifest = bundle_manifest(bundle)
    if current_manifest != manifest:
        raise ValueError("Figure 13 bundle manifest changed while being analyzed")
    if manifest.get("workflow") != STAGE_I_PRODUCTION_WORKFLOW:
        raise ValueError("Figure 13 requires the exact Stage I production workflow")
    manifest_digest = sha256_file(bundle / "manifest.json")
    matrix_provenance = figure_13_matrix_provenance()
    current_execution_authentication = figure_13_execution_authentication(
        bundle, manifest, cases
    )
    if execution_authentication is None:
        execution_authentication = current_execution_authentication
    elif execution_authentication != current_execution_authentication:
        raise ValueError(
            "Figure 13 execution authentication changed while being analyzed"
        )
    selected: dict[str, dict[str, object]] = {}
    for case_id, case in selected_cases.items():
        name = require_manifest_text(case, "name", f"Figure 13 case {case_id}")
        analyzed = analyzed_cases.get(name)
        if not isinstance(analyzed, dict):
            raise ValueError(f"Figure 13 case {case_id} lacks analyzed case output")
        selected[case_id] = figure_13_case_occupancy(
            bundle,
            case_id,
            case,
            analyzed,
            matrix_provenance,
            analysis_configuration,
            execution_authentication,
        )
    ordered = {
        case_id: selected[case_id]
        for case_id in FIGURE_13_FIREHOSE_CASES if case_id in selected
    }
    available = [
        case_id for case_id in FIGURE_13_FIREHOSE_CASES
        if case_id in ordered and ordered[case_id]["available"]
    ]
    missing = [
        case_id for case_id in FIGURE_13_FIREHOSE_CASES if case_id not in ordered
    ]
    unavailable = [
        case_id for case_id in FIGURE_13_FIREHOSE_CASES
        if case_id in ordered and not ordered[case_id]["available"]
    ]
    comparison_window: dict[str, object] | None = None
    comparison_grid: dict[str, object] | None = None
    if available:
        comparison_window = ordered[available[0]]["analysis_window"]
        comparison_grid = ordered[available[0]]["grid"]
        for case_id in available[1:]:
            if not figure_13_analysis_windows_agree(
                ordered[case_id]["analysis_window"], comparison_window
            ):
                raise ValueError(
                    "Figure 13 alternate-firehose cases use inconsistent windows"
                )
            if ordered[case_id]["grid"] != comparison_grid:
                raise ValueError(
                    "Figure 13 alternate-firehose cases use inconsistent grids"
                )
    execution_authenticated = execution_authentication.get("available") is True
    data_complete = not missing and not unavailable
    comparison_ready = data_complete and execution_authenticated
    if comparison_ready:
        revalidate_figure_13_snapshot_bytes(ordered, execution_authentication)
        validate_analysis_invocation_provenance(analysis_provenance)
        if bundle_manifest(bundle) != manifest:
            raise ValueError(
                "Figure 13 bundle manifest changed before comparison-ready publication"
            )
    blockers: list[str] = []
    if missing:
        blockers.append("required Figure 13 cases are missing")
    if unavailable:
        blockers.append("required Figure 13 cases lack exact analysis products")
    if not execution_authenticated:
        blockers.append(
            str(
                execution_authentication.get(
                    "reason", "canonical execution authentication is unavailable"
                )
            )
        )
    published_comparison: dict[str, object]
    if comparison_ready and "R15" in ordered and ordered["R15"]["available"]:
        r15 = ordered["R15"]
        published_comparison = {
            "available": True,
            "case_id": "R15",
            "primary_semantics": "strict_manuscript",
            "primary_temporal_aggregation": "trapezoidal_time_average",
            "parallel_total_unstable_volume_fraction": r15[
                "parallel_total_unstable_volume_fraction"
            ],
            "oblique_total_unstable_volume_fraction": r15[
                "oblique_total_unstable_volume_fraction"
            ],
            "oblique_only_reclassified_volume_fraction": r15[
                "oblique_only_reclassified_volume_fraction"
            ],
            "parallel_total_unstable_percent": r15[
                "parallel_total_unstable_percent"
            ],
            "oblique_total_unstable_percent": r15[
                "oblique_total_unstable_percent"
            ],
            "oblique_only_reclassified_percentage_points": r15[
                "oblique_only_reclassified_percentage_points"
            ],
        }
    else:
        published_comparison = {
            "available": False,
            "case_id": "R15",
            "reason": "comparison-ready R15 occupancy is unavailable",
        }
    product: dict[str, object] = {
        "schema_version": FIGURE_13_FIREHOSE_SCHEMA_VERSION,
        "product_id": "figure_13_alternate_firehose_occupancy",
        "definition": (
            "late-window total microinstability filling fraction, mirror OR "
            "firehose, under parallel and alternative oblique firehose definitions"
        ),
        "interpretation": (
            "both masks are evaluated on the same retained states; the oblique "
            "firehose mask is a postprocessing reclassification, not an "
            "oblique-policy rerun"
        ),
        "normalization": (
            "instantaneous equal-cell-volume fractions on consistent reconstructed "
            "snapshot grids, compared using a trapezoidal time average over the "
            "exact archived t=8..10 window at complete 0.25 cadence; equal-snapshot "
            "means are retained only as secondary diagnostics"
        ),
        "threshold_definitions": instability_threshold_definitions(),
        "primary_semantics": "strict_manuscript",
        "secondary_semantics": "inclusive_solver",
        "primary_temporal_aggregation": "trapezoidal_time_average",
        "secondary_temporal_aggregation": "equal_weight_retained_snapshot_mean",
        "required_case_ids": list(FIGURE_13_FIREHOSE_CASES),
        "available_case_ids": available,
        "missing_case_ids": missing,
        "unavailable_case_ids": unavailable,
        "complete": data_complete,
        "comparison_ready": comparison_ready,
        "comparison_ready_blockers": blockers,
        "comparison_window": comparison_window,
        "comparison_grid": comparison_grid,
        "published_quantitative_comparison": published_comparison,
        "analysis_provenance": analysis_provenance,
        "execution_authentication": execution_authentication,
        "matrix_provenance": matrix_provenance,
        "bundle_provenance": {
            "bundle": str(bundle),
            "bundle_manifest": str(bundle / "manifest.json"),
            "bundle_manifest_sha256_at_analysis": manifest_digest,
            "workflow": manifest.get("workflow"),
            "execution_epoch": manifest.get("execution_epoch"),
            "git_revision": manifest.get("git_revision"),
            "production_case_id": manifest.get("production_case_id"),
        },
        "cases": ordered,
    }
    product["product_sha256"] = canonical_json_sha256(product)
    return product


def write_figure_13_alternate_firehose_occupancy(
    product: dict[str, object], output_dir: Path,
) -> None:
    """Publish a standalone Figure 13 product through an immutable generation."""

    revalidate_figure_13_product_for_publication(product)
    write_figure_13_publication_generation(
        {"figure_13_alternate_firehose_occupancy": product}, output_dir
    )


def figure_13_publication_generation_sha256(result: dict[str, object]) -> str:
    """Hash one diagnostics/named-product pair without self-referential fields."""

    normalized = json.loads(json.dumps(result))
    normalized.pop("figure_13_publication_generation_sha256", None)
    product = normalized.get("figure_13_alternate_firehose_occupancy")
    if not isinstance(product, dict):
        raise ValueError("Figure 13 publication generation lacks its named product")
    product.pop("publication_generation_sha256", None)
    product["product_sha256"] = canonical_json_sha256({
        key: value for key, value in product.items() if key != "product_sha256"
    })
    return canonical_json_sha256({
        "diagnostics": normalized,
        "named_product": product,
    })


def write_figure_13_publication_generation(
    result: dict[str, object], output_dir: Path,
) -> None:
    """Publish diagnostics and its named Figure 13 product through one switch."""

    published_result = json.loads(json.dumps(result))
    product = published_result.get("figure_13_alternate_firehose_occupancy")
    if not isinstance(product, dict):
        raise ValueError("Figure 13 publication generation lacks its named product")
    generation_sha256 = figure_13_publication_generation_sha256(published_result)
    published_result["figure_13_publication_generation_sha256"] = generation_sha256
    product["publication_generation_sha256"] = generation_sha256
    product["product_sha256"] = canonical_json_sha256({
        key: value for key, value in product.items() if key != "product_sha256"
    })
    revalidate_figure_13_product_for_publication(product)
    if figure_13_publication_generation_sha256(published_result) != generation_sha256:
        raise ValueError("Figure 13 publication generation checksum changed")
    write_analysis_publication_generation(published_result, output_dir)


def publication_json_bytes(value: object) -> bytes:
    """Serialize one immutable publication member."""

    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def publication_file_profile(name: str, payload: bytes) -> dict[str, object]:
    """Return the complete immutable profile of one publication member."""

    return {
        "name": name,
        "size_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "mode": "0444",
    }


def ensure_publication_alias(output_dir: Path, name: str) -> None:
    """Install one stable alias that always resolves through analysis-current."""

    destination = output_dir / name
    target = f"{ANALYSIS_PUBLICATION_CURRENT_NAME}/{name}"
    if destination.is_symlink() and os.readlink(destination) == target:
        return
    if destination.exists() and not destination.is_symlink():
        raise ValueError(
            f"analysis publication refuses non-atomic legacy alias migration: "
            f"{destination}"
        )
    temporary = output_dir / f".{name}.alias.{os.getpid()}.{os.urandom(8).hex()}"
    try:
        os.symlink(target, temporary)
        os.replace(temporary, destination)
        fsync_directory(output_dir)
    finally:
        temporary.unlink(missing_ok=True)


def validate_publication_generation(
    generation: Path, expected_manifest_sha256: str | None = None,
) -> dict[str, object]:
    """Authenticate one complete immutable publication generation."""

    manifest_path = generation / ANALYSIS_PUBLICATION_MANIFEST_NAME
    manifest, manifest_sha256 = stable_json_file(
        manifest_path, "analysis publication generation manifest"
    )
    if (
        expected_manifest_sha256 is not None
        and manifest_sha256 != expected_manifest_sha256
    ):
        raise ValueError("analysis publication generation manifest checksum differs")
    if (
        manifest.get("schema_version") != ANALYSIS_PUBLICATION_SCHEMA_VERSION
        or manifest.get("generation") != generation.name
        or manifest.get("mode") not in {"diagnostics-only", "figure-13"}
        or not isinstance(manifest.get("files"), list)
    ):
        raise ValueError("analysis publication generation manifest is malformed")
    names: set[str] = set()
    for profile in manifest["files"]:
        if not isinstance(profile, dict) or set(profile) != {
            "name", "size_bytes", "sha256", "mode"
        }:
            raise ValueError("analysis publication file profile is malformed")
        name = str(profile["name"])
        path = generation / name
        if (
            name in names
            or Path(name).name != name
            or profile["mode"] != "0444"
            or not path.is_file()
            or path.is_symlink()
            or path.stat().st_mode & 0o777 != 0o444
            or path.stat().st_size != profile["size_bytes"]
            or stable_regular_file_sha256(path, "analysis publication member")
            != profile["sha256"]
        ):
            raise ValueError("analysis publication member profile differs")
        names.add(name)
    required = {"diagnostics.json"}
    if manifest["mode"] == "figure-13":
        required.add(FIGURE_13_FIREHOSE_PRODUCT_NAME)
    if names != required:
        raise ValueError("analysis publication generation file set is incomplete")
    if (
        generation.stat().st_mode & 0o777 != 0o555
        or manifest_path.stat().st_mode & 0o777 != 0o444
    ):
        raise ValueError("analysis publication generation modes differ")
    diagnostics, _ = stable_json_file(
        generation / "diagnostics.json", "analysis publication diagnostics"
    )
    generation_digest = generation.name.removeprefix("generation-")
    if manifest["mode"] == "figure-13":
        product, _ = stable_json_file(
            generation / FIGURE_13_FIREHOSE_PRODUCT_NAME,
            "analysis publication named product",
        )
        if (
            diagnostics.get("figure_13_alternate_firehose_occupancy") != product
            or diagnostics.get("figure_13_publication_generation_sha256")
            != generation_digest
            or figure_13_publication_generation_sha256(diagnostics)
            != generation_digest
        ):
            raise ValueError("Figure 13 publication generation content differs")
    elif canonical_json_sha256(diagnostics) != generation_digest:
        raise ValueError("diagnostics-only publication generation content differs")
    return manifest


def read_analysis_publication(output_dir: Path) -> dict[str, object]:
    """Resolve and authenticate the generation selected by analysis-current."""

    current = output_dir / ANALYSIS_PUBLICATION_CURRENT_NAME
    if not current.is_symlink():
        raise ValueError("analysis publication lacks its atomic current pointer")
    target = Path(os.readlink(current))
    if target.is_absolute() or ".." in target.parts:
        raise ValueError("analysis publication current pointer escapes its output")
    generation = (output_dir / target).resolve()
    generations = (output_dir / ANALYSIS_PUBLICATION_GENERATIONS_NAME).resolve()
    if not path_is_beneath(generation, generations) or not generation.is_dir():
        raise ValueError("analysis publication current pointer is invalid")
    manifest = validate_publication_generation(generation)
    return {
        "generation_dir": generation,
        "manifest": manifest,
        "diagnostics": generation / "diagnostics.json",
        "named_product": (
            generation / FIGURE_13_FIREHOSE_PRODUCT_NAME
            if manifest["mode"] == "figure-13" else None
        ),
    }


def write_analysis_publication_generation(
    result: dict[str, object], output_dir: Path,
) -> Path:
    """Publish one immutable generation with one atomic current-pointer switch."""

    published = json.loads(json.dumps(result))
    product = published.get("figure_13_alternate_firehose_occupancy")
    mode = "figure-13" if isinstance(product, dict) else "diagnostics-only"
    generation_digest = (
        str(published.get("figure_13_publication_generation_sha256"))
        if mode == "figure-13"
        else canonical_json_sha256(published)
    )
    if not re.fullmatch(r"[0-9a-f]{64}", generation_digest):
        raise ValueError("analysis publication generation digest is malformed")
    members: dict[str, bytes] = {
        "diagnostics.json": publication_json_bytes(published),
    }
    if mode == "figure-13":
        members[FIGURE_13_FIREHOSE_PRODUCT_NAME] = publication_json_bytes(product)
    profiles = [
        publication_file_profile(name, payload)
        for name, payload in sorted(members.items())
    ]
    generation_name = f"generation-{generation_digest}"
    manifest = {
        "schema_version": ANALYSIS_PUBLICATION_SCHEMA_VERSION,
        "generation": generation_name,
        "mode": mode,
        "files": profiles,
    }
    manifest_payload = publication_json_bytes(manifest)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in ("diagnostics.json", FIGURE_13_FIREHOSE_PRODUCT_NAME):
        destination = output_dir / name
        if destination.exists() and not destination.is_symlink():
            raise ValueError(
                f"analysis publication refuses non-atomic legacy alias migration: "
                f"{destination}"
            )
    generations = output_dir / ANALYSIS_PUBLICATION_GENERATIONS_NAME
    generations.mkdir(exist_ok=True)
    generation = generations / generation_name
    if generation.exists():
        if generation.is_symlink():
            raise ValueError("analysis publication generation must not be a symlink")
        validate_publication_generation(
            generation, hashlib.sha256(manifest_payload).hexdigest()
        )
    else:
        staged = Path(tempfile.mkdtemp(prefix=".stage-", dir=generations))
        try:
            for name, payload in {
                **members,
                ANALYSIS_PUBLICATION_MANIFEST_NAME: manifest_payload,
            }.items():
                path = staged / name
                with path.open("xb") as stream:
                    stream.write(payload)
                    stream.flush()
                    os.fsync(stream.fileno())
                path.chmod(0o444)
            fsync_directory(staged)
            staged.chmod(0o555)
            os.replace(staged, generation)
            fsync_directory(generations)
        finally:
            if staged.exists():
                staged.chmod(0o755)
                shutil.rmtree(staged)
    validate_publication_generation(
        generation, hashlib.sha256(manifest_payload).hexdigest()
    )

    for name in ("diagnostics.json", FIGURE_13_FIREHOSE_PRODUCT_NAME):
        ensure_publication_alias(output_dir, name)
    current = output_dir / ANALYSIS_PUBLICATION_CURRENT_NAME
    temporary = output_dir / (
        f".{ANALYSIS_PUBLICATION_CURRENT_NAME}.{os.getpid()}.{os.urandom(8).hex()}"
    )
    target = generation.relative_to(output_dir)
    try:
        os.symlink(str(target), temporary)
        os.replace(temporary, current)
        fsync_directory(output_dir)
    finally:
        temporary.unlink(missing_ok=True)
    for name in ("diagnostics.json", FIGURE_13_FIREHOSE_PRODUCT_NAME):
        ensure_publication_alias(output_dir, name)
    selected = read_analysis_publication(output_dir)
    if selected["generation_dir"] != generation:
        raise ValueError("analysis publication current pointer selected another generation")
    return generation


def revalidate_controller_artifact_files(record: dict[str, object]) -> None:
    """Revalidate live executable, approval, and restart bytes before publication."""

    context = f"Figure 13 controller artifact {record.get('case_id', '')}"
    executable = authenticate_regular_file(
        Path(str(record.get("executable_path", ""))),
        record.get("executable_sha256"),
        record.get("executable_size_bytes"),
        f"{context} executable",
    )
    if not os.access(Path(executable["path"]), os.X_OK):
        raise ValueError(f"{context} executable is not executable before publication")
    approval, approval_digest = stable_json_file(
        Path(str(record.get("qualification_approval_path", ""))),
        f"{context} qualification approval",
    )
    if (
        approval_digest != record.get("qualification_approval_sha256")
        or approval.get("approved_executable") != executable["path"]
        or approval.get("approved_executable_revision")
        != record.get("executable_revision")
        or approval.get("approved_executable_sha256") != executable["sha256"]
    ):
        raise ValueError(
            f"{context} executable approval changed before publication"
        )
    try:
        expected_ranks = int(record["expected_rank_count"])
        final_time = float(record["final_time"])
        case_root = Path(str(record["controller_case_root"]))
        canonical_root = Path(str(record["canonical_root"]))
        executable_revision = str(record["executable_revision"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{context} retained restart evidence is malformed") from error
    restart_binary_abi = qualified_restart_binary_abi(
        executable_revision, str(executable["sha256"]), context
    )
    restart_load_evidence = authenticate_qualified_restart_load_evidence(
        approval, canonical_root, executable, executable_revision, context
    )
    if (
        restart_binary_abi != record.get("restart_binary_abi")
        or restart_load_evidence != record.get("qualified_restart_load_evidence")
    ):
        raise ValueError(f"{context} restart qualification changed before publication")
    input_contract, runtime_contract, _ = retained_controller_restart_contract(
        record, context
    )
    terminal_records, terminal_loadability = authenticate_restart_rank_files(
        record.get("terminal_restart_rank_files"),
        expected_ranks,
        case_root,
        final_time,
        restart_binary_abi,
        input_contract,
        runtime_contract,
        f"{context} terminal restart",
    )
    if (
        terminal_records != record.get("terminal_restart_rank_files")
        or terminal_loadability != record.get("terminal_restart_loadability")
    ):
        raise ValueError(f"{context} terminal restart changed before publication")
    archive = record.get("restart_archive_rank_files")
    if archive is not None:
        archive_contract_record = {
            "path": record.get("restart_archive_parent_manifest"),
            "case_id": record.get("case_id"),
            "case_name": record.get("case_name"),
            "input_sha256": record.get("input_sha256"),
            "restart_parameter_contract": record.get(
                "restart_archive_parameter_contract"
            ),
        }
        archive_input_contract, archive_runtime_contract, _ = (
            retained_controller_restart_contract(
                archive_contract_record, f"{context} archived restart parent"
            )
        )
        archive_records, archive_loadability = authenticate_restart_rank_files(
            archive,
            expected_ranks,
            case_root,
            float(record.get("restart_archive_time", final_time)),
            restart_binary_abi,
            archive_input_contract,
            archive_runtime_contract,
            f"{context} archived restart",
        )
        if (
            archive_records != archive
            or archive_loadability != record.get("restart_archive_loadability")
        ):
            raise ValueError(f"{context} archived restart changed before publication")


def revalidate_figure_13_product_for_publication(
    product: dict[str, object],
) -> None:
    """Revalidate comparison-ready trust bindings immediately before publication."""

    if product.get("comparison_ready") is not True:
        return
    expected_product_digest = canonical_json_sha256({
        key: value for key, value in product.items() if key != "product_sha256"
    })
    if product.get("product_sha256") != expected_product_digest:
        raise ValueError("Figure 13 product checksum changed before publication")
    provenance = product.get("analysis_provenance")
    authentication = product.get("execution_authentication")
    cases = product.get("cases")
    bundle_provenance = product.get("bundle_provenance")
    if not all(
        isinstance(value, dict)
        for value in (provenance, authentication, cases, bundle_provenance)
    ):
        raise ValueError("Figure 13 comparison-ready product trust evidence is malformed")
    validate_analysis_invocation_provenance(provenance)
    revalidate_figure_13_snapshot_bytes(cases, authentication)
    manifest_path = Path(str(bundle_provenance.get("bundle_manifest", "")))
    _, manifest_digest = stable_json_file(
        manifest_path, "Figure 13 bundle manifest"
    )
    if manifest_digest != bundle_provenance.get("bundle_manifest_sha256_at_analysis"):
        raise ValueError("Figure 13 bundle manifest changed before publication")
    controller_artifacts = authentication.get("controller_artifacts")
    if not isinstance(controller_artifacts, dict):
        raise ValueError("Figure 13 execution authentication lacks controller artifacts")
    source_bundles: dict[str, str] = {}
    for records in controller_artifacts.values():
        if not isinstance(records, list):
            raise ValueError("Figure 13 controller artifact evidence is malformed")
        for record in records:
            if not isinstance(record, dict):
                raise ValueError("Figure 13 controller artifact record is malformed")
            _, digest = stable_json_file(
                Path(str(record.get("path", ""))),
                "Figure 13 controller artifact",
            )
            if digest != record.get("sha256"):
                raise ValueError("Figure 13 controller artifact changed before publication")
            revalidate_controller_artifact_files(record)
            source_path = str(record.get("source_bundle_path", ""))
            source_digest = str(record.get("source_bundle_sha256", ""))
            source_bundles[source_path] = source_digest
    for source_path, expected_digest in source_bundles.items():
        if stable_regular_file_sha256(
            Path(source_path), "Figure 13 source bundle"
        ) != expected_digest:
            raise ValueError("Figure 13 source bundle changed before publication")


def validate_digitized_source_figure(manifest_path: Path,
                                     record: dict[str, object],
                                     context: str) -> str:
    """Validate one source figure named by digitized reference provenance."""

    figure_name = require_manifest_text(record, "source_figure", context)
    figure_digest = require_manifest_text(record, "source_figure_sha256", context)
    if not re.fullmatch(r"[0-9a-f]{64}", figure_digest):
        raise ValueError("digitized source_figure_sha256 must be lowercase SHA-256")
    figure_path = (manifest_path.parent / figure_name).resolve()
    if not figure_path.is_file():
        raise ValueError(f"digitized source figure is missing: {figure_path}")
    if figure_digest != sha256_file(figure_path):
        raise ValueError("digitized source figure checksum does not match")
    return figure_name


def alignment_peak_curve(ensemble: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    """Return physical-k peak alignment cosine values from selected-shell PDFs."""

    alignment = ensemble.get("alignment", {})
    spectra = ensemble.get("spectra", {})
    if not isinstance(alignment, dict) or not isinstance(spectra, dict):
        raise ValueError("analyzed product is missing: alignment_peak.cos_theta")
    velocity = spectra.get("velocity", {})
    if not isinstance(velocity, dict) or "dk" not in velocity:
        raise ValueError("analyzed product is missing: alignment_peak.cos_theta")
    dk = float(velocity["dk"])
    peaks: list[tuple[float, float]] = []
    for shell, record in sorted(alignment.items(), key=lambda item: int(item[0])):
        edges = np.asarray(record["edges"], dtype=float)
        density = np.asarray(record["density"], dtype=float)
        centers = 0.5 * (edges[1:] + edges[:-1])
        peaks.append((float(shell) * dk, float(centers[np.argmax(density)])))
    if len(peaks) < 2:
        raise ValueError(
            "analyzed product requires at least two shells: "
            "alignment_peak.cos_theta"
        )
    return (
        np.asarray([peak[0] for peak in peaks], dtype=float),
        np.asarray([peak[1] for peak in peaks], dtype=float),
    )


def analyzed_product_curve(ensemble: dict[str, object], product: str,
                           histories: object = None
                           ) -> tuple[np.ndarray, np.ndarray]:
    """Return analyzed x/y arrays selected by a reference product identifier."""

    parts = product.split(".", maxsplit=1)
    if len(parts) != 2:
        raise ValueError(f"unsupported reference product: {product}")
    family, name = parts
    if family == "spectra":
        products = ensemble.get("spectra", {})
        if not isinstance(products, dict) or name not in products:
            raise ValueError(f"analyzed product is missing: {product}")
        record = products[name]
        return (
            np.asarray(record["k"], dtype=float),
            np.asarray(record["power_per_dk"], dtype=float),
        )
    if family == "pdf":
        products = ensemble.get("pdf", {})
        if not isinstance(products, dict) or name not in products:
            raise ValueError(f"analyzed product is missing: {product}")
        record = products[name]
        edges = np.asarray(record["edges"], dtype=float)
        return 0.5 * (edges[1:] + edges[:-1]), np.asarray(record["density"], dtype=float)
    if product in ("pressure_transfer.transfer", "pressure_transfer.signed_transfer"):
        record = ensemble.get("pressure_transfer", {})
        if not isinstance(record, dict):
            raise ValueError(f"analyzed product is missing: {product}")
        key = "signed_transfer" if product.endswith("signed_transfer") else "transfer"
        if key not in record:
            key = "transfer"
        return (
            np.asarray(record["k_perp"], dtype=float),
            np.asarray(record[key], dtype=float),
        )
    if product in (
        "pressure_transfer.transfer_normalized_by_total",
        "pressure_transfer.signed_transfer_normalized_by_total",
    ):
        record = ensemble.get("pressure_transfer", {})
        if (
            not isinstance(record, dict)
            or not record.get("normalization_available", False)
            or record.get("transfer_normalized_by_total") is None
        ):
            raise ValueError(f"analyzed product is missing: {product}")
        key = (
            "signed_transfer_normalized_by_total"
            if product.endswith("signed_transfer_normalized_by_total")
            else "transfer_normalized_by_total"
        )
        if record.get(key) is None:
            key = "transfer_normalized_by_total"
        return (
            np.asarray(record["k_perp"], dtype=float),
            np.asarray(record[key], dtype=float),
        )
    if family == "alignment":
        products = ensemble.get("alignment", {})
        if not isinstance(products, dict) or name not in products:
            raise ValueError(f"analyzed product is missing: {product}")
        record = products[name]
        edges = np.asarray(record["edges"], dtype=float)
        return 0.5 * (edges[1:] + edges[:-1]), np.asarray(record["density"], dtype=float)
    if family == "alignment_peak" and name == "cos_theta":
        return alignment_peak_curve(ensemble)
    if family == "eddy_anisotropy":
        products = ensemble.get("eddy_anisotropy", {})
        if not isinstance(products, dict) or name not in products:
            raise ValueError(f"analyzed product is missing: {product}")
        record = products[name]
        if not isinstance(record, dict) or not record.get("available", False):
            raise ValueError(f"analyzed product is missing: {product}")
        return (
            np.asarray(record["ell_perp_over_lperp"], dtype=float),
            np.asarray(record["ell_parallel_over_lperp"], dtype=float),
        )
    if family == "history":
        if not isinstance(histories, list) or len(histories) != 1:
            raise ValueError(
                f"analyzed product requires exactly one history: {product}"
            )
        time_series = histories[0].get("time_series", {})
        if not isinstance(time_series, dict) or name not in time_series:
            raise ValueError(f"analyzed product is missing: {product}")
        return (
            np.asarray(time_series["time"], dtype=float),
            np.asarray(time_series[name], dtype=float),
        )
    raise ValueError(f"unsupported reference product: {product}")


def interpolate_analysis_curve(x: np.ndarray, source_x: np.ndarray,
                               source_y: np.ndarray, method: str) -> np.ndarray:
    """Interpolate analyzed values at reference coordinates."""

    selected = np.isfinite(source_x) & np.isfinite(source_y)
    if method == "loglog":
        selected &= (source_x > 0.0) & (source_y > 0.0)
        if np.any(x <= 0.0):
            raise ValueError("loglog reference coordinates must be positive")
    source_x = source_x[selected]
    source_y = source_y[selected]
    if source_x.size < 2 or np.any(np.diff(source_x) <= 0.0):
        raise ValueError("analyzed reference product needs ordered finite samples")
    if np.any(x < source_x[0]) or np.any(x > source_x[-1]):
        raise ValueError("reference coordinates lie outside analyzed product range")
    if method == "linear":
        return np.interp(x, source_x, source_y)
    if method == "loglog":
        return np.exp(np.interp(np.log(x), np.log(source_x), np.log(source_y)))
    raise ValueError(f"unsupported reference interpolation: {method}")


def analyzed_product_surface(ensemble: dict[str, object], product: str
                             ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return analyzed x/y centers and density for a reference surface."""

    parts = product.split(".", maxsplit=1)
    if len(parts) != 2 or parts[0] != "pressure_density_joint":
        raise ValueError(f"unsupported reference surface product: {product}")
    products = ensemble.get("pressure_density_joint", {})
    if not isinstance(products, dict) or parts[1] not in products:
        raise ValueError(f"analyzed product is missing: {product}")
    record = products[parts[1]]
    x_edges = np.asarray(record["x_edges"], dtype=float)
    y_edges = np.asarray(record["y_edges"], dtype=float)
    density = np.asarray(record["density"], dtype=float)
    return (
        0.5 * (x_edges[1:] + x_edges[:-1]),
        0.5 * (y_edges[1:] + y_edges[:-1]),
        density,
    )


def interpolate_analysis_surface(
    x: np.ndarray, y: np.ndarray, source_x: np.ndarray, source_y: np.ndarray,
    source_z: np.ndarray, method: str,
) -> np.ndarray:
    """Interpolate one rectangular analyzed surface at reference samples."""

    if method != "bilinear":
        raise ValueError(f"unsupported reference surface interpolation: {method}")
    if (
        x.shape != y.shape or source_x.size < 2 or source_y.size < 2
        or source_z.shape != (source_x.size, source_y.size)
        or not np.isfinite(source_x).all() or not np.isfinite(source_y).all()
        or not np.isfinite(source_z).all()
        or np.any(np.diff(source_x) <= 0.0) or np.any(np.diff(source_y) <= 0.0)
    ):
        raise ValueError("analyzed reference surface needs a finite ordered grid")
    if (
        np.any(x < source_x[0]) or np.any(x > source_x[-1])
        or np.any(y < source_y[0]) or np.any(y > source_y[-1])
    ):
        raise ValueError("reference coordinates lie outside analyzed surface range")
    values = []
    for x_value, y_value in zip(x, y):
        row = np.asarray([
            np.interp(x_value, source_x, source_z[:, index])
            for index in range(source_y.size)
        ])
        values.append(np.interp(y_value, source_y, row))
    return np.asarray(values, dtype=float)


def reference_curve_comparisons(
    result: dict[str, object], manifest_path: Path,
    allow_missing_cases: bool = False,
    analysis_case_aliases: dict[str, str] | None = None,
    stage_i_reference_bindings: dict[str, dict[str, str]] | None = None,
) -> dict[str, object]:
    """Compare analysis products against provenance-qualified external data."""

    if analysis_case_aliases is None:
        analysis_case_aliases = {}
    manifest_sha256 = sha256_file(manifest_path)
    with manifest_path.open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    if manifest.get("schema_version") != REFERENCE_CURVE_SCHEMA_VERSION:
        raise ValueError(
            f"reference curves require schema_version={REFERENCE_CURVE_SCHEMA_VERSION}"
        )
    provenance = manifest.get("provenance", {})
    if not isinstance(provenance, dict):
        raise ValueError("reference curves require provenance metadata")
    method = require_manifest_text(provenance, "method", "reference provenance")
    if method not in ("machine_readable", "digitized"):
        raise ValueError(f"unsupported reference provenance method: {method}")
    require_manifest_text(provenance, "source_description", "reference provenance")
    require_manifest_text(
        provenance, "uncertainty_description", "reference provenance"
    )
    if method == "digitized":
        source_figures = provenance.get("source_figures")
        if source_figures is None:
            validate_digitized_source_figure(
                manifest_path, provenance, "digitized provenance"
            )
        else:
            if not isinstance(source_figures, list) or not source_figures:
                raise ValueError(
                    "digitized provenance source_figures must be a nonempty list"
                )
            figure_names: set[str] = set()
            for index, figure in enumerate(source_figures):
                if not isinstance(figure, dict):
                    raise ValueError(
                        "digitized provenance source_figures entries must be objects"
                    )
                name = validate_digitized_source_figure(
                    manifest_path, figure,
                    f"digitized provenance source_figures[{index}]",
                )
                if name in figure_names:
                    raise ValueError("digitized source figure is duplicated")
                figure_names.add(name)
        require_manifest_text(provenance, "digitization_tool", "digitized provenance")
    curves = manifest.get("curves", [])
    surfaces = manifest.get("surfaces", [])
    if not isinstance(curves, list) or not isinstance(surfaces, list):
        raise ValueError("reference manifest curves and surfaces must be lists")
    if not curves and not surfaces:
        raise ValueError("reference manifest must contain curves or surfaces")
    cases = result.get("cases", {})
    comparisons: dict[str, object] = {}
    omitted_products: list[dict[str, str]] = []
    product_ids: set[str] = set()
    for curve in curves:
        if not isinstance(curve, dict):
            raise ValueError("reference curve entries must be objects")
        curve_id = require_manifest_text(curve, "id", "reference curve")
        if curve_id in product_ids:
            raise ValueError(f"reference curve id is duplicated: {curve_id}")
        product_ids.add(curve_id)
        case = require_manifest_text(curve, "case", f"reference curve {curve_id}")
        product = require_manifest_text(curve, "product", f"reference curve {curve_id}")
        data_file = require_manifest_text(
            curve, "data_file", f"reference curve {curve_id}"
        )
        expected_digest = require_manifest_text(
            curve, "data_sha256", f"reference curve {curve_id}"
        )
        binding_validated = validate_stage_i_reference_binding(
            stage_i_reference_bindings, curve_id, "curve", case, product,
            data_file, expected_digest, manifest_sha256,
        )
        data_path = (manifest_path.parent / data_file).resolve()
        if not data_path.is_file():
            raise ValueError(f"reference curve data is missing: {data_path}")
        if expected_digest != sha256_file(data_path):
            raise ValueError(f"reference curve checksum does not match: {curve_id}")
        table = np.genfromtxt(data_path, delimiter=",", names=True, encoding="utf-8")
        rows = np.atleast_1d(table)
        if table.dtype.names is None or not {"x", "y", "y_uncertainty"}.issubset(
            table.dtype.names
        ):
            raise ValueError(
                f"reference curve {curve_id} requires x,y,y_uncertainty columns"
            )
        x = np.asarray(rows["x"], dtype=float)
        y = np.asarray(rows["y"], dtype=float)
        uncertainty = np.asarray(rows["y_uncertainty"], dtype=float)
        if (
            x.size < 2 or not np.isfinite(x).all() or not np.isfinite(y).all()
            or not np.isfinite(uncertainty).all() or np.any(uncertainty <= 0.0)
            or np.any(np.diff(x) <= 0.0)
        ):
            raise ValueError(
                f"reference curve {curve_id} needs ordered finite points and "
                "positive y_uncertainty"
            )
        analysis_case = analysis_case_aliases.get(case, case)
        if case == "direct":
            ensemble = result.get("snapshot_ensemble", {})
            histories = result.get("histories", [])
        elif isinstance(cases, dict) and analysis_case in cases:
            ensemble = cases[analysis_case].get("snapshot_ensemble", {})
            histories = cases[analysis_case].get("histories", [])
        else:
            if allow_missing_cases:
                omitted_products.append({
                    "id": curve_id,
                    "kind": "curve",
                    "case": case,
                    "analysis_case": analysis_case,
                    "product": product,
                    "stage_i_binding_validated": binding_validated,
                    "reason": "referenced case is absent from analyzed bundle",
                })
                continue
            raise ValueError(f"reference curve {curve_id} selects missing case {case}")
        if (
            not product.startswith("history.")
            and (
                not isinstance(ensemble, dict)
                or ensemble.get("snapshot_count", 0) == 0
            )
        ):
            raise ValueError(f"reference curve {curve_id} selects empty analysis")
        if not isinstance(ensemble, dict):
            ensemble = {}
        source_x, source_y = analyzed_product_curve(ensemble, product, histories)
        simulated = interpolate_analysis_curve(
            x, source_x, source_y, str(curve.get("interpolation", "linear"))
        )
        residual = simulated - y
        normalized = residual / uncertainty
        comparisons[curve_id] = {
            "available": True,
            "kind": "curve",
            "case": case,
            "analysis_case": analysis_case,
            "product": product,
            "stage_i_binding_validated": binding_validated,
            "reference_data_file": data_file,
            "reference_manifest_sha256": manifest_sha256,
            "data_file": str(data_path),
            "data_sha256": expected_digest,
            "interpolation": str(curve.get("interpolation", "linear")),
            "sample_count": int(x.size),
            "x": x.tolist(),
            "reference_y": y.tolist(),
            "reference_y_uncertainty": uncertainty.tolist(),
            "simulated_y": simulated.tolist(),
            "residual": residual.tolist(),
            "rms_residual": float(np.sqrt(np.mean(residual ** 2))),
            "maximum_absolute_residual": float(np.max(np.abs(residual))),
            "rms_normalized_by_reported_uncertainty": float(
                np.sqrt(np.mean(normalized ** 2))
            ),
        }
    surface_comparisons: dict[str, object] = {}
    for surface in surfaces:
        if not isinstance(surface, dict):
            raise ValueError("reference surface entries must be objects")
        surface_id = require_manifest_text(surface, "id", "reference surface")
        if surface_id in product_ids:
            raise ValueError(f"reference product id is duplicated: {surface_id}")
        product_ids.add(surface_id)
        case = require_manifest_text(surface, "case", f"reference surface {surface_id}")
        product = require_manifest_text(
            surface, "product", f"reference surface {surface_id}"
        )
        data_file = require_manifest_text(
            surface, "data_file", f"reference surface {surface_id}"
        )
        expected_digest = require_manifest_text(
            surface, "data_sha256", f"reference surface {surface_id}"
        )
        binding_validated = validate_stage_i_reference_binding(
            stage_i_reference_bindings, surface_id, "surface", case, product,
            data_file, expected_digest, manifest_sha256,
        )
        data_path = (manifest_path.parent / data_file).resolve()
        if not data_path.is_file():
            raise ValueError(f"reference surface data is missing: {data_path}")
        if expected_digest != sha256_file(data_path):
            raise ValueError(
                f"reference surface checksum does not match: {surface_id}"
            )
        table = np.genfromtxt(data_path, delimiter=",", names=True, encoding="utf-8")
        rows = np.atleast_1d(table)
        required = {"x", "y", "z", "z_uncertainty"}
        if table.dtype.names is None or not required.issubset(table.dtype.names):
            raise ValueError(
                f"reference surface {surface_id} requires "
                "x,y,z,z_uncertainty columns"
            )
        x = np.asarray(rows["x"], dtype=float)
        y = np.asarray(rows["y"], dtype=float)
        z = np.asarray(rows["z"], dtype=float)
        uncertainty = np.asarray(rows["z_uncertainty"], dtype=float)
        if (
            x.size < 1 or not np.isfinite(x).all() or not np.isfinite(y).all()
            or not np.isfinite(z).all() or not np.isfinite(uncertainty).all()
            or np.any(uncertainty <= 0.0)
        ):
            raise ValueError(
                f"reference surface {surface_id} needs finite samples and "
                "positive z_uncertainty"
            )
        analysis_case = analysis_case_aliases.get(case, case)
        if case == "direct":
            ensemble = result.get("snapshot_ensemble", {})
        elif isinstance(cases, dict) and analysis_case in cases:
            ensemble = cases[analysis_case].get("snapshot_ensemble", {})
        else:
            if allow_missing_cases:
                omitted_products.append({
                    "id": surface_id,
                    "kind": "surface",
                    "case": case,
                    "analysis_case": analysis_case,
                    "product": product,
                    "stage_i_binding_validated": binding_validated,
                    "reason": "referenced case is absent from analyzed bundle",
                })
                continue
            raise ValueError(
                f"reference surface {surface_id} selects missing case {case}"
            )
        if (
            not isinstance(ensemble, dict)
            or ensemble.get("snapshot_count", 0) == 0
        ):
            raise ValueError(f"reference surface {surface_id} selects empty analysis")
        source_x, source_y, source_z = analyzed_product_surface(ensemble, product)
        simulated = interpolate_analysis_surface(
            x, y, source_x, source_y, source_z,
            str(surface.get("interpolation", "bilinear")),
        )
        residual = simulated - z
        normalized = residual / uncertainty
        surface_comparisons[surface_id] = {
            "available": True,
            "kind": "surface",
            "case": case,
            "analysis_case": analysis_case,
            "product": product,
            "stage_i_binding_validated": binding_validated,
            "reference_data_file": data_file,
            "reference_manifest_sha256": manifest_sha256,
            "data_file": str(data_path),
            "data_sha256": expected_digest,
            "interpolation": str(surface.get("interpolation", "bilinear")),
            "sample_count": int(x.size),
            "x": x.tolist(),
            "y": y.tolist(),
            "reference_z": z.tolist(),
            "reference_z_uncertainty": uncertainty.tolist(),
            "simulated_z": simulated.tolist(),
            "residual": residual.tolist(),
            "rms_residual": float(np.sqrt(np.mean(residual ** 2))),
            "maximum_absolute_residual": float(np.max(np.abs(residual))),
            "rms_normalized_by_reported_uncertainty": float(
                np.sqrt(np.mean(normalized ** 2))
            ),
        }
    return {
        "available": True,
        "definition": (
            "analysis products interpolated onto provenance-qualified reference "
            "curve or surface coordinates"
        ),
        "manifest": str(manifest_path),
        "manifest_sha256": manifest_sha256,
        "provenance": provenance,
        "comparisons": comparisons,
        "surface_comparisons": surface_comparisons,
        "omitted_products": omitted_products,
        "allow_missing_cases": allow_missing_cases,
    }


def combined_reference_curve_comparisons(
    result: dict[str, object], manifest_paths: list[Path],
    allow_missing_cases: bool = False,
    analysis_case_aliases: dict[str, str] | None = None,
    stage_i_reference_bindings: dict[str, dict[str, str]] | None = None,
) -> dict[str, object]:
    """Combine distinct comparison products from qualified reference manifests."""

    manifests: list[dict[str, object]] = []
    comparisons: dict[str, object] = {}
    surface_comparisons: dict[str, object] = {}
    omitted_products: list[dict[str, str]] = []
    combined: dict[str, object] = {
        "available": True,
        "definition": (
            "analysis products interpolated onto provenance-qualified reference "
            "curve or surface coordinates"
        ),
        "manifests": manifests,
        "comparisons": comparisons,
        "surface_comparisons": surface_comparisons,
        "omitted_products": omitted_products,
        "allow_missing_cases": allow_missing_cases,
    }
    product_ids: set[str] = set()
    for manifest_path in manifest_paths:
        comparison = reference_curve_comparisons(
            result, manifest_path, allow_missing_cases=allow_missing_cases,
            analysis_case_aliases=analysis_case_aliases,
            stage_i_reference_bindings=stage_i_reference_bindings,
        )
        manifests.append({
            "manifest": comparison["manifest"],
            "manifest_sha256": comparison["manifest_sha256"],
            "provenance": comparison["provenance"],
        })
        all_ids = [
            *comparison["comparisons"].keys(),
            *comparison["surface_comparisons"].keys(),
            *[item["id"] for item in comparison["omitted_products"]],
        ]
        for product_id in all_ids:
            if product_id in product_ids:
                raise ValueError(
                    f"reference product id is duplicated across manifests: "
                    f"{product_id}"
                )
            product_ids.add(product_id)
        omitted_products.extend({
            **item,
            "manifest": comparison["manifest"],
        } for item in comparison["omitted_products"])
        for collection, retained in (
            ("comparisons", comparisons),
            ("surface_comparisons", surface_comparisons),
        ):
            products = comparison.get(collection, {})
            if not isinstance(products, dict):
                raise ValueError(f"reference comparison {collection} is not an object")
            for product_id, product in products.items():
                if not isinstance(product_id, str):
                    raise ValueError("reference product id is not text")
                retained[product_id] = product
    if len(manifest_paths) == 1:
        only = manifests[0]
        combined.update({
            "manifest": only["manifest"],
            "manifest_sha256": only["manifest_sha256"],
            "provenance": only["provenance"],
        })
    return combined


def require_text_list(record: dict[str, object], key: str, context: str,
                      required: bool = False) -> list[str]:
    """Return one unique list of nonempty manifest text values."""

    values = record.get(key, [])
    if not isinstance(values, list) or (
        required and not values
    ) or not all(isinstance(value, str) and value.strip() for value in values):
        qualifier = "nonempty " if required else ""
        raise ValueError(f"{context} requires a {qualifier}list of text {key}")
    retained = [str(value) for value in values]
    if len(set(retained)) != len(retained):
        raise ValueError(f"{context} {key} entries must be unique")
    return retained


def require_sha256(record: dict[str, object], key: str, context: str) -> str:
    """Return one required lowercase SHA-256 manifest value."""

    value = require_manifest_text(record, key, context)
    if not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError(f"{context} {key} must be lowercase SHA-256")
    return value


def validate_stage_i_reference_binding(
    bindings: dict[str, dict[str, str]] | None, product_id: str, kind: str,
    case: str, product: str, data_file: str, data_sha256: str,
    reference_manifest_sha256: str,
) -> bool:
    """Reject substituted semantics or provenance for one approved product id."""

    if bindings is None or product_id not in bindings:
        return False
    expected = bindings[product_id]
    observed = {
        "kind": kind,
        "case": case,
        "product": product,
        "data_file": data_file,
        "data_sha256": data_sha256,
        "reference_manifest_sha256": reference_manifest_sha256,
    }
    mismatches = sorted(
        key for key, value in observed.items() if expected.get(key) != value
    )
    if mismatches:
        raise ValueError(
            f"Stage I reference binding mismatch for {product_id}: "
            f"{', '.join(mismatches)}"
        )
    return True


def stage_i_panels_configuration(manifest_path: Path) -> dict[str, object]:
    """Read the versioned Stage I panel-gate configuration."""

    with manifest_path.open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    if not isinstance(manifest, dict):
        raise ValueError("Stage I manifest must be an object")
    configuration = manifest.get("panel_status")
    if configuration is None:
        configuration = {
            "schema_version": manifest.get("panels_schema_version"),
            "panels": manifest.get("panels"),
        }
    if not isinstance(configuration, dict):
        raise ValueError("Stage I manifest panel_status must be an object")
    if configuration.get("schema_version") != STAGE_I_PANEL_SCHEMA_VERSION:
        raise ValueError(
            "Stage I panels require panels_schema_version="
            f"{STAGE_I_PANEL_SCHEMA_VERSION}"
        )
    cases = manifest.get("cases", [])
    panels = configuration.get("panels", [])
    analysis_case_aliases = configuration.get("analysis_case_aliases", {})
    reference_manifests = configuration.get("reference_manifests", {})
    reference_product_bindings = configuration.get("reference_product_bindings", {})
    if not isinstance(cases, list) or not all(
        isinstance(case, dict) for case in cases
    ):
        raise ValueError("Stage I manifest cases must be a list of objects")
    if not isinstance(panels, list) or not panels or not all(
        isinstance(panel, dict) for panel in panels
    ):
        raise ValueError("Stage I manifest panels must be a nonempty list of objects")
    case_names: dict[str, str] = {}
    case_ids: set[str] = set()
    for case in cases:
        case_id = require_manifest_text(case, "id", "Stage I case")
        name = require_manifest_text(case, "name", f"Stage I case {case_id}")
        if case_id in case_ids or name in case_names:
            raise ValueError("Stage I case ids and names must be unique")
        case_ids.add(case_id)
        case_names[name] = case_id
    if not isinstance(analysis_case_aliases, dict) or not all(
        isinstance(alias, str) and alias.strip()
        and isinstance(canonical, str) and canonical.strip()
        for alias, canonical in analysis_case_aliases.items()
    ):
        raise ValueError("Stage I analysis_case_aliases must map text names to names")
    for alias, canonical in analysis_case_aliases.items():
        if alias in case_names:
            raise ValueError(f"Stage I analysis alias shadows a canonical case: {alias}")
        if canonical not in case_names:
            raise ValueError(
                f"Stage I analysis alias selects unknown canonical case: {canonical}"
            )
    if not isinstance(reference_manifests, dict) or not reference_manifests:
        raise ValueError("Stage I reference_manifests must be a nonempty object")
    normalized_manifests: dict[str, dict[str, str]] = {}
    for name, record in reference_manifests.items():
        if not isinstance(name, str) or not name.strip() or not isinstance(record, dict):
            raise ValueError("Stage I reference_manifests entries must be named objects")
        normalized_manifests[name] = {
            "path": require_manifest_text(record, "path", f"Stage I reference {name}"),
            "sha256": require_sha256(record, "sha256", f"Stage I reference {name}"),
        }
    if not isinstance(reference_product_bindings, dict):
        raise ValueError("Stage I reference_product_bindings must be an object")
    normalized_bindings: dict[str, dict[str, str]] = {}
    for product_id, binding in reference_product_bindings.items():
        context = f"Stage I reference product {product_id}"
        if (
            not isinstance(product_id, str) or not product_id.strip()
            or not isinstance(binding, dict)
        ):
            raise ValueError("Stage I reference_product_bindings entries are malformed")
        kind = require_manifest_text(binding, "kind", context)
        if kind not in ("curve", "surface"):
            raise ValueError(f"{context} has unsupported kind: {kind}")
        reference_manifest = require_manifest_text(
            binding, "reference_manifest", context
        )
        if reference_manifest not in normalized_manifests:
            raise ValueError(
                f"{context} selects unknown reference manifest: {reference_manifest}"
            )
        normalized_bindings[product_id] = {
            "kind": kind,
            "case": require_manifest_text(binding, "case", context),
            "product": require_manifest_text(binding, "product", context),
            "data_file": require_manifest_text(binding, "data_file", context),
            "data_sha256": require_sha256(binding, "data_sha256", context),
            "reference_manifest": reference_manifest,
            "reference_manifest_sha256": normalized_manifests[
                reference_manifest
            ]["sha256"],
        }
    normalized: list[dict[str, object]] = []
    panel_ids: set[str] = set()
    configured_products: set[str] = set()
    for panel in panels:
        panel_id = require_manifest_text(panel, "id", "Stage I panel")
        context = f"Stage I panel {panel_id}"
        if panel_id in panel_ids:
            raise ValueError(f"Stage I panel id is duplicated: {panel_id}")
        panel_ids.add(panel_id)
        disposition = require_manifest_text(panel, "disposition", context)
        if disposition not in ("comparison", "blocked_reference", "external_model"):
            raise ValueError(f"{context} has unsupported disposition: {disposition}")
        required_cases = require_text_list(
            panel, "required_cases", context, required=disposition == "comparison"
        )
        unknown_cases = sorted(set(required_cases) - case_ids)
        if unknown_cases:
            raise ValueError(f"{context} selects unknown cases: {unknown_cases}")
        products = require_text_list(
            panel, "reference_products", context,
            required=disposition == "comparison",
        )
        configured_products.update(products)
        criteria = panel.get("criteria", [])
        if not isinstance(criteria, list) or not all(
            isinstance(criterion, dict) for criterion in criteria
        ):
            raise ValueError(f"{context} criteria must be a list of objects")
        if disposition != "comparison" and criteria:
            raise ValueError(f"{context} static disposition cannot define criteria")
        criterion_state: str | None = None
        criterion_reason = str(panel.get("criterion_reason", "")).strip()
        if disposition == "comparison":
            criterion_state = require_manifest_text(panel, "criterion_state", context)
            if criterion_state not in ("pending_review", "reviewed"):
                raise ValueError(
                    f"{context} has unsupported criterion_state: {criterion_state}"
                )
            if criterion_state == "pending_review" and not criterion_reason:
                raise ValueError(
                    f"{context} pending_review requires criterion_reason"
                )
            if criterion_state == "reviewed" and not criteria:
                raise ValueError(f"{context} reviewed criteria cannot be empty")
        elif "criterion_state" in panel or criterion_reason:
            raise ValueError(
                f"{context} static disposition cannot define criterion lifecycle"
            )
        normalized_criteria: list[dict[str, object]] = []
        criterion_products: set[str] = set()
        for criterion in criteria:
            product = require_manifest_text(criterion, "product", context)
            metric = require_manifest_text(criterion, "metric", context)
            operator = require_manifest_text(criterion, "operator", context)
            if operator not in ("<", "<=", "==", ">=", ">"):
                raise ValueError(
                    f"{context} criterion has unsupported operator: {operator}"
                )
            try:
                limit = float(criterion["limit"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"{context} criterion requires numeric limit") from error
            if not math.isfinite(limit):
                raise ValueError(f"{context} criterion limit must be finite")
            if product not in products:
                raise ValueError(
                    f"{context} criterion selects unlisted product: {product}"
                )
            criterion_products.add(product)
            normalized_criteria.append({
                "product": product,
                "metric": metric,
                "operator": operator,
                "limit": limit,
            })
        if criterion_state == "reviewed" and criterion_products != set(products):
            raise ValueError(
                f"{context} requires reviewed criteria for every reference product"
            )
        reason = str(panel.get("reason", "")).strip()
        if disposition != "comparison" and not reason:
            raise ValueError(f"{context} static disposition requires reason")
        normalized.append({
            "id": panel_id,
            "figure": str(panel.get("figure", "")),
            "panel": str(panel.get("panel", "")),
            "description": str(panel.get("description", "")),
            "disposition": disposition,
            "required_cases": required_cases,
            "reference_products": products,
            "criteria": normalized_criteria,
            "criterion_state": criterion_state,
            "criterion_reason": criterion_reason,
            "configured_reason": reason,
        })
    if set(normalized_bindings) != configured_products:
        raise ValueError(
            "Stage I reference_product_bindings must exactly cover configured "
            "comparison products"
        )
    return {
        "schema_version": STAGE_I_PANEL_SCHEMA_VERSION,
        "manifest": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
        "case_names": case_names,
        "case_ids": sorted(case_ids),
        "analysis_case_aliases": analysis_case_aliases,
        "reference_manifests": normalized_manifests,
        "reference_product_bindings": normalized_bindings,
        "panels": normalized,
    }


def stage_i_bundle_case_ids(bundle: Path, manifest: dict[str, object],
                            configuration: dict[str, object]) -> list[str]:
    """Map accepted bundle cases to canonical Stage I ids."""

    admission = manifest.get("stage_i_admission_status")
    if admission is None and (
        manifest.get("workflow") == STAGE_I_PRODUCTION_WORKFLOW
        and manifest.get("status") == "accepted_for_analysis"
    ):
        # Explicit legacy accepted status is sufficient migration evidence.
        admission = "accepted_for_analysis"
    if admission != "accepted_for_analysis":
        raise ValueError(
            "Stage I panel status table requires an accepted_for_analysis bundle"
        )
    case_names = configuration["case_names"]
    if not isinstance(case_names, dict):
        raise ValueError("Stage I panel configuration case_names must be an object")
    cases = bundle_cases(bundle, manifest)
    accepted: list[str] = []
    for case in cases:
        case_id = case.get("case_id")
        if case_id is None and len(cases) == 1:
            case_id = manifest.get("production_case_id")
        if case_id is None:
            case_id = case_names.get(str(case.get("name", "")))
        if not isinstance(case_id, str) or case_id not in configuration["case_ids"]:
            raise ValueError(
                f"accepted bundle case does not map to Stage I inventory: "
                f"{case.get('name', 'unnamed_case')}"
            )
        accepted.append(case_id)
    if len(set(accepted)) != len(accepted):
        raise ValueError("accepted bundle maps more than one case to a Stage I id")
    return sorted(accepted)


def criterion_passes(observed: float, operator: str, limit: float) -> bool:
    """Apply one reviewed generic criterion without embedding scientific limits."""

    return {
        "<": observed < limit,
        "<=": observed <= limit,
        "==": observed == limit,
        ">=": observed >= limit,
        ">": observed > limit,
    }[operator]


def stage_i_panel_status_table(
    bundle: Path, bundle_metadata: dict[str, object],
    configuration: dict[str, object],
    comparisons: dict[str, object] | None,
) -> dict[str, object]:
    """Build one retained Stage I status row per configured paper panel."""

    accepted_cases = stage_i_bundle_case_ids(bundle, bundle_metadata, configuration)
    accepted = set(accepted_cases)
    available_products: dict[str, object] = {}
    omitted_products: dict[str, object] = {}
    if comparisons is not None:
        for collection in ("comparisons", "surface_comparisons"):
            products = comparisons.get(collection, {})
            if not isinstance(products, dict):
                raise ValueError(f"reference comparison {collection} is not an object")
            available_products.update(products)
        omissions = comparisons.get("omitted_products", [])
        if not isinstance(omissions, list):
            raise ValueError("reference comparison omitted_products is not a list")
        for omission in omissions:
            if not isinstance(omission, dict) or not isinstance(
                omission.get("id"), str
            ):
                raise ValueError("reference comparison omission is malformed")
            omitted_products[str(omission["id"])] = omission
    rows: list[dict[str, object]] = []
    counts = {status: 0 for status in STAGE_I_PANEL_STATUSES}
    for configured in configuration["panels"]:
        row = dict(configured)
        required_cases = set(row["required_cases"])
        required_products = set(row["reference_products"])
        missing_cases = sorted(required_cases - accepted)
        missing_products = sorted(required_products - set(available_products))
        omitted = sorted(required_products & set(omitted_products))
        row.update({
            "missing_cases": missing_cases,
            "available_reference_products": sorted(
                required_products & set(available_products)
            ),
            "missing_reference_products": missing_products,
            "omitted_reference_products": omitted,
            "criterion_results": [],
        })
        disposition = row["disposition"]
        if disposition == "external_model":
            status = "external_model"
            reason = row["configured_reason"]
        elif disposition == "blocked_reference":
            status = "blocked_reference"
            reason = row["configured_reason"]
        elif row["criterion_state"] == "pending_review":
            status = "not_run"
            reason = row["criterion_reason"]
        elif missing_cases or missing_products:
            status = "not_run"
            details = []
            if missing_cases:
                details.append(f"missing accepted cases: {', '.join(missing_cases)}")
            if missing_products:
                details.append(
                    f"missing reference products: {', '.join(missing_products)}"
                )
            reason = "; ".join(details)
        else:
            criterion_results = []
            for criterion in row["criteria"]:
                product = str(criterion["product"])
                record = available_products[product]
                if not isinstance(record, dict):
                    raise ValueError(f"reference product is not an object: {product}")
                if not record.get("stage_i_binding_validated", False):
                    raise ValueError(
                        f"reviewed Stage I reference product is not binding-validated: "
                        f"{product}"
                    )
                validate_stage_i_reference_binding(
                    configuration["reference_product_bindings"], product,
                    str(record.get("kind", "")), str(record.get("case", "")),
                    str(record.get("product", "")),
                    str(record.get("reference_data_file", "")),
                    str(record.get("data_sha256", "")),
                    str(record.get("reference_manifest_sha256", "")),
                )
                observed = record.get(str(criterion["metric"]))
                if isinstance(observed, bool) or not isinstance(
                    observed, (int, float)
                ) or not math.isfinite(float(observed)):
                    raise ValueError(
                        f"reference product {product} requires finite metric "
                        f"{criterion['metric']}"
                    )
                passed = criterion_passes(
                    float(observed), str(criterion["operator"]),
                    float(criterion["limit"]),
                )
                criterion_results.append({
                    **criterion,
                    "observed": float(observed),
                    "passed": passed,
                })
            row["criterion_results"] = criterion_results
            if all(bool(criterion["passed"]) for criterion in criterion_results):
                status = "passed"
                reason = "all reviewed criteria passed"
            else:
                status = "failed"
                reason = "one or more reviewed criteria failed"
        row.update({"status": status, "reason": reason})
        counts[status] += 1
        rows.append(row)
    return {
        "schema_version": STAGE_I_PANEL_SCHEMA_VERSION,
        "definition": (
            "per-panel Stage I comparison gate from accepted bundles, "
            "tracked dispositions, and reviewed data criteria"
        ),
        "stage_i_manifest": configuration["manifest"],
        "stage_i_manifest_sha256": configuration["manifest_sha256"],
        "bundle": str(bundle),
        "bundle_manifest": str(bundle / "manifest.json"),
        "bundle_manifest_sha256_at_analysis": sha256_file(bundle / "manifest.json"),
        "stage_i_admission_status": "accepted_for_analysis",
        "accepted_cases": accepted_cases,
        "analysis_case_aliases": configuration["analysis_case_aliases"],
        "status_counts": counts,
        "panels": rows,
    }


def markdown_cell(value: object) -> str:
    """Escape one compact Markdown table cell."""

    return str(value).replace("|", r"\|").replace("\n", " ")


def write_stage_i_panel_status_table(table: dict[str, object],
                                     output_dir: Path) -> None:
    """Retain machine-readable and concise human-readable Stage I gate tables."""

    json_path = output_dir / "stage_i_panel_status.json"
    markdown_path = output_dir / "stage_i_panel_status.md"
    json_path.write_text(json.dumps(table, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    lines = [
        "# CGL-LF MKS24 Stage I Panel Status",
        "",
        f"- Stage I manifest: `{table['stage_i_manifest']}`",
        f"- Accepted cases: `{', '.join(table['accepted_cases'])}`",
        "",
        "| Panel | Figure | Status | Required cases | Reference products | Reason |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for panel in table["panels"]:
        figure = str(panel["figure"])
        if panel["panel"]:
            figure += str(panel["panel"])
        lines.append(
            "| `{}` | `{}` | **{}** | `{}` | `{}` | {} |".format(
                markdown_cell(panel["id"]),
                markdown_cell(figure),
                markdown_cell(panel["status"]),
                markdown_cell(", ".join(panel["required_cases"])),
                markdown_cell(", ".join(panel["reference_products"])),
                markdown_cell(panel["reason"]),
            )
        )
    markdown_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    """Command-line entry point."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("snapshots", nargs="*", type=Path)
    command.add_argument("--bundle", type=Path)
    command.add_argument("--history", action="append", type=Path, default=[])
    command.add_argument("--lf-history", action="append", type=Path, default=[])
    command.add_argument("--output-dir", type=Path, default=Path("cgl_lf_paper_analysis"))
    command.add_argument("--pdf-bins", type=int, default=64)
    command.add_argument("--alignment-shells", default="1,2,3")
    command.add_argument("--eddy-samples", type=int, default=0)
    command.add_argument("--eddy-bins", type=int, default=20)
    command.add_argument("--eddy-seed", type=int, default=0)
    command.add_argument("--time-start", type=float)
    command.add_argument("--time-end", type=float)
    command.add_argument("--reference-curves", action="append", type=Path, default=[])
    command.add_argument(
        "--stage-i-manifest",
        type=Path,
        help=(
            "Tracked Stage I case and versioned panel configuration used to "
            "retain a per-panel status table for accepted bundles."
        ),
    )
    command.add_argument(
        "--allow-partial-reference-cases",
        action="store_true",
        help=(
            "Compare available bundle cases while retaining explicit omissions "
            "for other cases named by reference manifests."
        ),
    )
    command.add_argument(
        "--figure-13-only",
        action="store_true",
        help=(
            "Stream only the exact Figure 13 occupancy product from production "
            "Athena snapshots without constructing unrelated full-grid diagnostics."
        ),
    )
    command.add_argument("--synthetic-test", action="store_true")
    args = command.parse_args()
    if (
        args.time_start is not None
        and args.time_end is not None
        and args.time_end < args.time_start
    ):
        command.error("--time-end must be greater than or equal to --time-start")
    if args.eddy_samples < 0:
        command.error("--eddy-samples must be nonnegative")
    if args.eddy_bins < 2:
        command.error("--eddy-bins must be at least 2")
    if args.stage_i_manifest is not None and args.bundle is None:
        command.error("--stage-i-manifest requires --bundle")
    if args.figure_13_only and args.bundle is None:
        command.error("--figure-13-only requires --bundle")

    alignment_shells = [
        int(value) for value in args.alignment_shells.split(",") if value.strip()
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    analysis_configuration: dict[str, object] = {
        "snapshots": [str(path) for path in args.snapshots],
        "bundle": str(args.bundle) if args.bundle is not None else None,
        "history": [str(path) for path in args.history],
        "lf_history": [str(path) for path in args.lf_history],
        "output_dir": str(args.output_dir),
        "pdf_bins": args.pdf_bins,
        "alignment_shells": alignment_shells,
        "eddy_samples": args.eddy_samples,
        "eddy_bins": args.eddy_bins,
        "eddy_seed": args.eddy_seed,
        "time_start": args.time_start,
        "time_end": args.time_end,
        "reference_curves": [str(path) for path in args.reference_curves],
        "stage_i_manifest": (
            str(args.stage_i_manifest)
            if args.stage_i_manifest is not None else None
        ),
        "allow_partial_reference_cases": args.allow_partial_reference_cases,
        "figure_13_only": args.figure_13_only,
        "synthetic_test": args.synthetic_test,
    }
    analysis_provenance = analysis_invocation_provenance(analysis_configuration)
    result: dict[str, object] = {
        "histories": [],
        "lf_histories": [],
        "forcing_energy_budgets": [],
        "snapshots": {},
    }
    bundle_metadata: dict[str, object] | None = None
    retained_bundle_cases: list[dict[str, object]] | None = None
    figure_13_product: dict[str, object] | None = None
    figure_13_execution_authentication_record: dict[str, object] | None = None
    figure_13_provenance_cache: dict[str, dict[str, object]] = {}
    figure_13_occupancy_cache: dict[str, dict[str, object]] = {}
    figure_13_revalidated_cache: set[str] = set()
    if args.bundle is not None:
        bundle_metadata = bundle_manifest(args.bundle)
        retained_bundle_cases = bundle_cases(args.bundle, bundle_metadata)
        figure_13_execution_authentication_record = (
            figure_13_execution_authentication(
                args.bundle, bundle_metadata, retained_bundle_cases
            )
        )
        expected_ranks_by_path = figure_13_execution_authentication_record[
            "snapshot_expected_ranks"
        ]
        if not isinstance(expected_ranks_by_path, dict):
            raise ValueError("Figure 13 execution authentication rank map is malformed")
        result["cases"] = {}
        for case in retained_bundle_cases:
            name = str(case.get("name", "unnamed_case"))
            outputs = case.get("outputs", {})
            if not isinstance(outputs, dict):
                outputs = {}
            start, end = case_window(case, args.time_start, args.time_end)
            histories = []
            if "user_history" in outputs:
                histories.append(args.bundle / str(outputs["user_history"]))
            lf_histories = []
            if "mhd_history" in outputs:
                lf_histories.append(args.bundle / str(outputs["mhd_history"]))
            snapshots = [
                args.bundle / str(path)
                for path in outputs.get("snapshot_paths", [])
            ]
            model = case.get("model_choices", {})
            if not isinstance(model, dict):
                model = {}
            if args.figure_13_only:
                records, ensemble = analyze_figure_13_snapshot_paths(
                    snapshots,
                    start,
                    end,
                    expected_ranks_by_path,
                    figure_13_provenance_cache,
                    figure_13_occupancy_cache,
                    figure_13_revalidated_cache,
                )
            else:
                records, ensemble = analyze_snapshot_paths(
                    snapshots, args.pdf_bins, alignment_shells, start, end, model,
                    args.eddy_samples, args.eddy_bins, args.eddy_seed,
                    expected_ranks_by_path,
                )
            summaries = [summarize_history(path, start, end) for path in histories]
            lf_summaries = [
                summarize_lf_history(path, start, end) for path in lf_histories
            ]
            budgets = [
                summarize_forcing_energy_budget(
                    histories[0], lf_histories[0], model, start, end
                )
            ] if histories and lf_histories else []
            result["histories"].extend(summaries)
            result["lf_histories"].extend(lf_summaries)
            result["forcing_energy_budgets"].extend(budgets)
            result["snapshots"].update(records)
            result["cases"][name] = {
                "analysis_window": {"time_start": start, "time_end": end},
                "histories": summaries,
                "lf_histories": lf_summaries,
                "forcing_energy_budgets": budgets,
                "snapshots": records,
                "snapshot_ensemble": ensemble,
            }
        figure_13_product = figure_13_alternate_firehose_occupancy(
            args.bundle,
            bundle_metadata,
            retained_bundle_cases,
            result["cases"],
            analysis_provenance,
            figure_13_execution_authentication_record,
        )
        if figure_13_product is not None:
            result["figure_13_alternate_firehose_occupancy"] = figure_13_product
    direct_histories = [
        summarize_history(path, args.time_start, args.time_end)
        for path in args.history
    ]
    result["histories"].extend(direct_histories)
    result["lf_histories"].extend(
        summarize_lf_history(path, args.time_start, args.time_end)
        for path in args.lf_history
    )
    direct_records, direct_ensemble = analyze_snapshot_paths(
        list(args.snapshots), args.pdf_bins, alignment_shells,
        args.time_start, args.time_end, None, args.eddy_samples,
        args.eddy_bins, args.eddy_seed
    )
    result["snapshots"].update(direct_records)
    if args.snapshots:
        result["snapshot_ensemble"] = direct_ensemble
    if args.bundle is None and not args.snapshots and args.history:
        result["analysis_window"] = {
            "time_start": args.time_start,
            "time_end": args.time_end,
        }
    if args.synthetic_test:
        result["synthetic_test"] = synthetic_test()
    configuration: dict[str, object] | None = None
    if args.stage_i_manifest is not None:
        configuration = stage_i_panels_configuration(args.stage_i_manifest)
    if args.reference_curves:
        result["reference_curve_comparisons"] = combined_reference_curve_comparisons(
            result, args.reference_curves,
            allow_missing_cases=args.allow_partial_reference_cases,
            analysis_case_aliases=(
                configuration["analysis_case_aliases"]
                if configuration is not None else None
            ),
            stage_i_reference_bindings=(
                configuration["reference_product_bindings"]
                if configuration is not None else None
            ),
        )
    if args.stage_i_manifest is not None:
        assert args.bundle is not None and bundle_metadata is not None
        assert configuration is not None
        result["stage_i_panel_status"] = stage_i_panel_status_table(
            args.bundle, bundle_metadata, configuration,
            result.get("reference_curve_comparisons"),
        )
        write_stage_i_panel_status_table(
            result["stage_i_panel_status"], args.output_dir
        )
    if figure_13_product is not None:
        write_figure_13_publication_generation(result, args.output_dir)
    else:
        write_analysis_publication_generation(result, args.output_dir)
    destination = read_analysis_publication(args.output_dir)["diagnostics"]
    print(f"Wrote {destination}")
    if args.synthetic_test and not result["synthetic_test"]["passed"]:
        print("CGL-LF paper analysis synthetic test failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
