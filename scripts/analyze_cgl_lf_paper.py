#!/usr/bin/env python3
"""Generate MKS24-oriented diagnostics from archived CGL-LF snapshots."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR / "vis" / "python"))
import bin_convert  # noqa: E402


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
EDDY_ANGLE_DEGREES = 15.0


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
    if firehose_policy == "oblique":
        firehose_threshold = -0.7
    elif firehose_policy == "parallel":
        firehose_threshold = -1.0
    else:
        raise ValueError(f"invalid cgl_firehose_threshold={firehose_policy}")
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
    return {
        "grad_parallel_velocity_parallel": parallel_parallel,
        "grad_perp_velocity_parallel": perpendicular_parallel,
        "grad_parallel_velocity_perp": parallel_perp,
        "grad_perp_velocity_perp": perpendicular_perp,
        "bb_grad_velocity": strain_parallel,
    }


def pressure_work_decomposition(fields: dict[str, np.ndarray],
                                lengths: tuple[float, float, float],
                                strain_parallel: np.ndarray,
                                model: dict[str, object] | None
                                ) -> dict[str, object]:
    """Reconstruct cell-centered CGL pressure-force work from one snapshot.

    For P = p_perp I - Delta p b b, periodic integration of the pressure
    force contribution to kinetic energy is integral[P : grad(u)] dV.
    """

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
              ranges: tuple[tuple[float, float], tuple[float, float]] | None = None
              ) -> dict[str, object]:
    """Return a density-normalized two-dimensional histogram."""

    density, x_edges, y_edges = np.histogram2d(
        x_values.ravel(), y_values.ravel(), bins=bins, range=ranges, density=True
    )
    return {
        "x_edges": x_edges.tolist(),
        "y_edges": y_edges.tolist(),
        "density": density.tolist(),
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


def read_snapshot(path: Path
                  ) -> tuple[dict[str, np.ndarray], tuple[float, float, float], float]:
    """Read one current CGL primitive binary snapshot."""

    rank_local = path.parent.name == "rank_00000000"
    raw_reader = (
        bin_convert.read_all_ranks_binary if rank_local
        else bin_convert.read_binary
    )
    athdf_reader = (
        bin_convert.read_all_ranks_binary_as_athdf if rank_local
        else bin_convert.read_binary_as_athdf
    )
    raw = raw_reader(str(path))
    missing = sorted(set(REQUIRED_FIELDS) - set(raw["var_names"]))
    if missing:
        raise ValueError(
            f"{path} lacks CGL paper fields {missing}; regenerate with current "
            "mhd_w_bcc output (legacy eint is p_parallel)"
        )
    values = athdf_reader(
        str(path), quantities=list(REQUIRED_FIELDS), dtype=np.float64
    )
    lengths = (
        float(raw["x1max"] - raw["x1min"]),
        float(raw["x2max"] - raw["x2min"]),
        float(raw["x3max"] - raw["x3min"]),
    )
    return values, lengths, float(raw["time"])


def pdf_fields(fields: dict[str, np.ndarray],
               lengths: tuple[float, float, float] | None = None
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
        "beta_delta": 2.0 * (pperp - ppar)
        / np.maximum(bsqr, np.finfo(float).tiny),
    }
    if lengths is not None:
        magnetic = [fields["bcc1"], fields["bcc2"], fields["bcc3"]]
        velocity = [fields["velx"], fields["vely"], fields["velz"]]
        bhat = [component / np.sqrt(np.maximum(bsqr, np.finfo(float).tiny))
                for component in magnetic]
        values["bb_grad_velocity"] = velocity_gradient_products(
            velocity, bhat, lengths
        )["bb_grad_velocity"]
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
    transfer = pressure_transfer(rho, velocity, magnetic, delta_p, lengths, dk)
    pressure_work = pressure_work_decomposition(
        fields, lengths, velocity_products["bb_grad_velocity"], model_choices
    )
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
    pdf_values = pdf_fields(fields, lengths)
    pressure_density_values = pressure_density_fields(fields)
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
            [delta_p], lengths, dk, field_definition="pressure anisotropy Delta p"
        ),
        "grad_parallel_delta_p": shell_spectrum([gradient_parallel], lengths, dk),
        "grad_perp_delta_p": shell_spectrum([gradient_perp], lengths, dk),
    }
    spectra.update({
        name: shell_spectrum([values], lengths, dk)
        for name, values in velocity_products.items()
    })
    return {
        "time": time,
        "shape_z_y_x": list(rho.shape),
        "lengths_x_y_z": list(lengths),
        "pdf": {
            name: pdf(values, bins, None if pdf_ranges is None else pdf_ranges[name])
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
        "pressure_work_decomposition": pressure_work,
        "alignment": alignment_histograms(velocity, magnetic, lengths, dk,
                                          alignment_shells, bins),
        "eddy_anisotropy": local_field_eddy_anisotropy(
            velocity, magnetic, lengths, eddy_samples, eddy_bins, eddy_seed
        ),
        "heat_flux_transport_proxy": heat_flux_transport_proxy(
            fields, lengths, model_choices
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


def mean_joint_distribution(records: list[dict[str, object]]) -> dict[str, object]:
    """Average compatible two-dimensional histogram products."""

    return {
        "x_edges": records[0]["x_edges"],
        "y_edges": records[0]["y_edges"],
        "density": np.mean(
            [np.asarray(record["density"], dtype=float) for record in records], axis=0
        ).tolist(),
    }


def mean_spectrum(records: list[dict[str, object]]) -> dict[str, object]:
    """Average compatible shell-summed spectra."""

    result: dict[str, object] = {
        "dk": records[0]["dk"],
        "perpendicular": records[0]["perpendicular"],
        "k": records[0]["k"],
        "power_per_dk": np.mean(
            [np.asarray(record["power_per_dk"], dtype=float) for record in records],
            axis=0,
        ).tolist(),
    }
    for key in ("field_definition", "normalization_definition"):
        if key in records[0]:
            result[key] = records[0][key]
    return result


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
            pressure_integral[f"{name}_integral"] = float(np.trapz(
                [sample[name] for sample in pressure_work], times
            ))
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
                heat_flux_integral[f"{name}_integral"] = float(np.trapz(
                    [sample[name] for sample in heat_flux], heat_flux_times
                ))
        heat_flux_ensemble["time_integral_estimate"] = heat_flux_integral
    return {
        "snapshot_count": len(samples),
        "time_first": min(float(sample["time"]) for sample in samples),
        "time_last": max(float(sample["time"]) for sample in samples),
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
                ])
                for name in joint_names
            },
        },
        "spectra": {
            name: mean_spectrum([sample["spectra"][name] for sample in samples])
            for name in spectrum_names
        },
        "pressure_transfer": {
            "dk": transfer[0]["dk"],
            "k_perp": transfer[0]["k_perp"],
            "transfer": np.mean(
                [np.asarray(item["transfer"], dtype=float) for item in transfer], axis=0
            ).tolist(),
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
            "direct_real_space_mean": float(np.mean(
                [item["direct_real_space"] for item in transfer]
            )),
            "shell_sum_mean": float(np.mean([item["shell_sum"] for item in transfer])),
            "closure_error_mean": float(np.mean(
                [item["closure_error"] for item in transfer]
            )),
        },
        "pressure_work_decomposition": pressure_work_ensemble,
        "alignment": {
            name: mean_distribution([sample["alignment"][name] for sample in samples])
            for name in sorted(alignment_names)
        },
        "eddy_anisotropy": mean_eddy_anisotropy([
            sample["eddy_anisotropy"] for sample in samples
        ]),
        "heat_flux_transport_proxy": heat_flux_ensemble,
    }


def analyze_snapshot_paths(paths: list[Path], bins: int, alignment_shells: list[int],
                           time_start: float | None = None,
                           time_end: float | None = None,
                           model_choices: dict[str, object] | None = None,
                           eddy_samples: int = 0, eddy_bins: int = 20,
                           eddy_seed: int = 0
                           ) -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    """Analyze and time-average selected snapshots with common PDF bin edges."""

    selected: list[Path] = []
    extrema: dict[str, list[float]] = {}
    joint_extrema: dict[str, list[list[float]]] = {}
    for path in paths:
        fields, lengths, time = read_snapshot(path)
        if not time_mask(np.asarray([time]), time_start, time_end)[0]:
            continue
        selected.append(path)
        for name, values in pdf_fields(fields, lengths).items():
            finite = values[np.isfinite(values)]
            if finite.size == 0:
                continue
            low, high = float(np.min(finite)), float(np.max(finite))
            if name in extrema:
                extrema[name][0] = min(extrema[name][0], low)
                extrema[name][1] = max(extrema[name][1], high)
            else:
                extrema[name] = [low, high]
        for name, (x_values, y_values) in pressure_density_fields(fields).items():
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
        fields, lengths, time = read_snapshot(path)
        records[str(path)] = analyze_fields(
            fields, lengths, time, bins, alignment_shells, ranges, model_choices,
            eddy_samples, eddy_bins, eddy_seed, joint_ranges
        )
    ensemble = average_snapshot_records(records)
    ensemble["time_start"] = time_start
    ensemble["time_end"] = time_end
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


def bundle_cases(bundle: Path) -> list[dict[str, object]]:
    """Read case metadata used to preserve per-case analysis windows."""

    with (bundle / "manifest.json").open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    cases = manifest.get("cases", [])
    if not isinstance(cases, list):
        raise ValueError(f"manifest cases must be a list: {bundle}")
    return cases


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of one reference-data input."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_manifest_text(record: dict[str, object], key: str, context: str) -> str:
    """Return one required nonempty manifest text value."""

    value = record.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context} requires nonempty {key}")
    return value


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
    if product == "pressure_transfer.transfer":
        record = ensemble.get("pressure_transfer", {})
        if not isinstance(record, dict):
            raise ValueError(f"analyzed product is missing: {product}")
        return (
            np.asarray(record["k_perp"], dtype=float),
            np.asarray(record["transfer"], dtype=float),
        )
    if product == "pressure_transfer.transfer_normalized_by_total":
        record = ensemble.get("pressure_transfer", {})
        if (
            not isinstance(record, dict)
            or not record.get("normalization_available", False)
            or record.get("transfer_normalized_by_total") is None
        ):
            raise ValueError(f"analyzed product is missing: {product}")
        return (
            np.asarray(record["k_perp"], dtype=float),
            np.asarray(record["transfer_normalized_by_total"], dtype=float),
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


def reference_curve_comparisons(result: dict[str, object],
                                manifest_path: Path) -> dict[str, object]:
    """Compare analysis products against provenance-qualified external data."""

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
    for curve in curves:
        if not isinstance(curve, dict):
            raise ValueError("reference curve entries must be objects")
        curve_id = require_manifest_text(curve, "id", "reference curve")
        if curve_id in comparisons:
            raise ValueError(f"reference curve id is duplicated: {curve_id}")
        case = require_manifest_text(curve, "case", f"reference curve {curve_id}")
        product = require_manifest_text(curve, "product", f"reference curve {curve_id}")
        data_file = require_manifest_text(
            curve, "data_file", f"reference curve {curve_id}"
        )
        expected_digest = require_manifest_text(
            curve, "data_sha256", f"reference curve {curve_id}"
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
        if case == "direct":
            ensemble = result.get("snapshot_ensemble", {})
            histories = result.get("histories", [])
        elif isinstance(cases, dict) and case in cases:
            ensemble = cases[case].get("snapshot_ensemble", {})
            histories = cases[case].get("histories", [])
        else:
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
            "case": case,
            "product": product,
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
        if surface_id in comparisons or surface_id in surface_comparisons:
            raise ValueError(f"reference product id is duplicated: {surface_id}")
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
        if case == "direct":
            ensemble = result.get("snapshot_ensemble", {})
        elif isinstance(cases, dict) and case in cases:
            ensemble = cases[case].get("snapshot_ensemble", {})
        else:
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
            "case": case,
            "product": product,
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
        "manifest_sha256": sha256_file(manifest_path),
        "provenance": provenance,
        "comparisons": comparisons,
        "surface_comparisons": surface_comparisons,
    }


def combined_reference_curve_comparisons(
    result: dict[str, object], manifest_paths: list[Path],
) -> dict[str, object]:
    """Combine distinct comparison products from qualified reference manifests."""

    manifests: list[dict[str, object]] = []
    comparisons: dict[str, object] = {}
    surface_comparisons: dict[str, object] = {}
    combined: dict[str, object] = {
        "available": True,
        "definition": (
            "analysis products interpolated onto provenance-qualified reference "
            "curve or surface coordinates"
        ),
        "manifests": manifests,
        "comparisons": comparisons,
        "surface_comparisons": surface_comparisons,
    }
    product_ids: set[str] = set()
    for manifest_path in manifest_paths:
        comparison = reference_curve_comparisons(result, manifest_path)
        manifests.append({
            "manifest": comparison["manifest"],
            "manifest_sha256": comparison["manifest_sha256"],
            "provenance": comparison["provenance"],
        })
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
                if product_id in product_ids:
                    raise ValueError(
                        f"reference product id is duplicated across manifests: "
                        f"{product_id}"
                    )
                retained[product_id] = product
                product_ids.add(product_id)
    if len(manifest_paths) == 1:
        only = manifests[0]
        combined.update({
            "manifest": only["manifest"],
            "manifest_sha256": only["manifest_sha256"],
            "provenance": only["provenance"],
        })
    return combined


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

    alignment_shells = [
        int(value) for value in args.alignment_shells.split(",") if value.strip()
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, object] = {
        "histories": [],
        "lf_histories": [],
        "forcing_energy_budgets": [],
        "snapshots": {},
    }
    if args.bundle is not None:
        result["cases"] = {}
        for case in bundle_cases(args.bundle):
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
            records, ensemble = analyze_snapshot_paths(
                snapshots, args.pdf_bins, alignment_shells, start, end, model,
                args.eddy_samples, args.eddy_bins, args.eddy_seed
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
    if args.reference_curves:
        result["reference_curve_comparisons"] = combined_reference_curve_comparisons(
            result, args.reference_curves
        )
    destination = args.output_dir / "diagnostics.json"
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(f"Wrote {destination}")
    if args.synthetic_test and not result["synthetic_test"]["passed"]:
        print("CGL-LF paper analysis synthetic test failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
