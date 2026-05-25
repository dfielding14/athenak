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
                   dk: float, perpendicular: bool = True) -> dict[str, object]:
    """Compute an MKS24-style bin-summed spectrum for scalar or vector fields."""

    shell = shell_indices(fields[0].shape, lengths, dk, perpendicular)
    power = np.zeros(fields[0].shape, dtype=float)
    normalizer = float(fields[0].size)
    for field in fields:
        transformed = np.fft.fftn(field - np.mean(field)) / normalizer
        power += np.abs(transformed) ** 2
    binned = np.bincount(shell.ravel(), weights=power.ravel())
    return {
        "dk": dk,
        "perpendicular": perpendicular,
        "k": (np.arange(len(binned), dtype=float) * dk).tolist(),
        "power_per_dk": (binned / dk).tolist(),
    }


def pdf(values: np.ndarray, bins: int, value_range: tuple[float, float] | None = None
        ) -> dict[str, object]:
    """Return a density-normalized histogram."""

    counts, edges = np.histogram(
        values.ravel(), bins=bins, range=value_range, density=True
    )
    return {"edges": edges.tolist(), "density": counts.tolist()}


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
    """

    bsqr = sum(component * component for component in magnetic)
    safe_bsqr = np.maximum(bsqr, np.finfo(float).tiny)
    safe_rho = np.maximum(rho, np.finfo(float).tiny)
    root_rho = np.sqrt(safe_rho)
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
    return {
        "dk": dk,
        "k_perp": (np.arange(len(values), dtype=float) * dk).tolist(),
        "transfer": values,
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


def read_snapshot(path: Path
                  ) -> tuple[dict[str, np.ndarray], tuple[float, float, float], float]:
    """Read one current CGL primitive binary snapshot."""

    raw = bin_convert.read_binary(str(path))
    missing = sorted(set(REQUIRED_FIELDS) - set(raw["var_names"]))
    if missing:
        raise ValueError(
            f"{path} lacks CGL paper fields {missing}; regenerate with current "
            "mhd_w_bcc output (legacy eint is p_parallel)"
        )
    values = bin_convert.read_binary_as_athdf(
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


def analyze_fields(fields: dict[str, np.ndarray], lengths: tuple[float, float, float],
                   time: float, bins: int, alignment_shells: list[int],
                   pdf_ranges: dict[str, tuple[float, float]] | None = None,
                   model_choices: dict[str, object] | None = None
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
    spectra = {
        "velocity": shell_spectrum(velocity, lengths, dk),
        "magnetic_fluctuation": shell_spectrum(magnetic, lengths, dk),
        "density": shell_spectrum([rho], lengths, dk),
        "delta_p": shell_spectrum([delta_p], lengths, dk),
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
        "spectra": spectra,
        "pressure_transfer": transfer,
        "pressure_work_decomposition": pressure_work,
        "alignment": alignment_histograms(velocity, magnetic, lengths, dk,
                                          alignment_shells, bins),
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


def mean_spectrum(records: list[dict[str, object]]) -> dict[str, object]:
    """Average compatible shell-summed spectra."""

    return {
        "dk": records[0]["dk"],
        "perpendicular": records[0]["perpendicular"],
        "k": records[0]["k"],
        "power_per_dk": np.mean(
            [np.asarray(record["power_per_dk"], dtype=float) for record in records],
            axis=0,
        ).tolist(),
    }


def average_snapshot_records(records: dict[str, dict[str, object]]) -> dict[str, object]:
    """Average spatial products selected from one case and one time window."""

    if not records:
        return {"snapshot_count": 0}
    samples = sorted(records.values(), key=lambda sample: float(sample["time"]))
    times = np.asarray([float(sample["time"]) for sample in samples], dtype=float)
    can_integrate = bool(len(times) >= 2 and times[-1] > times[0])
    pdf_names = samples[0]["pdf"].keys()
    spectrum_names = samples[0]["spectra"].keys()
    alignment_names = set(samples[0]["alignment"].keys())
    for sample in samples[1:]:
        alignment_names.intersection_update(sample["alignment"].keys())
    transfer = [sample["pressure_transfer"] for sample in samples]
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
        "heat_flux_transport_proxy": heat_flux_ensemble,
    }


def analyze_snapshot_paths(paths: list[Path], bins: int, alignment_shells: list[int],
                           time_start: float | None = None,
                           time_end: float | None = None,
                           model_choices: dict[str, object] | None = None
                           ) -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    """Analyze and time-average selected snapshots with common PDF bin edges."""

    selected: list[Path] = []
    extrema: dict[str, list[float]] = {}
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
    ranges: dict[str, tuple[float, float]] = {}
    for name, (low, high) in extrema.items():
        if high <= low:
            delta = max(abs(low), 1.0) * 1.0e-12
            low -= delta
            high += delta
        ranges[name] = (low, high)
    records: dict[str, dict[str, object]] = {}
    for path in selected:
        fields, lengths, time = read_snapshot(path)
        records[str(path)] = analyze_fields(
            fields, lengths, time, bins, alignment_shells, ranges, model_choices
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
    x = (np.arange(shape[2]) + 0.5) * lengths[0] / shape[2]
    zz, _, xx = np.meshgrid(z, np.arange(shape[1]), x, indexing="ij")
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
    peak = int(np.argmax(velocity_power))
    derivative = periodic_gradient(fields["p_perp"] - fields["eint"], lengths)[2]
    exact = 0.1 * math.pi * np.cos(math.pi * zz)
    relative_gradient_error = float(
        np.max(np.abs(derivative - exact)) / np.max(np.abs(exact))
    )
    transfer = record["pressure_transfer"]
    pressure_work = record["pressure_work_decomposition"]
    alignment = record["alignment"].get("2", {})
    strain = np.asarray(record["pdf"]["bb_grad_velocity"]["density"])
    heat_flux = record["heat_flux_transport_proxy"]
    finite_alignment = bool(
        alignment and np.isfinite(np.asarray(alignment["density"])).all()
    )
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
    passed = bool(
        peak == 2
        and relative_gradient_error < 7.0e-3
        and abs(float(transfer["direct_real_space"])) < 1.0e-14
        and abs(float(transfer["closure_error"])) < 1.0e-14
        and finite_alignment
        and np.isfinite(strain).all()
        and heat_flux["available"]
        and heat_flux["regularized_perpendicular_power"] > 0.0
        and abs(float(heat_flux["regularized_parallel_power"])) < 1.0e-14
        and abs(float(pressure_work["anisotropic_stress_power"])) < 1.0e-14
        and correlated["anisotropic_stress_power"] < 0.0
        and abs(float(
            correlated["mks24_transfer_minus_anisotropic_stress_power"]
        )) < 1.0e-14
        and passive_work["applied_to_flow"] is False
        and correlated_quadrature["available"]
        and abs(float(correlated_quadrature["anisotropic_stress_power_integral"])
                - 2.0 * float(correlated["anisotropic_stress_power"])) < 1.0e-14
    )
    return {
        "passed": passed,
        "velocity_peak_kperp_bin": peak,
        "expected_velocity_peak_kperp_bin": 2,
        "parallel_gradient_relative_error": relative_gradient_error,
        "zero_transfer": transfer["direct_real_space"],
        "transfer_closure_error": transfer["closure_error"],
        "finite_alignment": finite_alignment,
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


def analyzed_product_curve(ensemble: dict[str, object],
                           product: str) -> tuple[np.ndarray, np.ndarray]:
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
    if family == "alignment":
        products = ensemble.get("alignment", {})
        if not isinstance(products, dict) or name not in products:
            raise ValueError(f"analyzed product is missing: {product}")
        record = products[name]
        edges = np.asarray(record["edges"], dtype=float)
        return 0.5 * (edges[1:] + edges[:-1]), np.asarray(record["density"], dtype=float)
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


def reference_curve_comparisons(result: dict[str, object],
                                manifest_path: Path) -> dict[str, object]:
    """Compare analysis products against provenance-qualified external curves."""

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
        figure_name = require_manifest_text(
            provenance, "source_figure", "digitized provenance"
        )
        figure_digest = require_manifest_text(
            provenance, "source_figure_sha256", "digitized provenance"
        )
        if not re.fullmatch(r"[0-9a-f]{64}", figure_digest):
            raise ValueError("digitized source_figure_sha256 must be lowercase SHA-256")
        figure_path = (manifest_path.parent / figure_name).resolve()
        if not figure_path.is_file():
            raise ValueError(f"digitized source figure is missing: {figure_path}")
        if figure_digest != sha256_file(figure_path):
            raise ValueError("digitized source figure checksum does not match")
        require_manifest_text(provenance, "digitization_tool", "digitized provenance")
    curves = manifest.get("curves", [])
    if not isinstance(curves, list) or not curves:
        raise ValueError("reference curves manifest must contain curves")
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
        elif isinstance(cases, dict) and case in cases:
            ensemble = cases[case].get("snapshot_ensemble", {})
        else:
            raise ValueError(f"reference curve {curve_id} selects missing case {case}")
        if not isinstance(ensemble, dict) or ensemble.get("snapshot_count", 0) == 0:
            raise ValueError(f"reference curve {curve_id} selects empty analysis")
        source_x, source_y = analyzed_product_curve(ensemble, product)
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
    return {
        "available": True,
        "definition": (
            "analysis products interpolated onto provenance-qualified reference "
            "curve coordinates"
        ),
        "manifest": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
        "provenance": provenance,
        "comparisons": comparisons,
    }


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
    command.add_argument("--time-start", type=float)
    command.add_argument("--time-end", type=float)
    command.add_argument("--reference-curves", type=Path)
    command.add_argument("--synthetic-test", action="store_true")
    args = command.parse_args()
    if (
        args.time_start is not None
        and args.time_end is not None
        and args.time_end < args.time_start
    ):
        command.error("--time-end must be greater than or equal to --time-start")

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
                snapshots, args.pdf_bins, alignment_shells, start, end, model
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
        args.time_start, args.time_end
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
    if args.reference_curves is not None:
        result["reference_curve_comparisons"] = reference_curve_comparisons(
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
