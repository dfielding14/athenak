#!/usr/bin/env python3
"""Reduce the compact no-CR/Hall-off/full-Hall parallel-shock pilots."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.publication.analyze_q011_section54_outputs import (  # noqa: E402
    compose_leaf_field,
    read_athenak_binary,
)
from tst.publication.pvtk_particles import read_particle_vtk  # noqa: E402


COLORS = ("#333333", "#31688e", "#d1495b", "#35b779", "#e69f00")
FIGURE_NOTE = "engineering_proxy | not_sun_bai_reproduction=true"
LONG_MORPHOLOGY_UPSTREAM_WIDTH = 1500.0
LONG_MORPHOLOGY_DOWNSTREAM_CELLS = 12.0
PVTK_HEADER = re.compile(
    rb"# AthenaK particle data at time= ([^ ]+)  nranks= ([0-9]+)"
    rb"  cycle=([0-9]+)  variables=prtcl_all"
)
HISTORY_LABEL = re.compile(r"\[[0-9]+\]=([^\s]+)")


def _files(run: Path, directory: str, pattern: str) -> list[Path]:
    paths = sorted((run / directory).glob(pattern))
    if not paths:
        paths = sorted(run.glob(pattern))
    return paths


def _parameter(
    parameters: dict[str, dict[str, str]], block: str, name: str, default: str
) -> str:
    return parameters.get(block, {}).get(name, default)


def _real_parameter(
    parameters: dict[str, dict[str, str]], block: str, name: str, default: float
) -> float:
    return float(_parameter(parameters, block, name, str(default)))


def _bool_parameter(
    parameters: dict[str, dict[str, str]], block: str, name: str, default: bool
) -> bool:
    raw = _parameter(parameters, block, name, str(default)).strip().lower()
    if raw in {"true", "1"}:
        return True
    if raw in {"false", "0"}:
        return False
    raise RuntimeError(f"invalid boolean <{block}>/{name}={raw!r}")


def _field(dataset: Any, name: str) -> np.ndarray:
    return np.asarray(compose_leaf_field(dataset, name).values, dtype=np.float64)


def _profile(values: np.ndarray) -> np.ndarray:
    return np.mean(values, axis=tuple(range(values.ndim - 1)))


def _json_vector(values: np.ndarray) -> list[float | None]:
    return [float(value) if math.isfinite(value) else None for value in values]


def _parallel_spectrum(
    by: np.ndarray, bz: np.ndarray, mask: np.ndarray, dx: float
) -> dict[str, np.ndarray | float] | None:
    """Return transverse-field power and helicity for modes parallel to x."""
    count = int(np.count_nonzero(mask))
    if count < 8:
        return None
    by_lines = np.reshape(by[..., mask], (-1, count))
    bz_lines = np.reshape(bz[..., mask], (-1, count))
    by_lines = by_lines - np.mean(by_lines, axis=1, keepdims=True)
    bz_lines = bz_lines - np.mean(bz_lines, axis=1, keepdims=True)
    window = np.hanning(count)[None, :]
    by_hat = np.fft.rfft(by_lines * window, axis=1)[:, 1:]
    bz_hat = np.fft.rfft(bz_lines * window, axis=1)[:, 1:]
    frequency = np.fft.rfftfreq(count, d=dx)[1:]
    power = np.mean(np.abs(by_hat) ** 2 + np.abs(bz_hat) ** 2, axis=0)
    cross = np.mean(by_hat * np.conj(bz_hat), axis=0)
    helicity = np.zeros_like(power)
    positive = power > np.finfo(float).tiny
    helicity[positive] = 2.0 * np.imag(cross[positive]) / power[positive]
    if not np.any(positive):
        return None
    peak = int(np.argmax(power))
    return {
        "wavenumber": 2.0 * math.pi * frequency,
        "wavelength": 1.0 / frequency,
        "relative_power": power / power[peak],
        "helicity": np.clip(helicity, -1.0, 1.0),
        "dominant_wavelength": float(1.0 / frequency[peak]),
        "dominant_helicity": float(np.clip(helicity[peak], -1.0, 1.0)),
    }


def _shock_front(
    x: np.ndarray, density: np.ndarray, domain: tuple[float, float]
) -> float | None:
    dx = float(np.mean(np.diff(x)))
    lower = domain[0] + max(4.0 * dx, 0.01 * (domain[1] - domain[0]))
    upper = domain[1] - max(4.0 * dx, 0.05 * (domain[1] - domain[0]))
    selected = np.flatnonzero((x >= lower) & (x <= upper))
    if selected.size == 0:
        return None
    gradient = np.gradient(density, x, edge_order=2)
    index = int(selected[np.argmin(gradient[selected])])
    scale = max(float(np.max(np.abs(density))), 1.0) / max(dx, 1.0)
    if gradient[index] >= -64.0 * np.finfo(float).eps * scale:
        return None
    return float(x[index])


def _weighted_quantile(
    values: np.ndarray, weights: np.ndarray, quantiles: tuple[float, ...]
) -> list[float]:
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights)
    if cumulative.size == 0 or cumulative[-1] <= 0.0:
        return [float("nan") for _ in quantiles]
    locations = np.asarray(quantiles) * cumulative[-1]
    return [float(value) for value in np.interp(locations, cumulative, sorted_values)]


def _parse_input_header(payload: bytes, source: Path) -> dict[str, dict[str, str]]:
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise RuntimeError(f"{source}: restart input header is not ASCII") from exc
    blocks: dict[str, dict[str, str]] = {}
    active: dict[str, str] | None = None
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            name = line[1:-1]
            if name == "par_end":
                break
            active = blocks.setdefault(name, {})
            continue
        if active is not None and "=" in line:
            name, value = line.split("=", 1)
            active[name.strip()] = value.strip()
    return blocks


def _restart_parameters(path: Path) -> dict[str, dict[str, str]]:
    marker = b"<par_end>\n"
    payload = bytearray()
    with path.open("rb") as stream:
        while len(payload) <= 64 * 1024 * 1024:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            payload.extend(chunk)
            end = payload.find(marker)
            if end >= 0:
                return _parse_input_header(bytes(payload[: end + len(marker)]), path)
    raise RuntimeError(f"{path}: restart has no bounded <par_end> header")


def _restart_time(parameters: dict[str, dict[str, str]]) -> float:
    for block in parameters.values():
        if block.get("file_type") == "rst" and "last_time" in block:
            return float(block["last_time"])
    raise RuntimeError("restart header has no restart-output last_time")


def _ledger_trace(run: Path) -> list[dict[str, float | None]]:
    rows: list[dict[str, float | None]] = []
    for path in _files(run, "rst", "*.rst"):
        parameters = _restart_parameters(path)
        problem = parameters.get("problem", {})
        if "ps_injected_cr_mass_global" not in problem:
            continue
        eta = float(problem.get("ps_eta", "0"))
        injected = float(problem["ps_injected_cr_mass_global"])
        reservoir = float(problem.get("ps_mass_reservoir_global", "0"))
        enabled = problem.get("ps_enable_injection", "true").lower() in {"true", "1"}
        swept = (injected + reservoir) / eta if enabled and eta > 0.0 else None
        rows.append(
            {
                "time": _restart_time(parameters),
                "injected_cr_mass": injected,
                "mass_reservoir": reservoir,
                "swept_mass": swept,
                "removed_cr_mass": float(problem.get("ps_removed_cr_mass_global", "0")),
            }
        )
    rows.sort(key=lambda row: row["time"])
    return rows


def _history(run: Path) -> dict[str, list[float]]:
    paths = sorted(run.glob("*.mhd.hst"))
    if len(paths) != 1:
        return {}
    labels: tuple[str, ...] | None = None
    rows: list[list[float]] = []
    for raw in paths[0].read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            found = tuple(HISTORY_LABEL.findall(line))
            if found:
                if labels is not None and found != labels:
                    raise RuntimeError(f"{paths[0]}: history labels changed")
                labels = found
            continue
        if labels is not None:
            values = [float(value) for value in line.split()]
            if len(values) != len(labels):
                raise RuntimeError(f"{paths[0]}: history row width changed")
            rows.append(values)
    if labels is None or not rows:
        return {}
    data = np.asarray(rows, dtype=np.float64)
    return {name: data[:, index].tolist() for index, name in enumerate(labels)}


def _pvtk_time(path: Path) -> tuple[float, int]:
    with path.open("rb") as stream:
        header = stream.read(512)
    match = PVTK_HEADER.search(header)
    if match is None:
        raise RuntimeError(f"{path}: malformed particle execution header")
    return float(match.group(1)), int(match.group(3))


def _particle_trace(
    run: Path,
    basename: str,
    parameters: dict[str, dict[str, str]],
    shock_speed: float,
) -> tuple[list[dict[str, float | int | None]], dict[str, np.ndarray] | None]:
    light_speed = _real_parameter(
        parameters, "particles", "pic_cr_light_speed", 1.0
    )
    u0 = _real_parameter(parameters, "problem", "ps_u0", 30.0)
    pinj = _real_parameter(parameters, "problem", "ps_vinj_over_u0", math.sqrt(10.0)) * u0
    birth_min = _real_parameter(parameters, "problem", "ps_inject_t_start", 0.0)
    trace: list[dict[str, float | int | None]] = []
    latest: dict[str, np.ndarray] | None = None
    pattern = f"{basename}.prtcl_all.*.part.vtk"
    for path in _files(run, "pvtk", pattern):
        time, cycle = _pvtk_time(path)
        data = read_particle_vtk(path)
        source = np.asarray(data.scalars.get("cr_source", []), dtype=np.int64)
        birth = np.asarray(data.scalars.get("birth_time", []), dtype=np.float64)
        mask = (source == 1) & (birth >= birth_min)
        if not np.any(mask):
            trace.append(
                {
                    "time": time,
                    "cycle": cycle,
                    "particle_count": 0,
                    "p50_over_pinj": None,
                    "p90_over_pinj": None,
                    "p99_over_pinj": None,
                    "pmax_over_pinj": None,
                    "weight_fraction_above_2pinj": None,
                }
            )
            continue
        velocity = np.asarray(data.vectors["vel"], dtype=np.float64)[mask]
        weights = np.asarray(data.scalars["macro_weight"], dtype=np.float64)[mask]
        denominator = 1.0 - shock_speed * velocity[:, 0] / light_speed**2
        gamma_surface = 1.0 / math.sqrt(1.0 - (shock_speed / light_speed) ** 2)
        relative = velocity.copy()
        relative[:, 0] = (velocity[:, 0] - shock_speed) / denominator
        relative[:, 1:] = velocity[:, 1:] / (gamma_surface * denominator[:, None])
        speed2 = np.sum(relative**2, axis=1)
        gamma = 1.0 / np.sqrt(1.0 - speed2 / light_speed**2)
        ratio = gamma * np.sqrt(speed2) / pinj
        p50, p90, p99 = _weighted_quantile(ratio, weights, (0.50, 0.90, 0.99))
        total_weight = float(np.sum(weights))
        above = float(np.sum(weights[ratio >= 2.0]) / total_weight)
        trace.append(
            {
                "time": time,
                "cycle": cycle,
                "particle_count": int(ratio.size),
                "p50_over_pinj": p50,
                "p90_over_pinj": p90,
                "p99_over_pinj": p99,
                "pmax_over_pinj": float(np.max(ratio)),
                "weight_fraction_above_2pinj": above,
            }
        )
        latest = {
            "p_over_pinj": ratio,
            "weight": weights,
            "x1": np.asarray(data.points, dtype=np.float64)[mask, 0],
        }
    return trace, latest


def _moment_datasets(run: Path, basename: str, token: str) -> dict[int, Any]:
    datasets: dict[int, Any] = {}
    for path in _files(run, "bin", f"{basename}.{token}.*.bin"):
        dataset = read_athenak_binary(path)
        datasets[dataset.cycle] = dataset
    return datasets


def _analyze_case(
    label: str,
    run: Path,
    upstream_width: float,
    upstream_gap_cells: float,
    collect_long_shock: bool = False,
    long_profile_bins: int = 256,
) -> tuple[
    dict[str, Any],
    dict[str, np.ndarray] | None,
    dict[str, Any] | None,
]:
    mhd_paths = _files(run, "bin", "*.mhd_w_bcc.*.bin")
    if not mhd_paths:
        raise RuntimeError(f"{run}: no mhd_w_bcc outputs")
    basenames = {path.name.split(".mhd_w_bcc.", 1)[0] for path in mhd_paths}
    if len(basenames) != 1:
        raise RuntimeError(f"{run}: multiple mhd_w_bcc basenames: {sorted(basenames)}")
    basename = basenames.pop()
    mhd = [read_athenak_binary(path) for path in mhd_paths]
    mhd.sort(key=lambda dataset: dataset.time)
    parameters = {
        block: dict(values) for block, values in mhd[0].input_parameters.items()
    }
    requested_terminal_time = _real_parameter(
        parameters, "time", "tlim", float(mhd[-1].time)
    )
    measured_terminal_time = float(mhd[-1].time)
    gamma_gas = _real_parameter(parameters, "mhd", "gamma", 5.0 / 3.0)
    u0 = _real_parameter(parameters, "problem", "ps_u0", 30.0)
    b0 = _real_parameter(parameters, "problem", "ps_b0", 1.0)
    hall_mode = _parameter(parameters, "particles", "pic_cr_hall_mode", "off")
    if collect_long_shock and hall_mode != "full":
        raise RuntimeError(
            f"{run}: long-shock diagnostics require pic_cr_hall_mode=full"
        )
    alpha_i = _real_parameter(
        parameters, "particles", "pic_background_ion_q_over_mc", 0.0
    )
    shock_model = _parameter(
        parameters, "problem", "ps_shock_speed_model", "finite_mach"
    )
    if shock_model != "ideal_surface":
        raise RuntimeError(f"{run}: analyzer expects ps_shock_speed_model=ideal_surface")
    model_speed = 0.5 * (gamma_gas - 1.0) * u0
    domain = (mhd[0].domain_bounds[0], mhd[0].domain_bounds[1])

    moments = {
        token: _moment_datasets(run, basename, token)
        for token in ("prtcl_rho", "prtcl_jx", "prtcl_jy", "prtcl_jz")
    }
    shock: list[dict[str, float | None]] = []
    upstream: list[dict[str, float | None]] = []
    long_profiles: list[dict[str, Any]] = []
    long_spectrum: dict[str, np.ndarray | float] | None = None
    long_morphology: dict[str, Any] | None = None
    injection_start = _real_parameter(
        parameters, "problem", "ps_inject_t_start", 0.0
    )
    for dataset in mhd:
        density_grid = compose_leaf_field(dataset, "dens")
        density_field = np.asarray(density_grid.values, dtype=np.float64)
        density_profile = _profile(density_field)
        x_faces = density_grid.x1_faces
        x = 0.5 * (x_faces[:-1] + x_faces[1:])
        front = _shock_front(x, density_profile, domain)
        if front is None or dataset.time < min(5.0, 0.2 * injection_start):
            continue
        dx = float(np.mean(np.diff(x)))
        gap = upstream_gap_cells*dx
        upstream_profile = (
            (x >= front + gap) & (x <= min(front + gap + upstream_width, domain[1]))
        )
        downstream_profile = (
            (x <= front - gap) &
            (x >= max(front - gap - upstream_width, domain[0]))
        )
        upstream_density = (
            float(np.mean(density_profile[upstream_profile]))
            if np.any(upstream_profile) else None
        )
        downstream_density = (
            float(np.mean(density_profile[downstream_profile]))
            if np.any(downstream_profile) else None
        )
        compression = (
            downstream_density/upstream_density
            if upstream_density is not None and downstream_density is not None and
            upstream_density > 0.0 else None
        )
        shock.append(
            {
                "time": float(dataset.time),
                "x1": front,
                "upstream_density": upstream_density,
                "downstream_density": downstream_density,
                "compression_ratio": compression,
            }
        )
        if not all(dataset.cycle in values for values in moments.values()):
            continue
        lower = front + upstream_gap_cells * dx
        upper = min(lower + upstream_width, domain[1] - upstream_gap_cells * dx)
        mask = (x >= lower) & (x <= upper)
        if not np.any(mask):
            continue
        ux = _field(dataset, "velx")
        bx = _field(dataset, "bcc1")
        by = _field(dataset, "bcc2")
        bz = _field(dataset, "bcc3")
        qcr = _field(moments["prtcl_rho"][dataset.cycle], "prtcl_rho")
        jx = _field(moments["prtcl_jx"][dataset.cycle], "prtcl_jx")
        jy = _field(moments["prtcl_jy"][dataset.cycle], "prtcl_jy")
        jz = _field(moments["prtcl_jz"][dataset.cycle], "prtcl_jz")
        region = (..., mask)
        mean_q = float(np.mean(qcr[region]))
        mean_density = float(np.mean(density_field[region]))
        mean_j = np.asarray(
            [np.mean(jx[region]), np.mean(jy[region]), np.mean(jz[region])]
        )
        drive_x = float(np.mean((jx - qcr * ux)[region]))
        mean_bx = float(np.mean(bx[region]))
        k_bell_off = (
            abs(drive_x) / (2.0 * abs(mean_bx)) if mean_bx != 0.0 else 0.0
        )
        electron_charge = alpha_i * mean_density + mean_q
        mean_r = (
            mean_q / electron_charge if electron_charge > 0.0 else None
        )
        va_parallel = (
            abs(mean_bx) / math.sqrt(mean_density)
            if mean_density > 0.0 else 0.0
        )
        lambda_parallel = (
            abs(drive_x) / (electron_charge * va_parallel)
            if electron_charge > 0.0 and va_parallel > 0.0 else None
        )
        hall_factor = (
            1.0 + (0.5 * lambda_parallel) ** 2
            if lambda_parallel is not None else 1.0
        )
        k_bell_hall = k_bell_off / hall_factor
        k_bell = k_bell_hall if hall_mode == "full" else k_bell_off
        wavelength_off = (
            2.0 * math.pi / k_bell_off if k_bell_off > 0.0 else None
        )
        wavelength_hall = (
            2.0 * math.pi / k_bell_hall if k_bell_hall > 0.0 else None
        )
        wavelength = 2.0 * math.pi / k_bell if k_bell > 0.0 else None
        bperp = float(np.sqrt(np.mean((by**2 + bz**2)[region])) / b0)
        upstream.append(
            {
                "time": float(dataset.time),
                "x1_lower": float(lower),
                "x1_upper": float(upper),
                "mean_qcr_over_c": mean_q,
                "mean_jcr_over_c_x": float(mean_j[0]),
                "mean_jcr_over_c_y": float(mean_j[1]),
                "mean_jcr_over_c_z": float(mean_j[2]),
                "mean_drive_current_x": drive_x,
                "mean_abs_R": abs(mean_r) if mean_r is not None else None,
                "mean_parallel_Lambda": lambda_parallel,
                "predicted_hall_factor": hall_factor,
                "predicted_bell_wavenumber_hall_off": k_bell_off,
                "predicted_bell_wavelength_hall_off": wavelength_off,
                "predicted_bell_wavenumber_full_hall": k_bell_hall,
                "predicted_bell_wavelength_full_hall": wavelength_hall,
                "predicted_bell_wavenumber": k_bell,
                "predicted_bell_wavelength": wavelength,
                "bperp_rms_over_b0": bperp,
            }
        )
        if collect_long_shock and dataset.time >= injection_start:
            distance = x[mask] - front
            density_x = density_profile[mask]
            ux_x = _profile(ux)[mask]
            bx_x = _profile(bx)[mask]
            qcr_x = _profile(qcr)[mask]
            drive_x_profile = _profile(jx - qcr * ux)[mask]
            bperp_x = np.sqrt(_profile(by**2 + bz**2)[mask]) / b0
            electron_x = alpha_i * density_x + qcr_x
            va_x = np.zeros_like(density_x)
            valid_density = density_x > 0.0
            va_x[valid_density] = (
                np.abs(bx_x[valid_density]) / np.sqrt(density_x[valid_density])
            )
            lambda_x = np.full_like(density_x, np.nan)
            valid_lambda = (electron_x > 0.0) & (va_x > 0.0)
            lambda_x[valid_lambda] = (
                np.abs(drive_x_profile[valid_lambda]) /
                (electron_x[valid_lambda] * va_x[valid_lambda])
            )
            k0_x = np.zeros_like(density_x)
            valid_bx = np.abs(bx_x) > 0.0
            k0_x[valid_bx] = (
                np.abs(drive_x_profile[valid_bx]) /
                (2.0 * np.abs(bx_x[valid_bx]))
            )
            hall_factor_x = np.ones_like(density_x)
            hall_factor_x[valid_lambda] += (0.5 * lambda_x[valid_lambda]) ** 2
            kbell_x = k0_x / hall_factor_x
            wavelength_x = np.full_like(density_x, np.nan)
            positive_k = kbell_x > 0.0
            wavelength_x[positive_k] = 2.0 * math.pi / kbell_x[positive_k]
            wavelength_cells_x = wavelength_x / dx

            current_weight = np.abs(drive_x_profile) * dx
            total_current = float(np.sum(current_weight))
            precursor_width: float | None = None
            if total_current > 0.0:
                index_90 = int(
                    np.searchsorted(np.cumsum(current_weight), 0.9 * total_current)
                )
                precursor_width = float(distance[min(index_90, distance.size - 1)])
            u_relative = np.abs(ux_x - model_speed)
            gamma_nominal = k0_x * va_x / np.sqrt(hall_factor_x)
            valid_growth = np.isfinite(gamma_nominal) & (u_relative > 0.0)
            n_bell = float(
                np.sum(gamma_nominal[valid_growth] / u_relative[valid_growth]) * dx
            )

            grid = np.linspace(
                upstream_gap_cells * dx,
                upstream_gap_cells * dx + upstream_width,
                long_profile_bins,
            )
            long_profiles.append(
                {
                    "time": float(dataset.time),
                    "shock_x1": front,
                    "distance": grid,
                    "drive_current_x": np.interp(
                        grid, distance, drive_x_profile, left=np.nan, right=np.nan
                    ),
                    "parallel_Lambda": np.interp(
                        grid, distance, lambda_x, left=np.nan, right=np.nan
                    ),
                    "bperp_rms_over_b0": np.interp(
                        grid, distance, bperp_x, left=np.nan, right=np.nan
                    ),
                    "bell_wavelength_cells": np.interp(
                        grid, distance, wavelength_cells_x,
                        left=np.nan, right=np.nan,
                    ),
                    "precursor_width_90_current": precursor_width,
                    "nominal_bell_efolds": n_bell,
                    "peak_abs_drive_current_x": float(
                        np.max(np.abs(drive_x_profile))
                    ),
                    "max_parallel_Lambda": (
                        float(np.nanmax(lambda_x))
                        if np.any(np.isfinite(lambda_x)) else None
                    ),
                    "minimum_bell_wavelength_cells": (
                        float(np.nanmin(wavelength_cells_x))
                        if np.any(np.isfinite(wavelength_cells_x)) else None
                    ),
                }
            )
            spectrum_mask = np.zeros_like(mask)
            if precursor_width is not None:
                spectrum_mask = (
                    (x >= lower) & (x <= front + precursor_width)
                )
            used_fallback = int(np.count_nonzero(spectrum_mask)) < 8
            if used_fallback:
                spectrum_mask = mask
            long_spectrum = _parallel_spectrum(by, bz, spectrum_mask, dx)
            if long_spectrum is not None:
                selected_x = x[spectrum_mask]
                long_spectrum.update(
                    {
                        "window_lower_distance": float(selected_x[0] - front),
                        "window_upper_distance": float(selected_x[-1] - front),
                        "window_cells": int(selected_x.size),
                        "used_configured_window_fallback": used_fallback,
                    }
                )
            if precursor_width is not None and precursor_width > 0.0:
                upstream_extent = min(
                    LONG_MORPHOLOGY_UPSTREAM_WIDTH, precursor_width
                )
                downstream_extent = LONG_MORPHOLOGY_DOWNSTREAM_CELLS * dx
                morphology_cells = np.flatnonzero(
                    (x >= front - downstream_extent) &
                    (x <= front + upstream_extent)
                )
                if morphology_cells.size >= 2:
                    first = int(morphology_cells[0])
                    last = int(morphology_cells[-1])
                    cell_slice = slice(first, last + 1)
                    density_2d = np.mean(
                        density_field[..., cell_slice], axis=0
                    )
                    density_mean = np.mean(density_2d, axis=0, keepdims=True)
                    density_contrast = np.zeros_like(density_2d)
                    np.divide(
                        density_2d,
                        density_mean,
                        out=density_contrast,
                        where=density_mean > 0.0,
                    )
                    density_contrast -= 1.0
                    bperp_2d = (
                        np.sqrt(
                            np.mean(
                                by[..., cell_slice] ** 2 +
                                bz[..., cell_slice] ** 2,
                                axis=0,
                            )
                        ) / b0
                    )
                    drive_2d = np.mean(
                        (jx - qcr * ux)[..., cell_slice], axis=0
                    )
                    long_morphology = {
                        "time": float(dataset.time),
                        "shock_x1": front,
                        "x_faces_shock_frame": (
                            x_faces[first:last + 2] - front
                        ),
                        "y_faces": np.asarray(
                            density_grid.x2_faces, dtype=np.float64
                        ),
                        "density_contrast": density_contrast,
                        "bperp_over_b0": bperp_2d,
                        "drive_current_x": drive_2d,
                        "configured_upstream_width": (
                            LONG_MORPHOLOGY_UPSTREAM_WIDTH
                        ),
                        "precursor_width_90_current": precursor_width,
                        "upstream_extent": upstream_extent,
                        "downstream_extent": downstream_extent,
                        "x1_cells": int(last - first + 1),
                        "x2_cells": int(density_2d.shape[0]),
                    }

    measured_speed: float | None = None
    if len(shock) >= 3:
        times = np.asarray([row["time"] for row in shock])
        positions = np.asarray([row["x1"] for row in shock])
        fit = times >= max(10.0, 0.5 * injection_start)
        if np.count_nonzero(fit) >= 3:
            measured_speed = float(np.polyfit(times[fit], positions[fit], 1)[0])

    particle_trace, latest_particles = _particle_trace(
        run, basename, parameters, model_speed
    )
    history = _history(run)
    hall: dict[str, Any] = {}
    if "hall_Rmax" in history and "hall_Lmax" in history:
        hall = {
            "time": history["time"],
            "max_abs_R": history["hall_Rmax"],
            "max_Lambda": history["hall_Lmax"],
        }
    result = {
        "label": label,
        "run_directory": str(run.resolve()),
        "basename": basename,
        "model": {
            "pic_cr_hall_mode": _parameter(
                parameters, "particles", "pic_cr_hall_mode", "off"
            ),
            "injection_enabled": _bool_parameter(
                parameters, "problem", "ps_enable_injection", True
            ),
            "injection_start": injection_start,
            "eta": _real_parameter(parameters, "problem", "ps_eta", 0.0),
            "rho0": _real_parameter(parameters, "problem", "ps_rho0", 1.0),
            "u0": u0,
            "b0": b0,
            "model_shock_speed": model_speed,
            "measured_shock_speed": measured_speed,
            "artificial_particle_light_speed": _real_parameter(
                parameters, "particles", "pic_cr_light_speed", 1.0
            ),
            "background_ion_q_over_mc": _real_parameter(
                parameters, "particles", "pic_background_ion_q_over_mc", 0.0
            ),
            "requested_terminal_time": requested_terminal_time,
            "measured_terminal_time": measured_terminal_time,
        },
        "shock_front": shock,
        "upstream": upstream,
        "injection_ledger": _ledger_trace(run),
        "particle_acceleration": particle_trace,
        "hall_history": hall,
    }
    long_shock: dict[str, Any] | None = None
    if collect_long_shock:
        long_shock = {
            "label": label,
            "run_directory": str(run.resolve()),
            "dx": (
                float(long_profiles[-1]["distance"][1] -
                      long_profiles[-1]["distance"][0])
                if len(long_profiles) > 0 and
                len(long_profiles[-1]["distance"]) > 1 else None
            ),
            "profiles": long_profiles,
            "spectrum": long_spectrum,
            "morphology": long_morphology,
        }
        if long_profiles:
            # The interpolation spacing above is for plotting; retain the actual
            # mesh spacing for wavelength-in-cell and spectrum diagnostics.
            x_faces = compose_leaf_field(mhd[-1], "dens").x1_faces
            long_shock["dx"] = float(np.mean(np.diff(x_faces)))
    return result, latest_particles, long_shock


def _style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.linewidth": 0.9,
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "legend.fontsize": 8,
            "lines.linewidth": 1.8,
        }
    )


def _plot_dynamics(cases: list[dict[str, Any]], path: Path, dpi: int) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(9.0, 9.0), constrained_layout=True)
    for index, case in enumerate(cases):
        color = COLORS[index % len(COLORS)]
        label = case["label"]
        shock = case["shock_front"]
        if shock:
            time = np.asarray([row["time"] for row in shock])
            x1 = np.asarray([row["x1"] for row in shock])
            axes[0, 0].plot(time, x1, color=color, label=label)
        ledger = case["injection_ledger"]
        if ledger and case["model"]["injection_enabled"]:
            time = np.asarray([row["time"] for row in ledger])
            injected = np.asarray([row["injected_cr_mass"] for row in ledger])
            reservoir = np.asarray([row["mass_reservoir"] for row in ledger])
            axes[0, 1].plot(time, injected, color=color, label=f"{label}: injected")
            axes[0, 1].plot(
                time,
                injected + reservoir,
                color=color,
                ls="--",
                label=rf"{label}: $\eta M_{{\rm swept}}$",
            )
        upstream = case["upstream"]
        if upstream:
            time = np.asarray([row["time"] for row in upstream])
            current = np.asarray([row["mean_drive_current_x"] for row in upstream])
            wavelength = np.asarray(
                [
                    np.nan
                    if row["predicted_bell_wavelength"] is None
                    else row["predicted_bell_wavelength"]
                    for row in upstream
                ]
            )
            bperp = np.asarray([row["bperp_rms_over_b0"] for row in upstream])
            axes[1, 0].plot(time, current, color=color, label=label)
            axes[1, 1].plot(time, wavelength, color=color, label=label)
            axes[2, 0].semilogy(
                time, np.maximum(bperp, 1.0e-16), color=color, label=label
            )
        hall = case["hall_history"]
        if hall:
            time = np.asarray(hall["time"])
            axes[2, 1].plot(
                time, hall["max_abs_R"], color=color, label=rf"{label}: $\max|R|$"
            )
            axes[2, 1].plot(
                time,
                hall["max_Lambda"],
                color=color,
                ls="--",
                label=rf"{label}: $\max\Lambda$",
            )

    longest = max(cases, key=lambda item: len(item["shock_front"]))
    if longest["shock_front"]:
        time = np.asarray([row["time"] for row in longest["shock_front"]])
        speed = longest["model"]["model_shock_speed"]
        axes[0, 0].plot(
            time, speed * time, color="#777777", ls=":", label="ideal surface"
        )
    axes[0, 0].set(
        xlabel="time", ylabel=r"$x_{\rm sh}$", title="Area-averaged shock front"
    )
    axes[0, 1].set(
        xlabel="time", ylabel="mass", title="Injection and swept-mass ledger"
    )
    axes[1, 0].axhline(0.0, color="#777777", lw=0.8)
    axes[1, 0].set(
        xlabel="time",
        ylabel=r"$\langle J_{\rm cr,x}/c-Q_{\rm cr}u_x\rangle$",
        title="Upstream driving current",
    )
    axes[1, 1].set(
        xlabel="time", ylabel=r"$\lambda_{\rm Bell}=2\pi/k_{\rm Bell}$",
        title="Predicted upstream Bell scale"
    )
    axes[1, 1].set_yscale("log")
    axes[2, 0].set(
        xlabel="time", ylabel=r"$B_{\perp,\rm rms}/B_0$",
        title="Upstream magnetic growth"
    )
    axes[2, 1].set(xlabel="time", ylabel="maximum", title="CR-Hall strength")
    for panel, axis in zip("abcdef", axes.flat):
        axis.text(
            0.02, 0.96, f"({panel})", transform=axis.transAxes,
            va="top", weight="bold"
        )
        axis.grid(alpha=0.18)
        handles, labels = axis.get_legend_handles_labels()
        if handles:
            axis.legend()
    fig.text(0.995, 0.002, FIGURE_NOTE, ha="right", va="bottom", fontsize=6.5)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _plot_acceleration(
    cases: list[dict[str, Any]], latest: list[dict[str, np.ndarray] | None],
    path: Path, dpi: int
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), constrained_layout=True)
    maxima = [float(np.max(item["p_over_pinj"])) for item in latest if item is not None]
    upper = max(2.0, 1.1 * max(maxima, default=2.0))
    bins = np.geomspace(0.5, upper, 72)
    for index, (case, particles) in enumerate(zip(cases, latest)):
        color = COLORS[index % len(COLORS)]
        label = case["label"]
        if particles is not None:
            ratio = particles["p_over_pinj"]
            weights = particles["weight"]
            counts, edges = np.histogram(ratio, bins=bins, weights=weights)
            counts = counts / max(float(np.sum(weights)), np.finfo(float).tiny)
            counts = counts / np.diff(np.log(edges))
            centers = np.sqrt(edges[:-1] * edges[1:])
            positive = counts > 0.0
            axes[0].loglog(centers[positive], counts[positive], color=color, label=label)
        trace = [
            row for row in case["particle_acceleration"]
            if row["p99_over_pinj"] is not None
        ]
        if trace:
            time = np.asarray([row["time"] for row in trace])
            p90 = np.asarray([row["p90_over_pinj"] for row in trace])
            p99 = np.asarray([row["p99_over_pinj"] for row in trace])
            axes[1].semilogy(time, p90, color=color, label=f"{label}: 90th")
            axes[1].semilogy(time, p99, color=color, ls="--", label=f"{label}: 99th")
    axes[0].axvline(1.0, color="#777777", lw=0.9, ls=":")
    axes[0].set(
        xlabel=r"shock-frame $p/p_{\rm inj}$",
        ylabel=r"$(1/N)\,dN/d\ln(p/p_{\rm inj})$",
        title="Final injected-CR spectrum",
    )
    axes[1].axhline(1.0, color="#777777", lw=0.9, ls=":")
    axes[1].set(
        xlabel="time", ylabel=r"shock-frame $p/p_{\rm inj}$",
        title="High-energy particle evolution"
    )
    for panel, axis in zip("ab", axes):
        axis.text(
            0.02, 0.96, f"({panel})", transform=axis.transAxes,
            va="top", weight="bold"
        )
        axis.grid(alpha=0.18, which="both")
        handles, labels = axis.get_legend_handles_labels()
        if handles:
            axis.legend()
    fig.text(0.995, 0.002, FIGURE_NOTE, ha="right", va="bottom", fontsize=6.5)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _positive_norm(values: np.ndarray) -> LogNorm | None:
    positive = values[np.isfinite(values) & (values > 0.0)]
    if positive.size == 0:
        return None
    lower, upper = np.percentile(positive, [2.0, 98.0])
    lower = max(float(lower), np.finfo(float).tiny)
    upper = max(float(upper), 1.01 * lower)
    return LogNorm(vmin=lower, vmax=upper)


def _plot_long_morphology(
    long_shock: dict[str, Any], path: Path, dpi: int
) -> None:
    morphology = long_shock.get("morphology")
    if morphology is None:
        return
    x_faces = np.asarray(
        morphology["x_faces_shock_frame"], dtype=np.float64
    )
    y_faces = np.asarray(morphology["y_faces"], dtype=np.float64)
    density = np.asarray(morphology["density_contrast"], dtype=np.float64)
    bperp = np.asarray(morphology["bperp_over_b0"], dtype=np.float64)
    drive = np.asarray(morphology["drive_current_x"], dtype=np.float64)
    density_limit = float(np.nanpercentile(np.abs(density), 99.0))
    drive_limit = float(np.nanpercentile(np.abs(drive), 99.0))
    density_limit = max(density_limit, np.finfo(float).eps)
    drive_limit = max(drive_limit, np.finfo(float).eps)

    fig, axes = plt.subplots(
        3, 1, figsize=(10.0, 7.2), sharex=True, constrained_layout=True
    )
    images = [
        axes[0].pcolormesh(
            x_faces,
            y_faces,
            density,
            shading="flat",
            cmap="RdBu_r",
            vmin=-density_limit,
            vmax=density_limit,
            rasterized=True,
        ),
        axes[1].pcolormesh(
            x_faces,
            y_faces,
            np.ma.masked_less_equal(bperp, 0.0),
            shading="flat",
            cmap="magma",
            norm=_positive_norm(bperp),
            rasterized=True,
        ),
        axes[2].pcolormesh(
            x_faces,
            y_faces,
            drive,
            shading="flat",
            cmap="RdBu_r",
            vmin=-drive_limit,
            vmax=drive_limit,
            rasterized=True,
        ),
    ]
    labels = (
        r"$\rho/\langle\rho\rangle_y-1$",
        r"$B_\perp/B_0$",
        r"$D_x=J_{{\rm cr},x}/c-Q_{\rm cr}u_x$",
    )
    titles = (
        "Transverse density contrast",
        "Transverse magnetic amplitude",
        "CR driving current",
    )
    for panel, axis, image, label, title in zip(
        "abc", axes, images, labels, titles
    ):
        axis.axvline(0.0, color="white", lw=0.8, ls=":", alpha=0.9)
        axis.set(ylabel="transverse y", title=title)
        axis.text(
            0.01,
            0.94,
            f"({panel})",
            transform=axis.transAxes,
            va="top",
            weight="bold",
            color="white",
        )
        fig.colorbar(image, ax=axis, label=label, pad=0.015)
    axes[-1].set_xlabel(r"shock-frame distance $x-x_{\rm sh}$")
    fig.suptitle(
        f"Full-Hall near-shock morphology at t={morphology['time']:g}"
    )
    fig.text(0.995, 0.002, FIGURE_NOTE, ha="right", va="bottom", fontsize=6.5)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _plot_long_profiles(long_shock: dict[str, Any], path: Path, dpi: int) -> None:
    profiles = long_shock["profiles"]
    if not profiles:
        return
    distance = np.asarray(profiles[-1]["distance"], dtype=np.float64)
    time = np.asarray([row["time"] for row in profiles], dtype=np.float64)
    fields = (
        ("drive_current_x", r"$D_x$", "Driving current", "RdBu_r"),
        ("parallel_Lambda", r"$\Lambda_\parallel$", "CR-Hall strength", "viridis"),
        (
            "bperp_rms_over_b0",
            r"$B_{\perp,{\rm rms}}/B_0$",
            "Transverse magnetic field",
            "magma",
        ),
    )
    fig, axes = plt.subplots(3, 2, figsize=(10.0, 9.0), constrained_layout=True)
    for row_index, (name, ylabel, title, cmap) in enumerate(fields):
        matrix = np.asarray([row[name] for row in profiles], dtype=np.float64)
        terminal = matrix[-1]
        axes[row_index, 0].plot(distance, terminal, color=COLORS[row_index + 1])
        if name != "drive_current_x" and np.any(terminal > 0.0):
            axes[row_index, 0].set_yscale("log")
        axes[row_index, 0].set(
            xlabel=r"distance upstream $x-x_{\rm sh}$",
            ylabel=ylabel,
            title=f"{title} at t={time[-1]:g}",
        )
        if name == "drive_current_x":
            limit = float(np.nanpercentile(np.abs(matrix), 99.0))
            image = axes[row_index, 1].pcolormesh(
                distance,
                time,
                matrix,
                shading="auto",
                cmap=cmap,
                vmin=-limit if limit > 0.0 else None,
                vmax=limit if limit > 0.0 else None,
            )
        else:
            norm = _positive_norm(matrix)
            image = axes[row_index, 1].pcolormesh(
                distance,
                time,
                np.ma.masked_less_equal(matrix, 0.0),
                shading="auto",
                cmap=cmap,
                norm=norm,
            )
        axes[row_index, 1].set(
            xlabel=r"distance upstream $x-x_{\rm sh}$",
            ylabel="time",
            title=f"{title} evolution",
        )
        fig.colorbar(image, ax=axes[row_index, 1], label=ylabel, pad=0.02)
    for panel, axis in zip("abcdef", axes.flat):
        axis.text(
            0.02,
            0.96,
            f"({panel})",
            transform=axis.transAxes,
            va="top",
            weight="bold",
        )
        axis.grid(alpha=0.15)
    fig.set_constrained_layout_pads(h_pad=0.1)
    fig.text(0.995, 0.002, FIGURE_NOTE, ha="right", va="bottom", fontsize=6.5)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _plot_long_bell(long_shock: dict[str, Any], path: Path, dpi: int) -> None:
    profiles = long_shock["profiles"]
    if not profiles:
        return
    time = np.asarray([row["time"] for row in profiles], dtype=np.float64)
    width = np.asarray(
        [
            np.nan if row["precursor_width_90_current"] is None
            else row["precursor_width_90_current"]
            for row in profiles
        ],
        dtype=np.float64,
    )
    efolds = np.asarray(
        [row["nominal_bell_efolds"] for row in profiles], dtype=np.float64
    )
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.4), constrained_layout=True)
    axes[0, 0].plot(time, width, color=COLORS[1])
    axes[0, 0].set(
        xlabel="time",
        ylabel="distance",
        title=r"Width containing 90% of $\int |D_x|dx$",
    )
    axes[0, 1].plot(time, efolds, color=COLORS[2])
    axes[0, 1].set(
        xlabel="time", ylabel=r"nominal $N_{\rm Bell}$",
        title="Local-growth opportunity"
    )

    spectrum = long_shock["spectrum"]
    if spectrum is not None:
        dx = float(long_shock["dx"])
        wavelength_cells = np.asarray(spectrum["wavelength"]) / dx
        order = np.argsort(wavelength_cells)
        power = np.asarray(spectrum["relative_power"])[order]
        helicity = np.asarray(spectrum["helicity"])[order]
        wavelength_cells = wavelength_cells[order]
        dominant = float(spectrum["dominant_wavelength"]) / dx
        axes[1, 0].loglog(wavelength_cells, power, color=COLORS[3])
        axes[1, 0].axvline(dominant, color="#777777", ls=":", lw=1.0)
        axes[1, 0].set(
            xlabel=r"parallel wavelength / $\Delta x$",
            ylabel="relative transverse power",
            title=rf"Final upstream spectrum; peak = {dominant:.1f} cells",
        )
        axes[1, 1].semilogx(wavelength_cells, helicity, color=COLORS[4])
        axes[1, 1].axhline(0.0, color="#777777", lw=0.8)
        axes[1, 1].axvline(dominant, color="#777777", ls=":", lw=1.0)
        axes[1, 1].set(
            xlabel=r"parallel wavelength / $\Delta x$",
            ylabel=r"$2\,\mathrm{Im}(\hat B_y\hat B_z^*)/P_\perp$",
            title="Circular-polarization diagnostic",
            ylim=(-1.05, 1.05),
        )
    else:
        for axis in axes[1]:
            axis.text(0.5, 0.5, "insufficient spectral signal", ha="center")
            axis.set_axis_off()
    for panel, axis in zip("abcd", axes.flat):
        if axis.axison:
            axis.text(
                0.02,
                0.96,
                f"({panel})",
                transform=axis.transAxes,
                va="top",
                weight="bold",
            )
            axis.grid(alpha=0.18, which="both")
    fig.set_constrained_layout_pads(h_pad=0.1)
    fig.text(0.995, 0.002, FIGURE_NOTE, ha="right", va="bottom", fontsize=6.5)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _write_long_summary(
    long_shock: dict[str, Any], figures: list[Path], output: Path
) -> Path:
    profiles = long_shock["profiles"]
    terminal: dict[str, Any] | None = None
    trace: list[dict[str, float | None]] = []
    for row in profiles:
        trace.append(
            {
                "time": row["time"],
                "shock_x1": row["shock_x1"],
                "precursor_width_90_current": row[
                    "precursor_width_90_current"
                ],
                "nominal_bell_efolds": row["nominal_bell_efolds"],
                "peak_abs_drive_current_x": row["peak_abs_drive_current_x"],
                "max_parallel_Lambda": row["max_parallel_Lambda"],
                "minimum_bell_wavelength_cells": row[
                    "minimum_bell_wavelength_cells"
                ],
            }
        )
    if profiles:
        row = profiles[-1]
        terminal = {
            "time": row["time"],
            "shock_x1": row["shock_x1"],
            "distance_upstream": _json_vector(row["distance"]),
            "drive_current_x": _json_vector(row["drive_current_x"]),
            "parallel_Lambda": _json_vector(row["parallel_Lambda"]),
            "bperp_rms_over_b0": _json_vector(row["bperp_rms_over_b0"]),
            "bell_wavelength_cells": _json_vector(
                row["bell_wavelength_cells"]
            ),
        }
    spectrum_json: dict[str, Any] | None = None
    spectrum = long_shock["spectrum"]
    if spectrum is not None:
        dx = float(long_shock["dx"])
        spectrum_json = {
            "wavenumber": _json_vector(np.asarray(spectrum["wavenumber"])),
            "wavelength_cells": _json_vector(
                np.asarray(spectrum["wavelength"]) / dx
            ),
            "relative_power": _json_vector(
                np.asarray(spectrum["relative_power"])
            ),
            "helicity": _json_vector(np.asarray(spectrum["helicity"])),
            "dominant_wavelength": spectrum["dominant_wavelength"],
            "dominant_wavelength_cells": spectrum["dominant_wavelength"] / dx,
            "dominant_helicity": spectrum["dominant_helicity"],
            "window_lower_distance": spectrum.get("window_lower_distance"),
            "window_upper_distance": spectrum.get("window_upper_distance"),
            "window_cells": spectrum.get("window_cells"),
            "used_configured_window_fallback": spectrum.get(
                "used_configured_window_fallback"
            ),
        }
    morphology_json: dict[str, Any] | None = None
    morphology = long_shock.get("morphology")
    if morphology is not None:
        x_faces = np.asarray(
            morphology["x_faces_shock_frame"], dtype=np.float64
        )
        y_faces = np.asarray(morphology["y_faces"], dtype=np.float64)
        morphology_json = {
            "time": morphology["time"],
            "shock_x1": morphology["shock_x1"],
            "shock_frame_x_lower": float(x_faces[0]),
            "shock_frame_x_upper": float(x_faces[-1]),
            "transverse_y_lower": float(y_faces[0]),
            "transverse_y_upper": float(y_faces[-1]),
            "configured_upstream_width": morphology[
                "configured_upstream_width"
            ],
            "precursor_width_90_current": morphology[
                "precursor_width_90_current"
            ],
            "upstream_extent": morphology["upstream_extent"],
            "downstream_extent": morphology["downstream_extent"],
            "x1_cells": morphology["x1_cells"],
            "x2_cells": morphology["x2_cells"],
            "density_contrast": "rho/<rho>_y-1 at each x",
        }
    summary = {
        "schema_version": 1,
        "analysis": "full_hall_long_parallel_shock",
        "evidence_class": "engineering_proxy",
        "not_sun_bai_reproduction": True,
        "case_label": long_shock["label"],
        "run_directory": long_shock["run_directory"],
        "assumptions": {
            "precursor_width": (
                "distance from the shock to the cell enclosing 90% of the "
                "configured upstream-window integral of abs(D_x)"
            ),
            "local_bell_scale": (
                "k0=abs(D_x)/(2*abs(Bx)); k_Bell=k0/[1+(Lambda/2)^2]; "
                "wavelength=2*pi/k_Bell"
            ),
            "nominal_bell_efolds": (
                "high-Mach, local-uniform approximation integral["
                "k0*v_A/sqrt(1+(Lambda/2)^2)/"
                "abs(u_x-v_shock)] dx over the configured upstream window; "
                "this applies Bai et al. (2015) equations 29-31 locally and is "
                "not an inhomogeneous linear-mode solution"
            ),
            "parallel_spectrum": (
                "final-snapshot Hann FFT along x from the shock gap through the "
                "90%-current precursor width; the configured upstream window is "
                "used only when that region has fewer than 8 cells; transverse-"
                "field power is averaged over transverse lines after removing "
                "each line mean"
            ),
            "polarization": (
                "helicity=2*Im(By_hat*conj(Bz_hat))/(abs(By_hat)^2+"
                "abs(Bz_hat)^2) using the numpy exp(-ikx) convention"
            ),
            "morphology_window": (
                "final analyzed snapshot from 12 cells downstream of the "
                "measured shock through min(1500 code units, the 90%-current "
                "precursor width) upstream; x3 is averaged"
            ),
        },
        "figures": [str(path) for path in figures if path.is_file()],
        "time_trace": trace,
        "terminal_profile": terminal,
        "final_upstream_parallel_spectrum": spectrum_json,
        "morphology_window": morphology_json,
    }
    path = output / "full_hall_long_shock_analysis.json"
    path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return path


def _parse_case(raw: str) -> tuple[str, Path]:
    if "=" not in raw:
        raise argparse.ArgumentTypeError("case must be LABEL=RUN_DIRECTORY")
    label, directory = raw.split("=", 1)
    if not label or not directory:
        raise argparse.ArgumentTypeError("case must be LABEL=RUN_DIRECTORY")
    path = Path(directory).expanduser().resolve()
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"case directory does not exist: {path}")
    return label, path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Analyze compact controlled CR-Hall parallel shocks."
    )
    parser.add_argument(
        "--case", action="append", type=_parse_case, required=True,
        help="repeat LABEL=RUN_DIRECTORY for no-CR, Hall-off, and full-Hall cases"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--upstream-width", type=float, default=300.0)
    parser.add_argument("--upstream-gap-cells", type=float, default=6.0)
    parser.add_argument("--dpi", type=int, default=320)
    parser.add_argument("--require-complete", action="store_true")
    parser.add_argument(
        "--long-full-hall-case",
        metavar="LABEL",
        help=(
            "add long-shock precursor, Bell-opportunity, and spectral diagnostics "
            "for this one full-Hall case"
        ),
    )
    parser.add_argument("--long-profile-bins", type=int, default=256)
    args = parser.parse_args()
    if args.upstream_width <= 0.0 or args.upstream_gap_cells < 0.0:
        parser.error("upstream width must be positive and gap must be non-negative")
    if args.long_profile_bins < 8:
        parser.error("--long-profile-bins must be at least 8")
    labels = [label for label, _ in args.case]
    if len(labels) != len(set(labels)):
        parser.error("--case labels must be unique")
    if (
        args.long_full_hall_case is not None and
        args.long_full_hall_case not in labels
    ):
        parser.error("--long-full-hall-case must match one --case label")

    output = args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    cases: list[dict[str, Any]] = []
    latest: list[dict[str, np.ndarray] | None] = []
    long_shock: dict[str, Any] | None = None
    for label, run in args.case:
        result, particles, long_result = _analyze_case(
            label,
            run,
            args.upstream_width,
            args.upstream_gap_cells,
            collect_long_shock=(label == args.long_full_hall_case),
            long_profile_bins=args.long_profile_bins,
        )
        cases.append(result)
        latest.append(particles)
        if long_result is not None:
            long_shock = long_result

    _style()
    dynamics = output / "controlled_shock_dynamics.png"
    acceleration = output / "controlled_shock_cr_acceleration.png"
    _plot_dynamics(cases, dynamics, args.dpi)
    _plot_acceleration(cases, latest, acceleration, args.dpi)
    long_figures: list[Path] = []
    long_summary: Path | None = None
    if long_shock is not None:
        long_profiles = output / "full_hall_long_shock_profiles.png"
        long_bell = output / "full_hall_long_shock_bell_diagnostics.png"
        long_morphology = output / "full_hall_long_shock_morphology.png"
        _plot_long_profiles(long_shock, long_profiles, args.dpi)
        _plot_long_bell(long_shock, long_bell, args.dpi)
        _plot_long_morphology(long_shock, long_morphology, args.dpi)
        long_figures = [long_profiles, long_bell, long_morphology]
        long_summary = _write_long_summary(long_shock, long_figures, output)
    completeness: dict[str, Any] = {"cases": {}}
    for case in cases:
        injected = case["model"]["injection_enabled"]
        injection_start = case["model"]["injection_start"]
        checks = {
            "terminal_output": (
                case["model"]["measured_terminal_time"] >=
                0.99 * case["model"]["requested_terminal_time"]
            ),
            "shock_track": len(case["shock_front"]) >= 3,
            "shock_speed_fit": case["model"]["measured_shock_speed"] is not None,
        }
        if injected:
            post_injection = [
                row for row in case["upstream"]
                if row["time"] >= injection_start
            ]
            checks.update(
                {
                    "injected_mass_ledger": any(
                        row["injected_cr_mass"] > 0.0
                        for row in case["injection_ledger"]
                    ),
                    "post_injection_current": any(
                        math.isfinite(row["mean_drive_current_x"]) and
                        abs(row["mean_drive_current_x"]) > 1.0e-12
                        for row in post_injection
                    ),
                    "post_injection_bell_scale": any(
                        row["predicted_bell_wavelength"] is not None and
                        math.isfinite(row["predicted_bell_wavelength"]) and
                        row["predicted_bell_wavelength"] > 0.0
                        for row in post_injection
                    ),
                    "post_injection_magnetic_measurement": any(
                        math.isfinite(row["bperp_rms_over_b0"]) and
                        row["bperp_rms_over_b0"] > 0.0
                        for row in post_injection
                    ),
                    "particle_acceleration": any(
                        row["particle_count"] > 0
                        for row in case["particle_acceleration"]
                    ),
                }
            )
        if case["model"]["pic_cr_hall_mode"] == "full":
            hall = case["hall_history"]
            hall_rows = zip(
                hall.get("time", []), hall.get("max_abs_R", []),
                hall.get("max_Lambda", [])
            )
            checks["post_injection_hall_history"] = any(
                time >= injection_start and math.isfinite(rmax) and
                math.isfinite(lambda_max) and lambda_max > 0.0
                for time, rmax, lambda_max in hall_rows
            )
        completeness["cases"][case["label"]] = checks
    completeness["passed"] = all(
        value
        for checks in completeness["cases"].values()
        for value in checks.values()
    )
    summary = {
        "schema_version": 1,
        "analysis": "compact_controlled_cr_hall_parallel_shock",
        "evidence_class": "engineering_proxy",
        "not_sun_bai_reproduction": True,
        "assumptions": {
            "shock_front": "strongest negative gradient of transverse-averaged density",
            "upstream_window": (
                f"measured front + {args.upstream_gap_cells} cells through "
                f"{args.upstream_width:g} code-length units farther upstream"
            ),
            "bell_scale": (
                "k0=abs(<Jcr_x/c-Qcr*u_x>)/(2*abs(<Bx>)); full-Hall cases "
                "use k_Bell=k0/[1+(Lambda_parallel/2)^2], while Hall-off uses k0; "
                "lambda_Bell=2*pi/k_Bell"
            ),
            "swept_mass": "(injected CR mass + fractional reservoir)/eta",
            "particle_momentum": "inverse artificial-C boost to the modeled shock frame",
        },
        "figures": [
            str(path) for path in [dynamics, acceleration, *long_figures]
            if path.is_file()
        ],
        "completeness": completeness,
        "cases": cases,
    }
    summary_path = output / "controlled_shock_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "summary": str(summary_path),
        "long_shock_summary": str(long_summary) if long_summary else None,
        "figures": summary["figures"],
        "completeness": completeness,
    }, indent=2))
    return 0 if (not args.require_complete or completeness["passed"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
