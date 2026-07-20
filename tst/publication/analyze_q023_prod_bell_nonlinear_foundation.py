#!/usr/bin/env python3
"""Source-local nonlinear Bell diagnostics with no science or launch authority."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
DECK = (
    REPO_ROOT
    / "inputs/tests/pic_q023_prod_bell_nonlinear_foundation_vl2_tsc.athinput"
)
CAMPAIGN_ID = "Q023-PROD-BELL-NONLINEAR"
CLAIM_ID = "CLAIM-PROD-BELL-NONLINEAR-NOHALL-001"
RECORD_TYPE = "q023_prod_bell_nonlinear_source_local_foundation_analysis"
K0 = 2.0 * math.pi
U_A = 1.0
B_G = 1.0
RHO0 = 1.0
EPSILON = 0.4
FIXED_CR_STREAM_PARALLEL_SPEED = U_A / EPSILON
EXPECTED_LINEAR_GROWTH = math.sqrt(1.0 - EPSILON * EPSILON)
LINEAR_WINDOW_TAU = (2.0, 6.0)
SATURATION_WINDOW_TAU = (14.0, 22.0)
NONLINEAR_ONSET_DELTA_B_OVER_BG = 0.1
SOURCE_LOCAL_FLATNESS_HEURISTIC = 0.05
MIN_WINDOW_SNAPSHOTS = 8
MAX_SNAPSHOTS = 128
SPECTRUM_DELTA_K_OVER_K0 = 0.125
LOW_K_BOUNDARY_OVER_K0 = 0.75
EXPECTED_EXTENT = (8.0 * math.sqrt(5.0), 8.0 * math.sqrt(1.25), 1.0)
EXPECTED_NX = (512, 256, 1)
EXPECTED_XMIN = (0.0, 0.0, 0.0)

_EXPECTED_DECK_VALUES = {
    ("time", "evolution"): "dynamic",
    ("time", "integrator"): "vl2",
    ("time", "cfl_number"): "0.1",
    ("time", "nlim"): "0",
    ("time", "tlim"): "4.0",
    ("mhd", "eos"): "ideal",
    ("mhd", "reconstruct"): "plm",
    ("mhd", "rsolver"): "llf",
    ("particles", "particle_type"): "cosmic_ray",
    ("particles", "ppc"): "16.0",
    ("particles", "pusher"): "boris_tsc",
    ("particles", "cr_distribution"): "center",
    ("particles", "deposit_moments"): "true",
    ("particles", "deposit_order"): "2",
    ("particles", "deposit_qscale"): "1.25e8",
    ("particles", "couple_moments_to_mhd"): "true",
    ("particles", "couple_j_to_efield_coeff"): "1.0",
    ("particles", "couple_moments_momentum_to_mhd"): "true",
    ("particles", "couple_moments_energy_to_mhd"): "true",
    ("particles", "pic_background_mode"): "coupled",
    ("particles", "pic_feedback_mode"): "coupled",
    ("particles", "pic_interp_scheme"): "tsc",
    ("particles", "pic_enable_2d3v"): "true",
    ("particles", "pic_cr_light_speed"): "2500.0",
    ("particles", "pic_cr_initial_state"): "velocity",
    ("particles", "pic_cr_hall_mode"): "off",
    ("particles", "pic_wave_damping_mode"): "off",
    ("particles", "pic_deltaf_mode"): "off",
    ("species0", "mass"): "1.0",
    ("species0", "charge"): "6.283185307179586e-6",
    ("problem", "pgen_name"): "q023_paper_bell_linear",
    ("q023_paper_bell_linear", "campaign_id"): "Q023-PAPER-BELL-LINEAR",
    ("q023_paper_bell_linear", "deck_role"):
        "source_local_runnable_preparation_not_authorized",
    ("q023_paper_bell_linear", "dimension"): "2",
    ("q023_paper_bell_linear", "epsilon"): "0.4",
    ("q023_paper_bell_linear", "rho"): "1.0",
    ("q023_paper_bell_linear", "pressure"): "1.0",
    ("q023_paper_bell_linear", "amplitude"): "1.0e-4",
    ("q023_paper_bell_linear", "b_g"): "1.0",
    ("q023_paper_bell_linear", "u_a"): "1.0",
    ("q023_paper_bell_linear", "wavelength"): "1.0",
    ("q023_paper_bell_linear", "k0"): "6.283185307179586",
    ("q023_paper_bell_linear", "c_over_v_cr"): "1000.0",
    ("q023_paper_bell_linear", "initial_eigenmode"):
        "section52_right_polarized_eigenmode",
    ("q023_prod_bell_nonlinear_foundation", "campaign_id"): CAMPAIGN_ID,
    ("q023_prod_bell_nonlinear_foundation", "claim_id"): CLAIM_ID,
    ("q023_prod_bell_nonlinear_foundation", "deck_role"):
        "source_local_foundation_not_authorized",
    ("q023_prod_bell_nonlinear_foundation", "carrier_wavelengths_x1"): "8",
    ("q023_prod_bell_nonlinear_foundation", "carrier_wavelengths_x2"): "8",
    ("q023_prod_bell_nonlinear_foundation", "linear_growth_window_tau"): "2.0,6.0",
    ("q023_prod_bell_nonlinear_foundation", "saturation_diagnostic_window_tau"):
        "14.0,22.0",
    ("q023_prod_bell_nonlinear_foundation",
     "nonlinear_onset_delta_b_rms_over_bg"): "0.1",
    ("q023_prod_bell_nonlinear_foundation", "spectrum_scope"):
        "mhd_spatial_magnetic_only",
    ("q023_prod_bell_nonlinear_foundation", "constant_current_diagnostic"):
        "riquelme_spitkovsky_2009_context_only",
    ("q023_prod_bell_nonlinear_foundation", "amplified_alfven_speed_estimator"):
        "sqrt_volume_mean_b_squared_over_volume_mean_density",
    ("q023_prod_bell_nonlinear_foundation",
     "anisotropic_pressure_saturation_relation"): "not_mapped_not_claimed",
    ("q023_prod_bell_nonlinear_foundation", "qualification_effect"):
        "none_no_launch_no_science_authority",
    ("output1", "file_type"): "bin",
    ("output1", "variable"): "mhd_w_bcc",
    ("output1", "id"): "mhd_w_bcc",
    ("output1", "dt"): "0.05",
    ("output1", "ghost_zones"): "false",
    ("output2", "file_type"): "rst",
    ("output2", "dt"): "0.5",
}


class ContractError(ValueError):
    """Raised when source-local nonlinear Bell diagnostics cannot be trusted."""


def _parse_athinput(path: Path) -> dict[str, dict[str, str]]:
    blocks: dict[str, dict[str, str]] = {}
    block: str | None = None
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), 1
    ):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            if not block or block in blocks:
                raise ContractError(f"{path}:{line_number}: invalid or duplicate block")
            blocks[block] = {}
            continue
        if block is None or "=" not in line:
            raise ContractError(f"{path}:{line_number}: malformed Athena input")
        name, value = (part.strip() for part in line.split("=", 1))
        if not name or not value or name in blocks[block]:
            raise ContractError(f"{path}:{line_number}: invalid or duplicate parameter")
        blocks[block][name] = value
    return blocks


def _require_close(label: str, measured: float, expected: float) -> None:
    if not math.isclose(measured, expected, rel_tol=1.0e-12, abs_tol=1.0e-12):
        raise ContractError(f"{label} does not match the source-local foundation")


def validate_production_deck(path: Path = DECK) -> dict[str, Any]:
    """Validate the exact source-local baseline deck without authorizing a run."""
    blocks = _parse_athinput(path)
    for (block, name), expected in _EXPECTED_DECK_VALUES.items():
        if blocks.get(block, {}).get(name) != expected:
            raise ContractError(f"{path}: <{block}>/{name} drifted")

    nx = tuple(int(blocks["mesh"][f"nx{axis}"]) for axis in (1, 2, 3))
    xmin = tuple(float(blocks["mesh"][f"x{axis}min"]) for axis in (1, 2, 3))
    xmax = tuple(float(blocks["mesh"][f"x{axis}max"]) for axis in (1, 2, 3))
    extent = tuple(hi - lo for lo, hi in zip(xmin, xmax))
    if nx != EXPECTED_NX:
        raise ContractError(f"{path}: production mesh cell counts drifted")
    for axis, (lo, expected) in enumerate(zip(xmin, EXPECTED_XMIN), 1):
        _require_close(f"{path}: x{axis}min", lo, expected)
    for axis, (value, expected) in enumerate(zip(extent, EXPECTED_EXTENT), 1):
        _require_close(f"{path}: x{axis} extent", value, expected)
    for axis in (1, 2, 3):
        if (
            blocks["mesh"][f"ix{axis}_bc"] != "periodic"
            or blocks["mesh"][f"ox{axis}_bc"] != "periodic"
        ):
            raise ContractError(f"{path}: all nonlinear Bell boundaries must be periodic")

    dx = tuple(extent[index] / nx[index] for index in range(2))
    _require_close(f"{path}: square-cell dx", dx[0], dx[1])
    parallel, _, _ = _mode_basis()
    periods = (
        K0 * parallel[0] * extent[0] / (2.0 * math.pi),
        K0 * parallel[1] * extent[1] / (2.0 * math.pi),
    )
    _require_close(f"{path}: x1 carrier periods", periods[0], 8.0)
    _require_close(f"{path}: x2 carrier periods", periods[1], 8.0)

    particles = blocks["particles"]
    stream = np.asarray(
        [float(particles[f"cr_v{axis}0"]) for axis in ("x", "y", "z")],
        dtype=float,
    )
    v_cr = float(np.linalg.norm(stream))
    ppc = float(particles["ppc"])
    qscale = float(particles["deposit_qscale"])
    charge = float(blocks["species0"]["charge"])
    light_speed = float(particles["pic_cr_light_speed"])
    _require_close(f"{path}: epsilon", U_A / v_cr, EPSILON)
    for axis, (measured, expected) in enumerate(zip(stream, v_cr * parallel), 1):
        _require_close(f"{path}: CR stream x{axis}", float(measured), float(expected))
    _require_close(
        f"{path}: j_CR",
        ppc * qscale * charge * v_cr,
        2.0 * B_G * light_speed * K0,
    )
    if float(blocks["time"]["tlim"]) * K0 * U_A < SATURATION_WINDOW_TAU[1]:
        raise ContractError(f"{path}: run ends before the fixed saturation window")

    return {
        "path": str(path.relative_to(REPO_ROOT)),
        "campaign_id": CAMPAIGN_ID,
        "analysis_scope": "source_local_foundation_only",
        "launch_authorized": False,
        "scientific_claim_authorized": False,
        "pressure_p0": 1.0,
        "epsilon_ua_over_vcr": EPSILON,
        "fixed_cr_stream_parallel_speed_over_ua":
            FIXED_CR_STREAM_PARALLEL_SPEED / U_A,
        "mesh_nx": list(nx),
        "mesh_extent": list(extent),
        "carrier_periods_by_periodic_axis": list(periods),
        "cells_per_carrier_period_by_axis": [nx[0] // 8, nx[1] // 8],
        "mhd_spatial_spectra_supported": True,
        "particle_spectra_supported": False,
    }


def _mode_basis() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    parallel = np.asarray([1.0, 2.0, 0.0], dtype=float)
    parallel /= np.linalg.norm(parallel)
    transverse_a = np.asarray([-parallel[1], parallel[0], 0.0], dtype=float)
    transverse_b = np.cross(parallel, transverse_a)
    return parallel, transverse_a, transverse_b


def _finite_time(dataset: Mapping[str, Any]) -> float:
    if "Time" not in dataset:
        raise ContractError("raw mhd_w_bcc snapshot is missing Time")
    time = float(dataset["Time"])
    if not math.isfinite(time):
        raise ContractError("raw mhd_w_bcc Time must be finite")
    return time


def _uniform_axis(
    values: Any, *, label: str, singleton_center: float | None = None
) -> tuple[np.ndarray, float, float]:
    coordinate = np.asarray(values, dtype=float)
    if (
        coordinate.ndim != 1
        or coordinate.size == 0
        or not np.all(np.isfinite(coordinate))
    ):
        raise ContractError(f"{label} must be a finite nonempty coordinate array")
    if singleton_center is not None:
        if coordinate.size != 1:
            raise ContractError(f"{label} must retain the exact singleton carrier")
        _require_close(label, float(coordinate[0]), singleton_center)
        return coordinate, 1.0, 0.0
    if coordinate.size < 8:
        raise ContractError(f"{label} requires at least eight cells")
    spacing = np.diff(coordinate)
    if np.any(spacing <= 0.0) or not np.allclose(
        spacing, spacing[0], rtol=1.0e-12, atol=1.0e-12
    ):
        raise ContractError(f"{label} must be a strictly increasing uniform grid")
    dx = float(spacing[0])
    xmin = float(coordinate[0] - 0.5 * dx)
    return coordinate, dx * coordinate.size, xmin


def _geometry(dataset: Mapping[str, Any]) -> dict[str, Any]:
    x1, l1, x1min = _uniform_axis(dataset.get("x1v"), label="x1v")
    x2, l2, x2min = _uniform_axis(dataset.get("x2v"), label="x2v")
    x3, l3, x3min = _uniform_axis(
        dataset.get("x3v"), label="x3v", singleton_center=0.5
    )
    extent = (l1, l2, l3)
    xmin = (x1min, x2min, x3min)
    for axis, (value, expected) in enumerate(zip(extent, EXPECTED_EXTENT), 1):
        _require_close(f"raw mhd_w_bcc x{axis} extent", value, expected)
    for axis, (value, expected) in enumerate(zip(xmin, EXPECTED_XMIN), 1):
        _require_close(f"raw mhd_w_bcc x{axis}min", value, expected)
    return {
        "coordinates": (x1, x2, x3),
        "nx": (x1.size, x2.size, x3.size),
        "extent": extent,
        "dx": (l1 / x1.size, l2 / x2.size, 1.0),
        "shape": (x3.size, x2.size, x1.size),
    }


def _field(
    dataset: Mapping[str, Any], name: str, geometry: Mapping[str, Any]
) -> np.ndarray:
    if name not in dataset:
        raise ContractError(f"raw mhd_w_bcc snapshot is missing {name}")
    values = np.asarray(dataset[name], dtype=float)
    if values.shape != geometry["shape"] or not np.all(np.isfinite(values)):
        raise ContractError(f"raw mhd_w_bcc {name} shape or values are invalid")
    return values


def _magnetic_spectrum(
    delta_b: np.ndarray, geometry: Mapping[str, Any]
) -> dict[str, Any]:
    field = delta_b[:, 0, :, :]
    count = field.shape[-2] * field.shape[-1]
    transformed = np.fft.fftn(field, axes=(-2, -1)) / count
    power = np.sum(np.abs(transformed) ** 2, axis=0)
    power[0, 0] = 0.0
    total_power = float(np.sum(power))
    if not math.isfinite(total_power) or total_power <= np.finfo(float).tiny:
        raise ContractError("magnetic spatial spectrum has no nonzero fluctuation power")

    nx1, nx2, _ = geometry["nx"]
    dx1, dx2, _ = geometry["dx"]
    k1 = 2.0 * math.pi * np.fft.fftfreq(nx1, d=dx1)
    k2 = 2.0 * math.pi * np.fft.fftfreq(nx2, d=dx2)
    k1_grid, k2_grid = np.meshgrid(k1, k2, indexing="xy")
    k_over_k0 = np.sqrt(k1_grid * k1_grid + k2_grid * k2_grid) / K0
    shell_index = np.rint(k_over_k0 / SPECTRUM_DELTA_K_OVER_K0).astype(int)
    shell_index[k_over_k0 > 0.0] = np.maximum(
        shell_index[k_over_k0 > 0.0], 1
    )
    shell_power = np.bincount(
        shell_index.ravel(), weights=power.ravel(), minlength=int(shell_index.max()) + 1
    )
    shell_power[0] = 0.0
    dominant_shell = int(np.argmax(shell_power))
    rows = [
        {
            "k_over_k0": float(index * SPECTRUM_DELTA_K_OVER_K0),
            "magnetic_power_fraction": float(value / total_power),
        }
        for index, value in enumerate(shell_power)
        if index > 0 and value > 0.0
    ]
    return {
        "dominant_k_over_k0": float(
            dominant_shell * SPECTRUM_DELTA_K_OVER_K0
        ),
        "dominant_wavelength_over_seed_wavelength": float(
            1.0 / (dominant_shell * SPECTRUM_DELTA_K_OVER_K0)
        ),
        "low_k_power_fraction_k_over_k0_lt_0p75": float(
            np.sum(power[k_over_k0 < LOW_K_BOUNDARY_OVER_K0]) / total_power
        ),
        "parseval_delta_b_squared": total_power,
        "shells": rows,
    }


def _snapshot_metrics(
    dataset: Mapping[str, Any], geometry: Mapping[str, Any]
) -> dict[str, Any]:
    magnetic = np.stack(
        [_field(dataset, name, geometry) for name in ("bcc1", "bcc2", "bcc3")]
    )
    velocity = np.stack(
        [_field(dataset, name, geometry) for name in ("velx", "vely", "velz")]
    )
    density = _field(dataset, "dens", geometry)
    if np.any(density <= 0.0):
        raise ContractError("raw mhd_w_bcc density must remain positive")

    mean_b = np.mean(magnetic, axis=(1, 2, 3))
    delta_b = magnetic - mean_b[:, None, None, None]
    density_sum = float(np.sum(density))
    mean_density = float(np.mean(density))
    mean_velocity = np.sum(
        density[None, ...] * velocity, axis=(1, 2, 3)
    ) / density_sum
    delta_velocity = velocity - mean_velocity[:, None, None, None]
    delta_b_squared = float(np.mean(np.sum(delta_b * delta_b, axis=0)))
    total_b_squared = float(np.mean(np.sum(magnetic * magnetic, axis=0)))
    magnetic_energy = 0.5 * delta_b_squared
    kinetic_energy = 0.5 * float(
        np.mean(density * np.sum(delta_velocity * delta_velocity, axis=0))
    )
    fluctuation_energy = magnetic_energy + kinetic_energy
    if delta_b_squared <= 0.0 or fluctuation_energy <= 0.0:
        raise ContractError("nonlinear Bell fluctuation energies must be positive")
    parallel, _, _ = _mode_basis()
    gas_velocity_parallel = float(np.dot(mean_velocity, parallel))
    relative_drift_parallel = (
        FIXED_CR_STREAM_PARALLEL_SPEED - gas_velocity_parallel
    )
    abs_relative_drift_parallel = abs(relative_drift_parallel)
    amplified_alfven_speed = math.sqrt(total_b_squared / mean_density)
    if not all(
        math.isfinite(value)
        for value in (
            mean_density,
            gas_velocity_parallel,
            relative_drift_parallel,
            amplified_alfven_speed,
        )
    ):
        raise ContractError("constant-current mechanism diagnostics must be finite")
    if mean_density <= 0.0:
        raise ContractError("mean density must remain positive")
    drift_floor = math.sqrt(np.finfo(float).tiny)
    alfven_to_drift = amplified_alfven_speed / max(
        abs_relative_drift_parallel, drift_floor
    )
    spectrum = _magnetic_spectrum(delta_b, geometry)
    return {
        "delta_b_rms_over_bg": math.sqrt(delta_b_squared) / B_G,
        "b_rms_over_bg": math.sqrt(total_b_squared) / B_G,
        "mean_b_parallel_over_bg": float(np.dot(mean_b, parallel) / B_G),
        "density_weighted_mean_gas_velocity_parallel_over_ua":
            gas_velocity_parallel / U_A,
        "fixed_cr_stream_parallel_speed_over_ua":
            FIXED_CR_STREAM_PARALLEL_SPEED / U_A,
        "configured_fixed_cr_minus_gas_relative_drift_parallel_over_ua":
            relative_drift_parallel / U_A,
        "abs_configured_fixed_cr_minus_gas_relative_drift_parallel_over_ua":
            abs_relative_drift_parallel / U_A,
        "mean_density_over_rho0": mean_density / RHO0,
        "amplified_total_field_alfven_speed_over_ua":
            amplified_alfven_speed / U_A,
        "v_a_amplified_over_abs_configured_fixed_cr_minus_gas_relative_drift":
            alfven_to_drift,
        "configured_fixed_cr_minus_gas_relative_drift_zero_within_"
        "machine_precision":
            abs_relative_drift_parallel < drift_floor,
        "mhd_fluctuation_magnetic_energy": magnetic_energy,
        "mhd_fluctuation_kinetic_energy": kinetic_energy,
        "mhd_fluctuation_magnetic_energy_fraction":
            magnetic_energy / fluctuation_energy,
        "dominant_k_over_k0": spectrum["dominant_k_over_k0"],
        "dominant_wavelength_over_seed_wavelength":
            spectrum["dominant_wavelength_over_seed_wavelength"],
        "low_k_power_fraction_k_over_k0_lt_0p75":
            spectrum["low_k_power_fraction_k_over_k0_lt_0p75"],
        "_spectrum": spectrum,
    }


def _window_fit(
    tau: np.ndarray,
    values: np.ndarray,
    window: tuple[float, float],
    *,
    label: str,
) -> tuple[np.ndarray, float, float]:
    mask = (tau >= window[0]) & (tau <= window[1])
    if np.count_nonzero(mask) < MIN_WINDOW_SNAPSHOTS:
        raise ContractError(f"{label} requires at least {MIN_WINDOW_SNAPSHOTS} snapshots")
    selected = values[mask]
    if np.any(selected <= 0.0):
        raise ContractError(f"{label} values must remain positive")
    coefficients = np.polyfit(tau[mask], np.log(selected), 1)
    fitted = np.polyval(coefficients, tau[mask])
    residual = np.log(selected) - fitted
    total = np.log(selected) - np.mean(np.log(selected))
    denominator = float(np.sum(total * total))
    r2 = 1.0 if denominator == 0.0 else 1.0 - float(
        np.sum(residual * residual)
    ) / denominator
    return mask, float(coefficients[0]), r2


def _window_summary(values: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    selected = np.asarray(values[mask], dtype=float)
    if selected.size < MIN_WINDOW_SNAPSHOTS or not np.all(np.isfinite(selected)):
        raise ContractError("fixed-window mechanism summary values are invalid")
    return {
        "start": float(selected[0]),
        "end": float(selected[-1]),
        "end_minus_start": float(selected[-1] - selected[0]),
        "median": float(np.median(selected)),
        "minimum": float(np.min(selected)),
        "maximum": float(np.max(selected)),
    }


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _raw_artifact_records(paths: Sequence[Path]) -> list[dict[str, str]]:
    records = []
    for path in paths:
        resolved = path.resolve(strict=True)
        if not resolved.is_file():
            raise ContractError(f"raw artifact is not a regular file: {resolved}")
        records.append(
            {
                "path": str(resolved),
                "sha256": _sha256_file(resolved),
            }
        )
    return records


def _validate_artifact_bindings(
    records: Sequence[Mapping[str, str]] | None,
    *,
    source_kind: str,
    snapshot_count: int,
) -> list[dict[str, str]]:
    if source_kind == "synthetic_contract_fixture":
        if records:
            raise ContractError("synthetic fixtures cannot carry raw artifact bindings")
        return []
    if records is None or len(records) != snapshot_count:
        raise ContractError("raw analysis requires one artifact binding per snapshot")
    validated = []
    for record in records:
        if set(record) != {"path", "sha256"}:
            raise ContractError("raw artifact binding keys do not match the contract")
        path = record["path"]
        sha256 = record["sha256"]
        if (
            not isinstance(path, str)
            or not path
            or not isinstance(sha256, str)
            or len(sha256) != 64
            or any(character not in "0123456789abcdef" for character in sha256)
        ):
            raise ContractError("raw artifact binding values are invalid")
        validated.append({"path": path, "sha256": sha256})
    return validated


def analyze_datasets(
    datasets: Sequence[Mapping[str, Any]],
    *,
    source_kind: str,
    raw_artifacts: Sequence[Mapping[str, str]] | None = None,
) -> dict[str, Any]:
    """Analyze decoded mhd_w_bcc snapshots under the source-local contract."""
    if source_kind not in {
        "source_local_raw_mhd_w_bcc",
        "synthetic_contract_fixture",
    }:
        raise ContractError("analysis source_kind is outside the source-local contract")
    if not datasets:
        raise ContractError("nonlinear Bell analysis requires raw snapshots")
    if len(datasets) > MAX_SNAPSHOTS:
        raise ContractError("nonlinear Bell analysis exceeds the bounded snapshot count")
    artifact_bindings = _validate_artifact_bindings(
        raw_artifacts, source_kind=source_kind, snapshot_count=len(datasets)
    )
    deck_contract = validate_production_deck()
    ordered = sorted(datasets, key=_finite_time)
    times = np.asarray([_finite_time(dataset) for dataset in ordered], dtype=float)
    if np.any(np.diff(times) <= 0.0):
        raise ContractError("raw mhd_w_bcc snapshot times must be unique and increasing")

    geometry = _geometry(ordered[0])
    metrics = []
    for dataset in ordered:
        current_geometry = _geometry(dataset)
        if (
            current_geometry["nx"] != geometry["nx"]
            or current_geometry["shape"] != geometry["shape"]
        ):
            raise ContractError("raw mhd_w_bcc geometry changed during the run")
        metrics.append(_snapshot_metrics(dataset, geometry))

    tau = times * K0 * U_A
    amplification = np.asarray(
        [metric["delta_b_rms_over_bg"] for metric in metrics], dtype=float
    )
    linear_mask, measured_growth, linear_r2 = _window_fit(
        tau, amplification, LINEAR_WINDOW_TAU, label="fixed linear-growth window"
    )
    saturation_mask, saturation_log_slope, saturation_r2 = _window_fit(
        tau,
        amplification,
        SATURATION_WINDOW_TAU,
        label="fixed saturation-diagnostic window",
    )
    saturation_values = amplification[saturation_mask]
    saturation_median = float(np.median(saturation_values))
    saturation_relative_span = float(
        (np.max(saturation_values) - np.min(saturation_values)) / saturation_median
    )
    mechanism_keys = (
        "density_weighted_mean_gas_velocity_parallel_over_ua",
        "fixed_cr_stream_parallel_speed_over_ua",
        "configured_fixed_cr_minus_gas_relative_drift_parallel_over_ua",
        "abs_configured_fixed_cr_minus_gas_relative_drift_parallel_over_ua",
        "mean_density_over_rho0",
        "amplified_total_field_alfven_speed_over_ua",
        "v_a_amplified_over_abs_configured_fixed_cr_minus_gas_relative_drift",
    )
    mechanism_summaries = {
        key: _window_summary(
            np.asarray([metric[key] for metric in metrics], dtype=float),
            saturation_mask,
        )
        for key in mechanism_keys
    }
    onset_rows = np.flatnonzero(amplification >= NONLINEAR_ONSET_DELTA_B_OVER_BG)
    onset_index = int(onset_rows[0]) if onset_rows.size else None

    selected_indices: dict[str, int] = {
        "initial": 0,
        "linear_window_end": int(np.argmin(np.abs(tau - LINEAR_WINDOW_TAU[1]))),
        "saturation_window_midpoint": int(
            np.argmin(np.abs(tau - sum(SATURATION_WINDOW_TAU) / 2.0))
        ),
        "final": len(metrics) - 1,
    }
    selected_spectra = []
    retained_indices = set()
    for label, index in selected_indices.items():
        if index in retained_indices:
            continue
        retained_indices.add(index)
        selected_spectra.append(
            {
                "label": label,
                "time": float(times[index]),
                "normalized_time_tau": float(tau[index]),
                **metrics[index]["_spectrum"],
            }
        )

    trace = []
    for time, normalized_time, metric in zip(times, tau, metrics):
        trace.append(
            {
                "time": float(time),
                "normalized_time_tau": float(normalized_time),
                **{key: value for key, value in metric.items() if key != "_spectrum"},
            }
        )

    return {
        "schema_version": 1,
        "record_type": RECORD_TYPE,
        "campaign_id": CAMPAIGN_ID,
        "claim_id": CLAIM_ID,
        "analysis_scope": "source_local_foundation_nonqualifying",
        "analysis_input_kind": source_kind,
        "qualification_effect":
            "none_source_local_diagnostics_do_not_authorize_launch_or_science_claims",
        "launch_authorized": False,
        "scientific_claim_authorized": False,
        "passed": False,
        "foundation_diagnostics_complete": True,
        "deck_contract": deck_contract,
        "raw_artifacts": artifact_bindings,
        "scientific_scope": {
            "linear_reference":
                "Sun_and_Bai_Section_5p2_no_Hall_right_polarized_mode_at_k0",
            "expected_growth_rate_over_k0_ua": EXPECTED_LINEAR_GROWTH,
            "pressure_scope":
                "p0_1_Bai_et_al_baseline_convention_not_pressure_independence_proof",
            "nonlinear_literature_scope":
                "Bell_mechanism_context_only_Q022_numeric_reference_mapping_open",
            "constant_current_mechanism_scope":
                "Riquelme_Spitkovsky_2009_context_only_configured_fixed_CR_"
                "stream_minus_density_weighted_mean_gas_velocity",
            "amplified_alfven_speed_estimator":
                "sqrt_volume_mean_total_B_squared_over_volume_mean_density",
            "constant_current_ratio_formula":
                "sqrt_volume_mean_total_B_squared_over_volume_mean_density_"
                "divided_by_abs_configured_CR_stream_minus_density_weighted_"
                "mean_gas_velocity_parallel",
            "evolving_cr_current_or_speed_measured": False,
            "anisotropic_pressure_saturation_relation":
                "not_mapped_not_claimed_normalization_and_CR_pressure_tensor_open",
            "geometry_scope":
                "deterministic_single_mode_periodic_2d_eight_wavelength_carrier",
            "spectra_scope":
                "mhd_spatial_magnetic_spectra_from_mhd_w_bcc_only",
            "particle_spectra_available": False,
            "three_dimensional_morphology_claimed": False,
        },
        "fixed_windows": {
            "linear_growth_tau": list(LINEAR_WINDOW_TAU),
            "saturation_diagnostic_tau": list(SATURATION_WINDOW_TAU),
            "window_status":
                "source_local_preregistered_engineering_windows_not_Q022_tolerances",
        },
        "linear_growth_diagnostic": {
            "snapshot_count": int(np.count_nonzero(linear_mask)),
            "measured_growth_rate_over_k0_ua": measured_growth,
            "expected_growth_rate_over_k0_ua": EXPECTED_LINEAR_GROWTH,
            "relative_difference_from_linear_reference": abs(
                measured_growth - EXPECTED_LINEAR_GROWTH
            ) / EXPECTED_LINEAR_GROWTH,
            "log_amplitude_fit_r2": linear_r2,
            "science_acceptance_threshold_frozen": False,
        },
        "nonlinear_onset_diagnostic": {
            "threshold_delta_b_rms_over_bg": NONLINEAR_ONSET_DELTA_B_OVER_BG,
            "threshold_status": "source_local_diagnostic_not_science_acceptance",
            "detected": onset_index is not None,
            "normalized_time_tau": (
                float(tau[onset_index]) if onset_index is not None else None
            ),
        },
        "saturation_window_diagnostic": {
            "snapshot_count": int(np.count_nonzero(saturation_mask)),
            "median_delta_b_rms_over_bg": saturation_median,
            "relative_span_delta_b_rms_over_bg": saturation_relative_span,
            "log_amplitude_slope_per_tau": saturation_log_slope,
            "log_amplitude_fit_r2": saturation_r2,
            "source_local_flatness_heuristic_abs_slope_limit":
                SOURCE_LOCAL_FLATNESS_HEURISTIC,
            "source_local_flatness_heuristic_met":
                abs(saturation_log_slope) <= SOURCE_LOCAL_FLATNESS_HEURISTIC,
            "heuristic_status":
                "diagnostic_only_not_reference_derived_or_science_acceptance",
            "constant_current_mechanism_diagnostic": {
                "interpretation":
                    "configured_fixed_CR_stream_reference_and_mhd_w_bcc_gas_"
                    "acceleration_only",
                "literature_scope":
                    "Riquelme_Spitkovsky_2009_constant_current_context_only",
                "science_acceptance_threshold_frozen": False,
                "anisotropic_pressure_saturation_relation":
                    "not_mapped_not_claimed_normalization_and_CR_pressure_"
                    "tensor_open",
                "fixed_window_summaries": mechanism_summaries,
            },
        },
        "snapshot_count": len(trace),
        "geometry": {
            "nx": list(geometry["nx"]),
            "extent": list(geometry["extent"]),
            "dx": list(geometry["dx"]),
        },
        "trace": trace,
        "selected_magnetic_spatial_spectra": selected_spectra,
    }


def read_binary_snapshots(
    paths: Sequence[Path],
    reader: Callable[[str], Mapping[str, Any]] | None = None,
) -> list[Mapping[str, Any]]:
    """Read source-local Athena binary snapshots without granting authority."""
    if not paths or len(paths) > MAX_SNAPSHOTS:
        raise ContractError("raw binary input count is outside the bounded contract")
    if reader is None:
        module_path = REPO_ROOT / "vis/python/bin_convert_new.py"
        spec = importlib.util.spec_from_file_location(
            "q023_prod_bell_nonlinear_bin_reader", module_path
        )
        if spec is None or spec.loader is None:
            raise ContractError("unable to load the Athena binary reader")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        reader = module.read_binary_as_athdf
    return [reader(str(path.resolve(strict=True))) for path in paths]


def analyze_binary_files(
    paths: Sequence[Path],
    reader: Callable[[str], Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Analyze raw source-local mhd_w_bcc files and bind their checksums."""
    if not paths or len(paths) > MAX_SNAPSHOTS:
        raise ContractError("raw binary input count is outside the bounded contract")
    artifacts = _raw_artifact_records(paths)
    return analyze_datasets(
        read_binary_snapshots(paths, reader=reader),
        source_kind="source_local_raw_mhd_w_bcc",
        raw_artifacts=artifacts,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", type=Path, nargs="+")
    args = parser.parse_args()
    report = analyze_binary_files(args.paths)
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
