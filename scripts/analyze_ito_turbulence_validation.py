#!/usr/bin/env python3
"""Analyze paired old/corrected Ito-2 turbulence runs."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import re
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_bin_convert():
    path = REPO_ROOT / "vis/python/bin_convert.py"
    spec = importlib.util.spec_from_file_location("athenak_bin_convert", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BIN_CONVERT = _load_bin_convert()


def _json_default(value: Any):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"cannot serialize {type(value).__name__}")


def discover_output(run_dir: Path, subdir: str, file_id: str, suffix: str) -> list[Path]:
    pattern = re.compile(rf"\.{re.escape(file_id)}\.(\d+)\.{re.escape(suffix)}$")
    matches = [path for path in (run_dir / subdir).glob("*") if pattern.search(path.name)]
    return sorted(matches, key=lambda path: int(pattern.search(path.name).group(1)))


def read_vtk_points(path: Path) -> tuple[np.ndarray, float]:
    """Read POINTS from AthenaK's big-endian legacy particle VTK output."""
    time_value = math.nan
    with path.open("rb") as stream:
        while True:
            line = stream.readline()
            if not line:
                raise ValueError(f"{path}: POINTS header not found")
            if b"AthenaK particle data at time=" in line:
                match = re.search(rb"time=\s*([+\-0-9.eE]+)", line)
                if match:
                    time_value = float(match.group(1))
            if line.startswith(b"POINTS "):
                fields = line.split()
                count = int(fields[1])
                dtype = fields[2].decode("ascii").lower()
                if dtype != "float":
                    raise ValueError(f"{path}: unsupported POINTS type {dtype}")
                points = np.fromfile(stream, dtype=">f4", count=3 * count)
                if points.size != 3 * count:
                    raise ValueError(f"{path}: truncated POINTS array")
                return points.astype(np.float64).reshape(count, 3), time_value


def read_vtk_time(path: Path) -> float:
    with path.open("rb") as stream:
        for line in stream:
            if b"AthenaK particle data at time=" in line:
                match = re.search(rb"time=\s*([+\-0-9.eE]+)", line)
                return float(match.group(1)) if match else math.nan
            if line.startswith(b"POINTS "):
                break
    return math.nan


def cic_deposit_periodic(
    points: np.ndarray,
    shape: tuple[int, int, int],
    bounds: tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
    chunk_size: int = 500_000,
) -> np.ndarray:
    """Deposit particles to cell centers with periodic 3D CIC weights."""
    nz, ny, nx = shape
    counts = np.zeros(nx * ny * nz, dtype=np.float64)
    mins = np.array([bounds[0][0], bounds[1][0], bounds[2][0]])
    widths = np.array(
        [
            bounds[0][1] - bounds[0][0],
            bounds[1][1] - bounds[1][0],
            bounds[2][1] - bounds[2][0],
        ]
    )
    dims_xyz = np.array([nx, ny, nz], dtype=np.int64)

    for start in range(0, points.shape[0], chunk_size):
        xyz = points[start: start + chunk_size]  # noqa: E203
        centered = (xyz - mins) / widths * dims_xyz - 0.5
        base = np.floor(centered).astype(np.int64)
        frac = centered - base

        for dz in (0, 1):
            iz = (base[:, 2] + dz) % nz
            wz = frac[:, 2] if dz else 1.0 - frac[:, 2]
            for dy in (0, 1):
                iy = (base[:, 1] + dy) % ny
                wy = frac[:, 1] if dy else 1.0 - frac[:, 1]
                for dx in (0, 1):
                    ix = (base[:, 0] + dx) % nx
                    wx = frac[:, 0] if dx else 1.0 - frac[:, 0]
                    flat = (iz * ny + iy) * nx + ix
                    counts += np.bincount(
                        flat, weights=wx * wy * wz, minlength=counts.size
                    )

    return counts.reshape(shape)


def normalize_density(field: np.ndarray) -> np.ndarray:
    mean = float(np.mean(field))
    if not np.isfinite(mean) or mean <= 0.0:
        raise ValueError(f"density mean must be positive and finite, got {mean}")
    return np.asarray(field, dtype=np.float64) / mean


def density_metrics(gas: np.ndarray, tracer: np.ndarray) -> dict[str, float]:
    x = gas.ravel()
    y = tracer.ravel()
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    if x.size < 2:
        raise ValueError("not enough finite cells for density metrics")

    pearson = float(np.corrcoef(x, y)[0, 1])
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept
    denom = float(np.sum((y - np.mean(y)) ** 2))
    r2_fit = (
        1.0 - float(np.sum((y - predicted) ** 2)) / denom if denom > 0.0 else math.nan
    )
    r2_one = 1.0 - float(np.sum((y - x) ** 2)) / denom if denom > 0.0 else math.nan

    positive = (x > 0.0) & (y > 0.0)
    log_ratio = np.log10(y[positive] / x[positive])
    return {
        "pearson_r": pearson,
        "pearson_r2": pearson * pearson,
        "fit_slope": float(slope),
        "fit_intercept": float(intercept),
        "r2_fit": r2_fit,
        "r2_one_to_one": r2_one,
        "normalized_l1": float(np.mean(np.abs(y - x))),
        "normalized_l2": float(np.sqrt(np.mean((y - x) ** 2))),
        "zero_tracer_fraction": float(np.mean(y <= 0.0)),
        "log10_ratio_mean": float(np.mean(log_ratio)),
        "log10_ratio_std": float(np.std(log_ratio)),
        "log10_ratio_p05": float(np.quantile(log_ratio, 0.05)),
        "log10_ratio_median": float(np.median(log_ratio)),
        "log10_ratio_p95": float(np.quantile(log_ratio, 0.95)),
    }


def column_metrics(
    gas: np.ndarray, tracer: np.ndarray
) -> tuple[dict[str, float], np.ndarray, np.ndarray]:
    gas_column = normalize_density(np.sum(gas, axis=0))
    tracer_column = normalize_density(np.sum(tracer, axis=0))
    ratio = tracer_column / gas_column
    metrics = {
        "column_pearson_r": float(
            np.corrcoef(gas_column.ravel(), tracer_column.ravel())[0, 1]
        ),
        "column_normalized_l1": float(np.mean(np.abs(tracer_column - gas_column))),
        "column_normalized_l2": float(
            np.sqrt(np.mean((tracer_column - gas_column) ** 2))
        ),
        "column_log10_ratio_std": float(np.std(np.log10(ratio))),
    }
    return metrics, gas_column, tracer_column


def isotropic_power_spectrum(field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    delta = normalize_density(field) - 1.0
    transform = np.fft.fftn(delta) / delta.size
    power = np.abs(transform) ** 2
    kz = np.fft.fftfreq(delta.shape[0]) * delta.shape[0]
    ky = np.fft.fftfreq(delta.shape[1]) * delta.shape[1]
    kx = np.fft.fftfreq(delta.shape[2]) * delta.shape[2]
    kmag = np.sqrt(
        kz[:, None, None] ** 2 + ky[None, :, None] ** 2 + kx[None, None, :] ** 2
    )
    shell = np.rint(kmag).astype(np.int64)
    max_shell = min(delta.shape) // 2
    flat_shell = shell.ravel()
    flat_power = power.ravel()
    sums = np.bincount(flat_shell, weights=flat_power, minlength=max_shell + 1)
    counts = np.bincount(flat_shell, minlength=max_shell + 1)
    valid = (np.arange(max_shell + 1) > 0) & (counts[: max_shell + 1] > 0)
    k = np.arange(max_shell + 1, dtype=np.float64)[valid]
    spectrum = sums[: max_shell + 1][valid] / counts[: max_shell + 1][valid]
    return k, spectrum


def spectrum_metrics(
    k: np.ndarray, gas_power: np.ndarray, tracer_power: np.ndarray, resolution: int
) -> dict[str, float]:
    ratio = np.divide(
        tracer_power,
        gas_power,
        out=np.full_like(tracer_power, np.nan),
        where=gas_power > 0.0,
    )

    def band_error(low: float, high: float) -> float:
        mask = (k >= low) & (k < high) & np.isfinite(ratio)
        return float(np.mean(np.abs(ratio[mask] - 1.0))) if np.any(mask) else math.nan

    return {
        "spectrum_relative_l1_large": band_error(1.0, 4.0),
        "spectrum_relative_l1_mid": band_error(4.0, max(5.0, resolution / 4.0)),
        "spectrum_relative_l1_small": band_error(
            max(5.0, resolution / 4.0), resolution / 2.0 + 1.0
        ),
    }


def jensen_shannon_divergence(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    a = a / np.sum(a)
    b = b / np.sum(b)
    middle = 0.5 * (a + b)

    def kl_divergence(p: np.ndarray, q: np.ndarray) -> float:
        mask = p > 0.0
        return float(np.sum(p[mask] * np.log2(p[mask] / q[mask])))

    return 0.5 * kl_divergence(a, middle) + 0.5 * kl_divergence(b, middle)


def _load_final_fields(
    run_dir: Path, deposition: str
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    gas_files = discover_output(run_dir, "bin", "gas", "bin")
    tracer_files = discover_output(run_dir, "bin", "tracer_ngp", "bin")
    if not gas_files or not tracer_files:
        raise FileNotFoundError(f"{run_dir}: gas or tracer binary output is missing")

    gas_data = BIN_CONVERT.read_binary_as_athdf(str(gas_files[-1]), quantities=["dens"])
    tracer_data = BIN_CONVERT.read_binary_as_athdf(
        str(tracer_files[-1]), quantities=["pdens"]
    )
    gas = np.asarray(gas_data["dens"], dtype=np.float64)
    tracer_ngp = np.asarray(tracer_data["pdens"], dtype=np.float64)
    gas_time = float(gas_data["Time"])
    tracer_time = float(tracer_data["Time"])
    if not np.isclose(gas_time, tracer_time, rtol=0.0, atol=1.0e-12):
        raise ValueError(f"{run_dir}: final gas and tracer times differ")

    selected = "ngp"
    tracer = tracer_ngp
    vtk_path = None
    if deposition in ("auto", "cic"):
        vtk_files = discover_output(run_dir, "pvtk", "particles", "part.vtk")
        candidates = []
        for path in vtk_files:
            vtk_time = read_vtk_time(path)
            if np.isfinite(vtk_time):
                candidates.append((abs(vtk_time - gas_time), vtk_time, path))
        if candidates:
            difference, vtk_time, vtk_path = min(candidates)
            if difference <= 1.0e-10:
                points, _ = read_vtk_points(vtk_path)
                bounds = (
                    (float(gas_data["x1f"][0]), float(gas_data["x1f"][-1])),
                    (float(gas_data["x2f"][0]), float(gas_data["x2f"][-1])),
                    (float(gas_data["x3f"][0]), float(gas_data["x3f"][-1])),
                )
                tracer = cic_deposit_periodic(points, gas.shape, bounds)
                selected = "cic"
            elif deposition == "cic":
                raise ValueError(
                    f"{run_dir}: nearest particle VTK time {vtk_time} "
                    f"!= gas time {gas_time}"
                )
        elif deposition == "cic":
            raise FileNotFoundError(
                f"{run_dir}: no particle VTK output for CIC deposition"
            )

    metadata = {
        "gas_file": str(gas_files[-1]),
        "tracer_ngp_file": str(tracer_files[-1]),
        "particle_vtk_file": str(vtk_path) if vtk_path else None,
        "time": gas_time,
        "cycle": int(gas_data["NumCycles"]),
        "deposition": selected,
        "particle_count": float(np.sum(tracer)),
        "shape": list(gas.shape),
    }
    return gas, tracer, metadata


def _histogram(log_ratio: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return np.histogram(log_ratio, bins=edges, density=True)[0]


def _write_metric_csv(
    path: Path, methods: dict[str, dict[str, float]], comparison: dict[str, float]
):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["method", "metric", "value"])
        for method, metrics in methods.items():
            for key, value in sorted(metrics.items()):
                writer.writerow([method, key, value])
        for key, value in sorted(comparison.items()):
            writer.writerow(["corrected_minus_old", key, value])


def _plot_columns(
    path: Path,
    gas_column: np.ndarray,
    old_column: np.ndarray,
    corrected_column: np.ndarray,
):
    log_gas = np.log10(gas_column)
    log_old = np.log10(old_column)
    log_corrected = np.log10(corrected_column)
    log_old_ratio = np.log10(old_column / gas_column)
    log_corrected_ratio = np.log10(corrected_column / gas_column)
    tracer_difference = corrected_column - old_column
    fields = [
        log_gas,
        log_old,
        log_corrected,
        log_old_ratio,
        log_corrected_ratio,
        tracer_difference,
    ]
    titles = [
        "log10 gas column",
        "log10 old tracer column",
        "log10 corrected tracer column",
        "log10(old/gas)",
        "log10(corrected/gas)",
        "corrected - old",
    ]
    density_min = min(np.min(log_gas), np.min(log_old), np.min(log_corrected))
    density_max = max(np.max(log_gas), np.max(log_old), np.max(log_corrected))
    ratio_limit = max(
        np.max(np.abs(log_old_ratio)), np.max(np.abs(log_corrected_ratio)), 1.0e-12
    )
    difference_limit = max(float(np.max(np.abs(tracer_difference))), 1.0e-12)
    figure, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for index, (axis, field, title) in enumerate(zip(axes.ravel(), fields, titles)):
        if index < 3:
            image = axis.imshow(
                field, origin="lower", cmap="viridis", vmin=density_min, vmax=density_max
            )
        else:
            limit = ratio_limit if index < 5 else difference_limit
            image = axis.imshow(
                field, origin="lower", cmap="coolwarm", vmin=-limit, vmax=limit
            )
        axis.set_title(title)
        axis.set_xticks([])
        axis.set_yticks([])
        figure.colorbar(image, ax=axis, shrink=0.82)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def analyze_pair(
    old_dir: Path, corrected_dir: Path, output_dir: Path, deposition: str
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    old_gas_raw, old_tracer_raw, old_meta = _load_final_fields(old_dir, deposition)
    new_gas_raw, new_tracer_raw, new_meta = _load_final_fields(corrected_dir, deposition)
    if old_gas_raw.shape != new_gas_raw.shape:
        raise ValueError("old and corrected gas grids have different shapes")

    gas_old = normalize_density(old_gas_raw)
    gas_new = normalize_density(new_gas_raw)
    tracer_old = normalize_density(old_tracer_raw)
    tracer_new = normalize_density(new_tracer_raw)

    gas_difference = gas_new - gas_old
    gas_reference = float(np.sqrt(np.mean(gas_old**2)))
    comparison = {
        "gas_relative_l2": float(np.sqrt(np.mean(gas_difference**2)) / gas_reference),
        "gas_max_abs": float(np.max(np.abs(gas_difference))),
        "tracer_old_corrected_l1": float(np.mean(np.abs(tracer_new - tracer_old))),
        "tracer_old_corrected_l2": float(
            np.sqrt(np.mean((tracer_new - tracer_old) ** 2))
        ),
    }

    method_metrics: dict[str, dict[str, float]] = {}
    columns: dict[str, np.ndarray] = {}
    spectra: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    gas_for_metrics = {"old": gas_old, "corrected": gas_new}
    tracer_for_metrics = {"old": tracer_old, "corrected": tracer_new}

    for method in ("old", "corrected"):
        gas = gas_for_metrics[method]
        tracer = tracer_for_metrics[method]
        metrics = density_metrics(gas, tracer)
        column_values, gas_column, tracer_column = column_metrics(gas, tracer)
        metrics.update(column_values)
        k, gas_power = isotropic_power_spectrum(gas)
        tracer_k, tracer_power = isotropic_power_spectrum(tracer)
        if not np.array_equal(k, tracer_k):
            raise ValueError("gas and tracer spectrum bins differ")
        metrics.update(spectrum_metrics(k, gas_power, tracer_power, gas.shape[-1]))
        method_metrics[method] = metrics
        columns[f"{method}_gas"] = gas_column
        columns[f"{method}_tracer"] = tracer_column
        spectra[method] = (k, gas_power, tracer_power)

    for key in method_metrics["old"]:
        comparison[f"delta_{key}"] = (
            method_metrics["corrected"][key] - method_metrics["old"][key]
        )

    edges = np.linspace(-3.0, 3.0, 121)
    centers = 0.5 * (edges[:-1] + edges[1:])
    log_ratio_old = np.log10(tracer_old[tracer_old > 0.0] / gas_old[tracer_old > 0.0])
    log_ratio_new = np.log10(tracer_new[tracer_new > 0.0] / gas_new[tracer_new > 0.0])
    pdf_old = _histogram(log_ratio_old, edges)
    pdf_new = _histogram(log_ratio_new, edges)
    comparison["ratio_pdf_js_divergence_bits"] = jensen_shannon_divergence(
        pdf_old + 1.0e-300, pdf_new + 1.0e-300
    )

    summary = {
        "old": {"snapshot": old_meta, "metrics": method_metrics["old"]},
        "corrected": {"snapshot": new_meta, "metrics": method_metrics["corrected"]},
        "comparison": comparison,
    }
    for label, directory in (("old", old_dir), ("corrected", corrected_dir)):
        metadata_path = directory / "run_metadata.json"
        if metadata_path.exists():
            summary[label]["runtime"] = json.loads(
                metadata_path.read_text(encoding="utf-8")
            )

    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    _write_metric_csv(output_dir / "metrics.csv", method_metrics, comparison)

    with (output_dir / "ratio_pdf.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["log10_tracer_over_gas", "old_pdf", "corrected_pdf"])
        writer.writerows(zip(centers, pdf_old, pdf_new))

    with (output_dir / "power_spectra.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "k",
                "old_gas_power",
                "old_tracer_power",
                "corrected_gas_power",
                "corrected_tracer_power",
            ]
        )
        old_k, old_gp, old_tp = spectra["old"]
        new_k, new_gp, new_tp = spectra["corrected"]
        if not np.array_equal(old_k, new_k):
            raise ValueError("old and corrected spectrum bins differ")
        writer.writerows(zip(old_k, old_gp, old_tp, new_gp, new_tp))

    np.savez_compressed(
        output_dir / "analysis_arrays.npz",
        gas_old=gas_old,
        gas_corrected=gas_new,
        tracer_old=tracer_old,
        tracer_corrected=tracer_new,
        gas_column=columns["old_gas"],
        old_tracer_column=columns["old_tracer"],
        corrected_tracer_column=columns["corrected_tracer"],
        k=old_k,
        old_gas_power=old_gp,
        old_tracer_power=old_tp,
        corrected_gas_power=new_gp,
        corrected_tracer_power=new_tp,
        ratio_pdf_centers=centers,
        old_ratio_pdf=pdf_old,
        corrected_ratio_pdf=pdf_new,
    )

    _plot_columns(
        output_dir / "column_density_comparison.png",
        columns["old_gas"],
        columns["old_tracer"],
        columns["corrected_tracer"],
    )

    low = min(np.quantile(log_ratio_old, 0.005), np.quantile(log_ratio_new, 0.005))
    high = max(np.quantile(log_ratio_old, 0.995), np.quantile(log_ratio_new, 0.995))
    margin = max(0.08 * (high - low), 1.0e-3)
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    axes[0].plot(centers, pdf_old, label="old Ito-2")
    axes[0].plot(centers, pdf_new, label="corrected Ito-2")
    axes[0].set_xlim(low - margin, high + margin)
    axes[0].set_xlabel(r"$\log_{10}(\rho_{\rm tracer}/\rho_{\rm gas})$")
    axes[0].set_ylabel("PDF")
    axes[0].legend()
    axes[1].plot(centers, pdf_new - pdf_old, color="black")
    axes[1].axhline(0.0, color="0.5", linewidth=0.8)
    axes[1].set_xlim(low - margin, high + margin)
    axes[1].set_xlabel(r"$\log_{10}(\rho_{\rm tracer}/\rho_{\rm gas})$")
    axes[1].set_ylabel("corrected PDF - old PDF")
    figure.savefig(output_dir / "tracer_gas_ratio_pdf.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    axes[0].loglog(old_k, old_gp, color="black", label="gas")
    axes[0].loglog(old_k, old_tp, label="old Ito-2")
    axes[0].loglog(new_k, new_tp, label="corrected Ito-2")
    axes[0].set_xlabel("k")
    axes[0].set_ylabel("P(k)")
    axes[0].legend()
    old_ratio = np.divide(
        old_tp, old_gp, out=np.full_like(old_tp, np.nan), where=old_gp > 0
    )
    new_ratio = np.divide(
        new_tp, new_gp, out=np.full_like(new_tp, np.nan), where=new_gp > 0
    )
    axes[1].semilogx(old_k, old_ratio - 1.0, label="old/gas - 1")
    axes[1].semilogx(new_k, new_ratio - 1.0, label="corrected/gas - 1")
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set_xlabel("k")
    axes[1].set_ylabel("relative power error")
    axes[1].legend()
    figure.savefig(output_dir / "density_power_spectra.png", dpi=180)
    plt.close(figure)

    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-run", type=Path, required=True)
    parser.add_argument("--corrected-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--deposition",
        choices=("auto", "cic", "ngp"),
        default="auto",
        help="Use final particle VTK CIC deposition when available, otherwise NGP.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    summary = analyze_pair(
        args.old_run, args.corrected_run, args.output_dir, args.deposition
    )
    print(json.dumps(summary["comparison"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
