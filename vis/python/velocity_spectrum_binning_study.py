#!/usr/bin/env python3
"""Compare velocity-spectrum binning choices for a single 2D snapshot."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from scalar_mixing_stream_run_analysis import (
    _canonical_field,
    _read_snapshot,
    _spectrum_kmax_from_field,
    _spectrum_shells_2d,
    _velocity_spectrum_2d,
)


LOW_K_MIN = 8.0
LOW_K_MAX = 40.0


@dataclass(frozen=True)
class BinningOption:
    key: str
    label: str
    family: str
    value: int


def _log_rebin(
    kvals: np.ndarray,
    spectrum: np.ndarray,
    counts: np.ndarray,
    *,
    kmax: int,
    nbins: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    edge_values = np.unique(np.rint(np.geomspace(1.0, kmax + 1.0, nbins + 1)).astype(int))
    if edge_values[0] != 1:
        edge_values = np.insert(edge_values, 0, 1)
    if edge_values[-1] != kmax + 1:
        edge_values = np.append(edge_values, kmax + 1)

    centers = []
    rebinned = []
    mode_counts = []
    widths = []
    for lo, hi in zip(edge_values[:-1], edge_values[1:]):
        mask = (kvals >= lo) & (kvals < hi)
        if not np.any(mask):
            continue
        centers.append(np.sqrt(lo * (hi - 1)))
        rebinned.append(np.sum(spectrum[mask]))
        mode_counts.append(int(np.sum(counts[mask])))
        widths.append(hi - lo)
    return (
        np.asarray(centers, dtype=np.float64),
        np.asarray(rebinned, dtype=np.float64),
        np.asarray(mode_counts, dtype=np.int64),
        np.asarray(widths, dtype=np.int64),
    )


def _linear_rebin(
    kvals: np.ndarray,
    spectrum: np.ndarray,
    counts: np.ndarray,
    *,
    kmax: int,
    width: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    centers = []
    rebinned = []
    mode_counts = []
    widths = []
    for lo in range(1, kmax + 1, width):
        hi = min(lo + width, kmax + 1)
        mask = (kvals >= lo) & (kvals < hi)
        if not np.any(mask):
            continue
        centers.append(np.sqrt(lo * (hi - 1)))
        rebinned.append(np.sum(spectrum[mask]))
        mode_counts.append(int(np.sum(counts[mask])))
        widths.append(hi - lo)
    return (
        np.asarray(centers, dtype=np.float64),
        np.asarray(rebinned, dtype=np.float64),
        np.asarray(mode_counts, dtype=np.int64),
        np.asarray(widths, dtype=np.int64),
    )


def _adaptive_rebin(
    kvals: np.ndarray,
    spectrum: np.ndarray,
    counts: np.ndarray,
    *,
    kmax: int,
    min_modes: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    centers = []
    rebinned = []
    mode_counts = []
    widths = []
    lo = 1
    while lo <= kmax:
        hi = lo + 1
        while hi <= kmax + 1:
            mask = (kvals >= lo) & (kvals < hi)
            if np.sum(counts[mask]) >= min_modes or hi == kmax + 1:
                break
            hi += 1
        mask = (kvals >= lo) & (kvals < hi)
        if not np.any(mask):
            lo = hi
            continue
        centers.append(np.sqrt(lo * (hi - 1)))
        rebinned.append(np.sum(spectrum[mask]))
        mode_counts.append(int(np.sum(counts[mask])))
        widths.append(hi - lo)
        lo = hi
    return (
        np.asarray(centers, dtype=np.float64),
        np.asarray(rebinned, dtype=np.float64),
        np.asarray(mode_counts, dtype=np.int64),
        np.asarray(widths, dtype=np.int64),
    )


def _local_curvature_metric(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 3:
        return float("inf")
    residuals = []
    for idx in range(1, x.size - 1):
        x0, x1, x2 = x[idx - 1], x[idx], x[idx + 1]
        y0, y1, y2 = y[idx - 1], y[idx], y[idx + 1]
        y_linear = y0 + (y2 - y0) * (x1 - x0) / (x2 - x0)
        residuals.append(y1 - y_linear)
    return float(np.sqrt(np.mean(np.square(residuals))))


def _evaluate_option(
    option: BinningOption,
    kvals: np.ndarray,
    spectrum: np.ndarray,
    counts: np.ndarray,
    *,
    kmax: int,
) -> dict[str, object]:
    if option.family == "log":
        centers, rebinned, mode_counts, widths = _log_rebin(
            kvals,
            spectrum,
            counts,
            kmax=kmax,
            nbins=option.value,
        )
    elif option.family == "linear":
        centers, rebinned, mode_counts, widths = _linear_rebin(
            kvals,
            spectrum,
            counts,
            kmax=kmax,
            width=option.value,
        )
    elif option.family == "adaptive":
        centers, rebinned, mode_counts, widths = _adaptive_rebin(
            kvals,
            spectrum,
            counts,
            kmax=kmax,
            min_modes=option.value,
        )
    else:
        raise ValueError(f"Unsupported family {option.family}")

    density = rebinned / widths
    low_k_mask = (centers >= LOW_K_MIN) & (centers <= LOW_K_MAX)
    x_low = np.log10(centers[low_k_mask])
    y_low = np.log10(density[low_k_mask])
    roughness = _local_curvature_metric(x_low, y_low)
    bins_lowk = int(np.sum(low_k_mask))
    mean_width_lowk = float(np.mean(widths[low_k_mask])) if bins_lowk else float("inf")
    median_modes_lowk = float(np.median(mode_counts[low_k_mask])) if bins_lowk else 0.0
    min_modes_lowk = int(np.min(mode_counts[low_k_mask])) if bins_lowk else 0

    return {
        "key": option.key,
        "label": option.label,
        "family": option.family,
        "value": option.value,
        "centers": centers,
        "rebinned": rebinned,
        "density": density,
        "mode_counts": mode_counts,
        "widths": widths,
        "roughness_lowk": roughness,
        "bins_lowk": bins_lowk,
        "mean_width_lowk": mean_width_lowk,
        "median_modes_lowk": median_modes_lowk,
        "min_modes_lowk": min_modes_lowk,
    }


def _pick_best(results: list[dict[str, object]]) -> dict[str, object]:
    candidates = [item for item in results if int(item["bins_lowk"]) >= 5]
    if not candidates:
        candidates = results
    return min(
        candidates,
        key=lambda item: (
            float(item["roughness_lowk"]),
            float(item["mean_width_lowk"]),
            -float(item["median_modes_lowk"]),
        ),
    )


def _save_summary_csv(path: Path, results: list[dict[str, object]], best_key: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "key",
                "label",
                "family",
                "value",
                "bins_lowk",
                "mean_width_lowk",
                "median_modes_lowk",
                "min_modes_lowk",
                "roughness_lowk",
                "recommended",
            ],
        )
        writer.writeheader()
        for item in results:
            writer.writerow(
                {
                    "key": item["key"],
                    "label": item["label"],
                    "family": item["family"],
                    "value": item["value"],
                    "bins_lowk": int(item["bins_lowk"]),
                    "mean_width_lowk": f"{float(item['mean_width_lowk']):.3f}",
                    "median_modes_lowk": f"{float(item['median_modes_lowk']):.1f}",
                    "min_modes_lowk": int(item["min_modes_lowk"]),
                    "roughness_lowk": f"{float(item['roughness_lowk']):.6f}",
                    "recommended": "yes" if item["key"] == best_key else "no",
                }
            )


def _make_plot(
    output: Path,
    raw_k: np.ndarray,
    raw_e: np.ndarray,
    results: list[dict[str, object]],
    best_key: str,
    *,
    kmax: int,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    ax_full, ax_zoom, ax_modes, ax_score = axes.ravel()

    ax_full.loglog(raw_k, raw_e, color="0.75", linewidth=0.9, alpha=0.9, label="raw shell spectrum")
    ax_zoom.loglog(raw_k, raw_e, color="0.85", linewidth=0.8, alpha=0.9, label="raw shell spectrum")

    colors = plt.get_cmap("tab10")(np.linspace(0.0, 0.9, len(results)))
    full_density_values = []
    zoom_density_values = []
    for color, item in zip(colors, results):
        centers = np.asarray(item["centers"])
        density = np.asarray(item["density"])
        modes = np.asarray(item["mode_counts"])
        linewidth = 3.0 if item["key"] == best_key else 1.7
        alpha = 1.0 if item["key"] == best_key else 0.9
        label = item["label"] + (" (recommended)" if item["key"] == best_key else "")
        ax_full.loglog(centers, density, color=color, linewidth=linewidth, alpha=alpha, label=label)
        ax_zoom.loglog(centers, density, color=color, linewidth=linewidth, alpha=alpha)
        low_k_mask = (centers >= LOW_K_MIN) & (centers <= LOW_K_MAX)
        full_mask = (centers >= 2.0) & (centers <= min(float(kmax), 512.0))
        full_density_values.extend(density[full_mask].tolist())
        zoom_density_values.extend(density[low_k_mask].tolist())
        ax_modes.semilogx(
            centers[low_k_mask],
            modes[low_k_mask],
            color=color,
            linewidth=linewidth,
            alpha=alpha,
            marker="o",
            markersize=3.5,
        )
        ax_score.scatter(
            float(item["median_modes_lowk"]),
            float(item["roughness_lowk"]),
            color=color,
            s=90 if item["key"] == best_key else 55,
            alpha=alpha,
        )
        ax_score.annotate(
            item["label"],
            (float(item["median_modes_lowk"]), float(item["roughness_lowk"])),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
        )

    if full_density_values:
        ylo = min(full_density_values) / 2.0
        yhi = max(full_density_values) * 1.5
        ax_full.set_ylim(ylo, yhi)
    if zoom_density_values:
        ylo = min(zoom_density_values) / 1.4
        yhi = max(zoom_density_values) * 1.2
        ax_zoom.set_ylim(ylo, yhi)

    ax_full.set_title("Width-Normalized Spectrum Comparison")
    ax_full.set_xlabel("k")
    ax_full.set_ylabel(r"$E_v(k) / \Delta k$")
    ax_full.set_xlim(1.0, float(kmax))
    ax_full.grid(True, which="both", alpha=0.22)
    ax_full.legend(loc="upper right", fontsize=8, frameon=False, ncol=2)

    ax_zoom.set_title(f"Low-k Zoom ({LOW_K_MIN:g} <= k <= {LOW_K_MAX:g})")
    ax_zoom.set_xlabel("k")
    ax_zoom.set_ylabel(r"$E_v(k) / \Delta k$")
    ax_zoom.set_xlim(LOW_K_MIN * 0.7, LOW_K_MAX * 1.2)
    ax_zoom.grid(True, which="both", alpha=0.22)

    ax_modes.set_title("Modes per Displayed Bin in the Low-k Range")
    ax_modes.set_xlabel("bin center k")
    ax_modes.set_ylabel("Fourier modes per bin")
    ax_modes.grid(True, which="both", alpha=0.22)

    ax_score.set_title("Smoothness vs Averaging Strength")
    ax_score.set_xlabel(f"median modes/bin for {LOW_K_MIN:g} <= k <= {LOW_K_MAX:g}")
    ax_score.set_ylabel("low-k roughness metric")
    ax_score.grid(True, alpha=0.22)

    fig.suptitle("Velocity Spectrum Binning Study", fontsize=16)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-file",
        type=Path,
        default=Path(
            "/lustre/orion/ast207/proj-shared/dfielding/autonomous_scalars/no_RHS/data/"
            "scalar_mix_stream_16384_step_kappa1e-7_noRHS/bin/"
            "scalar_mix_stream_16384_step_kappa1e-7_noRHS.hydro_w.00000.bin"
        ),
        help="Velocity snapshot to analyze.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "/lustre/orion/ast207/proj-shared/dfielding/autonomous_scalars/no_RHS/data/"
            "scalar_mix_stream_16384_step_kappa1e-7_noRHS/analysis/velocity_binning_study"
        ),
        help="Directory for comparison outputs.",
    )
    args = parser.parse_args()

    snapshot = _read_snapshot(args.input_file, ["velx", "vely", "velz"])
    vx = _canonical_field(np.asarray(snapshot["velx"]))
    vy = _canonical_field(np.asarray(snapshot["vely"]))
    vz = _canonical_field(np.asarray(snapshot["velz"]))

    kvals, spectrum = _velocity_spectrum_2d(vx, vy, vz)
    shell = _spectrum_shells_2d(vx.shape[1], vx.shape[0])
    counts = np.bincount(shell.ravel())[1 : len(spectrum) + 1]
    kmax = _spectrum_kmax_from_field(vx)

    options = [
        BinningOption("log44", "log nbins=44", "log", 44),
        BinningOption("log32", "log nbins=32", "log", 32),
        BinningOption("log24", "log nbins=24", "log", 24),
        BinningOption("log18", "log nbins=18", "log", 18),
        BinningOption("linear4", "linear dk=4", "linear", 4),
        BinningOption("linear8", "linear dk=8", "linear", 8),
        BinningOption("adaptive500", "adaptive min=500", "adaptive", 500),
        BinningOption("adaptive1000", "adaptive min=1000", "adaptive", 1000),
        BinningOption("adaptive2000", "adaptive min=2000", "adaptive", 2000),
    ]

    results = [_evaluate_option(option, kvals, spectrum, counts, kmax=kmax) for option in options]
    best = _pick_best(results)

    _save_summary_csv(args.output_dir / "velocity_binning_summary.csv", results, str(best["key"]))
    _make_plot(
        args.output_dir / "velocity_binning_comparison.png",
        kvals,
        spectrum,
        results,
        str(best["key"]),
        kmax=kmax,
    )

    print(f"recommended={best['key']}")
    for item in results:
        print(
            f"{item['key']}: bins_lowk={int(item['bins_lowk'])} "
            f"mean_width_lowk={float(item['mean_width_lowk']):.2f} "
            f"median_modes_lowk={float(item['median_modes_lowk']):.1f} "
            f"roughness_lowk={float(item['roughness_lowk']):.6f}"
        )


if __name__ == "__main__":
    main()
