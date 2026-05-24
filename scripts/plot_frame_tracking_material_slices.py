#!/usr/bin/env python3
"""Plot final tracer-density slices from cloud and TRML material validation runs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
import athena_read  # type: ignore[import-not-found]  # noqa: E402
import bin_convert as bin_reader  # type: ignore[import-not-found]  # noqa: E402


def final_snapshot(run_dir: Path) -> Path:
    paths = sorted((run_dir / "bin").glob("*.bin"))
    if not paths:
        raise FileNotFoundError(f"No binary snapshot found below {run_dir}.")
    return paths[-1]


def displacement(run_dir: Path, axis: str) -> float:
    paths = sorted(run_dir.glob("*.frame_tracker.hst"))
    if not paths:
        return 0.0
    data = athena_read.hst(str(paths[-1]))
    return float(np.asarray(data[f"ft_dx_{axis}"])[-1])


def read_tracer_density(run_dir: Path) -> dict[str, np.ndarray]:
    data = bin_reader.read_binary_as_athdf(
        str(final_snapshot(run_dir)), quantities=["dens", "s_00"], dtype=np.float64
    )
    return {
        "x1f": np.asarray(data["x1f"], dtype=float),
        "x2f": np.asarray(data["x2f"], dtype=float),
        "x3f": np.asarray(data["x3f"], dtype=float),
        "tracer": np.asarray(data["dens"], dtype=float)
        * np.maximum(np.asarray(data["s_00"], dtype=float), 0.0),
    }


def plot_pair(
    axes: tuple[plt.Axes, plt.Axes],
    off_dir: Path,
    on_dir: Path,
    problem: str,
    tracking_axis: str,
) -> None:
    off = read_tracer_density(off_dir)
    on = read_tracer_density(on_dir)
    if problem == "Cloud":
        off_slice = off["tracer"][off["tracer"].shape[0] // 2, :, :]
        on_slice = on["tracer"][on["tracer"].shape[0] // 2, :, :]
        off_extent = [off["x1f"][0], off["x1f"][-1], off["x2f"][0], off["x2f"][-1]]
        on_extent = [
            on["x1f"][0] + displacement(on_dir, tracking_axis),
            on["x1f"][-1] + displacement(on_dir, tracking_axis),
            on["x2f"][0],
            on["x2f"][-1],
        ]
        ylabel = "lab x2"
    else:
        off_slice = off["tracer"][:, off["tracer"].shape[1] // 2, :]
        on_slice = on["tracer"][:, on["tracer"].shape[1] // 2, :]
        off_extent = [off["x1f"][0], off["x1f"][-1], off["x3f"][0], off["x3f"][-1]]
        on_extent = [
            on["x1f"][0],
            on["x1f"][-1],
            on["x3f"][0] + displacement(on_dir, tracking_axis),
            on["x3f"][-1] + displacement(on_dir, tracking_axis),
        ]
        ylabel = "lab x3"
    vmax = max(float(np.max(off_slice)), float(np.max(on_slice)))
    for axis, image, extent, state in zip(
        axes, (off_slice, on_slice), (off_extent, on_extent), ("disabled", "enabled")
    ):
        plotted = axis.imshow(
            image,
            origin="lower",
            extent=extent,
            cmap="viridis",
            vmin=0.0,
            vmax=vmax,
            aspect="auto",
        )
        axis.set_title(f"{problem}: tracking {state}")
        axis.set_xlabel("lab x1")
        axis.set_ylabel(ylabel)
    plt.colorbar(plotted, ax=axes, label="tracer density", shrink=0.82)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cloud-off", type=Path, required=True)
    parser.add_argument("--cloud-on", type=Path, required=True)
    parser.add_argument("--trml-off", type=Path, required=True)
    parser.add_argument("--trml-on", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    figure, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
    plot_pair(tuple(axes[0]), args.cloud_off, args.cloud_on, "Cloud", "x1")
    plot_pair(tuple(axes[1]), args.trml_off, args.trml_on, "TRML", "x3")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
