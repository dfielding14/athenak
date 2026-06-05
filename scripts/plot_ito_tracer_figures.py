#!/usr/bin/env python3
"""Generate documentation figures from the Ito-2 tracer visualization runs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
import numpy as np  # noqa: E402

from read_prtcl_thermo_history import read_history  # noqa: E402


def _cycle_values(
    data: dict[str, np.ndarray], name: str
) -> tuple[np.ndarray, np.ndarray]:
    cycles = np.unique(data["cycle"])
    values = np.array([np.mean(data[name][data["cycle"] == cycle]) for cycle in cycles])
    return cycles, values


def _tag_history(
    data: dict[str, np.ndarray], tag: int, name: str
) -> tuple[np.ndarray, np.ndarray]:
    selected = data["tag"] == tag
    order = np.argsort(data["cycle"][selected])
    return data["time"][selected][order], data[name][selected][order]


def plot_1d(history_path: Path, output_path: Path) -> None:
    data = read_history(history_path)
    cycles = np.unique(data["cycle"])
    initial = data["cycle"] == cycles[0]
    final = data["cycle"] == cycles[-1]
    initial_by_tag = {
        int(tag): x1 for tag, x1 in zip(data["tag"][initial], data["x1"][initial])
    }
    displacement = np.array(
        [
            x1 - initial_by_tag[int(tag)]
            for tag, x1 in zip(data["tag"][final], data["x1"][final])
        ]
    )
    displacement = (displacement + 0.5) % 1.0 - 0.5

    _, times = _cycle_values(data, "time")
    _, velocity = _cycle_values(data, "v1")
    dt = np.diff(times)
    step_velocity = velocity[:-1]
    cell_width = 1.0 / 128.0
    courant = np.abs(step_velocity) * dt / cell_width
    expected_mean = np.concatenate(([0.0], np.cumsum(step_velocity * dt)))
    expected_variance = np.concatenate(
        ([0.0], np.cumsum(cell_width**2 * courant * (1.0 - courant)))
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    tags = np.unique(data["tag"])
    sample_tags = tags[np.linspace(0, tags.size - 1, 16, dtype=int)]
    initial_center = np.mean(data["x1"][initial])
    for tag in sample_tags:
        time, x1 = _tag_history(data, int(tag), "x1")
        axes[0].plot(time, x1, color="#1677b8", alpha=0.55, linewidth=0.8)
    axes[0].plot(
        times,
        initial_center + expected_mean,
        color="black",
        linestyle="--",
        linewidth=1.5,
        label="expected drift",
    )
    sigma = np.sqrt(expected_variance)
    axes[0].fill_between(
        times,
        initial_center + expected_mean - sigma,
        initial_center + expected_mean + sigma,
        color="#f2b134",
        alpha=0.3,
        label=r"expected $1\sigma$",
    )
    axes[0].set_xlabel("time")
    axes[0].set_ylabel(r"$x_1$")
    axes[0].set_title("Continuous sample trajectories")
    axes[0].legend(frameon=False, loc="upper left")

    bins = np.linspace(displacement.min(), displacement.max(), 38)
    axes[1].hist(
        displacement,
        bins=bins,
        density=True,
        color="#1677b8",
        alpha=0.65,
        label="Ito-2 particles",
    )
    final_mean = expected_mean[-1]
    final_sigma = sigma[-1]
    x = np.linspace(displacement.min(), displacement.max(), 400)
    normal = np.exp(-0.5 * ((x - final_mean) / final_sigma) ** 2)
    normal /= np.sqrt(2.0 * np.pi) * final_sigma
    axes[1].plot(
        x,
        normal,
        color="black",
        linestyle="--",
        linewidth=1.5,
        label="moment-matched normal guide",
    )
    axes[1].axvline(final_mean, color="#f2b134", linewidth=1.5, label="expected mean")
    axes[1].set_xlabel(r"$x_1(t_{\rm final})-x_1(0)$")
    axes[1].set_ylabel("probability density")
    axes[1].set_title("Final sheet displacement")
    axes[1].legend(frameon=False)

    fig.suptitle("1D-style uniform advection: Ito-2 drift and diffusion")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _set_grid(ax: plt.Axes, cell_width: float) -> None:
    ax.xaxis.set_minor_locator(MultipleLocator(cell_width))
    ax.yaxis.set_minor_locator(MultipleLocator(cell_width))
    ax.grid(which="minor", color="0.88", linewidth=0.45)
    ax.grid(which="major", color="0.75", linewidth=0.7)


def plot_2d(history_path: Path, output_path: Path) -> None:
    data = read_history(history_path)
    cycles = np.unique(data["cycle"])
    initial = data["cycle"] == cycles[0]
    final = data["cycle"] == cycles[-1]
    initial_center = np.array(
        [np.mean(data["x1"][initial]), np.mean(data["x2"][initial])]
    )
    final_center = np.array([np.mean(data["x1"][final]), np.mean(data["x2"][final])])

    all_x = np.concatenate((data["x1"][initial], data["x1"][final]))
    all_y = np.concatenate((data["x2"][initial], data["x2"][final]))
    margin = 0.025
    limits = (
        (all_x.min() - margin, all_x.max() + margin),
        (all_y.min() - margin, all_y.max() + margin),
    )
    cell_width = 1.0 / 64.0

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8), sharex=True, sharey=True)

    axes[0].scatter(
        data["x1"][initial],
        data["x2"][initial],
        s=12,
        color="0.35",
        alpha=0.55,
        linewidths=0,
        label="initial cloud",
    )
    axes[0].scatter(
        data["x1"][final],
        data["x2"][final],
        s=4,
        color="#1677b8",
        alpha=0.25,
        linewidths=0,
        label="final cloud",
    )
    axes[0].annotate(
        "",
        xy=final_center,
        xytext=initial_center,
        arrowprops={"arrowstyle": "->", "color": "#d95f02", "linewidth": 2.0},
    )
    axes[0].scatter(*initial_center, marker="x", color="black", s=35, zorder=3)
    axes[0].scatter(*final_center, marker="x", color="#d95f02", s=35, zorder=3)
    axes[0].set_title("Initial and final particle clouds")
    axes[0].legend(frameon=False, loc="upper left")

    tags = np.unique(data["tag"])
    sample_tags = tags[np.linspace(0, tags.size - 1, 48, dtype=int)]
    for tag in sample_tags:
        _, x1 = _tag_history(data, int(tag), "x1")
        _, x2 = _tag_history(data, int(tag), "x2")
        axes[1].plot(x1, x2, color="#1677b8", alpha=0.55, linewidth=0.75)
        axes[1].scatter(x1[0], x2[0], s=5, color="0.25", linewidths=0)
        axes[1].scatter(x1[-1], x2[-1], s=5, color="#d95f02", linewidths=0)
    axes[1].set_title("Continuous paths across cell boundaries")

    for ax in axes:
        ax.set_xlim(*limits[0])
        ax.set_ylim(*limits[1])
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x_1$")
        _set_grid(ax, cell_width)
    axes[0].set_ylabel(r"$x_2$")

    fig.suptitle("2D diagonal advection: continuous Ito-2 cloud transport")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--one-d", type=Path, required=True, help="1D-style .thp file")
    parser.add_argument("--two-d", type=Path, required=True, help="2D cloud .thp file")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("docs/source/modules/figures"),
        help="directory for generated PNG files",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_1d(args.one_d, args.output_dir / "ito_tracers_1d_sheet.png")
    plot_2d(args.two_d, args.output_dir / "ito_tracers_2d_cloud.png")


if __name__ == "__main__":
    main()
