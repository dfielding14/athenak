#!/usr/bin/env python3
"""Plot thermodynamic histories from the 2D Ito-tracer thermal-instability test."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as colors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from read_prtcl_thermo_history import read_history  # noqa: E402


def ism_cooling(temp: np.ndarray) -> np.ndarray:
    """Koyama-Inutsuka branch used by ISMCoolFn below 10^4.2 K."""
    return (
        2.0e-19 * np.exp(-1.184e5 / (temp + 1.0e3))
        + 2.8e-28 * np.sqrt(temp) * np.exp(-92.0 / temp)
    )


def tag_history(
    data: dict[str, np.ndarray], tag: int, field: str
) -> tuple[np.ndarray, np.ndarray]:
    selected = data["tag"] == tag
    order = np.argsort(data["time"][selected])
    return data["time"][selected][order], data[field][selected][order]


def select_sample_tags(
    data: dict[str, np.ndarray], final: np.ndarray, count: int
) -> np.ndarray:
    order = np.argsort(data["temperature"][final])
    tags = data["tag"][final][order]
    return tags[np.linspace(0, tags.size - 1, count, dtype=int)]


def plot_history(
    history_path: Path,
    output_path: Path,
    temperature_unit: float,
    density_unit: float,
    time_unit_myr: float,
    heating_rate: float,
    initial_pressure: float,
) -> None:
    data = read_history(history_path)
    cycles = np.unique(data["cycle"])
    final = data["cycle"] == cycles[-1]
    final_temperature = data["temperature"][final] * temperature_unit
    final_density = data["density"][final] * density_unit
    final_pressure = np.median(final_density * final_temperature)
    sample_tags = select_sample_tags(data, final, 28)

    times = np.array(
        [np.mean(data["time"][data["cycle"] == cycle]) for cycle in cycles]
    )
    temperature_percentiles = np.array(
        [
            np.percentile(
                data["temperature"][data["cycle"] == cycle] * temperature_unit,
                [10.0, 50.0, 90.0],
            )
            for cycle in cycles
        ]
    )

    norm = colors.LogNorm(vmin=30.0, vmax=8000.0)
    cmap = plt.get_cmap("coolwarm")
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.5))

    for tag in sample_tags:
        time, temperature = tag_history(data, int(tag), "temperature")
        tag_final = temperature[-1] * temperature_unit
        axes[0].plot(
            time * time_unit_myr,
            temperature * temperature_unit,
            color=cmap(norm(tag_final)),
            alpha=0.62,
            linewidth=0.75,
        )
    axes[0].fill_between(
        times * time_unit_myr,
        temperature_percentiles[:, 0],
        temperature_percentiles[:, 2],
        color="0.4",
        alpha=0.16,
        label="10th-90th percentile",
    )
    axes[0].plot(
        times * time_unit_myr,
        temperature_percentiles[:, 1],
        color="black",
        linewidth=1.6,
        label="median",
    )
    axes[0].axhline(temperature_unit, color="0.35", linestyle=":", linewidth=1.0)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("time [Myr]")
    axes[0].set_ylabel("particle temperature [K]")
    axes[0].set_title("Individual thermodynamic histories")
    axes[0].legend(frameon=False, loc="best")

    spatial = axes[1].scatter(
        data["x1"][final],
        data["x2"][final],
        c=final_temperature,
        cmap=cmap,
        norm=norm,
        s=5,
        linewidths=0,
        alpha=0.8,
    )
    axes[1].set_aspect("equal")
    axes[1].set_xlabel(r"$x_1$ [pc]")
    axes[1].set_ylabel(r"$x_2$ [pc]")
    axes[1].set_title("Final mass-tracer temperatures")
    colorbar = fig.colorbar(spatial, ax=axes[1], pad=0.02)
    colorbar.set_label("temperature [K]")

    axes[2].scatter(
        final_density,
        final_temperature,
        c=final_temperature,
        cmap=cmap,
        norm=norm,
        s=5,
        linewidths=0,
        alpha=0.16,
    )
    for tag in sample_tags:
        _, density = tag_history(data, int(tag), "density")
        _, temperature = tag_history(data, int(tag), "temperature")
        tag_final = temperature[-1] * temperature_unit
        axes[2].plot(
            density * density_unit,
            temperature * temperature_unit,
            color=cmap(norm(tag_final)),
            alpha=0.45,
            linewidth=0.65,
        )

    equilibrium_temperature = np.geomspace(20.0, 1.2e4, 600)
    equilibrium_density = heating_rate / ism_cooling(equilibrium_temperature)
    axes[2].plot(
        equilibrium_density,
        equilibrium_temperature,
        color="black",
        linewidth=1.5,
        label="heating-cooling equilibrium",
    )
    axes[2].plot(
        initial_pressure / equilibrium_temperature,
        equilibrium_temperature,
        color="0.35",
        linestyle="--",
        linewidth=1.0,
        label=rf"initial $P/k_B={initial_pressure:.0f}\ {{\rm K\,cm^{{-3}}}}$",
    )
    axes[2].plot(
        final_pressure / equilibrium_temperature,
        equilibrium_temperature,
        color="0.35",
        linestyle=":",
        linewidth=1.2,
        label=rf"final median $P/k_B={final_pressure:.0f}\ {{\rm K\,cm^{{-3}}}}$",
    )
    axes[2].set_xscale("log")
    axes[2].set_yscale("log")
    axes[2].set_xlim(0.15, 150.0)
    axes[2].set_ylim(20.0, 1.2e4)
    axes[2].set_xlabel(r"number density [cm$^{-3}$]")
    axes[2].set_ylabel("temperature [K]")
    axes[2].set_title("Tracer paths through phase space")
    axes[2].legend(frameon=False, fontsize=8, loc="lower left")

    fig.suptitle("2D ISM thermal instability traced by Ito-2 particles")
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("history", type=Path, help="particle thermo-history .thp file")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "docs/source/modules/figures/ito_tracers_thermal_instability_2d.png"
        ),
    )
    parser.add_argument("--temperature-unit", type=float, default=1573.006766288556)
    parser.add_argument("--density-unit", type=float, default=1.907175521614809)
    parser.add_argument("--time-unit-myr", type=float, default=0.3046955159499056)
    parser.add_argument("--heating-rate", type=float, default=1.9976230136434696e-26)
    parser.add_argument("--initial-pressure", type=float, default=3000.0)
    args = parser.parse_args()

    plot_history(
        args.history,
        args.output,
        args.temperature_unit,
        args.density_unit,
        args.time_unit_myr,
        args.heating_rate,
        args.initial_pressure,
    )


if __name__ == "__main__":
    main()
