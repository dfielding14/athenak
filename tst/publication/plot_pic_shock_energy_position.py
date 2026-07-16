#!/usr/bin/env python3
"""Plot the CR energy-position distribution for one shock snapshot."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import colors  # noqa: E402
import numpy as np  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.publication.analyze_q011_section54_outputs import (  # noqa: E402
    compose_leaf_field,
    read_athenak_binary,
)
from tst.publication.plot_pic_shock_cr_variants import (  # noqa: E402
    _particle_kinematics,
)
from tst.publication.plot_pic_shock_nine_panel import (  # noqa: E402
    _parameter,
    _particle_time,
    _shock_front,
    _style,
)
from tst.publication.pvtk_particles import read_particle_vtk  # noqa: E402


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mhd", type=Path, required=True)
    parser.add_argument("--particles", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--downstream-width", type=float, default=6000.0)
    parser.add_argument("--upstream-width", type=float, default=12000.0)
    parser.add_argument("--x-bins", type=int, default=1500)
    parser.add_argument("--energy-bins", type=int, default=120)
    parser.add_argument("--energy-min", type=float, default=0.3)
    parser.add_argument("--energy-max", type=float, default=1000.0)
    parser.add_argument(
        "--birth-cuts", type=float, nargs="*", default=None,
        help="also plot particles born at or after each supplied time",
    )
    parser.add_argument("--dpi", type=int, default=350)
    return parser.parse_args()


def _reduce(args: argparse.Namespace) -> dict[str, np.ndarray | float]:
    dataset = read_athenak_binary(args.mhd.resolve(strict=True))
    density_grid = compose_leaf_field(dataset, "dens")
    density = np.asarray(density_grid.values[0], dtype=np.float64)
    x1_faces = np.asarray(density_grid.x1_faces, dtype=np.float64)
    x1 = 0.5 * (x1_faces[:-1] + x1_faces[1:])
    shock = _shock_front(x1, density)
    parameters = {
        block: dict(values) for block, values in dataset.input_parameters.items()
    }
    rho0 = _parameter(parameters, "problem", "ps_rho0", 1.0)
    b0 = _parameter(parameters, "problem", "ps_b0", 1.0)
    u0 = _parameter(parameters, "problem", "ps_u0", 30.0)
    gamma_gas = _parameter(parameters, "mhd", "gamma", 5.0 / 3.0)
    light_speed = _parameter(parameters, "particles", "pic_cr_light_speed", 1.0)
    alpha_i = _parameter(
        parameters, "particles", "pic_background_ion_q_over_mc", 1.0
    )
    pinj = _parameter(
        parameters, "problem", "ps_vinj_over_u0", np.sqrt(10.0)
    ) * u0
    shock_speed = 0.5 * (gamma_gas - 1.0) * u0
    length_scale = alpha_i * np.sqrt(rho0)
    omega0 = alpha_i * b0

    particles = read_particle_vtk(args.particles.resolve(strict=True))
    points = np.asarray(particles.points, dtype=np.float64)
    velocity = np.asarray(particles.vectors["vel"], dtype=np.float64)
    weights = np.asarray(particles.scalars["macro_weight"], dtype=np.float64)
    birth_time = np.asarray(particles.scalars["birth_time"], dtype=np.float64)
    source = np.asarray(
        particles.scalars.get("cr_source", np.ones(len(points))), dtype=np.float64
    )
    _, kinetic_energy = _particle_kinematics(
        velocity, light_speed, shock_speed, pinj
    )
    energy_ratio = kinetic_energy / (0.5 * u0**2)
    x_relative = (points[:, 0] - shock) * length_scale
    selected = (
        (source == 1)
        & np.isfinite(energy_ratio)
        & (energy_ratio > 0.0)
        & (x_relative >= -args.downstream_width * length_scale)
        & (x_relative <= args.upstream_width * length_scale)
    )
    x_edges = np.linspace(
        -args.downstream_width * length_scale,
        args.upstream_width * length_scale,
        args.x_bins + 1,
    )
    log_energy_edges = np.linspace(
        np.log10(args.energy_min), np.log10(args.energy_max),
        args.energy_bins + 1,
    )
    dx = np.diff(x_edges)[:, None]
    dlog_energy = np.diff(log_energy_edges)[None, :]
    birth_cuts = np.asarray(
        [np.nan] + ([] if args.birth_cuts is None else args.birth_cuts),
        dtype=np.float64,
    )
    distributions = []
    selected_weights = []
    log_energy_ratio = np.log10(energy_ratio)
    for birth_cut in birth_cuts:
        cohort = selected if np.isnan(birth_cut) else selected & (birth_time >= birth_cut)
        histogram, _, _ = np.histogram2d(
            x_relative[cohort], log_energy_ratio[cohort],
            bins=(x_edges, log_energy_edges), weights=weights[cohort],
        )
        distributions.append(histogram / (dx * dlog_energy))
        selected_weights.append(float(np.sum(weights[cohort])))
    distribution = np.stack(distributions)
    if args.birth_cuts is None:
        distribution = distribution[0]
    return {
        "distribution": distribution,
        "birth_cuts": birth_cuts,
        "selected_weights": np.asarray(selected_weights),
        "x_edges": x_edges,
        "log_energy_edges": log_energy_edges,
        "time": _particle_time(args.particles) * omega0,
        "shock": shock * length_scale,
    }


def main() -> int:
    args = _arguments()
    cache = args.cache.resolve()
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.is_file():
        with np.load(cache) as saved:
            reduced = {name: np.asarray(saved[name]) for name in saved.files}
    else:
        reduced = _reduce(args)
        np.savez_compressed(cache, **reduced)

    distribution = np.asarray(reduced["distribution"], dtype=np.float64)
    x_edges = np.asarray(reduced["x_edges"], dtype=np.float64)
    energy_edges = np.asarray(reduced["log_energy_edges"], dtype=np.float64)
    time = float(np.asarray(reduced["time"]))
    normalization_distribution = (
        distribution[0] if distribution.ndim == 3 else distribution
    )
    positive = normalization_distribution[normalization_distribution > 0.0]
    norm = colors.LogNorm(
        vmin=max(float(np.percentile(positive, 1.0)), np.finfo(float).tiny),
        vmax=float(np.percentile(positive, 99.8)),
    )

    _style()
    if distribution.ndim == 2:
        panels = distribution[None, ...]
        labels = [None]
        figure, axes = plt.subplots(figsize=(10.2, 3.5), squeeze=False)
        axes = axes[0]
    else:
        panels = distribution
        birth_cuts = np.asarray(reduced["birth_cuts"], dtype=np.float64)
        labels = ["all particles"] + [
            rf"$t_{{\rm birth}}\geq {cut:g}$" for cut in birth_cuts[1:]
        ]
        figure, axes = plt.subplots(
            1, len(panels), figsize=(5.0 * len(panels), 3.5),
            sharex=True, sharey=True, squeeze=False,
        )
        axes = axes[0]

    for index, (axis, panel, label) in enumerate(zip(axes, panels, labels)):
        image = axis.pcolormesh(
            x_edges, energy_edges, panel.T, shading="flat",
            cmap="turbo", norm=norm, rasterized=True,
        )
        axis.axvline(0.0, color="white", lw=1.0, ls="--")
        axis.set_xlabel(r"$(x_1-x_{\rm sh})\,\omega_{pi}/c$")
        if index == 0:
            axis.set_ylabel(r"$\log_{10}(E/E_{\rm sh})$")
            axis.text(
                0.02, 0.96, rf"$\Omega_0 t={time:.1f}$", transform=axis.transAxes,
                ha="left", va="top", fontsize=10,
                bbox={"facecolor": "white", "alpha": 0.82, "pad": 2.0,
                      "edgecolor": "none"},
            )
        if label is not None:
            axis.set_title(label, fontsize=10)
        axis.text(0.02, 0.05, "downstream", transform=axis.transAxes,
                  ha="left", va="bottom", color="white", fontsize=8)
        axis.text(0.98, 0.05, "upstream", transform=axis.transAxes,
                  ha="right", va="bottom", color="white", fontsize=8)
    bar = figure.colorbar(image, ax=list(axes), pad=0.02)
    bar.set_label(r"$E f(E)$")
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(figure)
    print(f"time={time:g}")
    print(f"cache={cache}")
    print(f"figure={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
