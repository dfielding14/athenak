#!/usr/bin/env python3
"""Plot density, transverse magnetic field, and CR acceleration at three times."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import math
from pathlib import Path
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import colors  # noqa: E402
import cmasher as cmr  # noqa: E402
import numpy as np  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.publication.analyze_q011_section54_outputs import (  # noqa: E402
    compose_leaf_field,
    read_athenak_binary,
)
from tst.publication.pvtk_particles import read_particle_vtk  # noqa: E402


PVTK_HEADER = re.compile(
    rb"# AthenaK particle data at time= ([^ ]+)  nranks= ([0-9]+)"
    rb"  cycle=([0-9]+)  variables=prtcl_all"
)


def _binary_time(path: Path) -> float:
    with path.open("rb") as stream:
        for _ in range(8):
            line = stream.readline()
            if line.startswith(b"  time="):
                return float(line.split(b"=", 1)[1])
    raise RuntimeError(f"{path}: binary header has no time")


def _particle_time(path: Path) -> float:
    with path.open("rb") as stream:
        header = stream.read(512)
    match = PVTK_HEADER.search(header)
    if match is None:
        raise RuntimeError(f"{path}: malformed particle header")
    return float(match.group(1))


def _files(run: Path, directory: str, pattern: str) -> list[Path]:
    paths = sorted((run / directory).glob(pattern))
    if not paths:
        paths = sorted(run.glob(pattern))
    return paths


def _nearest(items: list[tuple[float, Path]], target: float) -> tuple[float, Path]:
    return min(items, key=lambda item: abs(item[0] - target))


def _parameter(
    parameters: dict[str, dict[str, str]], block: str, name: str, default: float
) -> float:
    return float(parameters.get(block, {}).get(name, str(default)))


def _shock_front(x1: np.ndarray, density: np.ndarray) -> float:
    profile = np.mean(density, axis=0)
    gradient = np.gradient(profile, x1)
    span = float(x1[-1] - x1[0])
    interior = np.flatnonzero(
        (x1 > x1[0] + 0.05 * span) & (x1 < x1[-1] - 0.05 * span)
    )
    return float(x1[interior[np.argmin(gradient[interior])]])


def _coarsened_edges(edges: np.ndarray, factor: int) -> np.ndarray:
    bins = math.ceil((len(edges) - 1) / factor)
    return np.linspace(float(edges[0]), float(edges[-1]), bins + 1)


def _momentum_ratio(
    velocity: np.ndarray, light_speed: float, shock_speed: float, pinj: float
) -> np.ndarray:
    denominator = 1.0 - shock_speed * velocity[:, 0] / light_speed**2
    gamma_surface = 1.0 / math.sqrt(1.0 - (shock_speed / light_speed) ** 2)
    relative = velocity.copy()
    relative[:, 0] = (velocity[:, 0] - shock_speed) / denominator
    relative[:, 1:] = velocity[:, 1:] / (gamma_surface * denominator[:, None])
    speed2 = np.sum(relative**2, axis=1)
    gamma = 1.0 / np.sqrt(1.0 - speed2 / light_speed**2)
    return gamma * np.sqrt(speed2) / pinj


def _acceleration_map(
    particle_path: Path,
    x2_edges: np.ndarray,
    x1_edges: np.ndarray,
    light_speed: float,
    shock_speed: float,
    pinj: float,
) -> np.ndarray:
    particles = read_particle_vtk(particle_path)
    points = np.asarray(particles.points, dtype=np.float64)
    velocity = np.asarray(particles.vectors["vel"], dtype=np.float64)
    weights = np.asarray(particles.scalars["macro_weight"], dtype=np.float64)
    source = np.asarray(particles.scalars.get("cr_source", np.ones(len(points))))
    selected = source == 1
    points = points[selected]
    weights = weights[selected]
    ratio = _momentum_ratio(
        velocity[selected], light_speed, shock_speed, pinj
    )

    dx2 = float(x2_edges[1] - x2_edges[0])
    dx1 = float(x1_edges[1] - x1_edges[0])
    if not np.allclose(np.diff(x2_edges), dx2) or not np.allclose(
        np.diff(x1_edges), dx1
    ):
        raise RuntimeError("fast CR binning requires uniform output coordinates")
    nx2 = len(x2_edges) - 1
    nx1 = len(x1_edges) - 1
    i2 = np.floor((points[:, 1] - x2_edges[0]) / dx2).astype(np.int64)
    i1 = np.floor((points[:, 0] - x1_edges[0]) / dx1).astype(np.int64)
    i2[points[:, 1] == x2_edges[-1]] = nx2 - 1
    i1[points[:, 0] == x1_edges[-1]] = nx1 - 1
    inside = (i2 >= 0) & (i2 < nx2) & (i1 >= 0) & (i1 < nx1)
    flat_index = i2[inside] * nx1 + i1[inside]
    weight_sum = np.bincount(
        flat_index, weights=weights[inside], minlength=nx2 * nx1
    ).reshape(nx2, nx1)
    momentum_sum = np.bincount(
        flat_index,
        weights=weights[inside] * ratio[inside],
        minlength=nx2 * nx1,
    ).reshape(nx2, nx1)
    mean_ratio = np.full_like(weight_sum, np.nan)
    np.divide(momentum_sum, weight_sum, out=mean_ratio, where=weight_sum > 0.0)
    return mean_ratio - 1.0


def _acceleration_task(arguments: tuple[object, ...]) -> np.ndarray:
    return _acceleration_map(*arguments)


def _style() -> None:
    plt.style.use("default")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "#eeeeee",
            "axes.edgecolor": "#202020",
            "axes.linewidth": 0.8,
            "axes.unicode_minus": False,
            "font.family": "serif",
            "font.serif": ["cmr10"],
            "mathtext.fontset": "cm",
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
        }
    )


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a horizontal three-time, nine-panel PIC shock figure."
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        action="append",
        required=True,
        help="Run output directory; repeat to combine a run and its continuation.",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fractions", type=float, nargs=3, default=(0.34, 0.66, 1.0))
    parser.add_argument("--upstream-width", type=float, default=1080.0)
    parser.add_argument("--downstream-width", type=float, default=720.0)
    parser.add_argument(
        "--x1-frame", choices=("shock", "absolute"), default="shock"
    )
    parser.add_argument("--cr-bin-factor", type=int, default=2)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--refresh-cache", action="store_true")
    parser.add_argument("--dpi", type=int, default=350)
    return parser.parse_args()


def main() -> int:
    args = _arguments()
    if args.cr_bin_factor < 1:
        raise RuntimeError("--cr-bin-factor must be positive")
    if args.workers < 1:
        raise RuntimeError("--workers must be positive")
    if any(value <= 0.0 or value > 1.0 for value in args.fractions):
        raise RuntimeError("--fractions must lie in (0, 1]")

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    cache_dir = (
        args.cache_dir.resolve()
        if args.cache_dir is not None
        else output.parent / "processed-particles"
    )
    cache_dir.mkdir(parents=True, exist_ok=True)

    runs = [path.resolve(strict=True) for path in args.run_dir]
    mhd_paths = sorted(
        {path for run in runs for path in _files(run, "bin", "*.mhd_w_bcc.*.bin")}
    )
    particle_paths = sorted(
        {path for run in runs for path in _files(run, "pvtk", "*.prtcl_all.*.part.vtk")}
    )
    if not mhd_paths or not particle_paths:
        raise RuntimeError(
            "run directories must contain MHD binary and particle VTK outputs"
        )

    mhd_timeline = sorted((_binary_time(path), path) for path in mhd_paths)
    particle_timeline = sorted((_particle_time(path), path) for path in particle_paths)
    latest = particle_timeline[-1][0]
    particle_selection = [
        _nearest(particle_timeline, fraction * latest) for fraction in args.fractions
    ]
    if len({path for _, path in particle_selection}) != 3:
        raise RuntimeError(
            "selected fractions do not resolve to three distinct snapshots"
        )
    mhd_selection = [_nearest(mhd_timeline, time) for time, _ in particle_selection]

    snapshots: list[dict[str, object]] = []
    acceleration_tasks: list[tuple[object, ...]] = []
    cache_paths: list[Path] = []
    for (particle_time, particle_path), (mhd_time, mhd_path) in zip(
        particle_selection, mhd_selection
    ):
        dataset = read_athenak_binary(mhd_path)
        density_grid = compose_leaf_field(dataset, "dens")
        by_grid = compose_leaf_field(dataset, "bcc2")
        bz_grid = compose_leaf_field(dataset, "bcc3")
        density = np.asarray(density_grid.values[0], dtype=np.float64)
        by = np.asarray(by_grid.values[0], dtype=np.float64)
        bz = np.asarray(bz_grid.values[0], dtype=np.float64)
        x1_faces = np.asarray(density_grid.x1_faces, dtype=np.float64)
        x2_faces = np.asarray(density_grid.x2_faces, dtype=np.float64)
        x1 = 0.5 * (x1_faces[:-1] + x1_faces[1:])
        shock = _shock_front(x1, density)
        selected = np.flatnonzero(
            (x1 >= shock - args.downstream_width)
            & (x1 <= shock + args.upstream_width)
        )
        first, last = int(selected[0]), int(selected[-1])
        crop_faces = x1_faces[first : last + 2]
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
        vinj_over_u0 = _parameter(
            parameters, "problem", "ps_vinj_over_u0", math.sqrt(10.0)
        )
        shock_speed = 0.5 * (gamma_gas - 1.0) * u0
        omega_pi_over_c = alpha_i * math.sqrt(rho0)
        omega0 = alpha_i * b0
        if args.x1_frame == "shock":
            x1_origin = shock
            shock_coordinate = 0.0
            display_downstream = -args.downstream_width * omega_pi_over_c
            display_upstream = args.upstream_width * omega_pi_over_c
        else:
            x1_origin = 0.0
            shock_coordinate = shock * omega_pi_over_c
            display_downstream = (
                shock - args.downstream_width
            ) * omega_pi_over_c
            display_upstream = (shock + args.upstream_width) * omega_pi_over_c

        cr_x1_edges = _coarsened_edges(crop_faces, args.cr_bin_factor)
        cr_x2_edges = _coarsened_edges(x2_faces, args.cr_bin_factor)
        acceleration_tasks.append(
            (
                particle_path,
                cr_x2_edges,
                cr_x1_edges,
                light_speed,
                shock_speed,
                vinj_over_u0 * u0,
            )
        )
        cache_name = (
            f"{particle_path.name.removesuffix('.part.vtk')}"
            f".t{particle_time:g}.u{args.upstream_width:g}"
            f".d{args.downstream_width:g}.b{args.cr_bin_factor}.npz"
        )
        cache_paths.append(cache_dir / cache_name)
        snapshots.append(
            {
                "time": particle_time * omega0,
                "mhd_time": mhd_time,
                "particle_path": particle_path,
                "mhd_path": mhd_path,
                "x2_faces": x2_faces * omega_pi_over_c,
                "x1_faces": (crop_faces - x1_origin) * omega_pi_over_c,
                "cr_x2_faces": cr_x2_edges * omega_pi_over_c,
                "cr_x1_faces": (cr_x1_edges - x1_origin) * omega_pi_over_c,
                "shock_coordinate": shock_coordinate,
                "display_downstream": display_downstream,
                "display_upstream": display_upstream,
                "density": density[:, first : last + 1] / rho0,
                "bperp": np.sqrt(
                    by[:, first : last + 1] ** 2 + bz[:, first : last + 1] ** 2
                ) / b0,
            }
        )

    acceleration_maps: list[np.ndarray | None] = [None] * len(acceleration_tasks)
    missing: list[int] = []
    for index, cache_path in enumerate(cache_paths):
        if cache_path.is_file() and not args.refresh_cache:
            with np.load(cache_path) as cached:
                acceleration_maps[index] = np.asarray(
                    cached["acceleration"], dtype=np.float64
                )
        else:
            missing.append(index)
    missing_tasks = [acceleration_tasks[index] for index in missing]
    if not missing_tasks:
        reduced = []
    elif args.workers == 1:
        reduced = [_acceleration_map(*task) for task in missing_tasks]
    else:
        with ProcessPoolExecutor(
            max_workers=min(args.workers, len(missing_tasks))
        ) as executor:
            reduced = list(executor.map(_acceleration_task, missing_tasks))
    for index, acceleration in zip(missing, reduced):
        acceleration_maps[index] = acceleration
        np.savez_compressed(
            cache_paths[index],
            acceleration=acceleration,
            particle_time=float(particle_selection[index][0]),
        )
    for snapshot, acceleration in zip(snapshots, acceleration_maps):
        if acceleration is None:
            raise RuntimeError("particle reduction did not produce an acceleration map")
        snapshot["acceleration"] = acceleration

    rho_values = np.concatenate([np.ravel(item["density"]) for item in snapshots])
    b_values = np.concatenate([np.ravel(item["bperp"]) for item in snapshots])
    acceleration_values = np.concatenate(
        [np.ravel(item["acceleration"]) for item in snapshots]
    )
    positive_b = b_values[np.isfinite(b_values) & (b_values > 0.0)]
    finite_acceleration = acceleration_values[np.isfinite(acceleration_values)]
    rho_norm = colors.Normalize(
        vmin=max(0.0, float(np.percentile(rho_values, 0.3))),
        vmax=float(np.percentile(rho_values, 99.7)),
    )
    b_norm = colors.LogNorm(
        vmin=max(float(np.percentile(positive_b, 1.0)), 1.0e-5),
        vmax=float(np.percentile(positive_b, 99.7)),
    )
    acceleration_limit = max(
        float(np.percentile(np.abs(finite_acceleration), 99.0)), 0.005
    )
    acceleration_norm = colors.TwoSlopeNorm(
        vmin=-acceleration_limit, vcenter=0.0, vmax=acceleration_limit
    )

    _style()
    figure = plt.figure(figsize=(17.5, 10.2))
    layout = figure.add_gridspec(2, 1, height_ratios=(0.025, 1.0), hspace=0.075)
    colorbar_grid = layout[0].subgridspec(1, 3, wspace=0.13)
    groups = layout[1].subgridspec(1, 3, wspace=0.13)
    mappables: list[object] | None = None
    for group_index, snapshot in enumerate(snapshots):
        panels = groups[group_index].subgridspec(1, 3, wspace=0.0)
        axes = [figure.add_subplot(panels[0, index]) for index in range(3)]
        images = [
            axes[0].pcolormesh(
                snapshot["x2_faces"],
                snapshot["x1_faces"],
                np.asarray(snapshot["density"]).T,
                shading="flat",
                cmap=cmr.rainforest,
                norm=rho_norm,
                rasterized=True,
            ),
            axes[1].pcolormesh(
                snapshot["x2_faces"],
                snapshot["x1_faces"],
                np.asarray(snapshot["bperp"]).T,
                shading="flat",
                cmap=cmr.chroma,
                norm=b_norm,
                rasterized=True,
            ),
            axes[2].pcolormesh(
                snapshot["cr_x2_faces"],
                snapshot["cr_x1_faces"],
                np.asarray(snapshot["acceleration"]).T,
                shading="flat",
                cmap="RdBu_r",
                norm=acceleration_norm,
                rasterized=True,
            ),
        ]
        if mappables is None:
            mappables = images
        for column, axis in enumerate(axes):
            axis.axhline(
                float(snapshot["shock_coordinate"]),
                color="white",
                lw=0.8,
                ls="--",
                alpha=0.9,
            )
            axis.set_xlim(
                float(np.asarray(snapshot["x2_faces"])[0]),
                float(np.asarray(snapshot["x2_faces"])[-1]),
            )
            axis.set_ylim(
                float(snapshot["display_upstream"]),
                float(snapshot["display_downstream"]),
            )
            axis.set_xlabel(r"$x_2\,\omega_{pi}/c$", fontsize=9)
            if column == 0:
                if group_index == 0:
                    ylabel = (
                        r"$(x_1-x_{\rm sh})\,\omega_{pi}/c$"
                        if args.x1_frame == "shock"
                        else r"$x_1\,\omega_{pi}/c$"
                    )
                    axis.set_ylabel(ylabel, fontsize=10)
            else:
                axis.tick_params(labelleft=False)
        axes[1].text(
            0.5,
            0.975,
            rf"$\Omega_0 t={float(snapshot['time']):.1f}$",
            transform=axes[1].transAxes,
            ha="center",
            va="top",
            fontsize=10,
            color="black",
            bbox={
                "facecolor": "white",
                "alpha": 0.82,
                "pad": 2.0,
                "edgecolor": "none",
            },
        )

    colorbar_labels = (
        r"$\rho/\rho_0$",
        r"$B_\perp/B_0$",
        r"$\langle p\rangle/p_{\rm inj}-1$",
    )
    if mappables is None:
        raise RuntimeError("no plot panels were created")
    for index, (image, label) in enumerate(zip(mappables, colorbar_labels)):
        axis = figure.add_subplot(colorbar_grid[0, index])
        bar = figure.colorbar(image, cax=axis, orientation="horizontal")
        bar.set_label(label, fontsize=9)
        bar.ax.tick_params(labelsize=7.5)
        bar.ax.xaxis.set_ticks_position("top")
        bar.ax.xaxis.set_label_position("top")

    figure.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(figure)

    print(f"latest particle time: {latest:g}")
    for fraction, snapshot in zip(args.fractions, snapshots):
        print(
            f"fraction={fraction:g} time={snapshot['time']:g} "
            f"mhd_time={snapshot['mhd_time']:g}\n"
            f"  mhd={snapshot['mhd_path']}\n"
            f"  particles={snapshot['particle_path']}"
        )
    print(f"figure: {output}")
    print(f"processed-particle cache: {cache_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
