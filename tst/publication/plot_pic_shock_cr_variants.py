#!/usr/bin/env python3
"""Create nine-panel shock figures with several cached CR diagnostics."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import math
from pathlib import Path
import sys

import cmasher as cmr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import colors  # noqa: E402
import numpy as np  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.publication.plot_pic_shock_nine_panel import (  # noqa: E402
    _binary_time,
    _coarsened_edges,
    _files,
    _nearest,
    _parameter,
    _particle_time,
    _shock_front,
    _style,
)
from tst.publication.analyze_q011_section54_outputs import (  # noqa: E402
    compose_leaf_field,
    read_athenak_binary,
)
from tst.publication.pvtk_particles import read_particle_vtk  # noqa: E402


P_EDGES = np.geomspace(0.1, 100.0, 97)
TAIL_THRESHOLD = 1.5
VARIANTS = (
    "tail_energy",
    "tail_fraction",
    "p90",
    "p99",
    "energy_density",
    "momentum_phase",
    "sun_bai_energy_phase",
)


def _particle_kinematics(
    velocity: np.ndarray,
    light_speed: float,
    shock_speed: float,
    pinj: float,
) -> tuple[np.ndarray, np.ndarray]:
    denominator = 1.0 - shock_speed * velocity[:, 0] / light_speed**2
    gamma_surface = 1.0 / math.sqrt(1.0 - (shock_speed / light_speed) ** 2)
    relative = velocity.copy()
    relative[:, 0] = (velocity[:, 0] - shock_speed) / denominator
    relative[:, 1:] = velocity[:, 1:] / (gamma_surface * denominator[:, None])
    speed2 = np.sum(relative**2, axis=1)
    gamma = 1.0 / np.sqrt(1.0 - speed2 / light_speed**2)
    momentum_ratio = gamma * np.sqrt(speed2) / pinj
    kinetic_energy = (gamma - 1.0) * light_speed**2
    return momentum_ratio, kinetic_energy


def _particle_distribution(
    arguments: tuple[object, ...],
) -> tuple[np.ndarray, np.ndarray]:
    (
        particle_path,
        x2_edges,
        x1_edges,
        light_speed,
        shock_speed,
        pinj,
    ) = arguments
    particles = read_particle_vtk(Path(particle_path))
    points = np.asarray(particles.points, dtype=np.float64)
    source = np.asarray(particles.scalars.get("cr_source", np.ones(len(points))))
    weights = np.asarray(particles.scalars["macro_weight"], dtype=np.float64)
    velocity = np.asarray(particles.vectors["vel"], dtype=np.float64)

    dx2 = float(x2_edges[1] - x2_edges[0])
    dx1 = float(x1_edges[1] - x1_edges[0])
    nx2 = len(x2_edges) - 1
    nx1 = len(x1_edges) - 1
    i2 = np.floor((points[:, 1] - x2_edges[0]) / dx2).astype(np.int64)
    i1 = np.floor((points[:, 0] - x1_edges[0]) / dx1).astype(np.int64)
    i2[points[:, 1] == x2_edges[-1]] = nx2 - 1
    i1[points[:, 0] == x1_edges[-1]] = nx1 - 1
    selected = (
        (source == 1)
        & (i2 >= 0)
        & (i2 < nx2)
        & (i1 >= 0)
        & (i1 < nx1)
    )
    i2 = i2[selected]
    i1 = i1[selected]
    weights = weights[selected]
    ratio, kinetic_energy = _particle_kinematics(
        velocity[selected], float(light_speed), float(shock_speed), float(pinj)
    )
    ip = np.searchsorted(P_EDGES, ratio, side="right") - 1
    inside = (ip >= 0) & (ip < len(P_EDGES) - 1)
    npbin = len(P_EDGES) - 1
    flat = ((i2[inside] * nx1 + i1[inside]) * npbin + ip[inside])
    size = nx2 * nx1 * npbin
    weight_hist = np.bincount(
        flat, weights=weights[inside], minlength=size
    ).reshape(nx2, nx1, npbin)
    energy_hist = np.bincount(
        flat,
        weights=weights[inside] * kinetic_energy[inside],
        minlength=size,
    ).reshape(nx2, nx1, npbin)
    return weight_hist.astype(np.float32), energy_hist.astype(np.float32)


def _weighted_quantile_map(weight_hist: np.ndarray, quantile: float) -> np.ndarray:
    cumulative = np.cumsum(weight_hist, axis=-1)
    total = cumulative[..., -1]
    target = quantile * total[..., None]
    index = np.argmax(cumulative >= target, axis=-1)
    centers = np.sqrt(P_EDGES[:-1] * P_EDGES[1:])
    result = centers[index]
    result[total <= 0.0] = np.nan
    return result


def _positive_norm(values: list[np.ndarray]) -> colors.LogNorm:
    joined = np.concatenate([np.ravel(value) for value in values])
    positive = joined[np.isfinite(joined) & (joined > 0.0)]
    lower = max(float(np.percentile(positive, 1.0)), np.finfo(float).tiny)
    upper = max(float(np.percentile(positive, 99.7)), 1.01 * lower)
    return colors.LogNorm(vmin=lower, vmax=upper)


def _metric(
    snapshot: dict[str, object],
    reduction: dict[str, np.ndarray],
    variant: str,
) -> dict[str, object]:
    weight_hist = reduction["weight_hist"]
    energy_hist = reduction["energy_hist"]
    p_centers = np.sqrt(P_EDGES[:-1] * P_EDGES[1:])
    tail = p_centers >= TAIL_THRESHOLD
    total_weight = np.sum(weight_hist, axis=-1)
    x1_edges = np.asarray(snapshot["cr_x1_faces"], dtype=np.float64)
    x2_edges = np.asarray(snapshot["cr_x2_faces"], dtype=np.float64)
    code_x1_edges = np.asarray(snapshot["code_cr_x1_faces"], dtype=np.float64)
    code_x2_edges = np.asarray(snapshot["code_cr_x2_faces"], dtype=np.float64)
    area = np.diff(code_x2_edges)[:, None] * np.diff(code_x1_edges)[None, :]
    rho0 = float(snapshot["rho0"])
    u0 = float(snapshot["u0"])

    if variant == "tail_energy":
        values = np.sum(energy_hist[..., tail], axis=-1) / (area * rho0 * u0**2)
        return {
            "values": values.T,
            "x_faces": x2_edges,
            "xlabel": r"$x_2\,\omega_{pi}/c$",
            "label": r"$E_{\rm cr}(p>1.5p_{\rm inj})/(\rho_0u_0^2)$",
            "kind": "log",
        }
    if variant == "tail_fraction":
        tail_weight = np.sum(weight_hist[..., tail], axis=-1)
        values = np.full_like(total_weight, np.nan, dtype=np.float64)
        np.divide(tail_weight, total_weight, out=values, where=total_weight > 0.0)
        return {
            "values": values.T,
            "x_faces": x2_edges,
            "xlabel": r"$x_2\,\omega_{pi}/c$",
            "label": r"$W(p>1.5p_{\rm inj})/W_{\rm cr}$",
            "kind": "fraction",
        }
    if variant == "p90":
        return {
            "values": _weighted_quantile_map(weight_hist, 0.90).T,
            "x_faces": x2_edges,
            "xlabel": r"$x_2\,\omega_{pi}/c$",
            "label": r"$p_{90}/p_{\rm inj}$",
            "kind": "linear",
        }
    if variant == "p99":
        return {
            "values": _weighted_quantile_map(weight_hist, 0.99).T,
            "x_faces": x2_edges,
            "xlabel": r"$x_2\,\omega_{pi}/c$",
            "label": r"$p_{99}/p_{\rm inj}$",
            "kind": "linear",
        }
    if variant == "energy_density":
        values = np.sum(energy_hist, axis=-1) / (area * rho0 * u0**2)
        return {
            "values": values.T,
            "x_faces": x2_edges,
            "xlabel": r"$x_2\,\omega_{pi}/c$",
            "label": r"$E_{\rm cr}/(\rho_0u_0^2)$",
            "kind": "log",
        }

    dx1 = np.diff(code_x1_edges)
    if variant == "momentum_phase":
        dlnp = np.diff(np.log(P_EDGES))
        values = np.sum(weight_hist, axis=0) / (dx1[:, None] * dlnp[None, :])
        return {
            "values": values,
            "x_faces": P_EDGES,
            "xlabel": r"$p/p_{\rm inj}$",
            "label": r"$dN_{\rm cr}/(dx_1\,d\ln p)$",
            "kind": "log",
            "xscale": "log",
            "xlim": (2.0e-1, 1.0e1),
        }
    if variant == "sun_bai_energy_phase":
        pinj = float(snapshot["pinj"])
        light_speed = float(snapshot["light_speed"])
        momentum_edges = P_EDGES * pinj
        gamma_edges = np.sqrt(1.0 + (momentum_edges / light_speed) ** 2)
        energy_edges = 2.0 * (gamma_edges - 1.0) * light_speed**2 / u0**2
        log_energy_edges = np.log10(energy_edges)
        dlog_energy = np.diff(log_energy_edges)
        values = np.sum(energy_hist, axis=0) / (
            dx1[:, None] * dlog_energy[None, :]
        )
        return {
            "values": values,
            "x_faces": log_energy_edges,
            "xlabel": r"$\log_{10}[2\epsilon/(m u_0^2)]$",
            "label": r"$\epsilon f(\epsilon)$",
            "kind": "log",
        }
    raise RuntimeError(f"unknown CR variant: {variant}")


def _normalization(metrics: list[dict[str, object]]) -> colors.Normalize:
    values = [np.asarray(metric["values"], dtype=np.float64) for metric in metrics]
    kind = str(metrics[0]["kind"])
    if kind == "log":
        return _positive_norm(values)
    joined = np.concatenate([np.ravel(value) for value in values])
    finite = joined[np.isfinite(joined)]
    if kind == "fraction":
        upper = max(float(np.percentile(finite, 99.5)), 0.01)
        return colors.Normalize(vmin=0.0, vmax=upper)
    lower, upper = np.percentile(finite, (0.5, 99.5))
    return colors.Normalize(vmin=float(lower), vmax=float(upper))


def _render(
    snapshots: list[dict[str, object]],
    metrics: list[dict[str, object]],
    variant: str,
    output: Path,
    dpi: int,
    x1_frame: str,
) -> None:
    rho_values = np.concatenate([np.ravel(item["density"]) for item in snapshots])
    b_values = np.concatenate([np.ravel(item["bperp"]) for item in snapshots])
    positive_b = b_values[np.isfinite(b_values) & (b_values > 0.0)]
    rho_norm = colors.Normalize(
        vmin=max(0.0, float(np.percentile(rho_values, 0.3))),
        vmax=float(np.percentile(rho_values, 99.7)),
    )
    b_norm = colors.LogNorm(
        vmin=max(float(np.percentile(positive_b, 1.0)), 1.0e-5),
        vmax=float(np.percentile(positive_b, 99.7)),
    )
    cr_norm = _normalization(metrics)
    cr_cmap = (
        cmr.voltage_r
        if variant in ("momentum_phase", "sun_bai_energy_phase")
        else cmr.voltage
    )

    _style()
    figure = plt.figure(figsize=(17.5, 10.2))
    layout = figure.add_gridspec(2, 1, height_ratios=(0.025, 1.0), hspace=0.075)
    colorbar_grid = layout[0].subgridspec(1, 3, wspace=0.13)
    groups = layout[1].subgridspec(1, 3, wspace=0.13)
    mappables: list[object] | None = None
    for group_index, (snapshot, metric) in enumerate(zip(snapshots, metrics)):
        panels = groups[group_index].subgridspec(1, 3, wspace=0.0)
        axes = [figure.add_subplot(panels[0, index]) for index in range(3)]
        images = [
            axes[0].pcolormesh(
                snapshot["x2_faces"], snapshot["x1_faces"],
                np.asarray(snapshot["density"]).T, shading="flat",
                cmap=cmr.rainforest, norm=rho_norm, rasterized=True,
            ),
            axes[1].pcolormesh(
                snapshot["x2_faces"], snapshot["x1_faces"],
                np.asarray(snapshot["bperp"]).T, shading="flat",
                cmap=cmr.chroma, norm=b_norm, rasterized=True,
            ),
            axes[2].pcolormesh(
                metric["x_faces"], snapshot["cr_x1_faces"], metric["values"],
                shading="flat", cmap=cr_cmap, norm=cr_norm, rasterized=True,
            ),
        ]
        if mappables is None:
            mappables = images
        for column, axis in enumerate(axes):
            axis.axhline(
                float(snapshot["shock_coordinate"]), color="white",
                lw=0.8, ls="--", alpha=0.9,
            )
            if column < 2:
                horizontal = np.asarray(snapshot["x2_faces"])
                axis.set_xlabel(r"$x_2\,\omega_{pi}/c$", fontsize=9)
            else:
                horizontal = np.asarray(metric["x_faces"])
                axis.set_xlabel(str(metric["xlabel"]), fontsize=9)
                if metric.get("xscale") == "log":
                    axis.set_xscale("log")
            limits = metric.get("xlim") if column == 2 else None
            if limits is None:
                limits = (float(horizontal[0]), float(horizontal[-1]))
            axis.set_xlim(*limits)
            axis.set_ylim(
                float(snapshot["display_upstream"]),
                float(snapshot["display_downstream"]),
            )
            if column == 0:
                if group_index == 0:
                    ylabel = (
                        r"$(x_1-x_{\rm sh})\,\omega_{pi}/c$"
                        if x1_frame == "shock" else r"$x_1\,\omega_{pi}/c$"
                    )
                    axis.set_ylabel(ylabel, fontsize=10)
            else:
                axis.tick_params(labelleft=False)
        axes[1].text(
            0.5, 0.975, rf"$\Omega_0 t={float(snapshot['time']):.1f}$",
            transform=axes[1].transAxes, ha="center", va="top", fontsize=10,
            color="black",
            bbox={"facecolor": "white", "alpha": 0.82, "pad": 2.0,
                  "edgecolor": "none"},
        )

    if mappables is None:
        raise RuntimeError("no plot panels were created")
    labels = (r"$\rho/\rho_0$", r"$B_\perp/B_0$", str(metrics[0]["label"]))
    for index, (image, label) in enumerate(zip(mappables, labels)):
        axis = figure.add_subplot(colorbar_grid[0, index])
        bar = figure.colorbar(image, cax=axis, orientation="horizontal")
        bar.set_label(label, fontsize=9)
        bar.ax.tick_params(labelsize=7.5)
        bar.ax.xaxis.set_ticks_position("top")
        bar.ax.xaxis.set_label_position("top")
    figure.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--fractions", type=float, nargs=3, default=(0.34, 0.66, 1.0))
    parser.add_argument("--upstream-width", type=float, default=1080.0)
    parser.add_argument("--downstream-width", type=float, default=720.0)
    parser.add_argument("--cr-bin-factor", type=int, default=2)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--x1-frame", choices=("shock", "absolute"), default="shock")
    parser.add_argument("--refresh-cache", action="store_true")
    parser.add_argument("--dpi", type=int, default=350)
    return parser.parse_args()


def main() -> int:
    args = _arguments()
    output_dir = args.output_dir.resolve()
    cache_dir = args.cache_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    runs = [path.resolve(strict=True) for path in args.run_dir]
    mhd_paths = sorted(
        {path for run in runs for path in _files(run, "bin", "*.mhd_w_bcc.*.bin")}
    )
    particle_paths = sorted(
        {path for run in runs for path in _files(run, "pvtk", "*.prtcl_all.*.part.vtk")}
    )
    mhd_timeline = sorted((_binary_time(path), path) for path in mhd_paths)
    particle_timeline = sorted((_particle_time(path), path) for path in particle_paths)
    latest = particle_timeline[-1][0]
    particle_selection = [
        _nearest(particle_timeline, fraction * latest) for fraction in args.fractions
    ]
    mhd_selection = [_nearest(mhd_timeline, time) for time, _ in particle_selection]

    snapshots: list[dict[str, object]] = []
    tasks: list[tuple[object, ...]] = []
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
        pinj = _parameter(
            parameters, "problem", "ps_vinj_over_u0", math.sqrt(10.0)
        ) * u0
        shock_speed = 0.5 * (gamma_gas - 1.0) * u0
        length_scale = alpha_i * math.sqrt(rho0)
        omega0 = alpha_i * b0
        if args.x1_frame == "shock":
            origin = shock
            shock_coordinate = 0.0
            display_downstream = -args.downstream_width * length_scale
            display_upstream = args.upstream_width * length_scale
        else:
            origin = 0.0
            shock_coordinate = shock * length_scale
            display_downstream = (shock - args.downstream_width) * length_scale
            display_upstream = (shock + args.upstream_width) * length_scale
        cr_x1_edges = _coarsened_edges(crop_faces, args.cr_bin_factor)
        cr_x2_edges = _coarsened_edges(x2_faces, args.cr_bin_factor)
        tasks.append(
            (particle_path, cr_x2_edges, cr_x1_edges, light_speed, shock_speed, pinj)
        )
        cache_name = (
            f"{particle_path.name.removesuffix('.part.vtk')}"
            f".t{particle_time:g}.u{args.upstream_width:g}"
            f".d{args.downstream_width:g}.b{args.cr_bin_factor}.mom96.npz"
        )
        cache_paths.append(cache_dir / cache_name)
        snapshots.append(
            {
                "time": particle_time * omega0,
                "mhd_time": mhd_time,
                "x2_faces": x2_faces * length_scale,
                "x1_faces": (crop_faces - origin) * length_scale,
                "cr_x1_faces": (cr_x1_edges - origin) * length_scale,
                "cr_x2_faces": cr_x2_edges * length_scale,
                "code_cr_x1_faces": cr_x1_edges,
                "code_cr_x2_faces": cr_x2_edges,
                "shock_coordinate": shock_coordinate,
                "display_downstream": display_downstream,
                "display_upstream": display_upstream,
                "density": density[:, first : last + 1] / rho0,
                "bperp": np.sqrt(
                    by[:, first : last + 1] ** 2 + bz[:, first : last + 1] ** 2
                ) / b0,
                "rho0": rho0,
                "u0": u0,
                "pinj": pinj,
                "light_speed": light_speed,
            }
        )

    reductions: list[dict[str, np.ndarray] | None] = [None] * len(tasks)
    missing: list[int] = []
    for index, path in enumerate(cache_paths):
        if path.is_file() and not args.refresh_cache:
            with np.load(path) as cached:
                reductions[index] = {
                    "weight_hist": np.asarray(cached["weight_hist"]),
                    "energy_hist": np.asarray(cached["energy_hist"]),
                }
        else:
            missing.append(index)
    missing_tasks = [tasks[index] for index in missing]
    if missing_tasks:
        with ProcessPoolExecutor(
            max_workers=min(args.workers, len(missing_tasks))
        ) as executor:
            completed = list(executor.map(_particle_distribution, missing_tasks))
        for index, (weight_hist, energy_hist) in zip(missing, completed):
            reductions[index] = {
                "weight_hist": weight_hist,
                "energy_hist": energy_hist,
            }
            np.savez_compressed(
                cache_paths[index], weight_hist=weight_hist, energy_hist=energy_hist,
                p_edges=P_EDGES, particle_time=particle_selection[index][0],
            )
    if any(reduction is None for reduction in reductions):
        raise RuntimeError("particle distribution reduction is incomplete")

    written: list[Path] = []
    for variant in VARIANTS:
        metrics = [
            _metric(snapshot, reduction, variant)
            for snapshot, reduction in zip(snapshots, reductions)
            if reduction is not None
        ]
        output = output_dir / f"shock_long_full_9panel_cr_{variant}.png"
        _render(snapshots, metrics, variant, output, args.dpi, args.x1_frame)
        written.append(output)
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
