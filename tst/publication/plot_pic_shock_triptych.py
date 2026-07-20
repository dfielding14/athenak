#!/usr/bin/env python3
"""Plot density, transverse field, and CR energy-position at one shock time."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import cmasher as cmr
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import colors, ticker  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402
import numpy as np  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.publication.analyze_q011_section54_outputs import (  # noqa: E402
    compose_leaf_field,
    read_athenak_binary,
)
from tst.publication.plot_pic_shock_energy_position import (  # noqa: E402
    _reduce as reduce_energy_position,
)
from tst.publication.plot_pic_shock_nine_panel import (  # noqa: E402
    _parameter,
    _shock_front,
    _style,
)


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot one publication-style rho/B/CR shock triptych."
    )
    parser.add_argument("--mhd", type=Path, required=True)
    parser.add_argument("--particles", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--downstream-width", type=float, default=6000.0)
    parser.add_argument("--upstream-width", type=float, default=9000.0)
    parser.add_argument("--x-bins", type=int, default=1500)
    parser.add_argument("--energy-bins", type=int, default=120)
    parser.add_argument("--energy-min", type=float, default=0.3)
    parser.add_argument("--energy-max", type=float, default=1000.0)
    parser.add_argument("--hide-y-label", action="store_true")
    parser.add_argument(
        "--rotated",
        action="store_true",
        help="Stack panels vertically with shock-relative x1 horizontal.",
    )
    parser.add_argument("--dpi", type=int, default=350)
    return parser.parse_args()


def _load_or_reduce_energy(args: argparse.Namespace) -> dict[str, np.ndarray]:
    cache = args.cache.resolve()
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.is_file():
        with np.load(cache) as saved:
            return {name: np.asarray(saved[name]) for name in saved.files}

    args.birth_cuts = None
    reduced = reduce_energy_position(args)
    np.savez_compressed(cache, **reduced)
    return {name: np.asarray(value) for name, value in reduced.items()}


def _positive_lognorm(values: np.ndarray, low: float, high: float) -> colors.LogNorm:
    positive = values[np.isfinite(values) & (values > 0.0)]
    if positive.size == 0:
        raise RuntimeError("logarithmic panel has no positive finite values")
    vmin = max(float(np.percentile(positive, low)), np.finfo(float).tiny)
    vmax = float(np.percentile(positive, high))
    if vmax <= vmin:
        vmax = float(np.max(positive))
    if vmax <= vmin:
        vmax = np.nextafter(vmin, np.inf)
    return colors.LogNorm(vmin=vmin, vmax=vmax)


def main() -> int:
    args = _arguments()
    if args.downstream_width <= 0.0 or args.upstream_width <= 0.0:
        raise RuntimeError("shock-window widths must be positive")

    reduced = _load_or_reduce_energy(args)
    distribution = np.asarray(reduced["distribution"], dtype=np.float64)
    if distribution.ndim != 2:
        raise RuntimeError("triptych requires one all-particle energy distribution")
    cr_x1_faces = np.asarray(reduced["x_edges"], dtype=np.float64)
    log_energy_faces = np.asarray(
        reduced["log_energy_edges"], dtype=np.float64
    )
    particle_time = float(np.asarray(reduced["time"]))

    dataset = read_athenak_binary(args.mhd.resolve(strict=True))
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

    parameters = {
        block: dict(values) for block, values in dataset.input_parameters.items()
    }
    rho0 = _parameter(parameters, "problem", "ps_rho0", 1.0)
    b0 = _parameter(parameters, "problem", "ps_b0", 1.0)
    alpha_i = _parameter(
        parameters, "particles", "pic_background_ion_q_over_mc", 1.0
    )
    length_scale = alpha_i * np.sqrt(rho0)
    selected = np.flatnonzero(
        (x1 >= shock - args.downstream_width)
        & (x1 <= shock + args.upstream_width)
    )
    if selected.size == 0:
        raise RuntimeError("requested shock window does not intersect the mesh")
    first, last = int(selected[0]), int(selected[-1])

    rho = density[:, first : last + 1] / rho0
    bperp = np.sqrt(
        by[:, first : last + 1] ** 2 + bz[:, first : last + 1] ** 2
    ) / b0
    plot_x1 = (x1_faces[first : last + 2] - shock) * length_scale
    plot_x2 = x2_faces * length_scale
    display_downstream = -args.downstream_width * length_scale
    display_upstream = args.upstream_width * length_scale

    rho_vmin = max(0.0, float(np.percentile(rho, 0.3)))
    rho_vmax = float(np.percentile(rho, 99.7))
    rho_norm = colors.Normalize(vmin=rho_vmin, vmax=rho_vmax)
    b_norm = _positive_lognorm(bperp, 1.0, 99.7)
    cr_norm = _positive_lognorm(distribution, 1.0, 99.8)

    _style()
    x2_span = float(plot_x2[-1] - plot_x2[0])
    if x2_span <= 0.0:
        raise RuntimeError("transverse plot extent must be positive")
    box_aspect = (display_upstream - display_downstream) / x2_span
    if args.rotated:
        rho_norm = colors.Normalize(vmin=1.0, vmax=8.0)
        b_norm = colors.LogNorm(vmin=0.5, vmax=12.0)
        cr_norm = colors.LogNorm(vmin=10.0**0.5, vmax=10.0**2.5)
        energy_faces = np.power(10.0, log_energy_faces)
        centered_x2 = plot_x2 - 0.5 * (plot_x2[0] + plot_x2[-1])
        panel_aspect = 1.0 / box_aspect
        panel_width = 4.42
        figure_height = 3.0 * panel_width * panel_aspect + 0.96
        figure = plt.figure(figsize=(6.5, figure_height))
        layout = figure.add_gridspec(
            3, 1, left=0.16, right=0.84, bottom=0.16, top=0.97, hspace=0.10
        )
        axes = [figure.add_subplot(layout[row, 0]) for row in range(3)]
        images = [
            axes[0].pcolormesh(
                plot_x1,
                centered_x2,
                rho,
                shading="flat",
                cmap=cmr.rainforest,
                norm=rho_norm,
                rasterized=True,
            ),
            axes[1].pcolormesh(
                plot_x1,
                centered_x2,
                bperp,
                shading="flat",
                cmap=cmr.chroma,
                norm=b_norm,
                rasterized=True,
            ),
            axes[2].pcolormesh(
                cr_x1_faces,
                energy_faces,
                distribution.T,
                shading="flat",
                cmap="turbo",
                norm=cr_norm,
                rasterized=True,
            ),
        ]
        axes[0].set_ylabel(r"$x_2\,\omega_{pi}/c$", fontsize=9)
        axes[1].set_ylabel(r"$x_2\,\omega_{pi}/c$", fontsize=9)
        axes[2].set_ylabel(r"$E/E_{\rm sh}$", fontsize=9)
        axes[2].set_yscale("log")
        axes[2].set_xlabel(
            r"$(x_1-x_{\rm sh})\,\omega_{pi}/c$", fontsize=10
        )
        for index, axis in enumerate(axes):
            axis.set_box_aspect(panel_aspect)
            axis.yaxis.set_label_position("left")
            axis.axvline(0.0, color="white", lw=0.5, ls="--", alpha=0.95)
            axis.set_xlim(display_downstream, display_upstream)
            axis.set_xticks(np.arange(-6000.0, 8000.1, 2000.0))
            axis.xaxis.set_major_formatter(
                ticker.FuncFormatter(
                    lambda value, _: (
                        r"$0$"
                        if abs(value) < 0.5
                        else rf"${value / 1000.0:g}\times10^3$"
                    )
                )
            )
            if index < 2:
                axis.set_ylim(float(centered_x2[0]), float(centered_x2[-1]))
                axis.set_yticks((-1000.0, 0.0, 1000.0))
                axis.set_yticklabels((r"$-10^3$", r"$0$", r"$10^3$"))
                axis.tick_params(labelbottom=False)
            else:
                axis.set_ylim(
                    float(energy_faces[0]), float(energy_faces[-1])
                )
        axes[1].text(
            0.985,
            0.94,
            rf"$\Omega_0 t={particle_time:.1f}$",
            transform=axes[1].transAxes,
            ha="right",
            va="top",
            fontsize=10,
            color="white",
        )
        labels = (r"$\rho/\rho_0$", r"$B_\perp/B_0$", r"$E f(E)$")
        for index, (axis, image, label) in enumerate(zip(axes, images, labels)):
            color_axis = axis.inset_axes([1.018, 0.0, 0.024, 1.0])
            bar = figure.colorbar(image, cax=color_axis, orientation="vertical")
            bar.set_label(label, fontsize=9)
            bar.ax.tick_params(labelsize=7.5)
            bar.ax.yaxis.set_ticks_position("right")
            bar.ax.yaxis.set_label_position("right")
            if index == 0:
                bar.set_ticks((2.0, 4.0, 6.0, 8.0))

        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        tight_pad = 0.03
        target_content_width = 6.5 - 2.0 * tight_pad
        for _ in range(4):
            figure.canvas.draw()
            tight_box = figure.get_tightbbox(figure.canvas.get_renderer())
            scale = target_content_width / tight_box.width
            if abs(scale - 1.0) < 1.0e-4:
                break
            width, height = figure.get_size_inches()
            figure.set_size_inches(width * scale, height * scale, forward=True)
        figure.canvas.draw()
        tight_box = figure.get_tightbbox(figure.canvas.get_renderer())
        center_x = 0.5 * (tight_box.x0 + tight_box.x1)
        output_box = Bbox.from_extents(
            center_x - 3.25,
            tight_box.y0 - tight_pad,
            center_x + 3.25,
            tight_box.y1 + tight_pad,
        )
        figure.savefig(output, dpi=args.dpi, bbox_inches=output_box, pad_inches=0.0)
        plt.close(figure)
        print(f"time={particle_time:g}")
        print(f"cache={args.cache.resolve()}")
        print(f"figure={output}")
        return 0

    figure = plt.figure(figsize=(10.8, 3.35 * box_aspect + 1.75))
    layout = figure.add_gridspec(
        2, 3, height_ratios=(0.035, 1.0), hspace=0.08, wspace=0.035
    )
    axes = [figure.add_subplot(layout[1, column]) for column in range(3)]
    images = [
        axes[0].pcolormesh(
            plot_x2,
            plot_x1,
            rho.T,
            shading="flat",
            cmap=cmr.rainforest,
            norm=rho_norm,
            rasterized=True,
        ),
        axes[1].pcolormesh(
            plot_x2,
            plot_x1,
            bperp.T,
            shading="flat",
            cmap=cmr.chroma,
            norm=b_norm,
            rasterized=True,
        ),
        axes[2].pcolormesh(
            log_energy_faces,
            cr_x1_faces,
            distribution,
            shading="flat",
            cmap="turbo",
            norm=cr_norm,
            rasterized=True,
        ),
    ]

    for index, axis in enumerate(axes):
        axis.set_box_aspect(box_aspect)
        axis.axhline(0.0, color="white", lw=0.85, ls="--", alpha=0.95)
        axis.set_ylim(display_upstream, display_downstream)
        if index == 0 and not args.hide_y_label:
            axis.set_ylabel(r"$(x_1-x_{\rm sh})\,\omega_{pi}/c$", fontsize=10)
        if index > 0 or args.hide_y_label:
            axis.tick_params(labelleft=False)
    for axis in axes[:2]:
        axis.set_xlim(float(plot_x2[0]), float(plot_x2[-1]))
        axis.set_xlabel(r"$x_2\,\omega_{pi}/c$", fontsize=9)
    axes[2].set_xlim(float(log_energy_faces[0]), float(log_energy_faces[-1]))
    axes[2].set_xlabel(r"$\log_{10}(E/E_{\rm sh})$", fontsize=9)
    axes[0].xaxis.set_major_locator(ticker.MaxNLocator(nbins=6, prune="upper"))
    axes[1].xaxis.set_major_locator(ticker.MaxNLocator(nbins=6, prune="both"))
    axes[2].xaxis.set_major_locator(ticker.MaxNLocator(nbins=7, prune="lower"))
    axes[1].text(
        0.5,
        0.975,
        rf"$\Omega_0 t={particle_time:.1f}$",
        transform=axes[1].transAxes,
        ha="center",
        va="top",
        fontsize=10,
        bbox={
            "facecolor": "white",
            "alpha": 0.82,
            "pad": 2.0,
            "edgecolor": "none",
        },
    )

    labels = (r"$\rho/\rho_0$", r"$B_\perp/B_0$", r"$E f(E)$")
    for column, (image, label) in enumerate(zip(images, labels)):
        color_axis = figure.add_subplot(layout[0, column])
        bar = figure.colorbar(image, cax=color_axis, orientation="horizontal")
        bar.set_label(label, fontsize=9)
        bar.ax.tick_params(labelsize=7.5)
        bar.ax.xaxis.set_ticks_position("top")
        bar.ax.xaxis.set_label_position("top")

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(figure)
    print(f"time={particle_time:g}")
    print(f"cache={args.cache.resolve()}")
    print(f"figure={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
