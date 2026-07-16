#!/usr/bin/env python3
"""Make a fast density and transverse-field plot from one MHD snapshot."""

from __future__ import annotations

import argparse
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

from tst.publication.analyze_q011_section54_outputs import (  # noqa: E402
    compose_leaf_field,
    read_athenak_binary,
)
from tst.publication.plot_pic_shock_nine_panel import (  # noqa: E402
    _binary_time,
    _files,
    _parameter,
    _shock_front,
    _style,
)


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Quick two-panel density and magnetic-field shock plot."
    )
    parser.add_argument("--run-dir", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--time", type=float, help="Nearest output time; default is the latest."
    )
    parser.add_argument("--upstream-width", type=float, default=1080.0)
    parser.add_argument("--downstream-width", type=float, default=720.0)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> int:
    args = _arguments()
    paths = sorted(
        {
            path
            for run in args.run_dir
            for path in _files(run.resolve(strict=True), "bin", "*.mhd_w_bcc.*.bin")
        }
    )
    if not paths:
        raise RuntimeError("no mhd_w_bcc binary outputs found")
    timeline = sorted((_binary_time(path), path) for path in paths)
    time, path = (
        timeline[-1]
        if args.time is None
        else min(timeline, key=lambda item: abs(item[0] - args.time))
    )

    dataset = read_athenak_binary(path)
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

    parameters = {
        block: dict(values) for block, values in dataset.input_parameters.items()
    }
    rho0 = _parameter(parameters, "problem", "ps_rho0", 1.0)
    b0 = _parameter(parameters, "problem", "ps_b0", 1.0)
    alpha_i = _parameter(
        parameters, "particles", "pic_background_ion_q_over_mc", 1.0
    )
    length_scale = alpha_i * np.sqrt(rho0)
    omega0 = alpha_i * b0
    rho = density[:, first : last + 1] / rho0
    bperp = np.sqrt(
        by[:, first : last + 1] ** 2 + bz[:, first : last + 1] ** 2
    ) / b0
    plot_x1 = (x1_faces[first : last + 2] - shock) * length_scale
    plot_x2 = x2_faces * length_scale

    rho_norm = colors.Normalize(
        vmin=max(0.0, float(np.percentile(rho, 0.3))),
        vmax=float(np.percentile(rho, 99.7)),
    )
    positive_b = bperp[np.isfinite(bperp) & (bperp > 0.0)]
    b_norm = colors.LogNorm(
        vmin=max(float(np.percentile(positive_b, 1.0)), 1.0e-5),
        vmax=float(np.percentile(positive_b, 99.7)),
    )

    _style()
    figure = plt.figure(figsize=(7.2, 8.2))
    layout = figure.add_gridspec(
        2, 2, height_ratios=(0.035, 1.0), hspace=0.08, wspace=0.10
    )
    axes = [figure.add_subplot(layout[1, column]) for column in range(2)]
    images = [
        axes[0].pcolormesh(
            plot_x2, plot_x1, rho.T, shading="flat", cmap=cmr.rainforest,
            norm=rho_norm, rasterized=True,
        ),
        axes[1].pcolormesh(
            plot_x2, plot_x1, bperp.T, shading="flat", cmap=cmr.chroma,
            norm=b_norm, rasterized=True,
        ),
    ]
    for index, axis in enumerate(axes):
        axis.axhline(0.0, color="white", lw=0.8, ls="--", alpha=0.9)
        axis.set_xlim(float(plot_x2[0]), float(plot_x2[-1]))
        axis.set_ylim(
            args.upstream_width * length_scale,
            -args.downstream_width * length_scale,
        )
        axis.set_xlabel(r"$x_2\,\omega_{pi}/c$", fontsize=9)
        if index == 0:
            axis.set_ylabel(r"$(x_1-x_{\rm sh})\,\omega_{pi}/c$", fontsize=10)
        else:
            axis.tick_params(labelleft=False)
    axes[1].text(
        0.5, 0.975, rf"$\Omega_0 t={time * omega0:.1f}$",
        transform=axes[1].transAxes, ha="center", va="top", fontsize=10,
        bbox={"facecolor": "white", "alpha": 0.82, "pad": 2.0,
              "edgecolor": "none"},
    )
    for column, (image, label) in enumerate(
        zip(images, (r"$\rho/\rho_0$", r"$B_\perp/B_0$"))
    ):
        axis = figure.add_subplot(layout[0, column])
        bar = figure.colorbar(image, cax=axis, orientation="horizontal")
        bar.set_label(label, fontsize=9)
        bar.ax.xaxis.set_ticks_position("top")
        bar.ax.xaxis.set_label_position("top")

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=args.dpi, bbox_inches="tight")
    plt.close(figure)
    print(f"time={time:g}")
    print(f"input={path}")
    print(f"figure={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
