"""Render documentation figures from the shipped projection-output examples."""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from projection import read_projection


def _axis_label(axis):
    return rf"$x_{axis[-1]}$"


def render_image(source, field, output, title, colorbar_label, cmap, log_scale=False,
                 level=None):
    """Render one two-dimensional reduction using its physical coordinate extent."""
    data = read_projection(source, level=level)
    image = data["fields"][field]
    extent = [
        float(data["metadata"]["xmin"]),
        float(data["metadata"]["xmax"]),
        float(data["metadata"]["ymin"]),
        float(data["metadata"]["ymax"]),
    ]
    norm = (
        colors.LogNorm(vmin=max(image[image > 0.0].min(), image.max() * 1.0e-6),
                       vmax=image.max())
        if log_scale else None
    )

    fig, ax = plt.subplots(figsize=(5.7, 4.8), constrained_layout=True)
    rendered = ax.imshow(
        image, origin="lower", extent=extent, cmap=cmap, norm=norm,
        interpolation="nearest"
    )
    ax.set_xlabel(_axis_label(data["metadata"]["image_axis0"]))
    ax.set_ylabel(_axis_label(data["metadata"]["image_axis1"]))
    ax.set_title(title)
    fig.colorbar(rendered, ax=ax, label=colorbar_label)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def render_profile(source, field, output, level=None):
    """Render the bounded two-axis projection as a one-dimensional profile."""
    data = read_projection(source, level=level)
    fig, ax = plt.subplots(figsize=(5.8, 4.0), constrained_layout=True)
    ax.plot(data["x"], data["fields"][field][0], color="tab:purple", linewidth=2)
    ax.set_xlabel(_axis_label(data["metadata"]["image_axis0"]))
    ax.set_ylabel(r"$\int \rho\,dx_1dx_2$")
    ax.set_title(r"Fixed-grid blast: $0\leq x_1\leq0.1$, $|x_2|\leq0.1$")
    ax.grid(alpha=0.25)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def render_amr_statistics(mean_source, stddev_source, output, level=None):
    """Render statistics plus the native leaf-patch resolution distribution."""
    mean = read_projection(mean_source, level=level)
    stddev = read_projection(stddev_source, level=level)
    extent = [
        float(mean["metadata"]["xmin"]),
        float(mean["metadata"]["xmax"]),
        float(mean["metadata"]["ymin"]),
        float(mean["metadata"]["ymax"]),
    ]
    level_map = np.full(mean["fields"]["eint"].shape, -1, dtype=np.int32)
    dx = (extent[1] - extent[0]) / level_map.shape[1]
    dy = (extent[3] - extent[2]) / level_map.shape[0]
    for patch in mean["patches"]:
        ix0 = int(round((patch["xmin"] - extent[0]) / dx))
        ix1 = int(round((patch["xmax"] - extent[0]) / dx))
        iy0 = int(round((patch["ymin"] - extent[2]) / dy))
        iy1 = int(round((patch["ymax"] - extent[2]) / dy))
        level_map[iy0:iy1, ix0:ix1] = np.maximum(
            level_map[iy0:iy1, ix0:ix1], patch["level"]
        )
    levels = sorted({patch["level"] for patch in mean["patches"]})
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.5), constrained_layout=True)
    panels = [
        (mean["fields"]["eint"], "Mass-weighted mean",
         r"$\langle e_{\rm int}\rangle_\rho$"),
        (stddev["fields"]["eint"], "Mass-weighted std. dev.",
         r"$\sigma_\rho(e_{\rm int})$"),
    ]
    for ax, (image, title, label) in zip(axes, panels):
        rendered = ax.imshow(
            image, origin="lower", extent=extent, cmap="inferno", interpolation="nearest"
        )
        ax.set_xlabel(_axis_label(mean["metadata"]["image_axis0"]))
        ax.set_ylabel(_axis_label(mean["metadata"]["image_axis1"]))
        ax.set_title(title)
        fig.colorbar(rendered, ax=ax, label=label)
    level_cmap = plt.get_cmap("viridis", max(levels) + 1)
    level_norm = colors.BoundaryNorm(
        np.arange(-0.5, max(levels) + 1.5), level_cmap.N
    )
    rendered = axes[2].imshow(
        level_map, origin="lower", extent=extent, cmap=level_cmap,
        norm=level_norm, interpolation="nearest"
    )
    for patch in mean["patches"]:
        if patch["level"] == max(levels):
            axes[2].add_patch(
                Rectangle(
                    (patch["xmin"], patch["ymin"]),
                    patch["xmax"] - patch["xmin"],
                    patch["ymax"] - patch["ymin"],
                    fill=False,
                    linewidth=0.18,
                    edgecolor="white",
                    alpha=0.45,
                )
            )
    axes[2].set_xlabel(_axis_label(mean["metadata"]["image_axis0"]))
    axes[2].set_ylabel(_axis_label(mean["metadata"]["image_axis1"]))
    axes[2].set_title("Finest contributing native patch")
    fig.colorbar(rendered, ax=axes[2], ticks=levels, label="leaf level")
    composition_level = mean["metadata"]["projection_level"]
    fig.suptitle(
        rf"AMR blast slab: $0\leq x_1\leq0.1$; native records composed at "
        rf"level {composition_level}"
    )
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = ArgumentParser()
    parser.add_argument("--fixed-map", required=True, type=Path)
    parser.add_argument("--fixed-profile", required=True, type=Path)
    parser.add_argument("--amr-mean", required=True, type=Path)
    parser.add_argument("--amr-stddev", required=True, type=Path)
    parser.add_argument("--mhd", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    render_image(
        args.fixed_map,
        "dens",
        args.output_dir / "projection_blast_column_density.png",
        "Fixed-grid blast: column density",
        r"$\int \rho\,dx_3$",
        "magma",
    )
    render_profile(
        args.fixed_profile,
        "dens",
        args.output_dir / "projection_blast_bounded_profile.png",
    )
    render_amr_statistics(
        args.amr_mean,
        args.amr_stddev,
        args.output_dir / "projection_blast_amr_slab_statistics.png",
    )
    render_image(
        args.mhd,
        "j2",
        args.output_dir / "projection_field_loop3d_current.png",
        r"3D field loop: bounded $\int |\mathbf{J}|^2\,dx_3$",
        r"$\int |\mathbf{J}|^2\,dx_3$",
        "viridis",
        log_scale=True,
    )


if __name__ == "__main__":
    main()
