"""Render documentation figures from the shipped projection-output examples."""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib.colors as colors
import matplotlib.pyplot as plt

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
    """Render narrow-slab AMR mean and dispersion panels."""
    mean = read_projection(mean_source, level=level)
    stddev = read_projection(stddev_source, level=level)
    extent = [
        float(mean["metadata"]["xmin"]),
        float(mean["metadata"]["xmax"]),
        float(mean["metadata"]["ymin"]),
        float(mean["metadata"]["ymax"]),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.3), constrained_layout=True)
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
    fig.suptitle(r"AMR blast slab: $0\leq x_1\leq0.1$")
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
        level=1,
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
