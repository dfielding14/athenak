#!/usr/bin/env python3
"""Plot full-resolution and coarsened Tiegan SGS snapshots.

Example:
  python inputs/hydro/tiegan_sgs/plot_sgs_slices.py \
      --run-dir /tmp/tiegan-sgs-m025-mpi4
"""

import argparse
from pathlib import Path
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("agg")
import matplotlib.pyplot as plt  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))

from bin_convert import (  # noqa: E402
    read_binary_as_athdf,
    read_coarsened_binary_as_athdf,
)


NUMBER_PATTERN = re.compile(r"\.(\d{5})\.(?:bin|cbin)$")
FACTOR_PATTERN = re.compile(r"cbin_sgs_f\d+_(\d+)$")


def file_number(path):
    """Return the five-digit output number encoded in an AthenaK filename."""
    match = NUMBER_PATTERN.search(path.name)
    if match is None:
        raise ValueError(f"cannot determine output number from {path}")
    return int(match.group(1))


def select_snapshot(run_dir, requested_number):
    """Find one full-resolution snapshot and all matching SGS snapshots."""
    full_files = sorted((run_dir / "bin").glob("*.prim.*.bin"))
    if not full_files:
        raise FileNotFoundError(
            f"no full-resolution *.prim.*.bin files under {run_dir}"
        )

    available = {file_number(path): path for path in full_files}
    number = max(available) if requested_number == "latest" else int(requested_number)
    if number not in available:
        raise FileNotFoundError(f"full-resolution output number {number:05d} not found")

    coarsened = []
    for directory in sorted(run_dir.glob("cbin_sgs_f*_*")):
        match = FACTOR_PATTERN.match(directory.name)
        if match is None:
            continue
        factor = int(match.group(1))
        matches = [
            path for path in directory.glob("*.cbin") if file_number(path) == number
        ]
        if matches:
            coarsened.append((factor, matches[0]))
    if not coarsened:
        raise FileNotFoundError(f"no coarsened output number {number:05d} found")

    return number, available[number], sorted(coarsened)


def plane(data, name):
    """Return the single x1-x2 plane from a fully two-dimensional output."""
    values = np.asarray(data[name])
    if values.shape[0] != 1:
        raise ValueError(f"{name} is not fully two-dimensional: shape={values.shape}")
    return values[0]


def robust_limits(values, signed=False, nonnegative=False):
    """Return color limits that suppress isolated floating-point outliers."""
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError("cannot plot a field without finite values")
    if signed:
        limit = np.percentile(np.abs(finite), 99.5)
        if limit == 0.0:
            limit = 1.0
        return -limit, limit
    lower = 0.0 if nonnegative else np.percentile(finite, 0.5)
    upper = np.percentile(finite, 99.5)
    if upper <= lower:
        scale = max(abs(lower), 1.0)
        lower -= 0.5 * scale
        upper += 0.5 * scale
    return lower, upper


def draw_panel(
    figure, axis, values, extent, title, cmap, signed=False, nonnegative=False
):
    """Draw one equal-aspect image with its own colorbar."""
    vmin, vmax = robust_limits(values, signed=signed, nonnegative=nonnegative)
    image = axis.imshow(
        values,
        origin="lower",
        extent=extent,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    axis.set_title(title)
    axis.set_xlabel(r"$x_1$")
    axis.set_ylabel(r"$x_2$")
    axis.set_aspect("equal")
    figure.colorbar(image, ax=axis, shrink=0.86)


def extent_from(data):
    """Return the Cartesian image extent from assembled AthenaK data."""
    return (data["x1f"][0], data["x1f"][-1], data["x2f"][0], data["x2f"][-1])


def plot_full_resolution(filename, output_path, dpi):
    """Write the requested one-row density and velocity figure."""
    data = read_binary_as_athdf(str(filename), quantities=["dens", "velx", "vely"])
    fields = [
        ("dens", r"$\rho$", "viridis", False),
        ("velx", r"$v_x$", "RdBu_r", True),
        ("vely", r"$v_y$", "RdBu_r", True),
    ]
    figure, axes = plt.subplots(
        1, 3, figsize=(15.5, 4.7), constrained_layout=True, squeeze=False
    )
    extent = extent_from(data)
    for axis, (name, title, cmap, signed) in zip(axes[0], fields):
        draw_panel(figure, axis, plane(data, name), extent, title, cmap, signed=signed)
    ny, nx = plane(data, "dens").shape
    figure.suptitle(
        f"Full resolution {nx} x {ny}: cycle {data['NumCycles']}, "
        f"t={data['Time']:.6g}"
    )
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)


def plot_coarsened(filename, factor, output_path, dpi):
    """Write the requested two-row coarsened state and SGS-stress figure."""
    names = ["dens", "velx", "vely", "tau_xx", "tau_xy", "tau_yy"]
    data = read_coarsened_binary_as_athdf(str(filename), quantities=names)
    fields = [
        ("dens", r"$\overline{\rho}$", "viridis", False, False),
        ("velx", r"$\widetilde{v}_x$", "RdBu_r", True, False),
        ("vely", r"$\widetilde{v}_y$", "RdBu_r", True, False),
        ("tau_xx", r"$\tau_{xx}$", "magma", False, True),
        ("tau_xy", r"$\tau_{xy}$", "RdBu_r", True, False),
        ("tau_yy", r"$\tau_{yy}$", "magma", False, True),
    ]
    figure, axes = plt.subplots(
        2, 3, figsize=(15.5, 9.0), constrained_layout=True, squeeze=False
    )
    extent = extent_from(data)
    for axis, (name, title, cmap, signed, nonnegative) in zip(axes.flat, fields):
        draw_panel(
            figure,
            axis,
            plane(data, name),
            extent,
            title,
            cmap,
            signed=signed,
            nonnegative=nonnegative,
        )
    ny, nx = plane(data, "dens").shape
    figure.suptitle(
        f"Square-filter factor {factor} ({nx} x {ny}): "
        f"cycle {data['NumCycles']}, t={data['Time']:.6g}"
    )
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument(
        "--number", default="latest", help="five-digit output number or 'latest'"
    )
    parser.add_argument("--output-dir", type=Path, help="default: <run-dir>/plots")
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    output_dir = (args.output_dir or run_dir / "plots").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    number, full_file, coarsened_files = select_snapshot(run_dir, args.number)

    full_output = output_dir / f"full_resolution_{number:05d}.png"
    plot_full_resolution(full_file, full_output, args.dpi)
    print(full_output)
    for factor, filename in coarsened_files:
        output = output_dir / f"sgs_f{factor:03d}_{number:05d}.png"
        plot_coarsened(filename, factor, output, args.dpi)
        print(output)


if __name__ == "__main__":
    main()
