#!/usr/bin/env python3
"""Make quick-look TRML slice plots from the latest AthenaK binary slice dump."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, SymLogNorm


GAMMA = 5.0 / 3.0
GM1 = GAMMA - 1.0


def add_repo_reader(repo_dir: Path) -> None:
    reader_dir = repo_dir / "vis" / "python"
    sys.path.insert(0, str(reader_dir))


def latest_slice(bin_dir: Path, basename: str, direction: str) -> Path:
    paths = sorted(bin_dir.glob(f"{basename}.slice_{direction}.*.bin"))
    if not paths:
        raise FileNotFoundError(f"No {direction}-normal slice files found in {bin_dir}")
    return paths[-1]


def sanitize(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def extent_to_start(min_value: float, global_min: float, dx: float) -> int:
    return int(round((min_value - global_min) / dx))


def assemble_slice(fd: dict, direction: str) -> dict[str, np.ndarray]:
    nx1, nx2, nx3 = fd["Nx1"], fd["Nx2"], fd["Nx3"]
    dx1 = (fd["x1max"] - fd["x1min"]) / nx1
    dx2 = (fd["x2max"] - fd["x2min"]) / nx2
    dx3 = (fd["x3max"] - fd["x3min"]) / nx3
    if direction == "x":
        shape = (nx3, nx2)
    elif direction == "y":
        shape = (nx3, nx1)
    elif direction == "z":
        shape = (nx2, nx1)
    else:
        raise ValueError(direction)

    out: dict[str, np.ndarray] = {}
    for var in fd["var_names"]:
        arr = np.full(shape, np.nan, dtype=float)
        for block_id, geom in enumerate(fd["mb_geometry"]):
            block = fd["mb_data"][var][block_id]
            n3, n2, n1 = block.shape
            x1min, _, x2min, _, x3min, _ = [float(v) for v in geom]
            i0 = extent_to_start(x1min, fd["x1min"], dx1)
            j0 = extent_to_start(x2min, fd["x2min"], dx2)
            k0 = extent_to_start(x3min, fd["x3min"], dx3)
            if direction == "x":
                arr[k0 : k0 + n3, j0 : j0 + n2] = block[:, :, 0]
            elif direction == "y":
                arr[k0 : k0 + n3, i0 : i0 + n1] = block[:, 0, :]
            else:
                arr[j0 : j0 + n2, i0 : i0 + n1] = block[0, :, :]
        out[var] = arr
    return out


def derived(fields: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    dens = fields["dens"]
    eint = fields["eint"]
    velx = fields["velx"]
    vely = fields["vely"]
    velz = fields["velz"]
    pressure = GM1 * eint
    temperature = pressure / np.maximum(dens, 1.0e-300)
    speed = np.sqrt(velx**2 + vely**2 + velz**2)
    cs = np.sqrt(np.maximum(GAMMA * pressure / np.maximum(dens, 1.0e-300), 1.0e-300))
    return {
        "dens": dens,
        "temperature": temperature,
        "pressure": pressure,
        "scalar_s00": fields.get("s_00", np.full_like(dens, np.nan)),
        "velx": velx,
        "vely": vely,
        "velz": velz,
        "speed": speed,
        "mach": speed / cs,
    }


def slice_extent(fd: dict, direction: str) -> tuple[tuple[float, float, float, float], tuple[str, str]]:
    if direction == "x":
        return (fd["x2min"], fd["x2max"], fd["x3min"], fd["x3max"]), ("y", "z")
    if direction == "y":
        return (fd["x1min"], fd["x1max"], fd["x3min"], fd["x3max"]), ("x", "z")
    return (fd["x1min"], fd["x1max"], fd["x2min"], fd["x2max"]), ("x", "y")


def read_last_lightdiag(run_dir: Path) -> dict[str, float]:
    path = run_dir / "trml_lightdiag.dat"
    if not path.exists():
        return {}
    with path.open("r") as handle:
        header = handle.readline().lstrip("#").split()
    data = np.genfromtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    return {name: float(data[-1, i]) for i, name in enumerate(header)}


def plot_norm(name: str, arr: np.ndarray):
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return None, None, None
    positive = finite[finite > 0.0]
    if name in {"dens", "temperature", "pressure"} and positive.size:
        vmin = max(np.nanpercentile(positive, 0.5), positive.min(), 1.0e-14)
        vmax = max(np.nanpercentile(positive, 99.7), vmin * 1.01)
        return LogNorm(vmin=vmin, vmax=vmax), None, None
    if name in {"velx", "vely", "velz"}:
        lim = max(np.nanpercentile(np.abs(finite), 99.0), 1.0e-8)
        return SymLogNorm(linthresh=max(0.02 * lim, 1.0e-6), vmin=-lim, vmax=lim), None, None
    vmin = np.nanpercentile(finite, 0.5)
    vmax = np.nanpercentile(finite, 99.5)
    if math.isclose(vmin, vmax):
        vmin, vmax = float(np.nanmin(finite)), float(np.nanmax(finite))
    return None, vmin, vmax


def cmap_for(name: str) -> str:
    return {
        "dens": "magma",
        "temperature": "inferno",
        "pressure": "viridis",
        "scalar_s00": "cividis",
        "speed": "plasma",
        "mach": "cubehelix",
        "velx": "coolwarm",
        "vely": "coolwarm",
        "velz": "coolwarm",
    }.get(name, "viridis")


def draw_field(ax: plt.Axes, arr: np.ndarray, extent, name: str, title: str):
    norm, vmin, vmax = plot_norm(name, arr)
    im = ax.imshow(
        arr,
        origin="lower",
        extent=extent,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap_for(name),
        norm=norm,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_title(title)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.035)
    cb.ax.tick_params(labelsize=8)
    return im


def add_tracking_lines(ax: plt.Axes, labels: tuple[str, str], diag: dict[str, float]) -> None:
    if labels[1] != "z":
        return
    for key, color, style in [
        ("ft_z_ctrl", "white", "-"),
        ("ft_z_cooling", "cyan", "--"),
        ("ft_z_peak", "lime", ":"),
    ]:
        value = diag.get(key)
        if value is not None and np.isfinite(value):
            ax.axhline(value, color=color, lw=1.1, ls=style, alpha=0.9)


def add_quiver(ax: plt.Axes, fields: dict[str, np.ndarray], extent, direction: str) -> None:
    if direction == "x":
        u, v = fields["vely"], fields["velz"]
    elif direction == "y":
        u, v = fields["velx"], fields["velz"]
    else:
        u, v = fields["velx"], fields["vely"]
    ny, nx = u.shape
    step_y = max(ny // 20, 1)
    step_x = max(nx // 20, 1)
    xs = np.linspace(extent[0], extent[1], nx)
    ys = np.linspace(extent[2], extent[3], ny)
    xx, yy = np.meshgrid(xs, ys)
    ax.quiver(
        xx[::step_y, ::step_x],
        yy[::step_y, ::step_x],
        u[::step_y, ::step_x],
        v[::step_y, ::step_x],
        color="white",
        alpha=0.55,
        scale=30.0,
        width=0.002,
    )


def write_summary(out_dir: Path, records: list[dict]) -> None:
    with (out_dir / "slice_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["direction", "time", "cycle", "quantity", "min", "p01", "mean", "median", "p99", "max"])
        for record in records:
            fields = record["fields"]
            for name, arr in fields.items():
                finite = arr[np.isfinite(arr)]
                if finite.size == 0:
                    continue
                writer.writerow(
                    [
                        record["direction"],
                        record["time"],
                        record["cycle"],
                        name,
                        np.nanmin(finite),
                        np.nanpercentile(finite, 1),
                        np.nanmean(finite),
                        np.nanmedian(finite),
                        np.nanpercentile(finite, 99),
                        np.nanmax(finite),
                    ]
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--repo-dir", type=Path, default=Path.cwd())
    parser.add_argument("--basename", default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()

    add_repo_reader(args.repo_dir)
    import bin_convert_new as bc  # noqa: PLC0415

    run_dir = args.run_dir.resolve()
    bin_dir = run_dir / "bin"
    basename = args.basename
    if basename is None:
        hydro = sorted(run_dir.glob("*.hydro.hst"))
        if not hydro:
            raise FileNotFoundError("Could not infer basename from *.hydro.hst")
        basename = hydro[0].name.removesuffix(".hydro.hst")
    out_dir = (args.out_dir or (run_dir / "status_slices_latest")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    diag = read_last_lightdiag(run_dir)

    records: list[dict] = []
    for direction in ("x", "y", "z"):
        path = latest_slice(bin_dir, basename, direction)
        fd = bc.read_binary(str(path))
        fields = derived(assemble_slice(fd, direction))
        extent, labels = slice_extent(fd, direction)
        records.append({"direction": direction, "time": fd["time"], "cycle": fd.get("cycle", np.nan), "fields": fields, "extent": extent, "labels": labels})

        overview = ["dens", "temperature", "scalar_s00", "pressure", "speed", "mach"]
        fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.2), sharex=True, sharey=True)
        for ax, name in zip(axes.ravel(), overview):
            draw_field(ax, fields[name], extent, name, name)
            if name in {"dens", "temperature", "speed"}:
                add_tracking_lines(ax, labels, diag)
            if name == "speed":
                add_quiver(ax, fields, extent, direction)
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])
        fig.suptitle(f"{direction}-normal TRML status slice, t={fd['time']:.3f}, cycle={fd.get('cycle', 'n/a')}")
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(out_dir / f"overview_slice_{direction}.png", dpi=180)
        plt.close(fig)

    for name in ["dens", "temperature", "scalar_s00", "pressure", "speed", "mach", "velx", "vely", "velz"]:
        fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
        for ax, record in zip(axes, records):
            draw_field(ax, record["fields"][name], record["extent"], name, f"{record['direction']}-normal")
            add_tracking_lines(ax, record["labels"], diag)
            ax.set_xlabel(record["labels"][0])
            ax.set_ylabel(record["labels"][1])
        time = records[0]["time"]
        fig.suptitle(f"{name} across latest orthogonal slices, t={time:.3f}")
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(out_dir / f"all_directions_{sanitize(name)}.png", dpi=180)
        plt.close(fig)

    write_summary(out_dir, records)
    with (out_dir / "README.md").open("w") as handle:
        handle.write("# TRML xi=3000 Status Slices\n\n")
        handle.write(f"Run directory: `{run_dir}`\n\n")
        handle.write(f"Slice dump time: `{records[0]['time']:.6g}`\n\n")
        handle.write("White/cyan/green horizontal lines on x/y-normal panels mark frame-tracking control, cooling, and peak z positions when available.\n")
        handle.write("The speed panels include sparse in-plane velocity arrows.\n")
    print(f"Wrote status slices to {out_dir}")


if __name__ == "__main__":
    main()
