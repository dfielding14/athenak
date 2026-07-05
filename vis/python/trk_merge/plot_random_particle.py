#!/usr/bin/env python3
"""Minimal example: plot one trajectory from a merged rich-track HDF5 file."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np


# Parse Athena history labels like "# [6]=B^2" into column indices.
HEADER_RE = re.compile(r"\[(\d+)\]=([^\s]+)")


def args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("tracks_h5")
    parser.add_argument("--out", default="particle_track.png")
    parser.add_argument("--row", type=int)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--b-rms", type=float)
    parser.add_argument("--history-file", type=Path)
    parser.add_argument("--min-mass", type=float)
    parser.add_argument("--mass-log-spacing", type=float)
    parser.add_argument("--smooth", type=int, default=256)
    parser.add_argument("--tmax", type=float, default=5.0)
    return parser.parse_args()


def h5_string(value) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def read_overrides(run_dir: Path) -> dict[str, str]:
    for path in sorted(run_dir.glob("*.runtime_overrides.txt")):
        values = {}
        for line in path.read_text().splitlines():
            if "=" in line and not line.lstrip().startswith("#"):
                key, value = line.split("=", 1)
                values[key.strip()] = value.strip()
        if values:
            return values
    return {}


def species_mass(species: int, run_dir: Path, opt: argparse.Namespace) -> float:
    overrides = read_overrides(run_dir)
    min_mass = opt.min_mass
    spacing = opt.mass_log_spacing
    if min_mass is None and "particles/min_mass" in overrides:
        min_mass = float(overrides["particles/min_mass"])
    if spacing is None and "particles/mass_log_spacing" in overrides:
        spacing = float(overrides["particles/mass_log_spacing"])
    if min_mass is None or spacing is None:
        raise SystemExit("missing mass ladder: pass --min-mass and --mass-log-spacing")
    return min_mass * spacing**species


def history_columns(path: Path) -> dict[str, int]:
    columns = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        for number, name in HEADER_RE.findall(line):
            columns[name] = int(number) - 1
    return columns


def brms_from_history(path: Path, target_time: float) -> tuple[float, Path, float]:
    columns = history_columns(path)
    data = np.loadtxt(path, comments="#")
    data = np.atleast_2d(data)
    row = data[np.argmin(np.abs(data[:, 0] - target_time))]
    if "B^2" in columns:
        return float(np.sqrt(row[columns["B^2"]])), path, float(row[0])
    if {"1-ME", "2-ME", "3-ME"} <= set(columns):
        magnetic_energy = (
            row[columns["1-ME"]] + row[columns["2-ME"]] + row[columns["3-ME"]]
        )
        return float(np.sqrt(2.0 * magnetic_energy)), path, float(row[0])
    raise SystemExit(f"could not infer B_rms from {path}")


def choose_brms(
    opt: argparse.Namespace, target_time: float
) -> tuple[float, Path | None, float | None]:
    if opt.b_rms is not None:
        return opt.b_rms, None, None
    if opt.history_file is not None:
        return brms_from_history(opt.history_file, target_time)
    print(
        "WARNING: no --history-file or --b-rms supplied; using B_rms=1.0. "
        "The curvature normalization is not physically scaled.",
        file=sys.stderr,
        flush=True,
    )
    return 1.0, None, None


def rolling_mean(y: np.ndarray, width: int) -> np.ndarray:
    if width <= 1:
        return y
    width = min(width, y.size)
    left = width // 2
    right = width - 1 - left
    kernel = np.ones(width) / width
    return np.convolve(np.pad(y, (left, right), mode="edge"), kernel, mode="valid")


def main() -> None:
    opt = args()

    with h5py.File(opt.tracks_h5, "r") as handle:
        fields = handle["values"].attrs["fields"]
        fields = fields.decode() if isinstance(fields, bytes) else fields
        idx = {name: i for i, name in enumerate(fields.split(","))}

        nprtcl = handle["values"].shape[0]
        row = opt.row if opt.row is not None else np.random.default_rng(opt.seed).integers(nprtcl)
        time_abs = handle["times"][:]
        time = time_abs - time_abs[0]
        keep = time <= opt.tmax
        data = handle["values"][row, keep, :]
        particle = handle["particles"][row]
        run_dir = Path(h5_string(handle.attrs.get("source_run_dir", Path(opt.tracks_h5).parent)))
        time = time[keep]

    v = data[:, [idx["vx"], idx["vy"], idx["vz"]]]
    b = data[:, [idx["bx"], idx["by"], idx["bz"]]]
    k = data[:, [idx["k1"], idx["k2"], idx["k3"]]]

    bmag = np.linalg.norm(b, axis=1)
    v2 = np.sum(v * v, axis=1)
    vpar = np.sum(v * b, axis=1) / np.maximum(bmag, 1.0e-30)
    mu_m = np.maximum(v2 - vpar * vpar, 1.0e-30) / np.maximum(2.0 * bmag, 1.0e-30)
    brms, history_path, history_time = choose_brms(opt, float(time_abs[0]))
    mass = species_mass(int(particle["species"]), run_dir, opt)
    gyro_period_length = 2.0 * np.pi * mass / brms
    kappa_scaled = np.maximum(np.linalg.norm(k, axis=1) * gyro_period_length, 1.0e-30)

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 13,
        "mathtext.fontset": "cm",
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })

    fig, axes = plt.subplots(2, 1, figsize=(5.1, 4.8), sharex=True,
                             gridspec_kw={"hspace": 0.08})
    axes[0].semilogy(time, mu_m, color="0.55", lw=0.5, alpha=0.35)
    axes[0].semilogy(time, rolling_mean(mu_m, opt.smooth), color="k", lw=1.4)
    axes[0].set_ylabel(r"$\mu_M = v_\perp^2/2B$")

    axes[1].semilogy(time, kappa_scaled, color="#1f77b4", lw=0.5, alpha=0.28)
    axes[1].semilogy(time, rolling_mean(kappa_scaled, opt.smooth), color="#1f77b4", lw=1.3)
    axes[1].axhline(1.0, color="k", ls="--", lw=0.9)
    axes[1].set_ylabel(r"$|\mathbf{K}|\,2\pi c/\Omega$")
    axes[1].set_xlabel(r"$tc/L$")

    fig.savefig(opt.out, dpi=220, bbox_inches="tight")
    print(f"{opt.out}  row={row} species={particle['species']} track_tag={particle['track_tag']}")
    print(f"mass={mass:.8g} b_rms={brms:.8g} 2pi_c_over_Omega={gyro_period_length:.8g}")
    if history_path is not None:
        print(f"history={history_path} history_time={history_time:.8g}")


if __name__ == "__main__":
    main()
