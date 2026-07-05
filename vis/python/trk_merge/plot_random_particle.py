#!/usr/bin/env python3
"""Minimal example: plot one trajectory from a merged rich-track HDF5 file."""

from __future__ import annotations

import argparse

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np


def args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("tracks_h5")
    parser.add_argument("--out", default="particle_track.png")
    parser.add_argument("--row", type=int)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--rg", type=float, default=1.0)
    parser.add_argument("--smooth", type=int, default=256)
    parser.add_argument("--tmax", type=float, default=5.0)
    return parser.parse_args()


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
        time = handle["times"][:] - handle["times"][0]
        keep = time <= opt.tmax
        data = handle["values"][row, keep, :]
        particle = handle["particles"][row]
        time = time[keep]

    v = data[:, [idx["vx"], idx["vy"], idx["vz"]]]
    b = data[:, [idx["bx"], idx["by"], idx["bz"]]]
    k = data[:, [idx["k1"], idx["k2"], idx["k3"]]]

    bmag = np.linalg.norm(b, axis=1)
    v2 = np.sum(v * v, axis=1)
    vpar = np.sum(v * b, axis=1) / np.maximum(bmag, 1.0e-30)
    mu_m = np.maximum(v2 - vpar * vpar, 1.0e-30) / np.maximum(2.0 * bmag, 1.0e-30)
    kappa_rg = np.maximum(np.linalg.norm(k, axis=1) * opt.rg, 1.0e-30)

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

    axes[1].semilogy(time, kappa_rg, color="#1f77b4", lw=0.5, alpha=0.28)
    axes[1].semilogy(time, rolling_mean(kappa_rg, opt.smooth), color="#1f77b4", lw=1.3)
    axes[1].axhline(1.0, color="k", ls="--", lw=0.9)
    axes[1].set_ylabel(r"$|\mathbf{K}| r_g$")
    axes[1].set_xlabel(r"$tc/L$")

    fig.savefig(opt.out, dpi=220, bbox_inches="tight")
    print(f"{opt.out}  row={row} species={particle['species']} track_tag={particle['track_tag']}")


if __name__ == "__main__":
    main()
