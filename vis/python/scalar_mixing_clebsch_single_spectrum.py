#!/usr/bin/env python3
"""Plot a single 3D Clebsch velocity spectrum with raw and compensated panels."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bin_convert_new import read_binary_as_athdf


def _shell_spectrum(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ntot = float(vx.size)
    uhat_x = np.fft.fftn(vx)
    uhat_y = np.fft.fftn(vy)
    uhat_z = np.fft.fftn(vz)
    energy = (np.abs(uhat_x) ** 2 + np.abs(uhat_y) ** 2 + np.abs(uhat_z) ** 2) / (ntot * ntot)

    nz, ny, nx = vx.shape
    kz = np.fft.fftfreq(nz) * nz if nz > 1 else np.array([0.0])
    ky = np.fft.fftfreq(ny) * ny if ny > 1 else np.array([0.0])
    kx = np.fft.fftfreq(nx) * nx
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    shell = np.ceil(np.sqrt(kx3 * kx3 + ky3 * ky3 + kz3 * kz3) - 1.0e-12).astype(np.int64)
    shell_energy = np.bincount(shell.ravel(), weights=energy.ravel())
    kvals = np.arange(shell_energy.size, dtype=np.float64)
    return kvals, shell_energy


def _fit_slope(kvals: np.ndarray, spectrum: np.ndarray, kmin: int, kmax: int) -> float:
    mask = (kvals >= kmin) & (kvals <= kmax) & (spectrum > 0.0)
    if np.count_nonzero(mask) < 2:
        return float("nan")
    coeff = np.polyfit(np.log(kvals[mask]), np.log(spectrum[mask]), 1)
    return float(coeff[0])


def _shape_normalize(spectrum: np.ndarray, kref: int) -> np.ndarray:
    normalized = np.full_like(spectrum, np.nan, dtype=np.float64)
    if kref >= spectrum.size or spectrum[kref] <= 0.0:
        return normalized
    normalized[:] = spectrum / spectrum[kref]
    return normalized


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hydro", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--nlow", type=int, required=True)
    parser.add_argument("--nhigh", type=int, required=True)
    parser.add_argument("--slope", type=float, required=True)
    parser.add_argument("--show-max", type=int, default=None)
    parser.add_argument("--kref", type=int, default=3)
    args = parser.parse_args()

    data = read_binary_as_athdf(
        str(args.hydro.resolve()),
        quantities=["velx", "vely", "velz"],
        dtype=np.float64,
    )
    vx = np.asarray(data["velx"], dtype=np.float64)
    vy = np.asarray(data["vely"], dtype=np.float64)
    vz = np.asarray(data["velz"], dtype=np.float64)
    kvals, shell_energy = _shell_spectrum(vx, vy, vz)
    fitted_slope = -_fit_slope(kvals, shell_energy, args.nlow, args.nhigh)

    show_max = int(args.show_max) if args.show_max is not None else int(min(kvals.size - 1, 2 * args.nhigh))
    show_max = min(show_max, int(kvals.size - 1))
    kref = int(args.kref)
    shape = _shape_normalize(shell_energy, kref)
    mask = kvals > 0.0
    reference = np.full_like(kvals, np.nan, dtype=np.float64)
    reference[mask] = (kvals[mask] / float(kref)) ** (-args.slope)
    if kref < reference.size and math.isfinite(reference[kref]) and reference[kref] != 0.0:
        reference *= shell_energy[kref] / reference[kref]

    compensated = np.full_like(shape, np.nan, dtype=np.float64)
    compensated[mask] = shape[mask] * (kvals[mask] / float(kref)) ** args.slope

    finite_raw = shell_energy[(kvals >= args.nlow) & (kvals <= show_max) & (shell_energy > 0.0)]
    ylo = float(np.min(finite_raw) * 0.6)
    yhi = float(np.max(finite_raw) * 1.7)

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True, constrained_layout=True)
    top, bottom = axes

    top.loglog(kvals[1:show_max + 1], shell_energy[1:show_max + 1], color="#275dad", lw=2.4, label=r"$E_v(k)$")
    top.loglog(
        kvals[max(args.nlow, 1):show_max + 1],
        reference[max(args.nlow, 1):show_max + 1],
        color="k",
        lw=1.3,
        ls="--",
        label=rf"$k^{{-{args.slope:.3f}}}$ reference",
    )
    top.axvline(args.nlow, color="0.6", lw=1.0, ls=":")
    top.axvline(args.nhigh, color="0.6", lw=1.0, ls=":")
    top.set_ylabel(r"$E_v(k)$")
    top.set_ylim(ylo, yhi)
    top.grid(True, which="both", alpha=0.25)
    top.legend(frameon=False, loc="lower left")
    top.set_title(
        rf"3D Clebsch Velocity Spectrum ($k_{{\min}}={args.nlow}$, $k_{{\max}}={args.nhigh}$, fit$={fitted_slope:.3f}$)"
    )

    bottom.axhline(1.0, color="k", lw=1.2, ls="--", label="target-compensated flat law")
    bottom.loglog(
        kvals[max(args.nlow, 1):show_max + 1],
        compensated[max(args.nlow, 1):show_max + 1],
        color="#bd632f",
        lw=2.4,
        marker="o",
        ms=2.8,
        label=rf"$E_v(k)\,k^{{{args.slope:.3f}}}$",
    )
    bottom.axvline(args.nlow, color="0.6", lw=1.0, ls=":")
    bottom.axvline(args.nhigh, color="0.6", lw=1.0, ls=":")
    bottom.set_xlabel("k")
    bottom.set_ylabel(r"$[E_v(k)/E_v(k_{\rm ref})]\,(k/k_{\rm ref})^s$")
    bottom.set_xlim(float(args.nlow), float(show_max))
    bottom.grid(True, which="both", alpha=0.25)
    bottom.legend(frameon=False, loc="lower left")

    args.output.resolve().parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output.resolve(), dpi=220)
    plt.close(fig)

    args.summary.resolve().write_text(
        json.dumps(
            {
                "plot": str(args.output.resolve()),
                "hydro": str(args.hydro.resolve()),
                "show_max": show_max,
                "ylim": [ylo, yhi],
                "target_velocity_slope": args.slope,
                "fitted_velocity_slope": fitted_slope,
                "kmin": args.nlow,
                "kmax": args.nhigh,
                "kref": kref,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
