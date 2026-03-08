#!/usr/bin/env python3
"""Post-process the 1024^2 scalar_mixing stream-function production run."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import cmasher as cmr
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

from bin_convert_new import read_binary_as_athdf


THETA_MEAN = 0.5
THETA_LINTHRESH = 1.0e-2
SPECTRUM_KMAX = 512


def _read_snapshot(path: Path) -> dict[str, np.ndarray | float]:
    data = read_binary_as_athdf(
        str(path),
        quantities=["velx", "vely", "velz", "s_00"],
        dtype=np.float64,
    )
    return {
        "time": float(data["Time"]),
        "x1v": np.asarray(data["x1v"], dtype=np.float64),
        "x2v": np.asarray(data["x2v"], dtype=np.float64),
        "velx": np.asarray(data["velx"], dtype=np.float64),
        "vely": np.asarray(data["vely"], dtype=np.float64),
        "velz": np.asarray(data["velz"], dtype=np.float64),
        "scalar": np.asarray(data["s_00"], dtype=np.float64),
    }


def _spectrum_shells_2d(nx: int, ny: int) -> np.ndarray:
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    ky2, kx2 = np.meshgrid(ky, kx, indexing="ij")
    return np.ceil(np.sqrt(kx2 * kx2 + ky2 * ky2) - 1.0e-12).astype(np.int64)


def _radial_scalar_spectrum(field2d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    fluctuating = field2d - THETA_MEAN
    shat = np.fft.fftn(fluctuating)
    energy = np.abs(shat) ** 2 / field2d.size**2
    shell = _spectrum_shells_2d(field2d.shape[1], field2d.shape[0])
    spectrum = np.bincount(shell.ravel(), weights=energy.ravel())
    kvals = np.arange(1, len(spectrum), dtype=np.float64)
    return kvals, spectrum[1:]


def _radial_velocity_spectrum(vx2d: np.ndarray, vy2d: np.ndarray,
                              vz2d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    uhat_x = np.fft.fftn(vx2d)
    uhat_y = np.fft.fftn(vy2d)
    uhat_z = np.fft.fftn(vz2d)
    energy = (np.abs(uhat_x) ** 2 + np.abs(uhat_y) ** 2 + np.abs(uhat_z) ** 2) / vx2d.size**2
    shell = _spectrum_shells_2d(vx2d.shape[1], vx2d.shape[0])
    spectrum = np.bincount(shell.ravel(), weights=energy.ravel())
    kvals = np.arange(1, len(spectrum), dtype=np.float64)
    return kvals, spectrum[1:]


def _log_rebin_spectrum(kvals: np.ndarray, spectrum: np.ndarray,
                        kmax: int = SPECTRUM_KMAX,
                        nbins: int = 44) -> tuple[np.ndarray, np.ndarray]:
    edge_values = np.unique(np.rint(np.geomspace(1.0, kmax + 1.0, nbins + 1)).astype(int))
    if edge_values[0] != 1:
        edge_values = np.insert(edge_values, 0, 1)
    if edge_values[-1] != kmax + 1:
        edge_values = np.append(edge_values, kmax + 1)

    centers = []
    rebinned = []
    for lo, hi in zip(edge_values[:-1], edge_values[1:]):
        mask = (kvals >= lo) & (kvals < hi)
        if not np.any(mask):
            continue
        centers.append(np.sqrt(lo * (hi - 1)))
        rebinned.append(np.sum(spectrum[mask]))

    return np.asarray(centers, dtype=np.float64), np.asarray(rebinned, dtype=np.float64)


def _make_spectra_figure(
    snapshots: list[dict[str, np.ndarray | float]],
    output: Path,
) -> None:
    velocity = snapshots[0]
    vx2d = np.asarray(velocity["velx"])[0]
    vy2d = np.asarray(velocity["vely"])[0]
    vz2d = np.asarray(velocity["velz"])[0]
    kvals_v, spec_v = _radial_velocity_spectrum(vx2d, vy2d, vz2d)
    kvals_v, spec_v = _log_rebin_spectrum(kvals_v, spec_v)

    fig, ax = plt.subplots(figsize=(9, 6))
    scalar_cmap = plt.get_cmap("viridis")
    times = [float(snapshot["time"]) for snapshot in snapshots]
    norm = plt.Normalize(vmin=min(times), vmax=max(times))

    for snapshot in snapshots:
        scalar2d = np.asarray(snapshot["scalar"])[0]
        kvals_s, spec_s = _radial_scalar_spectrum(scalar2d)
        kvals_s, spec_s = _log_rebin_spectrum(kvals_s, spec_s)
        color = scalar_cmap(norm(float(snapshot["time"])))
        scalar_floor = np.max(spec_s) * 1.0e-12
        mask = spec_s > scalar_floor
        ax.loglog(kvals_s[mask], spec_s[mask], color=color, linewidth=1.8, alpha=0.95)

    velocity_floor = np.max(spec_v) * 1.0e-12
    mask_v = spec_v > velocity_floor
    ax.loglog(kvals_v[mask_v], spec_v[mask_v], color="black", linewidth=2.6, linestyle="--",
              label="velocity spectrum")
    ax.set_xlabel("k")
    ax.set_ylabel(
        r"$E(k)=\sum_{k_-<|\mathbf{k}'|\leq k_+} P(\mathbf{k}')$" "\n"
        r"$P_\theta=|\hat{\theta}'|^2/(N_xN_y)^2,\ "
        r"P_v=(|\hat{v}_x|^2+|\hat{v}_y|^2+|\hat{v}_z|^2)/(N_xN_y)^2$"
    )
    ax.set_title("Scalar Power Spectra and Frozen Velocity Spectrum")
    ax.set_xlim(1.0, SPECTRUM_KMAX)
    ax.set_ylim(1.0e-6, 1.0)
    ax.grid(True, which="both", alpha=0.22)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=scalar_cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("scalar snapshot time")
    ax.legend(loc="lower left", frameon=False)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _make_scalar_map_figure(
    snapshots: list[dict[str, np.ndarray | float]],
    output: Path,
) -> None:
    nshots = len(snapshots)
    ncols = 4
    nrows = math.ceil(nshots / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.1 * ncols, 3.8 * nrows),
                             squeeze=False)

    theta_absmax = max(
        float(np.max(np.abs(np.asarray(snapshot["scalar"])[0] - THETA_MEAN)))
        for snapshot in snapshots
    )
    norm = colors.SymLogNorm(
        linthresh=THETA_LINTHRESH,
        vmin=-theta_absmax,
        vmax=theta_absmax,
        base=10.0,
    )

    image = None
    for ax, snapshot in zip(axes.ravel(), snapshots):
        scalar2d = np.asarray(snapshot["scalar"])[0] - THETA_MEAN
        x = np.asarray(snapshot["x1v"])
        y = np.asarray(snapshot["x2v"])
        image = ax.imshow(
            scalar2d,
            origin="lower",
            cmap=cmr.redshift,
            extent=[x[0], x[-1], y[0], y[-1]],
            norm=norm,
            aspect="equal",
        )
        ax.set_title(f"t = {float(snapshot['time']):.1f}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    for ax in axes.ravel()[nshots:]:
        ax.axis("off")

    if image is not None:
        cbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
        cbar.set_label(r"$\theta - \langle \theta \rangle$")

    fig.suptitle(r"Scalar Field Evolution: $\theta - \langle \theta \rangle$", fontsize=16)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("build-mpi/runs/scalar_mix_stream_1024_step/bin"),
        help="Directory containing scalar_mix_stream_1024_step hydro_w binary outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build-mpi/runs/scalar_mix_stream_1024_step/analysis"),
        help="Directory for analysis plots.",
    )
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob("scalar_mix_stream_1024_step.hydro_w.*.bin"))
    if not paths:
        raise FileNotFoundError(f"No hydro_w files found in {args.input_dir}")

    snapshots = [_read_snapshot(path) for path in paths]
    _make_spectra_figure(snapshots, args.output_dir / "scalar_velocity_spectra.png")
    _make_scalar_map_figure(snapshots, args.output_dir / "scalar_maps_redshift.png")


if __name__ == "__main__":
    main()
