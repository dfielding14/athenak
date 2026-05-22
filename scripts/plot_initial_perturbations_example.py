#!/usr/bin/env python3
"""Plot diagnostics for the initial Fourier perturbation example."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def _cell_dims(point_dims: tuple[int, int, int], ncell: int) -> tuple[int, int, int]:
    dims = tuple(dim - 1 if dim > 1 else dim for dim in point_dims)
    if np.prod(dims) == ncell:
        return dims
    if np.prod(point_dims) == ncell:
        return point_dims
    raise ValueError(
        f"VTK dimensions {point_dims} are inconsistent with CELL_DATA {ncell}"
    )


def read_vtk_scalars(filename: str | Path) -> dict[str, object]:
    """Read scalar fields from AthenaK legacy binary VTK output."""
    filename = Path(filename)
    with filename.open("rb") as stream:
        point_dims: tuple[int, int, int] | None = None
        ncell: int | None = None
        scalars: dict[str, np.ndarray] = {}

        while True:
            line = stream.readline()
            if not line:
                break
            text = line.decode("ascii", errors="replace").strip()
            if not text:
                continue
            if text.startswith("DIMENSIONS"):
                parts = text.split()
                point_dims = (int(parts[1]), int(parts[2]), int(parts[3]))
            elif text.startswith("CELL_DATA"):
                ncell = int(text.split()[1])
            elif text.startswith("SCALARS"):
                if point_dims is None or ncell is None:
                    raise ValueError("VTK scalar block appears before grid metadata")
                label = text.split()[1]
                lookup = stream.readline().decode("ascii", errors="replace").strip()
                if lookup != "LOOKUP_TABLE default":
                    raise ValueError(f"Unexpected VTK lookup table line: {lookup}")
                raw = stream.read(4 * ncell)
                if len(raw) != 4 * ncell:
                    raise ValueError(f"Unexpected end of file while reading {label}")
                nx, ny, nz = _cell_dims(point_dims, ncell)
                scalars[label] = np.frombuffer(raw, dtype=">f4").astype(
                    np.float64
                ).reshape((nz, ny, nx))

    if point_dims is None or ncell is None or not scalars:
        raise ValueError(f"No scalar cell data found in {filename}")
    nx, ny, nz = _cell_dims(point_dims, ncell)
    return {
        "cell_dims": (nx, ny, nz),
        "scalars": scalars,
    }


def density_contrast(
    fields: dict[str, np.ndarray], background: float = 1.0
) -> np.ndarray:
    """Return fractional density contrast from VTK scalar fields."""
    if "dens" not in fields:
        raise KeyError("VTK file does not contain the dens scalar")
    return fields["dens"] / background - 1.0


def shell_power(delta: np.ndarray) -> dict[str, np.ndarray | float]:
    """Return shell-integrated Fourier power and leakage diagnostics."""
    nz, ny, nx = delta.shape
    kx = np.fft.fftfreq(nx, d=1.0 / nx)
    ky = np.fft.fftfreq(ny, d=1.0 / ny)
    kz = np.fft.fftfreq(nz, d=1.0 / nz)
    kzz, kyy, kxx = np.meshgrid(kz, ky, kx, indexing="ij")
    kmag = np.sqrt(kxx**2 + kyy**2 + kzz**2)
    power = np.abs(np.fft.fftn(delta)) ** 2 / delta.size**2
    shell_index = np.floor(kmag).astype(int)
    total = np.bincount(shell_index.ravel(), weights=power.ravel())
    counts = np.bincount(shell_index.ravel())
    return {
        "k": np.arange(total.size),
        "power": total,
        "counts": counts,
        "mode_k": kmag,
        "mode_power": power,
    }


def density_metrics(delta: np.ndarray, nlow: int, nhigh: int) -> dict[str, float]:
    """Measure the RMS, mean, and out-of-band spectral leakage."""
    spectrum = shell_power(delta)
    kmag = np.asarray(spectrum["mode_k"])
    power = np.asarray(spectrum["mode_power"])
    selected = (kmag >= nlow - 1.0e-12) & (kmag <= nhigh + 1.0e-12)
    zero_mode = kmag < 1.0e-12
    total_power = float(power.sum())
    leakage = float(power[~(selected | zero_mode)].sum() / total_power)
    zero_fraction = float(power[zero_mode].sum() / total_power)
    return {
        "mean_delta": float(delta.mean()),
        "rms_delta": float(np.sqrt(np.mean(delta**2))),
        "leakage_fraction": leakage,
        "zero_mode_fraction": zero_fraction,
    }


def make_plots(
    delta: np.ndarray,
    output_dir: str | Path,
    nlow: int,
    nhigh: int,
    basename: str = "initial_perturbations_density",
) -> dict[str, Path]:
    """Create density slice and shell-power plots."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics = density_metrics(delta, nlow=nlow, nhigh=nhigh)

    slice_path = output_dir / f"{basename}_slice.png"
    spectrum_path = output_dir / f"{basename}_power_spectrum.png"

    mid = delta.shape[0] // 2
    vmax = np.max(np.abs(delta[mid]))
    fig, ax = plt.subplots(figsize=(6.0, 5.0), constrained_layout=True)
    image = ax.imshow(
        delta[mid],
        origin="lower",
        extent=(-0.5, 0.5, -0.5, 0.5),
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_xlabel("x1")
    ax.set_ylabel("x2")
    ax.set_title(
        "Initial density perturbation slice\n"
        f"RMS(delta rho / rho) = {metrics['rms_delta']:.6f}"
    )
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("delta rho / rho")
    fig.savefig(slice_path, dpi=180)
    plt.close(fig)

    spectrum = shell_power(delta)
    k = np.asarray(spectrum["k"])
    power = np.asarray(spectrum["power"])
    positive_power = power[1:].sum()
    normalized = np.zeros_like(power)
    if positive_power > 0.0:
        normalized = power / positive_power

    fig, ax = plt.subplots(figsize=(6.0, 4.4), constrained_layout=True)
    ax.semilogy(k, np.maximum(normalized, 1.0e-18), marker="o", color="#2a6f97")
    ax.axvspan(nlow, nhigh, color="#f2c14e", alpha=0.25, label="requested shell")
    ax.set_xlabel("integer Fourier shell k")
    ax.set_ylabel("fractional shell power")
    ax.set_title(
        "Initial density perturbation spectrum\n"
        f"out-of-band leakage = {metrics['leakage_fraction']:.2e}"
    )
    ax.set_xlim(0, min(k[-1], nhigh + 8))
    ax.set_ylim(1.0e-18, 2.0)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper right")
    fig.savefig(spectrum_path, dpi=180)
    plt.close(fig)

    return {"slice": slice_path, "spectrum": spectrum_path}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate docs figures from an AthenaK initial-perturbation VTK file."
    )
    parser.add_argument("vtk", type=Path, help="Path to InitialPerturbations VTK output")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("docs/source/_static"),
        help="Directory for generated PNG files",
    )
    parser.add_argument(
        "--basename",
        default="initial_perturbations_density",
        help="Base filename for generated PNG files",
    )
    parser.add_argument("--nlow", type=int, default=1, help="Lowest requested mode")
    parser.add_argument("--nhigh", type=int, default=4, help="Highest requested mode")
    args = parser.parse_args()

    data = read_vtk_scalars(args.vtk)
    fields = data["scalars"]
    if not isinstance(fields, dict):
        raise TypeError("VTK parser returned malformed scalar data")
    delta = density_contrast(fields)
    metrics = density_metrics(delta, nlow=args.nlow, nhigh=args.nhigh)
    paths = make_plots(
        delta,
        args.output_dir,
        nlow=args.nlow,
        nhigh=args.nhigh,
        basename=args.basename,
    )

    for key, value in metrics.items():
        print(f"{key} = {value:.16e}")
    for key, value in paths.items():
        print(f"{key}_path = {value}")


if __name__ == "__main__":
    main()
