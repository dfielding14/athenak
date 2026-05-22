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


def _fft_wavenumbers(shape: tuple[int, int, int]) -> tuple[np.ndarray, ...]:
    nz, ny, nx = shape
    kx = np.fft.fftfreq(nx, d=1.0 / nx)
    ky = np.fft.fftfreq(ny, d=1.0 / ny)
    kz = np.fft.fftfreq(nz, d=1.0 / nz)
    kzz, kyy, kxx = np.meshgrid(kz, ky, kx, indexing="ij")
    return kxx, kyy, kzz


def vector_decomposition(vector: np.ndarray) -> dict[str, np.ndarray | float]:
    """Split a vector field into spectral solenoidal and compressive parts."""
    vector = vector - np.mean(vector, axis=(1, 2, 3), keepdims=True)
    kxx, kyy, kzz = _fft_wavenumbers(vector.shape[1:])
    kvec = (kxx, kyy, kzz)
    k2 = kxx**2 + kyy**2 + kzz**2
    mask = k2 > 1.0e-24

    vhat = np.fft.fftn(vector, axes=(1, 2, 3))
    dot = sum(kvec[comp] * vhat[comp] for comp in range(3))
    compressive = np.zeros_like(vhat)
    for comp in range(3):
        compressive[comp][mask] = kvec[comp][mask] * dot[mask] / k2[mask]
    solenoidal = vhat - compressive

    sol_mode_power = np.sum(np.abs(solenoidal) ** 2, axis=0)
    comp_mode_power = np.sum(np.abs(compressive) ** 2, axis=0)
    total_mode_power = np.sum(np.abs(vhat) ** 2, axis=0)
    kmag = np.sqrt(k2)
    shell_index = np.floor(kmag).astype(int)

    sol_shell = np.bincount(shell_index.ravel(), weights=sol_mode_power.ravel())
    comp_shell = np.bincount(shell_index.ravel(), weights=comp_mode_power.ravel())
    total_shell = np.bincount(shell_index.ravel(), weights=total_mode_power.ravel())
    shell_count = np.bincount(shell_index.ravel())

    return {
        "k": np.arange(total_shell.size),
        "solenoidal_power": sol_shell,
        "compressive_power": comp_shell,
        "total_power": total_shell,
        "counts": shell_count,
        "solenoidal_energy": float(sol_mode_power[mask].sum()),
        "compressive_energy": float(comp_mode_power[mask].sum()),
        "total_energy": float(total_mode_power[mask].sum()),
    }


def _vector_metrics_from_decomposition(
    decomp: dict[str, np.ndarray | float],
) -> dict[str, float]:
    sol_energy = float(decomp["solenoidal_energy"])
    comp_energy = float(decomp["compressive_energy"])
    total_energy = float(decomp["total_energy"])
    if total_energy <= 0.0:
        raise ValueError("Cannot compute vector decomposition metrics for a zero field.")
    sol_amp = np.sqrt(sol_energy)
    comp_amp = np.sqrt(comp_energy)
    amp_sum = sol_amp + comp_amp
    if amp_sum <= 0.0:
        raise ValueError("Cannot compute vector decomposition metrics for a zero field.")
    return {
        "solenoidal_energy_fraction": sol_energy / total_energy,
        "compressive_energy_fraction": comp_energy / total_energy,
        "amplitude_solenoidal_fraction": sol_amp / amp_sum,
        "amplitude_compressive_fraction": comp_amp / amp_sum,
    }


def vector_decomposition_metrics(vector: np.ndarray) -> dict[str, float]:
    """Return scalar solenoidal/compressive fractions for a vector field."""
    return _vector_metrics_from_decomposition(vector_decomposition(vector))


def vector_rms(vector: np.ndarray) -> float:
    """Return RMS vector magnitude after removing component means."""
    centered = vector - np.mean(vector, axis=(1, 2, 3), keepdims=True)
    return float(np.sqrt(np.mean(np.sum(centered**2, axis=0))))


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


def _plot_vector_slice(
    magnitude: np.ndarray,
    output_path: Path,
    title: str,
    colorbar_label: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mid = magnitude.shape[0] // 2
    fig, ax = plt.subplots(figsize=(6.0, 5.0), constrained_layout=True)
    image = ax.imshow(
        magnitude[mid],
        origin="lower",
        extent=(-0.5, 0.5, -0.5, 0.5),
        cmap="viridis",
        interpolation="nearest",
    )
    ax.set_xlabel("x1")
    ax.set_ylabel("x2")
    ax.set_title(title)
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label(colorbar_label)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def make_velocity_plots(
    fields: dict[str, np.ndarray],
    output_dir: str | Path,
    basename: str,
    expected_f_solenoidal: float,
) -> dict[str, float | Path]:
    """Create velocity magnitude and Helmholtz decomposition plots."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    velocity = np.stack([fields["velx"], fields["vely"], fields["velz"]])
    decomp = vector_decomposition(velocity)
    metrics = _vector_metrics_from_decomposition(decomp)
    rms = vector_rms(velocity)

    slice_path = output_dir / f"{basename}_velocity_slice.png"
    spectrum_path = output_dir / f"{basename}_velocity_decomposition.png"
    speed = np.sqrt(np.sum(velocity**2, axis=0))
    title = (
        "2D initial velocity perturbation\n"
        f"RMS(|delta v|) = {rms:.6e}"
    )
    _plot_vector_slice(speed, slice_path, title, "|delta v|")

    k = np.asarray(decomp["k"])
    sol = np.asarray(decomp["solenoidal_power"])
    comp = np.asarray(decomp["compressive_power"])
    total = np.asarray(decomp["total_power"])
    norm = total[1:].sum()
    if norm > 0.0:
        sol = sol / norm
        comp = comp / norm

    fig, ax = plt.subplots(figsize=(6.0, 4.4), constrained_layout=True)
    ax.semilogy(k, np.maximum(sol, 1.0e-18), marker="o", label="solenoidal")
    ax.semilogy(k, np.maximum(comp, 1.0e-18), marker="s", label="compressive")
    ax.set_xlabel("integer Fourier shell k")
    ax.set_ylabel("fractional shell power")
    ax.set_title(
        "Velocity Helmholtz decomposition\n"
        f"amplitude solenoidal fraction = "
        f"{metrics['amplitude_solenoidal_fraction']:.4f} "
        f"(requested {expected_f_solenoidal:.2f})"
    )
    ax.set_xlim(0, min(k[-1], 16))
    ax.set_ylim(1.0e-18, 2.0)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper right")
    fig.savefig(spectrum_path, dpi=180)
    plt.close(fig)

    return {
        "velocity_rms": rms,
        "velocity_slice": slice_path,
        "velocity_decomposition": spectrum_path,
        **metrics,
    }


def make_magnetic_plots(
    fields: dict[str, np.ndarray],
    divb_fields: dict[str, np.ndarray],
    output_dir: str | Path,
    basename: str,
) -> dict[str, float | Path]:
    """Create magnetic perturbation and divB diagnostic plots."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    magnetic = np.stack([fields["bcc1"], fields["bcc2"], fields["bcc3"]])
    magnetic = magnetic - np.mean(magnetic, axis=(1, 2, 3), keepdims=True)
    brms = vector_rms(magnetic)
    bmag = np.sqrt(np.sum(magnetic**2, axis=0))
    divb = divb_fields["divb"]
    divb_max = float(np.max(np.abs(divb)))
    divb_rms = float(np.sqrt(np.mean(divb**2)))

    slice_path = output_dir / f"{basename}_magnetic_slice.png"
    divb_path = output_dir / f"{basename}_magnetic_divb.png"
    title = (
        "2D initial magnetic perturbation\n"
        f"RMS(|delta B|) = {brms:.6e}"
    )
    _plot_vector_slice(bmag, slice_path, title, "|delta B|")

    mid = divb.shape[0] // 2
    vmax = max(float(np.max(np.abs(divb[mid]))), 1.0e-30)
    fig, ax = plt.subplots(figsize=(6.0, 5.0), constrained_layout=True)
    image = ax.imshow(
        divb[mid],
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
        "Discrete div(B) from AthenaK mhd_divb\n"
        f"max |divB| = {divb_max:.3e}, RMS = {divb_rms:.3e}"
    )
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("divB")
    fig.savefig(divb_path, dpi=180)
    plt.close(fig)

    return {
        "magnetic_rms": brms,
        "divb_max_abs": divb_max,
        "divb_rms": divb_rms,
        "magnetic_slice": slice_path,
        "magnetic_divb": divb_path,
    }


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
    parser.add_argument(
        "--fields",
        default="density",
        help="Comma-separated diagnostics to plot: density, velocity, magnetic",
    )
    parser.add_argument(
        "--divb-vtk",
        type=Path,
        default=None,
        help="Path to mhd_divb VTK output, required for magnetic divB plots",
    )
    parser.add_argument(
        "--f-solenoidal",
        type=float,
        default=0.75,
        help="Requested amplitude-space solenoidal fraction for velocity plots",
    )
    parser.add_argument("--nlow", type=int, default=1, help="Lowest requested mode")
    parser.add_argument("--nhigh", type=int, default=4, help="Highest requested mode")
    args = parser.parse_args()

    data = read_vtk_scalars(args.vtk)
    fields = data["scalars"]
    if not isinstance(fields, dict):
        raise TypeError("VTK parser returned malformed scalar data")
    requested = {
        item.strip().lower()
        for item in args.fields.replace(";", ",").replace(" ", ",").split(",")
        if item.strip()
    }
    unknown = requested - {"density", "velocity", "magnetic"}
    if unknown:
        raise ValueError(f"Unknown diagnostic field(s): {', '.join(sorted(unknown))}")

    if "density" in requested:
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
            print(f"density_{key} = {value:.16e}")
        for key, value in paths.items():
            print(f"density_{key}_path = {value}")

    if "velocity" in requested:
        velocity = make_velocity_plots(
            fields,
            args.output_dir,
            basename=args.basename,
            expected_f_solenoidal=args.f_solenoidal,
        )
        for key, value in velocity.items():
            if isinstance(value, Path):
                print(f"{key}_path = {value}")
            else:
                print(f"{key} = {value:.16e}")

    if "magnetic" in requested:
        if args.divb_vtk is None:
            raise ValueError("--divb-vtk is required for magnetic diagnostics")
        divb_data = read_vtk_scalars(args.divb_vtk)
        divb_fields = divb_data["scalars"]
        if not isinstance(divb_fields, dict):
            raise TypeError("VTK parser returned malformed divB scalar data")
        magnetic = make_magnetic_plots(
            fields,
            divb_fields,
            args.output_dir,
            basename=args.basename,
        )
        for key, value in magnetic.items():
            if isinstance(value, Path):
                print(f"{key}_path = {value}")
            else:
                print(f"{key} = {value:.16e}")


if __name__ == "__main__":
    main()
