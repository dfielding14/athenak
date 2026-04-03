#!/usr/bin/env python3
"""Post-process scalar_mixing production runs in single-case or campaign mode."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

try:
    import cmasher as cmr

    SCALAR_CMAP = cmr.redshift
except ImportError:  # pragma: no cover - fallback for lighter Python environments.
    SCALAR_CMAP = "RdBu_r"

from bin_convert_new import read_binary_as_athdf


THETA_LINTHRESH = 1.0e-2
SPECTRUM_NBINS = 32
SCALAR_MAP_LIMIT = 1.5
SCALAR_MAP_DPI = 600
SCALAR_MAP_MAX_PIXELS = 4096
THREED_SPECTRUM_SLICE_FRACTIONS = (0.25, 0.5, 0.75)
THREED_MAP_SERIES_MAX_SNAPSHOTS = 9


def _warn(message: str) -> None:
    print(f"[WARN] {message}")


def _read_snapshot(path: Path, quantities: list[str]) -> dict[str, np.ndarray | float]:
    data = read_binary_as_athdf(str(path), quantities=quantities, dtype=np.float64)
    snapshot: dict[str, np.ndarray | float] = {
        "time": float(data["Time"]),
        "x1v": np.asarray(data["x1v"], dtype=np.float64),
        "x2v": np.asarray(data["x2v"], dtype=np.float64),
    }
    if "x3v" in data:
        snapshot["x3v"] = np.asarray(data["x3v"], dtype=np.float64)

    field_names = {"velx": "velx", "vely": "vely", "velz": "velz", "s_00": "scalar"}
    for quantity in quantities:
        snapshot[field_names[quantity]] = np.asarray(data[quantity], dtype=np.float64)
    return snapshot


def _canonical_field(field: np.ndarray) -> np.ndarray:
    array = np.asarray(field, dtype=np.float64)
    if array.ndim == 3 and array.shape[0] == 1:
        return array[0]
    return array


def _grid_dim(field: np.ndarray) -> int:
    array = _canonical_field(field)
    if array.ndim == 2:
        return 2
    if array.ndim == 3:
        return 3
    raise ValueError(f"Unsupported scalar field shape {array.shape}")


def _spectrum_kmax_from_field(field: np.ndarray) -> int:
    array = _canonical_field(field)
    return max(1, max(int(length) for length in array.shape) // 2)


def _spectrum_shells_2d(nx: int, ny: int) -> np.ndarray:
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    ky2, kx2 = np.meshgrid(ky, kx, indexing="ij")
    return np.ceil(np.sqrt(kx2 * kx2 + ky2 * ky2) - 1.0e-12).astype(np.int64)


def _spectrum_shells_3d(nx: int, ny: int, nz: int) -> np.ndarray:
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    kz = np.fft.fftfreq(nz) * nz
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    return np.ceil(np.sqrt(kx3 * kx3 + ky3 * ky3 + kz3 * kz3) - 1.0e-12).astype(np.int64)


def _theta_spectrum(field: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    array = _canonical_field(field)
    theta_mean = float(np.mean(array))
    fluctuating = array - theta_mean
    shat = np.fft.fftn(fluctuating)
    power = np.abs(shat) ** 2 / array.size**2

    if array.ndim == 2:
        shell = _spectrum_shells_2d(array.shape[1], array.shape[0])
    else:
        shell = _spectrum_shells_3d(array.shape[2], array.shape[1], array.shape[0])

    spectrum = np.bincount(shell.ravel(), weights=power.ravel())
    kvals = np.arange(1, len(spectrum), dtype=np.float64)
    return kvals, spectrum[1:], theta_mean


def _velocity_spectrum_2d(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vx2d = _canonical_field(vx)
    vy2d = _canonical_field(vy)
    vz2d = _canonical_field(vz)
    uhat_x = np.fft.fftn(vx2d)
    uhat_y = np.fft.fftn(vy2d)
    uhat_z = np.fft.fftn(vz2d)
    energy = (np.abs(uhat_x) ** 2 + np.abs(uhat_y) ** 2 + np.abs(uhat_z) ** 2) / vx2d.size**2
    shell = _spectrum_shells_2d(vx2d.shape[1], vx2d.shape[0])
    spectrum = np.bincount(shell.ravel(), weights=energy.ravel())
    kvals = np.arange(1, len(spectrum), dtype=np.float64)
    return kvals, spectrum[1:]


def _slice_indices(size: int) -> list[int]:
    indices = []
    for fraction in THREED_SPECTRUM_SLICE_FRACTIONS:
        index = int(round(fraction * (size - 1)))
        index = min(size - 1, max(0, index))
        if index not in indices:
            indices.append(index)
    return indices


def _average_raw_spectra(
    spectra: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray]:
    if not spectra:
        raise ValueError("No slice spectra available for averaging.")

    common_length = min(len(spec) for _, spec in spectra)
    raw_k = np.arange(1, common_length + 1, dtype=np.float64)
    raw_spec = np.mean(
        np.vstack([np.asarray(spec[:common_length], dtype=np.float64) for _, spec in spectra]),
        axis=0,
    )
    return raw_k, raw_spec


def _slice_averaged_theta_entry(
    field: np.ndarray,
    time_value: float,
    *,
    kmax: int,
    nbins: int,
) -> dict[str, np.ndarray | float | int | str]:
    array = _canonical_field(field)
    slice_spectra = []
    for axis in range(array.ndim):
        for index in _slice_indices(array.shape[axis]):
            plane = np.take(array, indices=index, axis=axis)
            kvals, spec, _ = _theta_spectrum(plane)
            slice_spectra.append((kvals, spec))

    entry = _build_rebinned_entry(
        time_value,
        *_average_raw_spectra(slice_spectra),
        kmax=kmax,
        nbins=nbins,
    )
    entry["method"] = "3d_slice_average_2dfft"
    entry["slice_count"] = len(slice_spectra)
    return entry


def _slice_averaged_velocity_curve(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    *,
    kmax: int,
    nbins: int,
) -> tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray], int, str]:
    vx3d = _canonical_field(vx)
    vy3d = _canonical_field(vy)
    vz3d = _canonical_field(vz)

    slice_spectra = []
    for axis in range(vx3d.ndim):
        for index in _slice_indices(vx3d.shape[axis]):
            plane_x = np.take(vx3d, indices=index, axis=axis)
            plane_y = np.take(vy3d, indices=index, axis=axis)
            plane_z = np.take(vz3d, indices=index, axis=axis)
            kvals, spec = _velocity_spectrum_2d(plane_x, plane_y, plane_z)
            slice_spectra.append((kvals, spec))

    raw_curve = _average_raw_spectra(slice_spectra)
    rebinned_curve = _log_rebin_spectrum(raw_curve[0], raw_curve[1], kmax=kmax, nbins=nbins)
    return raw_curve, rebinned_curve, len(slice_spectra), "3d_slice_average_2dfft"


def _log_rebin_spectrum(
    kvals: np.ndarray,
    spectrum: np.ndarray,
    kmax: int | None = None,
    nbins: int = SPECTRUM_NBINS,
) -> tuple[np.ndarray, np.ndarray]:
    if kvals.size == 0 or spectrum.size == 0:
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)
    if kmax is None:
        kmax = int(np.max(kvals))
    else:
        kmax = min(kmax, int(np.max(kvals)))
    if kmax < 1:
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)

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


def _save_theta_spectrum_npz(
    output_path: Path,
    metadata: dict[str, object],
    time_value: float,
    grid_shape: tuple[int, ...],
    kvals_raw: np.ndarray,
    spec_raw: np.ndarray,
    kvals_rebinned: np.ndarray,
    spec_rebinned: np.ndarray,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        time=np.float64(time_value),
        k_raw=np.asarray(kvals_raw, dtype=np.float64),
        theta_spectrum_raw=np.asarray(spec_raw, dtype=np.float64),
        k_rebinned=np.asarray(kvals_rebinned, dtype=np.float64),
        theta_spectrum_rebinned=np.asarray(spec_rebinned, dtype=np.float64),
        dim=np.asarray(str(metadata.get("dim", ""))),
        grid_dim=np.int64(int(metadata.get("grid_dim", 0))),
        grid_shape=np.asarray(grid_shape, dtype=np.int64),
        method=np.asarray(str(metadata.get("method", ""))),
        turb_expo=np.float64(float(metadata.get("turb_expo", 0.0))),
        scalar_diffusivity=np.float64(float(metadata.get("scalar_diffusivity", 0.0))),
        turb_rseed=np.int64(int(metadata.get("turb_rseed", 0))),
        spectrum_method=np.asarray(str(metadata.get("spectrum_method", ""))),
        spectrum_slice_count=np.int64(int(metadata.get("spectrum_slice_count", 0))),
        case_name=np.asarray(str(metadata.get("case_name", ""))),
        analyzed_dump=np.asarray(str(metadata.get("analyzed_dump", ""))),
    )


def _save_theta_series_npz(
    output_path: Path,
    metadata: dict[str, object],
    dump_paths: list[Path],
    grid_shape: tuple[int, ...],
    spectra: list[dict[str, np.ndarray | float]],
    velocity_raw: tuple[np.ndarray, np.ndarray],
    velocity_rebinned: tuple[np.ndarray, np.ndarray],
    *,
    kmax: int,
    nbins: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        times=np.asarray([float(entry["time"]) for entry in spectra], dtype=np.float64),
        dump_paths=np.asarray([str(path) for path in dump_paths]),
        theta_k_raw=np.asarray(spectra[0]["k_raw"], dtype=np.float64),
        theta_spectrum_raw=np.vstack(
            [np.asarray(entry["spec_raw"], dtype=np.float64) for entry in spectra]
        ),
        theta_k_rebinned=np.asarray(spectra[0]["k_rebinned"], dtype=np.float64),
        theta_spectrum_rebinned=np.vstack(
            [np.asarray(entry["spec_rebinned"], dtype=np.float64) for entry in spectra]
        ),
        velocity_k_raw=np.asarray(velocity_raw[0], dtype=np.float64),
        velocity_spectrum_raw=np.asarray(velocity_raw[1], dtype=np.float64),
        velocity_k_rebinned=np.asarray(velocity_rebinned[0], dtype=np.float64),
        velocity_spectrum_rebinned=np.asarray(velocity_rebinned[1], dtype=np.float64),
        kmax=np.int64(kmax),
        nbins=np.int64(nbins),
        dim=np.asarray(str(metadata.get("dim", ""))),
        grid_dim=np.int64(int(metadata.get("grid_dim", 0))),
        grid_shape=np.asarray(grid_shape, dtype=np.int64),
        method=np.asarray(str(metadata.get("method", ""))),
        turb_expo=np.float64(float(metadata.get("turb_expo", 0.0))),
        scalar_diffusivity=np.float64(float(metadata.get("scalar_diffusivity", 0.0))),
        turb_rseed=np.int64(int(metadata.get("turb_rseed", 0))),
        spectrum_method=np.asarray(str(metadata.get("spectrum_method", ""))),
        spectrum_slice_count=np.int64(int(metadata.get("spectrum_slice_count", 0))),
        case_name=np.asarray(str(metadata.get("case_name", ""))),
        analyzed_dump=np.asarray(str(metadata.get("analyzed_dump", ""))),
    )


def _plot_theta_spectrum(
    output: Path,
    title: str,
    spectra: list[dict[str, np.ndarray | float]],
    velocity_curve: tuple[np.ndarray, np.ndarray] | None = None,
    kmax: int | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    times = [float(entry["time"]) for entry in spectra]
    tmin = min(times)
    tmax = max(times)
    if tmax == tmin:
        tmax = tmin + 1.0
    norm = plt.Normalize(vmin=tmin, vmax=tmax)
    cmap = plt.get_cmap("viridis")

    for entry in spectra:
        kvals = np.asarray(entry["k_rebinned"], dtype=np.float64)
        spec = np.asarray(entry["spec_rebinned"], dtype=np.float64)
        if kvals.size == 0 or spec.size == 0:
            continue
        floor = np.max(spec) * 1.0e-12
        mask = spec > floor
        ax.loglog(
            kvals[mask],
            spec[mask],
            color=cmap(norm(float(entry["time"]))),
            linewidth=1.8,
            alpha=0.95,
        )

    if velocity_curve is not None:
        kvals_v, spec_v = velocity_curve
        if kvals_v.size and spec_v.size:
            floor = np.max(spec_v) * 1.0e-12
            mask = spec_v > floor
            ax.loglog(
                kvals_v[mask],
                spec_v[mask],
                color="black",
                linewidth=2.4,
                linestyle="--",
                label="velocity spectrum",
            )

    ax.set_xlabel("k")
    ax.set_ylabel(r"$E_\theta(k)$")
    ax.set_title(title)
    if kmax is not None:
        ax.set_xlim(1.0, float(kmax))
    ax.grid(True, which="both", alpha=0.22)
    if velocity_curve is not None:
        ax.legend(loc="lower left", frameon=False)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("snapshot time")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_scalar_map_series(paths: list[Path], output: Path) -> None:
    nshots = len(paths)
    ncols = 3 if nshots == 9 else 4
    nrows = math.ceil(nshots / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.1 * ncols, 3.8 * nrows),
        sharex=True,
        sharey=True,
        squeeze=False,
    )

    snapshots = []
    for path in paths:
        snapshot = _read_snapshot(path, ["s_00"])
        scalar = _canonical_field(np.asarray(snapshot["scalar"]))
        centered = scalar - float(np.mean(scalar))
        stride_y = max(1, math.ceil(centered.shape[0] / SCALAR_MAP_MAX_PIXELS))
        stride_x = max(1, math.ceil(centered.shape[1] / SCALAR_MAP_MAX_PIXELS))
        snapshots.append((snapshot, centered[::stride_y, ::stride_x]))

    norm = colors.SymLogNorm(
        linthresh=THETA_LINTHRESH,
        vmin=-SCALAR_MAP_LIMIT,
        vmax=SCALAR_MAP_LIMIT,
        base=10.0,
    )

    image = None
    active_axes = axes.ravel()[:nshots]
    for ax, (snapshot, centered) in zip(active_axes, snapshots):
        x = np.asarray(snapshot["x1v"])
        y = np.asarray(snapshot["x2v"])
        image = ax.imshow(
            centered,
            origin="lower",
            cmap=SCALAR_CMAP,
            extent=[x[0], x[-1], y[0], y[-1]],
            norm=norm,
            aspect="equal",
        )
        ax.set_title(f"t = {float(snapshot['time']):.1f}")
        ax.label_outer()

    for ax in axes.ravel()[nshots:]:
        ax.axis("off")

    if image is not None:
        cbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
        cbar.set_label(r"$\theta - \langle \theta \rangle$")

    fig.supxlabel("x")
    fig.supylabel("y")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=SCALAR_MAP_DPI, bbox_inches="tight")
    plt.close(fig)


def _plot_scalar_map_latest(snapshot: dict[str, np.ndarray | float], output: Path) -> None:
    scalar = _canonical_field(np.asarray(snapshot["scalar"]))
    centered = scalar - float(np.mean(scalar))
    theta_absmax = max(THETA_LINTHRESH, float(np.max(np.abs(centered))))
    norm = colors.SymLogNorm(
        linthresh=THETA_LINTHRESH,
        vmin=-theta_absmax,
        vmax=theta_absmax,
        base=10.0,
    )

    x = np.asarray(snapshot["x1v"])
    y = np.asarray(snapshot["x2v"])
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    image = ax.imshow(
        centered,
        origin="lower",
        cmap=SCALAR_CMAP,
        extent=[x[0], x[-1], y[0], y[-1]],
        norm=norm,
        aspect="equal",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Scalar map at t = {float(snapshot['time']):.3f}")
    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label(r"$\theta - \langle \theta \rangle$")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_scalar_midplanes(snapshot: dict[str, np.ndarray | float], output: Path) -> None:
    scalar = _canonical_field(np.asarray(snapshot["scalar"]))
    centered = scalar - float(np.mean(scalar))
    theta_absmax = max(THETA_LINTHRESH, float(np.max(np.abs(centered))))
    norm = colors.SymLogNorm(
        linthresh=THETA_LINTHRESH,
        vmin=-theta_absmax,
        vmax=theta_absmax,
        base=10.0,
    )

    x = np.asarray(snapshot["x1v"])
    y = np.asarray(snapshot["x2v"])
    z = np.asarray(snapshot["x3v"])
    iz = centered.shape[0] // 2
    iy = centered.shape[1] // 2
    ix = centered.shape[2] // 2

    slices = [
        ("x-y midplane", centered[iz, :, :], [x[0], x[-1], y[0], y[-1]], "x", "y"),
        ("x-z midplane", centered[:, iy, :], [x[0], x[-1], z[0], z[-1]], "x", "z"),
        ("y-z midplane", centered[:, :, ix], [y[0], y[-1], z[0], z[-1]], "y", "z"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), squeeze=False)
    image = None
    for ax, (title, plane, extent, xlabel, ylabel) in zip(axes.ravel(), slices):
        image = ax.imshow(
            plane,
            origin="lower",
            cmap=SCALAR_CMAP,
            extent=extent,
            norm=norm,
            aspect="equal",
        )
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    if image is not None:
        cbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
        cbar.set_label(r"$\theta - \langle \theta \rangle$")

    fig.suptitle(f"Scalar midplane slices at t = {float(snapshot['time']):.3f}", fontsize=15)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _select_series_paths(paths: list[Path], max_snapshots: int) -> list[Path]:
    if len(paths) <= max_snapshots:
        return list(paths)
    indices = []
    for value in np.linspace(0, len(paths) - 1, max_snapshots):
        index = int(round(float(value)))
        index = min(len(paths) - 1, max(0, index))
        if index not in indices:
            indices.append(index)
    return [paths[index] for index in indices]


def _downsample_plane(plane: np.ndarray) -> np.ndarray:
    stride_y = max(1, math.ceil(plane.shape[0] / SCALAR_MAP_MAX_PIXELS))
    stride_x = max(1, math.ceil(plane.shape[1] / SCALAR_MAP_MAX_PIXELS))
    return plane[::stride_y, ::stride_x]


def _plot_scalar_midplanes_series(paths: list[Path], output: Path) -> None:
    selected_paths = _select_series_paths(paths, THREED_MAP_SERIES_MAX_SNAPSHOTS)
    prepared_rows = []
    theta_absmax = THETA_LINTHRESH

    for path in selected_paths:
        snapshot = _read_snapshot(path, ["s_00"])
        scalar = _canonical_field(np.asarray(snapshot["scalar"]))
        centered = scalar - float(np.mean(scalar))

        x = np.asarray(snapshot["x1v"])
        y = np.asarray(snapshot["x2v"])
        z = np.asarray(snapshot["x3v"])
        iz = centered.shape[0] // 2
        iy = centered.shape[1] // 2
        ix = centered.shape[2] // 2

        slices = [
            ("x-y midplane", _downsample_plane(centered[iz, :, :]), [x[0], x[-1], y[0], y[-1]], "x", "y"),
            ("x-z midplane", _downsample_plane(centered[:, iy, :]), [x[0], x[-1], z[0], z[-1]], "x", "z"),
            ("y-z midplane", _downsample_plane(centered[:, :, ix]), [y[0], y[-1], z[0], z[-1]], "y", "z"),
        ]

        for _, plane, _, _, _ in slices:
            theta_absmax = max(theta_absmax, float(np.max(np.abs(plane))))
        prepared_rows.append((float(snapshot["time"]), slices))

    norm = colors.SymLogNorm(
        linthresh=THETA_LINTHRESH,
        vmin=-theta_absmax,
        vmax=theta_absmax,
        base=10.0,
    )

    nrows = len(prepared_rows)
    fig, axes = plt.subplots(
        nrows,
        3,
        figsize=(17, max(4.0 * nrows, 5.2)),
        squeeze=False,
    )
    image = None
    for row_axes, (time_value, slices) in zip(axes, prepared_rows, strict=True):
        for ax, (title, plane, extent, xlabel, ylabel) in zip(row_axes, slices, strict=True):
            image = ax.imshow(
                plane,
                origin="lower",
                cmap=SCALAR_CMAP,
                extent=extent,
                norm=norm,
                aspect="equal",
            )
            ax.set_title(f"{title}\nt = {time_value:.3f}")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

    if image is not None:
        cbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
        cbar.set_label(r"$\theta - \langle \theta \rangle$")

    fig.suptitle("Scalar midplane slices across time", fontsize=15)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _case_metadata_from_args(case_name: str, dump_path: Path) -> dict[str, object]:
    return {
        "case_name": case_name,
        "dim": "",
        "method": "",
        "turb_expo": 0.0,
        "scalar_diffusivity": 0.0,
        "turb_rseed": 0,
        "analyzed_dump": str(dump_path),
    }


def _build_rebinned_entry(
    time_value: float,
    kvals: np.ndarray,
    spectrum: np.ndarray,
    kmax: int | None = None,
    nbins: int = SPECTRUM_NBINS,
) -> dict[str, np.ndarray | float]:
    k_rebinned, spec_rebinned = _log_rebin_spectrum(kvals, spectrum, kmax=kmax, nbins=nbins)
    return {
        "time": time_value,
        "k_raw": np.asarray(kvals, dtype=np.float64),
        "spec_raw": np.asarray(spectrum, dtype=np.float64),
        "k_rebinned": k_rebinned,
        "spec_rebinned": spec_rebinned,
    }


def _analyze_case(
    case_name: str,
    dump_paths: list[Path],
    analysis_dir: Path,
    metadata: dict[str, object],
    latest_only: bool,
    save_spectrum_npz: bool,
) -> dict[str, object]:
    if not dump_paths:
        _warn(f"no hydro_w files found for {case_name}")
        return {
            "case_name": case_name,
            "status": "missing_outputs",
            "dump_count": 0,
            "analyzed_dump": None,
            "analysis_dir": str(analysis_dir),
        }

    analysis_dir.mkdir(parents=True, exist_ok=True)
    selected_paths = [dump_paths[-1]] if latest_only else dump_paths
    file_tag = f"_{case_name}"

    latest_snapshot = _read_snapshot(selected_paths[-1], ["s_00"])
    scalar_field = _canonical_field(np.asarray(latest_snapshot["scalar"]))
    dim = _grid_dim(scalar_field)
    spectrum_kmax = _spectrum_kmax_from_field(scalar_field)

    meta = dict(metadata)
    meta["case_name"] = case_name
    meta["dim"] = f"{dim}d"
    meta["grid_dim"] = dim
    meta["analyzed_dump"] = str(selected_paths[-1])
    meta["spectrum_method"] = "2d_full_fft" if dim == 2 else "3d_slice_average_2dfft"
    meta["spectrum_slice_count"] = 1 if dim == 2 else len(THREED_SPECTRUM_SLICE_FRACTIONS) * dim

    if dim == 2:
        latest_k, latest_spec, _ = _theta_spectrum(scalar_field)
        latest_entry = _build_rebinned_entry(
            float(latest_snapshot["time"]),
            latest_k,
            latest_spec,
            kmax=spectrum_kmax,
            nbins=SPECTRUM_NBINS,
        )
    else:
        _warn(
            "3D spectra use averaged 2D FFTs from multiple orthogonal slices "
            "instead of a full 3D FFT."
        )
        latest_entry = _slice_averaged_theta_entry(
            scalar_field,
            float(latest_snapshot["time"]),
            kmax=spectrum_kmax,
            nbins=SPECTRUM_NBINS,
        )

    if not latest_only:
        spectra = []
        for path in selected_paths:
            snapshot = _read_snapshot(path, ["s_00"])
            if dim == 2:
                kvals, spec, _ = _theta_spectrum(snapshot["scalar"])
                spectra.append(
                    _build_rebinned_entry(
                        float(snapshot["time"]),
                        kvals,
                        spec,
                        kmax=spectrum_kmax,
                        nbins=SPECTRUM_NBINS,
                    )
                )
            else:
                spectra.append(
                    _slice_averaged_theta_entry(
                        snapshot["scalar"],
                        float(snapshot["time"]),
                        kmax=spectrum_kmax,
                        nbins=SPECTRUM_NBINS,
                    )
                )

        velocity_snapshot = _read_snapshot(selected_paths[0], ["velx", "vely", "velz"])
        if dim == 2:
            kvals_v, spec_v = _velocity_spectrum_2d(
                velocity_snapshot["velx"],
                velocity_snapshot["vely"],
                velocity_snapshot["velz"],
            )
            velocity_raw_curve = (kvals_v, spec_v)
            velocity_curve = _log_rebin_spectrum(
                kvals_v,
                spec_v,
                kmax=spectrum_kmax,
                nbins=SPECTRUM_NBINS,
            )
        else:
            velocity_raw_curve, velocity_curve, velocity_slice_count, velocity_method = (
                _slice_averaged_velocity_curve(
                    velocity_snapshot["velx"],
                    velocity_snapshot["vely"],
                    velocity_snapshot["velz"],
                    kmax=spectrum_kmax,
                    nbins=SPECTRUM_NBINS,
                )
            )
            meta["spectrum_method"] = velocity_method
            meta["spectrum_slice_count"] = velocity_slice_count

        _save_theta_series_npz(
            analysis_dir / f"scalar_velocity_spectra{file_tag}.npz",
            meta,
            selected_paths,
            tuple(int(v) for v in scalar_field.shape),
            spectra,
            velocity_raw_curve,
            velocity_curve,
            kmax=spectrum_kmax,
            nbins=SPECTRUM_NBINS,
        )
        _plot_theta_spectrum(
            analysis_dir / f"scalar_velocity_spectra{file_tag}.png",
            (
                "Scalar Power Spectra and Frozen Velocity Spectrum"
                if dim == 2
                else "Slice-Averaged Scalar Power Spectra and Frozen Velocity Spectrum"
            ),
            spectra,
            velocity_curve=velocity_curve,
            kmax=spectrum_kmax,
        )
        if dim == 2:
            _plot_scalar_map_series(
                selected_paths,
                analysis_dir / f"scalar_maps_redshift{file_tag}.png",
            )
        else:
            _plot_scalar_midplanes(
                latest_snapshot,
                analysis_dir / f"scalar_midplanes{file_tag}.png",
            )
            _plot_scalar_midplanes_series(
                selected_paths,
                analysis_dir / f"scalar_midplanes_series{file_tag}.png",
            )
    else:
        _save_theta_spectrum_npz(
            analysis_dir / f"theta_power_spectrum{file_tag}.npz",
            meta,
            float(latest_snapshot["time"]),
            tuple(int(v) for v in scalar_field.shape),
            latest_entry["k_raw"],
            latest_entry["spec_raw"],
            latest_entry["k_rebinned"],
            latest_entry["spec_rebinned"],
        )
        title = (
            f"{case_name} theta power spectrum at t = {float(latest_snapshot['time']):.3f}"
            if dim == 2
            else (
                f"{case_name} slice-averaged theta power spectrum at "
                f"t = {float(latest_snapshot['time']):.3f}"
            )
        )
        _plot_theta_spectrum(
            analysis_dir / f"theta_power_spectrum{file_tag}.png",
            title,
            [latest_entry],
            kmax=spectrum_kmax,
        )
        if dim == 2:
            _plot_scalar_map_latest(
                latest_snapshot,
                analysis_dir / f"scalar_map{file_tag}.png",
            )
        else:
            _plot_scalar_midplanes(
                latest_snapshot,
                analysis_dir / f"scalar_midplanes{file_tag}.png",
            )

    return {
        "case_name": case_name,
        "status": "ok",
        "dump_count": len(dump_paths),
        "analyzed_dump": str(selected_paths[-1]),
        "analysis_dir": str(analysis_dir),
        "dim": dim,
    }


def _load_manifest(manifest_path: Path) -> list[dict[str, str]]:
    with manifest_path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def _discover_dump_paths(case_output_dir: Path) -> list[Path]:
    return sorted((case_output_dir / "bin").glob("*.hydro_w.*.bin"))


def _analyze_campaign(
    manifest_path: Path,
    analysis_root: Path,
    latest_only: bool,
    save_spectrum_npz: bool,
) -> list[dict[str, object]]:
    results = []
    for row in _load_manifest(manifest_path):
        case_name = row["case_name"]
        case_analysis_dir = analysis_root / case_name
        dump_paths = _discover_dump_paths(Path(row["output_dir"]))
        try:
            result = _analyze_case(
                case_name=case_name,
                dump_paths=dump_paths,
                analysis_dir=case_analysis_dir,
                metadata=row,
                latest_only=latest_only,
                save_spectrum_npz=save_spectrum_npz,
            )
        except Exception as exc:  # pragma: no cover - keep batch analysis alive.
            _warn(f"analysis failed for {case_name}: {exc}")
            result = {
                "case_name": case_name,
                "status": "error",
                "error": str(exc),
                "dump_count": len(dump_paths),
                "analyzed_dump": str(dump_paths[-1]) if dump_paths else None,
                "analysis_dir": str(case_analysis_dir),
            }
        results.append(result)

    analysis_root.mkdir(parents=True, exist_ok=True)
    summary_path = analysis_root / "campaign_summary.json"
    summary_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"Wrote campaign summary to {summary_path}")
    return results


def _infer_case_name(input_dir: Path) -> str:
    return input_dir.parent.name if input_dir.name == "bin" else input_dir.name


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("build-mpi/runs/scalar_mix_stream_1024_step/bin"),
        help="Directory containing hydro_w binary outputs for a single case.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build-mpi/runs/scalar_mix_stream_1024_step/analysis"),
        help="Directory for single-case analysis products.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Campaign manifest CSV for multi-case analysis.",
    )
    parser.add_argument(
        "--analysis-root",
        type=Path,
        help="Root directory for campaign analysis products.",
    )
    parser.add_argument(
        "--latest-only",
        action="store_true",
        help="Analyze only the latest available dump per case.",
    )
    parser.add_argument(
        "--save-spectrum-npz",
        action="store_true",
        help="Deprecated: spectrum NPZ files are now saved automatically.",
    )
    args = parser.parse_args()

    if args.manifest is not None:
        analysis_root = args.analysis_root if args.analysis_root is not None else args.output_dir
        _analyze_campaign(args.manifest, analysis_root, args.latest_only, args.save_spectrum_npz)
        return

    dump_paths = sorted(args.input_dir.glob("*.hydro_w.*.bin"))
    if not dump_paths:
        raise FileNotFoundError(f"No hydro_w files found in {args.input_dir}")

    case_name = _infer_case_name(args.input_dir)
    metadata = _case_metadata_from_args(case_name, dump_paths[-1])
    _analyze_case(
        case_name=case_name,
        dump_paths=dump_paths,
        analysis_dir=args.output_dir,
        metadata=metadata,
        latest_only=args.latest_only,
        save_spectrum_npz=args.save_spectrum_npz,
    )


if __name__ == "__main__":
    main()
