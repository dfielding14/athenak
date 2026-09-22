#!/usr/bin/env python3
"""Analyze the Tiegan SGS 1024^3 calibration run.

The full primitive cubes are intentionally not used for spectra.  Velocity spectra
are estimated from the three orthogonal full-resolution slice outputs and then
averaged.
"""

from __future__ import annotations

import argparse
import csv
import gc
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
import cmasher as cmr


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
from bin_convert import read_binary, read_coarsened_binary  # noqa: E402


RUN_NAME = "tiegan_sgs3d_m03333_1024_k2_dedt0022_mb256"
SGS_FACTORS = (4, 8, 16, 32)
SLICE_IDS = ("slice_x1", "slice_x2", "slice_x3")
VELOCITY_NAMES = ("velx", "vely", "velz")
SGS_NAMES = (
    "dens",
    "velx",
    "vely",
    "velz",
    "tau_xx",
    "tau_xy",
    "tau_xz",
    "tau_yy",
    "tau_yz",
    "tau_zz",
)
CMASHER_CMAPS = {
    "density": "cmr.rainforest",
    "speed": "cmr.chroma",
    "density_contrast": "cmr.fusion",
    "velocity": "cmr.wildfire",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=Path("/lustre/orion/ast207/proj-shared/dfielding/SGS/data") / RUN_NAME,
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("/lustre/orion/ast207/proj-shared/dfielding/SGS/plots") / RUN_NAME,
    )
    parser.add_argument("--force", action="store_true", help="ignore cached spectra")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--sample-points", type=int, default=400_000)
    return parser.parse_args()


def sorted_files(path: Path, pattern: str) -> list[Path]:
    files = sorted(path.glob(pattern))
    if not files:
        raise FileNotFoundError(f"no files match {path / pattern}")
    return files


def file_number(path: Path) -> int:
    match = re.search(r"\.(\d{5})\.(?:bin|cbin)$", path.name)
    if not match:
        raise ValueError(f"could not parse output number from {path}")
    return int(match.group(1))


def read_history(path: Path) -> dict[str, np.ndarray]:
    names: list[str] | None = None
    rows: list[list[float]] = []
    pattern = re.compile(r"\[(\d+)\]=([^\s]+)")
    with path.open() as stream:
        for raw in stream:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                matches = pattern.findall(line)
                if matches:
                    ordered = sorted((int(i), name) for i, name in matches)
                    names = [name for _i, name in ordered]
                continue
            rows.append([float(x) for x in line.split()])
    if names is None:
        raise ValueError(f"could not find history header in {path}")
    arr = np.asarray(rows, dtype=np.float64)
    if arr.shape[1] != len(names):
        raise ValueError(f"history column mismatch in {path}: {arr.shape[1]} vs {len(names)}")
    out = {name: arr[:, i] for i, name in enumerate(names)}
    ke = out["1-KE"] + out["2-KE"] + out["3-KE"]
    out["ke_total"] = ke
    out["mach"] = np.sqrt(2.0 * ke / out["mass"])
    out["mom_mag"] = np.sqrt(out["1-mom"] ** 2 + out["2-mom"] ** 2 + out["3-mom"] ** 2)
    return out


def plot_history(history: dict[str, np.ndarray], out_dir: Path) -> dict[str, float | None]:
    time = history["time"]
    mach = history["mach"]
    ke = history["ke_total"]
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.5), constrained_layout=True)
    axes = axes.ravel()

    axes[0].plot(time, mach, color="tab:blue", lw=1.8)
    axes[0].axhline(0.3333, color="k", ls="--", lw=1.0, alpha=0.7)
    axes[0].set_ylabel("rms Mach")
    axes[0].set_xlabel("time")
    axes[0].grid(alpha=0.25)

    axes[1].plot(time, ke, color="tab:orange", lw=1.8)
    axes[1].set_ylabel("total kinetic energy")
    axes[1].set_xlabel("time")
    axes[1].grid(alpha=0.25)

    axes[2].plot(time, history["1-KE"], label="x", lw=1.5)
    axes[2].plot(time, history["2-KE"], label="y", lw=1.5)
    axes[2].plot(time, history["3-KE"], label="z", lw=1.5)
    axes[2].set_ylabel("component KE")
    axes[2].set_xlabel("time")
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.25)

    axes[3].plot(time, history["dt"], color="tab:green", lw=1.8)
    axes[3].set_ylabel("dt")
    axes[3].set_xlabel("time")
    axes[3].grid(alpha=0.25)

    fig.savefig(out_dir / "history_mach_energy.png", dpi=180)
    plt.close(fig)

    crossed = np.flatnonzero(mach >= 0.3333)
    target_crossing = float(time[crossed[0]]) if crossed.size else None
    return {
        "final_time": float(time[-1]),
        "final_dt": float(history["dt"][-1]),
        "final_mach": float(mach[-1]),
        "final_ke": float(ke[-1]),
        "target_mach_crossing_time": target_crossing,
        "mean_mach_t_ge_2": float(np.mean(mach[time >= 2.0])),
        "std_mach_t_ge_2": float(np.std(mach[time >= 2.0])),
    }


def active_storage_axes(data: dict) -> list[int]:
    # Arrays are stored as (x3, x2, x1).  Keep only non-singleton output axes.
    out_block = [int(data["nx3_out_mb"]), int(data["nx2_out_mb"]), int(data["nx1_out_mb"])]
    return [axis for axis, size in enumerate(out_block) if size > 1]


def assemble_plane(data: dict, names: tuple[str, ...]) -> dict[str, np.ndarray]:
    axes = active_storage_axes(data)
    if len(axes) != 2:
        raise ValueError(f"expected a 2D slice, found active storage axes {axes}")

    logical_for_storage_axis = {0: 2, 1: 1, 2: 0}
    row_logical = logical_for_storage_axis[axes[0]]
    col_logical = logical_for_storage_axis[axes[1]]
    logical = np.asarray(data["mb_logical"])

    sample = np.squeeze(np.asarray(data["mb_data"][names[0]][0]))
    block_nrow, block_ncol = sample.shape
    nrow = (int(logical[:, row_logical].max()) + 1) * block_nrow
    ncol = (int(logical[:, col_logical].max()) + 1) * block_ncol
    out = {name: np.empty((nrow, ncol), dtype=np.float32) for name in names}

    for block, loc in enumerate(logical):
        row0 = int(loc[row_logical]) * block_nrow
        col0 = int(loc[col_logical]) * block_ncol
        for name in names:
            out[name][row0 : row0 + block_nrow, col0 : col0 + block_ncol] = np.squeeze(
                np.asarray(data["mb_data"][name][block])
            )
    return out


def assemble_3d(data: dict, names: tuple[str, ...]) -> dict[str, np.ndarray]:
    logical = np.asarray(data["mb_logical"])
    sample = np.asarray(data["mb_data"][names[0]][0])
    nz_block, ny_block, nx_block = sample.shape
    nx = (int(logical[:, 0].max()) + 1) * nx_block
    ny = (int(logical[:, 1].max()) + 1) * ny_block
    nz = (int(logical[:, 2].max()) + 1) * nz_block
    out = {name: np.empty((nz, ny, nx), dtype=np.float32) for name in names}

    for block, (lx1, lx2, lx3, _level) in enumerate(logical):
        i0 = int(lx1) * nx_block
        j0 = int(lx2) * ny_block
        k0 = int(lx3) * nz_block
        for name in names:
            out[name][k0 : k0 + nz_block, j0 : j0 + ny_block, i0 : i0 + nx_block] = (
                np.asarray(data["mb_data"][name][block])
            )
    return out


def block_reduce(a: np.ndarray, max_cells: int = 900) -> np.ndarray:
    if a.ndim != 2:
        raise ValueError("block_reduce expects a 2D array")
    ny, nx = a.shape
    factor = max(1, int(math.ceil(max(ny, nx) / max_cells)))
    if factor == 1:
        return a
    ny2 = (ny // factor) * factor
    nx2 = (nx // factor) * factor
    return a[:ny2, :nx2].reshape(ny2 // factor, factor, nx2 // factor, factor).mean(
        axis=(1, 3)
    )


def robust_limits(a: np.ndarray, symmetric: bool = False, pct: float = 99.5) -> tuple[float, float]:
    finite = np.asarray(a[np.isfinite(a)])
    if finite.size == 0:
        return -1.0, 1.0
    if symmetric:
        vmax = float(np.percentile(np.abs(finite), pct))
        if vmax == 0.0:
            vmax = float(np.max(np.abs(finite))) or 1.0
        return -vmax, vmax
    lo = float(np.percentile(finite, 100.0 - pct))
    hi = float(np.percentile(finite, pct))
    if lo == hi:
        hi = lo + 1.0
    return lo, hi


def add_image(ax, values: np.ndarray, title: str, *, signed: bool = False, symlog: bool = False):
    image = block_reduce(values)
    cmap = "RdBu_r" if signed else "viridis"
    vmin, vmax = robust_limits(image, symmetric=signed)
    norm = None
    if symlog:
        vmax_abs = max(abs(vmin), abs(vmax))
        norm = SymLogNorm(linthresh=max(vmax_abs * 1.0e-3, 1.0e-12), vmin=-vmax_abs, vmax=vmax_abs)
        vmin = vmax = None
    im = ax.imshow(image, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)


def positive_log_limits(a: np.ndarray, pct_lo: float = 0.5, pct_hi: float = 99.5) -> tuple[float, float]:
    finite = np.asarray(a[np.isfinite(a) & (a > 0.0)])
    if finite.size == 0:
        return 1.0e-12, 1.0
    lo = float(np.percentile(finite, pct_lo))
    hi = float(np.percentile(finite, pct_hi))
    if lo <= 0.0:
        lo = float(np.min(finite))
    if hi <= lo:
        hi = lo * 10.0
    return lo, hi


def add_latest_slice_image(
    ax,
    values: np.ndarray,
    title: str,
    *,
    cmap: str,
    signed: bool = False,
    log_scale: bool = False,
):
    image = block_reduce(values)
    norm = None
    if log_scale:
        vmin, vmax = positive_log_limits(image)
        norm = LogNorm(vmin=vmin, vmax=vmax)
        im = ax.imshow(image, origin="lower", cmap=cmap, norm=norm)
    else:
        vmin, vmax = robust_limits(image, symmetric=signed)
        im = ax.imshow(image, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)


def save_latest_slice_plots(run_dir: Path, out_dir: Path) -> dict[str, dict[str, float]]:
    summary: dict[str, dict[str, float]] = {}
    for slice_id in SLICE_IDS:
        path = sorted_files(run_dir / "bin", f"*.{slice_id}.*.bin")[-1]
        data = read_binary(str(path))
        fields = assemble_plane(data, ("dens", "velx", "vely", "velz"))
        speed = np.sqrt(fields["velx"] ** 2 + fields["vely"] ** 2 + fields["velz"] ** 2)
        density_contrast = fields["dens"] - np.mean(fields["dens"], dtype=np.float64)
        fields_to_plot = (
            ("density", fields["dens"], CMASHER_CMAPS["density"], False, False),
            ("velocity magnitude", speed, CMASHER_CMAPS["speed"], False, True),
            ("density contrast", density_contrast, CMASHER_CMAPS["density_contrast"], True, False),
            ("velx", fields["velx"], CMASHER_CMAPS["velocity"], True, False),
            ("vely", fields["vely"], CMASHER_CMAPS["velocity"], True, False),
            ("velz", fields["velz"], CMASHER_CMAPS["velocity"], True, False),
        )
        fig, axes = plt.subplots(2, 3, figsize=(12.4, 7.4), constrained_layout=True)
        for ax, (name, values, cmap, signed, log_scale) in zip(axes.ravel(), fields_to_plot):
            add_latest_slice_image(
                ax,
                values,
                f"{slice_id}: {name}",
                cmap=cmap,
                signed=signed,
                log_scale=log_scale,
            )
        fig.savefig(out_dir / f"latest_{slice_id}_fields.png", dpi=180)
        plt.close(fig)

        flat = {
            "dens": fields["dens"].ravel(),
            "velx": fields["velx"].ravel(),
            "vely": fields["vely"].ravel(),
            "velz": fields["velz"].ravel(),
            "speed": speed.ravel(),
        }
        fig, axes = plt.subplots(2, 3, figsize=(12.0, 6.8), constrained_layout=True)
        for ax, (name, values) in zip(axes.ravel(), flat.items()):
            ax.hist(values[np.isfinite(values)], bins=180, histtype="step", density=True, lw=1.4)
            ax.set_title(f"{slice_id}: {name}", fontsize=10)
            ax.grid(alpha=0.25)
        axes.ravel()[-1].axis("off")
        fig.savefig(out_dir / f"latest_{slice_id}_pdfs.png", dpi=180)
        plt.close(fig)

        summary[slice_id] = {
            "time": float(data["time"]),
            "mach_slice": float(np.sqrt(np.mean(speed * speed, dtype=np.float64))),
            "density_mean": float(np.mean(fields["dens"], dtype=np.float64)),
            "density_std": float(np.std(fields["dens"], dtype=np.float64)),
        }
        del data, fields
        gc.collect()
    return summary


def rfft2(array: np.ndarray, workers: int):
    try:
        from scipy import fft as spfft

        return spfft.rfft2(array, workers=workers)
    except Exception:
        return np.fft.rfft2(array)


def velocity_spectrum_2d(fields: dict[str, np.ndarray], workers: int) -> tuple[np.ndarray, np.ndarray]:
    ny, nx = fields["velx"].shape
    norm = float(nx * ny) ** 2
    power = np.zeros((ny, nx // 2 + 1), dtype=np.float64)
    for name in VELOCITY_NAMES:
        values = fields[name].astype(np.float32, copy=False)
        values = values - np.mean(values, dtype=np.float64)
        fhat = rfft2(values, workers=workers)
        power += np.abs(fhat) ** 2
    power *= 0.5 / norm
    if nx % 2 == 0:
        power[:, 1:-1] *= 2.0
    else:
        power[:, 1:] *= 2.0

    kx = np.fft.rfftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    k = np.sqrt(ky[:, None] ** 2 + kx[None, :] ** 2)
    kbin = np.rint(k).astype(np.int32)
    ncut = min(nx, ny) // 2
    valid = (kbin > 0) & (kbin <= ncut)
    spectrum = np.bincount(kbin[valid].ravel(), weights=power[valid].ravel(), minlength=ncut + 1)
    modes = np.arange(ncut + 1, dtype=np.int32)
    return modes[1:], spectrum[1:].astype(np.float64)


def spectrum_cache_path(out_dir: Path, slice_path: Path) -> Path:
    return out_dir / "cache" / "slice_spectra" / f"{slice_path.stem}.npz"


def compute_or_load_slice_spectrum(
    slice_path: Path, out_dir: Path, workers: int, force: bool
) -> dict[str, np.ndarray | float | int | str]:
    cache = spectrum_cache_path(out_dir, slice_path)
    if cache.exists() and not force:
        cached = np.load(cache)
        return {
            "file": str(slice_path),
            "time": float(cached["time"]),
            "cycle": int(cached["cycle"]),
            "modes": cached["modes"],
            "spectrum": cached["spectrum"],
            "mach_slice": float(cached["mach_slice"]),
            "kinetic_energy_slice": float(cached["kinetic_energy_slice"]),
        }

    data = read_binary(str(slice_path))
    fields = assemble_plane(data, ("velx", "vely", "velz"))
    speed2 = fields["velx"] ** 2 + fields["vely"] ** 2 + fields["velz"] ** 2
    kinetic_energy = 0.5 * float(np.mean(speed2, dtype=np.float64))
    mach = math.sqrt(float(np.mean(speed2, dtype=np.float64)))
    modes, spectrum = velocity_spectrum_2d(fields, workers)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache,
        time=float(data["time"]),
        cycle=int(data["cycle"]),
        modes=modes,
        spectrum=spectrum,
        mach_slice=mach,
        kinetic_energy_slice=kinetic_energy,
    )
    result = {
        "file": str(slice_path),
        "time": float(data["time"]),
        "cycle": int(data["cycle"]),
        "modes": modes,
        "spectrum": spectrum,
        "mach_slice": mach,
        "kinetic_energy_slice": kinetic_energy,
    }
    del data, fields
    gc.collect()
    return result


def rebin_spectrum(
    modes: np.ndarray, spectra: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    max_mode = int(modes[-1])
    ranges: list[tuple[int, int]] = []
    lo = 1
    while lo <= max_mode:
        if lo <= 16:
            width = 1
        elif lo <= 64:
            width = 2
        elif lo <= 128:
            width = 4
        elif lo <= 256:
            width = 8
        else:
            width = 16
        hi = min(max_mode, lo + width - 1)
        ranges.append((lo, hi))
        lo = hi + 1

    centers = np.array([math.sqrt(lo * hi) for lo, hi in ranges], dtype=np.float64)
    widths = np.array([hi - lo + 1 for lo, hi in ranges], dtype=np.float64)
    rebinned = np.full((spectra.shape[0], len(ranges)), np.nan, dtype=np.float64)
    for j, (lo, hi) in enumerate(ranges):
        mask = (modes >= lo) & (modes <= hi)
        rebinned[:, j] = np.nansum(spectra[:, mask], axis=1) / widths[j]
    edges = np.concatenate(([0.5], np.array([hi + 0.5 for _lo, hi in ranges], dtype=np.float64)))
    return centers, rebinned, edges, widths


def center_edges(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if values.size == 1:
        return np.array([values[0] - 0.5, values[0] + 0.5], dtype=np.float64)
    mid = 0.5 * (values[:-1] + values[1:])
    first = values[0] - (mid[0] - values[0])
    last = values[-1] + (values[-1] - mid[-1])
    return np.concatenate(([first], mid, [last]))


def fit_log_slope(modes: np.ndarray, spectrum: np.ndarray, lo: int, hi: int) -> float | None:
    mask = (modes >= lo) & (modes <= hi) & np.isfinite(spectrum) & (spectrum > 0.0)
    if np.count_nonzero(mask) < 4:
        return None
    x = np.log(modes[mask].astype(np.float64))
    y = np.log(spectrum[mask])
    return float(np.polyfit(x, y, 1)[0])


def analyze_slice_spectra(
    run_dir: Path, out_dir: Path, workers: int, force: bool
) -> tuple[dict, dict[str, float | None]]:
    records_by_slice: dict[str, list[dict]] = {}
    for slice_id in SLICE_IDS:
        paths = sorted_files(run_dir / "bin", f"*.{slice_id}.*.bin")
        records_by_slice[slice_id] = [
            compute_or_load_slice_spectrum(path, out_dir, workers, force) for path in paths
        ]

    numbers = sorted(set.intersection(*[set(file_number(Path(r["file"])) for r in records) for records in records_by_slice.values()]))
    combined_records = []
    for number in numbers:
        per_slice = []
        times = []
        machs = []
        for slice_id in SLICE_IDS:
            record = next(r for r in records_by_slice[slice_id] if file_number(Path(r["file"])) == number)
            per_slice.append(record)
            times.append(float(record["time"]))
            machs.append(float(record["mach_slice"]))
        min_n = min(len(r["modes"]) for r in per_slice)
        modes = per_slice[0]["modes"][:min_n]
        stack = np.vstack([r["spectrum"][:min_n] for r in per_slice])
        combined_records.append(
            {
                "number": number,
                "time": float(np.mean(times)),
                "modes": modes,
                "spectrum": np.mean(stack, axis=0),
                "spectrum_sum_three_slices": np.sum(stack, axis=0),
                "mach_slice_mean": float(np.mean(machs)),
                "mach_slice_std": float(np.std(machs)),
            }
        )

    modes = combined_records[0]["modes"]
    spectra = np.vstack([r["spectrum"] for r in combined_records])
    times = np.array([r["time"] for r in combined_records], dtype=np.float64)
    centers, binned, edges, _widths = rebin_spectrum(modes, spectra)

    np.savez_compressed(
        out_dir / "slice_velocity_spectra_over_time.npz",
        times=times,
        modes=modes,
        spectra=spectra,
        binned_modes=centers,
        binned_spectra=binned,
        slice_ids=np.array(SLICE_IDS),
    )

    with (out_dir / "slice_velocity_spectra_over_time.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["mode"] + [f"t={time:.9g}" for time in times])
        for j, mode in enumerate(modes):
            writer.writerow([int(mode)] + [f"{spectra[i, j]:.10e}" for i in range(spectra.shape[0])])

    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(times)))
    fig, ax = plt.subplots(figsize=(8.4, 5.6), constrained_layout=True)
    for i, time in enumerate(times):
        ax.loglog(centers, binned[i], color=colors[i], lw=1.15, alpha=0.9)
    ax.axvspan(1.0, 3.0, color="0.7", alpha=0.18, label="driving modes 1-3")
    ax.set_xlabel("2D slice Fourier mode n")
    ax.set_ylabel("mean slice velocity spectrum")
    ax.set_title("Velocity spectra from three orthogonal slices")
    ax.grid(True, which="both", alpha=0.25)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=times.min(), vmax=times.max()))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("time")
    fig.savefig(out_dir / "slice_velocity_spectra_lines.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.4, 5.6), constrained_layout=True)
    compensated = binned * centers[None, :] ** (5.0 / 3.0)
    for i, time in enumerate(times):
        ax.loglog(centers, compensated[i], color=colors[i], lw=1.15, alpha=0.9)
    ax.axvspan(1.0, 3.0, color="0.7", alpha=0.18)
    ax.set_xlabel("2D slice Fourier mode n")
    ax.set_ylabel(r"$n^{5/3} E_\mathrm{slice}(n)$")
    ax.set_title("Slice spectra compensated by n^(5/3)")
    ax.grid(True, which="both", alpha=0.25)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=times.min(), vmax=times.max()))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("time")
    fig.savefig(out_dir / "slice_velocity_spectra_k53_compensated.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.6), constrained_layout=True)
    positive = binned[np.isfinite(binned) & (binned > 0.0)]
    norm = LogNorm(vmin=max(float(np.percentile(positive, 1.0)), 1.0e-20), vmax=float(np.percentile(positive, 99.5)))
    mesh = ax.pcolormesh(edges, center_edges(times), binned, shading="auto", norm=norm, cmap="magma")
    ax.set_xscale("log")
    ax.set_xlabel("2D slice Fourier mode n")
    ax.set_ylabel("time")
    ax.set_title("Slice velocity spectrum evolution")
    fig.colorbar(mesh, ax=ax, label="mean slice spectrum")
    fig.savefig(out_dir / "slice_velocity_spectrum_heatmap.png", dpi=180)
    plt.close(fig)

    latest = spectra[-1]
    latest_binned = binned[-1]
    slope_summary = {
        "slope_4_16": fit_log_slope(centers, latest_binned, 4, 16),
        "slope_8_32": fit_log_slope(centers, latest_binned, 8, 32),
        "slope_16_64": fit_log_slope(centers, latest_binned, 16, 64),
        "slope_32_128": fit_log_slope(centers, latest_binned, 32, 128),
        "slope_64_256": fit_log_slope(centers, latest_binned, 64, 256),
        "slope_128_512": fit_log_slope(centers, latest_binned, 128, 512),
        "latest_slice_mach_mean": float(combined_records[-1]["mach_slice_mean"]),
        "latest_slice_mach_std": float(combined_records[-1]["mach_slice_std"]),
        "latest_time": float(times[-1]),
        "latest_energy_integral": float(np.sum(latest)),
    }
    return {"combined": combined_records, "by_slice": records_by_slice}, slope_summary


def accumulate_block_stats(data: dict) -> dict[str, float]:
    sums: dict[str, float] = {"count": 0.0}
    sums2: dict[str, float] = {}
    names = ("dens", "velx", "vely", "velz", "tau_xx", "tau_xy", "tau_xz", "tau_yy", "tau_yz", "tau_zz")
    for block in range(int(data["n_mbs"])):
        arrays = {name: np.asarray(data["mb_data"][name][block], dtype=np.float64) for name in names}
        derived = {
            "dens": arrays["dens"],
            "speed2": arrays["velx"] ** 2 + arrays["vely"] ** 2 + arrays["velz"] ** 2,
            "tau_trace": arrays["tau_xx"] + arrays["tau_yy"] + arrays["tau_zz"],
            "tau_norm": np.sqrt(
                arrays["tau_xx"] ** 2
                + arrays["tau_yy"] ** 2
                + arrays["tau_zz"] ** 2
                + 2.0 * (arrays["tau_xy"] ** 2 + arrays["tau_xz"] ** 2 + arrays["tau_yz"] ** 2)
            ),
            "tau_xx": arrays["tau_xx"],
            "tau_xy": arrays["tau_xy"],
            "tau_xz": arrays["tau_xz"],
            "tau_yy": arrays["tau_yy"],
            "tau_yz": arrays["tau_yz"],
            "tau_zz": arrays["tau_zz"],
        }
        n = float(arrays["dens"].size)
        sums["count"] += n
        for name, values in derived.items():
            sums[name] = sums.get(name, 0.0) + float(np.sum(values, dtype=np.float64))
            sums2[name] = sums2.get(name, 0.0) + float(np.sum(values * values, dtype=np.float64))
    count = sums["count"]
    out = {"count": count}
    for name, total in sums.items():
        if name == "count":
            continue
        mean = total / count
        var = max(sums2[name] / count - mean * mean, 0.0)
        out[f"{name}_mean"] = mean
        out[f"{name}_std"] = math.sqrt(var)
    out["mach_coarse"] = math.sqrt(out["speed2_mean"])
    return out


def analyze_sgs_time_stats(run_dir: Path, out_dir: Path) -> dict[int, list[dict[str, float]]]:
    records_by_factor: dict[int, list[dict[str, float]]] = {}
    for factor in SGS_FACTORS:
        subdir = run_dir / f"cbin_sgs_f{factor:03d}_{factor}"
        paths = sorted_files(subdir, f"*.sgs_f{factor:03d}.*.cbin")
        records: list[dict[str, float]] = []
        for path in paths:
            data = read_coarsened_binary(str(path))
            stats = accumulate_block_stats(data)
            stats.update({"time": float(data["time"]), "cycle": float(data["cycle"]), "factor": float(factor)})
            records.append(stats)
            del data
            gc.collect()
        records_by_factor[factor] = records

    fields = [
        ("mach_coarse", "coarse rms Mach"),
        ("tau_trace_mean", "mean tau trace"),
        ("tau_trace_std", "std tau trace"),
        ("tau_norm_mean", "mean tau norm"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.5), constrained_layout=True)
    for ax, (field, ylabel) in zip(axes.ravel(), fields):
        for factor, records in records_by_factor.items():
            time = [r["time"] for r in records]
            values = [r[field] for r in records]
            ax.plot(time, values, marker="o", ms=3, lw=1.4, label=f"f={factor}")
        ax.set_xlabel("time")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
    axes.ravel()[0].legend(fontsize=8)
    fig.savefig(out_dir / "sgs_time_stats.png", dpi=180)
    plt.close(fig)

    with (out_dir / "sgs_time_stats.csv").open("w", newline="") as stream:
        fieldnames = ["factor", "time", "cycle"] + sorted(
            k for k in records_by_factor[SGS_FACTORS[0]][0].keys() if k not in {"factor", "time", "cycle"}
        )
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for factor in SGS_FACTORS:
            for record in records_by_factor[factor]:
                writer.writerow(record)
    return records_by_factor


def central_diff(a: np.ndarray, axis: int, dx: float) -> np.ndarray:
    return (np.roll(a, -1, axis=axis) - np.roll(a, 1, axis=axis)) / (2.0 * dx)


def add_sgs_derived(fields: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    nz, ny, nx = fields["velx"].shape
    dx = 1.0 / float(nx)
    vx = fields["velx"]
    vy = fields["vely"]
    vz = fields["velz"]

    dudx = central_diff(vx, 2, dx)
    dudy = central_diff(vx, 1, dx)
    dudz = central_diff(vx, 0, dx)
    dvdx = central_diff(vy, 2, dx)
    dvdy = central_diff(vy, 1, dx)
    dvdz = central_diff(vy, 0, dx)
    dwdx = central_diff(vz, 2, dx)
    dwdy = central_diff(vz, 1, dx)
    dwdz = central_diff(vz, 0, dx)

    sxx = dudx
    syy = dvdy
    szz = dwdz
    sxy = 0.5 * (dudy + dvdx)
    sxz = 0.5 * (dudz + dwdx)
    syz = 0.5 * (dvdz + dwdy)

    derived = {
        "speed2": vx * vx + vy * vy + vz * vz,
        "speed": np.sqrt(vx * vx + vy * vy + vz * vz),
        "vxvy": vx * vy,
        "vxvz": vx * vz,
        "vyvz": vy * vz,
        "tau_trace": fields["tau_xx"] + fields["tau_yy"] + fields["tau_zz"],
        "tau_norm": np.sqrt(
            fields["tau_xx"] ** 2
            + fields["tau_yy"] ** 2
            + fields["tau_zz"] ** 2
            + 2.0 * (fields["tau_xy"] ** 2 + fields["tau_xz"] ** 2 + fields["tau_yz"] ** 2)
        ),
        "divergence": dudx + dvdy + dwdz,
        "strain_xx": sxx,
        "strain_yy": syy,
        "strain_zz": szz,
        "strain_xy": sxy,
        "strain_xz": sxz,
        "strain_yz": syz,
        "strain_mag": np.sqrt(2.0 * (sxx * sxx + syy * syy + szz * szz + 2.0 * (sxy * sxy + sxz * sxz + syz * syz))),
    }
    derived["Pi_sgs"] = -(
        fields["tau_xx"] * dudx
        + fields["tau_yy"] * dvdy
        + fields["tau_zz"] * dwdz
        + fields["tau_xy"] * (dudy + dvdx)
        + fields["tau_xz"] * (dudz + dwdx)
        + fields["tau_yz"] * (dvdz + dwdy)
    )
    return derived


def sample_arrays(arrays: dict[str, np.ndarray], max_points: int, seed: int) -> dict[str, np.ndarray]:
    first = next(iter(arrays.values()))
    n = first.size
    rng = np.random.default_rng(seed)
    if n > max_points:
        idx = rng.choice(n, size=max_points, replace=False)
    else:
        idx = slice(None)
    return {name: np.asarray(values).ravel()[idx] for name, values in arrays.items()}


def corrcoef(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(ok) < 3:
        return float("nan")
    return float(np.corrcoef(x[ok], y[ok])[0, 1])


def hist2d_panel(ax, x: np.ndarray, y: np.ndarray, xlabel: str, ylabel: str):
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if x.size < 3:
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return
    xlim = np.percentile(x, [0.5, 99.5])
    ylim = np.percentile(y, [0.5, 99.5])
    if xlim[0] == xlim[1]:
        xlim += [-1.0, 1.0]
    if ylim[0] == ylim[1]:
        ylim += [-1.0, 1.0]
    corr = corrcoef(x, y)
    ax.hist2d(x, y, bins=150, range=[xlim, ylim], norm=LogNorm(), cmap="magma")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"r={corr:+.3f}", fontsize=10)


def midplane_pi(fields: dict[str, np.ndarray], k: int) -> np.ndarray:
    nz, ny, nx = fields["velx"].shape
    dx = 1.0 / float(nx)
    kp = (k + 1) % nz
    km = (k - 1) % nz

    vx = fields["velx"]
    vy = fields["vely"]
    vz = fields["velz"]
    vxp = vx[k]
    vyp = vy[k]
    vzp = vz[k]

    dudx = (np.roll(vxp, -1, axis=1) - np.roll(vxp, 1, axis=1)) / (2.0 * dx)
    dudy = (np.roll(vxp, -1, axis=0) - np.roll(vxp, 1, axis=0)) / (2.0 * dx)
    dudz = (vx[kp] - vx[km]) / (2.0 * dx)

    dvdx = (np.roll(vyp, -1, axis=1) - np.roll(vyp, 1, axis=1)) / (2.0 * dx)
    dvdy = (np.roll(vyp, -1, axis=0) - np.roll(vyp, 1, axis=0)) / (2.0 * dx)
    dvdz = (vy[kp] - vy[km]) / (2.0 * dx)

    dwdx = (np.roll(vzp, -1, axis=1) - np.roll(vzp, 1, axis=1)) / (2.0 * dx)
    dwdy = (np.roll(vzp, -1, axis=0) - np.roll(vzp, 1, axis=0)) / (2.0 * dx)
    dwdz = (vz[kp] - vz[km]) / (2.0 * dx)

    return -(
        fields["tau_xx"][k] * dudx
        + fields["tau_yy"][k] * dvdy
        + fields["tau_zz"][k] * dwdz
        + fields["tau_xy"][k] * (dudy + dvdx)
        + fields["tau_xz"][k] * (dudz + dwdx)
        + fields["tau_yz"][k] * (dvdz + dwdy)
    )


def save_sgs_midplane_plot(factor: int, fields: dict[str, np.ndarray], out_dir: Path):
    k = fields["dens"].shape[0] // 2
    speed = np.sqrt(fields["velx"][k] ** 2 + fields["vely"][k] ** 2 + fields["velz"][k] ** 2)
    tau_trace = fields["tau_xx"][k] + fields["tau_yy"][k] + fields["tau_zz"][k]
    tau_norm = np.sqrt(
        fields["tau_xx"][k] ** 2
        + fields["tau_yy"][k] ** 2
        + fields["tau_zz"][k] ** 2
        + 2.0 * (fields["tau_xy"][k] ** 2 + fields["tau_xz"][k] ** 2 + fields["tau_yz"][k] ** 2)
    )
    panels = {
        "density": fields["dens"][k],
        "speed": speed,
        "tau trace": tau_trace,
        "tau norm": tau_norm,
        "tau_xy": fields["tau_xy"][k],
        "tau_xz": fields["tau_xz"][k],
        "tau_yz": fields["tau_yz"][k],
        "Pi_sgs": midplane_pi(fields, k),
    }
    fig, axes = plt.subplots(2, 4, figsize=(15.5, 7.2), constrained_layout=True)
    for ax, (name, values) in zip(axes.ravel(), panels.items()):
        signed = name.startswith("tau_") or name == "Pi_sgs"
        add_image(ax, values, f"f={factor}: {name}", signed=signed, symlog=signed)
    fig.savefig(out_dir / f"sgs_final_midplane_f{factor:03d}.png", dpi=180)
    plt.close(fig)


def save_sgs_pdf_plot(factor: int, sample: dict[str, np.ndarray], out_dir: Path):
    names = ("tau_trace", "tau_norm", "tau_xx", "tau_xy", "tau_xz", "tau_yy", "tau_yz", "tau_zz", "Pi_sgs")
    fig, axes = plt.subplots(3, 3, figsize=(12.0, 9.0), constrained_layout=True)
    for ax, name in zip(axes.ravel(), names):
        values = sample[name]
        values = values[np.isfinite(values)]
        ax.hist(values, bins=180, histtype="step", density=True, lw=1.3)
        ax.set_title(f"f={factor}: {name}", fontsize=10)
        ax.grid(alpha=0.25)
    fig.savefig(out_dir / f"sgs_tau_pdfs_f{factor:03d}_final.png", dpi=180)
    plt.close(fig)


def save_sgs_joint_plot(factor: int, sample: dict[str, np.ndarray], out_dir: Path):
    pairs = (
        ("speed2", "tau_trace"),
        ("vxvy", "tau_xy"),
        ("vxvz", "tau_xz"),
        ("vyvz", "tau_yz"),
        ("divergence", "tau_trace"),
        ("strain_mag", "Pi_sgs"),
    )
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 8.2), constrained_layout=True)
    for ax, (xname, yname) in zip(axes.ravel(), pairs):
        hist2d_panel(ax, sample[xname], sample[yname], xname, yname)
    fig.savefig(out_dir / f"sgs_joint_pdfs_f{factor:03d}_final.png", dpi=180)
    plt.close(fig)


def sample_derivative(a: np.ndarray, k: np.ndarray, j: np.ndarray, i: np.ndarray, axis: int) -> np.ndarray:
    nz, ny, nx = a.shape
    dx = 1.0 / float(nx)
    if axis == 0:
        return (a[(k + 1) % nz, j, i] - a[(k - 1) % nz, j, i]) / (2.0 * dx)
    if axis == 1:
        return (a[k, (j + 1) % ny, i] - a[k, (j - 1) % ny, i]) / (2.0 * dx)
    if axis == 2:
        return (a[k, j, (i + 1) % nx] - a[k, j, (i - 1) % nx]) / (2.0 * dx)
    raise ValueError(f"invalid axis {axis}")


def sample_sgs_fields(fields: dict[str, np.ndarray], max_points: int, seed: int) -> dict[str, np.ndarray]:
    nz, ny, nx = fields["velx"].shape
    total = nz * ny * nx
    rng = np.random.default_rng(seed)
    if total > max_points:
        flat = rng.choice(total, size=max_points, replace=False)
    else:
        flat = np.arange(total, dtype=np.int64)
    k, rem = np.divmod(flat, ny * nx)
    j, i = np.divmod(rem, nx)

    vx = fields["velx"][k, j, i]
    vy = fields["vely"][k, j, i]
    vz = fields["velz"][k, j, i]
    tau_xx = fields["tau_xx"][k, j, i]
    tau_xy = fields["tau_xy"][k, j, i]
    tau_xz = fields["tau_xz"][k, j, i]
    tau_yy = fields["tau_yy"][k, j, i]
    tau_yz = fields["tau_yz"][k, j, i]
    tau_zz = fields["tau_zz"][k, j, i]

    dudx = sample_derivative(fields["velx"], k, j, i, 2)
    dudy = sample_derivative(fields["velx"], k, j, i, 1)
    dudz = sample_derivative(fields["velx"], k, j, i, 0)
    dvdx = sample_derivative(fields["vely"], k, j, i, 2)
    dvdy = sample_derivative(fields["vely"], k, j, i, 1)
    dvdz = sample_derivative(fields["vely"], k, j, i, 0)
    dwdx = sample_derivative(fields["velz"], k, j, i, 2)
    dwdy = sample_derivative(fields["velz"], k, j, i, 1)
    dwdz = sample_derivative(fields["velz"], k, j, i, 0)

    sxx = dudx
    syy = dvdy
    szz = dwdz
    sxy = 0.5 * (dudy + dvdx)
    sxz = 0.5 * (dudz + dwdx)
    syz = 0.5 * (dvdz + dwdy)

    sample = {
        "dens": fields["dens"][k, j, i],
        "speed2": vx * vx + vy * vy + vz * vz,
        "speed": np.sqrt(vx * vx + vy * vy + vz * vz),
        "vxvy": vx * vy,
        "vxvz": vx * vz,
        "vyvz": vy * vz,
        "tau_xx": tau_xx,
        "tau_xy": tau_xy,
        "tau_xz": tau_xz,
        "tau_yy": tau_yy,
        "tau_yz": tau_yz,
        "tau_zz": tau_zz,
        "tau_trace": tau_xx + tau_yy + tau_zz,
        "tau_norm": np.sqrt(tau_xx * tau_xx + tau_yy * tau_yy + tau_zz * tau_zz + 2.0 * (tau_xy * tau_xy + tau_xz * tau_xz + tau_yz * tau_yz)),
        "divergence": dudx + dvdy + dwdz,
        "strain_xx": sxx,
        "strain_xy": sxy,
        "strain_xz": sxz,
        "strain_yy": syy,
        "strain_yz": syz,
        "strain_zz": szz,
        "strain_mag": np.sqrt(2.0 * (sxx * sxx + syy * syy + szz * szz + 2.0 * (sxy * sxy + sxz * sxz + syz * syz))),
    }
    sample["Pi_sgs"] = -(
        tau_xx * dudx
        + tau_yy * dvdy
        + tau_zz * dwdz
        + tau_xy * (dudy + dvdx)
        + tau_xz * (dudz + dwdx)
        + tau_yz * (dvdz + dwdy)
    )
    return sample


def analyze_sgs_final_correlations(
    run_dir: Path, out_dir: Path, sample_points: int
) -> tuple[dict[int, dict[str, float]], list[dict[str, float | int | str]]]:
    final_stats: dict[int, dict[str, float]] = {}
    corr_rows: list[dict[str, float | int | str]] = []
    pairs = (
        ("speed2", "tau_trace"),
        ("vxvy", "tau_xy"),
        ("vxvz", "tau_xz"),
        ("vyvz", "tau_yz"),
        ("divergence", "tau_trace"),
        ("strain_xx", "tau_xx"),
        ("strain_xy", "tau_xy"),
        ("strain_xz", "tau_xz"),
        ("strain_yy", "tau_yy"),
        ("strain_yz", "tau_yz"),
        ("strain_zz", "tau_zz"),
        ("strain_mag", "Pi_sgs"),
        ("tau_norm", "Pi_sgs"),
    )

    for factor in SGS_FACTORS:
        subdir = run_dir / f"cbin_sgs_f{factor:03d}_{factor}"
        path = sorted_files(subdir, f"*.sgs_f{factor:03d}.*.cbin")[-1]
        data = read_coarsened_binary(str(path))
        output_time = float(data["time"])
        fields = assemble_3d(data, SGS_NAMES)
        del data
        gc.collect()
        sample = sample_sgs_fields(fields, sample_points, seed=1000 + factor)

        stats: dict[str, float] = {"time": output_time, "sample_points": float(next(iter(sample.values())).size)}
        for name in (
            "dens",
            "speed",
            "speed2",
            "tau_trace",
            "tau_norm",
            "tau_xx",
            "tau_xy",
            "tau_xz",
            "tau_yy",
            "tau_yz",
            "tau_zz",
            "Pi_sgs",
            "divergence",
            "strain_mag",
        ):
            values = sample[name]
            stats[f"{name}_mean"] = float(np.mean(values, dtype=np.float64))
            stats[f"{name}_std"] = float(np.std(values, dtype=np.float64))
            stats[f"{name}_p01"] = float(np.percentile(values, 1.0))
            stats[f"{name}_p50"] = float(np.percentile(values, 50.0))
            stats[f"{name}_p99"] = float(np.percentile(values, 99.0))
        final_stats[factor] = stats
        for xname, yname in pairs:
            corr_rows.append(
                {
                    "factor": factor,
                    "time": stats["time"],
                    "x": xname,
                    "y": yname,
                    "corr": corrcoef(sample[xname], sample[yname]),
                }
            )

        save_sgs_midplane_plot(factor, fields, out_dir)
        save_sgs_pdf_plot(factor, sample, out_dir)
        save_sgs_joint_plot(factor, sample, out_dir)
        del fields, sample
        gc.collect()

    with (out_dir / "sgs_final_correlations.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["factor", "time", "x", "y", "corr"])
        writer.writeheader()
        writer.writerows(corr_rows)

    with (out_dir / "sgs_final_stats.csv").open("w", newline="") as stream:
        fieldnames = ["factor"] + sorted(next(iter(final_stats.values())).keys())
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for factor, stats in final_stats.items():
            writer.writerow({"factor": factor, **stats})

    return final_stats, corr_rows


def write_summary(
    out_dir: Path,
    inventory: dict,
    history_summary: dict,
    slice_summary: dict,
    spectra_summary: dict,
    sgs_time_stats: dict[int, list[dict[str, float]]],
    sgs_final_stats: dict[int, dict[str, float]],
    corr_rows: list[dict[str, float | int | str]],
):
    summary = {
        "inventory": inventory,
        "history": history_summary,
        "latest_slices": slice_summary,
        "slice_spectra": spectra_summary,
        "sgs_latest": sgs_final_stats,
        "sgs_correlations": corr_rows,
    }
    (out_dir / "analysis_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    strongest = sorted(
        corr_rows,
        key=lambda row: abs(float(row["corr"])) if np.isfinite(float(row["corr"])) else -1.0,
        reverse=True,
    )[:10]

    lines = [
        f"# Analysis summary: {RUN_NAME}",
        "",
        "## Run completion",
        f"- Final time: `{history_summary['final_time']:.6g}`",
        f"- Final rms Mach from history: `{history_summary['final_mach']:.6g}`",
        f"- Target Mach crossing time: `{history_summary['target_mach_crossing_time']}`",
        f"- Mean Mach for `t >= 2`: `{history_summary['mean_mach_t_ge_2']:.6g} +/- {history_summary['std_mach_t_ge_2']:.3g}`",
        "",
        "## Output inventory",
        f"- Full primitive cubes: `{inventory['prim_count']}`",
        f"- Full-resolution slices: `{inventory['slice_counts']}`",
        f"- SGS dumps per factor: `{inventory['sgs_counts']}`",
        "",
        "## Slice spectra",
        "- Spectra use all three velocity components on each orthogonal 2D slice.",
        "- The plotted combined spectrum is the mean of the three slice spectra.",
        f"- Latest mean slice Mach: `{spectra_summary['latest_slice_mach_mean']:.6g} +/- {spectra_summary['latest_slice_mach_std']:.3g}`",
        f"- Latest slopes: `4-16={spectra_summary['slope_4_16']}`, `8-32={spectra_summary['slope_8_32']}`, `16-64={spectra_summary['slope_16_64']}`, `32-128={spectra_summary['slope_32_128']}`",
        "",
        "## Latest SGS means",
    ]
    for factor in SGS_FACTORS:
        stats = sgs_final_stats[factor]
        lines.append(
            f"- f={factor}: mean tau_trace `{stats['tau_trace_mean']:.6e}`, "
            f"std tau_trace `{stats['tau_trace_std']:.6e}`, "
            f"mean Pi_sgs `{stats['Pi_sgs_mean']:.6e}`, std Pi_sgs `{stats['Pi_sgs_std']:.6e}`"
        )
    lines.extend(["", "## Strongest sampled SGS correlations"])
    for row in strongest:
        lines.append(
            f"- f={row['factor']} `{row['x']}` vs `{row['y']}`: r=`{float(row['corr']):+.4f}`"
        )
    lines.extend(
        [
            "",
            "## Key files",
            "- `history_mach_energy.png`",
            "- `slice_velocity_spectra_lines.png`",
            "- `slice_velocity_spectra_k53_compensated.png`",
            "- `slice_velocity_spectrum_heatmap.png`",
            "- `sgs_time_stats.png`",
            "- `sgs_final_correlations.csv`",
        ]
    )
    (out_dir / "analysis_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    history_path = args.run_dir / f"{args.run_dir.name}.hydro.hst"
    history = read_history(history_path)
    history_summary = plot_history(history, args.out_dir)

    inventory = {
        "run_dir": str(args.run_dir),
        "prim_count": len(list((args.run_dir / "bin").glob("*.prim.*.bin"))),
        "slice_counts": {
            slice_id: len(list((args.run_dir / "bin").glob(f"*.{slice_id}.*.bin"))) for slice_id in SLICE_IDS
        },
        "sgs_counts": {
            factor: len(list((args.run_dir / f"cbin_sgs_f{factor:03d}_{factor}").glob(f"*.sgs_f{factor:03d}.*.cbin")))
            for factor in SGS_FACTORS
        },
    }

    slice_summary = save_latest_slice_plots(args.run_dir, args.out_dir)
    _records, spectra_summary = analyze_slice_spectra(args.run_dir, args.out_dir, args.workers, args.force)
    sgs_time_stats = analyze_sgs_time_stats(args.run_dir, args.out_dir)
    sgs_final_stats, corr_rows = analyze_sgs_final_correlations(args.run_dir, args.out_dir, args.sample_points)

    write_summary(
        args.out_dir,
        inventory,
        history_summary,
        slice_summary,
        spectra_summary,
        sgs_time_stats,
        sgs_final_stats,
        corr_rows,
    )
    print(f"Wrote analysis to {args.out_dir}")


if __name__ == "__main__":
    main()
