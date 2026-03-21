#!/usr/bin/env python3
"""Generate the scalar-mixing velocity-method figures for the methods note."""

from __future__ import annotations

import argparse
import gc
import glob
import json
import math
import re
import shlex
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm

from bin_convert_new import read_binary_as_athdf


EXPO = 5.0 / 3.0
KMIN = 2
KMAX = 128
KREF = 3
TAIL_KMIN = 80
SECTOR_ORDER = (
    ("local_local", "local-local", "#234f7d"),
    ("low_high", "low-high", "#b24a2f"),
    ("high_low", "high-low", "#5c8c2b"),
    ("high_high_cancel", "high-high cancel", "#7a5fa1"),
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_athena_path() -> Path:
    for relpath in (
        Path("build-mpi/src/athena"),
        Path("build-codex-spectrum/src/athena"),
        Path("build/src/athena"),
    ):
        candidate = _repo_root() / relpath
        if candidate.exists():
            return candidate
    return _repo_root() / "build-mpi" / "src" / "athena"


def _figure_root() -> Path:
    return _repo_root() / "docs" / "figures" / "scalar_mixing_velocity_methods"


def _default_clebsch_audit_json() -> Path | None:
    candidate = _repo_root() / "tmp" / "clebsch_audit_report_prod256.json"
    return candidate if candidate.exists() else None


def _case_specs() -> list[dict[str, object]]:
    return [
        {
            "name": "2D Projection",
            "slug": "proj2d",
            "input": "inputs/hydro/scalar_mixing_methods_paper_2d.athinput",
            "basename": "scalar_mix_methods_proj2d",
            "overrides": [],
            "launcher": "mpirun -np 4",
            "dim": "2d",
            "method": "projection",
            "color": "#0b6e4f",
        },
        {
            "name": "2D Stream",
            "slug": "stream2d",
            "input": "inputs/hydro/scalar_mixing_methods_paper_2d.athinput",
            "basename": "scalar_mix_methods_stream2d",
            "overrides": ["problem/turb_use_stream_function=true"],
            "launcher": "mpirun -np 4",
            "dim": "2d",
            "method": "stream",
            "color": "#bd632f",
        },
        {
            "name": "3D Projection",
            "slug": "proj3d",
            "input": "inputs/hydro/scalar_mixing_methods_paper_3d.athinput",
            "basename": "scalar_mix_methods_proj3d",
            "overrides": [],
            "launcher": "mpirun -np 8",
            "dim": "3d",
            "method": "projection",
            "color": "#275dad",
        },
        {
            "name": "3D Stream",
            "slug": "stream3d",
            "input": "inputs/hydro/scalar_mixing_methods_paper_3d.athinput",
            "basename": "scalar_mix_methods_stream3d",
            "overrides": ["problem/turb_use_stream_function=true"],
            "launcher": "mpirun -np 8",
            "dim": "3d",
            "method": "stream",
            "color": "#8b3fb3",
        },
    ]


def _parse_init_seconds(log_text: str) -> float | None:
    match = re.search(r"velocity_init_wall_seconds=([0-9.eE+-]+)", log_text)
    if match is None:
        return None
    return float(match.group(1))


def _run_case(athena: Path, case: dict[str, object], force: bool) -> tuple[Path, Path, float | None]:
    athena = athena.resolve()
    run_dir = athena.parent.parent if athena.parent.name == "src" else athena.parent
    output_dir = run_dir / "bin"
    output_dir.mkdir(parents=True, exist_ok=True)
    basename = str(case["basename"])
    hydro_pattern = output_dir / f"{basename}.hydro_w.*.bin"
    diag_path = run_dir / f"{basename}.turb_init_diag.json"
    timing_path = run_dir / f"{basename}.velocity_init.json"

    hydro_matches = sorted(glob.glob(str(hydro_pattern)))
    if hydro_matches and diag_path.exists() and not force:
        init_seconds = None
        if timing_path.exists():
            timing = json.loads(timing_path.read_text(encoding="utf-8"))
            if timing.get("init_seconds") is not None:
                init_seconds = float(timing["init_seconds"])
        return Path(hydro_matches[-1]), diag_path, init_seconds

    for path in hydro_matches:
        Path(path).unlink()
    if diag_path.exists():
        diag_path.unlink()
    if timing_path.exists():
        timing_path.unlink()

    launcher = shlex.split(str(case["launcher"]))
    cmd = launcher + [
        str(athena),
        "-i",
        str(_repo_root() / str(case["input"])),
        f"job/basename={basename}",
    ] + list(case["overrides"])
    completed = subprocess.run(
        cmd,
        cwd=run_dir,
        capture_output=True,
        text=True,
        check=True,
    )
    log_text = completed.stdout + completed.stderr
    init_seconds = _parse_init_seconds(log_text)
    timing_path.write_text(
        json.dumps(
            {
                "case": str(case["name"]),
                "basename": basename,
                "command": cmd,
                "init_seconds": init_seconds,
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    hydro_matches = sorted(glob.glob(str(hydro_pattern)))
    if not hydro_matches:
        raise RuntimeError(f"Missing hydro_w output for {case['name']}")
    if not diag_path.exists():
        raise RuntimeError(f"Missing diagnostics JSON for {case['name']}")
    return Path(hydro_matches[-1]), diag_path, init_seconds


def _load_velocity(path: Path) -> dict[str, np.ndarray]:
    data = read_binary_as_athdf(
        str(path),
        quantities=["velx", "vely", "velz"],
        dtype=np.float64,
    )
    return {
        "velx": np.asarray(data["velx"], dtype=np.float64),
        "vely": np.asarray(data["vely"], dtype=np.float64),
        "velz": np.asarray(data["velz"], dtype=np.float64),
    }


def _fft_wave_numbers(n: int, xmin: float, xmax: float) -> np.ndarray:
    length = xmax - xmin
    dx = length / float(n)
    return 2.0 * np.pi * np.fft.fftfreq(n, d=dx)


def _integer_mode_numbers(n: int) -> np.ndarray:
    return np.fft.fftfreq(n) * n


def _midplane(field: np.ndarray) -> np.ndarray:
    if field.shape[0] == 1:
        return field[0]
    return field[field.shape[0] // 2]


def _shell_spectrum(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ntot = float(vx.size)
    uhat_x = np.fft.fftn(vx)
    uhat_y = np.fft.fftn(vy)
    uhat_z = np.fft.fftn(vz)
    energy = (np.abs(uhat_x) ** 2 + np.abs(uhat_y) ** 2 + np.abs(uhat_z) ** 2) / (ntot * ntot)

    nz, ny, nx = vx.shape
    kz = _integer_mode_numbers(nz) if nz > 1 else np.array([0.0])
    ky = _integer_mode_numbers(ny) if ny > 1 else np.array([0.0])
    kx = _integer_mode_numbers(nx)
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    shell = np.ceil(np.sqrt(kx3 * kx3 + ky3 * ky3 + kz3 * kz3) - 1.0e-12).astype(np.int64)
    shell_energy = np.bincount(shell.ravel(), weights=energy.ravel())
    kvals = np.arange(shell_energy.size, dtype=np.float64)
    return kvals, shell_energy


def _spectral_divergence_ratio(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray,
                               meta: dict[str, object]) -> float:
    uhat_x = np.fft.fftn(vx)
    uhat_y = np.fft.fftn(vy)
    uhat_z = np.fft.fftn(vz)
    mesh = meta["mesh"]
    nz, ny, nx = vx.shape
    kx = _fft_wave_numbers(nx, float(mesh["x1min"]), float(mesh["x1max"]))
    ky = _fft_wave_numbers(ny, float(mesh["x2min"]), float(mesh["x2max"])) if ny > 1 else np.array([0.0])
    kz = _fft_wave_numbers(nz, float(mesh["x3min"]), float(mesh["x3max"])) if nz > 1 else np.array([0.0])
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    div_hat = 1j * (kx3 * uhat_x + ky3 * uhat_y + kz3 * uhat_z)
    numer = np.sum(np.abs(div_hat) ** 2)
    denom = np.sum((kx3 * kx3 + ky3 * ky3 + kz3 * kz3) *
                   (np.abs(uhat_x) ** 2 + np.abs(uhat_y) ** 2 + np.abs(uhat_z) ** 2))
    if denom <= 0.0:
        return 0.0
    return float(np.sqrt(numer / denom))


def _fit_slope(kvals: np.ndarray, spectrum: np.ndarray, kmin: int, kmax: int) -> float:
    mask = (kvals >= kmin) & (kvals <= kmax) & (spectrum > 0.0)
    if np.count_nonzero(mask) < 2:
        return float("nan")
    x = np.log(kvals[mask])
    y = np.log(spectrum[mask])
    coeff = np.polyfit(x, y, 1)
    return float(coeff[0])


def _tail_smoothness(kvals: np.ndarray, spectrum: np.ndarray,
                     kmin: int = TAIL_KMIN, kmax: int = KMAX) -> float:
    usable = spectrum.copy()
    usable[usable <= 0.0] = np.nan
    residuals = []
    for k in range(max(kmin, 2), min(kmax, usable.size - 2) + 1):
        if not np.isfinite(usable[k - 1]) or not np.isfinite(usable[k]) or not np.isfinite(usable[k + 1]):
            continue
        local_mean = 0.5 * (np.log(usable[k - 1]) + np.log(usable[k + 1]))
        residuals.append(abs(np.log(usable[k]) - local_mean))
    if not residuals:
        return 0.0
    return float(np.max(residuals))


def _shape_normalize(spectrum: np.ndarray, kref: int = KREF) -> np.ndarray:
    normalized = np.full_like(spectrum, np.nan, dtype=np.float64)
    if kref >= spectrum.size or spectrum[kref] <= 0.0:
        return normalized
    normalized[:] = spectrum / spectrum[kref]
    return normalized


def _reference_law(kvals: np.ndarray, slope: float = -EXPO, kref: int = KREF) -> np.ndarray:
    ref = np.full_like(kvals, np.nan, dtype=np.float64)
    if kref <= 0:
        return ref
    mask = kvals > 0.0
    ref[mask] = (kvals[mask] / float(kref)) ** slope
    return ref


def _sector_fraction_series(diag: dict[str, object], response_key: str) -> dict[str, np.ndarray]:
    sector_block = diag["sector_response"][response_key]
    arrays = {
        key: np.asarray(sector_block[key]["shell_energy"], dtype=np.float64)
        for key, _, _ in SECTOR_ORDER
    }
    total = np.zeros_like(next(iter(arrays.values())))
    for values in arrays.values():
        total = total + values
    fractions = {}
    mask = total > 0.0
    for key, values in arrays.items():
        frac = np.full_like(values, np.nan, dtype=np.float64)
        frac[mask] = values[mask] / total[mask]
        fractions[key] = frac
    return fractions


def _plot_sector_fractions(ax: plt.Axes, kvals: np.ndarray, diag: dict[str, object],
                           response_key: str, title: str) -> None:
    fractions = _sector_fraction_series(diag, response_key)
    kstop = min(KMAX, len(kvals) - 1)
    for key, label, color in SECTOR_ORDER:
        ax.semilogx(kvals[1:kstop + 1], fractions[key][1:kstop + 1],
                    lw=2.0, color=color, label=label)
    ax.set_xlim(1.0, KMAX)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(r"$k$")
    ax.set_ylabel("fraction of shell response")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(fontsize=8, loc="best")


def _velocity_rms(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> float:
    return float(np.sqrt(np.mean(vx * vx + vy * vy + vz * vz)))


def _vorticity_2d(vx: np.ndarray, vy: np.ndarray, meta: dict[str, object]) -> np.ndarray:
    ny, nx = vx.shape
    mesh = meta["mesh"]
    kx = _fft_wave_numbers(nx, float(mesh["x1min"]), float(mesh["x1max"]))
    ky = _fft_wave_numbers(ny, float(mesh["x2min"]), float(mesh["x2max"]))
    ky2, kx2 = np.meshgrid(ky, kx, indexing="ij")
    vx_hat = np.fft.fft2(vx)
    vy_hat = np.fft.fft2(vy)
    omega_hat = 1j * kx2 * vy_hat - 1j * ky2 * vx_hat
    return np.fft.ifft2(omega_hat).real


def _extent(meta: dict[str, object]) -> list[float]:
    mesh = meta["mesh"]
    return [
        float(mesh["x1min"]),
        float(mesh["x1max"]),
        float(mesh["x2min"]),
        float(mesh["x2max"]),
    ]


def _scalar_mode_amplitude(diag: dict[str, object], coeff_key: str) -> np.ndarray:
    coeffs = diag["coefficients"][coeff_key]
    aka = np.asarray(coeffs["aka"], dtype=np.float64)
    akb = np.asarray(coeffs["akb"], dtype=np.float64)
    return np.sqrt(0.5 * (aka * aka + akb * akb))


def _vector_mode_amplitude(diag: dict[str, object], coeff_key: str) -> np.ndarray:
    coeffs = diag["coefficients"][coeff_key]
    aka0 = np.asarray(coeffs["aka0"], dtype=np.float64)
    aka1 = np.asarray(coeffs["aka1"], dtype=np.float64)
    aka2 = np.asarray(coeffs["aka2"], dtype=np.float64)
    akb0 = np.asarray(coeffs["akb0"], dtype=np.float64)
    akb1 = np.asarray(coeffs["akb1"], dtype=np.float64)
    akb2 = np.asarray(coeffs["akb2"], dtype=np.float64)
    return np.sqrt(0.5 * (aka0 * aka0 + aka1 * aka1 + aka2 * aka2 +
                          akb0 * akb0 + akb1 * akb1 + akb2 * akb2))


def _catalog_shell_counts(diag: dict[str, object], catalog_key: str, nmax: int) -> np.ndarray:
    shells = np.asarray(diag["catalogs"][catalog_key]["shell"], dtype=np.int64)
    counts = np.bincount(shells, minlength=nmax + 1)
    if counts.size < nmax + 1:
        counts = np.pad(counts, (0, nmax + 1 - counts.size))
    return counts[:nmax + 1]


def _reconstruct_slice(diag: dict[str, object], catalog_key: str,
                       coeff_key: str, slice_key: str) -> np.ndarray:
    meta = diag["metadata"]
    mesh = meta["mesh"]
    catalog = diag["catalogs"][catalog_key]
    coeffs = diag["coefficients"][coeff_key]
    slice_obj = diag["slices"][slice_key]
    ny, nx = slice_obj["shape"]
    x3_slice = float(slice_obj["x3_slice"])
    dkx = 2.0 * np.pi / (float(mesh["x1max"]) - float(mesh["x1min"]))
    dky = 2.0 * np.pi / (float(mesh["x2max"]) - float(mesh["x2min"]))
    dkz = 0.0
    if int(mesh["global_nx3"]) > 1:
        dkz = 2.0 * np.pi / (float(mesh["x3max"]) - float(mesh["x3min"]))
    nkx = np.asarray(catalog["nkx"], dtype=np.float64)
    nky = np.asarray(catalog["nky"], dtype=np.float64)
    nkz = np.asarray(catalog["nkz"], dtype=np.float64)
    aka = np.asarray(coeffs["aka"], dtype=np.float64)
    akb = np.asarray(coeffs["akb"], dtype=np.float64)
    x = np.linspace(float(mesh["x1min"]), float(mesh["x1max"]), nx, endpoint=False) + 0.5 * (
        float(mesh["x1max"]) - float(mesh["x1min"])
    ) / nx
    y = np.linspace(float(mesh["x2min"]), float(mesh["x2max"]), ny, endpoint=False) + 0.5 * (
        float(mesh["x2max"]) - float(mesh["x2min"])
    ) / ny
    yy, xx = np.meshgrid(y, x, indexing="ij")
    field = np.zeros((ny, nx), dtype=np.float64)
    chunk = 256
    for start in range(0, aka.size, chunk):
        stop = min(start + chunk, aka.size)
        phase = (
            dkx * nkx[start:stop, None, None] * xx[None, :, :]
            + dky * nky[start:stop, None, None] * yy[None, :, :]
            + dkz * nkz[start:stop, None, None] * x3_slice
        )
        field += np.sum(
            aka[start:stop, None, None] * np.cos(phase)
            - akb[start:stop, None, None] * np.sin(phase),
            axis=0,
        )
    return field


def _extract_case(velocity_path: Path, diag_path: Path, case: dict[str, object],
                  init_seconds: float | None) -> dict[str, object]:
    diag = json.loads(diag_path.read_text(encoding="utf-8"))
    velocity = _load_velocity(velocity_path)
    vx = velocity["velx"]
    vy = velocity["vely"]
    vz = velocity["velz"]
    kvals, shell_energy = _shell_spectrum(vx, vy, vz)
    meta = diag["metadata"]
    slice_vx = _midplane(vx)
    slice_vy = _midplane(vy)
    slice_vz = _midplane(vz)
    speed = np.sqrt(slice_vx * slice_vx + slice_vy * slice_vy + slice_vz * slice_vz)
    case_data: dict[str, object] = {
        "case": case,
        "velocity_path": velocity_path,
        "diag_path": diag_path,
        "diag": diag,
        "meta": meta,
        "extent": _extent(meta),
        "slice_vx": slice_vx,
        "slice_vy": slice_vy,
        "slice_vz": slice_vz,
        "slice_speed": speed,
        "kvals": kvals,
        "shell_energy": shell_energy,
        "shell_shape": _shape_normalize(shell_energy),
        "divergence_ratio": _spectral_divergence_ratio(vx, vy, vz, meta),
        "v_rms": _velocity_rms(vx, vy, vz),
        "slope": _fit_slope(kvals, shell_energy, KMIN, KMAX),
        "init_seconds": init_seconds if init_seconds is not None else meta.get("init_wall_seconds"),
    }
    if str(case["dim"]) == "2d":
        case_data["vorticity"] = _vorticity_2d(slice_vx, slice_vy, meta)
    if str(case["slug"]) == "stream3d":
        case_data["tail_smoothness"] = _tail_smoothness(kvals, shell_energy)
        quality_gate = diag["quality_gate"]
        if not bool(quality_gate["pass"]) or bool(quality_gate["forced_continue"]):
            raise RuntimeError(
                "3D stream figure generation requires a passing quality gate without forced continuation."
            )

    scale = float(meta["normalization"]["scale"])
    if str(case["method"]) == "projection":
        diag_velocity_shell = np.asarray(diag["spectra"]["velocity_shell_energy"], dtype=np.float64)
    elif str(case["slug"]) == "stream2d":
        diag_velocity_shell = np.asarray(diag["spectra"]["velocity_shell_energy"], dtype=np.float64)
    else:
        diag_velocity_shell = np.asarray(
            diag["spectra"]["realized_retained_velocity_shell_energy_raw"],
            dtype=np.float64,
        )
    diag_velocity_shell = diag_velocity_shell * (scale * scale)
    case_data["diag_velocity_shell"] = diag_velocity_shell
    case_data["diag_velocity_shape"] = _shape_normalize(diag_velocity_shell)
    valid_band = slice(KMIN, min(KMAX, len(shell_energy) - 1) + 1)
    band_fft = _shape_normalize(shell_energy)[valid_band]
    band_diag = _shape_normalize(diag_velocity_shell)[valid_band]
    mask = np.isfinite(band_fft) & np.isfinite(band_diag)
    if np.any(mask):
        case_data["velocity_diag_match_rms"] = float(
            np.sqrt(np.mean((np.log(band_fft[mask]) - np.log(band_diag[mask])) ** 2))
        )
    else:
        case_data["velocity_diag_match_rms"] = float("nan")

    if str(case["method"]) == "stream":
        if str(case["dim"]) == "2d":
            stored = np.asarray(diag["slices"]["psi_xy"]["data"], dtype=np.float64).reshape(
                diag["slices"]["psi_xy"]["shape"]
            )
            recon = _reconstruct_slice(diag, "psi", "psi", "psi_xy")
            case_data["psi"] = stored
            case_data["psi_reconstruct_max_abs"] = float(np.max(np.abs(stored - recon)))
            case_data["psi_shell_energy"] = np.asarray(diag["spectra"]["psi_shell_energy"], dtype=np.float64)
        else:
            phi1 = np.asarray(diag["slices"]["phi1_xy_mid"]["data"], dtype=np.float64).reshape(
                diag["slices"]["phi1_xy_mid"]["shape"]
            )
            phi2 = np.asarray(diag["slices"]["phi2_xy_mid"]["data"], dtype=np.float64).reshape(
                diag["slices"]["phi2_xy_mid"]["shape"]
            )
            recon1 = _reconstruct_slice(diag, "phi1", "phi1", "phi1_xy_mid")
            recon2 = _reconstruct_slice(diag, "phi2", "phi2", "phi2_xy_mid")
            case_data["phi1"] = phi1
            case_data["phi2"] = phi2
            case_data["phi1_reconstruct_max_abs"] = float(np.max(np.abs(phi1 - recon1)))
            case_data["phi2_reconstruct_max_abs"] = float(np.max(np.abs(phi2 - recon2)))
            case_data["phi1_shell_energy"] = np.asarray(diag["spectra"]["phi1_shell_energy"], dtype=np.float64)
            case_data["phi2_shell_energy"] = np.asarray(diag["spectra"]["phi2_shell_energy"], dtype=np.float64)

    del velocity, vx, vy, vz
    gc.collect()
    return case_data


def _overview_figure(results: dict[str, dict[str, object]], outdir: Path) -> None:
    ordered = [results[key] for key in ("proj2d", "stream2d", "proj3d", "stream3d")]
    vmax = max(float(np.nanmax(item["slice_speed"])) for item in ordered)
    fig = plt.figure(figsize=(16, 10), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.15])
    slice_axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
    ]
    spec_ax = fig.add_subplot(grid[:, 2])
    ref_k = np.arange(1, KMAX + 1, dtype=np.float64)
    ref = _reference_law(ref_k)
    image = None
    for ax, item in zip(slice_axes, ordered):
        image = ax.imshow(
            item["slice_speed"],
            origin="lower",
            extent=item["extent"],
            cmap="magma",
            vmin=0.0,
            vmax=vmax,
        )
        skip = max(1, item["slice_vx"].shape[0] // 20)
        x = np.linspace(item["extent"][0], item["extent"][1], item["slice_vx"].shape[1], endpoint=False)
        y = np.linspace(item["extent"][2], item["extent"][3], item["slice_vx"].shape[0], endpoint=False)
        xx, yy = np.meshgrid(x[::skip], y[::skip])
        ax.quiver(
            xx,
            yy,
            item["slice_vx"][::skip, ::skip],
            item["slice_vy"][::skip, ::skip],
            color="white",
            scale=20.0,
            width=0.0025,
            pivot="mid",
            alpha=0.8,
        )
        ax.set_title(
            f"{item['case']['name']}\n"
            f"t_init={float(item['init_seconds']):.2f} s, "
            f"div={float(item['divergence_ratio']):.2e}"
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    for item in ordered:
        shell = np.asarray(item["shell_shape"])
        spec_ax.loglog(
            item["kvals"][1:KMAX + 1],
            shell[1:KMAX + 1],
            lw=2.0,
            label=item["case"]["name"],
            color=item["case"]["color"],
        )
    spec_ax.loglog(ref_k, ref, color="black", ls="--", lw=1.5,
                   label=r"reference $E_u(k)/E_u(3) = (k/3)^{-5/3}$")
    spec_ax.set_xlim(1.0, KMAX)
    spec_ax.set_ylim(1.0e-4, 4.0)
    spec_ax.set_xlabel(r"$k$")
    spec_ax.set_ylabel(r"$E_u(k) / E_u(3)$")
    spec_ax.grid(True, which="both", alpha=0.2)
    spec_ax.legend(fontsize=8, loc="lower left")
    fig.colorbar(image, ax=slice_axes, shrink=0.92, location="right", label=r"$|\mathbf{v}|$")
    fig.suptitle(r"Prescribed velocity initialization in scalar_mixing: $512^2$ in 2D, $256^3$ in 3D")
    fig.savefig(outdir / "overview.png", dpi=220)
    plt.close(fig)


def _plot_spectrum(ax: plt.Axes, kvals: np.ndarray, spectra: list[tuple[str, np.ndarray, str]],
                   normalize: bool, ylabel: str, ylim: tuple[float, float]) -> None:
    for label, spec, color in spectra:
        curve = _shape_normalize(spec) if normalize else spec
        ax.loglog(kvals[1:KMAX + 1], curve[1:KMAX + 1], label=label, lw=2.0, color=color)
    if normalize:
        ref = _reference_law(kvals[1:KMAX + 1])
        ax.loglog(kvals[1:KMAX + 1], ref, color="black", ls="--", lw=1.2,
                  label=r"$k^{-5/3}$")
    ax.set_xlim(1.0, KMAX)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.2)


def _occupancy_scatter(ax: plt.Axes, diag: dict[str, object], catalog_key: str,
                       amplitude: np.ndarray, title: str) -> None:
    catalog = diag["catalogs"][catalog_key]
    nkx = np.asarray(catalog["nkx"], dtype=np.float64)
    nky = np.asarray(catalog["nky"], dtype=np.float64)
    scatter = ax.scatter(nkx, nky, c=amplitude, s=10.0, cmap="viridis")
    ax.set_xlabel(r"$n_x$")
    ax.set_ylabel(r"$n_y$")
    ax.set_title(title)
    plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04, label="mode amplitude")


def _shell_count_panel(ax: plt.Axes, diag: dict[str, object], catalog_key: str,
                       nmax: int, title: str) -> None:
    counts = _catalog_shell_counts(diag, catalog_key, nmax)
    ax.bar(np.arange(1, nmax + 1), counts[1:], width=0.85, color="#587498")
    ax.set_xlim(1, nmax)
    ax.set_xlabel(r"$k$")
    ax.set_ylabel("retained modes")
    ax.set_title(title)


def _methods_2d_figure(results: dict[str, dict[str, object]], outdir: Path) -> None:
    proj = results["proj2d"]
    stream = results["stream2d"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    vmax = max(float(np.nanmax(proj["slice_speed"])), float(np.nanmax(stream["slice_speed"])))
    wmax = max(float(np.nanmax(np.abs(proj["vorticity"]))), float(np.nanmax(np.abs(stream["vorticity"]))))
    im0 = axes[0, 0].imshow(proj["slice_speed"], origin="lower", extent=proj["extent"],
                            cmap="magma", vmin=0.0, vmax=vmax)
    axes[0, 0].set_title("2D projection: speed")
    im1 = axes[0, 1].imshow(proj["vorticity"], origin="lower", extent=proj["extent"],
                            cmap="coolwarm", norm=TwoSlopeNorm(vcenter=0.0, vmin=-wmax, vmax=wmax))
    axes[0, 1].set_title("2D projection: vorticity")
    _plot_spectrum(
        axes[0, 2],
        proj["kvals"],
        [("velocity", proj["shell_energy"], proj["case"]["color"])],
        normalize=True,
        ylabel=r"$E_u(k) / E_u(3)$",
        ylim=(1.0e-4, 4.0),
    )
    axes[0, 2].legend(fontsize=8, loc="lower left")

    im2 = axes[1, 0].imshow(stream["slice_speed"], origin="lower", extent=stream["extent"],
                            cmap="magma", vmin=0.0, vmax=vmax)
    axes[1, 0].set_title("2D stream: speed")
    psi = stream["psi"]
    psimax = float(np.nanmax(np.abs(psi)))
    im3 = axes[1, 1].imshow(psi, origin="lower", extent=stream["extent"],
                            cmap="cividis",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-psimax, vmax=psimax))
    axes[1, 1].set_title(r"2D stream: $\psi$")
    psi_shell = stream["psi_shell_energy"]
    ax = axes[1, 2]
    ax.loglog(stream["kvals"][1:KMAX + 1], _shape_normalize(stream["shell_energy"])[1:KMAX + 1],
              lw=2.0, color=stream["case"]["color"], label=r"$E_u(k)/E_u(3)$")
    psi_ref = (stream["kvals"][1:KMAX + 1] / float(KREF)) ** (-(EXPO + 2.0))
    ax.loglog(stream["kvals"][1:KMAX + 1], _shape_normalize(psi_shell)[1:KMAX + 1],
              lw=2.0, color="#4c8c2b", label=r"$E_\psi(k)/E_\psi(3)$")
    ax.loglog(stream["kvals"][1:KMAX + 1], _reference_law(stream["kvals"][1:KMAX + 1]),
              color="black", ls="--", lw=1.2, label=r"$k^{-5/3}$")
    ax.loglog(stream["kvals"][1:KMAX + 1], psi_ref,
              color="#3b5d25", ls=":", lw=1.2, label=r"$k^{-11/3}$")
    ax.set_xlim(1.0, KMAX)
    ax.set_ylim(1.0e-4, 4.0)
    ax.set_xlabel(r"$k$")
    ax.set_ylabel("shape-normalized spectra")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(fontsize=8, loc="lower left")

    for ax in axes.flat:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.colorbar(im0, ax=[axes[0, 0], axes[1, 0]], shrink=0.95, location="right",
                 label=r"$|\mathbf{v}|$")
    fig.colorbar(im1, ax=axes[0, 1], shrink=0.9, location="right", label=r"$\omega_z$")
    fig.colorbar(im3, ax=axes[1, 1], shrink=0.9, location="right", label=r"$\psi$")
    fig.savefig(outdir / "methods_2d.png", dpi=220)
    plt.close(fig)


def _methods_3d_figure(results: dict[str, dict[str, object]], outdir: Path) -> None:
    proj = results["proj3d"]
    stream = results["stream3d"]
    fig, axes = plt.subplots(2, 4, figsize=(18, 9), constrained_layout=True)
    vmax = max(float(np.nanmax(proj["slice_speed"])), float(np.nanmax(stream["slice_speed"])))
    vzmax = float(np.nanmax(np.abs(proj["slice_vz"])))
    im0 = axes[0, 0].imshow(proj["slice_speed"], origin="lower", extent=proj["extent"],
                            cmap="magma", vmin=0.0, vmax=vmax)
    axes[0, 0].set_title("3D projection: speed slice")
    im1 = axes[0, 1].imshow(proj["slice_vz"], origin="lower", extent=proj["extent"],
                            cmap="coolwarm",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-vzmax, vmax=vzmax))
    axes[0, 1].set_title("3D projection: $v_z$")
    _plot_spectrum(
        axes[0, 2],
        proj["kvals"],
        [("velocity", proj["shell_energy"], proj["case"]["color"])],
        normalize=True,
        ylabel=r"$E_u(k) / E_u(3)$",
        ylim=(1.0e-4, 4.0),
    )
    axes[0, 2].legend(fontsize=8, loc="lower left")
    _shell_count_panel(axes[0, 3], proj["diag"], "velocity", KMAX, "projection mode counts")

    im2 = axes[1, 0].imshow(stream["slice_speed"], origin="lower", extent=stream["extent"],
                            cmap="magma", vmin=0.0, vmax=vmax)
    axes[1, 0].set_title("3D stream: speed slice")
    phi1 = stream["phi1"]
    phi2 = stream["phi2"]
    pmax = max(float(np.nanmax(np.abs(phi1))), float(np.nanmax(np.abs(phi2))))
    im3 = axes[1, 1].imshow(phi1, origin="lower", extent=stream["extent"], cmap="cividis",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-pmax, vmax=pmax))
    axes[1, 1].set_title(r"3D stream: $\phi_1$")
    im4 = axes[1, 2].imshow(phi2, origin="lower", extent=stream["extent"], cmap="cividis",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-pmax, vmax=pmax))
    axes[1, 2].set_title(r"3D stream: $\phi_2$")
    ax = axes[1, 3]
    ax.loglog(stream["kvals"][1:KMAX + 1], _shape_normalize(stream["shell_energy"])[1:KMAX + 1],
              lw=2.0, color=stream["case"]["color"], label=r"$E_u(k)/E_u(3)$")
    ax.loglog(stream["kvals"][1:KMAX + 1], _shape_normalize(stream["phi1_shell_energy"])[1:KMAX + 1],
              lw=2.0, color="#3b7a57", label=r"$E_{\phi_1}(k)/E_{\phi_1}(3)$")
    ax.loglog(stream["kvals"][1:KMAX + 1], _shape_normalize(stream["phi2_shell_energy"])[1:KMAX + 1],
              lw=2.0, color="#8f5a3c", label=r"$E_{\phi_2}(k)/E_{\phi_2}(3)$")
    ax.loglog(stream["kvals"][1:KMAX + 1], _reference_law(stream["kvals"][1:KMAX + 1]),
              color="black", ls="--", lw=1.2, label=r"$k^{-5/3}$")
    ax.set_xlim(1.0, KMAX)
    ax.set_ylim(1.0e-4, 4.0)
    ax.set_xlabel(r"$k$")
    ax.set_ylabel("shape-normalized spectra")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(fontsize=8, loc="lower left")

    for ax in axes.flat:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.colorbar(im0, ax=[axes[0, 0], axes[1, 0]], shrink=0.95, location="right",
                 label=r"$|\mathbf{v}|$")
    fig.colorbar(im1, ax=axes[0, 1], shrink=0.9, location="right", label=r"$v_z$")
    fig.colorbar(im3, ax=[axes[1, 1], axes[1, 2]], shrink=0.9, location="right", label=r"$\phi$")
    fig.savefig(outdir / "methods_3d.png", dpi=220)
    plt.close(fig)


def _stream3d_realization_figure(results: dict[str, dict[str, object]], outdir: Path) -> None:
    stream = results["stream3d"]
    diag = stream["diag"]
    kvals = stream["kvals"]
    predicted = np.asarray(diag["spectra"]["predicted_velocity_shell_energy"], dtype=np.float64)
    realized = np.asarray(diag["spectra"]["realized_retained_velocity_shell_energy"], dtype=np.float64)
    realized_raw = np.asarray(diag["spectra"]["realized_retained_velocity_shell_energy_raw"], dtype=np.float64)
    weights1 = np.asarray(diag["spectra"]["shell_weights_phi1"], dtype=np.float64)
    weights2 = np.asarray(diag["spectra"]["shell_weights_phi2"], dtype=np.float64)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    shell = np.arange(1, len(weights1), dtype=np.float64)
    axes[0, 0].loglog(shell, weights1[1:], lw=2.0, color="#3b7a57", label=r"$A_s^{(1)}$")
    axes[0, 0].loglog(shell, weights2[1:], lw=2.0, color="#8f5a3c", label=r"$A_s^{(2)}$")
    axes[0, 0].set_xlabel(r"$s$")
    axes[0, 0].set_ylabel("shell weight")
    axes[0, 0].set_title("3D stream shell weights")
    axes[0, 0].grid(True, which="both", alpha=0.2)
    axes[0, 0].legend(fontsize=9)

    ax = axes[0, 1]
    ax.loglog(kvals[1:KMAX + 1], _shape_normalize(predicted)[1:KMAX + 1], lw=2.0,
              color="#555555", label="predicted fit")
    ax.loglog(kvals[1:KMAX + 1], _shape_normalize(realized)[1:KMAX + 1], lw=2.0,
              color="#c45b12", label="realized retained")
    ax.loglog(kvals[1:KMAX + 1], _shape_normalize(realized_raw)[1:KMAX + 1], lw=1.5,
              color="#e3a069", label="realized raw")
    ax.loglog(kvals[1:KMAX + 1], _shape_normalize(stream["shell_energy"])[1:KMAX + 1], lw=2.0,
              color=stream["case"]["color"], label="FFT of initialized field")
    ax.loglog(kvals[1:KMAX + 1], _reference_law(kvals[1:KMAX + 1]), color="black",
              ls="--", lw=1.2, label=r"$k^{-5/3}$")
    ax.set_xlim(1.0, KMAX)
    ax.set_ylim(1.0e-4, 4.0)
    ax.set_xlabel(r"$k$")
    ax.set_ylabel("shape-normalized spectra")
    ax.set_title("Predicted and realized velocity spectra")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(fontsize=8, loc="lower left")

    _plot_sector_fractions(axes[1, 0], kvals, diag, "fitted_response",
                           "Sector fractions of fitted response")

    ax = axes[1, 1]
    compensated = stream["shell_energy"] * np.power(np.maximum(kvals, 1.0), EXPO)
    target_comp = compensated[KREF]
    if target_comp > 0.0:
        compensated = compensated / target_comp
    ratio = np.full_like(stream["shell_energy"], np.nan)
    target = _reference_law(kvals, slope=-EXPO, kref=KREF)
    model = _shape_normalize(stream["shell_energy"])
    ratio[1:KMAX + 1] = model[1:KMAX + 1] / target[1:KMAX + 1]
    ax.loglog(kvals[1:KMAX + 1], compensated[1:KMAX + 1], color=stream["case"]["color"],
              lw=2.0, label=r"$k^{5/3} E_u(k)$")
    ax.loglog(kvals[1:KMAX + 1], ratio[1:KMAX + 1], color="#444444", lw=1.5,
              label=r"$E_u / k^{-5/3}$")
    ax.axvspan(TAIL_KMIN, KMAX, color="#d0d8e2", alpha=0.25, label="tail diagnostic band")
    ax.set_xlim(1.0, KMAX)
    ax.set_ylim(1.0e-2, 1.0e2)
    ax.set_xlabel(r"$k$")
    ax.set_ylabel("compensated response")
    ax.set_title(
        f"Tail smoothness = {float(stream['tail_smoothness']):.2e}, "
        f"leakage = {float(diag['fit']['realized_out_of_band_leakage']):.2e}"
    )
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(fontsize=8, loc="lower left")
    fig.savefig(outdir / "stream3d_realization.png", dpi=220)
    plt.close(fig)


def _load_stream3d_audit_records(audit_json: Path) -> list[dict[str, object]]:
    data = json.loads(audit_json.read_text(encoding="utf-8"))
    records = [record for record in data.get("records", [])
               if record.get("case") == "prod_thinned"]
    if not records:
      raise RuntimeError(f"No prod_thinned records found in {audit_json}")
    return records


def _load_stream3d_audit_diag(record: dict[str, object]) -> dict[str, object]:
    run_dir = Path(str(record["run_dir"]))
    matches = sorted(run_dir.glob("*.turb_init_diag.json"))
    if not matches:
        raise RuntimeError(f"Missing diagnostics JSON in {run_dir}")
    return json.loads(matches[-1].read_text(encoding="utf-8"))


def _stream3d_production_audit_figure(audit_json: Path, outdir: Path) -> None:
    records = _load_stream3d_audit_records(audit_json)
    diags = [_load_stream3d_audit_diag(record) for record in records]

    det_slopes = np.array([float(record["deterministic_slope"]) for record in records],
                          dtype=np.float64)
    real_slopes = np.array([float(record["realized_slope"]) for record in records],
                           dtype=np.float64)
    seeds = [str(record["seed"]) for record in records]

    predicted = np.vstack([
        np.asarray(diag["spectra"]["predicted_velocity_shell_energy"], dtype=np.float64)
        for diag in diags
    ])
    realized = np.vstack([
        np.asarray(diag["spectra"]["realized_retained_velocity_shell_energy"], dtype=np.float64)
        for diag in diags
    ])

    kvals = np.arange(predicted.shape[1], dtype=np.float64)
    mask = (
        (kvals >= KMIN)
        & (kvals <= KMAX)
        & (np.mean(predicted, axis=0) > 0.0)
        & (np.mean(realized, axis=0) > 0.0)
    )
    kvals_band = kvals[mask]
    predicted_band = predicted[:, mask]
    realized_band = realized[:, mask]

    pred_mean = np.mean(predicted_band, axis=0)
    real_mean = np.mean(realized_band, axis=0)
    pred_min = np.min(predicted_band, axis=0)
    pred_max = np.max(predicted_band, axis=0)
    real_min = np.min(realized_band, axis=0)
    real_max = np.max(realized_band, axis=0)

    kref = 16.0
    ref_idx = int(np.argmin(np.abs(kvals_band - kref)))
    ref_amp = pred_mean[ref_idx]
    target = ref_amp * np.power(kvals_band / kvals_band[ref_idx], -EXPO)

    ratio_mean = real_mean / pred_mean
    ratio_lo = real_min / np.maximum(pred_max, 1.0e-300)
    ratio_hi = real_max / np.maximum(pred_min, 1.0e-300)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10, 8),
        gridspec_kw={"height_ratios": [3.0, 1.35]},
        constrained_layout=True,
        sharex=True,
    )

    ax = axes[0]
    for seed, curve in zip(seeds, predicted_band):
        ax.loglog(kvals_band, curve, color="#2C7FB8", alpha=0.22, lw=1.0)
    for seed, curve in zip(seeds, realized_band):
        ax.loglog(kvals_band, curve, color="#D95F0E", alpha=0.22, lw=1.0)
    ax.fill_between(kvals_band, pred_min, pred_max, color="#2C7FB8", alpha=0.10)
    ax.fill_between(kvals_band, real_min, real_max, color="#D95F0E", alpha=0.10)
    ax.loglog(
        kvals_band,
        pred_mean,
        color="#2C7FB8",
        lw=2.6,
        label=f"deterministic mean (slope {np.mean(det_slopes):.3f})",
    )
    ax.loglog(
        kvals_band,
        real_mean,
        color="#D95F0E",
        lw=2.6,
        label=f"realized mean (slope {np.mean(real_slopes):.3f})",
    )
    ax.loglog(kvals_band, target, color="0.25", ls="--", lw=1.4, label=r"target $k^{-5/3}$")
    ax.set_xlim(1.0, KMAX)
    ax.set_ylim(1.0e-5, 4.0e-1)
    ax.set_ylabel(r"$E_u(k)$")
    ax.set_title("Production 3D Clebsch audit at $256^3$: deterministic vs realized spectra")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(fontsize=9, loc="lower left")
    ax.text(
        0.02,
        0.04,
        "Thin lines: individual seeds\nBands: seed-to-seed range\nAll audited seeds pass the quality gate",
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.85", alpha=0.92),
    )

    ax = axes[1]
    ax.axhline(1.0, color="0.35", lw=1.2, ls="--")
    ax.fill_between(kvals_band, ratio_lo, ratio_hi, color="#C994C7", alpha=0.20)
    ax.semilogx(kvals_band, ratio_mean, color="#7A0177", lw=2.1)
    ax.set_xlim(1.0, KMAX)
    ax.set_ylim(0.45, 1.65)
    ax.set_xlabel(r"$k$")
    ax.set_ylabel("realized / deterministic")
    ax.grid(True, which="both", alpha=0.2)

    fig.savefig(outdir / "stream3d_production_audit.png", dpi=220)
    plt.close(fig)


def _atlas_projection_2d(results: dict[str, dict[str, object]], outdir: Path) -> None:
    item = results["proj2d"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    vmax = float(np.nanmax(item["slice_speed"]))
    im0 = axes[0, 0].imshow(item["slice_speed"], origin="lower", extent=item["extent"],
                            cmap="magma", vmin=0.0, vmax=vmax)
    axes[0, 0].set_title("speed")
    vxmax = float(np.nanmax(np.abs(item["slice_vx"])))
    im1 = axes[0, 1].imshow(item["slice_vx"], origin="lower", extent=item["extent"],
                            cmap="coolwarm",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-vxmax, vmax=vxmax))
    axes[0, 1].set_title(r"$v_x$")
    vymax = float(np.nanmax(np.abs(item["slice_vy"])))
    im2 = axes[0, 2].imshow(item["slice_vy"], origin="lower", extent=item["extent"],
                            cmap="coolwarm",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-vymax, vmax=vymax))
    axes[0, 2].set_title(r"$v_y$")

    axes[1, 0].loglog(item["kvals"][1:KMAX + 1], item["shell_energy"][1:KMAX + 1],
                      color=item["case"]["color"], lw=2.0)
    axes[1, 0].set_xlim(1.0, KMAX)
    axes[1, 0].set_ylim(1.0e-8, 1.0e0)
    axes[1, 0].set_xlabel(r"$k$")
    axes[1, 0].set_ylabel(r"$E_u(k)=\sum_{k\in {\rm shell}} |\hat{\mathbf{u}}|^2/N^2$")
    axes[1, 0].grid(True, which="both", alpha=0.2)
    comp = item["shell_energy"] * np.power(np.maximum(item["kvals"], 1.0), EXPO)
    comp /= comp[KREF]
    axes[1, 1].loglog(item["kvals"][1:KMAX + 1], comp[1:KMAX + 1], color=item["case"]["color"], lw=2.0)
    axes[1, 1].set_xlim(1.0, KMAX)
    axes[1, 1].set_ylim(1.0e-2, 1.0e2)
    axes[1, 1].set_xlabel(r"$k$")
    axes[1, 1].set_ylabel(r"$k^{5/3} E_u(k)$")
    axes[1, 1].grid(True, which="both", alpha=0.2)
    amp = _vector_mode_amplitude(item["diag"], "velocity")
    _occupancy_scatter(axes[1, 2], item["diag"], "velocity", amp, "retained velocity modes")

    for ax in axes.flat[:3]:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.colorbar(im0, ax=axes[0, 0], shrink=0.9, location="right", label=r"$|\mathbf{v}|$")
    fig.colorbar(im1, ax=axes[0, 1], shrink=0.9, location="right", label=r"$v_x$")
    fig.colorbar(im2, ax=axes[0, 2], shrink=0.9, location="right", label=r"$v_y$")
    fig.savefig(outdir / "atlas_2d_projection.png", dpi=220)
    plt.close(fig)


def _atlas_stream_2d(results: dict[str, dict[str, object]], outdir: Path) -> None:
    item = results["stream2d"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    vmax = float(np.nanmax(item["slice_speed"]))
    im0 = axes[0, 0].imshow(item["slice_speed"], origin="lower", extent=item["extent"],
                            cmap="magma", vmin=0.0, vmax=vmax)
    axes[0, 0].set_title("speed")
    wmax = float(np.nanmax(np.abs(item["vorticity"])))
    im1 = axes[0, 1].imshow(item["vorticity"], origin="lower", extent=item["extent"],
                            cmap="coolwarm",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-wmax, vmax=wmax))
    axes[0, 1].set_title(r"$\omega_z$")
    psi = item["psi"]
    psimax = float(np.nanmax(np.abs(psi)))
    im2 = axes[0, 2].imshow(psi, origin="lower", extent=item["extent"], cmap="cividis",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-psimax, vmax=psimax))
    axes[0, 2].set_title(r"$\psi$")
    axes[1, 0].loglog(item["kvals"][1:KMAX + 1], item["shell_energy"][1:KMAX + 1],
                      color=item["case"]["color"], lw=2.0)
    axes[1, 0].set_xlim(1.0, KMAX)
    axes[1, 0].set_ylim(1.0e-8, 1.0e0)
    axes[1, 0].set_xlabel(r"$k$")
    axes[1, 0].set_ylabel(r"$E_u(k)$")
    axes[1, 0].grid(True, which="both", alpha=0.2)
    axes[1, 1].loglog(item["kvals"][1:KMAX + 1], item["psi_shell_energy"][1:KMAX + 1],
                      color="#4c8c2b", lw=2.0, label=r"$E_\psi(k)$")
    guide = item["psi_shell_energy"][KREF] * (
        item["kvals"][1:KMAX + 1] / float(KREF)
    ) ** (-(EXPO + 2.0))
    axes[1, 1].loglog(item["kvals"][1:KMAX + 1], guide,
                      color="black", ls="--", lw=1.2, label=r"$k^{-11/3}$")
    axes[1, 1].set_xlim(1.0, KMAX)
    axes[1, 1].set_ylim(1.0e-12, 1.0e0)
    axes[1, 1].set_xlabel(r"$k$")
    axes[1, 1].set_ylabel(r"$E_\psi(k)$")
    axes[1, 1].grid(True, which="both", alpha=0.2)
    axes[1, 1].legend(fontsize=8, loc="lower left")
    amp = _scalar_mode_amplitude(item["diag"], "psi")
    _occupancy_scatter(axes[1, 2], item["diag"], "psi", amp, "retained psi modes")
    for ax in axes.flat[:3]:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.colorbar(im0, ax=axes[0, 0], shrink=0.9, location="right", label=r"$|\mathbf{v}|$")
    fig.colorbar(im1, ax=axes[0, 1], shrink=0.9, location="right", label=r"$\omega_z$")
    fig.colorbar(im2, ax=axes[0, 2], shrink=0.9, location="right", label=r"$\psi$")
    fig.savefig(outdir / "atlas_2d_stream.png", dpi=220)
    plt.close(fig)


def _atlas_projection_3d(results: dict[str, dict[str, object]], outdir: Path) -> None:
    item = results["proj3d"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    vmax = float(np.nanmax(item["slice_speed"]))
    im0 = axes[0, 0].imshow(item["slice_speed"], origin="lower", extent=item["extent"],
                            cmap="magma", vmin=0.0, vmax=vmax)
    axes[0, 0].set_title("speed slice")
    vxmax = float(np.nanmax(np.abs(item["slice_vx"])))
    im1 = axes[0, 1].imshow(item["slice_vx"], origin="lower", extent=item["extent"],
                            cmap="coolwarm",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-vxmax, vmax=vxmax))
    axes[0, 1].set_title(r"$v_x$")
    vzmax = float(np.nanmax(np.abs(item["slice_vz"])))
    im2 = axes[0, 2].imshow(item["slice_vz"], origin="lower", extent=item["extent"],
                            cmap="coolwarm",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-vzmax, vmax=vzmax))
    axes[0, 2].set_title(r"$v_z$")
    axes[1, 0].loglog(item["kvals"][1:KMAX + 1], item["shell_energy"][1:KMAX + 1],
                      color=item["case"]["color"], lw=2.0)
    axes[1, 0].set_xlim(1.0, KMAX)
    axes[1, 0].set_ylim(1.0e-9, 1.0e0)
    axes[1, 0].set_xlabel(r"$k$")
    axes[1, 0].set_ylabel(r"$E_u(k)$")
    axes[1, 0].grid(True, which="both", alpha=0.2)
    comp = item["shell_energy"] * np.power(np.maximum(item["kvals"], 1.0), EXPO)
    comp /= comp[KREF]
    axes[1, 1].loglog(item["kvals"][1:KMAX + 1], comp[1:KMAX + 1], color=item["case"]["color"], lw=2.0)
    axes[1, 1].set_xlim(1.0, KMAX)
    axes[1, 1].set_ylim(1.0e-2, 1.0e2)
    axes[1, 1].set_xlabel(r"$k$")
    axes[1, 1].set_ylabel(r"$k^{5/3} E_u(k)$")
    axes[1, 1].grid(True, which="both", alpha=0.2)
    _shell_count_panel(axes[1, 2], item["diag"], "velocity", KMAX, "retained velocity modes")
    for ax in axes.flat[:3]:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.colorbar(im0, ax=axes[0, 0], shrink=0.9, location="right", label=r"$|\mathbf{v}|$")
    fig.colorbar(im1, ax=axes[0, 1], shrink=0.9, location="right", label=r"$v_x$")
    fig.colorbar(im2, ax=axes[0, 2], shrink=0.9, location="right", label=r"$v_z$")
    fig.savefig(outdir / "atlas_3d_projection.png", dpi=220)
    plt.close(fig)


def _atlas_stream_3d(results: dict[str, dict[str, object]], outdir: Path) -> None:
    item = results["stream3d"]
    diag = item["diag"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    vmax = float(np.nanmax(item["slice_speed"]))
    im0 = axes[0, 0].imshow(item["slice_speed"], origin="lower", extent=item["extent"],
                            cmap="magma", vmin=0.0, vmax=vmax)
    axes[0, 0].set_title("speed slice")
    pmax = max(float(np.nanmax(np.abs(item["phi1"]))), float(np.nanmax(np.abs(item["phi2"]))))
    im1 = axes[0, 1].imshow(item["phi1"], origin="lower", extent=item["extent"], cmap="cividis",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-pmax, vmax=pmax))
    axes[0, 1].set_title(r"$\phi_1$")
    im2 = axes[0, 2].imshow(item["phi2"], origin="lower", extent=item["extent"], cmap="cividis",
                            norm=TwoSlopeNorm(vcenter=0.0, vmin=-pmax, vmax=pmax))
    axes[0, 2].set_title(r"$\phi_2$")
    axes[1, 0].loglog(item["kvals"][1:KMAX + 1], item["shell_energy"][1:KMAX + 1],
                      color=item["case"]["color"], lw=2.0)
    axes[1, 0].set_xlim(1.0, KMAX)
    axes[1, 0].set_ylim(1.0e-9, 1.0e0)
    axes[1, 0].set_xlabel(r"$k$")
    axes[1, 0].set_ylabel(r"$E_u(k)$")
    axes[1, 0].grid(True, which="both", alpha=0.2)
    _plot_sector_fractions(axes[1, 1], item["kvals"], diag, "fitted_response",
                           "Fitted sector fractions")
    weights1 = np.asarray(diag["spectra"]["shell_weights_phi1"], dtype=np.float64)
    weights2 = np.asarray(diag["spectra"]["shell_weights_phi2"], dtype=np.float64)
    shell = np.arange(1, len(weights1))
    axes[1, 2].loglog(shell, weights1[1:], lw=2.0, color="#3b7a57", label=r"$A_s^{(1)}$")
    axes[1, 2].loglog(shell, weights2[1:], lw=2.0, color="#8f5a3c", label=r"$A_s^{(2)}$")
    axes[1, 2].set_xlim(1.0, len(weights1) - 1)
    axes[1, 2].set_xlabel(r"$s$")
    axes[1, 2].set_ylabel("shell weight")
    axes[1, 2].set_title("Clebsch shell weights")
    axes[1, 2].grid(True, which="both", alpha=0.2)
    axes[1, 2].legend(fontsize=8)
    for ax in axes.flat[:3]:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.colorbar(im0, ax=axes[0, 0], shrink=0.9, location="right", label=r"$|\mathbf{v}|$")
    fig.colorbar(im1, ax=[axes[0, 1], axes[0, 2]], shrink=0.9, location="right", label=r"$\phi$")
    fig.savefig(outdir / "atlas_3d_stream.png", dpi=220)
    plt.close(fig)


def _write_summary(results: dict[str, dict[str, object]], outdir: Path) -> None:
    summary = {}
    for slug, item in results.items():
        entry = {
            "name": item["case"]["name"],
            "init_seconds": float(item["init_seconds"]),
            "v_rms": float(item["v_rms"]),
            "target_v_rms": float(item["meta"]["target"]["v_rms"]),
            "v_rms_rel_error": abs(float(item["v_rms"]) - float(item["meta"]["target"]["v_rms"])) /
            float(item["meta"]["target"]["v_rms"]),
            "divergence_ratio": float(item["divergence_ratio"]),
            "slope": float(item["slope"]),
            "velocity_diag_match_rms": float(item["velocity_diag_match_rms"]),
        }
        if slug == "stream2d":
            entry["psi_reconstruct_max_abs"] = float(item["psi_reconstruct_max_abs"])
        if slug == "stream3d":
            entry["phi1_reconstruct_max_abs"] = float(item["phi1_reconstruct_max_abs"])
            entry["phi2_reconstruct_max_abs"] = float(item["phi2_reconstruct_max_abs"])
            entry["tail_smoothness"] = float(item["tail_smoothness"])
            entry["fit"] = item["diag"]["fit"]
        summary[slug] = entry
    (outdir / "campaign_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--athena", type=Path, default=_default_athena_path())
    parser.add_argument("--force", action="store_true",
                        help="rerun Athena even if outputs already exist")
    parser.add_argument(
        "--clebsch-audit-json",
        type=Path,
        default=_default_clebsch_audit_json(),
        help="optional production 3D Clebsch audit JSON used to build the multi-seed audit figure",
    )
    args = parser.parse_args()

    outdir = _figure_root()
    outdir.mkdir(parents=True, exist_ok=True)

    results: dict[str, dict[str, object]] = {}
    athena = args.athena.resolve()
    for case in _case_specs():
        print(f"[run] {case['name']}")
        hydro_path, diag_path, init_seconds = _run_case(athena, case, args.force)
        print(f"[analyze] {case['name']}")
        results[str(case["slug"])] = _extract_case(hydro_path, diag_path, case, init_seconds)

    print("[figure] overview")
    _overview_figure(results, outdir)
    print("[figure] methods_2d")
    _methods_2d_figure(results, outdir)
    print("[figure] methods_3d")
    _methods_3d_figure(results, outdir)
    print("[figure] stream3d_realization")
    _stream3d_realization_figure(results, outdir)
    print("[figure] atlas_2d_projection")
    _atlas_projection_2d(results, outdir)
    print("[figure] atlas_2d_stream")
    _atlas_stream_2d(results, outdir)
    print("[figure] atlas_3d_projection")
    _atlas_projection_3d(results, outdir)
    print("[figure] atlas_3d_stream")
    _atlas_stream_3d(results, outdir)
    if args.clebsch_audit_json is not None:
        print("[figure] stream3d_production_audit")
        _stream3d_production_audit_figure(args.clebsch_audit_json.resolve(), outdir)
    print("[write] campaign_summary.json")
    _write_summary(results, outdir)


if __name__ == "__main__":
    main()
