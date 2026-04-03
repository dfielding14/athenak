#!/usr/bin/env python3
"""Generate scalar-mixing velocity-method figures using the simplified Clebsch diagnostics."""

from __future__ import annotations

import argparse
import glob
import json
import math
import re
import shlex
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bin_convert_new import read_binary_as_athdf


EXPO = 5.0 / 3.0
ALPHA = 1.0 / 3.0
KMIN = 2
KMAX = 128
KREF = 3


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
    return _repo_root() / "build" / "src" / "athena"


def _figure_root() -> Path:
    return _repo_root() / "docs" / "figures" / "scalar_mixing_velocity_methods"


def _default_clebsch_audit_json() -> Path | None:
    for relpath in (
        Path("tmp/clebsch_audit_report_prod256.json"),
        Path("tmp/clebsch_audit_report.json"),
    ):
        candidate = _repo_root() / relpath
        if candidate.exists():
            return candidate
    return None


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
            "overrides": ["problem/turb_velocity_method=stream_2d"],
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
            "name": "3D Clebsch",
            "slug": "clebsch3d",
            "input": "inputs/hydro/scalar_mixing_methods_paper_3d.athinput",
            "basename": "scalar_mix_methods_clebsch3d",
            "overrides": [
                "problem/turb_velocity_method=clebsch",
                f"problem/turb_alpha={ALPHA:.16f}",
            ],
            "launcher": "mpirun -np 8",
            "dim": "3d",
            "method": "clebsch",
            "color": "#8b3fb3",
        },
    ]


def _parse_init_seconds(log_text: str) -> float | None:
    match = re.search(r"velocity_init_wall_seconds=([0-9.eE+-]+)", log_text)
    return None if match is None else float(match.group(1))


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
    init_seconds = _parse_init_seconds(completed.stdout + completed.stderr)
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
    return field[0] if field.shape[0] == 1 else field[field.shape[0] // 2]


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
    return 0.0 if denom <= 0.0 else float(np.sqrt(numer / denom))


def _fit_slope(kvals: np.ndarray, spectrum: np.ndarray, kmin: int, kmax: int) -> float:
    mask = (kvals >= kmin) & (kvals <= kmax) & (spectrum > 0.0)
    if np.count_nonzero(mask) < 2:
        return float("nan")
    coeff = np.polyfit(np.log(kvals[mask]), np.log(spectrum[mask]), 1)
    return float(coeff[0])


def _shape_normalize(spectrum: np.ndarray, kref: int = KREF) -> np.ndarray:
    normalized = np.full_like(spectrum, np.nan, dtype=np.float64)
    if kref >= spectrum.size or spectrum[kref] <= 0.0:
        return normalized
    normalized[:] = spectrum / spectrum[kref]
    return normalized


def _reference_law(kvals: np.ndarray, slope: float, kref: int = KREF) -> np.ndarray:
    ref = np.full_like(kvals, np.nan, dtype=np.float64)
    mask = kvals > 0.0
    ref[mask] = (kvals[mask] / float(kref)) ** slope
    return ref


def _extent(meta: dict[str, object]) -> list[float]:
    mesh = meta["mesh"]
    return [
        float(mesh["x1min"]),
        float(mesh["x1max"]),
        float(mesh["x2min"]),
        float(mesh["x2max"]),
    ]


def _vorticity_2d(vx: np.ndarray, vy: np.ndarray, meta: dict[str, object]) -> np.ndarray:
    ny, nx = vx.shape
    mesh = meta["mesh"]
    kx = _fft_wave_numbers(nx, float(mesh["x1min"]), float(mesh["x1max"]))
    ky = _fft_wave_numbers(ny, float(mesh["x2min"]), float(mesh["x2max"]))
    ky2, kx2 = np.meshgrid(ky, kx, indexing="ij")
    omega_hat = 1j * kx2 * np.fft.fft2(vy) - 1j * ky2 * np.fft.fft2(vx)
    return np.fft.ifft2(omega_hat).real


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
    dkz = 0.0 if int(mesh["global_nx3"]) == 1 else 2.0 * np.pi / (
        float(mesh["x3max"]) - float(mesh["x3min"])
    )
    nkx = np.asarray(catalog["nkx"], dtype=np.float64)
    nky = np.asarray(catalog["nky"], dtype=np.float64)
    nkz = np.asarray(catalog["nkz"], dtype=np.float64)
    aka = np.asarray(coeffs["aka"], dtype=np.float64)
    akb = np.asarray(coeffs["akb"], dtype=np.float64)
    x = np.linspace(float(mesh["x1min"]), float(mesh["x1max"]), nx, endpoint=False)
    y = np.linspace(float(mesh["x2min"]), float(mesh["x2max"]), ny, endpoint=False)
    x += 0.5 * (float(mesh["x1max"]) - float(mesh["x1min"])) / nx
    y += 0.5 * (float(mesh["x2max"]) - float(mesh["x2min"])) / ny
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


def _catalog_shell_counts(diag: dict[str, object], catalog_key: str, nmax: int) -> np.ndarray:
    shells = np.asarray(diag["catalogs"][catalog_key]["shell"], dtype=np.int64)
    counts = np.bincount(shells, minlength=nmax + 1)
    return counts[:nmax + 1]


def _extract_case(velocity_path: Path, diag_path: Path, case: dict[str, object],
                  init_seconds: float | None) -> dict[str, object]:
    diag = json.loads(diag_path.read_text(encoding="utf-8"))
    velocity = _load_velocity(velocity_path)
    vx = velocity["velx"]
    vy = velocity["vely"]
    vz = velocity["velz"]
    kvals, shell_energy = _shell_spectrum(vx, vy, vz)
    meta = diag["metadata"]
    diag_velocity_shell = None
    if "velocity_shell_energy" in diag["spectra"]:
        scale = float(meta["normalization"]["scale"])
        diag_velocity_shell = (
            np.asarray(diag["spectra"]["velocity_shell_energy"], dtype=np.float64) * (scale * scale)
        )

    item: dict[str, object] = {
        "case": case,
        "diag": diag,
        "meta": meta,
        "extent": _extent(meta),
        "slice_vx": _midplane(vx),
        "slice_vy": _midplane(vy),
        "slice_vz": _midplane(vz),
        "kvals": kvals,
        "shell_energy": shell_energy,
        "shell_shape": _shape_normalize(shell_energy),
        "diag_velocity_shell": diag_velocity_shell,
        "diag_velocity_shape": (
            _shape_normalize(diag_velocity_shell)
            if diag_velocity_shell is not None else None
        ),
        "divergence_ratio": _spectral_divergence_ratio(vx, vy, vz, meta),
        "v_rms": float(np.sqrt(np.mean(vx * vx + vy * vy + vz * vz))),
        "slope": _fit_slope(kvals, shell_energy, KMIN, min(KMAX, len(shell_energy) - 1)),
        "init_seconds": init_seconds if init_seconds is not None else meta.get("init_wall_seconds"),
    }
    item["slice_speed"] = np.sqrt(
        item["slice_vx"] * item["slice_vx"]
        + item["slice_vy"] * item["slice_vy"]
        + item["slice_vz"] * item["slice_vz"]
    )

    if str(case["slug"]) == "stream2d":
        psi = np.asarray(diag["slices"]["psi_xy"]["data"], dtype=np.float64).reshape(
            diag["slices"]["psi_xy"]["shape"]
        )
        recon = _reconstruct_slice(diag, "psi", "psi", "psi_xy")
        item["psi"] = psi
        item["psi_shell_energy"] = np.asarray(diag["spectra"]["psi_shell_energy"], dtype=np.float64)
        item["psi_reconstruct_max_abs"] = float(np.max(np.abs(psi - recon)))
    elif str(case["slug"]) == "clebsch3d":
        phi1 = np.asarray(diag["slices"]["phi1_xy_mid"]["data"], dtype=np.float64).reshape(
            diag["slices"]["phi1_xy_mid"]["shape"]
        )
        phi2 = np.asarray(diag["slices"]["phi2_xy_mid"]["data"], dtype=np.float64).reshape(
            diag["slices"]["phi2_xy_mid"]["shape"]
        )
        recon1 = _reconstruct_slice(diag, "phi", "phi1", "phi1_xy_mid")
        recon2 = _reconstruct_slice(diag, "phi", "phi2", "phi2_xy_mid")
        item["phi1"] = phi1
        item["phi2"] = phi2
        item["phi1_shell_energy"] = np.asarray(diag["spectra"]["phi1_shell_energy"], dtype=np.float64)
        item["phi2_shell_energy"] = np.asarray(diag["spectra"]["phi2_shell_energy"], dtype=np.float64)
        item["phi1_reconstruct_max_abs"] = float(np.max(np.abs(phi1 - recon1)))
        item["phi2_reconstruct_max_abs"] = float(np.max(np.abs(phi2 - recon2)))

    if str(case["dim"]) == "2d":
        item["vorticity"] = _vorticity_2d(item["slice_vx"], item["slice_vy"], meta)

    return item


def _image_limits(field: np.ndarray, symmetric: bool = True) -> tuple[float, float]:
    if symmetric:
        vmax = float(np.nanpercentile(np.abs(field), 99.0))
        return -vmax, vmax
    return 0.0, float(np.nanpercentile(field, 99.0))


def _plot_image(ax: plt.Axes, field: np.ndarray, extent: list[float], title: str,
                cmap: str, symmetric: bool = True) -> None:
    vmin, vmax = _image_limits(field, symmetric=symmetric)
    image = ax.imshow(field, origin="lower", extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax._image = image  # type: ignore[attr-defined]


def _plot_velocity_spectrum(ax: plt.Axes, item: dict[str, object], title: str,
                            color: str, kmax: int = KMAX) -> None:
    kvals = np.asarray(item["kvals"], dtype=np.float64)
    fft_shape = np.asarray(item["shell_shape"], dtype=np.float64)
    diag_shape = item["diag_velocity_shape"]
    ref = _reference_law(kvals, -float(item["meta"]["target"]["expo"]))
    stop = min(kmax, len(kvals) - 1)
    ax.loglog(kvals[1:stop + 1], ref[1:stop + 1], "k--", lw=1.5, label=r"target")
    ax.loglog(kvals[1:stop + 1], fft_shape[1:stop + 1], color=color, lw=2.0, label="FFT")
    if diag_shape is not None:
        diag_shape = np.asarray(diag_shape, dtype=np.float64)
        ax.loglog(kvals[1:stop + 1], diag_shape[1:stop + 1], color=color, lw=1.5,
                  ls=":", label="diag")
    ax.set_title(title)
    ax.set_xlabel("k")
    ax.set_ylabel(r"$E(k)/E(3)$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)


def _plot_scalar_spectrum(ax: plt.Axes, shell_energy: np.ndarray, slope: float,
                          label: str, color: str) -> None:
    kvals = np.arange(shell_energy.size, dtype=np.float64)
    shape = _shape_normalize(shell_energy)
    ref = _reference_law(kvals, -slope)
    stop = min(KMAX, len(kvals) - 1)
    ax.loglog(kvals[1:stop + 1], ref[1:stop + 1], "k--", lw=1.0, label="target")
    ax.loglog(kvals[1:stop + 1], shape[1:stop + 1], color=color, lw=2.0, label=label)
    ax.set_xlabel("k")
    ax.set_ylabel(r"$E(k)/E(3)$")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)


def _plot_shell_counts(ax: plt.Axes, counts: np.ndarray, title: str, color: str) -> None:
    kvals = np.arange(counts.size, dtype=np.float64)
    ax.bar(kvals[1:], counts[1:], width=0.9, color=color, alpha=0.8)
    ax.set_xlim(0.5, min(KMAX, len(kvals) - 1) + 0.5)
    ax.set_xlabel("shell")
    ax.set_ylabel("kept modes")
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.25)


def _overview_figure(results: dict[str, dict[str, object]], outdir: Path) -> None:
    ordered = [results[key] for key in ("proj2d", "stream2d", "proj3d", "clebsch3d")]
    vmax = max(float(np.nanmax(item["slice_speed"])) for item in ordered)
    fig = plt.figure(figsize=(16, 10), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.2])
    slice_axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
    ]
    spec_ax = fig.add_subplot(grid[:, 2])

    for ax, item in zip(slice_axes, ordered):
        field = np.asarray(item["slice_speed"], dtype=np.float64)
        image = ax.imshow(field, origin="lower", extent=item["extent"], cmap="magma",
                          vmin=0.0, vmax=vmax)
        ax.set_title(
            f"{item['case']['name']}\n"
            f"t_init={float(item['init_seconds']):.2f}s, "
            f"div={float(item['divergence_ratio']):.2e}"
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal")
        ax._image = image  # type: ignore[attr-defined]

    kvals = np.arange(1, KMAX + 1, dtype=np.float64)
    ref = (kvals / KREF) ** (-EXPO)
    spec_ax.loglog(kvals, ref, "k--", lw=2.0, label=r"target $k^{-5/3}$")
    for slug in ("proj2d", "stream2d", "proj3d", "clebsch3d"):
        item = results[slug]
        color = str(item["case"]["color"])
        shell_shape = np.asarray(item["shell_shape"], dtype=np.float64)
        stop = min(KMAX, len(shell_shape) - 1)
        spec_ax.loglog(np.arange(1, stop + 1), shell_shape[1:stop + 1],
                       lw=2.0, color=color, label=str(item["case"]["name"]))
    spec_ax.set_xlabel("k")
    spec_ax.set_ylabel(r"$E(k)/E(3)$")
    spec_ax.set_title("Velocity Spectra")
    spec_ax.grid(True, which="both", alpha=0.25)
    spec_ax.legend(frameon=False, loc="lower left")

    colorbar = fig.colorbar(slice_axes[0]._image, ax=slice_axes, fraction=0.02, pad=0.02)
    colorbar.set_label(r"$|\mathbf{v}|$")
    fig.savefig(outdir / "overview.png", dpi=220)
    plt.close(fig)


def _methods_2d_figure(results: dict[str, dict[str, object]], outdir: Path) -> None:
    proj = results["proj2d"]
    stream = results["stream2d"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)

    _plot_image(axes[0, 0], np.asarray(proj["slice_speed"]), proj["extent"],
                "2D Projection Speed", "magma", symmetric=False)
    _plot_image(axes[0, 1], np.asarray(proj["vorticity"]), proj["extent"],
                "2D Projection Vorticity", "coolwarm")
    _plot_velocity_spectrum(axes[0, 2], proj, "2D Projection Velocity Spectrum",
                            str(proj["case"]["color"]))

    _plot_image(axes[1, 0], np.asarray(stream["slice_speed"]), stream["extent"],
                "2D Stream Speed", "magma", symmetric=False)
    _plot_image(axes[1, 1], np.asarray(stream["psi"]), stream["extent"],
                r"2D Stream Function $\psi$", "coolwarm")
    _plot_velocity_spectrum(axes[1, 2], stream, "2D Stream Velocity Spectrum",
                            str(stream["case"]["color"]))
    _plot_scalar_spectrum(axes[1, 2], np.asarray(stream["psi_shell_energy"], dtype=np.float64),
                          EXPO + 2.0, r"$E_\psi(k)$", "#6a3d9a")
    fig.savefig(outdir / "methods_2d.png", dpi=220)
    plt.close(fig)


def _methods_3d_figure(results: dict[str, dict[str, object]], outdir: Path) -> None:
    proj = results["proj3d"]
    clebsch = results["clebsch3d"]
    fig, axes = plt.subplots(2, 4, figsize=(18, 9), constrained_layout=True)

    _plot_image(axes[0, 0], np.asarray(proj["slice_speed"]), proj["extent"],
                "3D Projection Speed", "magma", symmetric=False)
    _plot_image(axes[0, 1], np.asarray(proj["slice_vz"]), proj["extent"],
                r"3D Projection $v_z$", "coolwarm")
    _plot_velocity_spectrum(axes[0, 2], proj, "3D Projection Velocity Spectrum",
                            str(proj["case"]["color"]))
    _plot_shell_counts(
        axes[0, 3],
        _catalog_shell_counts(proj["diag"], "velocity", min(KMAX, len(proj["kvals"]) - 1)),
        "3D Projection Retained Modes",
        str(proj["case"]["color"]),
    )

    _plot_image(axes[1, 0], np.asarray(clebsch["slice_speed"]), clebsch["extent"],
                "3D Clebsch Speed", "magma", symmetric=False)
    _plot_image(axes[1, 1], np.asarray(clebsch["phi1"]), clebsch["extent"],
                r"3D Clebsch $\phi_1$", "coolwarm")
    _plot_image(axes[1, 2], np.asarray(clebsch["phi2"]), clebsch["extent"],
                r"3D Clebsch $\phi_2$", "coolwarm")
    _plot_velocity_spectrum(axes[1, 3], clebsch, "3D Clebsch Velocity Spectrum",
                            str(clebsch["case"]["color"]))
    _plot_scalar_spectrum(
        axes[1, 3],
        np.asarray(clebsch["phi1_shell_energy"], dtype=np.float64),
        float(clebsch["meta"]["target"]["phi_slope"]),
        r"$E_{\phi_1}(k)$",
        "#c04c4c",
    )
    _plot_scalar_spectrum(
        axes[1, 3],
        np.asarray(clebsch["phi2_shell_energy"], dtype=np.float64),
        float(clebsch["meta"]["target"]["phi_slope"]),
        r"$E_{\phi_2}(k)$",
        "#3b8f6b",
    )
    fig.savefig(outdir / "methods_3d.png", dpi=220)
    plt.close(fig)


def _clebsch_realization_figure(results: dict[str, dict[str, object]], outdir: Path) -> None:
    item = results["clebsch3d"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    _plot_image(axes[0, 0], np.asarray(item["phi1"]), item["extent"],
                r"$\phi_1$ midplane", "coolwarm")
    _plot_image(axes[0, 1], np.asarray(item["phi2"]), item["extent"],
                r"$\phi_2$ midplane", "coolwarm")

    ax = axes[1, 0]
    _plot_scalar_spectrum(
        ax, np.asarray(item["phi1_shell_energy"], dtype=np.float64),
        float(item["meta"]["target"]["phi_slope"]), r"$E_{\phi_1}(k)$", "#c04c4c"
    )
    _plot_scalar_spectrum(
        ax, np.asarray(item["phi2_shell_energy"], dtype=np.float64),
        float(item["meta"]["target"]["phi_slope"]), r"$E_{\phi_2}(k)$", "#3b8f6b"
    )
    ax.set_title("Clebsch Scalar Shell Spectra")

    ax = axes[1, 1]
    kvals = np.asarray(item["kvals"], dtype=np.float64)
    fft_shape = np.asarray(item["shell_shape"], dtype=np.float64)
    ref = _reference_law(kvals, -float(item["meta"]["target"]["velocity_slope"]))
    stop = min(KMAX, len(kvals) - 1)
    ax.loglog(kvals[1:stop + 1], ref[1:stop + 1], "k--", lw=1.5, label=r"target")
    ax.loglog(kvals[1:stop + 1], fft_shape[1:stop + 1], color="#275dad", lw=2.0,
              label="FFT")
    ax.set_xlabel("k")
    ax.set_ylabel(r"$E(k)/E(3)$")
    ax.set_title("Clebsch Velocity FFT Spectrum")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)

    fig.savefig(outdir / "clebsch3d_realization.png", dpi=220)
    plt.close(fig)


def _clebsch_seed_audit_figure(audit_json: Path, outdir: Path) -> None:
    report = json.loads(audit_json.read_text(encoding="utf-8"))
    records = report["records"]
    summary = report["summary"]
    fig, axes = plt.subplots(2, 1, figsize=(11, 10), constrained_layout=True)

    top = axes[0]
    color_map = {
        "small_unthinned": "#275dad",
        "small_thinned": "#bd632f",
        "prod_thinned": "#8b3fb3",
    }
    for record in records:
        spectrum = np.asarray(record["velocity_shell_energy"], dtype=np.float64)
        shape = _shape_normalize(spectrum)
        stop = min(int(record["analysis_nhigh"]), len(shape) - 1)
        kvals = np.arange(1, stop + 1, dtype=np.float64)
        top.loglog(kvals, shape[1:stop + 1], alpha=0.35,
                   color=color_map.get(str(record["case"]), "#444444"))
    for case_name, case_summary in summary.items():
        case_records = [r for r in records if str(r["case"]) == case_name]
        spectra = np.array([_shape_normalize(np.asarray(r["velocity_shell_energy"], dtype=np.float64))
                            for r in case_records])
        mean_shape = np.nanmean(spectra, axis=0)
        stop = min(int(case_summary["analysis_nhigh"]), len(mean_shape) - 1)
        kvals = np.arange(1, stop + 1, dtype=np.float64)
        top.loglog(kvals, mean_shape[1:stop + 1], lw=2.5,
                   color=color_map.get(case_name, "#444444"), label=case_name)
    ref_k = np.arange(1, KMAX + 1, dtype=np.float64)
    top.loglog(ref_k, (ref_k / KREF) ** (-EXPO), "k--", lw=1.5, label=r"target $k^{-5/3}$")
    top.set_xlabel("k")
    top.set_ylabel(r"$E(k)/E(3)$")
    top.set_title("Multi-seed Clebsch Velocity Spectra")
    top.grid(True, which="both", alpha=0.25)
    top.legend(frameon=False)

    bottom = axes[1]
    case_names = list(summary.keys())
    xpos = np.arange(len(case_names), dtype=np.float64)
    slopes = [-float(summary[name]["mean_retained_velocity_slope"]) for name in case_names]
    leakages = [float(summary[name]["mean_retained_velocity_leakage"]) for name in case_names]
    width = 0.38
    bottom.bar(xpos - 0.5 * width, slopes, width=width, color="#5c8c2b",
               label=r"mean slope magnitude")
    bottom.bar(xpos + 0.5 * width, leakages, width=width, color="#b24a2f",
               label="mean leakage")
    bottom.axhline(EXPO, color="k", linestyle="--", linewidth=1.2, label=r"target $5/3$")
    bottom.set_xticks(xpos, case_names)
    bottom.set_ylabel("diagnostic value")
    bottom.set_title("Audit Summary")
    bottom.grid(True, axis="y", alpha=0.25)
    bottom.legend(frameon=False)

    fig.savefig(outdir / "clebsch3d_seed_audit.png", dpi=220)
    plt.close(fig)


def _atlas_projection_2d(results: dict[str, dict[str, object]], outdir: Path) -> None:
    item = results["proj2d"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    _plot_image(axes[0, 0], np.asarray(item["slice_speed"]), item["extent"], "speed", "magma",
                symmetric=False)
    _plot_image(axes[0, 1], np.asarray(item["slice_vx"]), item["extent"], r"$v_x$", "coolwarm")
    _plot_image(axes[0, 2], np.asarray(item["slice_vy"]), item["extent"], r"$v_y$", "coolwarm")
    _plot_velocity_spectrum(axes[1, 0], item, "velocity spectrum", str(item["case"]["color"]))
    compensated = np.asarray(item["shell_shape"], dtype=np.float64) * (
        np.asarray(item["kvals"], dtype=np.float64) / KREF
    ) ** EXPO
    axes[1, 1].plot(item["kvals"][1:KMAX + 1], compensated[1:KMAX + 1], color=str(item["case"]["color"]))
    axes[1, 1].axhline(1.0, color="k", linestyle="--", linewidth=1.0)
    axes[1, 1].set_xlabel("k")
    axes[1, 1].set_ylabel(r"$E(k)k^{5/3}/E(3)$")
    axes[1, 1].set_title("compensated spectrum")
    axes[1, 1].grid(True, alpha=0.25)
    _plot_shell_counts(
        axes[1, 2],
        _catalog_shell_counts(item["diag"], "velocity", KMAX),
        "retained mode counts",
        str(item["case"]["color"]),
    )
    fig.savefig(outdir / "atlas_2d_projection.png", dpi=220)
    plt.close(fig)


def _atlas_stream_2d(results: dict[str, dict[str, object]], outdir: Path) -> None:
    item = results["stream2d"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    _plot_image(axes[0, 0], np.asarray(item["slice_speed"]), item["extent"], "speed", "magma",
                symmetric=False)
    _plot_image(axes[0, 1], np.asarray(item["psi"]), item["extent"], r"$\psi$", "coolwarm")
    _plot_image(axes[0, 2], np.asarray(item["vorticity"]), item["extent"], "vorticity",
                "coolwarm")
    _plot_velocity_spectrum(axes[1, 0], item, "velocity spectrum", str(item["case"]["color"]))
    ax = axes[1, 1]
    _plot_scalar_spectrum(ax, np.asarray(item["psi_shell_energy"], dtype=np.float64),
                          EXPO + 2.0, r"$E_\psi(k)$", "#6a3d9a")
    ax.set_title("stream-function spectrum")
    _plot_shell_counts(
        axes[1, 2],
        _catalog_shell_counts(item["diag"], "psi", KMAX),
        "retained mode counts",
        str(item["case"]["color"]),
    )
    fig.savefig(outdir / "atlas_2d_stream.png", dpi=220)
    plt.close(fig)


def _atlas_projection_3d(results: dict[str, dict[str, object]], outdir: Path) -> None:
    item = results["proj3d"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    _plot_image(axes[0, 0], np.asarray(item["slice_speed"]), item["extent"], "speed", "magma",
                symmetric=False)
    _plot_image(axes[0, 1], np.asarray(item["slice_vx"]), item["extent"], r"$v_x$", "coolwarm")
    _plot_image(axes[0, 2], np.asarray(item["slice_vz"]), item["extent"], r"$v_z$", "coolwarm")
    _plot_velocity_spectrum(axes[1, 0], item, "velocity spectrum", str(item["case"]["color"]))
    compensated = np.asarray(item["shell_shape"], dtype=np.float64) * (
        np.asarray(item["kvals"], dtype=np.float64) / KREF
    ) ** EXPO
    axes[1, 1].plot(item["kvals"][1:KMAX + 1], compensated[1:KMAX + 1], color=str(item["case"]["color"]))
    axes[1, 1].axhline(1.0, color="k", linestyle="--", linewidth=1.0)
    axes[1, 1].set_xlabel("k")
    axes[1, 1].set_ylabel(r"$E(k)k^{5/3}/E(3)$")
    axes[1, 1].set_title("compensated spectrum")
    axes[1, 1].grid(True, alpha=0.25)
    _plot_shell_counts(
        axes[1, 2],
        _catalog_shell_counts(item["diag"], "velocity", KMAX),
        "retained mode counts",
        str(item["case"]["color"]),
    )
    fig.savefig(outdir / "atlas_3d_projection.png", dpi=220)
    plt.close(fig)


def _atlas_clebsch_3d(results: dict[str, dict[str, object]], outdir: Path) -> None:
    item = results["clebsch3d"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    _plot_image(axes[0, 0], np.asarray(item["slice_speed"]), item["extent"], "speed", "magma",
                symmetric=False)
    _plot_image(axes[0, 1], np.asarray(item["phi1"]), item["extent"], r"$\phi_1$", "coolwarm")
    _plot_image(axes[0, 2], np.asarray(item["phi2"]), item["extent"], r"$\phi_2$", "coolwarm")
    _plot_velocity_spectrum(axes[1, 0], item, "velocity spectrum", str(item["case"]["color"]))
    ax = axes[1, 1]
    _plot_scalar_spectrum(ax, np.asarray(item["phi1_shell_energy"], dtype=np.float64),
                          float(item["meta"]["target"]["phi_slope"]),
                          r"$E_{\phi_1}(k)$", "#c04c4c")
    _plot_scalar_spectrum(ax, np.asarray(item["phi2_shell_energy"], dtype=np.float64),
                          float(item["meta"]["target"]["phi_slope"]),
                          r"$E_{\phi_2}(k)$", "#3b8f6b")
    ax.set_title("Clebsch scalar spectra")
    _plot_shell_counts(
        axes[1, 2],
        _catalog_shell_counts(item["diag"], "phi", KMAX),
        "retained scalar modes",
        str(item["case"]["color"]),
    )
    fig.savefig(outdir / "atlas_3d_clebsch.png", dpi=220)
    plt.close(fig)


def _write_summary(results: dict[str, dict[str, object]], outdir: Path) -> None:
    summary = {}
    for slug, item in results.items():
        entry = {
            "name": item["case"]["name"],
            "init_seconds": float(item["init_seconds"]),
            "v_rms": float(item["v_rms"]),
            "target_v_rms": float(item["meta"]["target"]["v_rms"]),
            "divergence_ratio": float(item["divergence_ratio"]),
            "fft_shell_slope": float(item["slope"]),
            "velocity_method": item["meta"]["velocity_method"],
        }
        if slug == "stream2d":
            entry["psi_reconstruct_max_abs"] = float(item["psi_reconstruct_max_abs"])
        if slug == "clebsch3d":
            entry["alpha"] = float(item["meta"]["target"]["alpha"])
            entry["phi_slope"] = float(item["meta"]["target"]["phi_slope"])
            entry["velocity_slope"] = float(item["meta"]["target"]["velocity_slope"])
            entry["phi1_reconstruct_max_abs"] = float(item["phi1_reconstruct_max_abs"])
            entry["phi2_reconstruct_max_abs"] = float(item["phi2_reconstruct_max_abs"])
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
        help="optional Clebsch audit JSON used to build the multi-seed audit figure",
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
    print("[figure] clebsch3d_realization")
    _clebsch_realization_figure(results, outdir)
    print("[figure] atlas_2d_projection")
    _atlas_projection_2d(results, outdir)
    print("[figure] atlas_2d_stream")
    _atlas_stream_2d(results, outdir)
    print("[figure] atlas_3d_projection")
    _atlas_projection_3d(results, outdir)
    print("[figure] atlas_3d_clebsch")
    _atlas_clebsch_3d(results, outdir)
    if args.clebsch_audit_json is not None and args.clebsch_audit_json.exists():
        print("[figure] clebsch3d_seed_audit")
        _clebsch_seed_audit_figure(args.clebsch_audit_json.resolve(), outdir)
    print("[write] campaign_summary.json")
    _write_summary(results, outdir)


if __name__ == "__main__":
    main()
