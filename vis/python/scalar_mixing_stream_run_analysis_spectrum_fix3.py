#!/usr/bin/env python3
"""Post-process scalar_mixing production runs with shell-integrated continuum spectra."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

import scalar_mixing_stream_run_analysis as base


SHELL_BINNING = "floor(kmag + 0.5)"
SPECTRUM_NORMALIZATION = "mean_power_times_continuous_shell"
PLOTTED_SPECTRUM = "continuum_surface_measure_times_shell_mean"
REBINNING_KIND = "log_bin_density_divided_by_bin_width"
DIAGNOSTIC_SHELL_SUM = "discrete_shell_sum"


def _spectrum_shells_2d(nx: int, ny: int) -> np.ndarray:
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    ky2, kx2 = np.meshgrid(ky, kx, indexing="ij")
    kmag = np.hypot(kx2, ky2)
    return np.floor(kmag + 0.5).astype(np.int64)


def _spectrum_shells_3d(nx: int, ny: int, nz: int) -> np.ndarray:
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    kz = np.fft.fftfreq(nz) * nz
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    kmag = np.sqrt(kx3 * kx3 + ky3 * ky3 + kz3 * kz3)
    return np.floor(kmag + 0.5).astype(np.int64)


def _shell_spectrum_stats(
    power: np.ndarray,
    shell: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    shell_sum = np.bincount(shell.ravel(), weights=power.ravel())
    shell_count = np.bincount(shell.ravel())
    kvals = np.arange(1, len(shell_sum), dtype=np.float64)
    return kvals, shell_sum[1:], shell_count[1:].astype(np.float64)


def _shell_mean(shell_sum: np.ndarray, shell_count: np.ndarray) -> np.ndarray:
    mean = np.zeros_like(shell_sum, dtype=np.float64)
    valid = shell_count > 0.0
    mean[valid] = shell_sum[valid] / shell_count[valid]
    return mean


def _continuum_shell_spectrum(
    kvals: np.ndarray,
    shell_sum: np.ndarray,
    shell_count: np.ndarray,
    *,
    spectrum_dim: int,
) -> np.ndarray:
    shell_mean = _shell_mean(shell_sum, shell_count)
    if spectrum_dim == 2:
        return 2.0 * np.pi * kvals * shell_mean
    if spectrum_dim == 3:
        return 4.0 * np.pi * kvals * kvals * shell_mean
    raise ValueError(f"Unsupported spectrum dimension {spectrum_dim}")


def _theta_spectrum(
    field: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    array = base._canonical_field(field)
    theta_mean = float(np.mean(array))
    fluctuating = array - theta_mean
    shat = np.fft.fftn(fluctuating)
    power = np.abs(shat) ** 2 / array.size**2

    if array.ndim == 2:
        shell = _spectrum_shells_2d(array.shape[1], array.shape[0])
        spectrum_dim = 2
    else:
        shell = _spectrum_shells_3d(array.shape[2], array.shape[1], array.shape[0])
        spectrum_dim = 3

    kvals, shell_sum, shell_count = _shell_spectrum_stats(power, shell)
    continuum = _continuum_shell_spectrum(
        kvals,
        shell_sum,
        shell_count,
        spectrum_dim=spectrum_dim,
    )
    return kvals, shell_sum, shell_count, continuum, theta_mean


def _velocity_spectrum_2d(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vx2d = base._canonical_field(vx)
    vy2d = base._canonical_field(vy)
    vz2d = base._canonical_field(vz)
    uhat_x = np.fft.fftn(vx2d)
    uhat_y = np.fft.fftn(vy2d)
    uhat_z = np.fft.fftn(vz2d)
    energy = (np.abs(uhat_x) ** 2 + np.abs(uhat_y) ** 2 + np.abs(uhat_z) ** 2) / vx2d.size**2
    shell = _spectrum_shells_2d(vx2d.shape[1], vx2d.shape[0])
    return _shell_spectrum_stats(energy, shell)


def _velocity_spectrum_3d(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vx3d = base._canonical_field(vx)
    vy3d = base._canonical_field(vy)
    vz3d = base._canonical_field(vz)
    uhat_x = np.fft.fftn(vx3d)
    uhat_y = np.fft.fftn(vy3d)
    uhat_z = np.fft.fftn(vz3d)
    energy = (np.abs(uhat_x) ** 2 + np.abs(uhat_y) ** 2 + np.abs(uhat_z) ** 2) / vx3d.size**2
    shell = _spectrum_shells_3d(vx3d.shape[2], vx3d.shape[1], vx3d.shape[0])
    return _shell_spectrum_stats(energy, shell)


def _log_rebin_spectrum_density(
    kvals: np.ndarray,
    spectrum: np.ndarray,
    kmax: int | None = None,
    nbins: int = base.SPECTRUM_NBINS,
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
        rebinned.append(np.sum(spectrum[mask]) / float(hi - lo))
    return np.asarray(centers, dtype=np.float64), np.asarray(rebinned, dtype=np.float64)


def _spectrum_payload(
    kvals: np.ndarray,
    shell_sum: np.ndarray,
    shell_count: np.ndarray,
    *,
    spectrum_dim: int,
    continuum_raw: np.ndarray | None = None,
    kmax: int | None = None,
    nbins: int = base.SPECTRUM_NBINS,
) -> dict[str, np.ndarray]:
    shell_mean = _shell_mean(shell_sum, shell_count)
    if continuum_raw is None:
        continuum = _continuum_shell_spectrum(
            kvals,
            shell_sum,
            shell_count,
            spectrum_dim=spectrum_dim,
        )
    else:
        continuum = np.asarray(continuum_raw, dtype=np.float64)
    k_rebinned, spec_rebinned = _log_rebin_spectrum_density(
        kvals,
        continuum,
        kmax=kmax,
        nbins=nbins,
    )
    return {
        "k_raw": np.asarray(kvals, dtype=np.float64),
        "spec_raw": continuum,
        "shell_sum_raw": np.asarray(shell_sum, dtype=np.float64),
        "shell_count_raw": np.asarray(shell_count, dtype=np.float64),
        "shell_mean_raw": shell_mean,
        "spec_continuum_raw": continuum,
        "k_rebinned": k_rebinned,
        "spec_rebinned": spec_rebinned,
    }


def _build_rebinned_entry(
    time_value: float,
    kvals: np.ndarray,
    shell_sum: np.ndarray,
    shell_count: np.ndarray,
    *,
    spectrum_dim: int,
    continuum_raw: np.ndarray | None = None,
    kmax: int | None = None,
    nbins: int = base.SPECTRUM_NBINS,
) -> dict[str, np.ndarray | float]:
    payload: dict[str, np.ndarray | float] = _spectrum_payload(
        kvals,
        shell_sum,
        shell_count,
        spectrum_dim=spectrum_dim,
        continuum_raw=continuum_raw,
        kmax=kmax,
        nbins=nbins,
    )
    payload["time"] = float(time_value)
    return payload


def _average_raw_shell_stats(
    spectra: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not spectra:
        raise ValueError("No slice spectra available for averaging.")

    common_length = min(len(shell_sum) for _, shell_sum, _ in spectra)
    raw_k = np.arange(1, common_length + 1, dtype=np.float64)
    raw_sum = np.mean(
        np.vstack([np.asarray(shell_sum[:common_length], dtype=np.float64) for _, shell_sum, _ in spectra]),
        axis=0,
    )
    raw_count = np.mean(
        np.vstack([np.asarray(shell_count[:common_length], dtype=np.float64) for _, _, shell_count in spectra]),
        axis=0,
    )
    return raw_k, raw_sum, raw_count


def _slice_averaged_theta_entry(
    field: np.ndarray,
    time_value: float,
    *,
    kmax: int,
    nbins: int,
) -> dict[str, np.ndarray | float | int | str]:
    array = base._canonical_field(field)
    slice_spectra: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for axis in range(array.ndim):
        for index in base._slice_indices(array.shape[axis]):
            plane = np.take(array, indices=index, axis=axis)
            kvals, shell_sum, shell_count, _, _ = _theta_spectrum(plane)
            slice_spectra.append((kvals, shell_sum, shell_count))

    raw_k, raw_sum, raw_count = _average_raw_shell_stats(slice_spectra)
    entry = _build_rebinned_entry(
        time_value,
        raw_k,
        raw_sum,
        raw_count,
        spectrum_dim=2,
        kmax=kmax,
        nbins=nbins,
    )
    entry["method"] = "3d_slice_average_2dfft"
    entry["slice_count"] = len(slice_spectra)
    return entry


def _slice_averaged_velocity_entry(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    *,
    kmax: int,
    nbins: int,
) -> tuple[dict[str, np.ndarray], int, str]:
    vx3d = base._canonical_field(vx)
    vy3d = base._canonical_field(vy)
    vz3d = base._canonical_field(vz)

    slice_spectra: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for axis in range(vx3d.ndim):
        for index in base._slice_indices(vx3d.shape[axis]):
            plane_x = np.take(vx3d, indices=index, axis=axis)
            plane_y = np.take(vy3d, indices=index, axis=axis)
            plane_z = np.take(vz3d, indices=index, axis=axis)
            kvals, shell_sum, shell_count = _velocity_spectrum_2d(plane_x, plane_y, plane_z)
            slice_spectra.append((kvals, shell_sum, shell_count))

    raw_k, raw_sum, raw_count = _average_raw_shell_stats(slice_spectra)
    payload = _spectrum_payload(
        raw_k,
        raw_sum,
        raw_count,
        spectrum_dim=2,
        kmax=kmax,
        nbins=nbins,
    )
    return payload, len(slice_spectra), "3d_slice_average_2dfft"


def _fixed_spectrum_metadata(metadata: dict[str, object]) -> dict[str, object]:
    updated = dict(metadata)
    updated["shell_binning"] = SHELL_BINNING
    updated["spectrum_normalization"] = SPECTRUM_NORMALIZATION
    updated["plotted_spectrum"] = PLOTTED_SPECTRUM
    updated["rebinning"] = REBINNING_KIND
    updated["diagnostic_shell_sum"] = DIAGNOSTIC_SHELL_SUM
    case_name = str(updated.get("case_name", ""))
    updated.setdefault(
        "problem_generator",
        "scalar_mixing_perfect_powerlaw" if "perfect_powerlaw" in case_name else "",
    )
    updated.setdefault(
        "turb_spectrum_contract",
        "exact_shell" if "perfect_powerlaw" in case_name else "",
    )
    return updated


def _save_theta_spectrum_npz(
    output_path: Path,
    metadata: dict[str, object],
    grid_shape: tuple[int, ...],
    entry: dict[str, np.ndarray | float],
) -> None:
    meta = _fixed_spectrum_metadata(metadata)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        time=np.float64(float(entry["time"])),
        k_raw=np.asarray(entry["k_raw"], dtype=np.float64),
        theta_spectrum_raw=np.asarray(entry["spec_raw"], dtype=np.float64),
        theta_shell_sum_raw=np.asarray(entry["shell_sum_raw"], dtype=np.float64),
        theta_shell_count_raw=np.asarray(entry["shell_count_raw"], dtype=np.float64),
        theta_shell_mean_raw=np.asarray(entry["shell_mean_raw"], dtype=np.float64),
        theta_spectrum_continuum_raw=np.asarray(entry["spec_continuum_raw"], dtype=np.float64),
        k_rebinned=np.asarray(entry["k_rebinned"], dtype=np.float64),
        theta_spectrum_rebinned=np.asarray(entry["spec_rebinned"], dtype=np.float64),
        dim=np.asarray(str(meta.get("dim", ""))),
        grid_dim=np.int64(int(meta.get("grid_dim", 0))),
        grid_shape=np.asarray(grid_shape, dtype=np.int64),
        method=np.asarray(str(meta.get("method", ""))),
        turb_expo=np.float64(float(meta.get("turb_expo", 0.0))),
        scalar_diffusivity=np.float64(float(meta.get("scalar_diffusivity", 0.0))),
        turb_rseed=np.int64(int(meta.get("turb_rseed", 0))),
        spectrum_method=np.asarray(str(meta.get("spectrum_method", ""))),
        spectrum_slice_count=np.int64(int(meta.get("spectrum_slice_count", 0))),
        theta_spectrum_method=np.asarray(str(meta.get("theta_spectrum_method", ""))),
        theta_spectrum_slice_count=np.int64(int(meta.get("theta_spectrum_slice_count", 0))),
        velocity_spectrum_method=np.asarray(str(meta.get("velocity_spectrum_method", ""))),
        velocity_spectrum_slice_count=np.int64(int(meta.get("velocity_spectrum_slice_count", 0))),
        case_name=np.asarray(str(meta.get("case_name", ""))),
        analyzed_dump=np.asarray(str(meta.get("analyzed_dump", ""))),
        problem_generator=np.asarray(str(meta.get("problem_generator", ""))),
        turb_spectrum_contract=np.asarray(str(meta.get("turb_spectrum_contract", ""))),
        shell_binning=np.asarray(str(meta["shell_binning"])),
        spectrum_normalization=np.asarray(str(meta["spectrum_normalization"])),
        plotted_spectrum=np.asarray(str(meta["plotted_spectrum"])),
        rebinning=np.asarray(str(meta["rebinning"])),
        diagnostic_shell_sum=np.asarray(str(meta["diagnostic_shell_sum"])),
    )


def _save_theta_series_npz(
    output_path: Path,
    metadata: dict[str, object],
    dump_paths: list[Path],
    grid_shape: tuple[int, ...],
    spectra: list[dict[str, np.ndarray | float]],
    velocity_entry: dict[str, np.ndarray],
    *,
    kmax: int,
    nbins: int,
) -> None:
    meta = _fixed_spectrum_metadata(metadata)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        times=np.asarray([float(entry["time"]) for entry in spectra], dtype=np.float64),
        dump_paths=np.asarray([str(path) for path in dump_paths]),
        theta_k_raw=np.asarray(spectra[0]["k_raw"], dtype=np.float64),
        theta_spectrum_raw=np.vstack([np.asarray(entry["spec_raw"], dtype=np.float64) for entry in spectra]),
        theta_shell_sum_raw=np.vstack(
            [np.asarray(entry["shell_sum_raw"], dtype=np.float64) for entry in spectra]
        ),
        theta_shell_count_raw=np.asarray(spectra[0]["shell_count_raw"], dtype=np.float64),
        theta_shell_mean_raw=np.vstack(
            [np.asarray(entry["shell_mean_raw"], dtype=np.float64) for entry in spectra]
        ),
        theta_spectrum_continuum_raw=np.vstack(
            [np.asarray(entry["spec_continuum_raw"], dtype=np.float64) for entry in spectra]
        ),
        theta_k_rebinned=np.asarray(spectra[0]["k_rebinned"], dtype=np.float64),
        theta_spectrum_rebinned=np.vstack(
            [np.asarray(entry["spec_rebinned"], dtype=np.float64) for entry in spectra]
        ),
        velocity_k_raw=np.asarray(velocity_entry["k_raw"], dtype=np.float64),
        velocity_spectrum_raw=np.asarray(velocity_entry["spec_raw"], dtype=np.float64),
        velocity_shell_sum_raw=np.asarray(velocity_entry["shell_sum_raw"], dtype=np.float64),
        velocity_shell_count_raw=np.asarray(velocity_entry["shell_count_raw"], dtype=np.float64),
        velocity_shell_mean_raw=np.asarray(velocity_entry["shell_mean_raw"], dtype=np.float64),
        velocity_spectrum_continuum_raw=np.asarray(velocity_entry["spec_continuum_raw"], dtype=np.float64),
        velocity_k_rebinned=np.asarray(velocity_entry["k_rebinned"], dtype=np.float64),
        velocity_spectrum_rebinned=np.asarray(velocity_entry["spec_rebinned"], dtype=np.float64),
        kmax=np.int64(kmax),
        nbins=np.int64(nbins),
        dim=np.asarray(str(meta.get("dim", ""))),
        grid_dim=np.int64(int(meta.get("grid_dim", 0))),
        grid_shape=np.asarray(grid_shape, dtype=np.int64),
        method=np.asarray(str(meta.get("method", ""))),
        turb_expo=np.float64(float(meta.get("turb_expo", 0.0))),
        scalar_diffusivity=np.float64(float(meta.get("scalar_diffusivity", 0.0))),
        turb_rseed=np.int64(int(meta.get("turb_rseed", 0))),
        spectrum_method=np.asarray(str(meta.get("spectrum_method", ""))),
        spectrum_slice_count=np.int64(int(meta.get("spectrum_slice_count", 0))),
        theta_spectrum_method=np.asarray(str(meta.get("theta_spectrum_method", ""))),
        theta_spectrum_slice_count=np.int64(int(meta.get("theta_spectrum_slice_count", 0))),
        velocity_spectrum_method=np.asarray(str(meta.get("velocity_spectrum_method", ""))),
        velocity_spectrum_slice_count=np.int64(int(meta.get("velocity_spectrum_slice_count", 0))),
        case_name=np.asarray(str(meta.get("case_name", ""))),
        analyzed_dump=np.asarray(str(meta.get("analyzed_dump", ""))),
        problem_generator=np.asarray(str(meta.get("problem_generator", ""))),
        turb_spectrum_contract=np.asarray(str(meta.get("turb_spectrum_contract", ""))),
        shell_binning=np.asarray(str(meta["shell_binning"])),
        spectrum_normalization=np.asarray(str(meta["spectrum_normalization"])),
        plotted_spectrum=np.asarray(str(meta["plotted_spectrum"])),
        rebinning=np.asarray(str(meta["rebinning"])),
        diagnostic_shell_sum=np.asarray(str(meta["diagnostic_shell_sum"])),
    )


def _analyze_case(
    case_name: str,
    dump_paths: list[Path],
    analysis_dir: Path,
    metadata: dict[str, object],
    latest_only: bool,
    save_spectrum_npz: bool,
) -> dict[str, object]:
    del save_spectrum_npz
    if not dump_paths:
        base._warn(f"no hydro_w files found for {case_name}")
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

    latest_snapshot = base._read_snapshot(selected_paths[-1], ["s_00"])
    scalar_field = base._canonical_field(np.asarray(latest_snapshot["scalar"]))
    dim = base._grid_dim(scalar_field)
    spectrum_kmax = base._spectrum_kmax_from_field(scalar_field)

    meta = dict(metadata)
    meta["case_name"] = case_name
    meta["dim"] = f"{dim}d"
    meta["grid_dim"] = dim
    meta["analyzed_dump"] = str(selected_paths[-1])
    meta["theta_spectrum_method"] = "2d_full_fft" if dim == 2 else "3d_slice_average_2dfft"
    meta["theta_spectrum_slice_count"] = (
        1 if dim == 2 else len(base.THREED_SPECTRUM_SLICE_FRACTIONS) * dim
    )
    meta["velocity_spectrum_method"] = "2d_full_fft" if dim == 2 else "3d_full_fft"
    meta["velocity_spectrum_slice_count"] = 1
    meta["spectrum_method"] = meta["theta_spectrum_method"]
    meta["spectrum_slice_count"] = meta["theta_spectrum_slice_count"]

    if dim == 2:
        latest_k, latest_sum, latest_count, latest_continuum, _ = _theta_spectrum(scalar_field)
        latest_entry = _build_rebinned_entry(
            float(latest_snapshot["time"]),
            latest_k,
            latest_sum,
            latest_count,
            spectrum_dim=2,
            continuum_raw=latest_continuum,
            kmax=spectrum_kmax,
            nbins=base.SPECTRUM_NBINS,
        )
    else:
        base._warn(
            "3D spectra use averaged 2D FFTs from multiple orthogonal slices "
            "instead of a full 3D FFT."
        )
        latest_entry = _slice_averaged_theta_entry(
            scalar_field,
            float(latest_snapshot["time"]),
            kmax=spectrum_kmax,
            nbins=base.SPECTRUM_NBINS,
        )

    if not latest_only:
        spectra: list[dict[str, np.ndarray | float]] = []
        for path in selected_paths:
            snapshot = base._read_snapshot(path, ["s_00"])
            if dim == 2:
                kvals, shell_sum, shell_count, continuum, _ = _theta_spectrum(snapshot["scalar"])
                spectra.append(
                    _build_rebinned_entry(
                        float(snapshot["time"]),
                        kvals,
                        shell_sum,
                        shell_count,
                        spectrum_dim=2,
                        continuum_raw=continuum,
                        kmax=spectrum_kmax,
                        nbins=base.SPECTRUM_NBINS,
                    )
                )
            else:
                spectra.append(
                    _slice_averaged_theta_entry(
                        snapshot["scalar"],
                        float(snapshot["time"]),
                        kmax=spectrum_kmax,
                        nbins=base.SPECTRUM_NBINS,
                    )
                )

        velocity_snapshot = base._read_snapshot(selected_paths[0], ["velx", "vely", "velz"])
        if dim == 2:
            kvals_v, shell_sum_v, shell_count_v = _velocity_spectrum_2d(
                velocity_snapshot["velx"],
                velocity_snapshot["vely"],
                velocity_snapshot["velz"],
            )
            velocity_entry = _spectrum_payload(
                kvals_v,
                shell_sum_v,
                shell_count_v,
                spectrum_dim=2,
                kmax=spectrum_kmax,
                nbins=base.SPECTRUM_NBINS,
            )
        else:
            kvals_v, shell_sum_v, shell_count_v = _velocity_spectrum_3d(
                velocity_snapshot["velx"],
                velocity_snapshot["vely"],
                velocity_snapshot["velz"],
            )
            velocity_entry = _spectrum_payload(
                kvals_v,
                shell_sum_v,
                shell_count_v,
                spectrum_dim=3,
                kmax=spectrum_kmax,
                nbins=base.SPECTRUM_NBINS,
            )
            meta["spectrum_method"] = meta["velocity_spectrum_method"]
            meta["spectrum_slice_count"] = meta["velocity_spectrum_slice_count"]

        _save_theta_series_npz(
            analysis_dir / f"scalar_velocity_spectra{file_tag}.npz",
            meta,
            selected_paths,
            tuple(int(v) for v in scalar_field.shape),
            spectra,
            velocity_entry,
            kmax=spectrum_kmax,
            nbins=base.SPECTRUM_NBINS,
        )
        base._plot_theta_spectrum(
            analysis_dir / f"scalar_velocity_spectra{file_tag}.png",
            (
                "Scalar Power Spectra and Frozen Velocity Spectrum"
                if dim == 2
                else "Slice-Averaged Scalar Power Spectra and 3D FFT Velocity Spectrum"
            ),
            spectra,
            velocity_curve=(
                np.asarray(velocity_entry["k_rebinned"], dtype=np.float64),
                np.asarray(velocity_entry["spec_rebinned"], dtype=np.float64),
            ),
            kmax=spectrum_kmax,
        )
        if dim == 2:
            base._plot_scalar_map_series(
                selected_paths,
                analysis_dir / f"scalar_maps_redshift{file_tag}.png",
            )
        else:
            base._plot_scalar_midplanes(
                latest_snapshot,
                analysis_dir / f"scalar_midplanes{file_tag}.png",
            )
            base._plot_scalar_midplanes_series(
                selected_paths,
                analysis_dir / f"scalar_midplanes_series{file_tag}.png",
            )
    else:
        _save_theta_spectrum_npz(
            analysis_dir / f"theta_power_spectrum{file_tag}.npz",
            meta,
            tuple(int(v) for v in scalar_field.shape),
            latest_entry,
        )
        title = (
            f"{case_name} theta power spectrum at t = {float(latest_snapshot['time']):.3f}"
            if dim == 2
            else (
                f"{case_name} slice-averaged theta power spectrum at "
                f"t = {float(latest_snapshot['time']):.3f}"
            )
        )
        base._plot_theta_spectrum(
            analysis_dir / f"theta_power_spectrum{file_tag}.png",
            title,
            [latest_entry],
            kmax=spectrum_kmax,
        )
        if dim == 2:
            base._plot_scalar_map_latest(
                latest_snapshot,
                analysis_dir / f"scalar_map{file_tag}.png",
            )
        else:
            base._plot_scalar_midplanes(
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


def _analyze_campaign(
    manifest_path: Path,
    analysis_root: Path,
    latest_only: bool,
    save_spectrum_npz: bool,
) -> list[dict[str, object]]:
    results = []
    with manifest_path.open(newline="", encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))

    for row in rows:
        case_name = row["case_name"]
        case_analysis_dir = analysis_root / case_name
        dump_paths = sorted((Path(row["output_dir"]) / "bin").glob("*.hydro_w.*.bin"))
        try:
            result = _analyze_case(
                case_name=case_name,
                dump_paths=dump_paths,
                analysis_dir=case_analysis_dir,
                metadata=row,
                latest_only=latest_only,
                save_spectrum_npz=save_spectrum_npz,
            )
        except Exception as exc:
            base._warn(f"analysis failed for {case_name}: {exc}")
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

    case_name = base._infer_case_name(args.input_dir)
    metadata = base._case_metadata_from_args(case_name, dump_paths[-1])
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
