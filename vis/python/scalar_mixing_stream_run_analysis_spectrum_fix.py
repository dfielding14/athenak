#!/usr/bin/env python3
"""Post-process scalar_mixing production runs with corrected radial spectra."""

from __future__ import annotations

from pathlib import Path

import numpy as np

import scalar_mixing_stream_run_analysis as base


SHELL_BINNING = "round"
SPECTRUM_NORMALIZATION = "mean_power_times_continuous_shell"


def _spectrum_shells_2d(nx: int, ny: int) -> np.ndarray:
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    ky2, kx2 = np.meshgrid(ky, kx, indexing="ij")
    return np.round(np.sqrt(kx2 * kx2 + ky2 * ky2)).astype(np.int64)


def _spectrum_shells_3d(nx: int, ny: int, nz: int) -> np.ndarray:
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    kz = np.fft.fftfreq(nz) * nz
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    return np.round(np.sqrt(kx3 * kx3 + ky3 * ky3 + kz3 * kz3)).astype(np.int64)


def _radial_shell_spectrum(
    power: np.ndarray,
    shell: np.ndarray,
    *,
    shell_measure_factor: float,
    shell_measure_exponent: float,
) -> tuple[np.ndarray, np.ndarray]:
    spectrum_sum = np.bincount(shell.ravel(), weights=power.ravel())
    mode_counts = np.bincount(shell.ravel())
    mean_power = np.zeros_like(spectrum_sum, dtype=np.float64)
    valid = mode_counts > 0
    mean_power[valid] = spectrum_sum[valid] / mode_counts[valid]

    kvals = np.arange(1, len(spectrum_sum), dtype=np.float64)
    spectrum = mean_power[1:] * shell_measure_factor * (kvals ** shell_measure_exponent)
    return kvals, spectrum


def _theta_spectrum(field: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    array = base._canonical_field(field)
    theta_mean = float(np.mean(array))
    fluctuating = array - theta_mean
    shat = np.fft.fftn(fluctuating)
    power = np.abs(shat) ** 2 / array.size**2

    if array.ndim == 2:
        shell = _spectrum_shells_2d(array.shape[1], array.shape[0])
        kvals, spectrum = _radial_shell_spectrum(
            power,
            shell,
            shell_measure_factor=2.0 * np.pi,
            shell_measure_exponent=1.0,
        )
    else:
        shell = _spectrum_shells_3d(array.shape[2], array.shape[1], array.shape[0])
        kvals, spectrum = _radial_shell_spectrum(
            power,
            shell,
            shell_measure_factor=4.0 * np.pi,
            shell_measure_exponent=2.0,
        )

    return kvals, spectrum, theta_mean


def _velocity_spectrum_2d(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    vx2d = base._canonical_field(vx)
    vy2d = base._canonical_field(vy)
    vz2d = base._canonical_field(vz)
    uhat_x = np.fft.fftn(vx2d)
    uhat_y = np.fft.fftn(vy2d)
    uhat_z = np.fft.fftn(vz2d)
    energy = (np.abs(uhat_x) ** 2 + np.abs(uhat_y) ** 2 + np.abs(uhat_z) ** 2) / vx2d.size**2
    shell = _spectrum_shells_2d(vx2d.shape[1], vx2d.shape[0])
    return _radial_shell_spectrum(
        energy,
        shell,
        shell_measure_factor=2.0 * np.pi,
        shell_measure_exponent=1.0,
    )


def _fixed_spectrum_metadata(metadata: dict[str, object]) -> dict[str, object]:
    updated = dict(metadata)
    updated["shell_binning"] = SHELL_BINNING
    updated["spectrum_normalization"] = SPECTRUM_NORMALIZATION
    return updated


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
    meta = _fixed_spectrum_metadata(metadata)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        time=np.float64(time_value),
        k_raw=np.asarray(kvals_raw, dtype=np.float64),
        theta_spectrum_raw=np.asarray(spec_raw, dtype=np.float64),
        k_rebinned=np.asarray(kvals_rebinned, dtype=np.float64),
        theta_spectrum_rebinned=np.asarray(spec_rebinned, dtype=np.float64),
        dim=np.asarray(str(meta.get("dim", ""))),
        grid_dim=np.int64(int(meta.get("grid_dim", 0))),
        grid_shape=np.asarray(grid_shape, dtype=np.int64),
        method=np.asarray(str(meta.get("method", ""))),
        turb_expo=np.float64(float(meta.get("turb_expo", 0.0))),
        scalar_diffusivity=np.float64(float(meta.get("scalar_diffusivity", 0.0))),
        turb_rseed=np.int64(int(meta.get("turb_rseed", 0))),
        spectrum_method=np.asarray(str(meta.get("spectrum_method", ""))),
        spectrum_slice_count=np.int64(int(meta.get("spectrum_slice_count", 0))),
        case_name=np.asarray(str(meta.get("case_name", ""))),
        analyzed_dump=np.asarray(str(meta.get("analyzed_dump", ""))),
        shell_binning=np.asarray(str(meta["shell_binning"])),
        spectrum_normalization=np.asarray(str(meta["spectrum_normalization"])),
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
    meta = _fixed_spectrum_metadata(metadata)
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
        dim=np.asarray(str(meta.get("dim", ""))),
        grid_dim=np.int64(int(meta.get("grid_dim", 0))),
        grid_shape=np.asarray(grid_shape, dtype=np.int64),
        method=np.asarray(str(meta.get("method", ""))),
        turb_expo=np.float64(float(meta.get("turb_expo", 0.0))),
        scalar_diffusivity=np.float64(float(meta.get("scalar_diffusivity", 0.0))),
        turb_rseed=np.int64(int(meta.get("turb_rseed", 0))),
        spectrum_method=np.asarray(str(meta.get("spectrum_method", ""))),
        spectrum_slice_count=np.int64(int(meta.get("spectrum_slice_count", 0))),
        case_name=np.asarray(str(meta.get("case_name", ""))),
        analyzed_dump=np.asarray(str(meta.get("analyzed_dump", ""))),
        shell_binning=np.asarray(str(meta["shell_binning"])),
        spectrum_normalization=np.asarray(str(meta["spectrum_normalization"])),
    )


base._spectrum_shells_2d = _spectrum_shells_2d
base._spectrum_shells_3d = _spectrum_shells_3d
base._theta_spectrum = _theta_spectrum
base._velocity_spectrum_2d = _velocity_spectrum_2d
base._save_theta_spectrum_npz = _save_theta_spectrum_npz
base._save_theta_series_npz = _save_theta_series_npz


def main() -> None:
    base.main()


if __name__ == "__main__":
    main()
