#!/usr/bin/env python3
"""Post-process scalar_mixing production runs with final spectrum fixes."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from functools import lru_cache
import json
import math
from pathlib import Path
import re

from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

import scalar_mixing_stream_run_analysis as base


THETA_LINTHRESH = base.THETA_LINTHRESH
SPECTRUM_NBINS = base.SPECTRUM_NBINS
SCALAR_MAP_LIMIT = base.SCALAR_MAP_LIMIT
SCALAR_MAP_DPI = base.SCALAR_MAP_DPI
SCALAR_MAP_MAX_PIXELS = base.SCALAR_MAP_MAX_PIXELS
SCALAR_CMAP = base.SCALAR_CMAP

SHELL_BINNING = "floor(kmag + 0.5)"
RAW_SPECTRUM_DEFINITION = "discrete integer-shell integrated power"
CONTINUUM_SPECTRUM_DEFINITION = "shell_mean multiplied by continuum shell measure"
REBINNING_DEFINITION = "log-bin mean of continuum spectrum over delta_k"
K_UNITS = "grid_mode"
VELOCITY_DEFINITION = "velocity_power"
SPECTRUM_METHOD_2D = "2d_full_fft"
SPECTRUM_METHOD_SLICE = "slice_averaged_2dfft"
DIRECT_SLICE_IDS = (
    "slice_x_m025",
    "slice_x",
    "slice_x_p025",
    "slice_y_m025",
    "slice_y",
    "slice_y_p025",
    "slice_z_m025",
    "slice_z",
    "slice_z_p025",
)
DIRECT_CENTER_SLICE_IDS = ("slice_x", "slice_y", "slice_z")
DIRECT_SLICE_FILENAME_RE = re.compile(r"^(?P<case>.+)\.(?P<slice_id>slice_[^.]+)\.(?P<dump>\d+)\.bin$")


@dataclass(frozen=True)
class SpectrumGeometry:
    shell_ids: np.ndarray
    shell_count_raw: np.ndarray
    k_shell_raw: np.ndarray
    k_continuum: np.ndarray
    continuum_shell_measure: np.ndarray
    continuum_valid_kmax: int


@dataclass(frozen=True)
class SpectrumResult:
    k_shell_raw: np.ndarray
    shell_sum_raw: np.ndarray
    shell_count_raw: np.ndarray
    shell_mean_raw: np.ndarray
    k_continuum: np.ndarray
    spectrum_continuum: np.ndarray
    continuum_valid_kmax: int


@dataclass(frozen=True)
class RebinSpec:
    starts: np.ndarray
    stops: np.ndarray
    widths: np.ndarray
    centers: np.ndarray


@dataclass(frozen=True)
class ScalarMapPayload2D:
    time: float
    extent: tuple[float, float, float, float]
    centered: np.ndarray


@dataclass(frozen=True)
class MidplanePayload3D:
    time: float
    slices: tuple[tuple[str, np.ndarray, tuple[float, float, float, float], str, str], ...]


@dataclass(frozen=True)
class DirectSliceSet3D:
    dump_token: str
    paths: dict[str, Path]


def _readonly(array: np.ndarray) -> np.ndarray:
    array.setflags(write=False)
    return array


def _discover_scalar_quantities(path: Path) -> tuple[str, ...]:
    raw = base.read_binary_as_athdf(str(path), raw=True)
    scalar_quantities = tuple(sorted(str(name) for name in raw["var_names"] if str(name).startswith("s_")))
    if not scalar_quantities:
        raise ValueError(f"No scalar fields found in {path}")
    return scalar_quantities


def _read_snapshot(path: Path, quantities: list[str] | tuple[str, ...]) -> dict[str, np.ndarray | float]:
    request = list(quantities)
    data = base.read_binary_as_athdf(str(path), quantities=request, dtype=np.float64)
    snapshot: dict[str, np.ndarray | float] = {
        "time": float(data["Time"]),
        "x1v": np.asarray(data["x1v"], dtype=np.float64),
        "x2v": np.asarray(data["x2v"], dtype=np.float64),
    }
    if "x3v" in data:
        snapshot["x3v"] = np.asarray(data["x3v"], dtype=np.float64)
    for quantity in request:
        snapshot[quantity] = np.asarray(data[quantity], dtype=np.float64)
    return snapshot


def _slice_dump_token(path: Path) -> str:
    match = DIRECT_SLICE_FILENAME_RE.match(path.name)
    if match is None:
        raise ValueError(f"Unrecognized slice filename format: {path.name}")
    return str(match.group("dump"))


def _slice_sort_key(path: Path) -> tuple[int, str]:
    dump_token = _slice_dump_token(path)
    if dump_token.isdigit():
        return (int(dump_token), dump_token)
    return (0, dump_token)


def _slice_plane_from_snapshot(
    snapshot: dict[str, np.ndarray | float],
    quantity: str,
    slice_id: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    field = np.asarray(snapshot[quantity], dtype=np.float64)
    x = np.asarray(snapshot["x1v"], dtype=np.float64)
    y = np.asarray(snapshot["x2v"], dtype=np.float64)
    z = np.asarray(snapshot["x3v"], dtype=np.float64)

    if slice_id.startswith("slice_x"):
        return np.asarray(field[:, :, 0], dtype=np.float64), y, z
    if slice_id.startswith("slice_y"):
        return np.asarray(field[:, 0, :], dtype=np.float64), x, z
    if slice_id.startswith("slice_z"):
        return np.asarray(field[0, :, :], dtype=np.float64), x, y
    raise ValueError(f"Unsupported direct slice identifier {slice_id}")


def _collect_direct_slice_sets(input_dir: Path, case_name: str) -> list[DirectSliceSet3D]:
    grouped: dict[str, dict[str, Path]] = {}
    for slice_id in DIRECT_SLICE_IDS:
        for path in sorted(input_dir.glob(f"{case_name}.{slice_id}.*.bin"), key=_slice_sort_key):
            dump_token = _slice_dump_token(path)
            grouped.setdefault(dump_token, {})[slice_id] = path

    direct_sets: list[DirectSliceSet3D] = []
    skipped_tokens: list[str] = []
    for dump_token in sorted(grouped, key=lambda token: (int(token), token) if token.isdigit() else (0, token)):
        paths = grouped[dump_token]
        missing = [slice_id for slice_id in DIRECT_SLICE_IDS if slice_id not in paths]
        if missing:
            skipped_tokens.append(f"{dump_token} missing {','.join(missing)}")
            continue
        direct_sets.append(
            DirectSliceSet3D(
                dump_token=dump_token,
                paths={slice_id: paths[slice_id] for slice_id in DIRECT_SLICE_IDS},
            )
        )

    if skipped_tokens:
        base._warn(
            "Skipped incomplete direct 3D slice sets for "
            f"{case_name}: {'; '.join(skipped_tokens[:6])}"
            + ("; ..." if len(skipped_tokens) > 6 else "")
        )
    return direct_sets


def _read_direct_slice_set(
    slice_set: DirectSliceSet3D,
    quantities: list[str] | tuple[str, ...],
) -> dict[str, dict[str, np.ndarray | float]]:
    return {
        slice_id: _read_snapshot(path, quantities)
        for slice_id, path in slice_set.paths.items()
    }


def _direct_slice_grid_shape(
    slice_snapshots: dict[str, dict[str, np.ndarray | float]]
) -> tuple[int, int, int]:
    x = np.asarray(slice_snapshots["slice_z"]["x1v"], dtype=np.float64)
    y = np.asarray(slice_snapshots["slice_z"]["x2v"], dtype=np.float64)
    z = np.asarray(slice_snapshots["slice_x"]["x3v"], dtype=np.float64)
    return (int(z.size), int(y.size), int(x.size))


def _direct_slice_mean(
    slice_snapshots: dict[str, dict[str, np.ndarray | float]],
    quantity: str,
) -> float:
    plane_means = []
    for slice_id in DIRECT_SLICE_IDS:
        plane, _, _ = _slice_plane_from_snapshot(slice_snapshots[slice_id], quantity, slice_id)
        plane_means.append(float(np.mean(plane)))
    return float(np.mean(plane_means))


def _direct_slice_theta_entry(
    slice_snapshots: dict[str, dict[str, np.ndarray | float]],
    quantity: str,
    time_value: float,
) -> dict[str, np.ndarray | float | int]:
    slice_results = []
    for slice_id in DIRECT_SLICE_IDS:
        plane, _, _ = _slice_plane_from_snapshot(slice_snapshots[slice_id], quantity, slice_id)
        slice_result, _ = _theta_spectrum(plane)
        slice_results.append(slice_result)
    entry = _series_entry(time_value, _average_shell_sums(slice_results))
    entry["spectrum_slice_count"] = len(slice_results)
    return entry


def _direct_slice_velocity_payload(
    slice_snapshots: dict[str, dict[str, np.ndarray | float]]
) -> dict[str, np.ndarray | int | str]:
    slice_results = []
    for slice_id in DIRECT_SLICE_IDS:
        plane_x, _, _ = _slice_plane_from_snapshot(slice_snapshots[slice_id], "velx", slice_id)
        plane_y, _, _ = _slice_plane_from_snapshot(slice_snapshots[slice_id], "vely", slice_id)
        plane_z, _, _ = _slice_plane_from_snapshot(slice_snapshots[slice_id], "velz", slice_id)
        slice_results.append(_velocity_spectrum_2d(plane_x, plane_y, plane_z))
    payload = _velocity_payload(_average_shell_sums(slice_results))
    payload["spectrum_slice_count"] = len(slice_results)
    return payload


def _prepare_direct_midplane_payload(
    slice_snapshots: dict[str, dict[str, np.ndarray | float]],
    quantity: str,
    time_value: float,
) -> MidplanePayload3D:
    mean_value = _direct_slice_mean(slice_snapshots, quantity)
    yz_plane, y, z = _slice_plane_from_snapshot(slice_snapshots["slice_x"], quantity, "slice_x")
    xz_plane, x, z_xz = _slice_plane_from_snapshot(slice_snapshots["slice_y"], quantity, "slice_y")
    xy_plane, x_xy, y_xy = _slice_plane_from_snapshot(slice_snapshots["slice_z"], quantity, "slice_z")

    slices = (
        (
            "x-y midplane",
            np.asarray(base._downsample_plane(xy_plane - mean_value), dtype=np.float32),
            (float(x_xy[0]), float(x_xy[-1]), float(y_xy[0]), float(y_xy[-1])),
            "x",
            "y",
        ),
        (
            "x-z midplane",
            np.asarray(base._downsample_plane(xz_plane - mean_value), dtype=np.float32),
            (float(x[0]), float(x[-1]), float(z_xz[0]), float(z_xz[-1])),
            "x",
            "z",
        ),
        (
            "y-z midplane",
            np.asarray(base._downsample_plane(yz_plane - mean_value), dtype=np.float32),
            (float(y[0]), float(y[-1]), float(z[0]), float(z[-1])),
            "y",
            "z",
        ),
    )
    return MidplanePayload3D(time=float(time_value), slices=slices)


def _scalar_file_tag(case_name: str, scalar_quantity: str, scalar_count: int) -> str:
    if scalar_count <= 1:
        return f"_{case_name}"
    return f"_{case_name}_{scalar_quantity}"


def _scalar_title_label(case_name: str, scalar_quantity: str, scalar_count: int) -> str:
    if scalar_count <= 1:
        return case_name
    return f"{case_name}\n{scalar_quantity}"


def _grid_size_from_field(field: np.ndarray) -> int:
    array = base._canonical_field(field)
    if array.ndim not in (2, 3):
        raise ValueError(f"Unsupported scalar field shape {array.shape}")
    if len(set(int(length) for length in array.shape)) != 1:
        raise ValueError(
            "scalar_mixing_stream_run_analysis_spectrum_final.py assumes square/cubic grids."
        )
    return int(array.shape[0])


def _safe_continuum_kmax(n: int) -> int:
    return max(1, int(math.floor(float(n) / 2.0 - 0.5)))


def _continuum_shell_measure(kvals: np.ndarray, spectrum_dim: int) -> np.ndarray:
    if spectrum_dim == 2:
        return 2.0 * np.pi * np.asarray(kvals, dtype=np.float64)
    if spectrum_dim == 3:
        kvals = np.asarray(kvals, dtype=np.float64)
        return (4.0 * np.pi / 3.0) * (((kvals + 0.5) ** 3) - ((kvals - 0.5) ** 3))
    raise ValueError(f"Unsupported spectrum dimension {spectrum_dim}")


@lru_cache(maxsize=8)
def _spectrum_geometry_2d(n: int) -> SpectrumGeometry:
    k = np.fft.fftfreq(n) * n
    kmag = np.hypot(k[:, None], k[None, :])
    shell_ids = np.floor(kmag + 0.5).astype(np.int32, copy=False).ravel()
    shell_count = np.bincount(shell_ids).astype(np.float64, copy=False)
    shell_count_raw = _readonly(shell_count[1:].copy())
    k_shell_raw = _readonly(np.arange(1, shell_count.size, dtype=np.float64))
    continuum_valid_kmax = min(_safe_continuum_kmax(n), int(k_shell_raw.size))
    k_continuum = _readonly(k_shell_raw[:continuum_valid_kmax].copy())
    continuum_shell_measure = _readonly(_continuum_shell_measure(k_continuum, 2))
    return SpectrumGeometry(
        shell_ids=_readonly(np.ascontiguousarray(shell_ids)),
        shell_count_raw=shell_count_raw,
        k_shell_raw=k_shell_raw,
        k_continuum=k_continuum,
        continuum_shell_measure=continuum_shell_measure,
        continuum_valid_kmax=continuum_valid_kmax,
    )


def _spectrum_geometry_3d(n: int) -> SpectrumGeometry:
    k = np.fft.fftfreq(n) * n
    kmag = np.sqrt(
        k[:, None, None] * k[:, None, None]
        + k[None, :, None] * k[None, :, None]
        + k[None, None, :] * k[None, None, :]
    )
    shell_ids = np.floor(kmag + 0.5).astype(np.int32, copy=False).ravel()
    shell_count = np.bincount(shell_ids).astype(np.float64, copy=False)
    shell_count_raw = shell_count[1:].copy()
    k_shell_raw = np.arange(1, shell_count.size, dtype=np.float64)
    continuum_valid_kmax = min(_safe_continuum_kmax(n), int(k_shell_raw.size))
    k_continuum = k_shell_raw[:continuum_valid_kmax].copy()
    continuum_shell_measure = _continuum_shell_measure(k_continuum, 3)
    return SpectrumGeometry(
        shell_ids=shell_ids,
        shell_count_raw=shell_count_raw,
        k_shell_raw=k_shell_raw,
        k_continuum=k_continuum,
        continuum_shell_measure=continuum_shell_measure,
        continuum_valid_kmax=continuum_valid_kmax,
    )


def _fft_power(transform: np.ndarray, norm_size: int) -> np.ndarray:
    real_part = np.asarray(transform.real, dtype=np.float64)
    imag_part = np.asarray(transform.imag, dtype=np.float64)
    return (real_part * real_part + imag_part * imag_part) / float(norm_size * norm_size)


def _shell_sum_from_power(power: np.ndarray, geometry: SpectrumGeometry) -> np.ndarray:
    shell_sum = np.bincount(
        geometry.shell_ids,
        weights=np.ravel(power),
        minlength=int(geometry.k_shell_raw.size) + 1,
    )
    return np.asarray(shell_sum[1:], dtype=np.float64)


def _make_spectrum_result(shell_sum_raw: np.ndarray, geometry: SpectrumGeometry) -> SpectrumResult:
    shell_mean_raw = np.zeros_like(shell_sum_raw, dtype=np.float64)
    valid = geometry.shell_count_raw > 0.0
    shell_mean_raw[valid] = shell_sum_raw[valid] / geometry.shell_count_raw[valid]
    continuum_size = geometry.k_continuum.size
    spectrum_continuum = (
        shell_mean_raw[:continuum_size] * geometry.continuum_shell_measure
    ).astype(np.float64, copy=False)
    return SpectrumResult(
        k_shell_raw=geometry.k_shell_raw,
        shell_sum_raw=np.asarray(shell_sum_raw, dtype=np.float64),
        shell_count_raw=geometry.shell_count_raw,
        shell_mean_raw=shell_mean_raw,
        k_continuum=geometry.k_continuum,
        spectrum_continuum=spectrum_continuum,
        continuum_valid_kmax=geometry.continuum_valid_kmax,
    )


def _theta_spectrum(field: np.ndarray) -> tuple[SpectrumResult, float]:
    array = base._canonical_field(field)
    n = _grid_size_from_field(array)
    theta_mean = float(np.mean(array))
    fluctuating = array - theta_mean
    shat = np.fft.fftn(fluctuating)
    power = _fft_power(shat, array.size)
    geometry = _spectrum_geometry_2d(n) if array.ndim == 2 else _spectrum_geometry_3d(n)
    return _make_spectrum_result(_shell_sum_from_power(power, geometry), geometry), theta_mean


def _velocity_spectrum_2d(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> SpectrumResult:
    vx2d = base._canonical_field(vx)
    n = _grid_size_from_field(vx2d)
    vy2d = base._canonical_field(vy)
    vz2d = base._canonical_field(vz)
    geometry = _spectrum_geometry_2d(n)
    energy = (
        _fft_power(np.fft.fftn(vx2d), vx2d.size)
        + _fft_power(np.fft.fftn(vy2d), vy2d.size)
        + _fft_power(np.fft.fftn(vz2d), vz2d.size)
    )
    return _make_spectrum_result(_shell_sum_from_power(energy, geometry), geometry)


@lru_cache(maxsize=32)
def _log_rebin_spec(continuum_valid_kmax: int, nbins: int) -> RebinSpec:
    edge_values = np.unique(
        np.rint(np.geomspace(1.0, continuum_valid_kmax + 1.0, nbins + 1)).astype(int)
    )
    if edge_values[0] != 1:
        edge_values = np.insert(edge_values, 0, 1)
    if edge_values[-1] != continuum_valid_kmax + 1:
        edge_values = np.append(edge_values, continuum_valid_kmax + 1)

    starts = []
    stops = []
    widths = []
    centers = []
    for lo, hi in zip(edge_values[:-1], edge_values[1:]):
        if hi <= lo:
            continue
        starts.append(lo - 1)
        stops.append(hi - 1)
        widths.append(hi - lo)
        centers.append(math.sqrt(lo * (hi - 1)))

    return RebinSpec(
        starts=_readonly(np.asarray(starts, dtype=np.int64)),
        stops=_readonly(np.asarray(stops, dtype=np.int64)),
        widths=_readonly(np.asarray(widths, dtype=np.float64)),
        centers=_readonly(np.asarray(centers, dtype=np.float64)),
    )


def _log_rebin_spectrum_density(
    continuum_spectrum: np.ndarray,
    continuum_valid_kmax: int,
    nbins: int = SPECTRUM_NBINS,
) -> tuple[np.ndarray, np.ndarray]:
    if continuum_valid_kmax < 1 or continuum_spectrum.size == 0:
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)

    spec = _log_rebin_spec(continuum_valid_kmax, nbins)
    rebinned = np.empty(spec.centers.size, dtype=np.float64)
    for idx, (start, stop, width) in enumerate(zip(spec.starts, spec.stops, spec.widths, strict=True)):
        rebinned[idx] = np.sum(continuum_spectrum[int(start) : int(stop)], dtype=np.float64) / width
    return spec.centers, rebinned


def _series_entry(time_value: float, result: SpectrumResult) -> dict[str, np.ndarray | float | int]:
    k_rebinned, spec_rebinned = _log_rebin_spectrum_density(
        result.spectrum_continuum,
        result.continuum_valid_kmax,
        nbins=SPECTRUM_NBINS,
    )
    return {
        "time": float(time_value),
        "theta_k_shell_raw": result.k_shell_raw,
        "theta_shell_sum_raw": result.shell_sum_raw,
        "theta_shell_count_raw": result.shell_count_raw,
        "theta_shell_mean_raw": result.shell_mean_raw,
        "theta_k_continuum": result.k_continuum,
        "theta_spectrum_continuum": result.spectrum_continuum,
        "theta_k_rebinned": k_rebinned,
        "theta_spectrum_rebinned_density": spec_rebinned,
        "continuum_valid_kmax": result.continuum_valid_kmax,
    }


def _velocity_payload(result: SpectrumResult) -> dict[str, np.ndarray | int]:
    k_rebinned, spec_rebinned = _log_rebin_spectrum_density(
        result.spectrum_continuum,
        result.continuum_valid_kmax,
        nbins=SPECTRUM_NBINS,
    )
    return {
        "velocity_k_shell_raw": result.k_shell_raw,
        "velocity_shell_sum_raw": result.shell_sum_raw,
        "velocity_shell_count_raw": result.shell_count_raw,
        "velocity_shell_mean_raw": result.shell_mean_raw,
        "velocity_k_continuum": result.k_continuum,
        "velocity_spectrum_continuum": result.spectrum_continuum,
        "velocity_k_rebinned": k_rebinned,
        "velocity_spectrum_rebinned_density": spec_rebinned,
        "continuum_valid_kmax": result.continuum_valid_kmax,
    }


def _average_shell_sums(results: list[SpectrumResult]) -> SpectrumResult:
    if not results:
        raise ValueError("No slice spectra available for averaging.")

    reference = results[0]
    stacked = np.vstack([np.asarray(result.shell_sum_raw, dtype=np.float64) for result in results])
    shell_sum_raw = np.mean(stacked, axis=0)
    geometry = SpectrumGeometry(
        shell_ids=np.asarray([], dtype=np.int32),
        shell_count_raw=reference.shell_count_raw,
        k_shell_raw=reference.k_shell_raw,
        k_continuum=reference.k_continuum,
        continuum_shell_measure=_continuum_shell_measure(reference.k_continuum, 2),
        continuum_valid_kmax=reference.continuum_valid_kmax,
    )
    return _make_spectrum_result(shell_sum_raw, geometry)


def _slice_averaged_theta_entry(
    field: np.ndarray,
    time_value: float,
) -> dict[str, np.ndarray | float | int | str]:
    array = base._canonical_field(field)
    slice_results = []
    for axis in range(array.ndim):
        for index in base._slice_indices(array.shape[axis]):
            plane = np.take(array, indices=index, axis=axis)
            slice_result, _ = _theta_spectrum(plane)
            slice_results.append(slice_result)

    entry = _series_entry(time_value, _average_shell_sums(slice_results))
    entry["spectrum_method"] = SPECTRUM_METHOD_SLICE
    entry["spectrum_slice_count"] = len(slice_results)
    return entry


def _slice_averaged_velocity_payload(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> dict[str, np.ndarray | int | str]:
    vx3d = base._canonical_field(vx)
    vy3d = base._canonical_field(vy)
    vz3d = base._canonical_field(vz)
    slice_results = []
    for axis in range(vx3d.ndim):
        for index in base._slice_indices(vx3d.shape[axis]):
            plane_x = np.take(vx3d, indices=index, axis=axis)
            plane_y = np.take(vy3d, indices=index, axis=axis)
            plane_z = np.take(vz3d, indices=index, axis=axis)
            slice_results.append(_velocity_spectrum_2d(plane_x, plane_y, plane_z))

    payload = _velocity_payload(_average_shell_sums(slice_results))
    payload["spectrum_method"] = SPECTRUM_METHOD_SLICE
    payload["spectrum_slice_count"] = len(slice_results)
    return payload


def _fixed_spectrum_metadata(metadata: dict[str, object], continuum_valid_kmax: int) -> dict[str, object]:
    updated = dict(metadata)
    updated["shell_binning"] = SHELL_BINNING
    updated["raw_spectrum_definition"] = RAW_SPECTRUM_DEFINITION
    updated["continuum_spectrum_definition"] = CONTINUUM_SPECTRUM_DEFINITION
    updated["rebinning_definition"] = REBINNING_DEFINITION
    updated["continuum_valid_kmax"] = continuum_valid_kmax
    updated["k_units"] = K_UNITS
    updated["velocity_definition"] = VELOCITY_DEFINITION
    return updated


def _metadata_npz_fields(
    metadata: dict[str, object],
    grid_shape: tuple[int, ...],
    continuum_valid_kmax: int,
) -> dict[str, np.ndarray]:
    meta = _fixed_spectrum_metadata(metadata, continuum_valid_kmax)
    return {
        "dim": np.asarray(str(meta.get("dim", ""))),
        "grid_dim": np.int64(int(meta.get("grid_dim", 0))),
        "grid_shape": np.asarray(grid_shape, dtype=np.int64),
        "case_method": np.asarray(str(meta.get("method", ""))),
        "turb_expo": np.float64(float(meta.get("turb_expo", 0.0))),
        "scalar_diffusivity": np.float64(float(meta.get("scalar_diffusivity", 0.0))),
        "turb_rseed": np.int64(int(meta.get("turb_rseed", 0))),
        "spectrum_method": np.asarray(str(meta.get("spectrum_method", ""))),
        "spectrum_slice_count": np.int64(int(meta.get("spectrum_slice_count", 0))),
        "case_name": np.asarray(str(meta.get("case_name", ""))),
        "analyzed_dump": np.asarray(str(meta.get("analyzed_dump", ""))),
        "scalar_quantity": np.asarray(str(meta.get("scalar_quantity", ""))),
        "scalar_count": np.int64(int(meta.get("scalar_count", 0))),
        "shell_binning": np.asarray(str(meta["shell_binning"])),
        "raw_spectrum_definition": np.asarray(str(meta["raw_spectrum_definition"])),
        "continuum_spectrum_definition": np.asarray(str(meta["continuum_spectrum_definition"])),
        "rebinning_definition": np.asarray(str(meta["rebinning_definition"])),
        "continuum_valid_kmax": np.int64(int(meta["continuum_valid_kmax"])),
        "k_units": np.asarray(str(meta["k_units"])),
        "velocity_definition": np.asarray(str(meta["velocity_definition"])),
    }


def _save_theta_spectrum_npz(
    output_path: Path,
    metadata: dict[str, object],
    grid_shape: tuple[int, ...],
    entry: dict[str, np.ndarray | float | int],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        time=np.float64(float(entry["time"])),
        theta_k_shell_raw=np.asarray(entry["theta_k_shell_raw"], dtype=np.float64),
        theta_shell_sum_raw=np.asarray(entry["theta_shell_sum_raw"], dtype=np.float64),
        theta_shell_count_raw=np.asarray(entry["theta_shell_count_raw"], dtype=np.float64),
        theta_shell_mean_raw=np.asarray(entry["theta_shell_mean_raw"], dtype=np.float64),
        theta_k_continuum=np.asarray(entry["theta_k_continuum"], dtype=np.float64),
        theta_spectrum_continuum=np.asarray(entry["theta_spectrum_continuum"], dtype=np.float64),
        theta_k_rebinned=np.asarray(entry["theta_k_rebinned"], dtype=np.float64),
        theta_spectrum_rebinned_density=np.asarray(
            entry["theta_spectrum_rebinned_density"],
            dtype=np.float64,
        ),
        **_metadata_npz_fields(metadata, grid_shape, int(entry["continuum_valid_kmax"])),
    )


def _save_theta_series_npz(
    output_path: Path,
    metadata: dict[str, object],
    dump_paths: list[Path],
    grid_shape: tuple[int, ...],
    spectra: list[dict[str, np.ndarray | float | int]],
    velocity_payload: dict[str, np.ndarray | int | str],
) -> None:
    continuum_valid_kmax = int(spectra[0]["continuum_valid_kmax"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        times=np.asarray([float(entry["time"]) for entry in spectra], dtype=np.float64),
        dump_paths=np.asarray([str(path) for path in dump_paths]),
        theta_k_shell_raw=np.asarray(spectra[0]["theta_k_shell_raw"], dtype=np.float64),
        theta_shell_sum_raw=np.vstack(
            [np.asarray(entry["theta_shell_sum_raw"], dtype=np.float64) for entry in spectra]
        ),
        theta_shell_count_raw=np.asarray(spectra[0]["theta_shell_count_raw"], dtype=np.float64),
        theta_shell_mean_raw=np.vstack(
            [np.asarray(entry["theta_shell_mean_raw"], dtype=np.float64) for entry in spectra]
        ),
        theta_k_continuum=np.asarray(spectra[0]["theta_k_continuum"], dtype=np.float64),
        theta_spectrum_continuum=np.vstack(
            [np.asarray(entry["theta_spectrum_continuum"], dtype=np.float64) for entry in spectra]
        ),
        theta_k_rebinned=np.asarray(spectra[0]["theta_k_rebinned"], dtype=np.float64),
        theta_spectrum_rebinned_density=np.vstack(
            [np.asarray(entry["theta_spectrum_rebinned_density"], dtype=np.float64) for entry in spectra]
        ),
        velocity_k_shell_raw=np.asarray(velocity_payload["velocity_k_shell_raw"], dtype=np.float64),
        velocity_shell_sum_raw=np.asarray(velocity_payload["velocity_shell_sum_raw"], dtype=np.float64),
        velocity_shell_count_raw=np.asarray(
            velocity_payload["velocity_shell_count_raw"],
            dtype=np.float64,
        ),
        velocity_shell_mean_raw=np.asarray(velocity_payload["velocity_shell_mean_raw"], dtype=np.float64),
        velocity_k_continuum=np.asarray(velocity_payload["velocity_k_continuum"], dtype=np.float64),
        velocity_spectrum_continuum=np.asarray(
            velocity_payload["velocity_spectrum_continuum"],
            dtype=np.float64,
        ),
        velocity_k_rebinned=np.asarray(velocity_payload["velocity_k_rebinned"], dtype=np.float64),
        velocity_spectrum_rebinned_density=np.asarray(
            velocity_payload["velocity_spectrum_rebinned_density"],
            dtype=np.float64,
        ),
        **_metadata_npz_fields(metadata, grid_shape, continuum_valid_kmax),
    )


def _format_two_line_spectrum_title(title: str) -> str:
    if "\n" in title:
        return title
    for token in (" Slice-Averaged ", " slice-averaged ", " and Frozen ", " at t = "):
        if token in title:
            return title.replace(token, f"\n{token.lstrip()}", 1)
    return title


def _plot_theta_spectrum(
    output: Path,
    title: str,
    spectra: list[dict[str, np.ndarray | float | int]],
    velocity_curve: tuple[np.ndarray, np.ndarray] | None = None,
    velocity_label: str = "velocity spectrum",
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
        kvals = np.asarray(entry["theta_k_rebinned"], dtype=np.float64)
        spec = np.asarray(entry["theta_spectrum_rebinned_density"], dtype=np.float64)
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
                label=velocity_label,
            )

    ax.set_xlabel("k")
    ax.set_ylabel(r"$E_\theta(k)$")
    ax.set_title(_format_two_line_spectrum_title(title))
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


def _prepare_scalar_map_payload(
    field: np.ndarray,
    time_value: float,
    x: np.ndarray,
    y: np.ndarray,
) -> ScalarMapPayload2D:
    scalar = base._canonical_field(np.asarray(field, dtype=np.float64))
    centered = scalar - float(np.mean(scalar))
    stride_y = max(1, math.ceil(centered.shape[0] / SCALAR_MAP_MAX_PIXELS))
    stride_x = max(1, math.ceil(centered.shape[1] / SCALAR_MAP_MAX_PIXELS))
    downsampled = np.asarray(centered[::stride_y, ::stride_x], dtype=np.float32)
    return ScalarMapPayload2D(
        time=float(time_value),
        extent=(float(x[0]), float(x[-1]), float(y[0]), float(y[-1])),
        centered=downsampled,
    )


def _plot_scalar_map_latest(payload: ScalarMapPayload2D, output: Path) -> None:
    norm = colors.SymLogNorm(
        linthresh=THETA_LINTHRESH,
        vmin=-SCALAR_MAP_LIMIT,
        vmax=SCALAR_MAP_LIMIT,
        base=10.0,
    )
    fig, ax = plt.subplots(figsize=(6.8, 6.0))
    image = ax.imshow(
        payload.centered,
        origin="lower",
        cmap=SCALAR_CMAP,
        extent=list(payload.extent),
        norm=norm,
        aspect="equal",
    )
    ax.set_title(f"t = {payload.time:.3f}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label(r"$\theta - \langle \theta \rangle$")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=SCALAR_MAP_DPI, bbox_inches="tight")
    plt.close(fig)


def _plot_scalar_map_series(payloads: list[ScalarMapPayload2D], output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    norm = colors.SymLogNorm(
        linthresh=THETA_LINTHRESH,
        vmin=-SCALAR_MAP_LIMIT,
        vmax=SCALAR_MAP_LIMIT,
        base=10.0,
    )
    for index, payload in enumerate(payloads):
        fig, ax = plt.subplots(figsize=(6.8, 6.0))
        image = ax.imshow(
            payload.centered,
            origin="lower",
            cmap=SCALAR_CMAP,
            extent=list(payload.extent),
            norm=norm,
            aspect="equal",
        )
        ax.set_title(f"t = {payload.time:.3f}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
        cbar.set_label(r"$\theta - \langle \theta \rangle$")
        fig.savefig(output / f"frame_{index:04d}.png", dpi=SCALAR_MAP_DPI, bbox_inches="tight")
        plt.close(fig)


def _prepare_midplane_payload(
    field: np.ndarray,
    time_value: float,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
) -> MidplanePayload3D:
    scalar = base._canonical_field(np.asarray(field, dtype=np.float64))
    centered = scalar - float(np.mean(scalar))
    iz = centered.shape[0] // 2
    iy = centered.shape[1] // 2
    ix = centered.shape[2] // 2

    slices = (
        (
            "x-y midplane",
            np.asarray(base._downsample_plane(centered[iz, :, :]), dtype=np.float32),
            (float(x[0]), float(x[-1]), float(y[0]), float(y[-1])),
            "x",
            "y",
        ),
        (
            "x-z midplane",
            np.asarray(base._downsample_plane(centered[:, iy, :]), dtype=np.float32),
            (float(x[0]), float(x[-1]), float(z[0]), float(z[-1])),
            "x",
            "z",
        ),
        (
            "y-z midplane",
            np.asarray(base._downsample_plane(centered[:, :, ix]), dtype=np.float32),
            (float(y[0]), float(y[-1]), float(z[0]), float(z[-1])),
            "y",
            "z",
        ),
    )
    return MidplanePayload3D(time=float(time_value), slices=slices)


def _plot_scalar_midplanes_latest(payload: MidplanePayload3D, output: Path) -> None:
    theta_absmax = THETA_LINTHRESH
    for _, plane, _, _, _ in payload.slices:
        theta_absmax = max(theta_absmax, float(np.max(np.abs(plane))))

    norm = colors.SymLogNorm(
        linthresh=THETA_LINTHRESH,
        vmin=-theta_absmax,
        vmax=theta_absmax,
        base=10.0,
    )

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), squeeze=False)
    image = None
    for ax, (title, plane, extent, xlabel, ylabel) in zip(axes[0], payload.slices, strict=True):
        image = ax.imshow(
            plane,
            origin="lower",
            cmap=SCALAR_CMAP,
            extent=list(extent),
            norm=norm,
            aspect="equal",
        )
        ax.set_title(f"{title}\nt = {payload.time:.3f}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    if image is not None:
        cbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
        cbar.set_label(r"$\theta - \langle \theta \rangle$")

    fig.suptitle("Scalar slice midplanes", fontsize=15)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_scalar_midplanes_series(payloads: list[MidplanePayload3D], output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    theta_absmax = THETA_LINTHRESH
    for payload in payloads:
        for _, plane, _, _, _ in payload.slices:
            theta_absmax = max(theta_absmax, float(np.max(np.abs(plane))))

    norm = colors.SymLogNorm(
        linthresh=THETA_LINTHRESH,
        vmin=-theta_absmax,
        vmax=theta_absmax,
        base=10.0,
    )

    for index, payload in enumerate(payloads):
        fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), squeeze=False)
        image = None
        for ax, (title, plane, extent, xlabel, ylabel) in zip(axes[0], payload.slices, strict=True):
            image = ax.imshow(
                plane,
                origin="lower",
                cmap=SCALAR_CMAP,
                extent=list(extent),
                norm=norm,
                aspect="equal",
            )
            ax.set_title(f"{title}\nt = {payload.time:.3f}")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

        if image is not None:
            cbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
            cbar.set_label(r"$\theta - \langle \theta \rangle$")

        fig.suptitle("Scalar slice midplanes", fontsize=15)
        fig.savefig(output / f"frame_{index:04d}.png", dpi=220, bbox_inches="tight")
        plt.close(fig)


def _spectrum_filenames(dim: int, latest_only: bool) -> tuple[str, str]:
    if dim == 2:
        return (
            "theta_power_spectrum" if latest_only else "scalar_velocity_spectra",
            "scalar_map" if latest_only else "scalar_maps_redshift",
        )
    return (
        "theta_slice_spectrum" if latest_only else "scalar_velocity_slice_spectra",
        "scalar_slice_midplanes" if latest_only else "scalar_slice_midplanes_series",
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
    input_dir = None
    if dump_paths:
        input_dir = dump_paths[0].parent
    elif "output_dir" in metadata:
        input_dir = Path(str(metadata["output_dir"])) / "bin"
    elif "input_dir" in metadata:
        input_dir = Path(str(metadata["input_dir"]))

    direct_slice_sets = _collect_direct_slice_sets(input_dir, case_name) if input_dir is not None else []
    use_direct_3d_slices = bool(direct_slice_sets)

    if not dump_paths and not use_direct_3d_slices:
        base._warn(f"no hydro_w or direct slice files found for {case_name}")
        return {
            "case_name": case_name,
            "status": "missing_outputs",
            "dump_count": 0,
            "analyzed_dump": None,
            "analysis_dir": str(analysis_dir),
        }

    analysis_dir.mkdir(parents=True, exist_ok=True)

    if use_direct_3d_slices:
        selected_slice_sets = [direct_slice_sets[-1]] if latest_only else direct_slice_sets
        selected_paths = [slice_set.paths["slice_z"] for slice_set in selected_slice_sets]
        scalar_quantities = _discover_scalar_quantities(selected_slice_sets[-1].paths["slice_z"])
        latest_slice_snapshots = _read_direct_slice_set(selected_slice_sets[-1], scalar_quantities)
        dim = 3
        grid_shape = _direct_slice_grid_shape(latest_slice_snapshots)
        scalar_count = len(scalar_quantities)
        latest_time = float(latest_slice_snapshots["slice_z"]["time"])
    else:
        selected_paths = [dump_paths[-1]] if latest_only else dump_paths
        scalar_quantities = _discover_scalar_quantities(selected_paths[-1])
        latest_snapshot = _read_snapshot(selected_paths[-1], scalar_quantities)
        first_scalar = scalar_quantities[0]
        scalar_field = base._canonical_field(np.asarray(latest_snapshot[first_scalar], dtype=np.float64))
        dim = base._grid_dim(scalar_field)
        grid_shape = tuple(int(v) for v in scalar_field.shape)
        scalar_count = len(scalar_quantities)

    meta_base = dict(metadata)
    meta_base["case_name"] = case_name
    meta_base["dim"] = f"{dim}d"
    meta_base["grid_dim"] = dim
    meta_base["analyzed_dump"] = str(selected_paths[-1])
    meta_base["spectrum_method"] = SPECTRUM_METHOD_2D if dim == 2 else SPECTRUM_METHOD_SLICE
    meta_base["spectrum_slice_count"] = 1 if dim == 2 else len(DIRECT_SLICE_IDS)
    meta_base["scalar_count"] = scalar_count

    spectrum_base, map_base = _spectrum_filenames(dim, latest_only)

    if dim == 3:
        if use_direct_3d_slices:
            base._warn(
                "3D outputs use the direct slice files for slice-averaged 2D FFT diagnostics."
            )
        else:
            base._warn(
                "3D outputs are slice-averaged 2D FFT diagnostics, not full 3D isotropic spectra."
            )

    if not latest_only:
        spectra_by_scalar: dict[str, list[dict[str, np.ndarray | float | int]]] = {
            scalar_quantity: [] for scalar_quantity in scalar_quantities
        }
        scalar_map_payloads: dict[str, list[ScalarMapPayload2D]] = {
            scalar_quantity: [] for scalar_quantity in scalar_quantities
        }
        midplane_payloads: dict[str, list[MidplanePayload3D]] = {
            scalar_quantity: [] for scalar_quantity in scalar_quantities
        }

        if use_direct_3d_slices:
            midplane_selection = {
                path
                for path in base._select_series_paths(
                    [slice_set.paths["slice_z"] for slice_set in selected_slice_sets],
                    base.THREED_MAP_SERIES_MAX_SNAPSHOTS,
                )
            }
            for slice_set in selected_slice_sets:
                slice_snapshots = _read_direct_slice_set(slice_set, scalar_quantities)
                time_value = float(slice_snapshots["slice_z"]["time"])
                for scalar_quantity in scalar_quantities:
                    spectra_by_scalar[scalar_quantity].append(
                        _direct_slice_theta_entry(slice_snapshots, scalar_quantity, time_value)
                    )
                    if slice_set.paths["slice_z"] in midplane_selection:
                        midplane_payloads[scalar_quantity].append(
                            _prepare_direct_midplane_payload(
                                slice_snapshots,
                                scalar_quantity,
                                time_value,
                            )
                        )
            velocity_slice_snapshots = _read_direct_slice_set(
                selected_slice_sets[0],
                ("velx", "vely", "velz"),
            )
            velocity_data = _direct_slice_velocity_payload(velocity_slice_snapshots)
            velocity_curve = (
                np.asarray(velocity_data["velocity_k_rebinned"], dtype=np.float64),
                np.asarray(velocity_data["velocity_spectrum_rebinned_density"], dtype=np.float64),
            )
            velocity_label = "slice-averaged 2D velocity spectrum"
        else:
            midplane_selection = {
                path
                for path in base._select_series_paths(
                    selected_paths,
                    base.THREED_MAP_SERIES_MAX_SNAPSHOTS,
                )
            } if dim == 3 else set()
            for path in selected_paths:
                snapshot = latest_snapshot if path == selected_paths[-1] else _read_snapshot(path, scalar_quantities)
                time_value = float(snapshot["time"])
                x_snapshot = np.asarray(snapshot["x1v"], dtype=np.float64)
                y_snapshot = np.asarray(snapshot["x2v"], dtype=np.float64)
                z_snapshot = np.asarray(snapshot["x3v"], dtype=np.float64) if dim == 3 else None
                for scalar_quantity in scalar_quantities:
                    field = np.asarray(snapshot[scalar_quantity], dtype=np.float64)
                    if dim == 2:
                        result, _ = _theta_spectrum(field)
                        spectra_by_scalar[scalar_quantity].append(_series_entry(time_value, result))
                        scalar_map_payloads[scalar_quantity].append(
                            _prepare_scalar_map_payload(field, time_value, x_snapshot, y_snapshot)
                        )
                    else:
                        spectra_by_scalar[scalar_quantity].append(
                            _slice_averaged_theta_entry(field, time_value)
                        )
                        if path in midplane_selection and z_snapshot is not None:
                            midplane_payloads[scalar_quantity].append(
                                _prepare_midplane_payload(
                                    field,
                                    time_value,
                                    x_snapshot,
                                    y_snapshot,
                                    z_snapshot,
                                )
                            )

            velocity_snapshot = _read_snapshot(selected_paths[0], ["velx", "vely", "velz"])
            if dim == 2:
                velocity_data = _velocity_payload(
                    _velocity_spectrum_2d(
                        velocity_snapshot["velx"],
                        velocity_snapshot["vely"],
                        velocity_snapshot["velz"],
                    )
                )
                velocity_curve = (
                    np.asarray(velocity_data["velocity_k_rebinned"], dtype=np.float64),
                    np.asarray(velocity_data["velocity_spectrum_rebinned_density"], dtype=np.float64),
                )
                velocity_label = "velocity spectrum"
            else:
                velocity_data = _slice_averaged_velocity_payload(
                    velocity_snapshot["velx"],
                    velocity_snapshot["vely"],
                    velocity_snapshot["velz"],
                )
                velocity_curve = (
                    np.asarray(velocity_data["velocity_k_rebinned"], dtype=np.float64),
                    np.asarray(velocity_data["velocity_spectrum_rebinned_density"], dtype=np.float64),
                )
                velocity_label = "slice-averaged 2D velocity spectrum"

        for scalar_quantity in scalar_quantities:
            scalar_meta = dict(meta_base)
            scalar_meta["scalar_quantity"] = scalar_quantity
            file_tag = _scalar_file_tag(case_name, scalar_quantity, scalar_count)
            title_label = _scalar_title_label(case_name, scalar_quantity, scalar_count)
            spectra = spectra_by_scalar[scalar_quantity]

            _save_theta_series_npz(
                analysis_dir / f"{spectrum_base}{file_tag}.npz",
                scalar_meta,
                selected_paths,
                grid_shape,
                spectra,
                velocity_data,
            )
            _plot_theta_spectrum(
                analysis_dir / f"{spectrum_base}{file_tag}.png",
                (
                    f"{title_label} Scalar Spectra and Frozen Velocity Spectrum"
                    if dim == 2
                    else (
                        f"{title_label} Slice-Averaged 2D Scalar Spectra and "
                        "Frozen Slice-Averaged 2D Velocity Spectrum"
                    )
                ),
                spectra,
                velocity_curve=velocity_curve,
                velocity_label=velocity_label,
                kmax=int(spectra[0]["continuum_valid_kmax"]),
            )
            if dim == 2:
                _plot_scalar_map_series(
                    scalar_map_payloads[scalar_quantity],
                    analysis_dir / f"{map_base}{file_tag}",
                )
            else:
                latest_midplane_payload = (
                    _prepare_direct_midplane_payload(
                        latest_slice_snapshots,
                        scalar_quantity,
                        latest_time,
                    )
                    if use_direct_3d_slices
                    else _prepare_midplane_payload(
                        latest_snapshot[scalar_quantity],
                        float(latest_snapshot["time"]),
                        np.asarray(latest_snapshot["x1v"], dtype=np.float64),
                        np.asarray(latest_snapshot["x2v"], dtype=np.float64),
                        np.asarray(latest_snapshot["x3v"], dtype=np.float64),
                    )
                )
                _plot_scalar_midplanes_latest(
                    latest_midplane_payload,
                    analysis_dir / f"scalar_slice_midplanes{file_tag}.png",
                )
                _plot_scalar_midplanes_series(
                    midplane_payloads[scalar_quantity],
                    analysis_dir / f"{map_base}{file_tag}",
                )
    else:
        for scalar_quantity in scalar_quantities:
            scalar_meta = dict(meta_base)
            scalar_meta["scalar_quantity"] = scalar_quantity
            file_tag = _scalar_file_tag(case_name, scalar_quantity, scalar_count)
            title_label = _scalar_title_label(case_name, scalar_quantity, scalar_count)
            if use_direct_3d_slices:
                latest_entry = _direct_slice_theta_entry(
                    latest_slice_snapshots,
                    scalar_quantity,
                    latest_time,
                )
            else:
                field = np.asarray(latest_snapshot[scalar_quantity], dtype=np.float64)
                if dim == 2:
                    latest_result, _ = _theta_spectrum(field)
                    latest_entry = _series_entry(float(latest_snapshot["time"]), latest_result)
                else:
                    latest_entry = _slice_averaged_theta_entry(
                        field,
                        float(latest_snapshot["time"]),
                    )

            _save_theta_spectrum_npz(
                analysis_dir / f"{spectrum_base}{file_tag}.npz",
                scalar_meta,
                grid_shape,
                latest_entry,
            )
            _plot_theta_spectrum(
                analysis_dir / f"{spectrum_base}{file_tag}.png",
                (
                    f"{title_label} theta spectrum at t = {float(latest_entry['time']):.3f}"
                    if dim == 2
                    else (
                        f"{title_label} slice-averaged 2D theta spectrum at "
                        f"t = {float(latest_entry['time']):.3f}"
                    )
                ),
                [latest_entry],
                kmax=int(latest_entry["continuum_valid_kmax"]),
            )
            if dim == 2:
                field = np.asarray(latest_snapshot[scalar_quantity], dtype=np.float64)
                _plot_scalar_map_latest(
                    _prepare_scalar_map_payload(
                        field,
                        float(latest_snapshot["time"]),
                        np.asarray(latest_snapshot["x1v"], dtype=np.float64),
                        np.asarray(latest_snapshot["x2v"], dtype=np.float64),
                    ),
                    analysis_dir / f"{map_base}{file_tag}.png",
                )
            else:
                latest_midplane_payload = (
                    _prepare_direct_midplane_payload(
                        latest_slice_snapshots,
                        scalar_quantity,
                        latest_time,
                    )
                    if use_direct_3d_slices
                    else _prepare_midplane_payload(
                        latest_snapshot[scalar_quantity],
                        float(latest_snapshot["time"]),
                        np.asarray(latest_snapshot["x1v"], dtype=np.float64),
                        np.asarray(latest_snapshot["x2v"], dtype=np.float64),
                        np.asarray(latest_snapshot["x3v"], dtype=np.float64),
                    )
                )
                _plot_scalar_midplanes_latest(
                    latest_midplane_payload,
                    analysis_dir / f"{map_base}{file_tag}.png",
                )

    available_count = len(direct_slice_sets) if use_direct_3d_slices else len(dump_paths)
    return {
        "case_name": case_name,
        "status": "ok",
        "dump_count": available_count,
        "analyzed_dump": str(selected_paths[-1]),
        "analysis_dir": str(analysis_dir),
        "dim": dim,
        "scalar_count": scalar_count,
        "scalar_quantities": list(scalar_quantities),
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
        dump_paths = sorted(args.input_dir.glob("*.slice_z.*.bin"), key=_slice_sort_key)
    if not dump_paths:
        raise FileNotFoundError(f"No hydro_w or slice_z files found in {args.input_dir}")

    case_name = base._infer_case_name(args.input_dir)
    metadata = base._case_metadata_from_args(case_name, dump_paths[-1])
    metadata["input_dir"] = str(args.input_dir)
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
