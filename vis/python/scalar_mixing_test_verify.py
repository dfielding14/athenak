"""Post-processing verifier for the scalar_mixing_test blob advection-diffusion problem.

The analytic solution is evaluated by taking the exact AthenaK cell-centered initial
condition, transforming it into periodic Fourier modes, and evolving each mode under
uniform advection plus diffusion. This keeps the reference solution aligned with the
discrete initialization used by the problem generator.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np

from bin_convert_new import athinput, read_binary_as_athdf


# Calibrated final-time thresholds for the shipped regression resolutions.
REFERENCE_THRESHOLDS: dict[tuple[str, tuple[int, int, int], float], dict[str, float]] = {
    ("scalar_mixing_blob_diffusion_2d.athinput", (1, 128, 128), 1.0): {
        "l1_max": 8.0e-4,
        "l2_max": 1.1e-3,
        "linf_max": 5.0e-3,
        "mass_rel_max": 1.0e-8,
        "centroid_max": 5.0e-4,
        "width_rel_max": 2.0e-2,
        "vel_uniform_max": 1.0e-12,
    },
    ("scalar_mixing_blob_diffusion_3d.athinput", (48, 48, 48), 1.0): {
        "l1_max": 9.0e-4,
        "l2_max": 2.5e-3,
        "linf_max": 3.5e-2,
        "mass_rel_max": 1.0e-8,
        "centroid_max": 1.0e-2,
        "width_rel_max": 4.0e-2,
        "vel_uniform_max": 1.0e-12,
    },
}
def _active_dimension(shape: tuple[int, int, int]) -> int:
    nz, ny, _ = shape
    return 3 if nz > 1 else (2 if ny > 1 else 1)


def _problem_defaults(config: dict[str, dict[str, Any]],
                      geometry: dict[str, float]) -> dict[str, float]:
    problem = config.get("problem", {})
    dim = geometry["dim"]
    active_min = min(geometry["lx"], geometry["ly"])
    if dim == 3:
        active_min = min(active_min, geometry["lz"])

    return {
        "scalar_inside": float(problem.get("scalar_inside", 1.0)),
        "scalar_outside": float(problem.get("scalar_outside", 0.0)),
        "blob_radius": float(problem.get("blob_radius", 0.125 * active_min)),
        "blob_center_x1": float(problem.get("blob_center_x1",
                                            geometry["x1min"] + 0.25*geometry["lx"])),
        "blob_center_x2": float(problem.get("blob_center_x2",
                                            geometry["x2min"] + 0.25*geometry["ly"])),
        "blob_center_x3": float(problem.get("blob_center_x3",
                                            geometry["x3min"] + 0.25*geometry["lz"])),
        "scalar_diffusivity": float(config.get("hydro", {}).get("scalar_diffusivity", 0.0)),
    }
def initial_scalar_field(config: dict[str, dict[str, Any]],
                         data: dict[str, Any]) -> np.ndarray:
    geometry = _geometry_from_dump(data)
    params = _problem_defaults(config, geometry)
    x1v = np.asarray(data["x1v"], dtype=np.float64)
    x2v = np.asarray(data["x2v"], dtype=np.float64)
    x3v = np.asarray(data["x3v"], dtype=np.float64)
    z3, y3, x3 = np.meshgrid(x3v, x2v, x1v, indexing="ij")

    dx = _periodic_delta(x3, params["blob_center_x1"], geometry["lx"])
    dy = _periodic_delta(y3, params["blob_center_x2"], geometry["ly"])
    if geometry["dim"] == 3:
        dz = _periodic_delta(z3, params["blob_center_x3"], geometry["lz"])
    else:
        dz = 0.0
    inside = dx*dx + dy*dy + dz*dz <= params["blob_radius"]*params["blob_radius"]
    return np.where(inside, params["scalar_inside"], params["scalar_outside"])


def _geometry_from_dump(data: dict[str, Any]) -> dict[str, float]:
    x1f = np.asarray(data["x1f"], dtype=np.float64)
    x2f = np.asarray(data["x2f"], dtype=np.float64)
    x3f = np.asarray(data["x3f"], dtype=np.float64)
    shape = tuple(int(n) for n in np.asarray(data["s_00"]).shape)
    nz, ny, nx = shape
    return {
        "shape": shape,
        "dim": _active_dimension(shape),
        "x1min": float(x1f[0]),
        "x2min": float(x2f[0]),
        "x3min": float(x3f[0]),
        "lx": float(x1f[-1] - x1f[0]),
        "ly": float(x2f[-1] - x2f[0]),
        "lz": float(x3f[-1] - x3f[0]),
        "dx1": float((x1f[-1] - x1f[0]) / nx),
        "dx2": float((x2f[-1] - x2f[0]) / ny) if ny > 1 else 1.0,
        "dx3": float((x3f[-1] - x3f[0]) / nz) if nz > 1 else 1.0,
        "dvol": float(((x1f[-1] - x1f[0]) / nx) *
                       (((x2f[-1] - x2f[0]) / ny) if ny > 1 else 1.0) *
                       (((x3f[-1] - x3f[0]) / nz) if nz > 1 else 1.0)),
    }


def _fft_mode_arrays(geometry: dict[str, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nz, ny, nx = geometry["shape"]
    n1 = np.fft.fftfreq(nx) * nx
    n2 = np.fft.fftfreq(ny) * ny if ny > 1 else np.array([0.0])
    n3 = np.fft.fftfreq(nz) * nz if nz > 1 else np.array([0.0])
    k1 = 2.0*np.pi*n1/geometry["lx"]
    k2 = 2.0*np.pi*n2/geometry["ly"] if ny > 1 else np.array([0.0])
    k3 = 2.0*np.pi*n3/geometry["lz"] if nz > 1 else np.array([0.0])
    return np.meshgrid(k3, k2, k1, indexing="ij")


def analytic_scalar_field(config: dict[str, dict[str, Any]],
                          data: dict[str, Any]) -> np.ndarray:
    geometry = _geometry_from_dump(data)
    params = _problem_defaults(config, geometry)
    dim = geometry["dim"]
    initial_field = initial_scalar_field(config, data)

    kz3, ky3, kx3 = _fft_mode_arrays(geometry)
    kmag2 = kx3*kx3 + ky3*ky3 + kz3*kz3

    velx = float(np.mean(np.asarray(data["velx"], dtype=np.float64)))
    vely = float(np.mean(np.asarray(data["vely"], dtype=np.float64)))
    velz = 0.0 if dim == 2 else float(np.mean(np.asarray(data["velz"], dtype=np.float64)))
    time = float(data["Time"])

    coeff = np.fft.fftn(initial_field).astype(np.complex128)
    coeff *= np.exp(-params["scalar_diffusivity"] * kmag2 * time)
    coeff *= np.exp(-1j * (kx3*velx + ky3*vely + kz3*velz) * time)
    return np.fft.ifftn(coeff).real


def _periodic_center(weights: np.ndarray, coords: np.ndarray,
                     xmin: float, length: float) -> float:
    theta = 2.0*np.pi*(coords - xmin)/length
    z = np.sum(weights * np.exp(1j*theta))
    angle = np.angle(z)
    if angle < 0.0:
        angle += 2.0*np.pi
    return xmin + length * angle / (2.0*np.pi)


def _periodic_delta(values: np.ndarray | float, center: float, length: float) -> np.ndarray:
    return (np.asarray(values, dtype=np.float64) - center + 0.5*length) % length - 0.5*length


def blob_moments(field: np.ndarray, data: dict[str, Any], scalar_outside: float) -> dict[str, Any]:
    weights = np.asarray(field, dtype=np.float64) - scalar_outside
    geometry = _geometry_from_dump(data)
    x1v = np.asarray(data["x1v"], dtype=np.float64)
    x2v = np.asarray(data["x2v"], dtype=np.float64)
    x3v = np.asarray(data["x3v"], dtype=np.float64)
    z3, y3, x3 = np.meshgrid(x3v, x2v, x1v, indexing="ij")

    mass = float(np.sum(weights) * geometry["dvol"])
    cx = _periodic_center(weights, x3, geometry["x1min"], geometry["lx"])
    cy = _periodic_center(weights, y3, geometry["x2min"], geometry["ly"])
    if geometry["dim"] == 3:
        cz = _periodic_center(weights, z3, geometry["x3min"], geometry["lz"])
        dz = _periodic_delta(z3, cz, geometry["lz"])
    else:
        cz = float(x3v[0])
        dz = 0.0

    dx = _periodic_delta(x3, cx, geometry["lx"])
    dy = _periodic_delta(y3, cy, geometry["ly"])
    second_moment = float(np.sum(weights * (dx*dx + dy*dy + dz*dz)) * geometry["dvol"] / mass)
    return {
        "mass": mass,
        "center": [float(cx), float(cy), float(cz)],
        "width_rms": math.sqrt(second_moment),
    }


def _uniform_velocity_stats(data: dict[str, Any]) -> dict[str, Any]:
    stats: dict[str, Any] = {}
    for key in ("velx", "vely", "velz"):
        arr = np.asarray(data[key], dtype=np.float64)
        mean = float(np.mean(arr))
        stats[key] = {"mean": mean, "max_dev": float(np.max(np.abs(arr - mean)))}
    return stats


def _evaluate_dump(input_path: str | os.PathLike[str],
                   dump_path: str | os.PathLike[str]) -> dict[str, Any]:
    config = athinput(str(input_path))
    data = read_binary_as_athdf(str(dump_path),
                                quantities=["s_00", "velx", "vely", "velz"],
                                dtype=np.float64)
    geometry = _geometry_from_dump(data)
    params = _problem_defaults(config, geometry)
    initial_field = initial_scalar_field(config, data)
    numeric = np.asarray(data["s_00"], dtype=np.float64)
    analytic = analytic_scalar_field(config, data)
    diff = numeric - analytic
    return {
        "config": config,
        "data": data,
        "geometry": geometry,
        "params": params,
        "initial_field": initial_field,
        "numeric": numeric,
        "analytic": analytic,
        "diff": diff,
    }


def analyze_dump(input_path: str | os.PathLike[str],
                 dump_path: str | os.PathLike[str]) -> dict[str, Any]:
    state = _evaluate_dump(input_path, dump_path)
    data = state["data"]
    geometry = state["geometry"]
    params = state["params"]
    initial_field = state["initial_field"]
    numeric = state["numeric"]
    analytic = state["analytic"]
    diff = state["diff"]

    l1 = float(np.mean(np.abs(diff)))
    l2 = float(np.sqrt(np.mean(diff*diff)))
    linf = float(np.max(np.abs(diff)))
    scalar_mass = float(np.sum(numeric) * geometry["dvol"])
    reference_mass = float(np.sum(initial_field) * geometry["dvol"])
    mass_abs_err = abs(scalar_mass - reference_mass)
    mass_rel_err = mass_abs_err / abs(reference_mass) if abs(reference_mass) > 0.0 else 0.0

    numeric_blob = blob_moments(numeric, data, params["scalar_outside"])
    analytic_blob = blob_moments(analytic, data, params["scalar_outside"])
    center_delta = np.array([
        _periodic_delta(numeric_blob["center"][0], analytic_blob["center"][0], geometry["lx"]),
        _periodic_delta(numeric_blob["center"][1], analytic_blob["center"][1], geometry["ly"]),
        _periodic_delta(numeric_blob["center"][2], analytic_blob["center"][2], geometry["lz"]),
    ], dtype=np.float64)
    if geometry["dim"] == 2:
        center_delta[2] = 0.0
    centroid_err = float(np.sqrt(np.sum(center_delta*center_delta)))
    width_abs_err = abs(numeric_blob["width_rms"] - analytic_blob["width_rms"])
    width_rel_err = width_abs_err / analytic_blob["width_rms"] \
        if analytic_blob["width_rms"] > 0.0 else 0.0

    metrics = {
        "dump": str(dump_path),
        "time": float(data["Time"]),
        "cycles": int(data["NumCycles"]),
        "shape": list(geometry["shape"]),
        "dim": geometry["dim"],
        "field_errors": {
            "l1": l1,
            "l2": l2,
            "linf": linf,
        },
        "mass": {
            "numeric": scalar_mass,
            "reference": reference_mass,
            "abs_error": mass_abs_err,
            "rel_error": mass_rel_err,
        },
        "blob": {
            "numeric_mass": numeric_blob["mass"],
            "analytic_mass": analytic_blob["mass"],
            "numeric_center": numeric_blob["center"],
            "analytic_center": analytic_blob["center"],
            "centroid_error": centroid_err,
            "numeric_width_rms": numeric_blob["width_rms"],
            "analytic_width_rms": analytic_blob["width_rms"],
            "width_abs_error": width_abs_err,
            "width_rel_error": width_rel_err,
        },
        "velocity": _uniform_velocity_stats(data),
    }

    thresholds = lookup_reference_thresholds(input_path, metrics)
    if thresholds is not None:
        metrics["thresholds"] = thresholds
        metrics["passes_thresholds"] = check_thresholds(metrics, thresholds)
    else:
        metrics["passes_thresholds"] = None
    return metrics


def analyze_outputs(input_path: str | os.PathLike[str],
                    dump_paths: list[str | os.PathLike[str]]) -> list[dict[str, Any]]:
    return [analyze_dump(input_path, dump_path) for dump_path in dump_paths]


def lookup_reference_thresholds(input_path: str | os.PathLike[str],
                                metrics: dict[str, Any]) -> dict[str, float] | None:
    time = float(metrics["time"])
    if abs(time - 1.0) > 1.0e-10:
        return None
    key = (Path(input_path).name, tuple(metrics["shape"]), time)
    return REFERENCE_THRESHOLDS.get(key)


def check_thresholds(metrics: dict[str, Any], thresholds: dict[str, float]) -> bool:
    return (
        metrics["field_errors"]["l1"] <= thresholds["l1_max"]
        and metrics["field_errors"]["l2"] <= thresholds["l2_max"]
        and metrics["field_errors"]["linf"] <= thresholds["linf_max"]
        and metrics["mass"]["rel_error"] <= thresholds["mass_rel_max"]
        and metrics["blob"]["centroid_error"] <= thresholds["centroid_max"]
        and metrics["blob"]["width_rel_error"] <= thresholds["width_rel_max"]
        and metrics["velocity"]["velx"]["max_dev"] <= thresholds["vel_uniform_max"]
        and metrics["velocity"]["vely"]["max_dev"] <= thresholds["vel_uniform_max"]
        and metrics["velocity"]["velz"]["max_dev"] <= thresholds["vel_uniform_max"]
    )


def _collect_dump_paths(dumps: list[str] | None, patterns: list[str] | None) -> list[str]:
    paths: list[str] = []
    if dumps:
        paths.extend(dumps)
    if patterns:
        for pattern in patterns:
            paths.extend(sorted(glob.glob(pattern)))
    deduped = sorted({str(Path(path)) for path in paths})
    if not deduped:
        raise ValueError("Specify at least one --dump or --glob path.")
    return deduped


def _import_matplotlib() -> tuple[Any, Any]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return matplotlib, plt


def _slice_plane_index(reference_field: np.ndarray, data: dict[str, Any]) -> tuple[int, str]:
    geometry = _geometry_from_dump(data)
    if geometry["dim"] == 3:
        plane_weights = np.sum(np.asarray(reference_field, dtype=np.float64), axis=(1, 2))
        k = int(np.argmax(plane_weights))
        slice_label = f"x3 = {float(np.asarray(data['x3v'])[k]):.3f}"
        return k, slice_label
    return 0, "true 2D"


def _slice_metadata(field: np.ndarray, k_index: int, slice_label: str) -> tuple[np.ndarray, str]:
    return np.asarray(field[k_index, :, :], dtype=np.float64), slice_label


def _slice_extent(data: dict[str, Any]) -> tuple[float, float, float, float]:
    x1f = np.asarray(data["x1f"], dtype=np.float64)
    x2f = np.asarray(data["x2f"], dtype=np.float64)
    return float(x1f[0]), float(x1f[-1]), float(x2f[0]), float(x2f[-1])


def _resolution_label(shape: list[int] | tuple[int, int, int]) -> str:
    nz, ny, nx = (int(val) for val in shape)
    if nz > 1:
        return f"{nx}^3"
    if ny > 1:
        return f"{nx}^2"
    return str(nx)


def _ensure_output_dir(plot_dir: str | os.PathLike[str] | None) -> Path:
    out_dir = Path(plot_dir) if plot_dir is not None else Path("scalar_mixing_test_plots")
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def plot_time_slices(input_path: str | os.PathLike[str],
                     dump_paths: list[str | os.PathLike[str]],
                     plot_dir: str | os.PathLike[str] | None = None) -> Path:
    if not dump_paths:
        raise ValueError("plot_time_slices requires at least one dump.")

    _, plt = _import_matplotlib()
    states = [_evaluate_dump(input_path, dump_path) for dump_path in dump_paths]
    ncols = len(states)
    fig, axes = plt.subplots(3, ncols, figsize=(4.6*ncols, 11.0), squeeze=False,
                             constrained_layout=True)

    scalar_min = min(float(np.min(state["numeric"])) for state in states)
    scalar_min = min(scalar_min, min(float(np.min(state["analytic"])) for state in states))
    scalar_max = max(float(np.max(state["numeric"])) for state in states)
    scalar_max = max(scalar_max, max(float(np.max(state["analytic"])) for state in states))
    err_max = max(float(np.max(np.abs(state["diff"]))) for state in states)
    err_max = max(err_max, 1.0e-16)

    scalar_im = None
    error_im = None
    for col, state in enumerate(states):
        data = state["data"]
        time = float(data["Time"])
        k_index, slice_label = _slice_plane_index(state["analytic"], data)
        numeric_slice, slice_label = _slice_metadata(state["numeric"], k_index, slice_label)
        analytic_slice, _ = _slice_metadata(state["analytic"], k_index, slice_label)
        error_slice = np.abs(numeric_slice - analytic_slice)
        extent = _slice_extent(data)

        ax = axes[0, col]
        scalar_im = ax.imshow(numeric_slice, origin="lower", extent=extent, aspect="equal",
                              cmap="viridis", vmin=scalar_min, vmax=scalar_max)
        ax.set_title(f"numeric, t={time:.3f}\n{slice_label}")
        ax.set_xlabel("x1")
        ax.set_ylabel("x2")

        ax = axes[1, col]
        ax.imshow(analytic_slice, origin="lower", extent=extent, aspect="equal",
                  cmap="viridis", vmin=scalar_min, vmax=scalar_max)
        ax.set_title(f"analytic, t={time:.3f}")
        ax.set_xlabel("x1")
        ax.set_ylabel("x2")

        ax = axes[2, col]
        error_im = ax.imshow(error_slice, origin="lower", extent=extent, aspect="equal",
                             cmap="magma", vmin=0.0, vmax=err_max)
        ax.set_title(f"|numeric - analytic|, t={time:.3f}")
        ax.set_xlabel("x1")
        ax.set_ylabel("x2")

    if scalar_im is not None:
        fig.colorbar(scalar_im, ax=axes[:2, :].ravel().tolist(), shrink=0.82, label="s_00")
    if error_im is not None:
        fig.colorbar(error_im, ax=axes[2, :].ravel().tolist(), shrink=0.82,
                     label="absolute error")

    input_name = Path(input_path).stem
    out_dir = _ensure_output_dir(plot_dir)
    out_path = out_dir / f"{input_name}_time_slices.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _resolve_case_dump(case_path: str) -> str:
    if glob.has_magic(case_path):
        matches = sorted(glob.glob(case_path))
        if not matches:
            raise ValueError(f"No dumps matched convergence pattern: {case_path}")
        return matches[-1]
    return case_path


def _convergence_groups(case_metrics: list[dict[str, Any]]) -> dict[tuple[int, float], list[dict[str, Any]]]:
    groups: dict[tuple[int, float], list[dict[str, Any]]] = {}
    for item in case_metrics:
        key = (int(item["dim"]), round(float(item["time"]), 10))
        groups.setdefault(key, []).append(item)
    for values in groups.values():
        values.sort(key=lambda item: max(int(val) for val in item["shape"]))
    return groups


def plot_convergence(case_metrics: list[dict[str, Any]],
                     plot_dir: str | os.PathLike[str] | None = None) -> list[Path]:
    if not case_metrics:
        raise ValueError("plot_convergence requires at least one analyzed case.")

    _, plt = _import_matplotlib()
    import matplotlib.ticker as mticker
    out_dir = _ensure_output_dir(plot_dir)
    output_paths: list[Path] = []

    for (dim, time), metrics_group in _convergence_groups(case_metrics).items():
        resolutions = np.asarray(
            [max(int(val) for val in item["shape"]) for item in metrics_group], dtype=np.float64
        )
        labels = [_resolution_label(item["shape"]) for item in metrics_group]
        l1 = np.asarray([item["field_errors"]["l1"] for item in metrics_group], dtype=np.float64)
        l2 = np.asarray([item["field_errors"]["l2"] for item in metrics_group], dtype=np.float64)
        linf = np.asarray([item["field_errors"]["linf"] for item in metrics_group],
                          dtype=np.float64)
        centroid = np.asarray([item["blob"]["centroid_error"] for item in metrics_group],
                              dtype=np.float64)
        width = np.asarray([item["blob"]["width_rel_error"] for item in metrics_group],
                           dtype=np.float64)
        mass = np.asarray([item["mass"]["rel_error"] for item in metrics_group], dtype=np.float64)

        fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), constrained_layout=True)

        axes[0].loglog(resolutions, l1, marker="o", linewidth=2.0, label="L1")
        axes[0].loglog(resolutions, l2, marker="s", linewidth=2.0, label="L2")
        axes[0].loglog(resolutions, linf, marker="^", linewidth=2.0, label="Linf")
        axes[0].set_xlabel("resolution")
        axes[0].set_ylabel("field error")
        axes[0].set_title(f"{dim}D field error at t={time:.3f}")
        axes[0].grid(True, which="both", alpha=0.3)
        axes[0].legend()

        axes[1].loglog(resolutions, centroid, marker="o", linewidth=2.0, label="centroid")
        axes[1].loglog(resolutions, width, marker="s", linewidth=2.0, label="width rel")
        axes[1].loglog(resolutions, np.maximum(mass, 1.0e-18), marker="^", linewidth=2.0,
                       label="mass rel")
        axes[1].set_xlabel("resolution")
        axes[1].set_ylabel("blob / conservation error")
        axes[1].set_title(f"{dim}D blob metrics at t={time:.3f}")
        axes[1].grid(True, which="both", alpha=0.3)
        axes[1].legend()

        for axis in axes:
            axis.set_xticks(resolutions)
            axis.set_xticklabels(labels)
            axis.xaxis.set_minor_formatter(mticker.NullFormatter())

        time_tag = str(time).replace(".", "p")
        out_path = out_dir / f"scalar_mixing_test_convergence_{dim}d_t{time_tag}.png"
        fig.savefig(out_path, dpi=160, bbox_inches="tight")
        plt.close(fig)
        output_paths.append(out_path)

    return output_paths


def _print_metrics(metrics: dict[str, Any]) -> None:
    print(Path(metrics["dump"]).name)
    print(f"  time={metrics['time']:.6f} cycles={metrics['cycles']} "
          f"shape={tuple(metrics['shape'])}")
    print("  field_errors:"
          f" l1={metrics['field_errors']['l1']:.6e}"
          f" l2={metrics['field_errors']['l2']:.6e}"
          f" linf={metrics['field_errors']['linf']:.6e}")
    print("  mass:"
          f" rel_error={metrics['mass']['rel_error']:.6e}"
          f" abs_error={metrics['mass']['abs_error']:.6e}")
    print("  blob:"
          f" centroid_error={metrics['blob']['centroid_error']:.6e}"
          f" width_rel_error={metrics['blob']['width_rel_error']:.6e}")
    print("  velocity_max_dev:"
          f" vx={metrics['velocity']['velx']['max_dev']:.3e}"
          f" vy={metrics['velocity']['vely']['max_dev']:.3e}"
          f" vz={metrics['velocity']['velz']['max_dev']:.3e}")
    if metrics["passes_thresholds"] is not None:
        print(f"  thresholds: {'PASS' if metrics['passes_thresholds'] else 'FAIL'}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",
                        help="AthenaK athinput file used for the run.")
    parser.add_argument("--dump", action="append",
                        help="Specific hydro_w binary dump to analyze.")
    parser.add_argument("--glob", action="append",
                        help="Glob pattern for one or more hydro_w binary dumps.")
    parser.add_argument("--json-out",
                        help="Optional JSON file for the full metric report.")
    parser.add_argument("--plot-dir",
                        help="Directory for generated PNG plots. Defaults to "
                             "./scalar_mixing_test_plots when plotting is requested.")
    parser.add_argument("--plot-slices", action="store_true",
                        help="Generate a time-sequence figure with scalar slices for the "
                             "selected --dump/--glob outputs.")
    parser.add_argument("--convergence-case", action="append", nargs=2, metavar=("INPUT", "DUMP"),
                        help="Add a case to a resolution/error convergence plot. "
                             "The second argument may be a specific dump or a glob pattern; "
                             "if it expands to multiple dumps, the latest one is used.")
    args = parser.parse_args(argv)

    metrics: list[dict[str, Any]] = []
    if args.dump or args.glob:
        if args.input is None:
            parser.error("--input is required when using --dump, --glob, or --plot-slices.")
        dump_paths = _collect_dump_paths(args.dump, args.glob)
        metrics = analyze_outputs(args.input, dump_paths)

        for item in metrics:
            _print_metrics(item)

        if args.plot_slices:
            slice_path = plot_time_slices(args.input, dump_paths, args.plot_dir)
            print(f"wrote slice plot: {slice_path}")
    elif args.plot_slices:
        parser.error("--plot-slices requires at least one --dump or --glob.")

    convergence_metrics: list[dict[str, Any]] = []
    if args.convergence_case:
        for input_path, case_path in args.convergence_case:
            resolved_dump = _resolve_case_dump(case_path)
            item = analyze_dump(input_path, resolved_dump)
            convergence_metrics.append(item)
        convergence_paths = plot_convergence(convergence_metrics, args.plot_dir)
        for out_path in convergence_paths:
            print(f"wrote convergence plot: {out_path}")

    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as fobj:
            json.dump(metrics if not convergence_metrics else {
                "selected_dumps": metrics,
                "convergence_cases": convergence_metrics,
            }, fobj, indent=2)
            fobj.write("\n")

    threshold_checks = [item["passes_thresholds"] for item in metrics + convergence_metrics
                        if item["passes_thresholds"] is not None]
    if threshold_checks and not all(threshold_checks):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
