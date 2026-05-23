#!/usr/bin/env python3
"""Compare frame-tracking runs in lab coordinates and write validation CSV rows."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
try:
    import bin_convert_new as bin_reader  # type: ignore[import-not-found]  # noqa: E402
except ModuleNotFoundError:
    import bin_convert as bin_reader  # type: ignore[import-not-found]  # noqa: E402
import athena_read  # type: ignore[import-not-found]  # noqa: E402


FIELDS = (
    "problem",
    "resolution",
    "tracking_mode",
    "restart_split",
    "ranks",
    "amr_mode",
    "diagnostic_metric",
    "comparison_reference",
    "absolute_error",
    "relative_error",
    "tolerance",
    "result",
)
QUANTITIES = ("dens", "eint", "velx", "vely", "velz")
AXIS_INDEX = {"x1": 0, "x2": 1, "x3": 2}
VELOCITY_NAME = {"x1": "velx", "x2": "vely", "x3": "velz"}


def final_binary(run_dir: Path) -> Path:
    paths = sorted((run_dir / "bin").glob("*.bin"))
    if not paths:
        paths = sorted((run_dir / "bin" / "rank_00000000").glob("*.bin"))
    if not paths:
        paths = sorted(run_dir.glob("*.bin"))
    if not paths:
        raise FileNotFoundError(f"No native binary snapshots found below {run_dir}.")
    return paths[-1]


def read_snapshot(run_dir: Path) -> dict[str, Any]:
    path = final_binary(run_dir)
    if "rank_00000000" in path.parts:
        data = bin_reader.read_all_ranks_binary_as_athdf(
            str(path), quantities=list(QUANTITIES)
        )
    else:
        data = bin_reader.read_binary_as_athdf(str(path), quantities=list(QUANTITIES))
    missing = [name for name in QUANTITIES if name not in data]
    if missing:
        raise ValueError(f"Snapshot {path} lacks required fields: {', '.join(missing)}.")
    return {
        "time": float(data["Time"]),
        "x1v": np.asarray(data["x1v"], dtype=float),
        "x2v": np.asarray(data["x2v"], dtype=float),
        "x3v": np.asarray(data["x3v"], dtype=float),
        "x1f": np.asarray(data["x1f"], dtype=float),
        "x2f": np.asarray(data["x2f"], dtype=float),
        "x3f": np.asarray(data["x3f"], dtype=float),
        "fields": {name: np.asarray(data[name], dtype=float) for name in QUANTITIES},
    }


def final_history(run_dir: Path, axis: str) -> dict[str, float]:
    paths = sorted(run_dir.glob("*.frame_tracker.hst"))
    if not paths:
        return {"displacement": 0.0, "velocity": 0.0}
    data = athena_read.hst(str(paths[0]))
    suffix = axis
    values: dict[str, float] = {
        "displacement": float(np.asarray(data[f"ft_dx_{suffix}"])[-1]),
        "velocity": float(np.asarray(data[f"ft_vf_{suffix}"])[-1]),
    }
    for field in ("ft_weight", "ft_misses", "ft_recov", "ft_limit", "ft_skip",
                  f"ft_vf_{suffix}", f"ft_dx_{suffix}", f"ft_pos_{suffix}",
                  f"ft_err_{suffix}", f"ft_dv_{suffix}"):
        if field in data:
            values[field] = float(np.asarray(data[field])[-1])
    values["max_misses"] = float(np.max(np.asarray(data["ft_misses"])))
    values["limit_events"] = float(np.sum(np.asarray(data["ft_limit"])))
    return values


def axis_centers(snapshot: dict[str, Any], axis: str, displacement: float) -> np.ndarray:
    return np.asarray(snapshot[f"{axis}v"]) + displacement


def selected_mask(snapshot: dict[str, Any], args: argparse.Namespace) -> np.ndarray:
    density = snapshot["fields"]["dens"]
    if args.selection == "density":
        value = density
    else:
        value = (args.gamma - 1.0) * snapshot["fields"]["eint"] / density
    return (value >= args.target_min) & (value <= args.target_max)


def selected_observables(snapshot: dict[str, Any], history: dict[str, float],
                         args: argparse.Namespace) -> tuple[float, float]:
    density = snapshot["fields"]["dens"]
    mask = selected_mask(snapshot, args)
    spacings = [float(np.mean(np.diff(snapshot[f"x{n}f"]))) for n in (1, 2, 3)]
    weights = np.where(mask, density * np.prod(spacings), 0.0)
    mass = float(np.sum(weights))
    if mass <= 0.0:
        return mass, math.nan
    coordinate = axis_centers(snapshot, args.axis, history["displacement"])
    shape = [1, 1, 1]
    shape[2 - AXIS_INDEX[args.axis]] = coordinate.size
    coordinate_grid = coordinate.reshape(shape)
    centroid = float(np.sum(weights * coordinate_grid) / mass)
    return mass, centroid


def transformed_field(snapshot: dict[str, Any], history: dict[str, float],
                      name: str, axis: str) -> np.ndarray:
    value = np.asarray(snapshot["fields"][name], dtype=float)
    if name == VELOCITY_NAME[axis]:
        value = value + history["velocity"]
    return value


def interpolate_to_reference(reference: dict[str, Any], candidate: dict[str, Any],
                             reference_history: dict[str, float],
                             candidate_history: dict[str, float], field: str,
                             axis: str) -> tuple[np.ndarray, np.ndarray]:
    ref_values = transformed_field(reference, reference_history, field, axis)
    cand_values = transformed_field(candidate, candidate_history, field, axis)
    dimension = 2 - AXIS_INDEX[axis]
    ref_x = axis_centers(reference, axis, reference_history["displacement"])
    cand_x = axis_centers(candidate, axis, candidate_history["displacement"])
    ref_lines = np.moveaxis(ref_values, dimension, -1).reshape(-1, ref_x.size)
    cand_lines = np.moveaxis(cand_values, dimension, -1).reshape(-1, cand_x.size)
    if ref_lines.shape[0] != cand_lines.shape[0]:
        raise ValueError("Reference and candidate transverse grids do not match.")
    interpolated = np.empty_like(ref_lines)
    for line, cand_line in enumerate(cand_lines):
        interpolated[line] = np.interp(ref_x, cand_x, cand_line,
                                       left=np.nan, right=np.nan)
    valid = np.isfinite(interpolated)
    return ref_lines[valid], interpolated[valid]


def scaled_error(candidate: float, reference: float) -> tuple[float, float]:
    absolute = abs(candidate - reference)
    relative = absolute / max(abs(reference), np.finfo(float).tiny)
    return absolute, relative


def norm_error(reference: np.ndarray, candidate: np.ndarray) -> tuple[float, float]:
    difference = candidate - reference
    absolute = float(np.linalg.norm(difference.ravel()))
    denominator = max(float(np.linalg.norm(reference.ravel())), np.finfo(float).tiny)
    return absolute, absolute / denominator


def add_row(rows: list[dict[str, object]], args: argparse.Namespace, metric: str,
            absolute: float, relative: float, tolerance: str,
            result: str) -> None:
    rows.append({
        "problem": args.problem,
        "resolution": args.resolution,
        "tracking_mode": args.tracking_mode,
        "restart_split": args.restart_split,
        "ranks": args.ranks,
        "amr_mode": args.amr_mode,
        "diagnostic_metric": metric,
        "comparison_reference": args.comparison_reference,
        "absolute_error": f"{absolute:.16e}",
        "relative_error": f"{relative:.16e}",
        "tolerance": tolerance,
        "result": result,
    })


def compare_physical(rows: list[dict[str, object]], reference: dict[str, Any],
                     candidate: dict[str, Any], reference_history: dict[str, float],
                     candidate_history: dict[str, float],
                     args: argparse.Namespace) -> None:
    reference_mass, reference_centroid = selected_observables(
        reference, reference_history, args
    )
    candidate_mass, candidate_centroid = selected_observables(
        candidate, candidate_history, args
    )
    absolute, relative = scaled_error(candidate_mass, reference_mass)
    add_row(rows, args, "selected_mass", absolute, relative, "relative <= 1.0e-2",
            "pass" if relative <= 1.0e-2 else "fail")
    absolute, relative = scaled_error(candidate_centroid, reference_centroid)
    cell = float(np.mean(np.diff(reference[f"{args.axis}f"])))
    add_row(rows, args, "tracked_centroid", absolute, relative,
            f"absolute <= {cell:.16e}", "pass" if absolute <= cell else "fail")
    ref_values, cand_values = interpolate_to_reference(
        reference, candidate, reference_history, candidate_history, "dens", args.axis
    )
    absolute, relative = norm_error(ref_values, cand_values)
    add_row(rows, args, "density_l2_lab_frame", absolute, relative,
            "trend against wiring baseline", "measured")


def compare_continuity(rows: list[dict[str, object]], reference: dict[str, Any],
                       candidate: dict[str, Any], reference_history: dict[str, float],
                       candidate_history: dict[str, float],
                       args: argparse.Namespace) -> None:
    tolerance = args.tolerance
    for field in QUANTITIES:
        ref_values, cand_values = interpolate_to_reference(
            reference, candidate, reference_history, candidate_history, field, args.axis
        )
        absolute, relative = norm_error(ref_values, cand_values)
        result = "pass" if (absolute <= tolerance or relative <= tolerance) else "fail"
        add_row(rows, args, f"{field}_l2_lab_frame", absolute, relative,
                f"absolute or relative <= {tolerance:.16e}", result)
    for key in (f"ft_vf_{args.axis}", f"ft_dx_{args.axis}", f"ft_pos_{args.axis}",
                f"ft_err_{args.axis}", f"ft_dv_{args.axis}"):
        if key not in reference_history or key not in candidate_history:
            continue
        absolute, relative = scaled_error(candidate_history[key], reference_history[key])
        scale = tolerance * max(1.0, abs(reference_history[key]))
        add_row(rows, args, key, absolute, relative,
                f"scaled absolute <= {scale:.16e}",
                "pass" if absolute <= scale else "fail")


def add_health_rows(rows: list[dict[str, object]], candidate: dict[str, Any],
                    history: dict[str, float], args: argparse.Namespace) -> None:
    finite = all(np.isfinite(value).all() for value in candidate["fields"].values())
    add_row(rows, args, "finite_fields", 0.0 if finite else 1.0,
            0.0 if finite else 1.0, "all finite", "pass" if finite else "fail")
    if "max_misses" in history:
        misses = history["max_misses"]
        add_row(rows, args, "max_missed_samples", misses, misses, "absolute <= 0",
                "pass" if misses == 0.0 else "fail")
    if "limit_events" in history:
        limited = history["limit_events"]
        add_row(rows, args, "limit_events", limited, limited, "absolute <= 0",
                "pass" if limited == 0.0 else "fail")


def write_rows(path: Path, rows: list[dict[str, object]], append: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not append or not path.exists()
    mode = "a" if append else "w"
    with path.open(mode, newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--reference-dir", type=Path, required=True)
    result.add_argument("--candidate-dir", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--append", action="store_true")
    result.add_argument("--problem", required=True)
    result.add_argument("--resolution", required=True)
    result.add_argument("--tracking-mode", required=True)
    result.add_argument("--restart-split", default="none")
    result.add_argument("--ranks", type=int, default=1)
    result.add_argument("--amr-mode", default="uniform")
    result.add_argument("--comparison-reference", required=True)
    result.add_argument("--comparison-kind",
                        choices=("physical", "restart", "mpi", "amr"),
                        default="physical")
    result.add_argument("--axis", choices=tuple(AXIS_INDEX), required=True)
    result.add_argument("--selection", choices=("density", "temperature"), required=True)
    result.add_argument("--target-min", type=float, required=True)
    result.add_argument("--target-max", type=float, required=True)
    result.add_argument("--gamma", type=float, default=5.0 / 3.0)
    result.add_argument("--tolerance", type=float, default=1.0e-10)
    return result


def main() -> None:
    args = parser().parse_args()
    reference = read_snapshot(args.reference_dir)
    candidate = read_snapshot(args.candidate_dir)
    reference_history = final_history(args.reference_dir, args.axis)
    candidate_history = final_history(args.candidate_dir, args.axis)
    rows: list[dict[str, object]] = []
    if args.comparison_kind == "physical":
        compare_physical(rows, reference, candidate, reference_history,
                         candidate_history, args)
    else:
        compare_continuity(rows, reference, candidate, reference_history,
                           candidate_history, args)
    add_health_rows(rows, candidate, candidate_history, args)
    write_rows(args.output, rows, args.append)
    print(f"Wrote {len(rows)} validation row(s) to {args.output}")


if __name__ == "__main__":
    main()
