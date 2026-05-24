#!/usr/bin/env python3
"""Compare frame-tracking runs and write validation CSV rows."""

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
PRIMITIVE_QUANTITIES = ("dens", "eint", "velx", "vely", "velz")
CONSERVED_QUANTITIES = ("dens", "mom1", "mom2", "mom3", "ener")
AXIS_INDEX = {"x1": 0, "x2": 1, "x3": 2}
VELOCITY_NAME = {"x1": "velx", "x2": "vely", "x3": "velz"}
BINARY_FIELD_TOLERANCE = 100.0 * np.finfo(np.float32).eps


def final_binary(run_dir: Path, binary_id: str | None) -> Path:
    pattern = f"*.{binary_id}.*.bin" if binary_id else "*.bin"
    paths = sorted((run_dir / "bin").glob(pattern))
    if not paths:
        paths = sorted((run_dir / "bin" / "rank_00000000").glob(pattern))
    if not paths:
        paths = sorted(run_dir.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No native binary snapshots found below {run_dir}.")
    return paths[-1]


def binary_variable_size(path: Path) -> int:
    with path.open("rb") as stream:
        for _ in range(8):
            line = stream.readline().decode("ascii")
            if "size of variable=" in line:
                return int(line.split("=")[-1])
    raise ValueError(f"Snapshot {path} lacks size-of-variable header metadata.")


def read_snapshot(run_dir: Path, quantities: tuple[str, ...],
                  binary_id: str | None) -> dict[str, Any]:
    path = final_binary(run_dir, binary_id)
    if "rank_00000000" in path.parts:
        data = bin_reader.read_all_ranks_binary_as_athdf(
            str(path), quantities=list(quantities), dtype=np.float64
        )
    else:
        data = bin_reader.read_binary_as_athdf(
            str(path), quantities=list(quantities), dtype=np.float64
        )
    missing = [name for name in quantities if name not in data]
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
        "variable_size": binary_variable_size(path),
        "fields": {name: np.asarray(data[name], dtype=float) for name in quantities},
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
    weight = np.asarray(data["ft_weight"])
    skipped = np.asarray(data["ft_skip"])
    sampled = np.flatnonzero(weight > 0.0)
    priming_index = int(sampled[0]) if sampled.size else -1
    priming_skips = skipped[priming_index] if priming_index >= 0 else 0.0
    values["max_misses"] = float(np.max(np.asarray(data["ft_misses"])))
    values["recovery_events"] = float(np.sum(np.asarray(data["ft_recov"])))
    values["limit_events"] = float(np.sum(np.asarray(data["ft_limit"])))
    values["priming_skip_events"] = float(priming_skips)
    values["unexpected_skip_events"] = float(np.sum(skipped) - priming_skips)
    return values


def axis_centers(snapshot: dict[str, Any], axis: str, displacement: float) -> np.ndarray:
    return np.asarray(snapshot[f"{axis}v"]) + displacement


def selected_mask(snapshot: dict[str, Any], args: argparse.Namespace) -> np.ndarray:
    density = snapshot["fields"]["dens"]
    if args.selection == "density":
        value = density
    elif args.selection == "temperature":
        value = (args.gamma - 1.0) * snapshot["fields"]["eint"] / density
    else:
        value = snapshot["fields"][args.scalar_field]
    return (value >= args.target_min) & (value <= args.target_max)


def selected_observables(snapshot: dict[str, Any], history: dict[str, float],
                         args: argparse.Namespace) -> tuple[float, float]:
    density = snapshot["fields"]["dens"]
    spacings = [float(np.mean(np.diff(snapshot[f"x{n}f"]))) for n in (1, 2, 3)]
    if args.selection == "scalar":
        scalar = np.maximum(snapshot["fields"][args.scalar_field], 0.0)
        weights = density * scalar * np.prod(spacings)
    else:
        mask = selected_mask(snapshot, args)
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


def hydro_conserved_lab_frame(snapshot: dict[str, Any], history: dict[str, float],
                              axis: str) -> dict[str, np.ndarray]:
    density = snapshot["fields"]["dens"]
    velocity = {
        name: transformed_field(snapshot, history, name, axis)
        for name in ("velx", "vely", "velz")
    }
    return {
        "dens": density,
        "mom1": density * velocity["velx"],
        "mom2": density * velocity["vely"],
        "mom3": density * velocity["velz"],
        "ener": (
            snapshot["fields"]["eint"]
            + 0.5
            * density
            * sum(velocity[name] * velocity[name] for name in velocity)
        ),
    }


def interpolate_to_reference(reference: dict[str, Any], candidate: dict[str, Any],
                             reference_history: dict[str, float],
                             candidate_history: dict[str, float], field: str,
                             axis: str) -> tuple[np.ndarray, np.ndarray]:
    ref_values = transformed_field(reference, reference_history, field, axis)
    cand_values = transformed_field(candidate, candidate_history, field, axis)
    return interpolate_values_to_reference(
        reference, candidate, reference_history, candidate_history,
        ref_values, cand_values, axis
    )


def interpolate_values_to_reference(
    reference: dict[str, Any],
    candidate: dict[str, Any],
    reference_history: dict[str, float],
    candidate_history: dict[str, float],
    ref_values: np.ndarray,
    cand_values: np.ndarray,
    axis: str,
) -> tuple[np.ndarray, np.ndarray]:
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


def grids_match(reference: dict[str, Any], candidate: dict[str, Any]) -> bool:
    return all(
        np.array_equal(reference[coordinate], candidate[coordinate])
        for coordinate in ("x1v", "x2v", "x3v", "x1f", "x2f", "x3f")
    )


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
    if args.selection == "scalar":
        scalar = candidate["fields"][args.scalar_field]
        density = candidate["fields"]["dens"]
        dvol = np.prod(
            [float(np.mean(np.diff(candidate[f"x{n}f"]))) for n in (1, 2, 3)]
        )
        excess = float(
            np.sum(
                density
                * (np.maximum(-scalar, 0.0) + np.maximum(scalar - 1.0, 0.0))
                * dvol
            )
        )
        relative_excess = excess / max(candidate_mass, np.finfo(float).tiny)
        add_row(rows, args, "tracer_out_of_range_mass", excess, relative_excess,
                "relative <= 1.0e-8",
                "pass" if relative_excess <= 1.0e-8 else "fail")
        ref_values, cand_values = interpolate_values_to_reference(
            reference,
            candidate,
            reference_history,
            candidate_history,
            reference["fields"]["dens"] * reference["fields"][args.scalar_field],
            candidate["fields"]["dens"] * candidate["fields"][args.scalar_field],
            args.axis,
        )
        absolute, relative = norm_error(ref_values, cand_values)
        add_row(rows, args, "tracer_density_l2_lab_frame", absolute, relative,
                "trend from low to medium resolution", "measured")
    reference_conserved = hydro_conserved_lab_frame(
        reference, reference_history, args.axis
    )
    candidate_conserved = hydro_conserved_lab_frame(
        candidate, candidate_history, args.axis
    )
    ref_state_values = []
    cand_state_values = []
    for name in ("dens", "mom1", "mom2", "mom3", "ener"):
        ref_values, cand_values = interpolate_values_to_reference(
            reference,
            candidate,
            reference_history,
            candidate_history,
            reference_conserved[name],
            candidate_conserved[name],
            args.axis,
        )
        ref_state_values.append(ref_values)
        cand_state_values.append(cand_values)
        absolute, relative = norm_error(ref_values, cand_values)
        add_row(rows, args, f"{name}_l2_lab_frame", absolute, relative,
                "trend from low to medium resolution", "measured")
    absolute, relative = norm_error(
        np.concatenate(ref_state_values), np.concatenate(cand_state_values)
    )
    add_row(rows, args, "hydro_conserved_l2_lab_frame", absolute, relative,
            "trend from low to medium resolution", "measured")


def compare_continuity(rows: list[dict[str, object]], reference: dict[str, Any],
                       candidate: dict[str, Any], reference_history: dict[str, float],
                       candidate_history: dict[str, float],
                       args: argparse.Namespace) -> None:
    field_tolerance = args.tolerance if args.require_real_output else args.field_tolerance
    common_grid = grids_match(reference, candidate)
    precision = "real" if candidate["variable_size"] == 8 else "float32"
    fields = list(args.state_quantities)
    if args.selection == "scalar":
        fields.append(args.scalar_field)
    for field in fields:
        if args.comparison_kind in ("restart", "mpi") and common_grid:
            ref_values = reference["fields"][field]
            cand_values = candidate["fields"][field]
            metric = f"{field}_l2_grid_frame_binary_{precision}"
        else:
            ref_values, cand_values = interpolate_to_reference(
                reference, candidate, reference_history, candidate_history,
                field, args.axis
            )
            metric = f"{field}_l2_lab_frame_binary_{precision}"
        absolute, relative = norm_error(ref_values, cand_values)
        result = "pass" if (absolute <= field_tolerance or
                            relative <= field_tolerance) else "fail"
        criterion = (
            f"native {precision} snapshot: absolute or relative <= {field_tolerance:.16e}"
        )
        add_row(rows, args, metric, absolute, relative, criterion, result)
    for key in (f"ft_vf_{args.axis}", f"ft_dx_{args.axis}", f"ft_pos_{args.axis}",
                f"ft_err_{args.axis}", f"ft_dv_{args.axis}"):
        if key not in reference_history or key not in candidate_history:
            continue
        absolute, relative = scaled_error(candidate_history[key], reference_history[key])
        scale = args.tolerance * max(1.0, abs(reference_history[key]))
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
    if "recovery_events" in history:
        recovered = history["recovery_events"]
        add_row(rows, args, "recovery_events", recovered, recovered, "absolute <= 0",
                "pass" if recovered == 0.0 else "fail")
    if "priming_skip_events" in history:
        priming = history["priming_skip_events"]
        add_row(rows, args, "priming_skip_events", priming, priming,
                "expected initialization event <= 1",
                "pass" if priming <= 1.0 else "fail")
    if "unexpected_skip_events" in history:
        skipped = history["unexpected_skip_events"]
        add_row(rows, args, "unexpected_skip_events", skipped, skipped,
                "absolute <= 0", "pass" if skipped == 0.0 else "fail")


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
    result.add_argument("--binary-id", default=None,
                        help="Optional binary output id, for example hydro_w or hydro_u.")
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
    result.add_argument("--selection", choices=("density", "temperature", "scalar"),
                        required=True)
    result.add_argument("--scalar-field", default="s_00")
    result.add_argument(
        "--state", choices=("primitive", "conserved"), default="primitive"
    )
    result.add_argument("--require-real-output", action="store_true",
                        help="Require native Real binary fields for strict comparisons.")
    result.add_argument("--target-min", type=float, required=True)
    result.add_argument("--target-max", type=float, required=True)
    result.add_argument("--gamma", type=float, default=5.0 / 3.0)
    result.add_argument("--tolerance", type=float, default=1.0e-10)
    result.add_argument(
        "--field-tolerance",
        type=float,
        default=BINARY_FIELD_TOLERANCE,
        help=("Tolerance for native binary snapshot fields; defaults to 100 times "
              "float32 epsilon because AthenaK binary outputs store fields as float."),
    )
    return result


def main() -> None:
    args = parser().parse_args()
    if args.require_real_output and args.comparison_kind == "physical":
        raise ValueError("--require-real-output is intended for continuity comparisons.")
    if (args.comparison_kind != "physical" and not args.require_real_output and
            args.field_tolerance < BINARY_FIELD_TOLERANCE):
        raise ValueError(
            "Native binary snapshots store fields as float32 and cannot certify a "
            f"field tolerance below {BINARY_FIELD_TOLERANCE:.16e}; use high-precision "
            "conserved-state output for strict restart certification."
        )
    args.state_quantities = (
        PRIMITIVE_QUANTITIES if args.state == "primitive" else CONSERVED_QUANTITIES
    )
    quantities = args.state_quantities + (
        (args.scalar_field,) if args.selection == "scalar" else ()
    )
    reference = read_snapshot(args.reference_dir, quantities, args.binary_id)
    candidate = read_snapshot(args.candidate_dir, quantities, args.binary_id)
    if reference["variable_size"] != candidate["variable_size"]:
        raise ValueError(
            "Reference and candidate binary outputs use different precisions."
        )
    if args.require_real_output and reference["variable_size"] != 8:
        raise ValueError(
            "Strict continuity comparison requires data_precision=real in a "
            "double-precision build."
        )
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
