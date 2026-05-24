#!/usr/bin/env python3
"""Localize lab-frame differences in tracked versus untracked cloud runs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import compare_frame_tracking_validation as validation


FIELDS = (
    "resolution",
    "diagnostic_metric",
    "region",
    "absolute_error",
    "relative_error",
    "squared_error_fraction",
    "notes",
)
STATE_FIELDS = ("dens", "mom1", "mom2", "mom3", "ener")


def aligned_values(
    reference: dict[str, object],
    candidate: dict[str, object],
    reference_history: dict[str, float],
    candidate_history: dict[str, float],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], np.ndarray]:
    reference_state = validation.hydro_conserved_lab_frame(
        reference, reference_history, "x1"
    )
    candidate_state = validation.hydro_conserved_lab_frame(
        candidate, candidate_history, "x1"
    )
    reference_x = validation.axis_centers(
        reference, "x1", reference_history["displacement"]
    )
    candidate_x = validation.axis_centers(
        candidate, "x1", candidate_history["displacement"]
    )
    valid = (reference_x >= candidate_x[0]) & (reference_x <= candidate_x[-1])

    def interpolate(values: np.ndarray) -> np.ndarray:
        lines = np.moveaxis(values, -1, -1).reshape(-1, candidate_x.size)
        interpolated = np.vstack(
            [
                np.interp(reference_x, candidate_x, line, left=np.nan, right=np.nan)
                for line in lines
            ]
        )
        shaped = interpolated.reshape(values.shape[:-1] + (reference_x.size,))
        return shaped[..., valid]

    reference_aligned = {
        name: reference_state[name][..., valid] for name in STATE_FIELDS
    }
    candidate_aligned = {
        name: interpolate(candidate_state[name]) for name in STATE_FIELDS
    }
    return reference_aligned, candidate_aligned, reference_x[valid]


def add_row(
    rows: list[dict[str, str]],
    resolution: str,
    metric: str,
    region: str,
    absolute: float,
    relative: float,
    fraction: float,
    notes: str,
) -> None:
    rows.append(
        {
            "resolution": resolution,
            "diagnostic_metric": metric,
            "region": region,
            "absolute_error": f"{absolute:.16e}",
            "relative_error": f"{relative:.16e}",
            "squared_error_fraction": f"{fraction:.16e}",
            "notes": notes,
        }
    )


def norm_row(
    rows: list[dict[str, str]],
    resolution: str,
    metric: str,
    region: str,
    reference: np.ndarray,
    candidate: np.ndarray,
    mask: np.ndarray,
    denominator_error: float,
    notes: str,
) -> None:
    difference = candidate[mask] - reference[mask]
    absolute = float(np.linalg.norm(difference))
    denominator = max(float(np.linalg.norm(reference[mask])), np.finfo(float).tiny)
    fraction = float(np.sum(difference * difference)) / max(
        denominator_error, np.finfo(float).tiny
    )
    add_row(
        rows, resolution, metric, region, absolute, absolute / denominator,
        fraction, notes
    )


def analyze_pair(
    resolution: str,
    reference_dir: Path,
    candidate_dir: Path,
    tracer_threshold: float,
) -> tuple[list[dict[str, str]], dict[str, np.ndarray]]:
    quantities = validation.PRIMITIVE_QUANTITIES + ("s_00",)
    reference = validation.read_snapshot(reference_dir, quantities, "hydro_w")
    candidate = validation.read_snapshot(candidate_dir, quantities, "hydro_w")
    reference_history = validation.final_history(reference_dir, "x1")
    candidate_history = validation.final_history(candidate_dir, "x1")
    reference_state, candidate_state, x1 = aligned_values(
        reference, candidate, reference_history, candidate_history
    )
    difference2 = sum(
        (candidate_state[name] - reference_state[name]) ** 2 for name in STATE_FIELDS
    )
    reference2 = sum(reference_state[name] ** 2 for name in STATE_FIELDS)
    total_error = float(np.sum(difference2))
    tracer_reference = (
        reference["fields"]["dens"] * np.maximum(reference["fields"]["s_00"], 0.0)
    )
    tracer_candidate = (
        candidate["fields"]["dens"] * np.maximum(candidate["fields"]["s_00"], 0.0)
    )
    reference_x = validation.axis_centers(
        reference, "x1", reference_history["displacement"]
    )
    candidate_x = validation.axis_centers(
        candidate, "x1", candidate_history["displacement"]
    )
    valid = (reference_x >= candidate_x[0]) & (reference_x <= candidate_x[-1])
    tracer_reference = tracer_reference[..., valid]
    tracer_lines = tracer_candidate.reshape(-1, candidate_x.size)
    tracer_interpolated = np.vstack(
        [
            np.interp(reference_x, candidate_x, line, left=np.nan, right=np.nan)
            for line in tracer_lines
        ]
    ).reshape(tracer_candidate.shape)[..., valid]

    xgrid = x1.reshape((1, 1, x1.size))
    regions = {
        "full_overlap": np.ones_like(difference2, dtype=bool),
        "inflow_x1_lt_0": np.broadcast_to(xgrid < 0.0, difference2.shape),
        "interaction_1_le_x1_lt_4": np.broadcast_to(
            (xgrid >= 1.0) & (xgrid < 4.0), difference2.shape
        ),
        "downstream_x1_ge_4": np.broadcast_to(xgrid >= 4.0, difference2.shape),
        "tracer_bearing": (tracer_reference > tracer_threshold)
        | (tracer_interpolated > tracer_threshold),
    }
    regions["nontracer"] = ~regions["tracer_bearing"]

    rows: list[dict[str, str]] = []
    for region, mask in regions.items():
        absolute = float(np.sqrt(np.sum(difference2[mask])))
        relative = absolute / max(float(np.sqrt(np.sum(reference2[mask]))),
                                  np.finfo(float).tiny)
        fraction = float(np.sum(difference2[mask])) / max(
            total_error, np.finfo(float).tiny
        )
        add_row(
            rows,
            resolution,
            "hydro_conserved_l2_lab_frame",
            region,
            absolute,
            relative,
            fraction,
            "Aggregate state norm; regions use fixed lab-frame x1 bounds.",
        )

    interaction_mask = regions["interaction_1_le_x1_lt_4"]
    for name in STATE_FIELDS:
        field_error = float(np.sum((candidate_state[name] - reference_state[name]) ** 2))
        norm_row(
            rows,
            resolution,
            f"{name}_l2_lab_frame",
            "interaction_1_le_x1_lt_4",
            reference_state[name],
            candidate_state[name],
            interaction_mask,
            field_error,
            "Per-field norm within the shock/cloud interaction zone.",
        )

    spacings = [
        float(np.mean(np.diff(reference[f"x{n}f"]))) for n in (1, 2, 3)
    ]
    dvol = float(np.prod(spacings))
    for name in STATE_FIELDS:
        reference_integral = float(np.sum(reference_state[name]) * dvol)
        difference = float(
            np.sum(candidate_state[name] - reference_state[name]) * dvol
        )
        relative = abs(difference) / max(abs(reference_integral), np.finfo(float).tiny)
        add_row(
            rows,
            resolution,
            f"integral_{name}_lab_frame",
            "full_overlap",
            abs(difference),
            relative,
            0.0,
            "Absolute and relative integrated-state difference over overlapping domain.",
        )

    reference_eint = reference["fields"]["eint"][..., valid]
    candidate_eint_lines = candidate["fields"]["eint"].reshape(-1, candidate_x.size)
    candidate_eint = np.vstack(
        [
            np.interp(reference_x, candidate_x, line, left=np.nan, right=np.nan)
            for line in candidate_eint_lines
        ]
    ).reshape(candidate["fields"]["eint"].shape)[..., valid]
    reference_profile = np.mean(reference_eint, axis=(0, 1))
    candidate_profile = np.mean(candidate_eint, axis=(0, 1))
    search = (x1 >= 1.0) & (x1 < 4.0)
    reference_gradient = np.abs(np.gradient(reference_profile, x1))
    candidate_gradient = np.abs(np.gradient(candidate_profile, x1))
    reference_peak = float(
        x1[np.flatnonzero(search)[np.argmax(reference_gradient[search])]]
    )
    candidate_peak = float(
        x1[np.flatnonzero(search)[np.argmax(candidate_gradient[search])]]
    )
    add_row(
        rows,
        resolution,
        "mean_eint_gradient_peak_position",
        "interaction_1_le_x1_lt_4",
        abs(candidate_peak - reference_peak),
        abs(candidate_peak - reference_peak) / spacings[0],
        0.0,
        f"disabled={reference_peak:.16e}; enabled={candidate_peak:.16e}; "
        "relative in cells.",
    )

    squared_by_x1 = np.sum(difference2, axis=(0, 1))
    return rows, {
        "x1": x1,
        "error_fraction_by_x1": squared_by_x1 / max(total_error, np.finfo(float).tiny),
        "reference_eint_profile": reference_profile,
        "candidate_eint_profile": candidate_profile,
    }


def write_plot(path: Path, data: dict[str, dict[str, np.ndarray]]) -> None:
    figure, axes = plt.subplots(
        2, 1, figsize=(8, 7), sharex=True, constrained_layout=True
    )
    for resolution, values in data.items():
        axes[0].semilogy(
            values["x1"],
            np.maximum(values["error_fraction_by_x1"], 1.0e-20),
            marker="o",
            label=resolution,
        )
        axes[1].plot(values["x1"], values["reference_eint_profile"],
                     label=f"{resolution}, disabled")
        axes[1].plot(values["x1"], values["candidate_eint_profile"], linestyle="--",
                     label=f"{resolution}, enabled")
    axes[0].axvspan(1.0, 4.0, color="0.9", label="interaction zone")
    axes[0].set_ylabel("squared error fraction per x1 column")
    axes[0].legend()
    axes[1].axvspan(1.0, 4.0, color="0.9")
    axes[1].set_xlabel("lab x1")
    axes[1].set_ylabel("transverse mean internal energy")
    axes[1].legend(ncol=2)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--low-reference-dir", type=Path, required=True)
    parser.add_argument("--low-candidate-dir", type=Path, required=True)
    parser.add_argument("--medium-reference-dir", type=Path, required=True)
    parser.add_argument("--medium-candidate-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot-output", type=Path, required=True)
    parser.add_argument("--tracer-threshold", type=float, default=1.0e-8)
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    plot_data: dict[str, dict[str, np.ndarray]] = {}
    for resolution, reference_dir, candidate_dir in (
        ("48x16x16", args.low_reference_dir, args.low_candidate_dir),
        ("96x32x32", args.medium_reference_dir, args.medium_candidate_dir),
    ):
        pair_rows, pair_plot = analyze_pair(
            resolution, reference_dir, candidate_dir, args.tracer_threshold
        )
        rows.extend(pair_rows)
        plot_data[resolution] = pair_plot

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_plot(args.plot_output, plot_data)
    print(f"Wrote {len(rows)} diagnostic row(s) to {args.output}")
    print(f"Wrote diagnostic figure to {args.plot_output}")


if __name__ == "__main__":
    main()
