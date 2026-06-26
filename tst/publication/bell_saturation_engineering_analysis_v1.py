#!/usr/bin/env python3
"""Reduce a compact 3D Bell run and apply a sustained-saturation gate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
from typing import Any, Sequence

import numpy as np

from tst.publication.analyze_q011_section54_outputs import (
    AthenaBinaryDataset,
    compose_leaf_field,
    read_athenak_binary,
)


K0 = 2.0 * math.pi
UA0 = 1.0
B0 = 1.0
INITIAL_JX = 4.0 * math.pi
ONSET_BPERP_OVER_B0 = 1.0
MINIMUM_CANDIDATE_SPAN_TAU = 4.0
MAXIMUM_CANDIDATE_SPAN_TAU = 8.0
MINIMUM_CONFIRMATION_SPAN_TAU = 8.0
MINIMUM_CONFIRMATION_SAMPLES = 12
MAXIMUM_CANDIDATE_LOG_ENERGY_SLOPE = 0.10
MAXIMUM_CONFIRMATION_LOG_ENERGY_SLOPE = 0.05
MAXIMUM_ENERGY_MAX_TO_MIN = 1.50
MINIMUM_CURRENT_RELAXATION_FRACTION = 0.10
MINIMUM_DENSITY_RMS_OVER_MEAN = 0.05
MINIMUM_FILAMENT_VOLUME_FRACTION = 1.0e-3
_INDEX_RE = re.compile(r"\.([0-9]{5})\.bin$")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _index(path: Path) -> int:
    match = _INDEX_RE.search(path.name)
    _require(match is not None, f"cannot parse output index: {path}")
    return int(match.group(1))


def _field(dataset: AthenaBinaryDataset, name: str) -> np.ndarray:
    return np.asarray(compose_leaf_field(dataset, name).values, dtype=np.float64)


def _load_product(path: Path, name: str) -> tuple[AthenaBinaryDataset, np.ndarray]:
    dataset = read_athenak_binary(path)
    _require(dataset.variable_names == (name,), f"unexpected field inventory in {path}")
    return dataset, _field(dataset, name)


def _slope(x: Sequence[float], y: Sequence[float]) -> float:
    coordinates = np.asarray(x, dtype=np.float64)
    values = np.asarray(y, dtype=np.float64)
    centered = coordinates - np.mean(coordinates)
    denominator = float(np.sum(centered * centered))
    _require(denominator > 0.0, "slope interval has no time variance")
    return float(np.sum(centered * (values - np.mean(values))) / denominator)


def _window_metrics(rows: Sequence[dict[str, float]], start: int, end: int) -> dict[str, float]:
    selected = rows[start : end + 1]
    tau = [row["tau"] for row in selected]
    energy = [row["bperp_rms_over_B0"] ** 2 for row in selected]
    slope = _slope(tau, np.log(energy))
    return {
        "start_index": start,
        "end_index": end,
        "start_tau": tau[0],
        "end_tau": tau[-1],
        "span_tau": tau[-1] - tau[0],
        "sample_count": len(selected),
        "log_energy_slope_per_tau": slope,
        "energy_max_to_min": max(energy) / min(energy),
        "median_bperp_rms_over_B0": float(
            np.median([row["bperp_rms_over_B0"] for row in selected])
        ),
        "terminal_energy_over_median": energy[-1] / float(np.median(energy)),
    }


def _candidate_window(rows: Sequence[dict[str, float]], onset: int) -> dict[str, float] | None:
    for start in range(onset, len(rows) - 4):
        for end in range(start + 4, len(rows)):
            span = rows[end]["tau"] - rows[start]["tau"]
            if span < MINIMUM_CANDIDATE_SPAN_TAU:
                continue
            if span > MAXIMUM_CANDIDATE_SPAN_TAU:
                break
            if rows[-1]["tau"] - rows[end]["tau"] < MINIMUM_CONFIRMATION_SPAN_TAU:
                continue
            metrics = _window_metrics(rows, start, end)
            if (
                abs(metrics["log_energy_slope_per_tau"])
                <= MAXIMUM_CANDIDATE_LOG_ENERGY_SLOPE
                and metrics["energy_max_to_min"] <= MAXIMUM_ENERGY_MAX_TO_MIN
            ):
                return metrics
    return None


def _rolling_positive_growth(rows: Sequence[dict[str, float]], start: int) -> list[float]:
    slopes: list[float] = []
    for left in range(start, len(rows) - 1):
        right = left + 1
        while right < len(rows) and rows[right]["tau"] - rows[left]["tau"] < 2.0:
            right += 1
        if right >= len(rows):
            break
        slopes.append(_window_metrics(rows, left, right)["log_energy_slope_per_tau"])
    return slopes


def _confirmation_window(
    rows: Sequence[dict[str, float]], candidate: dict[str, float]
) -> dict[str, Any] | None:
    candidate_start = int(candidate["start_index"])
    for start in range(candidate_start, len(rows) - MINIMUM_CONFIRMATION_SAMPLES + 1):
        metrics = _window_metrics(rows, start, len(rows) - 1)
        if metrics["span_tau"] < MINIMUM_CONFIRMATION_SPAN_TAU:
            continue
        slopes = _rolling_positive_growth(rows, start)
        consecutive_growth = any(
            left > MAXIMUM_CANDIDATE_LOG_ENERGY_SLOPE
            and right > MAXIMUM_CANDIDATE_LOG_ENERGY_SLOPE
            for left, right in zip(slopes, slopes[1:])
        )
        if (
            metrics["sample_count"] >= MINIMUM_CONFIRMATION_SAMPLES
            and abs(metrics["log_energy_slope_per_tau"])
            <= MAXIMUM_CONFIRMATION_LOG_ENERGY_SLOPE
            and metrics["energy_max_to_min"] <= MAXIMUM_ENERGY_MAX_TO_MIN
            and 2.0 / 3.0 <= metrics["terminal_energy_over_median"] <= 1.5
            and not consecutive_growth
        ):
            return {
                **metrics,
                "rolling_two_tau_log_energy_slopes": slopes,
                "maximum_rolling_two_tau_slope": max(slopes, default=0.0),
                "consecutive_positive_growth_windows": consecutive_growth,
            }
    return None


def analyze(run_root: Path, basename: str) -> dict[str, Any]:
    binary_root = run_root / "bin"
    mhd_paths = sorted(binary_root.glob(f"{basename}.mhd_w_bcc.*.bin"), key=_index)
    _require(len(mhd_paths) >= 20, f"too few MHD snapshots: {len(mhd_paths)}")
    products = ("prtcl_rho", "prtcl_jx", "prtcl_jy", "prtcl_jz", "prtcl_dedt")
    by_product = {
        product: {
            _index(path): path
            for path in binary_root.glob(f"{basename}.{product}.*.bin")
        }
        for product in products
    }
    indices = [_index(path) for path in mhd_paths]
    for product, inventory in by_product.items():
        _require(set(inventory) == set(indices), f"{product} output inventory is not matched")

    rows: list[dict[str, float]] = []
    cumulative_cr_to_mhd = 0.0
    prior_time: float | None = None
    prior_transfer_rate: float | None = None
    for path, index in zip(mhd_paths, indices):
        mhd = read_athenak_binary(path)
        _require(
            set(mhd.variable_names) == {
                "dens", "velx", "vely", "velz", "eint", "bcc1", "bcc2", "bcc3"
            },
            f"unexpected MHD inventory in {path}",
        )
        values = {name: _field(mhd, name) for name in mhd.variable_names}
        moment_values: dict[str, np.ndarray] = {}
        for product in products:
            dataset, moment = _load_product(by_product[product][index], product)
            _require(
                dataset.time == mhd.time and dataset.cycle == mhd.cycle,
                f"unmatched time/cycle for {product} output {index}",
            )
            moment_values[product] = moment

        density = values["dens"]
        bperp_squared = values["bcc2"] ** 2 + values["bcc3"] ** 2
        btotal_squared = values["bcc1"] ** 2 + bperp_squared
        jx = moment_values["prtcl_jx"]
        jy = moment_values["prtcl_jy"]
        jz = moment_values["prtcl_jz"]
        volume = (
            (mhd.domain_bounds[1] - mhd.domain_bounds[0])
            * (mhd.domain_bounds[3] - mhd.domain_bounds[2])
            * (mhd.domain_bounds[5] - mhd.domain_bounds[4])
        )
        cr_energy_gain_rate = float(np.mean(moment_values["prtcl_dedt"]) * volume)
        transfer_rate = -cr_energy_gain_rate
        if prior_time is not None and prior_transfer_rate is not None:
            cumulative_cr_to_mhd += 0.5 * (prior_transfer_rate + transfer_rate) * (
                float(mhd.time) - prior_time
            )
        prior_time = float(mhd.time)
        prior_transfer_rate = transfer_rate
        mean_density = float(np.mean(density))
        rows.append(
            {
                "output_index": index,
                "time": float(mhd.time),
                "tau": K0 * UA0 * float(mhd.time),
                "cycle": int(mhd.cycle),
                "bperp_rms_over_B0": float(np.sqrt(np.mean(bperp_squared)) / B0),
                "btotal_rms_over_B0": float(np.sqrt(np.mean(btotal_squared)) / B0),
                "density_rms_over_mean": float(np.std(density) / mean_density),
                "cavity_volume_fraction": float(np.mean(density < 0.5 * mean_density)),
                "filament_volume_fraction": float(
                    np.mean(bperp_squared > 4.0 * np.mean(bperp_squared))
                ),
                "mean_jx_over_initial": float(np.mean(jx) / INITIAL_JX),
                "jperp_rms_over_initial": float(
                    np.sqrt(np.mean(jy * jy + jz * jz)) / INITIAL_JX
                ),
                "cr_to_mhd_energy_transfer_rate": transfer_rate,
                "cumulative_cr_to_mhd_energy_transfer": cumulative_cr_to_mhd,
                "mhd_magnetic_energy": float(0.5 * np.mean(btotal_squared) * volume),
                "mhd_transverse_magnetic_energy": float(
                    0.5 * np.mean(bperp_squared) * volume
                ),
                "mhd_kinetic_energy": float(
                    0.5
                    * np.mean(
                        density
                        * (values["velx"] ** 2 + values["vely"] ** 2 + values["velz"] ** 2)
                    )
                    * volume
                ),
                "mhd_internal_energy": float(np.mean(values["eint"]) * volume),
            }
        )

    onset = next(
        (index for index, row in enumerate(rows) if row["bperp_rms_over_B0"] >= 1.0),
        None,
    )
    candidate = _candidate_window(rows, onset) if onset is not None else None
    confirmation = _confirmation_window(rows, candidate) if candidate is not None else None
    current_relaxation = 1.0 - rows[-1]["mean_jx_over_initial"]
    morphology_ready = bool(
        onset is not None
        and max(row["density_rms_over_mean"] for row in rows[onset:])
        >= MINIMUM_DENSITY_RMS_OVER_MEAN
        and max(row["filament_volume_fraction"] for row in rows[onset:])
        >= MINIMUM_FILAMENT_VOLUME_FRACTION
    )
    transfer_ready = rows[-1]["cumulative_cr_to_mhd_energy_transfer"] > 0.0
    current_relaxation_ready = current_relaxation >= MINIMUM_CURRENT_RELAXATION_FRACTION
    saturation_ready = bool(
        onset is not None
        and candidate is not None
        and confirmation is not None
        and morphology_ready
        and transfer_ready
        and current_relaxation_ready
    )
    return {
        "record_type": "q019_bell_saturation_engineering_analysis_v1",
        "run_root": str(run_root),
        "basename": basename,
        "snapshot_count": len(rows),
        "nonlinear_onset": (
            {
                "index": onset,
                "time": rows[onset]["time"],
                "tau": rows[onset]["tau"],
                "bperp_rms_over_B0": rows[onset]["bperp_rms_over_B0"],
            }
            if onset is not None
            else None
        ),
        "candidate_plateau": candidate,
        "sustained_confirmation": confirmation,
        "peak_bperp_rms_over_B0": max(row["bperp_rms_over_B0"] for row in rows),
        "final_bperp_rms_over_B0": rows[-1]["bperp_rms_over_B0"],
        "final_current_relaxation_fraction": current_relaxation,
        "morphology_ready": morphology_ready,
        "positive_cr_to_mhd_energy_transfer": transfer_ready,
        "current_relaxation_ready": current_relaxation_ready,
        "saturation_ready": saturation_ready,
        "criteria": {
            "onset_bperp_rms_over_B0": ONSET_BPERP_OVER_B0,
            "candidate_span_tau": [MINIMUM_CANDIDATE_SPAN_TAU, MAXIMUM_CANDIDATE_SPAN_TAU],
            "maximum_candidate_absolute_log_energy_slope_per_tau": (
                MAXIMUM_CANDIDATE_LOG_ENERGY_SLOPE
            ),
            "minimum_confirmation_span_tau": MINIMUM_CONFIRMATION_SPAN_TAU,
            "minimum_confirmation_samples": MINIMUM_CONFIRMATION_SAMPLES,
            "maximum_confirmation_absolute_log_energy_slope_per_tau": (
                MAXIMUM_CONFIRMATION_LOG_ENERGY_SLOPE
            ),
            "maximum_energy_max_to_min": MAXIMUM_ENERGY_MAX_TO_MIN,
            "minimum_current_relaxation_fraction": MINIMUM_CURRENT_RELAXATION_FRACTION,
            "minimum_density_rms_over_mean": MINIMUM_DENSITY_RMS_OVER_MEAN,
            "minimum_filament_volume_fraction": MINIMUM_FILAMENT_VOLUME_FRACTION,
        },
        "snapshots": rows,
        "claim_scope": (
            "compact_single_run_3d_finite_rigidity_engineering_evidence;"
            "not_convergence_box_independence_or_publication_qualification"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("basename")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-saturation", action="store_true")
    args = parser.parse_args()
    report = analyze(args.run_root.resolve(strict=True), args.basename)
    payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload, encoding="utf-8")
    return 0 if (not args.require_saturation or report["saturation_ready"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
