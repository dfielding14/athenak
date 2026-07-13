#!/usr/bin/env python3
"""Compare matched 128^3 full-CR-Hall and Hall-off turbulence reports."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


EVIDENCE_NOTE = "engineering_proxy | not_sun_bai_reproduction=true"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _load(path: Path) -> dict[str, Any]:
    report = json.loads(path.read_text(encoding="utf-8"))
    _require(report.get("schema_version") == 1, f"{path}: unsupported schema")
    _require(report.get("evidence_class") == "engineering_proxy",
             f"{path}: unexpected evidence class")
    _require(report.get("not_sun_bai_reproduction") is True,
             f"{path}: missing non-reproduction boundary")
    _require(report.get("acceptance", {}).get("passed") is True,
             f"{path}: source analysis did not pass")
    for key in ("model", "history", "snapshots"):
        _require(key in report, f"{path}: missing {key}")
    return report


def _array(report: dict[str, Any], name: str) -> np.ndarray:
    values = np.asarray(report["history"][name], dtype=np.float64)
    _require(values.ndim == 1 and values.size > 0 and np.all(np.isfinite(values)),
             f"invalid history array: {name}")
    return values


def _relative(numerator: float, denominator: float) -> float:
    return (numerator - denominator) / max(abs(denominator), 1.0e-30)


def _model_without_mode(report: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in report["model"].items()
            if key != "hall_mode"}


def _binding_hashes(path: Path) -> set[str]:
    hashes = set()
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split()
        if len(fields) >= 2 and len(fields[0]) == 64:
            hashes.add(fields[0])
    _require(hashes, f"{path}: no SHA-256 bindings found")
    return hashes


def _normalized_overrides(path: Path) -> tuple[set[str], str, str]:
    lines = {
        line.strip() for line in path.read_text(encoding="ascii").splitlines()
        if line.strip()
    }
    basenames = [line for line in lines if line.startswith("job/basename=")]
    hall_modes = [
        line for line in lines if line.startswith("particles/pic_cr_hall_mode=")
    ]
    _require(len(basenames) == 1 and len(hall_modes) == 1,
             f"{path}: incomplete runtime overrides")
    normalized = lines - {basenames[0], hall_modes[0]}
    return normalized, basenames[0], hall_modes[0]


def _series_summary(
    full_values: np.ndarray, off_values: np.ndarray,
    full_time: np.ndarray, off_time: np.ndarray, tail_start: float,
) -> dict[str, Any]:
    full_tail = full_values[full_time >= tail_start]
    off_tail = off_values[off_time >= tail_start]
    _require(full_tail.size >= 5 and off_tail.size >= 5,
             "comparison tail has fewer than five samples")
    full_median = float(np.median(full_tail))
    off_median = float(np.median(off_tail))
    return {
        "terminal": {
            "full_hall": float(full_values[-1]),
            "hall_off": float(off_values[-1]),
            "fractional_difference_full_minus_off": _relative(
                float(full_values[-1]), float(off_values[-1])
            ),
        },
        "final_30_percent_median": {
            "full_hall": full_median,
            "hall_off": off_median,
            "fractional_difference_full_minus_off": _relative(
                full_median, off_median
            ),
        },
    }


def _style() -> None:
    plt.rcParams.update({
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "axes.linewidth": 0.9,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.frameon": False,
        "lines.linewidth": 1.8,
    })


def _plot(full: dict[str, Any], off: dict[str, Any], path: Path) -> None:
    full_color = "#d1495b"
    off_color = "#3b4cc0"
    user_time = _array(full, "user_time")
    off_user_time = _array(off, "user_time")
    full_b = _array(full, "magnetic_rms")
    off_b = _array(off, "magnetic_rms")
    full_v = _array(full, "velocity_rms")
    off_v = _array(off, "velocity_rms")

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.4), constrained_layout=True)
    axes[0, 0].semilogy(user_time, full_b, color=full_color, label="full Hall")
    axes[0, 0].semilogy(off_user_time, off_b, color=off_color, linestyle="--",
                       label="Hall off")
    axes[0, 0].set(xlabel="time", ylabel=r"$B_{\rm rms}$",
                   title="Magnetic amplification")
    axes[0, 0].legend()

    axes[0, 1].plot(user_time, full_v, color=full_color, label="full Hall")
    axes[0, 1].plot(off_user_time, off_v, color=off_color, linestyle="--",
                    label="Hall off")
    axes[0, 1].set(xlabel="time", ylabel=r"$v_{\rm rms}$",
                   title="Turbulent velocity")
    axes[0, 1].legend()

    full_snapshot = full["snapshots"][-1]
    off_snapshot = off["snapshots"][-1]
    full_mode = np.asarray(full_snapshot["spectrum_mode"], dtype=np.float64)
    off_mode = np.asarray(off_snapshot["spectrum_mode"], dtype=np.float64)
    full_spectrum = np.asarray(
        full_snapshot["spectrum_magnetic_energy"], dtype=np.float64)
    off_spectrum = np.asarray(
        off_snapshot["spectrum_magnetic_energy"], dtype=np.float64)
    axes[1, 0].loglog(full_mode, full_spectrum, color=full_color,
                      label="full Hall")
    axes[1, 0].loglog(off_mode, off_spectrum, color=off_color,
                      linestyle="--", label="Hall off")
    axes[1, 0].set(xlabel=r"isotropic mode $kL/(2\pi)$",
                   ylabel=r"shell magnetic energy $E_B(k)$",
                   title="Final magnetic spectrum")
    axes[1, 0].legend()

    off_b_match = np.interp(user_time, off_user_time, off_b)
    off_v_match = np.interp(user_time, off_user_time, off_v)
    valid_b = off_b_match > max(float(np.max(off_b_match))*1.0e-12, 1.0e-30)
    valid_v = off_v_match > max(float(np.max(off_v_match))*1.0e-12, 1.0e-30)
    axes[1, 1].plot(user_time[valid_b], full_b[valid_b]/off_b_match[valid_b],
                    color=full_color, label=r"$B_{\rm rms}$ ratio")
    axes[1, 1].plot(user_time[valid_v], full_v[valid_v]/off_v_match[valid_v],
                    color="#2a9d8f", label=r"$v_{\rm rms}$ ratio")
    axes[1, 1].axhline(1.0, color="0.25", linestyle=":")
    axes[1, 1].set(xlabel="time", ylabel="full Hall / Hall off",
                   title="Matched response ratio")
    axes[1, 1].legend()

    for axis in axes.flat:
        axis.grid(alpha=0.18, which="both")
    fig.suptitle("128$^3$ turbulent-dynamo CR-Hall comparison")
    fig.text(0.995, 0.995, EVIDENCE_NOTE, ha="right", va="top",
             fontsize=7, color="0.35")
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-report", type=Path, required=True)
    parser.add_argument("--off-report", type=Path, required=True)
    parser.add_argument("--full-bindings", type=Path, required=True)
    parser.add_argument("--off-bindings", type=Path, required=True)
    parser.add_argument("--full-git-head", type=Path, required=True)
    parser.add_argument("--off-git-head", type=Path, required=True)
    parser.add_argument("--full-runtime-overrides", type=Path, required=True)
    parser.add_argument("--off-runtime-overrides", type=Path, required=True)
    parser.add_argument("--expected-snapshot-sha", required=True)
    parser.add_argument("--expected-input-sha", required=True)
    parser.add_argument("--expected-git-head", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    full = _load(args.full_report)
    off = _load(args.off_report)
    _require(full["model"].get("hall_mode") == "full",
             "full report does not select full Hall")
    _require(off["model"].get("hall_mode") == "off",
             "control report does not select Hall off")
    _require(_model_without_mode(full) == _model_without_mode(off),
             "full-Hall and Hall-off model settings differ")

    full_hashes = _binding_hashes(args.full_bindings)
    off_hashes = _binding_hashes(args.off_bindings)
    full_head = args.full_git_head.read_text(encoding="ascii").strip()
    off_head = args.off_git_head.read_text(encoding="ascii").strip()
    _require(args.expected_snapshot_sha in full_hashes and
             args.expected_snapshot_sha in off_hashes,
             "source-snapshot binding differs")
    _require(args.expected_input_sha in full_hashes and
             args.expected_input_sha in off_hashes, "input binding differs")
    _require(full_head == off_head == args.expected_git_head,
             "recorded git-head binding differs")
    full_overrides, _, full_mode_override = _normalized_overrides(
        args.full_runtime_overrides)
    off_overrides, _, off_mode_override = _normalized_overrides(
        args.off_runtime_overrides)
    _require(full_overrides == off_overrides,
             "runtime overrides differ beyond basename and Hall mode")
    _require(full_mode_override == "particles/pic_cr_hall_mode=full" and
             off_mode_override == "particles/pic_cr_hall_mode=off",
             "runtime Hall-mode overrides are not the matched full/off pair")

    time = _array(full, "time")
    off_time = _array(off, "time")
    user_time = _array(full, "user_time")
    off_user_time = _array(off, "user_time")
    _require(np.all(np.diff(time) > 0.0) and np.all(np.diff(off_time) > 0.0),
             "MHD history times are not monotone")
    _require(np.all(np.diff(user_time) > 0.0) and
             np.all(np.diff(off_user_time) > 0.0),
             "user history times are not monotone")
    _require(math.isclose(float(time[0]), float(off_time[0]), abs_tol=1.0e-12) and
             math.isclose(float(time[-1]), float(off_time[-1]), abs_tol=1.0e-12),
             "MHD history time ranges differ")
    _require(math.isclose(float(user_time[0]), float(off_user_time[0]),
                          abs_tol=1.0e-12) and
             math.isclose(float(user_time[-1]), float(off_user_time[-1]),
                          abs_tol=1.0e-12),
             "user history time ranges differ")

    full_snapshots = full["snapshots"]
    off_snapshots = off["snapshots"]
    _require(len(full_snapshots) == len(off_snapshots) >= 3,
             "snapshot inventories differ")
    full_snapshot_times = np.asarray(
        [item["time"] for item in full_snapshots], dtype=np.float64)
    off_snapshot_times = np.asarray(
        [item["time"] for item in off_snapshots], dtype=np.float64)
    _require(np.allclose(full_snapshot_times, off_snapshot_times,
                         rtol=0.0, atol=1.0e-2), "snapshot times differ")

    history_names = (
        "kinetic_energy", "magnetic_energy", "cr_energy",
        "velocity_rms", "magnetic_rms",
    )
    mhd_names = {"kinetic_energy", "magnetic_energy", "cr_energy"}
    tail_start = 0.7*float(full["model"]["tlim"])
    history_summary = {}
    for name in history_names:
        full_values = _array(full, name)
        off_values = _array(off, name)
        if name in mhd_names:
            full_series_time, off_series_time = time, off_time
        else:
            full_series_time, off_series_time = user_time, off_user_time
        history_summary[name] = _series_summary(
            full_values, off_values, full_series_time, off_series_time, tail_start
        )

    full_final = full_snapshots[-1]
    off_final = off_snapshots[-1]
    snapshot_names = (
        "density_mean", "density_fractional_rms", "velocity_rms", "magnetic_rms"
    )
    final_snapshot = {
        name: {
            "full_hall": float(full_final[name]),
            "hall_off": float(off_final[name]),
            "fractional_difference_full_minus_off": _relative(
                float(full_final[name]), float(off_final[name])
            ),
        }
        for name in snapshot_names
    }

    full_mode = np.asarray(full_final["spectrum_mode"], dtype=np.float64)
    off_mode = np.asarray(off_final["spectrum_mode"], dtype=np.float64)
    full_spectrum = np.asarray(
        full_final["spectrum_magnetic_energy"], dtype=np.float64)
    off_spectrum = np.asarray(
        off_final["spectrum_magnetic_energy"], dtype=np.float64)
    _require(np.all(np.isfinite(full_spectrum)) and
             np.all(np.isfinite(off_spectrum)), "non-finite final spectrum")
    common_mode, full_index, off_index = np.intersect1d(
        full_mode, off_mode, assume_unique=True, return_indices=True)
    _require(common_mode.size >= 10, "too few common final spectral modes")
    spectrum_l1 = float(
        np.sum(np.abs(full_spectrum[full_index] - off_spectrum[off_index])) /
        max(float(np.sum(np.abs(off_spectrum[off_index]))), 1.0e-30)
    )

    comparison = {
        "schema_version": 1,
        "evidence_class": "engineering_proxy",
        "comparison_role": "matched_engineering_control",
        "not_sun_bai_reproduction": True,
        "full_report": str(args.full_report.resolve()),
        "hall_off_report": str(args.off_report.resolve()),
        "provenance_match": {
            "source_snapshot_sha256": args.expected_snapshot_sha,
            "input_sha256": args.expected_input_sha,
            "git_head": args.expected_git_head,
            "matched": True,
        },
        "matched_model": _model_without_mode(full),
        "full_hall_tail": {
            key: full["hall_off_recommendation"][key]
            for key in (
                "tail_fraction_at_or_above_threshold", "tail_lambda_median",
                "tail_lambda_p90", "tail_lambda_peak",
            )
        },
        "final_snapshot": final_snapshot,
        "history_comparison": history_summary,
        "final_magnetic_spectrum": {
            "common_retained_mode_relative_l1_full_minus_off": spectrum_l1,
            "common_retained_mode_count": int(common_mode.size),
            "full_hall_peak_mode": float(full_mode[int(np.argmax(full_spectrum))]),
            "hall_off_peak_mode": float(off_mode[int(np.argmax(off_spectrum))]),
        },
        "acceptance": {
            "full_hall": full["acceptance"],
            "hall_off": off["acceptance"],
        },
        "completeness": {
            "matched_provenance": True,
            "matched_model": True,
            "compatible_history_ranges": True,
            "compatible_snapshot_times": True,
            "both_source_analyses_passed": True,
            "passed": True,
        },
        "interpretation": (
            "This matched control quantifies the finite-duration response to the "
            "CR-Hall closure; it is not a general turbulence-convergence claim."
        ),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    _style()
    _plot(full, off, args.output_dir / "turbulent_dynamo_hall_comparison.png")
    output = args.output_dir / "turbulent_dynamo_hall_comparison.json"
    output.write_text(json.dumps(comparison, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    print(json.dumps({
        "completeness": comparison["completeness"],
        "final_snapshot": final_snapshot,
        "final_magnetic_spectrum": comparison["final_magnetic_spectrum"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
