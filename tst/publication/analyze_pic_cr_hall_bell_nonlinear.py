#!/usr/bin/env python3
"""Reduce the compact full-CR-Hall nonlinear Bell pilot."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tst.publication import bell_saturation_engineering_analysis_v1 as base
from tst.publication.analyze_q011_section54_outputs import (
    compose_leaf_field,
    read_athenak_binary,
)


INITIAL_JX = 4.0 * math.pi
INITIAL_Q = 4.0 * math.pi / 100.0
INITIAL_R = 0.01
INITIAL_LAMBDA = 1.0
BACKGROUND_ALPHA = 12.44070690821558
_COLUMN_RE = re.compile(r"\[([0-9]+)\]=([^\s]+)")


def _history_path(run_root: Path) -> Path:
    paths = sorted(
        set(run_root.glob("*.mhd.hst")) | set((run_root / "bin").glob("*.mhd.hst"))
    )
    if len(paths) != 1:
        raise RuntimeError(f"expected one MHD history file below {run_root}")
    return paths[0]


def _hall_history(run_root: Path) -> dict[str, list[float]]:
    path = _history_path(run_root)
    header = ""
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#") and "[1]=" in line:
                header = line
    labels = {name: int(index) - 1 for index, name in _COLUMN_RE.findall(header)}
    for name in ("time", "hall_Rmax", "hall_Lmax"):
        if name not in labels:
            raise RuntimeError(f"{path} lacks required history column {name}")
    data = np.atleast_2d(np.loadtxt(path))
    return {
        "time": data[:, labels["time"]].tolist(),
        "Rmax": data[:, labels["hall_Rmax"]].tolist(),
        "Lambda_max": data[:, labels["hall_Lmax"]].tolist(),
    }


def _first_deposited_current(rows: list[dict[str, float]]) -> float:
    for row in rows:
        current = float(row["mean_jx_over_initial"])
        if current > 0.5:
            return current * INITIAL_JX
    raise RuntimeError("no deposited CR-current snapshot was found")


def _first_deposited_charge(run_root: Path, basename: str) -> float:
    paths = sorted((run_root / "bin").glob(f"{basename}.prtcl_rho.*.bin"))
    for path in paths:
        dataset = read_athenak_binary(path)
        charge = np.asarray(
            compose_leaf_field(dataset, "prtcl_rho").values, dtype=np.float64
        )
        mean_charge = float(np.mean(charge))
        if abs(mean_charge) > 0.5 * INITIAL_Q:
            return mean_charge
    raise RuntimeError("no deposited CR-charge snapshot was found")


def _requested_terminal_time(run_root: Path, basename: str) -> float:
    paths = sorted((run_root / "bin").glob(f"{basename}.mhd_w_bcc.*.bin"))
    if not paths:
        raise RuntimeError("no MHD snapshot was found")
    parameters = read_athenak_binary(paths[0]).input_parameters
    return float(parameters["time"]["tlim"])


def analyze(run_root: Path, basename: str) -> dict[str, Any]:
    report = base.analyze(run_root, basename)
    hall = _hall_history(run_root)
    time = np.asarray(hall["time"], dtype=np.float64)
    rmax = np.asarray(hall["Rmax"], dtype=np.float64)
    lambda_max = np.asarray(hall["Lambda_max"], dtype=np.float64)
    initialized = (time > 0.0) & np.isfinite(rmax) & np.isfinite(lambda_max)
    if not np.any(initialized):
        raise RuntimeError("Hall diagnostics contain no initialized cycle")
    first = int(np.flatnonzero(initialized)[0])
    measured_jx = _first_deposited_current(report["snapshots"])
    measured_q = _first_deposited_charge(run_root, basename)
    mean_electron_charge = BACKGROUND_ALPHA + measured_q
    measured_mean_r = measured_q / mean_electron_charge
    measured_mean_lambda = abs(measured_jx) / mean_electron_charge
    requested_terminal_time = _requested_terminal_time(run_root, basename)
    measured_terminal_time = float(report["snapshots"][-1]["time"])
    initial = {
        "expected_jcr_over_c": INITIAL_JX,
        "measured_jcr_over_c": measured_jx,
        "current_relative_error": abs(measured_jx - INITIAL_JX) / INITIAL_JX,
        "expected_qcr_over_c": INITIAL_Q,
        "measured_qcr_over_c": measured_q,
        "charge_relative_error": abs(measured_q - INITIAL_Q) / INITIAL_Q,
        "expected_R": INITIAL_R,
        "measured_mean_R": measured_mean_r,
        "first_history_Rmax": float(rmax[first]),
        "expected_Lambda": INITIAL_LAMBDA,
        "measured_mean_Lambda": measured_mean_lambda,
        "first_history_Lambda_max": float(lambda_max[first]),
    }
    checks = {
        "volume_aware_current": initial["current_relative_error"] < 0.01,
        "volume_aware_charge": initial["charge_relative_error"] < 0.01,
        "initial_mean_R": abs(initial["measured_mean_R"] - INITIAL_R) < 0.002,
        "initial_mean_Lambda": (
            abs(initial["measured_mean_Lambda"] - INITIAL_LAMBDA) < 0.05
        ),
        "entered_nonlinear_regime": report["peak_bperp_rms_over_B0"] >= 1.0,
        "reached_terminal_time": (
            measured_terminal_time >= 0.99 * requested_terminal_time
        ),
    }
    report.update(
        {
            "record_type": "pic_cr_hall_bell_nonlinear_pilot_v1",
            "hall_history": hall,
            "initial_normalization": initial,
            "pilot_checks": checks,
            "pilot_passed": all(checks.values()),
            "requested_terminal_time": requested_terminal_time,
            "measured_terminal_time": measured_terminal_time,
            "claim_scope": (
                "compact_3d_full_cr_hall_nonlinear_pilot;"
                "not_saturation_resolution_or_box_size_qualification"
            ),
        }
    )
    return report


def _plot(report: dict[str, Any], output: Path) -> None:
    rows = report["snapshots"]
    tau = np.asarray([row["tau"] for row in rows], dtype=np.float64)
    bperp = np.asarray(
        [row["bperp_rms_over_B0"] for row in rows], dtype=np.float64
    )
    current = np.asarray(
        [row["mean_jx_over_initial"] for row in rows], dtype=np.float64
    )
    if current.size > 1 and current[0] < 0.5 and current[1] > 0.9:
        current[0] = np.nan
    hall = report["hall_history"]

    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 0.9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6), constrained_layout=True)
    magnetic_line = axes[0].semilogy(
        tau, bperp, color="#c7432b", lw=2.0,
        label=r"$B_{\perp,\mathrm{rms}}/B_0$"
    )
    axes[0].axhline(1.0, color="#555555", lw=1.0, ls="--")
    current_axis = axes[0].twinx()
    current_line = current_axis.plot(
        tau, current, color="#2468a2", lw=1.6,
        label=r"$\langle J_{\rm cr,x}\rangle/J_{\rm cr,x}(0)$"
    )
    axes[0].set(xlabel=r"$k_0v_{A0}t$", ylabel=r"$B_{\perp,\mathrm{rms}}/B_0$")
    current_axis.set_ylabel(r"$\langle J_{\rm cr,x}\rangle/J_{\rm cr,x}(0)$")
    axes[0].legend(magnetic_line + current_line,
                   [line.get_label() for line in magnetic_line + current_line],
                   frameon=False, loc="best")
    axes[0].text(0.03, 0.94, "(a)", transform=axes[0].transAxes,
                 ha="left", va="top", fontweight="bold")

    axes[1].plot(hall["time"], hall["Rmax"], color="#4c78a8", lw=1.8,
                 label=r"$\max |R|$")
    axes[1].plot(hall["time"], hall["Lambda_max"], color="#d95f02", lw=1.8,
                 label=r"$\max \Lambda$")
    axes[1].set(xlabel="time", ylabel="CR-Hall strength")
    axes[1].legend(frameon=False)
    axes[1].text(0.03, 0.94, "(b)", transform=axes[1].transAxes,
                 ha="left", va="top", fontweight="bold")
    for axis in axes:
        axis.grid(alpha=0.18, lw=0.6)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("basename")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--require-nonlinear", action="store_true")
    args = parser.parse_args()
    report = analyze(args.run_root.resolve(strict=True), args.basename)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / "pic_cr_hall_bell_nonlinear_pilot.json"
    figure_path = args.output_dir / "pic_cr_hall_bell_nonlinear_pilot.png"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    _plot(report, figure_path)
    print(json.dumps({"checks": report["pilot_checks"],
                      "pilot_passed": report["pilot_passed"]}, indent=2))
    if args.require_nonlinear and not report["pilot_passed"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
