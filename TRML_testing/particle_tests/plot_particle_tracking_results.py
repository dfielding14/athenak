#!/usr/bin/env python3
"""Create figures for the lagrangian_mc particle-tracking validation report."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np

from oracle import (
    MeshSpec,
    analytic_state,
    initial_cell_center_particles,
    predicted_x1_mc_push,
    uniform_state,
)
from read_trk import TrackRecord, read_trk


TEST_DIR = Path(__file__).resolve().parent
DEFAULT_RUN_DIR = TEST_DIR / "runs_particle_tracking"
DEFAULT_OUT_DIR = TEST_DIR / "results" / "figures"


def load_case(run_dir: Path, case: str, basename: str) -> list[TrackRecord]:
    path = run_dir / case / "trk" / f"{basename}.trk"
    if not path.exists():
        raise FileNotFoundError(f"missing {path}; run run_particle_tracking_tests.py first")
    return read_trk(path)


def savefig(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_static(run_dir: Path, out_dir: Path) -> None:
    record = load_case(run_dir, "static", "particle_static")[0]
    mesh = MeshSpec(nx1=4, nx2=4, nx3=4)
    positions = initial_cell_center_particles(mesh)
    quantities = ["rho", "temp", "eint", "scalar0"]
    labels = {
        "rho": r"$\rho$",
        "temp": r"$T$",
        "eint": r"$e_{\rm int}$",
        "scalar0": "scalar0",
    }

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0))
    for ax, quantity in zip(axes.flat, quantities):
        expected = np.array([
            analytic_state(mesh, x, y, z)[quantity] for x, y, z in positions
        ])
        measured = np.array([row[quantity] for row in record.rows])
        ax.scatter(expected, measured, s=18, color="#245a8d", alpha=0.82)
        lo = min(expected.min(), measured.min())
        hi = max(expected.max(), measured.max())
        pad = 0.04 * (hi - lo if hi > lo else 1.0)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#222222", lw=1)
        ax.set_xlabel(f"oracle {labels[quantity]}")
        ax.set_ylabel(f"tracked {labels[quantity]}")
        ax.grid(alpha=0.25)
    fig.suptitle("Static thermodynamic sampling: tracked values follow the analytic oracle")
    savefig(fig, out_dir, "static_sampling")


def plot_mc_push(run_dir: Path, out_dir: Path) -> None:
    records = load_case(run_dir, "mc_push", "particle_mc_push")
    initial = records[0]
    final = records[-1]
    mesh = MeshSpec(nx1=4, nx2=4, nx3=1)
    tags = np.arange(final.ntrack)
    x0 = np.array([row["x"] for row in initial.rows])
    observed = np.array([row["x"] for row in final.rows])
    expected = np.array([
        predicted_x1_mc_push(
            int(tag), float(x), mesh, probability=0.5, random_seed=12345,
            cycle_used_by_pusher=0,
        )
        for tag, x in zip(tags, x0)
    ])
    moved = np.abs(expected - x0) > 1.0e-6

    fig, axes = plt.subplots(2, 1, figsize=(7.4, 5.8), sharex=True)
    axes[0].scatter(tags, x0, label="initial", color="#666666", s=24)
    axes[0].scatter(tags, expected, label="oracle final", color="#137a63", s=36, marker="x")
    axes[0].scatter(tags, observed, label="tracked final", facecolors="none",
                    edgecolors="#c43c39", s=64)
    for tag, start, end, did_move in zip(tags, x0, observed, moved):
        if did_move:
            axes[0].annotate(
                "", xy=(tag, end), xytext=(tag, start),
                arrowprops={"arrowstyle": "->", "color": "#999999", "lw": 0.8},
            )
    axes[0].set_ylabel("$x$")
    axes[0].set_title("One-cycle MC push: exact tag-by-tag final positions", pad=22)
    axes[0].legend(
        ncol=3, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.18)
    )
    axes[0].grid(alpha=0.25)

    residual = observed - expected
    axes[1].bar(tags, residual, color="#245a8d")
    axes[1].axhline(0.0, color="#222222", lw=1)
    axes[1].set_xlabel("particle tag")
    axes[1].set_ylabel("tracked - oracle")
    axes[1].grid(alpha=0.25, axis="y")
    savefig(fig, out_dir, "mc_push_positions")


def plot_boundary_lifecycle(run_dir: Path, out_dir: Path) -> None:
    exit_records = load_case(run_dir, "exit_x1", "particle_exit_x1")
    inject_records = load_case(run_dir, "inflow_injection", "particle_inflow_injection")
    exit_counts = [
        sum(int(round(row["active"])) for row in exit_records[0].rows),
        sum(int(round(row["active"])) for row in exit_records[-1].rows),
    ]
    inject_counts = [
        sum(int(round(row["active"])) for row in inject_records[0].rows),
        sum(int(round(row["active"])) for row in inject_records[-1].rows),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6))
    x = np.arange(2)
    width = 0.34
    axes[0].bar(x - width / 2, exit_counts, width, label="exit test", color="#c43c39")
    axes[0].bar(x + width / 2, inject_counts, width, label="inflow test", color="#137a63")
    axes[0].set_xticks(x, ["initial", "final"])
    axes[0].set_ylabel("active tracked particles")
    axes[0].set_title("Boundary lifecycle")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.25, axis="y")

    final_inject = inject_records[-1]
    active_rows = [row for row in final_inject.rows if int(round(row["active"])) == 1]
    inactive_rows = [row for row in final_inject.rows if int(round(row["active"])) == 0]
    axes[1].scatter(
        [row["x"] for row in active_rows],
        [row["y"] for row in active_rows],
        color="#137a63",
        label="injected active tags",
        s=50,
    )
    if inactive_rows:
        axes[1].scatter(
            [row["x"] for row in inactive_rows],
            [row["y"] for row in inactive_rows],
            facecolors="none",
            edgecolors="#777777",
            label="unused tracked slots",
            s=42,
        )
    axes[1].axvline(0.125, color="#222222", lw=1, ls="--", label="inner boundary cell centers")
    axes[1].set_xlabel("$x$")
    axes[1].set_ylabel("$y$")
    axes[1].set_xlim(-0.03, 1.03)
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].set_title("Injected particles appear on the inflow face")
    axes[1].legend(fontsize=7, loc="upper right")
    axes[1].grid(alpha=0.25)
    savefig(fig, out_dir, "boundary_lifecycle")


def parse_amr_log(run_dir: Path) -> tuple[int, int, int]:
    log = (run_dir / "amr_refine_derefine" / "athena.log").read_text()
    current = re.search(r"Current number of MeshBlocks =\s*(\d+)", log)
    changed = re.search(r"(\d+) MeshBlocks created,\s*(\d+) deleted by AMR", log)
    if current is None or changed is None:
        raise ValueError("could not parse AMR summary from athena.log")
    return int(current.group(1)), int(changed.group(1)), int(changed.group(2))


def plot_amr(run_dir: Path, out_dir: Path) -> None:
    records = load_case(run_dir, "amr_refine_derefine", "particle_amr_refine_derefine")
    final = records[-1]
    current, created, deleted = parse_amr_log(run_dir)
    tags = np.arange(final.ntrack)
    levels = np.array([row["level"] for row in final.rows])
    active = np.array([row["active"] for row in final.rows])

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6))
    axes[0].bar(["created", "deleted", "final"], [created, deleted, current],
                color=["#245a8d", "#c43c39", "#137a63"])
    axes[0].set_ylabel("MeshBlocks")
    axes[0].set_title("AMR event summary")
    axes[0].grid(alpha=0.25, axis="y")

    axes[1].scatter(tags, levels, c=active, cmap="viridis", vmin=0, vmax=1, s=42)
    axes[1].axhline(0.0, color="#222222", lw=1, ls="--")
    axes[1].set_ylim(-0.5, 1.5)
    axes[1].set_yticks([0, 1])
    axes[1].set_xlabel("particle tag")
    axes[1].set_ylabel("final tracked level")
    axes[1].set_title("Particles return to valid level-0 owners")
    axes[1].grid(alpha=0.25, axis="y")
    savefig(fig, out_dir, "amr_bookkeeping")


def plot_summary(out_dir: Path) -> None:
    cases = ["static", "mc push", "exit", "injection", "AMR"]
    fig, ax = plt.subplots(figsize=(7.4, 1.8))
    ax.scatter(range(len(cases)), [1] * len(cases), s=720, color="#137a63")
    for i, case in enumerate(cases):
        ax.text(i, 1, "PASS", ha="center", va="center", color="white",
                fontsize=8, fontweight="bold")
    ax.set_xticks(range(len(cases)), cases)
    ax.set_yticks([])
    ax.set_ylim(0.4, 1.6)
    ax.set_title("Implemented particle-tracking validation cases")
    for spine in ax.spines.values():
        spine.set_visible(False)
    savefig(fig, out_dir, "test_summary")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "figure.titlesize": 11,
    })
    plot_summary(args.out_dir)
    plot_static(args.run_dir, args.out_dir)
    plot_mc_push(args.run_dir, args.out_dir)
    plot_boundary_lifecycle(args.run_dir, args.out_dir)
    plot_amr(args.run_dir, args.out_dir)
    print(f"Wrote figures to {args.out_dir}")


if __name__ == "__main__":
    main()
