#!/usr/bin/env python3
"""Generate figures for docs/cgl_lf_validation.tex from quantitative CGL-LF runs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_metrics(path: Path) -> dict[str, float]:
    metrics: dict[str, float] = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            metrics[row["metric"]] = float(row["value"])
    return metrics


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def column(rows: list[dict[str, str]], name: str) -> np.ndarray:
    return np.array([float(row[name]) for row in rows], dtype=float)


def complex_column(rows: list[dict[str, str]], prefix: str) -> np.ndarray:
    real = column(rows, f"{prefix}_real")
    imag = column(rows, f"{prefix}_imag")
    return real + 1j*imag


def save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_decay(data_dir: Path, figure_dir: Path) -> None:
    cases = [
        ("parallel", data_dir / "cgl_lf_quant_parallel.decay.csv",
         r"$T_\parallel$ decay", r"$\sqrt{8/\pi}\,c_\parallel/|k_\parallel|$"),
        ("perp", data_dir / "cgl_lf_quant_perp.decay.csv",
         r"$T_\perp$ decay", r"$\sqrt{2/\pi}\,c_\parallel/|k_\parallel|$"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.0), sharey=False)
    for ax, (_, path, title, coeff_label) in zip(axes, cases):
        metrics = read_metrics(path)
        chi_key = "chi_parallel" if "parallel" in path.name else "chi_perp"
        chi = metrics.get("chi", metrics[chi_key])
        t_final = metrics["time"]
        time = np.linspace(0.0, t_final, 256)
        analytic = metrics["initial_amp"]*np.exp(
            -chi*metrics["k_wave"]**2*time)
        ax.plot(time, analytic, color="#1f77b4", lw=2.0, label="analytic")
        ax.scatter([t_final], [metrics["measured_amp"]], color="#d62728",
                   zorder=3, label="AthenaK")
        ax.set_title(title)
        ax.set_xlabel("time")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        ax.text(0.04, 0.08,
                f"{coeff_label}\nrel. error = {metrics['rel_err']:.2e}",
                transform=ax.transAxes, fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "0.8", "pad": 3})
        ax.grid(True, alpha=0.25)
    axes[0].set_ylabel("Fourier amplitude")
    axes[0].legend(loc="upper right", fontsize=8)
    fig.subplots_adjust(wspace=0.32)
    save(fig, figure_dir / "cgl_lf_decay.pdf")


def plot_grad_b(data_dir: Path, figure_dir: Path) -> None:
    rows = read_rows(data_dir / "cgl_lf_quant_grad_b.grad_b.csv")
    x = column(rows, "x")
    measured = column(rows, "measured_delta")
    expected = column(rows, "expected_delta")
    rel = np.linalg.norm(measured - expected)/np.linalg.norm(expected)

    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.plot(x, expected, color="#1f77b4", lw=2.0, label="offline FV RHS")
    ax.scatter(x, measured, s=14, color="#d62728", label="AthenaK")
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\Delta(p_\perp/|B|)$")
    ax.set_title(r"Isolated $\nabla_\parallel |B|$ response")
    ax.text(0.04, 0.08, f"RMS relative error = {rel:.2e}",
            transform=ax.transAxes, fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "0.8", "pad": 3})
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=8)
    save(fig, figure_dir / "cgl_lf_grad_b.pdf")


def plot_flux_limiter(data_dir: Path, figure_dir: Path) -> None:
    face_rows = read_rows(data_dir / "cgl_lf_flux_limiter.flux_limiter_faces.csv")
    cell_rows = read_rows(data_dir / "cgl_lf_flux_limiter.flux_limiter_cells.csv")

    x_face = column(face_rows, "x_face")
    q_unlimited = column(face_rows, "q_unlimited")
    q_limited = column(face_rows, "q_limited")
    q_max = column(face_rows, "q_max")
    x = column(cell_rows, "x")
    measured_delta = column(cell_rows, "measured_ppar") - column(cell_rows, "initial_ppar")
    expected_delta = column(cell_rows, "expected_ppar") - column(cell_rows, "initial_ppar")
    rel = np.linalg.norm(measured_delta - expected_delta)/np.linalg.norm(expected_delta)
    max_ratio = np.max(np.abs(q_unlimited)/np.maximum(q_max, 1.0e-300))

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1))
    ax = axes[0]
    ax.plot(x_face, q_unlimited, color="0.45", lw=1.2, label=r"unlimited $q_L$")
    ax.plot(x_face, q_limited, color="#d62728", lw=2.0, label="limited q")
    ax.plot(x_face, q_max, color="#1f77b4", lw=1.2, ls="--", label=r"$q_{\max}$")
    ax.plot(x_face, -q_max, color="#1f77b4", lw=1.2, ls="--")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlabel("face x")
    ax.set_ylabel(r"$q_\parallel$")
    ax.set_title("Flux cap")
    ax.text(0.04, 0.08, rf"max $|q_L|/q_{{max}}$ = {max_ratio:.1f}",
            transform=ax.transAxes, fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "0.8", "pad": 3})
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=7)

    ax = axes[1]
    ax.plot(x, expected_delta, color="#1f77b4", lw=2.0, label="offline FV update")
    ax.scatter(x, measured_delta, s=14, color="#d62728", label="AthenaK")
    ax.set_xlabel("cell x")
    ax.set_ylabel(r"$\Delta p_\parallel$")
    ax.set_title("Limited update")
    ax.text(0.04, 0.08, f"RMS relative error = {rel:.2e}",
            transform=ax.transAxes, fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "0.8", "pad": 3})
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=7)
    save(fig, figure_dir / "cgl_lf_flux_limiter.pdf")


def plot_paper_oblique(data_dir: Path, figure_dir: Path) -> None:
    cases = [
        ("Pure CGL", data_dir / "cgl_pure_paper_oblique_wave.paper_oblique_wave.csv"),
        ("CGL-LF", data_dir / "cgl_lf_paper_oblique_wave.paper_oblique_wave.csv"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.2), sharex=True)

    for row_idx, (title, path) in enumerate(cases):
        rows = read_rows(path)
        labels = [row["variable"] for row in rows]
        x = np.arange(len(labels), dtype=float)
        measured = complex_column(rows, "measured")
        reference = complex_column(rows, "reference")
        rel = column(rows, "rel_err")

        ax = axes[row_idx, 0]
        width = 0.34
        ax.bar(x - width/2, np.abs(reference), width, color="#1f77b4", label="reference")
        ax.bar(x + width/2, np.abs(measured), width, color="#d62728", label="AthenaK")
        ax.set_ylabel("mode amplitude")
        ax.set_title(title)
        ax.grid(True, axis="y", alpha=0.25)
        if row_idx == 0:
            ax.legend(loc="upper right", fontsize=7)

        ax = axes[row_idx, 1]
        ax.scatter(x, np.angle(reference), color="#1f77b4", marker="o", label="reference")
        ax.scatter(x, np.angle(measured), color="#d62728", marker="x", label="AthenaK")
        for xi, err in zip(x, rel):
            ax.text(xi, 0.02, f"{err:.1e}", rotation=90, va="bottom", ha="center",
                    fontsize=7, transform=ax.get_xaxis_transform())
        ax.set_ylabel("phase [rad]")
        ax.grid(True, alpha=0.25)
        if row_idx == 0:
            ax.legend(loc="upper right", fontsize=7)

        for ax in axes[row_idx, :]:
            ax.set_xticks(x)
            ax.set_xticklabels(labels)

    axes[1, 0].set_xlabel("Fourier component")
    axes[1, 1].set_xlabel("Fourier component")
    axes[0, 1].set_title("complex phase")
    save(fig, figure_dir / "cgl_lf_paper_oblique_wave.pdf")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="docs/figures/data", type=Path)
    parser.add_argument("--figure-dir", default="docs/figures", type=Path)
    args = parser.parse_args()

    args.figure_dir.mkdir(parents=True, exist_ok=True)
    plot_decay(args.data_dir, args.figure_dir)
    plot_grad_b(args.data_dir, args.figure_dir)
    plot_flux_limiter(args.data_dir, args.figure_dir)
    plot_paper_oblique(args.data_dir, args.figure_dir)
    print(f"Wrote CGL-LF validation figures to {args.figure_dir}")


if __name__ == "__main__":
    main()
