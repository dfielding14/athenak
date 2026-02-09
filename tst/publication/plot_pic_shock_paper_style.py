#!/usr/bin/env python3
"""Render an Athena++-style MHD-PIC shock figure from AthenaK outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys

import numpy as np
from matplotlib import colors

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError as exc:
    raise SystemExit("matplotlib is required: python3 -m pip install matplotlib") from exc

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
import bin_convert_new as bin_convert  # noqa: E402

_CYCLE_RE = re.compile(r"\.(\d+)\.bin$")


def _latest_bin(bin_dir: Path, basename: str, file_id: str) -> tuple[Path, int]:
    files = sorted(bin_dir.glob(f"{basename}.{file_id}.*.bin"))
    if not files:
        raise RuntimeError(f"No files found for {basename}.{file_id}")

    best = files[-1]
    match = _CYCLE_RE.search(best.name)
    cycle = int(match.group(1)) if match else -1
    return best, cycle


def _meshblock_lines(edges: np.ndarray, block_n: int) -> np.ndarray:
    if block_n <= 0:
        return np.asarray([])

    lines = edges[::block_n]
    if lines.size == 0 or lines[-1] != edges[-1]:
        lines = np.append(lines, edges[-1])
    return lines


def _draw_meshblock_grid(ax, x_lines: np.ndarray, y_lines: np.ndarray) -> None:
    for xv in x_lines:
        ax.axvline(xv, color="k", lw=0.25, alpha=0.55)
    for yv in y_lines:
        ax.axhline(yv, color="k", lw=0.25, alpha=0.55)


def _slice_xy(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 3:
        return np.asarray(arr[arr.shape[0] // 2, :, :], dtype=float)
    if arr.ndim == 2:
        return np.asarray(arr, dtype=float)
    raise RuntimeError(f"Expected 2D/3D field, got ndim={arr.ndim}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a paper-style MHD-PIC shock plot."
    )
    parser.add_argument("--basename", required=True)
    parser.add_argument(
        "--bin-dir",
        default=str(REPO_ROOT / "tst" / "build" / "src" / "bin"),
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--figure-name", default="shock_storyline_paper_style.png")
    parser.add_argument("--mb-nx1", type=int, default=8)
    parser.add_argument("--mb-nx2", type=int, default=8)
    parser.add_argument("--rho-vmax", type=float, default=4.0)
    parser.add_argument("--b-vmin", type=float, default=0.5)
    parser.add_argument("--b-vmax", type=float, default=32.0)
    parser.add_argument("--dpi", type=int, default=220)
    args = parser.parse_args()

    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    })

    bin_dir = Path(args.bin_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rho_file, rho_cycle = _latest_bin(bin_dir, args.basename, "rho")
    b_file, b_cycle = _latest_bin(bin_dir, args.basename, "bmag")
    cycle = min(rho_cycle, b_cycle)

    rho_data = bin_convert.read_binary_as_athdf(str(rho_file))
    b_data = bin_convert.read_binary_as_athdf(str(b_file))

    x1f = np.asarray(rho_data["x1f"], dtype=float)
    x2f = np.asarray(rho_data["x2f"], dtype=float)

    rho_xy = _slice_xy(np.asarray(rho_data["dens"], dtype=float))
    b_xy = _slice_xy(np.asarray(b_data["bmag"], dtype=float))

    # Normalize to upstream-like reference from left-most 10% in x.
    nx = rho_xy.shape[1]
    left_n = max(4, nx // 10)
    rho_ref = float(np.mean(rho_xy[:, :left_n]))
    b_ref = float(np.mean(b_xy[:, :left_n]))
    rho_ref = rho_ref if rho_ref > 0.0 else 1.0
    b_ref = b_ref if b_ref > 0.0 else 1.0

    rho_norm = rho_xy / rho_ref
    b_norm = b_xy / b_ref

    x_lines = _meshblock_lines(x1f, args.mb_nx1)
    y_lines = _meshblock_lines(x2f, args.mb_nx2)

    fig = plt.figure(figsize=(7.0, 8.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.14)

    ax1 = fig.add_subplot(gs[0, 0])
    m1 = ax1.pcolormesh(
        x1f,
        x2f,
        rho_norm,
        shading="auto",
        cmap="jet",
        vmin=0.0,
        vmax=max(args.rho_vmax, float(np.percentile(rho_norm, 99.5))),
    )
    _draw_meshblock_grid(ax1, x_lines, y_lines)
    ax1.set_aspect("equal", adjustable="box")
    ax1.set_xlim(float(x1f[0]), float(x1f[-1]))
    ax1.set_ylim(float(x2f[0]), float(x2f[-1]))
    ax1.set_ylabel(r"$y$", fontsize=14)
    tval = float(rho_data.get("Time", 0.0))
    ax1.set_title(rf"$t={tval:.3e}$", fontsize=18)
    c1 = fig.colorbar(m1, ax=ax1, pad=0.02)
    c1.set_label(r"$\rho/\rho_0$", fontsize=14)

    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    m2 = ax2.pcolormesh(
        x1f,
        x2f,
        np.clip(b_norm, args.b_vmin, None),
        shading="auto",
        cmap="jet",
        norm=colors.LogNorm(vmin=args.b_vmin, vmax=args.b_vmax),
    )
    _draw_meshblock_grid(ax2, x_lines, y_lines)
    ax2.set_aspect("equal", adjustable="box")
    ax2.set_ylabel(r"$y$", fontsize=14)
    ax2.set_xlabel(r"$x$", fontsize=18)
    c2 = fig.colorbar(m2, ax=ax2, pad=0.02)
    c2.set_label(r"$B/B_0$", fontsize=14)

    fig_path = output_dir / args.figure_name
    fig.savefig(fig_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "figure": str(fig_path),
        "rho_file": str(rho_file),
        "b_file": str(b_file),
        "cycle": cycle,
        "time": tval,
        "rho_ref": rho_ref,
        "b_ref": b_ref,
    }
    (output_dir / "paper_style_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("figure:", fig_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
