#!/usr/bin/env python3
"""Render an unqualified Section-5.4-inspired engineering quick-look figure."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
import bin_convert_new as bin_convert  # noqa: E402

from pvtk_particles import read_particle_vtk  # noqa: E402
from artifact_lineage import write_companion_manifest  # noqa: E402

_BIN_CYCLE_RE = re.compile(r"\.(\d+)\.bin$")
_PVTK_CYCLE_RE = re.compile(r"\.prtcl_all\.(\d+)\.part\.vtk$")
_PROXY_WATERMARK = "ENGINEERING PROXY - NOT SUN & BAI (2023) REPRODUCTION"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fobj:
        for chunk in iter(lambda: fobj.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _files_by_cycle(bin_dir: Path, basename: str, file_id: str) -> dict[int, Path]:
    out: dict[int, Path] = {}
    for path in sorted(bin_dir.glob(f"{basename}.{file_id}.*.bin")):
        match = _BIN_CYCLE_RE.search(path.name)
        if match is None:
            continue
        out[int(match.group(1))] = path
    return out


def _pvtk_by_cycle(pvtk_dir: Path, basename: str) -> dict[int, Path]:
    out: dict[int, Path] = {}
    for path in sorted(pvtk_dir.glob(f"{basename}.prtcl_all.*.part.vtk")):
        match = _PVTK_CYCLE_RE.search(path.name)
        if match is None:
            continue
        out[int(match.group(1))] = path
    return out


def _select_pvtk_for_cycle(
    pvtk_map: dict[int, Path], target_cycle: int
) -> tuple[Path | None, int | None]:
    if not pvtk_map:
        return None, None
    if target_cycle in pvtk_map:
        return pvtk_map[target_cycle], target_cycle

    lower = [cyc for cyc in pvtk_map if cyc <= target_cycle]
    if lower:
        cyc = max(lower)
        return pvtk_map[cyc], cyc

    cyc = min(pvtk_map)
    return pvtk_map[cyc], cyc


def _slice_xy(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 3:
        return np.asarray(arr[arr.shape[0] // 2, :, :], dtype=float)
    if arr.ndim == 2:
        return np.asarray(arr, dtype=float)
    raise RuntimeError(f"Expected 2D/3D field, got ndim={arr.ndim}")


def _meshblock_lines(edges: np.ndarray, block_n: int) -> np.ndarray:
    if block_n <= 0:
        return np.asarray([])
    lines = edges[::block_n]
    if lines.size == 0 or lines[-1] != edges[-1]:
        lines = np.append(lines, edges[-1])
    return lines


def _draw_meshblock_grid(ax, x_lines: np.ndarray, y_lines: np.ndarray) -> None:
    for xv in x_lines:
        ax.axvline(xv, color="k", lw=0.25, alpha=0.45)
    for yv in y_lines:
        ax.axhline(yv, color="k", lw=0.25, alpha=0.45)


def _shock_location(rho_xy: np.ndarray, x1v: np.ndarray) -> float:
    rho_x = np.mean(rho_xy, axis=0)
    grad = np.gradient(rho_x, x1v)
    ishock = int(np.argmax(np.abs(grad)))
    return float(x1v[ishock])


def _spectrum_proxy(
    points: np.ndarray,
    vel: np.ndarray,
    u0: float,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
) -> np.ndarray:
    nx = x_edges.size - 1
    ny = y_edges.size - 1
    if points.size == 0 or vel.size == 0:
        return np.zeros((nx, ny), dtype=float)

    x = np.asarray(points[:, 0], dtype=float)
    v2 = np.sum(np.asarray(vel, dtype=float) ** 2, axis=1)
    e_kin = 0.5 * np.maximum(v2, 0.0)
    y = np.log10(np.maximum(v2 / max(u0 * u0, 1.0e-30), 1.0e-30))

    hist_w, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges], weights=e_kin)
    dx = np.diff(x_edges)[:, None]
    dy = np.diff(y_edges)[None, :]
    cell = dx * dy
    proxy = hist_w / max(np.sum(hist_w), 1.0e-30)
    proxy = proxy / np.maximum(cell, 1.0e-30)
    proxy[~np.isfinite(proxy)] = 0.0
    return proxy


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a Section-5.4-inspired engineering quick-look figure."
    )
    parser.add_argument("--basename", required=True)
    parser.add_argument(
        "--bin-dir", default=str(REPO_ROOT / "tst" / "build" / "src" / "bin")
    )
    parser.add_argument(
        "--pvtk-dir", default=str(REPO_ROOT / "tst" / "build" / "src" / "pvtk")
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--figure-name", default="pic_parallel_shock_engineering_quicklook.png"
    )
    parser.add_argument("--mb-nx1", type=int, default=16)
    parser.add_argument("--mb-nx2", type=int, default=16)
    parser.add_argument("--rho-vmax", type=float, default=4.0)
    parser.add_argument("--b-vmin", type=float, default=0.5)
    parser.add_argument("--b-vmax", type=float, default=32.0)
    parser.add_argument("--spec-vmin", type=float, default=1.0e-6)
    parser.add_argument("--spec-vmax", type=float, default=1.0e-2)
    parser.add_argument("--u0", type=float, default=3.0)
    parser.add_argument("--spec-xbins", type=int, default=320)
    parser.add_argument("--spec-ybins", type=int, default=220)
    parser.add_argument("--spec-ymin", type=float, default=-0.5)
    parser.add_argument("--spec-ymax", type=float, default=3.0)
    parser.add_argument("--dpi", type=int, default=230)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
        }
    )

    bin_dir = Path(args.bin_dir).resolve()
    pvtk_dir = Path(args.pvtk_dir).resolve()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    rho_map = _files_by_cycle(bin_dir, args.basename, "rho")
    bmag_map = _files_by_cycle(bin_dir, args.basename, "bmag")
    shared = sorted(set(rho_map) & set(bmag_map))
    if not shared:
        raise RuntimeError("No shared rho/bmag cycles found for basename")
    cycle = shared[-1]
    rho_file = rho_map[cycle]
    bmag_file = bmag_map[cycle]

    rho_data = bin_convert.read_binary_as_athdf(str(rho_file))
    bmag_data = bin_convert.read_binary_as_athdf(str(bmag_file))

    x1f = np.asarray(rho_data["x1f"], dtype=float)
    x2f = np.asarray(rho_data["x2f"], dtype=float)
    x1v = np.asarray(rho_data["x1v"], dtype=float)
    rho_xy = _slice_xy(np.asarray(rho_data["dens"], dtype=float))
    b_xy = _slice_xy(np.asarray(bmag_data["bmag"], dtype=float))

    right_n = max(4, rho_xy.shape[1] // 10)
    rho_ref = float(np.mean(rho_xy[:, -right_n:]))
    b_ref = float(np.mean(b_xy[:, -right_n:]))
    rho_ref = rho_ref if rho_ref > 0.0 else 1.0
    b_ref = b_ref if b_ref > 0.0 else 1.0
    rho_norm = rho_xy / rho_ref
    b_norm = np.clip(b_xy / b_ref, args.b_vmin, None)
    shock_x = _shock_location(rho_xy, x1v)

    pvtk_map = _pvtk_by_cycle(pvtk_dir, args.basename)
    pvtk_file, pvtk_cycle = _select_pvtk_for_cycle(pvtk_map, cycle)
    if pvtk_file is None:
        points = np.zeros((0, 3), dtype=float)
        vel = np.zeros((0, 3), dtype=float)
    else:
        pdata = read_particle_vtk(pvtk_file)
        points = np.asarray(pdata.points, dtype=float)
        if "vel" in pdata.vectors:
            vel = np.asarray(pdata.vectors["vel"], dtype=float)
        elif pdata.vectors:
            vname = sorted(pdata.vectors.keys())[0]
            vel = np.asarray(pdata.vectors[vname], dtype=float)
        else:
            vel = np.zeros((0, 3), dtype=float)

    x_spec_edges = np.linspace(float(x1f[0]), float(x1f[-1]), args.spec_xbins + 1)
    y_spec_edges = np.linspace(args.spec_ymin, args.spec_ymax, args.spec_ybins + 1)
    spec_proxy = _spectrum_proxy(points, vel, args.u0, x_spec_edges, y_spec_edges)
    spec_plot = spec_proxy.T
    spec_plot[spec_plot <= 0.0] = np.nan

    x_lines = _meshblock_lines(x1f, args.mb_nx1)
    y_lines = _meshblock_lines(x2f, args.mb_nx2)

    fig = plt.figure(figsize=(13.8, 8.0))
    gs = fig.add_gridspec(
        3,
        2,
        width_ratios=[50.0, 0.8],
        height_ratios=[1.0, 1.0, 1.05],
        wspace=0.06,
        hspace=0.12,
    )

    ax1 = fig.add_subplot(gs[0, 0])
    cax1 = fig.add_subplot(gs[0, 1])
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
    ax1.axvline(shock_x, color="k", lw=0.8, ls="--", alpha=0.6)
    ax1.set_xlim(float(x1f[0]), float(x1f[-1]))
    ax1.set_ylim(float(x2f[0]), float(x2f[-1]))
    ax1.set_ylabel(r"$y$", fontsize=13)
    tval = float(rho_data.get("Time", 0.0))
    ax1.set_title(rf"$\Omega_0 t={tval:.3e}$", fontsize=20, pad=4)
    cbar1 = fig.colorbar(m1, cax=cax1)
    cbar1.set_label(r"$\rho/\rho_0$", fontsize=14)

    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    cax2 = fig.add_subplot(gs[1, 1])
    m2 = ax2.pcolormesh(
        x1f,
        x2f,
        b_norm,
        shading="auto",
        cmap="jet",
        norm=colors.LogNorm(vmin=args.b_vmin, vmax=args.b_vmax),
    )
    _draw_meshblock_grid(ax2, x_lines, y_lines)
    ax2.axvline(shock_x, color="k", lw=0.8, ls="--", alpha=0.6)
    ax2.set_ylabel(r"$y$", fontsize=13)
    cbar2 = fig.colorbar(m2, cax=cax2)
    cbar2.set_label(r"$B/B_0$", fontsize=14)

    ax3 = fig.add_subplot(gs[2, 0], sharex=ax1)
    cax3 = fig.add_subplot(gs[2, 1])
    m3 = ax3.pcolormesh(
        x_spec_edges,
        y_spec_edges,
        spec_plot,
        shading="auto",
        cmap="jet",
        norm=colors.LogNorm(vmin=args.spec_vmin, vmax=args.spec_vmax),
    )
    ax3.axvline(shock_x, color="k", lw=0.8, ls="--", alpha=0.6)
    ax3.set_ylim(args.spec_ymin, args.spec_ymax)
    ax3.set_ylabel(r"$\log_{10}\!\left(2\epsilon/(m u_0^2)\right)$", fontsize=14)
    ax3.set_xlabel(r"$x$", fontsize=22)
    if points.size == 0:
        ax3.text(
            0.5,
            0.5,
            "No particle data found",
            transform=ax3.transAxes,
            ha="center",
            va="center",
            fontsize=12,
        )
    cbar3 = fig.colorbar(m3, cax=cax3)
    cbar3.set_label(r"$\epsilon f(\epsilon)$ (proxy)", fontsize=14)

    fig_path = out_dir / args.figure_name
    fig.text(
        0.995, 0.005, _PROXY_WATERMARK, ha="right", va="bottom",
        fontsize=7, color="#8b0000",
    )
    fig.savefig(fig_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "figure": str(fig_path),
        "evidence_class": "engineering_proxy",
        "not_sun_bai_reproduction": True,
        "qualification_status": "unqualified_engineering_quicklook",
        "physical_model": "pic_parallel_shock_engineering_scaffold",
        "visible_watermark": _PROXY_WATERMARK,
        "basename": args.basename,
        "rho_file": str(rho_file),
        "rho_file_sha256": _sha256(rho_file),
        "bmag_file": str(bmag_file),
        "bmag_file_sha256": _sha256(bmag_file),
        "cycle": int(cycle),
        "pvtk_file": "" if pvtk_file is None else str(pvtk_file),
        "pvtk_file_sha256": "" if pvtk_file is None else _sha256(pvtk_file),
        "pvtk_cycle": None if pvtk_cycle is None else int(pvtk_cycle),
        "figure_sha256": _sha256(fig_path),
        "time": tval,
        "shock_x": float(shock_x),
        "rho_ref": float(rho_ref),
        "b_ref": float(b_ref),
        "n_particles": int(points.shape[0]),
        "u0": float(args.u0),
        "spec_xbins": int(args.spec_xbins),
        "spec_ybins": int(args.spec_ybins),
        "spec_ymin": float(args.spec_ymin),
        "spec_ymax": float(args.spec_ymax),
    }
    quicklook_summary = out_dir / "engineering_quicklook_summary.json"
    quicklook_summary.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    inputs = [rho_file, bmag_file]
    if pvtk_file is not None:
        inputs.append(pvtk_file)
    write_companion_manifest(
        out_dir / "engineering_quicklook_artifact_manifest.json",
        generator="plot_pic_parallel_shock_section54.py",
        physical_model="pic_parallel_shock_engineering_scaffold",
        inputs=inputs,
        outputs=(fig_path, quicklook_summary),
    )
    print(f"figure: {fig_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
