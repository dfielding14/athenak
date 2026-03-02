#!/usr/bin/env python3
"""Create diagnostics plots for a TRML run directory."""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from typing import List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np


GAMMA = 1.666666667
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BIN_READER_DIR = os.path.join(REPO_ROOT, "vis", "python")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--log-file", default="run.log")
    parser.add_argument("--out-dir", default="plots")
    parser.add_argument("--slice-glob", default="bin/TRML.slice_y.000*.bin")
    parser.add_argument("--full-glob", default="bin/TRML.hydro_w.000*.bin")
    return parser.parse_args()


def _parse_log_metrics(log_path: str) -> dict:
    current_time = None
    elapsed_re = re.compile(r"time=([0-9eE+\-\.]+)")
    metric_res = {
        "W_global": re.compile(r"W_global=([0-9eE+\-\.]+)"),
        "mean_z_ft": re.compile(r"mean_z_ft=([0-9eE+\-\.]+)"),
        "mean_velz_ft": re.compile(r"mean_velz_ft=([0-9eE+\-\.]+)"),
        "vel_z_shift": re.compile(r"vel_z_shift=([0-9eE+\-\.]+)"),
    }
    out = {k: {"t": [], "v": []} for k in metric_res}

    with open(log_path, "r", encoding="utf-8") as f:
        for line in f:
            m = elapsed_re.search(line)
            if m:
                current_time = float(m.group(1))
                continue
            if current_time is None:
                continue
            for name, rx in metric_res.items():
                mm = rx.search(line)
                if mm:
                    out[name]["t"].append(current_time)
                    out[name]["v"].append(float(mm.group(1)))
                    break
    return out


def _load_hst(path: str) -> np.ndarray:
    arr = np.loadtxt(path, comments="#")
    if arr.ndim == 1:
        arr = arr[None, :]
    return arr


def _clean_bad(vals: np.ndarray) -> np.ndarray:
    vals = vals.astype(float)
    vals[np.abs(vals) > 1.0e100] = np.nan
    return vals


def _snapshot_number(path: str) -> int:
    m = re.search(r"\.(\d+)\.bin$", os.path.basename(path))
    return int(m.group(1)) if m else -1


def _pick_three(files: List[str]) -> List[str]:
    if len(files) <= 3:
        return files
    idx = sorted(set([0, len(files) // 2, len(files) - 1]))
    return [files[i] for i in idx]


def _to_2d(arr: np.ndarray) -> np.ndarray:
    out = np.squeeze(np.asarray(arr, dtype=float))
    if out.ndim != 2:
        raise ValueError(f"Expected 2D array after squeeze, got shape={out.shape}")
    return out


def _get_vmin_vmax(
    fields: Sequence[np.ndarray], symmetric: bool = False
) -> Tuple[Optional[float], Optional[float]]:
    finite_vals = [f[np.isfinite(f)] for f in fields if np.any(np.isfinite(f))]
    if not finite_vals:
        return None, None
    all_vals = np.concatenate(finite_vals)
    if symmetric:
        amax = float(np.max(np.abs(all_vals)))
        if amax == 0.0:
            amax = 1.0
        return -amax, amax
    return float(np.min(all_vals)), float(np.max(all_vals))


def _plot_three_panel(
    out_path: str,
    fields: Sequence[np.ndarray],
    metas: Sequence[Tuple[np.ndarray, np.ndarray, float, str]],
    cmap: str,
    cbar_label: str,
    symmetric: bool = False,
) -> None:
    if not fields:
        return
    vmin, vmax = _get_vmin_vmax(fields, symmetric=symmetric)
    fig, axs = plt.subplots(
        1, len(fields), figsize=(5 * len(fields), 4.6), constrained_layout=True
    )
    if len(fields) == 1:
        axs = [axs]
    for ax, field, (x1f, x3f, tcur, bname) in zip(axs, fields, metas):
        im = ax.imshow(
            field,
            origin="lower",
            extent=[x1f[0], x1f[-1], x3f[0], x3f[-1]],
            aspect="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(f"{bname}\nt={tcur:.3f}")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
    cbar = fig.colorbar(im, ax=axs, shrink=0.95)
    cbar.set_label(cbar_label)
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


def _project_density_along_y(dens: np.ndarray, x2f: np.ndarray) -> np.ndarray:
    dens3d = np.asarray(dens, dtype=float)
    if dens3d.ndim != 3:
        dens3d = np.squeeze(dens3d)
    if dens3d.ndim != 3:
        raise ValueError(f"Expected 3D density for projection, got shape={dens3d.shape}")
    dy = np.diff(np.asarray(x2f, dtype=float))
    if dens3d.shape[1] != dy.size:
        raise ValueError(
            f"y-size mismatch: dens.shape[1]={dens3d.shape[1]} but len(dy)={dy.size}"
        )
    rho_col = np.sum(np.maximum(dens3d, 0.0) * dy[np.newaxis, :, np.newaxis], axis=1)
    return rho_col


def main() -> None:
    args = _parse_args()
    run_dir = os.path.abspath(args.run_dir)
    out_dir = os.path.join(run_dir, args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    log_path = os.path.join(run_dir, args.log_file)
    user_hst = os.path.join(run_dir, "history", "TRML.user.hst")
    hydro_hst = os.path.join(run_dir, "history", "TRML.hydro.hst")
    slice_glob = os.path.join(run_dir, args.slice_glob)
    full_glob = os.path.join(run_dir, args.full_glob)

    metrics = _parse_log_metrics(log_path)
    user = _load_hst(user_hst)
    hydro = _load_hst(hydro_hst)

    # Frame tracking timeseries
    fig, axs = plt.subplots(4, 1, figsize=(10, 10), sharex=True, constrained_layout=True)
    for ax, key in zip(
        axs, ["mean_z_ft", "mean_velz_ft", "vel_z_shift", "W_global"]
    ):
        t = np.asarray(metrics[key]["t"])
        v = np.asarray(metrics[key]["v"])
        if key == "W_global":
            ax.semilogy(t, v, lw=1.3)
            ax.set_ylabel("W_global")
        else:
            ax.plot(t, v, lw=1.3)
            ax.set_ylabel(key)
        ax.grid(alpha=0.3)
    axs[2].axhline(0.01, color="r", ls="--", lw=0.8, alpha=0.6)
    axs[2].axhline(-0.01, color="r", ls="--", lw=0.8, alpha=0.6)
    axs[0].set_title(f"Frame-tracking: {os.path.basename(run_dir)}")
    axs[-1].set_xlabel("time")
    fig.savefig(os.path.join(out_dir, "frame_tracking_timeseries.png"), dpi=170)
    plt.close(fig)

    # Interface z-ranges
    t = user[:, 0]
    zmin_shear = _clean_bad(user[:, 35])
    zmax_shear = _clean_bad(user[:, 36])
    zmin_peak = _clean_bad(user[:, 37])
    zmax_peak = _clean_bad(user[:, 38])
    zmin_vy = _clean_bad(user[:, 39])
    zmax_vy = _clean_bad(user[:, 40])

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    ax.fill_between(t, zmin_peak, zmax_peak, alpha=0.35, label="peak band z-range")
    ax.plot(t, zmin_shear, lw=1.0, label="zmin shear")
    ax.plot(t, zmax_shear, lw=1.0, label="zmax shear")
    ax.plot(t, zmin_vy, lw=1.0, ls="--", label="zmin vely")
    ax.plot(t, zmax_vy, lw=1.0, ls="--", label="zmax vely")
    ax.axhline(0.0, color="k", lw=0.8, alpha=0.6)
    ax.set_xlabel("time")
    ax.set_ylabel("z")
    ax.set_title("Interface position diagnostics")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=9)
    fig.savefig(os.path.join(out_dir, "interface_z_ranges.png"), dpi=170)
    plt.close(fig)

    # Hydro history
    th = hydro[:, 0]
    mass = hydro[:, 2]
    momz = hydro[:, 5]
    etot = hydro[:, 6]
    fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    ax[0].plot(th, mass, lw=1.2)
    ax[0].set_ylabel("mass")
    ax[0].grid(alpha=0.3)
    ax[1].plot(th, momz, lw=1.2)
    ax[1].set_ylabel("z-momentum")
    ax[1].grid(alpha=0.3)
    ax[2].plot(th, etot, lw=1.2)
    ax[2].set_ylabel("total E")
    ax[2].set_xlabel("time")
    ax[2].grid(alpha=0.3)
    fig.savefig(os.path.join(out_dir, "hydro_history.png"), dpi=170)
    plt.close(fig)

    # Three temperature slices, velocity slices, and 3D density projection
    sys.path.insert(0, BIN_READER_DIR)
    from bin_convert_new import read_binary_as_athdf  # pylint: disable=import-error

    slice_files = sorted(glob.glob(slice_glob), key=_snapshot_number)
    chosen = _pick_three(slice_files)
    if chosen:
        temps = []
        velx_fields = []
        vely_fields = []
        velz_fields = []
        meta = []
        for sf in chosen:
            d = read_binary_as_athdf(sf)
            dens = _to_2d(np.asarray(d["dens"], dtype=float))
            eint = _to_2d(np.asarray(d["eint"], dtype=float))
            temp = (GAMMA - 1.0) * eint / np.maximum(dens, 1.0e-30)
            temps.append(np.log10(np.maximum(temp, 1.0e-30)))
            velx_fields.append(_to_2d(np.asarray(d["velx"], dtype=float)))
            vely_fields.append(_to_2d(np.asarray(d["vely"], dtype=float)))
            velz_fields.append(_to_2d(np.asarray(d["velz"], dtype=float)))
            meta.append(
                (
                    np.asarray(d["x1f"], dtype=float),
                    np.asarray(d["x3f"], dtype=float),
                    float(d["Time"]),
                    os.path.basename(sf),
                )
            )
        _plot_three_panel(
            os.path.join(out_dir, "temperature_slices_3panel.png"),
            temps,
            meta,
            cmap="inferno",
            cbar_label="log10(T)",
            symmetric=False,
        )
        _plot_three_panel(
            os.path.join(out_dir, "velx_slices_3panel.png"),
            velx_fields,
            meta,
            cmap="coolwarm",
            cbar_label="v_x",
            symmetric=True,
        )
        _plot_three_panel(
            os.path.join(out_dir, "vely_slices_3panel.png"),
            vely_fields,
            meta,
            cmap="coolwarm",
            cbar_label="v_y",
            symmetric=True,
        )
        _plot_three_panel(
            os.path.join(out_dir, "velz_slices_3panel.png"),
            velz_fields,
            meta,
            cmap="coolwarm",
            cbar_label="v_z",
            symmetric=True,
        )

    full_files = sorted(glob.glob(full_glob), key=_snapshot_number)
    chosen_full = _pick_three(full_files)
    if chosen_full:
        rho_col_logs = []
        full_meta = []
        for ff in chosen_full:
            d = read_binary_as_athdf(ff)
            rho_col = _project_density_along_y(d["dens"], d["x2f"])
            rho_col_logs.append(np.log10(np.maximum(rho_col, 1.0e-30)))
            full_meta.append(
                (
                    np.asarray(d["x1f"], dtype=float),
                    np.asarray(d["x3f"], dtype=float),
                    float(d["Time"]),
                    os.path.basename(ff),
                )
            )
        _plot_three_panel(
            os.path.join(out_dir, "density_projection_y_3panel.png"),
            rho_col_logs,
            full_meta,
            cmap="magma",
            cbar_label="log10(∫ρ dy)",
            symmetric=False,
        )

    print(out_dir)


if __name__ == "__main__":
    main()
