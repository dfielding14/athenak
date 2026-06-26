#!/usr/bin/env python3
"""Render the linear-to-saturated Bell evolution from one accepted engineering run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

from tst.publication.analyze_q011_section54_outputs import (
    compose_leaf_field,
    read_athenak_binary,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "sha256": _sha256(path),
        "byte_count": path.stat().st_size,
    }


def _style() -> None:
    plt.style.use("default")
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#202020",
            "axes.linewidth": 0.9,
            "axes.labelsize": 10.5,
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
            "mathtext.fontset": "dejavuserif",
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "legend.fontsize": 8.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.10,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        clip_on=False,
    )


def _snapshot_path(run_root: Path, basename: str, index: int) -> Path:
    return run_root / "bin" / f"{basename}.mhd_w_bcc.{index:05d}.bin"


def _map(path: Path) -> dict[str, Any]:
    dataset = read_athenak_binary(path)
    by = np.asarray(compose_leaf_field(dataset, "bcc2").values, dtype=np.float64)
    bz = np.asarray(compose_leaf_field(dataset, "bcc3").values, dtype=np.float64)
    density = np.asarray(compose_leaf_field(dataset, "dens").values, dtype=np.float64)
    bperp = np.sqrt(by * by + bz * bz)
    return {
        "time": float(dataset.time),
        "bperp": bperp[bperp.shape[0] // 2],
        "density": density[density.shape[0] // 2],
        "x_faces": np.linspace(
            dataset.domain_bounds[0], dataset.domain_bounds[1], bperp.shape[2] + 1
        ),
        "y_faces": np.linspace(
            dataset.domain_bounds[2], dataset.domain_bounds[3], bperp.shape[1] + 1
        ),
    }


def _nearest(rows: list[dict[str, Any]], key: str, target: float) -> dict[str, Any]:
    return min(rows, key=lambda row: abs(float(row[key]) - target))


def make_figure(run_root: Path, output_root: Path, dpi: int) -> list[Path]:
    analysis_path = run_root / "saturation_analysis.json"
    report = json.loads(analysis_path.read_text(encoding="utf-8"))
    if report.get("saturation_ready") is not True:
        raise RuntimeError("refusing to label a run without a passing sustained-saturation gate")
    rows = report["snapshots"]
    basename = str(report["basename"])
    onset = report["nonlinear_onset"]
    confirmation = report["sustained_confirmation"]
    linear_row = _nearest(rows, "bperp_rms_over_B0", 0.05)
    onset_row = rows[int(onset["index"])]
    saturated_row = _nearest(
        rows,
        "tau",
        0.5 * (float(confirmation["start_tau"]) + float(confirmation["end_tau"])),
    )
    selected_rows = (linear_row, onset_row, saturated_row)
    selected_paths = [
        _snapshot_path(run_root, basename, int(row["output_index"]))
        for row in selected_rows
    ]
    maps = [_map(path) for path in selected_paths]

    tau = np.asarray([row["tau"] for row in rows], dtype=np.float64)
    bperp = np.asarray([row["bperp_rms_over_B0"] for row in rows], dtype=np.float64)
    current = np.asarray([row["mean_jx_over_initial"] for row in rows], dtype=np.float64)
    # The t=0 moment file precedes the first particle deposition.
    if current.size > 1 and current[0] < 0.5 and current[1] > 0.9:
        current[0] = np.nan
    _style()
    fig = plt.figure(figsize=(7.15, 4.70), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=(1.35, 0.62))
    ax = fig.add_subplot(grid[0, :])
    amplitude_line = ax.semilogy(
        tau, bperp, color="#D55E00", lw=2.0, label=r"$B_{\perp,\rm rms}/B_0$"
    )
    ax.axhline(1.0, color="#555555", ls="--", lw=1.0, label="nonlinear onset")
    ax.axvspan(
        float(confirmation["start_tau"]),
        float(confirmation["end_tau"]),
        color="#009E73",
        alpha=0.14,
        lw=0.0,
        label="sustained plateau",
    )
    for row, marker in zip(selected_rows, ("o", "s", "D")):
        ax.scatter(
            [row["tau"]],
            [row["bperp_rms_over_B0"]],
            marker=marker,
            s=30,
            color="#111111",
            edgecolor="white",
            linewidth=0.5,
            zorder=5,
        )
    ax.set_xlabel(r"$\tau=k_0u_{A0}t$")
    ax.set_ylabel(r"$B_{\perp,\rm rms}/B_0$", color="#D55E00")
    ax.tick_params(axis="y", labelcolor="#D55E00")
    ax.set_xlim(tau[0], tau[-1])
    current_ax = ax.twinx()
    current_line = current_ax.plot(
        tau, current, color="#0072B2", lw=1.7, label=r"$\langle J_{\rm CR,x}\rangle/J_{\rm CR,x}(0)$"
    )
    current_ax.set_ylabel(r"$\langle J_{\rm CR,x}\rangle/J_{\rm CR,x}(0)$", color="#0072B2")
    current_ax.tick_params(axis="y", labelcolor="#0072B2")
    current_ax.set_ylim(min(0.0, 1.05 * float(np.nanmin(current))), 1.05)
    handles = amplitude_line + current_line
    handles.extend(ax.get_legend_handles_labels()[0][1:])
    labels = [item.get_label() for item in handles]
    ax.legend(handles, labels, loc="lower left", ncol=2)
    ax.text(
        0.99,
        0.96,
        rf"peak $B_\perp/B_0={report['peak_bperp_rms_over_B0']:.2f}$"
        + "\n"
        + rf"final current loss $={100.0*report['final_current_relaxation_fraction']:.1f}\%$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.2,
    )
    _panel_label(ax, "(a)")

    positive = np.concatenate([item["bperp"].ravel() for item in maps])
    positive = positive[positive > 0.0]
    norm = colors.LogNorm(
        vmin=max(1.0e-3, float(np.percentile(positive, 0.1))),
        vmax=float(np.percentile(positive, 99.9)),
    )
    titles = ("linear", "nonlinear onset", "saturated")
    image = None
    for column, (data, row, title) in enumerate(zip(maps, selected_rows, titles)):
        ax = fig.add_subplot(grid[1, column])
        image = ax.pcolormesh(
            data["x_faces"],
            data["y_faces"],
            data["bperp"],
            shading="flat",
            cmap="magma",
            norm=norm,
            rasterized=True,
        )
        density = data["density"]
        if float(np.min(density)) < 0.8 or float(np.max(density)) > 1.2:
            x = 0.5 * (data["x_faces"][:-1] + data["x_faces"][1:])
            y = 0.5 * (data["y_faces"][:-1] + data["y_faces"][1:])
            levels = [level for level in (0.5, 1.5) if np.min(density) < level < np.max(density)]
            if levels:
                ax.contour(x, y, density, levels=levels, colors="white", linewidths=0.55)
        ax.set_title(
            title + "\n" + rf"$\tau={row['tau']:.1f}$, $B_{{\perp,\rm rms}}/B_0={row['bperp_rms_over_B0']:.2f}$",
            fontsize=9.0,
        )
        ax.set_xlabel(r"$x/\lambda_0$")
        if column == 0:
            ax.set_ylabel(r"$y/\lambda_0$")
        else:
            ax.set_yticklabels([])
        ax.set_aspect("equal")
        ax.text(
            0.025,
            0.95,
            f"({chr(ord('b') + column)})",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10.5,
            fontweight="bold",
            color="white",
            bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.45, "pad": 1.5},
        )
    if image is not None:
        figure_bar = fig.colorbar(image, ax=fig.axes[2:5], orientation="horizontal", shrink=0.72, pad=0.02)
        figure_bar.set_label(r"$B_\perp/B_0$ at $z=L_z/2$")

    output_root.mkdir(parents=True, exist_ok=True)
    outputs = [
        output_root / "figure_bell_nonlinear_saturation.pdf",
        output_root / "figure_bell_nonlinear_saturation.png",
    ]
    fig.savefig(outputs[0], bbox_inches="tight")
    fig.savefig(outputs[1], dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    caption = output_root / "figure_bell_nonlinear_saturation_caption.md"
    caption.write_text(
        "**Figure. Linear growth, nonlinear transition, and sustained saturation of "
        "the Bell instability in a compact three-dimensional finite-rigidity "
        "MHD-PIC calculation.** (a) Transverse rms magnetic amplitude and mean "
        "parallel cosmic-ray current. The green interval is the independently "
        "specified sustained-plateau window. (b-d) Midplane transverse-field "
        "amplitude during linear growth, at the first nonlinear crossing, and in "
        "the sustained saturated interval. White contours, where present, mark "
        "density contrasts of 0.5 and 1.5 relative to the initial density. The run "
        "uses volume-aware current normalization and finite particle rigidity. It "
        "is compact single-run engineering evidence, not a box-size or resolution "
        "convergence claim.\n",
        encoding="utf-8",
    )
    outputs.append(caption)

    raw_manifest = run_root / "raw_outputs.sha256"
    raw_inputs = []
    for line in raw_manifest.read_text(encoding="utf-8").splitlines():
        expected, name = line.split(maxsplit=1)
        path = Path(name)
        artifact = _artifact(path)
        if artifact["sha256"] != expected:
            raise RuntimeError(f"raw output hash mismatch: {path}")
        raw_inputs.append(artifact)
    manifest = {
        "record_type": "q019_bell_saturation_publication_figure_manifest_v1",
        "source_commit": (run_root / "source.commit").read_text(encoding="utf-8").strip(),
        "generator": _artifact(Path(__file__)),
        "inputs": [
            _artifact(analysis_path),
            _artifact(run_root / "bindings.sha256"),
            _artifact(raw_manifest),
            *raw_inputs,
        ],
        "outputs": [_artifact(path) for path in outputs],
        "metrics": {
            key: report[key]
            for key in (
                "snapshot_count",
                "nonlinear_onset",
                "candidate_plateau",
                "sustained_confirmation",
                "peak_bperp_rms_over_B0",
                "final_bperp_rms_over_B0",
                "final_current_relaxation_fraction",
                "saturation_ready",
            )
        },
        "scope": report["claim_scope"],
    }
    manifest_path = output_root / "figure_bell_nonlinear_saturation_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return [*outputs, manifest_path]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--dpi", type=int, default=400)
    args = parser.parse_args()
    paths = make_figure(
        args.run_root.resolve(strict=True), args.output_root.resolve(), args.dpi
    )
    print(json.dumps({"outputs": [str(path) for path in paths]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
