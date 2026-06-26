#!/usr/bin/env python3
"""Create publication-style Bell and shock figures from retained raw outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

from tst.publication import analyze_q011_section54_outputs as binary
from tst.publication import analyze_q023_paper_bell_linear_joverc as linear
from tst.publication import q019_registered_raw_reduction_v1 as nonlinear_raw
from tst.publication.pvtk_particles import read_particle_vtk
from tst.publication.q011_section54_particles import reconstruct_chi_from_pvtk_velocity


LINEAR_MEMBER_ID = "d1-fine-reference_x1-eps0p4"
NONLINEAR_CASE_ID = "q019-q023-carrier-s2-resolution-coarse-s0"
SHOCK_INJECTION_MOMENTUM = 30.0 * 3.16227766017
COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")


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
            "axes.titlesize": 10.5,
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
            "mathtext.fontset": "dejavuserif",
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.minor.visible": True,
            "ytick.minor.visible": True,
            "legend.frameon": False,
            "legend.fontsize": 8.0,
            "lines.linewidth": 1.8,
            "axes.unicode_minus": True,
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
        color="#111111",
        clip_on=False,
    )


def _write_figure(fig: plt.Figure, output_root: Path, stem: str, dpi: int) -> list[Path]:
    paths = [output_root / f"{stem}.pdf", output_root / f"{stem}.png"]
    fig.savefig(paths[0], bbox_inches="tight")
    fig.savefig(paths[1], dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return paths


def _linear_trace(run_root: Path) -> tuple[dict[str, Any], list[Path]]:
    member = next(
        item
        for item in linear.expected_deck_members()
        if item["member_id"] == LINEAR_MEMBER_ID
    )
    basename = "q023_joverc_" + LINEAR_MEMBER_ID.replace("-", "_")
    paths = sorted((run_root / "bin").glob(f"{basename}.mhd_w_bcc.*.bin"))
    datasets = [binary.read_athenak_binary(path) for path in paths]
    trace = linear.physics_trace_from_raw_datasets(datasets, member=member)
    report = linear._analyze_physics_trace(trace, float(member["epsilon"]))
    return {"member": member, "trace": trace, "report": report}, paths


def _nonlinear_snapshot(
    run_root: Path,
) -> tuple[dict[str, Any], dict[str, np.ndarray], tuple[np.ndarray, ...], list[Path]]:
    summary_path = run_root / "quicklook.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    peak = max(summary["snapshots"], key=lambda row: row["bperp_rms_over_B0"])
    index = int(peak["output_index"])
    products: dict[str, bytes] = {}
    paths: list[Path] = [summary_path]
    for product in nonlinear_raw.REQUIRED_BINARY_PRODUCTS:
        path = (
            run_root
            / "bin"
            / f"{NONLINEAR_CASE_ID}.{product}.{index:05d}.bin"
        )
        products[product] = path.read_bytes()
        paths.append(path)
    snapshot = nonlinear_raw.compose_snapshot(NONLINEAR_CASE_ID, products)
    fields = {
        name: np.asarray(values, dtype=np.float64)
        for name, values in snapshot["fields"].items()
    }
    faces = tuple(
        np.asarray(snapshot[name], dtype=np.float64)
        for name in ("x1_faces", "x2_faces", "x3_faces")
    )
    return {"summary": summary, "peak": peak}, fields, faces, paths


def make_bell_figure(
    linear_root: Path, nonlinear_root: Path, output_root: Path, dpi: int
) -> tuple[list[Path], dict[str, Any], list[Path]]:
    linear_data, linear_paths = _linear_trace(linear_root)
    nonlinear_data, fields, faces, nonlinear_paths = _nonlinear_snapshot(nonlinear_root)
    trace = linear_data["trace"]
    report = linear_data["report"]

    time = np.asarray(trace["normalized_time"], dtype=np.float64)
    right = np.asarray(trace["right_mode_real"]) + 1j * np.asarray(
        trace["right_mode_imag"]
    )
    amplitude = 2.0 * np.abs(right)
    phase = np.unwrap(np.angle(right))
    phase -= phase[0]
    expected_growth = float(report["expected_growth_rate_over_k0_ua"])
    expected_phase = float(report["expected_phase_frequency_over_k0_ua"])

    history = nonlinear_data["summary"]["snapshots"]
    nonlinear_time = np.asarray([row["time"] for row in history], dtype=np.float64)
    nonlinear_tau = 2.0 * np.pi * nonlinear_time
    bperp_history = np.asarray(
        [row["bperp_rms_over_B0"] for row in history], dtype=np.float64
    )
    btotal_history = np.asarray(
        [row["btotal_rms_over_B0"] for row in history], dtype=np.float64
    )
    peak = nonlinear_data["peak"]
    b0 = float(nonlinear_data["summary"]["initial_B0"])
    by_map = fields["bcc2"][0] / b0
    x1_faces, x2_faces, _ = faces

    fig, axes = plt.subplots(2, 2, figsize=(7.15, 6.15), constrained_layout=True)
    ax = axes[0, 0]
    ax.semilogy(time, amplitude, color=COLORS[0], label="PIC")
    ax.semilogy(
        time,
        amplitude[0] * np.exp(expected_growth * time),
        color="#222222",
        ls="--",
        label="linear theory",
    )
    ax.axvspan(
        1.0 / expected_growth,
        5.0 / expected_growth,
        color="#7f7f7f",
        alpha=0.13,
        lw=0.0,
        label="fit interval",
    )
    ax.set_xlabel(r"$k_0 u_{A0} t$")
    ax.set_ylabel(r"$2|\delta B_+(k_0)|/B_0$")
    ax.set_xlim(time[0], time[-1])
    ax.legend(loc="upper left")
    ax.text(
        0.98,
        0.05,
        rf"$\gamma_{{\rm PIC}}={report['measured_growth_rate_over_k0_ua']:.3f}$"
        + "\n"
        + rf"$\gamma_{{\rm th}}={expected_growth:.3f}$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
    )
    _panel_label(ax, "(a)")

    ax = axes[0, 1]
    ax.plot(time, phase, color=COLORS[0], label="PIC")
    ax.plot(time, expected_phase * time, color="#222222", ls="--", label="linear theory")
    ax.set_xlabel(r"$k_0 u_{A0} t$")
    ax.set_ylabel(r"$\arg B_+(t)-\arg B_+(0)$")
    ax.set_xlim(time[0], time[-1])
    ax.legend(loc="lower left")
    ax.text(
        0.98,
        0.95,
        rf"$\omega_{{\rm PIC}}={report['measured_phase_frequency_over_k0_ua']:.3f}$"
        + "\n"
        + rf"$\omega_{{\rm th}}={expected_phase:.3f}$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
    )
    _panel_label(ax, "(b)")

    ax = axes[1, 0]
    ax.semilogy(nonlinear_tau, bperp_history, color=COLORS[1], label=r"$B_{\perp,\rm rms}$")
    ax.semilogy(nonlinear_tau, btotal_history, color=COLORS[2], ls=":", label=r"$B_{\rm rms}$")
    ax.axhline(1.0, color="#555555", lw=1.0, ls="--")
    ax.scatter(
        [2.0 * np.pi * peak["time"]],
        [peak["bperp_rms_over_B0"]],
        s=26,
        color=COLORS[1],
        edgecolor="white",
        linewidth=0.5,
        zorder=4,
    )
    ax.set_xlabel(r"$\tau=k_0u_{A0}t$")
    ax.set_ylabel(r"$B/B_0$")
    ax.set_xlim(nonlinear_tau[0], nonlinear_tau[-1])
    ax.set_ylim(5.0e-5, 3.0)
    ax.legend(loc="upper left")
    ax.text(
        0.98,
        0.95,
        rf"peak $B_{{\perp,\rm rms}}/B_0={peak['bperp_rms_over_B0']:.2f}$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
    )
    _panel_label(ax, "(c)")

    ax = axes[1, 1]
    vmax = float(np.percentile(np.abs(by_map), 99.5))
    mesh = ax.pcolormesh(
        x1_faces,
        x2_faces,
        by_map,
        shading="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        rasterized=True,
    )
    ax.set_xlabel(r"$x/\lambda_0$")
    ax.set_ylabel(r"$y/\lambda_0$")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(rf"$\tau={2.0 * np.pi * peak['time']:.2f}$", pad=3)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$B_y/B_0$")
    cbar.ax.tick_params(labelsize=8)
    _panel_label(ax, "(d)")

    paths = _write_figure(fig, output_root, "figure_bell_instability", dpi)
    metrics = {
        "linear_snapshot_count": len(linear_paths),
        "linear_measured_growth": report["measured_growth_rate_over_k0_ua"],
        "linear_expected_growth": expected_growth,
        "linear_measured_phase": report["measured_phase_frequency_over_k0_ua"],
        "linear_expected_phase": expected_phase,
        "nonlinear_snapshot_count": len(history),
        "nonlinear_peak_time": peak["time"],
        "nonlinear_peak_tau": 2.0 * np.pi * peak["time"],
        "nonlinear_peak_bperp_rms_over_B0": peak["bperp_rms_over_B0"],
    }
    return paths, metrics, [*linear_paths, *nonlinear_paths]


def _shock_front(x: np.ndarray, density_profile: np.ndarray) -> float:
    gradient = np.gradient(density_profile, x)
    mask = (x > x[0] + 0.05 * (x[-1] - x[0])) & (
        x < x[0] + 0.90 * (x[-1] - x[0])
    )
    candidates = np.flatnonzero(mask)
    return float(x[candidates[np.argmin(gradient[candidates])]])


def _particle_state(run_root: Path, row: dict[str, Any]) -> tuple[dict[str, np.ndarray], Path]:
    path = run_root / row["path"]
    data = read_particle_vtk(path)
    mask = np.asarray(data.scalars["cr_source"]) == 1
    momentum = 30.0 * np.sqrt(
        reconstruct_chi_from_pvtk_velocity(np.asarray(data.vectors["vel"])[mask])
    )
    return {
        "x": np.asarray(data.points, dtype=np.float64)[mask, 0],
        "p": momentum,
        "weight": np.asarray(data.scalars["macro_weight"], dtype=np.float64)[mask],
        "ptag": np.asarray(data.scalars["ptag"], dtype=np.int64)[mask],
    }, path


def make_shock_figure(
    shock_root: Path, output_root: Path, dpi: int
) -> tuple[list[Path], dict[str, Any], list[Path]]:
    summary_path = shock_root / "quicklook.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    basename = str(summary["basename"])
    mhd_paths = sorted((shock_root / "bin").glob(f"{basename}.mhd_w_bcc.*.bin"))
    mhd_path = mhd_paths[-1]
    dataset = binary.read_athenak_binary(mhd_path)
    density_grid = binary.compose_leaf_field(dataset, "dens")
    density = np.asarray(density_grid.values[0], dtype=np.float64)
    x_faces = density_grid.x1_faces
    y_faces = density_grid.x2_faces
    x = 0.5 * (x_faces[:-1] + x_faces[1:])
    density_profile = np.mean(density, axis=0)
    shock_x = _shock_front(x, density_profile)

    nonempty_rows = [row for row in summary["snapshots"] if row["particle_count"] > 0]
    states_and_paths = [_particle_state(shock_root, row) for row in nonempty_rows]
    states = [item[0] for item in states_and_paths]
    particle_paths = [item[1] for item in states_and_paths]
    final = states[-1]

    fig, axes = plt.subplots(2, 2, figsize=(7.15, 6.15), constrained_layout=True)
    ax = axes[0, 0]
    mesh = ax.pcolormesh(
        x_faces,
        y_faces,
        density,
        shading="auto",
        cmap="viridis",
        vmin=0.9,
        vmax=float(np.percentile(density, 99.7)),
        rasterized=True,
    )
    ax.axvline(shock_x, color="white", lw=1.1, ls="--")
    ax.set_xlabel(r"$x/(c/\omega_{pi})$")
    ax.set_ylabel(r"$y/(c/\omega_{pi})$")
    ax.set_title(r"$\Omega_0 t=20$", pad=3)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$\rho/\rho_0$")
    cbar.ax.tick_params(labelsize=8)
    _panel_label(ax, "(a)")

    ax = axes[0, 1]
    x_edges = np.linspace(float(x_faces[0]), float(x_faces[-1]), 193)
    p_edges = np.linspace(0.88, 1.12, 121)
    histogram, _, _ = np.histogram2d(
        final["x"],
        final["p"] / SHOCK_INJECTION_MOMENTUM,
        bins=(x_edges, p_edges),
        weights=final["weight"],
    )
    positive = histogram[histogram > 0.0]
    phase = ax.pcolormesh(
        x_edges,
        p_edges,
        histogram.T,
        shading="auto",
        cmap="cividis",
        norm=colors.LogNorm(vmin=float(np.percentile(positive, 5.0)), vmax=float(positive.max())),
        rasterized=True,
    )
    ax.axvline(shock_x, color="white", lw=1.1, ls="--")
    ax.axhline(1.0, color="white", lw=0.8, ls=":")
    ax.set_xlabel(r"$x/(c/\omega_{pi})$")
    ax.set_ylabel(r"$p/p_{\rm inj}$")
    cbar = fig.colorbar(phase, ax=ax, pad=0.02)
    cbar.set_label("weighted counts")
    cbar.ax.tick_params(labelsize=8)
    _panel_label(ax, "(b)")

    ax = axes[1, 0]
    velocity_grid = binary.compose_leaf_field(dataset, "velx")
    pressure_grid = binary.compose_leaf_field(dataset, "eint")
    velocity_profile = np.mean(np.asarray(velocity_grid.values[0]), axis=0)
    pressure_profile = (5.0 / 3.0 - 1.0) * np.mean(
        np.asarray(pressure_grid.values[0]), axis=0
    )
    x_relative = x - shock_x
    shock_speed = 10.0
    upstream_shock_frame_speed = 40.0
    ax.plot(x_relative, density_profile, color=COLORS[0], label=r"$\rho/\rho_0$")
    ax.plot(
        x_relative,
        (shock_speed - velocity_profile) / upstream_shock_frame_speed,
        color=COLORS[1],
        label=r"$u_{\rm sh}/u_{\rm sh,0}$",
    )
    ax.plot(
        x_relative,
        pressure_profile / upstream_shock_frame_speed**2,
        color=COLORS[2],
        label=r"$P/(\rho_0u_{\rm sh,0}^2)$",
    )
    ax.axvline(0.0, color="#333333", lw=1.0, ls="--")
    ax.set_xlabel(r"$(x-x_{\rm sh})/(c/\omega_{pi})$")
    ax.set_ylabel("transverse average")
    ax.set_xlim(-120.0, 120.0)
    ax.set_ylim(-0.05, 4.25)
    ax.legend(loc="center right")
    _panel_label(ax, "(c)")

    ax = axes[1, 1]
    initial = states[0]
    initial_order = np.argsort(initial["ptag"])
    initial_tags = initial["ptag"][initial_order]
    initial_momentum = initial["p"][initial_order]
    cohort_time: list[float] = []
    cohort_survival: list[float] = []
    cohort_p10: list[float] = []
    cohort_p50: list[float] = []
    cohort_p90: list[float] = []
    cohort_p99: list[float] = []
    maximum_matched_gain = 0.0
    for row, state in zip(nonempty_rows, states):
        common, initial_index, state_index = np.intersect1d(
            initial_tags, state["ptag"], return_indices=True
        )
        delta_percent = 100.0 * (
            state["p"][state_index] / initial_momentum[initial_index] - 1.0
        )
        cohort_time.append(float(row["time"]))
        cohort_survival.append(float(common.size / initial_tags.size))
        p10, p50, p90, p99 = np.percentile(delta_percent, (10, 50, 90, 99))
        cohort_p10.append(float(p10))
        cohort_p50.append(float(p50))
        cohort_p90.append(float(p90))
        cohort_p99.append(float(p99))
        maximum_matched_gain = max(maximum_matched_gain, float(np.max(delta_percent)))
    cohort_time_array = np.asarray(cohort_time)
    ax.fill_between(
        cohort_time_array,
        cohort_p10,
        cohort_p90,
        color=COLORS[0],
        alpha=0.20,
        label="10th--90th percentile",
    )
    ax.plot(cohort_time_array, cohort_p50, color=COLORS[0], label="median")
    ax.plot(cohort_time_array, cohort_p99, color=COLORS[1], label="99th percentile")
    ax.axhline(0.0, color="#333333", lw=0.9, ls="--")
    ax.set_xlabel(r"$\Omega_0 t$")
    ax.set_ylabel(r"$100[p(t)/p(t_0)-1]$ (percent)")
    ax.set_xlim(cohort_time_array[0], cohort_time_array[-1])
    cohort_axis = ax.twinx()
    cohort_axis.plot(
        cohort_time_array,
        cohort_survival,
        color=COLORS[2],
        ls=":",
        label="surviving cohort",
    )
    cohort_axis.set_ylabel("surviving cohort fraction", color=COLORS[2])
    cohort_axis.tick_params(axis="y", colors=COLORS[2])
    cohort_axis.set_ylim(0.0, 1.05)
    lines, labels = ax.get_legend_handles_labels()
    other_lines, other_labels = cohort_axis.get_legend_handles_labels()
    ax.legend(lines + other_lines, labels + other_labels, loc="upper left", fontsize=7)
    ax.text(
        0.96,
        0.06,
        rf"maximum matched gain $={maximum_matched_gain:.3f}\%$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
    )
    _panel_label(ax, "(d)")

    # Shock-front trajectory is retained as a quantitative metric even though
    # the visible panel is reserved for exact particle-tag cohort retention.
    front_time: list[float] = []
    front_x: list[float] = []
    for path in mhd_paths:
        item = binary.read_athenak_binary(path)
        if item.time < 4.0:
            continue
        grid = binary.compose_leaf_field(item, "dens")
        centers = 0.5 * (grid.x1_faces[:-1] + grid.x1_faces[1:])
        profile = np.mean(np.asarray(grid.values[0], dtype=np.float64), axis=0)
        front_time.append(float(item.time))
        front_x.append(_shock_front(centers, profile))
    front_time_array = np.asarray(front_time)
    front_x_array = np.asarray(front_x)
    ideal_x = 10.0 * front_time_array
    maximum_front_error = float(np.max(np.abs(front_x_array - ideal_x)))

    downstream = (x_relative >= -100.0) & (x_relative <= -30.0)
    downstream_density = float(np.mean(density_profile[downstream]))
    downstream_velocity = float(np.mean(velocity_profile[downstream]))
    downstream_pressure = float(np.mean(pressure_profile[downstream]))

    paths = _write_figure(fig, output_root, "figure_nonrelativistic_shock", dpi)
    metrics = {
        "snapshot_count": summary["snapshot_count"],
        "final_time": summary["snapshots"][-1]["time"],
        "final_cycle": summary["snapshots"][-1]["cycle"],
        "final_particle_count": summary["snapshots"][-1]["particle_count"],
        "shock_x": shock_x,
        "maximum_shock_front_error": maximum_front_error,
        "downstream_density": downstream_density,
        "downstream_velocity_x": downstream_velocity,
        "downstream_pressure": downstream_pressure,
        "cohort_initial_count": int(initial_tags.size),
        "cohort_final_survival_fraction": cohort_survival[-1],
        "maximum_matched_momentum_gain_percent": maximum_matched_gain,
        "maximum_p99_over_injection": summary["maximum_p99_over_injection"],
        "final_p99_over_injection": summary["snapshots"][-1]["p99_over_injection"],
    }
    return paths, metrics, [summary_path, *mhd_paths, *particle_paths]


def _captions(bell: dict[str, Any], shock: dict[str, Any]) -> str:
    return f"""# Figure Captions

**Figure 1. Linear growth and nonlinear evolution of the Bell instability.**
(a) Amplitude of the right-handed magnetic eigenmode in the corrected-current
one-dimensional run (solid) compared with the linear dispersion relation
(dashed). The measured and predicted normalized growth rates are
{bell['linear_measured_growth']:.3f} and {bell['linear_expected_growth']:.3f}.
(b) Corresponding phase evolution; the measured and predicted frequencies are
{bell['linear_measured_phase']:.3f} and {bell['linear_expected_phase']:.3f}.
(c) Transverse and total rms magnetic-field amplitudes in the compact
two-dimensional nonlinear run. (d) The signed transverse component $B_y/B_0$
at the time of maximum rms amplification, when
$B_{{\\perp,\\mathrm{{rms}}}}/B_0={bell['nonlinear_peak_bperp_rms_over_B0']:.2f}$.
The simulations use the volume-aware deposited-current normalization.

**Figure 2. Structure and injected-particle transport in a compact coupled
non-relativistic shock.** (a) Density at $\\Omega_0 t=20$; the dashed line
marks the strongest transverse-averaged density gradient. (b) Weighted
particle phase space at the same time. (c) Transverse-averaged density,
shock-frame speed, and pressure profiles centered on the measured shock front.
(d) Momentum change and surviving fraction of the exact particle-tag cohort
present at the first nonempty output. The run contains
{shock['final_particle_count']:,} particles at the final output; the largest
matched-particle momentum increase is only
{shock['maximum_matched_momentum_gain_percent']:.3f} percent. This compact
uniform-grid calculation demonstrates shock formation, coupled injection,
feedback, and particle transport, but does not resolve diffusive shock
acceleration.
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--linear-root", type=Path, required=True)
    parser.add_argument("--nonlinear-root", type=Path, required=True)
    parser.add_argument("--shock-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=400)
    args = parser.parse_args()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    _style()

    bell_paths, bell_metrics, bell_inputs = make_bell_figure(
        args.linear_root.resolve(strict=True),
        args.nonlinear_root.resolve(strict=True),
        output_root,
        args.dpi,
    )
    shock_paths, shock_metrics, shock_inputs = make_shock_figure(
        args.shock_root.resolve(strict=True), output_root, args.dpi
    )
    caption_path = output_root / "figure_captions.md"
    caption_path.write_text(_captions(bell_metrics, shock_metrics), encoding="utf-8")

    repo_root = Path(__file__).resolve().parents[2]
    source_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()
    all_inputs = sorted({path.resolve() for path in [*bell_inputs, *shock_inputs]})
    all_outputs = [*bell_paths, *shock_paths, caption_path]
    manifest = {
        "record_type": "bell_shock_publication_figures_v1",
        "source_commit": source_commit,
        "generator": _artifact(Path(__file__)),
        "scope": (
            "publication_style_figures_from_compact_engineering_runs;"
            "not_convergence_or_publication_qualification"
        ),
        "bell_metrics": bell_metrics,
        "shock_metrics": shock_metrics,
        "inputs": [_artifact(path) for path in all_inputs],
        "outputs": [_artifact(path) for path in all_outputs],
    }
    manifest_path = output_root / "figure_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output_root": str(output_root), "figures": [str(p) for p in all_outputs]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
