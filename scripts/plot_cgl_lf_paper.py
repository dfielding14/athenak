#!/usr/bin/env python3
"""Render reduced MKS24-oriented plots from CGL-LF paper diagnostics JSON."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def save(fig: plt.Figure, path: Path, generated: list[str]) -> None:
    """Write one figure and retain its product name."""

    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    generated.append(path.name)


def curve(record: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    """Return histogram centers and density."""

    edges = np.asarray(record["edges"], dtype=float)
    return 0.5 * (edges[1:] + edges[:-1]), np.asarray(record["density"], dtype=float)


def spectral_curve(record: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    """Return positive spectrum coordinates suitable for log plotting."""

    k = np.asarray(record["k"], dtype=float)
    power = np.asarray(record["power_per_dk"], dtype=float)
    selected = (k > 0.0) & (power > 0.0) & np.isfinite(power)
    return k[selected], power[selected]


def case_ensembles(data: dict[str, object]) -> dict[str, dict[str, object]]:
    """Select cases for which snapshot-derived analysis exists."""

    selected: dict[str, dict[str, object]] = {}
    direct = data.get("snapshot_ensemble", {})
    if isinstance(direct, dict) and direct.get("snapshot_count", 0) > 0:
        selected["direct"] = direct
    cases = data.get("cases", {})
    if not isinstance(cases, dict):
        return selected
    for name, case in cases.items():
        ensemble = case.get("snapshot_ensemble", {})
        if isinstance(ensemble, dict) and ensemble.get("snapshot_count", 0) > 0:
            selected[str(name)] = ensemble
    return selected


def plot_history(data: dict[str, object], figure_dir: Path,
                 generated: list[str]) -> None:
    """Plot retained global history summaries when available."""

    histories = data.get("histories", [])
    if not isinstance(histories, list) or not histories:
        return
    labels = [Path(str(item["path"])).stem.replace(".user", "") for item in histories]
    cb2 = np.asarray([item.get("c_b2_final", np.nan) for item in histories], dtype=float)
    anis = np.asarray([item.get("abs_dp_mean_final", np.nan)
                       for item in histories], dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.2))
    axes[0].bar(labels, cb2, color="#1f77b4")
    axes[0].set_ylabel(r"$C_{B^2}$ at final output")
    axes[1].bar(labels, anis, color="#d62728")
    axes[1].set_ylabel(r"$\langle |\Delta p| \rangle$ at final output")
    for axis in axes:
        axis.tick_params(axis="x", labelrotation=30)
        axis.grid(True, axis="y", alpha=0.25)
    fig.suptitle("Archived CGL-LF reduced histories")
    save(fig, figure_dir / "paper_history_summary.pdf", generated)


def plot_pdfs(ensembles: dict[str, dict[str, object]], figure_dir: Path,
              generated: list[str]) -> None:
    """Plot anisotropy and rate-of-strain PDFs."""

    if not ensembles:
        return
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2))
    for name, ensemble in ensembles.items():
        pdfs = ensemble["pdf"]
        for axis, field, label in (
            (axes[0], "beta_delta", r"$\beta\Delta$"),
            (axes[1], "bb_grad_velocity", r"$\hat b\hat b:\nabla u$"),
        ):
            if field in pdfs:
                x, density = curve(pdfs[field])
                axis.plot(x, density, label=name)
                axis.set_xlabel(label)
    axes[0].set_ylabel("PDF")
    for axis in axes:
        axis.grid(True, alpha=0.25)
    axes[0].legend(fontsize=7)
    fig.suptitle("Pressure anisotropy and strain production")
    save(fig, figure_dir / "paper_pdfs.pdf", generated)


def plot_spatial_spectra(ensembles: dict[str, dict[str, object]], figure_dir: Path,
                         generated: list[str]) -> None:
    """Plot projected pressure and velocity gradient spectra."""

    if not ensembles:
        return
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.25))
    for name, ensemble in ensembles.items():
        spectra = ensemble["spectra"]
        for key, linestyle in (
            ("grad_parallel_delta_p", "-"),
            ("grad_perp_delta_p", "--"),
        ):
            k, power = spectral_curve(spectra[key])
            axes[0].loglog(k, power, linestyle=linestyle, label=f"{name} {key}")
        for key, linestyle in (
            ("grad_parallel_velocity_parallel", "-"),
            ("grad_perp_velocity_parallel", "--"),
        ):
            k, power = spectral_curve(spectra[key])
            axes[1].loglog(k, power, linestyle=linestyle, label=f"{name} {key}")
    axes[0].set_title(r"$\Delta p$ local gradients")
    axes[1].set_title(r"$u_\parallel$ local gradients")
    for axis in axes:
        axis.set_xlabel(r"$k_\perp$")
        axis.set_ylabel("spectral power")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=6)
    save(fig, figure_dir / "paper_projected_spectra.pdf", generated)


def plot_transfer_and_alignment(ensembles: dict[str, dict[str, object]],
                                figure_dir: Path, generated: list[str]) -> None:
    """Plot pressure-stress transfer and selected alignment peaks."""

    if not ensembles:
        return
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.25))
    for name, ensemble in ensembles.items():
        transfer = ensemble["pressure_transfer"]
        k = np.asarray(transfer["k_perp"], dtype=float)
        values = np.asarray(transfer["transfer"], dtype=float)
        selected = k > 0.0
        axes[0].plot(k[selected], values[selected], label=name)
        shells: list[float] = []
        peaks: list[float] = []
        for shell, distribution in sorted(
            ensemble["alignment"].items(), key=lambda item: int(item[0])
        ):
            centers, density = curve(distribution)
            shells.append(float(shell))
            peaks.append(float(centers[np.argmax(density)]))
        if shells:
            axes[1].plot(shells, peaks, marker="o", label=name)
    axes[0].axhline(0.0, color="0.3", lw=0.8)
    axes[0].set_xlabel(r"$k_\perp$")
    axes[0].set_ylabel(r"$\mathcal{T}_{\Delta p}$")
    axes[0].set_title("Anisotropic-pressure transfer")
    axes[1].set_xlabel(r"$k_\perp$ shell")
    axes[1].set_ylabel(r"peak $|\cos\theta|$")
    axes[1].set_ylim(0.0, 1.0)
    axes[1].set_title("Strain-field alignment")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize=7)
    save(fig, figure_dir / "paper_transfer_alignment.pdf", generated)


def plot_heat_flux_proxy(ensembles: dict[str, dict[str, object]], figure_dir: Path,
                         generated: list[str]) -> None:
    """Plot closure-reconstructed heat-flux temperature-smoothing proxies."""

    retained = {
        name: ensemble["heat_flux_transport_proxy"]
        for name, ensemble in ensembles.items()
        if ensemble.get("heat_flux_transport_proxy", {}).get("available", False)
    }
    if not retained:
        return
    labels = list(retained)
    locations = np.arange(len(labels), dtype=float)
    width = 0.25
    parallel = np.asarray(
        [retained[name]["regularized_parallel_power_mean"] for name in labels]
    )
    perpendicular = np.asarray(
        [retained[name]["regularized_perpendicular_power_mean"] for name in labels]
    )
    unlimited = np.asarray(
        [retained[name]["unlimited_total_power_mean"] for name in labels]
    )
    fig, axis = plt.subplots(figsize=(max(5.0, 1.1 * len(labels) + 2.5), 3.25))
    axis.bar(locations - width / 2.0, parallel, width, label=r"regularized $q_\parallel$")
    axis.bar(
        locations + width / 2.0,
        perpendicular,
        width,
        label=r"regularized $q_\perp$",
    )
    axis.plot(locations, unlimited, color="0.2", marker="o", label="unlimited total")
    axis.axhline(0.0, color="0.3", lw=0.8)
    axis.set_xticks(locations, labels, rotation=30, ha="right")
    axis.set_ylabel("integrated smoothing proxy")
    axis.set_title("Snapshot-reconstructed LF heat-flux proxy")
    axis.grid(True, axis="y", alpha=0.25)
    axis.legend(fontsize=7)
    save(fig, figure_dir / "paper_heat_flux_proxy.pdf", generated)


def plot_pressure_work(ensembles: dict[str, dict[str, object]], figure_dir: Path,
                       generated: list[str]) -> None:
    """Plot snapshot-reconstructed CGL pressure mechanical-power terms."""

    retained = {
        name: ensemble["pressure_work_decomposition"]
        for name, ensemble in ensembles.items()
        if ensemble.get("pressure_work_decomposition", {}).get("available", False)
    }
    if not retained:
        return
    labels = list(retained)
    locations = np.arange(len(labels), dtype=float)
    width = 0.25
    perpendicular = np.asarray([
        retained[name]["isotropic_perpendicular_pressure_power_mean"]
        for name in labels
    ])
    anisotropic = np.asarray([
        retained[name]["anisotropic_stress_power_mean"] for name in labels
    ])
    total = np.asarray([
        retained[name]["total_cgl_pressure_power_mean"] for name in labels
    ])
    fig, axis = plt.subplots(figsize=(max(5.0, 1.1 * len(labels) + 2.5), 3.25))
    axis.bar(
        locations - width / 2.0,
        perpendicular,
        width,
        label=r"$p_\perp\nabla\cdot u$",
    )
    axis.bar(
        locations + width / 2.0,
        anisotropic,
        width,
        label=r"$-\Delta p\,\hat b\hat b:\nabla u$",
    )
    axis.plot(locations, total, color="0.2", marker="o", label="total CGL pressure")
    axis.axhline(0.0, color="0.3", lw=0.8)
    axis.set_xticks(locations, labels, rotation=30, ha="right")
    axis.set_ylabel("integrated pressure power")
    axis.set_title("Snapshot-reconstructed CGL pressure work")
    axis.grid(True, axis="y", alpha=0.25)
    axis.legend(fontsize=7)
    save(fig, figure_dir / "paper_pressure_work.pdf", generated)


def main() -> int:
    """Command-line entry point."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics", required=True, type=Path)
    parser.add_argument("--figure-dir", required=True, type=Path)
    args = parser.parse_args()
    data = json.loads(args.diagnostics.read_text(encoding="utf-8"))
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    generated: list[str] = []
    ensembles = case_ensembles(data)
    plot_history(data, args.figure_dir, generated)
    plot_pdfs(ensembles, args.figure_dir, generated)
    plot_spatial_spectra(ensembles, args.figure_dir, generated)
    plot_transfer_and_alignment(ensembles, args.figure_dir, generated)
    plot_heat_flux_proxy(ensembles, args.figure_dir, generated)
    plot_pressure_work(ensembles, args.figure_dir, generated)
    (args.figure_dir / "paper_figures.json").write_text(
        json.dumps({"figures": generated}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(generated)} CGL-LF paper diagnostic figures to {args.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
