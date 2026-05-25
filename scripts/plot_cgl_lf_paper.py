#!/usr/bin/env python3
"""Render reduced MKS24-oriented plots from CGL-LF paper diagnostics JSON."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

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


def plot_pressure_density_joint(ensembles: dict[str, dict[str, object]],
                                figure_dir: Path, generated: list[str]) -> None:
    """Plot Figure 2(a)-style pressure-density joint distributions."""

    for name, ensemble in ensembles.items():
        product = ensemble.get("pressure_density_joint", {})
        if not isinstance(product, dict):
            continue
        if not all(key in product for key in ("parallel", "perpendicular")):
            continue
        fig, axes = plt.subplots(2, 1, figsize=(5.0, 6.0), sharex=True)
        for axis, key, ylabel in (
            (axes[0], "parallel", r"$\delta p_\parallel$"),
            (axes[1], "perpendicular", r"$\delta p_\perp$"),
        ):
            record = product[key]
            x_edges = np.asarray(record["x_edges"], dtype=float)
            y_edges = np.asarray(record["y_edges"], dtype=float)
            density = np.asarray(record["density"], dtype=float)
            image = axis.pcolormesh(
                x_edges, y_edges, density.T, shading="auto", cmap="viridis"
            )
            limits = np.asarray([x_edges[0], x_edges[-1]], dtype=float)
            axis.plot(limits, (5.0 / 3.0) * limits, "k:", lw=1.0)
            axis.set_ylabel(ylabel)
            fig.colorbar(image, ax=axis, label="joint PDF")
        axes[1].set_xlabel(r"$\langle p\rangle\,\delta\rho/\langle\rho\rangle$")
        fig.suptitle(f"Pressure-density joint PDF: {name}")
        suffix = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
        save(fig, figure_dir / f"paper_pressure_density_{suffix}.pdf", generated)


def plot_compressive_pressure_spectra(
    ensembles: dict[str, dict[str, object]], figure_dir: Path, generated: list[str]
) -> None:
    """Plot Figure 3/4/6-oriented compressive and pressure spectra."""

    if not ensembles:
        return
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.25))
    for name, ensemble in ensembles.items():
        spectra = ensemble["spectra"]
        for axis, key in (
            (axes[0], "compressive_velocity"),
            (axes[1], "density_fluctuation"),
        ):
            if key in spectra:
                k, power = spectral_curve(spectra[key])
                axis.loglog(k, power, label=name)
        for key, label, linestyle in (
            ("p_parallel", r"$p_\parallel$", "-"),
            ("p_perp", r"$p_\perp$", "--"),
            ("magnetic_pressure", r"$B^2/2$", ":"),
        ):
            if key in spectra:
                k, power = spectral_curve(spectra[key])
                axes[2].loglog(
                    k, power, linestyle=linestyle, label=f"{name} {label}"
                )
    axes[0].set_title(r"Compressive flow $\hat{k}\cdot u_k$")
    axes[1].set_title(r"Density fluctuation $\rho/\langle\rho\rangle-1$")
    axes[2].set_title("Thermal and magnetic pressure")
    for axis in axes:
        axis.set_xlabel(r"$k_\perp$")
        axis.set_ylabel("spectral power")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=6)
    save(fig, figure_dir / "paper_compressive_pressure_spectra.pdf", generated)


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
    normalized = all(
        ensemble.get("pressure_transfer", {}).get("normalization_available", False)
        for ensemble in ensembles.values()
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.25))
    for name, ensemble in ensembles.items():
        transfer = ensemble["pressure_transfer"]
        k = np.asarray(transfer["k_perp"], dtype=float)
        key = "transfer_normalized_by_total" if normalized else "transfer"
        values = np.asarray(transfer[key], dtype=float)
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
    axes[0].set_ylabel(
        r"$\mathcal{T}_{\Delta p}/\mathcal{T}_{\rm total}$"
        if normalized else r"$\mathcal{T}_{\Delta p}$"
    )
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


def plot_eddy_anisotropy(ensembles: dict[str, dict[str, object]], figure_dir: Path,
                         generated: list[str]) -> None:
    """Plot local-field-conditioned parallel/perpendicular eddy lengths."""

    retained = {
        name: ensemble["eddy_anisotropy"]
        for name, ensemble in ensembles.items()
        if ensemble.get("eddy_anisotropy", {}).get("available", False)
    }
    if not retained:
        return
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.25))
    for name, product in retained.items():
        for axis, key, title in (
            (axes[0], "velocity_perp", r"$u_\perp$"),
            (axes[1], "magnetic_perp", r"$B_\perp$"),
        ):
            curve_product = product.get(key, {})
            if not curve_product.get("available", False):
                continue
            axis.loglog(
                curve_product["ell_perp_over_lperp"],
                curve_product["ell_parallel_over_lperp"],
                label=name,
            )
            axis.set_title(title)
    for axis in axes:
        axis.set_xlabel(r"$\ell_\perp/L_\perp$")
        axis.set_ylabel(r"$\ell_\parallel/L_\perp$")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=7)
    fig.suptitle("Local-field-conditioned eddy anisotropy")
    save(fig, figure_dir / "paper_eddy_anisotropy.pdf", generated)


def plot_reference_comparisons(data: dict[str, object], figure_dir: Path,
                               generated: list[str]) -> None:
    """Plot optional provenance-qualified reference and simulated curves."""

    reference = data.get("reference_curve_comparisons", {})
    if not isinstance(reference, dict) or not reference.get("available", False):
        return
    comparisons = reference.get("comparisons", {})
    if not isinstance(comparisons, dict) or not comparisons:
        return
    names = list(comparisons)
    fig, axes = plt.subplots(
        len(names), 1, figsize=(6.2, max(3.0, 2.8 * len(names))), squeeze=False
    )
    for axis, name in zip(axes[:, 0], names):
        record = comparisons[name]
        x = np.asarray(record["x"], dtype=float)
        y = np.asarray(record["reference_y"], dtype=float)
        uncertainty = np.asarray(record["reference_y_uncertainty"], dtype=float)
        simulated = np.asarray(record["simulated_y"], dtype=float)
        axis.errorbar(x, y, yerr=uncertainty, fmt="o", ms=3, label="reference")
        axis.plot(x, simulated, label="analysis")
        if record.get("interpolation") == "loglog":
            axis.set_xscale("log")
            axis.set_yscale("log")
        axis.set_title(f"{name}: {record['product']}", fontsize=9)
        axis.set_xlabel("reference coordinate")
        axis.set_ylabel("diagnostic value")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=7)
    fig.suptitle("Reference-curve comparisons")
    save(fig, figure_dir / "paper_reference_comparisons.pdf", generated)


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
    plot_pressure_density_joint(ensembles, args.figure_dir, generated)
    plot_compressive_pressure_spectra(ensembles, args.figure_dir, generated)
    plot_spatial_spectra(ensembles, args.figure_dir, generated)
    plot_transfer_and_alignment(ensembles, args.figure_dir, generated)
    plot_heat_flux_proxy(ensembles, args.figure_dir, generated)
    plot_pressure_work(ensembles, args.figure_dir, generated)
    plot_eddy_anisotropy(ensembles, args.figure_dir, generated)
    plot_reference_comparisons(data, args.figure_dir, generated)
    (args.figure_dir / "paper_figures.json").write_text(
        json.dumps({"figures": generated}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(generated)} CGL-LF paper diagnostic figures to {args.figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
