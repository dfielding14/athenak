#!/usr/bin/env python3
"""Generate unqualified shock-rich engineering figures from diagnostics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np

from artifact_lineage import write_companion_manifest

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError as exc:
    raise SystemExit("matplotlib is required: python3 -m pip install matplotlib") from exc


_PROXY_WATERMARK = "ENGINEERING PROXY - NOT SUN & BAI (2023) REPRODUCTION"


def _parse_log_metrics(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")

    def _match_float(pattern: str, default: float = float("nan")) -> float:
        m = re.search(pattern, text)
        return float(m.group(1)) if m else default

    def _match_int(pattern: str, default: int = 0) -> int:
        m = re.search(pattern, text)
        return int(m.group(1)) if m else default

    return {
        "meshblocks": _match_int(r"Current number of MeshBlocks\s*=\s*(\d+)"),
        "created": _match_int(r"(\d+) MeshBlocks created,\s*\d+ deleted by AMR"),
        "deleted": _match_int(r"\d+ MeshBlocks created,\s*(\d+) deleted by AMR"),
        "communicated": _match_int(r"(\d+) communicated for load balancing"),
        "efficiency": _match_float(r"load balancing efficiency\s*=\s*([0-9eE+.-]+)"),
        "zone_rate": _match_float(r"zone-cycles/cpu_second\s*=\s*([0-9eE+.-]+)"),
        "particle_rate": _match_float(r"particle-updates/cpu_second\s*=\s*([0-9eE+.-]+)"),
    }


def _load_scan_runtime(scan_summary: Path) -> tuple[np.ndarray, np.ndarray]:
    data = json.loads(scan_summary.read_text(encoding="utf-8"))
    xs = []
    ys = []
    for case in data.get("cases", []):
        if case.get("status") != "passed":
            continue
        xs.append(float(case.get("nx", 0.0)))
        ys.append(float(case.get("runtime_s", 0.0)))
    if not xs:
        return np.array([]), np.array([])

    order = np.argsort(xs)
    return np.asarray(xs)[order], np.asarray(ys)[order]


def _panel_a(ax, profiles: dict) -> None:
    x = profiles["x"]
    rho = profiles["rho_x"]
    bmag = profiles["bmag_x"]
    j2 = profiles["j2_x"]
    shock_x = float(profiles["shock_x"][0])

    rho_n = rho / max(np.max(np.abs(rho)), 1.0e-30)
    bmag_n = bmag / max(np.max(np.abs(bmag)), 1.0e-30)
    j2_n = j2 / max(np.max(np.abs(j2)), 1.0e-30)

    ax.plot(x, rho_n, label=r"$\rho$", lw=1.8)
    ax.plot(x, bmag_n, label=r"$|B|$", lw=1.4)
    ax.plot(x, j2_n, label=r"$J^2$", lw=1.2)
    ax.axvline(shock_x, color="k", linestyle="--", linewidth=1.0, label="shock")
    ax.set_xlabel("x")
    ax.set_ylabel("normalized profile")
    ax.set_title("(A) Shock Structure: density, magnetic, current")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", fontsize=8)


def _panel_b(ax, xp: dict) -> None:
    hist = xp["hist"]
    x_edges = xp["x_edges"]
    p_edges = xp["p_edges"]

    # Avoid log(0) in display while preserving empties as dark bins.
    display = np.log10(hist.T + 1.0)
    mesh = ax.pcolormesh(x_edges, p_edges, display, shading="auto", cmap="magma")
    ax.set_yscale("log")
    ax.set_xlabel("x")
    ax.set_ylabel(r"$p_{proxy}=|v|$")
    ax.set_title("(B) CR Phase-Space Proxy: x-p histogram")
    cbar = plt.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$\log_{10}(N+1)$")


def _panel_c(ax, spectra: dict) -> None:
    p_edges = spectra["p_edges"]
    p_centers = np.sqrt(p_edges[:-1] * p_edges[1:])

    g = spectra["dndp_global"]
    u = spectra["dndp_upstream"]
    d = spectra["dndp_downstream"]

    ax.loglog(p_centers, g + 1.0, label="global", lw=1.8)
    ax.loglog(p_centers, u + 1.0, label="upstream", lw=1.4)
    ax.loglog(p_centers, d + 1.0, label="downstream", lw=1.4)
    ax.set_xlabel(r"$p_{proxy}$")
    ax.set_ylabel(r"$dN/dp$ (count proxy)")
    ax.set_title("(C) CR Spectrum Proxy")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)


def _panel_d(ax, local_metrics: dict, hpc_metrics: dict,
             scan_nx: np.ndarray, scan_runtime: np.ndarray) -> None:
    if scan_nx.size > 0:
        ax.plot(scan_nx, scan_runtime, marker="o", lw=1.5, color="tab:blue")
        ax.set_xlabel("resolution nx")
        ax.set_ylabel("runtime (s)")
    else:
        ax.set_xlabel("resolution nx")
        ax.set_ylabel("runtime (s)")
        ax.text(0.5, 0.5, "No scan runtime data", transform=ax.transAxes,
                ha="center", va="center")

    ax.set_title("(D) AMR/Load-Balance Context")
    ax.grid(alpha=0.25)

    txt = (
        "local: MB={l_mb}, AMR +{l_cr}/{l_del}, comm={l_comm}, "
        "eff={l_eff:.3f}\n"
        "hpc-plumb: MB={h_mb}, AMR +{h_cr}/{h_del}, comm={h_comm}, "
        "eff={h_eff:.3f}\n"
        "local zone/s={l_zone:.2e}, prtcl/s={l_part:.2e}\n"
        "hpc zone/s={h_zone:.2e}, prtcl/s={h_part:.2e}"
    ).format(
        l_mb=local_metrics["meshblocks"],
        l_cr=local_metrics["created"],
        l_del=local_metrics["deleted"],
        l_comm=local_metrics["communicated"],
        l_eff=local_metrics["efficiency"],
        h_mb=hpc_metrics["meshblocks"],
        h_cr=hpc_metrics["created"],
        h_del=hpc_metrics["deleted"],
        h_comm=hpc_metrics["communicated"],
        h_eff=hpc_metrics["efficiency"],
        l_zone=local_metrics["zone_rate"],
        l_part=local_metrics["particle_rate"],
        h_zone=hpc_metrics["zone_rate"],
        h_part=hpc_metrics["particle_rate"],
    )
    ax.text(
        0.02,
        0.98,
        txt,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.2", "fc": "white", "alpha": 0.8},
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot shock-rich engineering-stress multi-panel figure."
    )
    parser.add_argument("--diagnostics-dir", required=True)
    parser.add_argument("--step2-dir", required=True)
    parser.add_argument("--scan-summary", default="")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--figure-name", default="shock_storyline_multipanel.png")
    parser.add_argument("--dpi", type=int, default=220)
    args = parser.parse_args()

    diagnostics_dir = Path(args.diagnostics_dir).resolve()
    step2_dir = Path(args.step2_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    profiles = np.load(diagnostics_dir / "shock_profiles_latest.npz")
    xp = np.load(diagnostics_dir / "phase_space_xp.npz")
    spectra = np.load(diagnostics_dir / "spectra.npz")

    local_log = step2_dir / "local_np8.log"
    hpc_log = step2_dir / "hpc_plumb_np8.log"
    local_metrics = _parse_log_metrics(local_log)
    hpc_metrics = _parse_log_metrics(hpc_log)

    scan_nx = np.array([])
    scan_runtime = np.array([])
    if args.scan_summary:
        summary_path = Path(args.scan_summary).resolve()
        if summary_path.exists():
            scan_nx, scan_runtime = _load_scan_runtime(summary_path)

    fig, axs = plt.subplots(2, 2, figsize=(13.5, 9.0))
    _panel_a(axs[0, 0], profiles)
    _panel_b(axs[0, 1], xp)
    _panel_c(axs[1, 0], spectra)
    _panel_d(axs[1, 1], local_metrics, hpc_metrics, scan_nx, scan_runtime)

    fig.tight_layout()
    fig.text(
        0.995, 0.005, _PROXY_WATERMARK, ha="right", va="bottom",
        fontsize=7, color="#8b0000",
    )
    fig_path = output_dir / args.figure_name
    fig.savefig(fig_path, dpi=args.dpi)
    plt.close(fig)

    summary = {
        "figure": str(fig_path),
        "evidence_class": "engineering_proxy",
        "not_sun_bai_reproduction": True,
        "qualification_status": "unqualified_engineering_quicklook",
        "physical_model": "orszag_tang_shock_rich_engineering_stress",
        "visible_watermark": _PROXY_WATERMARK,
        "diagnostics_dir": str(diagnostics_dir),
        "step2_dir": str(step2_dir),
        "scan_summary": str(args.scan_summary),
        "local_metrics": local_metrics,
        "hpc_metrics": hpc_metrics,
        "scan_points": int(scan_nx.size),
    }
    figure_summary = output_dir / "figure_summary.json"
    figure_summary.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    inputs = [
        diagnostics_dir / "shock_profiles_latest.npz",
        diagnostics_dir / "phase_space_xp.npz",
        diagnostics_dir / "spectra.npz",
        local_log,
        hpc_log,
    ]
    if args.scan_summary:
        inputs.append(Path(args.scan_summary).resolve())
    write_companion_manifest(
        output_dir / "figure_artifact_manifest.json",
        generator="plot_pic_shock_storyline.py",
        physical_model="orszag_tang_shock_rich_engineering_stress",
        inputs=inputs,
        outputs=(fig_path, figure_summary),
    )

    print("figure:", fig_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
