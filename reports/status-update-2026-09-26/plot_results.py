"""Regenerate report figures from the saved, small analysis tables."""
import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FIGURES = ROOT / "figures/status_update"
FIGURES.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 180})

comparison = json.loads((DATA / "comparison.json").read_text())
fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
styles = [("current_limited", "Current limited", "C0"),
          ("constant", "Uniform resistivity", "C1"), ("ideal", "Ideal", "0.45")]
for model, label, color in styles:
    rows = sorted((r for r in comparison["flux_budgets"]
                   if r["model"] == model and r["cfl"] == 0.4),
                  key=lambda r: r["cells_per_di"])
    x = [r["cells_per_di"] for r in rows]
    axes[0].loglog(x, [abs(r["mean_flux_rate"]) for r in rows], "o-", color=color,
                   label=label)
    if model != "ideal":
        axes[0].loglog(x, [abs(r["mean_physical_EMF_difference"]) for r in rows],
                       ":", color=color)
    axes[1].loglog(x, [abs(r["mean_rate_discrepancy"]) for r in rows], "o-",
                   color=color, label=label)
for ax in axes:
    ax.set(xlabel=r"Upstream cells per $d_i$", xticks=[2, 4, 8])
    ax.set_xticklabels(["2", "4", "8"])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.grid(alpha=0.2)
axes[0].set(ylabel="Absolute interval-mean normalized rate",
            title=r"Startup interval: $0\leq t\leq0.02$")
axes[0].text(0.03, 0.15, "Solid: magnetic flux change\nDotted: physical X − reference EMF",
             transform=axes[0].transAxes, fontsize=9)
axes[0].legend(fontsize=9, loc="upper right")
axes[1].set(ylabel="Absolute difference\nflux rate − physical EMF",
            title="Ideal control identifies truncation error")
fig.savefig(FIGURES / "startup_convergence.png")
plt.close(fig)

rows = []
for name in ("cl-cpd2-pilot", "cl-cpd2-t1", "cl-cpd2-t5"):
    with (DATA / f"{name}.csv").open() as stream:
        rows.extend(csv.DictReader(stream))
rows.sort(key=lambda r: float(r["time"]))
t = np.array([float(r["time"]) for r in rows])
fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
axes[0].plot(t, [float(r["q_max_edge"]) for r in rows], label="Maximum edge q")
axes[0].axhline(1-np.sqrt(1e-6/0.005), color="0.4", ls="--", label=r"Transition $q_*$")
axes[0].set(xlabel=r"$t/(L/v_{A0})$", ylabel="q", title="Coarse onset pilot")
axes[0].legend(fontsize=9)
# Derivatives across all segments avoid treating each restart as a fresh experiment.
psi = np.array([float(r["psi_ref_minus_X"]) for r in rows])
rate = np.gradient(psi, t)
rate[[r["fixed_X_verified"] != "True" for r in rows]] = np.nan
axes[1].plot(t, rate, label="Fixed-reference flux derivative")
axes[1].plot(t, [float(r["physical_EMF_difference"]) for r in rows],
             label="Physical X − reference EMF")
axes[1].set(xlabel=r"$t/(L/v_{A0})$", ylabel="Signed normalized rate",
            title="Coarse grid: signed flux rates")
axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
axes[1].legend(fontsize=9)
for ax in axes:
    ax.grid(alpha=0.2)
fig.savefig(FIGURES / "coarse_evolution.png")
plt.close(fig)
