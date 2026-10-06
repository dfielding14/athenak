"""Read the preserved B4 runs; derive exact-advection errors and comparison figures."""
from pathlib import Path
import argparse
import json
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator


def final_state(path, family):
    if family == "squire":
        data = np.genfromtxt(sorted(path.glob("state_*_cycle_0.csv"))[-1],
                             delimiter=",", names=True)
        return dict(x=data["x"], rho=data["rho"], vx=data["mx"]/data["rho"],
                    ppar=data["ppar"], pperp=data["pperp"], B=data["By"],
                    time=float(data["time"][0]))
    file = sorted((path/"tab").glob("*mhd_w_bcc*.tab"))[-1]
    data = np.loadtxt(file)
    time = float(re.search(r"time=([\deE+.-]+)", file.read_text().splitlines()[0])[1])
    return dict(x=data[:, 2], rho=data[:, 3], vx=data[:, 4], ppar=data[:, 7],
                pperp=data[:, 8], B=data[:, 10], time=time)


def main(root, output):
    output.mkdir(parents=True, exist_ok=True)
    current = json.loads((root/"current/summary.json").read_text())
    budgets = []
    for r in current:
        if r["width"] == 0:
            residuals = [abs(c["sum_E"]/128 - 52.125 - .25*abs(r["velocity"])*c["time"])
                         for c in r["cycles"]]
            budgets.append(dict(case=Path(r["path"]).name,
                                max_abs_mean_energy_residual=max(residuals)))
    assert max(r["max_abs_mean_energy_residual"] for r in budgets) < 1e-12
    (output.parent/"energy-budget.json").write_text(json.dumps(budgets, indent=2)+"\n")
    matrix = []
    for r in current:
        if r["width"] == 0 and r["velocity"] == 10 and r["bfloor"] == "1e-10":
            matrix.append({"face": r["face"], "wall": r["wall"], "stop": r["stop"],
                           "first_failure": r["first_failure_ratio_outside_half_to_two"],
                           "final": r["final"], "minimum_pressure_over_run": r["min_pressure"]})
    smooth = []
    profiles = {}
    for family in ("current", "majeski", "squire"):
        for path in sorted((root/family/"runs").glob("smooth*")):
            s = final_state(path, family)
            assert abs(s["time"] - .012) < 1e-12, (path, s["time"])
            assert all(np.isfinite(s[k]).all() for k in ("rho", "vx", "ppar", "pperp", "B"))
            assert (s["ppar"] > 0).all() and (s["pperp"] > 0).all()
            exact_B = 1e-12+(1-1e-12)*.5*(1+np.tanh((s["x"]-10*s["time"]-.5)/.04))
            exact_p = 1.5-.5*exact_B**2
            ratio = s["pperp"]/s["ppar"]
            log_ratio = np.log(ratio)
            wall = "_fh1_" in path.name if family == "current" else None
            row = dict(family=family, wall=wall, case=path.name, nx=len(ratio), time=s["time"],
                       min_ratio=float(ratio.min()), max_ratio=float(ratio.max()),
                       l1_logratio=float(abs(log_ratio).mean()), linf_logratio=float(abs(log_ratio).max()),
                       l1_ratio_error=float(abs(ratio-1).mean()),
                       l1_B=float(abs(s["B"]-exact_B).mean()), linf_B=float(abs(s["B"]-exact_B).max()),
                       l1_rho=float(abs(s["rho"]-1).mean()), linf_rho=float(abs(s["rho"]-1).max()),
                       linf_velocity=float(abs(s["vx"]-10).max()),
                       l1_ppar=float(abs(s["ppar"]-exact_p).mean()),
                       l1_pperp=float(abs(s["pperp"]-exact_p).mean()))
            smooth.append(row)
            profiles[(family, wall, len(ratio))] = (s["x"], ratio)
    assert len(smooth) == 12, len(smooth)
    for row in smooth:
        previous = next((r for r in smooth if r["family"] == row["family"] and
                         r["wall"] == row["wall"] and r["nx"]*2 == row["nx"]), None)
        if previous:
            row["l1_logratio_order"] = float(np.log2(previous["l1_logratio"]/row["l1_logratio"]))
    (output/"derived-results.json").write_text(json.dumps(dict(matrix=matrix, smooth=smooth), indent=2)+"\n")

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    for face, color in [("original", "#a05000"), ("wo1", "#1764a0")]:
        for wall, style in [(True, "-"), (False, "--")]:
            r = next(r for r in current if r["face"] == face and r["wall"] == wall and
                     r["velocity"] == 10 and r["bfloor"] == "1e-10" and r["stop"] == "cycles")
            axes[0].plot([0]+[c["cycle"] for c in r["cycles"]],
                         [1]+[c["max_ratio"] for c in r["cycles"]], color=color, ls=style,
                         label=f"{'Original' if face == 'original' else 'WO1'} face; wall {'on' if wall else 'off'}")
    axes[0].axhline(2, color=".4", lw=1, ls=":")
    axes[0].set(yscale="log", xlabel="Completed cycle", ylabel=r"Maximum $p_\perp/p_\parallel$",
                title="Sharp contact: all four variants fail", ylim=(.8, 1e13))
    axes[0].legend(fontsize=9, loc="lower right")
    for family, wall, label, marker, color, size in [
        ("squire", None, "Squire archive", "s", "#777777", 10),
        ("majeski", None, "Majeski reference", "o", "#a05000", 7),
        ("current", False, "Current; wall off", "x", "#1764a0", 7),
        ("current", True, "Current; wall on", "+", "#152942", 11),
    ]:
        rows = sorted([r for r in smooth if r["family"] == family and r["wall"] == wall], key=lambda r:r["nx"])
        axes[1].loglog([r["nx"] for r in rows], [r["l1_logratio"] for r in rows],
                       marker=marker, ms=size, mfc="none", color=color, lw=1, label=label)
    axes[1].loglog([128, 512], [.13225, .13225/16], ":", color=".5", label=r"$N^{-2}$")
    axes[1].set(xlabel="Cells N", ylabel=r"Mean $|\ln(p_\perp/p_\parallel)|$",
                title="Resolved contact: reference agreement", xticks=[128,256,512])
    axes[1].set_xticklabels(["128","256","512"])
    axes[1].xaxis.set_minor_locator(NullLocator())
    axes[1].legend(fontsize=9)
    fig.savefig(output/"comparison.png", dpi=180)
    fig.savefig(output/"comparison.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
    for n, color in [(128,"#a05000"),(256,"#287855"),(512,"#1764a0")]:
        x, ratio = profiles[("current", True, n)]
        ax.plot(x, ratio, color=color, label=f"N={n}")
    ax.axhline(1, color=".4", lw=1, ls=":", label="Exact advected solution")
    ax.set(xlabel="x", ylabel=r"$p_\perp/p_\parallel$", title="Current code: smooth contact at t=0.012")
    ax.legend(fontsize=9)
    fig.savefig(output/"smooth-profiles.png", dpi=180)
    fig.savefig(output/"smooth-profiles.pdf")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs-root", type=Path, default=Path("/tmp/cgl-b4-compare-20261006"))
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent/"figures")
    args = parser.parse_args()
    main(args.runs_root, args.output)
