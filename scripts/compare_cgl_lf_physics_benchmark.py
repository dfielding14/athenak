#!/usr/bin/env python3
"""Compare retained, matched active/passive CGL-LF turbulence runs.

Example: compare_cgl_lf_physics_benchmark.py ACTIVE_RUN PASSIVE_RUN \
    --time-start 6 --time-end 18 --block-duration 2 --output-dir COMPARISON

Streams both primitive snapshot intervals to determine common exact PDF edges,
then reuses the single-run analyzer for all measurements and supplements. Never
launches a simulation or treats dependent time blocks as independent samples.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_cgl_lf_physics_benchmark as single  # noqa: E402


MODES = ("active", "passive")
MODE_KEYS = {("mhd", "passive"), ("problem", "passive_delta"),
             ("mhd", "passive_restart_encoding")}
ADMIN_KEYS = {("job", "basename"), ("time", "tlim"), ("time", "nlim"),
              ("time", "ndiag"), ("time", "ncycle_out"),
              ("time", "restart_time"),
              ("time", "walltime_limit"), ("time", "wall_time_limit")}
PDF_FIELDS = ("eint", "p_perp", "bcc1", "bcc2", "bcc3")


def canonical(value):
    """Compare numeric spellings exactly; no tolerance hides changed physics."""
    text = str(value).strip()
    if text.lower() in ("true", "false"):
        return ("boolean", text.lower() == "true")
    try:
        number = Decimal(text)
        if number.is_finite():
            return ("number", number)
    except InvalidOperation:
        pass
    return ("text", text)


def administrative(section, key):
    return section.startswith("output") or (section, key) in ADMIN_KEYS


def compare_parameters(active, passive):
    """Fail closed on any unlisted physics, integrator, mesh or forcing change."""
    for mode, header in zip(MODES, (active, passive)):
        actual = header.get("mhd", {}).get("passive", "").lower()
        if actual != str(mode == "passive").lower():
            raise ValueError(f"{mode} run must have mhd/passive={mode == 'passive'}")
        alias = header.get("problem", {}).get("passive_delta")
        if alias is not None and alias.lower() != actual:
            raise ValueError(f"{mode}: problem/passive_delta disagrees with mhd/passive")
        if header.get("mhd", {}).get("eos") != "cgl":
            raise ValueError(f"{mode} run must use mhd/eos=cgl")
        encoding = header.get("mhd", {}).get("passive_restart_encoding")
        if encoding is not None and encoding != "1":
            raise ValueError(f"{mode}: unsupported passive_restart_encoding={encoding}")
        sound = header.get("mhd", {}).get("iso_sound_speed")
        if sound is None or not np.isfinite(float(sound)) or float(sound) <= 0:
            raise ValueError(f"{mode}: explicitly retain positive mhd/iso_sound_speed in both inputs")
    mismatches, exempted = [], []
    for section in sorted(set(active) | set(passive)):
        for key in sorted(set(active.get(section, {})) | set(passive.get(section, {}))):
            a, p = active.get(section, {}).get(key), passive.get(section, {}).get(key)
            if a is not None and p is not None and canonical(a) == canonical(p):
                continue
            row = {"parameter": f"{section}/{key}", "active": a, "passive": p}
            if (section, key) in MODE_KEYS:
                row["reason"] = "declared active/passive model or restart encoding"
                exempted.append(row)
            elif administrative(section, key):
                row["reason"] = "output or runtime administration"
                exempted.append(row)
            else:
                mismatches.append(row)
    if mismatches:
        raise ValueError("Unmatched effective parameters: "+json.dumps(mismatches))
    return {"matched": True, "exempted_differences": exempted,
            "policy": "all embedded effective parameters agree except declared model flags, "
                      "basename, output blocks and enumerated runtime stop/diagnostic controls; "
                      "iso_sound_speed, mesh, reconstruction, integration, closure and forcing agree",
            "adaptive_timestep": "actual active/passive timesteps may differ with their dynamics; "
                                 "statistics use the same physical-time interval"}


def prescan(run, start, end, expected_mode):
    """Read physical primitive values, never estimate shared bins by rebinning."""
    run = Path(run).resolve()
    metadata_path = run/"benchmark_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    segments = single.load_segments(run, metadata, metadata_path)
    paths, _, discarded = single.segment_paths(
        segments, "snapshots", "**/*.mhd_w_bcc.*.bin", True)
    paths = single.bracketing_paths(paths, start, end)
    if not paths:
        raise ValueError(f"{expected_mode}: no snapshots bracket the requested interval")
    limits, records, reference, geometry, times = {}, [], None, None, []
    for path in paths:
        fields, lengths, header, info = single.read_uniform(path, PDF_FIELDS)
        model = single.infer_model(header, metadata)
        if model["passive"] != (expected_mode == "passive"):
            raise ValueError(f"{expected_mode}: retained snapshot has the wrong model flag")
        current_geometry = (info["shape_zyx"], list(lengths))
        if reference is None:
            reference, geometry = header, current_geometry
            for segment in segments:
                if segment["input"] is not None:
                    single.crosscheck_input(header, segment["input"])
        else:
            # Reuse the pair matcher to check all nonadministrative keys within
            # each run, after normalizing the deliberate mode difference.
            lhs, rhs = [dict((section, dict(values)) for section, values in h.items())
                        for h in (reference, header)]
            lhs["mhd"]["passive"], rhs["mhd"]["passive"] = "false", "true"
            for h, mode in ((lhs, "false"), (rhs, "true")):
                if "passive_delta" in h.get("problem", {}):
                    h["problem"]["passive_delta"] = mode
            compare_parameters(lhs, rhs)
            if current_geometry != geometry:
                raise ValueError(f"{expected_mode}: mesh geometry changed within the interval")
        b2 = sum(fields[name]**2 for name in ("bcc1", "bcc2", "bcc3"))
        if np.min(b2) <= 0 or min(np.min(fields[name]) for name in ("eint", "p_perp")) <= 0:
            raise ValueError(f"{path}: strictly positive field magnitude and pressures required")
        distributions = {"B_over_B0": np.sqrt(b2)/model["B0"],
                         "X": 2*(fields["p_perp"]-fields["eint"])/b2}
        for name, values in distributions.items():
            if not np.isfinite(values).all():
                raise ValueError(f"{path}: nonfinite distribution {name}")
            lo, hi = float(np.min(values)), float(np.max(values))
            old = limits.get(name, (lo, hi))
            limits[name] = min(lo, old[0]), max(hi, old[1])
        times.append(info["time"])
        records.append({"path": str(path), "info": info})
        print(f"prescan {expected_mode} t={info['time']:g}: {path.name}", flush=True)
    if min(times) > start or max(times) < end:
        raise ValueError(f"{expected_mode}: snapshots cover [{min(times)}, {max(times)}], "
                         f"not the full requested [{start}, {end}] interval; no extrapolation")
    return {"run": str(run), "metadata": single.retained_file(metadata_path),
            "header": reference, "geometry": geometry, "limits": limits,
            "snapshots": records, "discarded_snapshot_branches": discarded}


def shared_edges(scans, bins):
    edges, limits = {}, {}
    for name in ("B_over_B0", "X"):
        lo = min(scan["limits"][name][0] for scan in scans.values())
        hi = max(scan["limits"][name][1] for scan in scans.values())
        pad = max((hi-lo)*1e-8, 1e-12*max(1, abs(lo), abs(hi)))
        edges[name] = np.linspace(lo-pad, hi+pad, bins+1).tolist()
        limits[name] = [lo, hi]
    return {"schema_version": 1, "edges": edges, "observed_union_limits": limits,
            "method": "stream all cells from both retained, endpoint-bracketing intervals; "
                      "histogram the original fields with these identical linear edges"}


def numbers(values):
    return np.asarray(values, dtype=float)


def summary_difference(active, passive):
    """Retain paired window contrasts without inferential confidence intervals."""
    a, p = numbers(active["mean"]), numbers(passive["mean"])
    if a.shape != p.shape:
        raise ValueError("cannot compare summaries with different shapes")
    result = {"active_minus_passive": (a-p).tolist(),
              "active_over_passive": np.divide(a, p, out=np.full_like(a, np.nan),
                                               where=np.isfinite(p) & (p != 0)).tolist(),
              "block_difference_means": [], "block_difference_sd": None,
              "block_difference_min": None, "block_difference_max": None}
    if (active.get("block_starts") == passive.get("block_starts")
            and active.get("block_duration") == passive.get("block_duration")):
        aa, pp = numbers(active.get("block_means", [])), numbers(passive.get("block_means", []))
        if aa.shape != pp.shape:
            raise ValueError("matching time blocks have inconsistent data shapes")
        delta = aa-pp
        result["block_difference_means"] = delta
        if len(delta) >= 2:
            finite = np.all(np.isfinite(delta), axis=0)
            # Undefined ratios stay undefined; do not silently drop blocks.
            result["block_difference_sd"] = np.where(finite, np.std(delta, axis=0, ddof=1), np.nan).tolist()
            result["block_difference_min"] = np.where(finite, np.min(delta, axis=0), np.nan).tolist()
            result["block_difference_max"] = np.where(finite, np.max(delta, axis=0), np.nan).tolist()
    else:
        result["block_comparison_reason"] = "physical block boundaries differ"
    return single.json_value(result)


def pair_products(data):
    active, passive = data["active"], data["passive"]
    products = {}
    for category in ("PDFs", "spectra", "scalars"):
        products[category] = {}
        for key in sorted(set(active[category]) & set(passive[category])):
            a, p = active[category][key], passive[category][key]
            coordinate = "edges" if category == "PDFs" else "k" if category == "spectra" else None
            if coordinate and not np.array_equal(a[coordinate], p[coordinate]):
                raise ValueError(f"{category}.{key}: common bins or spectral shells differ")
            products[category][key] = {"active": a, "passive": p,
                                       "contrast": summary_difference(a, p)}
    products["pressure_balance_by_scale"] = {
        mode: data[mode]["pressure_balance_by_scale"] for mode in MODES}
    balance = products["pressure_balance_by_scale"]
    if not np.array_equal(balance["active"]["k"], balance["passive"]["k"]):
        raise ValueError("pressure-balance spectral shells differ")
    balance["contrast"] = {band: {key: summary_difference(balance["active"][band][key],
                                                          balance["passive"][band][key])
                                    for key in ("Paa", "Pbb", "Pab", "C", "R")}
                           for band in ("all_parallel", "full_k_resolved")}
    return products


def analysis_provenance(argv):
    """Identify this comparison independently of simulation/supplementary code."""
    source = single.retained_file(Path(__file__))
    root = Path(__file__).resolve().parents[1]
    try:
        revision = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
        status = subprocess.check_output(
            ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=normal"],
            text=True).splitlines()
        error = None
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        revision, status, error = None, None, str(exc)
    return {"analysis_revision": revision,
            "analysis_tree_dirty": None if status is None else bool(status),
            "analysis_tree_changes": status, "git_provenance_error": error,
            "comparison_script": source, "comparison_source_sha256": source["sha256"],
            "single_run_script": single.retained_file(Path(single.__file__)),
            "python": sys.version, "numpy": np.__version__,
            "argv": sys.argv if argv is None else argv}


def physical_evidence(data, integrity):
    """No automatic diagnostic check is a reviewed scientific conclusion."""
    reasons = ["Scientific interpretation has not been reviewed; parameter matching and "
               "successful analysis do not establish physical validation.",
               "One finite-time realization and dependent time blocks do not establish "
               "independent uncertainty, spatial convergence or a universal spectral interval."]
    if integrity["failed_runs"]:
        reasons.append("Failed simulation provenance: "+", ".join(integrity["failed_runs"])+".")
    for mode in MODES:
        purpose = data[mode].get("provenance", {}).get("retained_metadata", {}).get("classification")
        if purpose:
            reasons.append(f"{mode} retained run classification: {purpose}.")
        reasons.extend(f"{mode}: {reason}" for reason in data[mode].get("adequacy", {}).get("reasons", []))
    return {"classification": "inconclusive", "scientific_review_completed": False,
            "reasons": reasons}


def figure_status(data):
    """Carry retained simulation purpose and failure evidence onto every figure."""
    failed = [mode for mode in MODES if any(row.get("returncode") not in (0, None)
              for row in data[mode].get("simulation_integrity", {}).get("segments", []))]
    purposes = [data[mode].get("provenance", {}).get("retained_metadata", {}).get("classification", "")
                for mode in MODES]
    labels = []
    if failed:
        labels.append("FAILED RUN: "+", ".join(failed))
    if any("plumbing" in purpose.lower() for purpose in purposes):
        labels.append("PLUMBING ONLY")
    labels.append("Physical evidence inconclusive")
    return " | ".join(labels)


def make_figures(data, output):
    """Compact paired plots; all formulas and extended budgets stay supplementary."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import dbfplot

    audits = {}
    colors = dbfplot.colors("dark2", None)
    styles = {"active": "-", "passive": "--"}
    start, end = data["active"]["requested_window"]
    duration = data["active"]["sampling"]["block_duration"]
    band_note = f"Shading: range of {duration:g}-unit time averages."
    scale_note = (r"All $k_\parallel$; dotted: 8-cell transverse guide; "
                  r"gray: $k_\perp\geq0.75\,k_{\rm Ny}$ heuristic band.")
    status_note = f"Time [{start:g}, {end:g}] | "+figure_status(data)
    mode_handles = [Line2D([], [], color="0.25", ls=styles[mode], label=mode.capitalize())
                    for mode in MODES]

    def band(ax, x, rec, color, log=False, alpha=.10, first=0):
        if rec.get("block_min") is None:
            return
        x = x[first:]
        low, high = numbers(rec["block_min"])[first:], numbers(rec["block_max"])[first:]
        valid = np.isfinite(low) & np.isfinite(high)
        if log:
            valid &= (low > 0) & (high > 0)
        ax.fill_between(x, low, high, where=valid, color=color, alpha=alpha, linewidth=0)

    def mode_legend(ax, loc="upper right"):
        quantity_legend = ax.get_legend()
        if quantity_legend is not None:
            ax.add_artist(quantity_legend)
        legend = ax.legend(handles=mode_handles, loc=loc, frameon=False)
        ax.add_artist(legend)

    def finish(fig, name, note):
        fig.set_layout_engine("constrained", rect=(0, .15, 1, .745))
        fig.text(.5, .955, status_note, ha="center", va="top")
        fig.text(.5, .025, note, ha="center", va="bottom")
        for ax in fig.axes:
            for axis in ("x", "y"):
                if getattr(ax, f"get_{axis}scale")() == "log":
                    dbfplot.suppress_minor_labels(ax, axis=axis)
        exported = dbfplot.savefig(fig, output/name, profile="paper", strict=True,
                                   ignore=("figure.size",))
        audits[name] = {"ok": exported.report.ok,
                        "remaining_findings": list(exported.report.codes),
                        "figsize_inches": fig.get_size_inches().tolist(), "dpi": fig.dpi,
                        "exemptions": {"figure.size": "Multiple scientific panels need a larger canvas."},
                        "files": [p.name for p in exported.paths]}
        plt.close(fig)

    def spectrum(ax, key, color, label):
        for mode in MODES:
            rec = data[mode]["spectra"].get(key)
            if rec is None:
                continue
            k, mean = numbers(rec["k"]), numbers(rec["mean"])
            valid = (np.arange(len(k)) > 0) & (k > 0) & (mean > 0) & np.isfinite(mean)
            ax.loglog(k[valid], mean[valid], color=color, ls=styles[mode],
                      label=label if mode == "active" or key == "dynamic_isothermal_pressure" else "_nolegend_")
            band(ax, k, rec, color, log=True, first=1)
        ax.set_xlabel(r"$k_\perp$ [rad / length]")

    def scales(ax):
        reference = data["active"]
        ax.axvline(reference["scale_bands"]["resolved_kmax"], color="0.65", ls=":")
        ax.axvspan(reference["scale_bands"]["cutoff_band_kmin"],
                   reference["spectra"]["kinetic"]["k"][-1], color="0.7", alpha=.10)

    def guide(ax, key, exponent):
        reference = data["active"]
        rec = reference["spectra"][key]
        k, power = numbers(rec["k"]), numbers(rec["mean"])
        lo = max(1.4*reference["model"]["forcing_kmax"], k[0])
        hi = min(reference["scale_bands"]["resolved_kmax"], k[-1])
        valid = (np.arange(len(k)) > 0) & (k > 0) & (power > 0) & np.isfinite(power)
        if hi <= lo or np.count_nonzero(valid) < 2:
            return
        anchor = np.sqrt(lo*hi)
        amplitude = 2*np.exp(np.interp(np.log(anchor), np.log(k[valid]), np.log(power[valid])))
        x = np.geomspace(lo, hi, 32)
        ax.loglog(x, amplitude*(x/anchor)**exponent, color="0.35", ls="-.")
        label = r"$k_\perp^{-5/3}$" if exponent == -5/3 else r"$k_\perp^{1/3}$"
        ax.annotate(label, (anchor, amplitude), xytext=(0, 5), textcoords="offset points",
                    ha="center", va="bottom")

    with dbfplot.style_context():
        fig, axes = dbfplot.subplots(1, 2, figsize=(10, 4.2))
        for ax, key, color, xlabel in zip(axes, ("B_over_B0", "X"), (colors[2], colors[5]),
                (r"$|\mathbf{B}|/B_0$", r"$X=2(p_\perp-p_\parallel)/B^2$")):
            for mode in MODES:
                rec = data[mode]["PDFs"][key]
                edges = numbers(rec["edges"])
                centers = .5*(edges[:-1]+edges[1:])
                ax.plot(centers, np.ma.masked_less_equal(rec["mean"], 0),
                        color=color, ls=styles[mode], label=mode.capitalize())
                band(ax, centers, rec, color, log=True)
            ax.set(xscale="linear", yscale="log", xlabel=xlabel, ylabel="Volume PDF")
            ax.legend(frameon=False)
        for name in ("mirror", "firehose"):
            axes[1].axvline(data["active"]["model"][name+"_X"], color="0.6", ls=":")
        finish(fig, "pdfs", "Common linear bins. "+band_note)

        fig, axes = dbfplot.subplots(1, 2, figsize=(10, 4.2), sharey=True)
        width = data["active"]["sampling"]["near_threshold_halfwidth_X"]
        for ax, instability in zip(axes, ("mirror", "firehose")):
            threshold = data["active"]["model"][instability+"_X"]
            sign = ">" if instability == "mirror" else "<"
            difference = f"X-{threshold:g}" if threshold >= 0 else f"X+{-threshold:g}"
            for suffix, color, label in (("strict", colors[1], f"Strict: $X{sign}{threshold:g}$"),
                                         ("near", colors[0], rf"Near: $|{difference}|\leq {width:g}$")):
                for mode in MODES:
                    rows = data[mode]["snapshots"]
                    t = numbers([row["info"]["time"] for row in rows])
                    y = numbers([row["scalars"][instability+"_"+suffix] for row in rows])
                    ax.plot(t, y, color=color, ls=styles[mode],
                            label=label if mode == "active" else "_nolegend_")
            ax.text(.04, .94, instability.capitalize(), transform=ax.transAxes, va="top")
            ax.set(xlabel="Time [code units]", ylabel="Volume fraction", xlim=(start, end))
            ax.set_ylim(bottom=0)
            ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0, .82))
        mode_legend(axes[1])
        finish(fig, "occupancy", "Strict crossings and symmetric near bands are distinct; snapshot rounding matters.")

        fig, axes = dbfplot.subplots(1, 3, figsize=(13, 4.3))
        for key, color, label in (("kinetic", colors[0], r"$E_K$"),
                                   ("magnetic", colors[1], r"$E_B$")):
            spectrum(axes[0], key, color, label)
        for key, color, label in (("p_parallel", colors[2], r"$p_\parallel$"),
                                   ("p_perp", colors[3], r"$p_\perp$"),
                                   ("magnetic_pressure", colors[1], r"$p_B$"),
                                   ("dynamic_isothermal_pressure", colors[6], r"$c_{\rm iso}^2\rho$ (passive flow)")):
            spectrum(axes[1], key, color, label)
        for key, color, label in (("S_parallel", colors[2], r"$\hat b\cdot G\hat b$"),
                                   ("parallel_gradient_u_perp", colors[0], r"$PG\hat b$"),
                                   ("perp_gradient_u_parallel", colors[1], r"$\hat b\cdot GP$"),
                                   ("perp_gradient_u_perp", colors[3], r"$PGP$")):
            spectrum(axes[2], key, color, label)
        for ax, label in zip(axes, ("Energy spectrum", "Pressure spectrum", "Projected-gradient spectrum")):
            ax.set_ylabel(label)
            scales(ax)
            ax.legend(frameon=False, loc="lower left")
        guide(axes[0], "kinetic", -5/3)
        guide(axes[1], "p_perp", -5/3)
        guide(axes[2], "perp_gradient_u_perp", 1/3)
        mode_legend(axes[0])
        finish(fig, "spectra", scale_note+"\nSlopes are eye guides. "+band_note)

        fig, axes = dbfplot.subplots(1, 2, figsize=(10, 4.2))
        for ax, key, ylabel in zip(axes, ("C", "R"),
                                  ("Signed pressure correlation, $C$", "Pressure residual, $R$")):
            for mode in MODES:
                record = data[mode]["pressure_balance_by_scale"]
                rec = record["all_parallel"][key]
                k, mean = numbers(record["k"]), numbers(rec["mean"])
                valid = (np.arange(len(k)) > 0) & (k > 0) & np.isfinite(mean)
                ax.semilogx(k[valid], mean[valid], color=colors[3], ls=styles[mode],
                            label=mode.capitalize())
                band(ax, k, rec, colors[3], first=1)
            ax.set(xlabel=r"$k_\perp$ [rad / length]", ylabel=ylabel)
            ax.legend(frameon=False)
            scales(ax)
        axes[0].axhline(-1, color="0.6", ls=":")
        axes[0].axhline(0, color="0.6", ls="-.")
        axes[0].set_ylim(-1.05, 1.05)
        axes[1].axhline(1, color="0.6", ls=":")
        axes[1].set_ylim(bottom=0)
        finish(fig, "pressure_balance_scale", scale_note+"\n"+
               r"Thermal $p_\perp$, $p_B$: compensation at $C=-1$, $R=0$. "+band_note)
    return audits


def write_report(result, data, output):
    start, end = result["requested_window"]
    lines = ["# Matched active/passive CGL-LF comparison", "",
             f"Same physical interval **[{start:g}, {end:g}]**, with {result['block_duration']:g}-unit blocks. "
             "Effective physics, grid, reconstruction, integration and forcing parameters match; "
             "the active/passive model flags differ. Actual adaptive timesteps can differ.", "",
             "| Figure | Vector copy |", "| --- | --- |"]
    lines[2:2] = ["**Physical evidence: inconclusive.** These descriptive diagnostics require "
                  "scientific review; matching inputs and producing figures are not physical validation.", ""]
    integrity = result["simulation_integrity"]
    if integrity["failed_runs"]:
        lines[2:2] = ["**Failed simulation provenance: "+", ".join(integrity["failed_runs"])+
                      ". These retained-interval plots are diagnostic only; the pair is not a "
                      "successful benchmark or physical-validation result.**", ""]
    elif integrity["classification"] != "consistent":
        lines[2:2] = ["**Simulation integrity requires review.** See the retained numerical "
                      "health and completion records below; plots alone do not qualify the pair.", ""]
    for name in ("pdfs", "occupancy", "spectra", "pressure_balance_scale"):
        lines.append(f"| [{name.replace('_', ' ')}]({name}.png) | [PDF]({name}.pdf) |")
    lines += ["", "Colors identify quantities; solid curves are active and dashed curves passive. "
              "PDFs use identical linear edges determined from all retained cells in both intervals. "
              "Shading is the range of contiguous block means, not an independent-sample confidence interval. "
              "Spectral slopes are visual guides, not fitted acceptance criteria. The perpendicular zero "
              "shell (which contains pure-parallel modes) is omitted from logarithmic k plots and retained "
              "in metrics and Parseval sums; its bin midpoint is not a nonzero transverse mode.", "",
              "Main spectra and C/R curves sum **all parallel wavenumbers**. The dotted vertical line "
              "is the eight-cell transverse guide, k_perp = k_Ny/4; gray marks the heuristic "
              "k_perp >= 0.75 k_Ny band, not a measured numerical cutoff. Small k_perp alone does not "
              "guarantee a resolved full wavevector. Separately retained full_k_resolved C/R metrics "
              "select 0 < |k| <= min(k_Ny,x, k_Ny,y, k_Ny,z)/4 before shell aggregation; those filtered "
              "curves are not the curves displayed here. None of these guides establishes convergence.", "",
              "| Quantity | Active mean | Passive mean | Active − passive |", "| --- | ---: | ---: | ---: |"]
    for key in ("deltaB_rms_over_B0", "Mach_isotropic_proxy", "u_parallel_fraction", "S_parallel_rms",
                "mirror_strict", "firehose_strict", "mirror_near", "firehose_near",
                "pressure_correlation", "pressure_normalized_residual_variance",
                "thermal_density", "beta_volume_mean"):
        a = data["active"]["scalars"].get(key, {}).get("mean")
        p = data["passive"]["scalars"].get(key, {}).get("mean")
        av, pv = "unavailable" if a is None else f"{a:.6g}", "unavailable" if p is None else f"{p:.6g}"
        delta = "unavailable" if a is None or p is None else f"{a-p:.6g}"
        lines.append(f"| {key} | {av} | {pv} | {delta} |")
    passive_mach = data["passive"]["scalars"].get("Mach_isothermal", {}).get("mean")
    mach_text = "unavailable" if passive_mach is None else f"{passive_mach:.6g}"
    lines += ["", "Mach_isotropic_proxy is the common thermal sound-speed diagnostic. "
              f"The passive flow's separately defined isothermal Mach number is **{mach_text}**; "
              "its sound speed is the retained iso_sound_speed of the actual passive dynamics."]
    missing_pressure = {mode: [key for key in ("pressure_correlation", "pressure_normalized_residual_variance")
                               if key not in data[mode]["scalars"]] for mode in MODES}
    if any(missing_pressure.values()):
        lines += ["", "Unavailable global pressure ratios remain undefined when a required sampling or "
                  "endpoint-bracketing state has zero pressure variance. No undefined value is skipped "
                  "or interpolated across; detailed reasons and times remain in each supplementary report."]
    lines += ["", "Pressure-balance curves use signed cross powers of thermal perpendicular and magnetic "
              "pressure, averaged before forming C and R. The passive momentum pressure is instead "
              "c_iso² rho; its spectrum is labeled separately. Undefined spectral ratios are omitted, "
              "with null values retained in metrics.json.", "",
              "The common forcing parameters and OU seed define the same innovation prescription. "
              "The applied acceleration normalization can differ because it depends on the evolving flow. "
              "Passive invariant advection does not add irreversible shock heating, while the active "
              "calculation evolves total energy. Their thermal difference therefore is not an isolated "
              "measurement of anisotropy feedback.", "",
              "This is one finite-time paired realization. Block differences and variability are descriptive; "
              "no independence assumption, confidence interval or automatic tight physics tolerance is imposed.", "",
              "Interpretation limits:", ""]
    lines += ["- "+reason for reason in result["physical_evidence"]["reasons"]]
    lines += ["",
              "## Supplementary diagnostics", "",
              "[Active report](active/report.md) and [passive report](passive/report.md) retain forcing, "
              "conservation where applicable, curvature, Mach definitions, numerical health and longer formulas. "
              "The paired main figures do not replace those checks.", ""]
    for mode in MODES:
        rec = data[mode]
        health = rec.get("simulation_integrity", {})
        lines.append(f"- **{mode.capitalize()}**: {rec['sampling']['snapshots']} snapshots; "
                     f"{rec['sampling']['complete_blocks']} complete blocks; "
                     f"analysis normalization {rec['numerical_integrity']['classification']}; "
                     f"simulation health {health.get('classification', 'see supplementary metrics')}.")
        for caution in rec.get("temporal_dependence_cautions", []):
            lines.append(f"  {caution}.")
        for reason in rec.get("adequacy", {}).get("reasons", []):
            lines.append(f"  Sampling qualification: {reason}.")
    lines += ["", "[metrics.json](metrics.json) retains source/input hashes, exact matching exemptions, "
              "shared PDF edges, all numerical contrasts and block variability. "
              "[figure-audit.json](figure-audit.json) records strict dbfplot checks; only the multi-panel "
              "canvas size is exempted.", "",
              f"Comparison analysis revision: `{result['provenance']['analysis_revision']}`; "
              f"dirty worktree: `{result['provenance']['analysis_tree_dirty']}`; "
              f"comparison source SHA-256: `{result['provenance']['comparison_source_sha256']}`. "
              "These identify the comparison separately from simulation and supplementary analysis provenance.", ""]
    (output/"report.md").write_text("\n".join(lines))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("active_run", type=Path)
    parser.add_argument("passive_run", type=Path)
    parser.add_argument("--time-start", required=True, type=float)
    parser.add_argument("--time-end", required=True, type=float)
    parser.add_argument("--block-duration", default=2., type=float)
    parser.add_argument("--near-width", default=.05, type=float)
    parser.add_argument("--pdf-bins", default=128, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    if (not all(np.isfinite(x) for x in (args.time_start, args.time_end, args.block_duration, args.near_width))
            or args.time_end <= args.time_start or args.block_duration <= 0
            or args.near_width <= 0 or args.pdf_bins < 8):
        parser.error("require a finite ordered interval, positive block/band widths, and >=8 bins")
    if args.active_run.resolve() == args.passive_run.resolve():
        parser.error("active and passive run directories must differ")
    output = args.output_dir.resolve()
    scans = {mode: prescan(run, args.time_start, args.time_end, mode)
             for mode, run in zip(MODES, (args.active_run, args.passive_run))}
    match = compare_parameters(scans["active"]["header"], scans["passive"]["header"])
    if scans["active"]["geometry"] != scans["passive"]["geometry"]:
        raise ValueError("active/passive grids or domain lengths differ")
    edges = shared_edges(scans, args.pdf_bins)
    output.mkdir(parents=True, exist_ok=True)
    edges_path = output/"shared_pdf_edges.json"
    edges_path.write_text(json.dumps(edges, indent=2, allow_nan=False)+"\n")
    data = {}
    for mode in MODES:
        single.main([scans[mode]["run"], "--time-start", str(args.time_start),
                     "--time-end", str(args.time_end), "--block-duration", str(args.block_duration),
                     "--near-width", str(args.near_width), "--pdf-bins", str(args.pdf_bins),
                     "--pdf-edges", str(edges_path), "--output-dir", str(output/mode)])
        data[mode] = json.loads((output/mode/"metrics.json").read_text())
        if data[mode].get("requested_window") != [args.time_start, args.time_end]:
            raise ValueError(f"{mode}: unexpected supplementary analysis interval")
        for row in scans[mode]["snapshots"]:
            for item in row["info"]["files"]:
                if single.retained_file(Path(item["path"])) != item:
                    raise ValueError("snapshot changed between PDF prescan and analysis")
    products = pair_products(data)
    failed = [mode for mode in MODES if any(row.get("returncode") not in (0, None)
               for row in data[mode].get("simulation_integrity", {}).get("segments", []))]
    health_classes = [data[mode].get("simulation_integrity", {}).get("classification", "inconclusive")
                      for mode in MODES]
    integrity = {"classification": "concerning" if failed or "concerning" in health_classes else
                 "consistent" if all(value == "consistent" for value in health_classes) else "inconclusive",
                 "failed_runs": failed, "scope": "simulation completion and retained numerical health; "
                 "separate from parameter matching and descriptive contrast"}
    result = {"schema_version": 1, "classification": "descriptive matched comparison",
              "requested_window": [args.time_start, args.time_end], "block_duration": args.block_duration,
              "parameter_match": match, "shared_PDF_edges": edges, "products": products,
              "simulation_integrity": integrity,
              "physical_evidence": physical_evidence(data, integrity),
              "figure_status": figure_status(data),
              "main_spectral_plot_conventions": {
                  "parallel_modes": "all k_parallel; full_k_resolved comparisons remain in metrics only",
                  "zero_perpendicular_shell": "omitted from log-k plots; retained in metrics and Parseval sums",
                  "dotted_vertical": "k_perp = transverse inscribed Nyquist/4; eight-cell transverse guide",
                  "gray_band": "k_perp >= 0.75 transverse inscribed Nyquist; heuristic, not a measured cutoff",
                  "gradient_eye_guide_anchor": "perp_gradient_u_perp; illustrative slope +1/3",
                  "convergence": "No grid or spectral convergence is established by these guides."},
              "unavailable_scalar_comparisons": {mode: {key: reason for key, reason in
                  data[mode].get("unavailable_metric_reasons", {}).items() if key.startswith("scalars.")}
                  for mode in MODES},
              "block_interpretation": "paired contiguous physical-time means; dependent finite-time "
                                      "variability, not independent-sample confidence intervals",
              "runs": {mode: {"directory": scans[mode]["run"],
                              "metrics": single.retained_file(output/mode/"metrics.json"),
                              "simulation_revision": data[mode]["provenance"]["simulation_revision"],
                              "simulation_integrity": data[mode].get("simulation_integrity"),
                              "sampling": data[mode]["sampling"]} for mode in MODES},
              "provenance": {**analysis_provenance(argv), "PDF_prescans": scans}}
    result = single.json_value(result)
    (output/"metrics.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    audits = make_figures(data, output)
    (output/"figure-audit.json").write_text(json.dumps(audits, indent=2, allow_nan=False)+"\n")
    write_report(result, data, output)
    print(f"descriptive matched comparison: {output/'report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
