"""Plot small cached turbulence curves; never reread simulation fields.

python plot_comparison.py comparison.json --window 6 14 --output figures
The manifest lists cases with name, label, files (snapshot NPZs), and history
(one NPZ with ht/mirror/fire/backup), plus reference, kmax, and optional
dissipation_start and k_label. Snapshot times are canonical output times;
actual_time can preserve each writer's slightly offset simulation time.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np


SPECTRA = ("EK", "EB", "Pppar", "Ppperp", "PB2", "Pupar", "Pgrad", "PS", "Pdiv", "PQ",
           "gradient_parallel_parallel", "gradient_parallel_perp", "gradient_perp_parallel",
           "gradient_perp_perp", "gradient_total")


def time_weights(times):
    """Trapezoidal time weights, including half-weighted window endpoints."""
    times = np.asarray(times, dtype=float)
    assert times.ndim == 1 and len(times) > 0 and np.isfinite(times).all()
    if len(times) == 1:
        return np.ones(1)
    assert np.all(np.diff(times) > 0), "Duplicate or unordered snapshot times"
    return np.r_[times[1]-times[0], times[2:]-times[:-2], times[-1]-times[-2]] / (
        2*(times[-1]-times[0]))


def load_curves(manifest, directory, snapshot_time=None, window=None):
    curves, labels, selection = {}, {}, {}
    expected_times = None
    for case in manifest["cases"]:
        name = case["name"]
        assert name not in curves, name
        snapshots = []
        for filename in case["files"]:
            with np.load(directory / filename, allow_pickle=False) as data:
                nominal_time = float(data["time"])
                selected = (abs(nominal_time-snapshot_time) < 1e-9 if snapshot_time is not None
                            else window[0]-1e-9 <= nominal_time <= window[1]+1e-9)
                if selected:
                    snapshots.append((nominal_time, filename, dict(data)))
        snapshots.sort(key=lambda item: item[0])
        times = np.array([item[0] for item in snapshots])
        assert len(times), f"No selected snapshots for {name}"
        if snapshot_time is not None:
            assert len(times) == 1, f"Duplicate snapshot for {name} at {snapshot_time}"
        else:
            assert len(times) >= 2
            np.testing.assert_allclose(times[[0, -1]], window, rtol=0, atol=1e-9,
                                       err_msg=f"Missing averaging-window endpoint: {name}")
        if expected_times is None:
            expected_times = times
        np.testing.assert_array_equal(times, expected_times,
                                      err_msg="Cases must use matched snapshot times")
        weights = time_weights(times)
        data = {}
        for key in ("k", "b_centers", "beta_centers"):
            data[key] = snapshots[0][2][key]
            for _, _, snapshot in snapshots[1:]:
                np.testing.assert_array_equal(snapshot[key], data[key])
        for key in ("b_pdf", "beta_pdf", *SPECTRA):
            arrays = np.array([snapshot[key] for _, _, snapshot in snapshots])
            assert np.isfinite(arrays).all(), (name, key)
            data[key] = np.sum(weights[:, None]*arrays, axis=0)
        with np.load(directory / case["history"], allow_pickle=False) as history:
            visible = history["ht"] <= times[-1]+1e-9
            for key in ("ht", "mirror", "fire", "backup"):
                data[key] = history[key][visible]
        assert len(data["ht"]), f"No history rows before selected time: {name}"
        curves[name], labels[name] = data, case.get("label", name)
        selection[name] = {"files": [item[1] for item in snapshots], "times": times.tolist(),
                           "actual_times": [float(item[2].get("actual_time", item[0]))
                                            for item in snapshots], "weights": weights.tolist()}
    return curves, labels, selection


def make_plots(curves, output, labels=None, annotation="", reference="passive",
               kmax=None, dissipation_start=None, k_label=r"$kL_\perp/(2\pi)$"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path.home() / "dbfplot/src"))
    import dbfplot

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    dbfplot.use(latex=False, cycle="dark2")
    cases = list(curves)
    assert 2 <= len(cases) <= 3, "Expected two or three matched cases"
    labels = labels or {case: case for case in cases}
    assert reference in cases, "Specify the passive/reference case"
    colors = dict(zip(cases, dbfplot.colors("dark2", len(cases))))
    styles = dict(zip(cases, ("-", "--", "-.")))
    kmin = min(float(np.min(curves[case]["k"][curves[case]["k"] > 0])) for case in cases)
    kmax = kmax or min(float(np.max(curves[case]["k"])) for case in cases)
    assert kmax > kmin
    files = []

    def draw(ax, case, x, y, label=None, style=None, color=None):
        x, y = np.asarray(x), np.asarray(y)
        assert x.shape == y.shape and np.isfinite(x).all()
        assert np.isfinite(y).all(), (case, label)
        visible = y > 0 if ax.get_yscale() == "log" else np.ones(y.shape, bool)
        if ax.get_xscale() == "log":
            visible &= x > 0
        ax.plot(x, np.where(visible, y, np.nan), color=color or colors[case],
                linestyle=style or styles[case], label=label or labels[case])

    def all_cases(ax, xkey, ykey, factor=1):
        for case in cases:
            draw(ax, case, curves[case][xkey], factor*curves[case][ykey])

    def spectral_axes(ax, ylabel):
        ax.set(xscale="log", yscale="log", xlim=(kmin, kmax),
               xlabel=k_label, ylabel=ylabel)
        if dissipation_start is not None:
            assert kmin < dissipation_start < kmax
            ax.axvspan(dissipation_start, kmax, color="0.5", alpha=0.10, linewidth=0)

    def spectrum(ax, case, key, factor=1, **kwargs):
        data = curves[case]
        visible = (data["k"] >= kmin) & (data["k"] <= kmax)
        draw(ax, case, data["k"][visible], (factor*data[key])[visible], **kwargs)

    def ratio(ax, case, key, denominator, label, marker=None):
        data = curves[case]
        # Resolution comparisons share low-k shells, while their cutoffs differ.
        count = min(len(data["k"]), len(denominator))
        k, numerator, denominator = data["k"][:count], data[key][:count], denominator[:count]
        visible = (k >= kmin) & (k <= kmax) & (denominator > 0) & (numerator > 0)
        ax.plot(k[visible], numerator[visible]/denominator[visible],
                color=colors[case], linestyle=styles[case], label=label,
                marker=marker, markersize=3, markevery=12)

    def save(fig, name):
        if name != "02_threshold_history" and annotation:
            fig.axes[0].text(0.97, 0.97, annotation, ha="right", va="top",
                             transform=fig.axes[0].transAxes, fontsize=9)
        for ax in fig.axes:
            dbfplot.suppress_minor_labels(ax)
        exported = dbfplot.savefig(fig, output / name, profile="paper", strict=True,
                                   ignore=("figure.size",))
        files.extend(str(path) for path in exported.paths)
        plt.close(fig)

    fig, axes = dbfplot.subplots(1, 2, figsize=(10, 4), layout="constrained")
    axes[0].set(xscale="log", yscale="log", xlabel=r"$|B|/B_0$",
                ylabel=r"$p(|B|/B_0)$")
    all_cases(axes[0], "b_centers", "b_pdf")
    b_limits = (10**-0.5, 10**0.5)
    visible_pdf = []
    for case in cases:
        x, y = curves[case]["b_centers"], curves[case]["b_pdf"]
        positive = y > 0
        visible_pdf.extend(y[positive & (x >= b_limits[0]) & (x <= b_limits[1])])
        if positive.any():
            visible_pdf.extend(10**np.interp(np.log10(b_limits), np.log10(x[positive]),
                                            np.log10(y[positive])))
    assert visible_pdf, "No nonzero magnetic PDF bins"
    axes[0].set(xlim=b_limits, ylim=(min(visible_pdf)/1.5, max(visible_pdf)*1.5))
    axes[0].legend(frameon=False, fontsize=9)
    axes[1].set(yscale="log", xlabel=r"$\beta\Delta=2(p_\perp-p_\parallel)/B^2$",
                ylabel=r"$p(\beta\Delta)$", xlim=(-2.1, 2.1))
    all_cases(axes[1], "beta_centers", "beta_pdf")
    for value, label in ((-1.4, "Soft onset"), (1.0, "_nolegend_")):
        axes[1].axvline(value, color="0.35", linestyle="--", linewidth=1, label=label)
    for value, label in ((-2.0, "Backup onset"), (2.0, "_nolegend_")):
        axes[1].axvline(value, color="0.35", linestyle=":", linewidth=1, label=label)
    axes[1].set_ylim(bottom=2e-2)
    axes[1].legend(frameon=False, fontsize=9)
    save(fig, "01_pdfs")

    fig, axes = dbfplot.subplots(1, 3, figsize=(12, 3.8), layout="constrained")
    end = max(float(curves[case]["ht"][-1]) for case in cases)
    for ax, key, label in zip(axes, ("mirror", "fire", "backup"),
                              ("Mirror onset [%]", "Firehose onset [%]", "Backup region [%]")):
        all_cases(ax, "ht", key)
        ax.set(xlabel=r"$t$", ylabel=label, xlim=(0, end), ylim=(0, None))
    axes[0].legend(frameon=False, fontsize=9)
    save(fig, "02_threshold_history")

    fig, axes = dbfplot.subplots(2, 1, figsize=(6.5, 7), sharex=True, layout="constrained")
    spectral_axes(axes[0], r"$E(k)$")
    spectral_axes(axes[1], "$[" + k_label.strip("$") + r"]^{5/3}E(k)$")
    axes[0].set_xlabel("")
    for case in cases:
        for key, style, symbol in (("EK", "-", "K"), ("EB", "--", "B")):
            label = labels[case] + rf": $E_{{\rm {symbol}}}$"
            spectrum(axes[0], case, key, label=label, style=style)
            spectrum(axes[1], case, key, factor=curves[case]["k"]**(5/3),
                     label=label, style=style)
    # Guide amplitude follows the displayed power; its exponent is not fitted.
    data = curves[reference]
    band = (data["k"] >= 4) & (data["k"] <= min(16, kmax)) & (data["EK"] > 0)
    if band.any():
        amplitude = np.median(data["EK"][band]*data["k"][band]**(5/3))
        guide_k = np.geomspace(4, min(32, kmax), 100)
        axes[0].plot(guide_k, amplitude*guide_k**(-5/3), color="0.3",
                     linestyle=":", linewidth=1.2)
        at = np.sqrt(guide_k[0]*guide_k[-1])
        axes[0].text(at, 1.5*amplitude*at**(-5/3), r"$k^{-5/3}$", color="0.3", fontsize=10)
        axes[1].plot(guide_k, np.full_like(guide_k, amplitude), color="0.3",
                     linestyle=":", linewidth=1.2)
    handles, names = axes[0].get_legend_handles_labels()
    order = list(range(0, 2*len(cases), 2)) + list(range(1, 2*len(cases), 2))
    axes[0].legend([handles[i] for i in order], [names[i] for i in order],
                   frameon=False, fontsize=9, loc="lower left", ncol=2)
    save(fig, "03_energy_spectra")

    quantity_colors = dbfplot.colors("dark2", 4)
    fig, axes = dbfplot.subplots(1, len(cases), figsize=(4.1*len(cases), 3.9),
                                sharey=True, layout="constrained")
    for index, (ax, case) in enumerate(zip(axes, cases)):
        spectral_axes(ax, r"$E(k)$" if index == 0 else "")
        for key, quantity, color, style in zip(
                ("Pppar", "Ppperp", "PB2"),
                (r"$p_\parallel$", r"$p_\perp$", r"$B^2/2$"), quantity_colors, ("-", "--", "-.")):
            spectrum(ax, case, key, label=quantity, color=color, style=style)
        ax.text(0.97, 0.89, labels[case], ha="right", va="top", transform=ax.transAxes, fontsize=10)
    axes[0].legend(frameon=False, fontsize=9, loc="lower left")
    save(fig, "04_pressure_spectra")

    fig, axes = dbfplot.subplots(2, 3, figsize=(12, 7.3), layout="constrained")
    for ax, key, label in zip(axes.flat, ("PS", "PQ", "Pdiv", "Pupar"),
                              (r"$P_S(k),\quad S=\hat b_i\hat b_j\partial_j u_i$",
                               r"$P_Q(k),\quad Q=S-(\nabla\cdot u)/3$",
                               r"$P_{\nabla\cdot u}(k)$", r"$\frac{1}{2}P_{u_\parallel}(k)$")):
        spectral_axes(ax, label)
        for case in cases:
            spectrum(ax, case, key, factor=0.5 if key == "Pupar" else 1)
    axes[0, 0].legend(frameon=False, fontsize=9)
    axes[1, 0].text(0.03, 0.05, r"$u_\parallel=\mathbf{u}\cdot\hat{\mathbf{b}}$",
                    transform=axes[1, 0].transAxes, fontsize=10)
    spectral_axes(axes[1, 1], r"$P_S(k)/P_{\nabla u}(k)$")
    for case in cases:
        ratio(axes[1, 1], case, "PS", curves[case]["Pgrad"], labels[case])
    spectral_axes(axes[1, 2], "Power / " + labels[reference])
    axes[1, 2].axhline(1, color="0.5", linewidth=1, linestyle=":")
    for case in cases:
        if case == reference:
            continue
        count = min(len(curves[case]["k"]), len(curves[reference]["k"]))
        np.testing.assert_array_equal(curves[case]["k"][:count], curves[reference]["k"][:count])
        for key, symbol, marker in (("PS", "S", None), ("Pupar", r"u_\parallel", "o")):
            ratio(axes[1, 2], case, key, curves[reference][key],
                  labels[case] + rf": ${symbol}$", marker)
    axes[1, 2].legend(frameon=False, fontsize=8)
    save(fig, "05_strain_spectra")

    gradient_keys = ("parallel_parallel", "parallel_perp", "perp_parallel", "perp_perp")
    gradient_labels = (r"$\nabla_\parallel u_\parallel$", r"$\nabla_\parallel\mathbf{u}_\perp$",
                       r"$\nabla_\perp u_\parallel$", r"$\nabla_\perp\mathbf{u}_\perp$")
    fig, axes = dbfplot.subplots(1, len(cases), figsize=(4.1*len(cases), 4),
                                sharey=True, layout="constrained")
    for index, (ax, case) in enumerate(zip(axes, cases)):
        spectral_axes(ax, r"$E_{\nabla u}(k)$" if index == 0 else "")
        for key, label, color, style in zip(gradient_keys, gradient_labels, quantity_colors,
                                           ("-", "--", "-.", ":")):
            spectrum(ax, case, "gradient_" + key, label=label, color=color, style=style)
        spectrum(ax, case, "gradient_total", label=r"Total $\nabla\mathbf{u}$",
                 color="0.5", style="-")
        ax.text(0.97, 0.89, labels[case], ha="right", va="top", transform=ax.transAxes, fontsize=10)
    axes[0].legend(frameon=False, fontsize=9, loc="lower left")
    save(fig, "06_gradient_components")
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, nargs="?")
    parser.add_argument("--output", type=Path)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--time", type=float)
    selection.add_argument("--window", type=float, nargs=2)
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        np.testing.assert_allclose(time_weights([6, 7, 9, 14]), [1/16, 3/16, 7/16, 5/16])
        np.testing.assert_allclose(np.dot(time_weights([6, 7, 9, 14]), [6, 7, 9, 14]), 10)
        np.testing.assert_array_equal(time_weights([14]), [1])
        for invalid in ([1, 1], [2, 1]):
            try:
                time_weights(invalid)
            except AssertionError:
                continue
            raise AssertionError("Invalid time sequence accepted")
        print("Uneven-cadence trapezoidal averaging and duplicate-time checks passed.")
        return
    if args.manifest is None or args.output is None or (args.time is None and args.window is None):
        parser.error("Supply manifest, --output, and either --time or --window")
    if args.window is not None and args.window[0] >= args.window[1]:
        parser.error("--window requires increasing endpoints")
    manifest = json.loads(args.manifest.read_text())
    curves, labels, selected = load_curves(manifest, args.manifest.parent, args.time, args.window)
    annotation = (rf"$t={args.time:g}$" if args.time is not None else
                  rf"${args.window[0]:g}\leq t\leq {args.window[1]:g}$ average")
    files = make_plots(curves, args.output, labels, annotation=annotation,
                       **{key: manifest[key] for key in
                          ("reference", "kmax", "dissipation_start", "k_label")
                          if key in manifest})
    (args.output / "plot_selection.json").write_text(json.dumps(selected, indent=2) + "\n")
    print("Saved figures:\n" + "\n".join(files))


if __name__ == "__main__":
    main()
