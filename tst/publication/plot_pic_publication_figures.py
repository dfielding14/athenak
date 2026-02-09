#!/usr/bin/env python3
"""Generate proxy and publication-physics figure bundles from archived artifacts."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
import sys

import numpy as np

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except ModuleNotFoundError as exc:
    raise SystemExit(
        "matplotlib is required for plotting. "
        "Install with: python3 -m pip install matplotlib"
    ) from exc

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vis" / "python"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import bin_convert_new as bin_convert  # noqa: E402
from scripts.particles.pic_analysis_utils import (  # noqa: E402
    fit_exponential_growth_windowed,
)

from pvtk_particles import read_particle_vtk  # noqa: E402


_PUB_COLORS = [
    "#1b9e77",  # teal
    "#d95f02",  # orange
    "#7570b3",  # violet
    "#e7298a",  # magenta
    "#66a61e",  # green
    "#e6ab02",  # mustard
    "#a6761d",  # brown
    "#666666",  # gray
]


def _configure_publication_style() -> None:
    plt.style.use("default")
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#222222",
            "axes.linewidth": 1.0,
            "axes.labelsize": 14,
            "axes.titlesize": 17,
            "axes.titleweight": "semibold",
            "axes.labelweight": "normal",
            "axes.prop_cycle": mpl.cycler(color=_PUB_COLORS),
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
            "mathtext.fontset": "dejavuserif",
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 5.0,
            "ytick.major.size": 5.0,
            "xtick.minor.size": 3.0,
            "ytick.minor.size": 3.0,
            "lines.linewidth": 2.2,
            "lines.markersize": 5.5,
            "grid.color": "#7f7f7f",
            "grid.linewidth": 0.8,
            "grid.alpha": 0.25,
            "legend.frameon": True,
            "legend.fancybox": False,
            "legend.framealpha": 0.95,
            "legend.edgecolor": "#505050",
            "legend.facecolor": "white",
            "legend.fontsize": 11,
            "legend.title_fontsize": 11,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.autolayout": False,
        }
    )


def _flatten_axes(*axes):
    out = []
    for item in axes:
        if item is None:
            continue
        if isinstance(item, np.ndarray):
            out.extend(item.ravel().tolist())
        elif isinstance(item, (list, tuple)):
            out.extend(item)
        else:
            out.append(item)
    return out


def _style_axes(*axes, grid: bool = True, log_grid: bool = False) -> None:
    for ax in _flatten_axes(*axes):
        for spine in ax.spines.values():
            spine.set_linewidth(1.0)
            spine.set_color("#222222")
        ax.tick_params(
            axis="both",
            which="major",
            direction="in",
            top=True,
            right=True,
            width=0.9,
            length=5,
        )
        ax.tick_params(
            axis="both",
            which="minor",
            direction="in",
            top=True,
            right=True,
            width=0.7,
            length=3,
        )
        if grid:
            ax.grid(True, which="major", linestyle="-", alpha=0.22)
            ax.grid(True, which="minor", linestyle=":", alpha=0.18)
        if not log_grid:
            ax.minorticks_on()


def _add_legend(ax, **kwargs):
    defaults = {
        "frameon": True,
        "framealpha": 0.95,
        "facecolor": "white",
        "edgecolor": "#505050",
    }
    defaults.update(kwargs)
    return ax.legend(**defaults)


def _annotation_bbox(alpha: float = 0.92) -> dict[str, object]:
    return {
        "boxstyle": "round,pad=0.24",
        "fc": "white",
        "ec": "#505050",
        "alpha": alpha,
    }


def _style_colorbar(cbar, label: str | None = None) -> None:
    if label is not None:
        cbar.set_label(label)
    cbar.ax.tick_params(which="both", direction="in", width=0.8, length=3)


def _save_figure(fig, out: Path, dpi: int) -> None:
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")


def _case_dir(run_dir: Path, case_id: str) -> Path:
    return run_dir / "cases" / case_id / "files"


def _case_bin_dir(run_dir: Path, case_id: str) -> Path:
    return _case_dir(run_dir, case_id) / "tst" / "build" / "src" / "bin"


def _case_src_dir(run_dir: Path, case_id: str) -> Path:
    return _case_dir(run_dir, case_id) / "tst" / "build" / "src"


def _case_pvtk_dir(run_dir: Path, case_id: str) -> Path:
    return _case_dir(run_dir, case_id) / "tst" / "build" / "src" / "pvtk"


def _select_case_dir(
    run_dir: Path,
    case_ids: tuple[str, ...],
    subdir: str,
) -> tuple[Path, str]:
    for case_id in case_ids:
        if subdir == "bin":
            cand = _case_bin_dir(run_dir, case_id)
        elif subdir == "src":
            cand = _case_src_dir(run_dir, case_id)
        elif subdir == "pvtk":
            cand = _case_pvtk_dir(run_dir, case_id)
        else:
            raise ValueError("unsupported case subdir: " + subdir)
        if cand.is_dir():
            return cand, case_id

    fallback = case_ids[0]
    if subdir == "bin":
        return _case_bin_dir(run_dir, fallback), fallback
    if subdir == "src":
        return _case_src_dir(run_dir, fallback), fallback
    return _case_pvtk_dir(run_dir, fallback), fallback


def _latest_file(pattern: str, root: Path) -> Path | None:
    matches = sorted(root.glob(pattern))
    return matches[-1] if matches else None


def _integrate_quantity(dataset: dict, quantity: str) -> float:
    dx1 = np.diff(dataset["x1f"])
    dx2 = np.diff(dataset["x2f"])
    dx3 = np.diff(dataset["x3f"])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    return float(np.sum(dataset[quantity] * dvol))


def _mode1_complex(dataset: dict, field: str) -> complex:
    values = np.asarray(dataset[field], dtype=float)
    x_mode = np.mean(values, axis=(0, 1))
    fluc = x_mode - np.mean(x_mode)

    x1f = np.asarray(dataset["x1f"], dtype=float)
    x1v = np.asarray(dataset["x1v"], dtype=float)
    length = float(x1f[-1] - x1f[0])
    phase = np.exp(-2.0j * np.pi * (x1v - x1f[0]) / length)
    return np.sum(fluc * phase) / x_mode.size


def _series_from_bin(
    files: list[Path],
    compute_value,
) -> tuple[np.ndarray, np.ndarray]:
    times = []
    values = []
    for path in sorted(files):
        data = bin_convert.read_binary_as_athdf(str(path))
        times.append(float(data["Time"]))
        values.append(float(compute_value(data)))
    return np.asarray(times, dtype=float), np.asarray(values, dtype=float)


def _dominant_frequency_nonuniform(
    times: np.ndarray,
    signal: np.ndarray,
    fmin: float = 0.05,
    fmax: float = 0.40,
    nsample: int = 2000,
) -> tuple[float, float]:
    t = np.asarray(times, dtype=float)
    y = np.asarray(signal, dtype=float)
    if t.size != y.size or t.size < 8:
        return 0.0, 0.0

    y = y - np.mean(y)
    freqs = np.linspace(fmin, fmax, nsample)
    best_f = float(freqs[0])
    best_amp = -1.0

    for freq in freqs:
        omega_t = 2.0 * np.pi * freq * t
        design = np.column_stack((np.sin(omega_t), np.cos(omega_t)))
        coeff, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        amp = float(np.hypot(coeff[0], coeff[1]))
        if amp > best_amp:
            best_amp = amp
            best_f = float(freq)

    return best_f, best_amp


def _trim_startup_zeros(
    times: np.ndarray,
    values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(times, dtype=float)
    v = np.asarray(values, dtype=float)
    if t.size != v.size or v.size == 0:
        return t, v

    vmax = float(np.max(np.abs(v)))
    threshold = max(1.0e-30, 1.0e-12 * vmax)
    mask = np.abs(v) > threshold
    if np.count_nonzero(mask) < 3:
        return t, v
    first = int(np.argmax(mask))
    return t[first:], v[first:]


def _fit_log_growth(
    times: np.ndarray,
    values: np.ndarray,
    floor: float,
) -> tuple[float, float] | None:
    t = np.asarray(times, dtype=float)
    v = np.asarray(values, dtype=float)
    if t.size != v.size or t.size < 4:
        return None

    logv = np.log(np.maximum(v + floor, floor))
    coeff = np.polyfit(t, logv, 1)
    gamma = float(coeff[0])
    fit = coeff[0] * t + coeff[1]
    ss_res = float(np.sum((logv - fit) ** 2))
    ss_tot = float(np.sum((logv - np.mean(logv)) ** 2))
    r2 = 1.0 if ss_tot <= 1.0e-30 else (1.0 - ss_res / ss_tot)
    return gamma, r2


def _fit_sinusoid_fixed_freq(
    times: np.ndarray,
    signal: np.ndarray,
    freq: float,
) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(times, dtype=float)
    y = np.asarray(signal, dtype=float)
    omega_t = 2.0 * np.pi * freq * t
    design = np.column_stack((np.sin(omega_t), np.cos(omega_t), np.ones_like(t)))
    coeff, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    fit = design @ coeff
    resid = y - fit
    return fit, resid


def _linear_fit_with_uncertainty(
    x: np.ndarray,
    y: np.ndarray,
) -> tuple[float, float, float]:
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    if xx.size < 3 or yy.size != xx.size:
        return 0.0, 0.0, 0.0

    coeff = np.polyfit(xx, yy, 1)
    slope = float(coeff[0])
    intercept = float(coeff[1])
    fit = slope * xx + intercept

    dof = xx.size - 2
    if dof <= 0:
        return slope, intercept, 0.0

    s2 = float(np.sum((yy - fit) ** 2) / dof)
    x_centered = xx - np.mean(xx)
    sxx = float(np.sum(x_centered * x_centered))
    if sxx <= 1.0e-30:
        return slope, intercept, 0.0

    slope_stderr = float(np.sqrt(s2 / sxx))
    return slope, intercept, slope_stderr


def _plot_entity_deposit_maps(run_dir: Path, fig_dir: Path, dpi: int) -> Path | None:
    cases = [
        (("entity_deposit_mink_publication", "entity_deposit_mink"), "Minkowski"),
        (("entity_deposit_reflect_publication", "entity_deposit_reflect"), "Reflect"),
    ]
    loaded = []
    for candidates, label in cases:
        bin_dir, _ = _select_case_dir(run_dir, candidates, "bin")
        fpath = _latest_file("*_np1.prtcl_d.*.bin", bin_dir)
        if fpath is None:
            fpath = _latest_file("*.prtcl_d.*.bin", bin_dir)
        if fpath is None:
            continue
        data = bin_convert.read_binary_as_athdf(str(fpath))
        pdens = np.asarray(data["pdens"], dtype=float)

        rho_std = None
        rho_path = _latest_file("*.prtcl_rho.*.bin", bin_dir)
        if rho_path is not None:
            rho_data = bin_convert.read_binary_as_athdf(str(rho_path))
            rho_std = float(np.std(np.asarray(rho_data["prtcl_rho"], dtype=float)))

        rank_diffs = []
        for rank in (2, 3, 4):
            rpath = _latest_file(f"*_np{rank}.prtcl_d.*.bin", bin_dir)
            if rpath is None:
                continue
            rdata = bin_convert.read_binary_as_athdf(str(rpath))
            rpdens = np.asarray(rdata["pdens"], dtype=float)
            if rpdens.shape == pdens.shape:
                rank_diffs.append((rank, float(np.max(np.abs(rpdens - pdens)))))

        loaded.append((label, np.mean(pdens, axis=0), rho_std, rank_diffs))
    if not loaded:
        return None

    fig, axes = plt.subplots(1, len(loaded), figsize=(6.4 * len(loaded), 4.6))
    if len(loaded) == 1:
        axes = [axes]
    _style_axes(*axes, grid=False)

    for idx, (ax, (label, arr, rho_std, rank_diffs)) in enumerate(zip(axes, loaded)):
        im = ax.imshow(arr, origin="lower", aspect="auto", cmap="viridis")
        panel = "abcdefghijklmnopqrstuvwxyz"[idx]
        ax.set_title(f"({panel}) {label} final $n_{{prtcl}}$ support (np1)")
        ax.set_xlabel("x1 index")
        ax.set_ylabel("x2 index")
        cbar = fig.colorbar(im, ax=ax, shrink=0.9)
        _style_colorbar(cbar, r"$n_{prtcl}$")
        if rho_std is not None:
            ax.text(
                0.02,
                0.98,
                r"$\sigma(\rho_{prtcl})=$" + f"{rho_std:.1e}",
                transform=ax.transAxes,
                va="top",
                fontsize=10,
                bbox=_annotation_bbox(),
            )
        if rank_diffs:
            text = "max |npN-np1| " + r"$n_{prtcl}$" + ":\n" + "\n".join(
                [f"np{rank}: {diff:.1e}" for rank, diff in rank_diffs]
            )
            ax.text(
                0.98,
                0.98,
                text,
                transform=ax.transAxes,
                va="top",
                ha="right",
                fontsize=10,
                bbox=_annotation_bbox(),
            )
    out = fig_dir / "01_entity_deposit_maps.png"
    fig.tight_layout()
    _save_figure(fig, out, dpi)
    plt.close(fig)
    return out


def _plot_em_vacuum_convergence(run_dir: Path, fig_dir: Path, dpi: int) -> Path | None:
    src_dir, _ = _select_case_dir(
        run_dir, ("em_vacuum_wave_publication", "em_vacuum_wave"), "src"
    )
    files = sorted(src_dir.glob("pic_em_vacuum_*np*-errs.dat"))
    if not files:
        return None

    data_by_np: dict[int, list[tuple[int, float]]] = {}
    pat = re.compile(r"pic_em_vacuum_(?:pub_)?np(\d+)_n(\d+)-errs\.dat")
    for fpath in files:
        match = pat.match(fpath.name)
        if match is None:
            continue
        nproc = int(match.group(1))
        nres = int(match.group(2))
        table = np.loadtxt(fpath)
        row = table if table.ndim == 1 else table[-1, :]
        l1 = float(row[4])
        data_by_np.setdefault(nproc, []).append((nres, l1))
    if not data_by_np:
        return None

    fig, ax = plt.subplots(figsize=(6, 4))
    for nproc, values in sorted(data_by_np.items()):
        values.sort()
        x = [item[0] for item in values]
        y = [item[1] for item in values]
        ax.loglog(x, y, marker="o", label=f"np{nproc}")
    ax.set_xlabel("Resolution Nx1")
    ax.set_ylabel("L1 error")
    ax.set_title("EM vacuum proxy convergence")
    _style_axes(ax, grid=True, log_grid=True)
    _add_legend(ax, loc="upper right")
    out = fig_dir / "02_em_vacuum_convergence.png"
    fig.tight_layout()
    _save_figure(fig, out, dpi)
    plt.close(fig)
    return out


def _plot_langmuir_trace(run_dir: Path, fig_dir: Path, dpi: int) -> Path | None:
    bin_dir, _ = _select_case_dir(
        run_dir,
        ("langmuir_frequency_publication", "langmuir_frequency_proxy"),
        "bin",
    )

    def _rank_series(rank: int) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
        jx_files = sorted(
            bin_dir.glob(f"pic_langmuir_freq_proxy_np{rank}.prtcl_jx.*.bin")
        )
        jy_files = sorted(
            bin_dir.glob(f"pic_langmuir_freq_proxy_np{rank}.prtcl_jy.*.bin")
        )
        if not jx_files or not jy_files:
            return None

        cycles = {}
        for fpath in jx_files:
            cyc = fpath.name.split(".")[-2]
            cycles.setdefault(cyc, {})["jx"] = fpath
        for fpath in jy_files:
            cyc = fpath.name.split(".")[-2]
            cycles.setdefault(cyc, {})["jy"] = fpath

        times = []
        jxv = []
        jyv = []
        for cyc in sorted(cycles):
            pair = cycles[cyc]
            if "jx" not in pair or "jy" not in pair:
                continue
            jx_data = bin_convert.read_binary_as_athdf(str(pair["jx"]))
            jy_data = bin_convert.read_binary_as_athdf(str(pair["jy"]))
            times.append(float(jx_data["Time"]))
            jxv.append(_integrate_quantity(jx_data, "prtcl_jx"))
            jyv.append(_integrate_quantity(jy_data, "prtcl_jy"))
        if not times:
            return None
        return (
            np.asarray(times, dtype=float),
            np.asarray(jxv, dtype=float),
            np.asarray(jyv, dtype=float),
        )

    np1 = _rank_series(1)
    if np1 is None:
        return None
    np2 = _rank_series(2)
    times, jxv, jyv = np1

    f_exp = 1.0 / (2.0 * np.pi)
    fit_jx, resid_jx = _fit_sinusoid_fixed_freq(times, jxv, f_exp)
    fit_jy, resid_jy = _fit_sinusoid_fixed_freq(times, jyv, f_exp)

    fig, axes = plt.subplots(
        2, 1, figsize=(7.8, 6.0), sharex=True, gridspec_kw={"height_ratios": [3, 1.3]}
    )
    ax = axes[0]
    ax_res = axes[1]
    _style_axes(ax, grid=True)
    _style_axes(ax_res, grid=True)
    ax.plot(times, jxv, label="Jx np1", lw=2.0)
    ax.plot(times, jyv, label="Jy np1", lw=2.0)
    ax.plot(times, fit_jx, "k--", lw=1.2, alpha=0.8, label="Jx analytic fit")
    ax.plot(times, fit_jy, "k:", lw=1.2, alpha=0.8, label="Jy analytic fit")
    if np2 is not None:
        t2, jx2, jy2 = np2
        ax.plot(t2, jx2, "--", label="Jx np2", alpha=0.9)
        ax.plot(t2, jy2, "--", label="Jy np2", alpha=0.9)

    ax.set_ylabel("Integrated current")
    ax.set_title("(a) Langmuir proxy current oscillation")
    _add_legend(ax, loc="center", ncol=2)

    fx, ax_amp = _dominant_frequency_nonuniform(
        times, jxv
    )
    fy, ay_amp = _dominant_frequency_nonuniform(
        times, jyv
    )
    if ax_amp > 0.0 and ay_amp > 0.0:
        extra = ""
        if np2 is not None:
            t2, jx2, jy2 = np2
            fx2, _ = _dominant_frequency_nonuniform(t2, jx2)
            fy2, _ = _dominant_frequency_nonuniform(t2, jy2)
            extra = f"\nnp2: f(Jx)={fx2:.4f}, f(Jy)={fy2:.4f}"
        ax.text(
            0.02,
            0.98,
            f"np1: f(Jx)={fx:.4f}, f(Jy)={fy:.4f}, f_exp={f_exp:.4f}" + extra,
            transform=ax.transAxes,
            va="top",
            fontsize=10,
            bbox=_annotation_bbox(),
        )

    ax_res.plot(times, resid_jx, label="Jx residual (np1)")
    ax_res.plot(times, resid_jy, label="Jy residual (np1)")
    if np2 is not None:
        t2, jx2, jy2 = np2
        fit_jx2, resid_jx2 = _fit_sinusoid_fixed_freq(t2, jx2, f_exp)
        fit_jy2, resid_jy2 = _fit_sinusoid_fixed_freq(t2, jy2, f_exp)
        _ = fit_jx2, fit_jy2
        ax_res.plot(t2, resid_jx2, "--", alpha=0.8, label="Jx residual (np2)")
        ax_res.plot(t2, resid_jy2, "--", alpha=0.8, label="Jy residual (np2)")
    ax_res.axhline(0.0, color="k", lw=1.0, alpha=0.5)
    ax_res.set_xlabel("time")
    ax_res.set_ylabel("Residual")
    ax_res.set_title("(b) Residual to fixed-frequency fit")
    _add_legend(ax_res, loc="upper right", ncol=2, fontsize=9)

    out = fig_dir / "03_langmuir_currents.png"
    fig.tight_layout()
    _save_figure(fig, out, dpi)
    plt.close(fig)
    return out


def _plot_two_stream_weibel(run_dir: Path, fig_dir: Path, dpi: int) -> Path | None:
    ts_dir, _ = _select_case_dir(
        run_dir, ("two_stream_growth_publication", "two_stream_growth_proxy"), "bin"
    )
    wb_dir, _ = _select_case_dir(
        run_dir, ("weibel_growth_publication", "weibel_growth_proxy"), "bin"
    )

    ts_files_np1 = sorted(ts_dir.glob("pic_two_stream_pub_np1.mhd_bcc.*.bin"))
    ts_files_np2 = sorted(ts_dir.glob("pic_two_stream_pub_np2.mhd_bcc.*.bin"))
    ts_publication = bool(ts_files_np1 or ts_files_np2)
    if not ts_publication:
        ts_files_np1 = sorted(
            ts_dir.glob("pic_two_stream_growth_proxy_np1.prtcl_rho.*.bin")
        )
        ts_files_np2 = sorted(
            ts_dir.glob("pic_two_stream_growth_proxy_np2.prtcl_rho.*.bin")
        )

    wb_files_np1 = sorted(wb_dir.glob("pic_weibel_pub_np1.mhd_bcc.*.bin"))
    wb_files_np2 = sorted(wb_dir.glob("pic_weibel_pub_np2.mhd_bcc.*.bin"))
    wb_publication = bool(wb_files_np1 or wb_files_np2)
    if not wb_publication:
        wb_files_np1 = sorted(wb_dir.glob("pic_weibel_growth_proxy_np1.prtcl_jy.*.bin"))
        wb_files_np2 = sorted(wb_dir.glob("pic_weibel_growth_proxy_np2.prtcl_jy.*.bin"))

    if not (ts_files_np1 or ts_files_np2 or wb_files_np1 or wb_files_np2):
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    _style_axes(*axes, grid=True, log_grid=True)

    if ts_files_np1 or ts_files_np2:
        ann = []
        for rank, files, style in (
            (1, ts_files_np1, "-"),
            (2, ts_files_np2, "--"),
        ):
            if not files:
                continue
            t_ts, a_ts = _series_from_bin(
                files,
                (lambda data: np.abs(_mode1_complex(data, "bcc2")))
                if ts_publication
                else (lambda data: np.abs(_mode1_complex(data, "prtcl_rho"))),
            )
            t_plot, a_plot = _trim_startup_zeros(t_ts, a_ts)
            floor = max(1.0e-30, 1.0e-3 * float(np.max(a_plot)))
            a_safe = np.maximum(a_plot, floor)
            axes[0].semilogy(
                t_plot, a_safe, style, marker="o", ms=3, label=f"np{rank}"
            )
            fit = _fit_log_growth(t_plot, a_plot, floor=floor)
            if fit is not None:
                gamma, r2 = fit
                ann.append(f"np{rank}: g={gamma:.2e}, R2={r2:.2f}")
        if ann:
            axes[0].text(
                0.02,
                0.98,
                "\n".join(ann),
                transform=axes[0].transAxes,
                va="top",
                fontsize=10,
                bbox=_annotation_bbox(),
            )
        axes[0].set_title(
            "Two-stream publication growth (k=1)"
            if ts_publication
            else "Two-stream stability proxy (k=1)"
        )
        axes[0].set_xlabel("time")
        axes[0].set_ylabel("|B_y,k=1|" if ts_publication else "|rho_k=1|")
        _add_legend(axes[0], loc="lower right")
    else:
        axes[0].set_title("Two-stream data missing")

    if wb_files_np1 or wb_files_np2:
        ann = []
        for rank, files, style in (
            (1, wb_files_np1, "-"),
            (2, wb_files_np2, "--"),
        ):
            if not files:
                continue
            t_wb, a_wb = _series_from_bin(
                files,
                (lambda data: np.abs(_mode1_complex(data, "bcc2")))
                if wb_publication
                else (lambda data: np.abs(_mode1_complex(data, "prtcl_jy"))),
            )
            t_plot, a_plot = _trim_startup_zeros(t_wb, a_wb)
            floor = max(1.0e-30, 1.0e-3 * float(np.max(a_plot)))
            a_safe = np.maximum(a_plot, floor)
            axes[1].semilogy(
                t_plot, a_safe, style, marker="o", ms=3, label=f"np{rank}"
            )
            fit = _fit_log_growth(t_plot, a_plot, floor=floor)
            if fit is not None:
                gamma, r2 = fit
                ann.append(f"np{rank}: g={gamma:.2e}, R2={r2:.2f}")
        if ann:
            axes[1].text(
                0.02,
                0.98,
                "\n".join(ann),
                transform=axes[1].transAxes,
                va="top",
                fontsize=10,
                bbox=_annotation_bbox(),
            )
        axes[1].set_title(
            "Weibel publication growth (k=1)"
            if wb_publication
            else "Weibel stability proxy (k=1)"
        )
        axes[1].set_xlabel("time")
        axes[1].set_ylabel("|B_y,k=1|" if wb_publication else "|Jy_k=1|")
        _add_legend(axes[1], loc="lower right")
    else:
        axes[1].set_title("Weibel data missing")

    out = fig_dir / "04_two_stream_weibel_modes.png"
    fig.tight_layout()
    _save_figure(fig, out, dpi)
    plt.close(fig)
    return out


def _plot_bell_growth(run_dir: Path, fig_dir: Path, dpi: int) -> Path | None:
    bin_dir, _ = _select_case_dir(
        run_dir, ("bell_growth_publication", "bell_growth_proxy"), "bin"
    )
    pub_tags = [
        ("pic_bell_pub_serial_uncoupled", "serial uncoupled"),
        ("pic_bell_pub_serial_coupled", "serial coupled"),
        ("pic_bell_pub_mpi2_uncoupled", "mpi2 uncoupled"),
        ("pic_bell_pub_mpi2_coupled", "mpi2 coupled"),
    ]
    proxy_tags = [
        ("pic_bell_proxy_serial_uncoupled", "serial uncoupled"),
        ("pic_bell_proxy_serial_coupled", "serial coupled"),
        ("pic_bell_proxy_mpi2_uncoupled", "mpi2 uncoupled"),
        ("pic_bell_proxy_mpi2_coupled", "mpi2 coupled"),
    ]
    tags = (
        pub_tags
        if sorted(bin_dir.glob(pub_tags[0][0] + ".mhd_bcc.*.bin"))
        else proxy_tags
    )
    fig, ax = plt.subplots(figsize=(7, 4))
    _style_axes(ax, grid=True, log_grid=True)
    plotted = 0
    for basename, label in tags:
        files = sorted(bin_dir.glob(f"{basename}.mhd_bcc.*.bin"))
        if not files:
            continue
        times = []
        amps = []
        for fpath in files:
            data = bin_convert.read_binary_as_athdf(str(fpath))
            b2_k = _mode1_complex(data, "bcc2")
            b3_k = _mode1_complex(data, "bcc3")
            amp = np.sqrt(np.abs(b2_k) * np.abs(b2_k) + np.abs(b3_k) * np.abs(b3_k))
            times.append(float(data["Time"]))
            amps.append(float(amp))
        t = np.asarray(times, dtype=float)
        b = np.asarray(amps, dtype=float)
        floor = max(1.0e-30, 1.0e-6 * float(np.max(b)))
        growth = b[-1] / max(b[0], 1.0e-30)
        if "coupled" in basename and "uncoupled" not in basename:
            try:
                fit = fit_exponential_growth_windowed(
                    t, b, min_points=8, floor=floor, min_growth_factor=2.0
                )
            except RuntimeError:
                fit = None
            if fit is not None:
                curve_label = (
                    f"{label} (x{growth:.3f}, gamma={fit['gamma']:.3e}, "
                    f"R2={fit['r2']:.2f})"
                )
            else:
                curve_label = f"{label} (x{growth:.3f})"
        else:
            curve_label = f"{label} (x{growth:.3f})"
        style = "--" if "mpi2" in basename else "-"
        ax.semilogy(
            t, np.maximum(b, floor), style, marker="o", ms=3, label=curve_label
        )
        plotted += 1
    if plotted == 0:
        plt.close(fig)
        return None
    ax.set_xlabel("time")
    ax.set_ylabel(r"$|B_{\perp,k=1}|$")
    ax.set_title("Bell growth: coupled vs uncoupled (k=1 transverse mode)")
    _add_legend(ax, loc="center right")
    out = fig_dir / "05_bell_growth.png"
    fig.tight_layout()
    _save_figure(fig, out, dpi)
    plt.close(fig)
    return out


def _plot_multispecies_osc(run_dir: Path, fig_dir: Path, dpi: int) -> Path | None:
    bin_dir, _ = _select_case_dir(
        run_dir,
        (
            "multispecies_backreaction_publication",
            "multispecies_backreaction_oscillation",
        ),
        "bin",
    )
    pub_tags = ["uniform", "smr", "amr_publication"]
    use_publication = bool(sorted(bin_dir.glob("pic_mso_pub_uniform_np1.mhd_u_m2.*.bin")))
    tags = pub_tags if use_publication else ["uniform", "smr", "amr_proxy"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    _style_axes(axes[0], grid=True)
    _style_axes(axes[1], grid=True, log_grid=True)
    plotted = 0
    for tag in tags:
        prefix = f"pic_mso_pub_{tag}" if use_publication else f"pic_mso_{tag}"
        m2_files = sorted(bin_dir.glob(f"{prefix}_np1.mhd_u_m2.*.bin"))
        e_files = sorted(bin_dir.glob(f"{prefix}_np1.mhd_u_e.*.bin"))
        if not m2_files or not e_files:
            continue

        cycles = {}
        for fpath in m2_files:
            cyc = fpath.name.split(".")[-2]
            cycles.setdefault(cyc, {})["m2"] = fpath
        for fpath in e_files:
            cyc = fpath.name.split(".")[-2]
            cycles.setdefault(cyc, {})["e"] = fpath

        times = []
        m2_vals = []
        e_vals = []
        for cyc in sorted(cycles):
            pair = cycles[cyc]
            if "m2" not in pair or "e" not in pair:
                continue
            m2 = bin_convert.read_binary_as_athdf(str(pair["m2"]))
            ee = bin_convert.read_binary_as_athdf(str(pair["e"]))
            times.append(float(m2["Time"]))
            m2_vals.append(_integrate_quantity(m2, "mom2"))
            e_vals.append(_integrate_quantity(ee, "ener"))
        if not times:
            continue

        t = np.asarray(times, dtype=float)
        m2_arr = np.asarray(m2_vals, dtype=float)
        e_arr = np.asarray(e_vals, dtype=float)
        e_drift = np.abs(e_arr - e_arr[0]) / max(abs(e_arr[0]), 1.0)
        max_drift = float(np.max(e_drift))

        m2_centered = m2_arr - np.mean(m2_arr)
        m2_scale = max(float(np.max(np.abs(m2_centered))), 1.0e-30)
        m2_norm = m2_centered / m2_scale

        disp_tag = "amr" if tag == "amr_publication" else tag
        axes[0].plot(t, m2_norm, label=disp_tag + f" (drift={max_drift:.2e})")
        axes[1].plot(t, e_drift, label=disp_tag)
        plotted += 1

    if plotted == 0:
        plt.close(fig)
        return None

    axes[0].set_title("Multispecies oscillation")
    axes[0].set_xlabel("time")
    axes[0].set_ylabel("Normalized integrated mom2")
    _add_legend(axes[0], loc="upper right")

    axes[1].set_title("Relative energy drift")
    axes[1].set_xlabel("time")
    axes[1].set_ylabel(r"$|E-E_0|/|E_0|$")
    axes[1].set_yscale("log")
    _add_legend(axes[1], loc="lower right")

    out = fig_dir / "06_multispecies_oscillation.png"
    fig.tight_layout()
    _save_figure(fig, out, dpi)
    plt.close(fig)
    return out


def _polarization_split(
    by_mode: np.ndarray,
    bz_mode: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    right = 0.5 * (by_mode - 1.0j * bz_mode)
    left = 0.5 * (by_mode + 1.0j * bz_mode)
    return right, left


def _polarization_series(
    files: list[Path],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    times = []
    by = []
    bz = []
    for fpath in files:
        data = bin_convert.read_binary_as_athdf(str(fpath))
        times.append(float(data["Time"]))
        by.append(_mode1_complex(data, "bcc2"))
        bz.append(_mode1_complex(data, "bcc3"))
    t = np.asarray(times, dtype=float)
    by_mode = np.asarray(by, dtype=complex)
    bz_mode = np.asarray(bz, dtype=complex)
    right, left = _polarization_split(by_mode, bz_mode)
    return t, np.abs(right), np.abs(left)


def _plot_crsi_crpai(run_dir: Path, fig_dir: Path, dpi: int) -> list[Path]:
    outputs: list[Path] = []

    crsi_dir, _ = _select_case_dir(
        run_dir, ("crsi_deltaf_publication", "crsi_deltaf_proxy"), "bin"
    )
    crsi_serial = sorted(crsi_dir.glob("pic_crsi_pub_serial_on.mhd_bcc.*.bin"))
    crsi_mpi2 = sorted(crsi_dir.glob("pic_crsi_pub_mpi2_on.mhd_bcc.*.bin"))
    crsi_mpi4 = sorted(crsi_dir.glob("pic_crsi_pub_mpi4_on.mhd_bcc.*.bin"))
    if not crsi_serial:
        crsi_serial = sorted(crsi_dir.glob("pic_crsi_deltaf_serial_on.mhd_bcc.*.bin"))
        crsi_mpi2 = sorted(crsi_dir.glob("pic_crsi_deltaf_mpi2_on.mhd_bcc.*.bin"))
        crsi_mpi4 = sorted(crsi_dir.glob("pic_crsi_deltaf_mpi4_on.mhd_bcc.*.bin"))

    if crsi_serial:
        t, right, left = _polarization_series(crsi_serial)
        fig, ax = plt.subplots(figsize=(7, 4))
        _style_axes(ax, grid=True, log_grid=True)
        ax.semilogy(t, np.maximum(right, 1.0e-30), label="right np1")
        ax.semilogy(t, np.maximum(left, 1.0e-30), label="left np1")
        if crsi_mpi2:
            t2, r2, l2 = _polarization_series(crsi_mpi2)
            ax.semilogy(t2, np.maximum(r2, 1.0e-30), "--", label="right np2")
            ax.semilogy(t2, np.maximum(l2, 1.0e-30), "--", label="left np2")
        if crsi_mpi4:
            t4, r4, l4 = _polarization_series(crsi_mpi4)
            ax.semilogy(t4, np.maximum(r4, 1.0e-30), ":", label="right np4")
            ax.semilogy(t4, np.maximum(l4, 1.0e-30), ":", label="left np4")
        ax.set_title("CRSI delta-f polarization growth")
        ax.set_xlabel("time")
        ax.set_ylabel("mode amplitude")
        _add_legend(ax, loc="lower right")
        out = fig_dir / "07_crsi_polarization.png"
        fig.tight_layout()
        _save_figure(fig, out, dpi)
        plt.close(fig)
        outputs.append(out)

    crpai_dir, _ = _select_case_dir(
        run_dir, ("crpai_polarization_publication", "crpai_polarization_proxy"), "bin"
    )
    pairs = [
        ("prolate", "pic_crpai_pub_prolate_np", "pic_crpai_prolate_np"),
        ("oblate", "pic_crpai_pub_oblate_np", "pic_crpai_oblate_np"),
    ]
    loaded = []
    for label, pub_prefix, proxy_prefix in pairs:
        rank_series = []
        for rank, style in ((1, "-"), (2, "--"), (4, ":")):
            files = sorted(crpai_dir.glob(f"{pub_prefix}{rank}.mhd_bcc.*.bin"))
            if not files:
                files = sorted(crpai_dir.glob(f"{proxy_prefix}{rank}.mhd_bcc.*.bin"))
            if not files:
                continue
            rank_series.append((rank, style, _polarization_series(files)))
        if not rank_series:
            continue
        loaded.append((label, rank_series))
    if loaded:
        fig, axes = plt.subplots(1, len(loaded), figsize=(6 * len(loaded), 4))
        if len(loaded) == 1:
            axes = [axes]
        _style_axes(*axes, grid=True, log_grid=True)
        for ax, (label, rank_series) in zip(axes, loaded):
            for rank, style, (t, right, left) in rank_series:
                ax.semilogy(
                    t,
                    np.maximum(right, 1.0e-30),
                    style,
                    label=f"right np{rank}",
                )
                ax.semilogy(
                    t,
                    np.maximum(left, 1.0e-30),
                    style,
                    alpha=0.8,
                    label=f"left np{rank}",
                )
            ax.set_title("CRPAI " + label)
            ax.set_xlabel("time")
            ax.set_ylabel("mode amplitude")
            _add_legend(ax, fontsize=9, loc="lower right")
        out = fig_dir / "08_crpai_polarization.png"
        fig.tight_layout()
        _save_figure(fig, out, dpi)
        plt.close(fig)
        outputs.append(out)

    return outputs


def _plot_expanding_anisotropy(run_dir: Path, fig_dir: Path, dpi: int) -> Path | None:
    bin_dir, _ = _select_case_dir(
        run_dir,
        ("expanding_box_anisotropy_publication", "expanding_box_anisotropy_proxy"),
        "bin",
    )
    cases = [
        ("expanding", "pic_box_expanding_np1"),
        ("compressing", "pic_box_compressing_np1"),
    ]
    fig, ax = plt.subplots(figsize=(7, 4))
    _style_axes(ax, grid=True)
    plotted = 0
    for label, base in cases:
        jx_files = sorted(bin_dir.glob(f"{base}.prtcl_jx.*.bin"))
        jy_files = sorted(bin_dir.glob(f"{base}.prtcl_jy.*.bin"))
        if not jx_files or not jy_files:
            continue
        cyc_map = {}
        for fpath in jx_files:
            cyc = fpath.name.split(".")[-2]
            cyc_map.setdefault(cyc, {})["jx"] = fpath
        for fpath in jy_files:
            cyc = fpath.name.split(".")[-2]
            cyc_map.setdefault(cyc, {})["jy"] = fpath

        times = []
        ratio = []
        for cyc in sorted(cyc_map):
            pair = cyc_map[cyc]
            if "jx" not in pair or "jy" not in pair:
                continue
            jx = bin_convert.read_binary_as_athdf(str(pair["jx"]))
            jy = bin_convert.read_binary_as_athdf(str(pair["jy"]))
            jx_int = _integrate_quantity(jx, "prtcl_jx")
            jy_int = _integrate_quantity(jy, "prtcl_jy")
            if abs(jy_int) <= 1.0e-12:
                continue
            times.append(float(jx["Time"]))
            ratio.append(abs(jx_int) / max(abs(jy_int), 1.0e-30))
        if len(times) < 3:
            continue
        t_arr = np.asarray(times, dtype=float)
        r_arr = np.asarray(ratio, dtype=float)
        slope, intercept, slope_stderr = _linear_fit_with_uncertainty(t_arr, r_arr)
        fit = slope * t_arr + intercept
        ax.plot(t_arr, r_arr, "o", ms=3, label=label + " data")
        ax.plot(
            t_arr,
            fit,
            lw=1.8,
            label=label + f" fit (slope={slope:.3f}±{slope_stderr:.3f})",
        )
        plotted += 1

    if plotted == 0:
        plt.close(fig)
        return None
    ax.set_xlabel("time")
    ax.set_ylabel(r"$|J_x| / |J_y|$")
    ax.set_title("Expanding/compressing anisotropy trend")
    _add_legend(ax, loc="upper left")
    out = fig_dir / "09_expanding_box_anisotropy.png"
    fig.tight_layout()
    _save_figure(fig, out, dpi)
    plt.close(fig)
    return out


def _shock_dataset(
    run_dir: Path,
    case_id: str,
    prefix: str,
) -> tuple[dict, dict, dict, dict] | None:
    bin_dir = _case_bin_dir(run_dir, case_id)
    rho_file = _latest_file(f"{prefix}.rho.*.bin", bin_dir)
    bmag_file = _latest_file(f"{prefix}.bmag.*.bin", bin_dir)
    j2_file = _latest_file(f"{prefix}.j2.*.bin", bin_dir)
    if rho_file is None or bmag_file is None or j2_file is None:
        return None
    rho = bin_convert.read_binary_as_athdf(str(rho_file))
    bmag = bin_convert.read_binary_as_athdf(str(bmag_file))
    j2 = bin_convert.read_binary_as_athdf(str(j2_file))

    pvtk_dir = _case_pvtk_dir(run_dir, case_id)
    pfile = _latest_file(f"{prefix}.prtcl_all.*.part.vtk", pvtk_dir)
    pdata = read_particle_vtk(pfile) if pfile is not None else None
    return rho, bmag, j2, pdata


def _plot_shock_story(run_dir: Path, fig_dir: Path, dpi: int) -> Path | None:
    loaded = _shock_dataset(
        run_dir, "amr_shock_publication_local", "pic_amr_shock_pub_local"
    )
    title_suffix = "publication-local"
    if loaded is None:
        # fallback to smoke outputs (limited diagnostics, no particle spectrum)
        bin_dir = _case_bin_dir(run_dir, "amr_shock_lb_smoke")
        bfile = _latest_file("pic_amr_shock_lb_serial.mhd_bcc.*.bin", bin_dir)
        jx = _latest_file("pic_amr_shock_lb_serial.prtcl_jx.*.bin", bin_dir)
        jy = _latest_file("pic_amr_shock_lb_serial.prtcl_jy.*.bin", bin_dir)
        jz = _latest_file("pic_amr_shock_lb_serial.prtcl_jz.*.bin", bin_dir)
        if bfile is None or jx is None or jy is None or jz is None:
            return None
        bdat = bin_convert.read_binary_as_athdf(str(bfile))
        jxdat = bin_convert.read_binary_as_athdf(str(jx))
        jydat = bin_convert.read_binary_as_athdf(str(jy))
        jzdat = bin_convert.read_binary_as_athdf(str(jz))
        rho = {"dens": np.sqrt(np.asarray(bdat["bcc1"], dtype=float) ** 2)}
        bmag = {"bmag": np.sqrt(
            np.asarray(bdat["bcc1"], dtype=float) ** 2 +
            np.asarray(bdat["bcc2"], dtype=float) ** 2 +
            np.asarray(bdat["bcc3"], dtype=float) ** 2
        )}
        j2 = {"j2": (
            np.asarray(jxdat["prtcl_jx"], dtype=float) ** 2 +
            np.asarray(jydat["prtcl_jy"], dtype=float) ** 2 +
            np.asarray(jzdat["prtcl_jz"], dtype=float) ** 2
        )}
        pdata = None
        title_suffix = "smoke-fallback"
    else:
        rho, bmag, j2, pdata = loaded

    rho_2d = np.mean(np.asarray(rho.get("dens", rho.get("rho")), dtype=float), axis=0)
    b_2d = np.mean(np.asarray(bmag.get("bmag"), dtype=float), axis=0)
    j_2d = np.mean(np.sqrt(np.asarray(j2.get("j2"), dtype=float)), axis=0)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    _style_axes(axes, grid=False)
    im0 = axes[0, 0].imshow(rho_2d, origin="lower", aspect="auto", cmap="cividis")
    axes[0, 0].set_title("Density (z-avg)")
    axes[0, 0].set_xlabel("x1 index")
    axes[0, 0].set_ylabel("x2 index")
    cbar0 = fig.colorbar(im0, ax=axes[0, 0], shrink=0.85)
    _style_colorbar(cbar0, "density")

    im1 = axes[0, 1].imshow(b_2d, origin="lower", aspect="auto", cmap="magma")
    axes[0, 1].set_title("|B| (z-avg)")
    axes[0, 1].set_xlabel("x1 index")
    axes[0, 1].set_ylabel("x2 index")
    cbar1 = fig.colorbar(im1, ax=axes[0, 1], shrink=0.85)
    _style_colorbar(cbar1, "|B|")

    im2 = axes[1, 0].imshow(j_2d, origin="lower", aspect="auto", cmap="plasma")
    axes[1, 0].set_title("sqrt(J^2) (z-avg)")
    axes[1, 0].set_xlabel("x1 index")
    axes[1, 0].set_ylabel("x2 index")
    cbar2 = fig.colorbar(im2, ax=axes[1, 0], shrink=0.85)
    _style_colorbar(cbar2, r"$\sqrt{J^2}$")

    ax = axes[1, 1]
    if pdata is not None and "vel" in pdata.vectors:
        speed = np.linalg.norm(pdata.vectors["vel"], axis=1)
        ax.hist(speed, bins=80, density=True, alpha=0.85)
        ax.set_xlabel("speed")
        ax.set_ylabel("PDF")
        ax.set_title("Particle speed spectrum")
        _style_axes(ax, grid=True)
    else:
        ax.text(
            0.5,
            0.5,
            "Particle spectrum unavailable\n(no pvtk velocity dump)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox=_annotation_bbox(),
        )
        ax.set_axis_off()

    fig.suptitle("CR-shock storyline quick-look (" + title_suffix + ")")
    out = fig_dir / "10_shock_storyline.png"
    fig.tight_layout()
    _save_figure(fig, out, dpi)
    plt.close(fig)
    return out


def _run_bundle(bundle: str, run_dir: Path, fig_dir: Path, dpi: int) -> list[Path]:
    outputs: list[Path] = []

    if bundle == "proxy_regression":
        for fn in (
            _plot_entity_deposit_maps,
            _plot_em_vacuum_convergence,
            _plot_langmuir_trace,
            _plot_two_stream_weibel,
            _plot_bell_growth,
            _plot_multispecies_osc,
        ):
            out = fn(run_dir, fig_dir, dpi)
            if out is not None:
                outputs.append(out)

        outputs.extend(_plot_crsi_crpai(run_dir, fig_dir, dpi))

        out = _plot_expanding_anisotropy(run_dir, fig_dir, dpi)
        if out is not None:
            outputs.append(out)
        return outputs

    if bundle == "publication_physics":
        for fn in (
            _plot_em_vacuum_convergence,
            _plot_langmuir_trace,
            _plot_two_stream_weibel,
            _plot_bell_growth,
            _plot_multispecies_osc,
        ):
            out = fn(run_dir, fig_dir, dpi)
            if out is not None:
                outputs.append(out)

        outputs.extend(_plot_crsi_crpai(run_dir, fig_dir, dpi))

        out = _plot_expanding_anisotropy(run_dir, fig_dir, dpi)
        if out is not None:
            outputs.append(out)

        out = _plot_entity_deposit_maps(run_dir, fig_dir, dpi)
        if out is not None:
            outputs.append(out)

        out = _plot_shock_story(run_dir, fig_dir, dpi)
        if out is not None:
            outputs.append(out)
        return outputs

    raise ValueError(f"Unknown figure bundle '{bundle}'")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot publication-priority figures from a run directory."
    )
    parser.add_argument("--run-dir", required=True, help="Path to one run directory.")
    parser.add_argument(
        "--fig-dir",
        default="",
        help="Output directory root (default: <run-dir>/figures).",
    )
    parser.add_argument(
        "--bundle",
        choices=("proxy_regression", "publication_physics", "all"),
        default="all",
        help="Figure bundle selector.",
    )
    parser.add_argument("--dpi", type=int, default=320, help="PNG DPI.")
    args = parser.parse_args()

    _configure_publication_style()

    run_dir = Path(args.run_dir).resolve()
    fig_root = Path(args.fig_dir).resolve() if args.fig_dir else (run_dir / "figures")
    bundles = (
        ["proxy_regression", "publication_physics"]
        if args.bundle == "all"
        else [args.bundle]
    )

    all_outputs: list[Path] = []
    for bundle in bundles:
        fig_dir = fig_root / bundle
        fig_dir.mkdir(parents=True, exist_ok=True)
        all_outputs.extend(_run_bundle(bundle, run_dir, fig_dir, args.dpi))

    if not all_outputs:
        print("No figures generated (missing case artifacts).")
        return 2

    print("Generated figures:")
    for path in all_outputs:
        print(" -", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
