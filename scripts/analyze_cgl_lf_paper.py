#!/usr/bin/env python3
"""Analyze CGL-LF paper-validation histories and reduced diagnostics."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

try:
    import numpy as np
except ImportError:  # pragma: no cover - script still reports a useful error
    np = None


HEADER_RE = re.compile(r"\[(\d+)\]=([^\s]+)")


def parse_hst(path: Path) -> tuple[list[str], "np.ndarray"]:
    if np is None:
        raise RuntimeError("numpy is required to parse Athena history files")
    labels: list[str] | None = None
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                matches = HEADER_RE.findall(stripped)
                if matches:
                    labels = [name for _, name in sorted(matches, key=lambda item: int(item[0]))]
                continue
            rows.append([float(value) for value in stripped.split()])
    if labels is None:
        raise ValueError(f"{path} has no Athena history header")
    data = np.asarray(rows, dtype=float)
    if data.ndim != 2 or data.shape[1] != len(labels):
        raise ValueError(f"{path} has {data.shape[1] if data.ndim == 2 else 0} columns, "
                         f"expected {len(labels)} from header")
    return labels, data


def col(labels: list[str], data: "np.ndarray", name: str) -> "np.ndarray | None":
    return data[:, labels.index(name)] if name in labels else None


def last_scalar(labels: list[str], data: "np.ndarray", name: str,
                *, divide_by: str | None = None) -> float | None:
    values = col(labels, data, name)
    if values is None or len(values) == 0:
        return None
    if divide_by is not None:
        denom = col(labels, data, divide_by)
        if denom is None or denom[-1] == 0.0:
            return None
        return float(values[-1]/denom[-1])
    return float(values[-1])


def summarize_user_history(path: Path) -> dict[str, object]:
    labels, data = parse_hst(path)
    finite = bool(np.isfinite(data).all())
    summary: dict[str, object] = {
        "path": str(path),
        "n_rows": int(data.shape[0]),
        "finite": finite,
        "time_final": last_scalar(labels, data, "time"),
    }
    volume = col(labels, data, "vol")
    if volume is not None and volume[-1] != 0.0:
        for name in ("mass", "ekin", "emag", "ecgl", "beta", "dp", "absdp", "nueff", "hfpow"):
            value = last_scalar(labels, data, name, divide_by="vol")
            if value is not None:
                summary[f"{name}_mean_final"] = value
        for name in ("mir", "fire", "bmir", "bfire", "qpar1", "qperp1", "qpar10", "qperp10"):
            value = last_scalar(labels, data, name, divide_by="vol")
            if value is not None:
                summary[f"{name}_frac_final"] = value
        b2 = col(labels, data, "b2")
        b4 = col(labels, data, "b4")
        if b2 is not None and b4 is not None and b2[-1] != 0.0:
            mean_b2 = b2[-1]/volume[-1]
            mean_b4 = b4[-1]/volume[-1]
            summary["C_B2_final"] = float(mean_b4/(mean_b2*mean_b2) - 1.0)
    for key in ("mir_frac_final", "fire_frac_final", "qpar1_frac_final", "qperp1_frac_final"):
        if key in summary:
            value = float(summary[key])
            summary[f"{key}_in_range"] = bool(-1.0e-12 <= value <= 1.0 + 1.0e-12)
    return summary


def fit_damped_sine(t: "np.ndarray", y: "np.ndarray") -> dict[str, float]:
    y = np.asarray(y, dtype=float)
    signs = np.signbit(y)
    crossings = np.where(signs[1:] != signs[:-1])[0]
    if len(crossings) < 4:
        raise ValueError("need at least four zero crossings for frequency estimate")
    zero_times = []
    for idx in crossings:
        t0, t1 = t[idx], t[idx + 1]
        y0, y1 = y[idx], y[idx + 1]
        zero_times.append(float(t0 - y0*(t1 - t0)/(y1 - y0)))
    half_period = float(np.mean(np.diff(zero_times)))
    omega = math.pi/half_period

    abs_y = np.abs(y)
    peak_idx = []
    for idx in range(1, len(abs_y) - 1):
        if abs_y[idx] >= abs_y[idx - 1] and abs_y[idx] >= abs_y[idx + 1]:
            peak_idx.append(idx)
    if len(peak_idx) < 3:
        raise ValueError("need at least three peaks for damping estimate")
    peak_t = t[peak_idx]
    peak_a = np.maximum(abs_y[peak_idx], 1.0e-300)
    slope, intercept = np.polyfit(peak_t, np.log(peak_a), 1)
    return {"omega": float(omega), "gamma": float(-slope), "amplitude0": float(math.exp(intercept))}


def run_synthetic_wave_test() -> dict[str, object]:
    if np is None:
        raise RuntimeError("numpy is required for synthetic wave fitting")
    omega_true = 2.4
    gamma_true = 0.17
    t = np.linspace(0.0, 40.0, 4001)
    y = 0.6*np.exp(-gamma_true*t)*np.sin(omega_true*t + 0.35)
    fit = fit_damped_sine(t, y)
    fit["omega_true"] = omega_true
    fit["gamma_true"] = gamma_true
    fit["omega_rel_err"] = abs(fit["omega"] - omega_true)/omega_true
    fit["gamma_rel_err"] = abs(fit["gamma"] - gamma_true)/gamma_true
    fit["passed"] = bool(fit["omega_rel_err"] < 2.0e-3 and fit["gamma_rel_err"] < 2.0e-2)
    return fit


def plot_histories(summaries: list[dict[str, object]], output_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    for summary in summaries:
        path = Path(str(summary["path"]))
        labels, data = parse_hst(path)
        time = col(labels, data, "time")
        vol = col(labels, data, "vol")
        if time is None or vol is None:
            continue
        fig, ax = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
        for name in ("mir", "fire", "bmir", "bfire"):
            values = col(labels, data, name)
            if values is not None:
                ax[0].plot(time, values/vol, label=name)
        ax[0].set_ylabel("volume fraction")
        ax[0].legend(loc="best", fontsize=8)
        for name in ("qpar1", "qperp1", "qpar10", "qperp10"):
            values = col(labels, data, name)
            if values is not None:
                ax[1].plot(time, values/vol, label=name)
        ax[1].set_xlabel("time")
        ax[1].set_ylabel("cap activity")
        ax[1].legend(loc="best", fontsize=8)
        fig.tight_layout()
        fig.savefig(output_dir / f"{path.stem}_history.png", dpi=160)
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path,
                        help="History files or directories containing *.user.hst files")
    parser.add_argument("--output-dir", type=Path, default=Path("cgl_lf_paper_analysis"))
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--make-plots", action="store_true")
    parser.add_argument("--synthetic-wave-test", action="store_true")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    history_files: list[Path] = []
    for item in args.inputs:
        if item.is_dir():
            history_files.extend(sorted(item.glob("*.user.hst")))
        else:
            history_files.append(item)

    summaries = [summarize_user_history(path) for path in history_files]
    result: dict[str, object] = {
        "histories": summaries,
        "all_histories_finite": all(bool(item.get("finite", False)) for item in summaries),
    }
    if args.synthetic_wave_test:
        result["synthetic_wave_fit"] = run_synthetic_wave_test()
    if args.make_plots:
        plot_histories(summaries, args.output_dir)

    summary_path = args.summary_json or (args.output_dir / "summary.json")
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
    print(f"Wrote {summary_path}")
    if summaries:
        for item in summaries:
            print(f"{Path(str(item['path'])).name}: finite={item['finite']} "
                  f"t_final={item.get('time_final')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
