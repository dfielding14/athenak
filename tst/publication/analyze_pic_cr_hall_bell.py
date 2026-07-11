#!/usr/bin/env python3
"""Compact raw-output analysis for the full CR-Hall Bell validation."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[2]
K_BELL = 2.0 * math.pi
THEORY = {
    "off": {"k": K_BELL, "growth": 6.28287114006,
            "frequency": 0.0628318530718},
    "full": {"k": 0.8 * K_BELL, "growth": 5.55762565556,
             "frequency": 2.57547765741},
}


def _reader():
    path = REPO / "vis/python/bin_convert_new.py"
    spec = importlib.util.spec_from_file_location("pic_hall_bin", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("unable to load Athena binary reader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.read_binary_as_athdf


def _files(run: Path, token: str) -> list[Path]:
    paths = sorted(run.glob(f"*.{token}.*.bin"))
    if not paths:
        paths = sorted((run / "bin").glob(f"*.{token}.*.bin"))
    if not paths:
        raise RuntimeError(f"no {token} binary outputs in {run}")
    return paths


def _profile(dataset: dict, name: str) -> np.ndarray:
    values = np.asarray(dataset[name], dtype=float)
    return np.mean(values.reshape((-1, values.shape[-1])), axis=0)


def _mode(dataset: dict, wavenumber: float) -> tuple[complex, complex]:
    x = np.asarray(dataset["x1v"], dtype=float)
    by = _profile(dataset, "bcc2")
    bz = _profile(dataset, "bcc3")
    phase = np.exp(-1j * wavenumber * x)
    unstable = np.mean((by + 1j * bz) * phase)
    stable = np.mean((by - 1j * bz) * phase)
    return complex(unstable), complex(stable)


def _dominant_wavelength(dataset: dict) -> float:
    x = np.asarray(dataset["x1v"], dtype=float)
    transverse = _profile(dataset, "bcc2") + 1j * _profile(dataset, "bcc3")
    spectrum = np.abs(np.fft.fft(transverse - np.mean(transverse)))
    frequencies = np.fft.fftfreq(x.size, d=float(np.mean(np.diff(x))))
    spectrum[frequencies == 0.0] = 0.0
    index = int(np.argmax(spectrum))
    frequency = abs(float(frequencies[index]))
    if frequency == 0.0 or spectrum[index] == 0.0:
        raise RuntimeError("transverse magnetic field has no nonzero Fourier mode")
    return 1.0 / frequency


def _fit_trace(times: np.ndarray, mode: np.ndarray) -> dict:
    mask = (times >= 0.3) & (times <= 1.1)
    if np.count_nonzero(mask) < 8:
        raise RuntimeError("Bell fit requires at least eight snapshots in [0.3, 1.1]")
    growth_fit = np.polyfit(times[mask], np.log(np.abs(mode[mask])), 1)
    phase_fit = np.polyfit(times[mask], np.unwrap(np.angle(mode[mask])), 1)
    fitted = np.polyval(growth_fit, times[mask])
    observed = np.log(np.abs(mode[mask]))
    residual = observed - fitted
    variance = np.sum((observed - np.mean(observed)) ** 2)
    r2 = 1.0 if variance == 0.0 else 1.0 - np.sum(residual**2) / variance
    return {
        "growth": float(growth_fit[0]),
        "frequency": abs(float(phase_fit[0])),
        "growth_r2": float(r2),
        "mask": mask,
    }


def _volume_mean(dataset: dict, name: str) -> float:
    dx1 = np.diff(np.asarray(dataset["x1f"], dtype=float))
    dx2 = np.diff(np.asarray(dataset["x2f"], dtype=float))
    dx3 = np.diff(np.asarray(dataset["x3f"], dtype=float))
    volume = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    return float(np.sum(np.asarray(dataset[name], dtype=float) * volume) /
                 np.sum(volume))


def _history(run: Path) -> tuple[np.ndarray, np.ndarray]:
    paths = list(run.glob("*.mhd.hst"))
    if not paths:
        paths = list((run / "bin").glob("*.mhd.hst"))
    if len(paths) != 1:
        raise RuntimeError(f"expected one MHD history in {run}")
    data = np.atleast_2d(np.loadtxt(paths[0]))
    return data[:, 0], data[:, -2:]


def analyze_run(run: Path, model: str) -> dict:
    read = _reader()
    datasets = [read(str(path)) for path in _files(run, "mhd_w_bcc")]
    datasets.sort(key=lambda item: float(item["Time"]))
    times = np.asarray([float(item["Time"]) for item in datasets])
    modes = np.asarray([_mode(item, THEORY[model]["k"])[0]
                        for item in datasets])
    opposite = np.asarray([_mode(item, THEORY[model]["k"])[1]
                           for item in datasets])
    fit = _fit_trace(times, modes)
    fit_mask = fit.pop("mask")
    expected = THEORY[model]
    result = {
        "run": str(run.resolve()),
        "model": model,
        "times": times.tolist(),
        "mode_real": modes.real.tolist(),
        "mode_imag": modes.imag.tolist(),
        "amplitude": np.abs(modes).tolist(),
        "phase": np.unwrap(np.angle(modes)).tolist(),
        "fit_mask": fit_mask.tolist(),
        **fit,
        "expected_growth": expected["growth"],
        "expected_frequency": expected["frequency"],
        "growth_relative_error": abs(fit["growth"] - expected["growth"])
        / expected["growth"],
        "frequency_relative_error": abs(fit["frequency"] - expected["frequency"])
        / expected["frequency"],
        "polarization_ratio": float(np.median(
            np.abs(modes[fit_mask]) /
            np.maximum(np.abs(opposite[fit_mask]), np.finfo(float).tiny))),
        "dominant_wavelength": _dominant_wavelength(datasets[-1]),
    }
    if model == "full":
        htime, hall = _history(run)
        early = (htime > 0.0) & (htime <= 0.3)
        result["hall_history_time"] = htime.tolist()
        result["hall_rmax"] = hall[:, 0].tolist()
        result["hall_lambda_max"] = hall[:, 1].tolist()
        result["early_hall_rmax"] = float(np.median(hall[early, 0]))
        result["early_hall_lambda_max"] = float(np.median(hall[early, 1]))

        # The final raw snapshot is guaranteed to follow an actual moment deposit;
        # initialization-time outputs may precede the first staged deposit.
        rho = read(str(_files(run, "prtcl_rho")[-1]))
        jx = read(str(_files(run, "prtcl_jx")[-1]))
        jy = read(str(_files(run, "prtcl_jy")[-1]))
        jz = read(str(_files(run, "prtcl_jz")[-1]))
        result["raw_mean_qcr_over_c"] = _volume_mean(rho, "prtcl_rho")
        result["raw_mean_jcr_over_c"] = [
            _volume_mean(jx, "prtcl_jx"),
            _volume_mean(jy, "prtcl_jy"),
            _volume_mean(jz, "prtcl_jz"),
        ]
    return result


def _plot(results: list[dict], output: Path) -> None:
    colors = {"off": "#315a9b", "full": "#c33d32", "full-hi": "#e08b22"}
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)
    for label, result in zip(("off", "full", "full-hi"), results):
        time = np.asarray(result["times"])
        amplitude = np.asarray(result["amplitude"])
        phase = np.asarray(result["phase"])
        axes[0, 0].semilogy(time, amplitude, color=colors[label], label=label)
        axes[0, 1].plot(time, phase, color=colors[label], label=label)
        axes[1, 0].scatter(result["dominant_wavelength"], result["growth"],
                           color=colors[label], s=45, label=label)
    full = results[1]
    axes[1, 1].plot(full["hall_history_time"], full["hall_rmax"],
                    color="#4c78a8", label=r"$\max|R|$")
    axes[1, 1].plot(full["hall_history_time"], full["hall_lambda_max"],
                    color="#d95f02", label=r"$\max\Lambda$")
    axes[0, 0].set(xlabel="time", ylabel="Bell-unstable amplitude")
    axes[0, 1].set(xlabel="time", ylabel="unwrapped phase [rad]")
    axes[1, 0].set(xlabel="dominant wavelength", ylabel="growth rate")
    axes[1, 1].set(xlabel="time", ylabel="Hall strength")
    for axis in axes.flat:
        axis.grid(alpha=0.2)
        axis.legend(frameon=False)
    fig.suptitle("Full CR-Hall Bell validation")
    fig.savefig(output, dpi=240)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--off", type=Path, required=True)
    parser.add_argument("--full", type=Path, required=True)
    parser.add_argument("--full-hi", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    off = analyze_run(args.off, "off")
    full = analyze_run(args.full, "full")
    high = analyze_run(args.full_hi, "full")
    resolution_growth = (
        abs(high["growth"] - full["growth"])
        / max(abs(full["growth"]), np.finfo(float).tiny)
    )
    resolution_frequency = (
        abs(high["frequency"] - full["frequency"])
        / max(abs(full["frequency"]), np.finfo(float).tiny)
    )
    expected_q = 4.0 * math.pi / 100.0
    expected_j = 4.0 * math.pi
    checks = {
        "off_growth": off["growth_relative_error"] < 0.05,
        "off_frequency": off["frequency_relative_error"] < 0.05,
        "full_growth": full["growth_relative_error"] < 0.05,
        "full_frequency": full["frequency_relative_error"] < 0.05,
        "full_polarization": full["polarization_ratio"] > 10.0,
        "full_wavelength": abs(full["dominant_wavelength"] - 1.25) < 0.03,
        "resolution_growth": resolution_growth < 0.03,
        "resolution_frequency": resolution_frequency < 0.03,
        "raw_charge": abs(full["raw_mean_qcr_over_c"] - expected_q) / expected_q
        < 0.01,
        "raw_current": abs(full["raw_mean_jcr_over_c"][0] - expected_j) / expected_j
        < 0.01,
        "hall_R": abs(full["early_hall_rmax"] - 0.01) < 0.002,
        "hall_Lambda": abs(full["early_hall_lambda_max"] - 1.0) < 0.05,
    }
    report = {
        "off": off,
        "full": full,
        "full_high_resolution": high,
        "resolution_growth_difference": resolution_growth,
        "resolution_frequency_difference": resolution_frequency,
        "checks": checks,
        "passed": all(checks.values()),
    }
    (args.output_dir / "pic_cr_hall_bell_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _plot([off, full, high], args.output_dir / "pic_cr_hall_bell_validation.png")
    print(json.dumps({"checks": checks, "passed": report["passed"]}, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
