#!/usr/bin/env python3
"""Plot box-integrated energy diagnostics for the beta25 turbulence comparison."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


HEADER_RE = re.compile(r"\[(\d+)\]=(\S+)")


def read_history(path: Path) -> dict[str, np.ndarray]:
    header = ""
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") and "[1]=" in line:
                header = line
                break
    if not header:
        raise ValueError(f"Missing Athena history column header in {path}")
    indices = {name: int(index) - 1 for index, name in HEADER_RE.findall(header)}
    required = {"time", "mass", "tot-E", "1-KE", "2-KE", "3-KE", "1-ME", "2-ME", "3-ME"}
    missing = sorted(required - indices.keys())
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    values = np.loadtxt(path, comments="#", ndmin=2)
    if values.size == 0:
        raise ValueError(f"No history records in {path}")
    return {name: values[:, index] for name, index in indices.items()}


def diagnostics(history: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    mass = history["mass"]
    kinetic = history["1-KE"] + history["2-KE"] + history["3-KE"]
    magnetic = history["1-ME"] + history["2-ME"] + history["3-ME"]
    return {
        "time": history["time"],
        "total_energy": history["tot-E"],
        "kinetic_energy": kinetic,
        "magnetic_energy": magnetic,
        "velocity_rms": np.sqrt(2.0 * kinetic / mass),
        "alfven_rms": np.sqrt(2.0 * magnetic / mass),
    }


def save_single_plot(
    runs: list[tuple[str, dict[str, np.ndarray]]],
    key: str,
    ylabel: str,
    output: Path,
) -> None:
    fig, axis = plt.subplots(figsize=(6.6, 3.8))
    for label, values in runs:
        axis.plot(values["time"], values[key], label=label, linewidth=1.6)
    axis.set_xlabel("Time")
    axis.set_ylabel(ylabel)
    axis.grid(alpha=0.3)
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_csv(runs: list[tuple[str, dict[str, np.ndarray]]], output: Path) -> None:
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["model", "time", "total_energy", "kinetic_energy", "magnetic_energy",
             "velocity_rms", "alfven_rms"]
        )
        for label, values in runs:
            for row in zip(
                values["time"],
                values["total_energy"],
                values["kinetic_energy"],
                values["magnetic_energy"],
                values["velocity_rms"],
                values["alfven_rms"],
            ):
                writer.writerow([label, *(f"{item:.16e}" for item in row)])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mhd-history", type=Path, required=True)
    parser.add_argument("--cgl-history", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    runs = [
        ("MHD", diagnostics(read_history(args.mhd_history))),
        ("MHD-CGL", diagnostics(read_history(args.cgl_history))),
    ]

    save_single_plot(runs, "total_energy", r"Box-integrated total energy $E_{\rm tot}$",
                     args.output_dir / "total_energy_over_time.png")
    save_single_plot(runs, "velocity_rms", r"$v_{\rm rms} = (2 E_{\rm kin}/M)^{1/2}$",
                     args.output_dir / "velocity_rms_over_time.png")
    save_single_plot(runs, "alfven_rms",
                     r"$v_{\rm A,rms} = (2 E_{\rm mag}/M)^{1/2}$",
                     args.output_dir / "alfven_rms_over_time.png")

    fig, axes = plt.subplots(3, 1, figsize=(7.0, 8.6), sharex=True)
    panels = [
        ("total_energy", r"$E_{\rm tot}$"),
        ("velocity_rms", r"$v_{\rm rms}$"),
        ("alfven_rms", r"$v_{\rm A,rms}$"),
    ]
    for axis, (key, ylabel) in zip(axes, panels):
        for label, values in runs:
            axis.plot(values["time"], values[key], label=label, linewidth=1.6)
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.3)
    axes[0].legend(frameon=False)
    axes[-1].set_xlabel("Time")
    fig.tight_layout()
    fig.savefig(args.output_dir / "energy_velocity_comparison.png", dpi=180)
    plt.close(fig)

    write_csv(runs, args.output_dir / "energy_velocity_comparison.csv")
    for label, values in runs:
        print(
            f"{label}: rows={len(values['time'])} t_final={values['time'][-1]:.8g} "
            f"E_final={values['total_energy'][-1]:.8g} "
            f"v_rms_final={values['velocity_rms'][-1]:.8g} "
            f"v_A_rms_final={values['alfven_rms'][-1]:.8g}"
        )


if __name__ == "__main__":
    main()
