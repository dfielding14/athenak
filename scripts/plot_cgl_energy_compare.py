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
from matplotlib.lines import Line2D  # noqa: E402


HEADER_RE = re.compile(r"\[(\d+)\]=(\S+)")
RUN_COLORS = {"MHD": "tab:blue", "Pure CGL": "tab:orange"}
ENERGY_LINES = [
    ("total_energy", r"$E_{\rm tot}$", "-"),
    ("kinetic_energy", r"$E_{\rm kin}$", "--"),
    ("magnetic_energy", r"$E_{\rm mag}$", "-."),
    ("thermal_energy", r"$E_{\rm th}$", ":"),
]
ENERGY_CHANGE_LINES = [
    ("total_energy_change", r"$\Delta E_{\rm tot}$", "-"),
    ("kinetic_energy_change", r"$\Delta E_{\rm kin}$", "--"),
    ("magnetic_energy_change", r"$\Delta E_{\rm mag}$", "-."),
    ("thermal_energy_change", r"$\Delta E_{\rm th}$", ":"),
]


def read_history(path: Path) -> dict[str, np.ndarray]:
    header = ""
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") and "[1]=" in line:
                if header and HEADER_RE.findall(header) != HEADER_RE.findall(line):
                    raise ValueError(f"History columns change within {path}")
                header = line
    if not header:
        raise ValueError(f"Missing Athena history column header in {path}")
    indices = {name: int(index) - 1 for index, name in HEADER_RE.findall(header)}
    values = np.loadtxt(path, comments="#", ndmin=2)
    if values.size == 0:
        raise ValueError(f"No history records in {path}")
    if values.shape[1] != len(indices) or "time" not in indices:
        raise ValueError(f"History columns do not match header in {path}")
    times = values[:, indices["time"]]
    if not np.all(np.isfinite(times)):
        raise ValueError(f"Nonfinite history times in {path}")
    # A resumed branch replaces the abandoned rows at and after its first time.
    retained = []
    for row, time in enumerate(times):
        while retained and times[retained[-1]] >= time:
            retained.pop()
        retained.append(row)
    values = values[retained]
    return {name: values[:, index] for name, index in indices.items()}


def diagnostics(
    history: dict[str, np.ndarray], user_history: dict[str, np.ndarray] | None = None,
) -> dict[str, np.ndarray]:
    required = {"time", "mass", "tot-E", "1-KE", "2-KE", "3-KE"}
    missing = sorted(required - history.keys())
    if missing:
        raise ValueError(
            f"Energy history is missing required columns: {', '.join(missing)}")
    physical = history
    if user_history is not None:
        indices = np.searchsorted(user_history["time"], history["time"])
        if (np.any(indices >= len(user_history["time"])) or
                not np.array_equal(user_history["time"][indices], history["time"])):
            raise ValueError("MHD and user histories must contain matching sample times")
        physical = {key: values[indices] for key, values in user_history.items()}
    mass = history["mass"]
    kinetic = history["1-KE"] + history["2-KE"] + history["3-KE"]
    # Passive J/A heat is separate from the MHD total energy in the gas reservoir.
    thermal_key = "therm_mhd" if "cgl-J" in history else "therm_cgl"
    if "B^2" in physical:
        magnetic = 0.5*physical["B^2"]
    elif "magnetic" in physical:
        magnetic = physical["magnetic"]
    elif thermal_key in physical:
        magnetic = history["tot-E"] - kinetic - physical[thermal_key]
    else:
        raise ValueError("Exact thermal energy requires a companion .user.hst with "
                         "therm_cgl, therm_mhd, B^2, or magnetic; recover missing "
                         "thermal energy from snapshots")
    thermal = (physical[thermal_key] if thermal_key in physical else
               history["tot-E"] - kinetic - magnetic)
    return {
        "time": history["time"],
        "total_energy": history["tot-E"],
        "total_energy_change": history["tot-E"] - history["tot-E"][0],
        "kinetic_energy": kinetic,
        "kinetic_energy_change": kinetic - kinetic[0],
        "magnetic_energy": magnetic,
        "magnetic_energy_change": magnetic - magnetic[0],
        "thermal_energy": thermal,
        "thermal_energy_change": thermal - thermal[0],
        "velocity_rms": np.sqrt(2.0 * kinetic / mass),
        "alfven_rms": np.sqrt(2.0 * magnetic / mass),
    }


def plot_energy_curves(
    axis: plt.Axes,
    runs: list[tuple[str, dict[str, np.ndarray]]],
    curves: list[tuple[str, str, str]],
) -> None:
    for label, values in runs:
        for key, _, linestyle in curves:
            linewidth = 1.8 if key in {"total_energy", "total_energy_change"} else 1.4
            axis.plot(
                values["time"],
                values[key],
                color=RUN_COLORS[label],
                linestyle=linestyle,
                linewidth=linewidth,
            )
    model_handles = [
        Line2D([0], [0], color=RUN_COLORS[label], linewidth=1.8, label=label)
        for label, _ in runs
    ]
    component_handles = [
        Line2D([0], [0], color="black", linestyle=linestyle, linewidth=1.6, label=label)
        for _, label, linestyle in curves
    ]
    model_legend = axis.legend(handles=model_handles, loc="upper left", frameon=False)
    axis.add_artist(model_legend)
    axis.legend(handles=component_handles, loc="upper right", frameon=False)


def save_energy_plot(
    runs: list[tuple[str, dict[str, np.ndarray]]],
    output: Path,
) -> None:
    fig, axis = plt.subplots(figsize=(6.6, 3.8))
    plot_energy_curves(axis, runs, ENERGY_LINES)
    axis.set_xlabel("Time")
    axis.set_ylabel(r"Box-integrated energy")
    axis.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def save_energy_change_plot(
    runs: list[tuple[str, dict[str, np.ndarray]]],
    output: Path,
) -> None:
    fig, axis = plt.subplots(figsize=(6.6, 3.8))
    plot_energy_curves(axis, runs, ENERGY_CHANGE_LINES)
    axis.set_xlabel("Time")
    axis.set_ylabel(r"Change in box-integrated energy")
    axis.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


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
             "thermal_energy", "total_energy_change", "kinetic_energy_change",
             "magnetic_energy_change", "thermal_energy_change", "velocity_rms",
             "alfven_rms"]
        )
        for label, values in runs:
            for row in zip(
                values["time"],
                values["total_energy"],
                values["kinetic_energy"],
                values["magnetic_energy"],
                values["thermal_energy"],
                values["total_energy_change"],
                values["kinetic_energy_change"],
                values["magnetic_energy_change"],
                values["thermal_energy_change"],
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

    runs = []
    for label, path in (("MHD", args.mhd_history), ("Pure CGL", args.cgl_history)):
        user_path = path.with_name(path.name.removesuffix(".mhd.hst") + ".user.hst")
        user = read_history(user_path) if user_path.is_file() else None
        runs.append((label, diagnostics(read_history(path), user)))
    args.output_dir.mkdir(parents=True, exist_ok=True)

    save_energy_plot(runs, args.output_dir / "total_energy_over_time.png")
    save_energy_change_plot(runs, args.output_dir / "energy_change_over_time.png")
    save_single_plot(runs, "velocity_rms", r"$v_{\rm rms} = (2 E_{\rm kin}/M)^{1/2}$",
                     args.output_dir / "velocity_rms_over_time.png")
    save_single_plot(runs, "alfven_rms",
                     r"$v_{\rm A,rms} = (2 E_{\rm mag}/M)^{1/2}$",
                     args.output_dir / "alfven_rms_over_time.png")

    fig, axes = plt.subplots(3, 1, figsize=(7.0, 8.6), sharex=True)
    panels = [
        ("total_energy_change", r"$\Delta E_{\rm tot}$"),
        ("velocity_rms", r"$v_{\rm rms}$"),
        ("alfven_rms", r"$v_{\rm A,rms}$"),
    ]
    for axis, (key, ylabel) in zip(axes, panels):
        if key == "total_energy_change":
            plot_energy_curves(axis, runs, ENERGY_CHANGE_LINES)
        else:
            for label, values in runs:
                axis.plot(
                    values["time"],
                    values[key],
                    color=RUN_COLORS[label],
                    label=label,
                    linewidth=1.6,
                )
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.3)
    axes[-1].set_xlabel("Time")
    fig.tight_layout()
    fig.savefig(args.output_dir / "energy_velocity_comparison.png", dpi=180)
    plt.close(fig)

    write_csv(runs, args.output_dir / "energy_velocity_comparison.csv")
    for label, values in runs:
        print(
            f"{label}: rows={len(values['time'])} t_final={values['time'][-1]:.8g} "
            f"E_final={values['total_energy'][-1]:.8g} "
            f"delta_E_final={values['total_energy_change'][-1]:.8g} "
            f"v_rms_final={values['velocity_rms'][-1]:.8g} "
            f"v_A_rms_final={values['alfven_rms'][-1]:.8g}"
        )


if __name__ == "__main__":
    main()
