#!/usr/bin/env python3
"""Analyze one current-HEAD 128^3 full-CR-Hall turbulent-dynamo pilot."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
import re
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[2]
HISTORY_LABEL = re.compile(r"\[[0-9]+\]=([^\s]+)")
EVIDENCE_NOTE = "engineering_proxy | not_sun_bai_reproduction=true"
COLORS = ("#3b4cc0", "#2a9d8f", "#d1495b")
ENERGY_RESIDUAL_LIMIT = 1.0e-2
TOTAL_MOMENTUM_LIMIT = 1.0e-3
DIVB_LIMIT = 1.0e-10


def _binary_reader():
    path = REPO / "vis/python/bin_convert_new.py"
    spec = importlib.util.spec_from_file_location("pic_turb_bin", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("unable to load Athena binary reader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.read_binary_as_athdf


def _binary_header(path: Path) -> dict[str, Any]:
    """Read only the small ASCII preamble and embedded input deck."""
    with path.open("rb") as stream:
        if stream.readline().strip() != b"Athena binary output version=1.1":
            raise RuntimeError(f"{path}: unsupported Athena binary output")
        pheader_lines = int(stream.readline().split(b"=")[-1])
        preheader: dict[str, str] = {}
        for _ in range(pheader_lines - 1):
            key, value = stream.readline().decode("ascii").split("=", 1)
            preheader[key.strip()] = value.strip()
        nvars = int(stream.readline().split(b"=")[-1])
        variables = stream.readline().decode("ascii").split()[1:]
        if len(variables) != nvars:
            raise RuntimeError(f"{path}: inconsistent variable inventory")
        header_size = int(stream.readline().split(b"=")[-1])
        payload = stream.read(header_size)
        if len(payload) != header_size:
            raise RuntimeError(f"{path}: truncated input header")
    return {
        "time": float(preheader["time"]),
        "cycle": int(preheader["cycle"]),
        "variables": variables,
        "input": payload.decode("ascii").splitlines(),
    }


def _parameters(lines: list[str]) -> dict[str, dict[str, str]]:
    blocks: dict[str, dict[str, str]] = {}
    active: dict[str, str] | None = None
    for raw in lines:
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            name = line[1:-1]
            if name == "par_end":
                break
            active = blocks.setdefault(name, {})
        elif active is not None and "=" in line:
            name, value = line.split("=", 1)
            active[name.strip()] = value.strip()
    return blocks


def _history(path: Path) -> dict[str, np.ndarray]:
    labels: tuple[str, ...] | None = None
    rows: list[list[float]] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            found = tuple(HISTORY_LABEL.findall(line))
            if found:
                if labels is not None and found != labels:
                    raise RuntimeError(f"{path}: history labels changed")
                labels = found
            continue
        if labels is None:
            continue
        values = [float(value) for value in line.split()]
        if len(values) != len(labels):
            raise RuntimeError(f"{path}: history row width changed")
        rows.append(values)
    if labels is None or len(rows) < 3:
        raise RuntimeError(f"{path}: incomplete history")
    data = np.asarray(rows, dtype=np.float64)
    if not np.all(np.isfinite(data)) or np.any(np.diff(data[:, 0]) <= 0.0):
        raise RuntimeError(f"{path}: non-finite or non-monotone history")
    return {name: data[:, index] for index, name in enumerate(labels)}


def _one_file(run: Path, suffix: str) -> Path:
    paths = sorted(run.glob(f"*.{suffix}.hst"))
    if len(paths) != 1:
        raise RuntimeError(f"{run}: expected one *.{suffix}.hst")
    return paths[0]


def _real(
    parameters: dict[str, dict[str, str]], block: str, name: str
) -> float:
    try:
        return float(parameters[block][name])
    except KeyError as exc:
        raise RuntimeError(f"missing <{block}>/{name}") from exc


def _validate_model(
    parameters: dict[str, dict[str, str]], expected_hall_mode: str
) -> dict[str, Any]:
    particles = parameters.get("particles", {})
    if expected_hall_mode not in {"full", "off"}:
        raise RuntimeError(f"unsupported Hall mode: {expected_hall_mode}")
    if particles.get("pic_cr_hall_mode") != expected_hall_mode:
        raise RuntimeError(
            "analysis expected particles/pic_cr_hall_mode=" + expected_hall_mode
        )
    required_controls = {
        ("time", "integrator"): "vl2",
        ("particles", "particle_type"): "cosmic_ray",
        ("particles", "pusher"): "boris_tsc",
        ("particles", "deposit_moments"): "true",
        ("particles", "deposit_order"): "2",
        ("particles", "couple_moments_to_mhd"): "true",
        ("particles", "couple_moments_momentum_to_mhd"): "true",
        ("particles", "couple_moments_energy_to_mhd"): "true",
        ("particles", "pic_background_mode"): "coupled",
        ("particles", "pic_feedback_mode"): "coupled",
        ("particles", "pic_interp_scheme"): "tsc",
    }
    for (block, name), expected in required_controls.items():
        measured = parameters.get(block, {}).get(name)
        if measured != expected:
            raise RuntimeError(
                f"analysis requires {block}/{name}={expected}; measured {measured!r}"
            )
    if parameters.get("mesh_refinement", {}).get("refinement") != "none":
        raise RuntimeError("analysis requires the uniform 128^3 box")
    shape = tuple(int(parameters["mesh"][f"nx{axis}"]) for axis in (1, 2, 3))
    if shape != (128, 128, 128):
        raise RuntimeError(f"unexpected mesh: {shape}")
    alpha_i = _real(parameters, "particles", "pic_background_ion_q_over_mc")
    nspecies = int(particles["nspecies"])
    qom = []
    for species in range(nspecies):
        block = parameters[f"species{species}"]
        qom.append(float(block["charge"]) / float(block["mass"]))
    if not all(math.isclose(value, alpha_i, rel_tol=0.0, abs_tol=1.0e-12)
               for value in qom):
        raise RuntimeError("background-ion and CR species q/(mc) do not match")
    return {
        "mesh": list(shape),
        "hall_mode": expected_hall_mode,
        "integrator": parameters["time"]["integrator"],
        "coupling": "vl2_tsc",
        "background_ion_q_over_mc": alpha_i,
        "cr_species_q_over_mc": qom,
        "ppc": float(particles["ppc"]),
        "artificial_cr_light_speed": float(particles["pic_cr_light_speed"]),
        "tlim": _real(parameters, "time", "tlim"),
        "driving_energy_rate": _real(parameters, "turb_driving", "dedt"),
    }


def _spectrum(components: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    shape = components[0].shape
    if any(component.shape != shape for component in components):
        raise RuntimeError("spectrum components have inconsistent shapes")
    nz, ny, nx = shape
    kz = np.fft.fftfreq(nz) * nz
    ky = np.fft.fftfreq(ny) * ny
    kx = np.fft.rfftfreq(nx) * nx
    shell = np.rint(np.sqrt(kz[:, None, None] ** 2 +
                            ky[None, :, None] ** 2 +
                            kx[None, None, :] ** 2)).astype(np.int32)
    weight = np.ones(kx.size, dtype=np.float64)
    if kx.size > 1:
        weight[1:-1] = 2.0
        if nx % 2:
            weight[-1] = 2.0
    energy = np.zeros_like(shell, dtype=np.float64)
    norm = float(nx * ny * nz)
    for component in components:
        fluctuation = np.asarray(component, dtype=np.float64)
        fluctuation = fluctuation - np.mean(fluctuation)
        transform = np.fft.rfftn(fluctuation) / norm
        energy += 0.5 * np.abs(transform) ** 2 * weight[None, None, :]
    spectral_energy = np.bincount(shell.ravel(), weights=energy.ravel())
    modes = np.arange(spectral_energy.size, dtype=np.float64)
    display_floor = np.max(spectral_energy) * 1.0e-12
    valid = (modes > 0.0) & (spectral_energy > display_floor)
    return modes[valid], spectral_energy[valid]


def _snapshot(path: Path, read_binary, retain_slice: bool) -> dict[str, Any]:
    data = read_binary(str(path))
    names = ("dens", "velx", "vely", "velz", "bcc1", "bcc2", "bcc3")
    fields = {name: np.asarray(data[name], dtype=np.float64) for name in names}
    if not all(np.all(np.isfinite(value)) for value in fields.values()):
        raise RuntimeError(f"{path}: non-finite MHD snapshot")
    density = fields["dens"]
    velocity = [fields[name] for name in ("velx", "vely", "velz")]
    magnetic = [fields[name] for name in ("bcc1", "bcc2", "bcc3")]
    velocity_mean = np.asarray([np.mean(value) for value in velocity])
    velocity_rms = math.sqrt(max(
        float(sum(np.mean(value**2) for value in velocity) -
              np.sum(velocity_mean**2)), 0.0))
    magnetic_rms = math.sqrt(float(sum(np.mean(value**2) for value in magnetic)))
    mode, energy = _spectrum(magnetic)
    result: dict[str, Any] = {
        "path": str(path.resolve()),
        "time": float(data["Time"]),
        "density_mean": float(np.mean(density)),
        "density_fractional_rms": float(np.std(density) / np.mean(density)),
        "velocity_rms": velocity_rms,
        "magnetic_rms": magnetic_rms,
        "spectrum_mode": mode.tolist(),
        "spectrum_magnetic_energy": energy.tolist(),
    }
    if retain_slice:
        mid = density.shape[0] // 2
        bmag = np.sqrt(sum(value**2 for value in magnetic))
        result["density_slice"] = density[mid]
        result["magnetic_slice"] = bmag[mid]
        result["x1"] = np.asarray(data["x1v"], dtype=np.float64)
        result["x2"] = np.asarray(data["x2v"], dtype=np.float64)
    return result


def _history_metrics(
    mhd: dict[str, np.ndarray], user: dict[str, np.ndarray],
    volume: float, driving_rate: float, threshold: float,
    tail_fraction: float, hall_mode: str,
) -> dict[str, Any]:
    required_mhd = ("time", "tot-E", "1-KE", "2-KE", "3-KE",
                    "1-ME", "2-ME", "3-ME")
    if hall_mode == "full":
        required_mhd += ("hall_Rmax", "hall_Lmax")
    required_user = ("time", "vx_vol", "vy_vol", "vz_vol", "v2_vol",
                     "divB_max", "cr_Ekin", "cr_Px", "cr_Py", "cr_Pz")
    for name in required_mhd:
        if name not in mhd:
            raise RuntimeError(f"MHD history is missing {name}")
    for name in required_user:
        if name not in user:
            raise RuntimeError(f"user history is missing {name}")
    time = mhd["time"]
    hall_lambda = mhd["hall_Lmax"] if hall_mode == "full" else None
    hall_r = mhd["hall_Rmax"] if hall_mode == "full" else None
    if (hall_mode == "full" and
            (np.any(hall_lambda < 0.0) or np.any(hall_r < 0.0))):
        raise RuntimeError("Hall maxima must be non-negative")
    tail_start = time[0] + (1.0 - tail_fraction) * (time[-1] - time[0])
    tail = time >= tail_start
    if hall_mode == "full" and np.count_nonzero(tail) < 5:
        raise RuntimeError("Hall recommendation tail has fewer than five samples")
    sustained_fraction = (
        float(np.mean(hall_lambda[tail] >= threshold))
        if hall_mode == "full" else 0.0
    )
    recommend = (hall_mode == "full") and sustained_fraction >= 0.5

    kinetic = mhd["1-KE"] + mhd["2-KE"] + mhd["3-KE"]
    magnetic = mhd["1-ME"] + mhd["2-ME"] + mhd["3-ME"]
    cr_energy = np.interp(time, user["time"], user["cr_Ekin"])
    total = mhd["tot-E"] + cr_energy
    expected = total[0] + driving_rate * (time - time[0])
    energy_residual = (total - expected) / np.maximum(np.abs(expected), 1.0e-30)

    user_time = user["time"]
    velocity_mean2 = (user["vx_vol"]**2 + user["vy_vol"]**2 +
                      user["vz_vol"]**2) / volume**2
    velocity_rms = np.sqrt(np.maximum(user["v2_vol"] / volume -
                                      velocity_mean2, 0.0))
    user_magnetic = np.interp(user_time, time, magnetic)
    magnetic_rms = np.sqrt(np.maximum(2.0 * user_magnetic / volume, 0.0))
    gas_momentum = np.vstack([
        np.interp(user_time, time, mhd[name])
        for name in ("1-mom", "2-mom", "3-mom")
    ])
    cr_momentum = np.vstack([user[name] for name in ("cr_Px", "cr_Py", "cr_Pz")])
    total_momentum = np.sqrt(np.sum((gas_momentum + cr_momentum) ** 2, axis=0))
    return {
        "time": time.tolist(),
        "kinetic_energy": kinetic.tolist(),
        "magnetic_energy": magnetic.tolist(),
        "cr_energy": cr_energy.tolist(),
        "energy_residual_fraction": energy_residual.tolist(),
        "user_time": user_time.tolist(),
        "velocity_rms": velocity_rms.tolist(),
        "magnetic_rms": magnetic_rms.tolist(),
        "total_momentum_norm": total_momentum.tolist(),
        "divb_max": user["divB_max"].tolist(),
        "hall_rmax": hall_r.tolist() if hall_mode == "full" else None,
        "hall_lambda_max": hall_lambda.tolist() if hall_mode == "full" else None,
        "tail_start": float(tail_start),
        "tail_lambda_median": (
            float(np.median(hall_lambda[tail])) if hall_mode == "full" else None
        ),
        "tail_lambda_p90": (
            float(np.quantile(hall_lambda[tail], 0.9))
            if hall_mode == "full" else None
        ),
        "tail_lambda_peak": (
            float(np.max(hall_lambda[tail])) if hall_mode == "full" else None
        ),
        "tail_rmax_median": (
            float(np.median(hall_r[tail])) if hall_mode == "full" else None
        ),
        "tail_fraction_at_or_above_threshold": (
            sustained_fraction if hall_mode == "full" else None
        ),
        "recommend_full_size_hall_off": recommend,
    }


def _style() -> None:
    plt.rcParams.update({
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "axes.linewidth": 0.9,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.frameon": False,
        "lines.linewidth": 1.8,
    })


def _plot_history(
    metrics: dict[str, Any], threshold: float, hall_mode: str, path: Path
) -> None:
    time = np.asarray(metrics["time"])
    user_time = np.asarray(metrics["user_time"])
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2), constrained_layout=True)
    axes[0, 0].semilogy(time, metrics["kinetic_energy"], label=r"$E_{\rm kin}$")
    axes[0, 0].semilogy(time, metrics["magnetic_energy"], label=r"$E_B$")
    axes[0, 0].semilogy(time, metrics["cr_energy"], label=r"$E_{\rm CR}$")
    axes[0, 0].set(xlabel="time", ylabel="volume-integrated energy")
    axes[0, 0].legend()

    axes[0, 1].plot(user_time, metrics["velocity_rms"], label=r"$v_{\rm rms}$")
    axes[0, 1].plot(user_time, metrics["magnetic_rms"], label=r"$B_{\rm rms}$")
    axes[0, 1].set(xlabel="time", ylabel="RMS amplitude")
    axes[0, 1].legend()

    axes[1, 0].plot(time, metrics["energy_residual_fraction"],
                    color="#6f4e7c", label="driven energy residual")
    axes[1, 0].plot(user_time, metrics["total_momentum_norm"],
                    color="#2a9d8f", label="gas+CR momentum norm")
    axes[1, 0].set(xlabel="time", ylabel="conservation diagnostic")
    axes[1, 0].legend()

    if hall_mode == "full":
        axes[1, 1].semilogy(time, np.maximum(metrics["hall_rmax"], 1.0e-16),
                           label=r"$\max |R|$")
        axes[1, 1].semilogy(time, np.maximum(metrics["hall_lambda_max"], 1.0e-16),
                           label=r"$\max \Lambda$")
        axes[1, 1].axhline(threshold, color="0.25", linestyle="--",
                           label=rf"comparison threshold $\Lambda={threshold:g}$")
        axes[1, 1].axvspan(metrics["tail_start"], time[-1],
                           color="0.7", alpha=0.18)
        axes[1, 1].set(xlabel="time", ylabel="instantaneous domain maximum")
        axes[1, 1].legend()
    else:
        axes[1, 1].axis("off")
        axes[1, 1].text(
            0.5, 0.55, "CR-Hall closure disabled\nmatched control",
            transform=axes[1, 1].transAxes, ha="center", va="center", fontsize=12,
        )

    for axis in axes.flat:
        axis.grid(alpha=0.18)
    if hall_mode == "full":
        decision = "RUN Hall-off control" if metrics["recommend_full_size_hall_off"] \
            else "skip full-size Hall-off control"
        title = f"128$^3$ full-CR-Hall turbulent dynamo: {decision}"
    else:
        title = "128$^3$ Hall-off turbulent-dynamo control"
    fig.suptitle(title)
    fig.text(0.995, 0.995, EVIDENCE_NOTE, ha="right", va="top",
             fontsize=7, color="0.35")
    fig.savefig(path, dpi=300)
    plt.close(fig)


def _plot_slices(
    snapshots: list[dict[str, Any]], hall_mode: str, path: Path
) -> None:
    density = [item["density_slice"] / item["density_mean"] - 1.0
               for item in snapshots]
    log_b = [np.log10(np.maximum(item["magnetic_slice"] / item["magnetic_rms"],
                                 np.finfo(float).tiny)) for item in snapshots]
    density_limit = max(float(np.quantile(np.abs(item), 0.995)) for item in density)
    b_limits = (min(float(np.quantile(item, 0.01)) for item in log_b),
                max(float(np.quantile(item, 0.99)) for item in log_b))
    fig, axes = plt.subplots(2, len(snapshots), figsize=(12.0, 6.8),
                             constrained_layout=True)
    density_image = None
    magnetic_image = None
    for column, item in enumerate(snapshots):
        extent = [item["x1"][0], item["x1"][-1],
                  item["x2"][0], item["x2"][-1]]
        density_image = axes[0, column].imshow(
            density[column], origin="lower", extent=extent, cmap="RdBu_r",
            vmin=-density_limit, vmax=density_limit, interpolation="nearest")
        magnetic_image = axes[1, column].imshow(
            log_b[column], origin="lower", extent=extent, cmap="magma",
            vmin=b_limits[0], vmax=b_limits[1], interpolation="nearest")
        axes[0, column].set_title(f"$t={item['time']:.2f}$")
        axes[1, column].set_xlabel("$x$")
        if column == 0:
            axes[0, column].set_ylabel("$y$")
            axes[1, column].set_ylabel("$y$")
        for axis in axes[:, column]:
            axis.set_aspect("equal")
    fig.colorbar(density_image, ax=axes[0, :], shrink=0.82,
                 label=r"$\rho/\langle\rho\rangle-1$")
    fig.colorbar(magnetic_image, ax=axes[1, :], shrink=0.82,
                 label=r"$\log_{10}(|B|/B_{\rm rms})$")
    model = "full-CR-Hall" if hall_mode == "full" else "Hall-off"
    fig.suptitle(f"Midplane structure in the {model} turbulent box")
    fig.text(0.995, 0.995, EVIDENCE_NOTE, ha="right", va="top",
             fontsize=7, color="0.35")
    fig.savefig(path, dpi=300)
    plt.close(fig)


def _plot_spectra(
    snapshots: list[dict[str, Any]], hall_mode: str, path: Path
) -> None:
    fig, axis = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    for color, item in zip(COLORS, snapshots):
        mode = np.asarray(item["spectrum_mode"])
        energy = np.asarray(item["spectrum_magnetic_energy"])
        axis.loglog(mode, energy, color=color, label=f"$t={item['time']:.2f}$")
    final_mode = np.asarray(snapshots[-1]["spectrum_mode"])
    final_energy = np.asarray(snapshots[-1]["spectrum_magnetic_energy"])
    anchor_index = int(np.argmin(np.abs(final_mode - 8.0)))
    guide = final_energy[anchor_index] * (final_mode / final_mode[anchor_index]) ** (-5/3)
    axis.loglog(final_mode, guide, color="0.35", linestyle="--", linewidth=1.2,
                label=r"$k^{-5/3}$ guide")
    axis.set(xlabel=r"isotropic mode $kL/(2\pi)$",
             ylabel=r"shell magnetic energy $E_B(k)$")
    axis.grid(alpha=0.18, which="both")
    axis.legend()
    model = "full-CR-Hall" if hall_mode == "full" else "Hall-off"
    axis.set_title(f"{model} magnetic-energy spectrum")
    fig.text(0.995, 0.005, EVIDENCE_NOTE, ha="right", va="bottom",
             fontsize=7, color="0.35")
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--expected-hall-mode", choices=("full", "off"), default="full"
    )
    parser.add_argument("--lambda-threshold", type=float, default=0.1)
    parser.add_argument("--tail-fraction", type=float, default=0.3)
    args = parser.parse_args()
    if args.lambda_threshold <= 0.0:
        raise SystemExit("--lambda-threshold must be positive")
    if not 0.0 < args.tail_fraction <= 1.0:
        raise SystemExit("--tail-fraction must lie in (0, 1]")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    paths = sorted((args.run / "bin").glob("*.mhd_w_bcc.*.bin"))
    if len(paths) < 3:
        paths = sorted(args.run.glob("*.mhd_w_bcc.*.bin"))
    if len(paths) < 3:
        raise RuntimeError(f"{args.run}: fewer than three MHD snapshots")
    headers = [(_binary_header(path), path) for path in paths]
    headers.sort(key=lambda item: item[0]["time"])
    parameters = _parameters(headers[0][0]["input"])
    model = _validate_model(parameters, args.expected_hall_mode)
    if headers[-1][0]["time"] < 0.99 * model["tlim"]:
        raise RuntimeError("MHD snapshots do not reach the requested terminal time")
    selected_indices = sorted({0, len(headers) // 2, len(headers) - 1})
    selected_paths = [headers[index][1] for index in selected_indices]

    read_binary = _binary_reader()
    snapshots = [_snapshot(path, read_binary, True) for path in selected_paths]
    lengths = [
        snapshots[0]["x1"][-1] - snapshots[0]["x1"][0] +
        np.mean(np.diff(snapshots[0]["x1"])),
        snapshots[0]["x2"][-1] - snapshots[0]["x2"][0] +
        np.mean(np.diff(snapshots[0]["x2"])),
    ]
    x3min = _real(parameters, "mesh", "x3min")
    x3max = _real(parameters, "mesh", "x3max")
    volume = float(lengths[0] * lengths[1] * (x3max - x3min))

    mhd = _history(_one_file(args.run, "mhd"))
    user = _history(_one_file(args.run, "user"))
    metrics = _history_metrics(
        mhd, user, volume, model["driving_energy_rate"],
        args.lambda_threshold, args.tail_fraction, args.expected_hall_mode,
    )
    energy_residual = np.asarray(metrics["energy_residual_fraction"])
    total_momentum = np.asarray(metrics["total_momentum_norm"])
    divb = np.asarray(metrics["divb_max"])
    measured = {
        "mhd_terminal_time": float(mhd["time"][-1]),
        "user_terminal_time": float(user["time"][-1]),
        "max_abs_energy_residual": float(np.max(np.abs(energy_residual))),
        "max_total_momentum_norm": float(np.max(total_momentum)),
        "max_abs_divb": float(np.max(np.abs(divb))),
    }
    acceptance = {
        "limits": {
            "max_abs_energy_residual": ENERGY_RESIDUAL_LIMIT,
            "max_total_momentum_norm": TOTAL_MOMENTUM_LIMIT,
            "max_abs_divb": DIVB_LIMIT,
        },
        "measured": measured,
        "checks": {
            "mhd_history_terminal": (
                measured["mhd_terminal_time"] >= 0.99 * model["tlim"]
            ),
            "user_history_terminal": (
                measured["user_terminal_time"] >= 0.99 * model["tlim"]
            ),
            "driven_energy_residual": (
                measured["max_abs_energy_residual"] <= ENERGY_RESIDUAL_LIMIT
            ),
            "total_momentum": (
                measured["max_total_momentum_norm"] <= TOTAL_MOMENTUM_LIMIT
            ),
            "divb": measured["max_abs_divb"] <= DIVB_LIMIT,
        },
    }
    acceptance["passed"] = all(acceptance["checks"].values())
    recommendation = None
    if args.expected_hall_mode == "full":
        recommendation = {
            "run_full_size_hall_off": metrics["recommend_full_size_hall_off"],
            "criterion": (
                "Recommend the full-size Hall-off control when instantaneous max Lambda "
                f"is at least {args.lambda_threshold:g} in at least half of the final "
                f"{100.0 * args.tail_fraction:g}% of history samples."
            ),
            "physical_basis": (
                "Lambda=|v_H|/v_A measures the Hall drift relative to the Alfv\u00e9n "
                "speed; a sustained value at the selected threshold makes an "
                "induction-scale correction of that order plausible somewhere in "
                "the domain."
            ),
            "tail_fraction_at_or_above_threshold":
                metrics["tail_fraction_at_or_above_threshold"],
            "tail_lambda_median": metrics["tail_lambda_median"],
            "tail_lambda_p90": metrics["tail_lambda_p90"],
            "tail_lambda_peak": metrics["tail_lambda_peak"],
            "limitation": (
                "The history stores an instantaneous domain maximum, not a volume "
                "filling fraction; isolated spikes alone do not trigger the comparison."
            ),
        }
    clean_snapshots = []
    for item in snapshots:
        clean_snapshots.append({
            key: value for key, value in item.items()
            if key not in {"density_slice", "magnetic_slice", "x1", "x2"}
        })
    report = {
        "schema_version": 1,
        "evidence_class": "engineering_proxy",
        "not_sun_bai_reproduction": True,
        "run_directory": str(args.run.resolve()),
        "model": model,
        "history": metrics,
        "snapshots": clean_snapshots,
        "acceptance": acceptance,
        "hall_off_recommendation": recommendation,
    }

    _style()
    mode_slug = "full_hall" if args.expected_hall_mode == "full" else "hall_off"
    _plot_history(
        metrics, args.lambda_threshold, args.expected_hall_mode,
        args.output_dir / f"turbulent_dynamo_{mode_slug}_history.png",
    )
    _plot_slices(
        snapshots, args.expected_hall_mode,
        args.output_dir / f"turbulent_dynamo_{mode_slug}_slices.png",
    )
    _plot_spectra(
        snapshots, args.expected_hall_mode,
        args.output_dir / f"turbulent_dynamo_{mode_slug}_spectra.png",
    )
    (args.output_dir / f"turbulent_dynamo_{mode_slug}_analysis.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(
        {"acceptance": acceptance, "hall_off_recommendation": recommendation},
        indent=2, sort_keys=True
    ))
    return 0 if acceptance["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
