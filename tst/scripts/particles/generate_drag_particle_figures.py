#!/usr/bin/env python3
"""Generate documentation figures for drag-coupled particle validation."""

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
EXPECTED_GROWTH = 0.1424527496587504


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--athena", default="build/src/athena")
    parser.add_argument("--mpi-athena", default=None)
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--np", type=int, default=2)
    parser.add_argument(
        "--output-dir",
        default="docs/source/_static/particles",
        help="Directory for generated PNG and JSON documentation assets",
    )
    parser.add_argument(
        "--work-dir",
        default=None,
        help="Directory for AthenaK run products; temporary by default",
    )
    parser.add_argument(
        "--skip-runs",
        action="store_true",
        help="Reuse histories and VTK files already present in --work-dir",
    )
    parser.add_argument(
        "--keep-work",
        action="store_true",
        help="Do not remove the temporary run directory after generating figures",
    )
    return parser.parse_args()


def resolve_path(path):
    path = Path(path)
    return path if path.is_absolute() else REPO_ROOT / path


def read_hst(path):
    sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
    import athena_read  # noqa: WPS433

    return athena_read.hst(str(path))


def clean_case_dir(work_dir, name):
    case_dir = work_dir / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    return case_dir


def run_athena(binary, input_file, case_dir, overrides=None, mpi_np=1, mpirun="mpirun"):
    overrides = overrides or []
    binary = resolve_path(binary)
    input_file = resolve_path(input_file)
    command = [
        str(binary),
        "-i",
        str(input_file),
        "-d",
        str(case_dir),
    ] + overrides
    if mpi_np > 1:
        command = [mpirun, "-np", str(mpi_np)] + command

    process = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if process.returncode != 0:
        raise RuntimeError(
            "AthenaK run failed:\n"
            f"{' '.join(command)}\n\n"
            f"stdout:\n{process.stdout}\n\nstderr:\n{process.stderr}"
        )


def run_cases(args, work_dir):
    mpi_athena = args.mpi_athena or args.athena
    cases = {
        "relaxation": "inputs/tests/particle_drag_relaxation.athinput",
        "energy": "inputs/tests/particle_drag_relaxation.athinput",
        "uniform": "inputs/particles/streaming_instability.athinput",
        "amr_mpi": "inputs/particles/streaming_instability_amr.athinput",
        "dynamic_amr": "inputs/particles/streaming_instability_amr_dynamic.athinput",
        "visualization": "inputs/particles/streaming_instability_visualization.athinput",
    }
    case_dirs = {name: work_dir / name for name in cases}
    if args.skip_runs:
        return case_dirs

    for name in cases:
        case_dirs[name] = clean_case_dir(work_dir, name)

    run_athena(args.athena, cases["relaxation"], case_dirs["relaxation"])
    run_athena(
        args.athena,
        cases["energy"],
        case_dirs["energy"],
        overrides=[
            "job/basename=particle_drag_energy_docs",
            "time/tlim=0.1",
            "hydro/eos=ideal",
            "hydro/gamma=1.6666666666666667",
            "drag_particles/include_energy=true",
        ],
    )
    run_athena(args.athena, cases["uniform"], case_dirs["uniform"])
    run_athena(
        mpi_athena,
        cases["amr_mpi"],
        case_dirs["amr_mpi"],
        overrides=["time/tlim=4.0", "output1/dt=0.05"],
        mpi_np=args.np,
        mpirun=args.mpirun,
    )
    run_athena(
        mpi_athena,
        cases["dynamic_amr"],
        case_dirs["dynamic_amr"],
        mpi_np=args.np,
        mpirun=args.mpirun,
    )
    run_athena(args.athena, cases["visualization"], case_dirs["visualization"])
    return case_dirs


def mode_amplitude(data):
    gas_mode = np.hypot(data["gmodec"], data["gmodes"])
    particle_mode = np.hypot(data["pmodec"], data["pmodes"])
    return gas_mode + particle_mode


def fit_growth(data, t_min, t_max):
    amp = mode_amplitude(data)
    mask = (data["time"] >= t_min) & (data["time"] <= t_max) & (amp > 0.0)
    return float(np.polyfit(data["time"][mask], np.log(amp[mask]), 1)[0])


def relative_drift(values):
    scale = max(abs(float(values[0])), 1.0)
    return (values - values[0]) / scale


def floor_for_log(values, floor=1.0e-16):
    return np.maximum(np.abs(values), floor)


def latest_file(directory, pattern):
    files = sorted(directory.glob(pattern))
    if not files:
        raise FileNotFoundError(f"no files matching {pattern} in {directory}")
    return files[-1]


def read_vtk_scalar(path):
    dims = origin = spacing = None
    ncell = None
    scalar_name = None
    with open(path, "rb") as stream:
        while True:
            raw_line = stream.readline()
            if raw_line == b"":
                break
            line = raw_line.decode("ascii", errors="ignore").strip()
            if line.startswith("DIMENSIONS"):
                dims = tuple(int(x) for x in line.split()[1:4])
            elif line.startswith("ORIGIN"):
                origin = tuple(float(x) for x in line.split()[1:4])
            elif line.startswith("SPACING"):
                spacing = tuple(float(x) for x in line.split()[1:4])
            elif line.startswith("CELL_DATA"):
                ncell = int(line.split()[1])
            elif line.startswith("SCALARS"):
                scalar_name = line.split()[1]
            elif line.startswith("LOOKUP_TABLE"):
                raw = stream.read(4 * ncell)
                data = np.frombuffer(raw, dtype=">f4", count=ncell).astype(float)
                break
        else:
            raise ValueError(f"could not read scalar data from {path}")

    if dims is None or origin is None or spacing is None or ncell is None:
        raise ValueError(f"incomplete VTK header in {path}")
    nx, ny, nz = tuple(dim - 1 if dim > 1 else 1 for dim in dims)
    if nx * ny * nz != ncell:
        raise ValueError(f"VTK cell count mismatch in {path}")
    data = data.reshape((nz, ny, nx))
    centers = tuple(
        origin[axis] + (np.arange(n) + 0.5) * spacing[axis]
        for axis, n in enumerate((nx, ny, nz))
    )
    return {
        "name": scalar_name,
        "data": data,
        "x": centers[0],
        "y": centers[1],
        "z": centers[2],
        "spacing": spacing,
    }


def read_particle_points(path):
    with open(path, "rb") as stream:
        while True:
            raw_line = stream.readline()
            if raw_line == b"":
                break
            line = raw_line.decode("ascii", errors="ignore").strip()
            if line.startswith("POINTS"):
                npoint = int(line.split()[1])
                raw = stream.read(4 * 3 * npoint)
                return np.frombuffer(raw, dtype=">f4").reshape((npoint, 3)).astype(float)
    raise ValueError(f"could not read particle points from {path}")


def configure_matplotlib():
    plt.rcParams.update({
        "figure.dpi": 140,
        "savefig.dpi": 180,
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "legend.fontsize": 8,
        "axes.grid": True,
        "grid.alpha": 0.25,
    })


def plot_drag_relaxation(relaxation, energy, output_dir):
    gas_v = relaxation["gmom1"] / relaxation["gmass"]
    part_v = relaxation["pmom1"] / relaxation["pmass"]
    rel_v = np.abs(part_v - gas_v)
    total_mom = relaxation["gmom1"] + relaxation["pmom1"]
    total_energy = energy["gtotE"] + energy["pkinE"]

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.5))
    axes[0].semilogy(relaxation["time"], rel_v / rel_v[0], label="AthenaK")
    axes[0].semilogy(
        relaxation["time"],
        np.exp(-2.0 * relaxation["time"] / 0.1),
        "--",
        label=r"$\exp[-(1+\epsilon)t/t_s]$",
    )
    axes[0].set_xlabel("time")
    axes[0].set_ylabel(r"$|v_p-u_g|/|v_p-u_g|_0$")
    axes[0].set_title("Stopping-time drag relaxation")
    axes[0].legend()

    axes[1].plot(
        relaxation["time"],
        np.abs(relative_drift(total_mom)),
        label="total momentum",
    )
    axes[1].plot(
        energy["time"],
        np.abs(relative_drift(total_energy)),
        label="total energy, ideal gas",
    )
    axes[1].set_yscale("log")
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("relative drift")
    axes[1].set_title("Conserved exchange diagnostics")
    axes[1].legend()
    fig.tight_layout()
    path = output_dir / "drag_relaxation_conservation.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_streaming_growth(uniform, amr, output_dir):
    uniform_gamma = fit_growth(uniform, 0.0, 8.0)
    amr_gamma = fit_growth(amr, 0.0, 4.0)
    uniform_amp = mode_amplitude(uniform)
    amr_amp = mode_amplitude(amr)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.5))
    axes[0].semilogy(
        uniform["time"],
        uniform_amp / uniform_amp[0],
        color="tab:blue",
        label="uniform history",
    )
    axes[0].semilogy(
        amr["time"],
        amr_amp / amr_amp[0],
        color="tab:orange",
        label="AMR/MPI history",
    )
    axes[0].semilogy(
        uniform["time"],
        np.exp(uniform_gamma * (uniform["time"] - uniform["time"][0])),
        "--",
        color="tab:blue",
        label=fr"uniform fit $\gamma={uniform_gamma:.4f}$",
    )
    axes[0].semilogy(
        amr["time"],
        np.exp(amr_gamma * (amr["time"] - amr["time"][0])),
        "--",
        color="tab:orange",
        label=fr"AMR/MPI fit $\gamma={amr_gamma:.4f}$",
    )
    expected_t = np.linspace(0.0, 8.0, 200)
    axes[0].semilogy(
        expected_t,
        np.exp(EXPECTED_GROWTH * expected_t),
        "k--",
        label=fr"linear $\gamma={EXPECTED_GROWTH:.4f}$",
    )
    axes[0].set_xlabel("time")
    axes[0].set_ylabel(r"$A(t)/A(0)$")
    axes[0].set_title("Streaming-instability growth")
    axes[0].legend()

    axes[1].plot(
        uniform["time"],
        floor_for_log(relative_drift(uniform["pmass"])),
        label="uniform",
    )
    axes[1].plot(
        amr["time"],
        floor_for_log(relative_drift(amr["pmass"])),
        label="AMR/MPI",
    )
    axes[1].axhline(1.0e-12, color="k", ls=":", lw=1.0, label="test limit")
    axes[1].set_yscale("log")
    axes[1].set_ylim(1.0e-16, 1.0e-10)
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("absolute particle-mass drift")
    axes[1].set_title("Particle-mass conservation")
    axes[1].legend()
    fig.tight_layout()
    path = output_dir / "streaming_growth_uniform_amr_mpi.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path, uniform_gamma, amr_gamma


def plot_amr_remap(dynamic, output_dir):
    amp = mode_amplitude(dynamic)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.5))
    axes[0].plot(dynamic["time"], amp / amp[0], marker="o", ms=3)
    axes[0].set_xlabel("time")
    axes[0].set_ylabel(r"$A(t)/A(0)$")
    axes[0].set_title("Dynamic AMR remap history")

    mass_drift = np.abs(relative_drift(dynamic["pmass"]))
    axes[1].plot(dynamic["time"], floor_for_log(mass_drift), marker="o", ms=3)
    axes[1].axhline(1.0e-12, color="k", ls=":", lw=1.0, label="test limit")
    axes[1].set_yscale("log")
    axes[1].set_ylim(1.0e-16, 1.0e-10)
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("relative particle-mass drift")
    axes[1].set_title("AMR/MPI particle conservation")
    axes[1].legend()
    if np.max(mass_drift) == 0.0:
        axes[1].text(
            0.5,
            0.5,
            "all sampled values are exactly zero",
            ha="center",
            va="center",
            transform=axes[1].transAxes,
        )
    fig.tight_layout()
    path = output_dir / "streaming_dynamic_amr_remap.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def mesh_extent(vtk_data, j_index):
    x = vtk_data["x"]
    z = vtk_data["z"]
    dx = vtk_data["spacing"][0]
    dz = vtk_data["spacing"][2]
    return [x[0] - 0.5 * dx, x[-1] + 0.5 * dx, z[0] - 0.5 * dz, z[-1] + 0.5 * dz]


def plot_visual_slice(case_dir, output_dir):
    gas = read_vtk_scalar(latest_file(case_dir / "vtk", "*gas_density*.vtk"))
    dust = read_vtk_scalar(latest_file(case_dir / "vtk", "*particle_density*.vtk"))
    points = read_particle_points(latest_file(case_dir / "pvtk", "*particles*.vtk"))
    j_index = int(np.argmin(np.abs(gas["y"])))
    y0 = gas["y"][j_index]
    dy = gas["spacing"][1]
    point_mask = np.abs(points[:, 1] - y0) <= 0.5 * dy + 1.0e-12
    slice_points = points[point_mask]

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), sharex=True, sharey=True)
    for axis, data, title, cmap in [
        (axes[0], gas, "gas density", "viridis"),
        (axes[1], dust, "particle mass density", "magma"),
    ]:
        image = axis.imshow(
            data["data"][:, j_index, :],
            origin="lower",
            extent=mesh_extent(data, j_index),
            aspect="equal",
            cmap=cmap,
        )
        axis.scatter(
            slice_points[:, 0],
            slice_points[:, 2],
            s=5,
            c="white",
            edgecolors="black",
            linewidths=0.25,
            alpha=0.85,
        )
        axis.set_title(f"{title}, y={y0:.3f}")
        axis.set_xlabel("x")
        fig.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    axes[0].set_ylabel("z")
    fig.suptitle("Streaming visualization: mesh slice with particle markers")
    fig.tight_layout()
    path = output_dir / "streaming_slice_particles.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path, int(points.shape[0]), int(slice_points.shape[0])


def write_summary(output_dir, summary):
    path = output_dir / "drag_particle_validation_summary.json"
    with open(path, "w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return path


def main():
    args = parse_args()
    configure_matplotlib()
    output_dir = resolve_path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.skip_runs and args.work_dir is None:
        raise ValueError("--skip-runs requires --work-dir")

    owned_work_dir = args.work_dir is None
    work_dir = (
        Path(tempfile.mkdtemp(prefix="drag-particle-figures-"))
        if owned_work_dir else resolve_path(args.work_dir)
    )
    work_dir.mkdir(parents=True, exist_ok=True)

    try:
        case_dirs = run_cases(args, work_dir)
        relaxation = read_hst(
            case_dirs["relaxation"] / "particle_drag_relaxation.user.hst"
        )
        energy = read_hst(case_dirs["energy"] / "particle_drag_energy_docs.user.hst")
        uniform = read_hst(case_dirs["uniform"] / "streaming_instability.user.hst")
        amr = read_hst(case_dirs["amr_mpi"] / "streaming_instability_amr.user.hst")
        dynamic = read_hst(
            case_dirs["dynamic_amr"] / "streaming_instability_amr_dynamic.user.hst"
        )

        figure_paths = []
        figure_paths.append(plot_drag_relaxation(relaxation, energy, output_dir))
        growth_path, uniform_gamma, amr_gamma = plot_streaming_growth(
            uniform, amr, output_dir
        )
        figure_paths.append(growth_path)
        figure_paths.append(plot_amr_remap(dynamic, output_dir))
        slice_path, npoint, nslice = plot_visual_slice(
            case_dirs["visualization"], output_dir
        )
        figure_paths.append(slice_path)

        gas_v = relaxation["gmom1"] / relaxation["gmass"]
        part_v = relaxation["pmom1"] / relaxation["pmass"]
        rel_v = np.abs(part_v - gas_v)
        total_mom = relaxation["gmom1"] + relaxation["pmom1"]
        total_energy = energy["gtotE"] + energy["pkinE"]
        dynamic_mass_drift = np.max(np.abs(relative_drift(dynamic["pmass"])))
        summary = {
            "amr_mpi_growth": amr_gamma,
            "dynamic_amr_max_particle_mass_drift": float(dynamic_mass_drift),
            "expected_growth": EXPECTED_GROWTH,
            "relaxation_final_relative_velocity_ratio": float(rel_v[-1] / rel_v[0]),
            "relaxation_max_total_momentum_drift": float(
                np.max(np.abs(relative_drift(total_mom)))
            ),
            "energy_max_total_energy_drift": float(
                np.max(np.abs(relative_drift(total_energy)))
            ),
            "uniform_growth": uniform_gamma,
            "visualization_particles_total": npoint,
            "visualization_particles_in_slice": nslice,
        }
        summary_path = write_summary(output_dir, summary)
        print("Generated:")
        for path in figure_paths + [summary_path]:
            print(f"  {path.relative_to(REPO_ROOT)}")
    finally:
        if owned_work_dir and not args.keep_work:
            shutil.rmtree(work_dir)
        elif args.keep_work:
            print(f"Kept work directory: {work_dir}")


if __name__ == "__main__":
    main()
