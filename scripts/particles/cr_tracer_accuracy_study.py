#!/usr/bin/env python3
"""Run CR tracer accuracy studies and generate documentation figures."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict
from typing import Iterable
from typing import List
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import cr_tracer_inspect  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]
INPUTS = ROOT / "inputs" / "particles" / "accuracy"
DEFAULT_OUTPUT = ROOT / "docs" / "source" / "modules" / "figures" / (
    "cr_tracer_accuracy")


def _artifact_path(path: Path) -> str:
    """Return a repository-relative figure path when possible."""
    path = path.resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def _default_athena() -> Path:
    for candidate in (
        ROOT / "build-cr-robust-mpi" / "src" / "athena",
        ROOT / "build-cr-perf-mpi" / "src" / "athena",
        ROOT / "tst" / "build" / "src" / "athena",
    ):
        if candidate.exists():
            return candidate
    return ROOT / "build-cr-robust-mpi" / "src" / "athena"


def _run(args: argparse.Namespace, case: str, input_name: str,
         overrides: Sequence[str], ranks: int = 1) -> Path:
    run_dir = args.run_root / case
    shutil.rmtree(run_dir, ignore_errors=True)
    run_dir.mkdir(parents=True, exist_ok=True)
    if ranks > 1:
        command = [args.mpiexec, "-np", str(ranks), str(args.athena)]
    else:
        command = [str(args.athena)]
    command += ["-i", str(INPUTS / input_name), "-d", str(run_dir)]
    command += list(overrides)
    proc = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"{case} failed\nCOMMAND: {' '.join(command)}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    (run_dir / "stdout.txt").write_text(proc.stdout)
    (run_dir / "stderr.txt").write_text(proc.stderr)
    return run_dir


def _single_particle(summary: Dict) -> Dict:
    particles = [
        particle
        for restart in summary["restart"]
        for particle in restart["particles"]
    ]
    if len(particles) != 1:
        raise ValueError(f"Expected one particle, found {len(particles)}")
    return particles[0]


def _state(summary: Dict) -> Dict:
    particle = _single_particle(summary)
    return {
        "position": tuple(particle["reals"][0:3]),
        "velocity": tuple(particle["reals"][3:6]),
        "bfield": tuple(particle["reals"][7:10]),
    }


def _norm3(values: Sequence[float]) -> float:
    return math.sqrt(sum(value*value for value in values))


def _vector_error(a: Sequence[float], b: Sequence[float]) -> float:
    return _norm3([a[i] - b[i] for i in range(3)])


def _magnetic_moment_proxy(velocity: Sequence[float],
                           bfield: Sequence[float]) -> float:
    bmag = _norm3(bfield)
    speed = _norm3(velocity)
    vpar = sum(velocity[i]*bfield[i] for i in range(3))/bmag
    return (speed*speed - vpar*vpar)/bmag


def _phase_error(vx: float, vy: float, time: float) -> float:
    angle = math.atan2(vy, vx)
    target = -time
    return abs(math.atan2(math.sin(angle - target), math.cos(angle - target)))


def _fit_log_slope(x_values: Iterable[float], y_values: Iterable[float]) -> float:
    xs = [math.log(float(x)) for x in x_values]
    ys = [math.log(float(y)) for y in y_values]
    xbar = sum(xs)/len(xs)
    ybar = sum(ys)/len(ys)
    denom = sum((x - xbar)**2 for x in xs)
    return sum((xs[i] - xbar)*(ys[i] - ybar) for i in range(len(xs)))/denom


def _exact_bfield(profile: str, x: float, y: float, z: float,
                  b0: Sequence[float], bgrad: float = 1.0,
                  bamp: float = 0.05, bwave: float = 1.0) -> List[float]:
    bx, by, bz = b0
    if profile == "linear_cross":
        return [bx + bgrad*y, by + bgrad*x, bz]
    if profile in ("sinusoidal_divb_free", "turbulent"):
        k = 2.0*math.pi*bwave
        bx += bamp*k*math.sin(k*x)*(math.cos(k*y) - math.cos(k*z))
        by += bamp*k*math.sin(k*y)*(math.cos(k*z) - math.cos(k*x))
        bz += bamp*k*math.sin(k*z)*(math.cos(k*x) - math.cos(k*y))
        if profile == "turbulent":
            k2 = 2.0*k
            amp2 = 0.35*bamp
            bx += amp2*k2*math.sin(k2*(x + 0.13))*(
                math.cos(k2*(y - 0.07)) - math.cos(k2*(z + 0.11)))
            by += amp2*k2*math.sin(k2*(y - 0.07))*(
                math.cos(k2*(z + 0.11)) - math.cos(k2*(x + 0.13)))
            bz += amp2*k2*math.sin(k2*(z + 0.11))*(
                math.cos(k2*(x + 0.13)) - math.cos(k2*(y - 0.07)))
        return [bx, by, bz]
    if profile == "uniform":
        return list(b0)
    if profile == "mirror":
        bx += -bgrad*x*z
        by += -bgrad*y*z
        bz += bgrad*z*z
        return [bx, by, bz]
    if profile == "gradb":
        return [bx, by, bz + bgrad*x]
    raise ValueError(f"Unsupported profile {profile}")


def _sample_error(summary: Dict, profile: str, b0: Sequence[float],
                  bgrad: float = 1.0, bamp: float = 0.05,
                  bwave: float = 1.0) -> float:
    errors = []
    for restart in summary["restart"]:
        for particle in restart["particles"]:
            x, y, z = particle["reals"][0:3]
            sampled = particle["reals"][7:10]
            exact = _exact_bfield(profile, x, y, z, b0, bgrad, bamp, bwave)
            errors.append(_vector_error(sampled, exact))
    return math.sqrt(sum(error*error for error in errors)/len(errors))


def _field_plane(profile: str, b0: Sequence[float], bgrad: float,
                 bamp: float, bwave: float, plane: str, fixed: float,
                 limits: Sequence[float], npoint: int = 121):
    coords = np.linspace(limits[0], limits[1], npoint)
    du = np.zeros((npoint, npoint))
    dv = np.zeros((npoint, npoint))
    dbz = np.zeros((npoint, npoint))
    for j, vcoord in enumerate(coords):
        for i, ucoord in enumerate(coords):
            if plane == "xy":
                position = (ucoord, vcoord, fixed)
                components = (0, 1)
            elif plane == "xz":
                position = (ucoord, fixed, vcoord)
                components = (0, 2)
            else:
                raise ValueError(f"Unsupported qualitative plane {plane}")
            bfield = _exact_bfield(
                profile, *position, b0, bgrad=bgrad, bamp=bamp, bwave=bwave)
            du[j, i] = bfield[components[0]]
            dv[j, i] = bfield[components[1]]
            dbz[j, i] = bfield[2] - b0[2]
    return coords, du, dv, dbz


def _draw_field_plane(ax, config: Dict, limits: Sequence[float]):
    coords, du, dv, dbz = _field_plane(
        config["profile"], config["b0"], config["bgrad"], config["bamp"],
        config["bwave"], config["plane"], config["fixed"], limits)
    limit = max(float(np.max(np.abs(dbz))), 1.0e-12)
    image = ax.pcolormesh(
        coords, coords, dbz, cmap="RdBu_r", shading="auto",
        vmin=-limit, vmax=limit)
    transverse = np.hypot(du, dv)
    unit_u = np.divide(du, transverse, out=np.zeros_like(du),
                       where=transverse > 0.0)
    unit_v = np.divide(dv, transverse, out=np.zeros_like(dv),
                       where=transverse > 0.0)
    if float(np.max(transverse)) > 0.0:
        ax.streamplot(
            coords, coords, unit_u, unit_v, density=0.9,
            color="#262626", linewidth=0.55, arrowsize=0.65)
    ax.set_title(config["title"])
    ax.set_xlabel(config["plane"][0])
    ax.set_ylabel(config["plane"][1])
    ax.set_aspect("equal")
    ax.set_xlim(*limits)
    ax.set_ylim(*limits)
    return image


def _particle_tracks(run_dir: Path, per_species: int) -> Dict:
    snapshots = [
        cr_tracer_inspect.read_prst_file(path)
        for path in (run_dir / "prst").rglob("*.prst")
    ]
    snapshots.sort(key=lambda snapshot: snapshot["time"])
    initial_time = snapshots[0]["time"]
    first = [
        particle
        for snapshot in snapshots
        if math.isclose(snapshot["time"], initial_time)
        for particle in snapshot["particles"]
    ]
    selected = []
    for species in sorted({particle["species"] for particle in first}):
        species_tags = sorted(
            particle["tag"] for particle in first
            if particle["species"] == species)
        selected.extend((species, tag) for tag in species_tags[:per_species])

    tracks = {key: [] for key in selected}
    for snapshot in snapshots:
        for particle in snapshot["particles"]:
            key = (particle["species"], particle["tag"])
            if key in tracks:
                tracks[key].append(
                    (snapshot["time"], *particle["reals"][0:3]))
    return tracks


def _project_track(points: Sequence, plane: str, unwrap: bool) -> np.ndarray:
    values = np.asarray(points)
    indices = (1, 2) if plane == "xy" else (1, 3)
    projected = values[:, indices].copy()
    if unwrap:
        for axis in range(2):
            jumps = np.diff(projected[:, axis])
            offsets = np.where(jumps > 0.5, -1.0, 0.0)
            offsets += np.where(jumps < -0.5, 1.0, 0.0)
            projected[1:, axis] += np.cumsum(offsets)
    else:
        jumps = np.max(np.abs(np.diff(projected, axis=0)), axis=1)
        projected[1:][jumps > 0.5] = np.nan
    return projected


QUALITATIVE_CASES = [
    {
        "case": "uniform_gyro",
        "title": "1. Uniform-B gyro orbit",
        "input": "cr_uniform_gyro.athinput",
        "profile": "uniform",
        "b0": (0.0, 0.0, 1.0),
        "bgrad": 1.0,
        "bamp": 0.05,
        "bwave": 1.0,
        "plane": "xy",
        "fixed": 0.0,
        "unwrap": True,
        "limits": (-2.1, 1.1),
        "per_species": 1,
        "ranks": 1,
        "overrides": [
            "job/basename=cr_acc_qual_uniform_gyro",
            "particles/cfl_part=0.025",
            "time/nlim=512",
            "time/tlim=6.4",
            "output1/dt=0.04",
            "output2/dt=6.4",
        ],
    },
    {
        "case": "uniform_amr_mpi",
        "title": "2. Uniform-B AMR/MPI crossing",
        "input": "cr_uniform_gyro_amr.athinput",
        "profile": "uniform",
        "b0": (0.0, 0.0, 1.0),
        "bgrad": 1.0,
        "bamp": 0.05,
        "bwave": 1.0,
        "plane": "xy",
        "fixed": 0.0,
        "unwrap": False,
        "limits": (-0.5, 0.5),
        "per_species": 2,
        "ranks": 2,
        "overrides": [
            "job/basename=cr_acc_qual_uniform_amr_mpi",
            "time/nlim=512",
            "time/tlim=0.64",
            "output1/dt=0.02",
            "output2/dt=0.64",
            "output3/dt=0.64",
        ],
    },
    {
        "case": "linear_gather",
        "title": "3. Linear gather sample point",
        "input": "cr_linear_gather.athinput",
        "profile": "linear_cross",
        "b0": (0.3, -0.2, 1.0),
        "bgrad": 0.7,
        "bamp": 0.05,
        "bwave": 1.0,
        "plane": "xy",
        "fixed": 0.047,
        "unwrap": False,
        "limits": (-0.5, 0.5),
        "per_species": 1,
        "ranks": 1,
        "overrides": [
            "job/basename=cr_acc_qual_linear_gather",
            "particles/interpolation=trilinear",
            "time/nlim=1",
            "time/tlim=0.01",
            "output1/dt=0.01",
        ],
    },
    {
        "case": "manufactured_gather",
        "title": "4. Smooth divergence-free gather",
        "input": "cr_manufactured_divb_free.athinput",
        "profile": "sinusoidal_divb_free",
        "b0": (0.1, -0.2, 1.0),
        "bgrad": 1.0,
        "bamp": 0.025,
        "bwave": 1.0,
        "plane": "xy",
        "fixed": -0.11,
        "unwrap": False,
        "limits": (-0.5, 0.5),
        "per_species": 1,
        "ranks": 1,
        "overrides": [
            "job/basename=cr_acc_qual_manufactured_gather",
            "particles/interpolation=trilinear",
            "time/nlim=1",
            "time/tlim=0.01",
            "output1/dt=0.01",
        ],
    },
    {
        "case": "smooth_orbit_reference",
        "title": "5. Smooth-field orbit",
        "input": "cr_smooth_orbit_reference.athinput",
        "profile": "sinusoidal_divb_free",
        "b0": (0.1, -0.2, 1.0),
        "bgrad": 1.0,
        "bamp": 0.015,
        "bwave": 1.0,
        "plane": "xy",
        "fixed": 0.07,
        "unwrap": False,
        "limits": (-0.5, 0.5),
        "per_species": 1,
        "ranks": 1,
        "overrides": [
            "job/basename=cr_acc_qual_smooth_orbit",
            "time/nlim=256",
            "time/tlim=1.6",
            "output1/dt=0.02",
            "output2/dt=1.6",
        ],
    },
    {
        "case": "magnetic_mirror",
        "title": "6. Magnetic mirror",
        "input": "cr_magnetic_mirror.athinput",
        "profile": "mirror",
        "b0": (0.0, 0.0, 1.0),
        "bgrad": 0.8,
        "bamp": 0.05,
        "bwave": 1.0,
        "plane": "xz",
        "fixed": -0.03,
        "unwrap": False,
        "limits": (-0.5, 0.5),
        "per_species": 1,
        "ranks": 1,
        "overrides": [
            "job/basename=cr_acc_qual_mirror",
            "time/nlim=256",
            "time/tlim=1.6",
            "output1/dt=0.02",
            "output2/dt=1.6",
        ],
    },
    {
        "case": "gradb_drift",
        "title": "7. Grad-B drift",
        "input": "cr_gradb_curvature_drift.athinput",
        "profile": "gradb",
        "b0": (0.0, 0.0, 1.0),
        "bgrad": 1.0,
        "bamp": 0.05,
        "bwave": 1.0,
        "plane": "xy",
        "fixed": 0.0,
        "unwrap": False,
        "limits": (-0.5, 0.5),
        "per_species": 1,
        "ranks": 1,
        "overrides": [
            "job/basename=cr_acc_qual_gradb",
            "problem/Bgrad=1.0",
            "particles/cfl_part=0.01",
            "particles/subcycle_gyro_fraction=0.01",
            "time/nlim=1024",
            "time/tlim=1.6",
            "output1/dt=0.02",
            "output2/dt=1.6",
        ],
    },
    {
        "case": "amr_boundary",
        "title": "8. Smooth-field AMR boundary",
        "input": "cr_amr_boundary_convergence.athinput",
        "profile": "sinusoidal_divb_free",
        "b0": (0.1, -0.2, 1.0),
        "bgrad": 1.0,
        "bamp": 0.015,
        "bwave": 1.0,
        "plane": "xy",
        "fixed": 0.125,
        "unwrap": False,
        "limits": (-0.5, 0.5),
        "per_species": 2,
        "ranks": 2,
        "overrides": [
            "job/basename=cr_acc_qual_amr_boundary",
            "time/nlim=512",
            "time/tlim=0.64",
            "output1/dt=0.02",
            "output2/dt=0.64",
        ],
    },
    {
        "case": "mpi_decomposition",
        "title": "9. MPI decomposition crossing",
        "input": "cr_mpi_decomposition_invariance.athinput",
        "profile": "uniform",
        "b0": (0.0, 0.0, 1.0),
        "bgrad": 1.0,
        "bamp": 0.05,
        "bwave": 1.0,
        "plane": "xy",
        "fixed": 0.0,
        "unwrap": False,
        "limits": (-0.5, 0.5),
        "per_species": 2,
        "ranks": 2,
        "overrides": [
            "job/basename=cr_acc_qual_mpi_decomposition",
            "time/nlim=256",
            "time/tlim=1.0",
            "output1/dt=0.02",
            "output2/dt=1.0",
        ],
    },
    {
        "case": "isotropic_ensemble",
        "title": "10. Isotropic ensemble",
        "input": "cr_isotropic_ensemble.athinput",
        "profile": "uniform",
        "b0": (0.0, 0.0, 1.0),
        "bgrad": 1.0,
        "bamp": 0.05,
        "bwave": 1.0,
        "plane": "xy",
        "fixed": 0.0,
        "unwrap": False,
        "limits": (-0.5, 0.5),
        "per_species": 3,
        "ranks": 1,
        "overrides": [
            "job/basename=cr_acc_qual_isotropic",
            "time/nlim=128",
            "time/tlim=0.8",
            "output1/dt=0.02",
            "output2/dt=0.8",
            "output3/dt=0.8",
            "output4/dt=0.8",
            "output5/dt=0.8",
            "output6/dt=0.8",
        ],
    },
    {
        "case": "frozen_turbulent",
        "title": "11. Frozen multi-mode field",
        "input": "cr_frozen_turbulent_field.athinput",
        "profile": "turbulent",
        "b0": (0.0, 0.0, 1.0),
        "bgrad": 1.0,
        "bamp": 0.01,
        "bwave": 1.0,
        "plane": "xy",
        "fixed": 0.0,
        "unwrap": False,
        "limits": (-0.5, 0.5),
        "per_species": 3,
        "ranks": 1,
        "overrides": [
            "job/basename=cr_acc_qual_turbulent",
            "time/nlim=256",
            "time/tlim=1.28",
            "output1/dt=0.02",
            "output2/dt=1.28",
            "output3/dt=1.28",
        ],
    },
]


def study_qualitative_figures(args: argparse.Namespace) -> Dict:
    results = {}
    for config in QUALITATIVE_CASES:
        run_dir = _run(
            args, f"qualitative_{config['case']}", config["input"],
            config["overrides"], ranks=config["ranks"])
        tracks = _particle_tracks(run_dir, per_species=config["per_species"])
        if not tracks or any(len(points) < 2 for points in tracks.values()):
            raise ValueError(
                f"Qualitative case {config['case']} has an incomplete track")
        fig, ax = plt.subplots(figsize=(5.4, 4.5), constrained_layout=True)
        image = _draw_field_plane(ax, config, config["limits"])
        seen_species = set()
        for key, points in sorted(tracks.items()):
            projected = _project_track(points, config["plane"], config["unwrap"])
            indices = (1, 2) if config["plane"] == "xy" else (1, 3)
            raw_projected = np.asarray(points)[:, indices]
            color = plt.cm.tab10(key[0] % 10)
            label = f"species {key[0]}" if key[0] not in seen_species else None
            ax.plot(projected[:, 0], projected[:, 1], color=color,
                    linewidth=1.65, label=label)
            markers = projected if config["unwrap"] else raw_projected
            ax.scatter(markers[0, 0], markers[0, 1], color=color,
                       edgecolor="white", marker="o", s=27, zorder=4)
            ax.scatter(markers[-1, 0], markers[-1, 1], color=color,
                       edgecolor="white", marker="X", s=33, zorder=4)
            seen_species.add(key[0])
        final_time = max(points[-1][0] for points in tracks.values())
        ax.text(
            0.02, 0.98, f"t = 0 to {final_time:.2f}",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78})
        if config["profile"] in ("uniform", "linear_cross"):
            ax.text(
                0.02, 0.90, "Bz - B0z = 0",
                transform=ax.transAxes, va="top", fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "none",
                      "alpha": 0.78})
        if len(tracks) > 1:
            ax.legend(loc="lower right", fontsize=8, framealpha=0.85)
        if config["profile"] not in ("uniform", "linear_cross"):
            fig.colorbar(image, ax=ax, label="Bz - B0z")
        fig_path = args.output_dir / f"qualitative_{config['case']}.png"
        fig.savefig(fig_path, dpi=180)
        plt.close(fig)
        results[config["case"]] = {
            "figure": _artifact_path(fig_path),
            "input": config["input"],
            "plane": config["plane"],
            "fixed_coordinate": config["fixed"],
            "final_time": final_time,
            "track_count": len(tracks),
        }
    return results


def study_uniform_gyro(args: argparse.Namespace) -> Dict:
    cfl_values = [0.04, 0.02, 0.01]
    errors = []
    for cfl in cfl_values:
        run_dir = _run(
            args, f"uniform_gyro_cfl_{cfl}", "cr_uniform_gyro.athinput",
            [
                f"job/basename=cr_acc_gyro_cfl_{cfl}",
                f"particles/cfl_part={cfl}",
                "time/nlim=200",
                "time/tlim=0.32",
                "output1/dt=0.32",
                "output2/dt=0.32",
            ])
        summary = cr_tracer_inspect.summarize_restart(run_dir / "prst")
        particle = _single_particle(summary)
        vx, vy, _ = particle["reals"][3:6]
        errors.append(_phase_error(vx, vy, summary["restart"][0]["time"]))
    slope = _fit_log_slope(cfl_values, errors)

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.loglog(cfl_values, errors, marker="o")
    ax.set_xlabel("particle CFL")
    ax.set_ylabel("gyro phase error")
    ax.set_title(f"Uniform-B Boris phase convergence, slope={slope:.2f}")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig_path = args.output_dir / "uniform_gyro_phase_convergence.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "uniform_gyro",
        "cfl": cfl_values,
        "phase_error": errors,
        "slope": slope,
        "figure": _artifact_path(fig_path),
    }


def study_uniform_amr_mpi(args: argparse.Namespace) -> Dict:
    run_dir = _run(
        args, "uniform_amr_mpi", "cr_uniform_gyro_amr.athinput",
        [
            "job/basename=cr_acc_gyro_amr_mpi_docs",
            "particles/check_consistency_mode=full",
            "particles/validate_amr_lookup=true",
        ],
        ranks=2)
    output = (run_dir / "stdout.txt").read_text() + (run_dir / "stderr.txt").read_text()
    match = re.search(r"(\d+) MeshBlocks created, (\d+) deleted by AMR", output)
    created = int(match.group(1)) if match else 0
    deleted = int(match.group(2)) if match else 0
    expected = [512, 512]
    summary = cr_tracer_inspect.summarize_restart(run_dir / "prst")
    cr_tracer_inspect.validate_expected_counts(summary, sum(expected), expected)
    moments = cr_tracer_inspect.validate_moments(run_dir, expected)

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.bar(["created", "deleted"], [created, deleted])
    ax.set_ylabel("MeshBlocks")
    ax.set_title("Uniform-B AMR/MPI remap activity")
    fig.tight_layout()
    fig_path = args.output_dir / "uniform_amr_mpi_activity.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "uniform_amr_mpi",
        "meshblocks_created": created,
        "meshblocks_deleted": deleted,
        "rank_counts": summary["rank_counts"],
        "mean_speed2": [row["mean_speed2"] for row in moments],
        "figure": _artifact_path(fig_path),
    }


def study_linear_gather(args: argparse.Namespace) -> Dict:
    methods = ["lin", "trilinear", "tsc"]
    errors = []
    for method in methods:
        run_dir = _run(
            args, f"linear_gather_{method}", "cr_linear_gather.athinput",
            [
                f"job/basename=cr_acc_linear_{method}",
                f"particles/interpolation={method}",
                "time/nlim=1",
                "time/tlim=0.01",
            ])
        summary = cr_tracer_inspect.summarize_restart(run_dir / "prst")
        errors.append(_sample_error(
            summary, "linear_cross", (0.3, -0.2, 1.0), bgrad=0.7))

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.bar(methods, errors)
    ax.set_yscale("log")
    ax.set_ylabel("sampled-B error")
    ax.set_title("Linear-field gather error")
    fig.tight_layout()
    fig_path = args.output_dir / "linear_gather_error.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "linear_gather",
        "methods": methods,
        "sample_error": errors,
        "figure": _artifact_path(fig_path),
    }


def study_manufactured_gather(args: argparse.Namespace) -> Dict:
    nxs = [8, 16, 32]
    errors = []
    for nx in nxs:
        ppc = 1.0/(nx*nx*nx)
        run_dir = _run(
            args, f"manufactured_gather_{nx}",
            "cr_manufactured_divb_free.athinput",
            [
                f"job/basename=cr_acc_manufactured_{nx}",
                f"mesh/nx1={nx}",
                f"mesh/nx2={nx}",
                f"mesh/nx3={nx}",
                f"meshblock/nx1={nx}",
                f"meshblock/nx2={nx}",
                f"meshblock/nx3={nx}",
                f"particles/ppc={ppc:.17e}",
                "particles/interpolation=trilinear",
                "time/nlim=1",
                "time/tlim=0.01",
            ])
        summary = cr_tracer_inspect.summarize_restart(run_dir / "prst")
        errors.append(_sample_error(
            summary, "sinusoidal_divb_free", (0.1, -0.2, 1.0),
            bamp=0.025, bwave=1.0))
    slope = -_fit_log_slope(nxs, errors)

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.loglog(nxs, errors, marker="o")
    ax.set_xlabel("linear resolution")
    ax.set_ylabel("RMS sampled-B error")
    ax.set_title(f"Manufactured-field gather convergence, slope={slope:.2f}")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig_path = args.output_dir / "manufactured_gather_convergence.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "manufactured_gather",
        "nx": nxs,
        "sample_error": errors,
        "slope": slope,
        "figure": _artifact_path(fig_path),
    }


def study_smooth_orbit_reference(args: argparse.Namespace) -> Dict:
    resolutions = [16, 32, 64]
    states = []
    for nx in resolutions:
        ppc = 1.0/(nx*nx*nx)
        dt_fraction = 0.2/nx
        run_dir = _run(
            args, f"smooth_orbit_{nx}", "cr_smooth_orbit_reference.athinput",
            [
                f"job/basename=cr_acc_smooth_orbit_{nx}",
                f"mesh/nx1={nx}",
                f"mesh/nx2={nx}",
                f"mesh/nx3={nx}",
                f"meshblock/nx1={nx}",
                f"meshblock/nx2={nx}",
                f"meshblock/nx3={nx}",
                f"particles/ppc={ppc:.17e}",
                f"particles/cfl_part={dt_fraction}",
                f"particles/subcycle_gyro_fraction={dt_fraction}",
                "time/nlim=200",
                "time/tlim=0.1",
            ])
        states.append(_state(cr_tracer_inspect.summarize_restart(run_dir / "prst")))
    reference = states[-1]
    velocity_errors = [
        _vector_error(state["velocity"], reference["velocity"])
        for state in states[:-1]
    ]

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.loglog(resolutions[:-1], velocity_errors, marker="o")
    ax.set_xlabel("linear resolution")
    ax.set_ylabel("velocity error versus 64^3 reference")
    ax.set_title("Smooth-field orbit reference convergence")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig_path = args.output_dir / "smooth_orbit_reference_error.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "smooth_orbit_reference",
        "nx": resolutions,
        "velocity_error_to_reference": velocity_errors,
        "reference_nx": resolutions[-1],
        "figure": _artifact_path(fig_path),
    }


def study_magnetic_mirror(args: argparse.Namespace) -> Dict:
    run_dir = _run(
        args, "magnetic_mirror", "cr_magnetic_mirror.athinput",
        [
            "job/basename=cr_acc_mirror_docs",
            "time/nlim=100",
            "time/tlim=0.12",
        ])
    state = _state(cr_tracer_inspect.summarize_restart(run_dir / "prst"))
    initial_velocity = (0.7, 0.0, 0.25)
    initial_bfield = _exact_bfield(
        "mirror", 0.05, -0.03, -0.12, (0.0, 0.0, 1.0), bgrad=0.8)
    initial_moment = _magnetic_moment_proxy(initial_velocity, initial_bfield)
    final_moment = _magnetic_moment_proxy(state["velocity"], state["bfield"])
    relative_drift = abs(final_moment - initial_moment)/initial_moment

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.bar(["initial", "final"], [initial_moment, final_moment])
    ax.set_ylabel("magnetic moment proxy")
    ax.set_title(f"Mirror invariant drift = {relative_drift:.3e}")
    fig.tight_layout()
    fig_path = args.output_dir / "magnetic_mirror_moment.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "magnetic_mirror",
        "initial_moment": initial_moment,
        "final_moment": final_moment,
        "relative_moment_drift": relative_drift,
        "figure": _artifact_path(fig_path),
    }


def study_gradb_drift(args: argparse.Namespace) -> Dict:
    bgrads = [0.0, 0.5, 1.0]
    states = []
    for bgrad in bgrads:
        run_dir = _run(
            args, f"gradb_{bgrad}", "cr_gradb_curvature_drift.athinput",
            [
                f"job/basename=cr_acc_gradb_{bgrad}",
                f"problem/Bgrad={bgrad}",
                "particles/cfl_part=0.01",
                "particles/subcycle_gyro_fraction=0.01",
                "time/nlim=1000",
                "time/tlim=0.64",
            ])
        states.append(_state(cr_tracer_inspect.summarize_restart(run_dir / "prst")))
    velocity_delta = [
        _vector_error(state["velocity"], states[0]["velocity"])
        for state in states
    ]

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(bgrads, velocity_delta, marker="o")
    ax.set_xlabel("B gradient")
    ax.set_ylabel("velocity delta versus uniform field")
    ax.set_title("Grad-B response")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig_path = args.output_dir / "gradb_drift_response.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "gradb_drift",
        "bgrad": bgrads,
        "velocity_delta": velocity_delta,
        "figure": _artifact_path(fig_path),
    }


def study_amr_boundary(args: argparse.Namespace) -> Dict:
    run_dir = _run(
        args, "amr_boundary", "cr_amr_boundary_convergence.athinput",
        [
            "job/basename=cr_acc_amr_boundary_docs",
            "particles/check_consistency_mode=full",
            "particles/validate_amr_lookup=true",
        ],
        ranks=2)
    output = (run_dir / "stdout.txt").read_text() + (run_dir / "stderr.txt").read_text()
    match = re.search(r"(\d+) MeshBlocks created, (\d+) deleted by AMR", output)
    created = int(match.group(1)) if match else 0
    deleted = int(match.group(2)) if match else 0
    expected = [512, 512]
    summary = cr_tracer_inspect.summarize_restart(run_dir / "prst")
    cr_tracer_inspect.validate_expected_counts(summary, sum(expected), expected)
    moments = cr_tracer_inspect.validate_moments(run_dir, expected)

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.bar(["created", "deleted"], [created, deleted])
    ax.set_ylabel("MeshBlocks")
    ax.set_title("Smooth-field AMR boundary stress")
    fig.tight_layout()
    fig_path = args.output_dir / "amr_boundary_activity.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "amr_boundary",
        "meshblocks_created": created,
        "meshblocks_deleted": deleted,
        "rank_counts": summary["rank_counts"],
        "mean_speed2": [row["mean_speed2"] for row in moments],
        "figure": _artifact_path(fig_path),
    }


def study_mpi_decomposition(args: argparse.Namespace) -> Dict:
    run_dirs = []
    for ranks in [1, 2]:
        run_dirs.append(_run(
            args, f"mpi_decomposition_n{ranks}",
            "cr_mpi_decomposition_invariance.athinput",
            [f"job/basename=cr_acc_mpi_decomposition_n{ranks}"],
            ranks=ranks))
    expected = [2, 2]
    summaries = []
    moments = []
    for run_dir in run_dirs:
        summary = cr_tracer_inspect.summarize_restart(run_dir / "prst")
        cr_tracer_inspect.validate_expected_counts(summary, sum(expected), expected)
        summaries.append(summary)
        moments.append(cr_tracer_inspect.validate_moments(run_dir, expected))
    speed2_delta = [
        abs(moments[0][species]["mean_speed2"] -
            moments[1][species]["mean_speed2"])
        for species in range(2)
    ]
    particles = []
    for summary in summaries:
        particle_map = {}
        for restart in summary["restart"]:
            for particle in restart["particles"]:
                particle_map[(particle["species"], particle["tag"])] = particle
        particles.append(particle_map)
    per_tag_position_delta = []
    per_tag_velocity_delta = []
    for key in sorted(particles[0]):
        p0 = particles[0][key]
        p1 = particles[1][key]
        per_tag_position_delta.append(
            _vector_error(p0["reals"][0:3], p1["reals"][0:3]))
        per_tag_velocity_delta.append(
            _vector_error(p0["reals"][3:6], p1["reals"][3:6]))

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.bar(["position", "velocity"],
           [max(per_tag_position_delta), max(per_tag_velocity_delta)])
    ax.set_ylabel("max one-rank/two-rank per-tag delta")
    ax.set_title("MPI decomposition per-tag agreement")
    fig.tight_layout()
    fig_path = args.output_dir / "mpi_decomposition_delta.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "mpi_decomposition",
        "speed2_delta": speed2_delta,
        "max_position_delta": max(per_tag_position_delta),
        "max_velocity_delta": max(per_tag_velocity_delta),
        "figure": _artifact_path(fig_path),
    }


def study_isotropic_ensemble(args: argparse.Namespace) -> Dict:
    run_dir = _run(
        args, "isotropic_ensemble", "cr_isotropic_ensemble.athinput",
        [
            "job/basename=cr_acc_isotropic",
            "time/nlim=1",
            "time/tlim=0.02",
        ])
    expected = [1024, 1024]
    summary = cr_tracer_inspect.summarize_restart(run_dir / "prst")
    cr_tracer_inspect.validate_expected_counts(summary, sum(expected), expected)
    cr_tracer_inspect.validate_histograms(run_dir, expected, 16)
    moments = cr_tracer_inspect.validate_moments(run_dir, expected)
    rows = [{"species": row["species"], "mean_mu": row["mean_mu"],
             "mean_mu2": row["mean_mu2"]} for row in moments]

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    species = [row["species"] for row in rows]
    ax.bar([s - 0.18 for s in species], [row["mean_mu"] for row in rows],
           width=0.36, label="<mu>")
    ax.bar([s + 0.18 for s in species], [row["mean_mu2"] for row in rows],
           width=0.36, label="<mu^2>")
    ax.axhline(1.0/3.0, color="black", linestyle="--", linewidth=1.0)
    ax.set_xlabel("species")
    ax.set_ylabel("moment")
    ax.set_title("Isotropic ensemble pitch-angle moments")
    ax.legend()
    fig.tight_layout()
    fig_path = args.output_dir / "isotropic_ensemble_moments.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "isotropic_ensemble",
        "moments": rows,
        "figure": _artifact_path(fig_path),
    }


def study_frozen_turbulent(args: argparse.Namespace) -> Dict:
    run_dir = _run(
        args, "frozen_turbulent", "cr_frozen_turbulent_field.athinput",
        [
            "job/basename=cr_acc_turbulent_docs",
            "time/nlim=4",
            "time/tlim=0.08",
        ])
    expected = [128, 128]
    summary = cr_tracer_inspect.summarize_restart(run_dir / "prst")
    cr_tracer_inspect.validate_expected_counts(summary, sum(expected), expected)
    moments = cr_tracer_inspect.validate_moments(run_dir, expected)
    cr_tracer_inspect.validate_joint_spectra(run_dir, expected, 16, 16)
    rows = [{"species": row["species"], "mean_mu": row["mean_mu"],
             "mean_mu2": row["mean_mu2"]} for row in moments]

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    species = [row["species"] for row in rows]
    ax.bar([s - 0.18 for s in species], [row["mean_mu"] for row in rows],
           width=0.36, label="<mu>")
    ax.bar([s + 0.18 for s in species], [row["mean_mu2"] for row in rows],
           width=0.36, label="<mu^2>")
    ax.set_xlabel("species")
    ax.set_ylabel("moment")
    ax.set_title("Frozen turbulent-field pitch-angle moments")
    ax.legend()
    fig.tight_layout()
    fig_path = args.output_dir / "frozen_turbulent_moments.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "frozen_turbulent",
        "moments": rows,
        "figure": _artifact_path(fig_path),
    }


STUDIES = {
    "uniform_gyro": study_uniform_gyro,
    "uniform_amr_mpi": study_uniform_amr_mpi,
    "linear_gather": study_linear_gather,
    "manufactured_gather": study_manufactured_gather,
    "smooth_orbit_reference": study_smooth_orbit_reference,
    "magnetic_mirror": study_magnetic_mirror,
    "gradb_drift": study_gradb_drift,
    "amr_boundary": study_amr_boundary,
    "mpi_decomposition": study_mpi_decomposition,
    "isotropic_ensemble": study_isotropic_ensemble,
    "frozen_turbulent": study_frozen_turbulent,
}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run CR tracer accuracy studies and write docs artifacts.")
    parser.add_argument("--athena", type=Path, default=_default_athena())
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--run-root", type=Path,
                        default=Path("/tmp/cr_tracer_accuracy_runs"),
                        help="Temporary directory for raw Athena run outputs")
    parser.add_argument("--cases", default=",".join(STUDIES.keys()),
                        help="Comma-separated studies to run")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.run_root.mkdir(parents=True, exist_ok=True)

    selected = [case.strip() for case in args.cases.split(",") if case.strip()]
    results = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "athena": str(args.athena),
        "cases": [],
    }
    for case in selected:
        if case not in STUDIES:
            raise ValueError(f"Unknown accuracy study '{case}'")
        results["cases"].append(STUDIES[case](args))
    results["qualitative_figures"] = study_qualitative_figures(args)

    json_path = args.output_dir / "cr_tracer_accuracy_summary.json"
    json_path.write_text(json.dumps(results, indent=2) + "\n")
    md_path = args.output_dir / "cr_tracer_accuracy_results.md"
    lines = ["---", "orphan: true", "---", "", "# CR Tracer Accuracy Results", ""]
    for case in results["cases"]:
        lines.append(f"## {case['case']}")
        if "slope" in case:
            lines.append(f"- fitted slope: `{case['slope']:.3f}`")
        lines.append(f"- figure: `{case['figure']}`")
        lines.append("")
    lines.append("## qualitative_figures")
    for case, result in results["qualitative_figures"].items():
        lines.append(
            f"- {case}: `{result['figure']}` "
            f"(`t = 0` to `{result['final_time']:.2f}`, "
            f"`{result['track_count']}` tracks)")
    md_path.write_text("\n".join(lines).rstrip() + "\n")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
