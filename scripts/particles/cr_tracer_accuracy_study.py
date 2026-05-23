#!/usr/bin/env python3
"""Run CR tracer accuracy studies and generate documentation figures."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import time
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
    started = time.perf_counter()
    proc = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    wall_seconds = time.perf_counter() - started
    if proc.returncode != 0:
        raise RuntimeError(
            f"{case} failed\nCOMMAND: {' '.join(command)}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    (run_dir / "stdout.txt").write_text(proc.stdout)
    (run_dir / "stderr.txt").write_text(proc.stderr)
    args.executions.append({
        "case": case,
        "input": input_name,
        "ranks": ranks,
        "overrides": list(overrides),
        "wall_seconds": wall_seconds,
        "command": command,
    })
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


def _particle_map(summary: Dict) -> Dict:
    particles = {}
    for restart in summary["restart"]:
        for particle in restart["particles"]:
            key = (particle["species"], particle["tag"])
            if key in particles:
                raise ValueError(f"Duplicate particle key {key}")
            particles[key] = particle
    return particles


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


def _pitch_angle(particle: Dict) -> float:
    velocity = particle["reals"][3:6]
    bfield = particle["reals"][7:10]
    denominator = _norm3(velocity)*_norm3(bfield)
    return sum(velocity[i]*bfield[i] for i in range(3))/denominator


def _state_pitch_angle(state: Dict) -> float:
    denominator = _norm3(state["velocity"])*_norm3(state["bfield"])
    return sum(state["velocity"][i]*state["bfield"][i] for i in range(3))/denominator


def _state_magnetic_moment(state: Dict) -> float:
    return _magnetic_moment_proxy(state["velocity"], state["bfield"])


def _ensemble_mu_moments(particles: Dict) -> Dict:
    values = np.asarray([_pitch_angle(particle) for particle in particles.values()])
    return {
        "mean_mu": float(np.mean(values)),
        "mean_mu2": float(np.mean(values*values)),
    }


def _ensemble_max_speed2_error(sequence: Sequence[Dict], expected: float = 1.0) -> float:
    return max(
        abs(sum(value*value for value in particle["reals"][3:6]) - expected)
        for sample in sequence for particle in sample["particles"].values())


def _restart_sequence(run_dir: Path) -> List[Dict]:
    snapshots = [
        cr_tracer_inspect.read_prst_file(path)
        for path in (run_dir / "prst").rglob("*.prst")
    ]
    times = sorted({snapshot["time"] for snapshot in snapshots})
    sequence = []
    for sample_time in times:
        particles = {}
        for snapshot in snapshots:
            if snapshot["time"] != sample_time:
                continue
            for particle in snapshot["particles"]:
                key = (particle["species"], particle["tag"])
                if key in particles:
                    raise ValueError(f"Duplicate particle key {key}")
                particles[key] = particle
        sequence.append({"time": sample_time, "particles": particles})
    return sequence


def _pitch_angle_correlation(initial: Dict, current: Dict) -> float:
    if set(initial) != set(current):
        raise ValueError("Pitch-angle correlation requires matching particles")
    numerator = sum(_pitch_angle(initial[key])*_pitch_angle(current[key])
                    for key in initial)
    denominator = sum(_pitch_angle(initial[key])**2 for key in initial)
    return numerator/denominator


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


def _sample_errors(summary: Dict, profile: str, b0: Sequence[float],
                   bgrad: float = 1.0, bamp: float = 0.05,
                   bwave: float = 1.0) -> Dict:
    errors = []
    for restart in summary["restart"]:
        for particle in restart["particles"]:
            x, y, z = particle["reals"][0:3]
            sampled = particle["reals"][7:10]
            exact = _exact_bfield(profile, x, y, z, b0, bgrad, bamp, bwave)
            errors.append(_vector_error(sampled, exact))
    return {
        "l1": sum(errors)/len(errors),
        "l2": math.sqrt(sum(error*error for error in errors)/len(errors)),
        "linf": max(errors),
    }


def _sample_error(summary: Dict, profile: str, b0: Sequence[float],
                  bgrad: float = 1.0, bamp: float = 0.05,
                  bwave: float = 1.0) -> float:
    return _sample_errors(summary, profile, b0, bgrad, bamp, bwave)["l2"]


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
        "docs_mpi": True,
        "overrides": [
            "job/basename=cr_acc_qual_uniform_amr_mpi",
            "time/nlim=512",
            "time/tlim=0.64",
            "output1/dt=0.08",
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
        "bgrad": 4.0,
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
            "time/nlim=200000",
            "time/tlim=6.0",
            "output1/dt=0.04",
            "output2/dt=6.0",
        ],
        "docs_overrides": [
            "time/tlim=18.0",
            "output2/dt=18.0",
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
        "docs_mpi": True,
        "overrides": [
            "job/basename=cr_acc_qual_amr_boundary",
            "time/nlim=512",
            "time/tlim=0.64",
            "output1/dt=0.08",
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
        "docs_mpi": True,
        "overrides": [
            "job/basename=cr_acc_qual_mpi_decomposition",
            "time/nlim=256",
            "time/tlim=1.0",
            "output1/dt=0.08",
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
        "per_species": 1,
        "ranks": 1,
        "docs_mpi": True,
        "overrides": [
            "job/basename=cr_acc_qual_turbulent",
            "time/nlim=2000",
            "time/tlim=1.28",
            "output1/dt=0.08",
            "output2/dt=1.28",
            "output3/dt=1.28",
        ],
        "docs_overrides": [
            "time/tlim=3.20",
            "output1/dt=0.16",
            "output2/dt=3.20",
            "output3/dt=3.20",
        ],
    },
    {
        "case": "pitch_angle_decorrelation",
        "title": "12. Deterministic pitch-angle evolution",
        "input": "cr_pitch_angle_decorrelation.athinput",
        "profile": "turbulent",
        "b0": (0.0, 0.0, 1.0),
        "bgrad": 1.0,
        "bamp": 0.04,
        "bwave": 1.0,
        "plane": "xy",
        "fixed": 0.0,
        "unwrap": False,
        "limits": (-0.5, 0.5),
        "per_species": 3,
        "ranks": 2,
        "docs_mpi": True,
        "overrides": [
            "job/basename=cr_acc_qual_pitch_angle",
            "time/nlim=2000",
            "time/tlim=1.28",
            "output1/dt=0.08",
            "output4/dt=1.28",
        ],
    },
]


def study_qualitative_figures(args: argparse.Namespace) -> Dict:
    results = {}
    for config in QUALITATIVE_CASES:
        overrides = list(config["overrides"])
        if args.profile == "docs":
            overrides += config.get("docs_overrides", [])
        ranks = config["ranks"]
        if args.profile == "docs" and config.get("docs_mpi"):
            ranks = args.mpi_ranks
        run_dir = _run(
            args, f"qualitative_{config['case']}", config["input"],
            overrides, ranks=ranks)
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
            "ranks": ranks,
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
    ranks = args.mpi_ranks if args.profile == "docs" else min(2, args.mpi_ranks)
    run_dir = _run(
        args, "uniform_amr_mpi", "cr_uniform_gyro_amr.athinput",
        [
            "job/basename=cr_acc_gyro_amr_mpi_docs",
            "particles/check_consistency_mode=full",
            "particles/validate_amr_lookup=true",
        ],
        ranks=ranks)
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
        "ranks": ranks,
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
    nxs = [16, 32, 64, 128] if args.profile == "docs" else [8, 16, 32]
    error_rows = []
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
        error_rows.append(_sample_errors(
            summary, "sinusoidal_divb_free", (0.1, -0.2, 1.0),
            bamp=0.025, bwave=1.0))
    slopes = {
        key: -_fit_log_slope(nxs, [row[key] for row in error_rows])
        for key in ("l1", "l2", "linf")
    }

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for key, label in (("l1", "L1"), ("l2", "L2"), ("linf", "Linf")):
        ax.loglog(nxs, [row[key] for row in error_rows], marker="o",
                  label=f"{label}, slope={slopes[key]:.2f}")
    ax.set_xlabel("linear resolution")
    ax.set_ylabel("sampled-B error")
    ax.set_title("Manufactured-field gather convergence")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig_path = args.output_dir / "manufactured_gather_convergence.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "manufactured_gather",
        "nx": nxs,
        "sample_error": error_rows,
        "slope": slopes["l2"],
        "slopes": slopes,
        "figure": _artifact_path(fig_path),
    }


def study_smooth_orbit_reference(args: argparse.Namespace) -> Dict:
    resolutions = [16, 32, 64, 128] if args.profile == "docs" else [16, 32, 64]
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
    position_errors = [
        _vector_error(state["position"], reference["position"])
        for state in states[:-1]
    ]
    velocity_errors = [
        _vector_error(state["velocity"], reference["velocity"])
        for state in states[:-1]
    ]
    pitch_angle_errors = [
        abs(_state_pitch_angle(state) - _state_pitch_angle(reference))
        for state in states[:-1]
    ]
    magnetic_moment_errors = [
        abs(_state_magnetic_moment(state) - _state_magnetic_moment(reference))
        for state in states[:-1]
    ]

    temporal = None
    if args.profile == "docs":
        temporal_nx = 64
        temporal_cfl = [0.0125, 0.00625, 0.003125, 0.0015625]
        temporal_states = []
        for cfl in temporal_cfl:
            if math.isclose(cfl, 0.2/temporal_nx):
                temporal_states.append(states[resolutions.index(temporal_nx)])
                continue
            run_dir = _run(
                args, f"smooth_orbit_timestep_{cfl}",
                "cr_smooth_orbit_reference.athinput",
                [
                    f"job/basename=cr_acc_smooth_timestep_{cfl}",
                    f"mesh/nx1={temporal_nx}",
                    f"mesh/nx2={temporal_nx}",
                    f"mesh/nx3={temporal_nx}",
                    f"meshblock/nx1={temporal_nx}",
                    f"meshblock/nx2={temporal_nx}",
                    f"meshblock/nx3={temporal_nx}",
                    f"particles/ppc={1.0/(temporal_nx**3):.17e}",
                    f"particles/cfl_part={cfl}",
                    f"particles/subcycle_gyro_fraction={cfl}",
                    "time/nlim=200",
                    "time/tlim=0.1",
                ])
            temporal_states.append(
                _state(cr_tracer_inspect.summarize_restart(run_dir / "prst")))
        temporal_reference = temporal_states[-1]
        temporal = {
            "nx": temporal_nx,
            "cfl": temporal_cfl,
            "position_error_to_reference": [
                _vector_error(state["position"], temporal_reference["position"])
                for state in temporal_states[:-1]
            ],
            "velocity_error_to_reference": [
                _vector_error(state["velocity"], temporal_reference["velocity"])
                for state in temporal_states[:-1]
            ],
            "pitch_angle_error_to_reference": [
                abs(_state_pitch_angle(state) -
                    _state_pitch_angle(temporal_reference))
                for state in temporal_states[:-1]
            ],
            "magnetic_moment_error_to_reference": [
                abs(_state_magnetic_moment(state) -
                    _state_magnetic_moment(temporal_reference))
                for state in temporal_states[:-1]
            ],
        }

    columns = 2 if temporal else 1
    fig, axes = plt.subplots(1, columns, figsize=(5.4*columns, 3.8),
                             constrained_layout=True)
    axes = np.atleast_1d(axes)
    for values, label in (
        (position_errors, "position"),
        (velocity_errors, "velocity"),
        (pitch_angle_errors, "pitch angle"),
        (magnetic_moment_errors, "magnetic moment"),
    ):
        axes[0].loglog(resolutions[:-1], values, marker="o", label=label)
    axes[0].set_xlabel("linear resolution")
    axes[0].set_ylabel(f"error versus {resolutions[-1]}^3 reference")
    axes[0].set_title("Spatial convergence")
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=8)
    if temporal:
        for key, label in (
            ("position_error_to_reference", "position"),
            ("velocity_error_to_reference", "velocity"),
            ("pitch_angle_error_to_reference", "pitch angle"),
            ("magnetic_moment_error_to_reference", "magnetic moment"),
        ):
            axes[1].loglog(temporal["cfl"][:-1], temporal[key],
                           marker="o", label=label)
        axes[1].set_xlabel("particle CFL")
        axes[1].set_ylabel("error versus smallest-step reference")
        axes[1].set_title(f"Temporal convergence at {temporal['nx']}^3")
        axes[1].grid(True, which="both", alpha=0.3)
        axes[1].legend(fontsize=8)
    fig_path = args.output_dir / "smooth_orbit_reference_error.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    result = {
        "case": "smooth_orbit_reference",
        "nx": resolutions,
        "position_error_to_reference": position_errors,
        "velocity_error_to_reference": velocity_errors,
        "pitch_angle_error_to_reference": pitch_angle_errors,
        "magnetic_moment_error_to_reference": magnetic_moment_errors,
        "reference_nx": resolutions[-1],
        "requested_reference_nx": 256,
        "reference_status": "local reference; 256^3 reference pending",
        "figure": _artifact_path(fig_path),
    }
    if temporal:
        result["temporal_sweep"] = temporal
    return result


def study_magnetic_mirror(args: argparse.Namespace) -> Dict:
    final_time = 18.0 if args.profile == "docs" else 6.0
    run_dir = _run(
        args, "magnetic_mirror", "cr_magnetic_mirror.athinput",
        [
            "job/basename=cr_acc_mirror_docs",
            "time/nlim=200000",
            f"time/tlim={final_time}",
            "output1/dt=0.04",
            f"output2/dt={final_time}",
        ])
    sequence = _restart_sequence(run_dir)
    particles = [next(iter(row["particles"].values())) for row in sequence]
    times = np.asarray([row["time"] for row in sequence])
    z_values = np.asarray([particle["reals"][2] for particle in particles])
    vpar = np.asarray([
        sum(particle["reals"][3 + i]*particle["reals"][7 + i]
            for i in range(3))/_norm3(particle["reals"][7:10])
        for particle in particles
    ])
    moments = np.asarray([
        _magnetic_moment_proxy(particle["reals"][3:6], particle["reals"][7:10])
        for particle in particles
    ])
    reversal_indices = np.flatnonzero(vpar[1:]*vpar[:-1] < 0.0) + 1
    turning_times = times[reversal_indices]
    turning_locations = z_values[reversal_indices]
    periods = np.diff(turning_times[::2])
    final_relative_drift = abs(moments[-1] - moments[0])/moments[0]
    max_relative_excursion = float(np.max(np.abs(moments/moments[0] - 1.0)))

    passing_dir = _run(
        args, "magnetic_mirror_passing_control", "cr_magnetic_mirror.athinput",
        [
            "job/basename=cr_acc_mirror_passing",
            "problem/v0z=0.80",
            "time/nlim=200000",
            "time/tlim=6.0",
            "output1/dt=0.04",
            "output2/dt=6.0",
        ])
    passing_sequence = _restart_sequence(passing_dir)
    passing_particles = [
        next(iter(row["particles"].values())) for row in passing_sequence
    ]
    passing_vpar = np.asarray([
        sum(particle["reals"][3 + i]*particle["reals"][7 + i]
            for i in range(3))/_norm3(particle["reals"][7:10])
        for particle in passing_particles
    ])
    passing_reversals = int(np.count_nonzero(passing_vpar[1:]*passing_vpar[:-1] < 0.0))
    resolutions = [16, 32, 64, 128] if args.profile == "docs" else [16, 32, 64]
    resolution_states = []
    for nx in resolutions:
        run = _run(
            args, f"magnetic_mirror_convergence_{nx}",
            "cr_magnetic_mirror.athinput",
            [
                f"job/basename=cr_acc_mirror_convergence_{nx}",
                f"mesh/nx1={nx}",
                f"mesh/nx2={nx}",
                f"mesh/nx3={nx}",
                f"meshblock/nx1={nx}",
                f"meshblock/nx2={nx}",
                f"meshblock/nx3={nx}",
                f"particles/ppc={1.0/(nx**3):.17e}",
                "time/nlim=100",
                "time/tlim=0.12",
                "output1/dt=0.12",
                "output2/dt=0.12",
            ])
        resolution_states.append(
            _state(cr_tracer_inspect.summarize_restart(run / "prst")))
    convergence_reference = resolution_states[-1]
    position_errors = [
        _vector_error(row["position"], convergence_reference["position"])
        for row in resolution_states[:-1]
    ]
    velocity_errors = [
        _vector_error(row["velocity"], convergence_reference["velocity"])
        for row in resolution_states[:-1]
    ]

    fig, axes = plt.subplots(2, 1, figsize=(6.2, 5.0), sharex=True,
                             constrained_layout=True)
    axes[0].plot(times, z_values)
    axes[0].scatter(turning_times, turning_locations, color="tab:red", s=15)
    axes[0].set_ylabel("z")
    axes[0].set_title("Trapped magnetic-mirror trajectory")
    axes[1].plot(times, moments/moments[0] - 1.0)
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("relative moment change")
    axes[1].grid(True, alpha=0.3)
    fig_path = args.output_dir / "magnetic_mirror_moment.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.loglog(resolutions[:-1], position_errors, marker="o", label="position")
    ax.loglog(resolutions[:-1], velocity_errors, marker="o", label="velocity")
    ax.set_xlabel("linear resolution")
    ax.set_ylabel(f"error versus {resolutions[-1]}^3 reference")
    ax.set_title("Mirror short-interval spatial convergence")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    convergence_path = args.output_dir / "magnetic_mirror_convergence.png"
    fig.savefig(convergence_path, dpi=180)
    plt.close(fig)
    return {
        "case": "magnetic_mirror",
        "final_time": final_time,
        "z_range": [float(np.min(z_values)), float(np.max(z_values))],
        "parallel_velocity_reversals": int(len(reversal_indices)),
        "turning_locations": turning_locations.tolist(),
        "mean_bounce_period": float(np.mean(periods)) if len(periods) else None,
        "initial_moment": float(moments[0]),
        "final_moment": float(moments[-1]),
        "relative_moment_drift": float(final_relative_drift),
        "max_relative_moment_excursion": max_relative_excursion,
        "passing_control_reversals": passing_reversals,
        "nx": resolutions,
        "position_error_to_reference": position_errors,
        "velocity_error_to_reference": velocity_errors,
        "reference_nx": resolutions[-1],
        "requested_reference_nx": 256,
        "reference_status": "local reference; 256^3 reference pending",
        "figure": _artifact_path(fig_path),
        "convergence_figure": _artifact_path(convergence_path),
    }


def study_gradb_drift(args: argparse.Namespace) -> Dict:
    resolutions = [16, 32, 64, 128] if args.profile == "docs" else [16, 32, 64]
    states = []
    for nx in resolutions:
        run_dir = _run(
            args, f"gradb_{nx}", "cr_gradb_curvature_drift.athinput",
            [
                f"job/basename=cr_acc_gradb_{nx}",
                f"mesh/nx1={nx}",
                f"mesh/nx2={nx}",
                f"mesh/nx3={nx}",
                f"meshblock/nx1={nx}",
                f"meshblock/nx2={nx}",
                f"meshblock/nx3={nx}",
                f"particles/ppc={1.0/(nx**3):.17e}",
                "problem/Bgrad=1.0",
                "particles/cfl_part=0.01",
                "particles/subcycle_gyro_fraction=0.01",
                "time/nlim=1000",
                "time/tlim=0.64",
            ])
        states.append(_state(cr_tracer_inspect.summarize_restart(run_dir / "prst")))
    reference = states[-1]
    position_errors = [
        _vector_error(state["position"], reference["position"])
        for state in states[:-1]
    ]
    velocity_errors = [
        _vector_error(state["velocity"], reference["velocity"])
        for state in states[:-1]
    ]
    transverse_displacement = [state["position"][1] for state in states]

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.loglog(resolutions[:-1], position_errors, marker="o", label="position")
    ax.loglog(resolutions[:-1], velocity_errors, marker="o", label="velocity")
    ax.set_xlabel("linear resolution")
    ax.set_ylabel(f"error versus {resolutions[-1]}^3 reference")
    ax.set_title("Grad-B orbit convergence")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig_path = args.output_dir / "gradb_drift_response.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "gradb_drift",
        "nx": resolutions,
        "position_error_to_reference": position_errors,
        "velocity_error_to_reference": velocity_errors,
        "transverse_position": transverse_displacement,
        "reference_nx": resolutions[-1],
        "requested_reference_nx": 256,
        "reference_status": "local reference; 256^3 reference pending",
        "guiding_center_validity_ratio": 0.8/(0.9*0.9),
        "guiding_center_status": (
            "rho_L/L_B is order unity; use the full-orbit reference, "
            "not an asymptotic guiding-center estimate"),
        "figure": _artifact_path(fig_path),
    }


def study_amr_boundary(args: argparse.Namespace) -> Dict:
    ranks = args.mpi_ranks if args.profile == "docs" else min(2, args.mpi_ranks)
    base_resolutions = [16, 32, 64] if args.profile == "docs" else [16]
    run_dirs = []
    summaries = []
    expected = [512, 512]
    for nx in base_resolutions:
        run_dir = _run(
            args, f"amr_boundary_{nx}", "cr_amr_boundary_convergence.athinput",
            [
                f"job/basename=cr_acc_amr_boundary_{nx}",
                f"mesh/nx1={nx}",
                f"mesh/nx2={nx}",
                f"mesh/nx3={nx}",
                "meshblock/nx1=8",
                "meshblock/nx2=8",
                "meshblock/nx3=8",
                f"particles/ppc={512.0/(nx**3):.17e}",
                "mesh_refinement/max_nmb_per_rank=512",
                "particles/check_consistency_mode=full",
                "particles/validate_amr_lookup=true",
            ],
            ranks=ranks)
        summary = cr_tracer_inspect.summarize_restart(run_dir / "prst")
        cr_tracer_inspect.validate_expected_counts(summary, sum(expected), expected)
        run_dirs.append(run_dir)
        summaries.append(summary)
    run_dir = run_dirs[-1]
    summary = summaries[-1]
    output = (run_dir / "stdout.txt").read_text() + (run_dir / "stderr.txt").read_text()
    match = re.search(r"(\d+) MeshBlocks created, (\d+) deleted by AMR", output)
    created = int(match.group(1)) if match else 0
    deleted = int(match.group(2)) if match else 0
    moments = cr_tracer_inspect.validate_moments(run_dir, expected)
    uniform_grid_delta = None
    if args.profile == "docs":
        uniform_dir = _run(
            args, "amr_boundary_uniform_128",
            "cr_amr_boundary_convergence.athinput",
            [
                "job/basename=cr_acc_amr_boundary_uniform_128",
                "mesh/nx1=128",
                "mesh/nx2=128",
                "mesh/nx3=128",
                "meshblock/nx1=16",
                "meshblock/nx2=16",
                "meshblock/nx3=16",
                "mesh_refinement/refinement=none",
                f"particles/ppc={512.0/(128**3):.17e}",
                "mesh_refinement/max_nmb_per_rank=512",
                "particles/check_consistency_mode=full",
                "particles/validate_amr_lookup=true",
            ],
            ranks=ranks)
        uniform_summary = cr_tracer_inspect.summarize_restart(uniform_dir / "prst")
        cr_tracer_inspect.validate_expected_counts(
            uniform_summary, sum(expected), expected)
        amr_particles = _particle_map(summary)
        uniform_particles = _particle_map(uniform_summary)
        uniform_grid_delta = {
            "position": max(
                _vector_error(amr_particles[key]["reals"][0:3],
                              uniform_particles[key]["reals"][0:3])
                for key in amr_particles),
            "velocity": max(
                _vector_error(amr_particles[key]["reals"][3:6],
                              uniform_particles[key]["reals"][3:6])
                for key in amr_particles),
        }
    errors = []
    if len(summaries) > 1:
        reference = _particle_map(summaries[-1])
        for candidate in summaries[:-1]:
            particles = _particle_map(candidate)
            errors.append(max(
                _vector_error(particles[key]["reals"][0:3],
                              reference[key]["reals"][0:3])
                for key in reference))

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    if errors:
        ax.loglog([2*nx for nx in base_resolutions[:-1]], errors, marker="o")
        ax.set_xlabel("effective finest linear resolution")
        ax.set_ylabel(f"max position delta versus {2*base_resolutions[-1]}^3")
        ax.set_title("AMR boundary trajectory convergence")
        ax.grid(True, which="both", alpha=0.3)
    else:
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
        "ranks": ranks,
        "effective_nx": [2*nx for nx in base_resolutions],
        "max_position_delta_to_reference": errors,
        "uniform_grid_128_delta": uniform_grid_delta,
        "reference_nx": 2*base_resolutions[-1],
        "figure": _artifact_path(fig_path),
    }


def study_mpi_decomposition(args: argparse.Namespace) -> Dict:
    rank_values = [1, 2, 4, 8] if args.profile == "docs" else [1, 2]
    run_dirs = []
    for ranks in rank_values:
        run_dirs.append(_run(
            args, f"mpi_decomposition_n{ranks}",
            "cr_mpi_decomposition_invariance.athinput",
            [
                f"job/basename=cr_acc_mpi_decomposition_n{ranks}",
                "time/nlim=256",
                "time/tlim=0.40",
                "output1/dt=0.20",
                "output2/dt=0.20",
            ],
            ranks=ranks))
    expected = [8, 8]
    summaries = []
    moments = []
    for run_dir in run_dirs:
        summary = cr_tracer_inspect.summarize_restart(run_dir / "prst")
        cr_tracer_inspect.validate_expected_counts(summary, sum(expected), expected)
        summaries.append(summary)
        moments.append(cr_tracer_inspect.validate_moments(run_dir, expected))
    particles = [_particle_map(summary) for summary in summaries]
    speed2_delta = []
    max_position_delta = []
    max_velocity_delta = []
    for candidate_index in range(1, len(rank_values)):
        speed2_delta.append([
            abs(moments[0][species]["mean_speed2"] -
                moments[candidate_index][species]["mean_speed2"])
            for species in range(2)
        ])
        max_position_delta.append(max(
            _vector_error(particles[0][key]["reals"][0:3],
                          particles[candidate_index][key]["reals"][0:3])
            for key in particles[0]))
        max_velocity_delta.append(max(
            _vector_error(particles[0][key]["reals"][3:6],
                          particles[candidate_index][key]["reals"][3:6])
            for key in particles[0]))

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    xvalues = np.arange(len(rank_values) - 1)
    ax.bar(xvalues - 0.18, max_position_delta, width=0.36, label="position")
    ax.bar(xvalues + 0.18, max_velocity_delta, width=0.36, label="velocity")
    ax.set_xticks(xvalues, [str(rank) for rank in rank_values[1:]])
    ax.set_xlabel("MPI ranks compared with one-rank run")
    ax.set_ylabel("maximum per-tag delta")
    ax.set_title("MPI decomposition per-tag agreement")
    ax.legend()
    fig.tight_layout()
    fig_path = args.output_dir / "mpi_decomposition_delta.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "mpi_decomposition",
        "ranks": rank_values,
        "speed2_delta": speed2_delta,
        "max_position_delta": max_position_delta,
        "max_velocity_delta": max_velocity_delta,
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
    ranks = args.mpi_ranks if args.profile == "docs" else 1
    final_time = 1.28 if args.profile == "docs" else 0.08
    run_dir = _run(
        args, "frozen_turbulent", "cr_frozen_turbulent_field.athinput",
        [
            "job/basename=cr_acc_turbulent_docs",
            "time/nlim=2000",
            f"time/tlim={final_time}",
        ],
        ranks=ranks)
    expected = [128]*5
    summary = cr_tracer_inspect.summarize_restart(run_dir / "prst")
    cr_tracer_inspect.validate_expected_counts(summary, sum(expected), expected)
    moments = cr_tracer_inspect.validate_moments(run_dir, expected)
    spectra = cr_tracer_inspect.validate_joint_spectra(run_dir, expected, 16, 16)
    rows = [{"species": row["species"], "mean_mu": row["mean_mu"],
             "mean_mu2": row["mean_mu2"],
             "mean_speed2": row["mean_speed2"],
             "rms_displacement": math.sqrt(
                 row["mean_dx2"] + row["mean_dy2"] + row["mean_dz2"])}
            for row in moments]

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.7),
                             constrained_layout=True)
    species = [row["species"] for row in rows]
    axes[0].bar([s - 0.18 for s in species], [row["mean_mu"] for row in rows],
                width=0.36, label="<mu>")
    axes[0].bar([s + 0.18 for s in species], [row["mean_mu2"] for row in rows],
                width=0.36, label="<mu^2>")
    axes[0].set_xlabel("species")
    axes[0].set_ylabel("pitch-angle moment")
    axes[0].set_title("Pitch-angle diagnostics")
    axes[0].legend()
    axes[1].plot([0.25, 0.5, 1.0, 2.0, 4.0],
                 [row["rms_displacement"] for row in rows], marker="o")
    axes[1].set_xscale("log", base=2)
    axes[1].set_xlabel("mass / gyroradius parameter")
    axes[1].set_ylabel("RMS displacement")
    axes[1].set_title("Transport versus gyroradius")
    axes[1].grid(True, alpha=0.3)
    fig_path = args.output_dir / "frozen_turbulent_moments.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "frozen_turbulent",
        "moments": rows,
        "masses": [0.25, 0.5, 1.0, 2.0, 4.0],
        "joint_spectrum_totals": [
            sum(sum(values) for values in species_rows)
            for species_rows in spectra["histogram"]
        ],
        "ranks": ranks,
        "final_time": final_time,
        "figure": _artifact_path(fig_path),
    }


def study_pitch_angle_decorrelation(args: argparse.Namespace) -> Dict:
    ranks = args.mpi_ranks if args.profile == "docs" else min(2, args.mpi_ranks)
    resolutions = [32, 64, 128] if args.profile == "docs" else [16]
    final_time = 1.28 if args.profile == "docs" else 0.64
    structured_curves = []
    structured_final = []
    structured_final_moments = []
    structured_speed2_errors = []
    for nx in resolutions:
        ppc = 1024.0/(nx**3)
        meshblock_nx = 8 if nx == 16 else 16
        run_dir = _run(
            args, f"pitch_angle_structured_{nx}",
            "cr_pitch_angle_decorrelation.athinput",
            [
                f"job/basename=cr_acc_pitch_structured_{nx}",
                f"mesh/nx1={nx}",
                f"mesh/nx2={nx}",
                f"mesh/nx3={nx}",
                f"meshblock/nx1={meshblock_nx}",
                f"meshblock/nx2={meshblock_nx}",
                f"meshblock/nx3={meshblock_nx}",
                f"particles/ppc={ppc:.17e}",
                "time/nlim=200000",
                f"time/tlim={final_time}",
                "output1/dt=0.08",
                f"output4/dt={final_time}",
            ],
            ranks=ranks)
        sequence = _restart_sequence(run_dir)
        initial = sequence[0]["particles"]
        curve = [
            _pitch_angle_correlation(initial, sample["particles"])
            for sample in sequence
        ]
        structured_curves.append((sequence, curve))
        structured_final.append(curve[-1])
        structured_final_moments.append(_ensemble_mu_moments(sequence[-1]["particles"]))
        structured_speed2_errors.append(_ensemble_max_speed2_error(sequence))

    control_dir = _run(
        args, "pitch_angle_uniform_control",
        "cr_pitch_angle_uniform_control.athinput",
        [
            "job/basename=cr_acc_pitch_uniform_control",
            "mesh/nx1=32",
            "mesh/nx2=32",
            "mesh/nx3=32",
            "meshblock/nx1=16",
            "meshblock/nx2=16",
            "meshblock/nx3=16",
            f"particles/ppc={1024.0/(32**3):.17e}",
            "time/nlim=200000",
            f"time/tlim={final_time}",
            "output1/dt=0.08",
            f"output4/dt={final_time}",
        ],
        ranks=ranks)
    control_sequence = _restart_sequence(control_dir)
    control_initial = control_sequence[0]["particles"]
    control_curve = [
        _pitch_angle_correlation(control_initial, sample["particles"])
        for sample in control_sequence
    ]
    control_final_moments = _ensemble_mu_moments(control_sequence[-1]["particles"])
    control_speed2_error = _ensemble_max_speed2_error(control_sequence)
    finest_sequence = structured_curves[-1][0]
    initial_moments = _ensemble_mu_moments(finest_sequence[0]["particles"])
    finest_moment_curve = [
        _ensemble_mu_moments(sample["particles"]) for sample in finest_sequence
    ]
    bins = np.linspace(-1.0, 1.0, 17)
    centers = 0.5*(bins[:-1] + bins[1:])
    initial_histogram, _ = np.histogram(
        [_pitch_angle(particle) for particle in finest_sequence[0]["particles"].values()],
        bins=bins, density=True)
    final_particles = finest_sequence[-1]["particles"].values()
    final_histogram, _ = np.histogram(
        [_pitch_angle(particle) for particle in final_particles],
        bins=bins, density=True)
    control_particles = control_sequence[-1]["particles"].values()
    control_histogram, _ = np.histogram(
        [_pitch_angle(particle) for particle in control_particles],
        bins=bins, density=True)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.8),
                             constrained_layout=True)
    for nx, (sequence, curve) in zip(resolutions, structured_curves):
        axes[0].plot([sample["time"] for sample in sequence], curve,
                     label=f"structured {nx}^3")
    axes[0].plot([sample["time"] for sample in control_sequence], control_curve,
                 linestyle="--", color="black", label="uniform control")
    axes[0].set_xlabel("time")
    axes[0].set_ylabel("C_mu(t)")
    axes[0].set_title("Correlation")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=8)
    times = [sample["time"] for sample in finest_sequence]
    axes[1].plot(times, [row["mean_mu"] for row in finest_moment_curve],
                 label="<mu>")
    axes[1].plot(times, [row["mean_mu2"] for row in finest_moment_curve],
                 label="<mu^2>")
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("moment")
    axes[1].set_title("Structured field at 128^3")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=8)
    axes[2].step(centers, initial_histogram, where="mid", label="initial")
    axes[2].step(centers, final_histogram, where="mid",
                 label="structured final")
    axes[2].step(centers, control_histogram, where="mid",
                 label="uniform final", linestyle="--")
    axes[2].set_xlabel("mu")
    axes[2].set_ylabel("f(mu)")
    axes[2].set_title("Distribution evolution")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(fontsize=8)
    fig_path = args.output_dir / "pitch_angle_decorrelation.png"
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)
    return {
        "case": "pitch_angle_decorrelation",
        "nx": resolutions,
        "ranks": ranks,
        "final_time": final_time,
        "structured_final_correlation": structured_final,
        "initial_moments": initial_moments,
        "structured_final_moments": structured_final_moments,
        "structured_max_speed2_error": structured_speed2_errors,
        "uniform_control_final_correlation": control_curve[-1],
        "uniform_control_final_moments": control_final_moments,
        "uniform_control_max_speed2_error": control_speed2_error,
        "histogram_bin_centers": centers.tolist(),
        "initial_f_mu": initial_histogram.tolist(),
        "finest_structured_final_f_mu": final_histogram.tolist(),
        "uniform_control_final_f_mu": control_histogram.tolist(),
        "requested_reference_nx": 256,
        "reference_status": "256^3 structured reference pending",
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
    "pitch_angle_decorrelation": study_pitch_angle_decorrelation,
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
    parser.add_argument("--profile", choices=("ci", "docs"), default="ci",
                        help="Workload scale for regression or documentation")
    parser.add_argument("--mpi-ranks", type=int, default=8,
                        help="MPI ranks for distributed documentation cases")
    parser.add_argument("--cases", default=",".join(STUDIES.keys()),
                        help="Comma-separated studies to run")
    args = parser.parse_args()
    if args.mpi_ranks < 1 or args.mpi_ranks > 8:
        raise ValueError("--mpi-ranks must be in the inclusive range 1 to 8")
    args.executions = []
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.run_root.mkdir(parents=True, exist_ok=True)
    git_revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    git_dirty = bool(subprocess.check_output(
        ["git", "status", "--porcelain"], cwd=ROOT, text=True).strip())
    mpi_version = subprocess.check_output(
        [args.mpiexec, "--version"], text=True).splitlines()[0]

    selected = [case.strip() for case in args.cases.split(",") if case.strip()]
    results = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "athena": str(args.athena),
        "git_revision": git_revision,
        "git_worktree_dirty": git_dirty,
        "mpi_version": mpi_version,
        "profile": args.profile,
        "mpi_ranks": args.mpi_ranks,
        "cases": [],
    }
    for case in selected:
        if case not in STUDIES:
            raise ValueError(f"Unknown accuracy study '{case}'")
        results["cases"].append(STUDIES[case](args))
    results["qualitative_figures"] = study_qualitative_figures(args)
    results["executions"] = args.executions

    json_path = args.output_dir / "cr_tracer_accuracy_summary.json"
    json_path.write_text(json.dumps(results, indent=2) + "\n")
    md_path = args.output_dir / "cr_tracer_accuracy_results.md"
    lines = [
        "---", "orphan: true", "---", "", "# CR Tracer Accuracy Results", "",
        f"- profile: `{args.profile}`",
        f"- git revision: `{git_revision}`",
        f"- git worktree dirty during generation: `{git_dirty}`",
        f"- MPI: `{mpi_version}`",
        f"- distributed documentation ranks: `{args.mpi_ranks}`",
        "",
    ]
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
            f"`{result['track_count']}` tracks, `{result['ranks']}` ranks)")
    lines.extend(["", "## executions", "",
                  "| case | ranks | wall seconds |",
                  "|------|------:|-------------:|"])
    for execution in args.executions:
        lines.append(
            f"| {execution['case']} | {execution['ranks']} | "
            f"{execution['wall_seconds']:.3f} |")
    md_path.write_text("\n".join(lines).rstrip() + "\n")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
