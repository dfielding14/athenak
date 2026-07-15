"""Relativistic gyro timestep and uniform-field orbit regression."""

from __future__ import annotations

import glob
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena

_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
sys.path.insert(0, os.path.join(_SOURCE_ROOT, "tst", "publication"))
from pvtk_particles import read_particle_vtk  # noqa: E402

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_relativistic_gyro_timestep.athinput"
_MPIEXEC = shlex.split(os.environ.get("MPIEXEC", "mpiexec"))
_THETA_VALUES = (0.8, 0.4, 0.2)
_CFL = 0.5
_LIGHT_SPEED = 1.0
_GAMMA = 100.0
_Q_OVER_M = 100.0
_BMAG = 1.0
_OMEGA = _Q_OVER_M * _BMAG / _GAMMA
_SPEED = _LIGHT_SPEED * np.sqrt(1.0 - 1.0 / (_GAMMA * _GAMMA))
_RADIUS = _SPEED / _OMEGA
_TEND = 2.0 * np.pi / _OMEGA
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), "build", "src")


def _athena_input_path():
    return "../../" + athena.athena_rel_path + "inputs/" + _INPUT_DECK


def _athena_mpi_enabled():
    proc = subprocess.run(
        ["./athena", "-c"],
        cwd=_athena_exe_dir(),
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError("Unable to query Athena configuration with -c")
    output = (proc.stdout or "") + (proc.stderr or "")
    return "MPI parallelism:            ON" in output


def _case_tag(theta_max):
    return str(theta_max).replace(".", "p")


def _remove_outputs(basename):
    patterns = [
        os.path.join(_athena_exe_dir(), "pvtk", basename + ".*.part.vtk"),
        os.path.join(_athena_exe_dir(), basename + "-errs.dat"),
    ]
    for pattern in patterns:
        for path in glob.glob(pattern):
            os.remove(path)


def _extract_initial_dt(output):
    pattern = re.compile(
        r"cycle=\s*0\s+time=\s*[0-9eE+\-.]+\s+dt=\s*([0-9eE+\-.]+)"
    )
    match = pattern.search(output)
    if match is None:
        raise RuntimeError("Could not find the cycle-0 timestep in Athena output")
    return float(match.group(1))


def _snapshot_header(path):
    with open(path, "rb") as stream:
        header = stream.read(512)
    match = re.search(
        rb"# AthenaK particle data at time=\s*([^ ]+)\s+nranks=.*cycle=([0-9]+)",
        header,
    )
    if match is None:
        raise RuntimeError("Could not read particle VTK header in " + path)
    return float(match.group(1)), int(match.group(2))


def _load_snapshots(basename):
    pattern = os.path.join(
        _athena_exe_dir(), "pvtk", basename + ".prtcl_all.*.part.vtk"
    )
    paths = glob.glob(pattern)
    if not paths:
        raise RuntimeError("No particle VTK files found for " + basename)

    snapshots = []
    for path in paths:
        time, cycle = _snapshot_header(path)
        data = read_particle_vtk(path)
        if "ptag" not in data.scalars or "vel" not in data.vectors:
            raise RuntimeError("Particle VTK output lacks ptag or vel in " + path)
        order = np.argsort(data.scalars["ptag"])
        snapshots.append(
            {
                "time": time,
                "cycle": cycle,
                "tags": data.scalars["ptag"][order],
                "points": data.points[order],
                "velocity": data.vectors["vel"][order],
            }
        )
    snapshots.sort(key=lambda item: (item["cycle"], item["time"]))
    return snapshots


def _measure_orbit(basename, output, theta_max):
    snapshots = _load_snapshots(basename)
    times = np.array([item["time"] for item in snapshots])
    cycles = np.array([item["cycle"] for item in snapshots], dtype=int)
    points = np.stack([item["points"] for item in snapshots])
    velocity = np.stack([item["velocity"] for item in snapshots])
    tags = snapshots[0]["tags"]

    if any(not np.array_equal(item["tags"], tags) for item in snapshots):
        raise RuntimeError("Particle tags changed during " + basename)
    if cycles[0] != 0 or abs(times[0]) > 1.0e-14:
        raise RuntimeError("Orbit history does not start at cycle 0, t=0")

    phase = _OMEGA * times
    expected_velocity = np.zeros_like(velocity)
    expected_velocity[:, :, 0] = _SPEED * np.cos(phase)[:, None]
    expected_velocity[:, :, 1] = -_SPEED * np.sin(phase)[:, None]

    expected_displacement = np.zeros_like(points)
    expected_displacement[:, :, 0] = _RADIUS * np.sin(phase)[:, None]
    expected_displacement[:, :, 1] = (
        _RADIUS * (np.cos(phase) - 1.0)[:, None]
    )
    displacement = points - points[0][None, :, :]

    speed = np.linalg.norm(velocity, axis=2)
    gamma = 1.0 / np.sqrt(1.0 - np.minimum(speed * speed, 1.0 - 1.0e-15))
    kinetic_energy = (gamma - 1.0) * _LIGHT_SPEED * _LIGHT_SPEED

    measured_phase = np.arctan2(
        -np.mean(velocity[:, :, 1], axis=1),
        np.mean(velocity[:, :, 0], axis=1),
    )
    phase_delta = measured_phase - phase
    phase_error = np.abs(np.arctan2(np.sin(phase_delta), np.cos(phase_delta)))

    radius_vector = np.empty_like(points[:, :, :2])
    radius_vector[:, :, 0] = -_RADIUS * velocity[:, :, 1] / _SPEED
    radius_vector[:, :, 1] = _RADIUS * velocity[:, :, 0] / _SPEED
    inferred_centers = points[:, :, :2] - radius_vector
    initial_centers = inferred_centers[0]
    radius = np.linalg.norm(points[:, :, :2] - initial_centers[None, :, :], axis=2)

    return {
        "theta_max": theta_max,
        "initial_dt": _extract_initial_dt(output),
        "times": times,
        "cycles": cycles,
        "tags": tags,
        "points": points,
        "velocity": velocity,
        "nparticle": points.shape[1],
        "trajectory_error": float(
            np.max(np.linalg.norm(displacement - expected_displacement, axis=2))
        ),
        "velocity_error": float(
            np.max(np.linalg.norm(velocity - expected_velocity, axis=2))
        ),
        "final_phase_error": float(phase_error[-1]),
        "closure_error": float(
            np.max(np.linalg.norm(points[-1] - points[0], axis=1))
        ),
        "speed_rel_error": float(np.max(np.abs(speed - _SPEED)) / _SPEED),
        "gamma_rel_error": float(np.max(np.abs(gamma - _GAMMA)) / _GAMMA),
        "energy_rel_error": float(
            np.max(np.abs(kinetic_energy - (_GAMMA - 1.0))) / (_GAMMA - 1.0)
        ),
        "center_error": float(
            np.max(np.linalg.norm(inferred_centers - initial_centers[None, :, :], axis=2))
        ),
        "radius_error": float(np.max(np.abs(radius - _RADIUS))),
        "parallel_position_error": float(
            np.max(np.abs(points[:, :, 2] - points[0, :, 2][None, :]))
        ),
        "parallel_velocity_error": float(np.max(np.abs(velocity[:, :, 2]))),
        "translated_position_spread": float(
            np.max(np.abs(displacement - displacement[:, :1, :]))
        ),
        "velocity_spread": float(
            np.max(np.abs(velocity - velocity[:, :1, :]))
        ),
    }


def _run_case(label, basename, theta_max, nproc):
    _remove_outputs(basename)
    command = [
        "./athena",
        "-i",
        _athena_input_path(),
        "job/basename=" + basename,
        "particles/pic_theta_max=" + str(theta_max),
    ]
    if nproc > 1:
        command = [*_MPIEXEC, "-n", str(nproc), *command]
    logger.info("Executing %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command,
        cwd=_athena_exe_dir(),
        capture_output=True,
        text=True,
        check=False,
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Command failed for " + label + "\n" + output)
    return _measure_orbit(basename, output, theta_max)


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    _RESULTS["serial"] = []
    for theta_max in _THETA_VALUES:
        basename = "pic_rel_gyro_timestep_theta" + _case_tag(theta_max)
        _RESULTS["serial"].append(
            _run_case("serial theta=" + str(theta_max), basename, theta_max, 1)
        )

    launcher_available = bool(_MPIEXEC) and shutil.which(_MPIEXEC[0]) is not None
    if _athena_mpi_enabled() and launcher_available:
        theta_max = _THETA_VALUES[-1]
        _RESULTS["mpi2"] = _run_case(
            "mpi2 theta=" + str(theta_max),
            "pic_rel_gyro_timestep_mpi2",
            theta_max,
            2,
        )
    else:
        logger.info("Skipping MPI2 orbit parity: MPI build or launcher unavailable")


def _convergence_orders(cases, metric):
    return [
        np.log(cases[index][metric] / cases[index + 1][metric])
        / np.log(cases[index]["initial_dt"] / cases[index + 1]["initial_dt"])
        for index in range(len(cases) - 1)
    ]


def analyze():
    logger.debug("Analyzing test " + __name__)
    ok = True
    serial = _RESULTS["serial"]

    for case in serial:
        expected_dt = _CFL * case["theta_max"] * _GAMMA / (_Q_OVER_M * _BMAG)
        logger.info(
            "theta=%g dt=% .8e trajectory=% .8e phase=% .8e closure=% .8e "
            "speed=% .8e gamma=% .8e center=% .8e radius=% .8e",
            case["theta_max"],
            case["initial_dt"],
            case["trajectory_error"],
            case["final_phase_error"],
            case["closure_error"],
            case["speed_rel_error"],
            case["gamma_rel_error"],
            case["center_error"],
            case["radius_error"],
        )
        ok = abs(case["initial_dt"] - expected_dt) <= 1.0e-8 and ok
        ok = abs(case["times"][-1] - _TEND) <= 1.0e-12 and ok
        ok = case["nparticle"] == 8 and ok
        ok = np.unique(case["tags"]).size == 8 and ok
        ok = np.unique(case["points"][0], axis=0).shape[0] == 8 and ok
        ok = bool(np.all(np.isfinite(case["points"]))) and ok
        ok = bool(np.all(np.isfinite(case["velocity"]))) and ok
        ok = case["speed_rel_error"] <= 2.0e-7 and ok
        ok = case["gamma_rel_error"] <= 2.0e-3 and ok
        ok = case["energy_rel_error"] <= 2.0e-3 and ok
        ok = case["center_error"] <= 1.5e-5 and ok
        ok = case["radius_error"] <= 1.5e-5 and ok
        ok = case["parallel_position_error"] <= 1.0e-6 and ok
        ok = case["parallel_velocity_error"] <= 1.0e-7 and ok
        ok = case["translated_position_spread"] <= 1.0e-5 and ok
        ok = case["velocity_spread"] <= 1.0e-7 and ok

    for metric in ["trajectory_error", "final_phase_error", "closure_error"]:
        values = [case[metric] for case in serial]
        orders = _convergence_orders(serial, metric)
        logger.info("%s values=%s orders=%s", metric, values, orders)
        ok = bool(np.all(np.diff(values) < 0.0)) and ok
        ok = min(orders) >= 1.8 and ok

    if "mpi2" in _RESULTS:
        reference = serial[-1]
        mpi2 = _RESULTS["mpi2"]
        same_history = np.array_equal(reference["cycles"], mpi2["cycles"])
        same_history = same_history and np.allclose(
            reference["times"], mpi2["times"], atol=1.0e-14, rtol=0.0
        )
        point_error = float(np.max(np.abs(reference["points"] - mpi2["points"])))
        velocity_error = float(
            np.max(np.abs(reference["velocity"] - mpi2["velocity"]))
        )
        logger.info(
            "serial_vs_mpi2 same_history=%s point_error=% .8e "
            "velocity_error=% .8e",
            same_history,
            point_error,
            velocity_error,
        )
        ok = same_history and ok
        ok = np.array_equal(reference["tags"], mpi2["tags"]) and ok
        ok = abs(reference["initial_dt"] - mpi2["initial_dt"]) <= 1.0e-12 and ok
        ok = point_error <= 2.0e-6 and ok
        ok = velocity_error <= 2.0e-7 and ok

    return ok
