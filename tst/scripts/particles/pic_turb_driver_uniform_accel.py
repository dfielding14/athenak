"""Production-path one-cycle invariant for staged turbulence forcing."""

from __future__ import annotations

import glob
import logging
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import numpy as np


_SOURCE_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_SOURCE_ROOT / "vis" / "python"))
import bin_convert_new as bin_convert  # noqa: E402


logger = logging.getLogger("athena" + __name__[7:])

_MHD_INPUT = (
    _SOURCE_ROOT / "inputs" / "tests" /
    "turb_driver_uniform_accel_mhd.athinput"
)
_PIC_INPUT = (
    _SOURCE_ROOT / "inputs" / "tests" /
    "turb_driver_uniform_accel_pic_vl2.athinput"
)
_MPIEXEC = shlex.split(os.environ.get("MPIEXEC", "mpiexec"))

_RHO0 = 1.25
_PRESSURE0 = 1.0
_GAMMA = 5.0/3.0
_VELOCITY0 = np.array([0.25, 0.5, 0.0])
_ACCELERATION = np.array([4.0, -2.0, 1.0])

_RESULTS = {}
_SKIPPED = False


def _exe_dir() -> Path:
    override = os.environ.get("ATHENA_TURB_UNIFORM_ACCEL_EXE_DIR")
    return Path(override) if override else Path.cwd() / "build" / "src"


def _configuration() -> str:
    proc = subprocess.run(
        ["./athena", "-c"], cwd=_exe_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError("Unable to query Athena configuration")
    return (proc.stdout or "") + (proc.stderr or "")


def _configured_problem(config: str) -> str:
    for line in config.splitlines():
        if "Problem generator:" in line:
            return line.split(":", 1)[1].strip()
    return ""


def _mpi_launcher_available() -> bool:
    return bool(_MPIEXEC and shutil.which(_MPIEXEC[0]))


def _remove_outputs(basename: str) -> None:
    patterns = (
        _exe_dir() / "bin" / f"{basename}.*.bin",
        _exe_dir() / f"{basename}.mhd.hst",
    )
    for pattern in patterns:
        for path in glob.glob(str(pattern)):
            os.remove(path)


def _latest(pattern: Path) -> Path:
    matches = sorted(glob.glob(str(pattern)))
    if not matches:
        raise RuntimeError(f"Missing output matching {pattern}")
    return Path(matches[-1])


def _history_rows(basename: str) -> np.ndarray:
    path = _exe_dir() / f"{basename}.mhd.hst"
    rows = np.loadtxt(path)
    if rows.ndim == 1:
        rows = rows.reshape(1, -1)
    unique_rows = []
    for row in rows:
        if not unique_rows or row[0] != unique_rows[-1][0]:
            unique_rows.append(row)
        else:
            unique_rows[-1] = row
    if len(unique_rows) < 2:
        raise RuntimeError(f"Expected initial and final history rows in {path}")
    return np.stack(unique_rows)


def _expected_state(dt: float) -> tuple[np.ndarray, np.ndarray]:
    momentum0 = _RHO0*_VELOCITY0
    energy0 = (
        _PRESSURE0/(_GAMMA - 1.0) +
        0.5*_RHO0*np.dot(_VELOCITY0, _VELOCITY0)
    )
    velocity1 = _VELOCITY0 + _ACCELERATION*dt
    momentum1 = _RHO0*velocity1
    energy1 = (
        _PRESSURE0/(_GAMMA - 1.0) +
        0.5*_RHO0*np.dot(velocity1, velocity1)
    )
    return (
        np.concatenate(([_RHO0], momentum0, [energy0])),
        np.concatenate(([_RHO0], momentum1, [energy1])),
    )


def _measure(basename: str) -> dict:
    history = _history_rows(basename)
    initial = history[0]
    final = history[-1]
    dt = float(final[0] - initial[0])
    if dt <= 0.0:
        raise RuntimeError("Uniform-acceleration case did not advance in time")
    expected_initial, expected_final = _expected_state(dt)

    state_data = bin_convert.read_binary_as_athdf(_latest(
        _exe_dir() / "bin" / f"{basename}.mhd_u.*.bin"
    ))
    state = np.stack([
        np.asarray(state_data[name], dtype=float)
        for name in ("dens", "mom1", "mom2", "mom3", "ener")
    ])
    state_error = float(np.max(np.abs(
        state - expected_final[:, np.newaxis, np.newaxis, np.newaxis]
    )))

    force_data = bin_convert.read_binary_as_athdf(_latest(
        _exe_dir() / "bin" / f"{basename}.turb_force.*.bin"
    ))
    force = np.stack([
        np.asarray(force_data[name], dtype=float)
        for name in ("force1", "force2", "force3")
    ])
    force_error = float(np.max(np.abs(
        force - _ACCELERATION[:, np.newaxis, np.newaxis, np.newaxis]
    )))

    measured_power = float((final[6] - initial[6])/dt)
    expected_power = float(
        _RHO0*np.dot(_VELOCITY0, _ACCELERATION) +
        0.5*_RHO0*np.dot(_ACCELERATION, _ACCELERATION)*dt
    )
    return {
        "dt": dt,
        "initial_state": initial[2:7],
        "final_state": final[2:7],
        "expected_initial": expected_initial,
        "expected_final": expected_final,
        "history_error": float(max(
            np.max(np.abs(initial[2:7] - expected_initial)),
            np.max(np.abs(final[2:7] - expected_final)),
        )),
        "state_error": state_error,
        "force_error": force_error,
        "measured_power": measured_power,
        "expected_power": expected_power,
        "power_error": abs(measured_power - expected_power),
    }


def _run_case(label: str, input_path: Path, nproc: int, overrides=()) -> None:
    basename = "TurbUniformAccel_" + label
    _remove_outputs(basename)
    command = [
        "./athena", "-i", str(input_path), "job/basename=" + basename,
        *overrides,
    ]
    if nproc > 1:
        command = _MPIEXEC + ["-n", str(nproc)] + command
    proc = subprocess.run(
        command, cwd=_exe_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        output = (proc.stdout or "") + (proc.stderr or "")
        raise RuntimeError(f"Uniform-acceleration case {label} failed:\n{output}")
    _RESULTS[label] = _measure(basename)


def run(**kwargs):
    del kwargs
    global _SKIPPED
    logger.debug("Running test " + __name__)
    _SKIPPED = False
    _RESULTS.clear()
    config = _configuration()
    if _configured_problem(config) != "turb_uniform_accel_test":
        logger.info(
            "Skipping: configure with -DPROBLEM=turb_uniform_accel_test"
        )
        _SKIPPED = True
        return

    _run_case("mhd_rk1_np1", _MHD_INPUT, 1, ("time/integrator=rk1",))
    _run_case("mhd_rk2_np1", _MHD_INPUT, 1)
    _run_case("pic_vl2_np1", _PIC_INPUT, 1)

    if "MPI parallelism:            ON" in config and _mpi_launcher_available():
        _run_case("pic_vl2_np2", _PIC_INPUT, 2)


def analyze():
    if _SKIPPED:
        return True

    ok = True
    for label, result in _RESULTS.items():
        logger.info(
            "%s dt=% .8e history_err=% .8e field_err=% .8e force_err=% .8e",
            label, result["dt"], result["history_error"],
            result["state_error"], result["force_error"],
        )
        logger.info(
            "%s forcing_power=% .16e expected=% .16e residual=% .8e",
            label, result["measured_power"], result["expected_power"],
            result["power_error"],
        )
        ok = (result["history_error"] < 2.0e-12) and ok
        ok = (result["state_error"] < 2.0e-6) and ok
        ok = (result["force_error"] < 2.0e-6) and ok
        ok = (result["power_error"] < 2.0e-10) and ok

    heun = _RESULTS["mhd_rk2_np1"]
    midpoint = _RESULTS["pic_vl2_np1"]
    selector_error = float(np.max(np.abs(
        heun["final_state"] - midpoint["final_state"]
    )))
    selector_power_error = abs(
        heun["measured_power"] - midpoint["measured_power"]
    )
    logger.info(
        "Heun/PIC-midpoint state error=% .8e power error=% .8e",
        selector_error, selector_power_error,
    )
    ok = (selector_error < 2.0e-12) and ok
    ok = (selector_power_error < 2.0e-10) and ok

    if "pic_vl2_np2" in _RESULTS:
        decomp = _RESULTS["pic_vl2_np2"]
        decomp_error = float(np.max(np.abs(
            midpoint["final_state"] - decomp["final_state"]
        )))
        logger.info("PIC midpoint np1/np2 state error=% .8e", decomp_error)
        ok = (decomp_error < 2.0e-12) and ok

    return ok
