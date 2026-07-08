"""Initial-field and six-beam loading regression for the PIC dynamo pilot."""

from __future__ import annotations

import glob
import logging
import os
import shlex
import shutil
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena  # noqa: E402

_PUBLICATION_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "publication")
)
_VIS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "vis", "python")
)
sys.path.insert(0, _PUBLICATION_DIR)
sys.path.insert(0, _VIS_DIR)

import bin_convert_new as bin_convert  # noqa: E402
from pvtk_particles import read_particle_vtk  # noqa: E402

logger = logging.getLogger("athena" + __name__[7:])

_INPUT = "tests/pic_turbulent_dynamo_smoke.athinput"
_MHD_INPUT = "tests/turb_driver_mhd_smoke.athinput"
_MPIEXEC = shlex.split(os.environ.get("MPIEXEC", "mpiexec"))
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_RESULTS = {}
_MHD_RESULTS = {}
_MHD_RESTART = {}
_SKIPPED = False


def _exe_dir():
    return os.path.join(os.getcwd(), "build", "src")


def _configuration():
    proc = subprocess.run(
        ["./athena", "-c"], cwd=_exe_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError("Unable to query Athena configuration")
    return (proc.stdout or "") + (proc.stderr or "")


def _remove_outputs(basename):
    patterns = (
        os.path.join(_exe_dir(), "bin", basename + ".*.bin"),
        os.path.join(_exe_dir(), "pvtk", basename + ".*.vtk"),
        os.path.join(_exe_dir(), basename + ".*.hst"),
        os.path.join(_exe_dir(), "rst", basename + ".*"),
    )
    for pattern in patterns:
        for path in glob.glob(pattern):
            os.remove(path)


def _mpi_launcher_available():
    return bool(_MPIEXEC and shutil.which(_MPIEXEC[0]))


def _run(basename, nproc, overrides=()):
    input_path = os.path.join(_SOURCE_ROOT, "inputs", _INPUT)
    command = [
        "./athena", "-i", input_path, "job/basename=" + basename,
        "time/nlim=2", "time/tlim=1.0",
        *overrides,
    ]
    if nproc > 1:
        command = _MPIEXEC + ["-n", str(nproc)] + command
    proc = subprocess.run(
        command, cwd=_exe_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        output = (proc.stdout or "") + (proc.stderr or "")
        raise RuntimeError("PIC turbulent-dynamo smoke failed:\n" + output)


def _run_mhd(basename, nproc):
    input_path = os.path.join(_SOURCE_ROOT, "inputs", _MHD_INPUT)
    command = ["./athena", "-i", input_path, "job/basename=" + basename]
    if nproc > 1:
        command = _MPIEXEC + ["-n", str(nproc)] + command
    proc = subprocess.run(
        command, cwd=_exe_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        output = (proc.stdout or "") + (proc.stderr or "")
        raise RuntimeError("MHD turbulence-driver smoke failed:\n" + output)


def _expect_mhd_failure(parameter, message):
    input_path = os.path.join(_SOURCE_ROOT, "inputs", _MHD_INPUT)
    command = [
        "./athena", "-i", input_path,
        "job/basename=TurbDriverInvalid", parameter,
    ]
    proc = subprocess.run(
        command, cwd=_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode == 0 or message not in output:
        raise RuntimeError(
            "Invalid turbulence input did not fail as expected:\n" + output
        )


def _check_legacy_bounds_upgrade():
    input_path = os.path.join(_SOURCE_ROOT, "inputs", _MHD_INPUT)
    _remove_outputs("TurbDriverLegacyBounds")
    command = [
        "./athena", "-i", input_path,
        "job/basename=TurbDriverLegacyBounds",
        "time/nlim=0", "time/tlim=0.0",
        "turb_driving/min_kx=0", "turb_driving/min_ky=0",
        "turb_driving/min_kz=0",
    ]
    proc = subprocess.run(
        command, cwd=_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0 or "upgrading legacy nonnegative" not in output:
        raise RuntimeError("Legacy turbulence bounds were not upgraded:\n" + output)


def _latest(directory, pattern):
    paths = sorted(glob.glob(os.path.join(_exe_dir(), directory, pattern)))
    if not paths:
        raise RuntimeError("Missing output for " + pattern)
    return paths[-1]


def _earliest(directory, pattern):
    paths = sorted(glob.glob(os.path.join(_exe_dir(), directory, pattern)))
    if not paths:
        raise RuntimeError("Missing output for " + pattern)
    return paths[0]


def _first_history_row(basename):
    path = os.path.join(_exe_dir(), basename + ".user.hst")
    with open(path, "r", encoding="utf-8") as stream:
        rows = [line for line in stream if line.strip() and not line.startswith("#")]
    if not rows:
        raise RuntimeError("Missing user-history rows for " + basename)
    return np.fromstring(rows[0], sep=" ")


def _history_rows(basename, stream_name):
    path = os.path.join(_exe_dir(), basename + "." + stream_name + ".hst")
    with open(path, "r", encoding="utf-8") as stream:
        rows = [
            np.fromstring(line, sep=" ")
            for line in stream
            if line.strip() and not line.startswith("#")
        ]
    if not rows:
        raise RuntimeError("Missing " + stream_name + " history for " + basename)
    unique_rows = []
    for row in rows:
        if not unique_rows or row[0] != unique_rows[-1][0]:
            unique_rows.append(row)
    return np.stack(unique_rows)


def _force_metrics(force):
    shape = force.shape[1:]
    wave_axes = [np.fft.fftfreq(n)*n for n in shape]
    kz, ky, kx = np.meshgrid(*wave_axes, indexing="ij")
    force_hat = np.fft.fftn(force, axes=(1, 2, 3))
    wave_dot_force = (
        kx*force_hat[0] + ky*force_hat[1] + kz*force_hat[2]
    )
    wave_squared = kx*kx + ky*ky + kz*kz
    denominator = np.sqrt(np.sum(
        wave_squared*np.sum(np.abs(force_hat)**2, axis=0)
    ))
    divergence_ratio = float(
        np.sqrt(np.sum(np.abs(wave_dot_force)**2))/denominator
    )
    power = np.sum(np.abs(force_hat)**2, axis=0)
    shell = (wave_squared >= 1.0) & (wave_squared <= 9.0)
    outside_shell_fraction = float(np.sum(power[~shell])/np.sum(power))

    indices = [
        {int(round(wave)): n for n, wave in enumerate(axis)}
        for axis in wave_axes
    ]
    if shape[0] == 1:
        mixed_sign_modes = [(1, -1, 0), (2, -1, 0), (1, -2, 0)]
    else:
        mixed_sign_modes = [(1, -1, 0), (1, -1, 1), (2, -1, 0)]
    mixed_sign_power = []
    peak_power = float(np.max(power))
    for nx, ny, nz in mixed_sign_modes:
        mixed_sign_power.append(
            float(power[indices[0][nz], indices[1][ny], indices[2][nx]]/
                  peak_power)
        )
    return {
        "divergence_ratio": divergence_ratio,
        "outside_shell_fraction": outside_shell_fraction,
        "mixed_sign_relative_power": np.asarray(mixed_sign_power),
        "out_of_plane_fraction": float(
            np.linalg.norm(force[2])/np.linalg.norm(force)
        ),
    }


def _measure(basename):
    field_data = bin_convert.read_binary_as_athdf(
        _earliest("bin", basename + ".mhd_w_bcc.*.bin")
    )
    divb_data = bin_convert.read_binary_as_athdf(
        _earliest("bin", basename + ".mhd_divb.*.bin")
    )
    magnetic = np.stack([
        np.asarray(field_data[name], dtype=float)
        for name in ("bcc1", "bcc2", "bcc3")
    ])
    particles = read_particle_vtk(
        _earliest("pvtk", basename + ".prtcl_all.*.part.vtk")
    )
    species = particles.scalars["species"]
    states = particles.vectors["vel"]
    group_positions = particles.points.reshape(-1, 6, 3)
    history = _first_history_row(basename)
    mhd_history = _history_rows(basename, "mhd")
    user_history = _history_rows(basename, "user")
    if mhd_history.shape[0] < 3 or user_history.shape[0] < 3:
        raise RuntimeError("Expected initial, pre-driving, and driven history rows")
    before_energy = mhd_history[-2, 6] + user_history[-2, 10]
    after_energy = mhd_history[-1, 6] + user_history[-1, 10]
    driven_time = mhd_history[-1, 0] - mhd_history[-2, 0]
    expected_energy = 0.06*driven_time
    force_data = bin_convert.read_binary_as_athdf(
        _latest("bin", basename + ".turb_force.*.bin")
    )
    force = np.stack([
        np.asarray(force_data[name], dtype=float)
        for name in ("force1", "force2", "force3")
    ])
    return {
        "magnetic": magnetic,
        "magnetic_rms": float(np.sqrt(np.mean(np.sum(magnetic**2, axis=0)))),
        "stored_magnetic_mean": np.mean(
            magnetic, axis=tuple(range(1, magnetic.ndim))
        ),
        "history_magnetic_mean": history[6:9],
        "history_max_divb": float(history[9]),
        "stored_max_divb": float(np.max(np.abs(divb_data["divb"]))),
        "particle_count": int(states.shape[0]),
        "cell_count": int(np.prod(magnetic.shape[1:])),
        "species_counts": np.bincount(species, minlength=6),
        "position_spread": float(np.max(np.ptp(group_positions, axis=1))),
        "momentum_sum": np.sum(states, axis=0),
        "current_sum": 100.0*np.sum(states, axis=0),
        "history_cr_momentum": history[11:14],
        "forcing_energy_increment": float(after_energy - before_energy),
        "forcing_energy_expected": float(expected_energy),
        "forcing_energy_residual": float(
            after_energy - before_energy - expected_energy
        ),
        "force": force,
        **_force_metrics(force),
    }


def _measure_mhd(basename):
    mhd_history = _history_rows(basename, "mhd")
    if mhd_history.shape[0] < 3:
        raise RuntimeError("Expected initial, pre-driving, and driven MHD rows")
    energy_increment = mhd_history[-1, 6] - mhd_history[-2, 6]
    driven_time = mhd_history[-1, 0] - mhd_history[-2, 0]
    energy_expected = 0.06*driven_time
    force_data = bin_convert.read_binary_as_athdf(
        _latest("bin", basename + ".turb_force.*.bin")
    )
    force = np.stack([
        np.asarray(force_data[name], dtype=float)
        for name in ("force1", "force2", "force3")
    ])
    return {
        "forcing_energy_increment": float(energy_increment),
        "forcing_energy_expected": float(energy_expected),
        "forcing_energy_residual": float(energy_increment - energy_expected),
        "final_history": mhd_history[-1],
        "force": force,
        **_force_metrics(force),
    }


def _continue_mhd_from_cycle_one(source_basename):
    restart_path = os.path.join(
        _exe_dir(), "rst", source_basename + ".00001.rst"
    )
    if not os.path.isfile(restart_path):
        raise RuntimeError("Missing cycle-one restart: " + restart_path)
    continuation_basename = "TurbDriverMHD_restart"
    _remove_outputs(continuation_basename)
    command = [
        "./athena", "-r", restart_path,
        "job/basename=" + continuation_basename,
        "time/nlim=2", "time/tlim=1.0",
    ]
    proc = subprocess.run(
        command, cwd=_exe_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        output = (proc.stdout or "") + (proc.stderr or "")
        raise RuntimeError("MHD turbulence restart failed:\n" + output)
    force_data = bin_convert.read_binary_as_athdf(
        _latest("bin", continuation_basename + ".turb_force.*.bin")
    )
    return {
        "force": np.stack([
            np.asarray(force_data[name], dtype=float)
            for name in ("force1", "force2", "force3")
        ]),
        "final_history": _history_rows(continuation_basename, "mhd")[-1],
    }


def run(**kwargs):
    global _SKIPPED
    logger.debug("Running test " + __name__)
    _SKIPPED = False
    _RESULTS.clear()
    _MHD_RESULTS.clear()
    _MHD_RESTART.clear()
    config = _configuration()
    if "Problem generator:          turb" not in config:
        logger.info("Skipping: configure this test with -DPROBLEM=turb")
        _SKIPPED = True
        return

    mpi_available = (
        "MPI parallelism:            ON" in config
        and _mpi_launcher_available()
    )
    cases = [("np1", 1, ())]
    if mpi_available:
        cases.append(("np2", 2, ()))
    elif "MPI parallelism:            ON" in config:
        raise RuntimeError(
            "MPI-enabled turbulence regression requires an available MPIEXEC "
            "launcher"
        )

    cases.append((
        "2d_np1", 1,
        (
            "mesh/nx3=1", "meshblock/nx3=1",
            "particles/pic_enable_2d3v=true",
        ),
    ))
    for label, nproc, overrides in cases:
        basename = "PICTurbulentDynamoSmoke_" + label
        _remove_outputs(basename)
        _run(basename, nproc, overrides)
        _RESULTS[label] = _measure(basename)

    mhd_cases = [("np1", 1)]
    if mpi_available:
        mhd_cases.append(("np2", 2))
    for label, nproc in mhd_cases:
        basename = "TurbDriverMHD_" + label
        _remove_outputs(basename)
        _run_mhd(basename, nproc)
        _MHD_RESULTS[label] = _measure_mhd(basename)
    _MHD_RESTART.update(_continue_mhd_from_cycle_one("TurbDriverMHD_np1"))
    _check_legacy_bounds_upgrade()
    _expect_mhd_failure(
        "turb_driving/min_kx=0", "complete signed [-nhigh, nhigh] range"
    )
    _expect_mhd_failure(
        "turb_driving/dedt=-1", "require finite dedt >= 0"
    )
    _expect_mhd_failure(
        "turb_driving/driving_type=1", "only driving_type=0 is supported"
    )


def analyze():
    if _SKIPPED:
        return True

    ok = True
    for label, result in _RESULTS.items():
        max_history_bmean = float(np.max(np.abs(
            result["history_magnetic_mean"])))
        max_stored_bmean = float(np.max(np.abs(
            result["stored_magnetic_mean"])))
        max_momentum = float(np.max(np.abs(result["momentum_sum"])))
        max_current = float(np.max(np.abs(result["current_sum"])))
        max_history_cr_momentum = float(np.max(np.abs(
            result["history_cr_momentum"])))
        logger.info(
            "%s Brms=% .8e history_divB=% .8e stored_divB=% .8e",
            label, result["magnetic_rms"], result["history_max_divb"],
            result["stored_max_divb"],
        )
        logger.info(
            "%s particle_count=%d species_counts=%s position_spread=% .8e",
            label, result["particle_count"], result["species_counts"],
            result["position_spread"],
        )
        logger.info(
            "%s max_history_Bmean=% .8e max_stored_Bmean=% .8e "
            "max_momentum=% .8e max_current=% .8e max_history_cr_P=% .8e",
            label, max_history_bmean, max_stored_bmean, max_momentum,
            max_current, max_history_cr_momentum,
        )
        logger.info(
            "%s forcing_delta_E=% .16e expected=% .16e residual=% .8e "
            "spectral_div=% .8e outside_shell=% .8e mixed_sign=%s",
            label, result["forcing_energy_increment"],
            result["forcing_energy_expected"],
            result["forcing_energy_residual"], result["divergence_ratio"],
            result["outside_shell_fraction"],
            result["mixed_sign_relative_power"],
        )
        ok = (abs(result["magnetic_rms"] - 1.0e-3) < 1.0e-9) and ok
        ok = (max_history_bmean < 1.0e-12) and ok
        ok = (max_stored_bmean < 1.0e-9) and ok
        ok = (result["history_max_divb"] < 1.0e-12) and ok
        ok = (result["stored_max_divb"] < 1.0e-12) and ok
        expected_particles_per_species = result["cell_count"]
        ok = (result["particle_count"] ==
              expected_particles_per_species*6) and ok
        ok = np.array_equal(
            result["species_counts"],
            np.full(6, expected_particles_per_species),
        ) and ok
        ok = (result["position_spread"] < 1.0e-12) and ok
        ok = (max_momentum < 1.0e-12) and ok
        ok = (max_current < 1.0e-10) and ok
        ok = (max_history_cr_momentum < 1.0e-12) and ok
        energy_tolerance = max(
            2.0e-9, 5.0e-6*abs(result["forcing_energy_expected"])
        )
        ok = (abs(result["forcing_energy_residual"]) < energy_tolerance) and ok
        # Mesh binary output is float32, so FFT residuals are limited accordingly.
        ok = (result["divergence_ratio"] < 2.0e-7) and ok
        ok = (result["outside_shell_fraction"] < 2.0e-14) and ok
        ok = bool(np.all(result["mixed_sign_relative_power"] > 1.0e-8)) and ok
        if label == "2d_np1":
            ok = (result["out_of_plane_fraction"] > 0.05) and ok

    if "np2" in _RESULTS:
        error = float(np.max(np.abs(
            _RESULTS["np2"]["magnetic"] - _RESULTS["np1"]["magnetic"]
        )))
        logger.info("MPI magnetic decomposition error=% .8e", error)
        ok = (error < 2.0e-12) and ok
        force_error = float(np.max(np.abs(
            _RESULTS["np2"]["force"] - _RESULTS["np1"]["force"]
        )))
        energy_residual_error = abs(
            _RESULTS["np2"]["forcing_energy_residual"] -
            _RESULTS["np1"]["forcing_energy_residual"]
        )
        logger.info(
            "MPI force error=% .8e energy-residual error=% .8e",
            force_error, energy_residual_error,
        )
        ok = (force_error < 2.0e-11) and ok
        ok = (energy_residual_error < 2.0e-11) and ok

    for label, result in _MHD_RESULTS.items():
        logger.info(
            "MHD %s forcing_delta_E=% .16e expected=% .16e residual=% .8e "
            "spectral_div=% .8e outside_shell=% .8e mixed_sign=%s",
            label, result["forcing_energy_increment"],
            result["forcing_energy_expected"],
            result["forcing_energy_residual"], result["divergence_ratio"],
            result["outside_shell_fraction"],
            result["mixed_sign_relative_power"],
        )
        energy_tolerance = max(
            2.0e-9, 5.0e-6*abs(result["forcing_energy_expected"])
        )
        ok = (abs(result["forcing_energy_residual"]) < energy_tolerance) and ok
        ok = (result["divergence_ratio"] < 2.0e-7) and ok
        ok = (result["outside_shell_fraction"] < 2.0e-14) and ok
        ok = bool(np.all(result["mixed_sign_relative_power"] > 1.0e-8)) and ok

    if "np2" in _MHD_RESULTS:
        mhd_force_error = float(np.max(np.abs(
            _MHD_RESULTS["np2"]["force"] - _MHD_RESULTS["np1"]["force"]
        )))
        mhd_energy_error = abs(
            _MHD_RESULTS["np2"]["forcing_energy_residual"] -
            _MHD_RESULTS["np1"]["forcing_energy_residual"]
        )
        logger.info(
            "MHD MPI force error=% .8e energy-residual error=% .8e",
            mhd_force_error, mhd_energy_error,
        )
        ok = (mhd_force_error < 2.0e-11) and ok
        ok = (mhd_energy_error < 2.0e-11) and ok
    restart_force_error = float(np.max(np.abs(
        _MHD_RESTART["force"] - _MHD_RESULTS["np1"]["force"]
    )))
    restart_history_error = float(np.max(np.abs(
        _MHD_RESTART["final_history"] -
        _MHD_RESULTS["np1"]["final_history"]
    )))
    logger.info(
        "MHD restart force error=% .8e final-history error=% .8e",
        restart_force_error, restart_history_error,
    )
    ok = (restart_force_error == 0.0) and ok
    ok = (restart_history_error < 1.0e-14) and ok
    return ok
