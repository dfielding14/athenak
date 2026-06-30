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
import scripts.utils.athena as athena

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
_MPIEXEC = shlex.split(os.environ.get("MPIEXEC", "mpiexec"))
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_RESULTS = {}
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
    )
    for pattern in patterns:
        for path in glob.glob(pattern):
            os.remove(path)


def _mpi_launcher_available():
    return bool(_MPIEXEC and shutil.which(_MPIEXEC[0]))


def _run(basename, nproc):
    input_path = os.path.join(_SOURCE_ROOT, "inputs", _INPUT)
    command = [
        "./athena", "-i", input_path, "job/basename=" + basename,
        "time/nlim=0", "time/tlim=0.0",
    ]
    if nproc > 1:
        command = _MPIEXEC + ["-n", str(nproc)] + command
    proc = subprocess.run(
        command, cwd=_exe_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        output = (proc.stdout or "") + (proc.stderr or "")
        raise RuntimeError("PIC turbulent-dynamo smoke failed:\n" + output)


def _latest(directory, pattern):
    paths = sorted(glob.glob(os.path.join(_exe_dir(), directory, pattern)))
    if not paths:
        raise RuntimeError("Missing output for " + pattern)
    return paths[-1]


def _first_history_row(basename):
    path = os.path.join(_exe_dir(), basename + ".user.hst")
    with open(path, "r", encoding="utf-8") as stream:
        rows = [line for line in stream if line.strip() and not line.startswith("#")]
    if not rows:
        raise RuntimeError("Missing user-history rows for " + basename)
    return np.fromstring(rows[0], sep=" ")


def _measure(basename):
    field_data = bin_convert.read_binary_as_athdf(
        _latest("bin", basename + ".mhd_w_bcc.*.bin")
    )
    divb_data = bin_convert.read_binary_as_athdf(
        _latest("bin", basename + ".mhd_divb.*.bin")
    )
    magnetic = np.stack([
        np.asarray(field_data[name], dtype=float)
        for name in ("bcc1", "bcc2", "bcc3")
    ])
    particles = read_particle_vtk(
        _latest("pvtk", basename + ".prtcl_all.*.part.vtk")
    )
    species = particles.scalars["species"]
    states = particles.vectors["vel"]
    group_positions = particles.points.reshape(-1, 6, 3)
    history = _first_history_row(basename)
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
        "species_counts": np.bincount(species, minlength=6),
        "position_spread": float(np.max(np.ptp(group_positions, axis=1))),
        "momentum_sum": np.sum(states, axis=0),
        "current_sum": 100.0*np.sum(states, axis=0),
        "history_cr_momentum": history[11:14],
    }


def run(**kwargs):
    global _SKIPPED
    logger.debug("Running test " + __name__)
    config = _configuration()
    if "Problem generator:          turb" not in config:
        logger.info("Skipping: configure this test with -DPROBLEM=turb")
        _SKIPPED = True
        return

    cases = [("np1", 1)]
    if ("MPI parallelism:            ON" in config
            and _mpi_launcher_available()):
        cases.append(("np2", 2))
    elif "MPI parallelism:            ON" in config:
        logger.info("Skipping MPI case: launcher %s not found",
                    _MPIEXEC[0] if _MPIEXEC else "<empty>")

    for label, nproc in cases:
        basename = "PICTurbulentDynamoSmoke_" + label
        _remove_outputs(basename)
        _run(basename, nproc)
        _RESULTS[label] = _measure(basename)


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
        ok = (abs(result["magnetic_rms"] - 1.0e-3) < 1.0e-9) and ok
        ok = (max_history_bmean < 1.0e-12) and ok
        ok = (max_stored_bmean < 1.0e-9) and ok
        ok = (result["history_max_divb"] < 1.0e-12) and ok
        ok = (result["stored_max_divb"] < 1.0e-12) and ok
        ok = (result["particle_count"] == 16**3*6) and ok
        ok = np.array_equal(result["species_counts"], np.full(6, 16**3)) and ok
        ok = (result["position_spread"] < 1.0e-12) and ok
        ok = (max_momentum < 1.0e-12) and ok
        ok = (max_current < 1.0e-10) and ok
        ok = (max_history_cr_momentum < 1.0e-12) and ok

    if "np2" in _RESULTS:
        error = float(np.max(np.abs(
            _RESULTS["np2"]["magnetic"] - _RESULTS["np1"]["magnetic"]
        )))
        logger.info("MPI magnetic decomposition error=% .8e", error)
        ok = (error < 2.0e-12) and ok
    return ok
