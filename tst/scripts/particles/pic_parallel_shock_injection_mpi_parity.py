"""Bounded serial/MPI2 deterministic pic_parallel_shock injection parity."""

from __future__ import annotations

import glob
import json
import logging
import os
import shlex
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena

_PUBLICATION_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "publication")
)
sys.path.insert(0, _PUBLICATION_DIR)
from pvtk_particles import read_particle_vtk  # noqa: E402

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_parallel_shock_injection_mpi_parity.athinput"
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_MPI_LAUNCHER_ENV = "ATHENA_PIC_PARALLEL_SHOCK_MPI_PARITY_MPI_LAUNCHER"
_RESTART_FLOAT_FIELDS = [
    "ps_mass_reservoir_global",
    "ps_injected_cr_count_global",
    "ps_injected_cr_mass_global",
    "ps_injected_cr_momentum_x1_global",
    "ps_injected_cr_momentum_x2_global",
    "ps_injected_cr_momentum_x3_global",
    "ps_injected_cr_energy_global",
]
_LEDGER_FIELDS = _RESTART_FLOAT_FIELDS[1:]
_REQUIRED_PARTICLE_SCALARS = {
    "gid",
    "ptag",
    "species",
    "cr_source",
    "macro_weight",
    "birth_time",
}
_MACRO_MASS = 1.0e-3
_METADATA_ABS_TOL = 2.0e-12
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        "ATHENA_PIC_PARALLEL_SHOCK_MPI_PARITY_EXE_DIR",
        os.path.join(os.getcwd(), "build", "src"),
    )


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


def _mpi_launcher(nproc):
    raw_launcher = os.environ.get(
        _MPI_LAUNCHER_ENV, os.environ.get("MPIEXEC", "mpiexec")
    )
    command = shlex.split(raw_launcher)
    if not command:
        raise RuntimeError(_MPI_LAUNCHER_ENV + " must not be empty")
    replaced = False
    for index, value in enumerate(command):
        if "{nproc}" in value:
            command[index] = value.replace("{nproc}", str(nproc))
            replaced = True
    if not replaced:
        command += ["-n", str(nproc)]
    return command


def _remove_outputs(basename):
    for dirname in ["pvtk", "rst"]:
        pattern = os.path.join(_athena_exe_dir(), dirname, basename + ".*")
        for path in glob.glob(pattern):
            if os.path.isfile(path):
                os.remove(path)


def _latest_file(dirname, pattern, label):
    paths = sorted(glob.glob(os.path.join(_athena_exe_dir(), dirname, pattern)))
    if not paths:
        raise RuntimeError("No " + label + " files found for pattern: " + pattern)
    return paths[-1]


def _restart_problem_parameters(path):
    with open(path, "rb") as handle:
        text = handle.read(65536).decode("latin1", errors="ignore")
    if "<par_end>" not in text:
        raise RuntimeError("Restart parameter header is missing <par_end>: " + path)

    active_block = None
    parameters = {}
    for raw_line in text.split("<par_end>", 1)[0].splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("<") and line.endswith(">"):
            active_block = line[1:-1]
            continue
        if active_block != "problem" or "=" not in line:
            continue
        key, value = line.split("=", 1)
        parameters[key.strip()] = value.split("#", 1)[0].strip()
    return parameters


def _restart_metadata(basename):
    path = _latest_file("rst", basename + ".*.rst", "restart")
    parameters = _restart_problem_parameters(path)
    required = ["ps_next_tag", *_RESTART_FLOAT_FIELDS]
    missing = [name for name in required if name not in parameters]
    if missing:
        raise RuntimeError(
            "Restart metadata is missing injection-accounting keys "
            + repr(missing)
            + ": "
            + path
        )
    return {
        "ps_next_tag": int(parameters["ps_next_tag"]),
        **{name: float(parameters[name]) for name in _RESTART_FLOAT_FIELDS},
    }


def _particle_snapshot(basename):
    path = _latest_file(
        "pvtk", basename + ".prtcl_all.*.part.vtk", "particle VTK"
    )
    data = read_particle_vtk(path)
    missing = sorted(_REQUIRED_PARTICLE_SCALARS - set(data.scalars))
    if missing:
        raise RuntimeError("Particle VTK is missing scalars " + repr(missing))
    if "vel" not in data.vectors:
        raise RuntimeError("Particle VTK is missing vector: vel")

    tags = data.scalars["ptag"]
    if np.unique(tags).size != tags.size:
        raise RuntimeError("Particle VTK ptag values must be unique")
    order = np.argsort(tags)
    return {
        "path": path,
        "points": np.asarray(data.points)[order],
        "scalars": {
            name: np.asarray(values)[order] for name, values in data.scalars.items()
        },
        "vectors": {
            name: np.asarray(values)[order] for name, values in data.vectors.items()
        },
    }


def _run_case(name, nproc):
    basename = "pic_parallel_shock_injection_parity_" + name
    _remove_outputs(basename)
    command = ["./athena", "-i", _athena_input_path(), "job/basename=" + basename]
    if nproc > 1:
        command = _mpi_launcher(nproc) + command
    logger.info("Executing %s: %s", name, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Command failed for " + name + "\n" + output)
    return {
        "basename": basename,
        "metadata": _restart_metadata(basename),
        "particles": _particle_snapshot(basename),
    }


def _array_map_parity(left, right):
    if set(left) != set(right):
        return False
    return all(np.array_equal(left[name], right[name]) for name in left)


def _metadata_absolute_errors(serial, mpi2):
    return {
        name: abs(serial[name] - mpi2[name]) for name in _RESTART_FLOAT_FIELDS
    }


def _case_invariants(case):
    metadata = case["metadata"]
    scalars = case["particles"]["scalars"]
    count = int(scalars["ptag"].size)
    next_tag = metadata["ps_next_tag"]
    ledger = np.array([metadata[name] for name in _LEDGER_FIELDS])
    return {
        "particle_count": count,
        "next_tag": next_tag,
        "particle_count_matches_next_tag": count == next_tag,
        "tags_are_dense_from_zero": np.array_equal(
            scalars["ptag"], np.arange(next_tag, dtype=np.int64)
        ),
        "all_particles_are_shock_injected": bool(np.all(scalars["cr_source"] == 1)),
        "all_particles_are_species_zero": bool(np.all(scalars["species"] == 0)),
        "all_macro_weights_are_one": bool(np.all(scalars["macro_weight"] == 1.0)),
        "all_birth_times_are_zero": bool(np.all(scalars["birth_time"] == 0.0)),
        "ledger_is_finite": bool(np.all(np.isfinite(ledger))),
        "ledger_count_matches_next_tag": ledger[0] == next_tag,
        "ledger_mass_matches_macro_mass": (
            abs(ledger[1] - next_tag * _MACRO_MASS) <= _METADATA_ABS_TOL
        ),
        "reservoir_is_bounded": (
            0.0 <= metadata["ps_mass_reservoir_global"] < _MACRO_MASS
        ),
    }


def _summary():
    serial = _RESULTS["serial"]
    mpi2 = _RESULTS["mpi2"]
    serial_metadata = serial["metadata"]
    mpi2_metadata = mpi2["metadata"]
    serial_particles = serial["particles"]
    mpi2_particles = mpi2["particles"]
    metadata_errors = _metadata_absolute_errors(serial_metadata, mpi2_metadata)
    return {
        "evidence_class": "bounded_local_serial_vs_mpi2_host_regression",
        "not_qualification_evidence": True,
        "physical_setup": (
            "same_input_deck_same_two_meshblocks_zero_startup_particles_one_rk1_cycle"
        ),
        "serial": _case_invariants(serial),
        "mpi2": _case_invariants(mpi2),
        "next_tag_equal": serial_metadata["ps_next_tag"] == mpi2_metadata["ps_next_tag"],
        "metadata_absolute_errors": metadata_errors,
        "metadata_within_tolerance": all(
            error <= _METADATA_ABS_TOL for error in metadata_errors.values()
        ),
        "particle_points_equal": np.array_equal(
            serial_particles["points"], mpi2_particles["points"]
        ),
        "particle_scalars_equal": _array_map_parity(
            serial_particles["scalars"], mpi2_particles["scalars"]
        ),
        "particle_vectors_equal": _array_map_parity(
            serial_particles["vectors"], mpi2_particles["vectors"]
        ),
    }


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    _RESULTS["serial"] = _run_case("serial", 1)
    _RESULTS["mpi2"] = _run_case("mpi2", 2)


def analyze():
    logger.debug("Analyzing test " + __name__)
    summary = _summary()
    logger.info("PIC parallel-shock injection MPI parity metrics: %s", summary)
    return (
        summary["serial"]["particle_count"] > 0
        and all(
            value
            for name, value in summary["serial"].items()
            if name not in {"particle_count", "next_tag"}
        )
        and all(
            value
            for name, value in summary["mpi2"].items()
            if name not in {"particle_count", "next_tag"}
        )
        and summary["next_tag_equal"]
        and summary["metadata_within_tolerance"]
        and summary["particle_points_equal"]
        and summary["particle_scalars_equal"]
        and summary["particle_vectors_equal"]
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
