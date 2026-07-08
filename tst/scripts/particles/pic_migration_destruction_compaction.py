"""Serial/MPI survivor identity for simultaneous particle sends and destruction."""

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
import scripts.utils.athena as athena  # noqa: E402


_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_PUBLICATION_DIR = os.path.join(_SOURCE_ROOT, "tst", "publication")
sys.path.insert(0, _PUBLICATION_DIR)
from pvtk_particles import ParticleVTKData, read_particle_vtk  # noqa: E402

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_migration_destruction_compaction.athinput"
_MPIEXEC = shlex.split(os.environ.get("MPIEXEC", "mpiexec"))
_REMOVED_TAGS = {
    0, 1, 12, 13, 24, 25, 36, 37,
    155, 167, 179, 191,
}
_CROSS_RANK_TAGS = {
    59, 71, 83, 95,
    96, 97, 108, 109, 120, 121, 132, 133,
}
_INITIAL_PARTICLE_COUNT = 192
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        "ATHENA_PIC_MIGRATION_COMPACTION_EXE_DIR",
        os.path.join(os.getcwd(), "build", "src"),
    )


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


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


def _mpi_prefix(nproc):
    if not _MPIEXEC:
        raise RuntimeError("MPIEXEC must not be empty")
    return [*_MPIEXEC, "-n", str(nproc)]


def _remove_outputs(basename):
    pattern = os.path.join(
        _athena_exe_dir(), "pvtk", basename + ".*.part.vtk"
    )
    for path in glob.glob(pattern):
        if os.path.isfile(path):
            os.remove(path)


def _particle_cycles(basename):
    pattern = os.path.join(
        _athena_exe_dir(), "pvtk", basename + ".prtcl_all.*.part.vtk"
    )
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise RuntimeError("No particle VTK files found for " + basename)
    filename_pattern = re.compile(
        re.escape(basename)
        + r"\.prtcl_all\.(?:[0-9]+\.)?([0-9]{5})\.part\.vtk"
    )
    cycles = {}
    for path in paths:
        match = filename_pattern.fullmatch(os.path.basename(path))
        if match is None:
            raise RuntimeError("Malformed particle VTK filename: " + path)
        cycles.setdefault(int(match.group(1)), []).append(path)
    return cycles


def _merge_particle_slices(paths):
    slices = [read_particle_vtk(path) for path in sorted(paths)]
    scalar_names = set(slices[0].scalars)
    vector_names = set(slices[0].vectors)
    for data in slices[1:]:
        if set(data.scalars) != scalar_names or set(data.vectors) != vector_names:
            raise RuntimeError("Particle VTK schemas differ across GID slices")
    return ParticleVTKData(
        points=np.concatenate([data.points for data in slices], axis=0),
        scalars={
            name: np.concatenate([data.scalars[name] for data in slices])
            for name in scalar_names
        },
        vectors={
            name: np.concatenate([data.vectors[name] for data in slices], axis=0)
            for name in vector_names
        },
    )


def _snapshot(paths):
    data = _merge_particle_slices(paths)
    required_scalars = {
        "birth_time",
        "cr_source",
        "deltaf_f0",
        "deltaf_weight",
        "gid",
        "macro_weight",
        "ptag",
        "species",
    }
    missing = sorted(required_scalars - set(data.scalars))
    if missing:
        raise RuntimeError("Particle VTK is missing scalars " + repr(missing))
    if "vel" not in data.vectors:
        raise RuntimeError("Particle VTK is missing vector: vel")
    tags = data.scalars["ptag"]
    if np.unique(tags).size != tags.size:
        raise RuntimeError("Particle VTK tags are not unique")
    order = np.argsort(tags)
    return {
        "points": data.points[order],
        "vel": data.vectors["vel"][order],
        **{
            name: data.scalars[name][order]
            for name in sorted(required_scalars)
        },
    }


def _run_case(basename, nproc):
    _remove_outputs(basename)
    command = [
        "./athena",
        "-i",
        _athena_input_path(),
        "job/basename=" + basename,
    ]
    if nproc > 1:
        command = [*_mpi_prefix(nproc), *command]
    logger.info("Executing %s", " ".join(command))
    proc = subprocess.run(
        command,
        cwd=_athena_exe_dir(),
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "Command failed for " + basename + "\n"
            + (proc.stdout or "") + (proc.stderr or "")
        )
    cycles = _particle_cycles(basename)
    if len(cycles) < 2:
        raise RuntimeError("Expected initial and final particle outputs")
    return {
        "initial": _snapshot(cycles[min(cycles)]),
        "final": _snapshot(cycles[max(cycles)]),
    }


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    _RESULTS["serial"] = _run_case("pic_compaction_serial", 1)
    mpi_enabled = _athena_mpi_enabled()
    launcher_available = bool(_MPIEXEC) and shutil.which(_MPIEXEC[0]) is not None
    _RESULTS["mpi_enabled"] = mpi_enabled
    if mpi_enabled and not launcher_available:
        raise RuntimeError(
            "MPI-enabled compaction regression requires an available MPIEXEC "
            "launcher"
        )
    if mpi_enabled:
        _RESULTS["mpi2"] = _run_case("pic_compaction_mpi2", 2)
    else:
        logger.info("Skipping MPI2 compaction case: Athena MPI support is disabled")


def _tag_to_value(snapshot, field):
    return {
        int(tag): value
        for tag, value in zip(snapshot["ptag"], snapshot[field])
    }


def analyze():
    logger.debug("Analyzing test " + __name__)
    expected_initial = set(range(_INITIAL_PARTICLE_COUNT))
    expected_final = expected_initial - _REMOVED_TAGS
    serial = _RESULTS["serial"]
    serial_initial_tags = set(int(tag) for tag in serial["initial"]["ptag"])
    serial_final_tags = set(int(tag) for tag in serial["final"]["ptag"])
    serial_ok = (
        serial_initial_tags == expected_initial
        and serial_final_tags == expected_final
    )
    if not _RESULTS["mpi_enabled"]:
        return serial_ok

    mpi2 = _RESULTS["mpi2"]
    mpi_initial_tags = set(int(tag) for tag in mpi2["initial"]["ptag"])
    mpi_final_tags = set(int(tag) for tag in mpi2["final"]["ptag"])
    metadata_equal = all(
        np.array_equal(serial["final"][name], mpi2["final"][name])
        for name in ["ptag", "gid", "species", "cr_source"]
    )
    state_error = max(
        float(np.max(np.abs(serial["final"][name] - mpi2["final"][name])))
        for name in [
            "points",
            "vel",
            "birth_time",
            "deltaf_f0",
            "deltaf_weight",
            "macro_weight",
        ]
    )

    initial_gid = _tag_to_value(mpi2["initial"], "gid")
    final_gid = _tag_to_value(mpi2["final"], "gid")
    cross_rank_observed = all(
        (int(initial_gid[tag]) < 2) != (int(final_gid[tag]) < 2)
        for tag in _CROSS_RANK_TAGS
    )
    sends_per_rank = [0, 0]
    destructions_per_rank = [0, 0]
    for tag in _CROSS_RANK_TAGS:
        sends_per_rank[0 if int(initial_gid[tag]) < 2 else 1] += 1
    for tag in _REMOVED_TAGS:
        destructions_per_rank[0 if int(initial_gid[tag]) < 2 else 1] += 1
    simultaneous_on_each_rank = all(
        sends_per_rank[rank] > 0 and destructions_per_rank[rank] > 0
        for rank in range(2)
    )

    metrics = {
        "serial_survivors_exact": serial_ok,
        "mpi_initial_tags_exact": mpi_initial_tags == expected_initial,
        "mpi_survivors_exact": mpi_final_tags == expected_final,
        "serial_mpi_metadata_equal": metadata_equal,
        "serial_mpi_state_max_error": state_error,
        "cross_rank_tags": sorted(_CROSS_RANK_TAGS),
        "sends_per_rank": sends_per_rank,
        "destructions_per_rank": destructions_per_rank,
        "simultaneous_send_destroy_each_rank": simultaneous_on_each_rank,
    }
    logger.info("Particle compaction metrics: %s", metrics)
    return (
        serial_ok
        and mpi_initial_tags == expected_initial
        and mpi_final_tags == expected_final
        and metadata_equal
        and state_error <= 1.0e-6
        and cross_rank_observed
        and simultaneous_on_each_rank
    )
