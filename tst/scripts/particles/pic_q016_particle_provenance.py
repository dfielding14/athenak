"""Bounded Q-016 provenance, restart, migration, and spectra regression."""

from __future__ import annotations

import glob
import logging
import os
import re
import shlex
import struct
import subprocess
import sys

import numpy as np

_PUBLICATION_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "publication")
)
sys.path.insert(0, _PUBLICATION_DIR)
from pvtk_particles import ParticleVTKData, read_particle_vtk  # noqa: E402
from q016_particle_spectra import build_weighted_species_cohort_spectra  # noqa: E402

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_q016_particle_provenance.athinput"
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_PIC_RESTART_MAGIC = 0x5049435253543031
_EXPECTED_RESTART_SCHEMA = 7
_INT_SCALARS = ["gid", "ptag", "species", "cr_source"]
_FLOAT_SCALARS = [
    "macro_weight",
    "birth_time",
    "deltaf_f0",
    "deltaf_weight",
]
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        "ATHENA_Q016_EXE_DIR", os.path.join(os.getcwd(), "build", "src")
    )


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


def _rank_count():
    raw = os.environ.get("ATHENA_Q016_NPROC", "1")
    try:
        nproc = int(raw)
    except ValueError as error:
        raise RuntimeError("ATHENA_Q016_NPROC must be a positive integer") from error
    if nproc < 1:
        raise RuntimeError("ATHENA_Q016_NPROC must be a positive integer")
    return nproc


def _launcher_prefix():
    nproc = _rank_count()
    if nproc == 1:
        return []
    raw = os.environ.get(
        "ATHENA_Q016_LAUNCHER", os.environ.get("MPIEXEC", "mpiexec")
    )
    launcher = shlex.split(raw)
    if not launcher:
        raise RuntimeError("ATHENA_Q016_LAUNCHER must not be empty for multi-rank runs")
    return launcher + ["-n", str(nproc)]


def _athena_mpi_enabled():
    proc = subprocess.run(
        ["./athena", "-c"], cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError("Unable to query Athena configuration with -c")
    output = (proc.stdout or "") + (proc.stderr or "")
    return "MPI parallelism:            ON" in output


def _require_supported_rank_configuration():
    if _rank_count() > 1 and not _athena_mpi_enabled():
        raise RuntimeError(
            "ATHENA_Q016_NPROC > 1 requires an MPI-enabled Athena executable"
        )


def _remove_outputs(basename):
    exe_dir = _athena_exe_dir()
    for pattern in [
        os.path.join(exe_dir, "pvtk", basename + ".*.part.vtk"),
        os.path.join(exe_dir, "rst", basename + ".*.rst*"),
    ]:
        for fname in glob.glob(pattern):
            if os.path.isfile(fname):
                os.remove(fname)


def _run_athena(label, arguments, restart_file=None):
    _require_supported_rank_configuration()
    command = _launcher_prefix() + ["./athena"]
    if restart_file is None:
        command += ["-i", _athena_input_path()]
    else:
        command += ["-r", restart_file]
    command += list(arguments)
    logger.info("Executing %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "Command failed for " + label + "\n" + (proc.stdout or "") + (proc.stderr or "")
        )


def _matching_files(pattern, label):
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise RuntimeError("No " + label + " files found for pattern: " + pattern)
    return matches


def _particle_vtk_cycles(basename):
    paths = _matching_files(
        os.path.join(_athena_exe_dir(), "pvtk", basename + ".prtcl_all.*.part.vtk"),
        "particle VTK",
    )
    pattern = re.compile(
        re.escape(basename) + r"\.prtcl_all\.(?:[0-9]+\.)?([0-9]{5})\.part\.vtk"
    )
    cycles = {}
    for path in paths:
        match = pattern.fullmatch(os.path.basename(path))
        if match is None:
            raise RuntimeError("Malformed particle VTK filename: " + path)
        cycle = int(match.group(1))
        cycles.setdefault(cycle, []).append(path)
    return {cycle: sorted(partitions) for cycle, partitions in cycles.items()}


def _latest_particle_vtk_partitions(basename):
    cycles = _particle_vtk_cycles(basename)
    return cycles[max(cycles)]


def _latest_restart_file(basename):
    return _matching_files(
        os.path.join(_athena_exe_dir(), "rst", basename + ".*.rst"), "restart"
    )[-1]


def _restart_schema(path):
    with open(path, "rb") as fp:
        data = fp.read()
    for byte_order in ["<", ">"]:
        marker = struct.pack(byte_order + "Q", _PIC_RESTART_MAGIC)
        offset = data.find(marker)
        if offset >= 0:
            return struct.unpack_from(byte_order + "i", data, offset + len(marker))[0]
    raise RuntimeError("Particle restart marker not found in " + path)


def _merge_particle_vtk_partitions(paths):
    if not paths:
        raise RuntimeError("Particle VTK partition list is empty")
    partitions = [read_particle_vtk(path) for path in paths]
    scalar_schema = {
        name: values.dtype.str for name, values in partitions[0].scalars.items()
    }
    vector_schema = {
        name: values.dtype.str for name, values in partitions[0].vectors.items()
    }
    for data in partitions[1:]:
        if {
            name: values.dtype.str for name, values in data.scalars.items()
        } != scalar_schema:
            raise RuntimeError("Particle VTK scalar schemas differ across partitions")
        if {
            name: values.dtype.str for name, values in data.vectors.items()
        } != vector_schema:
            raise RuntimeError("Particle VTK vector schemas differ across partitions")
    return ParticleVTKData(
        points=np.concatenate([data.points for data in partitions], axis=0),
        scalars={
            name: np.concatenate([data.scalars[name] for data in partitions])
            for name in scalar_schema
        },
        vectors={
            name: np.concatenate([data.vectors[name] for data in partitions], axis=0)
            for name in vector_schema
        },
    )


def _snapshot(paths):
    data = _merge_particle_vtk_partitions(paths)
    for name in _INT_SCALARS + _FLOAT_SCALARS:
        if name not in data.scalars:
            raise RuntimeError("Missing particle VTK scalar: " + name)
    if "vel" not in data.vectors:
        raise RuntimeError("Missing particle VTK vector: vel")
    for name in _INT_SCALARS:
        if data.scalars[name].dtype.kind != "i":
            raise RuntimeError("Expected typed integer scalar: " + name)

    order = np.argsort(data.scalars["ptag"])
    if np.unique(data.scalars["ptag"]).size != data.scalars["ptag"].size:
        raise RuntimeError("Particle VTK ptag values must be unique across partitions")
    return {
        "points": data.points[order],
        "vel": data.vectors["vel"][order],
        **{name: data.scalars[name][order] for name in _INT_SCALARS + _FLOAT_SCALARS},
    }


def _manual_histogram(speed, weights, edges, mask):
    histogram = np.zeros(edges.size - 1, dtype=np.float64)
    for value, weight in zip(speed[mask], weights[mask]):
        for index in range(edges.size - 1):
            if edges[index] <= value < edges[index + 1]:
                histogram[index] += weight
                break
            if index == edges.size - 2 and value == edges[index + 1]:
                histogram[index] += weight
                break
    return histogram


def _spectrum_agreement(paths, semantics):
    data = _merge_particle_vtk_partitions(paths)
    speed = np.linalg.norm(data.vectors["vel"], axis=1)
    max_speed = max(1.0, float(np.max(speed)))
    edges = np.linspace(0.0, max_speed + 1.0, 9)
    birth_time = data.scalars["birth_time"]
    max_birth = max(1.0e-12, float(np.max(birth_time)))
    birth_edges = np.linspace(-1.0e-12, max_birth + 1.0e-12, 5)
    payload = build_weighted_species_cohort_spectra(
        data, edges, birth_edges, delta_f_semantics=semantics
    )

    weights = data.scalars["macro_weight"].astype(np.float64)
    if semantics == "delta_f_perturbation":
        weights *= data.scalars["deltaf_weight"]
    all_mask = np.ones(speed.size, dtype=bool)
    manual = _manual_histogram(speed, weights, edges, all_mask)
    errors = [
        float(
            np.max(
                np.abs(
                    manual
                    - np.asarray(payload["all_particles"]["weighted_sum_in_bins"])
                )
            )
        )
    ]
    grouped_count = 0
    for group in payload["groups"]:
        mask = data.scalars["species"] == group["species"]
        mask &= data.scalars["cr_source"] == group["cr_source"]
        mask &= birth_time >= group["birth_time_min"]
        if group["birth_cohort"] + 1 == len(birth_edges) - 1:
            mask &= birth_time <= group["birth_time_max"]
        else:
            mask &= birth_time < group["birth_time_max"]
        grouped_count += int(np.count_nonzero(mask))
        manual = _manual_histogram(speed, weights, edges, mask)
        errors.append(
            float(
                np.max(
                    np.abs(
                        manual
                        - np.asarray(group["spectrum"]["weighted_sum_in_bins"])
                    )
                )
            )
        )
    if grouped_count != speed.size:
        raise RuntimeError("Resolved spectra did not cover every particle")
    return max(errors), payload


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    _RESULTS["configured_ranks"] = _rank_count()
    full = "pic_q016_full"
    segment = "pic_q016_segment"
    restarted = "pic_q016_restart"
    for basename in [full, segment, restarted]:
        _remove_outputs(basename)

    _run_athena(
        "full",
        ["job/basename=" + full, "time/nlim=8", "output2/dcycle=0"],
    )
    full_cycles = _particle_vtk_cycles(full)
    _RESULTS["full_initial"] = _snapshot(full_cycles[min(full_cycles)])
    _RESULTS["full_final"] = _snapshot(full_cycles[max(full_cycles)])
    _RESULTS["full_final_paths"] = full_cycles[max(full_cycles)]

    _run_athena("segment", ["job/basename=" + segment, "time/nlim=4"])
    restart_path = _latest_restart_file(segment)
    _RESULTS["restart_schema"] = _restart_schema(restart_path)
    _run_athena(
        "restart",
        [
            "job/basename=" + restarted,
            "time/nlim=8",
            "output1/file_number=0",
            "output2/dcycle=0",
        ],
        restart_file=os.path.relpath(restart_path, _athena_exe_dir()),
    )
    _RESULTS["restart_final"] = _snapshot(_latest_particle_vtk_partitions(restarted))


def analyze():
    logger.debug("Analyzing test " + __name__)
    initial = _RESULTS["full_initial"]
    full = _RESULTS["full_final"]
    restarted = _RESULTS["restart_final"]

    int_equal = {
        name: bool(np.array_equal(full[name], restarted[name])) for name in _INT_SCALARS
    }
    float_errors = {
        name: float(np.max(np.abs(full[name] - restarted[name])))
        for name in ["points", "vel"] + _FLOAT_SCALARS
    }

    initial_gid = {int(tag): int(gid) for tag, gid in zip(initial["ptag"], initial["gid"])}
    migrated = any(
        tag in initial_gid and initial_gid[tag] != int(gid)
        for tag, gid in zip(full["ptag"], full["gid"])
    )
    sources = set(int(value) for value in full["cr_source"])
    shock_mask = full["cr_source"] == 1
    source_metadata_ok = (
        sources == {0, 1}
        and np.all(full["species"][shock_mask] == 1)
        and np.all(full["birth_time"][full["cr_source"] == 0] == 0.0)
        and np.all(full["birth_time"][shock_mask] >= 0.0)
    )
    full_f_error, full_f = _spectrum_agreement(_RESULTS["full_final_paths"], "full_f")
    delta_f_error, delta_f = _spectrum_agreement(
        _RESULTS["full_final_paths"], "delta_f_perturbation"
    )
    metrics = {
        "configured_ranks": _RESULTS["configured_ranks"],
        "restart_schema": _RESULTS["restart_schema"],
        "int_metadata_equal_after_restart": int_equal,
        "max_float_errors_after_restart": float_errors,
        "migration_observed": migrated,
        "source_metadata_ok": source_metadata_ok,
        "full_f_spectrum_agreement": full_f_error,
        "delta_f_spectrum_agreement": delta_f_error,
        "full_f_group_count": len(full_f["groups"]),
        "delta_f_definition": delta_f["delta_f_definition"],
    }
    logger.info("Q-016 provenance metrics: %s", metrics)
    _RESULTS["metrics"] = metrics
    return (
        _RESULTS["restart_schema"] == _EXPECTED_RESTART_SCHEMA
        and all(int_equal.values())
        and max(float_errors.values()) <= 1.0e-6
        and migrated
        and source_metadata_ok
        and full_f_error <= 1.0e-12
        and delta_f_error <= 1.0e-12
    )
