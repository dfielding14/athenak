"""Bounded pic_parallel_shock startup-cohort sink and restart regression."""

from __future__ import annotations

import glob
import json
import logging
import math
import os
import shutil
import struct
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

_INPUT_DECK = "tests/pic_parallel_shock_startup_cohort_sink_restart.athinput"
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_PARTICLE_INT_FIELDS = ["gid", "ptag", "species", "cr_source"]
_PARTICLE_FLOAT_FIELDS = ["points", "vel", "macro_weight", "birth_time"]
_REMOVED_FLOAT_FIELDS = [
    "ps_removed_cr_count_global",
    "ps_removed_cr_mass_global",
    "ps_removed_cr_momentum_x1_global",
    "ps_removed_cr_momentum_x2_global",
    "ps_removed_cr_momentum_x3_global",
    "ps_removed_cr_energy_global",
]
_MACRO_MASS = 1.0e-3
_PIC_RESTART_MAGIC = 0x5049435253543031
_EXPECTED_RESTART_SCHEMA = 8
_MODEL_INT_COUNT = 31
_MODEL_REAL_COUNT = 38
_REAL_BYTES = 8
_MPIEXEC = os.environ.get("MPIEXEC", "mpiexec")
_MPI_RELOAD_NPROC_ENV = "ATHENA_PIC_PARALLEL_SHOCK_MPI_RELOAD_NPROC"
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        "ATHENA_PIC_PARALLEL_SHOCK_STARTUP_COHORT_EXE_DIR",
        os.path.join(os.getcwd(), "build", "src"),
    )


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


def _remove_outputs(basename):
    exe_dir = _athena_exe_dir()
    for dirname in ["pvtk", "rst"]:
        for path in glob.glob(os.path.join(exe_dir, dirname, basename + ".*")):
            if os.path.isfile(path):
                os.remove(path)


def _latest_file(dirname, pattern, label):
    matches = sorted(glob.glob(os.path.join(_athena_exe_dir(), dirname, pattern)))
    if not matches:
        raise RuntimeError("No " + label + " files found for pattern: " + pattern)
    return matches[-1]


def _latest_restart(basename):
    return _latest_file("rst", basename + ".*.rst", "restart")


def _parse_boolean(value):
    normalized = value.strip().lower()
    if normalized in {"1", "true"}:
        return True
    if normalized in {"0", "false"}:
        return False
    raise RuntimeError("Unable to parse restart boolean: " + value)


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


def _restart_metadata(path):
    parameters = _restart_problem_parameters(path)
    required = [
        "ps_cr_ledger_schema",
        "ps_cr_ledger_complete",
        "ps_removed_excluded_early_cohort",
        *_REMOVED_FLOAT_FIELDS,
        "ps_tag_seeded",
        "ps_injection_tag_floor",
        "ps_next_tag",
        "ps_mass_reservoir_global",
        "ps_injected_cr_count_global",
    ]
    missing = [name for name in required if name not in parameters]
    if missing:
        raise RuntimeError(
            "Restart metadata is missing sink-accounting keys "
            + repr(missing)
            + ": "
            + path
        )
    return {
        "ps_cr_ledger_schema": int(parameters["ps_cr_ledger_schema"]),
        "ps_cr_ledger_complete": _parse_boolean(
            parameters["ps_cr_ledger_complete"]
        ),
        "ps_removed_excluded_early_cohort": _parse_boolean(
            parameters["ps_removed_excluded_early_cohort"]
        ),
        **{name: float(parameters[name]) for name in _REMOVED_FLOAT_FIELDS},
        "ps_tag_seeded": _parse_boolean(parameters["ps_tag_seeded"]),
        "ps_injection_tag_floor": int(parameters["ps_injection_tag_floor"]),
        "ps_next_tag": int(parameters["ps_next_tag"]),
        "ps_mass_reservoir_global": float(parameters["ps_mass_reservoir_global"]),
        "ps_injected_cr_count_global": float(
            parameters["ps_injected_cr_count_global"]
        ),
    }


def _run_athena(label, arguments, restart_path=None, nproc=1):
    command = ["./athena"]
    if restart_path is None:
        command += ["-i", _athena_input_path()]
    else:
        command += ["-r", os.path.relpath(restart_path, _athena_exe_dir())]
    command += list(arguments)
    if nproc > 1:
        command = [_MPIEXEC, "-n", str(nproc)] + command
    logger.info("Executing %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Command failed for " + label + "\n" + output)
    return output


def _run_athena_expect_fail(label, arguments, expected, restart_path=None):
    command = ["./athena"]
    if restart_path is None:
        command += ["-i", _athena_input_path()]
    else:
        command += ["-r", os.path.relpath(restart_path, _athena_exe_dir())]
    command += list(arguments)
    logger.info("Executing expected rejection %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode == 0:
        raise RuntimeError("Expected rejection passed for " + label)
    if expected not in output:
        raise RuntimeError(
            "Unexpected rejection reason for " + label + "\nExpected substring: "
            + expected + "\nOutput:\n" + output
        )


def _fnv1a64(path):
    value = 14695981039346656037
    with open(path, "rb") as handle:
        while True:
            payload = handle.read(1024 * 1024)
            if not payload:
                break
            for byte in payload:
                value ^= byte
                value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return value


def _write_completion_marker(path):
    with open(path + ".complete", "w", encoding="ascii") as handle:
        handle.write("ATHENAK_RESTART_COMPLETE_V1\n")
        handle.write("size=" + str(os.path.getsize(path)) + "\n")
        handle.write("fnv1a64=" + format(_fnv1a64(path), "016x") + "\n")


def _write_shared_restart_publication(path):
    _write_completion_marker(path)
    manifest_path = path + ".manifest"
    manifest = {
        "schema": "ATHENAK_RESTART_MANIFEST_V1",
        "members": [{
            "path": os.path.relpath(path, _athena_exe_dir()),
            "size": os.path.getsize(path),
            "fnv1a64": format(_fnv1a64(path), "016x"),
        }],
    }
    with open(manifest_path, "w", encoding="ascii") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    _write_completion_marker(manifest_path)


def _particle_restart_payload_offsets(data, restart_path):
    marker = struct.pack("<Q", _PIC_RESTART_MAGIC)
    section = data.find(marker)
    if section < 0:
        raise RuntimeError("Particle restart marker not found in " + restart_path)
    offset = section + struct.calcsize("<Q")
    meta_fmt = "<15i"
    (version, nmb_section, nrdata, nidata, *_unused) = struct.unpack_from(
        meta_fmt, data, offset
    )
    offset += struct.calcsize(meta_fmt)
    offset += _REAL_BYTES
    offset += _MODEL_INT_COUNT * struct.calcsize("<i")
    offset += _MODEL_REAL_COUNT * _REAL_BYTES
    npart_section = struct.unpack_from("<Q", data, offset)[0]
    offset += struct.calcsize("<Q")
    if (
        version != _EXPECTED_RESTART_SCHEMA
        or nmb_section <= 0
        or nrdata <= 0
        or nidata <= 3
        or npart_section <= 1
    ):
        raise RuntimeError("Unexpected particle restart metadata in " + restart_path)
    pr_real_offset = offset + nmb_section * struct.calcsize("<i")
    pr_int_offset = pr_real_offset + npart_section * nrdata * _REAL_BYTES
    return pr_real_offset, pr_int_offset, npart_section, nrdata, nidata


def _first_shock_particle_index(data, pr_int_offset, npart_section, nidata):
    for particle in range(npart_section):
        source = struct.unpack_from(
            "<i", data, pr_int_offset + (particle * nidata + 3) * 4
        )[0]
        if source == 1:
            return particle
    raise RuntimeError("No shock-injected particle found in restart payload")


def _crossed_particle_payload_mutations(tag_floor):
    return {
        "unknown_source": (
            lambda data, _pr, pi, _np, _nr, _ni: struct.pack_into(
                "<i", data, pi + 3 * 4, 99
            ),
            "pic_parallel_shock particle provenance payload is invalid",
        ),
        "baseline_tag_at_floor": (
            lambda data, _pr, pi, _np, _nr, _ni: struct.pack_into(
                "<i", data, pi + 1 * 4, tag_floor
            ),
            "pic_parallel_shock particle provenance payload is invalid",
        ),
        "non_finite_real_payload": (
            lambda data, pr, _pi, _np, _nr, _ni: struct.pack_into(
                "<d", data, pr, math.nan
            ),
            "pic_parallel_shock particle provenance payload is invalid",
        ),
        "duplicate_tag": (
            lambda data, _pr, pi, _np, _nr, nidata: struct.pack_into(
                "<i", data, pi + (nidata + 1) * 4,
                struct.unpack_from("<i", data, pi + 1 * 4)[0],
            ),
            "pic_parallel_shock particle tags are not globally unique",
        ),
    }


def _mutate_first_shock_real(field, value):
    def mutate(data, pr_real_offset, pr_int_offset, npart_section, nrdata, nidata):
        particle = _first_shock_particle_index(
            data, pr_int_offset, npart_section, nidata
        )
        struct.pack_into(
            "<d", data, pr_real_offset + (particle * nrdata + field) * _REAL_BYTES,
            value,
        )
    return mutate


def _mutate_first_shock_int(field, value):
    def mutate(data, _pr_real_offset, pr_int_offset, npart_section, _nrdata, nidata):
        particle = _first_shock_particle_index(
            data, pr_int_offset, npart_section, nidata
        )
        struct.pack_into(
            "<i", data, pr_int_offset + (particle * nidata + field) * 4, value
        )
    return mutate


def _shock_particle_payload_mutations():
    return {
        "shock_species": (
            _mutate_first_shock_int(2, 1),
            "pic_parallel_shock particle provenance payload is invalid",
        ),
        "shock_q_over_m": (
            _mutate_first_shock_real(6, 2.0),
            "pic_parallel_shock particle provenance payload is invalid",
        ),
        "shock_macro_weight": (
            _mutate_first_shock_real(22, 2.0),
            "pic_parallel_shock particle provenance payload is invalid",
        ),
        "shock_birth_time": (
            _mutate_first_shock_real(25, 1.0),
            "pic_parallel_shock particle provenance payload is invalid",
        ),
    }


def _run_corrupted_particle_restart_rejections(
    source_restart, restart_basename, mutations, continuation_nlim=3
):
    rejected = []
    for label, (mutate, expected) in mutations.items():
        corrupt_path = os.path.join(
            _athena_exe_dir(), "rst", restart_basename + "_" + label + ".00000.rst"
        )
        shutil.copyfile(source_restart, corrupt_path)
        with open(corrupt_path, "rb") as handle:
            data = bytearray(handle.read())
        pr_real_offset, pr_int_offset, npart_section, nrdata, nidata = (
            _particle_restart_payload_offsets(data, corrupt_path)
        )
        mutate(data, pr_real_offset, pr_int_offset, npart_section, nrdata, nidata)
        with open(corrupt_path, "wb") as handle:
            handle.write(data)
        _write_shared_restart_publication(corrupt_path)
        _run_athena_expect_fail(
            "corrupt_particle_" + label,
            [
                "job/basename=" + restart_basename + "_" + label + "_run",
                "time/nlim=" + str(continuation_nlim),
            ],
            expected,
            restart_path=corrupt_path,
        )
        rejected.append(label)
    return sorted(rejected)


def _particle_restart_payload_census(path):
    with open(path, "rb") as handle:
        data = bytearray(handle.read())
    _, pr_int_offset, npart_section, _, nidata = _particle_restart_payload_offsets(
        data, path
    )
    sources = [
        struct.unpack_from("<i", data, pr_int_offset + (particle * nidata + 3) * 4)[0]
        for particle in range(npart_section)
    ]
    gids = [
        struct.unpack_from("<i", data, pr_int_offset + particle * nidata * 4)[0]
        for particle in range(npart_section)
    ]
    return {
        "particle_count": int(npart_section),
        "shock_particle_count": sources.count(1),
        "distinct_meshblock_gids": len(set(gids)),
    }


def _particle_snapshot(basename):
    path = _latest_file(
        "pvtk", basename + ".prtcl_all.*.part.vtk", "particle VTK"
    )
    data = read_particle_vtk(path)
    for name in _PARTICLE_INT_FIELDS + _PARTICLE_FLOAT_FIELDS[2:]:
        if name not in data.scalars:
            raise RuntimeError("Missing particle VTK scalar: " + name)
    if "vel" not in data.vectors:
        raise RuntimeError("Missing particle VTK vector: vel")
    order = np.argsort(data.scalars["ptag"])
    return {
        "points": data.points[order],
        "vel": data.vectors["vel"][order],
        **{
            name: data.scalars[name][order]
            for name in _PARTICLE_INT_FIELDS + _PARTICLE_FLOAT_FIELDS[2:]
        },
    }


def _sink_from_pre_crossing_particles(snapshot):
    mask = snapshot["cr_source"] == 1
    weights = _MACRO_MASS * snapshot["macro_weight"][mask]
    velocity = snapshot["vel"][mask]
    momentum = np.sum(weights[:, np.newaxis] * velocity, axis=0)
    energy = np.sum(0.5 * weights * np.sum(velocity * velocity, axis=1))
    return {
        "ps_removed_cr_count_global": float(np.count_nonzero(mask)),
        "ps_removed_cr_mass_global": float(np.sum(weights)),
        "ps_removed_cr_momentum_x1_global": float(momentum[0]),
        "ps_removed_cr_momentum_x2_global": float(momentum[1]),
        "ps_removed_cr_momentum_x3_global": float(momentum[2]),
        "ps_removed_cr_energy_global": float(energy),
    }


def _max_abs(actual, expected):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    if actual.shape != expected.shape:
        return math.inf
    if actual.size == 0:
        return 0.0
    return float(np.max(np.abs(actual - expected)))


def _metadata_errors(actual, expected):
    return {
        name: abs(actual[name] - expected[name])
        for name in _REMOVED_FLOAT_FIELDS + ["ps_mass_reservoir_global"]
    }


def _run_optional_nonempty_mpi_restart_roundtrip():
    nproc = int(os.environ.get(_MPI_RELOAD_NPROC_ENV, "0"))
    if nproc <= 1:
        return "not_requested"
    segment = "pic_parallel_shock_startup_sink_mpi_seg"
    restart = "pic_parallel_shock_startup_sink_mpi_restart"
    _remove_outputs(segment)
    _remove_outputs(restart)
    mpi_layout = ["meshblock/nx1=4"]
    _run_athena(
        "nonempty_mpi_checkpoint",
        ["job/basename=" + segment, "time/nlim=1", *mpi_layout],
        nproc=nproc,
    )
    checkpoint = _latest_restart(segment)
    census = _particle_restart_payload_census(checkpoint)
    if (
        census["particle_count"] <= 0
        or census["shock_particle_count"] <= 0
        or census["distinct_meshblock_gids"] <= 1
    ):
        raise RuntimeError(
            "MPI restart checkpoint did not exercise distributed particles"
        )
    _run_athena(
        "nonempty_mpi_restart_continuation",
        ["job/basename=" + restart, "time/nlim=2"],
        restart_path=checkpoint,
        nproc=nproc,
    )
    return census


def _summary():
    pre = _RESULTS["pre_particles"]
    crossed = _RESULTS["crossed_particles"]
    full = _RESULTS["full_particles"]
    restarted = _RESULTS["restart_particles"]
    pre_metadata = _RESULTS["pre_metadata"]
    crossed_metadata = _RESULTS["crossed_metadata"]
    full_metadata = _RESULTS["full_metadata"]
    restart_metadata = _RESULTS["restart_metadata"]
    expected_sink = _sink_from_pre_crossing_particles(pre)
    sink_errors = {
        name: abs(crossed_metadata[name] - expected_sink[name])
        for name in _REMOVED_FLOAT_FIELDS
    }
    particle_int_equal = {
        name: bool(np.array_equal(restarted[name], full[name]))
        for name in _PARTICLE_INT_FIELDS
    }
    particle_float_errors = {
        name: _max_abs(restarted[name], full[name])
        for name in _PARTICLE_FLOAT_FIELDS
    }
    baseline_counts = {
        label: int(np.count_nonzero(snapshot["cr_source"] == 0))
        for label, snapshot in [
            ("pre", pre),
            ("crossed", crossed),
            ("full", full),
            ("restart", restarted),
        ]
    }
    shock_counts = {
        label: int(np.count_nonzero(snapshot["cr_source"] == 1))
        for label, snapshot in [
            ("pre", pre),
            ("crossed", crossed),
            ("full", full),
            ("restart", restarted),
        ]
    }
    return {
        "evidence_class": "bounded_local_serial_host_regression",
        "not_qualification_evidence": True,
        "cutoff_crossing_exercised": (
            not pre_metadata["ps_removed_excluded_early_cohort"]
            and crossed_metadata["ps_removed_excluded_early_cohort"]
        ),
        "sink_diagnostic_counts": {
            label: output.count("pic_parallel_shock removed_cr_sink:")
            for label, output in _RESULTS["outputs"].items()
        },
        "baseline_particle_counts": baseline_counts,
        "shock_injected_particle_counts": shock_counts,
        "expected_sink_from_pre_crossing_particles": expected_sink,
        "crossed_restart_metadata": crossed_metadata,
        "restart_schema_values": {
            label: metadata["ps_cr_ledger_schema"]
            for label, metadata in [
                ("pre", pre_metadata),
                ("crossed", crossed_metadata),
                ("full", full_metadata),
                ("restart", restart_metadata),
            ]
        },
        "restart_ledger_complete_values": {
            label: metadata["ps_cr_ledger_complete"]
            for label, metadata in [
                ("pre", pre_metadata),
                ("crossed", crossed_metadata),
                ("full", full_metadata),
                ("restart", restart_metadata),
            ]
        },
        "restart_tag_seeded_values": {
            label: metadata["ps_tag_seeded"]
            for label, metadata in [
                ("pre", pre_metadata),
                ("crossed", crossed_metadata),
                ("full", full_metadata),
                ("restart", restart_metadata),
            ]
        },
        "explicit_incomplete_restart_rejected":
            _RESULTS["explicit_incomplete_restart_rejected"],
        "invalid_restart_ledgers_rejected":
            _RESULTS["invalid_restart_ledgers_rejected"],
        "invalid_restart_particle_payloads_rejected":
            _RESULTS["invalid_restart_particle_payloads_rejected"],
        "nondecimal_ledger_roundtrip": _RESULTS["nondecimal_ledger_roundtrip"],
        "large_count_ledger_checkpoint": _RESULTS["large_count_ledger_checkpoint"],
        "optional_nonempty_mpi_restart_roundtrip":
            _RESULTS["optional_nonempty_mpi_restart_roundtrip"],
        "sink_vs_pre_crossing_particle_absolute_errors": sink_errors,
        "crossed_vs_full_metadata_absolute_errors": _metadata_errors(
            full_metadata, crossed_metadata
        ),
        "crossed_vs_restart_metadata_absolute_errors": _metadata_errors(
            restart_metadata, crossed_metadata
        ),
        "next_tag_values": {
            "crossed": crossed_metadata["ps_next_tag"],
            "full": full_metadata["ps_next_tag"],
            "restart": restart_metadata["ps_next_tag"],
        },
        "full_vs_restart_particle_int_payload_equal": particle_int_equal,
        "full_vs_restart_particle_float_payload_absolute_errors": particle_float_errors,
    }


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    basenames = {
        "pre": "pic_parallel_shock_startup_sink_pre",
        "crossed": "pic_parallel_shock_startup_sink_crossed",
        "full": "pic_parallel_shock_startup_sink_full",
        "restart": "pic_parallel_shock_startup_sink_restart",
    }
    for basename in basenames.values():
        _remove_outputs(basename)

    pre_output = _run_athena(
        "pre_crossing",
        ["job/basename=" + basenames["pre"], "time/nlim=1"],
    )
    pre_restart = _latest_restart(basenames["pre"])
    _RESULTS["pre_metadata"] = _restart_metadata(pre_restart)
    _RESULTS["pre_particles"] = _particle_snapshot(basenames["pre"])

    crossed_output = _run_athena(
        "crossing_checkpoint",
        ["job/basename=" + basenames["crossed"], "time/nlim=2"],
    )
    crossed_restart = _latest_restart(basenames["crossed"])
    _RESULTS["crossed_metadata"] = _restart_metadata(crossed_restart)
    _RESULTS["crossed_particles"] = _particle_snapshot(basenames["crossed"])
    _run_athena_expect_fail(
        "explicit_incomplete_restart",
        [
            "job/basename=" + basenames["restart"] + "_incomplete",
            "problem/ps_cr_ledger_complete=false",
        ],
        "pic_parallel_shock restart CR ledger is explicitly incomplete",
        restart_path=crossed_restart,
    )
    _RESULTS["explicit_incomplete_restart_rejected"] = True
    numeric_metadata_error = (
        "pic_parallel_shock restart CR ledger numeric metadata is invalid"
    )
    premature_sink_segment = basenames["restart"] + "_premature_sink_seg"
    _remove_outputs(premature_sink_segment)
    _run_athena(
        "premature_sink_checkpoint",
        [
            "job/basename=" + premature_sink_segment,
            "time/nlim=1",
            "problem/ps_inject_t_start=0.08",
            "problem/ps_inject_t_stop=0.16",
            "problem/ps_remove_birth_time_before=0.2",
        ],
    )
    premature_sink_restart = _latest_restart(premature_sink_segment)
    if _restart_metadata(premature_sink_restart)["ps_removed_excluded_early_cohort"]:
        raise RuntimeError("Premature-sink fixture crossed the sink threshold")
    invalid_ledger_overrides = {
        "inconsistent_injected_restart_mass": (
            ["problem/ps_injected_cr_mass_global=0.134"], numeric_metadata_error
        ),
        "inconsistent_removed_restart_mass": (
            ["problem/ps_removed_cr_mass_global=0.134"], numeric_metadata_error
        ),
        "non_finite_restart_ledger": (
            ["problem/ps_mass_reservoir_global=nan"], numeric_metadata_error
        ),
        "negative_restart_ledger": (
            ["problem/ps_removed_cr_count_global=-1"], numeric_metadata_error
        ),
        "non_integral_restart_count": (
            ["problem/ps_injected_cr_count_global=1.5"], numeric_metadata_error
        ),
        "out_of_range_restart_reservoir": (
            ["problem/ps_mass_reservoir_global=" + str(_MACRO_MASS)],
            numeric_metadata_error,
        ),
        "reset_restart_tag_floor": (
            ["problem/ps_injection_tag_floor=0"], numeric_metadata_error
        ),
        "reset_restart_tag_progression": (
            ["problem/ps_next_tag=0"], numeric_metadata_error
        ),
        "rewound_restart_tag_window": (
            [
                "time/nlim=4",
                "problem/ps_injection_tag_floor=0",
                "problem/ps_next_tag="
                + str(int(_RESULTS["crossed_metadata"]["ps_injected_cr_count_global"])),
            ],
            "pic_parallel_shock particle provenance payload is invalid",
        ),
        "unseeded_restart_tag_progression": (
            ["problem/ps_tag_seeded=false"], numeric_metadata_error
        ),
        "premature_sink_completion": (
            ["problem/ps_removed_excluded_early_cohort=true"], numeric_metadata_error
        ),
    }
    for label, (overrides, expected) in invalid_ledger_overrides.items():
        _run_athena_expect_fail(
            label,
            [
                "job/basename=" + basenames["restart"] + "_" + label,
                *overrides,
            ],
            expected,
            restart_path=(
                premature_sink_restart if label == "premature_sink_completion"
                else crossed_restart
            ),
        )
    _RESULTS["invalid_restart_ledgers_rejected"] = sorted(
        invalid_ledger_overrides
    )
    _RESULTS["invalid_restart_particle_payloads_rejected"] = sorted(
        _run_corrupted_particle_restart_rejections(
            crossed_restart,
            basenames["restart"] + "_corrupt",
            _crossed_particle_payload_mutations(
                _RESULTS["crossed_metadata"]["ps_injection_tag_floor"]
            ),
        )
        + _run_corrupted_particle_restart_rejections(
            pre_restart,
            basenames["restart"] + "_corrupt_pre",
            _shock_particle_payload_mutations(),
        )
        + _run_corrupted_particle_restart_rejections(
            pre_restart,
            basenames["restart"] + "_corrupt_zero_cycle",
            {
                "zero_cycle_shock_macro_weight": (
                    _mutate_first_shock_real(22, 2.0),
                    "pic_parallel_shock particle provenance payload is invalid",
                )
            },
            continuation_nlim=1,
        )
    )

    post_sink_segment = basenames["restart"] + "_post_sink_late_seg"
    _remove_outputs(post_sink_segment)
    _run_athena(
        "post_sink_late_particle_checkpoint",
        [
            "job/basename=" + post_sink_segment,
            "time/nlim=3",
            "problem/ps_inject_t_stop=1.0",
        ],
    )
    post_sink_restart = _latest_restart(post_sink_segment)
    if not _restart_metadata(post_sink_restart)["ps_removed_excluded_early_cohort"]:
        raise RuntimeError("Post-sink lifecycle fixture did not cross the sink threshold")
    _RESULTS["invalid_restart_particle_payloads_rejected"] += (
        _run_corrupted_particle_restart_rejections(
            post_sink_restart,
            basenames["restart"] + "_corrupt_post_sink",
            {
                "post_sink_early_birth_time": (
                    _mutate_first_shock_real(25, 0.0),
                    "pic_parallel_shock particle provenance payload is invalid",
                )
            },
            continuation_nlim=4,
        )
    )
    _RESULTS["invalid_restart_particle_payloads_rejected"].sort()

    precision_seg = basenames["restart"] + "_precision_seg"
    precision_rst = basenames["restart"] + "_precision_rst"
    _remove_outputs(precision_seg)
    _remove_outputs(precision_rst)
    _run_athena(
        "nondecimal_ledger_checkpoint",
        [
            "job/basename=" + precision_seg,
            "time/nlim=1",
            "particles/deposit_qscale=0.001234567890123",
        ],
    )
    _run_athena(
        "nondecimal_ledger_continuation",
        ["job/basename=" + precision_rst, "time/nlim=2"],
        restart_path=_latest_restart(precision_seg),
    )
    _RESULTS["nondecimal_ledger_roundtrip"] = True
    large_count_segment = basenames["restart"] + "_large_count_seg"
    _remove_outputs(large_count_segment)
    _run_athena(
        "large_count_ledger_checkpoint",
        [
            "job/basename=" + large_count_segment,
            "time/nlim=1",
            "problem/ps_eta=100.0",
        ],
    )
    large_count_restart = _latest_restart(large_count_segment)
    large_count_metadata = _restart_metadata(large_count_restart)
    if large_count_metadata["ps_injected_cr_count_global"] < 100000:
        raise RuntimeError("Large-count ledger fixture did not inject enough particles")
    _run_athena_expect_fail(
        "large_count_inconsistent_mass",
        [
            "job/basename=" + large_count_segment + "_bad_mass",
            "problem/ps_injected_cr_mass_global="
            + str(
                large_count_metadata["ps_injected_cr_count_global"] * _MACRO_MASS
                - _MACRO_MASS
            ),
        ],
        numeric_metadata_error,
        restart_path=large_count_restart,
    )
    large_count_restart_basename = large_count_segment + "_rst"
    _remove_outputs(large_count_restart_basename)
    _run_athena(
        "large_count_ledger_continuation",
        ["job/basename=" + large_count_restart_basename, "time/nlim=2"],
        restart_path=large_count_restart,
    )
    _RESULTS["large_count_ledger_checkpoint"] = {
        "injected_particle_count": large_count_metadata["ps_injected_cr_count_global"],
        "injected_mass": (
            large_count_metadata["ps_injected_cr_count_global"] * _MACRO_MASS
        ),
        "inconsistent_mass_rejected": True,
        "restart_continuation_passed": True,
    }
    _RESULTS["optional_nonempty_mpi_restart_roundtrip"] = (
        _run_optional_nonempty_mpi_restart_roundtrip()
    )

    full_output = _run_athena(
        "uninterrupted",
        ["job/basename=" + basenames["full"], "time/nlim=3"],
    )
    _RESULTS["full_metadata"] = _restart_metadata(_latest_restart(basenames["full"]))
    _RESULTS["full_particles"] = _particle_snapshot(basenames["full"])

    restart_output = _run_athena(
        "restart_continuation",
        [
            "job/basename=" + basenames["restart"],
            "time/nlim=3",
            "output1/file_number=0",
            "output2/file_number=0",
        ],
        restart_path=crossed_restart,
    )
    _RESULTS["restart_metadata"] = _restart_metadata(
        _latest_restart(basenames["restart"])
    )
    _RESULTS["restart_particles"] = _particle_snapshot(basenames["restart"])
    _RESULTS["outputs"] = {
        "pre": pre_output,
        "crossed": crossed_output,
        "full": full_output,
        "restart": restart_output,
    }


def analyze():
    logger.debug("Analyzing test " + __name__)
    summary = _summary()
    logger.info("PIC parallel-shock startup-cohort sink metrics: %s", summary)
    pre_metadata = _RESULTS["pre_metadata"]
    crossed_metadata = _RESULTS["crossed_metadata"]
    diagnostics = summary["sink_diagnostic_counts"]
    baseline_counts = summary["baseline_particle_counts"]
    shock_counts = summary["shock_injected_particle_counts"]
    sink_errors = summary["sink_vs_pre_crossing_particle_absolute_errors"]
    crossed_vs_full = summary["crossed_vs_full_metadata_absolute_errors"]
    crossed_vs_restart = summary["crossed_vs_restart_metadata_absolute_errors"]
    next_tags = summary["next_tag_values"]
    restart_schemas = summary["restart_schema_values"]
    restart_ledgers_complete = summary["restart_ledger_complete_values"]
    restart_tags_seeded = summary["restart_tag_seeded_values"]
    expected_sink = summary["expected_sink_from_pre_crossing_particles"]
    particle_errors = summary["full_vs_restart_particle_float_payload_absolute_errors"]
    momentum_scale = max(
        1.0,
        abs(expected_sink["ps_removed_cr_momentum_x1_global"]),
        abs(expected_sink["ps_removed_cr_momentum_x2_global"]),
        abs(expected_sink["ps_removed_cr_momentum_x3_global"]),
    )
    return (
        summary["cutoff_crossing_exercised"]
        and summary["explicit_incomplete_restart_rejected"]
        and summary["invalid_restart_ledgers_rejected"] == [
            "inconsistent_injected_restart_mass",
            "inconsistent_removed_restart_mass",
            "negative_restart_ledger",
            "non_finite_restart_ledger",
            "non_integral_restart_count",
            "out_of_range_restart_reservoir",
            "premature_sink_completion",
            "reset_restart_tag_floor",
            "reset_restart_tag_progression",
            "rewound_restart_tag_window",
            "unseeded_restart_tag_progression",
        ]
        and summary["invalid_restart_particle_payloads_rejected"] == [
            "baseline_tag_at_floor",
            "duplicate_tag",
            "non_finite_real_payload",
            "post_sink_early_birth_time",
            "shock_birth_time",
            "shock_macro_weight",
            "shock_q_over_m",
            "shock_species",
            "unknown_source",
            "zero_cycle_shock_macro_weight",
        ]
        and summary["nondecimal_ledger_roundtrip"]
        and summary["large_count_ledger_checkpoint"]["injected_particle_count"] >= 100000
        and summary["large_count_ledger_checkpoint"]["inconsistent_mass_rejected"]
        and summary["large_count_ledger_checkpoint"]["restart_continuation_passed"]
        and (
            summary["optional_nonempty_mpi_restart_roundtrip"] == "not_requested"
            or (
                summary["optional_nonempty_mpi_restart_roundtrip"]["particle_count"] > 0
                and summary["optional_nonempty_mpi_restart_roundtrip"][
                    "shock_particle_count"
                ] > 0
                and summary["optional_nonempty_mpi_restart_roundtrip"][
                    "distinct_meshblock_gids"
                ] > 1
            )
        )
        and all(schema == 3 for schema in restart_schemas.values())
        and all(restart_ledgers_complete.values())
        and all(restart_tags_seeded.values())
        and diagnostics == {"pre": 0, "crossed": 1, "full": 1, "restart": 0}
        and all(pre_metadata[name] == 0.0 for name in _REMOVED_FLOAT_FIELDS)
        and crossed_metadata["ps_removed_cr_count_global"] > 0.0
        and crossed_metadata["ps_removed_cr_energy_global"] > 0.0
        and all(math.isfinite(crossed_metadata[name]) for name in _REMOVED_FLOAT_FIELDS)
        and shock_counts["pre"] > 0
        and shock_counts["crossed"] == shock_counts["full"] == shock_counts["restart"] == 0
        and len(set(baseline_counts.values())) == 1
        and baseline_counts["pre"] > 0
        and sink_errors["ps_removed_cr_count_global"] == 0.0
        and sink_errors["ps_removed_cr_mass_global"] <= 1.0e-10
        and max(
            sink_errors["ps_removed_cr_momentum_x1_global"],
            sink_errors["ps_removed_cr_momentum_x2_global"],
            sink_errors["ps_removed_cr_momentum_x3_global"],
        )
        <= 5.0e-6 * momentum_scale
        and sink_errors["ps_removed_cr_energy_global"]
        <= 5.0e-6 * expected_sink["ps_removed_cr_energy_global"]
        and max(crossed_vs_full.values()) <= 1.0e-12
        and max(crossed_vs_restart.values()) <= 1.0e-12
        and next_tags["crossed"] == next_tags["full"] == next_tags["restart"]
        and all(summary["full_vs_restart_particle_int_payload_equal"].values())
        and max(particle_errors.values()) <= 1.0e-6
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
