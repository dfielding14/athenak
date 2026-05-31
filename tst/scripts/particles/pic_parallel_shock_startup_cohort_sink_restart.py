"""Bounded pic_parallel_shock startup-cohort sink and restart regression."""

from __future__ import annotations

import glob
import json
import logging
import math
import os
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
        "ps_next_tag",
        "ps_mass_reservoir_global",
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
        "ps_next_tag": int(parameters["ps_next_tag"]),
        "ps_mass_reservoir_global": float(parameters["ps_mass_reservoir_global"]),
    }


def _run_athena(label, arguments, restart_path=None):
    command = ["./athena"]
    if restart_path is None:
        command += ["-i", _athena_input_path()]
    else:
        command += ["-r", os.path.relpath(restart_path, _athena_exe_dir())]
    command += list(arguments)
    logger.info("Executing %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Command failed for " + label + "\n" + output)
    return output


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
    _RESULTS["pre_metadata"] = _restart_metadata(_latest_restart(basenames["pre"]))
    _RESULTS["pre_particles"] = _particle_snapshot(basenames["pre"])

    crossed_output = _run_athena(
        "crossing_checkpoint",
        ["job/basename=" + basenames["crossed"], "time/nlim=2"],
    )
    crossed_restart = _latest_restart(basenames["crossed"])
    _RESULTS["crossed_metadata"] = _restart_metadata(crossed_restart)
    _RESULTS["crossed_particles"] = _particle_snapshot(basenames["crossed"])

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
        and all(schema == 2 for schema in restart_schemas.values())
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
