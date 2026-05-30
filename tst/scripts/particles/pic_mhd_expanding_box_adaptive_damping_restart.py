"""Bounded Q-008/Q-033 active-MHD expanding-box restart continuity regression."""

from __future__ import annotations

import glob
import json
import logging
import os
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

_INPUT_DECK = "tests/pic_mhd_expanding_box_adaptive_damping_smoke.athinput"
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_MHD_FIELDS = {
    "mhd_u": ["dens", "mom1", "mom2", "mom3", "ener"],
    "mhd_w": ["dens", "velx", "vely", "velz", "eint"],
    "mhd_bcc": ["bcc1", "bcc2", "bcc3"],
}
_PARTICLE_INT_FIELDS = ["gid", "ptag", "species", "cr_source"]
_PARTICLE_FLOAT_FIELDS = [
    "points",
    "vel",
    "macro_weight",
    "birth_time",
    "deltaf_f0",
    "deltaf_weight",
]
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), "build", "src")


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


def _remove_outputs(basename):
    exe_dir = _athena_exe_dir()
    history = os.path.join(exe_dir, basename + ".mhd.hst")
    if os.path.exists(history):
        os.remove(history)
    for dirname in ["bin", "pvtk", "rst"]:
        for path in glob.glob(os.path.join(exe_dir, dirname, basename + ".*")):
            if os.path.isfile(path):
                os.remove(path)


def _latest_file(dirname, pattern, label):
    matches = sorted(glob.glob(os.path.join(_athena_exe_dir(), dirname, pattern)))
    if not matches:
        raise RuntimeError("No " + label + " files found for pattern: " + pattern)
    return matches[-1]


def _run_athena(label, arguments, restart_file=None):
    command = ["./athena"]
    if restart_file is None:
        command += ["-i", _athena_input_path()]
    else:
        command += ["-r", restart_file]
    command += list(arguments)
    logger.info("Executing %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Command failed for " + label + "\n" + output)
    for token in [
        "physical_mode=extended_mhd_pic",
        "background=coupled",
        "feedback=coupled",
        "deltaf_adapt=global_bikappa_moments_experimental",
        "expanding_box=on",
        "wave_damping=ion_neutral_friction",
    ]:
        if token not in output:
            raise RuntimeError(label + " missing runtime token: " + token)
    return output


def _mhd_snapshot(basename):
    snapshot = {}
    times = []
    for file_id, fields in _MHD_FIELDS.items():
        path = _latest_file("bin", basename + "." + file_id + ".*.bin", file_id)
        data = bin_convert.read_binary_as_athdf(path)
        times.append(float(data["Time"]))
        snapshot[file_id] = {
            field: np.asarray(data[field], dtype=np.float64) for field in fields
        }
    if max(times) != min(times):
        raise RuntimeError("MHD output times disagree for " + basename)
    return {"time": times[0], "fields": snapshot}


def _history_endpoint(basename):
    path = os.path.join(_athena_exe_dir(), basename + ".mhd.hst")
    rows = np.loadtxt(path)
    if rows.ndim == 1:
        rows = rows[np.newaxis, :]
    if rows.shape[0] == 0:
        raise RuntimeError("No history rows found for " + basename)
    return np.asarray(rows[-1], dtype=np.float64)


def _particle_snapshot(basename):
    path = _latest_file(
        "pvtk", basename + ".prtcl_all.*.part.vtk", "particle VTK"
    )
    data = read_particle_vtk(path)
    for field in _PARTICLE_INT_FIELDS + _PARTICLE_FLOAT_FIELDS[2:]:
        if field not in data.scalars:
            raise RuntimeError("Missing particle VTK scalar: " + field)
    if "vel" not in data.vectors:
        raise RuntimeError("Missing particle VTK vector: vel")
    order = np.argsort(data.scalars["ptag"])
    return {
        "points": data.points[order],
        "vel": data.vectors["vel"][order],
        **{
            field: data.scalars[field][order]
            for field in _PARTICLE_INT_FIELDS + _PARTICLE_FLOAT_FIELDS[2:]
        },
    }


def _max_abs(actual, expected):
    return float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))


def _summary():
    full_mhd = _RESULTS["full_mhd"]
    restart_mhd = _RESULTS["restart_mhd"]
    full_particles = _RESULTS["full_particles"]
    restart_particles = _RESULTS["restart_particles"]
    mhd_errors = {
        file_id + ":" + field: _max_abs(
            restart_mhd["fields"][file_id][field],
            full_mhd["fields"][file_id][field],
        )
        for file_id, fields in _MHD_FIELDS.items()
        for field in fields
    }
    particle_int_equal = {
        field: bool(np.array_equal(restart_particles[field], full_particles[field]))
        for field in _PARTICLE_INT_FIELDS
    }
    particle_float_errors = {
        field: _max_abs(restart_particles[field], full_particles[field])
        for field in _PARTICLE_FLOAT_FIELDS
    }
    return {
        "evidence_class": "bounded_local_serial_host_regression",
        "not_frontier_qualification_evidence": True,
        "checkpoint_continuations": 1,
        "full_time": full_mhd["time"],
        "restart_time": restart_mhd["time"],
        "mhd_time_absolute_error": abs(restart_mhd["time"] - full_mhd["time"]),
        "mhd_field_absolute_errors": mhd_errors,
        "history_endpoint_absolute_error": _max_abs(
            _RESULTS["restart_history"], _RESULTS["full_history"]
        ),
        "particle_count": int(full_particles["ptag"].size),
        "particle_int_payload_equal": particle_int_equal,
        "particle_float_payload_absolute_errors": particle_float_errors,
        "restart_refit_observed": _RESULTS["restart_refit_observed"],
    }


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    full = "pic_mhd_box_adaptive_damping_restart_full"
    segment = "pic_mhd_box_adaptive_damping_restart_segment"
    restarted = "pic_mhd_box_adaptive_damping_restart_continued"
    for basename in [full, segment, restarted]:
        _remove_outputs(basename)

    _run_athena(
        "full",
        ["job/basename=" + full, "time/nlim=2", "output6/dcycle=0"],
    )
    _RESULTS["full_mhd"] = _mhd_snapshot(full)
    _RESULTS["full_history"] = _history_endpoint(full)
    _RESULTS["full_particles"] = _particle_snapshot(full)

    _run_athena("segment", ["job/basename=" + segment, "time/nlim=1"])
    checkpoint = _latest_file("rst", segment + ".*.rst", "restart")
    restart_output = _run_athena(
        "restart",
        [
            "job/basename=" + restarted,
            "time/nlim=2",
            "output6/dcycle=0",
        ],
        restart_file=os.path.relpath(checkpoint, _athena_exe_dir()),
    )
    _RESULTS["restart_mhd"] = _mhd_snapshot(restarted)
    _RESULTS["restart_history"] = _history_endpoint(restarted)
    _RESULTS["restart_particles"] = _particle_snapshot(restarted)
    _RESULTS["restart_refit_observed"] = "PIC adaptive delta-f fit:" in restart_output


def analyze():
    logger.debug("Analyzing test " + __name__)
    summary = _summary()
    logger.info("Q-008/Q-033 expanding-box restart metrics: %s", summary)
    particle_errors = summary["particle_float_payload_absolute_errors"]
    return (
        summary["checkpoint_continuations"] == 1
        and summary["full_time"] > 0.0
        and summary["mhd_time_absolute_error"] <= 1.0e-14
        and max(summary["mhd_field_absolute_errors"].values()) <= 1.0e-12
        and summary["history_endpoint_absolute_error"] <= 1.0e-10
        and summary["particle_count"] == 64
        and all(summary["particle_int_payload_equal"].values())
        and max(particle_errors.values()) <= 1.0e-6
        and not summary["restart_refit_observed"]
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
