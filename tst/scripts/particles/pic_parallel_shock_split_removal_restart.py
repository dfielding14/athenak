"""Scaled Mignone startup-cohort removal timing and restart regression."""

from __future__ import annotations

import glob
import logging
import os
import subprocess
import sys

import numpy as np

_PUBLICATION_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "publication")
)
sys.path.insert(0, _PUBLICATION_DIR)
from pvtk_particles import read_particle_vtk  # noqa: E402

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_parallel_shock_split_removal_restart.athinput"
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_BIRTH_CUTOFF = 0.04
_REMOVAL_TRIGGER = 0.08
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        "ATHENA_PIC_PARALLEL_SHOCK_SPLIT_REMOVAL_EXE_DIR",
        os.path.join(os.getcwd(), "build", "src"),
    )


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


def _remove_outputs(basename):
    for dirname in ["pvtk", "rst"]:
        for path in glob.glob(
            os.path.join(_athena_exe_dir(), dirname, basename + ".*")
        ):
            if os.path.isfile(path):
                os.remove(path)


def _latest_file(dirname, pattern, label):
    matches = sorted(glob.glob(os.path.join(_athena_exe_dir(), dirname, pattern)))
    if not matches:
        raise RuntimeError("No " + label + " found for pattern " + pattern)
    return matches[-1]


def _latest_restart(basename):
    return _latest_file("rst", basename + ".*.rst", "restart")


def _restart_problem_parameters(path):
    with open(path, "rb") as handle:
        text = handle.read(65536).decode("latin1", errors="ignore")
    if "<par_end>" not in text:
        raise RuntimeError("Restart parameter header is incomplete: " + path)

    active_block = None
    parameters = {}
    for raw_line in text.split("<par_end>", 1)[0].splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("<") and line.endswith(">"):
            active_block = line[1:-1]
            continue
        if active_block == "problem" and "=" in line:
            name, value = line.split("=", 1)
            parameters[name.strip()] = value.split("#", 1)[0].strip()
    return parameters


def _strict_boolean(value):
    normalized = value.lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise RuntimeError("Invalid restart boolean " + value)


def _restart_state(path):
    parameters = _restart_problem_parameters(path)
    required = [
        "ps_remove_birth_time_before",
        "ps_remove_at_time",
        "ps_removed_excluded_early_cohort",
        "ps_removed_cr_count_global",
        "ps_removed_cr_mass_global",
        "ps_injected_cr_count_global",
    ]
    missing = [name for name in required if name not in parameters]
    if missing:
        raise RuntimeError("Restart metadata is missing " + repr(missing))
    return {
        "birth_cutoff": float(parameters["ps_remove_birth_time_before"]),
        "removal_trigger": float(parameters["ps_remove_at_time"]),
        "removed": _strict_boolean(
            parameters["ps_removed_excluded_early_cohort"]
        ),
        "removed_count": float(parameters["ps_removed_cr_count_global"]),
        "removed_mass": float(parameters["ps_removed_cr_mass_global"]),
        "injected_count": float(parameters["ps_injected_cr_count_global"]),
    }


def _particle_snapshot(basename):
    path = _latest_file(
        "pvtk", basename + ".prtcl_all.*.part.vtk", "particle VTK"
    )
    data = read_particle_vtk(path)
    required = {"ptag", "cr_source", "birth_time", "macro_weight"}
    missing = sorted(required - set(data.scalars))
    if missing:
        raise RuntimeError("Particle VTK is missing " + repr(missing))
    mask = data.scalars["cr_source"] == 1
    order = np.argsort(data.scalars["ptag"][mask])
    return {
        "tag": np.asarray(data.scalars["ptag"])[mask][order],
        "birth_time": np.asarray(data.scalars["birth_time"])[mask][order],
        "macro_weight": np.asarray(data.scalars["macro_weight"])[mask][order],
        "position": np.asarray(data.points)[mask][order],
        "velocity": np.asarray(data.vectors["vel"])[mask][order],
    }


def _run(label, basename, nlim, restart_path=None):
    _remove_outputs(basename)
    command = ["./athena"]
    if restart_path is None:
        command += ["-i", _athena_input_path()]
    else:
        command += ["-r", os.path.relpath(restart_path, _athena_exe_dir())]
    command += [
        "job/basename=" + basename,
        "time/tlim=0.12",
        "time/nlim=" + str(nlim),
    ]
    logger.info("Executing %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Command failed for " + label + "\n" + output)
    restart = _latest_restart(basename)
    return {
        "output": output,
        "restart": restart,
        "state": _restart_state(restart),
        "particles": _particle_snapshot(basename),
    }


def _max_particle_error(left, right):
    errors = []
    for name in ["birth_time", "macro_weight", "position", "velocity"]:
        if left[name].shape != right[name].shape:
            return float("inf")
        if left[name].size:
            errors.append(float(np.max(np.abs(left[name] - right[name]))))
    return max(errors, default=0.0)


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    pre = _run(
        "pre-removal checkpoint",
        "pic_parallel_shock_split_removal_pre",
        1,
    )
    full = _run(
        "uninterrupted removal",
        "pic_parallel_shock_split_removal_full",
        3,
    )
    continued = _run(
        "restart across removal trigger",
        "pic_parallel_shock_split_removal_continued",
        3,
        restart_path=pre["restart"],
    )
    _RESULTS.update({"pre": pre, "full": full, "continued": continued})


def analyze():
    logger.debug("Analyzing test " + __name__)
    pre = _RESULTS["pre"]
    full = _RESULTS["full"]
    continued = _RESULTS["continued"]

    early = pre["particles"]["birth_time"] < _BIRTH_CUTOFF
    full_late = full["particles"]["birth_time"] >= _BIRTH_CUTOFF
    states = [pre["state"], full["state"], continued["state"]]
    summary = {
        "pre_trigger_particle_count": int(pre["particles"]["tag"].size),
        "early_particle_count": int(np.count_nonzero(early)),
        "late_particle_count": int(np.count_nonzero(full_late)),
        "pre_trigger_removal_done": pre["state"]["removed"],
        "full_removal_done": full["state"]["removed"],
        "restart_removal_done": continued["state"]["removed"],
        "full_removed_count": full["state"]["removed_count"],
        "restart_removed_count": continued["state"]["removed_count"],
        "full_contains_only_late_particles": bool(np.all(full_late)),
        "restart_contains_only_late_particles": bool(
            np.all(continued["particles"]["birth_time"] >= _BIRTH_CUTOFF)
        ),
        "full_restart_tags_equal": bool(
            np.array_equal(
                full["particles"]["tag"], continued["particles"]["tag"]
            )
        ),
        "full_restart_max_particle_error": _max_particle_error(
            full["particles"], continued["particles"]
        ),
        "sink_diagnostic_counts": {
            name: case["output"].count("pic_parallel_shock removed_cr_sink:")
            for name, case in [("pre", pre), ("full", full),
                               ("continued", continued)]
        },
    }
    logger.info("Split removal/restart metrics: %s", summary)
    return (
        all(state["birth_cutoff"] == _BIRTH_CUTOFF for state in states)
        and all(state["removal_trigger"] == _REMOVAL_TRIGGER for state in states)
        and summary["early_particle_count"] > 0
        and summary["early_particle_count"] == summary["pre_trigger_particle_count"]
        and summary["late_particle_count"] > 0
        and not summary["pre_trigger_removal_done"]
        and pre["state"]["removed_count"] == 0.0
        and summary["full_removal_done"]
        and summary["restart_removal_done"]
        and summary["full_removed_count"] == summary["early_particle_count"]
        and summary["restart_removed_count"] == summary["early_particle_count"]
        and summary["full_contains_only_late_particles"]
        and summary["restart_contains_only_late_particles"]
        and summary["full_restart_tags_equal"]
        and summary["full_restart_max_particle_error"] <= 1.0e-6
        and summary["sink_diagnostic_counts"]
        == {"pre": 0, "full": 1, "continued": 1}
    )
