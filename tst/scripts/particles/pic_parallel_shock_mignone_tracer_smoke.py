"""Tiny 2D smoke for Mignone shock-tracer swept-mass injection."""

from __future__ import annotations

import glob
import logging
import math
import os
import re
import subprocess
import sys

import numpy as np

_PUBLICATION_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "publication")
)
sys.path.insert(0, _PUBLICATION_DIR)
from pvtk_particles import read_particle_vtk  # noqa: E402

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_parallel_shock_mignone_tracer_smoke.athinput"
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_BASENAME = "pic_parallel_shock_mignone_tracer_smoke"
_PRE_BASENAME = _BASENAME + "_pre"
_CONTINUED_BASENAME = _BASENAME + "_continued"
_SINGLE_BLOCK_BASENAME = _BASENAME + "_single_block"
_ALL_DEFER_BASENAME = _BASENAME + "_all_defer"
_DIAG_PREFIX = "pic_parallel_shock mignone_injection_diag:"
_MACRO_MASS = 0.05
_DX1 = 10.0
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        "ATHENA_PIC_PARALLEL_SHOCK_MIGNONE_SMOKE_EXE_DIR",
        os.path.join(os.getcwd(), "build", "src"),
    )


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


def _remove_outputs():
    for dirname in ["pvtk", "rst"]:
        basenames = [
            _BASENAME,
            _PRE_BASENAME,
            _CONTINUED_BASENAME,
            _SINGLE_BLOCK_BASENAME,
            _ALL_DEFER_BASENAME,
        ]
        for basename in basenames:
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


def _restart_state(basename):
    path = _latest_file("rst", basename + ".*.rst", "restart")
    parameters = _restart_problem_parameters(path)
    required = [
        "ps_injection_mode",
        "ps_mass_reservoir_global",
        "ps_injected_cr_count_global",
        "ps_injected_cr_mass_global",
        "ps_injected_cr_momentum_x1_global",
        "ps_injected_cr_momentum_x2_global",
        "ps_injected_cr_momentum_x3_global",
        "ps_injected_cr_energy_global",
        "ps_next_tag",
        "ps_mignone_accumulation_start_time",
        "ps_mignone_last_injection_time",
        "ps_mignone_last_accumulation_dt",
        "ps_mignone_last_swept_mass",
        "ps_mignone_injection_events",
    ]
    missing = [name for name in required if name not in parameters]
    if missing:
        raise RuntimeError("Restart metadata is missing " + repr(missing))
    return {
        "injection_mode": parameters["ps_injection_mode"],
        "mass_reservoir": float(parameters["ps_mass_reservoir_global"]),
        "injected_count": float(parameters["ps_injected_cr_count_global"]),
        "injected_mass": float(parameters["ps_injected_cr_mass_global"]),
        "injected_momentum": tuple(
            float(parameters[name])
            for name in (
                "ps_injected_cr_momentum_x1_global",
                "ps_injected_cr_momentum_x2_global",
                "ps_injected_cr_momentum_x3_global",
            )
        ),
        "injected_energy": float(parameters["ps_injected_cr_energy_global"]),
        "next_tag": int(parameters["ps_next_tag"]),
        "accumulation_start_time": float(
            parameters["ps_mignone_accumulation_start_time"]
        ),
        "last_injection_time": float(
            parameters["ps_mignone_last_injection_time"]
        ),
        "last_accumulation_dt": float(
            parameters["ps_mignone_last_accumulation_dt"]
        ),
        "last_swept_mass": float(parameters["ps_mignone_last_swept_mass"]),
        "injection_events": int(parameters["ps_mignone_injection_events"]),
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
    snapshot = {
        name: np.asarray(data.scalars[name])
        for name in sorted(required)
    }
    snapshot["position"] = np.asarray(data.points)
    snapshot["velocity"] = np.asarray(data.vectors["vel"])
    return snapshot


def _parse_diagnostics(output):
    integer_fields = {
        "event",
        "cycle",
        "requested_count",
        "injected_count",
        "deferred_cells",
    }
    required = {
        "event",
        "cycle",
        "time",
        "accumulation_dt",
        "tracer_mass",
        "swept_mass",
        "swept_fraction_trigger",
        "hot_swept_mass",
        "rejected_swept_fraction",
        "hot_pressure_min",
        "shock_x",
        "injected_x_min",
        "injected_x_max",
        "downstream_fraction",
        "requested_count",
        "injected_count",
        "deferred_cells",
        "deferred_swept_mass",
    }
    diagnostics = []
    for line in output.splitlines():
        if _DIAG_PREFIX not in line:
            continue
        payload = line.split(_DIAG_PREFIX, 1)[1]
        raw = dict(re.findall(r"([a-z_]+)=([^\s]+)", payload))
        missing = sorted(required - set(raw))
        if missing:
            raise RuntimeError("Mignone diagnostic is missing " + repr(missing))
        event = {
            name: int(raw[name]) if name in integer_fields else float(raw[name])
            for name in required
        }
        diagnostics.append(event)
    if not diagnostics:
        raise RuntimeError("No Mignone injection diagnostics were emitted")
    return diagnostics


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    _remove_outputs()
    command = ["./athena", "-i", _athena_input_path()]
    logger.info("Executing Mignone tracer smoke: %s", " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Mignone tracer smoke failed\n" + output)
    continuous = {
        "diagnostics": _parse_diagnostics(output),
        "restart": _restart_state(_BASENAME),
        "particles": _particle_snapshot(_BASENAME),
    }
    single_block_command = command + [
        "job/basename=" + _SINGLE_BLOCK_BASENAME,
        "meshblock/nx1=128",
        "mesh_refinement/max_nmb_per_rank=1",
    ]
    single_block = subprocess.run(
        single_block_command,
        cwd=_athena_exe_dir(),
        capture_output=True,
        text=True,
    )
    single_block_output = (single_block.stdout or "") + (
        single_block.stderr or ""
    )
    if single_block.returncode != 0:
        raise RuntimeError(
            "Mignone single-block reference failed\n" + single_block_output
        )

    pre_command = command + [
        "job/basename=" + _PRE_BASENAME,
        "time/nlim=30",
    ]
    pre = subprocess.run(
        pre_command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    pre_output = (pre.stdout or "") + (pre.stderr or "")
    if pre.returncode != 0:
        raise RuntimeError("Mignone accumulation checkpoint failed\n" + pre_output)
    pre_restart = _latest_file(
        "rst", _PRE_BASENAME + ".*.rst", "accumulation restart"
    )
    continued_command = [
        "./athena",
        "-r",
        os.path.relpath(pre_restart, _athena_exe_dir()),
        "job/basename=" + _CONTINUED_BASENAME,
        "time/nlim=250",
        "time/tlim=4.3",
    ]
    continued = subprocess.run(
        continued_command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    continued_output = (continued.stdout or "") + (continued.stderr or "")
    if continued.returncode != 0:
        raise RuntimeError("Mignone accumulation restart failed\n" + continued_output)

    all_defer_command = command + [
        "job/basename=" + _ALL_DEFER_BASENAME,
        "problem/ps_p_floor_frac=1.0e12",
        "time/nlim=60",
    ]
    logger.info(
        "Executing Mignone all-defer stress: %s", " ".join(all_defer_command)
    )
    all_defer = subprocess.run(
        all_defer_command,
        cwd=_athena_exe_dir(),
        capture_output=True,
        text=True,
    )
    all_defer_output = (all_defer.stdout or "") + (all_defer.stderr or "")
    if all_defer.returncode != 0:
        raise RuntimeError("Mignone all-defer stress failed\n" + all_defer_output)
    _RESULTS.update(
        {
            "continuous": continuous,
            "single_block": {
                "diagnostics": _parse_diagnostics(single_block_output),
                "restart": _restart_state(_SINGLE_BLOCK_BASENAME),
                "particles": _particle_snapshot(_SINGLE_BLOCK_BASENAME),
            },
            "pre_diagnostics": _parse_diagnostics(pre_output),
            "continued": {
                "diagnostics": _parse_diagnostics(continued_output),
                "restart": _restart_state(_CONTINUED_BASENAME),
                "particles": _particle_snapshot(_CONTINUED_BASENAME),
            },
            "all_defer": {
                "diagnostics": _parse_diagnostics(all_defer_output),
                "restart": _restart_state(_ALL_DEFER_BASENAME),
                "output": all_defer_output,
            },
        }
    )


def analyze():
    logger.debug("Analyzing test " + __name__)
    continuous = _RESULTS["continuous"]
    diagnostics = continuous["diagnostics"]
    restart = continuous["restart"]
    particles = continuous["particles"]
    continued = _RESULTS["continued"]
    restarted_diagnostics = (
        _RESULTS["pre_diagnostics"] + continued["diagnostics"]
    )
    restart_equivalent = (
        restarted_diagnostics == diagnostics
        and continued["restart"] == restart
        and all(
            np.array_equal(particles[name], continued["particles"][name])
            for name in particles
        )
    )
    checkpoint_during_accumulation = (
        _RESULTS["pre_diagnostics"][-1]["cycle"] < 30
        and continued["diagnostics"][0]["cycle"] > 30
    )
    single_block = _RESULTS["single_block"]
    decomposition_equivalent = (
        single_block["diagnostics"] == diagnostics
        and single_block["restart"] == restart
        and all(
            np.array_equal(particles[name], single_block["particles"][name])
            for name in particles
        )
    )

    finite_fields = [
        "time",
        "accumulation_dt",
        "tracer_mass",
        "swept_mass",
        "swept_fraction_trigger",
        "hot_swept_mass",
        "rejected_swept_fraction",
        "hot_pressure_min",
        "shock_x",
        "injected_x_min",
        "injected_x_max",
        "downstream_fraction",
        "deferred_swept_mass",
    ]
    finite = all(
        math.isfinite(event[name])
        for event in diagnostics
        for name in finite_fields
    )
    positive_budget = all(
        event["tracer_mass"] > 0.0
        and event["swept_mass"] > 0.0
        and event["swept_fraction_trigger"] == 0.8
        and event["hot_pressure_min"] == 30.0
        and event["hot_swept_mass"] > 0.0
        and event["hot_swept_mass"] <= event["swept_mass"]
        and event["hot_swept_mass"]
        > event["swept_fraction_trigger"] * event["tracer_mass"]
        and event["swept_mass"] <= (1.0 + 1.0e-12) * event["tracer_mass"]
        and 0.0 <= event["rejected_swept_fraction"] <= 1.0
        and abs(
            event["rejected_swept_fraction"]
            - (1.0 - event["hot_swept_mass"] / event["swept_mass"])
        ) <= 1.0e-12
        and event["requested_count"] == event["injected_count"]
        and event["injected_count"] > 0
        and event["deferred_cells"] == 0
        and event["deferred_swept_mass"] == 0.0
        for event in diagnostics
    )
    event_ids = [event["event"] for event in diagnostics]
    cycles = [event["cycle"] for event in diagnostics]
    event_times = np.array([event["time"] for event in diagnostics])
    accumulation_dt = np.array(
        [event["accumulation_dt"] for event in diagnostics]
    )
    local_downstream = all(
        event["downstream_fraction"] > 0.8
        and event["injected_x_min"] <= event["injected_x_max"]
        and event["injected_x_min"] >= event["shock_x"] - 4.0 * _DX1
        and event["injected_x_max"] <= event["shock_x"] + 2.0 * _DX1
        for event in diagnostics
    )
    cadence_is_sane = bool(
        np.all(accumulation_dt > 0.0)
        and all(right - left > 1 for left, right in zip(cycles, cycles[1:]))
        and np.all(np.diff(event_times) > 0.0)
    )
    injected_from_events = sum(event["injected_count"] for event in diagnostics)
    reservoir = 0.0
    requested_counts = []
    for event in diagnostics:
        budget = reservoir + 2.0e-3 * event["hot_swept_mass"]
        requested = math.floor(budget / _MACRO_MASS)
        requested_counts.append(requested)
        reservoir = budget - requested * _MACRO_MASS
    hot_swept_mass_budget = requested_counts == [
        event["injected_count"] for event in diagnostics
    ]
    particle_count = int(particles["ptag"].size)
    birth_times_match_events = all(
        np.any(np.abs(particles["birth_time"] - event["time"]) <= 1.0e-5)
        for event in diagnostics
    )
    all_defer = _RESULTS["all_defer"]
    all_defer_diagnostics = all_defer["diagnostics"]
    all_defer_restart = all_defer["restart"]
    repeated_defer_retries = any(
        right["cycle"] == left["cycle"] + 1
        for left, right in zip(
            all_defer_diagnostics, all_defer_diagnostics[1:]
        )
    )
    all_defer_events_are_clean = all(
        event["event"] == 1
        and event["requested_count"] > 0
        and event["injected_count"] == 0
        and event["deferred_cells"] > 0
        and event["deferred_swept_mass"] > 0.0
        and event["deferred_swept_mass"]
        <= (1.0 + 1.0e-12) * event["hot_swept_mass"]
        for event in all_defer_diagnostics
    )
    all_defer_state_unchanged = (
        all_defer_restart["injection_mode"] == "mignone_tracer"
        and all_defer_restart["mass_reservoir"] == 0.0
        and all_defer_restart["injected_count"] == 0.0
        and all_defer_restart["injected_mass"] == 0.0
        and all(value == 0.0 for value in all_defer_restart["injected_momentum"])
        and all_defer_restart["injected_energy"] == 0.0
        and all_defer_restart["next_tag"] == 0
        and all_defer_restart["injection_events"] == 0
        and all_defer_restart["accumulation_start_time"] == 0.0
        and all_defer_restart["last_injection_time"] == -1.0
        and all_defer_restart["last_accumulation_dt"] == 0.0
        and all_defer_restart["last_swept_mass"] == 0.0
    )
    no_fatal_floor_diagnostic = (
        "gas_subtraction_floor_diag:" not in all_defer["output"]
        and "gas subtraction would violate a fluid floor"
        not in all_defer["output"]
    )
    summary = {
        "event_count": len(diagnostics),
        "event_times": event_times.tolist(),
        "accumulation_dt": accumulation_dt.tolist(),
        "injected_from_events": injected_from_events,
        "hot_swept_mass_budget": hot_swept_mass_budget,
        "active_particle_count": particle_count,
        "local_downstream": local_downstream,
        "cadence_is_sane": cadence_is_sane,
        "restart_injection_events": restart["injection_events"],
        "restart_injected_count": restart["injected_count"],
        "restart_equivalent": restart_equivalent,
        "checkpoint_during_accumulation": checkpoint_during_accumulation,
        "decomposition_equivalent": decomposition_equivalent,
        "all_defer_attempts": len(all_defer_diagnostics),
        "repeated_defer_retries": repeated_defer_retries,
        "all_defer_events_are_clean": all_defer_events_are_clean,
        "all_defer_state_unchanged": all_defer_state_unchanged,
        "no_fatal_floor_diagnostic": no_fatal_floor_diagnostic,
    }
    logger.info("Mignone tracer smoke metrics: %s", summary)
    return (
        len(diagnostics) >= 2
        and finite
        and positive_budget
        and event_ids == list(range(1, len(diagnostics) + 1))
        and all(right > left for left, right in zip(cycles, cycles[1:]))
        and cadence_is_sane
        and local_downstream
        and restart["injection_mode"] == "mignone_tracer"
        and restart["injection_events"] == len(diagnostics)
        and restart["last_injection_time"] == diagnostics[-1]["time"]
        and restart["accumulation_start_time"] == diagnostics[-1]["time"]
        and restart["last_accumulation_dt"] == diagnostics[-1]["accumulation_dt"]
        and restart["last_swept_mass"] == diagnostics[-1]["swept_mass"]
        and injected_from_events == particle_count
        and hot_swept_mass_budget
        and restart["injected_count"] == particle_count
        and abs(restart["injected_mass"] - particle_count * _MACRO_MASS)
        <= 1.0e-12
        and restart["next_tag"] == particle_count
        and np.unique(particles["ptag"]).size == particle_count
        and np.all(particles["cr_source"] == 1)
        and np.all(particles["macro_weight"] == 1.0)
        and birth_times_match_events
        and checkpoint_during_accumulation
        and restart_equivalent
        and decomposition_equivalent
        and len(all_defer_diagnostics) >= 2
        and repeated_defer_retries
        and all_defer_events_are_clean
        and all_defer_state_unchanged
        and no_fatal_floor_diagnostic
    )
