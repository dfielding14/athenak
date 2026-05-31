"""Periodic AMR scheduled-intermediate restart ordering parity regression."""

from __future__ import annotations

import glob
import json
import logging
import os
import re
import subprocess
import sys

import numpy as np

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

_INPUT_DECK = "tests/pic_q009_periodic_amr_restart_ordering_parity.athinput"
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_BASENAMES = {
    "full": "pic_q009_periodic_amr_ordering_full",
    "seed": "pic_q009_periodic_amr_ordering_seed",
    "probe": "pic_q009_periodic_amr_ordering_probe",
    "resumed": "pic_q009_periodic_amr_ordering_resumed",
}
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        "ATHENA_Q009_PERIODIC_AMR_ORDERING_EXE_DIR",
        os.path.join(os.getcwd(), "build", "src"),
    )


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


def _remove_outputs(basename):
    for dirname in ["bin", "pvtk", "rst"]:
        pattern = os.path.join(_athena_exe_dir(), dirname, basename + ".*")
        for path in glob.glob(pattern):
            if os.path.isfile(path):
                os.remove(path)


def _matching_files(pattern, label):
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise RuntimeError("No " + label + " files found for pattern: " + pattern)
    return paths


def _latest_file(dirname, pattern, label):
    return _matching_files(os.path.join(_athena_exe_dir(), dirname, pattern), label)[-1]


def _restart_path(basename, number):
    return os.path.join(
        _athena_exe_dir(), "rst", f"{basename}.{number:05d}.rst"
    )


def _parse_amr_telemetry(output):
    transitions = re.findall(r"(\d+) MeshBlocks created,\s*(\d+) deleted by AMR", output)
    if not transitions:
        raise RuntimeError("Missing AMR telemetry in Athena output\n" + output)
    created, deleted = transitions[-1]
    return {"created": int(created), "deleted": int(deleted)}


def _run_athena(label, basename, nlim, restart_path=None):
    _remove_outputs(basename)
    command = ["./athena"]
    if restart_path is None:
        command += ["-i", _athena_input_path()]
    else:
        command += ["-r", os.path.relpath(restart_path, _athena_exe_dir())]
    command += [
        "job/basename=" + basename,
        "time/nlim=" + str(nlim),
        "output1/file_number=0",
        "output2/file_number=0",
        "output3/file_number=0",
    ]

    logger.info("Executing %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Command failed for " + label + "\n" + output)
    return output


def _logical_order(logical):
    return np.lexsort(tuple(logical[:, index] for index in reversed(range(4))))


def _mesh_snapshot_path(path):
    data = bin_convert.read_binary(path)
    logical = np.asarray(data["mb_logical"], dtype=np.int64)
    order = _logical_order(logical)
    state = {
        name: np.asarray(values)[order] for name, values in data["mb_data"].items()
    }
    return {
        "path": path,
        "time": float(data["time"]),
        "cycle": int(data["cycle"]),
        "nmeshblocks": int(data["n_mbs"]),
        "logical": logical[order],
        "geometry": np.asarray(data["mb_geometry"])[order],
        "state": state,
        "finite": all(np.all(np.isfinite(values)) for values in state.values()),
    }


def _mesh_snapshot(basename):
    path = _latest_file("bin", basename + ".mesh.*.bin", "mesh binary")
    return _mesh_snapshot_path(path)


def _particle_snapshot_path(path):
    data = read_particle_vtk(path)
    if "ptag" not in data.scalars:
        raise RuntimeError("Missing particle VTK scalar: ptag")
    order = np.argsort(data.scalars["ptag"])
    return {
        "path": path,
        "points": data.points[order],
        "scalars": {
            name: np.asarray(values)[order] for name, values in data.scalars.items()
        },
        "vectors": {
            name: np.asarray(values)[order] for name, values in data.vectors.items()
        },
    }


def _particle_snapshot(basename):
    path = _latest_file(
        "pvtk", basename + ".prtcl_all.*.part.vtk", "particle VTK"
    )
    return _particle_snapshot_path(path)


def _array_map_parity(left, right):
    if set(left) != set(right):
        return False
    return all(np.array_equal(left[name], right[name]) for name in left)


def _snapshot_parity(left, right):
    mesh_left = left["mesh"]
    mesh_right = right["mesh"]
    particles_left = left["particles"]
    particles_right = right["particles"]
    return {
        "time": mesh_left["time"] == mesh_right["time"],
        "cycle": mesh_left["cycle"] == mesh_right["cycle"],
        "meshblock_count": mesh_left["nmeshblocks"] == mesh_right["nmeshblocks"],
        "logical_topology": np.array_equal(
            mesh_left["logical"], mesh_right["logical"]
        ),
        "mesh_geometry": np.array_equal(
            mesh_left["geometry"], mesh_right["geometry"]
        ),
        "fluid_state": _array_map_parity(mesh_left["state"], mesh_right["state"]),
        "particle_points": np.array_equal(
            particles_left["points"], particles_right["points"]
        ),
        "particle_scalars": _array_map_parity(
            particles_left["scalars"], particles_right["scalars"]
        ),
        "particle_vectors": _array_map_parity(
            particles_left["vectors"], particles_right["vectors"]
        ),
    }


def _all_parity(checks):
    return all(checks.values())


def _snapshot(basename):
    return {
        "mesh": _mesh_snapshot(basename),
        "particles": _particle_snapshot(basename),
    }


def _snapshot_summary(snapshot):
    mesh = snapshot["mesh"]
    particles = snapshot["particles"]
    levels = mesh["logical"][:, 3]
    return {
        "time": mesh["time"],
        "cycle": mesh["cycle"],
        "meshblocks": mesh["nmeshblocks"],
        "min_level": int(np.min(levels)),
        "max_level": int(np.max(levels)),
        "finite_mesh_state": bool(mesh["finite"]),
        "particles": int(particles["points"].shape[0]),
    }


def _summary():
    return {
        "evidence_class": "targeted_local_serial_host_regression",
        "restart_selection": _RESULTS["restart_selection"],
        "seed_amr_telemetry": _RESULTS["seed_amr_telemetry"],
        "seed_initial": _snapshot_summary(_RESULTS["seed_initial"]),
        "seed_terminal": _snapshot_summary(_RESULTS["seed_terminal"]),
        "scheduled_reload_probe": _snapshot_summary(_RESULTS["probe"]),
        "uninterrupted_endpoint": _snapshot_summary(_RESULTS["full"]),
        "resumed_endpoint": _snapshot_summary(_RESULTS["resumed"]),
        "scheduled_reload_probe_parity": _RESULTS["probe_parity"],
        "resumed_endpoint_parity": _RESULTS["endpoint_parity"],
    }


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()

    _run_athena("uninterrupted", _BASENAMES["full"], 2)
    seed_output = _run_athena("transition_seed", _BASENAMES["seed"], 1)

    seed_restarts = _matching_files(
        os.path.join(_athena_exe_dir(), "rst", _BASENAMES["seed"] + ".*.rst"),
        "seed restart",
    )
    scheduled = _restart_path(_BASENAMES["seed"], 1)
    finalize = _restart_path(_BASENAMES["seed"], 2)
    expected = [_restart_path(_BASENAMES["seed"], index) for index in range(3)]
    if seed_restarts != expected:
        raise RuntimeError(
            "Expected initial, scheduled-intermediate, and Finalize checkpoints: "
            + repr(seed_restarts)
        )
    if scheduled == seed_restarts[-1]:
        raise RuntimeError(
            "Scheduled intermediate restart unexpectedly resolved to latest"
        )

    _run_athena("scheduled_reload_probe", _BASENAMES["probe"], 1, scheduled)
    _run_athena("scheduled_resume", _BASENAMES["resumed"], 2, scheduled)

    _RESULTS["restart_selection"] = {
        "selected_intermediate": scheduled,
        "unused_finalize_checkpoint": finalize,
        "latest_seed_checkpoint": seed_restarts[-1],
        "selected_is_latest": scheduled == seed_restarts[-1],
    }
    _RESULTS["seed_amr_telemetry"] = _parse_amr_telemetry(seed_output)
    _RESULTS["seed_initial"] = _mesh_and_particle_snapshot_numbers(
        _BASENAMES["seed"], 0
    )
    _RESULTS["seed_terminal"] = _snapshot(_BASENAMES["seed"])
    _RESULTS["probe"] = _snapshot(_BASENAMES["probe"])
    _RESULTS["full"] = _snapshot(_BASENAMES["full"])
    _RESULTS["resumed"] = _snapshot(_BASENAMES["resumed"])
    _RESULTS["probe_parity"] = _snapshot_parity(
        _RESULTS["probe"], _RESULTS["seed_terminal"]
    )
    _RESULTS["endpoint_parity"] = _snapshot_parity(
        _RESULTS["resumed"], _RESULTS["full"]
    )


def _mesh_and_particle_snapshot_numbers(basename, number):
    return {
        "mesh": _mesh_snapshot_path(
            os.path.join(
                _athena_exe_dir(), "bin", f"{basename}.mesh.{number:05d}.bin"
            )
        ),
        "particles": _particle_snapshot_path(
            os.path.join(
                _athena_exe_dir(),
                "pvtk",
                f"{basename}.prtcl_all.{number:05d}.part.vtk",
            )
        ),
    }


def analyze():
    logger.debug("Analyzing test " + __name__)
    summary = _summary()
    logger.info("Q-009 periodic AMR restart ordering parity metrics: %s", summary)
    selection = summary["restart_selection"]
    telemetry = summary["seed_amr_telemetry"]
    initial = summary["seed_initial"]
    terminal = summary["seed_terminal"]
    probe = summary["scheduled_reload_probe"]
    full = summary["uninterrupted_endpoint"]
    resumed = summary["resumed_endpoint"]

    return (
        not selection["selected_is_latest"]
        and selection["selected_intermediate"] != selection["unused_finalize_checkpoint"]
        and telemetry["created"] > 0
        and terminal["meshblocks"] > initial["meshblocks"]
        and terminal["max_level"] > initial["max_level"]
        and terminal["finite_mesh_state"]
        and probe["finite_mesh_state"]
        and full["finite_mesh_state"]
        and resumed["finite_mesh_state"]
        and initial["particles"] > 0
        and probe["particles"] == initial["particles"]
        and resumed["particles"] == full["particles"]
        and _all_parity(summary["scheduled_reload_probe_parity"])
        and _all_parity(summary["resumed_endpoint_parity"])
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
