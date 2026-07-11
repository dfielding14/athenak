"""Bounded serial-host Q-009 coupled-boundary particle lifetime regression."""

from __future__ import annotations

import glob
import json
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
_VIS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "vis", "python")
)
sys.path.insert(0, _PUBLICATION_DIR)
sys.path.insert(0, _VIS_DIR)

import bin_convert_new as bin_convert  # noqa: E402
from pvtk_particles import read_particle_vtk  # noqa: E402

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_q009_coupled_boundary_lifetime_vl2_tsc.athinput"
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_STAGES = [
    ("refine_1", "refine", 0.5, 2),
    ("derefine_1", "derefine", 1.0e6, 4),
    ("refine_2", "refine", 0.5, 6),
    ("derefine_2", "derefine", 1.0e6, 8),
    ("refine_3", "refine", 0.5, 10),
    ("derefine_3", "derefine", 1.0e6, 12),
]
_RUNTIME_IDENTITY_TOKENS = [
    "physical_mode=paper_mhd_pic_vl2_tsc",
    "background=coupled",
    "feedback=coupled",
    "induction=ideal_mhd_only",
    "restart_schema=8",
]
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        "ATHENA_Q009_COUPLED_BOUNDARY_EXE_DIR",
        os.path.join(os.getcwd(), "build", "src"),
    )


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


def _remove_outputs(basename):
    exe_dir = _athena_exe_dir()
    for dirname in ["bin", "pvtk", "rst"]:
        pattern = os.path.join(exe_dir, dirname, basename + ".*")
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


def _latest_restart(basename):
    return _latest_file("rst", basename + ".*.rst", "restart")


def _parse_telemetry(output):
    current = re.search(r"Current number of MeshBlocks =\s*(\d+)", output)
    transitions = re.search(
        r"(\d+) MeshBlocks created,\s*(\d+) deleted by AMR", output
    )
    if current is None or transitions is None:
        raise RuntimeError("Missing AMR telemetry in Athena output\n" + output)
    return {
        "meshblocks": int(current.group(1)),
        "created": int(transitions.group(1)),
        "deleted": int(transitions.group(2)),
    }


def _parse_runtime_identity(output):
    lines = [
        line for line in output.splitlines() if line.startswith("PIC runtime model:")
    ]
    if not lines:
        raise RuntimeError("Missing PIC runtime model identity in Athena output\n" + output)
    identity = lines[-1]
    missing = [token for token in _RUNTIME_IDENTITY_TOKENS if token not in identity]
    if missing:
        raise RuntimeError(
            "PIC runtime model identity missing expected tokens "
            + repr(missing)
            + "\n"
            + identity
        )
    return identity


def _mesh_snapshot(basename):
    path = _latest_file("bin", basename + ".mesh.*.bin", "mesh binary")
    data = bin_convert.read_binary(path)
    levels = np.asarray(data["mb_logical"][:, 3], dtype=np.int64)
    finite = all(
        np.all(np.isfinite(values))
        for arrays in data["mb_data"].values()
        for values in arrays
    )
    return {
        "path": path,
        "time": float(data["time"]),
        "nmeshblocks": int(data["n_mbs"]),
        "min_level": int(np.min(levels)),
        "max_level": int(np.max(levels)),
        "finite": bool(finite),
    }


def _particle_snapshot_path(path):
    data = read_particle_vtk(path)
    for name in ["gid", "ptag", "species"]:
        if name not in data.scalars:
            raise RuntimeError("Missing particle VTK scalar: " + name)
    if "vel" not in data.vectors:
        raise RuntimeError("Missing particle VTK vector: vel")
    order = np.argsort(data.scalars["ptag"])
    return {
        "path": path,
        "points": data.points[order],
        "velocity": data.vectors["vel"][order],
        "gid": data.scalars["gid"][order],
        "ptag": data.scalars["ptag"][order],
        "species": data.scalars["species"][order],
    }


def _particle_snapshot(basename):
    path = _latest_file(
        "pvtk", basename + ".prtcl_all.*.part.vtk", "particle VTK"
    )
    return _particle_snapshot_path(path)


def _initial_particle_snapshot(basename):
    paths = _matching_files(
        os.path.join(_athena_exe_dir(), "pvtk", basename + ".prtcl_all.*.part.vtk"),
        "particle VTK",
    )
    return _particle_snapshot_path(paths[0])


def _run_stage(stage, restart_path):
    label, transition, threshold, nlim = stage
    basename = "pic_q009_coupled_boundary_" + label
    _remove_outputs(basename)

    command = ["./athena"]
    if restart_path is None:
        command += ["-i", _athena_input_path()]
    else:
        command += ["-r", os.path.relpath(restart_path, _athena_exe_dir())]
    command += [
        "job/basename=" + basename,
        "mesh_refinement/dens_max=" + str(threshold),
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

    telemetry = _parse_telemetry(output)
    mesh = _mesh_snapshot(basename)
    particles = _particle_snapshot(basename)
    return {
        "label": label,
        "transition": transition,
        "threshold": threshold,
        "nlim": nlim,
        "runtime_identity": _parse_runtime_identity(output),
        "telemetry": telemetry,
        "mesh": mesh,
        "particles": particles,
        "restart": _latest_restart(basename),
    }


def _particle_integrity(initial, stages):
    snapshots = [initial] + [stage["particles"] for stage in stages]
    meshblock_limits = [min(stage["mesh"]["nmeshblocks"] for stage in stages)]
    meshblock_limits += [stage["mesh"]["nmeshblocks"] for stage in stages]
    first = snapshots[0]
    count = first["ptag"].size
    first_tags = first["ptag"]
    checks = {
        "count": int(count),
        "count_stable": True,
        "tags_unique": True,
        "tags_stable": True,
        "finite_points": True,
        "finite_velocity": True,
        "in_domain": True,
        "species_valid": True,
        "gid_valid": True,
    }
    transition_counts = []
    reflected_tags = set()
    lower = np.array([0.0, 0.0, 0.0])
    upper = np.array([1.0, 4.0, 1.0])
    for meshblock_limit, snapshot in zip(meshblock_limits, snapshots):
        checks["count_stable"] &= snapshot["ptag"].size == count
        checks["tags_unique"] &= np.unique(snapshot["ptag"]).size == count
        checks["tags_stable"] &= np.array_equal(snapshot["ptag"], first_tags)
        checks["finite_points"] &= bool(np.all(np.isfinite(snapshot["points"])))
        checks["finite_velocity"] &= bool(np.all(np.isfinite(snapshot["velocity"])))
        checks["in_domain"] &= bool(
            np.all((snapshot["points"] >= lower) & (snapshot["points"] <= upper))
        )
        checks["species_valid"] &= bool(np.all(snapshot["species"] == 0))
        checks["gid_valid"] &= bool(
            np.all((snapshot["gid"] >= 0) & (snapshot["gid"] < meshblock_limit))
        )
        reflected_tags.update(snapshot["ptag"][snapshot["velocity"][:, 0] < 0])

    for before, after in zip(snapshots[:-1], snapshots[1:]):
        transition_counts.append(int(np.count_nonzero(before["gid"] != after["gid"])))

    root_snapshots = [
        stage["particles"] for stage in stages if stage["transition"] == "derefine"
    ]
    root_gid_matches_position = all(
        np.array_equal(snapshot["gid"], (snapshot["points"][:, 0] >= 0.5).astype(int))
        for snapshot in root_snapshots
    )
    root_migration_count = int(
        np.count_nonzero(root_snapshots[0]["gid"] != root_snapshots[-1]["gid"])
    )
    point_displacement = float(
        np.max(np.abs(root_snapshots[0]["points"] - root_snapshots[-1]["points"]))
    )
    return {
        **{name: bool(value) if name != "count" else value for name, value in checks.items()},
        "ownership_changes_per_amr_transition": transition_counts,
        "all_amr_transitions_refresh_ownership": all(value > 0 for value in transition_counts),
        "root_gid_matches_position": root_gid_matches_position,
        "same_level_root_migration_count": root_migration_count,
        "same_level_root_max_point_displacement": point_displacement,
        "reflected_particle_count": len(reflected_tags),
        "reflect_boundary_observed": len(reflected_tags) > 0,
    }


def _summary():
    stages = _RESULTS["stages"]
    return {
        "evidence_class": "bounded_local_serial_host_regression",
        "not_qualification_evidence": True,
        "runtime_identity": stages[0]["runtime_identity"],
        "runtime_identity_stable": len(
            {stage["runtime_identity"] for stage in stages}
        )
        == 1,
        "supported_physical_boundaries": {
            "x1": "reflect",
            "x2": "outflow",
            "x3": "periodic",
        },
        "restart_continuations": len(stages) - 1,
        "stages": [
            {
                "label": stage["label"],
                "transition": stage["transition"],
                "time": stage["mesh"]["time"],
                "meshblocks": stage["mesh"]["nmeshblocks"],
                "min_level": stage["mesh"]["min_level"],
                "max_level": stage["mesh"]["max_level"],
                "finite_mesh_state": stage["mesh"]["finite"],
                "created": stage["telemetry"]["created"],
                "deleted": stage["telemetry"]["deleted"],
            }
            for stage in stages
        ],
        "particles": _RESULTS["particles"],
    }


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    stages = []
    restart_path = None
    for stage in _STAGES:
        result = _run_stage(stage, restart_path)
        stages.append(result)
        restart_path = result["restart"]
    _RESULTS["stages"] = stages
    initial = _initial_particle_snapshot("pic_q009_coupled_boundary_refine_1")
    _RESULTS["particles"] = _particle_integrity(initial, stages)


def analyze():
    logger.debug("Analyzing test " + __name__)
    summary = _summary()
    logger.info("Q-009 coupled-boundary lifetime metrics: %s", summary)
    particles = summary["particles"]
    stages = summary["stages"]
    levels = [stage["max_level"] for stage in stages]
    meshblocks = [stage["meshblocks"] for stage in stages]
    times = [stage["time"] for stage in stages]

    return (
        summary["runtime_identity_stable"]
        and summary["restart_continuations"] == 5
        and len(stages) == 6
        and all(after > before for before, after in zip(times[:-1], times[1:]))
        and all(stage["finite_mesh_state"] for stage in stages)
        and len(set(levels)) == 2
        and levels[0::2] == [max(levels)] * 3
        and levels[1::2] == [min(levels)] * 3
        and meshblocks[0::2] == [max(meshblocks)] * 3
        and meshblocks[1::2] == [min(meshblocks)] * 3
        and all(stage["created"] > 0 for stage in stages[0::2])
        and all(stage["deleted"] > 0 for stage in stages[1::2])
        and particles["count"] > 0
        and particles["count_stable"]
        and particles["tags_unique"]
        and particles["tags_stable"]
        and particles["finite_points"]
        and particles["finite_velocity"]
        and particles["in_domain"]
        and particles["species_valid"]
        and particles["gid_valid"]
        and particles["all_amr_transitions_refresh_ownership"]
        and particles["root_gid_matches_position"]
        and particles["same_level_root_migration_count"] > 0
        and math.isfinite(particles["same_level_root_max_point_displacement"])
        and particles["same_level_root_max_point_displacement"] > 0.0
        and particles["reflect_boundary_observed"]
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
