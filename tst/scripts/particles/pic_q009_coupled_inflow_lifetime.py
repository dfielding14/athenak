"""Bounded serial-host Q-009 coupled inflow lifetime sibling regression."""

from __future__ import annotations

from contextlib import contextmanager
import json
import logging
import math
import os
import subprocess

import numpy as np

from scripts.particles import pic_q009_coupled_boundary_lifetime as archived

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_q009_coupled_inflow_lifetime.athinput"
_EXE_ENV = "ATHENA_Q009_COUPLED_INFLOW_EXE_DIR"
_ARCHIVED_EXE_ENV = "ATHENA_Q009_COUPLED_BOUNDARY_EXE_DIR"
_BASENAME_PREFIX = "pic_q009_coupled_inflow_"
_INFLOW_CONTRACT = "zero_valued_moment_ghosts"
_ARCHIVED_SUMMARY = archived._summary
_RESULTS = archived._RESULTS


@contextmanager
def _adapt_archived_harness():
    original_deck = archived._INPUT_DECK
    original_exe_dir = os.environ.get(_ARCHIVED_EXE_ENV)
    archived._INPUT_DECK = _INPUT_DECK
    if _EXE_ENV in os.environ:
        os.environ[_ARCHIVED_EXE_ENV] = os.environ[_EXE_ENV]
    try:
        yield
    finally:
        archived._INPUT_DECK = original_deck
        if original_exe_dir is None:
            os.environ.pop(_ARCHIVED_EXE_ENV, None)
        else:
            os.environ[_ARCHIVED_EXE_ENV] = original_exe_dir


def _run_stage(stage, restart_path):
    label, transition, threshold, nlim = stage
    basename = _BASENAME_PREFIX + label
    archived._remove_outputs(basename)

    command = ["./athena"]
    if restart_path is None:
        command += ["-i", archived._athena_input_path()]
    else:
        command += ["-r", os.path.relpath(restart_path, archived._athena_exe_dir())]
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
        command, cwd=archived._athena_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Command failed for " + label + "\n" + output)

    telemetry = archived._parse_telemetry(output)
    mesh = archived._mesh_snapshot(basename)
    particles = archived._particle_snapshot(basename)
    return {
        "label": label,
        "transition": transition,
        "threshold": threshold,
        "nlim": nlim,
        "runtime_identity": archived._parse_runtime_identity(output),
        "telemetry": telemetry,
        "mesh": mesh,
        "particles": particles,
        "restart": archived._latest_restart(basename),
    }


def _summary():
    summary = _ARCHIVED_SUMMARY()
    summary["evidence_class"] = "bounded_local_serial_host_inflow_sibling_regression"
    summary["supported_physical_boundaries"]["x2"] = "inflow"
    summary["inflow_contract"] = _INFLOW_CONTRACT
    return summary


def _particle_integrity(initial, stages):
    snapshots = [initial] + [stage["particles"] for stage in stages]
    meshblock_limits = [min(stage["mesh"]["nmeshblocks"] for stage in stages)]
    meshblock_limits += [stage["mesh"]["nmeshblocks"] for stage in stages]
    lower = np.array([0.0, 0.0, 0.0])
    upper = np.array([1.0, 4.0, 1.0])
    checks = {
        "initial_count": int(initial["ptag"].size),
        "tags_unique": True,
        "retained_tags_stable": True,
        "finite_points": True,
        "finite_velocity": True,
        "in_domain": True,
        "species_valid": True,
        "gid_valid": True,
    }
    counts = []
    destruction_counts = []
    ownership_changes = []
    reflected_tags = set()
    prior_tags = None
    prior_gid_by_tag = None
    for meshblock_limit, snapshot in zip(meshblock_limits, snapshots):
        tags = snapshot["ptag"]
        gid_by_tag = dict(zip(tags.tolist(), snapshot["gid"].tolist()))
        tag_set = set(gid_by_tag)
        counts.append(int(tags.size))
        checks["tags_unique"] &= len(tag_set) == tags.size
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
        if prior_tags is not None:
            checks["retained_tags_stable"] &= tag_set <= prior_tags
            destruction_counts.append(len(prior_tags) - len(tag_set))
            ownership_changes.append(
                sum(
                    prior_gid_by_tag[tag] != gid_by_tag[tag]
                    for tag in tag_set
                )
            )
        prior_tags = tag_set
        prior_gid_by_tag = gid_by_tag

    empty_indices = [index for index, count in enumerate(counts[1:]) if count == 0]
    first_empty_stage = empty_indices[0] if empty_indices else None
    continuations_after_empty = (
        len(stages) - first_empty_stage - 1 if first_empty_stage is not None else 0
    )
    return {
        **{name: bool(value) if name != "initial_count" else value
           for name, value in checks.items()},
        "particle_counts": counts,
        "particle_count_nonincreasing": all(
            after <= before for before, after in zip(counts[:-1], counts[1:])
        ),
        "physical_boundary_destruction_counts": destruction_counts,
        "physical_boundary_destruction_observed": any(
            count > 0 for count in destruction_counts
        ),
        "ownership_changes_per_transition_for_retained_particles": ownership_changes,
        "empty_population_restart_continuations": continuations_after_empty,
        "reflected_particle_count": len(reflected_tags),
        "reflect_boundary_observed": len(reflected_tags) > 0,
    }


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    stages = []
    restart_path = None
    with _adapt_archived_harness():
        for stage in archived._STAGES:
            result = _run_stage(stage, restart_path)
            stages.append(result)
            restart_path = result["restart"]
        _RESULTS["stages"] = stages
        initial = archived._initial_particle_snapshot(_BASENAME_PREFIX + "refine_1")
        _RESULTS["particles"] = _particle_integrity(initial, stages)


def analyze():
    logger.debug("Analyzing test " + __name__)
    summary = _summary()
    logger.info("Q-009 coupled inflow lifetime metrics: %s", summary)
    particles = summary["particles"]
    stages = summary["stages"]
    levels = [stage["max_level"] for stage in stages]
    meshblocks = [stage["meshblocks"] for stage in stages]
    times = [stage["time"] for stage in stages]

    return (
        summary["runtime_identity_stable"]
        and summary["supported_physical_boundaries"]["x2"] == "inflow"
        and summary["inflow_contract"] == _INFLOW_CONTRACT
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
        and particles["initial_count"] > 0
        and particles["tags_unique"]
        and particles["retained_tags_stable"]
        and particles["finite_points"]
        and particles["finite_velocity"]
        and particles["in_domain"]
        and particles["species_valid"]
        and particles["gid_valid"]
        and particles["particle_count_nonincreasing"]
        and particles["physical_boundary_destruction_observed"]
        and particles["empty_population_restart_continuations"] > 0
        and particles["reflect_boundary_observed"]
        and math.isfinite(times[-1])
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
