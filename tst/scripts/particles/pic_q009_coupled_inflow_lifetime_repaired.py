"""Repaired bounded serial-host Q-009 coupled inflow lifetime regression."""

from __future__ import annotations

import json
import logging
import math

from scripts.particles import pic_q009_coupled_inflow_lifetime as characterized

logger = logging.getLogger("athena" + __name__[7:])

run = characterized.run
_summary = characterized._summary
_INFLOW_CONTRACT = characterized._INFLOW_CONTRACT


def analyze():
    logger.debug("Analyzing test " + __name__)
    summary = _summary()
    logger.info("Q-009 repaired coupled inflow lifetime metrics: %s", summary)
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
        and len(set(particles["particle_counts"])) == 1
        and not particles["physical_boundary_destruction_observed"]
        and particles["empty_population_restart_continuations"] == 0
        and all(
            count > 0
            for count in particles[
                "ownership_changes_per_transition_for_retained_particles"
            ]
        )
        and particles["reflect_boundary_observed"]
        and math.isfinite(times[-1])
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
