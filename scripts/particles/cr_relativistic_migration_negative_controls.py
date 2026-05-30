#!/usr/bin/env python3
"""Reject deliberate Phase-8 registration and migration-oracle weakenings."""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path
from typing import Callable

import cr_relativistic_migration_runtime_inspect as phase8


def _require_rejection(label: str, action: Callable[[], object], expected: str) -> None:
    try:
        action()
    except (RuntimeError, ValueError) as error:
        if expected not in str(error):
            raise RuntimeError(
                f"{label}: rejection {error!r} lacks expected substring {expected!r}"
            ) from error
        print(f"rejected {label}: {error}")
        return
    raise RuntimeError(f"{label}: mutation was accepted")


def _particle(
    tag: int, *, width: int = 22, value: float = 0.0, pgid: int = 0
) -> phase8.ParticleState:
    return phase8.ParticleState(
        pgid=pgid,
        tag=tag,
        species=0,
        reals=tuple(value for _ in range(width)),
    )


def _inspection(particles: list[phase8.ParticleState]) -> phase8.Inspection:
    checkpoint = phase8.Checkpoint(
        root=Path("."),
        number=0,
        manifest_relative=Path("prst/control.00000.prst.manifest"),
        mesh_relative=Path("rst/control.00000.rst"),
        witness_relative=Path("rst/control.00000.rst.rmeta"),
        shard_relatives=(Path("prst/rank_00000000/control.00000.prst"),),
        rank_counts=(len(particles),),
        topology_hash=1,
    )
    keyed = {particle.key: particle for particle in particles}
    return phase8.Inspection(checkpoint=checkpoint, particles=keyed)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="cr-rel-phase8-negative-") as temporary:
        root = Path(temporary)
        criteria = json.loads(phase8.DEFAULT_CRITERIA.read_text())

        tampered = json.loads(json.dumps(criteria))
        tampered["criteria"][0]["limit"] = 21
        tampered_path = root / "tampered-criteria.json"
        tampered_path.write_text(json.dumps(tampered, indent=2) + "\n")
        _require_rejection(
            "tampered criteria digest",
            lambda: phase8._validate_criteria_document(tampered_path),
            "criteria SHA-256",
        )

        reordered = json.loads(json.dumps(criteria))
        reordered["criteria"][0], reordered["criteria"][1] = (
            reordered["criteria"][1],
            reordered["criteria"][0],
        )
        reordered_path = root / "reordered-criteria.json"
        reordered_path.write_text(json.dumps(reordered, indent=2) + "\n")
        _require_rejection(
            "reordered criterion IDs",
            lambda: phase8._validate_criteria_document(
                reordered_path, enforce_digest=False),
            "criterion ID sequence",
        )

        correct = [_particle(tag) for tag in range(8)]
        _require_rejection(
            "typed-v2 width loss",
            lambda: phase8._validate_particle_states(
                [*correct[:-1], _particle(7, width=21)],
                expected_keys=phase8.EXPECTED_PARTICLE_KEYS,
            ),
            "typed-v2 real width",
        )
        _require_rejection(
            "non-finite typed-v2 state",
            lambda: phase8._validate_particle_states(
                [*correct[:-1], _particle(7, value=math.nan)],
                expected_keys=phase8.EXPECTED_PARTICLE_KEYS,
            ),
            "non-finite",
        )
        _require_rejection(
            "duplicate typed identity",
            lambda: phase8._validate_particle_states(
                [*correct[:-1], _particle(6)],
                expected_keys=phase8.EXPECTED_PARTICLE_KEYS,
            ),
            "duplicate particle identity",
        )
        _require_rejection(
            "lost typed identity",
            lambda: phase8._validate_particle_states(
                correct[:-1], expected_keys=phase8.EXPECTED_PARTICLE_KEYS),
            "particle identities",
        )

        left = _inspection(correct)
        shifted = [
            phase8.ParticleState(
                pgid=particle.pgid,
                tag=particle.tag,
                species=particle.species,
                reals=(1.0e-3, *particle.reals[1:]),
            )
            for particle in correct
        ]
        _require_rejection(
            "physical-state parity drift",
            lambda: phase8._compare(left, _inspection(shifted), "negative parity"),
            "IPX changed",
        )

        _require_rejection(
            "periodic-wrap witness loss",
            lambda: phase8._periodic_wrap_count(_inspection(correct)),
            "wrapped coordinate disagrees",
        )

        metrics = {
            criterion["metric"]: criterion["limit"]
            for criterion in criteria["criteria"]
        }
        metrics["typed_v2_real_field_count"] = 21
        report = {"metrics": metrics}
        _require_rejection(
            "runtime criterion weakening",
            lambda: phase8._validate_runtime_criteria(phase8.DEFAULT_CRITERIA, report),
            "P8-01",
        )

        for label, metric, value, criterion_id in (
            (
                "prescribed sampled-field witness weakening",
                "prescribed_sampled_field_max_real_error",
                1.0e-3,
                "P8-03",
            ),
            (
                "static-SMR multilevel exchange witness weakening",
                "static_smr_mpi_exchange_sent",
                0,
                "P8-19",
            ),
            (
                "static-SMR host-tree cross-level witness weakening",
                "static_smr_host_tree_cross_level_transition_count",
                0,
                "P8-61",
            ),
            (
                "adaptive off-rank remap witness weakening",
                "adaptive_device_off_rank_remap_sent",
                0,
                "P8-27",
            ),
        ):
            weakened = {
                criterion["metric"]: criterion["limit"]
                for criterion in criteria["criteria"]
            }
            weakened[metric] = value
            _require_rejection(
                label,
                lambda weakened=weakened: phase8._validate_runtime_criteria(
                    phase8.DEFAULT_CRITERIA, {"metrics": weakened}),
                criterion_id,
            )

        _require_rejection(
            "missing AMR-remap performance witness",
            lambda: phase8._performance_counters("", "amr_remap"),
            "lacks Particle performance",
        )
        _require_rejection(
            "initialization-only AMR-remap performance witness",
            lambda: phase8._performance_counters(
                "Particle performance amr_remap: sent=2 received=2\n",
                "amr_remap",
                after_setup=True,
            ),
            "setup-complete marker",
        )
        _require_rejection(
            "initialization-only AMR churn witness",
            lambda: phase8._require_runtime_amr_churn(
                "Total number of MeshBlocks = 8\n"
                "56 MeshBlocks created, 21 deleted by AMR\n"
                "Setup complete, executing task list(s)...\n",
                "negative AMR churn",
            ),
            "no post-setup MeshBlocks created",
        )
        _require_rejection(
            "missing static-SMR cross-level witness",
            lambda: phase8._performance_counters(
                "Setup complete, executing task list(s)...\n",
                "particle_multilevel_lookup",
                after_setup=True,
            ),
            "lacks Particle performance",
        )

        aggregated = phase8._performance_counters(
            "Particle performance amr_remap: sent=2 received=1\n"
            "Particle performance amr_remap: sent=3 received=4\n",
            "amr_remap",
        )
        if aggregated != {"sent": 5, "received": 5}:
            raise RuntimeError(f"performance aggregation changed: {aggregated!r}")

        mhd_ideal = phase8._mhd_ideal_source_text(phase8.DEFAULT_INPUT)
        if "cE0x" in mhd_ideal or "cE0y" in mhd_ideal or "cE0z" in mhd_ideal:
            raise RuntimeError("mhd_ideal negative-control input retained cE fields")
        if "relativistic_temporal_sampling = frozen_tn" not in mhd_ideal:
            raise RuntimeError("mhd_ideal negative-control input lacks frozen_tn")
        if "rsolver     = hlle" not in mhd_ideal:
            raise RuntimeError("mhd_ideal negative-control input lacks dynamic solver")

    print("CR relativistic Phase-8 mutation controls passed: 17 rejected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
