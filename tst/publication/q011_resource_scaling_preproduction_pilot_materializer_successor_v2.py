#!/usr/bin/env python3
"""Materialize exact-nine-output Q-011 scaling and I/O pilot review records."""

from __future__ import annotations

import hashlib
from typing import Mapping

from tst.publication import q011_section54_production_campaign_contract_successor_v2 as contract
from tst.publication import q011_section54_production_campaign_planner_successor_v2 as planner


RECORD_TYPE = "q011_section54_resource_pilot_materialization_successor_v2"
NODE_LADDERS = {
    "coarse_uniform_dx12": [4, 8, 16],
    "three_level_amr_root_dx12_finest_dx3": [8, 16, 32],
    "fine_uniform_dx3": [32, 64, 128],
}
PHASES = ("compute_scaling", "exact_output_io_scaling")
PHASE_SEEDS = {"compute_scaling": 24060601, "exact_output_io_scaling": 24060604}


class ResourcePilotMaterializationError(ValueError):
    """Reject a drifted or authority-bearing exact-deck pilot materialization."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ResourcePilotMaterializationError(message)


def _binding(path: str, payload: bytes) -> dict[str, str]:
    return {"path": path, "sha256": hashlib.sha256(payload).hexdigest()}


def _variant_overrides(variant: str) -> dict[tuple[str, str], str]:
    result: dict[tuple[str, str], str] = {}
    for token in contract.VARIANT_OVERRIDES[variant]:
        key, value = token.split("=", 1)
        block, name = key.split("/", 1)
        result[(block, name)] = value
    return result


def _overrides(phase: str, variant: str, basename: str) -> dict[tuple[str, str], str]:
    _require(phase in PHASES, f"unknown Q011 resource-pilot phase {phase!r}")
    _require(variant in contract.GRID_VARIANTS, f"unknown Q011 resource-pilot variant {variant!r}")
    seed = PHASE_SEEDS[phase]
    result = {
        ("job", "basename"): basename,
        ("problem", "ps_p0"): "1.0",
        ("particles", "pic_random_seed"): str(seed),
        ("problem", "ps_inject_seed"): str(seed),
        ("problem", "ps_seed_noise_seed"): str(seed),
        **_variant_overrides(variant),
    }
    if phase == "compute_scaling":
        result[("time", "nlim")] = "512"
        for index in range(1, 10):
            result[(f"output{index}", "dt")] = "-1.0"
    else:
        result[("mesh", "x2max")] = "240.0"
        result[("mesh", "nx2")] = (
            "80" if variant == "fine_uniform_dx3" else "20"
        )
    return result


def _validate_pilot_deck(payload: bytes, *, phase: str, variant: str, seed: int) -> None:
    blocks = contract.parse_deck(payload.decode("utf-8"))
    outputs = {name: values for name, values in blocks.items() if name.startswith("output")}
    _require(set(outputs) == set(contract.EXPECTED_OUTPUTS), "pilot deck output block inventory drifted")
    _require(blocks["problem"]["ps_p0"] == "1.0", "pilot deck pressure drifted")
    _require(
        {
            int(blocks["particles"]["pic_random_seed"]),
            int(blocks["problem"]["ps_inject_seed"]),
            int(blocks["problem"]["ps_seed_noise_seed"]),
        }
        == {seed},
        "pilot deck engineering seed drifted",
    )
    _require(seed in contract.ENGINEERING_SEEDS, "pilot deck does not use an engineering seed")
    _require(seed not in contract.CORE_QUALIFYING_SEEDS, "pilot deck uses a qualifying seed")
    if phase == "compute_scaling":
        _require(blocks["time"]["nlim"] == "512", "compute-scaling cycle limit drifted")
        _require(
            all(float(outputs[f"output{index}"]["dt"]) <= 0.0 for index in range(1, 10)),
            "compute-scaling deck did not suppress every successor output",
        )
        _require(blocks["mesh"]["x2max"] == "3120.0", "compute-scaling geometry drifted")
    else:
        _require(blocks["time"]["tlim"] == "1200.0", "I/O pilot terminal time drifted")
        _require(
            outputs == contract.EXPECTED_OUTPUTS,
            "I/O pilot does not retain the exact nine-output successor cadence",
        )
        _require(blocks["mesh"]["x2max"] == "240.0", "I/O pilot transverse extent drifted")
        _require(
            blocks["mesh"]["nx2"] == ("80" if variant == "fine_uniform_dx3" else "20"),
            "I/O pilot transverse resolution drifted",
        )


def build_materialization(
    plan: object | None = None,
) -> tuple[dict[str, object], dict[str, bytes]]:
    """Build deterministic exact-deck pilot records without writing or launching."""
    normalized_plan = planner.validate_plan(planner.build_plan() if plan is None else plan)
    stage = contract.build_stage_contract("resource_pilot_materializer")
    files: dict[str, bytes] = {}
    decks: list[dict[str, object]] = []
    cases: list[dict[str, object]] = []
    for phase in PHASES:
        seed = PHASE_SEEDS[phase]
        for variant in contract.GRID_VARIANTS:
            basename = f"q011-{phase.replace('_', '-')}-{variant}-seed-{seed}"
            payload = contract.render_successor_deck(_overrides(phase, variant, basename))
            _validate_pilot_deck(payload, phase=phase, variant=variant, seed=seed)
            relative = f"decks/{basename}.athinput"
            files[relative] = payload
            binding = _binding(relative, payload)
            decks.append(
                {
                    "phase": phase,
                    "variant": variant,
                    "engineering_seed": seed,
                    **binding,
                }
            )
            cases.append(
                {
                    "case_id": basename,
                    "phase": phase,
                    "variant": variant,
                    "engineering_seed": seed,
                    "input_deck": binding,
                    "node_candidates": (
                        NODE_LADDERS[variant] if phase == "compute_scaling" else []
                    ),
                    "status": (
                        "source_local_node_ladder_review_only"
                        if phase == "compute_scaling"
                        else "blocked_pending_compute_scaling_node_selection"
                    ),
                    "physical_output_scientific_inspection_authorized": False,
                    "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
                }
            )
    manifest = {
        "record_type": RECORD_TYPE,
        "schema_version": contract.SCHEMA_VERSION,
        "status": "source_local_materialization_complete_execution_blocked",
        "stage_contract": stage,
        "plan_sha256": contract.canonical_sha256(normalized_plan),
        "production_deck": stage["production_deck"],
        "successor_output_count": 9,
        "phases": list(PHASES),
        "engineering_seeds": dict(PHASE_SEEDS),
        "engineering_seeds_disjoint_from_qualifying_seeds": True,
        "maximum_consumed_node_hours": {
            "value": 500.0,
            "source_category": "engineering closure",
            "rationale": (
                "The cap is the preregistered Q011 engineering-pilot budget inside "
                "the project-wide node-hour ceiling and has no scientific literature "
                "origin."
            ),
        },
        "deck_count": len(decks),
        "decks": decks,
        "case_count": len(cases),
        "cases": cases,
        "gate_id": "exact_successor_deck_scaling_and_io_pilots",
        "gate_status": "blocked_pending_registered_runtime_evidence",
        "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
    }
    return manifest, files


def validate_materialization(
    manifest: object, files: Mapping[str, bytes], plan: object | None = None
) -> tuple[dict[str, object], dict[str, bytes]]:
    expected_manifest, expected_files = build_materialization(plan)
    _require(
        contract._strict_equal(manifest, expected_manifest),
        "Q011 resource-pilot manifest drifted or acquired authority",
    )
    _require(
        type(files) is dict
        and files.keys() == expected_files.keys()
        and all(type(files[path]) is bytes and files[path] == expected_files[path] for path in files),
        "Q011 resource-pilot deck bytes drifted",
    )
    return expected_manifest, expected_files


def main() -> None:
    manifest, _ = build_materialization()
    print(contract.canonical_json_bytes(manifest).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
