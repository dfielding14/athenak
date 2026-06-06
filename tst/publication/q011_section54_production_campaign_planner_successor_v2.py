#!/usr/bin/env python3
"""Source-local Q-011 nine-output production-campaign planner successor."""

from __future__ import annotations

from typing import Any

from tst.publication import q011_section54_production_campaign_contract_successor_v2 as contract


RECORD_TYPE = "q011_section54_production_campaign_plan_successor_v2"


class ProductionCampaignPlannerError(ValueError):
    """Reject a drifted or authority-bearing source-local plan."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ProductionCampaignPlannerError(message)


def _attempt_id(index: int, variant: str, seed: int) -> str:
    return f"baseline-{index:03d}-{variant}-seed-{seed}"


def build_plan() -> dict[str, object]:
    """Build the nine-attempt core plan without execution authority."""
    stage = contract.build_stage_contract("planner")
    attempts: list[dict[str, object]] = []
    index = 0
    for variant in contract.GRID_VARIANTS:
        for seed in contract.CORE_QUALIFYING_SEEDS:
            index += 1
            attempt_id = _attempt_id(index, variant, seed)
            attempts.append(
                {
                    "attempt_id": attempt_id,
                    "variant": variant,
                    "qualifying_seed": seed,
                    "production_deck": stage["production_deck"],
                    "argv": contract.expected_attempt_argv(attempt_id, variant, seed),
                    "status": "source_local_launch_prohibited_handoff",
                    "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
                }
            )
    return {
        "record_type": RECORD_TYPE,
        "schema_version": contract.SCHEMA_VERSION,
        "campaign_id": contract.CAMPAIGN_ID,
        "status": "source_local_plan_complete_execution_blocked",
        "stage_contract": stage,
        "selected_pressure": dict(contract.SELECTED_PRESSURE),
        "paired_seed_rule": "same seed across the complete coarse-AMR-fine triad",
        "attempt_count": len(attempts),
        "attempts": attempts,
        "blockers": contract.source_local_blockers(),
        "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
    }


def validate_plan(value: object) -> dict[str, object]:
    """Require the exact deterministic source-local core plan."""
    expected = build_plan()
    _require(
        contract._strict_equal(value, expected),
        "Q011 production-campaign planner successor drifted or acquired authority",
    )
    return expected


def main() -> None:
    print(contract.canonical_json_bytes(build_plan()).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
