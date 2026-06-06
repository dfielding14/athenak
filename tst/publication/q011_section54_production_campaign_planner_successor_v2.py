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


def build_plan(resource_only_freeze_record: object | None = None) -> dict[str, object]:
    """Build the core or valid resource-frozen plan without execution authority."""
    stage = contract.build_stage_contract("planner")
    normalized_freeze = (
        None
        if resource_only_freeze_record is None
        else contract.validate_resource_only_freeze_record(resource_only_freeze_record)
    )
    selected_seeds = (
        contract.CORE_QUALIFYING_SEEDS
        if normalized_freeze is None
        else tuple(normalized_freeze["selected_qualifying_seeds"])
    )
    planned_triad_count = len(selected_seeds)
    campaign_size_selection = {
        "status": (
            "mandatory_core_source_local_plan_pending_resource_only_freeze"
            if normalized_freeze is None
            else "resource_only_triad_count_frozen_execution_still_blocked"
        ),
        "planned_triad_count": planned_triad_count,
        "frozen_triad_count": (
            None if normalized_freeze is None else normalized_freeze["frozen_triad_count"]
        ),
        "selected_qualifying_seeds": list(selected_seeds),
        "unselected_preregistered_reserve_seeds": list(
            contract.QUALIFYING_SEED_POOL[planned_triad_count:]
        ),
        "complete_paired_triad_required": True,
        "ordered_preregistered_prefix_required": True,
        "maximum_preregistered_attempt_count_is_not_promised": True,
        "resource_only_freeze_record_sha256": (
            None
            if normalized_freeze is None
            else contract.canonical_sha256(normalized_freeze)
        ),
        "reporting_rule": contract.reporting_rule_for_frozen_n(planned_triad_count),
        "post_qualifying_execution_rule": contract.campaign_size_resource_design()[
            "post_qualifying_execution_rule"
        ],
    }
    attempts: list[dict[str, object]] = []
    index = 0
    for variant in contract.GRID_VARIANTS:
        for seed in selected_seeds:
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
        "campaign_size_resource_design_sha256": contract.canonical_sha256(
            contract.campaign_size_resource_design()
        ),
        "resource_only_freeze_record": normalized_freeze,
        "campaign_size_selection": campaign_size_selection,
        "attempt_count": len(attempts),
        "attempts": attempts,
        "blockers": contract.source_local_blockers(normalized_freeze),
        "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
    }


def validate_plan(value: object) -> dict[str, object]:
    """Require the exact deterministic core or resource-frozen plan."""
    embedded_freeze = (
        value.get("resource_only_freeze_record") if type(value) is dict else None
    )
    expected = build_plan(embedded_freeze)
    _require(
        contract._strict_equal(value, expected),
        "Q011 production-campaign planner successor drifted or acquired authority",
    )
    return expected


def main() -> None:
    print(contract.canonical_json_bytes(build_plan()).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
