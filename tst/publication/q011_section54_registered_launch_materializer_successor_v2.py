#!/usr/bin/env python3
"""Source-local registered-launch review successor for the Q-011 campaign."""

from __future__ import annotations

from tst.publication import q011_section54_production_campaign_contract_successor_v2 as contract
from tst.publication import q011_section54_production_campaign_planner_successor_v2 as planner
from tst.publication import (
    q011_resource_scaling_preproduction_pilot_materializer_successor_v2 as resource_pilots,
)


RECORD_TYPE = "q011_section54_registered_launch_review_candidates_successor_v2"


class RegisteredLaunchSuccessorError(ValueError):
    """Reject a drifted or authority-bearing registered-launch review record."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RegisteredLaunchSuccessorError(message)


def build_review_candidates(
    plan: object | None = None, pilot_materialization: object | None = None
) -> dict[str, object]:
    """Build launch-prohibited candidates from the exact successor plan."""
    normalized_plan = planner.validate_plan(planner.build_plan() if plan is None else plan)
    default_pilot, default_files = resource_pilots.build_materialization(normalized_plan)
    normalized_pilot, _ = resource_pilots.validate_materialization(
        default_pilot if pilot_materialization is None else pilot_materialization,
        default_files,
        normalized_plan,
    )
    stage = contract.build_stage_contract("registered_launch_materializer")
    candidates = [
        {
            "attempt_id": attempt["attempt_id"],
            "variant": attempt["variant"],
            "qualifying_seed": attempt["qualifying_seed"],
            "input_deck": stage["production_deck"],
            "argv": attempt["argv"],
            "registered_policy_slice_status": "not_materialized_not_authorized",
            "fresh_runtime_bindings_status": "absent",
            "status": "source_local_review_candidate_incomplete",
            "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
        }
        for attempt in normalized_plan["attempts"]
    ]
    return {
        "record_type": RECORD_TYPE,
        "schema_version": contract.SCHEMA_VERSION,
        "status": "source_local_review_candidates_complete_registration_blocked",
        "stage_contract": stage,
        "plan_sha256": contract.canonical_sha256(normalized_plan),
        "resource_pilot_materialization_sha256": contract.canonical_sha256(
            normalized_pilot
        ),
        "resource_only_freeze_record": normalized_plan["resource_only_freeze_record"],
        "campaign_size_selection": normalized_plan["campaign_size_selection"],
        "candidate_count": len(candidates),
        "candidates": candidates,
        "blockers": list(normalized_plan["blockers"]),
        "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
    }


def validate_review_candidates(
    value: object,
    plan: object | None = None,
    pilot_materialization: object | None = None,
) -> dict[str, object]:
    inferred_plan = plan
    if inferred_plan is None and type(value) is dict:
        inferred_plan = planner.build_plan(value.get("resource_only_freeze_record"))
    expected = build_review_candidates(inferred_plan, pilot_materialization)
    _require(
        contract._strict_equal(value, expected),
        "Q011 registered-launch successor drifted or acquired authority",
    )
    return expected


def main() -> None:
    print(contract.canonical_json_bytes(build_review_candidates()).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
