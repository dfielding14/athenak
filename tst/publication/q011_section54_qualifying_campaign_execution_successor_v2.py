#!/usr/bin/env python3
"""Launch-prohibited Q-011 production campaign execution-handoff successor."""

from __future__ import annotations

from tst.publication import q011_section54_production_campaign_contract_successor_v2 as contract
from tst.publication import q011_section54_production_campaign_planner_successor_v2 as planner
from tst.publication import q011_section54_registered_launch_materializer_successor_v2 as launch


RECORD_TYPE = "q011_section54_execution_handoffs_successor_v2"


class CampaignExecutionSuccessorError(ValueError):
    """Reject a drifted or authority-bearing source-local execution handoff."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise CampaignExecutionSuccessorError(message)


def build_execution_handoffs(
    plan: object | None = None, review_candidates: object | None = None
) -> dict[str, object]:
    """Build deterministic blocked handoffs for all core attempts."""
    normalized_plan = planner.validate_plan(planner.build_plan() if plan is None else plan)
    normalized_candidates = launch.validate_review_candidates(
        launch.build_review_candidates(normalized_plan)
        if review_candidates is None
        else review_candidates,
        normalized_plan,
    )
    stage = contract.build_stage_contract("campaign_execution")
    handoffs = [
        {
            "attempt_id": candidate["attempt_id"],
            "variant": candidate["variant"],
            "qualifying_seed": candidate["qualifying_seed"],
            "input_deck": stage["production_deck"],
            "argv": candidate["argv"],
            "registered_launch_candidate_sha256": contract.canonical_sha256(candidate),
            "raw_root_template": f"<registered-orion-raw-root>/{candidate['attempt_id']}",
            "execution_status": "blocked_no_registered_policy_or_fresh_runtime_bindings",
            "registered_execution_receipt": None,
            "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
        }
        for candidate in normalized_candidates["candidates"]
    ]
    return {
        "record_type": RECORD_TYPE,
        "schema_version": contract.SCHEMA_VERSION,
        "status": "source_local_handoffs_complete_execution_blocked",
        "stage_contract": stage,
        "plan_sha256": contract.canonical_sha256(normalized_plan),
        "registered_launch_review_sha256": contract.canonical_sha256(
            normalized_candidates
        ),
        "handoff_count": len(handoffs),
        "handoffs": handoffs,
        "blockers": contract.source_local_blockers(),
        "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
    }


def validate_execution_handoffs(
    value: object, plan: object | None = None, review_candidates: object | None = None
) -> dict[str, object]:
    expected = build_execution_handoffs(plan, review_candidates)
    _require(
        contract._strict_equal(value, expected),
        "Q011 execution-handoff successor drifted or acquired authority",
    )
    return expected


def main() -> None:
    print(contract.canonical_json_bytes(build_execution_handoffs()).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
