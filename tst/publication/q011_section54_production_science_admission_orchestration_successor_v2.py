#!/usr/bin/env python3
"""Fail-closed source-local admission successor for the Q-011 production path."""

from __future__ import annotations

from tst.publication import q011_section54_production_campaign_contract_successor_v2 as contract
from tst.publication import q011_section54_qualifying_campaign_execution_successor_v2 as execution


RECORD_TYPE = "q011_section54_source_local_admission_successor_v2"


class SourceLocalAdmissionSuccessorError(ValueError):
    """Reject a drifted or authority-bearing source-local admission record."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceLocalAdmissionSuccessorError(message)


def build_source_local_admission(handoffs: object | None = None) -> dict[str, object]:
    """Bind the complete handoff matrix while admitting no attempt as evidence."""
    normalized = execution.validate_execution_handoffs(
        execution.build_execution_handoffs() if handoffs is None else handoffs
    )
    stage = contract.build_stage_contract("admission")
    attempts = [
        {
            "attempt_id": handoff["attempt_id"],
            "variant": handoff["variant"],
            "qualifying_seed": handoff["qualifying_seed"],
            "input_deck": stage["production_deck"],
            "admitted": False,
            "raw_artifact_inventory_present": False,
            "registered_execution_receipt_present": False,
            "status": "blocked_pending_registered_execution_and_raw_artifact_bytes",
            "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
        }
        for handoff in normalized["handoffs"]
    ]
    return {
        "record_type": RECORD_TYPE,
        "schema_version": contract.SCHEMA_VERSION,
        "status": "source_local_admission_contract_complete_no_attempts_admitted",
        "stage_contract": stage,
        "execution_handoffs_sha256": contract.canonical_sha256(normalized),
        "runtime_source_closure": contract.runtime_source_closure(),
        "attempt_count": len(attempts),
        "admitted_attempt_count": 0,
        "attempts": attempts,
        "blockers": contract.source_local_blockers(),
        "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
    }


def validate_source_local_admission(
    value: object, handoffs: object | None = None
) -> dict[str, object]:
    expected = build_source_local_admission(handoffs)
    _require(
        contract._strict_equal(value, expected),
        "Q011 source-local admission successor drifted or acquired authority",
    )
    return expected


def main() -> None:
    print(contract.canonical_json_bytes(build_source_local_admission()).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
