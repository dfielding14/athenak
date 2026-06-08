#!/usr/bin/env python3
"""Focused tests for Q019 excluded-pilot engineering qualification."""

from __future__ import annotations

import copy
from pathlib import Path

from tst.publication import q019_excluded_pilot_campaign_driver_v1 as campaign
from tst.publication import q019_excluded_pilot_engineering_qualification_v1 as qual
from tst.publication.test_q019_excluded_pilot_campaign_driver_v1 import _attempts


INDEX_PATH = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/PIC/analysis/"
    "q019_nonlinear_bell_registered_successor_v1/"
    "q019_excluded_pilot_execution_index_v1.json"
)


def _build(attempts: list[dict[str, object]]) -> dict[str, object]:
    index = campaign.build_execution_index(attempts)
    return qual.build_qualification(
        index,
        execution_index_path=INDEX_PATH,
        execution_index_sha256="a" * 64,
    )


def test_fast_pairs_pass_engineering_only_without_science_authority() -> None:
    record = _build(_attempts())
    assert record["status"] == qual.STATUS_PASS
    assert record["decision"]["engineering_gate_pass"] is True
    assert record["decision"]["production_resource_freeze_recommended"] is True
    assert record["decision"]["production_resource_freeze_authorized"] is False
    assert record["saturation_evidence_eligible"] is False
    assert all(value is False for value in record["authorization"].values())
    assert all(
        item["pair_engineering_gate_pass"]
        for item in record["pair_qualifications"]
    )


def test_large_overhead_requires_redesign_even_with_valid_execution() -> None:
    attempts = _attempts()
    slow = copy.deepcopy(attempts)
    event = slow[1]["reconciliation_event"]
    event["elapsed_seconds"] = 150
    event["consumed_node_hours"] = event["billed_nodes"] * 150 / 3600.0
    event["event_sha256"] = campaign._ledger_event_sha256(event)
    slow[1]["admission_facts"]["registered_execution_identity"][
        "reconciliation_event_sha256"
    ] = event["event_sha256"]
    record = _build(slow)
    assert record["status"] == qual.STATUS_FAIL
    assert record["decision"]["engineering_gate_pass"] is False
    assert record["decision"]["selected_runtime_box_edge_monitor_dt"] is None
    assert record["pair_qualifications"][0]["overhead_gate_pass"] is False


def test_resource_headroom_is_required_independently_of_overhead() -> None:
    attempts = _attempts()
    near_limit = copy.deepcopy(attempts)
    for index in (2, 3):
        event = near_limit[index]["reconciliation_event"]
        event["elapsed_seconds"] = 3000
        event["consumed_node_hours"] = event["billed_nodes"] * 3000 / 3600.0
        event["event_sha256"] = campaign._ledger_event_sha256(event)
        near_limit[index]["admission_facts"]["registered_execution_identity"][
            "reconciliation_event_sha256"
        ] = event["event_sha256"]
    record = _build(near_limit)
    assert record["pair_qualifications"][1]["overhead_gate_pass"] is True
    assert record["pair_qualifications"][1]["resource_headroom_gate_pass"] is False
    assert record["decision"]["engineering_gate_pass"] is False
