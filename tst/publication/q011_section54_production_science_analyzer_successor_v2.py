#!/usr/bin/env python3
"""Fail-closed Q-011 production-science analyzer contract successor."""

from __future__ import annotations

from tst.publication import q011_section54_production_campaign_contract_successor_v2 as contract
from tst.publication import (
    q011_section54_production_science_admission_orchestration_successor_v2 as admission,
)


RECORD_TYPE = "q011_section54_source_local_analysis_gate_packet_successor_v2"


class AnalyzerContractSuccessorError(ValueError):
    """Reject a drifted or authority-bearing source-local analysis packet."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AnalyzerContractSuccessorError(message)


def build_source_local_analysis_packet(
    admission_record: object | None = None,
) -> dict[str, object]:
    """Bind all frozen gates while inspecting no qualifying output."""
    normalized = admission.validate_source_local_admission(
        admission.build_source_local_admission()
        if admission_record is None
        else admission_record
    )
    stage = contract.build_stage_contract("analyzer")
    gates = contract.validate_source_local_gate_records(contract.source_local_gate_records())
    return {
        "record_type": RECORD_TYPE,
        "schema_version": contract.SCHEMA_VERSION,
        "status": "source_local_gate_packet_complete_all_gates_blocked",
        "stage_contract": stage,
        "admission_sha256": contract.canonical_sha256(normalized),
        "production_deck": stage["production_deck"],
        "shock_front_detector": "unique_strongest_negative_density_gradient",
        "gate_count": len(gates),
        "passed_gate_count": 0,
        "measurements_inspected": False,
        "gates": gates,
        "blockers": contract.source_local_blockers(),
        "authorization": dict(contract.AUTHORIZATION_BOUNDARY),
    }


def validate_source_local_analysis_packet(
    value: object, admission_record: object | None = None
) -> dict[str, object]:
    expected = build_source_local_analysis_packet(admission_record)
    _require(
        contract._strict_equal(value, expected),
        "Q011 analyzer successor drifted or acquired authority/results",
    )
    return expected


def main() -> None:
    print(contract.canonical_json_bytes(build_source_local_analysis_packet()).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
