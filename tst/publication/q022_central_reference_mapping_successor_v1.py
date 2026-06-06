#!/usr/bin/env python3
"""Validation for the additive Q022 central-reference mapping successor."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_ROOT = REPO_ROOT / "tst/publication/readiness"

BELL_MAP_PATH = (
    READINESS_ROOT
    / "q022_xcmp_corrected_nonlinear_bell_reference_map_successor_v1_2026-06-06.json"
)
SHOCK_MAP_PATH = (
    READINESS_ROOT
    / "q022_xcmp_bai_sun_section54_shock_reference_map_successor_v1_2026-06-06.json"
)
LEDGER_PATH = (
    READINESS_ROOT
    / "q022_central_reference_discrepancy_open_input_ledger_successor_v1_2026-06-06.json"
)
READINESS_PATH = (
    READINESS_ROOT / "q022_central_reference_mapping_successor_v1_2026-06-06.json"
)

SUN_BAI_REFERENCE_ID = "sun_bai_2023_arxiv_2304.10568v1"
BAI_2015_REFERENCE_ID = "bai_2015_arxiv_1412.1087"
EXPECTED_CLOSURE = (
    "PPC*deposit_qscale*species_charge*v_CR/V_root_cell=J_CR/c=2*B0*k0"
)
EXPECTED_SUCCESSOR_ARTIFACTS = {
    "tst/publication/readiness/"
    "q022_xcmp_corrected_nonlinear_bell_reference_map_successor_v1_2026-06-06.json",
    "tst/publication/readiness/"
    "q022_xcmp_bai_sun_section54_shock_reference_map_successor_v1_2026-06-06.json",
    "tst/publication/readiness/"
    "q022_central_reference_discrepancy_open_input_ledger_successor_v1_2026-06-06.json",
    "tst/publication/q022_central_reference_mapping_successor_v1.py",
    "tst/publication/test_q022_central_reference_mapping_successor_v1.py",
    "tst/publication/readiness/q022_central_reference_mapping_successor_v1_2026-06-06.md",
}
EXPECTED_ATHENAK_CONTRACTS = {
    "tst/publication/readiness/q043_bell_current_normalization_supersession_2026-06-06.json",
    "tst/publication/readiness/"
    "q019_nonlinear_bell_volume_aware_campaign_successor_design_2026-06-06.json",
    "tst/publication/readiness/"
    "q011_section54_production_science_diagnostic_successor_v1_2026-06-06.json",
    "inputs/publication/"
    "pic_parallel_shock_section54_production_science_successor_v1_vl2_tsc.athinput",
    "tst/publication/readiness/"
    "q011_section54_pressure_selection_publication_candidate_successor_2026-06-05.json",
}


class ContractError(ValueError):
    """The bounded mapping successor drifted or was weakened."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    _require(isinstance(value, dict), f"{path}: expected JSON object")
    return value


def _validate_no_authority(record: dict[str, Any]) -> None:
    authority = record.get("authority")
    _require(isinstance(authority, dict) and authority, "authority inventory missing")
    _require(not any(authority.values()), "authority must remain false")


def _validate_reference_artifacts(record: dict[str, Any]) -> None:
    artifacts = record.get("reference_artifacts")
    _require(isinstance(artifacts, list) and artifacts, "reference artifacts missing")
    _require(
        {item["reference_id"] for item in artifacts} == set(record["reference_ids"]),
        "reference artifact inventory drifted",
    )
    for artifact in artifacts:
        path = Path(artifact["path"])
        _require(path.is_absolute(), "reference artifact path must be absolute")
        _require(path.is_file(), f"reference artifact missing: {path}")
        _require(sha256(path) == artifact["sha256"], f"reference PDF drifted: {path}")
        _require(
            artifact["artifact_role"].startswith("authoritative_locally_retained_"),
            "non-authoritative reference artifact admitted",
        )
        _require(
            isinstance(artifact["page_count"], int) and artifact["page_count"] > 0,
            "invalid PDF page count",
        )


def _validate_bound_artifacts(artifacts: list[dict[str, Any]]) -> None:
    for artifact in artifacts:
        path = Path(artifact["path"])
        if not path.is_absolute():
            path = REPO_ROOT / path
        _require(path.is_file(), f"bound artifact missing: {path}")
        _require(sha256(path) == artifact["sha256"], f"bound artifact drifted: {path}")


def _validate_provenance(record: dict[str, Any]) -> None:
    artifacts = {
        item["reference_id"]: item for item in record["reference_artifacts"]
    }
    entries = record.get("provenance_entries")
    _require(isinstance(entries, list) and entries, "page/equation provenance missing")
    ids = [entry["entry_id"] for entry in entries]
    _require(len(ids) == len(set(ids)), "duplicate provenance entry ID")
    for entry in entries:
        _require(entry["reference_id"] in artifacts, "unknown provenance reference")
        page = entry["physical_pdf_page"]
        _require(
            isinstance(page, int) and 1 <= page <= artifacts[entry["reference_id"]]["page_count"],
            "physical PDF page outside bound artifact",
        )
        for field in (
            "section",
            "equation",
            "mapping_class",
            "reference_statement",
            "athenak_mapping",
        ):
            _require(isinstance(entry[field], str) and entry[field], f"{field} missing")
        _require(
            entry["numeric_tolerance_supplied_by_reference"] is False,
            "reference tolerance was invented",
        )


def validate_common_map(record: dict[str, Any]) -> None:
    _require(record["schema_version"] == 2, "map schema drifted")
    _require(
        record["record_type"]
        == "q022_reference_specific_equation_normalization_map_successor",
        "map record type drifted",
    )
    _require(
        record["map_status"] == "bounded_map_recorded_external_review_open",
        "bounded map status drifted",
    )
    _require(
        record["qualification_effect"]
        == "source_local_reference_mapping_only_no_execution_no_claim_authority",
        "qualification boundary drifted",
    )
    _validate_no_authority(record)
    _validate_reference_artifacts(record)
    _validate_provenance(record)
    for field in (
        "matched_equations",
        "intentional_mismatches",
        "unit_map",
        "normalization_map",
        "parameter_overlap",
        "excluded_regimes",
        "source_records",
        "blocked_inputs",
    ):
        _require(bool(record[field]), f"{field} must be non-empty")
    _require(record["comparison_tolerances"] == [], "numeric tolerances were invented")
    _require(
        record["reviewer_disposition"] == "pending external review",
        "external review was falsely closed",
    )


def validate_bell_map(record: dict[str, Any]) -> None:
    validate_common_map(record)
    _require(
        record["comparison_id"] == "XCMP-CORRECTED-NONLINEAR-BELL-TARGET-PAPER",
        "Bell comparison ID drifted",
    )
    _require(
        set(record["reference_ids"]) == {SUN_BAI_REFERENCE_ID, BAI_2015_REFERENCE_ID},
        "Bell reference scope drifted",
    )
    normalization = record["normalization_map"]
    _require(
        normalization["sun_bai_reference_k0"] == "k0=j_CR/(2*B_g*c)",
        "Sun-Bai k0 normalization drifted",
    )
    _require(
        normalization["athenak_raw_prtcl_j"] == "prtcl_j=J_CR/c",
        "AthenaK raw current semantics drifted",
    )
    _require(
        normalization["athenak_corrected_general_closure"] == EXPECTED_CLOSURE,
        "corrected Bell closure drifted",
    )
    _require(
        normalization["artificial_C_multiplier_allowed"] is False
        and normalization["legacy_volume_blind_or_C_multiplied_normalization_allowed"]
        is False,
        "legacy Bell normalization was readmitted",
    )
    _require(
        record["parameter_overlap"]["standalone_nonlinear_reference_dataset_available"]
        is False,
        "target paper falsely acquired a standalone nonlinear dataset",
    )
    excluded = " ".join(record["excluded_regimes"])
    _require("CR-Hall" in excluded and "standalone nonlinear saturation" in excluded,
             "Bell excluded regimes drifted")


def validate_shock_map(record: dict[str, Any]) -> None:
    validate_common_map(record)
    _require(
        record["comparison_id"] == "XCMP-BAI-SUN-SECTION54-SHOCK",
        "shock comparison ID drifted",
    )
    _require(
        set(record["reference_ids"]) == {SUN_BAI_REFERENCE_ID, BAI_2015_REFERENCE_ID},
        "shock reference scope drifted",
    )
    overlap = record["parameter_overlap"]
    sun = overlap["exact_Sun_Bai_Section_5_4_overlap"]
    _require(sun["domain_c_over_omega_pi"] == [48000, 3120], "shock domain drifted")
    _require(sun["amr_cell_sizes_c_over_omega_pi"] == [12, 6, 3], "AMR ladder drifted")
    _require(sun["u0_over_U_A0"] == 30 and sun["M_A"] == 30, "Mach setup drifted")
    _require(sun["eta"] == "1e-3" and sun["C_over_U_A0"] == "1e4",
             "shock eta or C drifted")
    _require(sun["p_inj_over_m_u0"] == "sqrt(10)", "injection momentum drifted")
    _require(sun["remove_birth_time_before_Omega0_inverse"] == 45,
             "early-particle removal drifted")
    _require(sun["nominal_spectrum_times_Omega0_inverse"] == [500, 1200],
             "shock diagnostic times drifted")
    _require(
        overlap["not_explicit_in_Sun_Bai_Section_5_4"] == {"P0": 1, "T0": 1},
        "Sun-Bai pressure uncertainty was weakened",
    )
    _require(
        record["normalization_map"]["selected_pressure"]
        == "P0=T0=1_from_Bai_2015_not_explicitly_from_Sun_Bai_Section_5_4",
        "selected pressure provenance drifted",
    )
    _require(
        "existing Q011 engineering acceptance range" in " ".join(record["excluded_regimes"]),
        "engineering gate became a Q022 reference tolerance",
    )


REQUIRED_LEDGER_IDS = {
    "Q022-CENTRAL-BELL-OPEN-STANDALONE-NONLINEAR-DATASET",
    "Q022-CENTRAL-BELL-OPEN-EXACT-Q019-EXECUTION-CONTRACT",
    "Q022-CENTRAL-BELL-DISCREPANCY-HALL-SCOPE",
    "Q022-CENTRAL-BELL-OPEN-COMPARISON-TOLERANCES",
    "Q022-CENTRAL-SHOCK-OPEN-SUN-P0",
    "Q022-CENTRAL-SHOCK-DISCREPANCY-IDEAL-VS-DYNAMIC-SURFACE",
    "Q022-CENTRAL-SHOCK-DISCREPANCY-ETA-DOMAIN-HALL",
    "Q022-CENTRAL-SHOCK-OPEN-FIGURE-EXTRACTION",
    "Q022-CENTRAL-SHOCK-OPEN-COMPARISON-TOLERANCES",
    "Q022-CENTRAL-CROSSCUT-OPEN-EXTERNAL-REVIEW",
}


def validate_ledger(record: dict[str, Any]) -> None:
    _require(
        record["record_type"]
        == "q022_central_reference_discrepancy_open_input_ledger_successor",
        "ledger record type drifted",
    )
    _validate_no_authority(record)
    entries = record["entries"]
    ids = [entry["id"] for entry in entries]
    _require(len(ids) == len(set(ids)), "duplicate ledger ID")
    _require(set(ids) == REQUIRED_LEDGER_IDS, "ledger inventory drifted")
    for entry in entries:
        for field in ("comparison_id", "kind", "status", "finding", "required_resolution"):
            _require(isinstance(entry[field], str) and entry[field], f"{field} missing")
        _require(entry["blocks"], "ledger entry must block at least one overclaim")
    _require(record["numeric_tolerances_frozen"] is False, "tolerances falsely frozen")
    _require(
        record["qualifying_output_inspection_authorized"] is False,
        "qualifying output inspection falsely authorized",
    )


def validate_readiness(record: dict[str, Any]) -> None:
    _require(
        record["record_type"] == "q022_central_reference_mapping_successor_readiness",
        "readiness record type drifted",
    )
    _validate_no_authority(record)
    _require(
        record["central_comparisons"]
        == [
            "XCMP-CORRECTED-NONLINEAR-BELL-TARGET-PAPER",
            "XCMP-BAI-SUN-SECTION54-SHOCK",
        ],
        "central comparison scope drifted",
    )
    _require(
        record["reference_ids"] == [BAI_2015_REFERENCE_ID, SUN_BAI_REFERENCE_ID],
        "readiness reference scope drifted",
    )
    _require(
        {item["path"] for item in record["successor_artifacts"]}
        == EXPECTED_SUCCESSOR_ARTIFACTS,
        "successor artifact inventory drifted",
    )
    _require(
        {item["path"] for item in record["athenak_contract_bindings"]}
        == EXPECTED_ATHENAK_CONTRACTS,
        "AthenaK contract inventory drifted",
    )
    expected_historical = {
        str(path.relative_to(REPO_ROOT))
        for path in READINESS_ROOT.glob("q022_*_2026-05-30.json")
    } | {
        str(path.relative_to(REPO_ROOT))
        for path in (READINESS_ROOT / "schemas").glob("q022_*.schema.json")
    }
    _require(
        {item["path"] for item in record["historical_q022_bindings"]}
        == expected_historical,
        "historical Q022 inventory drifted",
    )
    _require(
        {item["reference_id"] for item in record["reference_pdf_bindings"]}
        == {BAI_2015_REFERENCE_ID, SUN_BAI_REFERENCE_ID},
        "readiness PDF inventory drifted",
    )
    for group in (
        "reference_pdf_bindings",
        "successor_artifacts",
        "historical_q022_bindings",
        "athenak_contract_bindings",
    ):
        _validate_bound_artifacts(record[group])
    _require(record["numeric_comparison_tolerances_frozen"] is False,
             "readiness falsely freezes numeric tolerances")
    _require(record["qualifying_output_inspection_authorized"] is False,
             "readiness falsely authorizes output inspection")
    _require(
        record["qualification_effect"]
        == "source_local_reference_mapping_successor_only_no_execution_no_claim_authority",
        "readiness qualification boundary drifted",
    )


def validate_all() -> None:
    validate_bell_map(load_json(BELL_MAP_PATH))
    validate_shock_map(load_json(SHOCK_MAP_PATH))
    validate_ledger(load_json(LEDGER_PATH))
    validate_readiness(load_json(READINESS_PATH))


if __name__ == "__main__":
    validate_all()
