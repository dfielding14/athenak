#!/usr/bin/env python3
"""Focused adversarial tests for the Q009 registered-pilot contract."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tempfile

import pytest

from tst.publication import q009_dynamic_amr_load_balance_registered_pilot_contract_v1 as contract


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _cases(value: dict[str, object]) -> dict[str, dict[str, object]]:
    return {case["case_id"]: case for case in value["pilot_matrix"]}


def _rejects(value: dict[str, object], match: str) -> None:
    with pytest.raises(contract.ContractError, match=match):
        contract.validate_contract(value, verify_repo_bindings=False)


@pytest.fixture(scope="module")
def frozen() -> dict[str, object]:
    return contract.load_contract()


def test_contract_is_exact_non_authorizing_five_case_frontier_hip_matrix(
    frozen: dict[str, object],
) -> None:
    assert frozen["record_type"] == contract.RECORD_TYPE
    assert frozen["contract_id"] == contract.CONTRACT_ID
    assert frozen["qualification_effect"] == contract.QUALIFICATION_EFFECT
    assert [case["case_id"] for case in frozen["pilot_matrix"]] == list(
        contract.EXPECTED_CASE_IDS
    )
    assert len(frozen["pilot_matrix"]) == 5
    assert all(case["resources"]["platform"] == "Frontier" for case in frozen["pilot_matrix"])
    assert all(case["resources"]["accelerator"] == "HIP" for case in frozen["pilot_matrix"])
    assert all(case["resources"]["mpi_ranks"] > 1 for case in frozen["pilot_matrix"])
    assert all(case["ceilings"]["maximum_attempts"] == 1 for case in frozen["pilot_matrix"])
    assert all(
        case["execution_boundary"]["launch_state"]
        == "blocked_pending_separate_admission"
        for case in frozen["pilot_matrix"]
    )
    assert all(value is False for value in frozen["authority_boundary"].values())


def test_exact_source_and_predecessor_bindings_are_current(
    frozen: dict[str, object],
) -> None:
    for group in ("predecessor_evidence", "source_templates"):
        for binding in frozen[group].values():
            assert _sha256(contract.REPO_ROOT / binding["path"]) == binding["sha256"]


def test_minimal_matrix_covers_required_dynamic_amr_load_balance_and_restart(
    frozen: dict[str, object],
) -> None:
    cases = _cases(frozen)
    migration = [
        cases["q009-amr-migration-lb-off-v1"],
        cases["q009-amr-migration-lb-on-v1"],
    ]
    restart = [
        cases["q009-shock-restart-post45-lb-off-v1"],
        cases["q009-shock-restart-post45-lb-on-v1"],
    ]
    assert {case["load_balance"]["pic_load_balance_cost_per_particle"] for case in migration} == {
        0.0,
        0.001,
    }
    assert {case["load_balance"]["pic_load_balance_cost_per_particle"] for case in restart} == {
        0.0,
        0.001,
    }
    assert all(
        case["dynamic_amr_profile_id"] == "repeated_refine_derefine" for case in migration
    )
    assert all(
        case["restart_profile_id"] == "checkpoint_before45_resume_after45"
        for case in restart
    )
    assert all(case["resources"]["scheduler_jobs_per_attempt"] == 2 for case in restart)
    assert all(
        case["resources"]["scheduler_jobs_per_attempt"] == 1
        for case in migration
    )
    profile = frozen["profile_definitions"]["restart_profiles"][
        "checkpoint_before45_resume_after45"
    ]
    assert profile["checkpoint_time_max_exclusive"] <= 45.0
    assert profile["resumed_terminal_time_min_exclusive"] == 45.0
    assert profile["resumed_terminal_time_max_inclusive"] > 45.0
    assert profile["restart_boundary_state_fingerprint_exact"]
    assert profile["restart_boundary_integer_counters_exact"]


def test_held_out_case_is_exact_full_geometry_and_cannot_be_retuned(
    frozen: dict[str, object],
) -> None:
    held = _cases(frozen)["q009-shock-heldout-full-geometry-lb-on-v1"]
    profile = frozen["profile_definitions"]["geometry_profiles"][held["geometry_profile_id"]]
    assert profile["held_out"]
    assert profile["exact_geometry"] == {
        "mesh/nx1": 4000,
        "mesh/nx2": 260,
        "mesh/x1max": 48000.0,
        "mesh/x2max": 3120.0,
        "meshblock/nx1": 20,
        "meshblock/nx2": 20,
    }
    assert profile["exact_overrides"] == {"time/nlim": 512}
    assert frozen["held_out_acceptance"]["same_contract_retuning_after_pair_telemetry"] == (
        "forbidden"
    )
    assert frozen["held_out_acceptance"]["allowed_postrun_read"] == (
        "aggregate_preregistered_telemetry_only"
    )


def test_budget_storage_queue_and_policy_boundaries_are_fail_closed(
    frozen: dict[str, object],
) -> None:
    envelope = frozen["budget_and_storage_envelope"]
    assert envelope["contract_hard_cap_node_hours"] == 100.0
    assert envelope["contract_hard_cap_artifact_storage_gib"] == 200.0
    assert sum(case["ceilings"]["maximum_node_hours"] for case in frozen["pilot_matrix"]) == 100.0
    assert (
        sum(case["ceilings"]["maximum_artifact_storage_gib"] for case in frozen["pilot_matrix"])
        == 200.0
    )
    assert envelope["overrun_action"] == "stop_fail_closed_no_retry"
    gates = set(frozen["registered_execution_gates"])
    assert "user_queue_empty_immediately_before_every_submission" in gates
    assert "separate_reviewed_registered_science_policy_slice_for_each_case" in gates
    assert "project_wide_node_hour_ledger_and_other_campaign_reserves_bound" in gates
    assert all(
        case["execution_boundary"]["policy_slice_state"]
        == "absent_not_authorized_by_this_contract"
        for case in frozen["pilot_matrix"]
    )


def test_uniform_exact_conservation_does_not_qualify_amr_and_telemetry_is_bounded(
    frozen: dict[str, object],
) -> None:
    boundary = frozen["scientific_boundary"]
    assert not boundary["uniform_mesh_exact_conservation_qualifies_dynamic_amr"]
    assert not boundary["paper_smooth_amr_individually_conservative_claimed"]
    assert "does not qualify AMR" in boundary["required_interpretation"]
    rules = frozen["required_telemetry"]["fail_closed_rules"]
    assert "missing_or_nonfinite_conservation_accounting_fails" in rules
    assert "unexplained_particle_loss_or_duplicate_count_must_equal_zero" in rules
    assert "no_scientific_conservation_tolerance_is_selected_by_this_pilot" in rules
    assert (
        frozen["inspection_policy"]["pre_preregistration_runtime_output_inspection"]
        == "forbidden"
    )
    assert "cell_or_particle_field_values" in frozen["inspection_policy"]["forbidden_inspection"]


def test_readiness_binds_all_new_contract_sources_with_exact_hashes(
    frozen: dict[str, object],
) -> None:
    readiness = json.loads(contract.READINESS_PATH.read_text(encoding="utf-8"))
    assert readiness["record_type"] == (
        "q009_dynamic_amr_load_balance_registered_pilot_contract_readiness"
    )
    assert readiness["status"] == "source_local_contract_validated_execution_blocked"
    assert readiness["contract_summary"]["case_count"] == len(frozen["pilot_matrix"])
    assert readiness["contract_summary"]["hard_cap_node_hours"] == 100.0
    assert readiness["contract_summary"]["hard_cap_artifact_storage_gib"] == 200.0
    assert readiness["source_constraints"]["new_files_only_under_tst_publication"]
    assert not readiness["authorization"]["launch_authorized"]
    assert not readiness["authorization"]["publication_authorized"]
    for binding in readiness["source_bindings"].values():
        assert _sha256(contract.REPO_ROOT / binding["path"]) == binding["sha256"]


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        (
            lambda value: value["pilot_matrix"].__delitem__(-1),
            "pilot matrix case ids or order drifted",
        ),
        (
            lambda value: value["pilot_matrix"][0]["resources"].update({"mpi_ranks": 8}),
            "summary drifted",
        ),
        (
            lambda value: value["pilot_matrix"][1]["load_balance"].update(
                {"pic_load_balance_cost_per_particle": 0.0}
            ),
            "summary drifted",
        ),
        (
            lambda value: value["pilot_matrix"][2].update({"restart_profile_id": "none"}),
            "summary drifted",
        ),
        (
            lambda value: value["profile_definitions"]["restart_profiles"][
                "checkpoint_before45_resume_after45"
            ].update({"resumed_terminal_time_min_exclusive": 44.0}),
            "post-45 restart profile drifted",
        ),
        (
            lambda value: value["profile_definitions"]["geometry_profiles"][
                "full_section54_geometry"
            ]["exact_geometry"].update({"mesh/nx2": 20}),
            "full held-out geometry drifted",
        ),
        (
            lambda value: value["budget_and_storage_envelope"].update(
                {"contract_hard_cap_node_hours": 101.0}
            ),
            "budget or storage envelope drifted",
        ),
        (
            lambda value: value["pilot_matrix"][4]["ceilings"].update(
                {"maximum_artifact_storage_gib": 256.0}
            ),
            "summary drifted",
        ),
        (
            lambda value: value["registered_execution_gates"].remove(
                "user_queue_empty_immediately_before_every_submission"
            ),
            "registered execution gates drifted",
        ),
        (
            lambda value: value["inspection_policy"]["post_preregistration_allowed_inspection"].append(
                "cell_or_particle_field_values"
            ),
            "inspection policy drifted",
        ),
        (
            lambda value: value["authority_boundary"].update({"publication_authorized": True}),
            "authority must remain false",
        ),
        (
            lambda value: value["scientific_boundary"].update(
                {"uniform_mesh_exact_conservation_qualifies_dynamic_amr": True}
            ),
            "scientific boundary drifted",
        ),
    ],
)
def test_adversarial_semantic_drift_fails_closed(
    frozen: dict[str, object], mutation, match: str
) -> None:
    value = deepcopy(frozen)
    result = mutation(value)
    if result is not None:
        # Mutation lambdas normally mutate in place; this catches accidental no-op replacements.
        raise AssertionError(f"mutation unexpectedly returned {result!r}")
    _rejects(value, match)


def test_drifted_bound_source_bytes_fail_closed(frozen: dict[str, object]) -> None:
    value = deepcopy(frozen)
    value["source_templates"]["q009_coupled_boundary_lifetime"]["sha256"] = "0" * 64
    with pytest.raises(contract.ContractError, match="SHA-256 drifted"):
        contract.validate_contract(value, verify_repo_bindings=True)


def test_duplicate_key_nonfinite_json_and_symlink_contract_fail_closed() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        duplicate = root / "duplicate.json"
        duplicate.write_text('{"record_type":"x","record_type":"y"}', encoding="utf-8")
        with pytest.raises(contract.ContractError, match="duplicate key"):
            contract.load_contract(duplicate)

        nonfinite = root / "nonfinite.json"
        nonfinite.write_text('{"value":NaN}', encoding="utf-8")
        with pytest.raises(contract.ContractError, match="forbidden JSON constant"):
            contract.load_contract(nonfinite)

        alias = root / "alias.json"
        alias.symlink_to(contract.CONTRACT_PATH)
        with pytest.raises(contract.ContractError, match="symlink or path alias forbidden"):
            contract.load_contract(alias)
