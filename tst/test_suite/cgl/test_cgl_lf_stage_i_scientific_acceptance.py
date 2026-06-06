"""Focused adversarial tests for standalone Stage I scientific acceptance."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
UTILITY = REPOSITORY / "scripts/frontier/cgl_lf_stage_i_scientific_acceptance.py"


def load_utility():
    spec = importlib.util.spec_from_file_location(
        "cgl_lf_stage_i_scientific_acceptance", UTILITY
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


acceptance = load_utility()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_history(path: Path, columns: dict[str, list[float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    labels = list(columns)
    rows = zip(*(columns[label] for label in labels))
    text = "# " + " ".join(
        f"[{index}]={label}" for index, label in enumerate(labels, start=1)
    ) + "\n"
    text += "\n".join(
        " ".join(format(value, ".17g") for value in row) for row in rows
    )
    path.write_text(text + "\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def policy():
    return acceptance.load_validated_policy(
        acceptance.DEFAULT_CRITERIA, acceptance.DEFAULT_CRITERIA_REVIEW
    )


@pytest.fixture
def fast_policy(policy):
    value = deepcopy(policy)
    value["criteria"]["statistics_policy"]["bootstrap_replicates"] = 80
    value["criteria"]["statistics_policy"]["gap_policy"]["expected_history_cadence"] = 0.25
    return value


def histories(root: Path, *, passive: bool = False) -> tuple[Path, Path]:
    times = [3.5 + 0.25 * index for index in range(27)]
    count = len(times)
    mhd = {
        "time": times,
        "lf_hwproj": [0.0] * count,
        "lf_cpwrk": [0.0 if passive else 1.0e-4 * index for index in range(count)],
        "lf_cawrk": [0.0 if passive else -2.0e-4 * index for index in range(count)],
        "lf_qface": [100.0 * index for index in range(count)],
        "lf_qprwrk": [1.0e-3 * index for index in range(count)],
        "lf_qpewrk": [-2.0e-3 * index for index in range(count)],
    }
    user = {
        "time": times,
        "kinetic": [1.0] * count,
        "magnetic": [1.5] * count,
        "beta": [100.0] * count,
        "abs_dp": [0.2] * count,
        "mirror_vol": [1.0e-3] * count,
        "fire_vol": [2.0e-3] * count,
        "hard_vol": [0.0] * count,
        "nu_eff": [20.0] * count,
        "force_pwr": [0.32] * count,
        "force_prp2": [1.0] * count,
        "force_prl2": [0.0] * count,
    }
    mhd_path = root / "case.mhd.hst"
    user_path = root / "case.user.hst"
    write_history(mhd_path, mhd)
    write_history(user_path, user)
    return mhd_path, user_path


def build_case_bundle(
    policy,
    root: Path,
    case_id: str,
    mhd_path: Path,
    user_path: Path,
) -> Path:
    case = acceptance.case_manifest_record(policy, case_id)
    executable = "a" * 64
    input_sha = "b" * 64
    segment_paths = [root / "segments" / "s00.json", root / "segments" / "s01.json"]
    for index, (path, final_time) in enumerate(zip(segment_paths, (0.1, 10.0))):
        parent = None
        if index:
            parent = {
                "manifest": str(segment_paths[index - 1].resolve()),
                "case_id": case_id,
                "result": "accepted",
                "final_time": 0.1,
                "executable_sha256": executable,
                "input_sha256": input_sha,
            }
        write_json(path, {
            "accounting": {
                "result": "accepted",
                "case_id": case_id,
                "case_name": case["name"],
                "segment": f"s0{index}_synthetic",
                "executable_sha256": executable,
            },
            "command": {
                "executable_sha256": executable,
                "input_sha256": input_sha,
                "matrix_sha256": policy["verified_sources"]["stage_i_manifest"]["sha256"],
                "parent_segment": parent,
                "restart_file": None if index == 0 else "synthetic-parent.rst",
                "restart_files": [] if index == 0 else ["synthetic-parent.rst"],
            },
            "scientific_inspection": {
                "accepted": True,
                "case_id": case_id,
                "manifest": str(path.resolve()),
                "final_time": final_time,
                "checks": {
                    "terminal_restart_physical_time_matches_final": index == 1,
                },
                "terminal_restart_time": final_time,
            },
        })
    bundle_path = root / "bundle.json"
    write_json(bundle_path, {
        "workflow": "paper-mks24-stage-i-production",
        "status": "accepted_for_analysis",
        "production_case_id": case_id,
        "required_final_time": 10.0,
        "accepted_final_time": 10.0,
        "cases": [{
            "name": case["name"],
            "input": case["input"],
            "status": "passed",
            "outputs": {
                "mhd_history": str(mhd_path.resolve().relative_to(root.resolve())),
                "user_history": str(user_path.resolve().relative_to(root.resolve())),
            },
            "model_choices": {"forcing_tcorr": "2.0"},
        }],
        "production_segment_manifests": [
            str(path.resolve()) for path in segment_paths
        ],
    })
    return bundle_path


def test_preregistered_criteria_bind_final_utility_and_completed_reviews(policy):
    assert policy["review_status"] == "approved"
    assert policy["approved"] is True
    assert policy["review"]["reviews"] == [
        {
            "role": "plasma_physics",
            "reviewer_id": "019e9a8d-253b-7010-b156-676866801f3c",
            "decision": "approved",
            "independent_of_implementation": True,
        },
        {
            "role": "statistical_methodology",
            "reviewer_id": "019e9ae1-d8c6-7861-9bbb-69e2cc97ba6f",
            "decision": "approved",
            "independent_of_implementation": True,
        },
    ]
    assert policy["review"]["method_revision"]["candidate_status"] == (
        "approved_by_independent_plasma_and_statistical_review"
    )
    assert policy["review"]["remaining_review_requirements"] == []
    utility_sha = acceptance.regular_file_binding(UTILITY, "utility")["sha256"]
    assert policy["criteria"]["source_bindings"]["acceptance_utility"]["sha256"] == utility_sha
    assert policy["review"]["acceptance_utility"]["sha256"] == utility_sha
    assert policy["review"]["criteria"]["sha256"] == policy["criteria_binding"]["sha256"]
    assert policy["criteria"]["family_gates"]["lf_strength"]["cases"] == [
        "R12", "R02", "R06", "R13"
    ]
    assert policy["criteria"]["statistics_policy"]["minimum_independent_time_blocks"] == {
        "full": 3.0,
        "comparison": 2.0,
    }
    assert policy["criteria"]["analysis_windows"] == {
        "full": [4.0, 10.0],
        "early": [4.0, 8.0],
        "late": [6.0, 10.0],
    }
    change = policy["criteria"]["criteria_change_record"]
    assert change["change_id"] == "physical-time-stationarity-r03-r17-v3"
    assert change["previous_policy"]["minimum_independent_time_blocks"] == {
        "full": 3.0,
        "half": 1.5,
    }
    assert change["review_disposition"] == (
        "approved_by_independent_plasma_and_statistical_review"
    )
    assert policy["criteria"]["extension_policy"]["forbidden_after_results"] == [
        "relax thresholds",
        "move the analysis windows",
        "move the convergence interval",
        "drop failed admitted products",
        "reclassify failed or inconclusive gates as passed",
        "use any family, panel, convergence, CT, or product gate as an extension trigger",
        "use a passed or failed stationarity gate as an extension trigger",
    ]
    assert policy["criteria"]["extension_policy"]["current_policy_authorizes_extension"] is False
    assert policy["criteria"]["extension_policy"]["prospective_extension_rule"][
        "analysis_windows"
    ] == {
        "full": [6.0, 12.0],
        "early": [6.0, 10.0],
        "late": [8.0, 12.0],
    }
    assert policy["criteria"]["extension_policy"]["prospective_extension_rule"][
        "maximum_extensions"
    ] == 1
    assert policy["criteria"]["extension_policy"]["prospective_extension_rule"][
        "fixed_combination_rule"
    ]["eligible_t10_gate_result"] == "inconclusive"
    assert policy["criteria"]["statistics_policy"]["stationarity_contrast"] == (
        "paired-physical-time-early-minus-late-effect-size-with-descriptive-bootstrap"
    )
    assert "median_cadence_upper_relative_tolerance" not in policy["criteria"][
        "statistics_policy"
    ]["gap_policy"]
    assert policy["criteria"]["scientific_products_policy"]["reviewed_generator_binding"][
        "status"
    ] == "unavailable_pending_companion_generator"
    assert "analyzer_contract" not in policy["criteria"]["source_bindings"]
    assert policy["verified_sources"]["current_source_authority_evidence"]["sha256"] == (
        "cb50beb064678a9446ac33801a0023547d06c432d8c59b6bd3fbe34b11cf0391"
    )
    assert policy["review"]["retained_ct_observation"][
        "approval_identity_available"
    ] is False
    evidence = acceptance.validate_criteria_evidence(policy)
    acceptance.verify_evidence_digest(evidence, "criteria validation")
    assert evidence["valid"] is True
    assert evidence["independent_review_complete"] is True


def test_approved_review_schema_requires_exact_completed_reviews(policy):
    review = deepcopy(policy["review"])
    review["reviews"][0]["independent_of_implementation"] = False
    with pytest.raises(acceptance.AcceptanceError, match="reviewer records are incoherent"):
        acceptance.validate_criteria_review(
            review,
            policy["review_binding"],
            policy["criteria"],
            policy["criteria_binding"],
            policy["verified_sources"]["acceptance_utility"],
        )


def test_approved_review_schema_rejects_nonminimal_reviewer_record(policy):
    review = deepcopy(policy["review"])
    review["reviews"][0]["review_note"] = "not part of the approved minimal record"
    with pytest.raises(acceptance.AcceptanceError, match="reviewer records are incoherent"):
        acceptance.validate_criteria_review(
            review,
            policy["review_binding"],
            policy["criteria"],
            policy["criteria_binding"],
            policy["verified_sources"]["acceptance_utility"],
        )


@pytest.mark.parametrize(
    ("field", "value", "match"),
    [
        ("review_status", "changes_required", "criteria disposition are incoherent"),
        (
            "candidate_status",
            "ready_for_independent_plasma_and_statistical_review",
            "method revision differs",
        ),
        (
            "remaining_review_requirements",
            ["review still required"],
            "remaining requirements are incoherent",
        ),
        (
            "required_review_roles",
            ["plasma_physics"],
            "required roles differ",
        ),
    ],
)
def test_approved_review_rejects_incoherent_status_dependent_state(
    policy, field, value, match
):
    review = deepcopy(policy["review"])
    if field == "candidate_status":
        review["method_revision"][field] = value
    else:
        review[field] = value
    with pytest.raises(acceptance.AcceptanceError, match=match):
        acceptance.validate_criteria_review(
            review,
            policy["review_binding"],
            policy["criteria"],
            policy["criteria_binding"],
            policy["verified_sources"]["acceptance_utility"],
        )


def test_endpoint_clipped_trapezoidal_statistics_are_deterministic():
    times = [7.5, 8.5, 9.5, 10.5]
    values = list(times)
    first = acceptance.window_statistics(
        times, values, 8.0, 10.0, replicates=100, seed_text="linear"
    )
    second = acceptance.window_statistics(
        times, values, 8.0, 10.0, replicates=100, seed_text="linear"
    )
    assert first == second
    assert first["mean"] == pytest.approx(9.0)
    assert first["sample_count"] == 4


def test_time_block_bootstrap_resamples_coherent_weighted_intervals():
    times = [0.0, 1.0, 9.0, 10.0]
    values = [0.0, 0.0, 10.0, 10.0]
    samples = acceptance.bootstrap_means(
        times, values, replicates=40, seed_text="full-time-block", block_duration=10.0
    )
    assert samples == pytest.approx([5.0] * 40)


def test_bootstrap_block_duration_honors_forcing_tcorr_floor():
    result = acceptance.window_statistics(
        [0.0, 0.5, 1.0, 1.5, 2.0],
        [0.0, 1.0, 0.0, 1.0, 0.0],
        0.0,
        2.0,
        replicates=40,
        seed_text="tcorr-floor",
        minimum_block_duration=1.5,
    )
    assert result["method"]["bootstrap_block_duration"] == pytest.approx(1.5)
    assert result["method"]["minimum_block_duration"] == pytest.approx(1.5)
    assert result["effective_sample_count"] <= 2.0 / 1.5 + 1.0e-12


def test_gap_aware_sampling_rejects_sparse_physical_time_and_caps_independence():
    dense = acceptance.window_statistics(
        [0.25 * index for index in range(25)],
        [(-1.0) ** index for index in range(25)],
        0.0,
        6.0,
        replicates=40,
        seed_text="dense-physical-cap",
        minimum_block_duration=2.0,
        expected_cadence=0.25,
    )
    assert dense["gap_adequacy"] == "pass"
    assert dense["physical_independence_cap"] == pytest.approx(3.0)
    assert dense["effective_sample_count"] <= 3.0

    sparse = acceptance.window_statistics(
        [0.0, 0.25, 0.5, 3.0, 3.25, 3.5, 6.0],
        [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
        0.0,
        6.0,
        replicates=40,
        seed_text="sparse-gap",
        minimum_block_duration=2.0,
        expected_cadence=0.25,
    )
    assert sparse["gap_adequacy"] == "inconclusive"
    assert sparse["maximum_gap"] > sparse["maximum_allowed_gap"]


def test_gap_adequacy_is_invariant_to_redundant_clustered_samples():
    policy = {
        "replicates": 40,
        "minimum_block_duration": 2.0,
        "expected_cadence": 0.4,
    }
    accepted_times = [float(index) for index in range(7)]
    accepted_values = [math.sin(time) for time in accepted_times]
    accepted_clustered_times = sorted({
        *accepted_times,
        *(2.0 + 0.001 * index for index in range(1, 500)),
    })
    accepted_clustered_values = [
        acceptance.interpolate_at(accepted_times, accepted_values, time)
        for time in accepted_clustered_times
    ]
    accepted = acceptance.window_statistics(
        accepted_times,
        accepted_values,
        0.0,
        6.0,
        seed_text="accepted-gap-base",
        **policy,
    )
    accepted_clustered = acceptance.window_statistics(
        accepted_clustered_times,
        accepted_clustered_values,
        0.0,
        6.0,
        seed_text="accepted-gap-clustered",
        **policy,
    )
    assert accepted["gap_adequacy"] == accepted_clustered["gap_adequacy"] == "pass"
    assert accepted["maximum_gap"] == accepted_clustered["maximum_gap"]
    assert accepted["physical_time_coverage_fraction"] == pytest.approx(
        accepted_clustered["physical_time_coverage_fraction"]
    )
    assert accepted["descriptive_median_cadence"] != accepted_clustered[
        "descriptive_median_cadence"
    ]

    sparse_times = [0.0, 0.5, 3.0, 3.5, 6.0]
    sparse_values = [math.cos(time) for time in sparse_times]
    sparse_clustered_times = sorted({
        *sparse_times,
        *(0.001 * index for index in range(1, 500)),
    })
    sparse_clustered_values = [
        acceptance.interpolate_at(sparse_times, sparse_values, time)
        for time in sparse_clustered_times
    ]
    sparse = acceptance.window_statistics(
        sparse_times,
        sparse_values,
        0.0,
        6.0,
        seed_text="sparse-gap-base",
        **policy,
    )
    sparse_clustered = acceptance.window_statistics(
        sparse_clustered_times,
        sparse_clustered_values,
        0.0,
        6.0,
        seed_text="sparse-gap-clustered",
        **policy,
    )
    assert sparse["gap_adequacy"] == sparse_clustered["gap_adequacy"] == "inconclusive"
    assert sparse["maximum_gap"] == sparse_clustered["maximum_gap"]
    assert sparse["physical_time_coverage_fraction"] == pytest.approx(
        sparse_clustered["physical_time_coverage_fraction"]
    )


def test_physical_time_effective_sample_estimator_detects_correlation():
    times = [0.25 * index for index in range(20)]
    independent = [(-1.0) ** index for index in range(20)]
    correlated = [0.0] * 10 + [1.0] * 10
    independent_neff, _ = acceptance.effective_sample_size(times, independent, 0.25)
    correlated_neff, _ = acceptance.effective_sample_size(times, correlated, 0.25)
    assert independent_neff > correlated_neff


def test_physical_time_ess_treats_roundoff_scale_constant_series_as_constant():
    times = [0.25 * index for index in range(25)]
    values = [0.2] * len(times)
    neff, duration = acceptance.effective_sample_size(times, values, 0.25)
    assert duration == pytest.approx(0.25)
    assert neff == pytest.approx(24.0)


def test_physical_time_ess_is_invariant_to_adversarial_clustered_sampling():
    base_times = [float(index) for index in range(7)]
    base_values = [0.0, 1.0, 0.5, -0.5, -1.0, 0.0, 0.5]
    clustered_times = sorted({
        *base_times,
        *(0.001 * index for index in range(1, 500)),
        *(4.0 + 0.001 * index for index in range(1, 500)),
    })
    clustered_values = [
        acceptance.interpolate_at(base_times, base_values, time)
        for time in clustered_times
    ]

    base_neff, base_duration = acceptance.effective_sample_size(
        base_times, base_values, 0.25
    )
    clustered_neff, clustered_duration = acceptance.effective_sample_size(
        clustered_times, clustered_values, 0.25
    )

    assert clustered_neff == pytest.approx(base_neff, rel=1.0e-11, abs=1.0e-11)
    assert clustered_duration == pytest.approx(
        base_duration, rel=1.0e-11, abs=1.0e-11
    )


def test_paired_physical_time_contrast_is_invariant_to_clustered_sampling():
    base_times = [float(index) for index in range(11)]
    base_values = [math.sin(0.7 * time) for time in base_times]
    clustered_times = sorted({
        *base_times,
        *(4.0 + 0.002 * index for index in range(1, 500)),
        *(6.0 + 0.002 * index for index in range(1, 500)),
    })
    clustered_values = [
        acceptance.interpolate_at(base_times, base_values, time)
        for time in clustered_times
    ]
    arguments = {
        "early_start": 4.0,
        "early_end": 8.0,
        "late_start": 6.0,
        "late_end": 10.0,
        "replicates": 120,
        "seed_text": "paired-cluster-invariance",
        "block_duration": 2.0,
    }

    base = acceptance.paired_window_contrast(base_times, base_values, **arguments)
    clustered = acceptance.paired_window_contrast(
        clustered_times, clustered_values, **arguments
    )

    assert clustered["signed_early_minus_late"] == pytest.approx(
        base["signed_early_minus_late"], rel=1.0e-12, abs=1.0e-12
    )
    assert clustered["standard_error"] == pytest.approx(
        base["standard_error"], rel=1.0e-11, abs=1.0e-11
    )


def test_stationarity_fails_large_resolved_half_window_drift():
    policy = {
        "decision_authority": (
            "paired physical-time early-minus-late effect size only; moving-block "
            "bootstrap uncertainty is descriptive and never pass/fail authority"
        ),
        "scalar_relative_change_lte": 0.25,
        "occupancy_absolute_change_lte": 0.002,
        "forcing_power_relative_change_lte": 0.1,
    }
    result = acceptance.stationarity_result(
        {"mean": 1.0},
        {"mean": 0.5, "standard_error": 0.01},
        {"mean": 1.5, "standard_error": 0.01},
        {
            "signed_early_minus_late": -1.0,
            "standard_error": 0.01,
            "confidence_interval_95": [-1.02, -0.98],
            "method": {"contrast": "fixture-paired"},
        },
        "scalar",
        policy,
    )
    assert result["result"] == "fail"
    assert result["descriptive_bootstrap"]["z_score"] > 3.0
    assert result["descriptive_bootstrap"]["inferential_authority"] is False


def test_forcing_power_stationarity_ignores_descriptive_bootstrap_significance():
    policy = {
        "decision_authority": (
            "paired physical-time early-minus-late effect size only; moving-block "
            "bootstrap uncertainty is descriptive and never pass/fail authority"
        ),
        "scalar_relative_change_lte": 0.25,
        "occupancy_absolute_change_lte": 0.002,
        "forcing_power_relative_change_lte": 0.1,
    }
    result = acceptance.stationarity_result(
        {"mean": 1.0},
        {"mean": 1.0},
        {"mean": 1.05},
        {
            "signed_early_minus_late": -0.05,
            "standard_error": 0.001,
            "confidence_interval_95": [-0.052, -0.048],
            "method": {"contrast": "fixture-paired"},
        },
        "forcing_power",
        policy,
    )
    assert result["relative_change"] < 0.1
    assert result["descriptive_bootstrap"]["z_score"] > 3.0
    assert result["result"] == "pass"

    failed = acceptance.stationarity_result(
        {"mean": 1.0},
        {"mean": 1.0},
        {"mean": 1.2},
        {
            "signed_early_minus_late": -0.2,
            "standard_error": 10.0,
            "confidence_interval_95": [-20.0, 20.0],
            "method": {"contrast": "fixture-paired"},
        },
        "forcing_power",
        policy,
    )
    assert failed["descriptive_bootstrap"]["z_score"] < 3.0
    assert failed["result"] == "fail"


def test_extension_nomination_uses_only_frozen_stationarity_triggers(policy):
    passing = acceptance.gate(
        "stationarity:kinetic",
        "pass",
        reason="fixture",
        observations={
            "sampling_adequacy": "pass",
            "stationarity": {"result": "pass"},
        },
    )
    failing = acceptance.gate(
        "stationarity:force_power",
        "fail",
        reason="fixture",
        observations={
            "sampling_adequacy": "pass",
            "stationarity": {"result": "fail"},
        },
    )
    inconclusive = acceptance.gate(
        "stationarity:magnetic",
        "inconclusive",
        reason="fixture",
        observations={
            "sampling_adequacy": "inconclusive",
            "stationarity": {"result": "pass"},
        },
    )
    assessment = acceptance.extension_assessment(
        policy["criteria"], [passing, failing, inconclusive], "inconclusive"
    )
    assert assessment["triggered"] is True
    assert assessment["triggered_by"] == [
        {
            "gate": "stationarity:magnetic",
            "trigger": "sampling_adequacy=inconclusive",
        }
    ]
    assert assessment["current_policy_authorizes_extension"] is False
    assert assessment["preserved_t10_result"] == "inconclusive"
    assert assessment["prospective_extension_rule"]["nominated_endpoint"] == 12.0
    terminal = acceptance.extension_assessment(
        policy["criteria"], [inconclusive], "fail"
    )
    assert terminal["triggered"] is False
    assert terminal["preserved_t10_result"] == "fail"


def test_prospective_extension_artifact_binds_replayed_preserved_t10_evidence(
    fast_policy, tmp_path
):
    mhd, user = histories(tmp_path / "history")
    retained_times = {3.5, 4.0, 4.25, 4.5, 6.0, 6.25, 6.5, 8.0, 8.25, 8.5, 10.0}
    for path, label in ((mhd, "MHD"), (user, "user")):
        history, _ = acceptance.load_history(path, f"fixture {label} history")
        indices = [
            index for index, time in enumerate(history["time"])
            if time in retained_times
        ]
        write_history(path, {
            name: [values[index] for index in indices]
            for name, values in history.items()
        })
    bundle = build_case_bundle(fast_policy, tmp_path, "R14", mhd, user)
    preserved = acceptance.evaluate_case(
        fast_policy, "R14", mhd, user, None, None, bundle
    )
    assert preserved["extension_assessment"]["triggered"] is True
    prospective = fast_policy["criteria"]["extension_policy"][
        "prospective_extension_rule"
    ]
    artifact = acceptance.seal_evidence({
        "schema_version": 1,
        "record_type": prospective["required_prospective_artifact"]["record_type"],
        "authority": "non-authorizing-prospective-scientific-extension-policy",
        "non_authorizing_statement": acceptance.NON_AUTHORIZING_STATEMENT,
        "release_authorizing": False,
        "approval": {
            "status": "approved",
            "approved_before_extension_execution": True,
            "independent_of_t10_assessment": True,
            "reviewer_id": "independent-statistical-reviewer",
        },
        "bindings": {
            "preserved_t10_case_evidence_sha256": preserved["evidence_digest"]["sha256"],
            "criteria_sha256": fast_policy["criteria_binding"]["sha256"],
            "acceptance_utility_sha256": fast_policy["verified_sources"][
                "acceptance_utility"
            ]["sha256"],
        },
        "nominated_stationarity_gates": [
            item["gate"] for item in preserved["extension_assessment"]["triggered_by"]
        ],
        "prospective_extension_rule": prospective,
    })
    assert acceptance.validate_prospective_extension_artifact(
        fast_policy, preserved, artifact
    ) == artifact

    forged = deepcopy(artifact)
    forged["bindings"]["preserved_t10_case_evidence_sha256"] = "0" * 64
    forged = acceptance.seal_evidence({key: value for key, value in forged.items()
                                       if key != "evidence_digest"})
    with pytest.raises(acceptance.AcceptanceError, match="exact preserved t<=10 evidence"):
        acceptance.validate_prospective_extension_artifact(
            fast_policy, preserved, forged
        )

    changed_rule = deepcopy(artifact)
    changed_rule["prospective_extension_rule"]["fixed_combination_rule"][
        "extension_pass"
    ] = "replace every t<=10 result"
    changed_rule = acceptance.seal_evidence({
        key: value for key, value in changed_rule.items() if key != "evidence_digest"
    })
    with pytest.raises(acceptance.AcceptanceError, match="contract differs"):
        acceptance.validate_prospective_extension_artifact(
            fast_policy, preserved, changed_rule
        )

    terminal = deepcopy(preserved)
    terminal["result"] = "fail"
    terminal["extension_assessment"]["preserved_t10_result"] = "fail"
    terminal = acceptance.seal_evidence({
        key: value for key, value in terminal.items() if key != "evidence_digest"
    })
    with pytest.raises(acceptance.AcceptanceError, match="semantic contents differ"):
        acceptance.validate_prospective_extension_artifact(
            fast_policy, terminal, artifact
        )


def test_prospective_extension_rejects_self_digested_minimal_case(policy):
    prospective = policy["criteria"]["extension_policy"]["prospective_extension_rule"]
    fabricated = acceptance.seal_evidence({
        "schema_version": 1,
        "record_type": "stage-i-scientific-case-evidence",
        "result": "inconclusive",
        "extension_assessment": {
            "decision_time": 10.0,
            "preserved_t10_result": "inconclusive",
            "triggered": True,
            "triggered_by": [{
                "gate": "stationarity:magnetic",
                "trigger": "sampling_adequacy=inconclusive",
            }],
            "current_policy_authorizes_extension": False,
            "prospective_extension_rule": prospective,
        },
    })
    artifact = acceptance.seal_evidence({
        "schema_version": 1,
        "record_type": prospective["required_prospective_artifact"]["record_type"],
        "authority": "non-authorizing-prospective-scientific-extension-policy",
        "non_authorizing_statement": acceptance.NON_AUTHORIZING_STATEMENT,
        "release_authorizing": False,
        "approval": {
            "status": "approved",
            "approved_before_extension_execution": True,
            "independent_of_t10_assessment": True,
            "reviewer_id": "fabricated-reviewer",
        },
        "bindings": {
            "preserved_t10_case_evidence_sha256": fabricated["evidence_digest"]["sha256"],
            "criteria_sha256": policy["criteria_binding"]["sha256"],
            "acceptance_utility_sha256": policy["verified_sources"][
                "acceptance_utility"
            ]["sha256"],
        },
        "nominated_stationarity_gates": ["stationarity:magnetic"],
        "prospective_extension_rule": prospective,
    })
    with pytest.raises(
        acceptance.AcceptanceError,
        match="provenance|evaluation_inputs",
    ):
        acceptance.validate_prospective_extension_artifact(
            policy, fabricated, artifact
        )


def test_case_authenticates_bundle_and_uses_forcing_tcorr_but_missing_products_never_pass(
    fast_policy, tmp_path
):
    mhd, user = histories(tmp_path / "history")
    bundle = build_case_bundle(fast_policy, tmp_path, "R14", mhd, user)
    evidence = acceptance.evaluate_case(
        fast_policy, "R14", mhd, user, None, None, bundle
    )
    gates = {item["name"]: item for item in evidence["gates"]}
    assert gates["accepted_case_bundle_lineage"]["result"] == "inconclusive"
    assert gates["accepted_case_bundle_lineage"]["observations"][
        "canonical_campaign_authority_eligible"
    ] is False
    assert gates["finite_limiter_semantics"]["result"] == "pass"
    assert gates["active_pressure_work_activity"]["result"] == "pass"
    assert gates["scientific_products_contract:full"]["result"] == "inconclusive"
    assert gates["stationarity:kinetic"]["result"] == "pass"
    assert evidence["metrics"]["kinetic"]["full"]["method"]["minimum_block_duration"] == 2.0
    assert evidence["metrics"]["kinetic"]["full"]["independent_time_block_count"] == 3.0
    assert evidence["metrics"]["kinetic"]["early"]["independent_time_block_count"] == 2.0
    assert evidence["metrics"]["kinetic"]["late"]["independent_time_block_count"] == 2.0
    assert evidence["metrics"]["force_power"]["stationarity"][
        "paired_contrast_method"
    ]["contrast"] == "paired-piecewise-linear-early-minus-late"
    assert evidence["extension_assessment"]["triggered"] is False
    assert evidence["extension_assessment"]["current_policy_authorizes_extension"] is False
    assert evidence["campaign_authority_eligible"] is False
    assert evidence["result"] == "inconclusive"


def test_passive_lf_case_enforces_exact_zero_pressure_work(fast_policy, tmp_path):
    mhd, user = histories(tmp_path / "history", passive=True)
    bundle = build_case_bundle(fast_policy, tmp_path, "R06", mhd, user)
    evidence = acceptance.evaluate_case(
        fast_policy, "R06", mhd, user, None, None, bundle
    )
    gates = {item["name"]: item for item in evidence["gates"]}
    assert gates["passive_pressure_work_exact_zero"]["result"] == "pass"
    assert gates["landau_fluid_activity"]["result"] == "pass"


def test_bundle_rejects_history_not_selected_by_accepted_case(fast_policy, tmp_path):
    mhd, user = histories(tmp_path / "history")
    bundle = build_case_bundle(fast_policy, tmp_path, "R14", mhd, user)
    wrong_mhd, _ = histories(tmp_path / "wrong")
    with pytest.raises(acceptance.AcceptanceError, match="accepted bundle product"):
        acceptance.authenticate_case_bundle(
            fast_policy,
            "R14",
            bundle,
            acceptance.regular_file_binding(wrong_mhd, "wrong MHD"),
            acceptance.regular_file_binding(user, "user"),
        )


def test_arbitrary_path_bundle_lineage_is_explicitly_campaign_ineligible(
    fast_policy, tmp_path
):
    mhd, user = histories(tmp_path / "history")
    bundle = build_case_bundle(fast_policy, tmp_path, "R14", mhd, user)
    authenticated, _, lineage_gate = acceptance.authenticate_case_bundle(
        fast_policy,
        "R14",
        bundle,
        acceptance.regular_file_binding(mhd, "MHD"),
        acceptance.regular_file_binding(user, "user"),
    )
    assert authenticated["canonical_campaign_authority_eligible"] is False
    assert lineage_gate["result"] == "inconclusive"
    assert "arbitrary noncanonical path" in lineage_gate["reason"]


def test_source_archive_catalog_parser_rejects_malformed_or_duplicate_rows(
    policy, tmp_path
):
    bad = tmp_path / "SHA256SUMS"
    bad.write_text(f"{'a' * 64}  one.bundle\n{'b' * 64}  one.bundle\n")
    forged = deepcopy(policy)
    forged["verified_sources"]["source_archive_catalog"] = (
        acceptance.regular_file_binding(bad, "bad source catalog")
    )
    with pytest.raises(acceptance.AcceptanceError, match="catalog is malformed"):
        acceptance.source_archive_catalog(forged)


def test_source_archive_catalog_reauthenticates_verified_bytes(policy, tmp_path):
    catalog = tmp_path / "SHA256SUMS"
    catalog.write_text(f"{'a' * 64}  first.bundle\n")
    forged = deepcopy(policy)
    forged["verified_sources"]["source_archive_catalog"] = (
        acceptance.regular_file_binding(catalog, "source catalog")
    )
    catalog.write_text(f"{'b' * 64}  second.bundle\n")
    with pytest.raises(acceptance.AcceptanceError, match="differs from the verified binding"):
        acceptance.source_archive_catalog(forged)


def test_sampling_feasibility_change_record_is_preregistered(policy):
    forged = deepcopy(policy["criteria"])
    forged["criteria_change_record"]["current_policy"]["analysis_windows"]["full"] = [
        5.0,
        10.0,
    ]
    with pytest.raises(acceptance.AcceptanceError, match="change record differs"):
        acceptance.validate_criteria_payload(forged, policy["criteria_binding"])


def diagnostics_contract(
    policy,
    case_id: str,
    window_name: str,
    bundle_binding: dict[str, object],
    mhd_binding: dict[str, object],
    user_binding: dict[str, object],
) -> dict[str, object]:
    window = policy["criteria"]["analysis_windows"][window_name]
    observed_window = {"time_start": float(window[0]), "time_end": float(window[1])}
    name = acceptance.case_name(policy, case_id)
    return {
        "scientific_acceptance_contract": {
            "schema_version": 2,
            "record_type": "stage-i-scientific-acceptance-reviewed-products",
            "case_id": case_id,
            "case_name": name,
            "analysis_window": observed_window,
            "accepted_bundle_manifest": bundle_binding,
            "generator": policy["criteria"]["scientific_products_policy"][
                "reviewed_generator_binding"
            ],
            "deterministic_replay_verification": (
                "exact-semantic-replay-from-bound-canonical-case-inputs"
            ),
            "stage_i_manifest_sha256": policy["verified_sources"]["stage_i_manifest"]["sha256"],
            "mhd_history": mhd_binding,
            "user_history": user_binding,
        },
        "cases": {name: {"analysis_window": observed_window}},
    }


def test_hand_authored_diagnostics_never_pass_without_reviewed_generator(policy, tmp_path):
    mhd, user = histories(tmp_path / "history")
    bundle_path = tmp_path / "bundle.json"
    write_json(bundle_path, {"binding_only": True})
    bundle = {"bundle_manifest": acceptance.regular_file_binding(bundle_path, "bundle")}
    mhd_binding = acceptance.regular_file_binding(mhd, "MHD")
    user_binding = acceptance.regular_file_binding(user, "user")
    diagnostics = diagnostics_contract(
        policy, "R02", "full", bundle["bundle_manifest"], mhd_binding, user_binding
    )
    trusted, result = acceptance.validate_diagnostics_contract(
        policy, "R02", "full", diagnostics, bundle, mhd_binding, user_binding
    )
    assert trusted is None
    assert result["result"] == "inconclusive"
    assert "hand-authored products cannot pass" in result["reason"]

    forged = deepcopy(diagnostics)
    forged["scientific_acceptance_contract"]["mhd_history"]["sha256"] = "0" * 64
    trusted, result = acceptance.validate_diagnostics_contract(
        policy, "R02", "full", forged, bundle, mhd_binding, user_binding
    )
    assert trusted is None
    assert result["result"] == "inconclusive"


def panel_comparison(policy, product_id: str, normalized_shift: float) -> dict[str, object]:
    bindings = policy["manifest"]["panel_status"]["reference_product_bindings"]
    product_binding = bindings[product_id]
    source = acceptance.verified_reference_product(policy, product_id, product_binding)
    reference = [float(value) for value in source["reference"]]
    uncertainty = [float(value) for value in source["uncertainty"]]
    simulated = [
        value + normalized_shift * error
        for value, error in zip(reference, uncertainty)
    ]
    residual = [value - ref for value, ref in zip(simulated, reference)]
    record = {
        "available": True,
        "kind": product_binding["kind"],
        "case": product_binding["case"],
        "analysis_case": product_binding["case"],
        "product": product_binding["product"],
        "stage_i_binding_validated": True,
        "reference_data_file": product_binding["data_file"],
        "reference_manifest_sha256": source["manifest_sha256"],
        "data_file": source["data_path"],
        "data_sha256": source["data_sha256"],
        "interpolation": source["interpolation"],
        "sample_count": len(reference),
        "x": source["x"],
        "residual": residual,
        "rms_residual": math.sqrt(sum(value * value for value in residual) / len(residual)),
        "maximum_absolute_residual": max(abs(value) for value in residual),
        "rms_normalized_by_reported_uncertainty": abs(normalized_shift),
    }
    if product_binding["kind"] == "surface":
        record.update({
            "y": source["y"],
            "reference_z": reference,
            "reference_z_uncertainty": uncertainty,
            "simulated_z": simulated,
        })
    else:
        record.update({
            "reference_y": reference,
            "reference_y_uncertainty": uncertainty,
            "simulated_y": simulated,
        })
    return record


def panel_diagnostics(product: str, full: dict, early: dict, late: dict) -> dict:
    return {
        window: {
            "reference_curve_comparisons": {
                "available": True,
                "comparisons": {product: record},
                "surface_comparisons": {},
            }
        }
        for window, record in (("full", full), ("early", early), ("late", late))
    }


def test_panel_metric_recomputation_is_available_but_hand_authored_panel_never_passes(policy):
    product = "fig9_alignment_active_alfvenic_beta10"
    diagnostics = panel_diagnostics(
        product,
        panel_comparison(policy, product, 0.2),
        panel_comparison(policy, product, 0.1),
        panel_comparison(policy, product, 0.2),
    )
    source_record = diagnostics["full"]["reference_curve_comparisons"]["comparisons"][product]
    binding = policy["manifest"]["panel_status"]["reference_product_bindings"][product]
    rms, maximum, _ = acceptance.normalized_reference_metrics(
        policy,
        product,
        source_record,
        binding,
        acceptance.case_name(policy, "R02"),
    )
    assert rms == pytest.approx(0.2)
    assert maximum == pytest.approx(0.2)
    records = acceptance.panel_product_assessments(policy, "R02", diagnostics)
    record = next(item for item in records if item.get("product_id") == product)
    assert record["result"] == "inconclusive"
    assert "reviewed scientific-products generator" in record["reason"]


def test_panel_gate_rejects_forged_summaries_and_never_passes_unavailable_product(policy):
    product = "fig9_alignment_active_alfvenic_beta10"
    full = panel_comparison(policy, product, 0.2)
    full["rms_normalized_by_reported_uncertainty"] = 0.0
    diagnostics = panel_diagnostics(
        product,
        full,
        panel_comparison(policy, product, 0.1),
        panel_comparison(policy, product, 0.2),
    )
    binding = policy["manifest"]["panel_status"]["reference_product_bindings"][product]
    with pytest.raises(acceptance.AcceptanceError, match="internally inconsistent"):
        acceptance.normalized_reference_metrics(
            policy,
            product,
            full,
            binding,
            acceptance.case_name(policy, "R02"),
        )

    unavailable = panel_comparison(policy, product, 0.2)
    unavailable["available"] = False
    diagnostics = panel_diagnostics(
        product,
        unavailable,
        panel_comparison(policy, product, 0.1),
        panel_comparison(policy, product, 0.2),
    )
    record = next(
        item for item in acceptance.panel_product_assessments(policy, "R02", diagnostics)
        if item.get("product_id") == product
    )
    assert record["result"] == "inconclusive"


def test_analyzer_scalar_metrics_are_recomputed_from_raw_samples(fast_policy):
    name = acceptance.case_name(fast_policy, "R02")
    times = [4.0 + 0.25 * index for index in range(25)]
    values = [0.2 + 0.01 * index for index in range(25)]
    stats = acceptance.window_statistics(
        times,
        values,
        4.0,
        10.0,
        replicates=80,
        seed_text=f"analyzer:peak_alignment:full:{fast_policy['criteria_binding']['sha256']}",
        minimum_block_duration=2.0,
        expected_cadence=0.25,
    )
    metric = {
        "sample_times": times,
        "sample_values": values,
        "mean": stats["mean"],
        "standard_error": stats["standard_error"],
        "standard_deviation": stats["standard_deviation"],
    }
    diagnostics = {
        "cases": {
            name: {"scientific_acceptance_metrics": {"peak_alignment": metric}}
        }
    }
    observed = acceptance.analyzer_metrics(diagnostics, name, fast_policy, 2.0)
    assert observed["peak_alignment"]["mean"] == pytest.approx(stats["mean"])
    forged = deepcopy(diagnostics)
    forged["cases"][name]["scientific_acceptance_metrics"]["peak_alignment"]["mean"] = 99.0
    with pytest.raises(acceptance.AcceptanceError, match="differs from raw samples"):
        acceptance.analyzer_metrics(forged, name, fast_policy, 2.0)


def case_scalar(mean: float, se: float = 0.01, deviation: float = 0.1):
    return {
        "full": {"mean": mean, "standard_error": se, "standard_deviation": deviation},
        "late": {"mean": mean + 1.0, "standard_error": se, "standard_deviation": deviation},
    }


def test_holm_step_down_stops_after_first_failure(fast_policy, monkeypatch):
    active = {
        "metrics": {
            "abs_dp": case_scalar(1.0, se=0.1),
            "mirror_occupancy": case_scalar(2.0, se=0.1),
            "firehose_occupancy": case_scalar(0.0, se=0.0),
        },
        "analyzer_metrics": {
            "peak_alignment": {
                "mean": 3.0,
                "standard_error": 0.1,
                "standard_deviation": 0.1,
            }
        },
    }
    passive = {
        "metrics": {
            "abs_dp": case_scalar(0.0, se=0.0),
            "mirror_occupancy": case_scalar(0.0, se=0.0),
            "firehose_occupancy": case_scalar(0.0, se=0.0),
        },
        "analyzer_metrics": {
            "peak_alignment": {
                "mean": 0.0,
                "standard_error": 0.0,
                "standard_deviation": 0.1,
            }
        },
    }
    p_values = {10.0: 0.001, 20.0: 0.03, 30.0: 0.031}
    monkeypatch.setattr(
        acceptance,
        "normal_two_sided_p",
        lambda z_score: p_values[round(z_score, 6)],
    )
    result = acceptance.pair_contrast(active, passive, fast_policy)
    ordered = sorted(result["metrics"], key=lambda item: item["two_sided_p"])
    assert [item["holm_significant"] for item in ordered] == [True, False, False]


def test_late_window_scalar_and_exact_convergence_support_are_enforced():
    case = {"metrics": {"nu_eff": case_scalar(4.0)}}
    assert acceptance.scalar_from_case(case, "nu_eff", "late")["mean"] == 5.0
    shells = [4.0, 6.0, 8.0, 12.0, 16.0, 24.0]
    physical = [value * math.pi for value in shells]
    complete = {"x": physical, "y": [1.0] * 6}
    missing_shell = {"x": [value * math.pi for value in shells if value != 16.0], "y": [1.0] * 5}
    with pytest.raises(acceptance.AcceptanceError, match="complete preregistered shell list"):
        acceptance.curve_distance(
            complete,
            missing_shell,
            (4.0 * math.pi, 24.0 * math.pi),
            "alignment",
            physical,
        )
    short_interval = {"x": physical[1:], "y": [1.0] * 5}
    with pytest.raises(acceptance.AcceptanceError, match="full preregistered interval"):
        acceptance.curve_distance(
            complete,
            short_interval,
            (4.0 * math.pi, 24.0 * math.pi),
            "spectrum",
        )


def test_curve_improvement_is_not_applicable_when_low_mid_is_unresolved():
    shells = [value * math.pi for value in (4, 6, 8, 12, 16, 24)]
    left = {"x": shells, "y": [1.0] * 6, "standard_error": [0.2] * 6}
    mid = {"x": shells, "y": [1.01] * 6, "standard_error": [0.2] * 6}
    result = acceptance.curve_distance_resolution(
        left,
        mid,
        (4.0 * math.pi, 24.0 * math.pi),
        "alignment",
        shells,
        2.0,
    )
    assert result["resolved"] is False
    assert result["resolution_applicability"] == "not_applicable_unresolved"


def complete_case(case_id: str) -> dict[str, object]:
    mandatory = [
        "accepted_case_bundle_lineage",
        "history_exact_window_coverage",
        "scientific_products_contract:full",
        "scientific_products_contract:early",
        "scientific_products_contract:late",
        "sampled_restart_ct_divb",
    ]
    return {
        "case_id": case_id,
        "result": "pass",
        "campaign_authority_eligible": False,
        "evaluation_inputs": {
            "mhd_history": {},
            "user_history": {},
            "accepted_bundle_manifest": {},
            "ct_evidence": {},
            "diagnostics": {"full": {}, "early": {}, "late": {}},
        },
        "gates": [
            acceptance.gate(name, "pass", reason="synthetic recomputed pass")
            for name in mandatory
        ],
        "metrics": {},
        "analyzer_metrics": {},
        "convergence_products": {},
        "panel_products": [],
    }


def test_campaign_approved_review_passes_but_other_authorities_still_block(
    fast_policy, tmp_path, monkeypatch
):
    paths: list[Path] = []
    for case_id in fast_policy["criteria"]["required_cases"]:
        path = tmp_path / f"{case_id}.json"
        write_json(path, {
            "record_type": "stage-i-scientific-case-evidence",
            "case_id": case_id,
        })
        paths.append(path)
    monkeypatch.setattr(acceptance, "verify_evidence_policy_binding", lambda *args, **kwargs: 0)
    monkeypatch.setattr(
        acceptance,
        "independently_recompute_case",
        lambda _policy, value, _label: complete_case(value["case_id"]),
    )
    campaign = acceptance.evaluate_campaign(fast_policy, paths)
    review_gate = next(
        gate for gate in campaign["gates"]
        if gate["name"] == "approved_independent_criteria_review"
    )
    lf_gate = next(
        gate for gate in campaign["gates"]
        if gate["name"] == "lf_strength_resolved_response"
    )
    assert review_gate["result"] == "pass"
    assert ["R06", "R02"] in [record["pair"] for record in lf_gate["observations"]]
    assert all(
        "missing_dependency" in record
        for record in lf_gate["observations"]
        if not record["available"]
    )
    assert campaign["result"] != "pass"
    canonical_gate = next(
        gate for gate in campaign["gates"]
        if gate["name"] == "canonical_campaign_authority"
    )
    assert canonical_gate["result"] == "inconclusive"
    products_gate = next(
        gate for gate in campaign["gates"]
        if gate["name"] == "reviewed_scientific_products_generator"
    )
    convergence_gate = next(
        gate for gate in campaign["gates"]
        if gate["name"] == "R16_R02_R17_resolution_convergence"
    )
    assert products_gate["result"] == "inconclusive"
    assert convergence_gate["result"] == "inconclusive"


def test_campaign_rejects_incomplete_or_voluntary_case_declaration(
    fast_policy, tmp_path
):
    case = acceptance.seal_evidence({
        "schema_version": 1,
        "record_type": "stage-i-scientific-case-evidence",
        "case_id": "R02",
        "result": "pass",
        "provenance": {
            "criteria": fast_policy["criteria_binding"],
            "criteria_review": fast_policy["review_binding"],
            "acceptance_utility": acceptance.regular_file_binding(UTILITY, "utility"),
            "inputs": [],
        },
    })
    path = tmp_path / "R02.json"
    write_json(path, case)
    with pytest.raises(acceptance.AcceptanceError, match="evaluation_inputs"):
        acceptance.evaluate_campaign(fast_policy, [path])


def test_verify_evidence_recomputes_semantics_not_only_self_digest(policy, tmp_path):
    evidence = acceptance.validate_criteria_evidence(policy)
    path = tmp_path / "criteria-validation.json"
    write_json(path, evidence)
    verification = acceptance.verify_evidence(policy, path, sha256(path))
    assert verification["verified"] is True

    forged = deepcopy(evidence)
    forged.pop("evidence_digest")
    forged["valid"] = False
    forged = acceptance.seal_evidence(forged)
    forged_path = tmp_path / "forged-validation.json"
    write_json(forged_path, forged)
    with pytest.raises(acceptance.AcceptanceError, match="semantic contents differ"):
        acceptance.verify_evidence(policy, forged_path, sha256(forged_path))


def test_verify_evidence_replays_truthful_incomplete_case_with_optional_inputs_absent(
    fast_policy, tmp_path
):
    mhd, user = histories(tmp_path / "history", passive=True)
    bundle = build_case_bundle(fast_policy, tmp_path, "R06", mhd, user)
    evidence = acceptance.evaluate_case(
        fast_policy, "R06", mhd, user, None, None, bundle
    )
    assert evidence["result"] == "inconclusive"
    path = tmp_path / "R06-incomplete.json"
    write_json(path, evidence)
    verification = acceptance.verify_evidence(fast_policy, path, sha256(path))
    assert verification["verified"] is True
    assert verification["evidence_result"] == "inconclusive"


def test_candidate_output_is_immutable_and_no_clobber(tmp_path):
    path = tmp_path / "candidate.json"
    acceptance.write_candidate(path, acceptance.seal_evidence({"value": 1}))
    assert path.stat().st_mode & 0o777 == 0o444
    with pytest.raises(acceptance.AcceptanceError, match="already exists"):
        acceptance.write_candidate(path, acceptance.seal_evidence({"value": 2}))
