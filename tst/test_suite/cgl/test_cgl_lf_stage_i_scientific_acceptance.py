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


def test_preregistered_criteria_bind_final_utility_but_remain_pending(policy):
    assert policy["review_status"] == "changes_required"
    assert policy["approved"] is False
    assert policy["review"]["reviews"] == []
    utility_sha = acceptance.regular_file_binding(UTILITY, "utility")["sha256"]
    assert policy["criteria"]["source_bindings"]["acceptance_utility"]["sha256"] == utility_sha
    assert policy["review"]["acceptance_utility"]["sha256"] == utility_sha
    assert policy["review"]["criteria"]["sha256"] == policy["criteria_binding"]["sha256"]
    assert policy["criteria"]["family_gates"]["lf_strength"]["cases"] == [
        "R12", "R02", "R06", "R13"
    ]
    assert policy["criteria"]["statistics_policy"]["minimum_independent_time_blocks"] == {
        "full": 3.0,
        "half": 1.5,
    }
    assert policy["criteria"]["analysis_windows"] == {
        "full": [4.0, 10.0],
        "early": [4.0, 7.0],
        "late": [7.0, 10.0],
    }
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
    assert evidence["independent_review_complete"] is False


def test_approved_review_schema_requires_declared_independence(policy):
    review = deepcopy(policy["review"])
    review["review_status"] = "approved"
    review["reviews"] = [
        {
            "role": "plasma_physics",
            "reviewer_id": "reviewer-a",
            "decision": "approved",
        },
        {
            "role": "statistical_methodology",
            "reviewer_id": "reviewer-b",
            "decision": "approved",
        },
    ]
    with pytest.raises(acceptance.AcceptanceError, match="reviewer is malformed"):
        acceptance.validate_criteria_review(
            review,
            policy["review_binding"],
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


def test_effective_sample_estimator_detects_correlation():
    independent = [(-1.0) ** index for index in range(20)]
    correlated = [0.0] * 10 + [1.0] * 10
    independent_neff, _ = acceptance.effective_sample_size(independent)
    correlated_neff, _ = acceptance.effective_sample_size(correlated)
    assert independent_neff > correlated_neff


def test_stationarity_fails_large_resolved_half_window_drift():
    policy = {
        "z_lte": 3.0,
        "scalar_relative_change_lte": 0.25,
        "occupancy_absolute_change_lte": 0.002,
        "forcing_power_relative_change_lte": 0.1,
    }
    result = acceptance.stationarity_result(
        {"mean": 1.0},
        {"mean": 0.5, "standard_error": 0.01},
        {"mean": 1.5, "standard_error": 0.01},
        "scalar",
        policy,
    )
    assert result["result"] == "fail"
    assert result["z_score"] > 3.0


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
    assert evidence["metrics"]["kinetic"]["early"]["independent_time_block_count"] == 1.5
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


def test_campaign_pending_review_cannot_pass_even_with_all_apparent_case_passes(
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
    assert review_gate["result"] == "inconclusive"
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
