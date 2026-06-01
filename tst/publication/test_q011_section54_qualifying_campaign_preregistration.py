#!/usr/bin/env python3
"""Fail-closed checks for the additive Q-011 Section 5.4 campaign policy."""

from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import unittest

try:
    from tst.publication import analyze_q011_section54_outputs as q011
except ModuleNotFoundError:
    import analyze_q011_section54_outputs as q011


POLICY = (
    Path(__file__).resolve().parent
    / "readiness"
    / "q011_section54_qualifying_campaign_preregistration_successor_v2_2026-06-01.json"
)


class PolicyError(ValueError):
    """Raised when the Q-011 preregistration schema or frozen policy drifts."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PolicyError(message)


def _keys(value: object, expected: set[str], label: str) -> dict[str, object]:
    _require(type(value) is dict, f"{label}: expected object")
    mapping = value
    _require(set(mapping) == expected, f"{label}: schema drift")
    return mapping


def _strict_equal(actual: object, expected: object, label: str) -> None:
    _require(type(actual) is type(expected), f"{label}: scalar type drift")
    if type(expected) is dict:
        actual_mapping = _keys(actual, set(expected), label)
        for key, expected_value in expected.items():
            _strict_equal(actual_mapping[key], expected_value, f"{label}/{key}")
    elif type(expected) is list:
        _require(len(actual) == len(expected), f"{label}: list length drift")
        for index, (actual_value, expected_value) in enumerate(zip(actual, expected)):
            _strict_equal(actual_value, expected_value, f"{label}[{index}]")
    elif type(expected) is float:
        _require(math.isfinite(actual), f"{label}: non-finite number")
        _require(actual == expected, f"{label}: numeric drift")
    else:
        _require(actual == expected, f"{label}: value drift")


def _finite_float(value: object, label: str) -> float:
    _require(type(value) is float, f"{label}: expected float")
    _require(math.isfinite(value), f"{label}: non-finite number")
    return value


def _reject_constant(value: str) -> None:
    raise PolicyError(f"JSON constant is forbidden: {value}")


def _decode_policy(payload: str) -> dict[str, object]:
    try:
        decoded = json.loads(payload, parse_constant=_reject_constant)
    except json.JSONDecodeError as exc:
        raise PolicyError("invalid JSON") from exc
    _require(type(decoded) is dict, "policy: expected object")
    return decoded


def _load_policy() -> dict[str, object]:
    return _decode_policy(POLICY.read_text(encoding="utf-8"))


def _validate_policy(policy: object) -> None:
    root = _keys(
        policy,
        {
            "record_type",
            "schema_version",
            "date",
            "gate",
            "claim_id",
            "qualification_effect",
            "predecessor_record",
            "scope",
            "schema_contract",
            "analysis_primitive_bindings",
            "paper_literal_observables",
            "athenak_selected_release_criteria",
            "artifact_retention_policy",
            "failure_policy",
            "qualifying_execution_bindings",
        },
        "policy",
    )
    _strict_equal(
        {
            key: root[key]
            for key in (
                "record_type",
                "schema_version",
                "date",
                "gate",
                "claim_id",
                "qualification_effect",
            )
        },
        {
            "record_type": "q011_section54_qualifying_campaign_preregistration",
            "schema_version": 1,
            "date": "2026-06-01",
            "gate": "Q-011",
            "claim_id": "CLAIM-PAPER-SHOCK-001",
            "qualification_effect": (
                "policy_freeze_only_no_execution_authorization_no_claim_closure"
            ),
        },
        "policy/identity",
    )
    _require(type(root["scope"]) is str and root["scope"], "policy/scope: expected text")
    _strict_equal(
        root["schema_contract"],
        {
            "schema_style": "self_contained_exact_key_policy",
            "unknown_keys": "reject",
            "missing_keys": "reject",
            "numeric_aliases": "reject_boolean_nonfinite_and_numeric_string_aliases",
            "policy_drift": (
                "requires_versioned_successor_before_qualifying_output_inspection"
            ),
        },
        "policy/schema_contract",
    )
    _strict_equal(
        root["analysis_primitive_bindings"],
        {
            "module": "tst/publication/analyze_q011_section54_outputs.py",
            "role": "policy_free_arithmetic_foundation_only",
            "fixed_histogram": "fixed_histogram",
            "slope_fit": "fit_fixed_bin_loglog_slope",
            "shock_front": "detect_shock_front",
            "upstream_magnetic_amplification": "estimate_upstream_amplification",
            "amr_vs_fine_uniform_residuals": "matched_grid_residual_metrics",
        },
        "policy/analysis_primitive_bindings",
    )

    literal = _keys(
        root["paper_literal_observables"],
        {"source_anchor", "tolerance_boundary", "observables"},
        "policy/paper_literal_observables",
    )
    _require("Section 5.4" in literal["source_anchor"], "paper source anchor drift")
    _require(
        "not a claimed manuscript tolerance" in literal["tolerance_boundary"],
        "paper tolerance boundary drift",
    )
    observables = literal["observables"]
    _require(type(observables) is list, "paper observables: expected list")
    expected_literal_ids = [
        "ideal_injection_surface",
        "startup_cohort_removal",
        "snapshot_morphology_at_t500",
        "upstream_magnetic_amplification_at_t500",
        "downstream_spectra_at_t500_and_t1200",
        "late_energy_tail_slope_at_t1200",
    ]
    _require(
        [item["id"] for item in observables] == expected_literal_ids,
        "paper observables: id drift",
    )
    for index, observable in enumerate(observables):
        item = _keys(
            observable,
            {"id", "classification", "statement"},
            f"paper observables[{index}]",
        )
        _require(
            item["classification"] == "paper_literal_observable",
            f"paper observables[{index}]: classification drift",
        )
        _require(
            type(item["statement"]) is str and item["statement"],
            f"paper observables[{index}]: expected statement",
        )

    criteria = _keys(
        root["athenak_selected_release_criteria"],
        {
            "criteria_provenance",
            "campaign_matrix",
            "snapshot_selection",
            "ideal_injection_surface_classifier",
            "particle_filter",
            "spectrum",
            "shock_front",
            "upstream_magnetic_amplification",
            "amr_vs_fine_uniform_residuals",
            "qualitative_figure_requirements",
        },
        "policy/athenak_selected_release_criteria",
    )
    _require(
        "not tolerances quoted or implied by the manuscript"
        in criteria["criteria_provenance"],
        "release-criteria provenance drift",
    )
    _strict_equal(
        criteria["campaign_matrix"],
        {
            "physical_mode": "paper_mhd_pic_vl2_tsc",
            "grid_variants": [
                "coarse_uniform_dx12",
                "three_level_amr_root_dx12_finest_dx3",
                "fine_uniform_dx3",
            ],
            "qualifying_seeds": [
                23050101,
                23050102,
                23050103,
                23050104,
                23050105,
                23050106,
                23050107,
                23050108,
            ],
            "expected_baseline_attempts": 24,
            "paired_seed_rule": (
                "Use the same qualifying seed for coarse-uniform, AMR and "
                "fine-uniform variants."
            ),
        },
        "criteria/campaign_matrix",
    )
    _strict_equal(
        criteria["snapshot_selection"],
        {
            "time_unit": "omega0_inverse",
            "required_times": [500.0, 1200.0],
            "absolute_match_tolerance": 1.0e-06,
            "missing_or_ambiguous_snapshot": "fail_endpoint",
        },
        "criteria/snapshot_selection",
    )
    _strict_equal(
        criteria["ideal_injection_surface_classifier"],
        {
            "formula": "x_ideal(t) = ((gamma - 1) * u0 / 2) * t",
            "gamma": 1.6666666666666667,
            "u0_over_u_a0": 30.0,
            "ideal_surface_speed_over_u_a0": 10.0,
            "downstream_rule": "x1 < x_ideal(t)",
            "upstream_rule": "x1 > x_ideal(t)",
            "on_surface_rule": (
                "x1 == x_ideal(t) belongs to neither open analysis region; "
                "injection remains assigned to the ideal surface"
            ),
            "detected_front_substitution": "forbidden",
        },
        "criteria/ideal_injection_surface_classifier",
    )
    _strict_equal(
        criteria["particle_filter"],
        {
            "source_scalar": "cr_source",
            "required_source": "shock_injected",
            "birth_time_scalar": "birth_time",
            "birth_time_operator": ">=",
            "birth_time_min_omega0_inverse": 45.0,
            "spatial_region": "downstream_by_ideal_injection_surface_classifier",
            "rejected_particle_policy": (
                "exclude_from_spectrum_and_archive_filter_counts"
            ),
        },
        "criteria/particle_filter",
    )

    spectrum = _keys(
        criteria["spectrum"],
        {
            "energy_variable",
            "energy_formula",
            "epsilon_formula",
            "histogram_weight",
            "distribution_for_fit",
            "plotted_quantity",
            "bin_spacing",
            "bin_edges",
            "underflow_policy",
            "overflow_policy",
            "max_overflow_macro_weight_fraction",
            "slope_snapshot_omega0_inverse",
            "slope_fit_window_chi",
            "slope_fit_rationale",
            "minimum_positive_fit_bins",
            "slope_target",
            "slope_absolute_tolerance",
            "slope_scope",
            "fit_failure_policy",
        },
        "criteria/spectrum",
    )
    _strict_equal(
        {key: spectrum[key] for key in spectrum if key != "bin_edges"},
        {
            "energy_variable": "chi",
            "energy_formula": "chi = 2 * epsilon / (m * u0^2)",
            "epsilon_formula": "epsilon = p_cr^2 / (2 * m)",
            "histogram_weight": "macro_weight",
            "distribution_for_fit": "f_chi = admitted_macro_weight / delta_chi",
            "plotted_quantity": (
                "chi * f_chi normalized by total admitted downstream macro weight"
            ),
            "bin_spacing": "logarithmic_fixed_quarter_octave",
            "underflow_policy": "archive_count_and_macro_weight",
            "overflow_policy": (
                "archive_count_and_macro_weight_and_fail_if_macro_weight_fraction_"
                "exceeds_limit"
            ),
            "max_overflow_macro_weight_fraction": 0.001,
            "slope_snapshot_omega0_inverse": 1200.0,
            "slope_fit_window_chi": [20.0, 160.0],
            "slope_fit_rationale": (
                "The fixed window starts above the mono-energetic injection "
                "feature at chi = 10 and ends below the visible high-energy "
                "rollover in the manuscript figure. It is selected before "
                "qualifying output inspection."
            ),
            "minimum_positive_fit_bins": 8,
            "slope_target": -1.5,
            "slope_absolute_tolerance": 0.2,
            "slope_scope": [
                "three_level_amr_root_dx12_finest_dx3",
                "fine_uniform_dx3",
            ],
            "fit_failure_policy": (
                "fail_endpoint_and_archive_fit_inputs; "
                "do_not_tune_bins_or_window_after_output_inspection"
            ),
        },
        "criteria/spectrum_without_edges",
    )
    edges = spectrum["bin_edges"]
    _require(type(edges) is list and len(edges) == 41, "spectrum edges: shape drift")
    for index, edge in enumerate(edges):
        _finite_float(edge, f"spectrum edges[{index}]")
    _require(edges[0] == 1.0 and edges[-1] == 1024.0, "spectrum edges: bounds drift")
    ratio = 2.0 ** 0.25
    for lower, upper in zip(edges, edges[1:]):
        _require(
            math.isclose(upper / lower, ratio, rel_tol=1.0e-14, abs_tol=0.0),
            "spectrum edges: logarithmic spacing drift",
        )

    _strict_equal(
        criteria["shock_front"],
        {
            "source_field": "rho",
            "profile": "y_area_weighted_mean",
            "detector": "unique_strongest_positive_density_gradient",
            "search_window_relative_to_x_ideal_c_over_omega_pi": [-1200.0, 1200.0],
            "max_absolute_offset_from_x_ideal_c_over_omega_pi": 600.0,
            "search_window_rationale": (
                "The fixed search window spans 100 root cells on either side of "
                "the ideal injection surface. It admits corrugation while "
                "excluding distant structures and may not be tuned after inspection."
            ),
        },
        "criteria/shock_front",
    )
    _strict_equal(
        criteria["upstream_magnetic_amplification"],
        {
            "source_field": "bmag",
            "observable": "area_weighted_mean_abs_b_over_b0",
            "snapshot_omega0_inverse": 500.0,
            "window_relative_to_x_ideal_c_over_omega_pi": [120.0, 1200.0],
            "reference_b0": 1.0,
            "acceptance_range": [1.2, 3.5],
            "rationale": (
                "The manuscript's approximate 2-4 factor refers to certain "
                "regions, not an area-weighted mean. The fixed mean-field interval "
                "is a separate AthenaK release criterion requiring resolved "
                "upstream amplification without claiming a manuscript tolerance."
            ),
        },
        "criteria/upstream_magnetic_amplification",
    )
    residuals = _keys(
        criteria["amr_vs_fine_uniform_residuals"],
        {
            "comparison_rule",
            "grid_matching_rule",
            "normalization_rule",
            "rows",
            "maximum_norm_policy",
            "tolerance_rationale",
            "threshold_boundary",
        },
        "criteria/amr_vs_fine_uniform_residuals",
    )
    expected_residual_rows = [
        {
            "observable": "shock_front_position_at_t500",
            "comparison_domain": "detected_front",
            "max_absolute_difference": 240.0,
            "max_relative_mean_absolute": None,
            "max_relative_root_mean_square": None,
        },
        {
            "observable": "upstream_magnetic_amplification_at_t500",
            "comparison_domain": "fixed_upstream_window",
            "max_absolute_difference": 0.35,
            "max_relative_mean_absolute": None,
            "max_relative_root_mean_square": None,
        },
        {
            "observable": "rho_y_average_at_t500",
            "comparison_domain": "x_ideal_plus_or_minus_1200_c_over_omega_pi",
            "max_absolute_difference": None,
            "max_relative_mean_absolute": 0.2,
            "max_relative_root_mean_square": 0.3,
        },
        {
            "observable": "bmag_y_average_at_t500",
            "comparison_domain": "x_ideal_plus_or_minus_1200_c_over_omega_pi",
            "max_absolute_difference": None,
            "max_relative_mean_absolute": 0.25,
            "max_relative_root_mean_square": 0.35,
        },
        {
            "observable": "normalized_downstream_chi_f_chi_at_t500",
            "comparison_domain": "fixed_chi_bins_10_to_256",
            "max_absolute_difference": None,
            "max_relative_mean_absolute": 0.2,
            "max_relative_root_mean_square": 0.3,
        },
        {
            "observable": "normalized_downstream_chi_f_chi_at_t1200",
            "comparison_domain": "fixed_chi_bins_10_to_512",
            "max_absolute_difference": None,
            "max_relative_mean_absolute": 0.2,
            "max_relative_root_mean_square": 0.3,
        },
    ]
    _strict_equal(
        residuals,
        {
            "comparison_rule": (
                "Compare each AMR qualifying seed only with the fine-uniform run "
                "carrying the same seed."
            ),
            "grid_matching_rule": (
                "For spatial profiles, conservatively area-restrict both variants "
                "to dx = 12 c/omega_pi and y-average before matched-grid residual "
                "metrics. For spectra, compare the same fixed chi bins after the "
                "particle filter."
            ),
            "normalization_rule": (
                "Normalize spectrum shapes by total admitted downstream macro "
                "weight before residual metrics."
            ),
            "rows": expected_residual_rows,
            "maximum_norm_policy": (
                "Archive maximum-absolute residuals but do not gate on a per-cell "
                "maximum norm because turbulent feature phase is not a stable "
                "cellwise oracle."
            ),
            "tolerance_rationale": (
                "The paired mean and RMS bounds are conservative AthenaK release "
                "screens: integrated spectra must remain close, while turbulent "
                "spatial profiles receive a modestly wider RMS allowance. They "
                "are frozen before qualifying output inspection."
            ),
            "threshold_boundary": (
                "These paired residual bounds are AthenaK release criteria, not "
                "manuscript tolerances."
            ),
        },
        "criteria/amr_vs_fine_uniform_residuals",
    )

    figures = criteria["qualitative_figure_requirements"]
    _strict_equal(
        figures,
        [
            {
                "id": "morphology_t500",
                "snapshot_omega0_inverse": 500.0,
                "required_panels": ["rho", "bmag", "filtered_cr_distribution"],
                "categorical_requirements": [
                    "show_shock_corrugation",
                    "show_upstream_density_filaments_and_cavities",
                    "show_regions_of_upstream_magnetic_amplification",
                    "show_meshblock_or_refinement_boundaries",
                ],
            },
            {
                "id": "profile_t500",
                "snapshot_omega0_inverse": 500.0,
                "required_panels": ["prtcl_jx_y_average", "bmag_y_average"],
                "categorical_requirements": [
                    "mark_ideal_injection_surface",
                    "mark_detected_shock_front_separately",
                    "overlay_coarse_uniform_amr_and_fine_uniform",
                ],
            },
            {
                "id": "spectra_t500_t1200",
                "snapshot_omega0_inverse": 1200.0,
                "required_panels": ["downstream_chi_f_chi"],
                "categorical_requirements": [
                    "use_fixed_chi_axis",
                    "overlay_coarse_uniform_amr_and_fine_uniform",
                    "distinguish_t500_and_t1200",
                    "show_late_slope_reference",
                    "record_external_reviewer_disposition",
                ],
            },
        ],
        "criteria/qualitative_figure_requirements",
    )

    _strict_equal(
        root["artifact_retention_policy"],
        {
            "retention_scope": (
                "retain every emitted raw artifact for every attempted baseline run, "
                "including failed attempts"
            ),
            "required_output_times_omega0_inverse": [
                0.0,
                100.0,
                200.0,
                300.0,
                400.0,
                500.0,
                600.0,
                700.0,
                800.0,
                900.0,
                1000.0,
                1100.0,
                1200.0,
            ],
            "raw_bin_ids": ["rho", "bmag", "prtcl_jx", "j2"],
            "raw_pvtk_ids": ["prtcl_all"],
            "raw_rst_policy": (
                "retain_each_emitted_restart_checkpoint_and_completion_metadata"
            ),
            "inventory_policy": (
                "write_root_relative_sha256_inventory_and_preserve_raw_bin_pvtk_rst_bytes"
            ),
            "derived_artifact_policy": (
                "retain_filter_counts_histograms_fit_inputs_residual_tables_figures_"
                "and_reviewer_disposition_without_replacing_raw_artifacts"
            ),
        },
        "policy/artifact_retention_policy",
    )
    _strict_equal(
        root["failure_policy"],
        {
            "archive_every_attempt": True,
            "archive_failed_attempts": True,
            "archive_finite_outliers": True,
            "delete_or_overwrite_attempt_artifacts": "forbidden",
            "failed_fit_or_missing_snapshot": "archive_available_artifacts_and_fail_endpoint",
            "post_inspection_threshold_or_window_tuning": (
                "forbidden_requires_versioned_successor_and_new_qualifying_dataset"
            ),
            "replacement_attempt_rule": (
                "archive_original_attempt_and_register_replacement_as_a_distinct_attempt"
            ),
        },
        "policy/failure_policy",
    )
    _strict_equal(
        root["qualifying_execution_bindings"],
        {
            "status": "open_run_blocked_until_immutable_campaign_record_binds_all_fields",
            "required_before_run": [
                "clean_candidate_git_commit",
                "clean_frontier_executable_sha256",
                "qualifying_input_deck_sha256_per_variant",
                "qualifying_analyzer_sha256",
                "authorized_orion_campaign_root",
                "registered_frontier_submission_policy",
                "independent_raw_artifact_recompute_plan",
            ],
            "campaign_results_inspected": False,
            "frontier_execution_authorized_by_this_record": False,
        },
        "policy/qualifying_execution_bindings",
    )


class Q011Section54QualifyingCampaignPreregistrationTests(unittest.TestCase):
    def test_checked_in_policy_is_self_contained_and_valid(self) -> None:
        _validate_policy(_load_policy())

    def test_bound_policy_free_analyzer_primitives_exist(self) -> None:
        bindings = _load_policy()["analysis_primitive_bindings"]
        for policy_key in (
            "fixed_histogram",
            "slope_fit",
            "shock_front",
            "upstream_magnetic_amplification",
            "amr_vs_fine_uniform_residuals",
        ):
            with self.subTest(policy_key=policy_key):
                self.assertTrue(callable(getattr(q011, bindings[policy_key])))

    def test_unknown_and_missing_keys_fail_closed(self) -> None:
        extra = _load_policy()
        extra["unexpected"] = "fail"
        with self.assertRaisesRegex(PolicyError, "schema drift"):
            _validate_policy(extra)

        missing = _load_policy()
        del missing["artifact_retention_policy"]["raw_rst_policy"]
        with self.assertRaisesRegex(PolicyError, "schema drift"):
            _validate_policy(missing)

        nested_extra = _load_policy()
        nested_extra["athenak_selected_release_criteria"]["spectrum"]["adaptive_bins"] = True
        with self.assertRaisesRegex(PolicyError, "schema drift"):
            _validate_policy(nested_extra)

    def test_scalar_aliases_and_nonfinite_values_fail_closed(self) -> None:
        for alias in (False, "1.0e-6", math.nan):
            with self.subTest(alias=alias):
                drift = _load_policy()
                drift["athenak_selected_release_criteria"]["snapshot_selection"][
                    "absolute_match_tolerance"
                ] = alias
                with self.assertRaises(PolicyError):
                    _validate_policy(drift)

        with self.assertRaisesRegex(PolicyError, "JSON constant is forbidden"):
            _decode_policy('{"value": NaN}')

    def test_snapshot_classifier_filter_and_energy_variable_drift_fail_closed(self) -> None:
        mutations = [
            (
                ("snapshot_selection", "required_times"),
                [500.0],
            ),
            (
                ("ideal_injection_surface_classifier", "detected_front_substitution"),
                "allowed",
            ),
            (
                ("particle_filter", "birth_time_operator"),
                ">",
            ),
            (
                ("spectrum", "energy_formula"),
                "chi = epsilon / (m * u0^2)",
            ),
            (
                ("spectrum", "slope_fit_window_chi"),
                [30.0, 200.0],
            ),
        ]
        for path, replacement in mutations:
            with self.subTest(path=path):
                drift = _load_policy()
                drift["athenak_selected_release_criteria"][path[0]][path[1]] = replacement
                with self.assertRaises(PolicyError):
                    _validate_policy(drift)

    def test_log_bins_and_release_threshold_drift_fail_closed(self) -> None:
        drift = _load_policy()
        drift["athenak_selected_release_criteria"]["spectrum"]["bin_edges"][20] = 33.0
        with self.assertRaisesRegex(PolicyError, "spacing drift"):
            _validate_policy(drift)

        drift = _load_policy()
        drift["athenak_selected_release_criteria"]["upstream_magnetic_amplification"][
            "acceptance_range"
        ][0] = 2.0
        with self.assertRaises(PolicyError):
            _validate_policy(drift)

        drift = _load_policy()
        drift["athenak_selected_release_criteria"]["amr_vs_fine_uniform_residuals"][
            "rows"
        ][2]["max_relative_mean_absolute"] = 0.5
        with self.assertRaises(PolicyError):
            _validate_policy(drift)

    def test_failure_archive_and_raw_retention_policy_drift_fail_closed(self) -> None:
        drift = _load_policy()
        drift["failure_policy"]["archive_failed_attempts"] = False
        with self.assertRaises(PolicyError):
            _validate_policy(drift)

        drift = _load_policy()
        drift["artifact_retention_policy"]["raw_pvtk_ids"] = []
        with self.assertRaises(PolicyError):
            _validate_policy(drift)

        drift = _load_policy()
        drift["artifact_retention_policy"]["raw_rst_policy"] = "retain_t1200_only"
        with self.assertRaises(PolicyError):
            _validate_policy(drift)


if __name__ == "__main__":
    unittest.main()
