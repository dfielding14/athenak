#!/usr/bin/env python3
"""Preregister excluded Q019 physical-window pilots and selection rules.

This module is intentionally non-authorizing.  It freezes which existing
preproduction rows may be inspected to choose finite-rigidity fit windows,
nonlinear onset/saturation windows, and a later production horizon.  The
future qualifying ensemble uses a disjoint seed inventory.
"""

from __future__ import annotations

import hashlib
import json
import math
import statistics
from typing import Mapping, Sequence

from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as design


SCHEMA_VERSION = 1
RECORD_TYPE = "q019_excluded_physical_window_preregistration_v1"
STATUS = "preregistered_excluded_physical_pilots_non_authorizing"

PREDECESSOR_CORE = (
    "q019-fr-predecessor-k8-rho3em06-s0",
    "q019-fr-predecessor-fiducial-ppc24-s0",
    "q019-fr-predecessor-fiducial-ppc48-s0",
    "q019-fr-predecessor-fiducial-resolution-coarse-s0",
    "q019-fr-predecessor-fiducial-resolution-fiducial-s0",
    "q019-fr-predecessor-fiducial-particle-step-small-s0",
    "q019-fr-predecessor-fiducial-noise-seeded-s0",
)
WINDOW_CORE = (
    "q019-hr-current-retention-s0",
    "q019-fr-grid-k8-rho3em06-s0",
    "q019-fr-3d-onset-small-s0",
    "q019-fr-3d-onset-large-s0",
    "q019-fr-3d-onset-small-ppc48-s0",
    "q019-fr-3d-onset-small-resolution-fine-s0",
    "q019-fr-3d-onset-small-particle-step-small-s0",
    "q019-fr-3d-onset-large-long-mode-sensitivity-s0",
)
WINDOW_REPLICATION = (
    "q019-fr-3d-onset-small-s1",
    "q019-fr-3d-onset-large-s1",
    "q019-fr-3d-onset-small-s2",
    "q019-fr-3d-onset-large-s2",
)

QUALIFYING_FIELD_SEEDS = (39019, 39023, 39041, 39043, 39047, 39079, 39089, 39103)
QUALIFYING_PARTICLE_SEEDS = (
    49019,
    49031,
    49033,
    49037,
    49043,
    49057,
    49069,
    49081,
)

ONSET_BPERP_RMS_OVER_B0 = 1.0
MINIMUM_PLATEAU_SPAN_TAU = 4.0
MAXIMUM_PLATEAU_SPAN_TAU = 8.0
MINIMUM_PLATEAU_SAMPLES = 5
MAXIMUM_ABSOLUTE_LOG_ENERGY_SLOPE_PER_TAU = 0.10
MAXIMUM_PLATEAU_ENERGY_MAX_TO_MIN = 1.50
MINIMUM_POST_PLATEAU_COVERAGE_TAU = 2.0
PRODUCTION_POST_PLATEAU_COVERAGE_TAU = 8.0
PRODUCTION_HORIZON_MULTIPLIER = 1.50

MAXIMUM_PAIRED_RELATIVE_DIFFERENCE = 0.15
MAXIMUM_PAIRED_TIME_DIFFERENCE_TAU = 2.0
MAXIMUM_CONSERVATION_FRACTIONAL_RESIDUAL = 1.0e-3
MINIMUM_CHARACTERISTIC_RL_OVER_DX = 8.0
MAXIMUM_BOX_EDGE_POWER_FRACTION = 0.10
MINIMUM_HIGH_RIGIDITY_CURRENT_RETENTION = 0.95
MAXIMUM_HIGH_RIGIDITY_CR_MOMENTUM_CHANGE = 0.05
MAXIMUM_HIGH_RIGIDITY_CR_ENERGY_CHANGE = 0.05
MAXIMUM_ABSOLUTE_LOCAL_HALL_PARAMETER = 0.10
MAXIMUM_LOCAL_BAI_CHARGE_FRACTION = 0.01
MAXIMUM_PILOT_NODE_HOURS = 500.0
MAXIMUM_UNRESERVED_BUDGET_FRACTION = 0.10
MAXIMUM_ATTEMPTS_PER_CASE = 1
MAXIMUM_RETRIES_PER_CASE = 0

AUTHORIZATION = {
    "launch_authorized": False,
    "scheduler_submission_authorized": False,
    "policy_mutation_authorized": False,
    "production_resource_freeze_authorized": False,
    "production_deck_freeze_authorized": False,
    "q019_qualification_authorized": False,
    "nonlinear_saturation_claim_authorized": False,
    "scientific_claim_authorized": False,
    "publication_authorized": False,
}


class PreregistrationError(ValueError):
    """Reject drifted pilot inventories or post-inspection rule changes."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PreregistrationError(message)


def _canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _case_bindings(case_ids: Sequence[str]) -> list[dict[str, object]]:
    cases = {str(case["case_id"]): case for case in design.expected_cases()}
    _require(len(cases) == 77, "Q019 source design inventory drifted")
    bindings = []
    for case_id in case_ids:
        _require(case_id in cases, f"unknown Q019 physical pilot: {case_id}")
        case = cases[case_id]
        bindings.append(
            {
                "case_id": case_id,
                "campaign_id": case["campaign_id"],
                "branch": case["branch"],
                "role": case["role"],
                "dimension": case["dimension"],
                "field_seed": case["field_seed"],
                "particle_seed": case["particle_seed"],
                "matrix_identity_fingerprint": case[
                    "matrix_identity_fingerprint"
                ],
                "saturation_candidate": False,
            }
        )
    return bindings


def build_preregistration() -> dict[str, object]:
    pilot_ids = PREDECESSOR_CORE + WINDOW_CORE + WINDOW_REPLICATION
    _require(len(pilot_ids) == len(set(pilot_ids)) == 19, "pilot inventory drifted")
    bindings = _case_bindings(pilot_ids)
    pilot_field_seeds = {int(binding["field_seed"]) for binding in bindings}
    pilot_particle_seeds = {int(binding["particle_seed"]) for binding in bindings}
    _require(
        pilot_field_seeds.isdisjoint(QUALIFYING_FIELD_SEEDS)
        and pilot_particle_seeds.isdisjoint(QUALIFYING_PARTICLE_SEEDS),
        "pilot and qualifying seeds overlap",
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "status": STATUS,
        "campaign": "q019_nonlinear_bell_registered_successor_v1",
        "purpose": (
            "select_physical_fit_onset_plateau_box_and_production_horizon_"
            "rules_from_excluded_pilots_before_qualifying_output"
        ),
        "prerequisites": [
            "complete_passing_registered_Q043_matrix_and_exact_retirement",
            "complete_passing_registered_Q023_linear_matrix_and_exact_retirement",
            "passing_Q019_controller_resource_engineering_gate",
            "fresh_clean_candidate_controller_storage_and_budget_bindings",
        ],
        "stages": [
            {
                "stage": 1,
                "name": "finite_rigidity_early_time_predecessor_core",
                "case_ids": list(PREDECESSOR_CORE),
                "advance_rule": (
                    "all_rows_complete_and_convergence_noise_current_and_"
                    "conservation_checks_pass_or_issue_reviewed_redesign"
                ),
            },
            {
                "stage": 2,
                "name": "onset_box_and_sensitivity_window_core",
                "case_ids": list(WINDOW_CORE),
                "advance_rule": (
                    "all_rows_complete_and_onset_plateau_resolution_box_"
                    "sensitivity_conservation_and_applicability_checks_pass"
                ),
            },
            {
                "stage": 3,
                "name": "paired_three_dimensional_replication",
                "case_ids": list(WINDOW_REPLICATION),
                "advance_rule": (
                    "both_small_large_pairs_complete_and_seed_level_box_"
                    "differences_pass_or_remove_three_dimensional_claim"
                ),
            },
        ],
        "resource_boundary": {
            "expected_case_count": len(pilot_ids),
            "maximum_attempts_per_case": MAXIMUM_ATTEMPTS_PER_CASE,
            "maximum_retries_per_case": MAXIMUM_RETRIES_PER_CASE,
            "maximum_campaign_node_hours": MAXIMUM_PILOT_NODE_HOURS,
            "maximum_fraction_of_then_unreserved_project_budget": (
                MAXIMUM_UNRESERVED_BUDGET_FRACTION
            ),
            "effective_node_hour_ceiling": (
                "min(500,0.10*then_unreserved_project_node_hours)"
            ),
            "stage_n_plus_one_submission_allowed_only_after_stage_n_pass": True,
            "launch_preparation_must_bind_measured_resource_model": True,
            "over_ceiling_disposition": (
                "versioned_reviewed_redesign_required_never_drop_required_"
                "paired_cases_after_inspection"
            ),
        },
        "case_bindings": bindings,
        "case_bindings_sha256": _canonical_sha256(bindings),
        "qualifying_seed_inventory": {
            "field_seeds": list(QUALIFYING_FIELD_SEEDS),
            "particle_seeds": list(QUALIFYING_PARTICLE_SEEDS),
            "paired": [
                [field, particle]
                for field, particle in zip(
                    QUALIFYING_FIELD_SEEDS, QUALIFYING_PARTICLE_SEEDS
                )
            ],
            "disjoint_from_every_physical_pilot": True,
            "use_before_separate_production_freeze": False,
        },
        "selection_rules": {
            "normalized_time": "tau=k0*U_A*t",
            "nonlinear_onset": {
                "definition": "first_retained_sample_at_or_above_threshold",
                "Bperp_rms_over_B0": ONSET_BPERP_RMS_OVER_B0,
            },
            "plateau": {
                "quantity": "E_Bperp_proportional_to_Bperp_rms_squared",
                "minimum_span_tau": MINIMUM_PLATEAU_SPAN_TAU,
                "maximum_span_tau": MAXIMUM_PLATEAU_SPAN_TAU,
                "minimum_samples": MINIMUM_PLATEAU_SAMPLES,
                "maximum_absolute_log_energy_slope_per_tau": (
                    MAXIMUM_ABSOLUTE_LOG_ENERGY_SLOPE_PER_TAU
                ),
                "maximum_energy_max_to_min": MAXIMUM_PLATEAU_ENERGY_MAX_TO_MIN,
                "minimum_post_window_coverage_tau": (
                    MINIMUM_POST_PLATEAU_COVERAGE_TAU
                ),
                "selection": "earliest_passing_contiguous_window_after_onset",
            },
            "production_horizon": {
                "tau": (
                    "max(1.5*plateau_end_tau,plateau_end_tau+8); round_up_to_"
                    "field_output_cadence"
                ),
                "multiplier": PRODUCTION_HORIZON_MULTIPLIER,
                "minimum_post_plateau_tau": PRODUCTION_POST_PLATEAU_COVERAGE_TAU,
            },
        },
        "numeric_gates": {
            "maximum_paired_relative_difference": (
                MAXIMUM_PAIRED_RELATIVE_DIFFERENCE
            ),
            "maximum_paired_time_difference_tau": (
                MAXIMUM_PAIRED_TIME_DIFFERENCE_TAU
            ),
            "maximum_conservation_fractional_residual": (
                MAXIMUM_CONSERVATION_FRACTIONAL_RESIDUAL
            ),
            "minimum_characteristic_rL_over_dx": MINIMUM_CHARACTERISTIC_RL_OVER_DX,
            "maximum_box_edge_power_fraction": MAXIMUM_BOX_EDGE_POWER_FRACTION,
            "minimum_high_rigidity_current_retention": (
                MINIMUM_HIGH_RIGIDITY_CURRENT_RETENTION
            ),
            "maximum_high_rigidity_CR_momentum_fractional_change": (
                MAXIMUM_HIGH_RIGIDITY_CR_MOMENTUM_CHANGE
            ),
            "maximum_high_rigidity_CR_energy_fractional_change": (
                MAXIMUM_HIGH_RIGIDITY_CR_ENERGY_CHANGE
            ),
            "maximum_absolute_local_Hall_parameter": (
                MAXIMUM_ABSOLUTE_LOCAL_HALL_PARAMETER
            ),
            "maximum_local_Bai_charge_fraction": (
                MAXIMUM_LOCAL_BAI_CHARGE_FRACTION
            ),
        },
        "threshold_classification": {
            "onset_Bperp_over_B0": "physics_definition_preregistered_by_project",
            "Hall_parameter_order_unity_context": "Bai_et_al_2015",
            "all_numeric_acceptance_thresholds": (
                "preregistered_project_numerical_and_scope_thresholds_not_"
                "universal_literature_tolerances"
            ),
        },
        "failure_policy": (
            "retain_all_attempts_and_issue_versioned_reviewed_redesign; never_"
            "change_windows_thresholds_seeds_or_claim_scope_after_inspection"
        ),
        "saturation_evidence_eligible": False,
        "authorization": dict(AUTHORIZATION),
    }


def validate_preregistration(value: object) -> dict[str, object]:
    expected = build_preregistration()
    _require(type(value) is dict and value == expected, "Q019 preregistration drifted")
    return expected


def _slope(x: Sequence[float], y: Sequence[float]) -> float:
    xbar = statistics.fmean(x)
    ybar = statistics.fmean(y)
    denominator = sum((value - xbar) ** 2 for value in x)
    _require(denominator > 0.0, "plateau fit has no time variance")
    return sum((a - xbar) * (b - ybar) for a, b in zip(x, y)) / denominator


def select_plateau_window(
    times: Sequence[float],
    bperp_rms_over_b0: Sequence[float],
    *,
    k0: float,
    u_a: float,
    field_output_dt: float,
) -> dict[str, object]:
    """Apply the frozen onset/plateau/horizon rule to one excluded trace."""
    _require(
        len(times) == len(bperp_rms_over_b0) and len(times) >= MINIMUM_PLATEAU_SAMPLES,
        "physical-window trace length is invalid",
    )
    _require(k0 > 0.0 and u_a > 0.0 and field_output_dt > 0.0, "normalization is invalid")
    _require(
        all(math.isfinite(value) for value in times)
        and all(
            math.isfinite(value) and value > 0.0 for value in bperp_rms_over_b0
        )
        and all(right > left for left, right in zip(times, times[1:])),
        "physical-window trace is nonfinite, nonpositive, or nonmonotonic",
    )
    tau = [k0 * u_a * value for value in times]
    onset = next(
        (
            index
            for index, value in enumerate(bperp_rms_over_b0)
            if value >= ONSET_BPERP_RMS_OVER_B0
        ),
        None,
    )
    _require(onset is not None, "nonlinear onset was not reached")
    log_energy = [math.log(value * value) for value in bperp_rms_over_b0]
    selected: tuple[int, int, float] | None = None
    for start in range(onset, len(times) - MINIMUM_PLATEAU_SAMPLES + 1):
        for end in range(start + MINIMUM_PLATEAU_SAMPLES - 1, len(times)):
            span = tau[end] - tau[start]
            if span < MINIMUM_PLATEAU_SPAN_TAU:
                continue
            if span > MAXIMUM_PLATEAU_SPAN_TAU:
                break
            if tau[-1] - tau[end] < MINIMUM_POST_PLATEAU_COVERAGE_TAU:
                continue
            slope = _slope(tau[start : end + 1], log_energy[start : end + 1])
            energies = [value * value for value in bperp_rms_over_b0[start : end + 1]]
            if (
                abs(slope) <= MAXIMUM_ABSOLUTE_LOG_ENERGY_SLOPE_PER_TAU
                and max(energies) / min(energies) <= MAXIMUM_PLATEAU_ENERGY_MAX_TO_MIN
            ):
                selected = (start, end, slope)
                break
        if selected is not None:
            break
    _require(selected is not None, "no preregistered sustained plateau was found")
    start, end, slope = selected
    amplitude = statistics.median(bperp_rms_over_b0[start : end + 1])
    horizon_tau = max(
        PRODUCTION_HORIZON_MULTIPLIER * tau[end],
        tau[end] + PRODUCTION_POST_PLATEAU_COVERAGE_TAU,
    )
    cadence_tau = k0 * u_a * field_output_dt
    horizon_tau = math.ceil(horizon_tau / cadence_tau) * cadence_tau
    return {
        "onset_index": onset,
        "onset_time": times[onset],
        "onset_tau": tau[onset],
        "plateau_start_index": start,
        "plateau_end_index": end,
        "plateau_start_time": times[start],
        "plateau_end_time": times[end],
        "plateau_start_tau": tau[start],
        "plateau_end_tau": tau[end],
        "plateau_log_energy_slope_per_tau": slope,
        "plateau_Bperp_rms_over_B0": amplitude,
        "recommended_production_terminal_tau": horizon_tau,
        "recommended_production_terminal_time": horizon_tau / (k0 * u_a),
        "production_freeze_authorized": False,
        "saturation_claim_authorized": False,
    }


build_preregistration()
