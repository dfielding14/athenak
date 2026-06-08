#!/usr/bin/env python3
"""Focused tests for cumulative Q019 excluded physical-window qualification."""

from __future__ import annotations

import copy
import math

from tst.publication import q019_excluded_physical_window_preregistration_v1 as prereg
from tst.publication import q019_excluded_physical_window_qualification_v1 as qual
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as design


def _complex_record(value: complex) -> dict[str, float]:
    return {
        "real": value.real,
        "imag": value.imag,
        "amplitude": abs(value),
        "phase": math.atan2(value.imag, value.real),
    }


def _amplitude(time: float) -> float:
    tau = 2.0 * math.pi * time
    if tau < 14.0:
        return 1.0e-3 * math.exp(0.5 * tau)
    return 2.5 * math.exp(0.01 * math.sin(tau))


def _report(case_id: str) -> dict[str, object]:
    case = next(row for row in design.expected_cases() if row["case_id"] == case_id)
    authority = {
        "launch_authorized": False,
        "policy_authorized": False,
        "qualification_authorized": False,
        "claim_authorized": False,
        "raw_production_authorized": False,
        "nonlinear_saturation_claim_authorized": False,
    }
    conservation = {
        "reference_bound": True,
        "energy_fractional_residual": 1.0e-8,
        "momentum_fractional_residual": [1.0e-8, 0.0, 0.0],
    }
    report: dict[str, object] = {
        "schema_version": 3,
        "record_type": "q019_physics_first_nonlinear_bell_successor_v2_analysis",
        "case_id": case_id,
        "campaign_id": case["campaign_id"],
        "branch": case["branch"],
        "analysis_input_kind": "raw_registered_bundle",
        "authority": authority,
        "structured_runtime_completion_gate": {
            "run_completion_status": "complete_terminal_time_reached",
            "incomplete_rejected": False,
            "evidence_disposition": "admitted_registered_raw_analysis",
            "saturation_evidence_eligible": False,
            "raw_science_admission_eligible": True,
            "trusted_execution_binding_present": True,
        },
        "independent_prerequisite_gate": {
            "q043_independent_raw_cycle_one_oracle_bound": True,
            "q023_independent_linear_predecessor_bound": True,
            "independent_prerequisites_complete": True,
            "nonlinear_execution_prerequisites_passed": True,
        },
        "particle_state_trace": [{"conservation": conservation}],
        "evolving_local_hall_applicability_diagnostics": {
            "trace": [
                {
                    "absolute_local_Lambda": {"maximum": 0.01},
                    "local_Bai_R": {"maximum": 1.0e-4},
                }
            ]
        },
        "deposited_grid_vs_reconstructed_particle_current_gate": {
            "maximum_relative_l2_error": 1.0e-10,
            "relative_l2_tolerance": 1.0e-8,
            "gate_passed": True,
        },
        "initial_rho_jx_noise_pair_gate": {
            "first_snapshot_measurement": {
                "deposited_charge_density_coefficient_of_variation": 0.0,
                "jx_over_c_coefficient_of_variation": 0.0,
                "transverse_current_rms_over_abs_parallel_mean": 0.0,
            }
        },
    }
    if case["branch"] == "high_rigidity_current_retention_candidate":
        report["evolving_finite_rigidity_resolution_stop_gate"] = {
            "applicable": False
        }
        report["high_rigidity_fixed_current_like_validity_gate"] = {
            "applicable": True,
            "trace": [
                {
                    "lab_parallel_current_retention": 0.99,
                    "CR_parallel_momentum_fractional_change": 0.01,
                    "CR_kinetic_energy_fractional_change": 0.01,
                }
            ],
        }
    else:
        report["evolving_finite_rigidity_resolution_stop_gate"] = {
            "applicable": True,
            "trace": [
                {"characteristic_median_particle_rl_over_dx": 12.0}
            ],
        }
        report["high_rigidity_fixed_current_like_validity_gate"] = {
            "applicable": False
        }

    if case["branch"] == "finite_rigidity_early_time_predecessor":
        times = [0.1 * index for index in range(1, 16)]
        modes = [
            1.0e-7
            * math.exp(0.8 * 2.0 * math.pi * time)
            * complex(
                math.cos(-0.2 * 2.0 * math.pi * time),
                math.sin(-0.2 * 2.0 * math.pi * time),
            )
            for time in times
        ]
        report["finite_rigidity_early_time_physics_predecessor"] = {
            "signed_complex_Bperp_k0_trace": [
                _complex_record(value) for value in modes
            ]
        }
        report["deposited_particle_moment_chronology"] = {
            "chronology": [
                {
                    "time": time,
                    "deposited_particle_moment_state": {"available": True},
                }
                for time in times
            ]
        }
        if case["finite_sampling_mode"] == "independent_position_noise_seeded":
            measurement = report["initial_rho_jx_noise_pair_gate"][
                "first_snapshot_measurement"
            ]
            measurement["deposited_charge_density_coefficient_of_variation"] = 1.0e-3
    else:
        times = [0.1 * index for index in range(1, 121)]
        trace = [
            {
                "time": time,
                "Bperp_rms_over_B0": _amplitude(time),
            }
            for time in times
        ]
        if case["branch"] == "high_rigidity_current_retention_candidate":
            report["high_rigidity_saturation_mechanism_discriminants"] = {
                "trace": trace
            }
        else:
            report["finite_nonlinear_onset_and_saturation_prerequisite_gate"] = {
                "trace": trace
            }
        report["nested_spectrum_trace"] = [
            {
                "time": time,
                "large_box_only_unavailable_mode_power_fraction": 0.02,
            }
            for time in times
        ]
    return report


def _reports(upto_stage: int) -> list[dict[str, object]]:
    groups = (
        prereg.PREDECESSOR_CORE,
        prereg.WINDOW_CORE,
        prereg.WINDOW_REPLICATION,
    )
    return [
        _report(case_id)
        for group in groups[:upto_stage]
        for case_id in group
    ]


def test_stage_one_passes_and_only_recommends_stage_two() -> None:
    result = qual.build_qualification(_reports(1), upto_stage=1)
    assert result["status"] == qual.STATUS_PASS
    assert result["decision"]["next_excluded_stage_submission_recommended"] is True
    assert result["decision"]["next_excluded_stage_submission_authorized"] is False
    assert result["decision"]["production_freeze_recommended"] is False
    assert not any(result["authorization"].values())


def test_full_three_stage_pass_recommends_conservative_production_freeze() -> None:
    result = qual.build_qualification(_reports(3), upto_stage=3)
    assert result["status"] == qual.STATUS_PASS
    assert result["decision"]["production_freeze_recommended"] is True
    assert result["decision"]["production_freeze_authorized"] is False
    assert result["decision"]["recommended_production_terminal_tau"] > 0.0
    assert (
        result["decision"]["recommended_qualifying_seed_inventory"]
        == prereg.build_preregistration()["qualifying_seed_inventory"]
    )
    assert result["excluded_pilot_saturation_evidence_eligible"] is False


def test_synthetic_or_incomplete_report_is_rejected_as_non_evidence() -> None:
    reports = _reports(1)
    reports[0]["analysis_input_kind"] = "synthetic_contract_fixture"
    try:
        qual.build_qualification(reports, upto_stage=1)
    except qual.QualificationError as error:
        assert "registered raw" in str(error)
    else:
        raise AssertionError("synthetic evidence was accepted")


def test_missing_onset_is_a_failed_physics_gate_not_an_authority_change() -> None:
    reports = _reports(2)
    target = next(
        row for row in reports if row["case_id"] == "q019-fr-3d-onset-small-s0"
    )
    for item in target["finite_nonlinear_onset_and_saturation_prerequisite_gate"][
        "trace"
    ]:
        item["Bperp_rms_over_B0"] = 0.1
    result = qual.build_qualification(reports, upto_stage=2)
    assert result["status"] == qual.STATUS_FAIL
    assert result["stage_qualifications"][1]["stage_gate_passed"] is False
    assert result["decision"]["failure_disposition"] is not None
    assert not any(result["authorization"].values())


def test_conservation_failure_is_independent_of_window_selection() -> None:
    reports = _reports(1)
    reports[0]["particle_state_trace"][0]["conservation"][
        "energy_fractional_residual"
    ] = 2.0e-3
    result = qual.build_qualification(reports, upto_stage=1)
    assert result["status"] == qual.STATUS_FAIL
    gate = result["stage_qualifications"][0]["case_gates"][0]["conservation"]
    assert gate["gate_passed"] is False


def test_validator_recomputes_all_derived_fields() -> None:
    reports = _reports(1)
    result = qual.build_qualification(reports, upto_stage=1)
    drifted = copy.deepcopy(result)
    drifted["decision"]["production_freeze_authorized"] = True
    try:
        qual.validate_qualification(drifted, reports)
    except qual.QualificationError as error:
        assert "derived fields" in str(error)
    else:
        raise AssertionError("authority-bearing mutation was accepted")

