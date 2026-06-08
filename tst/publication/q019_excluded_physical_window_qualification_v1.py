#!/usr/bin/env python3
"""Qualify cumulative Q019 excluded physical-window pilot evidence.

The qualification recomputes every numerical decision from registered raw
analysis reports.  It can recommend advancing to the next excluded stage or
freezing a later qualifying campaign, but it grants no launch, production,
scientific-claim, or publication authority.
"""

from __future__ import annotations

import hashlib
import json
import math
from typing import Mapping, Sequence

from tst.publication import analyze_q019_physics_first_nonlinear_bell_successor_v2 as analysis
from tst.publication import q019_excluded_physical_window_preregistration_v1 as prereg
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as design


SCHEMA_VERSION = 1
RECORD_TYPE = "q019_excluded_physical_window_qualification_v1"
STATUS_PASS = "passed_excluded_physical_window_gate_non_authorizing"
STATUS_FAIL = "failed_excluded_physical_window_gate_redesign_required_non_authorizing"

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


class QualificationError(ValueError):
    """Reject malformed, incomplete, synthetic, or inventory-drifted evidence."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise QualificationError(message)


def _canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _finite(value: object, label: str) -> float:
    _require(type(value) in {int, float}, f"{label} is not numeric")
    result = float(value)
    _require(math.isfinite(result), f"{label} is nonfinite")
    return result


def _mapping(value: object, label: str) -> Mapping[str, object]:
    _require(type(value) is dict, f"{label} is not an object")
    return value


def _sequence(value: object, label: str) -> Sequence[object]:
    _require(type(value) is list and value, f"{label} is not a nonempty list")
    return value


def _case_map() -> dict[str, dict[str, object]]:
    cases = {str(case["case_id"]): case for case in design.expected_cases()}
    _require(len(cases) == 77, "Q019 design inventory drifted")
    return cases


def _expected_ids(upto_stage: int) -> tuple[str, ...]:
    _require(upto_stage in {1, 2, 3}, "Q019 physical stage is invalid")
    groups = (
        prereg.PREDECESSOR_CORE,
        prereg.WINDOW_CORE,
        prereg.WINDOW_REPLICATION,
    )
    return tuple(case_id for group in groups[:upto_stage] for case_id in group)


def _validate_report(
    report: Mapping[str, object], case: Mapping[str, object]
) -> None:
    case_id = str(case["case_id"])
    _require(
        report.get("schema_version") == 3
        and report.get("record_type") == analysis.RECORD_TYPE
        and report.get("case_id") == case_id
        and report.get("campaign_id") == case["campaign_id"]
        and report.get("branch") == case["branch"],
        f"{case_id} analysis identity drifted",
    )
    _require(
        report.get("analysis_input_kind") == "raw_registered_bundle",
        f"{case_id} is not registered raw evidence",
    )
    authority = _mapping(report.get("authority"), f"{case_id} authority")
    _require(
        authority
        == {
            "launch_authorized": False,
            "policy_authorized": False,
            "qualification_authorized": False,
            "claim_authorized": False,
            "raw_production_authorized": False,
            "nonlinear_saturation_claim_authorized": False,
        },
        f"{case_id} analysis authority drifted",
    )
    completion = _mapping(
        report.get("structured_runtime_completion_gate"),
        f"{case_id} completion gate",
    )
    _require(
        completion.get("incomplete_rejected") is False
        and completion.get("raw_science_admission_eligible") is True
        and completion.get("trusted_execution_binding_present") is True
        and completion.get("evidence_disposition")
        == "admitted_registered_raw_analysis",
        f"{case_id} registered completion was not admitted",
    )
    prerequisites = _mapping(
        report.get("independent_prerequisite_gate"),
        f"{case_id} prerequisite gate",
    )
    _require(
        prerequisites.get("q043_independent_raw_cycle_one_oracle_bound") is True
        and prerequisites.get("q023_independent_linear_predecessor_bound") is True
        and prerequisites.get("independent_prerequisites_complete") is True
        and prerequisites.get("nonlinear_execution_prerequisites_passed") is True,
        f"{case_id} independent prerequisites are incomplete",
    )


def _conservation(report: Mapping[str, object]) -> dict[str, object]:
    case_id = str(report["case_id"])
    rows = _sequence(report.get("particle_state_trace"), f"{case_id} particle trace")
    maximum_energy = 0.0
    maximum_momentum = 0.0
    for row_value in rows:
        row = _mapping(row_value, f"{case_id} particle row")
        conservation = _mapping(
            row.get("conservation"), f"{case_id} conservation row"
        )
        _require(
            conservation.get("reference_bound") is True,
            f"{case_id} conservation reference is unbound",
        )
        maximum_energy = max(
            maximum_energy,
            abs(_finite(conservation.get("energy_fractional_residual"), "energy residual")),
        )
        momentum = _sequence(
            conservation.get("momentum_fractional_residual"),
            f"{case_id} momentum residual",
        )
        _require(len(momentum) == 3, f"{case_id} momentum residual shape drifted")
        maximum_momentum = max(
            maximum_momentum,
            *(abs(_finite(value, "momentum residual")) for value in momentum),
        )
    passed = (
        max(maximum_energy, maximum_momentum)
        <= prereg.MAXIMUM_CONSERVATION_FRACTIONAL_RESIDUAL
    )
    return {
        "maximum_absolute_energy_fractional_residual": maximum_energy,
        "maximum_absolute_momentum_fractional_residual": maximum_momentum,
        "maximum_allowed_fractional_residual": (
            prereg.MAXIMUM_CONSERVATION_FRACTIONAL_RESIDUAL
        ),
        "gate_passed": passed,
    }


def _resolution(report: Mapping[str, object]) -> dict[str, object]:
    case_id = str(report["case_id"])
    gate = _mapping(
        report.get("evolving_finite_rigidity_resolution_stop_gate"),
        f"{case_id} resolution gate",
    )
    if gate.get("applicable") is False:
        return {
            "applicable": False,
            "minimum_characteristic_rl_over_dx": None,
            "gate_passed": True,
        }
    rows = _sequence(gate.get("trace"), f"{case_id} resolution trace")
    minimum = min(
        _finite(
            _mapping(row, f"{case_id} resolution row").get(
                "characteristic_median_particle_rl_over_dx"
            ),
            "characteristic rL/dx",
        )
        for row in rows
    )
    return {
        "applicable": True,
        "minimum_characteristic_rl_over_dx": minimum,
        "minimum_allowed_characteristic_rl_over_dx": (
            prereg.MINIMUM_CHARACTERISTIC_RL_OVER_DX
        ),
        "gate_passed": minimum >= prereg.MINIMUM_CHARACTERISTIC_RL_OVER_DX,
    }


def _hall(report: Mapping[str, object]) -> dict[str, object]:
    case_id = str(report["case_id"])
    gate = _mapping(
        report.get("evolving_local_hall_applicability_diagnostics"),
        f"{case_id} Hall diagnostics",
    )
    rows = _sequence(gate.get("trace"), f"{case_id} Hall trace")
    maximum_lambda = 0.0
    maximum_bai_r = 0.0
    for row_value in rows:
        row = _mapping(row_value, f"{case_id} Hall row")
        absolute_lambda = _mapping(
            row.get("absolute_local_Lambda"), f"{case_id} absolute Lambda"
        )
        bai_r = _mapping(row.get("local_Bai_R"), f"{case_id} local Bai R")
        maximum_lambda = max(
            maximum_lambda,
            _finite(absolute_lambda.get("maximum"), "absolute local Lambda"),
        )
        maximum_bai_r = max(
            maximum_bai_r, _finite(bai_r.get("maximum"), "local Bai R")
        )
    return {
        "maximum_absolute_local_Hall_parameter": maximum_lambda,
        "maximum_allowed_absolute_local_Hall_parameter": (
            prereg.MAXIMUM_ABSOLUTE_LOCAL_HALL_PARAMETER
        ),
        "maximum_local_Bai_charge_fraction": maximum_bai_r,
        "maximum_allowed_local_Bai_charge_fraction": (
            prereg.MAXIMUM_LOCAL_BAI_CHARGE_FRACTION
        ),
        "gate_passed": (
            maximum_lambda <= prereg.MAXIMUM_ABSOLUTE_LOCAL_HALL_PARAMETER
            and maximum_bai_r <= prereg.MAXIMUM_LOCAL_BAI_CHARGE_FRACTION
        ),
    }


def _current_agreement(report: Mapping[str, object]) -> dict[str, object]:
    case_id = str(report["case_id"])
    gate = _mapping(
        report.get("deposited_grid_vs_reconstructed_particle_current_gate"),
        f"{case_id} current agreement",
    )
    maximum = _finite(gate.get("maximum_relative_l2_error"), "current error")
    tolerance = _finite(gate.get("relative_l2_tolerance"), "current tolerance")
    passed = gate.get("gate_passed") is True and maximum <= tolerance
    return {
        "maximum_relative_l2_error": maximum,
        "relative_l2_tolerance": tolerance,
        "gate_passed": passed,
    }


def _noise(report: Mapping[str, object], case: Mapping[str, object]) -> dict[str, object]:
    case_id = str(report["case_id"])
    gate = _mapping(
        report.get("initial_rho_jx_noise_pair_gate"), f"{case_id} noise gate"
    )
    measurement = _mapping(
        gate.get("first_snapshot_measurement"), f"{case_id} noise measurement"
    )
    names = (
        "deposited_charge_density_coefficient_of_variation",
        "jx_over_c_coefficient_of_variation",
        "transverse_current_rms_over_abs_parallel_mean",
    )
    maximum = max(_finite(measurement.get(name), name) for name in names)
    quiet = case["finite_sampling_mode"] != "independent_position_noise_seeded"
    passed = (
        maximum <= prereg.MAXIMUM_QUIET_START_RELATIVE_NOISE
        if quiet
        else maximum >= prereg.MINIMUM_NOISE_CONTROL_RELATIVE_NOISE
    )
    return {
        "sampling_class": "quiet_packet" if quiet else "independent_position_noise",
        "maximum_selected_relative_noise": maximum,
        "maximum_allowed_quiet_start_relative_noise": (
            prereg.MAXIMUM_QUIET_START_RELATIVE_NOISE
        ),
        "minimum_required_noise_control_relative_noise": (
            prereg.MINIMUM_NOISE_CONTROL_RELATIVE_NOISE
        ),
        "gate_passed": passed,
    }


def _complex_mode(value: object, label: str) -> complex:
    row = _mapping(value, label)
    real = _finite(row.get("real"), f"{label} real")
    imag = _finite(row.get("imag"), f"{label} imag")
    amplitude = _finite(row.get("amplitude"), f"{label} amplitude")
    result = complex(real, imag)
    _require(
        math.isclose(abs(result), amplitude, rel_tol=1.0e-12, abs_tol=1.0e-15),
        f"{label} amplitude drifted",
    )
    return result


def _predecessor_fit(
    report: Mapping[str, object], case: Mapping[str, object]
) -> dict[str, object]:
    case_id = str(report["case_id"])
    measured = _mapping(
        report.get("finite_rigidity_early_time_physics_predecessor"),
        f"{case_id} predecessor measurements",
    )
    modes = [
        _complex_mode(value, f"{case_id} Bperp mode")
        for value in _sequence(
            measured.get("signed_complex_Bperp_k0_trace"),
            f"{case_id} Bperp mode trace",
        )
    ]
    chronology = _mapping(
        report.get("deposited_particle_moment_chronology"),
        f"{case_id} moment chronology",
    )
    rows = _sequence(chronology.get("chronology"), f"{case_id} chronology")
    times = [
        _finite(_mapping(row, f"{case_id} chronology row").get("time"), "time")
        for row in rows
        if _mapping(
            _mapping(row, f"{case_id} chronology row").get(
                "deposited_particle_moment_state"
            ),
            f"{case_id} moment state",
        ).get("available")
        is True
    ]
    _require(len(times) == len(modes), f"{case_id} fit chronology drifted")
    try:
        selected = prereg.select_early_time_fit_window(
            times,
            modes,
            k0=float(case["k0"]),
            u_a=float(case["u_a"]),
            b0=float(case["b_g"]),
        )
    except prereg.PreregistrationError as error:
        return {
            "selection_passed": False,
            "failure_reason": str(error),
            "fit": None,
        }
    return {"selection_passed": True, "failure_reason": None, "fit": selected}


def _plateau(
    report: Mapping[str, object], case: Mapping[str, object]
) -> dict[str, object]:
    case_id = str(report["case_id"])
    if case["branch"] == "high_rigidity_current_retention_candidate":
        gate = _mapping(
            report.get("high_rigidity_saturation_mechanism_discriminants"),
            f"{case_id} high-rigidity discriminants",
        )
    else:
        gate = _mapping(
            report.get("finite_nonlinear_onset_and_saturation_prerequisite_gate"),
            f"{case_id} finite nonlinear trace",
        )
    rows = _sequence(gate.get("trace"), f"{case_id} nonlinear trace")
    times = [
        _finite(_mapping(row, f"{case_id} nonlinear row").get("time"), "time")
        for row in rows
    ]
    amplitudes = [
        _finite(
            _mapping(row, f"{case_id} nonlinear row").get("Bperp_rms_over_B0"),
            "Bperp_rms/B0",
        )
        for row in rows
    ]
    try:
        selected = prereg.select_plateau_window(
            times,
            amplitudes,
            k0=float(case["k0"]),
            u_a=float(case["u_a"]),
            field_output_dt=design.FIELD_OUTPUT_DT,
        )
    except prereg.PreregistrationError as error:
        return {
            "selection_passed": False,
            "failure_reason": str(error),
            "window": None,
        }
    return {"selection_passed": True, "failure_reason": None, "window": selected}


def _box_edge(
    report: Mapping[str, object], plateau: Mapping[str, object] | None
) -> dict[str, object]:
    case_id = str(report["case_id"])
    if plateau is None:
        return {
            "maximum_plateau_box_edge_power_fraction": None,
            "maximum_allowed_box_edge_power_fraction": (
                prereg.MAXIMUM_BOX_EDGE_POWER_FRACTION
            ),
            "gate_passed": False,
        }
    start = _finite(plateau.get("plateau_start_time"), "plateau start time")
    end = _finite(plateau.get("plateau_end_time"), "plateau end time")
    rows = [
        _mapping(value, f"{case_id} spectrum row")
        for value in _sequence(
            report.get("nested_spectrum_trace"), f"{case_id} spectrum trace"
        )
    ]
    selected = [
        row
        for row in rows
        if start
        <= _finite(row.get("time"), "spectrum time")
        <= end
    ]
    _require(selected, f"{case_id} has no spectrum samples in selected plateau")
    maximum = max(
        _finite(
            row.get("large_box_only_unavailable_mode_power_fraction"),
            "box-edge power fraction",
        )
        for row in selected
    )
    return {
        "maximum_plateau_box_edge_power_fraction": maximum,
        "maximum_allowed_box_edge_power_fraction": (
            prereg.MAXIMUM_BOX_EDGE_POWER_FRACTION
        ),
        "gate_passed": maximum <= prereg.MAXIMUM_BOX_EDGE_POWER_FRACTION,
    }


def _high_rigidity(report: Mapping[str, object]) -> dict[str, object]:
    case_id = str(report["case_id"])
    gate = _mapping(
        report.get("high_rigidity_fixed_current_like_validity_gate"),
        f"{case_id} high-rigidity validity",
    )
    if gate.get("applicable") is False:
        return {"applicable": False, "gate_passed": True}
    rows = _sequence(gate.get("trace"), f"{case_id} high-rigidity trace")
    minimum_current = min(
        _finite(
            _mapping(row, f"{case_id} high-rigidity row").get(
                "lab_parallel_current_retention"
            ),
            "current retention",
        )
        for row in rows
    )
    maximum_momentum = max(
        abs(
            _finite(
                _mapping(row, f"{case_id} high-rigidity row").get(
                    "CR_parallel_momentum_fractional_change"
                ),
                "CR momentum change",
            )
        )
        for row in rows
    )
    maximum_energy = max(
        abs(
            _finite(
                _mapping(row, f"{case_id} high-rigidity row").get(
                    "CR_kinetic_energy_fractional_change"
                ),
                "CR energy change",
            )
        )
        for row in rows
    )
    return {
        "applicable": True,
        "minimum_current_retention": minimum_current,
        "minimum_allowed_current_retention": (
            prereg.MINIMUM_HIGH_RIGIDITY_CURRENT_RETENTION
        ),
        "maximum_absolute_CR_momentum_fractional_change": maximum_momentum,
        "maximum_allowed_absolute_CR_momentum_fractional_change": (
            prereg.MAXIMUM_HIGH_RIGIDITY_CR_MOMENTUM_CHANGE
        ),
        "maximum_absolute_CR_energy_fractional_change": maximum_energy,
        "maximum_allowed_absolute_CR_energy_fractional_change": (
            prereg.MAXIMUM_HIGH_RIGIDITY_CR_ENERGY_CHANGE
        ),
        "gate_passed": (
            minimum_current >= prereg.MINIMUM_HIGH_RIGIDITY_CURRENT_RETENTION
            and maximum_momentum <= prereg.MAXIMUM_HIGH_RIGIDITY_CR_MOMENTUM_CHANGE
            and maximum_energy <= prereg.MAXIMUM_HIGH_RIGIDITY_CR_ENERGY_CHANGE
        ),
    }


def _relative_difference(left: float, right: float) -> float:
    scale = max(abs(left), abs(right))
    return 0.0 if scale == 0.0 else abs(left - right) / scale


def _predecessor_comparison(
    reference_id: str,
    control_id: str,
    fits: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    reference = _mapping(fits[reference_id].get("fit"), "reference fit")
    control = _mapping(fits[control_id].get("fit"), "control fit")
    growth_difference = _relative_difference(
        _finite(reference.get("normalized_growth_rate"), "reference growth"),
        _finite(control.get("normalized_growth_rate"), "control growth"),
    )
    frequency_difference = abs(
        _finite(reference.get("normalized_angular_frequency"), "reference frequency")
        - _finite(control.get("normalized_angular_frequency"), "control frequency")
    )
    passed = (
        growth_difference <= prereg.MAXIMUM_PAIRED_RELATIVE_DIFFERENCE
        and frequency_difference <= prereg.MAXIMUM_PAIRED_RELATIVE_DIFFERENCE
    )
    return {
        "reference_case_id": reference_id,
        "control_case_id": control_id,
        "normalized_growth_rate_relative_difference": growth_difference,
        "normalized_angular_frequency_absolute_difference": frequency_difference,
        "maximum_allowed_difference": prereg.MAXIMUM_PAIRED_RELATIVE_DIFFERENCE,
        "gate_passed": passed,
    }


def _nonlinear_comparison(
    left_id: str,
    right_id: str,
    plateaus: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    left = _mapping(plateaus[left_id].get("window"), "left plateau")
    right = _mapping(plateaus[right_id].get("window"), "right plateau")
    amplitude_difference = _relative_difference(
        _finite(left.get("plateau_Bperp_rms_over_B0"), "left plateau amplitude"),
        _finite(right.get("plateau_Bperp_rms_over_B0"), "right plateau amplitude"),
    )
    onset_difference = abs(
        _finite(left.get("onset_tau"), "left onset tau")
        - _finite(right.get("onset_tau"), "right onset tau")
    )
    end_difference = abs(
        _finite(left.get("plateau_end_tau"), "left plateau end tau")
        - _finite(right.get("plateau_end_tau"), "right plateau end tau")
    )
    return {
        "left_case_id": left_id,
        "right_case_id": right_id,
        "plateau_amplitude_relative_difference": amplitude_difference,
        "onset_absolute_difference_tau": onset_difference,
        "plateau_end_absolute_difference_tau": end_difference,
        "maximum_allowed_relative_difference": (
            prereg.MAXIMUM_PAIRED_RELATIVE_DIFFERENCE
        ),
        "maximum_allowed_time_difference_tau": (
            prereg.MAXIMUM_PAIRED_TIME_DIFFERENCE_TAU
        ),
        "gate_passed": (
            amplitude_difference <= prereg.MAXIMUM_PAIRED_RELATIVE_DIFFERENCE
            and onset_difference <= prereg.MAXIMUM_PAIRED_TIME_DIFFERENCE_TAU
            and end_difference <= prereg.MAXIMUM_PAIRED_TIME_DIFFERENCE_TAU
        ),
    }


def _replication_group(
    case_ids: Sequence[str],
    plateaus: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    comparisons = [
        _nonlinear_comparison(left, right, plateaus)
        for index, left in enumerate(case_ids)
        for right in case_ids[index + 1 :]
    ]
    return {
        "case_ids": list(case_ids),
        "pairwise_comparisons": comparisons,
        "gate_passed": all(item["gate_passed"] for item in comparisons),
    }


def build_qualification(
    reports: Sequence[Mapping[str, object]], *, upto_stage: int
) -> dict[str, object]:
    """Build a cumulative stage qualification from exact registered reports."""
    preregistration = prereg.build_preregistration()
    expected_ids = _expected_ids(upto_stage)
    _require(len(reports) == len(expected_ids), "Q019 report count drifted")
    by_id: dict[str, Mapping[str, object]] = {}
    for report in reports:
        _require(type(report) is dict, "Q019 report is not an object")
        case_id = report.get("case_id")
        _require(type(case_id) is str and case_id not in by_id, "duplicate report")
        by_id[case_id] = report
    _require(set(by_id) == set(expected_ids), "Q019 report inventory drifted")

    cases = _case_map()
    report_bindings = []
    common = {}
    for case_id in expected_ids:
        report = by_id[case_id]
        case = cases[case_id]
        _validate_report(report, case)
        conservation = _conservation(report)
        resolution = _resolution(report)
        hall = _hall(report)
        current = _current_agreement(report)
        common_pass = all(
            gate["gate_passed"]
            for gate in (conservation, resolution, hall, current)
        )
        common[case_id] = {
            "case_id": case_id,
            "conservation": conservation,
            "finite_rigidity_resolution": resolution,
            "local_Hall_and_Bai_R": hall,
            "deposited_vs_reconstructed_current": current,
            "common_gate_passed": common_pass,
        }
        report_bindings.append(
            {
                "case_id": case_id,
                "canonical_analysis_report_sha256": _canonical_sha256(report),
            }
        )

    predecessor_fits = {
        case_id: _predecessor_fit(by_id[case_id], cases[case_id])
        for case_id in prereg.PREDECESSOR_CORE
        if case_id in by_id
    }
    predecessor_noise = {
        case_id: _noise(by_id[case_id], cases[case_id])
        for case_id in prereg.PREDECESSOR_CORE
        if case_id in by_id
    }
    predecessor_comparisons = []
    if predecessor_fits:
        reference_id = prereg.PREDECESSOR_CORE[0]
        if all(item["selection_passed"] for item in predecessor_fits.values()):
            predecessor_comparisons = [
                _predecessor_comparison(reference_id, control_id, predecessor_fits)
                for control_id in prereg.PREDECESSOR_CORE[1:]
            ]
    stage1_pass = (
        all(common[case_id]["common_gate_passed"] for case_id in prereg.PREDECESSOR_CORE)
        and all(item["selection_passed"] for item in predecessor_fits.values())
        and all(item["gate_passed"] for item in predecessor_noise.values())
        and len(predecessor_comparisons) == len(prereg.PREDECESSOR_CORE) - 1
        and all(item["gate_passed"] for item in predecessor_comparisons)
    )
    stage_records = [
        {
            "stage": 1,
            "name": "finite_rigidity_early_time_predecessor_core",
            "case_ids": list(prereg.PREDECESSOR_CORE),
            "case_gates": [common[case_id] for case_id in prereg.PREDECESSOR_CORE],
            "fit_selections": [
                {"case_id": case_id, **predecessor_fits[case_id]}
                for case_id in prereg.PREDECESSOR_CORE
            ],
            "noise_gates": [
                {"case_id": case_id, **predecessor_noise[case_id]}
                for case_id in prereg.PREDECESSOR_CORE
            ],
            "convergence_comparisons": predecessor_comparisons,
            "stage_gate_passed": stage1_pass,
            "next_stage_submission_recommended": stage1_pass and upto_stage == 1,
            "next_stage_submission_authorized": False,
        }
    ]

    plateaus: dict[str, dict[str, object]] = {}
    if upto_stage >= 2:
        nonlinear_ids = prereg.WINDOW_CORE
        plateaus.update(
            {
                case_id: _plateau(by_id[case_id], cases[case_id])
                for case_id in nonlinear_ids
            }
        )
        box_gates = {
            case_id: _box_edge(
                by_id[case_id],
                plateaus[case_id]["window"]
                if plateaus[case_id]["selection_passed"]
                else None,
            )
            for case_id in nonlinear_ids
        }
        high_rigidity = {
            case_id: _high_rigidity(by_id[case_id]) for case_id in nonlinear_ids
        }
        comparisons = []
        if all(item["selection_passed"] for item in plateaus.values()):
            comparisons = [
                _nonlinear_comparison(left, right, plateaus)
                for left, right in preregistration["paired_comparison_rules"][
                    "stage_2_pairs"
                ]
            ]
        stage2_pass = (
            stage1_pass
            and all(common[case_id]["common_gate_passed"] for case_id in nonlinear_ids)
            and all(item["selection_passed"] for item in plateaus.values())
            and all(item["gate_passed"] for item in box_gates.values())
            and all(item["gate_passed"] for item in high_rigidity.values())
            and len(comparisons)
            == len(
                preregistration["paired_comparison_rules"]["stage_2_pairs"]
            )
            and all(item["gate_passed"] for item in comparisons)
        )
        stage_records.append(
            {
                "stage": 2,
                "name": "onset_box_and_sensitivity_window_core",
                "case_ids": list(nonlinear_ids),
                "case_gates": [common[case_id] for case_id in nonlinear_ids],
                "plateau_selections": [
                    {"case_id": case_id, **plateaus[case_id]}
                    for case_id in nonlinear_ids
                ],
                "box_edge_gates": [
                    {"case_id": case_id, **box_gates[case_id]}
                    for case_id in nonlinear_ids
                ],
                "high_rigidity_gates": [
                    {"case_id": case_id, **high_rigidity[case_id]}
                    for case_id in nonlinear_ids
                ],
                "paired_comparisons": comparisons,
                "stage_gate_passed": stage2_pass,
                "prior_stage_gate_passed": stage1_pass,
                "next_stage_submission_recommended": stage2_pass and upto_stage == 2,
                "next_stage_submission_authorized": False,
            }
        )
    else:
        stage2_pass = False

    if upto_stage == 3:
        nonlinear_ids = prereg.WINDOW_REPLICATION
        plateaus.update(
            {
                case_id: _plateau(by_id[case_id], cases[case_id])
                for case_id in nonlinear_ids
            }
        )
        box_gates = {
            case_id: _box_edge(
                by_id[case_id],
                plateaus[case_id]["window"]
                if plateaus[case_id]["selection_passed"]
                else None,
            )
            for case_id in nonlinear_ids
        }
        pair_rules = preregistration["paired_comparison_rules"]
        box_pairs = []
        seed_groups = []
        if all(item["selection_passed"] for item in plateaus.values()):
            box_pairs = [
                _nonlinear_comparison(left, right, plateaus)
                for left, right in pair_rules["stage_3_small_large_pairs"]
            ]
            seed_groups = [
                _replication_group(group, plateaus)
                for group in pair_rules["stage_3_seed_replication_groups"]
            ]
        stage3_pass = (
            stage2_pass
            and all(common[case_id]["common_gate_passed"] for case_id in nonlinear_ids)
            and all(plateaus[case_id]["selection_passed"] for case_id in nonlinear_ids)
            and all(item["gate_passed"] for item in box_gates.values())
            and len(box_pairs) == len(pair_rules["stage_3_small_large_pairs"])
            and all(item["gate_passed"] for item in box_pairs)
            and len(seed_groups) == len(pair_rules["stage_3_seed_replication_groups"])
            and all(item["gate_passed"] for item in seed_groups)
        )
        stage_records.append(
            {
                "stage": 3,
                "name": "paired_three_dimensional_replication",
                "case_ids": list(nonlinear_ids),
                "case_gates": [common[case_id] for case_id in nonlinear_ids],
                "plateau_selections": [
                    {"case_id": case_id, **plateaus[case_id]}
                    for case_id in nonlinear_ids
                ],
                "box_edge_gates": [
                    {"case_id": case_id, **box_gates[case_id]}
                    for case_id in nonlinear_ids
                ],
                "small_large_box_comparisons": box_pairs,
                "seed_replication_groups": seed_groups,
                "stage_gate_passed": stage3_pass,
                "prior_stage_gate_passed": stage2_pass,
                "next_stage_submission_recommended": False,
                "next_stage_submission_authorized": False,
            }
        )
    else:
        stage3_pass = False

    current_stage_pass = stage_records[-1]["stage_gate_passed"]
    full_pass = upto_stage == 3 and stage3_pass
    horizons = [
        _finite(item["window"]["recommended_production_terminal_tau"], "horizon")
        for case_id, item in plateaus.items()
        if item["selection_passed"]
        and case_id.startswith("q019-fr-3d-onset-")
    ]
    decision = {
        "cumulative_stage": upto_stage,
        "current_stage_gate_passed": current_stage_pass,
        "next_excluded_stage_submission_recommended": (
            current_stage_pass and upto_stage < 3
        ),
        "next_excluded_stage_submission_authorized": False,
        "production_freeze_recommended": full_pass,
        "production_freeze_authorized": False,
        "recommended_production_terminal_tau": max(horizons) if full_pass else None,
        "recommended_field_output_dt": design.FIELD_OUTPUT_DT if full_pass else None,
        "recommended_qualifying_seed_inventory": (
            preregistration["qualifying_seed_inventory"] if full_pass else None
        ),
        "scientific_claim_authorized": False,
        "publication_authorized": False,
        "failure_disposition": (
            None
            if current_stage_pass
            else "retain_all_attempts_and_prepare_versioned_reviewed_redesign"
        ),
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "status": STATUS_PASS if current_stage_pass else STATUS_FAIL,
        "campaign": "q019_nonlinear_bell_registered_successor_v1",
        "preregistration": {
            "record_type": prereg.RECORD_TYPE,
            "canonical_sha256": _canonical_sha256(preregistration),
        },
        "cumulative_stage": upto_stage,
        "expected_case_ids": list(expected_ids),
        "analysis_report_bindings": report_bindings,
        "analysis_report_bindings_sha256": _canonical_sha256(report_bindings),
        "stage_qualifications": stage_records,
        "decision": decision,
        "excluded_pilot_saturation_evidence_eligible": False,
        "authorization": dict(AUTHORIZATION),
    }


def validate_qualification(
    value: object, reports: Sequence[Mapping[str, object]]
) -> dict[str, object]:
    _require(
        type(value) is dict and value.get("record_type") == RECORD_TYPE,
        "Q019 physical qualification identity drifted",
    )
    upto_stage = value.get("cumulative_stage")
    _require(type(upto_stage) is int, "Q019 cumulative stage is invalid")
    rebuilt = build_qualification(reports, upto_stage=upto_stage)
    _require(value == rebuilt, "Q019 physical qualification derived fields drifted")
    return rebuilt

