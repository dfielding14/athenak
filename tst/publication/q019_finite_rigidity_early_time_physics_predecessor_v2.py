#!/usr/bin/env python3
"""Physics contract for the Q019 finite-rigidity early-time predecessor.

The exact drifting-isotropic-shell MHD-PIC system does not inherit a closed-form numeric
acceptance tolerance from the cited literature. This contract therefore
defines a measurable complex-mode reference and a convergence hierarchy while
leaving fit windows and tolerances to excluded pilots and external review.
"""

from __future__ import annotations

from typing import Mapping

from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as decks


SCHEMA_VERSION = 2
RECORD_TYPE = "q019_finite_rigidity_early_time_physics_predecessor_v2"
REFERENCE_ROLE = "finite_rigidity_early_time_physics_reference_grid"
CONTROL_ROLE = "finite_rigidity_early_time_convergence_control"
REQUIRED_MEASUREMENTS = (
    "signed_complex_Bperp_k0_trace",
    "signed_complex_deposited_Jperp_k0_trace",
    "complex_exponential_growth_and_frequency_fit",
    "polarization_handedness",
    "Bperp_Jperp_phase_relation",
    "parallel_current_retention",
    "transverse_force_noise",
    "deposited_grid_vs_reconstructed_particle_current_agreement",
    "initial_rho_jx_noise_pair_gate",
    "evolving_rl_over_dx",
)
REQUIRED_CONVERGENCE_AXES = (
    "resolution",
    "PPC",
    "timestep",
    "quiet_isotropic_shell_packet_vs_independent_positions",
)


class PredecessorContractError(ValueError):
    """Raised when the finite-rigidity predecessor contract drifts."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PredecessorContractError(message)


def _cases() -> tuple[dict[str, object], ...]:
    return decks.expected_cases()


def build_contract() -> dict[str, object]:
    references = [case for case in _cases() if case["role"] == REFERENCE_ROLE]
    controls = [case for case in _cases() if case["role"] == CONTROL_ROLE]
    _require(len(references) == 9, "finite predecessor reference grid drifted")
    _require(len(controls) == 6, "finite predecessor convergence controls drifted")
    _require(
        {
            (case["k0_rg0"], case["rho_cr_over_rho0"])
            for case in references
        }
        == {
            (k0_rg0, rho_cr)
            for k0_rg0 in decks.FINITE_K0_RG0_GRID
            for rho_cr in decks.FINITE_RHO_CR_GRID
        },
        "finite predecessor physical grid drifted",
    )
    _require(
        all(case["broadband_amplitude"] == 0.0 for case in references + controls),
        "finite predecessor must isolate the seeded k0 response",
    )
    _require(
        all(int(case["ppc"]) == max(decks.FINITE_PPC_LADDER) for case in references)
        and all(tuple(case["nx"]) == (1536, 768, 1) for case in references),
        "finite predecessor reference rows are not the highest-resolution/highest-PPC rows",
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "status": "source_local_physics_contract_complete_runtime_predecessor_pending",
        "reference_definition": {
            "kind": "measured_continuum_convergence_reference_for_exact_implementation",
            "reference_row": (
                "highest_resolution_highest_PPC_quiet_isotropic_shell_packet_member_at_each_"
                "physics_point"
            ),
            "fixed_current_Bell_growth_role": (
                "asymptotic_context_only_not_finite_rigidity_acceptance_target"
            ),
            "analytic_finite_shell_dispersion_relation_claimed": False,
        },
        "required_measurements": list(REQUIRED_MEASUREMENTS),
        "required_convergence_axes": list(REQUIRED_CONVERGENCE_AXES),
        "fit_window": None,
        "numeric_acceptance_tolerances": None,
        "threshold_source": (
            "freeze_after_excluded_predecessor_pilots_and_external_review_where_"
            "literature_does_not_supply_exact_implementation_tolerances"
        ),
        "reference_case_ids": [str(case["case_id"]) for case in references],
        "control_case_ids": [str(case["case_id"]) for case in controls],
        "runtime_predecessor_complete": False,
        "physical_pilot_authorized": False,
        "authority": {
            "launch_authorized": False,
            "policy_authorized": False,
            "qualification_authorized": False,
            "claim_authorized": False,
            "raw_production_authorized": False,
            "nonlinear_saturation_claim_authorized": False,
        },
    }


def validate_analysis_report(report: Mapping[str, object]) -> dict[str, object]:
    """Validate measurement coverage without granting a pass."""
    contract = build_contract()
    _require(type(report) is dict, "analysis report must be an object")
    _require(
        report.get("branch") == "finite_rigidity_early_time_predecessor",
        "analysis report is not a finite-rigidity predecessor row",
    )
    _require(
        report.get("authority") == contract["authority"],
        "analysis report authority drifted",
    )
    measured = report.get("finite_rigidity_early_time_physics_predecessor")
    _require(type(measured) is dict, "predecessor measurement block is missing")
    _require(
        all(name in measured for name in REQUIRED_MEASUREMENTS),
        "predecessor measurement inventory is incomplete",
    )
    _require(
        measured.get("fit_window") is None
        and measured.get("numeric_acceptance_tolerances") is None
        and measured.get("gate_passed") is None
        and measured.get("runtime_predecessor_complete") is False,
        "pilot-derived predecessor boundary drifted",
    )
    return {
        "analysis_contract_observed": True,
        "runtime_predecessor_complete": False,
        "physical_pilot_authorized": False,
    }


build_contract()
