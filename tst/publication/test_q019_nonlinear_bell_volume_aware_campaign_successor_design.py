#!/usr/bin/env python3
"""Adversarial checks for the corrected Q019 nonlinear Bell campaign design."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
RECORD_PATH = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q019_nonlinear_bell_volume_aware_campaign_successor_design_2026-06-06.json"
)
NOTE_PATH = RECORD_PATH.with_suffix(".md")

EXACT_GENERAL_CLOSURE = (
    "PPC*deposit_qscale*species_charge*v_CR/V_root_cell=J_CR/c=2*B0*k0"
)
UNIT_MASS_SPECIALIZED_CLOSURE = (
    "PPC*deposit_qscale*(q/(mc))*v_CR/V_root_cell=J_CR/c=2*B0*k0"
)
HIGH_RIGIDITY_ID = "Q019-HR-JOVERC-FIXED-CURRENT-LIKE-NOHALL"
FINITE_RIGIDITY_ID = "Q019-FR-JOVERC-SELF-CONSISTENT-NOHALL"
EXECUTABLE_LINEAGE = "q023_paper_bell_linear_joverc"
ORACLE_ID = "Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE"


class ContractError(ValueError):
    """The successor design was weakened or drifted."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_design(record: dict[str, object]) -> None:
    _require(
        record["record_type"]
        == "q019_nonlinear_bell_volume_aware_campaign_successor_design",
        "record type drifted",
    )
    _require(record["schema_version"] == 1, "schema version drifted")
    _require(
        record["qualification_effect"] == "none_source_local_non_authorizing_design_only",
        "qualification boundary drifted",
    )
    _require(
        record["authority"] and not any(record["authority"].values()),
        "authority must remain false",
    )

    lineage = record["authoritative_corrected_lineage"]
    _require(
        lineage["adversarial_audit_decision"]
        == "Q019_must_depend_first_on_the_passed_Q043_raw_cycle_one_oracle_and_"
        "then_on_a_separately_corrected_Q023_linear_Bell_predecessor",
        "authoritative predecessor order drifted",
    )
    first = lineage["first_predecessor"]
    _require(
        first["oracle_id"] == ORACLE_ID
        and first["required_snapshot"] == "raw_cycle_one"
        and first["required_result"] == "passed"
        and not first["q043_is_an_executable_generator"],
        "passed Q043 raw cycle-one predecessor drifted",
    )
    second = lineage["second_predecessor"]
    _require(
        second["campaign_id"] == "Q023-PAPER-BELL-LINEAR-JOVERC"
        and second["generator_name"] == EXECUTABLE_LINEAGE
        and "passed_Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE_raw_cycle_one_oracle"
        == second["required_after"],
        "separately corrected Q023 predecessor drifted",
    )
    _require(
        lineage["Q019_pgen_name"] == EXECUTABLE_LINEAGE,
        "Q019 pgen no longer binds the authoritative executable lineage",
    )
    _require(
        lineage["Q023_source_local_or_cycle_zero_oracle_as_Q043_substitute"]
        == "prohibited",
        "Q023 source-local or cycle-zero oracle became a Q043 substitute",
    )
    _require(
        "exact_passed_Q043_raw_cycle_one_oracle_receipt"
        in lineage["Q019_binding_rule"]
        and "separately_corrected_Q023_linear_predecessor"
        in lineage["Q019_binding_rule"],
        "Q019 authoritative-lineage binding rule drifted",
    )

    supersession = record["supersession"]
    _require(
        "all_physical_execution_qualification_and_publication_evidence"
        in supersession["effect"],
        "historical physical-use supersession drifted",
    )
    _require(
        set(supersession["invalid_inherited_contracts"])
        >= {
            "PPC_times_deposit_qscale_equals_2e9",
            "J_CR_over_2_B0_C_k0",
            "fixed_deposit_qscale_across_dimension_or_resolution",
            "legacy_q023_nonlinear_foundation_as_a_physical_predecessor",
        },
        "invalid inherited contracts were weakened",
    )
    for artifact in supersession["historical_artifacts"]:
        _require(
            _sha256(REPO_ROOT / artifact["path"]) == artifact["sha256"],
            f"historical artifact drifted: {artifact['path']}",
        )

    normalization = record["normalization_contract"]
    _require(
        normalization["exact_general_closure"] == EXACT_GENERAL_CLOSURE,
        "general closure drifted",
    )
    _require(
        normalization["deposited_quantity"] == "raw_prtcl_j_represents_J_CR_over_c",
        "deposited current semantics drifted",
    )
    _require(
        "V_root_cell" in normalization["general_deposit_qscale_derivation"],
        "qscale derivation lost root-cell volume",
    )
    specialization = normalization["q_over_mc_specialization"]
    _require(
        specialization["allowed_only_if"] == "species_mass_equals_1_is_explicitly_bound",
        "q/(mc) unit-mass specialization boundary drifted",
    )
    _require(
        specialization["specialized_closure"] == UNIT_MASS_SPECIALIZED_CLOSURE
        and not specialization["bare_q_over_mc_in_general_closure_allowed"],
        "bare q/(mc) was admitted into the general closure",
    )
    _require(
        not normalization["artificial_light_speed_in_target"]
        and not normalization["artificial_light_speed_in_qscale_derivation"],
        "artificial C returned to the current target",
    )
    _require(
        not normalization["fixed_qscale_across_dimension_resolution_or_PPC"],
        "fixed qscale was reintroduced",
    )

    bindings = record["source_local_evidence_bindings"]
    oracle = bindings["authoritative_Q043_raw_cycle_one_oracle_dependency"]
    _require(
        oracle["oracle_id"] == ORACLE_ID
        and oracle["required_snapshot"] == "raw_cycle_one"
        and oracle["required_result"] == "passed"
        and not oracle["Q023_source_local_or_cycle_zero_oracle_substitution_allowed"],
        "authoritative Q043 oracle dependency drifted",
    )

    predecessors = record["mandatory_predecessor_gates"]
    q043 = predecessors["Q043_raw_cycle_one_current_oracle"]
    _require(
        q043["required_before"]
        == "any_corrected_nonlinear_physical_pilot_or_production_execution"
        and q043["required_oracle_id"] == ORACLE_ID
        and q043["required_snapshot"] == "raw_cycle_one"
        and q043["required_result"] == "passed",
        "Q043 raw cycle-one oracle no longer blocks nonlinear pilots",
    )
    linear = predecessors["corrected_linear_predecessor"]
    _require(
        linear["required_campaign_id"] == "Q023-PAPER-BELL-LINEAR-JOVERC",
        "corrected linear predecessor identity drifted",
    )
    _require(
        linear["required_after"]
        == "passed_Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE_raw_cycle_one_oracle",
        "corrected Q023 predecessor no longer follows the Q043 oracle",
    )
    _require(
        linear["required_general_closure"] == EXACT_GENERAL_CLOSURE,
        "linear general closure drifted",
    )
    _require(
        not linear["source_local_preparation_alone_is_sufficient"],
        "source-local linear preparation falsely became sufficient",
    )
    _require(
        predecessors["gate_order"][:2]
        == [
            "passed_Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE_raw_cycle_one_oracle",
            "separately_corrected_Q023_linear_Bell_predecessor",
        ],
        "Q043-then-Q023 predecessor order drifted",
    )

    branches = record["branch_separation"]
    _require(
        branches["branches_must_not_share_a_campaign_id_or_be_pooled_statistically"],
        "branch separation was weakened",
    )
    high = branches["high_rigidity_fixed_current_like"]
    finite = branches["finite_rigidity_self_consistent"]
    _require(high["campaign_id"] == HIGH_RIGIDITY_ID, "high-rigidity ID drifted")
    _require(finite["campaign_id"] == FINITE_RIGIDITY_ID, "finite-rigidity ID drifted")
    _require(high["campaign_id"] != finite["campaign_id"], "branch IDs collided")
    centered = high["required_centered_beam_contract"]
    _require(centered["production_PPC"] == 1, "centered-beam production PPC drifted")
    _require(
        centered["PPC_greater_than_1_role"]
        == "output_level_deposition_oracle_and_debug_only_not_nonlinear_science_convergence",
        "duplicated centered particles became a convergence claim",
    )
    _require(
        not centered["PPC_convergence_claim_authorized"],
        "centered-beam PPC convergence was authorized",
    )
    _require(
        finite["PPC_role"]
        == "sampling_convergence_is_meaningful_and_requires_a_preregistered_PPC_ladder",
        "finite-rigidity sampling convergence drifted",
    )
    _require(
        "branch_specific_binding_to_the_passed_Q043_raw_cycle_one_oracle_and_"
        "then_a_reviewed_corrected_linear_or_early_growth_predecessor"
        in finite["required_new_source_features"],
        "finite-rigidity branch no longer binds the ordered predecessors",
    )

    geometry = record["geometry_and_box_contract"]
    pair = geometry["3D_required_pair"]
    _require(
        "L1_greater_than_L2_equal_L3"
        in pair["elongated_axis_aligned_fiducial"]["geometry"],
        "elongated axis-aligned 3D geometry drifted",
    )
    _require(
        "each_domain_extent_doubled" in pair["box_control"]["geometry"],
        "3D box control was weakened",
    )
    _require(
        pair["paired_seed_required"] and pair["same_cell_size_required"],
        "3D pairing was weakened",
    )
    _require(
        not geometry["3D_claim_without_complete_box_pair_authorized"],
        "unpaired 3D claim was authorized",
    )

    sensitivities = record["numerical_sensitivity_contract"]
    _require(
        {row["id"] for row in sensitivities["required_rows"]}
        == {"S-RIEMANN", "S-RECON", "S-RESOLUTION", "S-TIMESTEP"},
        "required numerical sensitivities drifted",
    )
    _require(
        any(row["change"] == "llf_to_hlld_only" for row in sensitivities["required_rows"]),
        "Riemann-solver sensitivity drifted",
    )
    _require(
        any(
            row["change"] == "plm_to_wenoz_only_with_required_ghost_zone_change"
            for row in sensitivities["required_rows"]
        ),
        "reconstruction sensitivity drifted",
    )

    diagnostics = record["required_diagnostics"]
    _require(
        "J_CR_gas=J_CR_lab-rho_CR*u_gas" in diagnostics["current_and_particle_state"],
        "gas-frame current was lost",
    )
    energy_residuals = set(diagnostics["conservation_residuals"]["energy"])
    _require(
        {
            "residual_over_absolute_cumulative_transferred_energy",
            "residual_over_absolute_MHD_energy_change",
            "residual_over_magnetic_energy_gain",
        }
        <= energy_residuals,
        "separate energy-residual normalizations were weakened",
    )
    _require(
        not diagnostics["conservation_residuals"]["single_denominator_only_is_sufficient"],
        "single-denominator residuals became sufficient",
    )

    seeds = record["seed_and_ensemble_contract"]
    _require(
        seeds["pilot_and_qualifying_seeds_must_be_disjoint"],
        "pilot and qualifying seeds can overlap",
    )
    _require(
        seeds["same_seed_paired_across_required_sensitivities_and_box_controls"],
        "paired-seed sensitivity contract drifted",
    )
    pilots = record["pilot_and_freeze_gates"]
    _require(pilots["pilot_outputs_are_never_qualifying_evidence"], "pilot scope drifted")
    _require(
        [stage["stage"] for stage in pilots["stages"]] == [0, 1, 2, 3, 4],
        "pilot/freeze stage order drifted",
    )
    _require(
        [stage["name"] for stage in pilots["stages"][:2]]
        == [
            "passed_Q043_raw_cycle_one_oracle_binding",
            "separately_corrected_Q023_linear_predecessor",
        ],
        "pilot/freeze predecessor order drifted",
    )
    _require(
        not any(stage["physical_execution_authorized_by_this_design"] for stage in pilots["stages"]),
        "a pilot stage gained execution authority",
    )

    resources = record["resource_contract"]
    _require(resources["project_wide_frontier_hard_cap_node_hours"] == 10000.0,
             "project budget cap drifted")
    _require(not resources["historical_Q019_cost_estimates_reusable"],
             "invalid historical resource estimates became reusable")
    _require(not resources["source_local_design_allocates_node_hours"],
             "source-local design allocated node hours")
    _require("never_silently_drop_seeds_sensitivities_or_the_3D_box_control"
             in resources["overrun_rule"], "resource overrun rule weakened")

    acceptance = set(record["fail_closed_acceptance"])
    _require(
        "corrected_linear_predecessor_closed_before_nonlinear_execution" in acceptance,
        "linear predecessor no longer blocks nonlinear execution",
    )
    _require(
        "passed_Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE_raw_cycle_one_oracle_"
        "receipt_and_artifacts_bound"
        in acceptance,
        "passed Q043 raw cycle-one acceptance gate was weakened",
    )
    _require(
        "required_resolution_timestep_solver_and_reconstruction_sensitivities_pass"
        in acceptance,
        "sensitivity acceptance gate was weakened",
    )
    _require(
        "3D_claim_requires_the_complete_paired_elongated_box_and_doubled_extent_control"
        in acceptance,
        "3D box acceptance gate was weakened",
    )


class Q019NonlinearBellVolumeAwareCampaignSuccessorDesignTests(unittest.TestCase):
    def setUp(self) -> None:
        self.record = json.loads(RECORD_PATH.read_text(encoding="utf-8"))

    def test_design_is_complete_non_authorizing_and_fail_closed(self) -> None:
        validate_design(self.record)

    def test_note_states_the_exact_closure_and_scientific_boundaries(self) -> None:
        note = NOTE_PATH.read_text(encoding="utf-8")
        for required in (
            "PPC * deposit_qscale * species_charge * v_CR / V_root_cell",
            "= J_CR/c",
            "= 2 B0 k0",
            "PPC=1",
            "J_CR,gas = J_CR,lab - rho_CR u_gas",
            "elongated, axis-aligned fiducial",
            "LLF -> HLLD",
            "PLM -> WENOZ",
            "10,000 node-hours",
            "not an exact locked-current experiment",
            "passed raw cycle-one",
            "Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE",
            "Q-023 source-local or",
            "cycle-zero oracle is not a substitute",
            "Q023-PAPER-BELL-LINEAR-JOVERC",
            "species_mass=1",
        ):
            with self.subTest(required=required):
                self.assertIn(required, note)

    def test_design_rejects_the_superseded_Q023_source_local_oracle_identity(self) -> None:
        combined = RECORD_PATH.read_text(encoding="utf-8") + NOTE_PATH.read_text(
            encoding="utf-8"
        )
        self.assertNotIn("Q023-BELL-JOVERC-DEPOSITED-CURRENT-ORACLE", combined)
        self.assertNotIn("q023_bell_joverc_deposited_current_oracle", combined)

    def test_rejects_volume_blind_or_artificial_C_closure(self) -> None:
        for closure in (
            "PPC*deposit_qscale*species_charge*v_CR=J_CR/c=2*B0*k0",
            "PPC*deposit_qscale*species_charge*v_CR/V_root_cell=2*B0*C*k0",
            UNIT_MASS_SPECIALIZED_CLOSURE,
        ):
            with self.subTest(closure=closure):
                candidate = copy.deepcopy(self.record)
                candidate["normalization_contract"]["exact_general_closure"] = closure
                with self.assertRaisesRegex(ContractError, "general closure drifted"):
                    validate_design(candidate)

    def test_rejects_bare_q_over_mc_without_explicit_unit_mass_specialization(self) -> None:
        candidate = copy.deepcopy(self.record)
        candidate["normalization_contract"]["q_over_mc_specialization"][
            "allowed_only_if"
        ] = "always"
        with self.assertRaisesRegex(ContractError, "unit-mass specialization"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["normalization_contract"]["q_over_mc_specialization"][
            "bare_q_over_mc_in_general_closure_allowed"
        ] = True
        with self.assertRaisesRegex(ContractError, "bare q/"):
            validate_design(candidate)

    def test_rejects_authority_or_weakened_predecessor_gate(self) -> None:
        candidate = copy.deepcopy(self.record)
        candidate["authority"]["physical_pilot_authorized"] = True
        with self.assertRaisesRegex(ContractError, "authority"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["mandatory_predecessor_gates"]["corrected_linear_predecessor"][
            "source_local_preparation_alone_is_sufficient"
        ] = True
        with self.assertRaisesRegex(ContractError, "falsely became sufficient"):
            validate_design(candidate)

    def test_rejects_wrong_Q043_or_drifted_predecessor_order(self) -> None:
        candidate = copy.deepcopy(self.record)
        candidate["authoritative_corrected_lineage"]["first_predecessor"][
            "required_snapshot"
        ] = "cycle_zero"
        with self.assertRaisesRegex(ContractError, "raw cycle-one"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["authoritative_corrected_lineage"]["first_predecessor"][
            "required_result"
        ] = "pending"
        with self.assertRaisesRegex(ContractError, "raw cycle-one"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["authoritative_corrected_lineage"]["second_predecessor"][
            "required_after"
        ] = "nothing"
        with self.assertRaisesRegex(ContractError, "corrected Q023"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["source_local_evidence_bindings"][
            "authoritative_Q043_raw_cycle_one_oracle_dependency"
        ]["Q023_source_local_or_cycle_zero_oracle_substitution_allowed"] = True
        with self.assertRaisesRegex(ContractError, "Q043 oracle dependency"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["mandatory_predecessor_gates"]["corrected_linear_predecessor"][
            "required_after"
        ] = "before_Q043"
        with self.assertRaisesRegex(ContractError, "no longer follows"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["authoritative_corrected_lineage"]["Q019_pgen_name"] = (
            "q023_paper_bell_linear"
        )
        with self.assertRaisesRegex(ContractError, "Q019 pgen"):
            validate_design(candidate)

    def test_rejects_centered_beam_PPC_branch_conflation_or_missing_box_control(self) -> None:
        candidate = copy.deepcopy(self.record)
        candidate["branch_separation"]["high_rigidity_fixed_current_like"][
            "required_centered_beam_contract"
        ]["production_PPC"] = 16
        with self.assertRaisesRegex(ContractError, "production PPC"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["branch_separation"]["finite_rigidity_self_consistent"]["campaign_id"] = (
            HIGH_RIGIDITY_ID
        )
        with self.assertRaisesRegex(ContractError, "finite-rigidity ID"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["geometry_and_box_contract"]["3D_required_pair"]["box_control"][
            "geometry"
        ] = "same_small_box"
        with self.assertRaisesRegex(ContractError, "box control"):
            validate_design(candidate)

    def test_rejects_missing_solver_diagnostic_seed_or_resource_gate(self) -> None:
        candidate = copy.deepcopy(self.record)
        candidate["numerical_sensitivity_contract"]["required_rows"].pop()
        with self.assertRaisesRegex(ContractError, "sensitivities"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["required_diagnostics"]["current_and_particle_state"].remove(
            "J_CR_gas=J_CR_lab-rho_CR*u_gas"
        )
        with self.assertRaisesRegex(ContractError, "gas-frame current"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["seed_and_ensemble_contract"][
            "pilot_and_qualifying_seeds_must_be_disjoint"
        ] = False
        with self.assertRaisesRegex(ContractError, "seeds can overlap"):
            validate_design(candidate)

        candidate = copy.deepcopy(self.record)
        candidate["resource_contract"]["historical_Q019_cost_estimates_reusable"] = True
        with self.assertRaisesRegex(ContractError, "resource estimates"):
            validate_design(candidate)


if __name__ == "__main__":
    unittest.main()
