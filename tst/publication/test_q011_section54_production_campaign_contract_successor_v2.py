#!/usr/bin/env python3
"""Adversarial tests for the fail-closed Q-011 production-campaign successor."""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest

from tst.publication import q011_resource_scaling_preproduction_pilot_materializer_successor_v2 as pilots
from tst.publication import q011_section54_production_campaign_contract_successor_v2 as contract
from tst.publication import q011_section54_production_campaign_planner_successor_v2 as planner
from tst.publication import (
    q011_section54_production_science_admission_orchestration_successor_v2 as admission,
)
from tst.publication import q011_section54_production_science_analyzer_successor_v2 as analyzer
from tst.publication import q011_section54_qualifying_campaign_execution_successor_v2 as execution
from tst.publication import q011_section54_registered_launch_materializer_successor_v2 as launch


def valid_resource_freeze_inputs() -> dict[str, object]:
    """Return resource-only inputs for a largest-affordable N=5 freeze."""
    return {
        "global_ledger_snapshot_binding": {
            "path": "registered/ledger/node_hours_snapshot.json",
            "sha256": "a" * 64,
        },
        "all_other_registered_campaign_reservations_manifest_binding": {
            "path": "registered/campaigns/all_other_reservations.json",
            "sha256": "b" * 64,
        },
        "exact_deck_scaling_io_pilot_completion_binding": {
            "path": "registered/q011/exact_deck_pilot_completion.json",
            "sha256": "c" * 64,
        },
        "global_ledger_cumulative_consumed_node_hours_after_excluded_pilots": 6000.0,
        "global_ledger_currently_reserved_node_hours_after_all_other_campaigns": 1000.0,
        "q011_unreserved_nonbaseline_obligations_node_hours": 500.0,
        "measured_conservative_complete_paired_triad_node_hours": 400.0,
        "exact_deck_scaling_io_pilots_passed": True,
        "exact_deck_scaling_io_pilots_excluded_from_qualifying_science": True,
        "all_other_registered_campaigns_reserved": True,
        "qualifying_production_reservations_present": False,
        "qualifying_launch_started": False,
        "qualifying_output_inspected": False,
    }


class Q011ProductionCampaignContractSuccessorV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.pipeline = contract.build_pipeline_contract()
        cls.plan = planner.build_plan()
        cls.pilot_manifest, cls.pilot_files = pilots.build_materialization()
        cls.launch_candidates = launch.build_review_candidates(
            cls.plan, cls.pilot_manifest
        )
        cls.execution_handoffs = execution.build_execution_handoffs(
            cls.plan, cls.launch_candidates
        )
        cls.admission = admission.build_source_local_admission(cls.execution_handoffs)
        cls.analysis = analyzer.build_source_local_analysis_packet(cls.admission)
        cls.resource_freeze = contract.build_resource_only_freeze_record(
            valid_resource_freeze_inputs()
        )
        cls.expanded_plan = planner.build_plan(cls.resource_freeze)
        cls.expanded_pilot_manifest, cls.expanded_pilot_files = pilots.build_materialization(
            cls.expanded_plan
        )
        cls.expanded_launch_candidates = launch.build_review_candidates(
            cls.expanded_plan, cls.expanded_pilot_manifest
        )
        cls.expanded_execution_handoffs = execution.build_execution_handoffs(
            cls.expanded_plan, cls.expanded_launch_candidates
        )
        cls.expanded_admission = admission.build_source_local_admission(
            cls.expanded_execution_handoffs
        )
        cls.expanded_analysis = analyzer.build_source_local_analysis_packet(
            cls.expanded_admission
        )

    def test_pipeline_binds_one_exact_nine_output_deck_end_to_end(self) -> None:
        deck = self.pipeline["production_deck"]
        self.assertEqual(deck["binding"]["path"], contract.SUCCESSOR_DECK)
        self.assertEqual(deck["binding"]["sha256"], contract.SUCCESSOR_DECK_SHA256)
        self.assertEqual(deck["launch_path"], contract.SUCCESSOR_DECK_LAUNCH_PATH)
        self.assertEqual(deck["output_count"], 9)
        self.assertEqual(deck["outputs"], contract.EXPECTED_OUTPUTS)
        self.assertEqual(self.pipeline["stage_order"], list(contract.STAGE_ORDER))
        for stage in self.pipeline["stages"]:
            self.assertEqual(stage["production_deck"], deck["binding"])
            self.assertEqual(
                stage["production_deck_launch_path"],
                contract.SUCCESSOR_DECK_LAUNCH_PATH,
            )
            self.assertEqual(stage["authorization"], contract.AUTHORIZATION_BOUNDARY)

    def test_stage_record_types_form_one_closed_handoff_chain(self) -> None:
        stages = self.pipeline["stages"]
        self.assertEqual(stages[0]["input_record_type"], "none_source_local_entrypoint")
        for previous, current in zip(stages, stages[1:]):
            self.assertEqual(previous["emitted_record_type"], current["input_record_type"])
            self.assertEqual(current["predecessor_stage"], previous["stage_id"])

    def test_campaign_size_design_preregisters_core_and_reserve_without_promising_24_runs(
        self,
    ) -> None:
        design = self.pipeline["campaign_size_resource_design"]
        self.assertEqual(
            design["mandatory_core_qualifying_seeds"],
            list(contract.CORE_QUALIFYING_SEEDS),
        )
        self.assertEqual(
            design["preregistered_reserve_qualifying_seeds"],
            list(contract.RESERVE_QUALIFYING_SEEDS),
        )
        self.assertEqual(
            design["ordered_qualifying_seed_pool"],
            list(contract.QUALIFYING_SEED_POOL),
        )
        numeric = {rule["rule_id"]: rule for rule in design["numeric_rules"]}
        self.assertEqual(numeric["mandatory_core_triad_count"]["value"], 3)
        self.assertEqual(numeric["mandatory_core_attempt_count"]["value"], 9)
        self.assertEqual(numeric["maximum_preregistered_triad_count"]["value"], 8)
        self.assertEqual(numeric["maximum_preregistered_attempt_count"]["value"], 24)
        self.assertIn(
            "not a promised run count",
            numeric["maximum_preregistered_attempt_count"]["rationale"],
        )
        self.assertEqual(self.plan["attempt_count"], 9)
        self.assertIsNone(
            self.plan["campaign_size_selection"]["resource_only_freeze_record_sha256"]
        )

    def test_resource_only_freeze_selects_largest_affordable_complete_prefix_end_to_end(
        self,
    ) -> None:
        freeze = contract.validate_resource_only_freeze_record(self.resource_freeze)
        self.assertEqual(freeze["frozen_triad_count"], 5)
        self.assertEqual(freeze["frozen_attempt_count"], 15)
        self.assertEqual(
            freeze["selected_qualifying_seeds"],
            list(contract.QUALIFYING_SEED_POOL[:5]),
        )
        evaluations = freeze["candidate_budget_evaluations"]
        self.assertEqual(
            [item["triad_count"] for item in evaluations if item["fits_global_ledger"]],
            [3, 4, 5],
        )
        self.assertEqual(evaluations[2]["projected_global_node_hours"], 9900.0)
        self.assertEqual(evaluations[3]["projected_global_node_hours"], 10300.0)
        maximum_inputs = valid_resource_freeze_inputs()
        maximum_inputs[
            "global_ledger_cumulative_consumed_node_hours_after_excluded_pilots"
        ] = 0.0
        maximum_inputs[
            "global_ledger_currently_reserved_node_hours_after_all_other_campaigns"
        ] = 0.0
        maximum_inputs["q011_unreserved_nonbaseline_obligations_node_hours"] = 0.0
        self.assertEqual(
            contract.build_resource_only_freeze_record(maximum_inputs)[
                "frozen_triad_count"
            ],
            8,
        )
        self.assertEqual(self.expanded_plan["attempt_count"], 15)
        self.assertNotIn(
            "resource_only_campaign_size_freeze_not_recorded_before_qualifying_launch",
            self.expanded_plan["blockers"],
        )
        selected = set(freeze["selected_qualifying_seeds"])
        for variant in contract.GRID_VARIANTS:
            self.assertEqual(
                {
                    attempt["qualifying_seed"]
                    for attempt in self.expanded_plan["attempts"]
                    if attempt["variant"] == variant
                },
                selected,
            )
        for value in (
            self.expanded_pilot_manifest,
            self.expanded_launch_candidates,
            self.expanded_execution_handoffs,
            self.expanded_admission,
            self.expanded_analysis,
        ):
            self.assertEqual(
                value["campaign_size_selection"],
                self.expanded_plan["campaign_size_selection"],
            )
            self.assertEqual(value["authorization"], contract.AUTHORIZATION_BOUNDARY)

    def test_resource_only_freeze_rejects_science_inputs_bad_timing_and_unaffordable_core(
        self,
    ) -> None:
        for key, value in (
            ("qualifying_output_inspected", True),
            ("qualifying_launch_started", True),
            ("qualifying_production_reservations_present", True),
            ("exact_deck_scaling_io_pilots_passed", False),
            ("exact_deck_scaling_io_pilots_excluded_from_qualifying_science", False),
            ("all_other_registered_campaigns_reserved", False),
        ):
            mutated = valid_resource_freeze_inputs()
            mutated[key] = value
            with self.assertRaises(contract.ProductionCampaignContractError):
                contract.build_resource_only_freeze_record(mutated)

        science_driven = valid_resource_freeze_inputs()
        science_driven["observed_acceleration_efficiency"] = 0.2
        with self.assertRaisesRegex(
            contract.ProductionCampaignContractError, "prohibited fields"
        ):
            contract.build_resource_only_freeze_record(science_driven)

        unaffordable = valid_resource_freeze_inputs()
        unaffordable[
            "global_ledger_cumulative_consumed_node_hours_after_excluded_pilots"
        ] = 9000.0
        with self.assertRaisesRegex(
            contract.ProductionCampaignContractError, "mandatory three-triad core"
        ):
            contract.build_resource_only_freeze_record(unaffordable)

    def test_frozen_prefix_seed_dropping_and_post_execution_mutation_fail_closed(
        self,
    ) -> None:
        mutated_freeze = copy.deepcopy(self.resource_freeze)
        mutated_freeze["selected_qualifying_seeds"].pop()
        with self.assertRaises(contract.ProductionCampaignContractError):
            contract.validate_resource_only_freeze_record(mutated_freeze)
        non_largest_freeze = copy.deepcopy(self.resource_freeze)
        non_largest_freeze["frozen_triad_count"] = 4
        with self.assertRaises(contract.ProductionCampaignContractError):
            contract.validate_resource_only_freeze_record(non_largest_freeze)

        mutated_plan = copy.deepcopy(self.expanded_plan)
        mutated_plan["attempts"] = [
            attempt
            for attempt in mutated_plan["attempts"]
            if attempt["qualifying_seed"] != contract.QUALIFYING_SEED_POOL[3]
        ]
        mutated_plan["attempt_count"] = len(mutated_plan["attempts"])
        with self.assertRaises(planner.ProductionCampaignPlannerError):
            planner.validate_plan(mutated_plan)

        mutated_execution = copy.deepcopy(self.expanded_execution_handoffs)
        mutated_execution["campaign_size_selection"]["planned_triad_count"] = 4
        with self.assertRaises(execution.CampaignExecutionSuccessorError):
            execution.validate_execution_handoffs(
                mutated_execution,
                self.expanded_plan,
                self.expanded_launch_candidates,
            )

    def test_reporting_rule_is_honest_for_n3_and_exactly_frozen_for_larger_n(
        self,
    ) -> None:
        n3 = contract.reporting_rule_for_frozen_n(3)
        self.assertFalse(n3["interval_estimator_frozen"])
        self.assertFalse(n3["broad_population_inference_authorized"])
        self.assertIn("descriptive", n3["rationale"])
        expected_coverage = {4: 0.875, 5: 0.9375, 6: 0.96875, 7: 0.875, 8: 0.9296875}
        for triad_count, coverage in expected_coverage.items():
            rule = contract.reporting_rule_for_frozen_n(triad_count)
            self.assertTrue(rule["interval_estimator_frozen"])
            self.assertEqual(rule["nominal_coverage"], 0.8)
            self.assertAlmostEqual(rule["achieved_finite_sample_coverage"], coverage)
            self.assertFalse(rule["broad_population_inference_authorized"])
        self.assertEqual(
            self.expanded_analysis["frozen_reporting_rule"],
            contract.reporting_rule_for_frozen_n(5),
        )

    def test_statistical_and_resource_numeric_rules_have_explicit_provenance(
        self,
    ) -> None:
        validated = contract.validate_statistical_resource_design()
        for section in (
            validated["campaign_size_resource_design"],
            validated["reporting_uncertainty_policy"],
        ):
            for rule in section["numeric_rules"]:
                self.assertIn(
                    rule["source_category"], contract.THRESHOLD_SOURCE_CATEGORIES
                )
                self.assertGreaterEqual(len(rule["rationale"].strip()), 24)
                self.assertNotEqual(rule["source_category"], "literature comparison")

    def test_runtime_source_closure_covers_requested_runtime_paths(self) -> None:
        closure = self.pipeline["runtime_source_closure"]
        members = {member["path"] for member in closure["members"]}
        required = {
            "src/particles/particles_pushers.cpp",
            "src/particles/particles_moments.cpp",
            "src/particles/particles_tasks.cpp",
            "src/bvals/bvals_part.cpp",
            "src/eos/eos.cpp",
            "src/eos/eos.hpp",
            "src/mhd/mhd_tasks.cpp",
            "src/mhd/mhd_update.cpp",
            "src/srcterms/srcterms.cpp",
            "src/srcterms/srcterms.hpp",
            "src/srcterms/srcterms_newdt.cpp",
            "src/mesh/mesh_refinement.cpp",
            "src/mesh/load_balance.cpp",
            "src/pgen/tests/pic_parallel_shock.cpp",
            "src/outputs/derived_variables.cpp",
            contract.SUCCESSOR_DECK,
        }
        self.assertTrue(required.issubset(members))
        self.assertEqual(closure["member_count"], len(closure["members"]))
        self.assertEqual(closure["sha256"], contract.canonical_sha256(closure["members"]))
        self.assertGreaterEqual(len(closure["groups"]), 5)

    def test_negative_gradient_preregistration_supersedes_future_rule_only(self) -> None:
        preregistration = contract.preregistration_binding()
        record = preregistration["record"]
        self.assertEqual(
            record["shock_front"]["detector"],
            "unique_strongest_negative_density_gradient",
        )
        self.assertEqual(
            record["shock_front"]["predecessor_positive_gradient_rule"],
            "superseded_for_future_successor_campaign_only",
        )
        self.assertFalse(record["shock_front"]["historical_records_modified"])
        self.assertEqual(
            self.analysis["shock_front_detector"],
            "unique_strongest_negative_density_gradient",
        )

    def test_fail_closed_science_and_engineering_gate_inventory_is_complete(self) -> None:
        self.assertEqual(
            [record["gate_id"] for record in self.analysis["gates"]],
            list(contract.GATE_DEFINITIONS),
        )
        self.assertEqual(self.analysis["gate_count"], 7)
        self.assertEqual(self.analysis["passed_gate_count"], 0)
        self.assertFalse(self.analysis["measurements_inspected"])
        for gate in self.analysis["gates"]:
            self.assertEqual(
                gate["status"], "blocked_pending_registered_runtime_evidence"
            )
            self.assertFalse(gate["passed"])
            self.assertFalse(gate["measurements_inspected"])
            self.assertFalse(gate["qualifying_evidence"])
            self.assertEqual(gate["authorization"], contract.AUTHORIZATION_BOUNDARY)

    def test_every_science_threshold_has_explicit_nonmisleading_provenance(self) -> None:
        definitions = contract.validate_gate_definition_threshold_provenance()
        self.assertEqual(
            self.pipeline["threshold_provenance_policy"]["allowed_source_categories"],
            sorted(contract.THRESHOLD_SOURCE_CATEGORIES),
        )
        for definition in definitions.values():
            self.assertNotIn("pass_rule", definition)
            for criterion in definition["acceptance_criteria"]:
                self.assertIn(
                    criterion["source_category"], contract.THRESHOLD_SOURCE_CATEGORIES
                )
                self.assertGreaterEqual(len(criterion["rationale"].strip()), 24)
                if criterion["source_category"] == "literature comparison":
                    rationale = criterion["rationale"].lower()
                    self.assertIn("comparison", rationale)
                    self.assertTrue(
                        "not" in rationale
                        or "no numeric" in rationale
                        or "no numerical" in rationale
                    )
            selection = definition["snapshot_selection"]
            self.assertIn(
                selection["source_category"], contract.THRESHOLD_SOURCE_CATEGORIES
            )
            self.assertGreaterEqual(len(selection["rationale"].strip()), 24)
        mutated = copy.deepcopy(contract.source_local_gate_records())
        del mutated[0]["definition"]["acceptance_criteria"][0]["source_category"]
        with self.assertRaises(contract.ProductionCampaignContractError):
            contract.validate_source_local_gate_records(mutated)

    def test_exact_deck_compute_and_io_pilots_are_materialized_fail_closed(self) -> None:
        self.assertEqual(self.pilot_manifest["production_deck"]["path"], contract.SUCCESSOR_DECK)
        self.assertEqual(self.pilot_manifest["successor_output_count"], 9)
        self.assertEqual(self.pilot_manifest["deck_count"], 6)
        self.assertEqual(self.pilot_manifest["case_count"], 6)
        self.assertEqual(
            self.pilot_manifest["gate_status"],
            "blocked_pending_registered_runtime_evidence",
        )
        budget = self.pilot_manifest["maximum_consumed_node_hours"]
        self.assertEqual(budget["value"], 500.0)
        self.assertEqual(budget["source_category"], "engineering closure")
        self.assertIn("no scientific literature origin", budget["rationale"])
        self.assertEqual(
            self.pilot_manifest["resource_only_freeze_definition"],
            contract.resource_only_freeze_definition(),
        )
        self.assertIn(
            "resource-only campaign sizing",
            self.pilot_manifest[
                "resource_only_freeze_must_use_excluded_pilot_measurements"
            ]["rationale"],
        )
        for deck_record in self.pilot_manifest["decks"]:
            blocks = contract.parse_deck(
                self.pilot_files[deck_record["path"]].decode("utf-8")
            )
            outputs = {
                name: values for name, values in blocks.items() if name.startswith("output")
            }
            self.assertEqual(set(outputs), set(contract.EXPECTED_OUTPUTS))
            if deck_record["phase"] == "compute_scaling":
                self.assertTrue(
                    all(float(outputs[f"output{index}"]["dt"]) <= 0.0 for index in range(1, 10))
                )
            else:
                self.assertEqual(outputs, contract.EXPECTED_OUTPUTS)

    def test_complete_source_local_chain_is_consistent_and_non_authorizing(self) -> None:
        self.assertEqual(self.plan["attempt_count"], 9)
        self.assertEqual(self.launch_candidates["candidate_count"], 9)
        self.assertEqual(self.execution_handoffs["handoff_count"], 9)
        self.assertEqual(self.admission["attempt_count"], 9)
        self.assertEqual(self.admission["admitted_attempt_count"], 0)
        self.assertEqual(self.analysis["passed_gate_count"], 0)
        for value in (
            self.pipeline,
            self.plan,
            self.pilot_manifest,
            self.launch_candidates,
            self.execution_handoffs,
            self.admission,
            self.analysis,
        ):
            self.assertEqual(value["authorization"], contract.AUTHORIZATION_BOUNDARY)
        for attempt in self.plan["attempts"]:
            self.assertEqual(attempt["argv"][1], contract.SUCCESSOR_DECK_LAUNCH_PATH)
            self.assertEqual(attempt["production_deck"]["path"], contract.SUCCESSOR_DECK)

    def test_protected_historical_artifacts_remain_byte_exact(self) -> None:
        observed = {
            record["path"]: record["sha256"]
            for record in contract.protected_historical_bindings()
        }
        self.assertEqual(observed, contract.PROTECTED_HISTORICAL_ARTIFACTS)

    def test_readiness_record_binds_exact_sources_and_deterministic_contracts(self) -> None:
        path = (
            contract.READINESS_ROOT
            / "q011_section54_production_campaign_contract_successor_v2_2026-06-06.json"
        )
        readiness = json.loads(path.read_text(encoding="utf-8"))
        for binding in readiness["source_bindings"].values():
            payload = (contract.REPO_ROOT / binding["path"]).read_bytes()
            self.assertEqual(hashlib.sha256(payload).hexdigest(), binding["sha256"])
        deterministic = readiness["deterministic_contract_bindings"]
        self.assertEqual(
            deterministic["pipeline_contract_sha256"],
            contract.canonical_sha256(contract.build_pipeline_contract()),
        )
        self.assertEqual(
            deterministic["runtime_source_closure_sha256"],
            contract.runtime_source_closure()["sha256"],
        )
        self.assertEqual(
            deterministic["source_local_gate_records_sha256"],
            contract.canonical_sha256(contract.source_local_gate_records()),
        )
        self.assertEqual(readiness["authorization"], contract.AUTHORIZATION_BOUNDARY)

    def test_mutated_or_six_output_deck_is_rejected(self) -> None:
        payload = (contract.REPO_ROOT / contract.SUCCESSOR_DECK).read_bytes()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory).resolve() / "mutated.athinput"
            path.write_bytes(payload.replace(b"<output9>", b"<output10>", 1))
            with self.assertRaisesRegex(
                contract.ProductionCampaignContractError, "SHA-256"
            ):
                contract.validate_successor_deck(path)
            path.write_bytes(
                (
                    contract.REPO_ROOT
                    / "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
                ).read_bytes()
            )
            with self.assertRaisesRegex(
                contract.ProductionCampaignContractError, "SHA-256"
            ):
                contract.validate_successor_deck(path)

    def test_authority_and_gate_result_injection_are_rejected(self) -> None:
        mutated_pipeline = copy.deepcopy(self.pipeline)
        mutated_pipeline["authorization"]["launch_authorized"] = True
        with self.assertRaisesRegex(
            contract.ProductionCampaignContractError, "acquired authority"
        ):
            contract.validate_pipeline_contract(mutated_pipeline)

        mutated_analysis = copy.deepcopy(self.analysis)
        mutated_analysis["gates"][0]["passed"] = True
        mutated_analysis["passed_gate_count"] = 1
        with self.assertRaisesRegex(analyzer.AnalyzerContractSuccessorError, "authority/results"):
            analyzer.validate_source_local_analysis_packet(
                mutated_analysis, self.admission
            )

    def test_source_closure_gradient_and_stage_lineage_tampering_are_rejected(self) -> None:
        for omitted in (
            "src/srcterms/srcterms.cpp",
            "src/srcterms/srcterms.hpp",
            "src/srcterms/srcterms_newdt.cpp",
            "src/eos/eos.cpp",
            "src/eos/eos.hpp",
        ):
            mutated = copy.deepcopy(self.pipeline)
            mutated["runtime_source_closure"]["members"] = [
                member
                for member in mutated["runtime_source_closure"]["members"]
                if member["path"] != omitted
            ]
            with self.assertRaises(contract.ProductionCampaignContractError):
                contract.validate_pipeline_contract(mutated)

        mutated = copy.deepcopy(self.pipeline)
        mutated["stages"][4]["production_deck"]["path"] = (
            "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
        )
        with self.assertRaises(contract.ProductionCampaignContractError):
            contract.validate_pipeline_contract(mutated)

        mutated = copy.deepcopy(self.analysis)
        mutated["shock_front_detector"] = "unique_strongest_positive_density_gradient"
        with self.assertRaises(analyzer.AnalyzerContractSuccessorError):
            analyzer.validate_source_local_analysis_packet(mutated, self.admission)

    def test_planner_launch_execution_admission_tampering_is_rejected(self) -> None:
        mutated_plan = copy.deepcopy(self.plan)
        mutated_plan["attempts"][0]["argv"][1] = (
            "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
        )
        with self.assertRaises(planner.ProductionCampaignPlannerError):
            planner.validate_plan(mutated_plan)

        mutated_launch = copy.deepcopy(self.launch_candidates)
        mutated_launch["candidates"][0]["registered_policy_slice_status"] = "authorized"
        with self.assertRaises(launch.RegisteredLaunchSuccessorError):
            launch.validate_review_candidates(
                mutated_launch, self.plan, self.pilot_manifest
            )

        mutated_execution = copy.deepcopy(self.execution_handoffs)
        mutated_execution["handoffs"][0]["registered_execution_receipt"] = {
            "forged": True
        }
        with self.assertRaises(execution.CampaignExecutionSuccessorError):
            execution.validate_execution_handoffs(
                mutated_execution, self.plan, self.launch_candidates
            )

        mutated_admission = copy.deepcopy(self.admission)
        mutated_admission["attempts"][0]["admitted"] = True
        mutated_admission["admitted_attempt_count"] = 1
        with self.assertRaises(admission.SourceLocalAdmissionSuccessorError):
            admission.validate_source_local_admission(
                mutated_admission, self.execution_handoffs
            )

    def test_resource_pilot_manifest_or_deck_byte_tampering_is_rejected(self) -> None:
        mutated_manifest = copy.deepcopy(self.pilot_manifest)
        mutated_manifest["authorization"]["scheduler_submission_authorized"] = True
        with self.assertRaises(pilots.ResourcePilotMaterializationError):
            pilots.validate_materialization(mutated_manifest, self.pilot_files)

        mutated_files = dict(self.pilot_files)
        path = next(iter(mutated_files))
        mutated_files[path] = mutated_files[path] + b"\n"
        with self.assertRaises(pilots.ResourcePilotMaterializationError):
            pilots.validate_materialization(self.pilot_manifest, mutated_files)

    def test_nonfinite_json_and_unknown_stage_fail_closed(self) -> None:
        with self.assertRaises(contract.ProductionCampaignContractError):
            contract.canonical_json_bytes({"bad": math.nan})
        with self.assertRaisesRegex(contract.ProductionCampaignContractError, "unknown"):
            contract.build_stage_contract("not_a_stage")


if __name__ == "__main__":
    unittest.main()
