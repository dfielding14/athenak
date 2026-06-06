#!/usr/bin/env python3
"""Focused source-local checks for the Q011 resource-scaling decision record."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
PLAN = (
    REPO_ROOT
    / "tst/publication/readiness"
    / "q011_resource_scaling_preproduction_plan_2026-06-06.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q011ResourceScalingPreproductionPlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.plan = json.loads(PLAN.read_text(encoding="utf-8"))

    def test_record_is_planning_only_and_forbids_qualifying_output_inspection(self) -> None:
        self.assertEqual(
            self.plan["record_type"], "q011_resource_scaling_preproduction_plan"
        )
        self.assertEqual(
            self.plan["qualification_effect"],
            "resource_planning_only_no_execution_authorization_no_claim_closure",
        )
        boundaries = self.plan["historical_and_evidence_boundaries"]
        self.assertEqual(boundaries["historical_preregistration_mutation"], "forbidden")
        self.assertEqual(boundaries["qualifying_output_inspection"], "forbidden")
        self.assertIs(boundaries["qualifying_campaign_results_inspected"], False)

    def test_source_local_bindings_are_exact(self) -> None:
        bindings = self.plan["source_and_evidence_bindings"]
        for key in [
            "paper_deck",
            "tracked_qualifying_preregistration_predecessor",
            "storage_planning_successor",
        ]:
            binding = bindings[key]
            self.assertEqual(_sha256(REPO_ROOT / binding["path"]), binding["sha256"])

    def test_budget_snapshot_arithmetic_is_exact(self) -> None:
        audit = self.plan["frontier_budget_audit_snapshot"]
        policy = audit["active_policy"]
        ledger = audit["ledger"]
        self.assertEqual(policy["maximum_node_hours"], 10000.0)
        self.assertEqual(audit["frontier_node_gpu_visible_gcd_count"], 8)
        self.assertEqual(ledger["currently_reserved_node_hours"], 0.0)
        self.assertAlmostEqual(
            ledger["remaining_unreserved_node_hours"],
            policy["maximum_node_hours"]
            - ledger["cumulative_consumed_node_hours"]
            - ledger["currently_reserved_node_hours"],
        )

    def test_geometry_and_work_ratios_recompute(self) -> None:
        geometry = self.plan["geometry_and_cadence"]
        pilot = geometry["pressure_pilot"]
        production = geometry["production"]
        model = self.plan["analytical_work_model"]
        exact = model["exact_inputs"]
        brackets = model["constant_rate_ideal_eight_gpu_node_brackets"]

        i60 = 0.5 * 45**2 + 0.5 * (60 - 45) ** 2
        i1200 = 0.5 * 45**2 + 0.5 * (1200 - 45) ** 2
        self.assertEqual(i60, exact["I_60"])
        self.assertEqual(i1200, exact["I_1200"])

        transverse_ratio = production["x2_extent"] / pilot["x2_extent"]
        coarse_particle = i1200 / i60 * transverse_ratio
        fine_particle = coarse_particle * 4
        coarse_zone = (
            production["coarse_uniform_dx12"]["cell_count"]
            / pilot["active_cells"]
            * production["coarse_uniform_dx12"][
                "projected_cycles_from_pressure_pilot_cadence"
            ]
            / pilot["cycles"]
        )
        fine_zone = (
            production["fine_uniform_dx3"]["cell_count"]
            / pilot["active_cells"]
            * production["fine_uniform_dx3"][
                "projected_cycles_from_pressure_pilot_cadence"
            ]
            / pilot["cycles"]
        )
        self.assertAlmostEqual(
            coarse_particle, brackets["coarse_uniform_dx12"]["particle_work_ratio"]
        )
        self.assertAlmostEqual(
            fine_particle, brackets["fine_uniform_dx3"]["particle_work_ratio"]
        )
        self.assertAlmostEqual(
            coarse_zone, brackets["coarse_uniform_dx12"]["zone_work_ratio"]
        )
        self.assertAlmostEqual(fine_zone, brackets["fine_uniform_dx3"]["zone_work_ratio"])

    def test_driver_brackets_recompute_from_selected_pressure_pilot(self) -> None:
        observations = self.plan["pressure_pilot_observations"]
        brackets = self.plan["analytical_work_model"][
            "constant_rate_ideal_eight_gpu_node_brackets"
        ]
        factor = observations["driver_seconds_rank_max"] / (8 * 3600)
        self.assertAlmostEqual(
            brackets["coarse_uniform_dx12"]["driver_node_hours_lower"],
            factor * brackets["coarse_uniform_dx12"]["particle_work_ratio"],
        )
        self.assertAlmostEqual(
            brackets["coarse_uniform_dx12"]["driver_node_hours_upper"],
            factor * brackets["coarse_uniform_dx12"]["zone_work_ratio"],
        )
        self.assertAlmostEqual(
            brackets["fine_uniform_dx3"]["driver_node_hours_lower"],
            factor * brackets["fine_uniform_dx3"]["particle_work_ratio"],
        )
        self.assertAlmostEqual(
            brackets["fine_uniform_dx3"]["driver_node_hours_upper"],
            factor * brackets["fine_uniform_dx3"]["zone_work_ratio"],
        )

    def test_scenario_envelope_recomputes_from_stated_assumptions(self) -> None:
        scenario = self.plan["scenario_only_likely_cost_envelope"]
        assumptions = scenario["assumptions"]
        brackets = self.plan["analytical_work_model"][
            "constant_rate_ideal_eight_gpu_node_brackets"
        ]

        def triad_cost(improvement: float, efficiency: float, amr_fraction: float) -> float:
            coarse = brackets["coarse_uniform_dx12"]
            amr = brackets["three_level_amr_root_dx12_finest_dx3"]
            fine = brackets["fine_uniform_dx3"]
            return sum(
                [
                    max(
                        coarse["driver_node_hours_lower"] / efficiency,
                        coarse["driver_node_hours_upper"]
                        / (improvement * efficiency),
                    ),
                    max(
                        amr["driver_node_hours_particle_floor"] / efficiency,
                        amr["driver_node_hours_full_fine_zone_value"]
                        * amr_fraction
                        / (improvement * efficiency),
                    ),
                    max(
                        fine["driver_node_hours_lower"] / efficiency,
                        fine["driver_node_hours_upper"] / (improvement * efficiency),
                    ),
                ]
            )

        improvement = assumptions[
            "large_grid_per_gpu_throughput_improvement_over_tiny_pressure_pilot"
        ]
        efficiency = assumptions["multi_node_efficiency"]
        amr_fraction = assumptions[
            "amr_time_weighted_active_cell_fraction_of_full_fine_domain"
        ]
        uplift = assumptions["output_wrapper_restart_and_failure_uplift_factor"]
        lower_triad = triad_cost(max(improvement), max(efficiency), min(amr_fraction))
        upper_triad = triad_cost(min(improvement), min(efficiency), max(amr_fraction))

        self.assertEqual(
            scenario["three_complete_paired_triads_total_node_hours"],
            [lower_triad * 3 * min(uplift), upper_triad * 3 * max(uplift)],
        )
        self.assertEqual(
            scenario["current_eight_complete_paired_triads_total_node_hours"],
            [lower_triad * 8 * min(uplift), upper_triad * 8 * max(uplift)],
        )

    def test_recommended_matrix_is_three_complete_paired_triads(self) -> None:
        current = self.plan["current_campaign_matrix"]
        recommended = self.plan["recommended_matrix"]
        self.assertEqual(len(current["grid_variants"]), 3)
        self.assertEqual(len(current["qualifying_seeds"]), 8)
        self.assertEqual(current["baseline_attempt_count"], 24)
        self.assertEqual(recommended["core_grid_variants"], current["grid_variants"])
        self.assertEqual(
            recommended["core_paired_qualifying_seeds"],
            current["qualifying_seeds"][:3],
        )
        self.assertEqual(recommended["core_baseline_attempt_count"], 9)
        self.assertEqual(
            recommended["expansion_order"],
            current["qualifying_seeds"][3:],
        )
        self.assertIs(recommended["matrix_successor_required"], True)

    def test_engineering_seeds_are_disjoint_and_pilot_cap_is_bounded(self) -> None:
        qualifying = set(self.plan["current_campaign_matrix"]["qualifying_seeds"])
        design = self.plan["preproduction_scaling_pilot_design"]
        engineering = set(design["engineering_seed_policy"]["seeds"])
        self.assertTrue(qualifying.isdisjoint(engineering))
        self.assertLessEqual(design["hard_consumed_node_hour_cap"], 500.0)
        self.assertIn(
            "Do not inspect or reduce engineering-pilot physical fields for scientific conclusions.",
            design["common_requirements"],
        )

    def test_project_budget_gate_retains_non_q011_reserve_and_failure_allowance(self) -> None:
        gates = {gate["gate_id"]: gate for gate in self.plan["decision_gates"]}
        condition = gates["Q011-RS-4"]["pass_condition"]
        definitions = self.plan["analytical_work_model"][
            "project_authorization_budget_term_definitions"
        ]
        self.assertIn("max(C_fine, 0.15 * C_baseline)", condition)
        self.assertIn("R_nonQ011", condition)
        self.assertIn("10000", condition)
        self.assertEqual(
            set(definitions),
            {
                "C_consumed",
                "C_reserved",
                "C_preproduction",
                "C_baseline",
                "C_restart",
                "C_fine",
                "R_nonQ011",
            },
        )
        self.assertIn("failed-attempt replacement reserve", definitions["C_fine"])
        self.assertIn("non-Q011", definitions["R_nonQ011"])
        self.assertIn(
            "protected non-Q011",
            self.plan["frontier_budget_audit_snapshot"]["budget_scope_warning"],
        )


if __name__ == "__main__":
    unittest.main()
