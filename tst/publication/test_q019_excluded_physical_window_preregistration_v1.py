#!/usr/bin/env python3
"""Tests for the non-authorizing Q019 physical-window preregistration."""

from __future__ import annotations

import copy
import math
import unittest

from tst.publication import q019_excluded_physical_window_preregistration_v1 as gate
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as design


class Q019ExcludedPhysicalWindowPreregistrationTests(unittest.TestCase):
    def test_inventory_is_exact_non_authorizing_and_seed_disjoint(self) -> None:
        value = gate.build_preregistration()
        self.assertEqual(len(value["case_bindings"]), 19)
        self.assertEqual(
            [len(stage["case_ids"]) for stage in value["stages"]], [7, 8, 4]
        )
        self.assertFalse(value["saturation_evidence_eligible"])
        self.assertTrue(not any(value["authorization"].values()))
        self.assertEqual(
            value["resource_boundary"]["effective_node_hour_ceiling"],
            "min(500,0.10*then_unreserved_project_node_hours)",
        )
        self.assertEqual(
            value["resource_boundary"]["maximum_attempts_per_case"], 1
        )
        self.assertEqual(value["resource_boundary"]["maximum_retries_per_case"], 0)
        pilot_field = {row["field_seed"] for row in value["case_bindings"]}
        pilot_particle = {row["particle_seed"] for row in value["case_bindings"]}
        qualifying = value["qualifying_seed_inventory"]
        self.assertTrue(pilot_field.isdisjoint(qualifying["field_seeds"]))
        self.assertTrue(pilot_particle.isdisjoint(qualifying["particle_seeds"]))
        known = {row["case_id"] for row in design.expected_cases()}
        self.assertEqual(
            {row["case_id"] for row in value["case_bindings"]} - known, set()
        )

    def test_preregistration_rejects_post_inspection_rule_change(self) -> None:
        value = gate.build_preregistration()
        drifted = copy.deepcopy(value)
        drifted["numeric_gates"]["maximum_paired_relative_difference"] = 0.20
        with self.assertRaisesRegex(gate.PreregistrationError, "drifted"):
            gate.validate_preregistration(drifted)

    def test_plateau_selection_uses_earliest_sustained_window(self) -> None:
        k0 = 2.0 * math.pi
        times = [0.1 * index for index in range(80)]
        tau = [k0 * value for value in times]
        amplitude = []
        for value in tau:
            if value < 12.0:
                amplitude.append(1.0e-3 * math.exp(0.65 * value))
            else:
                amplitude.append(2.5 * math.exp(0.015 * math.sin(value)))
        selected = gate.select_plateau_window(
            times, amplitude, k0=k0, u_a=1.0, field_output_dt=0.1
        )
        self.assertGreaterEqual(selected["onset_tau"], 10.0)
        self.assertLessEqual(
            abs(selected["plateau_log_energy_slope_per_tau"]),
            gate.MAXIMUM_ABSOLUTE_LOG_ENERGY_SLOPE_PER_TAU,
        )
        self.assertGreater(
            selected["recommended_production_terminal_tau"],
            selected["plateau_end_tau"],
        )
        self.assertFalse(selected["production_freeze_authorized"])

    def test_plateau_selection_fails_without_onset_or_plateau(self) -> None:
        times = [0.1 * index for index in range(80)]
        with self.assertRaisesRegex(gate.PreregistrationError, "onset"):
            gate.select_plateau_window(
                times,
                [0.1] * len(times),
                k0=2.0 * math.pi,
                u_a=1.0,
                field_output_dt=0.1,
            )
        with self.assertRaisesRegex(gate.PreregistrationError, "plateau"):
            gate.select_plateau_window(
                times,
                [math.exp(0.2 * index) for index in range(len(times))],
                k0=2.0 * math.pi,
                u_a=1.0,
                field_output_dt=0.1,
            )


if __name__ == "__main__":
    unittest.main()
