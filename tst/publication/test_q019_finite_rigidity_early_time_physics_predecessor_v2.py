#!/usr/bin/env python3
"""Focused tests for the Q019 finite-rigidity physics predecessor v2."""

from __future__ import annotations

import copy
import unittest

from tst.publication import analyze_q019_physics_first_nonlinear_bell_successor_v2 as analysis
from tst.publication import q019_finite_rigidity_early_time_physics_predecessor_v2 as predecessor
from tst.publication.test_q019_physics_first_nonlinear_bell_successor_v2 import _synthetic_evidence


class Q019FiniteRigidityPhysicsPredecessorV2Tests(unittest.TestCase):
    def test_reference_grid_and_convergence_contract_are_explicit(self) -> None:
        contract = predecessor.build_contract()
        self.assertEqual(len(contract["reference_case_ids"]), 9)
        self.assertEqual(len(contract["control_case_ids"]), 6)
        self.assertFalse(contract["reference_definition"]["analytic_finite_shell_dispersion_relation_claimed"])
        self.assertIsNone(contract["fit_window"])
        self.assertIsNone(contract["numeric_acceptance_tolerances"])
        self.assertFalse(contract["runtime_predecessor_complete"])
        self.assertFalse(any(contract["authority"].values()))

    def test_analysis_coverage_is_required_but_cannot_self_pass(self) -> None:
        case_id = "q019-fr-predecessor-k8-rho3em06-s0"
        snapshots, particles = _synthetic_evidence(case_id)
        report = analysis.analyze_snapshots(
            case_id, snapshots, particles, source_kind="synthetic_contract_fixture"
        )
        self.assertTrue(predecessor.validate_analysis_report(report)["analysis_contract_observed"])
        drifted = copy.deepcopy(report)
        drifted["finite_rigidity_early_time_physics_predecessor"].pop(
            "signed_complex_deposited_Jperp_k0_trace"
        )
        with self.assertRaisesRegex(predecessor.PredecessorContractError, "inventory"):
            predecessor.validate_analysis_report(drifted)


if __name__ == "__main__":
    unittest.main()
