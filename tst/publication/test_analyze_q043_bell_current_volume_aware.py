#!/usr/bin/env python3
"""Focused tests for the corrected Q-043 volume-aware Bell current preparation."""

from __future__ import annotations

import copy
import inspect
import json
import math
from pathlib import Path
import tempfile
import unittest

from tst.publication import analyze_q043_bell_current_volume_aware as bell
from tst.publication import (
    q043_bell_current_volume_aware_deposited_current_oracle as oracle,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q043_bell_current_volume_aware_successor_source_local_preparation_2026-06-06.json"
)


def _replace_once(path: Path, old: str, new: str) -> Path:
    text = path.read_text(encoding="utf-8")
    if text.count(old) != 1:
        raise AssertionError(f"expected exactly one occurrence of {old!r}")
    temporary = Path(tempfile.mkdtemp()) / path.name
    temporary.write_text(text.replace(old, new, 1), encoding="utf-8")
    return temporary


class Q043BellCurrentVolumeAwareTests(unittest.TestCase):
    def test_corrected_decks_close_deposited_j_over_c_and_remain_non_authorizing(
        self,
    ) -> None:
        decks = bell.validate_source_local_candidate_decks()
        self.assertEqual([deck["dimension"] for deck in decks], [1, 2, 3])
        expected = 2.0 * math.pi * 2.0
        qscale_per_root_volume = []
        for deck in decks:
            self.assertEqual(deck["campaign_id"], bell.CAMPAIGN_ID)
            self.assertEqual(deck["pgen_name"], bell.PGEN_NAME)
            self.assertEqual(
                deck["supersedes_campaign_id"], bell.SUPERSEDES_CAMPAIGN_ID
            )
            self.assertEqual(deck["current_normalization"], bell.CURRENT_NORMALIZATION)
            self.assertEqual(deck["supersession_effect"], bell.SUPERSESSION_EFFECT)
            self.assertEqual(deck["ppc"], 1.0)
            self.assertGreater(deck["root_cell_volume"], 0.0)
            self.assertAlmostEqual(
                deck["deposit_qscale"], deck["required_deposit_qscale"]
            )
            qscale_per_root_volume.append(
                deck["deposit_qscale"] / deck["root_cell_volume"]
            )
            self.assertAlmostEqual(deck["deposited_j_over_c"], expected)
            self.assertLessEqual(deck["deposited_j_over_c_absolute_residual"], 1.0e-13)
            self.assertAlmostEqual(
                deck["expected_deposited_j_over_c"],
                expected,
            )
            self.assertAlmostEqual(
                deck["historical_to_corrected_current_ratio"],
                deck["artificial_light_speed"] / deck["root_cell_volume"],
            )
            self.assertAlmostEqual(
                deck["unit_root_volume_qscale_to_corrected_current_ratio"],
                1.0 / deck["root_cell_volume"],
            )
            self.assertEqual(
                deck["deposited_current_output_variables"],
                ["prtcl_jx", "prtcl_jy", "prtcl_jz"],
            )
            self.assertTrue(
                deck["pgen_source_implementation_bound_by_this_preparation"]
            )
            self.assertTrue(deck["pgen_registration_bound_by_this_preparation"])
            self.assertEqual(deck["species_mass"], 1.0)
            self.assertTrue(
                deck["integration_binding"][
                    "integrated_raw_current_oracle_identity_compatible"
                ]
            )
            self.assertFalse(
                deck["integration_binding"][
                    "integrated_raw_current_oracle_accepted_as_qualification_evidence"
                ]
            )
            self.assertTrue(
                deck["integration_binding"]["integrated_raw_current_oracle_required"]
            )
            self.assertFalse(
                deck["integration_binding"][
                    "duplicate_q043_generator_or_oracle_required"
                ]
            )
            self.assertFalse(deck["launch_authorized"])
            self.assertFalse(deck["qualification_eligible"])
            self.assertEqual(deck["qualification_effect"], bell.QUALIFICATION_EFFECT)
        self.assertEqual(len({deck["deposit_qscale"] for deck in decks}), 3)
        for value in qscale_per_root_volume[1:]:
            self.assertAlmostEqual(value, qscale_per_root_volume[0])

    def test_current_formula_excludes_artificial_c_and_multic_observations_pass(
        self,
    ) -> None:
        signature = inspect.signature(bell.deposited_j_over_c_vector)
        self.assertNotIn("artificial_light_speed", signature.parameters)
        self.assertIn("root_cell_volume", signature.parameters)
        bundle = bell.synthetic_observed_current_bundle((250.0, 2500.0, 25000.0))
        report = bell.analyze_observed_current_bundle(bundle)
        self.assertTrue(report["normalization_contract_pass"])
        self.assertTrue(report["observed_deposited_current_contract_pass"])
        self.assertTrue(report["artificial_c_invariance_pass"])
        self.assertFalse(report["launch_authorized"])
        self.assertFalse(report["qualification_eligible"])
        self.assertFalse(report["scientific_claim_authorized"])
        self.assertFalse(report["passed"])
        self.assertEqual(report["qualification_effect"], bell.QUALIFICATION_EFFECT)
        for result in report["invariance_by_dimension"]:
            self.assertEqual(
                result["artificial_light_speeds"],
                [250.0, 2500.0, 25000.0],
            )
            self.assertTrue(result["artificial_c_invariance_pass"])
            self.assertAlmostEqual(result["maximum_absolute_spread"], 0.0)

    def test_current_formula_scales_with_root_cell_volume_and_not_decomposition(
        self,
    ) -> None:
        velocity = (2.5, 0.0, 0.0)
        charge = 2.0 * math.pi * 1.0e-6
        target = 4.0 * math.pi
        for root_cell_volume in (1.0, 0.125, 0.03125, 1.0e-4):
            qscale = bell.required_deposit_qscale(
                1.0, charge, 2.5, root_cell_volume, 1.0, 2.0 * math.pi
            )
            current = bell.deposited_j_over_c_vector(
                1.0, qscale, charge, velocity, root_cell_volume
            )
            self.assertAlmostEqual(current[0], target)
            self.assertEqual(current[1:], (0.0, 0.0))
            unit_volume_qscale = bell.required_deposit_qscale(
                1.0, charge, 2.5, 1.0, 1.0, 2.0 * math.pi
            )
            fixed_current = bell.deposited_j_over_c_vector(
                1.0, unit_volume_qscale, charge, velocity, root_cell_volume
            )
            self.assertAlmostEqual(fixed_current[0] / target, 1.0 / root_cell_volume)

    def test_historical_c_multiplied_current_observation_fails_closed(self) -> None:
        bundle = bell.synthetic_observed_current_bundle((250.0, 2500.0))
        bundle["records"][0]["volume_mean_deposited_j_over_c"][0] *= 2500.0
        with self.assertRaisesRegex(
            bell.ContractError, "observed-current dimension 1 x1"
        ):
            bell.analyze_observed_current_bundle(bundle)

    def test_observed_current_bundle_requires_complete_multic_dimension_grid(
        self,
    ) -> None:
        bundle = bell.synthetic_observed_current_bundle((250.0, 2500.0))
        bundle["records"] = [
            record for record in bundle["records"] if record["dimension"] != 3
        ]
        with self.assertRaisesRegex(
            bell.ContractError, "dimension 3 requires at least two"
        ):
            bell.analyze_observed_current_bundle(bundle)

        bundle = bell.synthetic_observed_current_bundle((250.0, 2500.0))
        bundle["records"].append(copy.deepcopy(bundle["records"][0]))
        with self.assertRaisesRegex(bell.ContractError, "duplicate observed-current"):
            bell.analyze_observed_current_bundle(bundle)

        bundle = bell.synthetic_observed_current_bundle((250.0, 2500.0))
        bundle["records"][0]["unexpected"] = "not allowed"
        with self.assertRaisesRegex(bell.ContractError, "record keys"):
            bell.analyze_observed_current_bundle(bundle)

    def test_corrected_deck_drift_and_old_deck_fail_closed(self) -> None:
        mutations = (
            (
                1,
                "integrator = vl2",
                "integrator = rk2",
                "time/integrator",
            ),
            (
                1,
                "pic_background_mode               = coupled",
                (
                    "pic_physical_mode                 = paper_mhd_pic_vl2_tsc\n"
                    "pic_background_mode               = coupled"
                ),
                "particles/pic_physical_mode",
            ),
            (
                1,
                "deposit_qscale                    = 6250.0",
                "deposit_qscale                    = 1.0e9",
                "deposit_qscale",
            ),
            (
                2,
                "ppc                               = 1.0",
                "ppc                               = 2.0",
                "particles/ppc",
            ),
            (
                3,
                "current_normalization = deposited_j_over_c_equals_2_b_g_k0",
                "current_normalization = deposited_j_over_c_equals_2_b_g_C_k0",
                "current_normalization",
            ),
            (
                1,
                "variable    = prtcl_jx",
                "variable    = prtcl_d",
                "output3/variable",
            ),
            (
                2,
                "pgen_name = q043_bell_current_volume_aware",
                "pgen_name = q023_paper_bell_linear",
                "problem/pgen_name",
            ),
            (
                1,
                "mass   = 1.0",
                "mass   = 2.0",
                "species0/mass",
            ),
        )
        for dimension, old, new, error in mutations:
            with self.subTest(dimension=dimension, mutation=error):
                path = _replace_once(bell.DECKS[dimension], old, new)
                with self.assertRaisesRegex(bell.ContractError, error):
                    bell.validate_candidate_deck(path, dimension)

        historical = (
            REPO_ROOT
            / "inputs/tests/pic_q023_paper_bell_linear_1d_candidate_vl2_tsc.athinput"
        )
        with self.assertRaisesRegex(bell.ContractError, "particles/ppc"):
            bell.validate_candidate_deck(historical, 1)

    def test_preparation_report_explicitly_invalidates_old_qualification(self) -> None:
        report = bell.build_preparation_report()
        self.assertTrue(report["normalization_contract_pass"])
        self.assertFalse(report["artificial_c_in_formula"])
        self.assertTrue(report["root_cell_volume_in_formula"])
        self.assertFalse(report["fixed_qscale_valid"])
        self.assertTrue(report["dimension_and_resolution_aware_qscale_required"])
        self.assertTrue(report["artificial_c_invariance_runtime_observation_pending"])
        self.assertEqual(report["supersedes_campaign_id"], bell.SUPERSEDES_CAMPAIGN_ID)
        disposition = report["historical_campaign_disposition"]
        self.assertEqual(disposition["campaign_id"], bell.SUPERSEDES_CAMPAIGN_ID)
        self.assertEqual(
            disposition["normalization"], bell.HISTORICAL_INVALID_NORMALIZATION
        )
        self.assertFalse(disposition["qualification_eligible"])
        self.assertIn("must_not_be_used_for_qualification", disposition["disposition"])
        self.assertFalse(report["launch_authorized"])
        self.assertFalse(report["qualification_eligible"])
        self.assertFalse(report["passed"])
        binding = report["integration_binding"]
        self.assertEqual(
            binding["selected_successor_campaign_id"],
            bell.CAMPAIGN_ID,
        )
        self.assertEqual(
            binding["selected_generator_name"],
            bell.PGEN_NAME,
        )
        self.assertEqual(
            binding["selected_readiness_record_prefix"],
            bell.READINESS_RECORD_PREFIX,
        )
        self.assertEqual(
            binding["integrated_raw_current_oracle_campaign_id"],
            oracle.CAMPAIGN_ID,
        )
        self.assertEqual(
            binding["integrated_raw_current_oracle_binding_status"],
            "source_local_ready_runtime_observation_pending",
        )
        self.assertTrue(binding["integrated_raw_current_oracle_identity_compatible"])
        self.assertTrue(binding["integrated_raw_current_oracle_required"])
        self.assertFalse(
            binding["integrated_raw_current_oracle_accepted_as_qualification_evidence"]
        )
        self.assertEqual(
            binding["supersession_identity_role"],
            "authoritative_foundational_current_lineage",
        )

    def test_selected_lineage_matches_integrated_raw_current_oracle(self) -> None:
        self.assertEqual(oracle.PGEN_CAMPAIGN_ID, bell.CAMPAIGN_ID)
        self.assertEqual(oracle.PGEN_NAME, bell.PGEN_NAME)
        self.assertEqual(oracle.CAMPAIGN_ID, bell.INTEGRATED_RAW_CURRENT_ORACLE_CAMPAIGN_ID)
        self.assertEqual(bell.PGEN_NAME, "q043_bell_current_volume_aware")
        self.assertEqual(oracle.CAMPAIGN_ID, bell.INTEGRATED_RAW_CURRENT_ORACLE_CAMPAIGN_ID)
        for deck_path in bell.DECKS.values():
            blocks = bell.parse_athinput(deck_path)
            self.assertEqual(blocks["problem"]["pgen_name"], oracle.PGEN_NAME)
            self.assertEqual(blocks[bell.PGEN_BLOCK]["campaign_id"], oracle.PGEN_CAMPAIGN_ID)

    def test_readiness_record_preserves_successor_preparation_metadata(self) -> None:
        readiness = json.loads(READINESS.read_text(encoding="utf-8"))
        self.assertEqual(readiness["gate"], "Q-043")
        self.assertEqual(readiness["campaign_id"], bell.CAMPAIGN_ID)
        self.assertEqual(
            readiness["supersedes_campaign_id"], bell.SUPERSEDES_CAMPAIGN_ID
        )
        self.assertEqual(
            readiness["current_normalization"], bell.CURRENT_NORMALIZATION
        )
        self.assertEqual(
            readiness["supersession_effect"], bell.SUPERSESSION_EFFECT
        )
        self.assertEqual(readiness["qualification_effect"], bell.QUALIFICATION_EFFECT)
        self.assertFalse(readiness["authority"]["launch_authorized"])
        self.assertFalse(readiness["authority"]["qualification_eligible"])
        self.assertFalse(readiness["authority"]["scientific_claim_authorized"])
        self.assertTrue(
            readiness["source_implementation_boundary"][
                "pgen_source_implementation_bound_by_this_preparation"
            ]
        )
        self.assertTrue(
            readiness["source_implementation_boundary"][
                "pgen_registration_bound_by_this_preparation"
            ]
        )
        self.assertEqual(
            readiness["historical_campaign_disposition"]["qualification_eligible"],
            False,
        )
        binding = readiness["integration_binding"]
        self.assertEqual(binding["selected_generator_name"], bell.PGEN_NAME)
        self.assertEqual(
            binding["selected_readiness_record_prefix"],
            bell.READINESS_RECORD_PREFIX,
        )
        self.assertEqual(
            binding["integrated_raw_current_oracle_campaign_id"],
            oracle.CAMPAIGN_ID,
        )
        self.assertEqual(
            binding["integrated_raw_current_oracle_identity"],
            bell.INTEGRATED_RAW_CURRENT_ORACLE_IDENTITY,
        )
        self.assertEqual(
            binding["integrated_raw_current_oracle_binding_status"],
            bell.INTEGRATED_RAW_CURRENT_ORACLE_BINDING_STATUS,
        )
        self.assertTrue(binding["integrated_raw_current_oracle_identity_compatible"])
        self.assertTrue(binding["integrated_raw_current_oracle_required"])
        self.assertFalse(
            binding["integrated_raw_current_oracle_accepted_as_qualification_evidence"]
        )
        self.assertEqual(
            binding["supersession_identity_role"],
            bell.SUPERSESSION_IDENTITY_ROLE,
        )

        # This dated record is historical metadata. Live deck byte integrity is
        # enforced by the generator-backed Q043 matrix validator.


if __name__ == "__main__":
    unittest.main()
