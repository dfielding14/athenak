#!/usr/bin/env python3
"""Focused tests for the Q-011 Section 5.4 normalization oracle."""

from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tst"))

from scripts.particles import pic_parallel_shock_section54_paper as q011  # noqa: E402


class PicQ011Section54NormalizationOracleTests(unittest.TestCase):
    def test_frozen_deck_derives_normalization_and_downstream_ppc(self) -> None:
        contract = q011.build_preparation_contract()
        self.assertEqual(contract["schema_version"], 2)
        calibration = contract["deck"][
            "normalization_and_macro_particle_calibration"
        ]

        self.assertAlmostEqual(calibration["alfven_speed_u_a0"], 1.0)
        self.assertAlmostEqual(calibration["normalized_q_over_mc"], 1.0)
        self.assertAlmostEqual(calibration["cyclotron_frequency_omega0"], 1.0)
        self.assertAlmostEqual(
            calibration["ion_inertial_length_c_over_omega_pi"],
            1.0,
        )
        self.assertAlmostEqual(
            calibration["ion_inertial_length_via_u_a0_over_omega0"],
            1.0,
        )
        self.assertAlmostEqual(calibration["numerical_light_speed_c"], 10000.0)
        self.assertAlmostEqual(
            calibration["numerical_light_speed_over_u_a0"],
            10000.0,
        )

        ambiguity = calibration["convention_ambiguity"]
        self.assertEqual(
            ambiguity["status"],
            "resolved_for_the_frozen_athenak_deck",
        )
        self.assertIn("normalized q/(m*c)", ambiguity["athenak_deck_mapping"])
        self.assertAlmostEqual(
            ambiguity["rejected_double_division_omega0"],
            1.0e-4,
        )

        cells = calibration["cell_sizes"]
        self.assertEqual([item["level"] for item in cells], [0, 1, 2])
        self.assertEqual(
            [item["dx_c_over_omega_pi"] for item in cells],
            [12.0, 6.0, 3.0],
        )
        self.assertEqual(
            [item["dy_c_over_omega_pi"] for item in cells],
            [12.0, 6.0, 3.0],
        )
        self.assertEqual(
            [item["effective_2d_cell_volume"] for item in cells],
            [144.0, 36.0, 9.0],
        )

        self.assertAlmostEqual(calibration["ideal_shock_speed_over_u_a0"], 10.0)
        self.assertAlmostEqual(
            calibration["upstream_relative_swept_speed_over_u_a0"],
            40.0,
        )
        self.assertAlmostEqual(calibration["macro_particle"]["mass"], 9.0e-4)

        downstream = calibration["ideal_downstream_calibration"]
        self.assertEqual(
            downstream["effective_dimension"],
            "2d_with_collapsed_x3_thickness",
        )
        self.assertAlmostEqual(downstream["compression_ratio"], 4.0)
        self.assertAlmostEqual(downstream["macro_particle_density"], 40.0 / 9.0)
        self.assertEqual(downstream["target_ppc_by_level"], [640.0, 160.0, 40.0])
        for measured, expected in zip(
            downstream["expected_ppc_by_level"],
            (640.0, 160.0, 40.0),
        ):
            self.assertAlmostEqual(measured, expected)
        self.assertAlmostEqual(downstream["expected_coarse_ppc"], 640.0)
        self.assertAlmostEqual(downstream["expected_fine_ppc"], 40.0)
        self.assertAlmostEqual(downstream["macro_mass_from_coarse_ppc"], 9.0e-4)
        self.assertAlmostEqual(downstream["macro_mass_from_fine_ppc"], 9.0e-4)

        self.assertNotIn(
            "downstream_40_640_ppc_macro_mass_calibration",
            contract["open_items"],
        )
        self.assertIn(
            "gas_pressure_thermodynamic_normalization_audit",
            contract["open_items"],
        )

    def test_qscale_drift_fails_closed(self) -> None:
        blocks = copy.deepcopy(q011.parse_athinput())
        blocks["particles"]["deposit_qscale"] = "1.0e-3"
        with self.assertRaisesRegex(q011.ContractError, "qscale to macro mass"):
            q011.derive_normalization_and_macro_particle_calibration(blocks)

    def test_artificial_light_speed_drift_fails_closed(self) -> None:
        blocks = copy.deepcopy(q011.parse_athinput())
        blocks["particles"]["pic_cr_light_speed"] = "1.0"
        with self.assertRaisesRegex(q011.ContractError, "C/U_A0"):
            q011.derive_normalization_and_macro_particle_calibration(blocks)


if __name__ == "__main__":
    unittest.main()
