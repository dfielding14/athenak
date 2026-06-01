#!/usr/bin/env python3
"""Focused tests for the Sun & Bai paper_smooth TSC interface oracle."""

from __future__ import annotations

from fractions import Fraction
import json
import os
from pathlib import Path
import subprocess
import sys
import unittest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.scripts.particles import pic_paper_smooth_tsc_oracle as oracle  # noqa: E402

ORACLE_PATH = (
    REPO_ROOT / "tst" / "scripts" / "particles" / "pic_paper_smooth_tsc_oracle.py"
)


class PicPaperSmoothTSCOracleTests(unittest.TestCase):
    def test_frozen_receiver_resolution_grid_preserves_raw_totals(self) -> None:
        expected_totals = {
            "a": Fraction(1),
            "b": Fraction(881, 800),
            "c": Fraction(37, 32),
            "d": Fraction(707, 800),
            "e": Fraction(49, 50),
            "f": Fraction(1),
        }
        self.assertEqual(oracle.EXPECTED_RAW_TOTALS, expected_totals)
        for label, coordinate in oracle.PARTICLES.items():
            with self.subTest(label=label):
                weights = oracle.receiver_raw_weights(coordinate)
                self.assertEqual(weights, oracle.EXPECTED_RAW_WEIGHTS[label])
                self.assertEqual(sum(weights, Fraction(0)), expected_totals[label])

        non_unit = tuple(
            label for label, total in expected_totals.items() if total != 1
        )
        self.assertEqual(non_unit, ("b", "c", "d", "e"))

    def test_receiver_resolution_changes_cross_interface_weights(self) -> None:
        particle_b = oracle.PARTICLES["b"]
        coarse_receiver = oracle.RECEIVERS[4]
        self.assertEqual(coarse_receiver.center, Fraction("0.5"))
        self.assertEqual(coarse_receiver.dx, Fraction(1))
        self.assertEqual(
            oracle.raw_tsc_weight(
                particle_b,
                coarse_receiver.center,
                coarse_receiver.dx,
            ),
            Fraction(81, 800),
        )
        self.assertEqual(
            oracle.raw_tsc_weight(particle_b, coarse_receiver.center, oracle.FINE_DX),
            0,
        )

    def test_support_bound_and_tensor_product(self) -> None:
        self.assertEqual(oracle.raw_tsc_weight(0, Fraction("0.75"), 0.5), 0)
        self.assertEqual(oracle.raw_tsc_weight(0, Fraction("0.76"), 0.5), 0)
        self.assertGreater(oracle.raw_tsc_weight(0, Fraction("0.74"), 0.5), 0)

        particle = (Fraction("-0.55"), Fraction("0.2"), Fraction("-0.4"))
        center = (Fraction("-0.25"), Fraction(0), Fraction(0))
        dx = (Fraction("0.5"), Fraction("0.5"), Fraction("0.5"))
        expected = (
            oracle.raw_tsc_weight(particle[0], center[0], dx[0])
            * oracle.raw_tsc_weight(particle[1], center[1], dx[1])
            * oracle.raw_tsc_weight(particle[2], center[2], dx[2])
        )
        self.assertEqual(expected, Fraction(234171, 4000000))
        self.assertEqual(oracle.tensor_product_weight(particle, center, dx), expected)
        self.assertEqual(
            oracle.tensor_product_weight(
                particle,
                (center[0], center[1], Fraction("0.35")),
                dx,
            ),
            0,
        )

    def test_standalone_executable_emits_checked_report(self) -> None:
        self.assertTrue(os.access(ORACLE_PATH, os.X_OK))
        proc = subprocess.run(
            [str(ORACLE_PATH), "--indent", "0"],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        report = json.loads(proc.stdout)
        self.assertEqual(report["status"], "pass")
        self.assertFalse(report["policy"]["renormalize_cross_interface_totals"])
        self.assertEqual(report["particles"]["c"]["raw_total"]["fraction"], "37/32")
        self.assertEqual(
            report["checks"],
            {
                "frozen_receiver_weight_grid": "pass",
                "raw_non_unit_totals_preserved": "pass",
                "support_bound": "pass",
                "tensor_product": "pass",
            },
        )


if __name__ == "__main__":
    unittest.main()
