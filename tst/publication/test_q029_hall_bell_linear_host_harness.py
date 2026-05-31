#!/usr/bin/env python3
"""Compiled host-contract tests for the Q-029 Hall-Bell preparation boundary."""

from __future__ import annotations

import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
HARNESS = REPO_ROOT / "tst/publication/q029_hall_bell_linear_host_harness.cpp"


class Q029HallBellLinearHostHarnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._temporary_directory = tempfile.TemporaryDirectory()
        binary = Path(cls._temporary_directory.name) / "q029-hall-bell-linear-host"
        subprocess.run(
            [
                os.environ.get("CXX", "c++"),
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-Werror",
                str(HARNESS),
                "-o",
                str(binary),
            ],
            cwd=REPO_ROOT,
            check=True,
        )
        cls.lines = subprocess.run(
            [str(binary)],
            cwd=REPO_ROOT,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        ).stdout.splitlines()

    @classmethod
    def tearDownClass(cls) -> None:
        cls._temporary_directory.cleanup()

    def test_current_density_contract(self) -> None:
        lines = [line.split() for line in self.lines if line.startswith("current_density ")]
        self.assertEqual(lines, [["current_density", "4", "20"], ["current_density", "16", "80"]])

    def test_chi_h_contract(self) -> None:
        lines = [line.split() for line in self.lines if line.startswith("chi_h ")]
        self.assertEqual(len(lines), 3)
        for fields in lines:
            alpha_h = float(fields[1])
            self.assertEqual(float(fields[2]), alpha_h*4.0/(2.0*1.0))

    def test_positive_prepared_whitelist_contract(self) -> None:
        lines = [line.split() for line in self.lines if line.startswith("positive_prepared ")]
        self.assertEqual(
            [fields[2] for fields in lines],
            ["0", "0", "1", "1", "0", "1", "1", "0"],
        )

    def test_q023_shared_carrier_reuse(self) -> None:
        lines = [line.split() for line in self.lines if line.startswith("q023_carrier ")]
        self.assertEqual(len(lines), 3)
        for fields in lines:
            dimension = int(fields[1])
            growth = float(fields[2])
            parallel = [float(value) for value in fields[3:6]]
            magnetic = [float(value) for value in fields[6:9]]
            velocity = [float(value) for value in fields[9:12]]
            vector_potential = [float(value) for value in fields[12:15]]
            self.assertAlmostEqual(growth, math.sqrt(1.0 - 0.4*0.4), places=14)
            self.assertAlmostEqual(sum(value*value for value in parallel), 1.0, places=14)
            self.assertAlmostEqual(sum(value*value for value in magnetic), 1.0, places=11)
            self.assertAlmostEqual(
                sum(value*value for value in velocity),
                1.0e-12,
                places=24,
            )
            self.assertTrue(all(math.isfinite(value) for value in vector_potential))
            if dimension == 1:
                self.assertEqual(parallel, [1.0, 0.0, 0.0])
            else:
                self.assertNotEqual(parallel[1], 0.0)


if __name__ == "__main__":
    unittest.main()
