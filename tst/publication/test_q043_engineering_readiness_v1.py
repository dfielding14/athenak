#!/usr/bin/env python3
"""Tests for the compact Q043 engineering-readiness workflow."""

from __future__ import annotations

import hashlib
from pathlib import Path
import tempfile
import unittest

from tst.publication import q043_engineering_readiness_v1 as readiness


class Q043EngineeringReadinessTests(unittest.TestCase):
    def test_selection_has_expected_orthogonal_coverage(self) -> None:
        cases = readiness.selected_cases()
        self.assertEqual(len(cases), 26)
        self.assertEqual(len({case["case_id"] for case in cases}), 26)

        baseline = {
            (case["dimension"], case["resolution"], case["ppc"])
            for case in cases
            if case["decomposition"] == "single"
            and case["artificial_c_over_v_cr"] == 1000
        }
        self.assertEqual(
            baseline,
            {
                (dimension, resolution, ppc)
                for dimension in (1, 2, 3)
                for resolution in ("coarse", "fine")
                for ppc in (1, 4)
            },
        )
        self.assertEqual(
            {
                (case["dimension"], case["artificial_c_over_v_cr"])
                for case in cases
                if case["resolution"] == "coarse"
                and case["ppc"] == 1
                and case["decomposition"] == "single"
                and case["artificial_c_over_v_cr"] != 1000
            },
            {
                (dimension, ratio)
                for dimension in (1, 2, 3)
                for ratio in (100, 10000)
            },
        )
        self.assertEqual(
            {
                (case["dimension"], case["decomposition"])
                for case in cases
                if case["decomposition"] != "single"
            },
            {
                (1, "split_x1"),
                (2, "split_x1"),
                (2, "split_x2"),
                (2, "split_x1x2"),
                (3, "split_x1"),
                (3, "split_x2"),
                (3, "split_x3"),
                (3, "split_xyz"),
            },
        )

    def test_manifest_binds_canonical_deck_bytes(self) -> None:
        manifest = readiness.build_manifest()
        self.assertEqual(manifest["engineering_case_count"], 26)
        self.assertEqual(manifest["canonical_q043_case_count"], 132)
        self.assertFalse(manifest["replaces_complete_q043_matrix"])
        self.assertFalse(manifest["scientific_claim_authorized"])
        for record in manifest["cases"]:
            path = readiness.oracle.REPO_ROOT / record["deck_path"]
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(),
                record["deck_sha256"],
            )

    def test_launch_command_uses_case_rank_count_and_raw_root(self) -> None:
        case = next(
            case
            for case in readiness.selected_cases()
            if case["case_id"]
            == "q043-current-oracle-d3-fine-ppc4-split_xyz-cvr1000"
        )
        command = readiness._command(
            case,
            executable=Path("/tmp/athena"),
            execution_deck=Path("/tmp/execution.athinput"),
            raw_root=Path("/tmp/raw"),
            launcher="/usr/bin/srun",
        )
        self.assertIn("--ntasks=8", command)
        self.assertIn("--gpu-bind=closest", command)
        self.assertIn("/tmp/execution.athinput", command)
        self.assertEqual(command[-3:], ["-d", "/tmp/raw", "time/nlim=1"])

    def test_execution_deck_adds_exact_species_velocities(self) -> None:
        case = next(
            case
            for case in readiness.selected_cases()
            if case["case_id"]
            == "q043-current-oracle-d2-coarse-ppc1-single-cvr100"
        )
        rendered = readiness._execution_deck_bytes(case).decode("utf-8")
        blocks = readiness.oracle.parse_athinput_text(rendered)
        for axis in ("x", "y", "z"):
            self.assertEqual(
                blocks["species0"][f"v{axis}0"],
                blocks["particles"][f"cr_v{axis}0"],
            )

    def test_dry_run_writes_all_commands_without_raw_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = root / "athena"
            executable.write_bytes(b"placeholder")
            output = root / "run"
            result = readiness.run_campaign(
                executable=executable,
                output_root=output,
                launcher="/usr/bin/srun",
                dry_run=True,
            )
            self.assertIsNone(result)
            self.assertTrue(
                (output / "engineering_readiness_manifest.json").is_file()
            )
            self.assertEqual(len(list(output.glob("*/command.json"))), 26)
            self.assertFalse(
                (output / "engineering_readiness_summary.json").exists()
            )


if __name__ == "__main__":
    unittest.main()
