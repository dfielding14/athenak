#!/usr/bin/env python3
"""Focused tests for Q-029 Hall-Bell source-local launch preparation."""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest

from tst.publication import analyze_q029_hall_bell_linear as hall


REPO_ROOT = Path(__file__).resolve().parents[2]
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q029_hall_bell_linear_source_local_preparation_2026-05-30.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q029HallBellLinearTests(unittest.TestCase):
    def test_source_local_decks_freeze_positive_grid_and_registration(self) -> None:
        decks = hall.validate_source_local_candidate_decks()
        self.assertEqual([deck["dimension"] for deck in decks], [1, 2, 3])
        self.assertEqual(
            [deck["launch_status"] for deck in decks],
            [hall.LAUNCH_STATUS] * 3,
        )
        self.assertEqual(
            decks[0]["carrier_semantics"],
            "physical_1d_transverse_invariant_thin_2d3v_nx2_4",
        )
        for deck in decks:
            self.assertEqual(deck["qualification_effect"], "none")
            self.assertFalse(deck["qualifying_evidence"])
            self.assertFalse(deck["hall_bell_qualification"])
            self.assertTrue(deck["pic_enable_2d3v"])
            self.assertEqual(deck["chi_h_grid"], list(hall.CHI_H_VALUES))
            self.assertTrue(all(value > 0.0 for value in deck["chi_h_grid"]))
            geometry = hall._EXPECTED_GEOMETRY[deck["dimension"]]
            self.assertEqual(deck["global_nx"], list(geometry["nx"]))
            self.assertEqual(deck["meshblock_nx"], list(geometry["meshblock_nx"]))
            self.assertEqual(deck["bounds"], [list(item) for item in geometry["bounds"]])

        registration = hall.validate_generator_registration()
        self.assertEqual(registration["pgen_name"], hall.PGEN_NAME)
        self.assertTrue(registration["fresh_dispatch"])
        self.assertTrue(registration["restart_dispatch"])
        self.assertTrue(registration["q023_seed_carrier_reused_from_guarded_header"])

    def test_complete_positive_launch_preparation_bundle_is_nonqualifying(self) -> None:
        report = hall.analyze_launch_preparation_bundle(
            hall.build_launch_preparation_bundle()
        )
        self.assertEqual(report["variant_count"], 45)
        self.assertEqual(report["launch_status"], hall.LAUNCH_STATUS)
        self.assertEqual(report["qualification_effect"], hall.QUALIFICATION_EFFECT)
        self.assertFalse(report["qualifying_evidence"])
        self.assertFalse(report["hall_bell_qualification"])
        self.assertFalse(report["linear_hall_bell_qualification"])
        self.assertFalse(report["nonlinear_hall_bell_qualification"])
        self.assertIn("not_qualifying_evidence", report["status"])
        self.assertEqual(
            {
                (item["dimension"], item["epsilon"], item["chi_h"])
                for item in report["variants"]
            },
            {
                (dimension, epsilon, chi_h)
                for dimension in hall.DIMENSIONS
                for epsilon in hall.EPSILON_VALUES
                for chi_h in hall.CHI_H_VALUES
            },
        )
        variants = {
            (item["dimension"], item["epsilon"], item["chi_h"]): item
            for item in report["variants"]
        }
        low_epsilon = variants[(1, 0.1, 0.25)]
        high_epsilon = variants[(1, 0.8, 0.25)]
        self.assertEqual(low_epsilon["stream_velocity"], [10.0, 0.0, 0.0])
        self.assertEqual(high_epsilon["stream_velocity"], [1.25, 0.0, 0.0])
        self.assertEqual(low_epsilon["pic_cr_light_speed"], 10000.0)
        self.assertEqual(high_epsilon["pic_cr_light_speed"], 1250.0)
        self.assertNotEqual(low_epsilon["j_cr"], high_epsilon["j_cr"])
        self.assertNotEqual(low_epsilon["alpha_h"], high_epsilon["alpha_h"])
        for variant in report["variants"]:
            speed = 1.0 / variant["epsilon"]
            self.assertAlmostEqual(
                math.sqrt(sum(value * value for value in variant["stream_velocity"])),
                speed,
            )
            self.assertAlmostEqual(variant["pic_cr_light_speed"], 1000.0 * speed)
            self.assertAlmostEqual(variant["alpha_h"], variant["chi_h"] / variant["j_cr"])
            self.assertEqual(len(variant["athena_overrides"]), 8)
            self.assertTrue(
                any(item.startswith("job/basename=") for item in variant["athena_overrides"])
            )
            self.assertTrue(
                any(
                    item.startswith("particles/couple_j_to_efield_coeff=")
                    for item in variant["athena_overrides"]
                )
            )

    def test_nonpositive_chi_h_variants_fail_closed(self) -> None:
        for chi_h in (0.0, -0.25):
            with self.subTest(chi_h=chi_h):
                bundle = hall.build_launch_preparation_bundle()
                bundle["variants"][0]["chi_h"] = chi_h
                with self.assertRaisesRegex(hall.ContractError, "prepared positive"):
                    hall.analyze_launch_preparation_bundle(bundle)

    def test_malformed_and_incomplete_launch_preparation_grids_fail_closed(
        self,
    ) -> None:
        bundle = hall.build_launch_preparation_bundle()
        bundle["variants"].pop()
        with self.assertRaisesRegex(hall.ContractError, "grid is incomplete"):
            hall.analyze_launch_preparation_bundle(bundle)

        bundle = hall.build_launch_preparation_bundle()
        bundle["variants"].append(copy.deepcopy(bundle["variants"][0]))
        with self.assertRaisesRegex(hall.ContractError, "duplicate"):
            hall.analyze_launch_preparation_bundle(bundle)

        bundle = hall.build_launch_preparation_bundle()
        bundle["variants"][0]["runtime_artifact_path"] = "/tmp/not-admitted"
        with self.assertRaisesRegex(hall.ContractError, "variant keys"):
            hall.analyze_launch_preparation_bundle(bundle)

        bundle = hall.build_launch_preparation_bundle()
        bundle["runtime_artifact_path"] = "/tmp/not-admitted"
        with self.assertRaisesRegex(hall.ContractError, "bundle keys"):
            hall.analyze_launch_preparation_bundle(bundle)

        with self.assertRaisesRegex(hall.ContractError, "bundle keys"):
            hall.analyze_launch_preparation_bundle([])  # type: ignore[arg-type]

    def test_materialized_variant_drift_and_nonfinite_numbers_fail_closed(self) -> None:
        for key in ("pic_cr_light_speed", "j_cr", "alpha_h"):
            with self.subTest(key=key):
                bundle = hall.build_launch_preparation_bundle()
                bundle["variants"][0][key] *= 2.0
                with self.assertRaisesRegex(hall.ContractError, "materialization"):
                    hall.analyze_launch_preparation_bundle(bundle)

        bundle = hall.build_launch_preparation_bundle()
        bundle["variants"][0]["stream_velocity"][0] *= 2.0
        with self.assertRaisesRegex(hall.ContractError, "stream velocity"):
            hall.analyze_launch_preparation_bundle(bundle)

        bundle = hall.build_launch_preparation_bundle()
        bundle["variants"][0]["athena_overrides"][0] = "job/basename=drifted"
        with self.assertRaisesRegex(hall.ContractError, "Athena overrides"):
            hall.analyze_launch_preparation_bundle(bundle)

        for key in ("epsilon", "chi_h", "pic_cr_light_speed", "j_cr", "alpha_h"):
            with self.subTest(nonfinite_key=key):
                bundle = hall.build_launch_preparation_bundle()
                bundle["variants"][0][key] = math.inf
                with self.assertRaisesRegex(hall.ContractError, "finite number"):
                    hall.analyze_launch_preparation_bundle(bundle)

    def test_candidate_deck_geometry_and_2d3v_drift_fail_closed(self) -> None:
        mutations = (
            (2, "x1min     = 0.0", "x1min     = 0.25", "x1min"),
            (
                2,
                "<meshblock>\nnx1       = 32",
                "<meshblock>\nnx1       = 16",
                "meshblock cell-count",
            ),
            (
                1,
                "pic_enable_2d3v                   = true",
                "pic_enable_2d3v                   = false",
                "pic_enable_2d3v",
            ),
            (2, "nx1       = 64", "nx1       = invalid", "expected an integer"),
            (
                1,
                "cr_vx0                            = 2.5",
                "cr_vx0                            = 0.0",
                "must be positive",
            ),
        )
        for dimension, old, new, error in mutations:
            with self.subTest(dimension=dimension, mutation=error):
                with tempfile.TemporaryDirectory() as temporary_directory:
                    path = Path(temporary_directory) / "candidate.athinput"
                    original = hall.DECKS[dimension].read_text(encoding="utf-8")
                    self.assertIn(old, original)
                    path.write_text(original.replace(old, new, 1), encoding="utf-8")
                    with self.assertRaisesRegex(hall.ContractError, error):
                        hall.validate_candidate_deck(path, dimension)

    def test_runtime_generator_freezes_guarded_header_geometry_and_finite_checks(
        self,
    ) -> None:
        generator = (
            REPO_ROOT / "src/pgen/tests/q029_hall_bell_linear.cpp"
        ).read_text(encoding="utf-8")
        self.assertIn('#include "q023_paper_bell_linear.hpp"', generator)
        self.assertNotIn('#include "q023_paper_bell_linear.cpp"', generator)
        self.assertIn("pmy_mesh_->mesh_indcs.nx1", generator)
        self.assertIn("pmy_mesh_->mb_indcs.nx1", generator)
        self.assertIn(
            'Q029RequireBoolean(pin, "particles", "pic_enable_2d3v", true);',
            generator,
        )
        self.assertGreaterEqual(generator.count("std::isfinite"), 2)

    def test_sidecar_binds_q029_local_artifacts_and_nonqualification_boundary(
        self,
    ) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-029")
        self.assertEqual(sidecar["campaign_id"], hall.CAMPAIGN_ID)
        self.assertEqual(sidecar["qualification_effect"], hall.QUALIFICATION_EFFECT)
        self.assertFalse(sidecar["claim_closure"])
        self.assertEqual(sidecar["frontier_authorization"], "not_bound")
        self.assertEqual(sidecar["launch_status"], hall.LAUNCH_STATUS)

        expected_paths = {
            "src/pgen/tests/q023_paper_bell_linear.hpp",
            "src/pgen/tests/q023_paper_bell_linear.cpp",
            "src/pgen/tests/q029_hall_bell_linear.cpp",
            "inputs/tests/pic_q029_hall_bell_linear_1d_candidate.athinput",
            "inputs/tests/pic_q029_hall_bell_linear_2d_candidate.athinput",
            "inputs/tests/pic_q029_hall_bell_linear_3d_candidate.athinput",
            "tst/publication/analyze_q029_hall_bell_linear.py",
            "tst/publication/test_analyze_q029_hall_bell_linear.py",
        }
        bindings = sidecar["artifact_bindings"]
        self.assertEqual(set(bindings), expected_paths)
        for relative, expected_sha256 in bindings.items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected_sha256)

        self.assertEqual(
            sidecar["shared_registration_contract"],
            {
                "src/CMakeLists.txt":
                    "includes pgen/tests/q029_hall_bell_linear.cpp",
                "src/pgen/pgen.hpp":
                    "declares ProblemGenerator::Q029HallBellLinear",
                "src/pgen/pgen.cpp":
                    "dispatches q029_hall_bell_linear for fresh and restart construction",
            },
        )
        boundary = sidecar["nonqualification_boundary"]
        self.assertEqual(boundary["artifact_role"], hall.ARTIFACT_ROLE)
        self.assertEqual(boundary["qualification_effect"], hall.QUALIFICATION_EFFECT)
        self.assertFalse(boundary["qualifying_evidence"])
        self.assertFalse(boundary["hall_bell_qualification"])
        self.assertTrue(
            any("Frontier" in item for item in sidecar["explicitly_not_claimed"])
        )
        self.assertTrue(
            any("Hall-Bell" in item for item in sidecar["explicitly_not_claimed"])
        )

    def test_sidecar_replays_retained_source_local_runtime_smoke(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        smoke = sidecar["source_local_runtime_smoke"]
        root = Path(smoke["artifact_root"])
        self.assertTrue(root.is_dir())
        files = [path for path in root.rglob("*") if path.is_file()]
        writable = [
            path
            for path in [root, *root.rglob("*")]
            if path.stat().st_mode & 0o222
        ]
        self.assertEqual(len(files), smoke["retention"]["file_count"])
        self.assertEqual(len(writable), smoke["retention"]["writable_entries"])
        self.assertEqual(smoke["retention"]["status"], "pass_recursively_read_only")
        self.assertEqual(_sha256(root / "src/athena"), smoke["executable_sha256"])
        for initialization in smoke["cycle_zero_initializations"]:
            dimension = initialization["dimension"]
            self.assertEqual(
                _sha256(root / f"{dimension}d_cycle_zero.stdout.txt"),
                initialization["stdout_sha256"],
            )
            self.assertEqual(
                _sha256(root / f"{dimension}d_cycle_zero.stderr.txt"),
                initialization["stderr_sha256"],
            )
            self.assertEqual(
                _sha256(root / initialization["selected_raw_mhd_bcc_path"]),
                initialization["selected_raw_mhd_bcc_sha256"],
            )
            self.assertEqual(initialization["result"], "pass")


if __name__ == "__main__":
    unittest.main()
