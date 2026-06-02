#!/usr/bin/env python3
"""Focused tests for the deterministic Q-011 storage estimator sidecar."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

if __package__:
    from . import q011_parallel_shock_storage_estimator as storage
else:
    import q011_parallel_shock_storage_estimator as storage


READINESS = (
    storage.REPO_ROOT
    / "tst/publication/readiness"
    / "q011_parallel_shock_storage_estimator_successor_v3_2026-06-02.json"
)


def _approx_quantity(value):
    if isinstance(value, int):
        return float(value)
    return value["approx"]


def _sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q011ParallelShockStorageEstimatorTests(unittest.TestCase):
    def _temporary_deck(self, transform) -> Path:
        source = storage.DEFAULT_DECK.read_text(encoding="utf-8")
        temporary = tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".athinput",
            encoding="utf-8",
            delete=False,
        )
        self.addCleanup(Path(temporary.name).unlink, missing_ok=True)
        with temporary:
            temporary.write(transform(source))
        return Path(temporary.name)

    def test_current_deck_projects_section54_particle_and_raw_pvtk_budget(self) -> None:
        report = storage.estimate_storage()
        self.assertEqual(storage.pvtk_bytes_per_particle(), 56)
        self.assertEqual(report["pvtk_layout"]["bytes_per_particle"], 56)
        self.assertEqual(report["per_run"]["pvtk_snapshot_count"], 13)
        self.assertEqual(
            [item["time"] for item in report["per_run"]["snapshots"]],
            list(range(0, 1201, 100)),
        )
        self.assertAlmostEqual(
            _approx_quantity(
                report["per_run"]["retained_particle_count_at_tlim_analytical"]
            ),
            1.60e8,
            delta=0.01e8,
        )
        self.assertAlmostEqual(
            report["per_run"]["raw_pvtk_payload_gb_decimal_before_overhead"],
            56.4,
            delta=0.1,
        )

    def test_campaign_multiplier_is_parameterized_by_grids_and_seeds(self) -> None:
        report = storage.estimate_storage(grid_count=2, seed_count=5)
        self.assertEqual(
            report["campaign_parameters"],
            {"grid_count": 2, "seed_count": 5, "run_count": 10},
        )
        self.assertAlmostEqual(
            _approx_quantity(
                report["campaign"]["raw_pvtk_payload_bytes_before_overhead"]
            ),
            10
            * _approx_quantity(
                report["per_run"]["raw_pvtk_payload_bytes_before_overhead"]
            ),
        )
        with self.assertRaisesRegex(storage.EstimatorError, "positive integer"):
            storage.estimate_storage(grid_count=0)

    def test_assumptions_are_explicit(self) -> None:
        report = storage.estimate_storage()
        assumptions = report["assumptions"]
        self.assertIn("no particle escape", assumptions["particle_retention"])
        self.assertIn("t=0", assumptions["cadence"])
        self.assertIn("timestep-edge", assumptions["quantization"])
        self.assertIn("raw-PVTK", assumptions["payload_scope"])
        self.assertIn("restart planning allowances", assumptions["payload_scope"])
        self.assertIn(
            "ASCII headers",
            report["planning_envelope"]["mesh_bin_policy"]["ascii_headers"],
        )
        self.assertIn(
            "not observed runtime bytes",
            report["planning_envelope"]["evidence_boundary"],
        )

    def test_mesh_bin_products_are_projected_for_each_required_variant(self) -> None:
        variants = storage.estimate_storage()["planning_envelope"]["variants"]
        self.assertEqual(
            [variant["variant"] for variant in variants],
            [
                "coarse_uniform_dx12",
                "three_level_amr_root_dx12_finest_dx3",
                "fine_uniform_dx3",
            ],
        )
        self.assertEqual(
            [variant["cell_count"] for variant in variants],
            [1_040_000, 16_640_000, 16_640_000],
        )
        self.assertEqual(
            [variant["meshblock_count"] for variant in variants],
            [2_600, 41_600, 41_600],
        )
        self.assertEqual(
            [
                variant["mesh_bin_products"][
                    "binary_payload_bytes_before_ascii_headers"
                ]
                for variant in variants
            ],
            [228_217_600, 3_651_481_600, 3_651_481_600],
        )
        for variant in variants:
            products = variant["mesh_bin_products"]["products"]
            self.assertEqual(
                [(product["id"], product["snapshot_count"]) for product in products],
                [("bmag", 13), ("j2", 13), ("prtcl_jx", 13), ("rho", 13)],
            )

    def test_restart_payload_is_an_explicit_planning_allowance(self) -> None:
        envelope = storage.estimate_storage()["planning_envelope"]
        restart_policy = envelope["restart_policy"]
        self.assertEqual(restart_policy["checkpoint_times"], list(range(100, 1201, 100)))
        self.assertEqual(
            restart_policy["particle_layout"],
            {
                "real_field_count": 26,
                "integer_field_count": 4,
                "real_bytes_planning_allowance": 8,
                "integer_bytes_planning_allowance": 4,
                "bytes_per_particle_planning_allowance": 224,
                "classification": (
                    "source-layout-derived allowance, not measured checkpoint bytes"
                ),
            },
        )
        self.assertIn("not measured", restart_policy["classification"])
        for variant in envelope["variants"]:
            self.assertEqual(variant["restart_payload"]["checkpoint_count"], 12)

    def test_full_envelope_applies_allowances_without_claiming_replication(self) -> None:
        envelope = storage.estimate_storage()["planning_envelope"]
        campaign = envelope["campaign"]
        logical = _approx_quantity(campaign["logical_bytes_before_filesystem_overhead"])
        filesystem = _approx_quantity(
            campaign["filesystem_allocation_overhead_bytes_allowance"]
        )
        replicated = _approx_quantity(campaign["bytes_after_replication_policy"])
        margin = _approx_quantity(campaign["safety_margin_bytes_allowance"])
        reservation = _approx_quantity(campaign["reservation_envelope_bytes"])
        self.assertAlmostEqual(filesystem, logical * 0.10)
        self.assertEqual(envelope["replication_policy"]["copy_count"], 1)
        self.assertIn(
            "does not provide",
            envelope["replication_policy"]["durability_risk"],
        )
        self.assertAlmostEqual(replicated, logical + filesystem, delta=0.01)
        self.assertAlmostEqual(margin, replicated * 0.25, delta=0.01)
        self.assertAlmostEqual(reservation, replicated + margin, delta=0.01)
        self.assertAlmostEqual(
            campaign["reservation_envelope_tb_decimal"],
            10.553535598859627,
        )

    def test_noncanonical_grid_count_uses_conservative_finest_equivalents(self) -> None:
        variants = storage.estimate_storage(grid_count=2)["planning_envelope"]["variants"]
        self.assertEqual(len(variants), 2)
        self.assertEqual(
            {variant["cell_count"] for variant in variants},
            {16_640_000},
        )
        self.assertTrue(
            all(
                "noncanonical_grid_count_conservative" in variant["projection_method"]
                for variant in variants
            )
        )

    def test_duplicate_parameter_fails_closed(self) -> None:
        deck = self._temporary_deck(
            lambda text: text.replace(
                "tlim       = 1200.0",
                "tlim       = 1200.0\n"
                "tlim       = 1200.0",
            )
        )
        with self.assertRaisesRegex(storage.EstimatorError, "duplicate time/tlim"):
            storage.estimate_storage(deck)

    def test_malformed_numeric_value_fails_closed(self) -> None:
        deck = self._temporary_deck(
            lambda text: text.replace(
                "ps_eta                        = 1.0e-3",
                "ps_eta                        = nan",
            )
        )
        with self.assertRaisesRegex(storage.EstimatorError, "finite decimal"):
            storage.estimate_storage(deck)

    def test_unsupported_initial_particles_fail_closed(self) -> None:
        deck = self._temporary_deck(
            lambda text: text.replace(
                "ppc                               = 0.0",
                "ppc                               = 1.0",
            )
        )
        with self.assertRaisesRegex(storage.EstimatorError, "initial particles"):
            storage.estimate_storage(deck)

    def test_nondivisible_pvtk_cadence_fails_closed(self) -> None:
        deck = self._temporary_deck(
            lambda text: text.replace(
                "<output5>\nfile_type   = pvtk\nvariable    = prtcl_all\n"
                "id          = prtcl_all\ndt          = 100.0",
                "<output5>\nfile_type   = pvtk\nvariable    = prtcl_all\n"
                "id          = prtcl_all\ndt          = 128.0",
            )
        )
        with self.assertRaisesRegex(storage.EstimatorError, "divide time/tlim"):
            storage.estimate_storage(deck)

    def test_nondivisible_restart_cadence_fails_closed(self) -> None:
        deck = self._temporary_deck(
            lambda text: text.replace(
                "<output6>\nfile_type   = rst\ndt          = 100.0",
                "<output6>\nfile_type   = rst\ndt          = 128.0",
            )
        )
        with self.assertRaisesRegex(storage.EstimatorError, "restart cadence"):
            storage.estimate_storage(deck)

    def test_missing_required_mesh_bin_product_fails_closed(self) -> None:
        deck = self._temporary_deck(
            lambda text: text.replace(
                "<output4>\nfile_type   = bin\nvariable    = mhd_j2\n"
                "id          = j2\ndt          = 100.0\nghost_zones = false\n",
                "",
            )
        )
        with self.assertRaisesRegex(storage.EstimatorError, "required Section 5.4"):
            storage.estimate_storage(deck)

    def test_readiness_sidecar_matches_current_projection(self) -> None:
        frozen = json.loads(READINESS.read_text(encoding="utf-8"))
        report = storage.estimate_storage()
        projection = frozen["current_projection"]
        self.assertEqual(projection["grid_count"], storage.DEFAULT_GRID_COUNT)
        self.assertEqual(projection["seed_count"], storage.DEFAULT_SEED_COUNT)
        self.assertAlmostEqual(
            projection["retained_particle_count_at_t1200_analytical_approx"],
            _approx_quantity(
                report["per_run"]["retained_particle_count_at_tlim_analytical"]
            ),
        )
        self.assertAlmostEqual(
            projection[
                "raw_pvtk_payload_bytes_per_run_before_overhead_analytical_approx"
            ],
            _approx_quantity(
                report["per_run"]["raw_pvtk_payload_bytes_before_overhead"]
            ),
        )
        self.assertAlmostEqual(
            projection[
                "raw_pvtk_payload_bytes_campaign_before_overhead_analytical_approx"
            ],
            _approx_quantity(
                report["campaign"]["raw_pvtk_payload_bytes_before_overhead"]
            ),
        )
        campaign = report["planning_envelope"]["campaign"]
        self.assertAlmostEqual(
            projection["logical_campaign_bytes_before_filesystem_overhead_approx"],
            _approx_quantity(campaign["logical_bytes_before_filesystem_overhead"]),
        )
        self.assertAlmostEqual(
            projection["reservation_envelope_bytes_approx"],
            _approx_quantity(campaign["reservation_envelope_bytes"]),
        )
        self.assertAlmostEqual(
            projection["reservation_envelope_tb_decimal"],
            campaign["reservation_envelope_tb_decimal"],
        )
        for binding in frozen["source_bindings"].values():
            self.assertEqual(
                binding["sha256"],
                _sha256_path(storage.REPO_ROOT / binding["path"]),
            )


if __name__ == "__main__":
    unittest.main()
