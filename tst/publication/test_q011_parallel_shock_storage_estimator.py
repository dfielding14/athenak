#!/usr/bin/env python3
"""Focused tests for the deterministic Q-011 storage estimator sidecar."""

from __future__ import annotations

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
    / "q011_parallel_shock_storage_estimator_successor_v2_2026-06-01.json"
)


def _approx_quantity(value):
    if isinstance(value, int):
        return float(value)
    return value["approx"]


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
        assumptions = storage.estimate_storage()["assumptions"]
        self.assertIn("no particle escape", assumptions["particle_retention"])
        self.assertIn("t=0", assumptions["cadence"])
        self.assertIn("timestep-edge", assumptions["quantization"])
        self.assertIn("ASCII headers", assumptions["payload_scope"])
        self.assertIn("restart files", assumptions["payload_scope"])

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


if __name__ == "__main__":
    unittest.main()
