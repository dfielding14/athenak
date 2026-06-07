#!/usr/bin/env python3
"""Unit tests for the Ito turbulence validation scripts."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from analyze_ito_turbulence_validation import (
    cic_deposit_periodic,
    density_metrics,
    isotropic_power_spectrum,
    jensen_shannon_divergence,
    read_vtk_points,
)
from aggregate_ito_turbulence_validation import aggregate_summaries, write_csv
from run_ito_turbulence_validation import command_for_run


class ItoTurbulenceAnalysisTests(unittest.TestCase):
    def test_cic_deposit_is_uniform_at_cell_centers(self):
        size = 4
        centers = (np.arange(size) + 0.5) / size
        z, y, x = np.meshgrid(centers, centers, centers, indexing="ij")
        points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
        counts = cic_deposit_periodic(
            points,
            (size, size, size),
            ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        )
        np.testing.assert_allclose(counts, 1.0, rtol=0.0, atol=1.0e-14)
        self.assertAlmostEqual(float(np.sum(counts)), float(points.shape[0]))

    def test_identical_density_metrics_are_exact(self):
        density = np.linspace(0.5, 1.5, 64).reshape(4, 4, 4)
        metrics = density_metrics(density, density.copy())
        self.assertAlmostEqual(metrics["pearson_r"], 1.0)
        self.assertAlmostEqual(metrics["r2_one_to_one"], 1.0)
        self.assertAlmostEqual(metrics["normalized_l2"], 0.0)
        self.assertAlmostEqual(metrics["log10_ratio_std"], 0.0)

    def test_power_spectrum_recovers_single_mode(self):
        size = 16
        x = np.arange(size)
        wave = 1.0 + 0.1 * np.sin(2.0 * np.pi * 2.0 * x / size)
        field = np.broadcast_to(wave[None, None, :], (size, size, size))
        k, power = isotropic_power_spectrum(field)
        self.assertEqual(int(k[np.argmax(power)]), 2)

    def test_jensen_shannon_is_zero_for_identical_histograms(self):
        distribution = np.array([0.1, 0.2, 0.7])
        self.assertAlmostEqual(jensen_shannon_divergence(distribution, distribution), 0.0)

    def test_particle_vtk_point_reader(self):
        points = np.array([[0.1, 0.2, 0.3], [0.9, 0.8, 0.7]], dtype=np.float32)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "particles.vtk"
            with path.open("wb") as stream:
                stream.write(b"# vtk DataFile Version 2.0\n")
                stream.write(
                    b"# AthenaK particle data at time= 0.125 nranks= 1 cycle=2\n"
                )
                stream.write(b"BINARY\nDATASET UNSTRUCTURED_GRID\n\n")
                stream.write(b"POINTS 2 float\n")
                stream.write(points.astype(">f4").tobytes())
            loaded, time_value = read_vtk_points(path)
        np.testing.assert_allclose(loaded, points)
        self.assertAlmostEqual(time_value, 0.125)

    def test_command_builder_preserves_overrides(self):
        command = command_for_run(
            Path("/tmp/athena"),
            Path("/tmp/test.athinput"),
            Path("/tmp/run"),
            "mpirun -np 2",
            ["time/cfl_number=0.32"],
        )
        self.assertEqual(command[:3], ["mpirun", "-np", "2"])
        self.assertEqual(command[-1], "time/cfl_number=0.32")

    def test_ensemble_aggregation_uses_paired_standard_error(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for index, value in enumerate((1.0, 2.0, 3.0)):
                path = Path(directory) / f"summary-{index}.json"
                path.write_text(
                    (
                        '{"comparison":{"delta_metric":'
                        f"{value}"
                        '},"old":{"runtime":{"elapsed_seconds":2.0}},'
                        '"corrected":{"runtime":{"elapsed_seconds":3.0}}}'
                    ),
                    encoding="utf-8",
                )
                paths.append(path)
            aggregate = aggregate_summaries(paths)

        self.assertEqual(aggregate["seed_count"], 3)
        self.assertAlmostEqual(
            aggregate["metrics"]["delta_metric"]["mean"], 2.0
        )
        self.assertAlmostEqual(
            aggregate["metrics"]["delta_metric"]["standard_error"],
            1.0 / np.sqrt(3.0),
        )
        self.assertAlmostEqual(aggregate["runtime_wall_ratio"]["mean"], 1.5)

    def test_ensemble_csv_uses_lf_records(self):
        aggregate = {
            "metrics": {
                "delta_metric": {
                    "mean": 1.0,
                    "standard_deviation": 0.5,
                    "standard_error": 0.25,
                }
            },
            "runtime_wall_ratio": {"mean": 1.5, "standard_error": 0.1},
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "aggregate.csv"
            write_csv(path, aggregate)
            payload = path.read_bytes()
        self.assertNotIn(b"\r\n", payload)
        self.assertEqual(payload.count(b"\n"), 3)


if __name__ == "__main__":
    unittest.main()
