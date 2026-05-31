#!/usr/bin/env python3
"""Focused formula tests for the Q008 expanding-box CPAW preparation matrix."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tst.scripts.particles import pic_mhd_expanding_box_cpaw_history_preparation as q008


_EXPECTED_AXES = {
    "x1": {
        "index": 0,
        "coordinate": "x1v",
        "phase_shape": (1, 1, -1),
        "parallel_b": "bcc1",
        "velocity_transverse": ("vely", "velz"),
        "magnetic_transverse": ("bcc2", "bcc3"),
        "history_parallel_me": 10,
        "history_transverse_ke": (8, 9),
        "history_transverse_me": (11, 12),
    },
    "x2": {
        "index": 1,
        "coordinate": "x2v",
        "phase_shape": (1, -1, 1),
        "parallel_b": "bcc2",
        "velocity_transverse": ("velz", "velx"),
        "magnetic_transverse": ("bcc3", "bcc1"),
        "history_parallel_me": 11,
        "history_transverse_ke": (9, 7),
        "history_transverse_me": (12, 10),
    },
    "x3": {
        "index": 2,
        "coordinate": "x3v",
        "phase_shape": (-1, 1, 1),
        "parallel_b": "bcc3",
        "velocity_transverse": ("velx", "vely"),
        "magnetic_transverse": ("bcc1", "bcc2"),
        "history_parallel_me": 12,
        "history_transverse_ke": (7, 8),
        "history_transverse_me": (10, 11),
    },
}


def _spec(
    axis: str = "x1",
    law: str = "exponential",
    rate: float = 0.01,
    resolution: int = 19,
    solver: str = "llf",
) -> dict[str, object]:
    return q008._spec("formula_probe", law, rate, axis, resolution, solver)


def _expected_scale_factor(
    law: str, rate: float, times: np.ndarray
) -> np.ndarray:
    if law == "linear":
        return 1.0 + rate * times
    if law == "reciprocal_linear":
        return 1.0 / (1.0 + rate * times)
    if law == "exponential":
        return np.exp(rate * times)
    raise AssertionError("unsupported test law " + law)


def _expected_geometry(
    law: str, rate: float, times: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    return (
        _expected_scale_factor(law, 2.0 * rate, times),
        _expected_scale_factor(law, rate, times),
    )


def _exact_history(spec: dict[str, object], times: np.ndarray) -> np.ndarray:
    parallel, transverse = _expected_geometry(
        str(spec["law"]), float(spec["rate"]), times
    )
    metadata = _EXPECTED_AXES[str(spec["axis"])]
    rows = np.zeros((times.size, 13))
    rows[:, 0] = times
    rows[:, 2] = 7.0

    transverse_ke = metadata["history_transverse_ke"]
    transverse_me = metadata["history_transverse_me"]
    parallel_me = metadata["history_parallel_me"]
    rows[:, transverse_ke[0]] = 3.0 * transverse**-2
    rows[:, transverse_ke[1]] = 5.0 * transverse**-2
    rows[:, transverse_me[0]] = 11.0 * parallel**-1
    rows[:, transverse_me[1]] = 13.0 * parallel**-1
    rows[:, parallel_me] = 17.0 * parallel / transverse**2
    return rows


class PicMHDExpandingBoxCPAWHistoryPreparationTests(unittest.TestCase):
    def test_axis_permutations_drive_dimensions_rates_and_command_flags(self) -> None:
        self.assertEqual(q008._AXES, _EXPECTED_AXES)
        expected = {
            "x1": ([19, 4, 4], [0.25, 0.125, 0.125]),
            "x2": ([4, 19, 4], [0.125, 0.25, 0.125]),
            "x3": ([4, 4, 19], [0.125, 0.125, 0.25]),
        }
        for axis, (dimensions, rates) in expected.items():
            with self.subTest(axis=axis):
                spec = _spec(axis=axis, law="linear", rate=0.125, solver="hlle")
                self.assertEqual(q008._dimensions(spec), dimensions)
                self.assertEqual(q008._rates(spec), rates)
                command = q008._command(spec)
                self.assertIn("mhd/rsolver=hlle", command)
                self.assertIn("particles/pic_expansion_law=linear", command)
                for candidate in _EXPECTED_AXES:
                    self.assertIn(
                        "problem/along_{}={}".format(
                            candidate, str(candidate == axis).lower()
                        ),
                        command,
                    )

    def test_mode_projection_uses_axis_specific_coordinate_and_transverse_pair(
        self,
    ) -> None:
        shapes = {"x1": (3, 2, 8), "x2": (3, 8, 2), "x3": (8, 3, 2)}
        coordinates = np.arange(8, dtype=np.float64) / 8.0
        for axis, shape in shapes.items():
            with self.subTest(axis=axis):
                metadata = _EXPECTED_AXES[axis]
                wave = np.exp(-1j * q008._WAVENUMBER * coordinates).reshape(
                    metadata["phase_shape"]
                )
                wave = np.broadcast_to(wave, shape)
                first, second = metadata["velocity_transverse"]
                snapshot = {
                    metadata["coordinate"]: coordinates,
                    first: wave.real,
                    second: wave.imag,
                }
                coefficient = q008._mode_coefficient(
                    snapshot, axis, metadata["velocity_transverse"]
                )
                self.assertAlmostEqual(coefficient.real, 1.0)
                self.assertAlmostEqual(coefficient.imag, 0.0)

    def test_geometry_scaling_uses_parallel_double_rate_for_each_law(self) -> None:
        times = np.asarray([0.0, 0.75, 2.0, 4.0])
        for law in ("linear", "reciprocal_linear", "exponential"):
            for rate in (-0.03, 0.03):
                with self.subTest(law=law, rate=rate):
                    spec = _spec(law=law, rate=rate)
                    parallel, transverse = q008._geometry(spec, times)
                    expected_parallel, expected_transverse = _expected_geometry(
                        law, rate, times
                    )
                    np.testing.assert_allclose(parallel, expected_parallel)
                    np.testing.assert_allclose(transverse, expected_transverse)
                    np.testing.assert_allclose(
                        q008._scale_factor(law, rate, times), expected_transverse
                    )
        with self.assertRaisesRegex(ValueError, "Unsupported expansion law"):
            q008._scale_factor("quadratic", 0.03, times)

    def test_phase_oracle_matches_independent_geometry_quadrature(self) -> None:
        time = 4.0
        sample_times = np.linspace(0.0, time, 100001)
        step = sample_times[1] - sample_times[0]
        for law in ("linear", "reciprocal_linear", "exponential"):
            for rate in (-0.03, 0.03):
                with self.subTest(law=law, rate=rate):
                    parallel, transverse = _expected_geometry(
                        law, rate, sample_times
                    )
                    integrand = 1.0 / (transverse * np.sqrt(parallel))
                    integral = step * (
                        0.5 * integrand[0]
                        + np.sum(integrand[1:-1])
                        + 0.5 * integrand[-1]
                    )
                    expected = q008._WAVENUMBER * q008._ALFVEN_SPEED * integral
                    self.assertAlmostEqual(
                        q008._phase_oracle(law, rate, time), expected, places=10
                    )
            self.assertAlmostEqual(
                q008._phase_oracle(law, 0.0, time),
                q008._WAVENUMBER * q008._ALFVEN_SPEED * time,
            )
        with self.assertRaisesRegex(ValueError, "Unsupported expansion law"):
            q008._phase_oracle("quadratic", 0.03, time)

    def test_history_energy_scaling_is_axis_permuted_and_law_independent(self) -> None:
        times = np.asarray([0.0, 0.75, 2.0, 4.0])
        for axis in _EXPECTED_AXES:
            for law in ("linear", "reciprocal_linear", "exponential"):
                with self.subTest(axis=axis, law=law):
                    spec = _spec(axis=axis, law=law, rate=0.03)
                    history = _exact_history(spec, times)
                    metrics = q008._history_oracle_metrics(history, spec)
                    for value in metrics.values():
                        self.assertAlmostEqual(value, 0.0)

                    transverse_me = _EXPECTED_AXES[axis]["history_transverse_me"][0]
                    history[-1, transverse_me] *= 1.05
                    metrics = q008._history_oracle_metrics(history, spec)
                    self.assertGreater(
                        metrics[
                            "transverse_magnetic_energy_scaled_relative_error_max"
                        ],
                        0.0,
                    )

    def test_roe_probe_requires_nonzero_rejection_text_and_no_outputs(self) -> None:
        def probe(
            returncode: int, output: str, produce_output: bool = False
        ) -> dict[str, object]:
            with tempfile.TemporaryDirectory() as directory:
                run_dir = Path(directory)

                def fake_run(command, cwd, capture_output, text):
                    self.assertTrue(capture_output)
                    self.assertTrue(text)
                    self.assertEqual(Path(cwd), run_dir)
                    self.assertIn("mhd/rsolver=roe", command)
                    self.assertIn("time/nlim=0", command)
                    if produce_output:
                        spec = q008._spec(
                            "roe_unsupported_probe", "exponential", 0.0, "x1", 32, "roe"
                        )
                        (run_dir / (q008._basename(spec) + ".mhd.hst")).write_text(
                            "unexpected\n", encoding="utf-8"
                        )
                    return subprocess.CompletedProcess(
                        command, returncode, stdout="", stderr=output
                    )

                with mock.patch.object(
                    q008, "_athena_run_dir", return_value=str(run_dir)
                ), mock.patch.object(q008.subprocess, "run", side_effect=fake_run):
                    result = q008._run_roe_rejection_probe()
                sidecar_dir = run_dir / "roe_unsupported_probe"
                self.assertTrue((sidecar_dir / "command.json").is_file())
                self.assertEqual(
                    (sidecar_dir / "returncode.txt").read_text(encoding="utf-8"),
                    str(returncode) + "\n",
                )
                return result

        accepted = probe(1, q008._ROE_REJECTION_TEXT)
        self.assertEqual(accepted["status"], "pass_fail_closed_unsupported")
        self.assertTrue(accepted["expected_rejection_text_observed"])
        self.assertEqual(accepted["outputs_after_rejection"], [])

        for returncode, output, produce_output in (
            (0, q008._ROE_REJECTION_TEXT, False),
            (1, "different failure", False),
            (1, q008._ROE_REJECTION_TEXT, True),
        ):
            with self.subTest(
                returncode=returncode,
                output=output,
                produce_output=produce_output,
            ):
                rejected = probe(returncode, output, produce_output)
                self.assertEqual(
                    rejected["status"], "fail_unexpected_roe_disposition"
                )

    def test_roe_is_not_a_successful_matrix_solver(self) -> None:
        self.assertEqual(q008._SOLVERS, ("llf", "hlle", "hlld"))
        self.assertNotIn("roe", {spec["solver"] for spec in q008._successful_specs()})


if __name__ == "__main__":
    unittest.main()
