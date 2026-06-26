#!/usr/bin/env python3

import math
import unittest

from tst.publication import bell_saturation_engineering_analysis_v1 as analysis


def _rows(amplitudes: list[float]) -> list[dict[str, float]]:
    return [
        {"tau": float(index), "bperp_rms_over_B0": amplitude}
        for index, amplitude in enumerate(amplitudes)
    ]


class SaturationWindowTests(unittest.TestCase):
    def test_sustained_plateau_passes(self) -> None:
        amplitudes = [1.0e-3 * math.exp(0.7 * index) for index in range(10)]
        amplitudes.extend(2.0 * (1.0 + 0.02 * math.sin(index)) for index in range(24))
        rows = _rows(amplitudes)
        onset = next(index for index, row in enumerate(rows) if row["bperp_rms_over_B0"] >= 1.0)
        candidate = analysis._candidate_window(rows, onset)
        self.assertIsNotNone(candidate)
        self.assertIsNotNone(analysis._confirmation_window(rows, candidate))

    def test_large_post_peak_oscillation_is_not_saturation(self) -> None:
        amplitudes = [1.0e-3 * math.exp(0.8 * index) for index in range(10)]
        amplitudes.extend([2.0, 1.8, 1.5, 0.8, 0.2, 0.08, 0.3, 0.9, 1.7, 2.1] * 3)
        rows = _rows(amplitudes)
        onset = next(index for index, row in enumerate(rows) if row["bperp_rms_over_B0"] >= 1.0)
        candidate = analysis._candidate_window(rows, onset)
        if candidate is not None:
            self.assertIsNone(analysis._confirmation_window(rows, candidate))


if __name__ == "__main__":
    unittest.main()
