#!/usr/bin/env python3
"""Synthetic tests for the strict Q-011 Section 5.4 analyzer foundation."""

from __future__ import annotations

import math
import struct
import unittest

import numpy as np

try:
    from tst.publication import analyze_q011_section54_outputs as q011
except ModuleNotFoundError:
    import analyze_q011_section54_outputs as q011


_DOMAIN_BOUNDS = (0.0, 4.0, 0.0, 2.0, 0.0, 1.0)
_ROOT_SHAPE = (4, 2, 1)
_BLOCK_SHAPE = (2, 2, 1)


def _geometry(
    logical_location: tuple[int, int, int],
    level: int,
    *,
    root_shape: tuple[int, int, int] = _ROOT_SHAPE,
    block_shape: tuple[int, int, int] = _BLOCK_SHAPE,
    domain_bounds: tuple[float, float, float, float, float, float] = _DOMAIN_BOUNDS,
) -> tuple[float, float, float, float, float, float]:
    geometry = []
    for axis, logical_index in enumerate(logical_location):
        root_blocks = root_shape[axis] // block_shape[axis]
        logical_extent = root_blocks * (2**level if root_shape[axis] > 1 else 1)
        lower = domain_bounds[2 * axis]
        upper = domain_bounds[2 * axis + 1]
        width = (upper - lower) / logical_extent
        geometry.extend((lower + logical_index * width, lower + (logical_index + 1) * width))
    return tuple(geometry)


def _block(
    logical_location: tuple[int, int, int],
    level: int,
    value: float,
    *,
    fields: tuple[str, ...] = ("dens",),
    root_shape: tuple[int, int, int] = _ROOT_SHAPE,
    block_shape: tuple[int, int, int] = _BLOCK_SHAPE,
    domain_bounds: tuple[float, float, float, float, float, float] = _DOMAIN_BOUNDS,
    geometry: tuple[float, float, float, float, float, float] | None = None,
) -> bytes:
    nx1, nx2, nx3 = block_shape
    index_and_logical = (
        0,
        nx1 - 1,
        0,
        nx2 - 1,
        0,
        nx3 - 1,
        *logical_location,
        level,
    )
    resolved_geometry = geometry or _geometry(
        logical_location,
        level,
        root_shape=root_shape,
        block_shape=block_shape,
        domain_bounds=domain_bounds,
    )
    cell_count = nx1 * nx2 * nx3
    values = np.concatenate(
        [
            np.full(cell_count, value + index, dtype="<f4")
            for index, _ in enumerate(fields)
        ]
    )
    return (
        struct.pack("<10i", *index_and_logical)
        + struct.pack("<6d", *resolved_geometry)
        + values.tobytes()
    )


def _payload(
    blocks: list[bytes],
    *,
    fields: tuple[str, ...] = ("dens",),
    root_shape: tuple[int, int, int] = _ROOT_SHAPE,
    block_shape: tuple[int, int, int] = _BLOCK_SHAPE,
    domain_bounds: tuple[float, float, float, float, float, float] = _DOMAIN_BOUNDS,
) -> bytes:
    x1min, x1max, x2min, x2max, x3min, x3max = domain_bounds
    parameter_header = (
        "<mesh>\n"
        f"nx1={root_shape[0]}\n"
        f"nx2={root_shape[1]}\n"
        f"nx3={root_shape[2]}\n"
        "nghost=0\n"
        f"x1min={x1min}\n"
        f"x1max={x1max}\n"
        f"x2min={x2min}\n"
        f"x2max={x2max}\n"
        f"x3min={x3min}\n"
        f"x3max={x3max}\n"
        "<meshblock>\n"
        f"nx1={block_shape[0]}\n"
        f"nx2={block_shape[1]}\n"
        f"nx3={block_shape[2]}\n"
    ).encode()
    header = (
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        b"  time=500.0\n"
        b"  cycle=42\n"
        b"  size of location=8\n"
        b"  size of variable=4\n"
        + f"  number of variables={len(fields)}\n".encode()
        + b"  variables:  "
        + b"  ".join(name.encode() for name in fields)
        + b"  \n"
        + f"  header offset={len(parameter_header)}\n".encode()
        + parameter_header
    )
    return header + b"".join(blocks)


def _uniform_payload() -> bytes:
    return _payload(
        [
            _block((0, 0, 0), 0, 1.0),
            _block((1, 0, 0), 0, 3.0),
        ]
    )


class Q011Section54BinaryTests(unittest.TestCase):
    def test_parser_reads_structural_dataset_and_composes_complete_grid(self) -> None:
        dataset = q011.parse_athenak_binary_bytes(_uniform_payload())
        self.assertEqual(dataset.time, 500.0)
        self.assertEqual(dataset.cycle, 42)
        self.assertEqual(dataset.root_block_shape, (2, 1, 1))
        self.assertEqual(len(dataset.blocks), 2)
        composite = q011.compose_leaf_field(dataset, "dens")
        self.assertEqual(composite.target_level, 0)
        np.testing.assert_array_equal(
            composite.values,
            np.array([[[1.0, 1.0, 3.0, 3.0], [1.0, 1.0, 3.0, 3.0]]]),
        )
        np.testing.assert_array_equal(composite.source_levels, np.zeros((1, 2, 4)))

    def test_parser_rejects_malformed_and_truncated_inputs(self) -> None:
        payload = _uniform_payload()
        with self.assertRaisesRegex(q011.AnalysisError, "bad magic"):
            q011.parse_athenak_binary_bytes(b"Not Athena\n" + payload.split(b"\n", 1)[1])
        with self.assertRaisesRegex(q011.AnalysisError, "truncated block 1 field data"):
            q011.parse_athenak_binary_bytes(payload[:-1])
        malformed = payload.replace(
            b"  number of variables=1\n", b"  number of variables=2\n", 1
        )
        with self.assertRaisesRegex(q011.AnalysisError, "variable count does not match"):
            q011.parse_athenak_binary_bytes(malformed)

    def test_parser_rejects_nonfinite_field_geometry_and_invalid_logical_level(self) -> None:
        nonfinite_field = _payload([_block((0, 0, 0), 0, math.nan)])
        with self.assertRaisesRegex(q011.AnalysisError, "field data is not finite"):
            q011.parse_athenak_binary_bytes(nonfinite_field)
        geometry = list(_geometry((0, 0, 0), 0))
        geometry[0] = math.inf
        nonfinite_geometry = _payload([_block((0, 0, 0), 0, 1.0, geometry=tuple(geometry))])
        with self.assertRaisesRegex(q011.AnalysisError, "geometry is not finite"):
            q011.parse_athenak_binary_bytes(nonfinite_geometry)
        negative_level = _payload([_block((0, 0, 0), -1, 1.0)])
        with self.assertRaisesRegex(q011.AnalysisError, "invalid physical refinement level"):
            q011.parse_athenak_binary_bytes(negative_level)
        excessive_level = _payload([_block((0, 0, 0), 31, 1.0)])
        with self.assertRaisesRegex(q011.AnalysisError, "invalid physical refinement level"):
            q011.parse_athenak_binary_bytes(excessive_level)

    def test_parser_rejects_duplicate_leaves(self) -> None:
        duplicate = _payload(
            [
                _block((0, 0, 0), 0, 1.0),
                _block((0, 0, 0), 0, 2.0),
            ]
        )
        with self.assertRaisesRegex(q011.AnalysisError, "duplicate leaf"):
            q011.parse_athenak_binary_bytes(duplicate)

    def test_composer_rejects_holes_and_coarse_fine_overlap(self) -> None:
        hole = q011.parse_athenak_binary_bytes(_payload([_block((0, 0, 0), 0, 1.0)]))
        with self.assertRaisesRegex(q011.AnalysisError, "holes"):
            q011.compose_leaf_field(hole, "dens")

        overlap_blocks = [_block((0, 0, 0), 0, 1.0)]
        overlap_blocks.extend(
            _block((lx1, lx2, 0), 1, 2.0)
            for lx2 in range(2)
            for lx1 in range(4)
        )
        overlap = q011.parse_athenak_binary_bytes(_payload(overlap_blocks))
        with self.assertRaisesRegex(q011.AnalysisError, "overlap"):
            q011.compose_leaf_field(overlap, "dens")

    def test_rank_datasets_merge_before_complete_composition(self) -> None:
        left = q011.parse_athenak_binary_bytes(
            _payload([_block((0, 0, 0), 0, 1.0)]), source="rank0"
        )
        right = q011.parse_athenak_binary_bytes(
            _payload([_block((1, 0, 0), 0, 3.0)]), source="rank1"
        )
        merged = q011.merge_athenak_binary_datasets([left, right])
        composite = q011.compose_leaf_field(merged, "dens")
        np.testing.assert_array_equal(
            composite.values,
            np.array([[[1.0, 1.0, 3.0, 3.0], [1.0, 1.0, 3.0, 3.0]]]),
        )


class Q011Section54ArithmeticTests(unittest.TestCase):
    def test_conservative_area_restriction_preserves_integral(self) -> None:
        values = np.array([[1.0, 3.0], [5.0, 7.0]])
        areas = np.array([[1.0, 1.0], [1.0, 3.0]])
        restricted = q011.restrict_area_average(values, cell_areas=areas)
        np.testing.assert_allclose(restricted.values, [[5.0]])
        np.testing.assert_allclose(restricted.areas, [[6.0]])
        self.assertAlmostEqual(
            float(np.sum(values * areas)),
            float(np.sum(restricted.values * restricted.areas)),
        )

    def test_fixed_histogram_tracks_bins_weights_and_out_of_range_samples(self) -> None:
        histogram = q011.fixed_histogram(
            [-1.0, 0.0, 0.5, 1.0, 2.0, 3.0],
            [0.0, 1.0, 2.0],
            weights=[10.0, 1.0, 2.0, 3.0, 4.0, 20.0],
        )
        np.testing.assert_array_equal(histogram.counts, [2, 2])
        np.testing.assert_allclose(histogram.weighted_counts, [3.0, 7.0])
        self.assertEqual(histogram.underflow_count, 1)
        self.assertEqual(histogram.overflow_count, 1)
        self.assertEqual(histogram.underflow_weight, 10.0)
        self.assertEqual(histogram.overflow_weight, 20.0)
        empty = q011.fixed_histogram([], [0.0, 1.0, 2.0])
        np.testing.assert_array_equal(empty.counts, [0, 0])
        np.testing.assert_array_equal(empty.weighted_counts, [0.0, 0.0])

    def test_fixed_bin_loglog_fit_recovers_energy_negative_three_halves(self) -> None:
        edges = np.geomspace(1.0, 256.0, 9)
        centers = np.sqrt(edges[:-1] * edges[1:])
        spectrum = 7.0 * centers ** (-1.5)
        fit = q011.fit_fixed_bin_loglog_slope(
            edges, spectrum, fit_window=(centers[0], centers[-1])
        )
        self.assertAlmostEqual(fit.slope, -1.5, places=12)
        self.assertEqual(fit.selected_bin_indices, tuple(range(8)))

    def test_fixed_bin_loglog_fit_rejects_insufficient_positive_bins(self) -> None:
        with self.assertRaisesRegex(q011.AnalysisError, "insufficient positive fit bins"):
            q011.fit_fixed_bin_loglog_slope(
                [1.0, 2.0, 4.0, 8.0],
                [1.0, 0.0, 0.0],
                fit_window=(1.0, 8.0),
            )

    def test_shock_front_search_ignores_stronger_distractor_outside_window(self) -> None:
        front = q011.detect_shock_front(
            np.arange(10.0),
            [1.0, 1.0, 1.0, 1.0, 2.0, 4.0, 4.0, 4.0, 20.0, 20.0],
            search_window=(3.0, 6.0),
            gradient_sign="positive",
        )
        self.assertEqual(front.index, 4)
        self.assertEqual(front.x1, 4.0)
        self.assertEqual(front.density_gradient, 1.5)

    def test_upstream_amplification_uses_only_explicit_window(self) -> None:
        estimate = q011.estimate_upstream_amplification(
            [0.0, 1.0, 2.0, 3.0],
            [[100.0, 4.0, 6.0, 100.0], [100.0, 2.0, 8.0, 100.0]],
            upstream_window=(1.0, 2.0),
            reference_magnetic_field=2.0,
        )
        self.assertEqual(estimate.selected_cell_count, 4)
        self.assertEqual(estimate.selected_area, 4.0)
        self.assertEqual(estimate.mean_magnetic_magnitude, 5.0)
        self.assertEqual(estimate.amplification, 2.5)

    def test_matched_grid_residual_metrics_cover_zero_and_nonzero_differences(self) -> None:
        zero = q011.matched_grid_residual_metrics([1.0, 2.0], [1.0, 2.0])
        self.assertEqual(zero.mean_absolute, 0.0)
        self.assertEqual(zero.root_mean_square, 0.0)
        self.assertEqual(zero.maximum_absolute, 0.0)

        residual = q011.matched_grid_residual_metrics(
            [1.0, 2.0], [2.0, 4.0], cell_areas=[1.0, 3.0]
        )
        self.assertEqual(residual.mean_absolute, 1.75)
        self.assertAlmostEqual(residual.root_mean_square, math.sqrt(3.25))
        self.assertEqual(residual.maximum_absolute, 2.0)
        self.assertEqual(residual.relative_mean_absolute, 1.0)
        self.assertEqual(residual.relative_root_mean_square, 1.0)
        self.assertEqual(residual.relative_maximum_absolute, 1.0)


if __name__ == "__main__":
    unittest.main()
