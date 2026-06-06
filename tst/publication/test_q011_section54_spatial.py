#!/usr/bin/env python3
"""Focused synthetic tests for the Q-011 Section 5.4 spatial reducers."""

from __future__ import annotations

from dataclasses import replace
import unittest

import numpy as np

try:
    from tst.publication import analyze_q011_section54_outputs as output_primitives
    from tst.publication import q011_section54_spatial as spatial
except ModuleNotFoundError:
    import analyze_q011_section54_outputs as output_primitives
    import q011_section54_spatial as spatial


def _dataset(
    field: str,
    values: object,
    *,
    time: float = 500.0,
    domain_bounds: tuple[float, float, float, float, float, float] = (
        0.0,
        2000.0,
        0.0,
        2.0,
        0.0,
        1.0,
    ),
) -> output_primitives.AthenaBinaryDataset:
    array = np.asarray(values, dtype=np.float64)
    nx3, nx2, nx1 = array.shape
    block = output_primitives.AthenaBinaryBlock(
        index_bounds=(0, nx1 - 1, 0, nx2 - 1, 0, nx3 - 1),
        logical_location=(0, 0, 0),
        level=0,
        geometry=domain_bounds,
        fields={field: array},
    )
    return output_primitives.AthenaBinaryDataset(
        source=f"synthetic-{field}",
        time=time,
        cycle=50,
        location_size=8,
        variable_size=8,
        variable_names=(field,),
        input_parameters={},
        root_grid_shape=(nx1, nx2, nx3),
        meshblock_shape=(nx1, nx2, nx3),
        nghost=0,
        domain_bounds=domain_bounds,
        blocks=(block,),
    )


def _raster(
    quantity: str,
    values: object,
    *,
    nominal_slot_time: float = 500.0,
    observed_committed_time: float | None = None,
    x1_faces: object = (0.0, 1.0, 2.0),
    x2_faces: object = (0.0, 1.0, 2.0),
    source_levels: object | None = None,
    target_level: int = 0,
) -> spatial.CartesianXYRaster:
    array = np.asarray(values, dtype=np.float64)
    if observed_committed_time is None:
        observed_committed_time = nominal_slot_time
    if source_levels is None:
        source_levels = np.zeros_like(array, dtype=np.int64)
    return spatial.CartesianXYRaster(
        quantity=quantity,
        source_field=spatial.mesh_field_name(quantity),
        nominal_slot_time=nominal_slot_time,
        observed_committed_time=observed_committed_time,
        x1_faces_c_over_omega_pi=np.asarray(x1_faces, dtype=np.float64),
        x2_faces_c_over_omega_pi=np.asarray(x2_faces, dtype=np.float64),
        collapsed_x3_faces=np.array([0.0, 1.0]),
        values_y_x=array,
        source_levels_y_x=source_levels,
        target_level=target_level,
    )


def _profile(
    x1: object,
    values: object,
) -> spatial.YAreaWeightedProfile:
    array = np.asarray(values, dtype=np.float64)
    return spatial.YAreaWeightedProfile(
        quantity="rho",
        source_field="dens",
        nominal_slot_time=500.0,
        observed_committed_time=500.0,
        x1_centers_c_over_omega_pi=np.asarray(x1, dtype=np.float64),
        values_x=array,
        column_areas=np.ones_like(array),
    )


def _snapshot_datasets(
    *, embedded_time: float = 500.0
) -> dict[str, output_primitives.AthenaBinaryDataset]:
    rho_x = np.array([1.0, 1.0, 2.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0])
    return {
        "rho": _dataset("dens", np.tile(rho_x, (1, 2, 1)), time=embedded_time),
        "bmag": _dataset("bmag", np.full((1, 2, 10), 2.0), time=embedded_time),
        "prtcl_jx": _dataset(
            "prtcl_jx", np.full((1, 2, 10), -0.25), time=embedded_time
        ),
        "j2": _dataset("j2", np.full((1, 2, 10), 0.125), time=embedded_time),
    }


class Q011Section54SpatialTests(unittest.TestCase):
    def test_exact_mesh_quantity_mapping_and_singleton_x3_collapse(self) -> None:
        self.assertEqual(
            dict(spatial.MESH_QUANTITY_FIELDS),
            {
                "rho": "dens",
                "bmag": "bmag",
                "prtcl_jx": "prtcl_jx",
                "j2": "j2",
            },
        )
        with self.assertRaisesRegex(spatial.AnalysisError, "unknown mesh quantity"):
            spatial.mesh_field_name("density")

        raster = spatial.compose_xy_quantity(
            _dataset("dens", np.arange(12.0).reshape(1, 3, 4)), "rho"
        )
        self.assertEqual(raster.values_y_x.shape, (3, 4))
        np.testing.assert_array_equal(
            raster.values_y_x, np.arange(12.0).reshape(3, 4)
        )
        np.testing.assert_array_equal(raster.source_levels_y_x, np.zeros((3, 4)))

        wrong_field = replace(
            _dataset("dens", np.ones((1, 2, 2))), variable_names=("rho",)
        )
        with self.assertRaisesRegex(spatial.AnalysisError, "exactly field 'dens'"):
            spatial.compose_xy_quantity(wrong_field, "rho")
        with self.assertRaisesRegex(spatial.AnalysisError, "x3 must contain exactly one"):
            spatial.compose_xy_quantity(
                _dataset("dens", np.ones((2, 2, 2))), "rho"
            )

    def test_snapshot_composition_requires_exact_matched_scalar_products(self) -> None:
        datasets = _snapshot_datasets()
        rasters = spatial.compose_xy_snapshot(datasets)
        self.assertEqual(tuple(rasters), spatial.REQUIRED_MESH_QUANTITIES)

        missing = dict(datasets)
        missing.pop("j2")
        with self.assertRaisesRegex(spatial.AnalysisError, "must be exactly"):
            spatial.compose_xy_snapshot(missing)

        mismatched = dict(datasets)
        mismatched["j2"] = _dataset(
            "j2",
            np.full((1, 2, 8), 0.125),
            domain_bounds=(0.0, 2000.0, 0.0, 2.0, 0.0, 1.0),
        )
        with self.assertRaisesRegex(spatial.AnalysisError, "Cartesian grids disagree"):
            spatial.compose_xy_snapshot(mismatched)

    def test_y_area_weighted_profile_uses_cartesian_cell_areas(self) -> None:
        raster = _raster(
            "rho",
            [[1.0, 2.0], [3.0, 6.0]],
            x1_faces=(0.0, 1.0, 3.0),
            x2_faces=(0.0, 1.0, 4.0),
        )
        profile = spatial.y_area_weighted_profile(raster)
        np.testing.assert_allclose(profile.values_x, [2.5, 5.0])
        np.testing.assert_allclose(profile.column_areas, [4.0, 8.0])
        np.testing.assert_allclose(profile.x1_centers_c_over_omega_pi, [0.5, 2.0])

    def test_detected_front_uses_fixed_positive_window_and_offset_bound(self) -> None:
        x1 = np.arange(100.0, 3300.0, 200.0)
        density = np.array(
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 5.0]
            + [5.0] * 6
            + [100.0]
        )
        front = spatial.detect_positive_gradient_front(
            _profile(x1, density), x_ideal_c_over_omega_pi=1500.0
        )
        self.assertEqual(front.search_window_c_over_omega_pi, (300.0, 2700.0))
        self.assertEqual(front.x_front_c_over_omega_pi, 1500.0)
        self.assertEqual(front.offset_from_x_ideal_c_over_omega_pi, 0.0)

        offset_density = np.array(
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 5.0, 5.0, 5.0]
        )
        with self.assertRaisesRegex(spatial.AnalysisError, "maximum absolute offset"):
            spatial.detect_positive_gradient_front(
                _profile(np.arange(0.0, 2200.0, 200.0), offset_density),
                x_ideal_c_over_omega_pi=100.0,
            )

    def test_upstream_b_amplification_uses_fixed_t500_area_weighted_gate(self) -> None:
        raster = _raster(
            "bmag",
            [[100.0, 1.0, 1.0, 1.0, 100.0], [100.0, 3.0, 3.0, 3.0, 100.0]],
            x1_faces=(0.0, 120.0, 240.0, 720.0, 1200.0, 1320.0),
            x2_faces=(0.0, 1.0, 4.0),
        )
        amplification = spatial.reduce_upstream_b_amplification_at_t500(
            raster, x_ideal_c_over_omega_pi=0.0
        )
        self.assertEqual(amplification.upstream_window_c_over_omega_pi, (120.0, 1200.0))
        self.assertEqual(amplification.selected_cell_count, 6)
        self.assertEqual(amplification.selected_area, 4320.0)
        self.assertEqual(amplification.mean_magnetic_magnitude, 2.5)
        self.assertEqual(amplification.amplification_over_b0, 2.5)
        self.assertTrue(amplification.passes_gate)

        low = _raster("bmag", np.ones((2, 5)), x1_faces=(0, 120, 240, 720, 1200, 1320))
        self.assertFalse(
            spatial.reduce_upstream_b_amplification_at_t500(
                low, x_ideal_c_over_omega_pi=0.0
            ).passes_gate
        )
        with self.assertRaisesRegex(spatial.AnalysisError, "only for nominal t=500"):
            spatial.reduce_upstream_b_amplification_at_t500(
                replace(raster, nominal_slot_time=400.0),
                x_ideal_c_over_omega_pi=0.0,
            )

    def test_morphology_record_carries_source_level_overlay_metadata(self) -> None:
        raster = _raster(
            "j2",
            [[1.0, 2.0], [3.0, 4.0]],
            source_levels=[[0, 1], [1, 1]],
            target_level=1,
        )
        record = spatial.morphology_raster_record(raster)
        self.assertEqual(record["shape_y_x"], [2, 2])
        self.assertEqual(
            record["source_level_overlay"],
            {
                "encoding": "physical_refinement_level_per_cell",
                "target_composite_level": 1,
                "levels_y_x": [[0, 1], [1, 1]],
                "level_cell_counts": [
                    {"source_level": 0, "cell_count": 1},
                    {"source_level": 1, "cell_count": 3},
                ],
            },
        )

    def test_t500_snapshot_reducer_emits_profiles_front_gate_and_rasters(self) -> None:
        record = spatial.reduce_t500_spatial_snapshot(
            _snapshot_datasets(),
            nominal_slot_time=500.0,
            observed_committed_time=500.0,
            x_ideal_c_over_omega_pi=500.0,
        )
        self.assertEqual(record["record_type"], "q011_section54_t500_spatial_reduction")
        self.assertEqual(record["detected_front"]["x_front_c_over_omega_pi"], 500.0)
        self.assertTrue(record["upstream_b_amplification"]["passes_gate"])
        self.assertEqual(
            tuple(record["y_area_weighted_profiles"]),
            spatial.REQUIRED_MESH_QUANTITIES,
        )
        self.assertEqual(
            tuple(record["morphology_rasters"]), spatial.REQUIRED_MESH_QUANTITIES
        )

    def test_t500_snapshot_reducer_preserves_full_observed_committed_time(self) -> None:
        observed = 500.053496123
        record = spatial.reduce_t500_spatial_snapshot(
            _snapshot_datasets(embedded_time=500.053),
            nominal_slot_time=500.0,
            observed_committed_time=observed,
            x_ideal_c_over_omega_pi=500.0,
        )
        self.assertEqual(record["nominal_slot_time"], 500.0)
        self.assertEqual(record["observed_committed_time"], observed)
        self.assertEqual(
            record["upstream_b_amplification"]["observed_committed_time"],
            observed,
        )
        self.assertEqual(
            record["y_area_weighted_profiles"]["rho"]["nominal_slot_time"],
            500.0,
        )
        with self.assertRaisesRegex(
            spatial.AnalysisError, "embedded mesh time differs"
        ):
            spatial.reduce_t500_spatial_snapshot(
                _snapshot_datasets(embedded_time=observed),
                nominal_slot_time=500.0,
                observed_committed_time=observed,
                x_ideal_c_over_omega_pi=500.0,
            )

    def test_malformed_shapes_and_nonfinite_values_fail_closed(self) -> None:
        with self.assertRaisesRegex(spatial.AnalysisError, "shape disagrees"):
            _raster("rho", [[1.0, 2.0, 3.0]])
        with self.assertRaisesRegex(spatial.AnalysisError, "must be finite"):
            _raster("rho", [[1.0, np.nan], [2.0, 3.0]])
        with self.assertRaisesRegex(spatial.AnalysisError, "integer array"):
            _raster("rho", [[1.0, 2.0], [3.0, 4.0]], source_levels=np.ones((2, 2)))
        with self.assertRaisesRegex(spatial.AnalysisError, "must be finite"):
            spatial.ideal_shock_search_window(np.inf)
        malformed_metadata = replace(
            _dataset("dens", np.ones((1, 2, 2))), root_grid_shape=(2, 2)
        )
        with self.assertRaisesRegex(spatial.AnalysisError, "exactly three axes"):
            spatial.compose_xy_quantity(malformed_metadata, "rho")
        with self.assertRaisesRegex(spatial.AnalysisError, "must be finite"):
            spatial.DetectedFrontRecord(
                nominal_slot_time=500.0,
                observed_committed_time=500.0,
                x_ideal_c_over_omega_pi=500.0,
                search_window_c_over_omega_pi=(-700.0, 1700.0),
                front_index=1,
                x_front_c_over_omega_pi=500.0,
                density_gradient=np.nan,
                offset_from_x_ideal_c_over_omega_pi=0.0,
                max_absolute_offset_c_over_omega_pi=600.0,
            )


if __name__ == "__main__":
    unittest.main()
