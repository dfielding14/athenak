#!/usr/bin/env python3
"""Focused synthetic tests for Q011 production-science successor v1."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import copy
import unittest
from unittest.mock import patch

import numpy as np

try:
    from tst.publication import analyze_q011_section54_outputs as output_primitives
    from tst.publication import q011_section54_production_science_successor_v1 as science
except ModuleNotFoundError:
    import analyze_q011_section54_outputs as output_primitives
    import q011_section54_production_science_successor_v1 as science


REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACT = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_production_science_diagnostic_successor_v1_2026-06-06.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _velocity_from_chi(values: object) -> np.ndarray:
    chi = np.asarray(values, dtype=np.float64)
    momentum = science.UPSTREAM_SPEED_U0 * np.sqrt(chi)
    velocity = momentum / np.sqrt(
        1.0 + (momentum / science.PARTICLE_LIGHT_SPEED) ** 2
    )
    return np.column_stack(
        (velocity, np.zeros_like(velocity), np.zeros_like(velocity))
    ).astype(np.float32).astype(np.float64)


def _particle_payload(maximum_chi: float, front_x: float) -> dict[str, object]:
    active_count = science.PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES
    active_chi = np.linspace(1.0, maximum_chi, active_count)
    active_x = np.where(np.arange(active_count) % 2 == 0, front_x - 500.0, front_x + 500.0)
    return {
        "points": np.column_stack(
            (
                np.concatenate((active_x, [front_x + 1500.0, front_x - 500.0])),
                np.zeros(active_count + 2),
                np.zeros(active_count + 2),
            )
        ),
        "cr_source": np.concatenate(
            (np.ones(active_count, dtype=np.int64), np.asarray([0, 1], dtype=np.int64))
        ),
        "birth_time": np.concatenate(
            (np.full(active_count, 45.0), np.asarray([50.0, 44.0]))
        ),
        "velocity": _velocity_from_chi(np.concatenate((active_chi, [64.0, 100.0]))),
        "macro_weight": np.ones(active_count + 2),
    }


def _mhd_dataset(time: float) -> output_primitives.AthenaBinaryDataset:
    nx1 = 72
    nx2 = 2
    dx1 = 200.0
    x1_centers = (np.arange(nx1) + 0.5) * dx1
    ideal = science.frozen_model.x_ideal(time)
    front_index = int(round((ideal - 100.0 - 0.5 * dx1) / dx1))
    front_x = float(x1_centers[front_index])

    density = np.ones(nx1)
    density[:front_index] = 4.0
    density[front_index] = 3.8
    density = np.tile(density, (1, nx2, 1))

    upstream = x1_centers > front_x + 120.0
    bcc1 = np.ones((1, nx2, nx1))
    bcc1[0, 0, upstream] = 2.0
    bcc1[0, 1, upstream] = 4.0
    zeros = np.zeros_like(bcc1)
    velocity_x = np.zeros_like(bcc1)
    velocity_x[:, :, upstream] = -30.0
    fields = {
        "dens": density,
        "velx": velocity_x,
        "vely": zeros,
        "velz": zeros,
        "eint": np.full_like(bcc1, 1.5),
        "bcc1": bcc1,
        "bcc2": zeros,
        "bcc3": zeros,
    }
    block = output_primitives.AthenaBinaryBlock(
        index_bounds=(0, nx1 - 1, 0, nx2 - 1, 0, 0),
        logical_location=(0, 0, 0),
        level=0,
        geometry=(0.0, nx1 * dx1, 0.0, 2.0, 0.0, 1.0),
        fields=fields,
    )
    return output_primitives.AthenaBinaryDataset(
        source=f"synthetic-mhd-w-bcc-t{time:g}",
        time=time,
        cycle=int(time),
        location_size=8,
        variable_size=8,
        variable_names=(
            "bcc3",
            "dens",
            "velz",
            "bcc1",
            "eint",
            "velx",
            "bcc2",
            "vely",
        ),
        input_parameters={},
        root_grid_shape=(nx1, nx2, 1),
        meshblock_shape=(nx1, nx2, 1),
        nghost=0,
        domain_bounds=(0.0, nx1 * dx1, 0.0, 2.0, 0.0, 1.0),
        blocks=(block,),
    )


def _current_datasets(
    mhd_dataset: output_primitives.AthenaBinaryDataset,
) -> dict[str, output_primitives.AthenaBinaryDataset]:
    shape = mhd_dataset.blocks[0].fields["dens"].shape
    values = {
        "prtcl_rho": np.full(shape, 0.1),
        "prtcl_jx": np.full(shape, 1.0),
        "prtcl_jy": np.full(shape, 2.0),
        "prtcl_jz": np.full(shape, 3.0),
        "mhd_j2": np.full(shape, 14.0),
    }
    datasets = {}
    reference_block = mhd_dataset.blocks[0]
    for product, field in science.CURRENT_PRODUCT_FIELDS.items():
        block = output_primitives.AthenaBinaryBlock(
            index_bounds=reference_block.index_bounds,
            logical_location=reference_block.logical_location,
            level=reference_block.level,
            geometry=reference_block.geometry,
            fields={field: values[product]},
        )
        datasets[product] = output_primitives.AthenaBinaryDataset(
            **{
                **mhd_dataset.__dict__,
                "source": f"synthetic-{product}",
                "variable_names": (field,),
                "blocks": (block,),
            }
        )
    return datasets


def _late_spectrum_payload() -> dict[str, object]:
    edges = np.asarray(science.particle_primitives.CHI_BIN_EDGES)
    centers = np.sqrt(edges[:-1] * edges[1:])
    selected = (centers >= 20.0) & (centers <= 160.0)
    selected_centers = centers[selected]
    selected_widths = np.diff(edges)[selected]
    repetitions = int(
        np.ceil(
            science.PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES
            / selected_centers.size
        )
    )
    chi = np.repeat(selected_centers, repetitions)
    weights = np.repeat(
        7.0 * selected_centers**-1.5 * selected_widths / repetitions,
        repetitions,
    )
    count = chi.size
    return {
        "points": np.column_stack(
            (np.full(count, 100.0), np.zeros(count), np.zeros(count))
        ),
        "cr_source": np.ones(count, dtype=np.int64),
        "birth_time": np.full(count, 45.0),
        "velocity": _velocity_from_chi(chi),
        "macro_weight": weights,
    }


def _production_snapshot(
    time: float,
    maximum_chi: float,
    particle_payload: dict[str, object] | None = None,
) -> dict[str, object]:
    dataset = _mhd_dataset(time)
    mhd = science.reduce_mhd_snapshot(
        dataset,
        nominal_slot_time=time,
        observed_committed_time=time,
    )
    particles = (
        _particle_payload(
            maximum_chi, mhd["detected_front"]["x_front_c_over_omega_pi"]
        )
        if particle_payload is None
        else particle_payload
    )
    return science.reduce_production_science_snapshot(
        dataset,
        _current_datasets(dataset),
        nominal_slot_time=time,
        observed_committed_time=time,
        **particles,
    )


class Q011ProductionScienceSuccessorV1Tests(unittest.TestCase):
    def test_successor_deck_preserves_physics_and_has_no_authority(self) -> None:
        record = science.validate_successor_deck()
        self.assertTrue(record["historical_physical_baseline_preserved"])
        self.assertEqual(
            record["output_ids"],
            [
                "mhd_w_bcc",
                "prtcl_rho",
                "prtcl_jx",
                "prtcl_jy",
                "prtcl_jz",
                "mhd_j2",
                "prtcl_all",
            ],
        )
        self.assertEqual(record["history_cadence_omega0_inverse"], 10.0)
        self.assertFalse(record["launch_authorized"])
        self.assertFalse(record["policy_mutation_authorized"])
        self.assertFalse(record["claim_closure_authorized"])

    def test_particle_maximum_energy_and_cr_energy_use_active_shock_population(self) -> None:
        payload = _particle_payload(9.0, 1900.0)
        record = science.reduce_particle_energy_snapshot(
            snapshot_time=200.0, **payload
        )
        self.assertEqual(
            record["particle_census"]["all_particle_count"],
            science.PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES + 2,
        )
        self.assertEqual(
            record["particle_census"]["active_particle_count"],
            science.PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES,
        )
        self.assertEqual(
            record["particle_census"]["active_positive_weight_particle_count"],
            science.PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES,
        )
        self.assertAlmostEqual(
            record["maximum_particle_energy"]["maximum_chi"], 9.0, places=5
        )
        self.assertGreater(record["active_cr_kinetic_energy"], 0.0)
        self.assertAlmostEqual(
            record["particle_census"]["active_cr_macro_mass"],
            science.PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES
            * science.PARTICLE_MACRO_MASS,
        )
        active_velocity = np.asarray(payload["velocity"])[
            : science.PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES
        ]
        speed_squared = np.sum(active_velocity * active_velocity, axis=1)
        momentum_squared = speed_squared / (
            1.0 - speed_squared / science.PARTICLE_LIGHT_SPEED**2
        )
        specific_energy = momentum_squared / (
            np.sqrt(1.0 + momentum_squared / science.PARTICLE_LIGHT_SPEED**2) + 1.0
        )
        expected_cr_energy = science.PARTICLE_MACRO_MASS * np.sum(specific_energy)
        self.assertAlmostEqual(record["active_cr_kinetic_energy"], expected_cr_energy)
        self.assertEqual(
            record["maximum_particle_energy"]["tail_quantile_sample_requirement"][
                "minimum_positive_weight_particle_count"
            ],
            science.PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES,
        )
        uncertainty = record["float32_pvtk_velocity_uncertainty"]
        max_chi_interval = uncertainty["maximum_particle_energy_intervals"][
            "maximum_chi"
        ]
        self.assertLessEqual(max_chi_interval["lower"], 9.0)
        self.assertGreaterEqual(max_chi_interval["upper"], 9.0)
        self.assertTrue(uncertainty["all_nominal_values_contained"])
        self.assertEqual(
            set(
                uncertainty["macro_weighted_tail_quantile_intervals"][
                    "specific_kinetic_energy"
                ]
            ),
            {"q990", "q999"},
        )

    def test_full_state_reducer_closes_shock_compression_and_local_b_statistics(
        self,
    ) -> None:
        record = science.reduce_mhd_snapshot(
            _mhd_dataset(200.0),
            nominal_slot_time=200.0,
            observed_committed_time=200.0,
        )
        self.assertEqual(
            record["detected_front"]["detector"],
            "unique_strongest_negative_density_gradient",
        )
        self.assertFalse(record["detected_front"]["qualifying_use_authorized"])
        self.assertIn(
            "requires_preregistration_supersession",
            record["detected_front"]["preregistration_status"],
        )
        self.assertLess(record["detected_front"]["density_gradient"], 0.0)
        self.assertAlmostEqual(
            record["shock_state"]["compression_ratio_downstream_over_upstream"],
            4.0,
        )
        amplification = record["upstream_magnetic_amplification"][
            "qualifying_ideal_surface_relative"
        ]
        self.assertIn("preregistered_qualifying_metric", amplification["qualification_role"])
        self.assertAlmostEqual(amplification["area_weighted_mean_abs_b_over_b0"], 3.0)
        self.assertAlmostEqual(amplification["local_max_abs_b_over_b0"], 4.0)
        self.assertEqual(
            amplification["area_weighted_quantiles_abs_b_over_b0"],
            {"q50": 2.0, "q90": 4.0, "q99": 4.0},
        )
        self.assertAlmostEqual(
            amplification["area_fractions"]["area_fraction_abs_b_over_b0_ge_3"],
            0.5,
        )
        supplemental = record["upstream_magnetic_amplification"][
            "supplemental_detected_front_relative"
        ]
        self.assertIn("supplemental_only", supplemental["qualification_role"])
        self.assertFalse(amplification["acceptance_gate"]["evaluated"])
        full_state = record["full_state"]
        self.assertEqual(
            full_state["raw_primitive_field_inventory"],
            list(science.MHD_PRIMITIVE_FIELDS),
        )
        self.assertEqual(
            set(full_state["y_area_weighted_profiles"]),
            {*science.MHD_PRIMITIVE_FIELDS, "pressure", "bmag"},
        )
        self.assertTrue(
            record["mhd_energy_components"]["closure_successor_required"]
        )

    def test_matched_current_reducer_reports_lab_gas_and_mhd_current_squared(self) -> None:
        dataset = _mhd_dataset(200.0)
        mhd = science.reduce_mhd_snapshot(dataset)
        record = science.reduce_cr_current_snapshot(
            dataset,
            _current_datasets(dataset),
            nominal_slot_time=200.0,
            observed_committed_time=200.0,
            detected_front_x1_c_over_omega_pi=mhd["detected_front"][
                "x_front_c_over_omega_pi"
            ],
        )
        full = record["statistics_by_region"]["full_domain"]
        self.assertEqual(full["j_cr_lab"]["area_weighted_mean_vector"], [1.0, 2.0, 3.0])
        ideal_upstream = record["statistics_by_region"][
            "qualifying_ideal_surface_relative_upstream"
        ]
        self.assertEqual(
            ideal_upstream["j_cr_gas"]["area_weighted_mean_vector"],
            [4.0, 2.0, 3.0],
        )
        self.assertEqual(
            ideal_upstream["mhd_current_squared"]["area_weighted_mean"], 14.0
        )
        self.assertTrue(
            record["mhd_current_squared_interpretation"]["not_particle_current_squared"]
        )
        self.assertEqual(
            record["source_representation"]["prtcl_j_components"],
            "deposited_J_CR_over_c",
        )
        self.assertEqual(
            record["source_representation"]["prtcl_rho"],
            "deposited_rho_CR_over_c",
        )
        self.assertEqual(
            record["frame_transform"]["formula"],
            "(J_CR/c)_gas = (J_CR/c)_lab - (rho_CR/c) * u_gas",
        )
        self.assertEqual(
            record["matched_metadata_contract"]["normalized_output_state_fields"],
            ["file_number", "last_time"],
        )
        self.assertEqual(
            set(record["y_area_weighted_profiles"]),
            {
                "x1_centers_c_over_omega_pi",
                "prtcl_rho",
                "j_cr_lab_x",
                "j_cr_lab_y",
                "j_cr_lab_z",
                "j_cr_lab_magnitude",
                "j_cr_gas_x",
                "j_cr_gas_y",
                "j_cr_gas_z",
                "j_cr_gas_magnitude",
                "mhd_current_squared",
            },
        )

    def test_matched_current_reducer_normalizes_only_output_state_counters(self) -> None:
        base = _mhd_dataset(200.0)
        parameters = {
            "particles": {"pic_cr_light_speed": "10000.0"},
            **{
                f"output{index}": {
                    "file_type": "bin",
                    "variable": f"field{index}",
                    "file_number": "2",
                    "last_time": "200",
                }
                for index in range(1, 7)
            },
        }
        dataset = output_primitives.AthenaBinaryDataset(
            **{**base.__dict__, "input_parameters": parameters}
        )
        currents = _current_datasets(dataset)
        for product_index, product in enumerate(science.CURRENT_PRODUCT_FIELDS, 1):
            product_parameters = copy.deepcopy(parameters)
            for output_index in range(1, 7):
                product_parameters[f"output{output_index}"]["file_number"] = str(
                    2 if output_index <= product_index else 1
                )
                product_parameters[f"output{output_index}"]["last_time"] = (
                    "200" if output_index <= product_index else "100"
                )
            currents[product] = output_primitives.AthenaBinaryDataset(
                **{
                    **currents[product].__dict__,
                    "input_parameters": product_parameters,
                }
            )
        mhd = science.reduce_mhd_snapshot(dataset)
        science.reduce_cr_current_snapshot(
            dataset,
            currents,
            nominal_slot_time=200.0,
            observed_committed_time=200.0,
            detected_front_x1_c_over_omega_pi=mhd["detected_front"][
                "x_front_c_over_omega_pi"
            ],
        )

        drifted = dict(currents)
        bad_parameters = copy.deepcopy(drifted["prtcl_jx"].input_parameters)
        bad_parameters["particles"]["pic_cr_light_speed"] = "9999.0"
        drifted["prtcl_jx"] = output_primitives.AthenaBinaryDataset(
            **{
                **drifted["prtcl_jx"].__dict__,
                "input_parameters": bad_parameters,
            }
        )
        with self.assertRaisesRegex(science.ProductionScienceError, "metadata disagrees"):
            science.reduce_cr_current_snapshot(
                dataset,
                drifted,
                nominal_slot_time=200.0,
                observed_committed_time=200.0,
                detected_front_x1_c_over_omega_pi=mhd["detected_front"][
                    "x_front_c_over_omega_pi"
                ],
            )

        malformed = dict(currents)
        bad_counter_parameters = copy.deepcopy(malformed["prtcl_jx"].input_parameters)
        bad_counter_parameters["output1"]["file_number"] = "not-an-integer"
        malformed["prtcl_jx"] = output_primitives.AthenaBinaryDataset(
            **{
                **malformed["prtcl_jx"].__dict__,
                "input_parameters": bad_counter_parameters,
            }
        )
        with self.assertRaisesRegex(science.ProductionScienceError, "file_number"):
            science.reduce_cr_current_snapshot(
                dataset,
                malformed,
                nominal_slot_time=200.0,
                observed_committed_time=200.0,
                detected_front_x1_c_over_omega_pi=mhd["detected_front"][
                    "x_front_c_over_omega_pi"
                ],
            )

    def test_preregistered_spectrum_reuses_particle_reducer_and_late_fit(self) -> None:
        record = science.reduce_particle_energy_snapshot(
            snapshot_time=1200.0,
            nominal_slot_time=1200.0,
            **_late_spectrum_payload(),
        )
        spectrum = record["preregistered_downstream_spectrum"]
        self.assertEqual(
            spectrum["reuse_contract"],
            "q011_section54_particles.reduce_particle_snapshot",
        )
        self.assertEqual(
            spectrum["minimum_admitted_particle_count"],
            science.DOWNSTREAM_SPECTRUM_MINIMUM_PARTICLE_SAMPLES,
        )
        self.assertEqual(
            spectrum["reduction"]["particle_filter"]["selection"]["spatial_rule"],
            "x1 < x_ideal(t)",
        )
        self.assertGreaterEqual(
            spectrum["reduction"]["late_slope"]["positive_fit_bin_count"],
            science.particle_primitives.LATE_SLOPE_MINIMUM_POSITIVE_BINS,
        )
        self.assertAlmostEqual(
            spectrum["reduction"]["late_slope"]["slope"], -1.5, places=12
        )

    def test_matched_snapshot_reports_full_and_shock_centered_energy_partitions(
        self,
    ) -> None:
        record = _production_snapshot(200.0, 9.0)
        for scope in ("full_domain", "shock_centered_window"):
            partition = record["energy_partition"][scope]
            self.assertEqual(
                partition["interpretation"],
                "instantaneous_wall_frame_partition_not_a_conserved_energy_budget",
            )
            self.assertFalse(
                partition["closure_status"]["conserved_energy_closure_evaluated"]
            )
            self.assertGreater(partition["gas_energy"], 0.0)
            self.assertGreater(partition["magnetic_energy"], 0.0)
            self.assertGreater(partition["cr_kinetic_energy"], 0.0)
            self.assertAlmostEqual(sum(partition["fractions"].values()), 1.0)
        self.assertEqual(
            record["energy_partition"]["shock_centered_window"][
                "active_cr_particle_count"
            ],
            science.PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES,
        )
        self.assertTrue(record["energy_partition"]["closure_successor_required"])
        self.assertIn("cr_current", record)

    def test_history_reports_acceleration_rate_shock_speed_and_time_series(self) -> None:
        history = science.reduce_production_science_history(
            [
                _production_snapshot(300.0, 16.0),
                _production_snapshot(200.0, 9.0),
                _production_snapshot(250.0, 12.5),
            ]
        )
        self.assertEqual(
            history["observed_committed_times_omega0_inverse"], [200.0, 250.0, 300.0]
        )
        self.assertAlmostEqual(
            history["particle_acceleration"]["maximum_chi_linear_fit"]["slope"],
            0.07,
        )
        self.assertAlmostEqual(
            history["shock_kinematics"][
                "measured_wall_frame_shock_speed_linear_fit"
            ]["slope"],
            10.0,
        )
        self.assertEqual(
            history["fit_requirements"]["linear_fit_minimum_snapshot_count"], 3
        )
        self.assertEqual(len(history["interval_rates"]), 2)
        self.assertEqual(len(history["energy_partition_time_series"]), 3)
        self.assertEqual(len(history["upstream_magnetic_amplification_time_series"]), 3)
        self.assertEqual(len(history["cr_current_statistics_time_series"]), 3)
        self.assertEqual(
            len(
                history["particle_acceleration"][
                    "float32_pvtk_velocity_uncertainty_time_series"
                ]
            ),
            3,
        )
        self.assertEqual(
            history["particle_acceleration"][
                "macro_weighted_chi_tail_quantile_history"
            ]["minimum_positive_weight_particle_count_per_snapshot"],
            science.PARTICLE_TAIL_MINIMUM_POSITIVE_WEIGHT_SAMPLES,
        )
        self.assertEqual(history["preregistered_downstream_spectrum_time_series"], [])
        self.assertEqual(json.loads(science.canonical_record_bytes(history)), history)

    def test_history_propagates_preregistered_t500_t1200_spectra_and_late_slope(
        self,
    ) -> None:
        history = science.reduce_production_science_history(
            [
                _production_snapshot(1200.0, 160.0, _late_spectrum_payload()),
                _production_snapshot(500.0, 160.0, _late_spectrum_payload()),
                _production_snapshot(700.0, 160.0),
            ]
        )
        spectra = history["preregistered_downstream_spectrum_time_series"]
        self.assertEqual([row["nominal_slot_time"] for row in spectra], [500.0, 1200.0])
        self.assertNotIn("late_slope", spectra[0]["reduction"])
        self.assertGreaterEqual(
            spectra[1]["reduction"]["late_slope"]["positive_fit_bin_count"],
            science.particle_primitives.LATE_SLOPE_MINIMUM_POSITIVE_BINS,
        )
        self.assertEqual(
            history["spectrum_requirements"][
                "minimum_admitted_particle_count_per_spectrum"
            ],
            science.DOWNSTREAM_SPECTRUM_MINIMUM_PARTICLE_SAMPLES,
        )

    def test_fail_closed_on_inventory_nonfinite_and_incomplete_windows(self) -> None:
        dataset = _mhd_dataset(200.0)
        missing = output_primitives.AthenaBinaryDataset(
            **{**dataset.__dict__, "variable_names": dataset.variable_names[:-1]}
        )
        with self.assertRaisesRegex(science.ProductionScienceError, "inventory drifted"):
            science.compose_full_mhd_state(missing)

        bad_fields = dict(dataset.blocks[0].fields)
        bad_fields["dens"] = np.array(bad_fields["dens"], copy=True)
        bad_fields["dens"][0, 0, 0] = np.nan
        bad_block = output_primitives.AthenaBinaryBlock(
            **{**dataset.blocks[0].__dict__, "fields": bad_fields}
        )
        bad = output_primitives.AthenaBinaryDataset(
            **{**dataset.__dict__, "blocks": (bad_block,)}
        )
        with self.assertRaisesRegex(science.ProductionScienceError, "not finite"):
            science.compose_full_mhd_state(bad)

        with self.assertRaisesRegex(science.ProductionScienceError, "escaped"):
            science.reduce_mhd_snapshot(
                _mhd_dataset(100.0),
                nominal_slot_time=100.0,
                observed_committed_time=100.0,
            )
        with self.assertRaisesRegex(science.ProductionScienceError, "begin at nominal t=200"):
            payload = _particle_payload(9.0, 900.0)
            science.reduce_production_science_snapshot(
                _mhd_dataset(100.0),
                _current_datasets(_mhd_dataset(100.0)),
                nominal_slot_time=100.0,
                observed_committed_time=100.0,
                **payload,
            )
        payload = _particle_payload(9.0, 1900.0)
        payload["macro_weight"] = np.zeros_like(payload["macro_weight"])
        with self.assertRaisesRegex(science.ProductionScienceError, "positive-weight"):
            science.reduce_particle_energy_snapshot(snapshot_time=200.0, **payload)

        too_few = _late_spectrum_payload()
        for key in too_few:
            too_few[key] = np.asarray(too_few[key])[:10]
        with self.assertRaisesRegex(
            science.ProductionScienceError, "tail quantiles require at least"
        ):
            science.reduce_particle_energy_snapshot(snapshot_time=200.0, **too_few)

        not_decoded_float32 = _particle_payload(9.0, 1900.0)
        not_decoded_float32["velocity"] = np.array(
            not_decoded_float32["velocity"], copy=True
        )
        not_decoded_float32["velocity"][0, 0] += 1.0e-10
        with self.assertRaisesRegex(
            science.ProductionScienceError, "exactly representable decoded float32"
        ):
            science.reduce_particle_energy_snapshot(
                snapshot_time=200.0, **not_decoded_float32
            )

        with self.assertRaisesRegex(
            science.ProductionScienceError, "fit requires at least 3 snapshots"
        ):
            science.reduce_production_science_history(
                [_production_snapshot(200.0, 9.0), _production_snapshot(300.0, 16.0)]
            )

        bad_uncertainty = _production_snapshot(250.0, 12.5)
        interval = bad_uncertainty["particles"]["float32_pvtk_velocity_uncertainty"][
            "maximum_particle_energy_intervals"
        ]["maximum_chi"]
        interval["lower"] = interval["upper"] + 1.0
        with self.assertRaisesRegex(
            science.ProductionScienceError, "nominal value escaped"
        ):
            science.reduce_production_science_history(
                [
                    _production_snapshot(200.0, 9.0),
                    bad_uncertainty,
                    _production_snapshot(300.0, 16.0),
                ]
            )

        current = _current_datasets(dataset)
        current["prtcl_jy"] = output_primitives.AthenaBinaryDataset(
            **{**current["prtcl_jy"].__dict__, "cycle": dataset.cycle + 1}
        )
        mhd = science.reduce_mhd_snapshot(dataset)
        with self.assertRaisesRegex(
            science.ProductionScienceError, "metadata disagrees"
        ):
            science.reduce_cr_current_snapshot(
                dataset,
                current,
                nominal_slot_time=200.0,
                observed_committed_time=200.0,
                detected_front_x1_c_over_omega_pi=mhd["detected_front"][
                    "x_front_c_over_omega_pi"
                ],
            )

    def test_public_contract_wraps_underlying_analysis_and_parser_exceptions(self) -> None:
        public_functions = (
            "validate_successor_deck",
            "compose_full_mhd_state",
            "reduce_mhd_snapshot",
            "reduce_cr_current_snapshot",
            "reduce_particle_energy_snapshot",
            "reduce_production_science_snapshot",
            "reduce_production_science_history",
            "canonical_record_bytes",
        )
        for name in public_functions:
            with self.subTest(public_function=name):
                self.assertTrue(hasattr(getattr(science, name), "__wrapped__"))

        nonfinite_error = output_primitives.AnalysisError(
            "synthetic nonfinite-data AnalysisError"
        )
        with patch.object(
            science.output_primitives,
            "compose_leaf_field",
            side_effect=nonfinite_error,
        ):
            with self.assertRaisesRegex(
                science.ProductionScienceError,
                "full MHD state composition failed: synthetic nonfinite-data AnalysisError",
            ) as caught:
                science.compose_full_mhd_state(_mhd_dataset(200.0))
        self.assertIs(caught.exception.__cause__, nonfinite_error)

        shock_error = output_primitives.AnalysisError(
            "synthetic shock arithmetic AnalysisError"
        )
        with patch.object(
            science.output_primitives,
            "detect_shock_front",
            side_effect=shock_error,
        ):
            with self.assertRaisesRegex(
                science.ProductionScienceError,
                "MHD snapshot reduction failed: synthetic shock arithmetic AnalysisError",
            ) as caught:
                science.reduce_mhd_snapshot(_mhd_dataset(200.0))
        self.assertIs(caught.exception.__cause__, shock_error)

        model_error = science.frozen_model.ModelContractError(
            "synthetic model parser failure"
        )
        with patch.object(
            science.frozen_model,
            "parse_deck_contract",
            side_effect=model_error,
        ):
            with self.assertRaisesRegex(
                science.ProductionScienceError,
                "successor deck validation failed: synthetic model parser failure",
            ) as caught:
                science.validate_successor_deck()
        self.assertIs(caught.exception.__cause__, model_error)

        particle_error = science.particle_primitives.ParticleReducerError(
            "synthetic downstream spectrum failure"
        )
        with patch.object(
            science.particle_primitives,
            "reduce_particle_snapshot",
            side_effect=particle_error,
        ):
            with self.assertRaisesRegex(
                science.ProductionScienceError,
                "particle-energy snapshot reduction failed: synthetic downstream spectrum failure",
            ) as caught:
                science.reduce_particle_energy_snapshot(
                    snapshot_time=1200.0,
                    nominal_slot_time=1200.0,
                    **_late_spectrum_payload(),
                )
        self.assertIs(caught.exception.__cause__, particle_error)

        missing = REPO_ROOT / "inputs/publication/absent-q011-successor.athinput"
        with self.assertRaisesRegex(
            science.ProductionScienceError, "successor deck validation failed"
        ) as caught:
            science.validate_successor_deck(missing)
        self.assertIsInstance(caught.exception.__cause__, FileNotFoundError)

        with self.assertRaisesRegex(
            science.ProductionScienceError,
            "production-science history reduction failed",
        ) as caught:
            science.reduce_production_science_history([{}, {}])
        self.assertIsInstance(caught.exception.__cause__, KeyError)

        cyclic: dict[str, object] = {}
        cyclic["self"] = cyclic
        with self.assertRaisesRegex(
            science.ProductionScienceError,
            "canonical diagnostic serialization failed",
        ) as caught:
            science.canonical_record_bytes(cyclic)
        self.assertIsInstance(caught.exception.__cause__, ValueError)

    def test_versioned_contract_binds_successor_sources_and_refuses_authority(self) -> None:
        contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        self.assertNotIn("REPLACE_", CONTRACT.read_text(encoding="utf-8"))
        self.assertEqual(contract["schema_version"], 1)
        self.assertEqual(contract["successor_id"], science.SUCCESSOR_ID)
        self.assertEqual(
            contract["qualification_effect"],
            "source_local_diagnostic_contract_only_no_launch_no_policy_authorization_no_claim_closure",
        )
        self.assertEqual(
            contract["authorization"],
            {
                "launch_authorized": False,
                "policy_mutation_authorized": False,
                "claim_closure_authorized": False,
                "qualifying_output_inspection_authorized": False,
            },
        )
        for binding in contract["source_bindings"].values():
            path = REPO_ROOT / binding["path"]
            self.assertTrue(path.is_file())
            self.assertEqual(binding["sha256"], _sha256(path))
        predecessor = contract["predecessor_contract"]
        self.assertEqual(
            predecessor["sha256"], _sha256(REPO_ROOT / predecessor["path"])
        )
        self.assertEqual(
            contract["output_contract"]["mhd_w_bcc_required_fields"],
            list(science.MHD_PRIMITIVE_FIELDS),
        )
        self.assertEqual(
            contract["diagnostic_contract"]["shock_front"]["gradient_sign"],
            "negative",
        )
        self.assertTrue(
            contract["diagnostic_contract"]["shock_front"][
                "preregistration_supersession_required"
            ]
        )
        self.assertEqual(
            contract["output_contract"]["deposited_cr_moment_products"],
            ["prtcl_rho", "prtcl_jx", "prtcl_jy", "prtcl_jz"],
        )
        self.assertEqual(
            contract["output_contract"]["mhd_j2_interpretation"],
            "MHD current magnitude squared |curl B|^2; not particle current squared",
        )
        self.assertEqual(
            contract["diagnostic_contract"]["upstream_magnetic_amplification"][
                "qualifying_metric"
            ]["window_reference"],
            "ideal_surface",
        )
        self.assertEqual(
            contract["diagnostic_contract"]["cr_current"]["gas_frame_formula"],
            "(J_CR/c)_gas = (J_CR/c)_lab - (rho_CR/c) * u_gas",
        )
        self.assertEqual(
            contract["diagnostic_contract"]["matched_snapshot"][
                "first_eligible_nominal_slot_omega0_inverse"
            ],
            science.MATCHED_SNAPSHOT_FIRST_ELIGIBLE_NOMINAL_TIME,
        )
        self.assertIn(
            "deposited J_CR/c",
            contract["diagnostic_contract"]["cr_current"][
                "raw_prtcl_j_representation"
            ],
        )
        self.assertEqual(
            contract["diagnostic_contract"]["downstream_spectrum"][
                "minimum_admitted_particle_count_per_spectrum"
            ],
            science.DOWNSTREAM_SPECTRUM_MINIMUM_PARTICLE_SAMPLES,
        )
        self.assertTrue(
            contract["diagnostic_contract"]["energy_partition"][
                "separate_conserved_energy_closure_successor_required"
            ]
        )
        self.assertEqual(
            contract["diagnostic_contract"]["particle_acceleration"][
                "linear_fit_minimum_snapshot_count"
            ],
            science.HISTORY_LINEAR_FIT_MINIMUM_SNAPSHOTS,
        )
        self.assertEqual(
            contract["diagnostic_contract"]["particle_acceleration"][
                "float32_pvtk_velocity_uncertainty"
            ]["bounded_observables"],
            [
                "maximum_chi",
                "maximum_specific_kinetic_energy",
                "macro_weighted_chi_q990_q999",
                "macro_weighted_specific_kinetic_energy_q990_q999",
            ],
        )


if __name__ == "__main__":
    unittest.main()
