#!/usr/bin/env python3
"""Focused adversarial tests for the Q011 physical-applicability successor."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

try:
    from tst.publication import analyze_q011_section54_outputs as output_primitives
    from tst.publication import q011_section54_physical_applicability_successor_v1 as app
    from tst.publication import q011_section54_production_science_successor_v1 as science
except ModuleNotFoundError:
    import analyze_q011_section54_outputs as output_primitives
    import q011_section54_physical_applicability_successor_v1 as app
    import q011_section54_production_science_successor_v1 as science


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_physical_applicability_successor_v1_2026-06-06.json"
)
DESIGN = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_physical_applicability_successor_v1_2026-06-06.md"
)
NX1 = 480
NX2 = 32
DX1 = 20.0
DX2 = 10.0
TIME = 500.0


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _mhd_dataset(
    *,
    no_fluctuation: bool = False,
    low_density_cell: bool = False,
    high_k_fluctuation: bool = False,
) -> output_primitives.AthenaBinaryDataset:
    x = (np.arange(NX1) + 0.5) * DX1
    y = (np.arange(NX2) + 0.5) * DX2
    front = 5000.0
    front_index = int(np.searchsorted(x, front))
    density_x = np.full(NX1, 0.04 if high_k_fluctuation else 1.0)
    density_x[:front_index] = 4.0
    density_x[front_index] = 3.8
    upstream = np.arange(NX1) > front_index
    density = np.broadcast_to(density_x[None, :], (NX2, NX1)).copy()
    if low_density_cell:
        density[0, -1] = 1.0e-4

    velx = np.broadcast_to(np.where(upstream, -30.0, 0.0)[None, :], (NX2, NX1)).copy()
    zeros = np.zeros((NX2, NX1), dtype=np.float64)
    bcc1 = np.ones((NX2, NX1), dtype=np.float64)
    if no_fluctuation:
        bcc2 = zeros.copy()
        bcc3 = zeros.copy()
    elif high_k_fluctuation:
        bcc2 = 0.05 * np.sin(2.0 * np.pi * x[None, :] / 40.0)
        bcc2 = np.broadcast_to(bcc2, (NX2, NX1)).copy()
        bcc3 = zeros.copy()
    else:
        bcc2 = 0.05 * np.sin(2.0 * np.pi * x[None, :] / 400.0)
        bcc2 = np.broadcast_to(bcc2, (NX2, NX1)).copy()
        bcc3 = 0.05 * np.cos(2.0 * np.pi * y[:, None] / 160.0)
        bcc3 = np.broadcast_to(bcc3, (NX2, NX1)).copy()
    fields = {
        "dens": density[None, :, :],
        "velx": velx[None, :, :],
        "vely": zeros[None, :, :],
        "velz": zeros[None, :, :],
        "eint": np.full((1, NX2, NX1), 1.5),
        "bcc1": bcc1[None, :, :],
        "bcc2": bcc2[None, :, :],
        "bcc3": bcc3[None, :, :],
    }
    block = output_primitives.AthenaBinaryBlock(
        index_bounds=(0, NX1 - 1, 0, NX2 - 1, 0, 0),
        logical_location=(0, 0, 0),
        level=0,
        geometry=(0.0, NX1 * DX1, 0.0, NX2 * DX2, 0.0, 1.0),
        fields=fields,
    )
    return output_primitives.AthenaBinaryDataset(
        source="synthetic-q011-applicability-mhd",
        time=TIME,
        cycle=500,
        location_size=8,
        variable_size=8,
        variable_names=tuple(reversed(science.MHD_PRIMITIVE_FIELDS)),
        input_parameters={},
        root_grid_shape=(NX1, NX2, 1),
        meshblock_shape=(NX1, NX2, 1),
        nghost=0,
        domain_bounds=(0.0, NX1 * DX1, 0.0, NX2 * DX2, 0.0, 1.0),
        blocks=(block,),
    )


def _current_datasets(
    mhd: output_primitives.AthenaBinaryDataset,
    *,
    rho_q: float = 0.005,
    gas_frame_current: float = 0.02,
    rare_lambda_failure: bool = False,
) -> dict[str, output_primitives.AthenaBinaryDataset]:
    velx = mhd.blocks[0].fields["velx"]
    shape = velx.shape
    jx = rho_q * velx + gas_frame_current
    if rare_lambda_failure:
        jx = jx.copy()
        jx[0, -1, -1] = rho_q * velx[0, -1, -1] + 0.2
    values = {
        "prtcl_rho": np.full(shape, rho_q),
        "prtcl_jx": jx,
        "prtcl_jy": np.zeros(shape),
        "prtcl_jz": np.zeros(shape),
        "mhd_j2": np.zeros(shape),
    }
    result: dict[str, output_primitives.AthenaBinaryDataset] = {}
    reference = mhd.blocks[0]
    for product, field in science.CURRENT_PRODUCT_FIELDS.items():
        block = output_primitives.AthenaBinaryBlock(
            index_bounds=reference.index_bounds,
            logical_location=reference.logical_location,
            level=reference.level,
            geometry=reference.geometry,
            fields={field: values[product]},
        )
        result[product] = output_primitives.AthenaBinaryDataset(
            **{
                **mhd.__dict__,
                "source": f"synthetic-q011-applicability-{product}",
                "variable_names": (field,),
                "blocks": (block,),
            }
        )
    return result


def _particles(*, velocity: object = 10.0, x: float = 5500.0) -> dict[str, object]:
    count = app.PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES
    y = ((np.arange(count) % NX2) + 0.5) * DX2
    points = np.column_stack((np.full(count, x), y, np.full(count, 0.5)))
    velocity_values = np.asarray(velocity, dtype=np.float32)
    if velocity_values.ndim == 0:
        velocity_values = np.full(count, velocity_values, dtype=np.float32)
    if velocity_values.shape != (count,):
        raise ValueError("synthetic particle velocity values must match particle count")
    velocities = np.column_stack(
        (
            velocity_values,
            np.zeros(count, dtype=np.float32),
            np.zeros(count, dtype=np.float32),
        )
    ).astype(np.float64)
    return {
        "points": points,
        "cr_source": np.ones(count, dtype=np.int64),
        "birth_time": np.full(count, 45.0),
        "velocity": velocities,
        "macro_weight": np.ones(count),
    }


def _snapshot(
    *,
    normalization: object | None = None,
    no_fluctuation: bool = False,
    low_density_cell: bool = False,
    high_k_fluctuation: bool = False,
    rho_q: float = 0.005,
    gas_frame_current: float = 0.02,
    rare_lambda_failure: bool = False,
    particle_velocity: object = 10.0,
    particle_x: float = 5500.0,
) -> app.ApplicabilitySnapshot:
    mhd = _mhd_dataset(
        no_fluctuation=no_fluctuation,
        low_density_cell=low_density_cell,
        high_k_fluctuation=high_k_fluctuation,
    )
    return app.reduce_physical_applicability_snapshot(
        mhd,
        _current_datasets(
            mhd,
            rho_q=rho_q,
            gas_frame_current=gas_frame_current,
            rare_lambda_failure=rare_lambda_failure,
        ),
        normalization=(
            dict(app.EXACT_NORMALIZATION) if normalization is None else normalization
        ),
        nominal_slot_time=TIME,
        observed_committed_time=TIME,
        **_particles(velocity=particle_velocity, x=particle_x),
    )


def _runtime(snapshot: app.ApplicabilitySnapshot) -> dict[str, object]:
    gates = snapshot.record["gates"]
    exposure = snapshot.record["particle_R_Lambda_exposure"]["populations"]["all_active"]
    return {
        "schema_version": app.SCHEMA_VERSION,
        "record_type": app.RUNTIME_RECORD_TYPE,
        "successor_id": app.SUCCESSOR_ID,
        "qualification_effect": app.QUALIFICATION_EFFECT,
        "authorization": dict(app.AUTHORIZATION),
        "exact_normalization": dict(app.EXACT_NORMALIZATION),
        "cycle_coverage": {
            "complete": True,
            "sampling_mode": "every_integrator_cycle_post_startup_removal_through_terminal_time",
            "post_startup_removal_start_time": app.STARTUP_REMOVAL_TIME,
            "terminal_time": 1200.0,
            "first_cycle": 100,
            "last_cycle": 200,
            "covered_cycle_count": 101,
            "expected_contiguous_cycle_count": 101,
            "gap_count": 0,
            "restart_segment_count": 1,
            "restart_segments_complete": True,
        },
        "all_cycle_extrema": {
            "R_maximum": gates["Q011-APP-R"]["observed_maximum"],
            "Lambda_maximum": gates["Q011-APP-LAMBDA"]["observed_maximum"],
            "S_delta_minimum_excluding_shock_transition": gates["Q011-APP-DI"][
                "observed_S_delta_minimum_excluding_shock_transition"
            ],
            "lambda_B_characteristic_over_local_di_maximum_minimum": gates[
                "Q011-APP-DI"
            ]["observed_lambda_B_characteristic_over_local_di_maximum"],
            "sub_10di_magnetic_power_fraction_maximum": gates["Q011-APP-DI"][
                "observed_sub_10di_magnetic_power_fraction"
            ],
            "macro_q999_over_Ly_maximum": gates["Q011-APP-RG"][
                "observed_macro_q999_over_Ly"
            ],
            "energy_q999_over_Ly_maximum": gates["Q011-APP-RG"][
                "observed_energy_q999_over_Ly"
            ],
            "particle_rg_maximum_over_Ly": gates["Q011-APP-RG"][
                "observed_maximum_over_Ly"
            ],
            "energy_fraction_with_rg_above_Ly_over_4_maximum": gates["Q011-APP-RG"][
                "observed_energy_fraction_with_rg_above_Ly_over_4"
            ],
        },
        "particle_exposure": {
            "complete": True,
            "method": "every_particle_update_and_pre_destruction_boundary_event",
            "active_particle_updates_included": True,
            "pre_destruction_boundary_events_included": True,
            "escaped_particles_included": True,
            "observation_count": 1000,
            "maximum_sampled_R": exposure["R"]["local_maximum"],
            "maximum_sampled_Lambda": exposure["Lambda"]["local_maximum"],
            "cumulative_macro_weighted_R_exceedance_fraction": exposure["R"][
                "macro_weighted"
            ]["exceedance_fraction"],
            "cumulative_CR_energy_weighted_R_exceedance_fraction": exposure["R"][
                "CR_kinetic_energy_weighted"
            ]["exceedance_fraction"],
            "cumulative_macro_weighted_Lambda_exceedance_fraction": exposure["Lambda"][
                "macro_weighted"
            ]["exceedance_fraction"],
            "cumulative_CR_energy_weighted_Lambda_exceedance_fraction": exposure[
                "Lambda"
            ]["CR_kinetic_energy_weighted"]["exceedance_fraction"],
        },
        "boundary_escape_ledger": {
            "complete": True,
            "scope": "all_Q011_shock_injected_particles_including_startup_removal",
            "nonperiodic_faces": {
                "ix1": {
                    "escaped_particle_count": 0,
                    "escaped_macro_weight": 0.0,
                    "escaped_kinetic_energy": 0.0,
                },
                "ox1": {
                    "escaped_particle_count": 0,
                    "escaped_macro_weight": 0.0,
                    "escaped_kinetic_energy": 0.0,
                },
            },
            "periodic_faces": ["ix2", "ox2"],
            "total_escaped_particle_count": 0,
            "total_escaped_macro_weight": 0.0,
            "total_escaped_kinetic_energy": 0.0,
            "unaccounted_particle_count": 0,
            "unaccounted_macro_weight": 0.0,
            "unaccounted_kinetic_energy": 0.0,
            "source_census_particle_count_residual": 0.0,
            "source_census_macro_weight_residual": 0.0,
            "escaped_particles_included_in_exposure": True,
        },
    }


class Q011PhysicalApplicabilitySuccessorV1Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.passing = _snapshot()

    def test_passing_snapshot_reports_required_maps_regions_gates_and_no_authority(self) -> None:
        record = self.passing.record
        self.assertTrue(record["snapshot_gate_pass_excluding_time_completeness"])
        self.assertFalse(record["gates"]["Q011-APP-TIME"]["pass"])
        self.assertFalse(record["authorization"]["launch_authorized"])
        self.assertFalse(record["authorization"]["claim_closure_authorized"])
        self.assertEqual(
            set(self.passing.cell_maps),
            {
                "R",
                "Lambda",
                "d_i",
                "S_delta",
                "actual_leaf_dx1",
                "actual_leaf_dx2",
                "gas_frame_current_magnitude",
            },
        )
        for values in self.passing.cell_maps.values():
            self.assertEqual(values.shape, (NX2, NX1))
            self.assertFalse(values.flags.writeable)
        self.assertEqual(
            set(record["regional_statistics"]), set(app.DETECTED_FRONT_REGIONS)
        )
        self.assertTrue(
            record["regional_statistics"]["detected_front_precursor"]["Lambda"][
                "gas_frame_current_weighted"
            ]["available"]
        )
        exposure = record["particle_R_Lambda_exposure"]["populations"]
        self.assertTrue(exposure["all_active"]["available"])
        self.assertTrue(exposure["detected_front_upstream"]["available"])
        self.assertTrue(exposure["high_energy_tail"]["available"])
        self.assertEqual(
            exposure["all_active"]["R"]["macro_weighted"]["exceedance_fraction"], 0.0
        )
        self.assertEqual(
            exposure["all_active"]["Lambda"]["CR_kinetic_energy_weighted"][
                "exceedance_fraction"
            ],
            0.0,
        )
        self.assertIn("all_history_and_global_applicability_claims", record["claim_rejections"])
        self.assertIn("microscopic_shock_structure_claim", record["claim_rejections"])

    def test_exact_normalization_rejects_inference_missing_fields_and_numeric_drift(self) -> None:
        missing = dict(app.EXACT_NORMALIZATION)
        missing.pop("prtcl_rho_representation")
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "keys drifted"):
            _snapshot(normalization=missing)
        drifted = dict(app.EXACT_NORMALIZATION)
        drifted["selected_species_abs_q_over_m"] = 2.0
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "drifted"):
            _snapshot(normalization=drifted)
        inferred = dict(app.EXACT_NORMALIZATION)
        inferred["selected_species_count"] = True
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "drifted"):
            _snapshot(normalization=inferred)

    def test_R_gate_uses_per_cell_maximum_and_rejects_all_physical_claims(self) -> None:
        failed = _snapshot(rho_q=0.02)
        self.assertFalse(failed.record["gates"]["Q011-APP-R"]["pass"])
        self.assertIn(
            "all_physical_MHD_PIC_Bell_shock_and_DSA_claims",
            failed.record["claim_rejections"],
        )

    def test_Lambda_gate_has_no_rare_cell_waiver(self) -> None:
        failed = _snapshot(rare_lambda_failure=True)
        gate = failed.record["gates"]["Q011-APP-LAMBDA"]
        regional = failed.record["regional_statistics"]["full_domain"]["Lambda"]
        self.assertFalse(gate["pass"])
        self.assertGreater(gate["observed_maximum"], app.LAMBDA_MAXIMUM)
        self.assertLess(
            regional["area_weighted"]["weighted_mean"], app.LAMBDA_MAXIMUM
        )
        self.assertIn("Hall_negligible_claim", failed.record["claim_rejections"])

    def test_actual_leaf_spacing_and_local_di_control_S_delta(self) -> None:
        failed = _snapshot(low_density_cell=True)
        gate = failed.record["gates"]["Q011-APP-DI"]
        self.assertFalse(gate["pass"])
        self.assertLess(
            gate["observed_S_delta_minimum_excluding_shock_transition"],
            app.S_DELTA_MINIMUM,
        )
        self.assertAlmostEqual(float(np.min(failed.cell_maps["actual_leaf_dx2"])), DX2)

    def test_sub_10di_power_and_characteristic_scale_fail_closed(self) -> None:
        failed = _snapshot(high_k_fluctuation=True)
        self.assertFalse(failed.record["gates"]["Q011-APP-DI"]["pass"])
        self.assertIn(
            "physical_precursor_turbulence_claim", failed.record["claim_rejections"]
        )

    def test_absent_precursor_magnetic_fluctuation_power_fails_closed(self) -> None:
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "no resolvable fluctuation power"
        ):
            _snapshot(no_fluctuation=True)

    def test_particle_gyroradius_gate_uses_macro_energy_tail_and_maximum(self) -> None:
        failed = _snapshot(particle_velocity=200.0)
        particle = failed.record["particle_gyroradius_containment"]
        self.assertFalse(failed.record["gates"]["Q011-APP-RG"]["pass"])
        self.assertGreater(particle["maximum_over_Ly"], app.RG_MAXIMUM_OVER_LY_MAXIMUM)
        self.assertIn("Emax_claim", failed.record["claim_rejections"])

    def test_particle_R_Lambda_exposure_reports_upstream_and_fixed_high_energy_tail(self) -> None:
        velocities = np.linspace(
            5.0, 15.0, app.PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES
        )
        snapshot = _snapshot(particle_velocity=velocities, gas_frame_current=0.2)
        populations = snapshot.record["particle_R_Lambda_exposure"]["populations"]
        self.assertEqual(
            populations["detected_front_upstream"]["particle_count"],
            app.PARTICLE_Q999_MINIMUM_POSITIVE_WEIGHT_SAMPLES,
        )
        self.assertGreater(populations["high_energy_tail"]["particle_count"], 0)
        self.assertLess(
            populations["high_energy_tail"]["particle_count"],
            populations["all_active"]["particle_count"],
        )
        self.assertEqual(
            populations["all_active"]["Lambda"]["macro_weighted"][
                "exceedance_fraction"
            ],
            1.0,
        )
        self.assertEqual(
            populations["all_active"]["Lambda"]["CR_kinetic_energy_weighted"][
                "exceedance_fraction"
            ],
            1.0,
        )

    def test_particle_TSC_nonperiodic_x1_stencil_must_remain_retained(self) -> None:
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "TSC stencil escaped"
        ):
            _snapshot(particle_x=1.0)

    def test_current_weighted_statistics_explicitly_report_zero_weight(self) -> None:
        snapshot = _snapshot(gas_frame_current=0.0)
        weighted = snapshot.record["regional_statistics"]["full_domain"]["Lambda"][
            "gas_frame_current_weighted"
        ]
        self.assertFalse(weighted["available"])
        self.assertIsNone(weighted["weighted_mean"])
        self.assertTrue(snapshot.record["gates"]["Q011-APP-LAMBDA"]["pass"])

    def test_complete_runtime_time_escape_interface_can_pass_without_authorizing(self) -> None:
        history = app.reduce_physical_applicability_history(
            [self.passing.record], _runtime(self.passing)
        )
        self.assertTrue(history["all_physical_applicability_gates_pass"])
        self.assertTrue(history["gates"]["Q011-APP-TIME"]["pass"])
        self.assertEqual(
            history["claim_rejections"],
            ["microscopic_shock_structure_claim", "self_consistent_injection_claim"],
        )
        self.assertEqual(history["authorization_effect"], "none_even_if_all_gates_pass")
        self.assertFalse(history["authorization"]["claim_closure_authorized"])

    def test_incomplete_time_coverage_and_unaccounted_escape_fail_APP_TIME(self) -> None:
        runtime = _runtime(self.passing)
        runtime["cycle_coverage"]["complete"] = False
        runtime["boundary_escape_ledger"]["unaccounted_particle_count"] = 1
        history = app.reduce_physical_applicability_history(
            [self.passing.record], runtime
        )
        self.assertFalse(history["gates"]["Q011-APP-TIME"]["pass"])
        self.assertIn(
            "all_history_and_global_applicability_claims", history["claim_rejections"]
        )

    def test_runtime_cycle_gap_and_escape_face_total_mismatch_fail_closed(self) -> None:
        runtime = _runtime(self.passing)
        runtime["cycle_coverage"]["gap_count"] = 1
        history = app.reduce_physical_applicability_history([self.passing.record], runtime)
        self.assertFalse(history["gates"]["Q011-APP-TIME"]["pass"])
        runtime = _runtime(self.passing)
        runtime["boundary_escape_ledger"]["total_escaped_particle_count"] = 1
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "do not close"):
            app.reduce_physical_applicability_history([self.passing.record], runtime)

    def test_all_cycle_Lambda_at_unity_rejects_target_plasma_claim(self) -> None:
        runtime = _runtime(self.passing)
        runtime["all_cycle_extrema"]["Lambda_maximum"] = 1.0
        history = app.reduce_physical_applicability_history(
            [self.passing.record], runtime
        )
        self.assertFalse(history["gates"]["Q011-APP-LAMBDA"]["pass"])
        self.assertIn(
            "no_Hall_run_approximates_target_plasma_claim",
            history["claim_rejections"],
        )

    def test_all_cycle_extrema_must_conservatively_dominate_supplied_snapshots(self) -> None:
        runtime = _runtime(self.passing)
        runtime["all_cycle_extrema"]["S_delta_minimum_excluding_shock_transition"] += 1.0
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "do not conservatively dominate"
        ):
            app.reduce_physical_applicability_history([self.passing.record], runtime)

    def test_history_recomputes_snapshot_pass_flags_from_bound_observables(self) -> None:
        failed = copy.deepcopy(dict(_snapshot(rho_q=0.02).record))
        failed["gates"]["Q011-APP-R"]["pass"] = True
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "pass flag disagrees"
        ):
            app.reduce_physical_applicability_history(
                [failed], _runtime(self.passing)
            )

    def test_runtime_schema_rejects_missing_fields_and_normalization_drift(self) -> None:
        runtime = _runtime(self.passing)
        runtime.pop("particle_exposure")
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "keys drifted"):
            app.reduce_physical_applicability_history([self.passing.record], runtime)
        runtime = _runtime(self.passing)
        runtime["exact_normalization"]["reference_b0"] = 2.0
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "drifted"):
            app.reduce_physical_applicability_history([self.passing.record], runtime)
        runtime = _runtime(self.passing)
        runtime["authorization"]["launch_authorized"] = 0
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "must be a boolean"):
            app.reduce_physical_applicability_history([self.passing.record], runtime)
        runtime = _runtime(self.passing)
        runtime["schema_version"] = True
        with self.assertRaisesRegex(app.PhysicalApplicabilityError, "schema version drifted"):
            app.reduce_physical_applicability_history([self.passing.record], runtime)
        runtime = _runtime(self.passing)
        runtime["particle_exposure"][
            "cumulative_macro_weighted_R_exceedance_fraction"
        ] = 0.1
        with self.assertRaisesRegex(
            app.PhysicalApplicabilityError, "reports exceedance below"
        ):
            app.reduce_physical_applicability_history([self.passing.record], runtime)

    def test_readiness_record_binds_only_new_successor_files_and_open_evidence(self) -> None:
        readiness = json.loads(READINESS.read_text(encoding="utf-8"))
        self.assertEqual(readiness["successor_id"], app.SUCCESSOR_ID)
        self.assertEqual(readiness["qualification_effect"], app.QUALIFICATION_EFFECT)
        self.assertFalse(readiness["authorization"]["launch_authorized"])
        self.assertFalse(readiness["authorization"]["claim_closure_authorized"])
        self.assertEqual(readiness["production_evidence_status"], "absent_not_fabricated")
        self.assertFalse(readiness["runtime_time_escape_telemetry_available"])
        for binding in readiness["source_bindings"].values():
            path = REPO_ROOT / binding["path"]
            self.assertEqual(binding["sha256"], _sha256(path))
        self.assertEqual(readiness["source_bindings"]["design"]["path"], str(DESIGN.relative_to(REPO_ROOT)))


if __name__ == "__main__":
    unittest.main()
