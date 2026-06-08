#!/usr/bin/env python3
"""Focused tests for the fail-closed Q019 physics-repair successor."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest
from unittest import mock

import numpy as np

from tst.publication import analyze_q019_physics_first_nonlinear_bell_successor_v2 as analysis
from tst.publication import analyze_q011_section54_outputs as binary_outputs
from tst.publication import q019_hardened_provenance_boundary_v2 as provenance
from tst.publication import q019_particle_state_analysis_bridge_v2 as bridge
from tst.publication import q019_physics_first_nonlinear_bell_successor_v2 as successor


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE = REPO_ROOT / "src/pgen/tests/q019_physics_first_nonlinear_bell_successor_v2.cpp"
HEADER = REPO_ROOT / "src/pgen/tests/q019_physics_first_nonlinear_bell_successor_v2.hpp"
MAIN = REPO_ROOT / "src/main.cpp"
HARNESS = (
    REPO_ROOT
    / "tst/publication/q019_physics_first_nonlinear_bell_successor_v2_host_harness.cpp"
)
REFERENCE_MAP = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q022_xcmp_corrected_nonlinear_bell_reference_map_successor_v2_2026-06-07.json"
)


def _synthetic_evidence(
    case_id: str, *, characteristic_gyroradius: float = 5.0, current_scale: float = 1.0
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    case = next(row for row in successor.expected_cases() if row["case_id"] == case_id)
    nx1, nx2, nx3 = (16, 8, 4) if int(case["dimension"]) == 3 else (16, 8, 1)
    shape = (nx3, nx2, nx1)
    x1_faces = np.linspace(0.0, float(case["extents"][0]), nx1 + 1)
    x2_faces = np.linspace(0.0, float(case["extents"][1]), nx2 + 1)
    x3_faces = np.linspace(0.0, float(case["extents"][2]), nx3 + 1)
    x1 = 0.5 * (x1_faces[:-1] + x1_faces[1:])
    snapshots: list[dict[str, object]] = []
    particles: list[dict[str, object]] = []
    jx = float(case["expected_j_over_c"])
    deposited_charge_density = (
        float(case["rho_cr_over_rho0"])
        * float(case["species_charge"])
        / float(case["species_mass"])
    )
    for index, time in enumerate((0.0, 0.2, 0.4)):
        deposited_moments_available = index > 0
        phase = float(case["k0"]) * x1
        by = np.broadcast_to((1.0e-5 * math.exp(time) * np.cos(phase))[None, None, :], shape)
        bz = np.broadcast_to((1.0e-5 * math.exp(time) * np.sin(phase))[None, None, :], shape)
        fields = {
            "dens": np.ones(shape),
            "eint": np.ones(shape),
            "velx": np.zeros(shape),
            "vely": np.zeros(shape),
            "velz": np.zeros(shape),
            "bcc1": np.full(shape, float(case["b_g"])),
            "bcc2": by,
            "bcc3": bz,
            "prtcl_rho": np.full(
                shape, deposited_charge_density if deposited_moments_available else 0.0
            ),
            "prtcl_jx": np.full(shape, jx if deposited_moments_available else 0.0),
            "prtcl_jy": np.zeros(shape),
            "prtcl_jz": np.zeros(shape),
            "prtcl_dedt": np.full(shape, 1.0e-5),
            "prtcl_dpxdt": np.full(shape, 1.0e-5),
            "prtcl_dpydt": np.zeros(shape),
            "prtcl_dpzdt": np.zeros(shape),
            "prtcl_ebdot": np.full(shape, 5.0e-6),
        }
        snapshots.append(
            {
                "cycle": index,
                "time": float(time),
                "x1_faces": x1_faces,
                "x2_faces": x2_faces,
                "x3_faces": x3_faces,
                "fields": fields,
            }
        )
        particles.append(
            {
                "schema_version": 2,
                "record_type": bridge.RECORD_TYPE,
                "case_id": case_id,
                "campaign_id": case["campaign_id"],
                "cycle": index,
                "time": float(time),
                "raw_restart_binding": None,
                "particle_bulk_velocity": [
                    float(case["guide_parallel_stream_speed"] - 0.02 * index),
                    0.0,
                    0.0,
                ],
                "particle_momentum": {
                    "volume_integrated_vector": [100.0 - 0.1 * index, 0.0, 0.0],
                    "momentum_flux_tensor": [
                        [2.0 - 0.1 * index, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0],
                    ],
                    "velocity_pressure_tensor": np.eye(3).tolist(),
                },
                "particle_kinetic_energy": 1000.0,
                "current": {
                    "lab_j_over_c": [current_scale * jx, 0.0, 0.0],
                    "gas_frame_j_over_c": [current_scale * jx - 0.01 * index, 0.0, 0.0],
                },
                "gyroradius": {
                    "definition": "fixture",
                    "minimum": 0.0,
                    "maximum": characteristic_gyroradius + 1.0,
                    "macro_mass_weighted_mean": characteristic_gyroradius + 0.1,
                    "macro_mass_weighted_median": characteristic_gyroradius,
                },
                "conservation": {
                    "reference_bound": True,
                    "energy_fractional_residual": 1.0e-8,
                    "momentum_fractional_residual": [1.0e-8, 0.0, 0.0],
                },
                "authority": {
                    "launch_authorized": False,
                    "policy_authorized": False,
                    "qualification_authorized": False,
                    "claim_authorized": False,
                    "raw_production_authorized": False,
                    "nonlinear_saturation_claim_authorized": False,
                },
                "limitations": ["fixture"],
            }
        )
    return snapshots, particles


class Q019PhysicsDesignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.cases = successor.expected_cases()
        cls.by_id = {str(case["case_id"]): case for case in cls.cases}
        cls.manifest, cls.rendered = successor.build_deck_manifest()

    def test_similarity_mapping_is_explicit_and_fail_closed(self) -> None:
        self.assertEqual(len(self.cases), 77)
        self.assertEqual(
            {float(case["species_charge"]) for case in self.cases}, {10000.0}
        )
        self.assertEqual(
            {float(case["background_q_over_mc_reference"]) for case in self.cases},
            {10000.0},
        )
        self.assertTrue(
            all(
                float(case["charge_density_ratio_equal_background_qom"]) <= 1.0e-5
                and float(case["guide_parallel_stream_speed"]) > float(case["u_a"])
                and float(case["minimum_active_dx_over_background_di"]) > 1.0
                and case["no_subion_cell_scale_envelope_satisfied"]
                and case["positive_finite_loading_accounting_satisfied"]
                and not case["no_hall_applicability_accepted"]
                and case["bounded_hall_omission_candidate"]
                and not case["bounded_hall_omission_review_complete"]
                and not case["mhd_resolved_scale_applicability_accepted"]
                and not case["R_much_less_than_one_applicability_accepted"]
                for case in self.cases
            )
        )
        fiducial = self.by_id["q019-fr-3d-onset-small-s0"]
        raw_eta = (
            fiducial["rho_cr_over_rho0"]
            * fiducial["species_charge"]
            / fiducial["background_q_over_mc_reference"]
        )
        self.assertAlmostEqual(
            fiducial["charge_density_ratio_equal_background_qom"],
            raw_eta / (1.0 + raw_eta),
        )
        self.assertAlmostEqual(fiducial["background_ion_inertial_length"], 1.0e-4)
        self.assertAlmostEqual(fiducial["k0_background_ion_inertial_length"], 2.0 * math.pi * 1.0e-4)
        self.assertAlmostEqual(
            fiducial["hall_parameter_from_current"],
            4.0
            * math.pi
            * 1.0e-4
            * (1.0 - fiducial["charge_density_ratio_equal_background_qom"]),
        )
        self.assertAlmostEqual(
            fiducial["hall_parameter_equal_background_qom"],
            fiducial["hall_parameter_from_current"],
        )
        self.assertAlmostEqual(
            fiducial["bai_hall_linear_factor"],
            1.0 + (fiducial["hall_parameter_from_current"] / 2.0) ** 2,
        )
        self.assertLess(fiducial["bai_hall_growth_rate_fractional_shift"], 0.0)
        self.assertLess(fiducial["bai_hall_wavenumber_fractional_shift"], 0.0)
        self.assertAlmostEqual(
            fiducial["bai_hall_growth_rate_fractional_shift"],
            fiducial["bai_hall_growth_rate_reduction_factor"] - 1.0,
        )
        self.assertAlmostEqual(
            fiducial["bai_hall_wavenumber_fractional_shift"],
            fiducial["bai_hall_wavenumber_reduction_factor"] - 1.0,
        )
        self.assertAlmostEqual(
            self.manifest["physics_grid"]["minimum_materialized_dx_over_di"],
            52.08333333333333,
        )
        self.assertFalse(
            self.manifest["literature_context"][
                "hall_order_unity_reference_is_accepted_threshold"
            ]
        )
        prerequisites = self.manifest["independent_prerequisite_boundary"]
        self.assertEqual(
            self.manifest["physical_pilot_gate_order"][:2],
            [
                f"independent_prerequisite:{successor.Q043_INDEPENDENT_ORACLE_ID}",
                f"independent_prerequisite:{successor.Q023_INDEPENDENT_PREDECESSOR_ID}",
            ],
        )
        self.assertFalse(prerequisites["q043_independent_raw_cycle_one_oracle_bound"])
        self.assertFalse(prerequisites["q023_independent_linear_predecessor_bound"])
        self.assertFalse(prerequisites["matrix_authorizes_execution"])
        self.assertEqual(
            self.manifest["diagnostic_contract"][
                "finite_packet_grouping_against_actual_initializer_proved_scope"
            ],
            "quiet_collocated_path_only",
        )

    def test_isotropic_shell_separates_current_from_shell_pressure(self) -> None:
        finite = [
            case
            for case in self.cases
            if case["branch"] != "high_rigidity_current_retention_candidate"
        ]
        self.assertEqual({int(case["ppc"]) % 6 for case in finite}, {0})
        for case in finite:
            moments = case["initial_shell_velocity_moments"]
            self.assertEqual(
                moments["mean"], [case["guide_parallel_stream_speed"], 0.0, 0.0]
            )
            covariance = np.asarray(moments["centered_covariance"])
            np.testing.assert_allclose(
                covariance,
                np.eye(3) * float(case["finite_shell_speed"]) ** 2 / 3.0,
            )
            tensor = np.asarray(case["initial_cr_momentum_flux_tensor"])
            self.assertAlmostEqual(
                tensor[0, 0] - 0.5 * (tensor[1, 1] + tensor[2, 2]),
                case["initial_anisotropic_momentum_flux_proxy"],
            )
            expected_anchors = (
                1
                if case["finite_sampling_mode"]
                == successor.FINITE_STRATIFIED_SAMPLING_MODE
                else int(case["ppc"])
            )
            self.assertEqual(
                case["finite_packet_independent_spatial_anchors_per_cell_expected"],
                expected_anchors,
            )
            quiet = (
                case["finite_sampling_mode"]
                == successor.FINITE_STRATIFIED_SAMPLING_MODE
            )
            self.assertEqual(
                case["finite_packet_grouping_against_actual_initializer_proved"],
                quiet,
            )
        fixed_rho = [
            case
            for case in finite
            if case["role"] == "finite_rigidity_density_drift_coupled_response_ensemble"
            and case["rho_cr_over_rho0"] == successor.FINITE_FIDUCIAL_RHO_CR
            and case["field_seed"] == successor.FIELD_SEEDS[0]
        ]
        self.assertEqual(len(fixed_rho), 3)
        self.assertEqual(
            {case["guide_parallel_stream_speed"] for case in fixed_rho},
            {fixed_rho[0]["guide_parallel_stream_speed"]},
        )
        self.assertEqual(
            {case["initial_anisotropic_momentum_flux_proxy"] for case in fixed_rho},
            {fixed_rho[0]["initial_anisotropic_momentum_flux_proxy"]},
        )
        self.assertEqual(len({case["isotropic_shell_pressure_proxy"] for case in fixed_rho}), 3)

    def test_current_3d_rows_are_onset_pilots_not_saturation_candidates(self) -> None:
        onset = [
            case
            for case in self.cases
            if str(case["role"]).startswith("finite_rigidity_3d_nonlinear_onset_box_")
        ]
        self.assertEqual(len(onset), 6)
        self.assertTrue(all(not case["saturation_candidate"] for case in self.cases))
        self.assertTrue(all(case["target_nonlinear_B_over_B0"] is None for case in self.cases))
        fiducial = self.by_id["q019-fr-3d-onset-small-s0"]
        self.assertAlmostEqual(fiducial["initial_nominal_rg0_over_max_active_dx"], 20.371832715762604)
        self.assertAlmostEqual(fiducial["minimum_active_extent_over_nominal_rg0"], 2.0 * math.pi)
        self.assertEqual(fiducial["required_initial_rg0_over_max_active_dx"], 8.0)
        self.assertGreaterEqual(
            fiducial["maximum_sampled_B_over_B0_before_resolution_stop"],
            fiducial["resolution_design_maximum_sampled_B_over_B0"],
        )
        all_three_d = [case for case in self.cases if case["dimension"] == 3]
        self.assertEqual(len(all_three_d), 11)
        self.assertTrue(
            all(
                not case["common_nonlinear_onset_resolution_reachable"]
                and case["common_nonlinear_onset_reachability_status"]
                == successor.COMMON_ONSET_REACHABILITY_STATUS
                for case in all_three_d
            )
        )
        self.assertTrue(
            self.by_id[
                "q019-fr-3d-onset-small-resolution-coarse-s0"
            ]["coarse_resolution_stop_or_intermittency_limitation"]
        )
        self.assertEqual(
            fiducial["characteristic_shell_p_iso_over_m"],
            fiducial["finite_shell_speed"],
        )
        self.assertGreater(
            fiducial["conservative_cr_kinetic_energy_density_upper_bound"],
            fiducial["initial_cr_kinetic_energy_density_rms_proxy"],
        )
        self.assertGreater(
            fiducial["initial_magnetic_energy_density_upper_bound"],
            0.5 * fiducial["b_g"] ** 2,
        )
        self.assertGreater(fiducial["absolute_total_energy_B_over_B0_bound"], 1.0)
        gate = self.manifest["finite_3d_nonlinear_onset_prerequisite_contract"]
        self.assertFalse(gate["current_rows_are_saturation_candidates"])
        self.assertIsNone(gate["target_B_over_B0"])
        self.assertIsNone(gate["numeric_thresholds"])
        self.assertEqual(gate["preregistered_nonlinear_onset_Bperp_rms_over_B0"], 1.0)
        self.assertEqual(gate["resolution_design_maximum_sampled_B_over_B0"], 2.0)
        self.assertFalse(gate["all_rows_reach_common_nonlinear_onset_resolution_envelope"])
        self.assertEqual(
            gate["common_nonlinear_onset_reachability_status"],
            successor.COMMON_ONSET_REACHABILITY_STATUS,
        )

    def test_box_edge_monitor_is_physical_passive_and_performance_blocked(self) -> None:
        gate = self.manifest["finite_resolution_envelope"]
        self.assertIn("physical-wave-number ball", gate["runtime_box_edge_monitor_definition"])
        self.assertNotIn("nearest nonzero", gate["runtime_box_edge_monitor_definition"])
        self.assertTrue(gate["runtime_box_edge_monitor_passive"])
        self.assertFalse(gate["runtime_box_edge_monitor_enabled"])
        self.assertEqual(gate["runtime_box_edge_monitor_enabled_case_ids"], [])
        self.assertFalse(gate["runtime_box_edge_monitor_can_request_stop"])
        self.assertTrue(gate["runtime_box_edge_monitor_restart_binds_prior_completed_time"])
        cost = self.manifest["staged_resource_boundary"]["box_edge_monitor_cost_model"]
        self.assertEqual(cost["physical_low_k_unique_modes"], {"2d": 7, "3d": 23})
        self.assertEqual(cost["full_grid_passes_per_valid_sample"], {"2d": 3, "3d": 6})
        self.assertEqual(
            cost["global_sum_reductions_per_valid_sample"], {"2d": 3, "3d": 6}
        )
        self.assertEqual(
            cost["complex_multiplications_per_cell_per_valid_sample"],
            {"2d": 28, "3d": 106},
        )
        self.assertGreater(cost["largest_case_upper_bound_full_grid_cell_visits"], 0)
        self.assertGreater(cost["largest_case_upper_bound_trigonometric_evaluations"], 0)
        self.assertEqual(cost["matrix_enabled_case_ids"], [])
        self.assertEqual(cost["matrix_enabled_case_count"], 0)
        self.assertEqual(cost["matrix_actual_full_grid_cell_visits"], 0)
        self.assertTrue(cost["excluded_pilot_benchmark_required"])
        self.assertTrue(cost["production_promotion_blocked_pending_benchmark"])
        self.assertIsNone(cost["performance_acceptance_thresholds"])
        for case in self.cases:
            dimension = int(case["dimension"])
            extents = [float(value) for value in case["extents"]]
            self.assertAlmostEqual(extents[0] / extents[1], 2.0)
            if dimension == 3:
                self.assertEqual(extents[1], extents[2])
            self.assertEqual(
                int(case["runtime_box_edge_monitor_unique_mode_count"]),
                successor.box_edge_unique_mode_count(dimension),
            )
            self.assertEqual(
                int(case["runtime_box_edge_monitor_cell_passes_per_sample"]),
                successor.box_edge_cell_passes_per_sample(dimension),
            )

    def test_decks_advance_but_are_launch_prohibited(self) -> None:
        self.assertEqual(len(self.rendered), 77)
        for case in self.cases:
            text = self.rendered[f"{case['case_id']}.athinput"]
            self.assertIn(f"nlim = {case['cycle_limit']}", text)
            self.assertIn("launch_authorized = false", text)
            self.assertIn("raw_production_authorized = false", text)
            self.assertIn("no_hall_applicability_accepted = false", text)
            self.assertIn("no_subion_cell_scale_envelope_satisfied = true", text)
            self.assertIn("positive_finite_loading_accounting_satisfied = true", text)
            self.assertIn("matrix_identity_fingerprint =", text)
            self.assertIn(
                "runtime_identity_checksum_status = "
                "sha256_cryptographic_immutable_runtime_semantics_bound_by_"
                "compiled_registry",
                text,
            )
            self.assertIn("external_sha256_execution_receipt_required = true", text)
            self.assertIn("external_sha256_execution_receipt_bound = false", text)
            self.assertIn("runtime_resolution_stop_controller_installed = false", text)
            self.assertIn("runtime_resolution_guard_pilot_qualified = false", text)
            self.assertIn("user_work_in_loop = true", text)
            self.assertIn(
                "deposited_prtcl_rho_semantics = single_species_charge_density", text
            )
            self.assertIn(
                "deposited_cr_mass_density_derivation = "
                "prtcl_rho_times_species_mass_over_species_charge",
                text,
            )
            self.assertIn(
                "species_mass_and_charge_bound_in_immutable_payload = true", text
            )
            self.assertIn("independent_prerequisites_complete = false", text)
            self.assertIn(
                "evolving_local_hall_applicability_diagnostics_complete = false", text
            )
            self.assertIn("saturation_candidate = false", text)
            self.assertIn("characteristic_shell_p_iso_over_m =", text)
            self.assertIn("runtime_box_edge_monitor_installed = true", text)
            self.assertIn("runtime_box_edge_monitor_enabled = false", text)
            self.assertIn("runtime_box_edge_stop_boundary_frozen = false", text)
            self.assertIn("runtime_box_edge_stop_armed = false", text)
            self.assertIn("runtime_box_edge_monitor_passive = true", text)
            self.assertIn("runtime_box_edge_monitor_schema = 2", text)
            self.assertIn(
                "runtime_box_edge_monitor_next_nominal_time = 0.10000000000000001",
                text,
            )
            self.assertIn("runtime_box_edge_monitor_last_prior_time = 0", text)
            self.assertIn("user_work_in_loop = true", text)
            self.assertIn("user_hist = true", text)
            self.assertIn("<mesh_refinement>\nrefinement = none", text)
            self.assertIn(
                "cr_rms_speed_kinetic_loading_proxy_to_background_magnetic_energy =",
                text,
            )
            self.assertIn("file_type = hst", text)
            self.assertIn("dcycle = 1", text)
            self.assertIn("data_format = %24.16e", text)
            self.assertEqual(text.count("ghost_zones = false"), 13)
            self.assertEqual(text.count("gid = -1"), 13)
            self.assertEqual(text.count("data_format = %12.5e"), 12)
            self.assertIn("user_hist_only = false", text)
            self.assertNotIn("hall_importance_boundary", text)
            self.assertNotIn("mhd_resolved_scale_envelope_satisfied", text)
            self.assertNotIn("loading_design_gate_passed", text)
            self.assertNotIn("gyrophase_ring", text)
            self.assertNotIn("target_nonlinear_B_over_B0", text)

    def test_checked_in_decks_match_closed_design(self) -> None:
        checked_manifest = json.loads(successor.CHECKED_IN_MANIFEST.read_text())
        self.assertEqual(checked_manifest, self.manifest)
        self.assertEqual(
            {path.name for path in successor.CHECKED_IN_DECK_ROOT.glob("*.athinput")},
            set(self.rendered),
        )
        for name, text in self.rendered.items():
            self.assertEqual((successor.CHECKED_IN_DECK_ROOT / name).read_text(), text)

    def test_runtime_resolution_guard_is_removed_and_identity_contract_is_exact(self) -> None:
        source = SOURCE.read_text()
        header = HEADER.read_text()
        self.assertNotIn("Q019RuntimeResolutionGuard", source)
        self.assertNotIn("Q019RequestIncompleteResolutionStop", source)
        self.assertNotIn("MPI_Allreduce(MPI_IN_PLACE, &local_particle_count", source)
        self.assertNotIn("MPI_Allreduce(MPI_IN_PLACE, &maximum_sampled_b_times_dx", source)
        self.assertNotIn("MPI_Allreduce(MPI_IN_PLACE, &maximum_sampled_b", source)
        self.assertIn("user_work_in_loop_func = Q019RuntimeDiagnostics", source)
        self.assertIn(
            'Q019RequireBoolean(pin, "problem", "user_work_in_loop", true)',
            source,
        )
        mutable_bookkeeping = source[
            source.index("bool Q019MutableRuntimeBookkeepingParameter") :
            source.index("std::string Q019MatrixIdentityPayload")
        ]
        identity_payload = source[
            source.index("std::string Q019MatrixIdentityPayload") :
            source.index("struct Q019CanonicalMatrixIdentity")
        ]
        self.assertIn("for (const auto &input_block : pin->block)", identity_payload)
        self.assertIn("for (const auto &input_line : input_block.line)", identity_payload)
        self.assertIn("input_line.param_value", identity_payload)
        self.assertIn("std::sort(entries.begin(), entries.end())", identity_payload)
        self.assertIn('input_block.block_name == "comment"', identity_payload)
        self.assertIn('input_line.param_name == "matrix_identity_fingerprint"', identity_payload)
        self.assertIn("Q019MutableRuntimeBookkeepingParameter", identity_payload)
        self.assertIn('parameter_name == "file_number"', mutable_bookkeeping)
        self.assertIn('parameter_name == "last_time"', mutable_bookkeeping)
        for name in successor.MUTABLE_BOX_EDGE_MONITOR_STATE_PARAMETERS:
            self.assertIn(f'parameter_name == "{name}"', mutable_bookkeeping)
        output_runtime_mutations = []
        for path in (successor.REPO_ROOT / "src/outputs").glob("*.cpp"):
            output_runtime_mutations.extend(
                re.findall(
                    r'pin->Set(?:Integer|Real|Boolean|String)\(\s*([^,]+),\s*"([^"]+)"',
                    path.read_text(),
                )
            )
        self.assertTrue(output_runtime_mutations)
        self.assertTrue(
            all(block.strip() == "out_params.block_name" for block, _ in output_runtime_mutations)
        )
        self.assertEqual(
            {name for _, name in output_runtime_mutations}, {"file_number", "last_time"}
        )
        self.assertNotIn("CharacteristicShellGyroradiusOverDx", source)
        self.assertIn("CharacteristicShellGyroradiusOverDx", header)
        self.assertIn("q019_box_edge_grouped_mode_amplitudes", source)
        self.assertIn("BoxEdgePhysicalLowKModeIsSelected", source)
        self.assertIn("BoxEdgeFirstCrossingChronologyIsValid", source)
        self.assertIn("runtime_box_edge_monitor_last_prior_time", source)
        self.assertIn("BoxEdgeDiagnosticStatus::zero_fluctuation", source)
        self.assertIn("BoxEdgeDiagnosticStatus::cadence_skipped", source)
        self.assertIn("BoxEdgeNextNominalTime", source)
        self.assertIn("BoxEdgeLastNominalTime", source)
        self.assertNotIn("q019_box_edge_mode_amplitude", source)
        self.assertNotIn("MPI_Barrier", source)
        self.assertNotIn("SampledGyroradiusOverDx", source)
        monitor = source[
            source.index("void Q019RuntimeBoxEdgeMonitor") :
            source.index("void Q019RuntimeDiagnostics")
        ]
        self.assertNotIn("Q019RequestIncompleteStop", monitor)
        self.assertNotIn("RequestUserStop", monitor)
        self.assertIn("NestedHaarOctahedralShellVelocityAtSample", source)
        self.assertIn("AxisAlignedMagneticFieldAt", source)
        self.assertIn("pgen_q019_initialize_particle_magnetic_field_cache", source)
        self.assertIn("GlobalCellLinearId", source)
        self.assertIn("cell_centered_nested_haar_octahedral_packets", source)
        self.assertNotIn("const int packet = p/finite_ppc", source)
        packet_identity = header[
            header.index("Q019_NLB_INLINE std::uint64_t NestedPacketIdentity") :
            header.index("Q019_NLB_INLINE double UniformOpen01")
        ]
        self.assertIn("global_cell_id", packet_identity)
        self.assertIn("packet_within_cell", packet_identity)
        self.assertIn("field_seed", packet_identity)
        self.assertIn("particle_seed", packet_identity)
        self.assertNotIn("PGID", packet_identity)
        self.assertNotIn("rank", packet_identity)
        self.assertIn("ppc % 6 == 0", header)
        self.assertIn("Q019_NLB_INLINE Vector3 AxisAlignedMagneticFieldAt", header)

    def test_restart_mutable_output_state_is_excluded_but_immutable_state_is_bound(self) -> None:
        case = self.by_id["q019-fr-runtime-initializer-ppc24-s0"]
        original = self.rendered[f"{case['case_id']}.athinput"]
        blocks = successor.parse_athinput_text(original)
        baseline = successor.deck_semantics_payload(blocks)
        output_name = next(name for name in blocks if name.startswith("output"))
        blocks[output_name]["file_number"] = "17"
        blocks[output_name]["last_time"] = "1.23456789"
        self.assertEqual(successor.deck_semantics_payload(blocks), baseline)
        for name in successor.MUTABLE_BOX_EDGE_MONITOR_STATE_PARAMETERS:
            original_value = blocks[successor.PGEN_BLOCK][name]
            blocks[successor.PGEN_BLOCK][name] = "17.25"
            self.assertEqual(successor.deck_semantics_payload(blocks), baseline)
            blocks[successor.PGEN_BLOCK][name] = original_value
        blocks[successor.PGEN_BLOCK]["runtime_box_edge_monitor_enabled"] = "true"
        self.assertNotEqual(successor.deck_semantics_payload(blocks), baseline)
        blocks[successor.PGEN_BLOCK]["runtime_box_edge_monitor_enabled"] = "false"
        blocks[output_name]["ghost_zones"] = "true"
        self.assertNotEqual(successor.deck_semantics_payload(blocks), baseline)
        blocks[output_name]["ghost_zones"] = "false"
        blocks["mhd"]["reconstruct"] = "wenoz"
        self.assertNotEqual(successor.deck_semantics_payload(blocks), baseline)
        validation = successor.validate_rendered_deck(case, original)
        self.assertFalse(validation["complete_runtime_deck_semantics_bound"])
        self.assertTrue(
            validation["exact_immutable_runtime_deck_semantics_payload_matches"]
        )
        self.assertTrue(
            validation[
                "immutable_runtime_deck_semantics_sha256_matches"
            ]
        )
        self.assertEqual(
            validation["immutable_runtime_deck_semantics_sha256"],
            hashlib.sha256(baseline.encode("utf-8")).hexdigest(),
        )
        self.assertTrue(
            validation["runtime_identity_checksum_is_cryptographic_integrity_binding"]
        )
        self.assertEqual(
            validation["mutable_output_bookkeeping_excluded_from_immutable_payload"],
            ["file_number", "last_time"],
        )
        self.assertEqual(
            validation[
                "mutable_box_edge_monitor_state_excluded_from_immutable_payload"
            ],
            list(successor.MUTABLE_BOX_EDGE_MONITOR_STATE_PARAMETERS),
        )
        self.assertTrue(validation["immutable_output_defaults_materialized_and_bound"])
        self.assertTrue(validation["external_sha256_execution_receipt_required"])
        self.assertFalse(validation["external_sha256_execution_receipt_bound"])
        self.assertFalse(
            validation["external_sha256_receipt_must_bind_exact_immutable_payload"]
        )
        self.assertTrue(
            validation["compiled_case_registry_binds_immutable_payload_sha256"]
        )

    def test_rendered_identity_mutations_fail_closed(self) -> None:
        case = self.by_id["q019-fr-runtime-initializer-ppc24-s0"]
        original = self.rendered[f"{case['case_id']}.athinput"]
        mutations = (
            original.replace(
                f"role = {case['role']}",
                "role = finite_rigidity_numerical_or_noise_control",
                1,
            ),
            original.replace("box_pair_id = not_applicable", "box_pair_id = drifted", 1),
            original.replace("saturation_candidate = false", "saturation_candidate = true", 1),
            original.replace(
                f"matrix_identity_fingerprint = {case['matrix_identity_fingerprint']}",
                f"matrix_identity_fingerprint = {'0' * 64}",
                1,
            ),
            original.replace("reconstruct = plm", "reconstruct = wenoz", 1),
            original.replace("pic_max_cell_cross = 2", "pic_max_cell_cross = 1", 1),
            original.replace(
                "eigenmode_amplitude = 9.9999999999999995e-07",
                "eigenmode_amplitude = 2.0000000000000001e-06",
                1,
            ),
            original.replace("x1max = 2", "x1max = 3", 1),
            original.replace("tlim = 0.0001", "tlim = 0.0002", 1),
            original.replace("dt = 0.10000000000000001", "dt = 0.2", 1),
            original.replace(
                "runtime_box_edge_monitor_dt = 0.10000000000000001",
                "runtime_box_edge_monitor_dt = 0.2",
                1,
            ),
            original.replace(
                "runtime_box_edge_monitor_next_nominal_time = 0.10000000000000001",
                "runtime_box_edge_monitor_next_nominal_time = 0.2",
                1,
            ),
            original.replace(
                "runtime_box_edge_monitor_last_prior_time = 0",
                "runtime_box_edge_monitor_last_prior_time = 0.01",
                1,
            ),
            original.replace(
                "runtime_box_edge_monitor_passive = true",
                "runtime_box_edge_monitor_passive = false",
                1,
            ),
            original.replace("refinement = none", "refinement = static", 1),
        )
        for text in mutations:
            with self.subTest(text=text.splitlines()[-1] if text else ""):
                with self.assertRaises(successor.ContractError):
                    successor.validate_rendered_deck(case, text)

    def test_coordinated_immutable_mutation_and_rehashed_deck_is_not_admitted(self) -> None:
        case = self.by_id["q019-fr-runtime-initializer-ppc24-s0"]
        original = self.rendered[f"{case['case_id']}.athinput"]
        mutated = original.replace("reconstruct = plm", "reconstruct = wenoz", 1)
        blocks = successor.parse_athinput_text(mutated)
        forged = hashlib.sha256(
            successor.deck_semantics_payload(blocks).encode("utf-8")
        ).hexdigest()
        mutated = mutated.replace(
            f"matrix_identity_fingerprint = {case['matrix_identity_fingerprint']}",
            f"matrix_identity_fingerprint = {forged}",
            1,
        )
        with self.assertRaisesRegex(
            successor.ContractError, "immutable runtime deck semantics drifted"
        ):
            successor.validate_rendered_deck(case, mutated)
        self.assertNotEqual(forged, case["matrix_identity_fingerprint"])

    def test_compiled_runtime_matrix_registry_matches_generated_matrix(self) -> None:
        source = SOURCE.read_text()
        compiled = dict(
            re.findall(
                r'\{"(q019-[^"]+)",\s*"([0-9a-f]{64})"\}',
                source[
                    source.index("q019_canonical_matrix_identities[]") :
                    source.index("Q019CanonicalMatrixFingerprint")
                ],
            )
        )
        expected = {
            str(case["case_id"]): str(case["matrix_identity_fingerprint"])
            for case in self.cases
        }
        self.assertEqual(compiled, expected)

    def test_exact_bai_literature_map_is_bound(self) -> None:
        reference = json.loads(REFERENCE_MAP.read_text())
        statements = "\n".join(reference["matched_equations"])
        self.assertIn("R=n_CR/(n_i+n_CR)", statements)
        self.assertIn("Lambda=R*(u_CR-v_g)/v_A", statements)
        self.assertIn("k_m=k0/f", statements)
        self.assertIn("growth-rate factor 1/sqrt(f)", statements)


class Q019AnalyzerAndAdmissionTests(unittest.TestCase):
    def test_transverse_large_box_fundamentals_are_unavailable_to_small_box(self) -> None:
        nx1, nx2, nx3 = 64, 32, 32
        dx = (0.5, 0.5, 0.5)
        x1 = np.arange(nx1) * dx[0]
        x2 = np.arange(nx2) * dx[1]
        x3 = np.arange(nx3) * dx[2]
        modes = (
            np.broadcast_to(np.cos(2.0 * math.pi * x1 / 32.0)[None, None, :], (nx3, nx2, nx1)),
            np.broadcast_to(np.cos(2.0 * math.pi * x2 / 16.0)[None, :, None], (nx3, nx2, nx1)),
            np.broadcast_to(np.cos(2.0 * math.pi * x3 / 16.0)[:, None, None], (nx3, nx2, nx1)),
        )
        for mode in modes:
            spectrum = analysis._spectrum(
                {"bcc2": mode, "bcc3": np.zeros_like(mode)},
                dx,
                k0=2.0 * math.pi,
                available_mode_extents=(16.0, 8.0, 8.0),
            )
            self.assertAlmostEqual(
                spectrum["large_box_only_unavailable_mode_power_fraction"], 1.0
            )
            self.assertAlmostEqual(spectrum["shared_small_box_mode_power_fraction"], 0.0)
        shared = np.broadcast_to(
            np.cos(2.0 * math.pi * x1 / 16.0)[None, None, :], (nx3, nx2, nx1)
        )
        spectrum = analysis._spectrum(
            {"bcc2": shared, "bcc3": np.zeros_like(shared)},
            dx,
            k0=2.0 * math.pi,
            available_mode_extents=(16.0, 8.0, 8.0),
        )
        self.assertAlmostEqual(spectrum["shared_small_box_mode_power_fraction"], 1.0)
        self.assertAlmostEqual(
            spectrum["large_box_only_unavailable_mode_power_fraction"], 0.0
        )

    def test_analyzer_retains_mechanism_data_and_remains_fail_closed(self) -> None:
        case_id = "q019-fr-grid-k8-rho3em06-s0"
        snapshots, particles = _synthetic_evidence(case_id)
        report = analysis.analyze_snapshots(
            case_id, snapshots, particles, source_kind="synthetic_contract_fixture"
        )
        resolution = report["evolving_finite_rigidity_resolution_stop_gate"]
        self.assertFalse(resolution["runtime_stop_controller_installed"])
        self.assertFalse(resolution["runtime_resolution_guard_pilot_qualified"])
        self.assertTrue(resolution["postprocessing_resolution_gate_required"])
        self.assertTrue(resolution["postprocessing_resolution_floor_satisfied"])
        self.assertFalse(resolution["gate_passed"])
        self.assertFalse(resolution["minimum_sampled_pitch_angle_guard_used"])
        self.assertEqual(
            resolution["trace"][0][
                "minimum_pitch_angle_sensitive_particle_rl_over_dx_report_only"
            ],
            0.0,
        )
        self.assertIn(
            "particle_momentum_flux_tensor", report["particle_state_trace"][0]
        )
        self.assertIn("conservation", report["particle_state_trace"][0])
        gate = report["finite_nonlinear_onset_and_saturation_prerequisite_gate"]
        self.assertFalse(gate["current_row_saturation_candidate"])
        self.assertIsNone(gate["target_B_over_B0"])
        self.assertIsNone(gate["numeric_thresholds"])
        self.assertFalse(gate["all_required_gates_passed"])
        mapping = report["mhd_scale_R_and_hall_applicability_gate"]
        self.assertTrue(mapping["no_subion_cell_scale_envelope_satisfied"])
        self.assertFalse(mapping["no_hall_applicability_accepted"])
        self.assertFalse(mapping["gate_passed"])
        local_hall = report["evolving_local_hall_applicability_diagnostics"]
        self.assertEqual(
            local_hall["derivation_scope"], "single_species_full_f_positive_charge"
        )
        self.assertTrue(local_hall["evolving_local_cr_concentration_and_R_available"])
        self.assertTrue(local_hall["local_cr_relative_drift_field_available"])
        self.assertTrue(local_hall["local_background_ion_inertial_length_available"])
        self.assertTrue(local_hall["local_dx_over_di_available"])
        self.assertTrue(local_hall["exact_local_Lambda_available"])
        self.assertFalse(local_hall["registered_evolving_local_hall_diagnostics_bound"])
        self.assertFalse(local_hall["gate_passed"])
        case = next(row for row in successor.expected_cases() if row["case_id"] == case_id)
        chronology = report["deposited_particle_moment_chronology"]
        self.assertTrue(chronology["natural_cycle_zero_pre_deposition_output_present"])
        self.assertEqual(chronology["first_valid_deposited_moment_cycle"], 1)
        self.assertEqual(
            chronology["deposited_moment_denominator_and_noise_baseline_cycle"], 1
        )
        self.assertFalse(
            report["grid_trace"][0]["deposited_particle_moment_state"]["available"]
        )
        self.assertAlmostEqual(
            report["grid_trace"][1]["deposited_grid_charge_density"],
            float(case["rho_cr_over_rho0"])
            * float(case["species_charge"])
            / float(case["species_mass"]),
        )
        self.assertAlmostEqual(
            report["grid_trace"][1]["deposited_grid_cr_mass_density"],
            float(case["rho_cr_over_rho0"]),
        )
        noise_gate = report["initial_rho_jx_noise_pair_gate"]
        self.assertTrue(
            noise_gate["first_snapshot_measurement_is_first_valid_deposited_moment"]
        )
        self.assertEqual(noise_gate["first_valid_deposited_moment_cycle"], 1)
        self.assertAlmostEqual(
            noise_gate["first_snapshot_measurement"]["deposited_charge_density_mean"],
            report["grid_trace"][1]["deposited_grid_charge_density"],
        )
        self.assertAlmostEqual(
            local_hall["trace"][0]["local_cr_mass_density"]["mean"],
            float(case["rho_cr_over_rho0"]),
        )
        self.assertAlmostEqual(
            local_hall["trace"][0]["local_Bai_R"]["mean"],
            float(case["charge_density_ratio_equal_background_qom"]),
        )
        self.assertAlmostEqual(
            local_hall["trace"][0]["signed_local_Lambda"]["mean"],
            float(case["hall_parameter_equal_background_qom"]),
            places=7,
        )
        prerequisite = report["independent_prerequisite_gate"]
        self.assertEqual(
            prerequisite["q043_independent_raw_cycle_one_oracle_id"],
            successor.Q043_INDEPENDENT_ORACLE_ID,
        )
        self.assertFalse(prerequisite["gate_passed"])
        timestep = report["actual_timestep_history_gate"]
        self.assertTrue(timestep["per_cycle_history_retention_configured"])
        self.assertFalse(timestep["registered_history_binding_present"])
        self.assertIsNone(timestep["actual_timestep_convergence_gate_passed"])

    def test_seeded_k0rg_rows_use_exact_case_k0_b0_and_mode_amplitude(self) -> None:
        for k0rg in (4, 16):
            case_id = f"q019-fr-rigidity-isolation-k{k0rg}-s0"
            case = next(row for row in successor.expected_cases() if row["case_id"] == case_id)
            snapshots, particles = _synthetic_evidence(case_id)
            report = analysis.analyze_snapshots(
                case_id,
                snapshots,
                particles,
                source_kind="synthetic_contract_fixture",
            )
            initial = report["grid_trace"][0]
            self.assertAlmostEqual(initial["spectrum"]["dominant_k_over_k0"], 1.0)
            self.assertAlmostEqual(
                initial["complex_Bperp_k0"]["amplitude"], 1.0e-5, places=12
            )
            self.assertAlmostEqual(
                initial["mean_Bperp_rms_over_B0"], 1.0e-5 / float(case["b_g"])
            )

    def test_characteristic_not_minimum_pitch_angle_fails_postprocessing_gate(self) -> None:
        case_id = "q019-fr-grid-k8-rho3em06-s0"
        snapshots, particles = _synthetic_evidence(case_id, characteristic_gyroradius=1.0)
        report = analysis.analyze_snapshots(
            case_id, snapshots, particles, source_kind="synthetic_contract_fixture"
        )
        gate = report["evolving_finite_rigidity_resolution_stop_gate"]
        self.assertFalse(gate["gate_passed"])
        self.assertFalse(gate["postprocessing_resolution_floor_satisfied"])
        self.assertEqual(gate["first_postprocessing_violation_time"], 0.0)

    def test_cycle_zero_empty_moments_are_natural_but_later_empty_moments_fail(self) -> None:
        case_id = "q019-fr-grid-k8-rho3em06-s0"
        snapshots, particles = _synthetic_evidence(case_id)
        for name in ("prtcl_rho", "prtcl_jx", "prtcl_jy", "prtcl_jz"):
            snapshots[1]["fields"][name] = np.zeros_like(snapshots[1]["fields"][name])
        with self.assertRaisesRegex(
            analysis.ContractError, "only at exact cycle-zero startup"
        ):
            analysis.analyze_snapshots(
                case_id,
                snapshots,
                particles,
                source_kind="synthetic_contract_fixture",
            )

    def test_cycle_zero_plus_one_valid_deposit_is_a_valid_minimal_chronology(self) -> None:
        case_id = "q019-fr-grid-k8-rho3em06-s0"
        snapshots, particles = _synthetic_evidence(case_id)
        report = analysis.analyze_snapshots(
            case_id,
            snapshots[:2],
            particles[:2],
            source_kind="synthetic_contract_fixture",
        )
        chronology = report["deposited_particle_moment_chronology"]
        self.assertEqual(chronology["valid_deposited_moment_snapshot_count"], 1)
        self.assertEqual(chronology["first_valid_deposited_moment_cycle"], 1)

    def test_grid_particle_cycle_time_chronology_is_exact(self) -> None:
        case_id = "q019-fr-grid-k8-rho3em06-s0"
        snapshots, particles = _synthetic_evidence(case_id)
        particles[1]["cycle"] = 7
        with self.assertRaisesRegex(analysis.ContractError, "cycle/time chronology"):
            analysis.analyze_snapshots(
                case_id,
                snapshots,
                particles,
                source_kind="synthetic_contract_fixture",
            )

    def test_local_diagnostic_availability_uses_physical_coverage_domains(self) -> None:
        case_id = "q019-fr-grid-k8-rho3em06-s0"
        snapshots, particles = _synthetic_evidence(case_id)
        for snapshot in snapshots[1:]:
            for name in ("prtcl_rho", "prtcl_jx", "prtcl_jy", "prtcl_jz"):
                snapshot["fields"][name][0, 0, 0] = 0.0
        report = analysis.analyze_snapshots(
            case_id, snapshots, particles, source_kind="synthetic_contract_fixture"
        )
        local = report["evolving_local_hall_applicability_diagnostics"]
        first = local["trace"][0]
        self.assertEqual(
            first["local_concentration_and_Bai_R_cell_coverage_fraction"], 1.0
        )
        self.assertEqual(
            first[
                "local_background_ion_inertial_length_and_dx_over_di_cell_coverage_fraction"
            ],
            1.0,
        )
        self.assertLess(first["local_cr_relative_drift_cell_coverage_fraction"], 1.0)
        self.assertLess(first["local_Lambda_cell_coverage_fraction"], 1.0)
        self.assertFalse(
            local["all_valid_snapshot_cells_have_local_diagnostics"]
        )
        self.assertTrue(local["unavailable_local_cr_drift_limited_to_zero_charge_cells"])
        self.assertTrue(
            local["unavailable_local_Lambda_limited_to_zero_charge_or_zero_B_cells"]
        )

    def test_structured_incomplete_stop_is_rejected_even_with_scheduler_completion(self) -> None:
        record = {
            "record_type": provenance.RUNTIME_COMPLETION_RECORD_TYPE,
            "run_completion_status": "incomplete_rejected_resolution_guard",
            "problem_stop_requested": True,
            "stop_reason_code": "characteristic_shell_resolution_floor",
            "process_exit_code": 1,
            "scheduler_terminal_state": "COMPLETED",
            "trusted_execution_binding_present": True,
        }
        classified = provenance.classify_runtime_completion(record)
        self.assertTrue(classified["incomplete_rejected"])
        self.assertEqual(classified["evidence_disposition"], "rejected_incomplete")
        self.assertFalse(classified["saturation_evidence_eligible"])
        bad = dict(record, process_exit_code=0)
        with self.assertRaisesRegex(provenance.ProvenanceBoundaryError, "nonzero"):
            provenance.classify_runtime_completion(bad)

    def test_raw_science_requires_hardened_registered_admission(self) -> None:
        contract = provenance.contract()
        self.assertIn(
            "structured_runtime_completion_status_and_problem_stop_reason",
            contract["required_bindings"],
        )
        self.assertTrue(contract["source_local_raw_science_analysis_enabled"])
        with self.assertRaisesRegex(
            analysis.ContractError, "requires a valid hardened"
        ):
            analysis.validate_raw_provenance({"self_attested": True})

    def test_raw_analysis_binds_exact_bundle_and_passed_independent_prerequisites(
        self,
    ) -> None:
        case_id = "q019-fr-grid-k8-rho3em06-s0"
        snapshots, particles = _synthetic_evidence(case_id)
        completion = {
            "record_type": provenance.RUNTIME_COMPLETION_RECORD_TYPE,
            "run_completion_status": "completed_not_acceptance_eligible",
            "problem_stop_requested": False,
            "stop_reason_code": "Terminating on time limit",
            "process_exit_code": 0,
            "scheduler_terminal_state": "COMPLETED",
            "trusted_execution_binding_present": True,
        }
        hardened = {
            "raw_science_admission_eligible": True,
            "saturation_evidence_eligible": False,
        }
        with mock.patch.object(
            provenance,
            "validate_raw_science_bundle",
            return_value=hardened,
        ) as validate:
            report = analysis.analyze_snapshots(
                case_id,
                snapshots,
                particles,
                source_kind="raw_registered_bundle",
                provenance={"record_type": provenance.REQUIRED_ADAPTER_RECORD_TYPE},
                completion_record=completion,
            )
        validate.assert_called_once_with(
            {"record_type": provenance.REQUIRED_ADAPTER_RECORD_TYPE},
            snapshots=snapshots,
            particle_states=particles,
            completion_record=completion,
        )
        gate = report["structured_runtime_completion_gate"]
        self.assertTrue(gate["raw_science_admission_eligible"])
        self.assertFalse(gate["saturation_evidence_eligible"])
        self.assertEqual(
            gate["evidence_disposition"], "admitted_registered_raw_analysis"
        )
        prerequisite = report["independent_prerequisite_gate"]
        self.assertTrue(prerequisite["q043_independent_raw_cycle_one_oracle_bound"])
        self.assertTrue(prerequisite["q023_independent_linear_predecessor_bound"])
        self.assertTrue(prerequisite["gate_passed"])


class Q019HostHarnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temp = tempfile.TemporaryDirectory()
        binary = Path(cls.temp.name) / "q019-host"
        subprocess.run(
            [
                os.environ.get("CXX", "c++"),
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-Werror",
                str(HARNESS),
                "-o",
                str(binary),
            ],
            cwd=REPO_ROOT,
            check=True,
        )
        cls.lines = subprocess.run(
            [str(binary)], cwd=REPO_ROOT, check=True, text=True, stdout=subprocess.PIPE
        ).stdout.splitlines()

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temp.cleanup()

    def test_mapping_packet_and_shell_oracles(self) -> None:
        mapping = next(line.split() for line in self.lines if line.startswith("mapping "))
        self.assertAlmostEqual(float(mapping[1]), 1.0e-4)
        self.assertAlmostEqual(float(mapping[3]), successor.FIDUCIAL_HALL_PARAMETER)
        self.assertAlmostEqual(float(mapping[7]), successor.FIDUCIAL_EXACT_BAI_R)
        self.assertAlmostEqual(float(mapping[8]), successor.FIDUCIAL_HALL_PARAMETER)
        self.assertLess(float(mapping[5]), 0.0)
        self.assertLess(float(mapping[6]), 0.0)
        self.assertAlmostEqual(float(mapping[5]), float(mapping[9]) - 1.0)
        self.assertAlmostEqual(float(mapping[6]), float(mapping[10]) - 1.0)
        current = next(line.split() for line in self.lines if line.startswith("current "))
        self.assertAlmostEqual(float(current[1]), successor.EXPECTED_J_OVER_C)
        scale = next(line.split() for line in self.lines if line.startswith("scale "))
        self.assertEqual(scale[1:3], ["1", "0"])
        packet = next(line.split() for line in self.lines if line.startswith("packet "))
        self.assertEqual(packet[1:], ["1", "0"])
        resolution = next(line.split() for line in self.lines if line.startswith("resolution "))
        self.assertEqual(resolution[1:3], ["1", "0"])
        selection = next(
            line.split() for line in self.lines if line.startswith("box_edge_selection ")
        )
        self.assertEqual(selection[1:3], ["7", "23"])
        self.assertEqual(selection[3:], ["1", "0", "1", "0", "1", "1"])
        numeric = next(
            line.split() for line in self.lines if line.startswith("box_edge_numeric ")
        )
        self.assertAlmostEqual(float(numeric[1]), 1.0, places=12)
        self.assertAlmostEqual(float(numeric[2]), 0.0, places=12)
        self.assertAlmostEqual(float(numeric[3]), 1.0, places=12)
        self.assertAlmostEqual(float(numeric[4]), 0.0, places=12)
        self.assertLess(float(numeric[5]), 1.0e-14)
        self.assertAlmostEqual(float(numeric[6]), 1024.0, places=10)
        self.assertAlmostEqual(float(numeric[7]), 1024.0**2 / 2.0, places=7)
        schedule = next(
            line.split() for line in self.lines if line.startswith("box_edge_schedule ")
        )
        self.assertEqual(float(schedule[1]), -1.0)
        self.assertEqual(int(schedule[2]), 2)
        self.assertEqual(schedule[3:5], ["1", "0"])
        self.assertAlmostEqual(float(schedule[5]), 0.4)
        self.assertAlmostEqual(float(schedule[6]), 0.5)
        self.assertEqual(schedule[7:], ["1", "0", "4", "16"])
        exact_cadence = next(
            line.split()
            for line in self.lines
            if line.startswith("box_edge_exact_cadence ")
        )
        self.assertEqual(exact_cadence[1:6], ["1", "1", "1", "1", "0"])
        self.assertAlmostEqual(float(exact_cadence[6]), 0.4)
        self.assertAlmostEqual(float(exact_cadence[7]), 0.5)
        analytic_field = next(
            line.split() for line in self.lines if line.startswith("analytic_field ")
        )
        self.assertLess(float(analytic_field[1]), 1.0e-9)
        self.assertLess(float(analytic_field[2]), 1.0e-9)
        self.assertGreater(float(analytic_field[3]), 0.999)
        self.assertGreater(float(analytic_field[4]), 0.999)
        nested = next(line.split() for line in self.lines if line.startswith("nested "))
        self.assertEqual(nested[1], "1")
        decomposition = next(
            line.split() for line in self.lines if line.startswith("decomposition ")
        )
        self.assertEqual(decomposition[1], "1")
        self.assertEqual(float(decomposition[2]), 0.0)
        moments = [
            line.split() for line in self.lines if line.startswith("packet_moments ")
        ]
        self.assertEqual([int(row[1]) for row in moments], list(successor.FINITE_PPC_LADDER))
        shell_speed = successor.finite_shell_speed(successor.FINITE_3D_FIDUCIAL_K0_RG0)
        eps = np.finfo(float).eps
        for row in moments:
            self.assertLessEqual(float(row[2]), 512.0 * eps * max(1.0, shell_speed))
            self.assertLessEqual(
                float(row[3]), 512.0 * eps * max(1.0, shell_speed * shell_speed)
            )
        fourth = [line.split() for line in self.lines if line.startswith("fourth ")]
        self.assertEqual([int(row[1]) for row in fourth], list(successor.FINITE_PPC_LADDER))
        rms = [float(row[2]) for row in fourth]
        self.assertGreater(rms[0], rms[1])
        self.assertGreater(rms[1], rms[2])


@unittest.skipUnless(
    os.environ.get("ATHENA_Q019_EXE_DIR"),
    "ATHENA_Q019_EXE_DIR is required for finite-deck runtime regression",
)
class Q019RuntimeInitializerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.binary = Path(os.environ["ATHENA_Q019_EXE_DIR"]) / "athena"
        if not cls.binary.is_file():
            raise RuntimeError(f"Q019 runtime executable does not exist: {cls.binary}")
        cls.temp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temp.cleanup()

    def _run(
        self,
        case_id: str,
        overrides: list[str],
        *,
        restart_file: Path | None = None,
    ) -> tuple[subprocess.CompletedProcess[str], Path]:
        deck = successor.CHECKED_IN_DECK_ROOT / f"{case_id}.athinput"
        run_dir = Path(self.temp.name) / f"{case_id}-{len(list(Path(self.temp.name).iterdir()))}"
        run_dir.mkdir()
        command = [str(self.binary)]
        command.extend(
            ["-r", str(restart_file)]
            if restart_file is not None
            else ["-i", str(deck)]
        )
        result = subprocess.run(
            [*command, *overrides],
            cwd=run_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=180,
            check=False,
        )
        return result, run_dir

    @staticmethod
    def _field_values(path: Path, field: str) -> tuple[int, np.ndarray]:
        dataset = binary_outputs.read_athenak_binary(path)
        values = np.concatenate(
            [np.asarray(block.fields[field], dtype=np.float64).ravel() for block in dataset.blocks]
        )
        return dataset.cycle, values

    def test_representative_finite_initializers_execute(self) -> None:
        for ppc in (24, 96):
            case_id = f"q019-fr-runtime-initializer-ppc{ppc}-s0"
            result, _ = self._run(case_id, [])
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn(
                "Q019_FINAL_EVIDENCE_STATUS=completed_not_acceptance_eligible",
                result.stdout,
            )
            self.assertIn("Q019_BOX_EDGE_MONITOR_COMPLETED_SLOTS=0", result.stdout)
            self.assertIn("Q019_BOX_EDGE_MONITOR_VALID_SAMPLES=0", result.stdout)

    def test_fresh_checkpoint_restart_and_hostile_immutable_override(self) -> None:
        case_id = "q019-fr-runtime-initializer-ppc24-s0"
        fresh, fresh_dir = self._run(case_id, [])
        self.assertEqual(fresh.returncode, 0, fresh.stdout)
        restart = fresh_dir / "rst" / f"{case_id}.00000.rst"
        self.assertTrue(restart.is_file())
        self.assertTrue(Path(f"{restart}.manifest").is_file())
        restarted, _ = self._run(case_id, [], restart_file=restart)
        self.assertEqual(restarted.returncode, 0, restarted.stdout)
        self.assertIn(
            "Q019_FINAL_EVIDENCE_STATUS=completed_not_acceptance_eligible",
            restarted.stdout,
        )
        hostile, _ = self._run(
            case_id,
            [
                f"{successor.PGEN_BLOCK}/"
                "role=finite_rigidity_numerical_or_noise_control"
            ],
            restart_file=restart,
        )
        self.assertNotEqual(hostile.returncode, 0, hostile.stdout)
        self.assertIn("immutable runtime deck semantics checksum mismatch", hostile.stdout)
        hostile_output, _ = self._run(
            case_id,
            ["output1/ghost_zones=true"],
            restart_file=restart,
        )
        self.assertNotEqual(hostile_output.returncode, 0, hostile_output.stdout)
        self.assertIn(
            "immutable runtime deck semantics checksum mismatch", hostile_output.stdout
        )

    def test_real_cycle_zero_empty_moments_and_charge_density_semantics(self) -> None:
        case_id = "q019-fr-runtime-initializer-ppc24-s0"
        case = next(
            row for row in successor.expected_cases() if row["case_id"] == case_id
        )
        result, run_dir = self._run(case_id, [])
        self.assertEqual(result.returncode, 0, result.stdout)
        rho_by_cycle = dict(
            self._field_values(path, "prtcl_rho")
            for path in sorted((run_dir / "bin").glob(f"{case_id}.prtcl_rho.*.bin"))
        )
        jx_by_cycle = dict(
            self._field_values(path, "prtcl_jx")
            for path in sorted((run_dir / "bin").glob(f"{case_id}.prtcl_jx.*.bin"))
        )
        self.assertIn(0, rho_by_cycle)
        self.assertIn(1, rho_by_cycle)
        self.assertTrue(np.all(rho_by_cycle[0] == 0.0))
        self.assertTrue(np.all(jx_by_cycle[0] == 0.0))
        expected_charge_density = (
            float(case["rho_cr_over_rho0"])
            * float(case["species_charge"])
            / float(case["species_mass"])
        )
        self.assertAlmostEqual(float(np.mean(rho_by_cycle[1])), expected_charge_density)
        self.assertAlmostEqual(
            float(np.mean(rho_by_cycle[1]))
            * float(case["species_mass"])
            / float(case["species_charge"]),
            float(case["rho_cr_over_rho0"]),
        )
        self.assertNotAlmostEqual(
            float(np.mean(rho_by_cycle[1])), float(case["rho_cr_over_rho0"])
        )
        self.assertAlmostEqual(
            float(np.mean(jx_by_cycle[1])), float(case["expected_j_over_c"])
        )

    def test_runtime_identity_mutations_are_rejected(self) -> None:
        block = successor.PGEN_BLOCK
        base = dict(
            next(
                case
                for case in successor.expected_cases()
                if case["case_id"] == "q019-fr-runtime-initializer-ppc24-s0"
            )
        )
        fabricated_case = dict(base)
        fabricated_case["case_id"] = "q019-fabricated-runtime-initializer-ppc24-s0"
        fabricated_role = dict(base)
        fabricated_role["role"] = "finite_rigidity_numerical_or_noise_control"
        cases = (
            [
                "job/basename=q019-mutated",
                f"{block}/case_id=q019-mutated",
            ],
            [
                f"job/basename={fabricated_case['case_id']}",
                f"{block}/case_id={fabricated_case['case_id']}",
                f"{block}/matrix_identity_fingerprint="
                f"{successor.matrix_identity_fingerprint(fabricated_case)}",
            ],
            [f"{block}/role=finite_rigidity_numerical_or_noise_control"],
            [
                f"{block}/role={fabricated_role['role']}",
                f"{block}/matrix_identity_fingerprint="
                f"{successor.matrix_identity_fingerprint(fabricated_role)}",
            ],
            [f"{block}/box_pair_id=drifted"],
            [f"{block}/saturation_candidate=true"],
            ["mhd/reconstruct=wenoz"],
            [f"{block}/eigenmode_amplitude=2e-6"],
            [
                "particles/pic_theta_max=1",
                f"{block}/pic_theta_max=1",
            ],
            [
                f"{block}/minimum_characteristic_shell_rl_over_dx=1",
                f"{block}/required_initial_rg0_over_max_active_dx=1",
            ],
            ["species0/mass=2"],
            ["species0/charge=5000"],
            [f"{block}/runtime_box_edge_monitor_dt=0.2"],
            [f"{block}/runtime_box_edge_monitor_enabled=true"],
            [f"{block}/runtime_box_edge_monitor_next_nominal_time=0.2"],
            [f"{block}/runtime_box_edge_monitor_last_prior_time=0.01"],
            [f"{block}/runtime_box_edge_monitor_passive=false"],
            ["mesh_refinement/refinement=static"],
        )
        for overrides in cases:
            with self.subTest(overrides=overrides):
                result, _ = self._run(
                    "q019-fr-runtime-initializer-ppc24-s0", overrides
                )
                self.assertNotEqual(result.returncode, 0, result.stdout)

    def test_restart_mutable_monitor_state_is_chronology_validated(self) -> None:
        case_id = "q019-fr-runtime-initializer-ppc24-s0"
        run_dir = Path(self.temp.name) / "restart-prior-time"
        run_dir.mkdir()
        first = subprocess.run(
            [str(self.binary), "-i", str(
                successor.CHECKED_IN_DECK_ROOT / f"{case_id}.athinput"
            )],
            cwd=run_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=180,
            check=False,
        )
        self.assertEqual(first.returncode, 0, first.stdout)
        restart = run_dir / "rst" / f"{case_id}.00000.rst"
        self.assertTrue(restart.is_file())
        continuation = subprocess.run(
            [str(self.binary), "-r", str(restart)],
            cwd=run_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=180,
            check=False,
        )
        self.assertEqual(continuation.returncode, 0, continuation.stdout)
        self.assertIn("Q019_BOX_EDGE_MONITOR_COMPLETED_SLOTS=0", continuation.stdout)
        hostile = subprocess.run(
            [
                str(self.binary),
                "-r",
                str(restart),
                f"{successor.PGEN_BLOCK}/runtime_box_edge_monitor_last_prior_time=0.01",
            ],
            cwd=run_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=180,
            check=False,
        )
        self.assertNotEqual(hostile.returncode, 0, hostile.stdout)
        self.assertIn("Q019 box-edge monitor restart chronology drifted", hostile.stdout)


if __name__ == "__main__":
    unittest.main()
