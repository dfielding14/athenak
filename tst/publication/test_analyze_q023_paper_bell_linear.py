#!/usr/bin/env python3
"""Focused tests for the Q-023 Section 5.2 Bell source-local candidates."""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from tst.publication import analyze_q023_paper_bell_linear as bell


REPO_ROOT = Path(__file__).resolve().parents[2]
PROFILES = (
    REPO_ROOT
    / "tst/publication/readiness/q023_local_preregistration_profiles_2026-05-30.json"
)
DRAFTS = (
    REPO_ROOT / "tst/publication/readiness/q023_campaign_drafts_2026-05-30.json"
)
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/q023_paper_bell_linear_source_local_implementation_successor_2026-05-31.json"
)
PAPER_INPUT_IDS = {
    "Q023-INPUT-PAPER-BELL-LINEAR-1D-CANDIDATE",
    "Q023-INPUT-PAPER-BELL-LINEAR-2D-CANDIDATE",
    "Q023-INPUT-PAPER-BELL-LINEAR-3D-CANDIDATE",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(dimension: int, epsilon: float) -> dict[str, object]:
    _, growth = bell.theoretical_dispersion(epsilon)
    time = np.linspace(0.0, 6.0 / growth, 49)
    right_amplitude = np.exp(growth * time)
    right = right_amplitude * np.exp(1.0j * epsilon * time)
    left = 1.0e-3 * right
    velocity_right = (-epsilon - 1.0j * growth) * right
    velocity_left = 1.0e-3 * velocity_right
    phase_interval, phase_change = bell._fixed_interval_phase_trace(time, right)
    velocity_phase_interval, velocity_phase_change = bell._fixed_interval_phase_trace(
        time, velocity_right
    )
    paper_phase_interval, paper_phase_change = bell._fixed_interval_phase_trace(
        time, velocity_right
    )
    return {
        "dimension": dimension,
        "epsilon": epsilon,
        "raw_provenance": bell.synthetic_contract_provenance(dimension, epsilon),
        "normalized_time": time.tolist(),
        "right_mode_real": right.real.tolist(),
        "right_mode_imag": right.imag.tolist(),
        "left_mode_real": left.real.tolist(),
        "left_mode_imag": left.imag.tolist(),
        "velocity_right_mode_real": velocity_right.real.tolist(),
        "velocity_right_mode_imag": velocity_right.imag.tolist(),
        "velocity_left_mode_real": velocity_left.real.tolist(),
        "velocity_left_mode_imag": velocity_left.imag.tolist(),
        "phase_interval": phase_interval.tolist(),
        "phase_change": phase_change.tolist(),
        "velocity_phase_interval": velocity_phase_interval.tolist(),
        "velocity_phase_change": velocity_phase_change.tolist(),
        "paper_delta_u_y_sine_fit_real": velocity_right.real.tolist(),
        "paper_delta_u_y_sine_fit_imag": velocity_right.imag.tolist(),
        "paper_delta_u_y_phase_interval": paper_phase_interval.tolist(),
        "paper_delta_u_y_phase_change": paper_phase_change.tolist(),
        "paper_volume_averaged_abs_delta_u": np.abs(velocity_right).tolist(),
    }


def _bundle() -> dict[str, object]:
    return {
        "schema_version": bell.TRACE_SCHEMA_VERSION,
        "campaign_id": bell.CAMPAIGN_ID,
        "qualifying_seed": bell.QUALIFYING_SEEDS[0],
        "records": [
            _record(dimension, epsilon)
            for dimension in bell.DIMENSIONS
            for epsilon in bell.EPSILON_VALUES
        ],
    }


def _raw_datasets(dimension: int, epsilon: float) -> list[dict[str, object]]:
    _, growth = bell.theoretical_dispersion(epsilon)
    parallel, transverse_a, transverse_b = bell._mode_basis(dimension)
    geometry = bell._EXPECTED_GEOMETRY[dimension]
    coordinates = [
        xmin + (np.arange(count, dtype=float) + 0.5) * extent / count
        for count, xmin, extent in zip(
            geometry["nx"], geometry["xmin"], geometry["extent"]
        )
    ]
    x3, x2, x1 = np.meshgrid(
        coordinates[2], coordinates[1], coordinates[0], indexing="ij"
    )
    spatial_phase = bell.K0 * (
        parallel[0]*x1 + parallel[1]*x2 + parallel[2]*x3
    )
    datasets = []
    for normalized_time in np.linspace(0.0, 7.0, 57):
        amplitude = 1.0e-6 * math.exp(growth*normalized_time)
        temporal_phase = spatial_phase + epsilon*normalized_time
        magnetic = (
            parallel[:, None, None, None]
            + amplitude*np.cos(temporal_phase)[None, ...]
            * transverse_a[:, None, None, None]
            - amplitude*np.sin(temporal_phase)[None, ...]
            * transverse_b[:, None, None, None]
        )
        velocity = amplitude * (
            (-epsilon*np.cos(temporal_phase) + growth*np.sin(temporal_phase))[None, ...]
            * transverse_a[:, None, None, None]
            + (growth*np.cos(temporal_phase) + epsilon*np.sin(temporal_phase))[None, ...]
            * transverse_b[:, None, None, None]
        )
        datasets.append(
            {
                "Time": normalized_time / (bell.K0*bell.U_A),
                "x1v": coordinates[0],
                "x2v": coordinates[1],
                "x3v": coordinates[2],
                "bcc1": magnetic[0],
                "bcc2": magnetic[1],
                "bcc3": magnetic[2],
                "velx": velocity[0],
                "vely": velocity[1],
                "velz": velocity[2],
            }
        )
    return datasets


class Q023PaperBellLinearTests(unittest.TestCase):
    def test_source_local_decks_freeze_paper_values_for_non_authorized_local_run(
        self,
    ) -> None:
        decks = bell.validate_source_local_candidate_decks()
        self.assertEqual([deck["dimension"] for deck in decks], [1, 2, 3])
        self.assertEqual(
            [deck["launch_status"] for deck in decks],
            [
                "source_local_runnable_thin_2d3v_carrier_preparation_only_not_authorized",
                "source_local_runnable_preparation_only_not_authorized",
                "source_local_runnable_preparation_only_not_authorized",
            ],
        )
        self.assertEqual(
            decks[0]["carrier_semantics"],
            "physical_1d_transverse_invariant_thin_2d3v_nx2_4",
        )
        self.assertEqual(bell._EXPECTED_GEOMETRY[1]["nx"], (32, 4, 1))
        self.assertEqual(decks[0]["active_dx"], [1.0 / 32.0])
        self.assertAlmostEqual(decks[1]["active_dx"][0], math.sqrt(5.0) / 64.0)
        self.assertAlmostEqual(decks[2]["active_dx"][2],
                               math.sqrt(1.3125) / 32.0)

    def test_complete_synthetic_analytical_grid_passes(self) -> None:
        report = bell.analyze_trace_bundle(_bundle())
        self.assertFalse(report["passed"])
        self.assertTrue(report["scientific_contract_pass"])
        self.assertEqual(report["record_count"], 15)
        self.assertEqual(report["analysis_input_kind"], "synthetic_contract_fixture")
        self.assertEqual(report["analysis_scope"], "synthetic_contract_fixture_test_only")
        self.assertTrue(report["synthetic_fixture_analysis"])
        self.assertFalse(report["materialized_source_local_candidate_pass"])
        self.assertFalse(report["section52_qualification_eligible"])
        self.assertIn("synthetic_contract_fixture_test_only",
                      report["qualification_effect"])

    def test_growth_fit_ignores_post_window_amplitude_change(self) -> None:
        bundle = _bundle()
        record = bundle["records"][0]
        _, growth = bell.theoretical_dispersion(record["epsilon"])
        time = np.asarray(record["normalized_time"])
        right = (
            np.asarray(record["right_mode_real"])
            + 1.0j * np.asarray(record["right_mode_imag"])
        )
        late = time > (5.0 / growth)
        right[late] *= np.exp(20.0 * (time[late] - 5.0 / growth))
        record["right_mode_real"] = right.real.tolist()
        record["right_mode_imag"] = right.imag.tolist()
        report = bell.analyze_trace_bundle(bundle)
        self.assertFalse(report["passed"])
        self.assertTrue(report["scientific_contract_pass"])

    def test_incomplete_or_duplicate_grid_fails_closed(self) -> None:
        bundle = _bundle()
        bundle["records"].pop()
        with self.assertRaisesRegex(bell.ContractError, "grid is incomplete"):
            bell.analyze_trace_bundle(bundle)

        bundle = _bundle()
        bundle["records"].append(copy.deepcopy(bundle["records"][0]))
        with self.assertRaisesRegex(bell.ContractError, "duplicate"):
            bell.analyze_trace_bundle(bundle)

    def test_epsilon_numeric_string_fails_closed(self) -> None:
        bundle = _bundle()
        bundle["records"][0]["epsilon"] = "0.1"
        with self.assertRaisesRegex(bell.ContractError, "JSON number"):
            bell.analyze_trace_bundle(bundle)

    def test_phase_interval_outside_frozen_rule_fails_closed(self) -> None:
        bundle = _bundle()
        bundle["records"][0]["phase_interval"][0] = 1.0
        with self.assertRaisesRegex(bell.ContractError, "phase intervals"):
            bell.analyze_trace_bundle(bundle)

    def test_phase_change_must_be_reproducible_from_retained_trace(self) -> None:
        bundle = _bundle()
        bundle["records"][0]["phase_change"][0] += 0.25
        with self.assertRaisesRegex(bell.ContractError, "retained mode trace"):
            bell.analyze_trace_bundle(bundle)

    def test_wrong_phase_propagation_sign_is_reported_as_scientific_failure(self) -> None:
        bundle = _bundle()
        record = bundle["records"][0]
        time = np.asarray(record["normalized_time"])
        right = np.conjugate(
            np.asarray(record["right_mode_real"])
            + 1.0j * np.asarray(record["right_mode_imag"])
        )
        left = np.conjugate(
            np.asarray(record["left_mode_real"])
            + 1.0j * np.asarray(record["left_mode_imag"])
        )
        record["right_mode_real"] = right.real.tolist()
        record["right_mode_imag"] = right.imag.tolist()
        record["left_mode_real"] = left.real.tolist()
        record["left_mode_imag"] = left.imag.tolist()
        interval, change = bell._fixed_interval_phase_trace(time, right)
        record["phase_interval"] = interval.tolist()
        record["phase_change"] = change.tolist()
        report = bell.analyze_trace_bundle(bundle)
        self.assertFalse(report["passed"])
        failed = report["records"][0]
        self.assertLess(failed["measured_phase_frequency_over_k0_ua"], 0.0)
        self.assertFalse(failed["phase_pass"])

    def test_nonqualifying_seed_fails_closed(self) -> None:
        bundle = _bundle()
        bundle["qualifying_seed"] = 23050091
        with self.assertRaisesRegex(bell.ContractError, "seed is not preregistered"):
            bell.analyze_trace_bundle(bundle)

    def test_wrong_polarization_is_reported_as_scientific_failure(self) -> None:
        bundle = _bundle()
        record = bundle["records"][0]
        record["left_mode_real"] = (
            2.0 * np.asarray(record["right_mode_real"])
        ).tolist()
        record["left_mode_imag"] = (
            2.0 * np.asarray(record["right_mode_imag"])
        ).tolist()
        report = bell.analyze_trace_bundle(bundle)
        self.assertFalse(report["passed"])
        failed = report["records"][0]
        self.assertFalse(failed["polarization_pass"])

    def test_raw_diagonal_mode_trace_extraction(self) -> None:
        dimension = 2
        epsilon = 0.4
        _, growth = bell.theoretical_dispersion(epsilon)
        datasets = _raw_datasets(dimension, epsilon)
        record = bell.extract_trace_record_from_datasets(
            dimension,
            epsilon,
            datasets,
            raw_provenance=bell.synthetic_contract_provenance(dimension, epsilon),
        )
        report = bell._analyze_record(record)
        self.assertTrue(report["passed"])
        self.assertFalse(report["section52_qualification_eligible"])
        self.assertAlmostEqual(report["measured_growth_rate_over_k0_ua"], growth)
        self.assertAlmostEqual(report["measured_phase_frequency_over_k0_ua"], epsilon)
        self.assertAlmostEqual(
            report["velocity_measured_growth_rate_over_k0_ua"], growth
        )
        self.assertAlmostEqual(
            report["velocity_measured_phase_frequency_over_k0_ua"], epsilon
        )
        self.assertTrue(report["velocity_magnetic_ratio_pass"])
        self.assertLess(report["velocity_magnetic_ratio_max_absolute_error"], 2.0e-12)
        self.assertAlmostEqual(
            report["paper_literal_measured_growth_rate_over_k0_ua"], growth
        )
        self.assertAlmostEqual(
            report["paper_literal_measured_phase_frequency_over_k0_ua"], epsilon
        )
        self.assertTrue(report["paper_literal_growth_pass"])
        self.assertTrue(report["paper_literal_phase_pass"])

        with tempfile.TemporaryDirectory() as directory:
            artifact_root = Path(directory).resolve()
            paths = []
            by_name = {}
            for index, dataset in enumerate(reversed(datasets)):
                path = artifact_root / f"snapshot-{index}.bin"
                path.write_bytes(f"snapshot-{index}\n".encode("ascii"))
                paths.append(path)
                by_name[str(path)] = dataset
            reread = bell.extract_trace_record_from_binary_files(
                dimension,
                epsilon,
                paths,
                "Q023-SOURCE-LOCAL-BASELINE-2D-EPSILON-0P4",
                reader=lambda name: by_name[name],
                artifact_root=artifact_root,
            )
            self.assertEqual(
                reread["raw_provenance"]["raw_artifacts"][0]["path"],
                "snapshot-0.bin",
            )
            self.assertTrue(
                bell._analyze_record(reread, artifact_root=artifact_root)["passed"]
            )
            with self.assertRaisesRegex(bell.ContractError, "explicit artifact root"):
                bell._analyze_record(reread)
            tampered = copy.deepcopy(reread)
            tampered["raw_provenance"]["raw_artifacts"][0]["sha256"] = "0" * 64
            with self.assertRaisesRegex(bell.ContractError, "artifact digest mismatch"):
                bell._analyze_record(tampered, artifact_root=artifact_root)
        self.assertEqual(reread["normalized_time"], record["normalized_time"])
        self.assertEqual(reread["right_mode_real"], record["right_mode_real"])
        self.assertEqual(reread["right_mode_imag"], record["right_mode_imag"])
        self.assertEqual(
            reread["raw_provenance"]["kind"], "source_local_materialized_variant"
        )
        self.assertEqual(len(reread["raw_provenance"]["raw_artifacts"]), len(datasets))

    def test_raw_geometry_outside_approved_variant_fails_closed(self) -> None:
        dimension = 2
        epsilon = 0.4
        datasets = _raw_datasets(dimension, epsilon)
        datasets[0]["x1v"] = 2.0 * np.asarray(datasets[0]["x1v"])
        with self.assertRaisesRegex(bell.ContractError, "geometry"):
            bell.extract_trace_record_from_datasets(
                dimension,
                epsilon,
                datasets,
                raw_provenance=bell.synthetic_contract_provenance(dimension, epsilon),
            )

    def test_dimension_and_retained_raw_geometry_numeric_aliases_fail_closed(
        self,
    ) -> None:
        for alias in (True, 1.0):
            with self.subTest(dimension_alias=alias):
                with self.assertRaisesRegex(bell.ContractError, "dimension"):
                    bell.synthetic_contract_provenance(alias, 0.4)
                with self.assertRaisesRegex(bell.ContractError, "dimension"):
                    bell._mode_basis(alias)
                provenance = bell.synthetic_contract_provenance(1, 0.4)
                with self.assertRaisesRegex(bell.ContractError, "dimension"):
                    bell._validate_raw_provenance(provenance, alias, 0.4)

        for key, index, alias in (
            ("nx", 0, 32.0),
            ("xmin", 0, False),
            ("extent", 0, 1),
        ):
            with self.subTest(geometry_key=key, geometry_alias=alias):
                provenance = bell.synthetic_contract_provenance(1, 0.4)
                provenance["raw_geometry"][key][index] = alias
                with self.assertRaisesRegex(bell.ContractError, "geometry"):
                    bell._validate_raw_provenance(provenance, 1, 0.4)

    def test_raw_float32_geometry_serialization_is_accepted(self) -> None:
        datasets = _raw_datasets(3, 0.4)
        for dataset in datasets:
            for axis in (1, 2, 3):
                dataset[f"x{axis}v"] = np.asarray(dataset[f"x{axis}v"]).astype(
                    np.float32
                )
        record = bell.extract_trace_record_from_datasets(
            3,
            0.4,
            datasets,
            raw_provenance=bell.synthetic_contract_provenance(3, 0.4),
        )
        self.assertTrue(bell._analyze_record(record)["passed"])

    def test_raw_combined_output_requires_velocity_fields(self) -> None:
        datasets = _raw_datasets(2, 0.4)
        datasets[0].pop("vely")
        with self.assertRaisesRegex(bell.ContractError, "velocity components"):
            bell.extract_trace_record_from_datasets(
                2,
                0.4,
                datasets,
                raw_provenance=bell.synthetic_contract_provenance(2, 0.4),
            )

    def test_velocity_magnetic_ratio_drift_is_a_scientific_failure(self) -> None:
        bundle = _bundle()
        record = bundle["records"][0]
        record["velocity_right_mode_real"] = (
            2.0 * np.asarray(record["velocity_right_mode_real"])
        ).tolist()
        record["velocity_right_mode_imag"] = (
            2.0 * np.asarray(record["velocity_right_mode_imag"])
        ).tolist()
        report = bell.analyze_trace_bundle(bundle)
        self.assertFalse(report["scientific_contract_pass"])
        failed = report["records"][0]
        self.assertFalse(failed["velocity_magnetic_ratio_pass"])

    def test_paper_literal_delta_u_y_phase_must_match_retained_sine_fit(self) -> None:
        bundle = _bundle()
        bundle["records"][0]["paper_delta_u_y_phase_change"][0] += 0.25
        with self.assertRaisesRegex(
            bell.ContractError, "paper-literal delta_u_y phase changes"
        ):
            bell.analyze_trace_bundle(bundle)

    def test_wrong_paper_literal_delta_u_y_propagation_sign_is_scientific_failure(
        self,
    ) -> None:
        bundle = _bundle()
        record = bundle["records"][0]
        time = np.asarray(record["normalized_time"])
        sine_fit = np.conjugate(
            np.asarray(record["paper_delta_u_y_sine_fit_real"])
            + 1.0j * np.asarray(record["paper_delta_u_y_sine_fit_imag"])
        )
        record["paper_delta_u_y_sine_fit_real"] = sine_fit.real.tolist()
        record["paper_delta_u_y_sine_fit_imag"] = sine_fit.imag.tolist()
        interval, change = bell._fixed_interval_phase_trace(time, sine_fit)
        record["paper_delta_u_y_phase_interval"] = interval.tolist()
        record["paper_delta_u_y_phase_change"] = change.tolist()
        report = bell.analyze_trace_bundle(bundle)
        self.assertFalse(report["scientific_contract_pass"])
        failed = report["records"][0]
        self.assertLess(
            failed["paper_literal_measured_phase_frequency_over_k0_ua"], 0.0
        )
        self.assertFalse(failed["paper_literal_phase_pass"])

    def test_paper_literal_volume_averaged_abs_delta_u_growth_is_required(self) -> None:
        bundle = _bundle()
        record = bundle["records"][0]
        record["paper_volume_averaged_abs_delta_u"] = (
            2.0 * np.asarray(record["paper_volume_averaged_abs_delta_u"])
        ).tolist()
        report = bell.analyze_trace_bundle(bundle)
        self.assertTrue(report["scientific_contract_pass"])
        self.assertTrue(report["records"][0]["paper_literal_growth_pass"])

        record["paper_volume_averaged_abs_delta_u"] = [1.0] * len(
            record["normalized_time"]
        )
        report = bell.analyze_trace_bundle(bundle)
        self.assertFalse(report["scientific_contract_pass"])
        self.assertFalse(report["records"][0]["paper_literal_growth_pass"])

    def test_binary_extraction_rejects_unapproved_variant(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            artifact_root = Path(directory).resolve()
            path = artifact_root / "snapshot.bin"
            path.write_bytes(b"snapshot\n")
            with self.assertRaisesRegex(bell.ContractError, "approved materialized"):
                bell.extract_trace_record_from_binary_files(
                    2,
                    0.4,
                    [path],
                    "Q023-SOURCE-LOCAL-UNREVIEWED",
                    reader=lambda name: {},
                    artifact_root=artifact_root,
                )

    def test_materialized_artifacts_require_normalized_authorized_root_paths(
        self,
    ) -> None:
        variant_id = "Q023-SOURCE-LOCAL-BASELINE-2D-EPSILON-0P4"
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory).resolve()
            artifact_root = directory_path / "authorized"
            artifact_root.mkdir()
            inside = artifact_root / "snapshot.bin"
            inside.write_bytes(b"inside\n")
            outside = directory_path / "outside.bin"
            outside.write_bytes(b"outside\n")

            provenance = bell._source_local_materialized_provenance(
                2,
                0.4,
                variant_id,
                [Path("snapshot.bin")],
                artifact_root=artifact_root,
            )
            self.assertEqual(provenance["raw_artifacts"][0]["path"], "snapshot.bin")
            bell._validate_raw_provenance(
                provenance, 2, 0.4, artifact_root=artifact_root
            )

            with self.assertRaisesRegex(bell.ContractError, "explicit artifact root"):
                bell._validate_raw_provenance(provenance, 2, 0.4)
            with self.assertRaisesRegex(bell.ContractError, "outside the authorized root"):
                bell._source_local_materialized_provenance(
                    2,
                    0.4,
                    variant_id,
                    [outside],
                    artifact_root=artifact_root,
                )
            with self.assertRaisesRegex(bell.ContractError, "root must be absolute"):
                bell._source_local_materialized_provenance(
                    2,
                    0.4,
                    variant_id,
                    [Path("snapshot.bin")],
                    artifact_root=Path("authorized"),
                )

            absolute = copy.deepcopy(provenance)
            absolute["raw_artifacts"][0]["path"] = str(inside)
            with self.assertRaisesRegex(bell.ContractError, "normalized root-relative"):
                bell._validate_raw_provenance(
                    absolute, 2, 0.4, artifact_root=artifact_root
                )

            tampered_deck = copy.deepcopy(provenance)
            tampered_deck["deck_sha256"] = "0" * 64
            with self.assertRaisesRegex(bell.ContractError, "deck digest mismatch"):
                bell._validate_raw_provenance(
                    tampered_deck, 2, 0.4, artifact_root=artifact_root
                )

            variant = bell._APPROVED_SOURCE_LOCAL_RAW_VARIANTS[variant_id]
            expected_deck_sha256 = variant["deck_sha256"]
            try:
                variant["deck_sha256"] = "0" * 64
                with self.assertRaisesRegex(bell.ContractError, "pinned value"):
                    bell._source_local_materialized_provenance(
                        2,
                        0.4,
                        variant_id,
                        [inside],
                        artifact_root=artifact_root,
                    )
            finally:
                variant["deck_sha256"] = expected_deck_sha256

    def test_readiness_bindings_distinguish_candidates_from_engineering_proxy(
        self,
    ) -> None:
        profiles = json.loads(PROFILES.read_text(encoding="utf-8"))
        drafts = json.loads(DRAFTS.read_text(encoding="utf-8"))
        inputs = {item["input_id"]: item for item in profiles["local_input_inventory"]}
        analyzers = {
            item["analyzer_id"]: item for item in profiles["local_analyzer_inventory"]
        }
        campaign = next(
            item for item in drafts["campaigns"]
            if item["campaign_id"] == bell.CAMPAIGN_ID
        )

        self.assertEqual(set(campaign["local_input_ids"]), PAPER_INPUT_IDS)
        self.assertNotIn("Q023-INPUT-BELL-PROXY", campaign["local_input_ids"])
        self.assertIn("Q023-INPUT-BELL-PROXY", inputs)
        self.assertEqual(
            campaign["local_analyzer_ids"],
            [
                "Q023-ANALYZER-PAPER-BELL-LINEAR-CANDIDATE",
                "Q023-MATERIALIZER-PAPER-BELL-LINEAR-VARIANTS",
            ],
        )
        analyzer = analyzers["Q023-ANALYZER-PAPER-BELL-LINEAR-CANDIDATE"]
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        sidecar_artifacts = {
            item["path"]: item["sha256"] for item in sidecar["source_local_artifacts"]
        }
        # The candidate analyzer digest remains consistent across both registries.
        self.assertEqual(sidecar_artifacts[analyzer["path"]], analyzer["sha256"])
        for input_id in PAPER_INPUT_IDS:
            candidate = inputs[input_id]
            self.assertEqual(
                _sha256(REPO_ROOT / candidate["path"]),
                candidate["sha256"],
            )

        bindings = campaign["candidate_bindings"]
        self.assertEqual(bindings["status"], "open_before_qualifying_run")
        self.assertIn("open_clean_candidate", bindings["git_commit"])
        self.assertIn("open_clean_frontier_executable",
                      bindings["executable_sha256"])
        self.assertIn("source_local", bindings["qualifying_input_checksums"])
        self.assertIn("source_local", bindings["production_analysis_checksums"])
        self.assertTrue(any("Section 5.2" in item
                            for item in campaign["claim_specific_exclusions"]))

    def test_sidecar_binds_approved_raw_variants_and_velocity_boundary(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        contract = sidecar["raw_trace_contract"]
        self.assertEqual(contract["signed_phase_frequency"],
                         "required_without_absolute_value")
        variants = {
            item["variant_id"]: item
            for item in contract["approved_source_local_materialized_variants"]
        }
        self.assertEqual(set(variants), set(bell._APPROVED_SOURCE_LOCAL_RAW_VARIANTS))
        for variant_id, expected in bell._APPROVED_SOURCE_LOCAL_RAW_VARIANTS.items():
            item = variants[variant_id]
            self.assertEqual(item["dimension"], expected["dimension"])
            self.assertEqual(item["epsilon"], expected["epsilon"])
            self.assertEqual(
                item["deck_path"],
                str(Path(expected["deck"]).relative_to(REPO_ROOT)),
            )
            self.assertEqual(item["deck_sha256"], expected["deck_sha256"])
            self.assertEqual(_sha256(Path(expected["deck"])), expected["deck_sha256"])
            self.assertEqual(
                item["raw_geometry"],
                bell._raw_geometry(expected["dimension"]),
            )
        boundary = contract["paper_velocity_observable_boundary"]
        self.assertFalse(boundary["section52_qualification_eligible"])
        self.assertEqual(
            boundary["current_retained_output"],
            "raw_mhd_w_bcc_combined_magnetic_and_fluid_velocity_modes",
        )
        self.assertEqual(
            boundary["paper_benchmark_observable"],
            "paper_literal_delta_u_y_spatial_sine_fit_phase_and_volume_averaged_"
            "absolute_delta_u_growth_estimator_frozen_source_local",
        )
        self.assertEqual(
            boundary["paper_literal_estimator_qualification_effect"],
            "source_local_contract_only_not_section52_qualification",
        )
        self.assertIn(
            "registered_frontier_execution",
            boundary["required_before_qualification"],
        )


if __name__ == "__main__":
    unittest.main()
