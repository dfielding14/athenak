#!/usr/bin/env python3
"""Focused tests for the Q-023 Section 5.2 Bell source-local candidates."""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
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
    phase_interval, phase_change = bell._fixed_interval_phase_trace(time, right)
    return {
        "dimension": dimension,
        "epsilon": epsilon,
        "normalized_time": time.tolist(),
        "right_mode_real": right.real.tolist(),
        "right_mode_imag": right.imag.tolist(),
        "left_mode_real": left.real.tolist(),
        "left_mode_imag": left.imag.tolist(),
        "phase_interval": phase_interval.tolist(),
        "phase_change": phase_change.tolist(),
    }


def _bundle() -> dict[str, object]:
    return {
        "schema_version": 1,
        "campaign_id": bell.CAMPAIGN_ID,
        "qualifying_seed": bell.QUALIFYING_SEEDS[0],
        "records": [
            _record(dimension, epsilon)
            for dimension in bell.DIMENSIONS
            for epsilon in bell.EPSILON_VALUES
        ],
    }


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
        self.assertTrue(report["passed"])
        self.assertEqual(report["record_count"], 15)
        self.assertIn("source_local_candidate_analysis_only",
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
        self.assertTrue(report["passed"])

    def test_incomplete_or_duplicate_grid_fails_closed(self) -> None:
        bundle = _bundle()
        bundle["records"].pop()
        with self.assertRaisesRegex(bell.ContractError, "grid is incomplete"):
            bell.analyze_trace_bundle(bundle)

        bundle = _bundle()
        bundle["records"].append(copy.deepcopy(bundle["records"][0]))
        with self.assertRaisesRegex(bell.ContractError, "duplicate"):
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
        dimension = 3
        epsilon = 0.4
        _, growth = bell.theoretical_dispersion(epsilon)
        parallel, transverse_a, transverse_b = bell._mode_basis(dimension)
        extents = bell._EXPECTED_GEOMETRY[dimension]["extent"]
        coordinates = [
            (np.arange(count, dtype=float) + 0.5) * extent / count
            for count, extent in zip((16, 8, 4), extents)
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
            datasets.append(
                {
                    "Time": normalized_time / (bell.K0*bell.U_A),
                    "x1v": coordinates[0],
                    "x2v": coordinates[1],
                    "x3v": coordinates[2],
                    "bcc1": magnetic[0],
                    "bcc2": magnetic[1],
                    "bcc3": magnetic[2],
                }
            )
        record = bell.extract_trace_record_from_datasets(dimension, epsilon, datasets)
        report = bell._analyze_record(record)
        self.assertTrue(report["passed"])
        self.assertAlmostEqual(report["measured_growth_rate_over_k0_ua"], growth)
        self.assertAlmostEqual(report["measured_phase_frequency_over_k0_ua"], epsilon)

        by_name = {f"snapshot-{index}.bin": dataset
                   for index, dataset in enumerate(reversed(datasets))}
        reread = bell.extract_trace_record_from_binary_files(
            dimension,
            epsilon,
            [Path(name) for name in by_name],
            reader=lambda name: by_name[name],
        )
        self.assertEqual(reread, record)

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
            ["Q023-ANALYZER-PAPER-BELL-LINEAR-CANDIDATE"],
        )
        analyzer = analyzers["Q023-ANALYZER-PAPER-BELL-LINEAR-CANDIDATE"]
        self.assertEqual(_sha256(REPO_ROOT / analyzer["path"]), analyzer["sha256"])
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


if __name__ == "__main__":
    unittest.main()
