#!/usr/bin/env python3
"""Focused tests for the bounded source-local Q-033 runtime successor."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest import mock

import numpy as np

from tst.publication import analyze_q033_crpai_transport_runtime_local as q033


REPO_ROOT = Path(__file__).resolve().parents[2]
SIDECAR = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q033_crpai_transport_runtime_local_2026-05-30.json"
)
BLOCKED_DECK = (
    REPO_ROOT / "inputs/tests/pic_q033_crpai_transport_calibration_candidate.athinput"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _seeded_fields() -> tuple[np.ndarray, np.ndarray]:
    descriptor = q033.build_static_descriptor()
    x1 = (np.arange(32, dtype=np.float64) + 0.5) / 32.0
    by = np.zeros_like(x1)
    bz = np.zeros_like(x1)
    for item in descriptor["wave_seed_contract"]["modes"]:
        phase = 2.0 * math.pi * item["mode"] * x1 + item["phase"]
        by += item["amplitude"] * np.cos(phase)
        handedness = 1.0 if item["handedness"] == "positive" else -1.0
        bz += handedness * item["amplitude"] * np.sin(phase)
    return np.tile(by, (1, 4, 1)), np.tile(bz, (1, 4, 1))


def _particle_payload(count: int = 8) -> dict[str, np.ndarray]:
    velocity = np.array([
        [1.0, 0.2, 0.1],
        [-1.0, -0.2, -0.1],
        [1.2, -0.1, 0.3],
        [-1.2, 0.1, -0.3],
        [0.8, 0.4, -0.2],
        [-0.8, -0.4, 0.2],
        [1.4, -0.3, -0.4],
        [-1.4, 0.3, 0.4],
    ])
    return {
        "tags": np.arange(count, dtype=np.int64),
        "velocity": velocity[:count],
        "macro_weight": np.ones(count),
        "deltaf_f0": np.full(count, 0.5),
        "deltaf_weight": np.zeros(count),
    }


class Q033CRPAITransportRuntimeLocalTests(unittest.TestCase):
    def test_additive_registration_preserves_reserved_blocked_open_name(self) -> None:
        registration = q033.validate_registration()
        self.assertEqual(registration["pgen_name"], q033.PGEN_NAME)
        self.assertTrue(registration["fresh_dispatch"])
        self.assertTrue(registration["restart_dispatch"])
        self.assertFalse(registration["reserved_blocked_pgen_registered"])
        dispatch = q033.DISPATCH.read_text(encoding="utf-8")
        self.assertNotIn(q033.RESERVED_BLOCKED_PGEN_NAME, dispatch)
        self.assertIn(
            "pgen_name = q033_crpai_transport_calibration_open",
            BLOCKED_DECK.read_text(encoding="utf-8"),
        )

    def test_deck_is_bounded_serial_thin_runtime_diagnostic_only(self) -> None:
        deck = q033.validate_deck()
        self.assertEqual(deck["campaign_id"], q033.CAMPAIGN_ID)
        self.assertEqual(deck["qualification_effect"], "none")
        self.assertFalse(deck["qualifying_evidence"])
        self.assertFalse(deck["physical_calibration_claimed"])
        self.assertEqual(deck["frontier_authorization"], "not_bound")
        self.assertEqual(deck["q022_closure"], "not_claimed")
        self.assertEqual(deck["external_review"], "not_claimed")

    def test_static_descriptor_freezes_seeded_local_formulas_without_claims(self) -> None:
        descriptor = q033.build_static_descriptor()
        self.assertEqual(descriptor, q033.build_static_descriptor())
        self.assertFalse(descriptor["qualifying_evidence"])
        self.assertFalse(descriptor["physical_calibration_claimed"])
        self.assertFalse(descriptor["q022_closure_claimed"])
        self.assertFalse(descriptor["external_review_claimed"])
        self.assertFalse(descriptor["frontier_qualification_claimed"])
        self.assertFalse(descriptor["launch_authorization_claimed"])
        modes = descriptor["wave_seed_contract"]["modes"]
        self.assertEqual([item["mode"] for item in modes], [1, 2, 3])
        self.assertEqual(
            [item["handedness"] for item in modes],
            ["positive", "negative", "positive"],
        )
        self.assertEqual(
            [item["amplitude"] for item in modes],
            [1.0e-4, 2.5e-5, 1.0e-4 / 9.0],
        )

    def test_cycle_local_diagnostics_recompute_wave_and_antipodal_payload(self) -> None:
        by, bz = _seeded_fields()
        particle = _particle_payload()
        report = q033.analyze_cycle_local_arrays(
            by,
            bz,
            particle["tags"],
            particle["velocity"],
            particle["macro_weight"],
            particle["deltaf_f0"],
            particle["deltaf_weight"],
        )
        self.assertFalse(report["qualifying_evidence"])
        self.assertFalse(report["physical_calibration_claimed"])
        self.assertEqual(len(report["wave_spectrum"]), 3)
        self.assertEqual(report["particle_payload"]["count"], 8)
        self.assertEqual(report["particle_payload"]["adjacent_antipodal_pair_count"], 4)
        self.assertEqual(report["particle_payload"]["maximum_velocity_pair_residual"], 0.0)
        self.assertEqual(
            report["interpretation"],
            "bounded_cycle_local_runtime_diagnostics_without_reference_thresholds",
        )

    def test_diagnostic_payload_fails_closed_on_shape_tag_and_nonfinite_drift(self) -> None:
        by, bz = _seeded_fields()
        particle = _particle_payload()
        cases = [
            ("thin-carrier", by[..., :-1], particle["tags"], particle["velocity"]),
            ("adjacent pairs", by, particle["tags"][:-1], particle["velocity"][:-1]),
            ("serial tag", by, np.array([0, 1, 2, 3, 4, 5, 6, 6]), particle["velocity"]),
        ]
        for expected, field, tags, velocity in cases:
            with self.subTest(expected=expected):
                with self.assertRaisesRegex(q033.ContractError, expected):
                    q033.analyze_cycle_local_arrays(
                        field,
                        bz,
                        tags,
                        velocity,
                        particle["macro_weight"][: tags.size],
                        particle["deltaf_f0"][: tags.size],
                        particle["deltaf_weight"][: tags.size],
                    )
        velocity = particle["velocity"].copy()
        velocity[0, 0] = np.nan
        with self.assertRaisesRegex(q033.ContractError, "finite"):
            q033.analyze_cycle_local_arrays(
                by,
                bz,
                particle["tags"],
                velocity,
                particle["macro_weight"],
                particle["deltaf_f0"],
                particle["deltaf_weight"],
            )

    def test_immutable_artifact_extractor_hashes_inputs_and_emits_nonclaims(self) -> None:
        by, bz = _seeded_fields()
        particle = _particle_payload()
        vtk = SimpleNamespace(
            scalars={
                "ptag": particle["tags"],
                "macro_weight": particle["macro_weight"],
                "deltaf_f0": particle["deltaf_f0"],
                "deltaf_weight": particle["deltaf_weight"],
            },
            vectors={"vel": particle["velocity"]},
        )
        with tempfile.TemporaryDirectory() as directory:
            mhd_path = Path(directory) / "snapshot.bin"
            vtk_path = Path(directory) / "snapshot.part.vtk"
            mhd_path.write_bytes(b"mhd")
            vtk_path.write_bytes(b"particles")
            with mock.patch.object(
                q033, "_read_mhd_bcc", return_value={"bcc2": by, "bcc3": bz, "Time": 0.0}
            ), mock.patch.object(q033, "_read_particle_vtk", return_value=vtk):
                report = q033.extract_runtime_artifacts(mhd_path, vtk_path)
        self.assertFalse(report["qualifying_evidence"])
        self.assertEqual(report["immutable_artifacts"]["mhd_bcc"]["sha256"],
                         hashlib.sha256(b"mhd").hexdigest())
        self.assertEqual(report["immutable_artifacts"]["prtcl_all"]["sha256"],
                         hashlib.sha256(b"particles").hexdigest())

    def test_strict_deck_parser_rejects_duplicate_parameter(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.athinput"
            path.write_text("<mesh>\nnx1 = 32\nnx1 = 32\n", encoding="utf-8")
            with self.assertRaisesRegex(q033.ContractError, "duplicate parameter"):
                q033.parse_athinput(path)

    def test_readiness_sidecar_binds_only_new_q033_runtime_local_artifacts(self) -> None:
        sidecar = json.loads(SIDECAR.read_text(encoding="utf-8"))
        self.assertEqual(sidecar["gate"], "Q-033")
        self.assertEqual(sidecar["qualification_effect"], "none")
        self.assertFalse(sidecar["claim_closure"])
        self.assertEqual(sidecar["frontier_authorization"], "not_bound")
        expected_paths = {
            "src/pgen/tests/q033_crpai_transport_runtime_local.cpp",
            "inputs/tests/pic_q033_crpai_transport_runtime_local.athinput",
            "tst/publication/analyze_q033_crpai_transport_runtime_local.py",
            "tst/publication/test_analyze_q033_crpai_transport_runtime_local.py",
        }
        self.assertEqual(set(sidecar["artifact_bindings"]), expected_paths)
        for relative, expected_sha256 in sidecar["artifact_bindings"].items():
            self.assertEqual(_sha256(REPO_ROOT / relative), expected_sha256)
        self.assertFalse(sidecar["reserved_candidate"]["registered"])
        self.assertEqual(
            sidecar["reserved_candidate"]["name"],
            q033.RESERVED_BLOCKED_PGEN_NAME,
        )
        for nonclaim in (
            "physical CRPAI transport calibration",
            "Q-022 closure",
            "external-review closure",
            "Frontier qualification or launch authorization",
        ):
            self.assertIn(nonclaim, sidecar["explicitly_not_claimed"])


if __name__ == "__main__":
    unittest.main()
