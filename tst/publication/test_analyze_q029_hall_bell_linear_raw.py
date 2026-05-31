#!/usr/bin/env python3
"""Focused regressions for the Q-029 source-local raw extraction scaffold."""

from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from tst.publication import analyze_q029_hall_bell_linear_raw as raw


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q029_hall_bell_linear_raw_extractor_source_local_2026-05-30.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Q029HallBellLinearRawTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        preparation = raw._load_preparation_record()
        cls.artifact_root = Path(
            preparation["source_local_runtime_smoke"]["artifact_root"]
        )
        cls.bundle = raw.extract_raw_trace_bundle(cls.artifact_root)

    def test_extracts_approved_combined_projected_traces_as_nonqualifying(self) -> None:
        bundle = self.bundle
        self.assertEqual(bundle["campaign_id"], raw.CAMPAIGN_ID)
        self.assertEqual(bundle["artifact_role"], raw.ARTIFACT_ROLE)
        self.assertEqual(bundle["qualification_effect"], raw.QUALIFICATION_EFFECT)
        self.assertFalse(bundle["qualifying_evidence"])
        self.assertFalse(bundle["hall_bell_qualification"])
        self.assertFalse(bundle["linear_hall_bell_qualification"])
        self.assertFalse(bundle["nonlinear_hall_bell_qualification"])
        self.assertEqual(bundle["projection_contract"], raw.PROJECTION_CONTRACT)
        self.assertEqual(bundle["provenance"]["artifact_writable_entries"], 0)
        self.assertEqual(
            (
                bundle["provenance"]["artifact_root_device"],
                bundle["provenance"]["artifact_root_inode"],
            ),
            raw.APPROVED_ARTIFACT_ROOT_IDENTITY,
        )
        self.assertEqual(
            bundle["provenance"]["binary_reader_sha256"],
            raw.BINARY_READER_SHA256,
        )
        self.assertEqual(
            [trace["trace_id"] for trace in bundle["traces"]],
            [source["trace_id"] for source in raw._APPROVED_TRACE_SOURCES],
        )
        self.assertEqual(
            [len(trace["raw_artifacts"]) for trace in bundle["traces"]],
            [1, 1, 1, 2],
        )
        for trace in bundle["traces"]:
            size = len(trace["raw_artifacts"])
            self.assertEqual(len(trace["time"]), size)
            self.assertEqual(len(trace["num_cycles"]), size)
            self.assertEqual(len(trace["magnetic_right_mode_real"]), size)
            self.assertEqual(len(trace["magnetic_left_mode_imag"]), size)
            self.assertEqual(len(trace["velocity_right_mode_real"]), size)
            self.assertEqual(len(trace["velocity_left_mode_imag"]), size)

    def test_exact_extracted_bundle_replays(self) -> None:
        reread = json.loads(json.dumps(self.bundle, allow_nan=False))
        self.assertEqual(
            raw.validate_raw_trace_bundle(reread, self.artifact_root),
            reread,
        )

    def test_raw_dataset_schema_geometry_and_velocity_drift_fail_closed(self) -> None:
        reader = raw._binary_reader()
        source = raw._APPROVED_TRACE_SOURCES[1]
        artifact = source["raw_artifacts"][0]
        dataset = reader(str(self.artifact_root / artifact["path"]))

        unexpected = copy.deepcopy(dataset)
        unexpected["not_admitted"] = np.zeros(1)
        with self.assertRaisesRegex(raw.ContractError, "dataset keys"):
            raw._project_dataset(unexpected, source["dimension"])

        geometry = copy.deepcopy(dataset)
        geometry["x1f"][0] = np.float32(0.25)
        with self.assertRaisesRegex(raw.ContractError, "geometry drift"):
            raw._project_dataset(geometry, source["dimension"])

        velocity = copy.deepcopy(dataset)
        velocity.pop("vely")
        with self.assertRaisesRegex(raw.ContractError, "dataset keys"):
            raw._project_dataset(velocity, source["dimension"])

    def test_bundle_schema_provenance_and_projected_value_drift_fail_closed(self) -> None:
        unexpected = copy.deepcopy(self.bundle)
        unexpected["not_admitted"] = True
        with self.assertRaisesRegex(raw.ContractError, "bundle keys"):
            raw.validate_raw_trace_bundle(unexpected, self.artifact_root)

        provenance = copy.deepcopy(self.bundle)
        provenance["provenance"]["artifact_inventory_sha256"] = "0" * 64
        with self.assertRaisesRegex(raw.ContractError, "raw provenance drift"):
            raw.validate_raw_trace_bundle(provenance, self.artifact_root)

        artifact = copy.deepcopy(self.bundle)
        artifact["traces"][0]["raw_artifacts"][0]["sha256"] = "0" * 64
        with self.assertRaisesRegex(raw.ContractError, "trace provenance drift"):
            raw.validate_raw_trace_bundle(artifact, self.artifact_root)

        projected = copy.deepcopy(self.bundle)
        projected["traces"][0]["velocity_right_mode_real"][0] += 1.0
        with self.assertRaisesRegex(raw.ContractError, "trace value drift"):
            raw.validate_raw_trace_bundle(projected, self.artifact_root)

    def test_dimension_geometry_and_provenance_numeric_aliases_fail_closed(
        self,
    ) -> None:
        for alias in (True, 1.0):
            with self.subTest(dimension_alias=alias):
                with self.assertRaisesRegex(raw.ContractError, "dimension"):
                    raw._raw_geometry(alias)
                with self.assertRaisesRegex(raw.ContractError, "dimension"):
                    raw._mode_basis(alias)
                bundle = copy.deepcopy(self.bundle)
                bundle["traces"][0]["dimension"] = alias
                with self.assertRaisesRegex(raw.ContractError, "dimension"):
                    raw.validate_raw_trace_bundle(bundle, self.artifact_root)

        for key, index, alias in (
            ("nx", 0, 32.0),
            ("xmin", 0, False),
            ("extent", 0, 1),
        ):
            with self.subTest(geometry_key=key, geometry_alias=alias):
                bundle = copy.deepcopy(self.bundle)
                bundle["traces"][0]["raw_geometry"][key][index] = alias
                with self.assertRaisesRegex(raw.ContractError, "geometry schema drift"):
                    raw.validate_raw_trace_bundle(bundle, self.artifact_root)

        for key, alias in (
            ("artifact_file_count", float(
                self.bundle["provenance"]["artifact_file_count"]
            )),
            ("artifact_writable_entries", False),
        ):
            with self.subTest(provenance_key=key, provenance_alias=alias):
                bundle = copy.deepcopy(self.bundle)
                bundle["provenance"][key] = alias
                with self.assertRaisesRegex(raw.ContractError, "identity schema drift"):
                    raw.validate_raw_trace_bundle(bundle, self.artifact_root)

    def test_preparation_record_digest_pin_rejects_schema_drift(self) -> None:
        preparation = json.loads(raw.PREPARATION_RECORD.read_text(encoding="utf-8"))
        preparation["not_admitted"] = True
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "preparation.json"
            path.write_text(json.dumps(preparation), encoding="utf-8")
            with self.assertRaisesRegex(raw.ContractError, "digest mismatch"):
                raw._load_preparation_record(path)

    def test_preparation_dimensions_and_retention_aliases_fail_closed(self) -> None:
        smoke = copy.deepcopy(
            raw._load_preparation_record()["source_local_runtime_smoke"]
        )
        smoke["cycle_zero_initializations"][0]["dimension"] = True
        with self.assertRaisesRegex(raw.ContractError, "dimension"):
            raw._validate_preparation_artifact_bindings(smoke)

        smoke = copy.deepcopy(
            raw._load_preparation_record()["source_local_runtime_smoke"]
        )
        smoke["restart_continuation"]["dimension"] = 2.0
        with self.assertRaisesRegex(raw.ContractError, "dimension"):
            raw._validate_preparation_artifact_bindings(smoke)

        for alias in (False, 0):
            with self.subTest(restart_parity_alias=alias):
                smoke = copy.deepcopy(
                    raw._load_preparation_record()["source_local_runtime_smoke"]
                )
                smoke["restart_continuation"]["max_absolute_field_difference"] = (
                    alias
                )
                with self.assertRaisesRegex(raw.ContractError, "restart-parity"):
                    raw._validate_preparation_artifact_bindings(smoke)

        preparation = json.loads(raw.PREPARATION_RECORD.read_text(encoding="utf-8"))
        preparation["source_local_runtime_smoke"]["retention"][
            "writable_entries"
        ] = False
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "preparation.json"
            path.write_text(json.dumps(preparation), encoding="utf-8")
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            with patch.object(raw, "PREPARATION_RECORD_SHA256", digest):
                with self.assertRaisesRegex(raw.ContractError, "retention provenance"):
                    raw._load_preparation_record(path)

    def test_artifact_root_must_be_the_exact_approved_absolute_root(self) -> None:
        smoke = raw._load_preparation_record()["source_local_runtime_smoke"]
        with self.assertRaisesRegex(raw.ContractError, "must be absolute"):
            raw._authorized_artifact_root(Path("relative"), smoke)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(raw.ContractError, "not the approved"):
                raw._authorized_artifact_root(Path(directory), smoke)
        for key, alias in (
            ("file_count", float(smoke["retention"]["file_count"])),
            ("writable_entries", False),
        ):
            with self.subTest(retention_key=key, retention_alias=alias):
                aliased = copy.deepcopy(smoke)
                aliased["retention"][key] = alias
                with self.assertRaisesRegex(raw.ContractError, "provenance drift"):
                    raw._authorized_artifact_root(self.artifact_root, aliased)

    def test_production_api_rejects_reader_injection_and_inventory_special_files(
        self,
    ) -> None:
        with self.assertRaises(TypeError):
            raw.extract_raw_trace_bundle(self.artifact_root, reader=lambda _: {})
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            os.mkfifo(root / "not-a-regular-file")
            with self.assertRaisesRegex(raw.ContractError, "special file"):
                raw._inventory_sha256(root)

    def test_inventory_counts_each_writable_directory_once(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            nested = root / "nested"
            nested.mkdir()
            payload = nested / "payload"
            payload.write_text("payload", encoding="utf-8")
            payload.chmod(0o444)
            file_count, writable_entries, _ = raw._inventory_sha256(root)
            self.assertEqual(file_count, 1)
            self.assertEqual(writable_entries, 2)

    def test_authorized_root_identity_and_inventory_buffers_fail_closed(self) -> None:
        state = os.stat(self.artifact_root, follow_symlinks=False)
        self.assertEqual(
            (state.st_dev, state.st_ino),
            raw.APPROVED_ARTIFACT_ROOT_IDENTITY,
        )
        with tempfile.TemporaryDirectory() as directory:
            root_fd = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
            try:
                with self.assertRaisesRegex(raw.ContractError, "identity drift"):
                    raw._revalidate_authorized_artifact_root(Path(directory), root_fd)
            finally:
                os.close(root_fd)
        retained = b"retained inventory bytes"
        digest = hashlib.sha256(retained).hexdigest()
        self.assertEqual(
            raw._normalized_artifact_bytes({"selected.bin": retained}, "selected.bin", digest),
            retained,
        )
        with self.assertRaisesRegex(raw.ContractError, "not an approved file"):
            raw._normalized_artifact_bytes({}, "selected.bin", digest)

    def test_decoder_bridge_passes_a_write_sealed_memfd(self) -> None:
        payload = b"sealed raw bytes"
        digest = hashlib.sha256(payload).hexdigest()
        observed = {}

        def reader(path: str) -> dict[str, str]:
            observed["path"] = path
            fd = os.open(path, os.O_WRONLY)
            try:
                with self.assertRaises(OSError):
                    os.write(fd, b"replacement")
            finally:
                os.close(fd)
            return {"status": "sealed"}

        self.assertEqual(
            raw._read_verified_dataset(payload, digest, reader),
            {"status": "sealed"},
        )
        self.assertRegex(observed["path"], r"^/proc/self/fd/[0-9]+$")

    def test_verified_decoder_source_and_runtime_versions_are_enforced(self) -> None:
        source = (
            REPO_ROOT / "tst/publication/analyze_q029_hall_bell_linear_raw.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("spec_from_file_location", source)
        self.assertNotIn("_normalized_artifact_file", source)
        self.assertIn('exec(compile(raw, str(BINARY_READER), "exec"), namespace)', source)
        self.assertIn(
            "root, root_fd, file_bytes = _authorized_artifact_root",
            source,
        )
        self.assertIn('reader(f"/proc/self/fd/{fd}")', source)
        self.assertIn('os.memfd_create("q029-raw-decode", os.MFD_ALLOW_SEALING)', source)
        self.assertNotIn("TemporaryFile", source)
        smoke = raw._load_preparation_record()["source_local_runtime_smoke"]
        with patch.object(raw.platform, "python_version", return_value="0.0.0"):
            with self.assertRaisesRegex(raw.ContractError, "runtime binding"):
                raw._provenance(self.artifact_root, smoke)
        with patch.object(raw.np, "__version__", "0.0.0"):
            with self.assertRaisesRegex(raw.ContractError, "runtime binding"):
                raw._provenance(self.artifact_root, smoke)

    def test_analyzer_contains_no_physical_oracle_grid_mapping_or_tolerance(self) -> None:
        source = (
            REPO_ROOT / "tst/publication/analyze_q029_hall_bell_linear_raw.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("theoretical_dispersion", source)
        self.assertNotIn("EPSILON_VALUES", source)
        self.assertNotIn("CHI_H_VALUES", source)
        self.assertNotIn("reference_mapping", source)
        self.assertNotIn("TOLERANCE", source)

    def test_readiness_sidecar_binds_only_new_q029_raw_extractor_files(self) -> None:
        readiness = json.loads(READINESS.read_text(encoding="utf-8"))
        self.assertEqual(readiness["gate"], "Q-029")
        self.assertEqual(readiness["campaign_id"], raw.CAMPAIGN_ID)
        self.assertFalse(readiness["claim_closure"])
        self.assertEqual(
            readiness["qualification_effect"], raw.QUALIFICATION_EFFECT
        )
        self.assertEqual(
            readiness["predecessor_binding"],
            {
                "path": str(raw.PREPARATION_RECORD.relative_to(REPO_ROOT)),
                "sha256": raw.PREPARATION_RECORD_SHA256,
            },
        )
        expected_paths = {
            "tst/publication/analyze_q029_hall_bell_linear_raw.py",
            "tst/publication/test_analyze_q029_hall_bell_linear_raw.py",
            "vis/python/bin_convert_new.py",
        }
        self.assertEqual(set(readiness["artifact_bindings"]), expected_paths)
        for relative, digest in readiness["artifact_bindings"].items():
            self.assertEqual(_sha256(REPO_ROOT / relative), digest)
        boundary = readiness["nonqualification_boundary"]
        self.assertFalse(boundary["qualifying_evidence"])
        self.assertFalse(boundary["hall_bell_qualification"])
        self.assertEqual(boundary["numeric_tolerance"], "intentionally_absent")
        self.assertEqual(boundary["physical_bai_dispersion_oracle"], "absent")
        self.assertEqual(boundary["physical_coefficient_grid"], "absent")
        self.assertEqual(boundary["extracted_reference_mapping"], "absent")


if __name__ == "__main__":
    unittest.main()
