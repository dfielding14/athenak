#!/usr/bin/env python3
"""Focused tests for the human Q-011 Section 5.4 pressure-selection boundary."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import unittest
from unittest import mock

PUBLICATION_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PUBLICATION_DIR / "frontier_control_plane"))

from . import publish_q011_section54_pressure_pilot_bundle as publisher
from . import q011_section54_pressure_selection as selection
from . import test_publish_q011_section54_pressure_pilot_bundle as publication_fixture


def _receipt(publication: dict[str, str], retained: dict[str, object]) -> dict[str, object]:
    descriptor_by_case = {
        item["case_id"]: item["descriptor_sha256"] for item in retained["raw_cases"]
    }
    return {
        "schema_version": 1,
        "record_type": selection.RECORD_TYPE,
        "selection_method": selection.SELECTION_METHOD,
        "published_pressure_pilot_receipt": {
            "path": publication["receipt_path"],
            "sha256": publication["receipt_sha256"],
        },
        "pilot_bundle_manifest_sha256": publication["manifest_sha256"],
        "aggregate_pilot_analysis_sha256": publication["analysis_result_sha256"],
        "case_descriptors": [
            {
                "case_id": case_id,
                "problem_ps_p0": problem_ps_p0,
                "descriptor_sha256": descriptor_by_case[case_id],
            }
            for case_id, problem_ps_p0 in selection.REGISTERED_CASES
        ],
        "selected_case": {
            "case_id": "ps_p0_0p10",
            "problem_ps_p0": 0.1,
        },
        "reviewer_identity": "Test-only Human Reviewer",
        "reviewed_utc": "2026-06-02T12:34:56Z",
        "rationale": "Test-only human rationale for selecting one registered pilot.",
    }


class Q011Section54PressureSelectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._fixture = publication_fixture._verified_raw_cases()
        cls.base, roots, digests = cls._fixture.__enter__()
        cls.bundle = cls.base / "bundle"
        cls.receipt_path = cls.base / "publication-receipt.json"
        cls.analysis_path = cls.base / "aggregate-analysis.json"
        cls.publication = publisher.publish_pressure_pilot_bundle(
            cls.bundle,
            receipt_path=cls.receipt_path,
            analysis_result_path=cls.analysis_path,
            case_artifact_dirs=roots,
            case_descriptor_sha256=digests,
            authorized_pic_root=cls.base.parent,
        )
        retained = json.loads(cls.receipt_path.read_text(encoding="utf-8"))
        cls.receipt_template = _receipt(cls.publication, retained)

    @classmethod
    def tearDownClass(cls) -> None:
        publication_fixture._make_writable(cls.bundle)
        for path in (cls.receipt_path, cls.analysis_path):
            path.chmod(0o600)
        cls._fixture.__exit__(None, None, None)

    def _receipt(self) -> dict[str, object]:
        return copy.deepcopy(self.receipt_template)

    def _validate(self, receipt: dict[str, object]) -> dict[str, object]:
        return selection.validate_pressure_selection_receipt(
            receipt,
            authorized_pic_root=self.base.parent,
        )

    def test_accepts_one_human_selection_from_verified_immutable_publication(self) -> None:
        receipt = self._receipt()
        with mock.patch.object(
            selection.pressure_pilot_publisher,
            "verify_published_pressure_pilot_receipt",
            wraps=publisher.verify_published_pressure_pilot_receipt,
        ) as verifier:
            self.assertEqual(self._validate(receipt), receipt)
        verifier.assert_called_once_with(
            str(self.receipt_path),
            authorized_pic_root=self.base.parent,
        )

    def test_rejects_self_asserted_arbitrary_publication_and_evidence_hashes(self) -> None:
        cases = []
        publication_receipt = self._receipt()
        publication_receipt["published_pressure_pilot_receipt"]["sha256"] = "0" * 64
        cases.append(publication_receipt)
        manifest = self._receipt()
        manifest["pilot_bundle_manifest_sha256"] = "1" * 64
        cases.append(manifest)
        analysis = self._receipt()
        analysis["aggregate_pilot_analysis_sha256"] = "2" * 64
        cases.append(analysis)
        descriptor = self._receipt()
        descriptor["case_descriptors"][0]["descriptor_sha256"] = "3" * 64
        cases.append(descriptor)
        for receipt in cases:
            with self.subTest(receipt=receipt), self.assertRaises(
                selection.PressureSelectionReceiptError
            ):
                self._validate(receipt)

    def test_requires_descriptor_stable_readonly_published_receipt(self) -> None:
        copied = self.base / "writable-self-asserted-receipt.json"
        copied.write_bytes(self.receipt_path.read_bytes())
        receipt = self._receipt()
        receipt["published_pressure_pilot_receipt"] = {
            "path": str(copied),
            "sha256": publisher._sha256(copied.read_bytes()),
        }
        try:
            with self.assertRaisesRegex(
                selection.PressureSelectionReceiptError,
                "expected read-only regular file",
            ):
                self._validate(receipt)
        finally:
            copied.unlink()

    def test_rejects_case_set_and_registered_problem_ps_p0_drift(self) -> None:
        missing = self._receipt()
        missing["case_descriptors"].pop()
        reordered = self._receipt()
        reordered["case_descriptors"][0], reordered["case_descriptors"][1] = (
            reordered["case_descriptors"][1],
            reordered["case_descriptors"][0],
        )
        changed = self._receipt()
        changed["case_descriptors"][2]["problem_ps_p0"] = 0.11
        integer_alias = self._receipt()
        integer_alias["case_descriptors"][0]["problem_ps_p0"] = 1
        repeated = self._receipt()
        repeated["case_descriptors"][3]["descriptor_sha256"] = repeated[
            "case_descriptors"
        ][0]["descriptor_sha256"]
        for receipt in (missing, reordered, changed, integer_alias, repeated):
            with self.subTest(receipt=receipt), self.assertRaises(
                selection.PressureSelectionReceiptError
            ):
                self._validate(receipt)

    def test_rejects_unregistered_mismatched_or_automated_selection(self) -> None:
        unknown = self._receipt()
        unknown["selected_case"]["case_id"] = "ps_p0_unknown"
        mismatched = self._receipt()
        mismatched["selected_case"]["problem_ps_p0"] = 0.2
        multiple = self._receipt()
        multiple["selected_case"] = [
            {"case_id": "ps_p0_0p10", "problem_ps_p0": 0.1},
            {"case_id": "ps_p0_0p20", "problem_ps_p0": 0.2},
        ]
        automated = {**self._receipt(), "selection_method": "automated_minimum_score"}
        extra_score = {**self._receipt(), "selection_score": 0.1}
        for receipt in (unknown, mismatched, multiple, automated, extra_score):
            with self.subTest(receipt=receipt), self.assertRaises(
                selection.PressureSelectionReceiptError
            ):
                self._validate(receipt)

    def test_requires_human_identity_canonical_utc_and_nonempty_rationale(self) -> None:
        cases = [
            {**self._receipt(), "reviewer_identity": "  "},
            {**self._receipt(), "reviewed_utc": "2026-06-02"},
            {**self._receipt(), "reviewed_utc": "2026-06-02T12:34:56+00:00"},
            {**self._receipt(), "rationale": ""},
            {**self._receipt(), "rationale": " \t "},
        ]
        for receipt in cases:
            with self.subTest(receipt=receipt), self.assertRaises(
                selection.PressureSelectionReceiptError
            ):
                self._validate(receipt)

    def test_canonical_bytes_reject_duplicate_nonfinite_and_reformatted_json(self) -> None:
        receipt = self._receipt()
        payload = selection.canonical_json_bytes(receipt)
        self.assertEqual(
            selection.validate_pressure_selection_receipt_bytes(
                payload,
                authorized_pic_root=self.base.parent,
            ),
            receipt,
        )
        duplicate = payload.replace(
            b'  "schema_version": 1,',
            b'  "schema_version": 1,\n  "schema_version": 1,',
            1,
        )
        nonfinite = payload.replace(b'"problem_ps_p0": 0.1', b'"problem_ps_p0": NaN', 1)
        reformatted = payload.replace(b"\n", b"", 1)
        for invalid in (duplicate, nonfinite, reformatted):
            with self.subTest(payload=invalid), self.assertRaises(
                selection.PressureSelectionReceiptError
            ):
                selection.validate_pressure_selection_receipt_bytes(
                    invalid,
                    authorized_pic_root=self.base.parent,
                )

    def test_rejects_numeric_boolean_schema_alias_and_input_mutation_isolated(self) -> None:
        receipt = self._receipt()
        alias = {**receipt, "schema_version": True}
        with self.assertRaises(selection.PressureSelectionReceiptError):
            self._validate(alias)

        validated = self._validate(receipt)
        receipt["selected_case"]["problem_ps_p0"] = 0.2
        self.assertEqual(validated["selected_case"]["problem_ps_p0"], 0.1)


if __name__ == "__main__":
    unittest.main()
