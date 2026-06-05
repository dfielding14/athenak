#!/usr/bin/env python3
"""Focused tests for the human Q-011 Section 5.4 pressure-selection boundary."""

from __future__ import annotations

import copy
from contextlib import contextmanager
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
        "schema_version": 3,
        "record_type": selection.RECORD_TYPE,
        "selection_method": selection.SELECTION_METHOD,
        "published_pressure_pilot_receipt": {
            "path": publication["receipt_path"],
            "sha256": publication["receipt_sha256"],
        },
        "published_pressure_pilot_review_packet_receipt": {
            "path": str(Path(publication["receipt_path"]).with_name("review-packet-receipt.json")),
            "sha256": "9" * 64,
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
        "authoritative_reanalysis_attestation": {
            "path": str(
                Path(publication["receipt_path"]).with_name(
                    "pressure-reanalysis-attestation.json"
                )
            ),
            "sha256": "7" * 64,
        },
        "reviewer_attestation": {
            "path": str(
                Path(publication["receipt_path"]).with_name(
                    "pressure-reviewer-attestation.json"
                )
            ),
            "sha256": "8" * 64,
        },
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

    def _packet_verification(self, receipt: dict[str, object]) -> dict[str, object]:
        aggregate_receipt = json.loads(self.receipt_path.read_text(encoding="utf-8"))
        return {
            "receipt_binding": copy.deepcopy(
                receipt.get(
                    "published_pressure_pilot_review_packet_receipt",
                    {"path": "/invalid/missing-review-packet-receipt", "sha256": "0" * 64},
                )
            ),
            "aggregate_receipt_binding": copy.deepcopy(
                receipt["published_pressure_pilot_receipt"]
            ),
            "packet_receipt": {},
            "aggregate_receipt": aggregate_receipt,
            "aggregate_bundle": copy.deepcopy(aggregate_receipt["aggregate_bundle"]),
            "aggregate_analysis": copy.deepcopy(aggregate_receipt["aggregate_analysis"]),
            "source_bindings": {},
            "inventory": {},
        }

    def _aggregate_verification(self) -> dict[str, object]:
        return {
            "packet_receipt_sha256": "9" * 64,
            "aggregate_receipt_sha256": self.publication["receipt_sha256"],
            "manifest_sha256": self.publication["manifest_sha256"],
            "analysis_result_sha256": self.publication["analysis_result_sha256"],
            "status": "pass_engineering_calibration_only",
        }

    def _reanalysis_verification(
        self, receipt: dict[str, object]
    ) -> dict[str, object]:
        return {
            "binding": copy.deepcopy(
                receipt.get(
                    "authoritative_reanalysis_attestation",
                    {"path": "/invalid/reanalysis-attestation", "sha256": "0" * 64},
                )
            ),
            "attestation": {},
            "operator_id": "test-only-reanalysis-operator",
            "recomputed_utc": "2026-06-02T12:00:00Z",
            "sealed_utc": "2026-06-02T12:01:00Z",
            "evidence": {},
            "source_authorization": {},
            "result": self._aggregate_verification(),
        }

    def _reviewer_verification(self, receipt: dict[str, object]) -> dict[str, object]:
        return {
            "binding": copy.deepcopy(
                receipt.get(
                    "reviewer_attestation",
                    {"path": "/invalid/reviewer-attestation", "sha256": "0" * 64},
                )
            ),
            "attestation": {},
            "reviewer_id": "test-only-human-reviewer",
            "reviewed_utc": "2026-06-02T12:02:00Z",
            "sealed_utc": "2026-06-02T12:03:00Z",
            "rationale": "Test-only human rationale.",
            "selected_case": copy.deepcopy(receipt["selected_case"]),
            "authoritative_reanalysis_attestation": copy.deepcopy(
                receipt.get(
                    "authoritative_reanalysis_attestation",
                    {"path": "/invalid/reanalysis-attestation", "sha256": "0" * 64},
                )
            ),
        }

    @contextmanager
    def _verification_context(self, receipt: dict[str, object]):
        reanalysis = self._reanalysis_verification(receipt)
        with mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=self._packet_verification(receipt),
        ), mock.patch.object(
            selection.historical_pressure_pilot_consumer,
            "consume_exact_historical_production_pressure_pilot",
            return_value=self._aggregate_verification(),
        ), mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_sealed_pressure_reanalysis_attestation",
            return_value=reanalysis,
        ), mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_sealed_pressure_reviewer_attestation",
            return_value=self._reviewer_verification(receipt),
        ):
            yield

    def _validate(self, receipt: dict[str, object]) -> dict[str, object]:
        with self._verification_context(receipt):
            return selection.validate_pressure_selection_receipt(
                receipt,
                authorized_pic_root=self.base.parent,
            )

    def test_accepts_one_human_selection_from_verified_immutable_publication(self) -> None:
        receipt = self._receipt()
        calls = []

        def verify_packet(*args: object, **kwargs: object) -> dict[str, object]:
            calls.append("packet")
            return self._packet_verification(receipt)

        def verify_aggregate(*args: object, **kwargs: object) -> dict[str, object]:
            calls.append("aggregate")
            return self._aggregate_verification()

        reanalysis_verification = self._reanalysis_verification(receipt)

        def verify_reanalysis(*args: object, **kwargs: object) -> dict[str, object]:
            calls.append("reanalysis")
            return reanalysis_verification

        def verify_reviewer(*args: object, **kwargs: object) -> dict[str, object]:
            calls.append("reviewer")
            return self._reviewer_verification(receipt)

        with mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            side_effect=verify_packet,
        ) as packet_verifier, mock.patch.object(
            selection.historical_pressure_pilot_consumer,
            "consume_exact_historical_production_pressure_pilot",
            side_effect=verify_aggregate,
        ) as aggregate_verifier, mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_sealed_pressure_reanalysis_attestation",
            side_effect=verify_reanalysis,
        ) as reanalysis_verifier, mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_sealed_pressure_reviewer_attestation",
            side_effect=verify_reviewer,
        ) as reviewer_verifier:
            self.assertEqual(
                selection.validate_pressure_selection_receipt(
                    receipt,
                    authorized_pic_root=self.base.parent,
                ),
                receipt,
            )
        self.assertEqual(calls, ["packet", "aggregate", "reanalysis", "reviewer"])
        packet_verifier.assert_called_once_with(
            receipt["published_pressure_pilot_review_packet_receipt"]["path"],
            aggregate_receipt_binding=receipt["published_pressure_pilot_receipt"],
            authorized_pic_root=self.base.parent,
        )
        aggregate_verifier.assert_called_once_with()
        reanalysis_verifier.assert_called_once_with(
            receipt["authoritative_reanalysis_attestation"],
            aggregate_receipt_binding=receipt["published_pressure_pilot_receipt"],
            packet_receipt_binding=receipt[
                "published_pressure_pilot_review_packet_receipt"
            ],
            pilot_bundle_manifest_sha256=receipt["pilot_bundle_manifest_sha256"],
            aggregate_pilot_analysis_sha256=receipt[
                "aggregate_pilot_analysis_sha256"
            ],
            authorized_pic_root=self.base.parent,
            expected_result=self._aggregate_verification(),
        )
        reviewer_verifier.assert_called_once_with(
            receipt["reviewer_attestation"],
            aggregate_receipt_binding=receipt["published_pressure_pilot_receipt"],
            packet_receipt_binding=receipt[
                "published_pressure_pilot_review_packet_receipt"
            ],
            reanalysis_verification=reanalysis_verification,
            selected_case=receipt["selected_case"],
            authorized_pic_root=self.base.parent,
        )

    def test_rejects_sealed_attestation_cross_binding(self) -> None:
        receipt = self._receipt()
        wrong_reanalysis = self._reanalysis_verification(receipt)
        wrong_reanalysis["binding"]["sha256"] = "0" * 64
        with mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=self._packet_verification(receipt),
        ), mock.patch.object(
            selection.historical_pressure_pilot_consumer,
            "consume_exact_historical_production_pressure_pilot",
            return_value=self._aggregate_verification(),
        ), mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_sealed_pressure_reanalysis_attestation",
            return_value=wrong_reanalysis,
        ), mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_sealed_pressure_reviewer_attestation",
        ) as reviewer_verifier, self.assertRaisesRegex(
            selection.PressureSelectionReceiptError,
            "authoritative reanalysis attestation binding drifted",
        ):
            selection.validate_pressure_selection_receipt(
                receipt,
                authorized_pic_root=self.base.parent,
            )
        reviewer_verifier.assert_not_called()

        wrong_reviewer = self._reviewer_verification(receipt)
        wrong_reviewer["binding"]["sha256"] = "1" * 64
        with mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=self._packet_verification(receipt),
        ), mock.patch.object(
            selection.historical_pressure_pilot_consumer,
            "consume_exact_historical_production_pressure_pilot",
            return_value=self._aggregate_verification(),
        ), mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_sealed_pressure_reanalysis_attestation",
            return_value=self._reanalysis_verification(receipt),
        ), mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_sealed_pressure_reviewer_attestation",
            return_value=wrong_reviewer,
        ), self.assertRaisesRegex(
            selection.PressureSelectionReceiptError,
            "human pressure-selection reviewer attestation binding drifted",
        ):
            selection.validate_pressure_selection_receipt(
                receipt,
                authorized_pic_root=self.base.parent,
            )

    def test_normalizes_authoritative_aggregate_publisher_failure(self) -> None:
        receipt = self._receipt()
        calls = []

        def verify_packet(*args: object, **kwargs: object) -> dict[str, object]:
            calls.append("packet")
            return self._packet_verification(receipt)

        def reject_aggregate(*args: object, **kwargs: object) -> dict[str, object]:
            calls.append("aggregate")
            raise selection.historical_pressure_pilot_consumer.HistoricalPressurePilotConsumerError(
                "test-only rejection"
            )

        with mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            side_effect=verify_packet,
        ), mock.patch.object(
            selection.historical_pressure_pilot_consumer,
            "consume_exact_historical_production_pressure_pilot",
            side_effect=reject_aggregate,
        ), self.assertRaisesRegex(
            selection.PressureSelectionReceiptError,
            "published pressure-pilot receipt failed authoritative verification",
        ) as raised:
            selection.validate_pressure_selection_receipt(
                receipt,
                authorized_pic_root=self.base.parent,
            )
        self.assertEqual(calls, ["packet", "aggregate"])
        self.assertIsInstance(
            raised.exception.__cause__,
            selection.historical_pressure_pilot_consumer.HistoricalPressurePilotConsumerError,
        )

    def test_rejects_authoritative_aggregate_publication_drift(self) -> None:
        cases = (
            (
                "packet_receipt_sha256",
                "0" * 64,
                "packet receipt SHA-256 differs from immutable publication binding",
            ),
            (
                "aggregate_receipt_sha256",
                "0" * 64,
                "receipt SHA-256 differs from immutable publication binding",
            ),
            (
                "manifest_sha256",
                "1" * 64,
                "manifest SHA-256 differs from immutable publication binding",
            ),
            (
                "analysis_result_sha256",
                "2" * 64,
                "analysis SHA-256 differs from immutable publication binding",
            ),
            (
                "status",
                "fail",
                "aggregate status is not pass_engineering_calibration_only",
            ),
        )
        for field, value, message in cases:
            receipt = self._receipt()
            aggregate_verification = self._aggregate_verification()
            aggregate_verification[field] = value
            with self.subTest(field=field), mock.patch.object(
                selection.pressure_review_packet_verifier,
                "consume_published_pressure_pilot_review_packet",
                return_value=self._packet_verification(receipt),
            ), mock.patch.object(
                selection.historical_pressure_pilot_consumer,
                "consume_exact_historical_production_pressure_pilot",
                return_value=aggregate_verification,
            ), self.assertRaisesRegex(
                selection.PressureSelectionReceiptError,
                message,
            ):
                selection.validate_pressure_selection_receipt(
                    receipt,
                    authorized_pic_root=self.base.parent,
                )

    def test_requires_exact_verified_review_packet_and_aggregate_cross_binding(self) -> None:
        schema_one = self._receipt()
        schema_one["schema_version"] = 1
        missing = self._receipt()
        del missing["published_pressure_pilot_review_packet_receipt"]
        malformed = self._receipt()
        malformed["published_pressure_pilot_review_packet_receipt"]["sha256"] = "0"
        for receipt in (schema_one, missing, malformed):
            with self.subTest(receipt=receipt), self.assertRaises(
                selection.PressureSelectionReceiptError
            ):
                self._validate(receipt)

        rejected = self._receipt()
        with mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            side_effect=selection.pressure_review_packet_verifier.PressureReviewPacketVerificationError(
                "test-only packet rejection"
            ),
        ), self.assertRaisesRegex(
            selection.PressureSelectionReceiptError,
            "review-packet receipt failed immutable verification",
        ):
            selection.validate_pressure_selection_receipt(
                rejected,
                authorized_pic_root=self.base.parent,
            )

        wrong_packet = self._receipt()
        packet_verification = self._packet_verification(wrong_packet)
        packet_verification["receipt_binding"]["sha256"] = "8" * 64
        with mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=packet_verification,
        ), self.assertRaisesRegex(
            selection.PressureSelectionReceiptError,
            "review-packet receipt binding drifted",
        ):
            selection.validate_pressure_selection_receipt(
                wrong_packet,
                authorized_pic_root=self.base.parent,
            )

        wrong_aggregate = self._receipt()
        packet_verification = self._packet_verification(wrong_aggregate)
        packet_verification["aggregate_receipt_binding"]["sha256"] = "7" * 64
        with mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=packet_verification,
        ), self.assertRaisesRegex(
            selection.PressureSelectionReceiptError,
            "different aggregate receipt",
        ):
            selection.validate_pressure_selection_receipt(
                wrong_aggregate,
                authorized_pic_root=self.base.parent,
            )

    def test_rejects_self_asserted_arbitrary_publication_and_evidence_hashes(self) -> None:
        cases = []
        publication_receipt = self._receipt()
        publication_receipt["published_pressure_pilot_receipt"]["sha256"] = "0" * 64
        with mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            side_effect=selection.pressure_review_packet_verifier.PressureReviewPacketVerificationError(
                "test-only forged aggregate binding rejection"
            ),
        ), self.assertRaises(selection.PressureSelectionReceiptError):
            selection.validate_pressure_selection_receipt(
                publication_receipt,
                authorized_pic_root=self.base.parent,
            )
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

    def test_does_not_consume_claimed_aggregate_before_packet_verification(self) -> None:
        receipt = self._receipt()
        receipt["published_pressure_pilot_receipt"] = {
            "path": str(self.base / "not-present-aggregate-receipt.json"),
            "sha256": receipt["published_pressure_pilot_receipt"]["sha256"],
        }
        with mock.patch.object(
            selection.pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            side_effect=selection.pressure_review_packet_verifier.PressureReviewPacketVerificationError(
                "test-only packet rejection"
            ),
        ), mock.patch.object(
            selection.historical_pressure_pilot_consumer,
            "consume_exact_historical_production_pressure_pilot",
            side_effect=AssertionError("aggregate publisher consumer must not be called"),
        ) as aggregate_verifier, self.assertRaisesRegex(
            selection.PressureSelectionReceiptError,
            "review-packet receipt failed immutable verification",
        ):
            selection.validate_pressure_selection_receipt(
                receipt,
                authorized_pic_root=self.base.parent,
            )
        aggregate_verifier.assert_not_called()

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

    def test_requires_exact_v3_sealed_attestation_bindings(self) -> None:
        schema_two = self._receipt()
        schema_two["schema_version"] = 2
        missing_reanalysis = self._receipt()
        del missing_reanalysis["authoritative_reanalysis_attestation"]
        missing_reviewer = self._receipt()
        del missing_reviewer["reviewer_attestation"]
        malformed_reanalysis = self._receipt()
        malformed_reanalysis["authoritative_reanalysis_attestation"]["sha256"] = "0"
        malformed_reviewer = self._receipt()
        malformed_reviewer["reviewer_attestation"]["path"] = ""
        legacy_inline_reviewer = {
            **self._receipt(),
            "reviewer_identity": "legacy inline reviewer",
        }
        cases = (
            schema_two,
            missing_reanalysis,
            missing_reviewer,
            malformed_reanalysis,
            malformed_reviewer,
            legacy_inline_reviewer,
        )
        for receipt in cases:
            with self.subTest(receipt=receipt), self.assertRaises(
                selection.PressureSelectionReceiptError
            ):
                self._validate(receipt)

    def test_canonical_bytes_reject_duplicate_nonfinite_and_reformatted_json(self) -> None:
        receipt = self._receipt()
        payload = selection.canonical_json_bytes(receipt)
        with self._verification_context(receipt):
            self.assertEqual(
                selection.validate_pressure_selection_receipt_bytes(
                    payload,
                    authorized_pic_root=self.base.parent,
                ),
                receipt,
            )
        duplicate = payload.replace(
            b'  "schema_version": 3,',
            b'  "schema_version": 3,\n  "schema_version": 3,',
            1,
        )
        nonfinite = payload.replace(b'"problem_ps_p0": 0.1', b'"problem_ps_p0": NaN', 1)
        reformatted = payload.replace(b"\n", b"", 1)
        with self._verification_context(receipt):
            for invalid in (duplicate, nonfinite, reformatted):
                with self.subTest(payload=invalid), self.assertRaises(
                    selection.PressureSelectionReceiptError
                ):
                    selection.validate_pressure_selection_receipt_bytes(
                        invalid,
                        authorized_pic_root=self.base.parent,
                    )

    def test_receipt_size_and_recursion_limits_fail_closed(self) -> None:
        oversized = b" " * (selection.MAX_PRESSURE_SELECTION_RECEIPT_BYTES + 1)
        deeply_nested = b'{"nested":' + (b"[" * 2000) + b"0" + (b"]" * 2000) + b"}\n"
        with self.assertRaisesRegex(
            selection.PressureSelectionReceiptError,
            "size limit",
        ):
            selection.validate_pressure_selection_receipt_bytes(
                oversized,
                authorized_pic_root=self.base.parent,
            )
        with self.assertRaisesRegex(
            selection.PressureSelectionReceiptError,
            "not valid UTF-8 JSON",
        ):
            selection.validate_pressure_selection_receipt_bytes(
                deeply_nested,
                authorized_pic_root=self.base.parent,
            )

        nested: object = 0
        for _ in range(2000):
            nested = [nested]
        with self.assertRaisesRegex(
            selection.PressureSelectionReceiptError,
            "not canonical JSON",
        ):
            selection.canonical_json_bytes(nested)

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
