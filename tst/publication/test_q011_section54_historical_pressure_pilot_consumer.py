#!/usr/bin/env python3
"""Focused adversarial tests for the exact historical pressure-pilot consumer."""

from __future__ import annotations

from contextlib import contextmanager
import copy
import hashlib
import inspect
import json
import os
from pathlib import Path
import stat
import tempfile
from typing import Iterator
import unittest
from unittest import mock

from tst.publication import q011_section54_historical_pressure_pilot_consumer as consumer


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


class Q011Section54HistoricalPressurePilotConsumerTests(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.root = Path(self._temporary.name)
        self.pic_root = self.root / "pic"
        self.publication_root = self.pic_root / "publication"
        self.bundle = self.publication_root / "bundle"
        self.nested = self.bundle / "nested"
        self.nested.mkdir(parents=True)
        self.manifest_payload = _canonical_json_bytes({"historical_manifest": True})
        self.member_payload = b"retained historical member\n"
        (self.bundle / consumer.AGGREGATE_MANIFEST_NAME).write_bytes(
            self.manifest_payload
        )
        (self.nested / "member.bin").write_bytes(self.member_payload)
        self.analysis_result = {
            "schema_version": 1,
            "status": "pass_engineering_calibration_only",
        }
        self.analysis_payload = _canonical_json_bytes(self.analysis_result)
        self.analysis_path = self.publication_root / "analysis.json"
        self.analysis_path.write_bytes(self.analysis_payload)
        self.packet_binding = {
            "path": str(self.publication_root / "packet-receipt.json"),
            "sha256": "1" * 64,
        }
        self.aggregate_receipt_binding = {
            "path": str(self.publication_root / "aggregate-receipt.json"),
            "sha256": "2" * 64,
        }
        self.bundle_binding = {
            "path": str(self.bundle),
            "manifest_sha256": _sha256(self.manifest_payload),
        }
        self.analysis_binding = {
            "path": str(self.analysis_path),
            "sha256": _sha256(self.analysis_payload),
        }
        self._freeze_bundle()

    def tearDown(self) -> None:
        for directory, _names, filenames in os.walk(self.root):
            base = Path(directory)
            base.chmod(0o755)
            for name in filenames:
                path = base / name
                if not path.is_symlink():
                    path.chmod(0o644)
        self._temporary.cleanup()

    def _freeze_bundle(self) -> None:
        for path in self.bundle.rglob("*"):
            if path.is_file():
                path.chmod(0o444)
        for path in sorted(self.bundle.rglob("*"), reverse=True):
            if path.is_dir():
                path.chmod(0o555)
        self.bundle.chmod(0o555)
        self.analysis_path.chmod(0o444)

    def _verification(self) -> dict[str, object]:
        return {
            "receipt_binding": copy.deepcopy(self.packet_binding),
            "aggregate_receipt_binding": copy.deepcopy(
                self.aggregate_receipt_binding
            ),
            "packet_receipt": {"record_type": "fixture_packet_receipt"},
            "aggregate_receipt": {"record_type": "fixture_aggregate_receipt"},
            "aggregate_bundle": copy.deepcopy(self.bundle_binding),
            "aggregate_analysis": copy.deepcopy(self.analysis_binding),
            "source_bindings": {"fixture": True},
            "inventory": {"fixture": True},
        }

    @contextmanager
    def _contract(self, **extra: object) -> Iterator[None]:
        values = {
            "AUTHORIZED_PRODUCTION_PIC_ROOT": self.pic_root,
            "AUTHORIZED_PRODUCTION_PUBLICATION_ROOT": self.publication_root,
            "AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING": self.packet_binding,
            "AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING": (
                self.aggregate_receipt_binding
            ),
            "AUTHORIZED_PRODUCTION_AGGREGATE_BUNDLE_BINDING": self.bundle_binding,
            "AUTHORIZED_PRODUCTION_AGGREGATE_ANALYSIS_BINDING": self.analysis_binding,
            **extra,
        }
        with mock.patch.multiple(consumer, **values):
            yield

    def test_public_api_exposes_no_paths_bypass_flags_or_publication_api(self) -> None:
        self.assertEqual(
            list(
                inspect.signature(
                    consumer.consume_exact_historical_production_pressure_pilot
                ).parameters
            ),
            [],
        )
        self.assertEqual(
            consumer.__all__,
            [
                "HistoricalPressurePilotConsumerError",
                "consume_exact_historical_production_pressure_pilot",
            ],
        )
        self.assertFalse(any(name.startswith("publish") for name in consumer.__all__))

    def test_accepts_exact_pair_and_recomputes_through_bounded_member_reader(
        self,
    ) -> None:
        calls: list[str] = []

        def analyze(
            root: Path,
            expected_manifest_sha256: str,
            **kwargs: object,
        ) -> dict[str, object]:
            calls.append("analyze")
            self.assertEqual(root, self.bundle)
            self.assertEqual(expected_manifest_sha256, self.bundle_binding["manifest_sha256"])
            self.assertEqual(
                kwargs["authorized_publication_root"],
                self.publication_root,
            )
            member_reader = kwargs["member_reader"]
            self.assertEqual(
                member_reader(consumer.AGGREGATE_MANIFEST_NAME),
                self.manifest_payload,
            )
            self.assertEqual(member_reader("nested/member.bin"), self.member_payload)
            self.assertEqual(
                kwargs["actual_files"],
                {consumer.AGGREGATE_MANIFEST_NAME, "nested/member.bin"},
            )
            return copy.deepcopy(self.analysis_result)

        with self._contract(), mock.patch.object(
            consumer._pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            side_effect=lambda *_args, **_kwargs: copy.deepcopy(self._verification()),
        ) as packet_verifier, mock.patch.object(
            consumer._pressure_pilot_analyzer,
            "analyze_exact_historical_production_pressure_pilot_bundle",
            side_effect=analyze,
            create=True,
        ) as analyzer:
            result = consumer.consume_exact_historical_production_pressure_pilot()

        self.assertEqual(calls, ["analyze"])
        self.assertEqual(
            result,
            {
                "packet_receipt_sha256": self.packet_binding["sha256"],
                "aggregate_receipt_sha256": self.aggregate_receipt_binding["sha256"],
                "manifest_sha256": self.bundle_binding["manifest_sha256"],
                "analysis_result_sha256": self.analysis_binding["sha256"],
                "status": self.analysis_result["status"],
            },
        )
        analyzer.assert_called_once()
        self.assertEqual(packet_verifier.call_count, 2)
        for invocation in packet_verifier.call_args_list:
            self.assertEqual(invocation.args, (self.packet_binding["path"],))
            self.assertEqual(
                invocation.kwargs,
                {
                    "aggregate_receipt_binding": self.aggregate_receipt_binding,
                    "authorized_pic_root": self.pic_root,
                },
            )

    def test_rejects_packet_verifier_result_bound_to_different_exact_pair(self) -> None:
        drifted = self._verification()
        drifted["receipt_binding"]["sha256"] = "f" * 64
        with self._contract(), mock.patch.object(
            consumer._pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=drifted,
        ), mock.patch.object(
            consumer._pressure_pilot_analyzer,
            "analyze_exact_historical_production_pressure_pilot_bundle",
            create=True,
        ) as analyzer, self.assertRaisesRegex(
            consumer.HistoricalPressurePilotConsumerError,
            "different packet receipt",
        ):
            consumer.consume_exact_historical_production_pressure_pilot()
        analyzer.assert_not_called()

    def test_rejects_exact_hash_bound_retained_analysis_byte_mismatch(self) -> None:
        different = _canonical_json_bytes(
            {
                "schema_version": 1,
                "status": "different_retained_status",
            }
        )
        self.analysis_path.chmod(0o600)
        self.analysis_path.write_bytes(different)
        self.analysis_path.chmod(0o444)
        self.analysis_binding["sha256"] = _sha256(different)
        with self._contract(), mock.patch.object(
            consumer._pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            side_effect=lambda *_args, **_kwargs: copy.deepcopy(self._verification()),
        ), mock.patch.object(
            consumer._pressure_pilot_analyzer,
            "analyze_exact_historical_production_pressure_pilot_bundle",
            return_value=copy.deepcopy(self.analysis_result),
            create=True,
        ), self.assertRaisesRegex(
            consumer.HistoricalPressurePilotConsumerError,
            "differs from retained analysis bytes",
        ):
            consumer.consume_exact_historical_production_pressure_pilot()

    def test_rejects_symlink_bundle_member_before_analyzer(self) -> None:
        self.bundle.chmod(0o755)
        (self.bundle / "linked-member").symlink_to("nested/member.bin")
        self.bundle.chmod(0o555)
        with self._contract(), mock.patch.object(
            consumer._pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=self._verification(),
        ), mock.patch.object(
            consumer._pressure_pilot_analyzer,
            "analyze_exact_historical_production_pressure_pilot_bundle",
            create=True,
        ) as analyzer, self.assertRaisesRegex(
            consumer.HistoricalPressurePilotConsumerError,
            "unsupported file type",
        ):
            consumer.consume_exact_historical_production_pressure_pilot()
        analyzer.assert_not_called()

    def test_rejects_actual_oversized_bundle_member_before_read(self) -> None:
        self.bundle.chmod(0o755)
        oversized = self.bundle / "oversized.bin"
        with oversized.open("wb") as stream:
            stream.truncate(consumer.MAX_RETAINED_FILE_BYTES + 1)
        oversized.chmod(0o444)
        self.bundle.chmod(0o555)
        with self._contract(), mock.patch.object(
            consumer._pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=self._verification(),
        ), mock.patch.object(
            consumer.os,
            "read",
            side_effect=AssertionError("oversized bundle member was read"),
        ), mock.patch.object(
            consumer._pressure_pilot_analyzer,
            "analyze_exact_historical_production_pressure_pilot_bundle",
            create=True,
        ) as analyzer, self.assertRaisesRegex(
            consumer.HistoricalPressurePilotConsumerError,
            "size limit",
        ):
            consumer.consume_exact_historical_production_pressure_pilot()
        analyzer.assert_not_called()

    def test_rejects_directory_entry_limit_before_analyzer(self) -> None:
        with self._contract(MAX_DIRECTORY_ENTRIES=1), mock.patch.object(
            consumer._pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=self._verification(),
        ), mock.patch.object(
            consumer._pressure_pilot_analyzer,
            "analyze_exact_historical_production_pressure_pilot_bundle",
            create=True,
        ) as analyzer, self.assertRaisesRegex(
            consumer.HistoricalPressurePilotConsumerError,
            "entry limit",
        ):
            consumer.consume_exact_historical_production_pressure_pilot()
        analyzer.assert_not_called()

    def test_rejects_same_account_member_replacement_during_recompute(self) -> None:
        member = self.nested / "member.bin"

        def replace_then_read(
            _root: Path,
            _manifest_sha256: str,
            **kwargs: object,
        ) -> dict[str, object]:
            self.nested.chmod(0o755)
            member.unlink()
            member.write_bytes(b"same-account replacement\n")
            member.chmod(0o444)
            self.nested.chmod(0o555)
            kwargs["member_reader"]("nested/member.bin")
            return copy.deepcopy(self.analysis_result)

        with self._contract(), mock.patch.object(
            consumer._pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=self._verification(),
        ), mock.patch.object(
            consumer._pressure_pilot_analyzer,
            "analyze_exact_historical_production_pressure_pilot_bundle",
            side_effect=replace_then_read,
            create=True,
        ), self.assertRaisesRegex(
            consumer.HistoricalPressurePilotConsumerError,
            "changed before",
        ):
            consumer.consume_exact_historical_production_pressure_pilot()

    def test_repeated_directory_identity_failures_do_not_leak_descriptors(self) -> None:
        root_fd = os.open(self.bundle, consumer._open_flags(directory=True))
        try:
            root_metadata = os.fstat(root_fd)
            reader = consumer._ImmutableBundleReader(
                root_fd,
                consumer._scan_bundle(root_fd, root_metadata),
            )
            original_fstat = os.fstat

            def drifted_identity(fd: int) -> os.stat_result:
                metadata = original_fstat(fd)
                values = list(metadata)
                values[stat.ST_INO] += 1
                return os.stat_result(values)

            baseline = len(os.listdir("/proc/self/fd"))
            with mock.patch.object(consumer.os, "fstat", side_effect=drifted_identity):
                for _index in range(25):
                    with self.assertRaisesRegex(
                        consumer.HistoricalPressurePilotConsumerError,
                        "directory changed while being opened: nested",
                    ):
                        reader.read("nested/member.bin")
                    self.assertEqual(len(os.listdir("/proc/self/fd")), baseline)
        finally:
            os.close(root_fd)

    def test_rejects_packet_verification_drift_after_recomputation(self) -> None:
        first = self._verification()
        second = self._verification()
        second["inventory"] = {"fixture": False}
        with self._contract(), mock.patch.object(
            consumer._pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            side_effect=[first, second],
        ), mock.patch.object(
            consumer._pressure_pilot_analyzer,
            "analyze_exact_historical_production_pressure_pilot_bundle",
            return_value=copy.deepcopy(self.analysis_result),
            create=True,
        ), self.assertRaisesRegex(
            consumer.HistoricalPressurePilotConsumerError,
            "verification changed",
        ):
            consumer.consume_exact_historical_production_pressure_pilot()

    def test_rejects_recursive_recomputed_analysis_fail_closed(self) -> None:
        deep: object = "terminal"
        for _index in range(2000):
            deep = [deep]
        with self._contract(), mock.patch.object(
            consumer._pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=self._verification(),
        ), mock.patch.object(
            consumer._pressure_pilot_analyzer,
            "analyze_exact_historical_production_pressure_pilot_bundle",
            return_value={"status": "pass_engineering_calibration_only", "deep": deep},
            create=True,
        ), self.assertRaisesRegex(
            consumer.HistoricalPressurePilotConsumerError,
            "not canonical JSON",
        ):
            consumer.consume_exact_historical_production_pressure_pilot()


if __name__ == "__main__":
    unittest.main()
