#!/usr/bin/env python3
"""Focused adversarial tests for Stage-4 pressure-selection publication."""

from __future__ import annotations

from contextlib import contextmanager
from contextlib import redirect_stderr, redirect_stdout
import copy
from datetime import datetime, timedelta, timezone
import fcntl
import hashlib
import io
import inspect
import json
import os
from pathlib import Path
import stat
import subprocess
import tarfile
import tempfile
from typing import Iterator
import unittest
from unittest import mock

from . import publish_q011_section54_pressure_selection as publisher


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
        "utf-8"
    )


def _git_style_archive(payload: bytes) -> bytes:
    stream = io.BytesIO()
    with tarfile.open(
        fileobj=stream,
        mode="w",
        format=tarfile.PAX_FORMAT,
        pax_headers={"comment": "e" * 40},
    ) as archive:
        archive.addfile(tarfile.TarInfo("unused"))
    header = bytearray(stream.getvalue()[: tarfile.BLOCKSIZE])
    name = b"pax_global_header"
    header[:100] = name + b"\0" * (100 - len(name))
    header[148:156] = b"        "
    header[148:156] = f"{sum(header):06o}\0 ".encode("ascii")
    return bytes(header) + stream.getvalue()[512:1024] + payload


class PressureSelectionPublisherTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.pic_root = self.root / "pic"
        self.project_home_root = self.root / "project-home"
        self.publication_root = self.pic_root / "publication"
        self.acceptance_root = self.pic_root / "publication_acceptance"
        self.archive_root = (
            self.pic_root / publisher.PRESSURE_GATE_ATTESTATION_ROOT_NAME
        )
        for path in (
            self.publication_root,
            self.acceptance_root,
            self.archive_root,
            self.pic_root / "policy",
            self.pic_root / "ledger",
            self.project_home_root / "policy",
        ):
            path.mkdir(parents=True, exist_ok=True)
        self.archive_root.chmod(0o700)
        self.version = "a" * 64
        self.manifest_path = (
            self.pic_root
            / "clean_candidates"
            / "candidate"
            / "clean_candidate_manifest.json"
        )
        self.clean_candidate_manifest: dict[str, object] = {
            "schema_version": 4,
            "source": {
                "git_commit": "e" * 40,
                "archive_sha256": "f" * 64,
            },
        }
        self._write_readonly(
            self.manifest_path, _canonical(self.clean_candidate_manifest)
        )
        self.policy: dict[str, object] = {
            "registered_science_slices": [],
            "frontier_admission_smoke": {
                "status": publisher.control_plane_common.CLOSED_ADMISSION_SMOKE_STATUS
            },
            "science_submission_freeze": {
                "status": publisher.control_plane_common.AUTHORIZED_CLEAN_CANDIDATE_FREEZE,
                "manifest_path": str(self.manifest_path),
                "manifest_sha256": _sha256(_canonical(self.clean_candidate_manifest)),
                "build_profile_control_plane_version": self.version,
            },
            "olcf_side_storage": {
                "installed_control_plane_version": self.version,
                "manual_accounting_authorizations": [
                    {
                        "authorization_id": "historical-accounting-only",
                        "path": str(
                            self.pic_root
                            / "policy"
                            / "manual_accounting_authorizations"
                            / "historical-accounting-only.json"
                        ),
                        "project_home_path": str(
                            self.project_home_root
                            / "policy"
                            / "manual_accounting_authorizations"
                            / "historical-accounting-only.json"
                        ),
                        "sha256": "b" * 64,
                    }
                ],
            },
        }
        self.promotion: dict[str, object] = {
            "schema_version": 2,
            "promotion_id": "d920ed7c-c54c-4adb-8510-5363efbe1932",
            "control_plane_version": self.version,
            "policy_path": str(self.pic_root / "policy" / "storage_policy.json"),
            "project_home_policy_path": str(
                self.project_home_root / "policy" / "storage_policy.json"
            ),
            "policy_sha256": _sha256(_canonical(self.policy)),
        }
        self.inventory: dict[str, object] = {
            "schema_version": 1,
            "version": self.version,
            "files": [],
        }
        self._write_active_state()
        self._write_installed_inventories()
        self.receipt_path = (
            self.publication_root / publisher.CANONICAL_RECEIPT_NAME
        )
        self.reanalysis_binding = {
            "path": str(
                self.archive_root
                / "20260605T120000Z-q011-section54-pressure-reanalysis-codex"
                / "attestation.json"
            ),
            "sha256": "c" * 64,
        }
        self.reviewer_binding = {
            "path": str(
                self.archive_root
                / "20260605T120100Z-q011-section54-pressure-selection-dfielding"
                / "attestation.json"
            ),
            "sha256": "d" * 64,
        }
        self.receipt = publisher.build_pressure_selection_receipt(
            published_pressure_pilot_receipt={
                "path": str(self.publication_root / "aggregate-receipt.json"),
                "sha256": "1" * 64,
            },
            published_pressure_pilot_review_packet_receipt={
                "path": str(self.publication_root / "review-packet-receipt.json"),
                "sha256": "2" * 64,
            },
            pilot_bundle_manifest_sha256="3" * 64,
            aggregate_pilot_analysis_sha256="4" * 64,
            case_descriptors=[
                {
                    "case_id": case_id,
                    "problem_ps_p0": p0,
                    "descriptor_sha256": f"{index + 5:064x}",
                }
                for index, (case_id, p0) in enumerate(
                    publisher.pressure_selection.REGISTERED_CASES
                )
            ],
            authoritative_reanalysis_attestation=self.reanalysis_binding,
            reviewer_attestation=self.reviewer_binding,
        )
        self.stage4_preparation_binding = {
            "path": str(
                self.archive_root
                / "20260605T120100Z-q011-section54-pressure-stage4-preparation-codex"
                / "attestation.json"
            ),
            "sha256": "9" * 64,
        }
        self.human_decision_binding = {
            "path": str(
                self.pic_root
                / publisher.PRESSURE_GATE_HUMAN_DECISION_ROOT_NAME
                / "dfielding-p0-1p0.json"
            ),
            "sha256": "",
        }
        self.reviewer_id = publisher.REVIEWED_REVIEWER_ID
        self.rationale = publisher.REVIEWED_RATIONALE
        self.reviewer_reviewed_utc = "2026-06-05T12:02:00Z"
        self.reanalysis_source_archive_sha256 = "f" * 64
        self.reanalysis_source_closure = [
            {"path": path, "sha256": f"{index + 10:064x}"}
            for index, path in enumerate(
                publisher.pressure_selection.pressure_review_packet_verifier.PRESSURE_REANALYSIS_SOURCE_PATHS
            )
        ]
        self.candidate_source_authorization = self._reanalysis_source_authorization()
        self.now = datetime(2026, 6, 5, 12, 2, 0, tzinfo=timezone.utc)
        human_decision_path = Path(self.human_decision_binding["path"])
        self._write_readonly(
            human_decision_path,
            _canonical(
                {
                    "schema_version": 1,
                    "record_type": publisher.HUMAN_DECISION_RECORD_TYPE,
                    "reviewer_id": publisher.REVIEWED_REVIEWER_ID,
                    "reviewed_utc": "2026-06-05T12:02:00Z",
                    "rationale": publisher.REVIEWED_RATIONALE,
                    "reviewer_statement": publisher.HUMAN_DECISION_STATEMENT,
                    "authoritative_reanalysis_attestation": self.reanalysis_binding,
                    "stage4_preparation_attestation": self.stage4_preparation_binding,
                    "selected_case": publisher.REVIEWED_SELECTED_CASE,
                }
            ),
        )
        human_decision_path.chmod(0o400)
        human_decision_path.parent.chmod(0o700)
        self.human_decision_binding["sha256"] = _sha256(
            human_decision_path.read_bytes()
        )

    def _write_readonly(self, path: Path, payload: bytes) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        path.chmod(0o444)

    def _write_active_state(self) -> None:
        policy_payload = _canonical(self.policy)
        self.promotion["policy_sha256"] = _sha256(policy_payload)
        promotion_payload = _canonical(self.promotion)
        for root in (self.pic_root, self.project_home_root):
            policy_path = root / "policy" / "storage_policy.json"
            promotion_path = root / "policy" / "active_promotion.json"
            for path in (policy_path, promotion_path):
                if path.exists():
                    path.chmod(0o600)
            self._write_readonly(policy_path, policy_payload)
            self._write_readonly(promotion_path, promotion_payload)

    def _write_installed_inventories(self) -> None:
        payload = _canonical(self.inventory)
        for root in (self.pic_root, self.project_home_root):
            self._write_readonly(
                root / "control_plane" / self.version / "inventory.json", payload
            )

    def _install_reanalysis_source_archive(self) -> Path:
        source_root = self.root / "source"
        source_payloads = {
            relative: f"# exact archived source: {relative}\n".encode("utf-8")
            for relative in (
                publisher.pressure_selection.pressure_review_packet_verifier.PRESSURE_REANALYSIS_SOURCE_PATHS
            )
        }
        archive_stream = io.BytesIO()
        with tarfile.open(fileobj=archive_stream, mode="w") as archive:
            for relative, payload in source_payloads.items():
                path = source_root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(payload)
                member = tarfile.TarInfo(relative)
                member.size = len(payload)
                member.mode = 0o644
                archive.addfile(member, io.BytesIO(payload))
        archive_payload = _git_style_archive(archive_stream.getvalue())
        archive_path = self.manifest_path.parent / "source.tar"
        self._write_readonly(archive_path, archive_payload)
        self.clean_candidate_manifest["source"]["archive_sha256"] = _sha256(
            archive_payload
        )
        self.manifest_path.chmod(0o600)
        self._write_readonly(
            self.manifest_path, _canonical(self.clean_candidate_manifest)
        )
        self.policy["science_submission_freeze"]["manifest_sha256"] = _sha256(
            _canonical(self.clean_candidate_manifest)
        )
        self._write_active_state()
        return source_root

    def _unlock_snapshot(self, **_kwargs: object) -> tuple[dict[str, object], dict[str, str]]:
        policy_payload = (self.pic_root / "policy" / "storage_policy.json").read_bytes()
        promotion_payload = (
            self.pic_root / "policy" / "active_promotion.json"
        ).read_bytes()
        return copy.deepcopy(self.policy), {
            "active_policy_sha256": _sha256(policy_payload),
            "active_promotion_sha256": _sha256(promotion_payload),
        }

    def _reanalysis(self, receipt: dict[str, object]) -> dict[str, object]:
        return {
            "binding": copy.deepcopy(receipt["authoritative_reanalysis_attestation"]),
            "sealed_utc": "2026-06-05T12:01:00Z",
            "source_authorization": self._reanalysis_source_authorization(),
        }

    def _reanalysis_source_authorization(self) -> dict[str, object]:
        closure = copy.deepcopy(self.reanalysis_source_closure)
        return {
            "execution_mode": (
                publisher.pressure_selection.pressure_review_packet_verifier.PRESSURE_REANALYSIS_EXECUTION_MODE
            ),
            "git_commit": "e" * 40,
            "source_archive_sha256": self.reanalysis_source_archive_sha256,
            "source_closure_sha256": _sha256(
                json.dumps(
                    closure,
                    separators=(",", ":"),
                    sort_keys=True,
                    allow_nan=False,
                ).encode("utf-8")
            ),
            "source_closure": closure,
            "historical_production_source_authorization": dict(
                publisher.pressure_selection.pressure_review_packet_verifier.AUTHORIZED_HISTORICAL_REANALYSIS_SOURCE_AUTHORIZATION
            ),
        }

    def _mirrored_ledger_state(self) -> dict[str, object]:
        return {
            "validator": "validated_read_only_mirrored_state_snapshot",
            "orion_ledger_path": str(self.pic_root / "ledger" / "node_hours.jsonl"),
            "orion_mirror_receipts_path": str(
                self.pic_root / "ledger" / "mirror_receipts.jsonl"
            ),
            "project_home_ledger_path": str(
                self.project_home_root / "ledger" / "node_hours.jsonl"
            ),
            "record_count": 5,
            "tail_event_sha256": "8" * 64,
            "active_reservation_ids": [],
            "currently_reserved_node_hours": 0.0,
            "cumulative_consumed_node_hours": 1.25,
            "pending_submission_marker": "absent",
            "pending_manual_accounting_markers": {
                "orion": "absent",
                "project_home": "absent",
            },
        }

    def _capture_mirrored_ledger_state(self, **_kwargs: object) -> dict[str, object]:
        markers = (
            (
                self.pic_root / "ledger" / "pending_submission.json",
                "pending scheduler submission blocks pressure selection",
            ),
            (
                self.pic_root / "ledger" / "pending_manual_accounting.json",
                "incomplete manual accounting blocks pressure selection",
            ),
            (
                self.project_home_root / "ledger" / "pending_manual_accounting.json",
                "mirrored incomplete manual accounting blocks pressure selection",
            ),
        )
        for path, message in markers:
            if path.exists():
                raise publisher.PressureSelectionPublicationError(message)
        return self._mirrored_ledger_state()

    def _reviewer(self, receipt: dict[str, object]) -> dict[str, object]:
        return {
            "binding": copy.deepcopy(receipt["reviewer_attestation"]),
            "reviewer_id": self.reviewer_id,
            "reviewed_utc": self.reviewer_reviewed_utc,
            "rationale": self.rationale,
            "selected_case": copy.deepcopy(receipt["selected_case"]),
            "authoritative_reanalysis_attestation": copy.deepcopy(
                receipt["authoritative_reanalysis_attestation"]
            ),
        }

    def _publisher_source_authentication(self) -> dict[str, object]:
        closure = [
            {"path": path, "sha256": "a" * 64}
            for path in publisher.PUBLISHER_SOURCE_PATHS
        ]
        return {
            "execution_mode": "worker_extracted_git_archive_expected_commit_verified",
            "git_commit": "e" * 40,
            "archive_sha256": "b" * 64,
            "source_closure_sha256": _sha256(
                json.dumps(
                    closure,
                    separators=(",", ":"),
                    sort_keys=True,
                    allow_nan=False,
                ).encode("utf-8")
            ),
            "source_closure": closure,
        }

    def _candidate_publication_authorization(
        self, receipt: dict[str, object] | None = None
    ) -> dict[str, object]:
        candidate = self.receipt if receipt is None else receipt
        return {
            "schema_version": 1,
            "record_type": publisher.CANDIDATE_AUTHORIZATION_RECORD_TYPE,
            "qualification_effect": publisher.QUALIFICATION_EFFECT,
            "sealed_utc": "2026-06-05T12:02:00Z",
            "publisher_source_authentication": self._publisher_source_authentication(),
            "stage4_preparation_attestation": copy.deepcopy(
                self.stage4_preparation_binding
            ),
            "human_decision": copy.deepcopy(self.human_decision_binding),
            "candidate_pressure_selection_receipt": {
                "path": str(
                    self.pic_root
                    / publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME
                    / "candidate.json"
                ),
                "sha256": _sha256(
                    publisher.pressure_selection.canonical_json_bytes(candidate)
                ),
            },
            "authoritative_reanalysis_attestation": copy.deepcopy(
                candidate["authoritative_reanalysis_attestation"]
            ),
            "reviewer_attestation": copy.deepcopy(candidate["reviewer_attestation"]),
        }

    @contextmanager
    def _verification_context(self) -> Iterator[None]:
        original_consume_preparation = publisher._consume_stage4_preparation_attestation

        def validate(
            receipt: object, *, authorized_pic_root: Path
        ) -> dict[str, object]:
            self.assertEqual(authorized_pic_root, self.pic_root)
            self.assertIsInstance(receipt, dict)
            return copy.deepcopy(receipt)

        def consume_reanalysis(
            binding: object, **_kwargs: object
        ) -> dict[str, object]:
            receipt = copy.deepcopy(self.receipt)
            receipt["authoritative_reanalysis_attestation"] = copy.deepcopy(binding)
            return self._reanalysis(receipt)

        def consume_reviewer(
            binding: object,
            *,
            selected_case: object,
            **_kwargs: object,
        ) -> dict[str, object]:
            receipt = copy.deepcopy(self.receipt)
            receipt["reviewer_attestation"] = copy.deepcopy(binding)
            receipt["selected_case"] = copy.deepcopy(selected_case)
            return self._reviewer(receipt)

        def consume_preparation(
            binding: object, *, authorized_pic_root: Path
        ) -> dict[str, object]:
            if isinstance(binding, dict) and Path(str(binding.get("path", ""))).exists():
                return original_consume_preparation(
                    binding, authorized_pic_root=authorized_pic_root
                )
            return {
                "binding": copy.deepcopy(binding),
                "attestation": {
                    "publisher_source_authentication": self._publisher_source_authentication(),
                    "authoritative_reanalysis_attestation": copy.deepcopy(
                        self.reanalysis_binding
                    ),
                },
            }

        with mock.patch.object(
            publisher,
            "_runtime_publisher_source_authentication",
            return_value=self._publisher_source_authentication(),
        ), mock.patch.object(
            publisher.pressure_selection,
            "validate_pressure_selection_receipt",
            side_effect=validate,
        ), mock.patch.object(
            publisher.pressure_selection.pressure_review_packet_verifier,
            "consume_sealed_pressure_reanalysis_attestation",
            side_effect=consume_reanalysis,
        ), mock.patch.object(
            publisher.pressure_selection.pressure_review_packet_verifier,
            "consume_sealed_pressure_reviewer_attestation",
            side_effect=consume_reviewer,
        ), mock.patch.object(
            publisher.pressure_selection.pressure_review_packet_verifier,
            "validate_pressure_reanalysis_source_snapshot",
            return_value=None,
        ), mock.patch.object(
            publisher,
            "_consume_stage4_preparation_attestation",
            side_effect=consume_preparation,
        ), mock.patch.object(
            publisher.control_plane_common,
            "require_storage_policy_unlock_snapshot",
            side_effect=self._unlock_snapshot,
        ), mock.patch.object(
            publisher.control_plane_common,
            "verify_installed_control_plane",
            return_value=copy.deepcopy(self.inventory),
        ), mock.patch.object(
            publisher.clean_candidate_revalidator,
            "revalidate_clean_candidate",
            return_value={
                "clean_candidate_manifest": {
                    "expected_sha256": self.policy["science_submission_freeze"][
                        "manifest_sha256"
                    ],
                    "path": str(self.manifest_path),
                    "sha256": self.policy["science_submission_freeze"][
                        "manifest_sha256"
                    ],
                },
                "current_control_plane_version": self.version,
                "source": {"git_commit": "e" * 40},
                "status": "passed",
            },
        ), mock.patch.object(
            publisher,
            "_capture_mirrored_ledger_state",
            side_effect=self._capture_mirrored_ledger_state,
        ), mock.patch.object(
            publisher,
            "_candidate_reanalysis_source_authorization",
            side_effect=lambda *_args, **_kwargs: (
                copy.deepcopy(self.candidate_source_authorization),
                {},
            ),
        ):
            yield

    def _publish(self) -> dict[str, object]:
        return publisher.publish_pressure_selection(
            copy.deepcopy(self.receipt),
            candidate_publication_authorization=self._candidate_publication_authorization(),
            controller_operator_id="codex",
            expected_git_commit="e" * 40,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
            now=self.now,
        )

    def _controller_binding(self) -> dict[str, str]:
        directories = list(self.archive_root.iterdir())
        self.assertEqual(len(directories), 1)
        path = directories[0] / "attestation.json"
        return {"path": str(path), "sha256": _sha256(path.read_bytes())}

    def _remove_human_decision_checkpoint(self) -> None:
        decision_path = Path(self.human_decision_binding["path"])
        decision_path.chmod(0o600)
        decision_path.unlink()
        decision_path.parent.rmdir()

    @contextmanager
    def _private_root_recovery_context(self) -> Iterator[None]:
        pic_status = self.pic_root.stat()
        archive_status = self.archive_root.stat()
        with self._verification_context(), mock.patch.object(
            publisher, "AUTHORIZED_PIC_ROOT", self.pic_root
        ), mock.patch.object(
            publisher, "AUTHORIZED_PROJECT_HOME_ROOT", self.project_home_root
        ), mock.patch.object(
            publisher,
            "REVIEWED_PRESSURE_GATE_PIC_ROOT_IDENTITY",
            (pic_status.st_dev, pic_status.st_ino),
        ), mock.patch.object(
            publisher,
            "REVIEWED_PRESSURE_GATE_ATTESTATION_ROOT_IDENTITY",
            (archive_status.st_dev, archive_status.st_ino),
        ), mock.patch.object(
            publisher, "REVIEWED_PRESSURE_GATE_ROOT_UID", archive_status.st_uid
        ), mock.patch.object(
            publisher, "REVIEWED_PRESSURE_GATE_ROOT_GID", archive_status.st_gid
        ), mock.patch.object(
            publisher,
            "REVIEWED_PRESSURE_GATE_PIC_ROOT_MODE",
            stat.S_IMODE(pic_status.st_mode),
        ), mock.patch.object(
            publisher,
            "REVIEWED_PRESSURE_GATE_PIC_ROOT_XATTR_BINDINGS",
            {
                name: _sha256(os.getxattr(self.pic_root, name))
                for name in os.listxattr(self.pic_root)
            },
        ), mock.patch.object(
            publisher,
            "REVIEWED_PRESSURE_GATE_ATTESTATION_ROOT_XATTR_BINDINGS",
            {
                name: _sha256(os.getxattr(self.archive_root, name))
                for name in os.listxattr(self.archive_root)
            },
        ):
            yield

    def _guard_path(self) -> Path:
        return self.publication_root / publisher.pilot_publisher._publication_guard_name(
            publisher.CANONICAL_RECEIPT_NAME
        )

    def _guard(self) -> dict[str, object]:
        return json.loads(self._guard_path().read_text(encoding="utf-8"))

    def _seal_path(self) -> Path:
        return self.acceptance_root / publisher._selection_success_seal_name(
            publisher.CANONICAL_RECEIPT_NAME
        )

    def test_publishes_exact_p0_1p00_receipt_and_controller_bound_success_seal(
        self,
    ) -> None:
        with self._verification_context():
            published = self._publish()
            consumed = publisher.consume_published_pressure_selection_receipt(
                self.receipt_path, authorized_pic_root=self.pic_root
            )

        self.assertEqual(published, consumed)
        self.assertEqual(consumed["selected_case"], publisher.REVIEWED_SELECTED_CASE)
        self.assertEqual(consumed["reviewer_id"], publisher.REVIEWED_REVIEWER_ID)
        self.assertEqual(consumed["rationale"], publisher.REVIEWED_RATIONALE)
        self.assertEqual(consumed["qualification_effect"], publisher.QUALIFICATION_EFFECT)
        self.assertEqual(stat.S_IMODE(self.receipt_path.stat().st_mode), 0o444)
        self.assertEqual(self.receipt_path.stat().st_nlink, 1)
        self.assertFalse(self._guard_path().exists())
        self.assertEqual(stat.S_IMODE(self._seal_path().stat().st_mode), 0o444)
        self.assertEqual(self._seal_path().stat().st_nlink, 1)
        controller = consumed["controller_state_attestation"]
        self.assertEqual(
            controller["attestation"]["publisher_source_authentication"],
            self._publisher_source_authentication(),
        )
        attestation_path = Path(controller["binding"]["path"])
        self.assertEqual(stat.S_IMODE(attestation_path.parent.stat().st_mode), 0o500)
        self.assertEqual(
            {path.name for path in attestation_path.parent.iterdir()},
            {"attestation.json", "active_policy.json", "active_promotion.json"},
        )
        for path in attestation_path.parent.iterdir():
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o400)
            self.assertEqual(path.stat().st_nlink, 1)
        state = controller["attestation"]["controller_state"]
        self.assertEqual(state["registered_science_slices"], [])
        self.assertEqual(
            state["frontier_admission_smoke"],
            {"status": publisher.control_plane_common.CLOSED_ADMISSION_SMOKE_STATUS},
        )
        self.assertEqual(state["pending_submission_marker"], "absent")
        self.assertEqual(
            state["manual_accounting_authorizations"],
            self.policy["olcf_side_storage"]["manual_accounting_authorizations"],
        )
        self.assertEqual(state["pending_manual_accounting_marker"], "absent")
        self.assertEqual(state["active_promotion_transaction"], "absent")
        self.assertEqual(state["mirrored_ledger_state"], self._mirrored_ledger_state())
        self.assertEqual(
            state["clean_candidate"],
            {
                "manifest_path": str(self.manifest_path),
                "manifest_sha256": self.policy["science_submission_freeze"][
                    "manifest_sha256"
                ],
                "git_commit": "e" * 40,
                "source_archive_sha256": "f" * 64,
            },
        )
        self.assertEqual(
            consumed["success_seal"]["controller_state_attestation"],
            controller["binding"],
        )

    def test_machine_reanalysis_cannot_create_reviewer_or_candidate_and_human_can_seal(
        self,
    ) -> None:
        self._remove_human_decision_checkpoint()
        self.archive_root.rmdir()
        self.pic_root.chmod(0o2755)
        result = {
            "packet_receipt_sha256": "2" * 64,
            "aggregate_receipt_sha256": "1" * 64,
            "manifest_sha256": "3" * 64,
            "analysis_result_sha256": "4" * 64,
            "status": "pass_engineering_calibration_only",
        }
        raw_cases = [
            {
                "case_id": case_id,
                "descriptor_sha256": f"{index + 5:064x}",
            }
            for index, (case_id, _pressure) in enumerate(
                publisher.pressure_selection.REGISTERED_CASES
            )
        ]
        aggregate_binding = {
            "path": str(self.publication_root / "aggregate-receipt.json"),
            "sha256": result["aggregate_receipt_sha256"],
        }
        packet_binding = {
            "path": str(self.publication_root / "packet-receipt.json"),
            "sha256": result["packet_receipt_sha256"],
        }
        packet = {
            "aggregate_receipt": {"raw_cases": raw_cases},
            "aggregate_bundle": {
                "path": str(self.publication_root / "bundle"),
                "manifest_sha256": result["manifest_sha256"],
            },
            "aggregate_analysis": {
                "path": str(self.publication_root / "analysis.json"),
                "sha256": result["analysis_result_sha256"],
            },
        }
        with self._verification_context(), mock.patch.object(
            publisher, "AUTHORIZED_PIC_ROOT", self.pic_root
        ), mock.patch.object(
            publisher, "AUTHORIZED_PROJECT_HOME_ROOT", self.project_home_root
        ), mock.patch.object(
            publisher.pressure_selection.historical_pressure_pilot_consumer,
            "AUTHORIZED_PRODUCTION_AGGREGATE_RECEIPT_BINDING",
            aggregate_binding,
        ), mock.patch.object(
            publisher.pressure_selection.historical_pressure_pilot_consumer,
            "AUTHORIZED_PRODUCTION_PACKET_RECEIPT_BINDING",
            packet_binding,
        ), mock.patch.object(
            publisher,
            "_execute_candidate_reanalysis",
            return_value=copy.deepcopy(result),
        ), mock.patch.object(
            publisher.pressure_selection.pressure_review_packet_verifier,
            "consume_published_pressure_pilot_review_packet",
            return_value=copy.deepcopy(packet),
        ):
            prepared = publisher.prepare_pressure_reanalysis(
                reanalysis_operator_id="codex",
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                now=self.now - timedelta(minutes=1),
            )
            self.assertNotIn("reviewer_attestation", prepared)
            self.assertNotIn("candidate_pressure_selection_receipt", prepared)
            self.assertFalse(
                (self.pic_root / publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME).exists()
            )
            reanalysis_path = Path(
                prepared["authoritative_reanalysis_attestation"]["path"]
            )
            preparation_path = Path(prepared["stage4_preparation_attestation"]["path"])
            self.assertEqual(stat.S_IMODE(reanalysis_path.stat().st_mode), 0o400)
            self.assertEqual(stat.S_IMODE(reanalysis_path.parent.stat().st_mode), 0o500)
            self.assertEqual(stat.S_IMODE(preparation_path.stat().st_mode), 0o400)
            self.assertEqual(stat.S_IMODE(preparation_path.parent.stat().st_mode), 0o500)
            decision_root = Path(prepared["human_decision_root"])
            self.assertEqual(list(decision_root.iterdir()), [])
            self.assertEqual(stat.S_IMODE(self.archive_root.stat().st_mode), 0o700)
            self.assertEqual(stat.S_IMODE(decision_root.stat().st_mode), 0o700)
            decision_path = decision_root / "dfielding-p0-1p0.json"
            self._write_readonly(
                decision_path,
                _canonical(
                    {
                        "schema_version": 1,
                        "record_type": publisher.HUMAN_DECISION_RECORD_TYPE,
                        "reviewer_id": publisher.REVIEWED_REVIEWER_ID,
                        "reviewed_utc": "2026-06-05T12:02:00Z",
                        "rationale": publisher.REVIEWED_RATIONALE,
                        "reviewer_statement": publisher.HUMAN_DECISION_STATEMENT,
                        "authoritative_reanalysis_attestation": prepared[
                            "authoritative_reanalysis_attestation"
                        ],
                        "stage4_preparation_attestation": prepared[
                            "stage4_preparation_attestation"
                        ],
                        "selected_case": publisher.REVIEWED_SELECTED_CASE,
                    }
                ),
            )
            decision_path.chmod(0o400)
            sealed = publisher.seal_human_pressure_selection(
                human_decision_path=decision_path,
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                now=self.now,
            )

        candidate = Path(sealed["candidate_pressure_selection_receipt"]["path"])
        candidate_authorization = Path(
            sealed["candidate_publication_authorization"]["path"]
        )
        self.assertEqual(
            stat.S_IMODE(
                (self.pic_root / publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME).stat().st_mode
            ),
            0o700,
        )
        self.assertEqual(stat.S_IMODE(candidate.stat().st_mode), 0o400)
        self.assertEqual(stat.S_IMODE(candidate_authorization.stat().st_mode), 0o400)
        receipt = json.loads(candidate.read_text(encoding="utf-8"))
        self.assertEqual(receipt["selected_case"], publisher.REVIEWED_SELECTED_CASE)
        self.assertEqual(
            receipt["authoritative_reanalysis_attestation"],
            prepared["authoritative_reanalysis_attestation"],
        )
        self.assertEqual(
            receipt["reviewer_attestation"], sealed["reviewer_attestation"]
        )
        for key in (
            "authoritative_reanalysis_attestation",
            "reviewer_attestation",
        ):
            path = Path((prepared if key.startswith("authoritative") else sealed)[key]["path"])
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o400)
            self.assertEqual(stat.S_IMODE(path.parent.stat().st_mode), 0o500)
            self.assertEqual([item.name for item in path.parent.iterdir()], ["attestation.json"])
        reviewer = json.loads(
            Path(sealed["reviewer_attestation"]["path"]).read_text(encoding="utf-8")
        )
        self.assertEqual(reviewer["schema_version"], 1)
        self.assertEqual(
            reviewer["record_type"],
            publisher.pressure_selection.pressure_review_packet_verifier.PRESSURE_REVIEWER_RECORD_TYPE,
        )
        self.assertEqual(reviewer["reviewer_id"], publisher.REVIEWED_REVIEWER_ID)
        self.assertEqual(reviewer["rationale"], publisher.REVIEWED_RATIONALE)
        self.assertEqual(reviewer["selected_case"], publisher.REVIEWED_SELECTED_CASE)
        self.assertEqual(
            reviewer["authoritative_reanalysis_attestation"],
            prepared["authoritative_reanalysis_attestation"],
        )
        self.assertEqual(reviewer["reviewed_utc"], "2026-06-05T12:02:00Z")
        self.assertEqual(
            sealed["human_decision"],
            {"path": str(decision_path), "sha256": _sha256(decision_path.read_bytes())},
        )

    def test_existing_inherited_setgid_private_root_is_not_implicitly_recovered(
        self,
    ) -> None:
        candidate_root = self.pic_root / publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME
        self.pic_root.chmod(0o2755)
        pic_descriptor = publisher.pilot_publisher._open_absolute_directory(
            self.pic_root
        )
        try:
            root, descriptor = publisher._open_or_create_private_root(
                self.pic_root,
                pic_descriptor,
                publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME,
                "pressure-selection candidate root",
            )
            os.close(descriptor)
        finally:
            os.close(pic_descriptor)
        self.assertEqual(root, candidate_root)
        retained_identity = (candidate_root.stat().st_dev, candidate_root.stat().st_ino)
        self.assertEqual(stat.S_IMODE(candidate_root.stat().st_mode), 0o700)
        candidate_root.chmod(0o2700)
        pic_descriptor = publisher.pilot_publisher._open_absolute_directory(
            self.pic_root
        )
        try:
            with self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "must be one private directory with mode 0700",
            ):
                publisher._open_or_create_private_root(
                    self.pic_root,
                    pic_descriptor,
                    publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME,
                    "pressure-selection candidate root",
                )
        finally:
            os.close(pic_descriptor)
        self.assertEqual(
            (candidate_root.stat().st_dev, candidate_root.stat().st_ino),
            retained_identity,
        )
        self.assertEqual(stat.S_IMODE(candidate_root.stat().st_mode), 0o2700)

    def test_private_root_creation_occurs_only_after_publication_lock(self) -> None:
        helper_source = inspect.getsource(publisher._open_or_create_private_root)
        self.assertIn("parent_descriptor: int", helper_source)
        self.assertNotIn("_open_absolute_directory(pic_root)", helper_source)
        for function in (
            publisher.prepare_pressure_reanalysis,
            publisher.seal_human_pressure_selection,
        ):
            with self.subTest(function=function.__name__):
                source = inspect.getsource(function)
                self.assertLess(
                    source.index("_lock_publication_transaction"),
                    source.index("_open_or_create_private_root"),
                )

    def test_private_root_creation_is_blocked_by_competing_publication_lock(
        self,
    ) -> None:
        decision_path = Path(self.human_decision_binding["path"])
        decision_payload = decision_path.read_bytes()
        self._remove_human_decision_checkpoint()
        self.archive_root.rmdir()
        self.pic_root.chmod(0o2755)
        anchor = publisher.pilot_publisher._open_absolute_directory(self.root)
        fcntl.flock(anchor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            with self._verification_context(), mock.patch.object(
                publisher, "AUTHORIZED_PIC_ROOT", self.pic_root
            ), mock.patch.object(
                publisher, "AUTHORIZED_PROJECT_HOME_ROOT", self.project_home_root
            ), self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "pressure-selection preparation failed closed",
            ):
                publisher.prepare_pressure_reanalysis(
                    reanalysis_operator_id="codex",
                    expected_git_commit="e" * 40,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    now=self.now,
                )
            self.assertFalse(self.archive_root.exists())
            self.assertFalse(decision_path.parent.exists())

            self.archive_root.mkdir(mode=0o700)
            decision_path.parent.mkdir(mode=0o700)
            decision_path.parent.chmod(0o700)
            self._write_readonly(decision_path, decision_payload)
            decision_path.chmod(0o400)
            with self._verification_context(), mock.patch.object(
                publisher, "AUTHORIZED_PIC_ROOT", self.pic_root
            ), mock.patch.object(
                publisher, "AUTHORIZED_PROJECT_HOME_ROOT", self.project_home_root
            ), self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "human pressure-selection sealing failed closed",
            ):
                publisher.seal_human_pressure_selection(
                    human_decision_path=decision_path,
                    expected_git_commit="e" * 40,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    now=self.now,
                )
            self.assertFalse(
                (self.pic_root / publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME).exists()
            )
        finally:
            fcntl.flock(anchor, fcntl.LOCK_UN)
            os.close(anchor)

    def test_private_root_creation_rejects_writable_pic_parent(self) -> None:
        candidate_root = self.pic_root / publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME
        for mode in (0o2775, 0o2777):
            with self.subTest(mode=f"{mode:04o}"):
                self.pic_root.chmod(mode)
                retained = publisher.pilot_publisher._open_absolute_directory(
                    self.pic_root
                )
                try:
                    with self.assertRaisesRegex(
                        ValueError,
                        "same-account isolated parent",
                    ):
                        publisher._open_or_create_private_root(
                            self.pic_root,
                            retained,
                            publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME,
                            "pressure-selection candidate root",
                        )
                finally:
                    os.close(retained)
                self.assertFalse(candidate_root.exists())

    def test_private_root_creation_strips_inherited_default_acl(self) -> None:
        candidate_root = self.pic_root / publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME
        self.pic_root.chmod(0o2755)
        parent_inode = self.pic_root.stat().st_ino
        xattrs: dict[int, dict[str, bytes]] = {
            parent_inode: {
                "lustre.lov": b"reviewed-layout",
                "system.posix_acl_default": b"parent-default-acl",
            }
        }
        removed: list[str] = []

        def inode_xattrs(descriptor: int) -> dict[str, bytes]:
            inode = os.fstat(descriptor).st_ino
            return xattrs.setdefault(
                inode,
                {
                    "lustre.lov": b"reviewed-layout",
                    "system.posix_acl_access": b"inherited-access-acl",
                    "system.posix_acl_default": b"inherited-default-acl",
                },
            )

        def remove_xattr(descriptor: int, name: str) -> None:
            removed.append(name)
            del inode_xattrs(descriptor)[name]

        retained = publisher.pilot_publisher._open_absolute_directory(self.pic_root)
        try:
            with mock.patch.object(
                publisher.os,
                "listxattr",
                side_effect=lambda descriptor: list(inode_xattrs(descriptor)),
            ), mock.patch.object(
                publisher.os,
                "getxattr",
                side_effect=lambda descriptor, name: inode_xattrs(descriptor)[name],
            ), mock.patch.object(
                publisher.os,
                "removexattr",
                side_effect=remove_xattr,
            ):
                root, descriptor = publisher._open_or_create_private_root(
                    self.pic_root,
                    retained,
                    publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME,
                    "pressure-selection candidate root",
                )
                os.close(descriptor)
        finally:
            os.close(retained)
        self.assertEqual(root, candidate_root)
        self.assertEqual(
            removed,
            ["system.posix_acl_access", "system.posix_acl_default"],
        )
        self.assertEqual(stat.S_IMODE(candidate_root.stat().st_mode), 0o700)

    def test_private_root_creation_rejects_parent_access_acl_before_mutation(
        self,
    ) -> None:
        candidate_root = self.pic_root / publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME
        retained = publisher.pilot_publisher._open_absolute_directory(self.pic_root)
        try:
            with mock.patch.object(
                publisher.os,
                "listxattr",
                return_value=["system.posix_acl_access"],
            ), self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "access ACL xattr",
            ):
                publisher._open_or_create_private_root(
                    self.pic_root,
                    retained,
                    publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME,
                    "pressure-selection candidate root",
                )
        finally:
            os.close(retained)
        self.assertFalse(candidate_root.exists())

    def test_private_root_creation_rejects_storage_xattr_drift(self) -> None:
        self.pic_root.chmod(0o2755)
        parent_inode = self.pic_root.stat().st_ino
        scenarios = (
            ("missing-layout", {}, "xattr bindings differ"),
            (
                "mismatched-layout",
                {"lustre.lov": b"different-layout"},
                "xattr bindings differ",
            ),
            (
                "unexpected-xattr",
                {"lustre.lov": b"reviewed-layout", "user.unexpected": b"value"},
                "unexpected xattrs",
            ),
        )
        for name, child_xattrs, message in scenarios:
            with self.subTest(name=name):
                private_name = f"private-root-{name}"
                retained = publisher.pilot_publisher._open_absolute_directory(
                    self.pic_root
                )

                def inode_xattrs(descriptor: int) -> dict[str, bytes]:
                    if os.fstat(descriptor).st_ino == parent_inode:
                        return {"lustre.lov": b"reviewed-layout"}
                    return child_xattrs

                try:
                    with mock.patch.object(
                        publisher.os,
                        "listxattr",
                        side_effect=lambda descriptor: list(inode_xattrs(descriptor)),
                    ), mock.patch.object(
                        publisher.os,
                        "getxattr",
                        side_effect=lambda descriptor, xattr: inode_xattrs(descriptor)[
                            xattr
                        ],
                    ), self.assertRaisesRegex(
                        publisher.PressureSelectionPublicationError,
                        message,
                    ):
                        publisher._open_or_create_private_root(
                            self.pic_root,
                            retained,
                            private_name,
                            "pressure-selection private root",
                        )
                finally:
                    os.close(retained)
                created = self.pic_root / private_name
                self.assertFalse(created.exists())

    def test_private_root_creation_removes_interrupted_fresh_checkpoint(self) -> None:
        self.pic_root.chmod(0o2755)
        original_parent_fsync = publisher.pilot_publisher._fsync_descriptor
        original_fsync = os.fsync
        original_open = os.open
        original_stat = os.stat
        for failure in (
            "descriptor-open",
            "identity-path-stat",
            "parent-fsync",
            "child-fsync",
        ):
            with self.subTest(failure=failure):
                private_name = f"private-root-{failure}"
                private_root = self.pic_root / private_name
                retained = publisher.pilot_publisher._open_absolute_directory(
                    self.pic_root
                )
                failed = False

                def fail_parent_fsync(descriptor: int) -> None:
                    nonlocal failed
                    if not failed:
                        failed = True
                        raise OSError("injected parent sync failure")
                    original_parent_fsync(descriptor)

                def fail_child_fsync(descriptor: int) -> None:
                    nonlocal failed
                    if descriptor != retained and not failed:
                        failed = True
                        raise OSError("injected child sync failure")
                    original_fsync(descriptor)

                def fail_identity_path_stat(
                    path: object, *args: object, **kwargs: object
                ) -> os.stat_result:
                    nonlocal failed
                    if (
                        path == private_name
                        and kwargs.get("dir_fd") == retained
                        and not failed
                    ):
                        failed = True
                        raise OSError("injected identity path-stat failure")
                    return original_stat(path, *args, **kwargs)

                def fail_descriptor_open(
                    path: object, flags: int, *args: object, **kwargs: object
                ) -> int:
                    nonlocal failed
                    if (
                        path == private_name
                        and kwargs.get("dir_fd") == retained
                        and not failed
                    ):
                        failed = True
                        raise OSError("injected descriptor open failure")
                    return original_open(path, flags, *args, **kwargs)

                patcher = (
                    mock.patch.object(
                        publisher.pilot_publisher,
                        "_fsync_descriptor",
                        side_effect=fail_parent_fsync,
                    )
                    if failure == "parent-fsync"
                    else (
                        mock.patch.object(
                            publisher.os,
                            "fsync",
                            side_effect=fail_child_fsync,
                        )
                        if failure == "child-fsync"
                        else mock.patch.object(
                            publisher.os,
                            "open" if failure == "descriptor-open" else "stat",
                            side_effect=(
                                fail_descriptor_open
                                if failure == "descriptor-open"
                                else fail_identity_path_stat
                            ),
                        )
                    )
                )
                try:
                    with patcher, self.assertRaises((OSError, ValueError)):
                        publisher._open_or_create_private_root(
                            self.pic_root,
                            retained,
                            private_name,
                            "pressure-selection private root",
                        )
                finally:
                    os.close(retained)
                self.assertTrue(failed)
                self.assertFalse(private_root.exists())

    def test_private_root_creation_rejects_replaced_pic_path(self) -> None:
        retained = publisher.pilot_publisher._open_absolute_directory(self.pic_root)
        displaced = self.root / "displaced-pic"
        self.pic_root.rename(displaced)
        self.pic_root.mkdir()
        try:
            with self.assertRaisesRegex(
                ValueError,
                "authorized PIC root changed",
            ):
                publisher._open_or_create_private_root(
                    self.pic_root,
                    retained,
                    publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME,
                    "pressure-selection candidate root",
                )
        finally:
            os.close(retained)
        self.assertFalse(
            (self.pic_root / publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME).exists()
        )

    def test_private_root_binding_rejects_replaced_child_path(self) -> None:
        retained = publisher.pilot_publisher._open_absolute_directory(self.pic_root)
        displaced = self.pic_root / "displaced-candidate-root"
        root, descriptor = publisher._open_or_create_private_root(
            self.pic_root,
            retained,
            publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME,
            "pressure-selection candidate root",
        )
        root.rename(displaced)
        root.mkdir(mode=0o700)
        try:
            with self.assertRaisesRegex(
                ValueError,
                "changed during publication",
            ):
                publisher._require_private_root_binding(
                    self.pic_root,
                    retained,
                    publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME,
                    descriptor,
                    "pressure-selection candidate root",
                )
        finally:
            os.close(descriptor)
            os.close(retained)

    def test_preparation_and_sealing_close_descriptors_when_pic_open_fails(
        self,
    ) -> None:
        original_open = publisher.pilot_publisher._open_absolute_directory
        decision_root = Path(self.human_decision_binding["path"]).parent
        operations = (
            (
                "prepare",
                lambda: publisher.prepare_pressure_reanalysis(
                    reanalysis_operator_id="codex",
                    expected_git_commit="e" * 40,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    now=self.now,
                ),
            ),
            (
                "seal",
                lambda: publisher.seal_human_pressure_selection(
                    human_decision_path=Path(self.human_decision_binding["path"]),
                    expected_git_commit="e" * 40,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    now=self.now,
                ),
            ),
        )
        for label, operation in operations:
            with self.subTest(operation=label):
                retained: list[int] = []

                def fail_pic_open(path: Path) -> int:
                    if Path(path) == self.pic_root:
                        raise OSError("injected PIC-root open failure")
                    descriptor = original_open(path)
                    if Path(path) != decision_root:
                        retained.append(descriptor)
                    return descriptor

                with self._verification_context(), mock.patch.object(
                    publisher, "AUTHORIZED_PIC_ROOT", self.pic_root
                ), mock.patch.object(
                    publisher, "AUTHORIZED_PROJECT_HOME_ROOT", self.project_home_root
                ), mock.patch.object(
                    publisher.pilot_publisher,
                    "_open_absolute_directory",
                    side_effect=fail_pic_open,
                ), self.assertRaisesRegex(
                    publisher.PressureSelectionPublicationError,
                    "failed closed",
                ):
                    operation()
                self.assertEqual(len(retained), 2)
                for descriptor in retained:
                    with self.assertRaises(OSError):
                        os.fstat(descriptor)

    def test_explicit_private_root_recovery_and_reconciliation_preserve_inode(
        self,
    ) -> None:
        self._remove_human_decision_checkpoint()
        self.archive_root.chmod(0o2700)
        retained_identity = (
            self.archive_root.stat().st_dev,
            self.archive_root.stat().st_ino,
        )
        with self._private_root_recovery_context():
            recovered = publisher.recover_pressure_gate_attestation_root(
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
            with self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "not the exact recoverable checkpoint",
            ):
                publisher.recover_pressure_gate_attestation_root(
                    expected_git_commit="e" * 40,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
            reconciled = publisher.recover_pressure_gate_attestation_root(
                expected_git_commit="e" * 40,
                reconcile_exact_empty_normalized_root=True,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

        self.assertEqual(
            recovered["action"], "recovered_exact_empty_inherited_setgid_root"
        )
        self.assertEqual(
            reconciled["action"], "reconciled_exact_empty_normalized_root"
        )
        self.assertEqual(recovered["authority"], "none")
        self.assertEqual(
            recovered["qualification_effect"],
            publisher.PRIVATE_ROOT_RECOVERY_QUALIFICATION_EFFECT,
        )
        self.assertEqual(recovered["before"]["mode"], "2700")
        self.assertEqual(recovered["after"]["mode"], "0700")
        self.assertEqual(
            (self.archive_root.stat().st_dev, self.archive_root.stat().st_ino),
            retained_identity,
        )
        self.assertEqual(stat.S_IMODE(self.archive_root.stat().st_mode), 0o700)
        self.assertEqual(list(self.archive_root.iterdir()), [])
        self.assertFalse(
            (self.pic_root / publisher.PRESSURE_GATE_HUMAN_DECISION_ROOT_NAME).exists()
        )
        self.assertFalse(
            (self.pic_root / publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME).exists()
        )
        self.assertFalse(self.receipt_path.exists())

    def test_explicit_private_root_recovery_rejects_nonempty_checkpoint(self) -> None:
        self._remove_human_decision_checkpoint()
        marker = self.archive_root / "unexpected"
        marker.write_text("not recoverable\n", encoding="utf-8")
        marker.chmod(0o400)
        self.archive_root.chmod(0o2700)
        with self._private_root_recovery_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "not the exact recoverable checkpoint",
        ):
            publisher.recover_pressure_gate_attestation_root(
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(stat.S_IMODE(self.archive_root.stat().st_mode), 0o2700)
        self.assertTrue(marker.is_file())

    def test_explicit_private_root_recovery_rejects_namespace_blocker(self) -> None:
        self._remove_human_decision_checkpoint()
        self.archive_root.chmod(0o2700)
        blockers = (
            (
                "human decision root",
                self.pic_root / publisher.PRESSURE_GATE_HUMAN_DECISION_ROOT_NAME,
                True,
            ),
            (
                "candidate root",
                self.pic_root / publisher.PRESSURE_GATE_CANDIDATE_ROOT_NAME,
                True,
            ),
            ("canonical receipt", self.receipt_path, False),
            ("publication guard", self._guard_path(), False),
            ("durable success seal", self._seal_path(), False),
        )
        for label, blocker, directory in blockers:
            with self.subTest(blocker=label):
                if directory:
                    blocker.mkdir(mode=0o700)
                else:
                    blocker.write_text("authority blocker\n", encoding="utf-8")
                    blocker.chmod(0o400)
                try:
                    with self._private_root_recovery_context(), self.assertRaisesRegex(
                        publisher.PressureSelectionPublicationError,
                        "pressure-gate private-root recovery failed closed",
                    ):
                        publisher.recover_pressure_gate_attestation_root(
                            expected_git_commit="e" * 40,
                            authorized_pic_root=self.pic_root,
                            authorized_project_home_root=self.project_home_root,
                        )
                    self.assertEqual(
                        stat.S_IMODE(self.archive_root.stat().st_mode), 0o2700
                    )
                finally:
                    if directory:
                        blocker.rmdir()
                    else:
                        blocker.chmod(0o600)
                        blocker.unlink()

    def test_explicit_private_root_recovery_rejects_competing_lock(self) -> None:
        self._remove_human_decision_checkpoint()
        self.archive_root.chmod(0o2700)
        anchor = publisher.pilot_publisher._open_absolute_directory(self.root)
        fcntl.flock(anchor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            with self._private_root_recovery_context(), self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "private-root recovery failed closed",
            ):
                publisher.recover_pressure_gate_attestation_root(
                    expected_git_commit="e" * 40,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        finally:
            fcntl.flock(anchor, fcntl.LOCK_UN)
            os.close(anchor)
        self.assertEqual(stat.S_IMODE(self.archive_root.stat().st_mode), 0o2700)

    def test_explicit_private_root_recovery_rejects_controller_state_drift(
        self,
    ) -> None:
        self._remove_human_decision_checkpoint()
        self.archive_root.chmod(0o2700)
        original_capture = publisher._capture_live_controller_state
        capture_count = 0

        def drift_second_capture(**kwargs: object) -> dict[str, object]:
            nonlocal capture_count
            capture_count += 1
            result = original_capture(**kwargs)
            if capture_count == 2:
                result = copy.deepcopy(result)
                result["controller_state"]["pending_submission_marker"] = "drifted"
            return result

        with self._private_root_recovery_context(), mock.patch.object(
            publisher,
            "_capture_live_controller_state",
            side_effect=drift_second_capture,
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "controller state changed",
        ):
            publisher.recover_pressure_gate_attestation_root(
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(stat.S_IMODE(self.archive_root.stat().st_mode), 0o700)

    def test_private_root_xattr_bindings_hash_values(self) -> None:
        with mock.patch.object(
            publisher.os, "listxattr", return_value=["lustre.lov"]
        ), mock.patch.object(
            publisher.os, "getxattr", return_value=b"reviewed-layout"
        ) as getxattr:
            bindings = publisher._private_root_xattr_bindings(
                -1, "pressure-gate attestation root"
            )
        self.assertEqual(bindings, {"lustre.lov": _sha256(b"reviewed-layout")})
        getxattr.assert_called_once_with(-1, "lustre.lov")

    def test_private_root_recovery_rejects_final_path_replacement(self) -> None:
        self._remove_human_decision_checkpoint()
        self.archive_root.chmod(0o2700)
        retained_identity = (
            self.archive_root.stat().st_dev,
            self.archive_root.stat().st_ino,
        )
        displaced = self.pic_root / "displaced-pressure-gate-attestations"
        original_capture = publisher._capture_live_controller_state
        capture_count = 0

        def replace_during_second_capture(**kwargs: object) -> dict[str, object]:
            nonlocal capture_count
            capture_count += 1
            result = original_capture(**kwargs)
            if capture_count == 2:
                self.archive_root.rename(displaced)
                self.archive_root.mkdir(mode=0o700)
            return result

        with self._private_root_recovery_context(), mock.patch.object(
            publisher,
            "_capture_live_controller_state",
            side_effect=replace_during_second_capture,
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "private-root recovery failed closed",
        ):
            publisher.recover_pressure_gate_attestation_root(
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertNotEqual(
            (self.archive_root.stat().st_dev, self.archive_root.stat().st_ino),
            retained_identity,
        )
        self.assertEqual(
            (displaced.stat().st_dev, displaced.stat().st_ino),
            retained_identity,
        )

    def test_private_root_recovery_rejects_final_pic_root_mode_drift(self) -> None:
        self._remove_human_decision_checkpoint()
        self.archive_root.chmod(0o2700)
        original_capture = publisher._capture_live_controller_state
        capture_count = 0

        def drift_during_second_capture(**kwargs: object) -> dict[str, object]:
            nonlocal capture_count
            capture_count += 1
            result = original_capture(**kwargs)
            if capture_count == 2:
                self.pic_root.chmod(0o2777)
            return result

        with self._private_root_recovery_context(), mock.patch.object(
            publisher,
            "_capture_live_controller_state",
            side_effect=drift_during_second_capture,
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "reviewed private-root checkpoint",
        ):
            publisher.recover_pressure_gate_attestation_root(
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_private_root_recovery_rebinds_path_after_final_state_read(self) -> None:
        self._remove_human_decision_checkpoint()
        self.archive_root.chmod(0o2700)
        retained_identity = (
            self.archive_root.stat().st_dev,
            self.archive_root.stat().st_ino,
        )
        displaced = self.pic_root / "displaced-after-final-state"
        original_state = publisher._private_root_state
        state_count = 0

        def replace_after_final_state(
            descriptor: int, label: str
        ) -> dict[str, object]:
            nonlocal state_count
            state_count += 1
            result = original_state(descriptor, label)
            if state_count == 3:
                self.archive_root.rename(displaced)
                self.archive_root.mkdir(mode=0o700)
            return result

        with self._private_root_recovery_context(), mock.patch.object(
            publisher,
            "_private_root_state",
            side_effect=replace_after_final_state,
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "private-root recovery failed closed",
        ):
            publisher.recover_pressure_gate_attestation_root(
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(
            (displaced.stat().st_dev, displaced.stat().st_ino),
            retained_identity,
        )

    def test_interrupted_private_root_recovery_requires_explicit_reconciliation(
        self,
    ) -> None:
        self._remove_human_decision_checkpoint()
        self.archive_root.chmod(0o2700)
        with self._private_root_recovery_context(), mock.patch.object(
            publisher.pilot_publisher,
            "_fsync_descriptor",
            side_effect=OSError("injected post-normalization sync failure"),
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "private-root recovery failed closed",
        ):
            publisher.recover_pressure_gate_attestation_root(
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(stat.S_IMODE(self.archive_root.stat().st_mode), 0o700)
        with self._private_root_recovery_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "not the exact recoverable checkpoint",
        ):
            publisher.recover_pressure_gate_attestation_root(
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with self._private_root_recovery_context():
            reconciled = publisher.recover_pressure_gate_attestation_root(
                expected_git_commit="e" * 40,
                reconcile_exact_empty_normalized_root=True,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(
            reconciled["action"], "reconciled_exact_empty_normalized_root"
        )

    def test_private_root_recovery_cli_never_dispatches_preparation(self) -> None:
        result = {
            "schema_version": 1,
            "record_type": publisher.PRIVATE_ROOT_RECOVERY_RECORD_TYPE,
            "action": "recovered_exact_empty_inherited_setgid_root",
        }
        for command, reconcile in (
            ("recover-pressure-gate-attestation-root", False),
            ("reconcile-pressure-gate-attestation-root", True),
        ):
            with self.subTest(command=command), mock.patch.object(
                publisher,
                "recover_pressure_gate_attestation_root",
                return_value=result,
            ) as recovery, mock.patch.object(
                publisher, "prepare_pressure_reanalysis"
            ) as preparation, redirect_stdout(io.StringIO()) as stdout:
                status = publisher.main(
                    [
                        command,
                        "--expected-git-commit",
                        "e" * 40,
                        "--authorized-pic-root",
                        str(self.pic_root),
                        "--authorized-project-home-root",
                        str(self.project_home_root),
                    ]
                )
            self.assertEqual(status, 0)
            preparation.assert_not_called()
            recovery.assert_called_once_with(
                expected_git_commit="e" * 40,
                reconcile_exact_empty_normalized_root=reconcile,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
            self.assertEqual(json.loads(stdout.getvalue()), result)

    def test_private_root_recovery_cli_enforces_exact_mode_matrix(self) -> None:
        self._remove_human_decision_checkpoint()
        self.archive_root.chmod(0o2700)
        arguments = [
            "--expected-git-commit",
            "e" * 40,
            "--authorized-pic-root",
            str(self.pic_root),
            "--authorized-project-home-root",
            str(self.project_home_root),
        ]
        cases = (
            (
                "recover-pressure-gate-attestation-root",
                0o2700,
                0,
                "recovered_exact_empty_inherited_setgid_root",
                0o700,
            ),
            (
                "recover-pressure-gate-attestation-root",
                0o700,
                2,
                None,
                0o700,
            ),
            (
                "reconcile-pressure-gate-attestation-root",
                0o700,
                0,
                "reconciled_exact_empty_normalized_root",
                0o700,
            ),
            (
                "reconcile-pressure-gate-attestation-root",
                0o2700,
                2,
                None,
                0o2700,
            ),
        )
        for command, initial_mode, status, action, final_mode in cases:
            with self.subTest(command=command, initial_mode=f"{initial_mode:04o}"):
                self.archive_root.chmod(initial_mode)
                stdout = io.StringIO()
                stderr = io.StringIO()
                with self._private_root_recovery_context(), redirect_stdout(
                    stdout
                ), redirect_stderr(stderr):
                    observed = publisher.main([command, *arguments])
                self.assertEqual(observed, status)
                self.assertEqual(stat.S_IMODE(self.archive_root.stat().st_mode), final_mode)
                if action is None:
                    self.assertEqual(stdout.getvalue(), "")
                    failure = json.loads(stderr.getvalue())
                    self.assertEqual(failure["status"], "failed_closed")
                    self.assertEqual(failure["command"], command)
                else:
                    self.assertEqual(stderr.getvalue(), "")
                    self.assertEqual(json.loads(stdout.getvalue())["action"], action)

    def test_private_root_recovery_cli_wraps_root_resolution_failure(self) -> None:
        stderr = io.StringIO()
        with self._verification_context(), mock.patch.object(
            publisher, "AUTHORIZED_PIC_ROOT", self.pic_root
        ), mock.patch.object(
            publisher, "AUTHORIZED_PROJECT_HOME_ROOT", self.project_home_root
        ), mock.patch.object(
            publisher.pilot_publisher,
            "_publication_root",
            side_effect=publisher.pilot_publisher.PressurePilotPublicationError(
                "injected root resolution failure"
            ),
        ), redirect_stderr(stderr):
            status = publisher.main(
                [
                    "recover-pressure-gate-attestation-root",
                    "--expected-git-commit",
                    "e" * 40,
                    "--authorized-pic-root",
                    str(self.pic_root),
                    "--authorized-project-home-root",
                    str(self.project_home_root),
                ]
            )
        self.assertEqual(status, 2)
        failure = json.loads(stderr.getvalue())
        self.assertEqual(failure["status"], "failed_closed")
        self.assertEqual(
            failure["message"], "pressure-gate private-root recovery failed closed"
        )

    def test_human_decision_must_be_strictly_later_than_reanalysis(self) -> None:
        decision = {
            "schema_version": 1,
            "record_type": publisher.HUMAN_DECISION_RECORD_TYPE,
            "reviewer_id": publisher.REVIEWED_REVIEWER_ID,
            "reviewed_utc": "2026-06-05T12:01:00Z",
            "rationale": publisher.REVIEWED_RATIONALE,
            "reviewer_statement": publisher.HUMAN_DECISION_STATEMENT,
            "authoritative_reanalysis_attestation": self.reanalysis_binding,
            "stage4_preparation_attestation": self.stage4_preparation_binding,
            "selected_case": publisher.REVIEWED_SELECTED_CASE,
        }
        with self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "strictly later than sealed reanalysis",
        ):
            publisher._validate_human_decision(
                decision,
                reanalysis_verification=self._reanalysis(self.receipt),
            )

    def test_live_verification_rejects_state_different_from_sealed_attestation(
        self,
    ) -> None:
        verified = {
            "controller_state_attestation": {
                "attestation": {"controller_state": {"state": "sealed"}}
            }
        }
        capture = {"controller_state": {"state": "changed"}}
        with mock.patch.object(
            publisher,
            "verify_published_pressure_selection_receipt",
            return_value=verified,
        ), mock.patch.object(
            publisher,
            "_capture_live_controller_state",
            return_value=capture,
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "live controller state differs",
        ):
            publisher.verify_published_pressure_selection_live_state(
                self.receipt_path,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_mutating_cli_source_authentication_requires_snapshot_environment(
        self,
    ) -> None:
        with mock.patch.dict(os.environ, {}, clear=True), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "requires an authenticated source snapshot",
        ):
            publisher._runtime_publisher_source_authentication("e" * 40)

    def test_runtime_source_authentication_rejects_archive_before_parsing(self) -> None:
        archive_path = self.root / "untrusted-runtime-source.tar"
        self._write_readonly(archive_path, b"attacker archive!!!")
        executing_root = Path(publisher.__file__).resolve().parents[2]
        environment = {
            "PIC_PRESSURE_PUBLICATION_SOURCE_ARCHIVE_PATH": str(archive_path),
            publisher.pilot_publisher.WORKER_SOURCE_SNAPSHOT_ROOT_ENV: str(
                executing_root
            ),
        }
        with mock.patch.dict(os.environ, environment, clear=True), mock.patch.object(
            publisher.pilot_publisher,
            "TRUSTED_SOURCE_REPOSITORY",
            self.root,
        ), mock.patch.object(
            publisher,
            "_trusted_publisher_source_archive",
            return_value=("e" * 40, b"trusted Git archive"),
        ), mock.patch.object(
            publisher.control_plane_common,
            "_source_archive_regular_files",
        ) as parser, self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "differs from the expected commit",
        ):
            publisher._runtime_publisher_source_authentication("e" * 40)
        parser.assert_not_called()

    def test_trusted_stage4_source_archive_excludes_ambient_git_poison(self) -> None:
        commit = "e" * 40
        archive = b"trusted Stage-4 archive"
        poison = {
            "GIT_DIR": "/tmp/poison.git",
            "GIT_WORK_TREE": "/tmp/poison-work-tree",
            "GIT_CONFIG_GLOBAL": "/tmp/poison-global.gitconfig",
            "GIT_CONFIG_SYSTEM": "/tmp/poison-system.gitconfig",
            "GIT_CONFIG_NOSYSTEM": "0",
            "GIT_OBJECT_DIRECTORY": "/tmp/poison-objects",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES": "/tmp/poison-alternates",
        }
        expected_environment = publisher.control_plane_common.trusted_git_environment()
        expected_environment["GIT_NO_REPLACE_OBJECTS"] = "1"
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(
            publisher.pilot_publisher, "TRUSTED_SOURCE_REPOSITORY", Path(directory)
        ), mock.patch.dict(os.environ, poison, clear=True), mock.patch.object(
            publisher.pilot_publisher.subprocess,
            "run",
            side_effect=[
                publisher.pilot_publisher.subprocess.CompletedProcess(
                    [], 0, commit.encode("ascii") + b"\n", b""
                ),
                publisher.pilot_publisher.subprocess.CompletedProcess(
                    [], 0, archive, b""
                ),
            ],
        ) as tracked:
            repository = str(Path(directory).resolve())
            self.assertEqual(
                publisher._trusted_publisher_source_archive(commit),
                (commit, archive),
            )
        expected_arguments = [
            publisher.control_plane_common.trusted_git_command(
                "--no-replace-objects",
                "-C",
                repository,
                "rev-parse",
                "--verify",
                f"{commit}^{{commit}}",
            ),
            publisher.control_plane_common.trusted_git_command(
                "--no-replace-objects",
                "-C",
                repository,
                "archive",
                "--format=tar",
                commit,
            ),
        ]
        self.assertEqual(len(tracked.call_args_list), 2)
        for call, arguments in zip(tracked.call_args_list, expected_arguments):
            self.assertEqual(call.args[0], arguments)
            self.assertEqual(call.kwargs["env"], expected_environment)
            for name, value in poison.items():
                self.assertNotEqual(call.kwargs["env"].get(name), value)

    def test_trusted_stage4_source_archive_ignores_replacement_refs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository = Path(directory)

            def git(*args: str) -> str:
                return subprocess.check_output(
                    ["/usr/bin/git", *args],
                    cwd=repository,
                    text=True,
                    env={"HOME": "/", "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin"},
                ).strip()

            git("init", "-q")
            tracked = repository / "tracked.txt"
            tracked.write_text("reviewed\n", encoding="utf-8")
            git("add", "tracked.txt")
            git(
                "-c",
                "user.name=Stage4",
                "-c",
                "user.email=stage4@example.invalid",
                "commit",
                "-qm",
                "reviewed",
            )
            reviewed_commit = git("rev-parse", "HEAD")
            tracked.write_text("replacement\n", encoding="utf-8")
            git("add", "tracked.txt")
            git(
                "-c",
                "user.name=Stage4",
                "-c",
                "user.email=stage4@example.invalid",
                "commit",
                "-qm",
                "replacement",
            )
            replacement_commit = git("rev-parse", "HEAD")
            git("replace", reviewed_commit, replacement_commit)

            with mock.patch.object(
                publisher.pilot_publisher,
                "TRUSTED_SOURCE_REPOSITORY",
                repository,
            ):
                commit, archive_payload = publisher._trusted_publisher_source_archive(
                    reviewed_commit
                )
            self.assertEqual(commit, reviewed_commit)
            with tarfile.open(fileobj=io.BytesIO(archive_payload), mode="r:") as archive:
                member = archive.extractfile("tracked.txt")
                self.assertIsNotNone(member)
                assert member is not None
                self.assertEqual(member.read(), b"reviewed\n")

    def test_rejects_choice_reviewer_and_rationale_drift(self) -> None:
        wrong_choice = copy.deepcopy(self.receipt)
        wrong_choice["selected_case"] = {
            "case_id": "ps_p0_0p10",
            "problem_ps_p0": 0.1,
        }
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "choice differs from the reviewed p0=1.0 baseline",
        ):
            publisher.publish_pressure_selection(
                wrong_choice,
                candidate_publication_authorization=self._candidate_publication_authorization(
                    wrong_choice
                ),
                controller_operator_id="codex",
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                now=self.now,
            )
        self.assertFalse(self.receipt_path.exists())

        self.reviewer_id = "another-reviewer"
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "reviewer differs from the reviewed human",
        ):
            self._publish()
        self.reviewer_id = publisher.REVIEWED_REVIEWER_ID

        self.rationale = "A different pressure-selection rationale."
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "rationale differs from the reviewed rationale",
        ):
            self._publish()
        self.assertFalse(self.receipt_path.exists())
        self.assertEqual(list(self.archive_root.iterdir()), [])

    def test_rejects_active_transaction_pending_submission_and_nonempty_slices(
        self,
    ) -> None:
        transaction = (
            self.pic_root
            / publisher.control_plane_common.ACTIVE_PROMOTION_TRANSACTION_RELATIVE
        )
        self._write_readonly(transaction, b"{}\n")
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "active-policy promotion transaction blocks pressure selection",
        ):
            self._publish()
        transaction.chmod(0o600)
        transaction.unlink()

        pending = self.pic_root / "ledger" / "pending_submission.json"
        pending.write_text("{}\n", encoding="utf-8")
        pending.chmod(0o444)
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "pending scheduler submission blocks pressure selection",
        ):
            self._publish()
        pending.chmod(0o600)
        pending.unlink()

        pending_accounting = self.pic_root / "ledger" / "pending_manual_accounting.json"
        pending_accounting.write_text("{}\n", encoding="utf-8")
        pending_accounting.chmod(0o444)
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "incomplete manual accounting blocks pressure selection",
        ):
            self._publish()
        pending_accounting.chmod(0o600)
        pending_accounting.unlink()

        self.policy["frontier_admission_smoke"] = {
            "status": publisher.control_plane_common.AUTHORIZED_ADMISSION_SMOKE_STATUS
        }
        self._write_active_state()
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "closed admission-smoke authority",
        ):
            self._publish()
        self.policy["frontier_admission_smoke"] = {
            "status": publisher.control_plane_common.CLOSED_ADMISSION_SMOKE_STATUS
        }
        self._write_active_state()

        self.policy["registered_science_slices"] = [{"authorization_id": "forbidden"}]
        self._write_active_state()
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "empty registered-science allowlist",
        ):
            self._publish()
        self.assertFalse(self.receipt_path.exists())
        self.assertEqual(list(self.archive_root.iterdir()), [])

    def test_rejects_reanalysis_source_different_from_clean_candidate(self) -> None:
        self.reanalysis_source_archive_sha256 = "0" * 64
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "reanalysis source differs from the authorized clean candidate",
        ):
            self._publish()
        self.assertFalse(self.receipt_path.exists())
        self.assertEqual(list(self.archive_root.iterdir()), [])

    def test_success_seal_failure_retains_guard_and_reconciles_exact_publication(
        self,
    ) -> None:
        with self._verification_context(), mock.patch.object(
            publisher,
            "_publish_selection_success_seal_at",
            side_effect=OSError("injected success-seal failure"),
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "retained under fail-closed guard",
        ):
            self._publish()

        self.assertTrue(self.receipt_path.is_file())
        self.assertTrue(self._guard_path().is_file())
        self.assertFalse(self._seal_path().exists())
        controller_binding = self._controller_binding()
        self.assertEqual(
            self._guard()["controller_state_attestation"], controller_binding
        )
        self.assertEqual(
            self._guard()["receipt"],
            {"path": str(self.receipt_path), "sha256": _sha256(self.receipt_path.read_bytes())},
        )

        with self._verification_context():
            reconciled = publisher.reconcile_pressure_selection_publication(
                self.receipt_path,
                controller_state_attestation=controller_binding,
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                now=self.now,
            )
        self.assertEqual(reconciled["selected_case"], publisher.REVIEWED_SELECTED_CASE)
        self.assertFalse(self._guard_path().exists())
        self.assertTrue(self._seal_path().is_file())

    def test_reconciliation_rejects_runtime_source_different_from_published_controller(
        self,
    ) -> None:
        with self._verification_context():
            self._publish()
            controller_binding = self._controller_binding()
            different_source = self._publisher_source_authentication()
            different_source["archive_sha256"] = "0" * 64
            with mock.patch.object(
                publisher,
                "_runtime_publisher_source_authentication",
                return_value=different_source,
            ), self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "differs from the runtime publisher source",
            ):
                publisher.reconcile_pressure_selection_publication(
                    self.receipt_path,
                    controller_state_attestation=controller_binding,
                    expected_git_commit="e" * 40,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    now=self.now,
                )

    def test_reconciliation_rejects_runtime_source_before_sealing_guarded_receipt(
        self,
    ) -> None:
        with self._verification_context(), mock.patch.object(
            publisher,
            "_publish_selection_success_seal_at",
            side_effect=OSError("injected success-seal failure"),
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "retained under fail-closed guard",
        ):
            self._publish()
        controller_binding = self._controller_binding()
        different_source = self._publisher_source_authentication()
        different_source["archive_sha256"] = "0" * 64

        with self._verification_context(), mock.patch.object(
            publisher,
            "_runtime_publisher_source_authentication",
            return_value=different_source,
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "differs from the runtime publisher source",
        ):
            publisher.reconcile_pressure_selection_publication(
                self.receipt_path,
                controller_state_attestation=controller_binding,
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                now=self.now,
            )

        self.assertTrue(self.receipt_path.is_file())
        self.assertTrue(self._guard_path().is_file())
        self.assertFalse(self._seal_path().exists())

    def test_post_commit_guard_disarm_wrapper_failure_reconciles_as_success(self) -> None:
        original = publisher.pilot_publisher._disarm_publication_guard_at

        def disarm_then_raise(*args: object, **kwargs: object) -> None:
            original(*args, **kwargs)
            raise OSError("injected post-commit guard-disarm wrapper failure")

        with self._verification_context(), mock.patch.object(
            publisher.pilot_publisher,
            "_disarm_publication_guard_at",
            side_effect=disarm_then_raise,
        ):
            result = self._publish()

        self.assertEqual(result["selected_case"], publisher.REVIEWED_SELECTED_CASE)
        self.assertFalse(self._guard_path().exists())
        self.assertTrue(self._seal_path().is_file())

    def test_publish_retains_guard_when_live_state_changes_after_success_seal(
        self,
    ) -> None:
        with self._verification_context():
            baseline = publisher._capture_live_controller_state(
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
            changed = copy.deepcopy(baseline)
            changed["controller_state"]["pending_submission_marker"] = "changed"
            with mock.patch.object(
                publisher,
                "_capture_live_controller_state",
                side_effect=[baseline, baseline, changed],
            ), self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "retained under fail-closed guard",
            ):
                self._publish()

        self.assertTrue(self.receipt_path.is_file())
        self.assertTrue(self._seal_path().is_file())
        self.assertTrue(self._guard_path().is_file())

    def test_reconciliation_retains_guard_when_state_changes_before_guard_removal(
        self,
    ) -> None:
        with self._verification_context(), mock.patch.object(
            publisher,
            "_publish_selection_success_seal_at",
            side_effect=OSError("injected success-seal failure"),
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "retained under fail-closed guard",
        ):
            self._publish()
        controller_binding = self._controller_binding()

        with self._verification_context():
            baseline = publisher._capture_live_controller_state(
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
            changed = copy.deepcopy(baseline)
            changed["controller_state"]["pending_submission_marker"] = "changed"
            with mock.patch.object(
                publisher,
                "_capture_live_controller_state",
                side_effect=[baseline, baseline, changed],
            ), self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "live controller state differs",
            ):
                publisher.reconcile_pressure_selection_publication(
                    self.receipt_path,
                    controller_state_attestation=controller_binding,
                    expected_git_commit="e" * 40,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    now=self.now,
                )

        self.assertTrue(self.receipt_path.is_file())
        self.assertTrue(self._seal_path().is_file())
        self.assertTrue(self._guard_path().is_file())

    def test_guard_arm_wrapper_failure_retains_exact_fail_closed_guard(self) -> None:
        original = publisher._arm_selection_recovery_guard_at

        def arm_then_raise(*args: object, **kwargs: object) -> object:
            original(*args, **kwargs)
            raise OSError("injected post-commit guard-arm wrapper failure")

        with self._verification_context(), mock.patch.object(
            publisher,
            "_arm_selection_recovery_guard_at",
            side_effect=arm_then_raise,
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "retained under fail-closed guard",
        ):
            self._publish()

        self.assertTrue(self._guard_path().is_file())
        self.assertFalse(self.receipt_path.exists())
        self.assertFalse(self._seal_path().exists())
        self.assertEqual(len(list(self.archive_root.iterdir())), 1)
        self.assertEqual(
            self._guard()["controller_state_attestation"], self._controller_binding()
        )
        with self._verification_context():
            recovered = publisher.reconcile_pressure_selection_publication(
                self.receipt_path,
                controller_state_attestation=self._controller_binding(),
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                now=self.now,
            )
        self.assertEqual(
            recovered["status"], "no_publication_exposed_guard_cleared"
        )
        self.assertFalse(self._guard_path().exists())

    def test_reconciliation_rejects_runtime_source_before_clearing_guard_only(
        self,
    ) -> None:
        original = publisher._arm_selection_recovery_guard_at

        def arm_then_raise(*args: object, **kwargs: object) -> object:
            original(*args, **kwargs)
            raise OSError("injected post-commit guard-arm wrapper failure")

        with self._verification_context(), mock.patch.object(
            publisher,
            "_arm_selection_recovery_guard_at",
            side_effect=arm_then_raise,
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "retained under fail-closed guard",
        ):
            self._publish()
        controller_binding = self._controller_binding()
        different_source = self._publisher_source_authentication()
        different_source["archive_sha256"] = "0" * 64

        with self._verification_context(), mock.patch.object(
            publisher,
            "_runtime_publisher_source_authentication",
            return_value=different_source,
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "differs from the runtime publisher source",
        ):
            publisher.reconcile_pressure_selection_publication(
                self.receipt_path,
                controller_state_attestation=controller_binding,
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                now=self.now,
            )

        self.assertFalse(self.receipt_path.exists())
        self.assertTrue(self._guard_path().is_file())
        self.assertFalse(self._seal_path().exists())

    def test_partial_controller_attestation_write_failure_is_removed(self) -> None:
        original = publisher._write_controller_member_at

        def fail_second_member(
            parent_descriptor: int, name: str, payload: bytes
        ) -> None:
            if name == "active_promotion.json":
                raise OSError("injected controller-state member failure")
            original(parent_descriptor, name, payload)

        with self._verification_context(), mock.patch.object(
            publisher,
            "_write_controller_member_at",
            side_effect=fail_second_member,
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "publication failed closed",
        ):
            self._publish()

        self.assertEqual(list(self.archive_root.iterdir()), [])
        self.assertFalse(self.receipt_path.exists())
        self.assertFalse(self._guard_path().exists())
        self.assertFalse(self._seal_path().exists())

    def test_reconciliation_rejects_a_different_controller_state_binding(self) -> None:
        with self._verification_context():
            self._publish()
            wrong_binding = self._controller_binding()
            wrong_binding["sha256"] = "0" * 64
            with self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "binds another controller state",
            ):
                publisher.reconcile_pressure_selection_publication(
                    self.receipt_path,
                    controller_state_attestation=wrong_binding,
                    expected_git_commit="e" * 40,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    now=self.now,
                )

    def test_published_receipt_rejects_guard_seal_and_controller_snapshot_tamper(
        self,
    ) -> None:
        with self._verification_context():
            self._publish()
            self._write_readonly(self._guard_path(), b"receipt publication is not authoritative\n")
            with self.assertRaises(ValueError):
                publisher.verify_published_pressure_selection_receipt(
                    self.receipt_path, authorized_pic_root=self.pic_root
                )
            self._guard_path().chmod(0o600)
            self._guard_path().unlink()

            self._seal_path().chmod(0o600)
            seal = json.loads(self._seal_path().read_text(encoding="utf-8"))
            seal["receipt_sha256"] = "0" * 64
            self._seal_path().write_bytes(_canonical(seal))
            self._seal_path().chmod(0o444)
            with self.assertRaises(ValueError):
                publisher.verify_published_pressure_selection_receipt(
                    self.receipt_path, authorized_pic_root=self.pic_root
                )

        controller_path = Path(self._controller_binding()["path"])
        controller_path.parent.chmod(0o700)
        controller_path.chmod(0o600)
        value = json.loads(controller_path.read_text(encoding="utf-8"))
        value["controller_state"]["registered_science_slices"] = [{}]
        controller_path.write_bytes(_canonical(value))
        controller_path.chmod(0o400)
        controller_path.parent.chmod(0o500)
        with self._verification_context(), self.assertRaises(ValueError):
            publisher.consume_sealed_controller_state_attestation(
                {"path": str(controller_path), "sha256": _sha256(controller_path.read_bytes())},
                authorized_pic_root=self.pic_root,
            )

    def test_module_has_no_scheduler_or_simulation_execution_path(self) -> None:
        source = inspect.getsource(publisher)
        self.assertNotIn("import subprocess", source)
        self.assertNotIn("TRUSTED_SBATCH", source)
        self.assertNotIn("os.system(", source)
        self.assertNotIn("pilot_publisher.subprocess.check_output", source)

    def test_executes_reanalysis_from_candidate_archive_modules(self) -> None:
        paths = (
            publisher.pressure_selection.pressure_review_packet_verifier.PRESSURE_REANALYSIS_SOURCE_PATHS
        )
        sources = {path: b"VALUE = 'unused'\n" for path in paths}
        sources["tst/publication/analyze_q011_section54_pressure_pilot.py"] = (
            b"def candidate_result():\n"
            b"    return {'status': 'archive-executed'}\n"
        )
        sources[
            "tst/publication/q011_section54_historical_pressure_pilot_consumer.py"
        ] = (
            b"from . import analyze_q011_section54_pressure_pilot as analyzer\n"
            b"def consume_exact_historical_production_pressure_pilot():\n"
            b"    return analyzer.candidate_result()\n"
        )
        self.assertEqual(
            publisher._execute_candidate_reanalysis(sources),
            {"status": "archive-executed"},
        )
        missing = dict(sources)
        del missing[paths[-1]]
        with self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "reanalysis archive omits",
        ):
            publisher._execute_candidate_reanalysis(missing)

    def test_candidate_reanalysis_absolute_imports_use_only_archive_modules(self) -> None:
        paths = (
            publisher.pressure_selection.pressure_review_packet_verifier.PRESSURE_REANALYSIS_SOURCE_PATHS
        )
        sources = {path: b"VALUE = 'unused'\n" for path in paths}
        sources["tst/publication/analyze_q011_section54_pressure_pilot.py"] = (
            b"def candidate_result():\n"
            b"    return {'status': 'archive-absolute-import'}\n"
        )
        sources[
            "tst/publication/q011_section54_historical_pressure_pilot_consumer.py"
        ] = (
            b"import tst.publication.analyze_q011_section54_pressure_pilot as analyzer\n"
            b"def consume_exact_historical_production_pressure_pilot():\n"
            b"    return analyzer.candidate_result()\n"
        )
        self.assertEqual(
            publisher._execute_candidate_reanalysis(sources),
            {"status": "archive-absolute-import"},
        )
        sources[
            "tst/publication/q011_section54_historical_pressure_pilot_consumer.py"
        ] = (
            b"import tst.publication.publish_q011_section54_pressure_selection\n"
            b"def consume_exact_historical_production_pressure_pilot():\n"
            b"    return {'status': 'escaped'}\n"
        )
        with self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "historical pressure-pilot recomputation failed",
        ):
            publisher._execute_candidate_reanalysis(sources)

    def test_candidate_archive_resource_limits_fail_closed(self) -> None:
        archive_stream = io.BytesIO()
        with tarfile.open(fileobj=archive_stream, mode="w") as archive:
            member = tarfile.TarInfo("one.py")
            member.size = 1
            archive.addfile(member, io.BytesIO(b"x"))
        with mock.patch.object(
            publisher, "MAX_CANDIDATE_ARCHIVE_REGULAR_FILES", 0
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "exceeds resource limits",
        ):
            publisher._validated_candidate_archive_regular_files(
                archive_stream.getvalue()
            )
        with mock.patch.object(
            publisher, "MAX_CANDIDATE_ARCHIVE_MEMBER_NAME_CHARACTERS", 5
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "exceeds resource limits",
        ):
            publisher._validated_candidate_archive_regular_files(
                archive_stream.getvalue()
            )

    def test_candidate_archive_raw_size_and_total_member_limits_fail_closed(self) -> None:
        archive_stream = io.BytesIO()
        with tarfile.open(fileobj=archive_stream, mode="w") as archive:
            archive.addfile(tarfile.TarInfo("directory/"))
            member = tarfile.TarInfo("directory/one.py")
            member.size = 1
            archive.addfile(member, io.BytesIO(b"x"))
        archive_payload = archive_stream.getvalue()
        archive_path = self.pic_root / "bounded-source.tar"
        self._write_readonly(archive_path, archive_payload)

        with self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "authorized clean-candidate source archive is unavailable",
        ):
            publisher._live_file(
                archive_path,
                self.pic_root,
                "authorized clean-candidate source archive",
                max_bytes=len(archive_payload) - 1,
            )

        with mock.patch.object(
            publisher, "MAX_CANDIDATE_ARCHIVE_RAW_BYTES", len(archive_payload) - 1
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "exceeds resource limits",
        ):
            publisher._validated_candidate_archive_regular_files(archive_payload)

        with mock.patch.object(
            publisher, "MAX_CANDIDATE_ARCHIVE_MEMBERS", 1
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "exceeds resource limits",
        ):
            publisher._validated_candidate_archive_regular_files(archive_payload)

    def test_candidate_archive_rejects_compression_and_wraps_malformed_tar(self) -> None:
        compressed_stream = io.BytesIO()
        with tarfile.open(fileobj=compressed_stream, mode="w:gz") as archive:
            member = tarfile.TarInfo("one.py")
            member.size = 1
            archive.addfile(member, io.BytesIO(b"x"))
        for payload in (compressed_stream.getvalue(), b"not a tar archive"):
            with self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "source archive is unreadable",
            ):
                publisher._validated_candidate_archive_regular_files(payload)
        with mock.patch.object(
            publisher.tarfile, "open", side_effect=RecursionError
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "source archive is unreadable",
        ):
            publisher._validated_candidate_archive_regular_files(b"plain tar probe")

    def test_live_capture_supplies_bounded_clean_candidate_reader(self) -> None:
        def revalidate(
            *_args: object,
            read_candidate_tree: object,
            **_kwargs: object,
        ) -> dict[str, object]:
            self.assertIs(
                read_candidate_tree, publisher._read_bounded_clean_candidate_tree
            )
            return {
                "clean_candidate_manifest": {
                    "expected_sha256": self.policy["science_submission_freeze"][
                        "manifest_sha256"
                    ],
                    "path": str(self.manifest_path),
                    "sha256": self.policy["science_submission_freeze"][
                        "manifest_sha256"
                    ],
                },
                "current_control_plane_version": self.version,
                "source": {"git_commit": "e" * 40},
                "status": "passed",
            }

        with self._verification_context(), mock.patch.object(
            publisher.clean_candidate_revalidator,
            "revalidate_clean_candidate",
            side_effect=revalidate,
        ):
            publisher._capture_live_controller_state(
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_bounded_clean_candidate_closure_rejects_archive_before_retention(
        self,
    ) -> None:
        directory = self.root / "bounded-candidate-probe"
        directory.mkdir()
        archive_stream = io.BytesIO()
        with tarfile.open(fileobj=archive_stream, mode="w") as archive:
            archive.addfile(tarfile.TarInfo("one/"))
            archive.addfile(tarfile.TarInfo("two/"))
        self._write_readonly(
            directory / "source.tar", _git_style_archive(archive_stream.getvalue())
        )

        closure = publisher._BoundedCleanCandidateClosure()
        descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
        try:
            with mock.patch.object(
                publisher, "MAX_CANDIDATE_ARCHIVE_MEMBERS", 1
            ), self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "exceeds resource limits",
            ):
                closure.read_regular_file_at(
                    descriptor,
                    "source.tar",
                    label="Clean-candidate source archive",
                )
        finally:
            closure.close()
            os.close(descriptor)

    def test_bounded_clean_candidate_closure_rejects_fifo_without_blocking(self) -> None:
        directory = self.root / "bounded-candidate-fifo-probe"
        directory.mkdir()
        os.mkfifo(directory / "source.tar", 0o400)

        closure = publisher._BoundedCleanCandidateClosure()
        descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
        try:
            with self.assertRaisesRegex(ValueError, "is not a regular file"):
                closure.read_regular_file_at(
                    descriptor,
                    "source.tar",
                    label="Clean-candidate source archive",
                )
        finally:
            closure.close()
            os.close(descriptor)

    def test_bounded_clean_candidate_closure_rejects_ambiguous_plain_tar_view(
        self,
    ) -> None:
        directory = self.root / "bounded-candidate-ambiguous-tar-probe"
        directory.mkdir()
        archive_stream = io.BytesIO()
        with tarfile.open(fileobj=archive_stream, mode="w") as archive:
            archive.addfile(tarfile.TarInfo("BZh9-polyglot-prefix"))
        self._write_readonly(directory / "source.tar", archive_stream.getvalue())

        closure = publisher._BoundedCleanCandidateClosure()
        descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
        try:
            with self.assertRaisesRegex(
                publisher.PressureSelectionPublicationError,
                "canonical uncompressed Git tar format",
            ):
                closure.read_regular_file_at(
                    descriptor,
                    "source.tar",
                    label="Clean-candidate source archive",
                )
        finally:
            closure.close()
            os.close(descriptor)

    def test_bounded_clean_candidate_closure_rejects_aggregate_archive_expansion(
        self,
    ) -> None:
        directory = self.root / "bounded-candidate-aggregate-probe"
        directory.mkdir()
        archive_stream = io.BytesIO()
        with tarfile.open(fileobj=archive_stream, mode="w") as archive:
            member = tarfile.TarInfo("one.py")
            member.size = 1
            archive.addfile(member, io.BytesIO(b"x"))
        for name in ("source.tar", "submodule.tar"):
            self._write_readonly(
                directory / name, _git_style_archive(archive_stream.getvalue())
            )

        closure = publisher._BoundedCleanCandidateClosure()
        descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
        try:
            with mock.patch.object(
                publisher, "MAX_CANDIDATE_TREE_ARCHIVE_REGULAR_BYTES", 1
            ):
                closure.read_regular_file_at(
                    descriptor,
                    "source.tar",
                    label="Clean-candidate source archive",
                )
                with self.assertRaisesRegex(
                    publisher.PressureSelectionPublicationError,
                    "tree exceeds resource limits",
                ):
                    closure.read_regular_file_at(
                        descriptor,
                        "submodule.tar",
                        label="Clean-candidate submodule archive",
                    )
        finally:
            closure.close()
            os.close(descriptor)

    def test_candidate_publication_rejects_stage4_source_discontinuity(self) -> None:
        authorization = self._candidate_publication_authorization()
        authorization["publisher_source_authentication"]["git_commit"] = "d" * 40
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "differs from the reviewed receipt or source",
        ):
            publisher.publish_pressure_selection(
                copy.deepcopy(self.receipt),
                candidate_publication_authorization=authorization,
                controller_operator_id="codex",
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                now=self.now,
            )

    def test_candidate_publication_reopens_and_rejects_equal_time_human_decision(
        self,
    ) -> None:
        decision_path = Path(self.human_decision_binding["path"])
        decision_path.chmod(0o600)
        decision = json.loads(decision_path.read_text(encoding="utf-8"))
        decision["reviewed_utc"] = "2026-06-05T12:01:00Z"
        decision_path.write_bytes(_canonical(decision))
        decision_path.chmod(0o400)
        authorization = self._candidate_publication_authorization()
        authorization["human_decision"]["sha256"] = _sha256(decision_path.read_bytes())
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "strictly later than sealed reanalysis",
        ):
            publisher.publish_pressure_selection(
                copy.deepcopy(self.receipt),
                candidate_publication_authorization=authorization,
                controller_operator_id="codex",
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                now=self.now,
            )

    def test_publication_and_consumption_reject_reviewer_older_than_human_decision(
        self,
    ) -> None:
        self.reviewer_reviewed_utc = "2026-06-05T12:01:00Z"
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "reviewer attestation differs from reopened human pressure-selection decision",
        ):
            self._publish()
        self.assertFalse(self.receipt_path.exists())

        self.reviewer_reviewed_utc = "2026-06-05T12:02:00Z"
        with self._verification_context():
            self._publish()
        self.reviewer_reviewed_utc = "2026-06-05T12:01:00Z"
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "reviewer attestation differs from reopened human pressure-selection decision",
        ):
            publisher.consume_published_pressure_selection_receipt(
                self.receipt_path, authorized_pic_root=self.pic_root
            )

    def test_rejects_forged_self_consistent_reanalysis_source_closure(self) -> None:
        self.reanalysis_source_closure[0]["sha256"] = "0" * 64
        receipt = copy.deepcopy(self.receipt)
        with self._verification_context(), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "source closure differs from the active clean-candidate archive",
        ):
            publisher.publish_pressure_selection(
                receipt,
                candidate_publication_authorization=self._candidate_publication_authorization(
                    receipt
                ),
                controller_operator_id="codex",
                expected_git_commit="e" * 40,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                now=self.now,
            )

    def test_locked_mirrored_ledger_capture_rejects_outstanding_reservation(self) -> None:
        @contextmanager
        def snapshot(*_args: object, **_kwargs: object) -> Iterator[list[dict[str, object]]]:
            yield [{"event_sha256": "8" * 64}]

        with mock.patch.object(
            publisher.ledger,
            "require_no_incomplete_manual_accounting_marker",
            return_value=None,
        ), mock.patch.object(
            publisher.ledger,
            "validated_read_only_mirrored_state_snapshot",
            side_effect=snapshot,
        ), mock.patch.object(
            publisher.ledger,
            "latest_reservations",
            return_value={"active": {"state": "submitted"}},
        ), mock.patch.object(
            publisher.ledger,
            "accounting",
            return_value={
                "currently_reserved_node_hours": 1.0,
                "cumulative_consumed_node_hours": 2.0,
            },
        ), self.assertRaisesRegex(
            publisher.PressureSelectionPublicationError,
            "outstanding reservation blocks pressure selection",
        ):
            publisher._capture_mirrored_ledger_state(
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_cli_failure_emits_canonical_structured_recovery_json(self) -> None:
        recovery = {
            "schema_version": 1,
            "record_type": publisher.RECOVERY_GUARD_RECORD_TYPE,
        }
        stdout = io.StringIO()
        stderr = io.StringIO()
        with mock.patch.object(
            publisher,
            "verify_published_pressure_selection_live_state",
            side_effect=publisher.PressureSelectionPublicationError(
                "injected failure", recovery=recovery
            ),
        ), redirect_stdout(stdout), redirect_stderr(stderr):
            status = publisher.main(
                [
                    "verify-published",
                    "--receipt",
                    str(self.receipt_path),
                    "--authorized-pic-root",
                    str(self.pic_root),
                    "--authorized-project-home-root",
                    str(self.project_home_root),
                ]
            )
        self.assertEqual(status, 2)
        self.assertEqual(stdout.getvalue(), "")
        failure = json.loads(stderr.getvalue())
        self.assertEqual(failure["status"], "failed_closed")
        self.assertEqual(failure["recovery"], recovery)
        self.assertEqual(stderr.getvalue().encode("utf-8"), _canonical(failure))


if __name__ == "__main__":
    unittest.main()
