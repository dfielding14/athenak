#!/usr/bin/env python3
"""Focused tests for source-local Q-011 Section 5.4 raw-attempt freezing."""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import stat
import tempfile
import threading
from typing import Any, Iterator
import unittest
from unittest import mock

from tst.publication import analyze_q011_section54_campaign as campaign
from tst.publication import immutable_orion_tree
from tst.publication import publish_q011_section54_campaign_attempt as publisher
from tst.publication import test_analyze_q011_section54_campaign as fixtures


_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _put(root: Path, relative: str, payload: bytes) -> dict[str, str]:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return {"path": relative, "sha256": _sha256(payload)}


def _write_manifest(source: Path, manifest: dict[str, Any]) -> None:
    (source / campaign.MANIFEST_NAME).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _tamper_read_only_member(path: Path) -> None:
    mode = stat.S_IMODE(path.stat().st_mode)
    path.chmod(mode | stat.S_IWUSR)
    path.write_bytes(path.read_bytes() + b"\npost-admission tamper\n")
    path.chmod(mode)


@contextmanager
def _concurrent_tamper_at_first_rename(
    destination_parent: Path, relative: str
) -> Iterator[mock.Mock]:
    original_rename = publisher._rename_no_replace_at
    requested = threading.Event()
    finished = threading.Event()
    member: list[Path] = []
    errors: list[BaseException] = []

    def tamper() -> None:
        try:
            if not requested.wait(timeout=5):
                raise AssertionError("publisher did not request concurrent retained-member tamper")
            if member:
                _tamper_read_only_member(member[0])
        except BaseException as error:
            errors.append(error)
        finally:
            finished.set()

    worker = threading.Thread(target=tamper)
    worker.start()
    tampered = False

    def tamper_then_rename(
        source_parent_fd: int,
        source_name: str,
        destination_parent_fd: int,
        destination_name: str,
    ) -> None:
        nonlocal tampered
        if not tampered:
            staging = Path("/proc/self/fd") / str(source_parent_fd) / source_name
            mode = stat.S_IMODE(staging.stat().st_mode)
            staging.chmod(mode | stat.S_IXUSR)
            member.append(staging / relative)
            requested.set()
            if not finished.wait(timeout=5):
                raise AssertionError("concurrent retained-member tamper did not finish")
            staging.chmod(mode)
            tampered = True
        original_rename(
            source_parent_fd, source_name, destination_parent_fd, destination_name
        )

    try:
        with mock.patch.object(
            publisher, "_rename_no_replace_at", side_effect=tamper_then_rename
        ) as rename:
            yield rename
    finally:
        requested.set()
        worker.join(timeout=5)
    if worker.is_alive():
        raise AssertionError("concurrent retained-member tamper worker did not finish")
    if errors:
        raise AssertionError("concurrent retained-member tamper failed") from errors[0]


def _omit_product(
    source: Path, manifest: dict[str, Any], *, kind: str, snapshot_time: float
) -> None:
    product = fixtures._find_product(manifest, kind, snapshot_time)
    manifest["products"].remove(product)
    (source / product["path"]).unlink()
    _write_manifest(source, manifest)


def _failed_manifest(
    source: Path,
    *,
    attempt_id: str = "q011-failed-attempt-001",
    logs: tuple[tuple[str, str, bytes], ...] = (
        ("stderr", "logs/stderr.txt", b"fixture failure\n"),
        ("stdout", "logs/stdout.txt", b"fixture startup\n"),
    ),
    retained: tuple[tuple[str, bytes], ...] = (),
) -> dict[str, Any]:
    failure_logs = []
    for kind, relative, payload in logs:
        failure_logs.append({"kind": kind, **_put(source, relative, payload)})
    retained_artifacts = [
        _put(source, relative, payload) for relative, payload in retained
    ]
    return {
        "schema_version": 1,
        "record_type": publisher.FAILED_MANIFEST_RECORD_TYPE,
        "attempt_status": "failed",
        "attempt_id": attempt_id,
        "failure_logs": failure_logs,
        "retained_artifacts": retained_artifacts,
    }


@contextmanager
def _raw_attempt() -> Iterator[tuple[Path, Path, Path, dict[str, Any]]]:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "raw"
        orion = root / "orion"
        source.mkdir()
        orion.mkdir()
        with mock.patch.object(campaign, "ORION_BULK_ROOT", orion):
            manifest = fixtures._manifest(source)
            contract_binding = manifest["artifact_bindings"]["attempt_contract"]
            contract = json.loads((source / contract_binding["path"]).read_text(encoding="utf-8"))
            destination_parent = Path(contract["authorized_orion_attempt_root"]).parent
            destination_parent.mkdir(parents=True)
            _write_manifest(source, manifest)
            try:
                yield root, source, destination_parent, manifest
            finally:
                for path in destination_parent.iterdir():
                    fixtures._make_writable_tree(path)


@contextmanager
def _failed_raw_attempt(
    *,
    attempt_id: str = "q011-failed-attempt-001",
    logs: tuple[tuple[str, str, bytes], ...] = (
        ("stderr", "logs/stderr.txt", b"fixture failure\n"),
        ("stdout", "logs/stdout.txt", b"fixture startup\n"),
    ),
    retained: tuple[tuple[str, bytes], ...] = (),
) -> Iterator[tuple[Path, Path, Path, dict[str, Any]]]:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "raw"
        destination_parent = root / "retained"
        source.mkdir()
        destination_parent.mkdir()
        manifest = _failed_manifest(
            source, attempt_id=attempt_id, logs=logs, retained=retained
        )
        _write_manifest(source, manifest)
        try:
            yield root, source, destination_parent, manifest
        finally:
            for path in destination_parent.iterdir():
                fixtures._make_writable_tree(path)


def _completed_external_closure() -> dict[str, object]:
    return {"validation": "fixture_external_candidate_closure"}


def _completed_pressure_pilot_publication(
    path: str | Path, *, authorized_pic_root: Path
) -> dict[str, str]:
    del authorized_pic_root
    payload = Path(path).read_bytes()
    parsed = json.loads(payload)
    return {
        "receipt_sha256": _sha256(payload),
        "manifest_sha256": parsed["aggregate_bundle"]["manifest_sha256"],
        "analysis_result_sha256": parsed["aggregate_analysis"]["sha256"],
    }


def _completed_registered_execution_ledger(
    receipt: dict[str, Any], **_kwargs: object
) -> dict[str, str]:
    return {
        "reservation_id": receipt["reservation_id"],
        "submission_id": receipt["submission_id"],
        "reconciliation_event_sha256": receipt["reconciliation_event_sha256"],
    }


@contextmanager
def _completed_registered_execution_ledger_snapshot(
    receipt: dict[str, Any], **kwargs: object
) -> Iterator[dict[str, str]]:
    yield _completed_registered_execution_ledger(receipt, **kwargs)


@contextmanager
def _completed_planner_tree(
    runtime_root: str | Path, expected_inventory_sha256: str, **kwargs: Any
) -> Iterator[tuple[dict[str, Any], immutable_orion_tree.VerifiedFrozenTree]]:
    kwargs["authorized_root"] = Path(runtime_root).parent
    with immutable_orion_tree.staged_verified_frozen_tree(
        runtime_root, expected_inventory_sha256, **kwargs
    ) as verified:
        yield verified


@contextmanager
def _completed_external_closures() -> Iterator[None]:
    with mock.patch.object(
        campaign,
        "_validate_external_clean_candidate_closure",
        return_value=_completed_external_closure(),
    ), mock.patch.object(
        campaign.pressure_selection.pressure_pilot_publisher,
        "verify_published_pressure_pilot_receipt",
        side_effect=_completed_pressure_pilot_publication,
    ), mock.patch.object(
        campaign,
        "staged_verified_frozen_tree",
        side_effect=_completed_planner_tree,
    ), mock.patch.object(
        campaign,
        "_validated_registered_execution_receipt_ledger_snapshot",
        side_effect=_completed_registered_execution_ledger_snapshot,
    ), mock.patch.object(
        campaign,
        "validate_planner_retention_binding",
        side_effect=lambda value, **_kwargs: value,
    ):
        yield


def _freeze_completed(
    source: Path, destination_parent: Path, **kwargs: Any
) -> dict[str, Any]:
    with _completed_external_closures():
        return publisher.freeze_campaign_attempt(
            source, destination_parent, attempt_status="completed", **kwargs
        )


class Q011Section54CampaignAttemptPublisherTests(unittest.TestCase):
    def test_completed_attempt_is_exclusively_frozen_with_analyzer_metadata(self) -> None:
        with _raw_attempt() as (root, source, destination_parent, manifest):
            with mock.patch.object(
                publisher, "_rename_no_replace_at", wraps=publisher._rename_no_replace_at
            ) as rename:
                result = _freeze_completed(source, destination_parent)
            rename.assert_called_once()
            destination = Path(result["campaign_root"])
            self.assertEqual(destination.name, manifest["run_identity"]["attempt_id"])
            self.assertTrue(result["prepublication_admission_passed"])
            self.assertTrue(result["admission_self_check_available"])
            self.assertIsNone(result["admission_self_check"])
            self.assertTrue((destination / campaign.MANIFEST_NAME).is_file())
            self.assertTrue((destination / immutable_orion_tree.INVENTORY_NAME).is_file())
            self.assertTrue((destination / immutable_orion_tree.FREEZE_RECEIPT_NAME).is_file())
            retained_manifest = json.loads(
                (destination / campaign.MANIFEST_NAME).read_text(encoding="utf-8")
            )
            self.assertEqual(
                retained_manifest["authorized_orion_campaign_root"], str(destination)
            )
            receipt = json.loads(
                (destination / immutable_orion_tree.FREEZE_RECEIPT_NAME).read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(receipt, campaign._EXPECTED_FREEZE_RECEIPT)
            self.assertFalse(
                any(
                    path.stat().st_mode & _WRITE_BITS
                    for path in [destination, *destination.rglob("*")]
                )
            )
            verified = immutable_orion_tree.verify_frozen_tree(
                destination,
                result["inventory_sha256"],
                authorized_root=destination_parent,
            )
            self.assertTrue(verified["recursively_read_only"])
            with _completed_external_closures():
                admitted = publisher.run_final_admission_self_check(
                    destination,
                    result["inventory_sha256"],
                    authorized_destination_root=destination_parent,
                    authorized_pic_root=root,
                )
            self.assertTrue(admitted["admitted_for_follow_on_numerical_qualification"])

    def test_completed_attempt_keeps_external_pic_root_distinct_from_destination_root(self) -> None:
        with _raw_attempt() as (root, source, destination_parent, _):
            with mock.patch.object(
                campaign,
                "_validated_retained_attempt_semantics_snapshot",
                wraps=campaign._validated_retained_attempt_semantics_snapshot,
            ) as retained_attempt_semantics:
                _freeze_completed(
                    source,
                    destination_parent,
                    authorized_pic_root=root,
                )
            self.assertEqual(
                retained_attempt_semantics.call_args.kwargs["authorized_pic_root"],
                root,
            )

    def test_failed_partial_attempt_is_frozen_without_admission_eligibility(self) -> None:
        with _failed_raw_attempt(
            retained=(("bin/q011.j2.partial.bin", b"partial emitted j2\n"),)
        ) as (_, source, destination_parent, _):
            result = publisher.freeze_campaign_attempt(
                source, destination_parent, attempt_status="failed"
            )
            destination = Path(result["campaign_root"])
            self.assertFalse(result["prepublication_admission_passed"])
            self.assertFalse(result["admission_self_check_available"])
            self.assertTrue((destination / "bin/q011.j2.partial.bin").is_file())
            receipt = json.loads(
                (destination / immutable_orion_tree.FREEZE_RECEIPT_NAME).read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(receipt, publisher._FAILED_RECEIPT)
            with self.assertRaisesRegex(
                publisher.PublicationError, "requires a completed attempt receipt"
            ):
                publisher.run_final_admission_self_check(
                    destination,
                    result["inventory_sha256"],
                    authorized_destination_root=destination_parent,
                )

    def test_failed_zero_scientific_product_attempt_retains_available_logs(self) -> None:
        with _failed_raw_attempt() as (_, source, destination_parent, manifest):
            self.assertEqual(manifest["retained_artifacts"], [])
            result = publisher.freeze_campaign_attempt(
                source, destination_parent, attempt_status="failed"
            )
            destination = Path(result["campaign_root"])
            self.assertEqual((destination / "logs/stderr.txt").read_bytes(), b"fixture failure\n")
            self.assertEqual((destination / "logs/stdout.txt").read_bytes(), b"fixture startup\n")

    def test_failed_attempt_requires_canonical_status_and_meaningful_failure_log(self) -> None:
        with _failed_raw_attempt(logs=()) as (_, source, destination_parent, _):
            with self.assertRaisesRegex(publisher.PublicationError, "at least one available log"):
                publisher.freeze_campaign_attempt(
                    source, destination_parent, attempt_status="failed"
                )
            self.assertEqual(list(destination_parent.iterdir()), [])

        with _failed_raw_attempt(logs=(("stderr", "logs/stderr.txt", b""),)) as (
            _,
            source,
            destination_parent,
            _,
        ):
            with self.assertRaisesRegex(publisher.PublicationError, "only empty failure logs"):
                publisher.freeze_campaign_attempt(
                    source, destination_parent, attempt_status="failed"
                )
            self.assertEqual(list(destination_parent.iterdir()), [])

        with _failed_raw_attempt() as (_, source, destination_parent, manifest):
            manifest["attempt_status"] = "completed"
            _write_manifest(source, manifest)
            with self.assertRaisesRegex(publisher.PublicationError, "must be 'failed'"):
                publisher.freeze_campaign_attempt(
                    source, destination_parent, attempt_status="failed"
                )
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_completed_attempt_missing_required_retained_product_is_rejected(self) -> None:
        with _raw_attempt() as (_, source, destination_parent, manifest):
            _omit_product(source, manifest, kind="j2", snapshot_time=100.0)
            with self.assertRaises(campaign.QualificationError):
                publisher.freeze_campaign_attempt(
                    source, destination_parent, attempt_status="completed"
                )
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_existing_destination_is_not_overwritten(self) -> None:
        with _raw_attempt() as (_, source, destination_parent, _):
            first = _freeze_completed(source, destination_parent)
            inventory = (
                Path(first["campaign_root"]) / immutable_orion_tree.INVENTORY_NAME
            ).read_bytes()
            with self.assertRaisesRegex(
                publisher.PublicationError, "destination already exists"
            ):
                _freeze_completed(source, destination_parent)
            self.assertEqual(
                (Path(first["campaign_root"]) / immutable_orion_tree.INVENTORY_NAME).read_bytes(),
                inventory,
            )

    def test_destination_appearing_at_atomic_rename_is_not_overwritten(self) -> None:
        with _raw_attempt() as (_, source, destination_parent, manifest):
            destination = destination_parent / manifest["run_identity"]["attempt_id"]
            original_rename = publisher._rename_no_replace_at

            def collide(
                source_parent_fd: int,
                source_name: str,
                destination_parent_fd: int,
                destination_name: str,
            ) -> None:
                destination.mkdir()
                (destination / "sentinel.txt").write_text("competing publication\n", encoding="utf-8")
                original_rename(
                    source_parent_fd, source_name, destination_parent_fd, destination_name
                )

            with mock.patch.object(publisher, "_rename_no_replace_at", side_effect=collide):
                with self.assertRaisesRegex(
                    publisher.PublicationError, "destination already exists"
                ):
                    _freeze_completed(source, destination_parent)
            self.assertEqual(
                (destination / "sentinel.txt").read_text(encoding="utf-8"),
                "competing publication\n",
            )
            self.assertEqual(
                [path.name for path in destination_parent.iterdir()], [destination.name]
            )

    def test_staging_path_swap_after_admission_is_not_published_or_removed(self) -> None:
        with _raw_attempt() as (_, source, destination_parent, manifest):
            destination = destination_parent / manifest["run_identity"]["attempt_id"]
            replacement: list[Path] = []

            def swap(staging: Path, *_: Any, **__: Any) -> dict[str, Any]:
                staging.rename(staging.with_name(f"{staging.name}.moved"))
                staging.mkdir()
                (staging / "sentinel.txt").write_text("replacement staging\n", encoding="utf-8")
                replacement.append(staging)
                return {}

            with mock.patch.object(
                publisher, "_validate_staged_completed_admission", side_effect=swap
            ):
                with self.assertRaisesRegex(publisher.PublicationError, "staging tree changed"):
                    _freeze_completed(source, destination_parent)
            self.assertFalse(destination.exists())
            self.assertEqual(
                (replacement[0] / "sentinel.txt").read_text(encoding="utf-8"),
                "replacement staging\n",
            )

    def test_private_container_substitution_never_exposes_replacement(self) -> None:
        with _failed_raw_attempt() as (_, source, destination_parent, manifest):
            destination = destination_parent / manifest["attempt_id"]
            original_rename = publisher._rename_no_replace_at
            substituted = False
            parked: list[Path] = []
            observed_public_replacement: list[bool] = []
            safe_public_inodes: list[int] = []
            observed_public_inodes: list[int] = []
            observed_public_modes: list[int] = []

            def substitute_then_rename(
                source_parent_fd: int,
                source_name: str,
                destination_parent_fd: int,
                destination_name: str,
            ) -> None:
                nonlocal substituted
                if not substituted:
                    private = next(
                        path
                        for path in destination_parent.iterdir()
                        if path.name.startswith(f".{destination.name}.staging-")
                    )
                    parked_root = private.with_name(f"{private.name}.parked")
                    private.rename(parked_root)
                    safe_public_inodes.append(
                        (parked_root / publisher._STAGING_ROOT_NAME).stat().st_ino
                    )
                    private.mkdir()
                    (private / publisher._STAGING_ROOT_NAME).mkdir()
                    (private / publisher._STAGING_ROOT_NAME / "sentinel.txt").write_text(
                        "attacker replacement\n", encoding="utf-8"
                    )
                    parked.append(parked_root)
                    substituted = True
                original_rename(
                    source_parent_fd, source_name, destination_parent_fd, destination_name
                )
                published = destination.stat()
                observed_public_inodes.append(published.st_ino)
                observed_public_modes.append(stat.S_IMODE(published.st_mode))
                try:
                    replacement_is_public = (destination / "sentinel.txt").exists()
                except PermissionError:
                    replacement_is_public = False
                observed_public_replacement.append(replacement_is_public)

            with mock.patch.object(
                publisher, "_rename_no_replace_at", side_effect=substitute_then_rename
            ):
                result = publisher.freeze_campaign_attempt(
                    source, destination_parent, attempt_status="failed"
                )
            self.assertEqual(Path(result["campaign_root"]), destination)
            self.assertTrue(destination.is_dir())
            self.assertEqual(observed_public_replacement, [False])
            self.assertEqual(observed_public_inodes, safe_public_inodes)
            self.assertEqual(observed_public_modes, [stat.S_IWUSR])
            self.assertTrue(parked)
            self.assertFalse((destination / "sentinel.txt").exists())

    def test_empty_directory_injected_at_atomic_publish_is_rolled_back(self) -> None:
        with _failed_raw_attempt() as (_, source, destination_parent, manifest):
            destination = destination_parent / manifest["attempt_id"]
            original_rename = publisher._rename_no_replace_at
            injected = False

            def inject_then_rename(
                source_parent_fd: int,
                source_name: str,
                destination_parent_fd: int,
                destination_name: str,
            ) -> None:
                nonlocal injected
                if not injected:
                    staging = Path("/proc/self/fd") / str(source_parent_fd) / source_name
                    mode = stat.S_IMODE(staging.stat().st_mode)
                    staging.chmod(mode | stat.S_IWUSR | stat.S_IXUSR)
                    (staging / "injected-empty").mkdir()
                    (staging / "injected-empty").chmod(0o500)
                    staging.chmod(mode)
                    injected = True
                original_rename(
                    source_parent_fd, source_name, destination_parent_fd, destination_name
                )

            with mock.patch.object(
                publisher, "_rename_no_replace_at", side_effect=inject_then_rename
            ) as rename:
                with self.assertRaisesRegex(
                    publisher.PublicationError, "directory membership drifted"
                ):
                    publisher.freeze_campaign_attempt(
                        source, destination_parent, attempt_status="failed"
                    )
            self.assertEqual(rename.call_count, 2)
            self.assertFalse(destination.exists())
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_cleanup_staging_does_not_delete_path_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            staging = parent / ".staging"
            staging.mkdir()
            (staging / "owned.txt").write_text("owned\n", encoding="utf-8")
            parent_fd = os.open(parent, publisher._DIRECTORY_FLAGS)
            staging_fd = os.open(staging.name, publisher._DIRECTORY_FLAGS, dir_fd=parent_fd)
            try:
                parked = parent / ".parked"
                staging.rename(parked)
                staging.mkdir()
                (staging / "sentinel.txt").write_text("replacement\n", encoding="utf-8")
                publisher._cleanup_staging(parent_fd, staging.name, staging_fd)
                self.assertEqual(
                    (staging / "sentinel.txt").read_text(encoding="utf-8"),
                    "replacement\n",
                )
                self.assertTrue((parked / "owned.txt").is_file())
            finally:
                os.close(staging_fd)
                os.close(parent_fd)

    def test_rollback_quarantines_moved_original_and_preserves_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            destination = parent / "retained-attempt"
            destination.mkdir()
            (destination / "owned.txt").write_text("owned\n", encoding="utf-8")
            moved = parent / "moved-public-attempt"
            parent_fd = os.open(parent, publisher._DIRECTORY_FLAGS)
            destination_fd = os.open(
                destination.name, publisher._DIRECTORY_FLAGS, dir_fd=parent_fd
            )
            original_rename = publisher._rename_no_replace_at
            raced = False

            def race_then_rename(*args: object, **kwargs: object) -> None:
                nonlocal raced
                if not raced:
                    destination.rename(moved)
                    destination.mkdir()
                    (destination / "sentinel.txt").write_text(
                        "replacement\n", encoding="utf-8"
                    )
                    raced = True
                original_rename(*args, **kwargs)

            try:
                with mock.patch.object(
                    publisher, "_rename_no_replace_at", side_effect=race_then_rename
                ), self.assertRaisesRegex(
                    publisher.PublicationError, "substituted public destination"
                ):
                    publisher._rollback_published_destination(
                        parent_fd, destination.name, destination_fd
                    )
                self.assertEqual(
                    (destination / "sentinel.txt").read_text(encoding="utf-8"),
                    "replacement\n",
                )
                self.assertFalse(moved.exists())
                self.assertFalse(
                    any(
                        path.name.startswith(f".{destination.name}.rollback-")
                        for path in parent.iterdir()
                    )
                )
            finally:
                os.close(destination_fd)
                os.close(parent_fd)

    def test_retained_member_tamper_after_admission_is_rejected_before_rename(self) -> None:
        with _raw_attempt() as (_, source, destination_parent, _):
            original_admission = publisher._validate_staged_completed_admission

            def admit_and_tamper(staging: Path, *args: Any, **kwargs: Any) -> dict[str, Any]:
                result = original_admission(staging, *args, **kwargs)
                _tamper_read_only_member(staging / campaign.MANIFEST_NAME)
                return result

            with mock.patch.object(
                publisher, "_validate_staged_completed_admission", side_effect=admit_and_tamper
            ), mock.patch.object(
                publisher, "_rename_no_replace_at", wraps=publisher._rename_no_replace_at
            ) as rename:
                with self.assertRaises(publisher.PublicationError):
                    _freeze_completed(source, destination_parent)
            rename.assert_not_called()
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_retained_member_tamper_at_rename_is_rolled_back(self) -> None:
        with _raw_attempt() as (_, source, destination_parent, manifest):
            destination = destination_parent / manifest["run_identity"]["attempt_id"]
            with _concurrent_tamper_at_first_rename(
                destination_parent, campaign.MANIFEST_NAME
            ) as rename:
                with self.assertRaises(publisher.PublicationError):
                    _freeze_completed(source, destination_parent)
            self.assertEqual(rename.call_count, 2)
            self.assertFalse(destination.exists())
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_failed_attempt_retained_member_tamper_at_rename_is_rolled_back(self) -> None:
        with _failed_raw_attempt() as (_, source, destination_parent, manifest):
            destination = destination_parent / manifest["attempt_id"]
            with mock.patch.object(
                publisher, "_validate_staged_completed_admission"
            ) as admission, _concurrent_tamper_at_first_rename(
                destination_parent, "logs/stderr.txt"
            ) as rename:
                with self.assertRaises(publisher.PublicationError):
                    publisher.freeze_campaign_attempt(
                        source, destination_parent, attempt_status="failed"
                    )
            admission.assert_not_called()
            self.assertEqual(rename.call_count, 2)
            self.assertFalse(destination.exists())
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_post_rename_verification_failure_with_cleanup_error_withdraws_destination(
        self,
    ) -> None:
        with _failed_raw_attempt() as (_, source, destination_parent, manifest):
            destination = destination_parent / manifest["attempt_id"]
            original_rename = publisher._rename_no_replace_at
            published = False

            def rename_then_tamper(
                source_parent_fd: int,
                source_name: str,
                destination_parent_fd: int,
                destination_name: str,
            ) -> None:
                nonlocal published
                original_rename(
                    source_parent_fd, source_name, destination_parent_fd, destination_name
                )
                if not published:
                    mode = stat.S_IMODE(destination.stat().st_mode)
                    destination.chmod(mode | stat.S_IXUSR)
                    _tamper_read_only_member(destination / "logs/stderr.txt")
                    destination.chmod(mode)
                    published = True

            with mock.patch.object(
                publisher, "_rename_no_replace_at", side_effect=rename_then_tamper
            ), mock.patch.object(
                publisher,
                "_remove_anchored_tree_at",
                side_effect=publisher.PublicationError("fixture cleanup failure"),
            ):
                with self.assertRaisesRegex(
                    publisher.PublicationError,
                    "cannot remove invalid retained attempt destination",
                ):
                    publisher.freeze_campaign_attempt(
                        source, destination_parent, attempt_status="failed"
                    )
            self.assertFalse(destination.exists())
            self.assertEqual(
                len(
                    [
                        path
                        for path in destination_parent.iterdir()
                        if path.name.startswith(f".{destination.name}.rollback-")
                    ]
                ),
                1,
            )

    def test_source_member_symlink_is_rejected(self) -> None:
        with _failed_raw_attempt() as (root, source, destination_parent, manifest):
            stdout = next(
                binding for binding in manifest["failure_logs"] if binding["kind"] == "stdout"
            )
            outside = root / "outside.txt"
            outside.write_text("outside\n", encoding="utf-8")
            (source / stdout["path"]).unlink()
            (source / stdout["path"]).symlink_to(outside)
            with self.assertRaisesRegex(publisher.PublicationError, "contains symlink"):
                publisher.freeze_campaign_attempt(
                    source, destination_parent, attempt_status="failed"
                )
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_source_member_hardlink_is_rejected(self) -> None:
        with _failed_raw_attempt() as (_, source, destination_parent, manifest):
            stderr = next(
                binding for binding in manifest["failure_logs"] if binding["kind"] == "stderr"
            )
            os.link(source / stderr["path"], source / "logs/stderr-alias.txt")
            with self.assertRaisesRegex(publisher.PublicationError, "hard-linked"):
                publisher.freeze_campaign_attempt(
                    source, destination_parent, attempt_status="failed"
                )
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_undeclared_source_member_is_rejected(self) -> None:
        with _failed_raw_attempt() as (_, source, destination_parent, _):
            (source / "undeclared.log").write_text("undeclared\n", encoding="utf-8")
            with self.assertRaisesRegex(
                publisher.PublicationError, "membership differs"
            ):
                publisher.freeze_campaign_attempt(
                    source, destination_parent, attempt_status="failed"
                )
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_attempt_id_path_escape_is_rejected_before_destination_creation(self) -> None:
        with _failed_raw_attempt() as (_, source, destination_parent, manifest):
            manifest["attempt_id"] = "../escape"
            _write_manifest(source, manifest)
            with self.assertRaises(publisher.PublicationError):
                publisher.freeze_campaign_attempt(
                    source, destination_parent, attempt_status="failed"
                )
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_declared_source_path_escape_is_rejected_before_destination_creation(self) -> None:
        with _failed_raw_attempt() as (_, source, destination_parent, manifest):
            manifest["failure_logs"][0]["path"] = "../outside.txt"
            _write_manifest(source, manifest)
            with self.assertRaises(publisher.PublicationError):
                publisher.freeze_campaign_attempt(
                    source, destination_parent, attempt_status="failed"
                )
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_completed_attempt_can_invoke_exposed_final_admission_self_check(self) -> None:
        with _raw_attempt() as (_, source, destination_parent, _):
            admitted = {"admitted_for_follow_on_numerical_qualification": True}
            with mock.patch.object(
                publisher, "run_final_admission_self_check", return_value=admitted
            ) as self_check:
                result = _freeze_completed(
                    source,
                    destination_parent,
                    run_admission_self_check=True,
                )
            self.assertEqual(result["admission_self_check"], admitted)
            self_check.assert_called_once_with(
                Path(result["campaign_root"]),
                result["inventory_sha256"],
                authorized_destination_root=destination_parent,
                authorized_pic_root=campaign.ORION_BULK_ROOT,
            )

    def test_final_admission_self_check_failure_rolls_back_publication(self) -> None:
        with _raw_attempt() as (_, source, destination_parent, manifest):
            destination = destination_parent / manifest["run_identity"]["attempt_id"]
            with mock.patch.object(
                publisher,
                "run_final_admission_self_check",
                side_effect=publisher.PublicationError("fixture final self-check failure"),
            ):
                with self.assertRaisesRegex(
                    publisher.PublicationError, "fixture final self-check failure"
                ):
                    _freeze_completed(
                        source,
                        destination_parent,
                        run_admission_self_check=True,
                    )
            self.assertFalse(destination.exists())
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_completed_attempt_rejects_nonplanner_destination(self) -> None:
        with _raw_attempt() as (root, source, _, _):
            wrong_parent = root / "wrong-retained-parent"
            wrong_parent.mkdir()
            with self.assertRaisesRegex(
                publisher.PublicationError, "differs from planner-authorized retained root"
            ):
                _freeze_completed(source, wrong_parent)
            self.assertEqual(list(wrong_parent.iterdir()), [])

    def test_failed_attempt_cannot_request_final_admission_self_check(self) -> None:
        with _failed_raw_attempt() as (_, source, destination_parent, _):
            with self.assertRaisesRegex(
                publisher.PublicationError, "only available for completed attempts"
            ):
                publisher.freeze_campaign_attempt(
                    source,
                    destination_parent,
                    attempt_status="failed",
                    run_admission_self_check=True,
                )
            self.assertEqual(list(destination_parent.iterdir()), [])

    def test_invalid_completed_attempt_never_becomes_visible(self) -> None:
        with _raw_attempt() as (_, source, destination_parent, manifest):
            stdout = fixtures._find_product(manifest, "stdout", None)
            fixtures._rewrite_product(source, stdout, b"invalid retained stdout\n")
            _write_manifest(source, manifest)
            with mock.patch.object(
                campaign,
                "_validate_external_clean_candidate_closure",
                return_value=_completed_external_closure(),
            ):
                with self.assertRaisesRegex(
                    publisher.PublicationError, "failed prepublication campaign admission"
                ):
                    publisher.freeze_campaign_attempt(
                        source, destination_parent, attempt_status="completed"
                    )
            self.assertEqual(list(destination_parent.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
