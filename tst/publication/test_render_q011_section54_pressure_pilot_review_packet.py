#!/usr/bin/env python3
"""Focused adversarial tests for Q011 pressure-pilot review-packet publication."""

from __future__ import annotations

from contextlib import contextmanager
import fcntl
import inspect
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Iterator
import unittest
from unittest.mock import patch


PUBLICATION_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PUBLICATION_DIR / "frontier_control_plane"))

from tst.publication import publish_q011_section54_pressure_pilot_bundle as publisher  # noqa: E402
from tst.publication import render_q011_section54_pressure_pilot_review_packet as renderer  # noqa: E402
from tst.publication import test_publish_q011_section54_pressure_pilot_bundle as fixture  # noqa: E402


@contextmanager
def _published_aggregate() -> Iterator[tuple[Path, Path]]:
    with fixture._verified_raw_cases() as (base, roots, digests):
        receipt = base / "aggregate-receipt.json"
        publisher.publish_pressure_pilot_bundle(
            base / "aggregate-bundle",
            receipt_path=receipt,
            analysis_result_path=base / "aggregate-analysis.json",
            case_artifact_dirs=roots,
            case_descriptor_sha256=digests,
            authorized_pic_root=base.parent,
        )
        yield base, receipt


@contextmanager
def _fast_figures() -> Iterator[None]:
    def qualitative(
        _bundle_root: Path,
        _manifest: dict[str, object],
        output: Path,
        *,
        member_reader: object,
    ) -> None:
        del member_reader
        output.write_bytes(b"deterministic qualitative PNG fixture\n")

    def profiles(_analysis: dict[str, object], output: Path) -> None:
        output.write_bytes(b"deterministic profile PNG fixture\n")

    with patch.object(renderer, "_qualitative_figure", side_effect=qualitative), patch.object(
        renderer, "_profile_figure", side_effect=profiles
    ):
        yield


class Q011Section54PressurePilotReviewPacketTests(unittest.TestCase):
    def test_public_receipt_verifier_exposes_no_marker_bypass_flags(self) -> None:
        self.assertEqual(
            list(inspect.signature(renderer.verify_published_review_packet_receipt).parameters),
            ["receipt_path", "authorized_pic_root"],
        )

    def test_renderer_preflights_source_bindings_before_exposure(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            guard = base / publisher._publication_guard_name(receipt.name)
            original = publisher._source_bindings
            observed_calls = 0

            def inspect_preflight(*, require_verified_archive: bool) -> dict[str, object]:
                nonlocal observed_calls
                self.assertFalse(packet.exists())
                self.assertFalse(receipt.exists())
                self.assertFalse(guard.exists())
                observed_calls += 1
                return original(require_verified_archive=require_verified_archive)

            with patch.object(
                publisher,
                "_source_bindings",
                side_effect=inspect_preflight,
            ):
                renderer.render_packet(
                    aggregate_receipt,
                    packet,
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertEqual(observed_calls, 1)

    def test_terminal_products_accept_dt_scheduled_binary_headers(self) -> None:
        with _published_aggregate() as (base, _aggregate_receipt):
            bundle = base / "aggregate-bundle"
            manifest = publisher.pilot._manifest_schema(
                (bundle / publisher.pilot.MANIFEST_NAME).read_bytes()
            )
            products = renderer._terminal_products(
                bundle,
                manifest["cases"][0],
                member_reader=lambda relative: (bundle / relative).read_bytes(),
            )
            self.assertEqual(len(products), 5)
            self.assertEqual(products[2].shape, products[3].shape)
            self.assertEqual(products[2].shape, products[4].shape)

    def test_review_packet_publishes_and_verifies_exact_closure(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            renderer.render_packet(
                aggregate_receipt,
                packet,
                receipt,
                authorized_pic_root=base.parent,
            )
            verification = renderer.verify_published_review_packet_receipt(
                receipt, authorized_pic_root=base.parent
            )
            record = json.loads(receipt.read_text(encoding="utf-8"))
            self.assertEqual(
                verification["inventory_sha256"], record["inventory_sha256"]
            )
            self.assertEqual(
                {
                    path.relative_to(packet).as_posix()
                    for path in packet.rglob("*")
                    if path.is_file()
                },
                {*renderer.PACKET_MEMBERS, renderer.INVENTORY_NAME},
            )
            self.assertEqual(
                record["source_bindings"]["runtime_source_archive"],
                {
                    "execution_mode": "direct_api_nonproduction_only",
                    "git_commit": None,
                    "archive_sha256": None,
                    "verified_source_closure_sha256": None,
                },
            )
            self.assertEqual(
                record["consumption_rule"], publisher.CONSUMPTION_RULE
            )
            for path in [packet, *packet.rglob("*"), receipt]:
                self.assertFalse(path.stat().st_mode & 0o222)
            success_seal = (
                base.parent
                / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY
                / publisher._publication_seal_name(receipt.name)
            )
            self.assertTrue(success_seal.is_file())
            self.assertFalse(success_seal.stat().st_mode & 0o222)

    def test_late_packet_failure_retains_guarded_public_packet_and_receipt(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            with patch.object(
                renderer,
                "_verify_published_review_packet_receipt",
                side_effect=renderer.PacketError("injected late packet failure"),
            ), self.assertRaisesRegex(renderer.PacketError, "reviewed reconciliation required"):
                renderer.render_packet(
                    aggregate_receipt,
                    packet,
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(packet.is_dir())
            self.assertTrue(receipt.is_file())
            self.assertTrue(
                (base / publisher._publication_guard_name(receipt.name)).is_file()
            )
            self.assertFalse(
                any(".staging-" in path.name for path in base.iterdir())
            )

    def test_retained_packet_verifier_rejects_payload_tamper(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            renderer.render_packet(
                aggregate_receipt,
                packet,
                receipt,
                authorized_pic_root=base.parent,
            )
            fixture._make_writable(packet)
            markdown = packet / "PRESSURE_REVIEW_PACKET.md"
            markdown.write_bytes(markdown.read_bytes() + b"forged\n")
            fixture._freeze_existing(packet)
            with self.assertRaisesRegex(renderer.PacketError, "checksum drifted"):
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )
            fixture._make_writable(packet)

    def test_public_build_substitution_is_not_followed_or_deleted(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            original = publisher._require_same_directory_at
            observed = []

            def substitute(
                parent_descriptor: int, name: str, descriptor: int, label: str
            ) -> None:
                if label == "review-packet public build tree" and not observed:
                    observed.append(name)
                    os.rename(
                        name,
                        name + ".descriptor-anchor",
                        src_dir_fd=parent_descriptor,
                        dst_dir_fd=parent_descriptor,
                    )
                    os.mkdir(name, dir_fd=parent_descriptor)
                original(parent_descriptor, name, descriptor, label)

            with patch.object(
                publisher, "_require_same_directory_at", side_effect=substitute
            ), self.assertRaisesRegex(
                renderer.PacketError, "reviewed reconciliation required"
            ):
                renderer.render_packet(
                    aggregate_receipt,
                    packet,
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(packet.is_dir())
            self.assertFalse(receipt.exists())
            self.assertEqual(len(observed), 1)
            substituted = base / observed[0]
            anchored = base / (observed[0] + ".descriptor-anchor")
            self.assertTrue(substituted.is_dir())
            self.assertTrue(anchored.is_dir())
            self.assertTrue(
                (base / publisher._publication_guard_name(receipt.name)).is_file()
            )
            fixture._make_writable(anchored)

    def test_analysis_substitution_fails_receipt_bound_digest_check(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            original = renderer.ImmutableReadonlyPublicationFile.read

            def forged_analysis(
                member: renderer.ImmutableReadonlyPublicationFile,
            ) -> bytes:
                payload = original(member)
                if member._label != "aggregate analysis":
                    return payload
                analysis = json.loads(payload)
                analysis["status"] = "forged"
                return renderer._json_bytes(analysis)

            with patch.object(
                renderer.ImmutableReadonlyPublicationFile,
                "read",
                side_effect=forged_analysis,
                autospec=True,
            ), self.assertRaisesRegex(
                renderer.PacketError,
                "aggregate analysis SHA-256 differs from aggregate receipt binding",
            ):
                renderer.render_packet(
                    aggregate_receipt,
                    packet,
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertFalse(packet.exists())
            self.assertFalse(receipt.exists())

    def test_publication_root_substitution_fails_closed(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            anchor = base.parent / (base.name + ".descriptor-anchor")
            original = publisher._require_same_directory
            calls = 0

            def substitute(path: Path, descriptor: int, label: str) -> None:
                nonlocal calls
                if label == "review-packet publication root":
                    calls += 1
                    if calls == 2:
                        os.rename(path, anchor)
                        path.mkdir()
                original(path, descriptor, label)

            try:
                with patch.object(
                    publisher, "_require_same_directory", side_effect=substitute
                ), self.assertRaisesRegex(
                    renderer.PacketError,
                    "publication root changed",
                ):
                    renderer.render_packet(
                        aggregate_receipt,
                        packet,
                        receipt,
                        authorized_pic_root=base.parent,
                    )
                self.assertTrue((anchor / packet.name).is_dir())
                self.assertFalse((anchor / receipt.name).exists())
                self.assertTrue(
                    (
                        anchor / publisher._publication_guard_name(receipt.name)
                    ).is_file()
                )
            finally:
                if base.exists():
                    base.rmdir()
                if anchor.exists():
                    os.rename(anchor, base)

    def test_analysis_path_substitution_is_detected_while_descriptor_is_retained(
        self,
    ) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            analysis = base / "aggregate-analysis.json"
            anchor = base / "aggregate-analysis.json.descriptor-anchor"
            original = renderer.ImmutableReadonlyPublicationFile.read
            substituted = False

            def substitute_path(
                member: renderer.ImmutableReadonlyPublicationFile,
            ) -> bytes:
                nonlocal substituted
                if member._label == "aggregate analysis" and not substituted:
                    substituted = True
                    os.rename(
                        member._name,
                        anchor.name,
                        src_dir_fd=member._parent_descriptor,
                        dst_dir_fd=member._parent_descriptor,
                    )
                    descriptor = os.open(
                        member._name,
                        os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                        0o444,
                        dir_fd=member._parent_descriptor,
                    )
                    try:
                        os.write(descriptor, b"forged replacement\n")
                    finally:
                        os.close(descriptor)
                return original(member)

            try:
                with patch.object(
                    renderer.ImmutableReadonlyPublicationFile,
                    "read",
                    side_effect=substitute_path,
                    autospec=True,
                ), self.assertRaisesRegex(renderer.PacketError, "path changed"):
                    renderer.render_packet(
                        aggregate_receipt,
                        packet,
                        receipt,
                        authorized_pic_root=base.parent,
                    )
                self.assertTrue(substituted)
                self.assertFalse(packet.exists())
                self.assertFalse(receipt.exists())
            finally:
                if analysis.exists():
                    analysis.unlink()
                if anchor.exists():
                    os.rename(anchor, analysis)

    def test_retained_file_reader_rejects_interrupted_two_link_commit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            destination = root / "destination"
            source.write_bytes(b"sealed\n")
            source.chmod(0o444)
            os.link(source, destination)
            descriptor = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
            try:
                with self.assertRaisesRegex(renderer.PacketError, "singly linked"):
                    with renderer.ImmutableReadonlyPublicationFile(
                        descriptor, destination.name, "interrupted commit"
                    ):
                        pass
            finally:
                os.close(descriptor)

    def test_packet_receipt_is_the_only_accepted_marker(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            original = publisher._require_same_directory_at
            observed_pre_receipt_window = False

            def inspect_pre_receipt_window(
                parent_descriptor: int, name: str, descriptor: int, label: str
            ) -> None:
                nonlocal observed_pre_receipt_window
                original(parent_descriptor, name, descriptor, label)
                if label != "review-packet public build tree" or observed_pre_receipt_window:
                    return
                observed_pre_receipt_window = True
                self.assertTrue(packet.is_dir())
                self.assertFalse(receipt.exists())
                self.assertTrue(
                    (base / publisher._publication_guard_name(receipt.name)).is_file()
                )
                with self.assertRaisesRegex(
                    publisher.PressurePilotPublicationError, "fail-closed guard"
                ):
                    renderer.verify_published_review_packet_receipt(
                        receipt, authorized_pic_root=base.parent
                    )

            with patch.object(
                publisher,
                "_require_same_directory_at",
                side_effect=inspect_pre_receipt_window,
            ):
                renderer.render_packet(
                    aggregate_receipt,
                    packet,
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_pre_receipt_window)
            renderer.verify_published_review_packet_receipt(
                receipt, authorized_pic_root=base.parent
            )

    def test_lustre_file_fallback_publishes_complete_review_packet(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures(), patch.object(
            publisher.ctypes, "CDLL", return_value=object()
        ):
            receipt = base / "review-packet-receipt.json"
            renderer.render_packet(
                aggregate_receipt,
                base / "review-packet",
                receipt,
                authorized_pic_root=base.parent,
            )
            renderer.verify_published_review_packet_receipt(
                receipt, authorized_pic_root=base.parent
            )
            self.assertFalse(any(".staging-" in path.name for path in base.iterdir()))

    def test_late_publication_root_clone_is_rejected_by_retained_identity(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            anchor = base.parent / (base.name + ".descriptor-anchor")
            original = renderer._verify_published_review_packet_receipt
            cloned = False

            def clone_root(*args: object, **kwargs: object) -> dict[str, str]:
                nonlocal cloned
                verification = original(*args, **kwargs)
                if not cloned:
                    cloned = True
                    os.rename(base, anchor)
                    shutil.copytree(anchor, base)
                return verification

            with patch.object(
                renderer,
                "_verify_published_review_packet_receipt",
                side_effect=clone_root,
            ), self.assertRaisesRegex(
                renderer.PacketError,
                "publication root changed",
            ):
                renderer.render_packet(
                    aggregate_receipt,
                    packet,
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(cloned)
            self.assertTrue(receipt.is_file())
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError,
                "fail-closed guard",
            ):
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )

    def test_failure_does_not_delete_substituted_public_packet(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            moved = base / "review-packet.moved-original"

            def substitute_before_failure(*_args: object, **_kwargs: object) -> object:
                os.rename(packet, moved)
                shutil.copytree(moved, packet)
                raise renderer.PacketError("injected post-substitution packet failure")

            with patch.object(
                renderer,
                "_verify_published_review_packet_receipt",
                side_effect=substitute_before_failure,
            ), self.assertRaisesRegex(
                renderer.PacketError,
                "reviewed reconciliation required",
            ):
                renderer.render_packet(
                    aggregate_receipt,
                    packet,
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(packet.is_dir())
            self.assertTrue(moved.is_dir())
            self.assertTrue(receipt.is_file())
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "fail-closed guard"
            ):
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )

    def test_success_seal_rejects_packet_receipt_substitution_before_publication(
        self,
    ) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            moved = base / "review-packet-receipt.json.moved-original"
            original = publisher._publish_publication_seal_at

            def substitute_before_seal(
                acceptance_descriptor: int,
                publication_descriptor: int,
                receipt_name: str,
                receipt_payload: bytes,
                receipt_identity: tuple[int, int],
            ) -> None:
                os.rename(receipt, moved)
                shutil.copy2(moved, receipt)
                return original(
                    acceptance_descriptor,
                    publication_descriptor,
                    receipt_name,
                    receipt_payload,
                    receipt_identity,
                )

            with patch.object(
                publisher,
                "_publish_publication_seal_at",
                side_effect=substitute_before_seal,
            ), self.assertRaisesRegex(renderer.PacketError, "reviewed reconciliation required"):
                renderer.render_packet(
                    aggregate_receipt,
                    packet,
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(receipt.is_file())
            self.assertFalse(
                (
                    base.parent
                    / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY
                    / publisher._publication_seal_name(receipt.name)
                ).exists()
            )
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "fail-closed guard"
            ):
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )

    def test_coordinated_replacement_failure_leaves_packet_receipt_invalidated(
        self,
    ) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"

            def replace_every_public_artifact(
                *_args: object, **_kwargs: object
            ) -> object:
                for path in (packet, receipt):
                    moved = base / f"{path.name}.moved-original"
                    os.rename(path, moved)
                    if moved.is_dir():
                        shutil.copytree(moved, path)
                    else:
                        shutil.copy2(moved, path)
                raise renderer.PacketError(
                    "injected coordinated review-packet replacement failure"
                )

            with patch.object(
                renderer,
                "_verify_published_review_packet_receipt",
                side_effect=replace_every_public_artifact,
            ), self.assertRaisesRegex(
                renderer.PacketError, "reviewed reconciliation required"
            ):
                renderer.render_packet(
                    aggregate_receipt,
                    packet,
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(receipt.is_file())
            self.assertTrue(
                (base / publisher._publication_guard_name(receipt.name)).is_file()
            )
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "fail-closed guard"
            ):
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )

    def test_success_seal_rename_commit_is_reconciled_after_wrapper_raise(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            receipt = base / "review-packet-receipt.json"
            original = publisher._rename_no_replace_at
            observed_committed_seal = False

            def raise_after_seal_commit(
                parent_descriptor: int, source_name: str, destination_name: str
            ) -> None:
                nonlocal observed_committed_seal
                original(parent_descriptor, source_name, destination_name)
                if (
                    destination_name
                    != publisher._publication_seal_name(receipt.name)
                    or observed_committed_seal
                ):
                    return
                observed_committed_seal = True
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )
                raise renderer.PacketError("injected post-commit seal-wrapper failure")

            with patch.object(
                publisher, "_rename_no_replace_at", side_effect=raise_after_seal_commit
            ):
                rendered = renderer.render_packet(
                    aggregate_receipt,
                    base / "review-packet",
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_committed_seal)
            self.assertEqual(
                rendered["record_type"],
                "q011_section54_pressure_pilot_review_packet_receipt",
            )
            renderer.verify_published_review_packet_receipt(
                receipt, authorized_pic_root=base.parent
            )

    def test_success_seal_helper_wrapper_failure_rearms_guard(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            receipt = base / "review-packet-receipt.json"
            original = publisher._publish_publication_seal_at
            observed_committed_seal = False

            def raise_after_seal_helper_commit(
                acceptance_descriptor: int,
                publication_descriptor: int,
                receipt_name: str,
                receipt_payload: bytes,
                receipt_identity: tuple[int, int],
            ) -> None:
                nonlocal observed_committed_seal
                original(
                    acceptance_descriptor,
                    publication_descriptor,
                    receipt_name,
                    receipt_payload,
                    receipt_identity,
                )
                observed_committed_seal = True
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )
                raise renderer.PacketError("injected post-commit seal-helper failure")

            with patch.object(
                publisher,
                "_publish_publication_seal_at",
                side_effect=raise_after_seal_helper_commit,
            ), self.assertRaisesRegex(
                renderer.PacketError,
                "reviewed reconciliation required",
            ):
                renderer.render_packet(
                    aggregate_receipt,
                    base / "review-packet",
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_committed_seal)
            self.assertTrue(
                (base / publisher._publication_guard_name(receipt.name)).is_file()
            )
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "fail-closed guard"
            ):
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )

    def test_success_seal_is_committed_while_guard_remains_armed(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            receipt = base / "review-packet-receipt.json"
            original = publisher._publish_publication_seal_at
            observed_guard = False

            def inspect_guard_before_seal(
                acceptance_descriptor: int,
                publication_descriptor: int,
                receipt_name: str,
                receipt_payload: bytes,
                receipt_identity: tuple[int, int],
            ) -> tuple[int, int]:
                nonlocal observed_guard
                publisher._file_identity_at(
                    publication_descriptor,
                    publisher._publication_guard_name(receipt_name),
                    "receipt publication guard",
                )
                observed_guard = True
                return original(
                    acceptance_descriptor,
                    publication_descriptor,
                    receipt_name,
                    receipt_payload,
                    receipt_identity,
                )

            with patch.object(
                publisher,
                "_publish_publication_seal_at",
                side_effect=inspect_guard_before_seal,
            ):
                renderer.render_packet(
                    aggregate_receipt,
                    base / "review-packet",
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_guard)
            renderer.verify_published_review_packet_receipt(
                receipt, authorized_pic_root=base.parent
            )

    def test_guard_disarm_commit_is_reconciled_after_wrapper_raise(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            receipt = base / "review-packet-receipt.json"
            original = publisher._disarm_publication_guard_at

            def raise_after_disarm(
                parent_descriptor: int,
                receipt_name: str,
                guard_identity: tuple[int, int],
            ) -> None:
                original(parent_descriptor, receipt_name, guard_identity)
                raise renderer.PacketError(
                    "injected post-commit guard-disarm wrapper failure"
                )

            with patch.object(
                publisher,
                "_disarm_publication_guard_at",
                side_effect=raise_after_disarm,
            ):
                rendered = renderer.render_packet(
                    aggregate_receipt,
                    base / "review-packet",
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertEqual(
                rendered["record_type"],
                "q011_section54_pressure_pilot_review_packet_receipt",
            )
            renderer.verify_published_review_packet_receipt(
                receipt, authorized_pic_root=base.parent
            )

    def test_guard_unlink_without_parent_sync_is_reconciled_durably(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            receipt = base / "review-packet-receipt.json"
            guard_unlinked = False
            observed_retry_sync = False
            publication_descriptor: int | None = None
            original_fsync = publisher._fsync_descriptor
            original_require_same_file = publisher._require_same_file_at
            original_require_guard_absent = publisher._require_publication_guard_absent_at

            def unlink_then_raise(
                parent_descriptor: int,
                receipt_name: str,
                guard_identity: tuple[int, int],
            ) -> None:
                nonlocal guard_unlinked, publication_descriptor
                publication_descriptor = parent_descriptor
                guard_name = publisher._publication_guard_name(receipt_name)
                publisher._require_same_file_at(
                    parent_descriptor,
                    guard_name,
                    guard_identity,
                    "receipt publication guard",
                )
                os.unlink(guard_name, dir_fd=parent_descriptor)
                guard_unlinked = True
                raise OSError("injected failure before parent sync")

            def record_fsync(descriptor: int) -> None:
                nonlocal observed_retry_sync
                if guard_unlinked and descriptor == publication_descriptor:
                    observed_retry_sync = True
                original_fsync(descriptor)

            def require_same_file_after_sync(
                parent_descriptor: int,
                name: str,
                identity: tuple[int, int],
                label: str,
            ) -> None:
                if guard_unlinked:
                    self.assertTrue(observed_retry_sync)
                original_require_same_file(parent_descriptor, name, identity, label)

            def require_guard_absent_after_sync(
                parent_descriptor: int, receipt_name: str, label: str
            ) -> None:
                if guard_unlinked:
                    self.assertTrue(observed_retry_sync)
                original_require_guard_absent(parent_descriptor, receipt_name, label)

            with patch.object(
                publisher,
                "_disarm_publication_guard_at",
                side_effect=unlink_then_raise,
            ), patch.object(
                publisher,
                "_fsync_descriptor",
                side_effect=record_fsync,
            ), patch.object(
                publisher,
                "_require_same_file_at",
                side_effect=require_same_file_after_sync,
            ), patch.object(
                publisher,
                "_require_publication_guard_absent_at",
                side_effect=require_guard_absent_after_sync,
            ):
                renderer.render_packet(
                    aggregate_receipt,
                    base / "review-packet",
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_retry_sync)
            self.assertIsNotNone(publication_descriptor)
            renderer.verify_published_review_packet_receipt(
                receipt, authorized_pic_root=base.parent
            )

    def test_packet_publication_holds_transaction_lock_through_guard_disarm(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            receipt = base / "review-packet-receipt.json"
            acceptance = base.parent / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY
            anchor = publisher._publication_transaction_anchor(base.parent)
            observed_lock = False
            original = publisher._disarm_publication_guard_at

            def inspect_lock_before_disarm(
                parent_descriptor: int,
                receipt_name: str,
                guard_identity: tuple[int, int],
            ) -> None:
                nonlocal observed_lock
                for path in (anchor, acceptance):
                    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
                    try:
                        with self.assertRaises(BlockingIOError):
                            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    finally:
                        os.close(descriptor)
                observed_lock = True
                original(parent_descriptor, receipt_name, guard_identity)

            with patch.object(
                publisher,
                "_disarm_publication_guard_at",
                side_effect=inspect_lock_before_disarm,
            ):
                renderer.render_packet(
                    aggregate_receipt,
                    base / "review-packet",
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(observed_lock)

    def test_post_disarm_packet_receipt_substitution_fails_closed(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            receipt = base / "review-packet-receipt.json"
            moved = base / "review-packet-receipt.json.moved-original"
            original = publisher._disarm_publication_guard_at

            def substitute_after_disarm(
                parent_descriptor: int,
                receipt_name: str,
                guard_identity: tuple[int, int],
            ) -> None:
                original(parent_descriptor, receipt_name, guard_identity)
                os.rename(receipt, moved)
                shutil.copy2(moved, receipt)

            with patch.object(
                publisher,
                "_disarm_publication_guard_at",
                side_effect=substitute_after_disarm,
            ), self.assertRaisesRegex(renderer.PacketError, "reviewed reconciliation required"):
                renderer.render_packet(
                    aggregate_receipt,
                    base / "review-packet",
                    receipt,
                    authorized_pic_root=base.parent,
                )
            self.assertTrue(receipt.is_file())
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "fail-closed guard"
            ):
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )

    def test_packet_receipt_schema_version_rejects_json_boolean_alias(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            renderer.render_packet(
                aggregate_receipt,
                packet,
                receipt,
                authorized_pic_root=base.parent,
            )
            receipt.chmod(0o600)
            record = json.loads(receipt.read_text(encoding="utf-8"))
            record["schema_version"] = True
            receipt.write_bytes(renderer._json_bytes(record))
            receipt.chmod(0o444)
            fixture._rewrite_publication_seal(base, receipt)
            with self.assertRaisesRegex(renderer.PacketError, "receipt schema drifted"):
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )

    def test_packet_inventory_schema_version_rejects_json_boolean_alias(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            packet = base / "review-packet"
            receipt = base / "review-packet-receipt.json"
            renderer.render_packet(
                aggregate_receipt,
                packet,
                receipt,
                authorized_pic_root=base.parent,
            )
            fixture._make_writable(packet)
            inventory = packet / renderer.INVENTORY_NAME
            inventory_record = json.loads(inventory.read_text(encoding="utf-8"))
            inventory_record["schema_version"] = True
            inventory_payload = renderer._json_bytes(inventory_record)
            inventory.write_bytes(inventory_payload)
            receipt.chmod(0o600)
            receipt_record = json.loads(receipt.read_text(encoding="utf-8"))
            receipt_record["inventory_sha256"] = renderer._sha256(inventory_payload)
            receipt.write_bytes(renderer._json_bytes(receipt_record))
            receipt.chmod(0o444)
            fixture._rewrite_publication_seal(base, receipt)
            fixture._freeze_existing(packet)
            with self.assertRaisesRegex(renderer.PacketError, "inventory schema drifted"):
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )
            fixture._make_writable(packet)

    def test_retained_packet_receipt_requires_durable_success_seal(self) -> None:
        with _published_aggregate() as (base, aggregate_receipt), _fast_figures():
            receipt = base / "review-packet-receipt.json"
            renderer.render_packet(
                aggregate_receipt,
                base / "review-packet",
                receipt,
                authorized_pic_root=base.parent,
            )
            (
                base.parent
                / publisher.PUBLICATION_ACCEPTANCE_DIRECTORY
                / publisher._publication_seal_name(receipt.name)
            ).unlink()
            fixture._write_forged_publication_namespace_seal(base, receipt)
            with self.assertRaisesRegex(
                publisher.PressurePilotPublicationError, "durable success seal"
            ):
                renderer.verify_published_review_packet_receipt(
                    receipt, authorized_pic_root=base.parent
                )


if __name__ == "__main__":
    unittest.main()
