#!/usr/bin/env python3
"""Focused tests for Q011 pressure-publication acceptance-root provisioning."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import stat
import tempfile
import unittest
from unittest.mock import patch

from tst.publication import provision_q011_pressure_publication_acceptance_root as root


COMMIT = "a" * 40
HELPER_SHA256 = "b" * 64


def _identity(path: Path) -> tuple[int, int]:
    status = path.stat()
    return status.st_dev, status.st_ino


class Q011PressurePublicationAcceptanceRootTests(unittest.TestCase):
    def _provision(
        self,
        pic_root: Path,
        *,
        expected_uid: int | None = None,
        expected_gid: int | None = None,
        recover: bool = False,
        reconcile: bool = False,
        recovery_identity: tuple[int, int] | None = None,
        recovery_xattrs: set[str] | None = None,
    ) -> dict[str, object]:
        return root.provision_or_verify(
            pic_root=pic_root,
            expected_uid=os.getuid() if expected_uid is None else expected_uid,
            expected_gid=os.getgid() if expected_gid is None else expected_gid,
            validated_source_commit=COMMIT,
            helper_sha256=HELPER_SHA256,
            expected_pic_root_identity=_identity(pic_root),
            expected_recovery_identity=recovery_identity,
            expected_recovery_xattrs=recovery_xattrs,
            recover_exact_empty_inherited_setgid_root=recover,
            reconcile_exact_empty_normalized_root=reconcile,
        )

    def test_fresh_creation_under_setgid_parent_normalizes_exact_mode(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            pic_root.mkdir(mode=0o2770)
            pic_root.chmod(0o2770)
            result = self._provision(pic_root)
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            self.assertEqual(result["action"], "fresh_created_and_normalized")
            self.assertEqual(stat.S_IMODE(acceptance.stat().st_mode), 0o700)
            self.assertEqual(result["after"]["mode"], "0700")
            self.assertEqual(result["after"]["entries"], [])

    def test_fresh_creation_removes_inherited_acl_before_closure_check(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            pic_root.mkdir()
            names = {"system.posix_acl_default"}

            def remove_xattr(_descriptor: int, name: str) -> None:
                names.remove(name)

            with patch.object(
                root, "_xattrs", side_effect=lambda _fd: set(names)
            ), patch.object(root.os, "removexattr", side_effect=remove_xattr):
                result = self._provision(pic_root)
            self.assertEqual(result["action"], "fresh_created_and_normalized")
            self.assertEqual(result["after"]["xattrs"], [])

    def test_explicit_recovery_preserves_reviewed_inode_and_changes_mode(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            pic_root.mkdir()
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            acceptance.chmod(0o2700)
            identity = _identity(acceptance)
            result = self._provision(
                pic_root,
                recover=True,
                recovery_identity=identity,
                recovery_xattrs=set(),
            )
            self.assertEqual(
                result["action"], "recovered_exact_empty_inherited_setgid_root"
            )
            self.assertEqual(_identity(acceptance), identity)
            self.assertEqual(stat.S_IMODE(acceptance.stat().st_mode), 0o700)

    def test_interrupted_recovery_reconciles_exact_normalized_inode(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            pic_root.mkdir()
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            acceptance.chmod(0o2700)
            identity = _identity(acceptance)
            with patch.object(root.os, "fsync", side_effect=OSError("stop")):
                with self.assertRaisesRegex(OSError, "stop"):
                    self._provision(
                        pic_root,
                        recover=True,
                        recovery_identity=identity,
                        recovery_xattrs=set(),
                    )
            self.assertEqual(stat.S_IMODE(acceptance.stat().st_mode), 0o700)
            result = self._provision(
                pic_root,
                reconcile=True,
                recovery_identity=identity,
                recovery_xattrs=set(),
            )
            self.assertEqual(result["action"], "reconciled_exact_empty_normalized_root")

    def test_default_verification_rejects_recovery_and_nonempty_checkpoints(self) -> None:
        for state in ("setgid", "nonempty"):
            with self.subTest(state=state), tempfile.TemporaryDirectory() as directory:
                pic_root = Path(directory) / "pic"
                pic_root.mkdir()
                acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
                acceptance.mkdir()
                acceptance.chmod(0o2700 if state == "setgid" else 0o700)
                if state == "nonempty":
                    (acceptance / "unexpected").write_text("stop\n", encoding="utf-8")
                with self.assertRaises(root.AcceptanceRootError):
                    self._provision(pic_root)

    def test_default_verification_rejects_substituted_bound_checkpoint(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            pic_root.mkdir()
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            acceptance.chmod(0o700)
            device, inode = _identity(acceptance)
            with self.assertRaisesRegex(root.AcceptanceRootError, "identity drifted"):
                self._provision(
                    pic_root,
                    recovery_identity=(device, inode + 1),
                )

    def test_default_verification_rejects_absent_bound_checkpoint(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            pic_root.mkdir()
            with self.assertRaisesRegex(root.AcceptanceRootError, "unavailable"):
                self._provision(pic_root, recovery_identity=(1, 1))

    def test_default_verification_requires_exact_bound_xattrs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            pic_root.mkdir()
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            acceptance.chmod(0o700)
            with patch.object(root, "_xattrs", return_value=set()):
                with self.assertRaisesRegex(
                    root.AcceptanceRootError, "xattr closure differs"
                ):
                    self._provision(
                        pic_root,
                        recovery_identity=_identity(acceptance),
                        recovery_xattrs={"lustre.lov"},
                    )

    def test_recovery_rejects_absent_normalized_nonempty_and_replaced_roots(self) -> None:
        for state in ("absent", "normalized", "nonempty", "replaced"):
            with self.subTest(state=state), tempfile.TemporaryDirectory() as directory:
                pic_root = Path(directory) / "pic"
                pic_root.mkdir()
                acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
                if state != "absent":
                    acceptance.mkdir()
                    acceptance.chmod(0o700 if state == "normalized" else 0o2700)
                identity = _identity(acceptance) if acceptance.exists() else (1, 1)
                if state == "nonempty":
                    (acceptance / "unexpected").write_text("stop\n", encoding="utf-8")
                if state == "replaced":
                    retained = Path(directory) / "retained"
                    acceptance.rename(retained)
                    acceptance.mkdir()
                    acceptance.chmod(0o2700)
                with self.assertRaises(root.AcceptanceRootError):
                    self._provision(
                        pic_root,
                        recover=True,
                        recovery_identity=identity,
                        recovery_xattrs=set(),
                    )

    def test_retained_open_rejects_stat_to_open_child_swap(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            pic_root.mkdir()
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            original = acceptance.stat()
            acceptance.rename(Path(directory) / "retained")
            acceptance.mkdir()
            descriptor = os.open(pic_root, root.DIRECTORY_FLAGS)
            try:
                with self.assertRaisesRegex(root.AcceptanceRootError, "before retained"):
                    root._open_retained_child(
                        descriptor, root.ACCEPTANCE_DIRECTORY, original
                    )
            finally:
                os.close(descriptor)

    def test_rejects_wrong_owner_group_mode_symlink_and_unexpected_xattr(self) -> None:
        for state in ("owner", "group", "mode", "symlink", "xattr"):
            with self.subTest(state=state), tempfile.TemporaryDirectory() as directory:
                pic_root = Path(directory) / "pic"
                pic_root.mkdir()
                acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
                if state == "symlink":
                    target = Path(directory) / "target"
                    target.mkdir()
                    acceptance.symlink_to(target, target_is_directory=True)
                else:
                    acceptance.mkdir()
                    acceptance.chmod(0o700 if state != "mode" else 0o750)
                expected_uid = os.getuid() + 1 if state == "owner" else os.getuid()
                expected_gid = os.getgid() + 1 if state == "group" else os.getgid()
                context = (
                    patch.object(root, "_xattrs", return_value={"user.unexpected"})
                    if state == "xattr"
                    else patch.object(root, "_xattrs", wraps=root._xattrs)
                )
                with context, self.assertRaises((root.AcceptanceRootError, OSError)):
                    self._provision(
                        pic_root,
                        expected_uid=expected_uid,
                        expected_gid=expected_gid,
                    )

    def test_recovery_rejects_acl_and_requires_exact_reviewed_xattrs(self) -> None:
        for names in ({"system.posix_acl_default"}, set()):
            with self.subTest(names=names), tempfile.TemporaryDirectory() as directory:
                pic_root = Path(directory) / "pic"
                pic_root.mkdir()
                acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
                acceptance.mkdir()
                acceptance.chmod(0o2700)
                with patch.object(root, "_xattrs", return_value=names):
                    with self.assertRaises(root.AcceptanceRootError):
                        self._provision(
                            pic_root,
                            recover=True,
                            recovery_identity=_identity(acceptance),
                            recovery_xattrs={"lustre.lov"},
                        )

    def test_rechecks_identity_and_syncs_child_and_parent(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            pic_root.mkdir()
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            acceptance.chmod(0o700)
            real_fsync = os.fsync
            synced = []

            def record_fsync(descriptor: int) -> None:
                synced.append(os.fstat(descriptor).st_ino)
                real_fsync(descriptor)

            with patch.object(root.os, "fsync", side_effect=record_fsync):
                result = self._provision(pic_root)
            self.assertEqual(result["action"], "verified_existing_empty_root")
            self.assertEqual(synced, [acceptance.stat().st_ino, pic_root.stat().st_ino])

    def test_canonical_receipt_validation_and_descriptor_relative_publication(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            policy = pic_root / "policy"
            policy.mkdir(parents=True)
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            acceptance.chmod(0o2700)
            identity = _identity(acceptance)
            with patch.multiple(
                root,
                AUTHORIZED_PIC_ROOT=pic_root,
                REVIEWED_PIC_ROOT_IDENTITY=_identity(pic_root),
                REVIEWED_POLICY_ROOT_IDENTITY=_identity(policy),
                REVIEWED_ACCEPTANCE_ROOT_IDENTITY=identity,
                REVIEWED_PIC_ROOT_MODE=stat.S_IMODE(pic_root.stat().st_mode),
                REVIEWED_POLICY_ROOT_MODE=stat.S_IMODE(policy.stat().st_mode),
                REVIEWED_PIC_ROOT_XATTRS=set(),
                REVIEWED_POLICY_ROOT_XATTRS=set(),
                AUTHORIZED_UID=os.getuid(),
                AUTHORIZED_GID=os.getgid(),
                REVIEWED_RECOVERY_XATTRS=set(),
            ):
                receipt = self._provision(
                    pic_root,
                    recover=True,
                    recovery_identity=identity,
                    recovery_xattrs=set(),
                )
                staging = policy / ".q011-acceptance-recovery.test"
                staging.write_bytes(root.canonical_json_bytes(receipt))
                root.publish_recovery_receipt(
                    staging,
                    validated_source_commit=COMMIT,
                    helper_sha256=HELPER_SHA256,
                )
                published = policy / root.RECOVERY_RECEIPT_NAME
                self.assertFalse(staging.exists())
                self.assertEqual(stat.S_IMODE(published.stat().st_mode), 0o400)
                root.verify_recovery_receipt(
                    published,
                    validated_source_commit=COMMIT,
                    helper_sha256=HELPER_SHA256,
                )
                published.chmod(0o600)
                with self.assertRaisesRegex(root.AcceptanceRootError, "mode"):
                    root.verify_recovery_receipt(
                        published,
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )
                published.chmod(0o400)
                with patch.object(root, "AUTHORIZED_UID", os.getuid() + 1):
                    with self.assertRaisesRegex(root.AcceptanceRootError, "owner"):
                        root.verify_recovery_receipt(
                            published,
                            validated_source_commit=COMMIT,
                            helper_sha256=HELPER_SHA256,
                        )
                forged = root.canonical_json_bytes({**receipt, "action": "forged"})
                with self.assertRaisesRegex(root.AcceptanceRootError, "action"):
                    root.validate_recovery_receipt(
                        forged,
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )

    def test_receipt_verification_durably_closes_post_link_interruption(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            policy = pic_root / "policy"
            policy.mkdir(parents=True)
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            acceptance.chmod(0o2700)
            identity = _identity(acceptance)
            with patch.multiple(
                root,
                AUTHORIZED_PIC_ROOT=pic_root,
                REVIEWED_PIC_ROOT_IDENTITY=_identity(pic_root),
                REVIEWED_POLICY_ROOT_IDENTITY=_identity(policy),
                REVIEWED_ACCEPTANCE_ROOT_IDENTITY=identity,
                REVIEWED_PIC_ROOT_MODE=stat.S_IMODE(pic_root.stat().st_mode),
                REVIEWED_POLICY_ROOT_MODE=stat.S_IMODE(policy.stat().st_mode),
                REVIEWED_PIC_ROOT_XATTRS=set(),
                REVIEWED_POLICY_ROOT_XATTRS=set(),
                AUTHORIZED_UID=os.getuid(),
                AUTHORIZED_GID=os.getgid(),
                REVIEWED_RECOVERY_XATTRS=set(),
            ):
                receipt = self._provision(
                    pic_root,
                    recover=True,
                    recovery_identity=identity,
                    recovery_xattrs=set(),
                )
                staging = policy / ".q011-acceptance-recovery.test"
                staging.write_bytes(root.canonical_json_bytes(receipt))
                real_fsync = os.fsync

                def stop_before_parent_sync(descriptor: int) -> None:
                    if os.fstat(descriptor).st_ino == policy.stat().st_ino:
                        raise OSError("stop before parent sync")
                    real_fsync(descriptor)

                with patch.object(root.os, "fsync", side_effect=stop_before_parent_sync):
                    with self.assertRaisesRegex(OSError, "stop before parent sync"):
                        root.publish_recovery_receipt(
                            staging,
                            validated_source_commit=COMMIT,
                            helper_sha256=HELPER_SHA256,
                        )
                published = policy / root.RECOVERY_RECEIPT_NAME
                self.assertTrue(published.exists())
                self.assertEqual(stat.S_IMODE(published.stat().st_mode), 0o400)
                root.verify_recovery_receipt(
                    published,
                    validated_source_commit=COMMIT,
                    helper_sha256=HELPER_SHA256,
                )

    def test_linked_receipt_alias_requires_explicit_reconciliation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            policy = pic_root / "policy"
            policy.mkdir(parents=True)
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            acceptance.chmod(0o2700)
            identity = _identity(acceptance)
            with patch.multiple(
                root,
                AUTHORIZED_PIC_ROOT=pic_root,
                REVIEWED_PIC_ROOT_IDENTITY=_identity(pic_root),
                REVIEWED_POLICY_ROOT_IDENTITY=_identity(policy),
                REVIEWED_ACCEPTANCE_ROOT_IDENTITY=identity,
                REVIEWED_PIC_ROOT_MODE=stat.S_IMODE(pic_root.stat().st_mode),
                REVIEWED_POLICY_ROOT_MODE=stat.S_IMODE(policy.stat().st_mode),
                REVIEWED_PIC_ROOT_XATTRS=set(),
                REVIEWED_POLICY_ROOT_XATTRS=set(),
                AUTHORIZED_UID=os.getuid(),
                AUTHORIZED_GID=os.getgid(),
                REVIEWED_RECOVERY_XATTRS=set(),
            ):
                receipt = self._provision(
                    pic_root,
                    recover=True,
                    recovery_identity=identity,
                    recovery_xattrs=set(),
                )
                staging = policy / ".q011-acceptance-recovery.test"
                staging.write_bytes(root.canonical_json_bytes(receipt))
                real_open = os.open

                def stop_after_link(path: str, *args: object, **kwargs: object) -> int:
                    if path == root.RECOVERY_RECEIPT_NAME:
                        raise OSError("stop after link")
                    return real_open(path, *args, **kwargs)

                with patch.object(root.os, "open", side_effect=stop_after_link):
                    with self.assertRaisesRegex(OSError, "stop after link"):
                        root.publish_recovery_receipt(
                            staging,
                            validated_source_commit=COMMIT,
                            helper_sha256=HELPER_SHA256,
                        )
                published = policy / root.RECOVERY_RECEIPT_NAME
                self.assertTrue(staging.exists())
                self.assertTrue(published.exists())
                with self.assertRaisesRegex(root.AcceptanceRootError, "namespace"):
                    root.verify_recovery_receipt(
                        published,
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )
                acceptance.chmod(0o777)
                with self.assertRaisesRegex(root.AcceptanceRootError, "mode"):
                    root.reconcile_linked_recovery_receipt(
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )
                acceptance.chmod(0o700)
                external = Path(directory) / "external-receipt-link"
                os.link(published, external)
                with self.assertRaisesRegex(root.AcceptanceRootError, "link count"):
                    root.reconcile_linked_recovery_receipt(
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )
                external.unlink()
                retained = Path(directory) / "retained-staging-link"
                staging.rename(retained)
                with self.assertRaisesRegex(root.AcceptanceRootError, "namespace"):
                    root.reconcile_linked_recovery_receipt(
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )
                wrong_alias = policy / ".q011-acceptance-recovery.wrong"
                wrong_alias.write_text("wrong\n", encoding="utf-8")
                with self.assertRaisesRegex(root.AcceptanceRootError, "retained staging alias"):
                    root.reconcile_linked_recovery_receipt(
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )
                wrong_alias.unlink()
                retained.rename(staging)
                root.reconcile_linked_recovery_receipt(
                    validated_source_commit=COMMIT,
                    helper_sha256=HELPER_SHA256,
                )
                self.assertFalse(staging.exists())
                root.verify_recovery_receipt(
                    published,
                    validated_source_commit=COMMIT,
                    helper_sha256=HELPER_SHA256,
                )

    def test_policy_receipt_operations_reject_unreviewed_pic_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            (pic_root / "policy").mkdir(parents=True)
            with patch.multiple(
                root,
                AUTHORIZED_PIC_ROOT=pic_root,
                REVIEWED_PIC_ROOT_IDENTITY=(1, 1),
            ):
                with self.assertRaisesRegex(root.AcceptanceRootError, "identity drifted"):
                    root._open_policy_descriptors()

    def test_policy_receipt_operations_reject_replaced_policy_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            policy = pic_root / "policy"
            policy.mkdir(parents=True)
            policy_identity = _identity(policy)
            policy.rename(pic_root / "retained-policy")
            policy.mkdir()
            with patch.multiple(
                root,
                AUTHORIZED_PIC_ROOT=pic_root,
                REVIEWED_PIC_ROOT_IDENTITY=_identity(pic_root),
                REVIEWED_POLICY_ROOT_IDENTITY=policy_identity,
                REVIEWED_PIC_ROOT_MODE=stat.S_IMODE(pic_root.stat().st_mode),
                REVIEWED_PIC_ROOT_XATTRS=set(),
                AUTHORIZED_UID=os.getuid(),
                AUTHORIZED_GID=os.getgid(),
            ):
                with self.assertRaisesRegex(root.AcceptanceRootError, "identity drifted"):
                    root._open_policy_descriptors()

    def test_receipt_operations_reject_live_acceptance_root_drift(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            policy = pic_root / "policy"
            policy.mkdir(parents=True)
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            acceptance.chmod(0o2700)
            identity = _identity(acceptance)
            with patch.multiple(
                root,
                AUTHORIZED_PIC_ROOT=pic_root,
                REVIEWED_PIC_ROOT_IDENTITY=_identity(pic_root),
                REVIEWED_POLICY_ROOT_IDENTITY=_identity(policy),
                REVIEWED_ACCEPTANCE_ROOT_IDENTITY=identity,
                REVIEWED_PIC_ROOT_MODE=stat.S_IMODE(pic_root.stat().st_mode),
                REVIEWED_POLICY_ROOT_MODE=stat.S_IMODE(policy.stat().st_mode),
                REVIEWED_PIC_ROOT_XATTRS=set(),
                REVIEWED_POLICY_ROOT_XATTRS=set(),
                AUTHORIZED_UID=os.getuid(),
                AUTHORIZED_GID=os.getgid(),
                REVIEWED_RECOVERY_XATTRS=set(),
            ):
                receipt = self._provision(
                    pic_root,
                    recover=True,
                    recovery_identity=identity,
                    recovery_xattrs=set(),
                )
                staging = policy / ".q011-acceptance-recovery.test"
                staging.write_bytes(root.canonical_json_bytes(receipt))
                acceptance.chmod(0o777)
                with self.assertRaisesRegex(root.AcceptanceRootError, "mode"):
                    root.publish_recovery_receipt(
                        staging,
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )
                acceptance.chmod(0o700)
                root.publish_recovery_receipt(
                    staging,
                    validated_source_commit=COMMIT,
                    helper_sha256=HELPER_SHA256,
                )
                published = policy / root.RECOVERY_RECEIPT_NAME
                acceptance.chmod(0o777)
                with self.assertRaisesRegex(root.AcceptanceRootError, "mode"):
                    root.verify_recovery_receipt(
                        published,
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )
                acceptance.chmod(0o700)
                (acceptance / "unexpected").write_text("stop\n", encoding="utf-8")
                with self.assertRaisesRegex(root.AcceptanceRootError, "not empty"):
                    root.verify_recovery_receipt(
                        published,
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )

    def test_receipt_operations_reject_weakened_policy_and_receipt_xattrs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            policy = pic_root / "policy"
            policy.mkdir(parents=True)
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            acceptance.chmod(0o2700)
            identity = _identity(acceptance)
            reviewed_policy_mode = stat.S_IMODE(policy.stat().st_mode)
            with patch.multiple(
                root,
                AUTHORIZED_PIC_ROOT=pic_root,
                REVIEWED_PIC_ROOT_IDENTITY=_identity(pic_root),
                REVIEWED_POLICY_ROOT_IDENTITY=_identity(policy),
                REVIEWED_ACCEPTANCE_ROOT_IDENTITY=identity,
                REVIEWED_PIC_ROOT_MODE=stat.S_IMODE(pic_root.stat().st_mode),
                REVIEWED_POLICY_ROOT_MODE=reviewed_policy_mode,
                REVIEWED_PIC_ROOT_XATTRS=set(),
                REVIEWED_POLICY_ROOT_XATTRS=set(),
                AUTHORIZED_UID=os.getuid(),
                AUTHORIZED_GID=os.getgid(),
                REVIEWED_RECOVERY_XATTRS=set(),
            ):
                receipt = self._provision(
                    pic_root,
                    recover=True,
                    recovery_identity=identity,
                    recovery_xattrs=set(),
                )
                staging = policy / ".q011-acceptance-recovery.test"
                staging.write_bytes(root.canonical_json_bytes(receipt))
                policy.chmod(0o777)
                with self.assertRaisesRegex(root.AcceptanceRootError, "PIC policy root mode"):
                    root.publish_recovery_receipt(
                        staging,
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )
                policy.chmod(reviewed_policy_mode)

                staging_inode = staging.stat().st_ino

                def receipt_acl(descriptor: int) -> set[str]:
                    if os.fstat(descriptor).st_ino == staging_inode:
                        return {"system.posix_acl_access"}
                    return set()

                with patch.object(root, "_xattrs", side_effect=receipt_acl):
                    with self.assertRaisesRegex(root.AcceptanceRootError, "retains ACL"):
                        root.publish_recovery_receipt(
                            staging,
                            validated_source_commit=COMMIT,
                            helper_sha256=HELPER_SHA256,
                        )

                def receipt_xattr(descriptor: int) -> set[str]:
                    if os.fstat(descriptor).st_ino == staging_inode:
                        return {"user.unexpected"}
                    return set()

                with patch.object(root, "_xattrs", side_effect=receipt_xattr):
                    with self.assertRaisesRegex(root.AcceptanceRootError, "unexpected xattrs"):
                        root.publish_recovery_receipt(
                            staging,
                            validated_source_commit=COMMIT,
                            helper_sha256=HELPER_SHA256,
                        )

    def test_recovery_namespace_preflight_rejects_orphan_staging_alias(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            policy = pic_root / "policy"
            policy.mkdir(parents=True)
            (policy / ".q011-acceptance-recovery.orphan").write_text(
                "orphan\n", encoding="utf-8"
            )
            with patch.multiple(
                root,
                AUTHORIZED_PIC_ROOT=pic_root,
                REVIEWED_PIC_ROOT_IDENTITY=_identity(pic_root),
                REVIEWED_POLICY_ROOT_IDENTITY=_identity(policy),
                REVIEWED_PIC_ROOT_MODE=stat.S_IMODE(pic_root.stat().st_mode),
                REVIEWED_POLICY_ROOT_MODE=stat.S_IMODE(policy.stat().st_mode),
                REVIEWED_PIC_ROOT_XATTRS=set(),
                REVIEWED_POLICY_ROOT_XATTRS=set(),
                AUTHORIZED_UID=os.getuid(),
                AUTHORIZED_GID=os.getgid(),
            ):
                with self.assertRaisesRegex(root.AcceptanceRootError, "not empty"):
                    root._require_empty_recovery_receipt_namespace()

    def test_publication_rejects_additional_orphan_staging_alias(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pic_root = Path(directory) / "pic"
            policy = pic_root / "policy"
            policy.mkdir(parents=True)
            acceptance = pic_root / root.ACCEPTANCE_DIRECTORY
            acceptance.mkdir()
            acceptance.chmod(0o2700)
            identity = _identity(acceptance)
            with patch.multiple(
                root,
                AUTHORIZED_PIC_ROOT=pic_root,
                REVIEWED_PIC_ROOT_IDENTITY=_identity(pic_root),
                REVIEWED_POLICY_ROOT_IDENTITY=_identity(policy),
                REVIEWED_ACCEPTANCE_ROOT_IDENTITY=identity,
                REVIEWED_PIC_ROOT_MODE=stat.S_IMODE(pic_root.stat().st_mode),
                REVIEWED_POLICY_ROOT_MODE=stat.S_IMODE(policy.stat().st_mode),
                REVIEWED_PIC_ROOT_XATTRS=set(),
                REVIEWED_POLICY_ROOT_XATTRS=set(),
                AUTHORIZED_UID=os.getuid(),
                AUTHORIZED_GID=os.getgid(),
                REVIEWED_RECOVERY_XATTRS=set(),
            ):
                receipt = self._provision(
                    pic_root,
                    recover=True,
                    recovery_identity=identity,
                    recovery_xattrs=set(),
                )
                staging = policy / ".q011-acceptance-recovery.test"
                staging.write_bytes(root.canonical_json_bytes(receipt))
                (policy / ".q011-acceptance-recovery.orphan").write_text(
                    "orphan\n", encoding="utf-8"
                )
                with self.assertRaisesRegex(root.AcceptanceRootError, "not exact"):
                    root.publish_recovery_receipt(
                        staging,
                        validated_source_commit=COMMIT,
                        helper_sha256=HELPER_SHA256,
                    )

    def test_self_authentication_rejects_malformed_bindings_and_missing_nofollow(
        self,
    ) -> None:
        with self.assertRaisesRegex(root.AcceptanceRootError, "Git SHA-1"):
            root._self_authenticate(
                validated_source_commit="A" * 40,
                expected_helper_sha256=HELPER_SHA256,
            )
        with self.assertRaisesRegex(root.AcceptanceRootError, "SHA-256"):
            root._self_authenticate(
                validated_source_commit=COMMIT,
                expected_helper_sha256="b" * 63,
            )
        with self.assertRaisesRegex(root.AcceptanceRootError, "SHA-256 drifted"):
            root._self_authenticate(
                validated_source_commit=COMMIT,
                expected_helper_sha256=HELPER_SHA256,
            )
        with patch.object(root, "NOFOLLOW_FLAG", None):
            with self.assertRaisesRegex(root.AcceptanceRootError, "O_NOFOLLOW"):
                root._open_absolute_directory(Path("/"))

    def test_self_authentication_and_cli_have_no_root_override(self) -> None:
        payload = Path(root.__file__).read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        self.assertEqual(
            root._self_authenticate(
                validated_source_commit=COMMIT,
                expected_helper_sha256=digest,
            ),
            digest,
        )
        with self.assertRaises(SystemExit):
            root.main(["--pic-root", "/tmp/not-authorized"])


if __name__ == "__main__":
    unittest.main()
