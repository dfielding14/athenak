#!/usr/bin/env python3
"""Tests for immutable Frontier PIC submission snapshots."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
import copy
from datetime import datetime, timedelta, timezone
import fcntl
import hashlib
import io
import json
import os
from pathlib import Path
import pwd
import inspect
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
from typing import Callable, Iterator
import unittest
from unittest.mock import call, patch
import uuid

import control_plane_common
import install_control_plane
import launch_trampoline
import ledger
import promote_active_policy
import q011_pressure_review_packet_verifier as pressure_packet_verifier
import reconcile_frontier_job
import reconcile_manual_frontier_allocations
import revalidate_clean_candidate
import terminal_recovery_handoff
import validate_and_reserve_frontier_job
from control_plane_common import atomic_write_bytes, durable_mkdir_parents
from control_plane_common import AUTHORIZED_STORAGE_PREFLIGHT_CAPTURE_SOURCE_BLOBS
from control_plane_common import AUTHORIZED_STORAGE_PREFLIGHT_OPERATIONS
from control_plane_common import (
    AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_SOURCE_AUTHENTICATION,
)
from control_plane_common import PinnedDirectoryAncestry
from control_plane_common import CONTROL_PLANE_FILES, inventory_digest, make_tree_read_only
from control_plane_common import durable_replace_tree
from control_plane_common import git_archive_commit_from_bytes
from control_plane_common import git_commit_tree_from_bytes
from control_plane_common import git_tree_sha1_from_archive
from control_plane_common import read_stable_regular_file_below, remove_tree
from control_plane_common import require_ledger_paths
from control_plane_common import launch_contract_sha256, record_for_role, sha256
from control_plane_common import prepared_artifact_manifest_from_source_archive
from control_plane_common import source_bundle_sha256
from control_plane_common import scheduler_account_matches_authorized
from control_plane_common import stable_serialization_anchor
from control_plane_common import trusted_git_command, trusted_git_environment
from control_plane_common import trusted_slurm_environment
from control_plane_common import require_storage_policy_unlock_snapshot
from control_plane_common import validate_clean_candidate_bundle
from control_plane_common import validate_planner_retention_binding
from control_plane_common import validate_launch_contract, verify_installed_control_plane
from control_plane_common import verify_historical_installed_control_plane
from control_plane_common import verify_snapshot_files
from control_plane_common import PRODUCTION_RUNTIME_MODULEPATH
from control_plane_common import TRUSTED_GIT, TRUSTED_PYTHON
from control_plane_common import TRUSTED_SACCT, TRUSTED_SBATCH, TRUSTED_SCANCEL
from control_plane_common import TRUSTED_SCONTROL, TRUSTED_SQUEUE
from create_clean_candidate_freeze import _authorized_source_path
from create_clean_candidate_freeze import _source_identity, _validated_submodules
from create_clean_candidate_freeze import create_freeze
from create_pre_submit_manifest import create_manifest
from initialize_frontier_ledger import initialize_from_policy
from install_control_plane import install
from launch_trampoline import _freeze_artifact_file_at
from launch_trampoline import _freeze_artifact_tree_at
from launch_trampoline import _capture_artifact_directory_identities_at
from launch_trampoline import _publish_frozen_artifact_inventory, _TASK_LOCAL_EXEC, launch
from ledger import accounting, append_primary_event, genesis_anchor_paths
from ledger import incomplete_manual_accounting_marker_paths
from ledger import validate_primary_chain
from operator_attestation import OPERATOR_STATEMENTS
from promote_active_policy import _promotion_lock, promote
from reconcile_frontier_job import reconcile
from reconcile_manual_frontier_allocations import reconcile_manual_allocations
from terminal_recovery_handoff import create_handoff
from terminal_recovery_handoff import require_closed_received_marker
from validate_and_reserve_frontier_job import _clear_matching_pending_marker
from validate_and_reserve_frontier_job import _require_current_reservation_marker
from validate_and_reserve_frontier_job import _require_scheduler_output_path
from validate_and_reserve_frontier_job import executable_reservation_bound_manifest
from validate_and_reserve_frontier_job import mark_dispatch_started, mark_submitted
from validate_and_reserve_frontier_job import repair_ledger_mirror
from validate_and_reserve_frontier_job import repair_reservation_attachments
from validate_and_reserve_frontier_job import reservation_bound_manifest
from validate_and_reserve_frontier_job import reserve, transition
from verify_compute_node_snapshot import verify
from write_orion_build_profile import _clone_exact_submodule
from write_orion_build_profile import _clone_fresh_source
from write_orion_build_profile import _require_exact_clone_capabilities
from write_orion_build_profile import _require_exact_standalone_clone
from write_orion_build_profile import _require_exact_standalone_full_clone
from write_orion_build_profile import _require_exact_standalone_shallow_clone
from write_orion_build_profile import _submodule_status
from write_orion_build_profile import build_profile as build_orion_profile
from write_orion_build_profile import write_profile


class SnapshotTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.pic_root = self.root / "pic"
        self.project_home_root = self.root / "project_home"
        self.sources = self.root / "sources"
        self.sources.mkdir()
        self.config = self.root / "config.json"
        self.policy = self.root / "storage_policy.json"
        self.submission_id = "804dca3d-f89f-4357-9407-e59804961ad7"
        self.authorized_clean_candidate_source_root: Path | None = None
        self.registered_science_slices: list[dict[str, object]] = []
        self.science_submission_freeze: dict[str, object] | None = None
        self.control_plane_dir = install(self.pic_root)
        self.control_plane_version = self.control_plane_dir.name
        self.project_home_control_plane_dir = install(self.project_home_root)
        self.assertEqual(
            self.project_home_control_plane_dir.name, self.control_plane_version
        )
        self._write(
            "job.sh",
            "#!/bin/bash\n#SBATCH -A AST207\n#SBATCH -p batch\n#SBATCH -q debug\n"
            f"#SBATCH -o {self.pic_root}/logs/slurm/%x.%j.log\n"
            "#SBATCH -N 1\n#SBATCH -t 00:10:00\n",
        )
        self._write("athena", "executable placeholder\n")
        self._write("input.athinput", "<job>\nbasename = pic\n")
        self._write(
            "environment.sh",
            Path(__file__).with_name("frontier_pic_environment.sh").read_text(
                encoding="utf-8"
            ),
        )
        self._write_timeout()
        self._write("queue.txt", "")
        real_check_output = subprocess.check_output

        def check_output(command: list[str], *args: object, **kwargs: object) -> str:
            if command[0] == TRUSTED_SQUEUE:
                return (self.sources / "queue.txt").read_text(encoding="utf-8")
            return real_check_output(command, *args, **kwargs)

        self.queue_output_patcher = patch(
            "validate_and_reserve_frontier_job.subprocess.check_output",
            side_effect=check_output,
        )
        self.queue_output_patcher.start()
        self.addCleanup(self.queue_output_patcher.stop)
        self._write("analysis.py", "print('analysis')\n")
        self._write_config()
        self._write_policy()
        self._promote_policy()
        self.ledger = self.pic_root / "ledger" / "node_hours.jsonl"
        self.csv = self.pic_root / "ledger" / "node_hours.csv"
        self.receipts = self.pic_root / "ledger" / "mirror_receipts.jsonl"
        self.mirror = self.project_home_root / "ledger" / "node_hours.jsonl"
        genesis = initialize_from_policy(
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            mirror_receipts=self.receipts,
            mirror_jsonl=self.mirror,
            mirror_transport="filesystem_copy",
            notes="snapshot control plane test",
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        receipt = json.loads(self.receipts.read_text(encoding="utf-8").splitlines()[0])
        self._closed_genesis = {
            "status": "initialized",
            "timestamp": genesis["timestamp"],
            "control_plane_version": genesis["control_plane_version"],
            "event_sha256": genesis["event_sha256"],
            "mirror_transport": "filesystem_copy",
            "mirror_ack_sha256": receipt["mirror_ack_sha256"],
        }
        self._write_policy()
        self._promote_policy()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write(self, name: str, text: str) -> Path:
        path = self.sources / name
        path.write_text(text, encoding="utf-8")
        return path

    def _utc(self, value: datetime) -> str:
        return value.replace(microsecond=0).isoformat().replace("+00:00", "Z")

    def _sealed_operator_attestation(
        self,
        authorization_id: str,
        phase: str,
        *,
        control_plane_version: str | None = None,
    ) -> Path:
        version = control_plane_version or self.control_plane_version
        archive = self.pic_root / "operator_attestations"
        archive.mkdir(mode=0o700, exist_ok=True)
        now = datetime.now(timezone.utc).replace(microsecond=0)
        for offset in range(10):
            timestamp = now - timedelta(seconds=offset)
            root = archive / (
                f"{timestamp.strftime('%Y%m%dT%H%M%SZ')}-{authorization_id}-{phase}"
            )
            try:
                root.mkdir(mode=0o700)
            except FileExistsError:
                existing = root / "attestation.json"
                if existing.is_file():
                    value = json.loads(existing.read_text(encoding="utf-8"))
                    if value.get("control_plane_version") == version:
                        return existing
                continue
            break
        else:
            raise AssertionError("could not allocate sealed-attestation fixture")
        manual_values = {
            str(self.pic_root / "ledger/pending_manual_accounting.json"): "absent",
            str(self.project_home_root / "ledger/pending_manual_accounting.json"): "absent",
        }
        counts = {
            str(self.pic_root / "ledger/node_hours.jsonl"): 1,
            str(self.pic_root / "ledger/mirror_receipts.jsonl"): 1,
            str(self.project_home_root / "ledger/node_hours.jsonl"): 1,
        }
        state = {
            "validation": "coherent",
            "validator": "existing_read_only_mirrored_state_snapshot",
            "ledger_record_count": 1,
            "ledger_tail_event_sha256": "0" * 64,
            "active_reservation_count": 0,
            "active_reservation_ids": [],
        }
        payloads = {
            "same_account_process_snapshot.txt": b"fixture same-account snapshot\n",
            "queue_snapshot.txt": b"",
            "pending_submission_marker.txt": b"absent\n",
            "pending_manual_accounting_marker.txt": "".join(
                f"{path} {value}\n" for path, value in manual_values.items()
            ).encode("utf-8"),
            "mirrored_ledger_line_counts.txt": "".join(
                f"{value} {path}\n" for path, value in counts.items()
            ).encode("utf-8"),
            "validated_mirrored_ledger_state.json": (
                json.dumps(state, indent=2, sort_keys=True) + "\n"
            ).encode("utf-8"),
        }
        payloads["capture_same_account_process_snapshot.txt"] = payloads[
            "same_account_process_snapshot.txt"
        ]
        for filename, payload in payloads.items():
            (root / filename).write_bytes(payload)
            if (
                filename != "same_account_process_snapshot.txt"
                and not filename.startswith("capture_")
            ):
                (root / f"capture_{filename}").write_bytes(payload)

        def record(filename: str, **values: object) -> dict[str, object]:
            return {
                "path": filename,
                "sha256": hashlib.sha256(payloads[filename]).hexdigest(),
                **values,
            }

        attestation = {
            "schema_version": 1,
            "record_type": (
                "q027_frontier_registered_science_same_account_isolation_attestation"
            ),
            "recorded_utc": self._utc(timestamp),
            "sealed_utc": self._utc(timestamp),
            "registered_science_authorization_id": authorization_id,
            "control_plane_version": version,
            "phase": phase,
            "same_account_process_snapshot": record(
                "same_account_process_snapshot.txt"
            ),
            "queue_snapshot": record("queue_snapshot.txt"),
            "pending_submission_marker": record(
                "pending_submission_marker.txt", value="absent"
            ),
            "pending_manual_accounting_marker": record(
                "pending_manual_accounting_marker.txt",
                value="absent",
                values=manual_values,
            ),
            "mirrored_ledger_line_counts": record(
                "mirrored_ledger_line_counts.txt", counts=counts
            ),
            "validated_mirrored_ledger_state": record(
                "validated_mirrored_ledger_state.json", state=state
            ),
            "captured_snapshots": {
                filename: {
                    "path": filename,
                    "sha256": hashlib.sha256((root / filename).read_bytes()).hexdigest(),
                }
                for filename in sorted(
                    {
                        "capture_same_account_process_snapshot.txt",
                        "capture_queue_snapshot.txt",
                        "capture_pending_submission_marker.txt",
                        "capture_pending_manual_accounting_marker.txt",
                        "capture_mirrored_ledger_line_counts.txt",
                        "capture_validated_mirrored_ledger_state.json",
                    }
                )
            },
            "operator_statement": OPERATOR_STATEMENTS[phase],
        }
        attestation_path = root / "attestation.json"
        attestation_path.write_text(
            json.dumps(attestation, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        for path in root.iterdir():
            path.chmod(0o400)
        root.chmod(0o500)
        return attestation_path

    def _publish_test_control_plane_successor(
        self, root: Path, *, schema_suffix: str = "\n"
    ) -> Path:
        staging = root / "control_plane" / "test-successor-staging"
        shutil.copytree(self.control_plane_dir, staging)
        for path in staging.iterdir():
            path.chmod(path.stat().st_mode | 0o200)
        schema = staging / "control_plane.schema.json"
        schema.write_text(
            schema.read_text(encoding="utf-8") + schema_suffix,
            encoding="utf-8",
        )
        records = [
            {"path": name, "sha256": sha256(staging / name)}
            for name in CONTROL_PLANE_FILES
        ]
        digest = inventory_digest(records)
        (staging / "inventory.json").write_text(
            json.dumps(
                {"schema_version": 1, "version": digest, "files": records},
                indent=2,
                sort_keys=True,
            ) + "\n",
            encoding="utf-8",
        )
        make_tree_read_only(
            staging,
            executable_names={
                path.name for path in staging.iterdir()
                if path.suffix in {".py", ".sh"}
            },
        )
        destination = staging.with_name(digest)
        staging.rename(destination)
        return destination

    def _revalidate_clean_candidate_with_test_source(
        self, candidate_manifest_path: Path, **kwargs: object
    ) -> dict[str, object]:
        if self.authorized_clean_candidate_source_root is None:
            raise AssertionError("Clean-candidate test source root was not injected")
        return revalidate_clean_candidate.revalidate_clean_candidate(
            candidate_manifest_path,
            **kwargs,
            authorized_source_root=self.authorized_clean_candidate_source_root,
        )

    def _promote_test_control_plane_successor(self, successor: Path) -> None:
        self._write_policy(
            admission_smoke_overrides=(
                {"status": "closed_after_pass"}
                if self.registered_science_slices
                else None
            ),
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )

        def invoke() -> None:
            promote(
                self.policy,
                **self._pre_policy_promotion_attestation_arguments(
                    control_plane_version=successor.name
                ),
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

        if self.authorized_clean_candidate_source_root is not None:
            with patch(
                "promote_active_policy.revalidate_clean_candidate",
                side_effect=self._revalidate_clean_candidate_with_test_source,
            ):
                invoke()
        else:
            invoke()

    def _create_test_terminal_recovery_handoff(
        self,
        successor: Path,
        *,
        authorize_purged_cancelled_zero_execution: bool = False,
    ) -> Path:
        return create_handoff(
            job_id="12345",
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            receipts_jsonl=self.receipts,
            mirror_jsonl=self.mirror,
            authorize_purged_cancelled_zero_execution=(
                authorize_purged_cancelled_zero_execution
            ),
            attest_reviewed_purged_reservation_job_binding=(
                authorize_purged_cancelled_zero_execution
            ),
            control_plane_dir=successor,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )

    def _write_timeout(self, *, expires_delta: timedelta = timedelta(hours=1)) -> None:
        now = datetime.now(timezone.utc)
        self._write(
            "timeout.json",
            json.dumps(
                {
                    "athena_walltime_seconds": 300,
                    "scheduler_walltime_seconds": 600,
                    "environment_profile_sha256": sha256(self.sources / "environment.sh"),
                    "measured_utc": self._utc(now - timedelta(minutes=1)),
                    "expires_utc": self._utc(now + expires_delta),
                }
            ),
        )

    def _write_manual_accounting_authorization(
        self,
        *,
        jobs: list[dict[str, str]] | None = None,
        authorization_id: str = "q016-login-host-direct-srun",
        accounting_scope: object = "manual_direct_srun_accounting_only",
        bind_policy: bool = True,
    ) -> Path:
        parent = self.pic_root / "policy" / "manual_accounting_authorizations"
        parent.mkdir(parents=True, exist_ok=True)
        path = parent / f"{authorization_id}.json"
        data = json.dumps(
            {
                "schema_version": 1,
                "authorization_id": authorization_id,
                "accounting_scope": accounting_scope,
                "scientific_evidence_eligible": False,
                "jobs": jobs
                or [
                    {"job_id": "4746332", "expected_qos": "normal"},
                    {"job_id": "4746335", "expected_qos": "normal"},
                ],
            }
        )
        path.write_text(data, encoding="utf-8")
        path.chmod(0o444)
        mirror_parent = (
            self.project_home_root / "policy" / "manual_accounting_authorizations"
        )
        mirror_parent.mkdir(parents=True, exist_ok=True)
        mirror_path = mirror_parent / path.name
        mirror_path.write_text(data, encoding="utf-8")
        mirror_path.chmod(0o444)
        if bind_policy:
            self._write_policy(
                manual_accounting_authorizations=[
                    self._manual_accounting_policy_binding(path)
                ]
            )
            self._promote_policy()
        return path

    def _manual_accounting_policy_binding(self, path: Path) -> dict[str, str]:
        return {
            "authorization_id": path.stem,
            "path": str(path),
            "project_home_path": str(
                self.project_home_root
                / "policy"
                / "manual_accounting_authorizations"
                / path.name
            ),
            "sha256": sha256(path),
        }

    def _manual_accounting_scheduler_output(
        self, command: list[str], *args: object, **kwargs: object
    ) -> str:
        del args, kwargs
        if command[0] == TRUSTED_SQUEUE:
            return ""
        self.assertEqual(command[0], TRUSTED_SACCT)
        self.assertIn("--allocations", command)
        self.assertIn("--clusters=frontier", command)
        return (
            "4746332|FAILED|5|1||ast207|batch|normal\n"
            "4746335|COMPLETED|7|1||ast207|batch|normal\n"
        )

    def _manual_accounting_arguments(self, authorization: Path) -> dict[str, Path]:
        return {
            "authorization": authorization,
            "ledger_jsonl": self.ledger,
            "ledger_csv": self.csv,
            "receipts_jsonl": self.receipts,
            "mirror_jsonl": self.mirror,
            "control_plane_dir": self.control_plane_dir,
            "authorized_pic_root": self.pic_root,
            "authorized_project_home_root": self.project_home_root,
        }

    def _append_registered_probe(
        self,
        *,
        final_event_type: str = "reservation_cancelled",
        job_id: str = "987654",
        include_operator_attestations: bool = False,
    ) -> dict[str, object]:
        reservation_id = str(uuid.uuid4())
        reservation = {
            "event_type": "reservation",
            "reservation_id": reservation_id,
            "submission_id": str(uuid.uuid4()),
            "control_plane_version": self.control_plane_version,
            "git_commit": "a" * 40,
            "campaign": "manual-accounting-fixture",
            "test_id": "registered-probe",
            "partition": "batch",
            "qos": "debug",
            "qos_selection_reason": "schema-valid manual-accounting fixture",
            "queue_snapshot_sha256": "b" * 64,
            "site_policy_checked_utc": "2026-05-31T00:00:00Z",
            "requested_nodes": 1,
            "requested_walltime": "00:01:00",
            "reserved_node_hours": 1.0 / 60.0,
            "artifact_dir": str(self.pic_root / "runs" / "registered-probe"),
            "state": "reserved",
            "reconciled": False,
        }
        if include_operator_attestations:
            authorization_id = "manual-accounting-registered-probe"
            reservation.update(
                {
                    "submission_scope": "registered_science",
                    "registered_science_authorization_id": authorization_id,
                    "clean_candidate_manifest_sha256": "c" * 64,
                    "manifest_path": str(
                        self.pic_root / "manifests/manual-accounting-fixture/manifest.json"
                    ),
                    "manifest_sha256": "d" * 64,
                    "job_script_sha256": "e" * 64,
                    "executable_sha256": "f" * 64,
                    "active_policy_sha256": "1" * 64,
                    "active_promotion_sha256": "2" * 64,
                }
            )
            for phase in ("pre_manifest", "pre_submit_wrapper"):
                attestation = self._sealed_operator_attestation(
                    authorization_id, phase
                )
                reservation[f"{phase}_attestation_path"] = str(attestation)
                reservation[f"{phase}_attestation_sha256"] = sha256(attestation)
        appended = append_primary_event(
            self.ledger,
            self.csv,
            self.receipts,
            self.mirror,
            reservation,
            mirror_transport="filesystem_copy",
        )
        if final_event_type == "reservation":
            return appended
        transition_record: dict[str, object] = {
            **ledger.transition_payload(appended),
            "event_type": (
                "job_id_attached"
                if final_event_type == "reconciliation"
                else final_event_type
            ),
            "state": (
                "cancelled"
                if final_event_type == "reservation_cancelled"
                else "submitted"
            ),
        }
        if final_event_type == "reservation_cancelled":
            transition_record["notes"] = "schema-valid manual-accounting fixture"
        elif final_event_type in {"job_id_attached", "reconciliation"}:
            transition_record["job_id"] = job_id
        else:
            raise ValueError(f"Unsupported registered probe event: {final_event_type}")
        appended = append_primary_event(
            self.ledger,
            self.csv,
            self.receipts,
            self.mirror,
            transition_record,
            mirror_transport="filesystem_copy",
        )
        if final_event_type != "reconciliation":
            return appended
        return append_primary_event(
            self.ledger,
            self.csv,
            self.receipts,
            self.mirror,
            {
                **ledger.transition_payload(appended),
                "event_type": "reconciliation",
                "job_id": job_id,
                "state": "COMPLETED",
                "reconciled": True,
                "scheduler_reported_allocated_nodes": 1,
                "billed_nodes": 1,
                "elapsed_seconds": 0,
                "consumed_node_hours": 0.0,
                "cumulative_consumed_node_hours": 0.0,
            },
            mirror_transport="filesystem_copy",
        )

    def _retry_manual_accounting_and_require_clean_markers(
        self, arguments: dict[str, Path]
    ) -> list[dict[str, object]]:
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            events = reconcile_manual_allocations(**arguments)
        local_marker, mirror_marker = incomplete_manual_accounting_marker_paths(
            self.ledger, self.mirror
        )
        self.assertFalse(local_marker.exists())
        assert mirror_marker is not None
        self.assertFalse(mirror_marker.exists())
        self.assertEqual(
            validate_primary_chain(self.ledger),
            ledger.validate_mirrored_state(self.ledger, self.receipts, self.mirror),
        )
        return events

    def _strand_manual_accounting_after_orion_append(
        self, arguments: dict[str, Path]
    ) -> None:
        real_append = ledger._append_jsonl
        interrupted = False

        def interrupt_before_mirror(
            path: Path, record: dict[str, object]
        ) -> None:
            nonlocal interrupted
            if (
                not interrupted
                and path == self.mirror
                and record.get("event_type") == "manual_allocation_reconciliation"
            ):
                interrupted = True
                raise RuntimeError("interrupted after Orion append")
            real_append(path, record)

        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ), patch("ledger._append_jsonl", side_effect=interrupt_before_mirror):
            with self.assertRaisesRegex(RuntimeError, "after Orion append"):
                reconcile_manual_allocations(**arguments)

    def _write_policy(
        self,
        *,
        science_submission_freeze: dict[str, object] | None = None,
        admission_smoke_overrides: dict[str, object] | None = None,
        registered_science_slices: list[dict[str, object]] | None = None,
        **storage_overrides: object,
    ) -> None:
        storage_preflight = self._write_storage_preflight_evidence()
        storage = {
            "installed_control_plane_version": self.control_plane_version,
            "staged_control_plane_candidate_version": self.control_plane_version,
            "installed_control_plane_lifecycle": "paired_installed_reviewed_generation",
            **storage_preflight,
            "project_home_mirror_root": str(self.project_home_root),
            "project_home_usage": [
                "small_append_only_ledger_and_control_plane_mirror",
            ],
            "project_home_retention_role": "operational_ledger_mirror_only",
            "project_home_ledger_mirror_transport": "filesystem_copy",
            "orion_bulk_evidence_root": str(self.pic_root),
            "orion_bulk_evidence_usage": [
                "simulation_outputs",
                "immutable_bulk_signoff_bundles",
                "private_reference_artifact_staging",
                "restore_drill_evidence",
            ],
            "orion_retention_role": (
                "user_selected_sole_bulk_evidence_root_with_documented_durability_risk"
            ),
            "manual_accounting_authorizations": [],
            "ledger_genesis_allowed": True,
        }
        if hasattr(self, "_closed_genesis"):
            storage.update(
                ledger_genesis_allowed=False,
                ledger_genesis=self._closed_genesis,
            )
        storage.update(storage_overrides)
        admission_smoke = {
            "status": "authorized_f0_parser_contract_only",
            "campaign": "f0_hipmpi_smoke",
            "test_id": "pic_parser_contract_guards",
            "evidence_class": "frontier_f0_admission_smoke_candidate",
            "physical_mode": "extended_mhd_pic_parser_contract",
            "selected_qos": "debug",
            "registered_short_nonproduction": True,
            "maximum_nodes": 1,
            "maximum_walltime_seconds": 15 * 60,
            "job_script_sha256": sha256(self.sources / "job.sh"),
            "input_deck_sha256": sha256(self.sources / "input.athinput"),
            "environment_profile_sha256": sha256(self.sources / "environment.sh"),
            "analysis_script_sha256": [sha256(self.sources / "analysis.py")],
            "executable_sha256": sha256(self.sources / "athena"),
            "launch_contract_sha256": launch_contract_sha256(
                self._launch_contract()
            ),
        }
        if (admission_smoke_overrides or {}).get("status") in {
            "closed_after_pass",
            "pending_exact_executable_binding",
        }:
            admission_smoke = {"status": admission_smoke_overrides["status"]}
        admission_smoke.update(admission_smoke_overrides or {})
        policy = {
            "schema_version": 1,
            "frontier": {
                "account": "AST207",
                "partition": "batch",
                "simulation_root": str(self.pic_root),
                "maximum_node_hours": 10000.0,
                "serial_pic_submissions": True,
            },
            "science_submission_freeze": (
                science_submission_freeze
                or self.science_submission_freeze
                or {"status": "pending_clean_candidate_freeze"}
            ),
            "registered_science_slices": (
                self.registered_science_slices
                if registered_science_slices is None
                else registered_science_slices
            ),
            "frontier_admission_smoke": admission_smoke,
            "olcf_side_storage": storage,
            "long_term_storage": {
                "status": "user_selected_orion_only_with_documented_durability_risk",
                "selected_destination": str(self.pic_root),
                "risk": (
                    "Orion-only retention is user-directed and does not provide an "
                    "institutional or approved off-site durable archive."
                ),
                "blocks": [
                    "terminal_durable_retention_signoff_pending_external_review",
                ],
            },
        }
        self.policy.write_text(json.dumps(policy), encoding="utf-8")

    def _write_storage_preflight_evidence(
        self,
        *,
        orion_root: Path | None = None,
        project_home_root: Path | None = None,
    ) -> dict[str, object]:
        orion_root = self.pic_root if orion_root is None else orion_root
        project_home_root = (
            self.project_home_root
            if project_home_root is None
            else project_home_root
        )
        probe_id = str(uuid.uuid4())
        completed_utc = "2026-06-03T00:00:00Z"
        relative = Path("policy/storage_preflight_evidence") / f"{probe_id}.json"
        orion_path = orion_root / relative
        project_home_path = project_home_root / relative
        artifact = {
            "completed_utc": completed_utc,
            "method": "local_create_write_sync_remove_probe",
            "probes": [
                {
                    "operations": AUTHORIZED_STORAGE_PREFLIGHT_OPERATIONS,
                    "path": str(root),
                    "payload_bytes": 32,
                    "payload_sha256": digest,
                    "role": role,
                    "st_dev": root.stat().st_dev,
                    "st_ino": root.stat().st_ino,
                    "status": "passed",
                }
                for role, root, digest in [
                    ("orion_simulation_root", orion_root, "1" * 64),
                    ("project_home_mirror_root", project_home_root, "2" * 64),
                ]
            ],
            "probe_id": probe_id,
            "publication": {
                "orion_path": str(orion_path),
                "project_home_path": str(project_home_path),
            },
            "record_type": "frontier_pic_storage_preflight_evidence",
            "schema_version": 2,
            "source_authentication": {
                **AUTHORIZED_STORAGE_PREFLIGHT_CAPTURE_SOURCE_BLOBS,
                "common_sha256": sha256(Path(control_plane_common.__file__)),
                "git_commit": "a" * 40,
                "tracked_clean_head_blobs": True,
            },
            "started_utc": completed_utc,
            "status": "passed",
        }
        payload = (
            json.dumps(
                artifact,
                allow_nan=False,
                ensure_ascii=True,
                separators=(",", ":"),
                sort_keys=True,
            )
            + "\n"
        ).encode("utf-8")
        for path in [orion_path, project_home_path]:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
            path.chmod(0o444)
        return {
            "last_preflight_utc": completed_utc,
            "orion_simulation_root_preflight": {
                "method": "local_create_write_sync_remove_probe",
                "path": str(orion_root),
                "status": "passed",
            },
            "project_home_preflight": {
                "method": "local_create_write_sync_remove_probe",
                "path": str(project_home_root),
                "status": "passed",
            },
            "storage_preflight_evidence": {
                "orion_path": str(orion_path),
                "probe_id": probe_id,
                "project_home_path": str(project_home_path),
                "sha256": hashlib.sha256(payload).hexdigest(),
            },
        }

    def _rewrite_reviewed_storage_preflight_artifact(
        self, mutate: Callable[[dict[str, object]], None]
    ) -> None:
        policy = json.loads(self.policy.read_text(encoding="utf-8"))
        binding = policy["olcf_side_storage"]["storage_preflight_evidence"]
        paths = [Path(str(binding[key])) for key in ["orion_path", "project_home_path"]]
        artifact = json.loads(paths[0].read_text(encoding="utf-8"))
        mutate(artifact)
        payload = (
            json.dumps(
                artifact,
                allow_nan=False,
                ensure_ascii=True,
                separators=(",", ":"),
                sort_keys=True,
            )
            + "\n"
        ).encode("utf-8")
        for path in paths:
            path.chmod(0o644)
            path.write_bytes(payload)
            path.chmod(0o444)
        binding["sha256"] = hashlib.sha256(payload).hexdigest()
        self.policy.write_text(json.dumps(policy), encoding="utf-8")

    def _rewrite_active_policy_as_historical_predecessor(self) -> None:
        policy_paths = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
        ]
        policy = json.loads(policy_paths[0].read_text(encoding="utf-8"))
        storage = policy["olcf_side_storage"]
        storage.pop("storage_preflight_evidence")
        storage["project_home_preflight"].pop("path")
        payload = json.dumps(policy).encode("utf-8")
        for path in policy_paths:
            path.chmod(0o644)
            path.write_bytes(payload)
            path.chmod(0o444)
        promotion_paths = [
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        promotion = json.loads(promotion_paths[0].read_text(encoding="utf-8"))
        promotion["policy_sha256"] = hashlib.sha256(payload).hexdigest()
        promotion_payload = (json.dumps(promotion, indent=2, sort_keys=True) + "\n").encode(
            "utf-8"
        )
        for path in promotion_paths:
            path.chmod(0o644)
            path.write_bytes(promotion_payload)
            path.chmod(0o444)

    def _rewrite_active_storage_preflight_artifact(
        self, mutate: Callable[[dict[str, object]], None]
    ) -> None:
        policy_paths = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
        ]
        policy = json.loads(policy_paths[0].read_text(encoding="utf-8"))
        binding = policy["olcf_side_storage"]["storage_preflight_evidence"]
        evidence_paths = [
            Path(str(binding[key])) for key in ["orion_path", "project_home_path"]
        ]
        artifact = json.loads(evidence_paths[0].read_text(encoding="utf-8"))
        mutate(artifact)
        policy["olcf_side_storage"]["last_preflight_utc"] = artifact["completed_utc"]
        evidence_payload = (
            json.dumps(
                artifact,
                allow_nan=False,
                ensure_ascii=True,
                separators=(",", ":"),
                sort_keys=True,
            )
            + "\n"
        ).encode("utf-8")
        for path in evidence_paths:
            path.chmod(0o644)
            path.write_bytes(evidence_payload)
            path.chmod(0o444)
        binding["sha256"] = hashlib.sha256(evidence_payload).hexdigest()
        policy_payload = json.dumps(policy).encode("utf-8")
        for path in policy_paths:
            path.chmod(0o644)
            path.write_bytes(policy_payload)
            path.chmod(0o444)
        promotion_paths = [
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        promotion = json.loads(promotion_paths[0].read_text(encoding="utf-8"))
        promotion["policy_sha256"] = hashlib.sha256(policy_payload).hexdigest()
        promotion_payload = (json.dumps(promotion, indent=2, sort_keys=True) + "\n").encode(
            "utf-8"
        )
        for path in promotion_paths:
            path.chmod(0o644)
            path.write_bytes(promotion_payload)
            path.chmod(0o444)

    def _install_deliberately_malformed_active_policy_fixture(self) -> None:
        """Fabricate malformed authorized state solely for downstream rejection tests."""
        policy_payload = self.policy.read_bytes()
        for path in [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
        ]:
            path.chmod(0o644)
            path.write_bytes(policy_payload)
            path.chmod(0o444)
        promotion_paths = [
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        promotion = json.loads(promotion_paths[0].read_text(encoding="utf-8"))
        promotion["policy_sha256"] = hashlib.sha256(policy_payload).hexdigest()
        promotion_payload = (json.dumps(promotion, indent=2, sort_keys=True) + "\n").encode(
            "utf-8"
        )
        for path in promotion_paths:
            path.chmod(0o644)
            path.write_bytes(promotion_payload)
            path.chmod(0o444)

    def _rewrite_active_policy_as_exact_reviewed_preflight_predecessor(self) -> None:
        def mutate(artifact: dict[str, object]) -> None:
            artifact.update(
                completed_utc="2026-06-02T00:00:00Z",
                schema_version=2,
                source_authentication={
                    **AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_SOURCE_AUTHENTICATION
                },
                started_utc="2026-06-02T00:00:00Z",
            )

        self._rewrite_active_storage_preflight_artifact(mutate)

    @contextmanager
    def _authorize_active_exact_reviewed_preflight_predecessor(
        self, *, control_plane_version: str | None = None
    ) -> Iterator[None]:
        policy_path = self.pic_root / "policy" / "storage_policy.json"
        promotion_path = self.pic_root / "policy" / "active_promotion.json"
        policy = json.loads(policy_path.read_text(encoding="utf-8"))
        binding = policy["olcf_side_storage"]["storage_preflight_evidence"]
        predecessor_control_plane_version = (
            control_plane_version or self.control_plane_version
        )
        with patch(
            "control_plane_common."
            "AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_POLICY_SHA256",
            sha256(policy_path),
        ), patch(
            "control_plane_common."
            "AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_PROMOTION_SHA256",
            sha256(promotion_path),
        ), patch(
            "control_plane_common."
            "AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_CONTROL_PLANE_VERSION",
            predecessor_control_plane_version,
        ), patch(
            "control_plane_common."
            "AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_PROBE_ID",
            binding["probe_id"],
        ), patch(
            "control_plane_common."
            "AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_EVIDENCE_SHA256",
            binding["sha256"],
        ):
            yield

    def _promote_policy(
        self, *, patch_clean_candidate_revalidation: bool = True
    ) -> None:
        def invoke() -> None:
            promote(
                self.policy,
                **self._pre_policy_promotion_attestation_arguments(),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

        if (
            patch_clean_candidate_revalidation
            and self.authorized_clean_candidate_source_root is not None
        ):
            with patch(
                "promote_active_policy.revalidate_clean_candidate",
                side_effect=self._revalidate_clean_candidate_with_test_source,
            ):
                invoke()
        else:
            invoke()

    def _pre_policy_promotion_attestation_arguments(
        self, *, control_plane_version: str | None = None
    ) -> dict[str, object]:
        policy = json.loads(self.policy.read_text(encoding="utf-8"))
        if not policy["registered_science_slices"]:
            return {}
        authorization_id = "reviewed-policy-promotion"
        return {
            "pre_policy_promotion_attestation": self._sealed_operator_attestation(
                authorization_id,
                "pre_policy_promotion",
                control_plane_version=control_plane_version,
            ),
            "pre_policy_promotion_authorization_id": authorization_id,
        }

    def _write_config(self, **overrides: object) -> None:
        now = datetime.now(timezone.utc)
        config = {
            "pic_root": str(self.pic_root),
            "campaign": "f0_hipmpi_smoke",
            "test_id": "pic_parser_contract_guards",
            "submission_id": self.submission_id,
            "job_script": str(self.sources / "job.sh"),
            "executable": str(self.sources / "athena"),
            "input_deck": str(self.sources / "input.athinput"),
            "environment_profile": str(self.sources / "environment.sh"),
            "timeout_margin_artifact": str(self.sources / "timeout.json"),
            "analysis_scripts": [str(self.sources / "analysis.py")],
            "queue_snapshot": str(self.sources / "queue.txt"),
            "submission_scope": "frontier_admission_smoke",
            "job_script_executable_env": "PIC_EXECUTABLE",
            "launch_contract": self._launch_contract(),
            "git_commit": "abc123",
            "evidence_class": "frontier_f0_admission_smoke_candidate",
            "physical_mode": "extended_mhd_pic_parser_contract",
            "selected_qos": "debug",
            "qos_selection_reason": "debug_available",
            "site_policy_checked_utc": self._utc(now),
            "registered_short_nonproduction": True,
            "artifact_dir": str(
                self.pic_root / "runs" / "f0_hipmpi_smoke" / self.submission_id
            ),
        }
        config.update(overrides)
        self.config.write_text(json.dumps(config), encoding="utf-8")

    def _launch_contract(self) -> dict[str, object]:
        return {
            "schema_version": 1,
            "executor": "trusted_trampoline_athena_argv_v1",
            "pre_actions": [],
            "actions": [
                {
                    "action_id": "athena-parser",
                    "kind": "athena",
                    "resources": {
                        "nodes": 1,
                        "tasks": 1,
                        "cpus_per_task": 1,
                        "gpus_per_task": 1,
                        "gpu_bind": "closest",
                    },
                    "arguments": [
                        {"literal": "-i"},
                        {"snapshot_role": "input-deck"},
                        {"literal": "-n"},
                    ],
                    "stdout_artifact": "athena_stdout.txt",
                    "stderr_artifact": "athena_stderr.txt",
                }
            ],
            "post_actions": [],
        }

    def _planner_clean_candidate(self) -> Path:
        import control_plane_common as common

        repo_root = Path(__file__).resolve().parents[3]
        source_root = self._clean_source("planner-candidate-source")
        self.authorized_clean_candidate_source_root = source_root
        inventory_relative = Path(self._prepared_artifact_inventory())
        inventory = json.loads((repo_root / inventory_relative).read_text(encoding="utf-8"))
        prepared_paths = {
            str(record["path"])
            for role in ("paper_decks", "analyzers")
            for record in inventory[role]
        }
        for stale in (
            "inputs/tests/pic_paper.athinput",
            "tst/publication/analyze_paper.py",
        ):
            if stale not in prepared_paths:
                (source_root / stale).unlink()
        relative_paths = {
            inventory_relative.as_posix(),
            *prepared_paths,
            *common.Q011_SECTION54_HELPER_SOURCES,
            *common.Q011_SECTION54_ARCHIVE_SOURCE_PATHS.values(),
        }
        for relative in relative_paths:
            destination = source_root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes((repo_root / relative).read_bytes())
        for role in ("paper_decks", "analyzers"):
            for record in inventory[role]:
                record["sha256"] = sha256(source_root / str(record["path"]))
        (source_root / inventory_relative).write_text(
            json.dumps(inventory, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        subprocess.run(["git", "-C", str(source_root), "add", "-A"], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "-c",
                "user.name=PIC Test",
                "-c",
                "user.email=pic-test@example.invalid",
                "commit",
                "-m",
                "add reviewed q011 planner sources",
            ],
            check=True,
            capture_output=True,
        )
        executable, profile = self._build_profile(
            source_root,
            self.pic_root / "planner-candidate-build",
            "hip-mpi-release-paper-pic",
        )
        return create_freeze(
            source_root=source_root,
            executable=executable,
            build_profile=profile,
            build_profile_id="hip-mpi-release-paper-pic",
            prepared_artifact_inventory=self._prepared_artifact_inventory(),
            freeze_id="b542ff53-0e5e-43d1-ae5b-988b9a94a92e",
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_source_root=source_root,
        )

    def _planner_retention(self) -> dict[str, object]:
        cached = getattr(self, "_planner_retention_binding", None)
        if cached is not None:
            return json.loads(json.dumps(cached))

        import control_plane_common as common

        repo_root = Path(__file__).resolve().parents[3]
        selected_case = {"case_id": "ps_p0_0p10", "problem_ps_p0": 0.1}
        published_pressure_receipt, published_packet_receipt = (
            self._planner_pressure_publication()
        )
        aggregate_receipt = json.loads(
            Path(published_pressure_receipt["path"]).read_text(encoding="utf-8")
        )
        aggregate_cases = aggregate_receipt["raw_cases"]
        environment_payload = (
            self.control_plane_dir / "frontier_pic_environment.sh"
        ).read_bytes()
        source_binding_paths = {
            name: path for name, path in common.Q011_SECTION54_SOURCE_BINDING_PATHS.items()
        }
        clean_candidate_path = self._planner_clean_candidate()
        clean_candidate_payload = clean_candidate_path.read_bytes()
        clean_candidate_manifest = json.loads(clean_candidate_payload)
        source = clean_candidate_manifest["source"]
        build = clean_candidate_manifest["build"]
        prepared = clean_candidate_manifest["prepared_artifacts"]
        freeze_id = str(clean_candidate_manifest["freeze_id"])
        candidate_root = clean_candidate_path.parent
        source_archive_path = Path(str(source["archive_path"]))
        source_archive_sha256 = hashlib.sha256(source_archive_path.read_bytes()).hexdigest()
        archive_members = {
            relative: b""
            for relative in {
                *common.Q011_SECTION54_HELPER_SOURCES,
                *common.Q011_SECTION54_ARCHIVE_SOURCE_PATHS.values(),
            }
        }
        with tarfile.open(fileobj=io.BytesIO(source_archive_path.read_bytes())) as archive:
            for relative in archive_members:
                member = archive.extractfile(relative)
                assert member is not None
                archive_members[relative] = member.read()
        helper_sources = [
            {
                "path": relative,
                "sha256": hashlib.sha256(archive_members[relative]).hexdigest(),
            }
            for relative in common.Q011_SECTION54_HELPER_SOURCES
        ]
        candidate_binding = {
            "clean_candidate_manifest": {
                "path": str(clean_candidate_path),
                "sha256": hashlib.sha256(clean_candidate_payload).hexdigest(),
            },
            "freeze_id": freeze_id,
            "git_commit": source["git_commit"],
            "git_tree": source["git_tree"],
            "source_archive_sha256": source_archive_sha256,
            "source_commit_sha256": source["commit_sha256"],
            "source_bundle_sha256": source["source_bundle_sha256"],
            "prepared_artifact_inventory_sha256": prepared["inventory_sha256"],
            "validated_submodules": [],
            "build_profile": {
                "path": build["profile_path"],
                "sha256": build["profile_sha256"],
            },
            "build_profile_receipt": {
                "path": build["profile_receipt_path"],
                "sha256": build["profile_receipt_sha256"],
            },
            "build_invocations_sha256": build["build_invocations_sha256"],
            "executable": {
                "path": build["executable_path"],
                "sha256": build["executable_sha256"],
            },
            "environment_profile": {
                "path": str(self.control_plane_dir / "frontier_pic_environment.sh"),
                "sha256": hashlib.sha256(environment_payload).hexdigest(),
                "control_plane_version": self.control_plane_version,
                "reviewed_source": {
                    "path": common.Q011_SECTION54_ARCHIVE_SOURCE_PATHS[
                        "environment_profile"
                    ],
                    "sha256": hashlib.sha256(environment_payload).hexdigest(),
                },
            },
        }
        reanalysis_attestation, reviewer_attestation = (
            self._planner_pressure_gate_attestations(
                published_pressure_receipt=published_pressure_receipt,
                published_packet_receipt=published_packet_receipt,
                pilot_bundle_manifest_sha256=aggregate_receipt["aggregate_bundle"][
                    "manifest_sha256"
                ],
                aggregate_pilot_analysis_sha256=aggregate_receipt[
                    "aggregate_analysis"
                ]["sha256"],
                selected_case=selected_case,
                git_commit=str(source["git_commit"]),
                source_archive_sha256=source_archive_sha256,
                helper_sources=helper_sources,
            )
        )
        pressure_receipt = {
            "schema_version": 3,
            "record_type": "q011_section54_pressure_selection_receipt",
            "selection_method": "human_review_only",
            "published_pressure_pilot_receipt": published_pressure_receipt,
            "published_pressure_pilot_review_packet_receipt": published_packet_receipt,
            "pilot_bundle_manifest_sha256": aggregate_receipt["aggregate_bundle"][
                "manifest_sha256"
            ],
            "aggregate_pilot_analysis_sha256": aggregate_receipt[
                "aggregate_analysis"
            ]["sha256"],
            "case_descriptors": [
                {
                    "case_id": aggregate_case["case_id"],
                    "problem_ps_p0": problem_ps_p0,
                    "descriptor_sha256": aggregate_case["descriptor_sha256"],
                }
                for aggregate_case, (_, problem_ps_p0) in zip(
                    aggregate_cases,
                    (
                        ("ps_p0_1p00", 1.0),
                        ("ps_p0_0p05", 0.05),
                        ("ps_p0_0p10", 0.1),
                        ("ps_p0_0p20", 0.2),
                    ),
                )
            ],
            "selected_case": selected_case,
            "authoritative_reanalysis_attestation": reanalysis_attestation,
            "reviewer_attestation": reviewer_attestation,
        }
        source_payloads = {
            "bindings/human_pressure_selection_receipt.json": (
                json.dumps(pressure_receipt, indent=2, sort_keys=True) + "\n"
            ).encode("utf-8"),
            "bindings/clean_candidate_manifest.json": clean_candidate_payload,
            "bindings/environment_profile.sh": environment_payload,
            "bindings/q011_section54_qualifying_campaign_preregistration.json": (
                repo_root
                / common.Q011_SECTION54_ARCHIVE_SOURCE_PATHS[
                    "qualifying_preregistration"
                ]
            ).read_bytes(),
            "bindings/q011_section54_restart_continuation_preregistration.json": (
                repo_root
                / common.Q011_SECTION54_ARCHIVE_SOURCE_PATHS["restart_preregistration"]
            ).read_bytes(),
            "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput": (
                repo_root / common.Q011_SECTION54_ARCHIVE_SOURCE_PATHS["paper_deck"]
            ).read_bytes(),
        }
        source_bindings = {
            name: {"path": path, "sha256": hashlib.sha256(source_payloads[path]).hexdigest()}
            for name, path in source_binding_paths.items()
        }
        campaign_matrix = common._planner_expected_matrix()
        basis = {
            "record_type": "q011_section54_qualifying_campaign_execution_plan",
            "schema_version": 1,
            "pressure_selection_receipt_sha256": source_bindings[
                "pressure_selection_receipt"
            ]["sha256"],
            "selected_case": selected_case,
            "candidate_binding": candidate_binding,
            "source_binding_sha256": {
                name: binding["sha256"] for name, binding in source_bindings.items()
            },
            "helper_source_closure": helper_sources,
            "campaign_matrix": campaign_matrix,
            "authorized_orion_root": str(self.pic_root),
        }
        plan_id = hashlib.sha256(
            json.dumps(
                basis, sort_keys=True, separators=(",", ":"), allow_nan=False
            ).encode("utf-8")
        ).hexdigest()
        planner_root = (
            self.pic_root / "plans" / f"q011-section54-qualifying-campaign-plan-{plan_id}"
        )
        campaign_root = self.pic_root / "campaigns" / f"q011-section54-{plan_id}"
        files = dict(source_payloads)
        helper_payload = (
            json.dumps(
                {
                    "record_type": "q011_section54_helper_source_closure",
                    "schema_version": 1,
                    "plan_id": plan_id,
                    "sources": helper_sources,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        ).encode("utf-8")
        files["helper_source_closure.json"] = helper_payload
        helper_binding = {
            "path": "helper_source_closure.json",
            "sha256": hashlib.sha256(helper_payload).hexdigest(),
        }
        descriptors = []
        selected_contract = None
        contract_bindings = []
        restart_source_attempt = None
        selected_attempt_id = common._planner_attempt_id(
            1, common.Q011_SECTION54_VARIANTS[0][0], common.Q011_SECTION54_SEEDS[0]
        )
        index = 0
        for variant, model_overrides in common.Q011_SECTION54_VARIANTS:
            for seed in common.Q011_SECTION54_SEEDS:
                index += 1
                attempt_id = common._planner_attempt_id(index, variant, seed)
                attempt_root = campaign_root / "baseline" / attempt_id
                contract_path = f"launch_contracts/baseline/{attempt_id}.json"
                contract = common._planner_expected_baseline_contract(
                    attempt_id=attempt_id,
                    variant=variant,
                    model_overrides=model_overrides,
                    seed=seed,
                    selected_ps_p0=selected_case["problem_ps_p0"],
                    candidate=candidate_binding,
                    paper_deck_binding=source_bindings["paper_deck"],
                    attempt_root=attempt_root,
                )
                contract_payload = (
                    json.dumps(contract, indent=2, sort_keys=True) + "\n"
                ).encode("utf-8")
                files[contract_path] = contract_payload
                contract_binding = {
                    "path": contract_path,
                    "sha256": hashlib.sha256(contract_payload).hexdigest(),
                }
                contract_bindings.append(contract_binding)
                descriptor_path = f"attempts/baseline/{attempt_id}.json"
                descriptor = common._planner_expected_baseline_descriptor(
                    index=index,
                    attempt_id=attempt_id,
                    variant=variant,
                    seed=seed,
                    selected_ps_p0=selected_case["problem_ps_p0"],
                    candidate=candidate_binding,
                    attempt_root=attempt_root,
                    contract_path=contract_path,
                    contract_payload=contract_payload,
                )
                descriptor_payload = (
                    json.dumps(descriptor, indent=2, sort_keys=True) + "\n"
                ).encode("utf-8")
                files[descriptor_path] = descriptor_payload
                descriptors.append(
                    {
                        "path": descriptor_path,
                        "sha256": hashlib.sha256(descriptor_payload).hexdigest(),
                    }
                )
                if index == 1:
                    selected_contract = contract
                if (
                    variant == "three_level_amr_root_dx12_finest_dx3"
                    and seed == common.Q011_SECTION54_SEEDS[0]
                ):
                    restart_source_attempt = descriptor
        assert restart_source_attempt is not None
        restart_preregistration = json.loads(
            source_payloads[
                "bindings/q011_section54_restart_continuation_preregistration.json"
            ]
        )
        carrier_id = common._planner_restart_carrier_id(
            restart_source_attempt["qualifying_seed"]
        )
        restart_root = campaign_root / "restart_continuation" / carrier_id
        restart_contract_path = f"launch_contracts/restart_continuation/{carrier_id}.json"
        restart_contract = common._planner_expected_restart_contract(
            carrier_id=carrier_id,
            source_attempt=restart_source_attempt,
            restart_preregistration=restart_preregistration,
            candidate=candidate_binding,
            paper_deck_binding=source_bindings["paper_deck"],
            attempt_root=restart_root,
        )
        restart_contract_payload = (
            json.dumps(restart_contract, indent=2, sort_keys=True) + "\n"
        ).encode("utf-8")
        files[restart_contract_path] = restart_contract_payload
        restart_contract_binding = {
            "path": restart_contract_path,
            "sha256": hashlib.sha256(restart_contract_payload).hexdigest(),
        }
        restart_carrier = common._planner_expected_restart_carrier(
            carrier_id=carrier_id,
            source_attempt=restart_source_attempt,
            restart_preregistration=restart_preregistration,
            restart_preregistration_binding=source_bindings["restart_preregistration"],
            attempt_root=restart_root,
            contract_path=restart_contract_path,
            contract_payload=restart_contract_payload,
        )
        restart_carrier_payload = (
            json.dumps(restart_carrier, indent=2, sort_keys=True) + "\n"
        ).encode("utf-8")
        restart_carrier_path = (
            "restart_continuation/amr_restart_continuation_carrier.json"
        )
        files[restart_carrier_path] = restart_carrier_payload
        restart_carrier_binding = {
            "path": restart_carrier_path,
            "sha256": hashlib.sha256(restart_carrier_payload).hexdigest(),
        }
        recompute = common._planner_expected_independent_recompute_plan(
            plan_id=plan_id,
            campaign_root=campaign_root,
            qualifying_preregistration_binding=source_bindings[
                "qualifying_preregistration"
            ],
        )
        recompute_payload = (
            json.dumps(recompute, indent=2, sort_keys=True) + "\n"
        ).encode("utf-8")
        recompute_path = "independent_raw_artifact_recompute_plan.json"
        files[recompute_path] = recompute_payload
        recompute_binding = {
            "path": recompute_path,
            "sha256": hashlib.sha256(recompute_payload).hexdigest(),
        }
        fragment = common._planner_expected_policy_fragment(
            plan_id=plan_id,
            pic_root=self.pic_root,
            campaign_root=campaign_root,
            candidate=candidate_binding,
            pressure_receipt_binding=source_bindings["pressure_selection_receipt"],
            contract_bindings=contract_bindings,
            restart_contract_binding=restart_contract_binding,
        )
        fragment_payload = (
            json.dumps(fragment, indent=2, sort_keys=True) + "\n"
        ).encode("utf-8")
        fragment_path = "nonauthorizing_policy_fragment.json"
        files[fragment_path] = fragment_payload
        fragment_binding = {
            "path": fragment_path,
            "sha256": hashlib.sha256(fragment_payload).hexdigest(),
        }
        qualifying = json.loads(
            source_payloads[
                "bindings/q011_section54_qualifying_campaign_preregistration.json"
            ]
        )
        plan_payload = (
            json.dumps(
                {
                    "record_type": "q011_section54_qualifying_campaign_execution_plan",
                    "schema_version": 1,
                    "plan_id": plan_id,
                    "artifact_role": (
                        "q011_section54_source_local_immutable_qualifying_campaign_plan"
                    ),
                    "qualification_effect": (
                        "plan_only_no_execution_authorization_no_claim_closure"
                    ),
                    "status": "source_local_immutable_review_plan_only",
                    "authorized_orion_root": str(self.pic_root),
                    "authorized_orion_campaign_root": str(campaign_root),
                    "selected_pressure": {
                        "selection_method": "human_review_only",
                        "selected_case": selected_case,
                        "receipt": source_bindings["pressure_selection_receipt"],
                    },
                    "candidate_binding": candidate_binding,
                    "source_bindings": source_bindings,
                    "helper_source_closure": helper_binding,
                    "campaign_matrix": campaign_matrix,
                    "baseline_attempt_count": 24,
                    "baseline_attempt_descriptors": descriptors,
                    "restart_continuation_carrier": restart_carrier_binding,
                    "independent_raw_artifact_recompute_plan": recompute_binding,
                    "nonauthorizing_policy_fragment": fragment_binding,
                    "execution_boundary": {
                        "mutates_live_policy": False,
                        "scheduler_calls": False,
                        "submits_jobs": False,
                        "infers_pressure_selection": False,
                        "launch_authorized": False,
                        "frontier_execution_authorized": False,
                        "claim_closure_authorized": False,
                    },
                    "preregistration_execution_boundary": qualifying[
                        "qualifying_execution_bindings"
                    ],
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        ).encode("utf-8")
        files["campaign_plan.json"] = plan_payload
        materialized_inventory_payload = "".join(
            f"{hashlib.sha256(payload).hexdigest()}  {path}\n"
            for path, payload in sorted(files.items())
        ).encode("utf-8")
        receipt_payload = (
            json.dumps(
                {
                    "record_type": (
                        "q011_section54_qualifying_campaign_plan_materialization_receipt"
                    ),
                    "schema_version": 1,
                    "plan_id": plan_id,
                    "campaign_plan": {
                        "path": "campaign_plan.json",
                        "sha256": hashlib.sha256(plan_payload).hexdigest(),
                    },
                    "helper_source_closure": helper_binding,
                    "tree_inventory": {
                        "algorithm": (
                            "sha256 of '<file_sha256>  <root-relative-path>\\n' "
                            "entries ordered lexically by root-relative path"
                        ),
                        "scope": (
                            "all materialized campaign-plan members before this "
                            "receipt and recursive-freeze metadata"
                        ),
                        "excludes": [
                            "materialization_receipt.json",
                            "freeze_receipt.json",
                            "artifact_inventory.sha256",
                        ],
                        "sha256": hashlib.sha256(
                            materialized_inventory_payload
                        ).hexdigest(),
                        "inventoried_file_count": len(files),
                    },
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        ).encode("utf-8")
        files["materialization_receipt.json"] = receipt_payload
        files["freeze_receipt.json"] = (
            json.dumps(
                {
                    "schema_version": 1,
                    "artifact_role": (
                        "q011_section54_source_local_immutable_qualifying_campaign_plan"
                    ),
                    "qualification_effect": (
                        "plan_only_no_execution_authorization_no_claim_closure"
                    ),
                    "inventory_excludes": "artifact_inventory.sha256",
                    "freeze_policy": (
                        "remove all owner, group and other write bits recursively"
                    ),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        ).encode("utf-8")
        inventory_payload = "".join(
            f"{hashlib.sha256(payload).hexdigest()}  {path}\n"
            for path, payload in sorted(files.items())
        ).encode("utf-8")
        planner_root.mkdir(parents=True)
        for relative, payload in files.items():
            path = planner_root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
        (planner_root / "artifact_inventory.sha256").write_bytes(inventory_payload)
        for path in sorted(planner_root.rglob("*"), reverse=True):
            path.chmod(0o555 if path.is_dir() else 0o444)
        planner_root.chmod(0o555)
        assert selected_contract is not None
        attempt_root = campaign_root / "baseline" / selected_attempt_id
        self._planner_retention_binding = {
            "schema_version": 1,
            "retention_role": "q011_section54_deterministic_retained_attempt",
            "planner_root": str(planner_root),
            "planner_inventory_sha256": hashlib.sha256(inventory_payload).hexdigest(),
            "planner_plan_id": plan_id,
            "planner_materialization_receipt": {
                "path": "materialization_receipt.json",
                "sha256": hashlib.sha256(receipt_payload).hexdigest(),
            },
            "attempt_id": selected_attempt_id,
            "authorized_orion_attempt_root": str(attempt_root),
            "authorized_orion_raw_root": str(attempt_root / "raw"),
            "argv": selected_contract["argv"],
        }
        self._planner_candidate_manifest = clean_candidate_path
        return json.loads(json.dumps(self._planner_retention_binding))

    def _planner_pressure_publication(
        self,
    ) -> tuple[dict[str, str], dict[str, str]]:
        publication = self.pic_root / "publication"
        acceptance = self.pic_root / "publication_acceptance"
        runs = self.pic_root / "runs"
        publication.mkdir(exist_ok=True)
        acceptance.mkdir(exist_ok=True)
        runs.mkdir(exist_ok=True)
        publication_identity = {
            "device": publication.stat().st_dev,
            "inode": publication.stat().st_ino,
        }
        source_bindings = {
            "postrun_aggregate_source_authorization": {
                "path": "tst/publication/readiness/source-authorization.json",
                "sha256": "1" * 64,
            },
            "registered_execution_preregistration": {
                "path": "tst/publication/readiness/registered-execution.json",
                "sha256": "2" * 64,
            },
            "historical_v2_execution_preregistration": {
                "path": "tst/publication/readiness/historical-v2.json",
                "sha256": "3" * 64,
            },
            "reviewed_source_closure": [
                {
                    "role": "aggregate_publisher",
                    "path": "tst/publication/publish_q011_section54_pressure_pilot_bundle.py",
                    "sha256": "4" * 64,
                },
                {
                    "role": "review_packet_renderer",
                    "path": "tst/publication/render_q011_section54_pressure_pilot_review_packet.py",
                    "sha256": "5" * 64,
                },
            ],
            "runtime_source_archive": {
                "execution_mode": "direct_api_nonproduction_only",
                "git_commit": None,
                "archive_sha256": None,
                "verified_source_closure_sha256": None,
            },
        }
        aggregate_bundle = publication / "pressure-pilot-bundle"
        aggregate_bundle.mkdir()
        raw_cases = []
        manifest_cases = []
        for case_id, pressure, argv_value in zip(
            pressure_packet_verifier.RAW_CASE_IDS,
            pressure_packet_verifier.RAW_CASE_PRESSURES,
            pressure_packet_verifier.RAW_CASE_ARGV_VALUES,
        ):
            artifact_dir = runs / f"raw-{case_id}" / "attempt"
            artifact_dir.mkdir(parents=True)
            raw_payloads: dict[str, bytes] = {}
            aggregate_payloads: dict[str, bytes] = {}
            raw_inventory: dict[str, dict[str, object]] = {}
            bundle_members: list[dict[str, object]] = []

            def bind_raw(source: str, payload: bytes) -> None:
                existing = raw_payloads.get(source)
                if existing is not None and existing != payload:
                    raise AssertionError(f"fixture payload collision: {source}")
                raw_payloads[source] = payload
                raw_inventory[source] = {
                    "path": source,
                    "sha256": hashlib.sha256(payload).hexdigest(),
                    "size": len(payload),
                }

            def bind_bundle(source: str, target: str, payload: bytes) -> dict[str, str]:
                bind_raw(source, payload)
                aggregate_payloads[target] = payload
                bundle_members.append(
                    {
                        "path": target,
                        "sha256": hashlib.sha256(payload).hexdigest(),
                        "size": len(payload),
                        "source_path": source,
                    }
                )
                return {"path": target, "sha256": hashlib.sha256(payload).hexdigest()}

            snapshots = []
            for snapshot_index, time in enumerate(
                pressure_packet_verifier.RAW_CASE_TIMES
            ):
                suffix = f"{snapshot_index:05d}"
                snapshot: dict[str, object] = {"time": time}
                for kind, directory, extension in (
                    ("mhd_w_bcc", "bin", "bin"),
                    ("bmag", "bin", "bin"),
                    ("prtcl_jx", "bin", "bin"),
                    ("j2", "bin", "bin"),
                    ("prtcl_all", "pvtk", "part.vtk"),
                ):
                    target = (
                        f"cases/{case_id}/{directory}/"
                        f"{case_id}.{kind}.{suffix}.{extension}"
                    )
                    source = (
                        f"output/{directory}/{case_id}.{kind}.{suffix}.{extension}"
                    )
                    snapshot[kind] = bind_bundle(
                        source,
                        target,
                        f"{case_id}:{kind}:{suffix}\n".encode("ascii"),
                    )
                snapshots.append(snapshot)

            stdout_payload = f"{case_id}: stdout\n".encode("ascii")
            stdout = bind_bundle(
                "athena_stdout.txt",
                f"cases/{case_id}/stdout.txt",
                stdout_payload,
            )
            restart_prefix = f"{case_id}.00004.rst"
            restart_bindings = {}
            for name, suffix in (
                ("manifest", ".manifest"),
                ("manifest_complete", ".manifest.complete"),
                ("artifact", ""),
                ("complete", ".complete"),
            ):
                filename = restart_prefix + suffix
                restart_bindings[name] = bind_bundle(
                    f"output/rst/{filename}",
                    f"cases/{case_id}/rst/{filename}",
                    f"{case_id}:{name}\n".encode("ascii"),
                )
            manifest_case = {
                "case_id": case_id,
                "ps_p0": pressure,
                "overrides": [
                    *pressure_packet_verifier.AUTHORIZED_COMMON_OVERRIDES,
                    f"problem/ps_p0={argv_value}",
                ],
                "snapshots": snapshots,
                "stdout": stdout,
                "terminal_restart": {
                    "time": pressure_packet_verifier.RAW_CASE_TIMES[-1],
                    "manifest": restart_bindings["manifest"],
                    "manifest_complete": restart_bindings["manifest_complete"],
                    "members": [
                        {
                            "artifact": restart_bindings["artifact"],
                            "complete": restart_bindings["complete"],
                        }
                    ],
                },
            }
            manifest_cases.append(manifest_case)

            allowlist = (
                f"q011-pressure-{case_id.replace('_', '-')}.environment.allowlist.txt"
            )
            runtime_payloads = {
                allowlist: b"PIC_FRONTIER_PROFILE=frontier_minimum_supported\n",
                "athena_stdout.txt": stdout_payload,
                "athena_stdout.sha256": (
                    hashlib.sha256(stdout_payload).hexdigest() + "\n"
                ).encode("ascii"),
                "athena_stderr.txt": b"frontier diagnostic stderr\n",
            }
            for path, payload in runtime_payloads.items():
                bind_raw(path, payload)
            runtime_artifacts = {
                path: hashlib.sha256(payload).hexdigest()
                for path, payload in runtime_payloads.items()
            }
            inventory_payload = self._planner_json_payload(
                {
                    "schema_version": 1,
                    "files": [raw_inventory[path] for path in sorted(raw_inventory)],
                }
            )
            self._write_planner_readonly(
                artifact_dir / pressure_packet_verifier.RAW_CASE_INVENTORY_NAME,
                inventory_payload,
            )
            descriptor = {
                "schema_version": 1,
                "record_type": pressure_packet_verifier.RAW_CASE_RECORD_TYPE,
                "evidence_class": pressure_packet_verifier.AGGREGATE_EVIDENCE_CLASS,
                "qualification_effect": (
                    pressure_packet_verifier.AGGREGATE_QUALIFICATION_EFFECT
                ),
                "launch_contract": "trusted_trampoline_athena_argv_v1",
                "case_id": case_id,
                "ps_p0": pressure,
                "argv_value": argv_value,
                "artifact_inventory_sha256": hashlib.sha256(
                    inventory_payload
                ).hexdigest(),
                "runtime_artifacts": runtime_artifacts,
                "runtime_profile": pressure_packet_verifier.AUTHORIZED_RUNTIME_PROFILE,
                "parallel_ranks": pressure_packet_verifier.AUTHORIZED_PARALLEL_RANKS,
                "rank_gpu_bindings": [
                    {
                        "host": "frontier00000",
                        "rank": 0,
                        "rocr_visible_device": (
                            pressure_packet_verifier.AUTHORIZED_ROCR_VISIBLE_DEVICE
                        ),
                    }
                ],
                "manifest_case": manifest_case,
                "bundle_members": sorted(
                    bundle_members, key=lambda member: member["path"]
                ),
            }
            descriptor_path = (
                artifact_dir / pressure_packet_verifier.RAW_CASE_DESCRIPTOR_PATH
            )
            descriptor_path.parent.mkdir()
            descriptor_payload = self._planner_json_payload(descriptor)
            self._write_planner_readonly(descriptor_path, descriptor_payload)
            for path, payload in raw_payloads.items():
                destination = artifact_dir / path
                destination.parent.mkdir(parents=True, exist_ok=True)
                self._write_planner_readonly(destination, payload)
            for path, payload in aggregate_payloads.items():
                destination = aggregate_bundle / path
                destination.parent.mkdir(parents=True, exist_ok=True)
                self._write_planner_readonly(destination, payload)
            for path in sorted(
                artifact_dir.rglob("*"),
                key=lambda candidate: len(candidate.parts),
                reverse=True,
            ):
                path.chmod(0o555 if path.is_dir() else 0o444)
            artifact_dir.chmod(0o555)
            (artifact_dir / "analysis").chmod(0o700)
            raw_cases.append(
                {
                    "case_id": case_id,
                    "artifact_dir": str(artifact_dir),
                    "descriptor_path": pressure_packet_verifier.RAW_CASE_DESCRIPTOR_PATH,
                    "descriptor_sha256": hashlib.sha256(descriptor_payload).hexdigest(),
                    "artifact_inventory_sha256": hashlib.sha256(
                        inventory_payload
                    ).hexdigest(),
                    "runtime_artifacts": runtime_artifacts,
                }
            )

        aggregate_manifest_payload = self._planner_json_payload(
            {
                "schema_version": 1,
                "record_type": pressure_packet_verifier.AGGREGATE_MANIFEST_RECORD_TYPE,
                "evidence_class": pressure_packet_verifier.AGGREGATE_EVIDENCE_CLASS,
                "qualification_effect": (
                    pressure_packet_verifier.AGGREGATE_QUALIFICATION_EFFECT
                ),
                "active_deck_binding": dict(
                    pressure_packet_verifier.AUTHORIZED_ACTIVE_DECK_BINDING
                ),
                "preregistration_binding": source_bindings[
                    "postrun_aggregate_source_authorization"
                ],
                "registered_execution_preregistration_binding": source_bindings[
                    "registered_execution_preregistration"
                ],
                "cases": manifest_cases,
            }
        )
        self._write_planner_readonly(
            aggregate_bundle / pressure_packet_verifier.AGGREGATE_MANIFEST_NAME,
            aggregate_manifest_payload,
        )
        for path in sorted(
            aggregate_bundle.rglob("*"),
            key=lambda candidate: len(candidate.parts),
            reverse=True,
        ):
            path.chmod(0o555 if path.is_dir() else 0o444)
        aggregate_bundle.chmod(0o555)
        aggregate_analysis = publication / "pressure-pilot-analysis.json"
        self._write_planner_readonly(
            aggregate_analysis,
            self._planner_json_payload({"status": "pass"}),
        )
        aggregate_receipt = publication / "pressure-pilot-receipt.json"
        aggregate_payload = self._planner_json_payload(
            {
                "schema_version": 1,
                "record_type": pressure_packet_verifier.AGGREGATE_RECEIPT_RECORD_TYPE,
                "evidence_class": pressure_packet_verifier.AGGREGATE_EVIDENCE_CLASS,
                "qualification_effect": (
                    pressure_packet_verifier.AGGREGATE_QUALIFICATION_EFFECT
                ),
                "consumption_rule": pressure_packet_verifier.CONSUMPTION_RULE,
                "publication_root_identity": publication_identity,
                "aggregate_bundle": {
                    "path": str(aggregate_bundle),
                    "manifest_sha256": hashlib.sha256(
                        aggregate_manifest_payload
                    ).hexdigest(),
                },
                "aggregate_analysis": {
                    "path": str(aggregate_analysis),
                    "sha256": hashlib.sha256(aggregate_analysis.read_bytes()).hexdigest(),
                },
                "source_bindings": source_bindings,
                "raw_cases": raw_cases,
            }
        )
        self._write_planner_readonly(aggregate_receipt, aggregate_payload)
        aggregate_binding = {
            "path": str(aggregate_receipt),
            "sha256": hashlib.sha256(aggregate_payload).hexdigest(),
        }
        packet_root = publication / "pressure-review-packet"
        figures = packet_root / "figures"
        figures.mkdir(parents=True)
        packet_payloads = {
            "PRESSURE_REVIEW_PACKET.md": b"test-only pressure review packet\n",
            "figures/terminal_mhd_pic_pressure_comparison.png": b"pressure comparison\n",
            "figures/terminal_profile_overlays.png": b"profile overlays\n",
            "pressure_review_metrics.json": self._planner_json_payload(
                {
                    "schema_version": 1,
                    "record_type": pressure_packet_verifier.REVIEW_METRICS_RECORD_TYPE,
                    "watermark": pressure_packet_verifier.WATERMARK,
                    "qualification_effect": pressure_packet_verifier.QUALIFICATION_EFFECT,
                    "aggregate_receipt": aggregate_binding,
                    "cases": [
                        {
                            "case_id": case_id,
                            "problem_ps_p0": float(problem_ps_p0),
                            "terminal_particle_count": 1,
                            "particle_efficiency": 1.0,
                            "zone_cycles_per_second": 1.0,
                            "particle_updates_per_second": 1.0,
                            "tracked_gpu_memory_high_water_bytes": 1.0,
                        }
                        for case_id, problem_ps_p0 in zip(
                            pressure_packet_verifier.RAW_CASE_IDS,
                            pressure_packet_verifier.RAW_CASE_PRESSURES,
                        )
                    ],
                }
            ),
        }
        for relative, payload in packet_payloads.items():
            self._write_planner_readonly(packet_root / relative, payload)
        inventory_payload = self._planner_json_payload(
            {
                "schema_version": 1,
                "record_type": pressure_packet_verifier.INVENTORY_RECORD_TYPE,
                "members": [
                    {
                        "path": relative,
                        "sha256": hashlib.sha256(packet_payloads[relative]).hexdigest(),
                        "size": len(packet_payloads[relative]),
                    }
                    for relative in sorted(packet_payloads)
                ],
            }
        )
        self._write_planner_readonly(
            packet_root / pressure_packet_verifier.INVENTORY_NAME,
            inventory_payload,
        )
        figures.chmod(0o555)
        packet_root.chmod(0o555)
        packet_receipt = publication / "pressure-review-packet-receipt.json"
        packet_receipt_payload = self._planner_json_payload(
            {
                "schema_version": 1,
                "record_type": pressure_packet_verifier.PACKET_RECEIPT_RECORD_TYPE,
                "watermark": pressure_packet_verifier.WATERMARK,
                "qualification_effect": pressure_packet_verifier.QUALIFICATION_EFFECT,
                "consumption_rule": pressure_packet_verifier.CONSUMPTION_RULE,
                "publication_root_identity": publication_identity,
                "aggregate_receipt": aggregate_binding,
                "packet_root": str(packet_root),
                "inventory_sha256": hashlib.sha256(inventory_payload).hexdigest(),
                "source_bindings": source_bindings,
            }
        )
        self._write_planner_readonly(packet_receipt, packet_receipt_payload)
        packet_binding = {
            "path": str(packet_receipt),
            "sha256": hashlib.sha256(packet_receipt_payload).hexdigest(),
        }
        self._seal_planner_pressure_receipt(aggregate_receipt)
        self._seal_planner_pressure_receipt(packet_receipt)
        self._planner_pressure_aggregate_receipt = aggregate_receipt
        return aggregate_binding, packet_binding

    @staticmethod
    def _write_planner_readonly(path: Path, payload: bytes) -> None:
        path.write_bytes(payload)
        path.chmod(0o444)

    def _seal_planner_pressure_receipt(self, receipt: Path) -> None:
        publication = self.pic_root / "publication"
        metadata = receipt.stat()
        payload = receipt.read_bytes()
        seal = self.pic_root / "publication_acceptance" / (
            f".{receipt.name}.publication-success"
        )
        self._write_planner_readonly(
            seal,
            self._planner_json_payload(
                {
                    "schema_version": 1,
                    "record_type": pressure_packet_verifier.SUCCESS_SEAL_RECORD_TYPE,
                    "publication_root_identity": {
                        "device": publication.stat().st_dev,
                        "inode": publication.stat().st_ino,
                    },
                    "receipt_name": receipt.name,
                    "receipt_sha256": hashlib.sha256(payload).hexdigest(),
                    "receipt_identity": {
                        "device": metadata.st_dev,
                        "inode": metadata.st_ino,
                    },
                }
            ),
        )

    def _alternate_planner_pressure_aggregate_binding(self) -> dict[str, str]:
        source = self._planner_pressure_aggregate_receipt
        alternate = source.with_name("alternate-pressure-pilot-receipt.json")
        payload = source.read_bytes()
        self._write_planner_readonly(alternate, payload)
        self._seal_planner_pressure_receipt(alternate)
        return {
            "path": str(alternate),
            "sha256": hashlib.sha256(payload).hexdigest(),
        }

    def _prepare_planner_retention_reconciliation(
        self,
    ) -> tuple[Path, dict[str, object], Path]:
        self._write_science_config(
            authorize=True, planner_retention=self._planner_retention()
        )
        manifest_path = self._create_manifest()
        with patch(
            "validate_and_reserve_frontier_job._verify_clean_candidate",
            return_value=sha256(self._planner_candidate_manifest),
        ):
            reservation = self._reserve(manifest_path)
        self._attach(str(reservation["reservation_id"]))
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        artifact_dir = Path(str(manifest["artifact_dir"]))
        (artifact_dir / "analysis").mkdir(parents=True, mode=0o700)
        return manifest_path, reservation, artifact_dir

    @staticmethod
    def _planner_json_payload(value: object) -> bytes:
        return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")

    def _planner_pressure_gate_attestations(
        self,
        *,
        published_pressure_receipt: dict[str, str],
        published_packet_receipt: dict[str, str],
        pilot_bundle_manifest_sha256: str,
        aggregate_pilot_analysis_sha256: str,
        selected_case: dict[str, object],
        git_commit: str,
        source_archive_sha256: str,
        helper_sources: list[dict[str, str]],
    ) -> tuple[dict[str, str], dict[str, str]]:
        helper_by_path = {
            source["path"]: source["sha256"] for source in helper_sources
        }
        self.assertEqual(len(helper_by_path), len(helper_sources))
        reanalysis_source_closure = [
            {"path": path, "sha256": helper_by_path[path]}
            for path in pressure_packet_verifier.PRESSURE_REANALYSIS_SOURCE_PATHS
        ]
        source_closure_sha256 = hashlib.sha256(
            json.dumps(
                reanalysis_source_closure,
                separators=(",", ":"),
                sort_keys=True,
                allow_nan=False,
            ).encode("utf-8")
        ).hexdigest()
        timestamp = self._utc(datetime.now(timezone.utc))
        compact_timestamp = datetime.strptime(
            timestamp, "%Y-%m-%dT%H:%M:%SZ"
        ).strftime("%Y%m%dT%H%M%SZ")
        operator_id = "planner-fixture-operator"
        reviewer_id = "planner-fixture-reviewer"
        evidence = {
            "published_pressure_pilot_receipt": published_pressure_receipt,
            "published_pressure_pilot_review_packet_receipt": published_packet_receipt,
            "pilot_bundle_manifest_sha256": pilot_bundle_manifest_sha256,
            "aggregate_pilot_analysis_sha256": aggregate_pilot_analysis_sha256,
        }
        result = {
            "packet_receipt_sha256": published_packet_receipt["sha256"],
            "aggregate_receipt_sha256": published_pressure_receipt["sha256"],
            "manifest_sha256": pilot_bundle_manifest_sha256,
            "analysis_result_sha256": aggregate_pilot_analysis_sha256,
            "status": "pass_engineering_calibration_only",
        }
        reanalysis_binding = self._write_planner_pressure_gate_attestation(
            (
                f"{compact_timestamp}-q011-section54-pressure-reanalysis-"
                f"{operator_id}"
            ),
            {
                "schema_version": 1,
                "record_type": pressure_packet_verifier.PRESSURE_REANALYSIS_RECORD_TYPE,
                "qualification_effect": (
                    pressure_packet_verifier.PRESSURE_REANALYSIS_QUALIFICATION_EFFECT
                ),
                "operator_id": operator_id,
                "recomputed_utc": timestamp,
                "sealed_utc": timestamp,
                "operator_statement": (
                    pressure_packet_verifier.PRESSURE_REANALYSIS_OPERATOR_STATEMENT
                ),
                "evidence": evidence,
                "source_authorization": {
                    "execution_mode": (
                        pressure_packet_verifier.PRESSURE_REANALYSIS_EXECUTION_MODE
                    ),
                    "git_commit": git_commit,
                    "source_archive_sha256": source_archive_sha256,
                    "source_closure_sha256": source_closure_sha256,
                    "source_closure": reanalysis_source_closure,
                    "historical_production_source_authorization": dict(
                        pressure_packet_verifier.AUTHORIZED_HISTORICAL_REANALYSIS_SOURCE_AUTHORIZATION
                    ),
                },
                "result": result,
            },
        )
        reviewer_binding = self._write_planner_pressure_gate_attestation(
            (
                f"{compact_timestamp}-q011-section54-pressure-selection-"
                f"{reviewer_id}"
            ),
            {
                "schema_version": 1,
                "record_type": pressure_packet_verifier.PRESSURE_REVIEWER_RECORD_TYPE,
                "qualification_effect": (
                    pressure_packet_verifier.PRESSURE_REVIEWER_QUALIFICATION_EFFECT
                ),
                "selection_method": "human_review_only",
                "reviewer_id": reviewer_id,
                "reviewed_utc": timestamp,
                "sealed_utc": timestamp,
                "rationale": "Test-only retained human pressure selection.",
                "reviewer_statement": pressure_packet_verifier.PRESSURE_REVIEWER_STATEMENT,
                "published_pressure_pilot_receipt": published_pressure_receipt,
                "published_pressure_pilot_review_packet_receipt": (
                    published_packet_receipt
                ),
                "authoritative_reanalysis_attestation": reanalysis_binding,
                "selected_case": selected_case,
            },
        )
        return reanalysis_binding, reviewer_binding

    def _write_planner_pressure_gate_attestation(
        self,
        directory_name: str,
        attestation: dict[str, object],
    ) -> dict[str, str]:
        archive = (
            self.pic_root / pressure_packet_verifier.PRESSURE_GATE_ATTESTATION_ROOT_NAME
        )
        archive.mkdir(mode=0o700, exist_ok=True)
        directory = archive / directory_name
        directory.mkdir(mode=0o700)
        path = directory / pressure_packet_verifier.PRESSURE_GATE_ATTESTATION_FILENAME
        payload = self._planner_json_payload(attestation)
        path.write_bytes(payload)
        path.chmod(0o400)
        directory.chmod(0o500)
        return {"path": str(path), "sha256": hashlib.sha256(payload).hexdigest()}

    def _rewrite_planner_pressure_gate_attestation(
        self,
        binding: dict[str, str],
        mutate: Callable[[dict[str, object]], None],
    ) -> dict[str, str]:
        path = Path(binding["path"])
        directory = path.parent
        attestation = json.loads(path.read_text(encoding="utf-8"))
        mutate(attestation)
        payload = self._planner_json_payload(attestation)
        directory.chmod(0o700)
        path.chmod(0o600)
        path.write_bytes(payload)
        path.chmod(0o400)
        directory.chmod(0o500)
        return {"path": str(path), "sha256": hashlib.sha256(payload).hexdigest()}

    def _rewrite_planner_pressure_reanalysis_attestation(
        self,
        binding: dict[str, object],
        mutate_reanalysis: Callable[[dict[str, object]], None],
    ) -> None:
        planner_root = Path(str(binding["planner_root"]))
        pressure_receipt = json.loads(
            (
                planner_root / "bindings/human_pressure_selection_receipt.json"
            ).read_text(encoding="utf-8")
        )
        reanalysis_binding = self._rewrite_planner_pressure_gate_attestation(
            pressure_receipt["authoritative_reanalysis_attestation"],
            mutate_reanalysis,
        )

        def rebind_reviewer(reviewer: dict[str, object]) -> None:
            reviewer["authoritative_reanalysis_attestation"] = reanalysis_binding

        reviewer_binding = self._rewrite_planner_pressure_gate_attestation(
            pressure_receipt["reviewer_attestation"],
            rebind_reviewer,
        )

        def rebind_receipt(receipt: dict[str, object]) -> None:
            receipt["authoritative_reanalysis_attestation"] = reanalysis_binding
            receipt["reviewer_attestation"] = reviewer_binding

        self._rewrite_planner_pressure_receipt(binding, rebind_receipt)

    def _rewrite_closed_planner_tree(
        self,
        binding: dict[str, object],
        mutate: Callable[[dict[str, bytes]], None],
    ) -> dict[str, object]:
        planner_root = Path(str(binding["planner_root"]))
        planner_root.chmod(0o755)
        for path in planner_root.rglob("*"):
            path.chmod(0o755 if path.is_dir() else 0o644)
        files = {
            path.relative_to(planner_root).as_posix(): path.read_bytes()
            for path in planner_root.rglob("*")
            if path.is_file() and path.name != "artifact_inventory.sha256"
        }
        mutate(files)
        receipt = json.loads(files["materialization_receipt.json"])
        pre_receipt_inventory = "".join(
            f"{hashlib.sha256(files[path]).hexdigest()}  {path}\n"
            for path in sorted(files)
            if path
            not in {
                "materialization_receipt.json",
                "freeze_receipt.json",
                "artifact_inventory.sha256",
            }
        ).encode("utf-8")
        receipt["tree_inventory"]["sha256"] = hashlib.sha256(
            pre_receipt_inventory
        ).hexdigest()
        receipt["tree_inventory"]["inventoried_file_count"] = len(
            [
                path
                for path in files
                if path not in {"materialization_receipt.json", "freeze_receipt.json"}
            ]
        )
        receipt_payload = self._planner_json_payload(receipt)
        files["materialization_receipt.json"] = receipt_payload
        inventory_payload = "".join(
            f"{hashlib.sha256(files[path]).hexdigest()}  {path}\n"
            for path in sorted(files)
        ).encode("utf-8")
        for relative, payload in files.items():
            (planner_root / relative).write_bytes(payload)
        (planner_root / "artifact_inventory.sha256").write_bytes(inventory_payload)
        for path in sorted(planner_root.rglob("*"), reverse=True):
            path.chmod(0o555 if path.is_dir() else 0o444)
        planner_root.chmod(0o555)
        binding["planner_inventory_sha256"] = hashlib.sha256(
            inventory_payload
        ).hexdigest()
        binding["planner_materialization_receipt"] = {
            "path": "materialization_receipt.json",
            "sha256": hashlib.sha256(receipt_payload).hexdigest(),
        }
        return binding

    def _rebind_planner_campaign_plan(self, files: dict[str, bytes]) -> None:
        receipt = json.loads(files["materialization_receipt.json"])
        receipt["campaign_plan"]["sha256"] = hashlib.sha256(
            files["campaign_plan.json"]
        ).hexdigest()
        files["materialization_receipt.json"] = self._planner_json_payload(receipt)

    def _rewrite_planner_pressure_receipt(
        self,
        binding: dict[str, object],
        mutate_receipt: Callable[[dict[str, object]], None],
    ) -> None:
        def mutate(files: dict[str, bytes]) -> None:
            path = "bindings/human_pressure_selection_receipt.json"
            pressure_receipt = json.loads(files[path])
            mutate_receipt(pressure_receipt)
            pressure_payload = self._planner_json_payload(pressure_receipt)
            files[path] = pressure_payload
            pressure_sha256 = hashlib.sha256(pressure_payload).hexdigest()
            plan = json.loads(files["campaign_plan.json"])
            plan["source_bindings"]["pressure_selection_receipt"][
                "sha256"
            ] = pressure_sha256
            plan["selected_pressure"]["receipt"]["sha256"] = pressure_sha256
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)

    def test_planner_retention_rejects_caller_selected_continuation_namespace(
        self,
    ) -> None:
        binding = self._planner_retention()
        self.assertEqual(
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            ),
            binding,
        )
        attempt_id = "q011-restart-fixture-001"
        attempt_root = (
            self.pic_root
            / "campaigns"
            / f"q011-section54-{binding['planner_plan_id']}"
            / "restart_continuation"
            / attempt_id
        )
        binding["attempt_id"] = attempt_id
        binding["authorized_orion_attempt_root"] = str(
            attempt_root
        )
        binding["authorized_orion_raw_root"] = str(attempt_root / "raw")
        binding["argv"] = ["-r", "/retained/source.rst", "-d", str(attempt_root / "raw")]
        with self.assertRaisesRegex(ValueError, "does not select one immutable descriptor"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_shaped_caller_overlay_substitution(self) -> None:
        accepted = self._planner_retention()
        for key, value in [
            ("authorized_orion_attempt_root", "/tmp/operator-selected-attempt"),
            ("authorized_orion_raw_root", "/tmp/operator-selected-attempt/raw"),
            ("argv", ["-i", "bindings/operator.athinput", "-d", "/tmp/raw"]),
        ]:
            forged = json.loads(json.dumps(accepted))
            forged[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(
                ValueError, "differs from immutable planner bytes"
            ):
                validate_planner_retention_binding(
                    forged, authorized_pic_root=self.pic_root
                )
        forged = json.loads(json.dumps(accepted))
        forged["planner_root"] = str(self.pic_root / "plans" / "operator-selected")
        with self.assertRaisesRegex(ValueError, "planner root is malformed"):
            validate_planner_retention_binding(
                forged, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_submission_clean_candidate_mismatch(self) -> None:
        binding = self._planner_retention()
        with self.assertRaisesRegex(
            ValueError, "differs from submission binding"
        ):
            validate_planner_retention_binding(
                binding,
                authorized_pic_root=self.pic_root,
                expected_clean_candidate_manifest_sha256="0" * 64,
            )

    def test_q011_helper_source_order_matches_execution_and_analyzer(self) -> None:
        import ast
        import control_plane_common as common

        repo_root = Path(__file__).resolve().parents[3]

        def source_tuple(relative: str, name: str) -> tuple[str, ...]:
            tree = ast.parse((repo_root / relative).read_text(encoding="utf-8"))
            for node in tree.body:
                if (
                    isinstance(node, ast.Assign)
                    and len(node.targets) == 1
                    and isinstance(node.targets[0], ast.Name)
                    and node.targets[0].id == name
                ):
                    value = ast.literal_eval(node.value)
                    self.assertIsInstance(value, tuple)
                    return value
            self.fail(f"Missing {name} in {relative}")

        expected = common.Q011_SECTION54_HELPER_SOURCES
        self.assertEqual(
            expected,
            source_tuple(
                "tst/publication/q011_section54_qualifying_campaign_execution.py",
                "_FIXED_HELPER_SOURCES",
            ),
        )
        self.assertEqual(
            expected,
            source_tuple(
                "tst/publication/analyze_q011_section54_campaign.py",
                "_EXPECTED_HELPER_SOURCE_PATHS",
            ),
        )

    def test_planner_retention_rejects_self_authored_frozen_helper_digest(self) -> None:
        binding = self._planner_retention()

        def mutate(files: dict[str, bytes]) -> None:
            helper = json.loads(files["helper_source_closure.json"])
            helper["sources"][0]["sha256"] = "f" * 64
            helper_payload = self._planner_json_payload(helper)
            files["helper_source_closure.json"] = helper_payload
            helper_sha256 = hashlib.sha256(helper_payload).hexdigest()
            plan = json.loads(files["campaign_plan.json"])
            plan["helper_source_closure"]["sha256"] = helper_sha256
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            receipt = json.loads(files["materialization_receipt.json"])
            receipt["helper_source_closure"]["sha256"] = helper_sha256
            files["materialization_receipt.json"] = self._planner_json_payload(receipt)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)
        with self.assertRaisesRegex(ValueError, "reviewed archive bytes"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_self_authored_arbitrary_matrix(self) -> None:
        binding = self._planner_retention()

        def mutate(files: dict[str, bytes]) -> None:
            plan = json.loads(files["campaign_plan.json"])
            plan["campaign_matrix"]["physical_mode"] = "operator_authored_mode"
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)
        with self.assertRaisesRegex(ValueError, "Section 5.4 matrix"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_self_authored_arbitrary_candidate(self) -> None:
        binding = self._planner_retention()

        def mutate(files: dict[str, bytes]) -> None:
            plan = json.loads(files["campaign_plan.json"])
            plan["candidate_binding"]["git_commit"] = "a" * 40
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)
        with self.assertRaisesRegex(ValueError, "candidate binding drifted"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_replayed_pressure_reanalysis_source_snapshot(
        self,
    ) -> None:
        binding = self._planner_retention()
        planner_root = Path(str(binding["planner_root"]))
        pressure_receipt = json.loads(
            (
                planner_root / "bindings/human_pressure_selection_receipt.json"
            ).read_text(encoding="utf-8")
        )
        original_reanalysis = json.loads(
            Path(
                pressure_receipt["authoritative_reanalysis_attestation"]["path"]
            ).read_text(encoding="utf-8")
        )

        def alternate(value: str, length: int) -> str:
            candidate = "0" * length
            return candidate if value != candidate else "f" * length

        def drift_commit(authorization: dict[str, object]) -> None:
            authorization["git_commit"] = alternate(
                str(authorization["git_commit"]), 40
            )

        def drift_archive(authorization: dict[str, object]) -> None:
            authorization["source_archive_sha256"] = alternate(
                str(authorization["source_archive_sha256"]), 64
            )

        def drift_closure(authorization: dict[str, object]) -> None:
            closure = authorization["source_closure"]
            closure[0]["sha256"] = alternate(str(closure[0]["sha256"]), 64)
            authorization["source_closure_sha256"] = hashlib.sha256(
                json.dumps(
                    closure,
                    separators=(",", ":"),
                    sort_keys=True,
                    allow_nan=False,
                ).encode("utf-8")
            ).hexdigest()

        for label, drift in (
            ("git commit", drift_commit),
            ("source archive", drift_archive),
            ("source closure", drift_closure),
        ):
            def replay_from_alternate_source(
                attestation: dict[str, object],
                *,
                drift: Callable[[dict[str, object]], None] = drift,
            ) -> None:
                attestation.clear()
                attestation.update(json.loads(json.dumps(original_reanalysis)))
                drift(attestation["source_authorization"])

            with self.subTest(binding=label):
                self._rewrite_planner_pressure_reanalysis_attestation(
                    binding, replay_from_alternate_source
                )
                with self.assertRaisesRegex(
                    ValueError,
                    "Planner pressure reanalysis source-snapshot binding failed",
                ):
                    validate_planner_retention_binding(
                        binding, authorized_pic_root=self.pic_root
                    )

    def test_planner_retention_rejects_self_authored_reduced_pressure_receipt(
        self,
    ) -> None:
        binding = self._planner_retention()

        def mutate(files: dict[str, bytes]) -> None:
            path = "bindings/human_pressure_selection_receipt.json"
            receipt = json.loads(files[path])
            receipt.pop("reviewer_attestation")
            receipt_payload = self._planner_json_payload(receipt)
            files[path] = receipt_payload
            receipt_sha256 = hashlib.sha256(receipt_payload).hexdigest()
            plan = json.loads(files["campaign_plan.json"])
            plan["source_bindings"]["pressure_selection_receipt"][
                "sha256"
            ] = receipt_sha256
            plan["selected_pressure"]["receipt"]["sha256"] = receipt_sha256
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)
        with self.assertRaisesRegex(ValueError, "receipt schema drifted"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_missing_pressure_packet_binding(self) -> None:
        binding = self._planner_retention()

        def mutate(files: dict[str, bytes]) -> None:
            path = "bindings/human_pressure_selection_receipt.json"
            receipt = json.loads(files[path])
            receipt.pop("published_pressure_pilot_review_packet_receipt")
            receipt_payload = self._planner_json_payload(receipt)
            files[path] = receipt_payload
            receipt_sha256 = hashlib.sha256(receipt_payload).hexdigest()
            plan = json.loads(files["campaign_plan.json"])
            plan["source_bindings"]["pressure_selection_receipt"][
                "sha256"
            ] = receipt_sha256
            plan["selected_pressure"]["receipt"]["sha256"] = receipt_sha256
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)
        with self.assertRaisesRegex(ValueError, "receipt schema drifted"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_pressure_packet_hash_drift(self) -> None:
        binding = self._planner_retention()

        def mutate(files: dict[str, bytes]) -> None:
            path = "bindings/human_pressure_selection_receipt.json"
            receipt = json.loads(files[path])
            receipt["published_pressure_pilot_review_packet_receipt"][
                "sha256"
            ] = "0" * 64
            receipt_payload = self._planner_json_payload(receipt)
            files[path] = receipt_payload
            receipt_sha256 = hashlib.sha256(receipt_payload).hexdigest()
            plan = json.loads(files["campaign_plan.json"])
            plan["source_bindings"]["pressure_selection_receipt"][
                "sha256"
            ] = receipt_sha256
            plan["selected_pressure"]["receipt"]["sha256"] = receipt_sha256
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)
        with self.assertRaisesRegex(
            ValueError,
            "Planner human pressure-selection review packet verifier result drifted",
        ):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_cross_bound_pressure_packet_aggregate(
        self,
    ) -> None:
        binding = self._planner_retention()
        alternate_aggregate = self._alternate_planner_pressure_aggregate_binding()

        def mutate(files: dict[str, bytes]) -> None:
            path = "bindings/human_pressure_selection_receipt.json"
            receipt = json.loads(files[path])
            receipt["published_pressure_pilot_receipt"] = alternate_aggregate
            receipt_payload = self._planner_json_payload(receipt)
            files[path] = receipt_payload
            receipt_sha256 = hashlib.sha256(receipt_payload).hexdigest()
            plan = json.loads(files["campaign_plan.json"])
            plan["source_bindings"]["pressure_selection_receipt"][
                "sha256"
            ] = receipt_sha256
            plan["selected_pressure"]["receipt"]["sha256"] = receipt_sha256
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)
        with self.assertRaisesRegex(
            ValueError, "does not bind the supplied aggregate receipt"
        ):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_pressure_bundle_manifest_hash_drift(
        self,
    ) -> None:
        binding = self._planner_retention()
        self._rewrite_planner_pressure_receipt(
            binding,
            lambda receipt: receipt.__setitem__(
                "pilot_bundle_manifest_sha256", "0" * 64
            ),
        )
        with self.assertRaisesRegex(ValueError, "aggregate evidence binding drifted"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_pressure_aggregate_analysis_hash_drift(
        self,
    ) -> None:
        binding = self._planner_retention()
        self._rewrite_planner_pressure_receipt(
            binding,
            lambda receipt: receipt.__setitem__(
                "aggregate_pilot_analysis_sha256", "0" * 64
            ),
        )
        with self.assertRaisesRegex(ValueError, "aggregate evidence binding drifted"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_pressure_aggregate_descriptor_drift(
        self,
    ) -> None:
        binding = self._planner_retention()

        def mutate(receipt: dict[str, object]) -> None:
            receipt["case_descriptors"][0]["descriptor_sha256"] = "0" * 64

        self._rewrite_planner_pressure_receipt(binding, mutate)
        with self.assertRaisesRegex(
            ValueError, "aggregate descriptor binding drifted"
        ):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_self_authored_extra_argv(self) -> None:
        binding = self._planner_retention()

        def mutate(files: dict[str, bytes]) -> None:
            plan = json.loads(files["campaign_plan.json"])
            descriptor_binding = plan["baseline_attempt_descriptors"][0]
            descriptor_path = descriptor_binding["path"]
            descriptor = json.loads(files[descriptor_path])
            contract_path = descriptor["launch_contract"]["path"]
            contract = json.loads(files[contract_path])
            contract["argv"].append("mesh/nx1=1")
            contract_payload = self._planner_json_payload(contract)
            files[contract_path] = contract_payload
            descriptor["launch_contract"]["sha256"] = hashlib.sha256(
                contract_payload
            ).hexdigest()
            descriptor_payload = self._planner_json_payload(descriptor)
            files[descriptor_path] = descriptor_payload
            descriptor_binding["sha256"] = hashlib.sha256(
                descriptor_payload
            ).hexdigest()
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)
        with self.assertRaisesRegex(ValueError, "launch contract.*drifted"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_self_authored_restart_carrier(self) -> None:
        binding = self._planner_retention()

        def mutate(files: dict[str, bytes]) -> None:
            path = "restart_continuation/amr_restart_continuation_carrier.json"
            carrier = json.loads(files[path])
            carrier["checkpoint_nominal_slot_omega0_inverse"] = 501.0
            carrier_payload = self._planner_json_payload(carrier)
            files[path] = carrier_payload
            plan = json.loads(files["campaign_plan.json"])
            plan["restart_continuation_carrier"]["sha256"] = hashlib.sha256(
                carrier_payload
            ).hexdigest()
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)
        with self.assertRaisesRegex(ValueError, "restart-continuation carrier drifted"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_self_authored_recompute_plan(self) -> None:
        binding = self._planner_retention()

        def mutate(files: dict[str, bytes]) -> None:
            path = "independent_raw_artifact_recompute_plan.json"
            recompute = json.loads(files[path])
            recompute["production_helper_imports_authorized"] = True
            recompute_payload = self._planner_json_payload(recompute)
            files[path] = recompute_payload
            plan = json.loads(files["campaign_plan.json"])
            plan["independent_raw_artifact_recompute_plan"]["sha256"] = hashlib.sha256(
                recompute_payload
            ).hexdigest()
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)
        with self.assertRaisesRegex(ValueError, "raw-artifact recompute plan"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def test_planner_retention_rejects_self_authored_policy_fragment(self) -> None:
        binding = self._planner_retention()

        def mutate(files: dict[str, bytes]) -> None:
            path = "nonauthorizing_policy_fragment.json"
            fragment = json.loads(files[path])
            fragment["launch_authorized"] = True
            fragment_payload = self._planner_json_payload(fragment)
            files[path] = fragment_payload
            plan = json.loads(files["campaign_plan.json"])
            plan["nonauthorizing_policy_fragment"]["sha256"] = hashlib.sha256(
                fragment_payload
            ).hexdigest()
            files["campaign_plan.json"] = self._planner_json_payload(plan)
            self._rebind_planner_campaign_plan(files)

        self._rewrite_closed_planner_tree(binding, mutate)
        with self.assertRaisesRegex(ValueError, "nonauthorizing policy fragment"):
            validate_planner_retention_binding(
                binding, authorized_pic_root=self.pic_root
            )

    def _clean_source(self, name: str) -> Path:
        source_root = self.root / name
        source_root.mkdir()
        subprocess.run(["git", "init", str(source_root)], check=True, capture_output=True)
        (source_root / "tracked.txt").write_text("tracked\n", encoding="utf-8")
        deck_root = source_root / "inputs/tests"
        deck_root.mkdir(parents=True)
        deck = deck_root / "pic_paper.athinput"
        deck.write_text("<job>\nbasename = prepared-paper\n", encoding="utf-8")
        publication_deck = (
            source_root
            / "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
        )
        publication_deck.parent.mkdir()
        publication_deck.write_text(
            "<job>\nbasename = prepared-section54-paper\n", encoding="utf-8"
        )
        analyzer_root = source_root / "tst/publication"
        analyzer_root.mkdir(parents=True)
        analyzer = analyzer_root / "analyze_paper.py"
        analyzer.write_text("print('prepared analysis')\n", encoding="utf-8")
        inventory = (
            analyzer_root / "frontier_control_plane/prepared_pic_artifact_inventory.json"
        )
        inventory.parent.mkdir()
        inventory.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "paper_decks": [
                        {
                            "path": (
                                "inputs/publication/"
                                "pic_parallel_shock_section54_paper_vl2_tsc.athinput"
                            ),
                            "sha256": sha256(publication_deck),
                        },
                        {"path": "inputs/tests/pic_paper.athinput", "sha256": sha256(deck)}
                    ],
                    "analyzers": [
                        {
                            "path": "tst/publication/analyze_paper.py",
                            "sha256": sha256(analyzer),
                        }
                    ],
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        subprocess.run(["git", "-C", str(source_root), "add", "."], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "-c",
                "user.name=PIC Test",
                "-c",
                "user.email=pic-test@example.invalid",
                "commit",
                "-m",
                "clean candidate",
            ],
            check=True,
            capture_output=True,
        )
        return source_root

    def _commit_source_change(
        self, source_root: Path, *, relative: str, content: str, message: str
    ) -> None:
        path = source_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        subprocess.run(["git", "-C", str(source_root), "add", "--", relative], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "-c",
                "user.name=PIC Test",
                "-c",
                "user.email=pic-test@example.invalid",
                "commit",
                "-m",
                message,
            ],
            check=True,
            capture_output=True,
        )

    def _tag_source_head(self, source_root: Path, tag: str) -> None:
        subprocess.run(["git", "-C", str(source_root), "tag", tag], check=True)

    def _prepared_artifact_inventory(self) -> str:
        return "tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json"

    def _profile_writer_arguments(
        self, source_root: Path, profile_id: str = "test-profile"
    ) -> dict[str, object]:
        commit = subprocess.check_output(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"], text=True
        ).strip()
        artifact_dir = self.pic_root / "bin" / commit[:12] / profile_id
        artifact_dir.mkdir(parents=True)
        log_dir = self.pic_root / "logs" / "build"
        log_dir.mkdir(parents=True, exist_ok=True)
        executable = artifact_dir / "athena"
        executable.write_text("built executable\n", encoding="utf-8")
        provenance = {
            "configure_log": log_dir / f"{commit[:12]}.{profile_id}.configure.log",
            "build_log": log_dir / f"{commit[:12]}.{profile_id}.build.log",
            "cmake_cache": artifact_dir / "CMakeCache.txt",
            "module_list": artifact_dir / "modules.txt",
            "toolchain_file": artifact_dir / "toolchain.txt",
            "build_invocations_file": artifact_dir / "build-invocations.json",
            "git_status_preconfigure_file": artifact_dir / "git_status.preconfigure.txt",
            "git_status_file": artifact_dir / "git_status.txt",
            "submodule_status_file": artifact_dir / "submodule_status.txt",
            "environment_allowlist_file": artifact_dir / "environment.allowlist.txt",
            "build_environment_file": artifact_dir / "build-environment.json",
        }
        for key, path in provenance.items():
            if key in {"git_status_preconfigure_file", "git_status_file"}:
                text = subprocess.check_output(
                    [
                        "git",
                        "-C",
                        str(source_root),
                        "status",
                        "--ignore-submodules=none",
                        "--porcelain",
                        "--untracked-files=all",
                    ],
                    text=True,
                )
            elif key == "submodule_status_file":
                text = subprocess.check_output(
                    ["git", "-C", str(source_root), "submodule", "status", "--recursive"],
                    text=True,
                )
            elif key == "build_invocations_file":
                text = json.dumps(
                    {"configure": ["/fake/cmake", "-S", "source"], "build": ["/fake/cmake", "--build", "build"]}
                )
            elif key == "build_environment_file":
                text = "{}"
            else:
                text = f"{key}=reviewed\n"
            path.write_text(text, encoding="utf-8")
        return {
            "source_root": source_root,
            "fresh_source_root": source_root,
            "executable": executable,
            "output": artifact_dir / "build_profile.json",
            "profile_id": profile_id,
            "expected_git_commit": commit,
            "configure_log": provenance["configure_log"],
            "build_log": provenance["build_log"],
            "cmake_cache": provenance["cmake_cache"],
            "module_list": provenance["module_list"],
            "toolchain_file": provenance["toolchain_file"],
            "build_invocations_file": provenance["build_invocations_file"],
            "git_status_preconfigure_file": provenance["git_status_preconfigure_file"],
            "git_status_file": provenance["git_status_file"],
            "submodule_status_file": provenance["submodule_status_file"],
            "environment_allowlist_file": provenance["environment_allowlist_file"],
            "build_environment_file": provenance["build_environment_file"],
            "control_plane_dir": self.control_plane_dir,
            "authorized_pic_root": self.pic_root,
            "authorized_source_root": source_root,
        }

    def _build_profile(self, source_root: Path, build: Path, profile_id: str) -> tuple[Path, Path]:
        del build
        commit = subprocess.check_output(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"], text=True
        ).strip()

        def execute(command: list[str], *, stream: object, environment: dict[str, str]) -> None:
            del environment
            stream.write(b"fake cmake invocation\n")
            if "--build" in command:
                cmake_dir = Path(command[command.index("--build") + 1])
                (cmake_dir / "src").mkdir(parents=True)
                (cmake_dir / "src" / "athena").write_text(
                    "built executable\n", encoding="utf-8"
                )
            else:
                cmake_dir = Path(command[command.index("-B") + 1])
                cmake_dir.mkdir(parents=True)
                (cmake_dir / "CMakeCache.txt").write_text(
                    "fixture cache\n", encoding="utf-8"
                )

        with patch("write_orion_build_profile._execute_logged_command", side_effect=execute):
            profile = build_orion_profile(
                source_root=source_root,
                expected_git_commit=commit,
                profile_id=profile_id,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_source_root=source_root,
            )
        executable = profile.with_name("athena")
        return executable, profile

    def _add_submodule(self, source_root: Path, name: str) -> Path:
        nested = self._clean_source(f"{source_root.name}-{name}-source")
        self._add_existing_submodule(source_root, nested, name)
        return source_root / name

    def _add_existing_submodule(
        self, source_root: Path, nested: Path, name: str
    ) -> None:
        subprocess.run(
            [
                "git",
                "-c",
                "protocol.file.allow=always",
                "-C",
                str(source_root),
                "submodule",
                "add",
                str(nested),
                name,
            ],
            check=True,
            capture_output=True,
        )
        subprocess.run(["git", "-C", str(source_root), "add", "."], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "-c",
                "user.name=PIC Test",
                "-c",
                "user.email=pic-test@example.invalid",
                "commit",
                "-m",
                f"add {name} source",
            ],
            check=True,
            capture_output=True,
        )

    def _clean_candidate(
        self,
        *,
        authorize: bool,
        source_name: str = "candidate-source",
        freeze_id: str = "03a7bd9a-7d4c-4e37-a12b-46de3817eff2",
        profile_id: str = "hip-mpi-release-paper-pic",
    ) -> tuple[Path, Path, str]:
        source_root = self._clean_source(source_name)
        self.authorized_clean_candidate_source_root = source_root
        self._add_submodule(source_root, "nested")
        executable, profile = self._build_profile(
            source_root, self.pic_root / "candidate-build", profile_id
        )
        manifest_path = create_freeze(
            source_root=source_root,
            executable=executable,
            build_profile=profile,
            build_profile_id=profile_id,
            prepared_artifact_inventory=self._prepared_artifact_inventory(),
            freeze_id=freeze_id,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_source_root=source_root,
        )
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        frozen_executable = Path(str(manifest["build"]["executable_path"]))
        git_commit = str(manifest["source"]["git_commit"])
        if authorize:
            self._write_policy(
                science_submission_freeze=self._authorized_science_freeze(
                    manifest_path
                )
            )
            self._promote_policy()
        return manifest_path, frozen_executable, git_commit

    def _authorized_science_freeze(
        self,
        manifest: Path,
        *,
        manifest_sha256: str | None = None,
        build_profile_control_plane_version: str | None = None,
    ) -> dict[str, object]:
        return {
            "status": "authorized",
            "manifest_path": str(manifest),
            "manifest_sha256": manifest_sha256 or sha256(manifest),
            "build_profile_control_plane_version": (
                build_profile_control_plane_version or self.control_plane_version
            ),
        }

    def _write_science_config(self, *, authorize: bool, **overrides: object) -> Path:
        planner_manifest = getattr(self, "_planner_candidate_manifest", None)
        if overrides.get("planner_retention") is None or planner_manifest is None:
            manifest, executable, git_commit = self._clean_candidate(authorize=authorize)
        else:
            manifest = planner_manifest
            planner_candidate = json.loads(manifest.read_text(encoding="utf-8"))
            executable = Path(str(planner_candidate["build"]["executable_path"]))
            git_commit = str(planner_candidate["source"]["git_commit"])
        config = {
            "campaign": "f1_gpu_gyro",
            "test_id": "pic_relativistic_gyro_paper",
            "submission_scope": "registered_science",
            "registered_science_authorization_id": "f1-clean-gyro-v1",
            "git_commit": git_commit,
            "evidence_class": "frontier_f1_registered_science",
            "physical_mode": "paper_test_particle",
            "executable": str(executable),
            "clean_candidate_manifest": str(manifest),
            "artifact_dir": str(
                self.pic_root / "runs" / "f1_gpu_gyro" / self.submission_id
            ),
        }
        config.update(overrides)
        self._write_config(**config)
        if authorize:
            self.science_submission_freeze = self._authorized_science_freeze(manifest)
            self.registered_science_slices = [
                {
                    "authorization_id": config["registered_science_authorization_id"],
                    "status": "authorized",
                    "campaign": config["campaign"],
                    "test_id": config["test_id"],
                    "evidence_class": config["evidence_class"],
                    "physical_mode": config["physical_mode"],
                    "runtime_profile": "frontier_minimum_supported",
                    "selected_qos": "debug",
                    "registered_short_nonproduction": True,
                    "maximum_nodes": 1,
                    "maximum_walltime_seconds": 10 * 60,
                    "maximum_attempts": 1,
                    "job_script_sha256": sha256(self.sources / "job.sh"),
                    "input_deck_sha256": sha256(self.sources / "input.athinput"),
                    "environment_profile_sha256": sha256(self.sources / "environment.sh"),
                    "analysis_script_sha256": [sha256(self.sources / "analysis.py")],
                    "executable_sha256": sha256(executable),
                    "launch_contract_sha256": launch_contract_sha256(
                        self._launch_contract()
                    ),
                    "clean_candidate_manifest_sha256": sha256(manifest),
                }
            ]
            self._write_policy(
                science_submission_freeze=self.science_submission_freeze,
                admission_smoke_overrides={"status": "closed_after_pass"},
            )
            self._promote_policy()
        return manifest

    def _update_registered_science_candidate_sha(self, candidate: Path) -> None:
        for record in self.registered_science_slices:
            record["clean_candidate_manifest_sha256"] = sha256(candidate)

    def _fresh_submission_manifest(self) -> Path:
        config = json.loads(self.config.read_text(encoding="utf-8"))
        submission_id = str(uuid.uuid4())
        config["submission_id"] = submission_id
        config["artifact_dir"] = str(
            self.pic_root / "runs" / str(config["campaign"]) / submission_id
        )
        self.config.write_text(json.dumps(config), encoding="utf-8")
        return self._create_manifest()

    def _rewrite_clean_candidate_profile(
        self, candidate: Path, profile: dict[str, object]
    ) -> None:
        candidate_value = json.loads(candidate.read_text(encoding="utf-8"))
        build = candidate_value["build"]
        self.assertIsInstance(build, dict)
        frozen_profile = Path(str(build["profile_path"]))
        candidate.parent.chmod(0o755)
        frozen_profile.chmod(0o644)
        frozen_profile.write_text(json.dumps(profile), encoding="utf-8")
        frozen_profile.chmod(0o444)
        build["profile_sha256"] = sha256(frozen_profile)
        candidate.chmod(0o644)
        candidate.write_text(json.dumps(candidate_value), encoding="utf-8")
        candidate.chmod(0o444)
        candidate.parent.chmod(0o555)
        self._update_registered_science_candidate_sha(candidate)
        self._write_policy(
            science_submission_freeze=self._authorized_science_freeze(candidate),
            admission_smoke_overrides={"status": "closed_after_pass"},
        )
        self._install_deliberately_malformed_active_policy_fixture()

    def _create_manifest(self, *, control_plane_dir: Path | None = None) -> Path:
        target_control_plane = control_plane_dir or self.control_plane_dir
        config = json.loads(self.config.read_text(encoding="utf-8"))
        if config.get("submission_scope") == "registered_science":
            config["pre_manifest_attestation"] = str(
                self._sealed_operator_attestation(
                    str(config["registered_science_authorization_id"]),
                    "pre_manifest",
                    control_plane_version=target_control_plane.name,
                )
            )
            self.config.write_text(json.dumps(config), encoding="utf-8")
        return create_manifest(
            self.config,
            control_plane_dir=target_control_plane,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )

    def _reserve(
        self,
        manifest_path: Path,
        cap: float = 10000.0,
        reservation_id: str = "89c76745-6c37-47f7-9847-800a98a47c9b",
        control_plane_dir: Path | None = None,
        *,
        patch_clean_candidate_bundle: bool = True,
        inject_wrapper_attestation: bool = True,
    ) -> dict[str, object]:
        def validate_with_test_roots(
            candidate: dict[str, object], **kwargs: object
        ) -> list[dict[str, str]]:
            if self.authorized_clean_candidate_source_root is None:
                raise AssertionError("Clean-candidate test source root was not injected")
            return validate_clean_candidate_bundle(
                candidate,
                **kwargs,
                authorized_pic_root=self.pic_root,
                authorized_source_root=self.authorized_clean_candidate_source_root,
            )

        def invoke() -> dict[str, object]:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            wrapper_attestation = None
            if (
                inject_wrapper_attestation
                and manifest.get("submission_scope") == "registered_science"
            ):
                wrapper_attestation = self._sealed_operator_attestation(
                    str(manifest["registered_science_authorization_id"]),
                    "pre_submit_wrapper",
                    control_plane_version=(control_plane_dir or self.control_plane_dir).name,
                )
            return reserve(
                manifest_path=manifest_path,
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                node_hour_cap=cap,
                reservation_id=reservation_id,
                pre_submit_wrapper_attestation=wrapper_attestation,
                control_plane_dir=control_plane_dir or self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

        if not patch_clean_candidate_bundle:
            return invoke()
        with patch(
            "validate_and_reserve_frontier_job.validate_clean_candidate_bundle",
            side_effect=validate_with_test_roots,
        ):
            return invoke()

    def _attach(self, reservation_id: str, job_id: str = "12345") -> None:
        scheduler = (
            f"JobId={job_id} JobState=PENDING Account=AST207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            mark_submitted(
                reservation_id=reservation_id,
                job_id=job_id,
                ledger_jsonl=self.ledger,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
            transition(
                reservation_id=reservation_id,
                job_id=job_id,
                event_type="job_id_attached",
                state="submitted",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def _mark_dispatch_started(self, reservation_id: str) -> None:
        mark_dispatch_started(
            reservation_id=reservation_id,
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            receipts_jsonl=self.receipts,
            mirror_jsonl=self.mirror,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )

    def _launch(
        self,
        manifest_path: Path,
        reservation: dict[str, object],
        *,
        runner: object = subprocess.run,
        environment_overrides: dict[str, str] | None = None,
        flock_error: OSError | None = None,
    ) -> None:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        job_script = record_for_role(manifest, "job-script")
        executable = record_for_role(manifest, "executable")
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=RUNNING Account=AST207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        self._attach(reservation_id)
        environment = {
            "PIC_MANIFEST_SHA256": str(reservation["manifest_sha256"]),
            "PIC_RESERVATION_ID": reservation_id,
            "PIC_SUBMISSION_ID": self.submission_id,
            "SLURM_JOB_ID": "12345",
        }
        environment.update(environment_overrides or {})
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            with patch.dict(
                os.environ,
                environment,
                clear=True,
            ):
                def execute() -> None:
                    launch(
                        manifest_path=manifest_path,
                        manifest_sha256=str(reservation["manifest_sha256"]),
                        job_script_sha256=str(job_script["sha256"]),
                        executable_sha256=str(executable["sha256"]),
                        reservation_id=reservation_id,
                        submission_id=self.submission_id,
                        ledger_jsonl=self.ledger,
                        receipts_jsonl=self.receipts,
                        mirror_jsonl=self.mirror,
                        runner=runner,
                        control_plane_dir=self.control_plane_dir,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                    )

                if flock_error is None:
                    execute()
                else:
                    with patch("ledger.fcntl.flock", side_effect=flock_error):
                        execute()

    def test_snapshot_verifies_and_detects_mutation(self) -> None:
        manifest_path = self._create_manifest()
        verify(manifest_path, submission_id=self.submission_id)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        executable = Path(str(record_for_role(manifest, "executable")["path"]))
        self.assertTrue(executable.stat().st_mode & 0o111)
        deck = Path(str(record_for_role(manifest, "input-deck")["path"]))
        deck.chmod(0o644)
        deck.write_text("mutated\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            verify(manifest_path)

    def test_snapshot_rejects_unresolved_placeholder(self) -> None:
        self._write("input.athinput", "value = REPLACE_INPUT\n")
        with self.assertRaises(ValueError):
            self._create_manifest()
        campaign_dir = self.pic_root / "manifests" / "f0_hipmpi_smoke"
        self.assertFalse(list(campaign_dir.glob(".tmp-*")))

    def test_snapshot_rejects_test_id_escape_before_overwriting_external_file(self) -> None:
        outside = self.root / "outside"
        outside.mkdir()
        victim = outside / "victim.athinput"
        victim.write_text("preserve external file\n", encoding="utf-8")
        self._write_config(test_id=str(outside / "victim"))
        with self.assertRaises(ValueError):
            self._create_manifest()
        self.assertEqual(victim.read_text(encoding="utf-8"),
                         "preserve external file\n")

    def test_snapshot_rename_failure_cleans_read_only_staging(self) -> None:
        campaign_dir = self.pic_root / "manifests" / "f0_hipmpi_smoke"
        with patch(
            "control_plane_common.os.replace", side_effect=OSError("rename failed")
        ):
            with self.assertRaises(OSError):
                self._create_manifest()
        self.assertFalse(list(campaign_dir.glob(".tmp-*")))

    def test_installed_control_plane_is_versioned_and_immutable(self) -> None:
        destination = self.control_plane_dir
        inventory = json.loads(
            (destination / "inventory.json").read_text(encoding="utf-8")
        )
        self.assertEqual(destination.name, inventory["version"])
        self.assertFalse(bool(destination.stat().st_mode & 0o222))
        with self.assertRaises(ValueError):
            install(self.pic_root)

    def test_historical_control_plane_accepts_closed_predecessor_inventory(self) -> None:
        staging = self.pic_root / "control_plane" / "historical-staging"
        shutil.copytree(self.control_plane_dir, staging)
        staging.chmod(0o755)
        for path in staging.iterdir():
            path.chmod(path.stat().st_mode | 0o200)
        (staging / "terminal_recovery_handoff.py").unlink()
        names = [
            name for name in CONTROL_PLANE_FILES
            if name != "terminal_recovery_handoff.py"
        ]
        records = [{"path": name, "sha256": sha256(staging / name)} for name in names]
        digest = inventory_digest(records)
        (staging / "inventory.json").write_text(
            json.dumps(
                {"schema_version": 1, "version": digest, "files": records},
                indent=2,
                sort_keys=True,
            ) + "\n",
            encoding="utf-8",
        )
        make_tree_read_only(
            staging,
            executable_names={
                path.name for path in staging.iterdir()
                if path.suffix in {".py", ".sh"}
            },
        )
        destination = staging.with_name(digest)
        staging.rename(destination)
        inventory = verify_historical_installed_control_plane(
            destination, authorized_pic_root=self.pic_root
        )
        self.assertEqual(inventory["version"], digest)
        with self.assertRaises(ValueError):
            verify_installed_control_plane(
                destination, authorized_pic_root=self.pic_root
            )

    def test_installed_control_plane_rename_failure_cleans_staging(self) -> None:
        target = self.root / "failed-install"
        with patch("control_plane_common.os.replace", side_effect=OSError("rename failed")):
            with self.assertRaises(OSError):
                install(target)
        self.assertFalse(list((target / "control_plane").glob(".tmp-*")))

    def test_durable_tree_publication_fsyncs_entries_and_parent(self) -> None:
        staging = self.root / "durable-tree-staging"
        nested = staging / "nested"
        nested.mkdir(parents=True)
        (nested / "evidence.txt").write_text("durable\n", encoding="utf-8")
        destination = self.root / "durable-tree-final"
        fsynced_modes: list[int] = []
        real_fsync = os.fsync

        def record_fsync(descriptor: int) -> None:
            fsynced_modes.append(os.fstat(descriptor).st_mode)
            real_fsync(descriptor)

        with patch("control_plane_common.os.fsync", side_effect=record_fsync):
            durable_replace_tree(staging, destination)
        self.assertEqual(
            (destination / "nested" / "evidence.txt").read_text(encoding="utf-8"),
            "durable\n",
        )
        self.assertEqual(sum(stat.S_ISREG(mode) for mode in fsynced_modes), 1)
        self.assertGreaterEqual(sum(stat.S_ISDIR(mode) for mode in fsynced_modes), 3)

    def test_durable_parent_creation_fsyncs_created_ancestors(self) -> None:
        target = self.root / "durable-parents" / "nested"
        fsynced_modes: list[int] = []
        real_fsync = os.fsync

        def record_fsync(descriptor: int) -> None:
            fsynced_modes.append(os.fstat(descriptor).st_mode)
            real_fsync(descriptor)

        with patch("control_plane_common.os.fsync", side_effect=record_fsync):
            durable_mkdir_parents(target)
        self.assertTrue(target.is_dir())
        self.assertGreaterEqual(sum(stat.S_ISDIR(mode) for mode in fsynced_modes), 4)

    def test_durable_tree_publication_parent_fsync_failure_rolls_back(self) -> None:
        staging = self.root / "failed-durable-tree-staging"
        staging.mkdir()
        (staging / "evidence.txt").write_text("durable\n", encoding="utf-8")
        destination = self.root / "failed-durable-tree-final"
        real_fsync = os.fsync
        directory_fsyncs = 0

        def fail_publish_parent_once(descriptor: int) -> None:
            nonlocal directory_fsyncs
            if stat.S_ISDIR(os.fstat(descriptor).st_mode):
                directory_fsyncs += 1
                if directory_fsyncs == 2:
                    raise OSError("publication parent fsync failed")
            real_fsync(descriptor)

        with patch("control_plane_common.os.fsync", side_effect=fail_publish_parent_once):
            with self.assertRaises(OSError):
                durable_replace_tree(staging, destination)
        self.assertFalse(staging.exists())
        self.assertFalse(destination.exists())

    def test_durable_tree_publication_parent_swap_cannot_redirect_publication(self) -> None:
        parent = self.root / "swapped-publication-parent"
        parent.mkdir()
        staging = parent / "staging"
        staging.mkdir()
        (staging / "evidence.txt").write_text("trusted\n", encoding="utf-8")
        destination = parent / "published"
        moved_parent = self.root / "swapped-publication-parent-original"
        outside = self.root / "outside-publication-parent"
        outside.mkdir()
        outside_staging = outside / staging.name
        outside_staging.mkdir()
        (outside_staging / "evidence.txt").write_text("outside\n", encoding="utf-8")
        real_replace = os.replace
        swapped = False

        def swap_parent_then_replace(*args: object, **kwargs: object) -> None:
            nonlocal swapped
            if not swapped:
                swapped = True
                parent.rename(moved_parent)
                parent.symlink_to(outside, target_is_directory=True)
            real_replace(*args, **kwargs)

        with patch("control_plane_common.os.replace", side_effect=swap_parent_then_replace):
            with self.assertRaises((OSError, ValueError)):
                durable_replace_tree(staging, destination)
        self.assertTrue((outside_staging / "evidence.txt").is_file())
        self.assertFalse((outside / destination.name).exists())
        self.assertFalse((moved_parent / destination.name).exists())

    def test_durable_tree_publication_rollback_never_deletes_outside_sentinel(
        self,
    ) -> None:
        parent = self.root / "swapped-rollback-parent"
        parent.mkdir()
        staging = parent / "staging"
        staging.mkdir()
        (staging / "evidence.txt").write_text("trusted\n", encoding="utf-8")
        destination = parent / "published"
        moved_parent = self.root / "swapped-rollback-parent-original"
        outside = self.root / "outside-rollback-parent"
        outside.mkdir()
        outside_destination = outside / destination.name
        outside_destination.mkdir()
        sentinel = outside_destination / "sentinel.txt"
        sentinel.write_text("preserve\n", encoding="utf-8")
        parent_inode = parent.stat().st_ino
        real_fsync = os.fsync
        failed = False

        def swap_parent_then_fail(descriptor: int) -> None:
            nonlocal failed
            descriptor_stat = os.fstat(descriptor)
            if stat.S_ISDIR(descriptor_stat.st_mode) and descriptor_stat.st_ino == parent_inode:
                if not failed:
                    failed = True
                    parent.rename(moved_parent)
                    parent.symlink_to(outside, target_is_directory=True)
                    raise OSError("publication parent fsync failed")
            real_fsync(descriptor)

        with patch("control_plane_common.os.fsync", side_effect=swap_parent_then_fail):
            with self.assertRaises(OSError):
                durable_replace_tree(staging, destination)
        self.assertTrue(sentinel.is_file(), "rollback deleted the outside sentinel")
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "preserve\n")
        self.assertFalse((moved_parent / destination.name).exists())

    def test_production_installer_guard_rejects_dirty_or_untracked_source(self) -> None:
        production_root = self.root / "production-pic"
        project_home_root = self.root / "production-project-home"
        repository = install_control_plane.SCRIPT_DIR.parent
        for status in [
            " M frontier_control_plane/README.md\n",
            "?? frontier_control_plane/untracked.py\n",
        ]:
            with self.subTest(status=status):
                with patch(
                    "install_control_plane.AUTHORIZED_PIC_ROOT", production_root
                ), patch(
                    "install_control_plane.AUTHORIZED_PROJECT_HOME_ROOT",
                    project_home_root,
                ), patch(
                    "install_control_plane.subprocess.check_output",
                    side_effect=[f"{repository}\n", status],
                ), patch("install_control_plane.subprocess.run") as tracked:
                    with self.assertRaisesRegex(
                        ValueError,
                        "requires clean tracked source files",
                    ):
                        install_control_plane._require_reviewed_source_for_production(
                            production_root,
                            expected_git_commit="a" * 40,
                        )
                tracked.assert_called_once()
                self.assertEqual(
                    tracked.call_args.kwargs["env"], trusted_git_environment()
                )

    def test_production_installer_requires_exact_expected_git_commit(self) -> None:
        production_root = self.root / "production-pic"
        project_home_root = self.root / "production-project-home"
        repository = install_control_plane.SCRIPT_DIR.parent
        with patch(
            "install_control_plane.AUTHORIZED_PIC_ROOT", production_root
        ), patch(
            "install_control_plane.AUTHORIZED_PROJECT_HOME_ROOT",
            project_home_root,
        ), self.assertRaisesRegex(
            ValueError,
            "requires an expected Git commit",
        ):
            install_control_plane._require_reviewed_source_for_production(
                production_root,
                expected_git_commit=None,
            )
        with patch(
            "install_control_plane.AUTHORIZED_PIC_ROOT", production_root
        ), patch(
            "install_control_plane.AUTHORIZED_PROJECT_HOME_ROOT",
            project_home_root,
        ), patch(
            "install_control_plane.subprocess.check_output",
            side_effect=[
                f"{repository}\n",
                "",
                "b" * 40 + "\n",
            ],
        ), patch(
            "install_control_plane.subprocess.run"
        ), self.assertRaisesRegex(
            ValueError,
            "differs from expected Git commit",
        ):
            install_control_plane._require_reviewed_source_for_production(
                production_root,
                expected_git_commit="a" * 40,
            )

    def test_common_git_hash_helpers_ignore_caller_configuration(self) -> None:
        commit = "1" * 40
        tree = "2" * 40
        commit_object = (
            f"tree {tree}\n"
            "author Test Author <test@example.com> 0 +0000\n"
            "committer Test Author <test@example.com> 0 +0000\n"
            "\nmessage\n"
        ).encode("ascii")
        with patch(
            "control_plane_common.subprocess.check_output",
            side_effect=[b"archive-commit\n", f"{commit}\n".encode("ascii")],
        ) as checked:
            self.assertEqual(
                git_archive_commit_from_bytes(b"archive"),
                "archive-commit",
            )
            self.assertEqual(
                git_commit_tree_from_bytes(
                    commit_object,
                    expected_commit=commit,
                ),
                tree,
            )
        self.assertEqual(len(checked.call_args_list), 2)
        for call in checked.call_args_list:
            self.assertEqual(
                call.kwargs["env"],
                trusted_git_environment(),
            )

    def test_repository_local_fsmonitor_cannot_execute_during_git_status(self) -> None:
        from create_clean_candidate_freeze import _git

        repository = self.root / "fsmonitor-repository"
        marker = self.root / "fsmonitor-executed"
        monitor = self.root / "fsmonitor.sh"
        monitor.write_text(
            f"#!/bin/sh\n: > {marker}\n",
            encoding="utf-8",
        )
        monitor.chmod(0o755)
        subprocess.run(["git", "init", str(repository)], check=True, capture_output=True)
        subprocess.run(
            ["git", "-C", str(repository), "config", "core.fsmonitor", str(monitor)],
            check=True,
        )
        self.assertEqual(_git(repository, "status", "--porcelain"), "")
        self.assertFalse(marker.exists())

    def test_pending_marker_clear_fsyncs_parent_directory(self) -> None:
        marker = self.root / "pending-marker" / "pending_submission.json"
        marker.parent.mkdir()
        marker.write_text(
            json.dumps({"reservation_id": "reservation-1"}),
            encoding="utf-8",
        )
        with patch(
            "validate_and_reserve_frontier_job.fsync_directory"
        ) as fsync_parent:
            _clear_matching_pending_marker(marker, "reservation-1")
        self.assertFalse(marker.exists())
        fsync_parent.assert_called_once_with(marker.parent)

    def test_installed_control_plane_rejects_symlinked_inventory(self) -> None:
        inventory = self.control_plane_dir / "inventory.json"
        outside = self.root / "outside-inventory.json"
        outside.write_bytes(inventory.read_bytes())
        outside.chmod(0o444)
        self.control_plane_dir.chmod(0o755)
        inventory.unlink()
        inventory.symlink_to(outside)
        with self.assertRaises(ValueError):
            verify_installed_control_plane(
                self.control_plane_dir, authorized_pic_root=self.pic_root
            )

    def test_installed_control_plane_rejects_symlinked_profile_launcher(self) -> None:
        launcher = self.control_plane_dir / "launch_with_frontier_profile.sh"
        outside = self.root / "outside-launcher.sh"
        outside.write_bytes(launcher.read_bytes())
        outside.chmod(0o555)
        self.control_plane_dir.chmod(0o755)
        launcher.unlink()
        launcher.symlink_to(outside)
        with self.assertRaises(ValueError):
            verify_installed_control_plane(
                self.control_plane_dir, authorized_pic_root=self.pic_root
            )

    def test_installed_control_plane_rejects_extra_adjacent_entry(self) -> None:
        self.control_plane_dir.chmod(0o755)
        unexpected = self.control_plane_dir / "pathlib.py"
        unexpected.write_text("raise RuntimeError('must not import')\n", encoding="utf-8")
        unexpected.chmod(0o444)
        self.control_plane_dir.chmod(0o555)
        with self.assertRaisesRegex(ValueError, "entries differ"):
            verify_installed_control_plane(
                self.control_plane_dir, authorized_pic_root=self.pic_root
            )

    def test_installed_control_plane_rejects_open_inventory_extensions(self) -> None:
        inventory_path = self.control_plane_dir / "inventory.json"
        original = inventory_path.read_bytes()
        for mutation in ("top-level", "record"):
            with self.subTest(mutation=mutation):
                inventory = json.loads(original)
                if mutation == "top-level":
                    inventory["unexpected"] = "must reject"
                else:
                    inventory["files"][0]["unexpected"] = "must reject"
                inventory_path.chmod(0o644)
                inventory_path.write_text(json.dumps(inventory), encoding="utf-8")
                inventory_path.chmod(0o444)
                with self.assertRaises(ValueError):
                    verify_installed_control_plane(
                        self.control_plane_dir, authorized_pic_root=self.pic_root
                    )
                inventory_path.chmod(0o644)
                inventory_path.write_bytes(original)
                inventory_path.chmod(0o444)

    def test_installed_control_plane_rejects_float_inventory_schema_version(self) -> None:
        inventory_path = self.control_plane_dir / "inventory.json"
        original = inventory_path.read_bytes()
        inventory = json.loads(original)
        inventory["schema_version"] = 1.0
        inventory_path.chmod(0o644)
        inventory_path.write_text(json.dumps(inventory), encoding="utf-8")
        inventory_path.chmod(0o444)
        try:
            with self.assertRaisesRegex(ValueError, "inventory schema"):
                verify_installed_control_plane(
                    self.control_plane_dir, authorized_pic_root=self.pic_root
                )
            result = subprocess.run(
                [
                    TRUSTED_PYTHON,
                    "-I",
                    str(self.control_plane_dir / "run_control_plane.py"),
                    "promote_active_policy.py",
                    "--help",
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "Unsupported installed control-plane inventory", result.stderr
            )
        finally:
            inventory_path.chmod(0o644)
            inventory_path.write_bytes(original)
            inventory_path.chmod(0o444)

    def test_captured_runner_imports_verified_sibling_bytes(self) -> None:
        import run_control_plane

        poison_root = self.root / "poison-import"
        poison_root.mkdir()
        marker = self.root / "captured-import-marker"
        poison_root.joinpath("sibling.py").write_text(
            "VALUE = 'mutable lexical import executed'\n", encoding="utf-8"
        )
        sources = {
            "entry.py": (
                "from pathlib import Path\n"
                "import sibling\n"
                f"Path({str(marker)!r}).write_text(sibling.VALUE, encoding='utf-8')\n"
            ).encode("utf-8"),
            "sibling.py": b"VALUE = 'captured verified bytes executed'\n",
        }
        sys.modules.pop("sibling", None)
        self.addCleanup(sys.modules.pop, "sibling", None)
        with patch.object(sys, "path", [str(poison_root), *sys.path]):
            run_control_plane._execute_captured(poison_root, "entry.py", sources)
        self.assertEqual(
            marker.read_text(encoding="utf-8"), "captured verified bytes executed"
        )

    def test_task_local_verifier_executes_pinned_input_deck_descriptor(self) -> None:
        task_root = self.root / "task-local-exec"
        task_root.mkdir()
        symbols = {
            "amdhip64": "test_amdhip64",
            "mpi_amd": "test_mpi_amd",
            "mpi_gtl_hsa": "test_mpi_gtl_hsa",
        }
        for library, symbol in symbols.items():
            source = task_root / f"{library}.c"
            source.write_text(f"void {symbol}(void) {{}}\n", encoding="utf-8")
            subprocess.run(
                [
                    "/usr/bin/cc",
                    "-shared",
                    "-fPIC",
                    "-o",
                    str(task_root / f"lib{library}.so"),
                    str(source),
                ],
                check=True,
            )
        executable = task_root / "cat-wrapper"
        source = task_root / "cat-wrapper.c"
        source.write_text(
            "#include <unistd.h>\n"
            + "".join(f"void {symbol}(void);\n" for symbol in symbols.values())
            + "int main(int argc, char **argv) {\n"
            + "".join(f"  {symbol}();\n" for symbol in symbols.values())
            + '  execl("/usr/bin/cat", "cat", argv[1], (char *)0);\n'
            + "  return 1;\n"
            + "}\n",
            encoding="utf-8",
        )
        subprocess.run(
            [
                "/usr/bin/cc",
                "-o",
                str(executable),
                str(source),
                "-L",
                str(task_root),
                "-Wl,-rpath,$ORIGIN",
                "-Wl,--no-as-needed",
                "-lamdhip64",
                "-lmpi_amd",
                "-lmpi_gtl_hsa",
            ],
            check=True,
        )
        executable.chmod(0o555)
        deck = task_root / "input.athinput"
        deck.write_text("verified task-local input\n", encoding="utf-8")
        deck.chmod(0o444)
        result = subprocess.run(
            [
                TRUSTED_PYTHON,
                "-I",
                "-c",
                _TASK_LOCAL_EXEC,
                str(task_root),
                str(executable),
                sha256(executable),
                str(deck),
                sha256(deck),
                "__PIC_INPUT_DECK_FD__",
            ],
            check=True,
            text=True,
            capture_output=True,
            env={
                **os.environ,
                "LD_LIBRARY_PATH": str(task_root),
                "SLURM_PROCID": "0",
                "ROCR_VISIBLE_DEVICES": "0",
            },
        )
        self.assertRegex(
            result.stdout,
            r"^PIC trusted GPU launch: rank=0 host=\S+ ROCR_VISIBLE_DEVICES=0 "
            r"linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa\n",
        )
        self.assertTrue(result.stdout.endswith("verified task-local input\n"))

    def test_task_local_verifier_rejects_gpu_library_prefix_lookalike(self) -> None:
        task_root = self.root / "task-local-fake-prefix"
        task_root.mkdir()
        symbols = {
            "amdhip64evil": "test_amdhip64evil",
            "mpi_amd": "test_mpi_amd",
            "mpi_gtl_hsa": "test_mpi_gtl_hsa",
        }
        for library, symbol in symbols.items():
            source = task_root / f"{library}.c"
            source.write_text(f"void {symbol}(void) {{}}\n", encoding="utf-8")
            subprocess.run(
                [
                    "/usr/bin/cc",
                    "-shared",
                    "-fPIC",
                    "-o",
                    str(task_root / f"lib{library}.so"),
                    str(source),
                ],
                check=True,
            )
        executable = task_root / "cat-wrapper"
        source = task_root / "cat-wrapper.c"
        source.write_text(
            "#include <unistd.h>\n"
            + "".join(f"void {symbol}(void);\n" for symbol in symbols.values())
            + "int main(int argc, char **argv) {\n"
            + "".join(f"  {symbol}();\n" for symbol in symbols.values())
            + '  execl("/usr/bin/cat", "cat", argv[1], (char *)0);\n'
            + "  return 1;\n"
            + "}\n",
            encoding="utf-8",
        )
        subprocess.run(
            [
                "/usr/bin/cc",
                "-o",
                str(executable),
                str(source),
                "-L",
                str(task_root),
                "-Wl,-rpath,$ORIGIN",
                "-Wl,--no-as-needed",
                "-lamdhip64evil",
                "-lmpi_amd",
                "-lmpi_gtl_hsa",
            ],
            check=True,
        )
        executable.chmod(0o555)
        deck = task_root / "input.athinput"
        deck.write_text("verified task-local input\n", encoding="utf-8")
        deck.chmod(0o444)
        result = subprocess.run(
            [
                TRUSTED_PYTHON,
                "-I",
                "-c",
                _TASK_LOCAL_EXEC,
                str(task_root),
                str(executable),
                sha256(executable),
                str(deck),
                sha256(deck),
                "__PIC_INPUT_DECK_FD__",
            ],
            text=True,
            capture_output=True,
            env={
                **os.environ,
                "LD_LIBRARY_PATH": str(task_root),
                "SLURM_PROCID": "0",
                "ROCR_VISIBLE_DEVICES": "0",
            },
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not linked against libamdhip64", result.stderr)

    def test_promoted_policy_anchor_is_mirrored_and_read_only(self) -> None:
        policy = self.pic_root / "policy" / "storage_policy.json"
        mirror_policy = self.project_home_root / "policy" / "storage_policy.json"
        promotion = self.pic_root / "policy" / "active_promotion.json"
        mirror_promotion = self.project_home_root / "policy" / "active_promotion.json"
        self.assertEqual(policy.read_bytes(), mirror_policy.read_bytes())
        self.assertEqual(promotion.read_bytes(), mirror_promotion.read_bytes())
        promotion_value = json.loads(promotion.read_text(encoding="utf-8"))
        self.assertEqual(promotion_value["schema_version"], 2)
        self.assertEqual(
            str(uuid.UUID(str(promotion_value["promotion_id"]))),
            promotion_value["promotion_id"],
        )
        for path in [policy, mirror_policy, promotion, mirror_promotion]:
            self.assertFalse(bool(path.stat().st_mode & 0o222))

    def test_active_policy_reader_fails_closed_during_promotion_transaction(
        self,
    ) -> None:
        marker = self.pic_root / "policy" / ".active_promotion_transaction.json"
        marker.write_text("{}\n", encoding="utf-8")
        marker.chmod(0o400)
        with self.assertRaisesRegex(ValueError, "requires locked recovery"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                allow_pending_genesis=True,
            )

    def test_active_policy_reader_fails_closed_on_reserved_rollback_anchor_entries(
        self,
    ) -> None:
        for root in [self.pic_root, self.project_home_root]:
            for prefix in control_plane_common.ACTIVE_PROMOTION_ROLLBACK_ANCHOR_PREFIXES:
                for kind in ["regular", "symlink", "directory"]:
                    with self.subTest(root=root, prefix=prefix, kind=kind):
                        path = root / "policy" / f"{prefix}malformed"
                        if kind == "regular":
                            path.write_text("ambiguous recovery evidence\n", encoding="utf-8")
                        elif kind == "symlink":
                            path.symlink_to("missing-rollback-anchor")
                        else:
                            path.mkdir()
                        try:
                            with self.assertRaisesRegex(
                                ValueError, "requires locked recovery"
                            ):
                                require_storage_policy_unlock_snapshot(
                                    control_plane_version=self.control_plane_version,
                                    authorized_pic_root=self.pic_root,
                                    authorized_project_home_root=self.project_home_root,
                                    allow_pending_genesis=True,
                                )
                        finally:
                            if kind == "directory":
                                path.rmdir()
                            else:
                                path.unlink()

    def test_policy_promotion_transaction_handles_each_publication_failure(
        self,
    ) -> None:
        anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        real_bytes = promote_active_policy.atomic_write_bytes_at
        real_json = promote_active_policy.atomic_write_json_at

        for fail_at in range(1, 9):
            with self.subTest(fail_at=fail_at):
                self._write_policy()
                before = {path: path.read_bytes() for path in anchors}
                calls = 0

                def maybe_fail(
                    implementation: Callable[..., None],
                    *args: object,
                    **kwargs: object,
                ) -> None:
                    nonlocal calls
                    calls += 1
                    if calls == fail_at:
                        raise RuntimeError(f"injected publication failure {fail_at}")
                    implementation(*args, **kwargs)

                with patch(
                    "promote_active_policy.atomic_write_bytes_at",
                    side_effect=lambda *args, **kwargs: maybe_fail(
                        real_bytes, *args, **kwargs
                    ),
                ), patch(
                    "promote_active_policy.atomic_write_json_at",
                    side_effect=lambda *args, **kwargs: maybe_fail(
                        real_json, *args, **kwargs
                    ),
                ):
                    if fail_at <= 6:
                        with self.assertRaisesRegex(
                            RuntimeError, "injected publication failure"
                        ):
                            self._promote_policy()
                    else:
                        self._promote_policy()
                current = {path: path.read_bytes() for path in anchors}
                if fail_at <= 6:
                    self.assertEqual(current, before)
                else:
                    self.assertNotEqual(current, before)
                    self.assertEqual(anchors[0].read_bytes(), self.policy.read_bytes())
                    require_storage_policy_unlock_snapshot(
                        control_plane_version=self.control_plane_version,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                        allow_pending_genesis=True,
                    )
                self.assertFalse(
                    (self.pic_root / "policy" / ".active_promotion_transaction.json").exists()
                )
                self.assertFalse(
                    (
                        self.project_home_root
                        / "policy"
                        / ".active_promotion_transaction.json"
                    ).exists()
                )
                for root in [self.pic_root, self.project_home_root]:
                    self.assertEqual(
                        list((root / "policy").glob(".*.transaction-rollback-*")),
                        [],
                    )

    def test_policy_promotion_transaction_cleans_each_setup_failure(self) -> None:
        anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        real_link = promote_active_policy.os.link

        for fail_at in range(1, 5):
            with self.subTest(fail_at=fail_at):
                self._write_policy()
                before = {path: path.read_bytes() for path in anchors}
                calls = 0

                def fail_link(*args: object, **kwargs: object) -> None:
                    nonlocal calls
                    calls += 1
                    if calls == fail_at:
                        raise RuntimeError(f"injected rollback-link failure {fail_at}")
                    real_link(*args, **kwargs)

                with patch(
                    "promote_active_policy.os.link", side_effect=fail_link
                ), self.assertRaisesRegex(RuntimeError, "injected rollback-link failure"):
                    self._promote_policy()
                self.assertEqual({path: path.read_bytes() for path in anchors}, before)
                for root in [self.pic_root, self.project_home_root]:
                    self.assertFalse(
                        (root / "policy" / ".active_promotion_transaction.json").exists()
                    )
                    self.assertEqual(
                        list((root / "policy").glob(".*.transaction-rollback-*")),
                        [],
                    )

    def test_policy_promotion_preserves_markerless_transaction_setup_links(
        self,
    ) -> None:
        orphan_id = str(uuid.uuid4())
        for root in [self.pic_root, self.project_home_root]:
            policy_root = root / "policy"
            for name in ["storage_policy.json", "active_promotion.json"]:
                os.link(
                    policy_root / name,
                    policy_root / f".{name}.transaction-rollback-{orphan_id}",
                )

        with self.assertRaisesRegex(ValueError, "requires locked recovery"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                allow_pending_genesis=True,
            )
        self._write_policy()
        with self.assertRaisesRegex(
            ValueError,
            "Markerless active-policy promotion rollback anchors require reviewed "
            "manual recovery",
        ):
            self._promote_policy()
        for root in [self.pic_root, self.project_home_root]:
            self.assertEqual(
                len(list((root / "policy").glob(".*.transaction-rollback-*"))),
                2,
            )

    def test_policy_promotion_preserves_markerless_complete_successor_evidence(
        self,
    ) -> None:
        anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        predecessor = {path: path.read_bytes() for path in anchors}
        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        self._promote_policy()
        successor = {path: path.read_bytes() for path in anchors}
        self.assertNotEqual(successor, predecessor)
        transaction_id = str(uuid.uuid4())
        rollback_paths = []
        for path in anchors:
            rollback_path = path.with_name(
                f".{path.name}.transaction-rollback-{transaction_id}"
            )
            rollback_path.write_bytes(predecessor[path])
            rollback_path.chmod(0o444)
            rollback_paths.append(rollback_path)

        with self.assertRaisesRegex(ValueError, "requires locked recovery"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                allow_pending_genesis=True,
            )
        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        with self.assertRaisesRegex(
            ValueError,
            "Markerless active-policy promotion rollback anchors require reviewed "
            "manual recovery",
        ):
            self._promote_policy()
        self.assertEqual({path: path.read_bytes() for path in anchors}, successor)
        self.assertTrue(all(path.exists() for path in rollback_paths))

    def test_policy_promotion_recovery_rejects_forged_null_predecessors_without_mutation(
        self,
    ) -> None:
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        baseline = {path: path.read_bytes() for _, path in anchors}

        for predecessor_state, null_indexes, error_pattern in [
            ("absent", {0, 1, 2, 3}, "explicit manual recovery"),
            ("complete", {0}, "predecessor state is malformed"),
        ]:
            with self.subTest(predecessor_state=predecessor_state):
                for _, path in anchors:
                    if path.exists():
                        path.chmod(0o600)
                    path.write_bytes(baseline[path])
                    path.chmod(0o444)
                transaction_id = str(uuid.uuid4())
                marker = {
                    "schema_version": (
                        promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION
                    ),
                    "record_type": (
                        promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE
                    ),
                    "transaction_id": transaction_id,
                    "state": "prepared",
                    "predecessor_state": predecessor_state,
                    "anchors": [
                        {
                            "root_role": root_role,
                            "name": path.name,
                            "rollback_name": (
                                f".{path.name}.transaction-rollback-{transaction_id}"
                            ),
                            "predecessor_sha256": (
                                None if index in null_indexes else sha256(path)
                            ),
                            "successor_sha256": sha256(path),
                        }
                        for index, (root_role, path) in enumerate(anchors)
                    ],
                }
                if predecessor_state == "absent":
                    for index in [1, 3]:
                        anchors[index][1].unlink()
                before = {
                    path: path.read_bytes() if path.exists() else None
                    for _, path in anchors
                }
                marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
                marker_paths = []
                for root in [self.pic_root, self.project_home_root]:
                    path = root / "policy" / ".active_promotion_transaction.json"
                    path.write_text(marker_payload, encoding="utf-8")
                    path.chmod(0o400)
                    marker_paths.append(path)
                self._write_policy()
                with self.assertRaisesRegex(ValueError, error_pattern):
                    self._promote_policy()
                self.assertEqual(
                    {
                        path: path.read_bytes() if path.exists() else None
                        for _, path in anchors
                    },
                    before,
                )
                self.assertTrue(all(path.exists() for path in marker_paths))
                for path in marker_paths:
                    path.chmod(0o600)
                    path.unlink()
        for _, path in anchors:
            if path.exists():
                path.chmod(0o600)
            path.write_bytes(baseline[path])
            path.chmod(0o444)

    def test_policy_promotion_recovery_prevalidates_all_rollback_anchors(self) -> None:
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        successor = {
            path.name: f"prepared successor {path.name}\n".encode("utf-8")
            for _, path in anchors
        }
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "prepared",
            "predecessor_state": "complete",
            "anchors": [],
        }
        for root_role, path in anchors:
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            os.link(path, path.with_name(rollback_name))
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": sha256(path),
                    "successor_sha256": hashlib.sha256(
                        successor[path.name]
                    ).hexdigest(),
                }
            )
        anchors[0][1].unlink()
        anchors[0][1].write_bytes(successor[anchors[0][1].name])
        anchors[0][1].chmod(0o444)
        before = {path: path.read_bytes() for _, path in anchors}
        corrupt = anchors[-1][1].with_name(
            f".{anchors[-1][1].name}.transaction-rollback-{transaction_id}"
        )
        corrupt.unlink()
        corrupt.write_text("corrupt rollback anchor\n", encoding="utf-8")
        corrupt.chmod(0o444)
        marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
        for root in [self.pic_root, self.project_home_root]:
            path = root / "policy" / ".active_promotion_transaction.json"
            path.write_text(marker_payload, encoding="utf-8")
            path.chmod(0o400)

        self._write_policy()
        with self.assertRaisesRegex(ValueError, "rollback anchor digest differs"):
            self._promote_policy()
        self.assertEqual({path: path.read_bytes() for _, path in anchors}, before)
        for root in [self.pic_root, self.project_home_root]:
            self.assertTrue(
                (root / "policy" / ".active_promotion_transaction.json").exists()
            )

    def test_policy_promotion_prepared_recovery_rejects_semantically_invalid_rollback(
        self,
    ) -> None:
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        before = {path: path.read_bytes() for _, path in anchors}
        invalid_predecessor = b"{}\n"
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "prepared",
            "predecessor_state": "complete",
            "anchors": [],
        }
        rollback_paths = []
        for index, (root_role, path) in enumerate(anchors):
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            rollback_path = path.with_name(rollback_name)
            rollback_path.write_bytes(invalid_predecessor)
            rollback_path.chmod(0o444)
            rollback_paths.append(rollback_path)
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": hashlib.sha256(
                        invalid_predecessor
                    ).hexdigest(),
                    "successor_sha256": sha256(path),
                }
            )
            if index != 0:
                path.unlink()
                path.write_bytes(invalid_predecessor)
                path.chmod(0o444)
        before = {path: path.read_bytes() for _, path in anchors}
        marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
        marker_paths = []
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            marker_path.write_text(marker_payload, encoding="utf-8")
            marker_path.chmod(0o400)
            marker_paths.append(marker_path)

        self._write_policy()
        with self.assertRaisesRegex(
            ValueError,
            "predecessor control-plane version is malformed",
        ):
            self._promote_policy()
        self.assertEqual({path: path.read_bytes() for _, path in anchors}, before)
        self.assertTrue(all(path.exists() for path in marker_paths))
        self.assertTrue(all(path.exists() for path in rollback_paths))
        self.assertEqual(
            {path: path.read_bytes() for path in rollback_paths},
            {path: invalid_predecessor for path in rollback_paths},
        )

    def test_policy_promotion_prepared_recovery_rejects_stale_unrelated_generation(
        self,
    ) -> None:
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        stale_predecessor = {path: path.read_bytes() for _, path in anchors}
        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        self._promote_policy()
        current = {path: path.read_bytes() for _, path in anchors}
        self.assertNotEqual(current, stale_predecessor)

        transaction_id = str(uuid.uuid4())
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "prepared",
            "predecessor_state": "complete",
            "anchors": [],
        }
        rollback_paths = []
        for index, (root_role, path) in enumerate(anchors):
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            rollback_path = path.with_name(rollback_name)
            rollback_path.write_bytes(stale_predecessor[path])
            rollback_path.chmod(0o444)
            rollback_paths.append(rollback_path)
            successor_payload = (
                b"unrelated stale-marker policy successor\n"
                if path.name == "storage_policy.json"
                else b"unrelated stale-marker promotion successor\n"
            )
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": hashlib.sha256(
                        stale_predecessor[path]
                    ).hexdigest(),
                    "successor_sha256": hashlib.sha256(
                        successor_payload
                    ).hexdigest(),
                }
            )
        marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
        marker_paths = []
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            marker_path.write_text(marker_payload, encoding="utf-8")
            marker_path.chmod(0o400)
            marker_paths.append(marker_path)

        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        with self.assertRaisesRegex(ValueError, "not one reachable publication-prefix"):
            self._promote_policy()
        self.assertEqual({path: path.read_bytes() for _, path in anchors}, current)
        self.assertTrue(all(path.exists() for path in marker_paths))
        self.assertTrue(all(path.exists() for path in rollback_paths))

    def test_policy_promotion_rollback_rejects_complete_successor(self) -> None:
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "prepared",
            "predecessor_state": "complete",
            "anchors": [],
        }
        rollback_paths = []
        for root_role, path in anchors:
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            rollback_path = path.with_name(rollback_name)
            os.link(path, rollback_path)
            rollback_paths.append(rollback_path)
            digest = sha256(path)
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": digest,
                    "successor_sha256": digest,
                }
            )
        before = {path: path.read_bytes() for _, path in anchors}
        policy_descriptor = os.open(
            self.pic_root / "policy",
            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
        )
        mirror_policy_descriptor = os.open(
            self.project_home_root / "policy",
            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
        )
        try:
            with self.assertRaisesRegex(
                ValueError,
                "Complete active-policy promotion successor cannot roll back",
            ):
                promote_active_policy._rollback_promotion_transaction(
                    marker,
                    policy_descriptor,
                    mirror_policy_descriptor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    post_publish_check=lambda: None,
                    validate_predecessor=lambda *_: None,
                )
        finally:
            os.close(mirror_policy_descriptor)
            os.close(policy_descriptor)
        self.assertEqual({path: path.read_bytes() for _, path in anchors}, before)
        self.assertTrue(all(path.exists() for path in rollback_paths))

    def test_policy_promotion_rollback_rejects_committed_marker_at_mutation_boundary(
        self,
    ) -> None:
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        successor_policy = b"{}\n"
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "prepared",
            "predecessor_state": "complete",
            "anchors": [],
        }
        rollback_paths = []
        for index, (root_role, path) in enumerate(anchors):
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            rollback_path = path.with_name(rollback_name)
            os.link(path, rollback_path)
            rollback_paths.append(rollback_path)
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": sha256(path),
                    "successor_sha256": (
                        hashlib.sha256(successor_policy).hexdigest()
                        if index in [0, 2]
                        else sha256(path)
                    ),
                }
            )
        anchors[0][1].unlink()
        anchors[0][1].write_bytes(successor_policy)
        anchors[0][1].chmod(0o444)
        before = {path: path.read_bytes() for _, path in anchors}
        marker_paths = []
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            marker_path.write_text(
                json.dumps({**marker, "state": "committed"}, indent=2, sort_keys=True)
                + "\n",
                encoding="utf-8",
            )
            marker_path.chmod(0o400)
            marker_paths.append(marker_path)

        expected_predecessor_snapshot = {
            "active_policy_sha256": marker["anchors"][0]["predecessor_sha256"],
            "active_promotion_sha256": marker["anchors"][1]["predecessor_sha256"],
        }
        policy_descriptor = os.open(
            self.pic_root / "policy",
            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
        )
        mirror_policy_descriptor = os.open(
            self.project_home_root / "policy",
            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
        )
        try:
            rollback_arguments = {
                "policy_descriptor": policy_descriptor,
                "mirror_policy_descriptor": mirror_policy_descriptor,
                "authorized_pic_root": self.pic_root,
                "authorized_project_home_root": self.project_home_root,
                "post_publish_check": lambda: None,
                "validate_predecessor": lambda *_: (
                    {},
                    expected_predecessor_snapshot,
                ),
            }
            with self.assertRaisesRegex(
                ValueError,
                "Committed active-policy promotion transaction cannot roll back",
            ):
                promote_active_policy._rollback_promotion_transaction(
                    {**marker, "state": "committed"},
                    **rollback_arguments,
                )
            with self.assertRaisesRegex(
                ValueError,
                "transaction marker changed before rollback",
            ):
                promote_active_policy._rollback_promotion_transaction(
                    marker,
                    **rollback_arguments,
                )
        finally:
            os.close(mirror_policy_descriptor)
            os.close(policy_descriptor)
        self.assertEqual({path: path.read_bytes() for _, path in anchors}, before)
        self.assertTrue(all(path.exists() for path in marker_paths))
        self.assertTrue(all(path.exists() for path in rollback_paths))

    def test_policy_promotion_complete_final_validation_failure_requires_recovery(
        self,
    ) -> None:
        anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        before = {path: path.read_bytes() for path in anchors}
        self._write_policy()
        with patch(
            "promote_active_policy.require_storage_policy_unlock_snapshot",
            side_effect=ValueError("injected final validation failure"),
        ), patch(
            "promote_active_policy._rollback_promotion_transaction",
            wraps=promote_active_policy._rollback_promotion_transaction,
        ) as rollback, self.assertRaisesRegex(
            RuntimeError,
            "Committed active-policy promotion requires locked recovery",
        ):
            self._promote_policy()
        rollback.assert_not_called()
        self.assertNotEqual({path: path.read_bytes() for path in anchors}, before)
        self.assertEqual(anchors[0].read_bytes(), self.policy.read_bytes())
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            self.assertTrue(marker_path.exists())
            self.assertEqual(
                json.loads(marker_path.read_text(encoding="utf-8"))["state"],
                "prepared",
            )
            self.assertEqual(
                len(list((root / "policy").glob(".*.transaction-rollback-*"))),
                2,
            )

    def test_policy_promotion_anchor_mutation_after_committed_markers_requires_recovery(
        self,
    ) -> None:
        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        real_atomic_write_json_at = promote_active_policy.atomic_write_json_at
        committed_markers = 0

        def mutate_after_second_committed_marker(
            parent_descriptor: int,
            name: str,
            value: object,
            **kwargs: object,
        ) -> None:
            nonlocal committed_markers
            real_atomic_write_json_at(parent_descriptor, name, value, **kwargs)
            if (
                name == ".active_promotion_transaction.json"
                and isinstance(value, dict)
                and value.get("state") == "committed"
            ):
                committed_markers += 1
                if committed_markers == 2:
                    active_policy = self.pic_root / "policy" / "storage_policy.json"
                    active_policy.chmod(0o644)
                    active_policy.write_text("{}\n", encoding="utf-8")
                    active_policy.chmod(0o444)

        with patch(
            "promote_active_policy.atomic_write_json_at",
            side_effect=mutate_after_second_committed_marker,
        ), self.assertRaisesRegex(
            RuntimeError,
            "Active-policy promotion commit state requires locked recovery",
        ):
            self._promote_policy()
        self.assertEqual(committed_markers, 2)
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            self.assertTrue(marker_path.exists())
            self.assertEqual(
                json.loads(marker_path.read_text(encoding="utf-8"))["state"],
                "committed",
            )
            self.assertEqual(
                len(list((root / "policy").glob(".*.transaction-rollback-*"))),
                2,
            )
        with self.assertRaisesRegex(ValueError, "requires locked recovery"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                allow_pending_genesis=True,
            )

    def test_policy_promotion_marker_mutation_after_committed_markers_rolls_forward(
        self,
    ) -> None:
        anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        before = {path: path.read_bytes() for path in anchors}
        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        real_atomic_write_json_at = promote_active_policy.atomic_write_json_at
        committed_markers = 0

        def replace_second_committed_marker_with_prepared(
            parent_descriptor: int,
            name: str,
            value: object,
            **kwargs: object,
        ) -> None:
            nonlocal committed_markers
            real_atomic_write_json_at(parent_descriptor, name, value, **kwargs)
            if (
                name == ".active_promotion_transaction.json"
                and isinstance(value, dict)
                and value.get("state") == "committed"
            ):
                committed_markers += 1
                if committed_markers == 2:
                    marker_path = (
                        self.project_home_root
                        / "policy"
                        / ".active_promotion_transaction.json"
                    )
                    marker_path.chmod(0o600)
                    marker_path.write_text(
                        json.dumps(
                            {**value, "state": "prepared"},
                            indent=2,
                            sort_keys=True,
                        )
                        + "\n",
                        encoding="utf-8",
                    )
                    marker_path.chmod(0o400)

        with patch(
            "promote_active_policy.atomic_write_json_at",
            side_effect=replace_second_committed_marker_with_prepared,
        ):
            self._promote_policy()
        self.assertEqual(committed_markers, 2)
        self.assertNotEqual({path: path.read_bytes() for path in anchors}, before)
        self.assertEqual(anchors[0].read_bytes(), self.policy.read_bytes())
        require_storage_policy_unlock_snapshot(
            control_plane_version=self.control_plane_version,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
            allow_pending_genesis=True,
        )
        for root in [self.pic_root, self.project_home_root]:
            self.assertFalse(
                (root / "policy" / ".active_promotion_transaction.json").exists()
            )
            self.assertEqual(
                list((root / "policy").glob(".*.transaction-rollback-*")),
                [],
            )

    def test_policy_promotion_unreadable_postcommit_marker_state_never_rolls_back(
        self,
    ) -> None:
        anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        before = {path: path.read_bytes() for path in anchors}
        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        real_read_marker = promote_active_policy._read_promotion_transaction_marker
        unreadable_committed_reads = 0

        def fail_while_both_committed(*args: object, **kwargs: object) -> object:
            nonlocal unreadable_committed_reads
            marker_paths = [
                root / "policy" / ".active_promotion_transaction.json"
                for root in [self.pic_root, self.project_home_root]
            ]
            if all(path.exists() for path in marker_paths) and all(
                json.loads(path.read_text(encoding="utf-8"))["state"] == "committed"
                for path in marker_paths
            ):
                unreadable_committed_reads += 1
                raise OSError("injected unreadable postcommit marker state")
            return real_read_marker(*args, **kwargs)

        with patch(
            "promote_active_policy._read_promotion_transaction_marker",
            side_effect=fail_while_both_committed,
        ), self.assertRaisesRegex(
            RuntimeError,
            "commit state requires locked recovery",
        ):
            self._promote_policy()
        self.assertEqual(unreadable_committed_reads, 2)
        self.assertNotEqual({path: path.read_bytes() for path in anchors}, before)
        self.assertEqual(anchors[0].read_bytes(), self.policy.read_bytes())
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            self.assertTrue(marker_path.exists())
            self.assertEqual(
                json.loads(marker_path.read_text(encoding="utf-8"))["state"],
                "committed",
            )
            self.assertEqual(
                len(list((root / "policy").glob(".*.transaction-rollback-*"))),
                2,
            )
        with self.assertRaisesRegex(ValueError, "requires locked recovery"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                allow_pending_genesis=True,
            )

    def test_policy_promotion_complete_invalid_successor_requires_locked_recovery(
        self,
    ) -> None:
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": (
                promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE
            ),
            "transaction_id": transaction_id,
            "state": "prepared",
            "predecessor_state": "complete",
            "anchors": [],
        }
        for root_role, path in anchors:
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            os.link(path, path.with_name(rollback_name))
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": sha256(path),
                    "successor_sha256": hashlib.sha256(b"{}\n").hexdigest(),
                }
            )
            path.unlink()
            path.write_text("{}\n", encoding="utf-8")
            path.chmod(0o444)
        marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
        for root in [self.pic_root, self.project_home_root]:
            path = root / "policy" / ".active_promotion_transaction.json"
            path.write_text(marker_payload, encoding="utf-8")
            path.chmod(0o400)

        self._write_policy()
        with self.assertRaisesRegex(
            ValueError,
            "successor controller is malformed",
        ):
            self._promote_policy()
        for root in [self.pic_root, self.project_home_root]:
            self.assertTrue(
                (root / "policy" / ".active_promotion_transaction.json").exists()
            )
            self.assertEqual(
                len(list((root / "policy").glob(".*.transaction-rollback-*"))),
                2,
            )
        with self.assertRaisesRegex(ValueError, "requires locked recovery"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                allow_pending_genesis=True,
            )

    def test_policy_promotion_complete_valid_prepared_or_mixed_successor_rolls_forward(
        self,
    ) -> None:
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        predecessor = {path: path.read_bytes() for _, path in anchors}
        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        self._promote_policy()
        successor = {path: path.read_bytes() for _, path in anchors}
        self.assertNotEqual(successor, predecessor)
        real_predecessor_validation = (
            promote_active_policy.require_policy_predecessor_snapshot_for_promotion
        )

        def fail_after_complete_successor_recovery(
            *args: object, **kwargs: object
        ) -> object:
            if kwargs.get("allow_active_promotion_transaction") is True:
                return real_predecessor_validation(*args, **kwargs)
            raise ValueError("injected after complete-successor recovery")

        for states in [("prepared", "prepared"), ("committed", "prepared")]:
            with self.subTest(states=states):
                transaction_id = str(uuid.uuid4())
                marker = {
                    "schema_version": (
                        promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION
                    ),
                    "record_type": (
                        promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE
                    ),
                    "transaction_id": transaction_id,
                    "state": "prepared",
                    "predecessor_state": "complete",
                    "anchors": [],
                }
                rollback_paths = []
                for root_role, path in anchors:
                    rollback_name = (
                        f".{path.name}.transaction-rollback-{transaction_id}"
                    )
                    rollback_path = path.with_name(rollback_name)
                    rollback_path.write_bytes(predecessor[path])
                    rollback_path.chmod(0o444)
                    rollback_paths.append(rollback_path)
                    marker["anchors"].append(
                        {
                            "root_role": root_role,
                            "name": path.name,
                            "rollback_name": rollback_name,
                            "predecessor_sha256": hashlib.sha256(
                                predecessor[path]
                            ).hexdigest(),
                            "successor_sha256": hashlib.sha256(
                                successor[path]
                            ).hexdigest(),
                        }
                    )
                marker_paths = []
                for state, root in zip(
                    states,
                    [self.pic_root, self.project_home_root],
                ):
                    marker_path = (
                        root / "policy" / ".active_promotion_transaction.json"
                    )
                    marker_path.write_text(
                        json.dumps(
                            {**marker, "state": state},
                            indent=2,
                            sort_keys=True,
                        )
                        + "\n",
                        encoding="utf-8",
                    )
                    marker_path.chmod(0o400)
                    marker_paths.append(marker_path)

                self._write_policy(
                    admission_smoke_overrides={"status": "closed_after_pass"}
                )
                with patch(
                    "promote_active_policy.require_policy_predecessor_snapshot_for_promotion",
                    side_effect=fail_after_complete_successor_recovery,
                ), self.assertRaisesRegex(
                    ValueError,
                    "injected after complete-successor recovery",
                ):
                    self._promote_policy()
                self.assertEqual(
                    {path: path.read_bytes() for _, path in anchors},
                    successor,
                )
                self.assertFalse(any(path.exists() for path in marker_paths))
                self.assertFalse(any(path.exists() for path in rollback_paths))

    def test_policy_promotion_recovers_reachable_partial_publication_prefix(
        self,
    ) -> None:
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        before = {path: path.read_bytes() for _, path in anchors}
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "prepared",
            "predecessor_state": "complete",
            "anchors": [],
        }
        rollback_paths = []
        successor = b"{}\n"
        for root_role, path in anchors:
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            rollback_path = path.with_name(rollback_name)
            os.link(path, rollback_path)
            rollback_paths.append(rollback_path)
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": sha256(path),
                    "successor_sha256": hashlib.sha256(successor).hexdigest(),
                }
            )
        for index in [0, 2]:
            path = anchors[index][1]
            path.unlink()
            path.write_bytes(successor)
            path.chmod(0o444)
        marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
        marker_paths = []
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            marker_path.write_text(marker_payload, encoding="utf-8")
            marker_path.chmod(0o400)
            marker_paths.append(marker_path)

        real_predecessor_validation = (
            promote_active_policy.require_policy_predecessor_snapshot_for_promotion
        )

        def fail_after_partial_recovery(*args: object, **kwargs: object) -> object:
            if kwargs.get("allow_active_promotion_transaction") is True:
                return real_predecessor_validation(*args, **kwargs)
            raise ValueError("injected after partial-prefix recovery")

        self._write_policy()
        with patch(
            "promote_active_policy.require_policy_predecessor_snapshot_for_promotion",
            side_effect=fail_after_partial_recovery,
        ), self.assertRaisesRegex(ValueError, "injected after partial-prefix recovery"):
            self._promote_policy()
        self.assertEqual({path: path.read_bytes() for _, path in anchors}, before)
        self.assertFalse(any(path.exists() for path in marker_paths))
        self.assertFalse(any(path.exists() for path in rollback_paths))

    def test_policy_promotion_interrupted_rollback_remains_reachable_and_recovers(
        self,
    ) -> None:
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        before = {path: path.read_bytes() for _, path in anchors}
        successor = b"{}\n"
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "prepared",
            "predecessor_state": "complete",
            "anchors": [],
        }
        rollback_paths = []
        for index, (root_role, path) in enumerate(anchors):
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            rollback_path = path.with_name(rollback_name)
            os.link(path, rollback_path)
            rollback_paths.append(rollback_path)
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": sha256(path),
                    "successor_sha256": hashlib.sha256(successor).hexdigest(),
                }
            )
            if index in [0, 1, 2]:
                path.unlink()
                path.write_bytes(successor)
                path.chmod(0o444)
        marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
        marker_paths = []
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            marker_path.write_text(marker_payload, encoding="utf-8")
            marker_path.chmod(0o400)
            marker_paths.append(marker_path)

        real_atomic_write = promote_active_policy.atomic_write_bytes_at
        calls = 0

        def interrupt_second_rollback_write(*args: object, **kwargs: object) -> None:
            nonlocal calls
            calls += 1
            if calls == 2:
                raise OSError("injected rollback publication interruption")
            real_atomic_write(*args, **kwargs)

        self._write_policy()
        with patch(
            "promote_active_policy.atomic_write_bytes_at",
            side_effect=interrupt_second_rollback_write,
        ), self.assertRaisesRegex(
            OSError,
            "injected rollback publication interruption",
        ):
            self._promote_policy()
        for index in [1, 3]:
            self.assertEqual(anchors[index][1].read_bytes(), before[anchors[index][1]])
        for index in [0, 2]:
            self.assertEqual(anchors[index][1].read_bytes(), successor)
        self.assertTrue(all(path.exists() for path in marker_paths))
        self.assertTrue(all(path.exists() for path in rollback_paths))

        real_predecessor_validation = (
            promote_active_policy.require_policy_predecessor_snapshot_for_promotion
        )

        def fail_after_rollback_retry(*args: object, **kwargs: object) -> object:
            if kwargs.get("allow_active_promotion_transaction") is True:
                return real_predecessor_validation(*args, **kwargs)
            raise ValueError("injected after rollback retry")

        self._write_policy()
        with patch(
            "promote_active_policy.require_policy_predecessor_snapshot_for_promotion",
            side_effect=fail_after_rollback_retry,
        ), self.assertRaisesRegex(ValueError, "injected after rollback retry"):
            self._promote_policy()
        self.assertEqual({path: path.read_bytes() for _, path in anchors}, before)
        self.assertFalse(any(path.exists() for path in marker_paths))
        self.assertFalse(any(path.exists() for path in rollback_paths))

    def test_policy_promotion_recovers_committed_marker_cleanup(self) -> None:
        transaction_id = str(uuid.uuid4())
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "committed",
            "predecessor_state": "complete",
            "anchors": [],
        }
        for root_role, path in [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]:
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            os.link(path, path.with_name(rollback_name))
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": sha256(path),
                    "successor_sha256": sha256(path),
                }
            )
        marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
        for root in [self.pic_root, self.project_home_root]:
            path = root / "policy" / ".active_promotion_transaction.json"
            path.write_text(marker_payload, encoding="utf-8")
            path.chmod(0o400)

        self._write_policy()
        self._promote_policy()
        require_storage_policy_unlock_snapshot(
            control_plane_version=self.control_plane_version,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
            allow_pending_genesis=True,
        )
        for root in [self.pic_root, self.project_home_root]:
            self.assertFalse(
                (root / "policy" / ".active_promotion_transaction.json").exists()
            )
            self.assertEqual(
                list((root / "policy").glob(".*.transaction-rollback-*")),
                [],
            )

    def test_policy_promotion_committed_recovery_requires_exact_successor_anchors(
        self,
    ) -> None:
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        baseline = {path: path.read_bytes() for _, path in anchors}

        for predecessor_state, mutation in [
            ("complete", "missing"),
            ("complete", "corrupt"),
            ("absent", "missing"),
            ("absent", "corrupt"),
        ]:
            with self.subTest(
                predecessor_state=predecessor_state,
                mutation=mutation,
            ):
                transaction_id = str(uuid.uuid4())
                marker = {
                    "schema_version": (
                        promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION
                    ),
                    "record_type": (
                        promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE
                    ),
                    "transaction_id": transaction_id,
                    "state": "committed",
                    "predecessor_state": predecessor_state,
                    "anchors": [],
                }
                rollback_paths = []
                for root_role, path in anchors:
                    rollback_name = (
                        f".{path.name}.transaction-rollback-{transaction_id}"
                    )
                    if predecessor_state == "complete":
                        rollback_path = path.with_name(rollback_name)
                        os.link(path, rollback_path)
                        rollback_paths.append(rollback_path)
                    marker["anchors"].append(
                        {
                            "root_role": root_role,
                            "name": path.name,
                            "rollback_name": rollback_name,
                            "predecessor_sha256": (
                                sha256(path)
                                if predecessor_state == "complete"
                                else None
                            ),
                            "successor_sha256": sha256(path),
                        }
                    )
                marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
                marker_paths = []
                for root in [self.pic_root, self.project_home_root]:
                    marker_path = (
                        root / "policy" / ".active_promotion_transaction.json"
                    )
                    marker_path.write_text(marker_payload, encoding="utf-8")
                    marker_path.chmod(0o400)
                    marker_paths.append(marker_path)

                corrupted = anchors[0][1]
                if mutation == "missing":
                    corrupted.unlink()
                else:
                    corrupted.chmod(0o600)
                    corrupted.write_text("corrupt committed successor\n", encoding="utf-8")
                    corrupted.chmod(0o400)
                self._write_policy()
                with self.assertRaises(
                    (FileNotFoundError, ValueError)
                ):
                    self._promote_policy()
                self.assertTrue(all(path.exists() for path in marker_paths))
                self.assertEqual(
                    [path.exists() for path in rollback_paths],
                    [True] * len(rollback_paths),
                )

                for path in marker_paths:
                    path.chmod(0o600)
                    path.unlink()
                for path in rollback_paths:
                    path.unlink()
                for _, path in anchors:
                    if path.exists():
                        path.chmod(0o600)
                        path.write_bytes(baseline[path])
                    else:
                        path.write_bytes(baseline[path])
                    path.chmod(0o444)

    def test_policy_promotion_committed_recovery_preserves_evidence_until_semantic_validation(
        self,
    ) -> None:
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        before = {path: path.read_bytes() for _, path in anchors}
        invalid_successor = b"{}\n"
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "committed",
            "predecessor_state": "complete",
            "anchors": [],
        }
        rollback_paths = []
        for index, (root_role, path) in enumerate(anchors):
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            rollback_path = path.with_name(rollback_name)
            os.link(path, rollback_path)
            rollback_paths.append(rollback_path)
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": hashlib.sha256(before[path]).hexdigest(),
                    "successor_sha256": hashlib.sha256(invalid_successor).hexdigest(),
                }
            )
            path.unlink()
            path.write_bytes(invalid_successor)
            path.chmod(0o444)
        marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
        marker_paths = []
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            marker_path.write_text(marker_payload, encoding="utf-8")
            marker_path.chmod(0o400)
            marker_paths.append(marker_path)

        self._write_policy()
        with self.assertRaisesRegex(
            ValueError,
            "successor controller is malformed",
        ):
            self._promote_policy()
        self.assertEqual(
            {path: path.read_bytes() for _, path in anchors},
            {path: invalid_successor for _, path in anchors},
        )
        self.assertTrue(all(path.exists() for path in marker_paths))
        self.assertTrue(all(path.exists() for path in rollback_paths))
        self.assertEqual(
            {path: path.read_bytes() for path in rollback_paths},
            {
                rollback_path: before[path]
                for (_, path), rollback_path in zip(anchors, rollback_paths)
            },
        )

    def test_policy_promotion_committed_recovery_rechecks_successor_after_candidate_validation(
        self,
    ) -> None:
        authorized_freeze = {
            "status": "authorized",
            "manifest_path": str(
                self.pic_root
                / "clean_candidates"
                / str(uuid.uuid4())
                / "clean_candidate_manifest.json"
            ),
            "manifest_sha256": "1" * 64,
            "build_profile_control_plane_version": self.control_plane_version,
        }
        self._write_policy(science_submission_freeze=authorized_freeze)
        passed_revalidation = {
            "status": "passed",
            "current_control_plane_version": self.control_plane_version,
            "build": {
                "receipt_control_plane_version": self.control_plane_version,
            },
        }
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            return_value=passed_revalidation,
        ):
            self._promote_policy(patch_clean_candidate_revalidation=False)
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "committed",
            "predecessor_state": "complete",
            "anchors": [],
        }
        rollback_paths = []
        for root_role, path in anchors:
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            rollback_path = path.with_name(rollback_name)
            os.link(path, rollback_path)
            rollback_paths.append(rollback_path)
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": sha256(path),
                    "successor_sha256": sha256(path),
                }
            )
        marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
        marker_paths = []
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            marker_path.write_text(marker_payload, encoding="utf-8")
            marker_path.chmod(0o400)
            marker_paths.append(marker_path)
        mutated = False

        def revalidate_then_mutate(*args: object, **kwargs: object) -> dict[str, object]:
            nonlocal mutated
            if not mutated:
                active_policy = self.pic_root / "policy" / "storage_policy.json"
                active_policy.chmod(0o644)
                active_policy.write_text("{}\n", encoding="utf-8")
                active_policy.chmod(0o444)
                mutated = True
            return passed_revalidation

        self._write_policy()
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=revalidate_then_mutate,
        ), self.assertRaisesRegex(ValueError, "successor anchor digest differs"):
            self._promote_policy(patch_clean_candidate_revalidation=False)
        self.assertTrue(mutated)
        self.assertTrue(all(path.exists() for path in marker_paths))
        self.assertTrue(all(path.exists() for path in rollback_paths))
        with self.assertRaisesRegex(ValueError, "requires locked recovery"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                allow_pending_genesis=True,
            )

    def test_policy_promotion_lone_committed_marker_requires_exact_successor(
        self,
    ) -> None:
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        baseline = {path: path.read_bytes() for _, path in anchors}

        for mutation in ["none", "missing", "corrupt"]:
            with self.subTest(mutation=mutation):
                transaction_id = str(uuid.uuid4())
                marker = {
                    "schema_version": (
                        promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION
                    ),
                    "record_type": (
                        promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE
                    ),
                    "transaction_id": transaction_id,
                    "state": "committed",
                    "predecessor_state": "complete",
                    "anchors": [],
                }
                rollback_paths = []
                for root_role, path in anchors:
                    rollback_name = (
                        f".{path.name}.transaction-rollback-{transaction_id}"
                    )
                    rollback_path = path.with_name(rollback_name)
                    os.link(path, rollback_path)
                    rollback_paths.append(rollback_path)
                    marker["anchors"].append(
                        {
                            "root_role": root_role,
                            "name": path.name,
                            "rollback_name": rollback_name,
                            "predecessor_sha256": sha256(path),
                            "successor_sha256": sha256(path),
                        }
                    )
                marker_path = (
                    self.pic_root
                    / "policy"
                    / ".active_promotion_transaction.json"
                )
                marker_path.write_text(
                    json.dumps(marker, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
                marker_path.chmod(0o400)

                mutated = anchors[0][1]
                if mutation == "missing":
                    mutated.unlink()
                elif mutation == "corrupt":
                    mutated.chmod(0o600)
                    mutated.write_text(
                        "corrupt committed successor\n",
                        encoding="utf-8",
                    )
                    mutated.chmod(0o400)

                self._write_policy()
                if mutation == "none":
                    with patch(
                        "promote_active_policy."
                        "require_policy_predecessor_snapshot_for_promotion",
                        side_effect=ValueError("injected after lone-marker cleanup"),
                    ), self.assertRaisesRegex(
                        ValueError,
                        "injected after lone-marker cleanup",
                    ):
                        self._promote_policy()
                    self.assertFalse(marker_path.exists())
                    self.assertEqual(
                        [path.exists() for path in rollback_paths],
                        [False] * len(rollback_paths),
                    )
                else:
                    with self.assertRaises((FileNotFoundError, ValueError)):
                        self._promote_policy()
                    self.assertTrue(marker_path.exists())
                    self.assertEqual(
                        [path.exists() for path in rollback_paths],
                        [True] * len(rollback_paths),
                    )

                if marker_path.exists():
                    marker_path.chmod(0o600)
                    marker_path.unlink()
                for path in rollback_paths:
                    if path.exists():
                        path.unlink()
                for _, path in anchors:
                    if path.exists():
                        path.chmod(0o600)
                        path.write_bytes(baseline[path])
                    else:
                        path.write_bytes(baseline[path])
                    path.chmod(0o444)

    def test_policy_promotion_mixed_commit_partial_state_retains_recovery_evidence(
        self,
    ) -> None:
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        before = {path: path.read_bytes() for _, path in anchors}
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "prepared",
            "predecessor_state": "complete",
            "anchors": [],
        }
        for index, (root_role, path) in enumerate(anchors):
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            os.link(path, path.with_name(rollback_name))
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": sha256(path),
                    "successor_sha256": hashlib.sha256(b"{}\n").hexdigest(),
                }
            )
            if index in [0, 2]:
                path.unlink()
                path.write_text("{}\n", encoding="utf-8")
                path.chmod(0o444)
        for root, state in [
            (self.pic_root, "committed"),
            (self.project_home_root, "prepared"),
        ]:
            path = root / "policy" / ".active_promotion_transaction.json"
            path.write_text(
                json.dumps({**marker, "state": state}, indent=2, sort_keys=True)
                + "\n",
                encoding="utf-8",
            )
            path.chmod(0o400)

        self._write_policy()
        with self.assertRaisesRegex(
            ValueError,
            "successor anchor digest differs",
        ):
            self._promote_policy()
        self.assertNotEqual({path: path.read_bytes() for _, path in anchors}, before)
        for root in [self.pic_root, self.project_home_root]:
            self.assertTrue(
                (root / "policy" / ".active_promotion_transaction.json").exists()
            )
            self.assertEqual(
                len(list((root / "policy").glob(".*.transaction-rollback-*"))),
                2,
            )

    def test_policy_promotion_prepared_rollback_interruption_remains_recoverable(
        self,
    ) -> None:
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        before = {path: path.read_bytes() for _, path in anchors}
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "prepared",
            "predecessor_state": "complete",
            "anchors": [],
        }
        rollback_paths = []
        for index, (root_role, path) in enumerate(anchors):
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            rollback_path = path.with_name(rollback_name)
            os.link(path, rollback_path)
            rollback_paths.append(rollback_path)
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": sha256(path),
                    "successor_sha256": hashlib.sha256(b"{}\n").hexdigest(),
                }
            )
            if index in [0, 2]:
                path.unlink()
                path.write_text("{}\n", encoding="utf-8")
                path.chmod(0o444)
        marker_paths = []
        for root, state in [
            (self.pic_root, "prepared"),
            (self.project_home_root, "prepared"),
        ]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            marker_path.write_text(
                json.dumps({**marker, "state": state}, indent=2, sort_keys=True)
                + "\n",
                encoding="utf-8",
            )
            marker_path.chmod(0o400)
            marker_paths.append(marker_path)

        real_unlink = promote_active_policy._unlink_if_exists_at
        real_fsync = promote_active_policy.os.fsync
        fail_next_fsync = False
        failed = False

        def unlink_first_marker_then_fail_fsync(
            parent_descriptor: int, name: str
        ) -> None:
            nonlocal fail_next_fsync
            is_first_marker = (
                name == ".active_promotion_transaction.json"
                and not fail_next_fsync
                and not failed
            )
            real_unlink(parent_descriptor, name)
            if is_first_marker:
                fail_next_fsync = True

        def fail_fsync_after_first_marker_unlink(descriptor: int) -> None:
            nonlocal fail_next_fsync, failed
            if fail_next_fsync:
                fail_next_fsync = False
                failed = True
                raise OSError("injected prepared rollback marker fsync failure")
            real_fsync(descriptor)

        self._write_policy()
        with patch(
            "promote_active_policy._unlink_if_exists_at",
            side_effect=unlink_first_marker_then_fail_fsync,
        ), patch(
            "promote_active_policy.os.fsync",
            side_effect=fail_fsync_after_first_marker_unlink,
        ), self.assertRaisesRegex(
            OSError,
            "injected prepared rollback marker fsync failure",
        ):
            self._promote_policy()
        self.assertTrue(failed)
        self.assertEqual({path: path.read_bytes() for _, path in anchors}, before)
        self.assertEqual(sum(path.exists() for path in marker_paths), 1)
        remaining_marker = next(path for path in marker_paths if path.exists())
        self.assertEqual(
            json.loads(remaining_marker.read_text(encoding="utf-8"))["state"],
            "prepared",
        )
        self.assertTrue(all(path.exists() for path in rollback_paths))

        real_predecessor_validation = (
            promote_active_policy.require_policy_predecessor_snapshot_for_promotion
        )

        def fail_after_interrupted_rollback_recovery(
            *args: object, **kwargs: object
        ) -> object:
            if kwargs.get("allow_active_promotion_transaction") is True:
                return real_predecessor_validation(*args, **kwargs)
            raise ValueError("injected after interrupted rollback recovery")

        self._write_policy()
        with patch(
            "promote_active_policy.require_policy_predecessor_snapshot_for_promotion",
            side_effect=fail_after_interrupted_rollback_recovery,
        ), self.assertRaisesRegex(
            ValueError,
            "injected after interrupted rollback recovery",
        ):
            self._promote_policy()
        self.assertEqual({path: path.read_bytes() for _, path in anchors}, before)
        self.assertFalse(any(path.exists() for path in marker_paths))
        self.assertFalse(any(path.exists() for path in rollback_paths))

    def test_policy_promotion_committed_absent_predecessor_marker_only_cleans_up(
        self,
    ) -> None:
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        before = {path: path.read_bytes() for _, path in anchors}
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "committed",
            "predecessor_state": "absent",
            "anchors": [
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": (
                        f".{path.name}.transaction-rollback-{transaction_id}"
                    ),
                    "predecessor_sha256": None,
                    "successor_sha256": sha256(path),
                }
                for root_role, path in anchors
            ],
        }
        marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
        for root in [self.pic_root, self.project_home_root]:
            path = root / "policy" / ".active_promotion_transaction.json"
            path.write_text(marker_payload, encoding="utf-8")
            path.chmod(0o400)

        self._write_policy()
        with patch(
            "promote_active_policy.require_policy_predecessor_snapshot_for_promotion",
            side_effect=ValueError("injected after committed cleanup"),
        ), self.assertRaisesRegex(ValueError, "injected after committed cleanup"):
            self._promote_policy()
        self.assertEqual({path: path.read_bytes() for _, path in anchors}, before)
        for root in [self.pic_root, self.project_home_root]:
            self.assertFalse(
                (root / "policy" / ".active_promotion_transaction.json").exists()
            )
            self.assertEqual(
                list((root / "policy").glob(".*.transaction-rollback-*")),
                [],
            )

    def test_policy_promotion_committed_rollback_anchor_cleanup_failure_retains_markers(
        self,
    ) -> None:
        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        reviewed_policy = self.policy.read_bytes()
        real_unlink = promote_active_policy._unlink_if_exists_at
        failed = False

        def fail_first_rollback_anchor_cleanup(
            parent_descriptor: int, name: str
        ) -> None:
            nonlocal failed
            if not failed and ".transaction-rollback-" in name:
                failed = True
                raise OSError("injected committed rollback-anchor cleanup failure")
            real_unlink(parent_descriptor, name)

        with patch(
            "promote_active_policy._unlink_if_exists_at",
            side_effect=fail_first_rollback_anchor_cleanup,
        ):
            self._promote_policy()
        self.assertTrue(failed)
        self.assertEqual(
            (self.pic_root / "policy" / "storage_policy.json").read_bytes(),
            reviewed_policy,
        )
        for root in [self.pic_root, self.project_home_root]:
            marker = root / "policy" / ".active_promotion_transaction.json"
            self.assertTrue(marker.exists())
            self.assertEqual(
                json.loads(marker.read_text(encoding="utf-8"))["state"],
                "committed",
            )
            self.assertNotEqual(
                list((root / "policy").glob(".*.transaction-rollback-*")),
                [],
            )
        with self.assertRaisesRegex(ValueError, "requires locked recovery"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                allow_pending_genesis=True,
            )

        real_predecessor_validation = (
            promote_active_policy.require_policy_predecessor_snapshot_for_promotion
        )

        def fail_after_committed_recovery(*args: object, **kwargs: object) -> object:
            if kwargs.get("allow_active_promotion_transaction") is True:
                return real_predecessor_validation(*args, **kwargs)
            raise ValueError("injected after committed anchor-cleanup recovery")

        with patch(
            "promote_active_policy.require_policy_predecessor_snapshot_for_promotion",
            side_effect=fail_after_committed_recovery,
        ), self.assertRaisesRegex(
            ValueError,
            "injected after committed anchor-cleanup recovery",
        ):
            self._promote_policy()
        for root in [self.pic_root, self.project_home_root]:
            self.assertFalse(
                (root / "policy" / ".active_promotion_transaction.json").exists()
            )
            self.assertEqual(
                list((root / "policy").glob(".*.transaction-rollback-*")),
                [],
            )

    def test_policy_promotion_committed_marker_cleanup_failure_reports_success(
        self,
    ) -> None:
        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        reviewed_policy = self.policy.read_bytes()
        real_unlink = promote_active_policy._unlink_if_exists_at
        failed = False

        def fail_first_committed_marker_cleanup(
            parent_descriptor: int, name: str
        ) -> None:
            nonlocal failed
            if (
                not failed
                and name == ".active_promotion_transaction.json"
            ):
                marker_path = (
                    Path("/proc")
                    / str(os.getpid())
                    / "fd"
                    / str(parent_descriptor)
                    / name
                )
                if marker_path.exists():
                    marker = json.loads(marker_path.read_text(encoding="utf-8"))
                    if marker.get("state") == "committed":
                        failed = True
                        raise OSError("injected committed marker cleanup failure")
            real_unlink(parent_descriptor, name)

        with patch(
            "promote_active_policy._unlink_if_exists_at",
            side_effect=fail_first_committed_marker_cleanup,
        ):
            self._promote_policy()
        self.assertTrue(failed)
        self.assertEqual(
            (self.pic_root / "policy" / "storage_policy.json").read_bytes(),
            reviewed_policy,
        )
        self.assertTrue(
            any(
                (root / "policy" / ".active_promotion_transaction.json").exists()
                for root in [self.pic_root, self.project_home_root]
            )
        )

        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        self._promote_policy()
        require_storage_policy_unlock_snapshot(
            control_plane_version=self.control_plane_version,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
            allow_pending_genesis=True,
        )

    def test_policy_promotion_committed_marker_fsync_failure_reports_success(
        self,
    ) -> None:
        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        reviewed_policy = self.policy.read_bytes()
        real_unlink = promote_active_policy._unlink_if_exists_at
        real_fsync = promote_active_policy.os.fsync
        fail_next_fsync = False
        failed = False

        def unlink_committed_marker_then_fail_fsync(
            parent_descriptor: int, name: str
        ) -> None:
            nonlocal fail_next_fsync
            marker_path = (
                Path("/proc")
                / str(os.getpid())
                / "fd"
                / str(parent_descriptor)
                / name
            )
            is_committed_marker = False
            if (
                name == ".active_promotion_transaction.json"
                and marker_path.exists()
            ):
                marker = json.loads(marker_path.read_text(encoding="utf-8"))
                is_committed_marker = marker.get("state") == "committed"
            real_unlink(parent_descriptor, name)
            if is_committed_marker:
                fail_next_fsync = True

        def fail_fsync_after_committed_marker_unlink(descriptor: int) -> None:
            nonlocal fail_next_fsync, failed
            if fail_next_fsync:
                fail_next_fsync = False
                failed = True
                raise OSError("injected post-unlink committed marker fsync failure")
            real_fsync(descriptor)

        with patch(
            "promote_active_policy._unlink_if_exists_at",
            side_effect=unlink_committed_marker_then_fail_fsync,
        ), patch(
            "promote_active_policy.os.fsync",
            side_effect=fail_fsync_after_committed_marker_unlink,
        ):
            self._promote_policy()
        self.assertTrue(failed)
        self.assertEqual(
            (self.pic_root / "policy" / "storage_policy.json").read_bytes(),
            reviewed_policy,
        )
        self.assertEqual(
            sum(
                (
                    root / "policy" / ".active_promotion_transaction.json"
                ).exists()
                for root in [self.pic_root, self.project_home_root]
            ),
            1,
        )

        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        self._promote_policy()
        require_storage_policy_unlock_snapshot(
            control_plane_version=self.control_plane_version,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
            allow_pending_genesis=True,
        )

    def test_policy_promotion_rejects_predecessor_swap_at_transaction_setup(
        self,
    ) -> None:
        active_promotion_paths = [
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        real_transaction = promote_active_policy._active_policy_transaction
        mutated: dict[Path, bytes] = {}

        @contextmanager
        def mutate_before_transaction(*args: object, **kwargs: object) -> Iterator[None]:
            promotion = json.loads(
                active_promotion_paths[0].read_text(encoding="utf-8")
            )
            promotion["promotion_id"] = str(uuid.uuid4())
            payload = (
                json.dumps(promotion, indent=2, sort_keys=True) + "\n"
            ).encode("utf-8")
            for path in active_promotion_paths:
                path.chmod(0o644)
                path.write_bytes(payload)
                path.chmod(0o444)
                mutated[path] = payload
            with real_transaction(*args, **kwargs):
                yield

        self._write_policy()
        with patch(
            "promote_active_policy._active_policy_transaction",
            side_effect=mutate_before_transaction,
        ), self.assertRaisesRegex(ValueError, "changed before transaction setup"):
            self._promote_policy()
        self.assertEqual(
            {path: path.read_bytes() for path in active_promotion_paths},
            mutated,
        )
        for root in [self.pic_root, self.project_home_root]:
            self.assertFalse(
                (root / "policy" / ".active_promotion_transaction.json").exists()
            )

    def test_policy_promotion_rejects_broken_pending_marker_symlink(self) -> None:
        marker = self.pic_root / "ledger" / "pending_submission.json"
        marker.symlink_to(self.root / "missing-pending-marker")
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_policy_promotion_rejects_paired_ledger_parent_symlink_transplant(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        self._reserve(manifest_path)
        ledger_parent = self.pic_root / "ledger"
        mirror_parent = self.project_home_root / "ledger"
        detached_ledger = self.root / "detached-ledger"
        detached_mirror = self.root / "detached-mirror"
        empty_ledger = self.root / "empty-ledger"
        empty_mirror = self.root / "empty-mirror"
        ledger_parent.rename(detached_ledger)
        mirror_parent.rename(detached_mirror)
        empty_ledger.mkdir()
        empty_mirror.mkdir()
        ledger_parent.symlink_to(empty_ledger, target_is_directory=True)
        mirror_parent.symlink_to(empty_mirror, target_is_directory=True)
        with self.assertRaises((NotADirectoryError, ValueError)):
            self._promote_policy()
        self.assertTrue((detached_ledger / "pending_submission.json").is_file())

    def test_policy_promotion_parent_swap_fails_without_writing_replacement(self) -> None:
        import promote_active_policy

        policy_parent = self.pic_root / "policy"
        displaced = self.pic_root / "displaced-policy"
        real_write = promote_active_policy.atomic_write_bytes_at
        swapped = False

        def swap_parent_then_write(*args: object, **kwargs: object) -> None:
            nonlocal swapped
            if not swapped:
                swapped = True
                policy_parent.rename(displaced)
                policy_parent.mkdir()
            real_write(*args, **kwargs)

        with patch(
            "promote_active_policy.atomic_write_bytes_at",
            side_effect=swap_parent_then_write,
        ):
            with self.assertRaises(ValueError):
                self._promote_policy()
        self.assertTrue(swapped)
        self.assertEqual(list(policy_parent.iterdir()), [])

    def test_policy_promotion_lock_replacement_fails_closed_while_held(self) -> None:
        lock = self.pic_root / ".promotion.lock"
        with self.assertRaisesRegex(ValueError, "lock path changed"):
            with _promotion_lock(self.pic_root):
                lock.unlink()
                lock.write_text("replacement lock\n", encoding="utf-8")

    def test_policy_promotion_lock_holds_stable_serialization_anchor(self) -> None:
        anchor = stable_serialization_anchor(self.pic_root)
        with _promotion_lock(self.pic_root):
            descriptor = os.open(anchor, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
            try:
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            finally:
                os.close(descriptor)

    def test_atomic_writer_fsyncs_parent_directory(self) -> None:
        output = self.root / "durable" / "output.json"
        modes = []
        real_fsync = os.fsync

        def record_fsync(descriptor: int) -> None:
            modes.append(os.fstat(descriptor).st_mode)
            real_fsync(descriptor)

        with patch("control_plane_common.os.fsync", side_effect=record_fsync):
            atomic_write_bytes(output, b'{"fixture": true}\n')
        self.assertEqual(output.read_bytes(), b'{"fixture": true}\n')
        self.assertTrue(any(stat.S_ISDIR(mode) for mode in modes))

    def test_atomic_replacement_parent_fsync_failure_rolls_back(self) -> None:
        output = self.root / "durable-replacement" / "output.json"
        output.parent.mkdir()
        output.write_bytes(b'{"generation": "old"}\n')
        real_fsync = os.fsync
        failed = False

        def fail_first_directory_fsync(descriptor: int) -> None:
            nonlocal failed
            if stat.S_ISDIR(os.fstat(descriptor).st_mode) and not failed:
                failed = True
                raise OSError("directory fsync failed")
            real_fsync(descriptor)

        with patch("control_plane_common.os.fsync", side_effect=fail_first_directory_fsync):
            with self.assertRaises(OSError):
                atomic_write_bytes(output, b'{"generation": "new"}\n')
        self.assertTrue(failed)
        self.assertEqual(output.read_bytes(), b'{"generation": "old"}\n')

    def test_atomic_writer_parent_swap_rolls_back_pinned_publication(self) -> None:
        parent = self.root / "swapped-atomic-parent"
        parent.mkdir()
        output = parent / "output.json"
        moved_parent = self.root / "swapped-atomic-parent-original"
        outside = self.root / "outside-atomic-parent"
        outside.mkdir()
        outside_output = outside / output.name
        outside_output.write_bytes(b'{"generation": "outside"}\n')
        parent_inode = parent.stat().st_ino
        real_fsync = os.fsync
        swapped = False

        def swap_parent_then_sync(descriptor: int) -> None:
            nonlocal swapped
            descriptor_stat = os.fstat(descriptor)
            if (
                stat.S_ISDIR(descriptor_stat.st_mode)
                and descriptor_stat.st_ino == parent_inode
                and not swapped
            ):
                swapped = True
                parent.rename(moved_parent)
                parent.symlink_to(outside, target_is_directory=True)
            real_fsync(descriptor)

        with patch("control_plane_common.os.fsync", side_effect=swap_parent_then_sync):
            with self.assertRaises((OSError, ValueError)):
                atomic_write_bytes(output, b'{"generation": "trusted"}\n')
        self.assertTrue(swapped)
        self.assertEqual(outside_output.read_bytes(), b'{"generation": "outside"}\n')
        self.assertFalse((moved_parent / output.name).exists())

    def test_policy_promotion_each_directory_fsync_failure_never_accepts_new_generation(
        self,
    ) -> None:
        baseline = self.policy.read_bytes()
        replacement = baseline + b"\n"
        active_anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        baseline_anchors = {path: path.read_bytes() for path in active_anchors}
        for failure_index in range(1, 5):
            with self.subTest(failure_index=failure_index):
                self.policy.write_bytes(replacement)
                real_fsync = os.fsync
                directory_fsyncs = 0

                def fail_selected_directory_fsync(descriptor: int) -> None:
                    nonlocal directory_fsyncs
                    if stat.S_ISDIR(os.fstat(descriptor).st_mode):
                        directory_fsyncs += 1
                        if directory_fsyncs == failure_index:
                            raise OSError("directory fsync failed")
                    real_fsync(descriptor)

                with patch(
                    "control_plane_common.os.fsync",
                    side_effect=fail_selected_directory_fsync,
                ):
                    with self.assertRaises(OSError):
                        self._promote_policy()
                try:
                    require_storage_policy_unlock_snapshot(
                        control_plane_version=self.control_plane_version,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                    )
                except ValueError:
                    pass
                else:
                    self.assertEqual(
                        (self.pic_root / "policy" / "storage_policy.json").read_bytes(),
                        baseline,
                    )
                for path, payload in baseline_anchors.items():
                    path.chmod(0o644)
                    path.write_bytes(payload)
                    path.chmod(0o444)
                self.policy.write_bytes(baseline)
                self._promote_policy()

    def test_trusted_project_home_mount_alias_preserves_lexical_ledger_path(self) -> None:
        real_parent = self.root / "real-project-home-parent"
        real_root = real_parent / "mirror"
        (real_root / "ledger").mkdir(parents=True)
        alias_parent = self.root / "project-home-alias"
        alias_parent.symlink_to(real_parent, target_is_directory=True)
        trusted_root = alias_parent / "mirror"
        mirror = trusted_root / "ledger" / "node_hours.jsonl"
        mirror.write_text("trusted mount alias\n", encoding="utf-8")
        require_ledger_paths(
            self.ledger,
            self.csv,
            self.receipts,
            mirror,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=trusted_root,
        )
        self.assertEqual(
            read_stable_regular_file_below(mirror, trusted_root / "ledger"),
            b"trusted mount alias\n",
        )
        durable_mkdir_parents(
            trusted_root / "new-ledger-parent" / "nested", root=trusted_root
        )
        self.assertTrue((real_root / "new-ledger-parent" / "nested").is_dir())
        output = trusted_root / "new-ledger-parent" / "output.json"
        atomic_write_bytes(
            output,
            b'{"trusted": true}\n',
            replace=False,
            root=trusted_root,
        )
        self.assertEqual(
            (real_root / "new-ledger-parent" / "output.json").read_bytes(),
            b'{"trusted": true}\n',
        )

    def test_durable_parent_creation_rejects_below_root_alias_without_side_effect(
        self,
    ) -> None:
        trusted_root = self.root / "trusted-root"
        trusted_root.mkdir()
        outside = self.root / "outside-root"
        outside.mkdir()
        (trusted_root / "aliased").symlink_to(outside, target_is_directory=True)
        with self.assertRaises(OSError):
            durable_mkdir_parents(
                trusted_root / "aliased" / "created", root=trusted_root
            )
        self.assertFalse((outside / "created").exists())

    def test_durable_parent_creation_rejects_path_outside_missing_trusted_root(
        self,
    ) -> None:
        missing_root = self.root / "missing-trusted-root"
        outside = self.root / "outside-missing-root" / "created"
        with self.assertRaises(ValueError):
            durable_mkdir_parents(outside, root=missing_root)
        self.assertFalse(outside.exists())

    def test_policy_promoter_rejects_duplicate_json_keys(self) -> None:
        text = self.policy.read_text(encoding="utf-8")
        self.policy.write_text(
            text.replace(
                '"schema_version": 1,',
                '"schema_version": 1, "schema_version": 1,',
                1,
            ),
            encoding="utf-8",
        )
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_policy_promoter_rejects_staged_candidate_version_mismatch(self) -> None:
        self._write_policy(staged_control_plane_candidate_version="0" * 64)
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_q011_policy_promoter_requires_active_same_controller_empty_baseline(
        self,
    ) -> None:
        policy = {
            "registered_science_slices": [
                {
                    "authorization_id": "q011-section54-pressure-ps-p0-1p00-v2",
                    "campaign": "q011_section54_pressure_ps_p0_1p00",
                }
            ]
        }
        kwargs = {
            "control_plane_version": self.control_plane_version,
            "authorized_pic_root": self.pic_root,
            "authorized_project_home_root": self.project_home_root,
            "authorized_account": "AST207",
        }
        with patch(
            "promote_active_policy.require_storage_policy_unlock_snapshot",
            return_value=({"registered_science_slices": [{}]}, {}),
        ), self.assertRaisesRegex(ValueError, "active launch-prohibited"):
            promote_active_policy._require_q011_launch_prohibited_baseline(
                policy, **kwargs
            )
        with patch(
            "promote_active_policy.require_storage_policy_unlock_snapshot",
            return_value=({"registered_science_slices": []}, {}),
        ):
            promote_active_policy._require_q011_launch_prohibited_baseline(
                policy, **kwargs
            )
        self.assertTrue(
            promote_active_policy._requires_q011_launch_prohibited_baseline(
                {
                    "registered_science_slices": [
                        {
                            "authorization_id": "renamed",
                            "campaign": "renamed",
                            "input_deck_sha256": promote_active_policy.Q011_INPUT_DECK_SHA256,
                        }
                    ]
                }
            )
        )

    def test_policy_promoter_rejects_uninstalled_candidate_lifecycle(self) -> None:
        self._write_policy(
            installed_control_plane_lifecycle=(
                "live_active_generation_successor_staged_not_installed"
            )
        )
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_reserve_attach_reconcile_and_compute_node_verify(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        verify(
            manifest_path,
            submission_id=self.submission_id,
            reservation_id=reservation_id,
            manifest_sha256=str(reservation["manifest_sha256"]),
        )
        self._attach(reservation_id)
        with patch(
            "reconcile_frontier_job._scheduler_result",
            return_value=("COMPLETED", 300, 1),
        ):
            result = reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertAlmostEqual(float(result["consumed_node_hours"]), 1.0 / 12.0)
        totals = accounting(validate_primary_chain(self.ledger))
        self.assertAlmostEqual(totals["cumulative_consumed_node_hours"], 1.0 / 12.0)
        self.assertEqual(totals["currently_reserved_node_hours"], 0.0)

    def test_trampoline_reverifies_reserved_snapshot_and_binds_executable(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        executable = record_for_role(manifest, "executable")
        calls: list[tuple[list[str], dict[str, object]]] = []

        def runner(command: list[str], **kwargs: object) -> None:
            calls.append((command, kwargs))

        self._launch(manifest_path, reservation, runner=runner)
        self.assertEqual(len(calls), 1)
        self.assertRegex(
            calls[0][0][0],
            r"^/proc/self/fd/[0-9]+/launch_with_frontier_profile[.]sh$",
        )
        self.assertEqual(calls[0][0][1], "/usr/bin/srun")
        self.assertEqual(calls[0][0][2], "--jobid=12345")
        self.assertEqual(calls[0][0][13], str(executable["path"]))
        self.assertEqual(
            calls[0][0][17:22],
            ["-t", "00:05:00", "-i", "__PIC_INPUT_DECK_FD__", "-n"],
        )
        self.assertTrue(calls[0][1]["check"])
        artifact_dir = Path(str(manifest["artifact_dir"]))
        inventory = json.loads(
            (artifact_dir / "artifact_inventory.json").read_text(encoding="utf-8")
        )
        self.assertEqual(inventory["schema_version"], 1)
        self.assertEqual(
            [record["path"] for record in inventory["files"]],
            [
                "athena-parser.environment.allowlist.txt",
                "athena_stderr.txt",
                "athena_stdout.txt",
            ],
        )
        self.assertEqual(stat.S_IMODE(artifact_dir.stat().st_mode), 0o555)
        self.assertEqual(stat.S_IMODE((artifact_dir / "analysis").stat().st_mode), 0o700)

    def test_trampoline_timeout_translation_rejects_malformed_alias_and_drift(
        self,
    ) -> None:
        canonical = json.loads((self.sources / "timeout.json").read_text(encoding="utf-8"))
        profile_sha256 = sha256(self.sources / "environment.sh")

        def translate(
            snapshot_margin: dict[str, object],
            *,
            manifest_margin: dict[str, object] | None = None,
        ) -> tuple[str, str]:
            directory = self.pic_root / "timeout-translation" / str(uuid.uuid4())
            directory.mkdir(parents=True)
            timeout_margin = directory / "timeout_margin.json"
            timeout_margin.write_text(json.dumps(snapshot_margin), encoding="utf-8")
            timeout_margin.chmod(0o444)
            return launch_trampoline._trusted_athena_timeout_arguments(
                {
                    "timeout_margin": (
                        snapshot_margin if manifest_margin is None else manifest_margin
                    ),
                    "snapshot_files": [
                        {
                            "role": "timeout-margin",
                            "path": str(timeout_margin),
                            "sha256": sha256(timeout_margin),
                        },
                        {
                            "role": "environment-profile",
                            "sha256": profile_sha256,
                        },
                    ],
                },
                root=self.pic_root,
            )

        self.assertEqual(translate(canonical), ("-t", "00:05:00"))
        self.assertEqual(
            translate(
                {
                    **canonical,
                    "athena_walltime_seconds": 3661,
                    "scheduler_walltime_seconds": 7200,
                }
            ),
            ("-t", "01:01:01"),
        )
        for now in [
            datetime.now(timezone.utc) - timedelta(hours=2),
            datetime.now(timezone.utc) + timedelta(hours=2),
        ]:
            with self.subTest(now=now):
                with patch("launch_trampoline._utc_now", return_value=now):
                    with self.assertRaisesRegex(ValueError, "stale or not yet valid"):
                        translate(canonical)
        unsafe = [
            (
                "exact integers",
                {**canonical, "athena_walltime_seconds": True},
                None,
            ),
            (
                "exact integers",
                {**canonical, "athena_walltime_seconds": "300"},
                None,
            ),
            (
                "exact integers",
                {**canonical, "scheduler_walltime_seconds": 600.0},
                None,
            ),
            (
                "below Slurm walltime",
                {**canonical, "athena_walltime_seconds": 600},
                None,
            ),
            (
                "differs from manifest",
                canonical,
                {**canonical, "athena_walltime_seconds": 301},
            ),
            (
                "differs from manifest",
                {**canonical, "athena_timeout_seconds": 300},
                None,
            ),
            (
                "environment profile",
                {**canonical, "environment_profile_sha256": "0" * 64},
                None,
            ),
            (
                "canonical RFC-3339",
                {**canonical, "measured_utc": "not-a-timestamp"},
                None,
            ),
        ]
        for message, snapshot_margin, manifest_margin in unsafe:
            with self.subTest(message=message, snapshot_margin=snapshot_margin):
                with self.assertRaisesRegex(ValueError, message):
                    translate(snapshot_margin, manifest_margin=manifest_margin)

    def test_trampoline_rejects_timeout_margin_stale_at_launch(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        calls: list[list[str]] = []
        stale_now = datetime.now(timezone.utc) + timedelta(hours=2)
        with patch("launch_trampoline._utc_now", return_value=stale_now):
            with self.assertRaisesRegex(ValueError, "stale or not yet valid"):
                self._launch(
                    manifest_path,
                    reservation,
                    runner=lambda command, **_: calls.append(command),
                )
        self.assertEqual(calls, [])

    def test_trampoline_rejects_action_level_timeout_override(self) -> None:
        for timeout_arguments in [
            [{"literal": "-t"}, {"literal": "00:00:01"}],
            [{"literal": "-t=00:00:01"}],
        ]:
            with self.subTest(timeout_arguments=timeout_arguments):
                contract = self._launch_contract()
                contract["actions"][0]["arguments"].extend(timeout_arguments)
                with self.assertRaisesRegex(ValueError, "timeout override"):
                    launch_trampoline._require_no_action_timeout_override(contract)
                with self.assertRaises(ValueError):
                    validate_launch_contract(contract)

    def test_trampoline_rejects_run_artifact_root_swap_during_launch(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        artifact_dir = Path(str(manifest["artifact_dir"]))

        def runner(_: list[str], **__: object) -> None:
            artifact_dir.rename(artifact_dir.with_name(f"{artifact_dir.name}.detached"))
            artifact_dir.mkdir()

        with self.assertRaisesRegex(ValueError, "path changed"):
            self._launch(manifest_path, reservation, runner=runner)

    def test_artifact_freeze_rejects_nested_directory_swap(self) -> None:
        root = self.root / "artifact-freeze"
        nested = root / "nested"
        detached = root / "nested.detached"
        nested.mkdir(parents=True)
        (nested / "artifact.txt").write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_fsync = os.fsync
        swapped = False

        def fsync(descriptor: int) -> None:
            nonlocal swapped
            metadata = os.fstat(descriptor)
            if not swapped:
                entry = nested.stat()
                if (metadata.st_dev, metadata.st_ino) == (entry.st_dev, entry.st_ino):
                    nested.rename(detached)
                    nested.mkdir()
                    swapped = True
            real_fsync(descriptor)

        try:
            with patch("launch_trampoline.os.fsync", side_effect=fsync):
                with self.assertRaisesRegex(ValueError, "directory changed"):
                    _freeze_artifact_tree_at(root_fd)
        finally:
            os.close(root_fd)

    def test_artifact_freeze_rejects_nested_directory_swap_before_open(self) -> None:
        root = self.root / "artifact-freeze-before-open"
        nested = root / "nested"
        detached = root / "nested.detached"
        nested.mkdir(parents=True)
        (nested / "artifact.txt").write_text("original\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_open = os.open
        swapped = False

        def open_with_replacement(*args: object, **kwargs: object) -> int:
            nonlocal swapped
            if not swapped and args[0] == "nested" and kwargs.get("dir_fd") == root_fd:
                nested.rename(detached)
                nested.mkdir()
                (nested / "artifact.txt").write_text("replacement\n", encoding="utf-8")
                swapped = True
            return real_open(*args, **kwargs)

        try:
            with patch("launch_trampoline.os.open", side_effect=open_with_replacement):
                with self.assertRaisesRegex(ValueError, "changed before freezing"):
                    _freeze_artifact_tree_at(root_fd)
        finally:
            os.close(root_fd)

    def test_artifact_freeze_rejects_regular_file_swap_before_open(self) -> None:
        root = self.root / "artifact-freeze-file-before-open"
        artifact = root / "artifact.txt"
        detached = root / "artifact.detached"
        root.mkdir()
        artifact.write_text("original\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_open = os.open
        swapped = False

        def open_with_replacement(*args: object, **kwargs: object) -> int:
            nonlocal swapped
            if not swapped and args[0] == artifact.name and kwargs.get("dir_fd") == root_fd:
                artifact.rename(detached)
                artifact.write_text("replacement\n", encoding="utf-8")
                swapped = True
            return real_open(*args, **kwargs)

        try:
            with patch("launch_trampoline.os.open", side_effect=open_with_replacement):
                with self.assertRaisesRegex(ValueError, "changed before freezing"):
                    _freeze_artifact_tree_at(root_fd)
        finally:
            os.close(root_fd)

    def test_artifact_inventory_freeze_never_reopens_created_inventory(self) -> None:
        root = self.root / "artifact-inventory-retained-descriptor"
        detached = root / "artifact_inventory.detached"
        root.mkdir()
        (root / "artifact.txt").write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_open = os.open
        swapped = False

        def open_with_replacement(*args: object, **kwargs: object) -> int:
            nonlocal swapped
            if (
                not swapped
                and args[0] == "artifact_inventory.json"
                and args[1] & os.O_ACCMODE == os.O_RDONLY
                and kwargs.get("dir_fd") == root_fd
                and root.stat().st_mode & 0o222
            ):
                inventory = root / "artifact_inventory.json"
                inventory.rename(detached)
                inventory.write_bytes(detached.read_bytes())
                swapped = True
            return real_open(*args, **kwargs)

        try:
            with patch("launch_trampoline.os.open", side_effect=open_with_replacement):
                _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)
        self.assertFalse(swapped)

    def test_artifact_inventory_freeze_rejects_replacement_after_freeze(self) -> None:
        root = self.root / "artifact-inventory-retained-identity"
        detached = root / "artifact_inventory.detached"
        root.mkdir()
        (root / "artifact.txt").write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_freeze = launch_trampoline._freeze_open_artifact_file_at
        swapped = False

        def freeze_with_replacement(
            *args: object, **kwargs: object
        ) -> dict[str, object]:
            nonlocal swapped
            record = real_freeze(*args, **kwargs)
            if not swapped and args[2] == "artifact_inventory.json":
                inventory = root / "artifact_inventory.json"
                inventory.rename(detached)
                inventory.write_bytes(detached.read_bytes())
                inventory.chmod(0o444)
                swapped = True
            return record

        try:
            with patch(
                "launch_trampoline._freeze_open_artifact_file_at",
                side_effect=freeze_with_replacement,
            ):
                with self.assertRaisesRegex(ValueError, "changed after freezing"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)
        self.assertTrue(swapped)

    def test_artifact_inventory_freeze_rejects_payload_replacement_after_freeze(
        self,
    ) -> None:
        root = self.root / "artifact-payload-retained-identity"
        artifact = root / "artifact.txt"
        detached = root.with_name(f"{root.name}.artifact-detached")
        root.mkdir()
        artifact.write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_freeze = launch_trampoline._freeze_open_artifact_file_at
        swapped = False

        def freeze_with_replacement(
            *args: object, **kwargs: object
        ) -> dict[str, object]:
            nonlocal swapped
            record = real_freeze(*args, **kwargs)
            if not swapped and args[2] == "artifact_inventory.json":
                artifact.rename(detached)
                artifact.write_bytes(detached.read_bytes())
                artifact.chmod(0o444)
                swapped = True
            return record

        try:
            with patch(
                "launch_trampoline._freeze_open_artifact_file_at",
                side_effect=freeze_with_replacement,
            ):
                with self.assertRaisesRegex(ValueError, "changed after freezing"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)
        self.assertTrue(swapped)

    def test_artifact_inventory_freeze_rejects_nested_payload_replacement_after_freeze(
        self,
    ) -> None:
        root = self.root / "artifact-nested-payload-retained-identity"
        artifact = root / "nested" / "artifact.txt"
        detached = root.with_name(f"{root.name}.nested-artifact-detached")
        artifact.parent.mkdir(parents=True)
        artifact.write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_freeze = launch_trampoline._freeze_open_artifact_file_at
        swapped = False

        def freeze_with_replacement(
            *args: object, **kwargs: object
        ) -> dict[str, object]:
            nonlocal swapped
            record = real_freeze(*args, **kwargs)
            if not swapped and args[2] == "artifact_inventory.json":
                artifact.parent.chmod(0o755)
                artifact.rename(detached)
                artifact.write_bytes(detached.read_bytes())
                artifact.chmod(0o444)
                artifact.parent.chmod(0o555)
                swapped = True
            return record

        try:
            with patch(
                "launch_trampoline._freeze_open_artifact_file_at",
                side_effect=freeze_with_replacement,
            ):
                with self.assertRaisesRegex(ValueError, "changed after freezing"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)
        self.assertTrue(swapped)

    def test_artifact_inventory_freeze_rejects_payload_replacement_after_first_final_check(
        self,
    ) -> None:
        root = self.root / "artifact-payload-retained-through-final-sweep"
        artifact = root / "artifact.txt"
        detached = root.with_name(f"{root.name}.artifact-detached")
        root.mkdir()
        artifact.write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_verify = launch_trampoline._verify_open_frozen_artifact_file_at
        swapped = False

        def verify_with_replacement(*args: object, **kwargs: object) -> None:
            nonlocal swapped
            real_verify(*args, **kwargs)
            if not swapped and args[2] == "artifact.txt":
                root.chmod(0o755)
                artifact.rename(detached)
                artifact.write_bytes(detached.read_bytes())
                artifact.chmod(0o444)
                root.chmod(0o555)
                swapped = True

        try:
            with patch(
                "launch_trampoline._verify_open_frozen_artifact_file_at",
                side_effect=verify_with_replacement,
            ):
                with self.assertRaisesRegex(ValueError, "changed after freezing"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)
        self.assertTrue(swapped)

    def test_artifact_inventory_freeze_rejects_nested_payload_replacement_after_first_final_check(
        self,
    ) -> None:
        root = self.root / "artifact-nested-payload-retained-through-final-sweep"
        artifact = root / "nested" / "artifact.txt"
        detached = root.with_name(f"{root.name}.nested-artifact-detached")
        artifact.parent.mkdir(parents=True)
        artifact.write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_verify = launch_trampoline._verify_open_frozen_artifact_file_at
        swapped = False

        def verify_with_replacement(*args: object, **kwargs: object) -> None:
            nonlocal swapped
            real_verify(*args, **kwargs)
            if not swapped and args[2] == "nested/artifact.txt":
                artifact.parent.chmod(0o755)
                artifact.rename(detached)
                artifact.write_bytes(detached.read_bytes())
                artifact.chmod(0o444)
                artifact.parent.chmod(0o555)
                swapped = True

        try:
            with patch(
                "launch_trampoline._verify_open_frozen_artifact_file_at",
                side_effect=verify_with_replacement,
            ):
                with self.assertRaisesRegex(ValueError, "changed after freezing"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)
        self.assertTrue(swapped)

    def test_artifact_inventory_freeze_rejects_inventory_replacement_during_final_tree_check(
        self,
    ) -> None:
        root = self.root / "artifact-inventory-retained-after-final-tree"
        artifact = root / "artifact.txt"
        detached = root.with_name(f"{root.name}.inventory-detached")
        root.mkdir()
        artifact.write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_verify = launch_trampoline._verify_open_frozen_artifact_file_at
        swapped = False

        def verify_with_replacement(*args: object, **kwargs: object) -> None:
            nonlocal swapped
            real_verify(*args, **kwargs)
            if not swapped and args[2] == "artifact.txt":
                inventory = root / "artifact_inventory.json"
                root.chmod(0o755)
                inventory.rename(detached)
                inventory.write_bytes(detached.read_bytes())
                inventory.chmod(0o444)
                root.chmod(0o555)
                swapped = True

        try:
            with patch(
                "launch_trampoline._verify_open_frozen_artifact_file_at",
                side_effect=verify_with_replacement,
            ):
                with self.assertRaisesRegex(ValueError, "changed after freezing"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)
        self.assertTrue(swapped)

    def test_artifact_inventory_freeze_rejects_nested_directory_transplant_before_closing_payload_check(
        self,
    ) -> None:
        root = self.root / "artifact-nested-directory-retained-after-payload-sweep"
        nested = root / "nested"
        detached = root / "nested.detached"
        nested.mkdir(parents=True)
        (nested / "artifact.txt").write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_verify = launch_trampoline._verify_open_frozen_artifact_file_at
        nested_payload_checks = 0
        swapped = False

        def verify_with_replacement(*args: object, **kwargs: object) -> None:
            nonlocal nested_payload_checks, swapped
            if args[2] == "nested/artifact.txt":
                nested_payload_checks += 1
                if nested_payload_checks == 2:
                    root.chmod(0o755)
                    nested.rename(detached)
                    nested.mkdir()
                    (nested / "artifact.txt").write_text(
                        "published replacement\n", encoding="utf-8"
                    )
                    (nested / "artifact.txt").chmod(0o444)
                    nested.chmod(0o555)
                    root.chmod(0o555)
                    swapped = True
            real_verify(*args, **kwargs)

        try:
            with patch(
                "launch_trampoline._verify_open_frozen_artifact_file_at",
                side_effect=verify_with_replacement,
            ):
                with self.assertRaisesRegex(ValueError, "directory changed after freezing"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)
        self.assertTrue(swapped)

    def test_pinned_directory_ancestry_rejects_real_component_relocation(self) -> None:
        anchor = self.root / "stable-anchor"
        pic_root = anchor / "project" / "PIC"
        leaf = pic_root / "runs" / "campaign"
        leaf.mkdir(parents=True)
        detached = anchor / "project" / "PIC.detached"
        with PinnedDirectoryAncestry(leaf, root=anchor) as ancestry:
            pic_root.rename(detached)
            pic_root.mkdir()
            (detached / "runs").rename(pic_root / "runs")
            with self.assertRaisesRegex(ValueError, "ancestry changed"):
                ancestry.require_same()

    def test_artifact_freeze_rejects_workload_descendant_replacement_after_capture(
        self,
    ) -> None:
        root = self.root / "artifact-freeze-workload-descendant"
        nested = root / "output" / "nested"
        nested.mkdir(parents=True)
        (nested / "artifact.txt").write_text("original\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        identities: dict[str, tuple[int, int]] = {}
        try:
            _capture_artifact_directory_identities_at(root_fd, identities)
            nested.rename(root / "output" / "nested.detached")
            nested.mkdir()
            (nested / "artifact.txt").write_text("replacement\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "changed before freezing"):
                _freeze_artifact_tree_at(
                    root_fd,
                    directory_identities=identities,
                )
        finally:
            os.close(root_fd)

    def test_artifact_freeze_hashes_final_read_only_bytes(self) -> None:
        root = self.root / "artifact-freeze-final-bytes"
        root.mkdir()
        artifact = root / "artifact.txt"
        artifact.write_text("old\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_fchmod = os.fchmod
        replaced = False

        def fchmod(descriptor: int, mode: int) -> None:
            nonlocal replaced
            metadata = os.fstat(descriptor)
            entry = artifact.stat()
            if (
                not replaced
                and (metadata.st_dev, metadata.st_ino) == (entry.st_dev, entry.st_ino)
            ):
                artifact.write_text("new\n", encoding="utf-8")
                replaced = True
            real_fchmod(descriptor, mode)

        try:
            with patch("launch_trampoline.os.fchmod", side_effect=fchmod):
                record = _freeze_artifact_file_at(root_fd, artifact.name, artifact.name)
        finally:
            os.close(root_fd)
        self.assertEqual(record["sha256"], hashlib.sha256(b"new\n").hexdigest())
        self.assertEqual(stat.S_IMODE(artifact.stat().st_mode), 0o444)

    def test_artifact_freeze_rejects_empty_directory(self) -> None:
        root = self.root / "artifact-freeze-empty-directory"
        (root / "empty").mkdir(parents=True)
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        try:
            with self.assertRaisesRegex(ValueError, "empty directory"):
                _freeze_artifact_tree_at(root_fd)
        finally:
            os.close(root_fd)

    def test_artifact_inventory_rejects_late_unlisted_file(self) -> None:
        root = self.root / "artifact-freeze-late-file"
        root.mkdir()
        (root / "artifact.txt").write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_freeze = launch_trampoline._freeze_open_artifact_file_at

        def freeze_with_late_file(*args: object, **kwargs: object) -> dict[str, object]:
            record = real_freeze(*args, **kwargs)
            if args[2] == "artifact_inventory.json":
                late = root / "late.txt"
                late.write_text("late\n", encoding="utf-8")
                late.chmod(0o444)
            return record

        try:
            with patch(
                "launch_trampoline._freeze_open_artifact_file_at",
                side_effect=freeze_with_late_file,
            ):
                with self.assertRaisesRegex(ValueError, "unlisted file"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)

    def test_artifact_inventory_rejects_late_nested_directory_swap(self) -> None:
        root = self.root / "artifact-freeze-late-directory"
        nested = root / "nested"
        detached = root / "nested.detached"
        nested.mkdir(parents=True)
        (nested / "artifact.txt").write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_freeze = launch_trampoline._freeze_open_artifact_file_at

        def freeze_with_late_directory(
            *args: object, **kwargs: object
        ) -> dict[str, object]:
            record = real_freeze(*args, **kwargs)
            if args[2] == "artifact_inventory.json":
                nested.rename(detached)
                nested.mkdir()
                replacement = nested / "artifact.txt"
                replacement.write_text("verified\n", encoding="utf-8")
                replacement.chmod(0o444)
                nested.chmod(0o555)
            return record

        try:
            with patch(
                "launch_trampoline._freeze_open_artifact_file_at",
                side_effect=freeze_with_late_directory,
            ):
                with self.assertRaisesRegex(ValueError, "directory changed"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)

    def test_artifact_inventory_rejects_late_analysis_directory_swap(self) -> None:
        root = self.root / "artifact-freeze-late-analysis-directory"
        detached = root / "analysis.detached"
        root.mkdir()
        (root / "artifact.txt").write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_freeze = launch_trampoline._freeze_open_artifact_file_at

        def freeze_with_late_analysis_directory(
            *args: object, **kwargs: object
        ) -> dict[str, object]:
            record = real_freeze(*args, **kwargs)
            if args[2] == "artifact_inventory.json":
                (root / "analysis").rename(detached)
                (root / "analysis").mkdir()
                (root / "analysis").chmod(0o777)
            return record

        try:
            with patch(
                "launch_trampoline._freeze_open_artifact_file_at",
                side_effect=freeze_with_late_analysis_directory,
            ):
                with self.assertRaisesRegex(ValueError, "analysis directory changed"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)

    def test_artifact_inventory_rejects_analysis_directory_publish_swap(self) -> None:
        root = self.root / "artifact-freeze-analysis-publish-swap"
        detached = root / "analysis.detached"
        root.mkdir()
        (root / "artifact.txt").write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_rename = os.rename

        def rename_with_replacement(*args: object, **kwargs: object) -> None:
            real_rename(*args, **kwargs)
            if len(args) >= 2 and args[1] == "analysis":
                (root / "analysis").rename(detached)
                (root / "analysis").mkdir(mode=0o700)

        try:
            with patch("launch_trampoline.os.rename", side_effect=rename_with_replacement):
                with self.assertRaisesRegex(ValueError, "analysis directory changed"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)

    def test_artifact_inventory_rejects_analysis_staging_swap_before_open(self) -> None:
        root = self.root / "artifact-freeze-analysis-staging-swap"
        root.mkdir()
        (root / "artifact.txt").write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_open = os.open
        swapped = False

        def open_with_replacement(*args: object, **kwargs: object) -> int:
            nonlocal swapped
            name = args[0]
            if (
                not swapped
                and isinstance(name, str)
                and name.startswith(".analysis.staging-")
                and kwargs.get("dir_fd") == root_fd
            ):
                staging = root / name
                staging.rename(root / f"{name}.detached")
                staging.mkdir(mode=0o700)
                swapped = True
            return real_open(*args, **kwargs)

        try:
            with patch("launch_trampoline.os.open", side_effect=open_with_replacement):
                with self.assertRaisesRegex(ValueError, "changed between creation and open"):
                    _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.close(root_fd)

    def test_artifact_root_rejects_swap_between_creation_and_open(self) -> None:
        pic_root = self.root / "pic-root-create-open"
        target = pic_root / "runs" / "campaign" / "submission"
        target.parent.mkdir(parents=True)
        detached = target.with_name("submission.detached")
        real_open = os.open
        swapped = False

        def open_with_replacement(*args: object, **kwargs: object) -> int:
            nonlocal swapped
            if not swapped and args[0] == target.name and target.exists():
                target.rename(detached)
                target.mkdir()
                swapped = True
            return real_open(*args, **kwargs)

        with patch("launch_trampoline.os.open", side_effect=open_with_replacement):
            with self.assertRaisesRegex(ValueError, "changed between creation and open"):
                launch_trampoline._create_artifact_directory(target, pic_root=pic_root)

    def test_artifact_root_rejects_pic_root_relocation_after_open(self) -> None:
        pic_root = self.root / "pic-root-after-open"
        target = pic_root / "runs" / "campaign" / "submission"
        target.parent.mkdir(parents=True)
        detached = self.root / "pic-root-after-open.detached"
        real_create = launch_trampoline._open_created_directory_at

        def create_with_relocation(*args: object, **kwargs: object) -> int:
            descriptor = real_create(*args, **kwargs)
            pic_root.rename(detached)
            pic_root.mkdir()
            (detached / "runs").rename(pic_root / "runs")
            return descriptor

        with patch(
            "launch_trampoline._open_created_directory_at",
            side_effect=create_with_relocation,
        ):
            with self.assertRaisesRegex(ValueError, "ancestry changed"):
                launch_trampoline._create_artifact_directory(target, pic_root=pic_root)

    def test_launch_directory_rejects_swap_between_creation_and_open(self) -> None:
        root = self.root / "launch-directory-create-open"
        root.mkdir()
        output = root / "output"
        detached = root / "output.detached"
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        real_open = os.open
        swapped = False

        def open_with_replacement(*args: object, **kwargs: object) -> int:
            nonlocal swapped
            if not swapped and args[0] == output.name and output.exists():
                output.rename(detached)
                output.mkdir()
                swapped = True
            return real_open(*args, **kwargs)

        try:
            with patch("launch_trampoline.os.open", side_effect=open_with_replacement):
                with self.assertRaisesRegex(ValueError, "changed between creation and open"):
                    launch_trampoline._mkdir_artifact_directory(
                        root_fd, root, output, {}
                    )
        finally:
            os.close(root_fd)

    def test_launch_directory_rejects_replacement_after_retention(self) -> None:
        root = self.root / "launch-directory-retained"
        root.mkdir()
        output = root / "output"
        detached = root / "output.detached"
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        identities: dict[str, tuple[int, int]] = {}
        try:
            launch_trampoline._mkdir_artifact_directory(
                root_fd, root, output, identities
            )
            output.rename(detached)
            output.mkdir()
            with self.assertRaisesRegex(ValueError, "changed during execution"):
                launch_trampoline._require_retained_artifact_directories_at(
                    root_fd, root, identities
                )
        finally:
            os.close(root_fd)

    def test_artifact_inventory_normalizes_restrictive_umask_for_analysis(self) -> None:
        root = self.root / "artifact-freeze-restrictive-umask"
        root.mkdir()
        (root / "artifact.txt").write_text("verified\n", encoding="utf-8")
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        inherited_umask = os.umask(0o777)
        try:
            _publish_frozen_artifact_inventory(root_fd, root)
        finally:
            os.umask(inherited_umask)
            os.close(root_fd)
        self.assertEqual(stat.S_IMODE((root / "analysis").stat().st_mode), 0o700)

    def test_trampoline_compute_snapshot_does_not_require_flock(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        calls: list[list[str]] = []

        def runner(command: list[str], **kwargs: object) -> None:
            calls.append(command)

        self._launch(
            manifest_path,
            reservation,
            runner=runner,
            flock_error=OSError(524, "Unknown error 524"),
        )
        self.assertEqual(len(calls), 1)

    def test_trampoline_rejects_wrong_executable_binding(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        job_script = record_for_role(manifest, "job-script")
        with patch.dict(
            os.environ,
            {
                "PIC_MANIFEST_SHA256": str(reservation["manifest_sha256"]),
                "PIC_RESERVATION_ID": str(reservation["reservation_id"]),
                "PIC_SUBMISSION_ID": self.submission_id,
            },
            clear=True,
        ):
            with self.assertRaises(ValueError):
                launch(
                    manifest_path=manifest_path,
                    manifest_sha256=str(reservation["manifest_sha256"]),
                    job_script_sha256=str(job_script["sha256"]),
                    executable_sha256="0" * 64,
                    reservation_id=str(reservation["reservation_id"]),
                    submission_id=self.submission_id,
                    ledger_jsonl=self.ledger,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )

    def test_trampoline_rejects_project_home_mirror_as_primary_ledger(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        self._attach(str(reservation["reservation_id"]))
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        job_script = record_for_role(manifest, "job-script")
        executable = record_for_role(manifest, "executable")
        with patch.dict(
            os.environ,
            {
                "PIC_MANIFEST_SHA256": str(reservation["manifest_sha256"]),
                "PIC_RESERVATION_ID": str(reservation["reservation_id"]),
                "PIC_SUBMISSION_ID": self.submission_id,
            },
            clear=True,
        ):
            with self.assertRaises(ValueError):
                launch(
                    manifest_path=manifest_path,
                    manifest_sha256=str(reservation["manifest_sha256"]),
                    job_script_sha256=str(job_script["sha256"]),
                    executable_sha256=str(executable["sha256"]),
                    reservation_id=str(reservation["reservation_id"]),
                    submission_id=self.submission_id,
                    ledger_jsonl=self.mirror,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    runner=lambda *_args, **_kwargs: None,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )

    def test_trampoline_ignores_job_template_shell_body(self) -> None:
        self._write(
            "job.sh",
            "#!/bin/bash\n#SBATCH -A AST207\n#SBATCH -p batch\n#SBATCH -q debug\n"
            f"#SBATCH -o {self.pic_root}/logs/slurm/%x.%j.log\n"
            "#SBATCH -N 1\n#SBATCH -t 00:10:00\n/bin/true\n",
        )
        self._write_policy()
        self._promote_policy()
        self._write_config()
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        commands: list[list[str]] = []

        def runner(command: list[str], **_: object) -> None:
            commands.append(command)

        self._launch(manifest_path, reservation, runner=runner)
        self.assertEqual(len(commands), 1)
        self.assertRegex(
            commands[0][0],
            r"^/proc/self/fd/[0-9]+/launch_with_frontier_profile[.]sh$",
        )
        self.assertEqual(commands[0][1], "/usr/bin/srun")
        self.assertEqual(commands[0][2], "--jobid=12345")
        self.assertNotIn("/bin/bash", commands[0])
        self.assertNotIn("/bin/true", commands[0])

    def test_profile_wrapper_closes_runtime_allowlist_fd_before_exec(self) -> None:
        wrapper_dir = self.root / "profile-wrapper"
        wrapper_dir.mkdir()
        wrapper = wrapper_dir / "launch_with_frontier_profile.sh"
        wrapper.write_text(
            Path(__file__).with_name("launch_with_frontier_profile.sh").read_text(
                encoding="utf-8"
            ),
            encoding="utf-8",
        )
        wrapper.chmod(0o755)
        (wrapper_dir / "frontier_pic_environment.sh").write_text(
            "record_pic_environment() { printf 'ALLOWLISTED=1\\n'; }\n",
            encoding="utf-8",
        )
        allowlist = wrapper_dir / "environment.allowlist.txt"
        with allowlist.open("wb") as stream:
            descriptor = stream.fileno()
            directory_descriptor = os.open(wrapper_dir, os.O_RDONLY | os.O_DIRECTORY)
            control_plane_descriptor = os.open(wrapper_dir, os.O_RDONLY | os.O_DIRECTORY)
            try:
                environment = dict(os.environ)
                environment["PIC_CONTROL_PLANE_DIR_FD"] = str(control_plane_descriptor)
                environment["PIC_RUNTIME_ALLOWLIST_FD"] = str(descriptor)
                environment["PIC_RUNTIME_ALLOWLIST_DIR_FD"] = str(directory_descriptor)
                environment["BASH_FUNC_module%%"] = "() {  :\n}"
                result = subprocess.run(
                    [
                        str(wrapper),
                        "/bin/bash",
                        "-c",
                        (
                            'test -z "${PIC_RUNTIME_ALLOWLIST_FD+x}" '
                            '&& test -z "${PIC_RUNTIME_ALLOWLIST_DIR_FD+x}" '
                            '&& ! printf "FORGED_CHILD_WRITE\\n" > "$1"'
                        ),
                        "bash",
                        str(allowlist),
                    ],
                    env=environment,
                    pass_fds=(descriptor, directory_descriptor, control_plane_descriptor),
                    check=False,
                )
            finally:
                os.close(directory_descriptor)
                os.close(control_plane_descriptor)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(allowlist.read_text(encoding="utf-8"), "ALLOWLISTED=1\n")
        self.assertEqual(stat.S_IMODE(allowlist.stat().st_mode), 0o400)

    def test_profile_wrapper_does_not_source_user_bashrc(self) -> None:
        wrapper = Path(__file__).with_name("launch_with_frontier_profile.sh")
        wrapper_text = wrapper.read_text(encoding="utf-8")
        self.assertNotIn("source /etc/profile", wrapper_text.splitlines())
        self.assertIn("source /opt/cray/pe/lmod/lmod/init/profile", wrapper_text)
        self.assertIn("export HOME=/", wrapper_text)

        wrapper_dir = self.root / "profile-wrapper"
        wrapper_dir.mkdir()
        copied_wrapper = wrapper_dir / wrapper.name
        copied_wrapper.write_text(wrapper_text, encoding="utf-8")
        copied_wrapper.chmod(0o755)
        (wrapper_dir / "frontier_pic_environment.sh").write_text(
            "record_pic_environment() { printf 'ALLOWLISTED=1\\n'; }\n",
            encoding="utf-8",
        )
        home = self.root / "forged-home"
        home.mkdir()
        marker = self.root / "user-bashrc-sourced"
        (home / ".bashrc").write_text(
            f": > {marker}\n",
            encoding="utf-8",
        )
        allowlist = wrapper_dir / "environment.allowlist.txt"
        with allowlist.open("wb") as stream:
            directory_descriptor = os.open(wrapper_dir, os.O_RDONLY | os.O_DIRECTORY)
            control_plane_descriptor = os.open(wrapper_dir, os.O_RDONLY | os.O_DIRECTORY)
            try:
                environment = {
                    **os.environ,
                    "HOME": str(home),
                    "PIC_CONTROL_PLANE_DIR_FD": str(control_plane_descriptor),
                    "PIC_RUNTIME_ALLOWLIST_FD": str(stream.fileno()),
                    "PIC_RUNTIME_ALLOWLIST_DIR_FD": str(directory_descriptor),
                }
                subprocess.run(
                    [str(copied_wrapper), "/bin/true"],
                    env=environment,
                    pass_fds=(
                        stream.fileno(),
                        directory_descriptor,
                        control_plane_descriptor,
                    ),
                    check=True,
                )
            finally:
                os.close(directory_descriptor)
                os.close(control_plane_descriptor)
        self.assertFalse(marker.exists())
        self.assertEqual(allowlist.read_text(encoding="utf-8"), "ALLOWLISTED=1\n")

    def test_production_profile_wrapper_emits_canonical_allowlist_from_stripped_environment(
        self,
    ) -> None:
        from launch_trampoline import _profile_environment

        repo_root = str(Path(__file__).resolve().parents[3])
        if repo_root not in sys.path:
            sys.path.insert(0, repo_root)
            self.addCleanup(sys.path.remove, repo_root)
        from tst.publication.pic_qualification_manifest import (
            _validate_environment_allowlist,
        )

        wrapper = Path(__file__).with_name("launch_with_frontier_profile.sh")
        control_plane_dir = wrapper.parent
        inherited = {**os.environ, **_profile_environment()}
        inherited.pop("BASH_ENV", None)
        inherited.pop("ENV", None)
        inherited.pop("OMP_NUM_THREADS", None)
        variants = {
            "stripped": _profile_environment(),
            "inherited": inherited,
            "poisoned": {
                **inherited,
                "MODULEPATH": "/tmp/caller-controlled-modulepath",
            },
        }
        payloads: list[bytes] = []
        for name, base_environment in variants.items():
            allowlist = self.root / f"{name}.environment.allowlist.txt"
            descriptor = os.open(allowlist, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            directory_descriptor = os.open(self.root, os.O_RDONLY | os.O_DIRECTORY)
            control_plane_descriptor = os.open(
                control_plane_dir, os.O_RDONLY | os.O_DIRECTORY
            )
            try:
                environment = {
                    **base_environment,
                    "PIC_CONTROL_PLANE_DIR_FD": str(control_plane_descriptor),
                    "PIC_RUNTIME_ALLOWLIST_FD": str(descriptor),
                    "PIC_RUNTIME_ALLOWLIST_DIR_FD": str(directory_descriptor),
                }
                subprocess.run(
                    ["/bin/bash", str(wrapper), "/bin/true"],
                    env=environment,
                    pass_fds=(
                        descriptor,
                        directory_descriptor,
                        control_plane_descriptor,
                    ),
                    check=True,
                )
            finally:
                os.close(descriptor)
                os.close(directory_descriptor)
                os.close(control_plane_descriptor)
            payload = allowlist.read_bytes()
            _validate_environment_allowlist(payload, require_frontier_values=True)
            payloads.append(payload)
        self.assertEqual(payloads, [payloads[0]] * len(payloads))

    def test_production_module_measurement_accepts_profile_and_rejects_late_mutation(
        self,
    ) -> None:
        profile = Path(__file__).with_name("frontier_pic_environment.sh")
        command = r"""
source /opt/cray/pe/lmod/lmod/init/profile
source "$1"
PYTHONPATH="$2" /opt/cray/pe/python/3.11.7/bin/python3 - <<'PY'
import os
from control_plane_common import measured_production_module_list_bytes
measured_production_module_list_bytes()
print("writer-precondition-canonical-ok")
os.environ["MODULEPATH"] += ":/tmp/post-activation-forged"
try:
    measured_production_module_list_bytes()
except ValueError:
    print("writer-precondition-post-mutation-rejected")
else:
    raise SystemExit("post-activation MODULEPATH mutation accepted")
PY
"""
        environment = {**os.environ, "MODULEPATH": "/tmp/caller-controlled-modulepath"}
        environment.pop("BASH_ENV", None)
        environment.pop("ENV", None)
        result = subprocess.run(
            ["/bin/bash", "-c", command, "bash", str(profile), str(profile.parent)],
            check=True,
            capture_output=True,
            text=True,
            env=environment,
        )
        self.assertEqual(
            result.stdout.splitlines(),
            [
                "writer-precondition-canonical-ok",
                "writer-precondition-post-mutation-rejected",
            ],
        )

    def test_trampoline_strips_bash_startup_hooks_before_profile_wrapper(self) -> None:
        hook = self.root / "bash-env-hook.sh"
        hook.write_text(
            'eval "exec 9>&${PIC_RUNTIME_ALLOWLIST_FD}"\n',
            encoding="utf-8",
        )
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        seen_environment: dict[str, str] = {}

        def runner(_: list[str], **kwargs: object) -> None:
            environment = dict(kwargs["env"])
            seen_environment.update(environment)
            result = subprocess.run(
                [
                    "/bin/bash",
                    "-c",
                    (
                        'printf "ALLOWLISTED=1\\n" >&"$PIC_RUNTIME_ALLOWLIST_FD"; '
                        'eval "exec ${PIC_RUNTIME_ALLOWLIST_FD}>&-"; '
                        'eval "exec ${PIC_RUNTIME_ALLOWLIST_DIR_FD}>&-"; '
                        'eval "exec ${PIC_CONTROL_PLANE_DIR_FD}>&-"; '
                        "unset PIC_CONTROL_PLANE_DIR_FD PIC_RUNTIME_ALLOWLIST_FD "
                        "PIC_RUNTIME_ALLOWLIST_DIR_FD; "
                        "exec /bin/bash -c '! printf \"FORGED_CHILD_WRITE\\\\n\" >&9'"
                    ),
                ],
                env=environment,
                pass_fds=kwargs["pass_fds"],
                check=False,
            )
            self.assertEqual(result.returncode, 0)

        injected = {
            "BASH_ENV": str(hook),
            "ENV": str(hook),
            "BASH_FUNC_injected%%": "() { :; }",
            "CDPATH": str(self.root),
            "LD_PRELOAD": str(hook),
            "LD_LIBRARY_PATH": str(self.root),
            "PYTHONPATH": str(self.root),
            "HSA_XNACK": "1",
            "MPICH_OFI_NIC_POLICY": "GPU",
            "FI_MR_CACHE_MONITOR": "forged",
        }
        self._launch(
            manifest_path,
            reservation,
            runner=runner,
            environment_overrides=injected,
        )
        for name in injected:
            self.assertNotIn(name, seen_environment)
        self.assertEqual(
            set(seen_environment),
            {
                "LC_ALL",
                "PATH",
                "PIC_FRONTIER_PROFILE",
                "PIC_RUNTIME_ALLOWLIST_FD",
                "PIC_RUNTIME_ALLOWLIST_DIR_FD",
                "PIC_CONTROL_PLANE_DIR_FD",
            },
        )
        self.assertEqual(seen_environment["LC_ALL"], "C")
        self.assertEqual(seen_environment["PATH"], "/usr/bin:/bin")
        self.assertEqual(
            seen_environment["PIC_FRONTIER_PROFILE"], "frontier_minimum_supported"
        )
        allowlist = (
            self.pic_root
            / "runs"
            / "f0_hipmpi_smoke"
            / self.submission_id
            / "athena-parser.environment.allowlist.txt"
        )
        self.assertEqual(allowlist.read_text(encoding="utf-8"), "ALLOWLISTED=1\n")

    def test_trampoline_requires_live_slurm_job_id(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        with self.assertRaises(ValueError):
            self._launch(
                manifest_path,
                reservation,
                runner=lambda *_args, **_kwargs: None,
                environment_overrides={"SLURM_JOB_ID": ""},
            )

    def test_trampoline_rejects_existing_run_artifact_directory(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        artifact_dir = (
            self.pic_root / "runs" / "f0_hipmpi_smoke" / self.submission_id
        )
        artifact_dir.mkdir(parents=True)
        (artifact_dir / "stale.txt").write_text("must not be reused\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            self._launch(
                manifest_path,
                reservation,
                runner=lambda *_args, **_kwargs: None,
            )

    def test_reservation_rejects_existing_run_artifact_directory(self) -> None:
        manifest_path = self._create_manifest()
        artifact_dir = (
            self.pic_root / "runs" / "f0_hipmpi_smoke" / self.submission_id
        )
        artifact_dir.mkdir(parents=True)
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_reservation_rejects_unbound_run_artifact_directory(self) -> None:
        self._write_config(
            artifact_dir=str(
                self.pic_root / "runs" / "f0_hipmpi_smoke" / str(uuid.uuid4())
            )
        )
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_reservation_rejects_operational_namespace_artifact_directory(self) -> None:
        self._write_config(artifact_dir=str(self.pic_root / "ledger" / "forged-run"))
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_executable_lookup_rejects_cancelled_reservation(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        transition(
            reservation_id=reservation_id,
            notes="cancel before scheduler submission",
            event_type="reservation_cancelled",
            state="cancelled",
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            receipts_jsonl=self.receipts,
            mirror_jsonl=self.mirror,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        scheduler = (
            f"JobId=12345 JobState=RUNNING Account=AST207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            with self.assertRaises(ValueError):
                reservation_bound_manifest(
                    manifest_path,
                    reservation_id,
                    ledger_jsonl=self.ledger,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    executable_job_id="12345",
                )

    def test_compute_node_snapshot_lookup_does_not_take_mutation_lock(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        self._attach(reservation_id)
        scheduler = (
            f"JobId=12345 JobState=RUNNING Account=AST207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            with patch(
                "validate_and_reserve_frontier_job.ledger_lock",
                side_effect=AssertionError("compute snapshot attempted mutation lock"),
            ):
                manifest, loaded = executable_reservation_bound_manifest(
                    manifest_path,
                    reservation_id,
                    ledger_jsonl=self.ledger,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    executable_job_id="12345",
                )
        self.assertEqual(manifest["submission_id"], self.submission_id)
        self.assertEqual(loaded["reservation_id"], reservation_id)

    def test_predispatch_snapshot_lookup_retains_mutation_lock(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        with patch(
            "validate_and_reserve_frontier_job.ledger_lock",
            side_effect=OSError(524, "Unknown error 524"),
        ):
            with self.assertRaises(OSError):
                reservation_bound_manifest(
                    manifest_path,
                    str(reservation["reservation_id"]),
                    ledger_jsonl=self.ledger,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    require_reserved=True,
                )

    def test_compute_node_snapshot_lookup_rejects_receipt_drift_during_use(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        self._attach(reservation_id)
        scheduler = (
            f"JobId=12345 JobState=RUNNING Account=AST207 "
            f"Comment=pic-reservation={reservation_id}"
        )

        def mutate_receipt(*args: object, **kwargs: object) -> str:
            self.receipts.write_text(
                self.receipts.read_text(encoding="utf-8") + "\n",
                encoding="utf-8",
            )
            return scheduler

        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            side_effect=mutate_receipt,
        ):
            with self.assertRaisesRegex(ValueError, "changed during snapshot use"):
                executable_reservation_bound_manifest(
                    manifest_path,
                    reservation_id,
                    ledger_jsonl=self.ledger,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    executable_job_id="12345",
                )

    def test_predispatch_lookup_rejects_cancelled_reservation(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        transition(
            reservation_id=reservation_id,
            notes="cancel before scheduler submission",
            event_type="reservation_cancelled",
            state="cancelled",
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            receipts_jsonl=self.receipts,
            mirror_jsonl=self.mirror,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        with self.assertRaises(ValueError):
            reservation_bound_manifest(
                manifest_path,
                reservation_id,
                ledger_jsonl=self.ledger,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                require_reserved=True,
            )

    def test_policy_promotion_rejects_reserved_submission_and_retains_accounting(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        self.assertRegex(str(reservation["active_policy_sha256"]), r"^[0-9a-f]{64}$")
        self.assertRegex(str(reservation["active_promotion_sha256"]), r"^[0-9a-f]{64}$")
        self.policy.write_text(
            self.policy.read_text(encoding="utf-8") + "\n",
            encoding="utf-8",
        )
        with self.assertRaises(ValueError):
            self._promote_policy()
        totals = accounting(validate_primary_chain(self.ledger))
        self.assertGreater(totals["currently_reserved_node_hours"], 0.0)
        self._mark_dispatch_started(reservation_id)

    def test_dispatch_started_reservation_cannot_be_cancelled_or_repaired(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        self._mark_dispatch_started(reservation_id)
        marker = json.loads(
            (self.pic_root / "ledger" / "pending_submission.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(marker["state"], "scheduler_dispatch_started")
        with self.assertRaises(ValueError):
            transition(
                reservation_id=reservation_id,
                notes="ambiguous sbatch result must retain accounting",
                event_type="reservation_cancelled",
                state="cancelled",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with self.assertRaises(ValueError):
            repair_reservation_attachments(
                reservation_id=reservation_id,
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_current_generation_paths_reject_marker_control_plane_mutation(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        marker_path = self.pic_root / "ledger" / "pending_submission.json"
        marker = json.loads(marker_path.read_text(encoding="utf-8"))
        marker["control_plane_version"] = "0" * 64
        marker_path.chmod(0o600)
        marker_path.write_text(json.dumps(marker), encoding="utf-8")
        marker_path.chmod(0o400)
        with self.assertRaises(ValueError):
            self._mark_dispatch_started(reservation_id)
        with self.assertRaises(ValueError):
            repair_reservation_attachments(
                reservation_id=reservation_id,
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        retained = json.loads(marker_path.read_text(encoding="utf-8"))
        self.assertEqual(retained["state"], "reserved_not_submitted")

    def test_successor_cannot_mark_prior_generation_submission(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        self._mark_dispatch_started(reservation_id)
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        scheduler = (
            f"JobId=12345 JobState=PENDING Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=successor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        marker = json.loads(
            (self.pic_root / "ledger" / "pending_submission.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(marker["state"], "scheduler_dispatch_started")

    def test_scheduler_output_is_restricted_to_dedicated_log_path(self) -> None:
        expected = self.pic_root / "logs" / "slurm" / "%x.%j.log"
        self.assertEqual(
            _require_scheduler_output_path(str(expected), self.pic_root),
            expected,
        )
        with self.assertRaises(ValueError):
            _require_scheduler_output_path(str(self.ledger), self.pic_root)

    def test_executable_lookup_rejects_pre_attachment_race(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=RUNNING Account=AST207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            mark_submitted(
                reservation_id=reservation_id,
                job_id="12345",
                ledger_jsonl=self.ledger,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
            with self.assertRaises(ValueError):
                reservation_bound_manifest(
                    manifest_path,
                    reservation_id,
                    ledger_jsonl=self.ledger,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    executable_job_id="12345",
                )

    def test_cancellation_rejects_submitted_not_attached_job(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=PENDING Account=AST207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            mark_submitted(
                reservation_id=reservation_id,
                job_id="12345",
                ledger_jsonl=self.ledger,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with self.assertRaises(ValueError):
            transition(
                reservation_id=reservation_id,
                notes="must not release accounting after sbatch",
                event_type="reservation_cancelled",
                state="cancelled",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_lookup_rejects_missing_paired_install_before_lock_creation(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        lock = self.ledger.with_suffix(self.ledger.suffix + ".lock")
        lock.unlink()
        remove_tree(self.project_home_control_plane_dir)
        with self.assertRaises(ValueError):
            reservation_bound_manifest(
                manifest_path,
                str(reservation["reservation_id"]),
                ledger_jsonl=self.ledger,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertFalse(lock.exists())

    def test_trampoline_does_not_open_shell_template_after_verification(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        template = Path(str(record_for_role(manifest, "job-script")["path"]))
        commands: list[list[str]] = []

        def runner(command: list[str], **_: object) -> None:
            template.chmod(0o644)
            template.write_text("#!/bin/bash\n/bin/true\n", encoding="utf-8")
            commands.append(command)

        self._launch(manifest_path, reservation, runner=runner)
        self.assertEqual(len(commands), 1)
        self.assertRegex(
            commands[0][0],
            r"^/proc/self/fd/[0-9]+/launch_with_frontier_profile[.]sh$",
        )
        self.assertEqual(commands[0][1], "/usr/bin/srun")
        self.assertEqual(commands[0][2], "--jobid=12345")
        self.assertNotIn(str(template), commands[0])

    def test_manifest_rejects_arbitrary_shell_launch_action(self) -> None:
        self._write_config()
        config = json.loads(self.config.read_text(encoding="utf-8"))
        config["launch_contract"]["actions"][0]["kind"] = "shell"
        config["launch_contract"]["actions"][0]["executable"] = "/bin/true"
        self.config.write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaises(ValueError):
            self._create_manifest()

    def test_launch_contract_rejects_untrusted_path_bearing_literals(self) -> None:
        unsafe_arguments = [
            [{"literal": "-i"}, {"literal": "/outside/input.athinput"}],
            [{"literal": "-d"}, {"literal": "/tmp"}],
            [{"literal": "-r"}, {"literal": "/outside/restart.rst"}],
            [{"literal": "-t"}, {"literal": "00:00:01"}],
            [{"literal": "-t=00:00:01"}],
            [{"literal": "job/basename=../../outside"}],
        ]
        for arguments in unsafe_arguments:
            with self.subTest(arguments=arguments):
                contract = self._launch_contract()
                contract["actions"][0]["arguments"] = arguments
                with self.assertRaises(ValueError):
                    validate_launch_contract(contract)

    def test_launch_contract_rejects_noninteger_schema_version(self) -> None:
        for value in (True, 1.0, "1"):
            with self.subTest(value=value):
                contract = self._launch_contract()
                contract["schema_version"] = value
                with self.assertRaisesRegex(ValueError, "launch-contract schema"):
                    validate_launch_contract(contract)

    def test_launch_contract_rejects_noncanonical_artifact_paths(self) -> None:
        for value in ("output/./stdout.txt", "output//stdout.txt", "output/"):
            with self.subTest(value=value):
                contract = self._launch_contract()
                contract["actions"][0]["stdout_artifact"] = value
                with self.assertRaisesRegex(ValueError, "relative artifact path"):
                    validate_launch_contract(contract)

    def test_admission_smoke_rejects_valid_but_unbound_launch_resource_drift(
        self,
    ) -> None:
        config = json.loads(self.config.read_text(encoding="utf-8"))
        config["launch_contract"]["actions"][0]["resources"]["tasks"] = 2
        self.config.write_text(json.dumps(config), encoding="utf-8")
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_snapshot_escape_rejects_before_reading_outside_root(self) -> None:
        outside = self.root / "outside-snapshot"
        outside.write_text("must not be read\n", encoding="utf-8")
        manifest = {
            "snapshot_files": [
                {
                    "role": "input-deck",
                    "path": str(outside),
                    "sha256": sha256(outside),
                    "source_path": str(outside),
                    "source_sha256": sha256(outside),
                }
            ]
        }
        with patch(
            "control_plane_common.read_stable_regular_file_below"
        ) as stable_read:
            with self.assertRaises(ValueError):
                verify_snapshot_files(manifest, root=self.pic_root)
        stable_read.assert_not_called()

    def test_trampoline_runs_only_declarative_bounded_hooks(self) -> None:
        self._write_config()
        config = json.loads(self.config.read_text(encoding="utf-8"))
        config["launch_contract"]["pre_actions"] = [
            {
                "action_id": "record-executable",
                "kind": "snapshot_sha256",
                "snapshot_role": "executable",
                "output_artifact": "checksums/executable.sha256",
            }
        ]
        config["launch_contract"]["post_actions"] = [
            {
                "action_id": "record-stdout",
                "kind": "artifact_sha256",
                "artifact": "athena_stdout.txt",
                "output_artifact": "checksums/athena_stdout.sha256",
            }
        ]
        self.config.write_text(json.dumps(config), encoding="utf-8")
        self._write_policy(
            admission_smoke_overrides={
                "launch_contract_sha256": launch_contract_sha256(
                    config["launch_contract"]
                ),
            }
        )
        self._promote_policy()
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        calls: list[tuple[list[str], dict[str, object]]] = []

        def runner(command: list[str], **kwargs: object) -> None:
            calls.append((command, kwargs))

        self._launch(manifest_path, reservation, runner=runner)
        self.assertEqual(len(calls), 1)
        self.assertRegex(
            calls[0][0][0],
            r"^/proc/self/fd/[0-9]+/launch_with_frontier_profile[.]sh$",
        )
        self.assertEqual(calls[0][0][1], "/usr/bin/srun")
        self.assertEqual(calls[0][0][2], "--jobid=12345")
        artifact_dir = self.pic_root / "runs" / "f0_hipmpi_smoke" / self.submission_id
        self.assertRegex(calls[0][1]["env"]["PIC_RUNTIME_ALLOWLIST_FD"], r"^[0-9]+$")
        self.assertRegex(
            calls[0][1]["env"]["PIC_RUNTIME_ALLOWLIST_DIR_FD"], r"^[0-9]+$"
        )
        self.assertRegex(calls[0][1]["env"]["PIC_CONTROL_PLANE_DIR_FD"], r"^[0-9]+$")
        self.assertEqual(len(calls[0][1]["pass_fds"]), 3)
        self.assertTrue(
            (artifact_dir / "athena-parser.environment.allowlist.txt").is_file()
        )
        self.assertEqual(
            (artifact_dir / "checksums" / "executable.sha256").read_text(
                encoding="utf-8"
            ),
            sha256(Path(str(record_for_role(
                json.loads(manifest_path.read_text(encoding="utf-8")),
                "executable",
            )["path"]))) + "\n",
        )
        self.assertEqual(
            (artifact_dir / "checksums" / "athena_stdout.sha256").read_text(
                encoding="utf-8"
            ),
            sha256(artifact_dir / "athena_stdout.txt") + "\n",
        )

    def test_manifest_rejects_user_python_hook(self) -> None:
        self._write_config()
        config = json.loads(self.config.read_text(encoding="utf-8"))
        config["launch_contract"]["post_actions"] = [
            {
                "action_id": "analysis",
                "kind": "analysis",
                "analysis_role": "analysis-script-000",
            }
        ]
        self.config.write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaises(ValueError):
            self._create_manifest()

    def test_manifest_preserves_authorized_analysis_support_module_basename(self) -> None:
        self._write(
            "analysis.py",
            "import importlib.util\n"
            "from pathlib import Path\n"
            "path = Path(__file__).with_name('frontier_f1_structured_artifacts.py')\n"
            "spec = importlib.util.spec_from_file_location('_verified_helper', path)\n"
            "module = importlib.util.module_from_spec(spec)\n"
            "spec.loader.exec_module(module)\n"
            "print(module.VALUE)\n",
        )
        helper = self._write(
            "frontier_f1_structured_artifacts.py", "VALUE = 'verified helper bytes'\n"
        )
        self._write_config(
            analysis_scripts=[str(self.sources / "analysis.py"), str(helper)]
        )
        manifest_path = self._create_manifest()
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        analysis_paths = {
            record["role"]: Path(str(record["path"])).name
            for record in manifest["snapshot_files"]
            if str(record["role"]).startswith("analysis-script-")
        }
        self.assertEqual(
            analysis_paths,
            {
                "analysis-script-000": "000-analysis.py",
                "analysis-script-001": "frontier_f1_structured_artifacts.py",
            },
        )
        snapshot_analysis = manifest_path.parent / "snapshot" / "analysis"
        self.assertEqual(
            subprocess.check_output(
                [TRUSTED_PYTHON, "-I", str(snapshot_analysis / "000-analysis.py")],
                text=True,
            ),
            "verified helper bytes\n",
        )

    def test_reserved_launch_rejects_self_consistent_snapshot_attachment_rewrite(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        executable = record_for_role(manifest, "executable")
        snapshot = Path(str(executable["path"]))
        snapshot.chmod(0o644)
        snapshot.write_text("forged after reservation\n", encoding="utf-8")
        snapshot.chmod(0o444)
        executable["sha256"] = sha256(snapshot)
        executable["source_sha256"] = sha256(snapshot)
        manifest_path.chmod(0o644)
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        manifest_path.chmod(0o444)
        attachment = manifest_path.parent / "manifest_sha256.txt"
        attachment.chmod(0o644)
        attachment.write_text(sha256(manifest_path) + "\n", encoding="utf-8")
        attachment.chmod(0o444)
        with self.assertRaises(ValueError):
            reservation_bound_manifest(
                manifest_path,
                str(reservation["reservation_id"]),
                ledger_jsonl=self.ledger,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_reserved_launch_rejects_private_ledger_substitution(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        fake = self.root / "private-ledger"
        with self.assertRaises(ValueError):
            reservation_bound_manifest(
                manifest_path,
                str(reservation["reservation_id"]),
                ledger_jsonl=fake / "node_hours.jsonl",
                receipts_jsonl=fake / "mirror_receipts.jsonl",
                mirror_jsonl=fake / "mirror.jsonl",
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_interrupted_reservation_attachment_write_is_repairable(self) -> None:
        manifest_path = self._create_manifest()
        with patch(
            "validate_and_reserve_frontier_job._write_reservation_attachments",
            side_effect=RuntimeError("simulated attachment interruption"),
        ):
            with self.assertRaises(RuntimeError):
                self._reserve(manifest_path)
        repair_reservation_attachments(
            reservation_id="89c76745-6c37-47f7-9847-800a98a47c9b",
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            receipts_jsonl=self.receipts,
            mirror_jsonl=self.mirror,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        marker = json.loads(
            (self.pic_root / "ledger" / "pending_submission.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(marker["state"], "reserved_not_submitted")
        for name in ["reservation_id.txt", "manifest_sha256.txt"]:
            attachment = manifest_path.parent / name
            self.assertTrue(attachment.is_file())
            self.assertFalse(bool(attachment.stat().st_mode & 0o222))

    def test_unappended_reservation_intent_is_recoverable(self) -> None:
        manifest_path = self._create_manifest()
        with patch(
            "validate_and_reserve_frontier_job._append_locked",
            side_effect=RuntimeError("simulated primary append interruption"),
        ):
            with self.assertRaises(RuntimeError):
                self._reserve(manifest_path)
        result = repair_reservation_attachments(
            reservation_id="89c76745-6c37-47f7-9847-800a98a47c9b",
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            receipts_jsonl=self.receipts,
            mirror_jsonl=self.mirror,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        self.assertEqual(result, "cleared_unappended_reservation_intent")
        self.assertFalse((self.pic_root / "ledger" / "pending_submission.json").exists())
        self.assertEqual(len(validate_primary_chain(self.ledger)), 1)

    def test_unappended_intent_repair_rejects_successor_control_plane(self) -> None:
        manifest_path = self._create_manifest()
        with patch(
            "validate_and_reserve_frontier_job._append_locked",
            side_effect=RuntimeError("simulated primary append interruption"),
        ):
            with self.assertRaises(RuntimeError):
                self._reserve(manifest_path)
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        with self.assertRaises(ValueError):
            repair_reservation_attachments(
                reservation_id="89c76745-6c37-47f7-9847-800a98a47c9b",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertTrue((self.pic_root / "ledger" / "pending_submission.json").is_file())

    def test_reserved_attachment_repair_rejects_successor_control_plane(self) -> None:
        manifest_path = self._create_manifest()
        with patch(
            "validate_and_reserve_frontier_job._write_reservation_attachments",
            side_effect=RuntimeError("simulated attachment interruption"),
        ):
            with self.assertRaises(RuntimeError):
                self._reserve(manifest_path)
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        with self.assertRaises(ValueError):
            repair_reservation_attachments(
                reservation_id="89c76745-6c37-47f7-9847-800a98a47c9b",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertTrue((self.pic_root / "ledger" / "pending_submission.json").is_file())

    def test_completed_cancellation_stranded_marker_is_recoverable(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        with patch(
            "validate_and_reserve_frontier_job._clear_matching_pending_marker",
            side_effect=RuntimeError("simulated marker cleanup interruption"),
        ):
            with self.assertRaises(RuntimeError):
                transition(
                    reservation_id=reservation_id,
                    event_type="reservation_cancelled",
                    state="cancelled",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        result = repair_reservation_attachments(
            reservation_id=reservation_id,
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            receipts_jsonl=self.receipts,
            mirror_jsonl=self.mirror,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        self.assertEqual(result, "cleared_completed_cancellation_pending_marker")
        self.assertFalse((self.pic_root / "ledger" / "pending_submission.json").exists())

    def test_attachment_repair_rejects_stale_policy_generation(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        transition(
            reservation_id=reservation_id,
            event_type="reservation_cancelled",
            state="cancelled",
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            receipts_jsonl=self.receipts,
            mirror_jsonl=self.mirror,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        self._promote_test_control_plane_successor(successor)
        marker_path = self.pic_root / "ledger" / "pending_submission.json"
        marker_path.write_text(
            json.dumps({"reservation_id": reservation_id}),
            encoding="utf-8",
        )
        with self.assertRaises(ValueError):
            repair_reservation_attachments(
                reservation_id=reservation_id,
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertTrue(marker_path.is_file())

    def test_completed_attachment_stranded_marker_is_recoverable(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            "JobId=12345 JobState=PENDING Account=AST207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            mark_submitted(
                reservation_id=reservation_id,
                job_id="12345",
                ledger_jsonl=self.ledger,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
            with patch(
                "validate_and_reserve_frontier_job._clear_matching_pending_marker",
                side_effect=RuntimeError("simulated marker cleanup interruption"),
            ):
                with self.assertRaises(RuntimeError):
                    transition(
                        reservation_id=reservation_id,
                        job_id="12345",
                        event_type="job_id_attached",
                        state="submitted",
                        ledger_jsonl=self.ledger,
                        ledger_csv=self.csv,
                        receipts_jsonl=self.receipts,
                        mirror_jsonl=self.mirror,
                        control_plane_dir=self.control_plane_dir,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                    )
        result = repair_reservation_attachments(
            reservation_id=reservation_id,
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            receipts_jsonl=self.receipts,
            mirror_jsonl=self.mirror,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        self.assertEqual(result, "cleared_completed_attachment_pending_marker")
        self.assertFalse((self.pic_root / "ledger" / "pending_submission.json").exists())

    def test_attachment_repair_rejects_falsey_terminal_recovery_provenance(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        synthetic = dict(reservation)
        synthetic.update(
            {
                "event_type": "job_id_attached",
                "state": "submitted",
                "terminal_recovery_handoff_path": "",
                "terminal_recovery_handoff_sha256": "",
                "terminal_recovery_mode": "",
            }
        )
        with patch(
            "validate_and_reserve_frontier_job.latest_reservations",
            return_value={reservation_id: synthetic},
        ):
            with self.assertRaises(ValueError):
                repair_reservation_attachments(
                    reservation_id=reservation_id,
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        self.assertTrue((self.pic_root / "ledger" / "pending_submission.json").is_file())

    def test_mutable_source_tree_cannot_run_mutating_ledger_repairs(self) -> None:
        common = {
            "ledger_jsonl": self.ledger,
            "ledger_csv": self.csv,
            "receipts_jsonl": self.receipts,
            "mirror_jsonl": self.mirror,
            "authorized_pic_root": self.pic_root,
            "authorized_project_home_root": self.project_home_root,
        }
        with patch(
            "validate_and_reserve_frontier_job.repair_mirrored_state_locked"
        ) as repair:
            with self.assertRaises(ValueError):
                repair_ledger_mirror(**common)
            with self.assertRaises(ValueError):
                repair_reservation_attachments(
                    reservation_id="89c76745-6c37-47f7-9847-800a98a47c9b",
                    **common,
                )
        repair.assert_not_called()

    def test_closed_policy_rejects_interrupted_looking_genesis_prefix(self) -> None:
        for path in [
            self.csv,
            self.receipts,
            self.mirror,
            self.pic_root / "ledger" / "genesis_anchor.json",
            self.project_home_root / "ledger" / "genesis_anchor.json",
        ]:
            if path.exists():
                path.chmod(0o600)
                path.unlink()
        with self.assertRaises(ValueError):
            initialize_from_policy(
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                mirror_receipts=self.receipts,
                mirror_jsonl=self.mirror,
                mirror_transport="filesystem_copy",
                notes="snapshot control plane test",
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_reservation_creation_rejects_broken_attachment_symlink_alias(self) -> None:
        manifest_path = self._create_manifest()
        outside = self.root / "outside-reservation-attachment"
        attachment = manifest_path.parent / "reservation_id.txt"
        attachment.symlink_to(outside)
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)
        self.assertTrue(attachment.is_symlink())
        self.assertFalse(outside.exists())

    def test_reservation_lookup_rejects_matching_attachment_symlink_alias(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        attachment = manifest_path.parent / "reservation_id.txt"
        outside = self.root / "outside-matching-reservation-id.txt"
        outside.write_bytes(attachment.read_bytes())
        outside.chmod(0o444)
        attachment.unlink()
        attachment.symlink_to(outside)
        with self.assertRaises(ValueError):
            reservation_bound_manifest(
                manifest_path,
                str(reservation["reservation_id"]),
                ledger_jsonl=self.ledger,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_compute_verify_rejects_matching_attachment_symlink_alias(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        attachment = manifest_path.parent / "manifest_sha256.txt"
        outside = self.root / "outside-matching-manifest-sha256.txt"
        outside.write_bytes(attachment.read_bytes())
        outside.chmod(0o444)
        attachment.unlink()
        attachment.symlink_to(outside)
        with self.assertRaises(ValueError):
            verify(
                manifest_path,
                reservation_id=str(reservation["reservation_id"]),
                manifest_sha256=str(reservation["manifest_sha256"]),
            )

    def test_reservation_lookup_rejects_noncanonical_lf_attachment_bytes(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        attachment = manifest_path.parent / "reservation_id.txt"
        attachment.chmod(0o644)
        attachment.write_bytes(attachment.read_bytes() + b"\n")
        attachment.chmod(0o444)
        with self.assertRaisesRegex(ValueError, "Reservation ID"):
            reservation_bound_manifest(
                manifest_path,
                str(reservation["reservation_id"]),
                ledger_jsonl=self.ledger,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_compute_verify_rejects_noncanonical_lf_attachment_bytes(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        attachment = manifest_path.parent / "manifest_sha256.txt"
        attachment.chmod(0o644)
        attachment.write_bytes(attachment.read_bytes() + b"\n")
        attachment.chmod(0o444)
        with self.assertRaisesRegex(ValueError, "Manifest checksum"):
            verify(
                manifest_path,
                reservation_id=str(reservation["reservation_id"]),
                manifest_sha256=str(reservation["manifest_sha256"]),
            )

    def test_attachment_repair_rejects_matching_symlink_alias(self) -> None:
        manifest_path = self._create_manifest()
        with patch(
            "validate_and_reserve_frontier_job._write_reservation_attachments",
            side_effect=RuntimeError("simulated attachment interruption"),
        ):
            with self.assertRaises(RuntimeError):
                self._reserve(manifest_path)
        outside = self.root / "outside-repair-reservation-id.txt"
        outside.write_text("89c76745-6c37-47f7-9847-800a98a47c9b\n", encoding="utf-8")
        outside.chmod(0o444)
        attachment = manifest_path.parent / "reservation_id.txt"
        attachment.symlink_to(outside)
        with self.assertRaises(ValueError):
            repair_reservation_attachments(
                reservation_id="89c76745-6c37-47f7-9847-800a98a47c9b",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(
            outside.read_text(encoding="utf-8"),
            "89c76745-6c37-47f7-9847-800a98a47c9b\n",
        )

    def test_unappended_intent_repair_rejects_broken_attachment_symlink_alias(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        with patch(
            "validate_and_reserve_frontier_job._append_locked",
            side_effect=RuntimeError("simulated primary append interruption"),
        ):
            with self.assertRaises(RuntimeError):
                self._reserve(manifest_path)
        outside = self.root / "outside-broken-repair-reservation-id.txt"
        attachment = manifest_path.parent / "reservation_id.txt"
        attachment.symlink_to(outside)
        with self.assertRaises(ValueError):
            repair_reservation_attachments(
                reservation_id="89c76745-6c37-47f7-9847-800a98a47c9b",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertTrue(attachment.is_symlink())
        self.assertFalse(outside.exists())

    def test_attach_rejects_scheduler_comment_or_state_mismatch(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value="JobId=12345 JobState=PENDING Account=AST207 Comment=wrong",
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=(
                f"JobId=12345 JobState=COMPLETED Account=AST207 "
                f"Comment=pic-reservation={reservation_id}"
            ),
        ):
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=(
                f"JobId=other JobState=PENDING Account=AST207 "
                f"Comment=pic-reservation={reservation_id}"
            ),
        ):
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )

    def test_attach_accepts_scheduler_canonical_lowercase_account(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=PENDING Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            mark_submitted(
                reservation_id=reservation_id,
                job_id="12345",
                ledger_jsonl=self.ledger,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_scheduler_account_matching_is_closed_to_canonical_spellings(self) -> None:
        self.assertTrue(scheduler_account_matches_authorized("AST207"))
        self.assertTrue(scheduler_account_matches_authorized("ast207"))
        for value in ["Ast207", "ast207-extra", "", None]:
            self.assertFalse(scheduler_account_matches_authorized(value))

    def test_submission_directive_rejects_lowercase_account(self) -> None:
        self._write(
            "job.sh",
            "#!/bin/bash\n#SBATCH -A ast207\n#SBATCH -p batch\n#SBATCH -q debug\n"
            f"#SBATCH -o {self.pic_root}/logs/slurm/%x.%j.log\n"
            "#SBATCH -N 1\n#SBATCH -t 00:10:00\n",
        )
        self._write_policy()
        self._promote_policy()
        with self.assertRaises(ValueError):
            self._reserve(self._create_manifest())

    def test_policy_rejects_lowercase_account(self) -> None:
        policy = json.loads(self.policy.read_text(encoding="utf-8"))
        policy["frontier"]["account"] = "ast207"
        self.policy.write_text(json.dumps(policy), encoding="utf-8")
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_reconcile_queries_slurm_and_rejects_comment_mismatch(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        self._attach(str(reservation["reservation_id"]))
        self.assertNotIn("state", inspect.signature(reconcile).parameters)
        self.assertNotIn("elapsed_seconds", inspect.signature(reconcile).parameters)
        self.assertNotIn("allocated_nodes", inspect.signature(reconcile).parameters)
        with patch(
            "reconcile_frontier_job.subprocess.check_output",
            return_value="12345|COMPLETED|300|1|wrong|AST207\n",
        ):
            with self.assertRaises(ValueError):
                reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )

    def test_scheduler_accounting_accepts_empty_comment_and_lowercase_account(self) -> None:
        import reconcile_frontier_job

        with patch(
            "reconcile_frontier_job._verify_scheduler_job_binding"
        ) as scheduler_binding:
            with patch.object(
                reconcile_frontier_job.subprocess,
                "check_output",
                return_value="12345|CANCELLED by 18664|0|0||ast207\n",
            ):
                self.assertEqual(
                    reconcile_frontier_job._scheduler_result("12345", "reservation"),
                    ("CANCELLED", 0, 0),
                )
        scheduler_binding.assert_called_once_with("12345", "reservation")
        with patch(
            "reconcile_frontier_job._verify_scheduler_job_binding",
            side_effect=ValueError("missing scheduler binding"),
        ):
            with patch.object(
                reconcile_frontier_job.subprocess,
                "check_output",
                return_value="12345|CANCELLED|0|0||ast207\n",
            ):
                with self.assertRaises(ValueError):
                    reconcile_frontier_job._scheduler_result("12345", "reservation")
        with patch.object(
            reconcile_frontier_job.subprocess,
            "check_output",
            return_value="12345|CANCELLED|0|0||wrong\n",
        ):
            with self.assertRaises(ValueError):
                reconcile_frontier_job._scheduler_result("12345", "reservation")
        for output in [
            "12345|CANCELLED|0|0|pic-reservation=reservation|ast207|extra\n",
            "12345|CANCELLED|0|0|pic-reservation=reservation|ast207\n"
            "12345|CANCELLED|0|0|pic-reservation=reservation|ast207\n",
        ]:
            with self.subTest(output=output):
                with patch.object(
                    reconcile_frontier_job.subprocess,
                    "check_output",
                    return_value=output,
                ):
                    with self.assertRaises(ValueError):
                        reconcile_frontier_job._scheduler_result(
                            "12345", "reservation"
                        )

    def test_purged_cancelled_zero_execution_snapshot_is_exact(self) -> None:
        row = (
            "12345|run_installed_control_plane_job.sh|CANCELLED by 18664|0|0||"
            "ast207|2026-05-30T15:47:31|None|2026-05-30T15:47:31|0:0\n"
        )

        def scheduler_output(command: list[str], **kwargs: object) -> str:
            if command[0] == TRUSTED_SCONTROL:
                raise subprocess.CalledProcessError(
                    1,
                    command,
                    stderr="slurm_load_jobs error: Invalid job id specified\n",
                )
            if command[0] == TRUSTED_SQUEUE:
                return ""
            self.assertEqual(command[0], TRUSTED_SACCT)
            return row

        with patch.object(
            terminal_recovery_handoff.subprocess,
            "check_output",
            side_effect=scheduler_output,
        ):
            snapshot = (
                terminal_recovery_handoff.require_purged_cancelled_zero_execution_snapshot(
                    "12345"
                )
            )
        self.assertEqual(
            snapshot["mode"],
            terminal_recovery_handoff.PURGED_CANCELLED_ZERO_EXECUTION_MODE,
        )
        self.assertEqual(snapshot["start"], "None")
        self.assertEqual(snapshot["elapsed_raw"], 0)
        self.assertEqual(snapshot["allocated_nodes"], 0)

    def test_purged_cancelled_zero_execution_snapshot_rejects_broader_shapes(
        self,
    ) -> None:
        fields = [
            "12345",
            "run_installed_control_plane_job.sh",
            "CANCELLED by 18664",
            "0",
            "0",
            "",
            "ast207",
            "2026-05-30T15:47:31",
            "None",
            "2026-05-30T15:47:31",
            "0:0",
        ]

        def purged_scontrol(command: list[str], **kwargs: object) -> str:
            if command[0] == TRUSTED_SCONTROL:
                raise subprocess.CalledProcessError(
                    1,
                    command,
                    stderr="slurm_load_jobs error: Invalid job id specified\n",
                )
            if command[0] == TRUSTED_SQUEUE:
                return ""
            return "|".join(fields) + "\n"

        for index, value in [
            (1, "other.sh"),
            (2, "COMPLETED"),
            (3, "1"),
            (4, "1"),
            (5, "pic-reservation=reservation"),
            (6, "wrong"),
            (8, "2026-05-30T15:47:31"),
            (9, "2026-05-30T15:47:32"),
            (10, "1:0"),
        ]:
            with self.subTest(index=index, value=value):
                original = fields[index]
                fields[index] = value
                try:
                    with patch.object(
                        terminal_recovery_handoff.subprocess,
                        "check_output",
                        side_effect=purged_scontrol,
                    ):
                        with self.assertRaises(ValueError):
                            terminal_recovery_handoff.require_purged_cancelled_zero_execution_snapshot(
                                "12345"
                            )
                finally:
                    fields[index] = original
        fields.append("unexpected")
        try:
            with patch.object(
                terminal_recovery_handoff.subprocess,
                "check_output",
                side_effect=purged_scontrol,
            ):
                with self.assertRaises(ValueError):
                    terminal_recovery_handoff.require_purged_cancelled_zero_execution_snapshot(
                        "12345"
                    )
        finally:
            fields.pop()
        with patch.object(
            terminal_recovery_handoff.subprocess,
            "check_output",
            side_effect=subprocess.CalledProcessError(
                2,
                [TRUSTED_SCONTROL],
                stderr="slurm_load_jobs error: Invalid job id specified\n",
            ),
        ):
            with self.assertRaises(ValueError):
                terminal_recovery_handoff.require_purged_cancelled_zero_execution_snapshot(
                    "12345"
                )
        with patch.object(
            terminal_recovery_handoff.subprocess,
            "check_output",
            side_effect=subprocess.CalledProcessError(
                1,
                [TRUSTED_SCONTROL],
                stderr="slurm_load_jobs error: Access/permission denied\n",
            ),
        ):
            with self.assertRaises(ValueError):
                terminal_recovery_handoff.require_purged_cancelled_zero_execution_snapshot(
                    "12345"
                )

    def test_purged_cancelled_zero_execution_snapshot_rejects_queue_and_duplicate_rows(
        self,
    ) -> None:
        row = (
            "12345|run_installed_control_plane_job.sh|CANCELLED by 18664|0|0||"
            "ast207|2026-05-30T15:47:31|None|2026-05-30T15:47:31|0:0\n"
        )
        queued = ""
        accounting = row

        def scheduler_output(command: list[str], **kwargs: object) -> str:
            if command[0] == TRUSTED_SCONTROL:
                raise subprocess.CalledProcessError(
                    1,
                    command,
                    stderr="slurm_load_jobs error: Invalid job id specified\n",
                )
            if command[0] == TRUSTED_SQUEUE:
                return queued
            self.assertEqual(command[0], TRUSTED_SACCT)
            return accounting

        with patch.object(
            terminal_recovery_handoff.subprocess,
            "check_output",
            side_effect=scheduler_output,
        ):
            queued = "12345\n"
            with self.assertRaises(ValueError):
                terminal_recovery_handoff.require_purged_cancelled_zero_execution_snapshot(
                    "12345"
                )
            queued = ""
            accounting = row + row
            with self.assertRaises(ValueError):
                terminal_recovery_handoff.require_purged_cancelled_zero_execution_snapshot(
                    "12345"
                )
            accounting = row + "12345|truncated\n"
            with self.assertRaises(ValueError):
                terminal_recovery_handoff.require_purged_cancelled_zero_execution_snapshot(
                    "12345"
                )

    def test_purged_cancelled_zero_execution_snapshot_accepts_exact_squeue_purge_only(
        self,
    ) -> None:
        row = (
            "12345|run_installed_control_plane_job.sh|CANCELLED by 18664|0|0||"
            "ast207|2026-05-30T15:47:31|None|2026-05-30T15:47:31|0:0\n"
        )
        queue_stderr = "slurm_load_jobs error: Invalid job id specified\n"
        squeue_returncode = 1

        def scheduler_output(command: list[str], **kwargs: object) -> str:
            if command[0] == TRUSTED_SCONTROL:
                raise subprocess.CalledProcessError(
                    1,
                    command,
                    stderr="slurm_load_jobs error: Invalid job id specified\n",
                )
            if command[0] == TRUSTED_SQUEUE:
                raise subprocess.CalledProcessError(
                    squeue_returncode,
                    command,
                    stderr=queue_stderr,
                )
            self.assertEqual(command[0], TRUSTED_SACCT)
            return row

        with patch.object(
            terminal_recovery_handoff.subprocess,
            "check_output",
            side_effect=scheduler_output,
        ):
            snapshot = (
                terminal_recovery_handoff.require_purged_cancelled_zero_execution_snapshot(
                    "12345"
                )
            )
            self.assertEqual(snapshot["job_id"], "12345")
            squeue_returncode = 2
            with self.assertRaises(ValueError):
                terminal_recovery_handoff.require_purged_cancelled_zero_execution_snapshot(
                    "12345"
                )
            squeue_returncode = 1
            queue_stderr = "slurm_load_jobs error: Access/permission denied\n"
            with self.assertRaises(ValueError):
                terminal_recovery_handoff.require_purged_cancelled_zero_execution_snapshot(
                    "12345"
                )

    def test_reconcile_recovers_terminal_submitted_not_attached_job(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=PENDING Account=AST207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            mark_submitted(
                reservation_id=reservation_id,
                job_id="12345",
                ledger_jsonl=self.ledger,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with patch(
            "reconcile_frontier_job._scheduler_result",
            return_value=("CANCELLED", 30, 1),
        ):
            result = reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(result["state"], "CANCELLED")
        self.assertTrue(result["reconciled"])
        self.assertFalse((self.pic_root / "ledger" / "pending_submission.json").exists())
        records = validate_primary_chain(self.ledger)
        self.assertEqual(records[-2]["event_type"], "job_id_attached")
        self.assertEqual(records[-1]["event_type"], "reconciliation")

    def test_reconcile_recovers_terminal_received_job_id(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=AST207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
            marker = json.loads(
                (self.pic_root / "ledger" / "pending_submission.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(marker["state"], "scheduler_job_id_received")
            with patch(
                "reconcile_frontier_job._scheduler_result",
                return_value=("CANCELLED", 0, 1),
            ):
                result = reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        self.assertEqual(result["state"], "CANCELLED")
        self.assertFalse((self.pic_root / "ledger" / "pending_submission.json").exists())

    def test_reconcile_received_job_id_rejects_scheduler_binding_mismatch(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=(
                f"JobId=12345 JobState=CANCELLED Account=AST207 "
                f"Comment=pic-reservation={reservation_id}"
            ),
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value="JobId=12345 JobState=CANCELLED Account=ast207 Comment=wrong",
        ):
            with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
                with self.assertRaises(ValueError):
                    reconcile(
                        job_id="12345",
                        ledger_jsonl=self.ledger,
                        ledger_csv=self.csv,
                        receipts_jsonl=self.receipts,
                        mirror_jsonl=self.mirror,
                        control_plane_dir=self.control_plane_dir,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                    )
        scheduler_result.assert_not_called()

    def test_successor_retry_requires_verified_prior_pair_before_marker_clear(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        handoff = self._create_test_terminal_recovery_handoff(successor)
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            with patch(
                "reconcile_frontier_job._scheduler_result",
                return_value=("CANCELLED", 0, 1),
            ):
                reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=successor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    terminal_recovery_handoff=handoff,
                )
        marker_path = self.pic_root / "ledger" / "pending_submission.json"
        marker_path.write_text(
            json.dumps({"reservation_id": reservation["reservation_id"]}),
            encoding="utf-8",
        )
        self.project_home_control_plane_dir.rename(
            self.project_home_control_plane_dir.with_name("missing-retry-prior-control-plane")
        )
        with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
            with self.assertRaises(ValueError):
                reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=successor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    terminal_recovery_handoff=handoff,
                )
        scheduler_result.assert_not_called()
        self.assertTrue(marker_path.is_file())

    def test_reconcile_terminal_retry_clears_stranded_marker_without_reaccounting(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        self._attach(reservation_id)
        with patch(
            "reconcile_frontier_job._scheduler_result",
            return_value=("COMPLETED", 300, 1),
        ):
            first = reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        marker_path = self.pic_root / "ledger" / "pending_submission.json"
        marker_path.write_text(
            json.dumps({"reservation_id": reservation_id}),
            encoding="utf-8",
        )
        with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
            second = reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        scheduler_result.assert_not_called()
        self.assertEqual(second["event_sha256"], first["event_sha256"])
        self.assertFalse(marker_path.exists())

    def test_reconcile_planner_receipt_is_derived_after_durable_mirror_append(
        self,
    ) -> None:
        _, _, artifact_dir = self._prepare_planner_retention_reconciliation()
        original_write = reconcile_frontier_job.atomic_write_bytes_at
        observed_post_mirror_authority = False

        def observe_post_mirror_authority(*args: object, **kwargs: object) -> None:
            nonlocal observed_post_mirror_authority
            records = ledger.validate_mirrored_state(
                self.ledger, self.receipts, self.mirror
            )
            self.assertEqual(records[-1]["event_type"], "reconciliation")
            receipts = ledger.validate_receipts(
                self.receipts,
                records,
                mirror_jsonl=self.mirror,
                mirror_transport="filesystem_copy",
            )
            self.assertEqual(
                receipts[-1]["mirrored_event_sha256"], records[-1]["event_sha256"]
            )
            observed_post_mirror_authority = True
            original_write(*args, **kwargs)

        with patch(
            "reconcile_frontier_job._scheduler_result",
            return_value=("COMPLETED", 60, 1),
        ), patch(
            "reconcile_frontier_job.atomic_write_bytes_at",
            side_effect=observe_post_mirror_authority,
        ):
            event = reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        receipt_path = (
            artifact_dir
            / "analysis"
            / reconcile_frontier_job.REGISTERED_EXECUTION_RECEIPT_NAME
        )
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        self.assertTrue(observed_post_mirror_authority)
        self.assertEqual(receipt["reconciliation_event_sha256"], event["event_sha256"])
        self.assertEqual(receipt["artifact_dir"], str(artifact_dir))
        self.assertEqual(receipt["planner_retention"], self._planner_retention())
        self.assertNotEqual(receipt["artifact_dir"], receipt["raw_output_root"])
        self.assertEqual(stat.S_IMODE(receipt_path.stat().st_mode), 0o444)

    def test_reconcile_planner_receipt_retry_is_idempotent(self) -> None:
        _, _, artifact_dir = self._prepare_planner_retention_reconciliation()
        with patch(
            "reconcile_frontier_job._scheduler_result",
            return_value=("COMPLETED", 60, 1),
        ):
            first = reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        receipt_path = (
            artifact_dir
            / "analysis"
            / reconcile_frontier_job.REGISTERED_EXECUTION_RECEIPT_NAME
        )
        receipt_bytes = receipt_path.read_bytes()
        with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
            second = reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        scheduler_result.assert_not_called()
        self.assertEqual(second["event_sha256"], first["event_sha256"])
        self.assertEqual(receipt_path.read_bytes(), receipt_bytes)

    def test_reconcile_planner_receipt_rejects_caller_selected_restart_namespace(
        self,
    ) -> None:
        retention = self._planner_retention()
        attempt_id = "q011-restart-fixture-001"
        attempt_root = (
            self.pic_root
            / "campaigns"
            / "q011-section54-fixture"
            / "restart_continuation"
            / attempt_id
        )
        retention.update(
            {
                "attempt_id": attempt_id,
                "authorized_orion_attempt_root": str(attempt_root),
                "authorized_orion_raw_root": str(attempt_root / "raw"),
                "argv": [
                    "-r",
                    "/retained/source.rst",
                    "-d",
                    str(attempt_root / "raw"),
                ],
            }
        )
        self._write_science_config(authorize=True, planner_retention=retention)
        with self.assertRaisesRegex(
            ValueError, "does not select one immutable descriptor"
        ):
            self._create_manifest()

    def test_reconcile_planner_receipt_rejects_analysis_directory_substitution(
        self,
    ) -> None:
        _, _, artifact_dir = self._prepare_planner_retention_reconciliation()
        analysis_dir = artifact_dir / "analysis"
        detached = artifact_dir / "analysis.detached"
        original_write = reconcile_frontier_job.atomic_write_bytes_at

        def substitute_then_write(*args: object, **kwargs: object) -> None:
            analysis_dir.rename(detached)
            analysis_dir.mkdir(mode=0o700)
            original_write(*args, **kwargs)

        with patch(
            "reconcile_frontier_job._scheduler_result",
            return_value=("COMPLETED", 60, 1),
        ), patch(
            "reconcile_frontier_job.atomic_write_bytes_at",
            side_effect=substitute_then_write,
        ), self.assertRaisesRegex(ValueError, "Directory ancestry changed"):
            reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        receipt_name = reconcile_frontier_job.REGISTERED_EXECUTION_RECEIPT_NAME
        self.assertFalse((analysis_dir / receipt_name).exists())
        self.assertFalse((detached / receipt_name).exists())
        analysis_dir.rmdir()
        detached.rename(analysis_dir)
        with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
            event = reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        scheduler_result.assert_not_called()
        self.assertEqual(event["event_type"], "reconciliation")
        self.assertTrue((analysis_dir / receipt_name).is_file())

    def test_reconcile_terminal_retry_rejects_stale_policy_generation(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        self._attach(reservation_id)
        with patch(
            "reconcile_frontier_job._scheduler_result",
            return_value=("COMPLETED", 300, 1),
        ):
            reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        self._promote_test_control_plane_successor(successor)
        marker_path = self.pic_root / "ledger" / "pending_submission.json"
        marker_path.write_text(
            json.dumps({"reservation_id": reservation_id}),
            encoding="utf-8",
        )
        with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
            with self.assertRaises(ValueError):
                reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        scheduler_result.assert_not_called()
        self.assertTrue(marker_path.is_file())

    def test_reconcile_requires_paired_installed_generation_before_scheduler_query(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        self._attach(str(reservation["reservation_id"]))
        self.project_home_control_plane_dir.rename(
            self.project_home_control_plane_dir.with_name("missing-control-plane")
        )
        with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
            with self.assertRaises(ValueError):
                reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        scheduler_result.assert_not_called()

    def test_successor_reconcile_accepts_verified_prior_installed_generation(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
            successor = self._publish_test_control_plane_successor(self.pic_root)
            project_home_successor = self._publish_test_control_plane_successor(
                self.project_home_root
            )
            self.assertEqual(successor.name, project_home_successor.name)
            handoff = self._create_test_terminal_recovery_handoff(successor)
            with patch(
                "reconcile_frontier_job._scheduler_result",
                return_value=("CANCELLED", 0, 1),
            ):
                result = reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=successor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    terminal_recovery_handoff=handoff,
                )
        self.assertEqual(result["state"], "CANCELLED")
        self.assertEqual(result["control_plane_version"], self.control_plane_version)
        self.assertEqual(result["reconciled_by_control_plane_version"], successor.name)
        records = validate_primary_chain(self.ledger)
        self.assertEqual(records[-2]["attached_by_control_plane_version"], successor.name)
        self.assertEqual(records[-1]["terminal_recovery_handoff_path"], str(handoff))
        self.assertIn("terminal_recovery_handoff_sha256", self.csv.read_text())
        self._promote_test_control_plane_successor(successor)

    def test_successor_reconcile_accepts_authorized_purged_zero_execution_cancellation(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        with self.assertRaises(ValueError):
            create_handoff(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                authorize_purged_cancelled_zero_execution=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with self.assertRaises(ValueError):
            create_handoff(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                attest_reviewed_purged_reservation_job_binding=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        full_row = (
            "12345|run_installed_control_plane_job.sh|CANCELLED by 18664|0|0||"
            "ast207|2026-05-30T15:47:31|None|2026-05-30T15:47:31|0:0\n"
        )
        short_row = "12345|CANCELLED by 18664|0|0||ast207\n"

        def scheduler_output(command: list[str], **kwargs: object) -> str:
            if command[0] == TRUSTED_SCONTROL:
                raise subprocess.CalledProcessError(
                    1,
                    command,
                    stderr="slurm_load_jobs error: Invalid job id specified\n",
                )
            if command[0] == TRUSTED_SQUEUE:
                return ""
            self.assertEqual(command[0], TRUSTED_SACCT)
            return short_row if "JobIDRaw,State,ElapsedRaw" in command[-1] else full_row

        with patch.object(
            terminal_recovery_handoff.subprocess,
            "check_output",
            side_effect=scheduler_output,
        ):
            handoff = self._create_test_terminal_recovery_handoff(
                successor,
                authorize_purged_cancelled_zero_execution=True,
            )
            result = reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                terminal_recovery_handoff=handoff,
            )
        self.assertEqual(result["state"], "CANCELLED")
        self.assertEqual(
            result["terminal_recovery_mode"],
            terminal_recovery_handoff.PURGED_CANCELLED_ZERO_EXECUTION_MODE,
        )
        self.assertFalse((self.pic_root / "ledger" / "pending_submission.json").exists())
        records = validate_primary_chain(self.ledger)
        self.assertEqual(records[-2]["terminal_recovery_mode"], result["terminal_recovery_mode"])
        self.assertIn("terminal_recovery_mode", self.csv.read_text())

    def test_successor_purged_reconcile_retry_uses_frozen_scheduler_snapshot(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        full_row = (
            "12345|run_installed_control_plane_job.sh|CANCELLED by 18664|0|0||"
            "ast207|2026-05-30T15:47:31|None|2026-05-30T15:47:31|0:0\n"
        )
        short_row = "12345|CANCELLED by 18664|0|0||ast207\n"

        def scheduler_output(command: list[str], **kwargs: object) -> str:
            if command[0] == TRUSTED_SCONTROL:
                raise subprocess.CalledProcessError(
                    1,
                    command,
                    stderr="slurm_load_jobs error: Invalid job id specified\n",
                )
            if command[0] == TRUSTED_SQUEUE:
                return ""
            self.assertEqual(command[0], TRUSTED_SACCT)
            return short_row if "JobIDRaw,State,ElapsedRaw" in command[-1] else full_row

        with patch.object(
            terminal_recovery_handoff.subprocess,
            "check_output",
            side_effect=scheduler_output,
        ):
            handoff = self._create_test_terminal_recovery_handoff(
                successor,
                authorize_purged_cancelled_zero_execution=True,
            )
            with patch(
                "reconcile_frontier_job._clear_matching_pending_marker",
                side_effect=RuntimeError("simulated marker cleanup interruption"),
            ):
                with self.assertRaises(RuntimeError):
                    reconcile(
                        job_id="12345",
                        ledger_jsonl=self.ledger,
                        ledger_csv=self.csv,
                        receipts_jsonl=self.receipts,
                        mirror_jsonl=self.mirror,
                        control_plane_dir=successor,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                        terminal_recovery_handoff=handoff,
                    )
        marker_path = self.pic_root / "ledger" / "pending_submission.json"
        self.assertTrue(marker_path.is_file())
        interrupted_records = validate_primary_chain(self.ledger)
        with patch.object(
            terminal_recovery_handoff.subprocess,
            "check_output",
            side_effect=AssertionError("retry must use the frozen scheduler snapshot"),
        ):
            with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
                result = reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=successor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    terminal_recovery_handoff=handoff,
                )
        scheduler_result.assert_not_called()
        self.assertEqual(result["event_type"], "reconciliation")
        self.assertFalse(marker_path.exists())
        self.assertEqual(validate_primary_chain(self.ledger), interrupted_records)

    def test_successor_purged_reconcile_resumes_after_terminal_append_failure(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        full_row = (
            "12345|run_installed_control_plane_job.sh|CANCELLED by 18664|0|0||"
            "ast207|2026-05-30T15:47:31|None|2026-05-30T15:47:31|0:0\n"
        )
        short_row = "12345|CANCELLED by 18664|0|0||ast207\n"

        def scheduler_output(command: list[str], **kwargs: object) -> str:
            if command[0] == TRUSTED_SCONTROL:
                raise subprocess.CalledProcessError(
                    1,
                    command,
                    stderr="slurm_load_jobs error: Invalid job id specified\n",
                )
            if command[0] == TRUSTED_SQUEUE:
                return ""
            self.assertEqual(command[0], TRUSTED_SACCT)
            return short_row if "JobIDRaw,State,ElapsedRaw" in command[-1] else full_row

        real_append = reconcile_frontier_job.append_primary_event_locked
        append_count = 0

        def fail_second_append(*args: object, **kwargs: object) -> dict[str, object]:
            nonlocal append_count
            append_count += 1
            if append_count == 2:
                raise RuntimeError("simulated terminal reconciliation interruption")
            return real_append(*args, **kwargs)

        with patch.object(
            terminal_recovery_handoff.subprocess,
            "check_output",
            side_effect=scheduler_output,
        ):
            handoff = self._create_test_terminal_recovery_handoff(
                successor,
                authorize_purged_cancelled_zero_execution=True,
            )
            with patch(
                "reconcile_frontier_job.append_primary_event_locked",
                side_effect=fail_second_append,
            ):
                with self.assertRaises(RuntimeError):
                    reconcile(
                        job_id="12345",
                        ledger_jsonl=self.ledger,
                        ledger_csv=self.csv,
                        receipts_jsonl=self.receipts,
                        mirror_jsonl=self.mirror,
                        control_plane_dir=successor,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                        terminal_recovery_handoff=handoff,
                    )
        marker_path = self.pic_root / "ledger" / "pending_submission.json"
        self.assertTrue(marker_path.is_file())
        interrupted_records = validate_primary_chain(self.ledger)
        self.assertEqual(interrupted_records[-1]["event_type"], "job_id_attached")
        with patch.object(
            terminal_recovery_handoff.subprocess,
            "check_output",
            side_effect=AssertionError("retry must use the frozen scheduler snapshot"),
        ):
            with patch(
                "reconcile_frontier_job._scheduler_result",
                side_effect=AssertionError("retry must not query scheduler accounting"),
            ):
                result = reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=successor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    terminal_recovery_handoff=handoff,
                )
        self.assertEqual(result["event_type"], "reconciliation")
        self.assertFalse(marker_path.exists())
        records = validate_primary_chain(self.ledger)
        self.assertEqual(
            [record["event_type"] for record in records].count("job_id_attached"),
            1,
        )
        self.assertEqual(
            [record["event_type"] for record in records].count("reconciliation"),
            1,
        )

    def test_successor_reconcile_resumes_after_attachment_append_failure(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        handoff = self._create_test_terminal_recovery_handoff(successor)
        real_append = reconcile_frontier_job.append_primary_event_locked
        append_count = 0

        def fail_second_append(*args: object, **kwargs: object) -> dict[str, object]:
            nonlocal append_count
            append_count += 1
            if append_count == 2:
                raise RuntimeError("simulated terminal reconciliation interruption")
            return real_append(*args, **kwargs)

        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            with patch(
                "reconcile_frontier_job._scheduler_result",
                return_value=("CANCELLED", 0, 1),
            ):
                with patch(
                    "reconcile_frontier_job.append_primary_event_locked",
                    side_effect=fail_second_append,
                ):
                    with self.assertRaises(RuntimeError):
                        reconcile(
                            job_id="12345",
                            ledger_jsonl=self.ledger,
                            ledger_csv=self.csv,
                            receipts_jsonl=self.receipts,
                            mirror_jsonl=self.mirror,
                            control_plane_dir=successor,
                            authorized_pic_root=self.pic_root,
                            authorized_project_home_root=self.project_home_root,
                            terminal_recovery_handoff=handoff,
                        )
        marker_path = self.pic_root / "ledger" / "pending_submission.json"
        self.assertTrue(marker_path.is_file())
        interrupted_records = validate_primary_chain(self.ledger)
        self.assertEqual(
            [record["event_type"] for record in interrupted_records].count(
                "job_id_attached"
            ),
            1,
        )
        with self.assertRaises(ValueError):
            repair_reservation_attachments(
                reservation_id=reservation_id,
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertTrue(marker_path.is_file())
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            with patch(
                "reconcile_frontier_job._scheduler_result",
                return_value=("CANCELLED", 0, 1),
            ):
                result = reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=successor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    terminal_recovery_handoff=handoff,
                )
        self.assertEqual(result["state"], "CANCELLED")
        self.assertFalse(marker_path.exists())
        records = validate_primary_chain(self.ledger)
        self.assertEqual(
            [record["event_type"] for record in records].count("job_id_attached"),
            1,
        )
        self.assertEqual(
            [record["event_type"] for record in records].count("reconciliation"),
            1,
        )

    def test_successor_reconcile_retry_clears_marker_after_terminal_append(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        handoff = self._create_test_terminal_recovery_handoff(successor)
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            with patch(
                "reconcile_frontier_job._scheduler_result",
                return_value=("CANCELLED", 0, 1),
            ):
                with patch(
                    "reconcile_frontier_job._clear_matching_pending_marker",
                    side_effect=RuntimeError("simulated marker cleanup interruption"),
                ):
                    with self.assertRaises(RuntimeError):
                        reconcile(
                            job_id="12345",
                            ledger_jsonl=self.ledger,
                            ledger_csv=self.csv,
                            receipts_jsonl=self.receipts,
                            mirror_jsonl=self.mirror,
                            control_plane_dir=successor,
                            authorized_pic_root=self.pic_root,
                            authorized_project_home_root=self.project_home_root,
                            terminal_recovery_handoff=handoff,
                        )
        marker_path = self.pic_root / "ledger" / "pending_submission.json"
        self.assertTrue(marker_path.is_file())
        interrupted_records = validate_primary_chain(self.ledger)
        with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
            result = reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                terminal_recovery_handoff=handoff,
            )
        scheduler_result.assert_not_called()
        self.assertFalse(marker_path.exists())
        self.assertEqual(result["event_type"], "reconciliation")
        self.assertEqual(validate_primary_chain(self.ledger), interrupted_records)

    def test_successor_reconcile_rejects_missing_prior_pair_before_scheduler_query(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        self._attach(str(reservation["reservation_id"]))
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        self.project_home_control_plane_dir.rename(
            self.project_home_control_plane_dir.with_name("missing-prior-control-plane")
        )
        with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
            with self.assertRaises(ValueError):
                reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=successor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                    terminal_recovery_handoff=(
                        self.pic_root / "policy" / "recovery_handoffs" / f"{uuid.uuid4()}.json"
                    ),
                )
        scheduler_result.assert_not_called()

    def test_successor_reconcile_requires_explicit_prior_version_before_scheduler_query(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
            with self.assertRaises(ValueError):
                reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=successor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        scheduler_result.assert_not_called()

    def test_policy_promotion_rejects_pending_terminal_recovery(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        with self.assertRaises(ValueError):
            self._promote_test_control_plane_successor(successor)

    def test_successor_reconcile_rejects_received_marker_field_mutation(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        handoff = self._create_test_terminal_recovery_handoff(successor)
        marker_path = self.pic_root / "ledger" / "pending_submission.json"
        canonical_marker = json.loads(marker_path.read_text(encoding="utf-8"))
        records_before = validate_primary_chain(self.ledger)
        marker_path.chmod(0o600)
        mutations = [
            ("schema_version", 3),
            ("state", "submitted_not_attached"),
            ("reservation_id", "wrong"),
            ("submission_id", "wrong"),
            ("manifest_path", "wrong"),
            ("manifest_sha256", "0" * 64),
            ("job_id", "wrong"),
            ("extra", "wrong"),
        ]
        for key, value in mutations:
            marker = dict(canonical_marker)
            marker[key] = value
            marker_path.write_text(json.dumps(marker), encoding="utf-8")
            with patch("reconcile_frontier_job._scheduler_result") as scheduler_result:
                with self.assertRaises(ValueError):
                    reconcile(
                        job_id="12345",
                        ledger_jsonl=self.ledger,
                        ledger_csv=self.csv,
                        receipts_jsonl=self.receipts,
                        mirror_jsonl=self.mirror,
                        control_plane_dir=successor,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                        terminal_recovery_handoff=handoff,
                    )
            scheduler_result.assert_not_called()
        marker_path.write_text(json.dumps(canonical_marker), encoding="utf-8")
        marker_path.chmod(0o400)
        self.assertEqual(validate_primary_chain(self.ledger), records_before)

    def test_successor_only_allows_terminal_reconciliation_for_prior_generation(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        records_before = validate_primary_chain(self.ledger)
        with self.assertRaises(ValueError):
            mark_dispatch_started(
                reservation_id=reservation_id,
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with self.assertRaises(ValueError):
            reservation_bound_manifest(
                manifest_path,
                reservation_id,
                ledger_jsonl=self.ledger,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                require_reserved=True,
            )
        self.assertEqual(validate_primary_chain(self.ledger), records_before)

    def test_successor_reconcile_rejects_active_prior_job_without_append(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        reservation_id = str(reservation["reservation_id"])
        scheduler = (
            f"JobId=12345 JobState=CANCELLED Account=ast207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            self._mark_dispatch_started(reservation_id)
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        handoff = self._create_test_terminal_recovery_handoff(successor)
        records_before = validate_primary_chain(self.ledger)
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            with patch(
                "reconcile_frontier_job._scheduler_result",
                return_value=("RUNNING", 0, 1),
            ):
                with self.assertRaises(ValueError):
                    reconcile(
                        job_id="12345",
                        ledger_jsonl=self.ledger,
                        ledger_csv=self.csv,
                        receipts_jsonl=self.receipts,
                        mirror_jsonl=self.mirror,
                        control_plane_dir=successor,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                        terminal_recovery_handoff=handoff,
                    )
        self.assertEqual(validate_primary_chain(self.ledger), records_before)

    def test_registered_science_rejects_pending_clean_candidate_freeze(self) -> None:
        self._write_science_config(authorize=False)
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_registered_science_accepts_exact_authorized_clean_candidate(self) -> None:
        candidate = self._write_science_config(authorize=True)
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        candidate_sha256 = str(reservation["clean_candidate_manifest_sha256"])
        self.assertEqual(len(candidate_sha256), 64)
        self.assertEqual(reservation["submission_scope"], "registered_science")
        prepared = json.loads(candidate.read_text(encoding="utf-8"))["prepared_artifacts"]
        self.assertEqual(prepared["inventory_path"], self._prepared_artifact_inventory())
        self.assertEqual(len(prepared["paper_decks"]), 2)
        self.assertEqual(len(prepared["analyzers"]), 1)
        self.assertIn("clean_candidate_manifest_sha256", self.csv.read_text())

    def test_q011_first_pressure_case_accepts_empty_ordered_closure_list(self) -> None:
        _, authorization_id, campaign, test_id = (
            validate_and_reserve_frontier_job.Q011_PRESSURE_CASES[0]
        )
        self._write_science_config(
            authorize=True,
            registered_science_authorization_id=authorization_id,
            campaign=campaign,
            test_id=test_id,
            artifact_dir=str(self.pic_root / "runs" / campaign / self.submission_id),
            planner_retention=self._planner_retention(),
        )
        manifest_path = self._create_manifest()
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["prior_case_closures"], [])
        with patch(
            "validate_and_reserve_frontier_job._verify_clean_candidate",
            return_value=sha256(self._planner_candidate_manifest),
        ), patch(
            "validate_and_reserve_frontier_job._verify_q011_repaired_clean_candidate"
        ), patch(
            "validate_and_reserve_frontier_job._q011_selected_pressure",
            return_value=validate_and_reserve_frontier_job.Q011_PRESSURE_BY_AUTHORIZATION[
                authorization_id
            ],
        ):
            reservation = self._reserve(manifest_path)
        self.assertEqual(
            reservation["registered_science_authorization_id"], authorization_id
        )

    def test_q011_later_pressure_case_rejects_handwritten_missing_closures(self) -> None:
        _, authorization_id, campaign, test_id = (
            validate_and_reserve_frontier_job.Q011_PRESSURE_CASES[1]
        )
        self._write_science_config(
            authorize=True,
            registered_science_authorization_id=authorization_id,
            campaign=campaign,
            test_id=test_id,
            artifact_dir=str(self.pic_root / "runs" / campaign / self.submission_id),
            planner_retention=self._planner_retention(),
        )
        manifest_path = self._create_manifest()
        with patch(
            "validate_and_reserve_frontier_job._verify_clean_candidate",
            return_value=sha256(self._planner_candidate_manifest),
        ), patch(
            "validate_and_reserve_frontier_job._verify_q011_repaired_clean_candidate"
        ), patch(
            "validate_and_reserve_frontier_job._q011_selected_pressure",
            return_value=validate_and_reserve_frontier_job.Q011_PRESSURE_BY_AUTHORIZATION[
                authorization_id
            ],
        ), self.assertRaisesRegex(ValueError, "exact ordered predecessor closures"):
            self._reserve(manifest_path)

    def test_q011_pressure_retry_rejects_failed_v1_and_unrepaired_carriers(self) -> None:
        _, authorization_id, campaign, test_id = (
            validate_and_reserve_frontier_job.Q011_PRESSURE_CASES[0]
        )
        repository = install_control_plane.SCRIPT_DIR.parents[2]
        manifest = {
            "registered_science_authorization_id": authorization_id,
            "campaign": campaign,
            "test_id": test_id,
            "launch_contract": json.loads(
                (
                    repository
                    / "tst/publication/readiness/"
                    "frontier_q011_section54_pressure_ps_p0_1p00_launch_contract.json"
                ).read_text(encoding="utf-8")
            ),
            "snapshot_files": [
                {
                    "role": "job-script",
                    "sha256": validate_and_reserve_frontier_job.Q011_JOB_SCRIPT_SHA256,
                },
                {
                    "role": "input-deck",
                    "sha256": validate_and_reserve_frontier_job.Q011_INPUT_DECK_SHA256,
                },
            ],
        }
        candidate = {"source": {"git_commit": "a" * 40}}
        with self.assertRaisesRegex(ValueError, "failed v1"):
            validate_and_reserve_frontier_job._verify_q011_repaired_clean_candidate(
                manifest,
                {"source": {"git_commit": validate_and_reserve_frontier_job.Q011_FAILED_V1_GIT_COMMIT}},
                candidate_sha256="a" * 64,
                source_archive=b"",
                executable_sha256="b" * 64,
            )
        with self.assertRaisesRegex(ValueError, "failed v1"):
            validate_and_reserve_frontier_job._verify_q011_repaired_clean_candidate(
                manifest,
                candidate,
                candidate_sha256=validate_and_reserve_frontier_job.Q011_FAILED_V1_CLEAN_CANDIDATE_MANIFEST_SHA256,
                source_archive=b"",
                executable_sha256="b" * 64,
            )
        with self.assertRaisesRegex(ValueError, "failed v1"):
            validate_and_reserve_frontier_job._verify_q011_repaired_clean_candidate(
                manifest,
                candidate,
                candidate_sha256="a" * 64,
                source_archive=b"",
                executable_sha256=validate_and_reserve_frontier_job.Q011_FAILED_V1_EXECUTABLE_SHA256,
            )
        archive_payload = io.BytesIO()
        with tarfile.open(fileobj=archive_payload, mode="w") as archive:
            payload = b"unrepaired\n"
            member = tarfile.TarInfo(
                validate_and_reserve_frontier_job.Q011_REPAIRED_GENERATOR_PATH
            )
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
        with self.assertRaisesRegex(ValueError, "repaired generator"):
            validate_and_reserve_frontier_job._verify_q011_repaired_clean_candidate(
                manifest,
                candidate,
                candidate_sha256="a" * 64,
                source_archive=archive_payload.getvalue(),
                executable_sha256="b" * 64,
            )

    def test_q011_pressure_retry_rejects_equivalent_authorization_alias(self) -> None:
        _, _, campaign, test_id = validate_and_reserve_frontier_job.Q011_PRESSURE_CASES[0]
        with self.assertRaisesRegex(ValueError, "aliases are forbidden"):
            validate_and_reserve_frontier_job._verify_q011_repaired_clean_candidate(
                {
                    "registered_science_authorization_id": "q011-equivalent-alias",
                    "campaign": campaign,
                    "test_id": test_id,
                    "launch_contract": {},
                    "snapshot_files": [],
                },
                {"source": {"git_commit": "a" * 40}},
                candidate_sha256="a" * 64,
                source_archive=b"",
                executable_sha256="b" * 64,
            )

    def test_q011_snapshot_analyzer_rejects_inherited_helper_override(self) -> None:
        with patch.dict(os.environ, {"PIC_F1_ANALYSIS_HELPER_FD": "7"}):
            with self.assertRaisesRegex(ValueError, "inherited helper overrides"):
                validate_and_reserve_frontier_job._q011_snapshot_analyzer({}, self.pic_root)

    def _q011_snapshot_analyzer_fixture(
        self,
    ) -> tuple[dict[str, object], Path, Path]:
        snapshot_root = self.root / "q011-snapshot" / "snapshot"
        analysis_root = snapshot_root / "analysis"
        analysis_root.mkdir(parents=True)
        repository = install_control_plane.SCRIPT_DIR.parents[2]
        records = []
        for role, name in (
            ("analysis-script-000", "analyze_q011_section54_pressure_pilot_case.py"),
            ("analysis-script-001", "frontier_f1_structured_artifacts.py"),
        ):
            source = repository / "tst" / "publication" / name
            target = analysis_root / (
                f"000-{name}" if role == "analysis-script-000" else name
            )
            shutil.copyfile(source, target)
            target.chmod(0o444)
            records.append({"role": role, "path": str(target), "sha256": sha256(target)})
        return {"snapshot_files": records}, snapshot_root.parent / "pre_submit_manifest.json", analysis_root

    def test_q011_snapshot_analyzer_executes_verified_source_without_bytecode(self) -> None:
        manifest, manifest_path, _ = self._q011_snapshot_analyzer_fixture()
        with patch.multiple(
            validate_and_reserve_frontier_job,
            Q011_RAW_ANALYZER_SHA256=record_for_role(
                manifest, "analysis-script-000"
            )["sha256"],
            Q011_STRUCTURED_HELPER_SHA256=record_for_role(
                manifest, "analysis-script-001"
            )["sha256"],
        ):
            module = validate_and_reserve_frontier_job._q011_snapshot_analyzer(
                manifest, manifest_path
            )
        self.assertTrue(callable(module.verify_published_case_descriptor))

    def test_q011_snapshot_analyzer_rejects_adjacent_cached_bytecode(self) -> None:
        manifest, manifest_path, analysis_root = self._q011_snapshot_analyzer_fixture()
        cache = analysis_root / "__pycache__"
        cache.mkdir()
        (cache / "injected.cpython-311.pyc").write_bytes(b"untrusted bytecode\n")
        with patch.multiple(
            validate_and_reserve_frontier_job,
            Q011_RAW_ANALYZER_SHA256=record_for_role(
                manifest, "analysis-script-000"
            )["sha256"],
            Q011_STRUCTURED_HELPER_SHA256=record_for_role(
                manifest, "analysis-script-001"
            )["sha256"],
        ), self.assertRaisesRegex(ValueError, "registered bytes"):
            validate_and_reserve_frontier_job._q011_snapshot_analyzer(
                manifest, manifest_path
            )

    def test_q011_snapshot_analyzer_rejects_post_hash_helper_inode_mutation(self) -> None:
        manifest, manifest_path, _ = self._q011_snapshot_analyzer_fixture()
        original = validate_and_reserve_frontier_job._q011_open_read_only_source

        def mutate_after_hash(path: Path, *, label: str) -> tuple[int, bytes]:
            descriptor, payload = original(path, label=label)
            if label == "Q011 predecessor helper snapshot":
                path.chmod(0o644)
                path.write_text(
                    "raise RuntimeError('MUTATED_HELPER_EXECUTED')\n",
                    encoding="utf-8",
                )
                path.chmod(0o444)
            return descriptor, payload

        with patch(
            "validate_and_reserve_frontier_job._q011_open_read_only_source",
            side_effect=mutate_after_hash,
        ), patch.multiple(
            validate_and_reserve_frontier_job,
            Q011_RAW_ANALYZER_SHA256=record_for_role(
                manifest, "analysis-script-000"
            )["sha256"],
            Q011_STRUCTURED_HELPER_SHA256=record_for_role(
                manifest, "analysis-script-001"
            )["sha256"],
        ), self.assertRaisesRegex(ValueError, "differs from registered bytes"):
            validate_and_reserve_frontier_job._q011_snapshot_analyzer(
                manifest, manifest_path
            )

    def test_q011_second_case_accepts_real_recomputed_predecessor_descriptor(self) -> None:
        from tst.publication import analyze_q011_section54_pressure_pilot_case as analyzer
        from tst.publication.test_publish_q011_section54_pressure_pilot_bundle import (
            _descriptor_sha256,
            _make_writable,
            _raw_tree,
        )

        first_case, first_authorization, first_campaign, first_test = (
            validate_and_reserve_frontier_job.Q011_PRESSURE_CASES[0]
        )
        _, second_authorization, second_campaign, second_test = (
            validate_and_reserve_frontier_job.Q011_PRESSURE_CASES[1]
        )
        first_submission = str(uuid.uuid4())
        artifact_dir = self.pic_root / "runs" / first_campaign / first_submission
        artifact_dir.parent.mkdir(parents=True)
        _raw_tree(artifact_dir, first_case)
        descriptor = analyzer.publish_case_descriptor(artifact_dir, first_case)
        descriptor_sha256 = _descriptor_sha256(descriptor)
        try:
            prior_manifest_path = (
                self.pic_root
                / "manifests"
                / first_campaign
                / first_submission
                / "pre_submit_manifest.json"
            )
            analysis_root = prior_manifest_path.parent / "snapshot" / "analysis"
            analysis_root.mkdir(parents=True)
            repository = install_control_plane.SCRIPT_DIR.parents[2]
            analyzer_path = (
                analysis_root / "000-analyze_q011_section54_pressure_pilot_case.py"
            )
            helper_path = analysis_root / "frontier_f1_structured_artifacts.py"
            shutil.copyfile(
                repository / "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
                analyzer_path,
            )
            shutil.copyfile(
                repository / "tst/publication/frontier_f1_structured_artifacts.py",
                helper_path,
            )
            analyzer_path.chmod(0o444)
            helper_path.chmod(0o444)
            candidate_sha256 = "c" * 64
            executable_sha256 = "e" * 64
            git_commit = "d" * 40
            prior_manifest = {
                "clean_candidate_manifest_sha256": candidate_sha256,
                "git_commit": git_commit,
                "snapshot_files": [
                    {
                        "role": "executable",
                        "path": str(self.pic_root / "fixture-athena"),
                        "sha256": executable_sha256,
                    },
                    {
                        "role": "analysis-script-000",
                        "path": str(analyzer_path),
                        "sha256": sha256(analyzer_path),
                    },
                    {
                        "role": "analysis-script-001",
                        "path": str(helper_path),
                        "sha256": sha256(helper_path),
                    },
                ],
            }
            prior_manifest_path.write_text(
                json.dumps(prior_manifest), encoding="utf-8"
            )
            prior_manifest_path.chmod(0o444)
            event_sha256 = "f" * 64
            closure = {
                "case_id": first_case,
                "submission_id": first_submission,
                "artifact_dir": str(artifact_dir),
                "descriptor_path": str(artifact_dir / "analysis" / "analysis.json"),
                "descriptor_sha256": descriptor_sha256,
                "reconciliation_event_sha256": event_sha256,
            }
            launch_contract = json.loads(
                (
                    repository
                    / "tst/publication/readiness/"
                    "frontier_q011_section54_pressure_ps_p0_0p05_launch_contract.json"
                ).read_text(encoding="utf-8")
            )
            manifest = {
                "registered_science_authorization_id": second_authorization,
                "campaign": second_campaign,
                "test_id": second_test,
                "launch_contract": launch_contract,
                "clean_candidate_manifest_sha256": candidate_sha256,
                "git_commit": git_commit,
                "prior_case_closures": [closure],
                "planner_retention": self._planner_retention(),
                "snapshot_files": [
                    {
                        "role": "job-script",
                        "sha256": validate_and_reserve_frontier_job.Q011_JOB_SCRIPT_SHA256,
                    },
                    {
                        "role": "input-deck",
                        "sha256": validate_and_reserve_frontier_job.Q011_INPUT_DECK_SHA256,
                    },
                    {
                        "role": "executable",
                        "sha256": executable_sha256,
                    },
                ],
            }
            records = [
                {
                    "event_sha256": event_sha256,
                    "event_type": "reconciliation",
                    "submission_scope": "registered_science",
                    "submission_id": first_submission,
                    "campaign": first_campaign,
                    "test_id": first_test,
                    "registered_science_authorization_id": first_authorization,
                    "artifact_dir": str(artifact_dir),
                    "clean_candidate_manifest_sha256": candidate_sha256,
                    "git_commit": git_commit,
                    "executable_sha256": executable_sha256,
                    "reconciled": True,
                    "state": "COMPLETED",
                    "manifest_path": str(prior_manifest_path),
                    "manifest_sha256": sha256(prior_manifest_path),
                }
            ]
            with patch.multiple(
                validate_and_reserve_frontier_job,
                Q011_RAW_ANALYZER_SHA256=sha256(analyzer_path),
                Q011_STRUCTURED_HELPER_SHA256=sha256(helper_path),
            ):
                validate_and_reserve_frontier_job._verify_q011_prior_case_closures(
                    manifest, records, authorized_pic_root=self.pic_root
                )
            records[0]["manifest_sha256"] = "0" * 64
            with self.assertRaisesRegex(ValueError, "manifest digest drifted"):
                validate_and_reserve_frontier_job._verify_q011_prior_case_closures(
                    manifest, records, authorized_pic_root=self.pic_root
                )
        finally:
            _make_writable(artifact_dir)

    def test_registered_science_policy_promotion_requires_pre_policy_attestation(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        with self.assertRaisesRegex(
            ValueError, "requires a sealed pre-policy-promotion attestation"
        ):
            promote(
                self.policy,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_registered_science_policy_promotion_retains_attestation_binding(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        promotion = json.loads(
            (self.pic_root / "policy" / "active_promotion.json").read_text(
                encoding="utf-8"
            )
        )
        path = Path(str(promotion["pre_policy_promotion_attestation_path"]))
        self.assertTrue(path.is_file())
        self.assertEqual(
            promotion["pre_policy_promotion_attestation_sha256"], sha256(path)
        )
        _, snapshot = require_storage_policy_unlock_snapshot(
            control_plane_version=self.control_plane_version,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        self.assertRegex(snapshot["active_promotion_sha256"], r"^[0-9a-f]{64}$")

    def test_registered_science_policy_unlock_revalidates_retained_attestation_tree(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        promotion = json.loads(
            (self.pic_root / "policy" / "active_promotion.json").read_text(
                encoding="utf-8"
            )
        )
        path = Path(str(promotion["pre_policy_promotion_attestation_path"]))
        member = path.parent / "same_account_process_snapshot.txt"
        path.parent.chmod(0o700)
        member.chmod(0o600)
        member.write_bytes(member.read_bytes() + b"forged\n")
        member.chmod(0o400)
        path.parent.chmod(0o500)
        with self.assertRaisesRegex(ValueError, "checksum differs"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_registered_science_policy_promotion_revalidates_attestation_inside_lock(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        arguments = self._pre_policy_promotion_attestation_arguments()
        path = Path(str(arguments["pre_policy_promotion_attestation"]))
        original_lock = promote_active_policy._promotion_lock

        @contextmanager
        def tampering_lock(root: Path):
            with original_lock(root) as descriptor:
                member = path.parent / "same_account_process_snapshot.txt"
                path.parent.chmod(0o700)
                member.chmod(0o600)
                member.write_bytes(member.read_bytes() + b"forged\n")
                member.chmod(0o400)
                path.parent.chmod(0o500)
                yield descriptor

        with patch.object(
            promote_active_policy, "_promotion_lock", tampering_lock
        ), self.assertRaisesRegex(ValueError, "checksum differs"):
            promote(
                self.policy,
                **arguments,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_registered_science_policy_promotion_revalidates_attestation_after_mirror_lock(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        arguments = self._pre_policy_promotion_attestation_arguments()
        path = Path(str(arguments["pre_policy_promotion_attestation"]))
        original_open = promote_active_policy.open_directory_below
        tampered = False

        def tampering_open(directory: Path, *, root: Path) -> int:
            nonlocal tampered
            descriptor = original_open(directory, root=root)
            if (
                not tampered
                and Path(os.path.abspath(directory))
                == self.project_home_root / "policy"
            ):
                tampered = True
                member = path.parent / "same_account_process_snapshot.txt"
                path.parent.chmod(0o700)
                member.chmod(0o600)
                member.write_bytes(member.read_bytes() + b"forged\n")
                member.chmod(0o400)
                path.parent.chmod(0o500)
            return descriptor

        with patch.object(
            promote_active_policy, "open_directory_below", tampering_open
        ), self.assertRaisesRegex(ValueError, "checksum differs"):
            promote(
                self.policy,
                **arguments,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertTrue(tampered)

    def test_registered_science_rejects_missing_wrapper_attestation(self) -> None:
        self._write_science_config(authorize=True)
        manifest_path = self._create_manifest()
        with self.assertRaisesRegex(
            ValueError, "requires a sealed pre-submit-wrapper attestation"
        ):
            self._reserve(manifest_path, inject_wrapper_attestation=False)

    def test_registered_science_rejects_legacy_two_field_manifest_attestation(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        config = json.loads(self.config.read_text(encoding="utf-8"))
        path = self._sealed_operator_attestation(
            str(config["registered_science_authorization_id"]),
            "pre_manifest",
        )
        path.parent.chmod(0o700)
        path.chmod(0o600)
        path.write_text(
            json.dumps(
                {
                    "phase": "pre_manifest",
                    "registered_science_authorization_id": (
                        config["registered_science_authorization_id"]
                    ),
                }
            )
            + "\n",
            encoding="utf-8",
        )
        path.chmod(0o400)
        path.parent.chmod(0o500)
        config["pre_manifest_attestation"] = str(path)
        self.config.write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "root schema"):
            create_manifest(
                self.config,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
            )

    def test_registered_science_rejects_prepared_artifact_manifest_tamper(self) -> None:
        candidate = self._write_science_config(authorize=True)
        value = json.loads(candidate.read_text(encoding="utf-8"))
        value["prepared_artifacts"]["paper_decks"][0]["sha256"] = "0" * 64
        candidate.parent.chmod(0o755)
        candidate.chmod(0o644)
        candidate.write_text(json.dumps(value), encoding="utf-8")
        candidate.chmod(0o444)
        candidate.parent.chmod(0o555)
        self._update_registered_science_candidate_sha(candidate)
        self._write_policy(
            science_submission_freeze=self._authorized_science_freeze(candidate),
            admission_smoke_overrides={"status": "closed_after_pass"},
        )
        self._install_deliberately_malformed_active_policy_fixture()
        manifest_path = self._create_manifest()
        with self.assertRaisesRegex(ValueError, "differs from archived source inventory"):
            self._reserve(manifest_path)

    def test_registered_science_accepts_historical_candidate_receipt_after_successor(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        self._promote_test_control_plane_successor(successor)
        manifest_path = self._create_manifest(control_plane_dir=successor)
        reservation = self._reserve(manifest_path, control_plane_dir=successor)
        self.assertEqual(
            reservation["registered_science_authorization_id"], "f1-clean-gyro-v1"
        )

    def test_registered_science_rejects_unclosed_authorized_admission_smoke(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        self._write_policy(
            admission_smoke_overrides={
                "status": "authorized_f0_parser_contract_only"
            }
        )
        with self.assertRaisesRegex(ValueError, "closed admission smoke"):
            self._promote_policy()

    def test_registered_science_rejects_unclosed_pending_admission_smoke(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        self._write_policy(
            admission_smoke_overrides={"status": "pending_exact_executable_binding"}
        )
        with self.assertRaisesRegex(ValueError, "closed admission smoke"):
            self._promote_policy()

    def test_registered_science_rejects_unbound_build_profile_control_plane(
        self,
    ) -> None:
        candidate = self._write_science_config(authorize=True)
        self._write_policy(
            science_submission_freeze=self._authorized_science_freeze(
                candidate, build_profile_control_plane_version="0" * 64
            ),
            admission_smoke_overrides={"status": "closed_after_pass"},
        )
        self._install_deliberately_malformed_active_policy_fixture()
        manifest_path = self._create_manifest()
        with self.assertRaisesRegex(ValueError, "unauthorized control-plane version"):
            self._reserve(manifest_path)

    def test_registered_science_rejects_missing_paired_historical_build_authority(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        predecessor = self.project_home_root / "control_plane" / self.control_plane_version
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(successor.name, project_home_successor.name)
        self._promote_test_control_plane_successor(successor)
        manifest_path = self._create_manifest(control_plane_dir=successor)
        predecessor.rename(predecessor.with_name(f"{predecessor.name}.missing"))
        with self.assertRaisesRegex(ValueError, "Missing historical"):
            self._reserve(manifest_path, control_plane_dir=successor)

    def test_registered_science_rejects_unknown_authorization_id(self) -> None:
        self._write_science_config(authorize=True)
        config = json.loads(self.config.read_text(encoding="utf-8"))
        config["registered_science_authorization_id"] = "unknown-slice"
        self.config.write_text(json.dumps(config), encoding="utf-8")
        manifest_path = self._create_manifest()
        with self.assertRaisesRegex(ValueError, "authorization ID is not active"):
            self._reserve(manifest_path)

    def test_registered_science_attempt_is_consumed_by_cancelled_reservation(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        reservation = self._reserve(self._create_manifest())
        transition(
            reservation_id=str(reservation["reservation_id"]),
            notes="cancel before scheduler submission",
            event_type="reservation_cancelled",
            state="cancelled",
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            receipts_jsonl=self.receipts,
            mirror_jsonl=self.mirror,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        with self.assertRaisesRegex(ValueError, "attempt ceiling"):
            self._reserve(
                self._fresh_submission_manifest(),
                reservation_id=str(uuid.uuid4()),
            )

    def test_registered_science_attempt_is_consumed_by_completed_reservation(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        reservation = self._reserve(self._create_manifest())
        self._attach(str(reservation["reservation_id"]))
        with patch(
            "reconcile_frontier_job._scheduler_result",
            return_value=("COMPLETED", 60, 1),
        ):
            reconcile(
                job_id="12345",
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with self.assertRaisesRegex(ValueError, "attempt ceiling"):
            self._reserve(
                self._fresh_submission_manifest(),
                reservation_id=str(uuid.uuid4()),
            )

    def test_registered_science_concurrent_first_attempt_admits_only_one(self) -> None:
        self._write_science_config(authorize=True)
        manifests = [self._create_manifest(), self._fresh_submission_manifest()]

        def validate_with_test_roots(
            candidate: dict[str, object], **kwargs: object
        ) -> list[dict[str, str]]:
            if self.authorized_clean_candidate_source_root is None:
                raise AssertionError("Clean-candidate test source root was not injected")
            return validate_clean_candidate_bundle(
                candidate,
                **kwargs,
                authorized_pic_root=self.pic_root,
                authorized_source_root=self.authorized_clean_candidate_source_root,
            )

        def reserve_one(index: int) -> dict[str, object] | ValueError:
            try:
                return self._reserve(
                    manifests[index],
                    reservation_id=str(uuid.uuid4()),
                    patch_clean_candidate_bundle=False,
                )
            except ValueError as error:
                return error

        with patch(
            "validate_and_reserve_frontier_job.validate_clean_candidate_bundle",
            side_effect=validate_with_test_roots,
        ):
            with ThreadPoolExecutor(max_workers=2) as executor:
                outcomes = list(executor.map(reserve_one, range(2)))
        self.assertEqual(sum(isinstance(value, dict) for value in outcomes), 1)
        self.assertEqual(sum(isinstance(value, ValueError) for value in outcomes), 1)

    def test_registered_science_rejects_launch_contract_drift(self) -> None:
        self._write_science_config(authorize=True)
        config = json.loads(self.config.read_text(encoding="utf-8"))
        config["launch_contract"]["actions"][0]["arguments"].append(
            {"literal": "time/nlim=2"}
        )
        self.config.write_text(json.dumps(config), encoding="utf-8")
        manifest_path = self._create_manifest()
        with self.assertRaisesRegex(ValueError, "launch contract"):
            self._reserve(manifest_path)

    def test_closed_admission_smoke_rejects_reuse(self) -> None:
        self._write_policy(admission_smoke_overrides={"status": "closed_after_pass"})
        self._promote_policy()
        manifest_path = self._create_manifest()
        with self.assertRaisesRegex(ValueError, "exemption is not authorized"):
            self._reserve(manifest_path)

    def test_manifest_schema_tracks_registered_science_authorization_id(self) -> None:
        schema = json.loads(
            (self.control_plane_dir / "control_plane.schema.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIn("registered_science_authorization_id", schema["properties"])
        branch = schema["allOf"][0]
        self.assertIn(
            "registered_science_authorization_id", branch["then"]["required"]
        )
        self.assertIn("pre_manifest_attestation_path", branch["then"]["required"])
        self.assertIn("pre_manifest_attestation_sha256", branch["then"]["required"])
        self.assertIn("prior_case_closures", branch["then"]["required"])
        prohibited = [
            value["required"][0]
            for value in branch["else"]["not"]["anyOf"]
        ]
        self.assertIn("registered_science_authorization_id", prohibited)
        self.assertIn("pre_manifest_attestation_path", prohibited)
        self.assertIn("pre_manifest_attestation_sha256", prohibited)
        self.assertIn("prior_case_closures", prohibited)

    def test_manifest_schema_closes_root_timeout_and_snapshot_records(self) -> None:
        schema = json.loads(
            (self.control_plane_dir / "control_plane.schema.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertFalse(schema["additionalProperties"])
        self.assertFalse(schema["properties"]["timeout_margin"]["additionalProperties"])
        self.assertFalse(
            schema["$defs"]["snapshot_file_record"]["additionalProperties"]
        )
        for key in [
            "control_plane_inventory",
            "campaign",
            "test_id",
            "git_commit",
            "evidence_class",
            "physical_mode",
            "artifact_dir",
            "queue_snapshot_sha256",
        ]:
            self.assertIn(key, schema["properties"])

    def test_reservation_rejects_extra_pre_submit_manifest_root_field(self) -> None:
        manifest_path = self._create_manifest()
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["unreviewed"] = True
        manifest_path.chmod(0o644)
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        manifest_path.chmod(0o444)
        with self.assertRaisesRegex(ValueError, "root schema"):
            self._reserve(manifest_path)

    def test_registered_science_rejects_policy_digest_mismatch(self) -> None:
        candidate = self._write_science_config(authorize=True)
        self._write_policy(
            science_submission_freeze=self._authorized_science_freeze(
                candidate, manifest_sha256="0" * 64
            ),
            admission_smoke_overrides={"status": "closed_after_pass"},
        )
        with self.assertRaisesRegex(ValueError, "another clean freeze"):
            self._promote_policy()

    def test_registered_science_rejects_missing_candidate_manifest(self) -> None:
        candidate = self._write_science_config(authorize=True)
        manifest_path = self._create_manifest()
        candidate.parent.chmod(0o755)
        candidate.unlink()
        candidate.parent.chmod(0o555)
        with self.assertRaises(FileNotFoundError):
            self._reserve(manifest_path)

    def test_registered_science_rejects_duplicate_candidate_keys(self) -> None:
        candidate = self._write_science_config(authorize=True)
        candidate.parent.chmod(0o755)
        candidate.chmod(0o644)
        candidate.write_text(
            '{"schema_version": 1, "schema_version": 1}\n',
            encoding="utf-8",
        )
        candidate.chmod(0o444)
        candidate.parent.chmod(0o555)
        self._update_registered_science_candidate_sha(candidate)
        self._write_policy(
            science_submission_freeze=self._authorized_science_freeze(candidate),
            admission_smoke_overrides={"status": "closed_after_pass"},
        )
        self._install_deliberately_malformed_active_policy_fixture()
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_registered_science_rejects_bound_executable_drift(self) -> None:
        self._write_science_config(authorize=True)
        manifest_path = self._create_manifest()
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        executable = Path(str(record_for_role(manifest, "executable")["source_path"]))
        executable.chmod(0o755)
        executable.write_text("mutated after freeze\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_registered_science_rejects_manifest_path_swap_after_stable_capture(
        self,
    ) -> None:
        candidate = self._write_science_config(authorize=True)
        manifest_path = self._create_manifest()
        from validate_and_reserve_frontier_job import _read_regular_file_at

        swapped = False

        def read_and_swap(directory_fd: int, name: str, *, label: str) -> bytes:
            nonlocal swapped
            data = _read_regular_file_at(directory_fd, name, label=label)
            if label == "Clean-candidate manifest" and not swapped:
                swapped = True
                candidate.parent.chmod(0o755)
                candidate.unlink()
                candidate.write_text('{"forged": true}\n', encoding="utf-8")
                candidate.chmod(0o444)
            return data

        with patch(
            "validate_and_reserve_frontier_job._read_regular_file_at",
            side_effect=read_and_swap,
        ):
            with self.assertRaises(ValueError):
                self._reserve(manifest_path)
        self.assertTrue(swapped)

    def test_registered_science_rejects_candidate_layout_widening(self) -> None:
        candidate = self._write_science_config(authorize=True)
        value = json.loads(candidate.read_text(encoding="utf-8"))
        alternate = self.pic_root / "alternate-athena"
        alternate.write_text("alternate\n", encoding="utf-8")
        value["build"]["executable_path"] = str(alternate)
        candidate.parent.chmod(0o755)
        candidate.chmod(0o644)
        candidate.write_text(json.dumps(value), encoding="utf-8")
        candidate.chmod(0o444)
        self._update_registered_science_candidate_sha(candidate)
        self._write_policy(
            science_submission_freeze=self._authorized_science_freeze(candidate),
            admission_smoke_overrides={"status": "closed_after_pass"},
        )
        self._install_deliberately_malformed_active_policy_fixture()
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_registered_science_rejects_forged_build_profile_source_root(self) -> None:
        candidate = self._write_science_config(authorize=True)
        candidate_value = json.loads(candidate.read_text(encoding="utf-8"))
        profile_path = Path(str(candidate_value["build"]["profile_path"]))
        profile = json.loads(profile_path.read_text(encoding="utf-8"))
        profile["authorized_source_root"] = str(self.root / "forged-source")
        self._rewrite_clean_candidate_profile(candidate, profile)
        manifest_path = self._create_manifest()
        with self.assertRaisesRegex(ValueError, "authorized source root"):
            self._reserve(manifest_path)

    def test_registered_science_rejects_forged_build_provenance_path(self) -> None:
        candidate = self._write_science_config(authorize=True)
        candidate_value = json.loads(candidate.read_text(encoding="utf-8"))
        profile_path = Path(str(candidate_value["build"]["profile_path"]))
        profile = json.loads(profile_path.read_text(encoding="utf-8"))
        profile["provenance_inputs"]["toolchain"]["path"] = str(
            self.pic_root / "bin" / "forged-toolchain.txt"
        )
        self._rewrite_clean_candidate_profile(candidate, profile)
        manifest_path = self._create_manifest()
        with self.assertRaisesRegex(ValueError, "documented Orion layout"):
            self._reserve(manifest_path)

    def test_registered_science_rejects_extra_candidate_file(self) -> None:
        candidate = self._write_science_config(authorize=True)
        candidate.parent.chmod(0o755)
        extra = candidate.parent / "unexpected.txt"
        extra.write_text("not part of the frozen layout\n", encoding="utf-8")
        extra.chmod(0o444)
        candidate.parent.chmod(0o555)
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_registered_science_rejects_self_consistent_snapshot_metadata_forgery(
        self,
    ) -> None:
        self._write_science_config(authorize=True)
        manifest_path = self._create_manifest()
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        executable = record_for_role(manifest, "executable")
        snapshot = Path(str(executable["path"]))
        snapshot.chmod(0o644)
        snapshot.write_text("forged executable snapshot\n", encoding="utf-8")
        executable["sha256"] = sha256(snapshot)
        executable["source_sha256"] = sha256(snapshot)
        manifest_path.chmod(0o644)
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        manifest_path.chmod(0o444)
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_admission_smoke_exemption_rejects_campaign_relabel(self) -> None:
        self._write_config(campaign="unit")
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_admission_smoke_exemption_rejects_unapproved_deck(self) -> None:
        self._write("input.athinput", "<job>\nbasename = relabeled-science\n")
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_admission_smoke_exemption_rejects_executable_drift(self) -> None:
        self._write("athena", "unreviewed smoke executable\n")
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_clean_candidate_creator_records_clean_git_attestation(self) -> None:
        source_root = self._clean_source("clean-source")
        executable, profile = self._build_profile(
            source_root, self.pic_root / "build", "test-profile"
        )
        manifest_path = create_freeze(
            source_root=source_root,
            executable=executable,
            build_profile=profile,
            build_profile_id="test-profile",
            prepared_artifact_inventory=self._prepared_artifact_inventory(),
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
        )
        candidate = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(candidate["source"]["worktree_status"], "clean")
        self.assertEqual(candidate["source"]["submodule_status"], "absent")
        self.assertEqual(candidate["build"]["executable_sha256"], sha256(executable))
        self.assertEqual(
            candidate["prepared_artifacts"],
            {
                "inventory_path": self._prepared_artifact_inventory(),
                "inventory_sha256": sha256(
                    source_root / self._prepared_artifact_inventory()
                ),
                "paper_decks": [
                    {
                        "path": (
                            "inputs/publication/"
                            "pic_parallel_shock_section54_paper_vl2_tsc.athinput"
                        ),
                        "sha256": sha256(
                            source_root
                            / "inputs/publication/"
                            "pic_parallel_shock_section54_paper_vl2_tsc.athinput"
                        ),
                    },
                    {
                        "path": "inputs/tests/pic_paper.athinput",
                        "sha256": sha256(source_root / "inputs/tests/pic_paper.athinput"),
                    }
                ],
                "analyzers": [
                    {
                        "path": "tst/publication/analyze_paper.py",
                        "sha256": sha256(source_root / "tst/publication/analyze_paper.py"),
                    }
                ],
            },
        )
        self.assertFalse(
            bool(Path(str(candidate["source"]["commit_path"])).stat().st_mode & 0o222)
        )
        self.assertFalse(bool(manifest_path.stat().st_mode & 0o222))
        self.assertFalse(
            bool(Path(str(candidate["build"]["executable_path"])).stat().st_mode & 0o222)
        )
        with Path(str(candidate["source"]["archive_path"])).open("rb") as archive:
            self.assertEqual(
                subprocess.check_output(
                    ["git", "get-tar-commit-id"], stdin=archive, text=True
                ).strip(),
                candidate["source"]["git_commit"],
            )
        self.assertFalse(list(manifest_path.parent.parent.glob(".tmp-*")))

    def test_clean_candidate_creator_rejects_omitted_required_publication_deck(
        self,
    ) -> None:
        source_root = self._clean_source("omitted-publication-deck-source")
        inventory = source_root / self._prepared_artifact_inventory()
        value = json.loads(inventory.read_text(encoding="utf-8"))
        value["paper_decks"] = [
            record
            for record in value["paper_decks"]
            if record["path"] != (
                "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
            )
        ]
        inventory.write_text(json.dumps(value), encoding="utf-8")
        subprocess.run(["git", "-C", str(source_root), "add", "."], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "-c",
                "user.name=PIC Test",
                "-c",
                "user.email=pic-test@example.invalid",
                "commit",
                "-m",
                "omit required publication deck",
            ],
            check=True,
            capture_output=True,
        )
        executable, profile = self._build_profile(
            source_root, self.pic_root / "omitted-publication-deck-build", "test-profile"
        )
        with self.assertRaisesRegex(ValueError, "must exactly cover archived"):
            create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="test-profile",
                prepared_artifact_inventory=self._prepared_artifact_inventory(),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
            )

    def test_clean_candidate_creator_rejects_missing_prepared_archive_member(self) -> None:
        source_root = self._clean_source("missing-prepared-member-source")
        inventory = source_root / self._prepared_artifact_inventory()
        value = json.loads(inventory.read_text(encoding="utf-8"))
        value["paper_decks"][0]["path"] = "inputs/tests/pic_missing.athinput"
        inventory.write_text(json.dumps(value), encoding="utf-8")
        subprocess.run(["git", "-C", str(source_root), "add", "."], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "-c",
                "user.name=PIC Test",
                "-c",
                "user.email=pic-test@example.invalid",
                "commit",
                "-m",
                "reference missing prepared deck",
            ],
            check=True,
            capture_output=True,
        )
        executable, profile = self._build_profile(
            source_root, self.pic_root / "missing-prepared-member-build", "test-profile"
        )
        with self.assertRaisesRegex(ValueError, "not an exact regular member"):
            create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="test-profile",
                prepared_artifact_inventory=self._prepared_artifact_inventory(),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
            )

    def test_prepared_artifact_inventory_rejects_noninteger_schema_version(self) -> None:
        source_root = self._clean_source("prepared-schema-source")
        inventory_path = source_root / self._prepared_artifact_inventory()
        original = json.loads(inventory_path.read_text(encoding="utf-8"))
        for value in (True, 1.0, "1"):
            with self.subTest(value=value):
                inventory = {**original, "schema_version": value}
                inventory_path.write_text(json.dumps(inventory), encoding="utf-8")
                subprocess.run(["git", "-C", str(source_root), "add", "."], check=True)
                subprocess.run(
                    [
                        "git",
                        "-C",
                        str(source_root),
                        "-c",
                        "user.name=PIC Test",
                        "-c",
                        "user.email=pic-test@example.invalid",
                        "commit",
                        "-m",
                        f"prepared schema {value!r}",
                    ],
                    check=True,
                    capture_output=True,
                )
                source_archive = subprocess.check_output(
                    ["git", "-C", str(source_root), "archive", "HEAD"]
                )
                with self.assertRaisesRegex(ValueError, "inventory schema"):
                    prepared_artifact_manifest_from_source_archive(
                        source_archive,
                        inventory_path=self._prepared_artifact_inventory(),
                    )

    def test_clean_candidate_creator_rejects_omitted_eligible_prepared_analyzer(self) -> None:
        source_root = self._clean_source("omitted-prepared-analyzer-source")
        omitted = source_root / "tst/publication/analyze_omitted.py"
        omitted.write_text("print('omitted analysis')\n", encoding="utf-8")
        subprocess.run(["git", "-C", str(source_root), "add", "."], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "-c",
                "user.name=PIC Test",
                "-c",
                "user.email=pic-test@example.invalid",
                "commit",
                "-m",
                "add omitted eligible analyzer",
            ],
            check=True,
            capture_output=True,
        )
        executable, profile = self._build_profile(
            source_root, self.pic_root / "omitted-prepared-analyzer-build", "test-profile"
        )
        with self.assertRaisesRegex(ValueError, "must exactly cover archived"):
            create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="test-profile",
                prepared_artifact_inventory=self._prepared_artifact_inventory(),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
            )

    def test_clean_candidate_creator_rejects_alternate_prepared_inventory_path(self) -> None:
        source_root = self._clean_source("alternate-prepared-inventory-source")
        canonical = source_root / self._prepared_artifact_inventory()
        alternate = source_root / "alternate/prepared_artifacts.json"
        alternate.parent.mkdir()
        alternate.write_bytes(canonical.read_bytes())
        subprocess.run(["git", "-C", str(source_root), "add", "."], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "-c",
                "user.name=PIC Test",
                "-c",
                "user.email=pic-test@example.invalid",
                "commit",
                "-m",
                "add alternate prepared inventory",
            ],
            check=True,
            capture_output=True,
        )
        executable, profile = self._build_profile(
            source_root, self.pic_root / "alternate-prepared-inventory-build", "test-profile"
        )
        with self.assertRaisesRegex(ValueError, "must use the canonical source-relative path"):
            create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="test-profile",
                prepared_artifact_inventory="alternate/prepared_artifacts.json",
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
            )

    def test_clean_candidate_creator_stages_captured_orion_input_bytes(self) -> None:
        source_root = self._clean_source("captured-input-source")
        executable, profile = self._build_profile(
            source_root, self.pic_root / "captured-input-build", "test-profile"
        )
        receipt = profile.with_name("profile_receipt.json")
        captured = {
            executable: executable.read_bytes(),
            profile: profile.read_bytes(),
            receipt: receipt.read_bytes(),
        }
        real_read = read_stable_regular_file_below
        mutated = False

        def mutate_after_capture(path: Path, root: Path, **kwargs: object) -> bytes:
            nonlocal mutated
            data = real_read(path, root, **kwargs)
            if not mutated and (Path(path) == receipt or Path(path) not in captured):
                mutated = True
                for artifact, payload, mode in [
                    (executable, b"forged executable\n", 0o555),
                    (profile, b'{"forged": true}\n', 0o444),
                    (receipt, b'{"forged": true}\n', 0o444),
                ]:
                    artifact.chmod(0o644)
                    artifact.write_bytes(payload)
                    artifact.chmod(mode)
            return data

        with patch(
            "create_clean_candidate_freeze.read_stable_regular_file_below",
            side_effect=mutate_after_capture,
        ):
            manifest_path = create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="test-profile",
                prepared_artifact_inventory=self._prepared_artifact_inventory(),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
            )
        self.assertTrue(mutated)
        candidate = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(
            Path(str(candidate["build"]["executable_path"])).read_bytes(),
            captured[executable],
        )
        self.assertEqual(
            Path(str(candidate["build"]["profile_path"])).read_bytes(),
            captured[profile],
        )
        self.assertEqual(
            Path(str(candidate["build"]["profile_receipt_path"])).read_bytes(),
            captured[receipt],
        )

    def test_orion_build_profile_writer_records_immutable_attestation(self) -> None:
        source_root = self._clean_source("profile-writer-source")
        arguments = self._profile_writer_arguments(source_root)
        executable = Path(str(arguments["executable"]))
        profile = Path(str(arguments["output"]))
        toolchain = Path(str(arguments["toolchain_file"]))
        toolchain.write_text("Frontier test toolchain\n", encoding="utf-8")
        profile = write_profile(**arguments)
        value = json.loads(profile.read_text(encoding="utf-8"))
        self.assertEqual(value["profile_id"], "test-profile")
        self.assertEqual(value["toolchain"], "Frontier test toolchain")
        self.assertEqual(
            value["build_invocations_sha256"],
            sha256(Path(str(arguments["build_invocations_file"]))),
        )
        self.assertEqual(value["executable_sha256"], sha256(executable))
        self.assertFalse(bool(profile.stat().st_mode & 0o222))
        self.assertFalse(list(profile.parent.glob(".build-profile-*")))
        with self.assertRaises(FileExistsError):
            write_profile(**arguments)

    def test_orion_build_profile_main_clone_is_bounded_and_submodules_preserve_status(
        self,
    ) -> None:
        child_source = self._clean_source("exact-clone-child-source")
        self._commit_source_change(
            child_source,
            relative="child-history.txt",
            content="second child revision\n",
            message="advance child source",
        )
        self._tag_source_head(child_source, "status-child")
        nested_source = self._clean_source("exact-clone-nested-source")
        self._add_existing_submodule(nested_source, child_source, "child")
        self._commit_source_change(
            nested_source,
            relative="nested-history.txt",
            content="second nested revision\n",
            message="advance nested source",
        )
        self._tag_source_head(nested_source, "status-nested")
        source_root = self._clean_source("exact-clone-source")
        self._add_existing_submodule(source_root, nested_source, "nested")
        self._commit_source_change(
            source_root,
            relative="top-history.txt",
            content="second top-level revision\n",
            message="advance top-level source",
        )
        subprocess.run(
            [
                "git",
                "-c",
                "protocol.file.allow=always",
                "-C",
                str(source_root),
                "submodule",
                "update",
                "--init",
                "--recursive",
            ],
            check=True,
            capture_output=True,
        )
        commit, tree, submodules = _source_identity(source_root)
        destination = self.root / "exact-clone-destination"
        clone_commands: list[list[str]] = []
        checkout_commands: list[list[str]] = []
        real_run = subprocess.run

        def track_clone_commands(
            command: list[str], *arguments: object, **keywords: object
        ) -> subprocess.CompletedProcess[bytes]:
            if "clone" in command:
                clone_commands.append(command)
            if "checkout" in command:
                checkout_commands.append(command)
            return real_run(command, *arguments, **keywords)

        with patch(
            "write_orion_build_profile.subprocess.run", side_effect=track_clone_commands
        ):
            _clone_fresh_source(
                source_root,
                destination,
                commit=commit,
                submodules=submodules,
            )

        ordered_submodules = sorted(
            submodules, key=lambda record: len(Path(record["path"]).parts)
        )
        self.assertEqual(len(clone_commands), 1 + len(ordered_submodules))
        main_clone = clone_commands[0]
        main_clone_index = main_clone.index("clone")
        self.assertEqual(
            main_clone[main_clone_index + 1 :],
            [
                "--no-local",
                "--no-tags",
                "--depth=1",
                f"--revision={commit}",
                str(source_root),
                str(destination),
            ],
        )
        for command, record in zip(
            clone_commands[1:], ordered_submodules, strict=True
        ):
            clone_index = command.index("clone")
            self.assertEqual(
                command[clone_index + 1 :],
                [
                    "--no-local",
                    "--reject-shallow",
                    "--no-checkout",
                    str(source_root / record["path"]),
                    str(destination / record["path"]),
                ],
            )
        self.assertEqual(len(checkout_commands), len(ordered_submodules))
        for command, record in zip(checkout_commands, ordered_submodules, strict=True):
            checkout_index = command.index("checkout")
            self.assertEqual(
                command[checkout_index + 1 :], ["--detach", record["git_commit"]]
            )

        self.assertEqual(_source_identity(destination), (commit, tree, submodules))
        authorized_submodule_status = subprocess.check_output(
            ["git", "-C", str(source_root), "submodule", "status", "--recursive"]
        )
        self.assertIn(b"(status-nested)", authorized_submodule_status)
        self.assertIn(b"(status-child)", authorized_submodule_status)
        self.assertEqual(
            subprocess.check_output(
                ["git", "-C", str(destination), "submodule", "status", "--recursive"]
            ),
            authorized_submodule_status,
        )
        repositories = [(destination, source_root, commit)] + [
            (
                destination / record["path"],
                source_root / record["path"],
                record["git_commit"],
            )
            for record in ordered_submodules
        ]
        for repository, original, revision in repositories:
            with self.subTest(repository=repository):
                self.assertEqual(
                    subprocess.check_output(
                        ["git", "-C", str(repository), "rev-parse", "HEAD"], text=True
                    ).strip(),
                    revision,
                )
                self.assertEqual(
                    subprocess.run(
                        ["git", "-C", str(repository), "symbolic-ref", "-q", "HEAD"],
                        capture_output=True,
                    ).returncode,
                    1,
                )
                self.assertEqual(
                    subprocess.check_output(
                        [
                            "git",
                            "-C",
                            str(repository),
                            "config",
                            "--get",
                            "remote.origin.url",
                        ],
                        text=True,
                    ).strip(),
                    str(original),
                )
        self.assertEqual(
            subprocess.check_output(
                ["git", "-C", str(destination), "rev-parse", "--is-shallow-repository"],
                text=True,
            ).strip(),
            "true",
        )
        self.assertEqual(
            subprocess.check_output(
                ["git", "-C", str(destination), "rev-list", "--count", "HEAD"],
                text=True,
            ).strip(),
            "1",
        )
        for repository, original, _ in repositories[1:]:
            with self.subTest(full_submodule=repository):
                self.assertEqual(
                    subprocess.check_output(
                        [
                            "git",
                            "-C",
                            str(repository),
                            "rev-parse",
                            "--is-shallow-repository",
                        ],
                        text=True,
                    ).strip(),
                    "false",
                )
                self.assertEqual(
                    subprocess.check_output(
                        ["git", "-C", str(repository), "rev-list", "--count", "HEAD"],
                        text=True,
                    ).strip(),
                    subprocess.check_output(
                        ["git", "-C", str(original), "rev-list", "--count", "HEAD"],
                        text=True,
                    ).strip(),
                )

        unavailable_source = source_root.with_name(f"{source_root.name}-unavailable")
        source_root.rename(unavailable_source)
        self.assertEqual(_source_identity(destination), (commit, tree, submodules))
        _require_exact_standalone_shallow_clone(destination, commit=commit)
        for repository, _, revision in repositories[1:]:
            _require_exact_standalone_full_clone(repository, commit=revision)
        for index, (repository, _, _) in enumerate(repositories):
            with self.subTest(standalone_repository=repository):
                subprocess.run(
                    [
                        "git",
                        "-C",
                        str(repository),
                        "archive",
                        "--format=tar",
                        f"--output={self.root / f'exact-clone-{index}.tar'}",
                        "HEAD",
                    ],
                    check=True,
                )
                subprocess.run(
                    ["git", "-C", str(repository), "cat-file", "commit", "HEAD"],
                    check=True,
                    capture_output=True,
                )

    def test_orion_build_profile_exact_clone_capabilities_precede_build_paths(
        self,
    ) -> None:
        _require_exact_clone_capabilities()
        source_root = self._clean_source("unsupported-exact-clone-capability-source")
        commit = subprocess.check_output(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"], text=True
        ).strip()
        profile_id = "unsupported-exact-clone-capability"
        build_root = self.pic_root / "build" / commit[:12] / profile_id
        artifact_root = self.pic_root / "bin" / commit[:12] / profile_id
        with patch(
            "write_orion_build_profile._require_exact_clone_capabilities",
            side_effect=ValueError("Trusted Git lacks required exact-clone capabilities"),
        ), patch(
            "write_orion_build_profile._execute_logged_command"
        ) as execute, self.assertRaisesRegex(
            ValueError, "lacks required exact-clone capabilities"
        ):
            build_orion_profile(
                source_root=source_root,
                expected_git_commit=commit,
                profile_id=profile_id,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_source_root=source_root,
            )
        execute.assert_not_called()
        self.assertFalse(build_root.exists())
        self.assertFalse(artifact_root.exists())

    def test_orion_build_profile_exact_clone_accepts_linked_authorized_worktree(
        self,
    ) -> None:
        source_root = self._clean_source("linked-authorized-source")
        linked_source = self.root / "linked-authorized-worktree"
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "worktree",
                "add",
                "--detach",
                str(linked_source),
                "HEAD",
            ],
            check=True,
            capture_output=True,
        )
        commit, tree, submodules = _source_identity(linked_source)
        destination = self.root / "linked-authorized-destination"
        _clone_fresh_source(
            linked_source,
            destination,
            commit=commit,
            submodules=submodules,
        )
        self.assertEqual(_source_identity(destination), (commit, tree, submodules))
        _require_exact_standalone_shallow_clone(destination, commit=commit)

    def test_orion_build_profile_rejects_recursive_status_drift_before_build(
        self,
    ) -> None:
        child_source = self._clean_source("status-drift-child-source")
        self._tag_source_head(child_source, "status-base")
        self._commit_source_change(
            child_source,
            relative="after-tag.txt",
            content="after tag\n",
            message="advance after tag",
        )
        source_root = self._clean_source("status-drift-source")
        self._add_existing_submodule(source_root, child_source, "child")
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root / "child"),
                "config",
                "--local",
                "core.abbrev",
                "12",
            ],
            check=True,
        )
        authorized_status = _submodule_status(source_root)
        self.assertIn(b"status-base-1-g", authorized_status)
        commit = subprocess.check_output(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"], text=True
        ).strip()
        with patch(
            "write_orion_build_profile._execute_logged_command"
        ) as execute, self.assertRaisesRegex(
            ValueError, "recursive submodule status differs"
        ):
            build_orion_profile(
                source_root=source_root,
                expected_git_commit=commit,
                profile_id="status-drift",
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_source_root=source_root,
            )
        execute.assert_not_called()

    def test_orion_build_profile_rejects_shallow_submodule_source(self) -> None:
        source_root = self._clean_source("full-submodule-source")
        self._commit_source_change(
            source_root,
            relative="history.txt",
            content="second revision\n",
            message="advance full source",
        )
        shallow_source = self.root / "shallow-submodule-source"
        subprocess.run(
            [
                "git",
                "clone",
                "--depth=1",
                f"file://{source_root}",
                str(shallow_source),
            ],
            check=True,
            capture_output=True,
        )
        commit = subprocess.check_output(
            ["git", "-C", str(shallow_source), "rev-parse", "HEAD"], text=True
        ).strip()
        self.assertEqual(
            subprocess.check_output(
                [
                    "git",
                    "-C",
                    str(shallow_source),
                    "rev-parse",
                    "--is-shallow-repository",
                ],
                text=True,
            ).strip(),
            "true",
        )
        with self.assertRaises(subprocess.CalledProcessError):
            _clone_exact_submodule(
                shallow_source,
                self.root / "rejected-shallow-submodule",
                commit=commit,
            )
        subprocess.run(
            ["git", "-C", str(shallow_source), "checkout", "--detach", commit],
            check=True,
            capture_output=True,
        )
        with self.assertRaisesRegex(ValueError, "submodule clone is shallow"):
            _require_exact_standalone_full_clone(shallow_source, commit=commit)

    def test_orion_build_profile_exact_revision_clone_postconditions_fail_closed(
        self,
    ) -> None:
        source_root = self._clean_source("adversarial-exact-clone-source")
        self._commit_source_change(
            source_root,
            relative="history.txt",
            content="second source revision\n",
            message="advance adversarial source",
        )
        commit = subprocess.check_output(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"], text=True
        ).strip()

        def fresh_clone(name: str) -> Path:
            destination = self.root / name
            _clone_fresh_source(source_root, destination, commit=commit, submodules=[])
            return destination

        wrong_revision = fresh_clone("wrong-revision-clone")
        with self.assertRaisesRegex(ValueError, "exact requested revision"):
            _require_exact_standalone_clone(wrong_revision, commit="0" * 40)

        attached = fresh_clone("attached-clone")
        subprocess.run(
            ["git", "-C", str(attached), "switch", "-c", "forged"],
            check=True,
            capture_output=True,
        )
        with self.assertRaisesRegex(ValueError, "not detached"):
            _require_exact_standalone_clone(attached, commit=commit)

        linked_source = fresh_clone("linked-source-clone")
        linked_worktree = self.root / "linked-worktree"
        subprocess.run(
            [
                "git",
                "-C",
                str(linked_source),
                "worktree",
                "add",
                "--detach",
                str(linked_worktree),
                "HEAD",
            ],
            check=True,
            capture_output=True,
        )
        with self.assertRaisesRegex(ValueError, "linked worktree"):
            _require_exact_standalone_clone(linked_worktree, commit=commit)

        unbounded = fresh_clone("unbounded-clone")
        unbounded_git_dir = Path(
            subprocess.check_output(
                ["git", "-C", str(unbounded), "rev-parse", "--absolute-git-dir"],
                text=True,
            ).strip()
        )
        (unbounded_git_dir / "shallow").unlink()
        with self.assertRaisesRegex(ValueError, "not shallow"):
            _require_exact_standalone_shallow_clone(unbounded, commit=commit)

        alternate = fresh_clone("alternate-clone")
        alternate_git_dir = Path(
            subprocess.check_output(
                ["git", "-C", str(alternate), "rev-parse", "--absolute-git-dir"],
                text=True,
            ).strip()
        )
        alternate_info = alternate_git_dir / "objects" / "info"
        alternate_info.mkdir(parents=True, exist_ok=True)
        (alternate_info / "alternates").write_text(
            f"{source_root / '.git' / 'objects'}\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError, "alternate object store"):
            _require_exact_standalone_clone(alternate, commit=commit)

        linked_objects = fresh_clone("linked-objects-clone")
        linked_objects_git_dir = Path(
            subprocess.check_output(
                ["git", "-C", str(linked_objects), "rev-parse", "--absolute-git-dir"],
                text=True,
            ).strip()
        )
        linked_pack = linked_objects_git_dir / "objects" / "pack"
        external_pack = self.root / "external-object-pack"
        linked_pack.rename(external_pack)
        linked_pack.symlink_to(external_pack, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "object store contains a symlink"):
            _require_exact_standalone_clone(linked_objects, commit=commit)

        hardlinked_objects = fresh_clone("hardlinked-objects-clone")
        hardlinked_objects_git_dir = Path(
            subprocess.check_output(
                [
                    "git",
                    "-C",
                    str(hardlinked_objects),
                    "rev-parse",
                    "--absolute-git-dir",
                ],
                text=True,
            ).strip()
        )
        object_file = next(
            path
            for path in (hardlinked_objects_git_dir / "objects").rglob("*")
            if path.is_file()
        )
        os.link(object_file, self.root / "external-hardlinked-object")
        with self.assertRaisesRegex(ValueError, "hard-linked object"):
            _require_exact_standalone_clone(hardlinked_objects, commit=commit)

        promisor_objects = fresh_clone("promisor-objects-clone")
        promisor_git_dir = Path(
            subprocess.check_output(
                ["git", "-C", str(promisor_objects), "rev-parse", "--absolute-git-dir"],
                text=True,
            ).strip()
        )
        (promisor_git_dir / "objects" / "pack" / "forged.promisor").touch()
        with self.assertRaisesRegex(ValueError, "promisor objects"):
            _require_exact_standalone_clone(promisor_objects, commit=commit)

        for key in [
            "remote.origin.promisor",
            "remote.origin.partialCloneFilter",
            "extensions.partialClone",
        ]:
            with self.subTest(configuration=key):
                promisor = fresh_clone(
                    f"promisor-clone-{key.lower().replace('.', '-')}"
                )
                subprocess.run(
                    ["git", "-C", str(promisor), "config", "--local", key, "true"],
                    check=True,
                )
                with self.assertRaisesRegex(ValueError, "promisor configuration"):
                    _require_exact_standalone_clone(promisor, commit=commit)

    def test_orion_build_profile_writer_rejects_dirty_source_tree(self) -> None:
        source_root = self._clean_source("dirty-profile-writer-source")
        (source_root / "untracked.txt").write_text("dirty\n", encoding="utf-8")
        arguments = self._profile_writer_arguments(source_root)
        profile = Path(str(arguments["output"]))
        with self.assertRaises(ValueError):
            write_profile(**arguments)
        self.assertFalse(profile.exists())

    def test_orion_build_profile_writer_rejects_invalid_metadata_and_aliases(self) -> None:
        source_root = self._clean_source("invalid-profile-writer-source")
        arguments = self._profile_writer_arguments(source_root)
        toolchain = Path(str(arguments["toolchain_file"]))
        with self.assertRaises(ValueError):
            write_profile(**{**arguments, "profile_id": " "})
        toolchain.write_bytes(b"\xff")
        with self.assertRaises(ValueError):
            write_profile(**arguments)
        for alias in (
            b" Frontier test toolchain\n",
            b"Frontier test toolchain \n",
            b"Frontier test toolchain",
            b"Frontier test toolchain\n\n",
        ):
            with self.subTest(alias=alias):
                toolchain.write_bytes(alias)
                with self.assertRaises(ValueError):
                    write_profile(**arguments)
        toolchain.write_text("Frontier test toolchain\n", encoding="utf-8")
        outside = self.root / "outside-profile-writer"
        outside.mkdir()
        alias = self.pic_root / "profile-writer-alias"
        alias.symlink_to(outside, target_is_directory=True)
        with self.assertRaises(ValueError):
            write_profile(**{**arguments, "output": alias / "build_profile.json"})
        self.assertEqual(list(outside.iterdir()), [])

    def test_orion_build_profile_writer_rejects_tracked_symlink_payload(self) -> None:
        source_root = self._clean_source("symlink-profile-writer-source")
        (source_root / "tracked-link").symlink_to("tracked.txt")
        subprocess.run(["git", "-C", str(source_root), "add", "tracked-link"], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "-c",
                "user.name=PIC Test",
                "-c",
                "user.email=pic-test@example.invalid",
                "commit",
                "-m",
                "add tracked symlink",
            ],
            check=True,
            capture_output=True,
        )
        commit = subprocess.check_output(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"], text=True
        ).strip()
        with patch("write_orion_build_profile._execute_logged_command") as execute:
            with self.assertRaises(ValueError):
                build_orion_profile(
                    source_root=source_root,
                    expected_git_commit=commit,
                    profile_id="test-profile",
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_source_root=source_root,
                )
        execute.assert_not_called()

    def test_orion_build_profile_api_has_no_caller_injection_hooks(self) -> None:
        parameters = inspect.signature(build_orion_profile).parameters
        self.assertNotIn("invocations_factory", parameters)
        self.assertNotIn("runner", parameters)
        self.assertNotIn("toolchain_description", parameters)

    def test_orion_build_profile_subprocess_environment_excludes_caller_poison(self) -> None:
        import write_orion_build_profile

        with patch.dict(
            os.environ,
            {
                "CMAKE_TOOLCHAIN_FILE": "/tmp/forged-cmake-toolchain",
                "GIT_CONFIG_GLOBAL": "/tmp/forged-git-config",
                "LD_PRELOAD": "/tmp/forged-loader.so",
                "PYTHONPATH": "/tmp/forged-python",
            },
        ):
            environment = write_orion_build_profile._production_build_environment()
        self.assertNotIn("CMAKE_TOOLCHAIN_FILE", environment)
        self.assertNotIn("GIT_CONFIG_GLOBAL", environment)
        self.assertNotIn("LD_PRELOAD", environment)
        self.assertNotIn("PYTHONPATH", environment)
        with patch.dict(
            os.environ,
            {
                "CRAY_FORGED": "apparently-safe-wrapper-option",
                "PE_ENV": "FORGED",
            },
        ):
            environment = write_orion_build_profile._production_build_environment()
        self.assertNotIn("CRAY_FORGED", environment)
        self.assertEqual(environment["PE_ENV"], "AMD")

    def test_orion_build_profile_rejects_symlinked_generated_outputs(self) -> None:
        for linked_output in ["athena", "CMakeCache.txt"]:
            with self.subTest(linked_output=linked_output):
                profile_id = f"linked-{linked_output.lower().replace('.', '-')}"
                source_root = self._clean_source(f"linked-{linked_output}-source")
                commit = subprocess.check_output(
                    ["git", "-C", str(source_root), "rev-parse", "HEAD"], text=True
                ).strip()
                outside = self.root / f"outside-{linked_output}"
                outside.write_text("outside payload\n", encoding="utf-8")

                def execute(
                    command: list[str], *, stream: object, environment: dict[str, str]
                ) -> None:
                    del environment
                    stream.write(b"fake cmake invocation\n")
                    if "--build" in command:
                        cmake_dir = Path(command[command.index("--build") + 1])
                        (cmake_dir / "src").mkdir(parents=True)
                        executable = cmake_dir / "src" / "athena"
                        if linked_output == "athena":
                            executable.symlink_to(outside)
                        else:
                            executable.write_text("built executable\n", encoding="utf-8")
                    else:
                        cmake_dir = Path(command[command.index("-B") + 1])
                        cmake_dir.mkdir(parents=True)
                        cache = cmake_dir / "CMakeCache.txt"
                        if linked_output == "CMakeCache.txt":
                            cache.symlink_to(outside)
                        else:
                            cache.write_text("fixture cache\n", encoding="utf-8")

                with patch(
                    "write_orion_build_profile._execute_logged_command",
                    side_effect=execute,
                ):
                    with self.assertRaises(OSError):
                        build_orion_profile(
                            source_root=source_root,
                            expected_git_commit=commit,
                            profile_id=profile_id,
                            control_plane_dir=self.control_plane_dir,
                            authorized_pic_root=self.pic_root,
                            authorized_source_root=source_root,
                        )

    def test_orion_build_profile_writer_rejects_executable_drift(self) -> None:
        source_root = self._clean_source("drifted-profile-writer-source")
        arguments = self._profile_writer_arguments(source_root)
        executable = Path(str(arguments["executable"]))
        from write_orion_build_profile import read_stable_regular_file_below

        executable_reads = 0

        def read_and_drift(path: Path, root: Path) -> bytes:
            nonlocal executable_reads
            if path == executable:
                executable_reads += 1
                if executable_reads == 2:
                    executable.write_text("drifted executable\n", encoding="utf-8")
            return read_stable_regular_file_below(path, root)

        profile = Path(str(arguments["output"]))
        with patch(
            "write_orion_build_profile.read_stable_regular_file_below",
            side_effect=read_and_drift,
        ):
            with self.assertRaises(ValueError):
                write_profile(**arguments)
        self.assertEqual(executable_reads, 2)
        self.assertFalse(profile.exists())

    def test_orion_build_profile_writer_rejects_dirty_recursive_submodule(self) -> None:
        source_root = self._clean_source("dirty-submodule-profile-writer-source")
        nested = self._add_submodule(source_root, "nested")
        (nested / "untracked.txt").write_text("dirty\n", encoding="utf-8")
        arguments = self._profile_writer_arguments(source_root)
        profile = Path(str(arguments["output"]))
        with self.assertRaises(ValueError):
            write_profile(**arguments)
        self.assertFalse(profile.exists())

    def test_orion_build_profile_writer_atomic_failure_cleans_staging(self) -> None:
        source_root = self._clean_source("failed-profile-writer-source")
        arguments = self._profile_writer_arguments(source_root)
        profile = Path(str(arguments["output"]))
        with patch("control_plane_common.os.link", side_effect=OSError("link failed")):
            with self.assertRaises(OSError):
                write_profile(**arguments)
        self.assertFalse(profile.exists())
        self.assertFalse(list(profile.parent.glob(".build-profile-*")))
        self.assertFalse(list(profile.parent.glob(".build_profile.json.tmp-*")))

    def test_orion_build_profile_writer_parent_fsync_failure_rolls_back(self) -> None:
        source_root = self._clean_source("fsync-failed-profile-writer-source")
        arguments = self._profile_writer_arguments(source_root)
        profile = Path(str(arguments["output"]))
        real_fsync = os.fsync
        directory_fsyncs = 0

        def fail_publish_parent_once(descriptor: int) -> None:
            nonlocal directory_fsyncs
            if stat.S_ISDIR(os.fstat(descriptor).st_mode):
                directory_fsyncs += 1
                if directory_fsyncs == 1:
                    raise OSError("publication parent fsync failed")
            real_fsync(descriptor)

        with patch("control_plane_common.os.fsync", side_effect=fail_publish_parent_once):
            with self.assertRaises(OSError):
                write_profile(**arguments)
        self.assertFalse(profile.exists())
        self.assertFalse(list(profile.parent.glob(".build-profile-*")))
        self.assertFalse(list(profile.parent.glob(".build_profile.json.tmp-*")))

    def test_orion_build_profile_writer_parent_swap_rejects_lexical_bytes(self) -> None:
        source_root = self._clean_source("swapped-profile-writer-source")
        arguments = self._profile_writer_arguments(source_root)
        profile = Path(str(arguments["output"]))
        parent = profile.parent
        moved_parent = parent.with_name(parent.name + "-original")
        outside = self.root / "outside-profile-parent"
        outside.mkdir()
        outside_profile = outside / profile.name
        outside_profile.write_text('{"generation": "outside"}\n', encoding="utf-8")
        parent_inode = parent.stat().st_ino
        real_fsync = os.fsync
        swapped = False

        def swap_parent_then_sync(descriptor: int) -> None:
            nonlocal swapped
            descriptor_stat = os.fstat(descriptor)
            if (
                stat.S_ISDIR(descriptor_stat.st_mode)
                and descriptor_stat.st_ino == parent_inode
                and not swapped
            ):
                swapped = True
                parent.rename(moved_parent)
                parent.symlink_to(outside, target_is_directory=True)
            real_fsync(descriptor)

        with patch("control_plane_common.os.fsync", side_effect=swap_parent_then_sync):
            with self.assertRaises((OSError, ValueError)):
                write_profile(**arguments)
        self.assertTrue(swapped)
        self.assertEqual(
            outside_profile.read_text(encoding="utf-8"),
            '{"generation": "outside"}\n',
        )
        self.assertFalse((moved_parent / profile.name).exists())

    def test_clean_candidate_creator_rejects_dirty_source_tree(self) -> None:
        source_root = self.root / "dirty-source"
        source_root.mkdir()
        subprocess.run(["git", "init", str(source_root)], check=True, capture_output=True)
        (source_root / "untracked.txt").write_text("untracked\n", encoding="utf-8")
        build = self.pic_root / "build"
        build.mkdir()
        executable = build / "athena"
        executable.write_text("built executable\n", encoding="utf-8")
        profile = build / "build_profile.json"
        profile.write_text("{}\n", encoding="utf-8")
        (build / "profile_receipt.json").write_text("{}\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="test-profile",
                prepared_artifact_inventory=self._prepared_artifact_inventory(),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
            )

    def test_clean_candidate_creator_rejects_structured_profile_mismatch(self) -> None:
        source_root = self._clean_source("profile-mismatch-source")
        executable, profile = self._build_profile(
            source_root, self.pic_root / "profile-mismatch-build", "test-profile"
        )
        value = json.loads(profile.read_text(encoding="utf-8"))
        value["source_archive_sha256"] = "0" * 64
        profile.chmod(0o644)
        profile.write_text(json.dumps(value), encoding="utf-8")
        with self.assertRaises(ValueError):
            create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="test-profile",
                prepared_artifact_inventory=self._prepared_artifact_inventory(),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
            )

    def test_clean_candidate_creator_rejects_noninteger_profile_schema_versions(self) -> None:
        source_root = self._clean_source("profile-schema-source")
        executable, profile = self._build_profile(
            source_root, self.pic_root / "profile-schema-build", "test-profile"
        )
        receipt = profile.with_name("profile_receipt.json")
        originals = {
            profile: json.loads(profile.read_text(encoding="utf-8")),
            receipt: json.loads(receipt.read_text(encoding="utf-8")),
        }
        for path, original in originals.items():
            for schema_version in (True, 1.0, "1"):
                with self.subTest(path=path.name, schema_version=schema_version):
                    path.chmod(0o644)
                    path.write_text(
                        json.dumps({**original, "schema_version": schema_version}),
                        encoding="utf-8",
                    )
                    path.chmod(0o444)
                    with self.assertRaisesRegex(ValueError, "schema version is invalid"):
                        create_freeze(
                            source_root=source_root,
                            executable=executable,
                            build_profile=profile,
                            build_profile_id="test-profile",
                            prepared_artifact_inventory=self._prepared_artifact_inventory(),
                            control_plane_dir=self.control_plane_dir,
                            authorized_pic_root=self.pic_root,
                        )
                    path.chmod(0o644)
                    path.write_text(json.dumps(original), encoding="utf-8")
                    path.chmod(0o444)

    def test_clean_candidate_rename_failure_cleans_read_only_staging(self) -> None:
        source_root = self._clean_source("candidate-rename-failure-source")
        executable, profile = self._build_profile(
            source_root, self.pic_root / "candidate-rename-failure-build", "test-profile"
        )
        candidate_root = self.pic_root / "clean_candidates"
        with patch(
            "control_plane_common.os.replace",
            side_effect=OSError("rename failed"),
        ):
            with self.assertRaises(OSError):
                create_freeze(
                    source_root=source_root,
                    executable=executable,
                    build_profile=profile,
                    build_profile_id="test-profile",
                    prepared_artifact_inventory=self._prepared_artifact_inventory(),
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                )
        self.assertFalse(list(candidate_root.glob(".tmp-*")))

    def test_clean_candidate_creator_archives_clean_pinned_submodules(self) -> None:
        source_root = self._clean_source("submodule-source")
        self._add_submodule(source_root, "nested")
        executable, profile = self._build_profile(
            source_root, self.pic_root / "submodule-build", "test-profile"
        )
        manifest_path = create_freeze(
            source_root=source_root,
            executable=executable,
            build_profile=profile,
            build_profile_id="test-profile",
            prepared_artifact_inventory=self._prepared_artifact_inventory(),
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_source_root=source_root,
        )
        candidate = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(candidate["source"]["submodule_status"], "clean_pinned_archived")
        self.assertEqual(len(candidate["source"]["submodules"]), 1)
        archive = Path(str(candidate["source"]["submodules"][0]["archive_path"]))
        self.assertEqual(archive, manifest_path.parent / "submodules" / "0000.tar")
        self.assertFalse(bool(archive.stat().st_mode & 0o222))
        commit_object = Path(str(candidate["source"]["submodules"][0]["commit_path"]))
        self.assertEqual(commit_object, manifest_path.parent / "submodules" / "0000.commit")
        self.assertFalse(bool(commit_object.stat().st_mode & 0o222))

    def test_authorized_source_path_normalizes_trusted_source_root_alias(self) -> None:
        source_root = self._clean_source("aliased-submodule-source")
        self._add_submodule(source_root, "nested")
        aliased_source_root = self.root / "aliased-submodule-source-root"
        aliased_source_root.symlink_to(source_root, target_is_directory=True)
        trusted_source_root = _authorized_source_path(
            aliased_source_root, aliased_source_root
        )
        self.assertEqual(trusted_source_root, source_root.resolve())
        records = _validated_submodules(trusted_source_root)
        self.assertEqual([record["path"] for record in records], ["nested"])

    def test_validated_submodules_rejects_source_root_alias(self) -> None:
        source_root = self._clean_source("untrusted-aliased-submodule-source")
        self._add_submodule(source_root, "nested")
        aliased_source_root = self.root / "untrusted-aliased-submodule-source-root"
        aliased_source_root.symlink_to(source_root, target_is_directory=True)
        with self.assertRaises(ValueError):
            _validated_submodules(aliased_source_root)

    def test_validated_submodules_rejects_empty_source_root_alias(self) -> None:
        source_root = self._clean_source("untrusted-empty-aliased-source")
        aliased_source_root = self.root / "untrusted-empty-aliased-source-root"
        aliased_source_root.symlink_to(source_root, target_is_directory=True)
        with self.assertRaises(ValueError):
            _validated_submodules(aliased_source_root)

    def test_validated_submodules_rejects_submodule_path_symlink(self) -> None:
        source_root = self.root / "symlinked-submodule-source"
        source_root.mkdir()
        nested_real = self.root / "symlinked-submodule-target"
        nested_real.mkdir()
        (source_root / "nested").symlink_to(nested_real, target_is_directory=True)
        with patch(
            "create_clean_candidate_freeze._git",
            return_value=f" {'0' * 40} nested",
        ):
            with self.assertRaisesRegex(ValueError, "traverses a symlink"):
                _validated_submodules(source_root)

    def test_orion_build_profile_writer_accepts_trusted_source_root_alias(self) -> None:
        source_root = self._clean_source("aliased-profile-writer-source")
        self._add_submodule(source_root, "nested")
        aliased_source_root = self.root / "aliased-profile-writer-source-root"
        aliased_source_root.symlink_to(source_root, target_is_directory=True)
        _, profile = self._build_profile(
            aliased_source_root,
            self.pic_root / "aliased-profile-writer-build",
            "aliased-profile-writer",
        )
        value = json.loads(profile.read_text(encoding="utf-8"))
        self.assertEqual(value["authorized_source_root"], str(aliased_source_root))

    def test_clean_candidate_creator_rejects_dirty_submodule(self) -> None:
        source_root = self._clean_source("dirty-submodule-source")
        nested = self._add_submodule(source_root, "nested")
        executable, profile = self._build_profile(
            source_root, self.pic_root / "dirty-submodule-build", "test-profile"
        )
        (nested / "untracked.txt").write_text("dirty\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="test-profile",
                prepared_artifact_inventory=self._prepared_artifact_inventory(),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_source_root=source_root,
            )

    def test_clean_candidate_creator_rejects_tracked_symlink_payload(self) -> None:
        source_root = self._clean_source("symlink-payload-source")
        (source_root / "tracked-link").symlink_to("tracked.txt")
        subprocess.run(["git", "-C", str(source_root), "add", "tracked-link"], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "-c",
                "user.name=PIC Test",
                "-c",
                "user.email=pic-test@example.invalid",
                "commit",
                "-m",
                "add tracked symlink",
            ],
            check=True,
            capture_output=True,
        )
        with self.assertRaises(ValueError):
            self._build_profile(
                source_root, self.pic_root / "symlink-payload-build", "test-profile"
            )

    def test_git_tree_reconstruction_rejects_missing_gitlink_placeholder(self) -> None:
        source_root = self._clean_source("missing-gitlink-placeholder-source")
        archive = self.root / "missing-gitlink-placeholder.tar"
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "archive",
                "--format=tar",
                f"--output={archive}",
                "HEAD",
            ],
            check=True,
        )
        with self.assertRaises(ValueError):
            git_tree_sha1_from_archive(
                archive, gitlinks={"missing": "0" * 40}, reject_symlinks=True
            )

    def test_clean_candidate_creator_archives_recursive_pinned_submodules(self) -> None:
        source_root = self._clean_source("recursive-submodule-source")
        child = self._clean_source("recursive-child-source")
        self._tag_source_head(child, "status-child")
        nested = self._clean_source("recursive-nested-source")
        self._add_existing_submodule(nested, child, "child")
        self._tag_source_head(nested, "status-nested")
        self._add_existing_submodule(source_root, nested, "nested")
        subprocess.run(
            [
                "git",
                "-c",
                "protocol.file.allow=always",
                "-C",
                str(source_root),
                "submodule",
                "update",
                "--init",
                "--recursive",
            ],
            check=True,
            capture_output=True,
        )
        authorized_submodule_status = subprocess.check_output(
            ["git", "-C", str(source_root), "submodule", "status", "--recursive"]
        )
        self.assertIn(b"(status-nested)", authorized_submodule_status)
        self.assertIn(b"(status-child)", authorized_submodule_status)
        executable, profile = self._build_profile(
            source_root, self.pic_root / "recursive-submodule-build", "test-profile"
        )
        self.assertEqual(
            profile.with_name("submodule_status.txt").read_bytes(),
            authorized_submodule_status,
        )
        manifest_path = create_freeze(
            source_root=source_root,
            executable=executable,
            build_profile=profile,
            build_profile_id="test-profile",
            prepared_artifact_inventory=self._prepared_artifact_inventory(),
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
        )
        candidate = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(
            [record["path"] for record in candidate["source"]["submodules"]],
            ["nested", "nested/child"],
        )

    def test_clean_candidate_creator_rejects_uninitialized_submodule(self) -> None:
        source_root = self._clean_source("uninitialized-submodule-source")
        self._add_submodule(source_root, "nested")
        subprocess.run(
            ["git", "-C", str(source_root), "submodule", "deinit", "-f", "nested"],
            check=True,
            capture_output=True,
        )
        with self.assertRaises(ValueError):
            _validated_submodules(source_root)

    def test_clean_candidate_creator_rejects_submodule_commit_drift(self) -> None:
        source_root = self._clean_source("drifted-submodule-source")
        nested = self._add_submodule(source_root, "nested")
        (nested / "tracked.txt").write_text("drifted\n", encoding="utf-8")
        subprocess.run(["git", "-C", str(nested), "add", "tracked.txt"], check=True)
        subprocess.run(
            [
                "git",
                "-C",
                str(nested),
                "-c",
                "user.name=PIC Test",
                "-c",
                "user.email=pic-test@example.invalid",
                "commit",
                "-m",
                "drift nested source",
            ],
            check=True,
            capture_output=True,
        )
        with self.assertRaises(ValueError):
            _validated_submodules(source_root)

    def test_registered_science_rejects_frozen_submodule_archive_drift(self) -> None:
        candidate_path = self._write_science_config(authorize=True)
        candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
        archive = Path(str(candidate["source"]["submodules"][0]["archive_path"]))
        archive.chmod(0o644)
        archive.write_bytes(archive.read_bytes() + b"drift")
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_reserved_manifest_digest_rejects_self_consistent_replacement(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest_path.chmod(0o644)
        manifest_path.write_text(json.dumps(manifest, indent=4), encoding="utf-8")
        with self.assertRaises(ValueError):
            verify(
                manifest_path,
                manifest_sha256=str(reservation["manifest_sha256"]),
            )

    def test_reservation_rejects_cap_overrun(self) -> None:
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path, cap=0.01)

    def test_reservation_rejects_cap_expansion(self) -> None:
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path, cap=10000.01)

    def test_normal_fallback_requires_claimed_debug_job(self) -> None:
        self._write(
            "job.sh",
            "#!/bin/bash\n#SBATCH -A AST207\n#SBATCH -p batch\n#SBATCH -q normal\n"
            f"#SBATCH -o {self.pic_root}/logs/slurm/%x.%j.log\n"
            "#SBATCH -N 1\n#SBATCH -t 00:10:00\n",
        )
        self._write_config(
            selected_qos="normal",
            qos_selection_reason="debug_slot_occupied",
        )
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_reservation_rejects_changed_queue_snapshot(self) -> None:
        manifest_path = self._create_manifest()
        self._write("queue.txt", "123|batch|normal|RUNNING|other-job|\n")
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_reservation_rejection_after_manifest_publication_leaves_no_ledger_intent(
        self,
    ) -> None:
        manifest_path = self._create_manifest()
        before = {
            path: path.read_bytes()
            for path in [self.ledger, self.csv, self.receipts, self.mirror]
        }
        self._write("queue.txt", "123|batch|normal|RUNNING|other-job|\n")
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)
        self.assertTrue(manifest_path.is_file())
        for name in ["reservation_id.txt", "manifest_sha256.txt"]:
            self.assertFalse((manifest_path.parent / name).exists())
        self.assertFalse((self.pic_root / "ledger" / "pending_submission.json").exists())
        for path, expected in before.items():
            self.assertEqual(path.read_bytes(), expected)

    def test_reservation_rejects_attestation_seven_field_queue_snapshot_without_ledger_intent(
        self,
    ) -> None:
        self._write("queue.txt", "123|AST207|batch|normal|RUNNING|other-job|\n")
        manifest_path = self._create_manifest()
        before = {
            path: path.read_bytes()
            for path in [self.ledger, self.csv, self.receipts, self.mirror]
        }
        self._write("queue.txt", "123|batch|normal|RUNNING|other-job|\n")
        with self.assertRaisesRegex(
            ValueError, "Fresh queue output differs from frozen queue snapshot"
        ):
            self._reserve(manifest_path)
        self.assertTrue(manifest_path.is_file())
        for name in ["reservation_id.txt", "manifest_sha256.txt"]:
            self.assertFalse((manifest_path.parent / name).exists())
        self.assertFalse((self.pic_root / "ledger" / "pending_submission.json").exists())
        for path, expected in before.items():
            self.assertEqual(path.read_bytes(), expected)

    def test_reservation_rejects_untrusted_environment_profile(self) -> None:
        self._write("environment.sh", "export MPICH_GPU_SUPPORT_ENABLED=0\n")
        self._write_timeout()
        self._write_policy()
        self._promote_policy()
        self._write_config()
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_initializer_rejects_locked_storage_policy(self) -> None:
        self._write_policy(ledger_genesis=None)
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_initializer_rejects_unapproved_control_plane_version(self) -> None:
        self._write_policy(installed_control_plane_version="0" * 64)
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_initializer_rejects_storage_policy_schema_and_boolean_aliases(self) -> None:
        for mutate in [
            lambda policy: policy.update(schema_version=True),
            lambda policy: policy["frontier"].update(maximum_node_hours=True),
            lambda policy: policy["frontier"].update(serial_pic_submissions=1),
            lambda policy: policy["frontier_admission_smoke"].update(maximum_nodes=True),
            lambda policy: policy["frontier_admission_smoke"].update(
                registered_short_nonproduction=1
            ),
        ]:
            with self.subTest(mutate=mutate):
                self._write_policy()
                policy = json.loads(self.policy.read_text(encoding="utf-8"))
                mutate(policy)
                self.policy.write_text(json.dumps(policy), encoding="utf-8")
                with self.assertRaises(ValueError):
                    self._promote_policy()

    def test_policy_accepts_closed_production_storage_metadata_shape(self) -> None:
        self._write_policy(
            status="passed_user_authorized_orion_only_storage",
            historical_project_home_bulk_artifacts=(
                "chronology_only_superseded_by_orion_policy_copies_do_not_add_new_bulk_artifacts"
            ),
            ledger_genesis_authorization=(
                "user_removed_kronos_dependency_and_selected_orion_only_bulk_evidence_root"
            ),
            orion_simulation_root_preflight={
                "status": "passed",
                "path": str(self.pic_root),
                "method": "local_create_write_sync_remove_probe",
            },
            project_home_preflight={
                "status": "passed",
                "path": str(self.project_home_root),
                "method": "local_create_write_sync_remove_probe",
            },
        )
        self._promote_policy()

    def test_storage_preflight_reviewed_source_digest_tuple_matches_current_bytes(
        self,
    ) -> None:
        expected = {
            "entrypoint_sha256": "capture_storage_preflight_evidence.py",
            "runner_sha256": "run_control_plane.py",
            "schema_sha256": "storage_preflight.schema.json",
        }
        control_plane = Path(__file__).parent
        self.assertEqual(
            AUTHORIZED_STORAGE_PREFLIGHT_CAPTURE_SOURCE_BLOBS,
            {
                key: hashlib.sha256((control_plane / filename).read_bytes()).hexdigest()
                for key, filename in expected.items()
            },
        )
        self.assertNotIn(
            "common_sha256",
            AUTHORIZED_STORAGE_PREFLIGHT_CAPTURE_SOURCE_BLOBS,
        )

    def test_exact_reviewed_storage_preflight_predecessor_authorization_is_literal(
        self,
    ) -> None:
        self.assertEqual(
            control_plane_common.
            AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_POLICY_SHA256,
            "48e74f3151b51ba84b72e3214ef4a66c03441c745912b833395f98f6f60a0bd1",
        )
        self.assertEqual(
            control_plane_common.
            AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_PROMOTION_SHA256,
            "4ecb45bb399cee750ce69198286afc9fb903dde92cffc7507600b0a53f1c547f",
        )
        self.assertEqual(
            control_plane_common.
            AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_CONTROL_PLANE_VERSION,
            "b56d96b40f2c666b9fa5b421fac589d6a4d6500a716f20b479c354a3d239cb48",
        )
        self.assertEqual(
            control_plane_common.
            AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_PROBE_ID,
            "bc399b56-8fbb-4b67-b1dc-5df1dbff62b6",
        )
        self.assertEqual(
            control_plane_common.
            AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_EVIDENCE_SHA256,
            "d4a289ce9f4cd7c406dbea8d457864790f3112488e47ea24766cb192beaafab2",
        )
        self.assertEqual(
            AUTHORIZED_STORAGE_PREFLIGHT_PREDECESSOR_MIGRATION_SOURCE_AUTHENTICATION,
            {
                "common_sha256": (
                    "8873527c19b46236ed66f1c8e3b0318cfcb106e017c5289eac70b482db0ef963"
                ),
                "entrypoint_sha256": (
                    "b6dae64b28dbcc7ce82877ad25d53bc4f0637016c4bd274431c1a4ba947ec94b"
                ),
                "git_commit": "f6471610a116ce5433550625c9dc61752315040f",
                "runner_sha256": (
                    "6053f190ed5bea093537ca5e6aef110212d54861f6726a6fa7294a1716eca2d5"
                ),
                "schema_sha256": (
                    "348b80f6b56fa56da57a4d30932b41784c939f5a2b1bf49c82d3a6acd024ee3f"
                ),
                "tracked_clean_head_blobs": True,
            },
        )

    def test_policy_rejects_missing_storage_preflight_binding_and_root_fields(
        self,
    ) -> None:
        for overrides in [
            {"storage_preflight_evidence": None},
            {
                "orion_simulation_root_preflight": {
                    "path": str(self.pic_root),
                    "status": "passed",
                }
            },
            {
                "project_home_preflight": {
                    "method": "local_create_write_sync_remove_probe",
                    "status": "passed",
                }
            },
        ]:
            with self.subTest(overrides=overrides):
                self._write_policy(**overrides)
                with self.assertRaises(ValueError):
                    self._promote_policy()

    def test_policy_rejects_storage_preflight_completion_mismatch(self) -> None:
        self._write_policy(last_preflight_utc="2026-06-03T00:00:01Z")
        with self.assertRaisesRegex(ValueError, "completion differs"):
            self._promote_policy()

    def test_policy_rejects_divergent_storage_preflight_mirror(self) -> None:
        self._write_policy()
        policy = json.loads(self.policy.read_text(encoding="utf-8"))
        mirror = Path(
            str(
                policy["olcf_side_storage"]["storage_preflight_evidence"][
                    "project_home_path"
                ]
            )
        )
        mirror.chmod(0o644)
        mirror.write_bytes(mirror.read_bytes() + b"divergent\n")
        mirror.chmod(0o444)
        with self.assertRaisesRegex(ValueError, "mirrored bytes differ"):
            self._promote_policy()

    def test_policy_rejects_unreviewed_storage_preflight_capture_source(self) -> None:
        self._write_policy()
        self._rewrite_reviewed_storage_preflight_artifact(
            lambda artifact: artifact["source_authentication"].update(
                runner_sha256="0" * 64
            )
        )
        with self.assertRaisesRegex(ValueError, "source authentication"):
            self._promote_policy()

    def test_historical_storage_preflight_retirement_is_promotion_only_and_one_time(
        self,
    ) -> None:
        self._rewrite_active_policy_as_historical_predecessor()
        historical_policy_sha256 = sha256(
            self.pic_root / "policy" / "storage_policy.json"
        )
        historical_promotion_sha256 = sha256(
            self.pic_root / "policy" / "active_promotion.json"
        )
        with self.assertRaisesRegex(ValueError, "lacks authenticated"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )
        with self.assertRaisesRegex(ValueError, "lacks authenticated"):
            promote(
                self.policy,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with patch(
            "control_plane_common."
            "AUTHORIZED_HISTORICAL_STORAGE_PREFLIGHT_RETIREMENT_POLICY_SHA256",
            historical_policy_sha256,
        ), patch(
            "control_plane_common."
            "AUTHORIZED_HISTORICAL_STORAGE_PREFLIGHT_RETIREMENT_PROMOTION_SHA256",
            historical_promotion_sha256,
        ):
            promote(
                self.policy,
                retire_historical_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        require_storage_policy_unlock_snapshot(
            control_plane_version=successor.name,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )
        with self.assertRaisesRegex(ValueError, "requires a legacy predecessor"):
            promote(
                self.policy,
                retire_historical_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_historical_storage_preflight_retirement_rejects_unreviewed_live_anchors(
        self,
    ) -> None:
        self._rewrite_active_policy_as_historical_predecessor()
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )
        with self.assertRaisesRegex(ValueError, "exact reviewed live anchors"):
            promote(
                self.policy,
                retire_historical_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_historical_storage_preflight_retirement_flag_requires_launch_prohibited_policy(
        self,
    ) -> None:
        self._write_policy(
            science_submission_freeze={
                "status": "authorized",
                "manifest_path": str(
                    self.pic_root
                    / "clean_candidates"
                    / str(uuid.uuid4())
                    / "clean_candidate_manifest.json"
                ),
                "manifest_sha256": "1" * 64,
                "build_profile_control_plane_version": self.control_plane_version,
            }
        )
        with self.assertRaisesRegex(ValueError, "launch-prohibited"):
            promote(
                self.policy,
                retire_historical_storage_preflight_predecessor=True,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_exact_reviewed_storage_preflight_predecessor_migration_is_one_time(
        self,
    ) -> None:
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )
        with self.assertRaisesRegex(ValueError, "source authentication"):
            promote(
                self.policy,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with self._authorize_active_exact_reviewed_preflight_predecessor():
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
            require_storage_policy_unlock_snapshot(
                control_plane_version=successor.name,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
            second_successor = self._publish_test_control_plane_successor(
                self.pic_root, schema_suffix="\n\n"
            )
            project_home_second_successor = self._publish_test_control_plane_successor(
                self.project_home_root, schema_suffix="\n\n"
            )
            self.assertEqual(project_home_second_successor.name, second_successor.name)
            self.assertNotEqual(second_successor.name, successor.name)
            self._write_policy(
                installed_control_plane_version=second_successor.name,
                staged_control_plane_candidate_version=second_successor.name,
            )
            with self.assertRaisesRegex(ValueError, "exact reviewed live anchors"):
                promote(
                    self.policy,
                    migrate_exact_reviewed_storage_preflight_predecessor=True,
                    control_plane_dir=second_successor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )

    def test_exact_reviewed_storage_preflight_migration_preserves_historical_authorized_freeze(
        self,
    ) -> None:
        candidate, _, _ = self._clean_candidate(
            authorize=True,
            source_name="migration-preserved-candidate-source",
            freeze_id=str(uuid.uuid4()),
            profile_id="hip-mpi-release-paper-pic-migration-preserved",
        )
        candidate_source_root = self.authorized_clean_candidate_source_root
        assert candidate_source_root is not None
        authorized_freeze = self._authorized_science_freeze(candidate)
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
            science_submission_freeze=authorized_freeze,
        )

        def revalidate_with_test_source(
            candidate_manifest_path: Path, **kwargs: object
        ) -> dict[str, object]:
            return revalidate_clean_candidate.revalidate_clean_candidate(
                candidate_manifest_path,
                **kwargs,
                authorized_source_root=candidate_source_root,
            )

        with self._authorize_active_exact_reviewed_preflight_predecessor(), patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=revalidate_with_test_source,
        ) as revalidate:
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(revalidate.call_count, 3)
        for invocation in revalidate.call_args_list:
            self.assertEqual(
                invocation.kwargs["expected_receipt_control_plane_version"],
                self.control_plane_version,
            )
            self.assertEqual(invocation.kwargs["control_plane_dir"], successor)
        active, _ = require_storage_policy_unlock_snapshot(
            control_plane_version=successor.name,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        self.assertEqual(active["science_submission_freeze"], authorized_freeze)

    def test_exact_reviewed_storage_preflight_migration_revalidates_preserved_freeze_before_publication(
        self,
    ) -> None:
        candidate, _, _ = self._clean_candidate(authorize=True)
        authorized_freeze = self._authorized_science_freeze(candidate)
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
            science_submission_freeze=authorized_freeze,
        )
        anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        before = {path: path.read_bytes() for path in anchors}

        with self._authorize_active_exact_reviewed_preflight_predecessor(), patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=ValueError("preserved candidate differs"),
        ), self.assertRaisesRegex(ValueError, "preserved candidate differs"):
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual({path: path.read_bytes() for path in anchors}, before)
        for root in [self.pic_root, self.project_home_root]:
            self.assertFalse(
                (root / "policy" / ".active_promotion_transaction.json").exists()
            )
            self.assertEqual(
                list((root / "policy").glob(".*.transaction-rollback-*")),
                [],
            )

    def test_exact_reviewed_storage_preflight_migration_committed_recovery_preserves_historical_authorized_freeze(
        self,
    ) -> None:
        candidate, _, _ = self._clean_candidate(
            authorize=True,
            source_name="migration-recovery-preserved-candidate-source",
            freeze_id=str(uuid.uuid4()),
            profile_id="hip-mpi-release-paper-pic-migration-recovery-preserved",
        )
        candidate_source_root = self.authorized_clean_candidate_source_root
        assert candidate_source_root is not None
        authorized_freeze = self._authorized_science_freeze(candidate)
        historical_build_controller = self.control_plane_version
        active_predecessor = self._publish_test_control_plane_successor(
            self.pic_root, schema_suffix="\n\n"
        )
        project_home_active_predecessor = self._publish_test_control_plane_successor(
            self.project_home_root, schema_suffix="\n\n"
        )
        self.assertEqual(
            project_home_active_predecessor.name,
            active_predecessor.name,
        )
        self.assertNotEqual(active_predecessor.name, historical_build_controller)
        self._write_policy(
            installed_control_plane_version=active_predecessor.name,
            staged_control_plane_candidate_version=active_predecessor.name,
            science_submission_freeze=authorized_freeze,
        )

        def revalidate_with_test_source(
            candidate_manifest_path: Path, **kwargs: object
        ) -> dict[str, object]:
            return revalidate_clean_candidate.revalidate_clean_candidate(
                candidate_manifest_path,
                **kwargs,
                authorized_source_root=candidate_source_root,
            )

        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=revalidate_with_test_source,
        ) as predecessor_revalidate:
            promote(
                self.policy,
                control_plane_dir=active_predecessor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(predecessor_revalidate.call_count, 2)
        for invocation in predecessor_revalidate.call_args_list:
            self.assertEqual(
                invocation.kwargs["expected_receipt_control_plane_version"],
                historical_build_controller,
            )
            self.assertEqual(invocation.kwargs["control_plane_dir"], active_predecessor)
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        successor = self._publish_test_control_plane_successor(
            self.pic_root, schema_suffix="\n\n\n"
        )
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root, schema_suffix="\n\n\n"
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self.assertNotEqual(successor.name, active_predecessor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
            science_submission_freeze=authorized_freeze,
        )
        real_unlink = promote_active_policy._unlink_if_exists_at
        cleanup_failed = False

        def fail_first_committed_marker_cleanup(
            parent_descriptor: int, name: str
        ) -> None:
            nonlocal cleanup_failed
            marker_path = (
                Path("/proc")
                / str(os.getpid())
                / "fd"
                / str(parent_descriptor)
                / name
            )
            if (
                not cleanup_failed
                and name == ".active_promotion_transaction.json"
                and marker_path.exists()
                and json.loads(marker_path.read_text(encoding="utf-8")).get("state")
                == "committed"
            ):
                cleanup_failed = True
                raise OSError("injected exact migration committed cleanup failure")
            real_unlink(parent_descriptor, name)

        with self._authorize_active_exact_reviewed_preflight_predecessor(
            control_plane_version=active_predecessor.name
        ), patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=revalidate_with_test_source,
        ) as migration_revalidate, patch(
            "promote_active_policy._unlink_if_exists_at",
            side_effect=fail_first_committed_marker_cleanup,
        ):
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertTrue(cleanup_failed)
        self.assertEqual(migration_revalidate.call_count, 3)
        marker_paths = [
            root / "policy" / ".active_promotion_transaction.json"
            for root in [self.pic_root, self.project_home_root]
        ]
        self.assertTrue(all(path.exists() for path in marker_paths))
        self.assertTrue(
            all(
                json.loads(path.read_text(encoding="utf-8"))["state"] == "committed"
                for path in marker_paths
            )
        )
        with self.assertRaisesRegex(ValueError, "requires locked recovery"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=successor.name,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

        real_predecessor_validation = (
            promote_active_policy.require_policy_predecessor_snapshot_for_promotion
        )

        def fail_after_committed_recovery(*args: object, **kwargs: object) -> object:
            if kwargs.get("allow_active_promotion_transaction") is True:
                return real_predecessor_validation(*args, **kwargs)
            raise ValueError("injected after historical-freeze committed recovery")

        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=revalidate_with_test_source,
        ) as recovery_revalidate, patch(
            "promote_active_policy.require_policy_predecessor_snapshot_for_promotion",
            side_effect=fail_after_committed_recovery,
        ), self.assertRaisesRegex(
            ValueError,
            "injected after historical-freeze committed recovery",
        ):
            promote(
                self.policy,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(recovery_revalidate.call_count, 1)
        recovery_invocation = recovery_revalidate.call_args
        assert recovery_invocation is not None
        self.assertEqual(
            recovery_invocation.kwargs["expected_receipt_control_plane_version"],
            historical_build_controller,
        )
        self.assertEqual(recovery_invocation.kwargs["control_plane_dir"], successor)
        for root in [self.pic_root, self.project_home_root]:
            self.assertFalse(
                (root / "policy" / ".active_promotion_transaction.json").exists()
            )
            self.assertEqual(
                list((root / "policy").glob(".*.transaction-rollback-*")),
                [],
            )
        active, _ = require_storage_policy_unlock_snapshot(
            control_plane_version=successor.name,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        self.assertEqual(active["science_submission_freeze"], authorized_freeze)
        self.assertEqual(
            active["science_submission_freeze"][
                "build_profile_control_plane_version"
            ],
            historical_build_controller,
        )
        self.assertEqual(
            active["olcf_side_storage"]["installed_control_plane_version"],
            successor.name,
        )
        self.assertEqual(active["registered_science_slices"], [])

    def test_exact_reviewed_storage_preflight_migration_complete_prepared_or_mixed_recovery_preserves_historical_authorized_freeze(
        self,
    ) -> None:
        candidate, _, _ = self._clean_candidate(
            authorize=True,
            source_name="migration-prepared-recovery-preserved-candidate-source",
            freeze_id=str(uuid.uuid4()),
            profile_id="hip-mpi-release-paper-pic-migration-prepared-recovery",
        )
        candidate_source_root = self.authorized_clean_candidate_source_root
        assert candidate_source_root is not None
        authorized_freeze = self._authorized_science_freeze(candidate)
        historical_build_controller = self.control_plane_version
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        predecessor = {path: path.read_bytes() for _, path in anchors}
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
            science_submission_freeze=authorized_freeze,
        )

        def revalidate_with_test_source(
            candidate_manifest_path: Path, **kwargs: object
        ) -> dict[str, object]:
            return revalidate_clean_candidate.revalidate_clean_candidate(
                candidate_manifest_path,
                **kwargs,
                authorized_source_root=candidate_source_root,
            )

        with self._authorize_active_exact_reviewed_preflight_predecessor(), patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=revalidate_with_test_source,
        ):
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        successor_anchors = {path: path.read_bytes() for _, path in anchors}
        self.assertNotEqual(successor_anchors, predecessor)
        real_predecessor_validation = (
            promote_active_policy.require_policy_predecessor_snapshot_for_promotion
        )

        def fail_after_complete_successor_recovery(
            *args: object, **kwargs: object
        ) -> object:
            if kwargs.get("allow_active_promotion_transaction") is True:
                return real_predecessor_validation(*args, **kwargs)
            raise ValueError(
                "injected after historical-freeze complete-successor recovery"
            )

        for states in [("prepared", "prepared"), ("committed", "prepared")]:
            with self.subTest(states=states):
                transaction_id = str(uuid.uuid4())
                marker = {
                    "schema_version": (
                        promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION
                    ),
                    "record_type": (
                        promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE
                    ),
                    "transaction_id": transaction_id,
                    "state": "prepared",
                    "predecessor_state": "complete",
                    "anchors": [],
                }
                rollback_paths = []
                for root_role, path in anchors:
                    rollback_name = (
                        f".{path.name}.transaction-rollback-{transaction_id}"
                    )
                    rollback_path = path.with_name(rollback_name)
                    rollback_path.write_bytes(predecessor[path])
                    rollback_path.chmod(0o444)
                    rollback_paths.append(rollback_path)
                    marker["anchors"].append(
                        {
                            "root_role": root_role,
                            "name": path.name,
                            "rollback_name": rollback_name,
                            "predecessor_sha256": hashlib.sha256(
                                predecessor[path]
                            ).hexdigest(),
                            "successor_sha256": hashlib.sha256(
                                successor_anchors[path]
                            ).hexdigest(),
                        }
                    )
                marker_paths = []
                for state, root in zip(
                    states,
                    [self.pic_root, self.project_home_root],
                ):
                    marker_path = (
                        root / "policy" / ".active_promotion_transaction.json"
                    )
                    marker_path.write_text(
                        json.dumps(
                            {**marker, "state": state},
                            indent=2,
                            sort_keys=True,
                        )
                        + "\n",
                        encoding="utf-8",
                    )
                    marker_path.chmod(0o400)
                    marker_paths.append(marker_path)

                with patch(
                    "promote_active_policy.revalidate_clean_candidate",
                    side_effect=revalidate_with_test_source,
                ) as recovery_revalidate, patch(
                    "promote_active_policy.require_policy_predecessor_snapshot_for_promotion",
                    side_effect=fail_after_complete_successor_recovery,
                ), self.assertRaisesRegex(
                    ValueError,
                    "injected after historical-freeze complete-successor recovery",
                ):
                    promote(
                        self.policy,
                        control_plane_dir=successor,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                    )
                self.assertEqual(recovery_revalidate.call_count, 1)
                recovery_invocation = recovery_revalidate.call_args
                assert recovery_invocation is not None
                self.assertEqual(
                    recovery_invocation.kwargs[
                        "expected_receipt_control_plane_version"
                    ],
                    historical_build_controller,
                )
                self.assertEqual(
                    recovery_invocation.kwargs["control_plane_dir"],
                    successor,
                )
                self.assertEqual(
                    {path: path.read_bytes() for _, path in anchors},
                    successor_anchors,
                )
                self.assertFalse(any(path.exists() for path in marker_paths))
                self.assertFalse(any(path.exists() for path in rollback_paths))
                active, _ = require_storage_policy_unlock_snapshot(
                    control_plane_version=successor.name,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
                self.assertEqual(active["science_submission_freeze"], authorized_freeze)
                self.assertEqual(active["registered_science_slices"], [])

    def test_exact_reviewed_storage_preflight_predecessor_rejects_source_drift(
        self,
    ) -> None:
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        self._rewrite_active_storage_preflight_artifact(
            lambda artifact: artifact["source_authentication"].update(
                runner_sha256="0" * 64
            )
        )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )
        with self._authorize_active_exact_reviewed_preflight_predecessor(), (
            self.assertRaisesRegex(ValueError, "source authentication")
        ):
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_exact_reviewed_storage_preflight_predecessor_rejects_wrong_anchor(
        self,
    ) -> None:
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )
        with self.assertRaisesRegex(ValueError, "exact reviewed live anchors"):
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_exact_reviewed_storage_preflight_predecessor_requires_new_controller(
        self,
    ) -> None:
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        self._write_policy()
        with self._authorize_active_exact_reviewed_preflight_predecessor(), (
            self.assertRaisesRegex(ValueError, "requires a new control plane")
        ):
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_exact_reviewed_storage_preflight_predecessor_rejects_freeze_drift(
        self,
    ) -> None:
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            science_submission_freeze={
                "status": "authorized",
                "manifest_path": str(
                    self.pic_root
                    / "clean_candidates"
                    / str(uuid.uuid4())
                    / "clean_candidate_manifest.json"
                ),
                "manifest_sha256": "1" * 64,
                "build_profile_control_plane_version": self.control_plane_version,
            },
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )
        with self._authorize_active_exact_reviewed_preflight_predecessor(), (
            self.assertRaisesRegex(ValueError, "exact authorized predecessor transformation")
        ):
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_exact_reviewed_storage_preflight_predecessor_rejects_metadata_drift(
        self,
    ) -> None:
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )
        policy = json.loads(self.policy.read_text(encoding="utf-8"))
        policy["reviewer"] = "unrelated successor drift"
        self.policy.write_text(json.dumps(policy), encoding="utf-8")
        with self._authorize_active_exact_reviewed_preflight_predecessor(), (
            self.assertRaisesRegex(ValueError, "exact authorized predecessor transformation")
        ):
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_exact_reviewed_storage_preflight_predecessor_rejects_numeric_type_drift(
        self,
    ) -> None:
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )
        policy = json.loads(self.policy.read_text(encoding="utf-8"))
        self.assertIsInstance(policy["frontier"]["maximum_node_hours"], float)
        policy["frontier"]["maximum_node_hours"] = 10000
        self.policy.write_text(json.dumps(policy), encoding="utf-8")
        with self._authorize_active_exact_reviewed_preflight_predecessor(), (
            self.assertRaisesRegex(ValueError, "exact authorized predecessor transformation")
        ):
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_exact_reviewed_storage_preflight_predecessor_requires_newer_preflight(
        self,
    ) -> None:
        self._rewrite_active_policy_as_exact_reviewed_preflight_predecessor()
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )
        self._rewrite_reviewed_storage_preflight_artifact(
            lambda artifact: artifact.update(
                completed_utc="2026-06-02T00:00:00Z",
                started_utc="2026-06-02T00:00:00Z",
            )
        )
        policy = json.loads(self.policy.read_text(encoding="utf-8"))
        policy["olcf_side_storage"]["last_preflight_utc"] = "2026-06-02T00:00:00Z"
        self.policy.write_text(json.dumps(policy), encoding="utf-8")
        with self._authorize_active_exact_reviewed_preflight_predecessor(), (
            self.assertRaisesRegex(ValueError, "requires a newer authenticated")
        ):
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=successor,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_exact_authorized_clean_candidate_freeze_replacement_is_compare_and_swap(
        self,
    ) -> None:
        first_freeze = {
            "status": "authorized",
            "manifest_path": str(
                self.pic_root
                / "clean_candidates"
                / "first"
                / "clean_candidate_manifest.json"
            ),
            "manifest_sha256": "1" * 64,
            "build_profile_control_plane_version": self.control_plane_version,
        }
        self._write_policy(science_submission_freeze=first_freeze)
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            return_value={
                "status": "passed",
                "current_control_plane_version": self.control_plane_version,
                "build": {
                    "receipt_control_plane_version": self.control_plane_version,
                },
            },
        ):
            self._promote_policy(patch_clean_candidate_revalidation=False)

        active_policy_path = self.pic_root / "policy" / "storage_policy.json"
        active_promotion_path = self.pic_root / "policy" / "active_promotion.json"
        expected_policy_sha256 = sha256(active_policy_path)
        expected_promotion_sha256 = sha256(active_promotion_path)
        predecessor = json.loads(active_policy_path.read_text(encoding="utf-8"))
        second = copy.deepcopy(predecessor)
        second["science_submission_freeze"] = {
            **first_freeze,
            "manifest_path": str(
                self.pic_root
                / "clean_candidates"
                / "second"
                / "clean_candidate_manifest.json"
            ),
            "manifest_sha256": "2" * 64,
        }
        stale = copy.deepcopy(predecessor)
        stale["science_submission_freeze"] = {
            **first_freeze,
            "manifest_path": str(
                self.pic_root
                / "clean_candidates"
                / "stale"
                / "clean_candidate_manifest.json"
            ),
            "manifest_sha256": "3" * 64,
        }
        self.policy.write_text(json.dumps(second), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "exact active-predecessor transition"):
            self._promote_policy()
        with self.assertRaisesRegex(ValueError, "requires both exact active predecessor"):
            promote(
                self.policy,
                replace_exact_authorized_clean_candidate_freeze=True,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            return_value={
                "status": "passed",
                "current_control_plane_version": self.control_plane_version,
                "build": {
                    "receipt_control_plane_version": self.control_plane_version,
                },
            },
        ) as revalidate:
            promote(
                self.policy,
                replace_exact_authorized_clean_candidate_freeze=True,
                expected_active_policy_sha256=expected_policy_sha256,
                expected_active_promotion_sha256=expected_promotion_sha256,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        expected_revalidation = call(
            Path(second["science_submission_freeze"]["manifest_path"]),
            expected_manifest_sha256=second["science_submission_freeze"][
                "manifest_sha256"
            ],
            expected_receipt_control_plane_version=self.control_plane_version,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        self.assertEqual(revalidate.call_args_list, [expected_revalidation] * 3)
        self.assertEqual(
            json.loads(active_policy_path.read_text(encoding="utf-8")),
            second,
        )

        second_policy_sha256 = sha256(active_policy_path)
        second_promotion_sha256 = sha256(active_promotion_path)
        self.policy.write_text(json.dumps(predecessor), encoding="utf-8")
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            return_value={
                "status": "passed",
                "current_control_plane_version": self.control_plane_version,
                "build": {
                    "receipt_control_plane_version": self.control_plane_version,
                },
            },
        ):
            promote(
                self.policy,
                replace_exact_authorized_clean_candidate_freeze=True,
                expected_active_policy_sha256=second_policy_sha256,
                expected_active_promotion_sha256=second_promotion_sha256,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(sha256(active_policy_path), expected_policy_sha256)
        self.assertNotEqual(sha256(active_promotion_path), expected_promotion_sha256)

        self.policy.write_text(json.dumps(stale), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "active predecessor hashes changed"):
            promote(
                self.policy,
                replace_exact_authorized_clean_candidate_freeze=True,
                expected_active_policy_sha256=expected_policy_sha256,
                expected_active_promotion_sha256=expected_promotion_sha256,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_authorized_clean_candidate_freeze_cannot_step_down_without_exact_mode(
        self,
    ) -> None:
        authorized_freeze = {
            "status": "authorized",
            "manifest_path": str(
                self.pic_root
                / "clean_candidates"
                / "first"
                / "clean_candidate_manifest.json"
            ),
            "manifest_sha256": "1" * 64,
            "build_profile_control_plane_version": self.control_plane_version,
        }
        self._write_policy(science_submission_freeze=authorized_freeze)
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            return_value={
                "status": "passed",
                "current_control_plane_version": self.control_plane_version,
                "build": {
                    "receipt_control_plane_version": self.control_plane_version,
                },
            },
        ):
            self._promote_policy(patch_clean_candidate_revalidation=False)
        anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        before = {path: path.read_bytes() for path in anchors}

        self._write_policy(
            science_submission_freeze={"status": "pending_clean_candidate_freeze"}
        )
        with self.assertRaisesRegex(ValueError, "exact active-predecessor transition"):
            self._promote_policy()
        self.assertEqual({path: path.read_bytes() for path in anchors}, before)

    def test_exact_authorized_clean_candidate_freeze_revalidation_failure_is_atomic(
        self,
    ) -> None:
        first_freeze = {
            "status": "authorized",
            "manifest_path": str(
                self.pic_root
                / "clean_candidates"
                / "first"
                / "clean_candidate_manifest.json"
            ),
            "manifest_sha256": "1" * 64,
            "build_profile_control_plane_version": self.control_plane_version,
        }
        self._write_policy(science_submission_freeze=first_freeze)
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            return_value={
                "status": "passed",
                "current_control_plane_version": self.control_plane_version,
                "build": {
                    "receipt_control_plane_version": self.control_plane_version,
                },
            },
        ):
            self._promote_policy(patch_clean_candidate_revalidation=False)
        anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        before = {path: path.read_bytes() for path in anchors}
        successor = json.loads(anchors[0].read_text(encoding="utf-8"))
        successor["science_submission_freeze"] = {
            **first_freeze,
            "manifest_path": str(
                self.pic_root
                / "clean_candidates"
                / "second"
                / "clean_candidate_manifest.json"
            ),
            "manifest_sha256": "2" * 64,
        }
        self.policy.write_text(json.dumps(successor), encoding="utf-8")
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=ValueError("candidate bundle differs"),
        ), self.assertRaisesRegex(ValueError, "candidate bundle differs"):
            promote(
                self.policy,
                replace_exact_authorized_clean_candidate_freeze=True,
                expected_active_policy_sha256=sha256(anchors[0]),
                expected_active_promotion_sha256=sha256(anchors[2]),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual({path: path.read_bytes() for path in anchors}, before)

    def test_authorized_clean_candidate_revalidation_uses_bound_historical_build_controller(
        self,
    ) -> None:
        historical_build_controller = "0" * 64
        manifest_path = (
            self.pic_root
            / "clean_candidates"
            / str(uuid.uuid4())
            / "clean_candidate_manifest.json"
        )
        policy = {
            "science_submission_freeze": {
                "status": "authorized",
                "manifest_path": str(manifest_path),
                "manifest_sha256": "1" * 64,
                "build_profile_control_plane_version": historical_build_controller,
            }
        }
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            return_value={
                "status": "passed",
                "current_control_plane_version": self.control_plane_version,
                "build": {
                    "receipt_control_plane_version": historical_build_controller,
                },
            },
        ) as revalidate:
            promote_active_policy._require_authorized_clean_candidate_freeze_revalidation(
                policy,
                control_plane_version=self.control_plane_version,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        revalidate.assert_called_once_with(
            manifest_path,
            expected_manifest_sha256="1" * 64,
            expected_receipt_control_plane_version=historical_build_controller,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )

    def test_active_launch_prohibited_generation_verifier_binds_exact_state(
        self,
    ) -> None:
        candidate, _, _ = self._clean_candidate(authorize=True)
        candidate_source_root = self.authorized_clean_candidate_source_root
        assert candidate_source_root is not None
        active_policy = self.pic_root / "policy" / "storage_policy.json"
        active_promotion = self.pic_root / "policy" / "active_promotion.json"
        freeze = self._authorized_science_freeze(candidate)
        serialization_lock_held = False

        def revalidate_with_test_source(
            candidate_manifest_path: Path, **kwargs: object
        ) -> dict[str, object]:
            nonlocal serialization_lock_held
            anchor = stable_serialization_anchor(self.pic_root)
            descriptor = os.open(
                anchor, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
            )
            try:
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                serialization_lock_held = True
            finally:
                os.close(descriptor)
            return revalidate_clean_candidate.revalidate_clean_candidate(
                candidate_manifest_path,
                **kwargs,
                authorized_source_root=candidate_source_root,
            )

        bindings = {
            "expected_control_plane_version": self.control_plane_version,
            "expected_active_policy_sha256": sha256(active_policy),
            "expected_active_promotion_sha256": sha256(active_promotion),
            "expected_authorized_freeze_manifest": candidate,
            "expected_authorized_freeze_manifest_sha256": freeze["manifest_sha256"],
            "expected_authorized_freeze_build_controller": self.control_plane_version,
            "authorized_pic_root": self.pic_root,
            "authorized_project_home_root": self.project_home_root,
        }
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=revalidate_with_test_source,
        ) as revalidate, patch.object(
            promote_active_policy,
            "SCRIPT_DIR",
            self.control_plane_dir,
        ):
            report = promote_active_policy.verify_active_launch_prohibited_generation(
                **bindings,
            )
        self.assertEqual(revalidate.call_count, 1)
        self.assertTrue(serialization_lock_held)
        self.assertEqual(
            report,
            {
                "active_policy_sha256": sha256(active_policy),
                "active_promotion_sha256": sha256(active_promotion),
                "control_plane_version": self.control_plane_version,
                "record_type": (
                    "frontier_pic_active_launch_prohibited_generation_verification"
                ),
                "schema_version": 1,
                "science_submission_freeze": freeze,
                "status": "passed",
            },
        )

        for label, override in [
            ("policy", {"expected_active_policy_sha256": "0" * 64}),
            ("promotion", {"expected_active_promotion_sha256": "0" * 64}),
            (
                "freeze",
                {"expected_authorized_freeze_manifest_sha256": "0" * 64},
            ),
        ]:
            with self.subTest(label=label), patch.object(
                promote_active_policy,
                "SCRIPT_DIR",
                self.control_plane_dir,
            ), self.assertRaisesRegex(
                ValueError,
                "differs from expected",
            ):
                promote_active_policy.verify_active_launch_prohibited_generation(
                    **{**bindings, **override},
                )

        rollback_anchor = (
            self.pic_root
            / "policy"
            / f".storage_policy.json.transaction-rollback-{uuid.uuid4()}"
        )
        rollback_anchor.write_bytes(active_policy.read_bytes())
        rollback_anchor.chmod(0o400)
        with patch.object(
            promote_active_policy,
            "SCRIPT_DIR",
            self.control_plane_dir,
        ), self.assertRaisesRegex(ValueError, "requires locked recovery"):
            promote_active_policy.verify_active_launch_prohibited_generation(
                **bindings,
            )

        with self.assertRaisesRegex(TypeError, "control_plane_dir"):
            promote_active_policy.verify_active_launch_prohibited_generation(
                **bindings,
                control_plane_dir=self.control_plane_dir,
            )

    def test_active_launch_prohibited_generation_verifier_cli_dispatches_exact_bindings(
        self,
    ) -> None:
        expected = {
            "schema_version": 1,
            "status": "passed",
        }
        manifest = Path("/tmp/clean_candidate_manifest.json")
        arguments = [
            "promote_active_policy.py",
            "--verify-active-launch-prohibited-generation",
            "--expected-control-plane-version",
            "a" * 64,
            "--expected-active-policy-sha256",
            "b" * 64,
            "--expected-active-promotion-sha256",
            "c" * 64,
            "--expected-authorized-freeze-manifest",
            str(manifest),
            "--expected-authorized-freeze-manifest-sha256",
            "d" * 64,
            "--expected-authorized-freeze-build-controller",
            "e" * 64,
        ]
        output = io.StringIO()
        with patch.object(
            promote_active_policy,
            "verify_active_launch_prohibited_generation",
            return_value=expected,
        ) as verify_generation, patch.object(
            sys, "argv", arguments
        ), patch.object(
            sys, "stdout", output
        ):
            promote_active_policy.main()
        verify_generation.assert_called_once_with(
            expected_active_policy_sha256="b" * 64,
            expected_active_promotion_sha256="c" * 64,
            expected_control_plane_version="a" * 64,
            expected_authorized_freeze_manifest=manifest,
            expected_authorized_freeze_manifest_sha256="d" * 64,
            expected_authorized_freeze_build_controller="e" * 64,
        )
        self.assertEqual(
            output.getvalue(),
            json.dumps(expected, separators=(",", ":"), sort_keys=True) + "\n",
        )

    def test_exact_authorized_clean_candidate_mutation_after_initial_revalidation_requires_recovery(
        self,
    ) -> None:
        self._clean_candidate(authorize=True)
        anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        before = {path: path.read_bytes() for path in anchors}
        candidate, _, _ = self._clean_candidate(
            authorize=False,
            source_name="replacement-candidate-toctou-precommit",
            freeze_id=str(uuid.uuid4()),
            profile_id="hip-mpi-release-paper-pic-toctou-precommit",
        )
        replacement_source_root = self.authorized_clean_candidate_source_root
        assert replacement_source_root is not None
        successor = json.loads(anchors[0].read_text(encoding="utf-8"))
        successor["science_submission_freeze"] = self._authorized_science_freeze(
            candidate
        )
        self.policy.write_text(json.dumps(successor), encoding="utf-8")
        reviewed_successor = copy.deepcopy(successor)
        validations = 0

        def revalidate_then_mutate(
            candidate_manifest_path: Path, **kwargs: object
        ) -> dict[str, object]:
            nonlocal validations
            result = revalidate_clean_candidate.revalidate_clean_candidate(
                candidate_manifest_path,
                **kwargs,
                authorized_source_root=replacement_source_root,
            )
            validations += 1
            if validations == 1:
                candidate.parent.chmod(0o755)
                candidate.chmod(0o644)
                candidate.write_bytes(candidate.read_bytes() + b"\n")
                candidate.chmod(0o444)
                candidate.parent.chmod(0o555)
            return result

        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=revalidate_then_mutate,
        ), self.assertRaisesRegex(
            RuntimeError,
            "Committed active-policy promotion requires locked recovery",
        ):
            promote(
                self.policy,
                replace_exact_authorized_clean_candidate_freeze=True,
                expected_active_policy_sha256=sha256(anchors[0]),
                expected_active_promotion_sha256=sha256(anchors[2]),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertEqual(validations, 1)
        self.assertNotEqual({path: path.read_bytes() for path in anchors}, before)
        self.assertEqual(
            json.loads(anchors[0].read_text(encoding="utf-8")),
            reviewed_successor,
        )
        self.assertEqual(anchors[0].read_bytes(), anchors[1].read_bytes())
        self.assertEqual(anchors[2].read_bytes(), anchors[3].read_bytes())
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            self.assertTrue(marker_path.exists())
            self.assertEqual(
                json.loads(marker_path.read_text(encoding="utf-8"))["state"],
                "prepared",
            )
            self.assertEqual(
                len(list((root / "policy").glob(".*.transaction-rollback-*"))),
                2,
            )

    def test_authorized_clean_candidate_mutation_after_committed_markers_requires_recovery(
        self,
    ) -> None:
        candidate, _, _ = self._clean_candidate(
            authorize=False,
            source_name="candidate-toctou-committed",
            freeze_id=str(uuid.uuid4()),
            profile_id="hip-mpi-release-paper-pic-toctou-committed",
        )
        candidate_source_root = self.authorized_clean_candidate_source_root
        assert candidate_source_root is not None
        self._write_policy(
            science_submission_freeze=self._authorized_science_freeze(candidate)
        )
        reviewed_policy = self.policy.read_bytes()
        real_atomic_write_json_at = promote_active_policy.atomic_write_json_at
        committed_markers = 0

        def mutate_after_second_committed_marker(
            parent_descriptor: int,
            name: str,
            value: object,
            **kwargs: object,
        ) -> None:
            nonlocal committed_markers
            real_atomic_write_json_at(parent_descriptor, name, value, **kwargs)
            if (
                name == ".active_promotion_transaction.json"
                and isinstance(value, dict)
                and value.get("state") == "committed"
            ):
                committed_markers += 1
                if committed_markers == 2:
                    candidate.parent.chmod(0o755)
                    candidate.chmod(0o644)
                    candidate.write_bytes(candidate.read_bytes() + b"\n")
                    candidate.chmod(0o444)
                    candidate.parent.chmod(0o555)

        def revalidate_with_test_source(
            candidate_manifest_path: Path, **kwargs: object
        ) -> dict[str, object]:
            return revalidate_clean_candidate.revalidate_clean_candidate(
                candidate_manifest_path,
                **kwargs,
                authorized_source_root=candidate_source_root,
            )

        with patch(
            "promote_active_policy.atomic_write_json_at",
            side_effect=mutate_after_second_committed_marker,
        ), patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=revalidate_with_test_source,
        ), self.assertRaisesRegex(
            RuntimeError,
            "Committed active-policy promotion requires locked recovery",
        ):
            self._promote_policy(patch_clean_candidate_revalidation=False)
        self.assertEqual(committed_markers, 2)
        self.assertEqual(
            (self.pic_root / "policy" / "storage_policy.json").read_bytes(),
            reviewed_policy,
        )
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            self.assertTrue(marker_path.exists())
            self.assertEqual(
                json.loads(marker_path.read_text(encoding="utf-8"))["state"],
                "committed",
            )
            self.assertEqual(
                len(list((root / "policy").glob(".*.transaction-rollback-*"))),
                2,
            )
        with self.assertRaisesRegex(ValueError, "requires locked recovery"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                allow_pending_genesis=True,
            )

    def test_committed_authorized_clean_candidate_recovery_preserves_evidence_when_candidate_is_invalid(
        self,
    ) -> None:
        candidate, _, _ = self._clean_candidate(
            authorize=True,
            source_name="candidate-toctou-recovery",
            freeze_id=str(uuid.uuid4()),
            profile_id="hip-mpi-release-paper-pic-toctou-recovery",
        )
        candidate_source_root = self.authorized_clean_candidate_source_root
        assert candidate_source_root is not None
        transaction_id = str(uuid.uuid4())
        anchors = [
            ("orion", self.pic_root / "policy" / "storage_policy.json"),
            ("orion", self.pic_root / "policy" / "active_promotion.json"),
            (
                "project_home",
                self.project_home_root / "policy" / "storage_policy.json",
            ),
            (
                "project_home",
                self.project_home_root / "policy" / "active_promotion.json",
            ),
        ]
        marker = {
            "schema_version": promote_active_policy.PROMOTION_TRANSACTION_SCHEMA_VERSION,
            "record_type": promote_active_policy.PROMOTION_TRANSACTION_RECORD_TYPE,
            "transaction_id": transaction_id,
            "state": "committed",
            "predecessor_state": "complete",
            "anchors": [],
        }
        rollback_paths = []
        for root_role, path in anchors:
            rollback_name = f".{path.name}.transaction-rollback-{transaction_id}"
            rollback_path = path.with_name(rollback_name)
            os.link(path, rollback_path)
            rollback_paths.append(rollback_path)
            marker["anchors"].append(
                {
                    "root_role": root_role,
                    "name": path.name,
                    "rollback_name": rollback_name,
                    "predecessor_sha256": sha256(path),
                    "successor_sha256": sha256(path),
                }
            )
        marker_payload = json.dumps(marker, indent=2, sort_keys=True) + "\n"
        marker_paths = []
        for root in [self.pic_root, self.project_home_root]:
            marker_path = root / "policy" / ".active_promotion_transaction.json"
            marker_path.write_text(marker_payload, encoding="utf-8")
            marker_path.chmod(0o400)
            marker_paths.append(marker_path)
        candidate.parent.chmod(0o755)
        candidate.chmod(0o644)
        candidate.write_bytes(candidate.read_bytes() + b"\n")
        candidate.chmod(0o444)
        candidate.parent.chmod(0o555)

        def revalidate_with_test_source(
            candidate_manifest_path: Path, **kwargs: object
        ) -> dict[str, object]:
            return revalidate_clean_candidate.revalidate_clean_candidate(
                candidate_manifest_path,
                **kwargs,
                authorized_source_root=candidate_source_root,
            )

        self._write_policy()
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            side_effect=revalidate_with_test_source,
        ), self.assertRaisesRegex(ValueError, "differs from expected binding"):
            self._promote_policy(patch_clean_candidate_revalidation=False)
        self.assertTrue(all(path.exists() for path in marker_paths))
        self.assertTrue(all(path.exists() for path in rollback_paths))
        with self.assertRaisesRegex(ValueError, "requires locked recovery"):
            require_storage_policy_unlock_snapshot(
                control_plane_version=self.control_plane_version,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
                allow_pending_genesis=True,
            )

    def test_exact_authorized_clean_candidate_freeze_replacement_revalidates_real_malformed_candidates(
        self,
    ) -> None:
        self._clean_candidate(authorize=True)
        anchors = [
            self.pic_root / "policy" / "storage_policy.json",
            self.project_home_root / "policy" / "storage_policy.json",
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        before = {path: path.read_bytes() for path in anchors}
        expected_policy_sha256 = sha256(anchors[0])
        expected_promotion_sha256 = sha256(anchors[2])

        def rewrite_manifest(
            candidate: Path, mutate: Callable[[dict[str, object]], None]
        ) -> None:
            value = json.loads(candidate.read_text(encoding="utf-8"))
            mutate(value)
            candidate.parent.chmod(0o755)
            candidate.chmod(0o644)
            candidate.write_text(json.dumps(value), encoding="utf-8")
            candidate.chmod(0o444)
            candidate.parent.chmod(0o555)

        def rewrite_profile(
            candidate: Path, mutate: Callable[[dict[str, object]], None]
        ) -> None:
            value = json.loads(candidate.read_text(encoding="utf-8"))
            build = value["build"]
            assert isinstance(build, dict)
            profile_path = Path(str(build["profile_path"]))
            profile = json.loads(profile_path.read_text(encoding="utf-8"))
            mutate(profile)
            candidate.parent.chmod(0o755)
            profile_path.chmod(0o644)
            profile_path.write_text(json.dumps(profile), encoding="utf-8")
            profile_path.chmod(0o444)
            build["profile_sha256"] = sha256(profile_path)
            candidate.chmod(0o644)
            candidate.write_text(json.dumps(value), encoding="utf-8")
            candidate.chmod(0o444)
            candidate.parent.chmod(0o555)

        def mutate_prepared_artifact(candidate: Path) -> None:
            def mutate(value: dict[str, object]) -> None:
                prepared = value["prepared_artifacts"]
                assert isinstance(prepared, dict)
                paper_decks = prepared["paper_decks"]
                assert isinstance(paper_decks, list)
                record = paper_decks[0]
                assert isinstance(record, dict)
                record["sha256"] = "0" * 64

            rewrite_manifest(candidate, mutate)

        def mutate_receipt_controller(candidate: Path) -> None:
            value = json.loads(candidate.read_text(encoding="utf-8"))
            build = value["build"]
            assert isinstance(build, dict)
            receipt_path = Path(str(build["profile_receipt_path"]))
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt["control_plane_version"] = "0" * 64
            candidate.parent.chmod(0o755)
            receipt_path.chmod(0o644)
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            receipt_path.chmod(0o444)
            candidate.parent.chmod(0o555)

        def mutate_duplicate_candidate_keys(candidate: Path) -> None:
            candidate.parent.chmod(0o755)
            candidate.chmod(0o644)
            candidate.write_text(
                '{"schema_version":4,"schema_version":4}\n',
                encoding="utf-8",
            )
            candidate.chmod(0o444)
            candidate.parent.chmod(0o555)

        def mutate_layout(candidate: Path) -> None:
            rewrite_manifest(
                candidate,
                lambda value: value["build"].update(
                    executable_path=str(self.pic_root / "alternate-athena")
                ),
            )

        def mutate_profile_source_root(candidate: Path) -> None:
            rewrite_profile(
                candidate,
                lambda profile: profile.update(
                    authorized_source_root=str(self.root / "forged-source")
                ),
            )

        def mutate_profile_provenance_path(candidate: Path) -> None:
            def mutate(profile: dict[str, object]) -> None:
                provenance = profile["provenance_inputs"]
                assert isinstance(provenance, dict)
                toolchain = provenance["toolchain"]
                assert isinstance(toolchain, dict)
                toolchain["path"] = str(self.pic_root / "bin" / "forged-toolchain.txt")

            rewrite_profile(candidate, mutate)

        cases: list[tuple[str, Callable[[Path], None], str]] = [
            (
                "prepared-artifact manifest tamper",
                mutate_prepared_artifact,
                "differs from archived source inventory",
            ),
            (
                "different build-receipt controller",
                mutate_receipt_controller,
                "different control-plane version",
            ),
            (
                "duplicate candidate JSON key",
                mutate_duplicate_candidate_keys,
                "Duplicate JSON object key",
            ),
            (
                "candidate layout widening",
                mutate_layout,
                "does not match the fixed layout",
            ),
            (
                "forged build-profile source root",
                mutate_profile_source_root,
                "authorized source root",
            ),
            (
                "forged build-provenance path",
                mutate_profile_provenance_path,
                "documented Orion layout",
            ),
        ]
        for index, (label, mutate, error_pattern) in enumerate(cases):
            with self.subTest(label=label):
                candidate, _, _ = self._clean_candidate(
                    authorize=False,
                    source_name=f"replacement-candidate-source-{index}",
                    freeze_id=str(uuid.uuid4()),
                    profile_id=f"hip-mpi-release-paper-pic-replacement-{index}",
                )
                replacement_source_root = self.authorized_clean_candidate_source_root
                assert replacement_source_root is not None
                mutate(candidate)
                successor = json.loads(anchors[0].read_text(encoding="utf-8"))
                successor["science_submission_freeze"] = self._authorized_science_freeze(
                    candidate
                )
                self.policy.write_text(json.dumps(successor), encoding="utf-8")

                def revalidate_with_test_source(
                    candidate_manifest_path: Path, **kwargs: object
                ) -> dict[str, object]:
                    return revalidate_clean_candidate.revalidate_clean_candidate(
                        candidate_manifest_path,
                        **kwargs,
                        authorized_source_root=replacement_source_root,
                    )

                with patch(
                    "promote_active_policy.revalidate_clean_candidate",
                    side_effect=revalidate_with_test_source,
                ), self.assertRaisesRegex(ValueError, error_pattern):
                    promote(
                        self.policy,
                        replace_exact_authorized_clean_candidate_freeze=True,
                        expected_active_policy_sha256=expected_policy_sha256,
                        expected_active_promotion_sha256=expected_promotion_sha256,
                        control_plane_dir=self.control_plane_dir,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                    )
                self.assertEqual({path: path.read_bytes() for path in anchors}, before)

    def test_exact_authorized_clean_candidate_freeze_replacement_rejects_policy_drift(
        self,
    ) -> None:
        first_freeze = {
            "status": "authorized",
            "manifest_path": str(
                self.pic_root
                / "clean_candidates"
                / "first"
                / "clean_candidate_manifest.json"
            ),
            "manifest_sha256": "1" * 64,
            "build_profile_control_plane_version": self.control_plane_version,
        }
        self._write_policy(science_submission_freeze=first_freeze)
        with patch(
            "promote_active_policy.revalidate_clean_candidate",
            return_value={
                "status": "passed",
                "current_control_plane_version": self.control_plane_version,
                "build": {
                    "receipt_control_plane_version": self.control_plane_version,
                },
            },
        ):
            self._promote_policy(patch_clean_candidate_revalidation=False)
        active_policy_path = self.pic_root / "policy" / "storage_policy.json"
        active_promotion_path = self.pic_root / "policy" / "active_promotion.json"
        successor = json.loads(active_policy_path.read_text(encoding="utf-8"))
        successor["science_submission_freeze"] = {
            **first_freeze,
            "manifest_path": str(
                self.pic_root
                / "clean_candidates"
                / "second"
                / "clean_candidate_manifest.json"
            ),
            "manifest_sha256": "2" * 64,
        }
        successor["frontier"]["maximum_node_hours"] = 10000
        self.policy.write_text(json.dumps(successor), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "changed unrelated policy fields"):
            promote(
                self.policy,
                replace_exact_authorized_clean_candidate_freeze=True,
                expected_active_policy_sha256=sha256(active_policy_path),
                expected_active_promotion_sha256=sha256(active_promotion_path),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_policy_predecessor_transition_modes_are_exclusive(self) -> None:
        with self.assertRaisesRegex(ValueError, "modes are exclusive"):
            promote(
                self.policy,
                retire_historical_storage_preflight_predecessor=True,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with self.assertRaisesRegex(ValueError, "modes are exclusive"):
            promote(
                self.policy,
                migrate_exact_reviewed_storage_preflight_predecessor=True,
                replace_exact_authorized_clean_candidate_freeze=True,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_policy_rejects_unknown_nested_storage_fields_and_boolean_metadata_aliases(
        self,
    ) -> None:
        for mutate in [
            lambda policy: policy["olcf_side_storage"].update(unreviewed=True),
            lambda policy: policy["long_term_storage"].update(unreviewed=True),
            lambda policy: policy["olcf_side_storage"][
                "orion_simulation_root_preflight"
            ].update(unreviewed=True),
            lambda policy: policy["olcf_side_storage"].update(status=True),
            lambda policy: policy.update(reviewer=True),
        ]:
            with self.subTest(mutate=mutate):
                self._write_policy()
                policy = json.loads(self.policy.read_text(encoding="utf-8"))
                mutate(policy)
                self.policy.write_text(json.dumps(policy), encoding="utf-8")
                with self.assertRaises(ValueError):
                    self._promote_policy()

    def test_policy_rejects_digest_and_registered_text_coercion_aliases(self) -> None:
        integer_digest = int("1" * 64)
        for mutate in [
            lambda policy: policy["olcf_side_storage"]["ledger_genesis"].update(
                event_sha256=integer_digest
            ),
            lambda policy: policy["frontier_admission_smoke"].update(
                job_script_sha256=integer_digest
            ),
            lambda policy: policy["frontier_admission_smoke"].update(
                analysis_script_sha256=[integer_digest]
            ),
        ]:
            with self.subTest(mutate=mutate):
                self._write_policy()
                policy = json.loads(self.policy.read_text(encoding="utf-8"))
                mutate(policy)
                self.policy.write_text(json.dumps(policy), encoding="utf-8")
                with self.assertRaises(ValueError):
                    self._promote_policy()
        self._write_policy(
            science_submission_freeze={
                "status": "authorized",
                "manifest_path": str(
                    self.pic_root
                    / "clean_candidates"
                    / str(uuid.uuid4())
                    / "clean_candidate_manifest.json"
                ),
                "manifest_sha256": integer_digest,
                "build_profile_control_plane_version": self.control_plane_version,
            }
        )
        with self.assertRaises(ValueError):
            self._promote_policy()
        self._write_science_config(authorize=True)
        for field, value in [
            ("authorization_id", 1),
            ("campaign", True),
            ("job_script_sha256", integer_digest),
            ("analysis_script_sha256", [integer_digest]),
        ]:
            with self.subTest(field=field):
                self._write_policy(
                    science_submission_freeze=self.science_submission_freeze,
                    admission_smoke_overrides={"status": "closed_after_pass"},
                )
                policy = json.loads(self.policy.read_text(encoding="utf-8"))
                policy["registered_science_slices"][0][field] = value
                self.policy.write_text(json.dumps(policy), encoding="utf-8")
                with self.assertRaises(ValueError):
                    self._promote_policy()

    def test_policy_promotion_rejects_retained_genesis_rollback_or_fabrication(self) -> None:
        self._write_policy(ledger_genesis_allowed=True, ledger_genesis=None)
        with self.assertRaisesRegex(ValueError, "retained ledger bindings"):
            self._promote_policy()
        fabricated = {**self._closed_genesis, "event_sha256": "1" * 64}
        self._write_policy(ledger_genesis_allowed=False, ledger_genesis=fabricated)
        with self.assertRaisesRegex(ValueError, "retained ledger bindings"):
            self._promote_policy()

    def test_initializer_rejects_non_filesystem_copy_argument(self) -> None:
        with self.assertRaises(ValueError):
            initialize_from_policy(
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                mirror_receipts=self.receipts,
                mirror_jsonl=self.mirror,
                mirror_transport="dtn_rsync",
                notes="must not append",
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_initializer_rejects_closed_policy_after_complete_ledger_deletion(
        self,
    ) -> None:
        for anchor in genesis_anchor_paths(self.ledger, self.mirror):
            anchor.chmod(0o600)
            anchor.unlink()
        for path in [self.ledger, self.csv, self.receipts, self.mirror]:
            path.unlink()
        with self.assertRaises(ValueError):
            initialize_from_policy(
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                mirror_receipts=self.receipts,
                mirror_jsonl=self.mirror,
                mirror_transport="filesystem_copy",
                notes="must not recreate closed genesis",
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertFalse(self.ledger.exists())

    def test_initializer_rejects_unreviewed_ledger_mirror_transport(self) -> None:
        self._write_policy(project_home_ledger_mirror_transport="dtn_rsync")
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_initializer_rejects_unauthorized_orion_bulk_evidence_root(self) -> None:
        self._write_policy(orion_bulk_evidence_root=str(self.root / "alternate"))
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_initializer_rejects_project_home_bulk_archive_role(self) -> None:
        self._write_policy(project_home_retention_role="permanent_archive")
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_initializer_rejects_missing_orion_only_durability_risk(self) -> None:
        policy = json.loads(self.policy.read_text(encoding="utf-8"))
        policy["long_term_storage"]["risk"] = ""
        self.policy.write_text(json.dumps(policy), encoding="utf-8")
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_mutable_source_tree_cannot_create_manifest(self) -> None:
        with self.assertRaises((FileNotFoundError, ValueError)):
            create_manifest(self.config, authorized_pic_root=self.pic_root)

    def test_manifest_rejects_unauthorized_root(self) -> None:
        self._write_config(pic_root=str(self.root / "alternate"))
        with self.assertRaises(ValueError):
            self._create_manifest()

    def test_policy_promotion_rejects_locked_storage_policy(self) -> None:
        self._write_policy(ledger_genesis_allowed=True, ledger_genesis=None)
        with self.assertRaisesRegex(ValueError, "retained ledger bindings"):
            self._promote_policy()

    def test_reservation_rejects_unreviewed_ledger_mirror_transport(self) -> None:
        self._write_policy(project_home_ledger_mirror_transport="dtn_rsync")
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_reservation_ignores_unpromoted_source_policy_override(self) -> None:
        manifest_path = self._create_manifest()
        self._write_policy(ledger_genesis_allowed=True, ledger_genesis=None)
        reservation = self._reserve(manifest_path)
        self.assertEqual(reservation["state"], "reserved")

    def test_active_policy_snapshot_rejects_configured_project_home_root_alias(
        self,
    ) -> None:
        project_home_alias = self.root / "project_home_alias"
        project_home_alias.symlink_to(self.project_home_root, target_is_directory=True)
        self._write_policy(
            project_home_mirror_root=str(project_home_alias),
            **self._write_storage_preflight_evidence(
                project_home_root=project_home_alias
            ),
        )
        with self.assertRaisesRegex(ValueError, "mirror root is not authorized"):
            promote(
                self.policy,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=project_home_alias,
            )

    def test_reservation_rejects_active_promotion_record_tamper(self) -> None:
        manifest_path = self._create_manifest()
        promotion_path = self.pic_root / "policy" / "active_promotion.json"
        promotion = json.loads(promotion_path.read_text(encoding="utf-8"))
        promotion["policy_sha256"] = "0" * 64
        promotion_path.chmod(0o644)
        promotion_path.write_text(json.dumps(promotion), encoding="utf-8")
        promotion_path.chmod(0o444)
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_active_policy_snapshot_rejects_noninteger_promotion_schema_version(self) -> None:
        promotion_paths = [
            self.pic_root / "policy" / "active_promotion.json",
            self.project_home_root / "policy" / "active_promotion.json",
        ]
        original = json.loads(promotion_paths[0].read_text(encoding="utf-8"))
        for schema_version in (True, 1.0, "1"):
            with self.subTest(schema_version=schema_version):
                payload = json.dumps({**original, "schema_version": schema_version})
                for path in promotion_paths:
                    path.chmod(0o644)
                    path.write_text(payload, encoding="utf-8")
                    path.chmod(0o444)
                with self.assertRaisesRegex(ValueError, "promotion record"):
                    require_storage_policy_unlock_snapshot(
                        control_plane_version=self.control_plane_version,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                    )

    def test_pending_markers_reject_noninteger_schema_versions(self) -> None:
        reservation = {
            "reservation_id": "reservation-1",
            "submission_id": "submission-1",
            "manifest_path": "/retained/manifest.json",
            "manifest_sha256": "1" * 64,
            "control_plane_version": self.control_plane_version,
        }
        current = {
            "schema_version": 2,
            "state": "scheduler_job_id_received",
            "reservation_id": reservation["reservation_id"],
            "submission_id": reservation["submission_id"],
            "manifest_path": reservation["manifest_path"],
            "manifest_sha256": reservation["manifest_sha256"],
            "control_plane_version": reservation["control_plane_version"],
            "job_id": "12345",
        }
        for schema_version in (True, 2.0, "2"):
            with self.subTest(current_schema_version=schema_version):
                with self.assertRaisesRegex(ValueError, "Pending marker"):
                    _require_current_reservation_marker(
                        {**current, "schema_version": schema_version},
                        reservation,
                        control_plane_version=self.control_plane_version,
                    )
        received = {**current, "schema_version": 1}
        received.pop("control_plane_version")
        for schema_version in (True, 1.0, "1"):
            with self.subTest(received_schema_version=schema_version):
                with self.assertRaisesRegex(ValueError, "marker schema"):
                    require_closed_received_marker(
                        {**received, "schema_version": schema_version},
                        reservation,
                        job_id="12345",
                    )

    def test_reserve_has_no_caller_selected_storage_policy(self) -> None:
        self.assertNotIn("storage_policy_path", inspect.signature(reserve).parameters)
        wrapper = Path(__file__).with_name("submit_frontier_job.sh")
        self.assertNotIn("STORAGE_POLICY", wrapper.read_text(encoding="utf-8"))

    def test_reservation_rejects_missing_genesis(self) -> None:
        manifest_path = self._create_manifest()
        self.ledger.write_text("", encoding="utf-8")
        self.receipts.write_text("", encoding="utf-8")
        self.mirror.write_text("", encoding="utf-8")
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_unauthorized_ledgers_reject_before_creating_lock_files(self) -> None:
        manifest_path = self._create_manifest()
        outside = self.root / "outside-ledger"
        for operation in ("reserve", "transition", "reconcile"):
            with self.subTest(operation=operation):
                ledger = outside / operation / "node_hours.jsonl"
                if operation == "reserve":
                    invoke = lambda: reserve(
                        manifest_path=manifest_path,
                        ledger_jsonl=ledger,
                        ledger_csv=self.csv,
                        receipts_jsonl=self.receipts,
                        mirror_jsonl=self.mirror,
                        node_hour_cap=10000.0,
                        control_plane_dir=self.control_plane_dir,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                    )
                elif operation == "transition":
                    invoke = lambda: transition(
                        reservation_id=str(uuid.uuid4()),
                        event_type="reservation_cancelled",
                        state="cancelled",
                        ledger_jsonl=ledger,
                        ledger_csv=self.csv,
                        receipts_jsonl=self.receipts,
                        mirror_jsonl=self.mirror,
                        control_plane_dir=self.control_plane_dir,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                    )
                else:
                    invoke = lambda: reconcile(
                        job_id="12345",
                        ledger_jsonl=ledger,
                        ledger_csv=self.csv,
                        receipts_jsonl=self.receipts,
                        mirror_jsonl=self.mirror,
                        control_plane_dir=self.control_plane_dir,
                        authorized_pic_root=self.pic_root,
                        authorized_project_home_root=self.project_home_root,
                    )
                with self.assertRaises(ValueError):
                    invoke()
                self.assertFalse(
                    ledger.with_suffix(ledger.suffix + ".lock").exists()
                )

    def test_symlink_ledger_aliases_reject_before_outside_lock_or_csv_write(self) -> None:
        manifest_path = self._create_manifest()
        outside = self.root / "outside-ledger-alias"
        outside.mkdir()
        alias = outside / "node_hours.jsonl"
        alias.symlink_to(self.ledger)
        with self.assertRaises(ValueError):
            reserve(
                manifest_path=manifest_path,
                ledger_jsonl=alias,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                node_hour_cap=10000.0,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertFalse(alias.with_suffix(".jsonl.lock").exists())

        csv_alias = outside / "node_hours.csv"
        csv_alias.symlink_to(self.csv)
        with self.assertRaises(ValueError):
            reserve(
                manifest_path=manifest_path,
                ledger_jsonl=self.ledger,
                ledger_csv=csv_alias,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                node_hour_cap=10000.0,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        self.assertTrue(csv_alias.is_symlink())

    def test_symlink_manifest_alias_rejects_before_outside_attachment_write(self) -> None:
        manifest_path = self._create_manifest()
        outside = self.root / "outside-manifest-alias"
        outside.mkdir()
        alias = outside / "pre_submit_manifest.json"
        alias.symlink_to(manifest_path)
        with self.assertRaises(ValueError):
            self._reserve(alias)
        self.assertFalse((outside / "reservation_id.txt").exists())
        self.assertFalse((outside / "manifest_sha256.txt").exists())

    def test_mutable_root_symlink_aliases_reject_before_outside_write(self) -> None:
        install_root = self.root / "install-root"
        install_root.mkdir()
        outside_install = self.root / "outside-install"
        outside_install.mkdir()
        (install_root / "control_plane").symlink_to(
            outside_install, target_is_directory=True
        )
        with self.assertRaises(ValueError):
            install(install_root)
        self.assertEqual(list(outside_install.iterdir()), [])

        outside_manifests = self.root / "outside-manifests"
        outside_manifests.mkdir()
        (self.pic_root / "manifests").symlink_to(
            outside_manifests, target_is_directory=True
        )
        with self.assertRaises(ValueError):
            self._create_manifest()
        self.assertEqual(list(outside_manifests.iterdir()), [])

    def test_clean_candidate_root_symlink_rejects_before_outside_write(self) -> None:
        source_root = self._clean_source("alias-candidate-source")
        executable, profile = self._build_profile(
            source_root,
            self.pic_root / "alias-candidate-build",
            "alias-profile",
        )
        outside = self.root / "outside-clean-candidates"
        outside.mkdir()
        (self.pic_root / "clean_candidates").symlink_to(
            outside, target_is_directory=True
        )
        with self.assertRaises(ValueError):
            create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="alias-profile",
                prepared_artifact_inventory=self._prepared_artifact_inventory(),
                freeze_id=str(uuid.uuid4()),
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
            )
        self.assertEqual(list(outside.iterdir()), [])

    def test_policy_root_symlink_rejects_before_outside_write(self) -> None:
        policy_root = self.pic_root / "policy"
        shutil.rmtree(policy_root)
        outside = self.root / "outside-policy"
        outside.mkdir()
        policy_root.symlink_to(outside, target_is_directory=True)
        with self.assertRaises((NotADirectoryError, ValueError)):
            self._promote_policy()
        self.assertEqual(list(outside.iterdir()), [])

    def test_installed_version_symlink_alias_is_not_a_control_plane(self) -> None:
        outside = self.root / "outside-control-plane-version"
        outside.mkdir()
        alias = self.pic_root / "control_plane" / "alias"
        alias.symlink_to(self.control_plane_dir, target_is_directory=True)
        from control_plane_common import verify_installed_control_plane

        with self.assertRaises(ValueError):
            verify_installed_control_plane(alias, authorized_pic_root=self.pic_root)
        self.assertEqual(list(outside.iterdir()), [])

    def test_manifest_creator_rejects_installed_version_symlink_alias(self) -> None:
        alias = self.pic_root / "control_plane" / "alias"
        alias.symlink_to(self.control_plane_dir, target_is_directory=True)
        with self.assertRaises(ValueError):
            create_manifest(
                self.config,
                control_plane_dir=alias,
                authorized_pic_root=self.pic_root,
            )

    def test_clean_candidate_creator_rejects_installed_version_symlink_alias(
        self,
    ) -> None:
        alias = self.pic_root / "control_plane" / "alias"
        alias.symlink_to(self.control_plane_dir, target_is_directory=True)
        build = self.pic_root / "alias-control-plane-build"
        build.mkdir()
        executable = build / "athena"
        executable.write_text("placeholder\n", encoding="utf-8")
        profile = build / "build_profile.json"
        profile.write_text("{}\n", encoding="utf-8")
        with patch(
            "create_clean_candidate_freeze._git",
            side_effect=AssertionError("alias reached source inspection"),
        ) as git:
            with self.assertRaises(ValueError):
                create_freeze(
                    source_root=self.sources,
                    executable=executable,
                    build_profile=profile,
                    build_profile_id="must-reject-alias",
                    prepared_artifact_inventory=self._prepared_artifact_inventory(),
                    control_plane_dir=alias,
                    authorized_pic_root=self.pic_root,
                )
        git.assert_not_called()

    def test_orion_build_profile_writer_rejects_installed_version_symlink_alias(
        self,
    ) -> None:
        alias = self.pic_root / "control_plane" / "alias"
        alias.symlink_to(self.control_plane_dir, target_is_directory=True)
        build = self.pic_root / "alias-profile-writer-build"
        build.mkdir()
        executable = build / "athena"
        executable.write_text("placeholder\n", encoding="utf-8")
        toolchain = build / "toolchain.txt"
        toolchain.write_text("Frontier test toolchain\n", encoding="utf-8")
        build_invocations = build / "build-invocations.json"
        build_invocations.write_text(
            '{"build":["/fake/cmake"],"configure":["/fake/cmake"]}\n',
            encoding="utf-8",
        )
        with patch(
            "write_orion_build_profile._source_identity",
            side_effect=AssertionError("alias reached source inspection"),
        ) as source_identity:
            with self.assertRaises(ValueError):
                write_profile(
                    source_root=self.sources,
                    fresh_source_root=self.sources,
                    executable=executable,
                    output=build / "build_profile.json",
                    profile_id="must-reject-alias",
                    expected_git_commit="0" * 40,
                    configure_log=toolchain,
                    build_log=toolchain,
                    cmake_cache=toolchain,
                    module_list=toolchain,
                    toolchain_file=toolchain,
                    build_invocations_file=build_invocations,
                    git_status_preconfigure_file=toolchain,
                    git_status_file=toolchain,
                    submodule_status_file=toolchain,
                    environment_allowlist_file=toolchain,
                    build_environment_file=toolchain,
                    control_plane_dir=alias,
                    authorized_pic_root=self.pic_root,
                )
        source_identity.assert_not_called()

    def test_policy_promoter_rejects_installed_version_symlink_alias(self) -> None:
        alias = self.pic_root / "control_plane" / "alias"
        alias.symlink_to(self.control_plane_dir, target_is_directory=True)
        with self.assertRaises(ValueError):
            promote(
                self.policy,
                control_plane_dir=alias,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_submit_wrapper_rejects_version_alias_before_adjacent_python(self) -> None:
        wrapper_root = self.root / "wrapper-pic"
        project_root = self.root / "wrapper-project-home"
        real = wrapper_root / "control_plane" / "real-version"
        real.mkdir(parents=True)
        marker = self.root / "adjacent-validator-executed"
        wrapper_text = Path(__file__).with_name("submit_frontier_job.sh").read_text(
            encoding="utf-8"
        )
        wrapper_text = wrapper_text.replace(
            "/lustre/orion/ast207/proj-shared/dfielding/PIC", str(wrapper_root)
        ).replace(
            "/ccs/proj/ast207/proj-shared/PIC", str(project_root)
        )
        wrapper = real / "submit_frontier_job.sh"
        wrapper.write_text(wrapper_text, encoding="utf-8")
        wrapper.chmod(0o755)
        validator = real / "validate_and_reserve_frontier_job.py"
        validator.write_text(
            "from pathlib import Path\n"
            f"Path({str(marker)!r}).write_text('executed', encoding='utf-8')\n",
            encoding="utf-8",
        )
        alias = wrapper_root / "control_plane" / "alias"
        alias.symlink_to(real, target_is_directory=True)
        result = subprocess.run(
            [
                str(alias / "submit_frontier_job.sh"),
                "unused-manifest",
                "unused-attestation",
            ],
            text=True,
            capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Refusing a symlink alias", result.stderr)
        self.assertFalse(marker.exists())

    def test_trampoline_rejects_symlink_primary_ledger_alias(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        job_script = record_for_role(manifest, "job-script")
        executable = record_for_role(manifest, "executable")
        outside = self.root / "outside-trampoline-ledger"
        outside.mkdir()
        alias = outside / "node_hours.jsonl"
        alias.symlink_to(self.ledger)
        with patch.dict(
            os.environ,
            {
                "PIC_MANIFEST_SHA256": str(reservation["manifest_sha256"]),
                "PIC_RESERVATION_ID": str(reservation["reservation_id"]),
                "PIC_SUBMISSION_ID": self.submission_id,
            },
            clear=True,
        ):
            with self.assertRaises(ValueError):
                launch(
                    manifest_path=manifest_path,
                    manifest_sha256=str(reservation["manifest_sha256"]),
                    job_script_sha256=str(job_script["sha256"]),
                    executable_sha256=str(executable["sha256"]),
                    reservation_id=str(reservation["reservation_id"]),
                    submission_id=self.submission_id,
                    ledger_jsonl=alias,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    runner=lambda *_args, **_kwargs: None,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        self.assertFalse(alias.with_suffix(".jsonl.lock").exists())

    def test_queue_snapshot_digest_is_bound_to_staged_copy(self) -> None:
        from create_pre_submit_manifest import snapshot_file as real_snapshot_file

        def mutate_source_after_snapshot(*args: object, **kwargs: object) -> object:
            record = real_snapshot_file(*args, **kwargs)
            if kwargs.get("role") == "queue-snapshot":
                self.sources.joinpath("queue.txt").write_text(
                    "mutated after snapshot\n", encoding="utf-8"
                )
            return record

        with patch(
            "create_pre_submit_manifest.snapshot_file",
            side_effect=mutate_source_after_snapshot,
        ):
            manifest_path = self._create_manifest()
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        queue_record = record_for_role(manifest, "queue-snapshot")
        self.assertEqual(manifest["queue_snapshot_sha256"], queue_record["sha256"])
        self.assertNotEqual(
            manifest["queue_snapshot_sha256"], sha256(self.sources / "queue.txt")
        )
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_reconcile_rejects_unknown_scheduler_terminal_state(self) -> None:
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        self._attach(str(reservation["reservation_id"]))
        with patch(
            "reconcile_frontier_job._scheduler_result",
            return_value=("NOT_A_SLURM_TERMINAL_STATE", 300, 1),
        ):
            with self.assertRaises(ValueError):
                reconcile(
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )

    def test_external_provenance_commands_are_absolute_and_path_independent(self) -> None:
        import create_clean_candidate_freeze
        import reconcile_frontier_job
        import validate_and_reserve_frontier_job

        self.assertEqual(TRUSTED_GIT, "/usr/bin/git")
        self.assertEqual(TRUSTED_SQUEUE, "/usr/bin/squeue")
        self.assertEqual(TRUSTED_SCONTROL, "/usr/bin/scontrol")
        self.assertEqual(TRUSTED_SACCT, "/usr/bin/sacct")
        self.assertEqual(TRUSTED_SBATCH, "/usr/bin/sbatch")
        self.assertEqual(TRUSTED_SCANCEL, "/usr/bin/scancel")
        self.assertEqual(TRUSTED_PYTHON, "/opt/cray/pe/python/3.11.7/bin/python3")
        with patch.dict(
            os.environ,
            {"LOGNAME": "forged-empty-queue-user", "USER": "forged-empty-queue-user"},
        ):
            with patch.object(
                validate_and_reserve_frontier_job.subprocess,
                "check_output",
                return_value="",
            ) as queue:
                validate_and_reserve_frontier_job._queue_output()
        self.assertEqual(queue.call_args.args[0][0], TRUSTED_SQUEUE)
        self.assertEqual(queue.call_args.args[0][2], pwd.getpwuid(os.getuid()).pw_name)
        self.assertEqual(queue.call_args.kwargs["env"], trusted_slurm_environment())
        with patch.object(
            validate_and_reserve_frontier_job.subprocess,
            "check_output",
            return_value="",
        ) as scheduler:
            validate_and_reserve_frontier_job._scheduler_job_output("123")
        self.assertEqual(scheduler.call_args.args[0][0], TRUSTED_SCONTROL)
        self.assertEqual(scheduler.call_args.kwargs["env"], trusted_slurm_environment())
        with patch.object(
            reconcile_frontier_job.subprocess,
            "check_output",
            return_value="",
        ) as accounting_call:
            with self.assertRaises(ValueError):
                reconcile_frontier_job._scheduler_result("123", "reservation")
        self.assertEqual(accounting_call.call_args.args[0][0], TRUSTED_SACCT)
        self.assertEqual(accounting_call.call_args.kwargs["env"], trusted_slurm_environment())
        self.assertIn(
            'trusted_git_command("-C", str(source_root), *arguments)',
            inspect.getsource(create_clean_candidate_freeze._git),
        )
        self.assertEqual(
            trusted_git_command("status"),
            [
                "/usr/bin/git",
                "-c",
                "core.fsmonitor=false",
                "-c",
                "core.hooksPath=/dev/null",
                "status",
            ],
        )
        wrapper = Path(__file__).with_name("submit_frontier_job.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn('PYTHON=(/opt/cray/pe/python/3.11.7/bin/python3 -I)', wrapper)
        self.assertIn('CONTROL_PLANE=("${PYTHON[@]}" "$RUNNER")', wrapper)
        self.assertIn(
            'SLURM_ENV=(/usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin '
            'SLURM_CLUSTERS=frontier)',
            wrapper,
        )
        self.assertIn(
            '"${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py verify-control-plane',
            wrapper,
        )
        self.assertIn('SBATCH="/usr/bin/sbatch"', wrapper)
        self.assertIn('SCANCEL="/usr/bin/scancel"', wrapper)
        self.assertIn('SCONTROL="/usr/bin/scontrol"', wrapper)
        self.assertIn('"${SLURM_ENV[@]}" "$SBATCH" --parsable --hold', wrapper)
        self.assertIn("--export=NIL", wrapper)
        self.assertNotIn("--get-user-env", wrapper)
        self.assertIn('"${SLURM_ENV[@]}" "$SCANCEL" "$job_id"', wrapper)
        self.assertIn('"${SLURM_ENV[@]}" "$SCONTROL" release "$job_id"', wrapper)
        launch_wrapper = Path(__file__).with_name("launch_with_frontier_profile.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn("/opt/cray/pe/python/3.11.7/bin/python3 -E -s -", launch_wrapper)

    def test_installed_python_entrypoints_are_isolated_from_pythonpath(self) -> None:
        entrypoints = [
            path
            for path in self.control_plane_dir.glob("*.py")
            if path.name
            not in {
                "q011_pressure_review_packet_verifier.py",
                "test_control_plane.py",
            }
        ]
        self.assertTrue(entrypoints)
        for entrypoint in entrypoints:
            self.assertEqual(
                entrypoint.read_text(encoding="utf-8").splitlines()[0],
                "#!/opt/cray/pe/python/3.11.7/bin/python3 -I",
            )
        poisoned = self.root / "poisoned-pythonpath"
        poisoned.mkdir()
        marker = self.root / "pythonpath-imported"
        (poisoned / "pathlib.py").write_text(
            f"open({str(marker)!r}, 'w').write('imported')\n",
            encoding="utf-8",
        )
        result = subprocess.run(
            [
                str(self.control_plane_dir / "run_control_plane.py"),
                "create_pre_submit_manifest.py",
                "--help",
            ],
            capture_output=True,
            text=True,
            env={**os.environ, "PYTHONPATH": str(poisoned)},
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(marker.exists())
        direct = subprocess.run(
            [str(self.control_plane_dir / "create_pre_submit_manifest.py"), "--help"],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(direct.returncode, 0)
        self.assertIn("run_control_plane.py", direct.stderr)
        self.control_plane_dir.chmod(0o755)
        (self.control_plane_dir / "pathlib.py").write_text(
            f"open({str(marker)!r}, 'w').write('imported')\n",
            encoding="utf-8",
        )
        self.control_plane_dir.chmod(0o555)
        adjacent = subprocess.run(
            [
                str(self.control_plane_dir / "run_control_plane.py"),
                "create_pre_submit_manifest.py",
                "--help",
            ],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(adjacent.returncode, 0)
        self.assertIn("entries differ", adjacent.stderr)
        self.assertFalse(marker.exists())

    def test_clean_candidate_schema_closes_production_profile_id(self) -> None:
        schema = json.loads(
            (self.control_plane_dir / "clean_candidate.schema.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            schema["properties"]["build"]["properties"]["profile_id"],
            {"const": "hip-mpi-release-paper-pic"},
        )

    def test_production_profile_explicitly_pins_double_precision(self) -> None:
        from control_plane_common import AUTHORIZED_PIC_ROOT
        from control_plane_common import PRODUCTION_BUILD_PROFILE
        from control_plane_common import production_build_invocations

        configure = production_build_invocations(
            authorized_pic_root=AUTHORIZED_PIC_ROOT,
            git_commit="0" * 40,
            profile_id=PRODUCTION_BUILD_PROFILE,
        )["configure"]
        self.assertEqual(configure.count("-DAthena_SINGLE_PRECISION=OFF"), 1)

    def test_clean_candidate_schema_requires_prepared_artifact_inventories(self) -> None:
        schema = json.loads(
            (self.control_plane_dir / "clean_candidate.schema.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(schema["properties"]["schema_version"], {"const": 4})
        self.assertIn("prepared_artifacts", schema["required"])
        prepared = schema["properties"]["prepared_artifacts"]
        self.assertEqual(
            set(prepared["required"]),
            {"inventory_path", "inventory_sha256", "paper_decks", "analyzers"},
        )
        self.assertEqual(prepared["properties"]["paper_decks"]["minItems"], 1)
        self.assertEqual(prepared["properties"]["analyzers"]["minItems"], 1)

    def test_production_semantics_reject_self_authored_profile_receipt_chain(self) -> None:
        from control_plane_common import AUTHORIZED_PIC_ROOT
        from control_plane_common import require_production_build_provenance

        with self.assertRaisesRegex(ValueError, "Unsupported installed Frontier build profile"):
            require_production_build_provenance(
                authorized_pic_root=AUTHORIZED_PIC_ROOT,
                git_commit="0" * 40,
                profile_id="self-authored-profile",
                toolchain="self-authored toolchain",
                invocations={"configure": ["/tmp/cmake"], "build": ["/tmp/cmake"]},
                module_list=b"self-authored/module\n",
                environment_allowlist=b"self-authored=1\n",
                build_environment=b'{"HOME": "/tmp"}\n',
            )

    def test_production_semantics_reject_extra_module_provenance(self) -> None:
        from control_plane_common import AUTHORIZED_PIC_ROOT
        from control_plane_common import PRODUCTION_BUILD_ENVIRONMENT
        from control_plane_common import PRODUCTION_BUILD_PROFILE
        from control_plane_common import PRODUCTION_TOOLCHAIN_DESCRIPTION
        from control_plane_common import production_build_invocations
        from control_plane_common import production_environment_allowlist_bytes
        from control_plane_common import production_module_list_bytes
        from control_plane_common import require_production_build_provenance

        commit = "0" * 40
        environment = (
            json.dumps(PRODUCTION_BUILD_ENVIRONMENT, indent=2, sort_keys=True) + "\n"
        ).encode("utf-8")
        with self.assertRaisesRegex(ValueError, "exact selection"):
            require_production_build_provenance(
                authorized_pic_root=AUTHORIZED_PIC_ROOT,
                git_commit=commit,
                profile_id=PRODUCTION_BUILD_PROFILE,
                toolchain=PRODUCTION_TOOLCHAIN_DESCRIPTION,
                invocations=production_build_invocations(
                    authorized_pic_root=AUTHORIZED_PIC_ROOT,
                    git_commit=commit,
                    profile_id=PRODUCTION_BUILD_PROFILE,
                ),
                module_list=production_module_list_bytes() + b"caller/module\n",
                environment_allowlist=production_environment_allowlist_bytes(),
                build_environment=environment,
            )

    def test_reserve_cli_rejects_caller_forged_queue_output_file(self) -> None:
        import validate_and_reserve_frontier_job

        arguments = [
            "validate_and_reserve_frontier_job.py",
            "reserve",
            "--manifest",
            "manifest.json",
            "--ledger-jsonl",
            "node_hours.jsonl",
            "--ledger-csv",
            "node_hours.csv",
            "--receipts-jsonl",
            "mirror_receipts.jsonl",
            "--mirror-jsonl",
            "mirror.jsonl",
            "--node-hour-cap",
            "10000",
            "--queue-output-file",
            "forged-queue.txt",
        ]
        with patch.object(sys, "argv", arguments), patch.object(
            validate_and_reserve_frontier_job, "reserve"
        ) as reserve_call:
            with self.assertRaises(SystemExit):
                validate_and_reserve_frontier_job.main()
        reserve_call.assert_not_called()

    def test_reserve_api_rejects_caller_forged_queue_output_file(self) -> None:
        self.assertNotIn(
            "_queue_output_file_for_test", inspect.signature(reserve).parameters
        )
        with self.assertRaises(TypeError):
            reserve(
                manifest_path=Path("manifest.json"),
                ledger_jsonl=Path("node_hours.jsonl"),
                ledger_csv=Path("node_hours.csv"),
                receipts_jsonl=Path("mirror_receipts.jsonl"),
                mirror_jsonl=Path("mirror.jsonl"),
                node_hour_cap=10000.0,
                _queue_output_file_for_test=Path("forged-queue.txt"),
            )

    def test_reservation_rejects_installed_checksum_drift(self) -> None:
        manifest_path = self._create_manifest()
        path = self.control_plane_dir / "frontier_pic_environment.sh"
        path.chmod(0o755)
        with path.open("a", encoding="utf-8") as stream:
            stream.write("\n# mutation\n")
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_frontier_environment_profiles_are_closed_and_redacted(self) -> None:
        path = Path(__file__).with_name("frontier_pic_environment.sh")
        command = 'module() { :; }\nsource "$1" || exit $?\nrecord_pic_environment\n'
        scrubbed = {
            "HSA_XNACK",
            "MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED",
            "MPICH_OFI_NIC_POLICY",
            "MPICH_GPU_IPC_CACHE_MAX_SIZE",
            "MPICH_MPIIO_HINTS",
            "MPICH_OFI_NUM_CQ_ENTRIES",
            "FI_MR_CACHE_MONITOR",
            "FI_CXI_RX_MATCH_MODE",
        }

        def capture(profile: str) -> dict[str, str]:
            environment = dict(os.environ)
            environment["PIC_FRONTIER_PROFILE"] = profile
            environment["MODULEPATH"] = "/tmp/caller-controlled-modulepath"
            for name in scrubbed:
                environment.pop(name, None)
            result = subprocess.run(
                ["/bin/bash", "-c", command, "bash", str(path)],
                check=True,
                capture_output=True,
                text=True,
                env=environment,
            )
            return dict(line.split("=", 1) for line in result.stdout.splitlines())

        baseline = capture("frontier_minimum_supported")
        self.assertEqual(baseline["PIC_FRONTIER_PROFILE"], "frontier_minimum_supported")
        self.assertEqual(baseline["HSA_XNACK"], "0")
        self.assertEqual(baseline["MPICH_GPU_SUPPORT_ENABLED"], "1")
        self.assertEqual(baseline["MPICH_OFI_NIC_POLICY"], "<unset>")
        self.assertEqual(baseline["SLURM_EXPORT_ENV"], "ALL")
        self.assertEqual(baseline["MODULEPATH"], PRODUCTION_RUNTIME_MODULEPATH)

        xnack = capture("frontier_xnack1_experimental")
        self.assertEqual(xnack["HSA_XNACK"], "1")
        self.assertEqual(xnack["MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED"], "1")
        self.assertEqual(xnack["MODULEPATH"], PRODUCTION_RUNTIME_MODULEPATH)

        ofi = capture("frontier_ofi_tuned_experimental")
        self.assertEqual(ofi["MPICH_OFI_NIC_POLICY"], "GPU")
        self.assertEqual(ofi["MPICH_GPU_IPC_CACHE_MAX_SIZE"], "1000")
        self.assertEqual(ofi["FI_CXI_RX_MATCH_MODE"], "software")
        self.assertEqual(ofi["MODULEPATH"], PRODUCTION_RUNTIME_MODULEPATH)

        result = subprocess.run(
            ["/bin/bash", "-c", command, "bash", str(path)],
            capture_output=True,
            text=True,
            env={**os.environ, "PIC_FRONTIER_PROFILE": "unsupported"},
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Unsupported PIC_FRONTIER_PROFILE", result.stderr)

        mutation_marker = self.root / "unsupported-profile-module-called"
        unsupported = subprocess.run(
            [
                "/bin/bash",
                "-c",
                'module() { : > "$MUTATION_MARKER"; }\nsource "$1"',
                "bash",
                str(path),
            ],
            capture_output=True,
            text=True,
            env={
                **os.environ,
                "PIC_FRONTIER_PROFILE": "unsupported",
                "MUTATION_MARKER": str(mutation_marker),
            },
        )
        self.assertNotEqual(unsupported.returncode, 0)
        self.assertFalse(mutation_marker.exists())

        failure_marker = self.root / "module-failure-called"
        failed_module = subprocess.run(
            [
                "/bin/bash",
                "-c",
                'module() { printf "%s\\n" "$*" >> "$FAILURE_MARKER"; return 1; }\nsource "$1"',
                "bash",
                str(path),
            ],
            capture_output=True,
            text=True,
            env={
                **os.environ,
                "PIC_FRONTIER_PROFILE": "frontier_minimum_supported",
                "FAILURE_MARKER": str(failure_marker),
            },
        )
        self.assertNotEqual(failed_module.returncode, 0)
        self.assertEqual(
            failure_marker.read_text(encoding="utf-8").splitlines(),
            ["--force purge"],
        )

    def test_reservation_rejects_wrong_account(self) -> None:
        self._write(
            "job.sh",
            "#!/bin/bash\n#SBATCH -A OTHER\n#SBATCH -p batch\n#SBATCH -q debug\n"
            f"#SBATCH -o {self.pic_root}/logs/slurm/%x.%j.log\n"
            "#SBATCH -N 1\n#SBATCH -t 00:10:00\n",
        )
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_reservation_rejects_stale_timeout_margin(self) -> None:
        self._write_timeout(expires_delta=timedelta(minutes=-1))
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_manifest_rejects_nonboolean_short_job_registration(self) -> None:
        self._write_config(registered_short_nonproduction="false")
        with self.assertRaises(ValueError):
            self._create_manifest()

    def test_manifest_rejects_nonstring_required_config_value(self) -> None:
        self._write_config(git_commit=True)
        with self.assertRaises(ValueError):
            self._create_manifest()

    def test_reservation_rejects_timeout_margin_numeric_aliases(self) -> None:
        for field, value in [
            ("athena_walltime_seconds", True),
            ("athena_walltime_seconds", "300"),
            ("athena_walltime_seconds", 300.9),
            ("scheduler_walltime_seconds", 600.0),
        ]:
            with self.subTest(field=field, value=value):
                self._write_timeout()
                margin_path = self.sources / "timeout.json"
                margin = json.loads(margin_path.read_text(encoding="utf-8"))
                margin[field] = value
                margin_path.write_text(json.dumps(margin), encoding="utf-8")
                submission_id = str(uuid.uuid4())
                self._write_config(
                    submission_id=submission_id,
                    artifact_dir=str(
                        self.pic_root / "runs" / "f0_hipmpi_smoke" / submission_id
                    ),
                )
                with self.assertRaises(ValueError):
                    self._reserve(self._create_manifest())

    def test_reservation_rejects_stale_site_policy_timestamp(self) -> None:
        stale = datetime.now(timezone.utc) - timedelta(days=2)
        self._write_config(site_policy_checked_utc=self._utc(stale))
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_concurrent_reservations_are_serialized(self) -> None:
        first = self._create_manifest()
        second_submission_id = str(uuid.uuid4())
        self._write_config(
            submission_id=second_submission_id,
            artifact_dir=str(
                self.pic_root / "runs" / "f0_hipmpi_smoke" / second_submission_id
            ),
        )
        second = self._create_manifest()
        identifiers = [str(uuid.uuid4()), str(uuid.uuid4())]

        def invoke(item: tuple[Path, str]) -> object:
            try:
                return self._reserve(item[0], reservation_id=item[1])
            except ValueError as error:
                return error

        with ThreadPoolExecutor(max_workers=2) as executor:
            results = list(executor.map(invoke, zip([first, second], identifiers)))
        self.assertEqual(sum(isinstance(result, dict) for result in results), 1)
        self.assertEqual(sum(isinstance(result, ValueError) for result in results), 1)

    def test_manual_direct_srun_accounting_appends_reviewed_events_once(self) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = {
            "authorization": authorization,
            "ledger_jsonl": self.ledger,
            "ledger_csv": self.csv,
            "receipts_jsonl": self.receipts,
            "mirror_jsonl": self.mirror,
            "control_plane_dir": self.control_plane_dir,
            "authorized_pic_root": self.pic_root,
            "authorized_project_home_root": self.project_home_root,
        }
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            first = reconcile_manual_allocations(**arguments)
            second = reconcile_manual_allocations(**arguments)
        self.assertEqual(first, second)
        records = validate_primary_chain(self.ledger)
        self.assertEqual(len(records), 3)
        self.assertEqual(records[-1]["event_type"], "manual_allocation_reconciliation")
        self.assertEqual(records[-1]["accounting_scope"], "manual_direct_srun_accounting_only")
        self.assertIs(records[-1]["scientific_evidence_eligible"], False)
        self.assertRegex(str(records[-1]["active_policy_sha256"]), r"^[0-9a-f]{64}$")
        self.assertRegex(str(records[-1]["active_promotion_sha256"]), r"^[0-9a-f]{64}$")
        self.assertAlmostEqual(
            accounting(records)["cumulative_consumed_node_hours"], 12.0 / 3600.0
        )
        self.assertEqual(
            json.loads(self.mirror.read_text(encoding="utf-8").splitlines()[-1]),
            records[-1],
        )
        self.assertIn("manual_accounting_authorization_sha256", self.csv.read_text())

    def test_manual_direct_scheduler_accounting_appends_nonqualifying_event(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization(
            accounting_scope="manual_direct_scheduler_accounting_only",
        )
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            reconcile_manual_allocations(
                **self._manual_accounting_arguments(authorization)
            )
        records = validate_primary_chain(self.ledger)
        self.assertEqual(
            records[-1]["accounting_scope"],
            "manual_direct_scheduler_accounting_only",
        )
        self.assertEqual(
            records[-1]["notes"],
            "Reviewed direct-scheduler accounting only; "
            "ineligible for scientific evidence.",
        )
        self.assertIs(records[-1]["scientific_evidence_eligible"], False)
        self.assertRegex(str(records[-1]["active_policy_sha256"]), r"^[0-9a-f]{64}$")
        self.assertRegex(str(records[-1]["active_promotion_sha256"]), r"^[0-9a-f]{64}$")
        self.assertAlmostEqual(
            accounting(records)["cumulative_consumed_node_hours"], 12.0 / 3600.0
        )
        self.assertEqual(
            json.loads(self.mirror.read_text(encoding="utf-8").splitlines()[-1]),
            records[-1],
        )
        self.assertIn("manual_accounting_authorization_sha256", self.csv.read_text())

    def test_manual_direct_srun_accounting_retry_appends_only_missing_suffix(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = {
            "authorization": authorization,
            "ledger_jsonl": self.ledger,
            "ledger_csv": self.csv,
            "receipts_jsonl": self.receipts,
            "mirror_jsonl": self.mirror,
            "control_plane_dir": self.control_plane_dir,
            "authorized_pic_root": self.pic_root,
            "authorized_project_home_root": self.project_home_root,
        }
        real_append = reconcile_manual_frontier_allocations.append_primary_event_locked
        append_count = 0

        def fail_second_append(*args: object, **kwargs: object) -> dict[str, object]:
            nonlocal append_count
            append_count += 1
            if append_count == 2:
                raise RuntimeError("interrupted")
            return real_append(*args, **kwargs)

        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ), patch(
            "reconcile_manual_frontier_allocations.append_primary_event_locked",
            side_effect=fail_second_append,
        ):
            with self.assertRaisesRegex(RuntimeError, "interrupted"):
                reconcile_manual_allocations(**arguments)
        prefix = self.ledger.read_bytes()
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            events = reconcile_manual_allocations(**arguments)
        self.assertEqual(len(events), 2)
        self.assertTrue(self.ledger.read_bytes().startswith(prefix))
        self.assertEqual(
            [record["job_id"] for record in validate_primary_chain(self.ledger)[1:]],
            ["4746332", "4746335"],
        )

    def test_manual_direct_srun_accounting_full_retry_uses_historical_prefix(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = {
            "authorization": authorization,
            "ledger_jsonl": self.ledger,
            "ledger_csv": self.csv,
            "receipts_jsonl": self.receipts,
            "mirror_jsonl": self.mirror,
            "control_plane_dir": self.control_plane_dir,
            "authorized_pic_root": self.pic_root,
            "authorized_project_home_root": self.project_home_root,
        }
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            first = reconcile_manual_allocations(**arguments)
        second_authorization = self._write_manual_accounting_authorization(
            authorization_id="later-reviewed-authorization",
            jobs=[{"job_id": "4746999", "expected_qos": "debug"}],
            bind_policy=False,
        )
        self._write_policy(
            manual_accounting_authorizations=[
                self._manual_accounting_policy_binding(authorization),
                self._manual_accounting_policy_binding(second_authorization),
            ]
        )
        self._promote_policy()
        self._append_registered_probe()
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            second = reconcile_manual_allocations(**arguments)
        self.assertEqual(first, second)
        self.assertAlmostEqual(
            accounting(validate_primary_chain(self.ledger))[
                "cumulative_consumed_node_hours"
            ],
            12.0 / 3600.0,
        )

    def test_manual_direct_srun_accounting_partial_retry_after_retained_successor_rotation(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = {
            "authorization": authorization,
            "ledger_jsonl": self.ledger,
            "ledger_csv": self.csv,
            "receipts_jsonl": self.receipts,
            "mirror_jsonl": self.mirror,
            "control_plane_dir": self.control_plane_dir,
            "authorized_pic_root": self.pic_root,
            "authorized_project_home_root": self.project_home_root,
        }
        real_append = reconcile_manual_frontier_allocations.append_primary_event_locked
        append_count = 0

        def fail_second_append(*args: object, **kwargs: object) -> dict[str, object]:
            nonlocal append_count
            append_count += 1
            if append_count == 2:
                raise RuntimeError("interrupted")
            return real_append(*args, **kwargs)

        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ), patch(
            "reconcile_manual_frontier_allocations.append_primary_event_locked",
            side_effect=fail_second_append,
        ):
            with self.assertRaisesRegex(RuntimeError, "interrupted"):
                reconcile_manual_allocations(**arguments)
        # Model one predecessor partial prefix created before paired markers existed.
        local_marker, mirror_marker = incomplete_manual_accounting_marker_paths(
            self.ledger, self.mirror
        )
        local_marker.unlink()
        assert mirror_marker is not None
        mirror_marker.unlink()
        second_authorization = self._write_manual_accounting_authorization(
            authorization_id="later-reviewed-authorization",
            jobs=[{"job_id": "4746999", "expected_qos": "debug"}],
            bind_policy=False,
        )
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        self._write_policy(
            manual_accounting_authorizations=[
                self._manual_accounting_policy_binding(authorization),
                self._manual_accounting_policy_binding(second_authorization),
            ],
            installed_control_plane_version=successor.name,
            staged_control_plane_candidate_version=successor.name,
        )
        promote(
            self.policy,
            control_plane_dir=successor,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        arguments["control_plane_dir"] = successor
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            events = reconcile_manual_allocations(**arguments)
        self.assertEqual([event["job_id"] for event in events], ["4746332", "4746335"])
        self.assertFalse(local_marker.exists())
        self.assertFalse(mirror_marker.exists())

    def test_manual_direct_srun_accounting_partial_retry_blocks_other_writers(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = {
            "authorization": authorization,
            "ledger_jsonl": self.ledger,
            "ledger_csv": self.csv,
            "receipts_jsonl": self.receipts,
            "mirror_jsonl": self.mirror,
            "control_plane_dir": self.control_plane_dir,
            "authorized_pic_root": self.pic_root,
            "authorized_project_home_root": self.project_home_root,
        }
        real_append = reconcile_manual_frontier_allocations.append_primary_event_locked
        append_count = 0

        def fail_second_append(*args: object, **kwargs: object) -> dict[str, object]:
            nonlocal append_count
            append_count += 1
            if append_count == 2:
                raise RuntimeError("interrupted")
            return real_append(*args, **kwargs)

        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ), patch(
            "reconcile_manual_frontier_allocations.append_primary_event_locked",
            side_effect=fail_second_append,
        ):
            with self.assertRaisesRegex(RuntimeError, "interrupted"):
                reconcile_manual_allocations(**arguments)
        local_marker, mirror_marker = incomplete_manual_accounting_marker_paths(
            self.ledger, self.mirror
        )
        self.assertTrue(local_marker.is_file())
        assert mirror_marker is not None
        self.assertEqual(local_marker.read_bytes(), mirror_marker.read_bytes())
        marker = json.loads(local_marker.read_text(encoding="utf-8"))
        self.assertEqual(marker["schema_version"], 3)
        self.assertEqual(marker["pre_tranche_sequence_number"], 1)
        self.assertEqual(marker["pre_tranche_authorized_job_count"], 0)
        self.assertEqual(
            marker["pre_tranche_chain_head"],
            validate_primary_chain(self.ledger)[0]["event_sha256"],
        )
        with self.assertRaisesRegex(ValueError, "incomplete manual accounting"):
            append_primary_event(
                self.ledger,
                self.csv,
                self.receipts,
                self.mirror,
                {"event_type": "historical_probe"},
                mirror_transport="filesystem_copy",
            )
        with self.assertRaisesRegex(ValueError, "incomplete manual accounting"):
            self._reserve(self._create_manifest())
        with self.assertRaisesRegex(ValueError, "incomplete manual accounting"):
            self._promote_policy()
        with self.assertRaisesRegex(ValueError, "incomplete manual accounting"):
            repair_ledger_mirror(
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            reconcile_manual_allocations(**arguments)
        self.assertFalse(local_marker.exists())
        self.assertFalse(mirror_marker.exists())

    def test_manual_direct_srun_accounting_recovers_after_orion_append(self) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        self._strand_manual_accounting_after_orion_append(arguments)
        self.assertEqual(len(validate_primary_chain(self.ledger)), 2)
        self.assertEqual(len(validate_primary_chain(self.mirror)), 1)
        events = self._retry_manual_accounting_and_require_clean_markers(arguments)
        self.assertEqual([event["job_id"] for event in events], ["4746332", "4746335"])

    def test_manual_direct_scheduler_accounting_recovers_after_orion_append(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization(
            accounting_scope="manual_direct_scheduler_accounting_only",
        )
        arguments = self._manual_accounting_arguments(authorization)
        self._strand_manual_accounting_after_orion_append(arguments)
        events = self._retry_manual_accounting_and_require_clean_markers(arguments)
        self.assertEqual([event["job_id"] for event in events], ["4746332", "4746335"])
        self.assertEqual(
            events[-1]["accounting_scope"],
            "manual_direct_scheduler_accounting_only",
        )

    def test_manual_accounting_recovery_rejects_retained_attestation_tamper_before_write(
        self,
    ) -> None:
        self._append_registered_probe(include_operator_attestations=True)
        reservation = validate_primary_chain(self.ledger)[1]
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        self._strand_manual_accounting_after_orion_append(arguments)
        attestation = Path(str(reservation["pre_manifest_attestation_path"]))
        member = attestation.parent / "same_account_process_snapshot.txt"
        attestation.parent.chmod(0o700)
        member.chmod(0o600)
        member.write_bytes(member.read_bytes() + b"forged\n")
        member.chmod(0o400)
        attestation.parent.chmod(0o500)
        mirror_before = self.mirror.read_bytes()
        receipts_before = self.receipts.read_bytes()
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ), self.assertRaisesRegex(ValueError, "checksum differs"):
            reconcile_manual_allocations(**arguments)
        self.assertEqual(self.mirror.read_bytes(), mirror_before)
        self.assertEqual(self.receipts.read_bytes(), receipts_before)

    def test_manual_direct_scheduler_recovery_rejects_cross_scope_suffix(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization(
            accounting_scope="manual_direct_scheduler_accounting_only",
        )
        arguments = self._manual_accounting_arguments(authorization)
        self._strand_manual_accounting_after_orion_append(arguments)
        lines = self.ledger.read_text(encoding="utf-8").splitlines()
        record = json.loads(lines[-1])
        record["accounting_scope"] = "manual_direct_srun_accounting_only"
        record["notes"] = (
            "Reviewed direct-srun accounting only; "
            "ineligible for scientific evidence."
        )
        record["event_sha256"] = ledger.record_sha256(record, "event_sha256")
        self.ledger.write_text(
            "".join(line + "\n" for line in lines[:-1])
            + ledger.canonical_json(record)
            + "\n",
            encoding="utf-8",
        )
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            with self.assertRaisesRegex(ValueError, "suffix event is invalid"):
                reconcile_manual_allocations(**arguments)

    def test_manual_direct_srun_accounting_recovery_rejects_pre_tranche_mirror_truncation(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        self._append_registered_probe()
        self._strand_manual_accounting_after_orion_append(arguments)
        local_marker, _ = incomplete_manual_accounting_marker_paths(
            self.ledger, self.mirror
        )
        self.assertEqual(
            json.loads(local_marker.read_text(encoding="utf-8"))[
                "pre_tranche_sequence_number"
            ],
            3,
        )
        lines = self.mirror.read_text(encoding="utf-8").splitlines()
        self.mirror.write_text(lines[0] + "\n", encoding="utf-8")
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            with self.assertRaisesRegex(ValueError, "pre-tranche publication loss"):
                reconcile_manual_allocations(**arguments)
        self.assertEqual(len(validate_primary_chain(self.mirror)), 1)
        self.assertTrue(local_marker.is_file())

    def test_manual_direct_srun_accounting_recovery_rejects_pre_tranche_receipt_truncation(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        self._append_registered_probe()
        self._strand_manual_accounting_after_orion_append(arguments)
        local_marker, _ = incomplete_manual_accounting_marker_paths(
            self.ledger, self.mirror
        )
        self.assertEqual(
            json.loads(local_marker.read_text(encoding="utf-8"))[
                "pre_tranche_sequence_number"
            ],
            3,
        )
        lines = self.receipts.read_text(encoding="utf-8").splitlines()
        self.receipts.write_text(lines[0] + "\n", encoding="utf-8")
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            with self.assertRaisesRegex(ValueError, "pre-tranche publication loss"):
                reconcile_manual_allocations(**arguments)
        self.assertEqual(
            len(self.receipts.read_text(encoding="utf-8").splitlines()), 1
        )
        self.assertTrue(local_marker.is_file())

    def test_manual_direct_srun_accounting_recovery_rejects_torn_orion_suffix(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        self._strand_manual_accounting_after_orion_append(arguments)
        with self.ledger.open("ab") as stream:
            stream.write(b'{"sequence_number":')
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
        ) as scheduler:
            with self.assertRaises(ValueError):
                reconcile_manual_allocations(**arguments)
        scheduler.assert_not_called()
        self.assertEqual(len(validate_primary_chain(self.mirror)), 1)

    def test_manual_direct_srun_accounting_recovery_rejects_corrupt_orion_suffix(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        self._strand_manual_accounting_after_orion_append(arguments)
        lines = self.ledger.read_text(encoding="utf-8").splitlines()
        record = json.loads(lines[-1])
        record["consumed_node_hours"] = 10.0
        record["event_sha256"] = ledger.record_sha256(record, "event_sha256")
        self.ledger.write_text(
            "".join(line + "\n" for line in lines[:-1])
            + ledger.canonical_json(record)
            + "\n",
            encoding="utf-8",
        )
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
        ) as scheduler:
            with self.assertRaisesRegex(ValueError, "ledger event usage differs"):
                reconcile_manual_allocations(**arguments)
        scheduler.assert_not_called()
        self.assertEqual(len(validate_primary_chain(self.mirror)), 1)

    def test_manual_direct_srun_accounting_recovery_rejects_unrelated_orion_suffix(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        self._strand_manual_accounting_after_orion_append(arguments)
        lines = self.ledger.read_text(encoding="utf-8").splitlines()
        record = json.loads(lines[-1])
        record["job_id"] = "4746335"
        record["event_sha256"] = ledger.record_sha256(record, "event_sha256")
        self.ledger.write_text(
            "".join(line + "\n" for line in lines[:-1])
            + ledger.canonical_json(record)
            + "\n",
            encoding="utf-8",
        )
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            with self.assertRaisesRegex(ValueError, "suffix event is invalid"):
                reconcile_manual_allocations(**arguments)
        self.assertEqual(len(validate_primary_chain(self.mirror)), 1)

    def test_manual_direct_srun_accounting_recovery_rejects_policy_binding_drift_before_mirror(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        self._strand_manual_accounting_after_orion_append(arguments)
        lines = self.ledger.read_text(encoding="utf-8").splitlines()
        record = json.loads(lines[-1])
        record["active_policy_sha256"] = "0" * 64
        record["event_sha256"] = ledger.record_sha256(record, "event_sha256")
        self.ledger.write_text(
            "".join(line + "\n" for line in lines[:-1])
            + ledger.canonical_json(record)
            + "\n",
            encoding="utf-8",
        )
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            with self.assertRaisesRegex(ValueError, "suffix event is invalid"):
                reconcile_manual_allocations(**arguments)
        self.assertEqual(len(validate_primary_chain(self.mirror)), 1)

    def test_manual_direct_srun_accounting_recovery_rejects_reviewed_qos_drift_before_mirror(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        self._strand_manual_accounting_after_orion_append(arguments)
        lines = self.ledger.read_text(encoding="utf-8").splitlines()
        record = json.loads(lines[-1])
        record["qos"] = "debug"
        record["event_sha256"] = ledger.record_sha256(record, "event_sha256")
        self.ledger.write_text(
            "".join(line + "\n" for line in lines[:-1])
            + ledger.canonical_json(record)
            + "\n",
            encoding="utf-8",
        )
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            with self.assertRaisesRegex(ValueError, "suffix event is invalid"):
                reconcile_manual_allocations(**arguments)
        self.assertEqual(len(validate_primary_chain(self.mirror)), 1)

    def test_manual_direct_srun_accounting_recovers_after_mirror_append(self) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        real_append_receipt = ledger._append_mirror_receipt
        interrupted = False

        def interrupt_before_receipt(
            receipts_jsonl: Path,
            event: dict[str, object],
            mirror_jsonl: Path,
            mirror_transport: str,
        ) -> dict[str, object]:
            nonlocal interrupted
            if (
                not interrupted
                and event.get("event_type") == "manual_allocation_reconciliation"
            ):
                interrupted = True
                raise RuntimeError("interrupted after mirror append")
            return real_append_receipt(
                receipts_jsonl, event, mirror_jsonl, mirror_transport
            )

        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ), patch("ledger._append_mirror_receipt", side_effect=interrupt_before_receipt):
            with self.assertRaisesRegex(RuntimeError, "after mirror append"):
                reconcile_manual_allocations(**arguments)
        self.assertEqual(len(validate_primary_chain(self.ledger)), 2)
        self.assertEqual(len(validate_primary_chain(self.mirror)), 2)
        self.assertEqual(len(self.receipts.read_text(encoding="utf-8").splitlines()), 1)
        self._retry_manual_accounting_and_require_clean_markers(arguments)

    def test_manual_direct_srun_accounting_recovers_after_receipt_before_csv(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ), patch("ledger.write_csv", side_effect=RuntimeError("interrupted before CSV")):
            with self.assertRaisesRegex(RuntimeError, "before CSV"):
                reconcile_manual_allocations(**arguments)
        self.assertEqual(len(validate_primary_chain(self.ledger)), 2)
        self.assertEqual(len(validate_primary_chain(self.mirror)), 2)
        self.assertEqual(len(self.receipts.read_text(encoding="utf-8").splitlines()), 2)
        self._retry_manual_accounting_and_require_clean_markers(arguments)

    def test_manual_direct_srun_accounting_historical_retry_repairs_csv_before_clear(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        self._retry_manual_accounting_and_require_clean_markers(arguments)
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ), patch(
            "reconcile_manual_frontier_allocations.write_csv",
            side_effect=RuntimeError("CSV publication failed"),
        ):
            with self.assertRaisesRegex(RuntimeError, "CSV publication failed"):
                reconcile_manual_allocations(**arguments)
        local_marker, mirror_marker = incomplete_manual_accounting_marker_paths(
            self.ledger, self.mirror
        )
        self.assertTrue(local_marker.is_file())
        assert mirror_marker is not None
        self.assertTrue(mirror_marker.is_file())
        self._retry_manual_accounting_and_require_clean_markers(arguments)

    def test_manual_direct_srun_accounting_recovers_one_sided_marker_publication(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        real_atomic_write = ledger.atomic_write_bytes_at
        marker_writes = 0

        def interrupt_second_marker(
            parent_descriptor: int,
            name: str,
            data: bytes,
            **kwargs: object,
        ) -> None:
            nonlocal marker_writes
            if name == ledger.INCOMPLETE_MANUAL_ACCOUNTING_MARKER_FILENAME:
                marker_writes += 1
                if marker_writes == 2:
                    raise RuntimeError("interrupted one-sided marker publication")
            real_atomic_write(parent_descriptor, name, data, **kwargs)

        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ), patch("ledger.atomic_write_bytes_at", side_effect=interrupt_second_marker):
            with self.assertRaisesRegex(RuntimeError, "one-sided marker"):
                reconcile_manual_allocations(**arguments)
        local_marker, mirror_marker = incomplete_manual_accounting_marker_paths(
            self.ledger, self.mirror
        )
        self.assertTrue(local_marker.is_file())
        assert mirror_marker is not None
        self.assertFalse(mirror_marker.exists())
        self._retry_manual_accounting_and_require_clean_markers(arguments)

    def test_manual_direct_srun_accounting_recovers_clear_one_marker_interruption(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        real_unlink = ledger.os.unlink
        marker_unlinks = 0

        def interrupt_second_marker_unlink(
            path: str, *args: object, **kwargs: object
        ) -> None:
            nonlocal marker_unlinks
            if path == ledger.INCOMPLETE_MANUAL_ACCOUNTING_MARKER_FILENAME:
                marker_unlinks += 1
                if marker_unlinks == 2:
                    raise RuntimeError("interrupted clearing marker pair")
            real_unlink(path, *args, **kwargs)

        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ), patch("ledger.os.unlink", side_effect=interrupt_second_marker_unlink):
            with self.assertRaisesRegex(RuntimeError, "clearing marker pair"):
                reconcile_manual_allocations(**arguments)
        local_marker, mirror_marker = incomplete_manual_accounting_marker_paths(
            self.ledger, self.mirror
        )
        self.assertFalse(local_marker.exists())
        assert mirror_marker is not None
        self.assertTrue(mirror_marker.is_file())
        self._retry_manual_accounting_and_require_clean_markers(arguments)

    def test_manual_direct_srun_accounting_legacy_partial_prefix_rejects_unrelated_suffix(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = self._manual_accounting_arguments(authorization)
        real_append = reconcile_manual_frontier_allocations.append_primary_event_locked
        append_count = 0

        def fail_second_append(*args: object, **kwargs: object) -> dict[str, object]:
            nonlocal append_count
            append_count += 1
            if append_count == 2:
                raise RuntimeError("interrupted")
            return real_append(*args, **kwargs)

        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ), patch(
            "reconcile_manual_frontier_allocations.append_primary_event_locked",
            side_effect=fail_second_append,
        ):
            with self.assertRaisesRegex(RuntimeError, "interrupted"):
                reconcile_manual_allocations(**arguments)
        local_marker, mirror_marker = incomplete_manual_accounting_marker_paths(
            self.ledger, self.mirror
        )
        local_marker.unlink()
        assert mirror_marker is not None
        mirror_marker.unlink()
        self._append_registered_probe()
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            with self.assertRaisesRegex(ValueError, "terminal ledger suffix"):
                reconcile_manual_allocations(**arguments)

    def test_manual_direct_srun_accounting_requires_frozen_authorization(self) -> None:
        authorization = self._write_manual_accounting_authorization()
        authorization.chmod(0o644)
        with self.assertRaisesRegex(ValueError, "not read-only"):
            reconcile_manual_allocations(
                authorization=authorization,
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_manual_accounting_rejects_nonstring_scope(self) -> None:
        authorization = self._write_manual_accounting_authorization(
            accounting_scope=[],
        )
        with self.assertRaisesRegex(
            ValueError, "Manual-accounting authorization scope is invalid"
        ):
            reconcile_manual_allocations(
                **self._manual_accounting_arguments(authorization)
            )

    def test_manual_accounting_rejects_unknown_scope(self) -> None:
        authorization = self._write_manual_accounting_authorization(
            accounting_scope="unknown_manual_scope",
        )
        with self.assertRaisesRegex(
            ValueError, "Manual-accounting authorization scope is invalid"
        ):
            reconcile_manual_allocations(
                **self._manual_accounting_arguments(authorization)
            )

    def test_manual_direct_srun_accounting_rejects_unbound_authorization(self) -> None:
        authorization = self._write_manual_accounting_authorization(
            bind_policy=False
        )
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
        ) as scheduler:
            with self.assertRaisesRegex(ValueError, "not bound by the promoted policy"):
                reconcile_manual_allocations(
                    authorization=authorization,
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        scheduler.assert_not_called()

    def test_policy_promotion_rejects_manual_accounting_authorization_hash_drift(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization(
            bind_policy=False
        )
        project_home_authorization = (
            self.project_home_root
            / "policy"
            / "manual_accounting_authorizations"
            / authorization.name
        )
        self._write_policy(
            manual_accounting_authorizations=[
                {
                    "authorization_id": authorization.stem,
                    "path": str(authorization),
                    "project_home_path": str(project_home_authorization),
                    "sha256": "0" * 64,
                }
            ]
        )
        with self.assertRaisesRegex(ValueError, "authorization bytes differ"):
            self._promote_policy()

    def test_manual_direct_srun_accounting_rejects_authorization_mirror_drift(
        self,
    ) -> None:
        authorization = self._write_manual_accounting_authorization()
        mirror = (
            self.project_home_root
            / "policy"
            / "manual_accounting_authorizations"
            / authorization.name
        )
        mirror.chmod(0o644)
        mirror.write_text("{}", encoding="utf-8")
        mirror.chmod(0o444)
        with self.assertRaisesRegex(ValueError, "mirror bytes differ"):
            reconcile_manual_allocations(
                authorization=authorization,
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

    def test_manual_direct_srun_accounting_rejects_nonempty_trusted_queue(self) -> None:
        authorization = self._write_manual_accounting_authorization()
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            return_value="4747000\n",
        ):
            with self.assertRaisesRegex(ValueError, "empty trusted Frontier queue"):
                reconcile_manual_allocations(
                    authorization=authorization,
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )

    def test_manual_direct_srun_accounting_rejects_pending_marker(self) -> None:
        authorization = self._write_manual_accounting_authorization()
        (self.pic_root / "ledger" / "pending_submission.json").write_text(
            "{}\n", encoding="utf-8"
        )
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
        ) as scheduler:
            with self.assertRaisesRegex(ValueError, "pending scheduler submission"):
                reconcile_manual_allocations(
                    authorization=authorization,
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        scheduler.assert_not_called()

    def test_manual_direct_srun_accounting_rejects_active_reservation(self) -> None:
        authorization = self._write_manual_accounting_authorization()
        self._append_registered_probe(final_event_type="reservation")
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
        ) as scheduler:
            with self.assertRaisesRegex(ValueError, "active reservation"):
                reconcile_manual_allocations(
                    authorization=authorization,
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        scheduler.assert_not_called()

    def test_manual_direct_srun_accounting_rejects_prior_job_event(self) -> None:
        authorization = self._write_manual_accounting_authorization()
        self._append_registered_probe(
            final_event_type="reconciliation",
            job_id="4746332",
        )
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
            side_effect=self._manual_accounting_scheduler_output,
        ):
            with self.assertRaisesRegex(ValueError, "prior ledger event"):
                reconcile_manual_allocations(
                    authorization=authorization,
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=self.control_plane_dir,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )

    def test_manual_direct_srun_accounting_requires_promoted_successor(self) -> None:
        authorization = self._write_manual_accounting_authorization()
        successor = self._publish_test_control_plane_successor(self.pic_root)
        project_home_successor = self._publish_test_control_plane_successor(
            self.project_home_root
        )
        self.assertEqual(project_home_successor.name, successor.name)
        with patch.object(
            reconcile_manual_frontier_allocations.subprocess,
            "check_output",
        ) as scheduler:
            with self.assertRaisesRegex(ValueError, "Active-policy promotion"):
                reconcile_manual_allocations(
                    authorization=authorization,
                    ledger_jsonl=self.ledger,
                    ledger_csv=self.csv,
                    receipts_jsonl=self.receipts,
                    mirror_jsonl=self.mirror,
                    control_plane_dir=successor,
                    authorized_pic_root=self.pic_root,
                    authorized_project_home_root=self.project_home_root,
                )
        scheduler.assert_not_called()

    def test_manual_direct_srun_accounting_rejects_scheduler_field_drift(self) -> None:
        authorization = self._write_manual_accounting_authorization()
        arguments = {
            "authorization": authorization,
            "ledger_jsonl": self.ledger,
            "ledger_csv": self.csv,
            "receipts_jsonl": self.receipts,
            "mirror_jsonl": self.mirror,
            "control_plane_dir": self.control_plane_dir,
            "authorized_pic_root": self.pic_root,
            "authorized_project_home_root": self.project_home_root,
        }
        scheduler_outputs = [
            "4746332|FAILED|5|1|unexpected|ast207|batch|normal\n"
            "4746335|COMPLETED|7|1||ast207|batch|normal\n",
            "4746332|FAILED|5|1||other|batch|normal\n"
            "4746335|COMPLETED|7|1||ast207|batch|normal\n",
            "4746332|FAILED|5|1||ast207|other|normal\n"
            "4746335|COMPLETED|7|1||ast207|batch|normal\n",
            "4746332|FAILED|5|1||ast207|batch|debug\n"
            "4746335|COMPLETED|7|1||ast207|batch|normal\n",
            "4746332|FAILED|-1|1||ast207|batch|normal\n"
            "4746335|COMPLETED|7|1||ast207|batch|normal\n",
            "4746332|FAILED|5|1||ast207|batch|normal\n",
            "4746332|FAILED|5|1||ast207|batch|normal\n"
            "4746332|FAILED|5|1||ast207|batch|normal\n"
            "4746335|COMPLETED|7|1||ast207|batch|normal\n",
            "4746332|FAILED|5|1||ast207|batch|normal\n"
            "4746335|COMPLETED|7|1||ast207|batch|normal\n"
            "4746999|COMPLETED|1|1||ast207|batch|normal\n",
        ]
        for scheduler_output in scheduler_outputs:
            with self.subTest(scheduler_output=scheduler_output):
                with patch.object(
                    reconcile_manual_frontier_allocations.subprocess,
                    "check_output",
                    side_effect=["", scheduler_output],
                ):
                    with self.assertRaises(ValueError):
                        reconcile_manual_allocations(**arguments)


if __name__ == "__main__":
    unittest.main()
