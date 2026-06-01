#!/usr/bin/env python3
"""Tests for immutable Frontier PIC submission snapshots."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import pwd
import inspect
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import uuid

import install_control_plane
import launch_trampoline
import ledger
import reconcile_frontier_job
import reconcile_manual_frontier_allocations
import terminal_recovery_handoff
from control_plane_common import atomic_write_bytes, durable_mkdir_parents
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
from control_plane_common import validate_launch_contract, verify_installed_control_plane
from control_plane_common import verify_historical_installed_control_plane
from control_plane_common import verify_snapshot_files
from control_plane_common import PRODUCTION_RUNTIME_MODULEPATH
from control_plane_common import TRUSTED_GIT, TRUSTED_PYTHON
from control_plane_common import TRUSTED_SACCT, TRUSTED_SBATCH, TRUSTED_SCANCEL
from control_plane_common import TRUSTED_SCONTROL, TRUSTED_SQUEUE
from create_clean_candidate_freeze import _authorized_source_path
from create_clean_candidate_freeze import create_freeze, _validated_submodules
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

    def _publish_test_control_plane_successor(self, root: Path) -> Path:
        staging = root / "control_plane" / "test-successor-staging"
        shutil.copytree(self.control_plane_dir, staging)
        for path in staging.iterdir():
            path.chmod(path.stat().st_mode | 0o200)
        schema = staging / "control_plane.schema.json"
        schema.write_text(
            schema.read_text(encoding="utf-8") + "\n",
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
        promote(
            self.policy,
            control_plane_dir=successor,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )

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
        bind_policy: bool = True,
    ) -> Path:
        parent = self.pic_root / "policy" / "manual_accounting_authorizations"
        parent.mkdir(parents=True, exist_ok=True)
        path = parent / f"{authorization_id}.json"
        data = json.dumps(
            {
                "schema_version": 1,
                "authorization_id": authorization_id,
                "accounting_scope": "manual_direct_srun_accounting_only",
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
        storage = {
            "installed_control_plane_version": self.control_plane_version,
            "staged_control_plane_candidate_version": self.control_plane_version,
            "installed_control_plane_lifecycle": "paired_installed_reviewed_generation",
            "orion_simulation_root_preflight": {"status": "passed"},
            "project_home_mirror_root": str(self.project_home_root),
            "project_home_usage": [
                "small_append_only_ledger_and_control_plane_mirror",
            ],
            "project_home_retention_role": "operational_ledger_mirror_only",
            "project_home_ledger_mirror_transport": "filesystem_copy",
            "project_home_preflight": {"status": "passed"},
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

    def _promote_policy(self) -> None:
        promote(
            self.policy,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )

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
            / "inputs/publication/pic_parallel_shock_section54_paper.athinput"
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
                                "pic_parallel_shock_section54_paper.athinput"
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
        self, *, authorize: bool
    ) -> tuple[Path, Path, str]:
        source_root = self._clean_source("candidate-source")
        self.authorized_clean_candidate_source_root = source_root
        self._add_submodule(source_root, "nested")
        executable, profile = self._build_profile(
            source_root, self.pic_root / "candidate-build", "hip-mpi-release-paper-pic"
        )
        manifest_path = create_freeze(
            source_root=source_root,
            executable=executable,
            build_profile=profile,
            build_profile_id="hip-mpi-release-paper-pic",
            prepared_artifact_inventory=self._prepared_artifact_inventory(),
            freeze_id="03a7bd9a-7d4c-4e37-a12b-46de3817eff2",
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
        manifest, executable, git_commit = self._clean_candidate(authorize=authorize)
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
        self._promote_policy()

    def _create_manifest(self, *, control_plane_dir: Path | None = None) -> Path:
        return create_manifest(
            self.config,
            control_plane_dir=control_plane_dir or self.control_plane_dir,
            authorized_pic_root=self.pic_root,
        )

    def _reserve(
        self,
        manifest_path: Path,
        cap: float = 10000.0,
        reservation_id: str = "89c76745-6c37-47f7-9847-800a98a47c9b",
        control_plane_dir: Path | None = None,
        *,
        patch_clean_candidate_bundle: bool = True,
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
            return reserve(
                manifest_path=manifest_path,
                ledger_jsonl=self.ledger,
                ledger_csv=self.csv,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                node_hour_cap=cap,
                reservation_id=reservation_id,
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
                            production_root
                        )
                tracked.assert_called_once()
                self.assertEqual(
                    tracked.call_args.kwargs["env"], trusted_git_environment()
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
        for path in [policy, mirror_policy, promotion, mirror_promotion]:
            self.assertFalse(bool(path.stat().st_mode & 0o222))

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
        self._promote_policy()
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
        self._promote_policy()
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
        prohibited = [
            value["required"][0]
            for value in branch["else"]["not"]["anyOf"]
        ]
        self.assertIn("registered_science_authorization_id", prohibited)

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
        self._promote_policy()
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

    def test_registered_science_manifest_path_swap_uses_stable_candidate_bytes(
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
            reservation = self._reserve(manifest_path)
        self.assertTrue(swapped)
        self.assertEqual(reservation["submission_scope"], "registered_science")

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
        self._promote_policy()
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
                            "pic_parallel_shock_section54_paper.athinput"
                        ),
                        "sha256": sha256(
                            source_root
                            / "inputs/publication/"
                            "pic_parallel_shock_section54_paper.athinput"
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
                "inputs/publication/pic_parallel_shock_section54_paper.athinput"
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
        nested = self._clean_source("recursive-nested-source")
        self._add_submodule(nested, "child")
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
        executable, profile = self._build_profile(
            source_root, self.pic_root / "recursive-submodule-build", "test-profile"
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
            last_preflight_utc="2026-05-30T22:52:40Z",
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
                "method": "local_create_write_sync_remove_probe",
            },
        )
        self._promote_policy()

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

    def test_active_policy_snapshot_accepts_configured_project_home_root_alias(
        self,
    ) -> None:
        project_home_alias = self.root / "project_home_alias"
        project_home_alias.symlink_to(self.project_home_root, target_is_directory=True)
        self._write_policy(project_home_mirror_root=str(project_home_alias))
        promote(
            self.policy,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=project_home_alias,
        )
        policy, snapshot = require_storage_policy_unlock_snapshot(
            control_plane_version=self.control_plane_version,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=project_home_alias,
        )
        self.assertEqual(
            policy["olcf_side_storage"]["project_home_mirror_root"],
            str(project_home_alias),
        )
        self.assertEqual(len(snapshot["active_policy_sha256"]), 64)
        self.assertEqual(len(snapshot["active_promotion_sha256"]), 64)

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
        for path in policy_root.iterdir():
            path.unlink()
        policy_root.rmdir()
        outside = self.root / "outside-policy"
        outside.mkdir()
        policy_root.symlink_to(outside, target_is_directory=True)
        with self.assertRaises(ValueError):
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
            [str(alias / "submit_frontier_job.sh"), "unused-manifest"],
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
            if path.name != "test_control_plane.py"
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
