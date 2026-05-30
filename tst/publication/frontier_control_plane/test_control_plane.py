#!/usr/bin/env python3
"""Tests for immutable Frontier PIC submission snapshots."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
import json
import os
from pathlib import Path
import pwd
import inspect
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import uuid

from control_plane_common import git_tree_sha1_from_archive, record_for_role, sha256
from control_plane_common import source_bundle_sha256
from control_plane_common import verify_installed_control_plane
from control_plane_common import TRUSTED_GIT, TRUSTED_PYTHON
from control_plane_common import TRUSTED_SACCT, TRUSTED_SBATCH, TRUSTED_SCANCEL
from control_plane_common import TRUSTED_SCONTROL, TRUSTED_SQUEUE
from create_clean_candidate_freeze import create_freeze, _validated_submodules
from create_pre_submit_manifest import create_manifest
from initialize_frontier_ledger import initialize_from_policy
from install_control_plane import install
from launch_trampoline import launch
from ledger import accounting, validate_primary_chain
from promote_active_policy import promote
from reconcile_frontier_job import reconcile
from validate_and_reserve_frontier_job import mark_submitted, repair_reservation_attachments
from validate_and_reserve_frontier_job import reservation_bound_manifest
from validate_and_reserve_frontier_job import reserve, transition
from verify_compute_node_snapshot import verify


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
        self._write_policy()
        self._promote_policy()
        self._write_config()
        self.ledger = self.pic_root / "ledger" / "node_hours.jsonl"
        self.csv = self.pic_root / "ledger" / "node_hours.csv"
        self.receipts = self.pic_root / "ledger" / "mirror_receipts.jsonl"
        self.mirror = self.project_home_root / "ledger" / "node_hours.jsonl"
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

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write(self, name: str, text: str) -> Path:
        path = self.sources / name
        path.write_text(text, encoding="utf-8")
        return path

    def _utc(self, value: datetime) -> str:
        return value.replace(microsecond=0).isoformat().replace("+00:00", "Z")

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

    def _write_policy(
        self,
        *,
        science_submission_freeze: dict[str, object] | None = None,
        **storage_overrides: object,
    ) -> None:
        storage = {
            "installed_control_plane_version": self.control_plane_version,
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
            "ledger_genesis_allowed": True,
        }
        storage.update(storage_overrides)
        policy = {
            "schema_version": 1,
            "frontier": {
                "account": "AST207",
                "partition": "batch",
                "simulation_root": str(self.pic_root),
                "maximum_node_hours": 10000.0,
                "serial_pic_submissions": True,
            },
            "science_submission_freeze": science_submission_freeze or {
                "status": "pending_clean_candidate_freeze",
            },
            "frontier_admission_smoke": {
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
            },
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
            "launch_contract": {
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
            },
            "git_commit": "abc123",
            "evidence_class": "frontier_f0_admission_smoke_candidate",
            "physical_mode": "extended_mhd_pic_parser_contract",
            "selected_qos": "debug",
            "qos_selection_reason": "debug_available",
            "site_policy_checked_utc": self._utc(now),
            "registered_short_nonproduction": True,
            "artifact_dir": str(self.pic_root / "runs" / "snapshot"),
        }
        config.update(overrides)
        self.config.write_text(json.dumps(config), encoding="utf-8")

    def _clean_source(self, name: str) -> Path:
        source_root = self.root / name
        source_root.mkdir()
        subprocess.run(["git", "init", str(source_root)], check=True, capture_output=True)
        (source_root / "tracked.txt").write_text("tracked\n", encoding="utf-8")
        subprocess.run(["git", "-C", str(source_root), "add", "tracked.txt"], check=True)
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

    def _build_profile(self, source_root: Path, build: Path, profile_id: str) -> tuple[Path, Path]:
        build.mkdir(parents=True)
        executable = build / "athena"
        executable.write_text("built executable\n", encoding="utf-8")
        declared_archive = build / "declared-source.tar"
        subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "archive",
                "--format=tar",
                f"--output={declared_archive}",
                "HEAD",
            ],
            check=True,
        )
        submodules = []
        declared_submodule_archives = []
        for index, record in enumerate(_validated_submodules(source_root)):
            declared_submodule_archive = build / f"declared-submodule-{index:04d}.tar"
            module_root = source_root.joinpath(*Path(record["path"]).parts)
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(module_root),
                    "archive",
                    "--format=tar",
                    f"--output={declared_submodule_archive}",
                    record["git_commit"],
                ],
                check=True,
            )
            declared_submodule_archives.append(declared_submodule_archive)
            submodules.append(
                {
                    "path": record["path"],
                    "archive_sha256": sha256(declared_submodule_archive),
                    "git_commit": record["git_commit"],
                    "git_tree": record["git_tree"],
                }
            )
        profile = build / "build_profile.json"
        archive_sha256 = sha256(declared_archive)
        profile.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "profile_id": profile_id,
                    "source_archive_sha256": archive_sha256,
                    "source_bundle_sha256": source_bundle_sha256(
                        archive_sha256, submodules
                    ),
                    "toolchain": "test-toolchain",
                    "build_command": "cmake --build build",
                    "executable_sha256": sha256(executable),
                    "submodules": submodules,
                }
            ),
            encoding="utf-8",
        )
        declared_archive.unlink()
        for declared_submodule_archive in declared_submodule_archives:
            declared_submodule_archive.unlink()
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
        self._add_submodule(source_root, "nested")
        executable, profile = self._build_profile(
            source_root, self.pic_root / "candidate-build", "hip-mpi-release-paper-pic"
        )
        manifest_path = create_freeze(
            source_root=source_root,
            executable=executable,
            build_profile=profile,
            build_profile_id="hip-mpi-release-paper-pic",
            freeze_id="03a7bd9a-7d4c-4e37-a12b-46de3817eff2",
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
        )
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        frozen_executable = Path(str(manifest["build"]["executable_path"]))
        git_commit = str(manifest["source"]["git_commit"])
        if authorize:
            self._write_policy(
                science_submission_freeze={
                    "status": "authorized",
                    "manifest_path": str(manifest_path),
                    "manifest_sha256": sha256(manifest_path),
                }
            )
            self._promote_policy()
        return manifest_path, frozen_executable, git_commit

    def _write_science_config(self, *, authorize: bool, **overrides: object) -> Path:
        manifest, executable, git_commit = self._clean_candidate(authorize=authorize)
        config = {
            "campaign": "f1_gpu_gyro",
            "test_id": "pic_relativistic_gyro_paper",
            "submission_scope": "registered_science",
            "git_commit": git_commit,
            "evidence_class": "frontier_f1_registered_science",
            "physical_mode": "paper_test_particle",
            "executable": str(executable),
            "clean_candidate_manifest": str(manifest),
        }
        config.update(overrides)
        self._write_config(**config)
        return manifest

    def _create_manifest(self) -> Path:
        return create_manifest(
            self.config,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
        )

    def _reserve(
        self,
        manifest_path: Path,
        cap: float = 10000.0,
        reservation_id: str = "89c76745-6c37-47f7-9847-800a98a47c9b",
    ) -> dict[str, object]:
        return reserve(
            manifest_path=manifest_path,
            ledger_jsonl=self.ledger,
            ledger_csv=self.csv,
            receipts_jsonl=self.receipts,
            mirror_jsonl=self.mirror,
            node_hour_cap=cap,
            reservation_id=reservation_id,
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )

    def _attach(self, reservation_id: str, job_id: str = "12345") -> None:
        scheduler = (
            f"JobId={job_id} JobState=PENDING Account=AST207 "
            f"Comment=pic-reservation={reservation_id}"
        )
        with patch(
            "validate_and_reserve_frontier_job._scheduler_job_output",
            return_value=scheduler,
        ):
            mark_submitted(
                reservation_id=reservation_id,
                job_id=job_id,
                ledger_jsonl=self.ledger,
                authorized_pic_root=self.pic_root,
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

    def _launch(
        self,
        manifest_path: Path,
        reservation: dict[str, object],
        *,
        runner: object = subprocess.run,
    ) -> None:
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
            launch(
                manifest_path=manifest_path,
                manifest_sha256=str(reservation["manifest_sha256"]),
                job_script_sha256=str(job_script["sha256"]),
                executable_sha256=str(executable["sha256"]),
                reservation_id=str(reservation["reservation_id"]),
                submission_id=self.submission_id,
                ledger_jsonl=self.ledger,
                receipts_jsonl=self.receipts,
                mirror_jsonl=self.mirror,
                runner=runner,
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
                authorized_project_home_root=self.project_home_root,
            )

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
            "create_pre_submit_manifest.os.replace", side_effect=OSError("rename failed")
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

    def test_installed_control_plane_rename_failure_cleans_staging(self) -> None:
        target = self.root / "failed-install"
        with patch("install_control_plane.os.replace", side_effect=OSError("rename failed")):
            with self.assertRaises(OSError):
                install(target)
        self.assertFalse(list((target / "control_plane").glob(".tmp-*")))

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

    def test_promoted_policy_anchor_is_mirrored_and_read_only(self) -> None:
        policy = self.pic_root / "policy" / "storage_policy.json"
        mirror_policy = self.project_home_root / "policy" / "storage_policy.json"
        promotion = self.pic_root / "policy" / "active_promotion.json"
        mirror_promotion = self.project_home_root / "policy" / "active_promotion.json"
        self.assertEqual(policy.read_bytes(), mirror_policy.read_bytes())
        self.assertEqual(promotion.read_bytes(), mirror_promotion.read_bytes())
        for path in [policy, mirror_policy, promotion, mirror_promotion]:
            self.assertFalse(bool(path.stat().st_mode & 0o222))

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
        self.assertEqual(
            calls[0][0][0],
            str(self.control_plane_dir / "launch_with_frontier_profile.sh"),
        )
        self.assertEqual(calls[0][0][1], "/usr/bin/srun")
        self.assertEqual(calls[0][0][7], str(executable["path"]))
        self.assertTrue(calls[0][1]["check"])

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
        self.assertEqual(
            commands[0][0],
            str(self.control_plane_dir / "launch_with_frontier_profile.sh"),
        )
        self.assertEqual(commands[0][1], "/usr/bin/srun")
        self.assertNotIn("/bin/bash", commands[0])
        self.assertNotIn("/bin/true", commands[0])

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
        self.assertEqual(
            commands[0][0],
            str(self.control_plane_dir / "launch_with_frontier_profile.sh"),
        )
        self.assertEqual(commands[0][1], "/usr/bin/srun")
        self.assertNotIn(str(template), commands[0])

    def test_manifest_rejects_arbitrary_shell_launch_action(self) -> None:
        self._write_config()
        config = json.loads(self.config.read_text(encoding="utf-8"))
        config["launch_contract"]["actions"][0]["kind"] = "shell"
        config["launch_contract"]["actions"][0]["executable"] = "/bin/true"
        self.config.write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaises(ValueError):
            self._create_manifest()

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
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        calls: list[tuple[list[str], dict[str, object]]] = []

        def runner(command: list[str], **kwargs: object) -> None:
            calls.append((command, kwargs))

        self._launch(manifest_path, reservation, runner=runner)
        self.assertEqual(len(calls), 1)
        self.assertEqual(
            calls[0][0][0],
            str(self.control_plane_dir / "launch_with_frontier_profile.sh"),
        )
        self.assertEqual(calls[0][0][1], "/usr/bin/srun")
        artifact_dir = self.pic_root / "runs" / "snapshot"
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
            authorized_pic_root=self.pic_root,
            authorized_project_home_root=self.project_home_root,
        )
        self.assertEqual(result, "cleared_unappended_reservation_intent")
        self.assertFalse((self.pic_root / "ledger" / "pending_submission.json").exists())
        self.assertEqual(len(validate_primary_chain(self.ledger)), 1)

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
            with self.assertRaises(ValueError):
                mark_submitted(
                    reservation_id=reservation_id,
                    job_id="12345",
                    ledger_jsonl=self.ledger,
                    authorized_pic_root=self.pic_root,
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
                    authorized_pic_root=self.pic_root,
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
                    authorized_pic_root=self.pic_root,
                )

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

    def test_registered_science_rejects_pending_clean_candidate_freeze(self) -> None:
        self._write_science_config(authorize=False)
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_registered_science_accepts_exact_authorized_clean_candidate(self) -> None:
        self._write_science_config(authorize=True)
        manifest_path = self._create_manifest()
        reservation = self._reserve(manifest_path)
        candidate_sha256 = str(reservation["clean_candidate_manifest_sha256"])
        self.assertEqual(len(candidate_sha256), 64)
        self.assertEqual(reservation["submission_scope"], "registered_science")
        self.assertIn("clean_candidate_manifest_sha256", self.csv.read_text())

    def test_registered_science_rejects_policy_digest_mismatch(self) -> None:
        candidate = self._write_science_config(authorize=True)
        self._write_policy(
            science_submission_freeze={
                "status": "authorized",
                "manifest_path": str(candidate),
                "manifest_sha256": "0" * 64,
            }
        )
        self._promote_policy()
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_registered_science_rejects_missing_candidate_manifest(self) -> None:
        candidate = self._write_science_config(authorize=True)
        manifest_path = self._create_manifest()
        candidate.parent.chmod(0o755)
        candidate.unlink()
        candidate.parent.chmod(0o555)
        with self.assertRaises(FileNotFoundError):
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
        self._write_policy(
            science_submission_freeze={
                "status": "authorized",
                "manifest_path": str(candidate),
                "manifest_sha256": sha256(candidate),
            }
        )
        self._promote_policy()
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
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
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
        )
        candidate = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(candidate["source"]["worktree_status"], "clean")
        self.assertEqual(candidate["source"]["submodule_status"], "absent")
        self.assertEqual(candidate["build"]["executable_sha256"], sha256(executable))
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
        with self.assertRaises(ValueError):
            create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="test-profile",
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
        profile.write_text(json.dumps(value), encoding="utf-8")
        with self.assertRaises(ValueError):
            create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="test-profile",
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
            )

    def test_clean_candidate_rename_failure_cleans_read_only_staging(self) -> None:
        source_root = self._clean_source("candidate-rename-failure-source")
        executable, profile = self._build_profile(
            source_root, self.pic_root / "candidate-rename-failure-build", "test-profile"
        )
        candidate_root = self.pic_root / "clean_candidates"
        with patch(
            "create_clean_candidate_freeze.os.replace",
            side_effect=OSError("rename failed"),
        ):
            with self.assertRaises(OSError):
                create_freeze(
                    source_root=source_root,
                    executable=executable,
                    build_profile=profile,
                    build_profile_id="test-profile",
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
            control_plane_dir=self.control_plane_dir,
            authorized_pic_root=self.pic_root,
        )
        candidate = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(candidate["source"]["submodule_status"], "clean_pinned_archived")
        self.assertEqual(len(candidate["source"]["submodules"]), 1)
        archive = Path(str(candidate["source"]["submodules"][0]["archive_path"]))
        self.assertEqual(archive, manifest_path.parent / "submodules" / "0000.tar")
        self.assertFalse(bool(archive.stat().st_mode & 0o222))

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
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
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
        executable, profile = self._build_profile(
            source_root, self.pic_root / "symlink-payload-build", "test-profile"
        )
        with self.assertRaises(ValueError):
            create_freeze(
                source_root=source_root,
                executable=executable,
                build_profile=profile,
                build_profile_id="test-profile",
                control_plane_dir=self.control_plane_dir,
                authorized_pic_root=self.pic_root,
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
        self._write_policy(ledger_genesis_allowed=False)
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_initializer_rejects_unapproved_control_plane_version(self) -> None:
        self._write_policy(installed_control_plane_version="0" * 64)
        with self.assertRaises(ValueError):
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

    def test_reservation_rejects_locked_storage_policy(self) -> None:
        self._write_policy(ledger_genesis_allowed=False)
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_reservation_rejects_unreviewed_ledger_mirror_transport(self) -> None:
        self._write_policy(project_home_ledger_mirror_transport="dtn_rsync")
        with self.assertRaises(ValueError):
            self._promote_policy()

    def test_reservation_ignores_unpromoted_source_policy_override(self) -> None:
        manifest_path = self._create_manifest()
        self._write_policy(ledger_genesis_allowed=False)
        reservation = self._reserve(manifest_path)
        self.assertEqual(reservation["state"], "reserved")

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
                    control_plane_dir=alias,
                    authorized_pic_root=self.pic_root,
                )
        git.assert_not_called()

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
        with patch.object(
            validate_and_reserve_frontier_job.subprocess,
            "check_output",
            return_value="",
        ) as scheduler:
            validate_and_reserve_frontier_job._scheduler_job_output("123")
        self.assertEqual(scheduler.call_args.args[0][0], TRUSTED_SCONTROL)
        with patch.object(
            reconcile_frontier_job.subprocess,
            "check_output",
            return_value="",
        ) as accounting_call:
            with self.assertRaises(ValueError):
                reconcile_frontier_job._scheduler_result("123", "reservation")
        self.assertEqual(accounting_call.call_args.args[0][0], TRUSTED_SACCT)
        self.assertIn(
            "[TRUSTED_GIT, \"-C\", str(source_root), *arguments]",
            inspect.getsource(create_clean_candidate_freeze._git),
        )
        wrapper = Path(__file__).with_name("submit_frontier_job.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn('PYTHON="/opt/cray/pe/python/3.11.7/bin/python3"', wrapper)
        self.assertIn('SBATCH="/usr/bin/sbatch"', wrapper)
        self.assertIn('SCANCEL="/usr/bin/scancel"', wrapper)

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

    def test_reservation_rejects_stale_site_policy_timestamp(self) -> None:
        stale = datetime.now(timezone.utc) - timedelta(days=2)
        self._write_config(site_policy_checked_utc=self._utc(stale))
        manifest_path = self._create_manifest()
        with self.assertRaises(ValueError):
            self._reserve(manifest_path)

    def test_concurrent_reservations_are_serialized(self) -> None:
        first = self._create_manifest()
        self._write_config(
            submission_id=str(uuid.uuid4()),
            artifact_dir=str(self.pic_root / "runs" / "snapshot-two"),
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


if __name__ == "__main__":
    unittest.main()
