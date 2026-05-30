#!/usr/bin/env python3
"""Regression tests for source-controlled PIC readiness registries."""

from __future__ import annotations

import base64
import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tst" / "publication" / "frontier_control_plane"))

from control_plane_common import CONTROL_PLANE_FILES
from control_plane_common import inventory_digest
from control_plane_common import launch_contract_sha256
from control_plane_common import validate_launch_contract
from tst.publication.pic_qualification_manifest import validate_schema


READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
SCHEMA_DIR = READINESS_DIR / "schemas"
CONTROL_PLANE_DIR = REPO_ROOT / "tst" / "publication" / "frontier_control_plane"
PAPER_TEX = (
    REPO_ROOT
    / "docs"
    / "reference_paper"
    / "arXiv-2304.10568v1"
    / "mnras_template.tex"
)

CLAIM_CLASSES = {
    "unit/regression",
    "engineering_proxy",
    "physics_validation",
    "sun_bai_2023_reproduction",
    "athenak_production_mode",
    "cross_code_comparison",
    "scoped_state_of_the_art",
    "unsupported",
}

REQUIRED_EXTENSION_CLAIMS = {
    "CLAIM-EXT-HALL-BELL-001",
    "CLAIM-EXT-CRSI-IN-DAMPING-001",
    "CLAIM-EXT-CRPAI-TRANSPORT-001",
    "CLAIM-STATEART-CRPAI-SCATTERING-001",
    "CLAIM-RELEASE-EXTENDED-MHD-PIC-001",
}


def _load(name: str) -> dict[str, object]:
    return json.loads((READINESS_DIR / name).read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _validation_manifest_schema() -> dict[str, object]:
    return json.loads(
        (SCHEMA_DIR / "validation_manifest.schema.json").read_text(
            encoding="utf-8"
        )
    )


def _minimum_validation_manifest() -> dict[str, object]:
    return {
        "schema_version": 1,
        "manifest_id": "host-scaffold-001",
        "created_utc": "2026-05-30T12:00:00Z",
        "claim_ids": ["CLAIM-PAPER-GYRO-001"],
        "test_id": "pic_relativistic_gyro_paper",
        "evidence_class": "unit/regression",
        "physical_mode": "paper_mhd_pic",
        "git": {
            "commit": "0" * 40,
            "tree": "0" * 40,
            "status": [],
            "source_archive": {
                "path": "source.tar",
                "sha256": "0" * 64,
            },
            "source_commit": {
                "path": "source.commit",
                "sha256": "0" * 64,
            },
            "source_bundle_sha256": "0" * 64,
            "submodule_status": "absent",
            "submodules": [],
            "clean_candidate_manifest": {
                "path": "clean_candidate_manifest.json",
                "sha256": "0" * 64,
            },
            "clean_candidate_build_profile": {
                "path": "build_profile.json",
                "sha256": "0" * 64,
            },
            "clean_candidate_build_profile_receipt": {
                "path": "profile_receipt.json",
                "sha256": "0" * 64,
            },
            "clean_candidate_build_provenance": {
                label: {
                    "path": f"build_provenance/{filename}",
                    "sha256": "0" * 64,
                }
                for label, filename in {
                    "configure_log": "configure.log",
                    "build_log": "build.log",
                    "cmake_cache": "CMakeCache.txt",
                    "module_list": "modules.txt",
                    "toolchain": "toolchain.txt",
                    "build_invocations": "build-invocations.json",
                    "git_status_preconfigure": "git_status.preconfigure.txt",
                    "git_status": "git_status.txt",
                    "submodule_status": "submodule_status.txt",
                    "environment_allowlist": "environment.allowlist.txt",
                    "build_environment": "build-environment.json",
                }.items()
            },
        },
        "authorization": {
            "control_plane_version": "0" * 64,
            "build_profile_control_plane_version": "0" * 64,
            "clean_candidate_manifest_path": (
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/"
                "clean_candidates/00000000-0000-0000-0000-000000000000/"
                "clean_candidate_manifest.json"
            ),
            "clean_candidate_manifest_sha256": "0" * 64,
            "active_policy": {
                "path": "active_policy.json",
                "sha256": "0" * 64,
            },
            "active_promotion": {
                "path": "active_promotion.json",
                "sha256": "0" * 64,
            },
        },
        "executable": {
            "path": "/tmp/athena",
            "sha256": "0" * 64,
            "cmake_cache": {"path": "/tmp/CMakeCache.txt", "sha256": "0" * 64},
            "modules": {"path": "/tmp/modules.txt", "sha256": "0" * 64},
            "environment_allowlist": {
                "path": "/tmp/environment.txt",
                "sha256": "0" * 64,
            },
        },
        "parameters": {},
        "oracle": {
            "kind": "analytic",
            "reference": "bounded host scaffold",
            "tolerances": {"relative_error": 1.0e-6},
        },
        "metrics": [{"name": "relative_error", "value": 0.0}],
        "resources": {
            "platform": "host",
            "artifact_root": "/tmp/pic-readiness",
        },
        "artifacts": [{"path": "metrics.json", "sha256": "0" * 64}],
        "review": {
            "reviewer": "pending external review",
            "disposition": "pending external review",
        },
    }


def _replace_nested(value: dict[str, object], path: tuple[object, ...],
                    replacement: object) -> None:
    target = value
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = replacement


class PicReadinessRegistryTests(unittest.TestCase):
    def test_json_documents_parse(self) -> None:
        paths = sorted(READINESS_DIR.glob("*.json"))
        paths += sorted(SCHEMA_DIR.glob("*.json"))
        self.assertTrue(paths)
        for path in paths:
            with self.subTest(path=path):
                json.loads(path.read_text(encoding="utf-8"))

    def test_storage_policy_records_authorized_frontier_boundary(self) -> None:
        policy = _load("storage_policy.json")
        frontier = policy["frontier"]
        storage = policy["olcf_side_storage"]
        long_term = policy["long_term_storage"]
        self.assertEqual(frontier["maximum_node_hours"], 10000)
        self.assertEqual(frontier["partition"], "batch")
        self.assertTrue(frontier["serial_pic_submissions"])
        self.assertEqual(
            frontier["simulation_root"],
            "/lustre/orion/ast207/proj-shared/dfielding/PIC",
        )
        self.assertFalse(storage["ledger_genesis_allowed"])
        self.assertEqual(storage["ledger_genesis"]["status"], "initialized")
        self.assertEqual(
            storage["ledger_genesis"]["mirror_transport"],
            "filesystem_copy",
        )
        self.assertEqual(
            storage["orion_bulk_evidence_root"],
            "/lustre/orion/ast207/proj-shared/dfielding/PIC",
        )
        self.assertEqual(
            storage["project_home_retention_role"],
            "operational_ledger_mirror_only",
        )
        source_alias_candidate = _load(
            "q027_control_plane_source_alias_hardening_candidate_2026-05-30.json"
        )
        recovery_candidate = _load(
            "q027_frontier_f0_purged_submission_recovery_activation_2026-05-30.json"
        )
        candidate = _load(
            "q027_frontier_f0_compute_snapshot_activation_2026-05-30.json"
        )
        f1_candidate = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        self.assertEqual(
            candidate["predecessor"]["control_plane_version"],
            recovery_candidate["active_successor"]["control_plane_version"],
        )
        self.assertEqual(
            f1_candidate["initial_live_predecessor"]["control_plane_version"],
            candidate["active_successor"]["control_plane_version"],
        )
        lifecycle = storage["installed_control_plane_lifecycle"]
        if lifecycle == "live_active_generation_successor_staged_not_installed":
            prior_active = f1_candidate.get(
                "active_policy_transition",
                f1_candidate.get("prior_active_policy_transition"),
            )
            expected_active_version = (
                prior_active["control_plane_version"]
                if prior_active is not None
                else candidate["active_successor"]["control_plane_version"]
            )
            self.assertEqual(
                storage["installed_control_plane_version"],
                expected_active_version,
            )
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                f1_candidate["staged_successor"]["control_plane_version"],
            )
            self.assertNotEqual(
                storage["installed_control_plane_version"],
                storage["staged_control_plane_candidate_version"],
            )
            self.assertEqual(
                f1_candidate["paired_install_transition"]["status"],
                "pending_clean_commit_and_paired_install",
            )
            self.assertEqual(
                f1_candidate["staged_successor"]["paired_install_status"],
                f1_candidate["paired_install_transition"]["status"],
            )
        elif lifecycle == "paired_installed_reviewed_generation":
            self.assertEqual(
                storage["installed_control_plane_version"],
                storage["staged_control_plane_candidate_version"],
            )
            paired = f1_candidate["paired_install_transition"]
            self.assertEqual(
                storage["installed_control_plane_version"],
                paired["control_plane_version"],
            )
            self.assertEqual(
                storage["ledger_genesis"]["event_sha256"],
                paired["genesis_event_sha256"],
            )
            self.assertEqual(
                storage["ledger_genesis"]["mirror_ack_sha256"],
                paired["genesis_mirror_ack_sha256"],
            )
            self.assertEqual(
                paired["orion_ledger_records"],
                paired["project_home_ledger_records"],
            )
            self.assertEqual(paired["active_reservations"], 0)
            active_transition = f1_candidate.get("active_policy_transition")
            if active_transition is None:
                self.assertEqual(
                    f1_candidate["staged_successor"]["active_policy_promotion_status"],
                    "pending",
                )
            else:
                self.assertEqual(active_transition["status"], "pass")
                self.assertEqual(
                    active_transition["control_plane_version"],
                    storage["installed_control_plane_version"],
                )
            science_freeze = policy["science_submission_freeze"]
            self.assertEqual(science_freeze["status"], "authorized")
            clean_candidate = candidate["inherited_clean_candidate_freeze"]
            clean_candidate_transition = source_alias_candidate[
                "clean_candidate_policy_transition"
            ]
            self.assertEqual(
                clean_candidate_transition["policy_sha256"],
                "e8909bf541d1c69d4d19c495ebfe983934291de37d6d1423ae4cbfca2e4bb155",
            )
            self.assertEqual(
                clean_candidate_transition["science_submission_freeze_status"],
                science_freeze["status"],
            )
            for key in ["manifest_path", "manifest_sha256"]:
                self.assertEqual(science_freeze[key], clean_candidate_transition[key])
                self.assertEqual(clean_candidate_transition[key], clean_candidate[key])
            self.assertEqual(clean_candidate_transition["orion_ledger_records"], 22)
            self.assertEqual(
                clean_candidate_transition["orion_ledger_records"],
                clean_candidate_transition["project_home_ledger_records"],
            )
            self.assertEqual(
                clean_candidate_transition["orion_ledger_records"],
                clean_candidate_transition["orion_receipt_records"],
            )

            self.assertEqual(clean_candidate_transition["active_reservations"], 0)
            admission_activation = _load(
                "q027_frontier_f0_clean_candidate_admission_policy_activation_2026-05-30.json"
            )
            self.assertEqual(
                admission_activation["control_plane_version"],
                recovery_candidate["predecessor"]["control_plane_version"],
            )
            self.assertEqual(
                admission_activation["policy_sha256"],
                recovery_candidate["predecessor"]["policy_sha256"],
            )
            self.assertEqual(
                admission_activation["science_submission_freeze"],
                {
                    key: science_freeze[key]
                    for key in ["status", "manifest_path", "manifest_sha256"]
                } | {
                    "executable_sha256": clean_candidate["executable_sha256"],
                },
            )
            self.assertEqual(
                admission_activation["frontier_admission_smoke"]["status"],
                "authorized_f0_parser_contract_only",
            )
            self.assertEqual(
                policy["frontier_admission_smoke"], {"status": "closed_after_pass"}
            )
            admission_ledger = admission_activation["ledger_validation"]
            self.assertEqual(admission_ledger["orion_ledger_records"], 22)
            self.assertEqual(
                admission_ledger["orion_ledger_records"],
                admission_ledger["project_home_ledger_records"],
            )
            self.assertEqual(
                admission_ledger["orion_ledger_records"],
                admission_ledger["orion_receipt_records"],
            )
            self.assertEqual(admission_ledger["active_reservations"], 0)
            self.assertEqual(admission_ledger["pending_submission_marker"], "absent")
            recovery_transition = candidate["active_policy_transition"]
            self.assertEqual(
                recovery_transition["control_plane_version"],
                candidate["active_successor"]["control_plane_version"],
            )
            self.assertEqual(
                recovery_transition["policy_sha256"],
                "cfe6610a4f7f54b153f2e20b30397d0417a2634b32997117f6d14a2eb4b771cc",
            )
            self.assertEqual(recovery_transition["orion_ledger_records"], 31)
            self.assertEqual(
                recovery_transition["orion_ledger_records"],
                recovery_transition["project_home_ledger_records"],
            )
            self.assertEqual(
                recovery_transition["orion_ledger_records"],
                recovery_transition["orion_receipt_records"],
            )
            self.assertEqual(recovery_transition["active_reservations"], 0)
            self.assertEqual(
                recovery_transition["pending_submission_marker"], "absent"
            )
        else:
            self.fail(f"Unknown installed-control-plane lifecycle: {lifecycle}")
        self.assertEqual(
            long_term["status"],
            "user_selected_orion_only_with_documented_durability_risk",
        )

    def test_registered_science_launch_contract_sidecars_match_policy(self) -> None:
        policy = _load("storage_policy.json")
        authorizations = {
            record["authorization_id"]: record
            for record in policy["registered_science_slices"]
        }
        sidecars = {
            "f1-clean-gyro-mpich-stderr-v3": "frontier_f1_clean_gyro_launch_contract.json",
            "f1-clean-paper-coupling-mpich-stderr-v2": (
                "frontier_f1_clean_paper_coupling_launch_contract.json"
            ),
        }
        self.assertEqual(set(authorizations), set(sidecars))
        for authorization_id, filename in sidecars.items():
            with self.subTest(authorization_id=authorization_id):
                contract = _load(filename)
                validate_launch_contract(contract)
                self.assertEqual(
                    launch_contract_sha256(contract),
                    authorizations[authorization_id]["launch_contract_sha256"],
                )

    def test_registered_science_staged_bindings_recompute_from_exact_files(self) -> None:
        policy = _load("storage_policy.json")
        storage = policy["olcf_side_storage"]
        staged_version = inventory_digest(
            [
                {"path": name, "sha256": _sha256(CONTROL_PLANE_DIR / name)}
                for name in CONTROL_PLANE_FILES
            ]
        )
        self.assertEqual(
            staged_version, storage["staged_control_plane_candidate_version"]
        )
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        self.assertEqual(
            staged_version, successor["staged_successor"]["control_plane_version"]
        )
        clean_manifest_path = Path(policy["science_submission_freeze"]["manifest_path"])
        clean_manifest = json.loads(clean_manifest_path.read_text(encoding="utf-8"))
        executable_path = Path(clean_manifest["build"]["executable_path"])
        binding_paths = {
            "f1-clean-gyro-mpich-stderr-v3": {
                "job_script_sha256": (
                    REPO_ROOT
                    / "tst/publication/frontier_f1_structured_gpu_relativistic_gyro_job.sh"
                ),
                "input_deck_sha256": (
                    REPO_ROOT / "inputs/tests/pic_relativistic_gyro_paper.athinput"
                ),
                "analysis_script_sha256": [
                    REPO_ROOT
                    / "tst/publication/frontier_f1_gpu_relativistic_gyro_analysis.py",
                    REPO_ROOT / "tst/publication/frontier_f1_structured_artifacts.py",
                ],
            },
            "f1-clean-paper-coupling-mpich-stderr-v2": {
                "job_script_sha256": (
                    REPO_ROOT
                    / "tst/publication/frontier_f1_structured_gpu_paper_coupling_job.sh"
                ),
                "input_deck_sha256": (
                    REPO_ROOT / "inputs/tests/pic_paper_coupling_conservation.athinput"
                ),
                "analysis_script_sha256": [
                    REPO_ROOT
                    / "tst/publication/frontier_f1_gpu_paper_coupling_analysis.py",
                    REPO_ROOT / "tst/publication/frontier_f1_structured_artifacts.py",
                ],
            },
        }
        environment_path = CONTROL_PLANE_DIR / "frontier_pic_environment.sh"
        successor_records = successor["registered_science_slices"]
        successors = {
            record["authorization_id"]: record
            for record in successor_records
        }
        self.assertEqual(len(successor_records), len(successors))
        self.assertEqual(set(successors), set(binding_paths))
        self.assertEqual(
            {
                record["authorization_id"]
                for record in policy["registered_science_slices"]
            },
            set(binding_paths),
        )
        for authorization in policy["registered_science_slices"]:
            authorization_id = authorization["authorization_id"]
            with self.subTest(authorization_id=authorization_id):
                paths = binding_paths[authorization_id]
                for key in (
                    "campaign",
                    "maximum_nodes",
                    "maximum_walltime_seconds",
                    "maximum_attempts",
                ):
                    self.assertEqual(
                        successors[authorization_id][key], authorization[key]
                    )
                self.assertEqual(
                    authorization["job_script_sha256"],
                    _sha256(paths["job_script_sha256"]),
                )
                self.assertEqual(
                    authorization["input_deck_sha256"],
                    _sha256(paths["input_deck_sha256"]),
                )
                self.assertEqual(
                    authorization["environment_profile_sha256"],
                    _sha256(environment_path),
                )
                analysis_sha256 = [
                    _sha256(path) for path in paths["analysis_script_sha256"]
                ]
                self.assertEqual(
                    authorization["analysis_script_sha256"], analysis_sha256
                )
                self.assertEqual(
                    successors[authorization_id]["analysis_script_sha256"],
                    analysis_sha256[0],
                )
                self.assertEqual(
                    successors[authorization_id]["analysis_support_sha256"],
                    analysis_sha256[1],
                )
                self.assertEqual(
                    authorization["clean_candidate_manifest_sha256"],
                    _sha256(clean_manifest_path),
                )
                self.assertEqual(
                    authorization["executable_sha256"], _sha256(executable_path)
                )

    def test_reviewed_mpich_stderr_fixture_matches_failed_attempt_provenance(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        superseded = successor["superseded_registered_execution"]
        provenance_path = (
            REPO_ROOT
            / superseded["reviewed_stderr_transcript_provenance"]
        )
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        failed = provenance["failed_attempt"]
        for key in (
            "job_id",
            "submission_id",
            "reservation_id",
            "registered_science_authorization_id",
            "run_artifact_dir",
        ):
            self.assertEqual(failed[key], superseded[key])
        manifest_path = Path(failed["pre_submit_manifest_path"])
        artifact_root = Path(failed["run_artifact_dir"])
        inventory_path = Path(failed["artifact_inventory_path"])
        self.assertEqual(_sha256(manifest_path), failed["pre_submit_manifest_sha256"])
        self.assertEqual(_sha256(inventory_path), failed["artifact_inventory_sha256"])
        self.assertEqual(f"{artifact_root.stat().st_mode & 0o777:04o}", failed["run_artifact_dir_mode"])
        self.assertEqual(f"{inventory_path.stat().st_mode & 0o777:04o}", failed["artifact_inventory_mode"])
        fixture = provenance["local_review_fixture"]
        fixture_path = REPO_ROOT / fixture["path"]
        decoded = base64.b64decode(
            fixture_path.read_bytes().replace(b"\n", b""),
            validate=True,
        )
        stderr_entry = failed["stderr_inventory_entry"]
        inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
        inventory_records = {
            record["path"]: record for record in inventory["files"]
        }
        stderr_path = artifact_root / stderr_entry["path"]
        self.assertEqual(_sha256(fixture_path), fixture["encoded_file_sha256"])
        self.assertEqual(hashlib.sha256(decoded).hexdigest(), fixture["decoded_sha256"])
        self.assertEqual(len(decoded), fixture["decoded_size"])
        self.assertEqual(fixture["decoded_sha256"], stderr_entry["sha256"])
        self.assertEqual(fixture["decoded_size"], stderr_entry["size"])
        self.assertEqual(inventory_records[stderr_entry["path"]], stderr_entry)
        self.assertEqual(stderr_path.read_bytes(), decoded)
        self.assertEqual(f"{stderr_path.stat().st_mode & 0o777:04o}", failed["stderr_mode"])
        self.assertTrue(fixture["verified_byte_identical_to_live_immutable_stderr"])
        binding = provenance["reviewed_validator_binding"]
        for key, relative in {
            "support_module_sha256": "tst/publication/frontier_f1_structured_artifacts.py",
            "gyro_analyzer_sha256": "tst/publication/frontier_f1_gpu_relativistic_gyro_analysis.py",
            "paper_coupling_analyzer_sha256": "tst/publication/frontier_f1_gpu_paper_coupling_analysis.py",
        }.items():
            self.assertEqual(binding[key], _sha256(REPO_ROOT / relative))
        self.assertEqual(
            provenance["disposition"],
            "pass_historical_transcript_bound_to_local_fixture_retry_requires_separate_v3_policy_activation",
        )

    def test_rejected_pre_reservation_manifest_chronology_is_bound(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        chronology = successor["rejected_pre_reservation_manifest_chronology"]
        fixture = _load(
            "q027_frontier_f1_rejected_pre_reservation_manifest_fixture_2026-05-30.json"
        )
        self.assertEqual(
            successor["rejected_pre_reservation_manifest_fixture"],
            "tst/publication/readiness/"
            "q027_frontier_f1_rejected_pre_reservation_manifest_fixture_2026-05-30.json",
        )
        self.assertEqual(chronology, fixture["chronology"])
        self.assertEqual(
            set(chronology),
            {
                "status",
                "submission_id",
                "registered_science_authorization_id",
                "control_plane_version",
                "manifest_path",
                "manifest_sha256",
                "manifest_analysis_support_sha256",
                "active_policy_analysis_support_sha256",
                "rejection",
                "pre_manifest_attestation_sha256",
                "pre_submit_wrapper_attestation_sha256",
                "reservation_attachments",
                "pending_submission_marker",
                "ledger_uuid_hits",
                "orion_ledger_records",
                "orion_receipt_records",
                "project_home_ledger_records",
                "active_reservations",
            },
        )
        self.assertEqual(
            chronology["status"],
            "fail_closed_before_reservation_intent_and_scheduler_submission",
        )
        self.assertEqual(
            chronology["active_policy_analysis_support_sha256"],
            "c5c77b6a952ed319498f08c91b9adc40101c090f91dc57d798dafdc455a38c01",
        )
        self.assertNotEqual(
            chronology["manifest_analysis_support_sha256"],
            chronology["active_policy_analysis_support_sha256"],
        )
        binding = fixture["manifest_binding"]
        self.assertEqual(binding["mode"], "0444")
        for key in (
            "submission_id",
            "control_plane_version",
            "registered_science_authorization_id",
            "manifest_sha256",
        ):
            self.assertEqual(binding[key], chronology[key])
        self.assertEqual(binding["analysis_support_role"], "analysis-script-001")
        self.assertEqual(
            binding["analysis_support_sha256"],
            chronology["manifest_analysis_support_sha256"],
        )
        self.assertEqual(chronology["reservation_attachments"], "absent")
        self.assertEqual(chronology["pending_submission_marker"], "absent")
        self.assertEqual(chronology["ledger_uuid_hits"], 0)
        self.assertEqual(chronology["orion_ledger_records"], 34)
        self.assertEqual(chronology["orion_receipt_records"], 34)
        self.assertEqual(chronology["project_home_ledger_records"], 34)
        self.assertEqual(chronology["active_reservations"], 0)
        template = _load(
            "q027_frontier_registered_science_same_account_isolation_attestation_template_2026-05-30.json"
        )
        attestations = fixture["attestations"]
        self.assertEqual(
            [record["contents"]["phase"] for record in attestations],
            ["pre_manifest", "pre_submit_wrapper"],
        )
        self.assertLess(
            attestations[0]["contents"]["recorded_utc"],
            attestations[1]["contents"]["recorded_utc"],
        )
        digest_keys = (
            "pre_manifest_attestation_sha256",
            "pre_submit_wrapper_attestation_sha256",
        )
        for record, digest_key in zip(attestations, digest_keys):
            contents = record["contents"]
            rendered = (json.dumps(contents, indent=2, sort_keys=True) + "\n").encode()
            self.assertEqual(
                hashlib.sha256(rendered).hexdigest(),
                chronology[digest_key],
            )
            self.assertEqual(record["attestation_sha256"], chronology[digest_key])
            self.assertEqual(contents["control_plane_version"], chronology["control_plane_version"])
            self.assertEqual(
                contents["registered_science_authorization_id"],
                chronology["registered_science_authorization_id"],
            )
            self.assertEqual(contents["operator_statement"], template["operator_statement"])
            self.assertEqual(contents["pending_submission_marker"]["value"], "absent")
            self.assertEqual(
                sorted(contents["mirrored_ledger_line_counts"]["counts"].values()),
                [34, 34, 34],
            )

    def test_rejected_pre_reservation_manifest_live_preflight_is_explicit(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        chronology = successor["rejected_pre_reservation_manifest_chronology"]
        preflight = successor["required_live_preflight"]
        self.assertEqual(
            preflight["status"],
            "required_immediately_before_policy_promotion_and_each_registered_science_submission_boundary",
        )
        queue_chronology = successor[
            "rejected_pre_reservation_operator_queue_format_chronology"
        ]
        self.assertEqual(
            preflight["historical_submission_ids_must_remain_absent_from_live_ledgers"],
            [chronology["submission_id"], queue_chronology["submission_id"]],
        )
        self.assertEqual(
            preflight["checks"],
            [
                "current Orion ledger, Orion mirror-receipt and Project Home mirror chains are coherent",
                "current Orion pending_submission.json marker is absent",
                "current mirrored ledger state has zero active reservations",
                "each historical rejected pre-reservation submission UUID has zero hits in current Orion JSONL, Orion CSV, Orion receipts and Project Home mirror streams",
                "same-account process and scheduler snapshots are reviewed for the current boundary",
            ],
        )
        self.assertIn("Legitimate later reservations", preflight["evidence_rule"])

    def test_rejected_operator_queue_format_manifest_chronology_is_bound(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        chronology = successor[
            "rejected_pre_reservation_operator_queue_format_chronology"
        ]
        fixture = _load(
            "q027_frontier_f1_rejected_operator_queue_format_manifest_fixture_2026-05-30.json"
        )
        self.assertEqual(
            successor["rejected_pre_reservation_operator_queue_format_fixture"],
            "tst/publication/readiness/"
            "q027_frontier_f1_rejected_operator_queue_format_manifest_fixture_2026-05-30.json",
        )
        self.assertEqual(chronology, fixture["chronology"])
        manifest_path = Path(chronology["manifest_path"])
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        binding = fixture["manifest_binding"]
        self.assertEqual(f"{manifest_path.stat().st_mode & 0o777:04o}", binding["mode"])
        self.assertEqual(_sha256(manifest_path), chronology["manifest_sha256"])
        self.assertEqual(manifest["submission_id"], chronology["submission_id"])
        self.assertEqual(
            manifest["queue_snapshot_sha256"], binding["queue_snapshot_sha256"]
        )
        self.assertEqual(
            binding["queue_snapshot_format"], "%i|%a|%P|%q|%T|%j|%k"
        )
        self.assertEqual(
            binding["validator_queue_snapshot_format"], "%i|%P|%q|%T|%j|%k"
        )
        self.assertNotEqual(
            binding["queue_snapshot_sha256"],
            binding["validator_queue_snapshot_sha256"],
        )
        self.assertEqual(chronology["reservation_attachments"], "absent")
        for name in ("reservation_id.txt", "manifest_sha256.txt"):
            self.assertFalse((manifest_path.parent / name).exists())
        for attestation, digest_key in zip(
            fixture["attestations"],
            (
                "pre_manifest_attestation_sha256",
                "pre_submit_wrapper_attestation_sha256",
            ),
        ):
            path = Path(attestation["path"])
            self.assertEqual(attestation["sha256"], chronology[digest_key])
            self.assertEqual(_sha256(path), chronology[digest_key])
            contents = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(contents["phase"], attestation["phase"])
            self.assertEqual(
                contents["queue_snapshot"]["sha256"],
                binding["queue_snapshot_sha256"],
            )
        live_surfaces = (
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl",
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.csv",
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/mirror_receipts.jsonl",
            "/ccs/proj/ast207/proj-shared/PIC/ledger/node_hours.jsonl",
        )
        for surface in live_surfaces:
            contents = Path(surface).read_text()
            for submission_id in successor["required_live_preflight"][
                "historical_submission_ids_must_remain_absent_from_live_ledgers"
            ]:
                self.assertNotIn(submission_id, contents)

    def test_reconciled_gyro_v2_analysis_rejection_chronology_is_bound(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        chronology = successor["superseded_gyro_v2_registered_execution"]
        fixture = _load(
            "q027_frontier_f1_gyro_v2_analysis_rejection_fixture_2026-05-30.json"
        )
        self.assertEqual(
            successor["superseded_gyro_v2_registered_execution_fixture"],
            "tst/publication/readiness/"
            "q027_frontier_f1_gyro_v2_analysis_rejection_fixture_2026-05-30.json",
        )
        self.assertEqual(chronology, fixture["chronology"])
        inventory_path = Path(fixture["artifact_inventory_path"])
        self.assertEqual(
            f"{inventory_path.stat().st_mode & 0o777:04o}",
            fixture["artifact_inventory_mode"],
        )
        self.assertEqual(_sha256(inventory_path), chronology["artifact_inventory_sha256"])
        self.assertEqual(
            json.loads(inventory_path.read_text(encoding="utf-8")),
            fixture["artifact_inventory"],
        )
        output_paths = [
            record["path"]
            for record in fixture["artifact_inventory"]["files"]
            if record["path"].startswith("output/")
        ]
        self.assertEqual(len(output_paths), 21)
        self.assertEqual(
            output_paths,
            sorted(output_paths),
        )
        artifact_dir = Path(chronology["run_artifact_dir"])
        self.assertEqual(list((artifact_dir / "analysis").iterdir()), [])
        ledger_events = [
            json.loads(line)
            for line in Path(
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl"
            ).read_text(encoding="utf-8").splitlines()
        ]
        self.assertIn(fixture["terminal_reconciliation_event"], ledger_events)
        self.assertEqual(chronology["analysis_result"], "absent")
        self.assertEqual(chronology["offline_analysis_receipt"], "absent")

    def test_registered_science_same_account_isolation_attestation_template(self) -> None:
        template = _load(
            "q027_frontier_registered_science_same_account_isolation_attestation_template_2026-05-30.json"
        )
        self.assertEqual(
            template["archive_root"],
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/operator_attestations",
        )
        self.assertEqual(
            template["status"],
            "required_before_each_registered_science_submission",
        )
        self.assertIn("same_account_process_snapshot", template["required_fields"])
        self.assertIn("operator_statement", template["required_fields"])
        self.assertIn(
            "test -e \"${PIC_ROOT}/ledger/pending_submission.json\" "
            "&& printf 'present\\n' || printf 'absent\\n'",
            template["required_snapshot_commands"],
        )

    def test_claim_classes_and_required_extensions(self) -> None:
        registry = _load("claims_registry.json")
        self.assertEqual(registry["default_reviewer"], "pending external review")
        claims = registry["claims"]
        ids = {claim["claim_id"] for claim in claims}
        self.assertEqual(len(ids), len(claims))
        self.assertTrue(REQUIRED_EXTENSION_CLAIMS.issubset(ids))
        for claim in claims:
            self.assertIn(claim["claim_class"], CLAIM_CLASSES)
            self.assertEqual(claim["disposition"], "open")
            self.assertTrue(claim["required_gates"])
            if claim["claim_id"] in REQUIRED_EXTENSION_CLAIMS:
                self.assertTrue(claim["authorized_extension_required"])

    def test_initial_findings_are_unique_and_have_valid_status(self) -> None:
        registry = _load("findings_registry.json")
        findings = registry["findings"]
        ids = {finding["finding_id"] for finding in findings}
        self.assertEqual(len(ids), len(findings))
        self.assertEqual(ids, {f"PIC-P0-00{i}" for i in range(1, 7)} | {
            f"PIC-P1-{i:03d}" for i in range(1, 11)
        })
        expected_verifying = {
            "PIC-P0-001", "PIC-P0-002", "PIC-P0-003", "PIC-P0-004",
            "PIC-P0-005", "PIC-P0-006", "PIC-P1-002", "PIC-P1-004",
            "PIC-P1-008", "PIC-P1-009", "PIC-P1-010",
        }
        for finding in findings:
            expected = "verifying" if finding["finding_id"] in expected_verifying else "open"
            self.assertEqual(finding["status"], expected)
            self.assertRegex(finding["verification_gate"], r"^Q-\d{3}$")

    def test_paper_source_checksum_is_frozen(self) -> None:
        inventory = _load("external_artifacts.json")
        artifacts = {
            artifact["artifact_id"]: artifact for artifact in inventory["artifacts"]
        }
        tex = artifacts["SUN_BAI_2023_ARXIV_V1_TEX"]
        self.assertEqual(tex["sha256"], _sha256(PAPER_TEX))
        entity = artifacts["ENTITY_TOOLKIT_REPLACEMENT_CANDIDATE"]
        self.assertEqual(
            entity["git_commit"],
            "512998c471bf3fdec292cb4a64150c4f0aeea539",
        )
        self.assertEqual(entity["historical_unavailable_commit"], "a59065fc")
        snapshot = _load("entity_snapshot.json")
        self.assertEqual(entity["snapshot_manifest"],
                         "tst/publication/readiness/entity_snapshot.json")
        self.assertEqual(snapshot["git_commit"], entity["git_commit"])
        self.assertEqual(snapshot["git_tree"], entity["git_tree"])
        self.assertEqual(snapshot["git_archive_sha256"],
                         entity["git_archive_sha256"])
        self.assertEqual(snapshot["worktree_status"], "clean")
        self.assertTrue(snapshot["files"])
        for source_file in snapshot["files"]:
            self.assertRegex(source_file["sha256"], r"^[0-9a-f]{64}$")

    def test_validation_manifest_schema_has_release_minimum(self) -> None:
        schema = _validation_manifest_schema()
        required = set(schema["required"])
        self.assertTrue(
            {
                "claim_ids",
                "evidence_class",
                "physical_mode",
                "git",
                "executable",
                "oracle",
                "metrics",
                "resources",
                "artifacts",
                "review",
            }.issubset(required)
        )
        classes = set(schema["properties"]["evidence_class"]["enum"])
        self.assertEqual(classes, CLAIM_CLASSES)

    def test_validation_manifest_schema_accepts_reviewable_scaffolds(self) -> None:
        schema = _validation_manifest_schema()
        pending = _minimum_validation_manifest()
        validate_schema(pending, schema)

        qualified = copy.deepcopy(pending)
        qualified["review"] = {
            "reviewer": "Named External Reviewer",
            "disposition": "qualified",
        }
        validate_schema(qualified, schema)

    def test_validation_manifest_schema_rejects_incomplete_evidence(self) -> None:
        schema = _validation_manifest_schema()
        cases = [
            ("empty metrics", ("metrics",), []),
            ("empty metric record", ("metrics",), [{}]),
            ("empty artifacts", ("artifacts",), []),
            ("empty tolerances", ("oracle", "tolerances"), {}),
            ("blank manifest ID", ("manifest_id",), " "),
            ("blank claim ID", ("claim_ids", 0), " "),
            ("blank test ID", ("test_id",), " "),
            ("blank physical mode", ("physical_mode",), " "),
            ("blank executable path", ("executable", "path"), " "),
            ("blank CMake cache path", ("executable", "cmake_cache", "path"), " "),
            ("blank modules path", ("executable", "modules", "path"), " "),
            (
                "blank environment allowlist path",
                ("executable", "environment_allowlist", "path"),
                " ",
            ),
            ("blank oracle kind", ("oracle", "kind"), " "),
            ("blank oracle reference", ("oracle", "reference"), " "),
            ("blank platform", ("resources", "platform"), " "),
            ("blank artifact root", ("resources", "artifact_root"), " "),
            ("blank artifact path", ("artifacts", 0, "path"), " "),
            ("blank reviewer", ("review", "reviewer"), " "),
            (
                "pending reviewer qualified disposition",
                ("review", "disposition"),
                "qualified",
            ),
        ]
        for label, path, replacement in cases:
            with self.subTest(label=label):
                manifest = _minimum_validation_manifest()
                _replace_nested(manifest, path, replacement)
                with self.assertRaises(ValueError):
                    validate_schema(manifest, schema)


if __name__ == "__main__":
    unittest.main()
