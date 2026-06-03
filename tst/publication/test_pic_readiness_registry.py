#!/usr/bin/env python3
"""Regression tests for source-controlled PIC readiness registries."""

from __future__ import annotations

import base64
from contextlib import contextmanager
import copy
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import unittest
from collections.abc import Iterator

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tst" / "publication" / "frontier_control_plane"))

from control_plane_common import CONTROL_PLANE_FILES
from control_plane_common import PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS
from control_plane_common import inventory_digest
from control_plane_common import launch_contract_sha256
from control_plane_common import validate_launch_contract
from ledger import record_sha256
from ledger import incomplete_manual_accounting_marker_paths
from ledger import validate_mirrored_state
from tst.publication import q011_section54_pressure_pilot_execution as pressure_execution
from tst.publication.pic_qualification_manifest import SCHEMA_PATH
from tst.publication.pic_qualification_manifest import validate_qualification_manifest
from tst.publication.pic_qualification_manifest import validate_schema


READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
SCHEMA_DIR = READINESS_DIR / "schemas"
CONTROL_PLANE_DIR = REPO_ROOT / "tst" / "publication" / "frontier_control_plane"
VALIDATION_MANIFEST_SCHEMA = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
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


def _descriptor_bytes(descriptor: int) -> bytes:
    os.lseek(descriptor, 0, os.SEEK_SET)
    chunks = []
    while chunk := os.read(descriptor, 1024 * 1024):
        chunks.append(chunk)
    return b"".join(chunks)


def _require_regular_namespace_identity(path: Path, identity: tuple[int, int]) -> None:
    observed = os.stat(path, follow_symlinks=False)
    if not stat.S_ISREG(observed.st_mode) or (observed.st_dev, observed.st_ino) != identity:
        raise ValueError(f"pinned regular file namespace changed: {path}")


@contextmanager
def _pinned_regular_bytes(path: Path) -> Iterator[bytes]:
    descriptor = os.open(
        path,
        os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ValueError(f"pinned path is not a regular file: {path}")
        identity = (before.st_dev, before.st_ino)
        data = _descriptor_bytes(descriptor)
        after = os.fstat(descriptor)
        if (
            identity != (after.st_dev, after.st_ino)
            or before.st_size != after.st_size
            or len(data) != after.st_size
        ):
            raise ValueError(f"pinned regular file changed during read: {path}")
        _require_regular_namespace_identity(path, identity)
        yield data
        if _descriptor_bytes(descriptor) != data:
            raise ValueError(f"pinned regular file bytes changed during validation: {path}")
        _require_regular_namespace_identity(path, identity)
    finally:
        os.close(descriptor)


def _git_blob_sha256(commit: str, relative_path: str) -> str:
    contents = subprocess.check_output(
        ["git", "show", f"{commit}:{relative_path}"],
        cwd=REPO_ROOT,
    )
    return hashlib.sha256(contents).hexdigest()


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
        manual_accounting_activation = _load(
            "q027_manual_frontier_accounting_activation_2026-05-30.json"
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
        science_freeze = policy["science_submission_freeze"]
        if science_freeze == {"status": "pending_clean_candidate_freeze"}:
            phase0_successor = _load(
                "phase0_curated_candidate_successor_v14_2026-06-01.json"
            )
            self.assertEqual(lifecycle, "paired_installed_reviewed_generation")
            self.assertEqual(
                storage["installed_control_plane_version"],
                phase0_successor["successor_source_control_plane_version"],
            )
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                phase0_successor["successor_source_control_plane_version"],
            )
            self.assertEqual(policy["registered_science_slices"], [])
            return
        d720_promotion = _load(
            "phase0_clean_candidate_freeze_and_policy_promotion_"
            "successor_v2_2026-06-02.json"
        )
        c83e_promotion = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v4_2026-06-02.json"
        )
        strict_q011_promotion = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v6_2026-06-02.json"
        )
        q011_pressure_promotion = _load(
            "phase0_curated_candidate_successor_v20_2026-06-02.json"
        )
        if (
            science_freeze
            == q011_pressure_promotion["policy_promotion"][
                "science_submission_freeze"
            ]
        ):
            self.assertEqual(lifecycle, "paired_installed_reviewed_generation")
            self.assertEqual(
                storage["installed_control_plane_version"],
                q011_pressure_promotion["live_paired_control_plane_version"],
            )
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                q011_pressure_promotion["live_paired_control_plane_version"],
            )
            self.assertEqual(
                [
                    record["authorization_id"]
                    for record in policy["registered_science_slices"]
                ],
                q011_pressure_promotion["policy_promotion"][
                    "registered_science_authorization_ids"
                ],
            )
            return
        if (
            science_freeze
            == d720_promotion["active_policy_promotion"]["science_submission_freeze"]
        ):
            expected_control_plane_version = d720_promotion["control_plane_version"]
            if (
                storage["installed_control_plane_version"]
                == c83e_promotion["control_plane_version"]
            ):
                expected_control_plane_version = c83e_promotion[
                    "control_plane_version"
                ]
                self.assertEqual(
                    science_freeze,
                    c83e_promotion["active_policy_promotion"][
                        "science_submission_freeze"
                    ],
                )
            if (
                storage["installed_control_plane_version"]
                == strict_q011_promotion["control_plane_version"]
            ):
                expected_control_plane_version = strict_q011_promotion[
                    "control_plane_version"
                ]
                self.assertEqual(
                    science_freeze,
                    strict_q011_promotion["active_policy_promotion"][
                        "science_submission_freeze"
                    ],
                )
            self.assertEqual(lifecycle, "paired_installed_reviewed_generation")
            self.assertEqual(
                storage["installed_control_plane_version"],
                expected_control_plane_version,
            )
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                expected_control_plane_version,
            )
            self.assertEqual(policy["registered_science_slices"], [])
            self.assertEqual(
                d720_promotion["active_policy_promotion"][
                    "frontier_launch_authorization"
                ],
                "none_no_registered_science_slices",
            )
            return
        replay = _load(
            "phase0_registered_prerequisite_replay_policy_promotion_2026-06-01.json"
        )
        if (
            storage["installed_control_plane_version"]
            == replay["control_plane_version"]
        ):
            accounting = _load(
                "phase0_scheduler_accounting_controller_successor_2026-06-01.json"
            )
            freeze = _load(
                "phase0_clean_candidate_freeze_and_policy_promotion_2026-06-01.json"
            )
            self.assertEqual(lifecycle, "paired_installed_reviewed_generation")
            self.assertEqual(
                storage["installed_control_plane_version"],
                storage["staged_control_plane_candidate_version"],
            )
            self.assertEqual(
                replay["predecessor_record"],
                "tst/publication/readiness/"
                "phase0_scheduler_accounting_controller_successor_2026-06-01.json",
            )
            self.assertEqual(
                accounting["predecessor_record"],
                "tst/publication/readiness/"
                "phase0_clean_candidate_freeze_and_policy_promotion_2026-06-01.json",
            )
            self.assertEqual(
                science_freeze["manifest_path"],
                freeze["clean_candidate"]["manifest_path"],
            )
            self.assertEqual(
                science_freeze["manifest_sha256"],
                freeze["clean_candidate"]["manifest_sha256"],
            )
            self.assertEqual(
                {
                    record["authorization_id"]: {
                        "campaign": record["campaign"],
                        "launch_contract_sha256": record["launch_contract_sha256"],
                    }
                    for record in policy["registered_science_slices"]
                },
                {
                    record["authorization_id"]: {
                        "campaign": record["campaign"],
                        "launch_contract_sha256": record["launch_contract_sha256"],
                    }
                    for record in replay["registered_science_slices"]
                },
            )
            terminal = replay["terminal_mirrored_ledger"]
            self.assertEqual(terminal["orion_ledger_records"], 66)
            self.assertEqual(
                terminal["orion_ledger_records"],
                terminal["project_home_ledger_records"],
            )
            self.assertEqual(
                terminal["orion_ledger_records"],
                terminal["orion_receipt_records"],
            )
            self.assertEqual(terminal["active_reservations"], 0)
            self.assertEqual(terminal["pending_submission_marker"], "absent")
            self.assertEqual(terminal["pending_manual_accounting_marker"], "absent")
            self.assertEqual(
                long_term["status"],
                "user_selected_orion_only_with_documented_durability_risk",
            )
            return
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
            current_transition = manual_accounting_activation[
                "control_plane_transition"
            ]
            self.assertEqual(
                storage["installed_control_plane_version"],
                current_transition["control_plane_version"],
            )
            terminal_ledger = manual_accounting_activation["terminal_ledger"]
            self.assertEqual(terminal_ledger["orion_records"], 56)
            self.assertEqual(
                terminal_ledger["orion_records"],
                terminal_ledger["project_home_records"],
            )
            self.assertEqual(
                terminal_ledger["orion_records"],
                terminal_ledger["mirror_receipts"],
            )
            self.assertEqual(terminal_ledger["currently_reserved_node_hours"], 0.0)
            self.assertEqual(
                manual_accounting_activation["authorization"][
                    "scientific_evidence_eligible"
                ],
                False,
            )
            paired = f1_candidate["paired_install_transition"]
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
            "f1_gpu_relativistic_gyro": "frontier_f1_clean_gyro_launch_contract.json",
            "f1_gpu_paper_coupling": (
                "frontier_f1_clean_paper_coupling_launch_contract.json"
            ),
            "f2_multirank_runtime_metadata": (
                "frontier_f2_multirank_runtime_metadata_launch_contract.json"
            ),
        }
        if not authorizations:
            self.assertIn(
                policy["science_submission_freeze"]["status"],
                ("pending_clean_candidate_freeze", "authorized"),
            )
            self.assertEqual(authorizations, {})
            for campaign, filename in sidecars.items():
                with self.subTest(campaign=campaign):
                    validate_launch_contract(_load(filename))
            return
        authorizations_by_campaign = {
            record["campaign"]: record
            for record in policy["registered_science_slices"]
        }
        q011_sidecars = {
            case.campaign: case.launch_contract_path.name
            for case in pressure_execution.CASES
        }
        if set(authorizations_by_campaign) == set(q011_sidecars):
            for campaign, filename in q011_sidecars.items():
                authorization = authorizations_by_campaign[campaign]
                with self.subTest(campaign=campaign):
                    contract = _load(filename)
                    self.assertEqual(
                        pressure_execution._launch_contract_sha256(contract),
                        authorization["launch_contract_sha256"],
                    )
            return
        self.assertEqual(set(authorizations_by_campaign), set(sidecars))
        for campaign, filename in sidecars.items():
            authorization = authorizations_by_campaign[campaign]
            authorization_id = authorization["authorization_id"]
            with self.subTest(authorization_id=authorization_id):
                contract = _load(filename)
                validate_launch_contract(contract)
                self.assertEqual(
                    launch_contract_sha256(contract),
                    authorization["launch_contract_sha256"],
                )

    def test_q011_pressure_pilot_registered_execution_source_tranche_is_frozen(
        self,
    ) -> None:
        record = _load(
            "q011_section54_pressure_pilot_registered_execution_preregistration_2026-06-02.json"
        )
        boundary = record["execution_boundary"]
        self.assertFalse(boundary["frontier_execution_authorized_by_this_record"])
        self.assertFalse(boundary["scheduler_commands_authorized_by_this_record"])
        self.assertFalse(boundary["storage_policy_mutation_authorized_by_this_record"])
        self.assertEqual(
            record["source_bindings"]["materializer"]["sha256"],
            "57f6300b4d9a7cdf2327dc92137c9dcb040eeb815d5c79f33f01e646703b6ca4",
        )
        self.assertEqual(
            [entry["authorization_id"] for entry in record["launch_matrix"]],
            [
                f"q011-section54-pressure-{case.case_id.replace('_', '-')}-v1"
                for case in pressure_execution.CASES
            ],
        )

    def test_q011_pressure_pilot_registered_execution_retry_successor_is_frozen(
        self,
    ) -> None:
        record = _load(
            "q011_section54_pressure_pilot_registered_execution_retry_"
            "successor_v2_2026-06-02.json"
        )
        self.assertEqual(
            record["predecessor_sha256"],
            _sha256(REPO_ROOT / record["predecessor_record"]),
        )
        chronology = record["failed_attempt_chronology"]
        self.assertEqual(chronology["job_id"], "4754211")
        self.assertEqual(chronology["terminal_state"], "FAILED")
        self.assertEqual(
            chronology["terminal_reconciliation_event_sha256"],
            "fc0082bef800733c433395d48f551559abe083367e5549cc4bffbe8a0ab48bfa",
        )
        boundary = record["execution_boundary"]
        self.assertFalse(boundary["frontier_execution_authorized_by_this_record"])
        self.assertFalse(boundary["scheduler_commands_authorized_by_this_record"])
        self.assertFalse(boundary["storage_policy_mutation_authorized_by_this_record"])
        self.assertEqual(
            record["source_bindings"],
            pressure_execution._historical_v2_preregistration()["source_bindings"],
        )
        status = pressure_execution.historical_v2_source_tranche_status()
        self.assertEqual(status["state"], "historical_consumed_slice_non_authorizing")
        self.assertFalse(status["source_bindings_match_current_checkout"])
        self.assertFalse(status["consumed_slice_reauthorization_allowed"])
        with self.assertRaisesRegex(
            pressure_execution.ContractError,
            "historical v2 registered-execution tranche is consumed",
        ):
            pressure_execution.validate_source_tranche()
        self.assertEqual(len(record["launch_matrix"]), len(pressure_execution.CASES))
        binding_by_case = {
            binding["case_id"]: binding
            for binding in record["source_bindings"]["launch_contracts"]
        }
        for case in pressure_execution.CASES:
            with self.subTest(case_id=case.case_id):
                contract = _load(case.launch_contract_path.name)
                self.assertEqual(
                    contract, pressure_execution.expected_launch_contract(case)
                )
                validate_launch_contract(contract)
                binding = binding_by_case[case.case_id]
                self.assertEqual(binding["file_sha256"], _sha256(case.launch_contract_path))
                self.assertEqual(
                    binding["launch_contract_sha256"],
                    launch_contract_sha256(contract),
                )

    def test_f2_multirank_runtime_metadata_candidate_resolves_source_commit(self) -> None:
        candidate = _load(
            "q027_frontier_f2_multirank_runtime_metadata_candidate_2026-05-30.json"
        )
        commit = candidate["implementation_source_commit"]
        self.assertEqual(
            subprocess.check_output(
                ["git", "cat-file", "-t", commit],
                cwd=REPO_ROOT,
                text=True,
            ).strip(),
            "commit",
        )
        binding = candidate["registered_science_slice"]
        self.assertEqual(
            candidate["staged_policy_sha256"],
            candidate["accepted_v2_execution"]["active_policy_sha256"],
        )
        for digest_key, relative_path in {
            "job_script_sha256":
                "tst/publication/frontier_f2_structured_multirank_runtime_metadata_job.sh",
            "input_deck_sha256": "inputs/tests/pic_parser_contract_guards.athinput",
            "analysis_script_sha256":
                "tst/publication/frontier_f2_multirank_runtime_metadata_analysis.py",
            "analysis_support_sha256":
                "tst/publication/frontier_f1_structured_artifacts.py",
        }.items():
            self.assertEqual(binding[digest_key], _git_blob_sha256(commit, relative_path))
        contract = _load("frontier_f2_multirank_runtime_metadata_launch_contract.json")
        self.assertEqual(
            binding["launch_contract_sha256"],
            launch_contract_sha256(contract),
        )

    def _assert_current_registered_replay_bindings(
        self, policy: dict[str, object], staged_version: str
    ) -> None:
        replay = _load(
            "phase0_registered_prerequisite_replay_policy_promotion_2026-06-01.json"
        )
        accounting = _load(
            "phase0_scheduler_accounting_controller_successor_2026-06-01.json"
        )
        closure = _load(
            "phase0_registered_prerequisite_replay_closure_2026-06-01.json"
        )
        storage = policy["olcf_side_storage"]
        self.assertEqual(staged_version, replay["control_plane_version"])
        self.assertEqual(
            storage["installed_control_plane_version"], staged_version
        )
        self.assertEqual(
            storage["staged_control_plane_candidate_version"], staged_version
        )
        promotion = replay["active_policy_promotion"]
        self.assertEqual(
            promotion["orion_policy_sha256"],
            _sha256(READINESS_DIR / "storage_policy.json"),
        )
        self.assertEqual(
            promotion["orion_policy_sha256"],
            _sha256(Path(promotion["orion_policy_path"])),
        )
        self.assertEqual(
            promotion["project_home_policy_sha256"],
            _sha256(Path(promotion["project_home_policy_path"])),
        )
        self.assertEqual(
            promotion["orion_promotion_sha256"],
            _sha256(
                Path(promotion["orion_policy_path"]).with_name(
                    "active_promotion.json"
                )
            ),
        )
        self.assertEqual(
            promotion["project_home_promotion_sha256"],
            _sha256(
                Path(promotion["project_home_policy_path"]).with_name(
                    "active_promotion.json"
                )
            ),
        )
        clean_manifest_path = Path(
            policy["science_submission_freeze"]["manifest_path"]
        )
        self.assertEqual(
            _sha256(clean_manifest_path),
            replay["clean_candidate"]["manifest_sha256"],
        )
        clean_manifest = json.loads(clean_manifest_path.read_text(encoding="utf-8"))
        executable_path = Path(clean_manifest["build"]["executable_path"])
        self.assertEqual(
            _sha256(executable_path),
            replay["clean_candidate"]["executable_sha256"],
        )
        binding_paths = {
            "f1_gpu_relativistic_gyro": {
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
            "f1_gpu_paper_coupling": {
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
            "f2_multirank_runtime_metadata": {
                "job_script_sha256": (
                    REPO_ROOT
                    / "tst/publication/frontier_f2_structured_multirank_runtime_metadata_job.sh"
                ),
                "input_deck_sha256": (
                    REPO_ROOT / "inputs/tests/pic_parser_contract_guards.athinput"
                ),
                "analysis_script_sha256": [
                    REPO_ROOT
                    / "tst/publication/frontier_f2_multirank_runtime_metadata_analysis.py",
                    REPO_ROOT / "tst/publication/frontier_f1_structured_artifacts.py",
                ],
            },
        }
        expected_metadata = {
            "f1_gpu_relativistic_gyro": {
                "test_id": "pic_relativistic_gyro_paper",
                "evidence_class": "frontier_f1_clean_candidate_gpu_pusher_oracle",
                "physical_mode": "paper_test_particle",
            },
            "f1_gpu_paper_coupling": {
                "test_id": "pic_paper_coupling_conservation",
                "evidence_class": "frontier_f1_clean_candidate_gpu_paper_coupling_oracle",
                "physical_mode": "paper_mhd_pic",
            },
            "f2_multirank_runtime_metadata": {
                "test_id": "pic_parser_contract_guards",
                "evidence_class": "frontier_f2_clean_candidate_multirank_runtime_metadata",
                "physical_mode": "extended_mhd_pic_parser_contract",
            },
        }
        authorizations = {
            record["campaign"]: record
            for record in policy["registered_science_slices"]
        }
        replay_slices = {
            record["campaign"]: record for record in replay["registered_science_slices"]
        }
        self.assertEqual(set(authorizations), set(binding_paths))
        self.assertEqual(set(replay_slices), set(binding_paths))
        environment_path = CONTROL_PLANE_DIR / "frontier_pic_environment.sh"
        for campaign, authorization in authorizations.items():
            with self.subTest(campaign=campaign):
                paths = binding_paths[campaign]
                for key, expected in expected_metadata[campaign].items():
                    self.assertEqual(authorization[key], expected)
                self.assertEqual(authorization["runtime_profile"], "frontier_minimum_supported")
                self.assertEqual(authorization["selected_qos"], "debug")
                self.assertIs(authorization["registered_short_nonproduction"], True)
                self.assertEqual(authorization["maximum_nodes"], 1)
                self.assertEqual(authorization["maximum_walltime_seconds"], 900)
                self.assertEqual(authorization["maximum_attempts"], 1)
                self.assertEqual(
                    authorization["authorization_id"],
                    replay_slices[campaign]["authorization_id"],
                )
                self.assertEqual(
                    authorization["launch_contract_sha256"],
                    replay_slices[campaign]["launch_contract_sha256"],
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
                self.assertEqual(
                    authorization["analysis_script_sha256"],
                    [_sha256(path) for path in paths["analysis_script_sha256"]],
                )
                self.assertEqual(
                    authorization["clean_candidate_manifest_sha256"],
                    _sha256(clean_manifest_path),
                )
                self.assertEqual(
                    authorization["executable_sha256"], _sha256(executable_path)
                )
        executions = {
            record["campaign"]: record for record in closure["registered_replays"]
        }
        self.assertEqual(set(executions), set(binding_paths))
        ledger_records = [
            json.loads(line)
            for line in Path(
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/"
                "node_hours.jsonl"
            )
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        terminal_events = {
            record["job_id"]: record
            for record in ledger_records
            if record.get("event_type") == "reconciliation"
            and record.get("job_id") in {
                execution["job_id"] for execution in executions.values()
            }
        }
        self.assertEqual(len(terminal_events), len(executions))
        for campaign, execution in executions.items():
            with self.subTest(campaign=campaign):
                self.assertEqual(
                    execution["authorization_id"],
                    authorizations[campaign]["authorization_id"],
                )
                manifest_path = Path(execution["manifest_path"])
                artifact_dir = Path(execution["artifact_dir"])
                self.assertEqual(_sha256(manifest_path), execution["manifest_sha256"])
                self.assertEqual(
                    _sha256(artifact_dir / "artifact_inventory.json"),
                    execution["artifact_inventory_sha256"],
                )
                result_path = artifact_dir / "analysis" / "analysis.json"
                self.assertEqual(
                    _sha256(result_path), execution["analysis_result_sha256"]
                )
                result = json.loads(result_path.read_text(encoding="utf-8"))
                self.assertEqual(result["status"], "pass")
                receipt_path = (
                    artifact_dir / "analysis" / "offline_analysis_receipt.json"
                )
                self.assertEqual(
                    _sha256(receipt_path),
                    execution["offline_analysis_receipt_sha256"],
                )
                receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
                self.assertEqual(
                    receipt["artifact_inventory"]["sha256"],
                    execution["artifact_inventory_sha256"],
                )
                self.assertEqual(
                    receipt["analysis_result"]["sha256"],
                    execution["analysis_result_sha256"],
                )
                for attestation_key in (
                    "pre_manifest_attestation",
                    "pre_submit_wrapper_attestation",
                ):
                    attestation = execution[attestation_key]
                    self.assertEqual(
                        _sha256(Path(attestation["path"])),
                        attestation["sha256"],
                    )
                terminal = terminal_events[execution["job_id"]]
                self.assertEqual(terminal["state"], "COMPLETED")
                self.assertEqual(
                    terminal["event_sha256"],
                    execution["terminal_ledger_event_sha256"],
                )
                self.assertEqual(
                    terminal["reservation_id"], execution["reservation_id"]
                )
                self.assertEqual(
                    terminal["submission_id"], execution["submission_id"]
                )
        ledger = closure["terminal_mirrored_ledger"]
        self.assertEqual(
            _sha256(Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl")),
            ledger["orion_node_hours_jsonl_sha256"],
        )
        self.assertEqual(
            _sha256(Path("/ccs/proj/ast207/proj-shared/PIC/ledger/node_hours.jsonl")),
            ledger["project_home_node_hours_jsonl_sha256"],
        )
        self.assertEqual(len(ledger_records), ledger["orion_ledger_records"])
        self.assertEqual(
            accounting["manual_accounting_activation"]["reviewed_job_count"], 10
        )

    def test_phase0_successor_v7_binds_completed_replay_boundary(self) -> None:
        successor = _load("phase0_curated_candidate_successor_v7_2026-06-01.json")
        self.assertEqual(
            successor["status"],
            "canonical_clean_candidate_frozen_accounting_closed_"
            "registered_prerequisite_replays_authorized",
        )
        self.assertEqual(successor["remaining_phase0_actions"], [])
        self.assertEqual(
            successor["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_curated_candidate_successor_v6_2026-06-01.json",
        )
        for key in (
            "clean_candidate_freeze_receipt",
            "scheduler_accounting_successor_receipt",
            "registered_prerequisite_replay_policy_promotion_receipt",
        ):
            receipt = successor[key]
            self.assertEqual(
                receipt["sha256"],
                _sha256(REPO_ROOT / receipt["path"]),
            )

    def test_phase0_successor_v8_binds_terminal_replay_closure(self) -> None:
        successor = _load("phase0_curated_candidate_successor_v8_2026-06-01.json")
        self.assertEqual(
            successor["status"],
            "canonical_clean_candidate_frozen_registered_prerequisite_replays_pass",
        )
        self.assertEqual(
            successor["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_curated_candidate_successor_v7_2026-06-01.json",
        )
        receipt = successor["registered_prerequisite_replay_closure_receipt"]
        self.assertEqual(receipt["sha256"], _sha256(REPO_ROOT / receipt["path"]))

    def test_phase0_successor_v15_binds_d720_policy_promotion(self) -> None:
        successor = _load("phase0_curated_candidate_successor_v15_2026-06-02.json")
        self.assertEqual(
            successor["status"],
            "canonical_vl2_tsc_clean_candidate_authorized_policy_promoted_"
            "launch_prohibited_pending_registered_science_slices",
        )
        self.assertEqual(
            successor["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_curated_candidate_successor_v14_2026-06-01.json",
        )
        self.assertEqual(
            successor["predecessor_sha256"],
            _sha256(REPO_ROOT / successor["predecessor_record"]),
        )
        receipt_binding = successor[
            "clean_candidate_freeze_and_policy_promotion_receipt"
        ]
        self.assertEqual(
            receipt_binding["sha256"],
            _sha256(REPO_ROOT / receipt_binding["path"]),
        )
        receipt = json.loads(
            (REPO_ROOT / receipt_binding["path"]).read_text(encoding="utf-8")
        )
        self.assertEqual(
            receipt["predecessor_sha256"],
            _sha256(REPO_ROOT / receipt["predecessor_record"]),
        )
        promotion = receipt["active_policy_promotion"]
        self.assertEqual(promotion["registered_science_slices"], [])
        self.assertEqual(
            promotion["frontier_launch_authorization"],
            "none_no_registered_science_slices",
        )
        self.assertEqual(
            promotion["repo_policy_sha256"],
            promotion["orion_policy_sha256"],
        )
        self.assertEqual(
            promotion["orion_policy_sha256"],
            promotion["project_home_policy_sha256"],
        )
        self.assertEqual(
            promotion["orion_promotion_sha256"],
            promotion["project_home_promotion_sha256"],
        )
        attestation = receipt["pre_promotion_operator_attestation"]
        self.assertEqual(_sha256(Path(attestation["path"])), attestation["sha256"])
        ledger = receipt["mirrored_ledger_invariant"]
        self.assertEqual(
            ledger["orion_node_hours_jsonl_sha256"],
            ledger["project_home_node_hours_jsonl_sha256"],
        )
        for key in (
            "orion_node_hours_jsonl_sha256",
            "project_home_node_hours_jsonl_sha256",
            "orion_node_hours_csv_sha256",
            "orion_mirror_receipts_jsonl_sha256",
        ):
            self.assertRegex(ledger[key], r"^[0-9a-f]{64}$")
        self.assertEqual(ledger["orion_ledger_records"], 108)
        self.assertEqual(ledger["active_reservations"], 0)
        self.assertEqual(receipt["external_review"]["reviewer"], "pending external review")

    def test_phase0_paired_successor_v4_binds_c83e_policy_promotion(self) -> None:
        receipt = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v4_2026-06-02.json"
        )
        self.assertEqual(
            receipt["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v3_2026-06-01.json",
        )
        self.assertEqual(
            receipt["predecessor_sha256"],
            _sha256(REPO_ROOT / receipt["predecessor_record"]),
        )
        installed = receipt["paired_install"]
        for key in ("orion", "project_home"):
            self.assertEqual(
                installed[f"{key}_inventory_sha256"],
                _sha256(Path(installed[f"{key}_root"]) / "inventory.json"),
            )
        promotion = receipt["active_policy_promotion"]
        self.assertEqual(promotion["registered_science_slices"], [])
        self.assertEqual(
            promotion["frontier_launch_authorization"],
            "none_no_registered_science_slices",
        )
        self.assertEqual(
            promotion["orion_policy_sha256"],
            promotion["project_home_policy_sha256"],
        )
        self.assertEqual(
            promotion["orion_promotion_sha256"],
            promotion["project_home_promotion_sha256"],
        )
        attestation = receipt["pre_promotion_operator_attestation"]
        self.assertEqual(_sha256(Path(attestation["path"])), attestation["sha256"])
        ledger = receipt["terminal_mirrored_ledger"]
        self.assertEqual(
            ledger["orion_node_hours_jsonl_sha256"],
            ledger["project_home_node_hours_jsonl_sha256"],
        )
        for key in (
            "orion_node_hours_jsonl_sha256",
            "project_home_node_hours_jsonl_sha256",
            "orion_node_hours_csv_sha256",
            "orion_mirror_receipts_jsonl_sha256",
        ):
            self.assertRegex(ledger[key], r"^[0-9a-f]{64}$")
        self.assertEqual(ledger["ledger_records"], 108)
        self.assertEqual(ledger["active_reservations"], 0)

    def test_phase0_successor_v16_binds_c83e_transition(self) -> None:
        successor = _load("phase0_curated_candidate_successor_v16_2026-06-02.json")
        self.assertEqual(
            successor["status"],
            "canonical_vl2_tsc_clean_candidate_authorized_successor_"
            "control_plane_promoted_launch_prohibited_pending_registered_"
            "science_slices",
        )
        self.assertEqual(
            successor["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_curated_candidate_successor_v15_2026-06-02.json",
        )
        self.assertEqual(
            successor["predecessor_sha256"],
            _sha256(REPO_ROOT / successor["predecessor_record"]),
        )
        for key in (
            "paired_control_plane_install_and_policy_promotion_receipt",
            "clean_candidate_freeze_and_policy_promotion_receipt",
        ):
            receipt = successor[key]
            self.assertEqual(
                receipt["sha256"],
                _sha256(REPO_ROOT / receipt["path"]),
            )
        baseline = successor["operational_baseline"]
        self.assertEqual(baseline["registered_science_slices"], [])
        self.assertEqual(
            successor["frontier_launch_authorization"],
            "none_no_registered_science_slices",
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
        current_repair = _load(
            "q011_section54_eleventh_frontier_sbatch_token_transition_"
            "2026-06-03.json"
        )
        repaired_staged = current_repair["repaired_staged_control_plane"]
        if staged_version == repaired_staged["version"]:
            self.assertEqual(
                current_repair["predecessor_sha256"],
                _sha256(REPO_ROOT / current_repair["predecessor_record"]),
            )
            self.assertEqual(
                repaired_staged["state"],
                "source_local_repair_staged_pending_clean_commit_worker_validation_"
                "install_build_and_freeze",
            )
            prepared = repaired_staged["prepared_artifacts"]
            prepared_path = REPO_ROOT / prepared["inventory_path"]
            prepared_inventory = json.loads(prepared_path.read_text(encoding="utf-8"))
            self.assertEqual(prepared["inventory_sha256"], _sha256(prepared_path))
            self.assertEqual(
                prepared["paper_deck_count"], len(prepared_inventory["paper_decks"])
            )
            self.assertEqual(
                prepared["publication_analyzer_count"],
                len(prepared_inventory["analyzers"]),
            )
            self.assertNotEqual(
                storage["installed_control_plane_version"], staged_version
            )
            self.assertEqual(
                storage["installed_control_plane_version"],
                storage["staged_control_plane_candidate_version"],
            )
            return
        strict_q011 = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v6_2026-06-02.json"
        )
        if (
            staged_version == strict_q011["control_plane_version"]
            and storage["installed_control_plane_version"]
            == strict_q011["control_plane_version"]
        ):
            self.assertEqual(
                strict_q011["predecessor_sha256"],
                _sha256(REPO_ROOT / strict_q011["predecessor_record"]),
            )
            installed = strict_q011["paired_install"]
            for key in ("orion", "project_home"):
                self.assertEqual(
                    installed[f"{key}_inventory_sha256"],
                    _sha256(Path(installed[f"{key}_root"]) / "inventory.json"),
                )
            self.assertEqual(
                installed["orion_inventory_sha256"],
                installed["project_home_inventory_sha256"],
            )
            promotion = _load(
                "phase0_curated_candidate_successor_v20_2026-06-02.json"
            )["policy_promotion"]
            self.assertEqual(
                promotion["repo_policy_sha256"],
                _sha256(READINESS_DIR / "storage_policy.json"),
            )
            for key in (
                "orion_policy",
                "project_home_policy",
                "orion_promotion",
                "project_home_promotion",
            ):
                self.assertEqual(
                    promotion[f"{key}_sha256"], _sha256(Path(promotion[f"{key}_path"]))
                )
            self.assertEqual(promotion["registered_science_slice_count"], 4)
            self.assertEqual(
                strict_q011["scheduler_isolation"],
                {
                    "status": "not_claimed_empty_allowlist_transition_only",
                    "reason": (
                        "unrelated_user_scheduler_job_4754394_active_"
                        "during_policy_transition"
                    ),
                },
            )
            terminal = strict_q011["terminal_mirrored_ledger"]
            self.assertEqual(
                terminal["orion_node_hours_jsonl_sha256"],
                "dfb6e24a431c4a44682f864be2dbd5520816488d58cada7f0672aeea39c204b0",
            )
            self.assertEqual(
                terminal["project_home_node_hours_jsonl_sha256"],
                terminal["orion_node_hours_jsonl_sha256"],
            )
            self.assertEqual(
                terminal["orion_node_hours_csv_sha256"],
                "92a833de7096c81954e9fed572e51c0706b61b54adc3d1704d7130652d47f946",
            )
            self.assertEqual(
                terminal["orion_mirror_receipts_jsonl_sha256"],
                "696e93105d386ea0cd77c14398552a0f17581a65e7714625e51d0b8d84b131d5",
            )
            self.assertEqual(terminal["ledger_records"], 111)
            self.assertEqual(terminal["latest_sequence_number"], 110)
            self.assertEqual(
                terminal["latest_event_sha256"],
                "fc0082bef800733c433395d48f551559abe083367e5549cc4bffbe8a0ab48bfa",
            )
            self.assertEqual(terminal["active_reservations"], 0)
            successor = _load("phase0_curated_candidate_successor_v18_2026-06-02.json")
            self.assertEqual(
                successor["predecessor_sha256"],
                _sha256(REPO_ROOT / successor["predecessor_record"]),
            )
            receipt = successor[
                "paired_control_plane_install_and_policy_promotion_receipt"
            ]
            self.assertEqual(receipt["sha256"], _sha256(REPO_ROOT / receipt["path"]))
            self.assertEqual(
                successor["live_paired_control_plane_version"], staged_version
            )
            self.assertEqual(
                successor["frontier_launch_authorization"],
                "none_no_registered_science_slices",
            )
            self.assertEqual(
                successor["operational_baseline"]["scheduler_isolation"],
                "not_claimed_unrelated_user_scheduler_job_4754394_active_"
                "empty_allowlist_transition_only",
            )
            return
        registered_pilot = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v5_2026-06-02.json"
        )
        if staged_version == registered_pilot["control_plane_version"]:
            self.assertEqual(
                registered_pilot["predecessor_sha256"],
                _sha256(REPO_ROOT / registered_pilot["predecessor_record"]),
            )
            installed = registered_pilot["paired_install"]
            for key in ("orion", "project_home"):
                self.assertEqual(
                    installed[f"{key}_inventory_sha256"],
                    _sha256(Path(installed[f"{key}_root"]) / "inventory.json"),
                )
            self.assertEqual(
                installed["orion_inventory_sha256"],
                installed["project_home_inventory_sha256"],
            )
            return
        q011_retry = _load("phase0_curated_candidate_successor_v17_2026-06-02.json")
        if staged_version == q011_retry["successor_source_control_plane_version"]:
            self.assertEqual(
                q011_retry["predecessor_record"],
                "tst/publication/readiness/"
                "phase0_curated_candidate_successor_v16_2026-06-02.json",
            )
            self.assertEqual(
                q011_retry["predecessor_sha256"],
                _sha256(REPO_ROOT / q011_retry["predecessor_record"]),
            )
            self.assertEqual(
                q011_retry["live_paired_control_plane_version"],
                registered_pilot["control_plane_version"],
            )
            registration = q011_retry["q011_retry_registration"]
            self.assertEqual(
                registration["sha256"], _sha256(REPO_ROOT / registration["path"])
            )
            prepared = q011_retry["prepared_artifacts"]
            prepared_path = REPO_ROOT / prepared["inventory_path"]
            prepared_inventory = json.loads(prepared_path.read_text(encoding="utf-8"))
            self.assertEqual(prepared["inventory_sha256"], _sha256(prepared_path))
            self.assertEqual(prepared["paper_deck_count"], len(prepared_inventory["paper_decks"]))
            self.assertEqual(
                prepared["publication_analyzer_count"], len(prepared_inventory["analyzers"])
            )
            staged = q011_retry["staged_control_plane"]
            self.assertEqual(staged["version"], staged_version)
            self.assertEqual(staged["inventoried_file_count"], len(CONTROL_PLANE_FILES))
            self.assertEqual(q011_retry["qualification_effect"], "none")
            self.assertEqual(
                q011_retry["frontier_launch_authorization"],
                "none_live_8f0a9d7f_empty_registered_science_allowlist",
            )
            return
        replay = _load(
            "phase0_registered_prerequisite_replay_policy_promotion_2026-06-01.json"
        )
        if (
            storage["installed_control_plane_version"] == replay["control_plane_version"]
            and staged_version == replay["control_plane_version"]
        ):
            self._assert_current_registered_replay_bindings(policy, staged_version)
            return
        if storage["installed_control_plane_version"] == replay["control_plane_version"]:
            successor = _load("phase0_curated_candidate_successor_v11_2026-06-01.json")
            self.assertEqual(
                successor["predecessor_record"],
                "tst/publication/readiness/"
                "phase0_curated_candidate_successor_v10_2026-06-01.json",
            )
            self.assertEqual(
                successor["predecessor_sha256"],
                _sha256(READINESS_DIR / successor["predecessor_record"].split("/")[-1]),
            )
            self.assertEqual(
                successor["status"],
                "local_vl2_tsc_successor_staged_validation_rereview_"
                "install_build_and_freeze_pending",
            )
            self.assertEqual(successor["qualification_effect"], "none")
            self.assertEqual(
                successor["candidate_freeze_source_commit"],
                "pending_final_clean_receipt_commit",
            )
            self.assertEqual(
                successor["frontier_launch_authorization"],
                "none_pending_validation_rereview_paired_install_"
                "clean_build_and_freeze",
            )
            self.assertEqual(
                successor["live_paired_control_plane_version"],
                storage["installed_control_plane_version"],
            )
            self.assertEqual(
                successor["successor_source_control_plane_version"], staged_version
            )
            predecessor_commit = successor["curated_source_predecessor_commit"]
            self.assertRegex(predecessor_commit, r"^[0-9a-f]{40}$")
            subprocess.check_call(
                ["git", "cat-file", "-e", f"{predecessor_commit}^{{commit}}"],
                cwd=REPO_ROOT,
            )
            prepared = successor["prepared_artifacts"]
            prepared_path = REPO_ROOT / prepared["inventory_path"]
            self.assertEqual(_sha256(prepared_path), prepared["inventory_sha256"])
            prepared_inventory = json.loads(
                prepared_path.read_text(encoding="utf-8")
            )
            expected_decks = sorted(
                [
                    *(REPO_ROOT / "inputs" / "tests").glob("pic*.athinput"),
                    *(
                        REPO_ROOT / path
                        for path in PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS
                    ),
                ]
            )
            expected_analyzers = sorted(
                (REPO_ROOT / "tst" / "publication").glob("analyze_*.py")
            )
            self.assertEqual(len(expected_decks), prepared["paper_deck_count"])
            self.assertEqual(
                len(expected_analyzers), prepared["publication_analyzer_count"]
            )
            for key, expected_paths in {
                "paper_decks": expected_decks,
                "analyzers": expected_analyzers,
            }.items():
                records = prepared_inventory[key]
                self.assertEqual(
                    [record["path"] for record in records],
                    [
                        path.relative_to(REPO_ROOT).as_posix()
                        for path in expected_paths
                    ],
                )
                for record in records:
                    self.assertEqual(
                        _sha256(REPO_ROOT / record["path"]), record["sha256"]
                    )
            baseline = successor["operational_baseline"]
            self.assertEqual(
                baseline["repo_policy_sha256"],
                _sha256(READINESS_DIR / "storage_policy.json"),
            )
            for key in (
                "orion_policy",
                "project_home_policy",
                "orion_promotion",
                "project_home_promotion",
            ):
                self.assertEqual(
                    baseline[f"{key}_sha256"], _sha256(Path(baseline[f"{key}_path"]))
            )
            return
        c83e_successor = _load(
            "phase0_curated_candidate_successor_v16_2026-06-02.json"
        )
        c83e_receipt_binding = c83e_successor[
            "paired_control_plane_install_and_policy_promotion_receipt"
        ]
        c83e_receipt = json.loads(
            (REPO_ROOT / c83e_receipt_binding["path"]).read_text(encoding="utf-8")
        )
        if (
            storage["installed_control_plane_version"]
            == c83e_receipt["control_plane_version"]
            and policy["science_submission_freeze"]
            == c83e_receipt["active_policy_promotion"]["science_submission_freeze"]
        ):
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                c83e_receipt["control_plane_version"],
            )
            self.assertEqual(staged_version, c83e_receipt["control_plane_version"])
            self.assertEqual(policy["registered_science_slices"], [])
            self.assertEqual(
                c83e_successor["predecessor_sha256"],
                _sha256(
                    READINESS_DIR
                    / c83e_successor["predecessor_record"].split("/")[-1]
                ),
            )
            self.assertEqual(
                c83e_receipt_binding["sha256"],
                _sha256(REPO_ROOT / c83e_receipt_binding["path"]),
            )
            baseline = c83e_successor["operational_baseline"]
            self.assertEqual(
                baseline["repo_policy_sha256"],
                _sha256(READINESS_DIR / "storage_policy.json"),
            )
            for key in (
                "orion_policy",
                "project_home_policy",
                "orion_promotion",
                "project_home_promotion",
            ):
                self.assertEqual(
                    baseline[f"{key}_sha256"],
                    _sha256(Path(baseline[f"{key}_path"])),
                )
            return
        d720_successor = _load(
            "phase0_curated_candidate_successor_v15_2026-06-02.json"
        )
        d720_receipt_binding = d720_successor[
            "clean_candidate_freeze_and_policy_promotion_receipt"
        ]
        d720_receipt = json.loads(
            (REPO_ROOT / d720_receipt_binding["path"]).read_text(encoding="utf-8")
        )
        if (
            storage["installed_control_plane_version"]
            == d720_receipt["control_plane_version"]
            and policy["science_submission_freeze"]
            == d720_receipt["active_policy_promotion"]["science_submission_freeze"]
        ):
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                d720_receipt["control_plane_version"],
            )
            self.assertEqual(policy["registered_science_slices"], [])
            self.assertEqual(
                d720_successor["predecessor_sha256"],
                _sha256(
                    READINESS_DIR
                    / d720_successor["predecessor_record"].split("/")[-1]
                ),
            )
            self.assertEqual(
                d720_receipt_binding["sha256"],
                _sha256(REPO_ROOT / d720_receipt_binding["path"]),
            )
            baseline = d720_successor["operational_baseline"]
            for key in (
                "orion_policy",
                "project_home_policy",
                "orion_promotion",
                "project_home_promotion",
            ):
                self.assertEqual(
                    baseline[f"{key}_sha256"],
                    _sha256(Path(baseline[f"{key}_path"])),
                )
            return
        transition = _load("phase0_curated_candidate_successor_v14_2026-06-01.json")
        if (
            storage["installed_control_plane_version"]
            == transition["successor_source_control_plane_version"]
        ):
            self.assertEqual(
                storage["staged_control_plane_candidate_version"], staged_version
            )
            self.assertEqual(
                transition["successor_source_control_plane_version"], staged_version
            )
            self.assertEqual(
                transition["predecessor_sha256"],
                _sha256(READINESS_DIR / transition["predecessor_record"].split("/")[-1]),
            )
            self.assertEqual(
                policy["science_submission_freeze"],
                {"status": "pending_clean_candidate_freeze"},
            )
            self.assertEqual(policy["registered_science_slices"], [])
            clean_candidate = transition["clean_candidate_freeze_receipt"]
            self.assertEqual(
                clean_candidate["sha256"],
                _sha256(REPO_ROOT / clean_candidate["path"]),
            )
            prepared = transition["prepared_artifacts"]
            self.assertEqual(
                prepared["inventory_sha256"],
                _sha256(REPO_ROOT / prepared["inventory_path"]),
            )
            return
        phase0_successor = _load(
            "phase0_curated_candidate_successor_v6_2026-06-01.json"
        )
        self.assertEqual(
            phase0_successor["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_curated_candidate_successor_v4_2026-05-31.json",
        )
        predecessor_commit = phase0_successor["curated_source_predecessor_commit"]
        self.assertRegex(predecessor_commit, r"^[0-9a-f]{40}$")
        subprocess.check_call(
            ["git", "cat-file", "-e", f"{predecessor_commit}^{{commit}}"],
            cwd=REPO_ROOT,
        )
        self.assertEqual(
            staged_version, phase0_successor["successor_source_control_plane_version"]
        )
        prepared = phase0_successor["prepared_artifacts"]
        prepared_path = REPO_ROOT / prepared["inventory_path"]
        self.assertEqual(_sha256(prepared_path), prepared["inventory_sha256"])
        prepared_inventory = json.loads(prepared_path.read_text(encoding="utf-8"))
        expected_decks = sorted(
            [
                *(REPO_ROOT / "inputs" / "tests").glob("pic*.athinput"),
                *(
                    REPO_ROOT / path
                    for path in PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS
                ),
            ]
        )
        expected_analyzers = sorted(
            (REPO_ROOT / "tst" / "publication").glob("analyze_*.py")
        )
        self.assertEqual(len(expected_decks), prepared["paper_deck_count"])
        self.assertEqual(
            len(expected_analyzers), prepared["publication_analyzer_count"]
        )
        for key, expected_paths in {
            "paper_decks": expected_decks,
            "analyzers": expected_analyzers,
        }.items():
            records = prepared_inventory[key]
            self.assertEqual(
                [record["path"] for record in records],
                [path.relative_to(REPO_ROOT).as_posix() for path in expected_paths],
            )
            for record in records:
                self.assertEqual(
                    _sha256(REPO_ROOT / record["path"]),
                    record["sha256"],
                )
        paired_binding = phase0_successor[
            "paired_install_and_policy_promotion_receipt"
        ]
        paired_path = REPO_ROOT / paired_binding["path"]
        self.assertEqual(_sha256(paired_path), paired_binding["sha256"])
        paired = json.loads(paired_path.read_text(encoding="utf-8"))
        self.assertEqual(
            paired["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_paired_control_plane_install_and_policy_promotion_successor_2026-05-31.json",
        )
        self.assertEqual(paired["control_plane_version"], staged_version)
        installed = paired["paired_install"]
        self.assertEqual(
            installed["orion_inventory_sha256"],
            installed["project_home_inventory_sha256"],
        )
        self.assertEqual(installed["inventoried_file_count"], len(CONTROL_PLANE_FILES))
        self.assertEqual(installed["generation_directory_mode"], "0555")
        self.assertTrue(installed["byte_identical_inventory"])
        self.assertEqual(installed["installed_pair_verify"], "pass")
        promotion = paired["active_policy_promotion"]
        self.assertEqual(promotion["maximum_node_hours"], 10000)
        self.assertEqual(promotion["registered_science_slices"], [])
        self.assertEqual(
            promotion["orion_policy_sha256"],
            promotion["project_home_policy_sha256"],
        )
        self.assertEqual(
            promotion["orion_policy_sha256"],
            _sha256(READINESS_DIR / "storage_policy.json"),
        )
        self.assertEqual(
            promotion["orion_promotion_sha256"],
            promotion["project_home_promotion_sha256"],
        )
        self.assertTrue(promotion["byte_identical_policy"])
        self.assertTrue(promotion["byte_identical_promotion"])
        self.assertEqual(
            storage["installed_control_plane_version"],
            storage["staged_control_plane_candidate_version"],
        )
        self.assertEqual(
            staged_version,
            storage["installed_control_plane_version"],
        )
        if policy["science_submission_freeze"] == {
            "status": "pending_clean_candidate_freeze"
        }:
            self.assertEqual(policy["registered_science_slices"], [])
            return
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
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
            "f2-parser-multirank-runtime-metadata-v2": {
                "job_script_sha256": (
                    REPO_ROOT
                    / "tst/publication/frontier_f2_structured_multirank_runtime_metadata_job.sh"
                ),
                "input_deck_sha256": (
                    REPO_ROOT / "inputs/tests/pic_parser_contract_guards.athinput"
                ),
                "analysis_script_sha256": [
                    REPO_ROOT
                    / "tst/publication/frontier_f2_multirank_runtime_metadata_analysis.py",
                    REPO_ROOT
                    / "tst/publication/frontier_f1_structured_artifacts.py",
                ],
            },
        }
        environment_path = CONTROL_PLANE_DIR / "frontier_pic_environment.sh"
        f2_candidate = _load(
            "q027_frontier_f2_multirank_runtime_metadata_candidate_2026-05-30.json"
        )
        successor_records = [
            *successor["registered_science_slices"],
            f2_candidate["registered_science_slice"],
        ]
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
        fixture = provenance["local_review_fixture"]
        fixture_path = REPO_ROOT / fixture["path"]
        decoded = base64.b64decode(
            fixture_path.read_bytes().replace(b"\n", b""),
            validate=True,
        )
        stderr_entry = failed["stderr_inventory_entry"]
        inventory_path = REPO_ROOT / failed["artifact_inventory_fixture_path"]
        self.assertEqual(_sha256(inventory_path), failed["artifact_inventory_sha256"])
        inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
        inventory_records = {
            record["path"]: record for record in inventory["files"]
        }
        self.assertEqual(_sha256(fixture_path), fixture["encoded_file_sha256"])
        self.assertEqual(hashlib.sha256(decoded).hexdigest(), fixture["decoded_sha256"])
        self.assertEqual(len(decoded), fixture["decoded_size"])
        self.assertEqual(fixture["decoded_sha256"], stderr_entry["sha256"])
        self.assertEqual(fixture["decoded_size"], stderr_entry["size"])
        self.assertEqual(inventory_records[stderr_entry["path"]], stderr_entry)
        self.assertTrue(fixture["verified_byte_identical_to_live_immutable_stderr"])
        binding = provenance["reviewed_validator_binding"]
        commit = binding["source_commit"]
        self.assertEqual(
            subprocess.check_output(
                ["git", "cat-file", "-t", commit],
                cwd=REPO_ROOT,
                text=True,
            ).strip(),
            "commit",
        )
        for key, relative_key in {
            "support_module_sha256": "support_module",
            "gyro_analyzer_sha256": "gyro_analyzer",
            "paper_coupling_analyzer_sha256": "paper_coupling_analyzer",
        }.items():
            self.assertEqual(binding[key], _git_blob_sha256(commit, binding[relative_key]))
        self.assertEqual(
            provenance["disposition"],
            "pass_historical_transcript_bound_to_local_fixture_retry_requires_separate_v2_policy_activation",
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

    def test_rejected_pre_reservation_manifest_live_preflight_contract_is_explicit(self) -> None:
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
            [
                chronology["submission_id"],
                queue_chronology["submission_id"],
                successor["coupling_v2_submission_artifact_scrub_transition"][
                    "failed_manifest_creation_submission_id"
                ],
            ],
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

    @unittest.skipUnless(
        os.environ.get("PIC_RUN_LIVE_PREFLIGHT") == "1",
        "set PIC_RUN_LIVE_PREFLIGHT=1 for Orion live-state checks",
    )
    def test_rejected_pre_reservation_manifest_live_preflight(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        preflight = successor["required_live_preflight"]
        ledger_path = Path(
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl"
        )
        live_surfaces = (
            ledger_path,
            Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.csv"),
            Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/mirror_receipts.jsonl"),
            Path("/ccs/proj/ast207/proj-shared/PIC/ledger/node_hours.jsonl"),
        )
        line_counts = [len(path.read_text().splitlines()) for path in live_surfaces]
        self.assertEqual(line_counts[0], line_counts[2])
        self.assertEqual(line_counts[0], line_counts[3])
        self.assertEqual(line_counts[0] + 1, line_counts[1])
        self.assertFalse(
            Path(
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/pending_submission.json"
            ).exists()
        )
        for marker in incomplete_manual_accounting_marker_paths(
            ledger_path, live_surfaces[3]
        ):
            self.assertFalse(marker.exists())
        ledger_events = validate_mirrored_state(
            ledger_path,
            live_surfaces[2],
            live_surfaces[3],
        )
        active_reservations = set()
        for event in ledger_events:
            if event["event_type"] == "reservation":
                active_reservations.add(event["reservation_id"])
            elif event["event_type"] in {"reconciliation", "reservation_cancelled"}:
                active_reservations.discard(event["reservation_id"])
        self.assertEqual(active_reservations, set())
        for surface in live_surfaces:
            contents = surface.read_text()
            for submission_id in preflight[
                "historical_submission_ids_must_remain_absent_from_live_ledgers"
            ]:
                self.assertNotIn(submission_id, contents)
        queue_fixture = _load(
            "q027_frontier_f1_rejected_operator_queue_format_manifest_fixture_2026-05-30.json"
        )
        rejected_manifest = Path(queue_fixture["chronology"]["manifest_path"])
        self.assertEqual(
            sorted(path.name for path in rejected_manifest.parent.iterdir()),
            ["pre_submit_manifest.json", "snapshot"],
        )
        gyro_fixture = _load(
            "q027_frontier_f1_gyro_v2_analysis_rejection_fixture_2026-05-30.json"
        )
        self.assertIn(gyro_fixture["terminal_reconciliation_event"], ledger_events)

    def test_accepted_registered_f1_source_local_closure_is_bound(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        closure_path = successor["source_local_accepted_closure_fixture"]
        closure = json.loads((REPO_ROOT / closure_path).read_text(encoding="utf-8"))
        executions = closure["registered_executions"]
        self.assertEqual(
            set(executions),
            {
                "accepted_gyro_v3_registered_execution",
                "accepted_paper_coupling_v2_registered_execution",
            },
        )
        for key, fixture in executions.items():
            execution = successor[key]
            documents = {}
            for name, record in fixture.items():
                path = REPO_ROOT / record["path"]
                self.assertEqual(_sha256(path), record["sha256"])
                documents[name] = json.loads(path.read_text(encoding="utf-8"))

            self.assertEqual(
                fixture["pre_submit_manifest"]["sha256"],
                execution["pre_submit_manifest_sha256"],
            )
            self.assertEqual(
                fixture["artifact_inventory"]["sha256"],
                execution["artifact_inventory_sha256"],
            )
            self.assertEqual(
                fixture["analysis_result"]["sha256"],
                execution["analysis_result_sha256"],
            )
            self.assertEqual(
                fixture["offline_analysis_receipt"]["sha256"],
                execution["offline_analysis_receipt_sha256"],
            )
            self.assertEqual(
                fixture["terminal_qualification_manifest"]["sha256"],
                execution["qualification_manifest_sha256"],
            )
            pre_submit_manifest = documents["pre_submit_manifest"]
            reconciliation = documents["terminal_reconciliation_event"]
            receipt = documents["terminal_reconciliation_mirror_receipt"]
            offline_receipt = documents["offline_analysis_receipt"]
            qualification = documents["terminal_qualification_manifest"]
            result = documents["analysis_result"]
            for field, expected in {
                "submission_id": execution["submission_id"],
                "registered_science_authorization_id":
                    execution["registered_science_authorization_id"],
                "artifact_dir": execution["run_artifact_dir"],
            }.items():
                self.assertEqual(pre_submit_manifest[field], expected)
            for field, expected in {
                "submission_id": execution["submission_id"],
                "reservation_id": execution["reservation_id"],
                "job_id": execution["job_id"],
                "registered_science_authorization_id":
                    execution["registered_science_authorization_id"],
                "manifest_sha256": execution["pre_submit_manifest_sha256"],
                "artifact_dir": execution["run_artifact_dir"],
                "consumed_node_hours": execution["consumed_node_hours"],
                "cumulative_consumed_node_hours":
                    execution["cumulative_consumed_node_hours"],
            }.items():
                self.assertEqual(reconciliation[field], expected)
            self.assertEqual(
                reconciliation["event_sha256"],
                record_sha256(reconciliation, "event_sha256"),
            )
            self.assertEqual(
                receipt["mirror_ack_sha256"],
                record_sha256(receipt, "mirror_ack_sha256"),
            )
            self.assertEqual(
                receipt["mirrored_event_sha256"],
                reconciliation["event_sha256"],
            )
            self.assertEqual(
                offline_receipt["artifact_inventory"]["sha256"],
                execution["artifact_inventory_sha256"],
            )
            self.assertEqual(
                offline_receipt["analysis_result"]["sha256"],
                execution["analysis_result_sha256"],
            )
            resources = qualification["resources"]
            for field, expected in {
                "submission_id": execution["submission_id"],
                "reservation_id": execution["reservation_id"],
                "job_id": execution["job_id"],
                "registered_science_authorization_id":
                    execution["registered_science_authorization_id"],
                "pre_submit_manifest_sha256": execution["pre_submit_manifest_sha256"],
                "artifact_inventory_sha256": execution["artifact_inventory_sha256"],
                "analysis_result_sha256": execution["analysis_result_sha256"],
                "offline_analysis_receipt_sha256":
                    execution["offline_analysis_receipt_sha256"],
                "node_hours": execution["consumed_node_hours"],
            }.items():
                self.assertEqual(resources[field], expected)
            self.assertEqual(
                qualification["authorization"]["active_policy"]["sha256"],
                execution["qualification_active_policy_sha256"],
            )
            self.assertEqual(
                qualification["authorization"]["active_promotion"]["sha256"],
                execution["qualification_active_promotion_sha256"],
            )
            if key == "accepted_gyro_v3_registered_execution":
                for field in [
                    "status",
                    "cycle",
                    "particle_count",
                    "max_abs_velocity_error",
                    "velocity_tolerance",
                ]:
                    self.assertEqual(result[field], execution["analysis"][field])
                self.assertEqual(
                    fixture["initial_qualification_manifest"]["sha256"],
                    execution["initial_qualification_manifest_sha256"],
                )
            else:
                coeff0 = result["cases"]["coeff0"]
                coeff7 = result["cases"]["coeff7"]
                self.assertEqual(
                    coeff0["max_abs_momentum_conservation_error"],
                    execution["analysis"]["coeff0_max_abs_momentum_conservation_error"],
                )
                self.assertEqual(
                    coeff0["abs_energy_conservation_error"],
                    execution["analysis"]["coeff0_abs_energy_conservation_error"],
                )
                self.assertEqual(
                    coeff7["max_abs_momentum_conservation_error"],
                    execution["analysis"]["coeff7_max_abs_momentum_conservation_error"],
                )
                self.assertEqual(
                    coeff7["abs_energy_conservation_error"],
                    execution["analysis"]["coeff7_abs_energy_conservation_error"],
                )
                self.assertEqual(
                    result["coefficient_invariance"]["max_abs_particle_momentum_error"],
                    execution["analysis"][
                        "max_abs_particle_momentum_coefficient_invariance_error"
                    ],
                )

    def test_accepted_registered_f2_source_local_closure_is_bound(self) -> None:
        candidate = _load(
            "q027_frontier_f2_multirank_runtime_metadata_candidate_2026-05-30.json"
        )
        execution = candidate["accepted_v2_execution"]
        closure_path = candidate["source_local_accepted_closure_fixture"]
        closure = json.loads((REPO_ROOT / closure_path).read_text(encoding="utf-8"))
        fixture = closure["accepted_v2_registered_execution"]
        self.assertEqual(fixture["submission_id"], execution["submission_id"])
        self.assertEqual(fixture["reservation_id"], execution["reservation_id"])
        self.assertEqual(fixture["job_id"], execution["job_id"])
        self.assertEqual(
            fixture["registered_science_authorization_id"],
            execution["authorization_id"],
        )
        self.assertEqual(
            set(fixture["documents"]),
            {
                "pre_submit_manifest",
                "artifact_inventory",
                "analysis_result",
                "offline_analysis_receipt",
                "terminal_qualification_manifest",
                "terminal_reconciliation_event",
                "terminal_reconciliation_mirror_receipt",
                "pre_policy_promotion_attestation",
                "pre_manifest_attestation",
                "pre_submit_wrapper_attestation",
            },
        )
        documents = {}
        for name, record in fixture["documents"].items():
            path = REPO_ROOT / record["path"]
            self.assertEqual(_sha256(path), record["sha256"])
            documents[name] = json.loads(path.read_text(encoding="utf-8"))

        for name, digest_key in {
            "pre_submit_manifest": "pre_submit_manifest_sha256",
            "artifact_inventory": "artifact_inventory_sha256",
            "analysis_result": "analysis_result_sha256",
            "offline_analysis_receipt": "offline_analysis_receipt_sha256",
            "terminal_qualification_manifest": "qualification_manifest_sha256",
            "pre_policy_promotion_attestation":
                "pre_policy_promotion_attestation_sha256",
            "pre_manifest_attestation": "pre_manifest_attestation_sha256",
            "pre_submit_wrapper_attestation":
                "pre_submit_wrapper_attestation_sha256",
        }.items():
            self.assertEqual(
                fixture["documents"][name]["sha256"],
                execution[digest_key],
            )

        self.assertEqual(
            fixture["contained_output_path"], execution["contained_output_path"]
        )
        self.assertEqual(
            fixture["contained_output_sha256"], execution["contained_output_sha256"]
        )
        manifest = documents["pre_submit_manifest"]
        self.assertEqual(manifest["submission_id"], execution["submission_id"])
        self.assertEqual(
            manifest["registered_science_authorization_id"],
            execution["authorization_id"],
        )
        self.assertEqual(manifest["artifact_dir"], fixture["run_artifact_dir"])
        inventory = {
            record["path"]: record
            for record in documents["artifact_inventory"]["files"]
        }
        contained_output = inventory[fixture["contained_output_path"]]
        self.assertEqual(
            contained_output["sha256"], fixture["contained_output_sha256"]
        )
        analysis = documents["analysis_result"]
        self.assertEqual(analysis["status"], "pass")
        self.assertEqual(analysis["parallel_ranks"], 8)
        self.assertEqual(len(analysis["hosts"]), 1)
        self.assertEqual(len(analysis["rank_gpu_bindings"]), 8)
        self.assertEqual(
            len({
                binding["rocr_visible_device"]
                for binding in analysis["rank_gpu_bindings"]
            }),
            8,
        )
        self.assertEqual(
            analysis["runtime_artifacts"][fixture["contained_output_path"]],
            fixture["contained_output_sha256"],
        )
        receipt = documents["offline_analysis_receipt"]
        self.assertEqual(
            receipt["artifact_inventory"]["sha256"],
            execution["artifact_inventory_sha256"],
        )
        self.assertEqual(
            receipt["analysis_result"]["sha256"],
            execution["analysis_result_sha256"],
        )
        qualification = documents["terminal_qualification_manifest"]
        resources = qualification["resources"]
        for field, expected in {
            "submission_id": execution["submission_id"],
            "reservation_id": execution["reservation_id"],
            "job_id": execution["job_id"],
            "registered_science_authorization_id": execution["authorization_id"],
            "pre_submit_manifest_sha256": execution["pre_submit_manifest_sha256"],
            "artifact_inventory_sha256": execution["artifact_inventory_sha256"],
            "analysis_result_sha256": execution["analysis_result_sha256"],
            "offline_analysis_receipt_sha256":
                execution["offline_analysis_receipt_sha256"],
            "node_hours": execution["consumed_node_hours"],
        }.items():
            self.assertEqual(resources[field], expected)
        self.assertEqual(
            qualification["authorization"]["active_policy"]["sha256"],
            execution["active_policy_sha256"],
        )
        self.assertEqual(
            qualification["authorization"]["active_promotion"]["sha256"],
            execution["active_promotion_sha256"],
        )
        reconciliation = documents["terminal_reconciliation_event"]
        for field, expected in {
            "event_type": "reconciliation",
            "state": execution["state"],
            "submission_id": execution["submission_id"],
            "reservation_id": execution["reservation_id"],
            "job_id": execution["job_id"],
            "registered_science_authorization_id": execution["authorization_id"],
            "manifest_sha256": execution["pre_submit_manifest_sha256"],
            "artifact_dir": fixture["run_artifact_dir"],
            "active_policy_sha256": execution["active_policy_sha256"],
            "active_promotion_sha256": execution["active_promotion_sha256"],
            "consumed_node_hours": execution["consumed_node_hours"],
            "cumulative_consumed_node_hours":
                execution["cumulative_consumed_node_hours"],
        }.items():
            self.assertEqual(reconciliation[field], expected)
        self.assertEqual(
            reconciliation["event_sha256"],
            record_sha256(reconciliation, "event_sha256"),
        )
        self.assertEqual(
            reconciliation["event_sha256"],
            execution["terminal_reconciliation_event_sha256"],
        )
        mirror_receipt = documents["terminal_reconciliation_mirror_receipt"]
        self.assertEqual(
            mirror_receipt["mirror_ack_sha256"],
            record_sha256(mirror_receipt, "mirror_ack_sha256"),
        )
        self.assertEqual(
            mirror_receipt["mirror_ack_sha256"],
            execution["terminal_reconciliation_mirror_ack_sha256"],
        )
        self.assertEqual(
            mirror_receipt["mirrored_event_sha256"],
            reconciliation["event_sha256"],
        )

        chronology = candidate["rejected_v2_pre_reservation_chronology"]
        rejected = closure["rejected_pre_reservation_chronology"]
        rejected_manifest = rejected["queue_format_manifest"]
        self.assertEqual(rejected_manifest["submission_id"], chronology["submission_id"])
        self.assertEqual(
            rejected_manifest["pre_submit_manifest"]["sha256"],
            chronology["pre_submit_manifest_sha256"],
        )
        self.assertEqual(
            rejected_manifest["attestations"][0]["sha256"],
            chronology["pre_manifest_attestation_sha256"],
        )
        self.assertEqual(
            rejected_manifest["attestations"][1]["sha256"],
            chronology["pre_submit_wrapper_attestation_sha256"],
        )
        queue_drift_attestation = rejected["transient_queue_drift"]["attestation"]
        self.assertEqual(
            queue_drift_attestation["sha256"],
            chronology["transient_queue_drift_retry_attestation_sha256"],
        )
        for record in [
            rejected_manifest["pre_submit_manifest"],
            *rejected_manifest["attestations"],
            queue_drift_attestation,
        ]:
            self.assertEqual(_sha256(REPO_ROOT / record["path"]), record["sha256"])
        self.assertEqual(chronology["reservation_attachments"], "absent")
        self.assertEqual(chronology["ledger_intent"], "absent")
        self.assertEqual(chronology["scheduler_submission"], "absent")
        terminal = closure["active_terminal_ledger"]
        self.assertEqual(terminal["records"], execution["terminal_ledger_records"])
        self.assertEqual(
            terminal["cumulative_consumed_node_hours"],
            execution["cumulative_consumed_node_hours"],
        )
        self.assertEqual(terminal["active_reservations"], 0)
        self.assertEqual(terminal["pending_submission_marker"], "absent")

    @unittest.skipUnless(
        os.environ.get("PIC_RUN_LIVE_PREFLIGHT") == "1",
        "set PIC_RUN_LIVE_PREFLIGHT=1 for Orion live-state checks",
    )
    def test_registered_terminal_review_manifests_replay_live(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        f2_candidate = _load(
            "q027_frontier_f2_multirank_runtime_metadata_candidate_2026-05-30.json"
        )
        manual_accounting_activation = _load(
            "q027_manual_frontier_accounting_activation_2026-05-30.json"
        )
        f2 = f2_candidate["accepted_v2_execution"]
        active_policy_path = Path(
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/policy/storage_policy.json"
        )
        active_promotion_path = Path(
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/policy/active_promotion.json"
        )
        current_policy_sha256 = _sha256(active_policy_path)
        current_promotion_sha256 = _sha256(active_promotion_path)
        current_transition = manual_accounting_activation["control_plane_transition"]
        paired_transition = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v3_2026-06-01.json"
        )
        paired_promotion = paired_transition["active_policy_promotion"]
        c83e_transition = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v4_2026-06-02.json"
        )
        c83e_promotion = c83e_transition["active_policy_promotion"]
        d720_transition = _load(
            "phase0_clean_candidate_freeze_and_policy_promotion_"
            "successor_v2_2026-06-02.json"
        )
        d720_promotion = d720_transition["active_policy_promotion"]
        if current_policy_sha256 == c83e_promotion["orion_policy_sha256"]:
            self.assertEqual(
                current_promotion_sha256, c83e_promotion["orion_promotion_sha256"]
            )
            terminal = c83e_transition["terminal_mirrored_ledger"]
        elif current_policy_sha256 == d720_promotion["orion_policy_sha256"]:
            self.assertEqual(
                current_promotion_sha256, d720_promotion["orion_promotion_sha256"]
            )
            terminal = d720_transition["mirrored_ledger_invariant"]
        elif current_policy_sha256 == paired_promotion["orion_policy_sha256"]:
            self.assertEqual(
                current_promotion_sha256, paired_promotion["orion_promotion_sha256"]
            )
            terminal = paired_transition["terminal_mirrored_ledger"]
        elif current_policy_sha256 == current_transition["active_policy_sha256"]:
            self.assertEqual(
                current_promotion_sha256, current_transition["active_promotion_sha256"]
            )
            terminal = manual_accounting_activation["terminal_ledger"]
        else:
            active_policy = json.loads(active_policy_path.read_text(encoding="utf-8"))
            active_promotion = json.loads(
                active_promotion_path.read_text(encoding="utf-8")
            )
            self.assertEqual(
                active_promotion["policy_sha256"],
                current_policy_sha256,
            )
            project_home_policy_path = Path(active_promotion["project_home_policy_path"])
            self.assertEqual(_sha256(project_home_policy_path), current_policy_sha256)
            self.assertEqual(
                json.loads(project_home_policy_path.read_text(encoding="utf-8")),
                active_policy,
            )
            validate_mirrored_state(
                Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl"),
                Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/mirror_receipts.jsonl"),
                Path("/ccs/proj/ast207/proj-shared/PIC/ledger/node_hours.jsonl"),
            )
            terminal = {
                "orion_node_hours_jsonl_sha256": _sha256(
                    Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl")
                ),
                "orion_node_hours_csv_sha256": _sha256(
                    Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.csv")
                ),
                "orion_mirror_receipts_jsonl_sha256": _sha256(
                    Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/mirror_receipts.jsonl")
                ),
                "project_home_node_hours_jsonl_sha256": _sha256(
                    Path("/ccs/proj/ast207/proj-shared/PIC/ledger/node_hours.jsonl")
                ),
            }
        for digest_key, path in {
            "orion_node_hours_jsonl_sha256": Path(
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl"
            ),
            "orion_node_hours_csv_sha256": Path(
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.csv"
            ),
            "orion_mirror_receipts_jsonl_sha256": Path(
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/"
                "mirror_receipts.jsonl"
            ),
            "project_home_node_hours_jsonl_sha256": Path(
                "/ccs/proj/ast207/proj-shared/PIC/ledger/node_hours.jsonl"
            ),
        }.items():
            expected_key = digest_key
            if terminal is manual_accounting_activation["terminal_ledger"]:
                expected_key = {
                    "orion_node_hours_jsonl_sha256": "orion_jsonl_sha256",
                    "orion_node_hours_csv_sha256": "orion_csv_sha256",
                    "orion_mirror_receipts_jsonl_sha256":
                        "orion_mirror_receipts_sha256",
                    "project_home_node_hours_jsonl_sha256":
                        "project_home_jsonl_sha256",
                }[digest_key]
            self.assertEqual(_sha256(path), terminal[expected_key])
        for key in [
            "accepted_gyro_v3_registered_execution",
            "accepted_paper_coupling_v2_registered_execution",
        ]:
            with self.subTest(execution=key):
                execution = successor[key]
                manifest_path = Path(execution["qualification_manifest_path"])
                with _pinned_regular_bytes(manifest_path) as manifest_bytes:
                    self.assertEqual(
                        hashlib.sha256(manifest_bytes).hexdigest(),
                        execution["qualification_manifest_sha256"],
                    )
                    manifest = json.loads(manifest_bytes)
                for path_key, sha256_key in {
                    "pre_submit_manifest_path": "pre_submit_manifest_sha256",
                    "artifact_inventory_path": "artifact_inventory_sha256",
                    "analysis_result_path": "analysis_result_sha256",
                    "offline_analysis_receipt_path":
                        "offline_analysis_receipt_sha256",
                }.items():
                    with _pinned_regular_bytes(
                        Path(manifest["resources"][path_key])
                    ) as resource_bytes:
                        self.assertEqual(
                            hashlib.sha256(resource_bytes).hexdigest(),
                            manifest["resources"][sha256_key],
                        )
                manifest_policy_sha256 = manifest["authorization"][
                    "active_policy"
                ]["sha256"]
                manifest_promotion_sha256 = manifest["authorization"][
                    "active_promotion"
                ]["sha256"]
                self.assertEqual(
                    manifest_policy_sha256,
                    execution["qualification_active_policy_sha256"],
                )
                self.assertEqual(
                    manifest_promotion_sha256,
                    execution["qualification_active_promotion_sha256"],
                )
                validate_schema(manifest, VALIDATION_MANIFEST_SCHEMA)
        with _pinned_regular_bytes(
            Path(f2["qualification_manifest_path"])
        ) as manifest_bytes:
            self.assertEqual(
                hashlib.sha256(manifest_bytes).hexdigest(),
                f2["qualification_manifest_sha256"],
            )
            validate_schema(json.loads(manifest_bytes), VALIDATION_MANIFEST_SCHEMA)

    def test_registered_parser_policy_transition_resolves_source_commit(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        transition = successor["gyro_v3_parser_policy_transition"]
        commit = transition["source_commit"]
        self.assertEqual(
            subprocess.check_output(
                ["git", "cat-file", "-t", commit],
                cwd=REPO_ROOT,
                text=True,
            ).strip(),
            "commit",
        )
        slices = {
            record["authorization_id"]: record
            for record in transition["registered_science_slices"]
        }
        for authorization_id, analyzer in {
            "f1-clean-gyro-mpich-stderr-v3":
                "tst/publication/frontier_f1_gpu_relativistic_gyro_analysis.py",
            "f1-clean-paper-coupling-mpich-stderr-v2":
                "tst/publication/frontier_f1_gpu_paper_coupling_analysis.py",
        }.items():
            self.assertEqual(
                slices[authorization_id]["analysis_script_sha256"],
                _git_blob_sha256(commit, analyzer),
            )
            self.assertEqual(
                slices[authorization_id]["analysis_support_sha256"],
                _git_blob_sha256(
                    commit,
                    "tst/publication/frontier_f1_structured_artifacts.py",
                ),
            )

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
        binding = fixture["manifest_binding"]
        manifest_path = REPO_ROOT / binding["local_fixture_path"]
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
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
        for attestation, digest_key in zip(
            fixture["attestations"],
            (
                "pre_manifest_attestation_sha256",
                "pre_submit_wrapper_attestation_sha256",
            ),
        ):
            path = REPO_ROOT / attestation["local_fixture_path"]
            self.assertEqual(attestation["sha256"], chronology[digest_key])
            self.assertEqual(_sha256(path), chronology[digest_key])
            contents = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(contents["phase"], attestation["phase"])
            self.assertEqual(contents["pending_submission_marker"]["value"], "absent")
            self.assertEqual(
                sorted(contents["mirrored_ledger_line_counts"]["counts"].values()),
                [
                    chronology["orion_ledger_records"],
                    chronology["orion_receipt_records"],
                    chronology["project_home_ledger_records"],
                ],
            )
            self.assertEqual(
                contents["queue_snapshot"]["sha256"],
                binding["queue_snapshot_sha256"],
            )

    def test_coupling_scrub_transition_resolves_source_commit(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        transition = successor["coupling_v2_submission_artifact_scrub_transition"]
        commit = transition["source_commit"]
        self.assertEqual(
            subprocess.check_output(
                ["git", "cat-file", "-t", commit],
                cwd=REPO_ROOT,
                text=True,
            ).strip(),
            "commit",
        )
        self.assertEqual(
            transition["scrub_safe_analysis_script_sha256"],
            _git_blob_sha256(
                commit,
                "tst/publication/frontier_f1_gpu_paper_coupling_analysis.py",
            ),
        )
        self.assertEqual(
            transition["active_policy_sha256"],
            _git_blob_sha256(
                commit,
                "tst/publication/readiness/storage_policy.json",
            ),
        )
        gyro = successor["accepted_gyro_v3_registered_execution"]
        coupling = successor["accepted_paper_coupling_v2_registered_execution"]
        f2 = _load(
            "q027_frontier_f2_multirank_runtime_metadata_candidate_2026-05-30.json"
        )["accepted_v2_execution"]
        self.assertNotEqual(
            gyro["initial_qualification_manifest_sha256"],
            gyro["qualification_manifest_sha256"],
        )
        self.assertEqual(
            coupling["active_policy_sha256"],
            transition["active_policy_sha256"],
        )
        self.assertEqual(
            coupling["active_promotion_sha256"],
            transition["active_promotion_sha256"],
        )
        for execution in (gyro, coupling):
            self.assertEqual(
                execution["qualification_active_policy_sha256"],
                f2["active_policy_sha256"],
            )
            self.assertEqual(
                execution["qualification_active_promotion_sha256"],
                f2["active_promotion_sha256"],
            )

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
        inventory_path = REPO_ROOT / fixture["artifact_inventory_fixture_path"]
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
        event = fixture["terminal_reconciliation_event"]
        self.assertEqual(event["submission_id"], chronology["submission_id"])
        self.assertEqual(event["reservation_id"], chronology["reservation_id"])
        self.assertEqual(event["job_id"], chronology["job_id"])
        self.assertEqual(event["manifest_sha256"], chronology["pre_submit_manifest_sha256"])
        self.assertEqual(event["consumed_node_hours"], chronology["consumed_node_hours"])
        self.assertEqual(
            event["event_sha256"],
            record_sha256(event, "event_sha256"),
        )
        receipt_fixture = fixture["terminal_reconciliation_mirror_receipt_fixture"]
        receipt_path = REPO_ROOT / receipt_fixture["path"]
        self.assertEqual(_sha256(receipt_path), receipt_fixture["sha256"])
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        self.assertEqual(
            receipt["mirror_ack_sha256"],
            record_sha256(receipt, "mirror_ack_sha256"),
        )
        self.assertEqual(receipt["mirrored_event_sha256"], event["event_sha256"])
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

    def test_pic_vl2_tsc_phase1_stage_ordering_review_receipt_is_bounded(self) -> None:
        receipt = _load(
            "pic_vl2_tsc_phase1_stage_ordering_review_receipt_successor_"
            "2026-06-02.json"
        )
        predecessor = receipt["predecessor"]
        predecessor_path = REPO_ROOT / predecessor["path"]
        self.assertEqual(predecessor["sha256"], _sha256(predecessor_path))
        registration = json.loads(predecessor_path.read_text(encoding="utf-8"))
        self.assertEqual(
            registration["identity_split"]["successor_mode"],
            receipt["physical_mode"],
        )
        self.assertEqual(receipt["phase"], "phase1")
        self.assertEqual(receipt["finding_id"], "PIC-P1-002")
        self.assertEqual(receipt["verification_gate"], "Q-004")
        self.assertIn("Bounded source-local", receipt["scope"])

        trace = receipt["reused_stage_trace"]
        self.assertEqual(trace["sha256"], _sha256(REPO_ROOT / trace["path"]))
        self.assertEqual(
            trace["sha256"],
            registration["active_successor_files"][trace["path"]],
        )
        self.assertIn(
            trace["listed_pass_result"],
            registration["local_validation"]["serial_oracles_passed"],
        )
        self.assertTrue(trace["reuse_only_no_new_dynamic_execution_or_artifact"])

        for binding_group in ("source_bindings", "script_bindings", "fixture_bindings"):
            for relative, expected_sha256 in receipt[binding_group].items():
                with self.subTest(binding_group=binding_group, relative=relative):
                    self.assertEqual(expected_sha256, _sha256(REPO_ROOT / relative))
                    predecessor_sha256 = registration["active_successor_files"].get(relative)
                    if predecessor_sha256 is not None:
                        self.assertEqual(expected_sha256, predecessor_sha256)

        review = receipt["stage_ordering_review"]
        expected_sequences = {
            "push": [
                "stage_1_inserted_push_is_predictor_no_op_at_x_ini",
                "stage_1_post_deposit_half_step_drift_reaches_x_mid",
                "stage_2_inserted_push_applies_midpoint_boris_kick_at_x_mid",
                "stage_2_post_deposit_half_step_drift_reaches_x_end",
            ],
            "deposition": [
                "Particles::Push",
                "Particles::SaveOldPositions",
                "Particles::ZeroMoments",
                "Particles::InitRecvMoments",
                "Particles::DepositMoments",
            ],
            "feedback_placement": [
                "MHD::RKUpdate",
                "paper_vl2_staged_particle_wrapper_chain",
                "MHD::MHDSrcTerms",
            ],
            "boundary_synchronization": [
                "Particles::DepositMoments",
                "Particles::RestrictMoments",
                "Particles::SendMoments",
                "Particles::RecvMoments",
                "Particles::ClearRecvMoments",
                "Particles::ClearSendMoments",
                "Particles::ApplyMomentPhysicalBCs",
                "Particles::ProlongateMoments",
            ],
            "migration_communication_ordering": [
                "Particles::DriftPaperCosmicRaysHalfStep",
                "after_stagen",
                "Particles::NewGID",
                "Particles::SendCnt",
                "Particles::InitRecv",
                "Particles::SendP",
                "Particles::RecvP",
                "Particles::ClearRecv",
                "Particles::ClearSend",
            ],
            "ct_ordering": [
                "MHD::CornerE",
                "MHD::EFieldSrc",
                "MHD::SendE",
                "MHD::RecvE",
                "MHD::CT",
            ],
        }
        self.assertEqual(
            set(review),
            set(expected_sequences),
        )
        for topic_name, topic in review.items():
            self.assertEqual(topic["status"], "covered_bounded_source_local_review")
            self.assertEqual(topic["sequence"], expected_sequences[topic_name])
        self.assertIn(
            "excluded from AddsCRCurrentToCT",
            review["ct_ordering"]["interpretation"],
        )

        boundary = receipt["qualification_boundary"]
        self.assertTrue(boundary["bounded_source_local_scope"])
        for key in (
            "adds_dynamic_evidence",
            "frontier_execution_authorized",
            "frontier_qualified",
            "mpi_qualified",
            "hip_qualified",
            "claim_closure",
        ):
            self.assertFalse(boundary[key])

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

    def test_q011_repaired_clean_candidate_freeze_receipt_is_frozen(self) -> None:
        receipt = _load("phase0_clean_candidate_freeze_successor_v3_2026-06-02.json")
        self.assertEqual(
            receipt["predecessor_sha256"],
            _sha256(REPO_ROOT / receipt["predecessor_record"]),
        )
        candidate = receipt["clean_candidate"]
        manifest_path = Path(candidate["manifest_path"])
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(candidate["manifest_sha256"], _sha256(manifest_path))
        self.assertEqual(candidate["manifest_sha256"], "ef527ed467995bd60fda07b5a3b09b56ea871595ace12fd64a948e246720dbe3")
        self.assertEqual(stat.S_IMODE(manifest_path.parent.stat().st_mode) & 0o222, 0)
        descendants = list(manifest_path.parent.rglob("*"))
        self.assertEqual(
            receipt["local_validation"]["expected_tree_entries"],
            len(descendants),
        )
        for descendant in descendants:
            self.assertFalse(descendant.is_symlink())
            self.assertEqual(stat.S_IMODE(descendant.stat().st_mode) & 0o222, 0)
        self.assertEqual(candidate["source_git_commit"], manifest["source"]["git_commit"])
        self.assertEqual(candidate["source_git_tree"], manifest["source"]["git_tree"])
        self.assertEqual(candidate["source_archive_sha256"], manifest["source"]["archive_sha256"])
        self.assertEqual(candidate["source_commit_sha256"], manifest["source"]["commit_sha256"])
        self.assertEqual(candidate["source_bundle_sha256"], manifest["source"]["source_bundle_sha256"])
        self.assertEqual(candidate["prepared_artifact_inventory_sha256"], manifest["prepared_artifacts"]["inventory_sha256"])
        self.assertEqual(candidate["paper_deck_count"], len(manifest["prepared_artifacts"]["paper_decks"]))
        self.assertEqual(candidate["publication_analyzer_count"], len(manifest["prepared_artifacts"]["analyzers"]))
        self.assertEqual(candidate["build_profile_sha256"], _sha256(Path(manifest["build"]["profile_path"])))
        self.assertEqual(candidate["profile_receipt_sha256"], _sha256(Path(manifest["build"]["profile_receipt_path"])))
        self.assertEqual(candidate["executable_sha256"], _sha256(Path(manifest["build"]["executable_path"])))
        fragment = receipt["post_freeze_policy_fragment"]
        self.assertEqual(fragment["sha256"], _sha256(Path(fragment["path"])))
        self.assertEqual(fragment["status"], "non_authorizing_review_fragment_only")
        fragment_payload = json.loads(Path(fragment["path"]).read_text(encoding="utf-8"))
        self.assertEqual(fragment["registered_science_slice_count"], len(fragment_payload["registered_science_slices"]))
        self.assertEqual(fragment["registered_science_slice_count"], 4)
        policy = receipt["policy_promotion"]
        self.assertEqual(policy["registered_science_slices"], [])
        self.assertEqual(
            policy["live_orion_policy_sha256"],
            "1a331137c7d83717a890fa64046890074101f0fec3de7b8b22ca41b8bb644d28",
        )
        self.assertEqual(
            policy["live_project_home_policy_sha256"],
            policy["live_orion_policy_sha256"],
        )
        self.assertEqual(
            policy["live_orion_promotion_sha256"],
            "e64db1ef1ca755e1fcf4ff86e078856e381ff31f6aab926f01658b7aa69cd73e",
        )
        self.assertEqual(
            policy["live_project_home_promotion_sha256"],
            policy["live_orion_promotion_sha256"],
        )
        successor = _load("phase0_curated_candidate_successor_v19_2026-06-02.json")
        self.assertEqual(successor["predecessor_sha256"], _sha256(REPO_ROOT / successor["predecessor_record"]))
        successor_receipt = successor["clean_candidate_freeze_receipt"]
        self.assertEqual(
            successor_receipt["sha256"],
            _sha256(REPO_ROOT / successor_receipt["path"]),
        )
        self.assertEqual(successor["frontier_launch_authorization"], "none_no_registered_science_slices")
        self.assertEqual(successor["operational_baseline"]["registered_science_slices"], [])

    def test_q011_pressure_pilot_four_slice_policy_promotion_is_current(self) -> None:
        successor = _load("phase0_curated_candidate_successor_v20_2026-06-02.json")
        self.assertEqual(
            successor["predecessor_sha256"],
            _sha256(REPO_ROOT / successor["predecessor_record"]),
        )
        promotion = successor["policy_promotion"]
        self.assertEqual(
            promotion["repo_policy_sha256"],
            _sha256(REPO_ROOT / promotion["repo_policy_path"]),
        )
        for key in (
            "orion_policy",
            "project_home_policy",
            "orion_promotion",
            "project_home_promotion",
        ):
            self.assertEqual(
                promotion[f"{key}_sha256"],
                _sha256(Path(promotion[f"{key}_path"])),
            )
        self.assertEqual(
            promotion["orion_policy_sha256"],
            promotion["project_home_policy_sha256"],
        )
        self.assertEqual(
            promotion["orion_promotion_sha256"],
            promotion["project_home_promotion_sha256"],
        )
        policy = _load("storage_policy.json")
        self.assertEqual(
            [
                record["authorization_id"]
                for record in policy["registered_science_slices"]
            ],
            promotion["registered_science_authorization_ids"],
        )
        self.assertEqual(promotion["registered_science_slice_count"], 4)
        attestation = successor["pre_policy_promotion_operator_attestation"]
        self.assertEqual(attestation["sha256"], _sha256(Path(attestation["path"])))
        self.assertEqual(
            attestation["queue_snapshot_sha256"],
            hashlib.sha256(b"").hexdigest(),
        )
        self.assertEqual(attestation["active_reservation_count"], 0)
        self.assertEqual(attestation["pending_submission_marker"], "absent")
        self.assertEqual(attestation["pending_manual_accounting_marker"], "absent")

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
