#!/usr/bin/env python3
"""Regression tests for source-controlled PIC readiness registries."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.publication.pic_qualification_manifest import validate_schema


READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
SCHEMA_DIR = READINESS_DIR / "schemas"
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
        active = _load("q027_active_control_plane_generation_2026-05-30.json")
        candidate = _load(
            "q027_control_plane_activation_lock_hardening_candidate_2026-05-30.json"
        )
        self.assertEqual(
            storage["staged_control_plane_candidate_version"],
            candidate["staged_successor"]["control_plane_version"],
        )
        lifecycle = storage["installed_control_plane_lifecycle"]
        if lifecycle == "live_active_generation_successor_staged_not_installed":
            self.assertEqual(
                storage["installed_control_plane_version"],
                active["active_generation"]["control_plane_version"],
            )
            self.assertNotEqual(
                storage["installed_control_plane_version"],
                storage["staged_control_plane_candidate_version"],
            )
        elif lifecycle == "paired_installed_reviewed_generation":
            self.assertEqual(
                storage["installed_control_plane_version"],
                storage["staged_control_plane_candidate_version"],
            )
            paired = _load("q027_outer_anchor_paired_activation_2026-05-30.json")
            self.assertEqual(
                storage["installed_control_plane_version"],
                paired["active_generation"]["control_plane_version"],
            )
            self.assertEqual(
                storage["ledger_genesis"]["event_sha256"],
                paired["anchor_migration"]["event_sha256"],
            )
            self.assertEqual(
                storage["ledger_genesis"]["mirror_ack_sha256"],
                paired["anchor_migration"]["mirror_ack_sha256"],
            )
            self.assertEqual(
                paired["mirrored_ledger_invariant"]["active_reservations"],
                0,
            )
        else:
            self.fail(f"Unknown installed-control-plane lifecycle: {lifecycle}")
        self.assertEqual(
            long_term["status"],
            "user_selected_orion_only_with_documented_durability_risk",
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
