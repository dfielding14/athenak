#!/usr/bin/env python3
"""Regression tests for the fail-closed PIC qualification-manifest gate."""

from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.publication import pic_qualification_manifest as qualification_manifest
from tst.publication.pic_qualification_manifest import (
    RUNTIME_ENVIRONMENT_ALLOWLIST_KEYS,
    _load_object_bytes,
    _require_frontier_ledger_binding,
    _source_bundle_sha256,
    _validate_environment_allowlist,
    _verify_frontier_offline_analysis,
    freeze_qualification_manifest,
    validate_qualification_manifest,
)
from control_plane_common import PRODUCTION_RUNTIME_LOADED_MODULES
from control_plane_common import PRODUCTION_RUNTIME_MODULEFILES
from control_plane_common import PRODUCTION_RUNTIME_MODULEPATH
from tst.publication.frontier_control_plane.control_plane_common import (
    BUILD_PROVENANCE_FILENAMES,
    git_commit_tree_from_bytes,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _runtime_allowlist() -> str:
    values = {key: "<unset>" for key in RUNTIME_ENVIRONMENT_ALLOWLIST_KEYS}
    values.update(
        PIC_FRONTIER_PROFILE="frontier_minimum_supported",
        HSA_XNACK="0",
        MPICH_ENV_DISPLAY="1",
        MPICH_VERSION_DISPLAY="1",
        MPICH_GPU_SUPPORT_ENABLED="1",
        SLURM_EXPORT_ENV="ALL",
        ROCM_PATH="/opt/rocm",
        LOADEDMODULES=":".join(PRODUCTION_RUNTIME_LOADED_MODULES),
        _LMFILES_=":".join(PRODUCTION_RUNTIME_MODULEFILES),
        MODULEPATH=PRODUCTION_RUNTIME_MODULEPATH,
    )
    return "".join(f"{key}={values[key]}\n" for key in RUNTIME_ENVIRONMENT_ALLOWLIST_KEYS)


class PicQualificationManifestTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        production_semantics = patch(
            "control_plane_common.require_production_build_provenance"
        )
        production_semantics.start()
        self.addCleanup(production_semantics.stop)
        environment_allowlist = _runtime_allowlist()
        for name, contents in (
            ("athena", "executable"),
            ("CMakeCache.txt", "cache"),
            ("modules.txt", "modules"),
            ("environment.txt", environment_allowlist),
            ("metrics.json", '{"relative_error": 0.0}\n'),
            ("active_policy.json", '{"fixture": "active-policy"}\n'),
            ("active_promotion.json", '{"fixture": "active-promotion"}\n'),
        ):
            (self.root / name).write_text(contents, encoding="utf-8")
        nested = self.root / "nested-repo"
        source = self.root / "source-repo"
        for repository in (nested, source):
            subprocess.run(["git", "init", str(repository)], check=True, capture_output=True)
            (repository / "tracked.txt").write_text(f"{repository.name}\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(repository), "add", "tracked.txt"], check=True)
            subprocess.run(
                [
                    "git", "-C", str(repository), "-c", "user.name=PIC Test",
                    "-c", "user.email=pic-test@example.invalid", "commit", "-m", "fixture",
                ],
                check=True, capture_output=True,
            )
        subprocess.run(
            [
                "git", "-c", "protocol.file.allow=always", "-C", str(source),
                "submodule", "add", str(nested), "kokkos",
            ],
            check=True, capture_output=True,
        )
        subprocess.run(["git", "-C", str(source), "add", "."], check=True)
        subprocess.run(
            [
                "git", "-C", str(source), "-c", "user.name=PIC Test",
                "-c", "user.email=pic-test@example.invalid", "commit", "-m", "add submodule",
            ],
            check=True, capture_output=True,
        )
        subprocess.run(
            ["git", "-C", str(source), "archive", "--format=tar",
             f"--output={self.root / 'source.tar'}", "HEAD"],
            check=True,
        )
        subprocess.run(
            ["git", "-C", str(source / "kokkos"), "archive", "--format=tar",
             f"--output={self.root / 'kokkos.tar'}", "HEAD"],
            check=True,
        )
        commit = subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
        ).strip()
        tree = subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "HEAD^{tree}"], text=True
        ).strip()
        kokkos_commit = subprocess.check_output(
            ["git", "-C", str(source / "kokkos"), "rev-parse", "HEAD"], text=True
        ).strip()
        kokkos_tree = subprocess.check_output(
            ["git", "-C", str(source / "kokkos"), "rev-parse", "HEAD^{tree}"], text=True
        ).strip()
        (self.root / "source.commit").write_bytes(
            subprocess.check_output(["git", "-C", str(source), "cat-file", "commit", "HEAD"])
        )
        (self.root / "kokkos.commit").write_bytes(
            subprocess.check_output(
                ["git", "-C", str(source / "kokkos"), "cat-file", "commit", "HEAD"]
            )
        )
        profile_submodules = [
            {
                "path": "kokkos",
                "archive_sha256": _sha256(self.root / "kokkos.tar"),
                "commit_sha256": _sha256(self.root / "kokkos.commit"),
                "git_commit": kokkos_commit,
                "git_tree": kokkos_tree,
            }
        ]
        source_bundle = _source_bundle_sha256(
            _sha256(self.root / "source.tar"),
            _sha256(self.root / "source.commit"),
            profile_submodules,
        )
        build_provenance_dir = self.root / "build_provenance"
        build_provenance_dir.mkdir()
        build_provenance_contents = {
            "configure_log": "configured\n",
            "build_log": "built\n",
            "cmake_cache": "cache\n",
            "module_list": "modules\n",
            "toolchain": "fixture\n",
            "build_invocations": (
                '{"configure": ["cmake", "-S", "."], '
                '"build": ["cmake", "--build", "."]}\n'
            ),
            "git_status_preconfigure": "",
            "git_status": "",
            "submodule_status": f" {kokkos_commit} kokkos\n",
            "environment_allowlist": environment_allowlist,
            "build_environment": "{}\n",
        }
        for label, filename in BUILD_PROVENANCE_FILENAMES.items():
            (build_provenance_dir / filename).write_text(
                build_provenance_contents[label], encoding="utf-8"
            )
        build_stem = f"{commit[:12]}.fixture"
        build_bin = (
            f"/lustre/orion/ast207/proj-shared/dfielding/PIC/bin/{commit[:12]}/fixture"
        )
        profile_provenance = {}
        for label, filename in BUILD_PROVENANCE_FILENAMES.items():
            if label in {"configure_log", "build_log"}:
                suffix = "configure.log" if label == "configure_log" else "build.log"
                path = (
                    "/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/build/"
                    f"{build_stem}.{suffix}"
                )
            else:
                path = f"{build_bin}/{filename}"
            profile_provenance[label] = {
                "path": path,
                "sha256": _sha256(build_provenance_dir / filename),
            }
        profile = {
            "schema_version": 3,
            "profile_id": "fixture",
            "authorized_source_root": "/ccs/home/dfielding/athenak-pic",
            "fresh_source_root": (
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/"
                f"build/{commit[:12]}/fixture/source"
            ),
            "git_commit": commit,
            "git_tree": tree,
            "source_archive_sha256": _sha256(self.root / "source.tar"),
            "source_commit_sha256": _sha256(self.root / "source.commit"),
            "source_bundle_sha256": source_bundle,
            "toolchain": "fixture",
            "build_invocations_sha256": profile_provenance["build_invocations"]["sha256"],
            "executable_sha256": _sha256(self.root / "athena"),
            "provenance_inputs": profile_provenance,
            "submodules": profile_submodules,
        }
        (self.root / "build_profile.json").write_text(
            json.dumps(profile), encoding="utf-8"
        )
        receipt = {
            "schema_version": 1,
            "control_plane_version": "a" * 64,
            "profile_path": f"{build_bin}/build_profile.json",
            "profile_sha256": _sha256(self.root / "build_profile.json"),
            "source_bundle_sha256": source_bundle,
            "fresh_source_root": profile["fresh_source_root"],
            "build_invocations_sha256": profile_provenance["build_invocations"]["sha256"],
            "git_status_preconfigure_sha256": profile_provenance[
                "git_status_preconfigure"
            ]["sha256"],
            "git_status_sha256": profile_provenance["git_status"]["sha256"],
            "configure_log_sha256": profile_provenance["configure_log"]["sha256"],
            "build_log_sha256": profile_provenance["build_log"]["sha256"],
            "executable_path": f"{build_bin}/athena",
            "executable_sha256": _sha256(self.root / "athena"),
        }
        (self.root / "profile_receipt.json").write_text(
            json.dumps(receipt), encoding="utf-8"
        )
        candidate = {
            "schema_version": 3,
            "freeze_id": "03a7bd9a-7d4c-4e37-a12b-46de3817eff2",
            "created_utc": "2026-05-30T12:00:00Z",
            "source": {
                "archive_path": "/original/03a7bd9a-7d4c-4e37-a12b-46de3817eff2/source.tar",
                "archive_sha256": _sha256(self.root / "source.tar"),
                "commit_path": "/original/03a7bd9a-7d4c-4e37-a12b-46de3817eff2/source.commit",
                "commit_sha256": _sha256(self.root / "source.commit"),
                "source_bundle_sha256": source_bundle,
                "git_commit": commit,
                "git_tree": tree,
                "worktree_status": "clean",
                "submodule_status": "clean_pinned_archived",
                "submodules": [
                    {
                        **profile_submodules[0],
                        "archive_path": "/original/03a7bd9a-7d4c-4e37-a12b-46de3817eff2/submodules/0000.tar",
                        "commit_path": "/original/03a7bd9a-7d4c-4e37-a12b-46de3817eff2/submodules/0000.commit",
                        "worktree_status": "clean",
                    }
                ],
            },
            "build": {
                "profile_id": "fixture",
                "profile_path": "/original/03a7bd9a-7d4c-4e37-a12b-46de3817eff2/build_profile.json",
                "profile_sha256": _sha256(self.root / "build_profile.json"),
                "profile_receipt_path": "/original/03a7bd9a-7d4c-4e37-a12b-46de3817eff2/profile_receipt.json",
                "profile_receipt_sha256": _sha256(self.root / "profile_receipt.json"),
                "source_archive_sha256": _sha256(self.root / "source.tar"),
                "source_commit_sha256": _sha256(self.root / "source.commit"),
                "source_bundle_sha256": source_bundle,
                "toolchain": "fixture",
                "build_invocations_sha256": profile_provenance["build_invocations"][
                    "sha256"
                ],
                "executable_path": "/original/03a7bd9a-7d4c-4e37-a12b-46de3817eff2/athena",
                "executable_sha256": _sha256(self.root / "athena"),
            },
        }
        (self.root / "clean_candidate_manifest.json").write_text(
            json.dumps(candidate), encoding="utf-8"
        )
        self.manifest = {
            "schema_version": 1,
            "manifest_id": "qualification-fixture-001",
            "created_utc": "2026-05-30T12:00:00Z",
            "claim_ids": ["CLAIM-PAPER-GYRO-001"],
            "test_id": "pic_relativistic_gyro_paper",
            "evidence_class": "sun_bai_2023_reproduction",
            "physical_mode": "paper_test_particle",
            "git": {
                "commit": commit,
                "tree": tree,
                "status": [],
                "source_archive": {
                    "path": "source.tar",
                    "sha256": _sha256(self.root / "source.tar"),
                },
                "source_commit": {
                    "path": "source.commit",
                    "sha256": _sha256(self.root / "source.commit"),
                },
                "source_bundle_sha256": "",
                "submodule_status": "clean_pinned_archived",
                "submodules": [
                    {
                        "path": "kokkos",
                        "archive_path": "kokkos.tar",
                        "archive_sha256": _sha256(self.root / "kokkos.tar"),
                        "commit_path": "kokkos.commit",
                        "commit_sha256": _sha256(self.root / "kokkos.commit"),
                        "git_commit": kokkos_commit,
                        "git_tree": kokkos_tree,
                        "worktree_status": "clean",
                    }
                ],
                "clean_candidate_manifest": {
                    "path": "clean_candidate_manifest.json",
                    "sha256": _sha256(self.root / "clean_candidate_manifest.json"),
                },
                "clean_candidate_build_profile": {
                    "path": "build_profile.json",
                    "sha256": _sha256(self.root / "build_profile.json"),
                },
                "clean_candidate_build_profile_receipt": {
                    "path": "profile_receipt.json",
                    "sha256": _sha256(self.root / "profile_receipt.json"),
                },
                "clean_candidate_build_provenance": {
                    label: {
                        "path": f"build_provenance/{filename}",
                        "sha256": _sha256(build_provenance_dir / filename),
                    }
                    for label, filename in BUILD_PROVENANCE_FILENAMES.items()
                },
            },
            "authorization": {
                "control_plane_version": "a" * 64,
                "build_profile_control_plane_version": "a" * 64,
                "clean_candidate_manifest_path": (
                    "/lustre/orion/ast207/proj-shared/dfielding/PIC/"
                    "clean_candidates/03a7bd9a-7d4c-4e37-a12b-46de3817eff2/"
                    "clean_candidate_manifest.json"
                ),
                "clean_candidate_manifest_sha256": _sha256(
                    self.root / "clean_candidate_manifest.json"
                ),
                "active_policy": {
                    "path": "active_policy.json",
                    "sha256": _sha256(self.root / "active_policy.json"),
                },
                "active_promotion": {
                    "path": "active_promotion.json",
                    "sha256": _sha256(self.root / "active_promotion.json"),
                },
            },
            "executable": {
                "path": "athena",
                "sha256": _sha256(self.root / "athena"),
                "cmake_cache": {
                    "path": "CMakeCache.txt",
                    "sha256": _sha256(self.root / "CMakeCache.txt"),
                },
                "modules": {
                    "path": "modules.txt",
                    "sha256": _sha256(self.root / "modules.txt"),
                },
                "environment_allowlist": {
                    "path": "environment.txt",
                    "sha256": _sha256(self.root / "environment.txt"),
                },
            },
            "parameters": {},
            "oracle": {
                "kind": "analytic",
                "reference": "bounded fixture",
                "tolerances": {"relative_error": 1.0e-6},
            },
            "metrics": [{"name": "relative_error", "value": 0.0}],
            "resources": {
                "platform": "host",
                "artifact_root": str(self.root),
            },
            "artifacts": [
                {
                    "path": "metrics.json",
                    "sha256": _sha256(self.root / "metrics.json"),
                }
            ],
            "review": {
                "reviewer": "pending external review",
                "disposition": "pending external review",
            },
        }
        self.manifest["git"]["source_bundle_sha256"] = _source_bundle_sha256(
            self.manifest["git"]["source_archive"]["sha256"],
            self.manifest["git"]["source_commit"]["sha256"],
            self.manifest["git"]["submodules"],
        )
        self.live_authorization = {
            "control_plane_version": self.manifest["authorization"]["control_plane_version"],
            "build_profile_control_plane_version": self.manifest["authorization"][
                "build_profile_control_plane_version"
            ],
            "clean_candidate_manifest_path": self.manifest["authorization"][
                "clean_candidate_manifest_path"
            ],
            "clean_candidate_manifest_sha256": self.manifest["authorization"][
                "clean_candidate_manifest_sha256"
            ],
            "active_policy_sha256": self.manifest["authorization"]["active_policy"][
                "sha256"
            ],
            "active_promotion_sha256": self.manifest["authorization"][
                "active_promotion"
            ]["sha256"],
        }
        self.authorization_patcher = patch(
            "tst.publication.pic_qualification_manifest._live_candidate_authorization",
            side_effect=self._live_authorization,
        )
        self.authorization_patcher.start()

    def tearDown(self) -> None:
        self.authorization_patcher.stop()
        self.temporary_directory.cleanup()

    def _live_authorization(self, candidate_sha256: str) -> dict[str, str]:
        if candidate_sha256 != self.live_authorization[
            "clean_candidate_manifest_sha256"
        ]:
            raise ValueError("fixture candidate is not policy-authorized")
        return self.live_authorization

    def _assert_rejected(self, manifest: dict[str, object]) -> None:
        with self.assertRaises(ValueError):
            validate_qualification_manifest(manifest)

    def test_reviewable_qualification_manifest_is_accepted(self) -> None:
        validate_qualification_manifest(self.manifest)

    def test_environment_allowlist_rejects_extra_duplicate_and_missing_keys(self) -> None:
        records = _runtime_allowlist().splitlines(keepends=True)
        modulepath = f"MODULEPATH={PRODUCTION_RUNTIME_MODULEPATH}\n"
        _validate_environment_allowlist(
            "".join(records).encode("utf-8"),
            require_frontier_values=True,
        )
        for data in [
            "".join(records + ["SECRET_TOKEN=redacted\n"]),
            "".join(records + [records[0]]),
            "".join(records[:-1]),
            "".join(
                "SLURM_EXPORT_ENV=NOT_ALL\n"
                if record.startswith("SLURM_EXPORT_ENV=")
                else record
                for record in records
            ),
            "".join(records).replace("\n", "\r\n"),
            "".join(records).replace(modulepath, f"MODULEPATH=/tmp/forged:{PRODUCTION_RUNTIME_MODULEPATH}\n"),
            "".join(records).replace(modulepath, f"MODULEPATH={PRODUCTION_RUNTIME_MODULEPATH}:/tmp/forged\n"),
            "".join(records).replace(
                modulepath,
                f"MODULEPATH={PRODUCTION_RUNTIME_MODULEPATH.split(':', 1)[1]}\n",
            ),
        ]:
            with self.subTest(data=data):
                with self.assertRaises(ValueError):
                    _validate_environment_allowlist(
                        data.encode("utf-8"),
                        require_frontier_values=True,
                    )

    def test_proxy_and_unit_evidence_are_rejected(self) -> None:
        for evidence_class in ("engineering_proxy", "unit/regression"):
            with self.subTest(evidence_class=evidence_class):
                manifest = copy.deepcopy(self.manifest)
                manifest["evidence_class"] = evidence_class
                self._assert_rejected(manifest)

    def test_dirty_source_and_unknown_claim_are_rejected(self) -> None:
        dirty = copy.deepcopy(self.manifest)
        dirty["git"]["status"] = [" M src/particles/particles.cpp"]
        self._assert_rejected(dirty)

        unknown = copy.deepcopy(self.manifest)
        unknown["claim_ids"] = ["CLAIM-DOES-NOT-EXIST"]
        self._assert_rejected(unknown)

    def test_submodule_bundle_and_clean_candidate_drift_are_rejected(self) -> None:
        dirty = copy.deepcopy(self.manifest)
        dirty["git"]["submodules"][0]["worktree_status"] = "dirty"
        self._assert_rejected(dirty)

        mismatched = copy.deepcopy(self.manifest)
        mismatched["git"]["source_bundle_sha256"] = "0" * 64
        self._assert_rejected(mismatched)

        candidate = copy.deepcopy(self.manifest)
        candidate["git"]["clean_candidate_manifest"]["sha256"] = "0" * 64
        self._assert_rejected(candidate)

        projected = copy.deepcopy(self.manifest)
        projected["git"]["tree"] = "0" * 40
        self._assert_rejected(projected)

    def test_forged_commit_profile_and_freeze_identity_are_rejected(self) -> None:
        forged = copy.deepcopy(self.manifest)
        candidate_path = self.root / "clean_candidate_manifest.json"
        candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
        candidate["source"]["git_commit"] = "0" * 40
        candidate_path.write_text(json.dumps(candidate), encoding="utf-8")
        forged["git"]["commit"] = "0" * 40
        forged["git"]["clean_candidate_manifest"]["sha256"] = _sha256(candidate_path)
        self._assert_rejected(forged)

        candidate["source"]["git_commit"] = self.manifest["git"]["commit"]
        candidate["build"]["profile_sha256"] = "4" * 64
        candidate_path.write_text(json.dumps(candidate), encoding="utf-8")
        fabricated_profile = copy.deepcopy(self.manifest)
        fabricated_profile["git"]["clean_candidate_manifest"]["sha256"] = _sha256(
            candidate_path
        )
        self._assert_rejected(fabricated_profile)

        candidate["build"]["profile_sha256"] = _sha256(self.root / "build_profile.json")
        candidate["freeze_id"] = "not-a-uuid"
        candidate_path.write_text(json.dumps(candidate), encoding="utf-8")
        malformed_identity = copy.deepcopy(self.manifest)
        malformed_identity["git"]["clean_candidate_manifest"]["sha256"] = _sha256(
            candidate_path
        )
        self._assert_rejected(malformed_identity)

    def test_minimal_pseudo_commit_object_is_rejected(self) -> None:
        candidate = json.loads(
            (self.root / "clean_candidate_manifest.json").read_text(encoding="utf-8")
        )
        pseudo_commit = f"tree {candidate['source']['git_tree']}\n\n".encode("ascii")
        digest = hashlib.sha1(
            f"commit {len(pseudo_commit)}\0".encode("ascii") + pseudo_commit
        ).hexdigest()
        with self.assertRaises(ValueError):
            git_commit_tree_from_bytes(pseudo_commit, expected_commit=digest)

        misplaced = (
            f"tree {candidate['source']['git_tree']}\n"
            "author PIC Test <pic-test@example.invalid> 1 +0000\n"
            f"parent {'0' * 40}\n"
            "committer PIC Test <pic-test@example.invalid> 1 +0000\n\n"
        ).encode("ascii")
        misplaced_digest = hashlib.sha1(
            f"commit {len(misplaced)}\0".encode("ascii") + misplaced
        ).hexdigest()
        with self.assertRaises(ValueError):
            git_commit_tree_from_bytes(misplaced, expected_commit=misplaced_digest)

        trailing = (
            f"tree {candidate['source']['git_tree']}\n"
            "author PIC Test <pic-test@example.invalid> 1 +0000\n"
            "committer PIC Test <pic-test@example.invalid> 1 +0000\n"
            f"parent {'0' * 40}\n\n"
        ).encode("ascii")
        trailing_digest = hashlib.sha1(
            f"commit {len(trailing)}\0".encode("ascii") + trailing
        ).hexdigest()
        with self.assertRaises(ValueError):
            git_commit_tree_from_bytes(trailing, expected_commit=trailing_digest)

    def test_schema_timestamp_and_extra_fields_are_rejected(self) -> None:
        for value in (
            "definitely-not-a-date",
            "2026-05-30 12:00Z",
            "2026-W22-6T12:00:00Z",
        ):
            with self.subTest(created_utc=value):
                malformed_time = copy.deepcopy(self.manifest)
                malformed_time["created_utc"] = value
                self._assert_rejected(malformed_time)

        extra = copy.deepcopy(self.manifest)
        extra["review"]["unexpected"] = "not allowed"
        self._assert_rejected(extra)

        nonfinite = copy.deepcopy(self.manifest)
        nonfinite["resources"]["node_hours"] = float("nan")
        self._assert_rejected(nonfinite)

        frontier = copy.deepcopy(self.manifest)
        frontier["resources"]["platform"] = "Frontier"
        self._assert_rejected(frontier)

        relabeled = copy.deepcopy(self.manifest)
        relabeled["resources"]["platform"] = "Frontier "
        self._assert_rejected(relabeled)

        host_with_frontier_claim = copy.deepcopy(self.manifest)
        host_with_frontier_claim["resources"]["submission_id"] = "submission-1"
        self._assert_rejected(host_with_frontier_claim)

    def test_duplicate_json_keys_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            _load_object_bytes(b'{"schema_version": 1, "schema_version": 2}',
                               label="duplicate fixture")

    def test_frontier_binding_requires_matching_completed_reconciliation(self) -> None:
        pic_root = self.root / "frontier-pic"
        project_home_root = self.root / "frontier-project-home"
        manifest_path = pic_root / "manifests" / "campaign" / "manifest.json"
        artifact_dir = pic_root / "runs" / "campaign"
        manifest_path.parent.mkdir(parents=True)
        artifact_dir.mkdir(parents=True)
        analysis_dir = artifact_dir / "analysis"
        analysis_dir.mkdir()
        snapshot_analysis_dir = manifest_path.parent / "snapshot" / "analysis"
        snapshot_analysis_dir.mkdir(parents=True)
        project_home_root.mkdir()
        analyzer_path = snapshot_analysis_dir / "000-analysis.py"
        helper_path = snapshot_analysis_dir / "frontier_f1_structured_artifacts.py"
        analyzer_path.write_text("pass\n", encoding="utf-8")
        helper_path.write_text("pass\n", encoding="utf-8")
        analyzer_path.chmod(0o444)
        helper_path.chmod(0o444)
        manifest_path.write_text(
            json.dumps(
                {
                    "registered_science_authorization_id": "f1-clean-gyro-v1",
                    "snapshot_files": [
                        {
                            "role": "analysis-script-000",
                            "path": str(analyzer_path),
                            "sha256": _sha256(analyzer_path),
                        },
                        {
                            "role": "analysis-script-001",
                            "path": str(helper_path),
                            "sha256": _sha256(helper_path),
                        },
                    ],
                    "artifact_dir": str(artifact_dir),
                    "submission_id": "submission-1",
                }
            )
            + "\n",
            encoding="utf-8",
        )
        manifest_path.chmod(0o444)
        manifest_sha256 = _sha256(manifest_path)
        inventory_path = artifact_dir / "artifact_inventory.json"
        inventory_path.write_text('{"files":[],"schema_version":1}\n', encoding="utf-8")
        inventory_path.chmod(0o444)
        result_path = analysis_dir / "analysis.json"
        result_path.write_text('{"schema_version":1,"status":"pass"}\n', encoding="utf-8")
        result_path.chmod(0o444)
        receipt_path = analysis_dir / "offline_analysis_receipt.json"
        receipt_path.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "runner": {
                        "python": "/opt/cray/pe/python/3.11.7/bin/python3",
                        "flags": ["-I", "-B"],
                    },
                    "analyzer": {
                        "path": "000-analysis.py",
                        "sha256": _sha256(analyzer_path),
                    },
                    "support_modules": [
                        {
                            "path": "frontier_f1_structured_artifacts.py",
                            "sha256": _sha256(helper_path),
                        }
                    ],
                    "artifact_inventory": {
                        "path": "artifact_inventory.json",
                        "sha256": _sha256(inventory_path),
                    },
                    "analysis_result": {
                        "path": "analysis/analysis.json",
                        "sha256": _sha256(result_path),
                    },
                }
            )
            + "\n",
            encoding="utf-8",
        )
        receipt_path.chmod(0o444)
        candidate_sha256 = "a" * 64
        control_plane_version = "b" * 64
        resources = {
            "platform": "Frontier",
            "submission_id": "submission-1",
            "reservation_id": "reservation-1",
            "job_id": "12345",
            "pre_submit_manifest_path": str(manifest_path),
            "pre_submit_manifest_sha256": manifest_sha256,
            "run_artifact_dir": str(artifact_dir),
            "artifact_root": str(artifact_dir),
            "node_hours": 1.5,
            "registered_science_authorization_id": "f1-clean-gyro-v1",
            "artifact_inventory_path": str(inventory_path),
            "artifact_inventory_sha256": _sha256(inventory_path),
            "analysis_result_path": str(result_path),
            "analysis_result_sha256": _sha256(result_path),
            "offline_analysis_receipt_path": str(receipt_path),
            "offline_analysis_receipt_sha256": _sha256(receipt_path),
        }
        record = {
            "event_type": "reconciliation",
            "reconciled": True,
            "state": "COMPLETED",
            "submission_scope": "registered_science",
            "registered_science_authorization_id": resources[
                "registered_science_authorization_id"
            ],
            "submission_id": resources["submission_id"],
            "reservation_id": resources["reservation_id"],
            "job_id": resources["job_id"],
            "control_plane_version": control_plane_version,
            "clean_candidate_manifest_sha256": candidate_sha256,
            "manifest_path": str(manifest_path),
            "manifest_sha256": manifest_sha256,
            "artifact_dir": str(artifact_dir),
            "consumed_node_hours": resources["node_hours"],
        }
        genesis = {
            "event_type": "genesis",
            "state": "initialized",
            "control_plane_version": "historical-genesis-control-plane",
        }
        with patch(
            "tst.publication.pic_qualification_manifest.validate_mirrored_state",
            return_value=[genesis, record],
        ), patch(
            "tst.publication.pic_qualification_manifest._verify_frontier_offline_analysis"
        ) as recompute:
            _require_frontier_ledger_binding(
                {"resources": resources},
                candidate_sha256=candidate_sha256,
                control_plane_version=control_plane_version,
                authorized_pic_root=pic_root,
                authorized_project_home_root=project_home_root,
            )
            record["state"] = "FAILED"
            with self.assertRaises(ValueError):
                _require_frontier_ledger_binding(
                    {"resources": resources},
                    candidate_sha256=candidate_sha256,
                    control_plane_version=control_plane_version,
                    authorized_pic_root=pic_root,
                    authorized_project_home_root=project_home_root,
                )
            record["state"] = "COMPLETED"
            resources["artifact_root"] = str(self.root / "qualification-bundle")
            _require_frontier_ledger_binding(
                {"resources": resources},
                candidate_sha256=candidate_sha256,
                control_plane_version=control_plane_version,
                authorized_pic_root=pic_root,
                authorized_project_home_root=project_home_root,
            )
            resources["artifact_root"] = str(artifact_dir)
            resources["node_hours"] = 2.0
            with self.assertRaises(ValueError):
                _require_frontier_ledger_binding(
                    {"resources": resources},
                    candidate_sha256=candidate_sha256,
                    control_plane_version=control_plane_version,
                    authorized_pic_root=pic_root,
                    authorized_project_home_root=project_home_root,
                )
            resources["node_hours"] = 1.5
            resources["registered_science_authorization_id"] = "f1-clean-paper-coupling-v1"
            with self.assertRaises(ValueError):
                _require_frontier_ledger_binding(
                    {"resources": resources},
                    candidate_sha256=candidate_sha256,
                    control_plane_version=control_plane_version,
                    authorized_pic_root=pic_root,
                    authorized_project_home_root=project_home_root,
                )
            resources["registered_science_authorization_id"] = "f1-clean-gyro-v1"
            resources["artifact_inventory_path"] = str(
                artifact_dir / "analysis" / ".." / "artifact_inventory.json"
            )
            with self.assertRaises(ValueError):
                _require_frontier_ledger_binding(
                    {"resources": resources},
                    candidate_sha256=candidate_sha256,
                    control_plane_version=control_plane_version,
                    authorized_pic_root=pic_root,
                    authorized_project_home_root=project_home_root,
                )
            resources["artifact_inventory_path"] = str(inventory_path)
            resources["artifact_inventory_sha256"] = "0" * 64
            with self.assertRaises(ValueError):
                _require_frontier_ledger_binding(
                    {"resources": resources},
                    candidate_sha256=candidate_sha256,
                    control_plane_version=control_plane_version,
                    authorized_pic_root=pic_root,
                    authorized_project_home_root=project_home_root,
                )
            resources["artifact_inventory_sha256"] = _sha256(inventory_path)
            resources["analysis_result_sha256"] = "0" * 64
            with self.assertRaises(ValueError):
                _require_frontier_ledger_binding(
                    {"resources": resources},
                    candidate_sha256=candidate_sha256,
                    control_plane_version=control_plane_version,
                    authorized_pic_root=pic_root,
                    authorized_project_home_root=project_home_root,
                )
            resources["analysis_result_sha256"] = _sha256(result_path)
            receipt_resource_path = resources.pop("offline_analysis_receipt_path")
            with self.assertRaises(ValueError):
                _require_frontier_ledger_binding(
                    {"resources": resources},
                    candidate_sha256=candidate_sha256,
                    control_plane_version=control_plane_version,
                    authorized_pic_root=pic_root,
                    authorized_project_home_root=project_home_root,
                )
            resources["offline_analysis_receipt_path"] = receipt_resource_path
            resources["offline_analysis_receipt_path"] = (
                f"{analysis_dir}/./offline_analysis_receipt.json"
            )
            with self.assertRaises(ValueError):
                _require_frontier_ledger_binding(
                    {"resources": resources},
                    candidate_sha256=candidate_sha256,
                    control_plane_version=control_plane_version,
                    authorized_pic_root=pic_root,
                    authorized_project_home_root=project_home_root,
                )
            resources["offline_analysis_receipt_path"] = receipt_resource_path
            resources["offline_analysis_receipt_sha256"] = "0" * 64
            with self.assertRaises(ValueError):
                _require_frontier_ledger_binding(
                    {"resources": resources},
                    candidate_sha256=candidate_sha256,
                    control_plane_version=control_plane_version,
                    authorized_pic_root=pic_root,
                    authorized_project_home_root=project_home_root,
                )
            resources["offline_analysis_receipt_sha256"] = _sha256(receipt_path)

            original_receipt = receipt_path.read_bytes()
            receipt_path.chmod(0o644)
            receipt_path.write_text('{"schema_version":1}\n', encoding="utf-8")
            receipt_path.chmod(0o444)
            resources["offline_analysis_receipt_sha256"] = _sha256(receipt_path)
            with self.assertRaisesRegex(ValueError, "receipt differs"):
                _require_frontier_ledger_binding(
                    {"resources": resources},
                    candidate_sha256=candidate_sha256,
                    control_plane_version=control_plane_version,
                    authorized_pic_root=pic_root,
                    authorized_project_home_root=project_home_root,
                )
            receipt_path.chmod(0o644)
            receipt_path.write_bytes(original_receipt)
            receipt_path.chmod(0o444)
            resources["offline_analysis_receipt_sha256"] = _sha256(receipt_path)

            for source_path in (analyzer_path, helper_path):
                with self.subTest(replaced_source=source_path.name):
                    detached = source_path.with_name(source_path.name + ".detached")

                    def replace_source(*_: object, **__: object) -> None:
                        source_path.rename(detached)
                        source_path.write_text("pass\n", encoding="utf-8")
                        source_path.chmod(0o444)

                    recompute.side_effect = replace_source
                    try:
                        with self.assertRaisesRegex(ValueError, "path changed"):
                            _require_frontier_ledger_binding(
                                {"resources": resources},
                                candidate_sha256=candidate_sha256,
                                control_plane_version=control_plane_version,
                                authorized_pic_root=pic_root,
                                authorized_project_home_root=project_home_root,
                            )
                    finally:
                        source_path.unlink(missing_ok=True)
                        detached.rename(source_path)
                    recompute.side_effect = None

            for evidence_path in (inventory_path, result_path, receipt_path):
                with self.subTest(replaced_evidence=evidence_path.name):
                    detached = evidence_path.with_name(evidence_path.name + ".detached")
                    original = evidence_path.read_bytes()

                    def replace_evidence(*_: object, **__: object) -> None:
                        evidence_path.rename(detached)
                        evidence_path.write_bytes(original)
                        evidence_path.chmod(0o444)

                    recompute.side_effect = replace_evidence
                    try:
                        with self.assertRaisesRegex(ValueError, "path changed"):
                            _require_frontier_ledger_binding(
                                {"resources": resources},
                                candidate_sha256=candidate_sha256,
                                control_plane_version=control_plane_version,
                                authorized_pic_root=pic_root,
                                authorized_project_home_root=project_home_root,
                            )
                    finally:
                        evidence_path.unlink(missing_ok=True)
                        detached.rename(evidence_path)
                    recompute.side_effect = None

            detached_receipt = receipt_path.with_name("offline_analysis_receipt.detached")

            def remove_receipt(*_: object, **__: object) -> None:
                receipt_path.rename(detached_receipt)

            recompute.side_effect = remove_receipt
            try:
                with self.assertRaisesRegex(ValueError, "path changed"):
                    _require_frontier_ledger_binding(
                        {"resources": resources},
                        candidate_sha256=candidate_sha256,
                        control_plane_version=control_plane_version,
                        authorized_pic_root=pic_root,
                        authorized_project_home_root=project_home_root,
                    )
            finally:
                detached_receipt.rename(receipt_path)
            recompute.side_effect = None

            detached_artifact_dir = artifact_dir.with_name("campaign-detached")
            replacement_artifact_dir = artifact_dir.with_name("campaign-replacement")
            replacement_artifact_dir.mkdir()

            def replace_artifact_dir(*_: object, **__: object) -> None:
                artifact_dir.rename(detached_artifact_dir)
                replacement_artifact_dir.rename(artifact_dir)

            recompute.side_effect = replace_artifact_dir
            try:
                with self.assertRaisesRegex(ValueError, "path changed"):
                    _require_frontier_ledger_binding(
                        {"resources": resources},
                        candidate_sha256=candidate_sha256,
                        control_plane_version=control_plane_version,
                        authorized_pic_root=pic_root,
                        authorized_project_home_root=project_home_root,
                    )
            finally:
                artifact_dir.rmdir()
                detached_artifact_dir.rename(artifact_dir)
            recompute.side_effect = None

            original_manifest = manifest_path.read_bytes()
            for replacement in (original_manifest, b'{"changed":true}\n'):
                with self.subTest(replaced_pre_submit_manifest=replacement):
                    detached_manifest = manifest_path.with_name("manifest.detached")

                    def replace_manifest(*_: object, **__: object) -> None:
                        manifest_path.rename(detached_manifest)
                        manifest_path.write_bytes(replacement)
                        manifest_path.chmod(0o444)

                    recompute.side_effect = replace_manifest
                    try:
                        with self.assertRaisesRegex(ValueError, "path changed"):
                            _require_frontier_ledger_binding(
                                {"resources": resources},
                                candidate_sha256=candidate_sha256,
                                control_plane_version=control_plane_version,
                                authorized_pic_root=pic_root,
                                authorized_project_home_root=project_home_root,
                            )
                    finally:
                        manifest_path.unlink(missing_ok=True)
                        detached_manifest.rename(manifest_path)
                    recompute.side_effect = None

            detached_transplanted_pic_root = self.root / "frontier-pic-transplanted-detached"

            def transplant_pic_root(*_: object, **__: object) -> None:
                pic_root.rename(detached_transplanted_pic_root)
                pic_root.mkdir()
                (detached_transplanted_pic_root / "runs").rename(pic_root / "runs")
                (detached_transplanted_pic_root / "manifests").rename(
                    pic_root / "manifests"
                )

            recompute.side_effect = transplant_pic_root
            try:
                with self.assertRaisesRegex(ValueError, "ancestry changed"):
                    _require_frontier_ledger_binding(
                        {"resources": resources},
                        candidate_sha256=candidate_sha256,
                        control_plane_version=control_plane_version,
                        authorized_pic_root=pic_root,
                        authorized_project_home_root=project_home_root,
                    )
            finally:
                (pic_root / "runs").rename(detached_transplanted_pic_root / "runs")
                (pic_root / "manifests").rename(
                    detached_transplanted_pic_root / "manifests"
                )
                pic_root.rmdir()
                detached_transplanted_pic_root.rename(pic_root)
            recompute.side_effect = None

            detached_pre_pin_pic_root = self.root / "frontier-pic-pre-pin-detached"
            real_open_pinned = (
                qualification_manifest._open_pinned_read_only_regular_file_at
            )
            transplanted = False

            def transplant_pic_root_after_manifest_open(
                *args: object, **kwargs: object
            ) -> tuple[int, bytes]:
                nonlocal transplanted
                result = real_open_pinned(*args, **kwargs)
                if (
                    not transplanted
                    and kwargs.get("label") == "Frontier pre-submit manifest"
                ):
                    pic_root.rename(detached_pre_pin_pic_root)
                    pic_root.mkdir()
                    for child in list(detached_pre_pin_pic_root.iterdir()):
                        child.rename(pic_root / child.name)
                    transplanted = True
                return result

            try:
                with patch(
                    "tst.publication.pic_qualification_manifest."
                    "_open_pinned_read_only_regular_file_at",
                    side_effect=transplant_pic_root_after_manifest_open,
                ):
                    with self.assertRaisesRegex(ValueError, "ancestry changed"):
                        _require_frontier_ledger_binding(
                            {"resources": resources},
                            candidate_sha256=candidate_sha256,
                            control_plane_version=control_plane_version,
                            authorized_pic_root=pic_root,
                            authorized_project_home_root=project_home_root,
                        )
            finally:
                for child in list(pic_root.iterdir()):
                    child.rename(detached_pre_pin_pic_root / child.name)
                pic_root.rmdir()
                detached_pre_pin_pic_root.rename(pic_root)

            runs_dir = pic_root / "runs"
            detached_runs_dir = pic_root / "runs-detached"
            runs_dir.rename(detached_runs_dir)
            runs_dir.symlink_to(detached_runs_dir, target_is_directory=True)
            try:
                with self.assertRaises((OSError, ValueError)):
                    _require_frontier_ledger_binding(
                        {"resources": resources},
                        candidate_sha256=candidate_sha256,
                        control_plane_version=control_plane_version,
                        authorized_pic_root=pic_root,
                        authorized_project_home_root=project_home_root,
                    )
            finally:
                runs_dir.unlink()
                detached_runs_dir.rename(runs_dir)

            snapshot_dir = manifest_path.parent / "snapshot"
            detached_snapshot_dir = manifest_path.parent / "snapshot-detached"
            snapshot_dir.rename(detached_snapshot_dir)
            snapshot_dir.symlink_to(detached_snapshot_dir, target_is_directory=True)
            try:
                with self.assertRaises((OSError, ValueError)):
                    _require_frontier_ledger_binding(
                        {"resources": resources},
                        candidate_sha256=candidate_sha256,
                        control_plane_version=control_plane_version,
                        authorized_pic_root=pic_root,
                        authorized_project_home_root=project_home_root,
                    )
            finally:
                snapshot_dir.unlink()
                detached_snapshot_dir.rename(snapshot_dir)

            detached_pic_root = self.root / "frontier-pic-detached"
            pic_root.rename(detached_pic_root)
            pic_root.symlink_to(detached_pic_root, target_is_directory=True)
            try:
                with self.assertRaises((OSError, ValueError)):
                    _require_frontier_ledger_binding(
                        {"resources": resources},
                        candidate_sha256=candidate_sha256,
                        control_plane_version=control_plane_version,
                        authorized_pic_root=pic_root,
                        authorized_project_home_root=project_home_root,
                    )
            finally:
                pic_root.unlink()
                detached_pic_root.rename(pic_root)

            _require_frontier_ledger_binding(
                {"resources": resources},
                candidate_sha256=candidate_sha256,
                control_plane_version=control_plane_version,
                authorized_pic_root=pic_root,
                authorized_project_home_root=project_home_root,
            )
            kwargs = recompute.call_args.kwargs
            self.assertIsInstance(recompute.call_args.args[0], int)
            self.assertIsInstance(kwargs["helper_fd"], int)
            self.assertIsInstance(kwargs["artifact_dir_fd"], int)
            self.assertEqual(kwargs["artifact_dir"], artifact_dir)
            self.assertEqual(
                kwargs["artifact_inventory_sha256"], _sha256(inventory_path)
            )
            self.assertEqual(kwargs["result_sha256"], _sha256(result_path))
        with patch(
            "tst.publication.pic_qualification_manifest.validate_mirrored_state",
            return_value=[record],
        ):
            with self.assertRaises(ValueError):
                _require_frontier_ledger_binding(
                    {"resources": resources},
                    candidate_sha256=candidate_sha256,
                    control_plane_version=control_plane_version,
                    authorized_pic_root=pic_root,
                    authorized_project_home_root=project_home_root,
                )

    def test_frontier_offline_analysis_recompute_uses_trusted_snapshot(self) -> None:
        completed = subprocess.CompletedProcess([], 0, stdout=b"", stderr=b"")
        with patch(
            "tst.publication.pic_qualification_manifest.subprocess.run",
            return_value=completed,
        ) as run:
            _verify_frontier_offline_analysis(
                10,
                helper_fd=11,
                artifact_dir_fd=12,
                artifact_dir=Path("/runs/f1"),
                artifact_inventory_sha256="b" * 64,
                result_sha256="a" * 64,
            )
        run.assert_called_once_with(
            [
                "/opt/cray/pe/python/3.11.7/bin/python3",
                "-I",
                "-B",
                "/proc/self/fd/10",
                "--artifact-dir",
                "/runs/f1",
                "--artifact-dir-fd",
                "12",
                "--verify-artifact-inventory-sha256",
                "b" * 64,
                "--verify-result-sha256",
                "a" * 64,
            ],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=300,
            pass_fds=(10, 11, 12),
            env={
                "HOME": "/",
                "LANG": "C",
                "LC_ALL": "C",
                "PATH": "/usr/bin:/bin",
                "PIC_F1_ANALYSIS_HELPER_FD": "11",
            },
            cwd="/",
        )

    def test_candidate_internal_path_alias_is_rejected(self) -> None:
        candidate_path = self.root / "clean_candidate_manifest.json"
        candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
        candidate["source"]["commit_path"] = (
            "/original/03a7bd9a-7d4c-4e37-a12b-46de3817eff2/./source.commit"
        )
        candidate_path.write_text(json.dumps(candidate), encoding="utf-8")
        manifest = copy.deepcopy(self.manifest)
        manifest["git"]["clean_candidate_manifest"]["sha256"] = _sha256(candidate_path)
        self._assert_rejected(manifest)

        candidate = json.loads(
            (self.root / "clean_candidate_manifest.json").read_text(encoding="utf-8")
        )
        for record in (
            candidate["source"],
            candidate["build"],
            *candidate["source"]["submodules"],
        ):
            for key in (
                "archive_path",
                "commit_path",
                "profile_path",
                "profile_receipt_path",
                "executable_path",
            ):
                if key in record:
                    record[key] = "/" + record[key].replace("/./", "/")
        candidate_path.write_text(json.dumps(candidate), encoding="utf-8")
        doubled = copy.deepcopy(self.manifest)
        doubled["git"]["clean_candidate_manifest"]["sha256"] = _sha256(candidate_path)
        self._assert_rejected(doubled)

        candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
        original = "/original/03a7bd9a-7d4c-4e37-a12b-46de3817eff2"
        alias = "/original/segment/../03a7bd9a-7d4c-4e37-a12b-46de3817eff2"
        for record in (
            candidate["source"],
            candidate["build"],
            *candidate["source"]["submodules"],
        ):
            for key in (
                "archive_path",
                "commit_path",
                "profile_path",
                "profile_receipt_path",
                "executable_path",
            ):
                if key in record:
                    value = record[key][1:] if record[key].startswith("//") else record[key]
                    record[key] = value.replace(original, alias)
        candidate_path.write_text(json.dumps(candidate), encoding="utf-8")
        parent_alias = copy.deepcopy(self.manifest)
        parent_alias["git"]["clean_candidate_manifest"]["sha256"] = _sha256(candidate_path)
        self._assert_rejected(parent_alias)

    def test_blank_clean_candidate_build_metadata_is_rejected(self) -> None:
        candidate_path = self.root / "clean_candidate_manifest.json"
        profile_path = self.root / "build_profile.json"
        candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
        profile = json.loads(profile_path.read_text(encoding="utf-8"))
        for field in ("profile_id", "toolchain"):
            candidate["build"][field] = ""
            profile[field] = ""
        profile_path.write_text(json.dumps(profile), encoding="utf-8")
        candidate["build"]["profile_sha256"] = _sha256(profile_path)
        candidate_path.write_text(json.dumps(candidate), encoding="utf-8")
        manifest = copy.deepcopy(self.manifest)
        manifest["git"]["clean_candidate_manifest"]["sha256"] = _sha256(candidate_path)
        manifest["git"]["clean_candidate_build_profile"]["sha256"] = _sha256(profile_path)
        self.live_authorization["clean_candidate_manifest_sha256"] = _sha256(
            candidate_path
        )
        self._assert_rejected(manifest)

    def test_escaped_missing_and_checksum_mismatched_files_are_rejected(
        self,
    ) -> None:
        escaped = copy.deepcopy(self.manifest)
        escaped["artifacts"][0]["path"] = "../outside.json"
        self._assert_rejected(escaped)

        missing = copy.deepcopy(self.manifest)
        missing["executable"]["modules"]["path"] = "missing-modules.txt"
        self._assert_rejected(missing)

        mismatched = copy.deepcopy(self.manifest)
        mismatched["artifacts"][0]["sha256"] = "0" * 64
        self._assert_rejected(mismatched)

        aliased = copy.deepcopy(self.manifest)
        aliased["git"]["source_archive"]["path"] = "./source.tar"
        self._assert_rejected(aliased)

        (self.root / "metrics-link.json").symlink_to("metrics.json")
        symlinked = copy.deepcopy(self.manifest)
        symlinked["artifacts"][0]["path"] = "metrics-link.json"
        self._assert_rejected(symlinked)

    def test_intermediate_directory_swap_is_rejected(self) -> None:
        from tst.publication import pic_qualification_manifest as validator

        inside = self.root / "inside"
        outside = self.root / "outside"
        inside.mkdir()
        outside.mkdir()
        (inside / "metrics.json").write_text('{"relative_error": 0.0}\n', encoding="utf-8")
        (outside / "metrics.json").write_text('{"relative_error": 0.0}\n', encoding="utf-8")
        manifest = copy.deepcopy(self.manifest)
        manifest["artifacts"][0] = {
            "path": "inside/metrics.json",
            "sha256": _sha256(inside / "metrics.json"),
        }
        original = validator._open_directory_component_at
        swapped = False

        def swap(directory_fd: int, name: str, *, label: str) -> int:
            nonlocal swapped
            if name == "inside" and not swapped:
                swapped = True
                inside.rename(self.root / "inside-original")
                inside.symlink_to(outside, target_is_directory=True)
            return original(directory_fd, name, label=label)

        with patch.object(validator, "_open_directory_component_at", side_effect=swap):
            self._assert_rejected(manifest)
        self.assertTrue(swapped)

    def test_artifact_root_component_swap_is_rejected(self) -> None:
        from tst.publication import pic_qualification_manifest as validator

        original_root = self.root.with_name(self.root.name + "-original")
        original = validator._open_directory_component_at
        swapped = False

        def swap(directory_fd: int, name: str, *, label: str) -> int:
            nonlocal swapped
            if name == self.root.name and not swapped:
                swapped = True
                self.root.rename(original_root)
                self.root.symlink_to(original_root, target_is_directory=True)
            return original(directory_fd, name, label=label)

        try:
            with patch.object(validator, "_open_directory_component_at", side_effect=swap):
                self._assert_rejected(self.manifest)
            self.assertTrue(swapped)
        finally:
            if self.root.is_symlink():
                self.root.unlink()
            if original_root.exists():
                original_root.rename(self.root)

    def test_freeze_writes_canonical_file_once(self) -> None:
        source = self.root / "prepared.json"
        output = self.root / "frozen.json"
        source.write_text(json.dumps(self.manifest), encoding="utf-8")
        freeze_qualification_manifest(source, output)
        self.assertEqual(json.loads(output.read_text(encoding="utf-8")),
                         self.manifest)
        self.assertFalse(output.stat().st_mode & 0o222)
        with self.assertRaises(FileExistsError):
            freeze_qualification_manifest(source, output)

    def test_freeze_chmod_failure_leaves_no_mutable_output(self) -> None:
        source = self.root / "prepared.json"
        output = self.root / "frozen.json"
        source.write_text(json.dumps(self.manifest), encoding="utf-8")
        with patch("tst.publication.pic_qualification_manifest.os.fchmod",
                   side_effect=OSError("chmod failed")):
            with self.assertRaises(OSError):
                freeze_qualification_manifest(source, output)
        self.assertFalse(output.exists())
        self.assertFalse(list(self.root.glob(".frozen.json.tmp-*")))

    def test_freeze_rejects_symlinked_output_parent(self) -> None:
        source = self.root / "prepared.json"
        outside = self.root / "outside"
        alias = self.root / "alias"
        outside.mkdir()
        alias.symlink_to(outside, target_is_directory=True)
        source.write_text(json.dumps(self.manifest), encoding="utf-8")
        with self.assertRaises(ValueError):
            freeze_qualification_manifest(source, alias / "frozen.json")
        self.assertEqual(list(outside.iterdir()), [])

    def test_freeze_directory_fsync_failure_rolls_back_output(self) -> None:
        source = self.root / "prepared.json"
        output = self.root / "frozen.json"
        source.write_text(json.dumps(self.manifest), encoding="utf-8")
        real_fsync = os.fsync

        def fail_directory_fsync(descriptor: int) -> None:
            if stat.S_ISDIR(os.fstat(descriptor).st_mode):
                raise OSError("directory fsync failed")
            real_fsync(descriptor)

        with patch(
            "tst.publication.pic_qualification_manifest.os.fsync",
            side_effect=fail_directory_fsync,
        ):
            with self.assertRaises(RuntimeError):
                freeze_qualification_manifest(source, output)
        self.assertFalse(output.exists())
        self.assertFalse(list(self.root.glob(".frozen.json.tmp-*")))

    def test_freeze_single_directory_fsync_failure_durably_rolls_back_output(self) -> None:
        source = self.root / "prepared.json"
        output = self.root / "frozen.json"
        source.write_text(json.dumps(self.manifest), encoding="utf-8")
        real_fsync = os.fsync
        failed = False

        def fail_first_directory_fsync(descriptor: int) -> None:
            nonlocal failed
            if stat.S_ISDIR(os.fstat(descriptor).st_mode) and not failed:
                failed = True
                raise OSError("directory fsync failed")
            real_fsync(descriptor)

        with patch(
            "tst.publication.pic_qualification_manifest.os.fsync",
            side_effect=fail_first_directory_fsync,
        ):
            with self.assertRaises(OSError):
                freeze_qualification_manifest(source, output)
        self.assertTrue(failed)
        self.assertFalse(output.exists())
        self.assertFalse(list(self.root.glob(".frozen.json.tmp-*")))


if __name__ == "__main__":
    unittest.main()
