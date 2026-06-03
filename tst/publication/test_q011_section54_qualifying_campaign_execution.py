#!/usr/bin/env python3
"""Focused tests for immutable Q-011 Section 5.4 campaign-plan materialization."""

from __future__ import annotations

import copy
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import tempfile
from typing import Iterator
import unittest
import uuid
from unittest.mock import patch

from tst.publication import q011_section54_pressure_selection as selection
from tst.publication import q011_section54_qualifying_campaign_execution as execution


_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _put(path: Path, payload: bytes, mode: int = 0o444) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    path.chmod(mode)
    return path


def _make_writable(root: Path) -> None:
    if not root.exists():
        return
    for path in [root, *root.rglob("*")]:
        if not path.is_symlink():
            path.chmod(path.stat().st_mode | stat.S_IWUSR)


def _receipt(*, case_id: str = "ps_p0_0p10", problem_ps_p0: float = 0.1) -> dict[str, object]:
    return {
        "schema_version": 1,
        "record_type": selection.RECORD_TYPE,
        "selection_method": selection.SELECTION_METHOD,
        "published_pressure_pilot_receipt": {
            "path": "/fixture/replaced-by-fixture.json",
            "sha256": "0" * 64,
        },
        "pilot_bundle_manifest_sha256": "a" * 64,
        "aggregate_pilot_analysis_sha256": "b" * 64,
        "case_descriptors": [
            {
                "case_id": registered_id,
                "problem_ps_p0": registered_ps_p0,
                "descriptor_sha256": character * 64,
            }
            for (registered_id, registered_ps_p0), character in zip(
                selection.REGISTERED_CASES, "cdef"
            )
        ],
        "selected_case": {
            "case_id": case_id,
            "problem_ps_p0": problem_ps_p0,
        },
        "reviewer_identity": "Focused Test Human Reviewer",
        "reviewed_utc": "2026-06-02T12:34:56Z",
        "rationale": "Focused test-only human pressure selection.",
    }


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _git(source: Path, *arguments: str) -> str:
    return subprocess.check_output(
        ["/usr/bin/git", "-C", str(source), *arguments],
        text=True,
    ).strip()


def _pressure_publication(root: Path) -> tuple[dict[str, str], dict[str, object]]:
    cases = [
        {
            "case_id": case_id,
            "artifact_dir": str(root / "runs" / case_id),
            "descriptor_path": "pressure_pilot_case_descriptor.json",
            "descriptor_sha256": character * 64,
            "artifact_inventory_sha256": "9" * 64,
            "runtime_artifacts": [],
        }
        for (case_id, _), character in zip(selection.REGISTERED_CASES, "cdef")
    ]
    payload = _json_bytes(
        {
            "aggregate_bundle": {
                "path": str(root / "publication" / "aggregate"),
                "manifest_sha256": "a" * 64,
            },
            "aggregate_analysis": {
                "path": str(root / "publication" / "aggregate-analysis.json"),
                "sha256": "b" * 64,
            },
            "raw_cases": cases,
        }
    )
    receipt = _put(root / "publication" / "pressure-pilot-receipt.json", payload)
    digest = _sha256(payload)
    return {"path": str(receipt), "sha256": digest}, {
        "receipt_sha256": digest,
        "manifest_sha256": "a" * 64,
        "analysis_result_sha256": "b" * 64,
        "status": "passed",
    }


@contextmanager
def _fixture(
    *,
    receipt: dict[str, object] | None = None,
    executable_payload: bytes | None = None,
) -> Iterator[dict[str, Path]]:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        orion = root / "orion"
        orion.mkdir()
        common = execution._load_control_plane_common()
        authorized_source = root / "source-authorized"
        authorized_source.mkdir()
        subprocess.run(["/usr/bin/git", "init", str(authorized_source)], check=True, capture_output=True)
        subprocess.run(
            ["/usr/bin/git", "-C", str(authorized_source), "config", "user.email", "focused@example.invalid"],
            check=True,
        )
        subprocess.run(
            ["/usr/bin/git", "-C", str(authorized_source), "config", "user.name", "Focused Test"],
            check=True,
        )
        prepared_sources = {
            relative: (execution.REPO_ROOT / relative).read_bytes()
            for relative in sorted(
                {
                    *common.Q011_SECTION54_ARCHIVE_SOURCE_PATHS.values(),
                    *common.Q011_SECTION54_HELPER_SOURCES,
                }
            )
        }
        inventory = {
            "schema_version": 1,
            "paper_decks": [
                {"path": path, "sha256": _sha256(prepared_sources[path])}
                for path in sorted(prepared_sources)
                if path in common.PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS
            ],
            "analyzers": [
                {"path": path, "sha256": _sha256(prepared_sources[path])}
                for path in sorted(prepared_sources)
                if path.startswith("tst/publication/analyze_")
                and path.endswith(".py")
                and "/" not in path[len("tst/publication/") :]
            ],
        }
        prepared_sources[common.PREPARED_ARTIFACT_INVENTORY_PATH] = _json_bytes(inventory)
        for relative, payload in prepared_sources.items():
            path = authorized_source / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
        subprocess.run(["/usr/bin/git", "-C", str(authorized_source), "add", "."], check=True)
        subprocess.run(
            ["/usr/bin/git", "-C", str(authorized_source), "commit", "-m", "focused fixture"],
            check=True,
            capture_output=True,
        )
        git_commit = _git(authorized_source, "rev-parse", "HEAD")
        git_tree = _git(authorized_source, "rev-parse", "HEAD^{tree}")
        freeze_id = str(uuid.uuid4())
        candidate_root = orion / "clean_candidates" / freeze_id
        candidate_root.mkdir(parents=True)
        archive = candidate_root / "source.tar"
        subprocess.run(
            ["/usr/bin/git", "-C", str(authorized_source), "archive", "--format=tar", f"--output={archive}", git_commit],
            check=True,
        )
        archive.chmod(0o444)
        commit = _put(
            candidate_root / "source.commit",
            subprocess.check_output(["/usr/bin/git", "-C", str(authorized_source), "cat-file", "commit", git_commit]),
        )
        executable = _put(
            candidate_root / "athena",
            Path("/bin/true").read_bytes() if executable_payload is None else executable_payload,
            0o555,
        )
        archive_payload = archive.read_bytes()
        commit_payload = commit.read_bytes()
        source_bundle = common.source_bundle_sha256(
            _sha256(archive_payload), _sha256(commit_payload), []
        )
        provenance_payloads = {
            "configure_log": b"focused configure log\n",
            "build_log": b"focused build log\n",
            "cmake_cache": b"focused cache\n",
            "module_list": common.production_module_list_bytes(),
            "toolchain": (common.PRODUCTION_TOOLCHAIN_DESCRIPTION + "\n").encode("utf-8"),
            "build_invocations": _json_bytes(
                common.production_build_invocations(
                    authorized_pic_root=orion,
                    git_commit=git_commit,
                    profile_id=common.PRODUCTION_BUILD_PROFILE,
                )
            ),
            "git_status_preconfigure": b"",
            "git_status": b"",
            "submodule_status": b"",
            "environment_allowlist": common.production_environment_allowlist_bytes(),
            "build_environment": _json_bytes(common.PRODUCTION_BUILD_ENVIRONMENT),
        }
        provenance_paths = common._documented_build_provenance_paths(
            authorized_pic_root=orion,
            git_commit=git_commit,
            profile_id=common.PRODUCTION_BUILD_PROFILE,
        )
        provenance_records = {}
        for label, filename in common.BUILD_PROVENANCE_FILENAMES.items():
            payload = provenance_payloads[label]
            _put(provenance_paths[label], payload)
            _put(candidate_root / "build_provenance" / filename, payload)
            provenance_records[label] = {
                "path": str(provenance_paths[label]),
                "sha256": _sha256(payload),
            }
        profile_value = {
            "schema_version": 3,
            "profile_id": common.PRODUCTION_BUILD_PROFILE,
            "authorized_source_root": str(authorized_source),
            "fresh_source_root": str(
                orion / "build" / git_commit[:12] / common.PRODUCTION_BUILD_PROFILE / "source"
            ),
            "git_commit": git_commit,
            "git_tree": git_tree,
            "source_archive_sha256": _sha256(archive_payload),
            "source_commit_sha256": _sha256(commit_payload),
            "source_bundle_sha256": source_bundle,
            "toolchain": common.PRODUCTION_TOOLCHAIN_DESCRIPTION,
            "build_invocations_sha256": _sha256(provenance_payloads["build_invocations"]),
            "executable_sha256": _sha256(executable.read_bytes()),
            "provenance_inputs": provenance_records,
            "submodules": [],
        }
        profile = _put(candidate_root / "build_profile.json", _json_bytes(profile_value))
        installed_control_plane_files = {
            name: (
                execution.REPO_ROOT / "tst/publication/frontier_control_plane" / name
            ).read_bytes()
            for name in common.CONTROL_PLANE_FILES
        }
        installed_control_plane_records = [
            {"path": name, "sha256": _sha256(installed_control_plane_files[name])}
            for name in common.CONTROL_PLANE_FILES
        ]
        control_plane_version = common.inventory_digest(installed_control_plane_records)
        installed_control_plane = orion / "control_plane" / control_plane_version
        for name, payload in installed_control_plane_files.items():
            _put(installed_control_plane / name, payload)
        _put(
            installed_control_plane / "inventory.json",
            _json_bytes(
                {
                    "schema_version": 1,
                    "version": control_plane_version,
                    "files": installed_control_plane_records,
                }
            ),
        )
        installed_control_plane.chmod(0o555)
        artifact_dir = orion / "bin" / git_commit[:12] / common.PRODUCTION_BUILD_PROFILE
        profile_receipt_value = {
            "schema_version": 1,
            "control_plane_version": control_plane_version,
            "profile_path": str(artifact_dir / "build_profile.json"),
            "profile_sha256": _sha256(profile.read_bytes()),
            "source_bundle_sha256": source_bundle,
            "fresh_source_root": profile_value["fresh_source_root"],
            "build_invocations_sha256": profile_value["build_invocations_sha256"],
            "git_status_preconfigure_sha256": provenance_records["git_status_preconfigure"]["sha256"],
            "git_status_sha256": provenance_records["git_status"]["sha256"],
            "configure_log_sha256": provenance_records["configure_log"]["sha256"],
            "build_log_sha256": provenance_records["build_log"]["sha256"],
            "executable_path": str(artifact_dir / "athena"),
            "executable_sha256": profile_value["executable_sha256"],
        }
        profile_receipt = _put(candidate_root / "profile_receipt.json", _json_bytes(profile_receipt_value))
        environment = installed_control_plane / "frontier_pic_environment.sh"
        manifest = {
            "schema_version": 4,
            "freeze_id": freeze_id,
            "created_utc": "2026-06-02T12:00:00Z",
            "prepared_artifacts": common.prepared_artifact_manifest_from_source_archive(
                archive_payload,
                inventory_path=common.PREPARED_ARTIFACT_INVENTORY_PATH,
            ),
            "source": {
                "archive_path": str(archive),
                "archive_sha256": _sha256(archive.read_bytes()),
                "commit_path": str(commit),
                "commit_sha256": _sha256(commit.read_bytes()),
                "source_bundle_sha256": source_bundle,
                "git_commit": git_commit,
                "git_tree": git_tree,
                "worktree_status": "clean",
                "submodule_status": "absent",
                "submodules": [],
            },
            "build": {
                "profile_id": common.PRODUCTION_BUILD_PROFILE,
                "profile_path": str(profile),
                "profile_sha256": _sha256(profile.read_bytes()),
                "profile_receipt_path": str(profile_receipt),
                "profile_receipt_sha256": _sha256(profile_receipt.read_bytes()),
                "source_archive_sha256": _sha256(archive.read_bytes()),
                "source_commit_sha256": _sha256(commit.read_bytes()),
                "source_bundle_sha256": source_bundle,
                "toolchain": common.PRODUCTION_TOOLCHAIN_DESCRIPTION,
                "build_invocations_sha256": profile_value["build_invocations_sha256"],
                "executable_path": str(executable),
                "executable_sha256": _sha256(executable.read_bytes()),
            },
        }
        candidate = _put(
            candidate_root / "clean_candidate_manifest.json",
            (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8"),
        )
        publication_binding, verified_publication = _pressure_publication(orion)
        receipt_value = copy.deepcopy(_receipt() if receipt is None else receipt)
        receipt_value["published_pressure_pilot_receipt"] = publication_binding
        pressure_receipt = _put(
            root / "human_pressure_selection_receipt.json",
            selection.canonical_json_bytes(receipt_value),
        )
        output_parent = root / "plans"
        output_parent.mkdir()
        with (
            patch.object(execution, "AUTHORIZED_ORION_ROOT", orion),
            patch.object(execution, "AUTHORIZED_SOURCE_ROOT", authorized_source),
            patch.object(
                selection.pressure_pilot_publisher,
                "verify_published_pressure_pilot_receipt",
                return_value=verified_publication,
            ),
        ):
            yield {
                "root": root,
                "orion": orion,
                "candidate": candidate,
                "executable": executable,
                "environment": environment,
                "pressure_receipt": pressure_receipt,
                "output_parent": output_parent,
            }
        for child in output_parent.iterdir():
            _make_writable(child)
            shutil.rmtree(child)
        _make_writable(installed_control_plane)


def _materialize(fixture: dict[str, Path], *, output_parent: Path | None = None) -> dict[str, object]:
    return execution.materialize_qualifying_campaign_plan(
        output_parent=fixture["output_parent"] if output_parent is None else output_parent,
        pressure_selection_receipt=fixture["pressure_receipt"],
        clean_candidate_manifest=fixture["candidate"],
        executable=fixture["executable"],
        environment_profile=fixture["environment"],
    )


@contextmanager
def _retained_planner(
    fixture: dict[str, Path],
) -> Iterator[tuple[dict[str, object], Path]]:
    planner_parent = fixture["orion"] / "plans"
    planner_parent.mkdir()
    result = _materialize(fixture, output_parent=planner_parent)
    root = Path(result["plan_root"])
    try:
        yield result, root
    finally:
        _make_writable(root)
        if root.exists():
            shutil.rmtree(root)
        planner_parent.rmdir()


def _planner_retention(
    fixture: dict[str, Path],
    result: dict[str, object],
    root: Path,
    attempt_id: str,
) -> dict[str, object]:
    return execution.materialize_planner_retention(
        planner_root=root,
        planner_inventory_sha256=str(result["inventory_sha256"]),
        attempt_id=attempt_id,
        authorized_pic_root=fixture["orion"],
    )


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


class Q011Section54QualifyingCampaignExecutionTests(unittest.TestCase):
    def test_planner_retention_overlay_is_exact_and_launch_prohibited(self) -> None:
        with _fixture() as fixture, _retained_planner(fixture) as (result, root):
            plan = _json(root / "campaign_plan.json")
            descriptor = _json(root / plan["baseline_attempt_descriptors"][8]["path"])
            contract = _json(root / descriptor["launch_contract"]["path"])
            attempt_root = descriptor["authorized_orion_attempt_root"]
            overlay = _planner_retention(
                fixture, result, root, str(descriptor["attempt_id"])
            )
            self.assertEqual(
                overlay,
                {
                    "schema_version": 1,
                    "retention_role": (
                        "q011_section54_deterministic_retained_attempt"
                    ),
                    "planner_root": str(root),
                    "planner_inventory_sha256": result["inventory_sha256"],
                    "planner_plan_id": result["plan_id"],
                    "planner_materialization_receipt": result[
                        "materialization_receipt"
                    ],
                    "attempt_id": descriptor["attempt_id"],
                    "authorized_orion_attempt_root": attempt_root,
                    "authorized_orion_raw_root": f"{attempt_root}/raw",
                    "argv": contract["argv"],
                },
            )
            self.assertEqual(
                set(overlay),
                {
                    "schema_version",
                    "retention_role",
                    "planner_root",
                    "planner_inventory_sha256",
                    "planner_plan_id",
                    "planner_materialization_receipt",
                    "attempt_id",
                    "authorized_orion_attempt_root",
                    "authorized_orion_raw_root",
                    "argv",
                },
            )
            self.assertFalse(contract["launch_authorized"])
            self.assertFalse(contract["scheduler_submission_authorized"])
            self.assertFalse(contract["live_policy_mutation_authorized"])

    def test_planner_retention_rejects_attempt_mismatch_and_nonbaseline_selection(
        self,
    ) -> None:
        with _fixture() as fixture, _retained_planner(fixture) as (result, root):
            plan = _json(root / "campaign_plan.json")
            carrier = _json(root / plan["restart_continuation_carrier"]["path"])
            for attempt_id in (
                "baseline-001-coarse_uniform_dx12-seed-99999999",
                str(carrier["carrier_id"]),
            ):
                with (
                    self.subTest(attempt_id=attempt_id),
                    self.assertRaisesRegex(
                        execution.CampaignPlanError,
                        "exactly one baseline planner descriptor",
                    ),
                ):
                    _planner_retention(fixture, result, root, attempt_id)

    def test_planner_retention_has_no_argv_or_destination_substitution_api(self) -> None:
        with _fixture() as fixture, _retained_planner(fixture) as (result, root):
            plan = _json(root / "campaign_plan.json")
            descriptor = _json(root / plan["baseline_attempt_descriptors"][0]["path"])
            arguments = {
                "planner_root": root,
                "planner_inventory_sha256": str(result["inventory_sha256"]),
                "attempt_id": str(descriptor["attempt_id"]),
                "authorized_pic_root": fixture["orion"],
            }
            for substitution in (
                {"argv": ["-d", "/tmp/operator-selected"]},
                {"authorized_orion_attempt_root": "/tmp/operator-selected"},
                {"authorized_orion_raw_root": "/tmp/operator-selected/raw"},
            ):
                with (
                    self.subTest(substitution=substitution),
                    self.assertRaises(TypeError),
                ):
                    execution.materialize_planner_retention(
                        **arguments,  # type: ignore[arg-type]
                        **substitution,
                    )

    def test_planner_retention_rejects_root_substitution_and_inventory_drift(
        self,
    ) -> None:
        with _fixture() as fixture, _retained_planner(fixture) as (result, root):
            plan = _json(root / "campaign_plan.json")
            descriptor = _json(root / plan["baseline_attempt_descriptors"][0]["path"])
            attempt_id = str(descriptor["attempt_id"])
            with self.assertRaises(execution.CampaignPlanError):
                execution.materialize_planner_retention(
                    planner_root=root,
                    planner_inventory_sha256="0" * 64,
                    attempt_id=attempt_id,
                    authorized_pic_root=fixture["orion"],
                )
            alias = fixture["orion"] / "operator-substituted-planner-root"
            alias.symlink_to(root, target_is_directory=True)
            try:
                with self.assertRaises(execution.CampaignPlanError):
                    execution.materialize_planner_retention(
                        planner_root=alias,
                        planner_inventory_sha256=str(result["inventory_sha256"]),
                        attempt_id=attempt_id,
                        authorized_pic_root=fixture["orion"],
                    )
            finally:
                alias.unlink()

    def test_materializes_exact_ordered_read_only_plan_and_one_restart_carrier(self) -> None:
        with _fixture() as fixture:
            result = _materialize(fixture)
            root = Path(result["plan_root"])
            plan = _json(root / "campaign_plan.json")
            self.assertEqual(result["baseline_attempt_count"], 24)
            self.assertEqual(result["restart_continuation_carrier_count"], 1)
            self.assertTrue(result["recursively_read_only"])
            receipt_binding = result["materialization_receipt"]
            receipt = _json(root / receipt_binding["path"])
            self.assertEqual(
                receipt_binding["sha256"],
                _sha256((root / receipt_binding["path"]).read_bytes()),
            )
            self.assertEqual(receipt["plan_id"], result["plan_id"])
            self.assertEqual(
                receipt["campaign_plan"],
                {
                    "path": "campaign_plan.json",
                    "sha256": result["campaign_plan_sha256"],
                },
            )
            self.assertEqual(
                receipt["tree_inventory"]["sha256"],
                result["materialized_member_inventory_sha256"],
            )
            self.assertEqual(
                receipt["tree_inventory"]["excludes"],
                [
                    execution.MATERIALIZATION_RECEIPT_NAME,
                    execution.immutable_orion_tree.FREEZE_RECEIPT_NAME,
                    execution.immutable_orion_tree.INVENTORY_NAME,
                ],
            )
            materialized_members = {
                path.relative_to(root).as_posix(): path.read_bytes()
                for path in root.rglob("*")
                if path.is_file()
                and path.name
                not in {
                    execution.MATERIALIZATION_RECEIPT_NAME,
                    execution.immutable_orion_tree.FREEZE_RECEIPT_NAME,
                    execution.immutable_orion_tree.INVENTORY_NAME,
                }
            }
            self.assertEqual(
                receipt["tree_inventory"]["sha256"],
                _sha256(execution._member_inventory_payload(materialized_members)),
            )
            self.assertEqual(plan["authorized_orion_root"], str(fixture["orion"]))
            self.assertEqual(
                plan["selected_pressure"]["selected_case"],
                {"case_id": "ps_p0_0p10", "problem_ps_p0": 0.1},
            )
            self.assertEqual(
                plan["candidate_binding"]["git_commit"],
                _json(fixture["candidate"])["source"]["git_commit"],
            )
            self.assertEqual(
                plan["candidate_binding"]["executable"]["sha256"],
                _sha256(fixture["executable"].read_bytes()),
            )
            self.assertEqual(
                plan["candidate_binding"]["environment_profile"]["sha256"],
                _sha256(fixture["environment"].read_bytes()),
            )
            self.assertEqual(
                plan["candidate_binding"]["environment_profile"]["reviewed_source"],
                {
                    "path": "tst/publication/frontier_control_plane/frontier_pic_environment.sh",
                    "sha256": _sha256(execution.REVIEWED_ENVIRONMENT_PROFILE_SOURCE.read_bytes()),
                },
            )

            descriptors = [
                _json(root / binding["path"])
                for binding in plan["baseline_attempt_descriptors"]
            ]
            self.assertEqual(
                [
                    (descriptor["variant"], descriptor["qualifying_seed"])
                    for descriptor in descriptors
                ],
                [
                    (variant, seed)
                    for variant in execution.CANONICAL_VARIANT_IDS
                    for seed in execution.QUALIFYING_SEEDS
                ],
            )
            for descriptor in descriptors:
                contract = _json(root / descriptor["launch_contract"]["path"])
                self.assertEqual(
                    descriptor["launch_contract"]["sha256"],
                    _sha256((root / descriptor["launch_contract"]["path"]).read_bytes()),
                )
                self.assertFalse(contract["launch_authorized"])
                self.assertFalse(contract["scheduler_submission_authorized"])
                self.assertIn("problem/ps_p0=0.1", contract["argv"])
            coarse_contract = _json(root / descriptors[0]["launch_contract"]["path"])
            amr_contract = _json(root / descriptors[8]["launch_contract"]["path"])
            self.assertTrue(
                set(execution.model.variant_binding("coarse_uniform_dx12").model_launch_overrides)
                <= set(coarse_contract["argv"])
            )
            self.assertFalse(
                set(amr_contract["argv"])
                & {
                    "mesh_refinement/refinement=adaptive",
                    "mesh_refinement/num_levels=3",
                    "problem/ps_enable_curvature_amr=true",
                }
            )

            carrier = _json(root / plan["restart_continuation_carrier"]["path"])
            self.assertEqual(carrier["variant"], "three_level_amr_root_dx12_finest_dx3")
            self.assertEqual(carrier["qualifying_seed"], execution.QUALIFYING_SEEDS[0])
            self.assertEqual(carrier["checkpoint_time_omega0_inverse"], 500.0)
            self.assertEqual(
                carrier["retained_output_schedule_after_checkpoint_omega0_inverse"],
                [600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0],
            )
            for path in [root, *root.rglob("*")]:
                self.assertFalse(path.stat().st_mode & _WRITE_BITS, path)

    def test_same_inputs_produce_byte_identical_plan_trees(self) -> None:
        with _fixture() as fixture:
            second_parent = fixture["root"] / "second-plans"
            second_parent.mkdir()
            first = _materialize(fixture)
            second = _materialize(fixture, output_parent=second_parent)
            first_root = Path(first["plan_root"])
            second_root = Path(second["plan_root"])
            try:
                self.assertEqual(first_root.name, second_root.name)
                first_files = {
                    path.relative_to(first_root).as_posix(): path.read_bytes()
                    for path in first_root.rglob("*")
                    if path.is_file()
                }
                second_files = {
                    path.relative_to(second_root).as_posix(): path.read_bytes()
                    for path in second_root.rglob("*")
                    if path.is_file()
                }
                self.assertEqual(first_files, second_files)
            finally:
                _make_writable(second_root)
                shutil.rmtree(second_root)

    def test_policy_fragment_and_independent_recompute_plan_remain_nonauthorizing(self) -> None:
        with _fixture() as fixture:
            result = _materialize(fixture)
            root = Path(result["plan_root"])
            plan = _json(root / "campaign_plan.json")
            fragment = _json(root / plan["nonauthorizing_policy_fragment"]["path"])
            recompute = _json(root / plan["independent_raw_artifact_recompute_plan"]["path"])
            self.assertFalse(fragment["mutates_live_policy"])
            self.assertFalse(fragment["scheduler_calls_authorized"])
            self.assertFalse(fragment["scheduler_submission_authorized"])
            self.assertFalse(fragment["frontier_execution_authorized"])
            self.assertFalse(fragment["launch_authorized"])
            self.assertFalse(recompute["production_helper_imports_authorized"])
            self.assertFalse(recompute["frontier_execution_authorized"])
            self.assertIn("reviewer-owned", recompute["implementation_rule"])
            closure = _json(root / plan["helper_source_closure"]["path"])
            self.assertIn(
                "tst/publication/q011_section54_restart.py",
                [record["path"] for record in closure["sources"]],
            )
            self.assertIn(
                "tst/publication/analyze_q011_section54_campaign.py",
                [record["path"] for record in closure["sources"]],
            )
            self.assertIn(
                "tst/publication/q011_section54_model.py",
                [record["path"] for record in closure["sources"]],
            )
            self.assertIn(
                "tst/publication/q011_section54_pressure_pilot_execution.py",
                [record["path"] for record in closure["sources"]],
            )
            for path in (
                "tst/publication/analyze_q011_section54_numerical_qualification.py",
                "tst/publication/q011_section54_particles.py",
                "tst/publication/q011_section54_spatial.py",
                "tst/publication/q011_section54_artifacts.py",
                "tst/publication/publish_q011_section54_pressure_pilot_bundle.py",
                "tst/publication/analyze_q011_section54_pressure_pilot.py",
                "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
                "tst/publication/frontier_f1_structured_artifacts.py",
                "tst/publication/q011_section54_attempt_manifest_materializer.py",
            ):
                with self.subTest(path=path):
                    self.assertIn(path, [record["path"] for record in closure["sources"]])

    def test_materializes_under_hidden_staging_before_atomic_no_replace_publish(self) -> None:
        with _fixture() as fixture:
            original_write = execution._write_new_file
            observed_staging: list[Path] = []

            def inspect_hidden_write(
                root: Path, descriptor: int, relative: str, payload: bytes
            ) -> None:
                observed_staging.append(root)
                self.assertTrue(
                    root.parent.name.startswith(".q011-section54-qualifying-campaign-plan-")
                )
                self.assertEqual(root.name, execution._STAGING_ROOT_NAME)
                self.assertFalse(
                    any(
                        path.name.startswith("q011-section54-qualifying-campaign-plan-")
                        for path in fixture["output_parent"].iterdir()
                    )
                )
                original_write(root, descriptor, relative, payload)

            with (
                patch.object(execution, "_write_new_file", side_effect=inspect_hidden_write),
                patch.object(
                    execution,
                    "_rename_no_replace_at",
                    wraps=execution._rename_no_replace_at,
                ) as rename,
            ):
                result = _materialize(fixture)
            self.assertTrue(observed_staging)
            rename.assert_called_once()
            self.assertTrue(Path(result["plan_root"]).is_dir())

    def test_destination_appearing_at_atomic_rename_is_not_overwritten(self) -> None:
        with _fixture() as fixture:
            original_rename = execution._rename_no_replace_at
            competing: list[Path] = []

            def collide(
                source_parent_descriptor: int,
                source_name: str,
                destination_parent_descriptor: int,
                destination_name: str,
            ) -> None:
                destination = fixture["output_parent"] / destination_name
                destination.mkdir()
                (destination / "sentinel.txt").write_text(
                    "competing publication\n", encoding="utf-8"
                )
                competing.append(destination)
                original_rename(
                    source_parent_descriptor,
                    source_name,
                    destination_parent_descriptor,
                    destination_name,
                )

            with patch.object(execution, "_rename_no_replace_at", side_effect=collide):
                with self.assertRaisesRegex(execution.CampaignPlanError, "already exists"):
                    _materialize(fixture)
            self.assertEqual(
                (competing[0] / "sentinel.txt").read_text(encoding="utf-8"),
                "competing publication\n",
            )
            self.assertEqual(list(fixture["output_parent"].iterdir()), competing)

    def test_publication_requires_atomic_no_replace_rename_support(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            (parent / "source").mkdir()
            descriptor = os.open(parent, execution._DIRECTORY_FLAGS)
            try:
                with patch.object(execution.ctypes, "CDLL", return_value=object()):
                    with self.assertRaisesRegex(
                        execution.CampaignPlanError,
                        "requires atomic no-replace rename support",
                    ):
                        execution._rename_no_replace_at(
                            descriptor, "source", descriptor, "destination"
                        )
                self.assertTrue((parent / "source").is_dir())
                self.assertFalse((parent / "destination").exists())
            finally:
                os.close(descriptor)

    def test_concurrent_injected_staging_member_fails_closed(self) -> None:
        with _fixture() as fixture:
            original_freeze = execution.immutable_orion_tree.freeze_tree_anchored

            def inject_member(root: Path, *args: object, **kwargs: object) -> dict[str, object]:
                (Path(root) / "injected-member.txt").write_text(
                    "not planner-owned\n", encoding="utf-8"
                )
                return original_freeze(root, *args, **kwargs)

            with patch.object(
                execution.immutable_orion_tree,
                "freeze_tree_anchored",
                side_effect=inject_member,
            ):
                with self.assertRaisesRegex(
                    execution.CampaignPlanError, "frozen inventory drifted"
                ):
                    _materialize(fixture)
            self.assertEqual(list(fixture["output_parent"].iterdir()), [])

    def test_member_injected_at_atomic_publish_is_rolled_back(self) -> None:
        with _fixture() as fixture:
            original_rename = execution._rename_no_replace_at
            injected = False

            def inject_then_rename(
                source_parent_descriptor: int,
                source_name: str,
                destination_parent_descriptor: int,
                destination_name: str,
            ) -> None:
                nonlocal injected
                if not injected:
                    staging = (
                        Path("/proc/self/fd") / str(source_parent_descriptor) / source_name
                    )
                    mode = stat.S_IMODE(staging.stat().st_mode)
                    staging.chmod(mode | stat.S_IWUSR | stat.S_IXUSR)
                    (staging / "injected-member.txt").write_text(
                        "not planner-owned\n", encoding="utf-8"
                    )
                    staging.chmod(mode)
                    injected = True
                original_rename(
                    source_parent_descriptor,
                    source_name,
                    destination_parent_descriptor,
                    destination_name,
                )

            with patch.object(
                execution, "_rename_no_replace_at", side_effect=inject_then_rename
            ) as rename:
                with self.assertRaises(execution.CampaignPlanError):
                    _materialize(fixture)
            self.assertEqual(rename.call_count, 2)
            self.assertEqual(list(fixture["output_parent"].iterdir()), [])

    def test_private_container_substitution_never_exposes_replacement(self) -> None:
        with _fixture() as fixture:
            original_rename = execution._rename_no_replace_at
            substituted = False
            parked: list[Path] = []
            observed_public_replacement: list[bool] = []
            safe_public_inodes: list[int] = []
            observed_public_inodes: list[int] = []
            observed_public_modes: list[int] = []

            def substitute_then_rename(
                source_parent_descriptor: int,
                source_name: str,
                destination_parent_descriptor: int,
                destination_name: str,
            ) -> None:
                nonlocal substituted
                if not substituted:
                    private = next(
                        path
                        for path in fixture["output_parent"].iterdir()
                        if path.name.startswith(
                            ".q011-section54-qualifying-campaign-plan-"
                        )
                    )
                    parked_root = private.with_name(f"{private.name}.parked")
                    private.rename(parked_root)
                    safe_public_inodes.append(
                        (parked_root / execution._STAGING_ROOT_NAME).stat().st_ino
                    )
                    private.mkdir()
                    (private / execution._STAGING_ROOT_NAME).mkdir()
                    (private / execution._STAGING_ROOT_NAME / "sentinel.txt").write_text(
                        "attacker replacement\n", encoding="utf-8"
                    )
                    parked.append(parked_root)
                    substituted = True
                original_rename(
                    source_parent_descriptor,
                    source_name,
                    destination_parent_descriptor,
                    destination_name,
                )
                destination = fixture["output_parent"] / destination_name
                published = destination.stat()
                observed_public_inodes.append(published.st_ino)
                observed_public_modes.append(stat.S_IMODE(published.st_mode))
                try:
                    replacement_is_public = (destination / "sentinel.txt").exists()
                except PermissionError:
                    replacement_is_public = False
                observed_public_replacement.append(replacement_is_public)

            with patch.object(
                execution, "_rename_no_replace_at", side_effect=substitute_then_rename
            ):
                result = _materialize(fixture)
            public = Path(result["plan_root"])
            self.assertTrue(public.is_dir())
            self.assertEqual(observed_public_replacement, [False])
            self.assertEqual(observed_public_inodes, safe_public_inodes)
            self.assertEqual(observed_public_modes, [stat.S_IWUSR])
            self.assertFalse((public / "sentinel.txt").exists())
            self.assertTrue(parked)

    def test_rollback_quarantines_moved_original_and_preserves_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            destination = parent / "campaign-plan"
            destination.mkdir()
            (destination / "owned.txt").write_text("owned\n", encoding="utf-8")
            moved = parent / "moved-public-campaign-plan"
            parent_descriptor = os.open(parent, execution._DIRECTORY_FLAGS)
            destination_descriptor = os.open(
                destination.name,
                execution._DIRECTORY_FLAGS,
                dir_fd=parent_descriptor,
            )
            original_rename = execution._rename_no_replace_at
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
                with patch.object(
                    execution, "_rename_no_replace_at", side_effect=race_then_rename
                ), self.assertRaisesRegex(
                    execution.CampaignPlanError, "substituted public destination"
                ):
                    execution._rollback_published_destination(
                        parent_descriptor,
                        destination.name,
                        destination_descriptor,
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
                os.close(destination_descriptor)
                os.close(parent_descriptor)

    def test_pressure_is_consumed_from_validated_human_receipt_and_never_inferred(self) -> None:
        with _fixture(receipt=_receipt(case_id="ps_p0_0p20", problem_ps_p0=0.2)) as fixture:
            result = _materialize(fixture)
            root = Path(result["plan_root"])
            plan = _json(root / "campaign_plan.json")
            self.assertEqual(plan["selected_pressure"]["selected_case"]["problem_ps_p0"], 0.2)
            descriptor = _json(root / plan["baseline_attempt_descriptors"][0]["path"])
            contract = _json(root / descriptor["launch_contract"]["path"])
            self.assertIn("problem/ps_p0=0.2", contract["argv"])

        automated = _receipt()
        automated["selection_method"] = "automated_minimum_score"
        with _fixture(receipt=automated) as fixture:
            with self.assertRaisesRegex(execution.CampaignPlanError, "human pressure-selection"):
                _materialize(fixture)

    def test_matrix_rejects_duplicate_boolean_and_reordered_aliases(self) -> None:
        matrix = {
            "physical_mode": execution.PHYSICAL_MODE,
            "grid_variants": list(execution.CANONICAL_VARIANT_IDS),
            "qualifying_seeds": list(execution.QUALIFYING_SEEDS),
            "expected_baseline_attempts": 24,
            "paired_seed_rule": (
                "Use the same qualifying seed for coarse-uniform, AMR and "
                "fine-uniform variants."
            ),
        }
        duplicate_variant = copy.deepcopy(matrix)
        duplicate_variant["grid_variants"][2] = duplicate_variant["grid_variants"][0]
        duplicate_seed = copy.deepcopy(matrix)
        duplicate_seed["qualifying_seeds"][7] = duplicate_seed["qualifying_seeds"][0]
        boolean_seed = copy.deepcopy(matrix)
        boolean_seed["qualifying_seeds"][0] = True
        boolean_count = copy.deepcopy(matrix)
        boolean_count["expected_baseline_attempts"] = True
        reordered = copy.deepcopy(matrix)
        reordered["grid_variants"].reverse()
        for invalid in (
            duplicate_variant,
            duplicate_seed,
            boolean_seed,
            boolean_count,
            reordered,
        ):
            with self.subTest(matrix=invalid), self.assertRaises(execution.CampaignPlanError):
                execution.validate_campaign_matrix(invalid)

    def test_hash_schema_symlink_escape_and_overwrite_fail_closed(self) -> None:
        with _fixture() as fixture:
            executable = fixture["executable"]
            executable.chmod(0o755)
            executable.write_bytes(b"drifted executable\n")
            executable.chmod(0o555)
            with self.assertRaisesRegex(execution.CampaignPlanError, "shared validation"):
                _materialize(fixture)

        with _fixture() as fixture:
            manifest = _json(fixture["candidate"])
            manifest["unexpected"] = False
            fixture["candidate"].chmod(0o644)
            fixture["candidate"].write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            fixture["candidate"].chmod(0o444)
            with self.assertRaisesRegex(execution.CampaignPlanError, "unexpected top-level fields"):
                _materialize(fixture)

        with _fixture() as fixture:
            alias = fixture["root"] / "plans-alias"
            alias.symlink_to(fixture["output_parent"], target_is_directory=True)
            with self.assertRaisesRegex(execution.CampaignPlanError, "symlink or path alias"):
                _materialize(fixture, output_parent=alias)

        with _fixture() as fixture:
            alias = fixture["root"] / "athena-alias"
            alias.symlink_to(fixture["executable"])
            with self.assertRaisesRegex(execution.CampaignPlanError, "supplied executable differs"):
                execution.materialize_qualifying_campaign_plan(
                    output_parent=fixture["output_parent"],
                    pressure_selection_receipt=fixture["pressure_receipt"],
                    clean_candidate_manifest=fixture["candidate"],
                    executable=alias,
                    environment_profile=fixture["environment"],
                )

        with _fixture() as fixture:
            _materialize(fixture)
            with self.assertRaisesRegex(execution.CampaignPlanError, "already exists"):
                _materialize(fixture)

    def test_rejects_non_elf_and_noncanonical_or_drifted_installed_environment(self) -> None:
        with _fixture(executable_payload=b"attested but not an ELF executable\n") as fixture:
            with self.assertRaisesRegex(execution.CampaignPlanError, "not an ELF binary"):
                _materialize(fixture)

        with _fixture() as fixture:
            noncanonical = _put(
                fixture["orion"] / "installed" / "frontier_pic_environment.sh",
                fixture["environment"].read_bytes(),
            )
            with self.assertRaisesRegex(execution.CampaignPlanError, "canonical control_plane"):
                execution.materialize_qualifying_campaign_plan(
                    output_parent=fixture["output_parent"],
                    pressure_selection_receipt=fixture["pressure_receipt"],
                    clean_candidate_manifest=fixture["candidate"],
                    executable=fixture["executable"],
                    environment_profile=noncanonical,
                )

        with _fixture() as fixture:
            fixture["environment"].chmod(0o644)
            fixture["environment"].write_bytes(b"#!/bin/bash\n# drifted reviewed environment\n")
            fixture["environment"].chmod(0o444)
            with self.assertRaisesRegex(execution.CampaignPlanError, "reviewed pressure-pilot source"):
                _materialize(fixture)

        with _fixture() as fixture:
            installed_helper = fixture["environment"].parent / "control_plane_common.py"
            installed_helper.chmod(0o644)
            installed_helper.write_bytes(installed_helper.read_bytes() + b"\n# drifted installed helper\n")
            installed_helper.chmod(0o444)
            with self.assertRaisesRegex(execution.CampaignPlanError, "Installed control-plane checksum"):
                _materialize(fixture)

    def test_rejects_frozen_build_provenance_drift_through_shared_validator(self) -> None:
        with _fixture() as fixture:
            modules = fixture["candidate"].parent / "build_provenance" / "modules.txt"
            modules.chmod(0o644)
            modules.write_bytes(modules.read_bytes() + b"drifted module\n")
            modules.chmod(0o444)
            with self.assertRaisesRegex(execution.CampaignPlanError, "provenance checksum mismatch"):
                _materialize(fixture)

    def test_descriptor_relative_output_write_rejects_nested_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "plan"
            outside = Path(directory) / "outside"
            root.mkdir()
            outside.mkdir()
            (root / "nested").symlink_to(outside, target_is_directory=True)
            descriptor = os.open(root, execution._DIRECTORY_FLAGS)
            try:
                with self.assertRaises(execution.CampaignPlanError):
                    execution._write_new_file(root, descriptor, "nested/member.json", b"{}\n")
            finally:
                os.close(descriptor)
            self.assertFalse((outside / "member.json").exists())

    def test_preregistration_and_deck_hash_drift_fail_closed(self) -> None:
        with _fixture() as fixture:
            drifted_preregistration = _put(
                fixture["root"] / "drifted-preregistration.json",
                execution.QUALIFYING_PREREGISTRATION.read_bytes() + b"\n",
            )
            with self.assertRaisesRegex(execution.CampaignPlanError, "preregistration SHA-256"):
                execution.materialize_qualifying_campaign_plan(
                    output_parent=fixture["output_parent"],
                    pressure_selection_receipt=fixture["pressure_receipt"],
                    clean_candidate_manifest=fixture["candidate"],
                    executable=fixture["executable"],
                    environment_profile=fixture["environment"],
                    qualifying_preregistration=drifted_preregistration,
                )
            drifted_deck = _put(
                fixture["root"] / "drifted.athinput",
                execution.PAPER_DECK.read_bytes() + b"\n",
            )
            with self.assertRaisesRegex(execution.CampaignPlanError, "paper deck SHA-256"):
                execution.materialize_qualifying_campaign_plan(
                    output_parent=fixture["output_parent"],
                    pressure_selection_receipt=fixture["pressure_receipt"],
                    clean_candidate_manifest=fixture["candidate"],
                    executable=fixture["executable"],
                    environment_profile=fixture["environment"],
                    paper_deck=drifted_deck,
                )


if __name__ == "__main__":
    unittest.main()
