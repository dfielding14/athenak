#!/usr/bin/env python3
"""Focused tests for the Q-011 registered-launch review materializer."""

from __future__ import annotations

import copy
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import stat
import tempfile
from typing import Any, Iterator
import unittest
from unittest.mock import patch

from tst.publication import q011_section54_registered_launch_materializer as materializer


_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH
_ATTEMPTS = (
    "baseline-001-coarse_uniform_dx12-seed-23050101",
    "baseline-002-coarse_uniform_dx12-seed-23050102",
)


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _put(path: Path, payload: bytes, mode: int = 0o444) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    path.chmod(mode)
    return path


def _binding(path: Path) -> dict[str, str]:
    return {"path": str(path), "sha256": _sha256(path.read_bytes())}


def _rewrite_json(path: Path, value: object) -> None:
    path.chmod(0o644)
    path.write_bytes(_json_bytes(value))
    path.chmod(0o444)


def _make_writable(root: Path) -> None:
    if not root.exists():
        return
    for path in [root, *root.rglob("*")]:
        if not path.is_symlink():
            path.chmod(path.stat().st_mode | stat.S_IWUSR)


@contextmanager
def _fixture(attempt_ids: tuple[str, ...] = _ATTEMPTS) -> Iterator[dict[str, Any]]:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory).resolve()
        pic_root = root / "pic"
        pic_root.mkdir()
        (pic_root / "plans").mkdir()
        (pic_root / "clean_candidates").mkdir()
        (pic_root / "control_plane").mkdir()
        (pic_root / "logs/slurm").mkdir(parents=True)
        plan_id = "a" * 64
        controller_version = "b" * 64
        git_commit = "c" * 40
        candidate_root = pic_root / "clean_candidates" / "fixture-freeze"
        candidate_root.mkdir()
        executable = _put(candidate_root / "athena", b"fixture-athena\n", 0o555)
        clean_manifest_value = {
            "schema_version": 4,
            "freeze_id": "fixture-freeze",
            "created_utc": "2026-06-06T12:00:00Z",
            "prepared_artifacts": {},
            "source": {"git_commit": git_commit},
            "build": {
                "executable_path": str(executable),
                "executable_sha256": _sha256(executable.read_bytes()),
            },
        }
        clean_manifest = _put(
            candidate_root / "clean_candidate_manifest.json",
            _json_bytes(clean_manifest_value),
        )
        controller = pic_root / "control_plane" / controller_version
        controller.mkdir()
        environment = _put(
            controller / "frontier_pic_environment.sh", b"fixture-environment\n"
        )
        controller.chmod(0o555)
        runtime = root / "runtime"
        runtime.mkdir()
        job_script = _put(
            runtime / "q011-job.sh",
            (
                "#!/bin/bash\n"
                "#SBATCH --account=AST207\n"
                "#SBATCH --partition=batch\n"
                "#SBATCH --qos=normal\n"
                "#SBATCH --nodes=2\n"
                "#SBATCH --time=01:00:00\n"
                f"#SBATCH --output={pic_root}/logs/slurm/%x.%j.log\n"
            ).encode("utf-8"),
        )
        input_deck = _put(runtime / "q011.athinput", b"<job>\nbasename=q011\n")
        analysis_one = _put(runtime / "analysis-one.py", b"print('one')\n")
        analysis_two = _put(runtime / "analysis-two.py", b"print('two')\n")
        candidate_binding = {
            "clean_candidate_manifest": _binding(clean_manifest),
            "executable": _binding(executable),
            "environment_profile": {
                **_binding(environment),
                "control_plane_version": controller_version,
                "reviewed_source": {
                    "path": "tst/publication/frontier_control_plane/frontier_pic_environment.sh",
                    "sha256": _sha256(environment.read_bytes()),
                },
            },
            "git_commit": git_commit,
        }
        planner_root = (
            pic_root / "plans" / f"q011-section54-qualifying-campaign-plan-{plan_id}"
        )
        planner_root.mkdir()
        plan = {
            "record_type": "q011_section54_qualifying_campaign_execution_plan",
            "schema_version": 1,
            "plan_id": plan_id,
            "status": "source_local_immutable_review_plan_only",
            "candidate_binding": candidate_binding,
            "source_bindings": {
                "paper_deck": {
                    "path": "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
                    "sha256": _sha256(input_deck.read_bytes()),
                }
            },
            "execution_boundary": {
                "mutates_live_policy": False,
                "scheduler_calls": False,
                "submits_jobs": False,
                "infers_pressure_selection": False,
                "launch_authorized": False,
                "frontier_execution_authorized": False,
                "claim_closure_authorized": False,
            },
        }
        plan_path = _put(planner_root / "campaign_plan.json", _json_bytes(plan))
        inventory_payload = (
            f"{_sha256(plan_path.read_bytes())}  campaign_plan.json\n".encode("utf-8")
        )
        inventory = _put(planner_root / "artifact_inventory.sha256", inventory_payload)
        planner_root.chmod(0o555)
        retention_root = root / "retention"
        retention_root.mkdir()
        retention_bindings = []
        for attempt_id in attempt_ids:
            retention = {
                "schema_version": 1,
                "retention_role": "q011_section54_deterministic_retained_attempt",
                "planner_root": str(planner_root),
                "planner_inventory_sha256": _sha256(inventory.read_bytes()),
                "planner_plan_id": plan_id,
                "planner_materialization_receipt": {
                    "path": "materialization_receipt.json",
                    "sha256": "d" * 64,
                },
                "attempt_id": attempt_id,
                "authorized_orion_attempt_root": str(
                    pic_root / "campaigns" / f"q011-section54-{plan_id}" / "baseline" / attempt_id
                ),
                "authorized_orion_raw_root": str(
                    pic_root
                    / "campaigns"
                    / f"q011-section54-{plan_id}"
                    / "baseline"
                    / attempt_id
                    / "raw"
                ),
                "argv": [
                    "-i",
                    "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
                    "-d",
                    str(
                        pic_root
                        / "campaigns"
                        / f"q011-section54-{plan_id}"
                        / "baseline"
                        / attempt_id
                        / "raw"
                    ),
                    f"job/basename={attempt_id}",
                    "problem/ps_p0=1.0",
                    "particles/pic_random_seed=23050101",
                ],
            }
            path = _put(retention_root / f"{attempt_id}.json", _json_bytes(retention))
            retention_bindings.append(
                {"attempt_id": attempt_id, "planner_retention": _binding(path)}
            )
        selected_value = {
            "record_type": materializer.SELECTED_BINDING_RECORD_TYPE,
            "schema_version": 1,
            "authorized_pic_root": str(pic_root),
            "planner": {
                "root": str(planner_root),
                "inventory_sha256": _sha256(inventory.read_bytes()),
                "plan_id": plan_id,
            },
            "clean_candidate": {
                "manifest": _binding(clean_manifest),
                "executable": _binding(executable),
                "git_commit": git_commit,
            },
            "controller": {
                "path": str(controller),
                "version": controller_version,
                "environment_profile": _binding(environment),
            },
            "runtime_artifacts": {
                "job_script": _binding(job_script),
                "input_deck": _binding(input_deck),
                "analysis_scripts": [_binding(analysis_one), _binding(analysis_two)],
            },
            "resource_profile": {
                "selected_qos": "normal",
                "qos_selection_reason": "normal_required_by_registered_campaign",
                "registered_short_nonproduction": False,
                "maximum_nodes": 2,
                "maximum_walltime_seconds": 3600,
                "maximum_attempts_per_slice": 1,
                "maximum_batch_node_hours": float(2 * len(attempt_ids)),
                "launch_resources": {
                    "nodes": 2,
                    "tasks": 16,
                    "cpus_per_task": 7,
                    "gpus_per_task": 1,
                    "gpu_bind": "closest",
                },
            },
            "attempts": retention_bindings,
        }
        selected = _put(root / "selected-bindings.json", _json_bytes(selected_value))
        output = root / "review-bundle"
        try:
            with (
                patch.object(
                    materializer,
                    "verify_installed_control_plane",
                    return_value={"version": controller_version},
                ) as verify_controller,
                patch.object(
                    materializer,
                    "validate_planner_retention_binding",
                    side_effect=lambda value, **_: copy.deepcopy(value),
                ) as validate_retention,
            ):
                yield {
                    "root": root,
                    "pic_root": pic_root,
                    "plan_id": plan_id,
                    "planner_root": planner_root,
                    "controller": controller,
                    "controller_version": controller_version,
                    "environment": environment,
                    "clean_manifest": clean_manifest,
                    "executable": executable,
                    "job_script": job_script,
                    "input_deck": input_deck,
                    "analysis_one": analysis_one,
                    "selected": selected,
                    "selected_value": selected_value,
                    "output": output,
                    "verify_controller": verify_controller,
                    "validate_retention": validate_retention,
                }
        finally:
            _make_writable(root)


def _materialize(fixture: dict[str, Any]) -> dict[str, Any]:
    return materializer.materialize_registered_launch_review_bundle(
        selected_bindings=fixture["selected"],
        output_root=fixture["output"],
        authorized_pic_root=fixture["pic_root"],
    )


class RegisteredLaunchMaterializerTests(unittest.TestCase):
    def test_materializes_complete_slice_candidates_and_incomplete_configs(self) -> None:
        with _fixture() as fixture:
            result = _materialize(fixture)
            self.assertFalse(result["launch_authorized"])
            self.assertTrue(result["recursively_read_only"])
            self.assertEqual(result["included_attempt_count"], 2)
            self.assertEqual(result["maximum_batch_node_hours"], 4.0)
            fragment = json.loads(
                (fixture["output"] / "registered_policy_slice_fragment.json").read_text()
            )
            self.assertEqual(len(fragment["registered_science_slices"]), 2)
            self.assertEqual(fragment["declared_batch_ceiling"]["maximum_registered_attempts"], 2)
            self.assertFalse(fragment["execution_boundary"]["complete_storage_policy"])
            for policy_slice in fragment["registered_science_slices"]:
                self.assertEqual(set(policy_slice), materializer.REGISTERED_SLICE_KEYS)
                self.assertEqual(policy_slice["status"], "authorized")
                self.assertEqual(policy_slice["maximum_attempts"], 1)
                self.assertEqual(policy_slice["maximum_nodes"], 2)
                self.assertEqual(policy_slice["maximum_walltime_seconds"], 3600)
                self.assertRegex(policy_slice["launch_contract_sha256"], r"^[0-9a-f]{64}$")
            for attempt_id in _ATTEMPTS:
                candidate = json.loads(
                    (
                        fixture["output"]
                        / "pre_submit_config_candidates"
                        / f"{attempt_id}.json"
                    ).read_text()
                )
                self.assertEqual(
                    candidate["status"],
                    "review_candidate_incomplete_not_pre_submit_config",
                )
                self.assertFalse(
                    candidate["execution_boundary"]["complete_pre_submit_config"]
                )
                static = candidate["static_pre_submit_config"]
                self.assertEqual(static["planner_retention"]["attempt_id"], attempt_id)
                self.assertNotIn("submission_id", static)
                self.assertNotIn("pre_manifest_attestation", static)
                self.assertNotIn("timeout_margin_artifact", static)
                self.assertNotIn("queue_snapshot", static)
                materializer.validate_launch_contract(static["launch_contract"])
            for path in [fixture["output"], *fixture["output"].rglob("*")]:
                self.assertFalse(path.stat().st_mode & _WRITE_BITS, path)

    def test_revalidates_every_included_retention_against_clean_candidate(self) -> None:
        with _fixture() as fixture:
            _materialize(fixture)
            validator = fixture["validate_retention"]
            self.assertEqual(validator.call_count, 2)
            for call in validator.call_args_list:
                self.assertEqual(call.kwargs["authorized_pic_root"], fixture["pic_root"])
                self.assertEqual(
                    call.kwargs["expected_clean_candidate_manifest_sha256"],
                    _sha256(fixture["clean_manifest"].read_bytes()),
                )

    def test_integrates_with_real_control_plane_retention_validator(self) -> None:
        from tst.publication.frontier_control_plane.test_control_plane import (
            SnapshotTests,
        )

        control_fixture = SnapshotTests(methodName="runTest")
        control_fixture.setUp()
        try:
            retention = control_fixture._planner_retention()
            pic_root = control_fixture.pic_root
            planner_root = Path(retention["planner_root"])
            plan = json.loads((planner_root / "campaign_plan.json").read_text())
            clean_manifest = Path(control_fixture._planner_candidate_manifest)
            manifest = json.loads(clean_manifest.read_text())
            executable = clean_manifest.parent / "athena"
            controller = control_fixture.control_plane_dir
            environment = controller / "frontier_pic_environment.sh"
            source_root = control_fixture.root / "registered-launch-integration"
            source_root.mkdir()
            paper_binding = plan["source_bindings"]["paper_deck"]
            input_deck = _put(
                source_root / "q011.athinput",
                (planner_root / paper_binding["path"]).read_bytes(),
            )
            job_script = _put(
                source_root / "q011-job.sh",
                (
                    "#!/bin/bash\n"
                    "#SBATCH --account=AST207\n"
                    "#SBATCH --partition=batch\n"
                    "#SBATCH --qos=normal\n"
                    "#SBATCH --nodes=1\n"
                    "#SBATCH --time=01:00:00\n"
                    f"#SBATCH --output={pic_root}/logs/slurm/%x.%j.log\n"
                ).encode("utf-8"),
            )
            analysis = _put(source_root / "analysis.py", b"print('integration')\n")
            retention_path = _put(
                source_root / "planner-retention.json", _json_bytes(retention)
            )
            selected_value = {
                "record_type": materializer.SELECTED_BINDING_RECORD_TYPE,
                "schema_version": 1,
                "authorized_pic_root": str(pic_root),
                "planner": {
                    "root": str(planner_root),
                    "inventory_sha256": retention["planner_inventory_sha256"],
                    "plan_id": retention["planner_plan_id"],
                },
                "clean_candidate": {
                    "manifest": _binding(clean_manifest),
                    "executable": _binding(executable),
                    "git_commit": manifest["source"]["git_commit"],
                },
                "controller": {
                    "path": str(controller),
                    "version": controller.name,
                    "environment_profile": _binding(environment),
                },
                "runtime_artifacts": {
                    "job_script": _binding(job_script),
                    "input_deck": _binding(input_deck),
                    "analysis_scripts": [_binding(analysis)],
                },
                "resource_profile": {
                    "selected_qos": "normal",
                    "qos_selection_reason": "normal_required_by_registered_campaign",
                    "registered_short_nonproduction": False,
                    "maximum_nodes": 1,
                    "maximum_walltime_seconds": 3600,
                    "maximum_attempts_per_slice": 1,
                    "maximum_batch_node_hours": 1.0,
                    "launch_resources": {
                        "nodes": 1,
                        "tasks": 8,
                        "cpus_per_task": 7,
                        "gpus_per_task": 1,
                        "gpu_bind": "closest",
                    },
                },
                "attempts": [
                    {
                        "attempt_id": retention["attempt_id"],
                        "planner_retention": _binding(retention_path),
                    }
                ],
            }
            selected = _put(
                source_root / "selected-bindings.json", _json_bytes(selected_value)
            )
            result = materializer.materialize_registered_launch_review_bundle(
                selected_bindings=selected,
                output_root=source_root / "review-bundle",
                authorized_pic_root=pic_root,
            )
            self.assertEqual(result["included_attempt_count"], 1)
            self.assertFalse(result["launch_authorized"])
        finally:
            _make_writable(control_fixture.root)
            control_fixture.tearDown()
            control_fixture.doCleanups()

    def test_rejects_unknown_selected_binding_data(self) -> None:
        with _fixture() as fixture:
            value = copy.deepcopy(fixture["selected_value"])
            value["unexpected"] = "not allowed"
            _rewrite_json(fixture["selected"], value)
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError, "schema is malformed"
            ):
                _materialize(fixture)

    def test_rejects_missing_retention_record(self) -> None:
        with _fixture() as fixture:
            retention = Path(
                fixture["selected_value"]["attempts"][0]["planner_retention"]["path"]
            )
            retention.unlink()
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError, "cannot be opened"
            ):
                _materialize(fixture)

    def test_rejects_retention_attempt_drift(self) -> None:
        with _fixture() as fixture:
            retention = Path(
                fixture["selected_value"]["attempts"][0]["planner_retention"]["path"]
            )
            value = json.loads(retention.read_text())
            value["attempt_id"] = _ATTEMPTS[1]
            _rewrite_json(retention, value)
            selected = copy.deepcopy(fixture["selected_value"])
            selected["attempts"][0]["planner_retention"] = _binding(retention)
            _rewrite_json(fixture["selected"], selected)
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "attempt binding drifted",
            ):
                _materialize(fixture)

    def test_rejects_duplicate_or_unordered_attempt_scope(self) -> None:
        with _fixture() as fixture:
            value = copy.deepcopy(fixture["selected_value"])
            value["attempts"] = list(reversed(value["attempts"]))
            _rewrite_json(fixture["selected"], value)
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "uniquely ordered",
            ):
                _materialize(fixture)

    def test_rejects_validator_drift_or_failure(self) -> None:
        with self.subTest("failure"), _fixture() as fixture:
            fixture["validate_retention"].side_effect = ValueError("fixture rejection")
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "planner retention is invalid",
            ):
                _materialize(fixture)
        with self.subTest("drift"), _fixture() as fixture:
            fixture["validate_retention"].side_effect = lambda value, **_: {
                **value,
                "argv": [*value["argv"], "problem/drift=true"],
            }
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "validator returned drifted data",
            ):
                _materialize(fixture)

    def test_rejects_selected_candidate_or_controller_drift(self) -> None:
        with self.subTest("candidate"), _fixture() as fixture:
            value = copy.deepcopy(fixture["selected_value"])
            value["clean_candidate"]["git_commit"] = "f" * 40
            _rewrite_json(fixture["selected"], value)
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "Git commit drifted",
            ):
                _materialize(fixture)
        with self.subTest("controller"), _fixture() as fixture:
            fixture["verify_controller"].return_value = {"version": "f" * 64}
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "version drifted",
            ):
                _materialize(fixture)

    def test_rejects_input_deck_and_analysis_drift(self) -> None:
        with self.subTest("deck"), _fixture() as fixture:
            value = copy.deepcopy(fixture["selected_value"])
            fixture["input_deck"].chmod(0o644)
            fixture["input_deck"].write_bytes(b"drifted deck\n")
            fixture["input_deck"].chmod(0o444)
            value["runtime_artifacts"]["input_deck"] = _binding(fixture["input_deck"])
            _rewrite_json(fixture["selected"], value)
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "differs from immutable planner paper deck",
            ):
                _materialize(fixture)
        with self.subTest("analysis"), _fixture() as fixture:
            fixture["analysis_one"].chmod(0o644)
            fixture["analysis_one"].write_bytes(b"drifted analysis\n")
            fixture["analysis_one"].chmod(0o444)
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError, "SHA-256 drifted"
            ):
                _materialize(fixture)

    def test_rejects_job_script_and_resource_ceiling_drift(self) -> None:
        with self.subTest("job nodes"), _fixture() as fixture:
            value = copy.deepcopy(fixture["selected_value"])
            value["resource_profile"]["maximum_nodes"] = 3
            value["resource_profile"]["maximum_batch_node_hours"] = 6.0
            _rewrite_json(fixture["selected"], value)
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "nodes differ from resource ceiling",
            ):
                _materialize(fixture)
        with self.subTest("batch ceiling"), _fixture() as fixture:
            value = copy.deepcopy(fixture["selected_value"])
            value["resource_profile"]["maximum_batch_node_hours"] = 3.0
            _rewrite_json(fixture["selected"], value)
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "differs from attempt/resource ceilings",
            ):
                _materialize(fixture)

    def test_rejects_noncanonical_or_writable_selected_bindings(self) -> None:
        with self.subTest("noncanonical"), _fixture() as fixture:
            fixture["selected"].chmod(0o644)
            fixture["selected"].write_text(
                json.dumps(fixture["selected_value"]), encoding="utf-8"
            )
            fixture["selected"].chmod(0o444)
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "canonical JSON",
            ):
                _materialize(fixture)
        with self.subTest("writable"), _fixture() as fixture:
            fixture["selected"].chmod(0o644)
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError, "must be read-only"
            ):
                _materialize(fixture)

    def test_rejects_protected_live_namespace_and_existing_output(self) -> None:
        with self.subTest("protected"), _fixture() as fixture:
            fixture["output"] = fixture["pic_root"] / "policy" / "candidate"
            fixture["output"].parent.mkdir()
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "live PIC control namespace",
            ):
                _materialize(fixture)
        with self.subTest("existing"), _fixture() as fixture:
            fixture["output"].mkdir()
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "already exists",
            ):
                _materialize(fixture)

    def test_readiness_contract_drift_fails_closed(self) -> None:
        with _fixture() as fixture:
            readiness = json.loads(materializer.READINESS_CONTRACT.read_text())
            readiness["status"] = "drifted"
            path = _put(fixture["root"] / "readiness.json", _json_bytes(readiness))
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError,
                "readiness contract drifted",
            ):
                materializer.materialize_registered_launch_review_bundle(
                    selected_bindings=fixture["selected"],
                    output_root=fixture["output"],
                    authorized_pic_root=fixture["pic_root"],
                    readiness_contract=path,
                )

    def test_rejects_writable_external_readiness_contract(self) -> None:
        with _fixture() as fixture:
            readiness = json.loads(materializer.READINESS_CONTRACT.read_text())
            path = _put(
                fixture["root"] / "readiness.json", _json_bytes(readiness), mode=0o644
            )
            with self.assertRaisesRegex(
                materializer.RegisteredLaunchMaterializationError, "must be read-only"
            ):
                materializer.materialize_registered_launch_review_bundle(
                    selected_bindings=fixture["selected"],
                    output_root=fixture["output"],
                    authorized_pic_root=fixture["pic_root"],
                    readiness_contract=path,
                )


if __name__ == "__main__":
    unittest.main()
