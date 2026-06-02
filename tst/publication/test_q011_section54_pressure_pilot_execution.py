#!/usr/bin/env python3
"""Fail-closed tests for the Q-011 pressure-pilot registered execution tranche."""

from __future__ import annotations

from contextlib import contextmanager, redirect_stderr
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
from typing import Iterator
import unittest
import uuid
from unittest.mock import patch

from tst.publication import q011_section54_pressure_pilot_execution as execution


CONTROL_PLANE = Path(__file__).resolve().parent / "frontier_control_plane"
sys.path.insert(0, str(CONTROL_PLANE))
from control_plane_common import launch_contract_sha256, validate_launch_contract  # noqa: E402


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _put(path: Path, payload: bytes, mode: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    path.chmod(mode)
    return path


class PressurePilotExecutionTest(unittest.TestCase):
    @contextmanager
    def _final_binding(
        self, root: Path, *, case: execution.PressureCase = execution.CASES[0]
    ) -> Iterator[dict[str, object]]:
        root.mkdir(parents=True, exist_ok=True)
        freeze = root / str(uuid.uuid4())
        freeze.mkdir()
        executable = _put(freeze / "athena", b"exact-clean-athena\n", 0o555)
        executable_sha256 = _sha256(executable.read_bytes())
        manifest = {
            "schema_version": 4,
            "freeze_id": freeze.name,
            "source": {"git_commit": "a" * 40},
            "build": {
                "executable_path": str(executable),
                "executable_sha256": executable_sha256,
            },
        }
        clean_manifest = _put(
            freeze / "clean_candidate_manifest.json",
            (json.dumps(manifest) + "\n").encode("utf-8"),
            0o444,
        )
        environment = _put(
            root / "installed" / "frontier_pic_environment.sh",
            execution.ENVIRONMENT_PROFILE_SOURCE.read_bytes(),
            0o444,
        )
        environment_sha256 = _sha256(environment.read_bytes())
        timeout = {
            "athena_walltime_seconds": 840,
            "scheduler_walltime_seconds": 900,
            "environment_profile_sha256": environment_sha256,
            "measured_utc": "2026-06-02T00:00:00Z",
            "expires_utc": "2026-06-03T00:00:00Z",
        }
        timeout_margin = _put(
            root / "timeout_margin.json",
            (json.dumps(timeout) + "\n").encode("utf-8"),
            0o644,
        )
        queue_snapshot = _put(root / "queue_snapshot.txt", b"", 0o644)
        pre_manifest_attestation = _put(
            root / "operator_attestations" / "pre_manifest" / "attestation.json",
            (
                json.dumps(
                    {
                        "phase": "pre_manifest",
                        "registered_science_authorization_id": case.authorization_id,
                    }
                )
                + "\n"
            ).encode("utf-8"),
            0o444,
        )
        yield {
            "case_id": case.case_id,
            "submission_id": str(uuid.uuid4()),
            "clean_candidate_manifest": clean_manifest,
            "executable": executable,
            "environment_profile": environment,
            "pre_manifest_attestation": pre_manifest_attestation,
            "timeout_margin_artifact": timeout_margin,
            "queue_snapshot": queue_snapshot,
            "site_policy_checked_utc": "2026-06-02T00:00:00Z",
        }

    def test_source_preregistration_and_shared_directive_only_template(self) -> None:
        preregistration = execution.validate_source_tranche()
        self.assertEqual(
            preregistration["record_type"],
            "q011_section54_pressure_pilot_registered_execution_preregistration",
        )
        boundary = preregistration["execution_boundary"]
        self.assertFalse(boundary["frontier_execution_authorized_by_this_record"])
        self.assertFalse(boundary["storage_policy_mutation_authorized_by_this_record"])
        self.assertFalse(boundary["prepared_inventory_mutation_authorized_by_this_record"])
        source_bindings = preregistration["source_bindings"]
        self.assertEqual(
            source_bindings["analysis_scripts"],
            [
                {
                    "path": "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
                    "sha256": "d9253fdd3b573bd87f4e94a84e88ac2d8110a0cbc608ce98fecd7b9eeece9693",
                },
                {
                    "path": "tst/publication/frontier_f1_structured_artifacts.py",
                    "sha256": "cf090115bcdfd143f67b12339115102b57e74cf3205b1ebb1521c144c3415a5a",
                },
            ],
        )
        self.assertEqual(
            source_bindings["materializer"],
            {
                "path": "tst/publication/q011_section54_pressure_pilot_execution.py",
                "sha256": "510ac6760231f39690f09afb6ec7855f0b03ba487938ea09f107fa1f68724675",
            },
        )
        self.assertNotIn(
            "q011_section54_pressure_pilot_registered_execution_preregistration_2026-06-02.json",
            json.dumps(source_bindings, sort_keys=True),
        )
        self.assertNotIn(
            "tst/publication/q011_section54_pressure_pilot_execution.py",
            json.dumps(source_bindings["analysis_scripts"], sort_keys=True),
        )
        template = execution.JOB_SCRIPT.read_text(encoding="utf-8").splitlines()
        self.assertTrue(template[0].startswith("#!"))
        self.assertTrue(all(not line or line.startswith("#") for line in template))
        directives = {
            line.removeprefix("#SBATCH ").split("=", 1)[0]: line.split("=", 1)[1]
            for line in template
            if line.startswith("#SBATCH ")
        }
        self.assertEqual(directives["--account"], "AST207")
        self.assertEqual(directives["--partition"], "batch")
        self.assertEqual(directives["--qos"], "debug")
        self.assertEqual(directives["--nodes"], "1")
        self.assertEqual(directives["--time"], "00:15:00")

    def test_each_case_is_one_trusted_action_with_exact_overrides_and_stdout_sha(
        self,
    ) -> None:
        for case in execution.CASES:
            with self.subTest(case_id=case.case_id):
                contract = json.loads(
                    case.launch_contract_path.read_text(encoding="utf-8")
                )
                validate_launch_contract(contract)
                self.assertEqual(contract, execution.expected_launch_contract(case))
                self.assertEqual(len(contract["actions"]), 1)
                action = contract["actions"][0]
                self.assertEqual(
                    action["resources"],
                    {
                        "nodes": 1,
                        "tasks": 1,
                        "cpus_per_task": 7,
                        "gpus_per_task": 1,
                        "gpu_bind": "closest",
                    },
                )
                literals = [
                    token["literal"]
                    for token in action["arguments"]
                    if "literal" in token
                ]
                self.assertEqual(literals[2], f"job/basename={case.case_id}")
                self.assertEqual(literals[3:-1], list(execution.FIXED_OVERRIDES))
                self.assertEqual(literals[-1], f"problem/ps_p0={case.argv_value}")
                self.assertEqual(
                    contract["post_actions"][-1],
                    {
                        "action_id": "sha-pressure-stdout",
                        "kind": "artifact_sha256",
                        "artifact": "athena_stdout.txt",
                        "output_artifact": "athena_stdout.sha256",
                    },
                )
                measured = next(
                    binding
                    for binding in execution.source_bindings()["launch_contracts"]
                    if binding["case_id"] == case.case_id
                )
                self.assertEqual(
                    measured["launch_contract_sha256"],
                    launch_contract_sha256(contract),
                )

    def test_materialized_policy_uses_four_separate_bounded_slices(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self._final_binding(Path(directory)) as binding:
                slices = execution.materialize_registered_science_slices(
                    clean_candidate_manifest=binding["clean_candidate_manifest"],
                    executable=binding["executable"],
                    environment_profile=binding["environment_profile"],
                )
        self.assertEqual(len(slices), 4)
        self.assertEqual(len({record["authorization_id"] for record in slices}), 4)
        self.assertEqual(len({record["campaign"] for record in slices}), 4)
        for record in slices:
            self.assertEqual(record["status"], "authorized")
            self.assertEqual(record["evidence_class"], "engineering_calibration_only")
            self.assertEqual(record["physical_mode"], "paper_mhd_pic_vl2_tsc")
            self.assertEqual(record["runtime_profile"], "frontier_minimum_supported")
            self.assertEqual(record["selected_qos"], "debug")
            self.assertTrue(record["registered_short_nonproduction"])
            self.assertEqual(record["maximum_nodes"], 1)
            self.assertEqual(record["maximum_walltime_seconds"], 900)
            self.assertEqual(record["maximum_attempts"], 1)
            for field in [
                "job_script_sha256",
                "input_deck_sha256",
                "environment_profile_sha256",
                "executable_sha256",
                "launch_contract_sha256",
                "clean_candidate_manifest_sha256",
            ]:
                self.assertRegex(record[field], r"^[0-9a-f]{64}$")
            self.assertEqual(
                record["analysis_script_sha256"],
                [
                    "d9253fdd3b573bd87f4e94a84e88ac2d8110a0cbc608ce98fecd7b9eeece9693",
                    "cf090115bcdfd143f67b12339115102b57e74cf3205b1ebb1521c144c3415a5a",
                ],
            )

    def test_reviewed_config_emits_exactly_one_selected_case_and_writes_exclusively(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for case in execution.CASES:
                with self.subTest(case_id=case.case_id):
                    with self._final_binding(root / case.case_id, case=case) as binding:
                        config = execution.materialize_reviewed_pre_submit_config(
                            **binding
                        )
                    self.assertEqual(config["campaign"], case.campaign)
                    self.assertEqual(config["test_id"], case.test_id)
                    self.assertEqual(
                        config["registered_science_authorization_id"],
                        case.authorization_id,
                    )
                    self.assertEqual(config["submission_scope"], "registered_science")
                    self.assertEqual(config["physical_mode"], execution.PHYSICAL_MODE)
                    self.assertEqual(
                        config["artifact_dir"],
                        str(
                            execution.AUTHORIZED_PIC_ROOT
                            / "runs"
                            / case.campaign
                            / binding["submission_id"]
                        ),
                    )
                    self.assertEqual(
                        config["launch_contract"],
                        execution.expected_launch_contract(case),
                    )
                    self.assertEqual(
                        config["analysis_scripts"],
                        [
                            str(
                                execution.REPO_ROOT
                                / "tst/publication/analyze_q011_section54_pressure_pilot_case.py"
                            ),
                            str(
                                execution.REPO_ROOT
                                / "tst/publication/frontier_f1_structured_artifacts.py"
                            ),
                        ],
                    )
                    output = root / f"{case.case_id}-reviewed-config"
                    manifest = execution.write_reviewed_pre_submit_config(
                        output, **binding
                    )
                    self.assertEqual(manifest["case_id"], case.case_id)
                    self.assertEqual(
                        set(path.name for path in output.iterdir()),
                        {"pre_submit_config.json", "materialization_manifest.json"},
                    )
                    self.assertEqual(output.stat().st_mode & 0o222, 0)
                    config_path = output / manifest["config"]["path"]
                    self.assertEqual(config_path.stat().st_mode & 0o222, 0)
                    self.assertEqual(
                        _sha256(config_path.read_bytes()),
                        manifest["config"]["sha256"],
                    )
                    self.assertEqual(
                        manifest["captured_inputs"]["six_field_queue_snapshot"]["format"],
                        "%i|%P|%q|%T|%j|%k",
                    )
                    self.assertEqual(len(manifest["required_next_steps"]), 3)
                    with self.assertRaisesRegex(
                        execution.ContractError, "already exists"
                    ):
                        execution.write_reviewed_pre_submit_config(output, **binding)

    def test_policy_fragment_write_is_additive_and_exclusive(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self._final_binding(root) as binding:
                output = root / "policy-fragment.json"
                fragment = execution.write_policy_fragment(
                    output,
                    clean_candidate_manifest=binding["clean_candidate_manifest"],
                    executable=binding["executable"],
                    environment_profile=binding["environment_profile"],
                )
                self.assertEqual(len(fragment["registered_science_slices"]), 4)
                self.assertEqual(output.stat().st_mode & 0o222, 0)
                with self.assertRaisesRegex(
                    execution.ContractError, "refusing to overwrite"
                ):
                    execution.write_policy_fragment(
                        output,
                        clean_candidate_manifest=binding["clean_candidate_manifest"],
                        executable=binding["executable"],
                        environment_profile=binding["environment_profile"],
                    )

    def test_final_bindings_fail_closed_on_environment_and_executable_drift(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self._final_binding(root) as binding:
                changed_environment = _put(
                    root / "changed" / "frontier_pic_environment.sh",
                    b"changed\n",
                    0o444,
                )
                with self.assertRaisesRegex(
                    execution.ContractError, "environment profile differs"
                ):
                    execution.materialize_registered_science_slices(
                        clean_candidate_manifest=binding["clean_candidate_manifest"],
                        executable=binding["executable"],
                        environment_profile=changed_environment,
                    )
                detached = _put(root / "detached-athena", b"exact-clean-athena\n", 0o555)
                with self.assertRaisesRegex(execution.ContractError, "not adjacent"):
                    execution.materialize_registered_science_slices(
                        clean_candidate_manifest=binding["clean_candidate_manifest"],
                        executable=detached,
                        environment_profile=binding["environment_profile"],
                    )

    def test_timeout_submission_id_and_source_drift_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self._final_binding(root) as binding:
                timeout = json.loads(
                    Path(binding["timeout_margin_artifact"]).read_text(encoding="utf-8")
                )
                timeout["scheduler_walltime_seconds"] = True
                Path(binding["timeout_margin_artifact"]).write_text(
                    json.dumps(timeout), encoding="utf-8"
                )
                with self.assertRaisesRegex(execution.ContractError, "exactly 900"):
                    execution.materialize_reviewed_pre_submit_config(**binding)
            with self._final_binding(root / "second") as binding:
                binding["case_id"] = "ps_p0_unknown"
                with self.assertRaisesRegex(execution.ContractError, "unknown"):
                    execution.materialize_reviewed_pre_submit_config(**binding)
            with self._final_binding(root / "third") as binding:
                Path(binding["queue_snapshot"]).write_text(
                    "1|AST207|batch|debug|RUNNING|job|comment\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(execution.ContractError, "six-field"):
                    execution.materialize_reviewed_pre_submit_config(**binding)
            with self.assertRaisesRegex(execution.ContractError, "required"):
                execution._selected_case(None)
        with patch.object(execution, "source_bindings", return_value={}):
            with self.assertRaisesRegex(execution.ContractError, "source bindings drifted"):
                execution.validate_source_tranche()

    def test_pre_submit_cli_requires_one_known_case_id(self) -> None:
        required = [
            "pre-submit-config",
            "--submission-id",
            str(uuid.uuid4()),
            "--clean-candidate-manifest",
            "/tmp/clean_candidate_manifest.json",
            "--executable",
            "/tmp/athena",
            "--environment-profile",
            "/tmp/frontier_pic_environment.sh",
            "--pre-manifest-attestation",
            "/tmp/pre_manifest_attestation.json",
            "--timeout-margin-artifact",
            "/tmp/timeout_margin.json",
            "--queue-snapshot",
            "/tmp/queue_snapshot.txt",
            "--site-policy-checked-utc",
            "2026-06-02T00:00:00Z",
            "--output-root",
            "/tmp/reviewed-config",
        ]
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                execution.build_parser().parse_args(required)
            with self.assertRaises(SystemExit):
                execution.build_parser().parse_args(
                    [*required, "--case-id", "ps_p0_unknown"]
                )


if __name__ == "__main__":
    unittest.main()
