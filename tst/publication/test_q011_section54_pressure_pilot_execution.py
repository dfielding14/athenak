#!/usr/bin/env python3
"""Fail-closed tests for the Q-011 pressure-pilot registered execution tranche."""

from __future__ import annotations

import copy
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
from tst.publication.frontier_control_plane.operator_attestation import (
    OPERATOR_STATEMENTS,
)


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


def _sealed_attestation(root: Path, case: execution.PressureCase) -> Path:
    now = execution.datetime.now(execution.timezone.utc).replace(microsecond=0)
    archive = root / "operator_attestations"
    archive.mkdir(mode=0o700)
    attestation_root = archive / (
        f"{now.strftime('%Y%m%dT%H%M%SZ')}-{case.authorization_id}-pre_manifest"
    )
    attestation_root.mkdir(mode=0o700)
    manual_values = {
        str(root / "ledger/pending_manual_accounting.json"): "absent",
        str(root / "project_home/ledger/pending_manual_accounting.json"): "absent",
    }
    counts = {
        str(root / "ledger/node_hours.jsonl"): 1,
        str(root / "ledger/mirror_receipts.jsonl"): 1,
        str(root / "project_home/ledger/node_hours.jsonl"): 1,
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
        "same_account_process_snapshot.txt": b"fixture process snapshot\n",
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
        _put(attestation_root / filename, payload, 0o400)
        if (
            filename != "same_account_process_snapshot.txt"
            and not filename.startswith("capture_")
        ):
            _put(attestation_root / f"capture_{filename}", payload, 0o400)

    def record(filename: str, **values: object) -> dict[str, object]:
        return {"path": filename, "sha256": _sha256(payloads[filename]), **values}

    utc = now.isoformat().replace("+00:00", "Z")
    attestation = {
        "schema_version": 1,
        "record_type": (
            "q027_frontier_registered_science_same_account_isolation_attestation"
        ),
        "recorded_utc": utc,
        "sealed_utc": utc,
        "registered_science_authorization_id": case.authorization_id,
        "control_plane_version": "a" * 64,
        "phase": "pre_manifest",
        "same_account_process_snapshot": record("same_account_process_snapshot.txt"),
        "queue_snapshot": record("queue_snapshot.txt"),
        "pending_submission_marker": record(
            "pending_submission_marker.txt", value="absent"
        ),
        "pending_manual_accounting_marker": record(
            "pending_manual_accounting_marker.txt", value="absent", values=manual_values
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
                "sha256": _sha256((attestation_root / filename).read_bytes()),
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
        "operator_statement": OPERATOR_STATEMENTS["pre_manifest"],
    }
    attestation_path = _put(
        attestation_root / "attestation.json",
        (json.dumps(attestation, indent=2, sort_keys=True) + "\n").encode("utf-8"),
        0o400,
    )
    attestation_root.chmod(0o500)
    return attestation_path


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
        queue_snapshot = _put(root / "queue_snapshot.txt", b"", 0o444)
        pre_manifest_attestation = _sealed_attestation(root, case)
        with patch.object(execution, "AUTHORIZED_PIC_ROOT", root), patch.object(
            execution, "AUTHORIZED_PROJECT_HOME_ROOT", root / "project_home"
        ):
            yield {
                "case_id": case.case_id,
                "submission_id": str(uuid.uuid4()),
                "clean_candidate_manifest": clean_manifest,
                "executable": executable,
                "environment_profile": environment,
                "control_plane_version": "a" * 64,
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
        bootstrap = preregistration["seed_timeout_margin_bootstrap"]
        self.assertEqual(
            bootstrap["classification"],
            "engineering_seed_margin_bootstrap_only_not_measurement",
        )
        self.assertEqual(bootstrap["athena_walltime_seconds"], 600)
        self.assertEqual(bootstrap["scheduler_walltime_seconds"], 900)
        source_bindings = preregistration["source_bindings"]
        self.assertEqual(
            source_bindings["analysis_scripts"],
            [
                {
                    "path": "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
                    "sha256": "7486c072d38c6c5ff78eebda0b913543327792f3aea38c7dca3b4a3d709b6d8f",
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
                "sha256": "da82f9d1e213f7725336312e0db2cabc7b63e747e46f4f9ea13d46d9379a0624",
            },
        )
        self.assertNotIn(
            "q011_section54_pressure_pilot_registered_execution_retry_successor_v2_2026-06-02.json",
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

    def test_source_preregistration_rejects_contract_drift(self) -> None:
        path, payload, preregistration = execution._read_json(
            execution.EXECUTION_PREREGISTRATION,
            label="Q-011 pressure-pilot registered-execution preregistration",
        )
        for variant in ("authorization", "record_type", "extra_field"):
            with self.subTest(variant=variant):
                changed = copy.deepcopy(preregistration)
                if variant == "authorization":
                    changed["execution_boundary"][
                        "frontier_execution_authorized_by_this_record"
                    ] = True
                elif variant == "record_type":
                    changed["record_type"] = "forged"
                else:
                    changed["unsupported_extra"] = "forged"
                with patch.object(
                    execution, "_read_json", return_value=(path, payload, changed)
                ), self.assertRaisesRegex(
                    execution.ContractError, "preregistration contract drifted"
                ):
                    execution.validate_source_tranche()

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

    def test_seed_timeout_margin_is_explicit_bootstrap_not_measurement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            environment = _put(
                root / "installed" / "frontier_pic_environment.sh",
                execution.ENVIRONMENT_PROFILE_SOURCE.read_bytes(),
                0o444,
            )
            output = root / "seed-timeout"
            rationale = execution.write_seed_timeout_margin(
                output,
                case_id=execution.CASES[0].case_id,
                environment_profile=environment,
                materialized_utc="2026-06-02T00:00:00Z",
                expires_utc="2026-06-02T01:00:00Z",
            )
            self.assertEqual(
                rationale["classification"],
                "engineering_seed_margin_bootstrap_only_not_measurement",
            )
            self.assertIn(
                "not an empirical runtime measurement",
                rationale["controller_field_semantics"]["measured_utc"],
            )
            timeout = json.loads(
                (output / execution.SEED_TIMEOUT_FILENAME).read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                set(timeout),
                execution._TIMEOUT_MARGIN_KEYS,
            )
            self.assertEqual(timeout["athena_walltime_seconds"], 600)
            self.assertEqual(timeout["scheduler_walltime_seconds"], 900)
            self.assertEqual(timeout["measured_utc"], "2026-06-02T00:00:00Z")
            self.assertEqual(
                rationale["timeout_margin_artifact"]["sha256"],
                _sha256((output / execution.SEED_TIMEOUT_FILENAME).read_bytes()),
            )
            for path in [output, *output.iterdir()]:
                self.assertEqual(path.stat().st_mode & 0o222, 0)
            with self.assertRaisesRegex(execution.ContractError, "already exists"):
                execution.write_seed_timeout_margin(
                    output,
                    case_id=execution.CASES[0].case_id,
                    environment_profile=environment,
                    materialized_utc="2026-06-02T00:00:00Z",
                    expires_utc="2026-06-02T01:00:00Z",
                )
            with self.assertRaisesRegex(execution.ContractError, "exceeds four hours"):
                execution.write_seed_timeout_margin(
                    root / "too-long",
                    case_id=execution.CASES[0].case_id,
                    environment_profile=environment,
                    materialized_utc="2026-06-02T00:00:00Z",
                    expires_utc="2026-06-02T04:00:01Z",
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
                    "7486c072d38c6c5ff78eebda0b913543327792f3aea38c7dca3b4a3d709b6d8f",
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
                    self.assertEqual(
                        config["pre_manifest_attestation"],
                        str(binding["pre_manifest_attestation"]),
                    )
                    self.assertEqual(config["physical_mode"], execution.PHYSICAL_MODE)
                    self.assertEqual(
                        config["artifact_dir"],
                        str(
                            root
                            / case.case_id
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
                    with patch.object(
                        execution, "AUTHORIZED_PIC_ROOT", root / case.case_id
                    ), patch.object(
                        execution,
                        "AUTHORIZED_PROJECT_HOME_ROOT",
                        root / case.case_id / "project_home",
                    ):
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
                    self.assertEqual(len(manifest["required_next_steps"]), 9)
                    with patch.object(
                        execution, "AUTHORIZED_PIC_ROOT", root / case.case_id
                    ), patch.object(
                        execution,
                        "AUTHORIZED_PROJECT_HOME_ROOT",
                        root / case.case_id / "project_home",
                    ):
                        with self.assertRaisesRegex(
                            execution.ContractError, "already exists"
                        ):
                            execution.write_reviewed_pre_submit_config(output, **binding)

    def test_pre_manifest_attestation_rejects_legacy_two_field_json(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self._final_binding(root) as binding:
                path = Path(binding["pre_manifest_attestation"])
                path.parent.chmod(0o700)
                path.chmod(0o600)
                path.write_text(
                    json.dumps(
                        {
                            "phase": "pre_manifest",
                            "registered_science_authorization_id": (
                                execution.CASES[0].authorization_id
                            ),
                        }
                    )
                    + "\n",
                    encoding="utf-8",
                )
                path.chmod(0o400)
                path.parent.chmod(0o500)
                with self.assertRaisesRegex(
                    execution.ContractError, "pre-manifest attestation is invalid"
                ):
                    execution.materialize_reviewed_pre_submit_config(**binding)

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

    def test_complete_policy_successors_preserve_baseline_and_replace_only_launch_fields(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline_value = {
                "schema_version": 1,
                "frontier": {"preserved": "frontier"},
                "science_submission_freeze": {
                    "status": "pending_clean_candidate_freeze"
                },
                "registered_science_slices": [],
                "frontier_admission_smoke": {"preserved": "smoke"},
                "olcf_side_storage": {
                    "installed_control_plane_version": "0" * 64,
                    "staged_control_plane_candidate_version": "0" * 64,
                    "last_preflight_utc": "2026-06-01T00:00:00Z",
                    "preserved": "storage",
                },
                "long_term_storage": {"preserved": "retention"},
            }
            baseline = _put(
                root / "baseline.json",
                (json.dumps(baseline_value) + "\n").encode("utf-8"),
                0o444,
            )
            successor = execution.materialize_baseline_policy_successor(
                baseline_policy=baseline,
                control_plane_version="a" * 64,
                last_preflight_utc="2026-06-02T01:02:03Z",
            )
            self.assertEqual(successor["frontier"], baseline_value["frontier"])
            self.assertEqual(
                successor["long_term_storage"], baseline_value["long_term_storage"]
            )
            self.assertEqual(successor["registered_science_slices"], [])
            self.assertEqual(
                successor["olcf_side_storage"]["installed_control_plane_version"],
                "a" * 64,
            )
            self.assertEqual(
                successor["olcf_side_storage"]["preserved"], "storage"
            )
            with self._final_binding(root / "binding") as binding:
                pilot = execution.materialize_pilot_policy_successor(
                    baseline_policy=baseline,
                    control_plane_version="a" * 64,
                    last_preflight_utc="2026-06-02T01:02:03Z",
                    clean_candidate_manifest=binding["clean_candidate_manifest"],
                    executable=binding["executable"],
                    environment_profile=binding["environment_profile"],
                )
            self.assertEqual(pilot["frontier"], baseline_value["frontier"])
            self.assertEqual(
                pilot["long_term_storage"], baseline_value["long_term_storage"]
            )
            self.assertEqual(pilot["science_submission_freeze"]["status"], "authorized")
            self.assertEqual(len(pilot["registered_science_slices"]), 4)

    def test_baseline_policy_successor_rejects_existing_registered_allowlist(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            baseline = _put(
                Path(directory) / "baseline.json",
                (
                    json.dumps(
                        {
                            "registered_science_slices": [{"unexpected": "slice"}],
                            "olcf_side_storage": {},
                        }
                    )
                    + "\n"
                ).encode("utf-8"),
                0o444,
            )
            with self.assertRaisesRegex(execution.ContractError, "empty"):
                execution.materialize_baseline_policy_successor(
                    baseline_policy=baseline,
                    control_plane_version="a" * 64,
                    last_preflight_utc="2026-06-02T01:02:03Z",
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
                queue_snapshot = Path(binding["queue_snapshot"])
                queue_snapshot.chmod(0o644)
                queue_snapshot.write_text(
                    "1|AST207|batch|debug|RUNNING|job|comment\n",
                    encoding="utf-8",
                )
                queue_snapshot.chmod(0o444)
                with self.assertRaisesRegex(execution.ContractError, "six-field"):
                    execution.materialize_reviewed_pre_submit_config(**binding)
            with self._final_binding(root / "fourth") as binding:
                Path(binding["queue_snapshot"]).chmod(0o644)
                with self.assertRaisesRegex(execution.ContractError, "must be read-only"):
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
