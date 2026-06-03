#!/usr/bin/env python3
"""Fail-closed tests for the Q-011 pressure-pilot registered execution tranche."""

from __future__ import annotations

import copy
from contextlib import contextmanager, redirect_stderr, redirect_stdout
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tarfile
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


def _storage_preflight_binding(root: Path, completed_utc: str) -> Path:
    probe_id = "12345678-1234-4234-8234-123456789abc"
    return _put(
        root
        / f"storage-preflight-{completed_utc.replace(':', '')}-{uuid.uuid4()}.json",
        (
            json.dumps(
                {
                    "last_preflight_utc": completed_utc,
                    "orion_simulation_root_preflight": {
                        "method": "local_create_write_sync_remove_probe",
                        "path": str(execution.AUTHORIZED_PIC_ROOT),
                        "status": "passed",
                    },
                    "project_home_preflight": {
                        "method": "local_create_write_sync_remove_probe",
                        "path": str(execution.AUTHORIZED_PROJECT_HOME_ROOT),
                        "status": "passed",
                    },
                    "storage_preflight_evidence": {
                        "orion_path": str(
                            execution.AUTHORIZED_PIC_ROOT
                            / "policy"
                            / "storage_preflight_evidence"
                            / f"{probe_id}.json"
                        ),
                        "probe_id": probe_id,
                        "project_home_path": str(
                            execution.AUTHORIZED_PROJECT_HOME_ROOT
                            / "policy"
                            / "storage_preflight_evidence"
                            / f"{probe_id}.json"
                        ),
                        "sha256": "a" * 64,
                    },
                }
            )
            + "\n"
        ).encode("utf-8"),
        0o444,
    )


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
        self,
        root: Path,
        *,
        case: execution.PressureCase = execution.CASES[0],
        git_commit: str = "a" * 40,
    ) -> Iterator[dict[str, object]]:
        root.mkdir(parents=True, exist_ok=True)
        freeze = root / str(uuid.uuid4())
        freeze.mkdir()
        executable = _put(freeze / "athena", b"exact-clean-athena\n", 0o555)
        executable_sha256 = _sha256(executable.read_bytes())
        archive_payload = io.BytesIO()
        generator_payload = execution.GENERATOR_SOURCE.read_bytes()
        with tarfile.open(fileobj=archive_payload, mode="w") as archive:
            member = tarfile.TarInfo("src/pgen/tests/pic_parallel_shock.cpp")
            member.size = len(generator_payload)
            archive.addfile(member, io.BytesIO(generator_payload))
        source_archive = _put(freeze / "source.tar", archive_payload.getvalue(), 0o444)
        source_archive_sha256 = _sha256(source_archive.read_bytes())
        manifest = {
            "schema_version": 4,
            "freeze_id": freeze.name,
            "source": {
                "archive_path": str(source_archive),
                "archive_sha256": source_archive_sha256,
                "git_commit": git_commit,
            },
            "build": {
                "executable_path": str(executable),
                "executable_sha256": executable_sha256,
                "source_archive_sha256": source_archive_sha256,
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
        now = execution.datetime.now(execution.timezone.utc).replace(microsecond=0)
        timeout = {
            "athena_walltime_seconds": 840,
            "scheduler_walltime_seconds": 900,
            "environment_profile_sha256": environment_sha256,
            "measured_utc": (
                now - execution.timedelta(seconds=1)
            ).isoformat().replace("+00:00", "Z"),
            "expires_utc": (
                now + execution.timedelta(hours=1)
            ).isoformat().replace("+00:00", "Z"),
        }
        timeout_margin = _put(
            root / "timeout_margin.json",
            (json.dumps(timeout) + "\n").encode("utf-8"),
            0o644,
        )
        queue_snapshot = _put(root / "queue_snapshot.txt", b"", 0o444)
        pre_manifest_attestation = _sealed_attestation(root, case)
        historical_sources = execution._historical_v2_preregistration()[
            "source_bindings"
        ]
        with patch.object(
            execution, "source_bindings", return_value=historical_sources
        ), patch.object(
            execution,
            "validate_source_tranche",
            return_value=execution._historical_v2_preregistration(),
        ), patch.object(execution, "AUTHORIZED_PIC_ROOT", root), patch.object(
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

    @staticmethod
    def _prior_case_bindings(case: execution.PressureCase) -> tuple[str, ...]:
        return tuple(
            f"{prior.case_id}={uuid.uuid4()}={'0' * 64}"
            for prior in execution.CASES[: execution.CASES.index(case)]
        )

    @staticmethod
    def _fake_prior_case_closure(
        case: execution.PressureCase,
        submission_id: str,
        descriptor_sha256: str,
        final: dict[str, str],
    ) -> dict[str, str]:
        return {
            "case_id": case.case_id,
            "submission_id": submission_id,
            "artifact_dir": f"/fixture/{case.campaign}/{submission_id}",
            "descriptor_path": f"/fixture/{case.campaign}/{submission_id}/analysis/analysis.json",
            "descriptor_sha256": descriptor_sha256,
            "reconciliation_event_sha256": "1" * 64,
        }

    def test_source_preregistration_and_shared_directive_only_template(self) -> None:
        preregistration = execution._historical_v2_preregistration()
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
            source_bindings["generator_source"],
            {
                "path": "src/pgen/tests/pic_parallel_shock.cpp",
                "sha256": "c0a01e4960f4ebb1a96bedc61fd918b4f7fc76addb59e0f9db3f9ab35eb8c9f9",
            },
        )
        self.assertEqual(
            source_bindings["analysis_scripts"],
            [
                {
                    "path": "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
                    "sha256": "8596d95c9b8952dcb760b10fbe00bf7ba8a2c895713cc3d36dd1efa1aac11ecc",
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
                "sha256": "774e61be728c7bdd90f6a9d264e7ce4e5413d5f2a898608ac726ba5b8baadd08",
            },
        )
        status = execution.historical_v2_source_tranche_status()
        self.assertEqual(status["state"], "historical_consumed_slice_non_authorizing")
        self.assertFalse(status["source_bindings_match_current_checkout"])
        self.assertEqual(status["launch_reauthorization_effect"], "none")
        self.assertFalse(status["consumed_slice_reauthorization_allowed"])
        with self.assertRaisesRegex(
            execution.ContractError, "historical v2 registered-execution tranche is consumed"
        ):
            execution.validate_source_tranche()
        with patch.object(
            execution, "source_bindings", return_value=source_bindings
        ), self.assertRaisesRegex(
            execution.ContractError, "historical v2 registered-execution tranche is consumed"
        ):
            execution.validate_source_tranche()
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
                    "8596d95c9b8952dcb760b10fbe00bf7ba8a2c895713cc3d36dd1efa1aac11ecc",
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
                        binding["prior_case_closures"] = self._prior_case_bindings(case)
                        with patch.object(
                            execution,
                            "_verified_prior_case_closure",
                            side_effect=self._fake_prior_case_closure,
                        ):
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
                    historical_sources = execution._historical_v2_preregistration()[
                        "source_bindings"
                    ]
                    with patch.object(
                        execution, "source_bindings", return_value=historical_sources
                    ), patch.object(
                        execution,
                        "validate_source_tranche",
                        return_value=execution._historical_v2_preregistration(),
                    ), patch.object(
                        execution, "AUTHORIZED_PIC_ROOT", root / case.case_id
                    ), patch.object(
                        execution,
                        "AUTHORIZED_PROJECT_HOME_ROOT",
                        root / case.case_id / "project_home",
                    ), patch.object(
                        execution,
                        "_verified_prior_case_closure",
                        side_effect=self._fake_prior_case_closure,
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
                        execution, "source_bindings", return_value=historical_sources
                    ), patch.object(
                        execution,
                        "validate_source_tranche",
                        return_value=execution._historical_v2_preregistration(),
                    ), patch.object(
                        execution, "AUTHORIZED_PIC_ROOT", root / case.case_id
                    ), patch.object(
                        execution,
                        "AUTHORIZED_PROJECT_HOME_ROOT",
                        root / case.case_id / "project_home",
                    ), patch.object(
                        execution,
                        "_verified_prior_case_closure",
                        side_effect=self._fake_prior_case_closure,
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
                storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T01:02:03Z"),
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
                    storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T01:02:03Z"),
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
            root = Path(directory)
            baseline = _put(
                root / "baseline.json",
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
                    storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T01:02:03Z"),
                )

    def test_policy_successor_rejects_malformed_or_writable_storage_preflight_binding(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = _put(
                root / "baseline.json",
                (
                    json.dumps(
                        {
                            "registered_science_slices": [],
                            "olcf_side_storage": {
                                "installed_control_plane_version": "0" * 64,
                                "staged_control_plane_candidate_version": "0" * 64,
                            },
                        }
                    )
                    + "\n"
                ).encode("utf-8"),
                0o444,
            )
            binding = _storage_preflight_binding(root, "2026-06-02T01:02:03Z")
            binding.chmod(0o644)
            with self.assertRaisesRegex(execution.ContractError, "must be read-only"):
                execution.materialize_baseline_policy_successor(
                    baseline_policy=baseline,
                    control_plane_version="a" * 64,
                    storage_preflight_binding=binding,
                )
            value = json.loads(binding.read_text(encoding="utf-8"))
            value["project_home_preflight"].pop("path")
            binding.write_text(json.dumps(value), encoding="utf-8")
            binding.chmod(0o444)
            with self.assertRaisesRegex(execution.ContractError, "is malformed"):
                execution.materialize_baseline_policy_successor(
                    baseline_policy=baseline,
                    control_plane_version="a" * 64,
                    storage_preflight_binding=binding,
                )

    def test_retire_consumed_slices_successor_is_exact_additive_and_exclusive(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            historical_value = json.loads(
                (execution.READINESS_ROOT / "storage_policy.json").read_text(
                    encoding="utf-8"
                )
            )
            baseline = _put(
                root / "historical-policy.json",
                (json.dumps(historical_value) + "\n").encode("utf-8"),
                0o444,
            )
            output = root / "retired-baseline-successor.json"
            successor = (
                execution.write_retire_consumed_slices_baseline_policy_successor(
                    output,
                    baseline_policy=baseline,
                    control_plane_version="a" * 64,
                    storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:36Z"),
                )
            )
            self.assertEqual(successor["frontier"], historical_value["frontier"])
            self.assertEqual(
                successor["long_term_storage"], historical_value["long_term_storage"]
            )
            self.assertEqual(successor["registered_science_slices"], [])
            self.assertEqual(
                successor["science_submission_freeze"],
                {"status": "pending_clean_candidate_freeze"},
            )
            self.assertEqual(
                successor["olcf_side_storage"]["installed_control_plane_version"],
                "a" * 64,
            )
            self.assertEqual(
                successor["olcf_side_storage"][
                    "staged_control_plane_candidate_version"
                ],
                "a" * 64,
            )
            self.assertEqual(
                successor["olcf_side_storage"]["last_preflight_utc"],
                "2026-06-02T05:08:36Z",
            )
            self.assertEqual(
                json.loads(baseline.read_text(encoding="utf-8")), historical_value
            )
            self.assertEqual(output.stat().st_mode & 0o222, 0)
            with self.assertRaisesRegex(execution.ContractError, "refusing to overwrite"):
                execution.write_retire_consumed_slices_baseline_policy_successor(
                    output,
                    baseline_policy=baseline,
                    control_plane_version="a" * 64,
                    storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:36Z"),
                )
            with self.assertRaisesRegex(execution.ContractError, "newer than baseline"):
                execution.materialize_retire_consumed_slices_baseline_policy_successor(
                    baseline_policy=baseline,
                    control_plane_version="a" * 64,
                    storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:35Z"),
                )
            with self.assertRaisesRegex(execution.ContractError, "must be new"):
                execution.materialize_retire_consumed_slices_baseline_policy_successor(
                    baseline_policy=baseline,
                    control_plane_version=historical_value["olcf_side_storage"][
                        "installed_control_plane_version"
                    ],
                    storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:36Z"),
                )

    def test_retire_consumed_slices_refuses_arbitrary_nonempty_policies(self) -> None:
        historical_value = json.loads(
            (execution.READINESS_ROOT / "storage_policy.json").read_text(
                encoding="utf-8"
            )
        )
        variants = {
            "arbitrary_allowlist": [{"unexpected": "slice"}],
            "drifted_consumed_allowlist": copy.deepcopy(
                historical_value["registered_science_slices"]
            ),
        }
        variants["drifted_consumed_allowlist"][0]["status"] = "consumed"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, slices in variants.items():
                with self.subTest(name=name):
                    changed = copy.deepcopy(historical_value)
                    changed["registered_science_slices"] = slices
                    baseline = _put(
                        root / f"{name}.json",
                        (json.dumps(changed) + "\n").encode("utf-8"),
                        0o444,
                    )
                    with self.assertRaisesRegex(
                        execution.ContractError, "exactly the consumed Q011"
                    ):
                        materialize = getattr(
                            execution,
                            "materialize_retire_consumed_slices_"
                            "baseline_policy_successor",
                        )
                        materialize(
                            baseline_policy=baseline,
                            control_plane_version="a" * 64,
                            storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:36Z"),
                        )
            changed = copy.deepcopy(historical_value)
            changed["science_submission_freeze"]["manifest_sha256"] = "0" * 64
            baseline = _put(
                root / "drifted-freeze.json",
                (json.dumps(changed) + "\n").encode("utf-8"),
                0o444,
            )
            with self.assertRaisesRegex(execution.ContractError, "exact consumed Q011"):
                execution.materialize_retire_consumed_slices_baseline_policy_successor(
                    baseline_policy=baseline,
                    control_plane_version="a" * 64,
                    storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:36Z"),
                )

    def test_candidate_only_successor_authorizes_freeze_without_launch_slices(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            historical_value = json.loads(
                (execution.READINESS_ROOT / "storage_policy.json").read_text(
                    encoding="utf-8"
                )
            )
            historical = _put(
                root / "historical-policy.json",
                (json.dumps(historical_value) + "\n").encode("utf-8"),
                0o444,
            )
            retired_value = (
                execution.materialize_retire_consumed_slices_baseline_policy_successor(
                    baseline_policy=historical,
                    control_plane_version="a" * 64,
                    storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:36Z"),
                )
            )
            retired = _put(
                root / "retired-baseline.json",
                (json.dumps(retired_value) + "\n").encode("utf-8"),
                0o444,
            )
            with self._final_binding(root / "binding") as binding:
                output = root / "candidate-only-successor.json"
                with patch.object(
                    execution,
                    "validate_source_tranche",
                    side_effect=AssertionError("candidate-only reopened consumed slices"),
                ):
                    successor = execution.write_candidate_only_policy_successor(
                        output,
                        baseline_policy=retired,
                        control_plane_version="a" * 64,
                        storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:37Z"),
                        clean_candidate_manifest=binding[
                            "clean_candidate_manifest"
                        ],
                        executable=binding["executable"],
                        environment_profile=binding["environment_profile"],
                    )
                with self.assertRaisesRegex(execution.ContractError, "must match"):
                    execution.materialize_candidate_only_policy_successor(
                        baseline_policy=retired,
                        control_plane_version="b" * 64,
                        storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:37Z"),
                        clean_candidate_manifest=binding[
                            "clean_candidate_manifest"
                        ],
                        executable=binding["executable"],
                        environment_profile=binding["environment_profile"],
                    )
                with patch.object(
                    execution,
                    "_bound_candidate_artifacts",
                    return_value={
                        "clean_candidate_manifest_sha256": (
                            execution
                            ._CONSUMED_HISTORICAL_V2_CLEAN_CANDIDATE_MANIFEST_SHA256
                        )
                    },
                ), self.assertRaisesRegex(execution.ContractError, "must be fresh"):
                    execution.materialize_candidate_only_policy_successor(
                        baseline_policy=retired,
                        control_plane_version="a" * 64,
                        storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:37Z"),
                        clean_candidate_manifest=binding[
                            "clean_candidate_manifest"
                        ],
                        executable=binding["executable"],
                        environment_profile=binding["environment_profile"],
                    )
                self.assertEqual(successor["registered_science_slices"], [])
                self.assertEqual(
                    successor["science_submission_freeze"]["status"], "authorized"
                )
                self.assertEqual(
                    successor["science_submission_freeze"]["manifest_path"],
                    str(binding["clean_candidate_manifest"]),
                )
                self.assertEqual(
                    successor["science_submission_freeze"]["manifest_sha256"],
                    _sha256(Path(binding["clean_candidate_manifest"]).read_bytes()),
                )
                self.assertEqual(
                    successor["science_submission_freeze"][
                        "build_profile_control_plane_version"
                    ],
                    "a" * 64,
                )
                self.assertEqual(
                    successor["olcf_side_storage"]["installed_control_plane_version"],
                    "a" * 64,
                )
                self.assertEqual(
                    successor["olcf_side_storage"][
                        "staged_control_plane_candidate_version"
                    ],
                    "a" * 64,
                )
                self.assertEqual(output.stat().st_mode & 0o222, 0)
                detached = _put(root / "detached-athena", b"exact-clean-athena\n", 0o555)
                with self.assertRaisesRegex(execution.ContractError, "not adjacent"):
                    execution.materialize_candidate_only_policy_successor(
                        baseline_policy=retired,
                        control_plane_version="a" * 64,
                        storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:37Z"),
                        clean_candidate_manifest=binding[
                            "clean_candidate_manifest"
                        ],
                        executable=detached,
                        environment_profile=binding["environment_profile"],
                    )
                equal_preflight = (
                    execution.materialize_candidate_only_policy_successor(
                        baseline_policy=retired,
                        control_plane_version="a" * 64,
                        storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:36Z"),
                        clean_candidate_manifest=binding[
                            "clean_candidate_manifest"
                        ],
                        executable=binding["executable"],
                        environment_profile=binding["environment_profile"],
                    )
                )
                self.assertEqual(
                    equal_preflight["olcf_side_storage"]["last_preflight_utc"],
                    "2026-06-02T05:08:36Z",
                )
                with self.assertRaisesRegex(
                    execution.ContractError, "equal to or newer"
                ):
                    execution.materialize_candidate_only_policy_successor(
                        baseline_policy=retired,
                        control_plane_version="a" * 64,
                        storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:35Z"),
                        clean_candidate_manifest=binding[
                            "clean_candidate_manifest"
                        ],
                        executable=binding["executable"],
                        environment_profile=binding["environment_profile"],
                    )

    def test_candidate_only_successor_requires_pending_empty_baseline(self) -> None:
        baseline_value = {
            "science_submission_freeze": {"status": "pending_clean_candidate_freeze"},
            "registered_science_slices": [],
            "olcf_side_storage": {"last_preflight_utc": "2026-06-02T05:08:36Z"},
        }
        variants = {
            "nonempty": {
                **baseline_value,
                "registered_science_slices": [{"unexpected": "slice"}],
            },
            "authorized": {
                **baseline_value,
                "science_submission_freeze": {"status": "authorized"},
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, value in variants.items():
                with self.subTest(name=name):
                    baseline = _put(
                        root / f"{name}.json",
                        (json.dumps(value) + "\n").encode("utf-8"),
                        0o444,
                    )
                    message = "empty" if name == "nonempty" else "pending"
                    with self.assertRaisesRegex(execution.ContractError, message):
                        execution.materialize_candidate_only_policy_successor(
                            baseline_policy=baseline,
                            control_plane_version="b" * 64,
                            storage_preflight_binding=_storage_preflight_binding(root, "2026-06-02T05:08:37Z"),
                            clean_candidate_manifest=root / "unused-manifest.json",
                            executable=root / "unused-athena",
                            environment_profile=root / "unused-environment.sh",
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

    def test_final_bindings_reject_failed_v1_and_unrepaired_source_archives(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self._final_binding(
                root / "failed", git_commit=execution._FAILED_V1_GIT_COMMIT
            ) as binding:
                with self.assertRaisesRegex(execution.ContractError, "failed v1 carrier"):
                    execution.materialize_registered_science_slices(
                        clean_candidate_manifest=binding["clean_candidate_manifest"],
                        executable=binding["executable"],
                        environment_profile=binding["environment_profile"],
                    )
            with self._final_binding(root / "unrepaired") as binding:
                archive_payload = io.BytesIO()
                with tarfile.open(fileobj=archive_payload, mode="w") as archive:
                    payload = b"unrepaired generator\n"
                    member = tarfile.TarInfo("src/pgen/tests/pic_parallel_shock.cpp")
                    member.size = len(payload)
                    archive.addfile(member, io.BytesIO(payload))
                source_archive = Path(binding["clean_candidate_manifest"]).parent / "source.tar"
                source_archive.chmod(0o644)
                source_archive.write_bytes(archive_payload.getvalue())
                source_archive.chmod(0o444)
                manifest_path = Path(binding["clean_candidate_manifest"])
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                digest = _sha256(source_archive.read_bytes())
                manifest["source"]["archive_sha256"] = digest
                manifest["build"]["source_archive_sha256"] = digest
                manifest_path.chmod(0o644)
                manifest_path.write_text(json.dumps(manifest) + "\n", encoding="utf-8")
                manifest_path.chmod(0o444)
                with self.assertRaisesRegex(execution.ContractError, "repaired generator"):
                    execution.materialize_registered_science_slices(
                        clean_candidate_manifest=manifest_path,
                        executable=binding["executable"],
                        environment_profile=binding["environment_profile"],
                    )

    def test_later_case_config_requires_exact_ordered_prior_case_closures(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case = execution.CASES[2]
            with self._final_binding(root, case=case) as binding:
                with self.assertRaisesRegex(execution.ContractError, "requires exactly 2"):
                    execution.materialize_reviewed_pre_submit_config(**binding)
                reversed_bindings = tuple(reversed(self._prior_case_bindings(case)))
                with self.assertRaisesRegex(execution.ContractError, "exact preregistered order"):
                    execution.materialize_reviewed_pre_submit_config(
                        **binding, prior_case_closures=reversed_bindings
                    )
                bindings = self._prior_case_bindings(case)
                with patch.object(
                    execution,
                    "_verified_prior_case_closure",
                    side_effect=self._fake_prior_case_closure,
                ):
                    config = execution.materialize_reviewed_pre_submit_config(
                        **binding, prior_case_closures=bindings
                    )
                self.assertEqual(
                    [record["case_id"] for record in config["prior_case_closures"]],
                    ["ps_p0_1p00", "ps_p0_0p05"],
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
            with self.assertRaisesRegex(
                execution.ContractError,
                "historical v2 registered-execution tranche is consumed",
            ):
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

    def test_policy_transition_cli_dispatches_reviewed_successor_writers(self) -> None:
        with patch.object(
            execution,
            "write_retire_consumed_slices_baseline_policy_successor",
            return_value={"transition": "retired"},
        ) as retirement, patch.object(
            sys,
            "argv",
            [
                "q011",
                "retire-consumed-slices-baseline-policy-successor",
                "--baseline-policy",
                "/tmp/historical-policy.json",
                "--control-plane-version",
                "a" * 64,
                "--storage-preflight-binding",
                "/tmp/storage-preflight.json",
                "--output",
                "/tmp/retired-policy.json",
            ],
        ), redirect_stdout(io.StringIO()):
            execution.main()
        retirement.assert_called_once_with(
            Path("/tmp/retired-policy.json"),
            baseline_policy=Path("/tmp/historical-policy.json"),
            control_plane_version="a" * 64,
            storage_preflight_binding=Path("/tmp/storage-preflight.json"),
        )
        with patch.object(
            execution,
            "write_candidate_only_policy_successor",
            return_value={"transition": "candidate-only"},
        ) as candidate, patch.object(
            sys,
            "argv",
            [
                "q011",
                "candidate-only-policy-successor",
                "--baseline-policy",
                "/tmp/retired-policy.json",
                "--control-plane-version",
                "b" * 64,
                "--storage-preflight-binding",
                "/tmp/storage-preflight.json",
                "--clean-candidate-manifest",
                "/tmp/clean_candidate_manifest.json",
                "--executable",
                "/tmp/athena",
                "--environment-profile",
                "/tmp/frontier_pic_environment.sh",
                "--output",
                "/tmp/candidate-only-policy.json",
            ],
        ), redirect_stdout(io.StringIO()):
            execution.main()
        candidate.assert_called_once_with(
            Path("/tmp/candidate-only-policy.json"),
            baseline_policy=Path("/tmp/retired-policy.json"),
            control_plane_version="b" * 64,
            storage_preflight_binding=Path("/tmp/storage-preflight.json"),
            clean_candidate_manifest=Path("/tmp/clean_candidate_manifest.json"),
            executable=Path("/tmp/athena"),
            environment_profile=Path("/tmp/frontier_pic_environment.sh"),
        )


if __name__ == "__main__":
    unittest.main()
