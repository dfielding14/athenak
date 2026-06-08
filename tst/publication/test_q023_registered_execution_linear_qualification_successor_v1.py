#!/usr/bin/env python3
"""Focused tests for registered Q023 aggregate qualification."""

from __future__ import annotations

import copy
from contextlib import contextmanager
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tst.publication import analyze_q023_paper_bell_linear_joverc as bell
from tst.publication import (
    q023_registered_execution_linear_qualification_successor_v1 as qualification,
)
from tst.publication.test_analyze_q023_paper_bell_linear_joverc import (
    _write_materialized_member,
    _write_registered_manifest_fixture,
)


def _synthetic_admissions() -> list[dict[str, object]]:
    admissions = []
    bundle = bell.synthetic_predecessor_bundle()
    for index, record in enumerate(bundle["records"], 1):
        member_id = str(record["member_id"])
        physics_record = {
            key: record[key] for key in bell._PHYSICS_MATRIX_RECORD_KEYS
        }
        report = bell._analyze_physics_matrix_record(
            physics_record, members=bell._manifest_members()
        )
        execution_identity = {
            "reservation_id": f"reservation-{index:03d}",
            "submission_id": f"submission-{index:03d}",
            "reconciliation_event_sha256": f"{index:064x}",
            "reconciliation_mirror_ack_sha256": f"{index + 100:064x}",
            "control_plane_version": "a" * 64,
            "registered_science_authorization_id": (
                f"q023-linear-{index:03d}-v1"
            ),
            "slurm_job_id": str(1000 + index),
            "pre_submit_manifest_path": f"/manifest/{index:03d}.json",
            "pre_submit_manifest_sha256": f"{index + 200:064x}",
        }
        admission = {
            "schema_version": 1,
            "record_type": qualification.CASE_RECORD_TYPE,
            "successor_id": qualification.SUCCESSOR_ID,
            "campaign_id": bell.CAMPAIGN_ID,
            "status": "registered_execution_case_admitted_non_authorizing",
            "qualification_effect": qualification.QUALIFICATION_EFFECT,
            "member_id": member_id,
            "artifact_root": f"/run/{index:03d}",
            "q043_dependency_sha256": "b" * 64,
            "source_record": record,
            "source_record_sha256": qualification.canonical_sha256(record),
            "physics_record": physics_record,
            "physics_record_sha256": qualification.canonical_sha256(
                physics_record
            ),
            "analysis_report": {
                **report,
                "provenance_kind": "registered_execution_trace",
                "provenance_gate_pass": True,
            },
            "analysis_report_sha256": "",
            "candidate_binding": {
                "source_commit": "c" * 40,
                "source_bundle_sha256": "d" * 64,
                "source_archive_sha256": "e" * 64,
                "clean_candidate_manifest_sha256": "f" * 64,
                "executable_sha256": "1" * 64,
                "environment_sha256": "2" * 64,
                "control_plane_version": "a" * 64,
            },
            "candidate_binding_sha256": "",
            "execution_identity": execution_identity,
            "paired_evidence": {},
            "authorization": dict(qualification.AUTHORIZATION_BOUNDARY),
        }
        admission["analysis_report_sha256"] = qualification.canonical_sha256(
            admission["analysis_report"]
        )
        admission["candidate_binding_sha256"] = qualification.canonical_sha256(
            admission["candidate_binding"]
        )
        admissions.append(admission)
    return admissions


def _write_production_case(
    orion_root: Path,
    project_home_root: Path,
    member: dict[str, object],
    dependency: dict[str, object],
    *,
    physical_times: list[float] | None = None,
) -> Path:
    submission_id = "11111111-1111-4111-8111-111111111111"
    case_root = (
        orion_root / qualification.prep.RUN_NAMESPACE / submission_id
    )
    cycles = [0, *(1000 + index for index in range(1, 89))]
    terminal_cycle = cycles[-1]
    provenance, _ = _write_materialized_member(
        case_root,
        member,
        dependency,
        cycles=cycles,
        physical_times=physical_times,
    )
    manifest_receipt, _ = _write_registered_manifest_fixture(
        orion_root, case_root, member
    )
    rank_count = math.prod(
        int(value) for value in member["decomposition_splits"]
    )
    stdout = case_root / "athena_stdout.txt"
    stdout.write_text(
        (
            f"Q023_REGISTERED_EXECUTION case_id={member['member_id']} "
            f"mpi_world_size={rank_count} "
            f"rank_ids={','.join(str(rank) for rank in range(rank_count))}\n"
            "Q023_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0\n"
            "Terminating on time limit\n"
            f"time=1.750000e+00 cycle={terminal_cycle}\n"
            "tlim=1.750000e+00 nlim=2000000\n"
        ),
        encoding="utf-8",
    )
    raw_inventory = [
        {
            "path": str(artifact["path"]).removeprefix("raw/"),
            "sha256": artifact["sha256"],
            "byte_count": (case_root / str(artifact["path"])).stat().st_size,
            "member_id": member["member_id"],
            "variable": artifact["variable"],
            "output_index": index // len(bell._RAW_OUTPUT_VARIABLES),
        }
        for index, artifact in enumerate(provenance["raw_artifacts"])
    ]
    producer = {
        "entrypoint": "reconcile_q023_registered_execution.py",
        "entrypoint_sha256": "1" * 64,
        "launch_trampoline_sha256": "2" * 64,
        "control_plane_version": manifest_receipt["control_plane_version"],
    }
    action = {
        "action_id": member["member_id"],
        "kind": "athena",
        "resources": {"tasks": rank_count},
        "arguments": [
            {"literal": "-i"},
            {"snapshot_role": "input-deck"},
            {"literal": "-d"},
            {"artifact_directory": "raw"},
        ],
    }
    wrapper = {
        "stdout_sha256": hashlib.sha256(stdout.read_bytes()).hexdigest(),
        "observed_world_size": rank_count,
        "observed_rank_ids": list(range(rank_count)),
        "terminal_cycle": terminal_cycle,
        "terminal_time": bell.LINEAR_RUNTIME_TLIM,
    }
    mirror_root = (
        project_home_root
        / qualification.prep.PROJECT_HOME_RECEIPT_NAMESPACE
        / submission_id
    )
    terminal_mirror = mirror_root / "q023_terminal_receipt.json"
    receipt_mirror = mirror_root / "q023_registered_execution_receipt.json"
    terminal = {
        "schema_version": 1,
        "record_type": "q023_registered_execution_terminal_receipt",
        "campaign_id": bell.CAMPAIGN_ID,
        "member_id": member["member_id"],
        "submission_id": submission_id,
        "slurm_job_id": "123",
        "slurm_terminal_state": "COMPLETED",
        "slurm_exit_code": "0:0",
        "terminal_cycle": terminal_cycle,
        "terminal_time": bell.LINEAR_RUNTIME_TLIM,
        "raw_inventory_sha256": bell._sha256_bytes(
            bell._canonical_json_bytes(raw_inventory)
        ),
    }
    terminal_payload = (
        json.dumps(terminal, sort_keys=True) + "\n"
    ).encode("utf-8")
    receipt = {
        "schema_version": 1,
        "record_type": "q023_reconciled_registered_execution_receipt",
        "receipt_role": "immutable_reconciled_registered_execution",
        "registration_scope": "registered_science",
        "reconciled": True,
        "campaign_id": bell.CAMPAIGN_ID,
        "member_id": member["member_id"],
        "reservation_id": "reservation-001",
        "submission_id": submission_id,
        "reconciliation_event_sha256": "3" * 64,
        "reconciliation_mirror_ack_sha256": "4" * 64,
        "control_plane_version": manifest_receipt["control_plane_version"],
        "project_home_mirrors": {
            "registered_execution_receipt_path": str(receipt_mirror),
            "terminal_receipt_path": str(terminal_mirror),
        },
        "producer": producer,
        "registered_science_authorization_id": (
            manifest_receipt["registered_science_authorization_id"]
        ),
        "source_commit": manifest_receipt["source_commit"],
        "source_bundle_sha256": manifest_receipt["source_bundle_sha256"],
        "source_archive_sha256": manifest_receipt["source_archive_sha256"],
        "clean_candidate_manifest_sha256": manifest_receipt[
            "clean_candidate_manifest_sha256"
        ],
        "executable_sha256": manifest_receipt["executable_sha256"],
        "environment_sha256": "7" * 64,
        "deck_sha256": member["deck_sha256"],
        "command_evidence": {
            "source": "trusted_pre_submit_manifest_and_installed_trampoline",
            "executor": "trusted_trampoline_athena_argv_v1",
            "action": action,
            "launch_trampoline_entrypoint": "launch_trampoline.py",
            "launch_trampoline_sha256": producer["launch_trampoline_sha256"],
            "trusted_wrapper_evidence": wrapper,
        },
        "mpi_evidence": {
            "tasks": rank_count,
            "observed_world_size": rank_count,
            "observed_rank_ids": list(range(rank_count)),
        },
        "slurm_job_id": "123",
        "slurm_terminal_state": "COMPLETED",
        "slurm_exit_code": "0:0",
        "terminal_cycle": terminal_cycle,
        "terminal_time": bell.LINEAR_RUNTIME_TLIM,
        "raw_output_root": str(case_root / "raw"),
        "artifact_dir": str(case_root),
        "artifact_inventory": {},
        "trampoline_completion_receipt": {},
        "terminal_receipt_sha256": hashlib.sha256(terminal_payload).hexdigest(),
        "pre_submit_manifest_path": manifest_receipt[
            "pre_submit_manifest_path"
        ],
        "pre_submit_manifest_sha256": manifest_receipt[
            "pre_submit_manifest_sha256"
        ],
        "raw_inventory": raw_inventory,
        "raw_inventory_sha256": terminal["raw_inventory_sha256"],
    }
    receipt_payload = (
        json.dumps(receipt, sort_keys=True) + "\n"
    ).encode("utf-8")
    analysis_root = case_root / "analysis"
    analysis_root.mkdir(exist_ok=True)
    (analysis_root / "q023_terminal_receipt.json").write_bytes(terminal_payload)
    (
        analysis_root / "q023_registered_execution_receipt.json"
    ).write_bytes(receipt_payload)
    mirror_root.mkdir(parents=True)
    terminal_mirror.write_bytes(terminal_payload)
    receipt_mirror.write_bytes(receipt_payload)
    for path in case_root.rglob("*"):
        path.chmod(0o555 if path.is_dir() else 0o444)
    case_root.chmod(0o555)
    for path in mirror_root.rglob("*"):
        path.chmod(0o555 if path.is_dir() else 0o444)
    mirror_root.chmod(0o555)
    return case_root


def _rederivation_fixture(case_root: Path) -> dict[str, object]:
    receipt_path = case_root / "analysis/q023_registered_execution_receipt.json"
    terminal_path = case_root / "analysis/q023_terminal_receipt.json"
    receipt_payload = receipt_path.read_bytes()
    receipt = json.loads(receipt_payload)
    return {
        "control_plane_version": receipt["control_plane_version"],
        "entrypoint": receipt["producer"]["entrypoint"],
        "entrypoint_sha256": receipt["producer"]["entrypoint_sha256"],
        "reconciliation_event_sha256": receipt[
            "reconciliation_event_sha256"
        ],
        "reconciliation_mirror_ack_sha256": receipt[
            "reconciliation_mirror_ack_sha256"
        ],
        "receipt_sha256": hashlib.sha256(receipt_payload).hexdigest(),
        "terminal_receipt_sha256": hashlib.sha256(
            terminal_path.read_bytes()
        ).hexdigest(),
        "exact_byte_rederivation_passed": True,
    }


class Q023RegisteredExecutionLinearQualificationTests(unittest.TestCase):
    def test_output_index_metadata_must_be_strictly_chronological(self) -> None:
        dependency = bell.synthetic_q043_registered_raw_oracle_dependency()
        dependency.update(
            {
                "binding_kind": "registered_matrix_qualification",
                "registered_admission_digest_bound": True,
                "registered_admission_schema_bound": True,
                "registered_execution_qualification_check_pass": True,
                "registered_raw_oracle_pass": True,
                "complete_foundational_raw_oracle_matrix_pass": True,
            }
        )
        member = next(iter(bell._manifest_members().values()))
        physical_times = [
            *(index * bell.LINEAR_OUTPUT_DT for index in range(88)),
            bell.LINEAR_RUNTIME_TLIM,
        ]
        physical_times[1], physical_times[2] = (
            physical_times[2],
            physical_times[1],
        )
        with (
            tempfile.TemporaryDirectory() as orion_directory,
            tempfile.TemporaryDirectory() as project_home_directory,
        ):
            orion_root = Path(orion_directory).resolve()
            project_home_root = Path(project_home_directory).resolve()
            case_root = _write_production_case(
                orion_root,
                project_home_root,
                member,
                dependency,
                physical_times=physical_times,
            )
            with (
                patch.object(
                    qualification.bell,
                    "validate_q043_dependency",
                    return_value=dependency,
                ),
                patch.object(
                    qualification,
                    "_trusted_reconciliation_and_rederivation",
                    return_value=_rederivation_fixture(case_root),
                ),
            ):
                with self.assertRaisesRegex(
                    qualification.QualificationError,
                    "output-index metadata is not strictly chronological",
                ):
                    qualification.build_case_admission(
                        member_id=str(member["member_id"]),
                        artifact_root=case_root,
                        q043_registered_raw_oracle_dependency=dependency,
                        q043_artifact_root=orion_root,
                        authorized_orion_root=orion_root,
                        authorized_project_home_root=project_home_root,
                    )

    def test_installed_producer_must_rederive_exact_receipt_bytes(self) -> None:
        with (
            tempfile.TemporaryDirectory() as orion_directory,
            tempfile.TemporaryDirectory() as project_home_directory,
        ):
            orion_root = Path(orion_directory).resolve()
            project_home_root = Path(project_home_directory).resolve()
            (orion_root / "ledger").mkdir()
            (project_home_root / "ledger").mkdir()
            receipt_path = orion_root / "receipt.json"
            terminal_path = orion_root / "terminal.json"
            receipt_mirror_path = project_home_root / "receipt.json"
            terminal_mirror_path = project_home_root / "terminal.json"
            receipt_payload = b'{"receipt":true}\n'
            terminal_payload = b'{"terminal":true}\n'
            version = "a" * 64
            producer_sha256 = "b" * 64
            trampoline_sha256 = "c" * 64
            receipt = {
                "control_plane_version": version,
                "producer": {
                    "entrypoint": "reconcile_q023_registered_execution.py",
                    "entrypoint_sha256": producer_sha256,
                    "launch_trampoline_sha256": trampoline_sha256,
                    "control_plane_version": version,
                },
                "reconciliation_event_sha256": "d" * 64,
                "reconciliation_mirror_ack_sha256": "e" * 64,
                "member_id": "member",
                "reservation_id": "reservation",
                "submission_id": "submission",
                "slurm_job_id": "123",
                "artifact_dir": str(orion_root / "run"),
                "source_commit": "f" * 40,
                "clean_candidate_manifest_sha256": "1" * 64,
                "executable_sha256": "2" * 64,
                "pre_submit_manifest_path": str(orion_root / "manifest.json"),
                "pre_submit_manifest_sha256": "3" * 64,
                "registered_science_authorization_id": "q023-linear-001-v1",
            }
            mirror_jsonl = project_home_root / "ledger/node_hours.jsonl"
            event = {
                "event_type": "reconciliation",
                "event_sha256": receipt["reconciliation_event_sha256"],
                "submission_scope": "registered_science",
                "campaign": qualification.prep.CAMPAIGN,
                "test_id": receipt["member_id"],
                "reservation_id": receipt["reservation_id"],
                "submission_id": receipt["submission_id"],
                "job_id": receipt["slurm_job_id"],
                "artifact_dir": receipt["artifact_dir"],
                "git_commit": receipt["source_commit"],
                "clean_candidate_manifest_sha256": receipt[
                    "clean_candidate_manifest_sha256"
                ],
                "executable_sha256": receipt["executable_sha256"],
                "manifest_path": receipt["pre_submit_manifest_path"],
                "manifest_sha256": receipt["pre_submit_manifest_sha256"],
                "control_plane_version": version,
                "reconciled_by_control_plane_version": version,
                "registered_science_authorization_id": receipt[
                    "registered_science_authorization_id"
                ],
                "reconciled": True,
                "state": "COMPLETED",
                "scheduler_exit_code": "0:0",
            }
            mirror_ack = {
                "mirrored_event_sha256": event["event_sha256"],
                "mirror_ack_sha256": receipt[
                    "reconciliation_mirror_ack_sha256"
                ],
                "mirror_transport": "filesystem_copy",
                "mirror_destination": str(mirror_jsonl),
            }

            class FakeLedger:
                @staticmethod
                @contextmanager
                def validated_read_only_mirrored_state_snapshot(*args, **kwargs):
                    del args, kwargs
                    yield [event]

                @staticmethod
                def require_explicit_genesis(records):
                    self.assertEqual(records, [event])

                @staticmethod
                def validate_receipts(*args, **kwargs):
                    del args, kwargs
                    return [mirror_ack]

            expected = (
                terminal_path,
                terminal_payload,
                receipt_path,
                receipt_payload,
                terminal_mirror_path,
                receipt_mirror_path,
            )

            class FakeProducer:
                derived = expected

                @classmethod
                def derive_q023_registered_execution_evidence(cls, *args, **kwargs):
                    del args, kwargs
                    return cls.derived

            pair = {
                "inventory": {"version": version},
                "digests": {
                    "reconcile_q023_registered_execution.py": producer_sha256,
                    "launch_trampoline.py": trampoline_sha256,
                },
            }

            @contextmanager
            def installed_modules(*args, **kwargs):
                del args, kwargs
                yield {
                    "ledger.py": FakeLedger,
                    "reconcile_q023_registered_execution.py": FakeProducer,
                }, pair

            with (
                patch.object(
                    qualification.prep, "AUTHORIZED_ORION_ROOT", orion_root
                ),
                patch.object(
                    qualification.prep,
                    "CANONICAL_PROJECT_HOME_ROOT",
                    project_home_root,
                ),
                patch.object(
                    qualification,
                    "_CANONICAL_PROJECT_HOME_LEDGER_ROOT",
                    project_home_root,
                ),
                patch.object(
                    qualification.bell.q043_registered,
                    "_installed_control_plane_modules",
                    side_effect=installed_modules,
                ),
            ):
                report = qualification._trusted_reconciliation_and_rederivation(
                    receipt=receipt,
                    receipt_path=receipt_path,
                    receipt_payload=receipt_payload,
                    terminal_path=terminal_path,
                    terminal_payload=terminal_payload,
                    receipt_mirror_path=receipt_mirror_path,
                    terminal_mirror_path=terminal_mirror_path,
                    authorized_orion_root=orion_root,
                    authorized_project_home_root=project_home_root,
                )
                self.assertTrue(report["exact_byte_rederivation_passed"])
                FakeProducer.derived = expected[:-1] + (
                    receipt_mirror_path.with_name("substituted.json"),
                )
                with self.assertRaisesRegex(
                    qualification.QualificationError, "different evidence bytes"
                ):
                    qualification._trusted_reconciliation_and_rederivation(
                        receipt=receipt,
                        receipt_path=receipt_path,
                        receipt_payload=receipt_payload,
                        terminal_path=terminal_path,
                        terminal_payload=terminal_payload,
                        receipt_mirror_path=receipt_mirror_path,
                        terminal_mirror_path=terminal_mirror_path,
                        authorized_orion_root=orion_root,
                        authorized_project_home_root=project_home_root,
                    )

    def test_case_admission_is_derived_from_445_retained_raw_files(self) -> None:
        dependency = bell.synthetic_q043_registered_raw_oracle_dependency()
        dependency.update(
            {
                "binding_kind": "registered_matrix_qualification",
                "registered_admission_digest_bound": True,
                "registered_admission_schema_bound": True,
                "registered_execution_qualification_check_pass": True,
                "registered_raw_oracle_pass": True,
                "complete_foundational_raw_oracle_matrix_pass": True,
            }
        )
        member = next(iter(bell._manifest_members().values()))
        with (
            tempfile.TemporaryDirectory() as orion_directory,
            tempfile.TemporaryDirectory() as project_home_directory,
        ):
            orion_root = Path(orion_directory).resolve()
            project_home_root = Path(project_home_directory).resolve()
            case_root = _write_production_case(
                orion_root, project_home_root, member, dependency
            )
            with (
                patch.object(
                    qualification.bell,
                    "validate_q043_dependency",
                    return_value=dependency,
                ),
                patch.object(
                    qualification,
                    "_trusted_reconciliation_and_rederivation",
                    return_value=_rederivation_fixture(case_root),
                ),
            ):
                admission = qualification.build_case_admission(
                    member_id=str(member["member_id"]),
                    artifact_root=case_root,
                    q043_registered_raw_oracle_dependency=dependency,
                    q043_artifact_root=orion_root,
                    authorized_orion_root=orion_root,
                    authorized_project_home_root=project_home_root,
                )
                self.assertEqual(
                    len(admission["source_record"]["provenance"]["raw_artifacts"]),
                    445,
                )
                self.assertEqual(
                    admission["source_record"]["provenance"]["raw_artifacts"][-1][
                        "time"
                    ],
                    bell.LINEAR_RUNTIME_TLIM,
                )
                self.assertEqual(
                    admission["source_record"]["provenance"]["raw_artifacts"][-1][
                        "cycle"
                    ],
                    1088,
                )
                altered = copy.deepcopy(admission)
                altered["source_record"]["physics_trace"]["right_mode_real"][10] += (
                    1.0e-8
                )
                with self.assertRaisesRegex(
                    qualification.QualificationError, "derived fields"
                ):
                    qualification.validate_case_admission(
                        altered,
                        q043_registered_raw_oracle_dependency=dependency,
                        q043_artifact_root=orion_root,
                        authorized_orion_root=orion_root,
                        authorized_project_home_root=project_home_root,
                    )

    def test_complete_registered_matrix_preserves_non_authority(self) -> None:
        admissions = _synthetic_admissions()
        dependency = {"binding_kind": "registered_matrix_qualification"}
        with (
            patch.object(
                qualification.bell,
                "validate_q043_dependency",
                return_value=dependency,
            ),
            patch.object(
                qualification.bell,
                "_dependency_digest",
                return_value="b" * 64,
            ),
            patch.object(
                qualification,
                "validate_case_admission",
                side_effect=lambda value, **_: value,
            ),
        ):
            matrix = qualification.build_matrix_qualification(
                case_admissions=admissions,
                q043_registered_raw_oracle_dependency=dependency,
                q043_artifact_root=qualification.prep.AUTHORIZED_ORION_ROOT,
            )
        self.assertEqual(matrix["case_count"], 55)
        self.assertTrue(matrix["registered_execution_qualification_check_pass"])
        self.assertTrue(matrix["registered_linear_qualification_pass"])
        self.assertFalse(matrix["authorization"]["scientific_claim_authorized"])
        self.assertFalse(matrix["authorization"]["publication_authorized"])

    def test_reused_execution_identity_fails_closed(self) -> None:
        admissions = _synthetic_admissions()
        admissions[1]["execution_identity"]["submission_id"] = admissions[0][
            "execution_identity"
        ]["submission_id"]
        dependency = {"binding_kind": "registered_matrix_qualification"}
        with (
            patch.object(
                qualification.bell,
                "validate_q043_dependency",
                return_value=dependency,
            ),
            patch.object(
                qualification,
                "validate_case_admission",
                side_effect=lambda value, **_: value,
            ),
            patch.object(
                qualification.bell,
                "_dependency_digest",
                return_value="b" * 64,
            ),
        ):
            with self.assertRaisesRegex(
                qualification.QualificationError, "identity reused: submission_id"
            ):
                qualification.build_matrix_qualification(
                    case_admissions=admissions,
                    q043_registered_raw_oracle_dependency=dependency,
                    q043_artifact_root=qualification.prep.AUTHORIZED_ORION_ROOT,
                )

    def test_complete_physics_failure_is_a_valid_failed_qualification(self) -> None:
        admissions = _synthetic_admissions()
        record = admissions[0]["physics_record"]
        trace = record["physics_trace"]
        trace["left_mode_real"][20] = trace["right_mode_real"][20]
        trace["left_mode_imag"][20] = trace["right_mode_imag"][20]
        admissions[0]["physics_record_sha256"] = qualification.canonical_sha256(
            record
        )
        dependency = {"binding_kind": "registered_matrix_qualification"}
        with (
            patch.object(
                qualification.bell,
                "validate_q043_dependency",
                return_value=dependency,
            ),
            patch.object(
                qualification.bell,
                "_dependency_digest",
                return_value="b" * 64,
            ),
            patch.object(
                qualification,
                "validate_case_admission",
                side_effect=lambda value, **_: value,
            ),
        ):
            matrix = qualification.build_matrix_qualification(
                case_admissions=admissions,
                q043_registered_raw_oracle_dependency=dependency,
                q043_artifact_root=qualification.prep.AUTHORIZED_ORION_ROOT,
            )
        self.assertTrue(matrix["registered_execution_qualification_check_pass"])
        self.assertFalse(matrix["registered_linear_qualification_pass"])
        self.assertIn("physics_fail", matrix["status"])

    def test_mixed_candidate_matrix_fails_closed(self) -> None:
        admissions = _synthetic_admissions()
        admissions[-1]["candidate_binding"]["executable_sha256"] = "9" * 64
        admissions[-1]["candidate_binding_sha256"] = (
            qualification.canonical_sha256(admissions[-1]["candidate_binding"])
        )
        dependency = {"binding_kind": "registered_matrix_qualification"}
        with (
            patch.object(
                qualification.bell,
                "validate_q043_dependency",
                return_value=dependency,
            ),
            patch.object(
                qualification.bell,
                "_dependency_digest",
                return_value="b" * 64,
            ),
            patch.object(
                qualification,
                "validate_case_admission",
                side_effect=lambda value, **_: value,
            ),
        ):
            with self.assertRaisesRegex(
                qualification.QualificationError, "mixes clean candidates"
            ):
                qualification.build_matrix_qualification(
                    case_admissions=admissions,
                    q043_registered_raw_oracle_dependency=dependency,
                    q043_artifact_root=qualification.prep.AUTHORIZED_ORION_ROOT,
                )

    def test_authority_drift_changes_canonical_admission(self) -> None:
        admission = _synthetic_admissions()[0]
        altered = copy.deepcopy(admission)
        altered["authorization"]["publication_authorized"] = True
        self.assertNotEqual(
            qualification.canonical_sha256(admission),
            qualification.canonical_sha256(altered),
        )


if __name__ == "__main__":
    unittest.main()
