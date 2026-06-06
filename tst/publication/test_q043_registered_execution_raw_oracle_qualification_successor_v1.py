#!/usr/bin/env python3
"""Filesystem-backed adversarial tests for Q043 registered admission."""

from __future__ import annotations

import copy
import hashlib
import io
import json
import math
import os
from pathlib import Path
import struct
import tarfile
import tempfile
import unittest
import uuid

import numpy as np

from tst.publication import (
    q043_bell_current_volume_aware_deposited_current_oracle as oracle,
)
from tst.publication import (
    q043_registered_execution_raw_oracle_qualification_successor_v1 as successor,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q043_registered_execution_raw_oracle_qualification_successor_v1_2026-06-06.json"
)


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _write(path: Path, payload: bytes, *, executable: bool = False) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    if executable:
        path.chmod(0o755)
    return {"path": str(path), "sha256": _sha(payload), "byte_count": len(payload)}


def _source_archive_payload() -> bytes:
    required = set(successor.required_candidate_source_paths())
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as archive:
        for relative in sorted(required):
            payload = (REPO_ROOT / relative).read_bytes()
            member = tarfile.TarInfo(relative)
            member.mode = 0o644
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
    return buffer.getvalue()


def _rewrite_json(binding: dict[str, object], value: object) -> None:
    payload = _json_bytes(value)
    Path(str(binding["path"])).write_bytes(payload)
    binding["sha256"] = _sha(payload)
    binding["byte_count"] = len(payload)


def _case(case_id: str) -> dict[str, object]:
    return next(case for case in oracle.expected_cases() if case["case_id"] == case_id)


def _block_payload(case: dict[str, object], rank: int, *, values: np.ndarray) -> bytes:
    nx1, nx2, nx3 = (int(value) for value in case["meshblock_nx"])
    bounds = tuple(tuple(float(value) for value in item) for item in case["bounds"])
    grid_x1, grid_x2, grid_x3 = (
        int(value) for value in case["meshblock_grid"]
    )
    if rank < 0 or rank >= grid_x1 * grid_x2 * grid_x3:
        raise AssertionError(f"rank {rank} is outside MeshBlock grid")
    logical_x1 = rank % grid_x1
    logical_x2 = (rank // grid_x1) % grid_x2
    logical_x3 = rank // (grid_x1 * grid_x2)
    logical = (logical_x1, logical_x2, logical_x3)
    geometry_values = []
    for axis in range(3):
        width = (bounds[axis][1] - bounds[axis][0]) / int(
            case["meshblock_grid"][axis]
        )
        geometry_values.extend(
            (
                bounds[axis][0] + logical[axis] * width,
                bounds[axis][0] + (logical[axis] + 1) * width,
            )
        )
    geometry = tuple(geometry_values)
    expected_shape = (nx3, nx2, nx1)
    if values.shape != expected_shape:
        raise AssertionError(f"expected shape {expected_shape}, got {values.shape}")
    index_and_logical = (
        0,
        nx1 - 1,
        0,
        nx2 - 1,
        0,
        nx3 - 1,
        logical_x1,
        logical_x2,
        logical_x3,
        0,
    )
    return (
        struct.pack("<10i", *index_and_logical)
        + struct.pack("<6d", *geometry)
        + np.asarray(values, dtype="<f4").tobytes()
    )


def _runtime_parameter_header(
    case: dict[str, object],
    *,
    field: str,
    cycle: int,
    common_overrides: dict[tuple[str, str], str] | None = None,
) -> bytes:
    parameters = copy.deepcopy(successor._expected_runtime_parameters(case))
    field_index = oracle.FIELDS.index(field) + 1
    for output_index in range(1, len(oracle.FIELDS) + 1):
        block = parameters[f"output{output_index}"]
        block["file_number"] = str(cycle if output_index <= field_index else cycle + 1)
        block["last_time"] = (
            "-1" if cycle == 0 and output_index <= field_index else "0"
        )
    for (block, key), value in (common_overrides or {}).items():
        if block not in parameters:
            parameters[block] = {}
        parameters[block][key] = value
    lines = []
    for block, values in parameters.items():
        lines.append(f"<{block}>")
        lines.extend(f"{key} = {value}" for key, value in values.items())
    return "\n".join(lines).encode("utf-8")


def _field_value(case: dict[str, object], field: str, cycle: int) -> float:
    if cycle == 0:
        return 0.0
    if field == "prtcl_rho":
        return oracle.EXPECTED_RHO
    basis = oracle._mode_basis(int(case["dimension"]))
    return oracle.EXPECTED_J_OVER_C * basis[oracle.FIELDS.index(field) - 1]


def _binary_payload(
    case: dict[str, object],
    *,
    field: str,
    cycle: int,
    rank: int,
    common_overrides: dict[tuple[str, str], str] | None = None,
) -> bytes:
    shape = tuple(reversed(tuple(int(value) for value in case["meshblock_nx"])))
    values = np.full(shape, _field_value(case, field, cycle), dtype=np.float64)
    parameter_header = _runtime_parameter_header(
        case, field=field, cycle=cycle, common_overrides=common_overrides
    )
    observed_time = "0.0" if cycle == 0 else "0.0025"
    return (
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        + f"  time={observed_time}\n".encode("ascii")
        + f"  cycle={cycle}\n".encode("ascii")
        + b"  size of location=8\n"
        b"  size of variable=4\n"
        b"  number of variables=1\n"
        + f"  variables:  {field}  \n".encode("ascii")
        + f"  header offset={len(parameter_header)}\n".encode("ascii")
        + parameter_header
        + _block_payload(case, rank, values=values)
    )


class Fixture:
    def __init__(self, root: Path):
        self.root = root
        self.authorized = root / "orion"
        self.authorized.mkdir()
        freeze = self.authorized / "clean_candidates" / str(
            uuid.uuid5(uuid.NAMESPACE_URL, "q043-test-freeze")
        )
        archive = _write(freeze / "source.tar", _source_archive_payload())
        executable = _write(freeze / "athena", b"\x7fELFq043-fixture\n", executable=True)
        environment = _write(
            freeze / "frontier_pic_environment.sh",
            b"#!/bin/bash\n# q043 fixture environment\n",
            executable=True,
        )
        git_commit = "1" * 40
        bundle = "2" * 64
        manifest_value = {
            "schema_version": 4,
            "freeze_id": str(uuid.uuid5(uuid.NAMESPACE_URL, "q043-test-freeze")),
            "created_utc": "2026-06-06T00:00:00Z",
            "prepared_artifacts": {},
            "source": {
                "archive_path": archive["path"],
                "archive_sha256": archive["sha256"],
                "source_bundle_sha256": bundle,
                "git_commit": git_commit,
                "worktree_status": "clean",
            },
            "build": {
                "source_archive_sha256": archive["sha256"],
                "source_bundle_sha256": bundle,
                "executable_path": executable["path"],
                "executable_sha256": executable["sha256"],
            },
        }
        manifest = _write(
            freeze / "clean_candidate_manifest.json", _json_bytes(manifest_value)
        )
        self.candidate = {
            "git_commit": git_commit,
            "source_bundle_sha256": bundle,
            "clean_candidate_manifest": manifest,
            "source_archive": archive,
            "executable": executable,
            "environment_profile": environment,
        }

    def inputs(
        self,
        case_id: str,
        *,
        common_overrides: dict[tuple[str, str], str] | None = None,
    ) -> tuple[dict[str, object], dict[str, object], list[dict[str, object]]]:
        case = _case(case_id)
        submission_id = str(uuid.uuid5(uuid.NAMESPACE_URL, f"submission:{case_id}"))
        reservation_id = str(uuid.uuid5(uuid.NAMESPACE_URL, f"reservation:{case_id}"))
        case_root = (
            self.authorized
            / successor.RUN_NAMESPACE
            / case_id
            / submission_id
        )
        raw_root = case_root / "raw"
        artifact_dir = case_root / "artifacts"
        raw_root.mkdir(parents=True)
        artifact_dir.mkdir()
        deck = successor._deck_binding(case, REPO_ROOT)
        normalized_candidate = successor._candidate_binding(
            self.candidate,
            self.authorized,
            source_root=REPO_ROOT,
            required_paths=successor.required_candidate_source_paths(),
        )
        launch_deck = _write(
            artifact_dir / "registration/input_deck.athinput",
            (REPO_ROOT / str(deck["path"])).read_bytes(),
        )
        launch_contract = successor._launch_contract(
            case=case,
            candidate=normalized_candidate,
            deck=deck,
            launch_deck=launch_deck,
            case_root=case_root,
            raw_root=raw_root,
            artifact_dir=artifact_dir,
        )
        launch_binding = _write(
            artifact_dir / "registration/launch_contract.json",
            _json_bytes(launch_contract),
        )
        ranks = int(case["mpi_ranks"])
        stdout_payload = (
            f"Q043_REGISTERED_EXECUTION case_id={case_id} mpi_world_size={ranks} "
            f"rank_ids={','.join(str(rank) for rank in range(ranks))}\n"
            "time=0.0025 cycle=1\n"
            "tlim=1 nlim=1\n"
            "Terminating on cycle limit\n"
            "Q043_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0\n"
        ).encode("utf-8")
        stdout_binding = _write(
            artifact_dir / "stdout/athena_stdout.txt", stdout_payload
        )
        job_id = str(int(hashlib.sha256(case_id.encode()).hexdigest()[:12], 16) + 1)
        terminal_value = {
            "schema_version": successor.SCHEMA_VERSION,
            "record_type": successor.TERMINAL_RECEIPT_RECORD_TYPE,
            "campaign_id": oracle.CAMPAIGN_ID,
            "case_id": case_id,
            "submission_id": submission_id,
            "slurm_job_id": job_id,
            "slurm_terminal_state": "COMPLETED",
            "slurm_exit_code": "0:0",
            "termination_reason": "cycle_limit",
            "terminal_cycle": 1,
            "observed_world_size": ranks,
            "stdout_sha256": stdout_binding["sha256"],
        }
        terminal_binding = _write(
            artifact_dir / "terminal/terminal_receipt.json",
            _json_bytes(terminal_value),
        )
        receipt_value = {
            "schema_version": successor.SCHEMA_VERSION,
            "record_type": successor.EXECUTION_RECEIPT_RECORD_TYPE,
            "receipt_role": "immutable_reconciled_registered_execution",
            "registration_scope": "registered_science",
            "reconciled": True,
            "campaign_id": oracle.CAMPAIGN_ID,
            "case_id": case_id,
            "reservation_id": reservation_id,
            "submission_id": submission_id,
            "reconciliation_event_sha256": hashlib.sha256(
                f"reconciliation:{case_id}".encode()
            ).hexdigest(),
            "source_commit": self.candidate["git_commit"],
            "source_bundle_sha256": self.candidate["source_bundle_sha256"],
            "source_archive_sha256": self.candidate["source_archive"]["sha256"],
            "clean_candidate_manifest_sha256": self.candidate[
                "clean_candidate_manifest"
            ]["sha256"],
            "executable_sha256": self.candidate["executable"]["sha256"],
            "environment_sha256": self.candidate["environment_profile"]["sha256"],
            "deck_sha256": deck["sha256"],
            "launch_contract_sha256": launch_binding["sha256"],
            "command": launch_contract["command"],
            "mpi_evidence": launch_contract["mpi_evidence"],
            "slurm_job_id": job_id,
            "slurm_terminal_state": "COMPLETED",
            "slurm_exit_code": "0:0",
            "raw_output_root": str(raw_root),
            "artifact_dir": str(artifact_dir),
            "stdout_sha256": stdout_binding["sha256"],
            "terminal_receipt_sha256": terminal_binding["sha256"],
            "pre_submit_manifest_sha256": hashlib.sha256(
                f"manifest:{case_id}".encode()
            ).hexdigest(),
        }
        receipt_binding = _write(
            artifact_dir / "registration/registered_execution_receipt.json",
            _json_bytes(receipt_value),
        )
        execution = {
            "case_root": str(case_root),
            "raw_output_root": str(raw_root),
            "artifact_dir": str(artifact_dir),
            "launch_deck": launch_deck,
            "launch_contract": launch_binding,
            "registered_execution_receipt": receipt_binding,
            "stdout_artifact": stdout_binding,
            "terminal_receipt": terminal_binding,
        }
        raw_artifacts = []
        for cycle in successor.REQUIRED_CYCLES:
            for field in successor.REQUIRED_FIELDS:
                for rank in range(ranks):
                    relative = successor._raw_relative_path(
                        case, field=field, cycle=cycle, rank=rank
                    )
                    payload = _binary_payload(
                        case,
                        field=field,
                        cycle=cycle,
                        rank=rank,
                        common_overrides=common_overrides,
                    )
                    binding = _write(raw_root / relative, payload)
                    raw_artifacts.append(
                        {
                            "path": relative,
                            "sha256": binding["sha256"],
                            "byte_count": binding["byte_count"],
                            "case_id": case_id,
                            "field": field,
                            "cycle": cycle,
                            "rank": rank,
                        }
                    )
        return case, execution, raw_artifacts

    def admission(
        self,
        case_id: str,
        *,
        common_overrides: dict[tuple[str, str], str] | None = None,
    ) -> dict[str, object]:
        _, execution, raw = self.inputs(case_id, common_overrides=common_overrides)
        return successor.build_case_admission(
            case_id=case_id,
            candidate_binding=self.candidate,
            execution_binding=execution,
            raw_artifacts=raw,
            source_root=REPO_ROOT,
            authorized_orion_root=self.authorized,
        )


class Q043RegisteredExecutionAdmissionTests(unittest.TestCase):
    def test_registered_dimension_valid_cases_rebuild_exactly_without_authority(
        self,
    ) -> None:
        selected = (
            "q043-current-oracle-d1-coarse-ppc1-single-cvr100",
            "q043-current-oracle-d1-coarse-ppc1-split_x1-cvr100",
            "q043-current-oracle-d2-coarse-ppc1-split_x2-cvr100",
            "q043-current-oracle-d2-coarse-ppc1-split_x1x2-cvr100",
            "q043-current-oracle-d3-coarse-ppc1-split_x3-cvr100",
            "q043-current-oracle-d3-coarse-ppc1-split_xyz-cvr100",
        )
        for case_id in selected:
            with self.subTest(case_id=case_id), tempfile.TemporaryDirectory() as temporary:
                fixture = Fixture(Path(temporary))
                record = fixture.admission(case_id)
                case = _case(case_id)
                self.assertEqual(
                    len(record["raw_artifacts"]),
                    2 * len(oracle.FIELDS) * int(case["mpi_ranks"]),
                )
                self.assertTrue(
                    record["hardened_raw_oracle_result"][
                        "registered_execution_raw_oracle_check_pass"
                    ]
                )
                self.assertFalse(
                    record["source_local_insufficiency"][
                        "source_local_or_synthetic_raw_analysis_sufficient"
                    ]
                )
                self.assertTrue(
                    all(
                        value is False
                        for key, value in record["authorization"].items()
                        if key.endswith("_authorized")
                    )
                )
                self.assertEqual(
                    successor.validate_case_admission(
                        record,
                        source_root=REPO_ROOT,
                        authorized_orion_root=fixture.authorized,
                    ),
                    record,
                )

    def test_relabeling_and_source_local_synthetic_substitution_are_insufficient(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        other = "q043-current-oracle-d1-coarse-ppc1-single-cvr1000"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            with self.assertRaises(successor.AdmissionError):
                successor.build_case_admission(
                    case_id=other,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )
            admitted = successor.build_case_admission(
                case_id=case_id,
                candidate_binding=fixture.candidate,
                execution_binding=execution,
                raw_artifacts=raw,
                source_root=REPO_ROOT,
                authorized_orion_root=fixture.authorized,
            )
            source_local = admitted["hardened_raw_oracle_result"]["source_local_oracle_result"]
            with self.assertRaisesRegex(
                successor.AdmissionError, "requires registered Q043 matrix"
            ):
                successor.validate_downstream_q023_q019_prerequisite(
                    source_local,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

    def test_source_local_result_schema_and_publication_authority_are_strict(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            admitted = fixture.admission(case_id)
            case = admitted["case_contract"]
            source_local = admitted["hardened_raw_oracle_result"][
                "source_local_oracle_result"
            ]
            self.assertEqual(
                successor._validate_source_local_case_result(source_local, case=case),
                source_local,
            )

            schema_drift = copy.deepcopy(source_local)
            schema_drift["unexpected_authority_surface"] = False
            with self.assertRaisesRegex(successor.AdmissionError, "keys drifted"):
                successor._validate_source_local_case_result(schema_drift, case=case)

            publication_drift = copy.deepcopy(source_local)
            publication_drift["publication_authorized"] = True
            with self.assertRaisesRegex(successor.AdmissionError, "authority drifted"):
                successor._validate_source_local_case_result(publication_drift, case=case)

    def test_split_without_exact_mpi_execution_evidence_fails(self) -> None:
        case_id = "q043-current-oracle-d3-coarse-ppc1-split_xyz-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            contract = json.loads(
                Path(str(execution["launch_contract"]["path"])).read_text()
            )
            contract["mpi_evidence"]["mpi_enabled"] = False
            _rewrite_json(execution["launch_contract"], contract)
            receipt = json.loads(
                Path(str(execution["registered_execution_receipt"]["path"])).read_text()
            )
            receipt["mpi_evidence"]["mpi_enabled"] = False
            receipt["launch_contract_sha256"] = execution["launch_contract"]["sha256"]
            _rewrite_json(execution["registered_execution_receipt"], receipt)
            with self.assertRaisesRegex(successor.AdmissionError, "launch contract"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

    def test_missing_duplicate_and_hardlinked_output_inventory_fails(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            with self.assertRaisesRegex(successor.AdmissionError, "exactly one cycle-zero"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw[:-1],
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )
            duplicate = copy.deepcopy(raw)
            duplicate.append(copy.deepcopy(duplicate[-1]))
            with self.assertRaisesRegex(successor.AdmissionError, "path reused"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=duplicate,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            Path(execution["raw_output_root"], "bin/undeclared.00002.bin").write_bytes(
                b"undeclared\n"
            )
            with self.assertRaisesRegex(successor.AdmissionError, "filesystem inventory"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            first = Path(execution["raw_output_root"], str(raw[0]["path"]))
            second = Path(execution["raw_output_root"], str(raw[1]["path"]))
            second.unlink()
            os.link(first, second)
            raw[1]["sha256"] = raw[0]["sha256"]
            raw[1]["byte_count"] = raw[0]["byte_count"]
            with self.assertRaisesRegex(successor.AdmissionError, "unaliased regular file"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

    def test_wrong_candidate_deck_and_executable_fail(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            bad_candidate = copy.deepcopy(fixture.candidate)
            bad_candidate["source_archive"]["sha256"] = "9" * 64
            with self.assertRaisesRegex(successor.AdmissionError, "SHA-256 drifted"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=bad_candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            bad_candidate = copy.deepcopy(fixture.candidate)
            malformed = b"not a source archive\n"
            Path(str(bad_candidate["source_archive"]["path"])).write_bytes(malformed)
            bad_candidate["source_archive"]["sha256"] = _sha(malformed)
            bad_candidate["source_archive"]["byte_count"] = len(malformed)
            manifest = json.loads(
                Path(str(bad_candidate["clean_candidate_manifest"]["path"])).read_text()
            )
            manifest["source"]["archive_sha256"] = bad_candidate["source_archive"]["sha256"]
            manifest["build"]["source_archive_sha256"] = bad_candidate["source_archive"][
                "sha256"
            ]
            _rewrite_json(bad_candidate["clean_candidate_manifest"], manifest)
            with self.assertRaisesRegex(successor.AdmissionError, "source archive"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=bad_candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            bad_candidate = copy.deepcopy(fixture.candidate)
            Path(str(bad_candidate["executable"]["path"])).write_bytes(b"wrong executable\n")
            with self.assertRaisesRegex(successor.AdmissionError, "SHA-256 drifted"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=bad_candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            launch_deck_path = Path(str(execution["launch_deck"]["path"]))
            bad_launch_deck = launch_deck_path.read_bytes() + b"# drift\n"
            launch_deck_path.write_bytes(bad_launch_deck)
            execution["launch_deck"]["sha256"] = _sha(bad_launch_deck)
            execution["launch_deck"]["byte_count"] = len(bad_launch_deck)
            with self.assertRaisesRegex(successor.AdmissionError, "immutable launch deck"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            contract = json.loads(
                Path(str(execution["launch_contract"]["path"])).read_text()
            )
            contract["deck"]["reviewed_checked_in_deck"]["sha256"] = "9" * 64
            _rewrite_json(execution["launch_contract"], contract)
            receipt = json.loads(
                Path(str(execution["registered_execution_receipt"]["path"])).read_text()
            )
            receipt["deck_sha256"] = "9" * 64
            receipt["launch_contract_sha256"] = execution["launch_contract"]["sha256"]
            _rewrite_json(execution["registered_execution_receipt"], receipt)
            with self.assertRaisesRegex(successor.AdmissionError, "launch contract"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

    def test_bad_terminal_state_and_stdout_fail(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            receipt = json.loads(
                Path(str(execution["registered_execution_receipt"]["path"])).read_text()
            )
            receipt["slurm_terminal_state"] = "FAILED"
            _rewrite_json(execution["registered_execution_receipt"], receipt)
            with self.assertRaisesRegex(successor.AdmissionError, "registered execution"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            stdout_path = Path(str(execution["stdout_artifact"]["path"]))
            bad_stdout = stdout_path.read_bytes().replace(
                b"Terminating on cycle limit", b"Terminating on wall clock limit"
            )
            stdout_path.write_bytes(bad_stdout)
            execution["stdout_artifact"]["sha256"] = _sha(bad_stdout)
            execution["stdout_artifact"]["byte_count"] = len(bad_stdout)
            terminal = json.loads(
                Path(str(execution["terminal_receipt"]["path"])).read_text()
            )
            terminal["stdout_sha256"] = execution["stdout_artifact"]["sha256"]
            _rewrite_json(execution["terminal_receipt"], terminal)
            receipt = json.loads(
                Path(str(execution["registered_execution_receipt"]["path"])).read_text()
            )
            receipt["stdout_sha256"] = execution["stdout_artifact"]["sha256"]
            receipt["terminal_receipt_sha256"] = execution["terminal_receipt"]["sha256"]
            _rewrite_json(execution["registered_execution_receipt"], receipt)
            with self.assertRaisesRegex(successor.AdmissionError, "stdout"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

    def test_path_escape_and_symlink_fail(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            escaped = copy.deepcopy(raw)
            escaped[0]["path"] = "../synthetic.bin"
            with self.assertRaisesRegex(successor.AdmissionError, "unsafe root-relative"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=escaped,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            path = Path(execution["raw_output_root"], str(raw[0]["path"]))
            retained = fixture.authorized / "retained.bin"
            retained.write_bytes(path.read_bytes())
            path.unlink()
            path.symlink_to(retained)
            with self.assertRaisesRegex(
                successor.AdmissionError, "escaped or is unavailable|symlink alias"
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

    def test_common_mode_runtime_metadata_drift_fails_registered_admission(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            case, execution, raw = fixture.inputs(
                case_id, common_overrides={("mhd", "gamma"): "9.0"}
            )
            cycle_one = {
                field: tuple(
                    Path(execution["raw_output_root"], str(item["path"]))
                    for item in raw
                    if item["cycle"] == 1 and item["field"] == field
                )
                for field in oracle.FIELDS
            }
            try:
                source_local = oracle.analyze_raw_case(case_id, cycle_one)
            except oracle.ContractError:
                source_local = None
            if source_local is not None:
                self.assertTrue(source_local["source_local_oracle_check_pass"])
            with self.assertRaisesRegex(successor.AdmissionError, "exact frozen deck"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )
            self.assertEqual(case["case_id"], case_id)

    def test_exact_132_case_dimension_valid_mpi_contract_is_required(
        self,
    ) -> None:
        records = [
            {
                "case_contract": copy.deepcopy(case),
                "execution_binding": {
                    "mpi_evidence": successor.expected_mpi_evidence(case)
                },
            }
            for case in oracle.expected_cases()
        ]
        x1_only = [
            record
            for record in records
            if record["case_contract"]["decomposition"] in ("single", "split_x1")
        ]
        with self.assertRaisesRegex(successor.AdmissionError, "exact 132-case"):
            successor.multidirectional_mpi_coverage_report(x1_only)
        report = successor.multidirectional_mpi_coverage_report(records)
        self.assertEqual(report["expected_case_count"], 132)
        self.assertEqual(report["registered_case_count"], 132)
        self.assertEqual(report["covered_partition_axes"], ["x1", "x2", "x3"])
        self.assertEqual(report["multi_axis_partition_count"], 24)
        self.assertEqual(
            report["registered_decompositions_by_dimension"],
            {
                "1": ["single", "split_x1"],
                "2": ["single", "split_x1", "split_x2", "split_x1x2"],
                "3": ["single", "split_x1", "split_x2", "split_x3", "split_xyz"],
            },
        )

        wrong_grid = copy.deepcopy(records)
        wrong_grid[-1]["execution_binding"]["mpi_evidence"]["meshblock_grid"] = [2, 2, 1]
        with self.assertRaisesRegex(successor.AdmissionError, "exact registered MPI"):
            successor.multidirectional_mpi_coverage_report(wrong_grid)

        noninteger_ppc = [copy.deepcopy(case) for case in oracle.expected_cases()]
        noninteger_ppc[0]["ppc"] = 1.5
        with self.assertRaisesRegex(successor.AdmissionError, "positive integer"):
            successor._validate_exact_case_matrix(noninteger_ppc)

    def test_incomplete_registered_matrix_fails_before_any_downstream_effect(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            admission = fixture.admission(case_id)
            with self.assertRaisesRegex(successor.AdmissionError, "matrix is incomplete"):
                successor.build_matrix_qualification(
                    case_admissions=[admission],
                    source_root=REPO_ROOT,
                    authorized_orion_root=fixture.authorized,
                )

    def test_readiness_binds_only_additive_non_authorizing_successor(self) -> None:
        record = json.loads(READINESS.read_text(encoding="utf-8"))
        for group in ("source_bindings", "foundation_bindings"):
            for binding in record[group].values():
                if isinstance(binding, dict):
                    path = REPO_ROOT / binding["path"]
                    self.assertEqual(_sha(path.read_bytes()), binding["sha256"])
        self.assertEqual(
            record["repaired_foundation_source_commit"],
            "5668ae65e8f02d9725da60b876221d3c441fdd95",
        )
        self.assertIn(
            "publication_authorized=false",
            record["registered_case_contract"]["raw_oracle_rule"],
        )
        self.assertTrue(record["source_local_insufficiency"]["explicit"])
        self.assertTrue(
            record["source_local_insufficiency"][
                "exact_132_case_dimension_valid_matrix_required"
            ]
        )
        self.assertTrue(record["mpi_coverage"]["multidirectional_required"])
        self.assertEqual(
            record["mpi_coverage"]["current_source_local_foundation_status"],
            "exact_132_case_dimension_valid_source_local_matrix_ready_nonqualifying",
        )
        self.assertEqual(record["mpi_coverage"]["exact_case_count"], 132)
        self.assertEqual(
            record["mpi_coverage"]["required_decompositions_by_dimension"],
            {
                "1": ["single", "split_x1"],
                "2": ["single", "split_x1", "split_x2", "split_x1x2"],
                "3": ["single", "split_x1", "split_x2", "split_x3", "split_xyz"],
            },
        )
        self.assertTrue(
            all(
                value is False
                for key, value in record["authorization"].items()
                if key.endswith("_authorized")
            )
        )
