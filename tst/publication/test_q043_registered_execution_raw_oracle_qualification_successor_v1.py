#!/usr/bin/env python3
"""Filesystem-backed adversarial tests for Q043 registered admission."""

from __future__ import annotations

import copy
from datetime import datetime, timedelta, timezone
import hashlib
import importlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import struct
import sys
import tarfile
import tempfile
from typing import Callable
import types
import unittest
from unittest.mock import patch
import uuid

import numpy as np

from tst.publication import (
    q043_bell_current_volume_aware_deposited_current_oracle as oracle,
)
from tst.publication import (
    q043_registered_execution_raw_oracle_qualification_successor_v1 as successor,
)
from tst.publication import (
    q043_registered_launch_policy_preparation_successor_v1 as preparation,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
CONTROL_PLANE_SOURCE = REPO_ROOT / "tst/publication/frontier_control_plane"
sys.path.insert(0, str(CONTROL_PLANE_SOURCE))
try:
    launch_trampoline = importlib.import_module("launch_trampoline")
    q043_producer = importlib.import_module("reconcile_q043_registered_execution")
finally:
    sys.path.pop(0)
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


def _binding(path: Path) -> dict[str, object]:
    payload = path.read_bytes()
    return {"path": str(path), "sha256": _sha(payload), "byte_count": len(payload)}


def _seal(path: Path, *, executable: bool = False) -> None:
    path.chmod(0o555 if executable else 0o444)


def _overwrite(path: Path, payload: bytes) -> None:
    path.chmod(0o644)
    path.write_bytes(payload)


REAL_TRUSTED_CANDIDATE_REVALIDATION = successor._trusted_clean_candidate_revalidation
REAL_TRUSTED_MIRRORED_LEDGER_STATE = successor._trusted_mirrored_ledger_state
REAL_TRUSTED_PRODUCER_REDERIVATION = successor._trusted_producer_rederivation


def _source_archive_payload(
    overrides: dict[str, bytes] | None = None,
) -> bytes:
    required = set(successor.required_candidate_source_paths())
    replacements = overrides or {}
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as archive:
        for relative in sorted(required):
            payload = replacements.get(relative, (REPO_ROOT / relative).read_bytes())
            member = tarfile.TarInfo(relative)
            member.mode = 0o644
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
    return buffer.getvalue()


def _rewrite_json(binding: dict[str, object], value: object) -> None:
    payload = _json_bytes(value)
    path = Path(str(binding["path"]))
    path.chmod(0o644)
    path.write_bytes(payload)
    _seal(path)
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
    value_override: float | None = None,
) -> bytes:
    shape = tuple(reversed(tuple(int(value) for value in case["meshblock_nx"])))
    value = _field_value(case, field, cycle) if value_override is None else value_override
    values = np.full(shape, value, dtype=np.float64)
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
        self.project_home = root / "project_home"
        self.project_home.mkdir()
        self.ledger_records: list[dict[str, object]] = []
        self.mirror_receipts: list[dict[str, object]] = []
        self.producer_sha256 = _sha(
            (
                REPO_ROOT
                / "tst/publication/frontier_control_plane/"
                "reconcile_q043_registered_execution.py"
            ).read_bytes()
        )
        self.trampoline_sha256 = _sha(
            (
                REPO_ROOT
                / "tst/publication/frontier_control_plane/launch_trampoline.py"
            ).read_bytes()
        )
        successor.AUTHORIZED_ORION_ROOT = self.authorized
        successor.AUTHORIZED_PROJECT_HOME_ROOT = self.project_home
        successor.AUTHORIZED_PROJECT_HOME_LEDGER_LEXICAL_ROOT = self.project_home
        successor._trusted_clean_candidate_revalidation = self._candidate_revalidation
        successor._trusted_mirrored_ledger_state = self._ledger_state
        successor._trusted_producer_rederivation = self._producer_rederivation
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
        for path, executable_mode in (
            (Path(str(archive["path"])), False),
            (Path(str(executable["path"])), True),
            (Path(str(environment["path"])), True),
            (Path(str(manifest["path"])), False),
        ):
            _seal(path, executable=executable_mode)
        self.candidate = {
            "git_commit": git_commit,
            "source_bundle_sha256": bundle,
            "clean_candidate_manifest": manifest,
            "source_archive": archive,
            "executable": executable,
            "environment_profile": environment,
        }

    def _ledger_state(self, version: str) -> dict[str, object]:
        if version != "a" * 64:
            raise successor.AdmissionError("test installed producer generation drifted")
        return {
            "records": copy.deepcopy(self.ledger_records),
            "mirror_receipts": copy.deepcopy(self.mirror_receipts),
            "installed_control_plane": {
                "version": version,
                "orion_path": str(self.authorized / "control_plane" / version),
                "project_home_path": str(self.project_home / "control_plane" / version),
                "producer_entrypoint": successor.Q043_PRODUCER_ENTRYPOINT,
                "producer_entrypoint_sha256": self.producer_sha256,
                "launch_trampoline_sha256": self.trampoline_sha256,
                "project_home_canonical_root": str(self.project_home),
                "project_home_ledger_root": str(self.project_home),
            },
        }

    def _candidate_revalidation(
        self, manifest_path: Path, *, manifest_sha256: str, git_commit: str
    ) -> dict[str, object]:
        if (
            manifest_path != Path(str(self.candidate["clean_candidate_manifest"]["path"]))
            or manifest_sha256 != self.candidate["clean_candidate_manifest"]["sha256"]
            or git_commit != self.candidate["git_commit"]
        ):
            raise successor.AdmissionError("test candidate revalidation binding drifted")
        return {
            "record_type": "frontier_pic_clean_candidate_read_only_revalidation",
            "schema_version": 1,
            "status": "passed",
            "clean_candidate_manifest": {
                "expected_sha256": manifest_sha256,
                "path": str(manifest_path),
                "sha256": manifest_sha256,
            },
            "source": {
                "git_commit": git_commit,
                "git_tree": "3" * 40,
                "source_bundle_sha256": self.candidate["source_bundle_sha256"],
            },
            "build": {
                "executable_sha256": self.candidate["executable"]["sha256"],
                "profile_id": "q043-test-profile",
                "receipt_control_plane_version": "a" * 64,
            },
            "validated_submodules": [],
        }

    def _producer_rederivation(self, **arguments: object) -> dict[str, object]:
        version = str(arguments["producer_version"])
        try:
            expected = q043_producer.derive_q043_registered_execution_evidence(
                copy.deepcopy(arguments["event"]),
                copy.deepcopy(arguments["mirror_ack"]),
                self._producer_inventory(version),
                authorized_pic_root=self.authorized,
                authorized_project_home_root=self.project_home,
            )
        except (OSError, TypeError, ValueError) as error:
            raise successor.AdmissionError(
                "installed Q043 producer could not re-derive execution evidence: "
                f"{error}"
            ) from error
        observed = (
            arguments["terminal_path"],
            arguments["terminal_payload"],
            arguments["receipt_path"],
            arguments["receipt_payload"],
            arguments["terminal_mirror_path"],
            arguments["receipt_mirror_path"],
        )
        if expected != observed:
            raise successor.AdmissionError(
                "installed Q043 producer re-derived different evidence bytes or paths"
            )
        return {
            "control_plane_version": version,
            "entrypoint": successor.Q043_PRODUCER_ENTRYPOINT,
            "entrypoint_sha256": self.producer_sha256,
            "launch_trampoline_sha256": self.trampoline_sha256,
            "receipt_sha256": _sha(arguments["receipt_payload"]),
            "terminal_receipt_sha256": _sha(arguments["terminal_payload"]),
            "exact_byte_rederivation_passed": True,
        }

    def _producer_inventory(self, version: str = "a" * 64) -> dict[str, object]:
        return {
            "schema_version": 1,
            "version": version,
            "files": [
                {
                    "path": successor.Q043_PRODUCER_ENTRYPOINT,
                    "sha256": self.producer_sha256,
                },
                {
                    "path": "launch_trampoline.py",
                    "sha256": self.trampoline_sha256,
                },
            ],
        }

    def _snapshot(
        self,
        destination: Path,
        *,
        role: str,
        source: Path,
        executable: bool = False,
    ) -> dict[str, object]:
        payload = source.read_bytes()
        binding = _write(destination, payload, executable=executable)
        _seal(destination, executable=executable)
        return {
            "path": binding["path"],
            "role": role,
            "sha256": binding["sha256"],
            "source_path": str(source),
            "source_sha256": binding["sha256"],
        }

    def _freeze_run(self, artifact_dir: Path) -> None:
        analysis = artifact_dir / "analysis"
        analysis.mkdir(mode=0o700)
        records = []
        for path in sorted(artifact_dir.rglob("*")):
            if path == analysis or analysis in path.parents or path.is_dir():
                continue
            relative = path.relative_to(artifact_dir).as_posix()
            payload = path.read_bytes()
            path.chmod(0o444)
            records.append(
                {"path": relative, "sha256": _sha(payload), "size": len(payload)}
            )
        inventory = artifact_dir / "artifact_inventory.json"
        inventory.write_bytes(_json_bytes({"schema_version": 1, "files": records}))
        inventory.chmod(0o444)
        for path in sorted(
            (item for item in artifact_dir.rglob("*") if item.is_dir() and item != analysis),
            key=lambda item: len(item.parts),
            reverse=True,
        ):
            path.chmod(0o555)
        artifact_dir.chmod(0o555)

    def inputs(
        self,
        case_id: str,
        *,
        common_overrides: dict[tuple[str, str], str] | None = None,
        raw_value_overrides: dict[tuple[int, str, int], float] | None = None,
        after_completion: Callable[["Fixture"], None] | None = None,
    ) -> tuple[dict[str, object], dict[str, object], list[dict[str, object]]]:
        case = _case(case_id)
        submission_id = str(uuid.uuid5(uuid.NAMESPACE_URL, f"submission:{case_id}"))
        reservation_id = str(uuid.uuid5(uuid.NAMESPACE_URL, f"reservation:{case_id}"))
        artifact_dir = (
            self.authorized
            / "runs"
            / successor.REGISTERED_CAMPAIGN
            / submission_id
        )
        raw_root = artifact_dir / "raw"
        raw_root.mkdir(parents=True)
        deck = successor._deck_binding(case, REPO_ROOT)
        ranks = int(case["mpi_ranks"])
        launch_contract = {
            "schema_version": 1,
            "executor": "trusted_trampoline_athena_argv_v1",
            "pre_actions": [],
            "actions": [
                {
                    "action_id": "q043-current-oracle",
                    "kind": "athena",
                    "resources": {
                        "nodes": 1,
                        "tasks": ranks,
                        "cpus_per_task": 1,
                        "gpus_per_task": 1,
                        "gpu_bind": "closest",
                    },
                    "arguments": [
                        {"literal": "-i"},
                        {"snapshot_role": "input-deck"},
                        {"literal": "-d"},
                        {"artifact_directory": "raw"},
                        {"literal": "time/nlim=1"},
                    ],
                    "stdout_artifact": "athena_stdout.txt",
                    "stderr_artifact": "athena_stderr.txt",
                }
            ],
            "post_actions": [],
        }
        stdout_payload = (
            "".join(
                "PIC trusted GPU launch: "
                f"rank={rank} host=frontier-test "
                f"ROCR_VISIBLE_DEVICES={rank} "
                "linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa\n"
                for rank in range(ranks)
            )
            + "time=0.0025 cycle=1\n"
            + "tlim=1 nlim=1\n"
            + "Terminating on cycle limit\n"
            + f"Q043_REGISTERED_EXECUTION case_id={case_id} "
            + f"mpi_world_size={ranks} "
            + f"rank_ids={','.join(str(rank) for rank in range(ranks))}\n"
            + "Q043_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0\n"
        ).encode("utf-8")
        _write(artifact_dir / "athena_stdout.txt", stdout_payload)
        _write(artifact_dir / "athena_stderr.txt", b"")
        job_id = str(int(hashlib.sha256(case_id.encode()).hexdigest()[:12], 16) + 1)
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
                        value_override=(raw_value_overrides or {}).get(
                            (cycle, field, rank)
                        ),
                    )
                    binding = _write(raw_root / relative, payload)
                    _seal(Path(str(binding["path"])))
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
        self._freeze_run(artifact_dir)
        manifest_root = (
            self.authorized
            / "manifests"
            / successor.REGISTERED_CAMPAIGN
            / submission_id
        )
        snapshot_root = manifest_root / "snapshot"
        job_script_source = _write(
            self.authorized / "test_sources" / "q043_job.sh",
            b"#!/bin/bash\n# trusted q043 fixture job\n",
            executable=True,
        )
        _seal(Path(str(job_script_source["path"])), executable=True)
        snapshot_files = [
            self._snapshot(
                snapshot_root / "clean_candidate_manifest.json",
                role="clean-candidate-manifest",
                source=Path(str(self.candidate["clean_candidate_manifest"]["path"])),
            ),
            self._snapshot(
                snapshot_root / "athena",
                role="executable",
                source=Path(str(self.candidate["executable"]["path"])),
                executable=True,
            ),
            self._snapshot(
                snapshot_root / f"{case_id}.athinput",
                role="input-deck",
                source=Path(str(deck["absolute_path"])),
            ),
            self._snapshot(
                snapshot_root / "frontier_pic_environment.sh",
                role="environment-profile",
                source=Path(str(self.candidate["environment_profile"]["path"])),
                executable=True,
            ),
            self._snapshot(
                snapshot_root / "q043_job.sh",
                role="job-script",
                source=Path(str(job_script_source["path"])),
                executable=True,
            ),
        ]
        control_plane_version = "a" * 64
        authorization_id = "q043-test-registered-raw-oracle"
        manifest_value = {
            "schema_version": 1,
            "control_plane_version": control_plane_version,
            "submission_id": submission_id,
            "pic_root": str(self.authorized),
            "campaign": successor.REGISTERED_CAMPAIGN,
            "test_id": case_id,
            "submission_scope": "registered_science",
            "registered_science_authorization_id": authorization_id,
            "git_commit": self.candidate["git_commit"],
            "launch_contract": launch_contract,
            "artifact_dir": str(artifact_dir),
            "clean_candidate_manifest_path": self.candidate[
                "clean_candidate_manifest"
            ]["path"],
            "clean_candidate_manifest_sha256": self.candidate[
                "clean_candidate_manifest"
            ]["sha256"],
            "snapshot_files": snapshot_files,
        }
        manifest_binding = _write(
            manifest_root / "pre_submit_manifest.json", _json_bytes(manifest_value)
        )
        _seal(Path(str(manifest_binding["path"])))
        artifact_dir_fd = os.open(artifact_dir, os.O_RDONLY | os.O_DIRECTORY)
        try:
            completion_publication = (
                launch_trampoline._publish_trampoline_completion_receipt_at(
                    artifact_dir_fd,
                    artifact_dir,
                    manifest_value,
                    append_anchor=lambda binding: {
                        "trampoline_completion": binding
                    },
                    manifest_path=Path(str(manifest_binding["path"])),
                    manifest_sha256=str(manifest_binding["sha256"]),
                    reservation_id=reservation_id,
                    submission_id=submission_id,
                    slurm_job_id=job_id,
                    authorized_pic_root=self.authorized,
                    authorized_project_home_root=self.project_home,
                )
            )
        finally:
            os.close(artifact_dir_fd)
        common_event = {
            "reservation_id": reservation_id,
            "submission_id": submission_id,
            "job_id": job_id,
            "control_plane_version": control_plane_version,
            "campaign": successor.REGISTERED_CAMPAIGN,
            "test_id": case_id,
            "submission_scope": "registered_science",
            "registered_science_authorization_id": authorization_id,
            "clean_candidate_manifest_sha256": self.candidate[
                "clean_candidate_manifest"
            ]["sha256"],
            "manifest_path": manifest_binding["path"],
            "manifest_sha256": manifest_binding["sha256"],
            "git_commit": self.candidate["git_commit"],
            "executable_sha256": self.candidate["executable"]["sha256"],
            "artifact_dir": str(artifact_dir),
            "trampoline_completion": completion_publication["ledger_binding"],
        }
        completion_event = {
            **common_event,
            "sequence_number": len(self.ledger_records) + 1,
            "event_type": "trampoline_completion",
            "state": "submitted",
            "reconciled": False,
        }
        completion_event["event_sha256"] = _sha(_json_bytes(completion_event))
        self.ledger_records.append(completion_event)
        completion_mirror_ack = {
            "mirrored_event_sha256": completion_event["event_sha256"],
            "mirror_destination": str(self.project_home / "ledger/node_hours.jsonl"),
            "mirror_transport": "filesystem_copy",
            "mirror_acknowledged_utc": "2026-06-06T00:00:00Z",
            "mirror_ack_sha256": "c" * 64,
        }
        self.mirror_receipts.append(completion_mirror_ack)
        self.last_completion_event = copy.deepcopy(completion_event)
        self.last_completion_mirror_ack = copy.deepcopy(completion_mirror_ack)
        if after_completion is not None:
            after_completion(self)
        event = {
            **common_event,
            "sequence_number": len(self.ledger_records) + 1,
            "event_type": "reconciliation",
            "state": "COMPLETED",
            "scheduler_exit_code": "0:0",
            "reconciled": True,
            "reconciled_by_control_plane_version": control_plane_version,
        }
        event["event_sha256"] = _sha(_json_bytes(event))
        self.ledger_records.append(event)
        mirror_ack = {
            "mirrored_event_sha256": event["event_sha256"],
            "mirror_destination": str(self.project_home / "ledger/node_hours.jsonl"),
            "mirror_transport": "filesystem_copy",
            "mirror_acknowledged_utc": "2026-06-06T00:00:00Z",
            "mirror_ack_sha256": "b" * 64,
        }
        self.mirror_receipts.append(mirror_ack)
        evidence = q043_producer.publish_q043_registered_execution_evidence(
            event,
            mirror_ack,
            self._producer_inventory(control_plane_version),
            authorized_pic_root=self.authorized,
            authorized_project_home_root=self.project_home,
        )
        produced = {"reconciliation": event, "evidence": evidence}
        self.last_event = copy.deepcopy(event)
        self.last_mirror_ack = copy.deepcopy(mirror_ack)
        self.last_produced = copy.deepcopy(produced)
        terminal_binding = _binding(
            Path(produced["evidence"]["terminal_receipt"]["orion_path"])
        )
        receipt_binding = _binding(
            Path(produced["evidence"]["registered_execution_receipt"]["orion_path"])
        )
        execution = {
            "artifact_dir": str(artifact_dir),
            "registered_execution_receipt": receipt_binding,
            "terminal_receipt": terminal_binding,
        }
        return case, execution, raw_artifacts

    def replace_completion_receipt_directories(self, *, malicious: bool = True) -> None:
        submission_id = str(self.last_completion_event["submission_id"])
        replacement_payload: bytes | None = None
        for root in (self.authorized, self.project_home):
            directory = (
                root
                / successor.TRAMPOLINE_COMPLETION_NAMESPACE
                / submission_id
            )
            detached = directory.with_name(f"{directory.name}.original")
            directory.rename(detached)
            shutil.copytree(detached, directory, copy_function=shutil.copy2)
            receipt = directory / successor.TRAMPOLINE_COMPLETION_NAME
            if malicious:
                if replacement_payload is None:
                    value = json.loads(receipt.read_text(encoding="utf-8"))
                    value["artifact_root_identity"]["inode"] += 1
                    replacement_payload = _json_bytes(value)
                receipt.chmod(0o644)
                receipt.write_bytes(replacement_payload)
                receipt.chmod(0o444)
            directory.chmod(0o500)

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
        )

    def reseal_receipt(
        self, execution: dict[str, object], receipt: dict[str, object]
    ) -> None:
        _rewrite_json(execution["registered_execution_receipt"], receipt)
        mirror_path = Path(
            receipt["project_home_mirrors"]["registered_execution_receipt_path"]
        )
        _overwrite(mirror_path, _json_bytes(receipt))
        _seal(mirror_path)

    def replace_raw_artifact(
        self,
        case: dict[str, object],
        execution: dict[str, object],
        artifact: dict[str, object],
        *,
        value: float,
    ) -> None:
        payload = _binary_payload(
            case,
            field=str(artifact["field"]),
            cycle=int(artifact["cycle"]),
            rank=int(artifact["rank"]),
            value_override=value,
        )
        path = Path(execution["artifact_dir"], "raw", str(artifact["path"]))
        _overwrite(path, payload)
        _seal(path)
        artifact["sha256"] = _sha(payload)
        artifact["byte_count"] = len(payload)


class Q043RegisteredExecutionAdmissionTests(unittest.TestCase):
    def test_captured_installed_generation_ignores_swap_after_verification(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            orion = root / "orion"
            project_home = root / "project_home"
            sources = {
                "captured_dependency.py": b'VALUE = "captured-original" \n',
                "captured_entrypoint.py": (
                    b"import captured_dependency\n"
                    b"RESULT = captured_dependency.VALUE\n"
                ),
            }
            records = [
                {"path": name, "sha256": _sha(payload)}
                for name, payload in sorted(sources.items())
            ]
            version = _sha(
                json.dumps(
                    records, sort_keys=True, separators=(",", ":")
                ).encode("utf-8")
            )
            for trusted_root in (orion, project_home):
                generation = trusted_root / "control_plane" / version
                generation.mkdir(parents=True)
                for name, payload in sources.items():
                    _write(generation / name, payload)
                    _seal(generation / name)
                _write(
                    generation / "inventory.json",
                    _json_bytes(
                        {"schema_version": 1, "version": version, "files": records}
                    ),
                )
                _seal(generation / "inventory.json")
                generation.chmod(0o555)
            original_roots = (
                successor.AUTHORIZED_ORION_ROOT,
                successor.AUTHORIZED_PROJECT_HOME_ROOT,
                successor.AUTHORIZED_PROJECT_HOME_LEDGER_LEXICAL_ROOT,
            )
            real_pair = successor._installed_control_plane_pair
            swapped = False

            def pair_then_swap(observed_version: str) -> dict[str, object]:
                nonlocal swapped
                pair = real_pair(observed_version)
                for trusted_root in (orion, project_home):
                    generation = trusted_root / "control_plane" / version
                    generation.chmod(0o755)
                    dependency = generation / "captured_dependency.py"
                    dependency.chmod(0o644)
                    dependency.write_bytes(b'VALUE = "malicious-replacement" \n')
                    dependency.chmod(0o444)
                    generation.chmod(0o555)
                swapped = True
                return pair

            try:
                successor.AUTHORIZED_ORION_ROOT = orion
                successor.AUTHORIZED_PROJECT_HOME_ROOT = project_home
                successor.AUTHORIZED_PROJECT_HOME_LEDGER_LEXICAL_ROOT = project_home
                with patch.object(
                    successor,
                    "_installed_control_plane_pair",
                    side_effect=pair_then_swap,
                ):
                    with successor._installed_control_plane_modules(
                        version, ("captured_entrypoint.py",)
                    ) as (modules, _):
                        self.assertTrue(swapped)
                        self.assertEqual(
                            modules["captured_entrypoint.py"].RESULT,
                            "captured-original",
                        )
            finally:
                (
                    successor.AUTHORIZED_ORION_ROOT,
                    successor.AUTHORIZED_PROJECT_HOME_ROOT,
                    successor.AUTHORIZED_PROJECT_HOME_LEDGER_LEXICAL_ROOT,
                ) = original_roots

    def test_authoritative_scientific_parser_ignores_ambient_import_poisoning(
        self,
    ) -> None:
        canonical_name = "tst.publication.analyze_q011_section54_outputs"
        poisoned = types.ModuleType(canonical_name)
        poisoned.parse_athenak_binary_bytes = lambda *_args, **_kwargs: "poisoned"
        previous = sys.modules.get(canonical_name)
        sys.modules[canonical_name] = poisoned
        try:
            exact = oracle._load_exact_binary_parser()
            self.assertIsNot(exact, poisoned)
            self.assertTrue(callable(exact.parse_athenak_binary_bytes))
            self.assertIs(successor.binary, oracle.binary)
            self.assertEqual(
                exact._Q043_CAPTURED_SOURCE_SHA256,
                oracle.BINARY_PARSER_SOURCE_SHA256,
            )
            self.assertIn(
                "tst/publication/analyze_q011_section54_outputs.py",
                successor.required_candidate_source_paths(),
            )
        finally:
            if previous is None:
                sys.modules.pop(canonical_name, None)
            else:
                sys.modules[canonical_name] = previous

    def test_all_132_cases_traverse_real_controller_path_to_admission(self) -> None:
        control_plane_path = str(CONTROL_PLANE_SOURCE)
        sys.path.insert(0, control_plane_path)
        try:
            import reconcile_frontier_job
            from test_control_plane import SnapshotTests
        finally:
            sys.path.pop(0)

        harness = SnapshotTests()
        original_roots = (
            successor.AUTHORIZED_ORION_ROOT,
            successor.AUTHORIZED_PROJECT_HOME_ROOT,
            successor.AUTHORIZED_PROJECT_HOME_LEDGER_LEXICAL_ROOT,
        )
        original_bootstrap = getattr(sys, "_pic_control_plane_bootstrapped", None)
        try:
            harness.setUp()
            cases = oracle.expected_cases()
            first_case = cases[0]
            first_case_id = str(first_case["case_id"])
            deck_path = Path(
                successor._deck_binding(first_case, REPO_ROOT)["absolute_path"]
            )
            (harness.sources / "input.athinput").write_bytes(deck_path.read_bytes())
            (harness.sources / "environment.sh").chmod(0o755)
            with patch.object(
                harness,
                "_launch_contract",
                return_value=preparation._launch_contract(first_case),
            ):
                candidate_path = harness._write_science_config(
                    authorize=False,
                    campaign=successor.REGISTERED_CAMPAIGN,
                    test_id=first_case_id,
                    registered_science_authorization_id="q043-controller-e2e",
                    input_deck=str(deck_path),
                    evidence_class="q043_registered_execution_test",
                    physical_mode="paper_mhd_pic_vl2_tsc",
                    artifact_dir=str(
                        harness.pic_root
                        / "runs"
                        / successor.REGISTERED_CAMPAIGN
                        / harness.submission_id
                    ),
                )
                candidate_root = candidate_path.parent
                candidate_manifest = json.loads(candidate_path.read_text(encoding="utf-8"))
                source_archive = candidate_root / "source.tar"
                executable = candidate_root / "athena"
                candidate_root.chmod(0o755)
                source_archive.chmod(0o644)
                source_archive.write_bytes(_source_archive_payload())
                source_archive.chmod(0o444)
                executable.chmod(0o755)
                executable.write_bytes(Path("/bin/true").read_bytes())
                executable.chmod(0o555)
                candidate_manifest["source"]["archive_sha256"] = _sha(
                    source_archive.read_bytes()
                )
                candidate_manifest["build"]["source_archive_sha256"] = _sha(
                    source_archive.read_bytes()
                )
                candidate_manifest["build"]["executable_sha256"] = _sha(
                    executable.read_bytes()
                )
                candidate_path.chmod(0o644)
                candidate_path.write_bytes(_json_bytes(candidate_manifest))
                candidate_path.chmod(0o444)
                candidate_root.chmod(0o555)
                harness.authorized_clean_candidate_source_root = None
                environment_profile = (
                    harness.control_plane_dir / "frontier_pic_environment.sh"
                )
                job_script = harness.control_plane_dir / "frontier_job.sh"
                reconciler = (
                    harness.control_plane_dir
                    / "reconcile_q043_registered_execution.py"
                )
                final_bindings = {
                    "record_type": preparation.FINAL_BINDING_RECORD_TYPE,
                    "schema_version": 1,
                    "source_commit": candidate_manifest["source"]["git_commit"],
                    "source_bundle_sha256": candidate_manifest["source"][
                        "source_bundle_sha256"
                    ],
                    "source_archive_path": str(source_archive),
                    "source_archive_sha256": _sha(source_archive.read_bytes()),
                    "clean_candidate_manifest_path": str(candidate_path),
                    "clean_candidate_manifest_sha256": _sha(
                        candidate_path.read_bytes()
                    ),
                    "executable_path": str(executable),
                    "executable_sha256": _sha(executable.read_bytes()),
                    "installed_control_plane_version": harness.control_plane_version,
                    "orion_installed_control_plane_root": str(
                        harness.control_plane_dir
                    ),
                    "project_home_installed_control_plane_root": str(
                        harness.project_home_control_plane_dir
                    ),
                    "environment_profile_path": str(environment_profile),
                    "environment_profile_sha256": _sha(
                        environment_profile.read_bytes()
                    ),
                    "job_script_path": str(job_script),
                    "job_script_sha256": _sha(job_script.read_bytes()),
                    "analysis_script_paths": [str(reconciler)],
                    "analysis_script_sha256": [_sha(reconciler.read_bytes())],
                    "reconcile_q043_registered_execution_path": str(reconciler),
                    "reconcile_q043_registered_execution_sha256": _sha(
                        reconciler.read_bytes()
                    ),
                }
                harness._write_policy(
                    registered_science_slices=[],
                    admission_smoke_overrides={"status": "closed_after_pass"},
                )
                harness._promote_policy()
                baseline_policy = json.loads(
                    harness.policy.read_text(encoding="utf-8")
                )
                now = datetime.now(timezone.utc).replace(microsecond=0)
                timeout_path = harness.sources / "q043-timeout.json"
                queue_path = harness.sources / "q043-empty-queue.txt"
                queue_path.write_bytes(b"")
                queue_path.chmod(0o444)
                with (
                    patch.object(
                        preparation, "AUTHORIZED_ORION_ROOT", harness.pic_root
                    ),
                    patch.object(
                        preparation,
                        "CANONICAL_PROJECT_HOME_ROOT",
                        harness.project_home_root,
                    ),
                ):
                    timeout_path.write_bytes(
                        _json_bytes(
                            preparation.materialize_q043_timeout_margin(
                                final_bindings=final_bindings,
                                measured_utc=harness._utc(now - timedelta(minutes=1)),
                                expires_utc=harness._utc(now + timedelta(hours=1)),
                            )
                        )
                    )
                    timeout_path.chmod(0o444)
                    policy = preparation.materialize_q043_promotable_policy(
                        baseline_policy=baseline_policy,
                        final_bindings=final_bindings,
                    )
                self.assertEqual(len(policy["registered_science_slices"]), 132)
                harness.policy.write_text(json.dumps(policy), encoding="utf-8")
                with patch(
                    "promote_active_policy.revalidate_clean_candidate",
                    return_value={
                        "status": "passed",
                        "current_control_plane_version": harness.control_plane_version,
                        "build": {
                            "receipt_control_plane_version": (
                                harness.control_plane_version
                            )
                        },
                    },
                ):
                    harness._promote_policy(patch_clean_candidate_revalidation=False)
            successor.AUTHORIZED_ORION_ROOT = harness.pic_root
            successor.AUTHORIZED_PROJECT_HOME_ROOT = harness.project_home_root
            successor.AUTHORIZED_PROJECT_HOME_LEDGER_LEXICAL_ROOT = (
                harness.project_home_root
            )
            candidate_manifest = json.loads(candidate_path.read_text(encoding="utf-8"))
            candidate_binding = {
                "git_commit": candidate_manifest["source"]["git_commit"],
                "source_bundle_sha256": candidate_manifest["source"][
                    "source_bundle_sha256"
                ],
                "clean_candidate_manifest": _binding(candidate_path),
                "source_archive": _binding(
                    Path(candidate_manifest["source"]["archive_path"])
                ),
                "executable": _binding(
                    Path(candidate_manifest["build"]["executable_path"])
                ),
                "environment_profile": _binding(environment_profile),
            }

            def revalidation_fixture(
                manifest: Path, *, manifest_sha256: str, git_commit: str
            ) -> dict[str, object]:
                self.assertEqual(manifest, candidate_path)
                self.assertEqual(manifest_sha256, _sha(candidate_path.read_bytes()))
                self.assertEqual(git_commit, candidate_manifest["source"]["git_commit"])
                return {
                    "record_type": "frontier_pic_clean_candidate_read_only_revalidation",
                    "schema_version": 1,
                    "status": "passed",
                    "clean_candidate_manifest": {
                        "expected_sha256": manifest_sha256,
                        "path": str(manifest),
                        "sha256": manifest_sha256,
                    },
                    "source": {
                        "git_commit": git_commit,
                        "source_bundle_sha256": candidate_manifest["source"][
                            "source_bundle_sha256"
                        ],
                    },
                    "build": {
                        "executable_sha256": candidate_manifest["build"][
                            "executable_sha256"
                        ],
                    },
                }

            authorizations = {
                str(item["test_id"]): str(item["authorization_id"])
                for item in policy["registered_science_slices"]
            }
            reserve_globals = harness._reserve.__func__.__globals__["reserve"].__globals__
            admissions = []
            with successor._installed_control_plane_modules(
                harness.control_plane_version,
                (successor.Q043_PRODUCER_ENTRYPOINT,),
            ) as (modules, _pair):
                installed_producer = modules[successor.Q043_PRODUCER_ENTRYPOINT]
                for index, case in enumerate(cases):
                    case_id = str(case["case_id"])
                    harness.submission_id = str(uuid.uuid4())
                    with (
                        self.subTest(case_id=case_id),
                        patch.object(
                            preparation, "AUTHORIZED_ORION_ROOT", harness.pic_root
                        ),
                        patch.object(
                            preparation,
                            "CANONICAL_PROJECT_HOME_ROOT",
                            harness.project_home_root,
                        ),
                    ):
                        config = preparation.materialize_q043_pre_submit_config(
                            case_id=case_id,
                            submission_id=harness.submission_id,
                            final_bindings=final_bindings,
                            pre_manifest_attestation=harness._sealed_operator_attestation(
                                authorizations[case_id], "pre_manifest"
                            ),
                            timeout_margin_artifact=timeout_path,
                            queue_snapshot=queue_path,
                            site_policy_checked_utc=harness._utc(now),
                            now=now,
                        )
                    self.assertEqual(config["test_id"], case_id)
                    self.assertEqual(config["job_script"], str(job_script))
                    self.assertEqual(config["analysis_scripts"], [str(reconciler)])
                    harness.config.write_text(json.dumps(config), encoding="utf-8")
                    manifest_path = harness._create_manifest()
                    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                    self.assertEqual(manifest["test_id"], case_id)
                    with patch.dict(
                        reserve_globals,
                        {
                            "validate_clean_candidate_bundle": (
                                lambda *_args, **_kwargs: []
                            ),
                            "_require_scheduler_output_path": (
                                lambda *_args, **_kwargs: None
                            ),
                        },
                    ):
                        reservation = harness._reserve(
                            manifest_path,
                            reservation_id=str(uuid.uuid4()),
                            patch_clean_candidate_bundle=False,
                        )

                    def runner(command: list[str], **kwargs: object) -> None:
                        raw_root = Path(command[command.index("-d") + 1])
                        task_lines = []
                        for rank in range(int(case["mpi_ranks"])):
                            task_lines.append(
                                "PIC trusted GPU launch: "
                                f"rank={rank} host=nid000001 "
                                f"ROCR_VISIBLE_DEVICES={rank} "
                                "linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa\n"
                            )
                            for cycle in successor.REQUIRED_CYCLES:
                                for field in successor.REQUIRED_FIELDS:
                                    relative = successor._raw_relative_path(
                                        case, field=field, cycle=cycle, rank=rank
                                    )
                                    path = raw_root / relative
                                    path.parent.mkdir(parents=True, exist_ok=True)
                                    path.write_bytes(
                                        _binary_payload(
                                            case,
                                            field=field,
                                            cycle=cycle,
                                            rank=rank,
                                        )
                                    )
                        kwargs["stdout"].write(
                            "".join(task_lines).encode("utf-8")
                            + b"time=0.0025 cycle=1\n"
                            b"tlim=1 nlim=1\n"
                            b"Terminating on cycle limit\n"
                        )

                    job_id = str(20000 + index)
                    harness._launch(
                        manifest_path, reservation, job_id=job_id, runner=runner
                    )
                    with patch.object(
                        reconcile_frontier_job,
                        "_scheduler_result",
                        return_value=("COMPLETED", 30, 1, "0:0"),
                    ):
                        generic_event = reconcile_frontier_job.reconcile(
                            job_id=job_id,
                            ledger_jsonl=harness.ledger,
                            ledger_csv=harness.csv,
                            receipts_jsonl=harness.receipts,
                            mirror_jsonl=harness.mirror,
                            control_plane_dir=harness.control_plane_dir,
                            authorized_pic_root=harness.pic_root,
                            authorized_project_home_root=harness.project_home_root,
                        )
                    self.assertEqual(generic_event["scheduler_exit_code"], "0:0")
                    if index == 0:
                        sys.__dict__.pop("_pic_control_plane_bootstrapped", None)
                        with self.assertRaisesRegex(
                            ValueError, "requires bootstrapped installed execution"
                        ):
                            installed_producer.reconcile_q043(
                                job_id=job_id,
                                ledger_jsonl=harness.ledger,
                                ledger_csv=harness.csv,
                                receipts_jsonl=harness.receipts,
                                mirror_jsonl=harness.mirror,
                                control_plane_dir=harness.control_plane_dir,
                                authorized_pic_root=harness.pic_root,
                                authorized_project_home_root=harness.project_home_root,
                            )
                        sys._pic_control_plane_bootstrapped = True
                    produced = installed_producer.reconcile_q043(
                        job_id=job_id,
                        ledger_jsonl=harness.ledger,
                        ledger_csv=harness.csv,
                        receipts_jsonl=harness.receipts,
                        mirror_jsonl=harness.mirror,
                        control_plane_dir=harness.control_plane_dir,
                        authorized_pic_root=harness.pic_root,
                        authorized_project_home_root=harness.project_home_root,
                    )
                    artifact_dir = Path(manifest["artifact_dir"])
                    inventory = json.loads(
                        (artifact_dir / "artifact_inventory.json").read_text(
                            encoding="utf-8"
                        )
                    )
                    inventory_paths = {
                        record["path"] for record in inventory["files"]
                    }
                    self.assertTrue(
                        any(path.startswith("raw/") for path in inventory_paths)
                    )
                    self.assertFalse(
                        any(path.startswith("analysis/") for path in inventory_paths)
                    )
                    execution = {
                        "artifact_dir": str(artifact_dir),
                        "registered_execution_receipt": _binding(
                            Path(
                                produced["evidence"]["registered_execution_receipt"][
                                    "orion_path"
                                ]
                            )
                        ),
                        "terminal_receipt": _binding(
                            Path(
                                produced["evidence"]["terminal_receipt"][
                                    "orion_path"
                                ]
                            )
                        ),
                    }
                    raw_artifacts = []
                    for cycle in successor.REQUIRED_CYCLES:
                        for field in successor.REQUIRED_FIELDS:
                            for rank in range(int(case["mpi_ranks"])):
                                relative = successor._raw_relative_path(
                                    case, field=field, cycle=cycle, rank=rank
                                )
                                binding = _binding(artifact_dir / "raw" / relative)
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
                    with (
                        patch.object(
                            successor,
                            "_trusted_clean_candidate_revalidation",
                            side_effect=revalidation_fixture,
                        ),
                        patch.object(
                            successor,
                            "_trusted_mirrored_ledger_state",
                            REAL_TRUSTED_MIRRORED_LEDGER_STATE,
                        ),
                        patch.object(
                            successor,
                            "_trusted_producer_rederivation",
                            REAL_TRUSTED_PRODUCER_REDERIVATION,
                        ),
                    ):
                        admitted = successor.build_case_admission(
                            case_id=case_id,
                            candidate_binding=candidate_binding,
                            execution_binding=execution,
                            raw_artifacts=raw_artifacts,
                        )
                    self.assertTrue(
                        admitted["hardened_raw_oracle_result"][
                            "registered_execution_raw_oracle_check_pass"
                        ]
                    )
                    self.assertFalse(
                        admitted["authorization"]["publication_authorized"]
                    )
                    admissions.append(admitted)
            self.assertEqual(
                [item["case_id"] for item in admissions],
                [str(case["case_id"]) for case in cases],
            )
            self.assertEqual(len(admissions), 132)
            self.assertTrue(
                all(
                    value is False
                    for admission in admissions
                    for key, value in admission["authorization"].items()
                    if key.endswith("_authorized")
                )
            )
        finally:
            if original_bootstrap is None:
                sys.__dict__.pop("_pic_control_plane_bootstrapped", None)
            else:
                sys._pic_control_plane_bootstrapped = original_bootstrap
            (
                successor.AUTHORIZED_ORION_ROOT,
                successor.AUTHORIZED_PROJECT_HOME_ROOT,
                successor.AUTHORIZED_PROJECT_HOME_LEDGER_LEXICAL_ROOT,
            ) = original_roots
            if hasattr(harness, "temporary"):
                harness.tearDown()

    def test_installed_q043_action_is_runner_reachable(self) -> None:
        control_plane = sys.modules["control_plane_common"]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "control_plane"
            staging = root / "staging"
            staging.mkdir(parents=True)
            records = []
            for name in control_plane.CONTROL_PLANE_FILES:
                payload = (CONTROL_PLANE_SOURCE / name).read_bytes()
                path = staging / name
                path.write_bytes(payload)
                path.chmod(0o555 if path.suffix in {".py", ".sh"} else 0o444)
                records.append({"path": name, "sha256": _sha(payload)})
            version = hashlib.sha256(
                json.dumps(records, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest()
            inventory = {
                "schema_version": 1,
                "version": version,
                "files": records,
            }
            inventory_path = staging / "inventory.json"
            inventory_path.write_bytes(_json_bytes(inventory))
            inventory_path.chmod(0o444)
            installed = staging.with_name(version)
            staging.rename(installed)
            installed.chmod(0o555)
            result = subprocess.run(
                [
                    sys.executable,
                    str(installed / "run_control_plane.py"),
                    successor.Q043_PRODUCER_ENTRYPOINT,
                    "--help",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("--job-id", result.stdout)

    def test_q043_producer_publishes_paired_consumable_evidence_and_retries(
        self,
    ) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            evidence = fixture.last_produced["evidence"]
            for name in ("terminal_receipt", "registered_execution_receipt"):
                item = evidence[name]
                orion = Path(item["orion_path"])
                mirror = Path(item["project_home_mirror_path"])
                self.assertEqual(orion.read_bytes(), mirror.read_bytes())
                self.assertEqual(_sha(orion.read_bytes()), item["sha256"])
                self.assertFalse(bool(mirror.parent.stat().st_mode & 0o222))
            retried = q043_producer.publish_q043_registered_execution_evidence(
                fixture.last_event,
                fixture.last_mirror_ack,
                fixture._producer_inventory(),
                authorized_pic_root=fixture.authorized,
                authorized_project_home_root=fixture.project_home,
            )
            self.assertEqual(retried, evidence)
            successor.build_case_admission(
                case_id=case_id,
                candidate_binding=fixture.candidate,
                execution_binding=execution,
                raw_artifacts=raw,
            )

    def test_generic_q043_reconciliation_without_completion_is_not_admissible(
        self,
    ) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.inputs(case_id)
            generic_event = copy.deepcopy(fixture.last_event)
            generic_event.pop("trampoline_completion")
            generic_event["event_sha256"] = _sha(_json_bytes(generic_event))
            mirror_ack = copy.deepcopy(fixture.last_mirror_ack)
            mirror_ack["mirrored_event_sha256"] = generic_event["event_sha256"]
            with self.assertRaisesRegex(
                ValueError,
                "canonical mirrored-ledger anchor",
            ):
                q043_producer.derive_q043_registered_execution_evidence(
                    generic_event,
                    mirror_ack,
                    fixture._producer_inventory(),
                    authorized_pic_root=fixture.authorized,
                    authorized_project_home_root=fixture.project_home,
                )

    def test_missing_different_replaced_and_symlinked_project_home_mirrors_fail(
        self,
    ) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        for artifact_name, execution_key in (
            ("registered_execution_receipt", "registered_execution_receipt"),
            ("terminal_receipt", "terminal_receipt"),
        ):
            for mutation, expected in (
                ("missing", "escaped or is unavailable"),
                ("different", "SHA-256 drifted"),
                ("replaced", "submission mirror directory is mutable"),
                ("symlinked", "escaped or is unavailable|symlink alias"),
            ):
                with (
                    self.subTest(artifact=artifact_name, mutation=mutation),
                    tempfile.TemporaryDirectory() as temporary,
                ):
                    fixture = Fixture(Path(temporary))
                    _, execution, raw = fixture.inputs(case_id)
                    mirror = Path(
                        fixture.last_produced["evidence"][artifact_name][
                            "project_home_mirror_path"
                        ]
                    )
                    mirror.parent.chmod(0o700)
                    if mutation == "missing":
                        mirror.unlink()
                        mirror.parent.chmod(0o500)
                    elif mutation == "different":
                        _overwrite(mirror, b"different trusted mirror bytes\n")
                        _seal(mirror)
                        mirror.parent.chmod(0o500)
                    elif mutation == "replaced":
                        payload = mirror.read_bytes()
                        mirror.unlink()
                        mirror.write_bytes(payload)
                        _seal(mirror)
                    else:
                        mirror.unlink()
                        mirror.symlink_to(
                            Path(str(execution[execution_key]["path"]))
                        )
                        mirror.parent.chmod(0o500)
                    with self.assertRaisesRegex(successor.AdmissionError, expected):
                        successor.build_case_admission(
                            case_id=case_id,
                            candidate_binding=fixture.candidate,
                            execution_binding=execution,
                            raw_artifacts=raw,
                        )

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
                )
            admitted = successor.build_case_admission(
                case_id=case_id,
                candidate_binding=fixture.candidate,
                execution_binding=execution,
                raw_artifacts=raw,
            )
            source_local = admitted["hardened_raw_oracle_result"]["source_local_oracle_result"]
            with self.assertRaisesRegex(
                successor.AdmissionError, "requires registered Q043 matrix"
            ):
                successor.validate_downstream_q023_q019_prerequisite(
                    source_local,
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
            receipt = json.loads(
                Path(str(execution["registered_execution_receipt"]["path"])).read_text()
            )
            receipt["mpi_evidence"]["tasks"] = 1
            fixture.reseal_receipt(execution, receipt)
            with self.assertRaisesRegex(
                successor.AdmissionError, "command or MPI|re-derived different evidence"
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
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
                )
            duplicate = copy.deepcopy(raw)
            duplicate.append(copy.deepcopy(duplicate[-1]))
            with self.assertRaisesRegex(successor.AdmissionError, "path reused"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=duplicate,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            raw_bin = Path(execution["artifact_dir"], "raw/bin")
            raw_bin.chmod(0o755)
            (raw_bin / "undeclared.00002.bin").write_bytes(b"undeclared\n")
            raw_bin.chmod(0o555)
            with self.assertRaisesRegex(
                successor.AdmissionError, "artifact tree differs|expected inventory"
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            first = Path(execution["artifact_dir"], "raw", str(raw[0]["path"]))
            second = Path(execution["artifact_dir"], "raw", str(raw[1]["path"]))
            second.parent.chmod(0o755)
            second.unlink()
            os.link(first, second)
            second.parent.chmod(0o555)
            raw[1]["sha256"] = raw[0]["sha256"]
            raw[1]["byte_count"] = raw[0]["byte_count"]
            with self.assertRaisesRegex(
                successor.AdmissionError, "not one read-only file|unaliased"
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

    def test_wrong_candidate_deck_and_executable_fail(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            bad_candidate = copy.deepcopy(fixture.candidate)
            bad_candidate["source_archive"]["sha256"] = "9" * 64
            with self.assertRaisesRegex(
                successor.AdmissionError, "SHA-256 drifted|lacks an execute bit"
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=bad_candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            bad_candidate = copy.deepcopy(fixture.candidate)
            malformed = b"not a source archive\n"
            _overwrite(Path(str(bad_candidate["source_archive"]["path"])), malformed)
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
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            bad_candidate = copy.deepcopy(fixture.candidate)
            _overwrite(
                Path(str(bad_candidate["executable"]["path"])), b"wrong executable\n"
            )
            with self.assertRaisesRegex(
                successor.AdmissionError, "SHA-256 drifted|lacks an execute bit"
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=bad_candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            launch_deck_path = (
                Path(str(fixture.last_event["manifest_path"])).parent
                / "snapshot"
                / f"{case_id}.athinput"
            )
            bad_launch_deck = launch_deck_path.read_bytes() + b"# drift\n"
            _overwrite(launch_deck_path, bad_launch_deck)
            with self.assertRaisesRegex(successor.AdmissionError, "input-deck|snapshot"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            receipt = json.loads(
                Path(str(execution["registered_execution_receipt"]["path"])).read_text()
            )
            receipt["deck_sha256"] = "9" * 64
            fixture.reseal_receipt(execution, receipt)
            with self.assertRaisesRegex(successor.AdmissionError, "cross-link|re-derived"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

    def test_self_issued_receipt_without_canonical_reconciliation_fails(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            fixture.ledger_records.clear()
            with self.assertRaisesRegex(
                successor.AdmissionError, "canonical mirrored-ledger reconciliation"
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

    def test_self_reported_fake_executable_fails_retained_elf_validation(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            executable = fixture.candidate["executable"]
            fake = b"self-reported fake executable\n"
            _overwrite(Path(str(executable["path"])), fake)
            _seal(Path(str(executable["path"])), executable=True)
            executable["sha256"] = _sha(fake)
            executable["byte_count"] = len(fake)
            manifest = json.loads(
                Path(str(fixture.candidate["clean_candidate_manifest"]["path"])).read_text()
            )
            manifest["build"]["executable_sha256"] = executable["sha256"]
            _rewrite_json(fixture.candidate["clean_candidate_manifest"], manifest)
            with self.assertRaisesRegex(successor.AdmissionError, "not an ELF binary"):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            executable = fixture.candidate["executable"]
            fake_elf = b"\x7fELFself-reported fake executable\n"
            _overwrite(Path(str(executable["path"])), fake_elf)
            _seal(Path(str(executable["path"])), executable=True)
            executable["sha256"] = _sha(fake_elf)
            executable["byte_count"] = len(fake_elf)
            manifest = json.loads(
                Path(str(fixture.candidate["clean_candidate_manifest"]["path"])).read_text()
            )
            manifest["build"]["executable_sha256"] = executable["sha256"]
            _rewrite_json(fixture.candidate["clean_candidate_manifest"], manifest)
            with patch.object(
                successor,
                "_trusted_clean_candidate_revalidation",
                REAL_TRUSTED_CANDIDATE_REVALIDATION,
            ), self.assertRaisesRegex(
                successor.AdmissionError, "failed trusted fixed-root build revalidation"
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

    def test_modified_deposition_and_output_sources_fail_archive_closure(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        for relative in (
            "src/particles/particles_moments.cpp",
            "src/particles/particles_pushers.cpp",
            "src/outputs/binary.cpp",
            "tst/publication/analyze_q011_section54_outputs.py",
        ):
            with self.subTest(relative=relative), tempfile.TemporaryDirectory() as temporary:
                fixture = Fixture(Path(temporary))
                _, execution, raw = fixture.inputs(case_id)
                archive = fixture.candidate["source_archive"]
                payload = _source_archive_payload({relative: b"review-bypassing drift\n"})
                _overwrite(Path(str(archive["path"])), payload)
                _seal(Path(str(archive["path"])))
                archive["sha256"] = _sha(payload)
                archive["byte_count"] = len(payload)
                manifest = json.loads(
                    Path(
                        str(fixture.candidate["clean_candidate_manifest"]["path"])
                    ).read_text()
                )
                manifest["source"]["archive_sha256"] = archive["sha256"]
                manifest["build"]["source_archive_sha256"] = archive["sha256"]
                _rewrite_json(fixture.candidate["clean_candidate_manifest"], manifest)
                with self.assertRaisesRegex(
                    successor.AdmissionError,
                    rf"source archive member drifted: {re.escape(relative)}",
                ):
                    successor.build_case_admission(
                        case_id=case_id,
                        candidate_binding=fixture.candidate,
                        execution_binding=execution,
                        raw_artifacts=raw,
                    )

    def test_post_receipt_raw_replacement_fails_inventory_seal(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            case, execution, raw = fixture.inputs(case_id)
            artifact = next(
                item
                for item in raw
                if item["cycle"] == 1 and item["field"] == "prtcl_jx"
            )
            fixture.replace_raw_artifact(
                case,
                execution,
                artifact,
                value=_field_value(case, "prtcl_jx", 1) + 0.25,
            )
            with self.assertRaisesRegex(
                successor.AdmissionError,
                "installed Q043 producer|Sealed launch artifact",
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

    def test_entire_run_root_substitution_fails_external_completion_anchor(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, _ = fixture.inputs(case_id)
            artifact_dir = Path(str(execution["artifact_dir"]))
            detached = artifact_dir.with_name(f"{artifact_dir.name}.original")
            artifact_dir.rename(detached)
            shutil.copytree(detached, artifact_dir, copy_function=shutil.copy2)
            with self.assertRaisesRegex(
                ValueError,
                "trampoline completion|Launch artifact root",
            ):
                q043_producer.publish_q043_registered_execution_evidence(
                    fixture.last_event,
                    fixture.last_mirror_ack,
                    fixture._producer_inventory(),
                    authorized_pic_root=fixture.authorized,
                    authorized_project_home_root=fixture.project_home,
                )

    def test_dual_completion_receipt_directory_replacement_fails_ledger_anchor(
        self,
    ) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with self.subTest(boundary="before-reconciliation"), tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            with self.assertRaisesRegex(
                ValueError,
                "canonical mirrored-ledger anchor",
            ):
                fixture.inputs(
                    case_id,
                    after_completion=lambda launched: (
                        launched.replace_completion_receipt_directories(
                            malicious=True
                        )
                    ),
                )

        with self.subTest(boundary="before-admission"), tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            fixture.replace_completion_receipt_directories(malicious=True)
            with self.assertRaisesRegex(
                successor.AdmissionError,
                "canonical mirrored-ledger anchor|installed Q043 producer",
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

    def test_raw_replacement_between_hash_and_parser_is_rejected(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            target = Path(
                str(execution["artifact_dir"]), "raw", str(raw[0]["path"])
            )
            original_parse = successor.binary.parse_athenak_binary_bytes
            swapped = False

            def parse_after_swap(payload: bytes, *, source: str):
                nonlocal swapped
                if not swapped:
                    replacement = bytearray(target.read_bytes())
                    replacement[-1] ^= 1
                    _overwrite(target, bytes(replacement))
                    _seal(target)
                    swapped = True
                return original_parse(payload, source=source)

            with patch.object(
                successor.binary,
                "parse_athenak_binary_bytes",
                side_effect=parse_after_swap,
            ), self.assertRaisesRegex(
                successor.AdmissionError,
                "retained raw bytes|exact-byte analysis|artifact",
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )
            self.assertTrue(swapped)

    def test_raw_replacement_before_exact_byte_oracle_is_rejected(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            target_record = next(item for item in raw if item["cycle"] == 1)
            target = Path(
                str(execution["artifact_dir"]), "raw", str(target_record["path"])
            )
            original_oracle = oracle.analyze_raw_case_bytes
            swapped = False

            def analyze_after_swap(
                observed_case_id: str,
                snapshots: dict[str, tuple[tuple[str, bytes], ...]],
            ) -> dict[str, object]:
                nonlocal swapped
                replacement = bytearray(target.read_bytes())
                replacement[-1] ^= 1
                _overwrite(target, bytes(replacement))
                _seal(target)
                swapped = True
                return original_oracle(observed_case_id, snapshots)

            with (
                patch.object(
                    oracle, "analyze_raw_case_bytes", side_effect=analyze_after_swap
                ),
                patch.object(
                    oracle,
                    "analyze_raw_case",
                    side_effect=AssertionError("path-reopening oracle must not execute"),
                ),
                self.assertRaisesRegex(
                    successor.AdmissionError,
                    "retained raw bytes|exact-byte analysis|artifact",
                ),
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )
            self.assertTrue(swapped)

    def test_invalid_cycle_zero_values_fail_scientific_validation(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        for value in (1.0, float("nan")):
            with self.subTest(value=value), tempfile.TemporaryDirectory() as temporary:
                fixture = Fixture(Path(temporary))
                _, execution, raw = fixture.inputs(
                    case_id,
                    raw_value_overrides={(0, "prtcl_rho", 0): value},
                )
                with self.assertRaisesRegex(
                    successor.AdmissionError,
                    "must be exactly finite zero|malformed AthenaK binary output",
                ):
                    successor.build_case_admission(
                        case_id=case_id,
                        candidate_binding=fixture.candidate,
                        execution_binding=execution,
                        raw_artifacts=raw,
                    )

    def test_bad_terminal_state_and_post_freeze_stdout_replacement_fail(self) -> None:
        case_id = "q043-current-oracle-d1-coarse-ppc1-single-cvr100"
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            receipt = json.loads(
                Path(str(execution["registered_execution_receipt"]["path"])).read_text()
            )
            receipt["slurm_terminal_state"] = "FAILED"
            _rewrite_json(execution["registered_execution_receipt"], receipt)
            with self.assertRaisesRegex(
                successor.AdmissionError, "registered.execution|receipt"
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            stdout_path = Path(str(execution["artifact_dir"])) / "athena_stdout.txt"
            bad_stdout = stdout_path.read_bytes().replace(
                b"Terminating on cycle limit", b"Terminating on wall clock limit"
            )
            _overwrite(stdout_path, bad_stdout)
            with self.assertRaisesRegex(
                successor.AdmissionError, "installed Q043 producer|artifact"
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
                )

    def test_installed_producer_requires_exact_trusted_rank_and_exit_stdout(self) -> None:
        case_id = "q043-current-oracle-d2-coarse-ppc1-split_x1-cvr100"
        valid = (
            "PIC trusted GPU launch: rank=0 host=nid000001 ROCR_VISIBLE_DEVICES=0 "
            "linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa\n"
            "PIC trusted GPU launch: rank=1 host=nid000001 ROCR_VISIBLE_DEVICES=1 "
            "linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa\n"
            "time=0.0025 cycle=1\n"
            "tlim=1 nlim=1\n"
            "Terminating on cycle limit\n"
            f"Q043_REGISTERED_EXECUTION case_id={case_id} "
            "mpi_world_size=2 rank_ids=0,1\n"
            "Q043_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0\n"
        ).encode("utf-8")
        evidence = q043_producer._trusted_wrapper_evidence(
            valid, case_id=case_id, ranks=2, stdout_sha256=_sha(valid)
        )
        self.assertEqual(evidence["observed_rank_ids"], [0, 1])
        self.assertEqual(evidence["exit_code"], 0)
        for bad in (
            valid.replace(b"rank_ids=0,1", b"rank_ids=0,0"),
            valid.replace(b"exit_code=0", b"exit_code=1"),
            valid.replace(
                b"PIC trusted GPU launch: rank=1",
                b"PIC trusted GPU launch: rank=0",
            ),
        ):
            with self.assertRaisesRegex(ValueError, "rank/exit|incomplete or duplicated"):
                q043_producer._trusted_wrapper_evidence(
                    bad, case_id=case_id, ranks=2, stdout_sha256=_sha(bad)
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
                )

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            _, execution, raw = fixture.inputs(case_id)
            path = Path(execution["artifact_dir"], "raw", str(raw[0]["path"]))
            retained = fixture.authorized / "retained.bin"
            retained.write_bytes(path.read_bytes())
            path.parent.chmod(0o755)
            path.unlink()
            path.symlink_to(retained)
            path.parent.chmod(0o555)
            with self.assertRaisesRegex(
                successor.AdmissionError,
                "escaped or is unavailable|symlink alias|unsupported entry",
            ):
                successor.build_case_admission(
                    case_id=case_id,
                    candidate_binding=fixture.candidate,
                    execution_binding=execution,
                    raw_artifacts=raw,
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
                    Path(execution["artifact_dir"], "raw", str(item["path"]))
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
