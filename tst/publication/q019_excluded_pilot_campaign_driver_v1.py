#!/opt/cray/pe/python/3.11.7/bin/python3
"""Run and admit the four registered Q019 runtime-controller pilots."""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import time
import uuid
from typing import Mapping, Sequence

from tst.publication import (
    q019_hardened_installed_control_plane_registered_admission_v1
    as registered_admission,
)


PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
PROJECT_HOME_LEDGER_ROOT = Path("/ccs/proj/ast207/proj-shared/PIC")
PYTHON = Path("/opt/cray/pe/python/3.11.7/bin/python3")
LEDGER_JSONL = PIC_ROOT / "ledger/node_hours.jsonl"
LEDGER_CSV = PIC_ROOT / "ledger/node_hours.csv"
RECEIPTS_JSONL = PIC_ROOT / "ledger/mirror_receipts.jsonl"
MIRROR_JSONL = PROJECT_HOME_LEDGER_ROOT / "ledger/node_hours.jsonl"
JOB_ID_PATTERN = re.compile(r"Submitted ([1-9][0-9]*) with reservation ")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
COMMIT_PATTERN = re.compile(r"[0-9a-f]{40}")
UUID_PATTERN = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-"
    r"[89ab][0-9a-f]{3}-[0-9a-f]{12}"
)
SCHEMA_VERSION = 1
EXECUTION_INDEX_RECORD_TYPE = "q019_excluded_pilot_execution_index_v1"
EXPECTED_ARTIFACTS = (
    "q019-controller-pilot-2d-baseline",
    "q019-controller-pilot-2d-instrumented",
    "q019-controller-pilot-3d-baseline",
    "q019-controller-pilot-3d-instrumented",
)
EXPECTED_CASES = {
    "q019-controller-pilot-2d-baseline": "q019-fr-grid-k8-rho1em05-s0",
    "q019-controller-pilot-2d-instrumented": "q019-fr-grid-k8-rho1em05-s0",
    "q019-controller-pilot-3d-baseline": "q019-fr-3d-onset-small-s0",
    "q019-controller-pilot-3d-instrumented": "q019-fr-3d-onset-small-s0",
}
EXPECTED_AUTHORIZATIONS = {
    "q019-controller-pilot-2d-baseline": "q019-pilot-01-2d-baseline-v1",
    "q019-controller-pilot-2d-instrumented": "q019-pilot-02-2d-instrumented-v1",
    "q019-controller-pilot-3d-baseline": "q019-pilot-03-3d-baseline-v1",
    "q019-controller-pilot-3d-instrumented": "q019-pilot-04-3d-instrumented-v1",
}
EXPECTED_NODES = {
    "q019-controller-pilot-2d-baseline": 4,
    "q019-controller-pilot-2d-instrumented": 4,
    "q019-controller-pilot-3d-baseline": 16,
    "q019-controller-pilot-3d-instrumented": 16,
}
MAXIMUM_ELAPSED_SECONDS = {
    "q019-controller-pilot-2d-baseline": 1800,
    "q019-controller-pilot-2d-instrumented": 1800,
    "q019-controller-pilot-3d-baseline": 3600,
    "q019-controller-pilot-3d-instrumented": 3600,
}
MAXIMUM_STORAGE_BYTES = {
    "q019-controller-pilot-2d-baseline": 64 * 1024**3,
    "q019-controller-pilot-2d-instrumented": 64 * 1024**3,
    "q019-controller-pilot-3d-baseline": 256 * 1024**3,
    "q019-controller-pilot-3d-instrumented": 256 * 1024**3,
}
AUTHORIZATION_BOUNDARY = {
    "launch_authorized": False,
    "scheduler_submission_authorized": False,
    "policy_mutation_authorized": False,
    "q019_qualification_authorized": False,
    "nonlinear_saturation_claim_authorized": False,
    "scientific_claim_authorized": False,
    "publication_authorized": False,
}
EXECUTION_PROVENANCE_KEYS = frozenset(
    {
        "source_commit",
        "source_bundle_sha256",
        "source_archive_sha256",
        "executable_sha256",
        "environment_sha256",
        "deck_sha256",
        "pre_submit_manifest",
        "clean_candidate_manifest",
        "source_archive",
        "executable_snapshot",
        "execution_receipt",
        "terminal_receipt",
        "installed_producer",
        "reconciliation_event_sha256",
        "reconciliation_mirror_ack_sha256",
        "raw_inventory_sha256",
        "retained_raw_bindings",
        "reduction_binding",
    }
)


class DriverError(RuntimeError):
    """Reject drifted live identities, ambiguous resume state, or bad pilots."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise DriverError(message)


def _utc_now() -> datetime:
    return datetime.now(timezone.utc).replace(microsecond=0)


def _utc(value: datetime) -> str:
    return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _canonical_sha256(value: object) -> str:
    payload = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _ledger_event_sha256(event: Mapping[str, object]) -> str:
    payload = {key: value for key, value in event.items() if key != "event_sha256"}
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _binding(path: Path) -> dict[str, object]:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and not before.st_mode & 0o222,
            f"Q019 pilot binding is not one read-only regular file: {path}",
        )
        digest = hashlib.sha256()
        byte_count = 0
        while chunk := os.read(descriptor, 1024 * 1024):
            digest.update(chunk)
            byte_count += len(chunk)
        after = os.fstat(descriptor)
        current = path.stat(follow_symlinks=False)
        identity = lambda item: (
            item.st_dev,
            item.st_ino,
            item.st_mode,
            item.st_nlink,
            item.st_size,
            item.st_mtime_ns,
            item.st_ctime_ns,
        )
        _require(
            identity(before) == identity(after)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and byte_count == after.st_size,
            f"Q019 pilot binding changed while reading: {path}",
        )
        return {
            "path": str(path),
            "sha256": digest.hexdigest(),
            "byte_count": byte_count,
        }
    finally:
        os.close(descriptor)


def _run(
    argv: Sequence[str | Path],
    *,
    cwd: Path | None = None,
    capture: bool = True,
) -> str:
    command = [str(item) for item in argv]
    result = subprocess.run(
        command,
        cwd=str(cwd) if cwd is not None else None,
        check=False,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
        env={
            **os.environ,
            "HOME": "/",
            "LANG": "C",
            "LC_ALL": "C",
            "SLURM_CLUSTERS": "frontier",
        },
    )
    if result.returncode != 0:
        raise DriverError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout:\n{result.stdout or ''}\nstderr:\n{result.stderr or ''}"
        )
    return (result.stdout or "").strip()


def _atomic_json(path: Path, value: object, *, mode: int = 0o444) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _json_bytes(value)
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        mode,
    )
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError(f"short write: {path}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def _replace_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _json_bytes(value)
    temporary = path.parent / f".{path.name}.tmp-{uuid.uuid4()}"
    descriptor = os.open(
        temporary,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o600,
    )
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError(f"short write: {temporary}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(temporary, path)
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def _write_empty_read_only(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o444,
    )
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _strict_equal(left: object, right: object) -> bool:
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return set(left) == set(right) and all(
            _strict_equal(left[key], right[key]) for key in left
        )
    if isinstance(left, list):
        return len(left) == len(right) and all(
            _strict_equal(a, b) for a, b in zip(left, right)
        )
    return left == right


def _artifact_payload_bytes(
    binding: Mapping[str, object], *, artifact_root: Path
) -> int:
    path = Path(str(binding.get("path", "")))
    _require(
        path == artifact_root / "artifact_inventory.json",
        "Q019 pilot artifact inventory path drifted",
    )
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        chunks = []
        while chunk := os.read(descriptor, 1024 * 1024):
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
        current = path.stat(follow_symlinks=False)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and not before.st_mode & 0o222
            and (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
            == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and len(payload) == binding.get("byte_count") == after.st_size
            and hashlib.sha256(payload).hexdigest() == binding.get("sha256"),
            "Q019 pilot artifact inventory bytes drifted",
        )
    finally:
        os.close(descriptor)
    try:
        inventory = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise DriverError("Q019 pilot artifact inventory is not UTF-8 JSON") from error
    records = inventory.get("files") if type(inventory) is dict else None
    _require(
        type(inventory) is dict
        and set(inventory) == {"schema_version", "files"}
        and inventory.get("schema_version") == 1
        and type(records) is list
        and all(
            type(record) is dict
            and set(record) == {"path", "sha256", "size"}
            and type(record["path"]) is str
            and bool(record["path"])
            and type(record["sha256"]) is str
            and SHA256_PATTERN.fullmatch(record["sha256"]) is not None
            and type(record["size"]) is int
            and record["size"] >= 0
            for record in records
        )
        and len({record["path"] for record in records}) == len(records),
        "Q019 pilot artifact inventory schema drifted",
    )
    return sum(record["size"] for record in records)


def _execution_provenance(admission: Mapping[str, object]) -> dict[str, object]:
    lineage = admission.get("execution_lineage")
    _require(
        type(lineage) is dict
        and set(lineage) == EXECUTION_PROVENANCE_KEYS
        and type(lineage.get("source_commit")) is str
        and COMMIT_PATTERN.fullmatch(str(lineage["source_commit"])) is not None
        and all(
            type(lineage.get(key)) is str
            and SHA256_PATTERN.fullmatch(str(lineage[key])) is not None
            for key in (
                "source_bundle_sha256",
                "source_archive_sha256",
                "executable_sha256",
                "environment_sha256",
                "deck_sha256",
                "reconciliation_event_sha256",
                "reconciliation_mirror_ack_sha256",
                "raw_inventory_sha256",
            )
        )
        and all(
            type(lineage.get(key)) is dict
            for key in (
                "pre_submit_manifest",
                "clean_candidate_manifest",
                "source_archive",
                "executable_snapshot",
                "execution_receipt",
                "terminal_receipt",
                "installed_producer",
                "reduction_binding",
            )
        )
        and type(lineage.get("retained_raw_bindings")) is list,
        "Q019 pilot admission execution provenance is malformed",
    )
    return dict(lineage)


def _reopen_attempt_admission(
    attempt: Mapping[str, object],
) -> dict[str, object]:
    binding = attempt.get("admission")
    _require(type(binding) is dict, "Q019 pilot admission binding is malformed")
    path = Path(str(binding.get("path", "")))
    _require(
        _binding(path) == binding,
        "Q019 pilot admission binding bytes or filesystem identity drifted",
    )
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise DriverError("Q019 pilot admission is not stable UTF-8 JSON") from error
    try:
        return registered_admission.validate_admission(value)
    except registered_admission.RegisteredAdmissionError as error:
        raise DriverError("Q019 pilot admission did not fully rederive") from error


def _validate_attempt_admission_projection(
    attempt: Mapping[str, object], admission: Mapping[str, object]
) -> None:
    expected_facts = {
        key: admission[key]
        for key in (
            "record_type",
            "case_id",
            "artifact_root",
            "registered_execution_identity",
            "execution_profile",
            "runtime_completion",
            "raw_science_admission_eligible",
            "problem_reported_saturation_evidence_eligible",
            "saturation_evidence_eligible",
            "authorization",
        )
    }
    _require(
        _strict_equal(attempt.get("admission_facts"), expected_facts)
        and _strict_equal(
            attempt.get("execution_provenance"), _execution_provenance(admission)
        ),
        "Q019 pilot attempt projection differs from rederived admission",
    )


def build_execution_index(
    attempts: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Validate and index four admitted excluded pilots without claim authority."""
    _require(
        type(attempts) in {list, tuple} and len(attempts) == len(EXPECTED_ARTIFACTS),
        "Q019 pilot execution index requires four attempts",
    )
    normalized: list[dict[str, object]] = []
    identities = {
        "submission_id": set(),
        "job_id": set(),
        "authorization_id": set(),
        "artifact_root": set(),
        "admission_sha256": set(),
        "reconciliation_event_sha256": set(),
    }
    for index, raw in enumerate(attempts):
        _require(type(raw) is dict, "Q019 pilot attempt record must be an object")
        attempt = dict(raw)
        artifact_id = EXPECTED_ARTIFACTS[index]
        case_id = EXPECTED_CASES[artifact_id]
        authorization_id = EXPECTED_AUTHORIZATIONS[artifact_id]
        expected_nodes = EXPECTED_NODES[artifact_id]
        _require(
            attempt.get("artifact_id") == artifact_id
            and attempt.get("source_case_id") == case_id
            and attempt.get("attempt_index") == index + 1
            and attempt.get("authorization_id") == authorization_id
            and type(attempt.get("submission_id")) is str
            and UUID_PATTERN.fullmatch(str(attempt["submission_id"])) is not None
            and type(attempt.get("job_id")) is str
            and str(attempt["job_id"]).isdigit()
            and attempt.get("artifact_root")
            == str(
                PIC_ROOT
                / "runs/q019_nonlinear_bell_registered_successor_v1"
                / str(attempt.get("submission_id"))
            ),
            f"{artifact_id}: Q019 pilot execution identity drifted",
        )
        admission = attempt.get("admission")
        admission_facts = attempt.get("admission_facts")
        execution_provenance = attempt.get("execution_provenance")
        inventory = attempt.get("artifact_inventory")
        payload_bytes = attempt.get("artifact_payload_bytes")
        event = attempt.get("reconciliation_event")
        _require(
            type(admission) is dict
            and set(admission) == {"path", "sha256", "byte_count"}
            and type(admission["sha256"]) is str
            and SHA256_PATTERN.fullmatch(admission["sha256"]) is not None
            and type(admission["byte_count"]) is int
            and admission["byte_count"] > 0
            and type(admission_facts) is dict
            and type(execution_provenance) is dict
            and set(execution_provenance) == EXECUTION_PROVENANCE_KEYS
            and type(inventory) is dict
            and set(inventory) == {"path", "sha256", "byte_count"}
            and type(inventory["sha256"]) is str
            and SHA256_PATTERN.fullmatch(inventory["sha256"]) is not None
            and type(inventory["byte_count"]) is int
            and inventory["byte_count"] > 0
            and type(payload_bytes) is int
            and 0 < payload_bytes <= MAXIMUM_STORAGE_BYTES[artifact_id]
            and type(event) is dict,
            f"{artifact_id}: Q019 pilot admission or storage binding drifted",
        )
        _require(
            execution_provenance.get("reconciliation_event_sha256")
            == admission_facts.get("registered_execution_identity", {}).get(
                "reconciliation_event_sha256"
            )
            and execution_provenance.get("reconciliation_mirror_ack_sha256")
            == admission_facts.get("registered_execution_identity", {}).get(
                "reconciliation_mirror_ack_sha256"
            ),
            f"{artifact_id}: Q019 pilot provenance identity drifted",
        )
        _execution_provenance({"execution_lineage": execution_provenance})
        registered_identity = admission_facts.get("registered_execution_identity")
        execution_profile = admission_facts.get("execution_profile")
        completion = admission_facts.get("runtime_completion")
        admission_authorization = admission_facts.get("authorization")
        _require(
            admission_facts.get("record_type")
            == "q019_hardened_installed_control_plane_registered_admission_v1"
            and admission_facts.get("case_id") == case_id
            and admission_facts.get("artifact_root") == attempt["artifact_root"]
            and type(registered_identity) is dict
            and registered_identity.get("submission_id") == attempt["submission_id"]
            and registered_identity.get("registered_science_authorization_id")
            == authorization_id
            and registered_identity.get("slurm_job_id") == attempt["job_id"]
            and type(registered_identity.get("reconciliation_event_sha256")) is str
            and type(execution_profile) is dict
            and execution_profile.get("kind") == "runtime_controller_overlay"
            and execution_profile.get("artifact_id") == artifact_id
            and execution_profile.get("source_case_id") == case_id
            and execution_profile.get("authority") == "excluded_pilot_only"
            and execution_profile.get("expected_stop_reason") == 1903
            and execution_profile.get("saturation_evidence_eligible") is False
            and type(completion) is dict
            and completion.get("run_completion_status")
            == "completed_not_acceptance_eligible"
            and completion.get("problem_stop_requested") is True
            and completion.get("stop_reason_code") == "1903"
            and completion.get("runtime_controller_trigger_cycle") == 20
            and completion.get("process_exit_code") == 0
            and completion.get("scheduler_terminal_state") == "COMPLETED"
            and completion.get("trusted_execution_binding_present") is True
            and admission_facts.get("raw_science_admission_eligible") is True
            and admission_facts.get("problem_reported_saturation_evidence_eligible")
            is False
            and admission_facts.get("saturation_evidence_eligible") is False
            and type(admission_authorization) is dict
            and bool(admission_authorization)
            and all(value is False for value in admission_authorization.values()),
            f"{artifact_id}: Q019 pilot admitted controller facts drifted",
        )
        _require(
            event.get("event_type") == "reconciliation"
            and event.get("campaign") == "q019_nonlinear_bell_registered_successor_v1"
            and event.get("test_id") == case_id
            and event.get("registered_science_authorization_id")
            == attempt["authorization_id"]
            and event.get("submission_id") == attempt["submission_id"]
            and event.get("job_id") == attempt["job_id"]
            and event.get("artifact_dir") == attempt.get("artifact_root")
            and event.get("state") == "COMPLETED"
            and event.get("reconciled") is True
            and event.get("scheduler_exit_code") == "0:0"
            and event.get("requested_nodes") == expected_nodes
            and event.get("scheduler_reported_allocated_nodes") == expected_nodes
            and type(event.get("elapsed_seconds")) is int
            and event["elapsed_seconds"] > 0
            and event["elapsed_seconds"] <= MAXIMUM_ELAPSED_SECONDS[artifact_id]
            and type(event.get("billed_nodes")) is int
            and event["billed_nodes"] == expected_nodes
            and type(event.get("consumed_node_hours")) in {int, float}
            and math.isclose(
                float(event["consumed_node_hours"]),
                event["billed_nodes"] * event["elapsed_seconds"] / 3600.0,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            and type(event.get("event_sha256")) is str
            and SHA256_PATTERN.fullmatch(event["event_sha256"]) is not None,
            f"{artifact_id}: Q019 pilot reconciliation event drifted",
        )
        _require(
            event["event_sha256"] == _ledger_event_sha256(event)
            and registered_identity["reconciliation_event_sha256"]
            == event["event_sha256"],
            f"{artifact_id}: Q019 pilot reconciliation hash binding drifted",
        )
        for key, value in (
            ("submission_id", attempt["submission_id"]),
            ("job_id", attempt["job_id"]),
            ("authorization_id", attempt["authorization_id"]),
            ("artifact_root", attempt["artifact_root"]),
            ("admission_sha256", admission["sha256"]),
            ("reconciliation_event_sha256", event["event_sha256"]),
        ):
            _require(
                value not in identities[key],
                f"{artifact_id}: Q019 pilot execution identity reused: {key}",
            )
            identities[key].add(value)
        normalized.append(attempt)
    pair_records = []
    for dimension, pair in (
        (2, normalized[:2]),
        (3, normalized[2:]),
    ):
        baseline, instrumented = pair
        _require(
            baseline["source_case_id"] == instrumented["source_case_id"]
            and baseline["artifact_id"].endswith("-baseline")
            and instrumented["artifact_id"].endswith("-instrumented"),
            f"Q019 {dimension}D pilot pair identity drifted",
        )
        base_elapsed = baseline["reconciliation_event"]["elapsed_seconds"]
        instrumented_elapsed = instrumented["reconciliation_event"]["elapsed_seconds"]
        pair_records.append(
            {
                "dimension": dimension,
                "source_case_id": baseline["source_case_id"],
                "baseline_artifact_id": baseline["artifact_id"],
                "instrumented_artifact_id": instrumented["artifact_id"],
                "baseline_elapsed_seconds": base_elapsed,
                "instrumented_elapsed_seconds": instrumented_elapsed,
                "instrumentation_elapsed_ratio": instrumented_elapsed / base_elapsed,
                "baseline_artifact_bytes": baseline["artifact_payload_bytes"],
                "instrumented_artifact_bytes": instrumented["artifact_payload_bytes"],
                "resource_measurement_only": True,
                "scientific_selection_authorized": False,
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": EXECUTION_INDEX_RECORD_TYPE,
        "status": "complete_registered_excluded_pilot_execution_index_non_authorizing",
        "campaign": "q019_nonlinear_bell_registered_successor_v1",
        "attempt_count": len(normalized),
        "attempts": normalized,
        "attempt_bindings_sha256": _canonical_sha256(normalized),
        "pairs": pair_records,
        "pair_bindings_sha256": _canonical_sha256(pair_records),
        "total_consumed_node_hours": sum(
            float(item["reconciliation_event"]["consumed_node_hours"])
            for item in normalized
        ),
        "total_artifact_bytes": sum(
            int(item["artifact_payload_bytes"]) for item in normalized
        ),
        "saturation_evidence_eligible": False,
        "production_resource_freeze_authorized": False,
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def validate_execution_index(value: object) -> dict[str, object]:
    _require(
        type(value) is dict
        and value.get("record_type") == EXECUTION_INDEX_RECORD_TYPE,
        "Q019 pilot execution index identity drifted",
    )
    rebuilt = build_execution_index(value.get("attempts"))
    _require(
        _strict_equal(value, rebuilt),
        "Q019 pilot execution index derived fields or authority drifted",
    )
    return rebuilt


def validate_execution_index_files(value: object) -> dict[str, object]:
    rebuilt = validate_execution_index(value)
    for attempt in rebuilt["attempts"]:
        admission = _reopen_attempt_admission(attempt)
        _validate_attempt_admission_projection(attempt, admission)
        _require(
            _artifact_payload_bytes(
                attempt["artifact_inventory"],
                artifact_root=Path(str(attempt["artifact_root"])),
            )
            == attempt["artifact_payload_bytes"],
            f"{attempt['artifact_id']}: Q019 artifact payload-byte total drifted",
        )
    return rebuilt


class Driver:
    def __init__(self, arguments: argparse.Namespace) -> None:
        self.source_root = arguments.source_root.resolve(strict=True)
        self.source_commit = arguments.source_commit
        self.control_plane_version = arguments.control_plane_version
        self.final_bindings_path = arguments.final_bindings.resolve(strict=True)
        self.final_bindings_sha256 = arguments.final_bindings_sha256
        self.reviewed_policy_path = arguments.reviewed_policy.resolve(strict=True)
        self.reviewed_policy_sha256 = arguments.reviewed_policy_sha256
        self.baseline_policy_sha256 = arguments.baseline_policy_sha256
        self.baseline_promotion_sha256 = arguments.baseline_promotion_sha256
        self.control_plane_root = PIC_ROOT / "control_plane" / self.control_plane_version
        self.run_control_plane = self.control_plane_root / "run_control_plane.py"
        self.attestation_helper = (
            self.source_root
            / "tst/publication/capture_frontier_pre_policy_promotion_attestation.py"
        )
        self._validate_cli_identities()
        sys.path.insert(0, str(self.source_root))
        from tst.publication import (  # noqa: PLC0415
            q019_excluded_pilot_launch_policy_preparation_v1 as preparation,
        )
        self.preparation = preparation
        self.admission = registered_admission
        self.final = self._read_bound_json(
            self.final_bindings_path,
            expected_sha256=self.final_bindings_sha256,
            label="Q019 pilot final bindings",
        )
        self.reviewed_policy = self._read_bound_json(
            self.reviewed_policy_path,
            expected_sha256=self.reviewed_policy_sha256,
            label="Q019 pilot reviewed policy",
        )
        self.preparation.validate_final_binding_files(self.final)
        _require(
            self.final["source_commit"] == self.source_commit
            and self.final["installed_control_plane_version"]
            == self.control_plane_version,
            "Q019 pilot final bindings differ from driver identities",
        )
        self.members = self.preparation._pilot_members()
        self.campaign = self.preparation.CAMPAIGN
        self.work_root = PIC_ROOT / "jobs" / self.campaign / "excluded_pilots_v1"
        self.admission_root = (
            PIC_ROOT / "analysis" / self.campaign / "excluded_pilot_admissions_v1"
        )
        self.index_path = (
            PIC_ROOT
            / "analysis"
            / self.campaign
            / "q019_excluded_pilot_execution_index_v1.json"
        )

    def _validate_cli_identities(self) -> None:
        for label, value, pattern in (
            ("source commit", self.source_commit, COMMIT_PATTERN),
            ("control-plane version", self.control_plane_version, SHA256_PATTERN),
            ("final bindings SHA-256", self.final_bindings_sha256, SHA256_PATTERN),
            ("reviewed policy SHA-256", self.reviewed_policy_sha256, SHA256_PATTERN),
            ("baseline policy SHA-256", self.baseline_policy_sha256, SHA256_PATTERN),
            (
                "baseline promotion SHA-256",
                self.baseline_promotion_sha256,
                SHA256_PATTERN,
            ),
        ):
            _require(pattern.fullmatch(value) is not None, f"malformed Q019 pilot {label}")
        head = _run(["/usr/bin/git", "-C", self.source_root, "rev-parse", "HEAD"])
        status = _run(
            [
                "/usr/bin/git",
                "-C",
                self.source_root,
                "status",
                "--porcelain=v1",
                "--untracked-files=all",
            ]
        )
        _require(
            head == self.source_commit and not status,
            "Q019 pilot source is not the exact clean selected commit",
        )
        _require(
            self.run_control_plane.is_file(),
            "Q019 pilot installed control-plane runner is absent",
        )

    @staticmethod
    def _read_bound_json(
        path: Path, *, expected_sha256: str, label: str
    ) -> dict[str, object]:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
        try:
            before = os.fstat(descriptor)
            _require(
                stat.S_ISREG(before.st_mode)
                and not before.st_mode & 0o222
                and before.st_nlink == 1,
                f"{label} is not one read-only regular file",
            )
            chunks = []
            while chunk := os.read(descriptor, 1024 * 1024):
                chunks.append(chunk)
            payload = b"".join(chunks)
            after = os.fstat(descriptor)
            _require(
                (
                    before.st_dev,
                    before.st_ino,
                    before.st_mode,
                    before.st_nlink,
                    before.st_size,
                    before.st_mtime_ns,
                    before.st_ctime_ns,
                )
                == (
                    after.st_dev,
                    after.st_ino,
                    after.st_mode,
                    after.st_nlink,
                    after.st_size,
                    after.st_mtime_ns,
                    after.st_ctime_ns,
                )
                and len(payload) == after.st_size
                and hashlib.sha256(payload).hexdigest() == expected_sha256,
                f"{label} changed or its digest drifted",
            )
        finally:
            os.close(descriptor)
        value = json.loads(payload)
        _require(type(value) is dict, f"{label} must be one JSON object")
        return value

    def _controller(self, entrypoint: str, *arguments: str | Path) -> str:
        return _run(
            [PYTHON, "-I", "-B", self.run_control_plane, entrypoint, *arguments]
        )

    @staticmethod
    def _queue() -> str:
        return _run(
            [
                "/usr/bin/squeue",
                "--clusters=frontier",
                "-u",
                os.environ.get("USER", "dfielding"),
                "-h",
                "-o",
                "%i|%a|%P|%q|%T|%j|%k",
            ]
        )

    def _wait_for_empty_queue(self) -> None:
        while True:
            queue = self._queue()
            if not queue:
                return
            print(f"{_utc(_utc_now())} waiting for empty queue:\n{queue}", flush=True)
            time.sleep(30)

    def _capture_attestation(self, authorization_id: str, phase: str) -> Path:
        staging = Path(
            _run(
                [
                    PYTHON,
                    "-B",
                    self.attestation_helper,
                    "capture",
                    "--authorization-id",
                    authorization_id,
                    "--control-plane-version",
                    self.control_plane_version,
                    "--phase",
                    phase,
                ],
                cwd=self.source_root,
            )
        )
        sealed = Path(
            _run(
                [
                    PYTHON,
                    "-B",
                    self.attestation_helper,
                    "seal",
                    "--staging-dir",
                    staging,
                    "--attest-reviewed",
                ],
                cwd=self.source_root,
            )
        )
        return sealed / "attestation.json"

    @staticmethod
    def _active_policy() -> tuple[dict[str, object], str, str]:
        policy_path = PIC_ROOT / "policy/storage_policy.json"
        promotion_path = PIC_ROOT / "policy/active_promotion.json"
        return (
            json.loads(policy_path.read_text(encoding="utf-8")),
            _sha256(policy_path),
            _sha256(promotion_path),
        )

    def _promote_policy_if_needed(self) -> None:
        policy, policy_sha256, promotion_sha256 = self._active_policy()
        if policy_sha256 == self.reviewed_policy_sha256:
            _require(
                policy == self.reviewed_policy
                and len(policy["registered_science_slices"]) == len(self.members),
                "active Q019 pilot policy differs from reviewed policy",
            )
            return
        _require(
            policy_sha256 == self.baseline_policy_sha256
            and promotion_sha256 == self.baseline_promotion_sha256
            and policy["registered_science_slices"] == [],
            "active policy is neither the exact Q019 pilot predecessor nor successor",
        )
        expected = self.preparation.materialize_q019_promotable_policy(
            baseline_policy=policy, final_bindings=self.final
        )
        _require(
            expected == self.reviewed_policy,
            "reviewed Q019 pilot policy is not the exact materialized successor",
        )
        self._wait_for_empty_queue()
        authorization_id = f"q019-pilot-policy-{self.source_commit[:8]}-v1"
        attestation = self._capture_attestation(
            authorization_id, "pre_policy_promotion"
        )
        policy, policy_sha256, promotion_sha256 = self._active_policy()
        _require(
            policy_sha256 == self.baseline_policy_sha256
            and promotion_sha256 == self.baseline_promotion_sha256
            and policy["registered_science_slices"] == [],
            "active Q019 pilot predecessor changed during promotion preparation",
        )
        self._controller(
            "promote_active_policy.py",
            "--reviewed-policy",
            self.reviewed_policy_path,
            "--pre-policy-promotion-attestation",
            attestation,
            "--pre-policy-promotion-authorization-id",
            authorization_id,
        )
        policy, policy_sha256, _ = self._active_policy()
        _require(
            policy_sha256 == self.reviewed_policy_sha256
            and len(policy["registered_science_slices"]) == len(self.members),
            "Q019 pilot policy promotion produced the wrong successor",
        )

    def _state_path(self, artifact_id: str) -> Path:
        _require(
            artifact_id in EXPECTED_ARTIFACTS,
            "Q019 pilot state key is not an exact artifact ID",
        )
        return self.work_root / "state" / f"{artifact_id}.json"

    def _load_state(self, artifact_id: str) -> dict[str, object] | None:
        path = self._state_path(artifact_id)
        if not path.exists():
            return None
        value = json.loads(path.read_text(encoding="utf-8"))
        _require(
            type(value) is dict and value.get("artifact_id") == artifact_id,
            f"{artifact_id}: Q019 pilot state identity drifted",
        )
        return value

    def _save_state(self, artifact_id: str, value: dict[str, object]) -> None:
        _replace_json(self._state_path(artifact_id), value)

    @staticmethod
    def _ledger_job_for_submission(submission_id: str) -> str | None:
        matches: list[str] = []
        with LEDGER_JSONL.open("r", encoding="utf-8") as stream:
            for line in stream:
                record = json.loads(line)
                if (
                    record.get("submission_id") == submission_id
                    and isinstance(record.get("job_id"), str)
                ):
                    matches.append(record["job_id"])
        unique = sorted(set(matches))
        _require(
            len(unique) <= 1,
            f"submission {submission_id} has multiple Q019 pilot job IDs",
        )
        return unique[0] if unique else None

    @staticmethod
    def _wait_for_terminal(job_id: str) -> None:
        while True:
            queue = _run(
                [
                    "/usr/bin/squeue",
                    "--clusters=frontier",
                    "-j",
                    job_id,
                    "-h",
                    "-o",
                    "%i|%T",
                ]
            )
            if queue:
                print(f"{_utc(_utc_now())} waiting for Q019 pilot job {queue}", flush=True)
                time.sleep(15)
                continue
            accounting = _run(
                [
                    "/usr/bin/sacct",
                    "--clusters=frontier",
                    "-X",
                    "-j",
                    job_id,
                    "--format=JobIDRaw,State,ExitCode",
                    "-n",
                    "-P",
                ]
            )
            rows = [line.split("|") for line in accounting.splitlines() if line]
            root = [row for row in rows if row[0] == job_id]
            if len(root) == 1 and root[0][1]:
                return
            time.sleep(10)

    def _materialize_attempt_inputs(
        self,
        *,
        artifact_id: str,
        submission_id: str,
        authorization_id: str,
    ) -> Path:
        root = self.work_root / "handoffs" / artifact_id / submission_id
        root.mkdir(parents=True, exist_ok=False)
        now = _utc_now()
        timeout = self.preparation.materialize_q019_timeout_margin(
            artifact_id=artifact_id,
            final_bindings=self.final,
            measured_utc=_utc(now - timedelta(minutes=1)),
            expires_utc=_utc(now + timedelta(hours=1)),
        )
        timeout_path = root / "timeout_margin.json"
        _atomic_json(timeout_path, timeout)
        queue_path = root / "queue_snapshot.txt"
        _require(not self._queue(), "Q019 pilot pre-manifest queue is not empty")
        _write_empty_read_only(queue_path)
        pre_manifest = self._capture_attestation(authorization_id, "pre_manifest")
        config = self.preparation.materialize_q019_pre_submit_config(
            artifact_id=artifact_id,
            submission_id=submission_id,
            final_bindings=self.final,
            pre_manifest_attestation=pre_manifest,
            timeout_margin_artifact=timeout_path,
            queue_snapshot=queue_path,
            site_policy_checked_utc=_utc(now),
            now=now,
        )
        config_path = root / "pre_submit_config.json"
        _atomic_json(config_path, config)
        return Path(
            self._controller("create_pre_submit_manifest.py", "--config", config_path)
        )

    @staticmethod
    def _reconciliation_event(
        *,
        submission_id: str,
        job_id: str,
        expected_event_sha256: str,
    ) -> dict[str, object]:
        matches = []
        with LEDGER_JSONL.open("r", encoding="utf-8") as stream:
            for line in stream:
                record = json.loads(line)
                if (
                    record.get("event_type") == "reconciliation"
                    and record.get("submission_id") == submission_id
                    and record.get("job_id") == job_id
                    and record.get("event_sha256") == expected_event_sha256
                ):
                    matches.append(record)
        _require(
            len(matches) == 1,
            "Q019 pilot lacks one exact authenticated reconciliation event",
        )
        return matches[0]

    def _attempt_record(
        self,
        *,
        index: int,
        member: Mapping[str, object],
        state: Mapping[str, object],
        admission_path: Path,
    ) -> dict[str, object]:
        artifact_root = (
            PIC_ROOT
            / self.preparation.RUN_NAMESPACE
            / str(state["submission_id"])
        )
        receipt_path = (
            artifact_root / "analysis" / self.admission.EXECUTION_RECEIPT_NAME
        )
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        admission_record = self.admission.validate_admission(
            json.loads(admission_path.read_text(encoding="utf-8"))
        )
        event = self._reconciliation_event(
            submission_id=str(state["submission_id"]),
            job_id=str(state["job_id"]),
            expected_event_sha256=str(receipt["reconciliation_event_sha256"]),
        )
        artifact_inventory = dict(receipt["artifact_inventory"])
        return {
            "attempt_index": index,
            "artifact_id": member["artifact_id"],
            "source_case_id": member["source_case_id"],
            "authorization_id": state["authorization_id"],
            "submission_id": state["submission_id"],
            "job_id": state["job_id"],
            "artifact_root": str(artifact_root),
            "admission": _binding(admission_path),
            "admission_facts": {
                key: admission_record[key]
                for key in (
                    "record_type",
                    "case_id",
                    "artifact_root",
                    "registered_execution_identity",
                    "execution_profile",
                    "runtime_completion",
                    "raw_science_admission_eligible",
                    "problem_reported_saturation_evidence_eligible",
                    "saturation_evidence_eligible",
                    "authorization",
                )
            },
            "execution_provenance": _execution_provenance(admission_record),
            "artifact_inventory": artifact_inventory,
            "artifact_payload_bytes": _artifact_payload_bytes(
                artifact_inventory, artifact_root=artifact_root
            ),
            "reconciliation_event": event,
        }

    @staticmethod
    def _validate_admitted_resume_state(
        *,
        index: int,
        member: Mapping[str, object],
        state: Mapping[str, object],
        admission_path: Path,
        admission_record: Mapping[str, object],
    ) -> None:
        identity = admission_record.get("registered_execution_identity")
        _require(
            type(identity) is dict
            and state.get("status") == "admitted"
            and state.get("attempt_index") == index
            and state.get("artifact_id") == member["artifact_id"]
            and state.get("source_case_id") == member["source_case_id"]
            and state.get("authorization_id")
            == identity.get("registered_science_authorization_id")
            and state.get("submission_id") == identity.get("submission_id")
            and state.get("job_id") == identity.get("slurm_job_id")
            and state.get("admission_path") == str(admission_path)
            and state.get("admission_sha256") == _sha256(admission_path),
            f"{member['artifact_id']}: admitted pilot resume state drifted",
        )

    def _run_attempt(
        self, index: int, member: Mapping[str, object]
    ) -> dict[str, object]:
        artifact_id = str(member["artifact_id"])
        admission_path = self.admission_root / f"{artifact_id}.json"
        state = self._load_state(artifact_id)
        if admission_path.exists():
            _require(state is not None, f"{artifact_id}: admitted pilot lacks resume state")
            record = json.loads(admission_path.read_text(encoding="utf-8"))
            record = self.admission.validate_admission(record)
            self._validate_admitted_resume_state(
                index=index,
                member=member,
                state=state,
                admission_path=admission_path,
                admission_record=record,
            )
            return self._attempt_record(
                index=index,
                member=member,
                state=state,
                admission_path=admission_path,
            )
        if state is None:
            self._wait_for_empty_queue()
            submission_id = str(uuid.uuid4())
            authorization_id = self.preparation._authorization_id(index, member)
            state = {
                "artifact_id": artifact_id,
                "source_case_id": member["source_case_id"],
                "attempt_index": index,
                "submission_id": submission_id,
                "authorization_id": authorization_id,
                "status": "materializing",
            }
            self._save_state(artifact_id, state)
            manifest = self._materialize_attempt_inputs(
                artifact_id=artifact_id,
                submission_id=submission_id,
                authorization_id=authorization_id,
            )
            state["manifest_path"] = str(manifest)
            state["status"] = "manifest_created"
            self._save_state(artifact_id, state)
        job_id = str(state.get("job_id", ""))
        if not job_id:
            recovered = self._ledger_job_for_submission(str(state["submission_id"]))
            if recovered:
                job_id = recovered
                state["job_id"] = job_id
                state["status"] = "submitted"
                self._save_state(artifact_id, state)
            else:
                _require(
                    state.get("status") == "manifest_created",
                    f"{artifact_id}: incomplete pre-submission state requires review",
                )
                self._wait_for_empty_queue()
                pre_submit = self._capture_attestation(
                    str(state["authorization_id"]), "pre_submit_wrapper"
                )
                output = _run(
                    [
                        self.control_plane_root / "submit_frontier_job.sh",
                        Path(str(state["manifest_path"])),
                        pre_submit,
                    ]
                )
                match = JOB_ID_PATTERN.search(output)
                _require(
                    match is not None,
                    f"{artifact_id}: could not parse submission output: {output}",
                )
                job_id = match.group(1)
                state["job_id"] = job_id
                state["status"] = "submitted"
                self._save_state(artifact_id, state)
                print(f"{artifact_id}: {output}", flush=True)
        self._wait_for_terminal(job_id)
        evidence = json.loads(
            self._controller(
                "reconcile_q019_registered_execution.py",
                "--job-id",
                job_id,
                "--ledger-jsonl",
                LEDGER_JSONL,
                "--ledger-csv",
                LEDGER_CSV,
                "--receipts-jsonl",
                RECEIPTS_JSONL,
                "--mirror-jsonl",
                MIRROR_JSONL,
            )
        )
        _require(type(evidence) is dict, f"{artifact_id}: malformed reconciler evidence")
        artifact_root = (
            PIC_ROOT
            / self.preparation.RUN_NAMESPACE
            / str(state["submission_id"])
        )
        admission, _ = self.admission.derive_case_bundle(
            case_id=str(member["source_case_id"]),
            artifact_root=artifact_root,
            q043_qualification_path=Path(
                str(self.final["q043_registered_matrix_path"])
            ),
            q043_artifact_root=PIC_ROOT,
            q023_qualification_path=Path(
                str(self.final["q023_registered_matrix_path"])
            ),
        )
        _require(
            admission["execution_profile"]["artifact_id"] == artifact_id
            and admission["runtime_completion"]["stop_reason_code"] == "1903"
            and admission["runtime_completion"][
                "runtime_controller_trigger_cycle"
            ]
            == 20
            and admission["saturation_evidence_eligible"] is False,
            f"{artifact_id}: admitted pilot stop or exclusion contract drifted",
        )
        _atomic_json(admission_path, admission)
        state["status"] = "admitted"
        state["admission_path"] = str(admission_path)
        state["admission_sha256"] = _sha256(admission_path)
        self._save_state(artifact_id, state)
        print(f"{artifact_id}: admitted {state['admission_sha256']}", flush=True)
        self._wait_for_empty_queue()
        return self._attempt_record(
            index=index,
            member=member,
            state=state,
            admission_path=admission_path,
        )

    def run(self) -> int:
        self._promote_policy_if_needed()
        self.admission_root.mkdir(parents=True, exist_ok=True)
        attempts = [
            self._run_attempt(index, member)
            for index, member in enumerate(self.members, 1)
        ]
        execution_index = build_execution_index(attempts)
        if self.index_path.exists():
            existing = json.loads(self.index_path.read_text(encoding="utf-8"))
            validate_execution_index_files(existing)
            _require(
                existing == execution_index,
                "existing Q019 pilot execution index differs",
            )
        else:
            _atomic_json(self.index_path, execution_index)
        print(
            json.dumps(
                {
                    "attempt_count": len(attempts),
                    "execution_index_path": str(self.index_path),
                    "execution_index_sha256": _sha256(self.index_path),
                    "total_consumed_node_hours": execution_index[
                        "total_consumed_node_hours"
                    ],
                    "saturation_evidence_eligible": False,
                },
                sort_keys=True,
            ),
            flush=True,
        )
        return 0


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--control-plane-version", required=True)
    parser.add_argument("--final-bindings", required=True, type=Path)
    parser.add_argument("--final-bindings-sha256", required=True)
    parser.add_argument("--reviewed-policy", required=True, type=Path)
    parser.add_argument("--reviewed-policy-sha256", required=True)
    parser.add_argument("--baseline-policy-sha256", required=True)
    parser.add_argument("--baseline-promotion-sha256", required=True)
    return parser


def main() -> int:
    return Driver(_parser().parse_args()).run()


if __name__ == "__main__":
    raise SystemExit(main())
