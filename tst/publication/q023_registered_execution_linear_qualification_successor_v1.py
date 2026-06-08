#!/usr/bin/env python3
"""Qualify the exact registered Q023 linear Bell matrix without granting authority."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import stat
from typing import Mapping, Sequence

from tst.publication import analyze_q023_paper_bell_linear_joverc as bell
from tst.publication import q023_registered_launch_policy_preparation_successor_v1 as prep


SCHEMA_VERSION = 1
SUCCESSOR_ID = "q023_registered_execution_linear_qualification_successor_v1"
CASE_RECORD_TYPE = "q023_registered_execution_linear_case_admission"
MATRIX_RECORD_TYPE = "q023_registered_execution_linear_matrix_qualification"
QUALIFICATION_EFFECT = (
    "registered_q023_linear_predecessor_only_no_launch_no_policy_no_q019_"
    "no_science_no_publication_authority"
)
AUTHORIZATION_BOUNDARY = {
    "launch_authorized": False,
    "scheduler_submission_authorized": False,
    "policy_mutation_authorized": False,
    "frontier_execution_authorized": False,
    "q019_qualification_authorized": False,
    "scientific_claim_authorized": False,
    "publication_authorized": False,
}
_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH
_CANONICAL_PROJECT_HOME_LEDGER_ROOT = Path(
    "/ccs/proj/ast207/proj-shared/PIC"
)


class QualificationError(ValueError):
    """Reject incomplete, substituted, mutable, or noncanonical registered evidence."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise QualificationError(message)


def _canonical_bytes(value: object) -> bytes:
    try:
        return (
            json.dumps(
                value,
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise QualificationError("Q023 qualification contains noncanonical JSON") from error


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


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


def _stable_read_only_bytes(path: Path, *, root: Path, label: str) -> tuple[Path, bytes]:
    root = Path(os.path.abspath(root)).resolve(strict=True)
    lexical = Path(os.path.abspath(path))
    try:
        resolved = lexical.resolve(strict=True)
    except OSError as error:
        raise QualificationError(f"{label} does not resolve to one file") from error
    _require(resolved == lexical, f"{label} path contains a symlink or alias")
    try:
        lexical.relative_to(root)
    except ValueError as error:
        raise QualificationError(f"{label} is outside the authorized root") from error
    try:
        descriptor = os.open(lexical, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    except OSError as error:
        raise QualificationError(f"{label} is not an openable regular file") from error
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and not before.st_mode & _WRITE_BITS,
            f"{label} is not one read-only regular file",
        )
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        current = os.stat(lexical, follow_symlinks=False)
        identity = lambda value: (
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_nlink,
            value.st_size,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
        _require(
            identity(before) == identity(after)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and len(payload) == after.st_size,
            f"{label} changed while reading",
        )
        return lexical, bytes(payload)
    finally:
        os.close(descriptor)


def _json_object(payload: bytes, *, label: str) -> dict[str, object]:
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise QualificationError(f"{label} is not UTF-8 JSON") from error
    _require(type(value) is dict, f"{label} is not an object")
    return value


def _expected_authorization_id(member_id: str) -> str:
    members = list(bell._manifest_members())
    _require(member_id in members, "Q023 case admission member is unknown")
    return f"q023-linear-{members.index(member_id) + 1:03d}-v1"


def _trusted_reconciliation_and_rederivation(
    *,
    receipt: Mapping[str, object],
    receipt_path: Path,
    receipt_payload: bytes,
    terminal_path: Path,
    terminal_payload: bytes,
    receipt_mirror_path: Path,
    terminal_mirror_path: Path,
    authorized_orion_root: Path,
    authorized_project_home_root: Path,
) -> dict[str, object]:
    """Rebuild exact Q023 receipts from the canonical ledger and installed producer."""
    _require(
        Path(authorized_orion_root).resolve(strict=True)
        == prep.AUTHORIZED_ORION_ROOT.resolve(strict=True)
        and Path(authorized_project_home_root).resolve(strict=True)
        == prep.CANONICAL_PROJECT_HOME_ROOT.resolve(strict=True),
        "Q023 installed-producer verification requires canonical production roots",
    )
    ledger_root = _CANONICAL_PROJECT_HOME_LEDGER_ROOT
    _require(
        ledger_root.resolve(strict=True)
        == prep.CANONICAL_PROJECT_HOME_ROOT.resolve(strict=True),
        "Q023 canonical Project Home ledger alias drifted",
    )
    version = receipt.get("control_plane_version")
    _require(
        type(version) is str and bell._SHA256.fullmatch(version) is not None,
        "Q023 receipt control-plane version is malformed",
    )
    filenames = ("ledger.py", "reconcile_q023_registered_execution.py")
    try:
        with bell.q043_registered._installed_control_plane_modules(
            version, filenames
        ) as (modules, pair):
            ledger = modules["ledger.py"]
            producer = modules["reconcile_q023_registered_execution.py"]
            ledger_jsonl = authorized_orion_root / "ledger/node_hours.jsonl"
            receipts_jsonl = authorized_orion_root / "ledger/mirror_receipts.jsonl"
            mirror_jsonl = ledger_root / "ledger/node_hours.jsonl"
            with ledger.validated_read_only_mirrored_state_snapshot(
                ledger_jsonl,
                receipts_jsonl,
                mirror_jsonl,
                ledger_root=authorized_orion_root,
                receipts_root=authorized_orion_root,
                mirror_root=ledger_root,
            ) as records:
                ledger.require_explicit_genesis(records)
                mirror_receipts = ledger.validate_receipts(
                    receipts_jsonl,
                    records,
                    mirror_jsonl=mirror_jsonl,
                    mirror_transport="filesystem_copy",
                    root=authorized_orion_root,
                )
                events = [dict(item) for item in records]
                acknowledgments = [dict(item) for item in mirror_receipts]
            matches = [
                item
                for item in events
                if item.get("event_type") == "reconciliation"
                and item.get("event_sha256")
                == receipt.get("reconciliation_event_sha256")
            ]
            _require(
                len(matches) == 1,
                "Q023 receipt lacks one canonical mirrored-ledger reconciliation",
            )
            event = matches[0]
            mirror_matches = [
                item
                for item in acknowledgments
                if item.get("mirrored_event_sha256") == event["event_sha256"]
            ]
            _require(
                len(mirror_matches) == 1,
                "Q023 receipt lacks one canonical Project Home mirror acknowledgment",
            )
            mirror_ack = mirror_matches[0]
            installed_producer_sha256 = pair["digests"].get(
                "reconcile_q023_registered_execution.py"
            )
            installed_trampoline_sha256 = pair["digests"].get(
                "launch_trampoline.py"
            )
            _require(
                receipt.get("producer")
                == {
                    "entrypoint": "reconcile_q023_registered_execution.py",
                    "entrypoint_sha256": installed_producer_sha256,
                    "launch_trampoline_sha256": installed_trampoline_sha256,
                    "control_plane_version": version,
                },
                "Q023 receipt producer differs from the paired installed generation",
            )
            _require(
                mirror_ack.get("mirror_ack_sha256")
                == receipt.get("reconciliation_mirror_ack_sha256")
                and mirror_ack.get("mirror_transport") == "filesystem_copy"
                and mirror_ack.get("mirror_destination") == str(mirror_jsonl),
                "Q023 receipt mirror acknowledgment differs from canonical ledger evidence",
            )
            _require(
                event.get("submission_scope") == "registered_science"
                and event.get("campaign") == prep.CAMPAIGN
                and event.get("test_id") == receipt.get("member_id")
                and event.get("reservation_id") == receipt.get("reservation_id")
                and event.get("submission_id") == receipt.get("submission_id")
                and event.get("job_id") == receipt.get("slurm_job_id")
                and event.get("artifact_dir") == receipt.get("artifact_dir")
                and event.get("git_commit") == receipt.get("source_commit")
                and event.get("clean_candidate_manifest_sha256")
                == receipt.get("clean_candidate_manifest_sha256")
                and event.get("executable_sha256")
                == receipt.get("executable_sha256")
                and event.get("manifest_path")
                == receipt.get("pre_submit_manifest_path")
                and event.get("manifest_sha256")
                == receipt.get("pre_submit_manifest_sha256")
                and event.get("control_plane_version") == version
                and event.get("reconciled_by_control_plane_version") == version
                and event.get("registered_science_authorization_id")
                == receipt.get("registered_science_authorization_id")
                and event.get("reconciled") is True
                and event.get("state") == "COMPLETED"
                and event.get("scheduler_exit_code") == "0:0",
                "Q023 receipt differs from canonical reconciliation evidence",
            )
            derived = producer.derive_q023_registered_execution_evidence(
                event,
                mirror_ack,
                pair["inventory"],
                authorized_pic_root=authorized_orion_root,
                authorized_project_home_root=authorized_project_home_root,
            )
    except (
        AttributeError,
        KeyError,
        OSError,
        TypeError,
        ValueError,
        bell.q043_registered.AdmissionError,
    ) as error:
        raise QualificationError(
            "Q023 installed producer could not re-derive registered evidence"
        ) from error
    expected = (
        terminal_path,
        terminal_payload,
        receipt_path,
        receipt_payload,
        terminal_mirror_path,
        receipt_mirror_path,
    )
    _require(
        derived == expected,
        "installed Q023 producer re-derived different evidence bytes or paths",
    )
    return {
        "control_plane_version": version,
        "entrypoint": "reconcile_q023_registered_execution.py",
        "entrypoint_sha256": receipt["producer"]["entrypoint_sha256"],
        "reconciliation_event_sha256": event["event_sha256"],
        "reconciliation_mirror_ack_sha256": mirror_ack["mirror_ack_sha256"],
        "receipt_sha256": hashlib.sha256(receipt_payload).hexdigest(),
        "terminal_receipt_sha256": hashlib.sha256(terminal_payload).hexdigest(),
        "exact_byte_rederivation_passed": True,
    }


def _receipt_and_terminal_evidence(
    *,
    artifact_root: Path,
    authorized_orion_root: Path,
    authorized_project_home_root: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    relative_receipt = "analysis/q023_registered_execution_receipt.json"
    receipt_path, receipt_payload = _stable_read_only_bytes(
        artifact_root / relative_receipt,
        root=authorized_orion_root,
        label="Q023 Orion registered execution receipt",
    )
    receipt = _json_object(receipt_payload, label="Q023 registered execution receipt")
    mirrors = receipt.get("project_home_mirrors")
    _require(
        type(mirrors) is dict
        and set(mirrors)
        == {"registered_execution_receipt_path", "terminal_receipt_path"},
        "Q023 Project Home mirror schema drifted",
    )
    submission_id = receipt.get("submission_id")
    expected_mirror_root = (
        authorized_project_home_root
        / prep.PROJECT_HOME_RECEIPT_NAMESPACE
        / str(submission_id)
    )
    expected_receipt_mirror = (
        expected_mirror_root / "q023_registered_execution_receipt.json"
    )
    expected_terminal_mirror = expected_mirror_root / "q023_terminal_receipt.json"
    _require(
        mirrors["registered_execution_receipt_path"] == str(expected_receipt_mirror)
        and mirrors["terminal_receipt_path"] == str(expected_terminal_mirror),
        "Q023 Project Home mirror paths drifted",
    )
    _, receipt_mirror_payload = _stable_read_only_bytes(
        expected_receipt_mirror,
        root=authorized_project_home_root,
        label="Q023 Project Home registered execution receipt",
    )
    _require(
        receipt_mirror_payload == receipt_payload,
        "Q023 Orion and Project Home execution receipts differ",
    )
    terminal_path = artifact_root / "analysis/q023_terminal_receipt.json"
    _, terminal_payload = _stable_read_only_bytes(
        terminal_path,
        root=authorized_orion_root,
        label="Q023 Orion terminal receipt",
    )
    _, terminal_mirror_payload = _stable_read_only_bytes(
        expected_terminal_mirror,
        root=authorized_project_home_root,
        label="Q023 Project Home terminal receipt",
    )
    _require(
        terminal_payload == terminal_mirror_payload
        and hashlib.sha256(terminal_payload).hexdigest()
        == receipt.get("terminal_receipt_sha256"),
        "Q023 paired terminal receipt bytes or digest drifted",
    )
    terminal = _json_object(terminal_payload, label="Q023 terminal receipt")
    _require(
        terminal.get("record_type")
        == "q023_registered_execution_terminal_receipt"
        and terminal.get("member_id") == receipt.get("member_id")
        and terminal.get("submission_id") == submission_id
        and terminal.get("slurm_job_id") == receipt.get("slurm_job_id")
        and terminal.get("slurm_terminal_state") == "COMPLETED"
        and terminal.get("slurm_exit_code") == "0:0"
        and terminal.get("terminal_time") == bell.LINEAR_RUNTIME_TLIM
        and terminal.get("raw_inventory_sha256")
        == receipt.get("raw_inventory_sha256"),
        "Q023 terminal and execution receipts do not cross-link",
    )
    installed_rederivation = _trusted_reconciliation_and_rederivation(
        receipt=receipt,
        receipt_path=receipt_path,
        receipt_payload=receipt_payload,
        terminal_path=terminal_path,
        terminal_payload=terminal_payload,
        receipt_mirror_path=expected_receipt_mirror,
        terminal_mirror_path=expected_terminal_mirror,
        authorized_orion_root=authorized_orion_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    return receipt, {
        "registered_execution_receipt": {
            "orion_path": str(receipt_path),
            "project_home_path": str(expected_receipt_mirror),
            "sha256": hashlib.sha256(receipt_payload).hexdigest(),
            "byte_count": len(receipt_payload),
        },
        "terminal_receipt": {
            "orion_path": str(terminal_path),
            "project_home_path": str(expected_terminal_mirror),
            "sha256": hashlib.sha256(terminal_payload).hexdigest(),
            "byte_count": len(terminal_payload),
        },
        "installed_reconciliation_rederivation": installed_rederivation,
    }


def _validated_case_root(
    artifact_root: Path, *, authorized_orion_root: Path
) -> Path:
    authorized_orion_root = Path(os.path.abspath(authorized_orion_root))
    artifact_root = Path(os.path.abspath(artifact_root))
    _require(
        artifact_root.resolve(strict=True) == artifact_root,
        "Q023 case root contains a symlink or alias",
    )
    expected_run_root = authorized_orion_root / prep.RUN_NAMESPACE
    try:
        relative_case_root = artifact_root.relative_to(expected_run_root)
    except ValueError as error:
        raise QualificationError(
            "Q023 case root is outside the registered run namespace"
        ) from error
    _require(
        len(relative_case_root.parts) == 1
        and artifact_root.is_dir()
        and not artifact_root.stat(follow_symlinks=False).st_mode & _WRITE_BITS,
        "Q023 case root is not one sealed registered submission directory",
    )
    return artifact_root


def _validate_output_index_metadata(
    values: list[tuple[int, list[int], list[float]]],
    *,
    terminal_cycle: object,
    terminal_time: object,
) -> None:
    previous_cycle: int | None = None
    previous_time: float | None = None
    for output_index, cycles, times in values:
        _require(
            len(cycles) == len(prep.REQUIRED_OUTPUT_FIELDS)
            and len(times) == len(prep.REQUIRED_OUTPUT_FIELDS)
            and len(set(cycles)) == 1
            and len(set(times)) == 1,
            f"Q023 output index {output_index} has inconsistent embedded metadata",
        )
        cycle = cycles[0]
        time = times[0]
        if output_index == prep.REQUIRED_OUTPUT_INDICES[0]:
            _require(
                cycle == 0 and time == 0.0,
                "Q023 output-index metadata does not begin at cycle/time zero",
            )
        else:
            _require(
                previous_cycle is not None
                and previous_time is not None
                and cycle > previous_cycle
                and time > previous_time,
                "Q023 output-index metadata is not strictly chronological",
            )
        previous_cycle = cycle
        previous_time = time
    _require(
        previous_cycle == terminal_cycle
        and type(terminal_time) in (int, float)
        and previous_time is not None
        and math.isclose(
            previous_time,
            float(terminal_time),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "Q023 final output-index metadata differs from the terminal receipt",
    )


def _derive_registered_case_evidence(
    *,
    member_id: str,
    artifact_root: Path,
    q043_registered_raw_oracle_dependency: Mapping[str, object],
    q043_artifact_root: Path,
    authorized_orion_root: Path,
    authorized_project_home_root: Path,
) -> tuple[
    dict[str, object],
    dict[str, object],
    dict[str, object],
    dict[str, object],
]:
    artifact_root = _validated_case_root(
        artifact_root, authorized_orion_root=authorized_orion_root
    )
    members = bell._manifest_members()
    _require(member_id in members, "Q023 retained case member is unknown")
    member = members[member_id]
    receipt, paired_evidence = _receipt_and_terminal_evidence(
        artifact_root=artifact_root,
        authorized_orion_root=authorized_orion_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    _require(
        receipt.get("member_id") == member_id,
        "Q023 reconciled receipt member differs from the requested case",
    )
    inventory = receipt.get("raw_inventory")
    _require(
        type(inventory) is list
        and len(inventory)
        == len(prep.REQUIRED_OUTPUT_INDICES) * len(prep.REQUIRED_OUTPUT_FIELDS),
        "Q023 reconciled receipt lacks the exact 445-member raw inventory",
    )
    retained_bindings: list[dict[str, object]] = []
    inventory_position = 0
    for output_index in prep.REQUIRED_OUTPUT_INDICES:
        for variable in prep.REQUIRED_OUTPUT_FIELDS:
            item = inventory[inventory_position]
            inventory_position += 1
            expected_path = prep._raw_relative_path(
                member_id, variable, output_index
            )
            _require(
                type(item) is dict
                and set(item)
                == {
                    "path",
                    "sha256",
                    "byte_count",
                    "member_id",
                    "variable",
                    "output_index",
                }
                and item["path"] == expected_path
                and item["member_id"] == member_id
                and item["variable"] == variable
                and type(item["output_index"]) is int
                and item["output_index"] == output_index
                and type(item["byte_count"]) is int
                and item["byte_count"] > 0
                and type(item["sha256"]) is str,
                "Q023 reconciled raw inventory ordering or identity drifted",
            )
            retained_bindings.append(
                {
                    "path": "raw/" + expected_path,
                    "sha256": item["sha256"],
                    "byte_count": item["byte_count"],
                }
            )
    try:
        retained_raw_batch = bell._RetainedRawBatch(
            artifact_root, retained_bindings
        )
    except bell.ContractError as error:
        raise QualificationError(
            "Q023 could not retain the exact reconciled raw inventory"
        ) from error
    try:
        raw_artifacts: list[dict[str, object]] = []
        mhd_datasets = []
        output_index_metadata: list[tuple[int, list[int], list[float]]] = []
        inventory_position = 0
        for output_index in prep.REQUIRED_OUTPUT_INDICES:
            output_cycles: list[int] = []
            output_times: list[float] = []
            for variable in prep.REQUIRED_OUTPUT_FIELDS:
                item = inventory[inventory_position]
                inventory_position += 1
                relative_path = "raw/" + str(item["path"])
                payload = retained_raw_batch.payload(relative_path)
                raw_path = retained_raw_batch.absolute_path(relative_path)
                try:
                    dataset = bell.binary.parse_athenak_binary_bytes(
                        payload, source=str(raw_path)
                    )
                except bell.binary.AnalysisError as error:
                    raise QualificationError(
                        "Q023 retained raw output failed strict AthenaK parsing"
                    ) from error
                raw_artifacts.append(
                    {
                        "path": relative_path,
                        "sha256": item["sha256"],
                        "variable": variable,
                        "cycle": dataset.cycle,
                        "time": dataset.time,
                    }
                )
                output_cycles.append(dataset.cycle)
                output_times.append(dataset.time)
                if variable == "mhd_w_bcc":
                    mhd_datasets.append(dataset)
            output_index_metadata.append(
                (output_index, output_cycles, output_times)
            )
        _validate_output_index_metadata(
            output_index_metadata,
            terminal_cycle=receipt.get("terminal_cycle"),
            terminal_time=receipt.get("terminal_time"),
        )
        executable = bell.registered_manifest_executable_binding(
            receipt,
            member=member,
            artifact_root=artifact_root,
            authorized_manifest_root=authorized_orion_root,
        )
        dependency = bell.validate_q043_dependency(
            q043_registered_raw_oracle_dependency,
            artifact_root=q043_artifact_root,
        )
        provenance = {
            "kind": "registered_execution_trace",
            "deck_path": (
                bell.DECK_ROOT.relative_to(bell.REPO_ROOT)
                / str(member["deck_path"])
            ).as_posix(),
            "deck_sha256": member["deck_sha256"],
            "source_path": bell.SOURCE_PATH.as_posix(),
            "source_sha256": executable["source_bindings"][
                bell.SOURCE_PATH.as_posix()
            ]["sha256"],
            "corrected_eigenmode_header_path": (
                bell.CORRECTED_EIGENMODE_HEADER_PATH.as_posix()
            ),
            "corrected_eigenmode_header_sha256": executable[
                "source_bindings"
            ][bell.CORRECTED_EIGENMODE_HEADER_PATH.as_posix()]["sha256"],
            "q043_registered_raw_oracle_dependency_sha256": (
                bell._dependency_digest(dependency)
            ),
            "authorized_artifact_root": str(artifact_root),
            "candidate_clean": True,
            "source_clean": True,
            "executable_clean": True,
            "executable_path": executable["path"],
            "executable_sha256": executable["sha256"],
            "registered_execution_receipt_path": (
                "analysis/q023_registered_execution_receipt.json"
            ),
            "registered_execution_receipt_sha256": paired_evidence[
                "registered_execution_receipt"
            ]["sha256"],
            "raw_artifacts": raw_artifacts,
        }
        record = {
            "member_id": member["member_id"],
            "dimension": member["dimension"],
            "epsilon": member["epsilon"],
            "resolution": member["resolution"],
            "decomposition": member["decomposition"],
            "decomposition_splits": member["decomposition_splits"],
            "physics_trace": bell.physics_trace_from_raw_datasets(
                mhd_datasets, member=member
            ),
            "provenance": provenance,
        }
        try:
            report = bell.analyze_rederived_predecessor_record(
                record,
                q043_registered_raw_oracle_dependency=dependency,
                q043_artifact_root=q043_artifact_root,
                artifact_root=artifact_root,
                authorized_manifest_root=authorized_orion_root,
                installed_reconciliation_rederivation=paired_evidence[
                    "installed_reconciliation_rederivation"
                ],
                retained_raw_batch=retained_raw_batch,
            )
            retained_raw_batch.revalidate()
        except bell.ContractError as error:
            raise QualificationError(
                "Q023 retained case failed independently reconstructed analyzer admission"
            ) from error
        return record, report, receipt, paired_evidence
    finally:
        retained_raw_batch.close()


def derive_predecessor_record_from_registered_execution(
    *,
    member_id: str,
    artifact_root: Path,
    q043_registered_raw_oracle_dependency: Mapping[str, object],
    q043_artifact_root: Path,
    authorized_orion_root: Path = prep.AUTHORIZED_ORION_ROOT,
    authorized_project_home_root: Path = prep.CANONICAL_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    """Derive one Q023 predecessor record solely from sealed registered evidence."""
    record, _, _, _ = _derive_registered_case_evidence(
        member_id=member_id,
        artifact_root=artifact_root,
        q043_registered_raw_oracle_dependency=(
            q043_registered_raw_oracle_dependency
        ),
        q043_artifact_root=q043_artifact_root,
        authorized_orion_root=authorized_orion_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    return record


def build_case_admission(
    *,
    member_id: str,
    artifact_root: Path,
    q043_registered_raw_oracle_dependency: Mapping[str, object],
    q043_artifact_root: Path,
    authorized_orion_root: Path = prep.AUTHORIZED_ORION_ROOT,
    authorized_project_home_root: Path = prep.CANONICAL_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    """Rebuild one registered Q023 case admission from immutable retained evidence."""
    artifact_root = _validated_case_root(
        artifact_root=artifact_root,
        authorized_orion_root=authorized_orion_root,
    )
    record, report, receipt, paired_evidence = _derive_registered_case_evidence(
        member_id=member_id,
        artifact_root=artifact_root,
        q043_registered_raw_oracle_dependency=(
            q043_registered_raw_oracle_dependency
        ),
        q043_artifact_root=q043_artifact_root,
        authorized_orion_root=authorized_orion_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    member_id = str(report["member_id"])
    _require(
        receipt.get("registered_science_authorization_id")
        == _expected_authorization_id(member_id),
        "Q023 registered authorization ID is not the deterministic policy slice",
    )
    dependency = bell.validate_q043_dependency(
        q043_registered_raw_oracle_dependency,
        artifact_root=q043_artifact_root,
    )
    physics_record = {
        key: record[key] for key in bell._PHYSICS_MATRIX_RECORD_KEYS
    }
    execution_identity = {
        key: receipt[key]
        for key in (
            "reservation_id",
            "submission_id",
            "reconciliation_event_sha256",
            "reconciliation_mirror_ack_sha256",
            "control_plane_version",
            "registered_science_authorization_id",
            "slurm_job_id",
            "pre_submit_manifest_path",
            "pre_submit_manifest_sha256",
        )
    }
    candidate_binding = {
        key: receipt[key]
        for key in (
            "source_commit",
            "source_bundle_sha256",
            "source_archive_sha256",
            "clean_candidate_manifest_sha256",
            "executable_sha256",
            "environment_sha256",
            "control_plane_version",
        )
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": CASE_RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "campaign_id": bell.CAMPAIGN_ID,
        "status": "registered_execution_case_admitted_non_authorizing",
        "qualification_effect": QUALIFICATION_EFFECT,
        "member_id": member_id,
        "artifact_root": str(artifact_root),
        "q043_dependency_sha256": bell._dependency_digest(dependency),
        "source_record": dict(record),
        "source_record_sha256": canonical_sha256(record),
        "physics_record": physics_record,
        "physics_record_sha256": canonical_sha256(physics_record),
        "analysis_report": report,
        "analysis_report_sha256": canonical_sha256(report),
        "candidate_binding": candidate_binding,
        "candidate_binding_sha256": canonical_sha256(candidate_binding),
        "execution_identity": execution_identity,
        "paired_evidence": paired_evidence,
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def validate_case_admission(
    value: object,
    *,
    q043_registered_raw_oracle_dependency: Mapping[str, object],
    q043_artifact_root: Path,
    authorized_orion_root: Path = prep.AUTHORIZED_ORION_ROOT,
    authorized_project_home_root: Path = prep.CANONICAL_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    _require(type(value) is dict, "Q023 case admission must be an object")
    rebuilt = build_case_admission(
        member_id=str(value.get("member_id", "")),
        artifact_root=Path(str(value.get("artifact_root", ""))),
        q043_registered_raw_oracle_dependency=(
            q043_registered_raw_oracle_dependency
        ),
        q043_artifact_root=q043_artifact_root,
        authorized_orion_root=authorized_orion_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    _require(
        _strict_equal(value, rebuilt),
        "Q023 case admission derived fields or authority drifted",
    )
    return rebuilt


def build_matrix_qualification(
    *,
    case_admissions: Sequence[Mapping[str, object]],
    q043_registered_raw_oracle_dependency: Mapping[str, object],
    q043_artifact_root: Path,
    authorized_orion_root: Path = prep.AUTHORIZED_ORION_ROOT,
    authorized_project_home_root: Path = prep.CANONICAL_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    """Build the exact complete registered Q023 linear qualification."""
    _require(type(case_admissions) is list, "Q023 case admissions must be a list")
    dependency = bell.validate_q043_dependency(
        q043_registered_raw_oracle_dependency,
        artifact_root=q043_artifact_root,
    )
    admissions = [
        validate_case_admission(
            item,
            q043_registered_raw_oracle_dependency=dependency,
            q043_artifact_root=q043_artifact_root,
            authorized_orion_root=authorized_orion_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        for item in case_admissions
    ]
    dependency_sha256 = bell._dependency_digest(dependency)
    _require(
        all(
            item["q043_dependency_sha256"] == dependency_sha256
            and item["authorization"] == AUTHORIZATION_BOUNDARY
            for item in admissions
        ),
        "registered Q023 case dependency or authority boundary drifted",
    )
    expected_ids = list(bell._manifest_members())
    _require(
        [item["member_id"] for item in admissions] == expected_ids,
        "registered Q023 matrix is incomplete, extra, or noncanonical",
    )
    identity_keys = (
        "reservation_id",
        "submission_id",
        "reconciliation_event_sha256",
        "reconciliation_mirror_ack_sha256",
        "registered_science_authorization_id",
        "slurm_job_id",
        "pre_submit_manifest_path",
        "pre_submit_manifest_sha256",
    )
    for key in identity_keys:
        values = [item["execution_identity"][key] for item in admissions]
        _require(
            len(values) == len(set(values)),
            f"registered Q023 execution identity reused: {key}",
        )
    artifact_roots = [item["artifact_root"] for item in admissions]
    _require(
        len(artifact_roots) == len(set(artifact_roots)),
        "registered Q023 artifact root reused",
    )
    _require(
        len({item["candidate_binding_sha256"] for item in admissions}) == 1,
        "registered Q023 matrix mixes clean candidates or controller generations",
    )
    physics = bell.analyze_physics_trace_matrix(
        [item["physics_record"] for item in admissions]
    )
    case_bindings = [
        {
            "member_id": item["member_id"],
            "case_admission_sha256": canonical_sha256(item),
            "source_record_sha256": item["source_record_sha256"],
            "physics_record_sha256": item["physics_record_sha256"],
            "analysis_report_sha256": item["analysis_report_sha256"],
            "candidate_binding_sha256": item["candidate_binding_sha256"],
        }
        for item in admissions
    ]
    registered_pass = bool(physics["predecessor_contract_pass"])
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": MATRIX_RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "campaign_id": bell.CAMPAIGN_ID,
        "status": (
            "complete_registered_Q023_linear_matrix_pass_non_authorizing"
            if registered_pass
            else "complete_registered_Q023_linear_matrix_physics_fail_non_authorizing"
        ),
        "qualification_effect": QUALIFICATION_EFFECT,
        "q043_dependency": dependency,
        "q043_dependency_sha256": dependency_sha256,
        "case_count": len(admissions),
        "case_admissions": admissions,
        "case_bindings": case_bindings,
        "case_bindings_sha256": canonical_sha256(case_bindings),
        "physics_analysis": physics,
        "physics_analysis_sha256": canonical_sha256(physics),
        "registered_execution_qualification_check_pass": True,
        "registered_linear_qualification_pass": registered_pass,
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def validate_matrix_qualification(
    value: object,
    *,
    q043_artifact_root: Path,
    authorized_orion_root: Path = prep.AUTHORIZED_ORION_ROOT,
    authorized_project_home_root: Path = prep.CANONICAL_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    _require(
        type(value) is dict and value.get("record_type") == MATRIX_RECORD_TYPE,
        "Q019 prerequisite requires a registered Q023 matrix record",
    )
    rebuilt = build_matrix_qualification(
        case_admissions=value.get("case_admissions"),
        q043_registered_raw_oracle_dependency=value.get("q043_dependency"),
        q043_artifact_root=q043_artifact_root,
        authorized_orion_root=authorized_orion_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    _require(
        _strict_equal(value, rebuilt),
        "Q023 matrix qualification derived fields or authority drifted",
    )
    return rebuilt


def validate_downstream_q019_prerequisite(
    value: object,
    *,
    q043_artifact_root: Path,
    authorized_orion_root: Path = prep.AUTHORIZED_ORION_ROOT,
    authorized_project_home_root: Path = prep.CANONICAL_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    record = validate_matrix_qualification(
        value,
        q043_artifact_root=q043_artifact_root,
        authorized_orion_root=authorized_orion_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    _require(
        record["case_count"] == prep.EXPECTED_CASE_COUNT
        and record["registered_execution_qualification_check_pass"] is True
        and record["registered_linear_qualification_pass"] is True
        and all(value is False for value in record["authorization"].values()),
        "Q019 prerequisite requires a complete passing non-authorizing Q023 matrix",
    )
    return record
